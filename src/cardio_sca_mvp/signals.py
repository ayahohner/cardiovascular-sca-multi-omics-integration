from __future__ import annotations

import csv
import gzip
import io
import tarfile
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread


def build_signal_ledger(
    root: Path,
    sources: dict[str, dict[str, Any]],
    config: dict[str, Any],
) -> pd.DataFrame:
    huvec = huvec_signature(
        root / sources["gse274340-huvec"]["local_path"], config
    )
    pbmc = pbmc_signature(
        root / sources["gse53441-pbmc"]["local_path"],
        root / sources["gpl570-annotation"]["local_path"],
        config,
    )
    ledger = pd.concat([huvec, pbmc], ignore_index=True)
    ledger.insert(0, "signal_id", [f"SIG-{index:04d}" for index in range(1, len(ledger) + 1)])
    return ledger.sort_values(["evidence_role", "direction", "gene"]).reset_index(drop=True)


def huvec_signature(path: Path, config: dict[str, Any]) -> pd.DataFrame:
    with tarfile.open(path, mode="r") as archive:
        control = _read_10x_tar_sample(archive, "GSM8447652_C1G")
        perturbed = _read_10x_tar_sample(archive, "GSM8447653_N1G")
    merged = control.merge(perturbed, on="gene", how="outer", suffixes=("_control", "_perturbed")).fillna(0)
    control_total = float(merged["counts_control"].sum())
    perturbed_total = float(merged["counts_perturbed"].sum())
    control_cpm = merged["counts_control"] / control_total * 1_000_000.0
    perturbed_cpm = merged["counts_perturbed"] / perturbed_total * 1_000_000.0
    merged["effect"] = np.log2(perturbed_cpm + 1.0) - np.log2(control_cpm + 1.0)
    merged["total_counts"] = merged["counts_control"] + merged["counts_perturbed"]
    eligible = merged.loc[
        (merged["total_counts"] >= float(config["huvec_min_total_counts"]))
        & (merged["effect"].abs() >= float(config["huvec_min_abs_log2_fold_change"]))
    ].copy()
    selected = _top_each_direction(
        eligible,
        effect_column="effect",
        maximum=int(config["max_genes_per_direction"]),
    )
    selected["source_id"] = "gse274340-huvec"
    selected["direction"] = np.where(selected["effect"] > 0, "higher_in_pdRBC_flow", "lower_in_pdRBC_flow")
    selected["fdr"] = np.nan
    selected["evidence_role"] = "endothelial_perturbational_prior"
    selected["replication"] = "one sample per condition; descriptive only"
    selected["transport_limit"] = "HUVEC microfluidic model, not adult heart or lung endothelium"
    selected["permitted_interpretation"] = "program compatibility with an in-vitro mechanical perturbation"
    return selected[
        [
            "source_id",
            "gene",
            "effect",
            "direction",
            "fdr",
            "evidence_role",
            "replication",
            "transport_limit",
            "permitted_interpretation",
        ]
    ]


def _read_10x_tar_sample(archive: tarfile.TarFile, prefix: str) -> pd.DataFrame:
    feature_member = archive.getmember(f"{prefix}_features.tsv.gz")
    matrix_member = archive.getmember(f"{prefix}_matrix.mtx.gz")
    feature_handle = archive.extractfile(feature_member)
    matrix_handle = archive.extractfile(matrix_member)
    if feature_handle is None or matrix_handle is None:
        raise ValueError(f"incomplete 10x sample in archive: {prefix}")
    with gzip.open(feature_handle, mode="rt") as features_stream:
        features = pd.read_csv(
            features_stream,
            sep="\t",
            header=None,
            names=["feature_id", "gene", "feature_type"],
        )
    with gzip.open(matrix_handle, mode="rb") as matrix_stream:
        matrix = mmread(matrix_stream, spmatrix=True)
    if not sparse.issparse(matrix):
        matrix = sparse.coo_matrix(matrix)
    counts = np.asarray(matrix.sum(axis=1)).ravel()
    if len(counts) != len(features):
        raise ValueError(f"feature/matrix mismatch for {prefix}")
    frame = pd.DataFrame(
        {
            "gene": features["gene"].astype(str).str.strip().str.upper(),
            "counts": counts,
        }
    )
    frame = frame.loc[frame["gene"].ne("") & frame["gene"].ne("NAN")]
    return frame.groupby("gene", as_index=False, observed=True)["counts"].sum()


def pbmc_signature(
    matrix_path: Path,
    annotation_path: Path,
    config: dict[str, Any],
) -> pd.DataFrame:
    expression, sample_metadata = _read_geo_series_matrix(matrix_path)
    annotation = _read_gpl_annotation(annotation_path)
    mapped = expression.join(annotation.set_index("probe_id"), how="inner")
    mapped = mapped.loc[mapped["gene"].notna()]
    sample_columns = list(expression.columns)
    gene_expression = mapped.groupby("gene", observed=True)[sample_columns].median()

    diagnosis = sample_metadata.set_index("accession")["diagnosis"].reindex(sample_columns)
    case_columns = diagnosis[diagnosis.eq("sickle cell anemia (SCA)")].index.tolist()
    control_columns = diagnosis[diagnosis.eq("normal")].index.tolist()
    if len(case_columns) != 24 or len(control_columns) != 10:
        raise ValueError(
            f"unexpected GSE53441 group sizes: {len(case_columns)} SCA, {len(control_columns)} control"
        )
    case = gene_expression[case_columns].to_numpy(dtype=float)
    control = gene_expression[control_columns].to_numpy(dtype=float)
    test = stats.ttest_ind(case, control, axis=1, equal_var=False, nan_policy="omit")
    effect = _hedges_g(case, control)
    fdr = benjamini_hochberg(np.asarray(test.pvalue, dtype=float))
    results = pd.DataFrame(
        {
            "gene": gene_expression.index.astype(str),
            "effect": effect,
            "fdr": fdr,
        }
    ).replace([np.inf, -np.inf], np.nan).dropna(subset=["effect", "fdr"])
    eligible = results.loc[
        (results["fdr"] <= float(config["pbmc_fdr_threshold"]))
        & (results["effect"].abs() >= float(config["pbmc_min_abs_standardized_effect"]))
    ].copy()
    selected = _top_each_direction(
        eligible,
        effect_column="effect",
        maximum=int(config["max_genes_per_direction"]),
    )
    selected["source_id"] = "gse53441-pbmc"
    selected["direction"] = np.where(selected["effect"] > 0, "higher_in_SCA_PBMC", "lower_in_SCA_PBMC")
    selected["evidence_role"] = "immune_association"
    selected["replication"] = "24 SCA and 10 control GEO biological-replicate samples"
    selected["transport_limit"] = "bulk PBMC association, not endothelial identity or a disease effect in heart or lung"
    selected["permitted_interpretation"] = "systemic immune program association compatible with a healthy endothelial-state contrast"
    return selected[
        [
            "source_id",
            "gene",
            "effect",
            "direction",
            "fdr",
            "evidence_role",
            "replication",
            "transport_limit",
            "permitted_interpretation",
        ]
    ]


def _read_geo_series_matrix(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    headers: dict[str, list[str]] = {}
    table_lines: list[str] = []
    in_table = False
    with gzip.open(path, mode="rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.rstrip("\n")
            if stripped == "!series_matrix_table_begin":
                in_table = True
                continue
            if stripped == "!series_matrix_table_end":
                break
            if in_table:
                table_lines.append(line)
            elif stripped.startswith("!Sample_"):
                row = next(csv.reader([stripped], delimiter="\t", quotechar='"'))
                headers.setdefault(row[0], []).append(row[1:])
    expression = pd.read_csv(io.StringIO("".join(table_lines)), sep="\t", index_col=0)
    expression.index = expression.index.astype(str)
    expression.columns = [column.strip('"') for column in expression.columns]
    accessions = headers["!Sample_geo_accession"][0]
    titles = headers["!Sample_title"][0]
    characteristic_rows = headers["!Sample_characteristics_ch1"]
    diagnoses = [value.removeprefix("diagnosis: ") for value in characteristic_rows[0]]
    metadata = pd.DataFrame(
        {"accession": accessions, "title": titles, "diagnosis": diagnoses}
    )
    return expression, metadata


def _read_gpl_annotation(path: Path) -> pd.DataFrame:
    table_lines: list[str] = []
    found_header = False
    with gzip.open(path, mode="rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not found_header and line.startswith("ID\t"):
                found_header = True
            if found_header:
                table_lines.append(line)
    annotation = pd.read_csv(io.StringIO("".join(table_lines)), sep="\t", dtype=str)
    annotation = annotation.rename(columns={"ID": "probe_id", "Gene symbol": "gene"})
    annotation["gene"] = annotation["gene"].map(_unambiguous_gene_symbol)
    return annotation[["probe_id", "gene"]]


def _unambiguous_gene_symbol(value: object) -> str | None:
    if not isinstance(value, str) or value in {"", "---", "nan"}:
        return None
    symbols = {part.strip().upper() for part in value.split("///") if part.strip() not in {"", "---"}}
    if len(symbols) != 1:
        return None
    return symbols.pop()


def _hedges_g(case: np.ndarray, control: np.ndarray) -> np.ndarray:
    n_case = case.shape[1]
    n_control = control.shape[1]
    mean_difference = np.nanmean(case, axis=1) - np.nanmean(control, axis=1)
    pooled_variance = (
        (n_case - 1) * np.nanvar(case, axis=1, ddof=1)
        + (n_control - 1) * np.nanvar(control, axis=1, ddof=1)
    ) / (n_case + n_control - 2)
    correction = 1.0 - 3.0 / (4.0 * (n_case + n_control) - 9.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return correction * mean_difference / np.sqrt(pooled_variance)


def benjamini_hochberg(p_values: np.ndarray) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    result = np.full(values.shape, np.nan, dtype=float)
    finite = np.isfinite(values)
    ordered_indices = np.argsort(values[finite])
    ordered = values[finite][ordered_indices]
    ranks = np.arange(1, len(ordered) + 1)
    adjusted = ordered * len(ordered) / ranks
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    finite_indices = np.flatnonzero(finite)[ordered_indices]
    result[finite_indices] = adjusted
    return result


def _top_each_direction(
    frame: pd.DataFrame,
    effect_column: str,
    maximum: int,
) -> pd.DataFrame:
    positive = frame.loc[frame[effect_column] > 0].nlargest(maximum, effect_column)
    negative = frame.loc[frame[effect_column] < 0].nsmallest(maximum, effect_column)
    return pd.concat([positive, negative], ignore_index=True)
