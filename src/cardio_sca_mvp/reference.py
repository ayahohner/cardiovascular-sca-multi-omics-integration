from __future__ import annotations

import hashlib
import json
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from .util import canonical_json_bytes, write_json, write_tsv


@dataclass
class ReferenceData:
    organ: str
    expression: pd.DataFrame
    groups: pd.DataFrame
    states: pd.DataFrame
    exclusions: pd.DataFrame
    gene_crosswalk: pd.DataFrame


def _config_digest(source: dict[str, Any], analysis: dict[str, Any]) -> str:
    payload = {
        "source_sha256": source["sha256"],
        "analysis": analysis,
        "reference_implementation_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    return hashlib.sha256(canonical_json_bytes(payload)).hexdigest()


def _cache_paths(cache_dir: Path, organ: str) -> dict[str, Path]:
    prefix = cache_dir / f"reference_{organ}"
    return {
        "matrix": prefix.with_suffix(".npz"),
        "groups": prefix.with_name(prefix.name + "_groups.tsv"),
        "states": prefix.with_name(prefix.name + "_states.tsv"),
        "exclusions": prefix.with_name(prefix.name + "_exclusions.tsv"),
        "crosswalk": prefix.with_name(prefix.name + "_gene_crosswalk.tsv"),
        "metadata": prefix.with_name(prefix.name + "_cache.json"),
    }


def load_or_build_reference(
    root: Path,
    source: dict[str, Any],
    reference_spec: dict[str, Any],
    pseudobulk_config: dict[str, Any],
    cache_dir: Path,
    force: bool = False,
) -> ReferenceData:
    cache = _cache_paths(cache_dir, reference_spec["organ"])
    digest = _config_digest(source, {"reference": reference_spec, "pseudobulk": pseudobulk_config})
    if not force and all(path.exists() for path in cache.values()):
        metadata = json.loads(cache["metadata"].read_text(encoding="utf-8"))
        if metadata.get("config_digest") == digest:
            arrays = np.load(cache["matrix"], allow_pickle=False)
            expression = pd.DataFrame(arrays["expression"], columns=arrays["genes"].astype(str))
            return ReferenceData(
                organ=reference_spec["organ"],
                expression=expression,
                groups=pd.read_csv(cache["groups"], sep="\t"),
                states=pd.read_csv(cache["states"], sep="\t"),
                exclusions=pd.read_csv(cache["exclusions"], sep="\t"),
                gene_crosswalk=pd.read_csv(cache["crosswalk"], sep="\t"),
            )

    result = build_reference(
        root / source["local_path"], reference_spec, pseudobulk_config
    )
    cache_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=cache_dir, suffix=".npz", delete=False) as handle:
        temporary = Path(handle.name)
    np.savez_compressed(
        temporary,
        expression=result.expression.to_numpy(dtype=np.float32),
        genes=result.expression.columns.to_numpy(dtype=str),
    )
    os.replace(temporary, cache["matrix"])
    write_tsv(cache["groups"], result.groups)
    write_tsv(cache["states"], result.states)
    write_tsv(cache["exclusions"], result.exclusions)
    write_tsv(cache["crosswalk"], result.gene_crosswalk)
    write_json(
        cache["metadata"],
        {
            "config_digest": digest,
            "source_id": source["id"],
            "source_sha256": source["sha256"],
            "matrix_shape": list(result.expression.shape),
        },
    )
    return result


def build_reference(
    path: Path,
    reference_spec: dict[str, Any],
    config: dict[str, Any],
) -> ReferenceData:
    dataset = ad.read_h5ad(path, backed="r")
    try:
        obs = dataset.obs
        required = {
            reference_spec["state_column"],
            reference_spec["donor_column"],
            reference_spec["sample_column"],
        }
        if reference_spec.get("study_column"):
            required.add(reference_spec["study_column"])
        missing = required - set(obs.columns)
        if missing:
            raise ValueError(f"reference source missing metadata: {', '.join(sorted(missing))}")
        if dataset.raw is None:
            raise ValueError("reference source does not contain raw.X")

        state_column = reference_spec["state_column"]
        source_state = obs[state_column].astype(str)
        mapped_state = source_state.map(reference_spec["state_map"])
        selected = mapped_state.notna().to_numpy()
        selected_positions = np.flatnonzero(selected)

        metadata = pd.DataFrame(
            {
                "obs_position": selected_positions,
                "source_state": source_state.iloc[selected_positions].to_numpy(),
                "state": mapped_state.iloc[selected_positions].to_numpy(),
                "donor": obs[reference_spec["donor_column"]].astype(str).iloc[selected_positions].to_numpy(),
                "sample": obs[reference_spec["sample_column"]].astype(str).iloc[selected_positions].to_numpy(),
                "study": (
                    obs[reference_spec["study_column"]].astype(str).iloc[selected_positions].to_numpy()
                    if reference_spec.get("study_column")
                    else np.repeat("not_applicable", len(selected_positions))
                ),
            }
        )
        metadata["organ"] = reference_spec["organ"]
        metadata["group_id"] = (
            metadata["organ"] + "::" + metadata["donor"] + "::" + metadata["state"]
        )

        group_sizes = metadata.groupby("group_id", observed=True).size()
        eligible_group_ids = set(
            group_sizes[group_sizes >= int(config["min_cells_per_donor_state"])].index
        )
        low_groups = metadata.loc[~metadata["group_id"].isin(eligible_group_ids)]
        metadata = metadata.loc[metadata["group_id"].isin(eligible_group_ids)].copy()

        state_donors = metadata.groupby("state", observed=True)["donor"].nunique()
        eligible_states = set(
            state_donors[state_donors >= int(config["min_donors_per_state"])].index
        )
        low_state_rows = metadata.loc[~metadata["state"].isin(eligible_states)]
        metadata = metadata.loc[metadata["state"].isin(eligible_states)].copy()
        if metadata.empty:
            raise ValueError("no donor-state groups passed reference eligibility")

        group_ids = sorted(metadata["group_id"].unique())
        group_code = {group_id: index for index, group_id in enumerate(group_ids)}
        codes = metadata["group_id"].map(group_code).to_numpy(dtype=np.int32)
        positions = metadata["obs_position"].to_numpy(dtype=np.int64)

        raw = dataset.raw
        n_features = raw.shape[1]
        accumulated = np.zeros((len(group_ids), n_features), dtype=np.float64)
        batch_cells = int(config["batch_cells"])
        for start in range(0, len(positions), batch_cells):
            stop = min(start + batch_cells, len(positions))
            batch_positions = positions[start:stop]
            batch_codes = codes[start:stop]
            batch = raw.X[batch_positions, :]
            if not sparse.issparse(batch):
                batch = sparse.csr_matrix(batch)
            for code in np.unique(batch_codes):
                total = batch[batch_codes == code].sum(axis=0)
                accumulated[code] += np.asarray(total).ravel()

        feature_symbols = (
            raw.var["feature_name"].astype(str).to_numpy()
            if "feature_name" in raw.var.columns
            else raw.var_names.astype(str).to_numpy()
        )
        feature_ids = raw.var_names.astype(str).to_numpy()
        accumulated, genes, crosswalk = collapse_gene_symbols(
            accumulated, feature_symbols, feature_ids
        )
        library_sizes = accumulated.sum(axis=1)
        if np.any(library_sizes <= 0):
            raise ValueError("one or more pseudobulk groups have zero library size")
        log_cpm = np.log1p(
            accumulated / library_sizes[:, np.newaxis] * float(config["normalization_target"])
        ).astype(np.float32)

        group_metadata = (
            metadata.groupby("group_id", observed=True)
            .agg(
                organ=("organ", "first"),
                state=("state", "first"),
                source_state=("source_state", "first"),
                donor=("donor", "first"),
                study=("study", "first"),
                cell_count=("obs_position", "size"),
                sample_count=("sample", "nunique"),
            )
            .reindex(group_ids)
            .reset_index()
        )
        state_metadata = (
            metadata.groupby(["organ", "state", "source_state"], observed=True)
            .agg(
                cell_count=("obs_position", "size"),
                donor_count=("donor", "nunique"),
                sample_count=("sample", "nunique"),
                study_count=("study", "nunique"),
            )
            .reset_index()
            .sort_values(["organ", "state"])
        )

        exclusions = build_exclusion_table(
            obs,
            source_state,
            reference_spec,
            low_groups,
            low_state_rows,
        )
        return ReferenceData(
            organ=reference_spec["organ"],
            expression=pd.DataFrame(log_cpm, columns=genes),
            groups=group_metadata,
            states=state_metadata,
            exclusions=exclusions,
            gene_crosswalk=crosswalk,
        )
    finally:
        dataset.file.close()


def collapse_gene_symbols(
    counts: np.ndarray,
    feature_symbols: np.ndarray,
    feature_ids: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    clean = np.array([symbol.strip().upper() for symbol in feature_symbols], dtype=object)
    valid = np.array([bool(symbol) and symbol not in {"NAN", "NONE"} for symbol in clean])
    clean = clean[valid]
    ids = feature_ids[valid]
    matrix = counts[:, valid]
    genes = np.array(sorted(set(clean)), dtype=str)
    gene_to_code = {gene: index for index, gene in enumerate(genes)}
    collapsed = np.zeros((counts.shape[0], len(genes)), dtype=np.float64)
    identifiers: dict[str, list[str]] = {gene: [] for gene in genes}
    for column, (gene, feature_id) in enumerate(zip(clean, ids, strict=True)):
        code = gene_to_code[str(gene)]
        collapsed[:, code] += matrix[:, column]
        identifiers[str(gene)].append(str(feature_id))
    crosswalk = pd.DataFrame(
        {
            "gene": genes,
            "feature_ids": [";".join(sorted(set(identifiers[gene]))) for gene in genes],
            "feature_count": [len(set(identifiers[gene])) for gene in genes],
        }
    )
    return collapsed, genes, crosswalk


def build_exclusion_table(
    obs: pd.DataFrame,
    source_state: pd.Series,
    reference_spec: dict[str, Any],
    low_groups: pd.DataFrame,
    low_state_rows: pd.DataFrame,
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    donor_column = reference_spec["donor_column"]
    for state, reason in reference_spec.get("excluded_states", {}).items():
        mask = source_state.eq(state)
        records.append(
            {
                "organ": reference_spec["organ"],
                "source_state": state,
                "level": "state",
                "cell_count": int(mask.sum()),
                "donor_count": int(obs.loc[mask, donor_column].nunique()),
                "reason": reason,
                "response": "excluded",
            }
        )
    for state, reason in reference_spec.get("unavailable_expected_states", {}).items():
        records.append(
            {
                "organ": reference_spec["organ"],
                "source_state": state,
                "level": "expected_state",
                "cell_count": 0,
                "donor_count": 0,
                "reason": reason,
                "response": "abstained",
            }
        )
    if not low_groups.empty:
        for group_id, group in low_groups.groupby("group_id", observed=True):
            records.append(
                {
                    "organ": reference_spec["organ"],
                    "source_state": group["source_state"].iloc[0],
                    "level": "donor_state",
                    "cell_count": int(len(group)),
                    "donor_count": 1,
                    "reason": "below min_cells_per_donor_state",
                    "response": "excluded",
                }
            )
    if not low_state_rows.empty:
        for state, group in low_state_rows.groupby("source_state", observed=True):
            records.append(
                {
                    "organ": reference_spec["organ"],
                    "source_state": state,
                    "level": "state",
                    "cell_count": int(len(group)),
                    "donor_count": int(group["donor"].nunique()),
                    "reason": "below min_donors_per_state",
                    "response": "excluded",
                }
            )
    columns = [
        "organ",
        "source_state",
        "level",
        "cell_count",
        "donor_count",
        "reason",
        "response",
    ]
    return pd.DataFrame.from_records(records, columns=columns)
