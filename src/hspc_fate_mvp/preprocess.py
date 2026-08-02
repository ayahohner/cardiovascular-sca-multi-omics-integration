from __future__ import annotations

import gzip
import json
import warnings
from dataclasses import replace
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
import scrublet as scr
from scipy import sparse

from .matrixio import DonorMatrix, read_discovery_zip, read_validation_sample
from .util import canonical_sha256, sha256_file, write_json


Reader = Callable[[], DonorMatrix]


def _knee_count(total_counts: np.ndarray) -> int:
    nonzero = np.sort(total_counts[total_counts > 0])[::-1]
    if len(nonzero) < 3:
        return 0
    ranks = np.arange(1, len(nonzero) + 1)
    x = np.log10(ranks)
    y = np.log10(nonzero)
    line = y[0] + (y[-1] - y[0]) * (x - x[0]) / (x[-1] - x[0])
    distance = y - line
    return int(nonzero[int(np.argmax(distance))])


def _scrublet_calls(
    counts: sparse.csr_matrix, expected_rate: float, seed: int
) -> tuple[np.ndarray, np.ndarray, float, str]:
    if counts.shape[0] < 100 or counts.shape[1] < 100:
        raise RuntimeError("Scrublet requires at least 100 cells and 100 genes")
    model = scr.Scrublet(counts, expected_doublet_rate=expected_rate, random_state=seed)
    with warnings.catch_warnings(record=True) as caught_warnings:
        warnings.simplefilter("always", RuntimeWarning)
        scores, calls = model.scrub_doublets(
            use_approx_neighbors=False,
            min_counts=3,
            min_cells=3,
            min_gene_variability_pctl=85,
            log_transform=True,
            n_prin_comps=min(30, counts.shape[0] - 2),
            verbose=False,
        )
    runtime_warning_count = sum(issubclass(item.category, RuntimeWarning) for item in caught_warnings)
    threshold_mode = "scrublet_automatic"
    if calls is None or len(np.unique(scores)) < 10:
        totals = np.asarray(counts.sum(axis=1)).ravel()
        genes = np.diff(counts.indptr)

        def robust_z(values: np.ndarray) -> np.ndarray:
            logged = np.log1p(values)
            median = np.median(logged)
            mad = np.median(np.abs(logged - median))
            return (logged - median) / max(1.4826 * mad, np.finfo(float).eps)

        joint_score = np.minimum(robust_z(totals), robust_z(genes))
        threshold = 3.5
        calls = joint_score > threshold
        scores = joint_score
        threshold_mode = "robust_joint_umi_gene_outlier"
    elif runtime_warning_count:
        threshold_mode = "scrublet_automatic_with_runtime_warning"
    if scores.shape != (counts.shape[0],) or calls is None or calls.shape != (counts.shape[0],):
        raise RuntimeError("Scrublet returned invalid result dimensions")
    if not np.isfinite(scores).all():
        raise RuntimeError("Scrublet returned non-finite scores")
    threshold = float(model.threshold_) if model.threshold_ is not None else float("nan")
    return scores, calls.astype(bool), threshold, threshold_mode


def qc_donor(
    donor: DonorMatrix,
    rules: dict[str, Any],
    seed: int,
    doublet_caller: Callable[
        [sparse.csr_matrix, float, int], tuple[np.ndarray, np.ndarray, float, str]
    ] = _scrublet_calls,
) -> tuple[DonorMatrix, list[dict[str, Any]], dict[str, Any]]:
    counts = donor.counts.tocsr()
    total = np.asarray(counts.sum(axis=1)).ravel()
    detected = np.diff(counts.indptr)
    mitochondrial = np.char.startswith(np.char.upper(donor.genes.astype(str)), "MT-")
    mito_counts = np.asarray(counts[:, mitochondrial].sum(axis=1)).ravel()
    mito_fraction = np.divide(mito_counts, total, out=np.ones_like(total, dtype=float), where=total > 0)
    cd34_index = np.flatnonzero(np.char.upper(donor.genes.astype(str)) == "CD34")
    cd34_counts = (
        np.asarray(counts[:, cd34_index].sum(axis=1)).ravel() if len(cd34_index) else np.zeros(len(total))
    )
    knee = _knee_count(total)

    masks: list[tuple[str, np.ndarray]] = [
        ("nonzero", total > 0),
        ("minimum_umis", total >= int(rules["minimum_umis"])),
        ("minimum_genes", detected >= int(rules["minimum_genes"])),
        ("maximum_mitochondrial_fraction", mito_fraction <= float(rules["maximum_mitochondrial_fraction"])),
    ]
    if rules["require_cd34_expression"]:
        if not len(cd34_index):
            raise RuntimeError(f"CD34 is absent from {donor.donor}")
        masks.append(("cd34_expression", cd34_counts > 0))

    cumulative = np.ones(counts.shape[0], dtype=bool)
    ledger: list[dict[str, Any]] = []
    ledger.append(
        {
            "cohort": donor.cohort,
            "donor": donor.donor,
            "condition": donor.condition,
            "control": "raw_barcodes",
            "threshold": "not_applicable",
            "input_cells": counts.shape[0],
            "retained_cells": counts.shape[0],
            "excluded_cells": 0,
        }
    )
    thresholds = {
        "nonzero": ">0 UMI",
        "minimum_umis": f">={rules['minimum_umis']} UMI",
        "minimum_genes": f">={rules['minimum_genes']} genes",
        "maximum_mitochondrial_fraction": f"<={rules['maximum_mitochondrial_fraction']}",
        "cd34_expression": ">0 CD34 UMI",
    }
    for control, mask in masks:
        before = int(cumulative.sum())
        cumulative &= mask
        after = int(cumulative.sum())
        ledger.append(
            {
                "cohort": donor.cohort,
                "donor": donor.donor,
                "condition": donor.condition,
                "control": control,
                "threshold": thresholds[control],
                "input_cells": before,
                "retained_cells": after,
                "excluded_cells": before - after,
            }
        )
    selected = counts[cumulative].tocsr()
    selected_barcodes = donor.barcodes[cumulative]
    scores, doublets, scrublet_threshold, scrublet_threshold_mode = doublet_caller(
        selected, float(rules["scrublet_expected_doublet_rate"]), seed
    )
    keep = ~doublets
    ledger.append(
        {
            "cohort": donor.cohort,
            "donor": donor.donor,
            "condition": donor.condition,
            "control": "scrublet_doublet",
            "threshold": scrublet_threshold,
            "input_cells": len(keep),
            "retained_cells": int(keep.sum()),
            "excluded_cells": int(doublets.sum()),
        }
    )
    retained = replace(donor, counts=selected[keep].tocsr(), barcodes=selected_barcodes[keep])
    diagnostics = {
        "raw_barcodes": int(counts.shape[0]),
        "raw_nonzero_barcodes": int((total > 0).sum()),
        "knee_umi_count": knee,
        "frozen_umi_threshold": int(rules["minimum_umis"]),
        "scrublet_threshold": scrublet_threshold,
        "scrublet_threshold_mode": scrublet_threshold_mode,
        "scrublet_score_median": float(np.median(scores)),
        "scrublet_score_maximum": float(np.max(scores)),
        "retained_cells": int(retained.counts.shape[0]),
    }
    return retained, ledger, diagnostics


def _cache_paths(root: Path, cohort: str, donor: str) -> tuple[Path, Path, Path, Path]:
    base = root / "data/derived/hspc-fate/preprocessed" / cohort / donor
    return (
        base.with_suffix(".npz"),
        base.with_suffix(".genes.txt.gz"),
        base.with_suffix(".barcodes.txt.gz"),
        base.with_suffix(".json"),
    )


def _save_cache(root: Path, donor: DonorMatrix, metadata: dict[str, Any]) -> None:
    matrix_path, genes_path, barcodes_path, metadata_path = _cache_paths(root, donor.cohort, donor.donor)
    matrix_path.parent.mkdir(parents=True, exist_ok=True)
    sparse.save_npz(matrix_path, donor.counts, compressed=True)
    with gzip.open(genes_path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(donor.genes) + "\n")
    with gzip.open(barcodes_path, "wt", encoding="utf-8") as handle:
        handle.write("\n".join(donor.barcodes) + "\n")
    write_json(metadata_path, metadata)


def _load_cache(root: Path, cohort: str, donor: str, condition: str) -> tuple[DonorMatrix, dict[str, Any]]:
    matrix_path, genes_path, barcodes_path, metadata_path = _cache_paths(root, cohort, donor)
    counts = sparse.load_npz(matrix_path).tocsr()
    with gzip.open(genes_path, "rt", encoding="utf-8") as handle:
        genes = np.asarray([line.rstrip("\n") for line in handle], dtype=str)
    with gzip.open(barcodes_path, "rt", encoding="utf-8") as handle:
        barcodes = np.asarray([line.rstrip("\n") for line in handle], dtype=str)
    with metadata_path.open(encoding="utf-8") as handle:
        metadata = json.load(handle)
    return DonorMatrix(donor, condition, cohort, counts, genes, barcodes), metadata


def prepare_cohort(
    root: Path,
    config: dict[str, Any],
    cohort: str,
    force: bool = False,
    doublet_caller: Callable[
        [sparse.csr_matrix, float, int], tuple[np.ndarray, np.ndarray, float, str]
    ] = _scrublet_calls,
) -> tuple[list[DonorMatrix], pd.DataFrame, list[dict[str, Any]]]:
    rules = config["qc"][cohort]
    config_hash = canonical_sha256(
        {
            "qc": rules,
            "seed": config["random_seed"],
            "preprocess_code_sha256": sha256_file(Path(__file__)),
            "matrixio_code_sha256": sha256_file(Path(__file__).with_name("matrixio.py")),
        }
    )
    readers: list[tuple[str, str, str, Reader]] = []
    if cohort == "discovery":
        for asset in config["required_assets"]:
            if asset["cohort"] != "discovery":
                continue
            path = root / asset["local_path"]
            readers.append(
                (
                    asset["donor"],
                    asset["condition"],
                    asset["sha256"],
                    lambda path=path, asset=asset: read_discovery_zip(path, asset["donor"], asset["condition"]),
                )
            )
    elif cohort == "validation":
        asset = next(item for item in config["required_assets"] if item["cohort"] == "validation")
        path = root / asset["local_path"]
        for sample in config["validation_samples"]:
            if not sample["include"]:
                continue
            readers.append(
                (
                    sample["donor"],
                    sample["condition"],
                    asset["sha256"],
                    lambda path=path, sample=sample: read_validation_sample(
                        path, sample["archive"], sample["donor"], sample["condition"]
                    ),
                )
            )
    else:
        raise ValueError(f"unknown cohort: {cohort}")

    donors: list[DonorMatrix] = []
    ledger_rows: list[dict[str, Any]] = []
    diagnostics: list[dict[str, Any]] = []
    for index, (donor_name, condition, source_hash, reader) in enumerate(readers):
        paths = _cache_paths(root, cohort, donor_name)
        if not force and all(path.is_file() for path in paths):
            cached, metadata = _load_cache(root, cohort, donor_name, condition)
            if metadata.get("source_sha256") == source_hash and metadata.get("config_sha256") == config_hash:
                donors.append(cached)
                ledger_rows.extend(metadata["qc_ledger"])
                diagnostics.append(metadata["diagnostics"])
                continue
        raw = reader()
        retained, rows, donor_diagnostics = qc_donor(
            raw, rules, int(config["random_seed"]) + index, doublet_caller=doublet_caller
        )
        metadata = {
            "source_sha256": source_hash,
            "config_sha256": config_hash,
            "qc_ledger": rows,
            "diagnostics": donor_diagnostics,
        }
        _save_cache(root, retained, metadata)
        donors.append(retained)
        ledger_rows.extend(rows)
        diagnostics.append(donor_diagnostics)

    retained_total = sum(donor.counts.shape[0] for donor in donors)
    published = int(rules["published_retained_cells"])
    tolerance = float(rules["published_count_tolerance_fraction"])
    if abs(retained_total - published) / published > tolerance:
        raise RuntimeError(
            f"{cohort} retained-cell total {retained_total} is outside {tolerance:.0%} of published {published}"
        )
    return donors, pd.DataFrame.from_records(ledger_rows), diagnostics
