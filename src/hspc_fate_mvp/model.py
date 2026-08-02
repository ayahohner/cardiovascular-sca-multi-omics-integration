from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Any

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.special import softmax
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

from .matrixio import DonorMatrix
from .util import canonical_sha256


LINEAGES = ("myeloid", "erythroid_megakaryocytic", "lymphoid")


@dataclass
class CombinedCohort:
    counts: sparse.csr_matrix
    genes: np.ndarray
    observations: pd.DataFrame


@dataclass
class FateModel:
    genes: np.ndarray
    variable_gene_indices: np.ndarray
    pca: PCA
    state_model: MiniBatchKMeans
    state_centroids: np.ndarray
    state_absorption: np.ndarray
    terminal_states: dict[str, int]
    root_state: int
    program_centers: dict[str, float]
    program_scales: dict[str, float]
    program_references: dict[str, np.ndarray]
    root_score_threshold: float
    root_maturity_threshold: float
    maturity_low: float
    maturity_high: float
    root_centroid: np.ndarray
    terminal_centroids: np.ndarray
    projection_temperature: float
    model_hashes: dict[str, str]
    state_definitions: dict[str, Any]
    topology_controls: dict[str, Any]


@dataclass
class FateResult:
    observations: pd.DataFrame
    root_mask: np.ndarray
    primary_probabilities: np.ndarray
    secondary_probabilities: np.ndarray
    program_scores: dict[str, np.ndarray]
    pseudotime: np.ndarray
    model_hashes_before: dict[str, str]
    model_hashes_after: dict[str, str]


def combine_donors(donors: list[DonorMatrix]) -> CombinedCohort:
    if not donors:
        raise ValueError("cannot combine an empty donor list")
    genes = donors[0].genes
    for donor in donors[1:]:
        if not np.array_equal(genes, donor.genes):
            raise ValueError(f"gene order differs for donor {donor.donor}")
    counts = sparse.vstack([donor.counts for donor in donors], format="csr")
    observations = pd.concat(
        [
            pd.DataFrame(
                {
                    "cohort": donor.cohort,
                    "donor": donor.donor,
                    "condition": donor.condition,
                    "cell_index": np.arange(donor.counts.shape[0], dtype=int),
                }
            )
            for donor in donors
        ],
        ignore_index=True,
    )
    if len(observations) != counts.shape[0]:
        raise RuntimeError("observation count does not match expression rows")
    return CombinedCohort(counts, genes, observations)


def log_normalize(counts: sparse.csr_matrix) -> sparse.csr_matrix:
    normalized = counts.astype(np.float32, copy=True)
    totals = np.asarray(normalized.sum(axis=1)).ravel()
    scale = np.divide(10_000.0, totals, out=np.zeros_like(totals, dtype=np.float32), where=totals > 0)
    normalized = sparse.diags(scale) @ normalized
    normalized.data = np.log1p(normalized.data)
    return normalized.tocsr()


def _gene_lookup(genes: np.ndarray) -> dict[str, int]:
    lookup: dict[str, int] = {}
    for index, gene in enumerate(genes):
        lookup.setdefault(str(gene).upper(), index)
    return lookup


def _raw_program_scores(
    expression: sparse.csr_matrix, genes: np.ndarray, marker_sets: dict[str, list[str]]
) -> tuple[dict[str, np.ndarray], dict[str, list[str]]]:
    lookup = _gene_lookup(genes)
    scores: dict[str, np.ndarray] = {}
    present: dict[str, list[str]] = {}
    for program, markers in marker_sets.items():
        available = [marker for marker in markers if marker.upper() in lookup]
        if len(available) < max(3, len(markers) // 2):
            raise RuntimeError(f"too few markers present for {program}: {available}")
        indices = [lookup[marker.upper()] for marker in available]
        scores[program] = np.asarray(expression[:, indices].mean(axis=1)).ravel()
        present[program] = available
    return scores, present


def _fit_program_standardization(
    raw_scores: dict[str, np.ndarray], healthy_mask: np.ndarray
) -> tuple[dict[str, float], dict[str, float]]:
    centers: dict[str, float] = {}
    scales: dict[str, float] = {}
    for program, values in raw_scores.items():
        reference = values[healthy_mask]
        center = float(np.median(reference))
        mad_scale = float(1.4826 * np.median(np.abs(reference - center)))
        standard_scale = float(np.std(reference))
        centers[program] = center
        scales[program] = max(mad_scale, standard_scale * 0.1, np.finfo(float).eps)
    return centers, scales


def _standardize_programs(
    raw_scores: dict[str, np.ndarray], centers: dict[str, float], scales: dict[str, float]
) -> dict[str, np.ndarray]:
    return {
        program: (values - centers[program]) / scales[program]
        for program, values in raw_scores.items()
    }


def _fit_program_references(
    standardized_scores: dict[str, np.ndarray], healthy_mask: np.ndarray
) -> dict[str, np.ndarray]:
    return {
        program: np.sort(values[healthy_mask].astype(float, copy=True))
        for program, values in standardized_scores.items()
    }


def _percentile_programs(
    standardized_scores: dict[str, np.ndarray], references: dict[str, np.ndarray]
) -> dict[str, np.ndarray]:
    return {
        program: (
            np.searchsorted(references[program], values, side="left")
            + np.searchsorted(references[program], values, side="right")
        )
        / (2 * len(references[program]))
        for program, values in standardized_scores.items()
    }


def _select_variable_genes(
    expression: sparse.csr_matrix,
    genes: np.ndarray,
    healthy_mask: np.ndarray,
    count: int,
) -> np.ndarray:
    healthy = expression[healthy_mask]
    means = np.asarray(healthy.mean(axis=0)).ravel()
    means_squared = np.asarray(healthy.power(2).mean(axis=0)).ravel()
    variances = np.maximum(means_squared - means**2, 0)
    upper = np.char.upper(genes.astype(str))
    eligible = ~(
        np.char.startswith(upper, "MT-")
        | np.char.startswith(upper, "RPL")
        | np.char.startswith(upper, "RPS")
    )
    eligible_indices = np.flatnonzero(eligible & (means > 0))
    if len(eligible_indices) < count:
        raise RuntimeError("too few expressed genes for the frozen variable-gene count")
    order = eligible_indices[np.argsort(variances[eligible_indices], kind="stable")]
    return np.sort(order[-count:])


def _empirical_pseudotime(maturity: np.ndarray, healthy_mask: np.ndarray) -> tuple[np.ndarray, float, float]:
    low, high = np.quantile(maturity[healthy_mask], [0.01, 0.99])
    if high <= low:
        raise RuntimeError("healthy maturity score has no usable range")
    pseudotime = np.clip((maturity - low) / (high - low), 0.0, 1.0)
    return pseudotime, float(low), float(high)


def _state_summary(
    state_labels: np.ndarray,
    healthy_mask: np.ndarray,
    observations: pd.DataFrame,
    pseudotime: np.ndarray,
    program_scores: dict[str, np.ndarray],
    state_count: int,
) -> pd.DataFrame:
    rows = []
    donor_values = observations["donor"].to_numpy()
    for state in range(state_count):
        mask = healthy_mask & (state_labels == state)
        if not mask.any():
            raise RuntimeError(f"empty healthy topology state: {state}")
        row: dict[str, Any] = {
            "state": state,
            "cells": int(mask.sum()),
            "donors": int(pd.Series(donor_values[mask]).nunique()),
            "pseudotime": float(np.mean(pseudotime[mask])),
        }
        for program, values in program_scores.items():
            row[f"{program}_score"] = float(np.mean(values[mask]))
        rows.append(row)
    return pd.DataFrame.from_records(rows)


def _choose_terminal_states(
    states: pd.DataFrame, minimum_donors: int
) -> dict[str, int]:
    candidates: dict[str, list[int]] = {}
    for lineage in LINEAGES:
        other_columns = [f"{other}_score" for other in LINEAGES if other != lineage]
        eligible = states.loc[
            (states["donors"] >= minimum_donors)
            & (states["pseudotime"] >= 0.5)
            & (states[f"{lineage}_score"] > states[other_columns].max(axis=1))
        ]
        if eligible.empty:
            raise RuntimeError(f"no donor-supported terminal state for {lineage}")
        candidates[lineage] = eligible.sort_values(f"{lineage}_score", ascending=False)["state"].tolist()

    best: tuple[float, tuple[int, int, int]] | None = None
    for selection in itertools.product(*(candidates[lineage] for lineage in LINEAGES)):
        if len(set(selection)) != len(LINEAGES):
            continue
        score = sum(
            float(states.loc[states["state"] == state, f"{lineage}_score"].iloc[0])
            for lineage, state in zip(LINEAGES, selection, strict=True)
        )
        if best is None or score > best[0]:
            best = (score, selection)
    if best is None:
        raise RuntimeError("terminal lineages cannot be assigned to distinct topology states")
    return {lineage: int(state) for lineage, state in zip(LINEAGES, best[1], strict=True)}


def _absorption_probabilities(
    centroids: np.ndarray,
    pseudotime: np.ndarray,
    terminal_states: dict[str, int],
    neighbors: int,
) -> tuple[np.ndarray, dict[int, list[tuple[int, float]]]]:
    terminal_lookup = {state: index for index, state in enumerate(terminal_states.values())}
    effective_time = pseudotime.copy()
    for state in terminal_lookup:
        effective_time[state] = 1.0
    graph: dict[int, list[tuple[int, float]]] = {}
    absorption = np.zeros((len(centroids), len(LINEAGES)), dtype=float)
    for state, lineage_index in terminal_lookup.items():
        absorption[state, lineage_index] = 1.0

    for state in np.argsort(effective_time)[::-1]:
        state = int(state)
        if state in terminal_lookup:
            graph[state] = [(state, 1.0)]
            continue
        forward = np.flatnonzero(effective_time > effective_time[state] + 1e-8)
        if not len(forward):
            forward = np.asarray(list(terminal_lookup), dtype=int)
        distances = np.linalg.norm(centroids[forward] - centroids[state], axis=1)
        chosen_order = np.argsort(distances, kind="stable")[:neighbors]
        chosen = forward[chosen_order]
        chosen_distances = distances[chosen_order]
        distance_scale = max(float(np.median(chosen_distances)), np.finfo(float).eps)
        weights = np.exp(-((chosen_distances / distance_scale) ** 2))
        weights *= np.exp(2.0 * (effective_time[chosen] - effective_time[state]))
        weights /= weights.sum()
        graph[state] = [(int(target), float(weight)) for target, weight in zip(chosen, weights, strict=True)]
        absorption[state] = weights @ absorption[chosen]
    row_sums = absorption.sum(axis=1)
    if not np.allclose(row_sums, 1.0, atol=1e-8):
        failed = np.flatnonzero(~np.isclose(row_sums, 1.0, atol=1e-8))
        raise RuntimeError(f"topology states lack terminal absorption: {failed.tolist()}")
    return absorption, graph


def _model_hashes(
    variable_indices: np.ndarray,
    pca: PCA,
    state_centroids: np.ndarray,
    state_absorption: np.ndarray,
    root_definition: dict[str, Any],
    terminal_definition: dict[str, Any],
    model_parameters: dict[str, Any],
) -> dict[str, str]:
    return {
        "transformation": canonical_sha256(
            {
                "variable_indices": variable_indices.tolist(),
                "pca_mean": pca.mean_.round(10).tolist(),
                "pca_components": pca.components_.round(10).tolist(),
            }
        ),
        "topology": canonical_sha256(
            {
                "state_centroids": state_centroids.round(10).tolist(),
                "state_absorption": state_absorption.round(10).tolist(),
            }
        ),
        "root_definition": canonical_sha256(root_definition),
        "terminal_definitions": canonical_sha256(terminal_definition),
        "model_parameters": canonical_sha256(model_parameters),
    }


def fit_healthy_model(
    cohort: CombinedCohort, config: dict[str, Any]
) -> tuple[FateModel, FateResult]:
    healthy_mask = cohort.observations["condition"].eq("healthy").to_numpy()
    if healthy_mask.sum() == 0:
        raise RuntimeError("healthy reference has no cells")
    expression = log_normalize(cohort.counts)
    raw_scores, present_markers = _raw_program_scores(
        expression, cohort.genes, config["lineage_markers"]
    )
    centers, scales = _fit_program_standardization(raw_scores, healthy_mask)
    scores = _standardize_programs(raw_scores, centers, scales)
    program_references = _fit_program_references(scores, healthy_mask)
    percentile_scores = _percentile_programs(scores, program_references)
    lineage_matrix = np.column_stack([percentile_scores[lineage] for lineage in LINEAGES])
    maturity = lineage_matrix.max(axis=1) - percentile_scores["root"]
    pseudotime, maturity_low, maturity_high = _empirical_pseudotime(maturity, healthy_mask)

    model_config = config["model"]
    root_score_threshold = float(model_config["root_score_quantile"])
    root_maturity_threshold = float(
        np.quantile(maturity[healthy_mask], model_config["root_maturity_quantile"])
    )
    root_mask = (percentile_scores["root"] >= root_score_threshold) & (
        maturity <= root_maturity_threshold
    )
    root_donor_count = int(cohort.observations.loc[root_mask & healthy_mask, "donor"].nunique())
    if root_donor_count < int(model_config["minimum_state_donors"]):
        raise RuntimeError("HSC/MPP root is not supported by enough healthy donors")

    variable_indices = _select_variable_genes(
        expression,
        cohort.genes,
        healthy_mask,
        int(model_config["highly_variable_genes"]),
    )
    dense = expression[:, variable_indices].toarray().astype(np.float32, copy=False)
    pca = PCA(
        n_components=int(model_config["pca_components"]),
        svd_solver="randomized",
        random_state=int(config["random_seed"]),
    )
    pca.fit(dense[healthy_mask])
    coordinates = pca.transform(dense)
    del dense

    state_count = int(model_config["topology_states"])
    state_model = MiniBatchKMeans(
        n_clusters=state_count,
        random_state=int(config["random_seed"]),
        batch_size=2048,
        n_init=10,
    )
    state_model.fit(coordinates[healthy_mask])
    labels = state_model.predict(coordinates)
    states = _state_summary(
        labels,
        healthy_mask,
        cohort.observations,
        pseudotime,
        percentile_scores,
        state_count,
    )
    terminal_states = _choose_terminal_states(states, int(model_config["minimum_state_donors"]))
    terminal_state_values = set(terminal_states.values())
    root_candidates = states.loc[
        (states["donors"] >= int(model_config["minimum_state_donors"]))
        & ~states["state"].isin(terminal_state_values)
    ].copy()
    root_candidates["root_advantage"] = root_candidates["root_score"] - root_candidates[
        [f"{lineage}_score" for lineage in LINEAGES]
    ].max(axis=1)
    root_state = int(root_candidates.sort_values("root_advantage", ascending=False).iloc[0]["state"])

    state_times = states.sort_values("state")["pseudotime"].to_numpy()
    absorption, graph = _absorption_probabilities(
        state_model.cluster_centers_,
        state_times,
        terminal_states,
        int(model_config["topology_neighbors"]),
    )
    sensitivity_absorption, _ = _absorption_probabilities(
        state_model.cluster_centers_,
        state_times,
        terminal_states,
        int(model_config["topology_neighbors"]) + 2,
    )
    stability_correlation = float(
        np.corrcoef(absorption.ravel(), sensitivity_absorption.ravel())[0, 1]
    )
    marker_order = {
        lineage: bool(
            states.loc[states["state"] == terminal_states[lineage], f"{lineage}_score"].iloc[0]
            > states.loc[states["state"] == root_state, f"{lineage}_score"].iloc[0]
        )
        for lineage in LINEAGES
    }
    root_absorption = absorption[root_state]
    topology_controls = {
        "root_to_all_terminals": bool(np.all(root_absorption > 0)),
        "marker_order": marker_order,
        "root_healthy_donors": root_donor_count,
        "terminal_healthy_donors": {
            lineage: int(states.loc[states["state"] == state, "donors"].iloc[0])
            for lineage, state in terminal_states.items()
        },
        "neighbor_sensitivity_correlation": stability_correlation,
        "minimum_neighbor_sensitivity_correlation": 0.8,
    }
    topology_controls["passed"] = bool(
        topology_controls["root_to_all_terminals"]
        and all(marker_order.values())
        and root_donor_count >= int(model_config["minimum_state_donors"])
        and all(
            count >= int(model_config["minimum_state_donors"])
            for count in topology_controls["terminal_healthy_donors"].values()
        )
        and stability_correlation >= topology_controls["minimum_neighbor_sensitivity_correlation"]
    )
    if not topology_controls["passed"]:
        raise RuntimeError(f"healthy topology controls failed: {topology_controls}")

    terminal_masks: dict[str, np.ndarray] = {}
    for lineage in LINEAGES:
        threshold = np.quantile(
            percentile_scores[lineage][healthy_mask], float(model_config["terminal_score_quantile"])
        )
        terminal_masks[lineage] = (
            healthy_mask
            & (lineage_matrix.argmax(axis=1) == LINEAGES.index(lineage))
            & (percentile_scores[lineage] >= threshold)
            & (pseudotime >= 0.5)
        )
        donors = cohort.observations.loc[terminal_masks[lineage], "donor"].nunique()
        if donors < int(model_config["minimum_state_donors"]):
            raise RuntimeError(f"terminal cells lack donor support for {lineage}")
    root_centroid = coordinates[root_mask & healthy_mask].mean(axis=0)
    terminal_centroids = np.vstack(
        [coordinates[terminal_masks[lineage]].mean(axis=0) for lineage in LINEAGES]
    )
    root_distances = np.linalg.norm(
        coordinates[root_mask & healthy_mask, None, :] - terminal_centroids[None, :, :], axis=2
    )
    projection_temperature = max(float(np.median(root_distances)), np.finfo(float).eps)

    root_definition = {
        "markers": present_markers["root"],
        "score_threshold": root_score_threshold,
        "maturity_threshold": root_maturity_threshold,
        "healthy_donor_count": root_donor_count,
        "reference_sha256": canonical_sha256(program_references["root"].round(10).tolist()),
    }
    terminal_definition = {
        lineage: {
            "markers": present_markers[lineage],
            "state": terminal_states[lineage],
            "healthy_donor_count": topology_controls["terminal_healthy_donors"][lineage],
            "reference_sha256": canonical_sha256(
                program_references[lineage].round(10).tolist()
            ),
        }
        for lineage in LINEAGES
    }
    hashes = _model_hashes(
        variable_indices,
        pca,
        state_model.cluster_centers_,
        absorption,
        root_definition,
        terminal_definition,
        model_config,
    )
    state_definitions = {
        "lineages": list(LINEAGES),
        "present_markers": present_markers,
        "root": root_definition,
        "terminals": terminal_definition,
        "states": states.to_dict(orient="records"),
        "graph": {str(state): targets for state, targets in graph.items()},
    }
    model = FateModel(
        genes=cohort.genes,
        variable_gene_indices=variable_indices,
        pca=pca,
        state_model=state_model,
        state_centroids=state_model.cluster_centers_.copy(),
        state_absorption=absorption,
        terminal_states=terminal_states,
        root_state=root_state,
        program_centers=centers,
        program_scales=scales,
        program_references=program_references,
        root_score_threshold=root_score_threshold,
        root_maturity_threshold=root_maturity_threshold,
        maturity_low=maturity_low,
        maturity_high=maturity_high,
        root_centroid=root_centroid,
        terminal_centroids=terminal_centroids,
        projection_temperature=projection_temperature,
        model_hashes=hashes,
        state_definitions=state_definitions,
        topology_controls=topology_controls,
    )
    result = map_fates(cohort, model, config)
    return model, result


def map_fates(cohort: CombinedCohort, model: FateModel, config: dict[str, Any]) -> FateResult:
    if not np.array_equal(cohort.genes, model.genes):
        raise ValueError("mapping gene order does not match frozen model")
    hashes_before = dict(model.model_hashes)
    expression = log_normalize(cohort.counts)
    raw_scores, _ = _raw_program_scores(expression, cohort.genes, config["lineage_markers"])
    scores = _standardize_programs(raw_scores, model.program_centers, model.program_scales)
    percentile_scores = _percentile_programs(scores, model.program_references)
    lineage_matrix = np.column_stack([percentile_scores[lineage] for lineage in LINEAGES])
    maturity = lineage_matrix.max(axis=1) - percentile_scores["root"]
    pseudotime = np.clip(
        (maturity - model.maturity_low) / (model.maturity_high - model.maturity_low), 0.0, 1.0
    )
    root_mask = (percentile_scores["root"] >= model.root_score_threshold) & (
        maturity <= model.root_maturity_threshold
    )
    dense = expression[:, model.variable_gene_indices].toarray().astype(np.float32, copy=False)
    coordinates = model.pca.transform(dense)
    del dense

    neighbor_count = int(config["model"]["mapping_neighbors"])
    neighbors = NearestNeighbors(n_neighbors=neighbor_count).fit(model.state_centroids)
    distances, indices = neighbors.kneighbors(coordinates)
    scale = np.maximum(np.median(distances, axis=1, keepdims=True), np.finfo(float).eps)
    weights = np.exp(-((distances / scale) ** 2))
    weights /= weights.sum(axis=1, keepdims=True)
    primary = np.einsum("nk,nkl->nl", weights, model.state_absorption[indices])

    terminal_distances = np.linalg.norm(
        coordinates[:, None, :] - model.terminal_centroids[None, :, :], axis=2
    )
    secondary = softmax(-terminal_distances / model.projection_temperature, axis=1)
    for probabilities in (primary, secondary):
        if not np.isfinite(probabilities).all() or np.any(probabilities < 0):
            raise RuntimeError("fate estimator returned invalid probabilities")
        if not np.allclose(probabilities.sum(axis=1), 1.0, atol=1e-6):
            raise RuntimeError("fate probabilities do not sum to one")

    hashes_after = _model_hashes(
        model.variable_gene_indices,
        model.pca,
        model.state_centroids,
        model.state_absorption,
        model.state_definitions["root"],
        model.state_definitions["terminals"],
        config["model"],
    )
    if hashes_before != hashes_after:
        raise RuntimeError("frozen healthy model changed during mapping")
    return FateResult(
        cohort.observations.copy(),
        root_mask,
        primary,
        secondary,
        percentile_scores,
        pseudotime,
        hashes_before,
        hashes_after,
    )
