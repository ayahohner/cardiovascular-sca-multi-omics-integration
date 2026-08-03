from __future__ import annotations

import hashlib
import warnings
from dataclasses import dataclass
from typing import Any

import anndata as ad
import numpy as np
import palantir
import pandas as pd
import scanpy as sc
from cellrank.estimators import GPCCA
from cellrank.kernels import CytoTRACEKernel
from scipy import sparse
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors

from .matrixio import DonorMatrix
from .util import canonical_sha256


LINEAGES = ("myeloid", "erythroid_megakaryocytic", "lymphoid")
ESTIMATORS = ("palantir", "cellrank_cytotrace")


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
    reference_coordinates: np.ndarray
    reference_cell_ids: np.ndarray
    reference_donors: np.ndarray
    primary_reference_probabilities: np.ndarray
    secondary_reference_probabilities: np.ndarray
    reference_pseudotime: np.ndarray
    program_centers: dict[str, float]
    program_scales: dict[str, float]
    program_references: dict[str, np.ndarray]
    root_score_threshold: float
    root_pseudotime_threshold: float
    ood_distance_threshold: float
    model_parameters: dict[str, Any]
    model_hashes: dict[str, str]
    state_definitions: dict[str, Any]
    topology_controls: dict[str, Any]


@dataclass
class FateResult:
    observations: pd.DataFrame
    root_mask: np.ndarray
    reference_supported: np.ndarray
    ood_distance: np.ndarray
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
    scale = np.divide(
        10_000.0,
        totals,
        out=np.zeros_like(totals, dtype=np.float32),
        where=totals > 0,
    )
    normalized = sparse.diags(scale) @ normalized
    normalized.data = np.log1p(normalized.data)
    return normalized.tocsr()


def _gene_lookup(genes: np.ndarray) -> dict[str, int]:
    lookup = {str(gene).upper(): index for index, gene in enumerate(genes)}
    if len(lookup) != len(genes):
        raise RuntimeError("duplicate gene symbols must be collapsed before modeling")
    return lookup


def _raw_program_scores(
    expression: sparse.csr_matrix,
    genes: np.ndarray,
    marker_sets: dict[str, list[str]],
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
    raw_scores: dict[str, np.ndarray], reference_mask: np.ndarray
) -> tuple[dict[str, float], dict[str, float]]:
    centers: dict[str, float] = {}
    scales: dict[str, float] = {}
    for program, values in raw_scores.items():
        reference = values[reference_mask]
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
    standardized_scores: dict[str, np.ndarray], reference_mask: np.ndarray
) -> dict[str, np.ndarray]:
    return {
        program: np.sort(values[reference_mask].astype(float, copy=True))
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


def _balanced_reference_indices(
    observations: pd.DataFrame,
    healthy_mask: np.ndarray,
    cells_per_donor: int,
    seed: int,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    indices: list[np.ndarray] = []
    for donor in sorted(observations.loc[healthy_mask, "donor"].unique()):
        available = np.flatnonzero(healthy_mask & observations["donor"].eq(donor).to_numpy())
        if len(available) < cells_per_donor:
            raise RuntimeError(
                f"healthy donor {donor} has {len(available)} cells; {cells_per_donor} required"
            )
        indices.append(np.sort(rng.choice(available, size=cells_per_donor, replace=False)))
    return np.sort(np.concatenate(indices))


def _select_variable_genes(
    expression: sparse.csr_matrix,
    genes: np.ndarray,
    reference_indices: np.ndarray,
    count: int,
) -> np.ndarray:
    reference = expression[reference_indices]
    means = np.asarray(reference.mean(axis=0)).ravel()
    means_squared = np.asarray(reference.power(2).mean(axis=0)).ravel()
    variances = np.maximum(means_squared - means**2, 0)
    upper = np.char.upper(genes.astype(str))
    eligible = ~(
        np.char.startswith(upper, "MT-")
        | np.char.startswith(upper, "RPL")
        | np.char.startswith(upper, "RPS")
    )
    eligible_indices = np.flatnonzero(eligible & (means > 0))
    if len(eligible_indices) < count:
        raise RuntimeError("too few expressed genes for the configured variable-gene count")
    order = eligible_indices[np.argsort(variances[eligible_indices], kind="stable")]
    return np.sort(order[-count:])


def _cell_ids(observations: pd.DataFrame) -> np.ndarray:
    return (
        observations["cohort"].astype(str)
        + "|"
        + observations["donor"].astype(str)
        + "|"
        + observations["cell_index"].astype(str)
    ).to_numpy()


def _terminal_candidates(
    percentile_scores: dict[str, np.ndarray],
    reference_mask: np.ndarray,
    observations: pd.DataFrame,
    quantile: float,
    minimum_donors: int,
) -> dict[str, np.ndarray]:
    lineage_matrix = np.column_stack([percentile_scores[lineage] for lineage in LINEAGES])
    winners = lineage_matrix.argmax(axis=1)
    candidates: dict[str, np.ndarray] = {}
    for lineage_index, lineage in enumerate(LINEAGES):
        threshold = float(np.quantile(percentile_scores[lineage][reference_mask], quantile))
        other = np.max(np.delete(lineage_matrix, lineage_index, axis=1), axis=1)
        advantage = percentile_scores[lineage] - other
        eligible = np.flatnonzero(
            reference_mask
            & (winners == lineage_index)
            & (percentile_scores[lineage] >= threshold)
            & (advantage > 0)
        )
        donor_count = observations.iloc[eligible]["donor"].nunique()
        if donor_count < minimum_donors:
            raise RuntimeError(f"terminal candidates lack donor support for {lineage}")
        candidates[lineage] = eligible[np.argsort(advantage[eligible], kind="stable")[::-1]]
    return candidates


def _balanced_terminal_cells(
    candidates: np.ndarray,
    observations: pd.DataFrame,
    cell_ids: np.ndarray,
    count: int,
) -> list[str]:
    donors = observations.iloc[candidates]["donor"].to_numpy()
    unique_donors = sorted(set(donors))
    per_donor = max(1, int(np.ceil(count / len(unique_donors))))
    selected: list[int] = []
    for donor in unique_donors:
        selected.extend(candidates[donors == donor][:per_donor].tolist())
    selected_lookup = set(selected)
    ordered = [index for index in candidates if index in selected_lookup][:count]
    return cell_ids[ordered].tolist()


def _ordered_probabilities(
    values: Any, names: np.ndarray, source: str
) -> tuple[np.ndarray, float]:
    matrix = values.to_numpy() if isinstance(values, pd.DataFrame) else np.asarray(values)
    order = []
    for lineage in LINEAGES:
        matches = np.flatnonzero(np.asarray(names, dtype=str) == lineage)
        if len(matches) != 1:
            raise RuntimeError(f"{source} did not return exactly one {lineage} probability")
        order.append(int(matches[0]))
    ordered = np.asarray(matrix[:, order], dtype=float)
    if not np.isfinite(ordered).all() or np.any(ordered < 0):
        raise RuntimeError(f"{source} returned invalid probabilities")
    row_sums = ordered.sum(axis=1)
    maximum_deviation = float(np.max(np.abs(row_sums - 1.0)))
    if maximum_deviation > 0.02:
        raise RuntimeError(
            f"{source} probabilities do not sum to one: "
            f"minimum={row_sums.min():.8g}; maximum={row_sums.max():.8g}"
        )
    ordered /= row_sums[:, None]
    return ordered, maximum_deviation


def _cross_donor_ood_threshold(
    coordinates: np.ndarray,
    donors: np.ndarray,
    neighbors: int,
    quantile: float,
) -> tuple[float, np.ndarray]:
    distances: list[np.ndarray] = []
    for donor in sorted(set(donors)):
        query = donors == donor
        reference = ~query
        count = min(neighbors, int(reference.sum()))
        model = NearestNeighbors(n_neighbors=count).fit(coordinates[reference])
        donor_distances, _ = model.kneighbors(coordinates[query])
        distances.append(np.median(donor_distances, axis=1))
    calibration = np.concatenate(distances)
    return float(np.quantile(calibration, quantile)), calibration


def _array_sha256(values: np.ndarray) -> str:
    normalized = np.ascontiguousarray(np.asarray(values, dtype=np.float64).round(10))
    return hashlib.sha256(normalized.tobytes()).hexdigest()


def _model_hashes(model: FateModel) -> dict[str, str]:
    return {
        "transformation": canonical_sha256(
            {
                "variable_indices": model.variable_gene_indices.tolist(),
                "pca_mean": _array_sha256(model.pca.mean_),
                "pca_components": _array_sha256(model.pca.components_),
            }
        ),
        "trajectory": canonical_sha256(
            {
                "reference_cell_ids": model.reference_cell_ids.tolist(),
                "reference_coordinates": _array_sha256(model.reference_coordinates),
                "palantir_probabilities": _array_sha256(model.primary_reference_probabilities),
                "cellrank_probabilities": _array_sha256(model.secondary_reference_probabilities),
                "pseudotime": _array_sha256(model.reference_pseudotime),
            }
        ),
        "root_definition": canonical_sha256(model.state_definitions["root"]),
        "terminal_definitions": canonical_sha256(model.state_definitions["terminals"]),
        "model_parameters": canonical_sha256(model.model_parameters),
    }


def _map_reference_field(
    coordinates: np.ndarray,
    model: FateModel,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    neighbor_count = min(
        int(model.model_parameters["mapping_neighbors"]), len(model.reference_coordinates)
    )
    neighbors = NearestNeighbors(n_neighbors=neighbor_count).fit(model.reference_coordinates)
    distances, indices = neighbors.kneighbors(coordinates)
    scale = np.maximum(np.median(distances, axis=1, keepdims=True), np.finfo(float).eps)
    weights = np.exp(-((distances / scale) ** 2))
    weights /= weights.sum(axis=1, keepdims=True)
    primary = np.einsum(
        "nk,nkl->nl", weights, model.primary_reference_probabilities[indices]
    )
    secondary = np.einsum(
        "nk,nkl->nl", weights, model.secondary_reference_probabilities[indices]
    )
    pseudotime = np.einsum("nk,nk->n", weights, model.reference_pseudotime[indices])
    ood_distance = np.median(distances, axis=1)
    return primary, secondary, pseudotime, ood_distance


def fit_healthy_model(
    cohort: CombinedCohort, config: dict[str, Any]
) -> tuple[FateModel, FateResult]:
    healthy_mask = cohort.observations["condition"].eq("healthy").to_numpy()
    if not healthy_mask.any():
        raise RuntimeError("healthy reference has no cells")
    model_config = dict(config["model"])
    seed = int(config["random_seed"])
    reference_indices = _balanced_reference_indices(
        cohort.observations,
        healthy_mask,
        int(model_config["reference_cells_per_healthy_donor"]),
        seed,
    )
    reference_mask = np.zeros(len(healthy_mask), dtype=bool)
    reference_mask[reference_indices] = True
    expression = log_normalize(cohort.counts)

    construction_markers = config["marker_roles"]["construction"]
    raw_scores, present_construction = _raw_program_scores(
        expression, cohort.genes, construction_markers
    )
    centers, scales = _fit_program_standardization(raw_scores, reference_mask)
    standardized = _standardize_programs(raw_scores, centers, scales)
    references = _fit_program_references(standardized, reference_mask)
    scores = _percentile_programs(standardized, references)

    variable_indices = _select_variable_genes(
        expression,
        cohort.genes,
        reference_indices,
        int(model_config["highly_variable_genes"]),
    )
    dense = expression[:, variable_indices].toarray().astype(np.float32, copy=False)
    pca = PCA(
        n_components=int(model_config["pca_components"]),
        svd_solver="randomized",
        random_state=seed,
    )
    pca.fit(dense[reference_indices])
    coordinates = pca.transform(dense)
    del dense

    cell_ids = _cell_ids(cohort.observations)
    reference_cell_ids = cell_ids[reference_indices]
    root_advantage = scores["root"] - np.column_stack(
        [scores[lineage] for lineage in LINEAGES]
    ).max(axis=1)
    root_seed_index = int(reference_indices[np.argmax(root_advantage[reference_indices])])
    terminal_candidates = _terminal_candidates(
        scores,
        reference_mask,
        cohort.observations,
        float(model_config["terminal_score_quantile"]),
        int(model_config["minimum_state_donors"]),
    )
    terminal_seed_indices = {
        lineage: int(candidates[0]) for lineage, candidates in terminal_candidates.items()
    }

    reference_expression = expression[reference_indices].copy()
    trajectory = ad.AnnData(reference_expression)
    trajectory.obs_names = reference_cell_ids
    trajectory.var_names = cohort.genes
    trajectory.obs["donor"] = cohort.observations.iloc[reference_indices]["donor"].to_numpy()
    trajectory.obsm["X_pca"] = coordinates[reference_indices]
    palantir.utils.run_diffusion_maps(
        trajectory,
        n_components=int(model_config["diffusion_components"]),
        knn=int(model_config["trajectory_neighbors"]),
        seed=seed,
        pca_key="X_pca",
    )
    palantir.utils.determine_multiscale_space(trajectory)
    terminal_seed_series = pd.Series(
        {
            cell_ids[index]: lineage
            for lineage, index in terminal_seed_indices.items()
        }
    )
    with warnings.catch_warnings(record=True) as trajectory_warnings:
        warnings.simplefilter("always", UserWarning)
        palantir.core.run_palantir(
            trajectory,
            early_cell=cell_ids[root_seed_index],
            terminal_states=terminal_seed_series,
            knn=int(model_config["trajectory_neighbors"]),
            num_waypoints=min(int(model_config["palantir_waypoints"]), len(reference_indices)),
            n_jobs=1,
            use_early_cell_as_start=True,
            seed=seed,
        )
    unreachable_warning_count = sum(
        "unreachable" in str(item.message).lower() for item in trajectory_warnings
    )
    palantir_values = trajectory.obsm["palantir_fate_probabilities"]
    if isinstance(palantir_values, pd.DataFrame):
        palantir_names = palantir_values.columns.to_numpy(dtype=str)
    else:
        palantir_names = np.asarray(
            trajectory.uns["palantir_fate_probabilities_columns"], dtype=str
        )
    primary_reference, palantir_normalization_deviation = _ordered_probabilities(
        palantir_values, palantir_names, "Palantir"
    )
    reference_pseudotime = trajectory.obs["palantir_pseudotime"].to_numpy(dtype=float)

    sc.pp.neighbors(
        trajectory,
        n_neighbors=int(model_config["trajectory_neighbors"]),
        use_rep="X_pca",
        random_state=seed,
    )
    cytotrace = CytoTRACEKernel(trajectory)
    cytotrace.compute_cytotrace(
        layer="X", n_genes=int(model_config["cellrank_cytotrace_genes"])
    )
    cytotrace.compute_transition_matrix(
        threshold_scheme="soft",
        frac_to_keep=0.3,
        n_jobs=1,
        backend="threading",
        show_progress_bar=False,
    )
    terminal_cells = {
        lineage: _balanced_terminal_cells(
            candidates,
            cohort.observations,
            cell_ids,
            int(model_config["cellrank_terminal_cells_per_lineage"]),
        )
        for lineage, candidates in terminal_candidates.items()
    }
    estimator = GPCCA(cytotrace)
    estimator.set_terminal_states(terminal_cells)
    estimator.compute_fate_probabilities(
        use_petsc=False,
        n_jobs=1,
        backend="threading",
        show_progress_bar=False,
    )
    secondary_reference, cellrank_normalization_deviation = _ordered_probabilities(
        estimator.fate_probabilities.X,
        np.asarray(estimator.fate_probabilities.names, dtype=str),
        "CellRank CytoTRACE",
    )

    root_pseudotime_threshold = float(
        np.quantile(reference_pseudotime, float(model_config["root_pseudotime_quantile"]))
    )
    root_reference = scores["root"][reference_indices] >= float(
        model_config["root_score_quantile"]
    )
    root_reference_indices = reference_indices[root_reference]
    root_donor_count = int(cohort.observations.iloc[root_reference_indices]["donor"].nunique())

    validation_raw, present_validation = _raw_program_scores(
        expression, cohort.genes, config["marker_roles"]["validation"]
    )
    validation_centers, validation_scales = _fit_program_standardization(
        validation_raw, reference_mask
    )
    validation_standardized = _standardize_programs(
        validation_raw, validation_centers, validation_scales
    )
    validation_references = _fit_program_references(
        validation_standardized, reference_mask
    )
    validation_scores = _percentile_programs(
        validation_standardized, validation_references
    )
    marker_order = {
        lineage: bool(
            np.median(validation_scores[lineage][terminal_candidates[lineage]])
            > np.median(validation_scores[lineage][root_reference_indices])
        )
        for lineage in LINEAGES
    }
    terminal_union = np.unique(np.concatenate(list(terminal_candidates.values())))
    root_marker_order = bool(
        np.median(validation_scores["root"][root_reference_indices])
        > np.median(validation_scores["root"][terminal_union])
    )

    reference_donors = cohort.observations.iloc[reference_indices]["donor"].to_numpy()
    ood_threshold, ood_calibration = _cross_donor_ood_threshold(
        coordinates[reference_indices],
        reference_donors,
        int(model_config["mapping_neighbors"]),
        float(model_config["ood_quantile"]),
    )
    root_to_all = bool(
        root_reference.any()
        and np.all(primary_reference[root_reference].mean(axis=0) > 0)
        and np.all(secondary_reference[root_reference].mean(axis=0) > 0)
    )
    terminal_donors = {
        lineage: int(cohort.observations.iloc[candidates]["donor"].nunique())
        for lineage, candidates in terminal_candidates.items()
    }
    estimator_correlations = {
        lineage: float(
            np.corrcoef(primary_reference[:, index], secondary_reference[:, index])[0, 1]
        )
        for index, lineage in enumerate(LINEAGES)
    }
    topology_controls = {
        "method": "palantir_primary_cellrank_cytotrace_sensitivity",
        "reference_cells": int(len(reference_indices)),
        "reference_cells_per_healthy_donor": int(
            model_config["reference_cells_per_healthy_donor"]
        ),
        "root_to_all_terminals": root_to_all,
        "held_out_marker_order": marker_order,
        "held_out_root_marker_order": root_marker_order,
        "root_healthy_donors": root_donor_count,
        "terminal_healthy_donors": terminal_donors,
        "estimator_reference_correlations": estimator_correlations,
        "maximum_probability_normalization_deviation": {
            "palantir": palantir_normalization_deviation,
            "cellrank_cytotrace": cellrank_normalization_deviation,
        },
        "palantir_unreachable_warning_count": unreachable_warning_count,
        "ood_distance_threshold": ood_threshold,
        "ood_calibration_median": float(np.median(ood_calibration)),
        "ood_calibration_99th_percentile": float(np.quantile(ood_calibration, 0.99)),
    }
    topology_controls["passed"] = bool(
        root_to_all
        and all(marker_order.values())
        and root_marker_order
        and root_donor_count >= int(model_config["minimum_state_donors"])
        and all(
            donor_count >= int(model_config["minimum_state_donors"])
            for donor_count in terminal_donors.values()
        )
    )
    state_definitions = {
        "estimand": "healthy-reference lineage priming",
        "construction_markers": present_construction,
        "validation_markers": present_validation,
        "root": {
            "seed_cell": cell_ids[root_seed_index],
            "score_threshold": float(model_config["root_score_quantile"]),
            "pseudotime_threshold_diagnostic_only": root_pseudotime_threshold,
            "healthy_donor_count": root_donor_count,
        },
        "terminals": {
            lineage: {
                "palantir_seed_cell": cell_ids[terminal_seed_indices[lineage]],
                "cellrank_terminal_cells": terminal_cells[lineage],
                "healthy_donor_count": terminal_donors[lineage],
            }
            for lineage in LINEAGES
        },
        "reference_cell_count": int(len(reference_indices)),
        "out_of_distribution": {
            "metric": "median Euclidean distance to frozen healthy PCA neighbors",
            "calibration": "leave-one-healthy-donor-out",
            "threshold_quantile": float(model_config["ood_quantile"]),
            "threshold": ood_threshold,
        },
    }
    model = FateModel(
        genes=cohort.genes,
        variable_gene_indices=variable_indices,
        pca=pca,
        reference_coordinates=coordinates[reference_indices],
        reference_cell_ids=reference_cell_ids,
        reference_donors=reference_donors,
        primary_reference_probabilities=primary_reference,
        secondary_reference_probabilities=secondary_reference,
        reference_pseudotime=reference_pseudotime,
        program_centers=centers,
        program_scales=scales,
        program_references=references,
        root_score_threshold=float(model_config["root_score_quantile"]),
        root_pseudotime_threshold=root_pseudotime_threshold,
        ood_distance_threshold=ood_threshold,
        model_parameters=model_config,
        model_hashes={},
        state_definitions=state_definitions,
        topology_controls=topology_controls,
    )
    model.model_hashes = _model_hashes(model)
    return model, map_fates(cohort, model, config)


def map_fates(cohort: CombinedCohort, model: FateModel, config: dict[str, Any]) -> FateResult:
    if not np.array_equal(cohort.genes, model.genes):
        raise ValueError("mapping gene order does not match frozen model")
    hashes_before = dict(model.model_hashes)
    expression = log_normalize(cohort.counts)
    raw_scores, _ = _raw_program_scores(
        expression, cohort.genes, config["marker_roles"]["construction"]
    )
    standardized = _standardize_programs(
        raw_scores, model.program_centers, model.program_scales
    )
    percentile_scores = _percentile_programs(standardized, model.program_references)
    dense = expression[:, model.variable_gene_indices].toarray().astype(np.float32, copy=False)
    coordinates = model.pca.transform(dense)
    del dense
    primary, secondary, pseudotime, ood_distance = _map_reference_field(coordinates, model)
    reference_supported = ood_distance <= model.ood_distance_threshold
    root_mask = (
        (percentile_scores["root"] >= model.root_score_threshold)
        & reference_supported
    )
    for source, probabilities in zip(ESTIMATORS, (primary, secondary), strict=True):
        if not np.isfinite(probabilities).all() or np.any(probabilities < 0):
            raise RuntimeError(f"{source} mapping returned invalid probabilities")
        if not np.allclose(probabilities.sum(axis=1), 1.0, atol=1e-5):
            raise RuntimeError(f"{source} mapping probabilities do not sum to one")
    hashes_after = _model_hashes(model)
    if hashes_before != hashes_after:
        raise RuntimeError("frozen healthy trajectory changed during query mapping")
    observations = cohort.observations.copy()
    observations["reference_supported"] = reference_supported
    observations["ood_distance"] = ood_distance
    return FateResult(
        observations=observations,
        root_mask=root_mask,
        reference_supported=reference_supported,
        ood_distance=ood_distance,
        primary_probabilities=primary,
        secondary_probabilities=secondary,
        program_scores=percentile_scores,
        pseudotime=pseudotime,
        model_hashes_before=hashes_before,
        model_hashes_after=hashes_after,
    )
