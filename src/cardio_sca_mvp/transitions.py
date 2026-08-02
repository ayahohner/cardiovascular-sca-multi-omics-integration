from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import pairwise_distances

from .reference import ReferenceData


@dataclass
class TransitionResult:
    edges: pd.DataFrame
    state_centroids: dict[str, pd.DataFrame]
    state_programs: pd.DataFrame
    diagnostics: pd.DataFrame


def generate_transitions(
    references: list[ReferenceData],
    config: dict[str, Any],
) -> TransitionResult:
    edge_frames: list[pd.DataFrame] = []
    program_frames: list[pd.DataFrame] = []
    diagnostic_frames: list[pd.DataFrame] = []
    centroids_by_organ: dict[str, pd.DataFrame] = {}
    for reference in references:
        edges, centroids, programs, diagnostics = _organ_transitions(reference, config)
        edge_frames.append(edges)
        program_frames.append(programs)
        diagnostic_frames.append(diagnostics)
        centroids_by_organ[reference.organ] = centroids
    return TransitionResult(
        edges=pd.concat(edge_frames, ignore_index=True),
        state_centroids=centroids_by_organ,
        state_programs=pd.concat(program_frames, ignore_index=True),
        diagnostics=pd.concat(diagnostic_frames, ignore_index=True),
    )


def _organ_transitions(
    reference: ReferenceData,
    config: dict[str, Any],
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    expression = reference.expression
    groups = reference.groups.reset_index(drop=True)
    if len(expression) != len(groups):
        raise ValueError(f"{reference.organ} expression and group metadata are misaligned")
    centroids = _state_centroids(expression, groups)
    feature_sets = _feature_sets(expression, centroids, config)
    main_features = feature_sets[max(feature_sets)]
    standardized = _standardize_state_matrix(centroids[main_features])
    neighbor_count = min(int(config["neighbors"]), len(standardized) - 1)
    base_edges, distances = mutual_knn_edges(standardized, neighbor_count)
    if not base_edges:
        raise ValueError(f"{reference.organ} produced no mutual-nearest-neighbor edges")

    donor_graphs = _leaveout_graphs(
        expression,
        groups,
        main_features,
        "donor",
        neighbor_count,
        required_states=set(centroids.index),
    )
    study_values = groups["study"].astype(str).unique()
    if len(study_values) > 1 and set(study_values) != {"not_applicable"}:
        study_graphs = _leaveout_graphs(
            expression,
            groups,
            main_features,
            "study",
            neighbor_count,
            required_states=set(centroids.index),
        )
    else:
        study_graphs = []

    parameter_graphs: list[set[tuple[str, str]]] = []
    candidate_neighbors = sorted(
        {
            max(1, neighbor_count - 1),
            neighbor_count,
            min(len(centroids) - 1, neighbor_count + 1),
        }
    )
    for feature_count, k in product(sorted(feature_sets), candidate_neighbors):
        matrix = _standardize_state_matrix(centroids[feature_sets[feature_count]])
        graph, _ = mutual_knn_edges(matrix, k)
        parameter_graphs.append(graph)

    state_metadata = reference.states.set_index("state")
    rows: list[dict[str, Any]] = []
    for state_a, state_b in sorted(base_edges):
        distance = float(distances.loc[state_a, state_b])
        donor_stability = _presence_fraction((state_a, state_b), donor_graphs)
        study_stability = (
            _presence_fraction((state_a, state_b), study_graphs)
            if study_graphs
            else np.nan
        )
        parameter_stability = _presence_fraction((state_a, state_b), parameter_graphs)
        donor_count = int(
            min(
                state_metadata.loc[state_a, "donor_count"],
                state_metadata.loc[state_b, "donor_count"],
            )
        )
        study_count = int(
            min(
                state_metadata.loc[state_a, "study_count"],
                state_metadata.loc[state_b, "study_count"],
            )
        )
        study_pass = (
            True
            if np.isnan(study_stability)
            else study_stability >= float(config["min_study_leaveout_stability"])
        )
        eligible = (
            donor_count >= 3
            and distance <= float(config["max_cosine_distance"])
            and donor_stability >= float(config["min_donor_leaveout_stability"])
            and parameter_stability >= float(config.get("min_parameter_stability", 0.5))
            and study_pass
        )
        failure_reasons: list[str] = []
        if donor_count < 3:
            failure_reasons.append("insufficient_donor_support")
        if distance > float(config["max_cosine_distance"]):
            failure_reasons.append("state_space_distance_exceeds_threshold")
        if donor_stability < float(config["min_donor_leaveout_stability"]):
            failure_reasons.append("donor_leaveout_instability")
        if parameter_stability < float(config.get("min_parameter_stability", 0.5)):
            failure_reasons.append("parameter_instability")
        if not study_pass:
            failure_reasons.append("study_leaveout_instability")
        rows.append(
            {
                "edge_id": f"{reference.organ}::{state_a}__{state_b}",
                "organ": reference.organ,
                "state_a": state_a,
                "state_b": state_b,
                "lineage": "endothelial",
                "adjacency_method": f"mutual_{neighbor_count}_nearest_neighbors",
                "distance": distance,
                "donor_count_min": donor_count,
                "study_count_min": study_count,
                "donor_leaveout_stability": donor_stability,
                "study_leaveout_stability": study_stability,
                "parameter_stability": parameter_stability,
                "direction_status": "undirected",
                "direction_evidence": "absent",
                "eligible": bool(eligible),
                "control_response": "retained" if eligible else "excluded",
                "failure_reasons": ";".join(failure_reasons) if failure_reasons else "none",
            }
        )
    edges = pd.DataFrame(rows)
    programs = _state_programs(reference.organ, standardized, maximum=20)
    diagnostics = pd.DataFrame(
        [
            {
                "organ": reference.organ,
                "group_count": len(groups),
                "state_count": len(centroids),
                "gene_count": expression.shape[1],
                "selected_feature_count": len(main_features),
                "base_edge_count": len(base_edges),
                "eligible_edge_count": int(edges["eligible"].sum()),
                "donor_leaveout_runs": len(donor_graphs),
                "study_leaveout_runs": len(study_graphs),
                "parameter_runs": len(parameter_graphs),
            }
        ]
    )
    return edges, centroids, programs, diagnostics


def _state_centroids(expression: pd.DataFrame, groups: pd.DataFrame) -> pd.DataFrame:
    frame = expression.copy()
    frame.index = groups["state"].to_numpy()
    return frame.groupby(level=0, observed=True).mean().sort_index()


def _feature_sets(
    expression: pd.DataFrame,
    centroids: pd.DataFrame,
    config: dict[str, Any],
) -> dict[int, list[str]]:
    detected = (np.expm1(expression) >= 1.0).mean(axis=0)
    eligible = detected[detected >= float(config["min_gene_detection_fraction"])].index
    prefixes = tuple(str(value).upper() for value in config["excluded_gene_prefixes"])
    eligible = [gene for gene in eligible if not str(gene).upper().startswith(prefixes)]
    variance = centroids[eligible].var(axis=0, ddof=0).sort_values(ascending=False)
    maximum = min(int(config["max_variable_genes"]), len(variance))
    if maximum < 10:
        raise ValueError("too few variable genes for adjacency")
    counts = sorted({maximum, max(10, maximum // 2)})
    return {count: variance.head(count).index.tolist() for count in counts}


def _standardize_state_matrix(matrix: pd.DataFrame) -> pd.DataFrame:
    centered = matrix - matrix.mean(axis=0)
    scale = matrix.std(axis=0, ddof=0)
    nonzero = scale > 0
    return centered.loc[:, nonzero] / scale.loc[nonzero]


def mutual_knn_edges(
    state_matrix: pd.DataFrame,
    neighbors: int,
) -> tuple[set[tuple[str, str]], pd.DataFrame]:
    states = list(state_matrix.index)
    if len(states) < 2:
        return set(), pd.DataFrame(index=states, columns=states, dtype=float)
    k = min(max(1, int(neighbors)), len(states) - 1)
    distances = pairwise_distances(state_matrix.to_numpy(), metric="cosine")
    distance_frame = pd.DataFrame(distances, index=states, columns=states)
    nearest: dict[str, set[str]] = {}
    for index, state in enumerate(states):
        order = np.argsort(distances[index])
        nearest[state] = set([states[position] for position in order if position != index][:k])
    edges: set[tuple[str, str]] = set()
    for state_a in states:
        for state_b in nearest[state_a]:
            if state_a in nearest[state_b]:
                edges.add(tuple(sorted((state_a, state_b))))
    return edges, distance_frame


def _leaveout_graphs(
    expression: pd.DataFrame,
    groups: pd.DataFrame,
    features: list[str],
    column: str,
    neighbors: int,
    required_states: set[str],
) -> list[set[tuple[str, str]]]:
    graphs: list[set[tuple[str, str]]] = []
    for value in sorted(groups[column].astype(str).unique()):
        keep = groups[column].astype(str).ne(value).to_numpy()
        subset_groups = groups.loc[keep].reset_index(drop=True)
        if set(subset_groups["state"]) != required_states:
            continue
        subset_expression = expression.loc[keep, features].reset_index(drop=True)
        centroids = _state_centroids(subset_expression, subset_groups)
        graph, _ = mutual_knn_edges(_standardize_state_matrix(centroids), neighbors)
        graphs.append(graph)
    return graphs


def _presence_fraction(
    edge: tuple[str, str],
    graphs: list[set[tuple[str, str]]],
) -> float:
    if not graphs:
        return float("nan")
    canonical = tuple(sorted(edge))
    return float(np.mean([canonical in graph for graph in graphs]))


def _state_programs(
    organ: str,
    standardized: pd.DataFrame,
    maximum: int,
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for state, row in standardized.iterrows():
        for rank, (gene, score) in enumerate(row.nlargest(maximum).items(), start=1):
            records.append(
                {
                    "organ": organ,
                    "state": state,
                    "rank": rank,
                    "gene": gene,
                    "state_specificity_z": float(score),
                }
            )
    return pd.DataFrame(records)


def map_signals_to_edges(
    edges: pd.DataFrame,
    centroids_by_organ: dict[str, pd.DataFrame],
    signals: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    records: list[dict[str, Any]] = []
    roles = sorted(signals["evidence_role"].unique())
    for edge in edges.itertuples(index=False):
        centroids = centroids_by_organ[edge.organ]
        contrast = centroids.loc[edge.state_b] - centroids.loc[edge.state_a]
        for role in roles:
            signature_frame = signals.loc[signals["evidence_role"].eq(role), ["gene", "effect"]]
            signature = signature_frame.groupby("gene", observed=True)["effect"].mean()
            genes = sorted(set(contrast.index) & set(signature.index))
            if len(genes) < 10:
                signed_alignment = np.nan
                absolute_alignment = np.nan
                orientation = "unresolved"
                response = "abstained"
            else:
                x = contrast.loc[genes].to_numpy(dtype=float)
                y = signature.loc[genes].to_numpy(dtype=float)
                denominator = np.linalg.norm(x) * np.linalg.norm(y)
                signed_alignment = float(np.dot(x, y) / denominator) if denominator else np.nan
                absolute_alignment = abs(signed_alignment)
                orientation = "toward_state_b" if signed_alignment >= 0 else "toward_state_a"
                response = "mapped"
            records.append(
                {
                    "edge_id": edge.edge_id,
                    "organ": edge.organ,
                    "state_a": edge.state_a,
                    "state_b": edge.state_b,
                    "evidence_role": role,
                    "overlap_gene_count": len(genes),
                    "signed_alignment_to_arbitrary_a_to_b_contrast": signed_alignment,
                    "absolute_alignment": absolute_alignment,
                    "contrast_orientation": orientation,
                    "direction_status": "undirected",
                    "permitted_interpretation": "molecular compatibility, not temporal direction or an SCA-induced shift",
                    "control_response": response,
                }
            )
    mapping = pd.DataFrame(records)
    compatibility = (
        mapping.groupby("edge_id", observed=True)["absolute_alignment"]
        .mean()
        .rename("sca_signal_compatibility")
        .reset_index()
    )
    return mapping, compatibility
