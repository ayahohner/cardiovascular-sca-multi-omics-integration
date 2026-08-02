from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd


COMPONENT_COLUMNS = [
    "reference_reproducibility",
    "directionality_evidence",
    "sca_signal_compatibility",
    "donor_study_coverage",
    "analysis_sensitivity",
    "transport_fitness",
]


def rank_hypotheses(
    edges: pd.DataFrame,
    compatibility: pd.DataFrame,
    config: dict[str, Any],
) -> pd.DataFrame:
    eligible = edges.loc[edges["eligible"]].merge(compatibility, on="edge_id", how="left").copy()
    if eligible.empty:
        return _empty_hypotheses()
    eligible["sca_signal_compatibility"] = eligible["sca_signal_compatibility"].fillna(0.0).clip(0.0, 1.0)
    eligible["reference_reproducibility"] = eligible.apply(_reference_reproducibility, axis=1)
    eligible["directionality_evidence"] = 0.0
    eligible["donor_study_coverage"] = eligible.apply(_coverage_score, axis=1)
    eligible["analysis_sensitivity"] = eligible["parameter_stability"].clip(0.0, 1.0)
    eligible["transport_fitness"] = float(config["transport_fitness"])
    weights = config["weights"]
    if set(weights) != set(COMPONENT_COLUMNS):
        raise ValueError("ranking weights must name every component exactly once")
    if not np.isclose(sum(float(value) for value in weights.values()), 1.0):
        raise ValueError("ranking weights must sum to 1")
    eligible["priority_score"] = sum(
        eligible[column] * float(weights[column]) for column in COMPONENT_COLUMNS
    )
    eligible = eligible.sort_values(
        ["priority_score", "reference_reproducibility", "edge_id"],
        ascending=[False, False, True],
    ).reset_index(drop=True)
    eligible.insert(0, "rank", np.arange(1, len(eligible) + 1))
    eligible["evidence_grade"] = "hypothesis_only"
    eligible["adjacency_confidence"] = np.where(
        (eligible["reference_reproducibility"] >= 0.8)
        & (eligible["donor_count_min"] >= 5)
        & ((eligible["study_count_min"] == 1) | (eligible["study_count_min"] >= 3)),
        "moderate",
        "low",
    )
    eligible["directionality_uncertainty"] = "unresolved; no eligible directional evidence in MVP"
    eligible["mapping_uncertainty"] = "indirect public HUVEC and bulk-PBMC transport"
    eligible["counterevidence"] = eligible.apply(_counterevidence, axis=1)
    eligible["permitted_interpretation"] = eligible.apply(
        lambda row: (
            f"{row.state_a} and {row.state_b} are reproducibly adjacent healthy {row.organ} "
            "endothelial states whose separating program is compatible with public SCA-relevant signals; "
            "no temporal direction or SCA-induced shift is inferred."
        ),
        axis=1,
    )
    eligible["next_data_requirement"] = eligible.apply(
        lambda row: (
            f"Independent adult {row.organ} endothelial data with donor-resolved SCA and comparator tissue, "
            "plus a time-resolved or direct perturbation design for this state pair."
        ),
        axis=1,
    )
    eligible["triggered_control_history"] = "none; all current eligibility gates passed"
    columns = [
        "rank",
        "edge_id",
        "organ",
        "state_a",
        "state_b",
        "priority_score",
        *COMPONENT_COLUMNS,
        "donor_count_min",
        "study_count_min",
        "donor_leaveout_stability",
        "study_leaveout_stability",
        "parameter_stability",
        "direction_status",
        "evidence_grade",
        "adjacency_confidence",
        "directionality_uncertainty",
        "mapping_uncertainty",
        "counterevidence",
        "permitted_interpretation",
        "next_data_requirement",
        "triggered_control_history",
    ]
    return eligible[columns]


def _reference_reproducibility(row: pd.Series) -> float:
    values = [float(row["donor_leaveout_stability"])]
    if pd.notna(row["study_leaveout_stability"]):
        values.append(float(row["study_leaveout_stability"]))
    return float(np.mean(values))


def _coverage_score(row: pd.Series) -> float:
    donor_score = min(1.0, float(row["donor_count_min"]) / 10.0)
    study_count = int(row["study_count_min"])
    study_score = min(1.0, study_count / 5.0) if study_count > 1 else 0.5
    return float((donor_score + study_score) / 2.0)


def _counterevidence(row: pd.Series) -> str:
    details = ["no eligible directional evidence"]
    if float(row["donor_leaveout_stability"]) < 0.8:
        details.append("donor leave-out stability below 0.8")
    if pd.notna(row["study_leaveout_stability"]) and float(row["study_leaveout_stability"]) < 0.8:
        details.append("study leave-out stability below 0.8")
    if float(row["parameter_stability"]) < 0.8:
        details.append("adjacency changes under some reasonable parameter choices")
    details.append("SCA evidence is transported from HUVEC and PBMC contexts")
    return "; ".join(details)


def _empty_hypotheses() -> pd.DataFrame:
    return pd.DataFrame(
        columns=[
            "rank",
            "edge_id",
            "organ",
            "state_a",
            "state_b",
            "priority_score",
            *COMPONENT_COLUMNS,
            "permitted_interpretation",
            "next_data_requirement",
        ]
    )
