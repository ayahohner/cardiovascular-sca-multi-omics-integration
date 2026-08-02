import numpy as np
import pandas as pd

from cardio_sca_mvp.ranking import COMPONENT_COLUMNS, rank_hypotheses
from cardio_sca_mvp.transitions import map_signals_to_edges


def edge_frame() -> pd.DataFrame:
    return pd.DataFrame(
        [
            {
                "edge_id": "heart::A__B",
                "organ": "heart",
                "state_a": "A",
                "state_b": "B",
                "donor_count_min": 8,
                "study_count_min": 1,
                "donor_leaveout_stability": 0.9,
                "study_leaveout_stability": np.nan,
                "parameter_stability": 0.8,
                "direction_status": "undirected",
                "eligible": True,
            },
            {
                "edge_id": "heart::B__C",
                "organ": "heart",
                "state_a": "B",
                "state_b": "C",
                "donor_count_min": 2,
                "study_count_min": 1,
                "donor_leaveout_stability": 0.2,
                "study_leaveout_stability": np.nan,
                "parameter_stability": 0.2,
                "direction_status": "undirected",
                "eligible": False,
            },
        ]
    )


def test_signal_alignment_is_compatibility_without_direction() -> None:
    genes = [f"G{i}" for i in range(12)]
    centroids = {
        "heart": pd.DataFrame(
            [np.zeros(12), np.ones(12), np.full(12, 2.0)],
            index=["A", "B", "C"],
            columns=genes,
        )
    }
    signals = pd.DataFrame(
        {
            "gene": genes * 2,
            "effect": [1.0] * 12 + [-1.0] * 12,
            "evidence_role": ["endothelial_perturbational_prior"] * 12 + ["immune_association"] * 12,
        }
    )
    mapping, compatibility = map_signals_to_edges(edge_frame(), centroids, signals)
    assert mapping["direction_status"].eq("undirected").all()
    np.testing.assert_allclose(mapping["absolute_alignment"].dropna(), 1.0, atol=1e-12)
    np.testing.assert_allclose(compatibility["sca_signal_compatibility"], 1.0, atol=1e-12)


def test_ranking_exposes_components_and_never_promotes_failed_edge() -> None:
    compatibility = pd.DataFrame(
        {"edge_id": ["heart::A__B", "heart::B__C"], "sca_signal_compatibility": [0.7, 1.0]}
    )
    config = {
        "weights": {
            "reference_reproducibility": 0.25,
            "directionality_evidence": 0.05,
            "sca_signal_compatibility": 0.25,
            "donor_study_coverage": 0.2,
            "analysis_sensitivity": 0.15,
            "transport_fitness": 0.1,
        },
        "transport_fitness": 0.25,
    }
    ranked = rank_hypotheses(edge_frame(), compatibility, config)
    assert ranked["edge_id"].tolist() == ["heart::A__B"]
    assert ranked["direction_status"].eq("undirected").all()
    assert ranked[COMPONENT_COLUMNS].notna().all().all()
    assert ((ranked[COMPONENT_COLUMNS] >= 0) & (ranked[COMPONENT_COLUMNS] <= 1)).all().all()
    assert "no temporal direction" in ranked["permitted_interpretation"].item()
