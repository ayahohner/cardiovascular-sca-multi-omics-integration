import numpy as np
import pandas as pd

from cardio_sca_mvp.ranking import COMPONENT_COLUMNS
from cardio_sca_mvp.validation import validate_analysis, validate_local_markdown_links


def valid_inputs():
    preflight = pd.DataFrame({"status": ["admitted", "admitted"]})
    states = pd.DataFrame(
        {"organ": ["heart", "lung"], "donor_count": [5, 5]}
    )
    edges = pd.DataFrame(
        {
            "edge_id": ["heart::A__B"],
            "eligible": [True],
            "donor_leaveout_stability": [0.9],
            "study_leaveout_stability": [0.8],
            "direction_status": ["undirected"],
        }
    )
    signals = pd.DataFrame(
        {
            "evidence_role": [
                "endothelial_perturbational_prior",
                "immune_association",
            ]
        }
    )
    mapping = pd.DataFrame(
        {
            "direction_status": ["undirected"],
            "control_response": ["mapped"],
            "overlap_gene_count": [20],
        }
    )
    hypotheses = pd.DataFrame(
        {
            "edge_id": ["heart::A__B"],
            **{column: [0.5] for column in COMPONENT_COLUMNS},
            "permitted_interpretation": [
                "Healthy states are adjacent; no temporal direction or disease shift is inferred."
            ],
        }
    )
    config = {
        "pseudobulk": {"min_donors_per_state": 3},
        "adjacency": {
            "min_donor_leaveout_stability": 0.6,
            "min_study_leaveout_stability": 0.5,
        },
    }
    return preflight, states, edges, signals, mapping, hypotheses, config


def test_validation_checks_scientific_not_just_file_shape() -> None:
    inputs = valid_inputs()
    report = validate_analysis(*inputs)
    assert report["status"] == "passed"
    assert report["failed_count"] == 0

    preflight, states, edges, signals, mapping, hypotheses, config = valid_inputs()
    edges.loc[0, "donor_leaveout_stability"] = 0.2
    edges.loc[0, "direction_status"] = "directed"
    mapping.loc[0, "overlap_gene_count"] = 3
    hypotheses.loc[0, "reference_reproducibility"] = 1.2
    hypotheses.loc[0, "permitted_interpretation"] = "SCA causes a transition and predicts patient prognosis."
    failed = validate_analysis(
        preflight, states, edges, signals, mapping, hypotheses, config
    )
    failed_names = {check["name"] for check in failed["checks"] if not check["passed"]}
    assert {
        "donor_stability_gate",
        "all_edges_undirected",
        "mapped_signals_have_gene_overlap",
        "component_scores_visible_and_bounded",
        "claim_boundaries",
    } <= failed_names


def test_local_markdown_link_check_detects_broken_release_navigation(tmp_path) -> None:
    (tmp_path / "README.md").write_text(
        "[missing](docs/result.md) [external](https://example.org)", encoding="utf-8"
    )
    failed = validate_local_markdown_links(tmp_path)
    assert not failed["passed"]
    (tmp_path / "docs").mkdir()
    (tmp_path / "docs/result.md").write_text("result", encoding="utf-8")
    assert validate_local_markdown_links(tmp_path)["passed"]
