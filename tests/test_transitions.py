import numpy as np
import pandas as pd

from cardio_sca_mvp.reference import ReferenceData
from cardio_sca_mvp.transitions import generate_transitions, mutual_knn_edges


def test_mutual_nearest_neighbors_recover_known_state_pairs() -> None:
    states = pd.DataFrame(
        [[1.0, 0.0], [0.95, 0.05], [0.0, 1.0], [0.05, 0.95]],
        index=["A", "B", "C", "D"],
        columns=["g1", "g2"],
    )
    edges, distances = mutual_knn_edges(states, neighbors=1)
    assert edges == {("A", "B"), ("C", "D")}
    assert distances.loc["A", "B"] < distances.loc["A", "C"]


def synthetic_reference() -> ReferenceData:
    rng = np.random.default_rng(3)
    patterns = {
        "A": np.r_[np.repeat(5.0, 10), np.repeat(1.0, 10)],
        "B": np.r_[np.repeat(4.8, 10), np.repeat(1.2, 10)],
        "C": np.r_[np.repeat(1.0, 10), np.repeat(5.0, 10)],
        "D": np.r_[np.repeat(1.2, 10), np.repeat(4.8, 10)],
    }
    expression_rows = []
    groups = []
    for donor_index in range(6):
        for state, pattern in patterns.items():
            expression_rows.append(pattern + rng.normal(0, 0.02, len(pattern)))
            groups.append(
                {
                    "group_id": f"g-{donor_index}-{state}",
                    "organ": "heart",
                    "state": state,
                    "source_state": state,
                    "donor": f"D{donor_index}",
                    "study": "not_applicable",
                    "cell_count": 100,
                    "sample_count": 1,
                }
            )
    state_metadata = pd.DataFrame(
        {
            "organ": "heart",
            "state": list(patterns),
            "source_state": list(patterns),
            "cell_count": [600] * 4,
            "donor_count": [6] * 4,
            "sample_count": [6] * 4,
            "study_count": [1] * 4,
        }
    )
    return ReferenceData(
        organ="heart",
        expression=pd.DataFrame(expression_rows, columns=[f"G{i}" for i in range(20)]),
        groups=pd.DataFrame(groups),
        states=state_metadata,
        exclusions=pd.DataFrame(),
        gene_crosswalk=pd.DataFrame(),
    )


def test_transition_gate_requires_stability_and_bounded_distance() -> None:
    config = {
        "neighbors": 1,
        "max_variable_genes": 20,
        "min_gene_detection_fraction": 0.0,
        "excluded_gene_prefixes": [],
        "max_cosine_distance": 0.8,
        "min_donor_leaveout_stability": 0.8,
        "min_study_leaveout_stability": 0.5,
        "min_parameter_stability": 0.5,
    }
    result = generate_transitions([synthetic_reference()], config)
    assert set(map(tuple, result.edges[["state_a", "state_b"]].to_numpy())) == {("A", "B"), ("C", "D")}
    assert result.edges["eligible"].all()
    assert result.edges["donor_leaveout_stability"].min() == 1.0

    strict = {**config, "max_cosine_distance": 0.0}
    rejected = generate_transitions([synthetic_reference()], strict).edges
    assert not rejected["eligible"].any()
    assert rejected["failure_reasons"].str.contains("state_space_distance_exceeds_threshold").all()
