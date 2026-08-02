import numpy as np
import pandas as pd
import pytest

from hspc_fate_mvp.model import FateResult, LINEAGES, _absorption_probabilities
from hspc_fate_mvp.statistics import (
    aggregate_donor_fates,
    classify_result,
    exact_permutation_tail,
    leave_one_out_effects,
)


def _fate_result() -> FateResult:
    observations = pd.DataFrame(
        {
            "cohort": ["discovery"] * 6,
            "donor": ["A", "A", "A", "A", "B", "B"],
            "condition": ["healthy"] * 4 + ["SCD"] * 2,
        }
    )
    primary = np.asarray(
        [
            [0.8, 0.1, 0.1],
            [0.7, 0.2, 0.1],
            [0.9, 0.05, 0.05],
            [0.0, 0.0, 1.0],
            [0.2, 0.3, 0.5],
            [0.0, 1.0, 0.0],
        ]
    )
    secondary = np.asarray(
        [
            [0.6, 0.3, 0.1],
            [0.5, 0.4, 0.1],
            [0.7, 0.2, 0.1],
            [0.2, 0.2, 0.6],
            [0.1, 0.2, 0.7],
            [0.1, 0.8, 0.1],
        ]
    )
    return FateResult(
        observations=observations,
        root_mask=np.asarray([True, True, True, False, True, False]),
        primary_probabilities=primary,
        secondary_probabilities=secondary,
        program_scores={},
        pseudotime=np.zeros(6),
        model_hashes_before={"model": "frozen"},
        model_hashes_after={"model": "frozen"},
    )


def test_donor_aggregation_uses_only_root_cells_and_preserves_donor_unit() -> None:
    frame = aggregate_donor_fates(_fate_result(), minimum_root_cells=2, enforce_minimum=False)
    donor_a = frame.loc[
        frame["donor"].eq("A") & frame["estimator"].eq("absorbing_markov")
    ].set_index("lineage")
    donor_b = frame.loc[frame["donor"].eq("B")]

    assert len(frame) == 2 * 2 * 3
    np.testing.assert_allclose(
        donor_a.loc[list(LINEAGES), "mean_fate_probability"], [0.8, 0.1166666667, 0.0833333333]
    )
    assert donor_a["root_cell_count"].eq(3).all()
    assert donor_a["root_cell_gate_passed"].all()
    assert donor_b["root_cell_count"].eq(1).all()
    assert not donor_b["root_cell_gate_passed"].any()
    assert frame.groupby(["donor", "estimator"])["mean_fate_probability"].sum().eq(1.0).all()

    with pytest.raises(RuntimeError, match="donor B has 1 HSC/MPP-like cells"):
        aggregate_donor_fates(_fate_result(), minimum_root_cells=2)


def test_exact_small_sample_controls_enumerate_all_assignments() -> None:
    values = np.asarray([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
    labels = np.asarray(["healthy"] * 3 + ["SCD"] * 3)
    tail, total = exact_permutation_tail(values, labels)
    leave_one_out = leave_one_out_effects(values, labels)

    assert total == 20
    assert tail == pytest.approx(0.1)
    assert np.all(leave_one_out > 0)


def test_absorbing_topology_reaches_every_terminal_with_normalized_probabilities() -> None:
    centroids = np.asarray(
        [[0.0, 0.0], [0.5, 0.4], [1.0, 0.0], [0.0, 1.0], [-1.0, 0.0]]
    )
    pseudotime = np.asarray([0.0, 0.4, 1.0, 1.0, 1.0])
    terminal_states = dict(zip(LINEAGES, [2, 3, 4], strict=True))
    absorption, graph = _absorption_probabilities(
        centroids, pseudotime, terminal_states, neighbors=3
    )

    np.testing.assert_allclose(absorption.sum(axis=1), 1.0)
    np.testing.assert_array_equal(absorption[[2, 3, 4]], np.eye(3))
    assert np.all(absorption[0] > 0)
    assert all(
        pseudotime[target] > pseudotime[source] or target in terminal_states.values()
        for source, edges in graph.items()
        for target, _ in edges
        if source not in terminal_states.values()
    )


def test_result_classification_distinguishes_null_specific_and_replicated() -> None:
    discovery = pd.DataFrame(
        {"lineage": list(LINEAGES), "stable": [False, False, False], "direction": [0, 0, 0]}
    )
    validation = discovery.copy()
    assert classify_result(discovery, validation) == (
        "no_stable_inferred_allocation_difference",
        [],
    )

    discovery.loc[0, ["stable", "direction"]] = [True, 1]
    assert classify_result(discovery, validation) == (
        "dataset_specific_inferred_allocation_difference",
        [],
    )

    validation.loc[0, ["stable", "direction"]] = [True, 1]
    assert classify_result(discovery, validation) == (
        "replicated_inferred_allocation_difference",
        ["myeloid"],
    )
