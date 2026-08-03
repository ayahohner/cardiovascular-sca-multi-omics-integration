from dataclasses import replace
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from hspc_fate_mvp.model import (
    FateResult,
    LINEAGES,
    _balanced_reference_indices,
    _cross_donor_ood_threshold,
    _map_reference_field,
)
from hspc_fate_mvp.statistics import (
    aggregate_donor_fates,
    classify_result,
    exact_permutation_tail,
    leave_one_out_effects,
)
from hspc_fate_mvp.robustness import attach_robustness


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
        reference_supported=np.ones(6, dtype=bool),
        ood_distance=np.zeros(6),
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
        frame["donor"].eq("A") & frame["estimator"].eq("palantir")
    ].set_index("lineage")
    donor_b = frame.loc[frame["donor"].eq("B")]

    assert len(frame) == 2 * 2 * 3
    np.testing.assert_allclose(
        donor_a.loc[list(LINEAGES), "mean_lineage_priming_probability"],
        [0.8, 0.1166666667, 0.0833333333],
    )
    assert donor_a["root_cell_count"].eq(3).all()
    assert donor_a["root_cell_gate_passed"].all()
    assert donor_b["root_cell_count"].eq(1).all()
    assert not donor_b["root_cell_gate_passed"].any()
    assert (
        frame.groupby(["donor", "estimator"])["mean_lineage_priming_probability"]
        .sum()
        .eq(1.0)
        .all()
    )

    with pytest.raises(RuntimeError, match="donor B has 1 HSC/MPP-like cells"):
        aggregate_donor_fates(_fate_result(), minimum_root_cells=2)


def test_donor_aggregation_rejects_an_unrepresented_donor_even_in_audit_mode() -> None:
    result = _fate_result()
    without_b = replace(result, root_mask=result.root_mask.copy())
    without_b.root_mask[without_b.observations["donor"].eq("B").to_numpy()] = False

    with pytest.raises(RuntimeError, match="donor B has no reference-supported"):
        aggregate_donor_fates(without_b, minimum_root_cells=2, enforce_minimum=False)


def test_exact_small_sample_controls_enumerate_all_assignments() -> None:
    values = np.asarray([0.0, 0.0, 0.0, 1.0, 1.0, 1.0])
    labels = np.asarray(["healthy"] * 3 + ["SCD"] * 3)
    tail, total = exact_permutation_tail(values, labels)
    leave_one_out = leave_one_out_effects(values, labels)

    assert total == 20
    assert tail == pytest.approx(0.1)
    assert np.all(leave_one_out > 0)


def test_reference_sampling_balances_donors_and_is_seed_reproducible() -> None:
    observations = pd.DataFrame(
        {
            "donor": ["A"] * 8 + ["B"] * 12 + ["C"] * 10,
            "condition": ["healthy"] * 30,
        }
    )
    healthy = np.ones(30, dtype=bool)
    first = _balanced_reference_indices(observations, healthy, 5, seed=11)
    second = _balanced_reference_indices(observations, healthy, 5, seed=11)

    np.testing.assert_array_equal(first, second)
    assert observations.iloc[first].groupby("donor").size().to_dict() == {
        "A": 5,
        "B": 5,
        "C": 5,
    }


def test_query_projection_normalizes_probabilities_and_exposes_far_ood_cells() -> None:
    coordinates = np.asarray(
        [[0.0, 0.0], [0.1, 0.0], [0.0, 0.1], [0.1, 0.1], [5.0, 5.0]]
    )
    donors = np.asarray(["A", "A", "B", "B"])
    threshold, calibration = _cross_donor_ood_threshold(
        coordinates[:4], donors, neighbors=2, quantile=0.99
    )
    model = SimpleNamespace(
        reference_coordinates=coordinates[:4],
        model_parameters={"mapping_neighbors": 3},
        primary_reference_probabilities=np.asarray(
            [[0.8, 0.1, 0.1], [0.7, 0.2, 0.1], [0.1, 0.8, 0.1], [0.1, 0.7, 0.2]]
        ),
        secondary_reference_probabilities=np.asarray(
            [[0.6, 0.2, 0.2], [0.5, 0.3, 0.2], [0.2, 0.6, 0.2], [0.2, 0.5, 0.3]]
        ),
        reference_pseudotime=np.asarray([0.1, 0.2, 0.7, 0.8]),
    )
    primary, secondary, pseudotime, distances = _map_reference_field(coordinates, model)

    assert calibration.shape == (4,)
    assert distances[-1] > threshold
    assert distances[0] <= threshold
    np.testing.assert_allclose(primary.sum(axis=1), 1.0)
    np.testing.assert_allclose(secondary.sum(axis=1), 1.0)
    assert pseudotime[0] < pseudotime[2]


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


def test_model_robustness_requires_both_estimators_and_all_bounded_refits() -> None:
    stability = pd.DataFrame({"lineage": ["myeloid"], "stable": [True]})
    rows = []
    for estimator in ("palantir", "cellrank_cytotrace"):
        rows.extend(
            [
                {
                    "lineage": "myeloid",
                    "estimator": estimator,
                    "perturbation": "baseline",
                    "direction": 1,
                    "status": "completed",
                },
                {
                    "lineage": "myeloid",
                    "estimator": estimator,
                    "perturbation": "donor_deletion_refit",
                    "direction": 1,
                    "status": "completed",
                },
                {
                    "lineage": "myeloid",
                    "estimator": estimator,
                    "perturbation": "hyperparameter_refit",
                    "direction": 1,
                    "status": "completed",
                },
            ]
        )
    robustness = pd.DataFrame(rows)
    config = {
        "model": {
            "discovery_leave_out_required": 1,
            "validation_leave_out_required": 1,
        }
    }

    passed = attach_robustness(stability, robustness, config, "discovery")
    assert bool(passed.loc[0, "model_robustness_passed"])
    assert bool(passed.loc[0, "stable"])

    failed_refit = robustness.copy()
    selector = (
        failed_refit["estimator"].eq("cellrank_cytotrace")
        & failed_refit["perturbation"].eq("hyperparameter_refit")
    )
    failed_refit.loc[selector, "direction"] = -1
    failed = attach_robustness(stability, failed_refit, config, "discovery")
    assert not bool(failed.loc[0, "model_robustness_passed"])
    assert not bool(failed.loc[0, "stable"])
