from __future__ import annotations

import itertools
import json
from typing import Any

import numpy as np
import pandas as pd

from .model import FateResult, LINEAGES


ESTIMATORS = ("absorbing_markov", "terminal_projection")


def aggregate_donor_fates(
    result: FateResult, minimum_root_cells: int, enforce_minimum: bool = True
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    donor_values = result.observations["donor"].to_numpy()
    condition_values = result.observations["condition"].to_numpy()
    cohort_values = result.observations["cohort"].to_numpy()
    probability_sets = {
        "absorbing_markov": result.primary_probabilities,
        "terminal_projection": result.secondary_probabilities,
    }
    for donor in pd.unique(donor_values):
        donor_mask = donor_values == donor
        root_mask = donor_mask & result.root_mask
        root_count = int(root_mask.sum())
        if enforce_minimum and root_count < minimum_root_cells:
            raise RuntimeError(
                f"donor {donor} has {root_count} HSC/MPP-like cells; minimum is {minimum_root_cells}"
            )
        condition = str(pd.unique(condition_values[donor_mask]).item())
        cohort = str(pd.unique(cohort_values[donor_mask]).item())
        for estimator, probabilities in probability_sets.items():
            means = probabilities[root_mask].mean(axis=0)
            for lineage, value in zip(LINEAGES, means, strict=True):
                rows.append(
                    {
                        "cohort": cohort,
                        "donor": donor,
                        "condition": condition,
                        "estimator": estimator,
                        "lineage": lineage,
                        "root_cell_count": root_count,
                        "root_cell_gate_passed": root_count >= minimum_root_cells,
                        "mean_fate_probability": float(value),
                    }
                )
    frame = pd.DataFrame.from_records(rows)
    expected = len(pd.unique(donor_values)) * len(ESTIMATORS) * len(LINEAGES)
    if len(frame) != expected:
        raise RuntimeError("donor aggregation did not produce one row per donor, estimator, and lineage")
    return frame


def _effect(values: np.ndarray, labels: np.ndarray) -> float:
    return float(values[labels == "SCD"].mean() - values[labels == "healthy"].mean())


def exact_permutation_tail(values: np.ndarray, labels: np.ndarray) -> tuple[float, int]:
    distribution = _exact_permutation_distribution(
        values, labels, np.asarray([str(index) for index in range(len(values))])
    )
    observed = abs(_effect(values, labels))
    tails = sum(abs(item["effect"]) >= observed - 1e-12 for item in distribution)
    return tails / len(distribution), len(distribution)


def _exact_permutation_distribution(
    values: np.ndarray, labels: np.ndarray, donors: np.ndarray
) -> list[dict[str, Any]]:
    scd_count = int((labels == "SCD").sum())
    distribution: list[dict[str, Any]] = []
    for scd_indices in itertools.combinations(range(len(values)), scd_count):
        permuted = np.full(len(values), "healthy", dtype=object)
        permuted[list(scd_indices)] = "SCD"
        distribution.append(
            {
                "scd_donors": [str(donors[index]) for index in scd_indices],
                "effect": _effect(values, permuted),
            }
        )
    return distribution


def leave_one_out_effects(values: np.ndarray, labels: np.ndarray) -> np.ndarray:
    effects = []
    for excluded in range(len(values)):
        keep = np.arange(len(values)) != excluded
        effects.append(_effect(values[keep], labels[keep]))
    return np.asarray(effects, dtype=float)


def bootstrap_interval(
    values: np.ndarray, labels: np.ndarray, seed: int, replicates: int = 2000
) -> tuple[float, float]:
    generator = np.random.default_rng(seed)
    healthy = values[labels == "healthy"]
    scd = values[labels == "SCD"]
    effects = np.empty(replicates, dtype=float)
    for index in range(replicates):
        effects[index] = generator.choice(scd, len(scd), replace=True).mean() - generator.choice(
            healthy, len(healthy), replace=True
        ).mean()
    low, high = np.quantile(effects, [0.025, 0.975])
    return float(low), float(high)


def analyze_stability(
    donor_fates: pd.DataFrame, config: dict[str, Any], cohort: str
) -> pd.DataFrame:
    required_lodo = int(
        config["model"][
            "discovery_leave_out_required" if cohort == "discovery" else "validation_leave_out_required"
        ]
    )
    estimator_results: dict[tuple[str, str], dict[str, Any]] = {}
    for (lineage, estimator), group in donor_fates.groupby(["lineage", "estimator"], sort=False):
        values = group["mean_fate_probability"].to_numpy(dtype=float)
        labels = group["condition"].to_numpy(dtype=str)
        donors = group["donor"].to_numpy(dtype=str)
        if set(labels) != {"healthy", "SCD"}:
            raise RuntimeError(f"{cohort} contrast lacks a required condition")
        effect = _effect(values, labels)
        permutation_distribution = _exact_permutation_distribution(values, labels, donors)
        permutation_count = len(permutation_distribution)
        permutation_tail = sum(
            abs(item["effect"]) >= abs(effect) - 1e-12 for item in permutation_distribution
        ) / permutation_count
        lodo = leave_one_out_effects(values, labels)
        sign = int(np.sign(effect))
        lodo_same = int((np.sign(lodo) == sign).sum()) if sign else 0
        interval_low, interval_high = bootstrap_interval(
            values,
            labels,
            int(config["random_seed"]) + LINEAGES.index(lineage) * 17 + ESTIMATORS.index(estimator),
        )
        estimator_results[(lineage, estimator)] = {
            "effect": effect,
            "direction": sign,
            "permutation_tail_fraction": permutation_tail,
            "permutation_count": permutation_count,
            "permutation_assignments_json": json.dumps(
                permutation_distribution, separators=(",", ":")
            ),
            "leave_one_out_same_direction": lodo_same,
            "leave_one_out_total": len(lodo),
            "leave_one_out_minimum_effect": float(lodo.min()),
            "leave_one_out_maximum_effect": float(lodo.max()),
            "leave_one_out_effects_json": json.dumps(
                [
                    {
                        "excluded_donor": str(donor),
                        "effect": float(lodo_effect),
                        "direction": int(np.sign(lodo_effect)),
                    }
                    for donor, lodo_effect in zip(donors, lodo, strict=True)
                ],
                separators=(",", ":"),
            ),
            "bootstrap_95_low": interval_low,
            "bootstrap_95_high": interval_high,
        }

    rows = []
    for lineage in LINEAGES:
        primary = estimator_results[(lineage, "absorbing_markov")]
        secondary = estimator_results[(lineage, "terminal_projection")]
        concordant = primary["direction"] != 0 and primary["direction"] == secondary["direction"]
        lodo_pass = (
            primary["leave_one_out_same_direction"] >= required_lodo
            and secondary["leave_one_out_same_direction"] >= required_lodo
        )
        null_pass = (
            primary["permutation_tail_fraction"] <= config["model"]["permutation_tail_maximum"]
            and secondary["permutation_tail_fraction"] <= config["model"]["permutation_tail_maximum"]
        )
        rows.append(
            {
                "cohort": cohort,
                "lineage": lineage,
                "absorbing_markov_effect": primary["effect"],
                "terminal_projection_effect": secondary["effect"],
                "direction": primary["direction"] if concordant else 0,
                "estimator_agreement": concordant,
                "absorbing_markov_permutation_tail_fraction": primary["permutation_tail_fraction"],
                "terminal_projection_permutation_tail_fraction": secondary[
                    "permutation_tail_fraction"
                ],
                "absorbing_markov_permutation_assignments_json": primary[
                    "permutation_assignments_json"
                ],
                "terminal_projection_permutation_assignments_json": secondary[
                    "permutation_assignments_json"
                ],
                "exact_permutations": primary["permutation_count"],
                "absorbing_markov_leave_one_out_same_direction": primary[
                    "leave_one_out_same_direction"
                ],
                "terminal_projection_leave_one_out_same_direction": secondary[
                    "leave_one_out_same_direction"
                ],
                "leave_one_out_total": primary["leave_one_out_total"],
                "absorbing_markov_leave_one_out_effects_json": primary[
                    "leave_one_out_effects_json"
                ],
                "terminal_projection_leave_one_out_effects_json": secondary[
                    "leave_one_out_effects_json"
                ],
                "absorbing_markov_bootstrap_95_low": primary["bootstrap_95_low"],
                "absorbing_markov_bootstrap_95_high": primary["bootstrap_95_high"],
                "terminal_projection_bootstrap_95_low": secondary["bootstrap_95_low"],
                "terminal_projection_bootstrap_95_high": secondary["bootstrap_95_high"],
                "leave_one_out_passed": lodo_pass,
                "null_control_passed": null_pass,
                "stable": bool(concordant and lodo_pass and null_pass),
            }
        )
    return pd.DataFrame.from_records(rows)


def build_replication_table(
    discovery: pd.DataFrame, validation: pd.DataFrame
) -> pd.DataFrame:
    rows = []
    transport = {
        "discovery": "plerixafor-mobilized HSPCs in a gene-therapy collection context",
        "validation": "pediatric bone-marrow CD34-positive cells collected in a distinct study",
    }
    for lineage in LINEAGES:
        discovery_row = discovery.loc[discovery["lineage"] == lineage].iloc[0]
        validation_row = validation.loc[validation["lineage"] == lineage].iloc[0]
        replicated = bool(
            discovery_row["stable"]
            and validation_row["stable"]
            and discovery_row["direction"] == validation_row["direction"]
        )
        for cohort, row in (("discovery", discovery_row), ("validation", validation_row)):
            rows.append(
                {
                    "cohort": cohort,
                    "lineage": lineage,
                    "direction": int(row["direction"]),
                    "magnitude": float(row["absorbing_markov_effect"]),
                    "uncertainty": (
                        f"donor bootstrap 95% [{row['absorbing_markov_bootstrap_95_low']:.6g}, "
                        f"{row['absorbing_markov_bootstrap_95_high']:.6g}]"
                    ),
                    "stability": bool(row["stable"]),
                    "estimator_agreement": bool(row["estimator_agreement"]),
                    "transport_limitations": transport[cohort],
                    "replication_concordant": replicated,
                }
            )
    return pd.DataFrame.from_records(rows)


def classify_result(discovery: pd.DataFrame, validation: pd.DataFrame) -> tuple[str, list[str]]:
    discovery_stable = discovery.loc[discovery["stable"]]
    if discovery_stable.empty:
        return "no_stable_inferred_allocation_difference", []
    replicated = []
    for lineage in discovery_stable["lineage"]:
        discovery_row = discovery.loc[discovery["lineage"] == lineage].iloc[0]
        validation_row = validation.loc[validation["lineage"] == lineage].iloc[0]
        if validation_row["stable"] and validation_row["direction"] == discovery_row["direction"]:
            replicated.append(str(lineage))
    if replicated:
        return "replicated_inferred_allocation_difference", replicated
    return "dataset_specific_inferred_allocation_difference", []
