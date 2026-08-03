from __future__ import annotations

from copy import deepcopy
from typing import Any

import numpy as np
import pandas as pd

from .model import ESTIMATORS, CombinedCohort, LINEAGES, FateResult, fit_healthy_model
from .statistics import aggregate_donor_fates


def _subset_cohort(cohort: CombinedCohort, keep: np.ndarray) -> CombinedCohort:
    return CombinedCohort(
        counts=cohort.counts[keep].tocsr(),
        genes=cohort.genes,
        observations=cohort.observations.loc[keep].reset_index(drop=True),
    )


def _effects(
    donor_fates: pd.DataFrame,
    cohort: str,
    perturbation: str,
    specification: str,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for lineage in LINEAGES:
        for estimator in ESTIMATORS:
            group = donor_fates.loc[
                donor_fates["lineage"].eq(lineage)
                & donor_fates["estimator"].eq(estimator)
            ]
            healthy = group.loc[
                group["condition"].eq("healthy"), "mean_lineage_priming_probability"
            ].to_numpy(dtype=float)
            scd = group.loc[
                group["condition"].eq("SCD"), "mean_lineage_priming_probability"
            ].to_numpy(dtype=float)
            if not len(healthy) or not len(scd):
                raise RuntimeError("robustness contrast lacks a required condition")
            effect = float(scd.mean() - healthy.mean())
            rows.append(
                {
                    "cohort": cohort,
                    "perturbation": perturbation,
                    "specification": specification,
                    "lineage": lineage,
                    "estimator": estimator,
                    "effect": effect,
                    "direction": int(np.sign(effect)),
                    "detail": "",
                    "status": "completed",
                }
            )
    return rows


def _failure_rows(
    cohort: str, perturbation: str, specification: str, error: Exception
) -> list[dict[str, Any]]:
    return [
        {
            "cohort": cohort,
            "perturbation": perturbation,
            "specification": specification,
            "lineage": lineage,
            "estimator": estimator,
            "effect": float("nan"),
            "direction": 0,
            "detail": f"{type(error).__name__}: {error}",
            "status": "failed",
        }
        for lineage in LINEAGES
        for estimator in ESTIMATORS
    ]


def run_model_robustness(
    cohort: CombinedCohort,
    baseline_result: FateResult,
    config: dict[str, Any],
    cohort_name: str,
) -> pd.DataFrame:
    minimum_root_cells = int(config["model"]["minimum_root_cells_per_donor"])
    baseline_fates = aggregate_donor_fates(
        baseline_result, minimum_root_cells, enforce_minimum=False
    )
    rows = _effects(baseline_fates, cohort_name, "baseline", "baseline")

    donor_conditions = cohort.observations[["donor", "condition"]].drop_duplicates()
    for donor, condition in donor_conditions.itertuples(index=False):
        keep = ~cohort.observations["donor"].eq(donor).to_numpy()
        specification = f"exclude_{donor}"
        if condition == "healthy":
            try:
                _, result = fit_healthy_model(_subset_cohort(cohort, keep), config)
                fates = aggregate_donor_fates(
                    result, minimum_root_cells, enforce_minimum=False
                )
                rows.extend(
                    _effects(fates, cohort_name, "donor_deletion_refit", specification)
                )
            except Exception as error:  # preserved as a scientific failure, not hidden
                rows.extend(
                    _failure_rows(
                        cohort_name, "donor_deletion_refit", specification, error
                    )
                )
        else:
            rows.extend(
                _effects(
                    baseline_fates.loc[~baseline_fates["donor"].eq(donor)],
                    cohort_name,
                    "donor_deletion_fixed_reference",
                    specification,
                )
            )

    for specification in config["model"].get("sensitivity_specs", []):
        varied = deepcopy(config)
        varied["random_seed"] = int(config["random_seed"]) + int(
            specification["random_seed_offset"]
        )
        for key in (
            "highly_variable_genes",
            "pca_components",
            "trajectory_neighbors",
            "reference_cells_per_healthy_donor",
        ):
            varied["model"][key] = specification[key]
        try:
            _, result = fit_healthy_model(cohort, varied)
            fates = aggregate_donor_fates(
                result, minimum_root_cells, enforce_minimum=False
            )
            rows.extend(
                _effects(
                    fates,
                    cohort_name,
                    "hyperparameter_refit",
                    str(specification["id"]),
                )
            )
        except Exception as error:  # preserved as a scientific failure, not hidden
            rows.extend(
                _failure_rows(
                    cohort_name,
                    "hyperparameter_refit",
                    str(specification["id"]),
                    error,
                )
            )
    return pd.DataFrame.from_records(rows)


def attach_robustness(
    stability: pd.DataFrame,
    robustness: pd.DataFrame,
    config: dict[str, Any],
    cohort: str,
) -> pd.DataFrame:
    updated = stability.copy()
    required_lodo = int(
        config["model"][
            "discovery_leave_out_required"
            if cohort == "discovery"
            else "validation_leave_out_required"
        ]
    )
    for row_index, row in updated.iterrows():
        lineage = str(row["lineage"])
        estimator_checks: list[bool] = []
        for estimator in ESTIMATORS:
            subset = robustness.loc[
                robustness["lineage"].eq(lineage)
                & robustness["estimator"].eq(estimator)
            ]
            baseline = subset.loc[subset["perturbation"].eq("baseline")].iloc[0]
            deletions = subset.loc[subset["perturbation"].str.startswith("donor_deletion")]
            specifications = subset.loc[subset["perturbation"].eq("hyperparameter_refit")]
            baseline_direction = int(baseline["direction"])
            deletion_same = int(
                (
                    deletions["status"].eq("completed")
                    & deletions["direction"].eq(baseline_direction)
                ).sum()
            )
            specification_passed = bool(
                not specifications.empty
                and specifications["status"].eq("completed").all()
                and specifications["direction"].eq(baseline_direction).all()
            )
            prefix = estimator
            updated.loc[row_index, f"{prefix}_model_refit_leave_out_same_direction"] = (
                deletion_same
            )
            updated.loc[row_index, f"{prefix}_hyperparameter_refit_passed"] = (
                specification_passed
            )
            estimator_checks.append(
                baseline_direction != 0
                and deletion_same >= required_lodo
                and specification_passed
            )
        robustness_passed = bool(all(estimator_checks))
        updated.loc[row_index, "model_robustness_passed"] = robustness_passed
        updated.loc[row_index, "stable"] = bool(row["stable"] and robustness_passed)
    return updated
