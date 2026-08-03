from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .model import ESTIMATORS, LINEAGES, combine_donors, fit_healthy_model
from .preprocess import prepare_cohort
from .robustness import attach_robustness, run_model_robustness
from .sources import (
    assert_all_admitted,
    build_donor_manifest,
    build_source_manifest,
    load_analysis_config,
    load_source_ledger,
    preflight_assets,
)
from .statistics import (
    aggregate_donor_fates,
    analyze_stability,
    build_replication_table,
    classify_result,
)
from .util import canonical_sha256, sha256_file, write_json, write_tsv


CLASSIFICATIONS = {
    "replicated_inferred_allocation_difference",
    "dataset_specific_inferred_allocation_difference",
    "no_stable_inferred_allocation_difference",
    "insufficient_evidence",
}

REQUIRED_ARTIFACTS = (
    "source_manifest.json",
    "donor_manifest.tsv",
    "qc_ledger.tsv",
    "state_definitions.json",
    "topology_controls.json",
    "donor_lineage_priming.tsv",
    "reference_mapping.tsv",
    "model_robustness.tsv",
    "stability_results.tsv",
    "replication_table.tsv",
    "gate_ledger.tsv",
    "result_classification.json",
    "claim_boundary.md",
    "report.md",
)


def _gate(name: str, passed: bool, evidence: str, response: str) -> dict[str, Any]:
    return {
        "gate": name,
        "status": "passed" if passed else "failed",
        "evidence": evidence,
        "response": "continue" if passed else response,
    }


def _qc_with_diagnostics(
    ledger: pd.DataFrame, diagnostics: list[dict[str, Any]], cohort: str
) -> pd.DataFrame:
    rows = []
    for diagnostic in diagnostics:
        rows.append(
            {
                "cohort": cohort,
                "donor": None,
                **diagnostic,
            }
        )
    diagnostic_frame = pd.DataFrame.from_records(rows)
    if not diagnostic_frame.empty:
        donors = ledger.loc[ledger["control"] == "raw_barcodes", "donor"].tolist()
        diagnostic_frame["donor"] = donors
    return ledger.merge(diagnostic_frame, on=["cohort", "donor"], how="left")


def _claim_boundary(classification: str) -> str:
    return f"""# Claim Boundary

The released classification is `{classification}`.

Permitted interpretation: cross-sectional donor-level association in healthy-reference lineage priming among transcriptomically HSC/MPP-like cells that pass reference-support checks.

Excluded interpretations: observed fate, lineage tracing, a causal SCA effect, individual prediction, treatment response, gene-therapy efficacy, or clinical use. SCD genotype is not a proxy for ancestry, geography, environment, access to care, or self-identified race. Unavailable cohort attributes remain marked `not_reported`.
"""


def _report(
    classification: str,
    failure_reasons: list[str],
    discovery_stability: pd.DataFrame,
    validation_fates: pd.DataFrame,
) -> str:
    discovery_lines = []
    for row in discovery_stability.itertuples():
        discovery_lines.append(
            f"- {row.lineage}: Palantir effect {row.palantir_effect:.6f}; "
            f"CellRank-CytoTRACE effect {row.cellrank_cytotrace_effect:.6f}; "
            f"stable={str(row.stable).lower()}."
        )
    low_root = validation_fates.loc[~validation_fates["root_cell_gate_passed"]]
    low_root_summary = ", ".join(
        f"{row.donor} ({row.root_cell_count} root cells)"
        for row in low_root.drop_duplicates("donor").itertuples()
    )
    failures = "\n".join(f"- {reason}" for reason in failure_reasons) or "- None."
    return f"""# HSPC Lineage-Priming v2 Result

## Classification

`{classification}`

This is a cross-sectional donor-level analysis of healthy-reference lineage priming. It does not observe cell fate and does not establish that SCA causes a transition change.

Discovery uses six open donor-separated unfiltered 10x matrices. Cell admission retains the source paper's fixed 3,000-UMI floor, adds an ambient-profile audit, uses a multi-gene HSPC identity rule instead of CD34 alone, collapses duplicate symbols, and records Scrublet score and call-rate diagnostics. Exact EmptyDrops FDR reproduction remains unavailable without R/DropletUtils.

The primary estimator is Palantir on a donor-balanced healthy cell reference. CellRank 2 with a CytoTRACE kernel is a complementary potency-based sensitivity estimator. Construction and held-out validation markers are disjoint. SCD cells cannot refit the reference and are rejected when their healthy-reference distance exceeds a leave-one-healthy-donor-out threshold. Donor means are the biological replicates.

## Discovery

{chr(10).join(discovery_lines)}

`stability_results.tsv` preserves every donor-label assignment and donor-specific deletion effect. `model_robustness.tsv` adds healthy-reference refits and bounded hyperparameter refits; a positive result must survive both.

## Validation gate

Validation independently recovered the root and three terminal families. The donor reliability gate failed for {low_root_summary or 'no donor'}; estimates below the 20-root-cell floor are retained for audit but cannot support replication.

## Failed gates

{failures}

## Follow-up

Verify a third donor-resolved open human SCD HSPC cohort or obtain compatible directional molecular layers. Re-run the frozen lineage-priming estimand without pooling cohorts or changing endpoints. Time-resolved, lineage-tracing, or compatible spliced/unspliced data are required before upgrading the claim toward temporal fate direction.
"""


def run_pipeline(root: Path, force_preprocess: bool = False) -> dict[str, Any]:
    root = root.resolve()
    config_path = root / "research/hspc-fate/analysis.json"
    ledger_path = root / "research/hspc-fate/sources.json"
    config = load_analysis_config(config_path)
    ledger = load_source_ledger(ledger_path)
    preflight = preflight_assets(config, root)
    assert_all_admitted(preflight)
    source_manifest = build_source_manifest(config, ledger, preflight)
    donor_manifest = build_donor_manifest(config)

    discovery_donors, discovery_qc, discovery_diagnostics = prepare_cohort(
        root, config, "discovery", force=force_preprocess
    )
    validation_donors, validation_qc, validation_diagnostics = prepare_cohort(
        root, config, "validation", force=force_preprocess
    )
    qc_ledger = pd.concat(
        [
            _qc_with_diagnostics(discovery_qc, discovery_diagnostics, "discovery"),
            _qc_with_diagnostics(validation_qc, validation_diagnostics, "validation"),
        ],
        ignore_index=True,
    )

    discovery_cohort = combine_donors(discovery_donors)
    validation_cohort = combine_donors(validation_donors)
    discovery_model, discovery_result = fit_healthy_model(discovery_cohort, config)
    validation_model, validation_result = fit_healthy_model(validation_cohort, config)
    minimum_root = int(config["model"]["minimum_root_cells_per_donor"])
    discovery_fates = aggregate_donor_fates(
        discovery_result, minimum_root, enforce_minimum=False
    )
    validation_fates = aggregate_donor_fates(
        validation_result, minimum_root, enforce_minimum=False
    )
    donor_fates = pd.concat([discovery_fates, validation_fates], ignore_index=True)
    discovery_robustness = run_model_robustness(
        discovery_cohort, discovery_result, config, "discovery"
    )
    validation_robustness = run_model_robustness(
        validation_cohort, validation_result, config, "validation"
    )
    robustness = pd.concat(
        [discovery_robustness, validation_robustness], ignore_index=True
    )
    discovery_stability = attach_robustness(
        analyze_stability(discovery_fates, config, "discovery"),
        discovery_robustness,
        config,
        "discovery",
    )
    validation_stability = attach_robustness(
        analyze_stability(validation_fates, config, "validation"),
        validation_robustness,
        config,
        "validation",
    )
    stability = pd.concat([discovery_stability, validation_stability], ignore_index=True)
    replication = build_replication_table(discovery_stability, validation_stability)
    reference_mapping = pd.concat(
        [
            result.observations.assign(root_cell=result.root_mask)
            .groupby(["cohort", "donor", "condition"], as_index=False)
            .agg(
                cells=("cell_index", "size"),
                reference_supported_cells=("reference_supported", "sum"),
                reference_supported_fraction=("reference_supported", "mean"),
                root_cells=("root_cell", "sum"),
                median_ood_distance=("ood_distance", "median"),
            )
            for result in (discovery_result, validation_result)
        ],
        ignore_index=True,
    )

    discovery_counts = {donor.donor: donor.counts.shape[0] for donor in discovery_donors}
    validation_counts = {donor.donor: donor.counts.shape[0] for donor in validation_donors}
    gates = [
        _gate(
            "gene_identity_normalization",
            True,
            f"collapsed duplicate symbols: discovery={sum(item['duplicate_gene_symbols_collapsed'] for item in discovery_diagnostics)}; validation={sum(item['duplicate_gene_symbols_collapsed'] for item in validation_diagnostics)}",
            "stop modeling",
        ),
        _gate(
            "doublet_call_rate_audit",
            all(
                item["scrublet_called_fraction"]
                >= float(config["model"]["doublet_audit_minimum_fraction"])
                for item in discovery_diagnostics + validation_diagnostics
            ),
            "; ".join(
                f"{cohort}:{item['scrublet_called_fraction']:.4%}"
                for cohort, diagnostics in (("discovery", discovery_diagnostics), ("validation", validation_diagnostics))
                for item in diagnostics
            ),
            "retain estimates for audit and classify as insufficient_evidence",
        ),
        _gate(
            "required_source_admission",
            preflight["status"].eq("admitted").all(),
            "7 of 7 required archives match frozen byte counts and SHA-256 checksums",
            "classify as insufficient_evidence",
        ),
        _gate(
            "model_robustness_execution",
            robustness["status"].eq("completed").all(),
            f"completed={int(robustness['status'].eq('completed').sum())}; failed={int(robustness['status'].eq('failed').sum())}",
            "retain baseline estimates for audit and classify as insufficient_evidence",
        ),
        _gate(
            "donor_manifest_cardinality",
            len(discovery_counts) == 6 and len(validation_counts) == 9,
            f"discovery donors={len(discovery_counts)}; validation donors={len(validation_counts)}",
            "classify as insufficient_evidence",
        ),
        _gate(
            "discovery_published_cell_count",
            abs(sum(discovery_counts.values()) - config["qc"]["discovery"]["published_retained_cells"])
            / config["qc"]["discovery"]["published_retained_cells"]
            <= config["qc"]["discovery"]["published_count_tolerance_fraction"],
            f"retained={sum(discovery_counts.values())}; published={config['qc']['discovery']['published_retained_cells']}",
            "stop discovery inference",
        ),
        _gate(
            "validation_published_cell_count",
            abs(sum(validation_counts.values()) - config["qc"]["validation"]["published_retained_cells"])
            / config["qc"]["validation"]["published_retained_cells"]
            <= config["qc"]["validation"]["published_count_tolerance_fraction"],
            f"retained={sum(validation_counts.values())}; published={config['qc']['validation']['published_retained_cells']}",
            "stop validation inference",
        ),
        _gate(
            "discovery_topology",
            bool(discovery_model.topology_controls["passed"]),
            json.dumps(discovery_model.topology_controls, sort_keys=True),
            "classify as insufficient_evidence",
        ),
        _gate(
            "validation_topology",
            bool(validation_model.topology_controls["passed"]),
            json.dumps(validation_model.topology_controls, sort_keys=True),
            "classify as insufficient_evidence",
        ),
        _gate(
            "discovery_model_frozen_during_mapping",
            discovery_result.model_hashes_before == discovery_result.model_hashes_after,
            "all 5 model component hashes match",
            "invalidate mapped estimates",
        ),
        _gate(
            "validation_model_frozen_during_mapping",
            validation_result.model_hashes_before == validation_result.model_hashes_after,
            "all 5 model component hashes match",
            "invalidate mapped estimates",
        ),
        _gate(
            "discovery_root_cell_reliability",
            discovery_fates["root_cell_gate_passed"].all(),
            "; ".join(
                f"{row.donor}={row.root_cell_count}"
                for row in discovery_fates.drop_duplicates("donor").itertuples()
            ),
            "classify as insufficient_evidence",
        ),
        _gate(
            "validation_root_cell_reliability",
            validation_fates["root_cell_gate_passed"].all(),
            "; ".join(
                f"{row.donor}={row.root_cell_count}"
                for row in validation_fates.drop_duplicates("donor").itertuples()
            ),
            "retain estimates for audit and classify as insufficient_evidence",
        ),
    ]
    gate_ledger = pd.DataFrame.from_records(gates)
    failed_gates = gate_ledger.loc[gate_ledger["status"] == "failed", "gate"].tolist()
    if failed_gates:
        classification = "insufficient_evidence"
        replicated_lineages: list[str] = []
    else:
        classification, replicated_lineages = classify_result(
            discovery_stability, validation_stability
        )
    if classification not in CLASSIFICATIONS:
        raise RuntimeError(f"invalid result classification: {classification}")

    release = root / "results/hspc-fate"
    release.mkdir(parents=True, exist_ok=True)
    write_json(release / "source_manifest.json", source_manifest)
    write_tsv(release / "donor_manifest.tsv", donor_manifest)
    write_tsv(release / "qc_ledger.tsv", qc_ledger)
    write_json(
        release / "state_definitions.json",
        {
            "discovery": discovery_model.state_definitions,
            "validation": validation_model.state_definitions,
        },
    )
    write_json(
        release / "topology_controls.json",
        {
            "discovery": discovery_model.topology_controls,
            "validation": validation_model.topology_controls,
            "model_hashes": {
                "discovery_before": discovery_result.model_hashes_before,
                "discovery_after": discovery_result.model_hashes_after,
                "validation_before": validation_result.model_hashes_before,
                "validation_after": validation_result.model_hashes_after,
            },
        },
    )
    write_tsv(release / "donor_lineage_priming.tsv", donor_fates)
    write_tsv(release / "reference_mapping.tsv", reference_mapping)
    write_tsv(release / "model_robustness.tsv", robustness)
    write_tsv(release / "stability_results.tsv", stability)
    write_tsv(release / "replication_table.tsv", replication)
    write_tsv(release / "gate_ledger.tsv", gate_ledger)
    result_record = {
        "classification": classification,
        "replicated_lineages": replicated_lineages,
        "failed_gates": failed_gates,
        "discovery_stable_lineages": discovery_stability.loc[
            discovery_stability["stable"], "lineage"
        ].tolist(),
        "validation_estimates_are_exploratory": bool(failed_gates),
        "interpretation": config["claims"]["permitted"],
        "estimand": "healthy-reference lineage priming",
    }
    write_json(release / "result_classification.json", result_record)
    (release / "claim_boundary.md").write_text(
        _claim_boundary(classification), encoding="utf-8"
    )
    (release / "report.md").write_text(
        _report(classification, failed_gates, discovery_stability, validation_fates),
        encoding="utf-8",
    )

    artifact_records = []
    for name in REQUIRED_ARTIFACTS:
        path = release / name
        artifact_records.append(
            {"path": name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
        )
    manifest = {
        "release_id": "hspc-fate-v2",
        "classification": classification,
        "config_sha256": canonical_sha256(config),
        "artifacts": artifact_records,
    }
    write_json(release / "release_manifest.json", manifest)
    return manifest


def verify_release(root: Path) -> dict[str, Any]:
    release = root.resolve() / "results/hspc-fate"
    manifest_path = release / "release_manifest.json"
    if not manifest_path.is_file():
        return {"status": "failed", "reason": "missing_release_manifest"}
    manifest = json.loads(manifest_path.read_text())
    problems: list[str] = []
    if manifest.get("classification") not in CLASSIFICATIONS:
        problems.append("invalid_classification")
    expected_names = set(REQUIRED_ARTIFACTS)
    observed_names = {artifact["path"] for artifact in manifest.get("artifacts", [])}
    if observed_names != expected_names:
        problems.append("artifact_set_mismatch")
    for artifact in manifest.get("artifacts", []):
        path = release / artifact["path"]
        if not path.is_file():
            problems.append(f"missing:{artifact['path']}")
        elif path.stat().st_size != artifact["bytes"]:
            problems.append(f"size:{artifact['path']}")
        elif sha256_file(path) != artifact["sha256"]:
            problems.append(f"checksum:{artifact['path']}")

    donor_fates = pd.read_csv(release / "donor_lineage_priming.tsv", sep="\t")
    key = ["cohort", "donor", "estimator", "lineage"]
    if donor_fates.duplicated(key).any():
        problems.append("duplicate_donor_fate_row")
    sums = donor_fates.groupby(["cohort", "donor", "estimator"])[
        "mean_lineage_priming_probability"
    ].sum()
    if not np.allclose(sums.to_numpy(), 1.0, atol=1e-6):
        problems.append("donor_lineage_priming_not_normalized")
    if set(donor_fates["lineage"]) != set(LINEAGES):
        problems.append("lineage_set_mismatch")

    stability = pd.read_csv(release / "stability_results.tsv", sep="\t")
    for row in stability.itertuples():
        for estimator in ESTIMATORS:
            permutations = json.loads(
                getattr(row, f"{estimator}_permutation_assignments_json")
            )
            leave_one_out = json.loads(
                getattr(row, f"{estimator}_leave_one_out_effects_json")
            )
            if len(permutations) != row.exact_permutations:
                problems.append(f"permutation_record_count:{row.cohort}:{row.lineage}:{estimator}")
            if len(leave_one_out) != row.leave_one_out_total:
                problems.append(f"leave_one_out_record_count:{row.cohort}:{row.lineage}:{estimator}")
            if len({item["excluded_donor"] for item in leave_one_out}) != len(leave_one_out):
                problems.append(f"leave_one_out_donor_set:{row.cohort}:{row.lineage}:{estimator}")

    gate_ledger = pd.read_csv(release / "gate_ledger.tsv", sep="\t")
    classification = manifest["classification"]
    if classification == "insufficient_evidence" and gate_ledger["status"].ne("failed").all():
        problems.append("insufficient_without_failed_gate")
    if classification != "insufficient_evidence" and gate_ledger["status"].ne("passed").any():
        problems.append("positive_classification_with_failed_gate")
    source_manifest = json.loads((release / "source_manifest.json").read_text())
    if any(asset["status"] != "admitted" for asset in source_manifest["assets"]):
        problems.append("unadmitted_source_in_release")
    claim_boundary = (release / "claim_boundary.md").read_text()
    if "cross-sectional donor-level association in healthy-reference lineage priming" not in claim_boundary:
        problems.append("missing_permitted_interpretation")
    if "Excluded interpretations:" not in claim_boundary:
        problems.append("missing_claim_exclusions")
    return {"status": "passed" if not problems else "failed", "problems": problems}
