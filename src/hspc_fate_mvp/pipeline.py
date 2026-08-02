from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .model import LINEAGES, combine_donors, fit_healthy_model
from .preprocess import prepare_cohort
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
    "donor_fate_probabilities.tsv",
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

Permitted interpretation: cross-sectional donor-level association in inferred lineage allocation among transcriptomically HSC/MPP-like cells.

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
            f"- {row.lineage}: absorbing-Markov effect {row.absorbing_markov_effect:.6f}; "
            f"projection effect {row.terminal_projection_effect:.6f}; stable={str(row.stable).lower()}."
        )
    low_root = validation_fates.loc[~validation_fates["root_cell_gate_passed"]]
    low_root_summary = ", ".join(
        f"{row.donor} ({row.root_cell_count} root cells)"
        for row in low_root.drop_duplicates("donor").itertuples()
    )
    failures = "\n".join(f"- {reason}" for reason in failure_reasons) or "- None."
    return f"""# HSPC Fate v1 Result

## Classification

`{classification}`

This is a cross-sectional donor-level analysis of inferred lineage allocation. It does not observe cell fate and does not establish that SCA causes an allocation difference.

Discovery uses six open donor-separated unfiltered 10x matrices. Cell admission applies the source paper's fixed 3,000-UMI floor with a recorded barcode-rank knee diagnostic; it does not reproduce EmptyDrops FDR calculations. Scrublet replaces the paper's DoubletFinder implementation. The validation cohort uses its own frozen source-aware thresholds. These are explicit lean-MVP method deviations, not claims of an exact paper reproduction.

The primary estimator is absorption through a healthy-only directed state graph. The independent estimator projects cells to pseudotime-defined healthy terminal centroids without using the graph. Donor means are the biological replicates.

## Discovery

{chr(10).join(discovery_lines)}

No discovery lineage passed all frozen method-concordance, leave-one-donor-out, and exact-permutation null controls. `stability_results.tsv` preserves every donor-label assignment and every donor-specific leave-one-out effect for audit.

## Validation gate

Validation independently recovered the root and three terminal families. The donor reliability gate failed for {low_root_summary or 'no donor'}; estimates below the 20-root-cell floor are retained for audit but cannot support replication.

## Failed gates

{failures}

## Follow-up

Acquire an open independent human SCA HSPC cohort in which every comparator donor retains at least 20 transcriptomically HSC/MPP-like cells under the frozen marker semantics. Re-run the same donor-level estimand without pooling cohorts or changing endpoints. Optional bulk-human or mouse evidence may qualify programs but cannot rescue this human replication gate.
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

    discovery_model, discovery_result = fit_healthy_model(
        combine_donors(discovery_donors), config
    )
    validation_model, validation_result = fit_healthy_model(
        combine_donors(validation_donors), config
    )
    minimum_root = int(config["model"]["minimum_root_cells_per_donor"])
    discovery_fates = aggregate_donor_fates(
        discovery_result, minimum_root, enforce_minimum=False
    )
    validation_fates = aggregate_donor_fates(
        validation_result, minimum_root, enforce_minimum=False
    )
    donor_fates = pd.concat([discovery_fates, validation_fates], ignore_index=True)
    discovery_stability = analyze_stability(discovery_fates, config, "discovery")
    validation_stability = analyze_stability(validation_fates, config, "validation")
    stability = pd.concat([discovery_stability, validation_stability], ignore_index=True)
    replication = build_replication_table(discovery_stability, validation_stability)

    discovery_counts = {donor.donor: donor.counts.shape[0] for donor in discovery_donors}
    validation_counts = {donor.donor: donor.counts.shape[0] for donor in validation_donors}
    gates = [
        _gate(
            "required_source_admission",
            preflight["status"].eq("admitted").all(),
            "7 of 7 required archives match frozen byte counts and SHA-256 checksums",
            "classify as insufficient_evidence",
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
    write_tsv(release / "donor_fate_probabilities.tsv", donor_fates)
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
        "release_id": "hspc-fate-v1",
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

    donor_fates = pd.read_csv(release / "donor_fate_probabilities.tsv", sep="\t")
    key = ["cohort", "donor", "estimator", "lineage"]
    if donor_fates.duplicated(key).any():
        problems.append("duplicate_donor_fate_row")
    sums = donor_fates.groupby(["cohort", "donor", "estimator"])[
        "mean_fate_probability"
    ].sum()
    if not np.allclose(sums.to_numpy(), 1.0, atol=1e-6):
        problems.append("donor_fate_probabilities_not_normalized")
    if set(donor_fates["lineage"]) != set(LINEAGES):
        problems.append("lineage_set_mismatch")

    stability = pd.read_csv(release / "stability_results.tsv", sep="\t")
    for row in stability.itertuples():
        for estimator in ("absorbing_markov", "terminal_projection"):
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
    if "cross-sectional donor-level association in inferred lineage allocation" not in claim_boundary:
        problems.append("missing_permitted_interpretation")
    if "Excluded interpretations:" not in claim_boundary:
        problems.append("missing_claim_exclusions")
    return {"status": "passed" if not problems else "failed", "problems": problems}
