from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd

from . import __version__
from .ranking import rank_hypotheses
from .reference import ReferenceData, load_or_build_reference
from .signals import build_signal_ledger
from .sources import assert_all_admitted, load_registry, preflight_sources, source_by_id
from .transitions import generate_transitions, map_signals_to_edges
from .util import (
    canonical_json_bytes,
    sha256_file,
    write_gzip_tsv,
    write_json,
    write_tsv,
)
from .validation import (
    validate_analysis,
    validate_compute_preflight,
    validate_git_boundary,
    validate_local_markdown_links,
)


RELEASE_FILES = [
    "source_registry.json",
    "state_metadata.tsv",
    "state_programs.tsv",
    "healthy_reference_state_centroids.tsv.gz",
    "candidate_adjacencies.tsv",
    "sca_signal_ledger.tsv",
    "sca_transition_mapping.tsv",
    "ranked_hypotheses.tsv",
    "decision_log.tsv",
    "validation_report.json",
    "compute_preflight.json",
    "follow_up_plan.md",
    "report.md",
    "provenance_manifest.json",
]


def run_pipeline(root: Path, force_reference: bool = False) -> dict[str, Any]:
    source_registry_path = root / "config/sources.json"
    analysis_config_path = root / "config/analysis.json"
    registry = load_registry(source_registry_path)
    analysis = json.loads(analysis_config_path.read_text(encoding="utf-8"))
    sources = source_by_id(registry)

    preflight = preflight_sources(registry, root)
    assert_all_admitted(preflight)
    signal_ledger = build_signal_ledger(root, sources, analysis["signatures"])

    references: list[ReferenceData] = []
    for reference_spec in analysis["reference_sources"]:
        references.append(
            load_or_build_reference(
                root=root,
                source=sources[reference_spec["source_id"]],
                reference_spec=reference_spec,
                pseudobulk_config=analysis["pseudobulk"],
                cache_dir=root / "data/derived",
                force=force_reference,
            )
        )
    transitions = generate_transitions(references, analysis["adjacency"])
    mapping, compatibility = map_signals_to_edges(
        transitions.edges, transitions.state_centroids, signal_ledger
    )
    hypotheses = rank_hypotheses(
        transitions.edges, compatibility, analysis["ranking"]
    )
    state_metadata = pd.concat([reference.states for reference in references], ignore_index=True)
    exclusions = pd.concat([reference.exclusions for reference in references], ignore_index=True)
    decision_log = build_decision_log(preflight, exclusions, transitions.edges, mapping)
    validation = validate_analysis(
        preflight,
        state_metadata,
        transitions.edges,
        signal_ledger,
        mapping,
        hypotheses,
        analysis,
    )
    validation["checks"].append(validate_git_boundary(root))
    validation["check_count"] += 1
    validation["passed_count"] += int(validation["checks"][-1]["passed"])
    validation["failed_count"] += int(not validation["checks"][-1]["passed"])
    validation["status"] = "passed" if validation["failed_count"] == 0 else "failed"
    compute_checks = validate_compute_preflight(
        analysis["compute_preflight"],
        sum(int(source["expected_bytes"]) for source in registry["sources"]),
    )
    validation["checks"].extend(compute_checks)
    validation["check_count"] += len(compute_checks)
    validation["passed_count"] += sum(check["passed"] for check in compute_checks)
    validation["failed_count"] += sum(not check["passed"] for check in compute_checks)
    validation["status"] = "passed" if validation["failed_count"] == 0 else "failed"
    if validation["status"] != "passed":
        failed = [check["name"] for check in validation["checks"] if not check["passed"]]
        raise RuntimeError(f"analysis validation failed: {', '.join(failed)}")

    output_dir = root / "results/mvp"
    output_dir.mkdir(parents=True, exist_ok=True)
    write_json(output_dir / "source_registry.json", _observed_registry(registry, preflight))
    write_tsv(output_dir / "state_metadata.tsv", state_metadata)
    write_tsv(output_dir / "state_programs.tsv", transitions.state_programs)
    write_gzip_tsv(
        output_dir / "healthy_reference_state_centroids.tsv.gz",
        _combined_centroids(transitions.state_centroids),
    )
    write_tsv(output_dir / "candidate_adjacencies.tsv", transitions.edges)
    write_tsv(output_dir / "sca_signal_ledger.tsv", signal_ledger)
    write_tsv(output_dir / "sca_transition_mapping.tsv", mapping)
    write_tsv(output_dir / "ranked_hypotheses.tsv", hypotheses)
    write_tsv(output_dir / "decision_log.tsv", decision_log)
    write_json(output_dir / "validation_report.json", validation)
    write_json(
        output_dir / "compute_preflight.json",
        {
            **analysis["compute_preflight"],
            "source_total_bytes": sum(
                int(source["expected_bytes"]) for source in registry["sources"]
            ),
            "interpretation": "observed local feasibility benchmark; not a performance guarantee on other hosts",
        },
    )
    (output_dir / "follow_up_plan.md").write_text(
        build_follow_up_plan(hypotheses), encoding="utf-8"
    )
    (output_dir / "report.md").write_text(
        build_report(state_metadata, transitions.edges, signal_ledger, hypotheses, validation),
        encoding="utf-8",
    )
    provenance = build_provenance_manifest(
        root,
        output_dir,
        source_registry_path,
        analysis_config_path,
        references,
    )
    write_json(output_dir / "provenance_manifest.json", provenance)
    release = build_release_manifest(output_dir, validation)
    write_json(output_dir / "release_manifest.json", release)
    return release


def verify_release(root: Path) -> dict[str, Any]:
    output_dir = root / "results/mvp"
    manifest_path = output_dir / "release_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError("release manifest is missing")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    checks: list[dict[str, Any]] = []
    for artifact in manifest["artifacts"]:
        path = output_dir / artifact["path"]
        exists = path.is_file()
        size_match = exists and path.stat().st_size == artifact["bytes"]
        hash_match = exists and sha256_file(path) == artifact["sha256"]
        checks.append(
            {
                "artifact": artifact["path"],
                "exists": exists,
                "size_match": size_match,
                "sha256_match": hash_match,
                "passed": bool(exists and size_match and hash_match),
            }
        )
    validation = json.loads((output_dir / "validation_report.json").read_text(encoding="utf-8"))
    git_check = validate_git_boundary(root)
    link_check = validate_local_markdown_links(root)
    passed = (
        all(check["passed"] for check in checks)
        and validation["status"] == "passed"
        and git_check["passed"]
        and link_check["passed"]
    )
    return {
        "status": "passed" if passed else "failed",
        "artifact_checks": checks,
        "analysis_validation": validation["status"],
        "git_boundary": git_check,
        "local_links": link_check,
    }


def _observed_registry(registry: dict[str, Any], preflight: pd.DataFrame) -> dict[str, Any]:
    observed = preflight.set_index("source_id").to_dict(orient="index")
    return {
        "schema_version": registry["schema_version"],
        "frozen_on": registry["frozen_on"],
        "sources": [
            {**source, "admission": observed[source["id"]]}
            for source in registry["sources"]
        ],
    }


def _combined_centroids(centroids_by_organ: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for organ, frame in sorted(centroids_by_organ.items()):
        part = frame.copy()
        part.insert(0, "state", part.index)
        part.insert(0, "organ", organ)
        part = part.reset_index(drop=True)
        rows.append(part)
    return pd.concat(rows, ignore_index=True, sort=True)


def build_decision_log(
    preflight: pd.DataFrame,
    exclusions: pd.DataFrame,
    edges: pd.DataFrame,
    mapping: pd.DataFrame,
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for row in preflight.itertuples(index=False):
        records.append(
            {
                "stage": "source_admission",
                "subject": row.source_id,
                "trigger": row.reason,
                "response": "admitted" if row.status == "admitted" else "quarantined",
                "detail": f"bytes={row.observed_bytes}; sha256={row.observed_sha256}",
            }
        )
    for row in exclusions.itertuples(index=False):
        records.append(
            {
                "stage": "reference_harmonization",
                "subject": f"{row.organ}:{row.source_state}",
                "trigger": row.reason,
                "response": row.response,
                "detail": f"level={row.level}; cells={row.cell_count}; donors={row.donor_count}",
            }
        )
    for row in edges.loc[~edges["eligible"]].itertuples(index=False):
        records.append(
            {
                "stage": "transition_validation",
                "subject": row.edge_id,
                "trigger": row.failure_reasons,
                "response": row.control_response,
                "detail": (
                    f"donor_stability={row.donor_leaveout_stability:.3f}; "
                    f"study_stability={row.study_leaveout_stability}; parameter_stability={row.parameter_stability:.3f}"
                ),
            }
        )
    for row in mapping.loc[mapping["control_response"].eq("abstained")].itertuples(index=False):
        records.append(
            {
                "stage": "sca_signal_mapping",
                "subject": f"{row.edge_id}:{row.evidence_role}",
                "trigger": "fewer than 10 overlapping signature genes",
                "response": "abstained",
                "detail": f"overlap_gene_count={row.overlap_gene_count}",
            }
        )
    return pd.DataFrame(records).sort_values(["stage", "subject"]).reset_index(drop=True)


def build_provenance_manifest(
    root: Path,
    output_dir: Path,
    source_registry_path: Path,
    analysis_config_path: Path,
    references: list[ReferenceData],
) -> dict[str, Any]:
    artifacts = []
    for name in RELEASE_FILES:
        if name == "provenance_manifest.json":
            continue
        path = output_dir / name
        artifacts.append(
            {
                "path": str(path.relative_to(root)),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "parents": ["config/sources.json", "config/analysis.json"],
            }
        )
    cache_artifacts = []
    for reference in references:
        path = root / "data/derived" / f"reference_{reference.organ}.npz"
        cache_artifacts.append(
            {
                "path": str(path.relative_to(root)),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "published": False,
                "rebuild_command": "cardio-sca-mvp run --force-reference",
            }
        )
    return {
        "schema_version": "1.0.0",
        "pipeline_version": __version__,
        "code": [
            {
                "path": str(path.relative_to(root)),
                "sha256": sha256_file(path),
            }
            for path in sorted((root / "src/cardio_sca_mvp").glob("*.py"))
        ],
        "configuration": [
            {
                "path": str(source_registry_path.relative_to(root)),
                "sha256": sha256_file(source_registry_path),
            },
            {
                "path": str(analysis_config_path.relative_to(root)),
                "sha256": sha256_file(analysis_config_path),
            },
        ],
        "published_artifacts": artifacts,
        "local_rebuildable_artifacts": cache_artifacts,
    }


def build_release_manifest(output_dir: Path, validation: dict[str, Any]) -> dict[str, Any]:
    artifacts = []
    for name in RELEASE_FILES:
        path = output_dir / name
        artifacts.append(
            {"path": name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
        )
    return {
        "release_id": "mvp-0.1.0",
        "pipeline_version": __version__,
        "scope": "open-data healthy heart/lung endothelial transition synthesis with indirect SCA-signal mapping",
        "claim_boundary": "hypothesis prioritization only; no patient fate, causal shift, or clinical prediction",
        "validation_status": validation["status"],
        "artifacts": artifacts,
    }


def build_report(
    states: pd.DataFrame,
    edges: pd.DataFrame,
    signals: pd.DataFrame,
    hypotheses: pd.DataFrame,
    validation: dict[str, Any],
) -> str:
    top = hypotheses.head(10)[
        ["rank", "organ", "state_a", "state_b", "priority_score", "adjacency_confidence"]
    ]
    top_table = _markdown_table(top) if not top.empty else "No edge passed eligibility."
    return f"""# MVP result: open-data endothelial transition synthesis

## Outcome

The working MVP assembled {len(states)} healthy endothelial states, evaluated {len(edges)} candidate adjacencies, retained {int(edges['eligible'].sum())}, and mapped {len(signals)} qualified SCA-relevant signature genes. It produced {len(hypotheses)} ranked hypotheses.

These are healthy-state adjacencies with indirect molecular compatibility. They are not observed fate, proof that SCA changes a transition, or clinical predictions. Every edge remains undirected because the MVP admitted no independent adult-organ directionality source.

## Highest-priority hypotheses

{top_table}

## Validation

{validation['passed_count']} of {validation['check_count']} release checks passed. See `validation_report.json` for the machine-readable evidence and `ranked_hypotheses.tsv` for all component scores, uncertainty, counterevidence, and next-data requirements.
"""


def build_follow_up_plan(hypotheses: pd.DataFrame) -> str:
    if hypotheses.empty:
        priorities = "No edge passed the current evidence gates. Acquire an independent healthy endothelial reference before adding SCA data."
    else:
        lines = []
        for row in hypotheses.head(5).itertuples(index=False):
            lines.append(
                f"- **{row.edge_id}:** {row.next_data_requirement}"
            )
        priorities = "\n".join(lines)
    return f"""# Follow-up data-acquisition plan

## Immediate priorities

{priorities}

## Minimum design requirements

- Keep heart and lung analyses organ-specific.
- Preserve donor, sample, treatment, genotype, anatomy, assay, and study identifiers.
- Include enough independent donors to reserve a donor or study holdout before discovery.
- Use direct perturbation, time-resolved sampling, or lineage-compatible measurements before assigning direction.
- Treat matched SCA tissue, longitudinal proof, and clinical endpoints as separate future scope decisions.
"""


def _markdown_table(frame: pd.DataFrame) -> str:
    headers = [str(column) for column in frame.columns]
    rows = []
    for row in frame.itertuples(index=False, name=None):
        rows.append(
            [f"{value:.3f}" if isinstance(value, float) else str(value) for value in row]
        )
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    lines.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(lines)
