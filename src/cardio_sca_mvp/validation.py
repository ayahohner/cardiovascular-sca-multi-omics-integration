from __future__ import annotations

import re
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .ranking import COMPONENT_COLUMNS
from .util import require_columns


FORBIDDEN_OUTPUT_PATTERNS = {
    "causal_sca_shift": re.compile(r"\bSCA (causes|caused|induces|increases|decreases) (a )?transition", re.I),
    "clinical_prediction": re.compile(r"\b(diagnos|prognos|predicts? patient|treatment recommendation)", re.I),
    "observed_fate": re.compile(r"\b(observed|estimated) (cellular )?fate\b", re.I),
}


def validate_analysis(
    preflight: pd.DataFrame,
    state_metadata: pd.DataFrame,
    edges: pd.DataFrame,
    signal_ledger: pd.DataFrame,
    mapping: pd.DataFrame,
    hypotheses: pd.DataFrame,
    analysis_config: dict[str, Any],
) -> dict[str, Any]:
    checks: list[dict[str, Any]] = []

    checks.append(_check("all_sources_admitted", bool((preflight["status"] == "admitted").all()), f"{len(preflight)} source assets"))
    organs = set(state_metadata["organ"])
    checks.append(_check("heart_and_lung_present", organs == {"heart", "lung"}, f"organs={sorted(organs)}"))
    checks.append(
        _check(
            "state_donor_support",
            bool((state_metadata["donor_count"] >= analysis_config["pseudobulk"]["min_donors_per_state"]).all()),
            f"minimum donors={int(state_metadata['donor_count'].min())}",
        )
    )

    eligible = edges.loc[edges["eligible"]]
    checks.append(_check("eligible_edges_exist", not eligible.empty, f"eligible={len(eligible)} total={len(edges)}"))
    checks.append(
        _check(
            "donor_stability_gate",
            bool(
                (
                    eligible["donor_leaveout_stability"]
                    >= analysis_config["adjacency"]["min_donor_leaveout_stability"]
                ).all()
            ),
            f"minimum={eligible['donor_leaveout_stability'].min() if not eligible.empty else 'n/a'}",
        )
    )
    study_values = eligible["study_leaveout_stability"].dropna()
    checks.append(
        _check(
            "study_stability_gate",
            bool(
                study_values.empty
                or (study_values >= analysis_config["adjacency"]["min_study_leaveout_stability"]).all()
            ),
            f"evaluated={len(study_values)}",
        )
    )
    checks.append(
        _check(
            "all_edges_undirected",
            bool(edges["direction_status"].eq("undirected").all()),
            "MVP has no eligible time, velocity, or direct adult-organ perturbation source",
        )
    )

    roles = set(signal_ledger["evidence_role"])
    checks.append(
        _check(
            "both_sca_evidence_roles_present",
            roles == {"endothelial_perturbational_prior", "immune_association"},
            f"roles={sorted(roles)}",
        )
    )
    checks.append(
        _check(
            "mapping_preserves_undirected_status",
            bool(mapping["direction_status"].eq("undirected").all()),
            f"mappings={len(mapping)}",
        )
    )
    mapped = mapping.loc[mapping["control_response"].eq("mapped")]
    checks.append(
        _check(
            "mapped_signals_have_gene_overlap",
            bool(not mapped.empty and (mapped["overlap_gene_count"] >= 10).all()),
            f"mapped={len(mapped)} minimum_overlap={mapped['overlap_gene_count'].min() if not mapped.empty else 'n/a'}",
        )
    )

    require_columns(hypotheses, ["edge_id", *COMPONENT_COLUMNS], "ranked hypotheses")
    component_values = hypotheses[COMPONENT_COLUMNS].to_numpy(dtype=float) if not hypotheses.empty else np.empty((0, len(COMPONENT_COLUMNS)))
    checks.append(
        _check(
            "component_scores_visible_and_bounded",
            bool(
                not hypotheses.empty
                and np.isfinite(component_values).all()
                and ((component_values >= 0) & (component_values <= 1)).all()
            ),
            f"hypotheses={len(hypotheses)} components={len(COMPONENT_COLUMNS)}",
        )
    )
    checks.append(
        _check(
            "ranked_edges_are_eligible",
            set(hypotheses["edge_id"]) <= set(eligible["edge_id"]),
            f"ranked={len(hypotheses)}",
        )
    )
    claim_text = "\n".join(hypotheses.get("permitted_interpretation", pd.Series(dtype=str)).astype(str))
    violations = [name for name, pattern in FORBIDDEN_OUTPUT_PATTERNS.items() if pattern.search(claim_text)]
    checks.append(_check("claim_boundaries", not violations, f"violations={violations}"))

    passed = all(check["passed"] for check in checks)
    return {
        "status": "passed" if passed else "failed",
        "check_count": len(checks),
        "passed_count": sum(check["passed"] for check in checks),
        "failed_count": sum(not check["passed"] for check in checks),
        "checks": checks,
    }


def validate_git_boundary(root: Path) -> dict[str, Any]:
    tracked = subprocess.run(
        ["git", "ls-files", "data/raw", "data/derived"],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    ).stdout.splitlines()
    return _check("raw_and_derived_data_not_tracked", not tracked, f"tracked={tracked}")


def validate_local_markdown_links(root: Path) -> dict[str, Any]:
    markdown_files = [root / "README.md"]
    markdown_files.extend(sorted((root / "docs").glob("*.md")))
    markdown_files.extend(sorted((root / "results/mvp").glob("*.md")))
    missing: list[str] = []
    pattern = re.compile(r"\[[^\]]+\]\(([^)]+)\)")
    for markdown in markdown_files:
        if not markdown.is_file():
            continue
        for target in pattern.findall(markdown.read_text(encoding="utf-8")):
            clean = target.strip("<>").split("#", maxsplit=1)[0]
            if not clean or clean.startswith(("https://", "http://", "mailto:")):
                continue
            if not (markdown.parent / clean).resolve().exists():
                missing.append(f"{markdown.relative_to(root)} -> {target}")
    return _check("local_markdown_links", not missing, f"missing={missing}")


def validate_compute_preflight(
    compute: dict[str, Any], source_total_bytes: int
) -> list[dict[str, Any]]:
    stages = [compute["heart_reference"], compute["lung_reference"], compute["cached_full_pipeline"]]
    peak_rss = max(
        compute["heart_reference"]["peak_rss_bytes"],
        compute["lung_reference"]["peak_rss_bytes"],
    )
    return [
        _check(
            "compute_memory_margin",
            peak_rss < compute["host_memory_bytes"],
            f"peak_rss_bytes={peak_rss}; host_memory_bytes={compute['host_memory_bytes']}",
        ),
        _check(
            "compute_disk_margin",
            source_total_bytes < compute["disk_available_bytes_before_analysis"],
            f"source_bytes={source_total_bytes}; available_bytes={compute['disk_available_bytes_before_analysis']}",
        ),
        _check(
            "compute_completion_and_repeatability",
            all(stage["completion_status"] == "passed" for stage in stages)
            and bool(compute["scientific_artifact_hashes_match_on_repeat"]),
            "heart, lung, cached full pipeline, and repeated scientific hashes",
        ),
    ]


def _check(name: str, passed: bool, detail: str) -> dict[str, Any]:
    return {"name": name, "passed": bool(passed), "detail": detail}
