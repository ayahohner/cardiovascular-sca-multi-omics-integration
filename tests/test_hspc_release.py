import json
import shutil
from pathlib import Path

import pandas as pd
import pytest

from hspc_fate_mvp.pipeline import run_pipeline, verify_release
from hspc_fate_mvp.util import sha256_file


ROOT = Path(__file__).resolve().parents[1]
RESULT = ROOT / "results/hspc-fate"
RAW_AVAILABLE = all(
    (ROOT / relative).is_file()
    for relative in (
        "data/raw/hspc-fate/discovery/HD175.zip",
        "data/raw/hspc-fate/discovery/HD177.zip",
        "data/raw/hspc-fate/discovery/HD192.zip",
        "data/raw/hspc-fate/discovery/SCD2.zip",
        "data/raw/hspc-fate/discovery/SCD3.zip",
        "data/raw/hspc-fate/discovery/SCD4.zip",
        "data/raw/hspc-fate/validation/GSE133181_RAW.tar",
    )
)


def _release_copy(tmp_path: Path) -> Path:
    destination = tmp_path / "results/hspc-fate"
    destination.parent.mkdir(parents=True)
    shutil.copytree(RESULT, destination)
    return destination


def test_release_verifier_detects_file_tampering(tmp_path: Path) -> None:
    release = _release_copy(tmp_path)
    with (release / "claim_boundary.md").open("a", encoding="utf-8") as handle:
        handle.write("tamper\n")

    verification = verify_release(tmp_path)
    assert verification["status"] == "failed"
    assert any(problem.endswith(":claim_boundary.md") for problem in verification["problems"])


def test_release_verifier_detects_semantic_probability_corruption(tmp_path: Path) -> None:
    release = _release_copy(tmp_path)
    path = release / "donor_lineage_priming.tsv"
    frame = pd.read_csv(path, sep="\t")
    frame.loc[0, "mean_lineage_priming_probability"] += 0.2
    frame.to_csv(path, sep="\t", index=False)

    manifest_path = release / "release_manifest.json"
    manifest = json.loads(manifest_path.read_text())
    artifact = next(item for item in manifest["artifacts"] if item["path"] == path.name)
    artifact["bytes"] = path.stat().st_size
    artifact["sha256"] = sha256_file(path)
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    verification = verify_release(tmp_path)
    assert verification["status"] == "failed"
    assert "donor_lineage_priming_not_normalized" in verification["problems"]


@pytest.mark.data
@pytest.mark.skipif(not RAW_AVAILABLE, reason="admitted public HSPC source files are not present")
def test_real_hspc_pipeline_is_reproducible_and_returns_a_valid_gated_result() -> None:
    first = run_pipeline(ROOT)
    second = run_pipeline(ROOT)

    assert first["classification"] == "insufficient_evidence"
    assert {item["path"]: item["sha256"] for item in first["artifacts"]} == {
        item["path"]: item["sha256"] for item in second["artifacts"]
    }
    assert verify_release(ROOT) == {"status": "passed", "problems": []}

    gates = pd.read_csv(RESULT / "gate_ledger.tsv", sep="\t")
    failed = gates.loc[gates["status"].eq("failed"), "gate"].tolist()
    assert failed == [
        "doublet_call_rate_audit",
        "model_robustness_execution",
        "discovery_topology",
        "validation_topology",
        "validation_root_cell_reliability",
    ]
    result = json.loads((RESULT / "result_classification.json").read_text())
    assert result["classification"] == "insufficient_evidence"
    assert result["discovery_stable_lineages"] == []
    assert result["validation_estimates_are_exploratory"] is True

    stability = pd.read_csv(RESULT / "stability_results.tsv", sep="\t")
    assert not stability["model_robustness_passed"].any()
    for row in stability.itertuples():
        expected = 20 if row.cohort == "discovery" else 126
        assert row.exact_permutations == expected
        assert len(json.loads(row.palantir_permutation_assignments_json)) == expected
        assert len(json.loads(row.cellrank_cytotrace_permutation_assignments_json)) == expected
        assert len(json.loads(row.palantir_leave_one_out_effects_json)) == row.leave_one_out_total
        assert len(json.loads(row.cellrank_cytotrace_leave_one_out_effects_json)) == row.leave_one_out_total

    robustness = pd.read_csv(RESULT / "model_robustness.tsv", sep="\t")
    assert {
        "baseline",
        "donor_deletion_refit",
        "donor_deletion_fixed_reference",
        "hyperparameter_refit",
    }.issubset(set(robustness["perturbation"]))
    assert robustness.loc[
        robustness["status"].eq("failed"), "detail"
    ].str.len().gt(0).all()
