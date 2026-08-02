import json
from pathlib import Path

import pytest

from cardio_sca_mvp.reference import load_or_build_reference
from cardio_sca_mvp.pipeline import run_pipeline, verify_release
from cardio_sca_mvp.signals import build_signal_ledger
from cardio_sca_mvp.sources import load_registry, source_by_id


ROOT = Path(__file__).resolve().parents[1]
DATA_AVAILABLE = (ROOT / "data/raw/heart/heart_vascular.h5ad").is_file() and (
    ROOT / "data/raw/lung/hlca_core.h5ad"
).is_file()


@pytest.mark.data
@pytest.mark.skipif(not DATA_AVAILABLE, reason="admitted public source files are not present")
def test_admitted_real_sources_produce_expected_evidence_units() -> None:
    registry = load_registry(ROOT / "config/sources.json")
    sources = source_by_id(registry)
    analysis = json.loads((ROOT / "config/analysis.json").read_text(encoding="utf-8"))
    ledger = build_signal_ledger(ROOT, sources, analysis["signatures"])
    assert set(ledger["evidence_role"]) == {
        "endothelial_perturbational_prior",
        "immune_association",
    }
    assert ledger.groupby("evidence_role").size().to_dict() == {
        "endothelial_perturbational_prior": 300,
        "immune_association": 300,
    }
    assert ledger.loc[ledger["evidence_role"].eq("immune_association"), "fdr"].notna().all()

    states = []
    for spec in analysis["reference_sources"]:
        reference = load_or_build_reference(
            ROOT,
            sources[spec["source_id"]],
            spec,
            analysis["pseudobulk"],
            ROOT / "data/derived",
        )
        states.append(reference.states)
    combined = __import__("pandas").concat(states, ignore_index=True)
    assert set(combined["organ"]) == {"heart", "lung"}
    assert len(combined) == 14
    assert combined["donor_count"].min() >= 3


@pytest.mark.data
@pytest.mark.skipif(not DATA_AVAILABLE, reason="admitted public source files are not present")
def test_real_pipeline_release_is_valid_and_deterministic() -> None:
    first = run_pipeline(ROOT)
    first_hashes = {artifact["path"]: artifact["sha256"] for artifact in first["artifacts"]}
    assert first["validation_status"] == "passed"

    second = run_pipeline(ROOT)
    second_hashes = {artifact["path"]: artifact["sha256"] for artifact in second["artifacts"]}
    assert first_hashes == second_hashes

    verification = verify_release(ROOT)
    assert verification["status"] == "passed"
    assert all(check["passed"] for check in verification["artifact_checks"])
