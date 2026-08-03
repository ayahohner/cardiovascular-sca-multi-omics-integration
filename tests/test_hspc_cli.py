import json
import sys
from pathlib import Path

import pandas as pd
import pytest

from hspc_fate_mvp import cli


def _configure(monkeypatch, tmp_path: Path, arguments: list[str]) -> None:
    monkeypatch.setattr(sys, "argv", ["hspc-fate-mvp", "--root", str(tmp_path), *arguments])
    monkeypatch.setattr(cli, "load_analysis_config", lambda path: {"config": str(path)})


def test_fetch_reports_only_newly_downloaded_assets(monkeypatch, tmp_path: Path, capsys) -> None:
    _configure(monkeypatch, tmp_path, ["fetch"])
    downloaded = tmp_path / "data/raw/hspc-fate/discovery/HD175.zip"
    monkeypatch.setattr(cli, "fetch_assets", lambda config, root: [downloaded])

    cli.main()

    assert json.loads(capsys.readouterr().out) == {"downloaded": [str(downloaded)]}


def test_preflight_exits_nonzero_and_preserves_quarantine_reason(
    monkeypatch, tmp_path: Path, capsys
) -> None:
    _configure(monkeypatch, tmp_path, ["preflight"])
    preflight = pd.DataFrame.from_records(
        [{"local_path": "data/raw/hspc-fate/discovery/HD175.zip", "status": "quarantined", "reason": "checksum_mismatch"}]
    )
    monkeypatch.setattr(cli, "preflight_assets", lambda config, root: preflight)

    with pytest.raises(SystemExit) as stopped:
        cli.main()

    assert stopped.value.code == 1
    output = json.loads(capsys.readouterr().out)
    assert output[0]["status"] == "quarantined"
    assert output[0]["reason"] == "checksum_mismatch"


def test_run_forwards_forced_preprocessing_and_reports_manifest(
    monkeypatch, tmp_path: Path, capsys
) -> None:
    _configure(monkeypatch, tmp_path, ["run", "--force-preprocess"])
    calls = []

    def fake_run(root, force_preprocess=False):
        calls.append((root, force_preprocess))
        return {"release_id": "hspc-fate-v2", "classification": "insufficient_evidence"}

    monkeypatch.setattr(cli, "run_pipeline", fake_run)
    cli.main()

    assert calls == [(tmp_path.resolve(), True)]
    assert json.loads(capsys.readouterr().out)["classification"] == "insufficient_evidence"


def test_verify_exits_nonzero_for_invalid_release(monkeypatch, tmp_path: Path, capsys) -> None:
    _configure(monkeypatch, tmp_path, ["verify"])
    monkeypatch.setattr(
        cli,
        "verify_release",
        lambda root: {"status": "failed", "problems": ["checksum:stability_results.tsv"]},
    )

    with pytest.raises(SystemExit) as stopped:
        cli.main()

    assert stopped.value.code == 1
    assert json.loads(capsys.readouterr().out)["problems"] == [
        "checksum:stability_results.tsv"
    ]
