import json
import sys
from pathlib import Path

import pandas as pd
import pytest

import cardio_sca_mvp.cli as cli


def invoke(monkeypatch, root: Path, *arguments: str) -> None:
    monkeypatch.setattr(sys, "argv", ["cardio-sca-mvp", "--root", str(root), *arguments])
    cli.main()


def test_cli_routes_fetch_preflight_run_and_verify(monkeypatch, tmp_path: Path, capsys) -> None:
    registry = {"schema_version": "1.0.0", "sources": []}
    monkeypatch.setattr(cli, "load_registry", lambda path: registry)
    monkeypatch.setattr(cli, "fetch_missing_sources", lambda value, root: [root / "downloaded.bin"])
    invoke(monkeypatch, tmp_path, "fetch")
    assert json.loads(capsys.readouterr().out)["downloaded"] == [str(tmp_path / "downloaded.bin")]

    admitted = pd.DataFrame({"source_id": ["s1"], "status": ["admitted"]})
    monkeypatch.setattr(cli, "preflight_sources", lambda value, root: admitted)
    invoke(monkeypatch, tmp_path, "preflight")
    assert json.loads(capsys.readouterr().out)[0]["status"] == "admitted"

    monkeypatch.setattr(cli, "run_pipeline", lambda root, force_reference: {"release_id": "mvp-test", "forced": force_reference})
    invoke(monkeypatch, tmp_path, "run", "--force-reference")
    assert json.loads(capsys.readouterr().out) == {"release_id": "mvp-test", "forced": True}

    monkeypatch.setattr(cli, "verify_release", lambda root: {"status": "passed"})
    invoke(monkeypatch, tmp_path, "verify")
    assert json.loads(capsys.readouterr().out)["status"] == "passed"


def test_cli_returns_failure_for_quarantine_or_invalid_release(monkeypatch, tmp_path: Path) -> None:
    monkeypatch.setattr(cli, "load_registry", lambda path: {"sources": []})
    quarantined = pd.DataFrame({"source_id": ["s1"], "status": ["quarantined"]})
    monkeypatch.setattr(cli, "preflight_sources", lambda value, root: quarantined)
    with pytest.raises(SystemExit) as preflight_exit:
        invoke(monkeypatch, tmp_path, "preflight")
    assert preflight_exit.value.code == 1

    monkeypatch.setattr(cli, "verify_release", lambda root: {"status": "failed"})
    with pytest.raises(SystemExit) as verify_exit:
        invoke(monkeypatch, tmp_path, "verify")
    assert verify_exit.value.code == 1
