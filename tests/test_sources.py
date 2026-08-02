import hashlib
import json
from pathlib import Path

import pytest

from cardio_sca_mvp.sources import fetch_missing_sources, load_registry, preflight_sources
import cardio_sca_mvp.sources as source_module


def source_record(path: str, payload: bytes) -> dict:
    return {
        "id": "source-1",
        "role": "healthy_reference",
        "accession": "TEST1",
        "release": "v1",
        "url": "https://example.org/source.bin",
        "local_path": path,
        "expected_bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "format": "test",
        "open_processed": True,
        "terms": "test terms",
        "terms_url": "https://example.org/terms",
        "citation": "test citation",
        "donor_count": 1,
        "sample_count": 1,
        "assay": "test assay",
        "anatomy": "test tissue",
        "covariates": ["donor"],
        "source_preprocessing": "test preprocessing",
        "replication_unit": "donor",
        "limitations": ["synthetic fixture"],
    }


def write_registry(root: Path, record: dict) -> Path:
    path = root / "sources.json"
    path.write_text(
        json.dumps({"schema_version": "1.0.0", "sources": [record]}),
        encoding="utf-8",
    )
    return path


def test_preflight_admits_only_exact_bytes_and_checksum(tmp_path: Path) -> None:
    payload = b"meaningful processed data"
    data_path = tmp_path / "data/source.bin"
    data_path.parent.mkdir()
    data_path.write_bytes(payload)
    registry = load_registry(write_registry(tmp_path, source_record("data/source.bin", payload)))

    admitted = preflight_sources(registry, tmp_path).iloc[0]
    assert admitted["status"] == "admitted"
    assert admitted["reason"] == "verified"

    data_path.write_bytes(b"different bytes but same len")
    mismatched = preflight_sources(registry, tmp_path).iloc[0]
    assert mismatched["status"] == "quarantined"
    assert mismatched["reason"] in {"size_mismatch", "checksum_mismatch"}


def test_registry_rejects_ineligible_or_untraceable_sources(tmp_path: Path) -> None:
    record = source_record("data/source.bin", b"x")
    record["open_processed"] = False
    with pytest.raises(ValueError, match="not eligible"):
        load_registry(write_registry(tmp_path, record))

    record = source_record("data/source.bin", b"x")
    record["terms_url"] = ""
    with pytest.raises(ValueError, match="empty fields"):
        load_registry(write_registry(tmp_path, record))


def test_fetch_replaces_same_size_corruption_only_after_checksum(monkeypatch, tmp_path: Path) -> None:
    payload = b"verified payload"
    record = source_record("data/source.bin", payload)
    registry = {"schema_version": "1.0.0", "sources": [record]}

    class Response:
        def __enter__(self):
            self.buffer = __import__("io").BytesIO(payload)
            return self.buffer

        def __exit__(self, *_):
            self.buffer.close()

    monkeypatch.setattr(source_module.urllib.request, "urlopen", lambda request: Response())
    downloaded = fetch_missing_sources(registry, tmp_path)
    target = tmp_path / record["local_path"]
    assert downloaded == [target]
    assert target.read_bytes() == payload

    target.write_bytes(b"corrupt payload!")
    assert target.stat().st_size == len(payload)
    downloaded = fetch_missing_sources(registry, tmp_path)
    assert downloaded == [target]
    assert target.read_bytes() == payload
