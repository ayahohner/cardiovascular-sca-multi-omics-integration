from __future__ import annotations

import json
import os
import re
import urllib.request
from pathlib import Path
from typing import Any

import pandas as pd

from .util import sha256_file


REQUIRED_SOURCE_FIELDS = {
    "id",
    "role",
    "accession",
    "release",
    "url",
    "local_path",
    "expected_bytes",
    "sha256",
    "format",
    "open_processed",
    "terms",
    "terms_url",
    "citation",
    "donor_count",
    "sample_count",
    "assay",
    "anatomy",
    "covariates",
    "source_preprocessing",
    "replication_unit",
    "limitations",
}


def load_registry(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        registry = json.load(handle)
    if registry.get("schema_version") != "1.0.0":
        raise ValueError("unsupported source-registry schema")
    sources = registry.get("sources")
    if not isinstance(sources, list) or not sources:
        raise ValueError("source registry must contain at least one source")
    identifiers: set[str] = set()
    for source in sources:
        validate_source_record(source)
        if source["id"] in identifiers:
            raise ValueError(f"duplicate source id: {source['id']}")
        identifiers.add(source["id"])
    return registry


def validate_source_record(source: dict[str, Any]) -> None:
    missing = REQUIRED_SOURCE_FIELDS - set(source)
    if missing:
        raise ValueError(f"source record missing fields: {', '.join(sorted(missing))}")
    empty = [
        key
        for key in REQUIRED_SOURCE_FIELDS
        if key not in {"donor_count", "sample_count"}
        and source[key] in (None, "", [], {})
    ]
    if empty:
        raise ValueError(f"source record has empty fields: {', '.join(sorted(empty))}")
    if not source["open_processed"]:
        raise ValueError(f"source {source['id']} is not eligible open processed data")
    if not isinstance(source["expected_bytes"], int) or source["expected_bytes"] <= 0:
        raise ValueError(f"source {source['id']} has invalid expected_bytes")
    if not re.fullmatch(r"[0-9a-f]{64}", source["sha256"]):
        raise ValueError(f"source {source['id']} has invalid sha256")
    if not str(source["url"]).startswith("https://"):
        raise ValueError(f"source {source['id']} must use HTTPS")
    if not str(source["terms_url"]).startswith("https://"):
        raise ValueError(f"source {source['id']} must identify HTTPS terms")
    if not isinstance(source["donor_count"], int) or source["donor_count"] < 0:
        raise ValueError(f"source {source['id']} has invalid donor_count")
    if not isinstance(source["sample_count"], int) or source["sample_count"] < 0:
        raise ValueError(f"source {source['id']} has invalid sample_count")


def source_by_id(registry: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {source["id"]: source for source in registry["sources"]}


def preflight_sources(registry: dict[str, Any], root: Path) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for source in registry["sources"]:
        path = root / source["local_path"]
        record = {
            "source_id": source["id"],
            "role": source["role"],
            "path": source["local_path"],
            "expected_bytes": source["expected_bytes"],
            "observed_bytes": None,
            "expected_sha256": source["sha256"],
            "observed_sha256": None,
            "status": "quarantined",
            "reason": "missing",
        }
        if path.is_file():
            record["observed_bytes"] = path.stat().st_size
            if path.stat().st_size != source["expected_bytes"]:
                record["reason"] = "size_mismatch"
            else:
                digest = sha256_file(path)
                record["observed_sha256"] = digest
                if digest == source["sha256"]:
                    record["status"] = "admitted"
                    record["reason"] = "verified"
                else:
                    record["reason"] = "checksum_mismatch"
        records.append(record)
    return pd.DataFrame.from_records(records)


def assert_all_admitted(preflight: pd.DataFrame) -> None:
    failed = preflight.loc[preflight["status"] != "admitted", ["source_id", "reason"]]
    if not failed.empty:
        summary = ", ".join(f"{row.source_id}:{row.reason}" for row in failed.itertuples())
        raise RuntimeError(f"source admission failed: {summary}")


def fetch_missing_sources(registry: dict[str, Any], root: Path) -> list[Path]:
    downloaded: list[Path] = []
    for source in registry["sources"]:
        destination = root / source["local_path"]
        if destination.is_file() and destination.stat().st_size == source["expected_bytes"]:
            if sha256_file(destination) == source["sha256"]:
                continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        partial = destination.with_suffix(destination.suffix + ".part")
        request = urllib.request.Request(source["url"], headers={"User-Agent": "cardio-sca-mvp/0.1"})
        with urllib.request.urlopen(request) as response, partial.open("wb") as output:
            while chunk := response.read(8 * 1024 * 1024):
                output.write(chunk)
        if partial.stat().st_size != source["expected_bytes"]:
            raise RuntimeError(f"download size mismatch for {source['id']}")
        if sha256_file(partial) != source["sha256"]:
            raise RuntimeError(f"download checksum mismatch for {source['id']}")
        os.replace(partial, destination)
        downloaded.append(destination)
    return downloaded
