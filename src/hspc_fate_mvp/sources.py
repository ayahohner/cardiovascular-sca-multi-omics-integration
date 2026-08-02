from __future__ import annotations

import json
import os
import re
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from .util import sha256_file


ASSET_FIELDS = {
    "source_id",
    "cohort",
    "donor",
    "condition",
    "url",
    "local_path",
    "expected_bytes",
    "sha256",
}


def load_analysis_config(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        config = json.load(handle)
    if config.get("schema_version") != "1.0.0":
        raise ValueError("unsupported HSPC analysis schema")
    assets = config.get("required_assets")
    if not isinstance(assets, list) or len(assets) != 7:
        raise ValueError("HSPC analysis requires six discovery archives and one validation archive")
    seen_paths: set[str] = set()
    for asset in assets:
        missing = ASSET_FIELDS - set(asset)
        if missing:
            raise ValueError(f"asset missing fields: {', '.join(sorted(missing))}")
        if asset["local_path"] in seen_paths:
            raise ValueError(f"duplicate asset path: {asset['local_path']}")
        seen_paths.add(asset["local_path"])
        if not str(asset["url"]).startswith("https://"):
            raise ValueError(f"asset URL is not HTTPS: {asset['source_id']}")
        if not str(asset["local_path"]).startswith("data/raw/hspc-fate/"):
            raise ValueError(f"asset escapes HSPC namespace: {asset['local_path']}")
        if not isinstance(asset["expected_bytes"], int) or asset["expected_bytes"] <= 0:
            raise ValueError(f"invalid asset size: {asset['source_id']}")
        if not re.fullmatch(r"[0-9a-f]{64}", asset["sha256"]):
            raise ValueError(f"invalid asset checksum: {asset['source_id']}")
    return config


def load_source_ledger(path: Path) -> dict[str, Any]:
    with path.open(encoding="utf-8") as handle:
        ledger = json.load(handle)
    if ledger.get("module") != "hspc-fate":
        raise ValueError("source ledger belongs to another module")
    sources = ledger.get("sources")
    if not isinstance(sources, list) or not sources:
        raise ValueError("source ledger is empty")
    identifiers = [source.get("id") for source in sources]
    if None in identifiers or len(set(identifiers)) != len(identifiers):
        raise ValueError("source ledger identifiers must be present and unique")
    return ledger


def preflight_assets(config: dict[str, Any], root: Path) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    for asset in config["required_assets"]:
        path = root / asset["local_path"]
        record = {
            **{key: asset[key] for key in ("source_id", "cohort", "donor", "condition", "local_path")},
            "expected_bytes": asset["expected_bytes"],
            "observed_bytes": None,
            "expected_sha256": asset["sha256"],
            "observed_sha256": None,
            "status": "quarantined",
            "reason": "missing",
        }
        if path.is_file():
            record["observed_bytes"] = path.stat().st_size
            record["retrieval_time_utc"] = datetime.fromtimestamp(
                path.stat().st_mtime, timezone.utc
            ).isoformat()
            if path.stat().st_size != asset["expected_bytes"]:
                record["reason"] = "size_mismatch"
            else:
                digest = sha256_file(path)
                record["observed_sha256"] = digest
                if digest == asset["sha256"]:
                    record["status"] = "admitted"
                    record["reason"] = "verified"
                else:
                    record["reason"] = "checksum_mismatch"
        records.append(record)
    return pd.DataFrame.from_records(records)


def assert_all_admitted(preflight: pd.DataFrame) -> None:
    failed = preflight.loc[preflight["status"] != "admitted", ["local_path", "reason"]]
    if not failed.empty:
        detail = ", ".join(f"{row.local_path}:{row.reason}" for row in failed.itertuples())
        raise RuntimeError(f"HSPC source admission failed: {detail}")


def fetch_assets(config: dict[str, Any], root: Path) -> list[Path]:
    downloaded: list[Path] = []
    for asset in config["required_assets"]:
        destination = root / asset["local_path"]
        if destination.is_file() and destination.stat().st_size == asset["expected_bytes"]:
            if sha256_file(destination) == asset["sha256"]:
                continue
        destination.parent.mkdir(parents=True, exist_ok=True)
        partial = destination.with_suffix(destination.suffix + ".part")
        request = urllib.request.Request(asset["url"], headers={"User-Agent": "hspc-fate-mvp/0.1"})
        with urllib.request.urlopen(request) as response, partial.open("wb") as output:
            while chunk := response.read(8 * 1024 * 1024):
                output.write(chunk)
        if partial.stat().st_size != asset["expected_bytes"]:
            raise RuntimeError(f"download size mismatch for {asset['local_path']}")
        if sha256_file(partial) != asset["sha256"]:
            raise RuntimeError(f"download checksum mismatch for {asset['local_path']}")
        os.replace(partial, destination)
        downloaded.append(destination)
    return downloaded


def build_donor_manifest(config: dict[str, Any]) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for asset in config["required_assets"]:
        if asset["cohort"] != "discovery":
            continue
        rows.append(
            {
                "source_id": asset["source_id"],
                "cohort": "discovery",
                "donor": asset["donor"],
                "condition": asset["condition"],
                "assay": "10x_scRNA_seq",
                "tissue": "plerixafor_mobilized_CD34_positive_HSPCs",
                "included": True,
                "decision": "required_discovery",
                "disease_genotype": "SCD" if asset["condition"] == "SCD" else "not_applicable",
                "ancestry": "not_reported",
                "geography": "not_reported_at_donor_level",
                "environment": "not_reported",
                "treatment_context": "plerixafor_mobilization_and_gene_therapy_collection",
                "access_to_care": "not_reported",
                "self_identified_race": "not_reported",
            }
        )
    for sample in config["validation_samples"]:
        rows.append(
            {
                "source_id": "gse133181-bm-cd34-scrna",
                "cohort": "validation",
                "donor": sample["donor"],
                "condition": sample["condition"],
                "assay": "10x_scRNA_seq",
                "tissue": "bone_marrow_CD34_positive_cells",
                "included": bool(sample["include"]),
                "decision": "required_validation" if sample["include"] else "predeclared_exclusion",
                "disease_genotype": "SCD" if sample["condition"] == "SCD" else "not_applicable",
                "ancestry": "not_reported",
                "geography": "not_reported_at_donor_level",
                "environment": "not_reported",
                "treatment_context": "pretransplant_bone_marrow_collection",
                "access_to_care": "not_reported",
                "self_identified_race": "not_reported",
            }
        )
    return pd.DataFrame.from_records(rows)


def build_source_manifest(
    config: dict[str, Any], ledger: dict[str, Any], preflight: pd.DataFrame
) -> dict[str, Any]:
    source_lookup = {source["id"]: source for source in ledger["sources"]}
    admitted_ids = {asset["source_id"] for asset in config["required_assets"]}
    asset_lookup = {asset["local_path"]: asset for asset in config["required_assets"]}
    tissue_lookup = {
        "s-bsst1258-hspc-scrna": "plerixafor-mobilized CD34-positive HSPCs",
        "gse133181-bm-cd34-scrna": "bone-marrow CD34-positive cells",
    }
    assets = []
    for record in preflight.to_dict(orient="records"):
        configured = asset_lookup[record["local_path"]]
        source = source_lookup[record["source_id"]]
        assets.append(
            {
                **record,
                "resolved_url": configured["url"],
                "terms": source["license"],
                "donor_label": configured["donor"],
                "assay": source["assay"],
                "tissue": tissue_lookup[record["source_id"]],
                "preprocessing_state": source["format"],
                "limitations": source["limitations"],
            }
        )
    return {
        "schema_version": "1.0.0",
        "module": "hspc-fate",
        "required_source_ids": sorted(admitted_ids),
        "sources": [source_lookup[source_id] for source_id in sorted(admitted_ids)],
        "assets": assets,
    }
