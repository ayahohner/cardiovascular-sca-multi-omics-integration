from __future__ import annotations

import hashlib
import gzip
import json
import os
import tempfile
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


def sha256_file(path: Path, chunk_bytes: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_bytes):
            digest.update(chunk)
    return digest.hexdigest()


def canonical_json_bytes(value: Any) -> bytes:
    return (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False) + "\n").encode()


def atomic_write_bytes(path: Path, payload: bytes) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        handle.write(payload)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def write_json(path: Path, value: Any) -> None:
    atomic_write_bytes(path, canonical_json_bytes(value))


def write_tsv(path: Path, frame: pd.DataFrame) -> None:
    payload = frame.to_csv(sep="\t", index=False, lineterminator="\n").encode()
    atomic_write_bytes(path, payload)


def write_gzip_tsv(path: Path, frame: pd.DataFrame) -> None:
    payload = frame.to_csv(sep="\t", index=False, lineterminator="\n").encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
    try:
        with temporary.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
                compressed.write(payload)
        os.replace(temporary, path)
    finally:
        if temporary.exists():
            temporary.unlink()


def require_columns(frame: pd.DataFrame, columns: Iterable[str], context: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{context} is missing columns: {', '.join(missing)}")


def minmax_score(values: pd.Series) -> pd.Series:
    numeric = values.astype(float)
    low = numeric.min()
    high = numeric.max()
    if high == low:
        return pd.Series(1.0, index=values.index)
    return (numeric - low) / (high - low)
