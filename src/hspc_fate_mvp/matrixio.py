from __future__ import annotations

import gzip
import io
import tarfile
import zipfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import sparse
from scipy.io import mmread


@dataclass
class DonorMatrix:
    donor: str
    condition: str
    cohort: str
    counts: sparse.csr_matrix
    genes: np.ndarray
    barcodes: np.ndarray
    gene_ids: np.ndarray | None = None


def _decode_lines(handle: io.BufferedIOBase) -> np.ndarray:
    return np.asarray([line.decode("utf-8").rstrip("\r\n") for line in handle], dtype=str)


def _features(lines: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    identifiers = []
    symbols = []
    for line in lines:
        fields = line.split("\t")
        identifiers.append(fields[0])
        symbols.append(fields[1] if len(fields) > 1 else fields[0])
    return np.asarray(identifiers, dtype=str), np.asarray(symbols, dtype=str)


def read_discovery_zip(path: Path, donor: str, condition: str) -> DonorMatrix:
    with zipfile.ZipFile(path) as archive:
        prefix = f"{donor}/"
        with archive.open(prefix + "features.tsv.gz") as compressed:
            with gzip.GzipFile(fileobj=compressed) as handle:
                gene_ids, genes = _features(_decode_lines(handle))
        with archive.open(prefix + "barcodes.tsv.gz") as compressed:
            with gzip.GzipFile(fileobj=compressed) as handle:
                barcodes = _decode_lines(handle)
        with archive.open(prefix + "matrix.mtx.gz") as compressed:
            with gzip.GzipFile(fileobj=compressed) as handle:
                matrix = mmread(handle, spmatrix=True).tocsr()
    if matrix.shape != (len(genes), len(barcodes)):
        raise ValueError(f"10x dimensions disagree for {donor}: {matrix.shape}")
    return DonorMatrix(donor, condition, "discovery", matrix.T.tocsr(), genes, barcodes, gene_ids)


def read_validation_sample(
    outer_path: Path, archive_name: str, donor: str, condition: str
) -> DonorMatrix:
    with tarfile.open(outer_path, mode="r") as outer:
        member = outer.getmember(archive_name)
        extracted = outer.extractfile(member)
        if extracted is None:
            raise ValueError(f"cannot read nested archive: {archive_name}")
        nested_bytes = extracted.read()
    with tarfile.open(fileobj=io.BytesIO(nested_bytes), mode="r:gz") as inner:
        feature_file = inner.extractfile("GRCh38/genes.tsv")
        barcode_file = inner.extractfile("GRCh38/barcodes.tsv")
        matrix_file = inner.extractfile("GRCh38/matrix.mtx")
        if feature_file is None or barcode_file is None or matrix_file is None:
            raise ValueError(f"incomplete validation archive: {archive_name}")
        gene_ids, genes = _features(_decode_lines(feature_file))
        barcodes = _decode_lines(barcode_file)
        matrix = mmread(matrix_file, spmatrix=True).tocsr()
    if matrix.shape != (len(genes), len(barcodes)):
        raise ValueError(f"10x dimensions disagree for {donor}: {matrix.shape}")
    return DonorMatrix(donor, condition, "validation", matrix.T.tocsr(), genes, barcodes, gene_ids)
