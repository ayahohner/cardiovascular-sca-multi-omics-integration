import gzip
import io
import json
import tarfile
import zipfile
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse
from scipy.io import mmwrite

from hspc_fate_mvp import preprocess
from hspc_fate_mvp.matrixio import DonorMatrix, read_discovery_zip, read_validation_sample
from hspc_fate_mvp.preprocess import _scrublet_calls, qc_donor
from hspc_fate_mvp.sources import load_analysis_config


ROOT = Path(__file__).resolve().parents[1]


def _matrix_market(matrix: sparse.csr_matrix) -> bytes:
    buffer = io.BytesIO()
    mmwrite(buffer, matrix)
    return buffer.getvalue()


def _add_tar_member(archive: tarfile.TarFile, name: str, payload: bytes) -> None:
    member = tarfile.TarInfo(name)
    member.size = len(payload)
    archive.addfile(member, io.BytesIO(payload))


def test_real_source_formats_round_trip_without_reordering(tmp_path: Path) -> None:
    genes = b"ENSG1\tCD34\nENSG2\tMPO\nENSG3\tGATA1\n"
    barcodes = b"cell-a\ncell-b\n"
    gene_by_cell = sparse.csr_matrix(np.asarray([[5, 0], [1, 3], [0, 7]], dtype=int))

    discovery = tmp_path / "HD175.zip"
    with zipfile.ZipFile(discovery, "w") as archive:
        archive.writestr("HD175/features.tsv.gz", gzip.compress(genes))
        archive.writestr("HD175/barcodes.tsv.gz", gzip.compress(barcodes))
        archive.writestr("HD175/matrix.mtx.gz", gzip.compress(_matrix_market(gene_by_cell)))
    observed_discovery = read_discovery_zip(discovery, "HD175", "healthy")

    nested_buffer = io.BytesIO()
    with tarfile.open(fileobj=nested_buffer, mode="w:gz") as nested:
        _add_tar_member(nested, "GRCh38/genes.tsv", genes)
        _add_tar_member(nested, "GRCh38/barcodes.tsv", barcodes)
        _add_tar_member(nested, "GRCh38/matrix.mtx", _matrix_market(gene_by_cell))
    validation = tmp_path / "GSE133181_RAW.tar"
    with tarfile.open(validation, mode="w") as outer:
        _add_tar_member(outer, "GSM3909256_BM1_filtered_gene_bc_matrices.tar.gz", nested_buffer.getvalue())
    observed_validation = read_validation_sample(
        validation, "GSM3909256_BM1_filtered_gene_bc_matrices.tar.gz", "BM1", "healthy"
    )

    for observed, cohort in ((observed_discovery, "discovery"), (observed_validation, "validation")):
        assert observed.cohort == cohort
        assert observed.genes.tolist() == ["CD34", "MPO", "GATA1"]
        assert observed.barcodes.tolist() == ["cell-a", "cell-b"]
        np.testing.assert_array_equal(observed.counts.toarray(), gene_by_cell.T.toarray())


def test_qc_ledger_accounts_for_each_exclusion_and_doublets() -> None:
    donor = DonorMatrix(
        donor="HD-test",
        condition="healthy",
        cohort="discovery",
        genes=np.asarray(["CD34", "MT-ND1", "A", "B", "C"]),
        barcodes=np.asarray([f"cell-{index}" for index in range(5)]),
        counts=sparse.csr_matrix(
            np.asarray(
                [
                    [1000, 0, 1000, 1000, 1000],
                    [250, 0, 250, 250, 250],
                    [1000, 1000, 1000, 1000, 1000],
                    [0, 0, 1500, 1500, 1500],
                    [1000, 0, 1000, 1000, 1000],
                ]
            )
        ),
    )
    rules = {
        "minimum_umis": 3000,
        "minimum_genes": 3,
        "maximum_mitochondrial_fraction": 0.15,
        "require_cd34_expression": True,
        "scrublet_expected_doublet_rate": 0.1,
    }

    def fake_doublets(counts, expected_rate, seed):
        assert counts.shape == (2, 5)
        assert expected_rate == 0.1
        assert seed == 17
        return np.asarray([0.1, 0.9]), np.asarray([False, True]), 0.5, "test"

    retained, ledger, diagnostics = qc_donor(donor, rules, seed=17, doublet_caller=fake_doublets)
    exclusions = {row["control"]: row["excluded_cells"] for row in ledger}

    assert retained.barcodes.tolist() == ["cell-0"]
    assert exclusions == {
        "raw_barcodes": 0,
        "nonzero": 0,
        "minimum_umis": 1,
        "minimum_genes": 0,
        "maximum_mitochondrial_fraction": 1,
        "cd34_expression": 1,
        "scrublet_doublet": 1,
    }
    assert sum(exclusions.values()) == donor.counts.shape[0] - retained.counts.shape[0]
    assert diagnostics["retained_cells"] == 1


def test_degenerate_scrublet_uses_auditable_complexity_fallback(monkeypatch) -> None:
    counts = np.zeros((100, 100), dtype=int)
    counts[:99, :10] = 10
    counts[99, :] = 100

    class DegenerateScrublet:
        threshold_ = None

        def __init__(self, matrix, **kwargs):
            assert matrix.shape == (100, 100)

        def scrub_doublets(self, **kwargs):
            return np.zeros(100), None

    monkeypatch.setattr(preprocess.scr, "Scrublet", DegenerateScrublet)
    scores, calls, threshold, mode = _scrublet_calls(sparse.csr_matrix(counts), 0.1, 13)

    assert mode == "robust_joint_umi_gene_outlier"
    assert np.isnan(threshold)
    assert np.isfinite(scores).all()
    assert np.flatnonzero(calls).tolist() == [99]


def test_analysis_config_rejects_assets_outside_module_namespace(tmp_path: Path) -> None:
    config = json.loads((ROOT / "research/hspc-fate/analysis.json").read_text())
    config["required_assets"][0]["local_path"] = "data/raw/endothelial/HD175.zip"
    path = tmp_path / "analysis.json"
    path.write_text(json.dumps(config), encoding="utf-8")

    with pytest.raises(ValueError, match="escapes HSPC namespace"):
        load_analysis_config(path)
