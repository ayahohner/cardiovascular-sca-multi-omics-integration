from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse

from cardio_sca_mvp.reference import build_reference, collapse_gene_symbols


def test_duplicate_gene_symbols_are_summed_with_stable_crosswalk() -> None:
    counts = np.array([[1.0, 2.0, 4.0], [3.0, 5.0, 7.0]])
    collapsed, genes, crosswalk = collapse_gene_symbols(
        counts,
        np.array(["G1", "G1", "G2"]),
        np.array(["ENSG1", "ENSG1_ALT", "ENSG2"]),
    )
    assert genes.tolist() == ["G1", "G2"]
    np.testing.assert_array_equal(collapsed, [[3.0, 4.0], [8.0, 7.0]])
    assert crosswalk.loc[crosswalk["gene"].eq("G1"), "feature_ids"].item() == "ENSG1;ENSG1_ALT"


def test_reference_build_uses_donor_state_pseudobulk_and_records_exclusions(tmp_path: Path) -> None:
    donors = [f"D{i}" for i in range(1, 5)]
    rows = []
    counts = []
    for donor_index, donor in enumerate(donors, start=1):
        for state, base in [("A", 10), ("B", 20)]:
            for cell in range(2):
                rows.append({"state": state, "donor": donor, "sample": f"{donor}-{state}"})
                counts.append([base + donor_index, cell + 1, 5, 2])
    for cell in range(2):
        rows.append({"state": "X", "donor": "D1", "sample": "D1-X"})
        counts.append([1, 1, 20, 3])
    obs = pd.DataFrame(rows, index=[f"cell-{i}" for i in range(len(rows))])
    var = pd.DataFrame(
        {"feature_name": ["G1", "G1", "G2", "MT-X"]},
        index=["ENSG1", "ENSG1_ALT", "ENSG2", "ENSGMT"],
    )
    fixture = ad.AnnData(X=sparse.csr_matrix(np.asarray(counts)), obs=obs, var=var)
    fixture.raw = fixture
    path = tmp_path / "reference.h5ad"
    fixture.write_h5ad(path)

    result = build_reference(
        path,
        {
            "organ": "heart",
            "state_column": "state",
            "donor_column": "donor",
            "sample_column": "sample",
            "study_column": None,
            "state_map": {"A": "state_a", "B": "state_b"},
            "excluded_states": {"X": "noncanonical fixture state"},
        },
        {
            "batch_cells": 3,
            "min_cells_per_donor_state": 2,
            "min_donors_per_state": 3,
            "normalization_target": 1_000_000.0,
        },
    )
    assert result.expression.shape == (8, 3)
    assert result.groups.groupby("state")["donor"].nunique().to_dict() == {"state_a": 4, "state_b": 4}
    assert set(result.expression.columns) == {"G1", "G2", "MT-X"}
    assert np.isfinite(result.expression.to_numpy()).all()
    exclusion = result.exclusions.loc[result.exclusions["source_state"].eq("X")].iloc[0]
    assert exclusion["response"] == "excluded"
    assert exclusion["cell_count"] == 2
