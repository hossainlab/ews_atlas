"""The same invariants, re-checked on the built atlas object when one exists.

Skipped automatically on a clean checkout. These read the object in backed mode
where possible so the suite stays runnable on a laptop.
"""

from __future__ import annotations

import numpy as np
import pytest

pytestmark = pytest.mark.slow


def _read(path):
    import anndata as ad

    return ad.read_h5ad(path)


def test_atlas_has_counts_layer_and_no_empty_cells(atlas_path):
    A = _read(atlas_path)
    assert "counts" in A.layers, "integration models need the raw count layer"
    totals = np.asarray(A.layers["counts"].sum(axis=1)).ravel()
    assert (totals > 0).all(), f"{int((totals == 0).sum())} cells carry zero counts"
    d = A.layers["counts"].data
    assert np.allclose(d, np.round(d))


def test_atlas_markers_present(atlas_path):
    from ewsatlas import genespace

    A = _read(atlas_path)
    missing = [g for g, ok in genespace.check_markers(A.var_names).items() if not ok]
    assert not missing, f"markers absent: {missing}"


def test_atlas_covariates_complete(atlas_path):
    A = _read(atlas_path)
    for col in ("study", "sample", "patient", "material_type", "assay"):
        assert col in A.obs.columns, f"missing covariate column {col}"
        assert A.obs[col].isna().sum() == 0, f"{col} has NaN; use the not_reported sentinel"


def test_atlas_mito_nonzero_in_whole_cell_studies(atlas_path):
    A = _read(atlas_path)
    if "pct_counts_mt" not in A.obs or "assay" not in A.obs:
        pytest.skip("QC metrics not stored on this object")
    for assay, sub in A.obs.groupby("assay", observed=True):
        if "sn" in str(assay):
            continue
        assert sub["pct_counts_mt"].median() > 0, f"{assay}: zero mito fraction"
