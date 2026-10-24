"""Shared fixtures.

The synthetic fixture reproduces the exact shape of the Aynaud deposit defect:
two samples whose feature files use different gene namings. Every invariant
test runs against it, so the suite is fast and needs no raw data, while
``test_atlas_object.py`` re-checks the same invariants on the real object when
one is present.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

GENES = ["CD99", "FLI1", "NKX2-2", "NR0B1", "PTPRC", "CD3E", "LYZ", "COL1A1",
         "PECAM1", "MKI67", "TOP2A", "EWSR1", "MT-CO1", "MT-ND1", "ACTB", "GAPDH"]


@pytest.fixture
def two_naming_studies():
    """One study with CellRanger v2 ``GRCh38_`` symbols, one with plain v3.

    Returns ``{study: AnnData}`` with raw integer counts in ``X``.
    """
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp

    rng = np.random.default_rng(0)
    out = {}
    for study, prefix, n in (("prefixed", "GRCh38_", 400), ("plain", "", 300)):
        X = sp.csr_matrix(rng.poisson(1.5, size=(n, len(GENES))).astype(np.float32))
        a = ad.AnnData(
            X=X,
            obs=pd.DataFrame({"sample": [f"{study}_S{i % 2}" for i in range(n)],
                              "study": study},
                             index=[f"{study}|{i}" for i in range(n)]),
            var=pd.DataFrame(index=pd.Index([prefix + g for g in GENES], name="gene")),
        )
        out[study] = a
    return out


@pytest.fixture
def harmonized(two_naming_studies):
    """The same pair after the reader's symbol harmonization."""
    from ewsatlas import io

    for a in two_naming_studies.values():
        a.var_names = io.harmonize_symbols(a.var_names)
        a.var_names_make_unique()
    return two_naming_studies


@pytest.fixture(scope="session")
def atlas_path():
    """The most advanced atlas object present.

    Order: the configured final object, then the pre-integration object from
    stage 02, then any EWSAtlas_v*.h5ad in data/processed. Skips on a clean
    checkout.
    """
    from ewsatlas import config

    cfg = config.load()
    final = config.ROOT / config.get(cfg, "paths.atlas_out")
    if final.is_file():
        return final
    interim = config.ROOT / config.get(cfg, "paths.interim") / "atlas_preintegration.h5ad"
    if interim.is_file():
        return interim
    processed = config.ROOT / config.get(cfg, "paths.processed")
    found = sorted(processed.glob("EWSAtlas_v*.h5ad")) if processed.is_dir() else []
    if found:
        return found[-1]
    pytest.skip("no built atlas object present")
