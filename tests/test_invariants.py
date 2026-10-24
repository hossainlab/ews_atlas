"""Invariants that each correspond to a defect actually observed in this project.

Every test names the failure it guards against; a test here is not a style
check, it is a tripwire on a bug that already cost a rebuild.
"""

from __future__ import annotations

import numpy as np
import pytest


# --- defect 1: GSE130024 ships two incompatible gene namings -----------------

def test_symbol_harmonization_strips_build_prefix():
    from ewsatlas import io

    assert io.harmonize_symbols(["GRCh38_CD99", "CD99", "hg19-MT-CO1"]) == ["CD99", "CD99", "MT-CO1"]
    assert io.harmonize_symbols(["ENSG00000002586__CD99"]) == ["CD99"]


def test_unharmonized_intersection_is_empty(two_naming_studies):
    """The bug's signature: without harmonization the two studies share no genes."""
    a, b = two_naming_studies["prefixed"], two_naming_studies["plain"]
    assert len(set(a.var_names) & set(b.var_names)) == 0


def test_harmonized_gene_space_is_complete(harmonized):
    from ewsatlas import genespace

    genes, report = genespace.build(harmonized, min_cells=3)
    assert report["n_annotation_intersection"] == len(harmonized["plain"].var_names)
    assert len(genes) == report["n_annotation_intersection"]


# --- defect 2: 13,064 of 15,057 cells silently zeroed by feature selection ---

def test_no_zero_count_cells_after_subsetting(harmonized):
    from ewsatlas import genespace

    genes, _ = genespace.build(harmonized, min_cells=3)
    for name, a in harmonized.items():
        sub = a[:, genes]
        totals = np.asarray(sub.X.sum(axis=1)).ravel()
        assert (totals > 0).all(), f"{name}: {int((totals == 0).sum())} cells zeroed"


def test_per_study_counts_reconcile(harmonized):
    from ewsatlas import genespace
    import anndata as ad

    genes, _ = genespace.build(harmonized, min_cells=3)
    parts = [a[:, genes].copy() for a in harmonized.values()]
    Z = ad.concat(parts, join="outer")
    expected = {k: v.n_obs for k, v in harmonized.items()}
    assert Z.obs["study"].value_counts().to_dict() == expected


# --- defect 3: marker panel lost, making annotation impossible ---------------

def test_marker_panel_survives(harmonized):
    from ewsatlas import genespace

    genes, _ = genespace.build(harmonized, min_cells=3)
    missing = [g for g, ok in genespace.check_markers(genes).items() if not ok]
    assert not missing, f"markers absent from feature space: {missing}"


# --- defect 4: 0.0% mitochondrial fraction because MT- was unrecognized ------

def test_mito_fraction_nonzero_for_whole_cell(harmonized):
    from ewsatlas import qc

    a = harmonized["prefixed"]
    qc.add_qc_metrics(a)
    assert a.var["mt"].sum() == 2, "MT- genes not recognized"
    assert a.obs["pct_counts_mt"].median() > 0, (
        "zero mito fraction across a droplet whole-cell study indicates the "
        "gene symbols were never matched")


# --- defect 5: non-integral values reaching a count-likelihood model ---------

def test_counts_layer_is_integral(harmonized):
    for name, a in harmonized.items():
        d = a.X.data if hasattr(a.X, "data") else np.asarray(a.X).ravel()
        assert np.allclose(d, np.round(d)), f"{name}: counts layer is not integral"


# --- defect 6: covariates merged from two incompatible schemas --------------

def test_covariate_table_covers_every_sample_key():
    import pandas as pd
    from ewsatlas import config

    cfg = config.load()
    p = config.ROOT / config.get(cfg, "paths.covariates")
    if not p.is_file():
        pytest.skip("covariate table not present")
    cov = pd.read_csv(p)
    for col in ("study", "sample"):
        assert col in cov.columns
    keys = cov["study"].astype(str) + "||" + cov["sample"].astype(str)
    assert keys.is_unique, "duplicate study||sample keys would multiply cells on join"
    assert cov.isna().sum().sum() == 0, "blanks must be the explicit not_reported sentinel"


# --- QC policy invariants ---------------------------------------------------

def test_mito_cap_is_assay_aware():
    """Nuclei are depleted of mitochondrial RNA, so their cap must be STRICTER."""
    from ewsatlas import config

    caps = config.get(config.load(), "qc.mito_cap_pct_by_assay")
    assert caps["10x_snRNA"] < caps["10x_scRNA"] <= caps["CELseq2_scRNA"], caps


def test_gene_floor_is_assay_aware():
    from ewsatlas import config

    cfg = config.load()
    floors = config.get(cfg, "qc.gene_floor_by_assay")
    assert floors["10x_snRNA"] < floors["10x_scRNA"], (
        "a nuclear arm held to the whole-cell floor loses its smallest nuclei "
        "preferentially, which is composition bias, not quality control")


def test_mad_outlier_handles_constant_metric():
    from ewsatlas import qc

    assert qc.mad_outlier(np.full(50, 7.0)).sum() == 0


def test_mad_outlier_flags_extremes():
    from ewsatlas import qc

    v = np.concatenate([np.random.default_rng(1).normal(10, 1, 200), [100.0]])
    assert qc.mad_outlier(v, n_mads=5)[-1]


def test_config_has_every_key_the_entrypoints_read():
    from ewsatlas import config

    cfg = config.load()
    for dotted in ("seed", "paths.raw", "paths.interim", "paths.covariates",
                   "qc.gene_floor_by_assay", "qc.mito_cap_pct_by_assay", "qc.mad_n",
                   "qc.doublet_threshold", "gene_space.min_cells_detected"):
        config.get(cfg, dotted)


def test_strict_config_get_raises_rather_than_defaulting():
    from ewsatlas import config

    with pytest.raises(KeyError):
        config.get(config.load(), "qc.nonexistent_threshold")


# --- provenance contract ----------------------------------------------------

def test_raw_directory_is_write_protected():
    from ewsatlas import config, provenance

    with pytest.raises(PermissionError):
        provenance.assert_raw_readonly([config.ROOT / "data" / "raw" / "anything.h5ad"])
    provenance.assert_raw_readonly([config.ROOT / "data" / "interim" / "ok.h5ad"])


def test_seeding_is_reproducible():
    from ewsatlas import provenance

    provenance.seed_everything(0)
    a = np.random.default_rng(0).normal(size=5)
    first = np.random.normal(size=5)
    provenance.seed_everything(0)
    assert np.allclose(first, np.random.normal(size=5))
    assert np.allclose(a, np.random.default_rng(0).normal(size=5))

# --- defect 7: CEL-seq2 tables are collision-corrected, hence fractional ------

def test_celseq_reader_rounds_to_counts(tmp_path):
    """A count likelihood is undefined on fractional input.

    Visser's deposited plate tables carry 1.8% fractional nonzero entries
    (max 29,569.7). The reader must round them and keep the original in
    ``layers['celseq_corrected']``.
    """
    import gzip
    import tarfile

    import pandas as pd
    from ewsatlas import io

    df = pd.DataFrame(
        [[1.0, 2.4, 0.0], [0.6, 3.5, 7.0], [0.0, 0.0, 1.2]],
        index=["CD99", "FLI1", "MT-CO1"], columns=["W1", "W2", "UNK"])
    plate = tmp_path / "GSM0000001_TM999.transcripts.txt.gz"
    with gzip.open(plate, "wt") as fh:
        df.to_csv(fh, sep="\t")
    tar = tmp_path / "plates.tar"
    with tarfile.open(tar, "w") as tf:
        tf.add(plate, arcname=plate.name)

    a = io.read_celseq_plates(tar, study="test")
    assert "UNK" not in [w for w in a.obs["well"]], "the UNK column must be dropped"
    d = a.X.data
    assert np.allclose(d, np.round(d)), "counts layer still fractional"
    assert "celseq_corrected" in a.layers, "original corrected values not retained"
    assert not np.allclose(a.layers["celseq_corrected"].data,
                           np.round(a.layers["celseq_corrected"].data))


# --- defect 8: a uniformly failed library survives per-cell QC ---------------

def test_library_gate_catches_a_failed_library(harmonized):
    """Per-sample MAD is computed *within* a library, so a bad one looks fine.

    Reproduces ScPCA SCPCL001111: median 44 detected genes while its neighbours
    sit above 237.
    """
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp
    from ewsatlas import qc

    good = harmonized["plain"]
    rng = np.random.default_rng(2)
    n = good.n_obs
    bad = good.copy()
    # knock all but a couple of genes out of every cell in this library
    M = sp.lil_matrix((n, good.n_vars), dtype=np.float32)
    for i in range(n):
        M[i, rng.choice(good.n_vars, size=2, replace=False)] = 1.0
    bad.X = M.tocsr()
    bad.obs["sample"] = "FAILED_LIB"
    merged = ad.concat([good, bad], join="outer", index_unique="-")
    qc.add_qc_metrics(merged)
    assert qc.failed_libraries(merged, gene_floor=5) == ["FAILED_LIB"]
    # and the per-sample MAD alone does NOT catch it
    flags = qc.flag_cells(merged, gene_floor=0, mito_cap=100.0, n_mads=5)
    in_bad = merged.obs["sample"].to_numpy() == "FAILED_LIB"
    assert flags["mad_outlier"].to_numpy()[in_bad].mean() < 0.5, (
        "if the MAD caught it, the library gate would be redundant")


# --- defect 9: a failed doublet caller was reported as zero doublets ---------

def test_doublet_failure_raises_instead_of_reporting_zero(harmonized, monkeypatch):
    """The first Waltner ingest marked all 82,496 nuclei 'failed' and the
    attrition table reported drop_doublet = 0. A failure must raise."""
    import scanpy as sc
    from ewsatlas import qc

    a = harmonized["plain"]
    qc.add_qc_metrics(a)

    def boom(*args, **kwargs):
        raise ValueError("simulated scrublet failure")

    monkeypatch.setattr(sc.pp, "scrublet", boom)
    with pytest.raises(RuntimeError, match="doublet detection failed"):
        qc.call_doublets(a, threshold=0.25, min_cells=10)

    table, errors = qc.call_doublets(a, threshold=0.25, min_cells=10, fatal=False)
    assert set(errors) == set(a.obs["sample"].unique().astype(str))
    assert (table["doublet_call"] == "failed").all()
    assert (table["doublet_call"] == "doublet").sum() == 0


def test_small_samples_are_not_assessed_not_cleaned(harmonized):
    """Under min_cells the call must be 'not_assessed', never 'singlet'."""
    from ewsatlas import qc

    a = harmonized["plain"]
    qc.add_qc_metrics(a)
    table, errors = qc.call_doublets(a, threshold=0.25, min_cells=10 ** 6)
    assert not errors
    assert (table["doublet_call"] == "not_assessed").all()


def test_expected_rate_scales_with_library_size():
    """0.8% per 1,000 cells: the cut must track library size, not a score."""
    import pandas as pd
    import scanpy as sc
    from ewsatlas import qc

    class _Fake:
        def __init__(self, n):
            self.n = n

    def fake_scrublet(sub, threshold=None, verbose=False):
        rng = np.random.default_rng(3)
        sub.obs["doublet_score"] = rng.random(sub.n_obs)
        sub.obs["predicted_doublet"] = False

    import anndata as ad
    import scipy.sparse as sp

    rng = np.random.default_rng(4)
    sizes = {"small": 2000, "big": 12000}
    parts = []
    for name, n in sizes.items():
        parts.append(ad.AnnData(
            X=sp.csr_matrix(rng.poisson(1.0, size=(n, 6)).astype(np.float32)),
            obs=pd.DataFrame({"sample": name}, index=[f"{name}-{i}" for i in range(n)])))
    A = ad.concat(parts, join="outer")
    import pytest as _pytest
    mp = _pytest.MonkeyPatch()
    mp.setattr(sc.pp, "scrublet", fake_scrublet)
    try:
        table, errors = qc.call_doublets(A, threshold=0.25, method="expected_rate",
                                         rate_per_1k=0.008)
    finally:
        mp.undo()
    assert not errors
    frac = {s: (table.loc[A.obs["sample"].to_numpy() == s, "doublet_call"] == "doublet").mean()
            for s in sizes}
    assert abs(frac["small"] - 0.016) < 0.002, frac
    assert abs(frac["big"] - 0.096) < 0.002, frac
    assert frac["big"] > frac["small"]


# --- defect 10: the "harmonized" table carried two spellings of one fact -----

def test_covariates_obey_the_controlled_vocabulary():
    """sex held both 'M' and 'male'; sort_strategy held 'cd45+' and
    'live_CD45pos_FACS'. A harmonized atlas whose own columns carry synonyms has
    harmonized nothing."""
    import pandas as pd
    from ewsatlas import config

    cfg = config.load()
    p = config.ROOT / config.get(cfg, "paths.covariates")
    if not p.is_file():
        pytest.skip("covariate table not present")
    cov = pd.read_csv(p)
    vocab = config.get(cfg, "controlled_vocabulary")
    offenders = {}
    for col, allowed in vocab.items():
        if col.startswith("_") or col not in cov.columns:
            continue
        bad = sorted(set(cov[col].astype(str)) - set(allowed))
        if bad:
            offenders[col] = bad
    assert not offenders, offenders


def test_by_assay_config_blocks_key_on_the_controlled_vocabulary():
    """config said 'CEL-seq2' while the covariate table said 'CELseq2_scRNA' —
    two spellings of one fact, which broke a figure that joined them."""
    from ewsatlas import config

    cfg = config.load()
    allowed = set(config.get(cfg, "controlled_vocabulary.assay"))
    for block in ("gene_floor_by_assay", "mito_cap_pct_by_assay"):
        keys = set(config.get(cfg, f"qc.{block}"))
        assert keys == allowed, (block, sorted(keys), sorted(allowed))
