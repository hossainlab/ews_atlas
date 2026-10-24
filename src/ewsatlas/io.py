"""Readers for the six source deposits, each a different on-disk shape.

One function per deposit format because the studies genuinely differ:

=====================  ====================================================
study                  shape
=====================  ====================================================
Aynaud GSE130024       10x MTX triplets; FOUR samples carry CellRanger v2
                       ``genes.tsv`` with ``GRCh38_``-prefixed symbols while
                       PDX-1058 ships plain v3 ``features.tsv``
Cillo GSE198896        per-sample nested tar.gz, each an uncompressed MTX dir
Goodspeed GSE277083    dense gene x cell TSV + separate cell metadata
Mao GSE250523          10x MTX triplets
Visser GSE243347       CEL-seq2 gene x well tables, ENSG__SYMBOL row names,
                       plus a per-plate annotation table
Waltner SCPCP000015    ScPCA per-library AnnData inside one zip
=====================  ====================================================

:func:`harmonize_symbols` is applied by every reader. Without it the gene-space
intersection silently matches only the unprefixed Aynaud sample and zeroes the
other four samples' cells — a failure that shows up downstream as hundreds of
uniform ~50-cell clusters from one study plus a 0.0% mitochondrial fraction.
"""

from __future__ import annotations

import gzip
import io as _io
import re
import tarfile
import zipfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.sparse as sp

#: barnyard / genome-build prefixes seen on multi-species CellRanger references
BARNYARD_RE = re.compile(r"^(GRCH38|HG19|MM10|GRCM38)[_-]")


def harmonize_symbols(symbols) -> list[str]:
    """Upper-case, drop any ``ENSG__`` prefix, strip barnyard build prefixes."""
    return [BARNYARD_RE.sub("", str(s).split("__")[-1].upper()) for s in symbols]


#: synthetic spike-in features that must never enter an expression matrix
SPIKEIN_RE = re.compile(r"^(ERCC-\d+|SIRV\d*|RGC-)", re.I)


def split_spikeins(adata):
    """Move spike-in features out of the matrix into per-cell QC covariates.

    Visser's CEL-seq2 tables ship 92 ERCC transcripts interleaved with the
    genes. Left in place they are counted as detected genes, inflate
    ``total_counts``, and dilute ``pct_counts_mt`` — in plate TM768 they are
    **60% of all counts**, which made a low-input plate look like a plate with
    few genes but many UMIs and tripped the library gate for the wrong reason.

    They are also a genuine quality metric: a high spike-in fraction means
    little endogenous RNA was captured. So they are removed from the matrix and
    retained as ``obs['pct_counts_spikein']``.
    """
    mask = pd.Index(adata.var_names).str.match(SPIKEIN_RE)
    if not mask.any():
        adata.obs["pct_counts_spikein"] = 0.0
        return adata
    total = np.asarray(adata.X.sum(axis=1)).ravel()
    spike = np.asarray(adata[:, mask].X.sum(axis=1)).ravel()
    out = adata[:, ~mask].copy()
    out.obs["pct_counts_spikein"] = np.where(total > 0, 100.0 * spike / total, 0.0)
    out.uns["n_spikein_features_removed"] = int(mask.sum())
    return out


def _anndata_from_mtx(mtx: bytes, barcodes: bytes, features: bytes, *, gzipped: bool, sample: str):
    opener = (lambda b: gzip.open(_io.BytesIO(b), "rb")) if gzipped else (lambda b: _io.BytesIO(b))
    X = sio.mmread(opener(mtx)).tocsr()
    bc = pd.read_csv(opener(barcodes), header=None, sep="\t")[0].astype(str).tolist()
    ft = pd.read_csv(opener(features), header=None, sep="\t")
    sym = harmonize_symbols(ft[1] if ft.shape[1] > 1 else ft[0])
    if X.shape[0] == len(sym):  # genes x cells -> cells x genes
        X = X.T.tocsr()
    a = ad.AnnData(
        X=X.astype(np.float32),
        obs=pd.DataFrame(index=[f"{sample}|{b}" for b in bc]),
        var=pd.DataFrame(index=pd.Index(sym, name="gene")),
    )
    a.var_names_make_unique()
    a.obs["sample"] = sample
    return a


def read_10x_tar(tar_path: Path, study: str) -> ad.AnnData:
    """Aynaud / Mao: flat MTX triplets inside one RAW.tar."""
    with tarfile.open(tar_path) as tf:
        members = tf.getnames()
        pattern = re.compile(r"GSM\d+_(.+?)_(barcodes|genes|features|matrix)\.")
        samples = sorted({m.group(1) for n in members if (m := pattern.match(n))})
        parts = []
        for s in samples:
            pick = lambda kind: [n for n in members if f"_{s}_{kind}" in n]
            f_feat = pick("genes") or pick("features")
            parts.append(_anndata_from_mtx(
                tf.extractfile(pick("matrix")[0]).read(),
                tf.extractfile(pick("barcodes")[0]).read(),
                tf.extractfile(f_feat[0]).read(),
                gzipped=True, sample=s))
    out = ad.concat(parts, join="outer")
    out.obs["study"] = study
    return out


def read_nested_tar(tar_path: Path, study: str, keep_prefixes: list[str]) -> ad.AnnData:
    """Cillo: one tar.gz per library, each containing an uncompressed MTX dir.

    ``keep_prefixes`` selects the in-scope libraries — GSE198896 also contains
    osteosarcoma and healthy-donor samples, which are out of scope for an
    EWSR1::FLI1 atlas.
    """
    parts = []
    with tarfile.open(tar_path) as tf:
        for prefix in keep_prefixes:
            member = next(n for n in tf.getnames() if n.startswith(prefix))
            blob = tf.extractfile(member).read()
            with tarfile.open(fileobj=_io.BytesIO(blob)) as inner:
                names = inner.getnames()
                g = lambda suffix: next(n for n in names if n.endswith(suffix))
                parts.append(_anndata_from_mtx(
                    inner.extractfile(g("matrix.mtx")).read(),
                    inner.extractfile(g("barcodes.tsv")).read(),
                    inner.extractfile(g("features.tsv")).read(),
                    gzipped=False, sample=prefix.split("_")[0]))
    out = ad.concat(parts, join="outer")
    out.obs["study"] = study
    return out


def read_dense_tsv(matrix_path: Path, metadata_path: Path, study: str, chunk_genes: int = 4000) -> ad.AnnData:
    """Goodspeed: dense gene x cell TSV, read in gene chunks to bound memory."""
    meta = pd.read_csv(metadata_path, sep="\t", index_col=0)
    blocks, genes, cells = [], [], None
    for ch in pd.read_csv(matrix_path, sep="\t", index_col=0, chunksize=chunk_genes):
        genes.extend(ch.index.tolist())
        blocks.append(sp.csr_matrix(ch.values.astype(np.float32)))
        cells = ch.columns.tolist()
    a = ad.AnnData(X=sp.vstack(blocks).T.tocsr(),
                   obs=pd.DataFrame(index=cells),
                   var=pd.DataFrame(index=pd.Index(harmonize_symbols(genes), name="gene")))
    a.var_names_make_unique()
    shared = a.obs_names.intersection(meta.index)
    a = a[shared].copy()
    a.obs["sample"] = meta.loc[shared, "Sample"].astype(str).values
    for col in ("scDblFinder.class", "Cell.Type", "Broad.Cell.Type", "Phase", "percent_mito"):
        if col in meta.columns:
            a.obs["src_" + col.replace(".", "_")] = meta.loc[shared, col].values
    a.obs_names = [f"{s}|{b}" for s, b in zip(a.obs["sample"], a.obs_names)]
    a.obs["study"] = study
    return a


def read_celseq_plates(tar_path: Path, study: str, annotation_path: Path | None = None) -> ad.AnnData:
    """Visser: one gene x well table per plate, resolved against the deposit's
    own annotation table rather than the publication.

    ``GSE243347_data_annotation.txt.gz`` is richer than the paper and is treated
    as authoritative for three things the matrices cannot carry:

    * **declared empty wells** - 16 per plate (column 1), 448 in total. Left in
      place they are counted as cells, and because the library gate uses a
      per-plate median they drag it down; they are dropped by name.
    * **pooled plates** - TM736 holds two patients, ES-039 in columns 2-3 and
      ES-030 in columns 4-18. The annotation lists the wells, so the plate is
      emitted as two samples (``TM736-1``, ``TM736-2``) instead of one
      placeholder. Wells on a pooled plate that are neither declared empty nor
      listed for a patient are unassigned and dropped.
    * **plate TM768** - patient ES-048, lymph-node metastasis, post-treatment:
      documented here and omitted from earlier builds entirely.

    Counts are not integers: 1.8% of nonzero entries are fractional because
    CEL-seq2 UMI counts are collision-corrected (the binomial correction that
    maps observed barcodes back to original molecules yields real numbers, and
    it grows with depth - the largest value here is 29,569.7). A
    negative-binomial likelihood is defined on counts, so the matrix is rounded
    at ingest: the median nonzero value is 1.0, so rounding moves 1.8% of
    entries by less than one molecule each, while passing fractional values to a
    count model is undefined. The original stays in
    ``layers['celseq_corrected']``.

    ERCC spike-ins are removed from the matrix and retained as
    ``obs['pct_counts_spikein']`` - see :func:`split_spikeins`.
    """
    ann = None
    if annotation_path is not None and Path(annotation_path).is_file():
        ann = pd.read_csv(annotation_path, sep="\t")
        ann["_plate"] = ann["plate_id"].astype(str).str.split("-").str[0]

    parts = []
    with tarfile.open(tar_path) as tf:
        for name in sorted(tf.getnames()):
            plate = re.match(r"GSM\d+_(.+?)\.transcripts", name).group(1)
            df = pd.read_csv(_io.BytesIO(tf.extractfile(name).read()), sep="\t",
                             index_col=0, compression="gzip")
            df = df.drop(columns=[c for c in ("UNK",) if c in df.columns])

            rows = ann[ann["_plate"] == plate] if ann is not None else None
            if rows is None or rows.empty:
                groups = [(plate, list(df.columns))]
            else:
                empty = {w.strip() for s in rows["empty_wells"].dropna()
                         for w in str(s).split(",") if w.strip()}
                groups = []
                for _, r in rows.iterrows():
                    if pd.isna(r.get("selected_wells")):
                        wells = [w for w in df.columns if w not in empty]
                    else:
                        sel = {w.strip() for w in str(r["selected_wells"]).split(",") if w.strip()}
                        wells = [w for w in df.columns if w in sel and w not in empty]
                    if wells:
                        groups.append((str(r["plate_id"]), wells))

            for sample_id, wells in groups:
                sub = df[wells]
                a = ad.AnnData(
                    X=sp.csr_matrix(sub.T.values.astype(np.float32)),
                    obs=pd.DataFrame(index=[f"{sample_id}|{w}" for w in wells]),
                    var=pd.DataFrame(index=pd.Index(harmonize_symbols(sub.index), name="gene")),
                )
                a.var_names_make_unique()
                a.obs["sample"] = sample_id
                a.obs["plate"] = plate
                a.obs["well"] = list(wells)
                a.layers["celseq_corrected"] = a.X.copy()
                a.X = a.X.copy()
                a.X.data = np.rint(a.X.data).astype(np.float32)
                a.X.eliminate_zeros()
                a = split_spikeins(a)
                parts.append(a)

    n_spike = max(int(a.uns.get("n_spikein_features_removed", 0)) for a in parts)
    out = ad.concat(parts, join="outer")
    out.obs["study"] = study
    out.uns["n_spikein_features_removed"] = n_spike
    return out


def read_scpca_zip(zip_path: Path, study: str, which: str = "filtered") -> ad.AnnData:
    """Waltner: ScPCA per-library AnnData inside one archive.

    ``which='filtered'`` takes ScPCA's emptyDrops-filtered matrices as the
    starting point and re-applies *our* QC on top, so the nuclei pass the same
    policy as the five GEO studies. ScPCA's own annotations are kept only as
    ``src_*`` provenance columns.
    """
    parts = []
    with zipfile.ZipFile(zip_path) as zf:
        targets = [n for n in zf.namelist() if n.endswith(f"_{which}_rna.h5ad")]
        for name in sorted(targets):
            with zf.open(name) as fh:
                a = ad.read_h5ad(_io.BytesIO(fh.read()))
            library = Path(name).name.split("_")[0]
            sample = Path(name).parent.name
            a.var_names = harmonize_symbols(
                a.var["gene_symbol"] if "gene_symbol" in a.var else a.var_names)
            a.var_names_make_unique()
            # ScPCA's own QC and annotation are retained as provenance only, never
            # as the atlas's labels: scDblFinder/miQC let us cross-check our own
            # doublet and compromised-cell calls, and is_xenograft cross-checks the
            # curated covariate table against the deposit.
            keep = {c: f"src_{c}" for c in a.obs.columns
                    if any(k in c for k in ("celltype", "cluster", "miQC", "scpca",
                                            "submitter", "scDblFinder", "prob_compromised",
                                            "is_xenograft", "participant_id",
                                            "subsets_mito_percent"))}
            a.obs = a.obs[list(keep)].rename(columns=keep)
            a.obs["sample"] = sample
            a.obs["library"] = library
            a.obs_names = [f"{sample}|{b}" for b in a.obs_names]
            parts.append(a)
    out = ad.concat(parts, join="outer")
    out.obs["study"] = study
    return out

def write_h5ad_atomic(adata, path, compression: str = "gzip") -> None:
    """Write an AnnData object without ever truncating the existing file.

    ``AnnData.write_h5ad`` hands the target to h5py with truncate-on-create, so
    a write that fails part way - a stale read handle elsewhere in the process,
    a full disk, an interrupt - leaves a zero-byte file where the input used to
    be. That happened to ``EWSAtlas_v5.h5ad`` during the stage 05 rewrite: the
    copy-number computation succeeded, the write failed on a lock held by the
    calling kernel, and the released object was destroyed.

    Writing to a sibling temporary file and then atomically replacing the target
    makes the failure mode "the old file is still there" instead.
    """
    import os
    from pathlib import Path

    path = Path(path)
    tmp = path.with_name(path.name + ".tmp")
    if tmp.exists():
        tmp.unlink()
    adata.write_h5ad(tmp, compression=compression)
    os.replace(tmp, path)
