"""Uniform quality control — one *criterion* for every study, not one number.

Scientific rationale for the assay-aware gene floor
---------------------------------------------------
A detected-gene floor is a proxy for "is this barcode a real cell with enough
signal to assign a state", and the mapping from that proxy to library quality
is assay-dependent. A nucleus holds roughly 10-20% of a cell's RNA and is
depleted of cytoplasmic and mitochondrial transcripts, so a nucleus with 150
detected genes is a normal nucleus, not a degraded cell. Applying one integer to
both assays is uniform only in arbitrary units: the *effective* stringency is
far harsher on the nuclear arm.

The failure mode is specific. A fixed floor does not remove nuclei at random —
it removes the smallest and least transcriptionally active ones, and because
transcriptional output tracks cell type and cycling state, the survivors skew
toward large, active, proliferating malignant cells. That converts a technical
threshold into a composition bias, which is precisely the artefact an atlas
exists to avoid. In this dataset a uniform 200-gene floor would discard ~45% of
the Waltner nuclei (median 231 genes, 296 UMIs).

So: the *criterion* is uniform (per-sample MAD outlier removal, applied
identically everywhere), and the absolute floor acts only as a backstop against
empty droplets, set per assay.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

from .compat import sequential_joblib


def mad_outlier(values, n_mads: int = 5) -> np.ndarray:
    """Median-absolute-deviation outlier mask (scuttle convention, 1.4826 scale).

    Adaptive rather than fixed so it transfers across libraries of different
    depth; returns all-False when MAD is zero (a constant metric cannot have
    outliers and must not nuke the sample).
    """
    v = np.asarray(values, dtype=float)
    med = np.median(v)
    mad = np.median(np.abs(v - med))
    if mad == 0:
        return np.zeros(v.shape[0], dtype=bool)
    return np.abs(v - med) > n_mads * 1.4826 * mad


def add_qc_metrics(adata) -> None:
    """Attach scanpy QC metrics plus mitochondrial / ribosomal flags in place."""
    adata.var["mt"] = pd.Index(adata.var_names).str.startswith(("MT-", "MT."))
    adata.var["ribo"] = pd.Index(adata.var_names).str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True,
                               percent_top=None, log1p=False)


def flag_cells(adata, *, gene_floor: int, mito_cap: float, n_mads: int) -> pd.DataFrame:
    """Return a per-cell flag table; callers decide what to drop.

    Returning flags rather than a filtered object keeps the attrition auditable:
    every stage can report how many cells each criterion removed, which is what
    the QC table in the paper is built from.
    """
    f = pd.DataFrame(index=adata.obs_names)
    f["below_gene_floor"] = adata.obs["n_genes_by_counts"].to_numpy() < gene_floor
    f["above_mito_cap"] = adata.obs["pct_counts_mt"].to_numpy() > mito_cap
    mad = np.zeros(adata.n_obs, dtype=bool)
    for _, idx in adata.obs.groupby("sample", observed=True).indices.items():
        mad[idx] = (
            mad_outlier(np.log1p(adata.obs["total_counts"].to_numpy()[idx]), n_mads)
            | mad_outlier(np.log1p(adata.obs["n_genes_by_counts"].to_numpy()[idx]), n_mads)
            | mad_outlier(adata.obs["pct_counts_mt"].to_numpy()[idx], n_mads)
        )
    f["mad_outlier"] = mad
    return f


def failed_libraries(adata, *, gene_floor: int) -> list[str]:
    """Libraries whose MEDIAN detected-gene count is below the assay floor.

    A per-cell floor plus a per-sample MAD cannot catch a uniformly failed
    library: the MAD is computed *within* that library, so a bad library looks
    internally consistent and survives. The median is the cheapest statistic
    that separates "nuclei are small" from "this library did not work".

    Measured case in this atlas: ScPCA library SCPCL001111 has median 44
    detected genes and 76 UMIs, with 98.5% of nuclei under 200 genes; ScPCA's
    own pipeline retained 89 of its 6,139 cells. Every other library in the
    deposit has a median of at least 237.
    """
    med = adata.obs.groupby("sample", observed=True)["n_genes_by_counts"].median()
    return sorted(med.index[med < gene_floor].astype(str))


def call_doublets(adata, *, threshold: float, min_cells: int = 100,
                  fatal: bool = True, method: str = "expected_rate",
                  rate_per_1k: float = 0.008) -> tuple[pd.DataFrame, dict[str, str]]:
    """Per-sample scrublet with an EXPLICIT threshold.

    scanpy's automatic threshold needs scikit-image, which is refused by host
    policy in this install, so the cut is supplied rather than inferred. Samples
    below ``min_cells`` are marked ``not_assessed`` instead of being forced
    through a simulation that cannot support it.

    Returns ``(table, errors)``. **A failure is never reported as zero
    doublets**: with ``fatal=True`` (the default) any per-sample exception
    raises, because the first run of this stage silently marked all 82,496
    Waltner nuclei ``failed`` and the attrition table dutifully reported
    ``drop_doublet = 0`` — a total failure that reads as a clean result.

    The matrix is forced to canonical CSR (deduplicated, sorted indices) before
    scrublet. ``ad.concat(join="outer")`` can leave non-canonical indices that
    h5ad serialization silently normalizes, which is why the same data failed
    in-memory and succeeded after a write/read round-trip.

    Two ways to turn scrublet scores into calls:

    ``method="expected_rate"`` (default)
        Flag the top-scoring cells per sample at the rate the chemistry
        predicts: 10x multiplet rate is ~0.8% per 1,000 cells recovered
        (``rate_per_1k``), so a 7,000-nucleus library expects ~5.6% and a
        14,000-nucleus library ~11%. The *criterion* is uniform — an instrument
        prior applied identically — while the resulting cut adapts to library
        size, which is what actually drives multiplet formation.

    ``method="fixed_threshold"``
        Apply ``threshold`` to the score directly. Retained because it is what
        earlier builds used, and because it is the right choice if a future
        environment can compute scrublet's own bimodal threshold.

    Why this matters here specifically: a malignant-plus-stromal doublet reads
    as a malignant cell with a mesenchymal programme, which is exactly the
    fusion-low state this atlas sets out to study. Under-calling doublets would
    manufacture the finding. On the Waltner arm the fixed 0.25 cut removed 2.9%
    of nuclei while ScPCA's independent scDblFinder flagged 9.9% of the same
    retained nuclei (score Spearman rho 0.41, median score 0.14 in their
    doublet class vs 0.09 in their singlets) — concordant in direction, far
    apart in stringency.
    """
    call = np.full(adata.n_obs, "not_assessed", dtype=object)
    score = np.full(adata.n_obs, np.nan)
    errors: dict[str, str] = {}
    for name, idx in adata.obs.groupby("sample", observed=True).indices.items():
        if len(idx) < min_cells:
            continue
        sub = adata[idx].copy()
        sub.X = sp.csr_matrix(sub.X)
        sub.X.sum_duplicates()
        sub.X.sort_indices()
        try:
            with sequential_joblib():  # see ewsatlas.compat for why
                sc.pp.scrublet(sub, threshold=threshold, verbose=False)
        except Exception as exc:  # noqa: BLE001 - recorded, then re-raised
            errors[str(name)] = f"{type(exc).__name__}: {exc}"
            call[idx] = "failed"
            continue
        s = sub.obs["doublet_score"].to_numpy()
        score[idx] = s
        if method == "expected_rate":
            k = int(round(min(0.5, rate_per_1k * len(idx) / 1000.0) * len(idx)))
            is_dbl = np.zeros(len(idx), dtype=bool)
            if k > 0:
                is_dbl[np.argsort(-s)[:k]] = True
        elif method == "fixed_threshold":
            is_dbl = s >= threshold
        else:
            raise ValueError(f"unknown doublet method: {method!r}")
        call[idx] = np.where(is_dbl, "doublet", "singlet")
    if errors and fatal:
        raise RuntimeError(
            f"doublet detection failed for {len(errors)} sample(s): {errors}. "
            "Refusing to continue: a failed caller would be recorded as zero doublets.")
    return pd.DataFrame({"doublet_score": score, "doublet_call": call}, index=adata.obs_names), errors


def apply_qc(adata, *, gene_floor: int, mito_cap: float, n_mads: int,
             doublet_threshold: float, library_gate: bool = True,
             doublet_method: str = "expected_rate", doublet_rate_per_1k: float = 0.008):
    """Full QC for one study. Returns ``(filtered, attrition_row)``."""
    add_qc_metrics(adata)
    n0 = adata.n_obs
    failed = failed_libraries(adata, gene_floor=gene_floor) if library_gate else []
    n_failed_cells = int(adata.obs["sample"].astype(str).isin(failed).sum())
    if failed:
        adata = adata[~adata.obs["sample"].astype(str).isin(failed)].copy()
    flags = flag_cells(adata, gene_floor=gene_floor, mito_cap=mito_cap, n_mads=n_mads)
    kept = adata[~flags.any(axis=1).to_numpy()].copy()
    dbl, dbl_errors = call_doublets(kept, threshold=doublet_threshold,
                                    method=doublet_method, rate_per_1k=doublet_rate_per_1k)
    kept.obs["doublet_score"] = dbl["doublet_score"].to_numpy()
    kept.obs["doublet_call"] = dbl["doublet_call"].to_numpy()
    out = kept[kept.obs["doublet_call"].to_numpy() != "doublet"].copy()
    row = {
        "study": str(adata.obs["study"].iloc[0]),
        "raw": n0,
        "gene_floor_applied": gene_floor,
        "mito_cap_applied": mito_cap,
        "failed_libraries": ";".join(failed) if failed else "",
        "drop_failed_library": n_failed_cells,
        "drop_below_gene_floor": int(flags["below_gene_floor"].sum()),
        "drop_above_mito_cap": int(flags["above_mito_cap"].sum()),
        "drop_mad_outlier": int((flags["mad_outlier"] & ~flags["below_gene_floor"]
                                 & ~flags["above_mito_cap"]).sum()),
        "doublet_method": doublet_method,
        "drop_doublet": int(kept.n_obs - out.n_obs),
        "doublet_not_assessed": int((kept.obs["doublet_call"] == "not_assessed").sum()),
        "doublet_errors": ";".join(dbl_errors) if dbl_errors else "",
        "kept": out.n_obs,
        "pct_kept": round(100 * out.n_obs / n0, 1),
    }
    return out, row
