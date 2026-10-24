"""Atlas as a deconvolution reference for 343 bulk Ewing tumours, with survival.

The atlas carries no outcome data, so the only way it can speak to prognosis is
as a cell-type reference projected into cohorts that do. 408 GEO expression records de-duplicate to 343 primary
records, one per patient, which are the samples deconvolved here; 129 unique patients (GSE63157
n=85, GSE17674 n=44) supply survival after the two accession-level duplications
documented in docs/07 are resolved.

Two choices that a reader should be able to check:

  * NNLS runs in LINEAR space. The first pass of this analysis (v4) was run on
    log2 microarray values and returned ~0% malignant content for 369 Ewing
    tumours. A mixture is additive on the linear scale and not after a log, so
    every bulk column is exponentiated before fitting.
  * Urothelial and Unassigned cells are excluded from the signature. Urothelium
    is Mao's bladder control and has no counterpart in bone or soft-tissue
    tumours; leaving it in gives it tumour signal to absorb.
"""
import sys
import pathlib
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
from scipy.optimize import nnls
import statsmodels.api as sm

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ewsatlas import config  # noqa: E402

DROP_TYPES = ["Urothelial", "Unassigned"]
N_MARKERS = 100
MIN_MEAN = 0.10
STROMAL = ["Fibroblast_CAF", "Endothelial"]
# Two reference profiles that differ only by cell-cycle genes are collinear, and
# NNLS resolves collinearity arbitrarily: with Ewing_sarcoma and
# Ewing_proliferating as separate columns the fit put 0.358 on the proliferating
# profile and 0.041 on the class that is eight times larger in the atlas. Bulk
# microarray cannot separate them, so they enter as one malignant profile.
MERGE = {"Ewing_sarcoma": "Malignant", "Ewing_proliferating": "Malignant"}
# The signature is built from WHOLE-CELL libraries only. 68,698 of the 97,586
# frozen malignant cells are single nuclei, whose transcriptome is depleted of
# cytoplasmic message and enriched for intronic and long-gene signal; a malignant
# profile dominated by nuclei does not resemble a bulk tissue extract, and the
# nucleus-based signature returned a mean malignant fraction of 0.33-0.53 against
# 0.64-0.78 for the whole-cell one.
WHOLE_CELL_ASSAYS = ["10x_scRNA", "CELseq2_scRNA"]


def signature(A: ad.AnnData, genes) -> pd.DataFrame:
    """Mean counts-per-10k per cell type, on the genes shared with the bulk."""
    # expm1 is applied per cell, before averaging: the mean of log values is not
    # the log of the mean. It preserves sparsity (expm1(0) = 0), so the matrix is
    # never densified - doing so needs 15.7 GiB. Per-type means are taken as a
    # single indicator-matrix product (types x cells @ cells x genes), which
    # allocates a 12 x 16,622 result instead of slicing out one group at a time.
    X = A.layers["lognorm"]
    if sp.issparse(X):
        X = X.tocsr(copy=True).astype(np.float32)
        X.data = np.expm1(X.data)
    else:
        X = np.expm1(np.asarray(X, dtype=np.float32))
    ct = pd.Series(A.obs["cell_type"].astype(str).to_numpy()).replace(MERGE).to_numpy()
    assay = A.obs["assay"].astype(str).to_numpy()
    ct = np.where(np.isin(assay, WHOLE_CELL_ASSAYS), ct, "__excluded__")
    types = [t for t in pd.unique(ct) if t not in DROP_TYPES and t != "__excluded__"]
    code = {t: i for i, t in enumerate(types)}
    keep = np.array([t in code for t in ct])
    rows = np.array([code[t] for t in ct[keep]])
    G = sp.csr_matrix((np.ones(keep.sum(), dtype=np.float32),
                       (rows, np.flatnonzero(keep))), shape=(len(types), A.n_obs))
    sums = np.asarray((G @ X).todense() if sp.issparse(X) else G @ X, dtype=np.float64)
    n = np.asarray(G.sum(1)).ravel()
    out = pd.DataFrame((sums / n[:, None]).T, index=list(A.var_names), columns=types)
    return out.loc[list(genes)]


def pick_markers(S: pd.DataFrame) -> list:
    """Specificity-selected markers: share of a gene's total signal in one type."""
    frac = S.div(S.sum(1).replace(0, np.nan), axis=0)
    out = []
    for t in S.columns:
        ok = S[t] >= MIN_MEAN
        out += list(frac.loc[ok, t].sort_values(ascending=False).head(N_MARKERS).index)
    return sorted(set(out))


def cox(df, cov, tcol, ecol, strata=None):
    d = df[[cov, tcol, ecol] + ([strata] if strata else [])].dropna()
    x = ((d[cov] - d[cov].mean()) / d[cov].std()).to_numpy()[:, None]
    kw = {"strata": d[strata].to_numpy()} if strata else {}
    m = sm.PHReg(d[tcol].to_numpy(), x, status=d[ecol].to_numpy(), **kw).fit()
    hr = float(np.exp(m.params[0]))
    lo, hi = np.exp(m.conf_int()[0])
    return dict(covariate=cov, endpoint=tcol.split("_")[0].upper(), n=len(d),
                events=int(d[ecol].sum()), hr=hr, ci_low=float(lo), ci_high=float(hi),
                p=float(m.pvalues[0]))


def main() -> int:
    cfg = config.load()
    PROJ = config.ROOT
    bulk_dir = PROJ / "data/external/bulk_cohorts"
    B = pd.read_parquet(bulk_dir / "ews_bulk_expression.parquet")
    C = pd.read_csv(bulk_dir / "ews_bulk_clinical.csv")
    A = ad.read_h5ad(PROJ / config.get(cfg, "paths.atlas_out"))

    shared = [g for g in A.var_names if g in B.index]
    S = signature(A, shared)
    print(f"{len(shared):,} shared genes, {S.shape[1]} cell types", flush=True)

    # Missing values in the harmonized bulk matrix are whole genes absent from a
    # platform (16.5% of cells, but by gene x cohort block, not scattered), so
    # each cohort is deconvolved against markers it actually measures. Fractions
    # are compared within cohort, never across.
    meta_all = C[C.included_in_matrix & C.is_primary_record].set_index("sample_id")
    cohort_of = meta_all.reindex(B.columns)["cohort"]
    frac, used = {}, []
    for ch, cols in cohort_of.groupby(cohort_of):
        sub = B[list(cols.index)]
        ok = [g for g in shared if sub.loc[g].notna().all()]
        markers = pick_markers(S.loc[ok])
        Sm = S.loc[markers]
        used.append(dict(cohort=ch, n_samples=len(cols.index), n_genes_measured=len(ok),
                         n_markers=len(markers)))
        for s in cols.index:
            y = np.power(2.0, sub[s].reindex(markers).to_numpy(dtype=float))
            w, _ = nnls(Sm.to_numpy(), y)
            frac[s] = pd.Series(w / w.sum() if w.sum() > 0 else w, index=Sm.columns)
    print(pd.DataFrame(used).to_string(index=False), flush=True)
    F = pd.DataFrame(frac).T
    F["malignant"] = F["Malignant"]
    F["stromal"] = F[[c for c in STROMAL if c in F]].sum(1)
    F["immune"] = F[[c for c in F.columns if c in
                     ["T_cell", "NK", "B_cell", "Plasma", "Monocyte_Macrophage", "DC",
                      "pDC", "Mast", "Platelet"]]].sum(1)
    F.index.name = "sample_id"
    F.to_csv(PROJ / "results/tables/f6_bulk_fractions.csv")

    meta = C[C.included_in_matrix & C.is_primary_record].set_index("sample_id")
    D = F.join(meta[["cohort", "os_time_days", "os_event", "efs_time_days", "efs_event"]])
    for c in ["os_time_days", "os_event", "efs_time_days", "efs_event"]:
        D[c] = pd.to_numeric(D[c], errors="coerce")
    print("mean malignant fraction by cohort:")
    print(D.groupby("cohort")["malignant"].agg(["mean", "size"]).round(3).to_string(), flush=True)

    surv = D[D.os_time_days.notna()].copy()
    res = []
    covs = ["malignant", "stromal", "immune", "Fibroblast_CAF"]
    for cov in covs:
        for tcol, ecol in [("os_time_days", "os_event"), ("efs_time_days", "efs_event")]:
            r = cox(surv, cov, tcol, ecol, strata="cohort")
            r["cohort"] = "both (stratified)"
            res.append(r)
            for ch, d in surv.groupby("cohort"):
                r = cox(d, cov, tcol, ecol)
                r["cohort"] = ch
                res.append(r)
    R = pd.DataFrame(res)
    R.to_csv(PROJ / "results/tables/f6_bulk_survival.csv", index=False)
    print(R[R.cohort == "both (stratified)"][
        ["covariate", "endpoint", "n", "events", "hr", "ci_low", "ci_high", "p"]]
        .round(4).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
