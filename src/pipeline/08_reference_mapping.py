"""Leave-one-study-out reference mapping: can EWSAtlas label data it never saw?

For each of the six studies in turn, a scVI reference is trained on the other
five, the held-out study is mapped onto it with scArches surgery, and cell-type
labels are transferred from reference to query by k-nearest neighbours in the
shared latent space.

Three declared departures from the atlas's own integration configuration, all in
the direction of making the experiment affordable on CPU (no GPU is available in
this sandbox) rather than of flattering the result:

  * 2,000 highly variable genes rather than 5,000, 100 epochs rather than 500,
    and the reference capped at ~50,000 cells per fold by trimming large
    (study, cell type) groups while keeping small ones whole. Six reference
    trainings at the full configuration need about six and a half hours on this
    CPU; the capped configuration needs well under one.
  * Highly variable genes are selected on the REFERENCE studies only. Selecting
    them on the full atlas would leak the query into the feature space and is the
    most common way this experiment is accidentally made easy.
  * batch_key='sample' with no extra categorical covariates. In the atlas the
    model also encodes study and assay, but those are nested within sample, and
    surgery onto unseen covariate levels adds a failure mode that is not the
    subject of this test.

What this measures, and what it does not: the reference labels are the atlas's
own annotation, which was produced with the query study present. The experiment
therefore measures whether the annotation is reproducible when a study is
removed from the reference - projection stability - and NOT generalization to an
external cohort. No external Ewing sarcoma single-cell cohort exists in any
public archive searched, so the stronger experiment cannot be run here.
"""
import json
import sys
import pathlib
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import scvi
from sklearn.neighbors import KNeighborsClassifier

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ewsatlas import config  # noqa: E402

N_HVG, MAX_EPOCHS, SURGERY_EPOCHS, K = 2000, 100, 50, 15
MIN_REF_CELLS_PER_TYPE = 20
# The reference is capped at roughly this many cells per fold by trimming large
# (study, cell type) groups while keeping every small one intact, so rare types
# are not lost. At full size a fold needs ~14 s/epoch and six folds do not finish
# in a session; the cap brings a fold to a few minutes without changing which
# cell types the reference can offer.
N_REF_MAX = 50000


def cap_reference(obs: pd.DataFrame, target: int, seed: int = 0) -> np.ndarray:
    """Indices of a per-group-capped subsample; small groups are kept whole."""
    groups = list(obs.groupby([obs["study"].astype(str), obs["cell_type"].astype(str)]).groups.items())
    sizes = np.array([len(ix) for _, ix in groups])
    if sizes.sum() <= target:
        return np.arange(len(obs))
    lo, hi = 1, int(sizes.max())
    while lo < hi:  # smallest per-group cap whose total still fits
        mid = (lo + hi) // 2
        if np.minimum(sizes, mid).sum() <= target:
            lo = mid + 1
        else:
            hi = mid
    cap = max(1, lo - 1)
    rng = np.random.default_rng(seed)
    pos = pd.Series(np.arange(len(obs)), index=obs.index)
    keep = []
    for _, ix in groups:
        q = pos.loc[ix].to_numpy()
        keep.append(q if len(q) <= cap else rng.choice(q, cap, replace=False))
    return np.sort(np.concatenate(keep))


def main() -> int:
    cfg = config.load()
    PROJ = config.ROOT
    scvi.settings.seed = 0
    A = ad.read_h5ad(PROJ / config.get(cfg, "paths.atlas_out"))
    A.X = A.layers["counts"].copy()
    studies = list(pd.unique(A.obs["study"].astype(str)))

    rows, per_type = [], []
    for held in studies:
        is_q = (A.obs["study"].astype(str) == held).to_numpy()
        ref = A[~is_q]
        ref = ref[cap_reference(ref.obs, N_REF_MAX)].copy()
        qry = A[is_q].copy()

        sc.pp.highly_variable_genes(ref, n_top_genes=N_HVG, flavor="seurat_v3",
                                    batch_key="study", subset=True)
        qry = qry[:, ref.var_names].copy()

        scvi.model.SCVI.setup_anndata(ref, batch_key="sample")
        model = scvi.model.SCVI(ref, n_latent=30, n_layers=2, gene_likelihood="nb",
                                encode_covariates=True, use_layer_norm="both",
                                use_batch_norm="none")
        model.train(max_epochs=MAX_EPOCHS, early_stopping=True,
                    early_stopping_patience=12, batch_size=512)

        scvi.model.SCVI.prepare_query_anndata(qry, model)
        qmodel = scvi.model.SCVI.load_query_data(qry, model)
        qmodel.train(max_epochs=SURGERY_EPOCHS, plan_kwargs={"weight_decay": 0.0},
                     batch_size=512)

        Zr = model.get_latent_representation()
        Zq = qmodel.get_latent_representation(qry)
        yr = ref.obs["cell_type"].astype(str).to_numpy()
        yq = qry.obs["cell_type"].astype(str).to_numpy()

        keep_t = pd.Series(yr).value_counts()
        keep_t = set(keep_t[keep_t >= MIN_REF_CELLS_PER_TYPE].index)
        m = np.isin(yr, list(keep_t))
        knn = KNeighborsClassifier(n_neighbors=K, weights="distance")
        knn.fit(Zr[m], yr[m])
        pred = knn.predict(Zq)
        prob = knn.predict_proba(Zq).max(1)

        # a query cell whose true label is absent from the reference cannot be
        # got right; report those separately rather than scoring them as errors
        answerable = np.isin(yq, list(keep_t))
        acc = float((pred[answerable] == yq[answerable]).mean())
        conf = prob >= 0.6
        rows.append(dict(held_out=held, n_query=int(is_q.sum()), n_ref=int((~is_q).sum()),
                         n_ref_used=int(ref.n_obs), n_hvg=int(ref.n_vars), accuracy=acc,
                         n_answerable=int(answerable.sum()),
                         types_absent_from_reference=";".join(
                             sorted(set(yq[~answerable]))) or "none",
                         accuracy_confident=float(
                             (pred[answerable & conf] == yq[answerable & conf]).mean())
                         if (answerable & conf).sum() else np.nan,
                         pct_confident=float(100 * conf[answerable].mean())))
        for ct in sorted(set(yq[answerable])):
            s = answerable & (yq == ct)
            per_type.append(dict(held_out=held, cell_type=ct, n=int(s.sum()),
                                 recall=float((pred[s] == ct).mean()),
                                 top_confusion=pd.Series(pred[s][pred[s] != ct]).value_counts()
                                 .index[0] if (pred[s] != ct).any() else ""))
        print(f"{held}: {is_q.sum():,} query cells, accuracy {acc:.3f}", flush=True)
        del ref, qry, model, qmodel

    out = PROJ / "results/tables"
    pd.DataFrame(rows).to_csv(out / "f6_loso_mapping.csv", index=False)
    PT = pd.DataFrame(per_type)
    PT.to_csv(out / "f6_loso_by_celltype.csv", index=False)
    # Urothelial cells exist in only one study (Mao 2,433; 20 elsewhere in total)
    # and Unassigned is not a biological type, so neither can be transferred when
    # its own study is the query. Accuracy is reported both ways.
    ex = ["Urothelial", "Unassigned"]
    dec = []
    for h, d in PT.groupby("held_out"):
        k = d[~d.cell_type.isin(ex)]
        dec.append(dict(held_out=h,
                        accuracy_all_types=float((d.n * d.recall).sum() / d.n.sum()),
                        accuracy_transferable=float((k.n * k.recall).sum() / k.n.sum()),
                        n_cells_excluded=int(d[d.cell_type.isin(ex)].n.sum())))
    pd.DataFrame(dec).to_csv(out / "f6_loso_accuracy_decomposition.csv", index=False)
    print(pd.DataFrame(rows)[["held_out", "n_query", "accuracy", "pct_confident"]]
          .round(3).to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
