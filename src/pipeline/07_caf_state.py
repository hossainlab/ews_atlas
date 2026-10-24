"""F4 analysis: is the CNA-rescued fibroblast-like population genuinely aneuploid?

Non-circular design. The frozen malignant set was built by correlating each cell's
genome-wide CNV profile to a per-sample consensus. Re-using that statistic to
validate the rescued cells would be circular. Instead:

  * windows are split into two disjoint halves by alternating chromosome;
  * candidate cells are selected on half A exactly as the pipeline selects them
    (correlation above the 95th percentile of that sample's reference cells);
  * every reported statistic is computed on half B, which played no part in
    selection;
  * the half-B statistic is own-patient specificity -- correlation to the cell's
    own sample consensus minus its mean correlation to the consensus of other
    samples in the same study. Ambient RNA, dissociation stress and technical
    noise raise correlation to any consensus; only a shared clonal karyotype
    raises it to one patient in particular.

Consensus profiles are seeded from marker-defined malignant cells only, so the
fibroblast-like candidates never contribute to the target they are scored against.
"""
import sys
import pathlib
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import infercnvpy as cnv
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ewsatlas import compat

PROJ = pathlib.Path(__file__).resolve().parents[2]
REFERENCE = ["T_cell", "B_cell", "Plasma", "Monocyte_Macrophage", "DC", "pDC",
             "Platelet", "Endothelial", "Urothelial"]
WINDOW, STEP, MIN_REF, MIN_MAL = 250, 10, 50, 20
FIBRO_LIKE = ["Fibroblast_CAF", "Unassigned"]


def corr_rows(mat, vec):
    """Pearson correlation of every row of mat against vec."""
    m = mat - mat.mean(1, keepdims=True)
    v = vec - vec.mean()
    d = np.sqrt((m ** 2).sum(1)) * np.sqrt((v ** 2).sum())
    return np.where(d > 0, (m @ v) / np.where(d > 0, d, 1.0), 0.0)


def main(downsample: int = 0) -> int:
    A = ad.read_h5ad(PROJ / "data/interim/caf_scored.h5ad")
    if downsample:
        # Own-patient specificity rises with sequencing depth, and the candidate
        # cells are 5-7x deeper than reference cells. Restricting to cells above
        # the target depth and downsampling every one of them to exactly that
        # depth removes the confound at the cost of cell number.
        import scanpy as sc
        A = A[np.asarray(A.layers["counts"].sum(1)).ravel() >= downsample].copy()
        A.X = A.layers["counts"].copy()
        sc.pp.downsample_counts(A, counts_per_cell=downsample)
        A.layers["counts"] = A.X.copy()
        sc.pp.normalize_total(A, target_sum=1e4)
        sc.pp.log1p(A)
        A.layers["lognorm"] = A.X.copy()
        print(f"depth-matched to {downsample} counts/cell: {A.n_obs:,} cells retained", flush=True)
    pos = pd.read_csv(PROJ / "data/external/references/gene_positions.csv", index_col=0)
    for c in ["chromosome", "start", "end"]:
        A.var[c] = pos.reindex(A.var_names)[c].values
    A = A[:, A.var["chromosome"].notna()].copy()
    A.X = A.layers["lognorm"].copy()

    rows = []
    for study in A.obs["study"].unique():
        sub = A[A.obs["study"] == study].copy()
        # infercnvpy indexes its gene list positionally; pandas 3 hands it an
        # arrow-backed string index, which raises ArrowInvalid. Same cast as 05_cnv.
        sub.var.index = pd.Index(list(sub.var_names), dtype=object, name="gene")
        ct = sub.obs["cell_type"].astype(str)
        ref = [c for c in REFERENCE if (ct == c).sum() >= 20]
        if not ref:
            continue
        with compat.sequential_process_map():
            cnv.tl.infercnv(sub, reference_key="cell_type", reference_cat=ref,
                            window_size=WINDOW, step=STEP, n_jobs=1)
        X = sub.obsm["X_cnv"]
        X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
        chr_pos = sub.uns["cnv"]["chr_pos"]
        order = sorted(chr_pos.items(), key=lambda kv: kv[1])
        wchr = np.empty(X.shape[1], dtype=object)
        for i, (c, s0) in enumerate(order):
            e = order[i + 1][1] if i + 1 < len(order) else X.shape[1]
            wchr[s0:e] = c
        chroms = [c for c, _ in order]
        halves = {"A": np.isin(wchr, chroms[0::2]), "B": np.isin(wchr, chroms[1::2])}

        ev = sub.obs["malignant_evidence"].astype(str)
        is_marker_mal = ev.str.startswith("marker").to_numpy()
        is_ref = ct.isin(REFERENCE).to_numpy()
        samp = sub.obs["sample"].astype(str).to_numpy()
        samples = [s for s in pd.unique(samp)
                   if (is_marker_mal & (samp == s)).sum() >= MIN_MAL
                   and (is_ref & (samp == s)).sum() >= MIN_REF]
        if len(samples) < 2:
            del sub, X
            continue
        cons = {s: X[is_marker_mal & (samp == s)].mean(0) for s in samples}

        for s in samples:
            inS = samp == s
            others = [o for o in samples if o != s]
            own_sp = {}
            for name, half in halves.items():
                own = corr_rows(X[inS][:, half], cons[s][half])
                oth = np.mean([corr_rows(X[inS][:, half], cons[o][half]) for o in others], axis=0)
                own_sp[name] = (own, own - oth)
            o = sub.obs.loc[inS].copy()
            o["own_A"], o["spec_A"] = own_sp["A"]
            o["own_B"], o["spec_B"] = own_sp["B"]
            cls = pd.Series("other", index=o.index)
            cls[is_marker_mal[inS]] = "marker_malignant"
            cls[is_ref[inS]] = "reference"
            fl = o["cell_type"].astype(str).isin(FIBRO_LIKE).to_numpy() & ~is_marker_mal[inS]
            cls[fl] = "fibro_like"
            thrA = np.percentile(o.loc[is_ref[inS], "own_A"], 95)
            selA = (o["own_A"] > thrA).to_numpy()
            for grp in ["marker_malignant", "reference", "fibro_like"]:
                for sel, tag in [(np.ones(len(o), bool), "all"), (selA, "selA")]:
                    m = (cls == grp).to_numpy() & sel
                    if m.sum() < 10:
                        continue
                    rows.append(dict(study=study, sample=s, group=grp, subset=tag,
                                     n=int(m.sum()),
                                     spec_B=float(o.loc[m, "spec_B"].mean()),
                                     own_B=float(o.loc[m, "own_B"].mean())))
            cand = (cls == "fibro_like").to_numpy() & selA
            refm = (cls == "reference").to_numpy()
            if cand.sum() >= 10 and refm.sum() >= 20:
                y = np.r_[np.ones(cand.sum()), np.zeros(refm.sum())]
                v = np.r_[o.loc[cand, "spec_B"].to_numpy(), o.loc[refm, "spec_B"].to_numpy()]
                rows.append(dict(study=study, sample=s, group="fibro_like_vs_ref",
                                 subset="aucB", n=int(cand.sum()),
                                 spec_B=float(roc_auc_score(y, v)), own_B=np.nan))
            # keep per-cell values for downstream figure panels
            for col in ["own_A", "spec_A", "own_B", "spec_B"]:
                A.obs.loc[o.index, col] = o[col].to_numpy()
            A.obs.loc[o.index, "selA"] = selA
            A.obs.loc[o.index, "caf_class"] = cls.to_numpy()
        del sub, X

    df = pd.DataFrame(rows)
    tag = f"_ds{downsample}" if downsample else ""
    out = PROJ / f"results/tables/f4_split_genome_specificity{tag}.csv"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    A.obs.to_csv(PROJ / f"data/interim/f4_cell_specificity{tag}.csv")
    print(df.groupby(["group", "subset"]).agg(
        n_samples=("sample", "nunique"), cells=("n", "sum"),
        spec_B_median=("spec_B", "median")).round(4).to_string())
    return 0


def tables() -> int:
    """Panel inputs for Figure 4, from the two cached per-cell runs."""
    T = PROJ / "results/tables"
    full = pd.read_csv(PROJ / "data/interim/f4_cell_specificity.csv",
                       index_col=0, low_memory=False)
    ds = pd.read_csv(PROJ / "data/interim/f4_cell_specificity_ds1500.csv",
                     index_col=0, low_memory=False)
    for d in (full, ds):
        d["selA"] = d["selA"] == True  # noqa: E712 - csv round-trip gives strings
        d["grp"] = np.where((d["caf_class"] == "fibro_like") & d["selA"], "candidate",
                            np.where(d["caf_class"] == "fibro_like", "fibro-like CNA-",
                                     d["caf_class"].astype(str)))

    # (a) cell-level specificity at matched depth, capped per group so the
    #     violins are not dominated by the largest group's point count
    pa = ds.loc[ds["caf_class"].notna(), ["study", "sample", "grp", "spec_B"]].dropna()
    pa = pa[pa.grp.isin(["reference", "fibro-like CNA-", "candidate", "marker_malignant"])]
    pd.concat([g.sample(min(len(g), 4000), random_state=0) for _, g in pa.groupby("grp")],
              ignore_index=True).to_csv(T / "f4_specificity_cells.csv", index=False)

    # (c) transcriptional profile of the same groups, at native depth
    rows = []
    keep = full.grp.isin(["reference", "fibro-like CNA-", "candidate", "marker_malignant"])
    for g, d in full[keep].groupby("grp"):
        rows.append(dict(group={"reference": "reference", "fibro-like CNA-": "fibro-like, CNA-",
                                "candidate": "candidate (CNA+ fibro-like)",
                                "marker_malignant": "marker malignant"}[g], n=len(d),
                         NT5E_pct=100 * (d["NT5E"] > 0).mean(), NT5E_mean=d["NT5E"].mean(),
                         caf_program=d["caf_program"].mean(),
                         ef_repressed=d["ef_repressed"].mean(),
                         ef_activated=d["ef_activated"].mean(), cycle=d["G2M_score"].mean(),
                         genes=d["n_genes_by_counts"].median(),
                         doublet=d["doublet_score"].mean()))
    pd.DataFrame(rows).to_csv(T / "f4_group_expression.csv", index=False)

    # (d) prevalence within each sample's malignant compartment
    pv = []
    for (st, s), d in full.groupby(["study", "sample"]):
        cls = d["caf_class"].astype(str)
        comp = int(((cls == "fibro_like") | (cls == "marker_malignant")).sum())
        if comp == 0:
            continue
        cand = int(((cls == "fibro_like") & d["selA"]).sum())
        pv.append(dict(study=st, sample=s, patient=str(d["patient"].iloc[0]),
                       compartment=comp, fibro_like=int((cls == "fibro_like").sum()),
                       candidate=cand, pct=100 * cand / comp,
                       assessed=bool(d["spec_B"].notna().any())))
    pd.DataFrame(pv).to_csv(T / "f4_prevalence.csv", index=False)
    print("wrote f4_specificity_cells.csv, f4_group_expression.csv, f4_prevalence.csv")
    return 0


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else "0"
    raise SystemExit(tables() if arg == "tables" else main(int(arg)))
