"""On-target / off-tumour profile of five Ewing sarcoma cell-surface targets.

STEAP1, IL1RAP, CD276 (B7-H3), CD99 and NT5E (CD73) all have active CAR-T,
CAR-NK or antibody-drug-conjugate programmes in this disease. What no single
study can supply is the denominator: how often the target is also carried by the
normal cells sitting next to the tumour, across every patient and study at once.

Two tables are written.

  * The landscape: per cent of cells positive and mean expression for each
    target in each annotated cell type, pooled. This gives absolute levels, which
    a selectivity margin alone cannot - a target with a healthy margin but 2%
    positivity is useless.
  * The margin: the same quantity computed WITHIN each sample and differenced
    against the malignant cells of that sample. Pooled comparisons are
    confounded, because 54% of atlas cells are Waltner nuclei with a median of
    ~1,100 genes detected, and cell types are not contributed evenly by the
    studies; a paired within-sample contrast holds assay, depth and patient
    fixed.

Detection is RNA, not protein. Dropout means the per-cent-positive figures are
lower bounds on surface presentation, and a negative margin is better evidence
than a positive one is.
"""
import sys
import pathlib
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from ewsatlas import config  # noqa: E402

TARGETS = ["STEAP1", "IL1RAP", "CD276", "CD99", "NT5E"]
NORMAL = ["T_cell", "NK", "Monocyte_Macrophage", "B_cell", "Endothelial", "Fibroblast_CAF",
          "DC", "pDC", "Mast", "Plasma", "Urothelial"]
MIN_MAL, MIN_NORM = 50, 20


def main() -> int:
    cfg = config.load()
    PROJ = config.ROOT
    A = ad.read_h5ad(PROJ / config.get(cfg, "paths.atlas_out"))
    missing = [g for g in TARGETS if g not in A.var_names]
    if missing:
        raise SystemExit(f"targets absent from the gene space: {missing}")

    X = A[:, TARGETS].layers["lognorm"]
    X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
    E = pd.DataFrame(X, index=A.obs_names, columns=TARGETS)
    for c in ["cell_type", "sample", "study", "assay", "patient"]:
        E[c] = A.obs[c].astype(str).to_numpy()
    E["malignant"] = A.obs["malignant_frozen"].to_numpy().astype(bool)

    # landscape, pooled
    g = E.groupby("cell_type")[TARGETS]
    land = pd.concat([g.apply(lambda d: 100 * (d > 0).mean()).add_suffix("_pct"),
                      g.mean().add_suffix("_mean")], axis=1)
    land["n"] = E.groupby("cell_type").size()
    land.to_csv(PROJ / "results/tables/f6_target_landscape.csv")

    # paired within-sample margin
    rows = []
    for s, d in E.groupby("sample"):
        m = d[d.malignant]
        if len(m) < MIN_MAL:
            continue
        for ct in NORMAL:
            n = d[(~d.malignant) & (d.cell_type == ct)]
            if len(n) < MIN_NORM:
                continue
            for gene in TARGETS:
                rows.append(dict(sample=s, study=d["study"].iloc[0], assay=d["assay"].iloc[0],
                                 patient=d["patient"].iloc[0], normal_type=ct, gene=gene,
                                 n_mal=len(m), n_norm=len(n),
                                 mal_pct=100 * (m[gene] > 0).mean(),
                                 norm_pct=100 * (n[gene] > 0).mean()))
    R = pd.DataFrame(rows)
    R["delta"] = R.mal_pct - R.norm_pct
    R.to_csv(PROJ / "results/tables/f6_target_margins.csv", index=False)

    print(f"{E.malignant.sum():,} malignant cells; margins from "
          f"{R['sample'].nunique()} samples")
    print("\nmedian within-sample margin (malignant minus normal, percentage points):")
    print(R.pivot_table(index="gene", columns="normal_type", values="delta",
                        aggfunc="median").reindex(TARGETS).round(1).to_string())
    print("\npooled per cent positive in malignant cells:")
    print((100 * (E.loc[E.malignant, TARGETS] > 0).mean()).round(1).to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
