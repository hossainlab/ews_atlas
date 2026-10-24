"""Stage 04 - cluster the primary embedding and assign cell types.

    python src/pipeline/04_annotate.py

Attaches the declared primary embedding (scVI, batch_key='sample') to the full
feature space, builds the neighbour graph, clusters with leiden, and labels
clusters by marker enrichment.

Labelling is at CLUSTER level, not per cell. Per-cell marker argmax was tried
twice on this data and failed both times in the same way: in v4 it produced
4,403 spurious platelets, and on v5 it produced 12,365 platelets, 22 Ewing cells
in a bladder-Ewing study, and urothelium inside a PDX-only study. Marker
expression is too sparse per cell to argmax directly.

Two guards encode knowledge the markers cannot carry:

* urothelial is callable only where the bladder study (Mao) supplies at least
  half the cluster - urothelium exists in exactly one deposit, so a urothelial
  call anywhere else is a marker artefact;
* human immune labels are refused for clusters that are predominantly
  orthotopic PDX, where those compartments are mouse-host derived and cannot be
  human immune cells in a human-reference alignment.

The leiden flavour is pinned to ``leidenalg``: scanpy's ``igraph`` flavour
defaults to CPM and over-fragmented v4 into 304 communities at resolution 1.0.
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ewsatlas import config, io, provenance  # noqa: E402
from ewsatlas.compat import sequential_joblib  # noqa: E402

PANELS = {
    "Ewing_sarcoma": ["CD99", "NKX2-2", "NR0B1", "FLI1", "CAV1", "PRKCB", "EZH2", "STEAP1", "LIPI"],
    "Ewing_proliferating": ["MKI67", "TOP2A", "CCNB1", "UBE2C", "BIRC5", "CENPF"],
    "T_cell": ["CD3D", "CD3E", "TRAC", "IL7R", "CD2", "CD247"],
    "NK_cell": ["NKG7", "GNLY", "KLRD1", "PRF1", "KLRF1"],
    "B_cell": ["MS4A1", "CD79A", "CD19", "BANK1"],
    "Plasma": ["MZB1", "JCHAIN", "DERL3"],
    "Monocyte_Macrophage": ["LYZ", "CD14", "FCGR3A", "CD68", "AIF1", "C1QA", "C1QB", "TYROBP"],
    "DC": ["CLEC9A", "CD1C", "FCER1A"],
    "pDC": ["LILRA4", "GZMB", "IRF7"],
    "Endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5", "EGFL7"],
    "Fibroblast_CAF": ["COL1A1", "COL1A2", "DCN", "LUM", "PDGFRB", "COL3A1"],
    "MSC_Osteoblast": ["RUNX2", "SP7", "SPP1", "IBSP"],
    "Urothelial": ["KRT18", "KRT8", "UPK1B", "PSCA", "KRT19"],
    "Platelet": ["PPBP", "PF4", "ITGA2B"],
}
IMMUNE = {"T_cell", "NK_cell", "B_cell", "Plasma", "Monocyte_Macrophage", "DC", "pDC", "Platelet"}
BLADDER_STUDY = "Mao_2024"

#: Obligatory markers: a cluster may only take a label if at least one of these
#: is expressed above `min_mean` (cluster-mean log-normalized expression). A
#: z-score says "this panel is relatively higher here than elsewhere"; it does
#: NOT say the panel is expressed at all, and on this data the difference
#: matters. Three calls in the first annotated build were wrong for exactly that
#: reason: a 7,739-cell cluster took NK on NKG7/GNLY while expressing CD3D/CD3E/
#: TRAC (cytotoxic T, which shares the cytotoxic panel), a 5,763-cell cluster
#: took DC with LYZ at 4.05 and CD1C at 0.20 (a monocyte), and four clusters
#: with essentially no lineage markers at all took immune labels on relative
#: enrichment alone.
REQUIRES = {
    "Ewing_sarcoma": (["CD99", "FLI1"], 0.3),
    "Ewing_proliferating": (["MKI67", "TOP2A"], 0.3),
    "T_cell": (["CD3D", "CD3E", "TRAC"], 0.25),
    "NK_cell": (["NKG7", "GNLY"], 0.5),
    "B_cell": (["MS4A1", "CD79A"], 0.25),
    "Plasma": (["MZB1", "JCHAIN"], 0.5),
    "Monocyte_Macrophage": (["LYZ", "CD68", "AIF1", "C1QA"], 0.5),
    "DC": (["CD1C", "CLEC9A", "FCER1A"], 0.3),
    "pDC": (["LILRA4", "IRF7"], 0.3),
    "Endothelial": (["PECAM1", "CDH5"], 0.3),
    "Fibroblast_CAF": (["COL1A1", "COL1A2", "DCN"], 0.5),
    # SPP1 is deliberately NOT obligatory here: macrophages express it, and
    # MSC/osteoblast vs SPP1+ macrophage is the known ambiguity in this atlas.
    "MSC_Osteoblast": (["RUNX2", "SP7", "IBSP"], 0.2),
    "Urothelial": (["KRT18", "KRT8"], 0.5),
    "Platelet": (["PPBP", "PF4"], 0.5),
}

#: Every NUCLEATED immune label additionally requires CD45. A cluster that does
#: not express PTPRC is not a leukocyte, whatever its relative panel score.
#: Platelets are exempt because they are anucleate and genuinely PTPRC-negative
#: - the first version of this rule rejected a clean PPBP/PF4 cluster and let it
#: fall through to Endothelial.
PTPRC_MIN = 0.2
PTPRC_EXEMPT = {"Platelet"}

#: A label reached only after other panels were refused still needs positive
#: evidence of its own. Below this cluster-level z the call is not made: four
#: clusters with no lineage markers at all otherwise inherited whichever panel
#: happened to rank first.
MIN_ACCEPT_Z = 0.5

#: Cytotoxic T cells express the NK panel. When a cluster scores NK highest but
#: carries comparable CD3, it is T.
CD3_FOR_NK_VETO = 0.5


def main() -> int:
    cfg = config.load()
    interim = config.path(cfg, "interim")
    src = interim / "atlas_preintegration.h5ad"
    emb = interim / "embeddings" / "latent_scvi_sample.npy"
    cells_f = interim / "embeddings" / "cell_index.npy"

    with provenance.start_run("04_annotate", cfg, [src, emb]) as run:
        A = ad.read_h5ad(src)
        order = np.load(cells_f, allow_pickle=True).astype(str)
        if list(A.obs_names) != list(order):
            raise AssertionError("cell order differs between the atlas object and the embedding")
        A.obsm["X_scvi"] = np.load(emb).astype(np.float32)

        res = float(config.get(cfg, "annotation.cluster.resolution"))
        with sequential_joblib():
            sc.pp.neighbors(A, use_rep="X_scvi", n_neighbors=20, random_state=int(cfg["seed"]))
            sc.tl.umap(A, min_dist=0.6, random_state=int(cfg["seed"]))
        sc.tl.leiden(A, resolution=res, key_added="leiden",
                     flavor=str(config.get(cfg, "annotation.cluster.flavor")),
                     random_state=int(cfg["seed"]), n_iterations=2)
        n_clusters = int(A.obs["leiden"].nunique())
        print(f"leiden: {n_clusters} clusters at resolution {res}", flush=True)

        # cluster-level marker scores on log-normalized expression
        A.X = A.layers["lognorm"].copy()
        for name, genes in PANELS.items():
            present = [g for g in genes if g in A.var_names]
            sc.tl.score_genes(A, present, score_name=f"sc_{name}",
                              ctrl_size=50, random_state=int(cfg["seed"]))
        score_cols = [f"sc_{k}" for k in PANELS]
        per_cluster = A.obs.groupby("leiden", observed=True)[score_cols].mean()
        z = (per_cluster - per_cluster.mean()) / per_cluster.std()

        comp = pd.crosstab(A.obs["leiden"], A.obs["study"], normalize="index")
        pdx = pd.crosstab(A.obs["leiden"], A.obs["material_type"], normalize="index").get(
            "O-PDX", pd.Series(0.0, index=per_cluster.index))

        # cluster-mean expression of every gene the guards consult
        guard_genes = sorted({g for gs, _ in REQUIRES.values() for g in gs}
                             | {"PTPRC", "CD3D", "CD3E", "TRAC"})
        guard_genes = [g for g in guard_genes if g in A.var_names]
        expr = pd.DataFrame(
            np.asarray(A[:, guard_genes].X.todense()), columns=guard_genes, index=A.obs_names)
        expr["leiden"] = A.obs["leiden"].to_numpy()
        gmean = expr.groupby("leiden", observed=True).mean()

        calls, notes = {}, {}
        for cl in per_cluster.index:
            ranked = z.loc[cl].sort_values(ascending=False)
            reasons = []
            for col in ranked.index:
                label = col[3:]
                if label == "Urothelial" and comp.loc[cl].get(BLADDER_STUDY, 0.0) < 0.5:
                    reasons.append("urothelial refused (bladder study <50% of cluster)")
                    continue
                if label in IMMUNE and pdx.get(cl, 0.0) > 0.8:
                    reasons.append("immune refused (cluster predominantly O-PDX)")
                    continue
                need, floor = REQUIRES[label]
                have = [g for g in need if g in gmean.columns]
                if have and gmean.loc[cl, have].max() < floor:
                    reasons.append(f"{label} refused (markers absent: "
                                   f"max {gmean.loc[cl, have].max():.2f} < {floor})")
                    continue
                if (label in IMMUNE and label not in PTPRC_EXEMPT
                        and gmean.loc[cl].get("PTPRC", 0.0) < PTPRC_MIN):
                    reasons.append(f"{label} refused (PTPRC {gmean.loc[cl].get('PTPRC', 0.0):.2f} "
                                   f"< {PTPRC_MIN})")
                    continue
                if label == "NK_cell":
                    cd3 = gmean.loc[cl, [g for g in ("CD3D", "CD3E", "TRAC") if g in gmean.columns]].max()
                    if cd3 >= CD3_FOR_NK_VETO:
                        reasons.append(f"NK refused (CD3 {cd3:.2f} - cytotoxic T, not NK)")
                        continue
                if ranked[col] < MIN_ACCEPT_Z:
                    reasons.append(f"{label} refused (cluster z {ranked[col]:.2f} "
                                   f"< {MIN_ACCEPT_Z}, no positive evidence)")
                    calls[cl] = "Unassigned"
                    break
                calls[cl] = label
                break
            else:
                calls[cl] = "Unassigned"
            if reasons:
                notes[cl] = "; ".join(reasons[:3])
        A.obs["cell_type"] = A.obs["leiden"].map(calls).astype(str)

        summary = pd.DataFrame({
            "leiden": per_cluster.index,
            "cells": A.obs["leiden"].value_counts().reindex(per_cluster.index).to_numpy(),
            "cell_type": [calls[c] for c in per_cluster.index],
            "top_z": z.max(axis=1).round(2).to_numpy(),
            "n_studies": (comp > 0.01).sum(axis=1).to_numpy(),
            "pct_pdx": (100 * pdx.reindex(per_cluster.index).fillna(0)).round(1).to_numpy(),
            "guard_note": [notes.get(c, "") for c in per_cluster.index],
        })
        out_tbl = config.ROOT / "results" / "tables" / "v5_cluster_annotation.csv"
        summary.to_csv(out_tbl, index=False)

        A.X = A.layers["counts"].copy()
        out = config.ROOT / config.get(cfg, "paths.atlas_out")
        out.parent.mkdir(parents=True, exist_ok=True)
        io.write_h5ad_atomic(A, out)
        run.add(n_clusters=n_clusters,
                cell_types=A.obs["cell_type"].value_counts().to_dict(),
                guards_fired=len(notes), out_bytes=int(out.stat().st_size))
        print(A.obs["cell_type"].value_counts().to_string(), flush=True)
        print(f"wrote {out.name} ({out.stat().st_size/1e6:.0f} MB)", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
