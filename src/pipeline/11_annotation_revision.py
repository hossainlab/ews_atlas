"""Stage 11 - apply the immune re-analysis corrections to the released object.

    python src/pipeline/11_annotation_revision.py

The atlas-level ``cell_type`` column was assigned at the level of 34 Leiden
communities over all 138,353 cells. Re-clustering the 29,114 immune cells on
their own neighbour graph (docs/10) found four errors that a coarse pass cannot
avoid, and until now they shipped on the released object:

  * all 715 cells labelled ``DC`` are MAST cells (TPSAB1 3.76, CPA3 5.80);
  * 1,267 NK cells sit inside ``T_cell`` (1,266) and ``Monocyte_Macrophage`` (1);
  * the real conventional dendritic cells - 1,137 cDC2 - sit inside
    ``Monocyte_Macrophage``, so the atlas had no dendritic-cell class at all;
  * 307 plasmacytoid dendritic cells sit inside ``Plasma``.

Twelve classes become fifteen. ``Low_quality`` immune cells (782) keep their
atlas label: a quality flag is not a cell type, and `immune_subset` already
records it.

The copy-number reference set is NOT changed by this. Mast and NK are normal
diploid classes and are added to ``05_cnv.REFERENCE``, so the same cells serve as
reference before and after; the stage asserts that set equality and refuses to
write if it fails, which is what makes it safe to leave the cnv_* columns and the
frozen malignant set untouched.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ewsatlas import config, io, provenance  # noqa: E402

# subset label -> revised atlas label, for subsets whose atlas label is wrong
REVISION = {"Mast": "Mast", "NK": "NK", "cDC2": "DC", "pDC": "pDC"}
REFERENCE_BEFORE = ["T_cell", "B_cell", "Plasma", "Monocyte_Macrophage", "DC", "pDC",
                    "Platelet", "Endothelial", "Urothelial"]
REFERENCE_AFTER = REFERENCE_BEFORE + ["Mast", "NK"]


def main() -> int:
    cfg = config.load()
    PROJ = config.ROOT
    atlas = PROJ / config.get(cfg, "paths.atlas_out")
    labels = PROJ / "results/tables/v6_immune_labels.csv"

    with provenance.start_run("11_annotation_revision", cfg, [atlas, labels]) as run:
        A = ad.read_h5ad(atlas)
        lab = pd.read_csv(labels, index_col=0)["immune_subset"]
        before = A.obs["cell_type"].astype(str).copy()

        if "cell_type_v5_original" not in A.obs:
            A.obs["cell_type_v5_original"] = before.to_numpy()
        new = before.copy()
        moved = {}
        for subset, target in REVISION.items():
            cells = lab.index[lab == subset]
            cells = [c for c in cells if c in A.obs_names]
            moved[f"{subset} -> {target}"] = len(cells)
            new.loc[cells] = target
        A.obs["cell_type"] = pd.Categorical(new.to_numpy())
        A.obs["immune_subset"] = pd.Categorical(
            lab.reindex(A.obs_names).fillna("not_immune").to_numpy())

        ref_before = set(A.obs_names[before.isin(REFERENCE_BEFORE).to_numpy()])
        ref_after = set(A.obs_names[new.isin(REFERENCE_AFTER).to_numpy()])
        if ref_before != ref_after:
            only_b, only_a = len(ref_before - ref_after), len(ref_after - ref_before)
            raise SystemExit(
                f"copy-number reference set changed ({only_b} lost, {only_a} gained); "
                "the cnv_* columns would have to be recomputed - refusing to write")

        A.uns["cell_type_revision"] = {
            "applied": "stage 11_annotation_revision",
            "source": "immune re-clustering, docs/10; labels in v6_immune_labels.csv",
            "moves": moved,
            "note": ("cell_type_v5_original preserves the labels the integration and "
                     "copy-number stages ran on; the copy-number reference set is "
                     "unchanged, so cnv_* and malignant_frozen are still valid."),
        }
        io.write_h5ad_atomic(A, atlas)

        counts = A.obs["cell_type"].value_counts()
        run.add(moves=moved, n_classes=int(counts.size),
                reference_set_unchanged=True, cells=int(A.n_obs))
        print("moves:", moved)
        print(f"\n{before.nunique()} classes -> {counts.size} classes")
        print(counts.to_string())
        print("\ncopy-number reference set unchanged:", len(ref_before), "cells")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
