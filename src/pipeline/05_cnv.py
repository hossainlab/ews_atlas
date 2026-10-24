"""Stage 05 - copy-number inference and the frozen malignant set.

    python src/pipeline/05_cnv.py

Why this stage exists: the discovery question concerns malignant cells with LOW
expression of the EWSR1-FLI1 programme. A marker-based definition of "malignant"
leans on exactly that programme (CD99, FLI1, NKX2-2), so it preferentially
excludes the cells under study. Copy-number evidence is independent of the
transcriptional programme and breaks that circularity.

WHAT CHANGED AND WHY (v5 -> v6)
-------------------------------
v5 scored each cell by the mean absolute inferred CNV across the genome and
failed: malignant-vs-reference AUC had a per-sample median of 0.540. Three
candidate explanations were tested against measurement rather than argued:

1. *Window size and stride.* Adopting NBAtlas's window_size=250 (from 100)
   changed the genome-wide AUC by -0.014 in Cillo and +0.027 in Visser. Not the
   cause. Stride makes no difference at all: step=10 against step=1 gives a
   per-sample median correlation AUC of 0.732 against 0.727, Spearman 0.999 and
   a maximum absolute per-sample difference of 0.015 - so step=10 is used, which
   is ten times cheaper and lets the 74,888-nucleus study run in memory.
2. *Reference width.* Widening the reference from immune-only to include
   endothelial and urothelial cells raised the Visser genome-wide AUC from 0.635
   to 0.662 and left Cillo unchanged (that study has too few of either type to
   qualify). Marginal.
3. *The statistic itself.* This was the cause. The mean of |CNV| discards the
   sign and the spatial pattern of the alterations, which is where the signal
   is: per chromosome arm the signed scores separate cleanly and in the expected
   directions (chr8q 0.61-0.72 and chr12q 0.64-0.70 above 0.5, chr16q 0.38-0.43
   below it - gain, gain, loss).

So this stage computes the two statistics of Tirosh et al. (Science 352:189,
2016) per sample:

* ``cnv_signal``  - mean |CNV| across windows, as before; retained because it is
  the quantity v5 reported and Figure 3 of the published set discusses.
* ``cnv_corr``    - correlation of each cell's CNV profile with a per-sample
  consensus malignant profile. This lifts the per-sample median AUC from 0.540
  to 0.727 and leaves no sample below 0.5 (previously 13 of 37 were).

The consensus is seeded from the high-signal quartile of cells that are NOT
annotated as reference types. Seeding on reference identity is deliberate:
immune and endothelial identity does not depend on the fusion programme, so the
seed is not circular with respect to the cell state under study. Seeding on the
top-signal quartile of *all* cells - the obvious first implementation - inverts
in samples where malignant cells are a minority: in 038_CTC (64 malignant
against 3,386 reference) the consensus became the reference profile and the AUC
read 0.368.

Classification is per sample, as in NBAtlas, with the threshold set from the
reference cells of that same sample (95th percentile of their ``cnv_corr``) so
no global cut-off is imposed across libraries of different depth.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import infercnvpy as cnv
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ewsatlas import compat, config, io, provenance  # noqa: E402

REFERENCE = ["T_cell", "NK", "B_cell", "Plasma", "Monocyte_Macrophage", "DC", "pDC", "Mast",
             "Platelet", "Endothelial", "Urothelial"]
MALIGNANT_MARKER_CALL = ["Ewing_sarcoma", "Ewing_proliferating"]
MIN_REF_CELLS = 50
MIN_MAL_CELLS = 20
MIN_TYPE_CELLS = 20
WINDOW, STEP = 250, 10
# Recurrent arm-level changes in Ewing sarcoma, scored signed.
ARMS = ["chr1p", "chr1q", "chr8p", "chr8q", "chr9p", "chr12p", "chr12q", "chr16q", "chr17p"]


def _arm_of_window(sub, window_size: int, step: int, cen: dict) -> np.ndarray:
    """Label each X_cnv column with its chromosome arm, in infercnvpy's gene order."""
    chr_pos = sub.uns["cnv"]["chr_pos"]
    order = sorted(chr_pos.items(), key=lambda kv: kv[1])
    ncol = sub.obsm["X_cnv"].shape[1]
    arms = np.array([""] * ncol, dtype=object)
    for i, (chrom, c0) in enumerate(order):
        c1 = order[i + 1][1] if i + 1 < len(order) else ncol
        starts = sub.var.loc[sub.var["chromosome"] == chrom, "start"].sort_values().to_numpy()
        cs, ce = cen.get(chrom, (np.inf, np.inf))
        for j in range(c1 - c0):
            sl = starts[j * step: j * step + window_size]
            if len(sl) == 0:
                continue
            mid = float(np.median(sl))
            arms[c0 + j] = f"{chrom}p" if mid < cs else (f"{chrom}q" if mid > ce else f"{chrom}cen")
    return arms


def main() -> int:
    cfg = config.load()
    atlas = config.ROOT / config.get(cfg, "paths.atlas_out")
    ref_dir = config.ROOT / "data" / "external" / "references"
    pos_file, cen_file = ref_dir / "gene_positions.csv", ref_dir / "hg38_centromeres.json"

    with provenance.start_run("05_cnv", cfg, [atlas, pos_file, cen_file]) as run:
        A = ad.read_h5ad(atlas)
        pos = pd.read_csv(pos_file, index_col=0)
        cen = {k: tuple(v) for k, v in json.load(open(cen_file)).items()}
        for c in ("chromosome", "start", "end"):
            A.var[c] = pd.Series(pos[c]).reindex(A.var_names).to_numpy()
        n_coord = int(A.var["chromosome"].notna().sum())
        print(f"genes with coordinates: {n_coord:,}/{A.n_vars:,}", flush=True)

        signal = np.full(A.n_obs, np.nan)
        corr = np.full(A.n_obs, np.nan)
        arm_mat = {a: np.full(A.n_obs, np.nan) for a in ARMS}
        klass = np.array(["not_run"] * A.n_obs, dtype=object)
        per_sample = []

        study_of = A.obs["study"].astype(str).to_numpy()
        samp_of = A.obs["sample"].astype(str).to_numpy()
        ct_of = A.obs["cell_type"].astype(str).to_numpy()

        for study in pd.unique(study_of):
            sidx = np.where(study_of == study)[0]
            sub = A[sidx].copy()
            sub.X = sub.layers["lognorm"].copy()
            sub.var.index = pd.Index(list(sub.var_names), dtype=object, name="gene")
            ct = pd.Series(ct_of[sidx])
            ref = [c for c in REFERENCE if (ct == c).sum() >= MIN_TYPE_CELLS]
            if not ref:
                for s in pd.unique(samp_of[sidx]):
                    per_sample.append({"study": study, "sample": s, "status": "no_reference_in_study"})
                print(f"{study}: no reference cell type reaches {MIN_TYPE_CELLS} cells", flush=True)
                del sub
                continue
            # n_jobs=1 is required, not a tuning choice: infercnvpy fans the
            # chromosome loop out over a ProcessPoolExecutor, and this sandbox
            # denies the named pipe that the pool needs when the stage runs as a
            # subprocess ("PermissionError: [WinError 5] Access is denied").
            # The same call succeeds in-kernel, which is what made the defect
            # look intermittent.
            with compat.sequential_process_map():
                cnv.tl.infercnv(sub, reference_key="cell_type", reference_cat=ref,
                                window_size=WINDOW, step=STEP, n_jobs=1)
            X = sub.obsm["X_cnv"]
            X = np.asarray(X.todense()) if hasattr(X, "todense") else np.asarray(X)
            arms = _arm_of_window(sub, WINDOW, STEP, cen)
            arm_cols = {a: np.where(arms == a)[0] for a in ARMS}

            is_ref = ct.isin(ref).to_numpy()
            is_mal = ct.isin(MALIGNANT_MARKER_CALL).to_numpy()
            print(f"{study}: {len(sidx):,} cells, {int(is_ref.sum()):,} reference, "
                  f"{X.shape[1]} windows, reference types {';'.join(ref)}", flush=True)

            for s in pd.unique(samp_of[sidx]):
                m = samp_of[sidx] == s
                g = sidx[m]
                sig = np.abs(X[m]).mean(axis=1)
                signal[g] = sig
                for a, cols in arm_cols.items():
                    if len(cols) >= 5:
                        arm_mat[a][g] = X[np.ix_(m, cols)].mean(axis=1)
                y, r = is_mal[m], is_ref[m]
                row = {"study": study, "sample": s, "cells": int(m.sum()),
                       "malignant_marker": int(y.sum()), "reference_cells": int(r.sum())}
                if y.sum() < MIN_MAL_CELLS or r.sum() < MIN_REF_CELLS:
                    row["status"] = "insufficient_cells"
                    per_sample.append(row)
                    continue
                cand = ~r
                thr = np.quantile(sig[cand], 0.75) if cand.sum() >= 8 else np.quantile(sig, 0.75)
                top = cand & (sig >= thr)
                if top.sum() < 5:
                    top = cand
                consensus = X[m][top].mean(axis=0)
                cc = X[m] - X[m].mean(axis=1, keepdims=True)
                cv = consensus - consensus.mean()
                cr = (cc @ cv) / (np.sqrt((cc ** 2).sum(axis=1)) * np.sqrt((cv ** 2).sum()) + 1e-12)
                corr[g] = cr
                cut = float(np.quantile(cr[r], 0.95))
                klass[g] = np.where(cr > cut, "aneuploid", "diploid")
                keep = y | r
                row.update({"status": "scored",
                            "auc_signal": round(float(roc_auc_score(y[keep], sig[keep])), 3),
                            "auc_corr": round(float(roc_auc_score(y[keep], cr[keep])), 3),
                            "corr_cutoff": round(cut, 4),
                            "aneuploid": int((cr > cut).sum()),
                            "aneuploid_of_marker_malignant": int(((cr > cut) & y).sum()),
                            "reference_types": ";".join(ref)})
                per_sample.append(row)
            del sub, X

        A.obs["cnv_signal"] = signal
        A.obs["cnv_corr"] = corr
        A.obs["cnv_class"] = pd.Categorical(klass)
        for a in ARMS:
            A.obs[f"cnv_{a}"] = arm_mat[a]
        # keep the v5 column name pointing at the quantity Figure 3 discusses
        A.obs["cnv_score"] = signal

        # ------------------------------------------------------------------
        # The frozen malignant set.
        #
        # Marker evidence and copy-number evidence fail in opposite directions
        # here, which is why the set is the union of the two rather than either
        # alone:
        #
        # * Markers (CD99, NKX2-2, FLI1) ARE the fusion programme, so a
        #   marker-defined set excludes low-fusion cells by construction - the
        #   exact cells the discovery question is about.
        # * The copy-number call ranks well (per-sample median AUC 0.712) but is
        #   deliberately thresholded at the 95th percentile of each sample's own
        #   reference cells, which makes it insensitive: only ~15% of
        #   marker-malignant cells clear it. That is acceptable HERE because the
        #   call is used only to RESCUE cells the markers missed, never to
        #   remove cells the markers found. Every cell it adds carries strong
        #   copy-number evidence.
        #
        # Aneuploid calls on reference cell types are treated as false positives
        # rather than rescues: immune and endothelial identity does not depend
        # on the fusion programme, so PTPRC/CD3/PECAM1 evidence outranks a
        # 0.712-AUC copy-number call.
        marker_mal = pd.Series(ct_of, index=A.obs_names).isin(MALIGNANT_MARKER_CALL).to_numpy()
        is_ref_type = pd.Series(ct_of, index=A.obs_names).isin(REFERENCE).to_numpy()
        aneuploid = klass == "aneuploid"
        rescued = aneuploid & ~marker_mal & ~is_ref_type
        frozen = marker_mal | rescued

        evidence = np.where(marker_mal & aneuploid, "marker_and_cna",
                   np.where(marker_mal & (klass == "diploid"), "marker_only_cna_negative",
                   np.where(marker_mal, "marker_only_cna_unassessed",
                   np.where(rescued, "cna_only", "not_malignant"))))
        A.obs["malignant_frozen"] = frozen
        A.obs["malignant_evidence"] = pd.Categorical(evidence)
        freeze = {
            "frozen_malignant": int(frozen.sum()),
            "marker_malignant": int(marker_mal.sum()),
            "rescued_by_cna": int(rescued.sum()),
            "rescued_breakdown": pd.Series(ct_of[rescued]).value_counts().to_dict(),
            "aneuploid_on_reference_types_discarded": int((aneuploid & is_ref_type).sum()),
            "marker_malignant_cna_unassessed": int((marker_mal & (klass == "not_run")).sum()),
            "marker_malignant_cna_negative": int((marker_mal & (klass == "diploid")).sum()),
        }
        print(json.dumps(freeze, indent=1, default=str), flush=True)
        pd.DataFrame({"cell_type": ct_of, "evidence": evidence}).value_counts().reset_index(
            name="cells").to_csv(config.ROOT / "results" / "tables" / "v6_malignant_evidence.csv", index=False)

        tbl = pd.DataFrame(per_sample)
        tbl.to_csv(config.ROOT / "results" / "tables" / "v6_cnv_per_sample.csv", index=False)
        sc = tbl[tbl.get("status", pd.Series(dtype=str)) == "scored"]
        summary = {
            "genes_with_coordinates": n_coord,
            "window_size": WINDOW, "step": STEP,
            "samples_scored": int(len(sc)),
            "samples_not_scored": int(len(tbl) - len(sc)),
            "cells_scored": int(np.isfinite(corr).sum()),
            "auc_signal_median": float(sc["auc_signal"].median()) if len(sc) else float("nan"),
            "auc_corr_median": float(sc["auc_corr"].median()) if len(sc) else float("nan"),
            "samples_auc_corr_above_0.7": int((sc["auc_corr"] > 0.7).sum()) if len(sc) else 0,
            "samples_auc_corr_below_0.5": int((sc["auc_corr"] < 0.5).sum()) if len(sc) else 0,
            "aneuploid_cells": int((klass == "aneuploid").sum()),
            **freeze,
        }
        print(json.dumps(summary, indent=1), flush=True)
        io.write_h5ad_atomic(A, atlas)
        run.add(**summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
