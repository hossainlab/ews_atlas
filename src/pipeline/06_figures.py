"""Stage 06 - the publication figure set.

    python src/pipeline/06_figures.py

Three main figures and four supplementary figures, composed on a 12-column grid
at 180 mm (double-column) width, 300 dpi. **No panel carries a title** - the
captions in results/reports/figure_captions.md carry the message, which is why
every panel must be readable from its axis labels, direct labels and legend
alone.

Panels are rendered individually to transparent PNGs at exact pixel sizes and
tiled by `compose_figure` from the `figure-composer` skill, which also stamps
the panel letters. Panel renderers live in `ewsatlas.figures` so a panel can be
regenerated on its own without redrawing a figure that already passes.

Every number drawn is read from the atlas object or a results/tables/ file
produced by an earlier stage; nothing is retyped.
"""

from __future__ import annotations

import json
import sys
import os
import shutil
import time
from pathlib import Path

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ewsatlas import config, figures as F, provenance  # noqa: E402

WIDTH_MM = 180.0
DPI = 300
GUTTER_MM = 4.0
LETTER_PT = 9.0   # bold panel letter, points
LETTER_BAND_MM = 3.8  # clean strip above each row for the letter
MM = DPI / 25.4


def _style():
    mpl.rcParams.update({
        "figure.dpi": DPI, "savefig.dpi": DPI, "font.size": 7,
        "axes.titlesize": 7, "axes.labelsize": 6.6, "legend.fontsize": 5.6,
        "xtick.labelsize": 5.6, "ytick.labelsize": 5.6,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.2, "ytick.major.size": 2.2,
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": "none", "savefig.transparent": True,
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })


def _panel_px(outline, letter):
    """Pixel size of one panel's grid slot."""
    ncol = outline["ncol"]
    total_w = WIDTH_MM * MM
    gut = GUTTER_MM * MM
    colw = (total_w - gut * (ncol - 1)) / ncol
    p = next(x for x in outline["panels"] if x["letter"] == letter)
    w = colw * p["colspan"] + gut * (p["colspan"] - 1)
    h = outline["row_heights_mm"][p["row"]] * MM
    return int(round(w)), int(round(h))


def render(outline, letter, draw, outdir: Path, **kw) -> Path:
    """Render one panel at exactly its grid-slot size, transparent, no bbox trim."""
    w, h = _panel_px(outline, letter)
    fig = plt.figure(figsize=(w / DPI, h / DPI), dpi=DPI)
    ax = fig.add_subplot(111)
    draw(ax, **kw)
    fig.tight_layout(pad=0.45)
    p = outdir / f"{outline['tag']}_{letter}.png"
    fig.savefig(p, dpi=DPI, transparent=True)
    plt.close(fig)
    return p


def compose(outline, paths, out: Path) -> Path:
    """Tile panels onto the grid and stamp bold panel letters."""
    from PIL import Image, ImageDraw, ImageFont

    ncol = outline["ncol"]
    gut = int(round(GUTTER_MM * MM))
    total_w = int(round(WIDTH_MM * MM))
    colw = (total_w - gut * (ncol - 1)) / ncol
    # Reserve a band above each row for the panel letter. Stamping the letter
    # over the panel collides with whatever is at its top-left corner (a tall
    # rotated y-label, for instance); adding the band keeps every letter in
    # clean space without shrinking any panel.
    band = int(round(LETTER_BAND_MM * MM)) if len(outline["panels"]) > 1 else 0
    row_y, y = [], 0
    for hmm in outline["row_heights_mm"]:
        y += band
        row_y.append(y)
        y += int(round(hmm * MM)) + gut
    total_h = y - gut
    canvas = Image.new("RGBA", (total_w, total_h), (255, 255, 255, 255))
    # ImageFont sizes are PIXELS, so the point size has to be scaled by the
    # output resolution: 14 px at 300 dpi is 3.4 pt, far below the ~8-9 pt
    # bold that journals ask for. Deriving it keeps the letters correct if
    # DPI or WIDTH_MM change.
    px = int(round(LETTER_PT / 72.0 * DPI))
    try:
        font = ImageFont.truetype("DejaVuSans-Bold.ttf", px)
    except OSError:
        font = ImageFont.load_default(size=px)
    drw = ImageDraw.Draw(canvas)
    for p in outline["panels"]:
        img = Image.open(paths[p["letter"]]).convert("RGBA")
        x = int(round(p["col"] * (colw + gut)))
        canvas.alpha_composite(img, (x, row_y[p["row"]]))
        if band:
            # Cell Press house style sets panel letters in bold UPPERCASE.
            drw.text((x, row_y[p["row"]] - band), p["letter"].upper(),
                     fill=(0, 0, 0, 255), font=font, anchor="la")
    # Write via a temporary file: a transient lock on the target (a viewer, an
    # indexer) otherwise aborts the whole stage part way through the set.
    tmp = out.with_name(out.stem + ".tmp.png")
    canvas.convert("RGB").save(tmp, dpi=(DPI, DPI))
    # A viewer or indexer holding the previous PNG makes os.replace fail with
    # WinError 5; retry briefly rather than abandoning a whole figure set.
    for attempt in range(4):
        try:
            os.replace(tmp, out)
            break
        except PermissionError:
            if attempt < 3:
                time.sleep(0.5)
                continue
            # A file open in a viewer can often still be written in place even
            # though it cannot be renamed over (Windows share modes allow
            # read/write but deny delete). Fall back to an in-place copy so one
            # open figure does not abort the whole set.
            shutil.copyfile(tmp, out)
            tmp.unlink(missing_ok=True)
    return out



def marker_dotplot(atlas_path) -> pd.DataFrame:
    """Fraction expressing and mean log-normalized expression per cell type."""
    import scipy.sparse as sp
    genes = [g for _, grp in F.MARKER_GROUPS for g in grp]
    B = ad.read_h5ad(atlas_path, backed="r")
    present = [g for g in genes if g in B.var_names]
    sub = B[:, present].to_memory()
    B.file.close()
    X = sub.layers["lognorm"]
    X = X.tocsr() if sp.issparse(X) else X
    ct = sub.obs["cell_type"].astype(str).to_numpy()
    rows = []
    for c in pd.unique(ct):
        m = ct == c
        Xi = X[m]
        frac = np.asarray((Xi > 0).sum(0)).ravel() / int(m.sum())
        mean = np.asarray(Xi.mean(0)).ravel()
        for g, fr, mn in zip(present, frac, mean):
            rows.append(dict(cell_type=c, gene=g, frac_expressing=float(fr),
                             mean_lognorm=float(mn)))
    return pd.DataFrame(rows)

def build_outlines():
    # Figure 1 follows the NBAtlas Figure 1 layout: cohort schematic, stacked
    # metadata embeddings, the ringed atlas, marker evidence.
    # Figure 1 was 235 mm tall, taller than the journal's type area leaves room
    # for a caption. The three stacked metadata embeddings were dropped: study
    # and assay are already a column header and a glyph row of panel (a) and a
    # band of the ringed atlas, and the sample embedding is Figure S6 and said
    # of itself "too many to read". The ringed atlas now spans the full width
    # with its keys beside it, so the row is set by the circle's HEIGHT.
    main1 = {"tag": "fig1", "ncol": 12, "row_heights_mm": [58, 86, 46], "panels": [
        {"letter": "a", "row": 0, "col": 0, "colspan": 12},
        {"letter": "b", "row": 1, "col": 0, "colspan": 12},
        {"letter": "c", "row": 2, "col": 0, "colspan": 12}]}
    main2 = {"tag": "fig2", "ncol": 12, "row_heights_mm": [54, 54], "panels": [
        {"letter": "a", "row": 0, "col": 0, "colspan": 5},
        {"letter": "b", "row": 0, "col": 5, "colspan": 7},
        {"letter": "c", "row": 1, "col": 0, "colspan": 6},
        {"letter": "d", "row": 1, "col": 6, "colspan": 6}]}
    main3 = {"tag": "fig3", "ncol": 12, "row_heights_mm": [52, 26], "panels": [
        {"letter": "a", "row": 0, "col": 0, "colspan": 5},
        {"letter": "b", "row": 0, "col": 5, "colspan": 7},
        {"letter": "c", "row": 1, "col": 0, "colspan": 12}]}
    main5 = {"tag": "fig5", "ncol": 12, "row_heights_mm": [62, 38], "panels": [
        {"letter": "a", "row": 0, "col": 0, "colspan": 5},
        {"letter": "b", "row": 0, "col": 5, "colspan": 7},
        {"letter": "c", "row": 1, "col": 0, "colspan": 12}]}
    # Figure 4 (CAF-like state): the hypothesis test, one panel per line of
    # evidence - genetic (a, b), transcriptional (c), prevalence (d).
    # retained as Figure S7 after the duplicate "included datasets" panel was dropped
    main4b = {"tag": "figS7caf", "ncol": 12, "row_heights_mm": [52, 48], "panels": [
        {"letter": "a", "row": 0, "col": 0, "colspan": 5},
        {"letter": "b", "row": 0, "col": 5, "colspan": 7},
        {"letter": "c", "row": 1, "col": 0, "colspan": 7},
        {"letter": "d", "row": 1, "col": 7, "colspan": 5}]}
    # Figure 6: three independent demonstrations that the atlas is usable as a
    # reference - surface-target selectivity (a, b), label transfer (c), bulk
    # deconvolution with outcome (d).
    main6 = {"tag": "fig6", "ncol": 12, "row_heights_mm": [43, 46], "panels": [
        {"letter": "a", "row": 0, "col": 0, "colspan": 7},
        {"letter": "b", "row": 0, "col": 7, "colspan": 5},
        {"letter": "c", "row": 1, "col": 0, "colspan": 5},
        {"letter": "d", "row": 1, "col": 5, "colspan": 7}]}
    main4 = {"tag": "fig4", "ncol": 12, "row_heights_mm": [60, 52], "panels": [
        {"letter": "a", "row": 0, "col": 0, "colspan": 12},
        {"letter": "b", "row": 1, "col": 0, "colspan": 12}]}
    s6 = {"tag": "figS6", "ncol": 12, "row_heights_mm": [150], "panels": [
        {"letter": "a", "row": 0, "col": 0, "colspan": 12}]}
    s5 = {"tag": "figS5", "ncol": 12, "row_heights_mm": [52, 46], "panels": [
        {"letter": "a", "row": 0, "col": 0, "colspan": 7},
        {"letter": "b", "row": 0, "col": 7, "colspan": 5},
        {"letter": "c", "row": 1, "col": 0, "colspan": 5},
        {"letter": "d", "row": 1, "col": 5, "colspan": 7}]}
    s1 = {"tag": "figS1", "ncol": 12, "row_heights_mm": [52, 44], "panels": [
        {"letter": "a", "row": 0, "col": 0, "colspan": 12},
        {"letter": "b", "row": 1, "col": 0, "colspan": 6},
        {"letter": "c", "row": 1, "col": 6, "colspan": 6}]}
    s2 = {"tag": "figS2", "ncol": 12, "row_heights_mm": [64], "panels": [
        {"letter": "a", "row": 0, "col": 0, "colspan": 6},
        {"letter": "b", "row": 0, "col": 6, "colspan": 6}]}
    # The three metadata embeddings moved here out of Figure 1, where study and
    # assay duplicated a column header, a glyph row and a ring band, and the
    # sample panel duplicated Figure S6. Integration detail is their subject.
    s3 = {"tag": "figS3", "ncol": 12, "row_heights_mm": [52, 52], "panels": [
        {"letter": "a", "row": 0, "col": 0, "colspan": 6},
        {"letter": "b", "row": 0, "col": 6, "colspan": 6},
        {"letter": "c", "row": 1, "col": 0, "colspan": 12}]}
    s4 = {"tag": "figS4", "ncol": 12, "row_heights_mm": [46], "panels": [
        {"letter": "a", "row": 0, "col": 0, "colspan": 12}]}
    return main1, main2, main3, main4, main4b, main5, main6, s1, s2, s3, s4, s5, s6


def main() -> int:
    cfg = config.load()
    atlas = config.ROOT / config.get(cfg, "paths.atlas_out")
    T = config.ROOT / "results" / "tables"
    # One home per kind of file: composed figures are deliverables and live in
    # results/pub_figs under their submission names; per-panel PNGs are
    # regenerable intermediates and live under data/interim.
    out = config.ROOT / "results" / "pub_figs"
    panels = config.ROOT / "data" / "interim" / "figure_panels"
    out.mkdir(parents=True, exist_ok=True); panels.mkdir(parents=True, exist_ok=True)

    with provenance.start_run("06_figures", cfg, [atlas, T / "v5_embedding_benchmark.csv"]) as run:
        _style()
        A = ad.read_h5ad(atlas, backed="r")
        cols = ["study", "sample", "patient", "cell_type", "material_type", "treatment_timing",
                "assay", "cnv_score", "n_genes_by_counts", "total_counts", "pct_counts_mt",
                "pct_counts_spikein", "malignant_frozen", "malignant_evidence",
                "cnv_corr", "cnv_signal", "cnv_chr1q", "cnv_chr8q", "cnv_chr12q",
                "cnv_chr16q", "cnv_chr9p", "leiden"]
        d = A.obs[[c for c in cols if c in A.obs.columns]].copy()
        d["umap1"] = A.obsm["X_umap"][:, 0]; d["umap2"] = A.obsm["X_umap"][:, 1]
        A.file.close()

        # The marker dot plot is derived from the released object at render time,
        # not read from a stored table. The stored v5_marker_dotplot.csv predated
        # the annotation revision: it carried twelve classes, so panel D silently
        # omitted NK, mast and plasmacytoid dendritic cells while panel C showed
        # them as sectors, and its "DC" row was the mast-cell profile
        # (FCER1A 41%, CD1C 0.4%) under a dendritic-cell label.
        dot = marker_dotplot(atlas)
        bench = pd.read_csv(T / "v5_embedding_benchmark.csv", index_col=0)
        mix = pd.read_csv(T / "v5_mixing_by_celltype.csv")
        knn = json.loads((T / "v5_knn_mixing.json").read_text(encoding="utf-8"))
        auc = pd.read_csv(T / "v5_cnv_auc_per_sample.csv")
        att = pd.read_csv(T / "qc_attrition.csv")
        prof = pd.read_csv(T / "scpca_library_profile.csv")
        cov = pd.read_csv(config.ROOT / config.get(cfg, "paths.covariates"))
        scib = pd.read_csv(T / "v5_scib_results.csv", index_col=0)
        dep = pd.read_csv(T / "v5_depositor_labels.csv", index_col=0).reindex(d.index)
        lab = dep["depositor_label"].astype(str)
        conf = pd.crosstab(d.cell_type[lab != "unlabelled"], lab[lab != "unlabelled"])
        hist = {k: pd.read_csv(T / f"v5_elbo_{k}.csv", index_col=0)
                for k in ("sample_train", "sample_val", "study_train", "study_val")}
        coh = (d.groupby("study", observed=True)
                 .agg(cells=("sample", "size"), samples=("sample", "nunique"),
                      patients=("patient", "nunique"))
                 .join(d.groupby("study", observed=True)["assay"]
                        .agg(lambda x: ";".join(sorted(set(x.astype(str)))))
                        .rename("assays"))
                 .sort_values("cells", ascending=False))
        totals = {"studies": d["study"].nunique(), "samples": d["sample"].nunique(),
                  "patients": d["patient"].nunique(), "cells": len(d)}
        cnv6 = T / "v6_cnv_per_sample.csv"
        evid = T / "v6_malignant_evidence.csv"
        occ = T / "v6_metaprogram_occupancy.csv"
        falsi = T / "v6_falsification.csv"
        byev = T / "v6_program_by_evidence.csv"
        imm_f = config.ROOT / "data" / "interim" / "immune_subset.h5ad"
        IM = ad.read_h5ad(imm_f)
        imm_mask = d.index.isin(IM.obs_names)
        xy_all = d[["umap1", "umap2"]].to_numpy()
        floors = config.get(cfg, "qc.gene_floor_by_assay")
        caps = config.get(cfg, "qc.mito_cap_pct_by_assay")

        o1, o2, o3, o4, o4b, o5, o6, s1, s2, s3, s4, s5, s6 = build_outlines()
        # match the single UMAP panel's height to the data aspect so the figure
        # is the plot and nothing else
        xr = d["umap1"].max() - d["umap1"].min()
        yr = d["umap2"].max() - d["umap2"].min()
        s6["row_heights_mm"] = [round(WIDTH_MM * float(yr / xr), 1)]
        jobs = [
            (o1, {"a": (F.panel_cohort_schematic, dict(g=coh, d=d, totals=totals)),
                  "b": (F.panel_circular_atlas, dict(d=d, aspect_wh=180.0 / 86.0)),
                  "c": (F.panel_marker_dotplot, dict(dot=dot))},
             "Figure_1_atlas_overview.png"),
            (o2, {"a": (F.panel_scib_scatter, dict(bench=bench)),
                  "b": (F.panel_knn_vs_expected, dict(mix=mix)),
                  "c": (F.panel_knn_by_study, dict(knn=knn)),
                  "d": (F.panel_depositor_confusion, dict(conf=conf))}, "Figure_2_integration_quality.png"),
            (o3, {"a": (F.panel_cnv_statistic_gain, dict(old_csv=T / "v5_cnv_auc_per_sample.csv",
                                                         new_csv=cnv6)),
                  "b": (F.panel_arm_scores, dict(csv=T / "v6_cnv_arm_auc.csv")),
                  "c": (F.panel_evidence_bar, dict(csv=evid))},
             "Figure_3_malignant_compartment.png"),
            (o5, {"a": (F.panel_mp_occupancy, dict(csv=occ)),
                  "b": (F.panel_confounding_heatmap, dict(csv=falsi)),
                  "c": (F.panel_mp24_origin, dict(csv=byev))},
             "Figure_5_recurrence_and_confounding.png"),
            (o4b, {"a": (F.panel_caf_specificity, dict(csv=T / "f4_specificity_cells.csv")),
                   "b": (F.panel_caf_per_sample, dict(csv=T / "f4_split_genome_specificity.csv",
                                                      csv_ds=T / "f4_split_genome_specificity_ds1500.csv",
                                                      prev_csv=T / "f4_prevalence.csv")),
                   "c": (F.panel_caf_expression, dict(csv=T / "f4_group_expression.csv")),
                   "d": (F.panel_caf_prevalence, dict(csv=T / "f4_prevalence.csv"))},
             "Figure_S7_caf_like_state.png"),
            (o4, {"a": (F.panel_immune_umap, dict(adata=IM, xy_all=xy_all, mask=imm_mask)),
                  "b": (F.panel_immune_markers, dict(adata=IM))},
             "Figure_4_immune_microenvironment.png"),
            (o6, {"a": (F.panel_target_landscape, dict(csv=T / "f6_target_landscape.csv")),
                   "b": (F.panel_target_margin, dict(csv=T / "f6_target_margins.csv")),
                   "c": (F.panel_loso_accuracy, dict(csv=T / "f6_loso_mapping.csv")),
                   "d": (F.panel_bulk_survival, dict(csv=T / "f6_bulk_survival.csv"))},
             "Figure_6_atlas_utility.png"),
            (s5, {"a": (F.panel_cnv_by_celltype, dict(d=d)), "b": (F.panel_cnv_auc, dict(auc=auc)),
                  "c": (F.panel_cnv_vs_depth, dict(d=d)),
                  "d": (F.panel_boundary_disagreement, dict(conf=conf))},
             "Figure_S5_copy_number_first_attempt.png"),
            (s6, {"a": (F.panel_umap_by_sample, dict(d=d))},
             "Figure_S6_umap_by_sample.png"),
            (s1, {"a": (F.panel_qc_attrition, dict(att=att)),
                  "b": (F.panel_genes_by_assay, dict(d=d, floors=floors)),
                  "c": (F.panel_mito_by_assay, dict(d=d, caps=caps))}, "Figure_S1_quality_control.png"),
            (s2, {"a": (F.panel_library_gate, dict(prof=prof)),
                  "b": (F.panel_spikein_by_plate, dict(d=d))}, "Figure_S2_deposit_specific_qc.png"),
            (s3, {"a": (F.panel_training_curves, dict(hist=hist)),
                  "b": (F.panel_scib_heatmap, dict(scib=scib)),
                  "c": (F.panel_umap_triple, dict(d=d, horizontal=True))},
             "Figure_S3_integration_detail.png"),
            (s4, {"a": (F.panel_covariate_coverage, dict(cov=cov))}, "Figure_S4_covariate_coverage.png"),
        ]
        made = []
        for outline, spec, fname in jobs:
            paths = {L: render(outline, L, fn, panels, **kw) for L, (fn, kw) in spec.items()}
            p = compose(outline, paths, out / fname)
            made.append(p)
            print(f"{p.name}  {p.stat().st_size/1e6:.2f} MB", flush=True)
        run.add(figures=[p.name for p in made])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
