"""Panel renderers for the EWSAtlas figure set.

Every function draws into a supplied axis and sets **no title** - captions carry
the message, per the project's figure convention. That raises the bar on
in-figure identifiability (figure-style §2.1): each panel must be readable from
its axis labels, direct labels and legend alone.

Shared vocabulary (colour threading, §4.1) lives at module scope so a cell type
keeps the same colour in every figure it appears in.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

CT_ORDER = ["Ewing_sarcoma", "Ewing_proliferating", "T_cell", "NK", "B_cell", "Plasma",
            "Monocyte_Macrophage", "DC", "pDC", "Mast", "Platelet", "Endothelial",
            "Fibroblast_CAF", "Urothelial", "Unassigned"]
PRETTY = {"Ewing_sarcoma": "Ewing sarcoma", "Ewing_proliferating": "Ewing, proliferating",
          "T_cell": "T cell", "B_cell": "B cell", "Plasma": "Plasma cell",
          "Monocyte_Macrophage": "Monocyte / macrophage", "DC": "Dendritic cell",
          "NK": "NK cell", "pDC": "Plasmacytoid DC", "Mast": "Mast cell",
          "Platelet": "Platelet", "Endothelial": "Endothelial",
          "Fibroblast_CAF": "Fibroblast / CAF", "Urothelial": "Urothelial",
          "Unassigned": "Unassigned"}
CT_COLOR = {"Ewing_sarcoma": "#5B2C83", "Ewing_proliferating": "#9A6FC0", "T_cell": "#1F6FB4",
            "B_cell": "#6BAED6", "Plasma": "#9ECAE1", "Monocyte_Macrophage": "#117733",
            "DC": "#66A61E", "NK": "#4393C3", "pDC": "#8DD3C7", "Mast": "#E7A33E", "Platelet": "#B2DF8A", "Endothelial": "#D95F02",
            "Fibroblast_CAF": "#A6761D", "Urothelial": "#E7298A", "Unassigned": "#BDBDBD"}
#: only spatially coherent classes are direct-labelled on an embedding; a
#: dispersed class labelled at its median names whatever it is scattered through
DIRECT = {"Ewing_sarcoma", "Ewing_proliferating", "T_cell", "B_cell",
          "Monocyte_Macrophage", "Endothelial", "Fibroblast_CAF", "Urothelial"}
STUDY_ORDER = ["Waltner_2025", "Mao_2024", "Aynaud_2020", "Cillo_2022",
               "Goodspeed_2025", "Visser_2023"]
SCOLOR = dict(zip(STUDY_ORDER, ["#332288", "#88CCEE", "#44AA99", "#DDCC77", "#CC6677", "#AA4499"]))
SPRETTY = {s: s.replace("_", " ") for s in STUDY_ORDER}
MAT = ["patient_tumor", "O-PDX", "CTC", "PBMC", "normal_tissue"]
MATP = {"patient_tumor": "Patient tumour", "O-PDX": "Orthotopic PDX",
        "CTC": "Circulating tumour cell", "PBMC": "Blood (PBMC)", "normal_tissue": "Normal tissue"}
MATC = dict(zip(MAT, ["#5B2C83", "#9A6FC0", "#CC6677", "#88CCEE", "#44AA99"]))
ASSAYP = {"10x_scRNA": "10x, whole cell", "10x_snRNA": "10x, nuclei", "CELseq2_scRNA": "CEL-seq2, plate"}
ASSAYC = {"10x_scRNA": "#1F6FB4", "10x_snRNA": "#882255", "CELseq2_scRNA": "#DDAA33"}
EPRETTY = {"scVI_batch_sample": "scVI (batch = sample)", "scVI_batch_study": "scVI (batch = study)",
           "Harmony_sample": "Harmony (batch = sample)", "Unintegrated_PCA": "Unintegrated PCA"}
FOCAL = "scVI_batch_sample"
REF_TYPES = ["T_cell", "B_cell", "Plasma", "Monocyte_Macrophage", "DC", "Platelet"]
MAL_TYPES = ["Ewing_sarcoma", "Ewing_proliferating"]
GREY = "#BDBDBD"

#: marker panel order for the dot plot, grouped by the compartment it identifies
MARKER_GROUPS = [
    ("Ewing", ["CD99", "NKX2-2", "NR0B1", "FLI1"]),
    ("Cycling", ["MKI67", "TOP2A"]),
    ("Leukocyte", ["PTPRC"]),
    ("T / NK", ["CD3D", "CD3E", "IL7R", "NKG7", "GNLY"]),
    ("B / plasma", ["MS4A1", "CD79A", "MZB1", "JCHAIN"]),
    ("Myeloid", ["LYZ", "CD14", "C1QA", "AIF1"]),
    # the annotation revision split mast, conventional and plasmacytoid
    # dendritic cells apart, so the panel needs a marker for each
    ("DC", ["CD1C", "FCER1A"]),
    ("pDC", ["LILRA4", "IRF7"]),
    ("Mast", ["TPSAB1", "TPSB2"]),
    ("Platelet", ["PPBP", "PF4"]),
    ("Endothelial", ["PECAM1", "CDH5"]),
    ("Fibroblast", ["COL1A1", "DCN", "LUM"]),
    ("Urothelial", ["KRT18", "KRT8", "UPK1B"]),
]


def _despine(ax, keep=("left", "bottom")):
    for sp, v in ax.spines.items():
        v.set_visible(sp in keep)


def _tier(name, meta):
    return f"{name}\n{meta}"


def _embedding_axes(ax, xlab="UMAP 1", ylab="UMAP 2"):
    ax.set_xticks([]); ax.set_yticks([])
    _despine(ax, keep=())
    ax.margins(0.04)
    x0, y0 = ax.get_xlim()[0], ax.get_ylim()[0]
    dx = 0.15 * (ax.get_xlim()[1] - x0); dy = 0.15 * (ax.get_ylim()[1] - y0)
    for xy, xytext in [((x0 + dx, y0), (x0, y0)), ((x0, y0 + dy), (x0, y0))]:
        ax.annotate("", xy=xy, xytext=xytext, arrowprops=dict(arrowstyle="-|>", color="0.35", lw=0.7))
    ax.text(x0 + dx * 1.06, y0, xlab, fontsize=5.5, color="0.35", va="center", ha="left")
    ax.text(x0, y0 + dy * 1.06, ylab, fontsize=5.5, color="0.35", ha="center", va="bottom", rotation=90)


# ---------------------------------------------------------------- Figure 1 ---

def panel_umap_celltype(ax, d):
    rng = np.random.default_rng(0)
    dd = d.iloc[rng.permutation(len(d))]
    for lvl in CT_ORDER:
        s = dd[dd.cell_type == lvl]
        ax.scatter(s.umap1, s.umap2, s=0.5, lw=0, c=CT_COLOR[lvl], rasterized=True, alpha=0.5)
    _embedding_axes(ax)
    for lvl in (c for c in CT_ORDER if c in DIRECT):
        s = d[d.cell_type == lvl]
        ax.text(s.umap1.median(), s.umap2.median(), PRETTY[lvl], fontsize=5.6, ha="center",
                va="center", bbox=dict(fc="white", ec="none", alpha=0.7, pad=0.6))
    rest = [c for c in CT_ORDER if c not in DIRECT]
    ax.legend([plt.Line2D([], [], marker="o", ls="", ms=3, mfc=CT_COLOR[c], mec="none") for c in rest],
              [f"{PRETTY[c]} ({int((d.cell_type == c).sum()):,})" for c in rest],
              frameon=False, fontsize=5.4, loc="lower right", handletextpad=0.25,
              labelspacing=0.2, borderpad=0.05)


def panel_umap_study(ax, d):
    rng = np.random.default_rng(0)
    dd = d.iloc[rng.permutation(len(d))]
    for lvl in STUDY_ORDER:
        s = dd[dd.study == lvl]
        ax.scatter(s.umap1, s.umap2, s=0.5, lw=0, c=SCOLOR[lvl], rasterized=True, alpha=0.5)
    _embedding_axes(ax)
    ax.legend([plt.Line2D([], [], marker="o", ls="", ms=3, mfc=SCOLOR[s], mec="none") for s in STUDY_ORDER],
              [f"{SPRETTY[s]} ({int((d.study == s).sum()):,})" for s in STUDY_ORDER],
              frameon=False, fontsize=5.4, loc="lower right", handletextpad=0.25,
              labelspacing=0.2, borderpad=0.05)


def panel_marker_dotplot(ax, dot):
    genes, boundaries, labels = [], [], []
    for name, gs in MARKER_GROUPS:
        present = [g for g in gs if g in set(dot.gene)]
        if not present:
            continue
        labels.append((name, len(genes) + len(present) / 2 - 0.5))
        genes += present
        boundaries.append(len(genes) - 0.5)
    cts = [c for c in CT_ORDER if c in set(dot.cell_type)][::-1]
    piv_f = dot.pivot(index="cell_type", columns="gene", values="frac_expressing").reindex(cts)[genes]
    piv_m = dot.pivot(index="cell_type", columns="gene", values="mean_lognorm").reindex(cts)[genes]
    X, Y = np.meshgrid(np.arange(len(genes)), np.arange(len(cts)))
    sizes = (piv_f.to_numpy() * 44.0) + 0.6
    vmax = float(np.nanpercentile(piv_m.to_numpy(), 98))
    sc = ax.scatter(X.ravel(), Y.ravel(), s=sizes.ravel(), c=piv_m.to_numpy().ravel(),
                    cmap="viridis", vmin=0, vmax=vmax, lw=0)
    for b in boundaries[:-1]:
        ax.axvline(b, color="0.85", lw=0.5, zorder=0)
    ax.set_xticks(np.arange(len(genes)))
    ax.set_xticklabels(genes, rotation=90, style="italic")
    ax.set_yticks(np.arange(len(cts)))
    ax.set_yticklabels([PRETTY[c] for c in cts])
    # group headers stagger across two tiers: adjacent short/long names collided
    for i, (name, x) in enumerate(labels):
        ax.text(x, len(cts) - 0.3 + (0.95 if i % 2 else 0.0), name, ha="center", va="bottom",
                fontsize=5.4, color="0.3")
    ax.set_ylim(-0.9, len(cts) + 1.3); ax.set_xlim(-0.8, len(genes) - 0.2)
    _despine(ax, keep=())
    ax.tick_params(length=0)
    # colorbar and size key sit in reserved insets so they cannot overlap
    cax = ax.inset_axes([1.012, 0.52, 0.012, 0.44])
    cb = ax.figure.colorbar(sc, cax=cax)
    cb.set_label("Mean expression\n(log-normalized)", fontsize=5.4, linespacing=1.2)
    cb.ax.tick_params(labelsize=5.0, length=1.2)
    cb.outline.set_visible(False)
    for f, lab in [(0.25, "25%"), (0.75, "75%")]:
        ax.scatter([], [], s=f * 44 + 0.6, c="0.45", lw=0, label=lab)
    ax.legend(frameon=False, fontsize=5.2, loc="upper left", bbox_to_anchor=(1.005, 0.40),
              title="Cells expressing", title_fontsize=5.4, handletextpad=0.3,
              labelspacing=0.55, borderpad=0.0)


def panel_material_by_study(ax, d):
    ypos = np.arange(len(STUDY_ORDER))
    mt = pd.crosstab(d.study, d.material_type).reindex(STUDY_ORDER).reindex(columns=MAT).fillna(0)[::-1]
    left = np.zeros(len(mt))
    for c in MAT:
        ax.barh(ypos, mt[c].to_numpy(), left=left, color=MATC[c], height=0.66, lw=0)
        left += mt[c].to_numpy()
    ax.set_yticks(ypos)
    ax.set_yticklabels([_tier(SPRETTY[s], f"{d[d.study == s]['sample'].nunique()} samples · "
                                        f"{d[d.study == s].patient.nunique()} patients")
                        for s in STUDY_ORDER[::-1]])
    ax.set_xlabel("Cells"); ax.set_xlim(0, 86000)
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, p: f"{v/1000:.0f}k" if v else "0"))
    ax.legend([plt.Rectangle((0, 0), 1, 1, fc=MATC[c], ec="none") for c in MAT],
              [MATP[c] for c in MAT], frameon=False, fontsize=5.4, loc="lower right",
              handlelength=0.8, handletextpad=0.35, labelspacing=0.2, borderpad=0.05)
    _despine(ax); ax.margins(y=0.04)


def panel_celltype_by_study(ax, d):
    ypos = np.arange(len(STUDY_ORDER))
    comp = pd.crosstab(d.study, d.cell_type, normalize="index").reindex(STUDY_ORDER)
    comp = comp.reindex(columns=CT_ORDER).fillna(0)[::-1]
    left = np.zeros(len(comp))
    for c in CT_ORDER:
        ax.barh(ypos, comp[c].to_numpy(), left=left, color=CT_COLOR[c], height=0.66, lw=0)
        left += comp[c].to_numpy()
    ax.set_yticks(ypos); ax.set_yticklabels([SPRETTY[s] for s in STUDY_ORDER[::-1]])
    ax.set_xlabel("Fraction of cells  (cell-type colours as in panel a)")
    ax.set_xlim(0, 1)
    _despine(ax); ax.margins(y=0.04)


# ---------------------------------------------------------------- Figure 2 ---

def panel_scib_scatter(ax, bench):
    for e, row in bench.iterrows():
        focal = e == FOCAL
        ax.scatter(row["Batch correction"], row["Bio conservation"], s=66 if focal else 40,
                   c=CT_COLOR["Ewing_sarcoma"] if focal else GREY,
                   ec="black" if focal else "0.45", lw=0.7 if focal else 0.4,
                   zorder=3 if focal else 2)
        ax.annotate(f"{EPRETTY[e]}\n{row['Total']:.3f}",
                    xy=(row["Batch correction"], row["Bio conservation"]),
                    xytext=(row["Batch correction"] + 0.013, row["Bio conservation"]),
                    ha="left", va="center", fontsize=5.6,
                    color="black" if focal else "0.3", linespacing=1.3)
    ax.set_xlabel("Batch correction (scIB aggregate)")
    ax.set_ylabel("Biological conservation\n(scIB aggregate)")
    ax.set_xlim(0.26, 0.82); ax.set_ylim(0.47, 0.64)
    ax.text(0.98, 0.03, "higher = better on both axes", transform=ax.transAxes,
            ha="right", fontsize=5.5, color="0.35")
    _despine(ax)


def panel_knn_vs_expected(ax, mix):
    m = mix.sort_values("cells", ascending=True)
    y = np.arange(len(m))
    for yi, (o, e) in enumerate(zip(m.observed, m.expected_given_composition)):
        ax.plot([e, o], [yi, yi], color="0.78", lw=1.2, zorder=1, solid_capstyle="round")
    ax.scatter(m.expected_given_composition, y, s=22, c="white", ec="0.35", lw=0.8, zorder=3)
    ax.scatter(m.observed, y, s=26, c=CT_COLOR["Ewing_sarcoma"], ec="none", zorder=4)
    top = len(m) - 1
    ax.annotate("observed", xy=(m.observed.iloc[top], top),
                xytext=(m.observed.iloc[top], top + 0.5), ha="center",
                fontsize=5.6, color=CT_COLOR["Ewing_sarcoma"])
    ax.annotate("expected given\ncell-type composition",
                xy=(m.expected_given_composition.iloc[top], top),
                xytext=(m.expected_given_composition.iloc[top] - 0.02, top + 0.55),
                ha="right", fontsize=5.6, color="0.35", linespacing=1.25)
    ax.set_yticks(y)
    ax.set_yticklabels([PRETTY.get(c, c) for c in m.cell_type])
    ax.set_xlabel("Share of 15 nearest neighbours from the same study")
    ax.set_xlim(0.1, 1.04); ax.set_ylim(-0.7, len(m) + 0.8)
    _despine(ax)


def panel_knn_by_study(ax, knn):
    y = np.arange(len(STUDY_ORDER))
    for e in ["Unintegrated_PCA", "Harmony_sample", "scVI_batch_study", FOCAL]:
        vals = [knn[e]["by_study"][s] for s in STUDY_ORDER]
        focal = e == FOCAL
        ax.scatter(vals, y, s=34 if focal else 22,
                   c=CT_COLOR["Ewing_sarcoma"] if focal else "none",
                   ec="black" if focal else "0.5", lw=0.7 if focal else 0.7,
                   marker="o" if focal else ["s", "^", "D"][
                       ["Unintegrated_PCA", "Harmony_sample", "scVI_batch_study"].index(e)],
                   zorder=4 if focal else 2, label=EPRETTY[e])
    ax.set_yticks(y); ax.set_yticklabels([SPRETTY[s] for s in STUDY_ORDER])
    ax.set_xlabel("Share of 15 nearest neighbours from the same study")
    ax.set_xlim(0.3, 1.02)
    ax.legend(frameon=False, fontsize=5.4, loc="lower left", handletextpad=0.25,
              labelspacing=0.25, borderpad=0.05)
    _despine(ax); ax.margins(y=0.08)


def panel_depositor_confusion(ax, conf):
    """Row-normalized agreement between our labels and the depositors'."""
    mat = conf.div(conf.sum(axis=1), axis=0)
    # order rows so each one's heaviest column walks down the matrix: the
    # correspondence is the point of the panel and alphabetical order hides it
    mat = mat.loc[sorted(mat.index, key=lambda r: (list(mat.columns).index(mat.loc[r].idxmax()),
                                                   -mat.loc[r].max()))]
    dep = {"B_plasma": "B / plasma", "Blood_other": "Blood, other", "T_NK": "T / NK",
           "pDC": "pDC", "Myeloid": "Myeloid", "Stromal": "Stromal",
           "Endothelial": "Endothelial", "Epithelial": "Epithelial", "Malignant": "Malignant"}
    im = ax.imshow(mat.to_numpy(), cmap="Purples", vmin=0, vmax=1, aspect="auto")
    ax.set_xticks(np.arange(mat.shape[1]))
    ax.set_xticklabels([dep.get(c, str(c).replace("_", " ")) for c in mat.columns], rotation=90)
    ax.set_yticks(np.arange(mat.shape[0]))
    ax.set_yticklabels([PRETTY.get(c, c) for c in mat.index])
    ax.set_xlabel("Depositor annotation")
    ax.set_ylabel("This atlas")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.iat[i, j]
            if v >= 0.02:
                ax.text(j, i, f"{100*v:.0f}", ha="center", va="center", fontsize=4.8,
                        color="white" if v > 0.55 else "0.25")
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("% of row", fontsize=5.6); cb.ax.tick_params(labelsize=5.2, length=1.5)
    cb.outline.set_visible(False)
    _despine(ax, keep=()); ax.tick_params(length=0)


# ---------------------------------------------------------------- Figure 3 ---

def panel_cnv_by_celltype(ax, d):
    cn = d[d.cnv_score.notna()]
    med = cn.groupby("cell_type", observed=True)["cnv_score"].median().sort_values()
    for i, ct in enumerate(med.index):
        v = cn.loc[cn.cell_type == ct, "cnv_score"].to_numpy()
        col = (CT_COLOR["Ewing_sarcoma"] if ct in MAL_TYPES
               else CT_COLOR["T_cell"] if ct in REF_TYPES else CT_COLOR["Fibroblast_CAF"])
        parts = ax.violinplot([v], positions=[i], widths=0.78, showextrema=False, vert=False)
        for b in parts["bodies"]:
            b.set_facecolor(col); b.set_alpha(0.55); b.set_edgecolor("none")
        ax.plot([np.median(v)], [i], marker="|", ms=7, mew=1.2, color="black", zorder=4)
    ax.set_yticks(range(len(med)))
    ax.set_yticklabels([PRETTY.get(c, c) for c in med.index])
    ax.set_xlabel("Inferred copy-number score (mean |CNV|)")
    ax.set_xlim(0, 0.045)
    ax.legend([plt.Rectangle((0, 0), 1, 1, fc=CT_COLOR["Ewing_sarcoma"], ec="none", alpha=0.55),
               plt.Rectangle((0, 0), 1, 1, fc=CT_COLOR["T_cell"], ec="none", alpha=0.55),
               plt.Rectangle((0, 0), 1, 1, fc=CT_COLOR["Fibroblast_CAF"], ec="none", alpha=0.55)],
              ["Malignant by markers", "Diploid reference (immune)", "Other"],
              frameon=False, fontsize=5.4, loc="lower right", handlelength=0.8,
              handletextpad=0.35, labelspacing=0.2, borderpad=0.05)
    _despine(ax)


def panel_cnv_auc(ax, auc):
    a = auc.sort_values("auc").reset_index(drop=True)
    ax.axvline(0.5, color="0.6", lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.scatter(a.auc, np.arange(len(a)), s=16, c=CT_COLOR["Ewing_sarcoma"], ec="none", zorder=3)
    ax.set_yticks([])
    ax.set_xlabel("Malignant-vs-reference AUC, per sample")
    ax.set_ylabel(f"{len(a)} samples")
    ax.set_xlim(0.28, 1.0)
    ax.annotate("no discrimination", xy=(0.5, len(a) * 0.96), xytext=(0.507, len(a) * 0.96),
                fontsize=5.6, color="0.35", ha="left", va="center")
    ax.annotate(f"median {a.auc.median():.3f}\n{int((a.auc > 0.7).sum())} of {len(a)} above 0.7",
                xy=(0.985, 0.06), xycoords="axes fraction", ha="right", fontsize=5.6,
                color="0.25", linespacing=1.3)
    _despine(ax); ax.margins(y=0.03)


def panel_cnv_vs_depth(ax, d):
    """How the score moves with sequencing depth, against how much it moves
    between malignant and reference cells.

    The global rank correlation with detected genes is only -0.08, so depth is
    not a simple monotonic driver. What the binned medians show instead is that
    the score swings further across depth bins than it ever separates malignant
    from reference cells, and that the two groups are indistinguishable in the
    shallowest bin - which is why the per-sample AUC in panel b sits at 0.54.
    """
    cn = d[d.cnv_score.notna()].copy()
    cn["grp"] = np.where(cn.cell_type.isin(MAL_TYPES), "malignant",
                         np.where(cn.cell_type.isin(REF_TYPES), "reference", "other"))
    bins = np.unique(np.quantile(cn.n_genes_by_counts, np.linspace(0, 1, 13)))
    cn["bin"] = pd.cut(cn.n_genes_by_counts, bins, include_lowest=True)
    curves = {}
    for grp, col, lab in [("reference", CT_COLOR["T_cell"], "Diploid reference cells"),
                          ("malignant", CT_COLOR["Ewing_sarcoma"], "Malignant cells")]:
        g = cn[cn.grp == grp].groupby("bin", observed=True)
        x, y, n = g.n_genes_by_counts.median(), g.cnv_score.median(), g.size()
        keep = (n >= 50).to_numpy()
        curves[grp] = y[keep]
        ax.plot(x[keep], y[keep], marker="o", ms=2.6, lw=1.1, color=col, label=lab)
    gap = (curves["malignant"] - curves["reference"]).median()
    allv = pd.concat(curves.values())
    swing = float(allv.max() - allv.min())
    ax.set_xscale("log")
    ticks = [100, 300, 1000, 3000, 10000]
    ax.set_xticks(ticks)
    ax.set_xticklabels(["100", "300", "1k", "3k", "10k"])
    ax.xaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    ax.set_xlabel("Detected genes per cell (median of bin)")
    ax.set_ylabel("Inferred copy-number score")
    ax.annotate(f"swing across depth {swing:.3f}\nmalignant − reference {gap:.3f}\n"
                f"rank correlation with depth −0.08",
                xy=(0.03, 0.06), xycoords="axes fraction", fontsize=5.6,
                va="bottom", color="0.25", linespacing=1.35)
    ax.legend(frameon=False, fontsize=5.4, loc="upper left", handletextpad=0.35,
              labelspacing=0.2, borderpad=0.05)
    _despine(ax)


def panel_boundary_disagreement(ax, conf):
    """Where our labels and the depositors' disagree, in cells."""
    pairs = [("Fibroblast_CAF", "Malignant"), ("Ewing_sarcoma", "Stromal"),
             ("T_cell", "Malignant"), ("Ewing_sarcoma", "T_NK"),
             ("Unassigned", "Epithelial"), ("Monocyte_Macrophage", "Malignant")]
    vals, labs, cols = [], [], []
    for ours, theirs in pairs:
        if ours in conf.index and theirs in conf.columns:
            vals.append(int(conf.at[ours, theirs]))
            labs.append(f"{PRETTY.get(ours, ours)}  ·  {theirs.replace('_', '/')}")
            cols.append(CT_COLOR.get(ours, GREY))
    o = np.argsort(vals)
    y = np.arange(len(vals))
    ax.barh(y, np.array(vals)[o], color=np.array(cols)[o], height=0.6, lw=0)
    for yi, v in zip(y, np.array(vals)[o]):
        ax.text(v + max(vals) * 0.02, yi, f"{v:,}", va="center", fontsize=5.6, color="0.25")
    ax.set_yticks(y); ax.set_yticklabels(np.array(labs)[o])
    ax.set_xlabel("Cells labelled differently (this atlas  ·  depositor)")
    ax.set_xlim(0, max(vals) * 1.22)
    _despine(ax); ax.margins(y=0.08)


# ------------------------------------------------------- Supplementary ------

FILTERS = [("drop_failed_library", "Failed library (whole sample)", "#882255"),
           ("drop_below_gene_floor", "Below detected-gene floor", "#CC6677"),
           ("drop_above_mito_cap", "Above mitochondrial cap", "#DDAA33"),
           ("drop_mad_outlier", "Per-sample MAD outlier", "#88CCEE"),
           ("drop_doublet", "Doublet", "#44AA99")]


def panel_qc_attrition(ax, att):
    """Where every discarded barcode went, per study."""
    a = att.set_index("study").reindex(STUDY_ORDER)
    y = np.arange(len(a))[::-1]
    left = np.zeros(len(a))
    for col, lab, colr in FILTERS:
        v = a[col].fillna(0).to_numpy()
        ax.barh(y, v, left=left, color=colr, height=0.62, lw=0, label=lab)
        left += v
    ax.barh(y, a["kept"].to_numpy(), left=left, color="0.82", height=0.62, lw=0, label="Retained")
    ax.set_yticks(y)
    ax.set_yticklabels([_tier(SPRETTY[s], f"{a.loc[s, 'pct_kept']:.1f}% retained") for s in STUDY_ORDER])
    ax.set_xlabel("Barcodes")
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, p: f"{v/1000:.0f}k" if v else "0"))
    ax.legend(frameon=False, fontsize=5.4, loc="lower right", handlelength=0.8,
              handletextpad=0.35, labelspacing=0.2, borderpad=0.05)
    _despine(ax); ax.margins(y=0.05)


def _assay_hist(ax, d, col, xlabel, logx=False, cap=None):
    for assay in ["10x_scRNA", "CELseq2_scRNA", "10x_snRNA"]:
        v = d.loc[d.assay == assay, col].dropna().to_numpy()
        if not len(v):
            continue
        if logx:
            v = v[v > 0]
            bins = np.logspace(np.log10(max(v.min(), 1)), np.log10(v.max()), 60)
        else:
            bins = np.linspace(0, cap if cap else v.max(), 60)
        ax.hist(v, bins=bins, density=True, histtype="step", lw=1.1,
                color=ASSAYC[assay], label=ASSAYP[assay])
    if logx:
        ax.set_xscale("log")
        ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(
            lambda v, p: f"{v/1000:.0f}k" if v >= 1000 else f"{v:.0f}"))
    ax.set_xlabel(xlabel); ax.set_ylabel("Density")
    ax.set_yticks([])
    _despine(ax)


def panel_genes_by_assay(ax, d, floors):
    _assay_hist(ax, d, "n_genes_by_counts", "Detected genes per cell (retained cells)", logx=True)
    for assay, style in [("10x_scRNA", "-"), ("10x_snRNA", (0, (4, 2)))]:
        ax.axvline(floors[assay], color=ASSAYC[assay], lw=0.8, ls=style, zorder=0)
    # rotated and set just right of each line: stacked horizontal labels collided
    # with the curves and with each other
    ax.annotate("floor 100 (nuclei)", xy=(106, 0.55), xycoords=("data", "axes fraction"),
                fontsize=5.2, color=ASSAYC["10x_snRNA"], ha="left", va="center", rotation=90)
    ax.annotate("floor 200 (whole cell, plate)", xy=(212, 0.55), xycoords=("data", "axes fraction"),
                fontsize=5.2, color=ASSAYC["10x_scRNA"], ha="left", va="center", rotation=90)
    ax.legend(frameon=False, fontsize=5.4, loc="upper right", handletextpad=0.5,
              labelspacing=0.25, borderpad=0.05)


def panel_mito_by_assay(ax, d, caps):
    _assay_hist(ax, d, "pct_counts_mt", "Mitochondrial content (%)", cap=45)
    for assay in ["10x_snRNA", "10x_scRNA", "CELseq2_scRNA"]:
        ax.axvline(caps[assay], color=ASSAYC[assay], lw=0.8, ls=(0, (4, 2)), zorder=0)
    ax.annotate("caps: 10 / 30 / 40", xy=(0.97, 0.93), xycoords="axes fraction",
                ha="right", va="top", fontsize=5.5, color="0.3")


def panel_library_gate(ax, prof):
    """Per-library median detected genes against its assay floor."""
    p = prof.sort_values("median_genes")
    y = np.arange(len(p))
    fail = p.median_genes < 100
    ax.scatter(p.median_genes[~fail], y[~fail.to_numpy()], s=18,
               c=ASSAYC["10x_snRNA"], ec="none", zorder=3)
    ax.scatter(p.median_genes[fail], y[fail.to_numpy()], s=30, facecolors="none",
               edgecolors="#882255", lw=1.0, zorder=4)
    ax.axvline(100, color="0.45", lw=0.8, ls=(0, (4, 3)), zorder=1)
    ax.set_yticks(y); ax.set_yticklabels(p.scpca_library_id, fontsize=4.8)
    ax.set_xscale("log")
    ax.set_xlabel("Median detected genes per nucleus")
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, p_: f"{v:.0f}"))
    ax.annotate("assay floor", xy=(100, len(p) - 0.4), xytext=(112, len(p) - 0.4),
                fontsize=5.4, color="0.4", va="center")
    for _, r in p[fail].iterrows():
        ax.annotate("excluded whole", xy=(r.median_genes, list(p.scpca_library_id).index(r.scpca_library_id)),
                    xytext=(r.median_genes * 1.25, list(p.scpca_library_id).index(r.scpca_library_id)),
                    fontsize=5.4, color="#882255", va="center")
    _despine(ax); ax.margins(y=0.05)


def panel_spikein_by_plate(ax, d):
    """CEL-seq2 spike-in fraction per retained plate."""
    v = (d[d.assay == "CELseq2_scRNA"].groupby("sample", observed=True)["pct_counts_spikein"]
         .median().sort_values())
    y = np.arange(len(v))
    ax.scatter(v.to_numpy(), y, s=16, c=ASSAYC["CELseq2_scRNA"], ec="none")
    ax.set_yticks(y); ax.set_yticklabels(v.index, fontsize=4.8)
    ax.set_xlabel("ERCC spike-in share of counts (%), median per plate")
    ax.set_xlim(0, max(60, float(v.max()) * 1.15))
    ax.annotate("four plates at 66–96% spike-in\nwere excluded before this step",
                xy=(0.97, 0.06), xycoords="axes fraction", ha="right", fontsize=5.5,
                color="0.3", linespacing=1.3)
    _despine(ax); ax.margins(y=0.04)


def panel_training_curves(ax, hist):
    """scVI ELBO for both batch-key arms; validation is what stops training."""
    for arm, col, dy in [("sample", CT_COLOR["Ewing_sarcoma"], -34), ("study", GREY, 30)]:
        tr, va = hist[f"{arm}_train"], hist[f"{arm}_val"]
        ax.plot(tr.index, tr.iloc[:, 0], lw=0.8, color=col, alpha=0.55)
        ax.plot(va.index, va.iloc[:, 0], lw=1.3, color=col)
        # single line, offset above or below its own curve so the two never meet
        ax.annotate(f"batch = {arm}, early stop at epoch {len(va)}",
                    xy=(len(va) - 1, va.iloc[-1, 0]),
                    xytext=(len(va) - 4, va.iloc[-1, 0] + dy),
                    fontsize=5.6, color=col, ha="right", va="center")
    ax.set_xlabel("Training epoch"); ax.set_ylabel("ELBO")
    ax.set_ylim(1450, 1760)
    ax.annotate("thin = training, thick = validation", xy=(0.97, 0.95),
                xycoords="axes fraction", ha="right", va="top", fontsize=5.5, color="0.35")
    _despine(ax)


def panel_scib_heatmap(ax, scib):
    s = scib[scib.index != "Metric Type"].astype(float)
    cols = [c for c in s.columns if c not in ("Batch correction", "Bio conservation", "Total")]
    s = s.loc[[FOCAL, "scVI_batch_study", "Harmony_sample", "Unintegrated_PCA"], cols]
    im = ax.imshow(s.to_numpy(), cmap="Purples", vmin=0, vmax=1, aspect="auto")
    # cLISI saturates at 1.00 and iLISI floors at 0.00 for every embedding here,
    # so those columns rank nothing; say so on the axis rather than letting a
    # reader weigh them alongside the metrics that do separate the methods
    flat = [c for c in cols if float(s[c].max() - s[c].min()) < 0.005]  # identical at 2 dp
    ax.set_xticks(np.arange(len(cols)))
    ax.set_xticklabels([f"{c} (constant)" if c in flat else c for c in cols], rotation=90)
    for t, c in zip(ax.get_xticklabels(), cols):
        if c in flat:
            t.set_color("0.55")
    ax.set_yticks(np.arange(s.shape[0])); ax.set_yticklabels([EPRETTY[i] for i in s.index])
    for i in range(s.shape[0]):
        for j in range(s.shape[1]):
            v = s.iat[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center", fontsize=4.9,
                    color="white" if v > 0.55 else "0.25")
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.025, pad=0.015)
    cb.set_label("Metric value", fontsize=5.6); cb.ax.tick_params(labelsize=5.2, length=1.5)
    cb.outline.set_visible(False)
    _despine(ax, keep=()); ax.tick_params(length=0)


def panel_covariate_coverage(ax, cov):
    fields = ["patient", "fusion_partner", "sex", "age_years", "tumor_site", "disease_state",
              "treatment_timing", "material_type", "assay", "sort_strategy"]
    fields = [f for f in fields if f in cov.columns]
    mat = pd.DataFrame(index=STUDY_ORDER, columns=fields, dtype=float)
    for s in STUDY_ORDER:
        sub = cov[cov.study == s]
        for f in fields:
            mat.loc[s, f] = 100.0 * (sub[f].astype(str) != "not_reported").mean()
    im = ax.imshow(mat.to_numpy(), cmap="Purples", vmin=0, vmax=100, aspect="auto")
    ax.set_xticks(np.arange(len(fields)))
    nice = {"patient": "Patient", "fusion_partner": "Fusion partner", "sex": "Sex",
            "age_years": "Age (years)", "tumor_site": "Tumour site",
            "disease_state": "Disease state", "treatment_timing": "Treatment timing",
            "material_type": "Material type", "assay": "Assay",
            "sort_strategy": "Sort strategy"}
    ax.set_xticklabels([nice.get(f, f.replace("_", " ")) for f in fields], rotation=90)
    ax.set_yticks(np.arange(len(STUDY_ORDER)))
    ax.set_yticklabels([SPRETTY[s] for s in STUDY_ORDER])
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.iat[i, j]
            ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=5.0,
                    color="white" if v > 55 else "0.25")
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.025, pad=0.015)
    cb.set_label("% of samples reported", fontsize=5.6)
    cb.ax.tick_params(labelsize=5.2, length=1.5); cb.outline.set_visible(False)
    _despine(ax, keep=()); ax.tick_params(length=0)

def locator_inset(ax, xy_all, mask, loc=(0.72, 0.72, 0.26, 0.26), size=0.6):
    """Small 'where are these cells' map: full atlas in grey, subset highlighted.

    NBAtlas puts one of these on every zoom panel, and it answers the question a
    reader always has about a re-clustered subset - which part of the atlas is
    this? - without a sentence of caption.
    """
    ins = ax.inset_axes(loc)
    ins.scatter(xy_all[~mask, 0], xy_all[~mask, 1], s=size, c="0.85", linewidths=0, rasterized=True)
    ins.scatter(xy_all[mask, 0], xy_all[mask, 1], s=size, c="#4C2A85", linewidths=0, rasterized=True)
    ins.set_xticks([]); ins.set_yticks([])
    for sp in ins.spines.values():
        sp.set_linewidth(0.4); sp.set_color("0.6")
    return ins


def density_panel(ax, xy, mask, cmap, levels=8, size=0.5):
    """Embedding density for one group, with the other cells as grey context.

    Preferred over two overplotted colours whenever two clinical groups are
    compared on the same embedding: with 10^5 cells the colour drawn second
    simply hides the first, and a density surface does not.
    """
    from scipy.stats import gaussian_kde
    import numpy as np

    ax.scatter(xy[:, 0], xy[:, 1], s=size, c="0.88", linewidths=0, rasterized=True)
    pts = xy[mask]
    if len(pts) > 20000:
        rng = np.random.default_rng(0)
        pts = pts[rng.choice(len(pts), 20000, replace=False)]
    k = gaussian_kde(pts.T)
    xg = np.linspace(xy[:, 0].min(), xy[:, 0].max(), 120)
    yg = np.linspace(xy[:, 1].min(), xy[:, 1].max(), 120)
    X, Y = np.meshgrid(xg, yg)
    Z = k(np.vstack([X.ravel(), Y.ravel()])).reshape(X.shape)
    ax.contourf(X, Y, Z, levels=levels, cmap=cmap, alpha=0.85)
    return ax

# ---------------------------------------------------------------------------
# Figures for the reframed Resource paper (v6): the malignant compartment, what
# recurs across six studies, the capture-confounding result, and the TME.
# ---------------------------------------------------------------------------

def panel_cnv_statistic_gain(ax, old_csv, new_csv):
    """Per-sample AUC, old genome-wide statistic against the correlation statistic."""
    import pandas as pd, numpy as np
    o = pd.read_csv(old_csv); n = pd.read_csv(new_csv)
    n = n[n["status"] == "scored"]
    data = [o["auc"].to_numpy(), n["auc_signal"].to_numpy(), n["auc_corr"].to_numpy()]
    names = ["v5 genome-wide\nmean |CNV|", "signal statistic\n(new parameters)", "correlation to\nsample consensus"]
    for i, (d, nm) in enumerate(zip(data, names)):
        xj = np.random.default_rng(0).normal(i, 0.055, len(d))
        ax.scatter(xj, d, s=7, c="#4C2A85" if i == 2 else "0.6", linewidths=0, zorder=3)
        ax.hlines(np.median(d), i - 0.26, i + 0.26, color="k", lw=1.4, zorder=4)
        ax.text(i, 1.035, f"{np.median(d):.3f}", ha="center", va="bottom", fontsize=6,
                color="#4C2A85" if i == 2 else "0.35")
    ax.axhline(0.5, color="0.5", ls="--", lw=0.6)
    ax.text(2.42, 0.5, "no discrimination", fontsize=5.4, color="0.45", va="bottom", ha="right")
    ax.set_xticks(range(3)); ax.set_xticklabels(names, fontsize=5.6)
    ax.set_ylim(0.25, 1.09); ax.set_ylabel("Malignant-vs-reference AUC,\nper sample")
    _despine(ax)
    return ax


def panel_arm_scores(ax, csv, arms=("chr1p", "chr1q", "chr8p", "chr8q", "chr9p",
                                    "chr12p", "chr12q", "chr16q", "chr17p")):
    """Per-sample AUC for each chromosome arm, malignant against diploid reference.

    Shown per sample rather than pooled per study on purpose: pooling two
    immune-rich studies gave 0.61-0.72 for chr8q and chr12q, which looked like a
    clean recovery of the recurrent Ewing pattern. Across all 45 scored samples
    the same arms sit at 0.574 and 0.526, and chr16q at 0.492 shows no loss
    signal at all. The arm scores are therefore NOT what rescued the caller.
    """
    import pandas as pd, numpy as np
    d = pd.read_csv(csv)
    rng = np.random.default_rng(0)
    med = [np.nanmedian(d[a]) for a in arms]
    best = int(np.argmax(np.abs(np.array(med) - 0.5)))
    for i, a in enumerate(arms):
        v = d[a].dropna().to_numpy()
        ax.scatter(rng.normal(i, 0.07, len(v)), v, s=4,
                   c="#4C2A85" if i == best else "0.68", linewidths=0, zorder=3)
        ax.hlines(med[i], i - 0.3, i + 0.3, color="k", lw=1.2, zorder=4)
    ax.axhline(0.5, color="0.5", ls="--", lw=0.6)
    ax.annotate(f"strongest arm {arms[best].replace('chr', '')}: median {med[best]:.3f}",
                xy=(best, med[best]), xytext=(0.97, 0.03), textcoords="axes fraction",
                ha="right", va="bottom", fontsize=5.6, color="#4C2A85",
                arrowprops=dict(arrowstyle="-", color="#4C2A85", lw=0.5,
                                shrinkA=0, shrinkB=2))
    ax.set_xticks(range(len(arms)))
    ax.set_xticklabels([a.replace("chr", "") for a in arms])
    ax.set_xlabel("Chromosome arm")
    ax.set_ylabel("Malignant-vs-reference AUC,\nper sample")
    ax.set_ylim(0.2, 0.95)
    _despine(ax)
    return ax


def panel_evidence_bar(ax, csv):
    """How each frozen-malignant cell earned its label."""
    import pandas as pd
    d = pd.read_csv(csv)
    order = ["marker_and_cna", "marker_only_cna_negative", "marker_only_cna_unassessed", "cna_only"]
    nice = {"marker_and_cna": "Markers and copy number",
            "marker_only_cna_negative": "Markers only, copy number negative",
            "marker_only_cna_unassessed": "Markers only, copy number unassessed",
            "cna_only": "Copy number only (rescued)"}
    cols = {"marker_and_cna": "#4C2A85", "marker_only_cna_negative": "#9B7FC0",
            "marker_only_cna_unassessed": "0.75", "cna_only": "#D7A32B"}
    tot = d[d.evidence.isin(order)].groupby("evidence")["cells"].sum()
    y = 0
    for k in order:
        ax.barh(0, tot[k], left=y, color=cols[k], height=0.55,
                edgecolor="w", linewidth=0.4, label=nice[k])
        if tot[k] > 6000:
            ax.text(y + tot[k] / 2, 0, f"{tot[k]:,}", ha="center", va="center",
                    fontsize=5.6, color="w" if k != "marker_only_cna_unassessed" else "0.25")
        else:  # the rescued slice is too narrow to hold its own label
            ax.annotate(f"{tot[k]:,}", xy=(y + tot[k] / 2, 0.29), xytext=(y + tot[k] / 2, 0.62),
                        ha="center", va="bottom", fontsize=5.6, color="0.25",
                        arrowprops=dict(arrowstyle="-", lw=0.4, color="0.5"))
        y += tot[k]
    ax.set_yticks([]); ax.set_ylim(-0.55, 0.9)
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f"{int(v):,}"))
    ax.set_xlabel(f"Cells   ({y:,} frozen malignant in total)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.12), frameon=False,
              fontsize=5.4, ncol=4, handlelength=0.9, columnspacing=0.9, handletextpad=0.4)
    for sp in ("left", "right", "top"):
        ax.spines[sp].set_visible(False)
    return ax


def panel_mp_occupancy(ax, csv):
    """Per metaprogram: how many samples and how many studies carry it."""
    import pandas as pd, numpy as np
    d = pd.read_csv(csv).sort_values("samples", ascending=True)
    y = np.arange(len(d))
    single = d["max_study_share"].to_numpy() == 1.0
    ax.barh(y, d["samples"], color=np.where(single, "#C4654E", "#4C2A85"), height=0.68)
    for yy, (s, st) in enumerate(zip(d["samples"], d["studies"])):
        ax.text(s + 0.25, yy, f"{st}", va="center", fontsize=5.2, color="0.35")
    ax.set_yticks(y); ax.set_yticklabels(d["MP"], fontsize=5.2)
    ax.set_xlabel("Samples carrying the program\n(number at bar end = studies)")
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color="#4C2A85", label="Spans >1 study"),
                       Patch(color="#C4654E", label="Confined to one study")],
              loc="lower right", frameon=False, fontsize=5.4, handlelength=0.9, borderpad=0.1)
    _despine(ax)
    return ax


def panel_confounding_heatmap(ax, csv):
    """Candidate programs against the technical drivers that could explain them."""
    import pandas as pd, numpy as np
    d = pd.read_csv(csv).set_index("program")
    nice_row = {"score_MP28": "MP28  NR0B1/PRKCB", "score_MP24": "MP24  fibroblast-like",
                "score_MP1": "MP1  cell cycle", "score_MP3": "MP3  hypoxia-like",
                "score_MP7": "MP7  mesenchymal", "fusion_clean": "Fusion axis (cycle-filtered)"}
    nice_col = {"n_genes_by_counts": "Detected genes", "total_counts": "Total counts",
                "pct_counts_mt": "Mitochondrial %", "doublet_score": "Doublet score",
                "S_score": "S phase", "G2M_score": "G2M phase",
                "warm_dissociation_up": "Warm dissociation",
                "nucleus_vs_cell_down": "Nucleus vs cell (down)",
                "nucleus_vs_cell_up": "Nucleus vs cell (up)"}
    d = d.rename(index=nice_row, columns=nice_col)
    V = d.to_numpy()
    im = ax.imshow(V, cmap="RdBu_r", vmin=-0.8, vmax=0.8, aspect="auto")
    ax.set_xticks(range(V.shape[1])); ax.set_xticklabels(d.columns, rotation=90)
    ax.set_yticks(range(V.shape[0])); ax.set_yticklabels(d.index)
    for i in range(V.shape[0]):
        for j in range(V.shape[1]):
            if abs(V[i, j]) >= 0.3:
                ax.text(j, i, f"{V[i, j]:.2f}", ha="center", va="center", fontsize=4.8,
                        color="w" if abs(V[i, j]) > 0.55 else "0.15")
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("Spearman $\\rho$", fontsize=5.6); cb.ax.tick_params(labelsize=5)
    cb.outline.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)
    return ax


def panel_mp24_origin(ax, csv):
    """MP24 scored by how each cell entered the frozen malignant set."""
    import pandas as pd, numpy as np
    d = pd.read_csv(csv, index_col=0)
    nice = {"cna_only": "Copy number only\n(rescued fibroblasts)", "marker_and_cna": "Markers and\ncopy number",
            "marker_only_cna_negative": "Markers only,\nCNA negative",
            "marker_only_cna_unassessed": "Markers only,\nCNA unassessed"}
    d = d.loc[[i for i in ["cna_only", "marker_and_cna", "marker_only_cna_negative",
                           "marker_only_cna_unassessed"] if i in d.index]]
    y = np.arange(len(d))
    ax.barh(y, d["score_MP24"], color=["#D7A32B" if i == "cna_only" else "#9B7FC0" for i in d.index],
            height=0.6)
    for yy, v in zip(y, d["score_MP24"]):
        ax.text(v + 0.012, yy, f"{v:.3f}", va="center", fontsize=5.6, color="0.3")
    ax.set_yticks(y); ax.set_yticklabels([nice[i] for i in d.index], fontsize=5.4)
    ax.set_xlabel("Mean MP24 score")
    ax.set_xlim(0, max(d["score_MP24"]) * 1.25)
    _despine(ax)
    return ax


def panel_immune_umap(ax, adata, key="immune_subset", xy_all=None, mask=None):
    """Immune subsets on their OWN embedding, with an atlas locator.

    This must use ``X_umap_immune`` - a UMAP computed on the immune subset's own
    neighbour graph - and not the atlas-wide ``X_umap`` restricted to these
    cells. The restricted version is a cut-out of Figure 1: subsets sit wherever
    they sat in a layout optimised for 138,353 cells, gaps appear where the
    malignant mass was removed, and the locator inset becomes a duplicate of the
    panel it sits on rather than a reference map.
    """
    import numpy as np
    xy = adata.obsm["X_umap_immune"]
    lab = adata.obs[key].astype(str).to_numpy()
    order = [l for l in ["T_CD4", "T_cytotoxic", "NK", "Treg", "B_cell", "Plasma",
                         "Monocyte", "Macrophage", "Macrophage_LAM", "cDC2", "pDC",
                         "Mast", "Platelet", "Low_quality"] if l in set(lab)]
    import matplotlib.cm as cm
    cols = dict(zip(order, list(cm.tab20.colors)[:len(order)]))
    cols["Low_quality"] = (0.82, 0.82, 0.82)
    for l in order:
        m = lab == l
        ax.scatter(xy[m, 0], xy[m, 1], s=0.7, c=[cols[l]], linewidths=0, rasterized=True,
                   label=f"{l.replace('_', ' ')} ({m.sum():,})")
    _embedding_axes(ax)
    ax.legend(loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=5,
              markerscale=7, handletextpad=0.3, labelspacing=0.28, borderpad=0.1)
    if xy_all is not None and mask is not None:
        w = 0.17
        h = w * ax.get_window_extent().width / max(ax.get_window_extent().height, 1)
        ins = locator_inset(ax, xy_all, mask, loc=(0.0, 1.0 - min(h, 0.3), w, min(h, 0.3)))
        ins.set_title("in the atlas", fontsize=5, color="0.45", pad=1.5)
    return ax


def panel_immune_markers(ax, adata, key="immune_subset"):
    """Marker evidence for the immune subsets, including the two corrected calls."""
    import numpy as np, pandas as pd
    import scipy.sparse as sp
    PAN = [("T", ["CD3D", "CD3E", "IL7R"]), ("Cytotoxic", ["CD8A", "GZMK"]),
           ("NK", ["KLRF1", "GNLY", "NKG7"]), ("Treg", ["FOXP3", "CTLA4"]),
           ("B / plasma", ["MS4A1", "CD79A", "MZB1"]),
           ("Mono", ["CD14", "FCN1", "VCAN"]), ("Mac", ["C1QA", "APOE", "MRC1"]),
           ("LAM", ["TREM2", "GPNMB", "SPP1"]), ("DC", ["CD1C", "FCER1A", "LILRA4"]),
           ("Mast", ["TPSAB1", "TPSB2"]), ("Plt", ["PPBP", "PF4"])]
    genes, groups = [], []
    for g, gs in PAN:
        for x in gs:
            if x in adata.var_names:
                genes.append(x); groups.append(g)
    order = [l for l in ["T_CD4", "T_cytotoxic", "NK", "Treg", "B_cell", "Plasma", "Monocyte",
                         "Macrophage", "Macrophage_LAM", "cDC2", "pDC", "Mast", "Platelet"]
             if l in set(adata.obs[key].astype(str))]
    X = adata[:, genes].X
    X = np.asarray(X.todense()) if sp.issparse(X) else np.asarray(X)
    lab = adata.obs[key].astype(str).to_numpy()
    mean = np.vstack([X[lab == l].mean(axis=0) for l in order])
    frac = np.vstack([(X[lab == l] > 0).mean(axis=0) for l in order])
    xs, ys = np.meshgrid(np.arange(len(genes)), np.arange(len(order)))
    sc_ = ax.scatter(xs.ravel(), ys.ravel(), s=frac.ravel() * 26,
                     c=mean.ravel(), cmap="viridis", linewidths=0)
    ax.set_xticks(range(len(genes)))
    ax.set_xticklabels(genes, rotation=90, style="italic")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([l.replace("_", " ") for l in order])
    # group headers: this used to call _tier(ax, groups), but _tier formats a
    # string and takes (name, meta) — the call returned a discarded string and
    # the panel shipped with no headers at all. Drawn directly now, staggered on
    # two tiers like panel_marker_dotplot because adjacent short names collide.
    spans = []
    for g in dict.fromkeys(groups):
        idx = [i for i, x in enumerate(groups) if x == g]
        spans.append((g, (idx[0] + idx[-1]) / 2.0))
    for i, (name, x) in enumerate(spans):
        ax.text(x, -1.45 - (0.75 if i % 2 else 0.0), name, ha="center", va="bottom",
                fontsize=5.4, color="0.3")
    cax = ax.inset_axes((1.02, 0.55, 0.016, 0.42))
    cb = ax.figure.colorbar(sc_, cax=cax); cb.set_label("Mean expression", fontsize=5.4)
    cb.ax.tick_params(labelsize=5); cb.outline.set_visible(False)
    # dot area encodes the fraction of cells expressing; it needs its own key
    for f_, lab in [(0.25, "25%"), (0.75, "75%")]:
        ax.scatter([], [], s=f_ * 26, c="0.45", lw=0, label=lab)
    ax.legend(frameon=False, fontsize=5.2, loc="upper left", bbox_to_anchor=(1.005, 0.44),
              title="Cells expressing", title_fontsize=5.4, handletextpad=0.3,
              labelspacing=0.5, borderpad=0.0)
    ax.set_ylim(-2.6, len(order) - 0.3); ax.set_xlim(-0.7, len(genes) - 0.3)
    ax.invert_yaxis()
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
    ax.tick_params(length=0)
    return ax

def panel_umap_by_sample(ax, d, seed=0):
    """All 65 samples on one embedding, deliberately without a legend.

    Hues are taken evenly around the colour wheel and then shuffled, so
    neighbouring samples are unlike and the panel reads as fine speckle wherever
    samples are interleaved. A legend for 65 categories would be unreadable and
    no colour could be matched back to a sample by eye, so this panel answers
    one question only - is any sample sitting apart from the rest - and the
    per-sample grid answers which.
    """
    import numpy as np
    import matplotlib as mpl

    samples = sorted(d["sample"].astype(str).unique())
    rng = np.random.default_rng(seed)
    hues = np.linspace(0, 1, len(samples), endpoint=False)
    rng.shuffle(hues)
    cmap = mpl.colormaps["hsv"]
    cols = {s: cmap(h) for s, h in zip(samples, hues)}

    lab = d["sample"].astype(str).to_numpy()
    idx = rng.permutation(len(d))          # no study paints over another
    ax.scatter(d["umap1"].to_numpy()[idx], d["umap2"].to_numpy()[idx], s=0.25,
               c=[cols[s] for s in lab[idx]], linewidths=0, rasterized=True)
    ax.set_aspect("equal")
    _embedding_axes(ax)
    ax.text(0.99, 0.01, f"{len(samples)} samples, {len(d):,} cells", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=5.6, color="0.35")
    return ax


ASSAY_ICON = {"10x_scRNA": "cell", "10x_snRNA": "nucleus", "CELseq2_scRNA": "well"}


def _assay_glyph(ax, x, y, kind, r=0.016, color="#2B6CB0"):
    """Droplet whole-cell, droplet nucleus, or plate well, drawn to scale."""
    from matplotlib.patches import Circle, Rectangle
    asp = ax.get_window_extent().width / max(ax.get_window_extent().height, 1)
    if kind == "cell":
        ax.add_patch(Circle((x, y), r, transform=ax.transAxes, facecolor="#EAF2FA",
                            edgecolor=color, lw=0.8, zorder=3, clip_on=False))
        ax.add_patch(Circle((x, y), r * 0.42, transform=ax.transAxes, facecolor="#E8A33D",
                            edgecolor="#9C6B16", lw=0.5, zorder=4, clip_on=False))
    elif kind == "nucleus":
        ax.add_patch(Circle((x, y), r * 0.52, transform=ax.transAxes, facecolor="#E8A33D",
                            edgecolor="#9C6B16", lw=0.6, zorder=3, clip_on=False))
    else:   # plate well
        ax.add_patch(Rectangle((x - r, y - r * asp), 2 * r, 2 * r * asp, transform=ax.transAxes,
                               facecolor="#FBF1DA", edgecolor="#9C6B16", lw=0.8, zorder=3,
                               clip_on=False))
        ax.add_patch(Circle((x, y), r * 0.34, transform=ax.transAxes, facecolor="#E8A33D",
                            edgecolor="#9C6B16", lw=0.4, zorder=4, clip_on=False))


def panel_cohort_schematic(ax, g, d, totals):
    """The 'included datasets' schematic: one column per study, one row per fact.

    Modelled on NBAtlas Figure 1A. Everything plotted here is read from the
    released object, so the panel cannot drift from the atlas it describes.
    """
    import numpy as np
    from matplotlib.patches import Rectangle, FancyArrow

    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

    L0, L1, B0, B1 = 0.015, 0.235, 0.245, 0.615          # label block, study block
    ax.add_patch(Rectangle((L0, 0.02), B1 - L0, 0.95, transform=ax.transAxes,
                           facecolor="none", edgecolor="0.25", lw=0.9, clip_on=False))
    ax.add_patch(Rectangle((L0, 0.02), L1 - L0, 0.95, transform=ax.transAxes,
                           facecolor="0.90", edgecolor="none", clip_on=False, zorder=0))

    studies = list(g.index)
    xs = np.linspace(B0, B1, len(studies) + 1)
    xc = (xs[:-1] + xs[1:]) / 2
    R_STUDY, R_ASSAY, R_SAMP, R_PAT, R_CELL = 0.90, 0.775, 0.56, 0.355, 0.175

    for lbl, y in (("Study", R_STUDY), ("Assay", R_ASSAY), ("Number of\nsamples", R_SAMP + 0.045),
                   ("Patients", R_PAT), ("Number of\ncells", R_CELL)):
        ax.text(L0 + 0.014, y, lbl, transform=ax.transAxes, ha="left", va="center",
                fontsize=6.0, fontweight="bold", color="0.15", linespacing=1.15)

    # The assay key used to be three glyph rows straddling the data row, so the
    # key entries read as row labels aligned to the icons. Each study now names
    # its own assay under its glyph and no key is needed.
    ASSAY_WORD = {"cell": "single-\ncell", "nucleus": "single-\nnucleus",
                  "well": "plate\nwell"}

    smax = float(g["samples"].max())
    base, height = R_SAMP - 0.075, 0.175
    ax.plot([B0 - 0.012, B0 - 0.012], [base, base + height], color="0.35", lw=0.7,
            transform=ax.transAxes, clip_on=False)
    for v in (0, int(smax // 2), int(smax)):
        yy = base + height * v / smax
        ax.plot([B0 - 0.019, B0 - 0.012], [yy, yy], color="0.35", lw=0.7,
                transform=ax.transAxes, clip_on=False)
        ax.text(B0 - 0.024, yy, str(v), transform=ax.transAxes, ha="right", va="center",
                fontsize=4.8, color="0.3")

    cmax = float(g["cells"].max())
    for i, (lab, v) in enumerate((("1K", 1_000), ("10K", 10_000), ("50K", 50_000))):
        yy = R_CELL + 0.072 - i * 0.058
        ax.scatter([L1 - 0.022], [yy], s=max(2.5, 150 * v / cmax), c="0.45",
                   transform=ax.transAxes, clip_on=False, linewidths=0)
        ax.text(L1 - 0.040, yy, lab, transform=ax.transAxes, ha="right", va="center",
                fontsize=4.6, color="0.3")

    for i, st in enumerate(studies):
        col = SCOLOR.get(st, "#4C2A85")
        row = g.loc[st]
        ax.text(xc[i], R_STUDY, SPRETTY.get(st, st).replace(" ", "\n"), transform=ax.transAxes,
                ha="center", va="center", fontsize=6.0, fontweight="bold", color="0.1",
                linespacing=1.1)
        kind = ASSAY_ICON.get(str(row["assays"]), "cell")
        _assay_glyph(ax, xc[i], R_ASSAY + 0.038, kind)
        ax.text(xc[i], R_ASSAY - 0.048, ASSAY_WORD[kind], transform=ax.transAxes,
                ha="center", va="center", fontsize=4.6, color="0.3", linespacing=1.15)
        w = (xs[1] - xs[0]) * 0.42
        h = height * row["samples"] / smax
        ax.add_patch(Rectangle((xc[i] - w / 2, base), w, h, transform=ax.transAxes,
                               facecolor=col, edgecolor="none", clip_on=False))
        ax.text(xc[i], base + 0.012, f"{int(row['samples'])}", transform=ax.transAxes,
                ha="center", va="bottom", fontsize=5.2, color="w" if h > 0.04 else "0.2")
        ax.text(xc[i], R_PAT, f"{int(row['patients'])}", transform=ax.transAxes, ha="center",
                va="center", fontsize=6.0, color="0.15")
        s_ = max(3.0, 150 * row["cells"] / cmax)
        ax.scatter([xc[i]], [R_CELL], s=s_, c=col, transform=ax.transAxes, clip_on=False,
                   linewidths=0, zorder=3)
        k = row["cells"] / 1000
        txt = f"{k:.0f}K" if k >= 10 else f"{k:.1f}K"
        if s_ > 110:
            ax.text(xc[i], R_CELL, txt, transform=ax.transAxes, ha="center", va="center",
                    fontsize=4.8, color="w", zorder=4)
        else:
            ax.text(xc[i], R_CELL - 0.055, txt, transform=ax.transAxes, ha="center",
                    va="center", fontsize=4.8, color="0.2", zorder=4)

    # integration box
    ax.add_patch(FancyArrow(0.645, 0.50, 0.045, 0, width=0.006, head_width=0.022,
                            head_length=0.018, transform=ax.transAxes, color="0.25",
                            length_includes_head=True, clip_on=False))
    ax.add_patch(Rectangle((0.715, 0.17), 0.275, 0.64, transform=ax.transAxes,
                           facecolor="none", edgecolor="0.25", lw=0.9, clip_on=False))
    ax.text(0.8525, 0.845, "Integration", transform=ax.transAxes, ha="center", va="bottom",
            fontsize=6.2, fontweight="bold", color="0.15")
    _assay_glyph(ax, 0.762, 0.745, "cell", r=0.013)
    ax.text(0.782, 0.745, "scvi-tools", transform=ax.transAxes, ha="left", va="center",
            fontsize=6.6, color="0.1")
    ins = ax.inset_axes((0.735, 0.30, 0.235, 0.36))
    rng = np.random.default_rng(0)
    sel = rng.choice(len(d), min(25000, len(d)), replace=False)
    dd = d.iloc[sel]
    ins.scatter(dd["umap1"], dd["umap2"], s=0.35, linewidths=0, rasterized=True,
                c=[SCOLOR.get(s, "#4C2A85") for s in dd["study"].astype(str)])
    ins.set_xticks([]); ins.set_yticks([])
    for sp in ins.spines.values():
        sp.set_linewidth(0.4); sp.set_color("0.7")
    ax.text(0.8525, 0.265, f"n = {totals['studies']} studies   "
                           f"n = {totals['samples']} samples\n"
                           f"n = {totals['patients']} patients   "
                           f"n = {totals['cells']:,} cells",
            transform=ax.transAxes, ha="center", va="top", fontsize=5.4, color="0.25",
            style="italic", linespacing=1.4)
    return ax

def _hsv_map(levels, seed=0):
    """Evenly spaced hues, shuffled so neighbouring levels are unlike."""
    import numpy as np
    import matplotlib as mpl

    rng = np.random.default_rng(seed)
    h = np.linspace(0, 1, len(levels), endpoint=False)
    rng.shuffle(h)
    cmap = mpl.colormaps["hsv"]
    return {lv: cmap(x) for lv, x in zip(levels, h)}


def panel_umap_triple(ax, d, horizontal=False):
    """Three embeddings coloured by study, assay and sample.

    ``horizontal`` lays them side by side, which is what a full-width
    supplementary row needs; stacked is the original narrow-column form.
    """
    import numpy as np

    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

    rng = np.random.default_rng(0)
    idx = rng.permutation(len(d))
    dd = d.iloc[idx]
    xy = dd[["umap1", "umap2"]].to_numpy()

    blocks = [("Study", "study", {s: SCOLOR[s] for s in STUDY_ORDER}, SPRETTY, True),
              ("Assay", "assay", ASSAYC, ASSAYP, True),
              ("Sample", "sample", _hsv_map(sorted(d["sample"].astype(str).unique())), None, False)]
    for i, (title, col, cmap, pretty, legend) in enumerate(blocks):
        if horizontal:
            a = ax.inset_axes((i / 3 + 0.012, 0.06, 1 / 3 * 0.60, 0.84))
        else:
            a = ax.inset_axes((0.0, 1 - (i + 1) / 3 + 0.035, 0.60, 1 / 3 - 0.07))
        lab = dd[col].astype(str).to_numpy()
        a.scatter(xy[:, 0], xy[:, 1], s=0.12, linewidths=0, rasterized=True,
                  c=[cmap.get(x, "0.7") for x in lab])
        a.set_aspect("equal"); a.set_xticks([]); a.set_yticks([])
        for sp in a.spines.values():
            sp.set_visible(False)
        a.set_title(title, fontsize=6.4, pad=1.5, color="0.1")
        if legend:
            order = [k for k in (STUDY_ORDER if col == "study" else ASSAYP) if k in set(lab)]
            a.legend([plt.Line2D([], [], marker="o", ls="", ms=2.6, mfc=cmap[k], mec="none")
                      for k in order],
                     [pretty.get(k, k) for k in order],
                     loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=4.6,
                     handletextpad=0.25, labelspacing=0.22, borderpad=0.05)
        else:
            a.text(1.04, 0.5, f"{d[col].nunique()} samples\n(no legend:\ntoo many to read)",
                   transform=a.transAxes, fontsize=4.6, color="0.35", va="center", ha="left",
                   linespacing=1.3)
    return ax


def panel_circular_atlas(ax, d, rings=("sample", "study"), aspect_wh=1.0):
    """Annotated embedding ringed by the metadata composition of each cell type.

    Follows the NBAtlas Figure 1C construction:

    * one arc per cell type, all arcs the SAME angular width. The arcs show
      composition, not abundance - sizing them by cell count would give the
      malignant compartment two thirds of the circle and squeeze every immune
      type into a sliver, which is the opposite of what the panel is for.
    * each arc carries a solid outer band in the cell type's colour, then one
      band per metadata variable. A type drawn from one sample or one study
      reads as a solid band; a pooled type reads as fine stripes.
    * cell-type names sit outside their arc, set tangentially.
    * inside, each type's island is outlined by a density contour. The contours
      are NOT labelled: the arc label already names every type and the embedding
      is coloured by the same scale, so on-plot text repeated eight of the
      fifteen names inside a panel that already printed them outside.

    Two bands by default, not four. Assay and cluster were dropped because assay
    is already a row of panel A and cluster is Figure S3's subject; five nested
    encodings in one panel is more than a reader will decode, and the two kept -
    sample and study - are the ones that show whether a cell type is pooled or
    comes from one deposit.

    ``aspect_wh`` is the axes' width/height ratio. When the panel is wider than
    it is tall the circle stays sized by the HEIGHT and the freed width carries
    the keys, so widening the panel does not force the figure taller.
    """
    import numpy as np
    from matplotlib.patches import Wedge
    from scipy.stats import gaussian_kde

    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_ylim(-1.78, 1.78)
    ax.set_xlim(-1.78, -1.78 + 3.56 * max(1.0, float(aspect_wh)))

    xy = d[["umap1", "umap2"]].to_numpy().astype(float)
    xy = xy - (xy.min(0) + xy.max(0)) / 2
    xy = xy / np.abs(xy).max() * 0.92
    ct = d["cell_type"].astype(str).to_numpy()
    rng = np.random.default_rng(0)
    idx = rng.permutation(len(d))
    ax.scatter(xy[idx, 0], xy[idx, 1], s=0.11, linewidths=0, rasterized=True,
               c=[CT_COLOR.get(c, GREY) for c in ct[idx]])

    order = [c for c in CT_ORDER if (ct == c).sum() > 0]
    gx, gy = np.mgrid[-1.0:1.0:160j, -1.0:1.0:160j]
    grid = np.vstack([gx.ravel(), gy.ravel()])
    # Only compact classes get a contour and an on-plot label. A kernel density
    # contour around a dispersed class (Unassigned, Plasma) encloses mostly
    # empty space and asserts an island that is not there; those types are
    # identified by their arc label instead.
    for c in [x for x in order if x in DIRECT]:
        pts = xy[ct == c]
        if len(pts) < 150:
            continue
        sub = pts[rng.choice(len(pts), min(2500, len(pts)), replace=False)]
        try:
            k = gaussian_kde(sub.T, bw_method=0.35)
        except Exception:
            continue
        z = k(grid).reshape(gx.shape)
        at_pts = k(sub.T)
        lvl = float(np.quantile(at_pts, 0.14))          # contour enclosing ~86% of the type
        ax.contour(gx, gy, z, levels=[lvl], colors=["#8A7F55"], linewidths=0.45, alpha=0.9)

    # --- arcs: equal width per cell type
    gap = 2.2
    span = (360 - gap * len(order)) / len(order)
    starts = 90.0 + np.arange(len(order)) * (span + gap)
    r0, dr = 1.015, 0.052
    cmaps = {}
    for k_, var in enumerate(rings):
        lv = sorted(d[var].astype(str).unique())
        cmaps[var] = ({s: SCOLOR[s] for s in STUDY_ORDER} if var == "study" else
                      ASSAYC if var == "assay" else _hsv_map(lv, seed=k_))
    for i, c in enumerate(order):
        sel = d[ct == c]
        a0 = starts[i]
        for k_, var in enumerate(rings):
            vc = sel[var].astype(str).value_counts()
            frac = (vc / vc.sum()).to_dict()
            off = 0.0
            for lv in sorted(frac, key=lambda x: -frac[x]):
                w = span * frac[lv]
                ax.add_patch(Wedge((0, 0), r0 + (k_ + 1) * dr, a0 + off, a0 + off + w,
                                   width=dr * 0.80,
                                   facecolor=cmaps[var].get(lv, "0.7"), edgecolor="none"))
                off += w
        ax.add_patch(Wedge((0, 0), r0 + (len(rings) + 1) * dr, a0, a0 + span,
                           width=dr * 1.05, facecolor=CT_COLOR.get(c, GREY), edgecolor="none"))
        mid = a0 + span / 2
        rr = r0 + (len(rings) + 1) * dr + 0.035
        th = np.deg2rad(mid)
        east = np.cos(th) >= 0
        ax.text(rr * np.cos(th), rr * np.sin(th), PRETTY.get(c, c), fontsize=4.8, color="0.15",
                ha="left" if east else "right", va="center",
                rotation=mid if east else mid - 180, rotation_mode="anchor")

    names = {"sample": "Sample", "assay": "Assay", "study": "Study", "leiden": "Cluster"}
    from matplotlib.patches import Rectangle as _Rect
    wide = float(aspect_wh) > 1.25
    kx = 2.05 if wide else -1.66
    ky = 1.25 if wide else -1.66
    kw, kh = 0.11, 0.105
    items = list(reversed([names[v] for v in rings] + ["Cell type"]))
    for k_, nm in enumerate(items):
        ax.add_patch(_Rect((kx, ky - k_ * kh), kw, kh * 0.72,
                           facecolor=CT_COLOR["Ewing_sarcoma"] if nm == "Cell type" else "0.82",
                           edgecolor="0.45", lw=0.3))
        ax.text(kx + kw + 0.05, ky - k_ * kh + kh * 0.36, nm, fontsize=5.0, color="0.25",
                ha="left", va="center")
    ax.text(kx, ky + kh * 1.05, "Bands, outer to inner", fontsize=5.0,
            color="0.3", ha="left", va="bottom")
    if wide:
        # The freed width carries the two keys the panel needs. The cell-type key
        # is not decoration: the contours inside carry no text any more, so
        # colour is the only way to read the embedding.
        sy = ky - (len(items) + 1.2) * kh
        ax.text(kx, sy + kh * 0.95, "Study", fontsize=5.0, color="0.3", ha="left", va="bottom")
        studies = [s for s in STUDY_ORDER if s in set(d["study"].astype(str))]
        for k_, s in enumerate(studies):
            ax.add_patch(_Rect((kx, sy - k_ * kh), kw, kh * 0.72,
                               facecolor=SCOLOR[s], edgecolor="none"))
            ax.text(kx + kw + 0.05, sy - k_ * kh + kh * 0.36, SPRETTY[s], fontsize=5.0,
                    color="0.25", ha="left", va="center")
        cy = sy - (len(studies) + 1.2) * kh
        ax.text(kx, cy + kh * 0.95, "Cell type", fontsize=5.0, color="0.3",
                ha="left", va="bottom")
        half = (len(order) + 1) // 2
        for k_, c in enumerate(order):
            col, rowi = divmod(k_, half)
            x = kx + col * 1.34
            ax.add_patch(_Rect((x, cy - rowi * kh), kw, kh * 0.72,
                               facecolor=CT_COLOR.get(c, GREY), edgecolor="none"))
            ax.text(x + kw + 0.05, cy - rowi * kh + kh * 0.36, PRETTY.get(c, c),
                    fontsize=5.0, color="0.25", ha="left", va="center")
    return ax


# --- Figure 4: the CD73-high fibroblast-like state -------------------------
CAF_GRP = ["reference", "fibro-like CNA-", "candidate", "marker_malignant"]
CAF_LAB = {"reference": "Reference\n(diploid)", "fibro-like CNA-": "Fibroblast-like\ncopy number \u2212",
           "candidate": "Fibroblast-like\ncopy number +", "marker_malignant": "Marker-defined\nmalignant"}
CAF_SHORT = {"reference": "Reference", "fibro-like CNA-": "Fibroblast-like, CN \u2212",
             "candidate": "Fibroblast-like, CN +", "marker_malignant": "Marker-defined malignant"}
CAF_COL = {"reference": "#1F6FB4", "fibro-like CNA-": "#8C6D31",
           "candidate": "#D7A32B", "marker_malignant": "#5B2C83"}


def panel_caf_specificity(ax, csv):
    """Half-B own-patient copy-number specificity, selection made on half A."""
    d = pd.read_csv(csv)
    for i, g in enumerate(CAF_GRP):
        v = d.loc[d.grp == g, "spec_B"].to_numpy()
        if not len(v):
            continue
        parts = ax.violinplot([v], positions=[i], widths=0.74, showextrema=False)
        for b in parts["bodies"]:
            b.set_facecolor(CAF_COL[g]); b.set_alpha(0.55); b.set_edgecolor("none")
        q1, med, q3 = np.percentile(v, [25, 50, 75])
        ax.vlines(i, q1, q3, color="0.2", lw=1.1)
        ax.plot(i, med, "o", ms=2.6, color="w", mec="0.2", mew=0.7, zorder=5)
        ax.text(i, 1.02, f"n={len(v):,}", ha="center", va="bottom", fontsize=5.0,
                color="0.35", transform=ax.get_xaxis_transform())
    ax.axhline(0, color="0.6", lw=0.5, ls=(0, (3, 2)))
    ax.set_xticks(range(len(CAF_GRP)))
    ax.set_xticklabels(["Reference\n(diploid)", "Fibroblast-\nlike, CN \u2212",
                        "Fibroblast-\nlike, CN +", "Marker-defined\nmalignant"],
                       fontsize=5.3, linespacing=1.15)
    ax.set_ylabel("Own-patient specificity\n(held-out half of genome)")
    ax.set_ylim(-0.75, 0.95)
    for sp in ("right", "top"):
        ax.spines[sp].set_visible(False)
    return ax


def panel_caf_per_sample(ax, csv, csv_ds=None, uc_patients=("UC_Patient",), prev_csv=None):
    """Per-sample position of the state between diploid reference and tumour."""
    def wide(path):
        d = pd.read_csv(path)
        w = d.pivot_table(index=["study", "sample"], columns=["group", "subset"], values="spec_B")
        w.columns = ["_".join(c) for c in w.columns]
        return w
    w = wide(csv)
    w = w[w["fibro_like_selA"].notna()].copy()
    ds = wide(csv_ds)["fibro_like_selA"] if csv_ds is not None else None
    uc = set()
    if prev_csv is not None:
        pv = pd.read_csv(prev_csv)
        uc = set(pv.loc[pv.patient.isin(uc_patients), "sample"])
    w = w.sort_values("fibro_like_selA")
    for i, (idx, r) in enumerate(w.iterrows()):
        xs = [r["reference_all"], r["fibro_like_selA"], r["marker_malignant_all"]]
        ax.plot([min(xs), max(xs)], [i, i], "-", color="0.78", lw=0.7, zorder=1)
        for x, g in zip(xs, ["reference", "candidate", "marker_malignant"]):
            ax.plot(x, i, "o", ms=3.4, color=CAF_COL[g], mec="w", mew=0.4, zorder=3)
        if ds is not None and idx in ds.index and pd.notna(ds.loc[idx]):
            ax.plot(ds.loc[idx], i, "o", ms=5.2, mfc="none", mec=CAF_COL["candidate"],
                    mew=0.9, zorder=4)
    labs = [f"{s}" + ("  \u2020" if s in uc else "") for _, s in w.index]
    ax.set_yticks(range(len(w))); ax.set_yticklabels(labs, fontsize=5.2)
    ax.set_ylim(-1.9, len(w) - 0.3)
    ax.axvline(0, color="0.6", lw=0.5, ls=(0, (3, 2)))
    ax.set_xlabel("Own-patient specificity (held-out half of genome)")
    hs = [plt.Line2D([], [], marker="o", ls="", ms=3.4, color=CAF_COL[g], mec="w", mew=0.4)
          for g in ["reference", "candidate", "marker_malignant"]]
    hs.append(plt.Line2D([], [], marker="o", ls="", ms=5.2, mfc="none",
                         mec=CAF_COL["candidate"], mew=0.9))
    ax.legend(hs, ["Diploid reference", "Fibroblast-like, CN +", "Marker-defined malignant",
                   "Fibroblast-like, CN +, depth-matched"],
              loc="lower center", bbox_to_anchor=(0.5, -0.02), frameon=False, fontsize=5.0,
              ncol=2, handletextpad=0.3, labelspacing=0.25, columnspacing=0.8)
    ax.text(0.015, 0.985, "\u2020 urothelial carcinoma control patient", transform=ax.transAxes,
            fontsize=5.0, color="0.35", va="top")
    for sp in ("right", "top"):
        ax.spines[sp].set_visible(False)
    return ax


def panel_caf_expression(ax, csv):
    """Transcriptional profile of the same four groups."""
    d = pd.read_csv(csv).set_index("group")
    key = {"reference": "reference", "fibro-like, CNA-": "fibro-like CNA-",
           "candidate (CNA+ fibro-like)": "candidate", "marker malignant": "marker_malignant"}
    d.index = [key.get(i, i) for i in d.index]
    feats = [("NT5E_pct", "CD73 (NT5E)\ndetected", 100.0),
             ("caf_program", "Matrisome /\nCAF program", 1.0),
             ("ef_repressed", "EWS::FLI1-repressed\ngenes", 1.0),
             ("cycle", "G2/M", 1.0),
             ("doublet", "Doublet score", 1.0)]
    width = 0.2
    for j, g in enumerate(CAF_GRP):
        if g not in d.index:
            continue
        vals = [d.loc[g, f] / s for f, _, s in feats]
        ax.bar(np.arange(len(feats)) + (j - 1.5) * width, vals, width=width * 0.92,
               color=CAF_COL[g], edgecolor="w", linewidth=0.3,
               label=CAF_SHORT[g])
    ax.axhline(0, color="0.5", lw=0.5)
    ax.set_xticks(range(len(feats)))
    ax.set_xticklabels([l for _, l, _ in feats], fontsize=5.4, linespacing=1.15)
    ax.set_ylabel("Mean score, or fraction of cells")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.18), frameon=False, fontsize=5.0,
              ncol=4, handlelength=0.7, columnspacing=0.6, handletextpad=0.3)
    ax.set_ylim(min(-0.2, ax.get_ylim()[0]), 1.18)
    for sp in ("right", "top"):
        ax.spines[sp].set_visible(False)
    return ax


def panel_caf_prevalence(ax, csv, uc_patients=("UC_Patient",)):
    """How much of the malignant compartment the state accounts for, per sample."""
    d = pd.read_csv(csv)
    d = d[d.compartment >= 50].copy()
    studies = [s for s in STUDY_ORDER if s in set(d.study)]
    rng = np.random.default_rng(1)
    for i, s in enumerate(studies):
        sub = d[d.study == s]
        x = i + rng.uniform(-0.18, 0.18, len(sub))
        isuc = sub.patient.isin(uc_patients).to_numpy()
        ax.scatter(x[~isuc], sub["pct"].to_numpy()[~isuc], s=13, color=SCOLOR[s],
                   edgecolor="w", linewidth=0.3, zorder=3)
        ax.scatter(x[isuc], sub["pct"].to_numpy()[isuc], s=15, facecolor="none",
                   edgecolor=SCOLOR[s], linewidth=0.8, zorder=3, marker="s")
        ax.text(i, 1.005, f"{int((sub['candidate'] >= 10).sum())}/{len(sub)}", ha="center",
                va="bottom", fontsize=5.0, color="0.35",
                transform=ax.get_xaxis_transform())
    ax.set_xticks(range(len(studies)))
    ax.set_xticklabels([SPRETTY[s] for s in studies], fontsize=5.0,
                       rotation=28, ha="right", rotation_mode="anchor")
    ax.set_xlim(-0.6, len(studies) - 0.4)
    ax.set_ylabel("Per cent of malignant\ncompartment")
    ax.text(0.5, 1.085, "samples with \u226510 cells of the state", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=5.0, color="0.35")
    ax.text(0.985, 0.90, "\u25a1 urothelial carcinoma\n    control patient",
            transform=ax.transAxes, ha="right", va="top", fontsize=5.0, color="0.35",
            linespacing=1.2)
    for sp in ("right", "top"):
        ax.spines[sp].set_visible(False)
    return ax

NORMAL_ORDER = ["T_cell", "NK", "B_cell", "Plasma", "Monocyte_Macrophage", "DC", "pDC",
                "Mast", "Endothelial", "Fibroblast_CAF", "Urothelial"]


# --- Figure 6: the atlas as a reference ------------------------------------
TARGET_ORDER = ["STEAP1", "IL1RAP", "CD276", "CD99", "NT5E"]
TARGET_LAB = {"STEAP1": "STEAP1", "IL1RAP": "IL1RAP", "CD276": "CD276\n(B7-H3)",
              "CD99": "CD99", "NT5E": "NT5E\n(CD73)"}


def panel_target_landscape(ax, csv):
    """Per-cell-type detection of five surface targets, pooled across the atlas."""
    d = pd.read_csv(csv, index_col=0)
    order = [c for c in CT_ORDER if c in d.index]
    pct = d.loc[order, [f"{g}_pct" for g in TARGET_ORDER]].to_numpy()
    mean = d.loc[order, [f"{g}_mean" for g in TARGET_ORDER]].to_numpy()
    xs, ys = np.meshgrid(np.arange(len(order)), np.arange(len(TARGET_ORDER)))
    sc = ax.scatter(xs.ravel(), ys.ravel(), s=(pct.T.ravel() / 100) * 95 + 2,
                    c=mean.T.ravel(), cmap="Purples", vmin=0,
                    vmax=float(np.nanpercentile(mean, 98)), edgecolor="0.4", linewidth=0.3)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([PRETTY.get(c, c) for c in order], fontsize=5.2,
                       rotation=32, ha="right", rotation_mode="anchor")
    ax.set_yticks(range(len(TARGET_ORDER)))
    ax.set_yticklabels([TARGET_LAB[g] for g in TARGET_ORDER], fontsize=5.6, linespacing=1.1)
    ax.set_xlim(-0.6, len(order) - 0.4); ax.set_ylim(-0.7, len(TARGET_ORDER) - 0.3)
    ax.invert_yaxis()
    cb = ax.figure.colorbar(sc, ax=ax, fraction=0.022, pad=0.012)
    cb.set_label("Mean log-normalized expression", fontsize=5.2)
    cb.ax.tick_params(labelsize=5.0, length=1.6, width=0.4)
    for frac in (25, 50, 75):
        ax.scatter([], [], s=(frac / 100) * 95 + 2, c="0.75", edgecolor="0.4",
                   linewidth=0.3, label=f"{frac}%")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, 1.17), frameon=False, fontsize=5.0,
              ncol=3, handletextpad=0.25, columnspacing=0.8, title="Cells positive",
              title_fontsize=5.0)
    for sp in ("right", "top"):
        ax.spines[sp].set_visible(False)
    return ax


def panel_target_margin(ax, csv):
    """Within-sample selectivity: malignant minus normal, in percentage points."""
    d = pd.read_csv(csv)
    piv = d.pivot_table(index="gene", columns="normal_type", values="delta",
                        aggfunc="median").reindex(TARGET_ORDER)
    cols = [c for c in NORMAL_ORDER if c in piv.columns]
    piv = piv[cols]
    v = float(np.nanmax(np.abs(piv.to_numpy())))
    im = ax.imshow(piv.to_numpy(), cmap="RdBu_r", vmin=-v, vmax=v, aspect="auto")
    nsamp = d.groupby("normal_type")["sample"].nunique()
    ax.set_xticks(range(len(cols)))
    short = {"T_cell": "T", "NK": "NK", "B_cell": "B", "Plasma": "Plasma",
             "Monocyte_Macrophage": "Mono", "DC": "DC", "pDC": "pDC", "Mast": "Mast",
             "Endothelial": "Endo", "Fibroblast_CAF": "Fib", "Urothelial": "Uro"}
    ax.set_xticklabels([f"{short.get(c, c)}\n({nsamp[c]})" for c in cols], fontsize=5.2,
                       rotation=0, ha="center", linespacing=1.1)
    ax.set_yticks(range(len(piv)))
    ax.set_yticklabels([TARGET_LAB[g] for g in piv.index], fontsize=5.6, linespacing=1.1)
    for i in range(piv.shape[0]):
        for j in range(piv.shape[1]):
            val = piv.to_numpy()[i, j]
            if np.isfinite(val):
                lab = "0" if abs(val) < 0.5 else f"{val:+.0f}"
                ax.text(j, i, lab, ha="center", va="center", fontsize=4.8,
                        color="w" if abs(val) > 0.6 * v else "0.2")
    cb = ax.figure.colorbar(im, ax=ax, fraction=0.028, pad=0.015)
    cb.set_label("Percentage points", fontsize=5.2)
    cb.ax.tick_params(labelsize=5.0, length=1.6, width=0.4)
    ax.tick_params(length=0)
    return ax


def panel_loso_accuracy(ax, csv, csv_type=None):
    """Label transfer onto a study the reference never saw."""
    d = pd.read_csv(csv).sort_values("accuracy")
    y = np.arange(len(d))
    ax.barh(y, d["accuracy"], height=0.52, color=[SCOLOR.get(s, GREY) for s in d.held_out],
            edgecolor="w", linewidth=0.4)
    for i, (_, r) in enumerate(d.iterrows()):
        ax.text(r["accuracy"] + 0.012, i, f"{r['accuracy']:.2f}", va="center",
                fontsize=5.2, color="0.25")
    ax.set_yticks(y)
    ax.set_yticklabels([f"{SPRETTY.get(s, s)}\n({n:,} cells)" for s, n in
                        zip(d.held_out, d.n_query)], fontsize=5.2, linespacing=1.15)
    ax.set_xlim(0, 1.08); ax.set_xlabel("Agreement with full-atlas annotation")
    med = float(d["accuracy"].median())
    ax.axvline(med, color="0.35", lw=0.7, ls=(0, (3, 2)))
    ax.text(med, len(d) - 0.35, f"median {med:.2f}", fontsize=5.0, color="0.35",
            ha="center", va="bottom")
    for sp in ("right", "top"):
        ax.spines[sp].set_visible(False)
    return ax


def panel_bulk_survival(ax, csv, frac_csv=None):
    """Hazard ratios per standard deviation of atlas-derived composition."""
    d = pd.read_csv(csv)
    d = d[d.cohort == "both (stratified)"].copy()
    lab = {"malignant": "Malignant fraction", "stromal": "Stromal fraction",
           "immune": "Immune fraction", "Fibroblast_CAF": "Fibroblast / CAF fraction"}
    d["key"] = d.covariate.map(lab)
    rows = [(k, e) for k in lab.values() for e in ["OS", "EFS"]]
    cmap = {"OS": "#4C2A85", "EFS": "#D7A32B"}
    for i, (k, e) in enumerate(rows):
        r = d[(d.key == k) & (d.endpoint == e)]
        if not len(r):
            continue
        r = r.iloc[0]
        ax.plot([r.ci_low, r.ci_high], [i, i], "-", color=cmap[e], lw=1.0)
        ax.plot(r.hr, i, "o", ms=3.6, color=cmap[e], mec="w", mew=0.4, zorder=3)
        ax.text(1.02, i, f"{r.hr:.2f} ({r.ci_low:.2f}\u2013{r.ci_high:.2f})  p={r.p:.3f}",
                transform=ax.get_yaxis_transform(), va="center", fontsize=4.9, color="0.25")
    ax.axvline(1.0, color="0.5", lw=0.6, ls=(0, (3, 2)))
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels([f"{k}, {e}" for k, e in rows], fontsize=5.2)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xticks([0.6, 0.8, 1.0, 1.25, 1.6])
    ax.set_xticklabels(["0.6", "0.8", "1.0", "1.25", "1.6"], fontsize=5.2)
    ax.set_xlabel("Hazard ratio per standard deviation")
    hs = [plt.Line2D([], [], marker="o", ls="-", ms=3.6, color=cmap[e], mec="w", mew=0.4)
          for e in ["OS", "EFS"]]
    ax.legend(hs, ["Overall survival", "Event-free survival"], loc="lower left",
              frameon=False, fontsize=5.0, handletextpad=0.3, labelspacing=0.25)
    for sp in ("right", "top"):
        ax.spines[sp].set_visible(False)
    return ax
