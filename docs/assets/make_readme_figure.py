"""Generate the README hero figure from the released atlas object.

Not a pipeline stage: reads only `obs` and `obsm` from
data/processed/EWSAtlas_v5.h5ad and writes docs/assets/README_atlas_umap.png.
Run from the repository root:  python docs/assets/make_readme_figure.py
"""
from pathlib import Path
import h5py, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

ROOT = Path(__file__).resolve().parents[2]
H5 = ROOT / "data/processed/EWSAtlas_v5.h5ad"
OUT = ROOT / "docs/assets/README_atlas_umap.png"

CT = {"Ewing_sarcoma": "#C0392B", "Ewing_proliferating": "#F5A58C",
      "Fibroblast_CAF": "#8C510A", "Endothelial": "#6A3D9A", "Urothelial": "#117733",
      "T_cell": "#1F78B4", "NK": "#7FC7E8", "B_cell": "#0B4F6C", "Plasma": "#44AA99",
      "Monocyte_Macrophage": "#DDAA33", "DC": "#AA4499", "pDC": "#882255",
      "Mast": "#999933", "Platelet": "#332288", "Unassigned": "#D9D9D9"}
PRETTY = {"Ewing_sarcoma": "Ewing sarcoma", "Ewing_proliferating": "Ewing, proliferating",
          "Fibroblast_CAF": "Fibroblast / CAF", "Monocyte_Macrophage": "Monocyte / macrophage",
          "T_cell": "T cell", "B_cell": "B cell", "pDC": "pDC", "DC": "DC", "NK": "NK",
          "Mast": "Mast", "Plasma": "Plasma", "Platelet": "Platelet",
          "Endothelial": "Endothelial", "Urothelial": "Urothelial", "Unassigned": "Unassigned"}
LAB = {"Ewing_sarcoma": (6.3, 4.2), "Ewing_proliferating": (15.2, 9.9),
       "Monocyte_Macrophage": (12.2, 15.8), "T_cell": (6.4, 10.9), "Endothelial": (4.4, 18.8),
       "NK": (3.4, 14.7), "Plasma": (-0.8, 16.2), "B_cell": (-4.8, 11.3), "DC": (8.3, 15.1),
       "Unassigned": (17.8, 13.9), "Urothelial": (15.6, 6.7), "Mast": (18.4, 3.0),
       "Fibroblast_CAF": (5.3, -5.5), "pDC": (-0.6, -3.0), "Platelet": (-5.8, 3.3)}


# Label anchors: densest 120x120 bin of each class on a grid common to the whole
# embedding, snapped to a real cell. Unassigned is scattered atlas-wide, so it is
# anchored on its one distinct island rather than on a central bin.
ANCHOR = {
    "B_cell": (-0.52, 11.72),
    "DC": (10.05, 11.92),
    "Endothelial": (4.36, 16.87),
    "Ewing_proliferating": (11.96, 8.71),
    "Ewing_sarcoma": (6.46, 5.27),
    "Fibroblast_CAF": (5.41, -3.78),
    "Mast": (19.56, 5.03),
    "Monocyte_Macrophage": (12.6, 12.78),
    "NK": (4.56, 13.44),
    "Plasma": (0.51, 14.07),
    "Platelet": (-4.73, 4.62),
    "T_cell": (4.77, 11.07),
    "Unassigned": (16.9, 12.5),
    "Urothelial": (14.92, 5.26),
    "pDC": (2.02, -1.83),
}

def load():
    with h5py.File(H5, "r") as f:
        um = f["obsm"]["X_umap"][:]
        out = {"umap1": um[:, 0], "umap2": um[:, 1]}
        g = f["obs"]["cell_type"]
        cats = [c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]]
        out["cell_type"] = pd.Categorical.from_codes(g["codes"][:], cats)
    return pd.DataFrame(out)


def main():
    plt.rcParams.update({"font.size": 9, "savefig.bbox": "tight",
                         "font.family": "sans-serif", "pdf.fonttype": 42})
    df = load()

    fig, ax = plt.subplots(figsize=(6.4, 6.2))
    for k in df.cell_type.value_counts().index:
        g = df[df.cell_type == k]
        ax.scatter(g.umap1, g.umap2, s=0.5, c=CT[k], alpha=.55, lw=0, rasterized=True)
    for k in CT:
        x, y = LAB[k]; ax_, ay = ANCHOR[k]
        if np.hypot(x - ax_, y - ay) > 1.4:
            ax.plot([ax_, x], [ay, y], lw=.45, color="0.35", zorder=3, solid_capstyle="round")
        ax.text(x, y, PRETTY[k], fontsize=8, ha="center", va="center", zorder=4,
                color="0.12", path_effects=[pe.withStroke(linewidth=2.6, foreground="white")])
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xlim(-7.4, 21.4); ax.set_ylim(-7.2, 20.8); ax.set_aspect("equal")
    x0, y0 = -6.6, -6.4
    for dx, dy in ((4.2, 0), (0, 4.2)):
        ax.annotate("", xy=(x0 + dx, y0 + dy), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="-|>", lw=.8, color="0.3"))
    ax.text(x0 + 4.4, y0, "UMAP 1", fontsize=7, va="center", color="0.3")
    ax.text(x0, y0 + 4.5, "UMAP 2", fontsize=7, ha="center", color="0.3")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=400)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
