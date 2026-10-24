"""Stage 03a - select HVGs and export the compact training object.

    python src/pipeline/03_integrate.py --export

scVI only ever sees the highly variable genes, so the object shipped to the GPU
carries those columns alone: 138,353 x 5,000 instead of 138,353 x 16,622, which
is the difference between a ~400 MB upload and a ~1.8 GB one. The full feature
space stays local and the latent spaces are joined back onto it by cell name.

HVG selection is at STUDY level rather than sample level: per-sample
``seurat_v3`` fails with a loess singularity on the small CEL-seq2 plates, and
study-level selection also stops any single deposit from dominating the gene
set. ``seurat_v3`` needs raw counts, which is why it runs on ``layers['counts']``
and not on the normalized matrix.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import anndata as ad
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ewsatlas import config, io, provenance  # noqa: E402
from ewsatlas.compat import sequential_joblib  # noqa: E402


def export_training_object(cfg: dict) -> Path:
    interim = config.path(cfg, "interim")
    src = interim / "atlas_preintegration.h5ad"
    with provenance.start_run("03a_hvg_export", cfg, [src]) as run:
        A = ad.read_h5ad(src)
        n_top = int(config.get(cfg, "integration.hvg.n_top_genes"))
        sc.pp.highly_variable_genes(
            A, n_top_genes=n_top,
            flavor=str(config.get(cfg, "integration.hvg.flavor")),
            batch_key=str(config.get(cfg, "integration.hvg.batch_key")),
            layer="counts", subset=False)
        hvg = A.var_names[A.var["highly_variable"].to_numpy()]
        keep_obs = ["sample", "study", "assay", "patient", "material_type", "disease_state",
                    "treatment_timing", "sort_strategy", "pct_counts_mt", "pct_counts_ribo",
                    "n_genes_by_counts", "total_counts"]
        T = ad.AnnData(
            X=A[:, hvg].layers["counts"].copy(),
            obs=A.obs[[c for c in keep_obs if c in A.obs.columns]].copy(),
            var=A.var.loc[hvg, [c for c in ("highly_variable", "means", "variances")
                                if c in A.var.columns]].copy(),
        )
        T.layers["counts"] = T.X.copy()
        out = interim / "atlas_hvg_for_scvi.h5ad"
        io.write_h5ad_atomic(T, out)
        # how much of the selection is shared across studies?
        run.add(n_hvg=int(len(hvg)), cells=int(T.n_obs),
                hvg_nbatches=A.var.loc[hvg, "highly_variable_nbatches"].value_counts().to_dict()
                if "highly_variable_nbatches" in A.var.columns else {},
                out_bytes=int(out.stat().st_size))
        print(f"HVG export: {T.n_obs:,} x {T.n_vars:,} -> {out.name} "
              f"({out.stat().st_size/1e6:.0f} MB)")
        if "highly_variable_nbatches" in A.var.columns:
            print("HVGs by number of studies selecting them:",
                  A.var.loc[hvg, "highly_variable_nbatches"].value_counts().sort_index().to_dict())
    return out


EMBEDDINGS = {
    "scVI_batch_sample": "latent_scvi_sample.npy",
    "scVI_batch_study": "latent_scvi_study.npy",
    "Harmony_sample": "latent_harmony_sample.npy",
    "Unintegrated_PCA": "pca_unintegrated.npy",
}


def benchmark(cfg: dict) -> Path:
    """Stage 03b - score the candidate embeddings on one footing.

    Batch mixing is measured on every cell as the share of each cell's 15
    nearest neighbours drawn from its own study, against the share expected if
    studies were mixed at random (sum of squared study proportions). This is the
    statistic already reported for v4, so the two are comparable.

    Biological conservation needs labels, and the labels used here are the
    DEPOSITORS' own annotations - OpenScPCA for the Waltner nuclei, the
    Goodspeed cell metadata - never labels derived from any of these
    embeddings. Scoring an embedding against clusters taken from that same
    embedding would be circular. Those labels cover 69,400 of 138,353 cells, so
    the bio metrics carry a smaller denominator than the batch metrics and both
    are reported.
    """
    import json

    import numpy as np
    import pandas as pd

    interim = config.path(cfg, "interim")
    emb_dir = interim / "embeddings"
    lab_path = config.ROOT / "results" / "tables" / "v5_depositor_labels.csv"
    src = interim / "atlas_preintegration.h5ad"

    with provenance.start_run("03b_benchmark", cfg, [src, lab_path]) as run:
        A = ad.read_h5ad(src, backed="r")
        obs = A.obs[["sample", "study"]].astype(str).copy()
        A.file.close()
        cells = np.load(emb_dir / "cell_index.npy", allow_pickle=True).astype(str)
        if list(obs.index) != list(cells):
            raise AssertionError("cell order differs between the atlas object and the GPU output")

        labels = pd.read_csv(lab_path, index_col=0).reindex(obs.index)
        B = ad.AnnData(X=np.zeros((len(obs), 1), dtype=np.float32), obs=obs.copy())
        B.obs["label"] = labels["depositor_label"].astype(str).to_numpy()
        for name, fname in EMBEDDINGS.items():
            B.obsm[name] = np.load(emb_dir / fname).astype(np.float32)

        st = B.obs["study"].to_numpy()
        shares = pd.Series(st).value_counts(normalize=True)
        expected = float((shares ** 2).sum())
        results = {"expected_same_study_knn": expected, "n_cells": int(B.n_obs)}
        for name in EMBEDDINGS:
            with sequential_joblib():  # see ewsatlas.compat
                sc.pp.neighbors(B, use_rep=name, n_neighbors=15, key_added=f"n_{name}")
            g = B.obsp[f"n_{name}_connectivities"].tocsr()
            same = np.empty(B.n_obs)
            for i in range(B.n_obs):
                nb = g.indices[g.indptr[i]:g.indptr[i + 1]]
                same[i] = np.mean(st[nb] == st[i]) if nb.size else np.nan
            results[name] = {
                "same_study_knn": float(np.nanmean(same)),
                "by_study": {s: float(np.nanmean(same[st == s])) for s in shares.index},
            }
            print(f"{name:20s} same-study kNN {results[name]['same_study_knn']:.3f} "
                  f"(expected {expected:.3f})", flush=True)

        try:
            from scib_metrics.benchmark import BatchCorrection, Benchmarker, BioConservation

            sub = B[B.obs["label"].to_numpy() != "unlabelled"].copy()
            bm = Benchmarker(
                sub, batch_key="study", label_key="label",
                embedding_obsm_keys=list(EMBEDDINGS),
                pre_integrated_embedding_obsm_key="Unintegrated_PCA",
                bio_conservation_metrics=BioConservation(
                    isolated_labels=True, nmi_ari_cluster_labels_leiden=True,
                    nmi_ari_cluster_labels_kmeans=False, silhouette_label=True, clisi_knn=True),
                batch_correction_metrics=BatchCorrection(
                    bras=True, ilisi_knn=True, kbet_per_label=True,
                    graph_connectivity=True, pcr_comparison=True, sbee=True),
                n_jobs=1)
            with sequential_joblib():
                bm.benchmark()
            tbl = bm.get_results(min_max_scale=False)
            out_tbl = config.ROOT / "results" / "tables" / "v5_scib_results.csv"
            tbl.to_csv(out_tbl)
            results["scib_n_cells"] = int(sub.n_obs)
            results["scib_table"] = str(out_tbl.name)
            print(tbl.to_string(), flush=True)
        except Exception as exc:  # recorded, never silently skipped
            results["scib_error"] = f"{type(exc).__name__}: {exc}"
            print("scib-metrics failed:", results["scib_error"], flush=True)

        out = config.ROOT / "results" / "tables" / "v5_knn_mixing.json"
        out.write_text(json.dumps(results, indent=1), encoding="utf-8")
        run.add(**{k: v for k, v in results.items() if not isinstance(v, dict)})
    return out


def mixing_by_celltype() -> "Path":
    """Per-cell-type same-study kNN mixing, read off the released atlas.

    This used to be computed ad hoc, so when the annotation was revised from 12
    to 15 classes the table kept its old labels and the figure built from it
    silently showed a mast-cell population under the dendritic-cell name. It is
    a pipeline step now so it cannot fall behind the object again.

    Observed = mean over cells of a type of the share of that cell's 15 nearest
    neighbours (in the released scVI space) drawn from its own study. Expected =
    the same quantity under random mixing within the type, i.e. the probability
    that two cells of the type drawn at random share a study.
    """
    import numpy as np
    import pandas as pd
    import scanpy as sc

    cfg = config.load()
    atlas = config.path(cfg, "processed") / "EWSAtlas_v5.h5ad"
    with provenance.start_run("03c_mixing_by_celltype", cfg, [atlas]) as run:
        A = ad.read_h5ad(atlas)
        st = A.obs["study"].astype(str).to_numpy()
        ct = A.obs["cell_type"].astype(str).to_numpy()
        with sequential_joblib():  # see ewsatlas.compat
            sc.pp.neighbors(A, use_rep="X_scvi", n_neighbors=15, key_added="n_mix")
        g = A.obsp["n_mix_connectivities"].tocsr()
        same = np.empty(A.n_obs)
        for i in range(A.n_obs):
            nb = g.indices[g.indptr[i]:g.indptr[i + 1]]
            same[i] = np.mean(st[nb] == st[i]) if nb.size else np.nan

        rows = []
        for c in pd.unique(ct):
            m = ct == c
            shares = pd.Series(st[m]).value_counts(normalize=True)
            rows.append({
                "cell_type": c,
                "cells": int(m.sum()),
                "observed": round(float(np.nanmean(same[m])), 3),
                "expected_given_composition": round(float((shares ** 2).sum()), 3),
                "n_studies": int(shares.size),
            })
        tbl = pd.DataFrame(rows).sort_values("cells", ascending=False)
        out = config.ROOT / "results" / "tables" / "v5_mixing_by_celltype.csv"
        tbl.to_csv(out, index=False)
        run.add(n_cell_types=int(len(tbl)), n_cells=int(A.n_obs))
        print(tbl.to_string(index=False), flush=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--export", action="store_true", help="write the HVG training object")
    ap.add_argument("--mixing-by-celltype", action="store_true",
                    help="per-cell-type same-study kNN mixing from the released atlas")
    ap.add_argument("--benchmark", action="store_true", help="score the candidate embeddings")
    args = ap.parse_args()
    cfg = config.load()
    if args.export:
        export_training_object(cfg)
        return 0
    if args.mixing_by_celltype:
        mixing_by_celltype()
        return 0
    if args.benchmark:
        benchmark(cfg)
        return 0
    ap.print_help()
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
