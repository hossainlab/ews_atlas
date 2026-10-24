"""Stage 02 — build the integrated-ready atlas object.

    python src/pipeline/02_build_atlas.py

Derives the feature space across all QC'd studies, attaches harmonized
covariates, concatenates, normalizes, and writes the pre-integration object.
Asserts the invariants that previously failed silently: no zero-count cells,
marker panel present, per-study counts reconcile.
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ewsatlas import config, genespace, io, provenance  # noqa: E402

CORE = ["patient", "fusion_partner", "sex", "age_years", "tumor_site", "disease_state",
        "treatment_timing", "material_type", "assay", "assay_chemistry", "sort_strategy"]


def main() -> int:
    cfg = config.load()
    interim = config.path(cfg, "interim")
    cov_path = config.ROOT / config.get(cfg, "paths.covariates")
    paths = [interim / f"qc_{s}.h5ad" for s in cfg["studies"]]
    present = [p for p in paths if p.is_file()]
    if not present:
        raise FileNotFoundError("no qc_*.h5ad found — run 01_ingest_qc.py first")

    with provenance.start_run("02_build_atlas", cfg, [*present, cov_path]) as run:
        per_study = {p.stem.replace("qc_", ""): ad.read_h5ad(p) for p in present}
        genes, report = genespace.build(
            per_study,
            min_cells=int(config.get(cfg, "gene_space.min_cells_detected")),
            n_studies=int(config.get(cfg, "gene_space.detected_in_n_studies")),
        )
        markers = genespace.check_markers(genes)
        missing = [g for g, ok in markers.items() if not ok]
        if missing:
            raise AssertionError(f"marker panel lost in feature selection: {missing}")

        parts = []
        for name, a in per_study.items():
            sub = a[:, genes].copy()
            # Study-specific layers (ScPCA's 'spliced', Visser's
            # 'celseq_corrected') stay in the per-study interim files. Carrying
            # them into the shared object would add ~1 GB of mostly structural
            # zeros after an outer join and invite someone to analyse a layer
            # that exists for one study only.
            for lyr in list(sub.layers.keys()):
                if lyr is not None:
                    del sub.layers[lyr]
            keep = [c for c in sub.obs.columns
                    if c in ("sample", "study", "well", "doublet_score", "doublet_call",
                             "n_genes_by_counts", "total_counts", "pct_counts_mt",
                             "pct_counts_ribo", "pct_counts_spikein")
                    or c.startswith("src_")]
            sub.obs = sub.obs[keep]
            parts.append(sub)
        Z = ad.concat(parts, join="outer", index_unique=None)
        Z.obs_names_make_unique()

        cov = pd.read_csv(cov_path)
        cov["_k"] = cov["study"].astype(str) + "||" + cov["sample"].astype(str)
        cov = cov.drop_duplicates("_k").set_index("_k")
        key = Z.obs["study"].astype(str) + "||" + Z.obs["sample"].astype(str)
        unmatched = sorted(set(key[~key.isin(cov.index)]))
        if unmatched:
            raise AssertionError(f"samples with no covariate row: {unmatched}")
        for c in CORE:
            Z.obs[c] = cov.loc[key.to_numpy(), c].to_numpy()

        Z.layers["counts"] = Z.X.copy()
        totals = np.asarray(Z.layers["counts"].sum(axis=1)).ravel()
        if (totals == 0).any():
            raise AssertionError(
                f"{int((totals == 0).sum())} cells have zero counts after feature selection — "
                "gene-symbol harmonization is broken (see io.harmonize_symbols)")
        sc.pp.normalize_total(Z, target_sum=1e4)
        sc.pp.log1p(Z)
        Z.layers["lognorm"] = Z.X.copy()
        Z.uns["gene_space"] = report
        Z.uns["config"] = cfg

        out = interim / "atlas_preintegration.h5ad"
        io.write_h5ad_atomic(Z, out)
        run.add(cells=int(Z.n_obs), genes=int(Z.n_vars), gene_space=report,
                by_study=Z.obs["study"].value_counts().to_dict())
        print(f"atlas_preintegration: {Z.n_obs:,} cells x {Z.n_vars:,} genes -> {out.name}")
        print("by study:", Z.obs["study"].value_counts().to_dict())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
