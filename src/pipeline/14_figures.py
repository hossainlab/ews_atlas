"""Stage 14 - the publication figure set.

    python src/pipeline/14_figures.py [--prepare-only]

Runs last: it draws on the atlas object and on tables written by stages 01-12.
Writes six main and seven supplementary figures to ``paths.figures``
(``results/figures/``), each as a 300 dpi PNG and a vector PDF (editable text,
dense scatter layers rasterized), plus ``figure_stats.json``. Captions are kept
by hand in ``results/figures/figure_captions.md``. Per-panel intermediates go
to ``paths.figure_panels``.

The sequence follows the arc shared by published harmonized atlases (GBmap,
the integrated breast-cancer atlas, HLCA, LuCA, NBAtlas): construction ->
validation -> composition -> compartments -> reference use -> clinical
association.

    1  Building the atlas              S1 QC policy
    2  Integration and annotation      S2 deposit-specific QC
    3  Composition across samples      S3 integration detail
    4  Malignant compartment           S4 covariate coverage
    5  Immune microenvironment         S5 copy-number detail
    6  The atlas as a reference        S6 embedding by sample
                                       S7 CAF-like malignant state

No panel carries a title; every number drawn is read from the atlas or a
results table. The malignant-only embedding and program scores are computed
once and cached at ``paths.malignant_embedding`` (``--prepare-only`` builds
just that cache, about three minutes on CPU).
"""

from __future__ import annotations

import argparse
import json
import pickle
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ewsatlas import config, provenance  # noqa: E402
from ewsatlas import figures as F  # noqa: E402
from ewsatlas import layout as L  # noqa: E402


def _p(cfg, key) -> Path:
    return config.ROOT / config.get(cfg, f"paths.{key}")


def prepare_malignant(cfg, atlas: Path, cache: Path) -> pd.DataFrame:
    """UMAP of the frozen malignant cells on their own scVI neighbour graph,
    plus per-cell scores for the recurrent programs shown in Figure 4."""
    import scanpy as sc
    from ewsatlas.compat import sequential_joblib

    if cache.is_file():
        return pd.read_parquet(cache)
    programs = config.get(cfg, "figures.malignant_programs")
    A = ad.read_h5ad(atlas, backed="r")
    mask = A.obs["malignant_frozen"].to_numpy().astype(bool)
    M = A[mask].to_memory()
    A.file.close()
    M.X = M.layers["lognorm"]
    seed = config.get(cfg, "seed")
    with sequential_joblib():
        sc.pp.neighbors(M, use_rep="X_scvi", random_state=seed,
                        n_neighbors=config.get(cfg, "figures.malignant_embedding_n_neighbors"))
    sc.tl.umap(M, random_state=seed)
    with open(_p(cfg, "metaprograms"), "rb") as fh:
        mps = pickle.load(fh)
    out = pd.DataFrame(M.obsm["X_umap"], index=M.obs_names, columns=["umap1", "umap2"])
    for mp in programs:
        genes = [g for g in mps[mp] if g in M.var_names]
        sc.tl.score_genes(M, genes, score_name=f"score_{mp}", random_state=seed)
        out[mp] = M.obs[f"score_{mp}"].to_numpy()
    out["study"] = M.obs["study"].astype(str).to_numpy()
    out["sample"] = M.obs["sample"].astype(str).to_numpy()
    cache.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(cache)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--prepare-only", action="store_true",
                    help="only build the malignant embedding cache")
    args = ap.parse_args()
    cfg = config.load()
    atlas = _p(cfg, "atlas_out")
    cache = _p(cfg, "malignant_embedding")
    if args.prepare_only:
        with provenance.start_run("14_figures_prepare", cfg, [atlas]) as run:
            emb = prepare_malignant(cfg, atlas, cache)
            run.add(n_malignant=len(emb))
        return 0
    return build(cfg, atlas, cache)


def outlines(umap_aspect: float) -> dict:
    """Grid layout of every figure: (letter, row, col, colspan) on 12 columns."""
    o = L.outline
    return {
        "Figure_1_atlas_construction": o("fig1", [22, 80, 56], [
            ("a", 0, 0, 12), ("b", 1, 0, 8), ("c", 1, 8, 4), ("d", 2, 0, 4), ("e", 2, 4, 8)]),
        "Figure_2_integration_and_annotation": o("fig2", [52, 62], [
            ("a", 0, 0, 5), ("b", 0, 5, 7), ("c", 1, 0, 6), ("d", 1, 6, 6)]),
        "Figure_3_composition_across_samples": o("fig3", [80, 46], [
            ("a", 0, 0, 12), ("b", 1, 0, 5), ("c", 1, 5, 7)]),
        "Figure_4_malignant_compartment": o("fig4", [60, 26, 68], [
            ("a", 0, 0, 4), ("b", 0, 4, 4), ("c", 0, 8, 4), ("d", 1, 0, 12),
            ("e", 2, 0, 5), ("f", 2, 5, 7)]),
        "Figure_5_immune_microenvironment": o("fig5", [62, 50, 54], [
            ("a", 0, 0, 12), ("b", 1, 0, 12), ("c", 2, 0, 12)]),
        "Figure_6_atlas_as_reference": o("fig6", [43, 50], [
            ("a", 0, 0, 7), ("b", 0, 7, 5), ("c", 1, 0, 3), ("d", 1, 3, 4), ("e", 1, 7, 5)]),
        "Figure_S1_quality_control": o("figS1", [52, 44], [
            ("a", 0, 0, 12), ("b", 1, 0, 6), ("c", 1, 6, 6)]),
        "Figure_S2_deposit_specific_qc": o("figS2", [64], [("a", 0, 0, 6), ("b", 0, 6, 6)]),
        "Figure_S3_integration_detail": o("figS3", [52, 52], [
            ("a", 0, 0, 6), ("b", 0, 6, 6), ("c", 1, 0, 12)]),
        "Figure_S4_covariate_coverage": o("figS4", [46], [("a", 0, 0, 12)]),
        "Figure_S5_copy_number_detail": o("figS5", [52, 46, 50], [
            ("a", 0, 0, 7), ("b", 0, 7, 5), ("c", 1, 0, 5), ("d", 1, 5, 7), ("e", 2, 0, 12)]),
        # one full-width embedding, row height set by the embedding's own aspect
        "Figure_S6_umap_by_sample": o("figS6", [round(L.WIDTH_MM * umap_aspect, 1)],
                                      [("a", 0, 0, 12)]),
        "Figure_S7_caf_like_state": o("figS7", [52, 48], [
            ("a", 0, 0, 5), ("b", 0, 5, 7), ("c", 1, 0, 7), ("d", 1, 7, 5)]),
    }


def build(cfg, atlas: Path, cache: Path) -> int:
    T = _p(cfg, "tables")
    out = _p(cfg, "figures")
    panels = _p(cfg, "figure_panels")
    out.mkdir(parents=True, exist_ok=True)
    panels.mkdir(parents=True, exist_ok=True)

    with provenance.start_run("14_figures", cfg, [atlas, cache]) as run:
        L.setup(cfg)
        A = ad.read_h5ad(atlas, backed="r")
        cols = ["study", "sample", "patient", "cell_type", "material_type", "treatment_timing",
                "disease_state", "assay", "cnv_score", "n_genes_by_counts", "total_counts",
                "pct_counts_mt", "pct_counts_spikein", "malignant_frozen", "malignant_evidence",
                "cnv_corr", "cnv_signal", "leiden"]
        d = A.obs[cols].copy()
        d["umap1"] = A.obsm["X_umap"][:, 0]
        d["umap2"] = A.obsm["X_umap"][:, 1]
        n_genes = A.n_vars
        A.file.close()

        emb = prepare_malignant(cfg, atlas, cache)
        programs = config.get(cfg, "figures.malignant_programs")
        dot = F.marker_table(atlas)
        bench = pd.read_csv(T / "v5_embedding_benchmark.csv", index_col=0)
        mix = pd.read_csv(T / "v5_mixing_by_celltype.csv")
        auc = pd.read_csv(T / "v5_cnv_auc_per_sample.csv")
        att = pd.read_csv(T / "qc_attrition.csv")
        prof = pd.read_csv(T / "scpca_library_profile.csv")
        cov = pd.read_csv(_p(cfg, "covariates"))
        scib = pd.read_csv(T / "v5_scib_results.csv", index_col=0)
        dep = pd.read_csv(T / "v5_depositor_labels.csv", index_col=0).reindex(d.index)
        lab = dep["depositor_label"].astype(str)
        conf = pd.crosstab(d.cell_type[lab != "unlabelled"], lab[lab != "unlabelled"])
        hist = {k: pd.read_csv(T / f"v5_elbo_{k}.csv", index_col=0)
                for k in ("sample_train", "sample_val", "study_train", "study_val")}
        loso = pd.read_csv(T / "f6_loso_mapping.csv")
        loso_ct = pd.read_csv(T / "f6_loso_by_celltype.csv")
        bulk = F.bulk_table(T / "f6_bulk_fractions.csv", _p(cfg, "bulk_clinical"))
        IM = ad.read_h5ad(_p(cfg, "immune_subset"))
        imm_mask = d.index.isin(IM.obs_names)
        xy_all = d[["umap1", "umap2"]].to_numpy()
        floors = config.get(cfg, "qc.gene_floor_by_assay")
        caps = config.get(cfg, "qc.mito_cap_pct_by_assay")

        coh = (d.groupby("study", observed=True)
                 .agg(cells=("sample", "size"), samples=("sample", "nunique"),
                      patients=("patient", "nunique"),
                      assay=("assay", lambda x: x.astype(str).mode().iat[0])))
        n_imm_types = IM.obs["immune_subset"].astype(str).replace("Low_quality", np.nan).nunique()
        steps = [
            ("Public deposits", f"{d['study'].nunique()} studies, GEO + ScPCA\n"
                                f"{int(att['raw'].sum()):,} barcodes\nEWSR1::FLI1 only"),
            ("Uniform QC", "assay-aware gene floor\nand mitochondrial cap\n"
                           "doublets, failed libraries"),
            ("Shared gene space", f"{n_genes:,} genes\nannotated in every arm\nraw counts kept"),
            ("scVI integration", "batch = sample\nbenchmarked against\nHarmony and study-level"),
            ("Annotation", f"{d['cell_type'].nunique()} cell types\n{n_imm_types} immune subsets\n"
                           f"{int(d['malignant_frozen'].sum()):,} malignant cells"),
            ("Applications", "reference mapping\nsurface targets\nbulk deconvolution"),
        ]

        xr = d["umap1"].max() - d["umap1"].min()
        yr = d["umap2"].max() - d["umap2"].min()
        lay = outlines(float(yr / xr))

        km = {}

        def _km(ax, **kw):
            km["p"], km["n"], km["events"] = F.panel_km(ax, **kw)

        spec = {
            "Figure_1_atlas_construction": {
                "a": (F.panel_flowchart, dict(steps=steps)),
                "b": (F.panel_umap_hero, dict(d=d)),
                "c": (F.panel_umap_meta, dict(d=d)),
                "d": (F.panel_cohort_bars, dict(coh=coh)),
                "e": (F.panel_marker_dotplot, dict(dot=dot, rotate_groups=True))},
            "Figure_2_integration_and_annotation": {
                "a": (F.panel_scib_scatter, dict(bench=bench)),
                "b": (F.panel_knn_vs_expected, dict(mix=mix)),
                "c": (F.panel_depositor_confusion, dict(conf=conf)),
                "d": (F.panel_loso_grid, dict(by_ct=loso_ct, overall=loso,
                                              min_n=config.get(cfg, "figures.loso_min_heldout_cells")))},
            "Figure_3_composition_across_samples": {
                "a": (F.panel_sample_composition, dict(d=d)),
                "b": (F.panel_composition_by_material, dict(d=d)),
                "c": (F.panel_study_celltype_grid, dict(d=d))},
            "Figure_4_malignant_compartment": {
                "a": (F.panel_cnv_statistic_gain, dict(
                    old_csv=T / "v5_cnv_auc_per_sample.csv", new_csv=T / "v6_cnv_per_sample.csv",
                    names=["Genome-wide\n|CNV|", "Signal\nstatistic", "Consensus\ncorr."])),
                "b": (F.panel_malignant_umap, dict(emb=emb)),
                "c": (F.panel_program_multiples, dict(emb=emb, programs=programs)),
                "d": (F.panel_evidence_bar, dict(csv=T / "v6_malignant_evidence.csv")),
                "e": (F.panel_mp_occupancy, dict(csv=T / "v6_metaprogram_occupancy.csv")),
                "f": (F.panel_confounding_heatmap, dict(csv=T / "v6_falsification.csv"))},
            "Figure_5_immune_microenvironment": {
                "a": (F.panel_immune_umap_threaded, dict(adata=IM, xy_all=xy_all, mask=imm_mask)),
                "b": (F.panel_immune_markers, dict(adata=IM)),
                "c": (F.panel_immune_by_sample, dict(
                    adata=IM, min_cells=config.get(cfg, "figures.immune_min_cells_per_sample")))},
            "Figure_6_atlas_as_reference": {
                "a": (F.panel_target_landscape, dict(csv=T / "f6_target_landscape.csv")),
                "b": (F.panel_target_margin, dict(csv=T / "f6_target_margins.csv")),
                "c": (F.panel_bulk_fractions, dict(D=bulk)),
                "d": (_km, dict(D=bulk)),
                "e": (F.panel_forest, dict(csv=T / "f6_bulk_survival.csv"))},
            "Figure_S1_quality_control": {
                "a": (F.panel_qc_attrition, dict(att=att)),
                "b": (F.panel_genes_by_assay, dict(d=d, floors=floors)),
                "c": (F.panel_mito_by_assay, dict(d=d, caps=caps))},
            "Figure_S2_deposit_specific_qc": {
                "a": (F.panel_library_gate, dict(prof=prof)),
                "b": (F.panel_spikein_by_plate, dict(d=d))},
            "Figure_S3_integration_detail": {
                "a": (F.panel_training_curves, dict(hist=hist)),
                "b": (F.panel_scib_heatmap, dict(scib=scib)),
                "c": (F.panel_umap_triple, dict(d=d, horizontal=True))},
            "Figure_S4_covariate_coverage": {
                "a": (F.panel_covariate_coverage, dict(cov=cov))},
            "Figure_S5_copy_number_detail": {
                "a": (F.panel_cnv_by_celltype, dict(d=d)),
                "b": (F.panel_cnv_auc, dict(auc=auc)),
                "c": (F.panel_cnv_vs_depth, dict(d=d)),
                "d": (F.panel_boundary_disagreement, dict(conf=conf)),
                "e": (F.panel_arm_scores, dict(csv=T / "v6_cnv_arm_auc.csv"))},
            "Figure_S6_umap_by_sample": {
                "a": (F.panel_umap_by_sample, dict(d=d))},
            "Figure_S7_caf_like_state": {
                "a": (F.panel_caf_specificity, dict(csv=T / "f4_specificity_cells.csv")),
                "b": (F.panel_caf_per_sample, dict(
                    csv=T / "f4_split_genome_specificity.csv",
                    csv_ds=T / "f4_split_genome_specificity_ds1500.csv",
                    prev_csv=T / "f4_prevalence.csv")),
                "c": (F.panel_caf_expression, dict(csv=T / "f4_group_expression.csv")),
                "d": (F.panel_caf_prevalence, dict(csv=T / "f4_prevalence.csv"))},
        }

        # The output directory holds exactly this set: drop any figure file a
        # renamed or removed figure left behind, so numbers never collide.
        keep = {f"{s}.{ext}" for s in spec for ext in ("png", "pdf")}
        for f in out.glob("Figure_*"):
            if f.name not in keep:
                f.unlink()

        made = []
        for stem, panels_spec in spec.items():
            o = lay[stem]
            paths = {k: L.render(o, k, fn, panels, **kw) for k, (fn, kw) in panels_spec.items()}
            made += [L.compose(o, paths, out / f"{stem}.png").name,
                     L.compose_pdf(o, paths, out / f"{stem}.pdf").name]
            print(f"{stem}  {L.WIDTH_MM:.0f} x {L.height_mm(o):.0f} mm", flush=True)

        stats = {"cells": len(d), "studies": int(d["study"].nunique()),
                 "samples": int(d["sample"].nunique()), "patients": int(d["patient"].nunique()),
                 "genes": int(n_genes), "raw_barcodes": int(att["raw"].sum()),
                 "malignant_cells": int(d["malignant_frozen"].sum()),
                 "immune_cells": int(IM.n_obs), "immune_subsets": int(n_imm_types),
                 "km_os_malignant": km}
        (out / "figure_stats.json").write_text(json.dumps(stats, indent=1) + "\n")
        run.add(figures=made, stats=stats)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
