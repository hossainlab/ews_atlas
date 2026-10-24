# EWSAtlas v5.0 — release notes

`data/processed/EWSAtlas_v5.h5ad` | **138,353 cells × 16,622 genes** |
1,899,286,803 bytes (1811 MiB), gzip |
6 studies · 65 samples · 38 patients

Every number here was read back from the object or from a stage manifest in
`logs/`. Build provenance: `src/pipeline/01_ingest_qc.py` →
`02_build_atlas.py` → `03_integrate.py` → `04_annotate.py` → `05_cnv.py`, all
parameterized from `src/config.json`.

## What changed from v4

| | v4 | v5 |
|---|---|---|
| studies | 5 | **6** (adds the Waltner snRNA arm) |
| cells | 64,702 | **138,353** |
| genes | 17,017 | 16,622 |
| patients | 33 | **38** |
| O-PDX studies | 1 (confounded) | **2, with 4 matched tumour/PDX patients** |
| scVI epochs | 120, early stopping never fired | **340, early-stopped** |
| HVGs | 2,000 | **5,000** |
| primary embedding | declared by construction | **declared by benchmark** |
| raw counts present | **no** (X was overwritten by lognorm) | yes, `layers['counts']` |

## Contents

| slot | contents |
|---|---|
| `X` | raw counts (integral; verified by test) |
| `layers['counts']` | the same counts, named |
| `layers['lognorm']` | library-size normalized (1e4) + log1p |
| `obsm['X_scvi']` | 30-d scVI latent — **the primary embedding** |
| `obsm['X_umap']` | UMAP on `X_scvi` |
| `obs` (55 cols) | controlled vocabulary + QC + `cell_type` + `cnv_score` |
| `uns['gene_space']` | feature-selection report |
| `uns['config']` | the full resolved config for this build |

Controlled vocabulary (enforced by `tests/test_invariants.py`): `sex` ∈ male /
female / not_reported; `material_type` ∈ patient_tumor / O-PDX / CTC / PBMC /
normal_tissue; `disease_state` ∈ primary / metastasis / recurrence / normal /
not_reported; `treatment_timing` ∈ treatment_naive / post_treatment /
not_reported; `sort_strategy` ∈ CD45pos_FACS / CD45neg_FACS / live_FACS /
microfiltration_CellSieve_8um / not_reported; `fusion_partner` ∈ FLI1 / ERG /
other / not_reported.

## Composition

Malignant (Ewing_sarcoma + Ewing_proliferating) 93,805 cells (67.8%), all 38
patients, all 6 studies. T_cell 16,025 · Monocyte_Macrophage 10,081 ·
Fibroblast_CAF 6,063 · Endothelial 4,109 · Unassigned 2,809 · Urothelial 2,453 ·
B_cell 1,679 · DC 715 · Plasma 396 · Platelet 218.

Covariate coverage: patient / material_type / assay 100%; disease_state 89.5%;
treatment_timing 84.6%; fusion_partner 82.9%; sex and age 44%.

## QC attrition

163,578 barcodes in → 138,353 cells (84.6%). Per study in
`results/tables/qc_attrition.csv`. Five samples were excluded whole by the
library gate (four Visser plates that are 66–96% ERCC spike-in, and ScPCA
library SCPCL001111 at median 44 genes).

## Validated

- **93.3% malignant agreement** with the depositors' independent annotations on
  the 69,400 cells they label.
- Composition sanity: Aynaud O-PDX 90.8% malignant, Cillo 82.5% immune, 2,433 of
  2,453 urothelial cells inside the bladder study.
- 26 invariant tests pass against this object.

## Known limitations — read before using

1. **Copy-number calling does not work as implemented** (AUC 0.677 overall,
   per-sample median 0.540, with stromal cells scoring above tumour). `cnv_score`
   is present on the object but **must not be used to call malignancy**. See
   `docs/05_annotation.md`.
2. **844 cells are `Fibroblast_CAF` here and `Malignant` to the depositors.**
   This is the fusion-low/mesenchymal boundary and it is unresolved.
3. **No MSC/osteoblast compartment is asserted** — the candidate clusters
   expressed RUNX2/SP7/IBSP at 0.02–0.05.
4. **Residual study structure in the immune compartments**: same-study kNN
   0.635 for T cells against 0.225 expected given composition.
5. **The treatment contrast is partly confounded with study** (Cramér's V 0.730
   in patient tumours); only Cillo and Visser contain both arms (19,336 cells).
6. **No survival data** in the atlas; outcome work is in the bulk cohorts and is
   not part of this object.
7. **No external query cohort exists** — Hong 2022 and Chen 2025 deposit no
   accession anywhere, so reference mapping can show projection stability but
   not generalization.
8. `sex` and `age_years` are missing for 56% of samples, mostly Aynaud, Goodspeed
   and Visser.

---

## v6 additions (this release)

New `obs` columns on the released object:

| column | meaning |
|---|---|
| `cnv_signal` | mean &#124;CNV&#124; per cell (the v5 quantity, retained) |
| `cnv_corr` | correlation to a per-sample consensus malignant profile — the statistic that works |
| `cnv_chr1p … cnv_chr17p` | signed arm-level scores for nine arms |
| `cnv_class` | `aneuploid` / `diploid` / `not_run`, thresholded per sample at the 95th percentile of that sample's reference cells |
| `malignant_frozen` | the frozen malignant set, 97,586 cells |
| `malignant_evidence` | how each cell earned its label: `marker_and_cna`, `marker_only_cna_negative`, `marker_only_cna_unassessed`, `cna_only` |

**Copy-number inference is now usable, and the fix was the statistic.** The v5
genome-wide mean of |CNV| gives a per-sample median malignant-vs-reference AUC
of 0.540. Adopting NBAtlas's window of 250 genes and a wider reference raises it
only to 0.597. Correlating each cell with a per-sample consensus profile reaches
0.712, with 23 of 45 scored samples above 0.7 and three still below 0.5
(PDX-184 0.435, SCPCS000490 0.452, CR-hPB120 0.459).

**Three things a user of this object must know.**

1. `cnv_class` is deliberately insensitive — only about 15% of marker-malignant
   cells clear the threshold. It is a rescue instrument, not a classifier.
2. `malignant_frozen` includes 3,781 cells rescued on copy-number evidence
   alone, all originally `Fibroblast_CAF` or `Unassigned`. The metaprogram MP24
   is carried almost entirely by these cells and is **not** a tumour state.
3. The `DC` label in `cell_type` is wrong: all 715 such cells are mast cells.
   The corrected immune labels are in `results/tables/v6_immune_labels.csv`
   (thirteen subsets, including 1,267 NK cells that sit inside `T_cell`).

**Negative results that ship with the object.** The ESCLA-derived fusion-activity
axis is a proliferation axis (Spearman 0.43 with S phase and 0.47 with G2M,
against 0.20 with NR0B1 and 0.23 with NKX2-2). Of the 28 recurrent expression
programs, 15 are confined to a single study, and those that look like hypoxia,
mesenchymal transition and cell cycle correlate with the
single-nucleus-versus-single-cell signature at 0.78, 0.69 and 0.64. Treating any
of these as biology without that correction would be a mistake.
