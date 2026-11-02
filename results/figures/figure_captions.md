# EWSAtlas figure captions

Produced by `src/pipeline/14_figures.py`. All figures are 180 mm wide,
written as 300 dpi PNG and vector PDF (editable text; dense scatter layers
rasterized at 300 dpi). Panel letters are bold uppercase; no panel carries a
title. n values are cells unless stated otherwise. Numbers below are read from
the released object or `results/tables/` at render time
(`figure_stats.json` holds the totals).

The sequence follows the structure shared by published harmonized atlases
(GBmap, the integrated breast-cancer atlas, HLCA, LuCA, NBAtlas):
construction → validation → composition → compartments → reference use →
clinical association.

---

## Figure 1. A harmonized single-cell reference for Ewing sarcoma from six studies and 38 patients

**(A)** Study design. 163,578 barcodes from six public deposits (GEO and ScPCA,
EWSR1::FLI1 Ewing sarcoma only) pass one QC policy, are placed in a shared
16,622-gene space, integrated with scVI using the sample as the batch, and
annotated into 15 cell types, 13 immune subsets and 97,586 malignant cells.
**(B)** The integrated embedding of 138,353 cells. Spatially coherent classes
are named on the map with their cell counts; classes dispersed through other
populations are keyed beside it rather than labelled at a misleading median.
**(C)** The same embedding coloured by study and by assay.
**(D)** Cells retained per study, coloured by assay (key in C), with patients
and samples per study.
**(E)** Mean log-normalized expression and fraction of expressing cells for the
marker panel, grouped by the compartment each gene marks. Every class carries
its own positive evidence.
See also Figures S1, S2 and S6.

## Figure 2. Integration removes study structure without erasing biology, and labels transfer to unseen studies

**(A)** scIB aggregate batch-correction against biological-conservation scores
for four candidate embeddings. scVI with `batch_key='sample'` leads on both
axes (total 0.606).
**(B)** Observed share of 15 nearest neighbours drawn from the same study,
against the share expected from cell-type composition alone.
**(C)** Agreement with depositor annotations, row-normalized percentage.
**(D)** Leave-one-study-out label transfer: each study is withheld and its
cells are labelled using a reference built from the other five.
Recall per cell type is shown as a percentage; overall agreement per study is
given beside its name (62–85%). Cell types with fewer than 20 held-out cells
are hatched, not coloured.
See also Figure S3.

## Figure 3. Composition across 65 samples varies most with source material

**(A)** Cell-type composition of every sample, grouped by study and ordered
within study by malignant share. Grey bars above give cells per sample (log
scale); strips below give study, material, disease state and treatment
timing. Not-reported covariates are drawn in light grey, never omitted.
**(B)** Mean per-sample composition by source material.
**(C)** Contribution of each study to each cell type. Dot size is the number
of cells (log scale); colour is the study's share of that type's cells.
See also Figure S4.

## Figure 4. The malignant compartment: a usable copy-number statistic, and recurrent programs that mostly track capture

**(A)** Malignant-vs-reference AUC per sample for three copy-number scoring
statistics: the genome-wide mean |CNV| of the first release (median 0.540), the
signal statistic (0.597) and correlation to a per-sample consensus malignant
profile (0.712).
**(B)** The 97,586 frozen malignant cells on their own scVI neighbour graph,
coloured by study.
**(C)** Scores of four recurrent programs on the same embedding (colour scale
clipped to the 1st–99th percentile of each score).
**(D)** How each malignant cell earned its label: markers and copy number,
markers alone, or copy number alone (3,781 rescued cells).
**(E)** Samples carrying each recurrent metaprogram; the number at the bar end
is the number of studies, and colour marks programs confined to one study.
**(F)** Spearman correlation of program scores with capture and technical
covariates.
See also Figures S5 and S7.

## Figure 5. The Ewing immune microenvironment resolves into thirteen subsets

**(A)** 29,114 immune cells on their own embedding, annotated as thirteen
subsets plus a low-quality class; the inset locates them in the atlas.
Lymphoid subsets are drawn in blues and purples and myeloid subsets in greens,
so the panel reads as a zoom into Figure 1.
**(B)** Marker evidence for each subset.
**(C)** Immune-subset share for the 49 samples with at least 50 immune cells
(median 185 per sample), grouped by study, with disease state and material
beneath.

## Figure 6. The atlas as a reference: surface targets and bulk-tumour composition associated with survival

**(A)** Detection of five surface targets per cell type, pooled across the
atlas.
**(B)** Within-sample selectivity: malignant minus normal detection, in
percentage points; the number of samples per comparison is under each column.
**(C)** Atlas-derived malignant, stromal and immune fractions deconvolved in
five bulk cohorts; bars mark cohort medians. Fractions are compared within
cohort, never across.
**(D)** Overall survival of 129 patients (51 deaths) split at the
within-cohort median of the atlas-derived malignant fraction; log-rank test
stratified by cohort, p = 0.047. The split is for display; the inference is
the continuous model in (E).
**(E)** Cox hazard ratios per standard deviation of each atlas-derived
fraction, stratified by cohort, for overall and event-free survival.

---

## Supplementary figures

**Figure S1. Barcode fate and assay-aware thresholds.** (A) Per-study barcode
fate through QC. (B) Detected genes per retained cell by assay, with each
assay's floor in the key. (C) Mitochondrial content by assay, with each assay's
cap in the key.

**Figure S2. Deposit-specific quality control.** (A) Median detected genes per
ScPCA nucleus library against the assay floor; one failed library is excluded
whole. (B) ERCC spike-in share per CEL-seq2 plate.

**Figure S3. Integration detail.** (A) scVI training and validation ELBO for
both batch-key choices. (B) Individual scIB metrics per embedding. (C) The
integrated embedding coloured by study, assay and sample.

**Figure S4. Covariate completeness.** Percentage of samples per study that
report each harmonized covariate.

**Figure S5. Copy-number detail.** (A) The superseded genome-wide score by cell
type. (B) Its per-sample AUC. (C) Its dependence on sequencing depth.
(D) Cells labelled differently by this atlas and by the depositors. (E)
Per-sample AUC for nine chromosome arms; arm-level scoring is not what made
the caller usable.

**Figure S6. The integrated embedding coloured by sample.**

**Figure S7. A CD73-high fibroblast-like malignant state.** (A) Own-patient
copy-number specificity by group. (B) The same per sample, including
depth-matched estimates. (C) Program and marker scores by group. (D)
Prevalence of the state per sample.
