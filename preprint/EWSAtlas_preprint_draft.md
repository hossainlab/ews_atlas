# EWSAtlas: a harmonized single-cell reference for Ewing sarcoma, and what it shows about recurrent tumour programs

**Draft v2.** Numbers in this draft were read back from the saved result tables,
not from notes.

---

## Summary

Ewing sarcoma is a fusion-driven paediatric bone and soft-tissue sarcoma, and it
is the conspicuous gap among sarcomas that already have an integrated single-cell
reference. We assembled EWSAtlas from all six public single-cell and
single-nucleus studies of Ewing sarcoma — 138,353 cells, 65 samples, 38 patients
— under one quality-control policy, one gene space, one controlled vocabulary and
one annotation pipeline, and defined the malignant compartment from copy number
rather than from the markers that name it. Two results follow that bear on how
such atlases are read. First, batch correction at the level of samples
outperforms correction at the level of studies on both the batch and the
biological axis, refuting our own pre-registered expectation. Second, and more
consequentially, most expression programs that recur across these six cohorts
track how the cells were captured rather than what they were doing: 15 of 28
recovered metaprograms are confined to a single study, and the three most
biologically appealing of them correlate with a single-nucleus-versus-single-cell
signature at 0.78, 0.70 and 0.64. A fusion-activity axis built from an external
knockdown experiment proves, on removal of cell-cycle genes, to track
proliferation as strongly as it tracks its own canonical targets. We therefore
present EWSAtlas as a reference rather than as a catalogue of new tumour states:
it transfers cell-type labels to a study it never saw at a median agreement of
0.82, it recovers tumour content from 343 archived bulk tumours well enough
to reproduce that cohort's published tumour–stroma survival association, and it
ranks five clinically pursued surface targets by a selectivity measure that
single cohorts cannot compute.

## Introduction

*(to be written: EWSR1::FLI1 biology; why an atlas; what NBAtlas established as a
template; what the six studies individually reported; the gap.)*

## Results

### A harmonized reference across six studies and 38 patients

Six public studies met inclusion (EWSR1::FLI1 or EWSR1::ERG-driven Ewing sarcoma;
DSRCT, clear cell sarcoma, FET–CREB tumours and myxoid liposarcoma excluded;
cell lines excluded): Aynaud 2020 (14,504 cells, 5 samples, 5 patients), Cillo
2022 (14,353 / 5 / 4), Goodspeed 2025 (12,639 / 15 / 11), Mao 2024 (15,394 / 4 /
2), Visser 2023 (6,575 / 24 / 11) and Waltner 2025 (74,888 / 12 / 5). Quality
control retained 138,353 of 163,578 barcodes (84.6%) under assay-aware gene
floors and mitochondrial caps, with a library-level gate that removes whole
libraries whose median complexity falls below their assay floor — a failure a
per-cell rule cannot see, because a uniformly failed library looks internally
consistent.

Sample count and cell count rank the studies in nearly opposite orders (Visser
contributes 24 samples and 6,575 cells; Waltner 12 samples and 74,888), which is
the practical argument for batching on sample rather than study (Figure 1A).
Fifteen cell types were assigned at the level of 34 Leiden communities with
obligatory-marker guards rather than per cell, leaving 2,809 cells (2.0%)
`Unassigned` rather than forcing a nearest label.

Two cautions belong in the reading of every figure that follows. Waltner supplies
54% of all cells and is the only single-nucleus study, so assay and cohort are
perfectly confounded for the nucleus arm. And the atlas is built entirely from
public data — which is the precedent set by NBAtlas, not a weakness, but it means
no external query cohort exists for reference-mapping validation.

### Batch correction at the level of samples outperforms correction at the level of studies

Four embeddings were scored with scIB metrics against an unintegrated PCA
baseline. scVI with `batch_key='sample'` leads on both axes (batch correction
0.611, biological conservation 0.603, total 0.606), ahead of scVI with
`batch_key='study'` (0.532), Harmony (0.483) and unintegrated PCA (0.437)
(Figure 2A). We had expected 65 sample-level batches to spend parameters on too
few cells each and over-correct; they do not.

Apparent batch structure must be read against composition. Observed same-study
15-nearest-neighbour share is 0.760 against 0.338 expected under random mixing,
but the studies sample different compartments, so roughly half of that residual
is composition rather than batch (Figure 2B).

Annotation agreement with the labels the depositors themselves released is 93.0%
over the 69,400 cells that carry them, and the disagreement is concentrated in
one row: cells deposited as `Stromal` are 76.1% malignant and 23.0% stromal under
our labelling, and cells deposited as `Epithelial` are 98.2% malignant
(Figure 2D). This boundary — malignant versus fibroblast — recurs throughout the
paper and is the subject of the next section.

### A correlation statistic, not a parameter change, makes copy-number inference usable

Defining the malignant compartment from markers alone is circular when the
question is which cells are malignant. We therefore scored copy number per study
with infercnvpy and classified cells per sample. Our first implementation — the
genome-wide mean of |CNV| — separated malignant from reference cells at a median
per-sample AUC of 0.540, which is not usable. Neither a wider window nor a wider
reference set fixed it (0.597). What fixed it was the statistic: correlating each
cell's profile with a per-sample consensus malignant profile, seeded from
high-signal non-reference cells, reaches a median AUC of 0.712, with 23 of 45
scored samples above 0.7 (Figure 3A). Three samples remain below 0.5 (PDX-184
0.435, SCPCS000490 0.452, CR-hPB120 0.459) and are reported as such.

Arm-level scoring, which we expected to be the fix for a disease with few focal
changes, is not: across all 45 samples chr8q is the strongest arm at a median AUC
of 0.574, chr12q reaches 0.526 and chr16q 0.492 (Figure 3B). The favourable
figures we first obtained (0.61–0.72) came from pooling two immune-rich studies
at study level and do not survive per-sample scoring.

The frozen malignant set is 97,586 cells, each carrying the evidence that put it
there: 12,114 with both markers and copy number, 67,587 with markers where copy
number was negative, 14,104 where copy number could not be assessed, and 3,781
admitted on copy number alone (Figure 3C). The copy-number threshold is the 95th
percentile of each sample's own reference cells and is used only to rescue cells
the markers missed, never to remove cells the markers found; 1,634 aneuploid
calls on immune and endothelial cells were discarded as false positives.

### Most recurrent expression programs track capture, not biology

Per-sample non-negative matrix factorisation across the 45 samples with at least
100 malignant cells (ranks 4–9, relative expression per Gavish et al.) yielded
1,755 programs, 864 robust, clustering into 28 metaprograms. **Fifteen of the 28
are confined to a single study** (Figure 4A); in NBAtlas, none were.

Testing the candidates against the technical variables that could explain them
(Figure 4B) shows the problem is not marginal. The programs that read as hypoxia,
mesenchymal transition and cell cycle correlate with a single-nucleus-versus-
single-cell signature at Spearman 0.781, 0.695 and 0.642, and the hypoxia-like
program correlates with detected genes at 0.524. Doublet score explains none of
them (|rho| <= 0.183).

The same test dissolved our leading biological hypothesis. An EWSR1::FLI1 activity
axis derived from an external knockdown experiment (ESCLA) — chosen precisely to
avoid defining fusion activity from atlas markers — correlates with S phase at
0.430 and G2/M at 0.468 against 0.200 for NR0B1 and 0.226 for NKX2-2. After
removing 133 cycle-associated genes it still tracks cycle (0.246, 0.256) as
strongly as its own canonical targets (0.229, 0.255). The mechanism is
unremarkable in hindsight: knocking down the fusion arrests growth in the source
experiment, so the signature carries proliferation by construction.

**Attribution.** Truong et al. (2023) reported warm-dissociation and gene-length
artifacts in 125,831 cells and nuclei from Ewing sarcoma, DSRCT and osteosarcoma
patient-derived xenografts, and showed they are computationally correctable. Our
contribution is the same phenomenon in patient tumours across six independent
cohorts, with the metaprogram decomposition that shows how it is mistaken for
biology. Because Waltner is our only single-nucleus study, assay and cohort
cannot be separated within this atlas; what supports assay as the cause is
mechanistic — the single-study programs are built from long, intron-rich genes —
together with Truong's matched cell/nucleus pairs, in which study is held
constant by design.

### The immune compartment resolves into thirteen subsets and corrects two atlas-level calls

Re-clustering 29,114 immune cells on their own neighbour graph resolves thirteen
subsets (T_CD4 7,532; T_cytotoxic 5,860; Monocyte 4,974; Macrophage 3,351; B_cell
1,681; NK 1,267; cDC2 1,137; Mast 715; Treg 670; Macrophage_LAM 532; pDC 307;
Platelet 218; Plasma 88) plus 782 low-quality cells (Figure 5).

The subset view corrected four atlas-level calls that a coarse pass cannot avoid,
and the released object has been revised accordingly: all 715 cells first labelled
`DC` are mast cells; 1,267 NK cells sat inside `T_cell`; the real conventional
dendritic cells (1,137 cDC2) sat inside `Monocyte_Macrophage`, so the atlas had no
dendritic-cell class at all; and 307 plasmacytoid dendritic cells sat inside
`Plasma`. Twelve classes became fifteen. The copy-number reference set is
unchanged by the revision, so the malignant calls above still stand.

No subset proportion differs by treatment arm after correction for multiple
testing: the smallest p across 13 subsets in 24 patients is 0.0339 (NK), and its
study-adjusted coefficient has the opposite sign to the raw group means (0.053
treatment-naive against 0.013 post-treatment) — a study-level reversal rather
than a biological shift.

### The atlas as a reference: target selectivity, label transfer and bulk deconvolution

*Surface targets.* STEAP1, IL1RAP, CD276 (B7-H3), CD99 and NT5E (CD73) all carry
active CAR-T, CAR-NK or antibody–drug-conjugate programmes in this disease. What
no single study can supply is the denominator. Pooled across the atlas, the
targets are detected in 56.6% (CD99), 28.1% (IL1RAP), 22.0% (STEAP1), 9.2%
(CD276) and 1.4% (NT5E) of malignant cells (Figure 6A). Because detection depends
on assay and the studies do not contribute cell types evenly, selectivity is
computed as a paired within-sample contrast against that sample's own normal
cells (Figure 6B). Across the nine non-urothelial microenvironment classes,
STEAP1 (+1.9 to +13.4 percentage points), IL1RAP (+2.9 to +16.7) and CD99 (+0.5
to +24.2) are positive against every class; CD276 is negative against fibroblasts
(−4.7); NT5E is positive against four classes, flat against two and negative
against three, including −15.0 against endothelium.

Margin and off-tumour burden are separate axes, and CD99 is the case that shows
why: its within-sample margin is positive everywhere, yet in absolute terms it is
carried by 74.2% of urothelial, 72.9% of fibroblast, 70.3% of endothelial and
68.3% of platelet cells against 56.1% of Ewing sarcoma cells — more normal cells
carry it than tumour cells do. Read together, the two axes favour IL1RAP and
STEAP1; both, however, are strongly negative against urothelium (−21.7, −20.9),
consistent with their documented normal expression and a reminder that this atlas
samples the tumour's neighbours, not every tissue a therapy will reach.

*Label transfer.* Training a scVI reference on five studies and mapping the
sixth by scArches surgery recovers the full-atlas annotation at a median
agreement of 0.819 across the six folds — 0.853 (Cillo), 0.844 (Visser), 0.822
(Waltner), 0.816 (Goodspeed), 0.787 (Aynaud) and 0.625 (Mao) — rising to
0.688–0.883 among the 81–95% of query cells where at least 60% of neighbours
agree (Figure 6C). Mao's value is structural rather than technical: its 2,433
urothelial cells have no counterpart in a reference built from the other five
studies, and excluding urothelial and unassigned cells it reaches 0.810, with a
median of 0.822 across folds. The four classes added by the annotation revision
transfer well — mast 0.97, NK 0.80, plasmacytoid dendritic 0.75, conventional
dendritic 0.71 by median recall — independent support for the revision.
The instructive failures are structural: a
cell type contributed by a single study cannot be transferred when that study is
the query, and the Ewing_sarcoma/Ewing_proliferating boundary is a cycle-phase
split rather than a discrete type. Because the reference labels were produced with
the query study present, this measures projection stability under study removal
and not generalization to an external cohort; no external Ewing sarcoma
single-cell cohort exists in any public archive we searched.

*Bulk deconvolution.* Using the atlas as a cell-type reference for 343 archived
bulk tumours — the primary record of each patient, de-duplicated from 408 GEO
expression records across six series — recovers mean malignant fractions of 0.73–0.83 across five cohorts.
In the 129 patients with survival, higher tumour content tracks better overall
survival (HR 0.76 per standard deviation, 95% CI 0.59–0.97) and higher stromal
content worse (1.32, 1.03–1.70) (Figure 6D). Neither survives Bonferroni
correction for the eight tests shown, and neither is new: this is the
tumour–stroma association GSE63157 itself published, recovered here as evidence
that the deconvolution carries real signal rather than as a finding.

Three corrections were required before these fractions were credible, and each
produced a plausible wrong answer first. NNLS on log-transformed expression
returns ~0% malignant content, because a mixture is additive only on the linear
scale. With `Ewing_sarcoma` and `Ewing_proliferating` as separate reference
profiles — collinear, differing mainly by cell-cycle genes — NNLS placed 0.358 on
the proliferating profile and 0.041 on the class 8.4-fold larger in the atlas.
And a malignant signature built from the atlas as a whole is nucleus-dominated
(68,698 of 97,586 malignant cells are nuclei) and does not resemble a bulk tissue
extract: it returned mean malignant fractions of 0.33–0.53, against 0.73–0.83
once the signature was restricted to whole-cell libraries.

### What does not survive removing the largest study

Waltner contributes 54% of cells, so every per-sample headline was recomputed
over the 53 non-Waltner samples (63,465 cells, 33 patients, 28,888 malignant
cells; Figure S9). The copy-number statistic holds (median AUC 0.724 against
0.712); 20 of 28 metaprograms survive, 8 of them single-study; and the number of
normal classes against which each target carries a positive margin is unchanged
for STEAP1 (9 of 10), IL1RAP (9 of 10), CD276 (9 of 10) and NT5E (4 of 10). CD99
is the exception, falling from 10 of 10 to 8 of 10 — its apparent universality is
the one selectivity result that depends on the single-nucleus arm. The headline
that cannot be re-tested this way is the capture result itself, since dropping
Waltner removes that arm entirely.

## Discussion

*(to be written. Argument to make: the value of a disease atlas built from
public data is not that it finds new states but that it supplies denominators,
controls and a mapping surface; and that the metaprogram result is a warning for
every atlas that mixes dissociation protocols. Place against Truong 2023, against
the EwS single-cell literature's individual claims, and against NBAtlas as
template.)*

## Limitations of the study

One study supplies 54% of cells and is the only single-nucleus dataset, so assay
and cohort cannot be separated within this atlas. Reference mapping is validated
against held-out studies rather than an external cohort, because none exists
publicly. Treatment contrasts rest on 28 patients from two studies that carry
both arms; the xenograft contrast on four matched patients from one study, where
0.125 is the smallest attainable two-sided p; circulating-tumour-cell claims are
not supportable at all, since 405 of 6,941 CTC-library cells are malignant. The
CD73-high fibroblast-like state we tested is genetically supported in three Ewing
patients of 38, too few for a prevalence claim. Proportional-hazards assumptions
for the survival models were not assessed. Target detection is RNA rather than
protein, so positivity figures are lower bounds on surface presentation and a
negative selectivity margin is stronger evidence than a positive one.

## STAR Methods

*(outline — each stage is a script in `src/pipeline/` and a document in `docs/`)*

1. Study inclusion and exclusion — docs/01
2. Quality control, assay-aware thresholds, library gate — docs/02, `01_ingest_qc.py`
3. Gene space — docs/03, `02_build_atlas.py`
4. Integration and benchmark — docs/04, `03_integrate.py`, `03b_train_scvi.py`
5. Annotation and its revision — docs/05, docs/10, `04_annotate.py`, `11_annotation_revision.py`
6. Copy number and the malignant set — docs/06, `05_cnv.py`
7. Metaprograms and falsification — docs/06, `07_caf_state.py`
8. Reference mapping — docs/12, `08_reference_mapping.py`
9. Bulk cohorts, deconvolution, survival — docs/07, docs/12, `09_bulk_deconvolution.py`
10. Surface targets — docs/12, `10_surface_targets.py`
11. Sensitivity analysis — `12_waltner_sensitivity.py`
12. Deposition — docs/14, `13_cellxgene_export.py`

## Figures

Main 1–6 and supplementary S1–S9 as listed in `results/pub_figs/figure_captions.md`,
which carries the full legends. Figures are 180 mm at 300 dpi with no in-figure
titles.

## Data and code availability

The analysis object, a CELLxGENE-schema deposit object (138,353 x 16,286 genes,
Ensembl identifiers) and the trained scVI model for scArches mapping are to be
deposited on Zenodo; `deposit/zenodo.json` is prepared and needs creators and a
licence. All code is in this repository. Source data are the six public studies
cited above.

One caveat belongs with the model. It was trained by the pipeline stage that the
released embedding predates, and the two agree globally but not locally
(Procrustes cosine 0.776, Spearman correlation of pairwise distances 0.904, but
34.7% overlap of 15-nearest-neighbour sets). The released `X_scvi` is retained as
the embedding every figure is built on, and the model's own latent is stored
beside it as `X_scvi_model_v5`; the query-mapping utility transfers labels in the
latter, because a cell projected by the model must find its reference neighbours
in the space the model produces. Users re-deriving the embedding should expect a
similar neighbourhood overlap and should not read it as a defect.
