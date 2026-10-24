# EWSAtlas figure captions

Figures carry no titles; each caption headline states the claim the figure is
meant to earn, in the style of Cell Press. Panel letters are bold uppercase.
All figures are 180 mm wide at 300 dpi with editable vector text
(`pdf.fonttype 42`). n values are cells unless stated otherwise.

---

## Figure 1. A harmonized single-cell reference for Ewing sarcoma spans six studies and 38 patients

**(A)** Cohort schematic. One column per study with its assay glyph (droplet
whole-cell, droplet nucleus, plate well), sample count, patient count and cell
count; every value is computed from the released object at render time, so the
panel cannot drift from the atlas. Sample count and cell count rank the
studies
in nearly opposite orders — Visser contributes 24 samples and 6,575 cells,
Waltner 12 samples and 74,888 — which is the practical argument for batch
correction at the level of samples (Figure 2A).
**(B)** The scVI-integrated embedding of 138,353 cells, coloured by cell type
and
ringed by the composition of each type. One arc per cell type, all arcs of
equal
angular width: the arcs show composition, not abundance, since sizing them by
cell count would give the malignant compartment two thirds of the circle. Each
arc carries a solid outer band in the type's colour and, inside it, one band
per
metadata variable — study, then sample. A type drawn from one deposit reads as
a
solid band; a pooled type reads as fine stripes. The Ewing sarcoma arc is
coarsely striped because Waltner's twelve nucleus libraries supply most
malignant
cells, whereas T cell, monocyte/macrophage and fibroblast arcs are finely
striped
across studies and samples; urothelial is a single solid block, being Mao's
bladder control alone. Contours inside enclose the spatially coherent types
and
carry no text — the arc labels already name every type. Fifteen cell types
were
assigned at the level of 34 Leiden communities (resolution 1.0), not per cell,
and each label additionally requires one of its own obligatory markers above a
floor; 2,809 cells (2.0%) in communities with no positive marker evidence are
left
`Unassigned` rather than given a nearest label.
**(C)** Mean log-normalized expression and the fraction of expressing cells
for
the marker panel, derived from the released object at render time and grouped
by
the compartment each gene marks. Every one of the fifteen classes carries its
own
positive evidence, including the four the immune re-analysis separated: NK on
NKG7 and GNLY, conventional dendritic cells on CD1C and FCER1A, plasmacytoid
dendritic cells on LILRA4 and IRF7, and mast cells on TPSAB1 and TPSB2. UPK1B
is
absent from the 16,622-gene space.
See also Figures S1, S2, S4 and S6.

## Figure 2. Sample-level batch correction outperforms study-level correction on both axes

**(A)** scIB aggregate batch-correction against biological-conservation scores
for four candidate embeddings. scVI with `batch_key='sample'` leads on both
axes
(total 0.606), which refuted the pre-registered expectation that study-level
batches would conserve more biology.
**(B)** Observed share of 15 nearest neighbours drawn from the same study,
against the share expected from cell-type composition alone. The residual is
real for T cells (0.639 observed against 0.220 expected) and monocytes (0.760
against 0.361), and close to what composition alone predicts for urothelial
cells, mast cells and platelets.
**(C)** Same-study neighbour share per study for all four embeddings.
**(D)** Agreement with depositor annotations on the 69,400 cells that carry
them, as a row-normalized percentage. The 77%/18% split on the fibroblast row
is
the unresolved boundary quantified.
See also Figure S3.

## Figure 3. A correlation statistic, not a parameter change, makes copy-number inference usable

**(A)** Malignant-vs-reference AUC per sample for three scoring statistics.
The
v5 genome-wide mean of |CNV| has a median of 0.540; re-running with NBAtlas's
window of 250 genes and a wider reference raises it only to 0.597; correlating
each cell's profile with a per-sample consensus malignant profile (Tirosh et
al.) reaches 0.712, with 23 of 45 samples above 0.7; three samples remain
below 0.5 (PDX-184 0.435, SCPCS000490 0.452, CR-hPB120 0.459).
**(B)** Per-sample AUC for nine chromosome arms scored signed. Arm-level
scoring
is **not** what rescued the caller: pooling two immune-rich studies gives
0.61-0.72 for chr8q and chr12q, but across all 45 scored samples chr8q is the
strongest arm at a median of 0.574, chr12q is 0.526, and chr16q at 0.492 shows
no
loss signal.
**(C)** How each of the 97,586 frozen malignant cells earned its label. The
copy-number threshold is the 95th percentile of each sample's own reference
cells, which is deliberately insensitive: it is used only to rescue cells the
markers missed (3,781, all originally `Fibroblast_CAF` or `Unassigned`), never
to remove cells the markers found. 1,634 aneuploid calls on immune and
endothelial cells were discarded as false positives.
See also Figure S5.

## Figure 4. The Ewing immune microenvironment resolves into thirteen subsets, correcting four atlas-level calls

**(A)** 29,114 immune cells re-clustered into 24 communities and annotated as
thirteen subsets plus a low-quality class, each carrying its own positive
marker evidence. The genes that separate the three myeloid-adjacent classes
are CD1C and LILRA4, not FCER1A: conventional dendritic cells reach CD1C 0.95
and plasmacytoid dendritic cells LILRA4 1.65 and IRF7 2.01, against CD1C 0.00,
LILRA4 0.00 and IRF7 0.06 in mast cells, whereas FCER1A is 1.45 in
conventional dendritic cells and still 0.52 in mast cells, which carry the
high-affinity IgE receptor themselves. Four of these corrections are applied
to the released object, taking its `cell_type` vocabulary from twelve classes
to
fifteen: all 715 cells first called `DC` are mast cells, 1,267 NK cells sat
inside `T_cell`, the 1,137 conventional dendritic cells sat inside
`Monocyte_Macrophage` and 307 plasmacytoid dendritic cells inside `Plasma`.
The
copy-number reference set is unchanged by the revision, so the malignant calls
of
Figure 3 still hold; `cell_type_v5_original` preserves the pre-revision
labels.
The embedding is on a UMAP computed from the immune
subset's own neighbour graph; the inset locates these cells in the full atlas
(grey). The subset embedding is necessary rather than cosmetic: the atlas-wide
UMAP restricted to immune cells is a cut-out of Figure 1, in which the subsets
sit where they sat in a layout optimised for 138,353 cells. Clustering used
the
existing scVI latent rather than a fresh integration of the subset (declared
deviation from NBAtlas, which re-integrated).
**(B)** Marker evidence for each subset, and the source of the four
corrections
now applied to the released object: all 715 cells first labelled `DC` are mast
cells (mean log-normalized TPSAB1 3.80, TPSB2 4.46, KIT 2.62, against CD1C
0.00;
CPA3 is absent from the 16,622-gene space), 1,267 NK cells sat within
`T_cell`,
the 1,137 conventional dendritic cells sat within `Monocyte_Macrophage`, and
307
plasmacytoid dendritic cells sat within `Plasma`. No subset proportion differs
between treatment arms after correction for multiple testing (smallest p =
0.0339 for NK across 13 subsets in 24 patients; the study-adjusted coefficient
for NK has the opposite sign to its raw group means, a study-level reversal).

---

## Figure 5. Most recurrent expression programs in a six-study atlas track capture, not biology

**(A)** The 28 metaprograms recovered by per-sample NMF (ranks 4-9, relative
expression, 45 samples with at least 100 frozen malignant cells), ordered by
the
number of samples carrying them. 15 of 28 are confined to a single study — in
NBAtlas, none were.
**(B)** Spearman correlation of each candidate program with the technical
drivers that could explain it. The programs that read as hypoxia, mesenchymal
transition and cell cycle correlate with the single-nucleus-versus-single-cell
signature at 0.78, 0.69 and 0.64, and the hypoxia-like program with detected
genes at 0.52. Doublet score is not a driver for any candidate (|rho| <=
0.18).
**(C)** MP24, the candidate fusion-low program, scored by how each cell
entered
the frozen malignant set: it is carried by the copy-number-rescued fibroblasts
(mean 0.503) and is near zero in every marker-defined class, so it is an
artefact of the rescue step rather than a tumour state.
The underlying artifact is not new: Truong et al. (*BMC Cancer* 2023;23:1105)
reported warm-dissociation and gene-length bias in 125,831 cells and nuclei
from sarcoma patient-derived xenografts. What is shown here is its effect on
the recurrence analysis of an integrated patient atlas, where it survives
every within-study check. See also Figure S5 and Table v6_metaprograms.

## Figure 6. The atlas is usable as a reference: it ranks surface targets by selectivity, labels studies it never saw, and recovers tumour content from bulk

**(A)** Detection of five cell-surface targets with active CAR-T, CAR-NK or
antibody-drug-conjugate programmes in Ewing sarcoma, by annotated cell type,
pooled across 138,353 cells. Dot area is the per cent of cells positive, fill
is
mean log-normalized expression. In malignant cells the targets are detected in
56.6% (CD99), 28.1% (IL1RAP), 22.0% (STEAP1), 9.2% (CD276) and 1.4% (NT5E) of
cells. Detection is RNA, so these are lower bounds on surface presentation.
The `Urothelial` class is the bladder control of Mao et al. rather than a
tissue at risk in most Ewing sarcoma.
**(B)** Selectivity as a paired, within-sample contrast: the per cent of
malignant cells positive minus the per cent of that sample's normal cells of
the
given type, median over the 47 samples with at least 50 malignant and 20
normal
cells (sample count per column in parentheses). Pooled comparisons are
confounded — 54% of the atlas is Waltner nuclei at ~1,100 detected genes, and
cell types are not contributed evenly by the studies — so the margin is
computed
with assay, depth and patient held fixed. Across the nine non-urothelial
classes,
STEAP1 (+1.9 to +13.4 percentage points), IL1RAP (+2.9 to +16.7) and CD99
(+0.5
to +24.2) are positive against every class; CD276 is positive against eight
but
negative against fibroblasts (−4.7); NT5E is positive against four, flat
against
two and negative against three, including −15.0 against endothelium. Margin
and
off-tumour burden are separate axes and must be read together with (A): CD99's
margin is positive everywhere, yet it is detected in 74.2% of urothelial,
72.9%
of fibroblast, 70.3% of endothelial and 68.3% of platelet cells against 56.1%
of
Ewing sarcoma cells — more normal cells carry it than tumour cells do. Read on
both axes the data favour IL1RAP and STEAP1, but both are strongly negative
against urothelium (−21.7, −20.9), consistent with their documented normal
expression. Column support is uneven: T cell 45 samples, monocyte/macrophage
24,
fibroblast 15, endothelial 14, dendritic 12, NK 10, B cell 8, urothelial 4,
mast
3, plasmacytoid dendritic 2 — the last three cannot rank targets. The plasma-
cell
column disappears after the annotation revision, because 307 of the 396 cells
first called plasma cells are plasmacytoid dendritic cells and no sample
retains
20.
**(C)** Leave-one-study-out label transfer. For each study a scVI reference is
trained on the other five (capped at ~50,000 cells by trimming large study x
cell
type groups, highly variable genes selected on the reference only), the held-
out
study is mapped by scArches surgery, and labels are transferred by 15-nearest
neighbours in the shared latent space. Agreement with the full-atlas
annotation is
0.853 (Cillo), 0.844 (Visser), 0.822 (Waltner), 0.816 (Goodspeed), 0.787
(Aynaud)
and 0.625 (Mao); median 0.819. Restricting to the 81.2–95.1% of query cells
where
at least 60% of neighbours agree gives 0.688–0.883. Mao's value is structural
rather than technical: its 2,433 urothelial cells have no counterpart in the
reference (20 urothelial cells exist across the other five studies combined),
so
they are necessarily called malignant; excluding urothelial and unassigned
cells
Mao reaches 0.810 and the median across folds is 0.822. The four classes added
by
the annotation revision transfer well — median recall 0.97 (mast), 0.80 (NK),
0.75 (plasmacytoid dendritic), 0.71 (conventional dendritic) — which is
independent support for the revision. Two failures are informative rather than
noise: Aynaud's T cells are recovered at 1% and called Ewing sarcoma, and
`Ewing_proliferating` at a median of 32%, a cycle-phase split rather than a
discrete type. The reference labels were produced with the query study
present,
so this measures projection stability under study removal, not generalization
to
an external cohort; no external Ewing sarcoma single-cell cohort exists in any
archive searched.
**(D)** Atlas-derived composition against outcome in 129 unique patients
(GSE63157 n=85, GSE17674 n=44), hazard ratio per standard deviation from a Cox
model stratified by cohort. Deconvolution is NNLS in linear space against a
signature built from whole-cell libraries only, with the two malignant classes
entered as one profile; mean malignant fraction is 0.73–0.79 across the five
cohorts. Higher tumour content tracks better overall survival (HR 0.77,
0.61–0.98) and higher stromal content worse (1.33, 1.04–1.71). Neither
survives
Bonferroni correction for the eight tests shown (threshold 0.0063), and
neither
is a new finding: this is the tumour–stroma association GSE63157 itself
published, recovered here as a positive control that the deconvolution carries
real signal.
See also Figures 1C and 4, and Tables f6_target_margins, f6_loso_by_celltype
and
f6_bulk_survival.

## Supplementary figures

**Figure S1. Barcode fate and assay-aware thresholds.** (A) Per-study
attrition
from 163,578 barcodes to 138,353 retained cells (84.6%), by exclusion reason.
(B) Detected genes and (C) mitochondrial content per assay, with the
assay-specific floors (200/200/100) and caps (30/40/10%) marked. The criterion
is uniform (per-sample 5-MAD); only the floors and caps vary by assay.

**Figure S2. Deposit-specific quality control.** (A) Median detected genes per
nucleus for the 13 ScPCA libraries; SCPCL001111 fails the library-level gate
whole. (B) ERCC spike-in share per CEL-seq2 plate; four plates at 66-96%
spike-in were excluded before this step.

**Figure S3. Integration detail.** (A) scVI training and validation ELBO for
both batch definitions, each early-stopped. (B) The individual scIB metrics
behind Figure 2A. (C) The integrated embedding coloured by study, by assay and
by
sample. These three panels were moved here out of Figure 1, where study and
assay
duplicated a column header of its schematic, a glyph row and a ring band, and
the
sample panel duplicated Figure S6. The sample panel carries no legend because
65
levels cannot be read as one; it is shown to demonstrate that no single
library
occupies a region of its own.

**Figure S4. Covariate completeness.** Percentage of samples per study
reporting
each harmonized covariate. Fusion partner, sex and age are the gaps.

**Figure S5. The first copy-number attempt, retained as the record of what
failed.** (A) v5 genome-wide score by cell type, (B) per-sample AUC, (C) score
against sequencing depth, and (D) the label disagreements this left
unresolved.
The score swings 0.011 across depth bins while malignant and reference cells
differ by 0.003. Superseded by Figure 3.

**Figure S6. The integrated embedding coloured by sample.** All 138,353 cells,
one hue per sample, hues spaced evenly around the colour wheel and shuffled so
that neighbouring samples are unlike; draw order is randomised so no study
paints over another. No legend is shown because 65 levels cannot be read as
one. The panel's purpose is negative: no single library occupies a region of
the embedding of its own.

**Figure S7. A CD73-high fibroblast-like malignant state is genetically real
but present in only a minority of tumours.** Test of the hypothesis that the
3,781 cells rescued into the frozen malignant set on copy-number evidence
alone are the CD73+ CAF-like tumour population reported by Wrenn et al. (Clin
Cancer Res 2023;29:5140-5154), rather than contaminating stroma.  **(A)** Own-
patient copy-number specificity: each cell's correlation to its own sample's
malignant consensus minus its mean correlation to the consensus of the other
samples in the same study. Ambient RNA and technical noise raise correlation
to any consensus; only a shared clonal karyotype raises it to one patient in
particular. Cells were selected on one half of the genome (alternating
chromosomes) exactly as the pipeline selects them, and are scored here on the
other half, which played no part in selection. All cells are downsampled to
exactly 1,500 counts, because the candidate cells are five to seven times
deeper than reference cells and specificity rises with depth (76,475 of
138,353 cells reach that depth). Cell-level means: diploid reference 0.083,
fibroblast-like copy-number-negative 0.121, fibroblast-like copy-number-
positive 0.229, marker-defined malignant 0.296 — the candidates sit between
stroma and tumour, nearer the tumour. **(B)** The same three quantities per
sample, for the 8 of 45 scored samples carrying at least ten candidate cells;
open rings repeat the candidate value after depth matching, available for 4
samples. Two of the eight samples are the urothelial carcinoma control patient
of Mao et al., in which the marker-defined malignant cells are urothelial
carcinoma rather than Ewing sarcoma; the two plate-based Visser samples are
negative. **(C)** Transcriptional profile of the same four groups at native
depth. CD73 (NT5E) is detected in 21.8% of candidates against 9.9% of copy-
number-negative fibroblast-like cells, 6.7% of reference cells and 0.8% of
marker-defined malignant cells, reproducing the published marker. The
matrisome/CAF program (0.861 against 0.962) and the EWS::FLI1-repressed gene
set (0.261 against 0.237) do **not** separate candidates from stroma, because
the EWS::FLI1-repressed set is itself mesenchymal — transcription alone cannot
adjudicate this question. Doublets are excluded: candidates carry the lowest
doublet score of any group (0.029 against 0.092 in marker-defined malignant
cells). **(D)** Share of each sample's malignant compartment held by the
state, by study. It reaches 8-21% in four Mao samples, 8.6% in one Waltner
sample and 9.3% in one Cillo sample, and is under 1% or absent in 37 of 45
samples; counts above each column are samples with at least ten cells of the
state. Five patients carry enough cells to test — Ewing_2 (Cillo),
Ewing_Patient (Mao), ES-025 (Visser), SJEWS001321 (Waltner) and the single
urothelial carcinoma control patient, which contributes two of the eight
samples. Specificity is positive and above the sample's own reference cells in
6 of those 8 samples and 4 of the 5 patients; ES-025 fails on both of its
samples. Excluding the urothelial carcinoma control, **three** Ewing sarcoma
patients of 38 carry positive evidence. See also Figure 3C and Tables
f4_split_genome_specificity and f4_prevalence.

