# EWSAtlas figure captions

Figures carry no titles; each caption headline states the claim the figure is
meant to earn, in the style of Cell Press. Panel letters are bold uppercase.
All figures are 180 mm wide at 300 dpi with editable vector text
(`pdf.fonttype 42`). n values are cells unless stated otherwise.

---

## Figure 1. A harmonized single-cell reference for Ewing sarcoma spans six studies, 65 samples and 38 patients

**(A)** The included deposits. Per study: assay chemistry (droplet whole-cell,
droplet nucleus, or plate well), number of samples, number of patients, and
number of cells retained after quality control, with the integrated totals at
right. Sample count and cell count rank the studies in nearly opposite orders —
Visser 2023 contributes the most samples (24) and the fewest cells (6,575),
Waltner 2025 the reverse (12 samples, 74,888 cells) — which is why integration
is batched on sample rather than on study.
**(B)** The integrated embedding coloured by study, by assay, and by sample. For
sample, hues are spread around the colour wheel and shuffled; with 65 categories
a legend is unreadable, so that panel shows only whether any sample sits apart.
**(C)** Annotated embedding of all 138,353 cells, ringed by the composition of
each cell type. Every cell type gets an arc of **equal angular width** — the
arcs show composition, not abundance, since sizing them by cell count would give
the malignant compartment two thirds of the circle and compress every immune
type into a sliver. Each arc carries a solid outer band in the cell type's
colour and, inside it, one band per metadata variable (inner to outer: sample,
assay, study, Leiden cluster). A type drawn from a single sample or study reads
as a solid band; a pooled type reads as fine stripes. Inside the circle, density
contours and on-plot labels are drawn only for compact classes — a contour
around a dispersed class such as Unassigned or Plasma would enclose mostly empty
space and assert an island that is not there — so those types are identified by
their arc label alone. Cell types were assigned at the level of 34 Leiden communities
(resolution 1.0), not per cell, and each label additionally requires one of its
own obligatory markers above a floor; 2,809 cells (2.0%) with no positive marker
evidence are left Unassigned. Only compact classes are direct-labelled.
**(D)** Mean log-normalized expression and the fraction of expressing cells for
the 31 panel markers present in the 16,622-gene feature space (UPK1B is absent),
grouped by the compartment they mark.
See also Figures S1, S2, S4 and S6.

## Figure 2. Sample-level batch correction outperforms study-level correction on both axes

**(A)** scIB aggregate batch-correction against biological-conservation scores
for four candidate embeddings. scVI with `batch_key='sample'` leads on both axes
(total 0.606), which refuted the pre-registered expectation that study-level
batches would conserve more biology.
**(B)** Observed share of 15 nearest neighbours drawn from the same study,
against the share expected from cell-type composition alone. The residual is
real for T cells (0.635 observed against 0.225 expected) and monocytes (0.736
against 0.347), and fully explained by composition for urothelial cells,
dendritic cells and platelets.
**(C)** Same-study neighbour share per study for all four embeddings.
**(D)** Agreement with depositor annotations on the 69,400 cells that carry
them, as a row-normalized percentage. The 77%/18% split on the fibroblast row is
the unresolved boundary quantified.
See also Figure S3.

## Figure 3. A correlation statistic, not a parameter change, makes copy-number inference usable

**(A)** Malignant-vs-reference AUC per sample for three scoring statistics. The
v5 genome-wide mean of |CNV| has a median of 0.540; re-running with NBAtlas's
window of 250 genes and a wider reference raises it only to 0.597; correlating
each cell's profile with a per-sample consensus malignant profile (Tirosh et
al.) reaches 0.712, with 23 of 45 samples above 0.7; three samples remain
below 0.5 (PDX-184 0.435, SCPCS000490 0.452, CR-hPB120 0.459).
**(B)** Per-sample AUC for nine chromosome arms scored signed. Arm-level scoring
is **not** what rescued the caller: pooling two immune-rich studies gives
0.61-0.72 for chr8q and chr12q, but across all 45 scored samples chr8q is the
strongest arm at a median of 0.574, chr12q is 0.526, and chr16q at 0.492 shows no
loss signal.
**(C)** How each of the 97,586 frozen malignant cells earned its label. The
copy-number threshold is the 95th percentile of each sample's own reference
cells, which is deliberately insensitive: it is used only to rescue cells the
markers missed (3,781, all originally `Fibroblast_CAF` or `Unassigned`), never
to remove cells the markers found. 1,634 aneuploid calls on immune and
endothelial cells were discarded as false positives.
See also Figure S5.

## Figure 4. Most recurrent expression programs in a six-study atlas track capture, not biology

**(A)** The 28 metaprograms recovered by per-sample NMF (ranks 4-9, relative
expression, 45 samples with at least 100 frozen malignant cells), ordered by the
number of samples carrying them. 15 of 28 are confined to a single study — in
NBAtlas, none were.
**(B)** Spearman correlation of each candidate program with the technical
drivers that could explain it. The programs that read as hypoxia, mesenchymal
transition and cell cycle correlate with the single-nucleus-versus-single-cell
signature at 0.78, 0.69 and 0.64, and the hypoxia-like program with detected
genes at 0.52. Doublet score is not a driver for any candidate (|rho| <= 0.18).
**(C)** MP24, the candidate fusion-low program, scored by how each cell entered
the frozen malignant set: it is carried by the copy-number-rescued fibroblasts
(mean 0.503) and is near zero in every marker-defined class, so it is an
artefact of the rescue step rather than a tumour state.
See also Figure S5 and Table v6_metaprograms.

## Figure 5. The Ewing immune microenvironment resolves into thirteen subsets, correcting two atlas-level calls

**(A)** 29,114 immune cells re-clustered into 24 communities and annotated as
thirteen subsets plus a low-quality class, on a UMAP computed from the immune
subset's own neighbour graph; the inset locates these cells in the full atlas
(grey). The subset embedding is necessary rather than cosmetic: the atlas-wide
UMAP restricted to immune cells is a cut-out of Figure 1, in which the subsets
sit where they sat in a layout optimised for 138,353 cells. Clustering used the
existing scVI latent rather than a fresh integration of the subset (declared
deviation from NBAtlas, which re-integrated).
**(B)** Marker evidence for each subset. Two Figure 1 calls are corrected here:
all 715 cells labelled `DC` in the atlas are mast cells (TPSAB1, CPA3), and
1,267 NK cells were contained within `T_cell`. No subset proportion differs
between treatment arms after correction for multiple testing (smallest p =
0.0339 for NK across 13 subsets in 24 patients; the study-adjusted coefficient
for NK has the opposite sign to its raw group means, a study-level reversal).

---

## Supplementary figures

**Figure S1. Barcode fate and assay-aware thresholds.** (A) Per-study attrition
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
behind Figure 2A.

**Figure S4. Covariate completeness.** Percentage of samples per study reporting
each harmonized covariate. Fusion partner, sex and age are the gaps.

**Figure S5. The first copy-number attempt, retained as the record of what
failed.** (A) v5 genome-wide score by cell type, (B) per-sample AUC, (C) score
against sequencing depth, and (D) the label disagreements this left unresolved.
The score swings 0.011 across depth bins while malignant and reference cells
differ by 0.003. Superseded by Figure 3.

**Figure S6. Every sample is distributed across the shared structure.** The
integrated embedding of all 138,353 cells coloured by sample, with hues taken
evenly around the colour wheel and shuffled so neighbouring samples are unlike.
With 65 categories a legend is unreadable and no colour can be matched back to a
sample by eye, so the panel answers one question: whether any sample sits apart
from the rest. The fine speckle throughout the shared structure is the point —
where samples interleave, colour mixes at the pixel level. Per-sample
neighbour-mixing statistics are in Figure 2C and
`results/tables/v5_knn_mixing.json`.

**Figure S7. The six studies that make up the atlas.** Per study: assay
chemistry (droplet whole-cell, droplet nucleus, or plate well), number of
samples, number of patients, and number of cells retained after quality
control, with the integrated totals at right. Every value is computed from the
released object at render time, so the panel cannot drift from the atlas it
describes. Waltner 2025 contributes 54% of all cells from 12 samples, while
Visser 2023 contributes 24 samples totalling 6,575 cells — sample count and
cell count rank the studies in nearly opposite orders, which is why the
integration is batched on sample rather than on study.
