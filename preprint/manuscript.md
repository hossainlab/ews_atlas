# EWSAtlas: a harmonized single-cell reference for Ewing sarcoma, and what it reveals about recurrent tumour programs

[[AUTHORS]]

## Summary

Ewing sarcoma is a fusion-driven paediatric bone and soft-tissue sarcoma, and the
conspicuous gap among sarcomas that already have an integrated single-cell
reference. We assembled EWSAtlas from all six public single-cell and
single-nucleus studies of Ewing sarcoma — 138,353 cells, 65 samples, 38 patients
— under one quality-control policy, one gene space, one controlled vocabulary and
one annotation pipeline, and defined the malignant compartment from copy number
rather than from the markers that name it. Two results follow that bear on how
such atlases are read. Batch correction at the level of samples outperforms
correction at the level of studies on both the batch and the biological axis,
refuting our pre-registered expectation. More consequentially, most expression
programs that recur across these six cohorts track how the cells were captured
rather than what they were doing: 15 of 28 recovered metaprograms are confined to
a single study, and the three most biologically appealing correlate with a
single-nucleus-versus-single-cell signature at 0.78, 0.70 and 0.64. A
fusion-activity axis built from an external knockdown experiment proves, once
cell-cycle genes are removed, to track proliferation as strongly as its own
canonical targets. We therefore present EWSAtlas as a reference rather than a
catalogue of new tumour states: it transfers cell-type labels to a study it never
saw at a median agreement of 0.82, recovers tumour content from 343 archived bulk
tumours well enough to reproduce a published tumour–stroma survival association,
and ranks five clinically pursued surface targets by a selectivity measure that
single cohorts cannot compute.

**Keywords:** Ewing sarcoma; single-cell RNA sequencing; single-nucleus RNA
sequencing; data integration; tumour microenvironment; EWSR1::FLI1; reference
mapping; CAR-T target selection

## Introduction

Ewing sarcoma is defined by a fusion between a FET-family gene and an ETS-family
transcription factor, most often EWSR1::FLI1, and arises in bone and soft tissue
predominantly in adolescents and young adults. The fusion is the disease's
founding lesion and its most studied feature, and a substantial literature
describes how its activity varies between cells, how that variation maps onto
migratory and proliferative phenotypes, and how the tumour's stroma and immune
infiltrate differ from those of other sarcomas. What the field has lacked is a
common substrate on which those claims can be compared.

Single-cell studies of Ewing sarcoma have been published by several groups, each
with its own dissociation protocol, chemistry, quality-control thresholds, gene
annotation and cell-type vocabulary. Aynaud et al. established fusion-activity
heterogeneity in patient-derived xenografts and cell lines [1]; Cillo et al.
contrasted the immune signatures of Ewing sarcoma and osteosarcoma [2]; Visser et
al. reported functionally impaired antigen-presenting cells [3]; Goodspeed et al.
profiled newly diagnosed tumours and circulating tumour cells [4]; Mao et al.
added a bladder primary with a urothelial control [5]; and Waltner et al.
contributed the largest single-nucleus cohort [6]. Individually each is a
five-to-fifteen-sample study. Collectively they are the raw material for a
reference, but only if the differences between them are removed deliberately
rather than averaged over.

Integrated disease references have changed how other tumours are studied.
NBAtlas, which harmonized seven neuroblastoma datasets from six publications,
provides the template this work follows: a resource paper that first establishes
that its dominant cluster is malignant, then interprets heterogeneity within it,
and finally demonstrates the atlas as a mapping surface for new data [7]. Among
sarcomas, osteosarcoma, chondrosarcoma, undifferentiated pleomorphic sarcoma and
rhabdomyosarcoma already have such references; Ewing sarcoma does not.

We set out to build one, and to test the kind of claim such an atlas is usually
made to support — that a recurrent tumour state exists across patients. The
resource behaved as intended. The claim did not survive the controls we built to
test it, and the reason it did not is itself the most transferable result here:
in a cohort that mixes whole-cell and nuclear capture, the programs that recur
most convincingly across studies are the ones that track the capture method.
Truong et al. showed this for dissociation protocols in sarcoma
patient-derived xenografts [8]; we find the same phenomenon in patient tumours
across six independent cohorts, where it is easily mistaken for biology because
it survives every within-study check. We report the atlas, the controls, and what
the atlas can and cannot be asked.

## Results

### A harmonized reference across six studies and 38 patients

Six public studies met inclusion — EWSR1::FLI1- or EWSR1::ERG-driven Ewing
sarcoma, with desmoplastic small round cell tumour, clear cell sarcoma,
FET–CREB-rearranged tumours and myxoid liposarcoma excluded, and cell lines
excluded from the atlas object. They contribute Aynaud 2020 (14,504 cells, 5
samples, 5 patients), Cillo 2022 (14,353 / 5 / 4), Goodspeed 2025 (12,639 / 15 /
11), Mao 2024 (15,394 / 4 / 2), Visser 2023 (6,575 / 24 / 11) and Waltner (74,888
/ 12 / 5) (Figure 1A). Quality control retained 138,353 of 163,578 barcodes
(84.6%) under assay-aware gene floors and mitochondrial caps, together with a
library-level gate that removes whole libraries whose median complexity falls
below their assay floor — a failure a per-cell rule cannot detect, because a
uniformly failed library looks internally consistent (Figure S1, Figure S2).

Sample count and cell count rank the studies in nearly opposite orders: Visser
contributes 24 samples and 6,575 cells, Waltner 12 samples and 74,888. Fifteen
cell types were assigned at the level of 34 Leiden communities rather than per
cell, each label additionally requiring one of its own obligatory markers above a
floor, leaving 2,809 cells (2.0%) `Unassigned` rather than forced to a nearest
label (Figure 1B, 1C). Coverage of the harmonized clinical covariates is given in
Figure S4, and the embedding coloured by individual sample in Figure S6.

Two cautions belong with every figure that follows. Waltner supplies 54% of all
cells and is the only single-nucleus study, so assay and cohort are perfectly
confounded for the nuclear arm. And the atlas is built entirely from public data,
which is the precedent NBAtlas set rather than a weakness, but it means no
external query cohort exists for reference-mapping validation.

[[FIGURE: Figure_1_atlas_overview.png]]

### Batch correction at the level of samples outperforms correction at the level of studies

Four embeddings were scored with scIB metrics against an unintegrated principal
component baseline [9]. scVI [10] with `batch_key='sample'` leads on both axes —
batch correction 0.611, biological conservation 0.603, total 0.606 — ahead of
scVI with `batch_key='study'` (0.532), Harmony [11] (0.483) and unintegrated PCA
(0.437) (Figure 2A). We had expected 65 sample-level batches to spend parameters
on too few cells each and over-correct; they do not.

Apparent batch structure must be read against composition. The observed
same-study share of 15 nearest neighbours is 0.760 against 0.338 expected under
random mixing, but the studies sample different compartments, so roughly half of
that residual is composition rather than batch (Figure 2B, 2C). Training curves for
both batch definitions, the individual scIB metrics and the embedding coloured by
study, assay and sample are in Figure S3.

Agreement with the labels the depositors themselves released is 93.0% over the
69,400 cells that carry them, and the disagreement is concentrated in one row:
cells deposited as `Stromal` are 76.1% malignant and 23.0% stromal under our
labelling, and cells deposited as `Epithelial` are 98.2% malignant (Figure 2D).
That boundary — malignant against fibroblast — recurs throughout this work.

[[FIGURE: Figure_2_integration_quality.png]]

### A correlation statistic, not a parameter change, makes copy-number inference usable

Defining the malignant compartment from markers alone is circular when the
question is which cells are malignant. We scored copy number per study and
classified cells per sample. Our first implementation, the genome-wide mean of
|CNV|, separated malignant from reference cells at a median per-sample AUC of
0.540, which is not usable. Neither a wider window nor a wider reference set
repaired it (0.597; Figure S5). What repaired it was the statistic: correlating each cell's
profile with a per-sample consensus malignant profile, seeded from high-signal
non-reference cells after Tirosh et al. [12], reaches a median AUC of 0.712, with
23 of 45 scored samples above 0.7 (Figure 3A). Three samples remain below 0.5
(PDX-184 0.435, SCPCS000490 0.452, CR-hPB120 0.459) and are reported as such.

Arm-level scoring, which we expected to be the repair for a disease with few
focal changes, is not. Across all 45 samples chr8q is the strongest arm at a
median AUC of 0.574, chr12q reaches 0.526 and chr16q 0.492 (Figure 3B). The
favourable figures we first obtained, 0.61–0.72, came from pooling two
immune-rich studies at study level and do not survive per-sample scoring.

The frozen malignant set is 97,586 cells, each carrying the evidence that placed
it there: 12,114 with both markers and copy number, 67,587 with markers where
copy number was negative, 14,104 where copy number could not be assessed, and
3,781 admitted on copy number alone (Figure 3C). The copy-number threshold is the
95th percentile of each sample's own reference cells and is used only to rescue
cells the markers missed, never to remove cells the markers found; 1,634
aneuploid calls on immune and endothelial cells were discarded as false
positives.

[[FIGURE: Figure_3_malignant_compartment.png]]

### The immune compartment resolves into thirteen subsets and corrects four atlas-level calls

Re-clustering 29,114 immune cells on their own neighbour graph resolves thirteen
subsets — T_CD4 7,532; T_cytotoxic 5,860; Monocyte 4,974; Macrophage 3,351;
B_cell 1,681; NK 1,267; cDC2 1,137; Mast 715; Treg 670; Macrophage_LAM 532; pDC
307; Platelet 218; Plasma 88 — plus 782 low-quality cells (Figure 4A, 4B).

The subset view corrected four atlas-level calls that a coarse pass cannot avoid,
and the released object has been revised accordingly: all 715 cells first labelled
`DC` are mast cells (mean log-normalized TPSAB1 3.80, TPSB2 4.46, KIT 2.62,
against CD1C 0.00); 1,267 NK cells sat inside `T_cell`; the 1,137 conventional
dendritic cells sat inside `Monocyte_Macrophage`, so the atlas had no
dendritic-cell class at all; and 307 plasmacytoid dendritic cells sat inside
`Plasma`. Twelve classes became fifteen. The copy-number reference set is
unchanged by the revision, so the malignant calls above still stand, and
`cell_type_v5_original` preserves the pre-revision labels.

No subset proportion differs by treatment arm after correction for multiple
testing: the smallest p across 13 subsets in 24 patients is 0.0339 (NK), and its
study-adjusted coefficient has the opposite sign to the raw group means — 0.053
treatment-naive against 0.013 post-treatment — a study-level reversal rather than
a biological shift.

[[FIGURE: Figure_4_immune_microenvironment.png]]

### Most recurrent expression programs track capture, not biology

Per-sample non-negative matrix factorisation across the 45 samples with at least
100 malignant cells, on relative expression as Gavish et al. specify [13],
yielded 1,755 programs, of which 864 were robust, clustering into 28
metaprograms. Fifteen of the 28 are confined to a single study (Figure 5A); in
NBAtlas, none were.

Testing the candidates against the technical variables that could explain them
shows the problem is not marginal (Figure 5B). The programs that read as hypoxia,
mesenchymal transition and cell cycle correlate with a
single-nucleus-versus-single-cell signature at Spearman 0.781, 0.695 and 0.642,
and the hypoxia-like program correlates with detected genes at 0.524. Doublet
score explains none of them (|rho| <= 0.183).

The same test dissolved our leading biological hypothesis. An EWSR1::FLI1
activity axis derived from the Ewing Sarcoma Cell Line Atlas [14] — an external
knockdown experiment chosen precisely so that fusion activity would not be
defined from atlas markers — correlates with S phase at 0.430 and G2/M at 0.468,
against 0.200 for NR0B1 and 0.226 for NKX2-2. After removing 133
cycle-associated genes it still tracks cycle (0.246, 0.256) as strongly as its
own canonical targets (0.229, 0.255). The mechanism is unremarkable in hindsight:
knocking down the fusion arrests growth in the source experiment, so the
signature carries proliferation by construction.

Truong et al. reported warm-dissociation and gene-length artifacts in 125,831
cells and nuclei from Ewing sarcoma, desmoplastic small round cell tumour and
osteosarcoma patient-derived xenografts, and showed they are computationally
correctable [8]. Our contribution is the same phenomenon in patient tumours across
six independent cohorts, with the metaprogram decomposition that shows how it is
mistaken for biology. Because Waltner is our only single-nucleus study, assay and
cohort cannot be separated within this atlas; what supports assay as the cause is
mechanistic — the single-study programs are built from long, intron-rich genes —
together with Truong's matched cell and nucleus pairs, in which study is held
constant by design.

[[FIGURE: Figure_5_recurrence_and_confounding.png]]

### The atlas as a reference: target selectivity, label transfer and bulk deconvolution

**Surface targets.** STEAP1, IL1RAP, CD276 (B7-H3), CD99 and NT5E (CD73) all
carry active chimeric antigen receptor or antibody–drug-conjugate programmes in
this disease. What no single study can supply is the denominator. Pooled across
the atlas the targets are detected in 56.6% (CD99), 28.1% (IL1RAP), 22.0%
(STEAP1), 9.2% (CD276) and 1.4% (NT5E) of malignant cells (Figure 6A). Because
detection depends on assay and the studies do not contribute cell types evenly,
selectivity is computed as a paired within-sample contrast against that sample's
own normal cells (Figure 6B). Across the nine non-urothelial microenvironment
classes, STEAP1 (+1.9 to +13.4 percentage points), IL1RAP (+2.9 to +16.7) and
CD99 (+0.5 to +24.2) are positive against every class; CD276 is positive against
eight but negative against fibroblasts (−4.7); NT5E is positive against four,
flat against two and negative against three, including −15.0 against endothelium.

Margin and off-tumour burden are separate axes, and CD99 is the case that shows
why: its within-sample margin is positive everywhere, yet in absolute terms it is
carried by 74.2% of urothelial, 72.9% of fibroblast, 70.3% of endothelial and
68.3% of platelet cells against 56.1% of Ewing sarcoma cells — more normal cells
carry it than tumour cells do. Read on both axes the data favour IL1RAP and
STEAP1; both, however, are strongly negative against urothelium (−21.7, −20.9),
consistent with their documented normal expression and a reminder that this atlas
samples the tumour's neighbours, not every tissue a therapy will reach.

**Label transfer.** Training a scVI reference on five studies and mapping the
sixth by scArches surgery [15] recovers the full-atlas annotation at a median
agreement of 0.819 — 0.853 (Cillo), 0.844 (Visser), 0.822 (Waltner), 0.816
(Goodspeed), 0.787 (Aynaud) and 0.625 (Mao) — rising to 0.688–0.883 among the
81–95% of query cells where at least 60% of neighbours agree (Figure 6C). Mao's
value is structural rather than technical: its 2,433 urothelial cells have no
counterpart in a reference built from the other five studies, and excluding
urothelial and unassigned cells it reaches 0.810, with a median of 0.822 across
folds. The four classes added by the annotation revision transfer well — median
recall 0.97 (mast), 0.80 (NK), 0.75 (plasmacytoid dendritic), 0.71 (conventional
dendritic) — which is independent support for the revision. Because the reference
labels were produced with the query study present, this measures projection
stability under study removal and not generalization to an external cohort.

**Bulk deconvolution.** Using the atlas as a cell-type reference for 343 archived
bulk tumours, de-duplicated from 408 expression records across six series,
recovers mean malignant fractions of 0.73–0.83 across five cohorts. In the 129
patients with survival, higher tumour content tracks better overall survival
(hazard ratio 0.76 per standard deviation, 95% CI 0.59–0.97) and higher stromal
content worse (1.32, 1.03–1.70) (Figure 6D). Neither survives Bonferroni
correction for the eight tests shown, and neither is new: this is the
tumour–stroma association one of these cohorts itself published, recovered here
as evidence that the deconvolution carries real signal rather than as a finding.

Three corrections were required before these fractions were credible, and each
produced a plausible wrong answer first. Non-negative least squares on
log-transformed expression returns approximately 0% malignant content, because a
mixture is additive only on the linear scale. With `Ewing_sarcoma` and
`Ewing_proliferating` as separate reference profiles — collinear, differing
mainly by cell-cycle genes — the fit placed 0.358 on the proliferating profile and
0.041 on the class 8.4-fold larger in the atlas. And a malignant signature built
from the atlas as a whole is nucleus-dominated, since 68,698 of 97,586 malignant
cells are nuclei, and does not resemble a bulk tissue extract: it returned mean
malignant fractions of 0.33–0.53, against 0.73–0.83 once the signature was
restricted to whole-cell libraries.

[[FIGURE: Figure_6_atlas_utility.png]]

### What survives removing the largest study

Waltner contributes 54% of cells, so every per-sample headline was recomputed over
the 53 non-Waltner samples (63,465 cells, 33 patients, 28,888 malignant cells;
Table S1). The copy-number statistic holds (median AUC 0.724 against 0.712); 20 of
28 metaprograms survive, 8 of them single-study; and the number of normal classes
against which each target carries a positive margin is unchanged for STEAP1 (9 of
10), IL1RAP (9 of 10), CD276 (9 of 10) and NT5E (4 of 10). CD99 is the exception,
falling from 10 of 10 to 8 of 10, so its apparent universality is the one
selectivity result that depends on the nuclear arm. The headline that cannot be
re-tested this way is the capture result itself, since dropping Waltner removes
that arm entirely.

## Discussion

EWSAtlas does for Ewing sarcoma what integrated references have done for other
tumours: it provides one latent space, one annotation pipeline, one controlled
vocabulary and one mapping surface for data that were previously comparable only
in prose. Its most useful outputs are denominators — quantities that require many
patients and cannot be computed from any single cohort. The surface-target
analysis is the clearest instance. Deciding whether STEAP1, IL1RAP, CD276, CD99
or NT5E is a defensible target requires knowing how often the normal cells beside
the tumour carry it, in enough patients that one tumour's stroma does not settle
the question. That table exists now, and it does not say what the single-cohort
literature would predict: CD99, the diagnostic marker of the disease, is carried
by more normal cells than tumour cells in absolute terms, while IL1RAP and STEAP1
combine a positive within-sample margin with a low normal burden — and both carry
a urothelial warning that only a cohort containing a bladder control could have
produced.

The result with the widest reach, however, is cautionary. Most programs that
recur across these six cohorts track capture rather than biology, and they do so
in a way that is almost designed to deceive: they are reproducible, they
correlate with plausible hallmarks, and they survive within-study checks. Only an
explicit confrontation with technical signatures — nuclear versus whole-cell,
dissociation, detection depth — separates them from biology, and in our hands
that confrontation also destroyed the fusion-activity axis we had built the
project around. That axis had been constructed from an external perturbation
experiment precisely to avoid circularity, and it still failed, because the
perturbation arrests growth and therefore encodes proliferation. We take two
lessons from this. The first is that a disease atlas assembled from
heterogeneously captured public data should treat capture as a primary
explanatory variable rather than a nuisance to be regressed out and forgotten.
The second is that signature provenance does not guarantee independence: a
signature derived from a perturbation inherits every downstream consequence of
that perturbation.

Our results also bear on two specific claims in the Ewing sarcoma literature. The
CD73-positive, matrisome-high malignant population described by Wrenn et al. [16]
is genetically supported here, but in three patients of 38 (Figure S7) — its cells
carry their own patient's copy-number profile on a held-out half of the genome,
yet the state is not recurrent at cohort scale. And the malignant-versus-fibroblast boundary
that this question turns on is the same boundary on which our labels most
disagree with the depositors' own. Any atlas-scale statement about
mesenchymal tumour states in this disease will be decided by how that boundary is
drawn, and we therefore ship the evidence class for every malignant cell rather
than a single binary call.

What EWSAtlas does not provide is a new mechanism. We looked for one, with
pre-specified tests, and the honest outcome is that the recurrent states we found
are explained by how the data were made. For the field, the value of this
resource lies in the controls, the denominators and the mapping surface it
supplies to the next study — and in the negative results that tell others which
questions these data cannot answer.

## Limitations of the study

One study supplies 54% of cells and is the only single-nucleus dataset, so assay
and cohort cannot be separated within this atlas. Reference mapping is validated
against held-out studies rather than an external cohort, because no external
Ewing sarcoma single-cell cohort exists in any public archive we searched.
Treatment contrasts rest on 28 patients from the two studies that carry both
arms; the xenograft contrast on four matched patients from one study, where 0.125
is the smallest attainable two-sided p; and circulating-tumour-cell claims are not
supportable at all, since 405 of 6,941 cells from such libraries are malignant.
The CD73-high fibroblast-like state is genetically supported in three Ewing
patients of 38, too few for a prevalence claim. Proportional-hazards assumptions
for the survival models were not assessed. Target detection is RNA rather than
protein, so positivity figures are lower bounds on surface presentation and a
negative selectivity margin is stronger evidence than a positive one. The saved
scVI model reproduces the released embedding globally but not locally, and the
two latent spaces are shipped side by side for that reason.

## STAR Methods

### Resource availability

**Lead contact.** [[LEAD_CONTACT]]

**Materials availability.** This study generated no new reagents.

**Data and code availability.** The analysis object, a CELLxGENE-schema deposit
object (138,353 cells x 16,286 Ensembl-identified genes) and the trained scVI
model for reference mapping are deposited at [[ZENODO_DOI]]. All analysis code is
at [[REPOSITORY_URL]]. Source data are the six public studies cited, under their
original accessions: GSE130024, GSE198896, GSE277083, GSE250523, GSE243347 and
ScPCA SCPCP000015. Any additional information is available from the lead contact
on request.

### Method details

**Study inclusion.** Candidate studies were required to profile EWSR1::FLI1- or
EWSR1::ERG-driven Ewing sarcoma at single-cell or single-nucleus resolution in
patient-derived material. Tumours driven by other FET fusions were excluded, as
were cell-line-only deposits. Eighteen candidates were triaged; six qualified.

**Quality control.** Thresholds are assay-aware in both directions: detected-gene
floors of 200 (10x whole-cell, CEL-seq2) and 100 (10x nuclei), and mitochondrial
caps of 30%, 40% and 10% respectively, the nuclear cap being strict because
nuclei are depleted of mitochondrial message. Within those bounds the criterion is
uniform, a per-sample 5-median-absolute-deviation rule. A library-level gate
excludes libraries whose median detected genes fall below their assay floor.
Doublet calls use the chemistry's expected multiplet rate rather than a fixed
score threshold.

**Gene space.** The feature space is the intersection of per-study annotations
restricted to genes detected in at least 3 cells in at least one study, giving
16,622 genes.

**Integration.** 5,000 highly variable genes were selected at study level
(`seurat_v3` on counts). scVI was trained with 30 latent dimensions, 2 layers, a
negative-binomial likelihood, covariate encoding, `batch_key='sample'` and study,
assay and mitochondrial fraction as covariates, with early stopping. Candidate
embeddings were scored with scIB.

**Annotation.** Labels were assigned per Leiden community with obligatory-marker
guards: each label requires one of its own markers above a floor, nucleated immune
labels require PTPRC, NK is refused when CD3 is comparable, and a community below
a z-score floor is left unassigned.

**Copy number.** inferCNV was run per study with a 250-gene window and step 10,
using immune, endothelial and urothelial cells as reference. Each cell was
correlated with a per-sample consensus malignant profile seeded from high-signal
non-reference cells; the aneuploidy threshold is the 95th percentile of each
sample's reference cells.

**Metaprograms.** Per-sample NMF at ranks 4–9 on relative expression (each gene
centred within sample, negatives clipped), with robust programs clustered into
metaprograms and each candidate correlated against detection depth, mitochondrial
fraction, doublet score, cell-cycle scores and dissociation and nuclear
signatures.

**Reference mapping.** For each held-out study a scVI reference was trained on the
other five with highly variable genes selected on the reference only; the query
was mapped by scArches surgery and labels transferred by 15 nearest neighbours.

**Bulk deconvolution.** Non-negative least squares in linear space against a
signature built from whole-cell libraries only, with the two malignant classes
entered as one profile and markers selected per cohort from the genes that cohort
measures.

### Quantification and statistical analysis

Patient-level analyses treat the patient, not the cell, as the unit and include
study as a covariate. Survival models are Cox proportional hazards stratified by
cohort, reporting hazard ratios per standard deviation. Multiple testing is
addressed by Bonferroni correction and stated per figure. Exact n values appear
in every figure legend.

[[CONTRIBUTIONS]]

## References

1. Aynaud M-M, Mirabeau O, Gruel N, et al. Transcriptional programs define
   intratumoral heterogeneity of Ewing sarcoma at single-cell resolution. Cell
   Reports. 2020;30(6):1767–1779.e6. doi:10.1016/j.celrep.2020.01.049
2. Cillo AR, Mukherjee E, Bailey NG, et al. Ewing sarcoma and osteosarcoma have
   distinct immune signatures and intercellular communication networks. Clinical
   Cancer Research. 2022;28(22):4968–4982. doi:10.1158/1078-0432.CCR-22-1471
3. Visser LL, Bleijs M, Margaritis T, et al. Ewing sarcoma single-cell
   transcriptome analysis reveals functionally impaired antigen-presenting cells.
   Cancer Research Communications. 2023;3(10):2158–2169.
   doi:10.1158/2767-9764.CRC-23-0027
4. Goodspeed A, Bodlak A, Duffy AB, et al. Single-cell RNA sequencing of Ewing
   sarcoma tumors demonstrates transcriptional heterogeneity and clonal
   evolution. Clinical Cancer Research. 2025;31(10):2010–2023.
   doi:10.1158/1078-0432.CCR-24-2040
5. Mao W, Ding J, Li Y, et al. Single-cell RNA sequencing and spatial
   transcriptomics of bladder Ewing sarcoma. iScience. 2024;27:110921.
   doi:10.1016/j.isci.2024.110921
6. Waltner OG, et al. Multimodal single-cell analyses reveal distinct
   fusion-regulated transcriptional programs in Ewing sarcoma. eLife.
   doi:10.7554/eLife.109124
7. Bonine N, Zanzani V, Van Haver A, et al. NBAtlas: a harmonized single-cell
   transcriptomic reference atlas of human neuroblastoma tumors. Cell Reports.
   2024;43:114804. doi:10.1016/j.celrep.2024.114804
8. Truong DD, Lamhamedi-Cherradi S-E, Ludwig JA, et al. Dissociation protocols
   used for sarcoma tissues bias the transcriptome observed in single-cell and
   single-nucleus RNA sequencing. BMC Cancer. 2023;23:1105.
   doi:10.1186/s12885-023-10977-1
9. Luecken MD, Büttner M, Chaichoompu K, et al. Benchmarking atlas-level data
   integration in single-cell genomics. Nature Methods. 2022;19(1):41–50.
   doi:10.1038/s41592-021-01336-8
10. Lopez R, Regier J, Cole MB, Jordan MI, Yosef N. Deep generative modeling for
    single-cell transcriptomics. Nature Methods. 2018;15(12):1053–1058.
    doi:10.1038/s41592-018-0229-2
11. Korsunsky I, Millard N, Fan J, et al. Fast, sensitive and accurate
    integration of single-cell data with Harmony. Nature Methods.
    2019;16(12):1289–1296. doi:10.1038/s41592-019-0619-0
12. Tirosh I, Izar B, Prakadan SM, et al. Dissecting the multicellular ecosystem
    of metastatic melanoma by single-cell RNA-seq. Science.
    2016;352(6282):189–196. doi:10.1126/science.aad0501
13. Gavish A, Tyler M, Greenwald AC, et al. Hallmarks of transcriptional
    intratumour heterogeneity across a thousand tumours. Nature.
    2023;618(7965):598–606. doi:10.1038/s41586-023-06130-4
14. Orth MF, Surdez D, Faehling T, et al. Systematic multi-omics cell line
    profiling uncovers principles of Ewing sarcoma fusion oncogene-mediated gene
    regulation. Cell Reports. 2022;41(10):111761.
    doi:10.1016/j.celrep.2022.111761
15. Lotfollahi M, Naghipourfar M, Luecken MD, et al. Mapping single-cell data to
    reference atlases by transfer learning. Nature Biotechnology.
    2022;40(1):121–130. doi:10.1038/s41587-021-01001-7
16. Wrenn ED, Apfelbaum AA, Rudzinski ER, et al. Cancer-associated
    fibroblast-like tumor cells remodel the Ewing sarcoma tumor microenvironment.
    Clinical Cancer Research. 2023;29(24):5140–5154.
    doi:10.1158/1078-0432.CCR-23-1111
17. Wolf FA, Angerer P, Theis FJ. SCANPY: large-scale single-cell gene expression
    data analysis. Genome Biology. 2018;19:15. doi:10.1186/s13059-017-1382-0
18. Bergsma FJ, et al. The immune landscape of pediatric (extra)cranial solid
    tumors: a systematic review and integrative analysis. European Journal of
    Cancer. 2025;228:115708. doi:10.1016/j.ejca.2025.115708
