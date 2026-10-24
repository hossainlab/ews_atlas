# 12. The atlas as a reference

Stages: `src/pipeline/08_reference_mapping.py`, `09_bulk_deconvolution.py`,
`10_surface_targets.py`
Outputs: `results/tables/f6_*.csv`, `results/pub_figs/Figure_6_atlas_utility.png`

A Resource paper is judged on whether anyone else can use the resource. This
stage answers that in three ways that do not depend on each other.

## 1. Surface targets, and why the comparison has to be paired

STEAP1, IL1RAP, CD276 (B7-H3), CD99 and NT5E (CD73) all carry active CAR-T,
CAR-NK or antibody-drug-conjugate programmes in Ewing sarcoma. What no single
study supplies is the denominator — how often the target is also carried by the
normal cells next to the tumour.

The pooled landscape is confounded and cannot be read as selectivity. 54% of
atlas cells are Waltner nuclei with a median of ~1,100 detected genes, and cell
types are not contributed evenly by the studies, so a cell type's detection rate
partly reports which assay produced it. The margin is therefore computed
**within each sample** and differenced against that sample's malignant cells,
holding assay, depth and patient fixed (47 samples qualify).

Median within-sample margin, in percentage points of cells positive. All ten
normal classes are shown with the number of contributing samples; the three
thinly sampled columns are not a basis for ranking targets. The plasma-cell
column present before the annotation revision has disappeared, because 307 of the
396 cells first called plasma cells are plasmacytoid dendritic cells and no sample
retains 20.

| target | T cell (45) | Mono/Mac (24) | Fibroblast (15) | Endothelial (14) | DC (12) | NK (10) | B cell (8) | Urothelial (4) | Mast (3) | pDC (2) |
|---|---|---|---|---|---|---|---|---|---|---|
| STEAP1 | +4.9 | +13.1 | +7.5 | +12.2 | +13.4 | +2.2 | +1.9 | −20.9 | +11.0 | +6.5 |
| IL1RAP | +7.7 | +4.7 | +13.1 | +6.3 | +16.7 | +6.6 | +2.9 | −21.7 | +10.0 | +10.3 |
| CD276 | +2.6 | +0.6 | −4.7 | +0.3 | +12.2 | +0.2 | +1.5 | +18.8 | +18.5 | +4.3 |
| CD99 | +3.5 | +19.8 | +2.0 | +2.3 | +24.2 | +0.5 | +6.4 | +7.9 | +6.8 | +8.7 |
| NT5E | +0.1 | +0.2 | 0.0 | −15.0 | +0.7 | 0.0 | −2.0 | −23.9 | +7.9 | −0.4 |

Reading, over the nine non-urothelial classes: STEAP1 (+1.9 to +13.4), IL1RAP
(+2.9 to +16.7) and CD99 (+0.5 to +24.2) are positive against every class. CD276
is positive against eight and negative against fibroblasts. NT5E is the weakest -
positive against four, flat against two, negative against three, worst against
endothelium. Both STEAP1 and IL1RAP are strongly negative against urothelium, a
normal-tissue warning that comes from Mao's bladder control rather than from any
Ewing sample.

**Margin and off-tumour burden are different axes, and CD99 is the case that
shows why.** Its within-sample margin is positive against all ten classes, but in
absolute terms it is detected in 74.2% of urothelial, 72.9% of fibroblast, 70.3%
of endothelial and 68.3% of platelet cells against 56.1% of Ewing sarcoma cells -
more normal cells carry it than tumour cells do. A positive margin says the
tumour carries the target more often than the normal cells *in that sample*; it
does not say the normal burden is low. Read on both axes the data favour IL1RAP
and STEAP1.

Waltner-exclusion (stage 12) changes one of these conclusions: the count of
classes with a positive margin is unchanged for STEAP1, IL1RAP, CD276 and NT5E,
but CD99 falls from 10 of 10 to 8 of 10, so its apparent universality is the one
selectivity result that depends on the single-nucleus arm.

Two caveats that belong next to any use of this table. Detection is RNA, so
per-cent-positive figures are lower bounds and a negative margin is stronger
evidence than a positive one. And the class labelled `DC` on the released object
is mast cells (docs/10, Figure 5), so that column must be read as mast.

## 2. Leave-one-study-out reference mapping

For each study in turn: train a scVI reference on the other five, map the
held-out study by scArches surgery, transfer labels by 15-nearest neighbours in
the shared latent space.

Three declared departures from the atlas's own integration configuration, all to
make six trainings affordable on CPU (no GPU is available in this sandbox):
2,000 highly variable genes rather than 5,000, 100 epochs rather than 500, and
the reference capped at ~50,000 cells per fold by trimming large
(study, cell type) groups while keeping small ones whole. Highly variable genes
are selected on the reference studies only — selecting them on the full atlas
leaks the query into the feature space and is the usual way this experiment is
accidentally made easy. `batch_key='sample'` with no extra categorical
covariates, since study and assay nest within sample and surgery onto unseen
covariate levels is a separate failure mode.

| held out | query cells | agreement | confident cells | agreement when confident |
|---|---|---|---|---|
| Cillo 2022 | 14,353 | 0.853 | 95.1% | 0.870 |
| Visser 2023 | 6,575 | 0.844 | 88.7% | 0.883 |
| Waltner 2025 | 74,888 | 0.822 | 81.2% | 0.876 |
| Goodspeed 2025 | 12,639 | 0.816 | 82.4% | 0.881 |
| Aynaud 2020 | 14,504 | 0.787 | 88.8% | 0.817 |
| Mao 2024 | 15,394 | 0.625 | 84.3% | 0.688 |

Mao's 0.625 is structural, not technical: its 2,433 urothelial cells have no
counterpart in the reference — 20 urothelial cells exist across the other five
studies combined — so they are necessarily called malignant. Excluding
urothelial and unassigned cells, Mao reaches 0.810 and the median across studies
is 0.822. The four classes added by the annotation revision transfer well - median
recall 0.97 (mast), 0.80 (NK), 0.75 (plasmacytoid dendritic), 0.71 (conventional
dendritic) - which is independent support for the revision. A cell type contributed by one study cannot be transferred when that
study is the query, and no amount of training fixes it.

Two other failures are worth reporting rather than smoothing over. Aynaud's 991
T cells are recovered at 1% and called Ewing sarcoma — the other five studies'
T cells do not support that annotation in an O-PDX sample. And
`Ewing_proliferating` is recovered at a median of 32%, which is what a
cycle-phase split of one population looks like when it is asked to behave like a
discrete type.

**What this measures.** The reference labels are the atlas's own annotation,
produced with the query study present. The experiment measures projection
stability under study removal, not generalization to an external cohort. No
external Ewing sarcoma single-cell cohort exists in any public archive searched
(docs/01), so the stronger experiment cannot be run here.

## 3. Bulk deconvolution and outcome

343 included primary records across five usable series supply expression; 129
unique patients supply survival after the accession-level duplications in
docs/07 are resolved.

Three corrections were needed before the numbers made sense, each recorded
because each produced a plausible-looking wrong answer first:

1. **Linear space.** NNLS on log2 microarray values returns ~0% malignant
   content. A mixture is additive on the linear scale only.
2. **Collinear malignant profiles.** With `Ewing_sarcoma` and
   `Ewing_proliferating` as separate reference columns, NNLS put 0.358 on the
   proliferating profile and 0.041 on the class eight times larger in the atlas.
   Bulk microarray cannot separate two profiles differing by cell-cycle genes;
   they enter as one.
3. **Whole-cell signature.** 68,698 of 97,586 frozen malignant cells are nuclei.
   A nucleus-dominated malignant profile does not resemble a bulk tissue extract:
   it returned mean malignant fractions of 0.33–0.53, against 0.73–0.79 once the
   signature was restricted to whole-cell libraries. Only the second is
   biologically plausible for this disease.
4. Per-platform missingness is by whole gene, so each cohort is deconvolved
   against markers it actually measures (896–900 of 1,000) and fractions are
   compared within cohort, never across.

Cox, stratified by cohort, hazard ratio per standard deviation, 129 patients:

| covariate | endpoint | HR (95% CI) | p |
|---|---|---|---|
| malignant fraction | OS | 0.77 (0.61–0.98) | 0.035 |
| malignant fraction | EFS | 0.82 (0.66–1.02) | 0.081 |
| stromal fraction | OS | 1.33 (1.04–1.71) | 0.025 |
| stromal fraction | EFS | 1.25 (0.98–1.58) | 0.068 |
| immune fraction | OS / EFS | 1.12 / 1.09 | 0.36 / 0.40 |
| fibroblast fraction | OS / EFS | 1.29 / 1.22 | 0.053 / 0.115 |

Nothing survives Bonferroni correction for the eight tests (threshold 0.0063),
and nothing here is new: this is the tumour–stroma association GSE63157 itself
published. It is reported as a positive control that the deconvolution carries
real signal, which is the only claim it can support.

Proportional-hazards assumptions were **not assessed**.

## Reproduce

```
python src/pipeline/10_surface_targets.py
python src/pipeline/09_bulk_deconvolution.py
python src/pipeline/08_reference_mapping.py   # ~1 h on 12 CPU cores
python src/pipeline/06_figures.py
```
