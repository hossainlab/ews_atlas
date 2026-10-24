# 06 · The discovery question

**Status: question fixed, power measured on v4, analysis pending v5.**

## Question

In cell lines, EWSR1::FLI1 activity fluctuates between a high state (proliferative,
fusion-target-high) and a low state (migratory, mesenchymal). Whether that axis
exists in patient tumours as a *patient-level* property — some tumours carrying
a larger fusion-low fraction than others — and whether that fraction shifts with
treatment, has not been tested across cohorts. A single study cannot test it:
the fraction has to be separable from study, assay and cell-cycle effects, which
needs several cohorts in one feature space.

Concretely:

1. Is per-patient fusion-low fraction reproducible across studies, or is it a
   study effect?
2. Does it separate from cell cycle, rather than relabelling it?
3. Does it differ between treatment-naive and post-treatment tumours?
4. Does a fusion-low signature carry prognostic information in bulk cohorts
   after adjusting for tumour purity?

## The axis, and its known confound

The fusion-activity signature was derived from ESCLA (GSE176339 + GSE212063,
107/108 Clariom D arrays usable — GSM5359664 is truncated in GEO): **648
EWSR1-ETS-high and 431 EWSR1-ETS-low genes** plus a graded axis
(`results/tables/fusion_activity_signatures.csv`).

Scored on v4, 566/648 high and 389/431 low genes mapped, and the axis behaved as
a **cell-cycle proxy**: `Ewing_proliferating` median +0.058 against
`Ewing_sarcoma` +0.003. That is the central methodological problem, not a side
note — the fusion-high programme and the proliferation programme overlap by
construction. Separating them (regressing cell-cycle scores, scoring within
cycle phase, and testing the axis in non-cycling malignant cells only) is a
prerequisite for any claim, which is why it is its own pipeline step.

A second caution: ESCLA's signature was re-derived here with RMA reimplemented
in Python rather than SST-RMA, because Bioconductor is unusable in this
environment. That compresses fold changes, and NKX2-2, CCND1 and FLI1 itself
fall below the cutoff as a result. The gene sets are usable as a direction, not
as a definitive list.

## What the data can and cannot support

Measured strata in v4 (malignant cells only):

| stratum | cells | patients | studies |
|---|---|---|---|
| patient tumour | 13,455 | 25 | 4 |
| orthotopic PDX | 14,468 | 5 | **1** |
| circulating tumour cells | 162 | 3 | 1 |
| normal tissue | 1,630 | — | 2 |
| primary | 13,191 | 23 | |
| metastasis | 169 | 3 | |
| recurrence | 257 | 2 | |
| treatment-naive | 9,196 | 17 | |
| post-treatment | 2,685 | 10 | |

Three consequences, stated before any analysis runs:

- **O-PDX fidelity is not answerable in v4.** The PDX stratum is one study
  (Aynaud), so "PDX versus patient" is perfectly confounded with "Aynaud versus
  everyone else". Waltner's matched material in v5 is what would make this
  testable; until then it is demoted to a secondary question.
- **CTC claims are unsupportable** at 162 cells from 3 patients.
- **The treatment contrast is the best-powered clinical axis** (9,196 vs 2,685
  cells; 17 vs 10 patients), which is why it is question 3 rather than an
  afterthought.

Also: 48.7% of v4 malignant cells (14,468/29,715) have `disease_state` =
`not_reported`. Any disease-state contrast is run on the annotated subset and
reports its own denominator.

## Falsification checks

A positive result is only reportable if it survives:

- **leave-one-study-out** — the patient-level effect must persist with each
  contributing study removed in turn;
- **random signatures** — size-matched random gene sets must not reproduce the
  effect;
- **variance partition** — patient must account for more latent variance than
  study for the cells in question.

## Novelty position

No integrated Ewing sarcoma single-cell atlas exists; EwS is the gap among
sarcomas that have one (osteosarcoma >110,000 cells / 17 samples,
10.1186/s12885-025-13714-y; chondrosarcoma; UPS; rhabdomyosarcoma; infantile
fibrosarcoma). Claims already occupied in the EwS literature, and therefore
*not* led with: fusion heterogeneity within a single PDX panel (Aynaud 2020);
CD73+ CAF-like populations; impaired antigen presentation (Visser); immune
signatures (Cillo); MIF–CD74 signalling and SLC40A1/C1QA macrophages; a T-cell
prognostic signature; neural-crest resemblance (Hong); and a pan-pediatric
immune meta-analysis spanning 26 studies and 272 samples.
