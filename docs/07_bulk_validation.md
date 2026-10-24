# 07 · Bulk cohorts as the outcome substrate

**Status: assembled and used for a v4-based deconvolution; to be re-run on v5.**

## Why bulk data is in this project at all

No source single-cell deposit carries survival. If the atlas is to say anything
about outcome, the association has to be made in cohorts that have it, by
projecting atlas-derived cell-state fractions into bulk expression. The bulk
data is **not part of the atlas object** and would be dropped if the project
stayed a pure resource.

## Cohorts

408 tumour samples × 21,360 genes, harmonized from 522 GEO records across seven
series (`results/tables/ews_bulk_expression.parquet`,
`ews_bulk_clinical.csv`).

Two accession-level duplications were found and are important for anyone reusing
these cohorts:

- **GSE63155 is a SubSeries of GSE63157** — identical GSMs. Treating them as two
  cohorts double-counts 46 patients.
- **GSE17618 and GSE17674 are two normalizations of the same 44 hybridizations**
  (median Pearson r = 0.994 between matched samples).

After de-duplication, survival is usable for **129 unique patients**: GSE63157
(n = 85, COG) and GSE17674 (n = 44, Savola). Savola's event-free survival events
are *derived*, not deposited.

## Deconvolution

NNLS against atlas-derived signatures, roughly 100 specificity-selected markers
per cell type, per cohort.

**A real error was caught here and is worth recording.** The first run was
performed on log-transformed expression and returned ~0% malignant content for
369 Ewing tumours — an impossible answer that should have been read as a bug
immediately. NNLS assumes the mixture is additive, which holds in linear space
and not after a log transform. Re-run in linear CPM space, mean malignant
fraction came to 0.64–0.78 per cohort, which is biologically sane for this
disease.

## Outcome association

Cox proportional hazards (`statsmodels.PHReg`), hazard ratio per standard
deviation, cohort as stratum, 129 patients:

| cohort | covariate | endpoint | HR (95% CI) | p |
|---|---|---|---|---|
| GSE63157 | malignant fraction | OS | 0.599 (0.422–0.849) | 0.0040 |
| GSE63157 | stromal fraction | OS | 1.647 (1.152–2.354) | 0.0062 |
| GSE63157 | Ewing_proliferating | OS | 1.488 (1.024–2.163) | 0.0372 |
| GSE63157 | malignant fraction | EFS | 0.705 (0.517–0.962) | 0.0274 |
| GSE17674 | (same directions) | OS/EFS | not significant | |

This **recapitulates GSE63157's own published tumour–stromal finding**, which is
the point: it is a positive control showing the deconvolution carries real
signal, not a new result. A fusion-low signature tested in these cohorts must be
adjusted for tumour purity, because purity is itself prognostic here.

Proportional-hazards assumptions were **not assessed**.
