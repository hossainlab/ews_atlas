# EWSAtlas v4.0 — validation report

Object: `data/processed/EWSAtlas_v4.h5ad` — 64,702 cells x 17,017 genes, 5 studies, 56 samples, 33 Ewing patients.
Built from the original GEO deposits; v3 was deleted at the user's instruction and contributed nothing.

## 1. Integration

scVI (n_latent=30, 2 layers, NB likelihood), batch = `sample` (56 batches), categorical covariates `study` and
`assay`, continuous covariate `pct_counts_mt`, 2,000 seurat_v3 HVGs selected at study level, 120 epochs on CPU
(12.4 min). Per-study latent standard deviation is 0.69-1.00, i.e. no study occupies a degenerate region.

**Benchmark not run.** The plan called for scVI to be compared against the retained Harmony embedding; that
embedding lived only in v3, which was deleted, so no comparison was possible. The primary embedding
(`X_scvi`) is therefore declared by construction, not by benchmark.

**Residual study structure is substantial.** Measured on the stored latent: each cell's 15 nearest neighbours
are drawn from its own study **78.1%** of the time, against **21.1%** expected if studies mixed perfectly
(1.000 would be no mixing at all). Per study: Mao 0.902, Cillo 0.824, Aynaud 0.749, Goodspeed 0.711,
Visser 0.617. Malignant cells mix somewhat better than the rest (0.739 vs 0.816). At the cluster level, 16 of 36
communities contain all five studies and 25 contain at least four; only one cluster (19 cells) is single-study.

This metric conflates batch effect with genuine biology — the cohorts differ in patients, anatomical site and
material type (O-PDX vs patient tumour vs PBMC vs CTC), and Mao is a single bladder case with a
urothelial-carcinoma comparator, so some separation is expected rather than artefactual. Even so, 0.781 is high,
and no claim in this report should be read as "the studies are well mixed". Any cross-study comparison of
composition or state frequency should be stratified by study, and a formal benchmark (scIB metrics against
Harmony, scANVI and an unintegrated baseline) remains outstanding.

**Diagnostic that caught a real defect.** The first integration produced 304 leiden communities, 279 of them
~50-cell knots consisting almost entirely of Aynaud cells. Root cause was not the integration but the ingest:
four of five Aynaud samples ship `GRCh38_`-prefixed gene symbols (CellRanger v2 barnyard naming) while
PDX-1058 ships plain v3 features, so the gene space matched only PDX-1058 and 13,064 of 15,057 Aynaud cells
carried zero counts. After stripping the prefix the object has 0 zero-count cells, Aynaud's mitochondrial
fraction reads 3.7% (it had read 0.0%, because its MT genes were never being recognised), and clustering
resolves to 36 communities, most spanning all five studies.

## 2. Annotation — PROVISIONAL

36 leiden communities (resolution 1.0 on `X_scvi`) assigned by per-cluster standardized signature z-argmax over
14 marker panels, with two guards: urothelium is allowed only in clusters where the bladder study contributes
>=50% of cells, and human immune identities are refused for clusters >80% O-PDX (an orthotopic xenograft carries
no human immune compartment). Clusters whose best lineage z is <1.0 default to `Ewing_sarcoma`.

Three independent sanity checks pass:

| Check | Expectation | Result |
|---|---|---|
| Aynaud (O-PDX, human tumour cells only) | ~100% malignant | **98.8%** |
| Cillo (TIL + PBMC libraries) | immune-dominated | **69.8% immune** |
| Urothelium | confined to the bladder study | **2,537 of 2,543 cells** |

Two earlier assignment schemes were rejected after inspection: fraction-of-cells thresholds misfired across
studies whose median gene counts differ ~10-fold (B cells called Ewing; urothelium called inside a PDX-only
study), and rank-percentile argmax produced 4,403 spurious platelets.

**Known caveats.** CNV inference was not performed, so malignant identity rests on markers alone.
`MSC_Osteoblast` in the immune-rich studies is most likely SPP1+ macrophages rather than mesenchymal cells,
since SPP1 is shared between the panels. Both are properly resolved by scANVI label transfer against a
validated reference — recommended before any publication use.

## 3. Fusion-activity axis (externally calibrated)

Scored against the ESCLA knockdown signatures derived in Phase 1 (566 of 648 high genes and 389 of 431 low genes
map into the v4 gene space); `fusion_activity` = high score - low score, stored per cell.

| Stratum | n | median fusion activity |
|---|---|---|
| `Ewing_proliferating` | 8,885 | +0.058 |
| `Ewing_sarcoma` | 20,830 | +0.003 |
| Malignant, recurrence | 257 | +0.066 |
| Malignant, primary | 13,191 | -0.003 |
| Malignant, metastasis | 169 | -0.056 |

Proliferating malignant cells carry higher EWSR1-ETS activity than non-cycling malignant cells. The
disease-state contrast rests on very few metastasis and recurrence cells (169 and 257) and on the 48.7% of
malignant cells (14,468 of 29,715) whose `disease_state` is `not_reported`, so it is reported as a description of the atlas,
not as an inference about progression.

The axis inherits the Phase 1 limitation: ESCLA was normalized with RMA reimplemented in Python rather than the
publication's SST-RMA, which compresses fold changes, and NKX2-2, CCND1 and FLI1 move in the expected direction
but fall below the signature's fold-change cutoff.

## 4. Single-cell vs single-nucleus comparison — NOT ASSESSED

All five included studies are whole-cell (10x 3'/5' and CEL-seq2). The only single-nucleus arm is
Waltner/SCPCP000015, absent pending a user-held ScPCA download token, so there is no nucleus stratum to compare
against. The GSE200529 calibration control (7 dissociation-bias gene sets from paired warm/cold/nuclei contrasts
across 9 sarcoma PDX) is built and waiting. Until that arm is added, v4.0 carries **no** sc-vs-sn confound —
which is a property of its reduced scope, not a resolved problem.

## 5. Outcome linkage

Deconvolution of the public bulk cohorts using v4 cell-type profiles as the reference: signature matrix built
in linear CPM space from raw counts, ~100 specificity-selected marker genes per type, bulk de-logged from log2,
non-negative least squares per tumour, run per cohort so each platform contributes only its own measured genes.
408 tumours across 6 cohorts. GSE63155 was dropped from the results: it is a SubSeries of GSE63157 with identical
GSM accessions and would double-count 46 patients.

An initial attempt was discarded: NNLS applied to log-space data returned ~0% malignant for 369 Ewing tumours,
which is impossible. The corrected linear-space estimates give mean malignant fraction 0.64-0.78 by cohort,
consistent with the high tumour purity expected in this disease.

Cox proportional hazards per cohort, HR per standard deviation of the estimated fraction, cohort treated as a
stratum rather than pooled (129 patients have usable survival: GSE63157 n=85, GSE17674 n=44):

| Endpoint | Cohort | n | events | Feature | HR per SD (95% CI) | p |
|---|---|---|---|---|---|---|
| OS | GSE63157 | 85 | 25 | malignant | 0.599 (0.422-0.849) | 0.0040 |
| OS | GSE63157 | 85 | 25 | stromal | 1.647 (1.152-2.354) | 0.0062 |
| OS | GSE63157 | 85 | 25 | Ewing_proliferating | 1.488 (1.024-2.163) | 0.0372 |
| EFS | GSE63157 | 85 | 34 | malignant | 0.705 (0.517-0.962) | 0.0274 |
| OS | GSE17674 | 44 | 26 | malignant | 0.737 (0.507-1.070) | 0.1087 |
| OS | GSE17674 | 44 | 26 | stromal | 1.211 (0.818-1.791) | 0.3386 |

Higher estimated stromal content and lower malignant content associate with worse overall survival in the
85-patient COG cohort; the 44-patient Savola cohort shows the same directions without reaching significance.

**This recapitulates a published result rather than extending it.** GSE63157's own report is titled for the
prognostic importance of tumour-stromal interactions, so the finding is best read as evidence that the v4
reference deconvolves sensibly, not as a new discovery. Two further cautions: deconvolved "malignant fraction"
is largely tumour purity, which is a known confounder of array-cohort survival analyses; and Savola's EFS
*events* are derived from status fields rather than deposited, as recorded in Phase 1.

Proportional-hazards assumptions and residual diagnostics: not assessed.

## 6. Reference mapping — round-trip validated

`src/map_query_to_ewsatlas.py` performs scArches surgery on the stored scVI model (reference weights frozen,
query batch embeddings trained) and transfers labels by 15-nearest-reference-neighbour vote in the shared latent
space.

Round-trip check: 2,000 Visser cells were held out of the reference, given unseen batch levels and re-projected.
After 5 surgery epochs, **1,680 of 2,000 (84.0%)** recovered their atlas cell type. Residual disagreement is
almost entirely `Ewing_sarcoma` <-> `Ewing_proliferating` (192 + 25 cells), a continuum rather than a discrete
boundary; immune, stromal and endothelial identities transfer cleanly (T cells 378/381, monocytes 65/66,
B cells 21/21). Urothelial transfer was **not** measured: the held-out query was drawn from Visser, which
contains no urothelial cells, so that category is absent from the confusion matrix entirely.

Caveat: the held-out cells came from a study that is itself in the reference, so this measures projection
stability, not generalization to a genuinely new cohort.
