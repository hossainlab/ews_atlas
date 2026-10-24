# 02 · Quality control

Implemented in `src/ewsatlas/qc.py`, driven by `config.qc`, executed by
`src/pipeline/01_ingest_qc.py`. Attrition is written to
`results/tables/qc_attrition.csv` and to the stage manifest.

## Policy

| filter | setting | applied |
|---|---|---|
| detected-gene floor | `qc.gene_floor_by_assay` — 200 (10x scRNA), 200 (CEL-seq2), 100 (10x snRNA) | per cell |
| mitochondrial cap | 30% | per cell |
| adaptive outliers | 5 MADs on log total counts, log detected genes, and % mito | **per sample** |
| doublets | scrublet, explicit threshold 0.25 | per sample, ≥100 cells |

## Why the gene floor is assay-aware — and the criterion is not

This is the one QC choice worth arguing in full, because it changes the atlas's
composition.

A detected-gene floor is a proxy for "is this barcode a real cell carrying
enough signal to assign a state". The map from that proxy to library quality is
assay-dependent. A nucleus holds roughly 10–20% of a cell's RNA and is depleted
of cytoplasmic and mitochondrial transcripts, so a nucleus with 150 detected
genes is a *normal nucleus*, not a degraded cell. One integer applied to both
assays is uniform only in arbitrary units; the effective stringency is far
harsher on the nuclear arm.

The failure mode is specific, and it is not noise. A fixed floor does not remove
nuclei at random — it removes the smallest and least transcriptionally active
ones. Because transcriptional output tracks cell type and cycling state, the
survivors skew toward large, active, proliferating malignant cells. A technical
threshold becomes a composition bias, which is precisely the artefact an atlas
exists to prevent.

**Measured cost — and a correction.** An earlier version of this page claimed
the Waltner nuclei have median 231 detected genes and 296 UMIs with 45.2% below
200, implying a uniform floor would cost ~30,000 nuclei. That figure was wrong.
Measured directly on the 13 ScPCA emptyDrops-filtered matrices (92,083 nuclei):

| | pooled | excluding the failed library |
|---|---|---|
| nuclei | 92,083 | 85,944 |
| median detected genes | 771 | 838 |
| median UMIs | 1,071 | 1,184 |
| median % mitochondrial | 2.14 | 1.99 |
| 95th percentile % mitochondrial | 15.8 | 10.2 |
| below 200 genes | 9,935 (10.8%) | 3,890 (4.53%) |
| below 100 genes | 6,123 (6.6%) | 476 (0.55%) |

So the assay-aware floor recovers **3,414 nuclei, about 4% of the arm**, not a
third of it. The directional argument above still holds — a whole-cell floor is
the wrong instrument for nuclei, and the median nuclear library here is indeed
shallower than a droplet whole-cell library — but the honest statement of its
magnitude is that it is a small correction, and this page previously oversold it.

The *criterion* stays uniform everywhere (per-sample 5-MAD outlier removal,
identical code path, identical parameters); the absolute floor is a per-assay
backstop against empty droplets. Guarded by `test_gene_floor_is_assay_aware`.

### Mitochondrial cap

The cap is assay-aware for the same reason as the floor, but in the opposite
direction. The four droplet whole-cell studies sit at median 2.8–3.9% MT with
p95 9–14%, so a 30% cap is effectively non-binding for them (a 20% cap would
cost 0–70 cells). Nuclei are *depleted* of mitochondrial transcripts — median
1.99%, p95 10.2% across the 85,944 nuclei that pass the library gate — so a
uniform 30% cap is lenient exactly where it should be strict: a nucleus above
10% MT carries cytoplasmic or ambient carryover, which inflates apparent
stromal signal in malignant nuclei, the same confounder as a doublet. The 10%
nuclear cap costs 4,441 of 85,944 nuclei (5.2%) against 366 (0.4%) at 30%.
CEL-seq2 plates run hotter once ERCC spike-ins are removed (per-plate medians
13–27%), so a 30% cap would cut around the median of several plates; 40% is the
backstop there.

(These figures were first computed against the full 92,083 nuclei, i.e. including
the failed library, which inflated them to 8,061 nuclei and 8.8% at a 10% cap.
The numbers above are measured on the gated population.)

## The library-level gate, which matters more than the floor

Profiling the ScPCA arm turned up something the per-cell policy cannot handle.
Library **SCPCL001111** has median 44 detected genes and 76 UMIs, with 98.5% of
its nuclei under 200 genes; ScPCA's own pipeline retained 89 of its 6,139 cells
(1.4%). Every other library in the deposit has a median of at least 237.

A per-sample MAD cannot catch this, because the MAD is computed *within* that
library — a uniformly failed library looks internally consistent and survives.
`qc.failed_libraries` therefore excludes any library whose **median** detected
genes falls below its assay floor, before per-cell filtering. One library fails
that test; 85,944 of 92,083 nuclei go forward to per-cell QC. The gate is
recorded per study in the attrition table (`drop_failed_library`,
`failed_libraries`) and guarded by `test_library_gate_catches_a_failed_library`,
which also asserts that the MAD alone does *not* catch it.

## Starting matrices for the ScPCA arm

ScPCA ships three matrices per library. The **processed** ones already enforce a
200-gene minimum (measured minimum per library 200–402, 0.0% below 200) — they
have already made the exact choice argued against above. Ingest therefore starts
from the **filtered** (emptyDrops) matrices, which retain nuclei down to 3
detected genes, and re-applies this project's policy on top. Starting from
`processed` would have made the assay-aware floor unobservable.

## Why outlier detection is per sample, not global

Library depth varies by an order of magnitude across these six deposits. A
global MAD cut would mostly separate deep libraries from shallow ones — that is,
it would detect the study, not the bad cells. Per-sample MADs make the criterion
depth-adaptive, which is the point of using MADs over fixed cuts at all.
`mad_outlier` returns all-False when the MAD is zero, so a constant metric
cannot wipe a sample.

## Doublets

**Why this is load-bearing here, not routine.** A malignant-plus-stromal
doublet reads as a malignant cell expressing a mesenchymal programme — which is
exactly the EWSR1-FLI1-low state this atlas sets out to study. Under-calling
doublets would manufacture the finding. That is why the doublet policy gets its
own section rather than a line in the table.

**Scores: scrublet, per sample, with an explicit cut.** scanpy's automatic
threshold requires `scikit-image`, whose install is refused by host policy here,
so the threshold is supplied rather than inferred. Samples under 100 cells are
`not_assessed` rather than forced through a simulation that cannot support them.

**Calls: the expected multiplet rate, not a score threshold.** The 10x multiplet
rate is ~0.8% per 1,000 cells recovered, so the expected doublet fraction scales
with library size — ~5.6% for a 7,000-nucleus library, ~11% for a 14,000-nucleus
one. A single score cut ignores that. `call_doublets(method="expected_rate")`
flags the top-scoring cells per sample at the rate the chemistry predicts; the
criterion (an instrument prior) is uniform, while the cut adapts to the quantity
that actually drives multiplet formation. `method="fixed_threshold"` is retained
for comparison and for an environment that can compute scrublet's own bimodal
threshold.

Measured on the Waltner arm (92,083 nuclei in, 12 libraries after the library
gate):

| doublet policy | nuclei removed | retention |
|---|---|---|
| fixed threshold 0.25 | 2,394 (2.9%) | 80,102 (87.0%) |
| expected rate, 0.8%/1k | 5,419 (6.6%) | 77,077 (83.7%) |

**Agreement with the depositor's independent caller is only partial, and that is
reported rather than tuned away.** ScPCA ran scDblFinder on the same nuclei. The
two scores correlate (Spearman rho 0.41) and ScPCA's doublet class has a higher
median scrublet score (0.142 vs 0.087), so they agree in direction — but after
the expected-rate cut, 6,753 of our 77,077 retained nuclei (8.8%) are still
scDblFinder doublets, ranging 5.8–10.7% per library. Raising our rate would not
fix this: the callers rank cells differently, and chasing agreement between two
heuristics is not evidence of correctness.

The honest mitigation is downstream and specific: the discovery step tests
directly whether fusion-low malignant cells are enriched for either caller's
doublet flag. If they are, the state is an artefact and the finding dies; if
they are not, the residual disagreement is not what explains it. That test is a
precondition of the discovery claim, not a robustness appendix.

### A failure is never reported as zero doublets

The first Waltner ingest marked **all 82,496 nuclei `failed`** — every
per-sample scrublet call raised, a bare `except` swallowed the exception, and
the attrition table dutifully reported `drop_doublet = 0`. A total failure read
as a clean result, and ScPCA's own calls were the only reason it surfaced.

The cause, once the exception was allowed to escape: joblib's threading backend
opens a named pipe for its worker pool, which this sandbox denies inside a
subprocess (`PermissionError: [WinError 5]`) while permitting it in a
long-lived kernel. The same code therefore passed interactively and failed under
`python src/pipeline/01_ingest_qc.py`. pynndescent passes an explicit `n_jobs`
that overrides `joblib.parallel_config`, so the backend itself is pinned to
`sequential` around the scrublet call. This is an environment workaround, not a
scientific choice — the result is identical, only slower.

Three changes followed:

- `call_doublets` returns `(table, errors)` and **raises** when any sample
  fails, unless the caller explicitly passes `fatal=False`.
- The attrition table carries `doublet_method`, `doublet_not_assessed` and
  `doublet_errors`, so assessed, not-assessed and failed cannot collapse into
  one number.
- The matrix is forced to canonical CSR (`sum_duplicates`, `sort_indices`)
  before scrublet, since `ad.concat(join="outer")` can leave non-canonical
  indices that h5ad serialization silently normalizes.

Guarded by `test_doublet_failure_raises_instead_of_reporting_zero`,
`test_small_samples_are_not_assessed_not_cleaned` and
`test_expected_rate_scales_with_library_size`.

### Cross-checking against the depositor's own calls

ScPCA's per-cell `scDblFinder` class and score, `miQC_pass` and
`prob_compromised` are retained as `src_*` provenance columns. They are never
used as this atlas's labels — two independent pipelines disagreeing is
information, and overwriting ours with theirs would discard it — but they are
the external reference that caught the failure above.

## Reporting flags, not filtered objects

`qc.flag_cells` returns a per-cell flag table and the caller decides what to
drop. This keeps attrition attributable per criterion, which is what the QC
table in a paper is actually built from, and makes "how many cells did the mito
cap cost?" answerable after the fact.

## What is not assessed

Ambient RNA is not corrected (no SoupX/CellBender step). The deposits do not
uniformly ship raw unfiltered droplet matrices, which ambient estimation needs,
so correcting some studies and not others would introduce exactly the kind of
study-dependent processing the atlas is meant to remove.
