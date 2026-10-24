# 04 · Integration

**Status: executed. The primary embedding is scVI with `batch_key='sample'`,
declared by the benchmark below.**

## Model

scVI is the primary integration method. It models raw counts with a negative
binomial likelihood, conditions on batch and covariates in the encoder, and
yields a latent space that downstream tools (scANVI label transfer, scArches
reference mapping) extend without retraining from scratch. Harmony runs as the
baseline because it is fast, linear, and the standard comparison in every
benchmark; it corrects a PCA rather than the count model, so when the two agree
the structure is unlikely to be an artefact of either.

scGPT and other zero-shot foundation embeddings are deliberately excluded:
published benchmarks place them below scVI on batch integration, and they add a
large dependency surface for no measured gain here.

## Configuration, and where it comes from

| parameter | value | source |
|---|---|---|
| HVGs | 5,000, `seurat_v3`, batch_key `study` | NBAtlas used 5,000 with the same flavour |
| `n_latent` | 30 | |
| `n_layers` | 2 | NBAtlas |
| `gene_likelihood` | `nb` | |
| `encode_covariates` | `True` | NBAtlas / scArches-recommended |
| `use_layer_norm` | `both` | NBAtlas / scArches-recommended |
| `use_batch_norm` | `none` | NBAtlas / scArches-recommended |
| batch key | two arms: `sample`, `study` | see below |
| categorical covariates | `study`, `assay` | |
| continuous covariates | `pct_counts_mt` | |
| `max_epochs` | 500 + early stopping (patience 12) | NBAtlas used 500 + early stopping |

NBAtlas's own values were read from its STAR Methods (scvi-tools v0.16.4): 5,000
HVGs, `sample` covariate, `n_layers=2`, `encode_covariates=True`,
`use_layer_norm='both'`, `use_batch_norm='none'`, 500 epochs with early
stopping; UMAP at `n_pcs=15`, `n_neighbors=20`, `min_dist=0.6`; Harmony compared;
CopyKAT + SCEVAN for CNA; scANVI (`n_layers=4`) → scArches with an uncertainty
threshold of 0.5.

HVG selection is at **study** level rather than sample level because per-sample
`seurat_v3` fails with a loess singularity on the small CEL-seq2 plates.

## Why two batch-key arms

v4 trained with `batch_key='sample'` (56 levels) and measured poorly: each
cell's 15 nearest neighbours came from its own study **78.1%** of the time
against **21.1%** expected under perfect mixing (per study: Mao 0.902, Cillo
0.824, Aynaud 0.749, Goodspeed 0.711, Visser 0.617; malignant cells 0.739 vs
0.816 for the rest). Spending batch-specific parameters on 56 sample-level
batches across ~100k cells is a plausible contributor, so v5 trains both arms
and the choice is decided by the benchmark rather than assumed.

v4 also stopped at 120 epochs with early stopping never triggering, i.e. it was
under-trained — the second reason to rerun rather than reuse.

## What the v4 numbers do and do not license

v4's 36 clusters spanned all five studies in 16 cases and ≥4 studies in 25; one
cluster of 19 cells was single-study. That is reasonable co-clustering
alongside poor neighbour mixing — the two measure different things, and the
honest summary is that v4 integrated well enough to annotate and not well
enough to claim a shared latent space.

## Benchmark and the declared primary embedding

Four embeddings on one footing, all from the same 5,000-HVG count matrix:

| embedding | batch correction | bio conservation | total | same-study kNN |
|---|---|---|---|---|
| **scVI, `batch_key='sample'`** | **0.611** | **0.603** | **0.606** | **0.760** |
| scVI, `batch_key='study'` | 0.567 | 0.509 | 0.532 | 0.950 |
| Harmony on `sample` | 0.403 | 0.537 | 0.483 | 0.825 |
| unintegrated PCA | 0.317 | 0.516 | 0.437 | 0.984 |

`scib-metrics` 0.6.1 aggregates (Leiden NMI/ARI, isolated labels, silhouette
label, cLISI; BRAS, iLISI, kBET, graph connectivity, PCR comparison, sBEE),
computed on the 69,400 cells carrying depositor labels, with the unintegrated
PCA as the `pre_integrated_embedding_obsm_key` reference. Same-study kNN is
computed on all 138,353 cells. Full tables:
`results/tables/v5_scib_results.csv`, `v5_knn_mixing.json`,
`v5_embedding_benchmark.csv`.

**`sample` beats `study`, which refutes the hypothesis that motivated the two
arms.** v5 was run in two arms because v4's poor mixing suggested that spending
batch-specific parameters on 56 sample-level batches was the problem. It is
not: the sample arm wins on *both* axes (Leiden ARI 0.384 vs 0.053), and the
study arm is barely better than no integration at all on bio conservation. The
finer batch definition is the right one here.

**Training actually converged this time.** Both arms early-stopped — 340 epochs
for `sample` (validation ELBO 1,488.4) and 380 for `study` (1,496.7) — against
v4's 120 epochs with early stopping never triggering.

### What the same-study kNN figure does and does not mean

0.760 against 0.338 expected looks like weak mixing, and the ranking is
informative, but the absolute level is **not** interpretable as batch effect in
this atlas, because the studies do not sample the same cell types. Aynaud is
98.8% malignant orthotopic PDX; Cillo is sorted PBMC and TIL; Mao carries the
only urothelium. Two studies that share no cell type cannot have mixed
neighbourhoods, and no integration method should force them to.

The per-study pattern under the primary embedding is consistent with exactly
that: Visser 0.358 and Aynaud 0.485 — the studies whose compartments are shared
— against Waltner 0.875 and Mao 0.853, which carry 54% of all cells and the
only urothelium respectively.

A composition-adjusted expectation (study shares computed within each cell type
rather than overall) is the right statistic, and it requires the cell-type
labels that the annotation stage produces. It is reported there, not here.

## Reproducibility note

Training ran on a Modal A100 (`config.compute`, image
`im-DdwhygaSGL9V5gKurM7ksW`, spec_sha `83a26f236d5ca71c`), 4,565 s wall for both
arms plus the PCA. GPU training reproduces **to tolerance, not bit-exactly**:
cuDNN selects kernels at runtime and atomic reductions are order-dependent. The
manifest records the device so a divergence is attributable.

Two failures worth recording, both environment rather than science:

- scanpy's `sc.external.pp.harmony_integrate` is broken against harmonypy 0.2.0
  — Harmony itself converged in 5 iterations, then the wrapper assigned a
  transposed array and raised. Harmony is therefore called directly through
  `harmonypy.run_harmony` with an explicit shape assertion, and because it needs
  only the PCA, it was finished locally in 77 s rather than paying for a second
  GPU run.
- Any pynndescent-backed call (`sc.pp.neighbors`, `sc.pp.scrublet`,
  scib-metrics' kNN graphs) fails with `PermissionError: [WinError 5]` when run
  as a script in this sandbox, because joblib's threading backend opens a named
  pipe that is denied in a subprocess but allowed in a long-lived kernel. The
  workaround is isolated in `ewsatlas.compat.sequential_joblib`.
