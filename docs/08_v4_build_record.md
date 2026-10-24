# 08 · Build record for EWSAtlas v4.0

v4 was built interactively, before the pipeline package existed. This page is
the retrospective record of what was actually run, so that v4's numbers remain
interpretable and the differences from the shipped scripts are explicit. It is a
record, not a reproduction recipe — **v5 is the reproducible object**.

## Result

| | |
|---|---|
| cells × genes | 64,702 × 17,017 |
| studies | 5 (no Waltner snRNA arm) |
| samples | 56 |
| EwS patients | 33 |
| cell types | 15 |
| malignant fraction | 45.9% |
| QC retention | 89.8% (64,702 / 72,023) |

Composition sanity checks that passed: Aynaud orthotopic PDX 98.8% malignant,
Cillo 69.8% immune, 2,537 of 2,543 urothelial cells inside the bladder study.

## Steps as executed

1. **Per-study ingest** — 10x MTX for Aynaud and Mao; nested tar.gz MTX for
   Cillo; gene × well plate tables for Visser; chunked dense TSV for Goodspeed.
   This is the logic now in `src/ewsatlas/io.py`.
2. **Uniform QC** — 200-gene floor, 30% mito cap, per-sample 5-MAD on log
   counts / log genes / mito, scrublet at an explicit 0.25.
3. **Feature space** — 17,017 genes: annotation intersection ∩ detected in ≥ 3
   cells in ≥ 1 study, barnyard prefixes stripped.
4. **scVI** → neighbours → UMAP → leiden (`flavor="leidenalg"`).
5. **Marker annotation** — standardized z-argmax with the urothelial and PDX
   guards (see [05](05_annotation.md)).
6. **Write h5ad.**

## scVI as actually trained in v4

Read back from `data/interim/scvi_model_v4/model.pt`:

| parameter | v4 value | v5 config |
|---|---|---|
| HVGs | 2,000 | 5,000 |
| `n_latent` | 30 | 30 |
| `n_layers` | 2 | 2 |
| `n_hidden` | 128 | 128 |
| dropout | 0.1 | 0.1 |
| likelihood | nb, dispersion per gene | nb |
| `use_observed_lib_size` | True | True |
| batch key | `sample` (56 levels) | two arms: `sample`, `study` |
| categorical covariates | `study`, `assay` | same |
| continuous covariates | `pct_counts_mt` | same |
| epochs | **120**, 12.4 min, final ELBO 696.8 | 500 + early stopping |

1,352,067 parameters. **Early stopping never triggered at 120 epochs** — v4 was
under-trained, which with the 2,000-HVG space is the main reason v5 retrains
rather than reuses.

## Known defects in v4, carried forward as v5 requirements

- **Integration mixing is poor.** Each cell's 15 nearest neighbours come from its
  own study 78.1% of the time against 21.1% expected under perfect mixing (Mao
  0.902, Cillo 0.824, Aynaud 0.749, Goodspeed 0.711, Visser 0.617; malignant
  0.739 vs 0.816 for the rest). Co-clustering is better than that figure
  suggests — 16 of 36 clusters span all five studies, 25 span ≥ 4, and only one
  cluster of 19 cells is single-study — but the two metrics measure different
  things. v4 integrated well enough to annotate, not well enough to claim a
  shared latent space.
- **No benchmark.** The primary embedding was declared by construction; the
  planned Harmony comparison lived only in the deleted v3 object.
- **Leiden flavour matters.** scanpy's `igraph` flavour defaults to CPM and
  over-fragmented to 304 communities at resolution 1.0; `leidenalg` is pinned in
  config for this reason.
- **The fusion axis tracks the cell cycle** (`Ewing_proliferating` +0.058 vs
  `Ewing_sarcoma` +0.003) — see [06](06_discovery.md).
- **48.7% of malignant cells** (14,468 / 29,715) have `disease_state` =
  `not_reported`.
- **O-PDX is one study**, so PDX-versus-patient is confounded in v4 and that
  question is not answerable from this object.

## Differences between v4 as run and the shipped scripts

The scripts in `src/pipeline/` are not a transcription of the v4 session. They
differ deliberately in four ways, each of which fixes something v4 got wrong or
left implicit:

1. The gene floor reads `qc.gene_floor_by_assay`, so the nuclear arm is not held
   to the whole-cell floor (v4 had no nuclear arm, so a single 200 was harmless
   then and would not be now).
2. Every stage writes a manifest; v4's parameters had to be recovered by reading
   `model.pt` after the fact, which is exactly the situation manifests prevent.
3. `02_build_atlas.py` asserts no zero-count cells, marker-panel survival, and a
   complete covariate join — the three failures that cost the earlier rebuilds.
4. Covariates are joined from `data/metadata/ews_covariates_harmonized.csv` on
   `study||sample` with an explicit uniqueness check, rather than merged from
   two schemas.

## Artefact findings (verified 2026-10-18)

Two defects in the released v4 object, both confirmed by reading the file.

**1 · It was stored uncompressed.** A rewrite had dropped
`compression="gzip"`, leaving 2,428,717,616 bytes on disk while the release
notes claimed 701 MB. Recompressed and re-verified at 64,702 × 17,017:
**747,720,744 bytes (713 MiB)**, 30.8% of the uncompressed size. The figure in
`results/reports/EWSAtlas_v4_release.md` now matches the byte count read from
disk.

**2 · It carries no raw counts.** `X` is identical to `layers["lognorm"]` — same
`data` and `indices` arrays across the first 3,000 cells — with non-integral
values capped at 7.822, and `expm1(X)` sums to exactly 10,000 per cell. The
normalization was applied in place and overwrote the count matrix; there is no
`layers["counts"]` and `raw` is `None`. The release notes' claim that `X` holds
"raw counts, int32" is wrong and has been corrected in place.

So v4 as released cannot retrain scVI or feed any count-likelihood model.
Nothing is permanently lost: counts survive in `data/interim/qc_<study>.h5ad`
(integral for Aynaud, Cillo, Goodspeed, Mao; Visser fractional by deposit — see
below), and v5 regenerates them.

This was found by
`tests/test_atlas_object.py::test_atlas_has_counts_layer_and_no_empty_cells` on
its **first** run against the object, which is the clearest available argument
for the invariant suite: the defect had been present through a release-notes
write and a full validation report without being noticed.

**3 · Visser's CEL-seq2 matrices are fractional.** 1.8% of nonzero entries are
non-integral (median nonzero value 1.0, maximum 29,569.7) because CEL-seq2 UMI
counts are collision-corrected. A negative-binomial likelihood is undefined on
fractional input, so `io.read_celseq_plates` now rounds at ingest and keeps the
corrected values in `layers["celseq_corrected"]`; the rounding moves 1.8% of
entries by less than one molecule each. How v4 handled this was not recorded —
another gap the per-stage manifests close.
