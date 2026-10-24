# 03 · Feature space

Implemented in `src/ewsatlas/genespace.py`, driven by `config.gene_space`.

## Rule

A gene is kept when

1. **every** study's annotation contains it (`intersection of annotations`), and
2. at least `detected_in_n_studies` (= 1) study detected it in at least
   `min_cells_detected` (= 3) cells.

## Why not the alternatives

**Strict detected-intersection** (detected in *every* study) penalises the
shallowest arm: a gene that is real in five cohorts gets dropped because the
sparsest one missed it. In v4 this cut the space from 17,017 to 13,923 features
— about 3,100 genes discarded for a reason that is technical, not biological.

**Union of annotations** admits genes absent from a study's reference, which
then appear as structural zeros. An integration model has no way to distinguish
a structural zero from a biological zero and will read the pattern as biology —
it is a direct route to study-driven latent structure. The union here is 71,917
features.

**Intersection of annotations ∩ any-detection** keeps only genes every study
could in principle have measured, then removes those nobody actually saw. In v4:
17,017 selected, of which 16,959 met the detection criterion.

## Symbol harmonization happens first

Every reader applies `io.harmonize_symbols`: upper-case, drop any `ENSG…__`
prefix, strip genome-build prefixes (`^(GRCH38|HG19|MM10|GRCM38)[_-]`). Without
this the Aynaud intersection collapses — see [01](01_data_sources.md).

## Marker panel is a hard gate

`genespace.MARKER_PANEL` (FLI1, EWSR1, NKX2-2, NR0B1, CD99, PTPRC, COL1A1,
PECAM1, CD3E, LYZ, MKI67, TOP2A) must survive feature selection.
`02_build_atlas.py` raises if any is absent, because their absence means symbol
harmonization broke, not that the biology changed.

## Highly variable genes

HVG selection for integration is separate from the atlas feature space: 5,000
genes, `seurat_v3` flavour, computed on counts with the study as the batch key
so that no single deposit dominates selection. NBAtlas used 5,000 HVGs with the
same flavour; v4 used only 2,000, which is one of the reasons it is being
rebuilt (see [04](04_integration.md)).
