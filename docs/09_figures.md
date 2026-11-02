# 09 · Figures

```bash
python src/pipeline/14_figures.py                  # full set, ~1 min once cached
python src/pipeline/14_figures.py --prepare-only   # malignant embedding cache only (~3 min, CPU)
```

The figure stage runs **last**: it reads the atlas object and tables written by
stages 01–12. It writes to `results/figures/`:

- `Figure_*.png` (300 dpi) and `Figure_*.pdf` (vector, editable text; dense
  scatter layers rasterized at 300 dpi), 180 mm wide;
- `figure_stats.json`, the totals quoted in the captions, read at render time;
- `figure_captions.md`, kept by hand.

Per-panel intermediates go to `data/interim/figure_panels/`. The stage deletes
any `Figure_*` file in `results/figures/` it did not write in the same run, so
a renamed figure never sits beside its old name.

Code: panel renderers in `ewsatlas.figures` (one axis in, no title);
grid, rendering and composition in `ewsatlas.layout`. Every geometric and
analytic parameter (width, dpi, gutter, letter size, malignant-embedding
neighbours, programs shown, minimum cells per held-out type and per sample) is
in the `figures` block of `src/config.json`; input and output locations are in
`paths`.

## The set

The order follows the arc shared by published harmonized atlases (GBmap, the
integrated breast-cancer atlas, HLCA, LuCA, NBAtlas): construction →
validation → composition → compartments → reference use → clinical
association. Supplementary figures are numbered in the order the main text
cites them.

| figure | subject |
|---|---|
| `Figure_1_atlas_construction` | study design, the integrated embedding, study and assay, cells per study, marker evidence |
| `Figure_2_integration_and_annotation` | embedding benchmark, residual study structure, depositor agreement, leave-one-study-out label transfer |
| `Figure_3_composition_across_samples` | composition of all 65 samples with covariate strips; by material; study contribution per cell type |
| `Figure_4_malignant_compartment` | the copy-number statistic, malignant-only embedding and program scores, malignant evidence, recurrence and its capture confounding |
| `Figure_5_immune_microenvironment` | thirteen immune subsets, marker evidence, subset share per sample |
| `Figure_6_atlas_as_reference` | surface-target selectivity; bulk deconvolution, Kaplan–Meier and Cox models |
| `Figure_S1_quality_control` | barcode fate per study; assay-aware floors and caps |
| `Figure_S2_deposit_specific_qc` | the failed ScPCA library; ERCC spike-in share per CEL-seq2 plate |
| `Figure_S3_integration_detail` | scVI training curves; individual scIB metrics; embedding by study, assay, sample |
| `Figure_S4_covariate_coverage` | covariate completeness per study |
| `Figure_S5_copy_number_detail` | the superseded genome-wide score, kept as the record of what failed; arm-level scores |
| `Figure_S6_umap_by_sample` | the integrated embedding coloured by sample |
| `Figure_S7_caf_like_state` | the CD73-high fibroblast-like malignant state |

The malignant-only embedding and program scores (Figure 4B–C) are cached at
`paths.malignant_embedding`. Delete the file to recompute it; it is seeded by
`config.seed` and regenerates identically on the same machine and package set.

## Conventions

- **No panel carries a title.** Captions carry the message, so every panel must
  be readable from its axis labels, direct labels and legend alone.
- **Panel letters** are bold uppercase (Cell Press house style), set in
  `figures.letter_pt`. All letters in one grid row share a baseline placed
  just above the highest drawn element in that row, read from the transparent
  panel canvases, so a letter neither floats in empty margin nor overlaps its
  panel. Single-panel figures carry no letter. The bold face is resolved
  through matplotlib's font manager (DejaVu Sans Bold, shipped with
  matplotlib), because PIL cannot find fonts by name on every platform and
  silently falls back to a regular-weight face.
- **Dispersed classes are never direct-labelled.** `Unassigned` cells are
  scattered through the malignant mass, so a label at their median names the
  wrong thing; only the classes in `figures.DIRECT` get on-plot text.
- **Colour is threaded.** A cell type, study or immune subset keeps one colour
  in every figure; immune subsets take the hue family of their parent class.
- **Not-reported covariates are drawn**, in light grey, never omitted.
- **Display versus inference.** The Kaplan–Meier split in Figure 6D (within-
  cohort median) is for display; the test is the continuous, cohort-stratified
  Cox model in Figure 6E.
