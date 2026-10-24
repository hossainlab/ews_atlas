# 09 · Figures

Produced by `src/pipeline/06_figures.py`; panel renderers live in
`ewsatlas.figures` so one panel can be regenerated without redrawing a figure
that already passes. Panels are rendered individually to transparent PNGs at
exact pixel sizes and tiled onto a 12-column grid at 180 mm (double-column)
width, 300 dpi, with bold panel letters stamped at composition time.

**No panel carries a title.** The captions in
`results/reports/figure_captions.md` carry the message, which raises the bar on
in-figure identifiability: every panel must be readable from its axis labels,
direct labels and legend alone.

## The set

| figure | subject |
|---|---|
| `Figure_1_atlas_overview.png` | the resource: embedding, marker evidence, cohort composition |
| `Figure_2_integration_quality.png` | embedding choice by benchmark, and how much residual study structure is composition |
| `Figure_3_malignant_compartment.png` | the statistic that makes copy-number inference usable, and the frozen malignant set |
| `Figure_4_recurrence_and_confounding.png` | what recurs across six studies, and why most of it is capture rather than biology |
| `Figure_5_immune_microenvironment.png` | thirteen immune subsets, correcting two atlas-level calls |
| `Figure_S1_quality_control.png` | barcode fate per study; assay-aware floors and caps |
| `Figure_S2_deposit_specific_qc.png` | the failed ScPCA library; ERCC spike-in share per CEL-seq2 plate |
| `Figure_S3_integration_detail.png` | scVI training curves; the individual scIB metrics |
| `Figure_S4_covariate_coverage.png` | covariate completeness per study |
| `Figure_S5_copy_number_first_attempt.png` | the superseded genome-wide score, kept as the record of what failed |

Panel letters are bold **uppercase** (Cell Press house style), 9 pt derived from
the output dpi, in a reserved band above each grid row; single-panel figures
carry no letter. Stage 06 writes the composed figures straight into `results/pub_figs/` under
submission filenames and copies `figure_captions.md` alongside them; per-panel
intermediates go to `data/interim/figure_panels/`.

## Conventions

- **Dispersed classes are never direct-labelled.** `Unassigned` cells are
  scattered through the malignant mass, so a label at their median names the
  wrong thing; only the classes in `figures.DIRECT` get on-plot text.
- **Colour is threaded.** A cell type keeps one colour across every figure, so
  Figure 1e needs no second legend.
- **Figure 3 reports a negative result.** `cnv_score` ships on the released
  object, so a reader will find it; the figure is what stops them using it.
- **Determinism is checked, not assumed.** Two consecutive runs of the stage
  produce byte-identical PNGs for all seven figures.

No discovery figures exist, because that analysis is blocked — see
[06_discovery.md](06_discovery.md) and the copy-number section of
[05_annotation.md](05_annotation.md). Drawing a fusion-activity axis over a
marker-defined malignant set would look convincing and be circular.
