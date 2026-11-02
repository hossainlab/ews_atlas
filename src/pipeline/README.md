# Pipeline entrypoints

Each stage reads parameters **only** from `src/config.json`, writes a run
manifest to `logs/`, and never mutates `data/raw`.

| stage | script | does |
|---|---|---|
| 01 | `01_ingest_qc.py` | ingest each deposit, uniform QC |
| 02 | `02_build_atlas.py` | shared gene space, integration-ready object |
| 03 | `03_integrate.py` | HVG selection, training object, scIB benchmark, mixing |
| 03b | `03b_train_scvi.py` | scVI training on GPU (Modal) |
| 04 | `04_annotate.py` | clustering and cell-type assignment |
| 05 | `05_cnv.py` | copy-number inference, frozen malignant set |
| 07 | `07_caf_state.py` | CAF-like malignant state |
| 08 | `08_reference_mapping.py` | leave-one-study-out label transfer |
| 09 | `09_bulk_deconvolution.py` | bulk deconvolution and survival |
| 10 | `10_surface_targets.py` | surface-target selectivity |
| 11 | `11_annotation_revision.py` | immune re-annotation applied to the object |
| 12 | `12_waltner_sensitivity.py` | results without the largest study |
| 13 | `13_cellxgene_export.py` | CELLxGENE-schema export |
| 14 | `14_figures.py` | the publication figure set (runs last) |

There is no stage 06: the figure stage that held that number was replaced by
stage 14, which has to run after the analyses it draws on (ADR-015).

Stages are committed only once executed and verified, so nothing in `src/` is
code that has never run. `run_all.py` executes the stages in its `SHIPPED` list
(currently 01–02); later stages are run individually, in the order above.

## Running

```bash
python src/pipeline/01_ingest_qc.py            # all studies in config.studies
python src/pipeline/01_ingest_qc.py --studies Visser_2023
python src/pipeline/02_build_atlas.py
python src/pipeline/run_all.py                 # shipped stages, in order
python src/pipeline/14_figures.py               # figures, after stages 03-12
```
