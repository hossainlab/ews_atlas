# Pipeline entrypoints

Each stage reads parameters **only** from `src/config.json`, writes a run
manifest to `logs/`, and never mutates `data/raw`.

| stage | script | status |
|---|---|---|
| 01 | `01_ingest_qc.py` | **shipped** — logic validated on all six deposits |
| 02 | `02_build_atlas.py` | **shipped** — includes the zero-count and marker-panel assertions |
| 03 | `03_integrate.py` | lands with the integration phase (scVI arms, Harmony, scIB benchmark) |
| 04 | `04_annotate.py` | lands with the annotation phase (leiden, scANVI, infercnvpy) |
| 05 | `05_discovery.py` | lands with the discovery phase (fusion axis, variance partition) |

Stages 03-05 are deliberately **not** pre-written as untested stubs: each is
committed when it is first executed and verified, so nothing in `src/` is code
that has never run. `run_all.py` executes the shipped stages in order.

## Running

```bash
python src/pipeline/01_ingest_qc.py            # all studies in config.studies
python src/pipeline/01_ingest_qc.py --studies Visser_2023
python src/pipeline/02_build_atlas.py
python src/pipeline/run_all.py                 # shipped stages, in order
```
