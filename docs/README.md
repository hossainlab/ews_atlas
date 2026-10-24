# EWSAtlas documentation

A harmonized single-cell / single-nucleus reference for **EWSR1::FLI1 Ewing
sarcoma**, assembled from every public Ewing sarcoma single-cell deposit we
could locate, and used to ask one discovery question (see
[06_discovery.md](06_discovery.md)).

Read in order:

1. [00_project_overview.md](00_project_overview.md) — scope, inclusion and exclusion criteria, what the atlas is and is not
2. [01_data_sources.md](01_data_sources.md) — the six deposits, their on-disk shapes, and the traps in each
3. [02_quality_control.md](02_quality_control.md) — the QC policy and why the gene floor is assay-aware
4. [03_gene_space.md](03_gene_space.md) — feature-space construction
5. [04_integration.md](04_integration.md) — scVI/Harmony configuration and how the primary embedding is chosen
6. [05_annotation.md](05_annotation.md) — cell-type labels and malignant-cell calling
7. [06_discovery.md](06_discovery.md) — the biological question, and what the data can and cannot answer
8. [07_bulk_validation.md](07_bulk_validation.md) — external outcome substrate
9. [reproduce.md](reproduce.md) — reproducing the atlas from a clean checkout
10. [08_v4_build_record.md](08_v4_build_record.md) — what the superseded v4 build actually ran, and the defects found in it
11. [09_figures.md](09_figures.md) — what each figure claims, and the conventions behind them
12. [decisions.md](decisions.md) — every non-obvious choice, with its scientific rationale

## Conventions used throughout

- Parameters live in `src/config.json`. No script defines a default; a missing
  key raises rather than silently substituting one.
- `data/raw` is immutable. Any stage that would write there raises
  `PermissionError` (`ewsatlas.provenance.assert_raw_readonly`).
- Every stage writes `logs/run_<stage>_<timestamp>.json` recording the seed,
  the full resolved config, each input's sha256, 15 package versions, and the
  host/device.
- Missing covariates are the string `not_reported`, never a blank or `NaN`,
  so "unknown" and "absent" stay distinguishable.
- A page that describes a stage not yet executed says so in its status line.
- `pytest -m "not slow"` must pass at all times (17 tests, no data needed).
  `pytest -m slow` runs the same invariants against the built object and
  **currently fails one test on v4 by design** — v4 carries no count matrix
  (see [08_v4_build_record.md](08_v4_build_record.md)). It is left failing
  rather than skipped, because that is the defect v5 must fix.
