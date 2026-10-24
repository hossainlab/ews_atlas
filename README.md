# EWSAtlas — project folder

Harmonized single-cell / single-nucleus reference for Ewing sarcoma (EWSR1::FLI1-driven),
built in the spirit of NBAtlas (Bonine et al., *Cell Reports* 43:114804, doi 10.1016/j.celrep.2024.114804).

## Scope boundary

EWSR1::FLI1 Ewing sarcoma only. Explicitly **out of scope**: DSRCT (EWSR1::WT1), clear cell sarcoma
(EWSR1::ATF1), intracranial FET-CREB tumours, FET-rearranged myoepithelial tumours, and myxoid liposarcoma
(FUS::DDIT3). These relatives dominate keyword searches for "Ewing" and "EWSR1" and must be filtered out.

## Layout

```
ews_atlas/
├── EWA_atlas_v3_integrated.h5ad   # INPUT, do not modify — 6.7 GB, 129,276 cells x 16,843 genes
├── ews_cohorts.csv                # original manifest (patient counts wrong for 4/6 studies)
├── NBAtlas.pdf                    # comparator paper
├── README.md                      # this file
├── data/
│   ├── raw/geo/                   # per-study downloads from GEO
│   ├── raw/scpca/                 # ScPCA / OpenScPCA downloads
│   ├── external/bulk_cohorts/     # bulk outcome cohorts (GSE34620, GSE63155/63157, ...)
│   ├── external/references/       # ESCLA, GSE332630, GSE200529 reference series
│   ├── interim/                   # intermediate objects (obs_v3.pkl, var_v3.pkl, h5ad_tree.txt)
│   └── processed/                 # EWSAtlas_v4.h5ad lands here
├── data/metadata/
│   ├── ews_cohorts_reconciled.csv     # manifest corrected against the v3 object
│   ├── ews_candidate_datasets.csv     # 26 datasets, include/defer/exclude decisions
│   └── curation_sources/              # raw GEO search dumps backing those decisions
├── results/
│   ├── reports/EWSAtlas_v3_review.md  # the deep review of v3
│   ├── tables/                        # coverage, composition, reconciliation tables
│   └── figures/
├── src/                           # reusable scripts (reference-mapping recipe)
├── docs/
└── logs/
```

## Path contract

Later phases write to fixed paths; nothing overwrites the v3 input.

| Output | Path |
|---|---|
| Harmonized covariate table | `data/metadata/ews_covariates_harmonized.csv` |
| Rebuilt atlas | `data/processed/EWSAtlas_v4.h5ad` |
| Validation report | `results/reports/EWSAtlas_v4_validation.md` |
| Outcome-linkage results | `results/tables/` |

## Environments

| Environment | Purpose | Key packages |
|---|---|---|
| `ews` | QC, metadata harmonization, annotation, plotting | scanpy 1.12.4, anndata 0.13.4, h5py 3.16.0, harmonypy 0.2.0, leidenalg 0.12.0 (Python 3.13.15) |
| `ews-scvi` | scVI / scANVI integration and label transfer | scvi-tools 1.5.1, torch 2.14.1+cpu, lightning 2.6.6, scanpy 1.12.4, anndata 0.13.4 (Python 3.12) |

scVI trains on **CPU** (torch reports `cuda: False`, 6 threads) — this machine has 12 cores, no GPU, and no remote
compute target is configured. Budget training time accordingly, or add a compute target via Customize → Compute.

Note on `ews-scvi`: scvi-tools resolves from a venv overlay inside the session workspace, not from the conda
prefix (an initial `scvi-tools` resolve pulled 0.12.0, which was incompatible with anndata 0.13 and had to be
removed from the conda layer to stop it shadowing 1.5.1). If an import of `scvi` fails in a later session,
reprovision with `pip install "scvi-tools>=1.2"` run from inside a python cell.

## Disk

82.0 GiB free on `C:` as of folder preparation. v3 is 6.7 GB; v4 is expected near 2 GB because counts will be
stored as int32 with a single normalized layer instead of v3's three float64 copies of a 178.2M-nnz matrix.

## State of v3 (why v4 is being built)

The v3 counts matrix is sound; its metadata layer is not. Two incompatible schemas were merged, so five of six
studies (63,458 cells, 49.1%) carry no sex, age or primary/metastasis annotation; `is_xenograft` misses the
14,702 PDX-derived Aynaud cells; 7,571 scDblFinder-labelled doublets were retained; and 17,453 cells are
annotation-review-flagged. Full detail and the ranked update roadmap: `results/reports/EWSAtlas_v3_review.md`.

## Repository layout

One home per kind of file. If a file could plausibly live in two places, it
lives in the one listed here and nowhere else.

```
ews_atlas/
├── src/
│   ├── config.json          single source of every parameter; stages read nothing else
│   ├── ewsatlas/            importable library (io, qc, genespace, figures, compat, provenance)
│   ├── pipeline/            numbered stages 01-06, each runnable alone
│   └── envs/                environment and lockfiles, plus the Modal GPU script
├── tests/                   28 fast + 4 slow invariant tests
├── docs/                    numbered documentation, decisions.md, reproduce.md
│   └── references/          source papers (NBAtlas)
├── data/
│   ├── raw/                 downloads, never written to by any stage
│   ├── metadata/            harmonized covariates and their curation sources
│   ├── external/            reference signatures and bulk cohorts
│   ├── interim/             regenerable intermediates (QC objects, embeddings, panels)
│   └── processed/           EWSAtlas_v5.h5ad, the release object
├── results/
│   ├── pub_figs/            THE figures - submission names, plus figure_captions.md
│   ├── tables/              numeric results, one CSV/JSON per analysis
│   └── reports/             written reports and release notes
└── logs/                    one run manifest per stage execution
```

**Rules that keep it this way.**

- Composed figures are written by stage 06 directly into `results/pub_figs`
  under submission names. There is no second copy under another naming scheme.
- Per-panel PNGs are intermediates and go to `data/interim/figure_panels`;
  nothing in `results/` is regenerable scratch.
- `data/interim` holds only what is expensive to recreate. The scVI latents in
  `data/interim/embeddings` are the output of a GPU job and are the reason the
  atlas can be rebuilt without re-running it.
- Nothing is written to the project root.
