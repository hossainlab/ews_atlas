# Reproducing EWSAtlas from a clean checkout

## 1 · Environment

```bash
conda env create -f src/envs/environment-local.yml
conda activate ewsatlas
pip install -e ".[dev]"
pytest -m "not slow"          # 16 invariant tests, no data needed
```

Python is floored at 3.12: anndata ≥ 0.13 requires it, and a 3.11 base resolves
to an older anndata that cannot read the h5ad objects written here.

Exact resolved versions for the reference run are in
`src/envs/requirements-lock-local.txt` (140 distributions) and in every
`logs/run_*.json`.

## 2 · Data

`data/raw` is **not** in version control and is immutable once populated. Fetch:

| target | source |
|---|---|
| `data/raw/geo/GSE130024/GSE130024_RAW.tar` | GEO |
| `data/raw/geo/GSE198896/GSE198896_RAW.tar` | GEO |
| `data/raw/geo/GSE277083/GSE277083_counts_matrix_{primary,blood}.txt.gz` + metadata | GEO |
| `data/raw/geo/GSE250523/GSE250523_RAW.tar` | GEO |
| `data/raw/geo/GSE243347/GSE243347_RAW.tar` + `GSE243347_data_annotation.txt.gz` | GEO |
| `data/raw/scpca/SCPCP000015_ann-data_*.zip` | ScPCA portal (requires accepting their terms and supplying an email — do this yourself) |

Checksums of the exact files used are recorded per stage in `logs/run_*.json`;
compare against those rather than against file sizes.

## 3 · Pipeline

```bash
python src/pipeline/01_ingest_qc.py        # -> data/interim/qc_<study>.h5ad
python src/pipeline/02_build_atlas.py      # -> data/interim/atlas_preintegration.h5ad
# or
python src/pipeline/run_all.py
```

Each stage writes `logs/run_<stage>_<timestamp>.json`. A stage that fails still
writes its manifest with `status: "error"` and the exception, so a partial run is
diagnosable from `logs/` rather than from terminal scrollback.

Later stages run individually, in the order listed in
`src/pipeline/README.md` (03, 03b on GPU, 04, 05, 07–13). The figure set is
built last:

```bash
python src/pipeline/14_figures.py   # -> results/figures/  (see docs/09_figures.md)
```

### Inputs no stage produces yet

The following files are read by a stage but written by none in `src/`. They
came from analyses run outside the committed pipeline. Until each has a
producing stage, a clean checkout cannot regenerate the figures that use them
or the stages listed:

| file | read by | analysis it records |
|---|---|---|
| `results/tables/v5_embedding_benchmark.csv` | 14 (Fig 2A) | aggregate scIB scores per embedding |
| `results/tables/v5_elbo_{sample,study}_{train,val}.csv` | 14 (Fig S3A) | scVI training history |
| `results/tables/v5_depositor_labels.csv` | 03, 14 (Fig 2C, S5D) | harmonized depositor annotations |
| `results/tables/v5_cnv_auc_per_sample.csv` | 14 (Fig 4A, S5B) | first-release genome-wide CNV AUC |
| `results/tables/v6_cnv_arm_auc.csv` | 14 (Fig S5E) | arm-level CNV AUC |
| `results/tables/scpca_library_profile.csv` | 14 (Fig S2A) | ScPCA per-library profile |
| `data/interim/metaprograms_relative.pkl` | 14 (Fig 4C) | per-sample NMF metaprograms |
| `results/tables/v6_metaprogram_occupancy.csv` | 12, 14 (Fig 4E) | metaprogram occupancy by study |
| `results/tables/v6_falsification.csv` | 14 (Fig 4F) | program vs technical-covariate correlations |
| `results/tables/v6_immune_labels.csv` | 11 | immune re-clustering labels |
| `data/interim/immune_subset.h5ad` | 14 (Fig 5) | immune sub-object |
| `data/interim/caf_scored.h5ad` | 07 | CAF-state scoring object |

## 4 · GPU training

Integration runs on Modal:

- image source of truth: `src/envs/ews_integration_gpu.py`
- built image: `im-DdwhygaSGL9V5gKurM7ksW`, spec_sha `83a26f236d5ca71c`
- hardware: A100, 8 CPU, 32 GiB

Two build traps are recorded there: `cell2location` 0.1.5 imports a pyro symbol
that current pyro removed (dropped — no spatial data here), and anndata ≥ 0.13
forces python 3.12 (a 3.11 base fails).

## 5 · What reproduces exactly, and what does not

| stage | reproducibility |
|---|---|
| ingest, QC, feature space, normalization | bit-reproducible under `config.seed` |
| scrublet | bit-reproducible (explicit threshold, seeded) |
| scVI / scANVI on GPU | **to tolerance only** — cuDNN kernel selection and atomic reductions are order-dependent |
| leiden, UMAP | seeded; stable given an identical embedding |
| NNLS deconvolution, Cox models | deterministic |

A re-run that differs should be diagnosed by comparing manifests first: seed,
config hash, input checksums and package versions are all there.

## 6 · Config is the only place to change a parameter

`src/config.json` is read by every stage through `ewsatlas.config.get`, which
raises on a missing key instead of substituting a default. Changing a threshold
in a script is a bug, not a customization.
