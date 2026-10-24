# EWSAtlas — Harmonized Single-Cell Transcriptomics Atlas of Ewing Sarcoma

## Scientific Goal

Build the first harmonized scRNA-seq atlas of human Ewing Sarcoma (EWS), focused on **tumor cell heterogeneity** — specifically the EWSR1-FLI1 activity spectrum and clonal architecture — using patient tumor samples only. Publication-ready, fully reproducible.

---

## Feasibility: CONDITIONAL GO

Sufficient data for a pilot atlas. One format blocker (Goodspeed 2025) must be resolved in parallel with pilot work.

---

## Datasets

| Dataset | GEO | Patients | Platform | Format | Status |
|---------|-----|----------|----------|--------|--------|
| Visser et al. 2023 | GSE243347 | 11 | CEL-Seq2 | .txt.gz | ✅ Downloaded |
| He et al. 2025 | — | 6 | 10x Chromium | MTX | ✅ Downloaded |
| Goodspeed et al. 2025 | GSE176190 | 7 | 10x (claimed) | .CEL.gz | 🚨 Format blocker |
| Tao et al. 2024 | TBD | ? | 10x | — | Accession needed |
| Aynaud et al. 2020 | TBD | ? | ? | — | Accession needed |

**Excluded (out of scope):** Franzetti 2020 GSE130024 (PDX), ESCLA GSE212063 (cell lines).

**Pilot can start** with He 2025 + Visser 2023 (~17 patients).

---

## Repository Structure

```
EWSAtlas/
├── data/
│   ├── raw/                          # Immutable raw downloads (never modified)
│   ├── processed/                    # Per-dataset QC-passed AnnData (.h5ad)
│   └── integrated/                   # Integrated + annotated atlas files
├── src/
│   └── ewsatlas/
│       ├── io/                       # Dataset-specific loaders
│       ├── qc/                       # MAD-based QC + doublet detection
│       ├── integration/              # Harmony + scVI wrappers
│       ├── annotation/               # CellHint annotation
│       └── utils/                    # Plotting, helpers
├── notebooks/
│   ├── 01_data_loading.ipynb
│   ├── 02_qc_per_dataset.ipynb
│   ├── 03_integration.ipynb
│   └── 04_annotation.ipynb
├── configs/
│   └── params.yaml                   # All tunable parameters
├── PLAN.md                           # This file
└── pyproject.toml
```

---

## Pipeline

### Phase 0 — Data Audit

**Goodspeed format issue**: Downloaded `.CEL.gz` files are Affymetrix microarray format, contradicting the paper's 10x Chromium description. Action: check GEO GSE176190 for separate scRNA-seq sub-accessions. If absent, contact authors.

Find missing GEO accessions for Tao 2024 (PMC10827204) and Aynaud 2020.

---

### Phase 1 — Data Loading

Each dataset loaded into `AnnData` with standardized `.obs` fields:

| Field | Values |
|-------|--------|
| `dataset` | `"Visser2023"`, `"He2025"`, `"Goodspeed2025"` |
| `patient_id` | Globally unique per patient |
| `sample_id` | Per-sample identifier |
| `treatment` | `"naive"`, `"neoadjuvant"`, `"relapsed"` |
| `platform` | `"CELSeq2"`, `"10x_Chromium"` |
| `tissue` | `"primary_tumor"`, `"CTC"` |

Raw counts preserved in `adata.layers["counts"]` before any normalization.

---

### Phase 2 — Automated QC (MAD-Based)

**No manual thresholds.** All cutoffs computed per-sample from data.

#### QC Metrics (via `sc.pp.calculate_qc_metrics`)
- `n_genes_by_counts`
- `total_counts`
- `pct_counts_mt`

#### MAD Outlier Detection
```python
def is_outlier(series, n_mads=5, log_transform=True):
    if log_transform:
        series = np.log1p(series)
    median = np.median(series)
    mad = np.median(np.abs(series - median))
    return np.abs(series - median) > n_mads * mad
```

Filters per sample:
- `n_genes_by_counts`: flag outlier low (empty droplets) + outlier high (potential doublets)
- `total_counts`: flag outlier low/high
- `pct_counts_mt`: flag outlier high only (dying cells)

**Scientific basis**: MAD is robust to non-normal distributions (Lun et al. 2016). Per-sample computation prevents batch-driven false positives. Recommended in Heumos et al. 2023 (*Nature Methods*). Default 5 MADs; configurable in `configs/params.yaml`.

All computed cutoff bounds logged to `data/processed/{dataset}_qc_report.csv`.

#### Doublet Detection
Run **per sample before merging**:
- Primary: `scDblFinder` (best benchmarked; Germain et al. 2022)
- Fallback: `scrublet` (pure Python)

Adds `adata.obs["predicted_doublet"]` and `adata.obs["doublet_score"]`.

---

### Phase 3 — Normalization & Feature Selection

```python
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(
    adata, n_top_genes=3000,
    batch_key="patient_id",    # per-patient HVG selection
    flavor="seurat_v3"
)
```

CEL-Seq2 and 10x use the same normalization; scVI handles remaining depth differences in latent space.

---

### Phase 4 — Integration

#### 4a. Harmony (sanity check, fast)
PCA → Harmony correction on `patient_id`. UMAP inspection to assess biological vs batch mixing.

#### 4b. scVI (final integration)
```python
scvi.model.SCVI.setup_anndata(adata, batch_key="patient_id", layer="counts")
model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
model.train(max_epochs=400, early_stopping=True)
adata.obsm["X_scVI"] = model.get_latent_representation()
```

**Batch key**: `patient_id` — prevents over-correction, preserves inter-patient biological variance.

Integration quality assessed with `scib-metrics`: iLISI, cLISI, NMI, ARI, kBET.

---

### Phase 5 — Cell Annotation

#### Automated cell typing: CellHint
```python
import cellhint
# Reference: HLCA or custom EwS marker panel
# adata.obs["cellhint_label"], adata.obs["cellhint_confidence"]
```

#### Tumor cell identification
- CNV inference: `infercnvpy` (malignant cells have genomic instability)
- EwS marker scoring: NKX2-2, CD99, VCAN, PRKCB, EWSR1

#### EWSR1-FLI1 activity scoring (primary heterogeneity axis)
Franzetti 2020 high/low activity gene signatures → `sc.tl.score_genes()`

Cell states to resolve on tumor cells:
- High EWSR1-FLI1 activity (proliferative)
- Low EWSR1-FLI1 activity (migratory/mesenchymal)
- Intermediate

---

### Phase 6 — Final Atlas Output

`data/integrated/ewsatlas_annotated.h5ad`:
- `adata.layers["counts"]` — raw counts
- `adata.X` — normalized log1p
- `adata.obsm["X_scVI"]` — scVI latent (30D)
- `adata.obsm["X_umap"]` — UMAP coordinates
- `adata.obs` — all metadata + annotations

---

## Configuration (`configs/params.yaml`)

```yaml
qc:
  mad_n_mads: 5
  mt_gene_prefix: "MT-"
  min_cells_per_gene: 3

integration:
  n_hvgs: 3000
  batch_key: "patient_id"
  scvi_n_latent: 30
  scvi_n_layers: 2
  scvi_max_epochs: 400

annotation:
  cellhint_reference: "hlca"
```

---

## Dependencies

Current + additions needed:

```toml
dependencies = [
    "cellhint>=1.0.0",
    "harmonypy>=0.0.9",
    "scanpy>=1.11.5",
    "scvi-tools>=1.4.2",
    "scrublet>=0.2.3",        # doublet detection fallback
    "infercnvpy>=0.4.3",      # CNV-based tumor cell ID
    "scib-metrics>=0.4.1",    # integration benchmarking
    "pyyaml>=6.0",
]
```

---

## Verification Checklist

- [ ] Load He 2025 TN.1 (11,571 cells) via `sc.read_10x_mtx()` — succeeds
- [ ] Compute QC metrics on He 2025 — MAD bounds biologically plausible
- [ ] Run doublet detection on one sample — score distribution looks right
- [ ] Parse one Visser 2023 sample (CEL-Seq2 .txt.gz) — gene × cell matrix correct
- [ ] Check GSE176190 on GEO — confirm whether scRNA-seq matrices exist
- [ ] Harmony integration on He 2025 alone — UMAP shows patient mixing without destroying biology
- [ ] scVI integration — scib-metrics pass threshold (cLISI > 0.6, iLISI > 0.6)
- [ ] CellHint annotation — confidence scores > 0.5 for majority of cells

---

## References

- Lun et al. (2016) *Genome Biology* — MAD-based QC
- Heumos et al. (2023) *Nature Methods* — Best practices scRNA-seq
- Germain et al. (2022) — scDblFinder benchmarking
- Franzetti et al. (2020) — EWSR1-FLI1 activity signatures
- Visser et al. (2023) *Nature Communications* — GSE243347
- He et al. (2025) — MIF-CD74 / macrophage polarization in EWS
- Goodspeed et al. (2025) — Clonal structure, TSPAN8, CTCs in EWS
