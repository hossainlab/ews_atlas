# EWSAtlas v4.0 — release notes and manifest

`data/processed/EWSAtlas_v4.h5ad` | 64,702 cells x 17,017 genes | 747,720,744 bytes (713 MiB), gzip-compressed | 5 studies, 56 samples, 33 Ewing patients

> Size corrected 2026-10-18. An earlier rewrite dropped `compression="gzip"`,
> leaving the file at 2,428,717,616 bytes while these notes still claimed
> 701 MB. The object has been recompressed and re-verified at 64,702 x 17,017;
> the figure above is the byte count read from disk. See
> `docs/08_v4_build_record.md`.

## Scope

EWSR1::FLI1-driven Ewing sarcoma only. **Excluded by design**: DSRCT (EWSR1::WT1), clear cell sarcoma
(EWSR1::ATF1), intracranial FET-CREB tumours, FET-rearranged myoepithelial tumours, myxoid liposarcoma
(FUS::DDIT3). These relatives dominate keyword searches for "Ewing"/"EWSR1" and were filtered out explicitly.
One urothelial-carcinoma patient is retained as the comparator shipped with the bladder study and is flagged
`patient == "UC_Patient"`, `fusion_partner == "other"`.

## Contents

| Slot | Contents |
|---|---|
| `X` | ~~raw counts, int32~~ — **incorrect, see the correction below** |
| `layers["lognorm"]` | library-size normalized (1e4) and log1p |
| `obsm["X_scvi"]` | 30-d scVI latent — **the primary embedding** |
| `obsm["X_umap"]` | UMAP on `X_scvi` |
| `obs` (26 cols) | controlled vocabulary + QC + `cell_type` + `fusion_activity` |
| `uns["atlas"]` | full build provenance |

> **Correction, 2026-10-18 — the released v4 object carries no raw counts.**
> Verified by reading the file: `X` is identical to `layers["lognorm"]` (same
> `data` and `indices` arrays over the first 3,000 cells), its values are
> non-integral with a maximum of 7.822, and `expm1(X)` sums to exactly 10,000
> per cell — i.e. the normalization was applied in place and the count matrix
> was overwritten. There is no `layers["counts"]` and `raw` is `None`.
>
> Consequence: v4 as released cannot be used to retrain scVI or for any other
> count-likelihood model. Nothing is permanently lost — the counts survive in
> `data/interim/qc_<study>.h5ad` (integral for Aynaud, Cillo, Goodspeed and Mao;
> Visser's CEL-seq2 tables are collision-corrected and therefore fractional by
> deposit), and the v5 rebuild regenerates them.
>
> Found by `tests/test_atlas_object.py::test_atlas_has_counts_layer_and_no_empty_cells`
> on its first run against the object.

Controlled vocabulary: `study`, `sample`, `patient`, `fusion_partner`, `sex`, `age_years`, `tumor_site`,
`disease_state`, `treatment_timing`, `material_type`, `assay`, `assay_chemistry`, `sort_strategy`.
Missing values are always the literal string `not_reported` — never blank, never imputed.

Cell-weighted covariate coverage: `material_type` 100%, `disease_state` 77%, `tumor_site` 77%,
`treatment_timing` 67%, `fusion_partner` 55%, `sex` 33%, `age_years` 33%.

## Sources

| Study | Accession | Assay | Cells kept |
|---|---|---|---|
| Aynaud 2020 | GSE130024 | 10x 3', O-PDX | 14,645 |
| Cillo 2022 | GSE198896 | 10x, TIL + PBMC | 14,523 |
| Goodspeed 2025 | GSE277083 | 10x, tumour + CTC | 12,712 |
| Mao 2024 | GSE250523 | 10x, bladder EwS + UC control | 15,664 |
| Visser 2023 | GSE243347 | CEL-seq2 plates | 7,158 |

Built from the deposits themselves. The superseded v3 object was deleted and contributed nothing.

## Decisions a reader needs to know

1. **QC is uniform across studies**: >=200 genes, <=30% mitochondrial reads, per-sample 5-MAD outlier removal on
   log counts / log genes / mito%, and scrublet at an explicit threshold of 0.25 (scikit-image is blocked in this
   environment, so the automatic threshold was unavailable). 89.8% of cells reach the final object (64,702 of 72,023); Visser
   retains 69.0% because its plates carry declared empty wells, 14 of which were removed after QC at assembly.
2. **Gene space** = intersection of the five per-study annotations, restricted to genes detected in >=3 cells in
   at least one study, after stripping `GRCh38_` barnyard prefixes. 17,017 genes. A strict detected-intersection
   would have been 13,923 and the union 71,917.
3. **Visser's pooled plate TM736 is resolved at well level** into patients ES-039 (columns 2-3) and ES-030
   (columns 4-18) using the deposited `selected_wells` assignments, rather than being left as an unassignable
   pseudo-patient.
4. **Plate TM768 (patient ES-048, lymph-node metastasis) is included.** It is present in the deposit and was
   absent from v3.
5. **Annotation is provisional** — marker-based, no CNV inference. See the validation report.

## Outputs alongside the object

- `results/reports/EWSAtlas_v4_validation.md` — integration, annotation, fusion axis, outcome linkage
- `results/tables/v4_bulk_deconvolution.csv` — estimated composition of 408 bulk tumours
- `results/tables/v4_outcome_association.csv` — Cox results per cohort
- `results/tables/v4_fusion_activity_per_cell.csv.gz` — per-cell EWSR1-ETS activity
- `metadata/ews_covariates_harmonized.csv` — per-sample covariates with per-field provenance
- `metadata/ews_candidate_datasets.csv` — the 26 datasets screened, with include/defer/exclude calls
- `src/map_query_to_ewsatlas.py` — reference-mapping recipe

## Open gaps

1. **Waltner / ScPCA SCPCP000015 is missing** (13 samples, ~65,000 nuclei, 5 St Jude patients, 6 O-PDX pairs).
   The portal requires a download token bound to an accepted terms-of-use agreement, which must be created by the
   data user. This is the single largest gap: it removes the only single-nucleus arm and the 100%-covered sex/age
   metadata, holding cell-weighted sex/age coverage at 33%.
2. **sc-vs-sn comparison not assessable** until (1) lands. The GSE200529 calibration control is already built.
3. **Hong 2022** (doi 10.1007/s00432-022-04073-3; 39,630 cells, 3 patients, includes meningeal EwS) and
   **Chen 2025** (doi 10.1186/s12967-025-06346-0; 13 immune-sorted bone tumours) have no locatable public
   accession. They remain the only route found to materially more patient cells.
4. **No patient spatial, snATAC/Multiome or TCR/BCR data exists for EwS** anywhere in the public archives
   searched. GSE332630 is lines/PDX and its deposited matrices are pooled 3-cell-line libraries with no per-cell
   demultiplexing, so it yielded no fusion signature either.
5. **Aynaud's clinical covariates remain `not_reported`** (sex, age, fusion, site, disease state): the paper has
   no PMC record and both the publisher and bioRxiv return 403 to non-browser clients. Cillo's per-patient
   demographics are behind a PMC proof-of-work challenge, which was not bypassed.
6. **`fusion_breakpoint_type` is `not_reported` for every sample** — no source reports it; all deposits record
   the partner gene only.
7. **No integration benchmark and no CNV inference**, both for reasons recorded in the validation report.
