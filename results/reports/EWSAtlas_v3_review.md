# EWSAtlas v3 — deep review: what is in it, what is missing, and what it can deliver

Scope: audit of `EWA_atlas_v3_integrated.h5ad` and `ews_cohorts.csv` (both in `C:\Users\User\Documents\ews_atlas`),
a systematic search of public archives and literature for Ewing sarcoma (EwS) data not yet included, and a
benchmark against NBAtlas. No files were modified, renamed or deleted.

---

## 1. What the object actually contains (verified)

| Property | Value |
|---|---|
| Cells x genes | 129,276 x 16,843 |
| Studies | 6 single-cell / single-nucleus studies |
| Donor-level ids | 39 (38 EwS donors + 1 urothelial-carcinoma control from Mao 2024) |
| Samples | 68 (`Sample` field) |
| Assay split | 10x snRNA 65,818 (50.9%) / 10x scRNA 56,681 (43.8%) / CEL-seq2 6,777 (5.2%) |
| Malignant cells | 89,817 (`Ewing_sarcoma` 77,145 + `Ewing_proliferating` 12,672) |
| Immune cells | 26,354 (20.4%) |
| Layers | `counts` (integer, nnz 178,215,370, max 1,377), `lognorm`; `X` = counts |
| Embeddings | `X_pca` (50), `X_scvi_v3` (30), `X_harmony_v3` (50), `X_umap`, `X_umap_scvi`, `X_umap_harmony`, `Unintegrated` |
| obs columns | 131 |

The matrix layer is sound. Every defect below is in the metadata and harmonization layer.

---

## 2. Manifest defects (`ews_cohorts.csv`)

Cell counts in the manifest match the object exactly for all six single-cell studies. Patient counts do not.

| Study | Manifest patients | Atlas patients | Atlas samples | Verdict |
|---|---|---|---|---|
| Aynaud 2020 | 10 | 5 | 5 | Wrong — and all five are `PDX-*` ids; no patient tumours from this study are present |
| Cillo 2022 | 7 | 4 | 5 | Wrong |
| Goodspeed 2025 | 5 | 11 | 15 | Wrong (understated) |
| Mao 2024 | 2 | 2 | 4 | Correct (1 EwS + 1 urothelial control) |
| Visser 2023 | 10 | 12 | 26 | Wrong; also conflicts with the publication (18 samples / 11 patients) |
| Waltner 2025 | 5 | 5 | 13 | Correct (7 patient tumours + 6 O-PDX) |

Further manifest issues:

- **Aynaud 2020 is only partially ingested.** The manifest cites GSE130024, which is the Chromium PDX arm only.
  The superseries GSE130025 also contains GSE130023 (PDX C1, 142 samples) and GSE130020 (96 marrow-derived MSC samples).
  The manifest's "O-PDX; Tumor" sample-type claim is not supported by the object.
- **Savola 2011 row understates the cohort.** GEO reports 55 samples for GSE17618 against 32 in the manifest, and the
  companion arm GSE17674 (62 samples, superseries GSE17679) is absent.
- **Volchenboum 2015 row is incomplete.** The companion COG series GSE63155 (46 samples) is absent.
- Both bulk rows are listed in the same sheet as the single-cell rows but are not integrated anywhere; `n_cells_post_qc` is empty for them.

---

## 3. Harmonization defects in the v3 object

**3.1 Two incompatible schemas were merged, not harmonized.**

| Field | Populated for | % of cells |
|---|---|---|
| `Fusion_Status` (`EWSR1-FLI1`) | Waltner only | 45.7 |
| `Fusion_Type` (`EWS-FLI1`) | Visser only | 5.1 |
| `Sex`, `Age`, `Primary_Metastasis`, `is_xenograft`, `disease_timing`, `Tissue_Location` | Waltner only | 50.9 |
| `Body_Site` | Visser only | 4.9 |
| `Tissue` | Mao only | 11.7 |
| `Sort` | Cillo + Visser | 16.3 |

Duplicate-case column pairs coexist (`Sex`/`sex`, `Age`/`age`, `Sample_Type`/`Sample_type`), and four parallel
annotation sets are carried side by side (`cellassign_*`, `singler_*`, `scimilarity_*`, `consensus_*`,
plus `openscpca_*`), all inherited from the ScPCA/OpenScPCA pipeline and therefore present for Waltner only.

**Consequence:** five of six studies — every study except Waltner — have no sex, age or primary/metastasis
annotation. That is 63,458 cells, 49.1% of the atlas. Fusion subtype is known for 50.8% of cells and is recorded
in two different columns with two different spellings. **No covariate in the object is populated for all six studies
except `Assay`, `Platform` and `Study` itself.** This is the single largest barrier to every clinically framed analysis.

**3.2 Xenograft content is undercounted.** `is_xenograft` is `True` for 38,283 Waltner cells and NA elsewhere — but
all 14,702 Aynaud cells are PDX-derived (`PDX-*` patient ids). True xenograft-derived content is at least
52,985 cells (41.0%), not 29.6%. Any statement about "patient tumours" currently includes PDX material.

**3.3 Doublet filtering was never unified.** `is_doublet` is `False` for all 129,276 cells, yet `scDblFinder_class`
labels 7,571 cells `doublet` and is NA for 63,458. Doublets from at least one study are retained in the integrated object.

**3.4 Annotation confidence.** `annotation_confidence_v3` is `review` for 17,453 cells (13.5%); `cell_type_v3`
contains an `Unresolved` class of 1,010 cells, 827 of them from Mao 2024.

**3.5 Composition is confounded with study and with assay.** Waltner contributes 70.5% of all malignant cells
(63,315 / 89,817) and is the only snRNA study (50.9% of the atlas); T cells come overwhelmingly from Cillo and
Goodspeed. Aynaud (PDX) contributes 9 cell types but only 143 T cells. Any compartment proportion computed across
the atlas is currently a statement about study design, not biology.

**3.6 Storage.** Three full copies of a 178.2M-nnz matrix at float64 (`X`, `counts`, `lognorm`) account for the
6.7 GB file. Casting counts to int32/float32 and dropping `X` as a duplicate of `counts` would cut this roughly 3-fold.

---

## 4. Data you are missing

Full table with per-dataset decisions: `ews_candidate_datasets.csv` (26 rows, 18 absent from the atlas).
Searches run: GEO (681 human EwS series; 25 EwS single-cell/spatial series triaged; 50 bulk series with n>=20),
ArrayExpress/BioStudies (0 EwS hits), PubMed, and the ScPCA portal.

### 4.1 Patient material absent from the atlas — highest impact

| Source | What it adds | Status |
|---|---|---|
| **Hong 2022**, doi 10.1007/s00432-022-04073-3 | 39,630 post-QC cells, 3 patients / 4 tumours, **including meningeal EwS** — a site absent from your atlas. Would add ~31% to atlas size | **No accession found** in GEO or ArrayExpress; likely GSA or author request |
| **Chen 2025**, doi 10.1186/s12967-025-06346-0 | 13 immune-sorted bone-tumour samples (EwS, osteosarcoma, giant cell tumour) — would deepen the thin immune compartment | **No accession found** |
| **GSE198896** (Cillo, already partially used) | Matched **PBMC** samples, plus per-sample `primary or recurrent`, `metastatic disease`, `tissue site`, `chemistry` fields in GEO | Public — covariates are retrievable today |
| **GSE130020 / GSE130023** (Aynaud siblings) | Patient-matched bone-marrow MSCs (96) as a cell-of-origin reference; C1 PDX arm (142) | Public |

### 4.2 Reference and control datasets that change what the atlas can claim

| Source | Why it matters |
|---|---|
| **GSE200529** (doi 10.1186/s12885-023-10977-1) | Measures, in sarcoma PDX, how warm vs cold dissociation vs nuclei isolation biases the transcriptome. Your atlas is split 49:51 scRNA:snRNA with the snRNA arm contributing 70.5% of malignant cells. Without this control, compartment-proportion claims are not defensible |
| **GSE332630** (bioRxiv 10.1101/2025.06.18.660457) | The only EwS single-cell **ATAC/Multiome** resource found (22 samples: cell lines, PDX, hMSC, shFLI1). Supplies fusion-regulated GRNs to score atlas cells against |
| **ESCLA** — GSE176339 + GSE212063 (doi 10.1016/j.celrep.2022.111761) | 18 EwS lines with inducible EWSR1-ETS knockdown, multi-omic. The reference EWSR1-ETS-high/low activity axis |
| **GSE291060 / GSE291136** | EWS::FLI1 fluctuation dynamics explicitly tied to prognosis (164 + 222 samples) |
| **GSE88826** | DNA-methylation disease spectrum across 188 EwS tumours — orthogonal axis for malignant-state mapping |

### 4.3 Bulk outcome cohorts absent from the manifest

GSE34620 (117 primary tumours — the largest public EwS expression cohort found), GSE17674 (62),
GSE63155 (46), GSE142162 (79), GSE37372 (78). Together with the two cohorts already listed this takes the
outcome-linked bulk base from 117 to roughly 500 tumours — the substrate for any survival-linked claim.

**Corrected after assembly (Phase 1).** That ~500 projection was wrong on two counts, both found by auditing the
accessions rather than trusting series titles. **GSE63155 is a SubSeries of GSE63157** — all 46 GSM accessions are
already inside it, so it adds nothing. **GSE17618 and GSE17674 are two normalizations of the same 44
hybridizations**, not two arms of one cohort (per-sample Pearson r median 0.994 over 54,675 probes), so the Savola
set does not double. Actual assembled substrate: **408 tumour samples / 343 unique patients across 6 cohorts, of
which only 129 carry usable survival** — GSE63157 (n=85, COG) and GSE17674 (n=44, Savola, whose EFS *events* are
derived from status/EFS<OVS rather than deposited). GSE34620 (117), GSE142162 (79) and GSE37372 (18 usable mRNA
tumours) carry no survival in GEO sample characteristics, and no GSM-keyed outcome table was recoverable from the
companion publications by the routes available here. Any survival-linked claim rests on 129 patients, not 500.

### 4.4 Explicitly out of scope

DSRCT (EWSR1::WT1: GSE263521/263522, GSE267803-807, GSE273438, GSE264741), clear cell sarcoma (EWSR1::ATF1,
GSE180187), intracranial FET-CREB tumours (GSE308071), FET-rearranged myoepithelial tumours (GSE300102),
myxoid liposarcoma (FUS::DDIT3, GSE191132). These are FET-fusion relatives, not EWSR1::FLI1 EwS. Recommend
stating this boundary explicitly in the methods — it is a recurring source of contamination in EwS searches.

### 4.5 Modalities with no EwS patient data found at all

No patient spatial transcriptomics (only GSE272958, an engineered hMSC model), no patient snATAC/Multiome
(GSE332630 is lines/PDX), no TCR/BCR repertoire. These are genuine field-level gaps, not oversights by you.

---

## 5. Benchmark: NBAtlas

NBAtlas (Bonine et al., *Cell Reports* 43:114804, 2024; doi 10.1016/j.celrep.2024.114804) integrated
**7 datasets, 362,991 cells, 61 patients** of neuroblastoma, with the same motivation you state for EwS.

| | NBAtlas | EWSAtlas v3 |
|---|---|---|
| Datasets | 7 | 6 |
| Cells | 362,991 | 129,276 |
| Patients | 61 | 38 EwS donors (14,702 cells PDX-only) |
| Integration | scVI benchmarked against alternatives | scVI + Harmony both stored, neither declared primary |
| sc vs sn comparison | Performed as a dedicated analysis | Not performed (and the confound is larger here) |
| Outcome linkage | Tumour-compartment profiles tied to clinical outcome | Not possible yet — no outcome covariates in the object |
| Reference mapping | Provided for query projection | Embeddings exist; no label-transfer interface |

The gap is not primarily cell count — it is that NBAtlas could stratify by patient and outcome, and v3 cannot.

---

## 6. Impactful updates EWSAtlas is positioned to deliver

Ranked by impact per unit of work. Each names what it additionally requires.

**U1. Harmonize the covariate layer (prerequisite for everything below).**
Collapse the duplicate schemas into one controlled vocabulary (`fusion_partner`, `sex`, `age_years`,
`tumor_site`, `disease_state`, `treatment_timing`, `material_type` with PDX/cell-line flags) and backfill the five
studies that have none. The covariates are recoverable: GSE198896 carries primary/recurrent, metastatic disease
and tissue site per sample in GEO; Goodspeed, Visser and Aynaud report per-patient tables in their publications.
*Requires: GEO sample-characteristic parsing plus manual curation from four papers. No new sequencing.*

**U2. The sc-vs-sn technical comparison, done properly.**
With a 49:51 split and the snRNA arm supplying 70.5% of malignant cells, this is the atlas's central validity
question and NBAtlas's precedent makes it an expected analysis. Quantify compartment bias per assay and report
composition only within assay strata or after explicit correction. *Requires: GSE200529 as the calibration control.*

**U3. A fusion-activity axis across all malignant cells.**
Score every malignant cell for EWSR1-FLI1 activity against an external reference (ESCLA knockdown series,
GSE332630 GRNs, GSE291060 fluctuation model) rather than against atlas-internal programs. This is the one axis
that is genuinely EwS-specific, and the existing `Merge.GEP.EWS` column is an internal construct with no external
calibration. *Requires: ESCLA + GSE332630 signature extraction.*

**U4. Deconvolution bridge to the outcome-linked bulk cohorts (129 patients with usable survival).**
Use harmonized atlas cell states as the reference to deconvolve GSE34620, GSE63155/63157, GSE17618/17674,
GSE142162 and GSE37372, then test state abundance against OS/EFS. This is how NBAtlas tied transcriptional state
to outcome, and it is reachable without any new single-cell data. *Requires: U1, U3, and survival tables from the bulk series.*

**U5. Publish as a queryable reference, not a file.**
A scVI/scANVI reference with a documented label-transfer path turns the atlas from a 6.7 GB download into
infrastructure others cite. *Requires: a frozen primary embedding (choose scVI or Harmony), unified doublet
filtering, and the annotation review backlog (17,453 cells) resolved.*

**U6. The meningeal / extraskeletal site question.**
Hong 2022's meningeal cases would let you ask whether anatomical site maps onto malignant cell state. Currently
`Body_Site` exists for Visser only. *Requires: access to Hong 2022 data (accession not public) and U1.*

---

## 7. Recommended sequence

1. Fix the manifest (`ews_cohorts_reconciled.csv` is the corrected starting point; patient counts now match the object).
2. Rebuild the metadata layer on the existing `counts` matrix — U1. The matrix does not need regenerating.
3. Decide and declare the primary embedding; unify doublet filtering; clear the 17,453-cell review backlog.
4. Add GSE200529 and run U2 before making any composition claim.
5. Add the public reference datasets (ESCLA, GSE332630) and run U3, then U4.
6. Pursue the two access-blocked patient cohorts (Hong 2022, Chen 2025) in parallel — they are the only route to
   materially more patient cells, and neither has a locatable public accession.

## Caveats

- The two highest-value missing cohorts could not be located in any public archive; "seek access" means contacting
  the authors or checking GSA/PRJCA, which this review could not do.
- Controlled-access repositories (EGA, dbGaP) were not searchable with the connectors available here; a controlled-access
  EwS single-cell cohort could exist and would not have appeared in this sweep.
- Per-study QC thresholds and the origin of the 16,843-gene feature space are not recorded in `uns`; whether that space
  is an intersection or a filtered union could not be determined from the object.
- Assumption/residual diagnostics for any downstream model: not assessed.
