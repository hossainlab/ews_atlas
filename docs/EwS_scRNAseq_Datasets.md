# Ewing Sarcoma (EwS) — Public scRNA-seq Dataset Inventory
> Compiled May 2026 | For NBAtlas-style atlas project planning
> Excludes spatial transcriptomics; scRNA-seq only

---

## TIER 1 — Primary Patient Tumor Datasets (scRNA-seq)

These are the core datasets for any EwS atlas project. All profiled from patient-derived primary tumors.

### 1. GSE243347 — Visser et al. 2023
| Field | Details |
|---|---|
| **GEO Accession** | GSE243347 ✅ Confirmed |
| **Authors** | L.L. Visser, M. Bleijs et al. |
| **Journal** | Cancer Research Communications |
| **Year** | 2023 |
| **Patients** | 11 |
| **Samples** | 18 (includes pre- and post-treatment) |
| **Platform** | CEL-Seq2 (plate-based, low throughput) |
| **Estimated cells** | ~5,000–15,000 (plate-based; lower yield than 10x) |
| **Focus** | Immune TME: NK cells, T cells, B cells, dendritic cells, immunosuppressive macrophages; functionally impaired antigen-presenting cells |
| **Key finding** | EwS tumors contain immunosuppressive APCs lacking costimulatory gene expression; T cells show dysfunction markers |
| **Data access** | Public (open GEO) |
| **Batch concern** | CEL-Seq2 vs 10x will require careful integration — different sensitivity and depth |

**Link:** https://doi.org/10.1158/2767-9764.crc-23-0027

---

### 2. GSE176190 — Goodspeed et al. 2025
| Field | Details |
|---|---|
| **GEO Accession** | GSE176190 ✅ Confirmed |
| **Authors** | Andrew Goodspeed, Avery Bodlak, et al. (Univ. Colorado Cancer Center) |
| **Journal** | Clinical Cancer Research |
| **Year** | 2025 (May) |
| **Patients** | 7 (untreated, treatment-naive) |
| **Samples** | 7 primary tumors + matched circulating tumor cells (CTCs) |
| **Platform** | 10x Chromium |
| **Cells** | ~6,228+ profiled; median 13,791 UMIs/cell; median 3,058 genes/cell |
| **Focus** | Transcriptional heterogeneity, clonal substructure, tumor microenvironment immunosuppression, CTC biology, therapeutic targeting |
| **Key finding** | Heterogeneous transcriptional programs in tumor cells correlate with patient outcome; TSPAN8 identified as therapeutic target via CTC profiling |
| **Data access** | Public (open GEO) |

**Link:** https://pubmed.ncbi.nlm.nih.gov/40029262/

---

### 3. He et al. 2025 — MIF-CD74 Axis Study
| Field | Details |
|---|---|
| **GEO Accession** | TBD — check paper Data Availability section |
| **Journal** | Cell Communication and Signaling |
| **Year** | 2025 |
| **Patients** | 6 |
| **Platform** | 10x Chromium (likely) |
| **Focus** | MIF-CD74 signaling axis in tumor-macrophage crosstalk; macrophage polarization; immunosuppression mechanisms |
| **Key finding** | MIF-CD74 pathway drives macrophage-mediated immunosuppression in EwS TME |
| **Validation cohorts** | Used bulk RNA-seq cohorts GSE63157 (n=85) and GSE17679 (n=88) for survival analysis |
| **Data access** | Likely public — confirm via paper |

---

### 4. Tao et al. 2024 — CTC Characterization
| Field | Details |
|---|---|
| **PMC ID** | PMC10827204 |
| **GEO Accession** | TBD — check Data Availability |
| **Year** | 2024 |
| **Material** | Primary tumor biopsies + peripheral blood (CTCs) |
| **Platform** | 10x Chromium |
| **Focus** | Transcriptional heterogeneity of primary tumors and circulating tumor cells; TSPAN8 as therapeutic target |
| **Data access** | Likely public — confirm via paper |

---

## TIER 2 — PDX / Model System Datasets (supporting scRNA-seq)

Useful for mechanistic validation, understanding EWSR1-FLI1 biology, and cross-referencing tumor cell programs to model systems.

### 5. GSE130024 — Franzetti et al. 2020
| Field | Details |
|---|---|
| **GEO Accession** | GSE130024 ✅ Confirmed |
| **Authors** | Franzetti et al. |
| **Journal** | Cell Reports |
| **Year** | 2020 |
| **Models** | 5 PDX lines + 2 MSC-derived models |
| **Platform** | 10x Chromium |
| **Focus** | EWSR1-FLI1 activity heterogeneity — high EWSR1-FLI1 = proliferative state; low EWSR1-FLI1 = migratory/invasive state |
| **Key finding** | First scRNA-seq demonstration that EWSR1-FLI1 activity fluctuation is the primary driver of intra-tumor heterogeneity (directly analogous to MYCN in NBAtlas) |
| **Data access** | Public (open GEO) |
| **Strategic value** | **CRITICAL REFERENCE** — defines the EwS-specific heterogeneity axis equivalent to MYCN in NB |

---

### 6. Aynaud et al. 2020 — Transcriptional Programs (ICA)
| Field | Details |
|---|---|
| **GEO Accession** | Check paper (Cell Reports, 2020) |
| **Journal** | Cell Reports |
| **Year** | 2020 |
| **Material** | Cell lines + model systems |
| **Method** | ICA (Independent Component Analysis) of scRNA-seq |
| **Focus** | Decomposing EwS transcriptional programs; time-resolved EWSR1-FLI1 binding site mapping |
| **Key finding** | Identified conserved transcriptional modules (ICA components) representing distinct cell states in EwS |
| **Strategic value** | Provides the ICA-based transcriptional program framework — analogous to NMF metaprograms in NBAtlas |

**Link:** https://doi.org/10.1016/j.celrep.2020.01.049

---

### 7. GSE212063 — ESCLA (Ewing Sarcoma Cell Line Atlas)
| Field | Details |
|---|---|
| **GEO Accession** | GSE212063 ✅ Confirmed |
| **Year** | 2022 |
| **Models** | 18 EwS cell lines |
| **Platform** | Multi-omic: scRNA-seq + scATAC-seq |
| **Focus** | Comprehensive cell line atlas; multi-omic characterization |
| **Data access** | Public (open GEO) |
| **Strategic value** | Reference for cell line validation; ATAC-seq enables chromatin accessibility analysis |

---

### 8. Khoogar et al. 2022 — EWSR1/FLI1 Downregulation
| Field | Details |
|---|---|
| **PMC ID** | PMC10959445 |
| **Journal** | Cellular Oncology |
| **Year** | 2022 |
| **Material** | EwS cell lines |
| **Focus** | Diverse cellular responses to EWSR1/FLI1 downregulation; dormancy programs; drug resistance |
| **Strategic value** | Characterizes cellular responses when the oncogenic driver is silenced — important for understanding metastatic/dormant states |

---

### 9. eLife Multimodal 2025 (Preprint/Reviewed)
| Field | Details |
|---|---|
| **Source** | eLife Reviewed Preprint / bioRxiv (2025) |
| **DOI** | 10.1101/2025.06.18.660457 |
| **Material** | Cell lines + primary patient tumors |
| **Platform** | scRNA + scATAC co-assay (10x Multiome or similar) |
| **Focus** | Fusion-regulated transcriptional programs; EWS::FLI1 regulatory modules; intratumoral heterogeneity |
| **Key finding** | Distinct regulatory modules are variably enriched across and within patient tumors; TGF-β modifies module usage |
| **Strategic value** | Most comprehensive multi-omic characterization to date; includes patient tumor data |

---

## TIER 3 — Bulk RNA-seq Validation Cohorts (NOT scRNA-seq)

Not suitable for atlas integration, but essential for survival analysis, gene signature validation, and deconvolution.

| Dataset | Patients | Year | Notes |
|---|---|---|---|
| GSE63157 | 85 | — | Large primary tumor cohort; commonly used for EwS survival analysis |
| GSE17679 | 88 | — | Another large patient cohort; frequently used alongside GSE63157 |
| ICGC EwS (EGAS00001000855) | ~100+ | 2014 | Genomic landscape study; WGS/WES + RNA; European Genome-phenome Archive (controlled access) |

---

## Summary: EwS scRNA-seq Landscape vs NBAtlas

| Metric | NBAtlas (NB) | EwS Public Data |
|---|---|---|
| Confirmed patient scRNA-seq datasets | 7 | ~3–4 |
| Total patients with scRNA-seq | 61 | ~24–30 (est.) |
| Total cells | 362,991 | ~20,000–40,000 (est.) |
| Primary tumor platform consistency | Mixed (dominated by 10x) | Mixed (10x + CEL-Seq2) |
| Existing harmonized atlas | ✅ Yes (NBAtlas) | ❌ None published |
| Controlled-access data | Yes (EGA) | Yes (EGA ICGC) |

**Gap assessment:** EwS has roughly 1/3 the patient numbers of NBAtlas and no existing atlas — making this both the challenge and the opportunity.

---

## Strategic Recommendations for EwS Atlas (No In-House Data)

### Option A — Public Data Only (Feasible Now)
- **Combine GSE243347 + GSE176190 + He et al. 2025** → ~24 patients, ~20,000–40,000 cells
- Integrate with scVI; use CEL-Seq2 as a batch covariate
- Scope: "First integrated single-cell landscape of EwS TME" — modest but publishable as a focused resource paper (not Cell Reports scale, but suitable for Cancer Research Communications, Frontiers in Oncology, or similar)

### Option B — Data Request Strategy (Recommended)
Contact authors directly for raw data not yet deposited:
1. **Tao et al. 2024 (PMC10827204)** — email corresponding author
2. **He et al. 2025** — corresponding author (Cell Commun Signal)
3. **EGA ICGC (EGAS00001000855)** — data access application; takes 4–8 weeks

### Option C — Collaboration (Highest Impact)
- Contact **Visser/Bleijs lab** (Utrecht) — largest patient cohort, active in EwS scRNA-seq
- Contact **Goodspeed lab** (Colorado) — 10x data, treatment-naive patients
- Propose: you lead computational integration; they provide validation/co-authorship

### Option D — Wait + Prepare (12–18 months)
- Prepare the full pipeline (QC, scVI integration, annotation, CNV, NMF) using NBAtlas as template
- Monitor GEO/EGA for new deposits (SARC028 trial data likely incoming)
- Position to rapidly publish the first atlas when sufficient data exists

---

## Key Papers to Read (Priority Order)

1. Goodspeed et al. 2025 (Clin Cancer Res) — most current patient scRNA-seq
2. Visser et al. 2023 (Cancer Res Commun) — largest patient cohort, immune TME
3. Franzetti et al. 2020 (Cell Reports) — EWSR1-FLI1 heterogeneity biology (your MYCN analog)
4. Aynaud et al. 2020 (Cell Reports) — ICA transcriptional programs framework
5. eLife Multimodal 2025 — most recent multi-omic view



PubMed Metadata — EwS scRNA-seq Dataset Papers                                   
                                                                                   
  1. Visser et al. 2023
                                                                                   
  ┌───────┬───────────────────────────────┐                                        
  │ Field │                                Value                                │  
  ├───────┼─────────────────────────────────────────────────────────────────────┤
  │ PMID  │ 37823774                                                            │
  ├─────────┼───────────────────────────────────────────────────────────────────┤
  │ PMCID   │ PMC10595530                                                       │  
  ├─────────┼───────────────────────────────────────────────────────────────────┤
  │ DOI     │ 10.1158/2767-9764.CRC-23-0027                                     │  
  ├─────────┼───────────────────────────────────────────────────────────────────┤  
  │ Title     │ Ewing Sarcoma Single-cell Transcriptome Analysis Reveals        │
  │           │ Functionally Impaired Antigen-presenting Cells                  │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Authors   │ Visser LL, Bleijs M, Margaritis T, van de Wetering M, Holstege  │  
  │           │ FCP, Clevers H                                                  │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Journal   │ Cancer Res Commun                                               │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ Published │ 2023 Oct 24; 3(10):2158–2169                                    │
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │           │ Sarcoma, Ewing/genetics; Single-Cell Gene Expression Analysis;  │
  │ MeSH      │ Antigen-Presenting Cells/metabolism; Tumor                      │  
  │           │ Microenvironment/genetics                                       │
  └───────────┴─────────────────────────────────────────────────────────────────┘  
           
  ---                                                                              
  2. Goodspeed et al. 2025
                                                                                   
  ┌───────┬───────────────────────────────┐
  │ Field │                                Value                                │
  ├───────┼─────────────────────────────────────────────────────────────────────┤
  │ PMID  │ 40029262                                                            │
  ├─────────┼───────────────────────────────────────────────────────────────────┤
  │ PMCID   │ PMC12081191                                                       │
  ├─────────┼───────────────────────────────────────────────────────────────────┤  
  │ DOI     │ 10.1158/1078-0432.CCR-24-2040                                     │
  ├─────────┼───────────────────────────────────────────────────────────────────┤  
  │ Title     │ Single-Cell RNA Sequencing of Ewing Sarcoma Tumors Demonstrates │
  │           │  Transcriptional Heterogeneity and Clonal Evolution             │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Authors   │ Goodspeed A, Bodlak A, Duffy AB, Nelson-Taylor S, Oike N,       │  
  │           │ Porfilio T, … Hayashi M                                         │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Journal   │ Clin Cancer Res                                                 │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ Published │ 2025 May 15; 31(10):2010–2023                                   │
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │           │ Sarcoma, Ewing/genetics/pathology/mortality; Clonal             │
  │ MeSH      │ Evolution/genetics; Genetic Heterogeneity; Neoplastic Cells,    │  
  │           │ Circulating; Tumor Microenvironment/genetics/immunology         │
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ Note      │ PMID 38293103 = bioRxiv preprint of same paper (PMC10827204)    │
  └───────────┴─────────────────────────────────────────────────────────────────┘  
   
  ---                                                                              
  3. He et al. 2025 (MIF-CD74)
                                                                                   
  ┌───────┬────────────────────────────┐
  │ Field │                                Value                                │  
  ├───────┼─────────────────────────────────────────────────────────────────────┤
  │ PMID  │ 39800691                                                            │
  ├─────────┼───────────────────────────────────────────────────────────────────┤
  │ PMCID   │ PMC11727170                                                       │  
  ├─────────┼───────────────────────────────────────────────────────────────────┤
  │ DOI     │ 10.1186/s12964-024-02020-y                                        │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ Title     │ Integrative analysis of Ewing's sarcoma reveals that the        │
  │           │ MIF-CD74 axis is a target for immunotherapy                     │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Authors   │ He F, Xu J, Zeng F, Wang B, Yang Y, Xu J, Sun X, Ren T, Tang X  │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ Journal   │ Cell Commun Signal                                              │
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ Published │ 2025 Jan 13; 23(1):23                                           │
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ MeSH      │ Sarcoma, Ewing/therapy/immunology; MIF/metabolism;              │
  │           │ Macrophages/immunology; CD8+ T cells; Immunotherapy             │  
  └───────────┴─────────────────────────────────────────────────────────────────┘
                                                                                   
  ---      
  4. Franzetti et al. 2017 (GSE130024 attribution note)
                                                                                   
  ┌───────┬──────────────────────┐
  │ Field │                                Value                                │  
  ├───────┼─────────────────────────────────────────────────────────────────────┤
  │ PMID  │ 28135250                                                            │
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Journal   │ Cell Commun Signal                                              │
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Published │ 2025 Jan 13; 23(1):23                                           │
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ MeSH      │ Sarcoma, Ewing/therapy/immunology; MIF/metabolism;              │
  │           │ Macrophages/immunology; CD8+ T cells; Immunotherapy             │
  └───────────┴─────────────────────────────────────────────────────────────────┘

  ---
  4. Franzetti et al. 2017 (GSE130024 attribution note)

  ┌───────────┬─────────────────────────────────────────────────────────────────┐
  │   Field   │                              Value                              │
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ PMID      │ 28135250                                                        │
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ PMCID     │ PMC5541267                                                      │
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ DOI       │ 10.1038/onc.2016.498                                            │
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Title     │ Cell-to-cell heterogeneity of EWSR1-FLI1 activity determines    │
  │           │ proliferation/migration choices in Ewing sarcoma cells          │
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Authors   │ Franzetti GA, Laud-Duval K, van der Ent W, Brisac A, Irondelle  │
  │           │ M, … Delattre O                                                 │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Journal   │ Oncogene                                                        │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Published │ 2017 Jun 22; 36(25):3505–3514                                   │
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ MeSH      │ EWSR1-FLI1 fusion; Cell Movement; Cell Proliferation; Neoplasm  │
  │           │ Invasiveness; Phenotypic Plasticity                             │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │           │ "Franzetti et al. 2020 Cell Reports" not found in PubMed —      │  
  │ Note      │ dataset doc may be misattributed. This 2017 Oncogene paper      │  
  │           │ covers the same biology (EWSR1-FLI1^high proliferative vs ^low  │
  │           │ migratory). Verify GEO GSE130024 submitter against this PMID.   │  
  └───────────┴─────────────────────────────────────────────────────────────────┘

  ---
  5. Aynaud et al. 2020
                                                                                   
  ┌───────────┬─────────────────────────────────────────────────────────────────┐
  │   Field   │                              Value                              │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ PMID      │ 32049009                                                        │
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ PMCID     │ not indexed (no PMC)                                            │
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ DOI       │ 10.1016/j.celrep.2020.01.049                                    │
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ Title     │ Transcriptional Programs Define Intratumoral Heterogeneity of   │
  │           │ Ewing Sarcoma at Single-Cell Resolution                         │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Authors   │ Aynaud MM, Mirabeau O, Gruel N, Grossetête S, Boeva V, Durand   │  
  │           │ S, … Delattre O, Zinovyev A                                     │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Journal   │ Cell Rep                                                        │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Published │ 2020 Feb 11; 30(6):1767–1779.e6                                 │
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ MeSH      │ Sarcoma, Ewing/genetics; RNA-Binding Protein EWS/metabolism;    │
  │           │ Signal Transduction; Transcription, Genetic                     │  
  └───────────┴─────────────────────────────────────────────────────────────────┘
                                                                                   
  ---      
  6. Khoogar et al. 2022
                                                                                   
  ┌───────────┬─────────────────────────────────────────────────────────────────┐
  │   Field   │                              Value                              │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ PMID      │ 34997546                                                        │
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ PMCID     │ PMC10959445                                                     │
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ DOI       │ 10.1007/s13402-021-00640-x                                      │
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ Title     │ Single-cell RNA profiling identifies diverse cellular responses │
  │           │  to EWSR1/FLI1 downregulation in Ewing sarcoma cells            │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Authors   │ Khoogar R, Li F, Chen Y, Ignatius M, Lawlor ER, Kitagawa K,     │  
  │           │ Huang TH, Phelps DA, Houghton PJ                                │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Journal   │ Cell Oncol (Dordr)                                              │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ Published │ 2022 Feb; 45(1):19–40                                           │
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ MeSH      │ EWSR1/FLI1 oncogene; Down-Regulation; Dormancy; Sarcoma,        │
  │           │ Ewing/genetics/metabolism                                       │  
  └───────────┴─────────────────────────────────────────────────────────────────┘
                                                                                   
  ---      
  7. Orth et al. 2022 (ESCLA)
                                                                                   
  ┌───────────┬─────────────────────────────────────────────────────────────────┐
  │   Field   │                              Value                              │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ PMID      │ 36476851                                                        │
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ PMCID     │ PMC10333306                                                     │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ DOI       │ 10.1016/j.celrep.2022.111761                                    │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ Title     │ Systematic multi-omics cell line profiling uncovers principles  │
  │           │ of Ewing sarcoma fusion oncogene-mediated gene regulation       │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Authors   │ Orth MF, Surdez D, Faehling T, Ehlers AC, Marchetto A,          │  
  │           │ Grossetête S, … Grünewald TGP                                   │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Journal   │ Cell Rep                                                        │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Published │ 2022 Dec 6; 41(10):111761                                       │
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ MeSH      │ Sarcoma, Ewing/genetics; Multiomics; Oncogenes; Transcription   │
  │           │ Factors                                                         │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Note      │ This IS the ESCLA paper (18 cell lines, GEO GSE212063 likely    │  
  │           │ linked here)                                                    │  
  └───────────┴─────────────────────────────────────────────────────────────────┘
                                                                                   
  ---      
  8. Waltner et al. 2025 (eLife multimodal, preprint)
                                                                                   
  ┌───────────┬─────────────────────────────────────────────────────────────────┐
  │   Field   │                              Value                              │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ PMID      │ 40667384                                                        │
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ PMCID     │ PMC12262694                                                     │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ DOI       │ 10.1101/2025.06.18.660457                                       │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ Title     │ Multimodal single-cell analyses reveal distinct                 │
  │           │ fusion-regulated transcriptional programs in Ewing sarcoma      │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Authors   │ Waltner OG, Apfelbaum AA, Wrenn ED, Bhise SS, Kanaan SB,        │  
  │           │ Gladden RG, … Furlan SN                                         │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤
  │ Journal   │ bioRxiv (preprint)                                              │  
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ Published │ 2025 Sep 3                                                      │
  ├───────────┼─────────────────────────────────────────────────────────────────┤  
  │ MeSH      │ Not yet indexed (preprint)                                      │
  └───────────┴─────────────────────────────────────────────────────────────────┘  
   
  ---                                                                              
  Key findings:
  - Tao 2024 (PMC10827204) = Goodspeed 2025 preprint — same paper, now published   
  - "Franzetti 2020 Cell Reports" not in PubMed; the EWSR1-FLI1 heterogeneity work
  is Franzetti 2017 Oncogene (PMID 28135250) — dataset doc has wrong year/journal. 
  Verify GSE130024 submitter.                                                      
  - ESCLA = Orth et al. 2022 (not originally named Franzetti)                      
  - 8/9 papers confirmed with full metadata; all PMIDs verified 