# Is EWSAtlas novel? A literature assessment

Scope: PubMed searches for Ewing sarcoma + atlas/integrated/harmonized/meta-analysis + single-cell (11 records),
sarcoma single-cell atlases (7), pediatric single-cell atlases (13), plus the 3CA and ScPCA resources and the
NBAtlas methods. Assessed against the three stated goals.

---

## Verdict in one line

**Goal 1 (harmonized resource): novel — no integrated Ewing sarcoma single-cell atlas exists, and EwS is now the
conspicuous gap among bone/soft-tissue sarcomas. Goal 2 (biological discovery): mostly occupied — the obvious
immune and fusion-heterogeneity claims are already published, so the discovery hook has to be chosen deliberately
rather than inherited from NBAtlas.**

---

## 1. Resource novelty — the gap is real

No multi-study integrated EwS single-cell reference has been published. Every EwS single-cell paper to date
analyses its own cohort, and the "integrative" EwS papers integrate *modalities* (scRNA + bulk + WGCNA), not
*cohorts*.

Meanwhile sibling sarcomas already have theirs:

| Disease | Resource | Scale |
|---|---|---|
| Osteosarcoma | Integrated OS cell atlas (BMC Cancer 2025, doi 10.1186/s12885-025-13714-y) | >110,000 cells, 17 samples |
| Chondrosarcoma | Single-cell atlas of conventional central chondrosarcoma (Commun Biol 2024) | — |
| Undifferentiated pleomorphic sarcoma | Single-cell landscape (Oncogene 2024) | — |
| Rhabdomyosarcoma | Immune suppression and outcome-predictive cell states (Nat Commun 2023) | — |
| Infantile fibrosarcoma | Malignant cells and microenvironment (Transl Pediatr 2024) | — |
| **Ewing sarcoma** | **none** | — |

So the format is established and publishable in this exact disease family, and EwS is the hole. That is the
strongest single argument for the project.

### The two resources a reviewer will raise

- **ScPCA / Single-cell Pediatric Cancer Atlas** (Cell Genomics 2026). Provides *uniformly processed* pediatric
  tumour data with cell-type annotations, including the EwS project SCPCP000015. It is per-sample processing and
  annotation, not cross-study integration, and it covers one EwS cohort.
- **3CA, Curated Cancer Cell Atlas** (Nat Cancer 2025): 124 datasets, 2,836 samples, 5,658,705 cells across 40+
  cancer types, and it does carry Ewing sarcoma. But its Sarcoma category holds only 2 studies / 26 samples /
  80,682 cells, and 3CA's design is per-dataset curation plus cross-dataset metaprograms — not a single integrated
  disease reference with harmonized clinical covariates.

**The defensible positioning** is therefore: cross-study *integration* (one latent space, one annotation
pipeline, one controlled vocabulary) plus a query-mapping interface — not "we collected the data", which ScPCA
and 3CA already did.

### Where the scale sits

| Atlas | Cells | Patients |
|---|---|---|
| NBAtlas | 362,991 | 61 |
| Osteosarcoma atlas | ~110,000 | 17 samples |
| **EWSAtlas v4.0 (current)** | **64,702** | **33** |
| EWSAtlas + Waltner (projected) | ~99,000-101,000 | 38 |

Mid-pack, not a scale leader. EwS is rarer than neuroblastoma and the public data simply does not exist — the
completeness review (`metadata/ews_integration_completeness.csv`) establishes that the integrable universe is
closed. That is an argument to make explicitly in the paper rather than a weakness to hide: the atlas is
*exhaustive* of public EwS single-cell data even though it is not large.

---

## 2. Discovery novelty — most of the obvious claims are taken

Already published for EwS, with the claim and its owner:

| Claim | Source |
|---|---|
| EWSR1-FLI1 activity heterogeneity defines intratumoral transcriptional programs (ICA) | Aynaud 2020, doi 10.1016/j.celrep.2020.01.049 |
| CAF-like malignant cells (CD73+) remodel the TME | Clin Cancer Res 2023, doi 10.1158/1078-0432.CCR-23-1111 |
| Functionally impaired antigen-presenting cells | Visser 2023, doi 10.1158/2767-9764.CRC-23-0027 |
| Distinct immune signatures and intercellular communication vs osteosarcoma | Cillo 2022, doi 10.1158/1078-0432.CCR-22-1471 |
| MIF-CD74 axis as immunotherapy target; SLC40A1/C1QA macrophages poor prognosis, CD8 T cells good prognosis | Cell Commun Signal 2025, doi 10.1186/s12964-024-02020-y |
| T-cell-associated prognostic signature from scRNA + WGCNA + bulk | Front Immunol 2025, doi 10.3389/fimmu.2025.1586544 |
| Neural-crest origin and developmental trajectory (incl. meningeal EwS) | Hong 2022, doi 10.1007/s00432-022-04073-3 |
| Transcriptional heterogeneity and clonal evolution incl. CTCs | Goodspeed 2025, doi 10.1158/1078-0432.CCR-24-2040 |
| Cross-study pediatric immune landscape (26 sc/snSeq studies, 272 samples, 35 tumour types) | EJC 2025, doi 10.1016/j.ejca.2025.115708 |

**Consequence: "the immune landscape of Ewing sarcoma" and "an immune prognostic signature" are both occupied,
the second by a paper that already combined single-cell with bulk outcome data.** Repeating either with more
cells is an incremental result, and it is the default thing an atlas paper does.

### What is genuinely open

1. **O-PDX fidelity.** 22.6% of v4 (14,645 Aynaud cells) is orthotopic xenograft, rising when Waltner's 6 O-PDX
   arrive. Nobody has asked how faithfully O-PDX reproduces patient-tumour cell states in EwS, and the
   Next-Generation Pediatric Cancer Models panel (GSE318841 / GSE319124 / GSE330883) provides patient tumours
   matched to their own PDX, organoid and 2D cultures as the control. NBAtlas never faced this because its
   datasets were patient tumours. **This is the clearest unoccupied, EwS-specific, atlas-dependent question.**
2. **Reproducibility of the malignant-state taxonomy across cohorts.** Aynaud's programs were derived in one
   study; no one has tested whether they recur in the other five. An atlas is the only instrument that can, and
   a negative result would be as informative as a positive one.
3. **Externally calibrated fusion activity.** v4 already scores every malignant cell against ESCLA knockdown
   signatures rather than atlas-internal programs — a different and more defensible construct than
   `Merge.GEP.EWS`-style internal scores.
4. **Single-cell vs single-nucleus composition bias in EwS**, with GSE200529 as the calibration control. NBAtlas
   did exactly this for neuroblastoma; it is unpublished for EwS and becomes available the moment the ScPCA arm
   lands.
5. **Anatomical-site biology** (bladder in v4, meningeal in Hong 2022) — currently anecdotal, never systematic.

---

## 3. Structural differences from NBAtlas that matter for reviewers

Verified from the NBAtlas methods:

- **Its reference is entirely public** — 7 datasets from 6 publications, 68 samples, 61 patients. So "an atlas
  built only from public data" is not an objection; NBAtlas is the precedent for it.
- **But it had in-house data for two jobs v4 cannot currently do.** A matched single-cell/single-nucleus pair
  from one patient (58,083 cells vs 63,703 nuclei) corroborated the assay-composition trend, and 3 single-cell +
  5 single-nucleus in-house samples (17,096 cells + 50,080 nuclei) served as the *query* demonstrating
  reference-based mapping. v4's mapping recipe is validated only by held-out projection of a study already inside
  the reference — stability, not generalization.
- **Its discovery hook was outcome-linked tumour states.** v4's current equivalent recapitulates a published
  tumour-stromal prognostic effect in 129 patients, which validates the deconvolution but is not a new finding.
- **It benchmarked integration** (Harmony, Seurat, ComBat) and called malignant cells by CNA (CopyKAT, SCEVAN).
  v4 does neither yet.

---

## 4. Recommendation

The resource is novel and worth building; the discovery claim needs choosing now, before more compute is spent.

1. **Lead with exhaustiveness, not size.** The completeness review proves the atlas covers every public EwS
   single-cell dataset that exists. That reframes 64,702-101,000 cells from "small" to "all there is".
2. **Pick the O-PDX fidelity question as the discovery hook.** It is unoccupied, EwS-specific, answerable with
   data already identified, and it turns v4's largest apparent weakness (22.6% xenograft content) into the
   subject.
3. **Do not lead with the immune landscape or an immune prognostic signature.** Both are published; at most
   confirm them and move on.
4. **Close the three NBAtlas-parity gaps** before submission: integration benchmark, CNA-based malignant calling,
   and a query dataset that is genuinely external to the reference.
5. **Secure an external query cohort.** Hong 2022 or Chen 2025 would serve — and if either becomes available it
   is worth more as held-out validation than as extra reference cells.


---

## 5. A discovery program that is both unoccupied and powered

Measured from `EWSAtlas_v4.h5ad` — these numbers decide what is answerable, and they correct the earlier
recommendation in §4.

**Malignant cells by material type**

| material_type | cells | patients | studies |
|---|---|---|---|
| patient_tumor | 13,455 | 25 | 4 |
| O-PDX | 14,468 | 5 | **1** |
| CTC | 162 | 3 | 1 |
| normal_tissue | 1,630 | 2 | 1 |

**Malignant cells by clinical stratum**

| stratum | cells | patients |
|---|---|---|
| primary | 13,191 | 23 |
| metastasis | 169 | 3 |
| recurrence | 257 | 2 |
| treatment_naive | 9,196 | 17 |
| post_treatment | 2,685 | 10 |

### Correction to §4: O-PDX fidelity is NOT answerable with v4.0 as it stands

All 14,468 O-PDX malignant cells come from Aynaud alone. "O-PDX vs patient tumour" is therefore perfectly
confounded with "Aynaud vs the other four studies" — any difference found would be indistinguishable from batch.
The question becomes answerable only with the Waltner arm, which contains 7 patient tumours and 6 O-PDX from the
same 5 patients, processed in one pipeline: those are the within-patient, within-study pairs the comparison needs.
Keep the question; it is still the best use of the ScPCA arm once it lands, but it cannot be the lead result now.

Likewise CTC claims are not supportable: only 162 of 6,845 CTC-library cells are malignant — the libraries are
overwhelmingly blood, as expected for circulating-tumour-cell enrichment.

### The lead question, answerable now

**Does the EWSR1-FLI1-low (mesenchymal/migratory) malignant state exist in patient tumours at appreciable
frequency, is its abundance a patient-level property rather than a cell-level one, and does it shift with
treatment exposure?**

Why this is the right question:

- **It is the field's central unresolved problem.** Ewing sarcoma has essentially one driver and few recurrent
  secondary mutations, yet outcomes vary widely. The leading explanation for that variability is fusion-activity
  state heterogeneity — the low-activity fraction being the proposed source of migration, metastasis and drug
  tolerance.
- **It is established in models and untested in patients.** Aynaud 2020 defined the axis in 5 PDX and cell lines;
  the EWS::FLI1-fluctuation work (GSE291060/GSE291136) modelled it in bulk. No one has quantified the
  low-activity fraction per patient across multiple single-cell cohorts.
- **v4 can answer it today.** Every malignant cell already carries an externally calibrated `fusion_activity`
  score built from ESCLA knockdown contrasts — not from atlas-internal programs, so the axis is not circular. The
  treatment contrast is the best-powered clinical stratum in the atlas: 17 treatment-naive vs 10 post-treatment
  patients, 9,196 vs 2,685 malignant cells, spanning 4 studies.
- **It extends to outcome.** A fusion-low signature derived here can be deconvolved into the 129 bulk tumours with
  usable survival, giving an independent test that does not depend on the atlas's own patients.
- **It fails informatively.** If the low-activity fraction is uniform across patients, the patient-level version
  of the model is wrong — which is itself worth publishing, and is the kind of negative result only an atlas can
  deliver.

Design notes that keep it honest: analyse per patient, not per cell (cells within a patient are not independent);
stratify by study, since residual study structure is substantial (78.1% same-study neighbours); exclude O-PDX from
the patient-level analysis rather than pooling it; and report the metastasis and recurrence strata as
underpowered descriptive context (3 and 2 patients) rather than as tests.

### Second question, for when the ScPCA arm lands

**Which malignant states does orthotopic xenografting preserve, and which does it lose?** Answerable with
Waltner's matched tumour/O-PDX pairs plus the Next-Generation Pediatric Cancer Models bulk panel
(GSE318841 / GSE319124 / GSE330883) as independent confirmation. If O-PDX systematically loses a state, every
PDX-based EwS preclinical study is blind to it — a result with immediate translational consequences and no
current answer in the literature.

### Third question, cheap and structural

**Do Aynaud's malignant programs reproduce in the other cohorts?** Score the published programs across all five
studies and ask whether they recur as coherent axes. This tests whether the field's working model of EwS
intratumoral heterogeneity generalizes beyond the study that defined it — the kind of cross-cohort check an
integrated atlas exists to make, and which no single-cohort paper can perform.

### What is not feasible with public data

Anything requiring germline or somatic genotype (GGAA microsatellite length, STAG2/TP53 co-mutation effects on
cell state), spatial context, or immune-receptor repertoire. None exists for EwS at single-cell resolution in any
archive searched.
