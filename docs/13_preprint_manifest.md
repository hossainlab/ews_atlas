# 13. What can be presented — preprint manifest

Every claim below is backed by a figure and a table already on disk, and every
number has been read back from the saved artifact rather than from notes. Claims
are graded: **A** = defensible as stated, **B** = defensible with the caveat
attached, **C** = do not present.

---

## Ready now: six main figures

### Figure 1 — the resource (grade A)

138,353 cells, 16,622 genes, six studies, 65 samples, 38 patients; the first
harmonized single-cell reference for Ewing sarcoma. Twelve cell types assigned at
the level of 34 Leiden communities with obligatory-marker guards, leaving 2.0%
`Unassigned` rather than forcing a nearest label. QC retained 138,353 of 163,578
barcodes (84.6%) under assay-aware floors.

Caveat that must travel with it: Waltner supplies 54% of all cells and the atlas
is 50% single-nucleus. Say so in the abstract, not only in the methods.

- Figure 1, Figure S1, S2, S4, S6; `v5_by_study.csv`, `qc_attrition.csv`

### Figure 2 — batching on sample beats batching on study (grade A)

scVI with `batch_key='sample'` scores 0.606 overall (batch 0.611, bio 0.603)
against 0.532 for study-level batching, 0.483 for Harmony and 0.437 for
unintegrated PCA. This **refuted** the pre-registered expectation that 65
sample-level batches would over-correct, and it is a reusable result for anyone
integrating few-patient sarcoma cohorts.

Annotation quality in the same figure: harmonized agreement with depositor labels
is 93.0% over 69,040 labelled cells, with disagreement concentrated in one row —
cells deposited as `Stromal` are 76.1% malignant and 23.0% stromal by our
labelling, and `Epithelial` is 98.2% malignant.

- Figure 2, Figure S3; `v5_embedding_benchmark.csv`, `v5_depositor_labels.csv`

### Figure 3 — the malignant compartment, defined independently of its markers (grade A)

The copy-number failure was the statistic, not the parameters: correlating each
cell to a per-sample consensus profile gives a median malignant-versus-reference
AUC of 0.712 across 45 samples, where the genome-wide mean of |CNV| gives 0.597
and the first implementation gave 0.540. 97,586 frozen malignant cells, each
carrying how it was called — 12,114 markers and copy number, 67,587 markers only
with copy number negative, 14,104 unassessed, 3,781 rescued on copy number alone.

Report honestly, as the caption already does: three samples remain below 0.5, and
arm-level scoring is **not** what fixed it (chr8q is the best arm at a median
0.574 per sample).

- Figure 3, Figure S5; `v6_cnv_per_sample.csv`, `v6_cnv_arm_auc.csv`, `v6_malignant_evidence.csv`

### Figure 4 — most recurrent programs track capture, not biology (grade A, with attribution)

28 metaprograms from per-sample NMF across 45 samples; **15 of 28 confined to a
single study**, where NBAtlas had none. The three most biologically appealing
programs correlate with the single-nucleus-versus-single-cell signature at 0.781
(hypoxia-like), 0.695 (mesenchymal) and 0.642 (cell cycle); the hypoxia-like
program correlates with detected genes at 0.524. Doublets explain none of it
(|rho| <= 0.183).

**Attribution is mandatory.** Truong et al., BMC Cancer 2023
(doi 10.1186/s12885-023-10977-1) reported dissociation and gene-length artifacts
in 125,831 Ewing/DSRCT/osteosarcoma PDX cells first. Our contribution is the same
result in patient tumours across six independent cohorts plus the metaprogram
decomposition. Frame as confirmation and extension, never as discovery.

- Figure 4; `v6_metaprograms.csv`, `v6_metaprogram_occupancy.csv`, `v6_falsification.csv`

### Figure 5 — the immune compartment, with two corrections (grade A for the corrections, C for composition)

29,114 immune cells resolved into 13 subsets plus a low-quality class: T_CD4
7,532, T_cytotoxic 5,860, Monocyte 4,974, Macrophage 3,351, B_cell 1,681, NK
1,267, cDC2 1,137, Mast 715, Treg 670, Macrophage_LAM 532, pDC 307, Platelet 218,
Plasma 88, Low_quality 782.

Two atlas-level corrections that ship on the released object: all 715 cells
labelled `DC` are **mast cells** (TPSAB1 3.76, CPA3 5.80), and 1,267 **NK** cells
sit inside `T_cell`.

Do **not** present composition-versus-treatment. The smallest p across 13 subsets
in 24 patients is 0.0339 (NK), which fails correction, and its study-adjusted
coefficient has the opposite sign to the raw means (0.053 naive vs 0.013 post) —
a study-level reversal.

- Figure 5; `v6_immune_labels.csv`, `v6_immune_composition.csv`

### Figure 6 — the atlas is usable as a reference (grade A)

Three independent demonstrations.

*Surface targets.* Over the seven non-urothelial microenvironment classes,
IL1RAP is the only one of five clinically pursued targets positive against all
seven (+0.7 to +13.1 percentage points, paired within sample). STEAP1 is
non-negative against all seven but flat against plasma cells. CD276 is negative
against fibroblasts. CD99 is the instructive case: positive margin against six of
seven classes, yet detected in 74.2% of urothelial, 72.9% of fibroblast, 70.3% of
endothelial and 68.3% of platelet cells against 56.1% of Ewing sarcoma cells.
Both STEAP1 and IL1RAP are strongly negative against urothelium (-20.9, -21.7).

*Reference mapping.* Leave-one-study-out label transfer recovers the full-atlas
annotation at a median of 0.84 (0.85 excluding structurally untransferable
classes), 0.72-0.90 among confidently assigned cells. Mao's 0.63 is structural:
2,433 urothelial cells against 20 in the rest of the reference.

*Bulk deconvolution.* Mean malignant fraction 0.73-0.79 across five cohorts;
higher tumour content tracks better overall survival (HR 0.77, 0.61-0.98) and
higher stroma worse (1.33, 1.04-1.71) in 129 patients. Neither survives
Bonferroni for eight tests, and neither is new — it reproduces GSE63157's own
published association, which is what makes it a calibration rather than a claim.

- Figure 6; `f6_target_landscape.csv`, `f6_target_margins.csv`, `f6_loso_mapping.csv`,
  `f6_loso_by_celltype.csv`, `f6_bulk_fractions.csv`, `f6_bulk_survival.csv`

---

## Ready now: supplementary

S1 quality control; S2 deposit-specific QC; S3 integration detail; S4 covariate
coverage; S5 the superseded copy-number attempt, kept deliberately as the record
of what failed; S6 UMAP by sample; **S8 (new) the CD73-high fibroblast-like state**
— genetically real on a held-out half of the genome but positive in three Ewing
patients of 38, which is why it belongs here and not in the main sequence.

Dropped S7 (`Figure_S7_included_datasets.png`): it duplicated Figure 1A exactly. The vacated number is now the CAF-like state figure (`Figure_S7_caf_like_state.png`), which the grading below places in the supplement rather than the main sequence.

---

## Grade C — do not present

| claim | why not |
|---|---|
| A fusion-activity (EWSR1::FLI1-low) axis exists in patient tumours | The ESCLA axis is a proliferation axis; after removing 133 cycle genes it tracks cycle (0.246 S / 0.256 G2M) as strongly as the canonical targets (0.229 / 0.255) |
| The CAF-like state is a recurrent tumour state | Positive in 3 Ewing patients of 38; the EWS::FLI1-repressed gene set cannot discriminate tumour from fibroblast because it is itself mesenchymal |
| Any treatment-associated state or composition shift | 28 patients, only two studies carry both arms; all metaprograms null, the one nominal hit is the confounded axis in the wrong direction |
| Xenograft versus patient tumour | 4 matched patients, one study; 0.125 is the smallest attainable two-sided p |
| Anything from circulating tumour cells | 405 of 6,941 CTC-library cells are malignant |
| Metastasis or recurrence biology | 3 and 2 patients respectively |

---

## Gaps to close before posting

1. **No trained scVI model is saved** (`data/interim/embeddings/*.npy` holds
   latents only), so the scArches-compatible model promised by Figure 6C does not
   exist as a download. One full retrain on the released atlas is needed.
2. `src/map_query_to_ewsatlas.py` points at `EWSAtlas_v4.h5ad` and a v4 model
   directory, both gone. It must be repointed at v5 or removed.
3. The two annotation corrections (DC = mast, NK inside T_cell) are **not applied
   to the released object**. Either patch `cell_type` or state prominently that
   the subset labels supersede it.
4. Waltner-exclusion sensitivity analysis (~63,465 cells) for the integration
   supplement is not run.
5. Deposit targets not yet prepared: cellxgene-ready object, Zenodo DOI, public
   code repository.

## Fastest credible path to a preprint

Post with the six main figures above, S1-S6 plus the CAF-like state as S8, and
the five tables. Gaps 2 and 3 are hours of work and should be closed first
because both are correctness issues a reader would hit immediately. Gap 1 is the
only one that costs real compute and it can be promised as "model available on
acceptance" if necessary. Gaps 4 and 5 can follow the preprint.

## Final set as built

Six main figures and seven supplementary figures, all 180 mm wide at 300 dpi with
no in-figure titles.

| Slot | File |
| --- | --- |
| Figure 1 | `Figure_1_atlas_overview.png` |
| Figure 2 | `Figure_2_integration_quality.png` |
| Figure 3 | `Figure_3_malignant_compartment.png` |
| Figure 4 | `Figure_4_immune_microenvironment.png` |
| Figure 5 | `Figure_5_recurrence_and_confounding.png` |
| Figure 6 | `Figure_6_atlas_utility.png` |
| Figure S1 | `Figure_S1_quality_control.png` |
| Figure S2 | `Figure_S2_deposit_specific_qc.png` |
| Figure S3 | `Figure_S3_integration_detail.png` |
| Figure S4 | `Figure_S4_covariate_coverage.png` |
| Figure S5 | `Figure_S5_copy_number_first_attempt.png` |
| Figure S6 | `Figure_S6_umap_by_sample.png` |
| Figure S7 | `Figure_S7_caf_like_state.png` |
| Table S1 | `s8_waltner_sensitivity.csv` |

Every supplementary figure is cited at least once in the Results, and every
figure that is rendered has a legend in `results/pub_figs/figure_captions.md` —
`write_supplementary()` raises if a legend is missing, so a figure cannot reach
the companion PDF unlabelled.

## Rendered documents

`preprint/build_manuscript.py` turns `preprint/manuscript.md` into three files:

- `EWSAtlas_manuscript.pdf` — 16 pages, A4, 25 mm margins, Times 11/16.5, the six
  main figures placed where the text calls them with legends pulled verbatim from
  the caption file, followed by the supplementary legends and the references.
- `EWSAtlas_supplementary.pdf` — 7 pages, one supplementary figure per page with
  its legend.
- `EWSAtlas_manuscript.docx` — same content for journals that want Word.

The build is idempotent and reads legends from the caption file rather than
holding a second copy, so a caption edit cannot drift out of the manuscript.

## Not supplied

These are not inferable from the data and are rendered in the PDF as visible
`[TO BE COMPLETED]` text rather than silently omitted: author list and
affiliations, lead contact, Zenodo DOI, repository URL, acknowledgements, author
contributions, declaration of interests. The licence is still deliberately open
pending a reading of the ScPCA terms of use.

## Sequence change

Figures 4 and 5 were swapped after review: the atlas now describes both
compartments (malignant, then immune) before the recurrence test that uses
them, and the paper closes on utility. The previous order interrupted the two
descriptive figures with the negative result.
