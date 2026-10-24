# 10 · Publication plan

Written after reading Bonine et al., *Cell Reports* 43:114804 (2024) — NBAtlas —
in full, including its STAR Methods, and inspecting its rendered figures. Every
cohort number below was recomputed from `EWSAtlas_v5.h5ad` for this document,
not carried from earlier notes.

---

## 1. What NBAtlas actually does, figure by figure

| | subject | panels |
|---|---|---|
| F1 | the resource | dataset overview; UMAP by dataset / assay / sample; annotated UMAP with metadata rings; marker expression |
| F2 | *neuroendocrine cells are the malignant cells* | genome-wide CNA score; arm-level CNA (1p loss, 17q gain); adrenergic vs mesenchymal signature scores; cell-cycle phase; SingleR mapping to normal sympathoadrenal development |
| F3 | malignant heterogeneity → risk groups | re-integrated malignant UMAP (7 clusters); marker dotplot; embedding **density** split high-risk vs low-risk; Jaccard heatmap of per-sample NMF programs clustered into 11 metaprograms; functional enrichment of metaprograms |
| F4 | the immune TME | immune UMAP; dotplot; monocyte/macrophage zoom; LAM signature; activation markers |
| F5 | a technical result | single-cell vs single-nucleus: neuroendocrine fraction, matched-sample UMAPs, cell-type proportions |
| F6 | utility | reference-based mapping of new data (scArches/scANVI, SingleR) vs manual annotation |

### The structural lesson

NBAtlas is a *Resource* paper that earns its discovery claims in **one figure of
definition before any figure of interpretation**. F2 exists solely to establish
that the big cluster is malignant — by three independent arguments (CNA,
lineage signature, developmental mapping) — and only then does F3 interpret
heterogeneity within it. Our Figure 3 is the mirror image: it reports that our
equivalent definitional evidence *failed*. That ordering is not optional, and it
is why the discovery work is still blocked.

Two further moves are worth copying directly:

- **The axis comes from outside the atlas.** Adrenergic/mesenchymal scores come
  from van Groningen et al., developmental labels from Jansky et al. Nothing
  that defines a cell state is derived from the atlas it is then used to
  describe. Our ESCLA-derived fusion signature is the same device.
- **Recurrence is proven per sample, not per cluster.** NMF is run on each
  sample separately, specifically to avoid batch correction, then programs are
  matched across samples; they explicitly checked no metaprogram came from a
  single study. For a six-study atlas with 38 patients this is the strongest
  available argument that a state is real.

### Their copy-number recipe, against ours

| | NBAtlas | EWSAtlas v5 (current) |
|---|---|---|
| window | 250 genes, step 1 | 100 genes |
| run unit | **re-run per study, merged** | one global run |
| reference | all types except neuroendocrine, fibroblast, Schwann | immune only |
| scored | genome-wide **and per chromosome arm** | genome-wide only |
| classification | **per sample**, ≥100 reference cells, else `not.run` | none |
| cross-check | CopyKAT **and** SCEVAN | none |

Four of those six differences plausibly contribute to our per-sample AUC of
0.540, and all four are fixable. Arm-level scoring matters most: Ewing sarcoma's
recurrent changes (1q gain, chr 8 gain, chr 12 gain, 16q loss) are few and
focal, so a genome-wide mean of |CNV| averages the signal away — exactly what
Figure 3c shows.

### Figure-style conventions to adopt

- Bold **uppercase** panel letters (Cell Press house style); ours are lowercase.
- A **locator inset** on every zoom panel — the full atlas in grey with the
  subset highlighted — so a reader always knows what was re-clustered.
- **Embedding density** split by clinical group, instead of overplotted
  coloured scatter, whenever groups are compared on a UMAP.
- The caption headline is a **claim**, not a label ("Neuroendocrine cells
  present with molecular hallmarks of malignant neuroblastoma cells"). This is
  compatible with our no-title rule and is how our captions should read.
- A `See also Figure SX` line closing each main caption.

---

## 2. Verified assets for the discovery

Recomputed from the released object:

| asset | value | consequence |
|---|---|---|
| atlas | 138,353 cells · 16,622 genes · 6 studies · 65 samples · 38 patients | |
| malignant by markers | 93,805 (67.8%) | NBAtlas's neuroendocrine share is 67.4% — a close structural parallel |
| fusion signature | 560 EF-high / 389 EF-low genes present, from 107 ESCLA arrays, 18 lines, dox-inducible EWSR1-ETS knockdown | the axis is **perturbation-defined**, not marker-defined — this is what removes the circularity |
| matched tumour + O-PDX | 4 patients (SJEWS001321, SJEWS049193, SJEWS063829, SJEWS063834), 62,719 malignant cells | tests whether the low state needs a human microenvironment |
| treatment arms, patient level | 15 treatment-naive vs 13 post-treatment patients (tumour material) | the clinical axis |
| within-study both arms | Visser 6 post / 5 naive · Cillo 2 / 2 (Mao has 1 naive and 1 not-reported patient, **no** post-treatment arm) | the only two places study is held constant |
| CTC | 6,941 cells, 8 patients, but **405 malignant** | too thin to carry a main claim |
| bulk | 408 harmonized tumours; survival for 129 patients (GSE63157 n=85, GSE17674 n=44) | external validation with purity adjustment |
| falsification controls | dissociation-bias signatures incl. EwS-specific nucleus-vs-cell and warm-dissociation sets | the fusion-low state must survive these |

### Confounds that constrain the claims

1. **Treatment is confounded with material type.** Of malignant cells labelled
   post-treatment, 41,102 are O-PDX and 29,810 are patient tumour; *no* O-PDX
   cell is treatment-naive. Any treatment contrast must therefore be run on
   patient-tumour material only, patient-level, with study as a covariate.
2. **The matched pairs are single-study.** All four tumour/PDX patients come
   from Waltner 2025, so that result cannot be replicated across studies.
3. **Prior art is partly occupied.** OpenScPCA already publishes EWS-high and
   EWS-low labels for the Waltner nuclei. The novel contribution is the
   cross-study, patient-level, treatment- and xenograft-linked framing with an
   externally calibrated axis — not the existence of a fusion-low state.

---

## 3. Proposed main figures

Figures 1 and the four supplementaries exist today. Everything else is work.

| | claim the figure must earn | status |
|---|---|---|
| **F1** | a harmonized six-study Ewing sarcoma reference exists and is annotated on marker evidence | **done** |
| **F2** | the malignant compartment is defined by evidence independent of the markers that named it — arm-level CNA, fusion-activity score, cell cycle | **blocked** on the CNA fix |
| **F3** | malignant heterogeneity is recurrent across patients and studies — re-integrated malignant UMAP, per-sample NMF metaprograms, cross-study occupancy | not started |
| **F4** | **the discovery** — a fusion-activity axis that is continuous, patient-recurrent, separable from cell cycle and depth, and associated with treatment and with xenografting | not started |
| **F5** | the Ewing immune microenvironment across studies and treatment arms | not started |
| **F6** | the fusion-low signature stratifies bulk cohorts after purity adjustment | substrate ready |

Supplementaries: S1–S4 as published today, plus S5 copy-number diagnostics and
caller agreement, S6 NMF robustness and study occupancy, S7 the falsification
panel (dissociation, depth, doublets, ambient), S8 immune detail.

---

## 4. Execution order

The order is forced by the dependency, not by convenience:

1. **Fix the copy-number caller** — window 250, re-run per study, arm-level
   scores for 1q / 8 / 12 / 16q, per-sample classification restricted to
   samples with ≥100 reference cells, and a second caller for agreement. If it
   still fails, Figure 3 stands as published and F2 rests on the fusion score
   and cell cycle alone, declared as such.
2. **Freeze the malignant set** on the union of CNA and marker evidence, and
   report the disagreement rather than resolving it silently.
3. **Score the fusion axis** with the ESCLA signature; regress out cell cycle
   and detection depth; show the axis survives.
4. **Per-sample NMF** on malignant cells, metaprograms by Jaccard clustering,
   with the explicit single-study check NBAtlas applied.
5. **Test the associations** — treatment (patient-level, tumour material only,
   study as covariate, confirmed within Visser and Cillo) and xenografting (the
   four matched patients).
6. **Falsify** — dissociation, ambient RNA, doublet enrichment, depth.
7. **Project into bulk** with purity adjustment; survival on 129 patients.
8. **Re-render all figures** in Cell Press style and write the paper.

Steps 1–2 gate everything after them. If step 1 fails and step 3's axis cannot
be separated from depth, the honest outcome is a Resource paper with a negative
methods result, not a discovery paper — and that decision is made at step 3,
on the measurement, not in advance.
