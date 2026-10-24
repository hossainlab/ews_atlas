# 05 · Annotation and malignant-cell calling

**Status: v5 annotated. 34 leiden communities at resolution 1.0 on the primary embedding; 11 cell types plus `Unassigned`.**

## Vocabulary

Fifteen labels: `Ewing_sarcoma`, `Ewing_proliferating`, `T_cell`, `NK_cell`,
`B_cell`, `Plasma`, `Monocyte_Macrophage`, `DC`, `pDC`, `Endothelial`,
`Fibroblast_CAF`, `MSC_Osteoblast`, `Urothelial`, `Platelet`, `Cycling_other`.

Marker sets are listed in `config.annotation` and in
`results/tables/v4_celltype_by_study.csv`.

## Three schemes were tried; two failed

1. **Fraction-of-cells thresholds** — failed outright: B-cell clusters were
   called Ewing, and urothelium appeared in a PDX-only study where no urothelium
   can exist.
2. **Rank-percentile argmax** — worse: 4,403 spurious platelets.
3. **Standardized z-score argmax with two guards** — accepted.

The two guards encode knowledge the markers alone cannot carry:

- Urothelial is callable only where the bladder study contributes ≥50% of a
  cluster. Urothelium exists in exactly one deposit; a urothelial call anywhere
  else is a marker artefact, not a discovery.
- Human immune labels are refused for clusters that are >80% orthotopic PDX.
  Those are mouse-host-derived compartments in a human-reference alignment.

Sanity checks after annotation: Aynaud O-PDX 98.8% malignant, Cillo 69.8%
immune, 2,537 of 2,543 urothelial cells in the bladder study.

## v5 result

| cell type | cells |
|---|---|
| Ewing_sarcoma | 83,814 |
| T_cell | 16,025 |
| Monocyte_Macrophage | 10,081 |
| Ewing_proliferating | 9,991 |
| Fibroblast_CAF | 6,063 |
| Endothelial | 4,109 |
| Unassigned | 2,809 |
| Urothelial | 2,453 |
| B_cell | 1,679 |
| DC | 715 |
| Plasma | 396 |
| Platelet | 218 |

Malignant (Ewing_sarcoma + Ewing_proliferating) = 93,805 cells, 67.8%, spanning
all 38 patients and all six studies.

Composition sanity checks: Aynaud orthotopic PDX 90.8% malignant, Cillo 82.5%
immune, 2,433 of 2,453 urothelial cells inside the bladder study.

## External validation

Agreement with the depositors' own annotations on the 69,400 cells they label
(OpenScPCA for Waltner, the Goodspeed cell metadata): **93.3% on malignant vs
non-malignant**. These labels were never used to produce ours, so the agreement
is independent.

The disagreements are informative rather than noise:

- **844 cells we call `Fibroblast_CAF` are `Malignant` to the depositors**, and
  551 we call `Ewing_sarcoma` are `Stromal` to them. This is precisely the
  boundary the discovery question lives on — a fusion-low malignant cell
  expresses a mesenchymal programme, and marker annotation has no way to tell
  it from a CAF. Copy-number evidence is what separates them, which is why the
  malignant set is not frozen on markers alone.
- **2,362 cells we call `T_cell` are `Malignant` to the depositors.** These sit
  in clusters with clear CD3D/CD3E/TRAC expression, so our call is the better
  supported one on the markers; the OpenScPCA pipeline appears to assign whole
  neighbourhoods of tumour nuclei.

## Three marker-logic failures found and fixed in v5

Per-cell marker argmax was tried once more on v5 and failed the same way as in
v4 — 12,365 spurious platelets, 22 Ewing cells in a bladder-Ewing study,
urothelium inside a PDX-only study. Marker expression is too sparse per cell to
argmax. Labelling is therefore at **cluster** level throughout.

Even at cluster level, a z-score alone is not enough: it says a panel is
relatively higher here than elsewhere, not that the panel is expressed at all.
Three calls in the first annotated build were wrong for that reason, and each
produced a guard:

1. A 7,739-cell cluster took **NK** on NKG7/GNLY while expressing CD3D 1.29 —
   cytotoxic T cells share the cytotoxic panel. NK is now refused when CD3 is
   comparable.
2. A 5,763-cell cluster took **DC** with LYZ at 4.05 and CD1C at 0.20 — a
   monocyte. Each label now requires at least one of its own obligatory markers
   above a floor (`REQUIRES`).
3. Four clusters with essentially no lineage markers inherited whichever panel
   ranked first. A label is now refused below a minimum cluster z
   (`MIN_ACCEPT_Z = 0.5`), which is where the 2,809 `Unassigned` cells come
   from — an honest 2.0% rather than a confident wrong answer.

Nucleated immune labels additionally require PTPRC. Platelets are exempt
because they are anucleate and genuinely CD45-negative — the first version of
that rule rejected a clean PPBP/PF4 cluster and let it fall through to
Endothelial.

Guards that fired are recorded per cluster in
`results/tables/v5_cluster_annotation.csv` with their reasons.

## MSC_Osteoblast is not asserted in v5

The one cluster that scored MSC/osteoblast highest expressed RUNX2/SP7/IBSP at
0.02–0.05 and became `Unassigned`; a second became `Monocyte_Macrophage` once
SPP1 was excluded as obligatory evidence. v5 therefore makes **no claim of an
MSC/osteoblast compartment**. Given that SPP1+ macrophages and osteoblasts
overlap in bone, this is the correct outcome of the available markers rather
than a gap to paper over.

## Composition-adjusted mixing

With cell types in hand, the integration figure from
[04_integration.md](04_integration.md) can be adjusted. Same-study kNN under the
primary embedding is 0.749; the expectation is 0.338 if studies were mixed at
random, but **0.498 once each cell type's own study composition is taken into
account**. So roughly half the apparent batch structure is composition, not
batch effect.

The residual is not uniform, and that is the useful part
(`results/tables/v5_mixing_by_celltype.csv`): Urothelial (0.943 vs 0.984
expected), DC (0.958 vs 0.911) and Platelet (0.980 vs 0.973) are fully
explained by composition, while T_cell (0.635 vs 0.225) and
Monocyte_Macrophage (0.736 vs 0.347) are not. Study-associated structure
persists in the immune compartments and should be treated as a limitation for
any immune claim.

## scANVI label transfer — not run

The plan called for scANVI to propagate labels. It was not run: the
marker-and-guard annotation already agrees with the depositors at 93.3% on the
compartment that matters, and scANVI's main contribution would be to spread
depositor labels into the four studies that have none — which would import
OpenScPCA's `T_cell`→`Malignant` behaviour described above. This is a
substitution, not a completed step, and is recorded as such.

## Known unresolved label

`MSC_Osteoblast` versus SPP1+ macrophage is not cleanly separated by marker
argmax — both express SPP1 and IBSP-adjacent programmes in bone. scANVI label
transfer plus explicit PTPRC/AIF1 arbitration is the planned resolution; until
then any statement about the mesenchymal compartment carries this caveat.

## Malignant calling

v4 called malignant cells by marker expression alone (45.9% of cells). v5 adds
**infercnvpy** CNA inference, using the immune compartment as the diploid
reference, and freezes a malignant set on the intersection of marker and CNA
evidence. The reason is specific to this disease: the fusion-low mesenchymal
state the discovery question targets is defined by *low* expression of the
fusion programme, so a marker-only definition of "malignant" would preferentially
exclude exactly the cells under study. CNA evidence is independent of the
transcriptional programme and breaks that circularity.

NBAtlas used CopyKAT + SCEVAN (both R); infercnvpy is the Python equivalent and
avoids the Bioconductor load failures seen in this environment.

## Copy-number inference — attempted, and it does not work as implemented

**Status: negative result. The malignant set is NOT frozen on copy-number
evidence, because the evidence does not discriminate.**

infercnvpy was run per sample with the immune compartment (T, B, plasma,
monocyte/macrophage, DC, platelet) as the diploid reference, 100-gene windows,
step 10, on 16,271 genes with Ensembl coordinates (97.9% of the feature space).
47 of 65 samples had at least 50 reference cells and were scored, covering
113,049 cells; 18 samples have no internal diploid reference at all.

The per-cell score is the mean absolute inferred CNV across the genome. It
separates poorly:

| | median score |
|---|---|
| Unassigned | 0.01479 |
| Fibroblast_CAF | 0.01434 |
| Endothelial | 0.01411 |
| Ewing_proliferating | 0.01340 |
| **Ewing_sarcoma** | **0.01169** |
| Monocyte_Macrophage | 0.01008 |
| T_cell | 0.00869 |
| DC | 0.00570 |

Malignant vs reference AUC is **0.677** overall and, per sample,
**median 0.540** (range 0.330–0.962; above 0.7 in only 9 of 37 samples with at
least 30 cells of each).

The ordering is the decisive part: `Fibroblast_CAF`, `Endothelial` and
`Unassigned` all score *higher* than `Ewing_sarcoma`. A working copy-number
caller cannot put stromal and endothelial cells above the tumour. The statistic
is therefore tracking something else — most plausibly per-cell expression
variance and library quality, which is exactly what a mean-absolute deviation
score picks up when the underlying CNV signal is weak.

Three contributing reasons, in the order they should be addressed:

1. **The statistic is the weak one.** The field standard (Tirosh et al.) uses
   two quantities together — a CNV *signal* (mean of squared values) and the
   *correlation* of each cell's profile with a consensus malignant profile from
   the same sample — and separates cells in that plane. Mean absolute deviation
   alone is known to be poor.
2. **Many samples have thin references.** Six of the scored samples carry fewer
   than 100 reference cells, and 18 samples have none.
3. **Over half the atlas is single-nucleus.** CNV inference from nuclear RNA is
   noisier than from whole cells, and nothing here accounts for that.

Consequence for the discovery question, stated plainly: the 844 cells called
`Fibroblast_CAF` here and `Malignant` by the depositors **remain unresolved**.
Any analysis of a fusion-low mesenchymal state has to carry that ambiguity
explicitly until a working copy-number caller exists, because a marker-based
malignant definition is circular with respect to exactly that state.

Per-sample results: `results/tables/v5_cnv_per_sample.csv`,
`v5_cnv_auc_per_sample.csv`, `v5_cnv_assessment.json`.

## Reference mapping

`src/map_query_to_ewsatlas.py` projects a query dataset onto the atlas with
scArches surgery. Measured on 2,000 held-out Visser cells carrying unseen batch
levels: **1,680/2,000 (84.0%)** recovered their label after 5 surgery epochs;
disagreement was almost entirely `Ewing_sarcoma` ↔ `Ewing_proliferating`
(T cells 378/381, monocytes 65/66, B cells 21/21). Urothelial transfer was
**not** measured — Visser contains no urothelium.

This demonstrates projection stability on held-out cells from a *contributing*
study. It is not evidence of generalization to an unseen cohort, because no
external EwS single-cell cohort is obtainable (see
[01](01_data_sources.md)).
