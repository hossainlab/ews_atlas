# Decision log

Each entry: the decision, the alternatives, and why. Ordered by when it was
taken, not by importance.

## ADR-001 · Scope is EWSR1::FLI1 Ewing sarcoma only

Other FET-rearranged tumours (DSRCT, clear cell sarcoma, FET-CREB, myxoid
liposarcoma) share the fusion grammar but not the lineage. Pooling them would
make every statement about "the malignant compartment" ambiguous. Cell lines are
excluded from the atlas object but used as external references.

## ADR-002 · v3 deleted and rebuilt from the original deposits

v3 (129,276 cells) merged two incompatible metadata schemas — `Fusion_Status`
populated only for Waltner (45.7% of cells), `Fusion_Type` only for Visser
(5.1%) — mislabelled xenograft status for 14,702 Aynaud cells (true xenograft
≥ 52,985, 41%), and carried `is_doublet = False` for every cell while
scDblFinder had flagged 7,571. 17,453 annotations were review-flagged. Patching
metadata on top of that object would have left the ingest bugs underneath it.

Consequence to state plainly: **v3's cell count is not comparable to v4/v5.**
Most of the difference is the Waltner nuclear arm plus the assay-aware floor,
not lost data.

## ADR-003 · Assay-aware gene floor, uniform criterion

See [02_quality_control.md](02_quality_control.md) for the full argument. Short
form: one integer across a whole-cell and a nuclear assay is uniform only in
arbitrary units, is harsher on nuclei, and removes the smallest nuclei
preferentially — converting a technical threshold into composition bias. The MAD
criterion is identical everywhere; the backstop floor is per assay.

**Amended 2026-10-18, after measuring it.** The decision stands but its stated
magnitude was wrong: the floor change recovers 3,414 nuclei (~4% of the ScPCA
arm), not ~30,000 (45.2%). Pooled medians on the filtered matrices are 771 genes
/ 1,071 UMIs, not 231 / 296. The justification is directional, not quantitative.

## ADR-003b · Library-level quality gate

Measuring the arm surfaced the decision that actually matters. SCPCL001111
(median 44 genes, 76 UMIs, 98.5% under 200 genes; ScPCA retained 1.4% of it) is
a failed library, and no per-cell rule catches it — the per-sample MAD is
computed within the library, so it looks internally consistent. Any library
whose median detected-gene count falls below its assay floor is excluded whole:
one library here, 6,139 of 92,083 nuclei.

## ADR-003c · Ingest ScPCA's `filtered`, not `processed`, matrices

ScPCA's processed matrices already enforce a 200-gene minimum (measured per
library: 200–402). Starting there would silently adopt the floor this project
argues against; the emptyDrops-filtered matrices retain nuclei down to 3 genes
and leave the choice to our own policy.

## ADR-004 · Feature space = intersection of annotations ∩ any-detection

Strict detected-intersection discards ~3,100 genes for the shallowest arm's
sake; the union admits structural zeros an integration model reads as biology.

## ADR-005 · scVI primary, Harmony baseline, scGPT excluded

Harmony corrects a PCA rather than the count model, so agreement between the two
is evidence against a method artefact. Zero-shot foundation embeddings
underperform scVI on batch integration in published benchmarks and add a large
dependency surface for no measured gain.

## ADR-006 · Two batch-key arms, decided by benchmark

v4 used `sample` (56 levels) and measured 78.1% same-study nearest neighbours
against 21.1% expected. v5 trains `sample` and `study` arms and declares the
primary embedding from `scib-metrics`, not by construction.

## ADR-007 · Doublets: explicit scrublet scores, expected-rate calls

`scikit-image` — required by scanpy's automatic threshold — is refused by host
policy here, so the cut is supplied and recorded rather than inferred. Calls are
made at the chemistry's expected multiplet rate (0.8% per 1,000 cells recovered)
rather than at one score threshold, because multiplet formation scales with
library size: the criterion is a uniform instrument prior, the cut adapts.

This matters more in this atlas than in most. A malignant-plus-stromal doublet
mimics the fusion-low mesenchymal state under study, so under-calling would
manufacture the finding. Measured on the Waltner arm, the fixed 0.25 cut removed
2.9% of nuclei against 6.6% for the expected-rate rule; ScPCA's independent
scDblFinder still flags 8.8% of the nuclei we retain (score rho 0.41). That
disagreement is reported, not tuned away, and the discovery step tests directly
whether fusion-low cells are enriched for either doublet flag.

## ADR-007b · A silent failure is a defect, not a result

The first Waltner ingest marked all 82,496 nuclei `failed` and reported
`drop_doublet = 0`; the cause was a sandbox denial of joblib's named pipe inside
a subprocess, masked by a bare `except`. Stage code now raises on any per-sample
failure and records `doublet_method`, `doublet_not_assessed` and
`doublet_errors` separately.

## ADR-008 · CNA-based malignant calling added in v5

The discovery question is about cells with *low* fusion-programme expression, so
a marker-only malignant definition would preferentially exclude the cells under
study. infercnvpy evidence is independent of the transcriptional programme and
breaks the circularity. Python equivalent of NBAtlas's CopyKAT + SCEVAN;
Bioconductor is unusable in this environment.

## ADR-009 · `not_reported` sentinel, never blank

Blanks conflate "not measured", "not applicable" and "we failed to curate it".
The covariate table has zero blanks by construction and the join is asserted to
cover every sample key.

## ADR-010 · Bulk cohorts justified only as the outcome substrate

The atlas has no survival data. Bulk is not part of the atlas object. The
de-duplication findings (GSE63155 ⊂ GSE63157; GSE17618 ≡ GSE17674) matter to
anyone reusing these cohorts and are documented in
[07_bulk_validation.md](07_bulk_validation.md).

## ADR-011 · Scripts drive the rebuild; stages ship when executed

Stages 03–05 are not pre-written as untested stubs. A script in `src/` has run
at least once. This is why `src/pipeline/README.md` carries a status column.

## ADR-012 · Reference mapping validates projection, not generalization

No external EwS single-cell cohort is obtainable — Hong 2022 and Chen 2025
deposit no accession anywhere. The 84.0% held-out label recovery is projection
stability on cells from a contributing study, and is reported as such.

## ADR-013 · CEL-seq2 values are rounded at ingest

Visser's collision-corrected matrices are 1.8% fractional. The options were to
keep the fractions (undefined under a count likelihood), re-derive raw UMI
counts (not deposited), or round. Rounding moves 1.8% of entries by less than
one molecule each and the median nonzero value is 1.0, so the distortion is
smaller than the correction it undoes. The corrected matrix is retained as
`layers['celseq_corrected']` so the choice is reversible.

## ADR-014 · A known-failing invariant stays failing, not skipped

`test_atlas_has_counts_layer_and_no_empty_cells` fails against v4 because v4 has
no count matrix. Marking it `xfail` would hide the one defect that blocks
retraining. It fails until v5 fixes it.

## Known gaps, not decisions

- Ambient RNA not corrected (deposits do not uniformly ship raw droplet matrices).
- Proportional-hazards assumptions not assessed.
- `MSC_Osteoblast` vs SPP1+ macrophage not cleanly resolved.
- `fusion_breakpoint_type` is 0/68 samples — not deposited by any study.
- The fusion-activity signature was derived with RMA reimplemented in Python,
  not SST-RMA; fold changes are compressed and NKX2-2, CCND1 and FLI1 fall below
  the cutoff.


## v6 · copy number, the fusion axis and the figure reframe

- **The copy-number statistic, not its parameters.** Window size, stride and
  reference width were each tested and rejected before the statistic was
  changed. Recorded because the parameter hypotheses were the plausible ones and
  measuring them cost less than arguing about them.
- **The consensus seed excludes reference cells.** Seeding from the top-signal
  quartile of all cells inverts the score wherever malignant cells are a
  minority (038_CTC: 64 malignant against 3,386 reference, AUC 0.368). Seeding
  on reference identity is not circular with respect to the fusion programme.
- **Arm-level scoring was reported honestly, not favourably.** Pooled per study
  it looked like a clean recovery of the recurrent Ewing pattern; per sample
  across all 45 samples it is weak. The figure shows the per-sample version.
- **The fusion-low discovery was abandoned on evidence.** The axis is
  proliferation-confounded and MP24 is the rescued fibroblasts. The paper leads
  on the capture-confounding result instead, which is measured and survives its
  own checks.
- **Atomic writes everywhere.** An in-place `write_h5ad` destroyed
  `EWSAtlas_v5.h5ad`; every stage now writes a temporary file and replaces the
  target, and the figure composer falls back to an in-place copy when a viewer
  holds the PNG open.

## Figure 4 / Figure S7: the planned swap was not carried out, and why

The approved architecture had the CAF-like malignant state as main Figure 4,
with the capture-confounding figure dropping to supplementary because Truong et
al. 2023 (*BMC Cancer* 23:1105, PMID 37254069) had published the underlying
artifact. Both halves of that swap were reversed once the results existed.

**The CAF-like state did not earn a main figure.** At matched depth on a held-out
half of the genome its cells carry their own patient's copy-number profile
(specificity 0.229 against 0.083 for known-diploid reference cells), so the state
is genetically real. But it clears its own sample's reference in only 6 of 45
samples and 4 patients, one of which is a non-Ewing control — three Ewing
patients of 38. A state present in 8% of patients is a report of a rare
observation, not a headline for a cohort-scale resource, and placing it first
would advertise a recurrence the data do not support.

**The capture figure is not the scooped result.** Truong et al. established the
artifact in matched cell/nucleus pairs from patient-derived xenografts and showed
it is computationally correctable; they are credited for that in the legend, in
the Introduction and in the Results. What this figure shows is different and is
ours: that in a six-study patient atlas the artifact becomes the dominant signal
in the recurrence analysis itself — 15 of 28 metaprograms confined to one study,
the three most convincing candidates correlating with the nucleus-versus-cell
signature at 0.78, 0.69 and 0.64 — so a reader who follows standard practice on
this kind of cohort will recover technical programs and read them as biology.
That is the caveat the resource exists to issue, which is what makes it main
text.

This reverses a decision the user approved, so it is theirs to overturn; the
figures are interchangeable by number alone (`06_figures.py` builds both
regardless) and the swap costs one render.
