# 11. The CD73-high fibroblast-like malignant state

Stage: `src/pipeline/07_caf_state.py`
Outputs: `results/tables/f4_*.csv`, `results/pub_figs/Figure_4_caf_like_state.png`

## Question

The frozen malignant set (docs/05) contains 3,781 cells admitted on copy-number
evidence alone — 2,220 originally annotated `Fibroblast_CAF`, 1,561
`Unassigned`. Figure 4 of the earlier figure set read these as an artefact of
the rescue step, because the candidate fusion-low metaprogram MP24 scores 0.503
in them and near zero in every marker-defined class.

Wrenn et al. (*Clin Cancer Res* 2023;29:5140-5154, PMID 37471463) report the opposite
reading: a CD73 (NT5E)-positive mesenchymal subpopulation of Ewing sarcoma cells
that upregulates EWS::FLI1-repressed genes and deposits extracellular matrix.
Under that reading the same cells are a real tumour state, and dismissing them
closed the falsification prematurely. This stage adjudicates.

## Why transcription cannot settle it

Scoring the four candidate groups on the ESCLA fusion-activity signatures gives
an EWS::FLI1-repressed score of 0.261 in the rescued cells and 0.237 in
copy-number-negative fibroblast-like cells, with matrisome/CAF scores of 0.861
and 0.962. The two are indistinguishable. This is expected rather than
surprising: the EWS::FLI1-repressed gene set *is* a mesenchymal gene set, so any
fibroblast scores high on it. "Expresses EWS::FLI1-repressed genes" therefore
has no discriminating power between a de-differentiated tumour cell and a
fibroblast, and the published criterion cannot be applied transcriptionally in a
mixed atlas.

CD73 itself does discriminate partially — detected in 21.8% of rescued cells
against 9.9% of copy-number-negative fibroblast-like cells and 6.7% of reference
cells — but 6.7% of diploid reference cells are also positive, so the marker
alone cannot assign a cell.

## The genetic test, and why it is not circular

The cells were *selected* by correlating their genome-wide copy-number profile
to a per-sample consensus. Re-using that statistic to validate them would be
circular. The design here breaks the circle three ways:

1. **Held-out genome.** Windows are split into two disjoint halves by
   alternating chromosome. Candidates are selected on half A, using the
   pipeline's own rule (correlation above the 95th percentile of that sample's
   reference cells). Every reported number is computed on half B.
2. **Own-patient specificity, not correlation.** The half-B statistic is the
   cell's correlation to its own sample's consensus *minus* its mean correlation
   to the consensus of the other samples in the same study. Ambient RNA,
   dissociation stress and low complexity raise correlation to any consensus;
   only a shared clonal karyotype raises it to one patient in particular.
3. **Consensus from marker-defined cells only.** The candidates never
   contribute to the profile they are scored against.

A fourth control was needed and changed the conclusion. Selecting *any* cells on
a correlated noise statistic inflates the held-out statistic: the half-A
selection raises half-B specificity by a median of +0.041 in the candidates and
by exactly the same +0.041 in known-diploid reference cells
(`f4_selection_carryover.csv`). The *gain* from selection is therefore
uninformative, and only the absolute level can be interpreted.

Absolute level is confounded by depth — candidates carry 12,472 median UMI
against 2,467 in reference cells — so the reported run downsamples every cell to
exactly 1,500 counts and discards cells below it (76,475 of 138,353 retained).

## Result

At matched depth, cell-level mean own-patient specificity on the held-out half:

| group | mean spec | n |
|---|---|---|
| diploid reference | 0.083 | 14,113 |
| fibroblast-like, copy number negative | 0.121 | 2,932 |
| fibroblast-like, copy number positive | 0.229 | 761 |
| marker-defined malignant | 0.296 | 39,707 |

The state is genetically real: the candidates sit between stroma and tumour and
nearer the tumour, they are not doublets (doublet score 0.029, the lowest of any
group, against 0.092 in marker-defined malignant cells), and they reproduce the
published CD73 marker.

It is also rare and unevenly distributed. Only 8 of 45 scored samples carry ten
or more candidate cells, in 5 patients: Ewing_2 (Cillo), Ewing_Patient (Mao),
ES-025 (Visser), SJEWS001321 (Waltner) and UC_Patient — the urothelial carcinoma
control of Mao et al., where the marker-defined malignant cells are urothelial
carcinoma rather than Ewing sarcoma, and which contributes two of the eight
samples (CR-hPB136, CR-hPB137). Four of the eight samples are Mao.

Carrying candidate cells is not the same as passing the test. Specificity is
positive and above the sample's own reference cells in 6 of the 8 samples and 4
of the 5 patients; the two plate-based Visser samples fail (TM416 −0.041, TM417
−0.083, both below their own reference), so ES-025 drops out. Excluding the
urothelial carcinoma control, **three** Ewing sarcoma patients of 38 carry
positive evidence: Ewing_2, Ewing_Patient and SJEWS001321.

## Conclusion

Both earlier readings were wrong in part. The rescued fibroblast-like cells are
not uniformly a rescue artefact — in a minority of samples they carry their own
patient's karyotype on a held-out half of the genome. But they are not a
recurrent tumour state across the cohort either, and the cohort cannot support a
claim about treatment, site or outcome from three patients. The marker-based
annotation and the frozen malignant set are unchanged by this result; what
changes is the interpretation of the `cna_only` evidence class, which should be
described as heterogeneous rather than as contamination.

## Reproduce

```
python src/pipeline/07_caf_state.py          # native depth
python src/pipeline/07_caf_state.py 1500     # depth-matched
python src/pipeline/07_caf_state.py tables   # panel inputs
python src/pipeline/06_figures.py
```
