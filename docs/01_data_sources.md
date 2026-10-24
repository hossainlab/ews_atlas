# 01 · Data sources

Six deposits, six different on-disk shapes. `src/ewsatlas/io.py` has one reader
per shape; `src/pipeline/01_ingest_qc.py` maps study → reader in `DEPOSITS`.

| study | accession | assay | material | on-disk shape |
|---|---|---|---|---|
| Aynaud 2020 | GSE130024 | 10x 3′ scRNA | 5 orthotopic PDX | MTX triplets in one RAW.tar |
| Cillo 2022 | GSE198896 | 10x scRNA | TIL + PBMC | per-library nested tar.gz |
| Goodspeed 2025 | GSE277083 | 10x scRNA | tumour + CTC | dense gene × cell TSV + metadata |
| Mao 2024 | GSE250523 | 10x scRNA | bladder EwS + urothelial control | MTX triplets |
| Visser 2023 | GSE243347 | CEL-seq2 | patient tumour | gene × well plate tables |
| Waltner 2025 | SCPCP000015 | 10x snRNA | patient tumour | per-library AnnData in one zip |

Primary citations: Aynaud 10.1016/j.celrep.2020.01.049 · Cillo
10.1158/1078-0432.CCR-22-1471 · Goodspeed 10.1158/1078-0432.CCR-24-2040 · Mao
10.1016/j.isci.2024.110921 · Visser 10.1158/2767-9764.CRC-23-0027 · Waltner
10.7554/eLife.109124.

## Traps, each of which cost a rebuild

**Aynaud GSE130024 ships two incompatible gene namings.** Four samples carry
CellRanger v2 `genes.tsv` with `GRCh38_`-prefixed symbols; PDX-1058 ships plain
v3 `features.tsv`. Intersecting the feature spaces without stripping the prefix
matches only PDX-1058 and silently zeroes 13,064 of 15,057 cells. The
downstream signature is unmistakable once you know it: hundreds of uniform
~50-cell leiden communities from a single study, plus a **0.0% mitochondrial
fraction** (the `MT-` genes were never matched). After stripping, Aynaud's
median mitochondrial content is 3.7%. Guarded by
`tests/test_invariants.py::test_symbol_harmonization_strips_build_prefix` and
four companion tests; removing the strip fails all five.

**Visser's deposit is richer than its paper.** `GSE243347_data_annotation.txt.gz`
resolves pooled plate **TM736** at well level (columns 2–3 = patient ES-039,
columns 4–18 = ES-030) and documents **plate TM768 (patient ES-048, lymph-node
metastasis)**, which earlier builds omitted entirely. Always reconcile plate
assignment against this file, not the manuscript tables.

**Visser's matrices are fractional, not counts.** 1.8% of nonzero entries are
non-integral (median nonzero 1.0, maximum 29,569.7) because CEL-seq2 UMI counts
are collision-corrected. `io.read_celseq_plates` rounds at ingest — a negative
binomial likelihood is undefined on fractional input — and keeps the corrected
values in `layers["celseq_corrected"]`. Guarded by
`test_celseq_reader_rounds_to_counts`.

**Cillo GSE198896 contains out-of-scope libraries.** The series includes
osteosarcoma and healthy-donor samples. Only the five Ewing libraries
(GSM5959129/130/139/140/141 → Ewing_1–4) are ingested; `GSM5959129` is PBMC,
not TIL.

**Goodspeed ships primary and blood matrices separately**; both are read and
concatenated, with the source cell metadata retained as `src_*` columns for
cross-checking but never as the atlas's own labels.

**Waltner is the only single-nucleus arm** and arrives through the ScPCA
portal, not GEO. ScPCA's emptyDrops-filtered matrices are the starting point;
our QC is then re-applied so the nuclei pass the same policy as the five GEO
studies. ScPCA's own annotations are kept only as `src_*` provenance.

**ScPCA's sample metadata has an internal ID conflict.** For two libraries the
`submitter_id` prefix disagrees with `participant_id`:

| library | submitter_id | participant_id | clinical block matches |
|---|---|---|---|
| SCPCL001114 | SJEWS063826_X1 | SJEWS063829 | `participant_id` (fusion unknown, initial diagnosis, primary, femur) |
| SCPCL001115 | SJEWS063829_X1 | SJEWS063834 | `participant_id` (EWSR1-FLI1, recurrence, metastasis, chest wall) |

The conflict is a one-position shift along the patient list. In both cases three
independent clinical fields travel with `participant_id`, not with the
`submitter_id` prefix, so `participant_id` is used and the conflict is flagged in
`results/tables/scpca_library_profile.csv` (`id_conflict`). This cannot be
settled from metadata alone — the decisive test is whether each PDX library
shares its putative patient's copy-number profile, which the CNA step can check
directly. Until then, no per-patient claim should rest on these two libraries.

**The ScPCA arm supplies the matched pairs v4 lacked.** 13 libraries from 5
patients, all single-nucleus 10x v3.1: 7 patient tumours and 6 orthotopic PDX,
with 4 patients contributing both. That is what makes the PDX-versus-patient
question answerable at all — in v4 the entire PDX stratum was one study.

## Integrity

Every input is fingerprinted by sha256 in the stage manifest at read time. The
ScPCA archive `SCPCP000015_ann-data_2026-10-17.zip` (2.67 GB) has sha256
beginning `322a6152364b7b88`.

## Deposits examined and rejected

Triaged as not integrable, with the reason: GSE318841 / GSE319124 / GSE330883
(bulk patient-matched PDX/organoid panels); GSE314150 / 314152 / 314154 (STAG2 ×
EWS-FLI1 in A673/TC71 — cell lines, ATAC/ChIP/bulk); GSE165973 (CITE-seq of
pooled A673 + CHLA10 — cell lines); GSE130020 (96 marrow MSCs); GSE130023
(142-cell C1 PDX arm). GSE261693 (3 human GEX libraries from a humanized-mouse
model) is the one remaining genuine single-cell candidate, deferred.

Two published cohorts could not be obtained at all: **Hong 2022**
(10.1007/s00432-022-04073-3, 39,630 cells, 3 patients) and **Chen 2025**
(10.1186/s12967-025-06346-0, 13 immune-sorted bone tumours) deposit no
accession in GEO, SRA, ENA, ArrayExpress, EGA or dbGaP. Their absence is the
reason the atlas has **no external query cohort**: reference mapping here can
demonstrate projection stability, not generalization to unseen data.
