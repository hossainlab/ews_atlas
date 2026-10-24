# 00 · Project overview

## What this is

EWSAtlas harmonizes the public single-cell and single-nucleus transcriptomic
data for Ewing sarcoma into one object with one QC policy, one feature space,
one batch-corrected embedding and one cell-type vocabulary. NBAtlas
(Bonine et al., *Cell Reports* 43:114804, 2024) is the template: a disease-focused
reference built mostly from already-public data, whose value is harmonization
and reuse rather than new sequencing.

## Scope

**Included:** EWSR1::FLI1 (and EWSR1::ERG) Ewing sarcoma — tumour, xenograft,
circulating tumour cells, and matched blood/normal tissue from the same studies.

**Excluded, deliberately.** Other FET-rearranged tumours are a different
disease despite the shared fusion grammar, and pooling them would make every
malignant-compartment statement ambiguous:

| entity | fusion | why excluded |
|---|---|---|
| DSRCT | EWSR1::WT1 | distinct lineage and transcriptional output |
| clear cell sarcoma | EWSR1::ATF1 | melanocytic differentiation |
| intracranial FET-CREB tumours | EWSR1/FUS::CREB family | distinct entity |
| FET-rearranged myoepithelial tumours | EWSR1::POU5F1 etc. | distinct entity |
| myxoid liposarcoma | FUS::DDIT3 | adipocytic lineage |

Cell lines and pooled cell-line experiments are excluded from the atlas object;
one cell-line resource (ESCLA) is used as an external *reference* for the
fusion-activity axis, not as atlas content (see [07](07_bulk_validation.md)).

## What the atlas is not

- **Not an outcome resource.** No source deposit carries survival. Outcome
  association is done in bulk cohorts and is clearly separated from the atlas
  object.
- **Not a size record.** The claim is coverage — exhaustive of locatable public
  EwS single-cell data — and cross-study harmonization, not cell count.
- **Not a cell-of-origin study.** That question is contested and the data here
  cannot settle it.

## Current status

| object | state |
|---|---|
| EWSAtlas v4 (5 studies, 64,702 cells × 17,017 genes) | built; superseded |
| EWSAtlas v5 (6 studies, adds the Waltner snRNA arm) | in progress |

v3 (129,276 cells) was deleted: it merged two incompatible metadata schemas,
mislabelled xenograft status for 14,702 cells, and carried an all-`False`
doublet column while the doublet caller had flagged 7,571 cells. Its cell count
is not comparable to v4/v5 — see [decisions.md](decisions.md), ADR-002.
