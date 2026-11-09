# CELLxGENE Discover — collection metadata and submission record

Prepared for submission of `data/processed/EWSAtlas_v5_cellxgene.h5ad`.
This file is the copy of record for what is entered into the Collection form;
the object itself is produced by `src/pipeline/13_cellxgene_export.py`.

## Validation status

```
cellxgene-schema 7.1.2
cellxgene-schema validate data/processed/EWSAtlas_v5_cellxgene.h5ad
Validation complete with status is_valid=True
```

| | |
|---|---|
| object | `EWSAtlas_v5_cellxgene.h5ad` |
| shape | 138,353 cells x 16,285 genes |
| size | 1,097,994,570 bytes (1.10 GB) |
| sha256 | `ecb8ba797ed998d89bb2a893e8b66640c190817c2cc16e3f4930f43fb2638891` |
| schema | CELLxGENE Discover dataset schema 7.x |

One warning is raised and is expected: *"Dataframe 'var' only has 16285 rows.
Features SHOULD NOT be filtered from expression matrix."* The atlas gene space
is the intersection of six per-study annotations (`docs/03_gene_space.md`), so
the unmeasured genes were never in the matrix to begin with and cannot be
restored with `feature_is_filtered`. This is a SHOULD, not a MUST.

## The process, per the official documentation

Source: <https://cellxgene.cziscience.com/docs/032__Contribute%20and%20Publish%20Data>

This is **not** a self-service upload. The documented order is:

1. Review the Data Eligibility criteria.
2. **Contact CZI first** with a description of the data, to confirm they will
   accept it.
3. Once confirmed, send files prepared to the submission requirements.
4. **CZI uploads to a private Collection** where you review it.
5. Revise as needed.
6. The data are made openly available when you are ready.

So the collection metadata below is what goes **in the email**, not into a form
you fill in yourself, and nothing becomes public until you say so.

### Eligibility, checked against the criteria

| criterion | EWSAtlas |
|---|---|
| not a drug screen | passes |
| species on the supported list | human only (`NCBITaxon:9606`) |
| assays on the accepted list | all three listed: `EFO:0009922` 10x 3' v3, `EFO:0009899` 10x 3' v2, `EFO:0010010` CEL-seq2 |
| max 50 GB per dataset | 1.10 GB |
| under 4.3M cells (else not visualized in Explorer) | 138,353 cells |
| single h5ad in **AnnData 0.10 format** | yes — and this is why stage 13 forces pre-0.11 encodings |

### Data Submission Policy — what you are agreeing to

Submitting asserts that you have authority to let CZI display, distribute and
create derivative works from the data; that **you** are responsible for ensuring
the data are not identifiable and for removing direct personal identifiers from
the metadata; that published data are available under **CC-BY 4.0**; and that
the Collection details, including your name and contact information, become
public. Published data can be deleted later, which does not undo prior
downloads.

De-identification was checked on the deposit object. Every identifier column
holds a study-internal code that is already public in the source deposit —
`donor_id` and `patient` are three-digit codes, `src_participant_id` holds the
SJEWS sample codes, `sample` and `src_submitter_id` combine those with a
suffix, and `age_years` is an integer year. No name, date, location or
free-text clinical note is carried. `src_submitter` retains one depositor
handle from the source metadata and may be dropped if you prefer.

## Collection metadata for the email

**Name**
EWSAtlas: a harmonized single-cell and single-nucleus reference for Ewing sarcoma

**Description**
EWSAtlas assembles the public single-cell and single-nucleus transcriptomic
record for EWSR1::FLI1 Ewing sarcoma into one reference: 138,353 cells from 65
samples and 38 donors across six published cohorts, processed under one
quality-control policy, one gene space, one controlled vocabulary and one
annotation pipeline. Cells are labelled with 14 cell types assigned per Leiden
community against obligatory marker panels, with communities meeting no label's
evidence left unassigned rather than forced to the nearest label. The malignant
compartment is defined by markers together with a per-sample consensus
copy-number statistic, and every malignant cell carries the evidence class that
placed it there. Batch correction is scVI with the sample as the batch, chosen
by an scIB benchmark against study-level scVI, Harmony and unintegrated
principal components. Covariates include donor, source material (patient
tumour, orthotopic xenograft, circulating tumour cells, blood, adjacent normal),
assay, treatment timing, disease state and fusion partner. The collection
includes one urothelial-carcinoma control patient from the bladder cohort,
carrying its own disease term so that it can be excluded.

**Contact**  Md. Jubayer Hossain — contact.jubayerhossain@gmail.com

**Publication / preprint DOI**  not yet posted; cite the data DOI below

**Links**
- Data source (archived release): https://doi.org/10.5281/zenodo.23622732
- Code: https://github.com/hossainlab/ews_atlas

**License**  CC BY 4.0

**Consortia**  none

## Dataset

Title (from `uns['title']`):
`EWSAtlas v5: a harmonized single-cell reference for Ewing sarcoma`

Default embedding `X_umap`; `X_scvi` and `X_scvi_model_v5` also present.
Queries must be mapped into `X_scvi_model_v5`, the saved model's own latent
space, not into the released annotation embedding — see
`docs/14_deposition.md`.

## Source deposits

| study | accession |
|---|---|
| Aynaud 2020 | GEO GSE130024 |
| Cillo 2022 | GEO GSE198896 |
| Visser 2023 | GEO GSE243347 |
| Goodspeed 2025 | GEO GSE277083 |
| Mao 2024 | GEO GSE250523 |
| Waltner 2025 | ScPCA SCPCP000015 |

Before signing the submission warranty, confirm redistribution terms for the
ScPCA arm specifically: a portal-hosted dataset carries its own terms of use,
and a journal's CC BY on an article does not automatically extend to data
hosted elsewhere.

## Curation decisions a reviewer should see

1. **`is_primary_data` is `False` for all 138,353 cells.** Every cohort is
   reused from another deposit. The consequence is that the dataset will not
   contribute primary cells to the CELLxGENE Census. Do not flip this flag.
2. **Orthotopic xenografts are `tissue_type = 'tissue'`.** The enum offers
   `primary cell culture`, `organoid`, `tissue` and `cell line`; a PDX is none
   of them exactly, and `tissue` is the closest, since those cells grew in a
   host tissue rather than in culture. 56,267 cells are affected. Flag this to
   the curation team and change it if they prefer another value.
3. **Circulating-tumour-cell and blood libraries are mapped to blood**
   (UBERON:0000178, 8,533 cells) rather than to the primary tumour's site,
   because the tissue field describes what was sampled.
4. **14,757 cells carry the generic term `UBERON:0000479` ("tissue")** because
   their deposits report no site (Aynaud's xenografts) or report one that does
   not resolve. Schema 7 forbids `unknown` when `tissue_type` is `tissue`, so
   the generic ancestor is used rather than a guessed specific site.
5. **`development_stage` and `self_reported_ethnicity` are `unknown`
   throughout.** Ethnicity is reported by no source deposit and must not be
   inferred. Development stage is partly recoverable — age is present for
   96,123 of 138,353 cells — and filling HsapDv terms is the single largest
   available improvement to this record.
6. **`sex_ontology_term_id` is `unknown` for 42,230 cells**, the same cells
   that carry no age: the deposits that omit one omit the other.
7. **`cell_type_ontology_term_id` is `unknown` for 2,809 cells**, the
   `Unassigned` class. These are communities that met no label's marker
   evidence, left unassigned by design.
8. **The two malignant classes map to `CL:0001064` (malignant cell).** The
   atlas distinction between cycling and non-cycling tumour cells is preserved
   in `obs['cell_type_atlas']`, which is not a schema field.

## Fields deliberately absent

CELLxGENE Discover annotates these itself, and submitting them fails
validation: `obs` `assay`, `cell_type`, `development_stage`, `disease`,
`organism`, `self_reported_ethnicity`, `sex`, `tissue`; `var` `feature_name`
and the other `feature_*` label columns; `uns` `schema_version` and
`schema_reference`. The atlas values survive as `obs['cell_type_atlas']`,
`obs['assay_atlas']`, `obs['sex_reported']` and `var['gene_symbol_atlas']`.

## Step-2 email: ask whether the data will be accepted

This is step 2 of the process: a description sent to ask whether the data will
be accepted. Do not attach the h5ad yet.

> To: cellxgene@chanzuckerberg.com
> Subject: Proposed contribution — EWSAtlas, a harmonized single-cell reference for Ewing sarcoma
>
> Dear CELLxGENE curation team,
>
> I would like to confirm whether you would accept a dataset for CZ CELLxGENE
> Discover.
>
> Title: EWSAtlas: a harmonized single-cell and single-nucleus reference for
> Ewing sarcoma
> Description: [the Description above]
> Contact: Md. Jubayer Hossain, contact.jubayerhossain@gmail.com
> Publication/preprint DOI: not yet posted; archived release at
> https://doi.org/10.5281/zenodo.23622732
> URLs: code at https://github.com/hossainlab/ews_atlas
> Consortia: none
>
> The dataset is one h5ad of 138,353 cells x 16,285 genes (1.10 GB), human
> only, assays 10x 3' v2, 10x 3' v3 and CEL-seq2, all on the accepted list. It
> is an integration of six already-published cohorts, so is_primary_data is
> False throughout. It passes cellxgene-schema 7.1.2 with is_valid=True; the
> only warning is that var has 16,285 rows, because the gene space is the
> intersection of the six source annotations and the unmeasured genes were
> never in the matrix.
>
> Two curation points I would value your guidance on before sending files:
> first, the tissue_type value for orthotopic patient-derived xenografts
> (56,267 cells), which fit none of tissue / organoid / primary cell culture /
> cell line exactly and are currently 'tissue'; second, the use of the generic
> UBERON:0000479 for 14,757 cells whose source deposit reports no anatomical
> site, since 'unknown' is not permitted when tissue_type is 'tissue'.
>
> Kind regards,
> Md. Jubayer Hossain
