# 14. Deposition and reuse

Stage: `src/pipeline/13_cellxgene_export.py`
Outputs: `data/processed/EWSAtlas_v5_cellxgene.h5ad`, `deposit/`

## What exists

| artefact | path | state |
|---|---|---|
| analysis object | `data/processed/EWSAtlas_v5.h5ad` | ready (138,353 x 16,622) |
| CELLxGENE deposit object | `data/processed/EWSAtlas_v5_cellxgene.h5ad` | **passes cellxgene-schema 7.1.2** (138,353 x 16,285, 1.10 GB) |
| trained scVI model | `data/interim/scvi_model_v5/model.pt` | ready (377 epochs, early stopping) — see the caveat below |
| query-mapping utility | `src/map_query_to_ewsatlas.py` | repointed at v5 and at the model's own latent |
| Zenodo metadata | `deposit/zenodo.json` | needs creators and license |
| citation file | `deposit/CITATION.cff` | needs authors, date, license |

## The saved model does not reproduce the released embedding

`03b_train_scvi.py` trains the configuration in `src/config.json` and stops by
early stopping at 377 of 500 epochs. Its latent agrees with the released
`X_scvi` globally but not locally: Procrustes-aligned cosine 0.776, Spearman
correlation of pairwise distances 0.904, but only **34.7%** of each cell's 15
nearest neighbours are shared (5,000 cells compared). The released embedding was
produced interactively before this stage existed, so the gap cannot be
reconciled after the fact.

Replacing `X_scvi` would invalidate the UMAP, the 34 Leiden communities and the
annotation built on them — that is, every figure. The retrained coordinates are
therefore stored alongside it as `obsm['X_scvi_model_v5']`, with the agreement
numbers in `uns['scvi_model_v5']`, and `map_query_to_ewsatlas.py` transfers labels
in **that** space. This is correct rather than a workaround: cell-type labels are
properties of cells, not of an embedding, so a query cell projected by the model
must find its reference neighbours in the space the model itself produces. Anyone
re-deriving the embedding from scratch should expect the same 0.3-0.4
neighbourhood overlap and should not read it as a defect in either run.

## The deposit object

Gene identifiers are Ensembl, fetched from Ensembl REST: 16,286 of 16,622
symbols (98.0%) resolve, and one of those, `ENSG00000285101`, is absent from the
GENCODE reference the validator uses, leaving **16,285** genes. All 337 dropped
symbols remain in the analysis object. `X` is log-normalized and `raw.X` is
counts, both `float32`, as the schema requires.

Ontology terms were resolved against EBI OLS and carry the label OLS returned;
none was written from memory. Cell types map to CL, assay to EFO by recorded
chemistry, sex to PATO, disease to MONDO (`MONDO:0012817` Ewing sarcoma;
`MONDO:0040679` urothelial carcinoma for the control patient; `PATO:0000461`
normal for tumour-adjacent bladder); organism is `NCBITaxon:9606`, carried in
`uns` as schema 7 requires.
`suspension_type` is `nucleus` for the Waltner libraries and `cell` elsewhere,
and `is_primary_data` is `False` throughout because all six studies are published
elsewhere.

### Remaining unknown fields

Three remain `unknown`, which schema 7 permits for these fields; the fourth,
`tissue_ontology_term_id`, no longer may be and has been resolved.

1. `development_stage_ontology_term_id` — 138,353 of 138,353. Needs one HsapDv
   term per age; 42,230 cells carry no age at all (the same 42,230 that carry no
   sex — deposits that omit one omit the other). This is the single largest
   available improvement to the record.
2. `self_reported_ethnicity_ontology_term_id` — 138,353 of 138,353. No source
   deposit reports it and it must not be inferred.
3. `sex_ontology_term_id` — 42,230 of 138,353, and `cell_type_ontology_term_id`
   for the 2,809 `Unassigned` cells. Both are genuinely missing rather than
   unmapped.
4. `tissue_ontology_term_id` — resolved for all cells. Circulating-tumour-cell
   and blood libraries map to blood rather than to the primary site, because the
   field describes what was sampled; chest wall and thorax to the thoracic
   segment of the trunk, plus axilla, ileum and brain; and 14,757 cells whose
   deposit reports no site carry the generic ancestor `UBERON:0000479`
   ("tissue"). Every term was checked against the ontology shipped with
   cellxgene-schema before use.

## Schema 7 validation

The deposit object passes `cellxgene-schema` 7.1.2 with `is_valid=True`
(138,353 x 16,285, 1.10 GB). Getting there required six changes beyond the
original export, all now in stage 13 and each one a rule that is easy to miss:

1. Discover annotates the human-readable mirror of every ontology field and the
   `feature_*` columns itself; submitting a copy of `assay`, `cell_type`, `sex`
   or `feature_name` fails validation. They are dropped and preserved as
   `cell_type_atlas`, `assay_atlas`, `sex_reported` and `gene_symbol_atlas`.
2. Schema 7 moved `organism_ontology_term_id` from `obs` to `uns`.
3. `var['feature_is_filtered']` is required; it is `False` throughout, because
   excluded genes were removed from the matrix rather than zeroed in it.
4. Both matrices must be `float32`; the analysis object holds int32 counts and a
   float64 normalized layer.
5. One Ensembl identifier the symbol map resolves, `ENSG00000285101`, is absent
   from the GENCODE reference the validator uses and is dropped.
6. `cellxgene-schema` pins an older anndata than this pipeline writes with, and
   that reader cannot open the nullable encodings anndata 0.13 emits by default.
   The stage sets `allow_write_nullable_strings = False`, casts nullable columns
   and strips `None` from `uns`. A file the validator cannot open cannot be
   ingested, so this is a blocker rather than a nicety.

Dropping the duplicated `layers` (`counts` is `raw.X`, `lognorm` is `X`) and
the float32 cast together took the object from 2.46 GB to 1.10 GB.

The collection form, the submission email and the curation decisions a reviewer
should see are in `deposit/cellxgene_collection.md`.

## Two steps that need your account

- **Zenodo DOI.** `deposit/zenodo.json` is complete except for creators and
  licence. Uploading needs a Zenodo token, which this project does not hold.
- **Public code repository.** The repository is ready to publish, but nothing has
  been pushed anywhere: publishing is a decision with no undo, so it waits for an
  explicit instruction, and a private repository first is the safer order.

A licence has deliberately **not** been chosen. CC-BY-4.0 for the data and a
permissive licence for the code is the usual pairing for a derived atlas, but the
source deposits' own terms should be checked first — ScPCA in particular.
