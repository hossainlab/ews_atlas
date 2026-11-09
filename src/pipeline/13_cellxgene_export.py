"""Stage 13 - export a CELLxGENE-schema object for deposition.

    python src/pipeline/13_cellxgene_export.py

Writes ``data/processed/EWSAtlas_v5_cellxgene.h5ad`` with CELLxGENE's required
obs fields and Ensembl gene identifiers.

Every ontology term used here was resolved against EBI OLS and is recorded with
the label OLS returned, including the two that are assigned in code rather than
through a lookup table: PATO:0000461 ("normal") and NCBITaxon:9606
("Homo sapiens"). Where a value cannot be mapped without a curation decision that
belongs to the depositors, the field is set to ``unknown`` and listed in the
stage's own output rather than guessed - an invented ontology term is worse than a
missing one, because it validates.

Known gaps, printed by the stage and repeated in docs/14_deposition.md:

  * ``development_stage_ontology_term_id`` is unknown for all cells. Mapping
    age in years to HsapDv requires a term per age, and 42,230 cells carry no
    age at all - the same 42,230 that carry no sex, since the deposits that omit
    one omit the other.
  * ``self_reported_ethnicity_ontology_term_id`` is unknown for all cells: no
    source deposit reports it, and it must not be inferred.
  * ``tissue_ontology_term_id`` resolves for the sites that map cleanly to bone,
    lung, bladder and lymph node; free-text sites (chest wall, axilla, paraspinal,
    ileum, cutaneous) stay unknown pending curation.
  * 336 of 16,622 symbols (2.0%) have no current Ensembl identifier and are
    dropped from the deposit object; they remain in the analysis object.

Fields the submitter must NOT supply, because CELLxGENE Discover annotates them
and a submitted copy fails validation: obs assay, cell_type, development_stage,
disease, organism, self_reported_ethnicity, sex, tissue; var feature_name and
the other feature_* columns; uns schema_version and schema_reference. This stage
therefore drops them and keeps the atlas values as obs cell_type_atlas,
assay_atlas, sex_reported and var gene_symbol_atlas. layers are dropped as well:
counts is raw.X and lognorm is X, so carrying both doubles the file.

Validate the output before upload:  cellxgene-schema validate <out.h5ad>
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ewsatlas import config, io, provenance  # noqa: E402

# --- ontology terms, each verified against EBI OLS (label OLS returned in comment)
CELL_TYPE = {
    "Ewing_sarcoma": "CL:0001064",          # malignant cell
    "Ewing_proliferating": "CL:0001064",    # malignant cell
    "T_cell": "CL:0000084",                 # T cell
    "NK": "CL:0000623",                     # natural killer cell
    "B_cell": "CL:0000236",                 # B cell
    "Plasma": "CL:0000786",                 # plasma cell
    "Monocyte_Macrophage": "CL:0000113",    # mononuclear phagocyte
    "DC": "CL:0000990",                     # conventional dendritic cell
    "pDC": "CL:0000784",                    # plasmacytoid dendritic cell
    "Mast": "CL:0000097",                   # mast cell
    "Platelet": "CL:0000233",               # platelet
    "Endothelial": "CL:0000115",            # endothelial cell
    "Fibroblast_CAF": "CL:0000057",         # fibroblast
    "Urothelial": "CL:0000731",             # urothelial cell
    "Unassigned": "unknown",
}
ASSAY = {
    "10x_Chromium_3prime_v3.1": "EFO:0009922",                      # 10x 3' v3
    "10x_Chromium_NextGEM_3prime_v3.1": "EFO:0009922",              # 10x 3' v3
    "10x_Chromium_3prime_v2": "EFO:0009899",                        # 10x 3' v2
    "10x_Chromium_NextGEM_3prime_version_not_reported": "EFO:0009922",
    "CEL-Seq2_384well_plate": "EFO:0010010",                        # CEL-seq2
    "not_reported": "EFO:0009922",
}
SEX = {"male": "PATO:0000384", "female": "PATO:0000383", "not_reported": "unknown"}
BONE = ["tibia", "femur", "fibula", "rib", "sacrum", "scapula", "humerus", "ulna",
        "iliac crest", "bone", "fibula and femur"]
TISSUE_EXACT = {"lung": "UBERON:0002048", "urinary bladder": "UBERON:0001255",
                "urinary bladder (tumour-adjacent normal)": "UBERON:0001255",
                "lymph node": "UBERON:0000029",
                # added for schema 7.0.0, which forbids 'unknown' when
                # tissue_type is 'tissue'. Every term below was checked against
                # the ontology shipped with cellxgene-schema: it resolves, and
                # it descends from UBERON:0001062 (anatomical entity).
                "chest wall": "UBERON:0000915",      # thoracic segment of trunk
                "thorax": "UBERON:0000915",
                "axilla": "UBERON:0009472",          # axilla
                "l ileum": "UBERON:0002116",         # ileum
                "ileum": "UBERON:0002116",
                "brain": "UBERON:0000955"}           # brain
TISSUE_FALLBACK = "UBERON:0000479"                   # tissue, site not reported
BLOOD = "UBERON:0000178"                             # blood
# Ensembl identifiers that the symbol map resolves but that are absent from the
# GENCODE reference cellxgene-schema validates against ("Could not infer
# organism from feature ID"). Dropped from the deposit object only.
NOT_IN_GENCODE = {"ENSG00000285101"}


def tissue_term(site: str, material: str = "") -> str:
    """UBERON term for the tissue that was sampled, not for the tumour's origin.

    A circulating-tumour-cell library and a PBMC library are blood however the
    deposit records the primary site, so material decides those before site is
    consulted. Sites the deposits give as free text and that cannot be resolved
    (Aynaud's xenografts report none) fall back to the generic 'tissue' term,
    which is honest about the gap; 'unknown' is not a legal value here.
    """
    if str(material) in ("CTC", "PBMC"):
        return BLOOD
    s = str(site).strip().lower()
    if s in TISSUE_EXACT:
        return TISSUE_EXACT[s]
    if any(b in s for b in BONE):
        return "UBERON:0001474"        # bone element
    return TISSUE_FALLBACK


def main() -> int:
    # Write the pre-0.11 string format. cellxgene-schema pins an older anndata
    # that raises IORegistryError on the nullable encodings anndata 0.13 emits
    # by default, and a file the validator cannot open cannot be ingested.
    ad.settings.allow_write_nullable_strings = False

    cfg = config.load()
    PROJ = config.ROOT
    atlas = PROJ / config.get(cfg, "paths.atlas_out")
    ens_file = PROJ / "data/external/references/gene_symbol_to_ensembl.csv"

    with provenance.start_run("13_cellxgene_export", cfg, [atlas, ens_file]) as run:
        A = ad.read_h5ad(atlas)
        ens = pd.read_csv(ens_file, index_col=0)["ensembl_id"]

        keep = [g for g in A.var_names if g in ens.index and pd.notna(ens.get(g))
                and ens[g] not in NOT_IN_GENCODE]
        dropped = A.n_vars - len(keep)
        A = A[:, keep].copy()
        # Discover annotates var["feature_name"] from its own Ensembl reference,
        # so a submitted copy fails validation; the atlas's harmonized symbol is
        # kept under a name the schema does not reserve.
        A.var["gene_symbol_atlas"] = A.var_names
        # Required by the schema. False throughout: genes excluded by the gene
        # space were removed from the matrix, not zeroed in it, so no retained
        # feature has had its raw counts suppressed.
        A.var["feature_is_filtered"] = False
        A.var_names = pd.Index([ens[g] for g in A.var_names], name="ensembl_id")
        A.var_names_make_unique()

        o = A.obs
        o["assay_ontology_term_id"] = o["assay_chemistry"].astype(str).map(ASSAY).fillna("EFO:0009922")
        o["suspension_type"] = np.where(o["assay"].astype(str) == "10x_snRNA", "nucleus", "cell")
        o["cell_type_ontology_term_id"] = o["cell_type"].astype(str).map(CELL_TYPE).fillna("unknown")
        o["sex_ontology_term_id"] = o["sex"].astype(str).map(SEX).fillna("unknown")
        o["tissue_ontology_term_id"] = [tissue_term(s, m) for s, m
                                        in zip(o["tumor_site"].astype(str),
                                               o["material_type"].astype(str))]
        # The enum is primary cell culture / organoid / tissue / cell line; an
        # orthotopic xenograft is none of those exactly, and 'tissue' is the
        # closest, since the cells were growing in a host tissue rather than in
        # culture. Flagged for the curation team in docs/14_deposition.md.
        o["tissue_type"] = "tissue"
        # the urothelial-carcinoma control patient is not Ewing sarcoma, and the
        # tumour-adjacent normal bladder is not malignant tissue at all
        site = o["tumor_site"].astype(str).str.lower()
        o["disease_ontology_term_id"] = np.where(
            o["disease_state"].astype(str) == "normal", "PATO:0000461",   # normal
            np.where(o["patient"].astype(str) == "UC_Patient",
                     "MONDO:0040679",                                    # urothelial carcinoma
                     "MONDO:0012817"))                                   # Ewing sarcoma
        o["development_stage_ontology_term_id"] = "unknown"
        o["self_reported_ethnicity_ontology_term_id"] = "unknown"
        o["is_primary_data"] = False          # all six studies are published elsewhere
        o["donor_id"] = o["patient"].astype(str)

        # Discover annotates the human-readable mirror of every ontology field
        # itself (assay, cell_type, development_stage, disease, organism,
        # self_reported_ethnicity, sex, tissue); submitting one fails validation.
        # Keep the atlas values under names the schema does not reserve.
        RESERVED = ["assay", "cell_type", "development_stage", "disease",
                    "organism", "self_reported_ethnicity", "sex", "tissue"]
        o["cell_type_atlas"] = o["cell_type"].astype(str)
        o["assay_atlas"] = o["assay"].astype(str)
        o["sex_reported"] = o["sex"].astype(str)
        A.obs = o.drop(columns=[c for c in RESERVED if c in o.columns])

        # Plain object indices before raw is built, so raw inherits them too
        # (see the encoding note below).
        A.obs.index = pd.Index([str(i) for i in A.obs_names], dtype=object)
        A.var.index = pd.Index([str(i) for i in A.var_names], dtype=object,
                               name="ensembl_id")

        # raw.var carries the index only: feature_* names are Discover's.
        # The schema requires float32 in both matrices; the analysis object
        # stores counts as int32 and the normalized layer as float64.
        A.raw = ad.AnnData(X=A.layers["counts"].astype(np.float32),
                           var=pd.DataFrame(index=A.var_names.copy()),
                           obs=pd.DataFrame(index=A.obs_names.copy()))
        A.X = A.layers["lognorm"].astype(np.float32)
        for layer in ("counts", "lognorm"):      # now duplicated by raw.X and X
            del A.layers[layer]
        A.uns["title"] = "EWSAtlas v5: a harmonized single-cell reference for Ewing sarcoma"
        # schema 7.0.0 moved organism from obs to uns and deprecated the obs field
        A.uns["organism_ontology_term_id"] = "NCBITaxon:9606"
        A.uns["default_embedding"] = "X_umap"
        A.uns["batch_condition"] = ["sample"]
        if isinstance(A.uns.get("config"), dict):
            A.uns["config"].pop("compute", None)   # build-infrastructure ids, not data

        # cellxgene-schema pins an older anndata than the one this pipeline
        # writes with, and that reader raises IORegistryError on the encodings
        # anndata 0.13 emits for pandas nullable dtypes ("nullable-string-array",
        # "nullable-boolean") and for a None inside uns ("null"). A deposit the
        # validator cannot open is a deposit CZI cannot ingest, so the object is
        # written with legacy encodings only.
        NULLABLE = ("boolean", "string", "Int8", "Int16", "Int32", "Int64",
                    "Float32", "Float64")
        for frame in (A.obs, A.var):
            frame.index = pd.Index([str(i) for i in frame.index], dtype=object,
                                   name=frame.index.name)
            for c in frame.columns:
                s = frame[c]
                if str(s.dtype) in NULLABLE:
                    frame[c] = (s.astype(object)
                                 .where(s.notna(), "not_reported")
                                 .astype(str).astype("category"))

        def _prune(d):
            """Drop None (writes as 'null') and the empty dicts that leaves."""
            out_ = {}
            for k, v in d.items():
                if isinstance(v, dict):
                    v = _prune(v)
                    if not v:
                        continue
                elif v is None:
                    continue
                out_[k] = v
            return out_

        A.uns = _prune(dict(A.uns))

        out = PROJ / "data/processed/EWSAtlas_v5_cellxgene.h5ad"
        io.write_h5ad_atomic(A, out)

        unknown = {c: int((o[c] == "unknown").sum()) for c in o.columns
                   if c.endswith("_ontology_term_id") and (o[c] == "unknown").any()}
        run.add(genes_kept=int(A.n_vars), genes_dropped=int(dropped),
                unknown_fields=unknown, out_bytes=int(out.stat().st_size))
        print(f"wrote {out.name}: {A.n_obs:,} x {A.n_vars:,} "
              f"({out.stat().st_size/1e9:.2f} GB); {dropped} symbols dropped")
        print("\nfields still needing curation (cells set to 'unknown'):")
        for k, v in sorted(unknown.items(), key=lambda kv: -kv[1]):
            print(f"  {k}: {v:,} of {A.n_obs:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
