"""Stage 01 — ingest each deposit and apply uniform QC.

    python src/pipeline/01_ingest_qc.py [--studies Aynaud_2020 ...]

Reads every raw deposit named in ``config.studies``, harmonizes gene symbols,
applies the assay-aware QC policy, and writes ``data/interim/qc_<study>.h5ad``
plus ``results/tables/qc_attrition.csv``. ``data/raw`` is never written to.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ewsatlas import config, io, provenance, qc  # noqa: E402

#: deposit layout per study: (reader, kwargs relative to data/raw)
DEPOSITS = {
    "Aynaud_2020":    ("read_10x_tar",       {"tar_path": "geo/GSE130024/GSE130024_RAW.tar"}),
    "Mao_2024":       ("read_10x_tar",       {"tar_path": "geo/GSE250523/GSE250523_RAW.tar"}),
    "Cillo_2022":     ("read_nested_tar",    {"tar_path": "geo/GSE198896/GSE198896_RAW.tar",
                                              "keep_prefixes": ["GSM5959129_Ewing_1_PBMC",
                                                                "GSM5959130_Ewing_1_TIL",
                                                                "GSM5959139_Ewing_2_PBMC_TIL",
                                                                "GSM5959140_Ewing_3_PBMC_TIL",
                                                                "GSM5959141_Ewing_4_PBMC_TIL"]}),
    "Goodspeed_2025": ("read_dense_tsv",     {"matrix_path": "geo/GSE277083/GSE277083_counts_matrix_primary.txt.gz",
                                              "metadata_path": "geo/GSE277083/GSE277083_cell_metadata_primary.txt.gz"}),
    "Visser_2023":    ("read_celseq_plates", {"tar_path": "geo/GSE243347/GSE243347_RAW.tar",
                                              "annotation_path": "geo/GSE243347/GSE243347_data_annotation.txt.gz"}),
    "Waltner_2025":   ("read_scpca_zip",     {"zip_path": "scpca/SCPCP000015_ann-data_2026-10-17.zip"}),
}

#: Goodspeed ships primary and blood matrices separately; both are ingested.
EXTRA_PARTS = {
    "Goodspeed_2025": [{"matrix_path": "geo/GSE277083/GSE277083_counts_matrix_blood.txt.gz",
                        "metadata_path": "geo/GSE277083/GSE277083_cell_metadata_blood.txt.gz"}],
}

#: which absolute gene floor applies, by the assay the study used
ASSAY_OF = {
    "Aynaud_2020": "10x_scRNA", "Cillo_2022": "10x_scRNA", "Goodspeed_2025": "10x_scRNA",
    "Mao_2024": "10x_scRNA", "Visser_2023": "CELseq2_scRNA", "Waltner_2025": "10x_snRNA",
}


def main(argv: list[str] | None = None) -> int:
    cfg = config.load()
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--studies", nargs="*", default=cfg["studies"])
    args = ap.parse_args(argv)

    raw_root = config.path(cfg, "raw")
    interim = config.path(cfg, "interim")
    interim.mkdir(parents=True, exist_ok=True)
    provenance.assert_raw_readonly([interim])

    floors = config.get(cfg, "qc.gene_floor_by_assay")
    mito_caps = config.get(cfg, "qc.mito_cap_pct_by_assay")
    inputs = []
    for s in args.studies:
        _, kw = DEPOSITS[s]
        inputs += [str(raw_root / v) for k, v in kw.items() if k.endswith("path")]

    rows = []
    with provenance.start_run("01_ingest_qc", cfg, inputs) as run:
        for study in args.studies:
            reader_name, kwargs = DEPOSITS[study]
            reader = getattr(io, reader_name)
            parts = [{k: (raw_root / v if k.endswith("path") else v) for k, v in kwargs.items()}]
            for extra in EXTRA_PARTS.get(study, []):
                merged = dict(kwargs)
                merged.update(extra)
                parts.append({k: (raw_root / v if k.endswith("path") else v) for k, v in merged.items()})

            adatas = [reader(study=study, **p) for p in parts]
            if len(adatas) == 1:
                adata = adatas[0]
            else:
                import anndata as ad
                adata = ad.concat(adatas, join="outer")
                adata.obs_names_make_unique()
                adata.obs["study"] = study

            filtered, row = qc.apply_qc(
                adata,
                gene_floor=int(floors[ASSAY_OF[study]]),
                mito_cap=float(mito_caps[ASSAY_OF[study]]),
                n_mads=int(config.get(cfg, "qc.mad_n")),
                doublet_threshold=float(config.get(cfg, "qc.doublet_threshold")),
                library_gate=bool(config.get(cfg, "qc.library_gate_median_genes_below_floor")),
                doublet_method=str(config.get(cfg, "qc.doublet_call_method")),
                doublet_rate_per_1k=float(config.get(cfg, "qc.doublet_rate_per_1k_cells")),
            )
            io.write_h5ad_atomic(filtered, interim / f"qc_{study}.h5ad")
            rows.append(row)
            print(f"{study:16s} raw={row['raw']:7,d} kept={row['kept']:7,d} ({row['pct_kept']:.1f}%)"
                  f"  floor={row['gene_floor_applied']} mito_cap={row['mito_cap_applied']}")

        # Upsert by study: running one study must not erase the others' rows.
        # A --studies Visser_2023 re-run previously left a one-row table and
        # broke every downstream figure that reads it.
        out = config.ROOT / "results" / "tables" / "qc_attrition.csv"
        out.parent.mkdir(parents=True, exist_ok=True)
        table = pd.DataFrame(rows)
        if out.is_file():
            prior = pd.read_csv(out)
            prior = prior[~prior["study"].isin(table["study"])]
            table = pd.concat([prior, table], ignore_index=True)
        order = {s: i for i, s in enumerate(cfg["studies"])}
        table = table.sort_values("study", key=lambda c: c.map(order)).reset_index(drop=True)
        table.to_csv(out, index=False)
        run.add(attrition=rows, total_kept=int(table["kept"].sum()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
