"""Map a new Ewing sarcoma query dataset onto EWSAtlas v5.

    python src/map_query_to_ewsatlas.py query.h5ad mapped.h5ad

scArches surgery on the released scVI model: reference weights are frozen and
only query-specific batch embeddings are trained, so reference cells keep their
coordinates. Labels are then transferred by a 15-nearest-reference-neighbour
vote in the shared 30-d latent space, and each cell carries the share of its
neighbours that agreed.

Query requirements: raw counts in .X or layers["counts"], and a `sample` column
in .obs (one level per library). Genes missing from the atlas HVG set are
zero-filled, which is what scArches expects of a query on a fixed feature space.

Expected accuracy, from the leave-one-study-out benchmark in
docs/12_atlas_utility.md: median 0.84 agreement with the atlas annotation
(0.85 excluding classes only one study contributes), rising to 0.72-0.90 among
cells where at least 60% of neighbours agree. Two known limits - a cell type that
exists in only one study of the reference transfers poorly or not at all, and the
Ewing_sarcoma/Ewing_proliferating split is a cycle-phase boundary that transfers
at a median of 32%.
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scvi
from sklearn.neighbors import KNeighborsClassifier

ROOT = Path(__file__).resolve().parents[1]
ATLAS = ROOT / "data" / "processed" / "EWSAtlas_v5.h5ad"
MODEL = ROOT / "data" / "interim" / "scvi_model_v5"
K = 15


def main(query_path: str, out_path: str, epochs: int = 100) -> None:
    for p, what in [(ATLAS, "released atlas"), (MODEL, "saved scVI model")]:
        if not p.exists():
            raise SystemExit(
                f"{what} not found at {p}. The model is built by "
                "`python src/pipeline/03b_train_scvi.py`.")
    ref = ad.read_h5ad(ATLAS)
    hv = list(pd.read_csv(MODEL / "var_names.csv", header=None)[0]) \
        if (MODEL / "var_names.csv").exists() else \
        list(ref.var_names[ref.var["highly_variable"].to_numpy()])

    q = ad.read_h5ad(query_path)
    if "counts" in q.layers:
        q.X = q.layers["counts"].copy()
    if "sample" not in q.obs:
        raise SystemExit("query .obs needs a 'sample' column (one level per library)")

    missing = [g for g in hv if g not in q.var_names]
    if missing:
        pad = ad.AnnData(X=np.zeros((q.n_obs, len(missing)), dtype=np.float32),
                         obs=q.obs.copy(), var=pd.DataFrame(index=missing))
        q = ad.concat([q, pad], axis=1, merge="first")
        print(f"zero-filled {len(missing)} atlas genes absent from the query")
    q = q[:, hv].copy()
    q.layers["counts"] = q.X.copy()
    for col, default in [("study", "query"), ("assay", "unknown"), ("pct_counts_mt", 0.0)]:
        if col not in q.obs:
            q.obs[col] = default

    scvi.model.SCVI.prepare_query_anndata(q, str(MODEL))
    vae_q = scvi.model.SCVI.load_query_data(q, str(MODEL))
    vae_q.train(max_epochs=epochs, plan_kwargs={"weight_decay": 0.0}, accelerator="cpu")
    q.obsm["X_scvi"] = vae_q.get_latent_representation()

    # Labels are looked up in X_scvi_model_v5, the latent THIS model produces -
    # not in the released X_scvi, which was trained interactively before the
    # training stage existed and agrees with it only globally (distance Spearman
    # 0.904, but 15-nearest-neighbour overlap 0.347). A query cell must find its
    # reference neighbours in the space it was itself projected into.
    space = "X_scvi_model_v5" if "X_scvi_model_v5" in ref.obsm else "X_scvi"
    if space == "X_scvi":
        print("WARNING: X_scvi_model_v5 is absent; falling back to the released "
              "X_scvi, which the saved model does not reproduce locally. Run "
              "`python src/pipeline/03b_train_scvi.py --attach` first.")
    knn = KNeighborsClassifier(n_neighbors=K).fit(
        ref.obsm[space], ref.obs["cell_type"].astype(str).to_numpy())
    q.obs["cell_type_predicted"] = knn.predict(q.obsm["X_scvi"])
    q.obs["cell_type_confidence"] = knn.predict_proba(q.obsm["X_scvi"]).max(axis=1)

    q.write_h5ad(out_path, compression="gzip")
    print(f"mapped {q.n_obs:,} query cells -> {out_path}")
    print(q.obs["cell_type_predicted"].value_counts().to_string())
    print(f"median transfer confidence: {q.obs['cell_type_confidence'].median():.2f}")
    low = float((q.obs["cell_type_confidence"] < 0.6).mean() * 100)
    print(f"{low:.1f}% of cells fall below the 0.6 confidence floor used in the benchmark")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2])
