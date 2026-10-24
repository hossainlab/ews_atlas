"""Stage 03b - train the released scVI model and SAVE it.

    python src/pipeline/03b_train_scvi.py [--max-epochs N]

This stage did not exist. The v5 latent spaces in
``data/interim/embeddings/*.npy`` were produced in an interactive session and
only the coordinates were kept, so (a) the integration step was the one stage of
the pipeline that was not reproducible from the repository, and (b) no model
existed for anyone to map new data onto - which the Resource claim in Figure 6
requires.

The configuration is read from ``integration.scvi`` in src/config.json, i.e. the
same settings the released embedding was built with. The run writes the model to
``data/interim/scvi_model_v5/`` and the latent to
``embeddings/latent_scvi_sample_retrained.npy`` under a NEW name, then compares
it against the released ``latent_scvi_sample.npy`` rather than overwriting it.
Agreement is reported three ways because an exact match is not expected across
library versions: mean cosine of per-cell vectors after orthogonal Procrustes
alignment, correlation of pairwise distances on a 5,000-cell subsample, and
k-nearest-neighbour overlap. Whether the saved model may be shipped alongside the
released coordinates is a judgement on those numbers, and the numbers are printed
rather than asserted.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import anndata as ad
import scvi
from scipy.linalg import orthogonal_procrustes
from scipy.stats import spearmanr
from sklearn.neighbors import NearestNeighbors

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from ewsatlas import config, io, provenance  # noqa: E402


def compare(new: np.ndarray, old: np.ndarray, seed: int = 0, k: int = 15) -> dict:
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(new), size=min(5000, len(new)), replace=False)
    R, _ = orthogonal_procrustes(new[idx], old[idx])
    aligned = new[idx] @ R
    num = (aligned * old[idx]).sum(1)
    den = np.linalg.norm(aligned, axis=1) * np.linalg.norm(old[idx], axis=1)
    cos = float(np.mean(np.where(den > 0, num / np.where(den > 0, den, 1), 0.0)))

    def pdist_sub(X):
        d = np.linalg.norm(X[:1000, None, :] - X[None, :1000, :], axis=2)
        return d[np.triu_indices(1000, 1)]
    rho = float(spearmanr(pdist_sub(new[idx]), pdist_sub(old[idx])).statistic)

    nn_new = NearestNeighbors(n_neighbors=k + 1).fit(new[idx])
    nn_old = NearestNeighbors(n_neighbors=k + 1).fit(old[idx])
    a = nn_new.kneighbors(new[idx], return_distance=False)[:, 1:]
    b = nn_old.kneighbors(old[idx], return_distance=False)[:, 1:]
    overlap = float(np.mean([len(set(x) & set(y)) / k for x, y in zip(a, b)]))
    return {"procrustes_cosine": cos, "distance_spearman": rho, f"knn{k}_overlap": overlap,
            "n_compared": int(len(idx))}


def main(max_epochs: int | None = None) -> int:
    cfg = config.load()
    PROJ = config.ROOT
    interim = config.path(cfg, "interim")
    src = interim / "atlas_hvg_for_scvi.h5ad"
    s = config.get(cfg, "integration.scvi")
    scvi.settings.seed = 0

    with provenance.start_run("03b_train_scvi", cfg, [src]) as run:
        T = ad.read_h5ad(src)
        scvi.model.SCVI.setup_anndata(
            T, layer="counts", batch_key=str(s["batch_key"]),
            categorical_covariate_keys=list(s["categorical_covariate_keys"]),
            continuous_covariate_keys=list(s["continuous_covariate_keys"]))
        model = scvi.model.SCVI(
            T, n_latent=int(s["n_latent"]), n_layers=int(s["n_layers"]),
            gene_likelihood=str(s["gene_likelihood"]),
            encode_covariates=bool(s["encode_covariates"]),
            use_layer_norm=str(s["use_layer_norm"]),
            use_batch_norm=str(s["use_batch_norm"]))
        model.train(max_epochs=int(max_epochs or s["max_epochs"]),
                    early_stopping=bool(s["early_stopping"]),
                    early_stopping_patience=int(s["early_stopping_patience"]),
                    batch_size=int(s["batch_size"]),
                    plan_kwargs={"lr": float(s["lr"])})

        model_dir = interim / "scvi_model_v5"
        model.save(str(model_dir), overwrite=True, save_anndata=False)
        Z = model.get_latent_representation()
        emb = interim / "embeddings"
        np.save(emb / "latent_scvi_sample_retrained.npy", Z)

        old = np.load(emb / "latent_scvi_sample.npy")
        stats = compare(Z, old) if old.shape == Z.shape else {"shape_mismatch": True}
        n_epochs = int(len(model.history["elbo_train"]))
        run.add(epochs_run=n_epochs, model_dir=str(model_dir), **stats)
        print(f"trained {n_epochs} epochs; model -> {model_dir}")
        for k, v in stats.items():
            print(f"  {k}: {v}")
    return 0


def attach_model_latent() -> int:
    """Store the retrained latent on the atlas as the model's own reference space.

    The retrained model does not reproduce the released embedding locally: global
    structure agrees (Spearman 0.904 on pairwise distances) but only 34.7% of each
    cell's 15 nearest neighbours are shared. The released ``X_scvi`` was trained
    interactively before this stage existed, so the difference cannot be
    reconciled, only handled.

    Replacing ``X_scvi`` would invalidate the UMAP, the Leiden communities and the
    annotation built on it. Instead the retrained coordinates are stored next to
    it as ``X_scvi_model_v5``, and the query-mapping utility looks labels up in
    THAT space. This is correct because cell-type labels are properties of cells,
    not of an embedding: a query cell mapped by the model must find its reference
    neighbours in the space the model itself produces.
    """
    import numpy as np
    import anndata as ad
    cfg = config.load()
    PROJ = config.ROOT
    atlas = PROJ / config.get(cfg, "paths.atlas_out")
    Z = np.load(PROJ / "data/interim/embeddings/latent_scvi_sample_retrained.npy")
    A = ad.read_h5ad(atlas)
    if Z.shape[0] != A.n_obs:
        raise SystemExit(f"latent has {Z.shape[0]} rows, atlas has {A.n_obs}")
    A.obsm["X_scvi_model_v5"] = Z
    A.uns["scvi_model_v5"] = {
        "model_dir": "data/interim/scvi_model_v5",
        "epochs_run": 377,
        "agreement_with_released_X_scvi": {
            "procrustes_cosine": 0.776, "distance_spearman": 0.904,
            "knn15_overlap": 0.347, "n_compared": 5000},
        "note": ("X_scvi is the released embedding that the UMAP, Leiden communities "
                 "and annotation were built on. X_scvi_model_v5 is the latent the "
                 "saved model produces; map queries into THAT space, not X_scvi."),
    }
    io.write_h5ad_atomic(A, atlas)
    print(f"attached X_scvi_model_v5 {Z.shape} to {atlas.name}")
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-epochs", type=int, default=None)
    ap.add_argument("--attach", action="store_true",
                    help="store the retrained latent on the atlas as X_scvi_model_v5")
    a = ap.parse_args()
    raise SystemExit(attach_model_latent() if a.attach else main(a.max_epochs))
