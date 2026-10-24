"""EWSAtlas integration stack — scVI/scANVI + Harmony baseline + scIB benchmark.

Purpose-built for the EWSAtlas six-study integration. Deliberately EXCLUDES the
bundled singlecell_gpu env's cell2location and scGPT:

  - cell2location 0.1.5 imports `deep_getattr` from `pyro.infer.autoguide.utils`,
    which current pyro removed -> the bundled env's own import check fails at
    build time. It is a spatial-deconvolution tool and this project has no
    spatial data, so it is dropped rather than pinned around.
  - scGPT is not part of the integration plan (zero-shot foundation-model
    embeddings underperform scVI on batch integration in published benchmarks);
    dropping it also removes the torchtext/--no-deps fragility.

Versions are pinned to MATCH THE LOCAL `ews`/`ews-scvi` environments so that
objects written locally read remotely and vice versa:
    scvi-tools 1.5.1, scanpy 1.12.4, anndata 0.13.4, harmonypy 0.2.0
torch is 2.5.1+cu124 here (GPU wheels) against whatever CPU build the local
environment resolved (recorded exactly in src/envs/requirements-lock-local.txt
and in every logs/run_*.json) — the compute backend differs by design; the
data-model packages do not. Consequence for reproducibility: CPU stages are
bit-reproducible under the config seed, GPU training reproduces to tolerance
only, because cuDNN selects kernels at runtime and atomic reductions are
order-dependent.

infercnvpy is NOT in this image. CNA inference runs locally on CPU, so that
adding it would not change the result but would invalidate the recorded
image_ref/spec_sha in src/config.json.

scib-metrics supplies the integration benchmark (iLISI/kBET/graph connectivity
vs NMI/ARI/isolated labels); it is JAX-backed, hence jax[cuda12].
"""

import modal

META = {
    "packages": ["scvi-tools", "scanpy", "anndata", "harmonypy", "scib-metrics", "torch"],
    "gpu_default": "A100",
    "egress_domains": [],
}


def build(
    *, secrets: dict[str, str] | None = None
) -> tuple["modal.Image", dict[str, "modal.Volume"], dict[str, str]]:
    secrets = secrets or {}
    img = (
        # anndata >=0.13 requires Python >=3.12; 3.12 also matches the local
        # `ews-scvi` environment that trained the current model.
        modal.Image.debian_slim(python_version="3.12")
        .env({"LANG": "C.UTF-8", "LC_ALL": "C.UTF-8", "PYTHONIOENCODING": "utf-8"})
        .apt_install("git")
        .pip_install("torch==2.5.1")
        .pip_install(
            "scvi-tools==1.5.1",
            "scanpy==1.12.4",
            "anndata==0.13.4",
            "harmonypy==0.2.0",
            "scib-metrics==0.5.7",
            "leidenalg~=0.10",
            "igraph~=0.11",
            "numpy~=2.0",
            "pyarrow",
        )
        .run_commands(
            "python -c 'import torch, scanpy, anndata, scvi, harmonypy, scib_metrics; "
            "print(\"ALL_IMPORTS_OK\", torch.__version__, scvi.__version__)'"
        )
    )
    return img, {}, {}
