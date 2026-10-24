"""Run manifests, deterministic seeding, and the raw-data contract.

Each entrypoint opens with :func:`start_run` and closes with :func:`finish_run`,
producing ``logs/run_<stage>_<timestamp>.json`` that records the seed, resolved
config, input checksums, package versions and host. That file is what lets a
later reader establish *which* code and inputs produced a given object.

On determinism: seeding fixes the CPU path exactly, but GPU training reproduces
to tolerance rather than bit-exactness — cuDNN picks kernels at runtime and
atomic reductions are order-dependent. The manifest records the device so a
mismatch is attributable instead of mysterious.
"""

from __future__ import annotations

import datetime as _dt
import hashlib
import json
import os
import platform
import socket
import sys
from pathlib import Path
from typing import Any

from .config import ROOT

_PKGS = (
    "anndata", "scanpy", "scvi-tools", "torch", "harmonypy", "scib-metrics",
    "infercnvpy", "numpy", "scipy", "pandas", "scikit-learn", "h5py",
    "statsmodels", "leidenalg", "igraph",
)


def sha256(path: str | Path, chunk: int = 1 << 22) -> str:
    """Streaming sha256 of a file — used to pin raw inputs by content."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def package_versions() -> dict[str, str | None]:
    import importlib.metadata as im

    out: dict[str, str | None] = {}
    for name in _PKGS:
        try:
            out[name] = im.version(name)
        except Exception:
            out[name] = None
    return out


def seed_everything(seed: int) -> None:
    """Seed numpy, torch and scvi-tools from the single config seed."""
    import random

    import numpy as np

    random.seed(seed)
    np.random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    try:
        import torch

        torch.manual_seed(seed)
        if torch.cuda.is_available():
            torch.cuda.manual_seed_all(seed)
    except ModuleNotFoundError:
        pass
    try:
        import scvi

        scvi.settings.seed = seed
    except ModuleNotFoundError:
        pass


def device_info() -> dict[str, Any]:
    info: dict[str, Any] = {"platform": platform.platform(), "python": sys.version.split()[0],
                            "hostname_hashed": hashlib.sha256(socket.gethostname().encode()).hexdigest()[:12]}
    try:
        import torch

        info["torch"] = torch.__version__
        info["cuda_available"] = bool(torch.cuda.is_available())
        info["device"] = torch.cuda.get_device_name(0) if torch.cuda.is_available() else "cpu"
    except ModuleNotFoundError:
        info["torch"] = None
    return info


def assert_raw_readonly(paths: list[str | Path]) -> None:
    """Fail loudly if a stage is about to treat data/raw as writable.

    The contract is: ``data/raw`` is immutable. Re-deriving an intermediate is
    always allowed; mutating a downloaded deposit is never allowed, because it
    breaks the checksum chain that makes the rest of the record meaningful.
    """
    raw = (ROOT / "data" / "raw").resolve()
    for p in paths:
        rp = Path(p).resolve()
        if raw == rp or raw in rp.parents:
            raise PermissionError(f"refusing to write inside data/raw: {rp}")


class Run:
    """Context for one pipeline stage; writes its manifest on exit."""

    def __init__(self, stage: str, cfg: dict[str, Any], inputs: list[str | Path]):
        self.stage = stage
        self.cfg = cfg
        self.inputs = [str(Path(i)) for i in inputs]
        self.started = _dt.datetime.now()
        self.record: dict[str, Any] = {}

    def __enter__(self) -> "Run":
        seed = int(self.cfg["seed"])
        seed_everything(seed)
        self.record = {
            "stage": self.stage,
            "started": self.started.isoformat(timespec="seconds"),
            "seed": seed,
            "config": self.cfg,
            "inputs": [
                {"path": p, "sha256": sha256(p) if Path(p).is_file() else None,
                 "bytes": Path(p).stat().st_size if Path(p).is_file() else None}
                for p in self.inputs
            ],
            "packages": package_versions(),
            "device": device_info(),
            "determinism_note": (
                "CPU stages are bit-reproducible under this seed. GPU training is "
                "reproducible to tolerance only: cuDNN kernel selection and atomic "
                "reductions are nondeterministic."
            ),
        }
        return self

    def add(self, **kw: Any) -> None:
        """Attach stage-specific numbers (cell counts, thresholds applied)."""
        self.record.setdefault("results", {}).update(kw)

    def __exit__(self, exc_type, exc, tb) -> None:
        self.record["finished"] = _dt.datetime.now().isoformat(timespec="seconds")
        self.record["wall_s"] = round((_dt.datetime.now() - self.started).total_seconds(), 1)
        self.record["status"] = "error" if exc_type else "ok"
        if exc_type:
            self.record["error"] = f"{exc_type.__name__}: {exc}"
        out = ROOT / "logs" / f"run_{self.stage}_{self.started:%Y%m%dT%H%M%S}.json"
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(self.record, indent=1, default=str), encoding="utf-8")


def start_run(stage: str, cfg: dict[str, Any], inputs: list[str | Path]) -> Run:
    return Run(stage, cfg, inputs)
