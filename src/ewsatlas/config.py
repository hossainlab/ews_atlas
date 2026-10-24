"""Single source of truth for pipeline parameters.

Rationale: thresholds scattered across scripts are the usual reason a
single-cell pipeline cannot be re-run years later. Every entrypoint reads
``src/config.json`` through :func:`load` and never defines a default of its own;
a missing key is an error, not a silent fallback.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

#: project root = two levels above this file (src/ewsatlas/config.py)
ROOT = Path(__file__).resolve().parents[2]
CONFIG_PATH = ROOT / "src" / "config.json"


def load(path: Path | None = None) -> dict[str, Any]:
    """Return the parsed config. Raises if the file is absent."""
    p = Path(path) if path else CONFIG_PATH
    if not p.is_file():
        raise FileNotFoundError(f"config not found: {p}")
    with p.open(encoding="utf-8") as fh:
        return json.load(fh)


def get(cfg: dict[str, Any], dotted: str) -> Any:
    """Fetch ``cfg['a']['b']`` as ``get(cfg, 'a.b')``; KeyError names the path.

    Deliberately strict: an analysis that silently substitutes a default for a
    missing QC threshold is not reproducible.
    """
    node: Any = cfg
    for part in dotted.split("."):
        if not isinstance(node, dict) or part not in node:
            raise KeyError(f"missing config key: {dotted!r} (at {part!r})")
        node = node[part]
    return node


def path(cfg: dict[str, Any], key: str) -> Path:
    """Resolve a ``paths.*`` entry to an absolute path under the project root."""
    return ROOT / get(cfg, f"paths.{key}")
