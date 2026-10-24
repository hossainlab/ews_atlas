"""Feature-space construction.

Rule (recorded in ``config.gene_space.rule``): the intersection of the per-study
*annotations*, restricted to genes detected in at least three cells in at least
one study.

Why not the stricter detected-intersection: requiring detection in *every*
study penalises the shallowest arm, so a gene real in five cohorts is dropped
because the sparsest one missed it. Why not the union: genes absent from a
study's reference become structural zeros that an integration model reads as
biology. The intersection-of-annotations keeps only genes every study could in
principle measure, and the detection filter removes those nobody saw.
"""

from __future__ import annotations

import anndata as ad
import numpy as np


def detected_genes(adata, min_cells: int = 3) -> set[str]:
    counts = (adata.X > 0).sum(axis=0)
    counts = np.asarray(counts).ravel()
    return set(np.asarray(adata.var_names)[counts >= min_cells])


def build(per_study: dict[str, "ad.AnnData"], min_cells: int = 3,
          n_studies: int = 1) -> tuple[list[str], dict]:
    """Return ``(genes, report)`` for the six-study feature space.

    ``min_cells`` and ``n_studies`` are ``config.gene_space.min_cells_detected``
    and ``config.gene_space.detected_in_n_studies``: a gene is kept when every
    study's annotation contains it and at least ``n_studies`` studies detected
    it in at least ``min_cells`` cells.
    """
    annotations = {k: set(v.var_names) for k, v in per_study.items()}
    detected = {k: detected_genes(v, min_cells) for k, v in per_study.items()}
    inter_annot = set.intersection(*annotations.values())
    hits: dict[str, int] = {}
    for s in detected.values():
        for g in s:
            hits[g] = hits.get(g, 0) + 1
    genes = sorted(inter_annot & {g for g, n in hits.items() if n >= n_studies})
    report = {
        "n_selected": len(genes),
        "n_annotation_intersection": len(inter_annot),
        "n_detected_intersection_strict": len(set.intersection(*detected.values())),
        "n_union": len(set.union(*annotations.values())),
        "per_study_annotation": {k: len(v) for k, v in annotations.items()},
        "per_study_detected": {k: len(v) for k, v in detected.items()},
        "min_cells_detected": min_cells,
        "detected_in_n_studies": n_studies,
    }
    return genes, report


#: Canonical markers that MUST survive feature selection; absence means the
#: symbol harmonization broke, not that the biology changed.
MARKER_PANEL = ("FLI1", "EWSR1", "NKX2-2", "NR0B1", "CD99", "PTPRC", "COL1A1",
                "PECAM1", "CD3E", "LYZ", "MKI67", "TOP2A")


def check_markers(genes) -> dict[str, bool]:
    present = set(genes)
    return {g: (g in present) for g in MARKER_PANEL}
