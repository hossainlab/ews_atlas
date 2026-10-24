"""Environment workarounds, isolated so they are visible and removable.

Nothing here changes a scientific result. Each entry records why it exists so a
future reader can delete it when the underlying constraint is gone.
"""

from __future__ import annotations

import contextlib


@contextlib.contextmanager
def sequential_joblib():
    """Force joblib onto its sequential backend for the enclosed block.

    Why: joblib's threading backend opens a named pipe for its worker pool, and
    this Windows sandbox denies ``CreateNamedPipe`` inside a *subprocess* while
    permitting it in a long-lived kernel. Any pynndescent-backed call -
    ``sc.pp.neighbors``, ``sc.pp.scrublet``, scib-metrics' kNN graphs -
    therefore succeeds interactively and fails with ``PermissionError: [WinError
    5]`` when the same code runs as ``python src/pipeline/...``. pynndescent
    passes an explicit ``n_jobs`` that overrides ``joblib.parallel_config``, so
    the backend has to be pinned rather than the worker count.

    Cost: single-threaded nearest-neighbour construction, i.e. slower. The
    result is identical - this is a scheduling constraint, not a numerical one.
    """
    import joblib

    with joblib.parallel_backend("sequential"):
        yield

@contextlib.contextmanager
def sequential_process_map():
    """Run infercnvpy's chromosome fan-out in-process.

    ``infercnvpy.tl._infercnv`` calls ``tqdm.contrib.concurrent.process_map``
    unconditionally - ``n_jobs=1`` does not bypass it - and the pool's named
    pipe is denied in this sandbox whenever the stage runs as a subprocess
    ("PermissionError: [WinError 5] Access is denied"). The same call succeeds
    in-kernel, which is what made the defect look intermittent. Swapping in a
    sequential map keeps the stage runnable from the command line, which is
    what reproducibility requires; the chromosome loop is not the bottleneck.
    """
    import infercnvpy.tl._infercnv as _m

    original = _m.process_map

    def _seq(fn, *iterables, **kwargs):
        for key in ("max_workers", "chunksize", "tqdm_class", "desc", "unit", "total"):
            kwargs.pop(key, None)
        return [fn(*args) for args in zip(*iterables)]

    _m.process_map = _seq
    try:
        yield
    finally:
        _m.process_map = original
