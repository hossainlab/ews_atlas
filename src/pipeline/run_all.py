"""Run the shipped pipeline stages in order.

    python src/pipeline/run_all.py [--from 01] [--to 02]

Stops at the first failure and prints the manifest path, so a partial run is
diagnosable from logs/ rather than from scrollback.
"""

from __future__ import annotations

import argparse
import runpy
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SHIPPED = ["01_ingest_qc.py", "02_build_atlas.py"]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--from", dest="start", default="01")
    ap.add_argument("--to", dest="end", default="99")
    args = ap.parse_args()
    for script in SHIPPED:
        stage = script[:2]
        if not (args.start <= stage <= args.end):
            continue
        print(f"\n=== {script} ===", flush=True)
        sys.argv = [script]
        try:
            runpy.run_path(str(HERE / script), run_name="__main__")
        except SystemExit as e:
            if e.code not in (0, None):
                print(f"stage {script} failed with exit {e.code}; see logs/", file=sys.stderr)
                return int(e.code)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
