"""EWSAtlas — a harmonized single-cell reference for EWSR1::FLI1 Ewing sarcoma.

Every pipeline parameter lives in ``src/config.json``; nothing in this package
hard-codes a threshold. Import order mirrors the pipeline:
``config -> io -> qc -> genespace -> integrate -> annotate -> discovery``.
"""

__version__ = "0.5.0.dev0"
