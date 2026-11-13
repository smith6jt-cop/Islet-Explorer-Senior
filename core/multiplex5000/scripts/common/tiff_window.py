"""Compatibility import for source-tree scripts using the windowed TIFF reader."""
from pathlib import Path
import sys

try:
    from multiplex5000.tiff_window import read_window, selftest_bit_exact
except ModuleNotFoundError as exc:
    if exc.name not in {"multiplex5000", "multiplex5000.tiff_window"}:
        raise
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from multiplex5000.tiff_window import read_window, selftest_bit_exact
