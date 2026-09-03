"""
Export clusters / phenotypes from an ``.h5ad`` to Xenium Explorer.

Usage
-----
    python scripts/export_xenium_explorer_groups.py \
        --input data/processed/THYHDL065_annotated.h5ad \
        --output data/xenium_explorer \
        --sample THYHDL065 \
        --groups leiden_0.5 celltype phenotype \
        --cell-id-key cell_id

This produces:

    data/xenium_explorer/THYHDL065_cell_groups.csv
    data/xenium_explorer/THYHDL065/analysis.zarr.zip
    data/xenium_explorer/THYHDL065_color_palette.json

Load the CSV in Xenium Explorer via ``Cell`` -> ``Add cell categorization``,
or drop the ``analysis.zarr.zip`` next to the ``experiment.xenium`` file so
the groupings appear in the ``Clusters`` dropdown automatically.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import scanpy as sc  # noqa: E402

from utils.xenium_explorer_export import export_for_xenium_explorer  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", "-i", type=Path, required=True, help="Input .h5ad file")
    p.add_argument("--output", "-o", type=Path, required=True, help="Output directory")
    p.add_argument("--sample", "-s", required=True, help="Sample name used as file prefix")
    p.add_argument(
        "--groups", "-g", nargs="+", required=True,
        help="One or more .obs columns (e.g. leiden_0.5 celltype phenotype).",
    )
    p.add_argument(
        "--cell-id-key", default=None,
        help="Column in .obs holding the Xenium cell_id. Defaults to common "
             "candidates, otherwise the AnnData index is used.",
    )
    p.add_argument(
        "--palette", default="tab20",
        help="Matplotlib colormap name used for default category colors.",
    )
    p.add_argument("--no-csv", action="store_true", help="Skip CSV export.")
    p.add_argument("--no-zarr", action="store_true", help="Skip analysis.zarr.zip export.")
    return p.parse_args()


def main() -> int:
    args = parse_args()

    if not args.input.exists():
        raise FileNotFoundError(f"Input file not found: {args.input}")

    adata = sc.read_h5ad(args.input)

    missing = [g for g in args.groups if g not in adata.obs.columns]
    if missing:
        raise KeyError(
            f"Groups not found in adata.obs: {missing}. "
            f"Available columns: {list(adata.obs.columns)}"
        )

    artifacts = export_for_xenium_explorer(
        adata,
        group_keys=args.groups,
        output_dir=args.output,
        sample_name=args.sample,
        cell_id_key=args.cell_id_key,
        palette=args.palette,
        write_csv=not args.no_csv,
        write_zarr=not args.no_zarr,
    )

    for kind, path in artifacts.items():
        print(f"  {kind:>8s}: {path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
