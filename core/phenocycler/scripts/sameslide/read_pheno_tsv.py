#!/usr/bin/env python3
"""
QuPath measurement TSV -> a contract-valid PhenoCycler cell table, for a same-slide pair.

Why this does not go through ``phenocycler/cells_parquet.py``: that path is built for the
serial-section cohort and three of its assumptions do not hold for donor 6533's 70-plex run.

  1. Its header readers use ``csv.reader`` with the default comma dialect
     (``cells_parquet.header_cols``, ``artifacts.QuPathImageManifest.create``). On a
     tab-delimited export that yields ONE column and exits with "missing required canonical
     columns". DuckDB's ``read_csv`` further down sniffs the delimiter and would cope; the
     Python header readers never get that far.
  2. It keeps only whole-cell and ``Nucleus:`` means. A same-slide pair is the one case where
     the protein compartment breakdown is worth carrying, because it describes the SAME cell
     the RNA came from -- so all four compartments are kept here.
  3. Everything downstream of it (``redsea``, ``expression``) validates against
     ``marker_registry.json``, which has 60 markers to this panel's 70 and would refuse.

This reader therefore does the narrow job -- geometry + per-compartment means, on the
contract -- and leaves REDSEA/expression to a later pass once the registry is extended.

**INS/GCG/SST are deliberately KEPT.** ``export_pheno._IDENTITY_ONLY_MARKERS`` strips them,
and on serial sections that is right: a hormone measured on a different cell cannot be that
cell's identity. On one physical section it is a per-cell measurement of the same cell, and
it is the whole reason this pair exists -- none of the three is on either Xenium panel.
That exclusion is a serial-section policy, not a general one, so it does not apply here.

Usage:
    python scripts/sameslide/read_pheno_tsv.py \
        --tsv  /blue/.../data/phenocycler/6533.tsv \
        --donor 6533 --roi panc --serial 0080870
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.csv as pacsv

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from phenocycler.integration.contract import (  # noqa: E402
    CONTRACT_PROFILE_GEOMETRY,
    as_feature,
    coerce_cell_table,
    write_cell_table,
)

#: QuPath's four measurement compartments. All are kept -- see the module docstring.
COMPARTMENTS = ("Nucleus", "Cell", "Membrane", "Cytoplasm")

#: Only the mean enters the feature space. Min/Max/Std.Dev. exist in the export but are
#: single-cell dispersion statistics, not intensities, and would quadruple the width.
DEFAULT_STAT = "Mean"

#: QuPath writes cells with no enclosing annotation under this parent. Detection is run
#: inside annotations (``scripts/groovy/0_cell_detection.groovy`` calls ``selectAnnotations``),
#: so these are strays rather than a tissue region, and they are dropped by default.
UNANNOTATED_PARENT_PREFIX = "Root object"


def sanitize(name: str) -> str:
    """Marker name -> column-safe token. Matches ``cells_parquet.sanitize`` exactly.

    Kept identical on purpose: ``redsea.SANITIZE`` replaces only space, ``/`` and ``-``, so it
    renders ``Histone H3 (p Ser28)`` as ``Histone_H3_(p_Ser28)`` where this gives
    ``Histone_H3_p_Ser28``. Two spellings of one marker will not join, so the mapping this
    reader used is written to a sidecar rather than left to be re-derived.
    """
    return re.sub(r"[^0-9A-Za-z]+", "_", name).strip("_")


def read_header(tsv: Path) -> list[str]:
    """Column names only, without reading the body."""
    opts = pacsv.ReadOptions(block_size=1 << 20)
    with pacsv.open_csv(tsv, read_options=opts,
                        parse_options=pacsv.ParseOptions(delimiter="\t")) as rdr:
        return list(rdr.schema.names)


def plan_columns(header: list[str], stat: str = DEFAULT_STAT) -> tuple[dict, dict, list[str]]:
    """Split the header into the metadata we keep, the marker means, and what to request.

    Returns ``(meta_map, feature_map, include)`` where ``meta_map`` is
    ``{source_column: output_column}``, ``feature_map`` is
    ``{source_column: feat__<compartment>__<marker>}``.
    """
    meta: dict[str, str] = {}

    def find(pattern: str, out: str, *, required: bool = True) -> None:
        hits = [c for c in header if re.fullmatch(pattern, c)]
        if not hits:
            if required:
                raise SystemExit(f"required column matching {pattern!r} not found in the TSV")
            return
        meta[hits[0]] = out

    # Matched by pattern, not by literal: the centroid headers carry U+00B5 MICRO SIGN and
    # the area header a literal caret-2, both easy to mistype and silently miss.
    find(r"Object ID", "cell_id")
    find(r"Centroid X .?m", "x_um")
    find(r"Centroid Y .?m", "y_um")
    find(r"Cell: Area .?m\^2", "area_um2")
    find(r"Image", "image", required=False)
    find(r"Parent", "parent_raw", required=False)
    find(r"Object type", "object_type", required=False)
    find(r"Signed distance to annotation with Islet .?m", "dist_islet_um", required=False)
    find(r"Signed distance to annotation with Tissue .?m", "dist_tissue_um", required=False)

    features: dict[str, str] = {}
    pat = re.compile(rf"^({'|'.join(COMPARTMENTS)}): (.+): {re.escape(stat)}$")
    for col in header:
        m = pat.match(col)
        if m:
            comp, marker = m.group(1), m.group(2)
            features[col] = as_feature(f"{comp}__{sanitize(marker)}")

    return meta, features, list(meta) + list(features)


def read_table(tsv: Path, meta: dict, features: dict, *, block_mb: int = 256,
               keep_unannotated: bool = False) -> tuple[pd.DataFrame, dict]:
    """Stream the TSV, keeping only the planned columns."""
    ro = pacsv.ReadOptions(block_size=block_mb << 20)
    po = pacsv.ParseOptions(delimiter="\t")
    co = pacsv.ConvertOptions(include_columns=list(meta) + list(features))

    frames, n_read = [], 0
    t0 = time.time()
    with pacsv.open_csv(tsv, read_options=ro, parse_options=po, convert_options=co) as rdr:
        for batch in rdr:
            n_read += batch.num_rows
            frames.append(batch.to_pandas())
    df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    df = df.rename(columns={**meta, **features})

    stats = {"n_rows_tsv": int(n_read), "read_seconds": round(time.time() - t0, 1)}

    if "object_type" in df.columns:
        kinds = df["object_type"].value_counts().to_dict()
        stats["object_type"] = {str(k): int(v) for k, v in kinds.items()}
        df = df[df["object_type"] == "Cell"]

    if "parent_raw" in df.columns:
        parent = df["parent_raw"].fillna("")
        unannotated = parent.str.startswith(UNANNOTATED_PARENT_PREFIX)
        stats["n_unannotated"] = int(unannotated.sum())
        if not keep_unannotated:
            df = df[~unannotated]
        # Islet provenance, as audit context only. structures.py re-derives islets by DBSCAN
        # on both modalities so the comparison is not confounded by method; this column must
        # never be promoted to a structure_id.
        num = parent.str.extract(r"Islet_(\d+)", expand=False)
        df = df.assign(qupath_islet_num=pd.to_numeric(num, errors="coerce"))

    for c in ("x_um", "y_um", "area_um2", "dist_islet_um", "dist_tissue_um"):
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    feat_cols = [c for c in df.columns if c.startswith("feat__")]
    df[feat_cols] = df[feat_cols].astype(np.float32)

    stats["n_cells"] = int(len(df))
    stats["n_features"] = len(feat_cols)
    return df.reset_index(drop=True), stats


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", required=True, type=Path)
    ap.add_argument("--donor", required=True)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--serial", required=True, help="Xenium slide serial this section pairs with")
    ap.add_argument("--tissue", default="pancreas")
    ap.add_argument("--out", type=Path, default=None,
                    help="output root (default: data/integration/sameslide)")
    ap.add_argument("--stat", default=DEFAULT_STAT)
    ap.add_argument("--block-mb", type=int, default=256)
    ap.add_argument("--keep-unannotated", action="store_true")
    ap.add_argument("--limit-preview", type=int, default=0,
                    help="plan columns and exit without reading the body")
    args = ap.parse_args(argv)

    header = read_header(args.tsv)
    meta, features, include = plan_columns(header, args.stat)
    print(f"[read_pheno] header {len(header)} cols -> keeping {len(include)} "
          f"({len(meta)} meta + {len(features)} features)")
    if args.limit_preview:
        for src, out in list(features.items())[: args.limit_preview]:
            print(f"    {src!r} -> {out}")
        return 0

    df, stats = read_table(args.tsv, meta, features, block_mb=args.block_mb,
                           keep_unannotated=args.keep_unannotated)
    print(f"[read_pheno] {stats['n_rows_tsv']:,} rows read in {stats['read_seconds']}s "
          f"-> {stats['n_cells']:,} cells, {stats['n_features']} features")

    df = coerce_cell_table(df, modality="phenocycler", donor_id=args.donor,
                           profile=CONTRACT_PROFILE_GEOMETRY, roi=args.roi,
                           tissue=args.tissue)

    # Section-scoped output. The shared contract dirs key on (donor_id, roi) only, and donor
    # 6533 has two pancreas sections -- 0041393 and 0080870 -- so writing there would have
    # one silently replace the other. This lane keeps its own serial-qualified tree instead
    # of re-keying a layout 20 other donors already occupy.
    root = args.out or (Path(__file__).resolve().parents[2] / "data" / "integration" / "sameslide")
    out_dir = root / f"donor_id={args.donor}" / f"roi={args.roi}" / f"serial={args.serial}"
    path = write_cell_table(df, out_dir / "cells_pheno.parquet")

    sidecar = {
        "donor_id": args.donor, "roi": args.roi, "xenium_serial": args.serial,
        "tsv": str(args.tsv), "stat": args.stat,
        "marker_columns": {v: k for k, v in features.items()},
        "meta_columns": meta,
        "sanitizer": "cells_parquet.sanitize (non-alphanumeric runs -> _)",
        "ins_gcg_sst_kept": True,
        **stats,
    }
    (out_dir / "read_pheno_tsv.json").write_text(json.dumps(sidecar, indent=2) + "\n")
    print(f"[read_pheno] wrote {path}")
    print(f"[read_pheno] sidecar {out_dir / 'read_pheno_tsv.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
