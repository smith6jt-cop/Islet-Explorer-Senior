#!/usr/bin/env python3
"""
QuPath measurement export -> a contract-valid PhenoCycler cell table, for a same-slide pair.

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

**v3 export shape: two Parquet files, not one TSV.** The re-run ships per image
``results/v3/cells/<tag>.parquet`` -- all measurements, QuPath column names, all four
compartments -- plus ``results/v3/context/<tag>.parquet`` -- ``Object ID``, ``Parent``,
centroids, ``Signed distance to annotation with <Class> um``. A TSV of the cells table is
still kept for 6533 and 0112, so all three input shapes must produce the identical
contract-valid table:

  * ``--tsv``                          one tab-delimited export, everything in one file.
  * ``--parquet``                      one Parquet export, everything in one file (same
                                        shape as ``--tsv``, just columnar).
  * ``--parquet --context-parquet``    measurements and context split across two Parquet
                                        files, joined here on the QuPath ``Object ID``. A
                                        column that could come from either file (``Area``
                                        is a measurement; centroids and ``Parent`` are
                                        context) is looked for in the context file FIRST,
                                        falling back to the cells file -- the context export
                                        does not always carry every meta column, and when it
                                        does, it is authoritative over a coordinate rounded
                                        differently in an older cells export.

A joined read fails loudly, before any row leaves this module, when the two files' Object ID
sets do not match exactly: a silent inner join would quietly drop the very cells the mismatch
is telling you about.

Usage:
    python scripts/sameslide/read_pheno_tsv.py \
        --tsv  /blue/.../data/phenocycler/6533.tsv \
        --donor 6533 --roi panc --serial 0080870

    python scripts/sameslide/read_pheno_tsv.py \
        --parquet results/v3/cells/6533.parquet \
        --context-parquet results/v3/context/6533.parquet \
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
import pyarrow.parquet as pq

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

#: (pattern, output column, required). Matched by pattern, not by literal: the centroid
#: headers carry U+00B5 MICRO SIGN and the area header a literal caret-2, both easy to
#: mistype and silently miss.
META_PATTERNS: tuple[tuple[str, str, bool], ...] = (
    (r"Object ID", "cell_id", True),
    (r"Centroid X .?m", "x_um", True),
    (r"Centroid Y .?m", "y_um", True),
    (r"Cell: Area .?m\^2", "area_um2", True),
    (r"Image", "image", False),
    (r"Parent", "parent_raw", False),
    (r"Object type", "object_type", False),
    (r"Signed distance to annotation with Islet .?m", "dist_islet_um", False),
    (r"Signed distance to annotation with Tissue .?m", "dist_tissue_um", False),
)


def sanitize(name: str) -> str:
    """Marker name -> column-safe token. Matches ``cells_parquet.sanitize`` exactly.

    Kept identical on purpose: ``redsea.SANITIZE`` replaces only space, ``/`` and ``-``, so it
    renders ``Histone H3 (p Ser28)`` as ``Histone_H3_(p_Ser28)`` where this gives
    ``Histone_H3_p_Ser28``. Two spellings of one marker will not join, so the mapping this
    reader used is written to a sidecar rather than left to be re-derived.
    """
    return re.sub(r"[^0-9A-Za-z]+", "_", name).strip("_")


def read_header(tsv: Path) -> list[str]:
    """Column names only, without reading the body -- TSV."""
    opts = pacsv.ReadOptions(block_size=1 << 20)
    with pacsv.open_csv(tsv, read_options=opts,
                        parse_options=pacsv.ParseOptions(delimiter="\t")) as rdr:
        return list(rdr.schema.names)


def read_header_parquet(path: Path) -> list[str]:
    """Column names only, without reading the body -- Parquet."""
    return list(pq.ParquetFile(path).schema_arrow.names)


def _match_meta(header: list[str], patterns) -> tuple[dict, list]:
    """Match as many ``patterns`` against ``header`` as possible.

    Returns ``(matched, unmatched)`` where ``matched`` is ``{source_column: out_column}`` and
    ``unmatched`` is the leftover pattern tuples, so a caller can try a second header for
    whatever the first one did not have.
    """
    matched: dict[str, str] = {}
    unmatched = []
    for pattern, out, required in patterns:
        hits = [c for c in header if re.fullmatch(pattern, c)]
        if hits:
            matched[hits[0]] = out
        else:
            unmatched.append((pattern, out, required))
    return matched, unmatched


def find_meta_columns(*headers: list[str]) -> dict[str, str]:
    """Meta columns for a single header, or (source_column -> out_column) across several.

    Headers are tried in order and the first one to carry a given field wins -- used for the
    TSV/single-Parquet case (one header) and, unchanged in behaviour, for the split
    cells+context case (context first, cells file as fallback for e.g. ``Area``, which is a
    measurement carried only in the cells export).
    """
    patterns = list(META_PATTERNS)
    matched: dict[str, str] = {}
    for header in headers:
        hit, patterns = _match_meta(header, patterns)
        matched.update(hit)
    missing_required = [out for _, out, required in patterns if required]
    if missing_required:
        raise SystemExit(
            f"required column(s) for {missing_required} not found in any of the "
            f"{len(headers)} header(s) given")
    return matched


def find_feature_columns(header: list[str], stat: str = DEFAULT_STAT) -> dict[str, str]:
    """``{source_column: feat__<compartment>__<marker>}`` for every compartment mean."""
    features: dict[str, str] = {}
    pat = re.compile(rf"^({'|'.join(COMPARTMENTS)}): (.+): {re.escape(stat)}$")
    for col in header:
        m = pat.match(col)
        if m:
            comp, marker = m.group(1), m.group(2)
            features[col] = as_feature(f"{comp}__{sanitize(marker)}")
    return features


def plan_columns(header: list[str], stat: str = DEFAULT_STAT) -> tuple[dict, dict, list[str]]:
    """Split a single header into the metadata we keep, the marker means, and what to request.

    Returns ``(meta_map, feature_map, include)`` where ``meta_map`` is
    ``{source_column: output_column}``, ``feature_map`` is
    ``{source_column: feat__<compartment>__<marker>}``.
    """
    meta = find_meta_columns(header)
    features = find_feature_columns(header, stat)
    return meta, features, list(meta) + list(features)


def _find_object_id_column(header: list[str], *, source: str) -> str:
    hits = [c for c in header if re.fullmatch(r"Object ID", c)]
    if not hits:
        raise SystemExit(f"required column matching 'Object ID' not found in {source}")
    return hits[0]


def _postprocess(df: pd.DataFrame, *, keep_unannotated: bool, stats: dict) -> tuple[pd.DataFrame, dict]:
    """Everything after the raw columns are read and renamed -- shared by every input shape."""
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

    stats = {"input_format": "tsv", "n_rows_tsv": int(n_read), "read_seconds": round(time.time() - t0, 1)}
    return _postprocess(df, keep_unannotated=keep_unannotated, stats=stats)


def read_table_parquet(path: Path, meta: dict, features: dict, *,
                       keep_unannotated: bool = False) -> tuple[pd.DataFrame, dict]:
    """One Parquet file carrying both metadata and marker means (the ``--parquet``-only case)."""
    t0 = time.time()
    columns = {**meta, **features}
    table = pq.read_table(path, columns=list(columns))
    df = table.to_pandas().rename(columns=columns)

    stats = {"input_format": "parquet", "n_rows_parquet": int(len(df)),
             "read_seconds": round(time.time() - t0, 1)}
    return _postprocess(df, keep_unannotated=keep_unannotated, stats=stats)


def read_table_parquet_joined(cells_path: Path, context_path: Path,
                              meta: dict, features: dict, *,
                              keep_unannotated: bool = False) -> tuple[pd.DataFrame, dict]:
    """Cells + context Parquet, joined on the QuPath ``Object ID``.

    ``meta`` maps each SOURCE column (from whichever file it was found in by
    :func:`find_meta_columns`) to its output name; this function figures out, per source
    column, which of the two files to read it from.
    """
    t0 = time.time()
    cells_header = read_header_parquet(cells_path)
    context_header = read_header_parquet(context_path)
    cells_oid = _find_object_id_column(cells_header, source=f"cells parquet {cells_path.name!r}")
    context_oid = _find_object_id_column(context_header, source=f"context parquet {context_path.name!r}")

    cells_header_set = set(cells_header)
    context_header_set = set(context_header)
    meta_from_cells = {c: out for c, out in meta.items() if c not in context_header_set and c in cells_header_set}
    meta_from_context = {c: out for c, out in meta.items() if c in context_header_set}

    cells_columns = dict(features)
    cells_columns.update(meta_from_cells)
    cells_columns.setdefault(cells_oid, "cell_id")
    cells_df = pq.read_table(cells_path, columns=list(cells_columns)).to_pandas()
    cells_df = cells_df.rename(columns=cells_columns)

    context_columns = dict(meta_from_context)
    context_columns.setdefault(context_oid, "cell_id")
    context_df = pq.read_table(context_path, columns=list(context_columns)).to_pandas()
    context_df = context_df.rename(columns=context_columns)

    cells_ids = set(cells_df["cell_id"].astype(str))
    context_ids = set(context_df["cell_id"].astype(str))
    if cells_ids != context_ids:
        only_cells = len(cells_ids - context_ids)
        only_context = len(context_ids - cells_ids)
        raise SystemExit(
            f"Object ID mismatch joining {cells_path.name!r} ({len(cells_ids):,} ids) against "
            f"{context_path.name!r} ({len(context_ids):,} ids): {only_cells:,} id(s) only in "
            f"the cells file, {only_context:,} only in the context file -- refusing a silent "
            f"inner join")

    df = cells_df.merge(context_df, on="cell_id", how="inner", validate="one_to_one")

    stats = {"input_format": "parquet+context", "n_rows_cells_parquet": int(len(cells_df)),
             "n_rows_context_parquet": int(len(context_df)),
             "read_seconds": round(time.time() - t0, 1)}
    return _postprocess(df, keep_unannotated=keep_unannotated, stats=stats)


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--tsv", type=Path, help="one tab-delimited QuPath measurement export")
    src.add_argument("--parquet", type=Path,
                     help="one Parquet QuPath measurement export, or (with --context-parquet) "
                          "just the per-cell measurements half of a split v3 export")
    ap.add_argument("--context-parquet", type=Path, default=None,
                    help="v3 context export (Object ID, Parent, centroids, signed distances); "
                         "only valid together with --parquet")
    ap.add_argument("--donor", required=True)
    ap.add_argument("--roi", default="panc")
    ap.add_argument("--serial", required=True, help="Xenium slide serial this section pairs with")
    ap.add_argument("--tissue", default="pancreas")
    ap.add_argument("--out", type=Path, default=None,
                    help="output root (default: data/integration/sameslide)")
    ap.add_argument("--stat", default=DEFAULT_STAT)
    ap.add_argument("--block-mb", type=int, default=256, help="TSV read block size (--tsv only)")
    ap.add_argument("--keep-unannotated", action="store_true")
    ap.add_argument("--limit-preview", type=int, default=0,
                    help="plan columns and exit without reading the body")
    args = ap.parse_args(argv)

    if args.context_parquet is not None and args.parquet is None:
        ap.error("--context-parquet requires --parquet")

    if args.tsv is not None:
        header = read_header(args.tsv)
        meta, features, include = plan_columns(header, args.stat)
        print(f"[read_pheno] header {len(header)} cols -> keeping {len(include)} "
              f"({len(meta)} meta + {len(features)} features)")
        if args.limit_preview:
            for src_col, out in list(features.items())[: args.limit_preview]:
                print(f"    {src_col!r} -> {out}")
            return 0
        df, stats = read_table(args.tsv, meta, features, block_mb=args.block_mb,
                               keep_unannotated=args.keep_unannotated)
        print(f"[read_pheno] {stats['n_rows_tsv']:,} rows read in {stats['read_seconds']}s "
              f"-> {stats['n_cells']:,} cells, {stats['n_features']} features")
    elif args.context_parquet is not None:
        cells_header = read_header_parquet(args.parquet)
        context_header = read_header_parquet(args.context_parquet)
        meta = find_meta_columns(context_header, cells_header)
        features = find_feature_columns(cells_header, args.stat)
        print(f"[read_pheno] cells parquet {len(cells_header)} cols, context parquet "
              f"{len(context_header)} cols -> keeping {len(meta) + len(features)} "
              f"({len(meta)} meta + {len(features)} features)")
        if args.limit_preview:
            for src_col, out in list(features.items())[: args.limit_preview]:
                print(f"    {src_col!r} -> {out}")
            return 0
        df, stats = read_table_parquet_joined(args.parquet, args.context_parquet, meta, features,
                                              keep_unannotated=args.keep_unannotated)
        print(f"[read_pheno] {stats['n_rows_cells_parquet']:,} cells rows x "
              f"{stats['n_rows_context_parquet']:,} context rows joined in "
              f"{stats['read_seconds']}s -> {stats['n_cells']:,} cells, {stats['n_features']} features")
    else:
        header = read_header_parquet(args.parquet)
        meta, features, include = plan_columns(header, args.stat)
        print(f"[read_pheno] header {len(header)} cols -> keeping {len(include)} "
              f"({len(meta)} meta + {len(features)} features)")
        if args.limit_preview:
            for src_col, out in list(features.items())[: args.limit_preview]:
                print(f"    {src_col!r} -> {out}")
            return 0
        df, stats = read_table_parquet(args.parquet, meta, features,
                                       keep_unannotated=args.keep_unannotated)
        print(f"[read_pheno] {stats['n_rows_parquet']:,} rows read in {stats['read_seconds']}s "
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
        "tsv": str(args.tsv) if args.tsv else None,
        "parquet": str(args.parquet) if args.parquet else None,
        "context_parquet": str(args.context_parquet) if args.context_parquet else None,
        "stat": args.stat,
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
