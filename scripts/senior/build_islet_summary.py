#!/usr/bin/env python3
"""Build the small, donor-stratified islet size/composition layer (spec §1.5).

Writes islets.parquet (one core per donor/islet), composition.parquet (observed
phenotype counts), donors.parquet (size strata), and availability.parquet under
$ISLET_LAYERS_DIR/islet_summary. No cell hive is scanned. The upstream tables
already exclude detections without nuclei; this builder does not reclassify cells.

These are aggregate tables: their explicit (donor_id, islet_key) identifies an
islet, not a cell. No donor is parsed from combined_islet_id, and no cell join or
synthetic object_id is introduced. Cell-level joins belong to their UUID-keyed
producers. The size table's n_cells and composition's region_total are different
denominators and are retained separately.

typing/geometry contains endocrine CELL UNIONS, not the islet border. Its WKB is
in pixels. The producer records area_um2 = WKB.area * image_scale**2, so the
recorded pair recovers that outline's isotropic calibration without a constant
or an image-name lookup. Conversion to phenocycler_um happens here, never in R.
Circularity and solidity describe this union, not the enclosing islet.

The converted outline itself is NOT written. Until 2026-09-23 it went out as a
``wkb`` column: 147.3 of islets.parquet's 148.7 MB, and no card read it
(card_s1_size_composition.R reads the counts, the states and the strata). What
it measured is kept -- endocrine_circularity, endocrine_solidity and
outline_um_per_px -- and the outline stays in typing/geometry, where it came from.

    .venv/bin/python scripts/senior/build_islet_summary.py --layers-dir /tmp/layers
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
from shapely import from_wkb
from shapely.affinity import scale

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from scripts.senior.build_sections_layer import filter_eligible_donors, norm_donor  # noqa: E402

DEFAULT_SOURCE = REPO / "data/development/typing"
KEY = ["donor_id", "islet_key"]
STRATA = ("0 cells", "1 cell", "2–10 cells", ">10 cells")
GEOMETRY_COLUMNS = [
    *KEY,
    "combined_islet_id",
    "region",
    "region_um2",
    "n_cells",
    "islet_diam_um",
    "centroid_x_um",
    "centroid_y_um",
]
COMPOSITION_COLUMNS = [
    *KEY,
    "combined_islet_id",
    "phenotype",
    "n_cells",
    "region",
    "region_total",
    "prop",
]
OUTLINE_COLUMNS = [*KEY, "endocrine_um2_from_outlines", "n_endocrine_outlines", "wkb"]


def refuse_hormone_composition(path: Path, frame: pd.DataFrame | None = None) -> None:
    """Refuse the known argmax file by path and its schema even when renamed."""
    hormone_columns = {"Ins_any", "Ins_single", "Multi_Pos", "Triple_Neg"}
    if path.name == "islet_hormone_composition.parquet" or (
        frame is not None and hormone_columns.intersection(frame.columns)
    ):
        raise ValueError(
            f"Refusing degenerate hormone composition source {path}: typing/app's "
            "7,634 rows have Ins_any == Ins_single and Multi_Pos = Triple_Neg = 0. "
            "Argmax phenotype counts are not hormone positivity. Use "
            "islet_composition_long.parquet for descriptive phenotype composition. "
            "data/app_data/parquet/uns_composition.parquet contains genuine "
            "per-donor thresholded positivity, a different domain and grain; it "
            "must not be substituted for per-islet composition."
        )


def _read(path: Path, columns: list[str], *, composition: bool = False):
    if composition:
        refuse_hormone_composition(path)
    if not path.exists():
        return None, f"source_missing:{path}"
    # The geometry directory also holds GeoJSON. Only its parquet is this source.
    files = sorted(path.rglob("*.parquet")) if path.is_dir() else [path]
    if not files:
        return None, f"no_parquet_files:{path}"
    table = ds.dataset([str(p) for p in files], format="parquet", partitioning=None)
    if composition:
        refuse_hormone_composition(path, pd.DataFrame(columns=table.schema.names))
    missing = set(columns) - set(table.schema.names)
    if missing:
        raise ValueError(f"{path}: missing columns {sorted(missing)}")
    frame = table.to_table(columns=columns).to_pandas()
    if frame.empty:
        return None, f"source_has_no_rows:{path}"
    for column in KEY:
        if (
            frame[column].isna().any()
            or not frame[column].map(lambda v: isinstance(v, str) and bool(v.strip())).all()
        ):
            raise ValueError(f"{path}: {column} must contain nonempty strings")
    frame["donor_id"] = frame["donor_id"].map(norm_donor)
    eligible = filter_eligible_donors(frame["donor_id"].unique())
    frame = frame[frame["donor_id"].isin(eligible)].copy()
    if frame.empty:
        return None, f"no_eligible_donors:{path}"
    return frame, ""


def _unique(frame: pd.DataFrame, keys: list[str], source: str) -> None:
    n = int(frame.duplicated(keys, keep=False).sum())
    if n:
        raise ValueError(f"{source}: {n} duplicate rows for {keys}")


def _is_integer(values: pd.Series) -> pd.Series:
    return values.notna() & np.isfinite(values) & values.ge(0) & values.mod(1).eq(0)


def _stratum(n: int) -> str:
    return STRATA[0 if n == 0 else 1 if n == 1 else 2 if n <= 10 else 3]


def _outlines(frame: pd.DataFrame | None) -> pd.DataFrame | None:
    if frame is None:
        return None
    _unique(frame, KEY, "endocrine outlines")
    rows = []
    for row in frame.to_dict("records"):
        out = {key: row[key] for key in KEY}
        out.update(
            outline_state="unavailable",
            outline_reason="missing_endocrine_outline_measurement",
            endocrine_um2_from_outlines=None,
            n_endocrine_outlines=None,
            outline_um_per_px=None,
            endocrine_circularity=None,
            endocrine_solidity=None,
        )
        area = row["endocrine_um2_from_outlines"]
        n = row["n_endocrine_outlines"]
        if row["wkb"] is None or pd.isna(area) or pd.isna(n):
            rows.append(out)
            continue
        try:
            outline = from_wkb(row["wkb"])
        except (ValueError, TypeError) as exc:
            raise ValueError(f"Invalid WKB for {out}: {exc}") from exc
        if (
            outline.geom_type not in {"Polygon", "MultiPolygon"}
            or outline.is_empty
            or not outline.is_valid
            or not np.isfinite(area)
            or area <= 0
            or not np.isfinite(n)
            or n <= 0
            or n % 1
            or outline.area <= 0
        ):
            out["outline_reason"] = "invalid_endocrine_outline_measurement"
        else:
            calibration = float(np.sqrt(area / outline.area))
            outline_um = scale(outline, xfact=calibration, yfact=calibration, origin=(0, 0))
            out.update(
                outline_state="measured",
                outline_reason="",
                endocrine_um2_from_outlines=float(area),
                n_endocrine_outlines=int(n),
                outline_um_per_px=calibration,
                endocrine_circularity=float(4 * np.pi * area / outline_um.length**2),
                endocrine_solidity=float(area / outline_um.convex_hull.area),
            )
        rows.append(out)
    return pd.DataFrame(rows)


def build_islet_summary(
    source_root: Path = DEFAULT_SOURCE,
    *,
    composition_source: Path | None = None,
) -> dict[str, pd.DataFrame]:
    """Read only the three small sources; return tables ready for serialization."""
    source_root = Path(source_root)
    geometry_path = source_root / "app/islet_geometry.parquet"
    composition_path = (
        Path(composition_source)
        if composition_source
        else (source_root / "app/islet_composition_long.parquet")
    )
    # Refuse the dangerous input even when another input happens to be missing.
    refuse_hormone_composition(composition_path)
    geometry, geometry_reason = _read(geometry_path, GEOMETRY_COLUMNS)
    composition, composition_reason = _read(
        composition_path,
        COMPOSITION_COLUMNS,
        composition=True,
    )
    outlines, outline_reason = _read(source_root / "geometry", OUTLINE_COLUMNS)
    availability = []
    for capability, frame, reason, path in (
        ("islet_size", geometry, geometry_reason, geometry_path),
        ("phenotype_composition", composition, composition_reason, composition_path),
        ("endocrine_outline", outlines, outline_reason, source_root / "geometry"),
    ):
        availability.append(
            dict(
                capability=capability,
                state="measured" if frame is not None else "unavailable",
                reason=reason,
                n=len(frame) if frame is not None else None,
                source=str(path.resolve()),
                source_domain_root=str((source_root / "app").resolve()),
            )
        )
    tables = {"availability": pd.DataFrame(availability)}
    tables["availability"]["n"] = tables["availability"]["n"].astype("Int64")
    if geometry is None:
        return tables
    islets = geometry[geometry["region"].eq("core")].copy()
    _unique(islets, KEY, "core size")
    if islets.empty:
        tables["availability"].loc[0, ["state", "reason", "n"]] = (
            "unavailable",
            "no_core_rows_in_islet_geometry",
            pd.NA,
        )
        return tables
    count_valid = _is_integer(islets["n_cells"])
    islets["state"] = np.where(count_valid, "measured", "unavailable")
    islets["reason"] = np.where(count_valid, "", "missing_or_invalid_core_n_cells")
    islets["n_cells"] = islets["n_cells"].where(count_valid).astype("Int64")
    islets["size_stratum"] = islets["n_cells"].map(
        lambda n: "unavailable" if pd.isna(n) else _stratum(n)
    )
    # Coordinate and area availability is independent of a measured cell count.
    spatial = ["region_um2", "islet_diam_um", "centroid_x_um", "centroid_y_um"]
    spatial_valid = np.isfinite(islets[spatial]).all(axis=1)
    spatial_valid &= islets["region_um2"].gt(0) & islets["islet_diam_um"].gt(0)
    islets["geometry_state"] = np.where(spatial_valid, "measured", "unavailable")
    islets["geometry_reason"] = np.where(spatial_valid, "", "missing_or_invalid_core_geometry")
    islets.loc[~spatial_valid, spatial] = np.nan
    islets["frame"] = "phenocycler_um:" + islets["donor_id"] + ":panc"
    islets["composition_state"] = "unavailable"
    islets["composition_reason"] = composition_reason or "no_core_composition_exported"
    islets["region_total"] = pd.array([None] * len(islets), dtype="Int64")
    if composition is not None:
        core = composition[composition["region"].eq("core")].copy()
        _unique(core, KEY + ["phenotype"], "core composition")
        if not core.empty:
            valid = _is_integer(core["n_cells"]) & _is_integer(core["region_total"])
            valid &= core["region_total"].gt(0) & core["phenotype"].notna()
            valid &= core["phenotype"].astype(str).str.strip().ne("")
            valid &= np.isfinite(core["prop"]) & np.isclose(
                core["prop"],
                core["n_cells"] / core["region_total"],
                rtol=1e-9,
                atol=1e-12,
            )
            if not valid.all():
                raise ValueError(
                    f"core composition: {int((~valid).sum())} invalid count/proportion rows"
                )
            grouped = core.groupby(KEY, sort=False)
            counts = grouped["n_cells"].sum()
            totals = grouped["region_total"].first()
            if not grouped["region_total"].nunique().eq(1).all() or not counts.eq(totals).all():
                raise ValueError(
                    "core composition: phenotype counts do not sum to one recorded region_total"
                )
            core = core.rename(columns={"n_cells": "phenotype_n_cells"})
            core = core.merge(
                islets[KEY + ["n_cells", "size_stratum", "state", "reason"]].rename(
                    columns={"state": "n_cells_state", "reason": "n_cells_reason"}
                ),
                on=KEY,
                how="left",
                validate="many_to_one",
                indicator=True,
            )
            if core["_merge"].ne("both").any():
                raise ValueError("core composition references islets absent from core geometry")
            core = core.drop(columns="_merge")
            core["state"] = "measured"
            core["reason"] = ""
            core["composition_kind"] = "phenotype_counts"
            core["frame"] = "none"  # counts/proportions are not tissue coordinates
            tables["composition"] = core.sort_values(KEY + ["phenotype"]).reset_index(drop=True)
            islets = islets.drop(columns="region_total").merge(
                totals.rename("region_total").reset_index(),
                on=KEY,
                how="left",
                validate="one_to_one",
            )
            has_composition = islets["region_total"].notna()
            islets.loc[has_composition, ["composition_state", "composition_reason"]] = [
                "measured",
                "",
            ]
            islets["region_total"] = islets["region_total"].astype("Int64")
    outline_table = _outlines(outlines)
    if outline_table is not None:
        islets = islets.merge(outline_table, on=KEY, how="left", validate="one_to_one")
    else:
        for column in (
            "endocrine_um2_from_outlines",
            "n_endocrine_outlines",
            "outline_um_per_px",
            "endocrine_circularity",
            "endocrine_solidity",
        ):
            islets[column] = None
        islets["outline_state"] = "unavailable"
        islets["outline_reason"] = outline_reason
    islets["outline_state"] = islets["outline_state"].fillna("unavailable")
    islets["outline_reason"] = islets["outline_reason"].fillna("no_endocrine_outline_exported")
    islets["n_endocrine_outlines"] = islets["n_endocrine_outlines"].astype("Int64")
    islets["outline_kind"] = "endocrine_cell_union"
    islets["outline_calibration"] = "sqrt(endocrine_um2_from_outlines / source_wkb_area_px2)"
    islets["positivity_state"] = "unavailable"
    islets["positivity_reason"] = "phenotype_argmax_is_not_thresholded_hormone_positivity"
    tables["islets"] = islets.sort_values(KEY).reset_index(drop=True)
    donors = []
    for donor, group in islets.groupby("donor_id", sort=True):
        for stratum in STRATA + ("unavailable",):
            here = group[group["size_stratum"].eq(stratum)]
            if stratum == "unavailable" and here.empty:
                continue
            measured = stratum != "unavailable"
            donors.append(
                dict(
                    donor_id=donor,
                    size_stratum=stratum,
                    n_islets=len(here),
                    n_cells=int(here["n_cells"].sum()) if measured else None,
                    state="measured" if measured else "unavailable",
                    reason="" if measured else "missing_or_invalid_core_n_cells",
                    n_composition_measured=int(here["composition_state"].eq("measured").sum()),
                    n_outline_measured=int(here["outline_state"].eq("measured").sum()),
                    frame="none",
                )
            )
    tables["donors"] = pd.DataFrame(donors)
    tables["donors"]["n_cells"] = tables["donors"]["n_cells"].astype("Int64")
    # Availability counts the islets actually emitted, not all four source regions.
    for index, column in enumerate(("state", "composition_state", "outline_state")):
        n = int(islets[column].eq("measured").sum())
        if n:
            tables["availability"].loc[index, "n"] = n
        elif tables["availability"].loc[index, "state"] == "measured":
            tables["availability"].loc[index, ["state", "reason", "n"]] = (
                "unavailable",
                f"no_measured_core_islets:{column}",
                pd.NA,
            )
    return tables


def write_summary(tables: dict[str, pd.DataFrame], layers_dir: Path) -> Path:
    """Replace this builder's four artifacts; remove a previously measured stale table."""
    target = Path(layers_dir) / "islet_summary"
    target.mkdir(parents=True, exist_ok=True)
    for name in ("islets", "composition", "donors", "availability"):
        path = target / f"{name}.parquet"
        if name not in tables:
            path.unlink(missing_ok=True)
            continue
        table = pa.Table.from_pandas(tables[name], preserve_index=False)
        pq.write_table(table, path)
    return target


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--composition-source", type=Path)
    parser.add_argument(
        "--layers-dir",
        type=Path,
        default=Path(os.environ.get("ISLET_LAYERS_DIR", "/data/islet_layers")),
    )
    args = parser.parse_args(argv)
    tables = build_islet_summary(args.source_root, composition_source=args.composition_source)
    target = write_summary(tables, args.layers_dir)
    print(f"Wrote {target}")
    if "islets" in tables:
        islets = tables["islets"]
        print(
            f"{len(islets):,} islets; {islets.donor_id.nunique()} donors kept separate; "
            f"{islets.composition_state.eq('measured').sum():,} with core composition; "
            f"{islets.outline_state.eq('measured').sum():,} with endocrine outlines"
        )
        print(tables["donors"].to_string(index=False))
    else:
        print(tables["availability"].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
