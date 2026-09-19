#!/usr/bin/env python3
"""Build trajectory/data.parquet: rank islets separately inside each donor.

This consumes the measured single-cell summaries; it does not fit a new trajectory,
orient it by disease, or impute the unsampled cells. The upstream typing/build_cells.py
excludes objects without nuclei before trajectory_singlecell_features.py selects cells.
No cell-level joins, spatial coordinates, or phenotype-based tuning occur here.

The core median determines the within-donor rank (average ranks for ties). The paired
peri-minus-core difference retains the source scale. Batch eta squared is a descriptive,
cell-weighted stain_group effect, computed separately inside each donor_status on the
sampled-cell pseudotime, never an acceptance gate or a correction. It is not the historical
0.988 diagnostic of the deployed cohort trajectory.

WHEN THE SOURCE FIT IS WITHIN DONOR (2026-09-17), that batch diagnostic stops being a
diagnostic and is reported Unavailable with its reason instead of a number. The producer now
fits one diffusion map per donor and hands over the WITHIN-DONOR QUANTILE, so every donor is
rank-uniform on (0, 1) by construction; pooling donors inside a status then gives every stain
group the same mean and an eta squared of 0. A zero that arithmetic forces is not evidence
that batch was removed, and printing it beside the word "batch" would read as exactly that.
The cross-donor relation is measured where it can be measured, with a permutation null, by
core/phenocycler/scripts/typing/trajectory_crossdonor_relation.py.

Run with .venv/bin/python. --check builds and round-trips the artifact in memory without
creating files. Otherwise only $ISLET_LAYERS_DIR/trajectory/data.parquet is written.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
from scripts.senior.availability_states import AVAILABILITY_STATES  # noqa: E402
from scripts.senior.build_sections_layer import filter_eligible_donors, norm_donor

DEFAULT_SOURCE = REPO / "data/development/typing/trajectory"
SUMMARY_FILE = "singlecell_islet_summary.parquet"
CELLS_FILE = "singlecell_pseudotime.parquet"
METADATA_FILE = "singlecell_trajectory.json"
#: Derived, never restated -- see `scripts/senior/availability_states.py`. Kept as a SET here
#: because this module only ever asks membership of it; the ordered tuple is the authority and
#: its order is not this module's to carry.
STATES = set(AVAILABILITY_STATES)
REGION_FIELDS = [
    f"sc_{field}_{region}"
    for region in ("core", "peri")
    for field in ("pseudotime", "pseudotime_mean", "n")
]
DELTA = "sc_pseudotime_peri_minus_core"
SUMMARY_COLUMNS = ["case_id", "islet_key", *REGION_FIELDS, DELTA]
CELL_COLUMNS = ["donor_id", "object_id", "donor_status", "stain_group", "pseudotime"]

SCHEMA = pa.schema([
    ("donor_id", pa.string()), ("donor_status", pa.string()),
    ("islet_key", pa.string()),
    *[(c, pa.int64() if c.startswith("sc_n_") else pa.float64()) for c in REGION_FIELDS],
    (DELTA, pa.float64()), ("within_donor_rank", pa.float64()),
    ("n_islets_ordered", pa.int64()),
    *[(c, pa.string()) for c in (
        "state", "reason", "core_state", "core_reason", "peri_state", "peri_reason",
        "delta_state", "delta_reason",
    )],
    ("batch_eta_sq_within_status", pa.float64()),
    ("batch_eta_sq_state", pa.string()), ("batch_eta_sq_reason", pa.string()),
    ("batch_n_cells", pa.int64()), ("batch_n_groups", pa.int64()),
    ("batch_factor", pa.string()),
    ("n_ordered_cells", pa.int64()), ("n_population_cells", pa.int64()),
    ("sampling_state", pa.string()), ("sampling_reason", pa.string()),
    ("frame", pa.string()), ("source", pa.string()), ("domain_data_root", pa.string()),
])


def _donors(values: pd.Series) -> pd.Series:
    # Reject missing/floating identifiers before normalization; never invent a donor
    # from a file name or truncate a numeric identifier through int(). case_id is int64.
    text = values.astype("string").str.strip()
    if text.isna().any() or not text.str.fullmatch(r"[0-9]+").all():
        raise ValueError("donor IDs must be nonmissing integer strings")
    return text.map(norm_donor)


def _require_columns(frame: pd.DataFrame, columns: list[str], source: str) -> None:
    missing = sorted(set(columns) - set(frame.columns))
    if missing:
        raise ValueError(f"{source}: missing columns {missing}")


def _unavailable(reason: str, *, source: str, domain_root: str) -> dict:
    return {
        **dict.fromkeys(SCHEMA.names),
        **{c: "unavailable" for c in SCHEMA.names if c.endswith("_state") or c == "state"},
        **{c: reason for c in SCHEMA.names if c.endswith("_reason") or c == "reason"},
        "frame": "none", "batch_factor": "stain_group",
        "source": source, "domain_data_root": domain_root,
    }


def batch_eta_sq_within_status(cells: pd.DataFrame, *, fit_scope: str | None = None) -> dict[str, dict]:
    """SS_between(stain_group) / SS_total, independently within each status.

    A single observed batch or constant pseudotime cannot estimate this diagnostic.
    Nonfinite pseudotime is not measured; missing batch labels invalidate the audit
    rather than silently dropping that batch. Disease labels are audit context only.

    A within-donor fit hands over a per-donor quantile, which makes this quantity zero by
    arithmetic rather than by measurement, so it is refused rather than reported. See the
    module docstring.
    """
    within_donor = fit_scope == "within_donor"
    result = {}
    for status, group in cells.groupby("donor_status", observed=True, sort=True):
        if within_donor:
            result[str(status)] = {
                "batch_eta_sq_within_status": None, "batch_eta_sq_state": "unavailable",
                "batch_eta_sq_reason": "within_donor_fit_makes_pooled_batch_eta_sq_arithmetic"
                                       "(each_donor_is_rank_uniform_on_its_own_axis)",
                "batch_n_cells": int(len(group)), "batch_n_groups": None,
            }
            continue
        group = group.loc[np.isfinite(group["pseudotime"])].copy()
        diagnostic = {
            "batch_eta_sq_within_status": None, "batch_eta_sq_state": "unavailable",
            "batch_eta_sq_reason": "no_finite_sampled_pseudotime",
            "batch_n_cells": None, "batch_n_groups": None,
        }
        batches = group["stain_group"].astype("string").str.strip()
        if group.empty:
            pass
        elif batches.isna().any() or batches.isin(["", "unknown", "Unknown"]).any():
            diagnostic["batch_eta_sq_reason"] = "sampled_cells_missing_stain_group"
        elif batches.nunique() < 2:
            diagnostic["batch_eta_sq_reason"] = "fewer_than_two_stain_groups_within_status"
        else:
            group["stain_group"] = batches
            values = group["pseudotime"]
            grand = values.mean()
            total = float(((values - grand) ** 2).sum())
            if total == 0:
                diagnostic["batch_eta_sq_reason"] = "constant_pseudotime_within_status"
            else:
                between = sum(
                    len(batch) * (batch["pseudotime"].mean() - grand) ** 2
                    for _, batch in group.groupby("stain_group", observed=True)
                )
                diagnostic.update(
                    batch_eta_sq_within_status=float(between / total),
                    batch_eta_sq_state="measured", batch_eta_sq_reason="",
                    batch_n_cells=len(group), batch_n_groups=int(batches.nunique()),
                )
        result[str(status)] = diagnostic
    return result


def _sampling(metadata: dict | None, n_rows: int) -> dict:
    result = {
        "n_ordered_cells": n_rows, "n_population_cells": None,
        "sampling_state": "unavailable", "sampling_reason": "trajectory_metadata_not_found",
    }
    if metadata is None:
        return result
    ordered, population = metadata.get("n_ordered"), metadata.get("n_population")
    if (type(ordered) is not int or type(population) is not int
            or ordered != n_rows or not 0 < ordered <= population):
        raise ValueError("trajectory sampling metadata disagrees with the pseudotime table")
    result.update(n_population_cells=population, sampling_state="measured", sampling_reason="")
    return result


def _region(row: dict, region: str) -> tuple[str, str]:
    count = row[f"sc_n_{region}"]
    median, mean = row[f"sc_pseudotime_{region}"], row[f"sc_pseudotime_mean_{region}"]
    if not pd.isna(count) and (not np.isfinite(count) or count <= 0 or count % 1):
        raise ValueError(f"sc_n_{region} must be a positive integer or absent")
    if pd.isna(count):
        if np.isfinite(median) or np.isfinite(mean):
            raise ValueError(f"{region} pseudotime has no sampled-cell count")
        return "unavailable", f"no_sampled_{region}_cells"
    if np.isfinite(median) != np.isfinite(mean):
        raise ValueError(f"{region} median and mean have inconsistent availability")
    if np.isfinite(median):
        return "measured", ""
    # Respect the producer's reporting policy; do not introduce a new cell-count gate.
    return "unavailable", f"{region}_pseudotime_not_reported(n_sampled={int(count)})"


def trajectory_table(
    summary: pd.DataFrame,
    cells: pd.DataFrame,
    metadata: dict | None,
    *,
    source: str = SUMMARY_FILE,
    domain_root: str = "",
) -> pa.Table:
    """Pure, testable transformation. All islet ordering is grouped by donor_id."""
    _require_columns(summary, SUMMARY_COLUMNS, SUMMARY_FILE)
    _require_columns(cells, CELL_COLUMNS, CELLS_FILE)
    summary, cells = summary.copy(), cells.copy()
    sampling = _sampling(metadata, len(cells))
    summary["donor_id"] = _donors(summary["case_id"])
    cells["donor_id"] = _donors(cells["donor_id"])
    eligible = filter_eligible_donors(set(summary["donor_id"]) | set(cells["donor_id"]))
    summary = summary.loc[summary["donor_id"].isin(eligible)].copy()
    cells = cells.loc[cells["donor_id"].isin(eligible)].copy()
    if summary["islet_key"].isna().any() or summary["islet_key"].eq("").any():
        raise ValueError("summary contains an absent islet_key")
    if summary.duplicated(["donor_id", "islet_key"]).any():
        raise ValueError("duplicate (donor_id, islet_key) after donor normalization")
    if (cells["object_id"].isna().any() or cells["object_id"].eq("").any()
            or cells["object_id"].duplicated().any()):
        raise ValueError("sampled cells must have unique, nonmissing object_id UUIDs")

    statuses = {}
    for donor, group in cells.groupby("donor_id", observed=True):
        labels = group["donor_status"].dropna().astype(str).str.strip().unique()
        if len(labels) > 1:
            raise ValueError(f"donor {donor} has conflicting donor_status values")
        statuses[donor] = labels[0] if len(labels) and labels[0] else None
    diagnostics = batch_eta_sq_within_status(
        cells, fit_scope=(metadata or {}).get("fit_scope"))
    rows = []
    base = _unavailable("donor_status_not_recorded", source=source, domain_root=domain_root)
    for donor in sorted(eligible):
        donor_rows = summary.loc[summary["donor_id"] == donor]
        for item in donor_rows.to_dict("records") or [{"donor_id": donor}]:
            row = {**base, **sampling, "donor_id": donor, "donor_status": statuses.get(donor)}
            row.update(diagnostics.get(statuses.get(donor), {}))
            if "islet_key" not in item:
                for prefix in ("core_", "peri_", "delta_", ""):
                    row[f"{prefix}reason"] = "no_islet_summary_rows"
                rows.append(row)
                continue
            row.update({c: item[c] for c in ["islet_key", *REGION_FIELDS, DELTA]})
            for region in ("core", "peri"):
                state, reason = _region(item, region)
                row[f"{region}_state"], row[f"{region}_reason"] = state, reason
                if state != "measured":
                    row[f"sc_pseudotime_{region}"] = None
                    row[f"sc_pseudotime_mean_{region}"] = None
            row["state"], row["reason"] = row["core_state"], row["core_reason"]
            if donor not in statuses:
                row["state"], row["reason"] = "unavailable", "no_sampled_cells_for_donor"
            paired = row["core_state"] == row["peri_state"] == "measured"
            if paired and np.isfinite(item[DELTA]):
                if not np.isclose(item[DELTA], item["sc_pseudotime_peri"] - item["sc_pseudotime_core"],
                                  rtol=1e-10, atol=1e-12):
                    raise ValueError("peri-minus-core differs from the paired source medians")
                row["delta_state"], row["delta_reason"] = "measured", ""
            else:
                row[DELTA] = None
                row["delta_state"] = "unavailable"
                row["delta_reason"] = (
                    "peri_minus_core_not_reported" if paired else ";".join(
                        row[f"{r}_reason"] for r in ("core", "peri") if row[f"{r}_reason"]
                    )
                )
            rows.append(row)
    if not rows:
        rows = [{**_unavailable("no_eligible_trajectory_rows", source=source,
                               domain_root=domain_root), **sampling}]
    out = pd.DataFrame(rows)
    # Keeping donor_id in the grouping is the acceptance contract. No cohort rank exists.
    measured = out["state"].eq("measured")
    ranked = out.loc[measured].groupby("donor_id", observed=True)["sc_pseudotime_core"]
    out.loc[measured, "within_donor_rank"] = ranked.rank(method="average")
    counts = out.loc[measured].groupby("donor_id", observed=True).size()
    out["n_islets_ordered"] = out["donor_id"].map(counts)
    out = out.sort_values(["donor_id", "within_donor_rank", "islet_key"], na_position="last")
    return pa.Table.from_pandas(out, schema=SCHEMA, preserve_index=False)


def build_trajectory_layer(source_dir: Path = DEFAULT_SOURCE, *, domain_root: Path | None = None) -> pa.Table:
    source_dir = Path(source_dir)
    domain = str((domain_root or source_dir.parent / "app").resolve())
    paths = [source_dir / SUMMARY_FILE, source_dir / CELLS_FILE]
    absent = [p.name for p in paths if not p.is_file()]
    if absent:
        return pa.Table.from_pylist([
            _unavailable("source_not_found:" + ",".join(absent),
                         source=str(source_dir.resolve()), domain_root=domain)
        ], schema=SCHEMA)
    summary, cells = [
        ds.dataset(path, format="parquet", partitioning=None).to_table(columns=columns).to_pandas()
        for path, columns in zip(paths, (SUMMARY_COLUMNS, CELL_COLUMNS), strict=True)
    ]
    metadata_path = source_dir / METADATA_FILE
    metadata = json.loads(metadata_path.read_text()) if metadata_path.is_file() else None
    return trajectory_table(summary, cells, metadata, source=str(paths[0].resolve()), domain_root=domain)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--domain-data-root", type=Path)
    parser.add_argument("--layers-dir", type=Path,
                        default=Path(os.environ.get("ISLET_LAYERS_DIR", "/data/islet_layers")))
    parser.add_argument("--check", action="store_true", help="round-trip parquet in memory; write no files")
    args = parser.parse_args(argv)
    table = build_trajectory_layer(args.source_dir, domain_root=args.domain_data_root)
    if args.check:
        buffer = pa.BufferOutputStream()
        pq.write_table(table, buffer, compression="zstd")
        restored = pq.ParquetFile(pa.BufferReader(buffer.getvalue())).read()
        if not table.equals(restored):
            raise AssertionError("trajectory parquet did not round-trip")
    else:
        target = args.layers_dir / "trajectory/data.parquet"
        target.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(table, target, compression="zstd")
        print(f"Wrote {target}")
    frame = table.to_pandas()
    print(f"{len(frame):,} rows; {frame.donor_id.nunique()} donors; "
          f"{frame.state.eq('measured').sum():,} ordered islets; "
          f"{frame.delta_state.eq('measured').sum():,} paired differences")
    print(frame[["donor_status", "batch_eta_sq_within_status", "batch_eta_sq_state",
                 "batch_eta_sq_reason", "batch_n_cells", "batch_n_groups"]]
          .drop_duplicates().to_string(index=False))
    print(frame[["n_ordered_cells", "n_population_cells", "sampling_state", "sampling_reason"]]
          .drop_duplicates().to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
