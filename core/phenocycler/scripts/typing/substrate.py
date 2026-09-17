#!/usr/bin/env python3
"""Build the one per-cell table every typing step reads.

One row per QuPath object, keyed by ``object_id`` (the QuPath detection UUID), carrying for each
usable marker its value in that marker's own compartment (see :mod:`compartments`) and its
within-donor rank.

Three rules the table enforces, because getting any of them wrong has silently broken previous
attempts:

* **A missing measurement stays missing.**  Cells that fail geometry QC keep their row and carry
  NaN plus the QC reason; they are never filled with zero and never dropped, so any denominator
  can be stated over the population it was actually computed on.
* **Ranks are within donor.**  Donors differ in distribution *shape*, not just level, so no
  cohort-wide cutoff is valid.  Every downstream comparison is between cells of one donor.
* **The join is the UUID.**  Never positional: the run's ``redsea`` and ``geometry_qc``
  partitions are not written in ``expression`` row order.

Written to ``data/development/typing/substrate/donor_id=<donor>/data.parquet``.
"""

from __future__ import annotations

import argparse
import glob
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import compartments as C  # noqa: E402

REPO = Path(__file__).resolve().parents[4]
RUN = REPO / "core/phenocycler/data/runs/30f6689e4d27d295f890"
OUT = REPO / "data/development/typing/substrate"
CONTEXT = ("cell_region", "islet_num", "X_centroid", "Y_centroid")


def _one(pattern: str) -> str:
    hits = sorted(glob.glob(pattern))
    if not hits:
        raise FileNotFoundError(pattern)
    return hits[0]


def registry_markers() -> list[str]:
    reg = json.loads((REPO / "core/phenocycler/phenocycler/marker_registry.json").read_text())
    markers = reg["markers"] if isinstance(reg, dict) and "markers" in reg else reg
    return [m["name"] for m in markers]


def within_donor_rank(values: np.ndarray) -> np.ndarray:
    """Rank in [0, 1] over the finite entries; NaN stays NaN.

    Ties get their average rank, so a marker with a large clamped-to-zero mass does not order
    those cells arbitrarily among themselves.
    """
    out = np.full(values.shape, np.nan, dtype=np.float32)
    ok = np.isfinite(values)
    n = int(ok.sum())
    if n == 0:
        return out
    out[ok] = (pd.Series(values[ok]).rank(method="average").to_numpy() - 0.5) / n
    return out


def build(donor: str) -> Path:
    markers = C.typing_markers(registry_markers())

    geom = pd.read_parquet(
        _one(f"{RUN}/geometry_qc/donor_id={donor}/*.parquet"),
        columns=["object_id", "analysis_eligible", "analysis_reason"],
    )
    # Which compartments a cell actually HAS.  An object with no nucleus is not typed at all; an
    # object whose cell polygon equals its nucleus has no cytoplasm and no membrane, so a surface
    # or cytoplasmic marker measured on it is measured in a compartment that does not exist there
    # and only its nuclear markers may be read.  Measured on 6442: 225,345 objects have no
    # nucleus, and 127,803 of the analysis-eligible ones have a nucleus but no cytoplasm.
    areas = pd.read_parquet(
        _one(f"{RUN}/redsea/donor_id={donor}/*.parquet"),
        columns=["object_id", "nucleus_area_px", "cytoplasm_area_px"],
    )
    areas["object_id"] = areas["object_id"].astype(str)
    expr = pd.read_parquet(
        _one(f"{RUN}/expression/donor_id={donor}/*.parquet"),
        columns=["object_id", *CONTEXT],
    )
    red_cols = ["object_id"] + [C.column_for(m) for m in markers]
    red = pd.read_parquet(_one(f"{RUN}/redsea/donor_id={donor}/*.parquet"), columns=red_cols)

    for frame in (geom, expr, red):
        frame["object_id"] = frame["object_id"].astype(str)

    df = expr.merge(geom, on="object_id", how="left", validate="one_to_one")
    df = df.merge(red, on="object_id", how="left", validate="one_to_one")
    df = df.merge(areas, on="object_id", how="left", validate="one_to_one")
    if len(df) != len(expr):
        raise AssertionError(f"join changed row count: {len(expr)} -> {len(df)}")

    eligible = df["analysis_eligible"].fillna(False).to_numpy(bool)
    has_nucleus = df["nucleus_area_px"].fillna(0).to_numpy() > 0
    has_cytoplasm = df["cytoplasm_area_px"].fillna(0).to_numpy() > 0
    out = pd.DataFrame({
        "donor_id": donor,
        "object_id": df["object_id"].to_numpy(),
        "eligible": eligible,
        "has_nucleus": has_nucleus,
        "has_cytoplasm": has_cytoplasm,
        "unavailable_reason": np.where(eligible, "", df["analysis_reason"].fillna("no_qc_record")),
    })
    for col in CONTEXT:
        out[col] = df[col].to_numpy()

    coverage = {}
    for m in markers:
        v = df[C.column_for(m)].to_numpy(np.float64)
        v = np.where(eligible, v, np.nan)          # ineligible cells hold no measurement
        out[f"v__{m}"] = v.astype(np.float32)
        out[f"r__{m}"] = within_donor_rank(v)
        coverage[m] = {
            "compartment": C.compartment_for(m),
            "finite_of_eligible": float(np.isfinite(v[eligible]).mean()) if eligible.any() else 0.0,
            "zero_of_finite": float((v[np.isfinite(v)] == 0).mean()) if np.isfinite(v).any() else 0.0,
        }

    dest = OUT / f"donor_id={donor}"
    dest.mkdir(parents=True, exist_ok=True)
    path = dest / "data.parquet"
    out.to_parquet(path, index=False)
    (dest / "coverage.json").write_text(json.dumps({
        "donor_id": donor,
        "n_objects": int(len(out)),
        "n_eligible": int(eligible.sum()),
        "markers": coverage,
    }, indent=2))
    return path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--donors", nargs="+", required=True)
    args = ap.parse_args()
    for donor in args.donors:
        path = build(donor)
        size = path.stat().st_size / 1e9
        print(f"[{donor}] wrote {path} ({size:.2f} GB)")


if __name__ == "__main__":
    main()
