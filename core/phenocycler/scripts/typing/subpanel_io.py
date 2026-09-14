"""Inputs the subpanel classification reads, and the scale it decides on.

Every number here was reproduced against the classification already on disk before this
module was written, so it is a reconstruction rather than a reimplementation:

    celltype/uniform/donor_id=<d>/data.parquet      UniFORM-aligned per-cell values
      -> apply_scale(floor, span)                   broadtype/scratch/cosine_scale.parquet
      -> reference level = pooled p99 of max(subpanel) over cells in NO classified object
                           and outside every islet

Checked 2026-09-12 against ``subpanel/subpanel_log.json``, all four to four decimals:

    hormone 0.1010 · Nerve 0.2272 · Lymphatic 0.1776 · Blood vessel 0.3073

The scale is cohort-wide by the celltype lane's own design (``run_step3.build_scale``):
one floor per marker at the aligned background mode, one top at the cohort p99.9.  The
floor is the half that UniFORM did not flatten for twelve markers -- see the plan's
"Raised, not resolved here" -- and it is read, not recomputed, so this lane cannot drift
from the values the classification was made on.
"""
from __future__ import annotations

import pathlib

import numpy as np
import pandas as pd

REPO = pathlib.Path(__file__).resolve().parents[4]
RUN_ID = "30f6689e4d27d295f890"

RUN = REPO / "core/phenocycler/data/runs" / RUN_ID
ALIGNED = REPO / "data/development/celltype/uniform"
SCALE_CSV = REPO / "data/development/broadtype/scratch/cosine_scale.parquet"
OUT = REPO / "data/development/typing"
GEOJSON = pathlib.Path.home() / "IO60panc2nd/results/Geojson"
ISLETS = REPO / "data/Islets"

#: context -> the markers that answer that context's question.  Vessel is Caveolin/CD34 and
#: NOT CD31: the objects were drawn by `Cavdil1200` / `CD34dil1430`, so CD31 was never what
#: made them.  Exocrine is three markers because EpCAM (0.9-3.1x) and Keratin-5 (1.0-1.3x)
#: were measured not to separate the bulk from any removed class.
SUBPANEL: dict[str, tuple[str, ...]] = {
    "islet": ("INS", "GCG", "SST", "CD99"),
    "Nerve": ("B3TUBB", "CD56"),
    "Lymphatic": ("Podoplanin",),
    "Blood vessel": ("Caveolin", "CD34"),
    "exocrine": ("Ker8_18", "E_cadherin", "Pan_Cytokeratin"),
}

#: The three hormones, read together against the non-islet reference level.  CD99 is in the
#: islet subpanel but cannot hold a cell endocrine: when it wins the argmax the cell is not
#: endocrine and continues down the ordinary chain.
HORMONES: tuple[str, ...] = ("INS", "GCG", "SST")

#: The quantile of the reference population that sets every structure's level.
REFERENCE_QUANTILE = 99.0


def load_scale() -> pd.DataFrame:
    """The cohort floor/top/span per marker, as the classification decided on it."""
    return pd.read_parquet(SCALE_CSV).set_index("marker")


def apply_scale(values: np.ndarray, markers: list[str], scale: pd.DataFrame) -> np.ndarray:
    """(cells x markers) -> [0, 1], clipped at the measured floor and the cohort top.

    Identical to ``celltype/score.apply_scale``; duplicated rather than imported so this
    lane's inputs cannot change under it when that lane is edited.
    """
    lo = scale.loc[markers, "floor"].to_numpy(float)
    span = scale.loc[markers, "span"].to_numpy(float)
    span = np.where(np.isfinite(span) & (span > 0), span, 1.0)
    return np.clip((np.asarray(values, dtype=np.float32) - lo) / span, 0.0, 1.0)


def reference_level(pooled: np.ndarray, quantile: float = REFERENCE_QUANTILE) -> float:
    """The level a structure's own markers must clear, from the tissue it sits in.

    A non-finite level would silently classify nothing -- ``value >= nan`` is False for every
    cell -- which is exactly how eight NaN cells in 10.4 million disabled a whole tier and
    left every downstream number looking plausible.  It fails loudly instead.
    """
    level = float(np.nanpercentile(pooled, quantile))
    if not np.isfinite(level):
        raise ValueError(f"reference level is not finite ({level!r}) -- refusing to classify")
    return level
