"""Arm B-aligned ``som-kintsugi-aligned``: the SOM arm on the KINTSUGI lane's **arm-A block** -- the
aligned intensities themselves (features.py ``load_block(..., arm="A")``: no rank, StandardScaler
inside ``som_cluster``), restricted to the declared type markers.  This is the plan's literal
specification for arm B ("scripts/kintsugi/features.py block restricted to registry type markers").

Declared 2026-09-06, AFTER ``som-kintsugi`` (within-donor-rank features) had been screened on both
donors and awarded no lineage cluster on either (its clusters separate by overall brightness; every
label contrast sat below the random-marker null).  Per the protocol a configuration change after
results is a NEW arm name, so the rank configuration's row stands as recorded and this row is judged
on its own.  Everything except the feature transform is byte-identical to ``som_kintsugi``: the same
node grid, seeds, null-gated cluster annotation, node choice by seed stability, and false-award
measurement on a shuffled twin.
"""
from __future__ import annotations

import pathlib
import sys

import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
import arm_base  # noqa: E402
import pins  # noqa: E402
import som_kintsugi as base  # noqa: E402

ARM = "som-kintsugi-aligned"
STOCHASTIC = base.STOCHASTIC


class _AlignedFeatures:
    """``arm_base`` with ``within_donor_rank`` replaced by the identity transform on aligned values."""

    @staticmethod
    def within_donor_rank(X: pd.DataFrame) -> pd.DataFrame:
        return X.astype(float).copy()

    def __getattr__(self, name):
        return getattr(arm_base, name)


DEFAULT_PARAMS = dict(base.DEFAULT_PARAMS)
DEFAULT_PARAMS["features"] = ("aligned intensity as loaded (UniFORM Step-2 values; KINTSUGI features.py arm A), declared type "
                              "markers only, no rank; StandardScaler inside som_cluster weights markers; a cell with any "
                              "non-finite feature is Unavailable, never imputed")
DEFAULT_PARAMS["declared_after"] = ("som-kintsugi (rank features) screened 2026-09-06 with no lineage cluster awarded on either "
                                    "donor; this configuration is the plan's literal arm-B specification and a new arm name")


def run(ctx: arm_base.ArmContext) -> pd.DataFrame:
    base.ARM = ARM
    base.SCRATCH_ROOT = pins.OUT / "scratch" / ARM
    base.DEFAULT_PARAMS = DEFAULT_PARAMS
    base.arm_base = _AlignedFeatures()
    try:
        return base.run(ctx)
    finally:
        base.arm_base = arm_base


if __name__ == "__main__":
    raise SystemExit(arm_base.cli(ARM, run, env=str(base.SOM_PYTHON), default_params=DEFAULT_PARAMS))
