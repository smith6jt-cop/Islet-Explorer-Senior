"""Computation parity: on one CODEX image, the ported method reproduces the pinned commit's tables exactly.

It shows the code moved intact, nothing about biology. The pinned CODEX run (Islet-Explorer 578767e, the v19 typing
outputs) is on the workstation, not on HiPerGator, so this skips unless a run directory is named:

    MPLX5K_CODEX_PARITY_RUN=/data/panc_CODEX/.../uniform-20donor-20261027 [MPLX5K_CODEX_PARITY_IMAGE=ND_0112] \\
        pytest tests/current/test_positivity_parity.py

The run must hold typing/{redsea,context,coverage}, cells/, typing/empty_area_nulls.csv and the tables levels.py,
canary.py and identity.py wrote. The option each marker uses is the run's own choice (levels_selection.csv,
levels_differentiation_selection.csv): the choice is made across the cohort, so one image reproduces the per-image
levels given it, and canary and identity given those levels.
"""

from __future__ import annotations

import os
from itertools import chain
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from multiplex5000.positivity import canary as CN
from multiplex5000.positivity import identity as ID
from multiplex5000.positivity import levels as LV
from multiplex5000.positivity import registry as RG
from multiplex5000.positivity import sources

RUN = os.environ.get("MPLX5K_CODEX_PARITY_RUN")
pytestmark = pytest.mark.skipif(
    not RUN or not Path(RUN, "typing", "levels.csv").is_file(),
    reason="the pinned CODEX typing run (Islet-Explorer 578767e) is only on the workstation; "
           "set MPLX5K_CODEX_PARITY_RUN to its run directory")
PANCREAS = ("codex-uniform-32", "pancreas")


@pytest.fixture(scope="module")
def run():
    return Path(RUN)


@pytest.fixture(scope="module")
def image(run):
    return os.environ.get("MPLX5K_CODEX_PARITY_IMAGE") or str(
        pd.read_csv(run / "typing" / "levels.csv")["image"].iloc[0])


@pytest.fixture(scope="module")
def tables(run):
    t = run / "typing"
    return {"levels": pd.read_csv(t / "levels.csv"),
            "levels_differentiation": pd.read_csv(t / "levels_differentiation.csv"),
            "selection": pd.read_csv(t / "levels_selection.csv"),
            "differentiation_selection": pd.read_csv(t / "levels_differentiation_selection.csv"),
            "canary_levels": pd.read_csv(t / "canary_levels.csv")}


def _same(got: pd.DataFrame, pinned: pd.DataFrame, key: list[str]) -> None:
    got = got[list(pinned.columns)].sort_values(key).reset_index(drop=True)
    pinned = pinned.sort_values(key).reset_index(drop=True)
    pd.testing.assert_frame_equal(got, pinned, check_exact=True, check_dtype=False)


def test_levels_are_reproduced(run, image, tables):
    reg = RG.load(*PANCREAS)
    diff = tables["differentiation_selection"]
    chosen = {r.marker: (r.measurement, r.method) for r in chain(
        tables["selection"][~tables["selection"]["marker"].isin(diff["marker"])].itertuples(), diff.itertuples())}
    positivity = {r.marker: (r.measurement, r.method)
                  for r in tables["selection"][tables["selection"]["marker"].isin(diff["marker"])].itertuples()}
    source = sources.CodexRunSource(run)
    rows = pd.DataFrame(LV.image_pass(source, image, reg.spec, "final", chosen, LV.empty_nulls(source, reg),
                                      positivity))
    for reading, name in (("positivity", "levels"), ("differentiation", "levels_differentiation")):
        pinned = tables[name][tables[name]["image"] == image]
        _same(rows[rows["reading"] == reading], pinned, ["marker"])


def test_canary_levels_are_reproduced(run, image, tables):
    reg = RG.load(*PANCREAS)
    got = pd.DataFrame(CN.image_levels(sources.CodexRunSource(run), image, reg.spec, tables))
    _same(got, tables["canary_levels"][tables["canary_levels"]["image"] == image], ["marker"])


def test_identity_is_reproduced(run, image, tables):
    reg = RG.load(*PANCREAS)
    got, _ = ID.classify_image(sources.CodexRunSource(run), image, reg.spec, tables)
    pinned = pd.read_parquet(run / "typing" / "identity" / f"{image}.parquet")
    assert np.array_equal(got["object_id"].to_numpy(), pinned["object_id"].to_numpy())
    pd.testing.assert_frame_equal(got[list(pinned.columns)].reset_index(drop=True), pinned.reset_index(drop=True),
                                  check_exact=True, check_dtype=False)
