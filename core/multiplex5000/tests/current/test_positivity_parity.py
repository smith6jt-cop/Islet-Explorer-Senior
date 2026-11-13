"""Computation parity: on one CODEX image, the ported method reproduces the pinned commit's tables exactly.

It shows the code moved intact, nothing about biology. The run is the CODEX "current candidate generation" (source
docs/typing_handoff_20261030.md), on the workstation only, so this skips when it is absent (as on HiPerGator):

    MPLX5K_CODEX_PARITY_RUN   default /data/panc_CODEX/Islet-Explorer/data/development/uniform-20donor-20261029-v19-qc-work
    MPLX5K_CODEX_PARITY_IMAGE default the first image of its levels.csv

The run must hold typing/{redsea,context,coverage}, cells/, typing/empty_area_nulls.csv and the tables levels.py,
canary.py and identity.py wrote. The option each marker uses is the run's own choice (levels_selection.csv,
levels_differentiation_selection.csv): the choice is made across the cohort, so one image reproduces the per-image
levels given it, and canary and identity given those levels.

Which code made the run. The typing method at 578767e is the one 2075d31 committed (2026-10-30 17:35:28 -04:00: whole-
cell positivity, membrane only for differentiating; no later commit up to 578767e touches levels, canary, identity,
area_coverage, redsea_cells or the registry). A run typed before 2075d31 is expected to differ on the ruled markers and
must be regenerated with 578767e first; the handoff of 2026-10-30 lists that rerun as still to do. The CODEX typing
scripts record no commit in levels.csv, canary_levels.csv or identity/, so the run's age is read from what it does
record: typing/coverage/<image>.json holds the sha256 of the registry, area_coverage.py, redsea_cells.py and the
levels.csv it read. The test fails loudly (it does not skip) when that record is missing, when the registry or code
hashes are not 578767e's (2075d31 changed the registry), when coverage read another levels.csv than the one present, or
when the levels, canary or identity tables are older than 2075d31's commit time (a weaker check: file times).
"""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime
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

DEFAULT_RUN = "/data/panc_CODEX/Islet-Explorer/data/development/uniform-20donor-20261029-v19-qc-work"
RUN = os.environ.get("MPLX5K_CODEX_PARITY_RUN") or DEFAULT_RUN
pytestmark = pytest.mark.skipif(
    not Path(RUN, "typing", "levels.csv").is_file(),
    reason=f"the CODEX typing run {RUN} is not here (workstation only); set MPLX5K_CODEX_PARITY_RUN to override")
#: 2075d31, the last commit before 578767e to change the typing method
RULINGS_COMMIT, RULINGS_TIME = "2075d31", datetime.fromisoformat("2026-10-30T17:35:28-04:00")
#: what typing/coverage/<image>.json records for a run made with 578767e's files
PINNED_SHA256 = {
    "registry": "aa3dfd2c34902acb6b7af049811ccd1cde3c212b46c5907b7534f10926f4d657",
    "coverage_code": "0684db0e5475f253776f47de37b372de065fc2f309cf23ba0724ad274f5099d2",
    "pixel_processing_code": "af1ae8577a19c9189e5672a997e6e8c84ad919bd4abb80670c8607d4e187eeef",
}
STALE = (f"regenerate the run's typing with Islet-Explorer 578767e first: tables made before {RULINGS_COMMIT} "
         f"({RULINGS_TIME.isoformat()}) predate the whole-cell positivity rulings and are expected to differ")
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


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


@pytest.fixture(scope="module", autouse=True)
def provenance(run, image):
    """Fail loudly when the run was not made with 578767e's typing method (see the module docstring)."""
    record = run / "typing" / "coverage" / f"{image}.json"
    if not record.is_file():
        pytest.fail(f"{record} is missing: the run records no provenance for {image}, so its commit is unknown; "
                    + STALE)
    recorded = json.loads(record.read_text()).get("inputs", {}).get("sha256", {})
    wrong = {k: recorded.get(k) for k, v in PINNED_SHA256.items() if recorded.get(k) != v}
    if wrong:
        pytest.fail(f"{record}: coverage was made with other files than 578767e's ({sorted(wrong)}); " + STALE)
    if recorded.get("typing/levels.csv") != _sha256(run / "typing" / "levels.csv"):
        pytest.fail(f"{record}: coverage read another levels.csv than the one in the run; " + STALE)
    tables = [run / "typing" / n for n in ("levels.csv", "levels_differentiation.csv", "levels_selection.csv",
                                           "canary_levels.csv")] + [run / "typing" / "identity" / f"{image}.parquet"]
    old = [p.name for p in tables if datetime.fromtimestamp(p.stat().st_mtime).astimezone() < RULINGS_TIME]
    if old:
        pytest.fail(f"{old} are older than {RULINGS_COMMIT}; " + STALE)


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


def test_identity_outside_spatial_conflicts_is_reproduced(run, image, tables):
    reg = RG.load(*PANCREAS)
    got, _ = ID.classify_image(sources.CodexRunSource(run), image, reg.spec, tables)
    pinned = pd.read_parquet(run / "typing" / "identity" / f"{image}.parquet")
    assert np.array_equal(got["object_id"].to_numpy(), pinned["object_id"].to_numpy())
    # The original independent areas settled some equal-depth rows even where the new method cannot
    # validate a scale. Equal-depth conflicts and tied hormones may legitimately propagate new labels.
    touched = (got["spatial_conflict_evidence"].ne("") |
               pinned["tags"].str.contains("tie_brightness_coverage", na=False) |
               pinned["status"].eq("tied_highest_hormones"))
    pd.testing.assert_frame_equal(got.loc[~touched, list(pinned.columns)].reset_index(drop=True),
                                  pinned.loc[~touched].reset_index(drop=True),
                                  check_exact=True, check_dtype=False)
