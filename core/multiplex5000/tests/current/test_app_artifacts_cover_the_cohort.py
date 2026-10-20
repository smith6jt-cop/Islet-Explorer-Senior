"""The app's per-islet artifacts must cover the cohort their inputs reach, and agree with it.

The app showed ten donors while the typing covered twenty, and every test passed the whole time.
Nothing was broken: each artifact was internally consistent, every path resolved and every column
existed. They had simply been built on an older, smaller donor list and never rebuilt, and no
check compared one artifact's coverage against another's.

So this compares them to each other and to `islet_geometry`, and it derives the expected donor set
from `cohort.resolve()` rather than naming it. A test that pins the cohort as a literal has to be
edited every time a donor joins, which is the same failure in a different file -- the number would
be updated to match whatever was built, and would then agree with a stale build forever.

Covering the cohort is necessary and not sufficient: an artifact can carry every donor and still
have been built before an input it reads, in which case it carries the previous generation's
numbers for all twenty of them.  So the coverage assertions below are joined by a staleness one,
measured against the graph `scripts/senior/build_layers.py` declares -- for the same artifacts,
found through the declaration rather than through a second list kept here.  The layers and the
column families get the same treatment in `test_layer_staleness.py`.

Marked `data`: these are built artifacts, not fixtures, so the suite skips when they are absent
rather than failing a synthetic CI run.

    pytest core/multiplex5000/tests/current/test_app_artifacts_cover_the_cohort.py
"""

from __future__ import annotations

import importlib.util
import os
import sys
import time
from pathlib import Path

import pandas as pd
import pytest

REPO = Path(__file__).resolve().parents[4]
TYPING = REPO / "data/development/typing"
APP = TYPING / "app"
DRIVER = REPO / "scripts/senior/build_layers.py"

sys.path.insert(0, str(REPO / "core/multiplex5000/scripts/typing"))

pytestmark = pytest.mark.data


#: artifact -> the column that carries the donor. Two spellings because the app's CSV readers
#: type-convert `case_id` per column and the parquet writers keep `donor_id`; both are read here
#: as strings so `115` and `0115` cannot compare unequal.
ARTIFACTS = {
    "islet_composition_long.parquet": "donor_id",
    "islet_hormone_composition.parquet": "donor_id",
    "islet_celltype_proportions.csv": "case_id",
    "islet_spatial_lookup.csv": "case_id",
    "islet_marker_means.parquet": "case_id",
    "neighborhood_metrics.csv": "case_id",
    "parquet/islets.parquet": "case_id",
}


def _read(name: str) -> pd.DataFrame:
    path = APP / name
    if not path.is_file():
        pytest.skip(f"{path} is not built")
    return pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)


def _donors(frame: pd.DataFrame, column: str) -> set[str]:
    return set(frame[column].astype(str))


@pytest.fixture(scope="module")
def expected() -> set[str]:
    """The donors the builders would resolve, straight from the module they all call."""
    if not (APP / "islet_geometry.parquet").is_file():
        pytest.skip("islet_geometry.parquet is not built")
    import cohort as COHORT

    return set(COHORT.resolve(None, include_marginal=True))


@pytest.fixture(scope="module")
def geometry() -> pd.DataFrame:
    return _read("islet_geometry.parquet")


@pytest.mark.parametrize("name,column", sorted(ARTIFACTS.items()))
def test_artifact_covers_every_donor_the_geometry_has(name, column, geometry, expected):
    """No artifact may stop short of the islets the geometry already describes.

    `islet_geometry` is the floor because it is what the app draws an islet FROM; an artifact
    that covers fewer donors leaves the app with regions it can outline and nothing to say
    about them, which is exactly how ten donors disappeared.
    """
    built = _donors(_read(name), column)
    geometry_donors = _donors(geometry, "donor_id") & expected
    missing = sorted(geometry_donors - built, key=int)
    assert not missing, (
        f"{name} covers {len(built)} donors; {len(geometry_donors)} have islet geometry. "
        f"Never built or rebuilt for: {missing}")


def test_the_geometry_itself_covers_the_resolved_cohort(geometry, expected):
    missing = sorted(expected - _donors(geometry, "donor_id"), key=int)
    assert not missing, f"islet_geometry is missing donors cohort.resolve() returns: {missing}"


def test_no_artifact_carries_a_donor_the_cohort_excludes(expected):
    """A donor nothing resolves is a donor no builder was asked for -- including 6457.

    Read from `DONOR_EXCLUSIONS` rather than named, so the day an exclusion is lifted this
    keeps testing the rule instead of failing on the decision.
    """
    from phenocycler.cohort import DONOR_EXCLUSIONS

    for name, column in sorted(ARTIFACTS.items()):
        built = _donors(_read(name), column)
        assert not (built & set(DONOR_EXCLUSIONS)), (
            f"{name} carries centrally excluded donor(s): "
            f"{sorted(built & set(DONOR_EXCLUSIONS))}")
        assert not (built - expected), (
            f"{name} carries donor(s) no builder resolves: {sorted(built - expected)}")


def test_the_islet_artifacts_describe_the_same_islets():
    """Internal coherence: the per-islet tables are one islet set seen four ways.

    Not a reproduction of a previous result -- a statement that these files can be joined. The
    app keys them all on `combined_islet_id`, so a key in one and not another silently drops
    rows out of whichever merge reaches it first.
    """
    composition = set(_read("islet_composition_long.parquet")["combined_islet_id"])
    keyed = {
        "islet_marker_means.parquet": set(_read("islet_marker_means.parquet")["combined_islet_id"]),
        "islet_celltype_proportions.csv":
            set(_read("islet_celltype_proportions.csv")["combined_islet_id"]),
        "parquet/islets.parquet": set(_read("parquet/islets.parquet")["combined_islet_id"]),
    }
    for name, keys in keyed.items():
        assert keys == composition, (
            f"{name} and islet_composition_long disagree on {len(keys ^ composition):,} islets "
            f"({len(keys - composition):,} only here, {len(composition - keys):,} only there)")


def test_marker_means_use_the_membership_the_other_artifacts_use():
    """The islet a marker mean is filed under must be the islet that key names elsewhere.

    The aligned per-cell table carries the ORIGINAL export's `islet_num`, which the batch-2
    re-export renumbered. Deriving `islet_key` from it agreed with the rebuilt membership on
    0.00-0.15% of the shared cells for those donors, so the means were filed under the wrong
    islets while every key still looked valid. The check is a join, not a stored number.
    """
    means = _read("islet_marker_means.parquet")
    rebuilt = means["case_id"].astype(str) + "_" + means["islet_key"].astype(str)
    assert (rebuilt == means["combined_islet_id"].astype(str)).all(), (
        "combined_islet_id does not agree with case_id + islet_key in islet_marker_means")

    geometry = _read("islet_geometry.parquet")
    known = set(geometry["combined_islet_id"])
    unknown = sorted(set(means["combined_islet_id"]) - known)[:10]
    assert not unknown, f"marker means filed under islets the geometry does not have: {unknown}"


def test_proportions_are_over_a_denominator_that_is_present_not_assumed():
    """A proportion over an empty region is unavailable, never 0.0.

    Writing 0.0 there states "0% of this islet's peri cells are Alpha" about an islet with no
    peri cells to be 0% of, and pulls every pooled median toward zero.
    """
    frame = _read("islet_celltype_proportions.csv")
    for total, prefix in (("total_cells_core", "prop_"), ("total_cells_peri", "peri_prop_")):
        if total not in frame.columns:
            continue
        empty = frame[frame[total].fillna(0) == 0]
        columns = [c for c in frame.columns if c.startswith(prefix)]
        if empty.empty or not columns:
            continue
        assert empty[columns].isna().all().all(), (
            f"{prefix}* is filled in on {len(empty)} islet(s) with no cells in that region")


# ---------------------------------------------------------------------------------------------
# Staleness: covering the cohort is not the same as having been built from what is on disk now
# ---------------------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def driver():
    """`scripts/senior/build_layers.py`, loaded by path -- it is a script, not a package module."""
    if not DRIVER.is_file():
        pytest.skip(f"{DRIVER} is not present; there is no declared graph to check against")
    spec = importlib.util.spec_from_file_location("build_layers", DRIVER)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def builders(driver):
    """artifact path -> the node(s) declared to write it, and one measured status per node.

    Derived from the declaration, so an artifact that changes builder, or gains one, is followed
    here without an edit. Nothing in this fixture opens an artifact; it reads mtimes only.
    """
    if not APP.is_dir():
        pytest.skip(f"{APP} is not built")
    layers = Path(os.environ.get("ISLET_LAYERS_DIR", "/data/islet_layers"))
    nodes = driver.build_graph(REPO, layers)
    statuses = {status.name: status for status in driver.check(nodes)}
    return driver.producers(nodes), statuses


def _when(mtime):
    return time.strftime("%Y-%m-%d %H:%M", time.localtime(mtime)) if mtime else "unknown"


def _relative(path):
    try:
        return str(Path(path).relative_to(REPO))
    except (ValueError, TypeError):
        return str(path)


@pytest.mark.parametrize("name", sorted(ARTIFACTS) + ["islet_geometry.parquet"])
def test_artifact_was_built_after_every_input_its_builder_reads(name, driver, builders):
    """An artifact older than an input it reads describes a generation that is no longer there.

    This is the check that was missing while five layers and three column families sat hours
    behind their own sources and every other test passed. It asserts nothing about the numbers in
    the artifact -- only that they were computed from the files that are on disk now.
    """
    if not (APP / name).exists():
        pytest.skip(f"{APP / name} is not built")
    by_output, statuses = builders
    path = os.path.normpath(str(APP / name))
    writers = by_output.get(path)
    assert writers, (
        f"{name} is covered by this suite but no node in {_relative(DRIVER)} declares it as an "
        f"output, so nothing can tell whether it is stale. Declare its builder there.")

    offenders = []
    for writer in writers:
        status = statuses[writer]
        if status.state == driver.FRESH:
            continue
        if status.state == driver.MISSING_INPUT:
            offenders.append(f"{writer} declares inputs that do not exist: "
                             + ", ".join(_relative(p) for p in status.missing))
        elif status.state == driver.MISSING_OUTPUT:
            offenders.append(f"{writer} has not written: "
                             + ", ".join(_relative(p) for p in status.missing))
        else:
            offenders.append(
                f"{writer} wrote {_relative(status.oldest_output)} "
                f"({_when(status.oldest_output_mtime)}) BEFORE it read "
                f"{_relative(status.newest_input)} ({_when(status.newest_input_mtime)})")
    assert not offenders, f"{name} is stale:\n  " + "\n  ".join(offenders)
