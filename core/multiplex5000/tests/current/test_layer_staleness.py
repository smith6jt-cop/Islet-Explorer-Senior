"""An artifact built before an input it reads is stale, and nothing said so.

Five of the nine artifacts under ``/data/islet_layers`` were built hours before inputs they read,
and three column families in ``app/parquet/islets.parquet`` sat at a smaller donor set than the
files they are merged from.  Every existing test passed the whole time: each artifact was
internally consistent, every path resolved, every column existed.  What nobody compared was an
output's mtime against its inputs' mtimes, and one artifact's donor coverage against the coverage
of the file it was built out of.

This is that comparison.  It reads the graph ``scripts/senior/build_layers.py`` declares -- the
same declaration ``--check`` reads -- so a builder that gains an input gains it here too, and a
node added there is checked here without an edit.  Nothing in this file asserts a biological
number, reruns a builder, or writes anything.

Two things it deliberately does NOT do:

* It never names a donor count.  Both sides of every coverage assertion are derived: the family's
  coverage from ``islets.parquet``, the expected coverage from the source file that family is
  merged in from, and that source's own floor from ``islet_geometry`` intersected with
  ``cohort.resolve()``.  A test that pinned "20" would have to be edited the day a donor joins,
  which is the same failure in a different file -- the literal would be updated to match whatever
  was built, and would then agree with a stale build forever.
* It never treats a passing check as a statement about the numbers in an artifact.  Fresh means
  the artifact was computed from the generation of the inputs that is on disk now.  It is not a
  claim that what it computed is right.

Marked ``data`` like the other artifact tests, so a synthetic CI run without the data tree skips
rather than fails.

    .venv/bin/python -m pytest core/multiplex5000/tests/current/test_layer_staleness.py -q
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
DRIVER = REPO / "scripts/senior/build_layers.py"
TYPING = REPO / "data/development/typing"
APP = TYPING / "app"
#: Resolved the way every card and the driver resolve it, so the test reads what the app reads.
LAYERS = Path(os.environ.get("ISLET_LAYERS_DIR", "/data/islet_layers"))

sys.path.insert(0, str(REPO / "core/multiplex5000/scripts/typing"))

pytestmark = pytest.mark.data


# ---------------------------------------------------------------------------------------------
# The declared graph, loaded from the file that declares it
# ---------------------------------------------------------------------------------------------

def _load_driver():
    """Import ``build_layers.py`` by path. It is a script, not a package module."""
    spec = importlib.util.spec_from_file_location("build_layers", DRIVER)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    # Registered before exec so the frozen dataclasses in it resolve their own module.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _layer_names() -> tuple[str, ...]:
    """The layer node names, for parametrisation at collection time.

    Empty when the driver is absent, which pytest reports as a skip rather than as a silent pass.
    """
    if not DRIVER.is_file():
        return ()
    try:
        return tuple(_load_driver().LAYER_NODES)
    except Exception:                                    # pragma: no cover - driver is broken
        return ()


LAYER_NAMES = _layer_names()


@pytest.fixture(scope="module")
def driver():
    if not DRIVER.is_file():
        pytest.skip(f"{DRIVER} is not present; there is no declared graph to check")
    return _load_driver()


@pytest.fixture(scope="module")
def nodes(driver):
    if not LAYERS.is_dir():
        pytest.skip(f"{LAYERS} is not mounted")
    if not APP.is_dir():
        pytest.skip(f"{APP} is not built")
    return driver.build_graph(REPO, LAYERS)


@pytest.fixture(scope="module")
def statuses(driver, nodes):
    """One NodeStatus per node, measured once. Reads mtimes only; opens no artifact."""
    return {status.name: status for status in driver.check(nodes)}


def _stamp(mtime: float | None) -> str:
    return time.strftime("%Y-%m-%d %H:%M", time.localtime(mtime)) if mtime else "unknown"


def _short(path: str | None) -> str:
    if not path:
        return "-"
    try:
        return str(Path(path).relative_to(REPO))
    except ValueError:
        return str(path)


def _pair(status) -> str:
    """The offending pair, named: which output, which input, and how far apart."""
    gap = ""
    if status.newest_input_mtime and status.oldest_output_mtime:
        hours = (status.newest_input_mtime - status.oldest_output_mtime) / 3600.0
        gap = f", built {hours:,.1f} h before it"
    return (f"{status.name}: output {_short(status.oldest_output)} "
            f"({_stamp(status.oldest_output_mtime)}) is OLDER than input "
            f"{_short(status.newest_input)} ({_stamp(status.newest_input_mtime)}){gap}")


# ---------------------------------------------------------------------------------------------
# Every declared output is newer than every declared input
# ---------------------------------------------------------------------------------------------

def test_no_declared_output_is_older_than_a_declared_input(driver, statuses):
    """The whole graph, app artifacts and layers alike.

    ``check_node`` compares each node's OLDEST output against its NEWEST input, which is exactly
    "every output is newer than every input" measured in two numbers.  A builder's own source file
    is a declared input, so a builder edited after its artifact was written reads as stale here
    too -- the artifact no longer describes that builder.
    """
    stale = [s for s in statuses.values() if s.state == driver.STALE]
    assert not stale, (
        f"{len(stale)} of {len(statuses)} declared nodes were built before an input they read:\n  "
        + "\n  ".join(_pair(s) for s in stale)
        + "\n\nStale is not a claim about the numbers in these artifacts. It is a claim that they "
          "were computed from a different generation of the inputs than the ones on disk now.")


def test_every_declared_input_and_output_is_present(driver, statuses):
    """A gate that passed on an absent artifact would pass before the build had happened."""
    absent = [s for s in statuses.values()
              if s.state in (driver.MISSING_INPUT, driver.MISSING_OUTPUT)]
    assert not absent, (
        "declared paths that do not exist:\n  "
        + "\n  ".join(f"{s.name} [{s.state}]: " + ", ".join(_short(p) for p in s.missing)
                      for s in absent))


# ---------------------------------------------------------------------------------------------
# The same, for the nine artifacts under /data/islet_layers
# ---------------------------------------------------------------------------------------------

@pytest.mark.parametrize("layer", LAYER_NAMES)
def test_the_layer_artifact_is_newer_than_every_input_it_reads(layer, driver, nodes, statuses):
    """One case per layer, so a stale layer fails under its own name.

    These are the artifacts the six layer cards draw from directly.  A card showing a layer built
    before the typing it describes is showing the previous generation's numbers under this
    generation's caption, and nothing in the app can tell the difference.
    """
    status = statuses[layer]
    if status.state == driver.MISSING_OUTPUT:
        pytest.fail(f"{layer} has no artifact on disk: "
                    + ", ".join(_short(p) for p in status.missing))
    if status.state == driver.MISSING_INPUT:
        pytest.fail(f"{layer} declares inputs that do not exist: "
                    + ", ".join(_short(p) for p in status.missing))
    assert status.state == driver.FRESH, _pair(status)


def test_every_layer_the_driver_declares_has_a_root_on_disk(nodes):
    """The nine layer roots, derived from the declaration rather than listed here.

    The root is the first path component under ``ISLET_LAYERS_DIR``, so a declared output with a
    glob in its directory part (``pln/donor_id=*/...``) still names the directory it lives under.
    """
    missing = []
    for name in LAYER_NAMES:
        roots = set()
        for output in nodes[name].outputs:
            try:
                parts = Path(output).relative_to(LAYERS).parts
            except ValueError:
                continue                      # an output this layer writes outside the tree
            if parts:
                roots.add(LAYERS / parts[0])
        for root in sorted(roots):
            if not root.is_dir():
                missing.append(f"{name}: {root}")
    assert not missing, "layer roots the driver declares and the tree does not have:\n  " + \
        "\n  ".join(missing)


# ---------------------------------------------------------------------------------------------
# The column families in islets.parquet, against the files they are merged from
# ---------------------------------------------------------------------------------------------

#: The join keys of every per-islet artifact; never part of a column family.
KEYS = ("case_id", "islet_key", "combined_islet_id", "donor_id", "region")


def _passthrough_columns(source: pd.DataFrame, islets: pd.DataFrame) -> list[str]:
    """Columns the merge carries across under their own names."""
    return [c for c in source.columns if c not in KEYS and c in islets.columns]


def _region_prefixed_columns(source: pd.DataFrame, islets: pd.DataFrame) -> list[str]:
    """The galectin merge splits one `region` column into two prefixes.

    Derived exactly the way ``build_parquet_views.py`` derives it (the class list comes off the
    source's own ``prop_*`` columns), so a new galectin class is picked up here without an edit.
    """
    classes = [c[len("prop_"):] for c in source.columns if c.startswith("prop_")]
    names = [f"{prefix}{cls}" for cls in classes for prefix in ("prop_", "peri_prop_")]
    return [c for c in names if c in islets.columns]


#: family -> (the file islets.parquet merges those columns in from, how they are spelled there).
#: The mapping is a fact about the schema; every count on either side of an assertion is measured.
FAMILIES = {
    "leiden_* (islet-level Leiden and its UMAP)":
        (APP / "islet_leiden.parquet", _passthrough_columns),
    "sc_pseudotime* / sc_n_* (the single-cell trajectory summary)":
        (TYPING / "trajectory/singlecell_islet_summary.parquet", _passthrough_columns),
    "prop_Gal* / peri_prop_* (the galectin layer)":
        (APP / "islet_galectin.parquet", _region_prefixed_columns),
}


def _read(path: Path) -> pd.DataFrame:
    if not path.is_file():
        pytest.skip(f"{path} is not built")
    return pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path)


def _donors(frame: pd.DataFrame, column: str = "case_id") -> set[str]:
    """Read as strings so `115` and `0115` cannot compare unequal."""
    return set(frame[column].astype(str))


@pytest.fixture(scope="module")
def islets() -> pd.DataFrame:
    return _read(APP / "parquet/islets.parquet")


@pytest.fixture(scope="module")
def geometry_floor() -> set[str]:
    """The donors an artifact over these islets must reach: what the geometry has, that the
    cohort resolves.  Both sides derived -- `islet_geometry` is what the app draws an islet FROM,
    and `cohort.resolve()` is the module every builder calls."""
    geometry = _read(APP / "islet_geometry.parquet")
    import cohort as COHORT

    return _donors(geometry, "donor_id") & set(COHORT.resolve(None, include_marginal=True))


@pytest.mark.parametrize("family", sorted(FAMILIES))
def test_a_column_family_reaches_every_donor_its_source_carries(family, islets):
    """The merge must not stop short of the file it merged from.

    This is the half of the defect that is a merge problem rather than a build problem: the source
    was rebuilt at the full cohort and the merge step that carries it into `islets.parquet` was
    never re-run, so the column sits at the previous generation's donor set while the file beside
    it carries the current one.
    """
    path, spell = FAMILIES[family]
    source = _read(path)
    columns = spell(source, islets)
    assert columns, f"{family}: none of {path.name}'s columns are in islets.parquet at all"

    carried = islets.loc[islets[columns].notna().any(axis=1)]
    present = _donors(carried)
    expected = _donors(source)
    missing = sorted(expected - present, key=int)
    assert not missing, (
        f"{family}: {len(columns)} column(s) filled in for {len(present)} donor(s) in "
        f"islets.parquet, while {_short(str(path))} carries {len(expected)}. "
        f"Never merged for: {missing}. "
        f"islets.parquet written {_stamp((APP / 'parquet/islets.parquet').stat().st_mtime)}, "
        f"{path.name} written {_stamp(path.stat().st_mtime)}.")


@pytest.mark.parametrize("family", sorted(FAMILIES))
def test_the_source_of_a_column_family_covers_the_islets_the_geometry_describes(
        family, geometry_floor):
    """The other half: a source that was itself built on a smaller cohort.

    Equality against the source is not enough on its own -- a family merged faithfully from a
    ten-donor file agrees with it perfectly and is still ten donors in the app.  `islet_geometry`
    is the floor because it is what the app draws an islet from; a per-islet artifact covering
    fewer donors leaves the app with regions it can outline and nothing to say about them.
    """
    path, _ = FAMILIES[family]
    source = _read(path)
    missing = sorted(geometry_floor - _donors(source), key=int)
    assert not missing, (
        f"{family}: {_short(str(path))} covers {len(_donors(source))} donor(s); "
        f"{len(geometry_floor)} have islet geometry. Never built or rebuilt for: {missing}. "
        f"Written {_stamp(path.stat().st_mtime)}.")
