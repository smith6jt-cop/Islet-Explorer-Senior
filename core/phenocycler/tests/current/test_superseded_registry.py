"""Superseded artifacts must refuse to be read, and say so with the numbers.

The defect this guards is not a wrong calculation.  It is a RIGHT calculation over the wrong
generation of the same table -- ``typing/subpanel/subpanel_calls.parquet`` read as the current
endocrine typing, which is 79,572 cells and three materially different endocrine counts away
from the build on disk.  Nothing about that read looks wrong while it is happening.  Prose said
the file was superseded; prose does not open files, so prose did not stop it.

So what is asserted here is refusals, and the CONTENT of the refusals:

  * every registry entry marked ``REFUSE`` raises when a reader touches it, before any existence
    or schema check, and by path, by basename, and through a symlink pointing at it;
  * the message names what was refused, when, by what, why, and what to read instead -- a
    refusal a reader cannot act on just sends them back to the old file;
  * the measured counts are IN the message, both generations of them, so the reader does not
    have to trust the word "superseded";
  * the current artifacts beside them stay readable -- ``typing/cells``,
    ``subpanel/subpanel_log.json``, ``broadtype/scratch/cosine_scale.parquet``.  A guard that
    takes the live inputs with it is not a guard, it is an outage;
  * a count nobody measured says UNAVAILABLE.  It is never 0, and never an estimate;
  * the one artifact that is historical but is the ONLY copy of what the lane needs -- the
    pre-2026-09-16 ``data/runs/`` tree -- announces instead of raising.  Loud, not silent, and
    not a method decision taken in a guard.

Nothing here asserts against a frozen copy of a result, and nothing passes by reproducing a
biological finding.  The 20-donor counts appear only as the reason a refusal exists.

Run with
    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python -m pytest -p no:cacheprovider \
      core/phenocycler/tests/current/test_superseded_registry.py -q
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pyarrow.parquet as pq
import pytest

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO / "core/phenocycler/scripts/typing"))

import subpanel_io as S      # noqa: E402
import superseded as SUP     # noqa: E402

MODULE = REPO / "core/phenocycler/scripts/typing/superseded.py"

REFUSED = [r for r in SUP.REGISTRY if r.severity == SUP.REFUSE]
ANNOUNCED = [r for r in SUP.REGISTRY if r.severity == SUP.ANNOUNCE]


@pytest.fixture(autouse=True)
def _forget_announcements():
    """Each test gets a fresh process-memory, so "once per process" is testable."""
    SUP._announced.clear()
    yield
    SUP._announced.clear()


# --------------------------------------------------------------------------------------------
# The refusal itself
# --------------------------------------------------------------------------------------------

@pytest.mark.parametrize("record", REFUSED, ids=[r.key for r in REFUSED])
def test_every_refused_fragment_raises_by_path(record):
    """Not one of them may be openable, and the artifact need not exist to be refused."""
    for fragment in record.fragments:
        target = REPO / fragment
        with pytest.raises(SUP.SupersededArtifactError):
            SUP.refuse_if_superseded(target)
        # A tree's children are refused too, existing or not.
        if record.kind == "tree":
            with pytest.raises(SUP.SupersededArtifactError):
                SUP.refuse_if_superseded(target / "donor_id=115" / "data.parquet")


@pytest.mark.parametrize("record", REFUSED, ids=[r.key for r in REFUSED])
def test_the_message_names_the_replacement_and_where_the_table_is(record):
    """A refusal a reader cannot act on sends them straight back to the superseded file."""
    with pytest.raises(SUP.SupersededArtifactError) as exc:
        SUP.refuse_if_superseded(REPO / record.fragments[0])
    message = str(exc.value)
    assert record.fragments[0] in message                   # what was refused
    assert record.superseded_on in message                  # when
    assert record.superseded_by in message                  # by what
    assert record.read_instead in message                   # what to read instead
    assert "superseded.py" in message                       # where the one table lives
    assert record.key in message


@pytest.mark.parametrize("record", SUP.REGISTRY, ids=[r.key for r in SUP.REGISTRY])
def test_a_count_nobody_measured_says_unavailable_not_zero(record):
    """A missing measurement is unavailable. It is never 0 and never an estimate."""
    counts = record.counts
    assert counts.strip(), f"{record.key} carries no measurement at all"
    has_number = any(ch.isdigit() for ch in counts)
    assert has_number or "UNAVAILABLE" in counts, record.key
    assert "approximately" not in counts.lower()
    assert "roughly" not in counts.lower()


# --------------------------------------------------------------------------------------------
# subpanel_calls.parquet -- the artifact this whole module exists for
# --------------------------------------------------------------------------------------------

SUBPANEL_CALLS = "data/development/typing/subpanel/subpanel_calls.parquet"


def test_subpanel_calls_refusal_carries_both_generations_of_counts():
    with pytest.raises(SUP.SupersededArtifactError) as exc:
        SUP.refuse_if_superseded(REPO / SUBPANEL_CALLS)
    message = str(exc.value)
    # The previous generation, and the build that replaced it, side by side.
    assert "25,150,537" in message and "25,230,109" in message
    assert "79,572" in message
    assert "Beta 191,485 / Alpha 246,055 / Delta 62,797" in message
    assert "Beta 226,183 / Alpha 279,686 / Delta 90,658" in message
    assert "2026-09-11" in message and "2026-09-15" in message
    assert "data/development/typing/cells/" in message


def test_the_counts_in_the_message_are_the_counts_on_disk():
    """The refusal quotes a measurement; the measurement has to be true of this checkout.

    A refusal carrying a stale number is a refusal nobody will believe the second time.
    """
    calls = REPO / SUBPANEL_CALLS
    cells = REPO / "data/development/typing/cells"
    if not calls.exists() or not cells.exists():
        pytest.skip("the superseded artifact or its replacement is not in this checkout")
    superseded_rows = pq.ParquetFile(calls).metadata.num_rows
    partitions = sorted(cells.glob("donor_id=*/*.parquet"))
    current_rows = sum(pq.ParquetFile(p).metadata.num_rows for p in partitions)
    record = SUP._BY_KEY["subpanel_calls"]
    assert f"{superseded_rows:,}" in record.counts
    assert f"{current_rows:,}" in record.counts
    assert f"{current_rows - superseded_rows:,}" in record.counts


def test_an_alias_or_symlink_to_it_is_the_same_file(tmp_path):
    """Copying a superseded artifact under another name is how one gets read as current."""
    target = tmp_path / "subpanel_calls.parquet"
    target.touch()
    alias = tmp_path / "cells.parquet"
    alias.symlink_to(target)
    with pytest.raises(SUP.SupersededArtifactError, match="Beta 191,485"):
        SUP.refuse_if_superseded(alias)
    # And the bare basename, anywhere at all.
    with pytest.raises(SUP.SupersededArtifactError):
        SUP.refuse_if_superseded(tmp_path / "somewhere" / "subpanel_calls.parquet")


# --------------------------------------------------------------------------------------------
# refinement_log.json -- superseded AND unreproducible, which is worse
# --------------------------------------------------------------------------------------------

def test_refinement_log_refusal_says_no_producer_and_names_what_to_read():
    with pytest.raises(SUP.SupersededArtifactError) as exc:
        SUP.refuse_if_superseded(REPO / "data/development/typing/umap/refinement_log.json")
    message = str(exc.value)
    assert "NO PRODUCER" in message
    assert "2026-09-11" in message
    assert "subpanel_log.json" in message          # the current file to read instead
    assert "UNAVAILABLE" in message                # there is no current own-marker table
    assert "NOT WIRED" in message                  # refine.py is not a replacement either
    assert "54, not 21" in message                 # and its marker set is a different one


def test_the_current_classification_log_beside_it_stays_readable():
    """The replacement the refusal points at must not itself be refused."""
    current = S.OUT / "subpanel/subpanel_log.json"
    assert SUP.match(current) is None
    if current.exists():
        assert json.loads(current.read_text())


# --------------------------------------------------------------------------------------------
# The guard is attached to the paths, so a reader gets it without asking
# --------------------------------------------------------------------------------------------

def test_subpanel_io_hands_out_guarded_paths():
    assert isinstance(S.OUT, SUP.GuardedPath)
    assert isinstance(S.RUN, SUP.GuardedPath)
    assert isinstance(S.OUT / "umap" / "refinement_log.json", SUP.GuardedPath)


def test_a_reader_that_never_heard_of_the_registry_is_still_refused():
    """`S.OUT / ...` is what all thirty-three modules in this lane compose paths from."""
    path = S.OUT / "umap/refinement_log.json"
    with pytest.raises(SUP.SupersededArtifactError):
        path.read_text()
    with pytest.raises(SUP.SupersededArtifactError):
        os.fspath(path)                       # pyarrow / pandas take this door
    with pytest.raises(SUP.SupersededArtifactError):
        path.open("rb")
    with pytest.raises(SUP.SupersededArtifactError):
        pq.read_table(S.OUT / "subpanel/subpanel_calls.parquet")


def test_the_current_artifacts_are_untouched():
    """A guard that takes the live inputs with it is an outage, not a guard."""
    for live in (S.OUT / "cells",
                 S.OUT / "cells/donor_id=115/data.parquet",
                 S.OUT / "subpanel/subpanel_log.json",
                 S.OUT / "subpanel_rebuilt",
                 S.OUT / "islet_membership",
                 S.OUT / "immune",
                 S.SCALE_CSV):
        assert SUP.match(live) is None, live
        os.fspath(live)                        # must not raise


def test_broadtype_scratch_is_exempt_because_the_scale_is_read_from_it():
    """`subpanel_io.load_scale()` reads the cohort floor/top the classification was made on."""
    assert "broadtype/scratch" in " ".join(SUP.EXEMPT)
    assert SUP.match(S.SCALE_CSV) is None
    assert SUP.match(REPO / "data/development/broadtype/scratch/anything.parquet") is None
    # An exemption beats the basename fallback too, or a cache file that happened to be called
    # subpanel_calls.parquet would take the live scale root down with it.
    assert SUP.match(REPO / "data/development/broadtype/scratch/subpanel_calls.parquet") is None
    # its siblings in the withdrawn comparison are still registered
    assert SUP.match(REPO / "data/development/broadtype/release") is not None


def test_a_segment_is_matched_not_a_substring():
    """`data/restore_gated_redsea` must not swallow a differently named neighbour."""
    assert SUP.match(REPO / "data/restore_gated_redsea_v2/x.parquet") is None
    assert SUP.match(REPO / "data/development/celltype/calls") is None
    assert SUP.match(REPO / "data/development/typing/calls") is not None


# --------------------------------------------------------------------------------------------
# ANNOUNCE: loud, not silent, and not a method decision taken inside a guard
# --------------------------------------------------------------------------------------------

@pytest.mark.parametrize("record", ANNOUNCED, ids=[r.key for r in ANNOUNCED])
def test_announced_artifacts_read_but_say_so(record, capsys):
    target = REPO / record.fragments[0] / "expression"
    assert SUP.refuse_if_superseded(target) == target
    printed = capsys.readouterr().err
    assert record.superseded_on in printed
    assert record.read_instead in printed
    assert record.counts in printed


def test_the_historical_run_tree_announces_once_and_keeps_working(capsys):
    """The 20-donor cohort output lives nowhere else; refusing it would stop the lane."""
    first = S.RUN / "expression/donor_id=115/data.parquet"
    os.fspath(first)
    assert "Reading historical (not current)" in capsys.readouterr().err
    os.fspath(S.RUN / "expression/donor_id=6442/data.parquet")
    assert capsys.readouterr().err == ""          # once per process, not once per donor
    path, which = S.cell_source("115")
    assert which in ("run", "revised")


# --------------------------------------------------------------------------------------------
# The one allowed reader
# --------------------------------------------------------------------------------------------

def _call_as(filename: str, function: str, path):
    """Run `refuse_if_superseded(path)` from a frame that really is <filename>::<function>."""
    source = f"def {function}(p, f):\n    return f(p)\n"
    namespace: dict = {}
    exec(compile(source, filename, "exec"), namespace)
    return namespace[function](path, SUP.refuse_if_superseded)


def test_only_subpanel_gate_may_read_the_previous_generation(capsys):
    """CLAUDE.md keeps the pooled levels alive so `subpanel.py --gate` can reproduce the table.

    Allowed is not the same as silent: the gate still gets the whole notice on stderr.
    """
    target = REPO / SUBPANEL_CALLS
    gate_file = str(REPO / "core/phenocycler/scripts/typing/subpanel.py")
    assert _call_as(gate_file, "gate", target) == target
    printed = capsys.readouterr().err
    assert "allowed for" in printed
    assert "Beta 191,485 / Alpha 246,055 / Delta 62,797" in printed

    # The same file, a different function: refused.
    with pytest.raises(SUP.SupersededArtifactError):
        _call_as(gate_file, "main", target)
    # The same function name, a different file: refused.
    with pytest.raises(SUP.SupersededArtifactError):
        _call_as(str(REPO / "scripts/senior/build_layers.py"), "gate", target)


def test_nothing_else_in_the_registry_has_an_allowed_reader():
    allowed = {r.key: r.allowed_readers for r in SUP.REGISTRY if r.allowed_readers}
    assert set(allowed) == {"subpanel_calls"}


# --------------------------------------------------------------------------------------------
# The senior readers that used to open these without a word
# --------------------------------------------------------------------------------------------

def test_the_distance_builder_refuses_through_the_shared_table():
    sys.path.insert(0, str(REPO))
    from scripts.senior import build_distance_profiles as builder

    with pytest.raises(ValueError) as exc:
        builder.build_distance_profiles(cells_root=Path(SUBPANEL_CALLS))
    message = str(exc.value)
    assert "Beta 191,485 / Alpha 246,055 / Delta 62,797" in message
    assert "Beta 226,183 / Alpha 279,686 / Delta 90,658" in message
    # It refuses the other registered roots now too, which it never knew about before.
    with pytest.raises(ValueError):
        builder.refuse_superseded(REPO / "data/restore_gated_redsea")


@pytest.mark.parametrize("script,args", [
    ("scripts/senior/build_cell_distance_parquet.py", []),
    ("scripts/senior/build_phenotyped_h5ad.py", []),
])
def test_the_retired_senior_producers_refuse_at_the_command_line(script, args):
    """Their correction banners are prose. This is the part that stops the run."""
    result = subprocess.run([sys.executable, "-B", str(REPO / script), *args],
                            capture_output=True, text=True, check=False, cwd=REPO)
    assert result.returncode != 0, result.stdout
    combined = result.stdout + result.stderr
    assert "superseded.py" in combined
    assert "read instead" in combined


def test_the_module_answers_on_its_own_at_the_command_line():
    """One place to ask what is superseded, and it prints the counts."""
    listing = subprocess.run([sys.executable, "-B", str(MODULE)],
                             capture_output=True, text=True, check=False, cwd=REPO)
    assert listing.returncode == 0
    assert "subpanel_calls" in listing.stdout
    assert "25,150,537" in listing.stdout

    checked = subprocess.run([sys.executable, "-B", str(MODULE), SUBPANEL_CALLS,
                              "data/development/typing/cells"],
                             capture_output=True, text=True, check=False, cwd=REPO)
    assert checked.returncode == 1
    assert "REFUSED" in checked.stdout
    assert "ok        data/development/typing/cells" in checked.stdout
