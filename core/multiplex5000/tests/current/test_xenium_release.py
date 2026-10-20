"""Contracts for the one Xenium release resolver.

Five builders used to name five different Xenium sources and only one of them
read the current pair.  Everything asserted here is a way that could happen
again: a directory mtime deciding which release is current, a superseded release
winning on a timestamp, a tree with no ``calls/`` silently losing a comparison to
an older one, an absent release reaching a count as a zero, a donor id inferred
into an integer, or the two nucleus absences collapsing into one.

Nothing here checks that a number matches a previous biological result.  The
quantitative checks are internal coherence -- the parts of a tally sum to the rows
it read, and a release's own book-keeping of its own file agrees with reading that
file -- which are statements about the build, not about the biology.
"""

from __future__ import annotations

import ast
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq
import pytest

REPO = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(REPO))

from scripts.senior import xenium_release as xr  # noqa: E402

#: build_availability.py left this list on 2026-09-23: it dropped xenium_panc and xenium_pln, which
#: no card read, and reads no release at all now (test_declared_producers.py holds it to that).
BUILDERS = {
    "build_compartment_gradient.py": REPO / "scripts/senior/build_compartment_gradient.py",
    "modality_concordance.py": REPO / "scripts/senior/modality_concordance.py",
    "build_grid_composition.py": REPO / "scripts/senior/build_grid_composition.py",
}


# ---------------------------------------------------------------------------
# synthetic releases
# ---------------------------------------------------------------------------
def write_release(root: Path, *, release_id: str, when: str, tissue: str | None = "pln",
                  roi="pln_1", donors=("112", "0112"), supersedes: str | None = None,
                  calls: bool = True, labels=None, partial: bool = False,
                  donor_column: bool = False, time_key: str = "generated_at_utc",
                  nucleus: bool = True, calls_dir: str | None = None,
                  sidecars: dict | None = None) -> Path:
    """One release on disk, shaped like the real ones.

    ``donor_column`` reproduces the pLN releases, which carry ``donor_id`` and
    ``roi`` as physical ``large_string`` columns as well as in the hive path.
    """
    root.mkdir(parents=True, exist_ok=True)
    declared: dict = {"schema_version": 1, "release_id": release_id, "modality": "xenium",
                      "roi": roi, time_key: when, "donors": list(donors)}
    if tissue is not None:
        declared["tissue"] = tissue
    if calls_dir is not None:
        declared["calls_dir"] = calls_dir
    if supersedes is not None:
        declared["producer"] = {"lane": release_id, "supersedes": supersedes}
    if partial:
        declared["scope"] = {"is_partial": True, "slug": "one_section"}
    (root / "release.json").write_text(json.dumps(declared))
    for name, text in (sidecars or {}).items():
        (root / name).write_text(text)
    if not calls:
        return root
    labels = labels or ["Immune_T", "Non_immune", "Uncalled_ambiguous", "Unavailable_no_nucleus"]
    for donor in donors:
        part = root / (calls_dir or "calls") / f"donor_id={donor}" / f"roi={roi}"
        part.mkdir(parents=True, exist_ok=True)
        table = {"cell_id": [f"{donor}-{i}" for i in range(len(labels))], "lineage": labels}
        if donor_column:
            table["donor_id"] = pa.array([donor] * len(labels), type=pa.large_string())
            table["roi"] = pa.array([roi] * len(labels), type=pa.large_string())
        if nucleus:
            table["nucleus_count"] = [1] * (len(labels) - 1) + [0]
            table["nucleus_measured"] = [True] * len(labels)
        pq.write_table(pa.table(table), part / "data.parquet")
    return root


@pytest.fixture
def roots(tmp_path):
    a, b = tmp_path / "repo_tree", tmp_path / "external"
    a.mkdir()
    b.mkdir()
    return [a, b]


def iso(day: int, hour: int = 0) -> str:
    return f"2026-09-{day:02d}T{hour:02d}:00:00Z"


# ---------------------------------------------------------------------------
# selection is on the declaration, never on an mtime
# ---------------------------------------------------------------------------
def test_the_module_reads_no_filesystem_mtime_at_all():
    """A touch must not be able to change which release the cohort is built from.

    The scan this replaces ordered candidates by ``manifest.stat().st_mtime``.
    """
    source = (REPO / "scripts/senior/xenium_release.py").read_text()
    for forbidden in ("st_mtime", "getmtime", "st_ctime", "os.stat("):
        assert forbidden not in source, f"{forbidden} reappeared in the resolver"


def test_a_newer_mtime_does_not_make_an_older_release_current(roots):
    old = write_release(roots[0] / "old", release_id="old", when=iso(10))
    write_release(roots[0] / "new", release_id="new", when=iso(17))
    # Reorder the mtimes: the OLD release becomes the newest thing on disk.
    future = datetime(2027, 1, 1, tzinfo=timezone.utc).timestamp()
    for path in (old, old / "release.json"):
        os.utime(path, (future, future))
    found = xr.resolve("pln", roots=roots)
    assert found.available
    assert found.release_id == "new"
    assert found.generated_at.isoformat() == "2026-09-17T00:00:00+00:00"


def test_both_roots_are_searched(roots):
    write_release(roots[0] / "in_repo", release_id="in_repo", when=iso(10))
    write_release(roots[1] / "outside", release_id="outside", when=iso(17))
    found = xr.resolve("pln", roots=roots)
    assert found.available and found.release_id == "outside"
    assert found.root == roots[1] / "outside"


def test_a_tissue_is_read_from_roi_when_no_tissue_key_is_declared(roots):
    """The pancreas release declares no ``tissue`` at all, only ``roi: "panc"``."""
    write_release(roots[0] / "panc", release_id="panc", when=iso(17), tissue=None, roi="panc")
    found = xr.resolve("pancreas", roots=roots)
    assert found.available and found.tissue == "panc"


def test_the_tissue_aliases_the_other_lanes_use_all_resolve(roots):
    write_release(roots[0] / "pln", release_id="pln", when=iso(17))
    for alias in ("pln", "PLN", "pancreatic_lymph_node", "lymph_node", "pln_1"):
        assert xr.resolve(alias, roots=roots).release_id == "pln"
    with pytest.raises(xr.XeniumReleaseError, match="unknown tissue"):
        xr.resolve("spleen", roots=roots)


def test_a_release_declaring_no_parsable_time_is_rejected_not_mtime_ordered(roots):
    write_release(roots[0] / "nostamp", release_id="nostamp", when="not a date")
    found = xr.resolve("pln", roots=roots)
    assert not found.available
    reasons = " ".join(c["rejected"] or "" for c in found.rejected)
    assert "no parsable write time" in reasons
    assert "mtime" in reasons


# ---------------------------------------------------------------------------
# supersession
# ---------------------------------------------------------------------------
def test_a_superseded_release_loses_even_when_its_timestamp_is_newer(roots):
    """A re-run of a retired lane must not resurrect it."""
    v2 = write_release(roots[0] / "release_pln_immune_v2", release_id="v2", when=iso(18))
    write_release(roots[0] / "release_pln_immune_v3", release_id="v3", when=iso(17),
                  supersedes=f"{v2} (its partition contained nucleus_count == 0 cells)")
    found = xr.resolve("pln", roots=roots)
    assert found.available and found.release_id == "v3"
    # The loser really was the newer one on paper, so the win came from the
    # declaration and not from the ordering.
    retired = next(c for c in xr.scan(roots) if c.release_id == "v2")
    assert retired.generated_at > found.generated_at


def test_supersession_matches_on_the_resolved_path_not_on_a_substring(roots):
    """``release_pln_immune`` is a PREFIX of ``release_pln_immune_v2``.

    A substring test would retire the successor along with its predecessor.
    """
    draft = write_release(roots[0] / "release_pln_immune", release_id="draft", when=iso(17, 1))
    write_release(roots[0] / "release_pln_immune_v2", release_id="v2", when=iso(17, 2),
                  supersedes=f"{draft} (draft; its stage-2 labels were refuted)")
    found = xr.resolve("pln", roots=roots)
    assert found.available and found.release_id == "v2"


def test_supersession_may_name_the_release_id_instead_of_a_path(roots):
    write_release(roots[0] / "a", release_id="lane-a-20260918", when=iso(18))
    write_release(roots[0] / "b", release_id="lane-b-20260917", when=iso(17),
                  supersedes="lane-a-20260918")
    assert xr.resolve("pln", roots=roots).release_id == "lane-b-20260917"


def test_a_relative_supersedes_path_resolves_against_the_repository(roots, monkeypatch):
    """The real v3 writes ``data/development/xenium_typing/release_pln_immune_v2``."""
    monkeypatch.setattr(xr, "REPO", roots[0])
    write_release(roots[0] / "old", release_id="old", when=iso(18))
    write_release(roots[0] / "new", release_id="new", when=iso(17), supersedes="old (superseded)")
    assert xr.resolve("pln", roots=roots).release_id == "new"


# ---------------------------------------------------------------------------
# a release must have a real calls tree behind it
# ---------------------------------------------------------------------------
def test_a_release_without_a_calls_tree_is_rejected_by_name(roots):
    """``release_panc_v1`` holds only ``_work/``; it must not win by being newest."""
    write_release(roots[0] / "no_calls", release_id="no_calls", when=iso(18), calls=False)
    write_release(roots[0] / "with_calls", release_id="with_calls", when=iso(17))
    found = xr.resolve("pln", roots=roots)
    assert found.available and found.release_id == "with_calls"


def test_an_empty_calls_tree_is_not_a_calls_tree(roots):
    root = write_release(roots[0] / "hollow", release_id="hollow", when=iso(18), calls=False)
    (root / "calls" / "donor_id=112").mkdir(parents=True)
    found = xr.resolve("pln", roots=roots)
    assert not found.available
    assert "no parquet file" in " ".join(c["rejected"] or "" for c in found.rejected)


def test_a_partial_run_is_not_the_cohort_release(roots):
    """A partial run writes the SAME release_id under ``_partial/<kind>s=<slug>/``."""
    full = write_release(roots[0] / "full", release_id="same-id", when=iso(17))
    write_release(full / "_partial" / "sections=6347", release_id="same-id",
                  when=iso(18), partial=True)
    found = xr.resolve("pln", roots=roots)
    assert found.available and found.root == full


def test_a_nested_release_is_never_a_candidate(roots):
    """The walk stops at the first manifest, so a release inside a release is not one."""
    full = write_release(roots[0] / "full", release_id="outer", when=iso(17))
    write_release(full / "compare" / "inner", release_id="inner", when=iso(18))
    found = xr.resolve("pln", roots=roots)
    assert found.available and found.release_id == "outer"
    assert [c.release_id for c in xr.scan(roots)] == ["outer"]


def test_a_declared_calls_dir_is_honoured(roots):
    write_release(roots[0] / "named", release_id="named", when=iso(17), calls_dir="assignments")
    found = xr.resolve("pln", roots=roots)
    assert found.available and found.calls_dir.name == "assignments"


# ---------------------------------------------------------------------------
# an absence is named, and cannot read as a zero
# ---------------------------------------------------------------------------
def test_nothing_current_returns_a_named_unavailable_not_none(roots):
    found = xr.resolve("panc", roots=roots)
    assert found is not None
    assert found.available is False
    assert found.reason.strip()
    assert [str(r) for r in found.roots] == [str(r) for r in roots]


def test_an_unavailable_carries_no_root_so_it_cannot_be_read_as_one(roots):
    """Code that assumed a release must fail loudly, not count an absence as zero."""
    found = xr.resolve("panc", roots=roots)
    for attribute in ("root", "calls_dir", "release_id", "donor_partitions"):
        with pytest.raises(AttributeError):
            getattr(found, attribute)
    with pytest.raises(xr.XeniumReleaseError, match="no current Xenium panc release"):
        xr.require("panc", roots=roots)


def test_an_unavailable_names_every_candidate_it_rejected(roots):
    write_release(roots[0] / "no_calls", release_id="no_calls", when=iso(18), calls=False)
    found = xr.resolve("pln", roots=roots)
    assert not found.available
    described = found.describe()
    assert described["available"] is False
    assert described["candidates_rejected"][0]["release_id"] == "no_calls"
    assert "calls" in described["candidates_rejected"][0]["rejected"]


def test_a_foreign_repo_root_never_reaches_the_external_machine_root(tmp_path, monkeypatch):
    """A synthetic tree must not silently read a release nobody put there.

    A builder run against a fixture repo root resolves releases under that root; if
    the external root travelled with it, that fixture would be measuring this
    machine.
    """
    monkeypatch.delenv(xr.ROOTS_ENV, raising=False)
    monkeypatch.delenv(xr.EXTERNAL_ROOT_ENV, raising=False)
    assert xr.default_roots(tmp_path) == (tmp_path / xr.REPO_RELATIVE_ROOT,)
    assert xr.EXTERNAL_ROOT_DEFAULT in xr.default_roots(REPO)
    monkeypatch.setenv(xr.EXTERNAL_ROOT_ENV, str(tmp_path / "named"))
    assert xr.default_roots(tmp_path)[-1] == tmp_path / "named"
    monkeypatch.setenv(xr.ROOTS_ENV, os.pathsep.join([str(tmp_path / "a"), str(tmp_path / "b")]))
    assert xr.default_roots(REPO) == (tmp_path / "a", tmp_path / "b")


# ---------------------------------------------------------------------------
# donor_id is a string, everywhere, always
# ---------------------------------------------------------------------------
def test_a_donor_id_is_read_as_a_string_and_keeps_its_leading_zero(roots):
    write_release(roots[0] / "r", release_id="r", when=iso(17), donors=("112", "0112"))
    release = xr.resolve("pln", roots=roots)
    assert set(release.donor_partitions()) == {"112", "0112"}
    assert all(isinstance(d, str) for d in release.observed_donors())
    table = release.read_calls(columns=["donor_id"], donors=["0112"])
    assert table.num_rows and set(table["donor_id"].to_pylist()) == {"0112"}


def test_inferred_partitioning_is_exactly_what_this_module_refuses_to_do(roots):
    """The contrast, measured: an inferring reader turns "112" into the integer 112.

    The pancreas release warns about this in its own ``partition_key_note``.
    """
    write_release(roots[0] / "r", release_id="r", when=iso(17), donors=("112",))
    release = xr.resolve("pln", roots=roots)
    inferred = ds.dataset(release.calls_dir, format="parquet", partitioning="hive")
    assert pa.types.is_integer(inferred.schema.field("donor_id").type)
    assert pa.types.is_string(release.dataset().schema.field("donor_id").type)


def test_a_physical_large_string_donor_column_does_not_break_the_read(roots):
    """The pLN releases carry donor_id and roi as physical large_string columns.

    Declaring them as ``string`` raises ArrowTypeError: incompatible types.
    """
    write_release(roots[0] / "r", release_id="r", when=iso(17), donors=("112",),
                  donor_column=True)
    release = xr.resolve("pln", roots=roots)
    with pytest.raises(pa.ArrowTypeError):
        ds.dataset(release.calls_dir, format="parquet",
                   partitioning=ds.partitioning(
                       pa.schema([("donor_id", pa.string()), ("roi", pa.string())]),
                       flavor="hive")).to_table(columns=["donor_id"])
    table = release.read_calls(columns=["donor_id"], donors=["112"])
    assert table.num_rows and set(table["donor_id"].to_pylist()) == {"112"}


def test_a_partition_key_that_is_not_a_string_is_refused(roots):
    root = write_release(roots[0] / "r", release_id="r", when=iso(17), donors=("112",))
    part = root / "calls" / "donor_id=112" / "roi=pln_1"
    pq.write_table(pa.table({"cell_id": ["a"], "lineage": ["Immune_T"],
                             "donor_id": pa.array([112], type=pa.int64())}),
                   part / "data.parquet")
    with pytest.raises(xr.XeniumReleaseError, match="stops matching the donor it names"):
        xr.resolve("pln", roots=roots).partitioning()


# ---------------------------------------------------------------------------
# columns
# ---------------------------------------------------------------------------
def test_a_column_role_falls_back_when_the_release_declares_no_map(roots):
    """The pancreas release declares ``columns``; neither pLN release does."""
    write_release(roots[0] / "r", release_id="r", when=iso(17))
    release = xr.resolve("pln", roots=roots)
    assert release.column("label") == "lineage"
    assert release.column("cell_key") == "cell_id"
    assert release.column("common_lineage") is None
    with pytest.raises(xr.XeniumReleaseError, match="no column fills the role"):
        release.require_column("common_lineage")


def test_a_declared_column_map_wins_over_the_fallback(roots):
    root = write_release(roots[0] / "r", release_id="r", when=iso(17))
    declared = json.loads((root / "release.json").read_text())
    declared["columns"] = {"label": "cell_id"}
    (root / "release.json").write_text(json.dumps(declared))
    assert xr.resolve("pln", roots=roots).column("label") == "cell_id"


def test_a_sidecar_is_found_or_named(roots):
    write_release(roots[0] / "r", release_id="r", when=iso(17),
                  sidecars={"per_donor_summary.csv": "donor_id,cells\n112,4\n"})
    release = xr.resolve("pln", roots=roots)
    assert release.sidecar("per_donor_summary.csv").name == "per_donor_summary.csv"
    assert release.sidecar("controls.csv") is None
    with pytest.raises(xr.XeniumReleaseError, match="not on disk"):
        release.require_sidecar("controls.csv")


# ---------------------------------------------------------------------------
# labels and the two nucleus absences
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("label,state", [
    ("Immune_T", xr.CALL),
    ("Non_immune", xr.CALL),               # a measured negative is an ANSWER
    ("Endocrine", xr.CALL),
    ("Uncalled_ambiguous", xr.UNCALLED),
    ("Uncalled_no_enrichment", xr.UNCALLED),
    ("Immune_uncalled", xr.UNCALLED),      # the other spelling, in the same file
    ("Unavailable_no_nucleus", xr.LABEL_UNAVAILABLE),
    ("", xr.BLANK),
    (None, xr.BLANK),
])
def test_a_label_is_classified_by_what_it_says(label, state):
    assert xr.classify_label(label) == state


def test_the_parts_of_a_tally_sum_to_the_rows_it_read(roots):
    write_release(roots[0] / "r", release_id="r", when=iso(17), donors=("112",))
    tally = xr.donor_tallies(xr.resolve("pln", roots=roots))["112"]
    assert (tally.n_called + tally.n_uncalled + tally.n_label_unavailable + tally.n_blank
            + tally.n_no_nucleus_excluded + tally.n_nucleus_not_measured) == tally.n_rows


def test_a_measured_absent_nucleus_excludes_the_cell(roots):
    write_release(roots[0] / "r", release_id="r", when=iso(17), donors=("112",),
                  labels=["Immune_T", "Immune_T"])
    root = roots[0] / "r" / "calls" / "donor_id=112" / "roi=pln_1"
    pq.write_table(pa.table({"cell_id": ["a", "b"], "lineage": ["Immune_T", "Immune_T"],
                             "nucleus_count": [1, 0], "nucleus_measured": [True, True]}),
                   root / "data.parquet")
    tally = xr.donor_tallies(xr.resolve("pln", roots=roots))["112"]
    assert tally.n_called == 1 and tally.n_no_nucleus_excluded == 1
    assert tally.n_nucleus_not_measured == 0


def test_a_never_measured_nucleus_is_not_an_absent_one(roots):
    """Measured 2026-09-17: every cell of 6347, 6579 and 6580 in the pancreas release.

    Counting those as absent nuclei would silently drop three whole donors under a
    rule that never applied to them.
    """
    write_release(roots[0] / "r", release_id="r", when=iso(17), donors=("112",))
    root = roots[0] / "r" / "calls" / "donor_id=112" / "roi=pln_1"
    pq.write_table(pa.table({"cell_id": ["a", "b"], "lineage": ["Immune_T", "Immune_T"],
                             "nucleus_count": [-1, -1], "nucleus_measured": [False, False]}),
                   root / "data.parquet")
    tally = xr.donor_tallies(xr.resolve("pln", roots=roots))["112"]
    assert tally.n_nucleus_not_measured == 2
    assert tally.n_no_nucleus_excluded == 0
    assert tally.n_called == 0
    assert "nucleus_not_measured=2" in tally.reason()


def test_a_release_with_no_nucleus_column_says_so_rather_than_passing(roots):
    write_release(roots[0] / "r", release_id="r", when=iso(17), donors=("112",), nucleus=False)
    tally = xr.donor_tallies(xr.resolve("pln", roots=roots))["112"]
    assert tally.nucleus_rule_applied is False
    assert "UNMEASURABLE" in tally.reason() or "could NOT be applied" in tally.reason()
    assert tally.n_no_nucleus_excluded == 0


# ---------------------------------------------------------------------------
# the builders hold no pin any more
# ---------------------------------------------------------------------------
PINS = ("release_pln_immune", "release_panc_v1", "panc_2026-09-17", "/data/xenium_typing")


def executable_strings(path: Path) -> list[str]:
    """Every string literal the module EXECUTES, with docstrings left out.

    A comment or a docstring may name a retired release -- explaining which pin
    was removed and why is the point.  What must not survive is a literal the code
    reads a path out of.
    """
    tree = ast.parse(path.read_text())
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            body = getattr(node, "body", [])
            if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
                    and isinstance(body[0].value.value, str):
                docstrings.add(id(body[0].value))
    return [n.value for n in ast.walk(tree)
            if isinstance(n, ast.Constant) and isinstance(n.value, str)
            and id(n) not in docstrings]


@pytest.mark.parametrize("name", sorted(BUILDERS))
def test_no_builder_pins_a_release_path_any_more(name):
    """Five pins in five files is the defect; one resolver is the fix."""
    for literal in executable_strings(BUILDERS[name]):
        for pin in PINS:
            assert pin not in literal, f"{name} still pins {pin!r} in the literal {literal!r}"


def test_the_pin_check_would_catch_a_pin_that_came_back(tmp_path):
    """The guard above is only worth having if it fails on the thing it forbids."""
    decoy = tmp_path / "decoy.py"
    decoy.write_text('"""A docstring naming release_panc_v1 is fine."""\n'
                     'ROOT = "data/development/xenium_typing/release_panc_v1"\n')
    literals = executable_strings(decoy)
    assert any("release_panc_v1" in v for v in literals)
    assert not any("A docstring" in v for v in literals)


@pytest.mark.parametrize("name", sorted(BUILDERS))
def test_every_builder_resolves_through_the_one_module(name):
    source = BUILDERS[name].read_text()
    assert "xenium_release" in source, f"{name} does not import the resolver"


# ---------------------------------------------------------------------------
# the real trees
# ---------------------------------------------------------------------------
@pytest.mark.data
def test_both_current_releases_resolve_on_this_checkout():
    for tissue in xr.TISSUES:
        found = xr.resolve(tissue)
        if not found.available:
            pytest.skip(f"no {tissue} release on this machine: {found.reason}")
        assert found.calls_dir.is_dir()
        assert found.observed_donors()
        assert found.generated_at.tzinfo is not None
        assert found.column("label")


@pytest.mark.data
def test_the_grid_layer_reads_the_roots_the_resolver_returns():
    """The grid layer is the reference implementation; its inputs must not move."""
    manifest = Path(os.environ.get("ISLET_LAYERS_DIR", "/data/islet_layers")) \
        / "grid_composition" / "manifest.json"
    if not manifest.is_file():
        pytest.skip(f"no grid_composition layer at {manifest}")
    built = {row["tissue"] + "/" + row["modality"]: row["root"]
             for row in json.loads(manifest.read_text())["grid"]}
    for tissue, key in (("panc", "pancreas/RNA"), ("pln", "pln/RNA")):
        found = xr.resolve(tissue)
        if not found.available:
            pytest.skip(f"no {tissue} release on this machine")
        assert str(found.root) == built[key], (
            f"the resolver returns {found.root} where the built layer read {built[key]}; "
            "the grid layer's numbers were computed on the latter"
        )


@pytest.mark.data
def test_a_tally_agrees_with_the_releases_own_book_keeping():
    """Internal coherence only: the release counted its own file, and so did we.

    This is not a check against any prior result, biological or otherwise. It
    compares two counts of the SAME rows -- the producer's, recorded in
    release.json, and this module's, read back off the parquet.
    """
    found = xr.resolve("panc")
    if not found.available:
        pytest.skip(f"no panc release on this machine: {found.reason}")
    declared_total = found.metadata.get("n_cells_total")
    declared_no_nucleus = found.metadata.get("n_cells_no_nucleus_excluded")
    if declared_total is None:
        pytest.skip("this release declares no cell totals")
    tallies = xr.donor_tallies(found)
    assert sum(t.n_rows for t in tallies.values()) == declared_total
    if declared_no_nucleus is not None:
        assert sum(t.n_no_nucleus_excluded for t in tallies.values()) == declared_no_nucleus
    for donor, tally in tallies.items():
        assert (tally.n_called + tally.n_uncalled + tally.n_label_unavailable
                + tally.n_blank + tally.n_no_nucleus_excluded
                + tally.n_nucleus_not_measured) == tally.n_rows, donor
