"""The code tier of `scripts/retire/`: what it lists, what it refuses, and what it cannot touch.

Two things this suite is deliberately built around.

**The exclusion is a property of the code, not of a paragraph.** `merge/hpg-sameslide` is being
developed on right now, and the guarantee that a retirement run cannot take one of its files is
only worth anything if something checks it. `dispositions.KEEP_BESIDE_CODE` and
`KEEP_BESIDE_CODE_TREES` are that guarantee and these tests are what hold them in place -- a path
added to a retirement tier that happens to sit under `scripts/sameslide/` fails here rather than
in a quarantine directory.

**Resolution is by path.** Two files in this repo are called `aggregate_islets.py` and three are
called `pins.py`. A retirement check that matched basenames would have taken
`scripts/senior/aggregate_islets.py` off the list (because the typing lane's file of the same name
is a live node in `build_layers.py`'s DAG) and put `core/phenocycler/scripts/broadtype/` on it
(because `scripts/report/bakeoff.py` imports *a* `pins`). Both are asserted below.
"""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[4]
RETIRE = REPO / "scripts" / "retire"
sys.path.insert(0, str(RETIRE))

import build_manifest as BM  # noqa: E402
import dispositions as D  # noqa: E402
import retire as R  # noqa: E402

CODE_ENTRIES = D.CODE_TIER1 + D.CODE_TIER2
ALL_ENTRIES = CODE_ENTRIES + D.ARTIFACT_TIER2


@pytest.fixture(scope="module")
def manifest(tmp_path_factory) -> dict:
    """The code manifest, built once. Written to tmp_path -- never under data/."""
    out = tmp_path_factory.mktemp("retire") / "MANIFEST.csv"
    subprocess.run([sys.executable, str(RETIRE / "build_manifest.py"),
                    "--kind", "code", "--out", str(out), "--quiet"],
                   check=True, cwd=REPO)
    with out.open() as fh:
        rows = list(csv.DictReader(fh))
    return {"path": out, "rows": rows, "by_path": {r["path"]: r for r in rows}}


# ----------------------------------------------------------- the exclusion ----
def test_the_hpg_trees_the_plan_names_are_all_guarded():
    """sameslide/, xenium_celltype/ (with its hpg/) and run_field_refine.py, by name."""
    guarded = {p for p, _ in D.KEEP_BESIDE_CODE} | set(D.KEEP_BESIDE_CODE_TREES)
    for required in (f"{D.REPO}/core/phenocycler/scripts/sameslide",
                     f"{D.REPO}/core/phenocycler/scripts/xenium_celltype",
                     f"{D.REPO}/core/phenocycler/scripts/restain/run_field_refine.py"):
        assert required in guarded, f"{required} is not on the keep-beside list"
    assert f"{D.REPO}/core/phenocycler/scripts/xenium_celltype" in D.KEEP_BESIDE_CODE_TREES, \
        "the hpg/ subdir is only covered if xenium_celltype is guarded as a TREE"


def test_no_retirement_entry_is_a_protected_path():
    protected = {p for p, _ in D.KEEP_BESIDE_CODE}
    for entry in ALL_ENTRIES:
        assert entry["path"] not in protected, entry["path"]
        assert R.under_protected_tree(entry["path"]) is None, \
            f"{entry['path']} sits under a tree merge/hpg-sameslide is developing on"


def test_retire_refuses_a_protected_path_even_if_it_reaches_a_tier(tmp_path):
    """The belt-and-braces case: a path put on tier 1 by mistake, with a PASSing check."""
    victim = f"{D.REPO}/core/phenocycler/scripts/sameslide/lift_calls.py"
    assert R.under_protected_tree(victim), "precondition: the tree is guarded"
    manifest = tmp_path / "M.csv"
    with manifest.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=BM.FIELDS)
        w.writeheader()
        w.writerow(BM.row(1, victim, kind="code"))
    out = subprocess.run([sys.executable, str(RETIRE / "retire.py"), "--kind", "code",
                          "--tier", "1", "--manifest", str(manifest)],
                         capture_output=True, text=True, cwd=REPO)
    assert out.returncode == 0, out.stderr
    assert "SKIP" in out.stdout and "lift_calls.py" in out.stdout
    assert "0 entries" in out.stdout
    assert D.HPG_BRANCH in out.stdout


def test_code_is_never_deleted_by_this_script():
    out = subprocess.run([sys.executable, str(RETIRE / "retire.py"), "--kind", "code",
                          "--tier", "1", "--mode", "delete"],
                         capture_output=True, text=True, cwd=REPO)
    assert out.returncode != 0
    assert "refusing --mode delete" in (out.stderr + out.stdout)


# --------------------------------------------------- resolution is by path ----
def test_a_basename_shared_with_a_live_file_does_not_block_the_orphan(manifest):
    """`scripts/senior/aggregate_islets.py` is an orphan; the typing lane's file of the same
    name is a live node in build_layers.py. Only the first is listed, and it passes."""
    senior = f"{D.REPO}/scripts/senior/aggregate_islets.py"
    typing = f"{D.REPO}/core/phenocycler/scripts/typing/aggregate_islets.py"
    assert senior in manifest["by_path"], "the senior orphan should be listed"
    assert typing not in manifest["by_path"], "the live typing producer must never be listed"
    assert manifest["by_path"][senior]["verified"] == "PASS", \
        manifest["by_path"][senior]["verify_detail"]
    assert Path(typing).exists(), "precondition: the live file of the same name is still there"


def test_an_import_that_resolves_to_another_lane_is_not_a_live_reference():
    """`scripts/report/{bakeoff,identity_io}.py` import `pins`, meaning their OWN pins.py."""
    report = Path(D.REPO) / "core/phenocycler/scripts/report"
    broadtype = Path(D.REPO) / "core/phenocycler/scripts/broadtype"
    assert (report / "pins.py").exists() and (broadtype / "pins.py").exists()
    assert not BM._import_resolves_here(report / "bakeoff.py", "pins", broadtype)
    assert BM._import_resolves_here(report / "bakeoff.py", "pins", report)


def test_an_import_reached_through_sys_path_is_a_live_reference():
    """kronos2_doc/run.py puts scripts/broadtype on sys.path, so its `import pins` is broadtype's."""
    run = Path(D.REPO) / "core/phenocycler/scripts/kronos2_doc/run.py"
    broadtype = Path(D.REPO) / "core/phenocycler/scripts/broadtype"
    assert run.exists(), "precondition"
    assert BM._import_resolves_here(run, "pins", broadtype)


def test_a_mention_inside_a_python_docstring_is_not_a_call(tmp_path):
    """Otherwise this lane's own rebuild_adata_ins_root.py, whose opening paragraph explains why
    two legacy producers are not its ancestors, would read as a live caller of both."""
    f = tmp_path / "m.py"
    f.write_text('"""Explains scripts/reaggregate_islets.py at length."""\nimport os\n')
    assert 1 in BM._py_doc_lines(f)
    assert BM._py_imports(f) == frozenset({"os"})


def test_imports_are_read_from_the_ast_not_the_text(tmp_path):
    f = tmp_path / "m.py"
    f.write_text('MSG = "import substrate"\n# import positivity\nimport json\n')
    assert BM._py_imports(f) == frozenset({"json"})


# ------------------------------------------------------- the table itself ----
def test_every_listed_code_path_still_exists():
    """A table that has drifted from disk must be corrected, not acted on."""
    for entry in CODE_ENTRIES:
        assert Path(entry["path"]).exists(), f"{entry['path']} is on a code tier but is gone"


def test_no_path_is_listed_twice():
    seen = [e["path"] for e in ALL_ENTRIES] + [p for p, _ in D.CODE_TIER3] \
        + [p for p, _ in D.CODE_TIER4] + [p for p, _ in D.KEEP_BESIDE_CODE]
    assert len(seen) == len(set(seen)), \
        sorted({p for p in seen if seen.count(p) > 1})


def test_code_and_data_tiers_do_not_overlap():
    data = {e["path"] for e in D.TIER1 + D.TIER2}
    assert not data & {e["path"] for e in ALL_ENTRIES}


def test_every_blocker_is_itself_on_a_code_tier():
    listed = {e["path"] for e in ALL_ENTRIES}
    for entry in ALL_ENTRIES:
        for blocker in entry["blocked_by"]:
            assert blocker in listed, f"{entry['path']} is blocked by unlisted {blocker}"


def test_a_blocked_entry_is_held_when_its_blocker_stays(manifest):
    """fixed_islet_aggregation.py is imported by two files on the list; alone, it is a broken
    import, so the tier-1 pass holds it while reaggregate_islets.py (tier 2) is still there."""
    out = subprocess.run([sys.executable, str(RETIRE / "retire.py"), "--kind", "code",
                          "--tier", "1", "--manifest", str(manifest["path"])],
                         capture_output=True, text=True, cwd=REPO)
    assert out.returncode == 0, out.stderr
    held = [ln for ln in out.stdout.splitlines()
            if "fixed_islet_aggregation.py" in ln and "SKIP" in ln]
    assert held, out.stdout
    assert "blocked_by" in held[0]


def test_a_dry_run_moves_nothing(manifest):
    before = {e["path"]: Path(e["path"]).exists() for e in ALL_ENTRIES}
    subprocess.run([sys.executable, str(RETIRE / "retire.py"), "--kind", "code", "--tier", "1",
                    "--manifest", str(manifest["path"])], check=True, capture_output=True, cwd=REPO)
    subprocess.run([sys.executable, str(RETIRE / "retire.py"), "--kind", "code", "--tier", "2",
                    "--manifest", str(manifest["path"])], check=True, capture_output=True, cwd=REPO)
    assert {e["path"]: Path(e["path"]).exists() for e in ALL_ENTRIES} == before


def test_the_overturned_entries_carry_their_measurement():
    overturned = dict(D.CODE_TIER3)
    for path in (f"{D.REPO}/core/phenocycler/scripts/broadtype",
                 f"{D.REPO}/core/phenocycler/scripts/dct",
                 f"{D.REPO}/core/phenocycler/scripts/typing/call.py"):
        assert path in overturned, f"{path} must be recorded as overturned, with why"
        assert "import" in overturned[path].lower()


# -------------------------------------- adata_ins_root.h5ad has one producer ----
TRAJECTORY = REPO / "core/phenocycler/scripts/typing/trajectory_islet.py"
REBUILD = REPO / "scripts/rebuild_adata_ins_root.py"


def test_trajectory_islet_does_not_target_the_deployed_h5ad():
    source = TRAJECTORY.read_text()
    assert 'out_dir / "adata_ins_root.h5ad"' not in source
    assert 'adata_ins_root.rejected.h5ad"' not in source
    assert "REFUSED_H5AD_STEMS" in source


def test_trajectory_islet_refuses_that_name_on_the_command_line():
    out = subprocess.run([sys.executable, str(TRAJECTORY), "--h5ad-name", "adata_ins_root.h5ad"],
                         capture_output=True, text=True, cwd=REPO)
    assert out.returncode != 0
    assert "rebuild_adata_ins_root.py" in (out.stderr + out.stdout)


def test_trajectory_islet_still_writes_the_clustering_before_the_gate():
    """islet_leiden.parquet is what build_layers.py declares as this node's output, and the
    pseudotime gate does not judge it -- so it must stay on the early side of the gate."""
    lines = TRAJECTORY.read_text().splitlines()
    leiden = next(i for i, ln in enumerate(lines, 1) if "to_parquet(out_dir / args.leiden_name" in ln)
    gate = next(i for i, ln in enumerate(lines, 1) if "if not all(gates.values()):" in ln)
    assert leiden < gate, f"leiden write at {leiden} must precede the gate at {gate}"


def test_rebuild_is_the_single_producer():
    source = REBUILD.read_text()
    assert 'OUT_NAME = "adata_ins_root.h5ad"' in source
    assert "SINGLE PRODUCER" in source


def test_the_frame_is_resolved_from_an_explicit_ordered_list_not_an_mtime():
    """`build_distance_profiles.py` picks the newest release.json mtime and a `touch` reorders it.
    The frame this writes the deployed file from must not be choosable that way."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("rebuild_ins", REBUILD)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert len(module.FRAME_CANDIDATES) >= 2
    assert all(isinstance(c, Path) for c in module.FRAME_CANDIDATES)
    import inspect
    chooser = inspect.getsource(module.resolve_frame)
    for forbidden in ("mtime", "glob", "sorted", "max("):
        assert forbidden not in chooser, f"resolve_frame() chooses with {forbidden!r}"


def test_only_one_file_in_the_tree_writes_the_deployed_path():
    """A second producer is the state root CLAUDE.md forbids. The legacy pair target
    <repo>/data/adata_ins_root.h5ad, a different path in an unrunnable lane, and are on the
    retirement list; nothing else may name the app root's copy as an output."""
    writers = []
    for path in REPO.rglob("*.py"):
        if any(part in {".venv", "external", "node_modules", "data", "tests"} for part in path.parts):
            continue
        text = path.read_text(errors="replace")
        if "APP / OUT_NAME" in text or 'APP / "adata_ins_root.h5ad"' in text:
            writers.append(str(path.relative_to(REPO)))
    assert writers == ["scripts/rebuild_adata_ins_root.py"], writers


# --------------------------------------------------------------- ie_holds ----
def test_ie_holds_is_listed_with_its_prerequisites():
    entry = next(e for e in D.CODE_TIER1 if e["path"].endswith("ie_holds.R"))
    assert "test-guards.R" in entry["prereq"] and "test-no-domain-pixel.R" in entry["prereq"]
    assert "NAMESPACE" in entry["prereq"]


def test_ie_holds_still_has_no_caller():
    """The decision was 'retire', not 'wire it up'. If a caller appears, that decision is stale."""
    callers = []
    for path in list(REPO.rglob("*.R")) + list(REPO.rglob("*.py")):
        if any(part in {".venv", "external", "node_modules", "data", "tests", "retire"}
               for part in path.parts):
            continue          # the disposition table and this suite name the two functions
        if path.name == "ie_holds.R":
            continue
        text = path.read_text(errors="replace")
        for fn in ("ie_assert_hold_lifted(", "ie_hold_is_lifted("):
            if fn in text:
                callers.append(f"{path.relative_to(REPO)}: {fn}")
    assert callers == [], callers
