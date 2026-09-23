"""The code tier of `scripts/retire/`: what it lists, what it refuses, and what it cannot touch.

Two things this suite is deliberately built around.

**The guard is a property of the code, not of a paragraph.** A module something live imports
must not be movable by a retirement run, and that guarantee is only worth anything if something
checks it. `dispositions.CODE_LIBRARY` -- every path the 2026-09-22 reference walk measured as
reached by live code, importers named -- is that guarantee and these tests are what hold it in
place: a path put on a retirement tier that is a declared library, or sits under a library
directory, fails here rather than in a quarantine directory. `CODE_CANDIDATE` is the other half
of that walk: paths with NO live reference, recorded for the maintainer and neither refused nor
acted on. (Both replaced the `merge/hpg-sameslide` branch protection, which lapsed at the merge.)

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


# ------------------------------------------------------ the declared libraries ----
LIBRARY = {p for p, _ in D.CODE_LIBRARY}
CANDIDATE = {p for p, _ in D.CODE_CANDIDATE}
SCRIPTS = f"{D.REPO}/core/phenocycler/scripts"

#: What the 2026-09-22 walk measured LIVE, by name. The list is the guard; this pins it.
MEASURED_LIBRARIES = (
    *(f"{SCRIPTS}/celltype/{m}.py" for m in ("panel", "rules", "score", "usability",
                                              "run_step3", "run_islet")),
    *(f"{SCRIPTS}/restain/{m}.py" for m in ("align", "registration_inputs", "run_align",
                                             "run_intensity_refine", "run_register")),
    f"{D.REPO}/core/phenocycler/phenocycler/integration/transform.py",
    *(f"{SCRIPTS}/xenium_celltype/{m}.py" for m in ("xenio", "nullgate", "build_phenotyped_base")),
    f"{SCRIPTS}/xenium_celltype/hpg/nullgate_fast.py",
    f"{SCRIPTS}/xenium_celltype/tests",
    f"{SCRIPTS}/xenium_celltype/hpg/tests",
    f"{D.REPO}/analysis/xenium/scripts/tissue_markers/pancreas.py",
    f"{D.REPO}/core/phenocycler/tests/current/test_no_duplicate_tracked_data.py",
    f"{D.REPO}/docs/XENIUM_PIPELINE.md",
)


def test_every_library_the_walk_measured_is_declared_with_its_importers():
    for required in MEASURED_LIBRARIES:
        assert required in LIBRARY, f"{required} is not on CODE_LIBRARY"
    for path, why in D.CODE_LIBRARY:
        assert ":" in why, f"{path}: a library row names its importers as path:line"


def test_every_declared_library_exists():
    """The branch that once justified tolerating absence has merged; a missing library is a
    drifted table, and the manifest check must FAIL on it rather than guard a ghost."""
    for path in LIBRARY:
        assert Path(path).exists(), f"{path} is declared a library but is gone"
    ok, detail = BM.check_library(Path(f"{D.REPO}/no/such/library.py"), None)
    assert not ok and "MISSING" in detail


def test_no_retirement_entry_is_a_library_or_under_one():
    for entry in ALL_ENTRIES:
        assert entry["path"] not in LIBRARY, entry["path"]
        assert R.under_protected_tree(entry["path"]) is None, \
            f"{entry['path']} sits under a declared library directory"


def test_a_library_directory_protects_what_is_under_it():
    inside = f"{SCRIPTS}/xenium_celltype/tests/test_nullgate.py"
    assert R.under_protected_tree(inside) == f"{SCRIPTS}/xenium_celltype/tests"
    assert R.under_protected_tree(f"{SCRIPTS}/xenium_celltype/hpg/sync_code.sh") is None, \
        "a file beside a library directory, not under it, is not protected by it"


@pytest.mark.parametrize("victim", [
    f"{SCRIPTS}/celltype/panel.py",                       # a library file, refused by name
    f"{SCRIPTS}/xenium_celltype/tests/test_nullgate.py",  # under a library directory
])
def test_retire_refuses_a_declared_library_even_if_it_reaches_a_tier(tmp_path, victim):
    """The belt-and-braces case: a library put on tier 1 by mistake, with a PASSing check."""
    manifest = tmp_path / "M.csv"
    with manifest.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=BM.FIELDS)
        w.writeheader()
        w.writerow(BM.row(1, victim, kind="code"))
    out = subprocess.run([sys.executable, str(RETIRE / "retire.py"), "--kind", "code",
                          "--tier", "1", "--manifest", str(manifest)],
                         capture_output=True, text=True, cwd=REPO)
    assert out.returncode == 0, out.stderr
    assert "SKIP" in out.stdout and Path(victim).name in out.stdout
    assert "0 entries" in out.stdout
    assert "declared library" in out.stdout


# ------------------------------------------------------------ the candidates ----
def test_the_candidates_the_walk_found_were_ruled_on_and_carry_their_evidence():
    """The 2026-09-22 walk recorded 23 candidates; the maintainer ruled on them 2026-09-23.
    The same-slide lane is a library WITH the ruling on its row; the rest are tier 2 with the
    walk's evidence carried unchanged; nothing is a candidate any more."""
    libraries = dict(D.CODE_LIBRARY)
    for kept in (f"{SCRIPTS}/sameslide",
                 f"{SCRIPTS}/restain/run_field_refine.py",
                 f"{SCRIPTS}/celltype/run_sameslide_protein.py"):
        assert kept in libraries, f"{kept} must be a declared library (ruled kept 2026-09-23)"
        assert "RULED KEPT 2026-09-23" in libraries[kept], f"{kept}: the ruling is on its row"
        assert "2026-09-22" in libraries[kept], f"{kept}: the walk's evidence stays on its row"
    tier2 = {e["path"]: e for e in D.CODE_TIER2}
    for moved in (f"{SCRIPTS}/typing/compare_clustering.py",
                  f"{SCRIPTS}/typing/call.py",
                  f"{SCRIPTS}/xenium_celltype/eval_detection.py",
                  f"{SCRIPTS}/xenium_celltype/hpg/fetch_results.sh",
                  f"{SCRIPTS}/xenium_celltype/hpg/sync_code.sh",
                  f"{SCRIPTS}/xenium_celltype/hpg/slurm"):
        assert moved in tier2, f"{moved} must be on CODE_TIER2 (ruled 2026-09-23)"
        assert tier2[moved]["evidence"], f"{moved}: the walk's evidence is carried"
    assert tier2[f"{SCRIPTS}/typing/call.py"]["blocked_by"] == (
        f"{SCRIPTS}/typing/compare_clustering.py",), "call.py moves with its one importer"
    assert tier2[f"{SCRIPTS}/xenium_celltype/hpg/stage_prep.py"]["blocked_by"] == (
        f"{SCRIPTS}/xenium_celltype/hpg/stage_call.py",)
    assert D.CODE_CANDIDATE == [], "every 2026-09-22 candidate was ruled on 2026-09-23"
    for path, why in D.CODE_CANDIDATE:
        assert Path(path).exists(), f"{path} is a candidate but is gone"
        assert why, f"{path}: a candidate carries the walk's evidence"


def test_a_candidate_is_on_no_tier_and_in_no_refusal_set():
    tiered = {e["path"] for e in ALL_ENTRIES} | {p for p, _ in D.CODE_TIER3 + D.CODE_TIER4}
    for path in CANDIDATE:
        assert path not in tiered, f"{path} is both a candidate and on a tier"
        assert path not in LIBRARY, f"{path} is both a candidate and a library"
        assert path not in R.REFUSE, f"{path}: a candidate is not refused, it is recorded"
        assert R.under_protected_tree(path) is None, f"{path} sits under a library directory"


def test_a_candidate_is_neither_refused_nor_acted_on(manifest):
    """The manifest carries every candidate as tier CANDIDATE with the orphan walk run; a
    retire run over tiers 1 and 2 never mentions one, in either the acting or the skip list."""
    for path in CANDIDATE:
        row = manifest["by_path"][path]
        assert row["tier"] == "CANDIDATE" and row["check"] == "candidate_orphan", row
        assert row["verified"] in ("PASS", "FAIL"), "a candidate row carries a measurement"
    # context.py was the worked example of "a reference from another candidate is internal to
    # the question": its only importer was eval_detection.py, a candidate too. Both were ruled
    # to tier 2 and moved together on 2026-09-23, which the log records; the manifest now
    # reports the path absent rather than measuring a file that is not there.
    context = f"{SCRIPTS}/xenium_celltype/context.py"
    assert context in _moved_by_the_ledger(), "context.py moved with eval_detection.py"
    assert manifest["by_path"][context]["verify_detail"] == "MISSING"
    for tier in ("1", "2"):
        out = subprocess.run([sys.executable, str(RETIRE / "retire.py"), "--kind", "code",
                              "--tier", tier, "--manifest", str(manifest["path"])],
                             capture_output=True, text=True, cwd=REPO)
        assert out.returncode == 0, out.stderr
        for path in CANDIDATE:
            assert path.replace(f"{D.REPO}/", "") not in out.stdout, \
                f"tier {tier} run names the candidate {path}"


def test_every_library_row_is_in_the_manifest_and_verified(manifest):
    for path in LIBRARY:
        row = manifest["by_path"][path]
        assert row["tier"] == "LIBRARY" and row["check"] == "library", row
        assert row["verified"] == "PASS", row["verify_detail"]


# ---------------------------------------------- the alternative typers, declared ----
def test_the_alternative_typers_are_declared_and_not_movable():
    """Each runnable alternative typer the producers census named is a decision with its
    measured live use, in the tier the census proposed. None is a fifth producer."""
    libraries = dict(D.CODE_LIBRARY)
    for path, importer in ((f"{D.REPO}/core/phenocycler/phenocycler/simultaneous_typing.py",
                            "production_typing.py:21"),
                           (f"{SCRIPTS}/pln/cluster_pln.py", "build_pln_immune_pheno.py:166"),
                           (f"{SCRIPTS}/typing/refine.py", "vocab.py:206")):
        assert path in libraries, f"{path} must be a declared library"
        assert importer in libraries[path], f"{path}: the live use {importer} must be named"
    tier2 = {e["path"]: e for e in D.CODE_TIER2}
    for name in ("release", "nullgate_typer", "baseline_analyzer", "scanvi_replicates"):
        entry = tier2[f"{SCRIPTS}/xenium_typing/{name}.py"]
        assert "config.ini:153" in entry["prereq"], f"{name}: the Step 5 prerequisite is named"
    hier = tier2[f"{D.REPO}/analysis/xenium/scripts/hierarchical_phenotype.py"]
    assert "nullgate_typer.py" in hier["prereq"] and "run_pipeline.sh" in hier["prereq"]
    assert f"{SCRIPTS}/xenium_typing/nullgate_typer.py" in hier["blocked_by"]
    overturned = dict(D.CODE_TIER3)
    assert "xen_identity.py" in overturned[f"{SCRIPTS}/xenium_typing/joint_model.py"]


def test_the_apparatus_moved_once_its_binding_was_blank(manifest):
    """Until 2026-09-23 core/phenocycler/config.ini:153 bound the release directory the X1/X3
    arms and release.py write, and the manifest FAILED all four by name on that binding. Step 5
    blanked it; the check then passed on all four (baseline_analyzer once its lane test was
    removed with it) and the maintainer's ruling executed the tier the same day. The log
    records the move; the key is still read blank; a value there would be the binding back."""
    moved = _moved_by_the_ledger()
    for name in ("release", "nullgate_typer", "scanvi_replicates", "baseline_analyzer"):
        path = f"{SCRIPTS}/xenium_typing/{name}.py"
        row = manifest["by_path"][path]
        assert row["check"] == "code_orphan_and_config_blank", row["check"]
        assert path in moved, f"{name}: quarantined 2026-09-23, the log must say so"
        assert not Path(path).exists() and row["verify_detail"] == "MISSING", row
    assert BM._config_binding(BM.CONFIG_INI, "xenium_identity_release") == (153, "")
    assert BM._config_binding(BM.CONFIG_INI, "phenocycler_identity_release") == (152, "")
    row = manifest["by_path"][f"{D.REPO}/analysis/xenium/scripts/hierarchical_phenotype.py"]
    assert row["verified"] == "FAIL" and "run_pipeline.sh" in row["verify_detail"], row


def test_the_config_binding_check_reads_the_key_by_name(tmp_path):
    ini = tmp_path / "config.ini"
    ini.write_text("[integration]\nphenocycler_identity_release = ../../x\n"
                   "xenium_identity_release = \n")
    assert BM._config_binding(ini, "phenocycler_identity_release") == (2, "../../x")
    assert BM._config_binding(ini, "xenium_identity_release") == (3, "")
    assert BM._config_binding(ini, "no_such_key") is None
    live = BM._config_binding(BM.CONFIG_INI, "xenium_identity_release")
    assert live is not None and live[0] == 153, "the prerequisite cites config.ini:153"


# ------------------------------------------------------------ the log migration ----
V1_HEADER = "retired_utc,tier,path,bytes,action,superseded_by,evidence"
V1_ROWS = [
    ("2026-09-16T14:11:24+00:00", "1", "/data/tmp_b2r", "79883166872",
     "/data/ARCHIVE_2026-09-15_quarantine/tmp_b2r", "/data/islet_ome_tiff/6593.ome.tiff",
     "Single orphaned mktemp dir; the conversion completed, both 2026-09-14 11:06."),
    ("2026-09-16T14:11:24+00:00", "2", "/x/data/cells", "1", "/x/quarantine/data/cells",
     "/x/runs/{cells,geometry_qc}", "donor_id=6457 (728,698 rows) exist here, \"quoted\"."),
]


def _v1_log(tmp_path) -> Path:
    log = tmp_path / "RETIRED.csv"
    with log.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(V1_HEADER.split(","))
        w.writerows(V1_ROWS)
    return log


def _code_row(path="/x/scripts/old.py"):
    return dict(retired_utc="2026-09-22T00:00:00+00:00", kind="code", tier=1, path=path,
                bytes="10", action="/x/quarantine/scripts/old.py", superseded_by="-",
                evidence="orphan")


def test_the_first_code_append_migrates_a_v1_log_once(tmp_path):
    """One log stays: header gains `kind` second, old rows are stamped kind=data, the new row
    lands under the right columns, and a second append does not migrate again."""
    log = _v1_log(tmp_path)
    R.append_log(log, [_code_row()])
    with log.open(newline="") as fh:
        reader = csv.DictReader(fh)
        assert reader.fieldnames == R.LOG_FIELDS
        rows = list(reader)
    assert len(rows) == len(V1_ROWS) + 1
    for old, migrated in zip(V1_ROWS, rows, strict=False):
        assert migrated["kind"] == "data"
        assert [migrated[k] for k in R.LOG_FIELDS_V1] == list(old), "no column shifted"
    assert rows[-1]["kind"] == "code" and rows[-1]["path"] == "/x/scripts/old.py"
    assert rows[-1]["evidence"] == "orphan"
    assert not R.migrate_log(log), "already migrated"
    R.append_log(log, [_code_row("/x/scripts/other.py")])
    with log.open(newline="") as fh:
        again = list(csv.DictReader(fh))
    assert len(again) == len(V1_ROWS) + 2 and again[0]["kind"] == "data"
    assert not (tmp_path / "RETIRED.csv.migrating").exists()


def test_a_fresh_log_is_written_in_the_current_shape(tmp_path):
    log = tmp_path / "RETIRED.csv"
    assert not R.migrate_log(log)
    R.append_log(log, [_code_row()])
    with log.open(newline="") as fh:
        reader = csv.DictReader(fh)
        assert reader.fieldnames == R.LOG_FIELDS
        assert [r["kind"] for r in reader] == ["code"]


def test_the_real_log_is_never_touched_by_this_suite():
    """The migration runs against the live log only inside `retire.py --execute`."""
    assert R.LOG == Path(D.REPO) / "data" / "ARCHIVE_2026-09-15" / "RETIRED.csv"


def test_code_is_never_deleted_by_this_script():
    out = subprocess.run([sys.executable, str(RETIRE / "retire.py"), "--kind", "code",
                          "--tier", "1", "--mode", "delete"],
                         capture_output=True, text=True, cwd=REPO)
    assert out.returncode != 0
    assert "refusing --mode delete" in (out.stderr + out.stdout)



def _moved_by_the_ledger() -> set[str]:
    """Paths the retirement log records as quarantined (kind=code). Read from the real log,
    never written: an entry that is absent from disk AND recorded here is a retirement that
    happened, not a table that drifted."""
    if not R.LOG.exists():
        return set()
    with R.LOG.open(newline="") as handle:
        return {row["path"] for row in csv.DictReader(handle) if row.get("kind") == "code"}


# --------------------------------------------------- resolution is by path ----
def test_a_basename_shared_with_a_live_file_does_not_block_the_orphan(manifest):
    """`scripts/senior/aggregate_islets.py` is an orphan; the typing lane's file of the same
    name is a live node in build_layers.py. Only the first is listed, and it passes."""
    senior = f"{D.REPO}/scripts/senior/aggregate_islets.py"
    typing = f"{D.REPO}/core/phenocycler/scripts/typing/aggregate_islets.py"
    assert senior in manifest["by_path"], "the senior orphan should be listed"
    assert typing not in manifest["by_path"], "the live typing producer must never be listed"
    if senior in _moved_by_the_ledger():
        assert not Path(senior).exists(), "recorded as quarantined, so it must be gone"
    else:
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
    # A fake name: a real retiring path here would count as a LIVE reference to that entry.
    f.write_text('"""Explains scripts/no_such_producer.py at length."""\nimport os\n')
    assert 1 in BM._py_doc_lines(f)
    assert BM._py_imports(f) == frozenset({"os"})


def test_imports_are_read_from_the_ast_not_the_text(tmp_path):
    f = tmp_path / "m.py"
    f.write_text('MSG = "import substrate"\n# import positivity\nimport json\n')
    assert BM._py_imports(f) == frozenset({"json"})


# ------------------------------------------------------- the table itself ----
def test_every_listed_code_path_still_exists():
    """A table that has drifted from disk must be corrected, not acted on."""
    moved = _moved_by_the_ledger()
    for entry in CODE_ENTRIES:
        if entry["path"] in moved:
            assert not Path(entry["path"]).exists(), \
                f"{entry['path']} is recorded as quarantined but is still on disk"
            continue
        assert Path(entry["path"]).exists(), f"{entry['path']} is on a code tier but is gone"


def test_no_path_is_listed_twice():
    seen = [e["path"] for e in ALL_ENTRIES] + [p for p, _ in D.CODE_TIER3] \
        + [p for p, _ in D.CODE_TIER4] + [p for p, _ in D.CODE_LIBRARY] \
        + [p for p, _ in D.CODE_CANDIDATE]
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
    if f"{D.REPO}/scripts/senior/fixed_islet_aggregation.py" in _moved_by_the_ledger():
        assert "already gone" in held[0], held[0]
    else:
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
                 f"{D.REPO}/core/phenocycler/scripts/dct"):
        assert path in overturned, f"{path} must be recorded as overturned, with why"
        assert "import" in overturned[path].lower()
    # call.py was overturned while compare_clustering.py imported it; ruled retired 2026-09-23.
    assert f"{D.REPO}/core/phenocycler/scripts/typing/call.py" not in overturned


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
        if any(part in {".venv", "external", "node_modules", "data", "tests"}
               for part in path.relative_to(REPO).parts):
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
               for part in path.relative_to(REPO).parts):
            continue          # the disposition table and this suite name the two functions
        if path.name == "ie_holds.R":
            continue
        text = path.read_text(errors="replace")
        for fn in ("ie_assert_hold_lifted(", "ie_hold_is_lifted("):
            if fn in text:
                callers.append(f"{path.relative_to(REPO)}: {fn}")
    assert callers == [], callers
