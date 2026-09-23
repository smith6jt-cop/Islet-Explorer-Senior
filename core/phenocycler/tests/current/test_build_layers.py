"""The build driver: its declared graph, its ordering guard, and its staleness gate.

    .venv/bin/python -m pytest core/phenocycler/tests/current/test_build_layers.py -q

Nothing here asserts a biological number, and nothing here reruns a builder.  The mtime tests run
against files this test creates in tmp_path; the tests that read the real repo are marked `data`
and assert structure -- which nodes exist, which paths they name -- never what is in an artifact.
"""
from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[4]
DRIVER = REPO / "scripts/senior/build_layers.py"
WRAPPER = REPO / "scripts/senior/run_pipeline.sh"


def _load():
    spec = importlib.util.spec_from_file_location("build_layers", DRIVER)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    # Registered before exec so the frozen dataclasses in it can resolve their own module.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


bl = _load()


def _touch(path: Path, when: float) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("x")
    os.utime(path, (when, when))
    return path


# ---------------------------------------------------------------------------------------------
# The declared graph
# ---------------------------------------------------------------------------------------------

def test_the_graph_builds_and_is_acyclic(tmp_path):
    nodes = bl.build_graph(REPO, tmp_path / "layers")
    order = bl.topological_order(nodes)
    assert set(order) == set(nodes)
    seen: set[str] = set()
    dependencies = bl.edges(nodes)
    for name in order:
        assert dependencies[name] <= seen, f"{name} runs before {dependencies[name] - seen}"
        seen.add(name)


def test_all_nine_layers_are_nodes(tmp_path):
    nodes = bl.build_graph(REPO, tmp_path / "layers")
    missing = sorted(set(bl.LAYER_NODES) - set(nodes))
    assert not missing, f"the driver does not cover {missing}"
    assert len(bl.LAYER_NODES) == 9


def test_every_declared_builder_exists_on_disk(tmp_path):
    """A driver that names a script that is not there is prose, not a driver."""
    nodes = bl.build_graph(REPO, tmp_path / "layers")
    for node in nodes.values():
        script = Path(node.command[1])
        assert script.is_file(), f"{node.name} names a builder that does not exist: {script}"


def test_every_node_runs_the_repo_venv_python(tmp_path):
    """`python` on PATH is not this checkout's environment; /usr/bin/python3 has no pyarrow."""
    nodes = bl.build_graph(REPO, tmp_path / "layers")
    for node in nodes.values():
        assert node.command[0] == str(bl.VENV_PYTHON), (
            f"{node.name} does not use the repo venv: {node.command[0]}")


def test_the_layers_dir_is_honoured_rather_than_hardcoded(tmp_path):
    nodes = bl.build_graph(REPO, tmp_path / "elsewhere")
    outputs = nodes["sections"].outputs
    assert outputs == (str(tmp_path / "elsewhere/sections/data.parquet"),)


def test_no_node_declares_a_relative_path(tmp_path):
    """Every path is absolute after expansion, so a run's cwd cannot change what it means."""
    nodes = bl.build_graph(REPO, tmp_path / "layers")
    for node in nodes.values():
        for path in node.inputs + node.outputs + node.soft_inputs:
            assert Path(path).is_absolute(), f"{node.name}: {path}"


# ---------------------------------------------------------------------------------------------
# The islets.parquet ordering hazard
# ---------------------------------------------------------------------------------------------

def test_parquet_views_is_the_one_writer_of_islets_parquet_and_runs_after_both_rollups(tmp_path):
    """Since 2026-09-22 the two rollup scripts write only their summary files, which
    build_parquet_views.py merges; the dependency is the files, not a declared `after`."""
    nodes = bl.build_graph(REPO, tmp_path / "layers")
    islets = str(REPO / "data/development/typing/app/parquet/islets.parquet")
    assert bl.producers(nodes)[islets] == ["parquet_views"]
    order = bl.topological_order(nodes)
    assert order.index("crossed_fraction") < order.index("parquet_views")
    assert order.index("singlecell_trajectory") < order.index("parquet_views")
    for name in ("crossed_fraction", "singlecell_trajectory"):
        assert islets not in nodes[name].inputs and islets not in nodes[name].outputs


def test_two_writers_of_one_path_with_no_order_between_them_is_refused():
    """The guard guards: drop the `after` edge and the graph refuses to build."""
    a = bl.Node(name="a", command=("python", "a.py"), inputs=(), outputs=("/tmp/shared",))
    b = bl.Node(name="b", command=("python", "b.py"), inputs=(), outputs=("/tmp/shared",))
    with pytest.raises(ValueError, match="both write /tmp/shared"):
        bl._assert_shared_outputs_are_ordered({"a": a, "b": b})


def test_the_same_two_writers_are_accepted_once_an_order_is_declared():
    a = bl.Node(name="a", command=("python", "a.py"), inputs=(), outputs=("/tmp/shared",))
    b = bl.Node(name="b", command=("python", "b.py"), inputs=(), outputs=("/tmp/shared",),
                after=("a",))
    bl._assert_shared_outputs_are_ordered({"a": a, "b": b})
    assert bl.topological_order({"a": a, "b": b}) == ["a", "b"]


def test_the_two_rollups_are_independent_and_both_feed_parquet_views():
    nodes = bl.build_graph(REPO, Path("/tmp/layers-for-this-test"))
    dependencies = bl.edges(nodes)
    assert "singlecell_trajectory" not in dependencies["crossed_fraction"]
    assert "crossed_fraction" not in dependencies["singlecell_trajectory"]
    assert {"crossed_fraction", "singlecell_trajectory"} <= set(dependencies["parquet_views"])


def test_nonendocrine_outline_is_not_a_node(tmp_path):
    """It was one until 2026-09-21, and the cycle-break it needed went with it.

    The node streamed `data/redsea_scratch/geojson` -- 61 GB, up to 5.5 GB per donor, parsed one
    feature at a time -- to union each islet's endocrine cell outlines. The numbers that union
    fed are already measured per cell: `build_islet_geometry.py` takes the non-endocrine AREA
    from `cells_um2 = subset["cell_area"].sum()` and says so in its own docstring, and the
    STRUCTURE areas are the same evidence as `typing/crossed`'s per-cell `crossed_frac_*`. Only
    the drawn outline needs the polygon, so the script runs on demand instead.

    The cycle break this test used to assert -- a dead read of islet_geometry.parquet declared as
    a soft_input rather than an edge -- existed only to keep that node acyclic. With the node
    gone there is no cycle to break, and asserting on the vanished declaration would pin the
    scaffolding rather than the behaviour.
    """
    nodes = bl.build_graph(REPO, tmp_path / "layers")
    assert "nonendocrine_outline" not in nodes
    # islet_geometry still READS the unions; they are a source now, with a named producer.
    union = str(REPO / "data/development/typing/geometry")
    assert any(union in i for i in nodes["islet_geometry"].inputs), (
        "islet_geometry stopped reading the endocrine unions -- if that is intended, this test "
        "and the SOURCE_PRODUCERS entry should go too")
    assert any("typing/geometry" in k for k in bl.SOURCE_PRODUCERS), (
        "the unions are an input no node produces, so the check table must name what does")


# ---------------------------------------------------------------------------------------------
# The staleness gate
# ---------------------------------------------------------------------------------------------

def _node(tmp_path, inputs, outputs):
    return bl.Node(name="n", command=("python", "x.py"),
                   inputs=tuple(str(tmp_path / p) for p in inputs),
                   outputs=tuple(str(tmp_path / p) for p in outputs))


def test_an_output_newer_than_every_input_is_fresh(tmp_path):
    _touch(tmp_path / "in.csv", 1000)
    _touch(tmp_path / "out.parquet", 2000)
    status = bl.check_node(_node(tmp_path, ["in.csv"], ["out.parquet"]))
    assert status.state == bl.FRESH
    assert not status.is_red


def test_an_output_older_than_an_input_is_stale_and_names_the_input(tmp_path):
    _touch(tmp_path / "in.csv", 3000)
    _touch(tmp_path / "out.parquet", 2000)
    status = bl.check_node(_node(tmp_path, ["in.csv"], ["out.parquet"]))
    assert status.state == bl.STALE
    assert status.is_red
    assert status.newest_input == str(tmp_path / "in.csv")


def test_the_oldest_output_decides_not_the_newest(tmp_path):
    """One rewritten output must not hide a sibling that predates the inputs."""
    _touch(tmp_path / "in.csv", 3000)
    _touch(tmp_path / "old.parquet", 2000)
    _touch(tmp_path / "new.parquet", 4000)
    status = bl.check_node(_node(tmp_path, ["in.csv"], ["old.parquet", "new.parquet"]))
    assert status.state == bl.STALE
    assert status.oldest_output == str(tmp_path / "old.parquet")


def test_a_missing_output_is_named_rather_than_treated_as_fresh(tmp_path):
    _touch(tmp_path / "in.csv", 1000)
    status = bl.check_node(_node(tmp_path, ["in.csv"], ["never_built.parquet"]))
    assert status.state == bl.MISSING_OUTPUT
    assert status.is_red
    assert status.missing == (str(tmp_path / "never_built.parquet"),)


def test_a_missing_input_is_named_rather_than_treated_as_fresh(tmp_path):
    _touch(tmp_path / "out.parquet", 2000)
    status = bl.check_node(_node(tmp_path, ["gone.csv"], ["out.parquet"]))
    assert status.state == bl.MISSING_INPUT
    assert status.is_red
    assert status.missing == (str(tmp_path / "gone.csv"),)


def test_a_directory_input_takes_the_newest_file_in_the_tree(tmp_path):
    _touch(tmp_path / "tree/a/old.parquet", 1000)
    _touch(tmp_path / "tree/b/c/new.parquet", 5000)
    _touch(tmp_path / "out.parquet", 3000)
    status = bl.check_node(_node(tmp_path, ["tree"], ["out.parquet"]))
    assert status.state == bl.STALE
    assert status.newest_input == str(tmp_path / "tree/b/c/new.parquet")


def test_a_glob_input_only_matches_its_own_family(tmp_path):
    """typing/geometry holds two families from two builders; a glob keeps them apart."""
    _touch(tmp_path / "geometry/endocrine_union_115.parquet", 1000)
    _touch(tmp_path / "geometry/nonendocrine_115.geojson", 9000)
    _touch(tmp_path / "out.parquet", 3000)
    status = bl.check_node(_node(tmp_path, ["geometry/endocrine_union_*.parquet"],
                                 ["out.parquet"]))
    assert status.state == bl.FRESH


def test_a_glob_in_the_directory_part_is_resolved(tmp_path):
    """The pln layer's per-donor partitions are declared as donor_id=*/<file>."""
    _touch(tmp_path / "pln/donor_id=115/pln_node_membership.parquet", 7000)
    _touch(tmp_path / "in.csv", 1000)
    status = bl.check_node(_node(tmp_path, ["in.csv"],
                                 ["pln/donor_id=*/pln_node_membership.parquet"]))
    assert status.state == bl.FRESH
    assert status.oldest_output.endswith("donor_id=115/pln_node_membership.parquet")


def test_a_half_written_temporary_is_not_counted(tmp_path):
    """Every layer builder writes .parquet.tmp and renames; the temp must not read as an input."""
    _touch(tmp_path / "tree/data.parquet", 1000)
    _touch(tmp_path / "tree/data.parquet.tmp", 9000)
    _touch(tmp_path / "out.parquet", 3000)
    assert bl.check_node(_node(tmp_path, ["tree"], ["out.parquet"])).state == bl.FRESH


def test_a_pycache_directory_is_not_counted(tmp_path):
    _touch(tmp_path / "tree/data.parquet", 1000)
    _touch(tmp_path / "tree/__pycache__/x.pyc", 9000)
    _touch(tmp_path / "out.parquet", 3000)
    assert bl.check_node(_node(tmp_path, ["tree"], ["out.parquet"])).state == bl.FRESH


def test_check_reports_every_node_in_topological_order(tmp_path):
    nodes = bl.build_graph(REPO, tmp_path / "layers")
    statuses = bl.check(nodes)
    assert [s.name for s in statuses] == bl.topological_order(nodes)


def test_check_honours_only(tmp_path):
    nodes = bl.build_graph(REPO, tmp_path / "layers")
    statuses = bl.check(nodes, only=["sections", "trajectory"])
    assert [s.name for s in statuses] == ["sections", "trajectory"]


# ---------------------------------------------------------------------------------------------
# The command line
# ---------------------------------------------------------------------------------------------

def _driver(*args, layers: Path | None = None):
    argv = [sys.executable, str(DRIVER), *args]
    if layers is not None:
        argv += ["--layers-dir", str(layers)]
    return subprocess.run(argv, cwd=str(REPO), capture_output=True, text=True)


def test_check_exits_non_zero_when_an_artifact_is_missing(tmp_path):
    """Pointed at an empty layers root, every layer is absent and the gate must say so."""
    result = _driver("--check", layers=tmp_path / "empty")
    assert result.returncode != 0
    assert bl.MISSING_OUTPUT in result.stdout
    assert not (tmp_path / "empty").exists(), "--check must create nothing"


def test_plan_runs_nothing_and_writes_nothing(tmp_path):
    result = _driver("--plan", layers=tmp_path / "empty")
    assert result.returncode == 0
    assert "build_parquet_views.py" in result.stdout
    assert not (tmp_path / "empty").exists()


def test_plan_prints_both_rollups_before_parquet_views(tmp_path):
    result = _driver("--plan", layers=tmp_path / "empty")
    text = result.stdout
    assert text.index("summarise_crossed_fraction.py") < text.index("build_parquet_views.py")
    assert text.index("summarise_singlecell_trajectory.py") < text.index("build_parquet_views.py")


def test_an_unknown_node_name_is_refused_rather_than_silently_skipped(tmp_path):
    result = _driver("--check", "--only", "not_a_node", layers=tmp_path / "empty")
    assert result.returncode == 2
    assert "unknown node" in result.stderr


@pytest.mark.data
def test_run_refuses_a_node_whose_own_input_is_stale_and_builds_nothing(tmp_path):
    """The refusal is the point: a fresh artifact built on a stale input is the failure mode.

    Pointed at an empty layers root, islet_summary's upstream geometry nodes are not fresh, so the
    driver must refuse before it starts the builder and must leave the root uncreated.
    """
    target = tmp_path / "probe_layers"
    result = _driver("--run", "--only", "islet_summary", layers=target)
    assert result.returncode == 1
    assert "REFUSED" in result.stdout
    assert "islet_geometry" in result.stdout
    assert not target.exists(), "a refused run must write nothing"


def test_a_mode_is_required(tmp_path):
    result = _driver(layers=tmp_path / "empty")
    assert result.returncode == 2


# ---------------------------------------------------------------------------------------------
# run_pipeline.sh -- the wrapper that could not run
# ---------------------------------------------------------------------------------------------

def test_the_wrapper_does_not_default_to_a_repo_root_config():
    """Root CLAUDE.md: a repo-root config.ini is a bug. core/phenocycler owns the only one."""
    code = "\n".join(line for line in WRAPPER.read_text().splitlines()
                     if not line.lstrip().startswith("#"))
    assert '${REPO}/config.ini' not in code
    assert '${PACKAGE_DIR}/config.ini' in code
    assert not (REPO / "config.ini").exists(), "a repo-root config.ini reappeared"
    assert (REPO / "core/phenocycler/config.ini").is_file()


def test_the_wrapper_does_not_call_a_bare_python():
    code = "\n".join(line for line in WRAPPER.read_text().splitlines()
                     if not line.lstrip().startswith("#"))
    assert "exec python -m phenocycler" not in code
    assert '${REPO}/.venv/bin/python' in code


def test_the_wrapper_runs_from_core_phenocycler():
    """The workflow runs from core/phenocycler, which owns the config and the data roots."""
    assert 'cd "${PACKAGE_DIR}"' in WRAPPER.read_text()


def test_the_wrapper_is_syntactically_valid():
    assert subprocess.run(["bash", "-n", str(WRAPPER)]).returncode == 0


def test_the_wrapper_refuses_a_configuration_that_is_not_there(tmp_path):
    result = subprocess.run(["bash", str(WRAPPER), "status"], cwd=str(REPO),
                            env={**os.environ, "PHENOCYCLER_CONFIG": str(tmp_path / "nope.ini")},
                            capture_output=True, text=True)
    assert result.returncode == 2
    assert "no configuration at" in result.stderr


def test_there_is_exactly_one_runnable_wrapper_for_the_source_workflow():
    """Root CLAUDE.md: replacing a workflow removes the obsolete entrypoint, never leaves two."""
    wrappers = sorted(p for p in (REPO / "scripts").rglob("run_pipeline.sh"))
    assert wrappers == [WRAPPER]


# ---------------------------------------------------------------------------------------------
# Against the artifacts that are actually on disk
# ---------------------------------------------------------------------------------------------

@pytest.mark.data
def test_the_real_check_names_which_of_the_nine_layers_are_not_fresh():
    """No expected count is asserted: the gate's job is to name, and the number moves as lanes
    land.  What is pinned is that it reaches a verdict for every layer and says so out loud."""
    nodes = bl.build_graph()
    statuses = {s.name: s for s in bl.check(nodes)}
    for layer in bl.LAYER_NODES:
        assert layer in statuses
        assert statuses[layer].state in (bl.FRESH, bl.STALE, bl.MISSING_OUTPUT, bl.MISSING_INPUT)
    red = [n for n in bl.LAYER_NODES if statuses[n].is_red]
    if red:
        assert any(statuses[n].newest_input for n in red if statuses[n].state == bl.STALE)


@pytest.mark.data
def test_every_declared_input_that_is_not_produced_here_has_a_named_producer():
    """A source with no node must still be attributable, or `stale against X` is a dead end."""
    nodes = bl.build_graph()
    produced = set(bl.producers(nodes))
    unattributed = []
    for node in nodes.values():
        for path in node.inputs:
            if path in produced or path.endswith((".py", ".R", ".sh")):
                continue
            relative = str(Path(path)).replace(str(REPO) + "/", "")
            if not any(relative.startswith(k) or k in relative for k in bl.SOURCE_PRODUCERS):
                unattributed.append(relative)
    # Inventories and manifests are tracked files a person edits; they need no producer.
    unattributed = [p for p in unattributed
                    if not p.endswith((".csv", ".json", ".tsv", ".ini"))]
    assert not unattributed, f"sources with no node and no named producer: {sorted(set(unattributed))}"
