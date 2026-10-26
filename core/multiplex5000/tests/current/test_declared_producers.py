"""The declared producers -- one per tissue x modality -- and the guard that a class column
reaches the artifact graph from nowhere else.

    .venv/bin/python -m pytest core/multiplex5000/tests/current/test_declared_producers.py -q

``docs/PRODUCERS.md`` is the ledger this file holds to.  Four producers are declared there:

    pancreas x PhenoCycler   core/multiplex5000/scripts/typing/subpanel.py + build_cells.py
                             -> data/development/typing/cells            (cell_class)
    pancreas x Xenium        multiplex5000.tissues.pancreas_xenium_worker, run by
                             panc_xenium_release.py -> /data/xenium_typing/panc_2026-09-17/release
                                                                          (lineage, common_lineage)
    pLN x Xenium             multiplex5000.tissues.lymph_node_xenium_worker
                             -> data/development/xenium_typing/release_pln_immune_v3
                                                                          (lineage, common_lineage)
    pLN x PhenoCycler        core/multiplex5000/scripts/typing/build_pln_immune_pheno.py
                             -> data/development/typing/pln/pheno          (immune_class, ...)

The two Xenium roots are never spelled by a builder: ``scripts/senior/xenium_release.py``
resolves them from each release's own ``release.json``, searching its two roots
(``REPO_RELATIVE_ROOT`` inside the repo, ``EXTERNAL_ROOT_DEFAULT`` outside it).  A declared
input under either root is therefore a declared-producer root here; which release inside it is
current is the resolver's answer, and that answer is checked in ``test_xenium_release.py`` and,
against the real manifests, in the one ``data``-marked test at the bottom of this file.

WHAT "READS A CLASS COLUMN FROM A ROOT" MEANS HERE.  The graph declares inputs, not columns, so
the reading is measured two ways at once and both must hold: the node declares the root as an
input (``build_layers.py``), and the node's builder names one of the class columns that root
carries as a string constant in its CODE -- docstrings excluded, comments never parsed.  The
signature is coarse on purpose: a builder that opens a root and spells its class column is
reading it, whatever the surrounding control flow says about when.

THE LEAK THAT WAS MEASURED (2026-09-21) AND CLOSED (2026-09-23).  ``config.ini:153`` bound
``data/development/xenium_typing/release`` -- the scANVI two-donor release that the declared
pancreas producer supersedes (``panc_xenium_release.py:86-96``) -- as the integration lane's
Xenium identity release; ``phenospace.py:532-533`` wrote its ``lineage_native`` and
``lineage_common`` onto ``cells_xen``, and two nodes read them as a fallback when no typed
release resolved: ``availability`` (``xenium_ingest_availability``) and ``distance_profiles``
(``_pln_arm``).  Step 5 blanked both bindings at ``config.ini:152-153``, removed both readers,
and declared the release roots as inputs of the two nodes.  The guard below is now an ordinary
assertion; the binding test at the end pins the two keys BLANK so the leak cannot come back
unannounced.  Nothing here asserts a biological number, and nothing here reruns a builder.
"""
from __future__ import annotations

import ast
import importlib.util
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[4]
DRIVER = REPO / "scripts/senior/build_layers.py"
PRODUCERS_MD = REPO / "docs/PRODUCERS.md"
CONFIG_INI = REPO / "core/multiplex5000/config.ini"

sys.path.insert(0, str(REPO))
from scripts.senior import xenium_release as xr  # noqa: E402


def _load_driver():
    """``build_layers.py`` by path, the way ``test_build_layers.py`` loads it."""
    spec = importlib.util.spec_from_file_location("build_layers", DRIVER)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    # Registered before exec so the frozen dataclasses in it can resolve their own module.
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


bl = _load_driver()


# ---------------------------------------------------------------------------------------------
# The ledger, as docs/PRODUCERS.md declares it
# ---------------------------------------------------------------------------------------------

#: (tissue, modality) -> the scripts that ARE the producer.  Paths, not names: two files in
#: this tree are called aggregate_islets.py and three pins.py, so a basename is not an identity.
DECLARED_PRODUCERS: dict[tuple[str, str], tuple[str, ...]] = {
    ("pancreas", "PhenoCycler"): (
        "core/multiplex5000/scripts/typing/subpanel.py",
        "core/multiplex5000/scripts/typing/build_cells.py",
        # The refinement pass, joined back by `build_cells.py --refine` (ruled in 2026-09-23).
        "core/multiplex5000/scripts/typing/refine.py",
    ),
    ("pancreas", "Xenium"): (
        "core/multiplex5000/multiplex5000/tissues/pancreas_xenium_worker.py",
        "analysis/xenium/scripts/panc_xenium_release.py",
    ),
    ("pLN", "Xenium"): (
        "core/multiplex5000/multiplex5000/tissues/lymph_node_xenium_worker.py",
    ),
    ("pLN", "PhenoCycler"): (
        "core/multiplex5000/scripts/typing/build_pln_immune_pheno.py",
    ),
}

#: The output roots the four producers write, spelled as PRODUCERS.md spells them, so the doc
#: can be checked for naming each one.  The two Xenium spellings are the CURRENT releases as
#: measured 2026-09-22; the graph itself is checked against the resolver's roots, not these.
DECLARED_OUTPUT_ROOTS_AS_WRITTEN = (
    "data/development/typing/cells",
    "/data/xenium_typing/panc_2026-09-17/release",
    "data/development/xenium_typing/release_pln_immune_v3",
    "data/development/typing/pln/pheno",
)

#: The release ids the two current Xenium releases declare, read from their release.json on
#: 2026-09-22 (panc_xenium_release.py writes the first; pln_immune_typing_v3.py:103 the second).
DECLARED_RELEASE_IDS = {
    xr.PANC: "xenium-panc-cluster-name-20260917",
    xr.PLN: "xenium-pln-immune-tnk-split-dev-20260917",
}

#: Every root the graph declares that carries a class column, and which columns.  Measured
#: 2026-09-21 and re-read 2026-09-22:
#:   typing/cells      cell_class (subpanel.py:1214; build_cells.py:88), and refined_class
#:                     joined beside it by `build_cells.py --refine` (build_cells.py:148),
#:                     the refinement's verdict since the 2026-09-23 ruling
#:   typing/immune     the same lane's second pass over Immune (immune_subtype.py:98, :577)
#:   typing/pln/pheno  build_pln_immune_pheno.py:1732-1746
#:   the two resolver roots   lineage / common_lineage (panc_xenium_release.py:552;
#:                     pln_immune_typing_v3.py:1128-1129); the superseded releases under the
#:                     repo root spell lineage_native (xenium_typing/release.py:149)
#:   cells_xen         lineage_native / lineage_common, attached via phenospace.py:532-533
#:                     while config.ini:153 was bound (until 2026-09-23) -- NOT a producer, and
#:                     kept here so a builder that starts naming them again is caught
#:   cells_pheno       the same two columns from config.ini:152 via export_pheno.py:381-386
#:                     -- NOT a producer, and declared by no node
CLASS_COLUMNS_BY_ROOT: dict[Path, frozenset[str]] = {
    REPO / "data/development/typing/cells": frozenset({"cell_class", "refined_class"}),
    REPO / "data/development/typing/immune": frozenset({"immune_class", "immune_subtype"}),
    REPO / "data/development/typing/pln/pheno": frozenset(
        {"is_immune", "immune_class", "immune_subtype", "t_nk_class"}),
    REPO / xr.REPO_RELATIVE_ROOT: frozenset({"lineage", "common_lineage", "lineage_native"}),
    xr.EXTERNAL_ROOT_DEFAULT: frozenset({"lineage", "common_lineage"}),
    REPO / "core/phenocycler/data/integration/cells_xen": frozenset(
        {"lineage_native", "lineage_common"}),
    REPO / "core/phenocycler/data/integration/cells_pheno": frozenset(
        {"lineage_native", "lineage_common"}),
}

#: The roots a class column may come from: the four producers' outputs, with typing/immune as
#: producer 1's own second pass and the two resolver roots standing for producers 2 and 3.
DECLARED_PRODUCER_ROOTS: frozenset[Path] = frozenset({
    REPO / "data/development/typing/cells",
    REPO / "data/development/typing/immune",
    REPO / "data/development/typing/pln/pheno",
    REPO / xr.REPO_RELATIVE_ROOT,
    xr.EXTERNAL_ROOT_DEFAULT,
})

CELLS_XEN = REPO / "core/phenocycler/data/integration/cells_xen"
CELLS_PHENO = REPO / "core/phenocycler/data/integration/cells_pheno"

# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------

def _code_strings(path: Path) -> set[str]:
    """Every string constant in the module's code.  Docstrings are left out; comments are not
    in the AST at all.  A column named in a docstring is prose about the column, not a read."""
    tree = ast.parse(path.read_text(), filename=str(path))
    docstrings: set[int] = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            body = node.body
            if (body and isinstance(body[0], ast.Expr)
                    and isinstance(body[0].value, ast.Constant)
                    and isinstance(body[0].value.value, str)):
                docstrings.add(id(body[0].value))
    return {n.value for n in ast.walk(tree)
            if isinstance(n, ast.Constant) and isinstance(n.value, str)
            and id(n) not in docstrings}


def _under(path: str, root: Path) -> bool:
    p = Path(path)
    return p == root or root in p.parents


def _builder(node) -> Path:
    """The script a node runs.  Every node runs the repo venv on a script (test_build_layers
    asserts the interpreter); the script is the second token."""
    assert len(node.command) >= 2, f"{node.name}: command has no script: {node.command}"
    return Path(node.command[1])


def _class_reads(nodes) -> set[tuple[str, Path]]:
    """(node, root) for every declared input that carries a class column the node's builder
    names in code.  Both conditions, as the module docstring says."""
    out: set[tuple[str, Path]] = set()
    for node in nodes.values():
        names = _code_strings(_builder(node))
        for declared in node.inputs:
            for root, columns in CLASS_COLUMNS_BY_ROOT.items():
                if _under(declared, root) and names & columns:
                    out.add((node.name, root))
    return out


def _leaks(nodes) -> set[tuple[str, Path]]:
    return {(name, root) for name, root in _class_reads(nodes)
            if root not in DECLARED_PRODUCER_ROOTS}


def _rel(root: Path) -> str:
    try:
        return str(root.relative_to(REPO))
    except ValueError:
        return str(root)


@pytest.fixture(scope="module")
def nodes(tmp_path_factory):
    return bl.build_graph(REPO, tmp_path_factory.mktemp("layers"))


# ---------------------------------------------------------------------------------------------
# The producers exist, and the ledger names them
# ---------------------------------------------------------------------------------------------

def test_the_four_declared_producers_exist_on_disk():
    missing = [p for scripts in DECLARED_PRODUCERS.values() for p in scripts
               if not (REPO / p).is_file()]
    assert not missing, f"declared producers that are not on disk: {missing}"
    assert len(DECLARED_PRODUCERS) == 4
    assert {t for t, _ in DECLARED_PRODUCERS} == {"pancreas", "pLN"}
    assert {m for _, m in DECLARED_PRODUCERS} == {"PhenoCycler", "Xenium"}


def test_producers_md_names_each_declared_producer_with_its_path():
    assert PRODUCERS_MD.is_file(), "docs/PRODUCERS.md is the ledger this file holds to"
    text = PRODUCERS_MD.read_text()
    unnamed = [p for scripts in DECLARED_PRODUCERS.values() for p in scripts if p not in text]
    assert not unnamed, f"docs/PRODUCERS.md does not name {unnamed} by path"
    roots_unnamed = [r for r in DECLARED_OUTPUT_ROOTS_AS_WRITTEN if r not in text]
    assert not roots_unnamed, f"docs/PRODUCERS.md does not name the output root(s) {roots_unnamed}"
    for key in ("phenocycler_identity_release", "xenium_identity_release", "config.ini:152",
                "config.ini:153"):
        assert key in text, f"docs/PRODUCERS.md must record the measured leak by name: {key}"


def test_every_declared_producer_root_is_spelled_by_its_producer_or_its_resolver():
    """The roots this file checks the graph against are the producers' own spellings, read out
    of their modules rather than restated: build_cells.OUT, build_pln_immune_pheno.OUT, and the
    resolver's two search roots.  If a producer moves its output, this is the test that says
    the table above has drifted from it."""
    typing = REPO / "core/multiplex5000/scripts/typing"
    text_cells = (typing / "build_cells.py").read_text()
    assert 'OUT = S.OUT / "cells"' in text_cells, "build_cells.py:88 no longer writes S.OUT/cells"
    text_pln = (typing / "build_pln_immune_pheno.py").read_text()
    assert 'OUT = REPO / "data/development/typing/pln/pheno"' in text_pln, (
        "build_pln_immune_pheno.py:168 no longer writes data/development/typing/pln/pheno")
    assert REPO / xr.REPO_RELATIVE_ROOT == REPO / "data/development/xenium_typing"
    assert xr.EXTERNAL_ROOT_DEFAULT == Path("/data/xenium_typing")
    assert set(xr.default_roots(REPO)) == {REPO / xr.REPO_RELATIVE_ROOT, xr.EXTERNAL_ROOT_DEFAULT}


# ---------------------------------------------------------------------------------------------
# The guard
# ---------------------------------------------------------------------------------------------

def test_class_columns_reach_the_graph_only_from_declared_producer_roots(nodes):
    """Every node reading a class column reads it only from the four declared producers'
    roots.  cells_xen, whose lineage columns were attached from config.ini:153 until
    2026-09-23, is read for a class column by no node: availability and distance_profiles
    still declare it (IF pairing; a pancreas cell's coordinates) and neither builder names
    lineage_native or lineage_common in code.  A third node declaring cells_xen and naming
    those columns, a node declaring cells_pheno, or a release spelled outside the resolver's
    roots all land here."""
    leaks = _leaks(nodes)
    assert not leaks, "class columns read from a root that is not a declared producer: " + \
        ", ".join(f"{name} <- {_rel(root)}" for name, root in sorted(leaks, key=str))


def test_the_two_nodes_that_read_the_leak_still_declare_the_ingest_and_the_releases(nodes):
    """What replaced the fallback: distance_profiles declares the two declared release roots
    as inputs (so a new release makes it stale); both nodes keep cells_xen, declared for the
    columns only the ingest carries, and neither builder names its lineage columns.

    availability declares NO release root since 2026-09-23 (Step 5c): it dropped xenium_panc
    and xenium_pln, which no card read, and with them its only read of a release.  A release
    root declared on a node that does not read it would make it stale for nothing."""
    for name in ("availability", "distance_profiles"):
        node = nodes[name]
        assert any(_under(i, CELLS_XEN) for i in node.inputs), f"{name} no longer declares cells_xen"
        named = _code_strings(_builder(node)) & CLASS_COLUMNS_BY_ROOT[CELLS_XEN]
        assert not named, f"{name}'s builder names a cells_xen lineage column in code: {named}"
    profiles = nodes["distance_profiles"]
    assert any(_under(i, xr.EXTERNAL_ROOT_DEFAULT) for i in profiles.inputs), (
        "distance_profiles does not declare the external release root as an input")
    assert any(_under(i, REPO / xr.REPO_RELATIVE_ROOT) for i in profiles.inputs), (
        "distance_profiles does not declare the repo release root as an input")
    availability = nodes["availability"]
    imported = set()
    for node in ast.walk(ast.parse(_builder(availability).read_text())):
        if isinstance(node, ast.Import):
            imported |= {alias.name for alias in node.names}
        elif isinstance(node, ast.ImportFrom):
            imported |= {f"{node.module}.{alias.name}" for alias in node.names}
    assert not {name for name in imported if name.endswith("xenium_release")}, (
        "build_availability.py imports the release resolver again")
    assert not any(_under(i, xr.EXTERNAL_ROOT_DEFAULT) or _under(i, REPO / xr.REPO_RELATIVE_ROOT)
                   for i in availability.inputs), (
        "availability declares a release root, and its builder reads no release")


def test_cells_pheno_is_a_declared_input_of_no_node(nodes):
    """config.ini:152's twin.  export_pheno.py:381-386 attaches the withdrawn broad-typing
    release onto cells_pheno, and the census measured that no build_layers node declares it.
    The report lane reads it (identity_io.py:19); the graph must not."""
    declaring = sorted(n.name for n in nodes.values()
                       if any(_under(i, CELLS_PHENO) for i in n.inputs))
    assert not declaring, f"cells_pheno has become an input of the graph: {declaring}"


def test_every_declared_xenium_release_input_sits_under_a_resolver_root(nodes):
    """No node spells a Xenium release outside the two roots xenium_release.py searches.  A
    release declared elsewhere would be one the resolver can never select, read by a builder
    that pinned it -- the 2026-09-18 defect xenium_release.py exists to end."""
    roots = (REPO / xr.REPO_RELATIVE_ROOT, xr.EXTERNAL_ROOT_DEFAULT)
    stray = sorted((n.name, i) for n in nodes.values() for i in n.inputs
                   if "xenium_typing" in i and not any(_under(i, r) for r in roots))
    assert not stray, f"Xenium release inputs outside the resolver's roots: {stray}"


def test_no_node_declares_the_superseded_pln_release(nodes):
    """release_pln_immune_v3 declares that it supersedes _v2 (pln_immune_typing_v3.py:1484),
    and the resolver honours that, so no builder reads _v2 any more. A node still declaring
    _v2 has the staleness check watching a tree its builder never opens -- which is what
    `compartments` did until 2026-09-23. The maintainer's go to fix it came the same day."""
    stale = sorted((n.name, i) for n in nodes.values() for i in n.inputs
                   if i.rstrip("/").endswith("release_pln_immune_v2"))
    assert not stale, f"nodes declaring the superseded pLN release: {stale}"
    assert any(i.rstrip("/").endswith("release_pln_immune_v3")
               for i in nodes["compartments"].inputs), (
        "compartments resolves the pLN release (build_compartment_gradient.py:870) and must "
        "declare the v3 root it reads")


def test_both_identity_release_bindings_are_blank():
    """config.ini:152-153 are the two bindings that carried a second method into the ingest
    trees.  Both keys stay in the file (the ledger cites them by line) and both are BLANK,
    which config.ini:151 names as the fail-closed state.  A value here is the leak coming
    back: phenospace.load_release_or_none attaches whatever the key points at."""
    import configparser
    cfg = configparser.ConfigParser(interpolation=None)
    cfg.read(CONFIG_INI)
    values = {k: v for section in cfg.sections() for k, v in cfg.items(section)
              if k in ("xenium_identity_release", "phenocycler_identity_release")}
    assert set(values) == {"xenium_identity_release", "phenocycler_identity_release"}, (
        f"config.ini no longer carries both identity-release keys: {sorted(values)}")
    for key, value in values.items():
        assert value.strip() == "", (
            f"{key} = {value!r}: a declared identity release attaches a class column to the "
            "ingest that no declared producer wrote. Blank is the fail-closed state "
            "(config.ini:151; docs/PRODUCERS.md)")


# ---------------------------------------------------------------------------------------------
# Against the real releases (deselected by default; `make test-data`)
# ---------------------------------------------------------------------------------------------

@pytest.mark.data
def test_the_resolver_selects_the_declared_producers_releases_on_this_checkout():
    """The two Xenium roots the graph is checked against are search roots; this is the one
    place the selection inside them is checked against the ledger: the resolver's answer for
    each tissue is the release the declared producer wrote, by its own release_id."""
    expected_root = {
        xr.PANC: Path("/data/xenium_typing/panc_2026-09-17/release"),
        xr.PLN: REPO / "data/development/xenium_typing/release_pln_immune_v3",
    }
    for tissue, release_id in DECLARED_RELEASE_IDS.items():
        found = xr.resolve(tissue, repo=REPO)
        if not found.available:
            # The pattern test_xenium_release.py uses: a worktree carries no data/, and the
            # pLN release lives only under the main checkout's. A skip is reported, not hidden.
            pytest.skip(f"no {tissue} release on this checkout: {found.reason}")
        assert found.release_id == release_id, (
            f"{tissue}: the resolver selects {found.release_id} at {found.root}, but the "
            f"ledger declares {release_id}; docs/PRODUCERS.md has drifted from the tree")
        assert found.root == expected_root[tissue], (
            f"{tissue}: {release_id} now lives at {found.root}; docs/PRODUCERS.md spells "
            f"{expected_root[tissue]}")
    for root in (REPO / "data/development/typing/cells", REPO / "data/development/typing/pln/pheno"):
        if not root.is_dir():
            pytest.skip(f"{root} is not on this checkout")
        assert any(root.glob("donor_id=*")), f"{root} holds no donor partition"
