"""No leak: every class the ONE approved producer can emit has to survive to the far end.

The maintainer's ruling is that the islet chain (CD99 as the key differentiator, elongation in
the larger islets only, per-donor AND per-stratum reference levels) and the recovered refinement
pass both live in the approved implementation, and that nothing downstream keeps the old
behaviour.  "Nothing downstream" is the part a producer test cannot see: ``subpanel.py`` and
``refine.py`` can be perfectly correct while a consumer three files away maps their output
through a nine-key dictionary written when there were eight classes, and the cells fall out of
the numerator AND the denominator without a line of output.

So this file asserts about the CONSUMERS, and it DERIVES the vocabulary rather than restating it:

  * the class list comes from ``subpanel.STRUCTURES``, the string literals ``subpanel.py``
    actually assigns into ``cell_class``, and ``refine.OWN_MARKERS`` / ``UNTOUCHED`` /
    ``UNRESOLVED`` -- cross-checked against each other so neither module can drift alone, and
    against ``subpanel_log.json``'s hormone names when that file is on disk;
  * the consumers come from reading the tree, not from a list kept here.  Each test names every
    file it inspected in its failure message, so the coverage is legible and a new builder that
    nobody told this test about still gets checked.

This guard is EXPECTED TO BE RED until the rerun and the wiring that follows it.  A guard that
passes on the tree it was written against has not been shown to catch anything.  What each test
catches today is written in its own docstring, with the counts, so a reader can tell a leak that
is still open from a regression that arrived later.

Nothing here reads a result.  No test asserts a cell count, a proportion or a demotion rate: a
check that passed by reproducing a previous number would be measuring the expectation.  The two
tests that touch ``data/`` read a LOG of what the method decided (marker names, class names),
never an outcome, and skip cleanly when it is absent.

Run from the repository root:

    .venv/bin/python -m pytest core/phenocycler/tests/current/test_no_typing_leak.py -v
"""

from __future__ import annotations

import ast
import importlib.util
import json
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest
import yaml

REPO = pathlib.Path(__file__).resolve().parents[4]
TYPING = REPO / "core/phenocycler/scripts/typing"
SENIOR = REPO / "scripts/senior"
RUN20 = REPO / "apps/islet/config/run20_default.yml"

#: The modules that DECIDE a class.  Everything else that reads one is a consumer, and it is
#: the consumers this file is about.  Kept as filenames rather than as a list of "the eight
#: classes" on purpose: a name of a file is a fact about the tree, a list of classes would be
#: the very restatement this test exists to catch.
PRODUCERS = ("subpanel.py", "subpanel_io.py", "context.py", "refine.py")

#: A file naming this many of the producer's base classes as literals is not naming a subgroup
#: -- the two real subgroups are three long (``STRUCTURAL`` is Blood vessel / Nerve / Lymphatic,
#: the hormones are Beta / Alpha / Delta) -- it is enumerating the vocabulary, and an
#: enumeration of the vocabulary has to be complete.
ENUMERATION_THRESHOLD = 6

#: Aggregating calls that turn a class column into counts.  A class missing from one of these
#: is missing from both the numerator and the denominator it would have been part of.
AGGREGATORS = ("groupby", "value_counts", "pivot_table", "pivot", "crosstab", "size")
CLASS_COLUMNS = ("phenotype", "cell_class", "refined_class")


# --------------------------------------------------------------------------- loading

def _import(name: str):
    """Import a typing module by path, with the two search paths its own CLI sets up."""
    for extra in (TYPING, TYPING.parent / "celltype"):
        if str(extra) not in sys.path:
            sys.path.insert(0, str(extra))
    spec = importlib.util.spec_from_file_location(name, TYPING / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _tree(path: pathlib.Path) -> ast.Module:
    return ast.parse(path.read_text(), filename=str(path))


def _consumer_paths() -> list[pathlib.Path]:
    """Every builder under the two typing lanes that reads a cell class, producers excluded."""
    out = []
    for directory in (TYPING, SENIOR):
        for path in sorted(directory.glob("*.py")):
            if path.name in PRODUCERS:
                continue
            text = path.read_text()
            if any(col in text for col in CLASS_COLUMNS):
                out.append(path)
    return out


#: How a file says it reads the run cohort's per-cell table: it either goes through the lane's
#: own path module or it names the stage root.  Derived rather than listed, so a new builder is
#: in scope the moment it opens the table -- and the CODEX-legacy scripts that share the
#: ``scripts/senior/`` directory but read the other domain's h5ad files stay out of scope,
#: because this ruling is about the PhenoCycler run cohort and they are not part of it.
RUN_COHORT_MARKERS = ("import subpanel_io", "typing/cells", '"typing", "cells"',
                      "'typing', 'cells'")


def _run_cohort_consumers() -> list[pathlib.Path]:
    """The consumers that actually read the table the producer writes."""
    return [p for p in _consumer_paths()
            if any(m in p.read_text() for m in RUN_COHORT_MARKERS)]


def _rel(paths) -> list[str]:
    return [str(p.relative_to(REPO)) for p in paths]


def _string_constants(node: ast.AST) -> set[str]:
    return {n.value for n in ast.walk(node)
            if isinstance(n, ast.Constant) and isinstance(n.value, str)}


# ----------------------------------------------------------------- the vocabulary

def _chain_literals() -> set[str]:
    """Class names ``subpanel.py`` assigns into ``cell_class`` -- read out of the chain itself."""
    found = set()
    for node in ast.walk(_tree(TYPING / "subpanel.py")):
        if not isinstance(node, ast.Assign):
            continue
        for target in node.targets:
            if (isinstance(target, ast.Subscript) and isinstance(target.value, ast.Name)
                    and target.value.id == "cell_class"
                    and isinstance(node.value, ast.Constant)
                    and isinstance(node.value.value, str)):
                found.add(node.value.value)
    return found


@pytest.fixture(scope="module")
def vocabulary() -> dict:
    """Every class the producer can emit, derived from the two modules and cross-checked.

    ``refine.OWN_MARKERS`` is the only place the eight base classes are enumerated in one
    object, so it is the spine -- but it is only trusted because the structures and the chain
    literals read independently out of ``subpanel.py`` are a subset of it.  If either module
    grows a class the other does not know about, the cross-check fails here rather than
    surfacing as a class with no cells in it six steps later.
    """
    subpanel = _import("subpanel")
    refine = _import("refine")

    spine = set(refine.OWN_MARKERS) | set(refine.UNTOUCHED) | {refine.UNRESOLVED}
    structures = set(subpanel.STRUCTURES)
    literals = _chain_literals()
    hormones = spine - structures - literals - {refine.UNRESOLVED}
    return {"all": spine, "base": set(refine.OWN_MARKERS), "structures": structures,
            "chain_literals": literals, "hormones": hormones,
            "unresolved": refine.UNRESOLVED, "untouched": set(refine.UNTOUCHED),
            "refine": refine, "subpanel": subpanel}


def test_the_two_producer_modules_agree_on_the_vocabulary(vocabulary):
    """Foundation for every other test here: the class list can be derived, and it is one list.

    GREEN today, and it has to stay green -- if this one goes red the others are asserting
    against a vocabulary that no longer describes the producer.
    """
    spine = vocabulary["all"]
    drifted = (vocabulary["structures"] | vocabulary["chain_literals"]) - spine
    assert not drifted, (
        f"subpanel.py emits {sorted(drifted)}, which refine.py's OWN_MARKERS / UNTOUCHED / "
        f"UNRESOLVED do not know about. The refinement would then be asked to keep or demote a "
        f"class it has no own-marker test for. Producer vocabulary derived here: {sorted(spine)}")
    assert vocabulary["hormones"], (
        "no hormone class names survive the derivation -- subpanel.py's structures plus its "
        "cell_class literals already account for every class refine.py declares, which means "
        "Beta/Alpha/Delta have gone missing from one side or the other")

    # The hormone names are decided by the method and recorded in its log. Read as INPUT --
    # which classes exist -- never as a result.
    log = TYPING.parent.parent.parent.parent / "data/development/typing/subpanel/subpanel_log.json"
    if log.exists():
        declared = set(json.loads(log.read_text()).get("hormones", {}))
        assert declared == vocabulary["hormones"], (
            f"the hormone classes derived from code {sorted(vocabulary['hormones'])} are not the "
            f"ones the method's own log declares {sorted(declared)}")


# --------------------------------------------------- 1. the app's declared vocabulary

def _explicitly_removed(consumers) -> set[str]:
    """Classes a consumer drops BY NAME -- ``cell_class != "Excluded"`` and friends.

    A class removed by name is accounted for: somebody wrote down that it goes and can be asked
    why. That is the only reason a producer class is allowed to be absent from the app's
    display vocabulary, and the set is read out of the consumers rather than listed here.
    """
    removed = set()
    for path in consumers:
        for node in ast.walk(_tree(path)):
            if not isinstance(node, ast.Compare):
                continue
            if "cell_class" not in _string_constants(node.left):
                continue
            for op, comp in zip(node.ops, node.comparators):
                # ONLY ``!=``. ``cell_class == "Immune"`` is a selection -- immune_subtype.py
                # reads the Immune cells and writes nothing back -- and a selection says
                # nothing about whether the class survives to the app.
                if isinstance(op, ast.NotEq) and isinstance(comp, ast.Constant) \
                        and isinstance(comp.value, str):
                    removed.add(comp.value)
    return removed


def test_the_app_declares_every_class_the_producer_can_emit(vocabulary):
    """A ninth class must not be able to exist in the data and be invisible in the app.

    WAS RED, WENT GREEN DURING THIS SESSION as the app lane landed. When this test was written,
    ``build_cells.DISPLAY`` had nine keys and ``run20_default.yml`` ``phenotypes.display_order``
    had eight entries, and ``Unresolved (no lineage evidence)`` was in neither -- so the class
    the refinement can emit had nowhere to go. ``build_cells`` raises ``KeyError`` on an unmapped
    class, so the failure mode was never a silent hole in the app: it was that the per-cell table
    would refuse to build at all the moment the refinement was wired in. Both halves have moved.

    It also checks the restated constant. ``build_cells`` deliberately does NOT import
    ``refine`` (that would pull in the whole celltype scoring stack) and restates the label
    instead, so two spellings of the same class now exist in the tree and only an assertion
    keeps them one class.
    """
    build_cells = _import("build_cells")
    display = dict(build_cells.DISPLAY)
    config = yaml.safe_load(RUN20.read_text())
    order = list(config["phenotypes"]["display_order"])
    colors = dict(config["phenotypes"]["colors"])

    removed = _explicitly_removed(_run_cohort_consumers())
    must_display = sorted(vocabulary["all"] - removed)

    unmapped = [c for c in must_display if c not in display]
    assert not unmapped, (
        f"{unmapped} can be emitted by the producer and have no entry in build_cells.DISPLAY, "
        f"which is the one crosswalk from the producer's words to the app's. build_cells raises "
        f"KeyError on an unmapped class, so the per-cell table will not build at all. "
        f"DISPLAY covers {sorted(display)}; classes dropped by name elsewhere and therefore not "
        f"required to display: {sorted(removed)}")

    missing = [display[c] for c in must_display if display[c] not in order]
    assert not missing, (
        f"{missing} would be written into the cells table and never appear in the app: "
        f"run20_default.yml phenotypes.display_order is {order}. A class that exists in the "
        f"data and not in the vocabulary is counted by every builder and shown by nothing")

    uncoloured = [display[c] for c in must_display if display[c] not in colors]
    assert not uncoloured, (
        f"{uncoloured} are missing from phenotypes.colors; ieCore::domain_config refuses a "
        f"domain whose colors do not cover display_order, so the app will not start")

    restated = getattr(build_cells, "UNRESOLVED", None)
    if restated is not None:
        assert restated == vocabulary["unresolved"], (
            f"build_cells restates the unresolved label as {restated!r} while refine.py produces "
            f"{vocabulary['unresolved']!r}. Two spellings are two classes: the join would map "
            "one of them to nothing and the cells would leave the table without a word")


# ------------------------------------------- 2. no consumer filters to the old eight

def test_no_consumer_filters_cell_class_to_the_old_eight(vocabulary):
    """A builder that enumerates the vocabulary must enumerate all of it.

    Caught three of the 21 run-cohort consumers when it was written; ``build_cells.py``
    (``DISPLAY``, then 8 base classes + ``Excluded``) was fixed during this session and has
    dropped out. RED ON TWO, for two different reasons:

      * ``scripts/senior/build_grid_composition.py`` (``PANC_CLASS_TO_COMMON``, 7 keys -- Immune
        is routed by its definition instead) predates the refinement and omits only the new
        class. Its loop is ``if value in PANC_CLASS_TO_COMMON``, so an unmapped class is skipped
        with no branch and no count.
      * ``scripts/senior/modality_concordance.py`` (``PANC_PHENO_CLASS_TO_COMMON``) omits the new
        class AND ``Lymphatic`` AND ``Excluded``, and those two halves predate this work
        entirely: ``common.value_counts(dropna=True)`` drops the unmapped rows, so Lymphatic and
        Excluded cells leave the lineage counts while ``n_cells`` still counts them and the
        per-lineage numbers no longer sum to the total they are reported against.

    The threshold is six of eight because the two real subgroups are three long -- the
    structural classes and the hormones -- so a file naming six is naming the vocabulary and an
    ordinary ``STRUCTURAL``/hormone tuple is not mistaken for one.
    """
    consumers = _run_cohort_consumers()
    base = vocabulary["base"]
    offenders = {}
    for path in consumers:
        literals = _string_constants(_tree(path)) & vocabulary["all"]
        if len(literals & base) < ENUMERATION_THRESHOLD:
            continue
        missing = sorted(vocabulary["all"] - literals)
        if missing:
            offenders[str(path.relative_to(REPO))] = missing

    assert not offenders, (
        "these builders enumerate the producer's classes as literals and the enumeration is "
        "incomplete, so the missing classes are dropped or mapped to nothing without a word:\n"
        + "\n".join(f"  {f}: missing {m}" for f, m in sorted(offenders.items()))
        + f"\n\nfiles inspected ({len(consumers)}): {_rel(consumers)}")


# -------------------------------------------- 3. Unresolved keeps its row and its denominator

def _aggregation_sites(path: pathlib.Path) -> list[tuple[int, str, bool]]:
    """(line, calls, drops_nulls) for every place this file counts by a class column.

    Sites are collapsed per line because the fix lands on one call of a chain and has to count
    for the chain: ``frame.groupby([...], dropna=False).size()`` keeps the null group, but the
    ``dropna`` sits on the ``groupby`` and the ``.size()`` would otherwise be reported as still
    dropping it -- and a guard that cannot be satisfied by the correct fix is not a guard.
    """
    per_line: dict[int, dict] = {}
    for node in ast.walk(_tree(path)):
        if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr in AGGREGATORS):
            continue
        named = set()
        for arg in list(node.args) + [kw.value for kw in node.keywords]:
            named |= _string_constants(arg)
        named |= _string_constants(node.func.value)
        if not named & set(CLASS_COLUMNS):
            continue
        keeps_nulls = any(kw.arg == "dropna" and isinstance(kw.value, ast.Constant)
                          and kw.value.value is False for kw in node.keywords)
        entry = per_line.setdefault(node.lineno, {"calls": set(), "keeps": False})
        entry["calls"].add(node.func.attr)
        entry["keeps"] |= keeps_nulls
    return [(line, ".".join(sorted(e["calls"])), not e["keeps"])
            for line, e in sorted(per_line.items())]


def test_unresolved_is_never_a_zero_and_never_leaves_a_denominator(vocabulary):
    """A cell the refinement could not resolve must be carried, not erased.

    RED TODAY on all three composition builders that read the run cohort's per-cell table --
    ``aggregate_islets.py``, ``build_workbook.py`` and ``neighborhood.py``. None of them names
    the unresolved class, and all of them aggregate with pandas' default. The concrete path:
    ``build_cells`` maps ``cell_class`` through ``DISPLAY`` into ``phenotype``; an unmapped
    class becomes NaN; ``aggregate_islets.build_donor`` groups by ``phenotype`` and then takes
    ``region_total`` as the sum of those same groups -- so an unresolved cell leaves the
    numerator and the denominator in one step, and every remaining proportion silently rises.

    The requirement is not that a builder must report the class; it is that the builder must
    either carry it as its own row (``dropna=False``) or name it, so that its absence is a
    decision on the page rather than a default in pandas.

    Out of scope, and deliberately: ``scripts/senior/fixed_islet_aggregation.py`` and
    ``scripts/senior/build_phenotyped_h5ad.py`` also aggregate a ``phenotype`` column, but they
    read the legacy CODEX h5ad lane rather than this cohort's per-cell table. The two datasets
    are never merged, so a class this producer can emit is not their problem.
    """
    unresolved = vocabulary["unresolved"]
    token = unresolved.split(" (")[0]          # "Unresolved" -- the name, not the parenthetical
    consumers = _run_cohort_consumers()
    checked, offenders = [], {}
    for path in consumers:
        sites = _aggregation_sites(path)
        if not sites:
            continue
        checked.append(path)
        text = path.read_text()
        names_it = token in text or "refine" in text
        silent = [f"line {ln}: {calls}()" for ln, calls, drops in sites if drops]
        if silent and not names_it:
            offenders[str(path.relative_to(REPO))] = silent

    assert checked, (
        "no composition builder was found at all, which means the discovery in this test is "
        f"broken rather than that the tree is clean; files read: {_rel(consumers)}")
    assert not offenders, (
        f"{unresolved!r} is dropped from these counts by pandas' default and named in none of "
        "them, so it leaves the denominator without a line of output:\n"
        + "\n".join(f"  {f}:\n" + "\n".join(f"      {s}" for s in sites)
                    for f, sites in sorted(offenders.items()))
        + f"\n\ncomposition builders inspected ({len(checked)}): {_rel(checked)}")


# ------------------------------------------------------- 4. the refinement is single-pass

ITERATION_WORDS = ("pass", "passes", "iter", "iterations", "step", "steps",
                   "round", "rounds", "repeat", "converge")


def test_the_refinement_exposes_no_way_to_iterate_it(vocabulary):
    """No argument, anywhere in refine.py, lets a caller ask for a second pass.

    GREEN today. It is asserted because the ratchet it prevents is invisible in the output: a
    second pass re-derives its bound from a population the first pass already cut, so it always
    finds another five per cent to cut and the result looks tidier every time.
    """
    source = _tree(TYPING / "refine.py")
    offenders = []
    for node in ast.walk(source):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            args = node.args
            for arg in (list(args.posonlyargs) + list(args.args) + list(args.kwonlyargs)):
                if any(w == arg.arg or arg.arg.endswith(f"_{w}") or arg.arg.startswith(f"{w}_")
                       for w in ITERATION_WORDS):
                    offenders.append(f"{node.name}(..., {arg.arg}, ...) line {node.lineno}")
        if (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "add_argument"):
            for arg in node.args:
                if isinstance(arg, ast.Constant) and isinstance(arg.value, str):
                    flag = arg.value.lstrip("-").replace("-", "_")
                    if any(w == flag for w in ITERATION_WORDS):
                        offenders.append(f"argparse {arg.value} line {node.lineno}")

    assert not offenders, (
        "refine.py takes an iteration count, so a caller can run the pass to convergence "
        f"without editing it: {offenders}. The pass is single by measurement, not by "
        "convention -- the second pass was measured and rejected (REFINEMENT.md, "
        "'Step 2, tested and rejected'); a caller that wants one has to answer that first")


def test_no_caller_loops_the_refinement():
    """Nothing in the tree calls the pass from inside a loop.

    GREEN today -- and green for the wrong reason, which the next test is about: nothing calls
    it at all. It will start doing real work the moment the pass is wired into a builder.
    """
    roots = (REPO / "core/phenocycler/scripts", REPO / "scripts",
             REPO / "core/phenocycler/tests", REPO / "apps")
    entrypoints = {"apply_once", "refine_donor", "main"}
    offenders = []
    for root in roots:
        for path in sorted(root.rglob("*.py")):
            if path.name == "refine.py":
                continue
            try:
                tree = _tree(path)
            except SyntaxError:
                continue
            loops = []
            for node in ast.walk(tree):
                if isinstance(node, (ast.For, ast.AsyncFor, ast.While)):
                    loops.append(node)
            for loop in loops:
                for node in ast.walk(loop):
                    if not isinstance(node, ast.Call):
                        continue
                    name = (node.func.attr if isinstance(node.func, ast.Attribute)
                            else getattr(node.func, "id", None))
                    if name not in entrypoints:
                        continue
                    receiver = (node.func.value.id if isinstance(node.func, ast.Attribute)
                                and isinstance(node.func.value, ast.Name) else "")
                    if receiver == "refine" or (name != "main" and "refine" in path.read_text()):
                        offenders.append(f"{path.relative_to(REPO)}:{node.lineno} {name}()")

    assert not offenders, (
        f"the refinement is called from inside a loop: {offenders}. Looping it is the ratchet "
        "the single-pass guard exists to stop -- it does not raise on the second iteration of a "
        "loop that rebuilds the labels array each time")


def test_the_single_pass_guards_actually_raise(vocabulary):
    """The guards are exercised, not just read: both entry points refuse an already-refined input.

    GREEN today. Synthetic inputs only -- three rows and a fake neighbour index; nothing here
    depends on a cohort number.
    """
    refine = vocabulary["refine"]
    n = 3
    Z = np.zeros((n, 2), dtype=float)
    labels = np.array(["Exocrine", refine.UNRESOLVED, "Immune"], dtype=object)
    strata = np.array(["", "", ""], dtype=object)
    leaf = np.zeros(n, dtype=bool)
    nbr = np.zeros((n, 1), dtype=int)

    with pytest.raises(refine.RefinementAlreadyApplied):
        refine.apply_once(Z, ["INS", "GCG"], labels, strata, leaf, nbr=nbr)

    with pytest.raises(refine.RefinementAlreadyApplied):
        refine._refuse_second_pass(pd.DataFrame({"cell_class": [refine.UNRESOLVED]}))

    with pytest.raises(refine.RefinementAlreadyApplied):
        refine._refuse_second_pass(
            pd.DataFrame({"cell_class": ["Exocrine"], "refined_class": ["Exocrine"]}))

    # And it does NOT refuse an ordinary input, or the guard would be a blanket refusal.
    refine._refuse_second_pass(pd.DataFrame({"cell_class": ["Exocrine", "Immune"]}))


# ------------------------------------------- 5. the refinement has to reach the table

def test_the_refinement_reaches_the_table_the_app_reads(vocabulary):
    """The pass has to be joined back, and a run that skips the join has to say so.

    WAS RED, WENT GREEN DURING THIS SESSION on its first half. When this test was written,
    ``refine.py`` wrote a sidecar keyed on ``object_id`` -- the right shape for a pass that must
    not be run twice -- and no builder joined it back, so ``data/development/typing/cells/``,
    the table ``build_app_data.py`` and every senior script read, carried the pre-refinement
    ``cell_class`` for every cell. 136 builder and app sources were searched and none mentioned
    it. ``build_cells.py`` has since grown ``--refine`` and ``_join_refinement``.

    STILL RED ON THE SECOND HALF. ``--refine`` defaults to ``None`` and the join sits behind a
    bare ``if refine is not None:`` with no ``else``, so a build run without the flag writes a
    table with no ``refined_class`` column, no ``refine_pass`` column and nothing in the ledger
    that says the pass was skipped. That table is indistinguishable from one where the pass ran
    and demoted nothing, and it is the one every consumer would then serve. Either the flag is
    required, or the absent case is recorded on the way past.
    """
    needles = ("refined_class", "import refine", "from refine", "refine_pass")
    searched, hits = [], []
    roots = [(TYPING, "*.py"), (SENIOR, "*.py"),
             (REPO / "apps/islet/R", "*.R"), (REPO / "apps/R/ieCore/R", "*.R")]
    for root, glob in roots:
        if not root.exists():
            continue
        for path in sorted(root.rglob(glob)):
            if path.name == "refine.py":
                continue
            searched.append(path)
            text = path.read_text()
            if any(needle in text for needle in needles):
                hits.append(str(path.relative_to(REPO)))

    assert hits, (
        f"nothing consumes the refinement. {len(searched)} builder and app sources were "
        f"searched for {list(needles)} and none of them mentions it, so the refined class is "
        "written to a sidecar and read by no one: the app, the composition builders and the "
        "senior scripts all still serve the pre-refinement cell_class. "
        f"Sources searched under: {[str(r.relative_to(REPO)) for r, _ in roots if r.exists()]}")

    # And the join must not be skippable in silence.
    builder = TYPING / "build_cells.py"
    tree = _tree(builder)
    required = any(
        isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
        and node.func.attr == "add_argument"
        and any(isinstance(a, ast.Constant) and a.value == "--refine" for a in node.args)
        and any(kw.arg == "required" and isinstance(kw.value, ast.Constant)
                and kw.value.value is True for kw in node.keywords)
        for node in ast.walk(tree))
    unhandled = []
    if not required:
        for node in ast.walk(tree):
            if not isinstance(node, ast.If):
                continue
            test_names = {n.id for n in ast.walk(node.test) if isinstance(n, ast.Name)}
            if "refine" not in test_names:
                continue
            body = ast.dump(ast.Module(body=node.body, type_ignores=[]))
            if "_join_refinement" not in body and "refined_class" not in body:
                continue
            if not node.orelse:
                unhandled.append(f"{builder.relative_to(REPO)}:{node.lineno}")
    assert not unhandled, (
        f"the refinement join is optional and its absence is not recorded: {unhandled}. "
        "`--refine` defaults to None, so a run without it writes a cells table with no "
        "refined_class and no refine_pass, which reads exactly like a run where the pass "
        "demoted nothing. Either make the flag required, or write the 'not applied' case onto "
        "the rows and into the ledger so a reader can tell the two apart")


# ------------------------------------ 6. no switch defaults to the superseded behaviour

def test_no_producer_switch_defaults_to_the_superseded_behaviour(vocabulary):
    """A diagnostic switch that can restore old behaviour must default to the current method.

    RED TODAY on one flag. ``subpanel.py``'s three "isolate the change in a diff" switches are
    right: ``--no-top-marker`` and ``--no-islet-chain`` are ``store_true`` (so the new chain is
    on by default) and ``--levels`` defaults to ``per-donor`` with help that says in words which
    choice is the method and which "exists only to reproduce the old table".

    ``--membership`` is the exception. It offers ``parent`` and ``geometry``, defaults to
    ``parent``, and carries NO help at all -- so nothing in the file says which of the two is
    the method, and this guard cannot confirm the default is the current one. That matters
    because ``in_islet`` decides whether a cell is asked the islet question at all: a cell the
    membership does not place inside an islet never reaches the CD99 / elongation chain the
    maintainer just ruled on, and falls through to the complement instead.

    The assertion is deliberately about the help text, not about which value is right: this test
    does not get to rule on the method, only to insist that the file says which value is it.
    """
    findings = []
    for module in ("subpanel.py", "refine.py"):
        for node in ast.walk(_tree(TYPING / module)):
            if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                    and node.func.attr == "add_argument"):
                continue
            flag = next((a.value for a in node.args
                         if isinstance(a, ast.Constant) and isinstance(a.value, str)), None)
            if flag is None or not flag.startswith("--"):
                continue
            kw = {k.arg: k.value for k in node.keywords}
            choices = kw.get("choices")
            if not isinstance(choices, (ast.Tuple, ast.List)):
                continue
            values = [c.value for c in choices.elts if isinstance(c, ast.Constant)]
            default = kw.get("default")
            default = default.value if isinstance(default, ast.Constant) else None
            help_text = kw.get("help")
            help_text = help_text.value if isinstance(help_text, ast.Constant) else (
                "".join(v.value for v in ast.walk(help_text)
                        if isinstance(v, ast.Constant) and isinstance(v.value, str))
                if help_text is not None else "")
            if default is not None and default in help_text:
                continue
            findings.append(
                f"{module} {flag} choices={values} default={default!r} "
                f"help={help_text or '(none)'!r}")

    assert not findings, (
        "a switch that selects between the current method and a superseded behaviour does not "
        "name its default as the method, so nothing in the file says which value is the method "
        "and a run that keeps the old behaviour looks like an ordinary run:\n"
        + "\n".join(f"  {f}" for f in findings))
