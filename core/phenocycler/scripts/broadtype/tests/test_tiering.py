"""The tiering rule, asserted from the code and from the records: identity = type markers only.

Three layers, because each catches what the others cannot:

1. **The loader.**  ``arm_base.load_context`` is the one place an arm receives values, so the
   rule holds for every arm if ``ArmContext.X`` carries the declared identity markers and
   nothing else -- no state marker, no audit column.  Tested on a synthetic uniform parquet.

2. **The arm source (static, heuristic).**  For every ``arms/*.py`` that exists at test time,
   import it by path (it must expose ``ARM`` and ``run``), then walk its AST.  The heuristic:
   * FORBIDDEN identifiers, as a ``Name`` or an ``Attribute``: ``state_markers_for`` and
     ``state_columns`` (the state tier's two entry points);
   * FORBIDDEN column keys -- ``cell_region``, ``islet_num``, ``Parent`` -- as a string
     constant used as a subscript, inside a list/tuple subscript, as a call argument, or as
     an attribute name (``ctx.context.cell_region``); these are QuPath parentage, audit only;
   * anything QuPath-derived as an INPUT: an identifier containing ``qupath`` (case-insensitive,
     e.g. ``pins.QUPATH_MANIFEST``) or a string constant containing ``qupath`` used as a
     subscript or call argument (``"qupath:Nucleus__..."`` columns).  A string mentioning
     qupath in a message or a docstring is not an input and is not flagged;
   * another arm's output as an input: a string constant containing ``records.jsonl`` or
     ``calls__`` used as a subscript / call argument;
   * REQUIRED: the module reads ``.X`` or ``.markers`` off some object (the context).
   Docstrings are skipped.  This is a lexical check on names, not a proof: an arm that
   smuggled a state marker through an alias would pass it, which is why layer 3 exists.

3. **The records.**  For every ``records.jsonl`` under ``data/development/broadtype/arms/``,
   ``identity_markers`` (what actually entered identity) must be a subset of the declared
   identity markers and disjoint from the state markers.  Skipped with a reason when no arm
   has written a record yet.
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

import arm_base
import declarations
import pins
from criteria import CRITERIA

BROADTYPE = pathlib.Path(__file__).resolve().parents[1]
ARMS_DIR = BROADTYPE / "arms"
ARM_MODULES = sorted(p for p in ARMS_DIR.glob("*.py") if not p.name.startswith("_") and not p.name.endswith("_worker.py")) if ARMS_DIR.exists() else []
#: QuPath-named identifiers that read IMAGES or GEOMETRY (channel map, qptiff path, the GeoJSON
#: feature stream that turns a mask label into a UUID): provenance and geometry, not identity.
GEOMETRY_IDENTIFIERS = {"QUPATH_MANIFEST", "manifest_image", "qptiff_path", "geojson_path",
                        "iter_qupath_geojson_features", "qptiff_channels", "qupath_manifest"}
RECORDS = sorted((pins.OUT / "arms").glob("**/records.jsonl")) if (pins.OUT / "arms").exists() else []

FORBIDDEN_IDENTIFIERS = {"state_markers_for", "state_columns"}
FORBIDDEN_COLUMNS = {"cell_region", "islet_num", "Parent"}
FORBIDDEN_INPUT_SUBSTRINGS = ("qupath:", "records.jsonl", "calls__")   # the QuPath GeoJSON *stream* is geometry (mask label -> UUID), not identity; cell_region/islet_num/Parent are the identity-shaped QuPath fields and stay forbidden


# --------------------------------------------------------------------------- 1. the loader
@pytest.fixture
def synthetic_uniform(tmp_path, dec):
    """A uniform (Step-2) table for donor 6436 with identity, state, structural and audit columns."""
    rng = np.random.default_rng(0)
    n = 400
    import panel
    have = panel.donor_markers("6436")
    tab = pd.DataFrame({"donor_id": "6436", "object_id": [f"cell-{i}" for i in range(n)],
                        "qc_analysis_eligible": True, "X_centroid": rng.uniform(0, 1e4, n),
                        "Y_centroid": rng.uniform(0, 1e4, n),
                        "cell_region": rng.choice(["core", "peri", "tissue"], n), "islet_num": rng.integers(-1, 5, n)})
    for m in have:
        tab[m] = rng.normal(size=n)
    p = tmp_path / "uniform.parquet"
    tab.to_parquet(p, index=False)
    fc = pd.DataFrame({"object_id": tab["object_id"].iloc[:250], "region_idx": np.arange(250) % 16,
                       "draw": "blind", "in_expression": True})
    fc.loc[fc.index[:5], "in_expression"] = False                   # GeoJSON-only cells never load
    fp = tmp_path / "frame_cells.parquet"
    fc.to_parquet(fp, index=False)
    return {"uniform": p, "frame_cells": fp, "n": n, "have": have}


def test_load_context_hands_arms_type_markers_only(synthetic_uniform, dec, monkeypatch):
    monkeypatch.setattr(pins, "uniform_parquet", lambda d: synthetic_uniform["uniform"])
    monkeypatch.setattr(arm_base.pd, "read_parquet", arm_base.pd.read_parquet)
    monkeypatch.setattr(arm_base, "frame_cells", lambda d: pd.read_parquet(synthetic_uniform["frame_cells"]))
    identity = arm_base.identity_markers_for("6436")
    state = arm_base.state_markers_for("6436")
    assert set(identity) <= set(dec["identity_markers"]) and set(state) <= set(dec["state_markers"])
    assert not set(identity) & set(state)
    assert identity == [m for m in dec["identity_markers"] if m in synthetic_uniform["have"]]

    ctx = arm_base.load_context("6436", "donor")
    assert ctx.markers == identity and list(ctx.X.columns) == identity
    assert not set(ctx.X.columns) & set(dec["state_markers"])
    assert not set(ctx.X.columns) & set(dec["structural_markers"])
    assert not set(ctx.X.columns) & set(arm_base.CONTEXT_COLUMNS), "audit columns never enter X"
    assert len(ctx.ids) == synthetic_uniform["n"] and ctx.config_id == "primary"
    assert ctx.region_idx is None and ctx.scope == "donor"
    # the audit columns live in ctx.context, aligned to ids, and nowhere else
    assert {"cell_region", "islet_num", "X_centroid", "Y_centroid"} <= set(ctx.context.columns)
    assert len(ctx.context) == len(ctx.ids)

    fr = arm_base.load_context("6436", "frame")
    assert len(fr.ids) == 245, "frame scope = frame cells that are in expression"
    assert fr.region_idx is not None and len(fr.region_idx) == 245
    assert list(fr.X.columns) == identity

    ab = arm_base.load_context("6436", "frame", ablate="CD3e")
    assert "CD3e" not in ab.markers and "CD3e" not in ab.X.columns
    assert ab.markers == [m for m in identity if m != "CD3e"] and ab.config_id == "ablate:CD3e"
    with pytest.raises(KeyError, match="is not an identity marker"):
        arm_base.load_context("6436", "frame", ablate="Ki67")            # a state marker cannot be 'ablated'
    with pytest.raises(KeyError):
        arm_base.load_context("6436", "frame", ablate="cell_region")

    sub = arm_base.load_context("6436", "frame", subsample_frac=CRITERIA["subsample_frac"], seed=1,
                                replicate="subsample:1")
    assert len(sub.ids) == int(round(CRITERIA["subsample_frac"] * 245)) and sub.config_id == "subsample:1"
    assert set(sub.ids) <= set(fr.ids) and len(sub.region_idx) == len(sub.ids)
    again = arm_base.load_context("6436", "frame", subsample_frac=CRITERIA["subsample_frac"], seed=1,
                                  replicate="subsample:1")
    assert np.array_equal(sub.ids, again.ids), "a subsample replicate is a deterministic hash draw"


def test_state_columns_are_computed_after_identity_and_never_thresholded(synthetic_uniform, dec, monkeypatch):
    monkeypatch.setattr(pins, "uniform_parquet", lambda d: synthetic_uniform["uniform"])
    ids = np.array([f"cell-{i}" for i in (5, 2, 9)], dtype=object)
    st = arm_base.state_columns("6436", ids)
    assert st["object_id"].tolist() == list(ids)
    cols = [c for c in st.columns if c != "object_id"]
    assert cols == [f"state_{m}_donor_rank" for m in arm_base.state_markers_for("6436")]
    assert not {c.removeprefix("state_").removesuffix("_donor_rank") for c in cols} & set(dec["identity_markers"])
    vals = st[cols].to_numpy(float)
    assert ((vals > 0) & (vals < 1)).all(), "continuous within-donor ranks, never a positive/negative call"


def test_within_donor_rank_is_in_unit_interval_and_keeps_nan():
    X = pd.DataFrame({"a": [3.0, np.nan, 1.0, 2.0], "b": [np.nan] * 4})
    r = arm_base.within_donor_rank(X)
    assert r["a"].tolist()[1] != r["a"].tolist()[1]                        # NaN stays NaN
    assert r["a"].dropna().tolist() == pytest.approx([2.5 / 3, 0.5 / 3, 1.5 / 3])
    assert r["b"].isna().all()


# --------------------------------------------------------------------------- 2. the arm source
def _string_used_as_input(node: ast.AST, parents: dict) -> bool:
    """Is this string constant a subscript key, an element of a subscript list/tuple, or a call
    argument?  Those are the positions where a column name or a path becomes an INPUT."""
    p = parents.get(node)
    if isinstance(p, (ast.List, ast.Tuple)):
        p = parents.get(p)
    if isinstance(p, ast.Subscript):
        return True
    if isinstance(p, ast.Call):
        return node in p.args or any(kw.value is node for kw in p.keywords)
    if isinstance(p, ast.keyword):
        return True
    if isinstance(p, ast.JoinedStr):
        return _string_used_as_input(p, parents)
    if isinstance(p, ast.BinOp):                                           # pins.OUT / "arms" / ...
        return True
    return False


def _docstring_nodes(tree: ast.AST) -> set:
    out = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            body = node.body
            if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) \
                    and isinstance(body[0].value.value, str):
                out.add(body[0].value)
    return out


def audit_arm_source(path: pathlib.Path) -> dict:
    """The heuristic in the module docstring, applied to one arm module.  Returns the findings."""
    tree = ast.parse(path.read_text(), filename=str(path))
    parents = {}
    for node in ast.walk(tree):
        for child in ast.iter_child_nodes(node):
            parents[child] = node
    docstrings = _docstring_nodes(tree)
    hits, reads_ctx = [], False
    for node in ast.walk(tree):
        if isinstance(node, ast.Name):
            name = node.id
        elif isinstance(node, ast.Attribute):
            name = node.attr
            if name in ("X", "markers"):
                reads_ctx = True
        else:
            name = None
        if name is not None:
            if name in FORBIDDEN_IDENTIFIERS:
                hits.append((node.lineno, f"identifier {name!r} (state tier entry point)"))
            elif name in FORBIDDEN_COLUMNS:
                hits.append((node.lineno, f"identifier {name!r} (QuPath parentage, audit only)"))
            elif "qupath" in name.lower() and name not in GEOMETRY_IDENTIFIERS:
                hits.append((node.lineno, f"identifier {name!r} (QuPath-derived input)"))
            continue
        if isinstance(node, ast.Constant) and isinstance(node.value, str) and node not in docstrings:
            s = node.value
            if not _string_used_as_input(node, parents):
                continue
            if s in FORBIDDEN_COLUMNS:
                hits.append((node.lineno, f"column {s!r} used as an input (QuPath parentage, audit only)"))
            for tok in FORBIDDEN_INPUT_SUBSTRINGS:
                if tok in s.lower():
                    hits.append((node.lineno, f"string {s!r} used as an input ({tok})"))
    return {"hits": sorted(hits), "reads_ctx": reads_ctx}


def _import_by_path(path: pathlib.Path):
    if str(BROADTYPE) not in sys.path:
        sys.path.insert(0, str(BROADTYPE))
    spec = importlib.util.spec_from_file_location(f"broadtype_arm_{path.stem}", path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod   # dataclasses and typing look the module up by name during exec
    spec.loader.exec_module(mod)
    return mod


@pytest.mark.parametrize("path", ARM_MODULES or [None], ids=lambda p: p.stem if p else "no-arms-yet")
def test_no_state_marker_in_identity_inputs(path, dec):
    if path is None:
        pytest.skip(f"no arm module exists yet under {ARMS_DIR.relative_to(pins.REPO)}/")
    mod = _import_by_path(path)
    assert isinstance(getattr(mod, "ARM", None), str) and mod.ARM, f"{path.name}: no ARM name"
    assert mod.ARM.replace("-", "_") == path.stem, f"{path.name}: ARM {mod.ARM!r} does not name the module"
    assert callable(getattr(mod, "run", None)), f"{path.name}: no run(ctx)"
    audit = audit_arm_source(path)
    assert audit["reads_ctx"], f"{path.name}: never reads ctx.X / ctx.markers -- where does its identity input come from?"
    assert not audit["hits"], f"{path.name}: forbidden identity inputs:\n" + "\n".join(
        f"  line {ln}: {why}" for ln, why in audit["hits"])


def test_the_source_heuristic_catches_what_it_claims(tmp_path):
    """The static check is only as good as its own test: plant each forbidden form and one clean arm."""
    bad = tmp_path / "bad_arm.py"
    bad.write_text('''"""An arm that mentions cell_region in its docstring, which is fine."""
import arm_base
ARM = "bad-arm"
def run(ctx):
    a = ctx.context["cell_region"]                     # audit column as input
    b = ctx.context[["X_centroid", "islet_num"]]      # inside a list subscript
    c = ctx.context.get("Parent")                     # as a call argument
    d = arm_base.state_columns(ctx.donor_id, ctx.ids)  # the state tier
    e = tab["qupath:Nucleus__CD3e"]                   # a QuPath-derived column
    f = tab.qupath_classification                     # a QuPath class as an input
    g = pd.read_parquet(pins.OUT / "arms" / "x" / "calls__primary.parquet")   # another arm's calls
    note = "this arm ignores cell_region"             # a message: NOT flagged
    return ctx.X
''')
    audit = audit_arm_source(bad)
    reasons = " | ".join(why for _, why in audit["hits"])
    assert audit["reads_ctx"]
    assert "column 'cell_region'" in reasons and "column 'islet_num'" in reasons and "column 'Parent'" in reasons
    assert "identifier 'state_columns'" in reasons
    assert "qupath:Nucleus__CD3e" in reasons and "identifier 'qupath_classification'" in reasons
    assert "calls__primary.parquet" in reasons
    assert "ignores cell_region" not in reasons, "a message string is not an input"
    assert len(audit["hits"]) == 7

    good = tmp_path / "good_arm.py"
    good.write_text('''"""Uses ctx.X and the centroids only; cell_region is named here in prose only."""
ARM = "good-arm"
def run(ctx):
    xy = ctx.context[["X_centroid", "Y_centroid"]].to_numpy()
    Z = ctx.X[ctx.markers]
    return Z
''')
    audit = audit_arm_source(good)
    assert audit["hits"] == [] and audit["reads_ctx"]


# --------------------------------------------------------------------------- 3. the records
@pytest.mark.parametrize("path", RECORDS or [None],
                         ids=lambda p: str(p.relative_to(pins.OUT)) if p else "no-records-yet")
def test_recorded_identity_markers_are_declared_type_markers_only(path, dec):
    if path is None:
        pytest.skip(f"no records.jsonl exists yet under {(pins.OUT / 'arms').relative_to(pins.REPO)}/")
    identity, state = set(dec["identity_markers"]), set(dec["state_markers"])
    lines = [ln for ln in path.read_text().splitlines() if ln.strip()]
    assert lines, f"{path}: empty record file"
    for ln in lines:
        rec = json.loads(ln)
        used = set(rec["identity_markers"])
        assert used <= identity, f"{path} {rec['arm']}/{rec['config_id']}: undeclared identity input {sorted(used - identity)}"
        assert not used & state, f"{path} {rec['arm']}/{rec['config_id']}: state marker in identity {sorted(used & state)}"
        assert rec["n_called"] + rec["n_ambiguous"] + rec["n_unavailable"] == rec["n_cells"]
        if rec["config_id"].startswith("ablate:"):
            assert rec["config_id"].split(":", 1)[1] not in used, "the withheld marker re-entered identity"


def test_declared_identity_and_state_are_disjoint_registry_kinds(dec):
    kinds = pins.markers_by_kind()
    assert set(dec["identity_markers"]) <= set(kinds[pins.KIND_TYPE])
    assert set(dec["state_markers"]) == set(kinds[pins.KIND_STATE])
    assert not set(dec["identity_markers"]) & set(dec["state_markers"])
    for L, spec in dec["labels"].items():
        assert set(spec["markers"]) <= set(dec["identity_markers"]), L
    assert declarations.load() == dec
