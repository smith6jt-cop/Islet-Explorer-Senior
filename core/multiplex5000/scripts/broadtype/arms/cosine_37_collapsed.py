"""Arm A' ``cosine-37-collapsed``: the 37-reference lane's calls, collapsed to the 7 labels.

The continuity row.  It runs EXACTLY the machinery of ``cosine_broad8`` -- ``signed`` encoding,
OR-group weights per declaring row, cohort-wide scale, ``require_min`` 0.10 with negatives
enforced, ceiling-normalised cosine inside the gate -- on the maintainer's 37-reference table
``phenotype_workflow.csv`` (with the 2026-09-05 fixes to ``rules``/``score``), then maps each
winning reference to a broad label by walking the table's own hierarchy:

    subtree of Epi-endocrine                                  -> Endocrine
    every other Epithelial descendant, and Epithelial itself   -> Epithelial_nonendocrine
    Imm-T-cell / Imm-NK-cell subtrees                          -> T_NK
    Imm-B-cell / Imm-plasma-cell                               -> B_Plasma
    Imm-macrophage / Imm-dendritic-cell / Imm-neutrophil       -> Myeloid
    Endothelial subtree                                        -> Vascular
    Mesenchymal subtree                                        -> Stromal
    Immune, undivided                                          -> Ambiguous(immune_branch_undivided)
    Neural and its descendants                                 -> Unavailable(label_not_declared)
    AMBIGUOUS-CD11b-pos-epithelial                             -> Ambiguous(hybrid_call)

``collapse_map`` derives the map from the table and refuses a table with an unmapped row, so a
new reference is placed deliberately or the arm stops.  ``score`` is the winning reference's
cosine; a cell meeting no reference is ``Ambiguous(no_reference_requirements_met)``.

**Flagged: this arm's identity input is the 37 table's marker set, which includes registry
``state`` markers** (Ki67, PCNA, HLA-DR, Granzyme B, ...), because that is what the 37 lane
scores -- its leaves are state-gated (``T-exhausted``, ``Exo-acinar-ductal-proliferating``).  The
tiering rule of the protocol is therefore violated on purpose, once, for continuity with the
numbers the lane has been producing; ``params["flag"]`` says so and the record's
``identity_markers`` lists the markers actually used (the module sets ``ctx.markers`` to that
list before ``arm_base.write`` reads it, because the contract's ``cli`` leaves no other hook).
The extra columns are read from this donor's aligned parquet for ``ctx.ids`` only; nothing
else about the context is touched, and ``cell_region``/``islet_num`` are never read.

Marker availability follows arm A, not ``run_step3``: a table marker this donor never acquired
(``b_Catenin1`` on every donor, ``CD38`` on 7) or the ablated marker is dropped from the matrix
and the gate rather than making every reference that names it unscorable.  That is the
"same machinery" the protocol asks for; where it differs from the 37 lane's own driver (which
leaves such a reference unscored on that donor) is recorded here.  A reference left with no
positive marker is unavailable; a LABEL is unavailable only when every reference collapsing
to it is, and only then does a no-rule cell become ``Unavailable(ablation_removes_required_input)``.

The nucleus/DAPI ratio is not applied, as in arm A: ``rules.flatten`` composes every leaf with
its ancestors, so no reference of this table constrains only nucleus-measured markers
(``tests/test_usability.py`` asserts it) and the arm fails loud if one ever does.

Deterministic; ``ctx.seed`` is ignored.

    cd core/multiplex5000 && ../../.venv/bin/python scripts/broadtype/arms/cosine_37_collapsed.py --donor 6436 --scope frame
"""
from __future__ import annotations

import pathlib
import sys

import pandas as pd
import pyarrow.parquet as pq

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE.parent) not in sys.path:
    sys.path.insert(0, str(HERE.parent))            # scripts/broadtype
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import arm_base  # noqa: E402
import declarations  # noqa: E402
import pins  # noqa: E402
import records  # noqa: E402
import panel  # noqa: E402  (scripts/celltype, on sys.path via pins)
import rules  # noqa: E402
import cosine_broad8 as A  # noqa: E402

ARM = "cosine-37-collapsed"
RULES_CSV = pins.RULES_37

#: Subtree root -> broad label.  Order matters: a more specific root claims its subtree before
#: an enclosing one (Epi-endocrine before Epithelial; the immune branches before Immune).
SUBTREE_LABEL = (
    ("Epi-endocrine", "Endocrine"),
    ("Imm-T-cell", "T_NK"), ("Imm-NK-cell", "T_NK"),
    ("Imm-B-cell", "B_Plasma"), ("Imm-plasma-cell", "B_Plasma"),
    ("Imm-macrophage", "Myeloid"), ("Imm-dendritic-cell", "Myeloid"), ("Imm-neutrophil", "Myeloid"),
    ("Endothelial", "Vascular"),
    ("Mesenchymal", "Stromal"),
    ("Epithelial", "Epithelial_nonendocrine"),
)
#: Rows (and, for Neural, a subtree) that collapse to a non-label with a recorded reason.
SPECIAL = {
    "Immune": (records.AMBIGUOUS, "immune_branch_undivided"),
    "AMBIGUOUS-CD11b-pos-epithelial": (records.AMBIGUOUS, "hybrid_call"),
}
NEURAL_ROOT = "Neural"
NEURAL = (records.UNAVAILABLE, "label_not_declared")

PARAMS = {
    **A.PARAMS,
    "rules": str(RULES_CSV.relative_to(pins.REPO)),
    "rules_sha256": A._sha(RULES_CSV),
    "flag": "inherits state-gated rule table",
    "identity_input": "all markers of the 37-reference table this donor carries, state markers included",
    "collapse": "rule-table hierarchy walk (SUBTREE_LABEL, SPECIAL, Neural -> Unavailable)",
    "absent_marker_policy": ("as arm A: dropped from the matrix and the gate (differs from run_step3, which "
                             "leaves a reference naming an unacquired marker unscored on that donor)"),
}


# --------------------------------------------------------------------------- the collapse
def _children(table: pd.DataFrame) -> dict:
    kids: dict = {}
    for p, c in zip(table["parent"].astype(str), table["phenotype"].astype(str)):
        kids.setdefault(p, []).append(c)
    return kids


def _subtree(root: str, kids: dict) -> list:
    out, stack = [], [root]
    while stack:
        n = stack.pop()
        out.append(n)
        stack.extend(kids.get(n, []))
    return out


def collapse_map(table: pd.DataFrame) -> dict:
    """reference -> (label, reason), derived from the table's parent column; fails on a gap."""
    kids = _children(table)
    names = set(table["phenotype"].astype(str))
    out: dict = {}
    for ph, val in SPECIAL.items():
        if ph in names:
            out[ph] = val
    if NEURAL_ROOT in names:
        for ph in _subtree(NEURAL_ROOT, kids):
            out.setdefault(ph, NEURAL)
    for root, label in SUBTREE_LABEL:
        if root not in names:
            continue
        for ph in _subtree(root, kids):
            out.setdefault(ph, (label, ""))
    unmapped = sorted(names - set(out))
    if unmapped:
        raise RuntimeError(f"{ARM}: references with no broad label: {unmapped} -- place them deliberately")
    declared = set(declarations.load()["labels"])
    bad = sorted({v[0] for v in out.values()} - declared - {records.AMBIGUOUS, records.UNAVAILABLE})
    if bad:
        raise RuntimeError(f"{ARM}: collapse targets undeclared labels {bad}")
    return out


def unavailable_labels(refs: A.References, collapse: dict) -> list:
    """Declared labels none of whose references survived the marker restriction."""
    reachable = {collapse[n][0] for n in refs.names if n in collapse}
    lost = {collapse[n][0] for n in refs.unavailable if n in collapse}
    return sorted(L for L in lost if L not in reachable
                  and L not in (records.AMBIGUOUS, records.UNAVAILABLE))


# --------------------------------------------------------------------------- the input
def markers_used(ctx: arm_base.ArmContext, table: pd.DataFrame) -> list:
    """Every table marker this donor carries and the aligned file holds, minus the ablated one."""
    cols = [c for c in table.columns if c not in ("parent", "phenotype")]
    if isinstance(ctx.extra.get("X_full"), pd.DataFrame):
        on_disk = set(ctx.extra["X_full"].columns)
        carried = on_disk
    else:
        on_disk = {f.name for f in pq.ParquetFile(pins.uniform_parquet(ctx.donor_id)).schema_arrow}
        carried = set(panel.donor_markers(ctx.donor_id))
    return [m for m in cols if m in on_disk and m in carried and m != ctx.ablated]


def values_for(ctx: arm_base.ArmContext, used: list) -> pd.DataFrame:
    """``ctx.X`` widened with the table's other markers for ``ctx.ids``, in ``ctx.ids`` order."""
    if isinstance(ctx.extra.get("X_full"), pd.DataFrame):
        full = ctx.extra["X_full"]
        return full[used].reset_index(drop=True)
    extra = [m for m in used if m not in ctx.X.columns]
    X = ctx.X.copy()
    if extra:
        tab = pd.read_parquet(pins.uniform_parquet(ctx.donor_id), columns=["object_id", *extra])
        tab["object_id"] = tab["object_id"].astype(str)
        add = tab.set_index("object_id").reindex(pd.Index(ctx.ids.astype(str)))
        for m in extra:
            X[m] = add[m].to_numpy(float)
    return X[used]


def run(ctx: arm_base.ArmContext) -> pd.DataFrame:
    table = rules.load_rules(pathlib.Path(ctx.extra.get("rules_csv", RULES_CSV)))
    collapse = collapse_map(table)
    used = markers_used(ctx, table)
    X = values_for(ctx, used)
    refs = A.restrict(table, used)
    scale = ctx.extra.get("scale")
    if scale is None:
        scale = A.cohort_scale(refs.markers)
    lost = unavailable_labels(refs, collapse)
    calls = A.call_cells(ctx, refs, scale, X=X, collapse=collapse, label_unavailable=bool(lost))
    # The record must say what actually entered identity (state markers included, flagged):
    # arm_base.cli hands ctx straight to write(), which reads ctx.markers, so it is set here.
    ctx.markers = list(refs.markers)
    ctx.extra["unavailable_labels"] = lost
    print(f"{ARM} {ctx.donor_id} {ctx.config_id}: {A.summarise(calls, refs)}; "
          f"unavailable labels {lost or 'none'}; {len(refs.markers)} markers used", flush=True)
    return calls


if __name__ == "__main__":
    raise SystemExit(arm_base.cli(ARM, run, env=sys.executable, default_params=PARAMS))
