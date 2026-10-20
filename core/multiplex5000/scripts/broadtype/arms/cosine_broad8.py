"""Arm A ``cosine-broad8``: the celltype lane's cosine + requirement gate on the 7 declared labels.

Method (frozen 2026-09-05; README arm table): ``scripts/celltype`` ``rules``/``score`` on
``rules_broad8.csv`` (display names; ``rules.load_rules`` maps them through the frozen
crosswalk), the ``signed`` encoding with OR-group weights shared per DECLARING ROW, the
cohort-wide [aligned background mode, p99.9] scale of ``run_step3.build_scale`` (DAPI ratio
off), ``require_min`` 0.10 with the negative constraints ENFORCED (``none_of``), cosine ranking
inside the gate, ``score.call_from_scores``.  The rows of the table ARE the labels, so no
collapse is needed; ``score`` is the winning cosine (ceiling-normalised, as ``run_step3`` scores
it); a cell meeting no row is ``Ambiguous(no_reference_requirements_met)``.

What the arm does with what ``ArmContext`` hands it -- each a recorded decision, none a tuning:

* ``ctx.markers`` is the declared identity set minus the ablated marker.  A marker the table
  names and the context lacks (the ablated one, or one this donor never acquired) is dropped
  from the reference matrix AND from the requirement gate, positives and negatives alike.  The
  ablation asks whether a label still lands on the withheld marker's tail through its OTHER
  markers, so the gate must neither demand nor forbid a marker nobody can see.  A label left
  with no positive marker cannot be scored: it leaves the candidate set, and a cell that then
  meets no remaining label is ``Unavailable(ablation_removes_required_input)`` rather than
  ``Ambiguous`` -- the removed label was a possible answer that cannot be evaluated, and
  Ambiguous is never a parking bucket.  With every label carrying >= 2 declared markers no
  single ablation empties one on this table; the path exists for a donor missing a marker.
* The nucleus/DAPI ratio is NOT applied (``Xs_nuc=None``).  ``score.nuclear_only_references``
  admits a reference only when EVERY marker it constrains is nucleus-measured, and no row of
  ``rules_broad8`` is (each mixes whole-cell markers).  Rather than carry an inert branch, the
  arm asserts that and fails loud if a table ever makes a row qualify, so it can never score a
  reference on the wrong quantity silently.
* ``qc_analysis_eligible`` is the geometry-QC flag carried on the context (not a QuPath class).
  ``run_step3`` scores none of the ineligible cells; there they land in Ambiguous by a quirk of
  the gate, here they are ``Unavailable(qc_analysis_ineligible)``.  A cell with a NaN in any
  scored marker is ``Unavailable(missing_measurement)``: a missing measurement is not zero.

The scale is cohort-wide by design -- per donor would re-impose the distribution matching Step 2
avoided -- and reads all 20 donors' aligned parquets once; it is cached at
``data/development/broadtype/scratch/cosine_scale.parquet`` with its marker list and a key over
its inputs (marker list, donors, ``build_scale`` parameters, the expression content id, the
Step-2 ``shifts.json`` and the aligned files' size/mtime), so repeated runs and the sibling arm
``cosine-37-collapsed`` reuse it and a Step-2 re-run invalidates it.

Deterministic; ``ctx.seed`` is ignored.  Uses ``ctx.X`` and the QC flag only -- no state marker,
no ``cell_region``/``islet_num``, no other arm's calls.

    cd core/multiplex5000 && ../../.venv/bin/python scripts/broadtype/arms/cosine_broad8.py --donor 6436 --scope frame
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import sys
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

HERE = pathlib.Path(__file__).resolve().parent
if str(HERE.parent) not in sys.path:
    sys.path.insert(0, str(HERE.parent))            # scripts/broadtype
import arm_base  # noqa: E402
import pins  # noqa: E402
import records  # noqa: E402
import panel  # noqa: E402  (scripts/celltype, on sys.path via pins)
import rules  # noqa: E402
import run_step3 as S3  # noqa: E402
import score as SC  # noqa: E402

ARM = "cosine-broad8"
RULES_CSV = HERE.parent / "rules_broad8.csv"

#: Frozen configuration.  A change to any of these after a result is seen is a new arm name.
ENCODING = "signed"                 # run_step3.ENCODING: resolves every reference; zero cannot
SPLIT_OR_GROUPS = True              # unit weight shared per (declaring row, token) group
REQUIRE_MIN = 0.10                  # run_step3.REQUIRE_MIN: distance above the marker's own background
CEILING = True                      # divide by each reference's attainable maximum (run_step3 default)
NUCLEAR_RATIO = False               # see the module docstring: no nucleus-only reference exists
SCALE_N_PER_DONOR = 200_000         # build_scale default
SCALE_SEED = 0                      # build_scale default
SCALE_CACHE = pins.OUT / "scratch" / "cosine_scale.parquet"
SCALE_META_KEY = b"broadtype_cosine_scale"

REASON_NO_RULE = "no_reference_requirements_met"
REASON_ABLATED = "ablation_removes_required_input"
REASON_QC = "qc_analysis_ineligible"
REASON_NAN = "missing_measurement"


def _sha(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


PARAMS = {
    "rules": str(RULES_CSV.relative_to(pins.REPO)),
    "rules_sha256": _sha(RULES_CSV),
    "encoding": ENCODING,
    "split_or_groups": SPLIT_OR_GROUPS,
    "or_group_key": "(declaring_row, token)",
    "require_min": REQUIRE_MIN,
    "negatives_enforced": True,
    "ceiling_normalised": CEILING,
    "nuclear_ratio": NUCLEAR_RATIO,
    "scale": {"builder": "scripts/celltype/run_step3.build_scale", "cohort_wide": True,
              "n_per_donor": SCALE_N_PER_DONOR, "seed": SCALE_SEED, "top_quantile": SC.SCALE_Q,
              "floor": "aligned background mode", "cache": str(SCALE_CACHE.relative_to(pins.REPO))},
    "absent_marker_policy": "dropped from the reference matrix and the gate; a label with no positive left is unavailable",
    "qc_ineligible": REASON_QC,
    "nan_input": REASON_NAN,
}


# --------------------------------------------------------------------------- references
@dataclass
class References:
    """The scorable part of a rule table on exactly the markers a context carries."""
    names: list                                  # scorable references (matrix columns)
    markers: list                                # markers some scorable reference uses (rows)
    R: np.ndarray                                # markers x names, ENCODING weights
    reqs: dict                                   # name -> (all_of, any_groups, none_of), restricted
    dropped_markers: list = field(default_factory=list)   # named by the table, absent from the context
    unavailable: dict = field(default_factory=dict)       # reference -> reason (no positive marker left)

    def frame(self) -> pd.DataFrame:
        return pd.DataFrame(self.R, index=self.markers, columns=self.names)


def restrict(table: pd.DataFrame, markers: list) -> References:
    """Reference matrix and gate for ``markers`` only, recording what could not be kept.

    Absent markers leave both the matrix (``reference_matrix(strict=False)`` drops the row and
    reports it) and the gate (filtered out of ``all_of``/``none_of``; an OR-group loses the
    member and an emptied group is dropped rather than made unsatisfiable).  A reference whose
    positive evidence is gone entirely is unavailable, not scored against nothing.  Markers no
    remaining reference uses leave the cell vector too, as ``run_step3`` does, so an unused
    marker cannot inflate the norm the cosine divides by.
    """
    have = set(markers)
    flat = rules.flatten(table)
    Rm, dropped = rules.reference_matrix(flat, list(markers), ENCODING, strict=False,
                                         split_or_groups=SPLIT_OR_GROUPS)
    reqs, unavailable = {}, {}
    for name, (all_of, groups, none_of) in rules.requirements_with_negatives(table).items():
        a = [m for m in all_of if m in have]
        g = [[m for m in grp if m in have] for grp in groups]
        g = [grp for grp in g if grp]
        n = [m for m in none_of if m in have]
        if not a and not g:
            unavailable[name] = REASON_ABLATED
            continue
        if not (Rm[name] > 0).any():                       # the matrix must agree with the gate
            raise RuntimeError(f"{name}: gate keeps positives {a}/{g} but the matrix has none")
        reqs[name] = (a, g, n)
    names = [c for c in Rm.columns if c in reqs]
    Rm = Rm[names]
    Rm = Rm.loc[(Rm != 0).any(axis=1)]
    return References(names=names, markers=list(Rm.index), R=Rm.to_numpy("float32"), reqs=reqs,
                      dropped_markers=sorted(dropped), unavailable=unavailable)


# --------------------------------------------------------------------------- the cohort scale
def _scale_universe() -> list:
    """Every marker the aligned files can hold, so one build serves both cosine arms."""
    return sorted(panel.markers(include_variant=True))


def _scale_key(markers: list, donors: list) -> str:
    stats = {}
    for d in donors:
        p = pins.uniform_parquet(d)
        st = p.stat()
        stats[d] = [st.st_size, st.st_mtime_ns]
    payload = {
        "markers": list(markers), "donors": list(donors),
        "n_per_donor": SCALE_N_PER_DONOR, "seed": SCALE_SEED, "ratio": NUCLEAR_RATIO, "q": SC.SCALE_Q,
        "expression_content_id": pins.EXPRESSION_CONTENT_ID,
        "shifts_sha256": _sha(pins.SHIFTS) if pins.SHIFTS.exists() else "",
        "uniform_files": stats,
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def _load_cached_scale(donors: list) -> pd.DataFrame | None:
    if not SCALE_CACHE.exists():
        return None
    t = pq.read_table(SCALE_CACHE)
    meta = json.loads((t.schema.metadata or {}).get(SCALE_META_KEY, b"{}"))
    if not meta or meta.get("sha256") != _scale_key(meta.get("markers", []), donors):
        return None                                        # inputs moved: rebuild
    return t.to_pandas().set_index("marker")


def cohort_scale(markers: list) -> pd.DataFrame:
    """floor / top / span per marker, cohort-wide, cached with its marker list and input key."""
    donors = panel.cohort()["donor_id"].tolist()
    scale = _load_cached_scale(donors)
    if scale is None or not set(markers) <= set(scale.index):
        universe = sorted(set(_scale_universe()) | set(markers))
        scale = S3.build_scale(donors, universe, n_per_donor=SCALE_N_PER_DONOR, seed=SCALE_SEED,
                               ratio=NUCLEAR_RATIO)
        meta = {"markers": universe, "sha256": _scale_key(universe, donors), "donors": donors,
                "n_per_donor": SCALE_N_PER_DONOR, "seed": SCALE_SEED, "ratio": NUCLEAR_RATIO,
                "top_quantile": SC.SCALE_Q, "builder": "scripts/celltype/run_step3.build_scale"}
        SCALE_CACHE.parent.mkdir(parents=True, exist_ok=True)
        tab = pa.Table.from_pandas(scale.reset_index())
        pq.write_table(tab.replace_schema_metadata({**(tab.schema.metadata or {}),
                                                    SCALE_META_KEY: json.dumps(meta).encode()}), SCALE_CACHE)
    missing = [m for m in markers if m not in scale.index]
    if missing:
        raise KeyError(f"no cohort scale for {missing}: the aligned files hold no values for them")
    return scale.loc[list(markers)]


# --------------------------------------------------------------------------- calling
def call_cells(ctx: arm_base.ArmContext, refs: References, scale: pd.DataFrame, *,
               X: pd.DataFrame | None = None, collapse: dict | None = None,
               label_unavailable: bool = False) -> pd.DataFrame:
    """Score, gate and call every cell in ``ctx.ids``; returns ``records.CALL_COLUMNS`` (+ margin).

    ``collapse`` maps a reference name to ``(label, reason)`` for a table whose rows are not the
    labels (arm A'); without it a row is its own label.  ``label_unavailable`` says a declared
    label could not be scored at all on this input, which turns a no-rule cell from Ambiguous
    into ``Unavailable(ablation_removes_required_input)``.
    """
    X = ctx.X if X is None else X
    n = len(ctx.ids)
    label = np.full(n, records.UNAVAILABLE, dtype=object)
    reason = np.full(n, "", dtype=object)
    sc = np.full(n, np.nan, dtype=float)
    margin = np.full(n, np.nan, dtype=float)
    elig = (ctx.context["qc_analysis_eligible"].to_numpy(bool)
            if "qc_analysis_eligible" in ctx.context.columns else np.ones(n, dtype=bool))
    if refs.markers:
        Xv = X[refs.markers].to_numpy("float32")
        finite = np.isfinite(Xv).all(axis=1)
    else:
        Xv = np.zeros((n, 0), dtype="float32")
        finite = np.ones(n, dtype=bool)
    reason[~elig] = REASON_QC
    reason[elig & ~finite] = REASON_NAN
    ok = elig & finite
    if not refs.names:
        reason[ok] = REASON_ABLATED
    elif ok.any():
        nuc = SC.nuclear_only_references(refs.frame(), panel.NUCLEUS_MEASURED)
        if nuc.any():
            raise RuntimeError(f"{ARM}: nucleus-only references {list(np.array(refs.names)[nuc])} -- the "
                               "DAPI ratio is deliberately not implemented in this arm; a new table "
                               "that needs it is a new arm")
        Xs = SC.apply_scale(Xv[ok], refs.markers, scale)
        S, allowed = SC.score_and_gate(Xs, None, refs.R, refs.markers, refs.reqs, refs.names,
                                       REQUIRE_MIN, nuc, CEILING)
        out = SC.call_from_scores(S, refs.names, allowed=allowed)
        call = out["call"].to_numpy(dtype=object)
        if (call == "Unscored").any():                     # a hot marker gives a norm; cannot happen
            raise RuntimeError(f"{ARM}: {int((call == 'Unscored').sum())} cells passed a gate with no angle")
        idx = np.flatnonzero(ok)
        amb = call == "Ambiguous"
        sc[idx] = out["score"].to_numpy(float)
        margin[idx] = out["margin"].to_numpy(float)
        if collapse is None:
            label[idx] = call
        else:
            mapped = [collapse.get(c, (c, "")) for c in call]
            label[idx] = [m[0] for m in mapped]
            reason[idx] = [m[1] for m in mapped]
        label[idx[amb]] = records.UNAVAILABLE if label_unavailable else records.AMBIGUOUS
        reason[idx[amb]] = REASON_ABLATED if label_unavailable else REASON_NO_RULE
    return pd.DataFrame({"object_id": pd.Series(ctx.ids).astype(str), "label": label,
                         "score": sc, "reason": reason, "margin": margin})


def summarise(calls: pd.DataFrame, refs: References) -> str:
    counts = calls["label"].value_counts().to_dict()
    return (f"labels {counts}; dropped markers {refs.dropped_markers or 'none'}; "
            f"unavailable references {sorted(refs.unavailable) or 'none'}")


def run(ctx: arm_base.ArmContext) -> pd.DataFrame:
    table = rules.load_rules(pathlib.Path(ctx.extra.get("rules_csv", RULES_CSV)))
    refs = restrict(table, list(ctx.markers))
    scale = ctx.extra.get("scale")
    if scale is None:
        scale = cohort_scale(refs.markers)
    calls = call_cells(ctx, refs, scale, label_unavailable=bool(refs.unavailable))
    print(f"{ARM} {ctx.donor_id} {ctx.config_id}: {summarise(calls, refs)}", flush=True)
    return calls


if __name__ == "__main__":
    raise SystemExit(arm_base.cli(ARM, run, env=sys.executable, default_params=PARAMS))
