"""The concordance gate: does the model agree with this workflow's own measurements?

Pre-registered 2026-08-24, before the confound factorial runs.  It answers one question
with rank statistics and nothing else: for each marker the marker registry already calls
a *type* marker of a broad compartment, does the model's probability for that compartment
track the measured intensity of that marker, better than it tracks any other compartment?

These are acceptance bars on **rank agreement between the model and the run's own REDSEA
expression**.  They are not marker thresholds, they derive no cutoff, they need no
target/reference pair, and nothing here may be tuned toward an expected cell-type
frequency, cross-donor uniformity or disease ordering.

The anchors are DERIVED, not chosen (see :func:`derive_anchors`): every registry marker
with ``kind == "type"`` and a non-null ``broad_compartment`` that resolves into the
DeepCellTypes vocabulary and is acquired on every cohort donor.  Hand-picking anchors after
seeing a result would be the tuning this whole lane exists to avoid, so the rule is code
and ``anchors.json`` is its frozen output.

Two calibration ends are reported beside every statistic, because a number with neither
is unreadable:

* the **permutation null** -- the same statistic on a probability matrix shuffled across
  cells, which is what "no information" looks like through this join;
* the **measured ceiling** -- ``plane_identity_*.txt``'s median self-rho of 0.920, what
  agreement looks like when the two quantities really are the same thing measured twice.
"""
from __future__ import annotations

import json
import pathlib
import sys
from dataclasses import dataclass, field

import numpy as np
from scipy.stats import rankdata, spearmanr

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "common"))  # gates/tiff_window/window_mask/label_space live in scripts/common (Track E.3)

from gates import Verdict, assess_anchor_concordance, require_evidence  # noqa: E402
from label_space import load_label_space  # noqa: E402

__all__ = ["derive_anchors", "GateInputs", "GateResult", "run_gate", "CRITERIA", "CEILING"]

ANCHORS = HERE / "anchors.json"
CROSSWALK = HERE / "broad_class_crosswalk.json"

#: Measured on donor 6539, 2026-08-21: the median Spearman rho between a plane's own
#: mean intensity over the mask and that marker's column in the run's expression parquet,
#: through this same object_id join.  The practical ceiling for any agreement statistic
#: computed this way.  See data/development/deepcell_types_review/6539/plane_identity_6539.txt
CEILING = 0.920

#: Pre-registered thresholds.  Fixed 2026-08-24, before the factorial ran.  Counts are
#: expressed as FRACTIONS of the derived anchor set, because that set is derived rather
#: than enumerated and its size is a property of the panel, not of this file.
CRITERIA = {
    "sign_min_rho": 0.0,           # (1) every anchor must be positive against its compartment
    "specificity_rank1_frac": 11 / 15,   # (2) this fraction must rank its compartment 1st ...
    "specificity_max_rank": 3,           #     ... and none may rank worse than this
    "auroc_immune": 0.75,          # (3) the compartment this study consumes
    "auroc_other": 0.70,
    "coherence_frac": 11 / 15,     # (4) the argmax must agree with the model's own probability
    "impossible_mass_mean": 0.10,  # (5) mean probability on classes the tissue cannot hold
    "impossible_mass_outweighed": 0.10,  #     fraction of cells where it beats the survivor
    "confidence_frac": 11 / 15,    # (7) does NOT block release; licenses stratification
    "min_obs": 50,                 # below this an anchor is unjudged, never passed
    "null_draws": 200,
    "null_quantile": 0.95,
}


# --------------------------------------------------------------------- anchors
def derive_anchors(repo_root: pathlib.Path | None = None) -> dict:
    """Every registry type-marker with a compartment that the model can actually see.

    Mechanical on purpose.  The inputs -- ``marker_registry.json``'s ``kind`` and
    ``broad_compartment``, and the DCT crosswalk -- all predate any inference output.
    """
    root = pathlib.Path(repo_root) if repo_root else HERE.parents[1]
    from deepcell_types.utils import resolve_supported_marker
    from phenocycler.cohort import ensure_eligible_donors

    reg = json.loads((root / "phenocycler" / "marker_registry.json").read_text())["markers"]
    reg = reg if isinstance(reg, list) else [dict(name=k, **v) for k, v in reg.items()]
    xw = json.loads((HERE / "deepcell_types_crosswalk.json").read_text())["entries"]
    man = json.loads((root / "data" / "qupath_manifest.json").read_text())
    donors = ensure_eligible_donors(
        [str(im["donor_id"]).lstrip("0") for im in man["images"]],
        context="DCT anchor cohort",
    )

    per_donor, marker_to_dct = [], {}
    for im in man["images"]:
        got = set()
        for c in im["channel_map"]:
            nm = c["channel_name"]
            t = xw[nm]["target"] if nm in xw else resolve_supported_marker(nm)
            if t:
                got.add(t)
                marker_to_dct[c["marker"]] = t
        per_donor.append(got)
    intersection = set.intersection(*per_donor)

    anchors, dropped = [], []
    for m in reg:
        if m["kind"] != "type" or not m.get("broad_compartment"):
            continue
        dct = marker_to_dct.get(m["name"])
        if dct is None or dct not in intersection:
            dropped.append({"marker": m["name"], "compartment": m["broad_compartment"],
                            "reason": "no DeepCellTypes vocabulary entry" if dct is None
                                      else "not acquired on every donor"})
            continue
        anchors.append({"marker": m["name"], "dct_marker": dct,
                        "compartment": m["broad_compartment"]})
    anchors.sort(key=lambda a: (a["compartment"], a["marker"]))
    return {
        "schema_version": 1,
        "method_version": "dct-anchors-2026-08-24",
        "rule": (
            "Every marker in phenocycler/marker_registry.json with kind == 'type' and a "
            "non-null broad_compartment, that resolves into the DeepCellTypes 278-marker "
            f"vocabulary and is acquired on all {len(set(donors))} cohort donors in "
            f"qupath_manifest.json (the {len(intersection)}-marker cross-donor intersection). "
            "Derived, never hand-picked: choosing anchors after seeing a "
            "result would be the tuning this lane exists to avoid."
        ),
        "reads": "the run's REDSEA expression parquet, one authoritative value per marker",
        "anchors": anchors,
        "excluded": dropped,
    }


def load_anchors(path: pathlib.Path = ANCHORS) -> list:
    d = json.loads(pathlib.Path(path).read_text())
    return [(a["marker"], a["compartment"]) for a in d["anchors"]]


# ----------------------------------------------------------------- statistics
def _auroc(score, positive) -> float:
    positive = np.asarray(positive, dtype=bool)
    n1, n0 = int(positive.sum()), int((~positive).sum())
    if n1 == 0 or n0 == 0:
        return np.nan
    r = rankdata(score)
    return float((r[positive].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def _decile_contrast(x):
    """Top vs bottom decile of a measured marker -- an evaluation contrast, not a call.

    Deciles are used to *score* the model, never to label a cell.  Nothing downstream
    consumes this split and no threshold is derived from it.
    """
    ok = np.isfinite(x)
    if ok.sum() < 20:
        return np.zeros(len(x), bool), np.zeros(len(x), bool)
    lo, hi = np.quantile(x[ok], [0.10, 0.90])
    return ok & (x >= hi), ok & (x <= lo)


@dataclass(frozen=True)
class GateInputs:
    """One arm on one region: the model's probabilities beside our own measurements."""
    probabilities: np.ndarray    # (n, 51) softmax in the model's own ct2idx column order
    classes: tuple               # those 51 column names
    measured: dict               # registry marker name -> (n,) measured intensity, NaN = unavailable
    label: str = ""              # donor/region/arm, for reporting


@dataclass
class GateResult:
    label: str
    restricted: bool
    anchor_verdict: Verdict
    rho: dict = field(default_factory=dict)          # marker -> {compartment: rho}
    auroc: dict = field(default_factory=dict)        # marker -> AUROC of its compartment
    rank: dict = field(default_factory=dict)
    n_obs: dict = field(default_factory=dict)
    call_share: dict = field(default_factory=dict)
    coherence: dict = field(default_factory=dict)    # marker -> (share called in top decile, bottom)
    impossible_mass_mean: float = np.nan
    impossible_outweighed: float = np.nan
    confidence_auroc: dict = field(default_factory=dict)   # marker -> (low, high)
    null: dict = field(default_factory=dict)
    criteria: dict = field(default_factory=dict)     # name -> (passed, detail)

    @property
    def passed(self) -> bool:
        blocking = [v[0] for k, v in self.criteria.items() if k != "confidence_informativeness"]
        return bool(blocking) and all(blocking)


def run_gate(inp: GateInputs, *, restricted: bool = True, anchors=None,
             criteria: dict | None = None, rng_seed: int = 0) -> GateResult:
    """Score one arm.  Every criterion fails closed on absent or under-powered evidence."""
    C = dict(CRITERIA if criteria is None else criteria)
    anchors = anchors if anchors is not None else load_anchors()
    ls = load_label_space(inp.classes)

    res = ls.restrict(inp.probabilities, inp.classes)
    if restricted:
        P, cols = res.probabilities, ls.kept
        call = res.restricted_call
    else:
        P, cols = np.asarray(inp.probabilities, float), ls.classes
        call = res.unrestricted_call
    groups_map = ls.broad_groups(restricted=restricted)
    groups = sorted(groups_map)
    gidx = {g: [cols.index(c) for c in groups_map[g]] for g in groups}

    # Cells the restriction could not call carry no compartment probability; excluding
    # them here is the same fail-closed rule as everywhere else, and it is recorded.
    usable = ~np.isnan(P).any(axis=1) if restricted else np.ones(len(P), bool)
    G = np.column_stack([np.where(usable, P[:, gidx[g]].sum(axis=1), np.nan) for g in groups])

    out = GateResult(label=inp.label, restricted=restricted,
                     anchor_verdict=Verdict(False, [], "not assessed"))

    S = np.full((len(anchors), len(groups)), np.nan)
    counts = np.zeros(len(anchors), dtype=int)
    for i, (marker, comp) in enumerate(anchors):
        x = np.asarray(inp.measured.get(marker, np.full(len(P), np.nan)), dtype=float)
        for j, g in enumerate(groups):
            ok = np.isfinite(x) & np.isfinite(G[:, j])
            if ok.sum() >= C["min_obs"]:
                S[i, j] = spearmanr(x[ok], G[ok, j]).statistic
        counts[i] = int((np.isfinite(x) & usable).sum())
        out.rho[marker] = {g: S[i, j] for j, g in enumerate(groups)}
        out.n_obs[marker] = counts[i]
        if comp in groups:
            j = groups.index(comp)
            finite = np.isfinite(S[i])
            out.rank[marker] = int(1 + (S[i][finite] > S[i, j]).sum()) if np.isfinite(S[i, j]) else len(groups)
            hi, lo = _decile_contrast(x)
            sel = (hi | lo) & usable
            out.auroc[marker] = _auroc(G[sel, j], hi[sel]) if sel.sum() >= C["min_obs"] else np.nan

    out.anchor_verdict = assess_anchor_concordance(
        anchors, groups, S, counts, min_obs=C["min_obs"],
        min_self=C["sign_min_rho"], max_rank=C["specificity_max_rank"])

    # ---- call composition and what the restriction removed
    called = call[usable] if restricted else call
    uniq, cnt = np.unique([ls.broad_class.get(c, c) for c in called], return_counts=True)
    out.call_share = dict(zip(uniq.tolist(), (cnt / max(cnt.sum(), 1)).tolist()))

    # Does the ARGMAX carry what the probability vector knows?  Reported per anchor as
    # (share called this compartment in the marker's top decile, same in its bottom).
    call_comp = np.array([ls.broad_class.get(c, c) for c in call])
    for marker, comp in anchors:
        if comp not in groups:
            continue
        x = np.asarray(inp.measured.get(marker, np.full(len(P), np.nan)), dtype=float)
        hi, lo = _decile_contrast(x)
        hi, lo = hi & usable, lo & usable
        out.coherence[marker] = (
            float((call_comp[hi] == comp).mean()) if hi.sum() >= C["min_obs"] else np.nan,
            float((call_comp[lo] == comp).mean()) if lo.sum() >= C["min_obs"] else np.nan)
    out.impossible_mass_mean = float(np.mean(res.impossible_mass))
    out.impossible_outweighed = float(np.mean(res.ambiguous_reason == "mass_on_impossible_classes"))

    # ---- does max-softmax predict agreement?  (licenses stratification; never blocks)
    conf = np.where(usable, np.nanmax(P, axis=1), np.nan)
    fin = np.isfinite(conf)
    if fin.sum() >= 3 * C["min_obs"]:
        b1, b2 = np.quantile(conf[fin], [1 / 3, 2 / 3])
        for i, (marker, comp) in enumerate(anchors):
            if comp not in groups:
                continue
            j = groups.index(comp)
            x = np.asarray(inp.measured.get(marker, np.full(len(P), np.nan)), dtype=float)
            hi, lo = _decile_contrast(x)
            pair = []
            for m in (fin & (conf <= b1), fin & (conf >= b2)):
                sel = (hi | lo) & m
                pair.append(_auroc(G[sel, j], hi[sel]) if sel.sum() >= C["min_obs"] else np.nan)
            out.confidence_auroc[marker] = tuple(pair)

    # ---- permutation null: the same statistics with the model's rows shuffled
    rng = np.random.default_rng(rng_seed)
    nr, na = [], []
    for _ in range(int(C["null_draws"])):
        perm = rng.permutation(len(P))
        Gp = G[perm]
        marker, comp = anchors[rng.integers(len(anchors))]
        if comp not in groups:
            continue
        j = groups.index(comp)
        x = np.asarray(inp.measured.get(marker, np.full(len(P), np.nan)), dtype=float)
        ok = np.isfinite(x) & np.isfinite(Gp[:, j])
        if ok.sum() >= C["min_obs"]:
            nr.append(abs(spearmanr(x[ok], Gp[ok, j]).statistic))
            hi, lo = _decile_contrast(x)
            sel = hi | lo
            na.append(abs(_auroc(Gp[sel, j], hi[sel]) - 0.5))
    q = float(C["null_quantile"])
    out.null = {"rho_abs": float(np.quantile(nr, q)) if nr else np.nan,
                "auroc_dev": float(np.quantile(na, q)) if na else np.nan,
                "draws": len(nr), "quantile": q}

    out.criteria = _judge(out, anchors, groups, C)
    return out


def _judge(out: GateResult, anchors, groups, C) -> dict:
    n = len(anchors)
    crit = {}

    bad_sign = [m for m, _ in anchors
                if not (np.isfinite(out.rho.get(m, {}).get(_comp(anchors, m), np.nan))
                        and out.rho[m][_comp(anchors, m)] > C["sign_min_rho"])]
    crit["sign"] = (not bad_sign, f"{n - len(bad_sign)}/{n} positive"
                                  + (f"; negative or unjudged: {', '.join(bad_sign)}" if bad_sign else ""))

    rank1 = [m for m, _ in anchors if out.rank.get(m, 99) == 1]
    worse = [m for m, _ in anchors if out.rank.get(m, 99) > C["specificity_max_rank"]]
    crit["specificity"] = (
        len(rank1) >= C["specificity_rank1_frac"] * n and not worse,
        f"{len(rank1)}/{n} rank 1 (need {int(np.ceil(C['specificity_rank1_frac'] * n))}); "
        f"{len(worse)} worse than rank {C['specificity_max_rank']}"
        + (f": {', '.join(worse)}" if worse else ""))

    fails = []
    for m, comp in anchors:
        bar = C["auroc_immune"] if comp == "Immune" else C["auroc_other"]
        a = out.auroc.get(m, np.nan)
        if not (np.isfinite(a) and a >= bar):
            fails.append(f"{m}={a:.3f}" if np.isfinite(a) else f"{m}=unjudged")
    crit["detectability"] = (not fails, f"{n - len(fails)}/{n} clear their bar"
                                        + (f"; below: {', '.join(fails)}" if fails else ""))

    # (4) Argmax coherence.  An earlier draft of this criterion capped and floored each
    # compartment's SHARE of calls.  That was invalid and is recorded here rather than
    # quietly dropped: ~85 % of exocrine pancreas is acinar, so a *correct* model would
    # have failed a 45 % cap, and demanding >=0.5 % Neural calls is an expected-frequency
    # assumption -- the thing CLAUDE.md's first non-negotiable forbids.  Replaced before
    # any verdict was recorded, and the replacement is strictly internal: if a
    # compartment's probability separates a marker's top decile from its bottom, the
    # argmax must separate them the same way.  No external frequency is consulted.
    pairs = [(m, out.coherence[m]) for m, _ in anchors
             if m in out.coherence and all(np.isfinite(v) for v in out.coherence[m])]
    incoherent = [f"{m}({t:.1%} vs {b:.1%})" for m, (t, b) in pairs if not t > b]
    crit["argmax_coherence"] = (
        bool(pairs) and len(pairs) - len(incoherent) >= C["coherence_frac"] * len(pairs),
        f"{len(pairs) - len(incoherent)}/{len(pairs)} anchors are called their own compartment "
        f"more often in the marker's top decile than its bottom "
        f"(need {int(np.ceil(C['coherence_frac'] * max(len(pairs), 1)))})"
        + (f"; not: {', '.join(incoherent)}" if incoherent else ""))

    crit["impossible_mass"] = (
        out.impossible_mass_mean <= C["impossible_mass_mean"]
        and out.impossible_outweighed <= C["impossible_mass_outweighed"],
        f"mean {out.impossible_mass_mean:.3f} (bar {C['impossible_mass_mean']}); "
        f"outweighed {out.impossible_outweighed:.1%} (bar {C['impossible_mass_outweighed']:.0%})")

    nr, na = out.null.get("rho_abs", np.nan), out.null.get("auroc_dev", np.nan)
    if not (np.isfinite(nr) and np.isfinite(na)):
        crit["null_separation"] = (False, "permutation null not measured -- cannot judge")
    else:
        weak = [m for m, comp in anchors
                if np.isfinite(out.auroc.get(m, np.nan)) and abs(out.auroc[m] - 0.5) <= na]
        crit["null_separation"] = (not weak,
                                   f"null |rho|<={nr:.3f}, |AUROC-0.5|<={na:.3f} at q{out.null['quantile']}; "
                                   + (f"inside the null: {', '.join(weak)}" if weak else "all anchors outside"))

    pairs = [v for v in out.confidence_auroc.values()
             if np.isfinite(v[0]) and np.isfinite(v[1])]
    better = sum(1 for lo, hi in pairs if hi >= lo)
    crit["confidence_informativeness"] = (
        bool(pairs) and better >= C["confidence_frac"] * len(pairs),
        f"high >= low on {better}/{len(pairs)} judged anchors "
        f"(need {int(np.ceil(C['confidence_frac'] * max(len(pairs), 1)))}); "
        "does NOT block release -- it licenses confidence stratification")
    return crit


def _comp(anchors, marker):
    return next(c for m, c in anchors if m == marker)
