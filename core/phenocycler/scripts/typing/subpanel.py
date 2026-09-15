#!/usr/bin/env python3
"""The subpanel classification: context from geometry, then one question per context.

WHY THIS FILE EXISTS. `subpanel_calls.parquet` -- 25,150,537 rows, the typing the whole app and
every downstream lane reads -- had **no producer**. It was made in a session and the script was
never saved: an exhaustive search of the working tree, every branch, the full git history and the
other checkouts finds only `subpanel_io.py` (its inputs), `build_cells.py` (a reader) and prose.
So the typing could not be re-run, which meant no correction to it could ever be applied.

This is that producer, written from `RECONSTRUCTION.md` -- which was itself recovered by
reproducing the calls from their inputs and checked against them, not inferred from prose.

THE METHOD, unchanged from what was recovered:

    context comes from geometry; each context is asked only what its own markers can answer;
    every level is the 99th percentile of a REFERENCE POPULATION, never a constant.

  reference population   cells in no classified object and outside every islet -- the ordinary
                         tissue a structure sits in. It excludes the objects and the islets so it
                         cannot contain the thing being detected.
  levels (pooled, p99)   hormone 0.1010 · Nerve 0.2272 · Lymphatic 0.1776 · Blood vessel 0.3073
  scale                  READ from the celltype lane, never recomputed -- `build_scale` draws a
                         seeded 200,000-cell sample per donor, so recomputing it would drift.

TWO THINGS THIS CHANGES FROM THE TABLE ON DISK, both deliberate, both behind a flag:

  --membership geometry  Islet membership from the islet OUTLINES, not QuPath's `Parent`. The
                         Parent labels cells `peri` that sit up to 74.6 um inside a core outline;
                         because a false `in_islet` means the cell is never asked the islet
                         question, 11,556 islet-core cells fell through to Exocrine, 5,250 of them
                         hormone-positive. `--membership parent` reproduces the old behaviour and
                         is the DEFAULT, because stage 1 is a reconstruction gate.
  --emit-detail          Carries `immune_definition` and `resolved_by` as columns. The table on
                         disk collapsed both into the parent label and threw them away.

THE GATE COMES FIRST. `--gate` rebuilds with `--membership parent` and compares row-for-row on
(donor_id, object_id, cell_class) against the table on disk. Until that matches, a difference
cannot be attributed: it is either a reconstruction bug or a correction working, and there is no
way to tell which. Two runs were lost to exactly that confusion before.

WHAT IS NOT RECOVERED, and is reported rather than hidden: which cells the supporting tier is
eligible to claim. Applying it to every unclaimed cell reproduces 99.998 % on 6442 and 115 but
99.59 % on 6436, so an eligibility condition remains unrecovered. This runs the permissive form
and prints the per-donor residual every time, because a reconstruction that quietly papers over
its own known gap is worth less than one that names it.

    python core/phenocycler/scripts/typing/subpanel.py --gate --donors 6442 6436 115
    python core/phenocycler/scripts/typing/subpanel.py --membership geometry --emit-detail
"""
from __future__ import annotations

import argparse
import json
import pathlib
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "celltype"))
import context as CX
import islet_membership as MEMBERSHIP
import lineages as LIN
import rules as RULES
import subpanel_io as S

OUT = S.OUT / "subpanel"
WORK = S.OUT / "scratch/subpanel_context"

#: Strongest first. Recovered from the calls rather than assumed -- see RECONSTRUCTION.md.
STRUCTURES = ("Blood vessel", "Nerve", "Lymphatic")

#: The marker each pixel classifier was actually BUILT ON -- what drew the object.
#:
#: Deliberately not the same as `S.SUBPANEL`, which is what ANSWERS the question once a cell is
#: in the object. Nerve's subpanel reads CD56 as well, but the object came from
#: `B3Tub_Ero0.5_Thr600` -- B3TUBB alone. Vessel objects came from `Cavdil1200` and `CD34dil1430`,
#: so Caveolin and CD34 and NOT CD31, which never drew anything. Lymphatic from
#: `PDPNmorphClos0-5_Thr300`.
CLASSIFIER_MARKER: dict[str, tuple[str, ...]] = {
    "Nerve": ("B3TUBB",),
    "Lymphatic": ("Podoplanin",),
    "Blood vessel": ("Caveolin", "CD34"),
}

#: Markers that can carry an identity, from `lineages.LINEAGE_MARKERS` -- "a statement of what
#: each antibody means, made before any result is seen". State markers (proliferation,
#: exhaustion, activation) are absent by that module's own design: they say what a cell is doing,
#: not what it is, so they must not be able to win a cell's top-marker argument.
IDENTITY_MARKERS: tuple[str, ...] = tuple(sorted(
    {m for markers in LIN.LINEAGE_MARKERS.values() for m in markers}))

#: Scaled value an immune rule's markers must clear. Recovered, not chosen.
IMMUNE_REQUIRE_MIN = 0.10

#: The supporting tier's level over `bottom_tier_immune`, cohort-wide. Recovers independently on
#: 6442 (0.8715) and 6436 (0.8716).
SUPPORTING_LEVEL = 0.8715

#: A cell inside an islet that falls through to Exocrine is not called Exocrine -- the maintainer
#: ruled this set out of analysis. It keeps its own reason so the count stays auditable.
ISLET_EXOCRINE_REASON = "exocrine_inside_islet_border"


def _log(): return json.loads((OUT / "subpanel_log.json").read_text())


def scaled_markers(donor: str, wanted: set[str], scale: pd.DataFrame,
                   object_ids: pd.Series) -> pd.DataFrame:
    """Scaled values for the markers this donor actually carries, joined to `object_ids`.

    By NAME, never by position: three panels exist and a marker absent from one donor is absent,
    not shifted. An absent marker yields a column of NaN, and `NaN >= level` is False, so a rule
    needing it simply cannot fire -- which is the correct behaviour and not zero-filling.
    """
    path = S.ALIGNED / f"donor_id={donor}/data.parquet"
    have = set(pq.ParquetFile(path).schema_arrow.names)
    present = sorted(wanted & have)
    frame = pq.read_table(path, columns=["object_id"] + present).to_pandas()
    frame = object_ids.to_frame("object_id").merge(frame, on="object_id", how="left")
    values = pd.DataFrame(S.apply_scale(frame[present].to_numpy(), present, scale),
                          columns=present, index=frame.index)
    for marker in sorted(wanted - have):
        values[marker] = np.nan
    return values


#: Columns the geometry half of the classification needs, whichever source supplies them.
GEOMETRY_COLUMNS = ["object_id", "X_centroid", "Y_centroid", "islet_num", "cell_region",
                    "qc_analysis_eligible"]


def donor_context(donor: str, membership: str, pixel: dict, stem: dict) -> pd.DataFrame:
    """object_id, context, cell_region for every analysis-eligible cell of one donor."""
    path, _ = S.cell_source(donor)
    run = pq.read_table(path, columns=GEOMETRY_COLUMNS).to_pandas()
    run = run[run["qc_analysis_eligible"]].reset_index(drop=True)

    if membership == "parent":
        parent = run["islet_num"].fillna("").astype(str).str.strip()
        in_islet = ((parent != "") & run["cell_region"].eq("core")).to_numpy()
    elif membership == "geometry":
        geo = MEMBERSHIP.membership(donor, run["X_centroid"].to_numpy(),
                                    run["Y_centroid"].to_numpy(), pixel[donor], stem[donor])
        in_islet = geo["in_islet"].to_numpy()
    else:
        raise ValueError(f"unknown membership {membership!r}")

    # Ask the object question of EVERY cell, islet included -- correction C. `context` keeps the
    # precedence answer; `object_context` keeps what the object actually said, which the original
    # table discarded and which is why an earlier analysis of this same lane went wrong.
    object_context = CX.assign(run["X_centroid"].to_numpy(), run["Y_centroid"].to_numpy(),
                               np.zeros(len(run), dtype=bool), stem[donor], pixel[donor])
    context = np.where(in_islet, "islet", object_context)
    frame = pd.DataFrame({"object_id": run["object_id"], "context": context,
                          "object_context": object_context, "in_islet": in_islet,
                          "cell_region": run["cell_region"]})
    # Per-structure membership, NOT the collapsed label: a cell inside both a vessel and a
    # lymphatic must be able to answer for both. See `context.memberships`.
    inside = CX.memberships(run["X_centroid"].to_numpy(), run["Y_centroid"].to_numpy(),
                            stem[donor], pixel[donor])
    for name, mask in inside.items():
        frame[f"in_{name}"] = mask
    return frame


LEVEL_CACHE = WORK / "reference_levels.json"


def fit_levels(donors: list[str], scale: pd.DataFrame, log: dict,
               membership: str = "geometry") -> dict:
    """Every structure's level, from the reference population pooled over the given donors.

    POOLING IS NOT THE DEFAULT. It is kept only so `--gate` can reproduce the table on disk, which
    was built this way. Measured 2026-09-15, the same reference population fitted per donor gives:

        hormone       0.0253 - 0.4936   19.5x spread   (pooled 0.1010)
        Blood vessel  0.1031 - 0.5387    5.2x          (pooled 0.3073)
        Nerve         0.0913 - 0.6631    7.3x          (pooled 0.2272)
        Lymphatic     0.0592 - 0.4175    7.0x          (pooled 0.1776)

    and the hormone level separates almost perfectly by acquisition batch -- batch 1 median
    0.0418, batch 2 median 0.3167, with only 6436 crossing. So the pooled level is a blend of two
    batch populations: 2.4x too strict for batch 1 and 3x too permissive for batch 2. Acquisition
    batch is aliased with disease status at the coarse level, so a level built that way carries
    the confound directly into the calls.

    Fitting per donor also makes the level mean what it says -- "above the ordinary tissue THIS
    structure sits in" -- since the reference population is that donor's own tissue.
    """
    pools: dict[str, list] = {name: [] for name in ("hormone", *STRUCTURES)}
    for donor in donors:
        ctx = pd.read_parquet(WORK / f"{donor}__{membership}.parquet")
        tissue = ctx["context"].to_numpy() == "tissue"
        if not tissue.any():
            continue
        wanted = set(S.HORMONES) | {m for name in STRUCTURES for m in log["subpanels"][name]}
        values = scaled_markers(donor, wanted, scale, ctx["object_id"])
        for name, markers in [("hormone", list(S.HORMONES))] + \
                             [(n, log["subpanels"][n]) for n in STRUCTURES]:
            present = [m for m in markers if values[m].notna().any()]
            if not present:
                continue
            pools[name].append(np.nanmax(values.loc[tissue, present].to_numpy(), axis=1))
    levels = {}
    for name, chunks in pools.items():
        if not chunks:
            raise ValueError(f"{name}: reference population is empty -- refusing to classify")
        levels[name] = S.reference_level(np.concatenate(chunks))
    return levels


def immune_rules() -> dict:
    """phenotype -> (all_of, [any_of groups], none_of), immune phenotypes only."""
    table = RULES.requirements(RULES.load_rules(), with_negatives=True)
    return {k: v for k, v in table.items() if k.startswith(("Imm-", "T-"))}


def fires(values: pd.DataFrame, rule: tuple) -> np.ndarray:
    """Does each cell satisfy one immune definition, with its declared negatives enforced.

    A marker the donor does not carry is NaN, and every comparison against NaN is False. That
    makes an all_of containing an absent marker unsatisfiable (correct -- the evidence is not
    there) and an any_of fall back to whichever members ARE present (also correct).
    """
    all_of, any_of, none_of = rule
    ok = np.ones(len(values), dtype=bool)
    for marker in all_of:
        ok &= (values[marker].to_numpy() >= IMMUNE_REQUIRE_MIN)
    for group in any_of:
        members = [m for m in group if m in values.columns]
        if not members:
            continue
        ok &= (values[members].to_numpy() >= IMMUNE_REQUIRE_MIN).any(axis=1)
    for marker in none_of:
        # A negative must be ABSENT. NaN is "not measured", which cannot demonstrate presence,
        # so it does not veto -- the alternative silently kills every rule on a partial panel.
        v = values[marker].to_numpy()
        ok &= ~(v >= IMMUNE_REQUIRE_MIN)
    return ok


def classify(donor: str, ctx: pd.DataFrame, levels: dict, scale: pd.DataFrame,
             log: dict, detail: bool, require_top_marker: bool = True) -> pd.DataFrame:
    rules = immune_rules()
    wanted = set(S.HORMONES) | {"CD99"}
    for markers in log["subpanels"].values():
        wanted |= set(markers)
    for all_of, any_of, none_of in rules.values():
        wanted |= set(all_of) | set(none_of) | {m for g in any_of for m in g}
    wanted |= set(log["bottom_tier_immune"])
    values = scaled_markers(donor, wanted, scale, ctx["object_id"])

    n = len(ctx)
    cell_class = np.full(n, None, dtype=object)
    definition = np.full(n, "", dtype=object)
    resolved_by = np.full(n, "", dtype=object)
    reason = np.full(n, "", dtype=object)
    obj = ctx["object_context"].to_numpy()
    passed_through: dict[str, np.ndarray] = {}

    # --- 1. structures, weakest first so a stronger class overwrites -------------------------
    # Asked of EVERY cell, islet included: the object evidence is real inside an islet too.
    #
    # THREE conditions, not one. A centroid inside the object is necessary and nowhere near
    # sufficient -- an object is a region of tissue and most of what sits in it is whatever the
    # region runs through. So the cell must also clear that structure's level AND its own
    # top-expressed identity marker must be one of the markers that DREW the object. A cell whose
    # brightest identity marker is Ker8_18 is not a blood vessel because it happens to lie inside
    # one; it is an exocrine cell the vessel passes through, and it keeps its own identity and
    # carries the crossing as a feature instead.
    identity = [m for m in IDENTITY_MARKERS if m in values.columns]
    identity_values = values[identity].to_numpy()
    top_marker = np.array(identity, dtype=object)[
        np.nanargmax(np.nan_to_num(identity_values, nan=-1.0), axis=1)]

    for name in reversed(STRUCTURES):
        markers = [m for m in log["subpanels"][name] if m in values.columns]
        if not markers:
            continue
        inside = ctx[f"in_{name}"].to_numpy()
        clears = np.nanmax(values[markers].to_numpy(), axis=1) >= levels[name]
        owns = (np.isin(top_marker, CLASSIFIER_MARKER[name]) if require_top_marker
                else np.ones(n, dtype=bool))
        hit = inside & clears & owns
        cell_class[hit] = name
        resolved_by[hit] = "structure"
        # Everything else the object contains is passed through, and says so per cell.
        passed_through[name] = inside & ~hit

    # --- 2. inside an islet: the hormone question ---------------------------------------------
    in_islet = ctx["in_islet"].to_numpy()
    horm = [h for h in S.HORMONES if h in values.columns]
    top = np.nanmax(values[horm].to_numpy(), axis=1)
    which = np.array(horm, dtype=object)[np.nanargmax(np.nan_to_num(values[horm].to_numpy(),
                                                                   nan=-1.0), axis=1)]
    name_of = {v: k for k, v in log["hormones"].items()}
    endocrine = in_islet & (cell_class == None) & (top >= levels["hormone"])  # noqa: E711
    cell_class[endocrine] = [name_of[w] for w in which[endocrine]]
    resolved_by[endocrine] = "hormone"

    # --- 3. the immune definitions ------------------------------------------------------------
    unclaimed = cell_class == None  # noqa: E711
    for phenotype, rule in sorted(rules.items()):
        hit = unclaimed & fires(values, rule)
        cell_class[hit] = "Immune"
        definition[hit] = phenotype
        resolved_by[hit] = "definition"
        unclaimed &= ~hit

    # --- 4. the supporting tier ---------------------------------------------------------------
    # A cell is eligible when its immune evidence CLEARS the level and its exocrine evidence does
    # NOT -- the same level on both sides, so the tier never overrules a cell the exocrine markers
    # speak for as loudly as the immune ones do.
    #
    # This was recovered, not chosen, and the recovery is worth recording because it is the one
    # part of the method RECONSTRUCTION.md left open. All 35,332 cells the original table claims
    # by this tier on 6436 carry exocrine evidence below 0.8715; none is above it. Applying the
    # tier to every unclaimed cell instead -- the permissive form -- wrongly claims 5,671 strongly
    # exocrine cells there (exocrine max p50 0.98). Restricting to cells with NO exocrine evidence
    # was already measured too narrow: it keeps only 19,142 of the 35,332.
    #
    # It is not exact. 1,223 cells on 6436 remain wrongly claimed, 98.4 % of them driven by PD_1
    # on cells whose exocrine evidence sits just under the level (p50 0.843 against 0.044 for the
    # cells the tier should claim). Closing that would mean picking a threshold because it fits
    # the answer, which is not recovery, so it is left open and reported.
    tier = [m for m in log["bottom_tier_immune"] if m in values.columns]
    exocrine = [m for m in S.SUBPANEL["exocrine"] if m in values.columns]
    if tier:
        clears = np.nanmax(values[tier].to_numpy(), axis=1) >= SUPPORTING_LEVEL
        speaks_exocrine = (np.nanmax(values[exocrine].to_numpy(), axis=1) >= SUPPORTING_LEVEL
                           if exocrine else np.zeros(n, dtype=bool))
        hit = unclaimed & clears & ~speaks_exocrine
        cell_class[hit] = "Immune"
        resolved_by[hit] = "supporting"
        unclaimed &= ~hit

    # --- 5. fall through ----------------------------------------------------------------------
    cell_class[unclaimed] = "Exocrine"
    resolved_by[unclaimed] = "fallthrough"

    # A cell inside an islet that fell through is not Exocrine -- ruled out of analysis, by name.
    ruled_out = in_islet & (cell_class == "Exocrine")
    cell_class[ruled_out] = "Excluded"
    reason[ruled_out] = ISLET_EXOCRINE_REASON

    frame = pd.DataFrame({
        "donor_id": donor, "object_id": ctx["object_id"].to_numpy(),
        "cell_class": cell_class.astype(str), "context": ctx["context"].to_numpy(),
        "cell_region": ctx["cell_region"].to_numpy(), "excluded_reason": reason.astype(str)})
    if detail:
        frame["immune_definition"] = definition.astype(str)
        frame["resolved_by"] = resolved_by.astype(str)
        frame["object_context"] = obj
        frame["top_identity_marker"] = top_marker
        for name, mask in passed_through.items():
            frame[f"passed_through_{name.replace(' ', '_')}"] = mask
    return frame


def gate(rebuilt: pd.DataFrame, donor: str) -> dict:
    """Row-for-row against the table on disk, on (donor_id, object_id, cell_class)."""
    old = pq.read_table(OUT / "subpanel_calls.parquet",
                        columns=["object_id", "cell_class", "context"],
                        filters=[("donor_id", "=", donor)]).to_pandas()
    merged = old.merge(rebuilt[["object_id", "cell_class", "context"]], on="object_id",
                       how="outer", suffixes=("_old", "_new"), indicator=True)
    both = merged["_merge"] == "both"
    agree = (merged.loc[both, "cell_class_old"] == merged.loc[both, "cell_class_new"])
    ctx_agree = (merged.loc[both, "context_old"] == merged.loc[both, "context_new"])
    confusion = (merged.loc[both & ~agree.reindex(merged.index, fill_value=False)]
                 .groupby(["cell_class_old", "cell_class_new"]).size()
                 .sort_values(ascending=False).head(6))
    return {"donor_id": donor, "rows_old": int((merged["_merge"] != "right_only").sum()),
            "rows_new": int((merged["_merge"] != "left_only").sum()),
            "only_old": int((merged["_merge"] == "left_only").sum()),
            "only_new": int((merged["_merge"] == "right_only").sum()),
            "class_match_pct": 100.0 * float(agree.mean()) if len(agree) else float("nan"),
            "context_match_pct": 100.0 * float(ctx_agree.mean()) if len(ctx_agree) else float("nan"),
            "confusion": confusion}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--membership", choices=("parent", "geometry"), default="parent")
    ap.add_argument("--emit-detail", action="store_true")
    ap.add_argument("--no-top-marker", action="store_true",
                    help="drop correction C's requirement, to isolate its effect in a diff")
    ap.add_argument("--levels", choices=("per-donor", "pooled"), default="per-donor",
                    help="per-donor is the method; pooled exists only to reproduce the old table")
    ap.add_argument("--refit", action="store_true", help="recompute context and levels")
    ap.add_argument("--gate", action="store_true",
                    help="compare against the table on disk instead of writing")
    ap.add_argument("--out", default=None, help="output directory (default: subpanel_rebuilt)")
    args = ap.parse_args()

    log = _log()
    scale = S.load_scale()
    pixel, stem = MEMBERSHIP._manifest()
    every = sorted(p.name.split("=")[1] for p in (S.RUN / "expression").glob("donor_id=*"))
    donors = args.donors or every
    WORK.mkdir(parents=True, exist_ok=True)

    print(f"membership: {args.membership}   classifying {len(donors)} of {len(every)} donors")
    t0 = time.time()
    # Per donor, a donor needs only its own context. Pooled levels cannot be fitted without the
    # whole cohort, so that mode -- and only that mode -- pays for every donor.
    scope = every if args.levels == "pooled" else donors
    # Keyed by membership as well as donor: the two modes give different context, and a cache
    # that ignored that would silently serve one mode's answer to the other.
    key = lambda d: WORK / f"{d}__{args.membership}.parquet"
    stale = [d for d in scope if not key(d).exists() or args.refit]
    for donor in stale:
        ctx = donor_context(donor, args.membership, pixel, stem)
        ctx.to_parquet(key(donor), index=False)
        print(f"  context {donor}: {len(ctx):,} eligible cells, "
              f"{int(ctx['in_islet'].sum()):,} in an islet  [{S.cell_source(donor)[1]}]"
              f"   [{time.time()-t0:.0f}s]", flush=True)

    if args.levels == "pooled":
        if LEVEL_CACHE.exists() and not args.refit and not stale:
            levels = json.loads(LEVEL_CACHE.read_text())
            print("   (pooled levels read from cache; --refit to recompute)")
        else:
            levels = fit_levels(every, scale, log, args.membership)
            LEVEL_CACHE.write_text(json.dumps(levels, indent=1))
        print("\nPOOLED reference levels -- reproduces the old table; NOT the method:")
        for name, value in levels.items():
            recorded = log["bars"].get(name) or (0.1010 if name == "hormone" else None)
            print(f"   {name:14s} {value:.4f}" +
                  (f"   recorded {recorded:.4f}" if recorded else ""))
        per_donor = {d: levels for d in donors}
    else:
        per_donor = {d: fit_levels([d], scale, log, args.membership) for d in donors}
        print("\nPER-DONOR reference levels (p99 of that donor's OWN ordinary tissue):")
        frame = pd.DataFrame(per_donor).T
        frame.index.name = "donor"
        print(frame.to_string(float_format=lambda v: f"{v:.4f}"))

    target = pathlib.Path(args.out) if args.out else S.OUT / "subpanel_rebuilt"
    if not args.gate:
        target.mkdir(parents=True, exist_ok=True)

    if args.gate and args.levels != "pooled":
        print("\nNOTE: --gate against per-donor levels is NOT a reconstruction gate. The table on\n"
              "disk was built with pooled levels, so what follows is the CORRECTION'S EFFECT, not\n"
              "a reproduction failure. Run `--gate --levels pooled` for the gate itself.")

    gates = []
    for donor in donors:
        ctx = pd.read_parquet(key(donor))
        frame = classify(donor, ctx, per_donor[donor], scale, log, args.emit_detail,
                         not args.no_top_marker)
        if args.gate:
            g = gate(frame, donor)
            gates.append(g)
            print(f"\n{donor}: class {g['class_match_pct']:.4f}%  context "
                  f"{g['context_match_pct']:.4f}%   rows old {g['rows_old']:,} / "
                  f"new {g['rows_new']:,}  (only-old {g['only_old']:,}, only-new {g['only_new']:,})")
            if len(g["confusion"]):
                for (a, b), n in g["confusion"].items():
                    print(f"      {a} -> {b}: {n:,}")
        else:
            part = target / f"donor_id={donor}"
            part.mkdir(parents=True, exist_ok=True)
            frame.drop(columns=["donor_id"]).to_parquet(part / "data.parquet", index=False)
            print(f"  {donor}: {len(frame):,} cells  " +
                  "  ".join(f"{k} {v:,}" for k, v in
                            frame['cell_class'].value_counts().head(5).items()), flush=True)

    if args.gate and gates:
        table = pd.DataFrame([{k: v for k, v in g.items() if k != "confusion"} for g in gates])
        pd.set_option("display.width", 200)
        print("\n" + table.to_string(index=False))
        worst = table["class_match_pct"].min()
        print(f"\nGATE: worst donor {worst:.4f}% on cell_class. "
              f"{'PASS' if worst >= 99.99 else 'NOT YET EXACT -- do not apply corrections'}")
    elif not args.gate:
        print(f"\n-> {target}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
