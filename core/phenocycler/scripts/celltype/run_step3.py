"""Step 3: score every cell against the rule table, by cosine, on the aligned values.

    python scripts/celltype/run_step3.py                  # whole cohort
    python scripts/celltype/run_step3.py --donors 6450    # one donor
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

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel                                     # noqa: E402
import rules as R                                # noqa: E402
import score as SC                               # noqa: E402

ALIGNED = panel.OUT / "uniform"
OUT = panel.OUT / "calls"
ENCODING = "signed"          # measured: resolves 37 of 37 references; zero leaves 34

#: Allow the DAPI ratio where -- and only where -- it is a like-for-like comparison.
#:
#: The eight raw-nucleus markers reach ``expression`` as uncompensated nucleus means while the
#: other 52 are REDSEA whole-cell means, and dividing the group by DAPI was meant to put it on
#: a common internal reference.  Measured, that helps in exactly one shape and harms in every
#: other: Ki67 & PCNA survives at 375.7x -> 147.0x because both sides carry the same
#: denominator, while cross-compartment partner tails fall 53.6x -> 0.7x and the TOX &
#: E-cadherin negative control rises 0.0x -> 10.5x -- an association the ratio manufactures.
#:
#: So this is no longer a global on/off.  ``score.nuclear_only_references`` decides per
#: reference, and a reference qualifies only if every marker it constrains is nucleus-measured.
#: ``--no-nuclear-ratio`` still switches it off outright.
NUCLEAR_RATIO = True

#: How far above a marker's own measured background a cell must sit for that marker to count
#: as shown, on the [background mode, robust top] scale.  One number for the whole panel,
#: swept and reported -- not a cut per marker, and deliberately not a quantile, which would
#: assume how common each cell type is.
REQUIRE_MIN = 0.10


def _read_donor(donor: str, markers: list, extra: list, ratio: bool):
    """One donor's aligned values, with the nuclear group put on DAPI if asked.

    Availability is intersected with what the ALIGNED file actually holds, not only with what
    the cohort table says the donor carries.  The two differ: Steps 1 and 2 ran on the cohort
    intersection, so CD38 and b_Catenin1 were never aligned even though 13 and 4 donors
    respectively acquired them.  A reference needing one is therefore scored on no donor at
    present -- a gap to close by re-running Steps 1 and 2 with ``include_variant=True``, not
    something to paper over by requesting a column that is not there.
    """
    on_disk = {f.name for f in
               pq.ParquetFile(ALIGNED / f"donor_id={donor}" / "data.parquet").schema_arrow}
    have = [m for m in markers if m in set(panel.donor_markers(donor)) & on_disk]
    need = list(dict.fromkeys(have + extra + ["DAPI"]))
    need = [c for c in need if c not in ("donor_id",)]
    df = pq.read_table(ALIGNED / f"donor_id={donor}" / "data.parquet",
                       columns=need).to_pandas()
    changed = []
    if ratio:
        df, changed = SC.nuclear_ratio(df, have, "DAPI", panel.NUCLEUS_MEASURED)
    return df, have, changed


def build_scale(donors: list, markers: list, n_per_donor: int = 200_000,
                seed: int = 0, ratio: bool = NUCLEAR_RATIO) -> pd.DataFrame:
    """Cohort-wide floor and top per marker, on the ALIGNED values.

    Cohort-wide on purpose.  A per-donor scale would re-impose exactly the distribution
    matching Step 2 avoided, and the shift has already put the donors on one scale -- so one
    floor and one top per marker is the whole point of having done it.
    """
    rng = np.random.default_rng(seed)
    pool: dict = {m: [] for m in markers}
    for d in donors:
        t, have, _ = _read_donor(d, markers, ["qc_analysis_eligible"], ratio)
        t = t[t["qc_analysis_eligible"]]
        if len(t) > n_per_donor:
            t = t.iloc[rng.choice(len(t), n_per_donor, replace=False)]
        for m in have:
            pool[m].append(t[m].to_numpy("float64"))
    floors = {}
    for m, chunks in pool.items():
        if not chunks:
            continue
        v = np.concatenate(chunks)
        v = v[np.isfinite(v) & (v > 0)]
        if len(v) < 1000:
            continue
        h, e = np.histogram(np.log1p(v), bins=512)
        floors[m] = float(np.expm1(e[1:][int(np.argmax(h[1:])) + 1]))   # the aligned mode
    return SC.marker_scale({m: np.concatenate(c) for m, c in pool.items() if c}, floors)


def _nuclear_scale(scale: pd.DataFrame, donors: list, markers: list) -> pd.DataFrame:
    """``scale`` with the nucleus-measured rows recomputed on the DAPI ratio.

    Dividing changes those markers' distribution, so their floor and top must be derived from
    the ratio rather than inherited from the untouched values.  Built lazily and only over
    that group, because on the current rule table nothing asks for it.
    """
    nuc = [m for m in markers if m in panel.NUCLEUS_MEASURED and m != "DAPI"]
    if not nuc:
        return scale
    out = scale.copy()
    ratio_scale = build_scale(donors, nuc, ratio=True)
    for m in ratio_scale.index:
        if m in out.index:
            out.loc[m] = ratio_scale.loc[m]
    return out


def score_donor(donor: str, Rm: pd.DataFrame, scale: pd.DataFrame,
                ceiling: bool = True, ratio: bool = NUCLEAR_RATIO,
                reqs: dict | None = None, require_min: float = REQUIRE_MIN) -> pd.DataFrame:
    """Cosine-score one donor, over the markers that donor actually carries.

    A reference that needs a marker this donor lacks is not scored here -- it is scored on
    the donors that have it.  Zeroing the missing row instead would quietly turn a
    13-donor comparison into a 20-donor one with weaker evidence on 7 of them.
    """
    missing = [m for m in Rm.index if m not in set(panel.donor_markers(donor))]
    usable = [p for p in Rm.columns if not (Rm.loc[missing, p] != 0).any()]

    keep = ["object_id", "qc_analysis_eligible", "cell_region", "islet_num"]
    df, have, _ = _read_donor(donor, list(Rm.index), keep, ratio=False)
    elig = df["qc_analysis_eligible"].to_numpy()
    X = SC.apply_scale(df[have].to_numpy("float32"), have, scale)

    # The ratio is a per-reference decision, not a per-cell one -- see NUCLEAR_RATIO.
    nuc_ref = (SC.nuclear_only_references(Rm[usable], panel.NUCLEUS_MEASURED)
               if ratio else np.zeros(len(usable), dtype=bool))
    X_nuc = None
    if nuc_ref.any():
        dfr, _ = SC.nuclear_ratio(df.copy(), have, "DAPI", panel.NUCLEUS_MEASURED)
        X_nuc = SC.apply_scale(dfr[have].to_numpy("float32"), have,
                               _nuclear_scale(scale, [donor], have))

    S, allowed = SC.score_and_gate(X, X_nuc, Rm.loc[have, usable].to_numpy("float32"),
                                   have, reqs, usable, require_min, nuc_ref, ceiling)
    S[~elig] = np.nan
    if allowed is not None:
        allowed &= elig[:, None]
    out = SC.call_from_scores(S, usable, allowed=allowed)
    out.insert(0, "object_id", df["object_id"].to_numpy())
    out.insert(0, "donor_id", donor)
    out["cell_region"] = df["cell_region"].to_numpy()
    out["islet_num"] = df["islet_num"].to_numpy()
    out["n_references_scored"] = len(usable)
    out["n_markers_scored"] = len(have)
    return out


def islet_concordance(calls: pd.DataFrame) -> pd.DataFrame:
    """Enrichment of each call inside QuPath islet cores, over that donor's own base rate.

    `cell_region` comes from polygon containment and cannot be predicted from intensity, so
    it is the one check here that owes nothing to the markers the call was made from.
    """
    rows = []
    for d, g in calls.groupby("donor_id"):
        g = g[~g["call"].isin(["Unscored", "Ambiguous"])]
        base = (g["cell_region"] == "core").mean()
        for c, gg in g.groupby("call"):
            if len(gg) < 100:
                continue
            frac = (gg["cell_region"] == "core").mean()
            rows.append({"donor_id": d, "call": c, "n": len(gg), "core_frac": frac,
                         "base_rate": base,
                         "enrichment": frac / base if base > 0 else np.nan})
    return pd.DataFrame(rows)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--encoding", default=ENCODING, choices=("zero", "signed"))
    ap.add_argument("--no-nuclear-ratio", action="store_true",
                    help="score the raw-nucleus markers as-is instead of dividing by DAPI")
    ap.add_argument("--no-gate", action="store_true",
                    help="rank by cosine alone, without enforcing what each rule demands")
    ap.add_argument("--require-min", type=float, default=REQUIRE_MIN)
    ap.add_argument("--no-ceiling", action="store_true",
                    help="do not divide by each reference's attainable maximum")
    args = ap.parse_args()

    donors = args.donors or panel.cohort()["donor_id"].tolist()
    flat = R.flatten(R.load_rules())
    Rm, dropped = R.reference_matrix(flat, panel.markers(include_variant=True),
                                     args.encoding, strict=False)
    Rm = Rm.loc[(Rm != 0).any(axis=1)]                     # markers no reference uses
    print(f"Step 3 — {R.RULES_CSV.name}: {len(Rm.columns)} references over {len(Rm)} markers, "
          f"encoding {args.encoding}")
    if dropped:
        print(f"  unencodable: {dropped}")
    ev = R.evaluable_donors(flat)
    for _, r in ev[ev["n_donors"] < len(donors)].iterrows():
        print(f"  {r['phenotype']}: scored on {r['n_donors']} donors (needs "
              f"{r['variant_markers']})")

    t0 = time.time()
    print("  building the cohort-wide marker scale ...", flush=True)
    ratio = not args.no_nuclear_ratio
    print(f"  nuclear markers on DAPI: {'yes' if ratio else 'no'} "
          f"({', '.join(m for m in panel.NUCLEUS_MEASURED if m != 'DAPI')})")
    scale = build_scale(donors, list(Rm.index), ratio=ratio)
    reqs = None if args.no_gate else R.requirements(R.load_rules())
    if reqs is not None:
        n_all = int(np.median([len(v[0]) for v in reqs.values()]))
        print(f"  requirement gate ON: median {n_all} markers must be shown per reference, "
              f"at >= {args.require_min} above their own background")
    print(f"    {time.time() - t0:.0f}s", flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    frames = []
    for d in donors:
        t = time.time()
        c = score_donor(d, Rm, scale, ceiling=not args.no_ceiling, ratio=ratio,
                        reqs=reqs, require_min=args.require_min)
        p = OUT / f"donor_id={d}"
        p.mkdir(parents=True, exist_ok=True)
        c.to_parquet(p / "data.parquet", index=False)
        frames.append(c)
        n = (c["call"] != "Unscored").sum()
        print(f"  {d}: {len(c):,} cells, {n:,} scored, "
              f"{c['n_references_scored'].iloc[0]} references   {time.time() - t:.0f}s",
              flush=True)
    calls = pd.concat(frames, ignore_index=True)

    print("\ncall composition, cohort (scored cells only):")
    v = calls.loc[~calls["call"].isin(["Unscored", "Ambiguous"]), "call"].value_counts(
        normalize=True)
    for c, f in v.items():
        print(f"    {c:<38}{f:7.2%}")
    print(f"    {'Ambiguous (no rule satisfied)':<38}"
          f"{(calls['call'] == 'Ambiguous').mean():7.2%}")
    print(f"    {'Unscored (no signal at all)':<38}"
          f"{(calls['call'] == 'Unscored').mean():7.2%}")

    conc = islet_concordance(calls)
    conc.to_parquet(panel.OUT / "islet_concordance.parquet", index=False)
    print("\nislet-core enrichment over each donor's own base rate, pooled (>=100 cells):")
    agg = (conc.groupby("call")
                .agg(n=("n", "sum"), donors=("donor_id", "nunique"),
                     enrichment=("enrichment", "median"))
                .sort_values("enrichment", ascending=False))
    print(agg.to_string(float_format=lambda x: f"{x:8.2f}"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
