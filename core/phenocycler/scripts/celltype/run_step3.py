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

#: Divide the eight raw-nucleus markers by DAPI before scoring, so that group is compared on
#: a common internal reference rather than against REDSEA whole-cell means.
NUCLEAR_RATIO = True


def _read_donor(donor: str, markers: list, extra: list, ratio: bool):
    """One donor's aligned values, with the nuclear group put on DAPI if asked."""
    have = [m for m in markers if m in set(panel.donor_markers(donor))]
    need = list(dict.fromkeys(have + extra + (["DAPI"] if ratio else [])))
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


def score_donor(donor: str, Rm: pd.DataFrame, scale: pd.DataFrame,
                ceiling: bool = True, ratio: bool = NUCLEAR_RATIO) -> pd.DataFrame:
    """Cosine-score one donor, over the markers that donor actually carries.

    A reference that needs a marker this donor lacks is not scored here -- it is scored on
    the donors that have it.  Zeroing the missing row instead would quietly turn a
    13-donor comparison into a 20-donor one with weaker evidence on 7 of them.
    """
    missing = [m for m in Rm.index if m not in set(panel.donor_markers(donor))]
    usable = [p for p in Rm.columns if not (Rm.loc[missing, p] != 0).any()]

    keep = ["object_id", "qc_analysis_eligible", "cell_region", "islet_num"]
    df, have, _ = _read_donor(donor, list(Rm.index), keep, ratio)
    elig = df["qc_analysis_eligible"].to_numpy()
    X = SC.apply_scale(df[have].to_numpy("float32"), have, scale)
    S = SC.cosine_scores(X, Rm.loc[have, usable].to_numpy("float32"),
                         normalise_ceiling=ceiling)
    S[~elig] = np.nan
    out = SC.call_from_scores(S, usable)
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
        g = g[g["call"] != "Unscored"]
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
    print(f"    {time.time() - t0:.0f}s", flush=True)

    OUT.mkdir(parents=True, exist_ok=True)
    frames = []
    for d in donors:
        t = time.time()
        c = score_donor(d, Rm, scale, ceiling=not args.no_ceiling, ratio=ratio)
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
    v = calls.loc[calls["call"] != "Unscored", "call"].value_counts(normalize=True)
    for c, f in v.items():
        print(f"    {c:<38}{f:7.2%}")
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
