"""Who is inside the islets, given that QuPath already told us where they are.

The islets are a **prior, not a finding**.  Every image in this cohort was run through a
QuPath pixel classifier on CD99, INS, GCG and SST, and its output is the ``Islet_N`` polygons
and their 20 um expansion ring -- already joined per cell as ``cell_region``.  Rediscovering
them is not analysis: UTAG recovers them at 25x enrichment, which is a good check that the
data is sound and no new information at all.

That prior is strong and it is measurable.  Median within-donor percentile inside ``Islet_N``
minus outside it, over six donors: **CD99 +0.49 to +0.50 in every donor** -- more consistent
than INS, which runs +0.04 to +0.51 because T1D donors have lost their beta cells -- with
GCG +0.49 to +0.50 and CD56 +0.35 to +0.45.  Immune markers go the other way (CD3e -0.07 to
-0.29, CD11b -0.09 to -0.33), which is the point: inside an islet the non-endocrine cells are
the minority, and they are what this study is about.

So the work is not to call 25 million cells.  It is:

* **inside ``Islet_N``** -- take endocrine as the default and move a cell off it only on
  positive evidence of something else.  Islet cores are 1.5-6.2% of cells per donor.
* **in the 20 um ring** -- no default.  This is the immune environment around islets, the
  thing the study exists to measure, so every reference competes.

Endocrine SUBTYPE is a separate question and is answered here directly from the hormones,
because inside the prior we already know the cell is endocrine and the only question is which.

    python scripts/celltype/run_islet.py                 # whole cohort
    python scripts/celltype/run_islet.py --donors 6450
"""
from __future__ import annotations

import argparse
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
import run_step3 as S3                           # noqa: E402

OUT = panel.OUT / "islet_residents"

#: References that describe an endocrine cell.  Inside the prior these are not scored -- the
#: pixel classifier already made that call, and re-deriving it is what produced Endo-beta at
#: 0.00% while the same cells sat inside polygons drawn from INS.
ENDOCRINE_PREFIX = ("Epi-endocrine", "Endo-alpha", "Endo-beta", "Endo-delta",
                    "Endo-hormone-negative")


def intruder_references(rules: pd.DataFrame, columns) -> list:
    """References that can legitimately displace the endocrine prior inside an islet.

    The question inside the prior is immune / vascular / lymphatic / other -- not what kind
    of epithelium.  So the whole epithelial branch is out of contention there, ancestors and
    siblings alike, because an endocrine cell satisfies those rules by construction: it is
    E-cadherin+ and carries Ker8/18 and Pan-CK. Left in, the generic ``Epithelial`` parent
    took 42% of core cells off the prior and ``Epi-exocrine`` another 26%.

    That is a scope decision, not a claim that no duct ever crosses an islet outline.  A
    Ker8/18-high cell inside a polygon drawn by an INS/GCG/SST/CD99 pixel classifier is a
    question about the outline, and this lane does not answer it -- the ring, where there is
    no prior, scores every reference including the epithelial ones.
    """
    par = dict(zip(rules["phenotype"], rules["parent"]))

    def branch(ph):
        cur, seen = ph, []
        while cur and cur != "all":
            seen.append(cur)
            cur = par.get(cur)
            if cur in seen:
                break
        return seen

    endo_branch = set()
    for ph in rules["phenotype"]:
        if str(ph).startswith(ENDOCRINE_PREFIX):
            endo_branch |= set(branch(ph))
    implied = {c for c in rules["phenotype"] if endo_branch & set(branch(c))}
    return [c for c in columns if c not in implied]

#: Hormones, for subtyping inside the prior.  Read directly: within an islet the question is
#: not "is this endocrine" but "which hormone", and that is what these markers answer.
HORMONE = {"Endo-beta": "INS", "Endo-alpha": "GCG", "Endo-delta": "SST"}


def subtype_endocrine(Xs: np.ndarray, markers: list, min_scaled: float) -> pd.DataFrame:
    """Which hormone, for cells the prior already says are endocrine.

    A cell showing none is ``Endo-hormone-negative`` -- a real category in T1D islets, not a
    failure to call.  A cell showing more than one keeps the strongest and records the
    runner-up, because hormone spillover between touching endocrine cells is expected and
    hiding it would misreport a measurement as a certainty.
    """
    idx = {m: i for i, m in enumerate(markers)}
    names = [k for k, m in HORMONE.items() if m in idx]
    if not names:
        return pd.DataFrame({"call": ["Endo-unresolved"] * len(Xs), "score": np.nan,
                             "runner_up": "", "margin": np.nan})
    V = np.column_stack([Xs[:, idx[HORMONE[k]]] for k in names])
    order = np.argsort(-V, axis=1)
    top, second = order[:, 0], order[:, 1] if V.shape[1] > 1 else order[:, 0]
    tv = V[np.arange(len(V)), top]
    sv = V[np.arange(len(V)), second] if V.shape[1] > 1 else np.zeros(len(V))
    call = np.where(tv >= min_scaled, np.asarray(names, dtype=object)[top],
                    "Endo-hormone-negative")
    return pd.DataFrame({
        "call": call, "score": tv,
        "runner_up": np.where(tv >= min_scaled, np.asarray(names, dtype=object)[second], ""),
        "margin": np.where(tv >= min_scaled, tv - sv, np.nan)})


def run_donor(donor: str, Rm: pd.DataFrame, scale: pd.DataFrame, reqs: dict,
              require_min: float, ratio: bool, intruders: list) -> pd.DataFrame:
    missing = [m for m in Rm.index if m not in set(panel.donor_markers(donor))]
    usable = [p for p in Rm.columns if not (Rm.loc[missing, p] != 0).any()]
    intruder = [p for p in usable if p in intruders]

    keep = ["object_id", "qc_analysis_eligible", "cell_region", "islet_num"]
    df, have, _ = S3._read_donor(donor, list(Rm.index), keep, ratio=False)
    df = df[df["qc_analysis_eligible"] & df["cell_region"].isin(("core", "peri"))]
    if df.empty:
        return pd.DataFrame()
    df = df.reset_index(drop=True)
    Xs = SC.apply_scale(df[have].to_numpy("float32"), have, scale)
    # The DAPI ratio is a per-reference decision -- see run_step3.NUCLEAR_RATIO.
    nuc_all = (SC.nuclear_only_references(Rm, panel.NUCLEUS_MEASURED)
               if ratio else np.zeros(len(Rm.columns), dtype=bool))
    Xs_nuc = None
    if nuc_all.any():
        dfr, _ = SC.nuclear_ratio(df.copy(), have, "DAPI", panel.NUCLEUS_MEASURED)
        Xs_nuc = SC.apply_scale(dfr[have].to_numpy("float32"), have,
                                S3._nuclear_scale(scale, [donor], have))
    nuc_by_ref = dict(zip(Rm.columns, nuc_all))
    core = (df["cell_region"] == "core").to_numpy()

    # The two regions are different questions and get different reference sets.  Inside
    # Islet_N the prior settles the epithelial one, so only intruders compete.  The ring has no
    # prior -- it is ordinary tissue, mostly acinar -- so every reference does.  Scoring the
    # ring against intruders alone left 49% of it Ambiguous: the acinar and ductal cells had
    # nowhere to go.
    def _call(rows, refs):
        if not rows.any():
            return None
        Xr = Xs[rows]
        Xr_nuc = Xs_nuc[rows] if Xs_nuc is not None else None
        nuc_ref = np.array([nuc_by_ref.get(r, False) for r in refs], dtype=bool)
        S, allowed = SC.score_and_gate(Xr, Xr_nuc, Rm.loc[have, refs].to_numpy("float32"),
                                       have, reqs, refs, require_min, nuc_ref)
        o = SC.call_from_scores(S, refs, allowed=allowed)
        o.index = np.flatnonzero(rows)
        return o

    parts = [x for x in (_call(core, intruder), _call(~core, usable)) if x is not None]
    out = pd.concat(parts).sort_index().reset_index(drop=True)

    # Inside Islet_N, a cell that shows no intruder is endocrine BY THE PRIOR, and is
    # subtyped by hormone rather than left Ambiguous.
    sub = subtype_endocrine(Xs, have, require_min)
    fell_through = core & (out["call"] == "Ambiguous").to_numpy()
    for c in ("call", "runner_up"):
        out[c] = np.where(fell_through, sub[c].to_numpy(dtype=object),
                          out[c].to_numpy(dtype=object))
    for c in ("score", "margin"):
        out[c] = np.where(fell_through, sub[c].to_numpy("float64"),
                          out[c].to_numpy("float64"))
    out["reason"] = np.where(fell_through, "endocrine_by_islet_prior",
                             out.get("reason", pd.Series([""] * len(out))).to_numpy(dtype=object))

    out.insert(0, "object_id", df["object_id"].to_numpy())
    out.insert(0, "donor_id", donor)
    out["cell_region"] = df["cell_region"].to_numpy()
    out["islet_num"] = df["islet_num"].to_numpy()
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--require-min", type=float, default=S3.REQUIRE_MIN)
    ap.add_argument("--no-nuclear-ratio", action="store_true")
    args = ap.parse_args()

    donors = args.donors or panel.cohort()["donor_id"].tolist()
    flat = R.flatten(R.load_rules())
    Rm, _ = R.reference_matrix(flat, panel.markers(include_variant=True), S3.ENCODING,
                               strict=False)
    Rm = Rm.loc[(Rm != 0).any(axis=1)]
    reqs = R.requirements(R.load_rules())
    ratio = not args.no_nuclear_ratio
    intruders = intruder_references(R.load_rules(), Rm.columns)
    excluded = [c for c in Rm.columns
                if c not in intruders and not str(c).startswith(ENDOCRINE_PREFIX)]
    print(f"Islet residents — {len(intruders)} references may displace the endocrine prior")
    print(f"  excluded as implied by it: {excluded}")

    scale = S3.build_scale(donors, list(Rm.index), ratio=False)
    OUT.mkdir(parents=True, exist_ok=True)
    frames = []
    for d in donors:
        t = time.time()
        c = run_donor(d, Rm, scale, reqs, args.require_min, ratio, intruders)
        if c.empty:
            print(f"  {d}: no islet cells"); continue
        p = OUT / f"donor_id={d}"; p.mkdir(parents=True, exist_ok=True)
        c.to_parquet(p / "data.parquet", index=False)
        frames.append(c)
        print(f"  {d}: {len(c):,} islet-associated cells "
              f"({(c.cell_region=='core').sum():,} core, {(c.cell_region=='peri').sum():,} "
              f"peri)   {time.time()-t:.0f}s", flush=True)
    calls = pd.concat(frames, ignore_index=True)

    for region in ("core", "peri"):
        g = calls[calls.cell_region == region]
        print(f"\n--- inside Islet_N ({region}), {len(g):,} cells ---"
              if region == "core" else
              f"\n--- the 20 um ring ({region}), {len(g):,} cells — the study's question ---")
        v = g["call"].value_counts(normalize=True)
        for c, f in v.head(14).items():
            print(f"    {c:<38}{f:7.2%}")
    calls.to_parquet(panel.OUT / "islet_residents_all.parquet", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
