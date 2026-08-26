"""Why do a quarter of the peri-islet cells satisfy no rule?

``Ambiguous`` here means one specific thing: the cell showed none of the marker combinations
any reference demands.  That is honest, but it is not an answer, and at ~25% of the 20 um ring
-- the tissue this study exists to measure -- it is too large to leave unexamined.

Four questions, one panel each:

* **Where are they?**  A donor effect, a region effect, or uniform.
* **Are they simply dim?**  If the ambiguous cells are just cells with nothing above
  background, that is a segmentation and staining story, not a rule story.
* **What DO they express?**  Their marker profile against the cells that were called.
* **What blocked them?**  For every reference, how many ambiguous cells met all its
  requirements but one, and which marker was the single blocker.  This is the actionable
  panel: it says "N cells would be Imm-macrophage if CD68 cleared" rather than "N cells are
  unexplained".

    python scripts/celltype/figures_ambiguous.py
"""
from __future__ import annotations

import pathlib
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel                             # noqa: E402
import rules as R                        # noqa: E402
import score as SC                       # noqa: E402
import run_step3 as S3                   # noqa: E402
import run_islet as RI                   # noqa: E402

AMB, CORE = "Ambiguous", "core"


def near_miss(donor: str, Rm: pd.DataFrame, scale: pd.DataFrame, reqs: dict,
              require_min: float, ratio: bool, intruders: list) -> pd.DataFrame:
    """For ambiguous cells: which single requirement, if met, would have admitted them.

    A cell short by one marker is a different object from one short by five.  The first is a
    threshold question about that marker; the second is a cell the panel has nothing to say
    about, and reporting them together would hide which is which.
    """
    keep = ["object_id", "qc_analysis_eligible", "cell_region", "islet_num"]
    df, have, _ = S3._read_donor(donor, list(Rm.index), keep, ratio)
    df = df[df["qc_analysis_eligible"] & df["cell_region"].isin(("core", "peri"))]
    if df.empty:
        return pd.DataFrame()
    df = df.reset_index(drop=True)
    Xs = SC.apply_scale(df[have].to_numpy("float32"), have, scale)
    hot = Xs >= require_min
    idx = {m: i for i, m in enumerate(have)}
    core = (df["cell_region"] == CORE).to_numpy()

    rows = []
    for j, ph in enumerate(Rm.columns):
        if ph not in reqs:
            continue
        # Intruders compete everywhere; the rest only in the ring, where there is no prior.
        in_play = np.ones(len(Xs), bool) if ph in intruders else ~core
        all_of, groups = reqs[ph]
        blockers = [m for m in all_of if m in idx]
        if len(blockers) != len(all_of):
            continue                              # needs a marker this donor lacks
        met_groups = np.ones(len(Xs), bool)
        for g in groups:
            cols = [idx[m] for m in g if m in idx]
            if not cols:
                met_groups[:] = False
                break
            met_groups &= hot[:, cols].any(axis=1)
        M = hot[:, [idx[m] for m in blockers]] if blockers else np.ones((len(Xs), 0), bool)
        n_short = (~M).sum(axis=1) + (~met_groups).astype(int)
        one_short = (n_short == 1) & in_play
        if not one_short.any():
            continue
        which = np.where((~M).any(axis=1), np.asarray(blockers, dtype=object)[
            np.argmax(~M, axis=1)] if blockers else "", "<any-of group>")
        for m in pd.unique(which[one_short]):
            sel = one_short & (which == m)
            rows.append({"donor_id": donor, "reference": ph, "blocker": str(m),
                         "n": int(sel.sum()),
                         "n_core": int((sel & core).sum()),
                         "n_peri": int((sel & ~core).sum())})
    return pd.DataFrame(rows)


def figure(calls: pd.DataFrame, prof: pd.DataFrame, miss: pd.DataFrame,
           out: pathlib.Path) -> pathlib.Path:
    fig, ax = plt.subplots(2, 2, figsize=(17, 11), dpi=140)

    # (A) where they are
    a = ax[0, 0]
    t = (calls.assign(amb=calls["call"].eq(AMB))
              .groupby(["donor_id", "cell_region"])["amb"].mean().unstack() * 100)
    t = t.sort_values("peri")
    x = np.arange(len(t))
    a.barh(x - 0.2, t["core"], height=0.4, color="#1f77b4", label="inside Islet_N")
    a.barh(x + 0.2, t["peri"], height=0.4, color="#ff7f0e", label="20 um ring")
    a.set_yticks(x); a.set_yticklabels(t.index, fontsize=7)
    a.set_xlabel("% of cells satisfying no rule")
    a.set_title("A. Where the ambiguous cells are", fontweight="bold")
    a.legend(fontsize=8, frameon=False)

    # (B) are they simply dim?
    b = ax[0, 1]
    for lab, grp, colour in (("called", "called", "#2ca02c"),
                             ("ambiguous", "amb", "#d62728")):
        v = prof.loc[prof["group"] == grp, "max_scaled"].to_numpy()
        if len(v):
            b.hist(v, bins=60, range=(0, 1), histtype="step", lw=2, color=colour, label=lab,
                   density=True)
    b.axvline(S3.REQUIRE_MIN, color="k", ls="--", lw=1)
    b.text(S3.REQUIRE_MIN, b.get_ylim()[1] * 0.9, "  require_min", fontsize=8)
    b.set_xlabel("brightest marker in the cell, on the [background, robust top] scale")
    b.set_ylabel("density")
    b.set_title("B. Are they simply dim?", fontweight="bold")
    b.legend(fontsize=8, frameon=False)

    # (C) what do they express?
    c = ax[1, 0]
    m = prof.drop(columns=["max_scaled"]).groupby("group").median(numeric_only=True).T
    m = m.dropna().sort_values("amb")
    y = np.arange(len(m))
    c.barh(y - 0.2, m["called"], height=0.4, color="#2ca02c", label="called")
    c.barh(y + 0.2, m["amb"], height=0.4, color="#d62728", label="ambiguous")
    c.set_yticks(y); c.set_yticklabels(m.index, fontsize=6)
    c.set_xlabel("median scaled value")
    c.set_title("C. What the ambiguous cells do express", fontweight="bold")
    c.legend(fontsize=8, frameon=False)

    # (D) the single blocker
    d = ax[1, 1]
    if len(miss):
        g = (miss.groupby(["reference", "blocker"])["n"].sum()
                 .sort_values(ascending=False).head(18)[::-1])
        lbl = [f"{r}  ←  {b}" for r, b in g.index]
        d.barh(np.arange(len(g)), g.to_numpy(), color="#9467bd")
        d.set_yticks(np.arange(len(g))); d.set_yticklabels(lbl, fontsize=7)
        d.set_xlabel("ambiguous cells that met every requirement but this one")
    else:
        d.text(0.5, 0.5, "no cell was short by exactly one requirement",
               ha="center", transform=d.transAxes)
    d.set_title("D. What blocked them — one requirement short", fontweight="bold")

    fig.suptitle(f"Cells satisfying no rule — {len(calls):,} islet-associated cells, "
                 f"{calls['call'].eq(AMB).mean():.1%} ambiguous", fontweight="bold", y=1.00)
    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, bbox_inches="tight")
    plt.close(fig)
    return out


def main() -> int:
    calls = pd.read_parquet(panel.OUT / "islet_residents_all.parquet")
    flat = R.flatten(R.load_rules())
    Rm, _ = R.reference_matrix(flat, panel.markers(include_variant=True), S3.ENCODING,
                               strict=False)
    Rm = Rm.loc[(Rm != 0).any(axis=1)]
    reqs = R.requirements(R.load_rules())
    intruders = RI.intruder_references(R.load_rules(), Rm.columns)
    donors = panel.cohort()["donor_id"].tolist()
    scale = S3.build_scale(donors, list(Rm.index), ratio=S3.NUCLEAR_RATIO)

    prof, miss = [], []
    for d in donors:
        keep = ["object_id", "qc_analysis_eligible", "cell_region"]
        df, have, _ = S3._read_donor(d, list(Rm.index), keep, S3.NUCLEAR_RATIO)
        df = df[df["qc_analysis_eligible"] & df["cell_region"].isin(("core", "peri"))]
        if df.empty:
            continue
        df = df.reset_index(drop=True)
        Xs = SC.apply_scale(df[have].to_numpy("float32"), have, scale)
        lab = calls.loc[calls.donor_id == d].set_index("object_id")["call"]
        g = np.where(df["object_id"].map(lab).eq(AMB).to_numpy(), "amb", "called")
        t = pd.DataFrame(Xs, columns=have)
        t["group"] = g
        t["max_scaled"] = Xs.max(axis=1)
        prof.append(t)
        miss.append(near_miss(d, Rm, scale, reqs, S3.REQUIRE_MIN, S3.NUCLEAR_RATIO, intruders))
        print(f"  {d}: {(g == 'amb').sum():,} ambiguous of {len(g):,}", flush=True)

    prof = pd.concat(prof, ignore_index=True)
    miss = pd.concat([m for m in miss if len(m)], ignore_index=True) if miss else pd.DataFrame()
    miss.to_parquet(panel.OUT / "ambiguous_near_miss.parquet", index=False)
    out = figure(calls, prof, miss, panel.OUT / "figures_ambiguous" / "ambiguous.png")
    print(f"\nwrote {out}")

    if len(miss):
        print("\nambiguous cells one requirement short, pooled over donors:")
        g = miss.groupby(["reference", "blocker"])[["n", "n_core", "n_peri"]].sum()
        print(g.sort_values("n", ascending=False).head(15).to_string())
        print(f"\n{miss['n'].sum():,} of {int(calls['call'].eq(AMB).sum()):,} ambiguous cells "
              f"({miss['n'].sum() / max(calls['call'].eq(AMB).sum(), 1):.1%}) were short by "
              f"exactly one requirement")
    return 0


if __name__ == "__main__":
    sys.exit(main())
