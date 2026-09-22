#!/usr/bin/env python3
"""Phase 1: is the undecided islet-core population two cell types, and which cells are which?

WHY THE TWO OBVIOUS TESTS BOTH FAILED, AND WHAT REPLACED THEM.

Attempt 1 -- compare a 1- and a 2-component Gaussian mixture by BIC. It said "two components"
for the undecided cells and ALSO for the reference population, which is endocrine by
construction. A statistic that fires on its own control discriminates nothing.

Attempt 2 -- a parametric bootstrap likelihood-ratio test, simulating the null from the fitted
1-component model. It rejected on the control too, and the reason is measurable: no feature here
is remotely Gaussian. In donor 115's reference, skew runs to 2.61 and excess kurtosis to 18.7;
Shapiro rejects normality at p < 1e-17 on every one. So the observed LR was large because the
data are skewed, while the simulated null was drawn from a genuinely Gaussian model and could
never be. The test had power against "not Gaussian", not against "two populations".

THE FIX IS TO STOP SIMULATING THE NULL AND START SAMPLING IT.

The null here is available as real data: the reference population is one endocrine population,
and it carries the same skew, the same bounds, the same feature correlations, and the same
measurement noise as the target. So the null distribution of the statistic is built by drawing
SUBSAMPLES OF THE REFERENCE at the target's own size and computing the identical statistic on
each. Non-Gaussianity then appears on both sides of the comparison and cancels.

  statistic   LR = 2 * (loglik of a 2-component fit - loglik of a 1-component fit), the same
              quantity as before -- what changes is only what it is compared against.
  null        B subsamples of reference HALF A, each of size n.
  control     one subsample of reference HALF B, size n. Disjoint from the half the null was
              built on, so the control is not scoring itself. It MUST NOT reject.
  target      the undecided cells, size n.
  p           fraction of the null at or above the observed statistic.

The reference pools alpha, beta and delta cells, which need not share a nucleus shape. That
makes the null WIDER than a single cell type would be, so the test is conservative: it is
harder to reject, and a rejection is not explained by endocrine subtype heterogeneity.

Matched n throughout, because LR grows with sample size -- comparing a statistic computed on
9,000 reference cells against one computed on 1,500 undecided cells would measure nothing but
the count.

PER-CELL OUTPUT. Where the target rejects and the control does not, every undecided cell gets
the posterior probability of the component whose centroid lies nearer the reference centroid --
P(endocrine-like). Written per donor. It is a probability, not a call: nothing here assigns a
cell_class and the latch does not move.

ISLET SIZE IS A STRATUM, NOT A NUISANCE. The endocrine shape difference is reported to appear
in islets of more than about ten cells; a two-cell classifier fragment has no architecture to
impose one. Pooling both dilutes exactly the thing being looked for, so the test runs SEPARATELY
in each size stratum, and each stratum is calibrated against ITS OWN hormone-positive cells --
a large islet's endocrine cells are the only fair yardstick for a large islet's uncertain ones.
The small stratum is not discarded: it is the contrast that says whether a positive result in
the large one is about islet architecture or about the whole cohort.

RUN FOR THE FIRST TIME 2026-09-18, over 20 donors, to `data/development/typing/islet_identity/`
(`mixture_test.parquet`, `mixture_test_per_feature.parquet`, `undecided_posterior.parquet`).
In the LARGE stratum, of the 20 donor rows: 10 `displaced subpopulation`, 3 `two populations`,
1 `one population`, 5 UNFIT because the CONTROL rejected, 1 too few cells -- so 14 testable.  Per
feature, against a Bonferroni bar of 0.0083 over the six features, CD99 separates more strongly
than any shape feature and `nucleus_elongation` is well behind it.

WHAT THAT DOES AND DOES NOT SETTLE.  It is evidence about ONE threshold against a cohort-pooled
null within each stratum.  The maintainer's standing point is that a single stratification is not
expected to hold across donors whose islets differ in shape by disease, so this is not the last
word on the elongation route in `subpanel.py`, whose weight remains OPEN.  Five donors' large
strata are UNFIT because their own control rejected, and a statistic that fires on its own
control discriminates nothing -- those five are not evidence either way.

DROPPED DRAWS ARE COUNTED (2026-09-21).  A null draw whose mixture fit raised, or returned a
non-finite statistic, is left out of the null, and the p-value is the fraction of the SURVIVING
draws at or above the observed statistic -- that has not changed.  What the 2026-09-18 outputs
do not carry is how many draws each p-value lost: a null of 300 requested draws could be 30 real
ones and read the same.  `n_draws_requested` and `n_dropped_draws` now sit beside every p-value
in `mixture_test.parquet`, `mixture_test_per_feature.parquet` and `undecided_posterior.parquet`.
The 2026-09-18 files predate the columns and have NOT been re-run.

AN EMPTY NULL IS NOT A VERDICT.  When every draw drops, the null has no survivors, every p-value
is NaN, and NaN compares False at each step of the verdict chain -- which used to fall through
to `one population`.  Such a stratum is now `UNFIT: null empty`, beside the existing
`UNFIT: control rejects`: it is not testable and is counted with the other UNFITs, not as
evidence of anything.  The statistic and the p-value are unchanged; only the verdict is.

A STATISTIC THAT IS NOT FINITE IS NOT A VERDICT EITHER (2026-09-22).  The null can be intact and
the CONTROL's or the TARGET's own fit still return a non-finite statistic -- a zero pooled SD
gives a NaN separation -- for which `pvalue` returns NaN, and that NaN fell through the same
chain to `one population` (or, on the control side, past the `control rejects` check it could
not fail).  Such a stratum is now `UNFIT: statistic not finite`, counted with the other UNFITs.
The statistic and the p-value are reported as they are, NaN; only the verdict changes.

    python core/phenocycler/scripts/typing/islet_identity_mixture.py [--donors 115] [--null 300]
"""
from __future__ import annotations

import argparse
import pathlib
import sys
from typing import NamedTuple

import numpy as np
import pandas as pd
from sklearn.mixture import GaussianMixture

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import subpanel_io as S
from islet_identity import (ISLET_IDENTITY_IN_FLUX, MARKER, NUCLEUS_SHAPE, load)


def fit_statistics(x: np.ndarray, seed: int) -> tuple[float, float]:
    """(LR, separation) for a sample of n rows.

    LR = 2 * (loglik_2 - loglik_1), in nats.

    SEPARATION is the distance between the two fitted component centroids in pooled-SD units.
    It exists because LR alone cannot tell a mixture from a skew: a two-component model fitted
    to one skewed population puts a wide component and a narrow one at nearly the SAME place
    and still gains a large likelihood. That is a curve-fit, not two cell types. Two genuinely
    distinct populations show up as two centroids that are APART. Both statistics are calibrated
    against the same reference null, so "apart" is measured against how far apart a fit to one
    real endocrine population drifts, not against zero.
    """
    if x.ndim == 1:
        x = x[:, None]
    g1 = GaussianMixture(1, random_state=seed, covariance_type="full").fit(x)
    g2 = GaussianMixture(2, random_state=seed, n_init=3, covariance_type="full",
                         reg_covar=1e-6).fit(x)
    lr = float(2.0 * (g2.score(x) - g1.score(x)) * len(x))
    d = float(np.linalg.norm(g2.means_[0] - g2.means_[1]))
    pooled = float(np.sqrt(sum(w * np.trace(c) for w, c in zip(g2.weights_, g2.covariances_))))
    return lr, (d / pooled if pooled > 0 else np.nan)


class NullDraws(NamedTuple):
    """What `reference_null` hands back: the surviving draws, and the ledger of what was lost."""
    values: np.ndarray            #: (n_surviving, 2) of [LR, separation]
    n_draws_requested: int
    n_dropped_draws: int          #: fits that raised, or returned a non-finite statistic


def reference_null(pool: np.ndarray, n: int, draws: int, seed: int, fit=None) -> NullDraws:
    """[LR, separation] on ``draws`` subsamples of size n WITHOUT replacement, minus the drops.

    Without replacement, deliberately: resampling with replacement duplicates points, and a
    mixture fit will happily put a tight component on a duplicated cluster, inflating the null's
    upper tail in a way that has nothing to do with real structure.

    A draw is DROPPED when its fit raises or returns a non-finite statistic (a zero pooled SD
    gives a NaN separation).  The p-value has always been taken over the survivors and still
    is; since 2026-09-21 the drops are COUNTED and carried beside every p-value, because a
    silently thinned null reads exactly like a full one.  ``fit`` is the statistic function --
    `fit_statistics` unless a test hands in a stand-in.
    """
    fit = fit_statistics if fit is None else fit
    rng = np.random.default_rng(seed)
    out = np.full((draws, 2), np.nan)
    for b in range(draws):
        idx = rng.choice(len(pool), n, replace=False)
        try:
            out[b] = fit(pool[idx], seed + b + 1)
        except Exception:
            pass
    kept = np.isfinite(out).all(axis=1)
    return NullDraws(out[kept], int(draws), int(draws - kept.sum()))


def pvalue(null: np.ndarray, observed: float) -> float:
    """One-sided: how often one real endocrine population reaches this statistic by itself."""
    if not len(null) or not np.isfinite(observed):
        return np.nan
    return float((np.sum(null >= observed) + 1) / (len(null) + 1))


def test_population(ref: pd.DataFrame, und: pd.DataFrame, space: list, donor: str,
                    stratum: str, draws: int, draws_1d: int, seed: int):
    """Run the reference-calibrated test on one donor x stratum. Returns (row, per-feature, U, g2)."""
    rng = np.random.default_rng(seed)
    half = rng.permutation(len(ref))
    A, B = half[: len(half) // 2], half[len(half) // 2:]
    n = int(min(len(und), len(A), len(B)))
    base = {"donor_id": donor, "stratum": stratum, "n_reference": len(ref),
            "n_undecided": len(und), "n_matched": n}
    if n < 200:
        return {**base, "verdict": f"too few (n={n})"}, [], None, None

    # Centre and scale by the REFERENCE's moments, so the target cannot set the scale of the
    # very spread the test is looking for.
    mu = ref[space].to_numpy(float).mean(axis=0)
    sd = ref[space].to_numpy(float).std(axis=0)
    sd[sd <= 0] = 1.0
    R = (ref[space].to_numpy(float) - mu) / sd
    U = (und[space].to_numpy(float) - mu) / sd
    ctrl = R[B][rng.choice(len(B), n, replace=False)]
    targ = U[rng.choice(len(U), n, replace=False)] if len(U) > n else U

    null_draws = reference_null(R[A], n, draws, seed)
    null = null_draws.values
    (lr_c, sep_c) = fit_statistics(ctrl, seed)
    (lr_t, sep_t) = fit_statistics(targ, seed)
    p_lr_c, p_lr_t = pvalue(null[:, 0], lr_c), pvalue(null[:, 0], lr_t)
    p_sep_c, p_sep_t = pvalue(null[:, 1], sep_c), pvalue(null[:, 1], sep_t)

    # Three outcomes, not two, because the two statistics can and do disagree. Collapsing them
    # into one yes/no would hide which kind of evidence is actually present.
    #
    # And three UNFITs before them.  A null with no surviving draw gives four NaN p-values; a
    # control or target whose own fit returned a non-finite statistic gives NaN for its two.
    # A NaN is False at every comparison below -- the fall-through would be "one population"
    # for a stratum that was never tested.
    if not len(null):
        verdict = "UNFIT: null empty"
    elif not np.isfinite([p_lr_c, p_sep_c, p_lr_t, p_sep_t]).all():
        verdict = "UNFIT: statistic not finite"
    elif p_lr_c < 0.05 or p_sep_c < 0.05:
        verdict = "UNFIT: control rejects"
    elif p_lr_t < 0.05 and p_sep_t < 0.05:
        verdict = "two populations"
    elif p_sep_t < 0.05:
        verdict = "displaced subpopulation"
    else:
        verdict = "one population"

    row = {**base, "n_draws_requested": null_draws.n_draws_requested,
           "n_dropped_draws": null_draws.n_dropped_draws,
           "null_lr_median": float(np.median(null[:, 0])) if len(null) else np.nan,
           "null_sep_median": float(np.median(null[:, 1])) if len(null) else np.nan,
           "lr_control": lr_c, "p_lr_control": p_lr_c,
           "sep_control": sep_c, "p_sep_control": p_sep_c,
           "lr_target": lr_t, "p_lr_target": p_lr_t,
           "sep_target": sep_t, "p_sep_target": p_sep_t, "verdict": verdict}

    # WHICH FEATURE carries it, same null, one dimension at a time.
    feats = []
    for j, feat in enumerate(space):
        nj_draws = reference_null(R[A][:, j], n, draws_1d, seed)
        nj = nj_draws.values
        lr_tj, sep_tj = fit_statistics(targ[:, j], seed)
        lr_cj, sep_cj = fit_statistics(ctrl[:, j], seed)
        feats.append({"donor_id": donor, "stratum": stratum, "feature": feat,
                      "n_draws_requested": nj_draws.n_draws_requested,
                      "n_dropped_draws": nj_draws.n_dropped_draws,
                      "p_lr_control": pvalue(nj[:, 0], lr_cj),
                      "p_sep_control": pvalue(nj[:, 1], sep_cj),
                      "p_lr_target": pvalue(nj[:, 0], lr_tj),
                      "p_sep_target": pvalue(nj[:, 1], sep_tj),
                      "sep_target": sep_tj, "sep_control": sep_cj,
                      "null_sep_median": float(np.median(nj[:, 1])) if len(nj) else np.nan})

    g2 = None
    if verdict in ("two populations", "displaced subpopulation"):
        g2 = GaussianMixture(2, random_state=seed, n_init=5, covariance_type="full").fit(U)
    return row, feats, U, g2


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--donors", nargs="*", default=None)
    ap.add_argument("--null", type=int, default=300, help="null subsamples, joint space")
    ap.add_argument("--null-1d", type=int, default=1000, help="null subsamples, per feature")
    ap.add_argument("--min-islet-cells", type=int, default=10,
                    help="islets STRICTLY LARGER than this form the large stratum")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    available = sorted(p.name.split("=")[1] for p in (S.OUT / "cells").glob("donor_id=*"))
    donors = [d for d in (args.donors or available) if d not in ISLET_IDENTITY_IN_FLUX]
    out_dir = S.OUT / "islet_identity"
    out_dir.mkdir(parents=True, exist_ok=True)
    k = args.min_islet_cells

    print("null = subsamples of the reference itself; control = a DISJOINT reference half")
    print(f"joint space {args.null} draws, per feature {args.null_1d} draws")
    print(f"strata: islets of >{k} cells, and fragments of <={k} cells, each with its own "
          f"reference\n")

    summary, per_feature, posteriors = [], [], []
    for donor in donors:
        df = load(donor)
        if df is None:
            print(f"  {donor}: inputs incomplete"); continue
        space = [MARKER] + [f for f in NUCLEUS_SHAPE if f in df.columns]
        df = df.dropna(subset=space)
        print(f"{donor}  islet-core cells {len(df):,}  "
              f"({100 * float((df['islet_cells'] > k).mean()):.1f}% in islets of >{k} cells)")

        for stratum, mask in (("large", df["islet_cells"] > k), ("fragment", df["islet_cells"] <= k)):
            sub = df[mask]
            ref, und = sub[sub["reference"]], sub[sub["undecided"]]
            row, feats, U, g2 = test_population(ref, und, space, donor, stratum,
                                                args.null, args.null_1d, args.seed)
            summary.append(row); per_feature.extend(feats)
            if "n_matched" in row and row.get("verdict", "").startswith("too few"):
                print(f"   {stratum:9s} reference {len(ref):,} | undecided {len(und):,}"
                      f"  -> {row['verdict']}")
                continue
            print(f"   {stratum:9s} reference {len(ref):,} | undecided {len(und):,} | n {row['n_matched']:,}")
            print(f"      control  LR p {row['p_lr_control']:.3f}  separation {row['sep_control']:.3f} "
                  f"p {row['p_sep_control']:.3f}")
            print(f"      target   LR p {row['p_lr_target']:.3f}  separation {row['sep_target']:.3f} "
                  f"p {row['p_sep_target']:.3f}   -> {row['verdict']}")
            lost_1d = sum(f["n_dropped_draws"] for f in feats)
            if row["n_dropped_draws"] or lost_1d:
                print(f"      dropped  {row['n_dropped_draws']} of {row['n_draws_requested']} joint "
                      f"null draws; {lost_1d} of {sum(f['n_draws_requested'] for f in feats)} "
                      f"per-feature draws (fit raised or non-finite)")
            if g2 is not None:
                # After standardising by the reference, the reference centroid IS the origin, so
                # the endocrine-like component is the one with the smaller centroid norm.
                endo = int(np.argmin(np.linalg.norm(g2.means_, axis=1)))
                post = g2.predict_proba(U)[:, endo]
                posteriors.append(pd.DataFrame({
                    "donor_id": donor, "stratum": stratum,
                    "object_id": und["object_id"].to_numpy(),
                    "islet_num": und["islet_num"].to_numpy(),
                    "p_endocrine_like": post,
                    "evidence": row["verdict"],
                    "n_draws_requested": row["n_draws_requested"],
                    "n_dropped_draws": row["n_dropped_draws"]}))
                print(f"      -> {float((post > 0.5).mean()):.1%} of undecided cells score "
                      f"endocrine-like (component weight {g2.weights_[endo]:.3f})")

    if not summary:
        return 0
    frame, feats = pd.DataFrame(summary), pd.DataFrame(per_feature)
    # A "too few" row asked for no draws: its ledger is unavailable, not 0, so the columns stay
    # nullable integers rather than floats with a NaN in them.
    for col in ("n_draws_requested", "n_dropped_draws"):
        if col in frame.columns:
            frame[col] = frame[col].astype("Int64")
    frame.to_parquet(out_dir / "mixture_test.parquet", index=False)
    feats.to_parquet(out_dir / "mixture_test_per_feature.parquet", index=False)
    if posteriors:
        pd.concat(posteriors, ignore_index=True).to_parquet(
            out_dir / "undecided_posterior.parquet", index=False)

    pd.set_option("display.width", 200)
    print("\n" + "=" * 80)
    for stratum in ("large", "fragment"):
        sl = frame[frame["stratum"] == stratum]
        tested = sl[sl["verdict"].isin(["two populations", "displaced subpopulation",
                                        "one population"])]
        print(f"\nISLETS OF >{k} CELLS" if stratum == "large" else
              f"\nFRAGMENTS OF <={k} CELLS  (the contrast)")
        print(f"  testable {len(tested)} of {len(sl)} donors; " +
              ", ".join(f"{v}: {c}" for v, c in sl['verdict'].value_counts().items()))
        cols = ["donor_id", "n_matched", "n_dropped_draws", "p_lr_control", "p_sep_control",
                "p_lr_target", "sep_target", "p_sep_target", "verdict"]
        print(sl[[c for c in cols if c in sl.columns]].to_string(index=False))

    fl = feats[feats["stratum"] == "large"]
    if len(fl):
        bon = 0.05 / fl["feature"].nunique()
        fl = fl[(fl["p_lr_control"] >= 0.05) & (fl["p_sep_control"] >= 0.05)]
        tab = (fl.assign(fires=(fl["p_sep_target"] < bon))
               .groupby("feature")
               .agg(donors_firing=("fires", "sum"), donors_fit=("fires", "size"),
                    median_sep_target=("sep_target", "median"),
                    median_sep_null=("null_sep_median", "median"),
                    median_p_sep=("p_sep_target", "median"),
                    median_p_lr=("p_lr_target", "median"))
               .sort_values(["donors_firing", "median_p_sep"], ascending=[False, True]))
        print(f"\nPER FEATURE, islets of >{k} cells, donors whose control cleared both statistics")
        print(f"  Bonferroni over {feats['feature'].nunique()} features: p < {bon:.4f}")
        print(tab.to_string())
    print(f"\n-> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
