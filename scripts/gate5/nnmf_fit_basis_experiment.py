#!/usr/bin/env python
"""Fit the NNMF on the DATA instead of on the pre-thresholded arms, and measure what changes.

THE QUESTION
------------
The manuscript's design is for the clustering to FIND the two groups. Lemma 1: for a mutually
exclusive pair, {Y1, Y2, Y3} cannot be written as positive combinations of one another, so a
2-component non-negative factorisation of the pair's intensity profile recovers the split *wherever it
lies in that donor's distribution*. That adaptivity is why the method clusters rather than thresholds,
and it is what should let one method serve markers with very different frequencies and intensity ranges.

This implementation instead hands the model the answer: it picks a training set by per-marker raw Otsu
floors (target-high AND reference-low vs reference-high), samples the two arms to EQUAL SIZE, fits on
that, and predicts the rest. The clustering can then only ratify the split the floors already drew.
Donor 6566 INS is the clean failure: 412,459 cells supplied as "INS-positive examples", of which a few
hundred are (maintainer, from the image).

So: does fitting on the data recover the adaptivity?

WHAT IT RUNS
------------
The Gate-3 representative subset -- 7 donors x the 7 ACCEPTED_PAIRS, 49 donor-pairs -- through the
production `evaluate_locked_pair` three times, changing exactly one thing:

    arms         the shipped basis: balanced, floor-selected, equal-size exclusive arms
    qc_retained  a deterministic sample of cells above at least ONE floor. The floors still drop
                 double-negatives -- the role the published `idx_select` plays -- but no longer decide
                 which cells are which class.
    all          a deterministic sample of every valid cell; floors unused at fit time.

Nothing is written to canonical `data/`. No policy changes: `nnmf_fit_basis` defaults to "arms" and
all 254 tests pass on the default.

USAGE
    ~/miniconda3/envs/scvi-env/bin/python scripts/gate5/nnmf_fit_basis_experiment.py
    ... --donors 6566 --pairs "CD68<-E_cadherin"
"""
from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "Phenocycler_Analysis"))

from phenocycler.config import load_config                                     # noqa: E402
from phenocycler.restore_validation import (                                   # noqa: E402
    ACCEPTED_PAIRS,
    ExpressionReview,
    PairReview,
    evaluate_locked_pair,
    load_expression_reviews,
    load_pair_reviews,
    load_validation_sample,
)

# The Gate-3 full-cell stress pilot subset (`gate3_pilot_v2`), i.e. the maintainer's representative
# donors: 6539 and 6476 plus the weakest-case donor per pair.
PILOT_DONORS = ("6442", "6476", "6505", "6539", "6551", "6566", "6591")
# (nnmf_fit_basis, control_membership). Changing the fit basis alone does not hand the decision to the
# clustering -- the raw floors are ANDed into the membership predicates after it -- so both are swept.
ARMS = (
    ("arms",        "floor_and_component"),   # shipped
    ("qc_retained", "floor_and_component"),
    ("all",         "floor_and_component"),
    ("arms",        "component"),
    ("qc_retained", "component"),             # closest to the manuscript
    ("all",         "component"),
)


def run(donor_df, donor, target, reference, config, reviews, pairs, arms) -> list[dict]:
    rows = []
    for basis, membership, k, loss in arms:
        solver = "cd" if loss == "frobenius" else "mu"
        cfg = replace(config, nnmf_fit_basis=basis, control_membership=membership,
                      n_components=k, nnmf_beta_loss=loss, nnmf_solver=solver)
        try:
            e = evaluate_locked_pair(
                donor_df, donor, target, reference, config=cfg,
                expression_review=reviews.get((donor, target), ExpressionReview.UNREVIEWED),
                pair_review=pairs.get((target, reference), PairReview.UNREVIEWED),
            )
        except Exception as exc:                       # a basis that cannot fit is a result, not a crash
            rows.append({"donor": donor, "target": target, "reference": reference, "basis": basis,
                         "membership": membership, "k": k, "beta_loss": loss,
                         "state": "ERROR", "state_reason": str(exc)[:160]})
            continue
        stats = e.get("target_full_stats") or {}
        n = int(e.get("n_assigned") or 0)
        rows.append({
            "donor": donor, "target": target, "reference": reference, "basis": basis,
            "membership": membership, "k": k, "beta_loss": loss,
            "state": str(e["state"]), "state_reason": str(e.get("state_reason", ""))[:160],
            "divisor": stats.get("divisor"), "manuscript_max": stats.get("maximum"),
            "n_assigned": n,
            "threshold_positive_n": e.get("threshold_positive_n"),
            "call_rate": (e.get("threshold_positive_n") or 0) / n if n else None,
            "target_component_n": e.get("target_component_n"),
            "reference_component_n": e.get("reference_component_n"),
            "reference_control_n": e.get("reference_control_n"),
            "target_anchor_n": e.get("target_anchor_n"),
            "target_fold": e.get("target_fold"),
            "reference_fold": e.get("reference_fold"),
            "divisor_reproducibility": e.get("divisor_reproducibility"),
            "anchor_recovery": e.get("anchor_recovery"),
            "double_high_fraction": e.get("double_high_fraction"),
            "candidate_target_n": e.get("candidate_target_n"),
            "candidate_reference_n": e.get("candidate_reference_n"),
            "target_input_floor": e.get("target_input_floor"),
            "feature_compartments": e.get("feature_compartments"),
        })
    return rows


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--config", type=Path, default=REPO / "config.ini")
    p.add_argument("--screen", type=Path,
                   default=REPO / "data/restore_pair_validation/screen_v12_production"
                                  "/ctrl-reference_only_div-p99")
    p.add_argument("--sample", type=Path, default=None)
    p.add_argument("--donors", nargs="*", default=list(PILOT_DONORS))
    p.add_argument("--pairs", nargs="*", default=None, help='e.g. "CD68<-E_cadherin"')
    p.add_argument("--out", type=Path, default=REPO / "data" / "nnmf_fit_basis" / "results.csv")
    p.add_argument("--basis", nargs="*", default=["arms"], choices=["arms", "qc_retained", "all"])
    p.add_argument("--membership", nargs="*", default=["floor_and_component"],
                   choices=["floor_and_component", "component"])
    p.add_argument("--n-components", nargs="*", type=int, default=[2])
    p.add_argument("--beta-loss", nargs="*", default=["frobenius"],
                   choices=["frobenius", "kullback-leibler"])
    p.add_argument("--expression-reviews", type=Path,
                   default=REPO / "data/restore_pair_validation/expression_reviews.csv",
                   help="defaults to the CANONICAL file; the screen's own copy still carries the "
                        "withdrawn 6566 INS CONFIRMED_ABSENT row")
    a = p.parse_args(argv)

    cfg = load_config(a.config)
    sample = a.sample or (cfg.restore_pair_validation_dir / "qupath_compartment_full_20donor")
    spec = json.loads((a.screen / "method_spec.json").read_text())
    import phenocycler.restore_pair_review as rpr
    config = rpr._config_from_specification(spec)
    reviews, _ = load_expression_reviews(a.expression_reviews)
    pairs, _ = load_pair_reviews(a.screen / "pair_reviews.csv")

    wanted = ACCEPTED_PAIRS if not a.pairs else [
        tuple(s.split("<-")) for s in a.pairs]
    wanted = [(t.strip(), r.strip()) for t, r in wanted]

    arms = [(b, m, k, l) for b in a.basis for m in a.membership
            for k in a.n_components for l in a.beta_loss]
    print(f"[basis] {len(arms)} arms x {len(a.donors)} donors", flush=True)

    rows = []
    for donor in a.donors:
        print(f"[basis] loading donor {donor} …", flush=True)
        donor_df, _ = load_validation_sample(sample, cfg.cells_dir, cfg.cells_redsea_dir, donor)
        for target, reference in wanted:
            if f"raw__{target}" not in donor_df or f"raw__{reference}" not in donor_df:
                continue
            rows += run(donor_df, donor, target, reference, config, reviews, pairs, arms)
            print(f"[basis]   {donor} {target} <- {reference}", flush=True)
        del donor_df

    d = pd.DataFrame(rows)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    d.to_csv(a.out, index=False)

    print(f"\n=== states by basis ===")
    d["arm"] = (d.basis + " / " + d.membership + " / k=" + d.k.astype(str)
                + " / " + d.beta_loss)
    print(pd.crosstab(d.arm, d.state).to_string())
    ok = d[d.state == "VALID"]
    if len(ok):
        print(f"\n=== median over VALID donor-pairs ===")
        print(ok.groupby("arm")[["call_rate", "divisor", "target_fold", "reference_fold",
                                 "divisor_reproducibility"]].median().round(4).to_string())
    print(f"\n-> {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
