#!/usr/bin/env python3
"""
Score the pipeline against the manual gold standard.

Reports three things, in this order of authority:

  1. INTER-ANNOTATOR agreement (Cohen's kappa, per-compartment agreement, unsure rate). This bounds
     everything else — the pipeline cannot be shown to be more accurate than two humans agree the
     truth is. If kappa is low the label definitions are the problem, not the pipeline, and the
     accuracy numbers below should not be quoted.
  2. PER-COMPARTMENT precision / recall / F1 of the pipeline against the consensus label, both raw
     (over the sample) and POPULATION-WEIGHTED. The sample is stratified and over-samples immune, so
     the raw numbers describe the sample and the weighted numbers describe the cohort. Quote the
     weighted ones.
  3. The confusion matrix, which is what says where the errors go — an Epithelial cell called Other
     and an Epithelial cell called Immune are very different failures for this study.

Consensus rule: the two annotators must agree. Disagreements and "Unsure" are EXCLUDED from accuracy
and reported as a separate count — adjudicating them by majority with two annotators is just picking
one, and adjudicating them against the pipeline would make the measurement circular.

    python scripts/gold_standard/score.py
"""

from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "Phenocycler_Analysis"))
from phenocycler.config import COMPARTMENT_ORDER, OTHER_LABEL, UNRESOLVED_LABEL   # noqa: E402

LABELS = COMPARTMENT_ORDER + [OTHER_LABEL]
VOCAB = set(LABELS) | {"Unsure"}


def load_annotations(completed: Path) -> pd.DataFrame:
    files = sorted(glob.glob(str(completed / "annotations_*.csv")))
    if not files:
        raise SystemExit(f"[err] no annotations_*.csv under {completed} — nothing has been annotated yet")
    d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    bad = set(d["compartment"].dropna().unique()) - VOCAB - {""}
    if bad:
        raise SystemExit(f"[err] labels outside the vocabulary: {sorted(bad)}")
    return d


def consensus(d: pd.DataFrame) -> tuple[pd.DataFrame, dict]:
    """One row per review_id with the agreed label; plus the agreement diagnostics."""
    wide = d.pivot_table(index=["gs_review_id", "object_id"], columns="annotator",
                         values="compartment", aggfunc="first").reset_index()
    who = [c for c in wide.columns if c not in ("gs_review_id", "object_id")]
    if len(who) < 2:
        raise SystemExit(f"[err] need two annotators, found {who}")
    a, b = wide[who[0]].fillna(""), wide[who[1]].fillna("")
    both = (a != "") & (b != "")
    unsure = both & ((a == "Unsure") | (b == "Unsure"))
    scorable = both & ~unsure
    agree = scorable & (a == b)

    diag = {
        "n_rows": int(len(wide)),
        "n_both_annotated": int(both.sum()),
        "n_unsure_either": int(unsure.sum()),
        "n_scorable": int(scorable.sum()),
        "n_agree": int(agree.sum()),
        "raw_agreement": float(agree.sum() / max(scorable.sum(), 1)),
        "cohens_kappa": _kappa(a[scorable], b[scorable]),
        "annotators": who,
    }
    out = wide.loc[agree, ["gs_review_id", "object_id"]].copy()
    out["truth"] = a[agree].to_numpy()
    return out, diag


def _kappa(a: pd.Series, b: pd.Series) -> float:
    cats = sorted(set(a) | set(b))
    if not cats:
        return float("nan")
    m = pd.crosstab(pd.Categorical(a, cats), pd.Categorical(b, cats), dropna=False).to_numpy(float)
    n = m.sum()
    po = np.trace(m) / n
    pe = float((m.sum(0) / n * (m.sum(1) / n)).sum())
    return float((po - pe) / (1 - pe)) if pe < 1 else float("nan")


def pipeline_calls(gated_dir: Path, cells_dir: Path, epithelial_gate: str, object_ids) -> pd.DataFrame:
    """Recompute the pipeline compartment for the reviewed cells only, from a gated-parquet directory.

    Defaults to the PRODUCTION gated output and the gate set in ``config.COMPARTMENT_GATES`` as it
    actually ships. That matters: the anatomy-check arms under ``data/_arm/`` were built with the old
    7-marker gate set, so scoring against one measures a configuration that no longer exists — and
    since ``lineage`` fails loudly on a missing ``{marker}_state``, it does not even run."""
    from phenocycler import lineage
    from phenocycler.config import COMPARTMENT_GATES
    gates = dict(COMPARTMENT_GATES)
    if epithelial_gate == "anatomy-union":
        # The anatomy check's 4-marker comparator. NOT what production ships: the shipped Epithelial
        # gate also carries Keratin_5 and the three hormones, and in islet core INS/GCG do most of the
        # work (docs/GATE5_ROUND_2026-07-27.md, Finding 5). Kept only to reproduce that round.
        gates["Epithelial"] = ["E_cadherin", "Pan_Cytokeratin", "Ker8_18", "EpCAM"]
    want = set(map(str, object_ids))
    files = sorted(glob.glob(str(gated_dir / "donor_id=*" / "*.parquet")))
    if not files:
        raise SystemExit(f"[err] no gated parquets under {gated_dir}")
    rows = []
    for f in files:
        donor = Path(f).parent.name.split("=", 1)[1]
        lineage.COMPARTMENT_GATES = gates
        try:
            out = lineage.assign_donor(donor, f, cells_dir)
        finally:
            lineage.COMPARTMENT_GATES = dict(COMPARTMENT_GATES)
        rows.append(out.loc[out["object_id"].isin(want), ["object_id", "compartment"]])
    return pd.concat(rows, ignore_index=True).rename(columns={"compartment": "predicted"})


def rates(truth, pred, weight, labels) -> pd.DataFrame:
    rows = []
    w = np.asarray(weight, float)
    for lab in labels:
        t, p = (np.asarray(truth) == lab), (np.asarray(pred) == lab)
        tp, fp, fn = float(w[t & p].sum()), float(w[~t & p].sum()), float(w[t & ~p].sum())
        prec = tp / (tp + fp) if tp + fp else np.nan
        rec = tp / (tp + fn) if tp + fn else np.nan
        f1 = 2 * prec * rec / (prec + rec) if prec and rec and (prec + rec) else np.nan
        rows.append(dict(compartment=lab, n_truth=int(t.sum()), n_pred=int(p.sum()),
                         precision=prec, recall=rec, f1=f1))
    return pd.DataFrame(rows)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--bundle", type=Path, default=REPO / "data" / "gold_standard_v1")
    ap.add_argument("--config", type=Path, default=REPO / "config.ini")
    ap.add_argument("--gated-dir", type=Path, default=None,
                    help="gated parquet dir to score (default: the production restore_gated_redsea). "
                         "Point at data/_arm/<arm>/restore_gated_redsea to score an anatomy-check arm "
                         "instead — but those were built with the old 7-marker gate set.")
    ap.add_argument("--cells-dir", type=Path, default=None, help="default: the production cells dir")
    ap.add_argument("--epithelial-gate", choices=("config", "anatomy-union"), default="config",
                    help="'config' = the Epithelial gate as it ships (default); 'anatomy-union' = the "
                         "anatomy check's 4-marker comparator, for reproducing that round only")
    a = ap.parse_args(argv)

    from phenocycler.config import load_config
    cfg = load_config(a.config)
    gated_dir = Path(a.gated_dir) if a.gated_dir else cfg.restore_gated_dir
    cells_dir = Path(a.cells_dir) if a.cells_dir else cfg.cells_dir
    scored = "production" if a.gated_dir is None else str(gated_dir)

    key = pd.read_csv(a.bundle / "sample_key.csv", dtype={"object_id": str})
    ann = load_annotations(a.bundle / "completed")
    truth, diag = consensus(ann)

    print("=" * 72)
    print("1. INTER-ANNOTATOR AGREEMENT (this bounds every number below)")
    print("=" * 72)
    for k, v in diag.items():
        print(f"  {k:20s} {v}")
    if diag["n_scorable"] and diag["cohens_kappa"] < 0.6:
        print("\n  ** kappa < 0.6: the label definitions disagree, not just the calls. Fix the")
        print("     protocol before quoting any accuracy number below. **")

    m = truth.merge(key[["object_id", "donor", "stratum", "weight", "cell_region"]],
                    on="object_id", how="left", validate="1:1")
    pred = pipeline_calls(gated_dir, cells_dir, a.epithelial_gate, m["object_id"])
    m = m.merge(pred, on="object_id", how="left")
    missing = int(m["predicted"].isna().sum())
    if missing:
        print(f"\n  ** {missing} of {len(m)} reviewed cells have no pipeline call in {gated_dir} — "
              "they are counted as Unresolved, which makes recall pessimistic. Check the donor set. **")
    m["predicted"] = m["predicted"].fillna(UNRESOLVED_LABEL)

    print("\n" + "=" * 72)
    print(f"2. PIPELINE ACCURACY — scoring {scored}, Epithelial gate = {a.epithelial_gate}")
    print("=" * 72)
    print("\n-- over the sample (stratified: over-samples immune, describes the SAMPLE) --")
    print(rates(m.truth, m.predicted, np.ones(len(m)), LABELS).round(3).to_string(index=False))
    print("\n-- population-weighted (describes the COHORT; quote these) --")
    weighted = rates(m.truth, m.predicted, m.weight, LABELS)
    print(weighted.round(3).to_string(index=False))

    print("\n" + "=" * 72)
    print("3. CONFUSION (rows = human truth, cols = pipeline)")
    print("=" * 72)
    cm = pd.crosstab(m.truth, m.predicted, dropna=False)
    print(cm.to_string())

    out = a.bundle / "scores"
    out.mkdir(exist_ok=True)
    tag = ("production" if a.gated_dir is None else Path(a.gated_dir).parent.name) + f"_{a.epithelial_gate}"
    weighted.to_csv(out / f"accuracy_{tag}_weighted.csv", index=False)
    cm.to_csv(out / f"confusion_{tag}.csv")
    pd.Series(diag).to_csv(out / "inter_annotator.csv")
    print(f"\n[done] {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
