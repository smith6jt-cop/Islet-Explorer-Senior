"""Verify that every gene named in a tissue's marker module exists in
the corresponding `panel_audit_{tissue}.csv` and is `exclude_recommended=='no'`.

For each panel (MARKER_DICT, ANTI_OVERRIDE_MARKERS, LINEAGE_CANARY_PANELS,
PANEL_COVERAGE_PANELS), report which genes are absent / excluded, and
warn if a MARKER_DICT entry would drop below 3 panel-present genes
after filtering (the ≥3-genes rule from CLAUDE.md decision #12).
"""
import argparse
import importlib
import os
from pathlib import Path

import pandas as pd

# Two roots, deliberately separate.  ROOT is where the CODE and the tracked tables live;
# DATA_ROOT is where the multi-GB h5ads live and is NEVER inside the repo.  Collapsing them
# -- as the single HiPerGator ROOT did -- would put data/processed/ under analysis/xenium/,
# and the repo-root .gitignore pattern "/data/" is ANCHORED, so it would not be ignored.
ROOT = Path(__file__).resolve().parents[2]                     # analysis/xenium
REPO = Path(__file__).resolve().parents[4]                     # monorepo root
DATA_ROOT = Path(os.environ.get("XENIUM_DATA_ROOT", "/data/xenium"))


def load_audit(tissue):
    p = ROOT / f"panel_audit_{tissue}.csv"
    if not p.exists() and tissue == "pancreas":
        # Pancreas uses unsuffixed filename.
        p = ROOT / "panel_audit.csv"
    df = pd.read_csv(p)
    return p, df


def check_panel(name, panel_dict, audit_df, min_genes=3, kind="scoring"):
    """Report which genes are absent or excluded per cell-type entry.

    `kind` ∈ {scoring, override, canary, coverage} controls the
    minimum-genes warning behavior:
      - scoring: every entry must have ≥ min_genes panel-present.
      - others: warn but don't fail; small panels are OK for canary /
                override / coverage.
    """
    in_panel = set(audit_df["gene"])
    excluded = set(audit_df.loc[audit_df["exclude_recommended"] == "yes", "gene"])
    print(f"\n--- {name} ---")
    issues_critical = 0
    for ct, genes in panel_dict.items():
        absent = [g for g in genes if g not in in_panel]
        flagged = [g for g in genes if g in excluded]
        present = [g for g in genes
                   if g in in_panel and g not in excluded]
        status = "ok"
        if absent:
            status = "absent"
        if flagged:
            status = "flagged"
        if kind == "scoring" and len(present) < min_genes:
            status = "FAIL"
            issues_critical += 1
        if absent or flagged or status != "ok":
            print(f"  {ct:<18} status={status}  present={len(present)}/{len(genes)}")
            if absent:
                print(f"      absent from panel:  {absent}")
            if flagged:
                print(f"      flagged as exclude: {flagged}")
            if kind == "scoring" and len(present) < min_genes:
                print(f"      → DROP from score_genes argmax "
                      f"(only {len(present)} panel-present, need ≥{min_genes})")
    if issues_critical:
        print(f"  *** {issues_critical} panels would be dropped from scoring ***")
    return issues_critical


def check_seed_markers(mod, audit_df, min_genes=2):
    """Validate the scANVI broad-classifier SEED_MARKERS (redesign 2026-06-08):
    keys must equal BROAD_TO_FINE keys; each lineage must have >= min_genes
    panel-present, non-excluded seed genes (decision #12 — no thin seeds, which
    would starve a scANVI class and make a lineage unpredictable)."""
    seeds = getattr(mod, "SEED_MARKERS", None)
    if not seeds:
        print("\n--- SEED_MARKERS: none defined (will fall back to MARKER_DICT union) ---")
        return 0
    in_panel = set(audit_df["gene"])
    excluded = set(audit_df.loc[audit_df["exclude_recommended"] == "yes", "gene"])
    broad = set(getattr(mod, "BROAD_TO_FINE", {}))
    print(f"\n--- SEED_MARKERS (broad-classifier seeds, >= {min_genes}/lineage) ---")
    n_fail = 0
    missing_keys = broad - set(seeds)
    extra_keys = set(seeds) - broad
    if missing_keys:
        print(f"  *** SEED_MARKERS missing broad lineages: {sorted(missing_keys)} ***")
        n_fail += 1
    if extra_keys:
        print(f"  *** SEED_MARKERS has non-BROAD_TO_FINE keys: {sorted(extra_keys)} ***")
        n_fail += 1
    for L, genes in seeds.items():
        present = [g for g in genes if g in in_panel and g not in excluded]
        absent = [g for g in genes if g not in in_panel]
        flagged = [g for g in genes if g in excluded]
        status = "ok" if len(present) >= min_genes else "FAIL"
        if status != "ok" or absent or flagged:
            print(f"  {L:<12} status={status}  present={len(present)}/{len(genes)}")
            if absent:
                print(f"      absent from panel:  {absent}")
            if flagged:
                print(f"      flagged as exclude: {flagged}")
        if status == "FAIL":
            n_fail += 1
    if n_fail == 0:
        print(f"  all {len(seeds)} broad lineages have >= {min_genes} panel-present seeds.")
    return n_fail


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", choices=("pancreas", "thymus"),
                    default="thymus")
    ap.add_argument("--audit-csv", default=None,
                    help="Override audit CSV path")
    args = ap.parse_args()

    mod = importlib.import_module(f"scripts.tissue_markers.{args.tissue}")
    if args.audit_csv:
        audit_path = Path(args.audit_csv)
        audit_df = pd.read_csv(audit_path)
    else:
        audit_path, audit_df = load_audit(args.tissue)

    print(f"=== verifying {args.tissue} markers against {audit_path.name} ===")
    print(f"audit rows: {len(audit_df):,}; "
          f"excluded: {(audit_df['exclude_recommended']=='yes').sum()}")

    n_fail = 0
    n_fail += check_panel("MARKER_DICT (scoring, ≥3)", mod.MARKER_DICT,
                          audit_df, min_genes=3, kind="scoring")
    n_fail += check_panel("ANTI_OVERRIDE_MARKERS (override)",
                          mod.ANTI_OVERRIDE_MARKERS, audit_df,
                          min_genes=1, kind="override")
    n_fail += check_panel("LINEAGE_CANARY_PANELS (canary)",
                          mod.LINEAGE_CANARY_PANELS, audit_df,
                          min_genes=1, kind="canary")
    n_fail += check_panel("PANEL_COVERAGE_PANELS (coverage)",
                          mod.PANEL_COVERAGE_PANELS, audit_df,
                          min_genes=1, kind="coverage")
    n_fail += check_seed_markers(mod, audit_df, min_genes=2)

    if n_fail == 0:
        print("\nALL panels meet the ≥3-genes-for-scoring rule.")
    else:
        print(f"\n*** {n_fail} scoring panels would be dropped — review.")


if __name__ == "__main__":
    main()
