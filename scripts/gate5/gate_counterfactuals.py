#!/usr/bin/env python
"""Regenerate the three Gate-5 tables that backed tracked conclusions but had no producing script.

`data/gate5/counterfactuals.csv`, `followup.csv` and `mes_loss.csv` back Findings 2, 4 and 8 of
`docs/GATE5_ROUND_2026-07-27.md` and the `cd11b-cd34-gate-defects` memory. They were produced ad-hoc on
2026-07-27 and a repo-wide search for their filenames returns nothing, so the numbers in those documents
could not be checked, refreshed after a gate edit, or reproduced after a re-run. This script is the
missing instrument.

It does NOT reimplement the lineage tree. It calls `phenocycler.lineage.assign_donor` with
`COMPARTMENT_GATES` temporarily swapped, which is the same technique `scripts/gold_standard/score.py`
uses, so a counterfactual arm cannot drift from what production would actually do.

    counterfactuals.csv  20 donors x 4 arms (production, no_CD11b, no_CD34, no_both): compartment
                         percentages tissue-wide and in islet core
    followup.csv         per donor: the CD11b and CD34 overlap criteria behind Findings 2 and 4 --
                         CD11b+ that are epithelial-marker+, outside islets, or carry no other Immune
                         gate marker; CD34/CD31 agreement both ways; SMA-only cells lost to higher gates
    mes_loss.csv         per donor x claiming compartment: where Vimentin+/SMA+ cells actually end up

USAGE
    ~/miniconda3/envs/scvi-env/bin/python scripts/gate5/gate_counterfactuals.py --config config.ini
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

from phenocycler import lineage                                                # noqa: E402
from phenocycler.config import COMPARTMENT_GATES, COMPARTMENT_ORDER, OTHER_LABEL, load_config  # noqa: E402

ARMS = {
    "production":  {},
    "no_CD11b":    {"Immune": [m for m in COMPARTMENT_GATES["Immune"] if m != "CD11b"]},
    "no_CD34":     {"Endothelial": [m for m in COMPARTMENT_GATES["Endothelial"] if m != "CD34"]},
    "no_both":     {"Immune": [m for m in COMPARTMENT_GATES["Immune"] if m != "CD11b"],
                    "Endothelial": [m for m in COMPARTMENT_GATES["Endothelial"] if m != "CD34"]},
}
IMMUNE = COMPARTMENT_GATES["Immune"]
MESENCHYMAL = COMPARTMENT_GATES["Mesenchymal"]

# "Epithelial-marker positive" is reported under TWO definitions, because the ad-hoc 2026-07-27 tables
# used one and `COMPARTMENT_GATES["Epithelial"]` is the other, and nothing recorded the difference.
#
#   structural  the shipped Epithelial gate MINUS the three hormones. Reproduces the published numbers
#               EXACTLY (verified on donor 6591: cd11b_epi 13,454 and cd34_not_cd31_epi 26,123 both to
#               the cell). It is also the defensible definition for these particular questions -- "is
#               this CD11b+ cell also a cytokeratin+ epithelial cell" is a different claim from "is it
#               INS+", and a hormone+ CD11b+ cell is an islet-spillover question, not an epithelial one.
#   gate        COMPARTMENT_GATES["Epithelial"] as it ships, hormones included, for continuity with
#               the lineage assignment.
#
# The two differ little for CD11b (13,494 vs 13,454, 0.3%) so the published "81.1% of CD11b+ cells are
# epithelial-marker+" stands; they differ by 5% for cd34_not_cd31_epi.
EPITHELIAL_STRUCTURAL = [m for m in COMPARTMENT_GATES["Epithelial"]
                         if m not in ("INS", "GCG", "SST")]
EPITHELIAL_GATE = COMPARTMENT_GATES["Epithelial"]


def _assign(donor: str, gated_f: str, cells_dir: Path, gate_override: dict) -> pd.DataFrame:
    """One arm, through production's own assignment function."""
    gates = dict(COMPARTMENT_GATES)
    gates.update(gate_override)
    lineage.COMPARTMENT_GATES = gates
    try:
        return lineage.assign_donor(donor, gated_f, cells_dir)
    finally:
        lineage.COMPARTMENT_GATES = dict(COMPARTMENT_GATES)


def _core_mask(out: pd.DataFrame) -> pd.Series:
    """Islet-core cells. `cell_region` is QuPath-derived and its value is `core` (not `islet_core`);
    getting this wrong silently yields an all-False mask and empty core percentages."""
    if "cell_region" not in out:
        return pd.Series(False, index=out.index)
    return out["cell_region"].astype(str) == "core"


def _pct(assigned: pd.Series) -> dict:
    n = len(assigned)
    return {c: 100.0 * float((assigned == c).sum()) / n for c in COMPARTMENT_ORDER + [OTHER_LABEL]} \
        if n else {}


def donor_tables(donor: str, gated_f: str, cells_dir: Path):
    gated = pd.read_parquet(gated_f)
    pos = {m: gated[f"{m}_pos"].to_numpy(bool) for m in
           sorted({m for g in COMPARTMENT_GATES.values() for m in g}) if f"{m}_pos" in gated.columns}

    cf_rows, base = [], None
    for arm, override in ARMS.items():
        out = _assign(donor, gated_f, cells_dir, override)
        if arm == "production":
            base = out
        core = _core_mask(out)
        row = {"donor": donor, "arm": arm, "n": len(out), **_pct(out["compartment"])}
        row.update({f"core_{k}": v for k, v in _pct(out.loc[core, "compartment"]).items()})
        cf_rows.append(row)

    # ---- followup: the CD11b and CD34 criteria, on production's assignment
    epi_any = np.logical_or.reduce([pos[m] for m in EPITHELIAL_STRUCTURAL if m in pos])
    epi_gate = np.logical_or.reduce([pos[m] for m in EPITHELIAL_GATE if m in pos])
    imm_other = np.logical_or.reduce([pos[m] for m in IMMUNE if m in pos and m != "CD11b"])
    mes_any = np.logical_or.reduce([pos[m] for m in MESENCHYMAL if m in pos])
    core = _core_mask(base).to_numpy()
    cd11b, cd31, cd34 = pos.get("CD11b"), pos.get("CD31"), pos.get("CD34")
    vim, sma = pos.get("Vimentin"), pos.get("SMA")
    fu = {
        "donor": donor, "n": len(base),
        "cd11b_n": int(cd11b.sum()),
        "cd11b_epi": int((cd11b & epi_any).sum()),
        "cd11b_epi_outside": int((cd11b & epi_any & ~core).sum()),
        "cd11b_no_other_immune": int((cd11b & ~imm_other).sum()),
        "cd11b_epi_no_other_immune": int((cd11b & epi_any & ~imm_other).sum()),
        "cd31_n": int(cd31.sum()), "cd34_n": int(cd34.sum()),
        "cd31_and_cd34": int((cd31 & cd34).sum()),
        "cd34_not_cd31": int((cd34 & ~cd31).sum()),
        "cd34_not_cd31_vim": int((cd34 & ~cd31 & vim).sum()),
        "cd34_not_cd31_epi": int((cd34 & ~cd31 & epi_any).sum()),
        "cd34_not_cd31_core": int((cd34 & ~cd31 & core).sum()),
        "mes_elig": int(mes_any.sum()),
        "sma_only": int((sma & ~vim).sum()),
        "sma_only_lost": int((sma & ~vim & (base["compartment"] != "Mesenchymal").to_numpy()).sum()),
        # the same three quantities under the SHIPPED Epithelial gate, hormones included
        "cd11b_epi_gate": int((cd11b & epi_gate).sum()),
        "cd11b_epi_outside_gate": int((cd11b & epi_gate & ~core).sum()),
        "cd34_not_cd31_epi_gate": int((cd34 & ~cd31 & epi_gate).sum()),
    }

    # ---- mes_loss: where the Vimentin+/SMA+ cells actually go
    comp = base["compartment"].to_numpy()
    ml = [{"donor": donor, "claimed_by": c, "n": int((mes_any & (comp == c)).sum()),
           "mes_n": int(mes_any.sum())} for c in COMPARTMENT_ORDER + [OTHER_LABEL]]
    # the CD34-only route into Endothelial, called out separately because it is the actionable lever
    if cd34 is not None and cd31 is not None:
        pod = pos.get("Podoplanin", np.zeros(len(base), bool))
        ml.append({"donor": donor, "claimed_by": "Endothelial_via_CD34_only",
                   "n": int((mes_any & (comp == "Endothelial") & cd34 & ~cd31 & ~pod).sum()),
                   "mes_n": int(mes_any.sum())})
    print(f"[cf] {donor}: {len(base):,} cells, {len(ARMS)} arms", flush=True)
    return cf_rows, fu, ml


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--config", type=Path, default=REPO / "config.ini")
    p.add_argument("--outdir", type=Path, default=REPO / "data" / "gate5")
    p.add_argument("--donors", nargs="*", default=None)
    a = p.parse_args(argv)

    cfg = load_config(a.config)
    files = sorted(glob.glob(str(cfg.restore_gated_dir / "donor_id=*" / "*.parquet")))
    if not files:
        raise SystemExit(f"[err] no gated parquets under {cfg.restore_gated_dir}")
    if a.donors:
        files = [f for f in files if Path(f).parent.name.split("=", 1)[1] in set(a.donors)]

    cf, fu, ml = [], [], []
    for f in files:
        donor = Path(f).parent.name.split("=", 1)[1]
        c, u, m = donor_tables(donor, f, cfg.cells_dir)
        cf += c
        fu.append(u)
        ml += m

    a.outdir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(cf).to_csv(a.outdir / "counterfactuals.csv", index=False)
    pd.DataFrame(fu).to_csv(a.outdir / "followup.csv", index=False)
    pd.DataFrame(ml).to_csv(a.outdir / "mes_loss.csv", index=False)

    d = pd.DataFrame(cf)
    print("\n=== cohort mean compartment % by arm (the Finding 4 table) ===")
    cols = [c for c in COMPARTMENT_ORDER + [OTHER_LABEL] if c in d.columns]
    cols += [c for c in ("core_Endothelial", "core_Immune", "core_Epithelial") if c in d.columns]
    print(d.groupby("arm")[cols].mean().round(1).to_string())
    print(f"\n-> {a.outdir}/counterfactuals.csv, followup.csv, mes_loss.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
