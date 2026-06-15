#!/usr/bin/env python3
"""
Convert the scimap hierarchical rules CSV into an Astir marker YAML.

Astir (Geuijen/Campbell, Cell Systems 2021) is a reference-free probabilistic classifier:
given a cell_type -> POSITIVE-markers map, it learns each type's marker distribution and
assigns per-cell posterior probabilities — NO per-marker gating. That is exactly why it is
the right tool here: empirically the per-marker GMM gates on this 59-plex panel are unreliable
in BOTH directions (SOX2/EpCAM/Ker8_18 over-called ~90% positive; PanCK under-called ~1-5%),
so neither univariate gating nor relative-z cluster annotation is robust. Astir learns marker
informativeness instead of thresholding.

Mapping (terminal cell types only — Astir assigns leaves):
  leaf positive markers = the leaf's OWN positive markers (the discriminating ones), plus the
  positive markers of any INTERMEDIATE parent (e.g. 'T cell' -> CD3e, 'Macrophage' ->
  CD68/CD163) so subtypes carry their lineage marker. Top lineages (Endocrine/Epithelial/...)
  are NOT merged, so Beta stays {INS} not {INS,GCG,SST}.

Usage:
  python scripts/senior/rules_to_astir_yaml.py --out scripts/senior/astir_markers.yaml
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd
import yaml

RULES = Path(__file__).resolve().parent / "phenotype_rules_senior.csv"
POS_KW = {"pos", "anypos", "allpos"}
# Intermediate parents whose own positive markers should be inherited by their subtypes.
INTERMEDIATE = {"T cell", "Macrophage"}


def parse(rules_path):
    df = pd.read_csv(rules_path)
    markers = list(df.columns[2:])
    own_pos, parent = {}, {}
    is_parent = set()
    for _, r in df.iterrows():
        name = str(r["phenotype_name"]).strip()
        parent[name] = str(r["parent_class"]).strip()
        is_parent.add(parent[name])
        pos = [m for m in markers if not pd.isna(r[m]) and str(r[m]).strip().lower() in POS_KW]
        own_pos[name] = pos
    leaves = [n for n in own_pos if n not in is_parent]
    return own_pos, parent, leaves


def leaf_markers(leaf, own_pos, parent):
    mk = list(own_pos.get(leaf, []))
    p = parent.get(leaf)
    # inherit only from intermediate parents (carry the lineage marker, e.g. CD3e for T subtypes)
    while p and p != "all":
        if p in INTERMEDIATE:
            mk += [m for m in own_pos.get(p, []) if m not in mk]
        p = parent.get(p)
    return mk


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--rules", type=Path, default=RULES)
    ap.add_argument("--exclude-markers", nargs="*", default=None,
                    help="markers to drop from all cell-type definitions (e.g. audit DROP_background like SOX2)")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args(argv)
    drop = set(a.exclude_markers or [])
    own_pos, parent, leaves = parse(a.rules)
    cell_types = {}
    for lf in leaves:
        mk = [m for m in leaf_markers(lf, own_pos, parent) if m not in drop]
        if mk:
            cell_types[lf] = mk
        else:
            print(f"[warn] leaf '{lf}' has no positive markers — skipped")
    doc = {"cell_types": cell_types, "cell_states": {}}
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with open(a.out, "w") as f:
        yaml.safe_dump(doc, f, sort_keys=False, default_flow_style=False)
    print(f"[done] wrote {a.out} with {len(cell_types)} cell types:")
    for k, v in cell_types.items():
        print(f"    {k:24s} <- {', '.join(v)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
