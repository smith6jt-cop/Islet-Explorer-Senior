#!/usr/bin/env python3
"""Phase 2 -- name the immune cells, inside the class they are already in.

1,479,904 cells carry the single label ``Immune`` and nothing finer, which is the population the
study question is about.  This resolves them against the established rule table
(``phenotype_workflow.csv``, the phenotyping guide) and **re-decides nothing outside them**: a cell
that is not ``Immune`` is not read, not scored and not written.  ``cell_class`` is never modified.

The hierarchy is the table's own, three levels deep:

    Immune -> Imm-T-cell | Imm-B-cell | Imm-macrophage | Imm-dendritic-cell
              Imm-NK-cell | Imm-neutrophil | Imm-plasma-cell
      Imm-T-cell     -> T-regulatory | T-exhausted | T-progenitor-memory
                        T-CD8-cytotoxic | T-CD8 | T-CD4-helper
      Imm-macrophage -> Mac-M1-like | Mac-M2-like

**CD4 and CD8 sit where the literature puts them, because the table already does.** CD8 appears
only beneath ``Imm-T-cell`` -- never as a top-level discriminator -- and CD4 is carried in the
panel but named by no rule at all: ``T-CD4-helper`` is ``CD8 = allneg, FOXP3 = allneg``. So CD8
positivity is decisive for the T lineage and CD4 alone names nothing, which is the instruction,
and it needed no new rule to honour.

**Panel variants are availability, not absence.** ``Imm-plasma-cell`` needs CD38, which seven
donors do not carry; ``Imm-NK-cell`` names CD57, which is at background across batch 1 (scaled p99
0.0035-0.0266 against a 0.10 level). A rule whose markers a donor lacks is ``unavailable`` for that
donor and recorded as such -- never scored as a negative, which would quietly call every cell
"not a plasma cell" on exactly the donors that cannot be asked.

**Multiple rules can fire.** That is not resolved by a silent precedence: every match is counted,
the most specific rule wins (most required markers, ties broken by name for determinism), and both
the match count and the losing rules are written per cell so the choice is auditable.

    python core/phenocycler/scripts/typing/immune_subtype.py [--donors 6539]
"""
from __future__ import annotations

import argparse
import json
import sys
import time

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = __import__("pathlib").Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "celltype"))
import cohort as COHORT
import rules as RULES
import subpanel_io as S

OUT = S.OUT / "immune"

#: The same level the coarse immune definitions were evaluated at (RECONSTRUCTION.md).
REQUIRE_MIN = 0.10

LEVEL2_PARENT = "Immune"
LEVEL3_PARENTS = ("Imm-T-cell", "Imm-macrophage")


def immune_rules() -> tuple[pd.DataFrame, dict]:
    table = RULES.load_rules()
    reqs = RULES.requirements(table, with_negatives=True)
    return table, reqs


def children_of(table: pd.DataFrame, parent: str) -> list[str]:
    return [p for p in table.loc[table["parent"] == parent, "phenotype"].astype(str).unique()]


def evaluate(scaled: pd.DataFrame, requirement, available: set[str]):
    """(fires, unavailable) for one rule over every row.

    ``requirement`` is (all_of, [any_of groups], none_of) as the table declares it.
    A rule naming a marker this donor does not carry is unavailable, not false.
    """
    all_of, any_groups, none_of = requirement
    needed = set(all_of) | set(none_of) | {m for g in any_groups for m in g}
    if not needed <= available:
        return None, sorted(needed - available)

    fires = np.ones(len(scaled), dtype=bool)
    for marker in all_of:
        fires &= scaled[marker].to_numpy() >= REQUIRE_MIN
    for group in any_groups:
        if not group:
            continue
        hit = np.zeros(len(scaled), dtype=bool)
        for marker in group:
            hit |= scaled[marker].to_numpy() >= REQUIRE_MIN
        fires &= hit
    for marker in none_of:
        fires &= scaled[marker].to_numpy() < REQUIRE_MIN
    return fires, []


def resolve(frame: pd.DataFrame, candidates: list[str], reqs, scaled, available, specificity):
    """Assign the most specific firing rule; record the count and the runners-up."""
    hits, skipped = {}, {}
    for name in candidates:
        fires, missing = evaluate(scaled, reqs[name], available)
        if fires is None:
            skipped[name] = missing
            continue
        hits[name] = fires
    if not hits:
        return (pd.Series([None] * len(frame), index=frame.index),
                pd.Series(0, index=frame.index), pd.Series([""] * len(frame), index=frame.index),
                skipped)

    order = sorted(hits, key=lambda n: (-specificity[n], n))
    stacked = np.vstack([hits[n] for n in order])
    n_match = stacked.sum(axis=0)
    first = np.argmax(stacked, axis=0)
    call = np.where(n_match > 0, np.array(order)[first], None)
    others = [
        ";".join(n for k, n in enumerate(order) if stacked[k, i] and k != first[i])
        if n_match[i] > 1 else ""
        for i in range(len(frame))
    ]
    return (pd.Series(call, index=frame.index), pd.Series(n_match, index=frame.index),
            pd.Series(others, index=frame.index), skipped)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donors", nargs="*", default=None)
    parser.add_argument("--include-marginal", action="store_true")
    parser.add_argument("--all", action="store_true", help="every manifest donor, not just clean")
    args = parser.parse_args()

    if args.all:
        manifest = json.loads((S.REPO / "core/phenocycler/data/qupath_manifest.json").read_text())
        sys.path.insert(0, str(S.REPO / "core/phenocycler"))
        from phenocycler import cohort as pc
        donors = sorted({r["donor_id"] for r in manifest["images"]} - set(pc.DONOR_EXCLUSIONS))
    else:
        donors = COHORT.resolve(args.donors, include_marginal=args.include_marginal)

    table, reqs = immune_rules()
    level2 = children_of(table, LEVEL2_PARENT)
    level3 = {p: children_of(table, p) for p in LEVEL3_PARENTS}
    specificity = {n: len(reqs[n][0]) + sum(len(g) for g in reqs[n][1]) + len(reqs[n][2])
                   for n in reqs}
    print(f"level 2 ({len(level2)}): {', '.join(level2)}")
    for parent, kids in level3.items():
        print(f"  {parent} ({len(kids)}): {', '.join(kids)}")
    print()

    scale = S.load_scale()
    OUT.mkdir(parents=True, exist_ok=True)
    summary, skipped_report = [], {}

    for donor in donors:
        started = time.time()
        cells = pq.read_table(S.OUT / f"cells/donor_id={donor}/data.parquet",
                              columns=["object_id", "cell_class"]).to_pandas()
        immune = cells[cells["cell_class"] == "Immune"]
        if immune.empty:
            print(f"  {donor}: no Immune cells"); continue

        available_all = set(pq.ParquetFile(
            S.ALIGNED / f"donor_id={donor}/data.parquet").schema_arrow.names)
        markers = sorted((available_all & set(scale.index)) - {"DAPI"})
        values = pq.read_table(S.ALIGNED / f"donor_id={donor}/data.parquet",
                               columns=["object_id"] + markers).to_pandas()
        values = immune.merge(values, on="object_id", how="left", validate="one_to_one")
        scaled = pd.DataFrame(
            S.apply_scale(values[markers].to_numpy(), markers, scale),
            columns=markers, index=values.index)

        call2, n2, other2, skip2 = resolve(values, level2, reqs, scaled, set(markers), specificity)
        out = pd.DataFrame({
            "donor_id": donor, "object_id": values["object_id"],
            "immune_class": call2, "immune_class_matches": n2, "immune_class_also": other2,
            "immune_subtype": None, "immune_subtype_matches": 0, "immune_subtype_also": "",
        })
        for parent, kids in level3.items():
            mask = (out["immune_class"] == parent).to_numpy()
            if not mask.any():
                continue
            sub, n3, other3, skip3 = resolve(values[mask], kids, reqs,
                                             scaled[mask], set(markers), specificity)
            out.loc[mask, "immune_subtype"] = sub.to_numpy()
            out.loc[mask, "immune_subtype_matches"] = n3.to_numpy()
            out.loc[mask, "immune_subtype_also"] = other3.to_numpy()
            skip2.update(skip3)

        skipped_report[donor] = skip2
        target = OUT / f"donor_id={donor}"
        target.mkdir(parents=True, exist_ok=True)
        out.to_parquet(target / "data.parquet", index=False)

        named = out["immune_class"].notna().sum()
        multi = int((out["immune_class_matches"] > 1).sum())
        print(f"  {donor}: {len(out):>9,} Immune | named {100*named/len(out):5.1f}% | "
              f"multi-match {100*multi/len(out):4.1f}% | "
              f"unavailable rules: {sorted(skip2) or 'none'} | {time.time()-started:.0f}s",
              flush=True)
        summary.append(out.assign(donor=donor))

    if not summary:
        return 1
    everything = pd.concat(summary, ignore_index=True)
    print(f"\n{len(everything):,} Immune cells over {len(donors)} donors")
    print("\nimmune_class:")
    print(everything["immune_class"].fillna("(no rule fired)").value_counts().to_string())
    sub = everything[everything["immune_subtype"].notna()]
    if len(sub):
        print("\nimmune_subtype:")
        print(sub["immune_subtype"].value_counts().to_string())
    (OUT / "unavailable_rules.json").write_text(json.dumps(skipped_report, indent=2))
    print(f"\n-> {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
