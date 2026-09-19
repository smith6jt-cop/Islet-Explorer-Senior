#!/usr/bin/env python3
"""``phenotype_rules.csv`` -- what the app's rules modal and every inline hint read.

The app's schema is a marker gate: per phenotype, which markers must be positive, which must be
negative, and which groups need at least one.  This classification is not a marker gate.  It asks a
different question per context and compares against a **reference level** measured from the ordinary
tissue, and two of its eight calls are not expressible as marker requirements at all -- ``Immune`` is
a union of six named combinations (the permissive flat OR was measured to remove 9,310,209 cells
carrying *more* keratin than the exocrine bulk), and ``Exocrine`` is the complement.

So each row carries both: the markers that answer that phenotype's question, in the columns, and the
rule itself in a ``description`` column.  ``format_phenotype_rule()`` prefers the description when
one is present, so a hint states the method rather than implying a flat OR, while the modal's marker
columns still say which channels are involved.

Every number here comes from ``RECONSTRUCTION.md``, which was verified against the calls on disk.

WHAT THIS FILE DESCRIBES IS THE TABLE ON DISK, NOT THE 2026-09-18 ISLET CHAIN.  The descriptions
below state the hormone level 0.1010 and say CD99 "gates nothing", and both are true of
``typing/cells`` as written on 2026-09-15 -- which is the build the app serves, because no donor
has been re-run.  When the chain does run, the hormone sentence becomes wrong on two counts: the
levels are per donor rather than pooled, and CD99 at its own donor x stratum level is what decides
whether a cell is endocrine.  Rebuild this file in the same pass as the typing.

    python core/phenocycler/scripts/typing/build_phenotype_rules.py
"""
from __future__ import annotations

import csv
import sys

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent))
import subpanel_io as S

APP = S.OUT / "app"

#: Column order: the hormones, then each structure's own markers, then the immune channels the
#: named combinations are built from.  CD31 is deliberately absent -- the vessel objects were drawn
#: by `Cavdil1200` / `CD34dil1430`, so CD31 was never what made them.
MARKERS = [
    "INS", "GCG", "SST",
    "Caveolin", "CD34", "B3TUBB", "CD56", "Podoplanin",
    "CD3e", "CD68", "Iba1", "CD20", "CD79a", "MPO", "CD66", "CD11b",
    "CD11c", "CD209", "IDO1", "HLA_DR", "Granzyme_B",
    "Ker8_18", "E_cadherin", "Pan_Cytokeratin",
]

_HORMONE = ("Inside an islet, where no structure object claimed the cell: max(INS, GCG, SST) on the "
            "scaled values clears the hormone reference level 0.1010, and {marker} wins the argmax "
            "of the three. CD99 is in the islet subpanel -- it is one of the markers the islet "
            "outline's pixel classifier was built on -- but it is not in the argmax and gates "
            "nothing here: three variants were tested against the calls and only this one fits.")

_STRUCTURE = ("A structure is a region and a cell type: the centroid is inside a {classifier} "
              "classified object AND that object's own markers clear its reference level {level}. A "
              "cell the object merely crosses keeps its own identity. This test is asked inside "
              "islets too, with the same object evidence.")

#: (parent, phenotype, {marker: pos|anypos|allneg}, description)
RULES = [
    ("all", "Beta_cell",  {"INS": "pos"},  _HORMONE.format(marker="INS")),
    ("all", "Alpha_cell", {"GCG": "pos"},  _HORMONE.format(marker="GCG")),
    ("all", "Delta_cell", {"SST": "pos"},  _HORMONE.format(marker="SST")),
    ("all", "Blood_Vessel",
     {"Caveolin": "anypos", "CD34": "anypos"},
     _STRUCTURE.format(classifier="Cavdil1200 / CD34dil1430", level="0.3073")
     + " CD31 is not used: it is not what the objects were drawn from."),
    ("all", "Neural",
     {"B3TUBB": "anypos", "CD56": "anypos"},
     _STRUCTURE.format(classifier="B3Tub_Ero0.5_Thr600", level="0.2272")),
    ("all", "Lymphatic",
     {"Podoplanin": "pos"},
     _STRUCTURE.format(classifier="PDPNmorphClos0-5_Thr300", level="0.1776")),
    ("all", "Immune",
     {m: "anypos" for m in ("CD3e", "CD68", "Iba1", "CD20", "CD79a", "MPO", "CD66", "CD11b",
                            "CD11c", "CD209", "IDO1", "HLA_DR", "CD56", "Granzyme_B")},
     "Not a flat OR over these markers: the union of six named combinations, each with its declared "
     "negatives enforced, evaluated at require_min = 0.10 on the scaled values -- CD3e; CD68 or "
     "Iba1; CD20 and CD79a; CD66 and (CD11b or MPO); CD11c / CD209 / IDO1 / HLA-DR; CD56 with "
     "Granzyme B or CD57. A supporting tier at 0.8715 settles cells the definitions leave disputed. "
     "On this ten-donor cohort CD57 sits at background (scaled p99 0.0035-0.0266 against a 0.10 "
     "level), so the NK definition reduces to CD56 with Granzyme B. Immune subtypes are not yet "
     "named -- this is the resolved coarse class."),
    ("all", "Exocrine",
     {"Ker8_18": "anypos", "E_cadherin": "anypos", "Pan_Cytokeratin": "anypos"},
     "The complement, not a gate: a cell in no classified object, outside every islet, that no "
     "structure and no immune definition claimed. Ker8-18, E-cadherin and Pan-Cytokeratin separate "
     "this bulk from every removed class by 5-15x, which is why they are listed, but the call is "
     "the fall-through. Never read as acinar or ductal: the panel carries no amylase, trypsin, "
     "CPA1, elastase, CK19, CA19-9, SOX9, HNF1B, CFTR or MUC1, so exocrine is one label."),
]


def main() -> int:
    target = APP / "phenotype_rules.csv"
    with open(target, "w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["", ""] + MARKERS + ["description"])
        for parent, phenotype, calls, description in RULES:
            unknown = set(calls) - set(MARKERS)
            if unknown:
                raise ValueError(f"{phenotype}: {sorted(unknown)} not in the marker columns")
            writer.writerow([parent, phenotype]
                            + [calls.get(m, "") for m in MARKERS] + [description])

    print(f"{len(RULES)} phenotypes x {len(MARKERS)} markers -> {target}")
    for _, phenotype, calls, _ in RULES:
        print(f"  {phenotype:<14s} {', '.join(f'{m}:{v}' for m, v in calls.items())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
