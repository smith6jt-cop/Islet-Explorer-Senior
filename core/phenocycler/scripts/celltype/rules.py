"""Turn `pancreas_phenotype_rules_v4.csv` into reference vectors pSAM's scoring can use.

pSAM scores a cell by the cosine of its marker vector against a set of reference
pseudospectra -- one column per cell type, one row per marker, values 0-1.  The maintainer's
rule table is the same object in a different encoding, with two structural differences that
have to be resolved before a cosine means anything:

**It is a hierarchy, not a set of references.**  ``Proliferating_beta``'s row says only
``Ki67 pos``, because its parent supplies "is a beta cell".  Read as an independent
reference that is indistinguishable from every other proliferating phenotype -- measured,
ten of them collapse to one.  :func:`flatten` composes each phenotype with its ancestors, so
``Proliferating_beta`` becomes ``INS allpos, SST neg, GCG neg, CD68 neg, CD20 neg, CD3e neg,
Ki67 pos``.  That takes the median reference from 1 constraint to 8 and the number of
distinct references from 66 to 100 of 104.

**Cosine cannot express a negative constraint.**  A 0 in a reference means "does not
contribute", not "must be absent", and 84 of the table's 208 constraints are negative.
Measured after flattening, dropping them costs 4 further collapses -- 100 distinct references
becomes 96 -- and every phenotype lost is either a vacuous child that restates its parent
(``Ductal_cell``, ``Neuron``, ``Stellate_cell``) or one defined purely by negation
(``DN_T_cell``, ``Fibroblast``, ``Quiescent_stellate``, ``Pericyte``, ``Stressed_beta``).
So both encodings are provided and neither is assumed:

* ``zero``   -- pSAM's own semantics: positive constraints 1, everything else 0.
* ``signed`` -- a small extension: negative constraints -1, so a cell bright in a marker the
  reference forbids is pushed away from it rather than merely not pulled toward it.  Still
  magnitude-invariant, which is the property this lane adopted cosine for.
"""
from __future__ import annotations

import json
import pathlib
import re
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel                                     # noqa: E402

RULES_CSV = panel.REPO / "pancreas_phenotype_rules_v4.csv"
CROSSWALK = HERE / "marker_name_crosswalk.json"

POSITIVE = ("pos", "allpos", "anypos")
NEGATIVE = ("neg", "allneg")


def load_crosswalk() -> dict:
    """Display name -> registry name, from the frozen file.

    A frozen map rather than a regex: a regex that silently misses one name turns that rule's
    constraint into no constraint, which fails open.  :func:`build_crosswalk` regenerates it
    and the test suite asserts the two still agree, so a v5 rename fails a test instead of
    dropping a constraint.
    """
    return json.loads(CROSSWALK.read_text())["map"]


def build_crosswalk() -> dict:
    cols = list(pd.read_csv(RULES_CSV, nrows=0).columns)[2:]
    known = set(panel.markers(include_variant=True)) | set(panel.REGISTRY_EXCLUDED)
    out, missing = {}, []
    for c in cols:
        t = re.sub(r"[^A-Za-z0-9]", "_", c)
        (out if t in known else missing).__setitem__(c, t) if t in known else missing.append(c)
    if missing:
        raise ValueError(f"no registry marker for: {missing}")
    return out


def load_rules() -> pd.DataFrame:
    d = pd.read_csv(RULES_CSV)
    d.columns = ["parent", "phenotype"] + [load_crosswalk()[c] for c in d.columns[2:]]
    return d


def flatten(rules: pd.DataFrame) -> dict:
    """phenotype -> {marker: token}, composed with every ancestor.  The child wins."""
    mk = [c for c in rules.columns if c not in ("parent", "phenotype")]
    row = {r["phenotype"]: {c: r[c] for c in mk if pd.notna(r[c])}
           for _, r in rules.iterrows()}
    par = dict(zip(rules["phenotype"], rules["parent"]))
    out = {}
    for ph in rules["phenotype"]:
        chain, cur = [], ph
        while cur and cur != "all":
            chain.append(cur)
            cur = par.get(cur)
            if cur in chain:          # a cycle would loop forever; refuse rather than guess
                raise ValueError(f"cycle in the rule hierarchy at {ph}")
        ref = {}
        for a in reversed(chain):
            ref.update(row.get(a, {}))
        out[ph] = ref
    return out


def reference_matrix(flat: dict, markers: list, encoding: str = "zero",
                     strict: bool = True) -> tuple:
    """markers x phenotypes in the encoding named, plus what could not be encoded.

    A constraint on a marker outside ``markers`` cannot be represented, and dropping one
    silently is a fail-open: the rule simply stops constraining, and two phenotypes that
    differ only by that marker become indistinguishable with nothing to show it.  That is not
    hypothetical -- CD38 is absent from 7 donors, and its loss is exactly what makes
    ``B_cell`` and ``CD38_high_B`` identical.  So every dropped constraint is returned, and
    with ``strict`` the caller has to look at them.
    """
    if encoding not in ("zero", "signed"):
        raise ValueError(f"unknown encoding {encoding!r}")
    M = pd.DataFrame(0.0, index=markers, columns=list(flat))
    dropped: dict = {}
    for ph, ref in flat.items():
        for m, tok in ref.items():
            if m not in M.index:
                dropped.setdefault(m, []).append(ph)
                continue
            if tok in POSITIVE:
                M.loc[m, ph] = 1.0
            elif tok in NEGATIVE and encoding == "signed":
                M.loc[m, ph] = -1.0
    if strict and dropped:
        raise ValueError(
            "constraints on markers outside the panel would be dropped silently: "
            + "; ".join(f"{m} ({len(v)} phenotypes)" for m, v in sorted(dropped.items()))
            + " -- pass strict=False and record them, or widen the marker set")
    return M, dropped


def distinct_references(M: pd.DataFrame) -> dict:
    """Which phenotypes share an identical reference vector, and are therefore unresolvable."""
    groups: dict = {}
    for ph in M.columns:
        groups.setdefault(tuple(M[ph].to_numpy()), []).append(ph)
    return {k: v for k, v in groups.items() if len(v) > 1}


def main() -> int:
    cw = build_crosswalk()
    CROSSWALK.write_text(json.dumps(
        {"source": RULES_CSV.name, "n": len(cw),
         "note": "Frozen so a v5 rename fails a test rather than silently dropping a "
                 "constraint. Regenerate with `python scripts/celltype/rules.py`.",
         "map": dict(sorted(cw.items()))}, indent=2))
    print(f"wrote {CROSSWALK}  ({len(cw)} marker columns)")

    rules = load_rules()
    flat = flatten(rules)
    markers = panel.markers()
    print(f"{len(flat)} phenotypes, {len(markers)} markers in this lane's panel")
    for enc in ("zero", "signed"):
        M, dropped = reference_matrix(flat, markers, enc, strict=False)
        if dropped:
            print(f"    NOT ENCODED (marker outside the panel): "
                  + "; ".join(f"{m} -> {len(v)} phenotypes" for m, v in sorted(dropped.items())))
        dup = distinct_references(M)
        n_dup = sum(len(v) for v in dup.values())
        print(f"  encoding {enc:<7} distinct references {len(M.columns) - n_dup + len(dup):3d}"
              f" of {len(M.columns)}   ({n_dup} phenotypes in {len(dup)} unresolvable groups)")
        for v in sorted(dup.values(), key=len, reverse=True)[:6]:
            print(f"      {' == '.join(v)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
