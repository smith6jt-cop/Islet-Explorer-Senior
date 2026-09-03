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
import os
import pathlib
import re
import sys

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panel                                     # noqa: E402

#: The maintainer's rule table.  Overridable so a new version can be checked before it is
#: adopted; the crosswalk is regenerated from whichever file is named, and the test suite
#: asserts the frozen map still matches, so a rename in a v5 fails a test rather than
#: silently turning that rule's constraint into no constraint.
#: Candidates in preference order; the first that exists wins.  Override with
#: ``CELLTYPE_RULES_CSV`` to check a new table before adopting it.
RULES_CANDIDATES = ("phenotype_workflow.csv", "phenotype_workflow_59plex_v2.csv",
                    "pancreas_phenotype_rules_v4.csv")


def _default_rules_csv() -> pathlib.Path:
    for name in RULES_CANDIDATES:
        if (panel.REPO / name).exists():
            return panel.REPO / name
    return panel.REPO / RULES_CANDIDATES[0]


RULES_CSV = pathlib.Path(os.environ["CELLTYPE_RULES_CSV"]) \
    if os.environ.get("CELLTYPE_RULES_CSV") else _default_rules_csv()
CROSSWALK = HERE / "marker_name_crosswalk.json"

#: ``anypos`` is an OR-group: at least one member.  A row may need TWO of them -- T-exhausted
#: is "(CD4 or CD8) AND (VISTA or LAG3 or TOX)" -- and a bare ``anypos`` on both cannot say
#: which marker belongs to which, so they merge into one OR that CD68 alone could satisfy.
#: ``anypos2``, ``anypos3`` ... open further groups on the same row.  Backward compatible:
#: a table using only ``anypos`` behaves exactly as before.
OR_POS = re.compile(r"^anypos(\d*)$")
OR_NEG = re.compile(r"^anyneg(\d*)$")

POSITIVE = ("pos", "allpos", "anypos")
NEGATIVE = ("neg", "allneg", "anyneg")


def is_positive(token: str) -> bool:
    return token in ("pos", "allpos") or bool(OR_POS.match(str(token)))


def is_negative(token: str) -> bool:
    return token in ("neg", "allneg") or bool(OR_NEG.match(str(token)))

#: How a table names a phenotype's parent.  ``pancreas_phenotype_rules_v4.csv`` used an
#: explicit ``parent`` column; ``phenotype_workflow_59plex_v2.csv`` puts the hierarchy in the
#: name -- "Immune - T cell - CD4+ helper" under "Immune - T cell" under "Immune".  Both are
#: read, because either way a row is a partial statement: in the 59-plex table a child drops
#: between 2 and 8 of its parent's constraints and overrides up to 4, so reading a row alone
#: loses real evidence.
NAME_SEPARATOR = " - "


def _marker_columns(path: pathlib.Path | None = None) -> list:
    cols = list(pd.read_csv(path or RULES_CSV, nrows=0).columns)
    return [c for c in cols if not str(c).startswith("Unnamed:")]


def load_crosswalk(path: pathlib.Path | None = None) -> dict:
    """Display name -> registry name, from the frozen file.

    A frozen map rather than a regex: a regex that silently misses one name turns that rule's
    constraint into no constraint, which fails open.  :func:`build_crosswalk` regenerates it
    and the test suite asserts the two still agree, so a v5 rename fails a test instead of
    dropping a constraint.
    """
    cw = json.loads(CROSSWALK.read_text())["map"]
    missing = [c for c in _marker_columns(path) if c not in cw]
    if missing:
        raise ValueError(
            f"{len(missing)} marker column(s) are not in the frozen crosswalk: {missing}. "
            "Regenerate it with `python scripts/celltype/rules.py` and review the diff -- an "
            "unmapped column would silently turn its rule's constraint into no constraint.")
    return cw


def build_crosswalk(path: pathlib.Path | None = None) -> dict:
    cols = _marker_columns(path)
    known = set(panel.markers(include_variant=True)) | set(panel.REGISTRY_EXCLUDED)
    out, missing = {}, []
    for c in cols:
        t = re.sub(r"[^A-Za-z0-9]", "_", c)
        (out if t in known else missing).__setitem__(c, t) if t in known else missing.append(c)
    if missing:
        raise ValueError(f"no registry marker for: {missing}")
    return out


def load_rules(path: pathlib.Path | None = None) -> pd.DataFrame:
    """Any supported rule table, as ``parent, phenotype, <marker columns>``.

    Two layouts are read.  A table whose first two columns are unnamed is the v4 layout,
    ``parent, phenotype``.  A table with one unnamed column is the 59-plex layout, where the
    hierarchy lives in the name and the parent is the longest proper prefix that is itself a
    row -- so "Immune - T cell - CD4+ helper" sits under "Immune - T cell", and a leaf whose
    prefix is not a row (``Endothelial - vascular``) is top level.
    """
    d = pd.read_csv(path or RULES_CSV)
    unnamed = [c for c in d.columns if str(c).startswith("Unnamed:")]
    cw = load_crosswalk(path)
    if len(unnamed) >= 2:
        d.columns = ["parent", "phenotype"] + [cw[c] for c in d.columns[2:]]
        return d
    d = d.rename(columns={unnamed[0]: "phenotype"})
    d.columns = ["phenotype"] + [cw[c] for c in d.columns[1:]]
    names = set(d["phenotype"])

    def parent(n: str) -> str:
        parts = n.split(NAME_SEPARATOR)
        for k in range(len(parts) - 1, 0, -1):
            cand = NAME_SEPARATOR.join(parts[:k])
            if cand in names:
                return cand
        return "all"

    d.insert(0, "parent", [parent(n) for n in d["phenotype"]])
    return d


def flatten(rules: pd.DataFrame, with_source: bool = False) -> dict:
    """phenotype -> {marker: token}, composed with every ancestor.  The child wins.

    With ``with_source`` the value is ``(token, origin_phenotype)``.  That provenance is not
    decoration: an ``anypos`` group means "at least one of THESE", and a flattened reference
    can carry several such groups from different levels -- "is an immune cell" from Immune and
    "is CD4+ or CD8+" from a T subset.  Merged into one OR they would be satisfied by CD68
    alone, which is the opposite of what the subset rule asks.
    """
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
            ref.update({m: (t, a) for m, t in row.get(a, {}).items()})
        out[ph] = ref if with_source else {m: t for m, (t, _) in ref.items()}
    return out


def reference_matrix(flat: dict, markers: list, encoding: str = "zero",
                     strict: bool = True, split_or_groups: bool = True) -> tuple:
    """markers x phenotypes in the encoding named, plus what could not be encoded.

    Build this over the panel UNION (``panel.markers(include_variant=True)``), not the
    cohort intersection.  The cohort is three panel variants, so marker availability is a
    per-donor property: a rule naming CD38 is evaluable on the 13 donors that carry it and
    ``unavailable`` on the other 7.  Restricting the reference to the intersection instead
    throws that rule away for everyone, and does it silently -- which is exactly what made
    ``B_cell`` and ``CD38_high_B`` come out identical.  Use :func:`evaluable_donors` to get
    the n a phenotype is actually measured on.

    Anything still unrepresentable is returned rather than skipped, and ``strict`` makes the
    caller look at it.

    ``split_or_groups`` shares unit weight across an ``anypos`` / ``anyneg`` group rather than
    giving every member 1.  A group says "expect elevation in ONE of these", and encoding it
    as k separate 1s says "expect all k", which also makes a parent's OR-group identical to a
    child that names one of its members: on the 59-plex table that alone left Immune,
    Imm-T-cell, Imm-B-cell and Imm-neutrophil sharing one vector.  Sharing the weight keeps
    the group's contribution to the norm the same while letting a child that commits to one
    member out-align it.
    """
    if encoding not in ("zero", "signed"):
        raise ValueError(f"unknown encoding {encoding!r}")
    M = pd.DataFrame(0.0, index=markers, columns=list(flat))
    dropped: dict = {}
    for ph, ref in flat.items():
        usable = {m: t for m, t in ref.items() if m in M.index}
        for m in set(ref) - set(usable):
            dropped.setdefault(m, []).append(ph)
        n_any: dict = {}
        for t in usable.values():
            if OR_POS.match(str(t)) or OR_NEG.match(str(t)):
                n_any[t] = n_any.get(t, 0) + 1
        for m, tok in usable.items():
            w = 1.0
            if split_or_groups and n_any.get(tok, 0) > 1:
                w = 1.0 / n_any[tok]
            if is_positive(tok):
                M.loc[m, ph] = w
            elif is_negative(tok) and encoding == "signed":
                M.loc[m, ph] = -w
    if strict and dropped:
        raise ValueError(
            "constraints on markers outside the panel would be dropped silently: "
            + "; ".join(f"{m} ({len(v)} phenotypes)" for m, v in sorted(dropped.items()))
            + " -- pass strict=False and record them, or widen the marker set")
    return M, dropped


def evaluable_donors(flat: dict) -> pd.DataFrame:
    """Per phenotype: the markers its reference needs, and the donors that carry them all."""
    rows = []
    for ph, ref in flat.items():
        need = sorted(ref)
        d = panel.donors_with(need)
        variant = sorted(set(need) & set(panel.PANEL_VARIANT_MARKERS))
        rows.append({"phenotype": ph, "n_markers": len(need), "n_donors": len(d),
                     "variant_markers": ",".join(variant), "donors": ",".join(d)})
    return pd.DataFrame(rows).sort_values(["n_donors", "phenotype"])


def requirements(rules: pd.DataFrame) -> dict:
    """phenotype -> (all_of, [any_of groups]) -- what the rule says a cell MUST show.

    A cosine is a similarity, not a requirement: ``allpos`` reads as "must be positive" but
    scoring it only adds a direction component, so a cell can win a reference without meeting
    it.  Measured on donor 6450 before this gate existed, ``T-regulatory`` requires CD4 and
    FOXP3 and 10.4% of the cells it was given met both; ``Imm-B-cell`` requires CD20 and 50.3%
    met it.  References with one dominant requirement held up -- T-CD4-helper 98.0%,
    T-CD8-cytotoxic 93.1% -- which is the pattern you would expect if nothing is enforced and
    dilution decides.  So requirements are extracted here and enforced separately from the
    ranking.

    OR-groups are read from the row that DECLARED them, not from the flattened result, and
    that distinction is load-bearing.  ``Epi-endocrine`` declares "any of B3TUBB|GCG|INS|SST";
    ``Endo-alpha`` then commits GCG to ``allpos`` and INS/SST to ``allneg``.  Read off the
    flattened tokens, only B3TUBB still carries ``anypos`` and the group collapses into a hard
    requirement for it -- which took alpha cells from thousands to 79.  Read off the declaring
    row, GCG is seen to satisfy the group and it is dropped.

    So, per group: satisfied and dropped if the phenotype commits any member to a positive;
    otherwise members it commits to a negative are removed; a group emptied that way is
    dropped rather than made unsatisfiable.
    """
    mk = [c for c in rules.columns if c not in ("parent", "phenotype")]
    row = {r["phenotype"]: {c: r[c] for c in mk if pd.notna(r[c])}
           for _, r in rules.iterrows()}
    par = dict(zip(rules["phenotype"], rules["parent"]))
    flat = flatten(rules)
    out = {}
    for ph in rules["phenotype"]:
        final = flat[ph]
        all_of = sorted(m for m, t in final.items() if t in ("pos", "allpos"))
        chain, cur = [], ph
        while cur and cur != "all":
            chain.append(cur)
            cur = par.get(cur)
            if cur in chain:
                break
        groups = []
        for anc in chain:                       # every level that declared an OR-group
            by_token: dict = {}
            for m, t in row.get(anc, {}).items():
                if OR_POS.match(str(t)):
                    by_token.setdefault(t, []).append(m)
            for _, members in sorted(by_token.items()):
                if any(final.get(m) in ("pos", "allpos") for m in members):
                    continue                    # the phenotype already commits to one
                members = [m for m in members if not is_negative(str(final.get(m, "")))]
                if members:
                    groups.append(sorted(members))
        out[ph] = (all_of, groups)
    return out


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
    markers = panel.markers(include_variant=True)      # the union: availability is per donor
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

    ev = evaluable_donors(flat)
    short = ev[ev["n_donors"] < 20]
    print(f"\nphenotypes measured on fewer than all 20 donors: {len(short)} of {len(ev)}")
    for _, r in short.iterrows():
        print(f"    {r['phenotype']:<22} n={r['n_donors']:2d} donors   "
              f"needs {r['variant_markers']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
