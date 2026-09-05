"""The broad label space, derived by code from two pre-existing declarations -- never from
Step-1 values, never from a result.

Sources
  * ``phenocycler/marker_registry.json``: ``kind`` (type / state / excluded / structural) and
    ``broad_compartment`` (Immune / Epithelial / Endothelial / Mesenchymal / Neural).  The
    registry is the marker authority (root CLAUDE.md).
  * ``phenotype_workflow.csv`` (the 37-reference rule table): used for exactly two things --
    the split of the registry's single Immune compartment into T_NK / B_Plasma / Myeloid by
    the branch in which a marker carries a positive token, and the endocrine hormones (the
    ``allpos`` markers of the Endo-beta / Endo-alpha / Endo-delta leaves).

Tiering rule (the PhenoCycler mirror of Xenium decision #17): identity markers are registry
``type`` markers only.  ``state`` markers never enter identity; they are appended after typing as
continuous per-cell columns.  Conflicts between the two sources are recorded, and a marker is
EXCLUDED when the rule table makes a definitional (``allpos``/``pos``) identity claim for it in a
different top-level branch than its registry compartment (CD66: registry Epithelial, rule-table
neutrophil).  OR-membership (``anypos``) in another branch is a note, not an exclusion
(B3TUBB: registry Neural, one of five OR-members of Epi-endocrine).

    cd core/phenocycler && ../../.venv/bin/python scripts/broadtype/declarations.py
writes ``scripts/broadtype/declarations.json`` and ``scripts/broadtype/rules_broad8.csv``;
``tests/test_declarations.py`` regenerates both and fails on any hand edit.
"""
from __future__ import annotations

import hashlib
import json
import pathlib
import sys

import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pins  # noqa: E402
import rules  # noqa: E402  (scripts/celltype/rules.py)

DECLARATIONS_JSON = HERE / "declarations.json"
RULES_BROAD8_CSV = HERE / "rules_broad8.csv"

#: lineage_native -> common lineage (phenospace.COMMON_LINEAGES).  The only hand-written map.
COMMON = {
    "Endocrine": "Endocrine",
    "Epithelial_nonendocrine": "Exocrine",
    "Vascular": "Vascular",
    "Stromal": "Stromal",
    "T_NK": "T_NK",
    "B_Plasma": "B_Plasma",
    "Myeloid": "Myeloid",
    "Neural": "Neural",
}
LABELS = tuple(COMMON)

#: registry compartment -> label(s); Epithelial splits by the hormone rule, Immune by branch.
COMPARTMENT_LABEL = {"Endothelial": "Vascular", "Mesenchymal": "Stromal", "Neural": "Neural"}
#: rule-table subtree roots that define the immune split.
IMMUNE_BRANCHES = {
    "T_NK": ("Imm-T-cell", "Imm-NK-cell"),
    "B_Plasma": ("Imm-B-cell", "Imm-plasma-cell"),
    "Myeloid": ("Imm-macrophage", "Imm-dendritic-cell", "Imm-neutrophil"),
}
ENDOCRINE_LEAVES = ("Endo-beta", "Endo-alpha", "Endo-delta")
#: top-level rule-table rows -> the label group they speak for (for the conflict check).
TOP_GROUP = {"Epithelial": {"Epithelial_nonendocrine", "Endocrine"},
             "Immune": {"T_NK", "B_Plasma", "Myeloid"},
             "Endothelial": {"Vascular"}, "Mesenchymal": {"Stromal"}, "Neural": {"Neural"}}
REASONS = {
    "Endocrine": "hormones are mutually exclusive subtype markers (panel.DECLARED_CONTROLS); held-out ablation would test the neighbourhood of a removed subtype, so Endocrine is gated by spatial coherence and stability",
    "Epithelial_nonendocrine": "no acinar- or ductal-specific marker is declared on this panel; mapped to the common lineage Exocrine, which the integration vocabulary already records as PhenoCycler's default sink",
    "Vascular": "registry Endothelial compartment; Podoplanin marks lymphatic endothelium and stays declared",
    "Stromal": "registry Mesenchymal compartment",
    "T_NK": "registry Immune markers whose positive tokens sit in the T-cell / NK branch of the rule table",
    "B_Plasma": "registry Immune markers whose positive tokens sit in the B / plasma branch",
    "Myeloid": "registry Immune markers whose positive tokens sit in the macrophage / dendritic / neutrophil branch; granulocytes folded in as the Xenium typer does",
    "Neural": "registry Neural compartment; single declared marker, so reported and not held-out-gated",
}
#: Negative constraints for rules_broad8.csv (registry names).  Taken from the top-level rows of
#: the rule table (Epithelial / Immune / Endothelial / Neural / Mesenchymal) plus the hormones
#: for the non-endocrine epithelium and the cross-immune negatives the table uses at leaf level.
NEGATIVES = {
    "Endocrine": ("CD68", "CD20", "CD3e"),
    "Epithelial_nonendocrine": ("CD68", "CD20", "CD3e", "INS", "GCG", "SST"),
    "T_NK": ("E_cadherin", "Ker8_18", "Pan_Cytokeratin", "EpCAM", "CD20", "CD68"),
    "B_Plasma": ("E_cadherin", "Ker8_18", "Pan_Cytokeratin", "EpCAM", "CD3e", "CD68"),
    "Myeloid": ("E_cadherin", "Ker8_18", "Pan_Cytokeratin", "EpCAM", "CD3e", "CD20"),
    "Vascular": ("E_cadherin", "Pan_Cytokeratin", "CD3e"),
    "Neural": ("CD31", "E_cadherin", "CD3e"),
    "Stromal": ("CD31", "CD68", "E_cadherin", "Pan_Cytokeratin", "CD3e"),
}


def _sha(path: pathlib.Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _children(table: pd.DataFrame) -> dict:
    kids: dict = {}
    for p, c in zip(table["parent"], table["phenotype"]):
        kids.setdefault(str(p), []).append(str(c))
    return kids


def _subtree(root: str, kids: dict) -> list:
    out, stack = [], [root]
    while stack:
        n = stack.pop()
        out.append(n)
        stack.extend(kids.get(n, []))
    return out


def _top_ancestor(name: str, parent_of: dict) -> str:
    while parent_of.get(name, "all") != "all":
        name = parent_of[name]
    return name


def derive() -> dict:
    reg = {m["name"]: m for m in pins.registry()["markers"]}
    kind = {n: m["kind"] for n, m in reg.items()}
    comp = {n: m.get("broad_compartment") for n, m in reg.items()}
    table = rules.load_rules(pins.RULES_37)
    # load_rules already renames the display columns to registry names via the frozen crosswalk
    marker_cols = [c for c in table.columns if c not in ("parent", "phenotype")]
    kids = _children(table)
    parent_of = dict(zip(table["phenotype"].astype(str), table["parent"].astype(str)))
    rows = {str(r["phenotype"]): r for _, r in table.iterrows()}

    def positive_rows(marker: str, tokens=("pos", "allpos", "anypos")) -> list:
        if marker not in marker_cols:
            return []
        return [ph for ph, r in rows.items()
                if isinstance(r.get(marker), str) and any(str(r[marker]).startswith(t) for t in tokens)]

    notes, excluded = [], {}
    declared: dict = {L: [] for L in LABELS}

    # 1. Endocrine: the hormones, from the Endo-* leaves' definitional tokens.
    hormones = sorted({m for m in marker_cols
                       if any(str(rows[leaf].get(m)) == "allpos" for leaf in ENDOCRINE_LEAVES if leaf in rows)})
    for m in hormones:
        if kind.get(m) != pins.KIND_TYPE:
            excluded[m] = f"hormone row marker is registry kind={kind.get(m)}"
            continue
        declared["Endocrine"].append(m)
    # 2. Registry compartments.
    for m, c in sorted(comp.items()):
        if kind[m] != pins.KIND_TYPE or not c or m in declared["Endocrine"]:
            continue
        # conflict check: a definitional claim in a different top-level branch excludes the marker
        strong = positive_rows(m, tokens=("pos", "allpos"))
        strong_groups = {g for ph in strong for g in TOP_GROUP.get(_top_ancestor(ph, parent_of), set())}
        weak = positive_rows(m, tokens=("anypos",))
        weak_groups = {g for ph in weak for g in TOP_GROUP.get(_top_ancestor(ph, parent_of), set())}
        if c == "Immune":
            pos = set(positive_rows(m))
            branches = [L for L, roots in IMMUNE_BRANCHES.items()
                        if any(ph in pos for root in roots if root in rows for ph in _subtree(root, kids))]
            if len(branches) == 1:
                declared[branches[0]].append(m)
            elif not branches:
                excluded[m] = "registry Immune marker with no positive token in any immune branch of the rule table"
            else:
                excluded[m] = f"registry Immune marker with positive tokens in more than one branch: {branches}"
            continue
        if c == "Epithelial":
            target = "Epithelial_nonendocrine"
        else:
            target = COMPARTMENT_LABEL[c]
        own = TOP_GROUP.get(c, set())
        foreign_strong = strong_groups - own
        if foreign_strong:
            excluded[m] = f"registry {c} but the rule table makes a definitional claim for it under {sorted(foreign_strong)} (rows {strong})"
            continue
        foreign_weak = weak_groups - own
        if foreign_weak:
            notes.append(f"{m}: registry {c}; also an OR-member under {sorted(foreign_weak)} (rows {weak}) -- kept, noted")
        declared[target].append(m)
    # 3. Type markers with no compartment and not a hormone.
    for m in sorted(kind):
        if kind[m] == pins.KIND_TYPE and not comp.get(m) and m not in declared["Endocrine"]:
            excluded[m] = "registry type marker with no broad_compartment" + (
                f"; rule-table positive rows {positive_rows(m)}" if positive_rows(m) else "")
    for m in sorted(kind):
        if kind[m] == pins.KIND_EXCLUDED:
            excluded[m] = "registry kind=excluded"

    labels, undeclared = {}, {}
    for L in LABELS:
        ms = sorted(declared[L])
        if not ms:
            undeclared[L] = ("no registry type marker survives the conflict rule for this label; "
                             "the label is not declared and no arm may emit it")
            continue
        labels[L] = {
            "common_lineage": COMMON[L],
            "markers": ms,
            "n_markers": len(ms),
            "held_out_gated": (len(ms) >= 2 and L != "Endocrine"),
            "negatives": list(NEGATIVES[L]),
            "reason": REASONS[L],
        }
    state = sorted(m for m, k in kind.items() if k == pins.KIND_STATE)
    structural = sorted(m for m, k in kind.items() if k == pins.KIND_STRUCTURAL)
    all_identity = sorted({m for L in labels for m in labels[L]["markers"]})
    bad = [m for m in all_identity if kind[m] != pins.KIND_TYPE]
    if bad:
        raise RuntimeError(f"non-type markers reached identity: {bad}")
    return {
        "method_version": "broadtype-declarations-2026-09-05",
        "rule": ("identity = registry kind=='type'; compartments from the registry; Immune split by the "
                 "rule-table branch holding a positive token; Endocrine = allpos markers of the Endo-* "
                 "leaves; definitional (pos/allpos) cross-branch claims exclude; OR-membership notes"),
        "labels": labels,
        "undeclared_labels": undeclared,
        "identity_markers": all_identity,
        "state_markers": state,
        "structural_markers": structural,
        "excluded": dict(sorted(excluded.items())),
        "notes": notes,
        "sources": {"registry": str(pins.REGISTRY.relative_to(pins.REPO)), "registry_sha256": _sha(pins.REGISTRY),
                    "rules": str(pins.RULES_37.relative_to(pins.REPO)), "rules_sha256": _sha(pins.RULES_37)},
    }


def rules_broad8(dec: dict) -> pd.DataFrame:
    """The 8-row scimap-format rule table (display names) for the cosine arm."""
    cw = rules.load_crosswalk(pins.RULES_37)             # display -> registry
    inv = {v: k for k, v in cw.items()}                   # registry -> display
    cols = [c for c in pd.read_csv(pins.RULES_37, nrows=0).columns if not str(c).startswith("Unnamed:")]
    out = []
    for L, spec in dec["labels"].items():
        row = {"parent": "all", "phenotype": L}
        for c in cols:
            row[c] = ""
        ms = spec["markers"]
        tok = "allpos" if len(ms) == 1 else "anypos"
        for m in ms:
            if m not in inv:
                raise RuntimeError(f"{L}: declared marker {m} has no column in the rule table's crosswalk")
            row[inv[m]] = tok
        for m in spec["negatives"]:
            if m in ms:
                raise RuntimeError(f"{L}: {m} is both declared and negative")
            row[inv[m]] = "allneg"
        out.append(row)
    return pd.DataFrame(out, columns=["parent", "phenotype", *cols])


def write() -> dict:
    dec = derive()
    DECLARATIONS_JSON.write_text(json.dumps(dec, indent=1, sort_keys=False) + "\n")
    df = rules_broad8(dec)
    # v4 layout the celltype loader reads: two unnamed leading columns
    df = df.rename(columns={"parent": "", "phenotype": ""})
    df.to_csv(RULES_BROAD8_CSV, index=False)
    return dec


def load() -> dict:
    return json.loads(DECLARATIONS_JSON.read_text())


if __name__ == "__main__":
    d = write()
    for L, s in d["labels"].items():
        print(f"{L:26s} -> {s['common_lineage']:10s} {s['n_markers']} markers: {s['markers']}  held_out_gated={s['held_out_gated']}")
    print("state markers:", len(d["state_markers"]))
    print("excluded:")
    for m, r in d["excluded"].items():
        print(f"   {m}: {r}")
    for n in d["notes"]:
        print("note:", n)
