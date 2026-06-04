#!/usr/bin/env python3
"""
Generate phenotype_rules_senior.csv (Junior's hierarchical format) for the
IO60 59-plex PhenoCycler panel. DRAFT for review.

Gating vocab (per Junior): pos = positive, allneg = all listed must be negative,
anypos = at least one of the listed positive, blank = no constraint. Positivity
itself is decided by per-marker, per-image thresholds in the phenotyper
(percentile/Otsu on the QuPath Cell:<marker>:Mean scale) — NOT hard-coded here.

Hierarchy is resolved most-specific-first; cells matching no leaf -> Unknown.
"""
import csv
from pathlib import Path

# Lineage markers used for gating (subset of the 59; checkpoint/state markers
# like PD-1/LAG3/TOX/Ki67/PCNA are overlays, not lineage gates).
MARKERS = ["INS","GCG","SST","Pan_Cytokeratin","EpCAM","Ker8_18",
           "Keratin_5","TP63","SOX2","CD31","CD34","Podoplanin","SMA","Vimentin",
           "Collagen_IV","B3TUBB","CD56","CD3e","CD8","CD4","FOXP3","CD20",
           "CD79a","CD68","CD163","Iba1","CD206","M2Gal3","CD11c","CD209",
           "CD11b","MPO","HLA_DR"]

# (parent_class, phenotype_name, {marker: gate})
RULES = [
    ("all", "Endocrine", {"INS":"anypos","GCG":"anypos","SST":"anypos"}),  # IAPP dropped (unreliable: GMM over-called it)
    ("Endocrine", "Beta cell",  {"INS":"pos"}),
    ("Endocrine", "Alpha cell", {"GCG":"pos"}),
    ("Endocrine", "Delta cell", {"SST":"pos"}),

    ("all", "Epithelial", {"INS":"allneg","GCG":"allneg","SST":"allneg",
                           "Pan_Cytokeratin":"anypos","EpCAM":"anypos","Ker8_18":"anypos",
                           "CD3e":"allneg","CD68":"allneg"}),
    ("Epithelial", "Ductal", {"Keratin_5":"anypos","TP63":"anypos","SOX2":"anypos"}),
    ("Epithelial", "Acinar", {"Pan_Cytokeratin":"pos","Keratin_5":"allneg","TP63":"allneg"}),

    ("all", "Immune", {"INS":"allneg","GCG":"allneg","SST":"allneg",
                       "CD3e":"anypos","CD20":"anypos","CD79a":"anypos","CD68":"anypos",
                       "CD163":"anypos","Iba1":"anypos","CD11c":"anypos","CD11b":"anypos",
                       "MPO":"anypos","CD56":"anypos"}),
    ("Immune", "T cell", {"CD3e":"pos"}),
    ("T cell", "CD8 T cell", {"CD8":"pos","CD4":"allneg"}),
    ("T cell", "CD4 T cell", {"CD4":"pos","CD8":"allneg"}),
    ("T cell", "Regulatory T cell", {"CD4":"pos","FOXP3":"pos"}),
    ("Immune", "B cell", {"CD3e":"allneg","CD20":"anypos","CD79a":"anypos"}),
    ("Immune", "Macrophage", {"CD3e":"allneg","CD68":"anypos","CD163":"anypos","Iba1":"anypos"}),
    ("Macrophage", "M2 macrophage", {"CD163":"anypos","CD206":"anypos","M2Gal3":"anypos"}),
    ("Immune", "Dendritic cell", {"CD3e":"allneg","CD11c":"anypos","CD209":"anypos","HLA_DR":"pos"}),
    ("Immune", "Granulocyte", {"CD3e":"allneg","CD68":"allneg","CD11b":"anypos","MPO":"anypos"}),
    ("Immune", "NK cell", {"CD3e":"allneg","CD56":"pos"}),

    ("all", "Endothelial", {"INS":"allneg","GCG":"allneg","SST":"allneg","CD3e":"allneg",
                            "Pan_Cytokeratin":"allneg","CD31":"anypos","CD34":"anypos",
                            "Podoplanin":"anypos"}),
    ("Endothelial", "Blood vessel", {"CD31":"anypos","CD34":"anypos"}),
    ("Endothelial", "Lymphatic", {"Podoplanin":"pos"}),

    ("all", "Stromal", {"INS":"allneg","GCG":"allneg","SST":"allneg","CD3e":"allneg",
                        "CD31":"allneg","Pan_Cytokeratin":"allneg",
                        "SMA":"anypos","Vimentin":"anypos","Collagen_IV":"anypos"}),
    ("Stromal", "Smooth muscle", {"SMA":"pos"}),
    ("Stromal", "Fibroblast", {"SMA":"allneg","Vimentin":"anypos","Collagen_IV":"anypos"}),

    ("all", "Neural", {"CD3e":"allneg","B3TUBB":"pos"}),
]


def main():
    out = Path(__file__).resolve().parent / "phenotype_rules_senior.csv"
    with open(out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["parent_class", "phenotype_name"] + MARKERS)
        for parent, name, gates in RULES:
            bad = set(gates) - set(MARKERS)
            if bad:
                raise SystemExit(f"{name}: unknown markers {bad}")
            w.writerow([parent, name] + [gates.get(m, "") for m in MARKERS])
    print(f"[done] wrote {out} ({len(RULES)} phenotype rules, {len(MARKERS)} gating markers)")


if __name__ == "__main__":
    main()
