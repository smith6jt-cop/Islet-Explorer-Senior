**Eight gates, plus DAPI for segmentation.**

| Marker | Compartment |
|---|---|
| CD3e (C46) | Immune — T |
| CD20 (C30) | Immune — B |
| CD68 (C7) | Immune — macrophage |
| CD11b (C57) | Immune — broad myeloid (granulocyte, DC, NK) |
| E-cadherin (C11) | Epithelial (endocrine + exocrine) |
| CD31 (C5) | Endothelial |
| B3TUBB (C22) | Neural |
| Vimentin (C9) | Mesenchymal |

Four for immune is not padding — with no CD45, immune can't be one gate, so it's the minimum union that nets the populations actually present in pancreas. SMA is dropped because every SMA+ cell is also Vimentin+, so it adds nothing at compartment level; it's the first thing to add back if you want a positive mesenchymal call rather than a residual one. Hormones are dropped because endocrine sits under Epithelial.

**Gate order — immune moves first, reversing what I said earlier:**

1. CD3e / CD20 / CD68 / CD11b any-positive → **Immune**
2. E-cadherin+ → **Epithelial**
3. CD31+ → **Endothelial**
4. B3TUBB+ → **Neural**
5. Vimentin+ → **Mesenchymal**
6. → **Other**

Two constraints force this. CD31/PECAM1 is expressed on monocytes, neutrophils, NK, and naive T cells, so gating endothelium before immune miscalls leukocytes as endothelial — immune must precede CD31. And B3TUBB must follow E-cadherin, since endocrine cells are TUBB3+. My earlier ordering put immune late because CD56 needed endocrine and neural removed first; CD56 isn't in this set, so that constraint disappears and the CD31 one governs. Bonus: removing immune first also clears the rare E-cadherin+ DC out of the epithelial gate.

**RESTORE pairs — 8 targets.** Counterparts draw on the full panel; they don't have to be gates themselves, which is why EpCAM is available for B3TUBB (EpCAM-high selects exocrine, excluding the TUBB3+ endocrine confound).

```python
first_pass_pairs = [
    ["E-cadherin", "CD31"],
    ["CD31",       "E-cadherin"],
    ["CD3e",       "E-cadherin"],
    ["CD20",       "E-cadherin"],
    ["CD68",       "E-cadherin"],
    ["CD11b",      "E-cadherin"],
    ["B3TUBB",     "EpCAM"],
    ["Vimentin",   "E-cadherin"],
]
```

Every target here has real prevalence in pancreas, so the near-zero-prevalence trap doesn't apply to this set — the two remaining runtime issues are the `>50` scale floor and acinar autofluorescence. Eight pairs also runs fast enough to validate scaling and thresholds before committing to the full list.

**Falls to Other:** mast cells, plasma cells (CD20−; CD79a recovers them), Schwann cells, adipocytes, and vimentin-dim quiescent stellate.

Add INS (C14), GCG (C47), SST (C36) — three markers — if you want the endocrine/exocrine split resolved in the same pass rather than a second one. IAPP is redundant with INS at this level.