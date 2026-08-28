# Ground-truth anchor markers — for biological review

These are the **curated canonical anchors**: the genes that confidently label a cell as a
lineage/subtype. In the corrected design these define the scANVI **anchor cells** (treated
as ground truth) and override CL annotation gaps. A cell is anchored to a lineage when it
cleanly expresses that lineage's anchors (absolute expression). Everything else is
propagated by scANVI. **So the accuracy of cell typing rests on these lists — please vet.**

Critical panel constraint (T1D add-on 5K panel): the canonical hormones **INS, GCG, SST,
PPY, TTR** and most acinar zymogens (PRSS1/2, CPA1/2, CTRB1/2, CELA, REG) are **NOT on the
panel** → surrogates are used and flagged with ⚠. Genes I consider **weak/questionable
anchors are flagged 🚩** — these are where I most want your judgment.

---

## PANCREAS — broad-lineage anchors

| Lineage | Anchor genes | Notes / flags |
|---|---|---|
| **Endocrine** | CHGA, INSM1, ISL1, NEUROD1, PAX6, SCG2, SCGN, IAPP, MAFA, NKX6-1, PDX1, SLC30A8, ARX, MAFB, HHEX | pan-islet TFs + granins (solid) + β/α/δ surrogates. ⚠ INS/GCG/SST off-panel. 🚩 MAFB also myeloid; PDX1 also duct progenitor |
| **Exocrine** | AMY1A, CUZD1, KRT19, SOX9, HNF1B, CFTR, MUC1, TSPAN8 | ⚠ AMY1A/CUZD1 are the only acinar markers on-panel (zymogens off). 🚩 SOX9/KRT19 broad-epithelial (also ductal-progenitor) |
| **Immune_T** | CD3E, CD3G, CD2, CD5, CD7, PTPRC, TRAC, TRBC2, CD8A, CD8B, CD4 | strong T set. 🚩 PTPRC is pan-leukocyte (shouldn't anchor T alone); CD4 also low in myeloid |
| **Immune_B** | MS4A1, CD79A, CD79B, CD19, CD22, MZB1, TNFRSF17, PRDM1 | MS4A1/CD79A/CD19 = B; MZB1/TNFRSF17/PRDM1 = plasma. solid |
| **Myeloid** | CD68, CD163, CSF1R, MARCO, CD14, ITGAM, ITGAX, CD1C | macrophage/mono + DC. ⚠ LYZ/C1QA/C1QB off-panel. solid |
| **Stromal** | ACTA2, PDGFRB, RGS5, PDGFRA, THY1, CDH11 | pericyte/SMC + fibroblast/stellate. ⚠ COL1A1/DCN/LUM off-panel. 🚩 ACTA2 also myoepithelial; THY1 broad |
| **Vascular** | PECAM1, CDH5, CLDN5, PLVAP, KDR, ENG, ERG, TIE1, CD93 | strong endothelial. ⚠ VWF/FLT1 off-panel |
| **Neural** | MPZ, PMP22, NGFR, SOX10 | Schwann (peripheral glia). ⚠ PLP1/S100B off-panel. NES/S100A1 dropped (non-specific) |

### PANCREAS — fine-subtype anchors
| Subtype | Anchors | Flags |
|---|---|---|
| **Beta** | IAPP, MAFA, NKX6-1, PDX1, SLC30A8, ABCC8, KCNJ11, GLP1R, PCSK1 | ⚠ INS off-panel → IAPP primary. reasonable depth |
| **Alpha** | ARX, MAFB, GC, IRX1 | ⚠ GCG off-panel. 🚩 only 4 genes; MAFB cross-lineage (myeloid/δ) — **thin** |
| **Delta** | HHEX, LEPR, PCSK1, POU2F2 | ⚠ SST off-panel. 🚩🚩 **weakest anchor**: HHEX ok, but LEPR weak, PCSK1 shared w/ β, POU2F2 also B-cells |
| **Endocrine (residual)** | CHGA, INSM1, ISL1, NEUROD1, FEV, PAX6, SCG5, SCG2, SCGN | generic islet — for cells +endocrine but no β/α/δ subtype |
| **Gamma / Epsilon** | *(rule-based, not here)* | ⚠ PPY/GHRL: Gamma via AQP3+ETV1+pan-endocrine; Epsilon via GHRL. 🚩 ultra-rare, rule-driven |
| **Acinar** | AMY1A, CUZD1 | 🚩 only 2 on-panel (panel limitation) |
| **Ductal** | KRT19, SOX9, HNF1B, MUC1, CFTR, TSPAN8 | solid |
| **Stellate** | ACTA2, PDGFRB, RGS5, PDGFRA, THY1, CDH11 | (= broad Stromal; pancreas has only this stromal fine type) |
| **T_cells / NK / B / Plasma / Myeloid / DC / Mast / Endothelial / Schwann** | as broad + NK: NKG7,KLRD1,KLRF1,NCR3,KLRK1,PRF1,GZMA,KLRB1,FCGR3A; DC: ITGAX,CD1C,CLEC9A,CLEC10A,BATF3,IRF8,FLT3,LAMP3; Mast: KIT,CMA1,HDC,MS4A2,FCER1A | 🚩 NK shares GZMA/PRF1/NKG7 with CD8-effector T (NK↔CD8 ambiguity) |

---

## THYMUS — broad-lineage anchors

| Lineage | Anchor genes | Notes / flags |
|---|---|---|
| **Thymocyte** | CD3E, CD3G, TRAC, TRBC2, CD7, CD8A, CD8B, CD4, BCL11B, LEF1, TCF7 | ⚠ CD3D off-panel. strong T-lineage |
| **TEC** | EPCAM, KRT14, PSMB11, PRSS16, CCL25, FOXN1, AIRE, FEZF2, CTSL, LY75 | cTEC (PSMB11/PRSS16/CCL25/LY75) + mTEC (AIRE/FEZF2/KRT14). AIRE/FEZF2 promoted from upstream state panels. solid |
| **Myeloid** | CD68, CD163, CSF1R, MARCO, CD14, ITGAM, ITGAX, CD1C | solid |
| **B_lymphoid** | MS4A1, CD79A, CD79B, CD19, BANK1, PAX5, MZB1, TNFRSF17 | solid |
| **NK** | NKG7, KLRD1, KLRF1, NCR3, KLRK1, PRF1, GZMA, FCGR3A | 🚩 NK↔CD8 thymocyte ambiguity (shared cytotoxic genes) |
| **Vascular** | PECAM1, CDH5, CLDN5, PLVAP, KDR, ENG, TIE1, ERG, CD93, LYVE1, PROX1, PDPN | endothelial + lymphatic (LYVE1/PROX1/PDPN) |
| **Stromal** | COL5A1, COL5A2, COL6A1, CDH11, THY1, FAP, PDGFRA, ACTA2, MCAM, PDGFRB, RGS5 | ⚠ COL1A1/DCN off-panel → COL5/6. fibroblast + SMC/pericyte |
| **Neural** | MPZ, PMP22, NGFR, GLDN, PRPH | Schwann. ⚠ PLP1/S100B off-panel; NES/S100A1 dropped |

### THYMUS — fine-subtype anchors
| Subtype | Anchors | Flags |
|---|---|---|
| **DN_thymocyte** | CD34, NOTCH1, HES1, IL7R, BCL11B, CD7 | early/DN. solid |
| **DP_thymocyte** | CD8B, CCR9, LEF1, TCF7 | 🚩 DP defined biologically as CD4+CD8+; here via CD8B+DP-TFs |
| **SP4_thymocyte** | CD4, GATA3, CCR4, CD40LG, ICOS, KLF2, CCR7, FOXP1 | CD4 SP |
| **SP8_thymocyte** | CD8A, RUNX3, EOMES, CRTAM | CD8 SP. 🚩 thin (4) |
| **Treg** | FOXP3, IL2RA, CTLA4, IKZF2, TIGIT, TNFRSF18 | FOXP3 = definitive |
| **cTEC** | PSMB11, CTSL, PRSS16, CCL25, LY75, FOXN1, CXCL12, EPCAM | solid |
| **mTEC** | KRT14, AIRE, FEZF2, DSG3 | 🚩 thin (4); AIRE/FEZF2 low-detection in FFPE |
| **Fibroblast / Smooth_muscle / Endothelial / Lymphatic_EC / B / Plasma / Myeloid / DC / Mast / NK / Schwann** | as listed in MARKER_DICT | Mast/Lymphatic_EC/Schwann are rare; thin panels |

---

## My top concerns for your review
1. **Delta (pancreas)** 🚩🚩 — weakest anchor (SST off-panel; HHEX/LEPR/PCSK1/POU2F2 each have issues). Likely under/over-calls δ.
2. **Alpha (pancreas)** 🚩 — only ARX/MAFB/GC/IRX1 (GCG off-panel; MAFB cross-lineage).
3. **Acinar = AMY1A+CUZD1 only** — hard panel limit; everything else acinar is off-panel.
4. **PTPRC in the T anchor** — pan-leukocyte; should not anchor T on its own.
5. **NK ↔ CD8 cytotoxic-T overlap** (both tissues) — shared GZMA/PRF1/NKG7/KLR genes.
6. **mTEC / DP_thymocyte** — thin/indirect definitions; AIRE/FEZF2 low-detection.
7. **Gamma/Epsilon** — rule-based (AQP3/ETV1/GHRL), ultra-rare, not in the main anchor lists.

Off-panel canonicals replaced by surrogates throughout: INS→IAPP/MAFA, GCG→ARX/MAFB,
SST→HHEX, PPY→(AQP3/ETV1 rule), acinar zymogens→AMY1A/CUZD1, COL1A1/DCN→COL5/6,
VWF→PECAM1/CDH5, PLP1/S100B→MPZ/PMP22.
