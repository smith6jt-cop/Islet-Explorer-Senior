# Thymus Panel Audit — Xenium Prime 5K + 100 Custom (5K+100) vs. hIO 380-plex

This audit uses the canonical thymus-curated input (`panel_audit_thymus.csv`,
5,062 rows) and computes detection-weighted metrics against two thymus FFPE
Xenium runs (sample IDs `0041435` and `0041534`).

The structure mirrors the prior pancreas audit:
1. separate the 5K reference from the 100-gene custom add-on,
2. build cell-type / pathway subpanels of the 5K,
3. assess what is lost from immune analysis if hIO were used instead.

The previous heuristic reinstatement of pancreas-flagged exclusions is dropped:
the `exclude_recommended` column in `panel_audit_thymus.csv` is canonical.

---

## 1. Inputs separated: 5K reference vs. custom 100

Reconstructed against the canonical 5K reference:

| Set | Count |
|---|---:|
| 5K reference (canonical) | 5,001 |
| 5K-shared (5K ∩ audit) | 4,962 |
| 5K-pre-removed by user (5K − audit) | 39 |
| Custom add-on (audit − 5K reference) | 100 |
| Audit total | 5,062 |
| Audit-excluded (`exclude_recommended=yes`) | 106 |
| 5K-kept (in audit) | 4,864 |
| Custom-kept | 92 |
| **Working pool** | **4,956** |

The 39 pre-removed 5K genes include the 9 removed in the pancreas audit
(CMTM5, CYP2A6, F2, KLK4, KLK6, MSLN, MUC4, PAEP, PCDH11Y) plus 30 additional
removals that are tissue-restricted to lineages absent from thymus
(MAGEA/C cancer-testis antigens, serotonin receptors HTR1A/B, neural genes
ARR3/CRYGC/GLRA2/GUCY2D/NTSR2/NYX/POU3F4/PRDM13/SIM2/TRPC5, liver/adrenal
markers AHSG/AMH/CYP11B2/LBP/SULT2A1/TDO2, and several others).

### Custom 100 in thymus context

The 100-gene custom add-on was originally a T1D-GWAS / pancreas autoimmune set.
In thymus, 92 of those 100 genes survive the audit, 61 detect ≥1% in the
thymus runs, and the 6 TCR-pattern genes (TRAC, TRAV1-2, TRAV24, TRBC2,
TRBV6-4, TRGC2) are the most consequential additions for repertoire readout.
Together with TRDC from the 5K reference, the audit carries 7 TCR genes total
— the only TCR coverage available across either panel under comparison.

---

## 2. User's canonical exclusion calls (106 genes)

| Category | n | Note |
|---|---:|---|
| spermatogenesis | 18 | Cancer-testis / germline; not in thymus |
| hepatocyte | 16 | Tissue-restricted; not in thymus |
| cardiac_muscle | 13 | Tissue-restricted; not in thymus |
| photoreceptor_retinal | 12 | Tissue-restricted; not in thymus |
| skeletal_muscle | 10 | Tissue-restricted; not in thymus |
| oocyte_germline | 7 | Germ cell-restricted |
| pancreatic_islet_TF | 6 | ARX, MAFA, MAFB, NKX6-1, PDX1, SLC30A8 — β-cell-restricted; pan-neuroendocrine TFs (NEUROD1, ISL1, INSM1, FEV, CHGA, CHGB, SCG2, SCGN) deliberately retained for thymic medullary neuroendocrine cells |
| sex_chromosome | 6 | Y-linked / sex-restricted |
| nephron_specific | 5 | Renal tubule-restricted |
| bone_cartilage_tooth | 5 | Skeletal lineage-restricted |
| embryonic_pluripotency | 3 | Pluripotency-restricted |
| pancreatic_exocrine | 2 | AMY1A, CUZD1 — acinar/ductal-restricted |
| pancreatic_endocrine | 2 | GHRL, IAPP — islet hormone-restricted |
| chemosensory_receptor | 1 | TAS2R38 — gustatory receptor |

These exclusions remove the 5K's tissue-restricted antigen reporters wholesale.
A consequence: thymic AIRE-driven promiscuous TRA expression and thymic myoid
cells (TMCs, which express skeletal-muscle markers) are not readable from the
audit panel — those cells will appear in clusters but cannot be specifically
identified by lineage markers. This is the user's explicit curation decision
and is honored throughout.

---

## 3. Thymus detection profile of the working pool

`detection_pct_0041435` and `detection_pct_0041534` are derived from two
thymus runs. The max-of-two metric (NaN → 0) gives a per-gene detection
ceiling across the two samples.

| Threshold | Genes (kept pool, n=4,956) | % of pool |
|---|---:|---:|
| ≥0.5% | 3,268 | 65.9% |
| ≥1% | 2,721 | 54.9% |
| ≥5% | 701 | 14.1% |
| ≥10% | 189 | 3.8% |
| ≥20% | 40 | 0.8% |

For reference, the pancreas v3 audit reported 3,151 / 1,601 / 767 genes at
≥1% / ≥5% / ≥10% in the kept pool. Thymus has a markedly shallower
distribution at higher thresholds — only ~25% of the genes that crossed 10%
in pancreas do so in thymus. This is consistent with thymus being dominated
by thymocytes (small cells, lower RNA content per cell) versus the larger
acinar/ductal/endocrine cells in pancreas.

---

## 4. Twenty curated thymus subpanels

Built from the 4,956-gene working pool, using 10x `cell_type` token matching +
10x `cellchat_pathway` matching + curated marker overrides, intersected with
the kept pool. Each per-subpanel CSV carries the thymus detection columns.

| Subpanel | n_genes | ≥1% | ≥5% | ≥10% | median det_max |
|---|---:|---:|---:|---:|---:|
| 01_thymic_epithelial_cortical | 49 | 40 | 16 | 2 | 3.28% |
| 02_thymic_epithelial_medullary | 52 | 39 | 12 | 2 | 1.66% |
| 03_hassall_keratinized | 5 | 3 | 2 | 1 | 4.49% |
| 04_thymic_specialized_TEC | 152 | 58 | 12 | 3 | 0.46% |
| 05_thymic_myoid | 3 | 2 | 0 | 0 | 1.11% |
| 06_thymocyte_DN | 163 | 107 | 45 | 20 | 2.65% |
| 07_thymocyte_DP | 61 | 36 | 17 | 9 | 2.46% |
| 08_thymocyte_SP | 230 | 146 | 60 | 23 | 2.01% |
| 09_thymic_treg | 83 | 69 | 31 | 11 | 3.30% |
| 10_TCR_VDJ_machinery | 15 | 6 | 4 | 3 | 0.88% |
| 11_positive_negative_selection | 22 | 17 | 9 | 3 | 3.41% |
| 12_thymic_DC_myeloid | 470 | 334 | 113 | 37 | 2.30% |
| 13_thymic_B_plasma | 343 | 221 | 68 | 21 | 1.77% |
| 14_thymic_NK_ILC | 209 | 138 | 53 | 19 | 1.94% |
| 15_thymic_stroma_fibroblast | 400 | 173 | 35 | 11 | 0.73% |
| 16_thymic_endothelial_pericyte | 549 | 265 | 60 | 20 | 0.95% |
| 17_notch_il7_signaling | 27 | 17 | 8 | 2 | 1.64% |
| 18_crosstalk_RANKL_LTBR | 27 | 17 | 4 | 1 | 1.46% |
| 19_AIRE_TRA_program | 7 | 2 | 1 | 0 | 0.44% |
| 20_thymic_egress_S1P | 8 | 7 | 1 | 0 | 2.98% |
| 99_uncategorized | 3,216 | 1,736 | 408 | 100 | 1.18% |

Notes on specific subpanels:

- **05_thymic_myoid (3 genes)** has collapsed because the user excluded 9 of
  the 16 canonical TMC markers (MYOD1, MYOG, MYH1/2/4/13, ACTN2, DMD, RYR1)
  under `skeletal_muscle`. Only ACHE and CHRNA1 remain plus one cell_type
  token match. TMCs cannot be resolved from this panel.
- **19_AIRE_TRA_program (7 genes)** is narrow because most of the cross-tissue
  TRA reporters are excluded. The retained set is the autoimmune-disease
  canonical autoantigens (CHRNA1, GAD1, TG, MOG) plus AIRE, FEZF2, CIITA.
  Detection is low (median 0.44%) because mTECs are rare.
- **03_hassall_keratinized (5 genes)** is small but high-signal — KRT1, KRT14,
  TGM3, DSG1, DSG3 are all retained (not in the exclusion list) and at least
  three detect ≥1% in thymus, anchoring Hassall corpuscle identity.
- **04_thymic_specialized_TEC (152 genes)** is broad due to 10x cell_type
  tokens for tuft/neuroendocrine/ionocyte matching across many tissues; only
  1 gene detects ≥1% in thymus, reflecting the rarity of these subtypes.
- **10_TCR_VDJ_machinery (15 genes)** has 6 detecting ≥1%; TRAC at 40.7% and
  TRBC2 at 36.5% are the highest. TRAV/TRBV variable segments fire poorly
  because each is one rearrangement target among many.
- **The high-density thymocyte/T-cell subpanels (06–09)** dominate the
  panel's thymus-relevant signal: ≥1% counts of 107, 36, 146, and 69 genes
  in DN, DP, SP, and Treg subpanels respectively. This is the panel's
  strongest readout.
- **MHC class II machinery** is essentially absent in both the 5K+100 and
  hIO panels (see §6); thymic positive selection of CD4+ thymocytes cannot
  be inferred from class II expression.

---


---

## 5. Pathway re-annotation of the uncategorized residual (HALLMARK / Reactome / KEGG)

Mirroring the pancreas v3 audit, the 3,216-gene `99_uncategorized` pool was
re-annotated against three pathway libraries fetched from Enrichr:

| Library | n_sets | Source |
|---|---:|---|
| MSigDB Hallmark 2020 | 50 | maayanlab.cloud/Enrichr |
| Reactome 2022 | 1,818 | maayanlab.cloud/Enrichr |
| KEGG 2021 Human | 320 | maayanlab.cloud/Enrichr |

Per-gene HALLMARK / REACTOME / KEGG tags were written into
`subpanel_99_uncategorized_REANNOTATED.csv`. Of the 3,216 uncategorized genes,
2,082 carry at least one pathway tag and 1,134 carry none (lncRNAs, MHC-region
genes without annotated function, RNA-processing helpers, and recently
characterized loci).

### 5.1 Pathway subpanels (21–49)

Twenty-nine pathway-derived subpanels were built from the kept pool. Selection
strategy: rank candidate pathways by overlap size with the pathway-eligible
pool (uncategorized genes), then take the highest-overlap **thymus-relevant**
terms first. Thymus-relevance is determined by keyword priors (TCR, T cell,
MHC, antigen processing, NF-κB, interferon, JAK-STAT, IL-7, Notch, apoptosis,
p53, mTOR, PI3K, MAPK, cytokine, chemokine, immune, cell cycle, proteasome,
ubiquitin, RANK, complement, hematopoiesis, etc.). Counts mirror the pancreas
v3 audit: 8 HALLMARK + 12 Reactome + 9 KEGG = 29 subpanels.

| # | Source | Pathway | n_in_pool | ≥1% | ≥5% | ≥10% |
|---:|---|---|---:|---:|---:|---:|
| 21 | HALLMARK | mTORC1 Signaling | 94 | 74 | 23 | 5 |
| 22 | HALLMARK | p53 Pathway | 92 | 67 | 16 | 5 |
| 23 | HALLMARK | Apoptosis | 105 | 72 | 33 | 16 |
| 24 | HALLMARK | Hypoxia | 98 | 58 | 21 | 9 |
| 25 | HALLMARK | PI3K/AKT/mTOR Signaling | 68 | 56 | 23 | 9 |
| 26 | HALLMARK | Interferon Gamma Response | 120 | **102** | 35 | 13 |
| 27 | HALLMARK | TNF-alpha Signaling via NF-κB | 115 | 68 | 22 | 8 |
| 28 | HALLMARK | IL-2/STAT5 Signaling | 125 | 86 | 19 | 5 |
| 29 | Reactome | Immune System | 829 | **553** | 206 | **71** |
| 30 | Reactome | Innate Immune System | 418 | 304 | 116 | 38 |
| 31 | Reactome | Cytokine Signaling In Immune System | 417 | 287 | 108 | 40 |
| 32 | Reactome | GPCR Downstream Signaling | 301 | 93 | 28 | 9 |
| 33 | Reactome | Cell Cycle | 209 | 152 | 58 | 14 |
| 34 | Reactome | Adaptive Immune System | 293 | 212 | 86 | 24 |
| 35 | Reactome | Signaling By Interleukins | 286 | 192 | 75 | 26 |
| 36 | Reactome | Transcriptional Regulation By TP53 | 144 | 113 | 38 | 11 |
| 37 | Reactome | Cell Cycle, Mitotic | 164 | 119 | 46 | 10 |
| 38 | Reactome | Deubiquitination | 98 | 77 | 35 | 7 |
| 39 | Reactome | MAPK Family Signaling Cascades | 135 | 88 | 36 | 7 |
| 40 | Reactome | Class I MHC Antigen Processing & Presentation | 99 | 81 | 29 | 5 |
| 41 | KEGG | PI3K-Akt signaling pathway | 213 | 131 | 40 | 10 |
| 42 | KEGG | MAPK signaling pathway | 180 | 113 | 33 | 7 |
| 43 | KEGG | Cytokine-cytokine receptor interaction | 232 | 81 | 17 | 7 |
| 44 | KEGG | Chemokine signaling pathway | 128 | 88 | 35 | 12 |
| 45 | KEGG | JAK-STAT signaling pathway | 128 | 68 | 28 | 9 |
| 46 | KEGG | Apoptosis | 97 | 78 | 17 | 3 |
| 47 | KEGG | mTOR signaling pathway | 77 | 57 | 17 | 4 |
| 48 | KEGG | Ubiquitin mediated proteolysis | 62 | 53 | 13 | 0 |
| 49 | KEGG | Cell cycle | 79 | 55 | 21 | 4 |

Notes:

- **Member genes are intersected with the full thymus working pool**, not
  restricted to genes that were uncategorized. The `was_uncategorized` column
  in each per-subpanel CSV flags which members came from the residual vs.
  which were already assigned to curated subpanels 01–20.
- **Pathway subpanels overlap heavily** by design (many genes belong to
  multiple pathways). For example, mTORC1 / PI3K-Akt / mTOR / p53 / apoptosis
  share several effectors; the Reactome Immune System set is a superset of
  Innate, Adaptive, Cytokine-Signaling, and Interleukin-Signaling subsets.
  This is the same redundancy pattern as in the pancreas v3 audit.
- **Highest thymus-detection density**: Reactome Class I MHC (median 3.35%),
  KEGG Ubiquitin-mediated proteolysis (3.61%), Reactome Deubiquitination
  (3.11%), HALLMARK PI3K/AKT/mTOR (2.96%), HALLMARK IFNγ Response (2.93%).
  Class I MHC machinery and IFNγ-response programs fire strongly in thymus,
  consistent with active antigen presentation by TECs.
- **Reactome Immune System (829 genes; 553 ≥1% in thymus; 71 ≥10%)** is the
  single largest pathway subpanel and serves as the broad immune backbone.
  For more focused analysis, Adaptive Immune System (293, 212 ≥1%) is the
  thymus-relevant subset; Innate Immune (418) covers the thymic DC/macrophage
  component and the IFN response machinery.

### 5.2 Post-pathway residual

- **Picked up by pathway subpanels (21–49)**: 972 genes from the residual
- **Still residual after pathway pass (99c)**: 2,244 genes
  (`subpanel_99c_residual_after_pathway_pass.csv`; sorted by `n_pathways`)
- **Truly unannotated (99d)**: 1,134 genes with no HALLMARK / Reactome / KEGG
  membership (`subpanel_99d_truly_unannotated.csv`)

A quick look at the top-detected 99d entries — SOX2-OT (45.7%), LENG8 (32.5%),
NONO (24.2%), PRRC2A (20.1% — sits in the MHC region as BAT2), RNF187 (18.8%),
H3F3B (18.6%), NOTCH2NLA (13.5%) — confirms that the truly-unannotated set is
dominated by lncRNAs, RNA-processing factors, ubiquitin-machinery genes, and
recently characterized loci that are simply not in the public pathway
libraries. None of these omissions reflect panel-design failures; they reflect
documentation gaps in HALLMARK/Reactome/KEGG.

---

## 6. 5K+100 vs. hIO — set summary

| | n |
|---|---:|
| 5K+100 thymus working pool | 4,956 |
| hIO | 380 |
| Shared | 294 |
| Only 5K+100 | 4,662 |
| Only hIO | 86 |
| hIO genes detected ≥1% in thymus (shared subset of n=294) | 176 |

Immune classifier (curated lineage roster + TCR/IG regex + hIO annotation
tokens, applied to the 5K+100 ∪ hIO union):

| | n |
|---|---:|
| Immune-classified union | 380 |
| Shared (immune) | 219 |
| Only 5K+100 (immune) | 115 |
| Only hIO (immune) | 46 |
| 5K+100 immune detected ≥1% in thymus | 186 / 334 |
| 5K+100 immune detected ≥5% in thymus | 67 / 334 |

---

## 7. Module-by-module thymus coverage (with detection)

`canonical_thymus_markers_side_by_side.csv` and
`canonical_thymus_module_summary.csv` give the full 23-module table.

**Modules dominated by 5K+100:**

| Module | n | 5K+100 kept | hIO | ≥1% in thymus |
|---|---:|---:|---:|---:|
| TCR variable/constant segments | 9 | 7 | **0** | 3 |
| V(D)J recombination machinery | 9 | 7 | **0** | 3 |
| Cortical TEC (cTEC) | 9 | 9 | 3 | 7 |
| Medullary TEC (mTEC) | 12 | 11 | 4 | 9 |
| Notch signaling | 13 | 12 | 3 | 9 |
| Specialized TEC: tuft / neuroendocrine | 14 | 11 | 0 | 1 |
| Crosstalk (RANKL/LTβR/CD40) | 15 | 15 | 3 | 10 |
| AIRE/FEZF2 & TRAs | 9 | 7 | 1 | 2 |
| Negative-selection apoptosis | 11 | 9 | 1 | 5 |
| MHC class I peptide-loading complex | 15 | 7 | 2 | 7 |
| Thymocyte transcription factors | 12 | 12 | 5 | 12 |
| Treg lineage | 10 | 9 | 4 | 7 |
| Thymic DCs | 18 | 18 | 9 | 9 |

**Modules where hIO contributes:**

| Module | n | 5K+100 | hIO | hIO-unique |
|---|---:|---:|---:|---|
| Thymic B cells / plasma cells | 25 | 18 | 15 | 5 (IGHG1, IGHM, IGKC, IGLC3, JCHAIN) |
| Thymic NK / ILC | 15 | 14 | 12 | 1 (GNLY) |
| Thymic macrophages | 12 | 7 | 7 | 2 (APOE, C1QB) |
| CD3 / CD4 / CD8 coreceptors | 11 | 10 | 6 | 1 (CD3D) |

**Modules where both panels fail:**

| Module | n | 5K+100 | hIO | Shared gap |
|---|---:|---:|---:|---|
| MHC class II | 12 | 2 | 1 | All α/β chains missing: HLA-DRA, HLA-DRB1, HLA-DPA1/B1, HLA-DQA1/B1, HLA-DMA/B, HLA-DOA/B |
| Hassall corpuscle | 10 | 5 | 1 | IVL, KRT10, SPRR2A, LOR, FLG, S100A8 absent from both |
| Thymic myoid cells (TMC) | 16 | 2 | 0 | 9 audit-excluded; 5 not in either panel |

---

## 8. The "what is lost from hIO" picture — by thymus detection

Of 115 immune-classified genes only in 5K+100:

| Detection threshold | Lost from hIO |
|---|---:|
| ≥10% in thymus | 6 |
| ≥5% | 26 |
| ≥1% | 63 |
| ≥0.5% | 85 |

Top 15 by thymus detection:

| Gene | det_max (%) | Function |
|---|---:|---|
| B2M | 63.1 | MHC class I light chain (essential for positive selection) |
| TRAC | 40.7 | TCR α constant region |
| TRBC2 | 36.5 | TCR β constant region |
| IRF9 | 15.1 | Type I IFN signaling |
| LEF1 | 14.0 | T-cell transcription factor |
| IFI6 | 12.2 | IFN-induced effector |
| LAT | 9.7 | TCR signaling adaptor |
| CIITA | 9.5 | MHC class II master regulator |
| IFITM1 | 9.1 | IFN-induced |
| TYK2 | 8.8 | JAK family kinase |
| PSMB9 | 8.7 | Immunoproteasome subunit |
| TAP1 | 8.4 | MHC I peptide transporter |
| LY75 | 8.3 | DC endocytic receptor (DEC-205) |
| PSMB8 | 8.1 | Immunoproteasome subunit |
| ZAP70 | 7.7 | TCR signaling kinase |

Nearly all top-ranked losses are central to T-cell development and tolerance.
A panel restricted to hIO would not detect TCR repertoire (TRAC/TRBC2 at >35%
in thymus), MHC class I assembly (B2M at 63.1%, TAP1/2, PSMB8/9, ERAP1), or
the MHC class II master regulator (CIITA at 9.5%).

---

## 9. What hIO uniquely contributes (46 immune genes)

| Functional group | n | Genes |
|---|---:|---|
| Myeloid / phagocyte effector | 11 | AIF1, APOE, AREG, CORO1A, DGKA, MIS18BP1, MPEG1, PLD4, PRDX4, SSR4, VCAN |
| Immunoglobulin constant region | 8 | IGHG1, IGHG2, IGHG3, IGHG4, IGHGP, IGHM, IGKC, IGLC3 |
| Chemokine ligand | 6 | CCL2, CCL3, CCL3L1, CCL4, CXCL1, CXCL14 |
| Complement classical pathway | 3 | C1QB, C1R, C1S |
| GPCR / T-cell migration | 3 | GPR171, GPR183, SAMD3 |
| TCR/CD3 complex | 2 | CD3D, CD37 |
| Mast cell | 2 | CPA3, MCEMP1 |
| S100/alarmin | 2 | S100A9, S100B |
| Cell-cycle inhibitor | 2 | CDKN2C, CDKN2D |
| Cytokine receptor | 2 | FLT1, IL22RA2 |
| Singletons | 5 | HLA-B, JCHAIN, GNLY, IRF2, ARG1 |

For thymus medullary biology specifically, the consequential additions are:
**HLA-B** (the only classical class I heavy chain in either panel),
**CD3D** (closing the CD3 complex), the **Ig constant-region set** (medullary
plasma-cell isotyping and AID/class-switch tracking), and **C1QB/C1R/C1S**
(complement deposition on apoptotic thymocytes during negative selection).

These genes are in hIO but not in the audit, so their thymus detection cannot
be computed from these files. Their behavior in thymus would have to be
inferred from the user's downstream dataset (5,101-gene panel, three thymus
samples) directly.

---

## 10. Output files

```
/mnt/user-data/outputs/thymus_audit/
├── README.md                                              (this file)
├── thymus_working_panel.csv                               4,956-gene working pool
├── subpanel_summary.csv                                   all 50 subpanels × detection bins
├── pathway_subpanel_picks.csv                             ranking of pathway picks (subpanels 21-49)
├── canonical_thymus_markers_side_by_side.csv              module × gene × panel × thymus detection
├── canonical_thymus_module_summary.csv                    module-level counts + ≥1% in thymus
├── immune_ledger_5K100_vs_hIO_thymus.csv                  380-gene immune union with thymus detection
├── hIO_unique_for_thymus_immune_grouped.csv               46 hIO-only immune genes by function
└── subpanels/
    ├── subpanel_01_thymic_epithelial_cortical.csv         curated 01-20 (cell type / biology)
    ├── ...
    ├── subpanel_20_thymic_egress_S1P.csv
    ├── subpanel_21_HALLMARK_mTORC1_Signaling.csv          pathway-derived 21-49
    ├── ...
    ├── subpanel_49_KEGG_Cell_cycle.csv
    ├── subpanel_99_uncategorized.csv                      original residual (no pathway tags)
    ├── subpanel_99_uncategorized_REANNOTATED.csv          residual + HALLMARK/REACTOME/KEGG cols
    ├── subpanel_99c_residual_after_pathway_pass.csv       residual after pathway subpanels (2,244)
    └── subpanel_99d_truly_unannotated.csv                 no pathway membership at all (1,134)
```

---

## 11. Bottom line

Same conclusion as in the pancreas audit, with a wider margin for thymus:
**no off-the-shelf Xenium panel replaces the 5K+100 for thymus.** The
content-based gap is large and the detection-weighted gap is consistent —
the 5K+100 carries 63 immune-classified genes ≥1% in thymus that hIO does
not, and the top 15 of those losses are precisely the central T-cell
development / tolerance / MHC-I-assembly machinery.

hIO contributes a small, well-defined set of additions (HLA-B, CD3D, Ig
constant regions, C1Q complement) that would be candidate add-on targets
for a future panel redesign — alongside the persistent shared gaps in
classical MHC class II α/β chains, the rest of the class I heavy chains
(HLA-A/C/E/F/G), and broader TCR V-segment coverage.

Two thymus subpanels in this audit are degraded relative to a fully
informed panel design: **05_thymic_myoid** (3 genes; the user excluded the
skeletal-muscle markers) and **19_AIRE_TRA_program** (7 genes; cross-tissue
TRA reporters excluded). These are the operational consequences of the
exclusion calls and not panel-content failures per se.
