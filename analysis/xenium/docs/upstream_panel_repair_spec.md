# Spec: repair the XeniumPanelExplorer subpanels for *specific* cell-type identity

**Hand-off document** for an agent working in `smith6jt-cop/XeniumPanelExplorer`.
Written 2026-06-08 from the consumer side (`Xenium_Analysis`) after the over-inclusive
panels caused catch-all mis-phenotyping. The goal is to make each cell-type subpanel
ship a clean, gene-level **identity / subclass / state** classification so a consumer
can score *specific* identity genes instead of a noisy panel mean.

## The problem (evidence)
The cell-type subpanels (`01–07`, `11–16` and the pancreas/thymus equivalents) are
**feature pools**, not identity-marker sets: each is padded with subclass markers,
activation-state genes, cell-cycle genes, and broadly-expressed/housekeeping genes.
Consequences measured downstream:
- The 226-gene `04_immune_B_plasma` panel contains cell-cycle genes (MKI67, TOP2A,
  CCNA2, CDK1, AURKB, BIRC5…). A proliferating **acinar** cell therefore scores high
  on the "B-cell" panel. Scoring the panel **mean** mislabeled ~11k non-B cells as B
  (≈0.9% expressed any canonical B marker).
- `panel_gene_roles.csv` only tier-classifies the **mixed** panels (08/09/10/14 +
  thymus 10/11/17/19/20). The **cell-type identity panels were never gene-tiered**, so
  their identity genes and their state/cycle/housekeeping padding are indistinguishable.

Net: "is this gene useful for *identity*?" is unanswered per-gene for the panels that
matter most. This spec fixes that.

## Goal
Extend the gene-level tier classification (`identity` / `subclass` / `state` /
`non_specific` / `epithelial_general`) to **every gene in every cell-type subpanel**,
grounded in the **Cell Ontology (CL)** hierarchy + the CellxGene `cell_type`
annotation already in the metadata, plus a curated canonical override. Emit a per-panel
**`identity_core`** (the `identity`-tier subset) that a consumer scores directly.

## Method — per gene, per panel (must be hierarchy-aware, not string matching)

1. **Use the Cell Ontology hierarchy.** The CellxGene `cell_type` terms ARE CL terms and
   form an IS-A graph. Parse CL (`cl-basic.obo`, OBO Foundry); map each gene's
   `cell_type` terms → CL IDs (by name + synonym). Each panel's target lineage has a CL
   **anchor** term; a term *supports* lineage L iff it is a **descendant of (or equal to)**
   L's anchor. This makes subclasses roll up correctly and places general terms right:
   - `pancreatic A cell` ⊂ `pancreatic endocrine cell` ⊂ `endocrine cell` (verified in CL).
   - acinar / ductal / islet / **keratinocyte** ⊂ `epithelial cell`. **The pancreas and
     thymus are epithelial organs** — do NOT discard epithelial/keratinocyte terms (ducts
     undergo squamous metaplasia; thymic Hassall's corpuscles are keratinized TEC).
   Suggested anchors: Endocrine=`pancreatic endocrine cell`; Acinar=`pancreatic acinar cell`;
   Ductal=`pancreatic ductal cell`; T=`T cell`; B=`B cell`(+`plasma cell`);
   NK=`natural killer cell`(+`innate lymphoid cell`); Myeloid=`myeloid leukocyte`
   (+`dendritic cell`,`mast cell`); Vascular=`endothelial cell`;
   Stromal=`fibroblast`/`smooth muscle cell`/`pericyte`/`mesenchymal cell`;
   Neural(peripheral)=`Schwann cell`/`peripheral nervous system neuron` (NOT generic CNS
   `neuron`); TEC=`thymic epithelial cell` (mTEC=`medullary…`, cTEC=`cortical…`).

2. **Assign a tier per gene** (dominant over its mapped terms; suggest ≥70% share):
   - **`identity`** — terms dominantly under the lineage's CL anchor → marks the lineage.
   - **`subclass`** — dominantly under a *subtype* anchor inside the lineage (e.g.
     `pancreatic A cell`→alpha; `CD8-positive, alpha-beta T cell`→CD8; `medullary thymic
     epithelial cell`→mTEC). Carry a `subtype_hint`.
   - **`state`** — primarily a pathway/program gene: present in a `cellchat_pathway`, or a
     member of a HALLMARK/REACTOME/KEGG/cell-cycle/IFN/stress/metabolism program panel, or
     annotated only to a *process* not a cell type. Cross-lineage; never an identity gene.
   - **`non_specific`** — `cell_type` maps to many lineages (broad terms: `leukocyte`,
     `hematopoietic cell`, `epithelial cell` [generic], `muscle cell`, `connective tissue
     cell`, `professional antigen presenting cell`), OR unannotated AND not canonical.
   - **`epithelial_general`** — under `epithelial cell` but NOT under a specific epithelial
     lineage anchor → shared epithelium (an embedding-extra, not a per-lineage mask).

3. **Canonical override (fills CL's gaps — REQUIRED).** ~43% of the 5K genes have *no*
   `cell_type` annotation, **including textbook markers** (IAPP, MAFA, NKX6-1, SLC30A8,
   AMY1A, CUZD1, KRT19, CHGA, TRAC/TRBC2, PTPRC). Maintain a curated, literature-backed
   canonical marker list per cell type that **forces `identity`** regardless of CL, and
   resolves multi-lineage mis-maps (MAFB = alpha-cell TF though also myeloid; SOX9 =
   ductal though CL lists other epithelia; PMP22 = Schwann though listed broadly). Flag
   each gene's `evidence ∈ {cellxgene, canonical, both}`.

4. **Specificity validation (data-backed).** The fat CSVs already carry
   `detection_pct_<run>`. For each panel's `identity_core`, confirm the genes are detected
   *predominantly* in that lineage's cells on the reference runs (or at least are not
   ubiquitously detected). Report genes that fail (candidate `non_specific`).

## Deliverables
1. **Extend `data/panel_gene_roles.csv` to cover ALL cell-type panels**, columns:
   `source, subpanel, gene, gene_tier {identity|subclass|state|non_specific|epithelial_general},
   lineage_hint, subtype_hint, cl_anchor, cl_terms, evidence {cellxgene|canonical|both}, note`.
2. **Per-panel `identity_core` list** (the `identity`-tier genes) — the consumer scores THESE
   for the lineage, not the whole panel. Expose in `panel_roles.csv` (`n_identity_genes`) and
   as `data/<...>/identity_core/<panel>.csv`.
3. **Vendor the inputs for reproducibility**: the CL release used (`cl-basic.obo` + version/date)
   and the CellxGene annotation source/version.
4. Keep the existing tier framework for the mixed/pathway panels (don't regress
   `panel_classification.md`); this spec *extends* it to the cell-type panels.

## Why this is the right fix
`panel_classification.md` already establishes identity/subclass/state (Morris 2019,
*Development*; Heumos 2023, *Nat Rev Genet*) and applied it to the disputed panels. This
spec applies the *same principle, gene-by-gene, to the identity panels*, using the CL
ontology so it is reproducible and not hand-wavy. Once shipped, a consumer scores each
panel's `identity_core` (absolute expression — NOT a z-scored panel mean) and gets clean,
specific, no-catch-all cell typing; the rich full panels remain available for embedding /
feature learning.

## Concrete errors found in the current annotation (fix these specifically)
Measured on the consumer side by deriving fine anchors from the CellxGene `cell_type`
annotation alone — these are why fine subtyping fails and must be corrected at the source:
- **Mis-filed master TF**: `ARX` (the **alpha**-cell master TF) is annotated dominantly to
  `pancreatic PP cell` → lands in **Gamma**, not Alpha. Pan-endocrine TFs `INSM1`/`ISL1` get
  forced into single subtypes. Fine endocrine subtype annotation is unreliable.
- **Gut/islet conflation**: anchoring Endocrine on `enteroendocrine cell` (intestine) pulls in
  gut genes — `CDX1/CDX2`, `FABP2`, `TPH1`, `DEFA6`, `MGAM`, `GUCY2C`, `SLC15A1`. Keep
  pancreatic islet (`pancreatic endocrine cell` + its A/B/D/PP/epsilon descendants) SEPARATE
  from gut enteroendocrine.
- **Housekeeping leakage**: ubiquitous genes are tagged to specific islet cells — `ATP5F1A`
  (ATP synthase), `U2AF1` (splicing), `TUBA4A` (tubulin), `KPNB1`, `CALM3`, `GLUD1`, `AAMP`.
  A true marker is not pan-expressed.
- **Uneven coverage / shared-gene capture**: NK accrues 121 "specific" genes (incl. shared
  cytotoxic GZMA/PRF1/EOMES/TBX21/ZAP70 that are also T) while generic `T cell` yields ~0
  unique → T and thymocyte stages come out empty. Resolve T vs NK and assign shared cytotoxic
  genes to the correct level.
- **Specificity can't be rescued by *overall* detection** in a dominant-cell-type tissue
  (acinar-dominated pancreas): even housekeeping shows low overall detection. Use **per-lineage
  detection enrichment** (detected in the target lineage ≫ others), not overall detection.

A literature-grounded, adversarially-verified DRAFT of the corrected per-subtype `identity_core`
(both tissues) accompanies this spec (consumer artifact `data/reference/identity_cores_verified.csv`)
— use it as a starting point, but the upstream agent + a human should sign off the biology.

## Reference implementation
A working consumer-side prototype of this exact method is in `Xenium_Analysis`:
`scripts/derive_identity_panels.py` (parses `cl-basic.obo`, maps terms via descendancy
under lineage anchors, merges a canonical core, emits `identity_panels_<tissue>.csv`,
`identity_audit_<tissue>.csv`, `identity_epithelial_general_<tissue>.csv`). Port/improve
that logic upstream so the classification ships *with the panels* rather than being
re-derived by every consumer.
