#!/usr/bin/env python3
"""Gene-level identity classification for the cell-type subpanels.

Background
----------
``scripts/build_panel_roles.py`` tier-classifies the *mixed* panels (the ones
that bundle program genes among a few lineage markers: pancreas 08/09/10/14,
thymus 11/17/19/20).  The cell-type *identity* panels themselves (pancreas
01-07/11-13/15-16, thymus 01-09/12-16) were never gene-tiered, so a consumer
could not tell a panel's identity markers from its subclass / state / padding
genes -- and scoring the panel *mean* mislabels off-lineage cells (e.g. a
proliferating acinar cell scoring high on the cell-cycle-padded B-cell panel).

This script fixes that.  Per gene, per cell-type panel, it assigns a tier::

    identity            terms dominantly under the panel lineage's CL anchor
    subclass            terms dominantly under one subtype anchor in the lineage
    state               a program gene (cellchat pathway / pathway-panel member)
                        that does not mark the lineage
    epithelial_general  under `epithelial cell` but not a specific epi lineage
    non_specific        maps to many lineages, or unannotated and not canonical

Method (hierarchy-aware, not string matching)
---------------------------------------------
The CellxGene ``cell_type`` annotations in ``reference_5k/xenium5k_genes.csv``
are Cell Ontology (CL) terms forming an IS-A graph.  We vendor a CL release
(``data/ontology/cl-basic.obo``) and, for each gene, map its terms to CL IDs
and test descendancy: a term *supports* lineage L iff it is a descendant of (or
equal to) one of L's CL anchors.  Subclasses roll up correctly because e.g.
`pancreatic A cell` is-a `pancreatic endocrine cell` is-a `endocrine cell`.

CL is silent for ~43% of the panel (incl. textbook markers IAPP, MAFA, KRT19,
CHGA, TRAC ...), so a curated, literature-backed canonical list per lineage
forces ``identity`` regardless of CL and resolves multi-lineage mis-maps
(MAFB = alpha-cell TF though CL also lists myeloid; SOX9 = ductal; PMP22 =
Schwann).  Each gene's ``evidence`` records whether the call came from
``cellxgene`` (CL), ``canonical`` (override), or ``both``.

Outputs (written by ``run()``; ``build_panel_roles.py`` merges the gene rows
into the unified ``data/panel_gene_roles.csv`` and adds ``n_identity_genes``
to ``data/panel_roles.csv``):

  * cell-type gene rows (returned to the orchestrator)             -> panel_gene_roles.csv
  * data/identity_core/<panel>.csv          (shared panels)
    data/tissues/<t>/identity_core/<panel>.csv (tissue panels)    identity-tier genes only
  * data/identity_audit_<tissue>.csv        specificity check on each identity_core
  * data/identity_epithelial_general_<tissue>.csv   shared-epithelium genes

Run `python3 scripts/build_panel_roles.py` to regenerate everything (it merges
the mixed-panel rows into the unified file). Running this module standalone
(`python3 scripts/derive_identity_panels.py`) writes only the identity_core /
audit / epithelial-general outputs and prints a summary — useful while iterating
on the classifier without rebuilding panel_roles.

Stdlib only; the classification is curated metadata, not a numeric recompute.
"""
import csv
import os
import re
from functools import lru_cache

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
SHARED = os.path.join(DATA, "subpanels_shared")
OBO = os.path.join(DATA, "ontology", "cl-basic.obo")
REFERENCE = os.path.join(DATA, "reference_5k", "xenium5k_genes.csv")

TISSUE_DIR = {
    "pancreas": os.path.join(DATA, "tissues", "pancreas"),
    "thymus":   os.path.join(DATA, "tissues", "thymus"),
}

# Tier vocabulary.
IDENTITY, SUBCLASS, STATE, NONSPEC, EPIGEN = (
    "identity", "subclass", "state", "non_specific", "epithelial_general")

# Decision thresholds (share of a gene's *mapped* CL terms). The spec suggests
# a dominant >=70% share; kept as named constants so they are tunable.
SUPPORT_MIN = 0.70   # >= this share supporting the lineage anchor -> marks it
SUBTYPE_MIN = 0.70   # >= this share of supporting terms on one subtype -> subclass
EPI_MIN     = 0.70   # >= this share under `epithelial cell` (off-lineage) -> epi_general
UBIQUITOUS_PCT = 50.0  # identity gene detected in >= this % of all cells = ubiquitous flag


# ---------------------------------------------------------------------------
# Cell Ontology parsing
# ---------------------------------------------------------------------------
def parse_obo(path):
    """Parse an OBO file into (parents, name2id, id2name).

    Only CL terms and their is_a edges are kept; obsolete terms are dropped.
    name2id maps each term's primary name plus its EXACT/NARROW synonyms to the
    CL id (primary names win on collision -- they are loaded last).
    """
    parents, name2id, id2name = {}, {}, {}
    names_first = {}  # primary names, applied after synonyms so they win
    with open(path, encoding="utf-8") as fh:
        blocks = fh.read().split("\n\n")
    for block in blocks:
        if not block.startswith("[Term]"):
            continue
        cid = nm = None
        pars, syns, obsolete = set(), [], False
        for ln in block.splitlines():
            if ln.startswith("id: CL:"):
                cid = ln[4:].strip()
            elif ln.startswith("name: "):
                nm = ln[6:].strip()
            elif ln.startswith("is_a: CL:"):
                pars.add(ln.split()[1])
            elif ln.startswith("synonym: "):
                m = re.match(r'synonym: "(.*?)" (EXACT|NARROW)', ln)
                if m:
                    syns.append(m.group(1).strip())
            elif ln.startswith("is_obsolete: true"):
                obsolete = True
        if not cid or obsolete:
            continue
        parents[cid] = pars
        id2name[cid] = nm
        for s in syns:
            name2id.setdefault(s.lower(), cid)
        if nm:
            names_first[nm.lower()] = cid
    name2id.update(names_first)  # primary names override synonyms
    return parents, name2id, id2name


PARENTS, NAME2ID, ID2NAME = parse_obo(OBO)


@lru_cache(maxsize=None)
def ancestors(cid):
    """Transitive is_a closure of cid, inclusive of cid itself."""
    out = {cid}
    for p in PARENTS.get(cid, ()):
        out |= ancestors(p)
    return out


def anchor(name):
    """Resolve a CL term name to its id, failing loudly if absent.

    Anchors are the backbone of every lineage definition, so a typo or an
    ontology rename must surface immediately rather than silently widen a mask.
    """
    cid = NAME2ID.get(name.lower())
    if cid is None:
        raise KeyError("CL term not found in %s: %r" % (os.path.basename(OBO), name))
    return cid


def cl_data_version():
    with open(OBO, encoding="utf-8") as fh:
        for ln in fh:
            if ln.startswith("data-version:"):
                return ln.split(":", 1)[1].strip()
            if ln.startswith("[Term]"):
                break
    return "unknown"


# ---------------------------------------------------------------------------
# Lineage definitions: CL anchors + subtype anchors per cell-type lineage.
# `anchors`   : a term supports the lineage iff descendant-or-equal of any.
# `subtypes`  : {hint: [anchor names]} used to call `subclass` and tag subtype_hint.
# `epithelial`: lineage sits under `epithelial cell` -> off-lineage epithelial
#               terms route to `epithelial_general` (the spec's "embedding-extra").
# `epi_general_panel`: the dedicated pan-epithelial panel (16); generic-epithelial
#               genes are epithelial_general, only canonical pan-epi -> identity.
# ---------------------------------------------------------------------------
class Lineage:
    def __init__(self, key, anchor_names, subtypes=None, epithelial=False,
                 epi_general_panel=False):
        self.key = key
        self.anchors = [anchor(n) for n in anchor_names]
        self.anchor_names = anchor_names
        self.subtypes = {h: [anchor(n) for n in names]
                         for h, names in (subtypes or {}).items()}
        self.epithelial = epithelial
        self.epi_general_panel = epi_general_panel


LINEAGES = {
    "Endocrine": Lineage("Endocrine", ["pancreatic endocrine cell"], epithelial=True,
        subtypes={"alpha": ["pancreatic A cell"], "beta": ["type B pancreatic cell"],
                  "delta": ["pancreatic D cell"], "PP": ["pancreatic PP cell"],
                  "epsilon": ["pancreatic epsilon cell"]}),
    "Exocrine": Lineage("Exocrine", ["pancreatic acinar cell", "pancreatic ductal cell"],
        epithelial=True,
        subtypes={"acinar": ["pancreatic acinar cell"], "ductal": ["pancreatic ductal cell"]}),
    "T": Lineage("T", ["T cell"],
        subtypes={"CD4": ["CD4-positive, alpha-beta T cell"],
                  "CD8": ["CD8-positive, alpha-beta T cell"],
                  "Treg": ["regulatory T cell"], "gdT": ["gamma-delta T cell"],
                  "NKT": ["mature NK T cell"]}),
    "B": Lineage("B", ["B cell", "plasma cell"],
        subtypes={"plasma": ["plasma cell"], "memory": ["memory B cell"],
                  "naive": ["naive B cell"]}),
    "NK": Lineage("NK", ["natural killer cell", "innate lymphoid cell"],
        subtypes={"ILC2": ["group 2 innate lymphoid cell"]}),
    "Myeloid": Lineage("Myeloid", ["myeloid leukocyte", "dendritic cell", "mast cell"],
        subtypes={"macrophage": ["macrophage"], "monocyte": ["monocyte"],
                  "DC": ["conventional dendritic cell"],
                  "pDC": ["plasmacytoid dendritic cell"], "mast": ["mast cell"]}),
    "Granulocyte": Lineage("Granulocyte", ["granulocyte"],
        subtypes={"neutrophil": ["neutrophil"], "eosinophil": ["eosinophil"],
                  "basophil": ["basophil"]}),
    "Vascular": Lineage("Vascular", ["endothelial cell"],
        subtypes={"lymphatic": ["endothelial cell of lymphatic vessel"],
                  "vascular": ["endothelial cell of vascular tree"]}),
    "Mural": Lineage("Mural", ["pericyte", "smooth muscle cell"],
        subtypes={"pericyte": ["pericyte"], "SMC": ["smooth muscle cell"]}),
    "Fibroblast": Lineage("Fibroblast", ["fibroblast", "mesenchymal cell"],
        subtypes={"stellate": ["pancreatic stellate cell"],
                  "myofibroblast": ["myofibroblast cell"]}),
    "Neural": Lineage("Neural", ["Schwann cell", "peripheral nervous system neuron"],
        subtypes={"neuron": ["peripheral nervous system neuron", "enteric neuron"],
                  "glia": ["Schwann cell"]}),
    "Epithelial": Lineage("Epithelial", ["epithelial cell"], epithelial=True,
        epi_general_panel=True),
    # Thymus
    "TEC": Lineage("TEC", ["epithelial cell of thymus"], epithelial=True,
        subtypes={"cTEC": ["cortical thymic epithelial cell"],
                  "mTEC": ["medullary thymic epithelial cell"],
                  "tuft": ["tuft cell"], "neuroendocrine": ["neuroendocrine cell"],
                  "ionocyte": ["ionocyte"]}),
    "TEC_Hassall": Lineage("TEC_Hassall", ["keratinocyte", "epithelial cell of thymus"],
        epithelial=True, subtypes={"keratinocyte": ["keratinocyte"]}),
    "Myoid": Lineage("Myoid", ["muscle cell"],
        subtypes={"skeletal": ["skeletal muscle cell"]}),
    "Thymocyte": Lineage("Thymocyte", ["T cell", "thymocyte"],
        subtypes={"DN": ["double negative thymocyte"],
                  "DP": ["double-positive, alpha-beta thymocyte"],
                  "Treg": ["regulatory T cell"], "CD8": ["CD8-positive, alpha-beta T cell"],
                  "CD4": ["CD4-positive, alpha-beta T cell"]}),
}

# (source, stem) -> lineage key.  Only cell-type *identity* panels appear here;
# the mixed / pathway / residual panels are owned by build_panel_roles.py.
PANEL_LINEAGE = {
    ("shared", "03_immune_T_cell"): "T",
    ("shared", "04_immune_B_plasma"): "B",
    ("shared", "05_immune_NK_ILC"): "NK",
    ("shared", "06_immune_myeloid"): "Myeloid",
    ("shared", "07_immune_granulocyte"): "Granulocyte",
    ("shared", "11_endothelial_vascular"): "Vascular",
    ("shared", "12_pericyte_smooth_muscle"): "Mural",
    ("shared", "15_neural_glial"): "Neural",
    ("shared", "16_epithelial_general"): "Epithelial",
    ("pancreas", "01_pancreas_endocrine"): "Endocrine",
    ("pancreas", "02_pancreas_exocrine"): "Exocrine",
    ("pancreas", "13_fibroblast_stellate"): "Fibroblast",
    ("thymus", "01_thymic_epithelial_cortical"): "TEC",
    ("thymus", "02_thymic_epithelial_medullary"): "TEC",
    ("thymus", "03_hassall_keratinized"): "TEC_Hassall",
    ("thymus", "04_thymic_specialized_TEC"): "TEC",
    ("thymus", "05_thymic_myoid"): "Myoid",
    ("thymus", "06_thymocyte_DN"): "Thymocyte",
    ("thymus", "07_thymocyte_DP"): "Thymocyte",
    ("thymus", "08_thymocyte_SP"): "Thymocyte",
    ("thymus", "09_thymic_treg"): "Thymocyte",
    ("thymus", "12_thymic_DC_myeloid"): "Myeloid",
    ("thymus", "13_thymic_B_plasma"): "B",
    ("thymus", "14_thymic_NK_ILC"): "NK",
    ("thymus", "15_thymic_stroma_fibroblast"): "Fibroblast",
    ("thymus", "16_thymic_endothelial_pericyte"): "Vascular",
}


# ---------------------------------------------------------------------------
# Canonical override.  CANONICAL[lineage_key] = {gene: subtype_hint}.
# Every listed gene is forced to `identity`; subtype_hint ("" = pan-lineage)
# is carried through so a consumer can still subtype an identity call.
# Literature-backed textbook markers, deliberately conservative; a gene only
# takes effect if it is actually present in the panel being classified.
# ---------------------------------------------------------------------------
def _mk(pan, **subs):
    d = {g: "" for g in pan}
    for hint, genes in subs.items():
        for g in genes:
            d[g] = hint
    return d


CANONICAL = {
    "Endocrine": _mk(
        ["CHGA", "CHGB", "SCG2", "SCG3", "SCG5", "PCSK1", "PCSK2", "PCSK1N",
         "PAX6", "NKX2-2", "NEUROD1", "ISL1", "INSM1", "FEV", "RFX6", "PAX4"],
        alpha=["GCG", "ARX", "MAFB", "IRX2", "TTR"],
        beta=["INS", "IAPP", "MAFA", "NKX6-1", "PDX1", "SLC30A8", "DLK1", "UCN3",
              "G6PC2", "ADCYAP1", "HADH", "ERO1B", "ABCC8", "KCNJ11"],
        delta=["SST", "HHEX", "RBP4"], PP=["PPY"], epsilon=["GHRL"]),
    "Exocrine": _mk(
        [],
        acinar=["PRSS1", "PRSS2", "CPA1", "CPA2", "CPB1", "CELA2A", "CELA3A",
                "CELA3B", "CTRB1", "CTRB2", "CTRC", "CTRL", "PNLIP", "PNLIPRP1",
                "CLPS", "PLA2G1B", "AMY2A", "AMY2B", "AMY1A", "REG1A", "REG1B",
                "REG3A", "REG3G", "SPINK1", "CEL", "RBPJL", "PTF1A", "BHLHA15",
                "CUZD1", "GP2", "SYCN"],
        ductal=["KRT19", "KRT7", "SOX9", "CFTR", "MMP7", "TFF1", "TFF2", "TFF3",
                "CEACAM6", "MUC1", "SPP1", "ANXA4", "HNF1B", "ONECUT2", "CLDN1",
                "CLDN4", "AQP1", "SLC4A4", "KRT23", "VTCN1"]),
    "T": _mk(
        ["CD3D", "CD3E", "CD3G", "CD2", "CD5", "CD6", "CD7", "TRAC", "TRBC1",
         "TRBC2", "LCK", "LAT", "ZAP70", "THEMIS", "CD28", "TRAT1", "ITK",
         "CD247", "LCP2", "SKAP1"],
        CD4=["CD4"], CD8=["CD8A", "CD8B"],
        Treg=["FOXP3", "IL2RA", "CTLA4", "IKZF2"]),
    "B": _mk(
        ["CD19", "MS4A1", "CD79A", "CD79B", "CD22", "PAX5", "EBF1", "TCL1A",
         "FCRL1", "FCRL5", "VPREB3", "BANK1", "BLK", "CD72", "CR2", "FCER2",
         "TNFRSF13B", "SPIB"],
        plasma=["SDC1", "PRDM1", "XBP1", "MZB1", "JCHAIN", "TNFRSF17", "DERL3",
                "TXNDC5", "FKBP11", "CD38", "IGHG1", "IGKC", "IGHM", "SLAMF7",
                "POU2AF1"]),
    "NK": _mk(
        ["NCR1", "NCR3", "KLRD1", "KLRF1", "KLRB1", "KLRC1", "KLRC2", "KLRC3",
         "KLRK1", "KIR2DL1", "KIR2DL3", "KIR2DL4", "KIR3DL1", "KIR3DL2", "GNLY",
         "NKG7", "PRF1", "NCAM1", "FCGR3A", "EOMES", "TBX21", "TYROBP", "KLRG1",
         "CD244"],
        ILC2=["IL7R", "RORC", "GATA3", "KIT", "IL1RL1"]),
    "Myeloid": _mk(
        ["CD68", "CD14", "CD163", "LYZ", "CSF1R", "ITGAM", "ITGAX", "FCGR1A",
         "AIF1", "FCER1G", "TYROBP", "CD33", "CEBPA", "CEBPB", "SPI1", "LST1",
         "SIRPA"],
        macrophage=["MRC1", "MARCO", "MSR1", "C1QA", "C1QB", "C1QC", "TREM2",
                    "APOE", "MERTK"],
        monocyte=["FCN1", "VCAN", "CCR2"],
        DC=["FLT3", "ZBTB46", "CLEC9A", "XCR1", "BATF3", "IRF8", "CD1C",
            "FCER1A", "CLEC10A"],
        pDC=["LILRA4", "IL3RA", "CLEC4C", "TCF4", "IRF7", "BCL11A"],
        mast=["TPSAB1", "TPSB2", "CPA3", "MS4A2", "HDC", "CMA1"]),
    "Granulocyte": _mk(
        ["MPO", "ELANE", "PRTN3", "AZU1", "CTSG", "LCN2", "LTF", "CAMP"],
        neutrophil=["FCGR3B", "CSF3R", "FUT4", "CEACAM3", "S100A8", "S100A9",
                    "S100A12", "CXCR1", "CXCR2", "FFAR2", "MMP8", "MMP9"],
        eosinophil=["SIGLEC8", "IL5RA", "CCR3", "PRG2", "EPX", "CLC", "RNASE2",
                    "RNASE3"],
        basophil=["ENPP3", "GATA2"]),
    "Vascular": _mk(
        ["PECAM1", "CDH5", "VWF", "CLDN5", "KDR", "FLT1", "TEK", "ERG", "FLI1",
         "EGFL7", "CD34", "ENG", "EMCN", "PLVAP", "CLEC14A", "ROBO4", "TIE1",
         "ESAM", "SOX17", "SOX18", "NOTCH4", "CALCRL", "RAMP2", "CD93", "ICAM2"],
        lymphatic=["LYVE1", "PROX1", "PDPN", "FLT4", "CCL21", "MMRN1"]),
    "Mural": _mk(
        [],
        pericyte=["RGS5", "KCNJ8", "ABCC9", "NOTCH3", "CSPG4", "PDGFRB",
                  "HIGD1B", "COX4I2", "MCAM", "ANPEP"],
        SMC=["ACTA2", "MYH11", "TAGLN", "CNN1", "DES", "MYLK", "MYL9", "LMOD1",
             "PLN", "ACTG2", "MYOCD"]),
    "Fibroblast": _mk(
        ["PDGFRA", "COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "COL6A3",
         "LUM", "DCN", "THY1", "FAP", "POSTN", "FN1", "COL5A1", "COL5A2", "DPT",
         "MMP2", "SFRP2", "GSN", "COL14A1", "FBLN1", "MFAP5", "PDPN"],
        stellate=["RGS5", "GFAP", "NES", "SPARC"]),
    "Neural": _mk(
        ["PLP1", "MPZ", "MBP", "S100B", "PMP22", "SOX10", "GFAP", "NGFR",
         "CDH19", "ERBB3", "FOXD3", "PRX", "MAL", "PMP2", "EGR2", "CRYAB",
         "SCN7A"],
        neuron=["ELAVL4", "ELAVL3", "ELAVL2", "PHOX2B", "PHOX2A", "RET",
                "STMN2", "PRPH", "UCHL1", "NEFL", "NEFM", "NEFH", "TH", "DBH",
                "CHAT", "SNAP25", "SYP", "ASCL1", "TUBB3", "NRXN1", "NRXN3"]),
    "Epithelial": _mk(
        ["EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "CLDN3", "CLDN4", "CLDN7",
         "TACSTD2", "KRT7", "ELF3", "MUC1", "CEACAM5"]),
    "TEC": _mk(
        ["EPCAM", "KRT8", "KRT18", "KRT5", "KRT14", "KRT17", "KRT15", "KRT19",
         "CDH1", "TP63", "PAX1", "PAX9", "FOXN1", "DLL4"],
        cTEC=["PSMB11", "PRSS16", "LY75", "CCL25", "SLC46A2", "ENPEP"],
        mTEC=["AIRE", "FEZF2", "CCL21", "SPIB", "TNFRSF11A", "CD40", "IVL",
              "KRT1", "KRT10"],
        tuft=["POU2F3", "TRPM5", "GNAT3", "DCLK1", "AVIL", "GNG13", "SOX9"],
        neuroendocrine=["CHGA", "CHGB", "NEUROD1", "INSM1", "ASCL1", "SYP"],
        ionocyte=["FOXI1", "CFTR", "ASCL3", "ATP6V0A4"],
        microfold=["GP2", "CCL20", "SPIB", "TNFRSF11A"]),
    "TEC_Hassall": _mk(
        ["EPCAM", "KRT5", "KRT14"],
        keratinocyte=["KRT1", "KRT10", "IVL", "SPINK5", "FLG", "DSG1", "DSG3",
                      "KRT6A", "KRT16", "SBSN", "LOR", "TGM1", "CALML5",
                      "KRTDAP", "CDSN", "DSP"]),
    "Myoid": _mk(
        [],
        skeletal=["DES", "MYOG", "MYH1", "MYH2", "MYH3", "TTN", "ACTN2",
                  "CHRNA1", "CHRNG", "RYR1", "MYL1", "ACTA1", "TNNT3", "TNNC2",
                  "NEB", "MYF5", "MYOD1"]),
    "Thymocyte": _mk(
        ["CD3D", "CD3E", "CD3G", "CD2", "CD7", "PTCRA", "TRAC", "LCK", "ZAP70",
         "TRBC2"],
        DN=["CD34", "KIT", "IL7R", "CD44", "RAG1", "RAG2", "DNTT", "LYL1"],
        DP=["CD1A", "CD1B", "CD1E", "RORC", "AQP3"],
        CD8=["CD8A", "CD8B"], CD4=["CD4"],
        Treg=["FOXP3", "IL2RA", "CTLA4", "IKZF2"]),
}

# Negative overrides: gene -> set(lineage keys) where CL may map it but it is
# NOT that lineage's identity (resolves the multi-lineage mis-maps the spec
# calls out).  Applied only if the gene is not canonical-positive for the panel.
CANONICAL_NEGATIVE = {
    "MAFB": {"Myeloid"},   # alpha-cell TF; CL also lists myeloid
    "SOX9": {"Neural"},    # ductal/chondro TF; not neural identity here
    "GATA3": {"T", "Thymocyte"},  # also Th2/ILC; keep as NK-ILC/subclass elsewhere
}

EPI_CELL = anchor("epithelial cell")

# Stems whose genes form the cross-lineage program union (cell-cycle / apoptosis
# / stress / metabolism / signaling + HALLMARK/REACTOME/KEGG).  Used to call
# `state` for genes that do not mark the panel lineage.
PROGRAM_RE = re.compile(r"^(1[7-9]|20|2[1-9]|3[0-9]|4[0-9])_")

# Cell-cycle / proliferation is the textbook cross-lineage contaminant the spec
# calls out (a proliferating cell of ANY lineage expresses MKI67/TOP2A/...), and
# CellxGene sometimes annotates such genes to a single proliferating subtype
# (e.g. CCNB2/CEP55 -> "DN4 thymocyte"), which would otherwise read as lineage
# identity.  Membership in a cell-cycle panel OR this curated marker set forces
# `state`, overriding CL -- but never a curated canonical identity marker.
CELL_CYCLE_RE = re.compile(r"cell_cycle", re.IGNORECASE)
CELL_CYCLE_MARKERS = {
    "MKI67", "TOP2A", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7",
    "CCNA2", "CCNB1", "CCNB2", "CCNE1", "CCNE2", "CCND1", "CDK1", "CDK2", "CDK4",
    "CDK6", "CDC20", "CDC45", "CDC6", "CDT1", "AURKA", "AURKB", "BUB1", "BUB1B",
    "BIRC5", "CENPA", "CENPE", "CENPF", "CEP55", "KIF11", "KIF23", "NUSAP1",
    "UBE2C", "TYMS", "RRM2", "TK1", "PLK1", "FOXM1", "E2F1", "ASPM", "NUF2",
    "HMMR", "PBK", "TPX2", "PCLAF", "GINS2", "EXO1", "ANLN", "PTTG1", "CKS1B",
    "CKS2", "NDC80", "SPC25", "DLGAP5", "CDKN3",
}


# ---------------------------------------------------------------------------
# IO helpers
# ---------------------------------------------------------------------------
def read_genes(path):
    """Ordered gene symbols from a subpanel CSV (the `gene` column)."""
    out = []
    with open(path, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            g = (r.get("gene") or "").strip()
            if g:
                out.append(g)
    return out


def load_reference():
    """gene -> (cell_type terms list, cellchat nonempty bool)."""
    ct, cc = {}, {}
    with open(REFERENCE, newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            g = r["gene"].strip()
            ct[g] = [t.strip() for t in (r.get("cell_type") or "").split(";") if t.strip()]
            cc[g] = bool((r.get("cellchat_pathway") or "").strip())
    return ct, cc


def load_detection(tissue):
    """gene -> (det_r1, det_r2, det_max) from the tissue audit (or None)."""
    path = os.path.join(TISSUE_DIR[tissue], "audit", "xenium5k_audit.csv")
    det = {}
    with open(path, newline="", encoding="utf-8") as fh:
        rd = csv.DictReader(fh)
        dcols = [c for c in rd.fieldnames if c.startswith("detection_pct_")]
        for r in rd:
            vals = []
            for c in dcols:
                try:
                    vals.append(float(r[c]))
                except (TypeError, ValueError):
                    vals.append(0.0)  # blank pct = not detected; keep run alignment
            if vals:
                det[r["gene"].strip()] = (vals[0], vals[-1], max(vals))
    return det, dcols


def panel_path(source, stem):
    if source == "shared":
        return os.path.join(SHARED, stem + ".csv")
    return os.path.join(TISSUE_DIR[source], "subpanels", stem + ".csv")


def _panel_dirs():
    return [SHARED] + [os.path.join(TISSUE_DIR[t], "subpanels") for t in TISSUE_DIR]


def build_program_union():
    """Union of gene symbols across all program/pathway panels (stems 17-49)."""
    union = set()
    for d in _panel_dirs():
        if not os.path.isdir(d):
            continue
        for f in os.listdir(d):
            if f.endswith(".csv") and PROGRAM_RE.match(f):
                union.update(read_genes(os.path.join(d, f)))
    return union


def build_cellcycle_union():
    """Genes in any cell-cycle panel, plus the curated proliferation markers."""
    union = set(CELL_CYCLE_MARKERS)
    for d in _panel_dirs():
        if not os.path.isdir(d):
            continue
        for f in os.listdir(d):
            if f.endswith(".csv") and CELL_CYCLE_RE.search(f):
                union.update(read_genes(os.path.join(d, f)))
    return union


# ---------------------------------------------------------------------------
# Core classification
# ---------------------------------------------------------------------------
def map_terms(terms):
    """cell_type term names -> (mapped CL ids in order, unmapped names)."""
    ids, unmapped = [], []
    for t in terms:
        cid = NAME2ID.get(t.lower())
        if cid:
            ids.append(cid)
        else:
            unmapped.append(t)
    return ids, unmapped


def _term_names(ids, cap=6):
    seen, out = set(), []
    for cid in ids:
        nm = ID2NAME.get(cid, cid)
        if nm not in seen:
            seen.add(nm)
            out.append(nm)
    return ";".join(out[:cap])


def classify_gene(gene, lin, ct_terms, is_program, force_state=False):
    """Return (tier, lineage_hint, subtype_hint, cl_anchor, cl_terms, evidence, note)."""
    canon = CANONICAL.get(lin.key, {})
    term_ids, _ = map_terms(ct_terms)
    supports = [t for t in term_ids if any(a in ancestors(t) for a in lin.anchors)]
    n_mapped = len(term_ids)

    # 1. Canonical override -> identity (carry subtype_hint; CL agreement -> both).
    if gene in canon:
        sub = canon[gene]
        evidence = "both" if supports else "canonical"
        cl_terms = _term_names(supports) if supports else ""
        cl_anc = lin.anchor_names[0]
        if sub and sub in lin.subtypes:
            cl_anc = ID2NAME.get(lin.subtypes[sub][0], cl_anc)
        note = "canonical identity marker" + ("" if supports else " (CL-unannotated)")
        return (IDENTITY, lin.key, sub, cl_anc, cl_terms, evidence, note)

    # 1b. Cell-cycle / proliferation program -> state, overriding any CL lineage
    # annotation (these are cross-lineage; a CellxGene mono-annotation to one
    # proliferating subtype must not read as identity).
    if force_state:
        ev = "cellxgene" if n_mapped else "none"
        return (STATE, lin.key, "", "", _term_names(term_ids),
                ev, "cell-cycle/proliferation program (cross-lineage)")

    # 2. Unannotated and not canonical -> non_specific.
    if n_mapped == 0:
        note = "no cell_type annotation; not a canonical marker"
        if is_program:
            return (STATE, lin.key, "", "", "", "none", "program gene; " + note)
        return (NONSPEC, lin.key, "", "", "", "none", note)

    # 2b. Dedicated pan-epithelial panel (16): under `epithelial cell` ->
    # epithelial_general (this panel is the shared-epithelium embedding-extra,
    # never a per-lineage mask).  Only canonical pan-epi markers (step 1) are
    # identity here; the generic support->identity path (step 3) must not fire.
    if lin.epi_general_panel:
        epi = [t for t in term_ids if EPI_CELL in ancestors(t)]
        if epi and len(epi) / n_mapped >= EPI_MIN:
            return (EPIGEN, lin.key, "", ID2NAME[EPI_CELL], _term_names(epi),
                    "cellxgene", "shared epithelium (pan-epithelial panel)")
        if is_program:
            return (STATE, lin.key, "", "", _term_names(term_ids),
                    "cellxgene", "program gene; cross-lineage")
        return (NONSPEC, lin.key, "", "", _term_names(term_ids),
                "cellxgene", "not epithelial; maps off-lineage")

    neg = lin.key in CANONICAL_NEGATIVE.get(gene, set())
    support_share = len(supports) / n_mapped

    # 3. CL dominantly supports the lineage anchor (and not a negative override).
    if support_share >= SUPPORT_MIN and not neg:
        # subclass iff supporting terms concentrate on a single subtype and none
        # is a lineage-general marker (supports the anchor but no subtype).
        sub_hits = {}            # hint -> count
        general = 0
        for t in supports:
            anc_t = ancestors(t)
            hit = [h for h, sa in lin.subtypes.items() if any(a in anc_t for a in sa)]
            if hit:
                for h in hit:
                    sub_hits[h] = sub_hits.get(h, 0) + 1
            else:
                general += 1
        if sub_hits and general == 0:
            top, n_top = max(sub_hits.items(), key=lambda kv: kv[1])
            if n_top / len(supports) >= SUBTYPE_MIN:
                cl_anc = ID2NAME.get(lin.subtypes[top][0], lin.anchor_names[0])
                return (SUBCLASS, lin.key, top, cl_anc, _term_names(supports),
                        "cellxgene", "subtype-restricted (CL)")
        return (IDENTITY, lin.key, "", lin.anchor_names[0], _term_names(supports),
                "cellxgene", "marks lineage (CL)")

    # 4. Off-lineage but under `epithelial cell` -> epithelial_general
    # (specific epithelial lineages; the pan-epithelial panel 16 returned at 2b).
    if lin.epithelial:
        epi = [t for t in term_ids if EPI_CELL in ancestors(t) and t not in supports]
        if epi and len(epi) / n_mapped >= EPI_MIN:
            return (EPIGEN, lin.key, "", ID2NAME[EPI_CELL], _term_names(epi),
                    "cellxgene", "shared epithelium; not lineage-specific")

    # 5. Program / state gene that does not mark the lineage.
    if is_program:
        return (STATE, lin.key, "", "", _term_names(term_ids),
                "cellxgene", "program gene; cross-lineage")

    # 6. Maps elsewhere / too broad -> non_specific.
    note = "negative override; not lineage identity" if neg else \
           "maps off-lineage or to broad terms"
    return (NONSPEC, lin.key, "", "", _term_names(term_ids), "cellxgene", note)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------
def classify_celltype_panels():
    """Classify every gene of every cell-type identity panel.

    Returns a dict with:
      rows          list of [source, subpanel, gene, gene_tier, lineage_hint,
                    subtype_hint, cl_anchor, cl_terms, evidence, note]
      identity_core {(source, stem): [genes]}    identity-tier genes per panel
      n_identity    {(source, stem): int}
      audit         {tissue: [rows]}             specificity check
      epi_general   {tissue: [rows]}             shared-epithelium genes
    """
    ct_terms, cellchat = load_reference()
    program = build_program_union()
    cellcycle = build_cellcycle_union()
    det = {t: load_detection(t) for t in TISSUE_DIR}
    # Per-tissue audit covers each tissue's own panels plus the shared panels it
    # actually masks with (pancreas uses the shared immune/vascular/neural set;
    # thymus appends only shared neural + epithelial-general).
    SHARED_BY_TISSUE = {
        "pancreas": {"03_immune_T_cell", "04_immune_B_plasma", "05_immune_NK_ILC",
                     "06_immune_myeloid", "07_immune_granulocyte",
                     "11_endothelial_vascular", "12_pericyte_smooth_muscle",
                     "15_neural_glial", "16_epithelial_general"},
        "thymus": {"15_neural_glial", "16_epithelial_general"},
    }

    rows, identity_core, n_identity = [], {}, {}
    epi_general = {t: [] for t in TISSUE_DIR}

    for (source, stem), lin_key in sorted(PANEL_LINEAGE.items()):
        lin = LINEAGES[lin_key]
        genes = read_genes(panel_path(source, stem))
        core = []
        for g in genes:
            is_prog = cellchat.get(g, False) or (g in program)
            tier, lh, sh, ca, clt, ev, note = classify_gene(
                g, lin, ct_terms.get(g, []), is_prog, force_state=g in cellcycle)
            rows.append([source, stem, g, tier, lh, sh, ca, clt, ev, note])
            if tier == IDENTITY:
                core.append((g, sh, ev, note))
            if tier == EPIGEN:
                for t in TISSUE_DIR:
                    if (source != "shared" and source == t) or \
                       (source == "shared" and stem in SHARED_BY_TISSUE[t]):
                        epi_general[t].append([source, stem, g, lin_key, clt])
        identity_core[(source, stem)] = core
        n_identity[(source, stem)] = len(core)

    # Specificity audit: identity_core genes vs per-tissue detection.
    audit = {t: [] for t in TISSUE_DIR}
    for (source, stem), core in identity_core.items():
        for tissue in TISSUE_DIR:
            applies = (source == tissue) or \
                      (source == "shared" and stem in SHARED_BY_TISSUE[tissue])
            if not applies:
                continue
            d, _ = det[tissue]
            for g, sh, ev, _note in core:
                rec = d.get(g)
                det_max = rec[2] if rec else None
                flag = ""
                if det_max is None:
                    flag = "not_in_audit"
                elif det_max >= UBIQUITOUS_PCT:
                    flag = "ubiquitous_candidate_non_specific"
                audit[tissue].append([
                    source, stem, g, sh, ev,
                    "" if rec is None else round(rec[0], 3),
                    "" if rec is None else round(rec[1], 3),
                    "" if det_max is None else round(det_max, 3), flag])

    return {"rows": rows, "identity_core": identity_core,
            "n_identity": n_identity, "audit": audit, "epi_general": epi_general}


# ---------------------------------------------------------------------------
# Output writers
# ---------------------------------------------------------------------------
GENE_ROLE_HEADER = ["source", "subpanel", "gene", "gene_tier", "lineage_hint",
                    "subtype_hint", "cl_anchor", "cl_terms", "evidence", "note"]


def _write(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def identity_core_path(source, stem):
    if source == "shared":
        return os.path.join(DATA, "identity_core", stem + ".csv")
    return os.path.join(TISSUE_DIR[source], "identity_core", stem + ".csv")


def write_identity_outputs(result):
    """Write identity_core/, identity_audit_<t>.csv, identity_epithelial_general_<t>.csv."""
    for (source, stem), core in result["identity_core"].items():
        _write(identity_core_path(source, stem),
               ["gene", "subtype_hint", "evidence", "note"], core)
    for tissue, audit in result["audit"].items():
        _write(os.path.join(DATA, "identity_audit_%s.csv" % tissue),
               ["source", "subpanel", "gene", "subtype_hint", "evidence",
                "detection_pct_r1", "detection_pct_r2", "det_max",
                "specificity_flag"], audit)
    for tissue, eg in result["epi_general"].items():
        _write(os.path.join(DATA, "identity_epithelial_general_%s.csv" % tissue),
               ["source", "subpanel", "gene", "lineage_hint", "cl_terms"], eg)


def write_panel_gene_roles(celltype_rows, mixed_rows, path=None):
    """Unified data/panel_gene_roles.csv: mixed-panel rows + cell-type rows."""
    path = path or os.path.join(DATA, "panel_gene_roles.csv")
    _write(path, GENE_ROLE_HEADER, list(mixed_rows) + list(celltype_rows))


def run(write_outputs=True):
    result = classify_celltype_panels()
    if write_outputs:
        write_identity_outputs(result)
    return result


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # Standalone: classify and write the identity_core / audit / epithelial-general
    # outputs, and print a summary. The unified data/panel_gene_roles.csv (which
    # also needs the hand-curated mixed-panel rows) is written by the orchestrator
    # `python3 scripts/build_panel_roles.py`; this module never imports it.
    from collections import Counter
    res = run(write_outputs=True)

    print("CL release:", cl_data_version())
    print("cell-type gene rows: %d across %d panels"
          % (len(res["rows"]), len(res["identity_core"])))
    tiers = Counter(r[3] for r in res["rows"])
    print("tier counts:", dict(tiers))
    print("\nidentity_core sizes (n_identity / n_genes):")
    for (src, stem), n in sorted(res["n_identity"].items()):
        ng = len(read_genes(panel_path(src, stem)))
        print("  %-9s %-32s %4d / %4d" % (src, stem, n, ng))
    for t in TISSUE_DIR:
        flagged = sum(1 for r in res["audit"][t] if r[-1])
        print("audit %s: %d identity rows, %d flagged; epi_general: %d"
              % (t, len(res["audit"][t]), flagged, len(res["epi_general"][t])))
