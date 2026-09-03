"""Tier-driven config from the vendored XeniumPanelExplorer manifests.

Reads ``data/reference/subpanels/panel_roles.csv`` (one row per (tissue, subpanel):
tier, phenotyping_use, lineage_or_group, n_genes, phenotyping_source, note) and
``panel_gene_roles.csv`` (per-gene curation for the mixed panels) and exposes the
per-tissue config the phenotyper consumes:

  identity_allowlist(tissue)  {local_rel: lineage}      phenotyping_use == broad_mask
  embedding_extra(tissue)     [local_rel]               phenotyping_use == embedding_extra
  state_panels(tissue)        {group: [local_rel]}       phenotyping_use == module_score
  fine_tier_panels(tissue)    {group: [local_rel]}       phenotyping_use == fine_tier
  gene_roles()                DataFrame                  panel_gene_roles.csv

The 3-tier framework is documented upstream in docs/panel_classification.md and
docs/xenium_analysis_recommendations.md (see data/reference/subpanels/SOURCE.txt).
``assert_consistent(mod, tissue)`` fail-loud-diffs the derived identity allowlist
against a tissue module's hand-maintained SUBPANEL_ALLOWLIST so an upstream
re-pin that silently re-classifies a panel stops the run instead of re-mapping
every cell (CLAUDE.md #6 reproducibility, applied at the config layer).
"""
import os
from pathlib import Path

import pandas as pd

# Two roots, deliberately separate.  ROOT is where the CODE and the tracked tables live;
# DATA_ROOT is where the multi-GB h5ads live and is NEVER inside the repo.  Collapsing them
# -- as the single HiPerGator ROOT did -- would put data/processed/ under analysis/xenium/,
# and the repo-root .gitignore pattern "/data/" is ANCHORED, so it would not be ignored.
ROOT = Path(__file__).resolve().parents[2]                     # analysis/xenium
REPO = Path(__file__).resolve().parents[4]                     # monorepo root
DATA_ROOT = Path(os.environ.get("XENIUM_DATA_ROOT", "/data/xenium"))
SUBPANEL_DIR = REPO / "data/development/xenium_panels/subpanels"
ROLE_CSV = SUBPANEL_DIR / "panel_roles.csv"
GENE_ROLE_CSV = SUBPANEL_DIR / "panel_gene_roles.csv"

# Pathway rows all carry lineage_or_group == "Pathway"; key those state groups by
# file stem so 40 distinct pathways don't collapse into one score. Semantically
# named groups (IFN_response, Cytokine_signaling, …) keep their lineage_or_group.
_GENERIC_STATE_GROUP = "Pathway"


def _local_rel(phenotyping_source, stem):
    """Map a panel_roles `phenotyping_source` to the vendored local relative path
    (mirrors scripts/fetch_subpanels.py::_local_rel)."""
    s = str(phenotyping_source).replace("\\", "/")
    if "subpanels_shared/" in s:
        return f"shared/{stem}.csv"
    if "tissues/pancreas/" in s:
        return f"pancreas/{stem}.csv"
    if "tissues/thymus/" in s:
        return f"thymus/{stem}.csv"
    return f"misc/{stem}.csv"


def load_roles(tissue=None):
    if not ROLE_CSV.exists():
        raise FileNotFoundError(
            f"{ROLE_CSV} missing — vendor it with `python scripts/fetch_subpanels.py`")
    df = pd.read_csv(ROLE_CSV)
    if tissue is not None:
        df = df[df["tissue"] == tissue].copy()
    df["local_rel"] = [_local_rel(s, st)
                       for s, st in zip(df["phenotyping_source"], df["subpanel"])]
    return df


def _use(tissue, use):
    df = load_roles(tissue)
    return df[df["phenotyping_use"] == use]


def identity_allowlist(tissue):
    """{local_rel: lineage} for the broad-mask identity panels."""
    df = _use(tissue, "broad_mask")
    return {r.local_rel: r.lineage_or_group for r in df.itertuples()}


def embedding_extra(tissue):
    """[local_rel] for the embedding-extra panels (16_epithelial_general)."""
    return list(_use(tissue, "embedding_extra")["local_rel"])


def _grouped(df):
    out = {}
    for r in df.itertuples():
        group = r.subpanel if r.lineage_or_group == _GENERIC_STATE_GROUP \
            else r.lineage_or_group
        out.setdefault(group, []).append(r.local_rel)
    return out


def state_panels(tissue):
    """{group: [local_rel]} for the module_score (Tier-3 state) panels."""
    return _grouped(_use(tissue, "module_score"))


def fine_tier_panels(tissue):
    """{group: [local_rel]} for the fine_tier subclass panels (thymus 10/17/19/20)."""
    return _grouped(_use(tissue, "fine_tier"))


def gene_roles(source=None, subpanel=None):
    """panel_gene_roles.csv (mixed-panel per-gene curation). Optional filters."""
    if not GENE_ROLE_CSV.exists():
        raise FileNotFoundError(
            f"{GENE_ROLE_CSV} missing — vendor it with `python scripts/fetch_subpanels.py`")
    df = pd.read_csv(GENE_ROLE_CSV)
    if source is not None:
        df = df[df["source"] == source]
    if subpanel is not None:
        df = df[df["subpanel"] == subpanel]
    return df


def panel_genes(tissue, uses, reference_dir=SUBPANEL_DIR):
    """Union of the `gene` column over all vendored panels for `tissue` whose
    phenotyping_use is in `uses` (e.g. {'module_score','fine_tier'} for the
    state-side genes used by the identity GUARD)."""
    df = load_roles(tissue)
    df = df[df["phenotyping_use"].isin(set(uses))]
    root = Path(reference_dir)
    genes = set()
    for rel in df["local_rel"]:
        p = root / rel
        if p.exists():
            col = pd.read_csv(p)["gene"].dropna().astype(str)
            genes |= {g for g in col if g}
    return genes


def promoted_identity_genes(tissue):
    """Genes curated as `identity` in a NON-mask (state/fine_tier) panel that
    should be promoted into a broad mask — {gene: lineage}. Currently the thymus
    mTEC TFs AIRE/FEZF2 (upstream Rec 3): in 11_selection (state) + 19_AIRE_TRA
    (fine_tier), gene_tier == identity, lineage_hint == mTEC → promote to TEC."""
    gr = gene_roles()
    gr = gr[(gr["gene_tier"] == "identity")]
    # Only consider rows for THIS tissue's panels (source == tissue) — the shared
    # 08-panel identity genes (KIR/KLR) are NK and already in the 05 mask.
    gr = gr[gr["source"] == tissue]
    hint_to_lineage = {"mTEC": "TEC"}
    out = {}
    for r in gr.itertuples():
        lineage = hint_to_lineage.get(str(r.lineage_hint))
        if lineage:
            out[r.gene] = lineage
    return out


def assert_consistent(mod, tissue):
    """Fail-loud if the vendored panel_roles.csv disagrees with the tissue module's
    hand-maintained SUBPANEL_ALLOWLIST / SUBPANEL_EMBEDDING_EXTRA (drift guard)."""
    derived = identity_allowlist(tissue)
    hand = dict(getattr(mod, "SUBPANEL_ALLOWLIST", {}))
    if derived != hand:
        only_csv = {k: v for k, v in derived.items() if hand.get(k) != v}
        only_mod = {k: v for k, v in hand.items() if derived.get(k) != v}
        raise SystemExit(
            f"[panel_roles] SUBPANEL_ALLOWLIST drift for {tissue}:\n"
            f"  in panel_roles.csv but not matching module: {only_csv}\n"
            f"  in module but not matching panel_roles.csv: {only_mod}\n"
            f"  Reconcile scripts/tissue_markers/{tissue}.py with the vendored "
            f"panel_roles.csv (or re-pin/re-fetch).")
    d_extra = set(embedding_extra(tissue))
    h_extra = set(getattr(mod, "SUBPANEL_EMBEDDING_EXTRA", []))
    if d_extra != h_extra:
        raise SystemExit(
            f"[panel_roles] SUBPANEL_EMBEDDING_EXTRA drift for {tissue}: "
            f"panel_roles={sorted(d_extra)} vs module={sorted(h_extra)}")
    print(f"  [panel_roles] {tissue}: allowlist + embedding_extra consistent "
          f"({len(derived)} identity masks, {len(d_extra)} embedding-extra)")
