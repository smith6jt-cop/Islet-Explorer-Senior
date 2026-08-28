"""Load the upstream IDENTITY_CORE gene sets (scANVI-redesign 2026-06-08, @076d1192).

The identity_core CSVs
(``data/reference/subpanels/identity_core/{shared,<tissue>}/<stem>.csv``; columns
``gene,subtype_hint,evidence,note``) are the CLEAN, CL-ontology-curated +
specificity-audited identity gene set per cell-type panel. They REPLACE the fat
subpanels in the broad CALL — a fat-panel mean is not an identity score, and
z-scoring it is the catch-all's silent killer (HANDOFF 2026-06-08 lessons 1+2) —
and drive fine subtyping via ``subtype_hint``.

The per-tissue file→lineage mapping reuses the tissue module's ``SUBPANEL_ALLOWLIST``
(broad_mask) + ``SUBPANEL_EMBEDDING_EXTRA`` (16_epithelial_general). Fine subtyping
uses the module's ``IC_HINT_FINE`` (``subtype_hint`` → fine label: ``beta`` → ``Beta``,
``CD4`` → ``SP4_thymocyte``, ``plasma`` → ``Plasma_cells``, …); genes whose hint is
unmapped do NOT enter any fine score panel (those fine types fall back to the tissue
module's MARKER_DICT — the plan's sanctioned thin-fine fallback).

Only the ``gene`` + ``subtype_hint`` columns are consumed. ``broad_genes`` /
``embedding_extra_genes`` raise ``FileNotFoundError`` on a missing core
(reproducibility: CLAUDE.md #6 — the vendored snapshot must be complete).
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
DEFAULT_SUBPANEL_DIR = REPO / "data/development/xenium_panels/subpanels"


def _ic_rel(allowlist_key):
    """``'shared/03_immune_T_cell.csv'`` → ``'identity_core/shared/03_immune_T_cell.csv'``;
    ``'pancreas/01_...csv'`` → ``'identity_core/pancreas/01_...csv'``.

    The allowlist key's leading dir (``shared`` / ``<tissue>``) is preserved under
    ``identity_core/`` — matching the layout fetch_subpanels.py vendors."""
    sub, stem = allowlist_key.split("/", 1)
    return f"identity_core/{sub}/{stem}"


def load_core_table(tissue, mod, reference_dir=DEFAULT_SUBPANEL_DIR):
    """Concatenated identity_core rows for all of `tissue`'s broad-mask +
    embedding-extra panels. Adds columns: ``lineage`` (None for embedding-extra),
    ``file_stem``, ``role`` ∈ {broad_mask, embedding_extra}."""
    root = Path(reference_dir)
    frames = []
    for key, lineage in mod.SUBPANEL_ALLOWLIST.items():
        df = pd.read_csv(root / _ic_rel(key))
        df = df.assign(lineage=lineage, file_stem=key.split("/", 1)[1],
                       role="broad_mask")
        frames.append(df)
    for key in getattr(mod, "SUBPANEL_EMBEDDING_EXTRA", []):
        df = pd.read_csv(root / _ic_rel(key))
        df = df.assign(lineage=None, file_stem=key.split("/", 1)[1],
                       role="embedding_extra")
        frames.append(df)
    out = pd.concat(frames, ignore_index=True)
    out["gene"] = out["gene"].astype(str)
    if "subtype_hint" not in out.columns:
        out["subtype_hint"] = ""
    out["subtype_hint"] = out["subtype_hint"].fillna("").astype(str).str.strip()
    return out


def broad_genes(tissue, mod, reference_dir=DEFAULT_SUBPANEL_DIR):
    """``{lineage: set(gene)}`` — broad-mask identity_core unioned per lineage
    (folded files, e.g. shared/05_immune_NK_ILC → Immune_T, are merged)."""
    tbl = load_core_table(tissue, mod, reference_dir)
    tbl = tbl[tbl["role"] == "broad_mask"]
    return {L: set(sub["gene"]) for L, sub in tbl.groupby("lineage")}


def embedding_extra_genes(tissue, mod, reference_dir=DEFAULT_SUBPANEL_DIR):
    """``set(gene)`` for the embedding-extra cores (16_epithelial_general) — they
    join the scANVI feature union but seed no per-lineage mask."""
    tbl = load_core_table(tissue, mod, reference_dir)
    return set(tbl[tbl["role"] == "embedding_extra"]["gene"])


def hint_fine_panels(tissue, mod, reference_dir=DEFAULT_SUBPANEL_DIR,
                     live_genes=None):
    """``{fine_type: sorted([gene])}`` from the broad-mask cores' ``subtype_hint``
    mapped through the tissue module's ``IC_HINT_FINE``. Genes with an unmapped /
    empty hint are dropped (their fine type falls back to MARKER_DICT). Restricted
    to `live_genes` when given."""
    hint_fine = getattr(mod, "IC_HINT_FINE", {})
    if not hint_fine:
        return {}
    tbl = load_core_table(tissue, mod, reference_dir)
    tbl = tbl[tbl["role"] == "broad_mask"]
    out = {}
    for r in tbl.itertuples():
        ft = hint_fine.get(r.subtype_hint)
        if not ft:
            continue
        if live_genes is not None and r.gene not in live_genes:
            continue
        out.setdefault(ft, set()).add(r.gene)
    return {ft: sorted(g) for ft, g in out.items()}
