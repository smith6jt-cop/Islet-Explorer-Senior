"""Tier-3 STATE programs — score as orthogonal obs module scores, never a lineage axis.

Upstream `panel_roles.csv` tags each pathway/program panel `phenotyping_use ==
module_score` (IFN response, cytokine signaling, ECM remodeling, proliferation,
apoptosis, ER-stress, the HALLMARK/REACTOME/KEGG pathways, and the pancreas
T1D-insulin-biology panel). Per the user-directed redesign (decision #3:
"score + report + guard identity") these are scored into `obs['state_<group>_score']`
and reported, but NEVER used to seed/argmax a lineage — that is the GUARD against
an IFN-high / stressed / cycling cell being mis-lineaged.

This supersedes scripts/program_panels.py (which loaded the pancreas pathway
snapshot for the fig_program_coverage figure only); state_panels is tissue-aware
(reads the per-tissue vendored state panels via _panel_roles) and computes obs
scores. Only the `gene` column of each panel CSV is read.
"""
import os
from pathlib import Path

import pandas as pd
import scanpy as sc

from scripts.tissue_markers import _panel_roles

# Two roots, deliberately separate.  ROOT is where the CODE and the tracked tables live;
# DATA_ROOT is where the multi-GB h5ads live and is NEVER inside the repo.  Collapsing them
# -- as the single HiPerGator ROOT did -- would put data/processed/ under analysis/xenium/,
# and the repo-root .gitignore pattern "/data/" is ANCHORED, so it would not be ignored.
ROOT = Path(__file__).resolve().parents[1]                     # analysis/xenium
REPO = Path(__file__).resolve().parents[3]                     # monorepo root
DATA_ROOT = Path(os.environ.get("XENIUM_DATA_ROOT", "/data/xenium"))
DEFAULT_SUBPANEL_DIR = REPO / "data/development/xenium_panels/subpanels"

MIN_STATE_GENES = 5  # skip (report) a group with fewer panel-present genes
RANDOM_STATE = 0


def load_state_panels(tissue, reference_dir=DEFAULT_SUBPANEL_DIR):
    """{group: set(gene)} for the tissue's module_score state panels (gene column only)."""
    groups = _panel_roles.state_panels(tissue)   # {group: [local_rel]}
    root = Path(reference_dir)
    out = {}
    for group, rels in groups.items():
        genes = set()
        for rel in rels:
            df = pd.read_csv(root / rel)
            genes |= {g for g in df["gene"].dropna().astype(str) if g}
        if genes:
            out[group] = genes
    return out


def score_state_panels(adata, tissue, reference_dir=DEFAULT_SUBPANEL_DIR,
                       layer="lognorm"):
    """Compute obs['state_<group>_score'] via sc.tl.score_genes on `layer`.

    Panel-present genes only; groups with < MIN_STATE_GENES on this sample's panel
    are skipped (reported, non-fatal). Returns {group: n_genes_used}.
    """
    panels = load_state_panels(tissue, reference_dir)
    used = {}
    for group in sorted(panels):
        present = [g for g in panels[group] if g in adata.var_names]
        if len(present) < MIN_STATE_GENES:
            print(f"    [state] {group}: {len(present)} on-panel (<{MIN_STATE_GENES}) — skipped")
            continue
        sc.tl.score_genes(adata, present, score_name=f"state_{group}_score",
                          layer=layer, use_raw=False, random_state=RANDOM_STATE)
        used[group] = len(present)
    print(f"    [state] scored {len(used)} programs: "
          f"{', '.join(f'{g}({n})' for g, n in sorted(used.items()))}")
    return used
