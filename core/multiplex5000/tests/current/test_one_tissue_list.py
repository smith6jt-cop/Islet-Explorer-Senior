"""One authority for tissue names: the tissue profiles, served by ``integration/tissues.py``.

    .venv/bin/python -m pytest core/multiplex5000/tests/current/test_one_tissue_list.py -q

Measured 2026-10-19: the tissue names were written out in many places, each a little different
(``panc`` / ``pancreas``, ``pln`` / ``pLN`` / ``pancreatic_lymph_node``, a ``pln_1`` here and a
``lymphnode`` there). Since D2 (2026-10-20) ``profiles/tissues/*.yaml`` is the authority and
``integration/tissues.py`` reads it. A module-level constant that lists two or more tissue
spellings, or binds a ``*TISSUE*`` name to a tissue id, is a second list; it is named in KEPT with
the reason it stays and the stage that takes it away, or it fails here.

Tests are not scanned (a fixture may spell a tissue any way it likes), nor ``core/panels`` (the
panel authority keeps its own tissue folders), nor the Xenium lane's tissue plugins, nor the
recovered files (their bytes are recorded).
"""

from __future__ import annotations

import ast
import subprocess
from pathlib import Path

from multiplex5000.integration import tissues

REPO = Path(__file__).resolve().parents[4]
REGISTRY = "core/multiplex5000/multiplex5000/integration/tissues.py"
SKIPPED = ("analysis/recovered_20261008/", "core/panels/", "analysis/xenium/scripts/tissue_markers/")

#: Every spelling of every profiled tissue, and the two tissues the panel authority also covers.
SPELLINGS = {s for t, prof in tissues.PROFILES.items() for s in (t, *prof.aliases)} | {
    "spleen", "thymus"}

XENIUM = ("the Xenium lane runs in scvi-env, which cannot import multiplex5000; its sample->tissue "
          "maps are acquisition identity and move to the acquisition table at D8")
SAMESLIDE = ("the HPG same-slide lane's own labels, which its outputs carry; changed with that "
             "session at D10")

#: Tissue lists that stay for now, and why.
KEPT = {
    ("analysis/xenium/scripts/benchmark/bench_common.py", "SAMPLE_TISSUE"): XENIUM,
    ("analysis/xenium/scripts/refresh_broad_umap.py", "DEFAULT_TISSUE"): XENIUM,
    ("analysis/xenium/scripts/render_stage02_phenotyping.py", "DEFAULT_TISSUE"): XENIUM,
    ("analysis/xenium/scripts/nucleus_from_bundle.py", "PANC_TOKENS"): XENIUM,
    ("analysis/xenium/scripts/nucleus_from_bundle.py", "PLN_TOKENS"): XENIUM,
    ("analysis/xenium/scripts/fetch_subpanels.py", "TISSUES"):
        "the tissue folders core/panels ships (pancreas, thymus), which that authority owns",
    ("core/multiplex5000/multiplex5000/sections.py", "REGION_TO_ROI"):
        "section keys read a tissue token out of the image name; recorded identities replace that "
        "parsing at D4",
    ("core/multiplex5000/scripts/sameslide/compare_sections.py", "IDENTITY"): SAMESLIDE,
    ("core/multiplex5000/scripts/sameslide/completion_four.py", "TISSUE_ONLY"): SAMESLIDE,
    ("core/multiplex5000/scripts/sameslide/figures_four.py", "COVER"): SAMESLIDE,
    ("core/multiplex5000/scripts/sameslide/join_four.py", "TISSUE"): SAMESLIDE,
    ("core/multiplex5000/scripts/typing/plot_protein_literature_endpoints.py", "COLORS"):
        "a figure palette keyed by its input's own labels; palettes move to the figure stages (D11)",
    ("scripts/senior/xenium_release.py", "TISSUE_ALIASES"):
        "the Xenium release reader's own tokens (panc/pln) and two spellings the profiles do not "
        "carry (pln_1, lymphnode); unified with the RNA releases at D8",
}


def _norm(value: str) -> str:
    return value.strip().lower().replace(" ", "_").replace("-", "_")


def _strings(node: ast.AST) -> list[str]:
    if isinstance(node, (ast.Tuple, ast.List, ast.Set)):
        items = node.elts
    elif isinstance(node, ast.Dict):
        items = [*node.keys, *node.values]
    elif isinstance(node, ast.Call) and node.args:
        return _strings(node.args[0])
    else:
        return []
    return [i.value for i in items if isinstance(i, ast.Constant) and isinstance(i.value, str)]


def _lists_in(tree: ast.Module) -> list[str]:
    """Module-level names bound to a tissue list, or to a tissue id under a *TISSUE* name."""
    hits = []
    for node in tree.body:
        if not isinstance(node, (ast.Assign, ast.AnnAssign)) or node.value is None:
            continue
        targets = node.targets if isinstance(node, ast.Assign) else [node.target]
        names = [t.id for t in targets if isinstance(t, ast.Name)]
        if not names:
            continue
        spelled = {_norm(s) for s in _strings(node.value)} & SPELLINGS
        restated = (isinstance(node.value, ast.Constant) and isinstance(node.value.value, str)
                    and _norm(node.value.value) in tissues.TISSUES
                    and any("TISSUE" in n.upper() for n in names))
        if len(spelled) >= 2 or restated:
            hits.append(names[0])
    return hits


def _tissue_lists() -> set[tuple[str, str]]:
    files = subprocess.run(["git", "ls-files", "*.py"], cwd=REPO, capture_output=True, text=True,
                           check=True).stdout.split()
    found = set()
    for rel in files:
        if ("/tests/" in rel or rel.startswith("tests/") or rel.startswith(SKIPPED)
                or rel == REGISTRY):
            continue
        try:
            tree = ast.parse((REPO / rel).read_text())
        except (SyntaxError, UnicodeDecodeError, FileNotFoundError):
            continue
        found.update((rel, name) for name in _lists_in(tree))
    return found


def test_no_tissue_list_outside_the_profiles_but_the_named_ones():
    found = _tissue_lists()
    unnamed = found - set(KEPT)
    assert not unnamed, (
        f"a second list of tissue names: {sorted(unnamed)}. Read them from "
        f"multiplex5000.integration.tissues (its PROFILES, TISSUES, PANCREAS, LYMPH_NODE, "
        f"normalize), or name the list in KEPT with the reason it stays and when it goes")
    vanished = set(KEPT) - found
    assert not vanished, f"named tissue lists that are gone; take them off KEPT: {sorted(vanished)}"


def test_the_detector_sees_a_list_and_a_restated_name():
    """The scan is only as good as what it can see; hold it to both shapes it looks for."""
    tree = ast.parse('TISSUES = ("pancreas", "pln")\n'
                     'ONE = ("pancreas",)\n'
                     'PLN_TISSUE = "pancreatic_lymph_node"\n'
                     'ROI = "pancreas"\n'
                     'BY_TISSUE: dict = {"Pancreas": 1, "lymph node": 2}\n')
    assert _lists_in(tree) == ["TISSUES", "PLN_TISSUE", "BY_TISSUE"]
