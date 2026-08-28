"""Shared constants + path helpers for the spatial-method benchmark.

Imported by prep_benchmark_input.py, evaluate.py, run_baseline_niche.py.
Pure-stdlib + pathlib so it imports cleanly in any of the benchmark envs.
"""
from __future__ import annotations

import re
import os
from pathlib import Path

# Two roots, deliberately separate.  ROOT is where the CODE and the tracked tables live;
# DATA_ROOT is where the multi-GB h5ads live and is NEVER inside the repo.  Collapsing them
# -- as the single HiPerGator ROOT did -- would put data/processed/ under analysis/xenium/,
# and the repo-root .gitignore pattern "/data/" is ANCHORED, so it would not be ignored.
ROOT = Path(__file__).resolve().parents[2]                     # analysis/xenium
REPO = Path(__file__).resolve().parents[4]                     # monorepo root
DATA_ROOT = Path(os.environ.get("XENIUM_DATA_ROOT", "/data/xenium"))

# Sample -> tissue (drives the tissue_markers module + biological-recovery track).
SAMPLE_TISSUE = {
    "0041323": "pancreas",
    "0041326": "pancreas",
    "0041435": "thymus",
    "0041534": "thymus",
}

# Reproducibility: every stochastic step pins this.
SEED = 0
# Shared evaluation spatial graph (k-NN) and matched domain count (= stage-03
# baseline n_clusters). Methods are ALSO reported at their natural settings.
EVAL_SPATIAL_K = 10
N_DOMAINS_TARGET = 12

# Domain-producing methods (head-to-head) + the interaction method.
DOMAIN_METHODS = ["banksy", "utag", "kodama", "baseline"]
INTERACTION_METHODS = ["misty"]
ALL_METHODS = DOMAIN_METHODS + INTERACTION_METHODS

# Subsample variants. crop = compact contiguous region (preserves the micron
# neighbourhood graph); thin = random (sampling-robustness control).
HEADLINE_VARIANT = "crop50k"
SWEEP_VARIANTS = ["crop25k", "crop50k", "crop100k", "crop200k", "full"]
ALL_VARIANTS = ["crop25k", "crop50k", "thin50k", "crop100k", "crop200k", "full"]


def variant_to_n(variant: str) -> int | None:
    """crop50k -> 50000 ; thin50k -> 50000 ; full -> None (all cells)."""
    if variant == "full":
        return None
    m = re.fullmatch(r"(crop|thin)(\d+)k", variant)
    if not m:
        raise ValueError(f"unrecognised variant {variant!r}")
    return int(m.group(2)) * 1000


def variant_mode(variant: str) -> str:
    """'crop' | 'thin' | 'full'."""
    if variant == "full":
        return "full"
    return variant[:4] if variant.startswith("thin") else "crop"


# ---- paths ------------------------------------------------------------------
def bench_dir(sample: str) -> Path:
    return DATA_ROOT / "processed" / "benchmark" / sample


def results_dir(sample: str) -> Path:
    return bench_dir(sample) / "results"


def fig_dir(sample: str) -> Path:
    return DATA_ROOT / "figures" / sample / "benchmark"


def _stem(sample: str, variant: str) -> str:
    return f"{sample}_bench_{variant}"


def input_paths(sample: str, variant: str) -> dict[str, Path]:
    d = bench_dir(sample)
    s = _stem(sample, variant)
    return {
        "h5ad": d / f"{s}_input.h5ad",
        "expr_mtx": d / f"{s}_expr.mtx",
        "genes": d / f"{s}_genes.csv",
        "coords": d / f"{s}_coords.csv",
        "meta": d / f"{s}_meta.csv",
        "barcodes": d / f"{s}_barcodes.csv",
        "misty_markers": d / f"{s}_misty_markers.txt",
        "cropinfo": d / f"{s}_cropinfo.json",
    }


def labels_path(sample: str, variant: str, method: str) -> Path:
    return results_dir(sample) / f"{method}_{variant}_labels.csv"


def metrics_path(sample: str, variant: str, method: str) -> Path:
    return results_dir(sample) / f"{method}_{variant}_metrics.json"


def misty_results_dir(sample: str, variant: str) -> Path:
    return results_dir(sample) / f"misty_{variant}"
