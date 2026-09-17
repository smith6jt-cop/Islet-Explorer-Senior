"""Arm X1 ``scanvi-per-donor`` -- the incumbent per-donor scANVI typer, read from its sidecar.

Nothing is re-run here.  The sidecar ``<serial>_hierarchical.h5ad`` (analysis/xenium
``hierarchical_phenotype.py --output slim``, decisions #16/#17) already carries
``celltype_broad`` (8 lineages), ``celltype`` (18 fine), ``broad_confidence``, ``broad_lowsignal``,
``cycling`` and the ``state_*_score`` columns; this module writes them in the lane's one record
shape so the same gate (``gate.judge``) judges X1 beside X2 and X3.

Honest limits, recorded rather than smoothed over: the typer is stochastic (scVI/scANVI) and was run
once per donor with ``random_state 0``, so it supplies NO replicates -- the gate's stability
criterion is unjudged for X1, and unjudgeable is never a pass.  ``reason`` is ``lowsignal`` where the
typer flagged the cell (no eligible lineage at the floor; the label is its absolute-score argmax).

    cd core/phenocycler && ../../.venv/bin/python scripts/xenium_typing/nullgate_typer.py --donors 6436 6442 [--judge]
"""
from __future__ import annotations

import argparse
import pathlib
import sys
import time

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import gate  # noqa: E402
import pins  # noqa: E402

ARM = pins.ARM_SCANVI


def read_sidecar(donor: str) -> pd.DataFrame:
    import anndata as ad
    p = pins.sidecar_path(donor)
    if not p.exists():
        raise FileNotFoundError(f"donor {donor}: no typer sidecar at {p}")
    a = ad.read_h5ad(p)
    obs = a.obs
    out = pd.DataFrame({
        "cell_id": np.asarray(a.obs_names).astype(str),
        "lineage_native": obs["celltype_broad"].astype(str).to_numpy(),
        "celltype_fine": obs["celltype"].astype(str).to_numpy(),
        "confidence": obs["broad_confidence"].astype(float).to_numpy(),
        "reason": np.where(obs["broad_lowsignal"].astype(bool).to_numpy(), "lowsignal", ""),
        "cycling": obs["cycling"].astype(bool).to_numpy(),
    })
    for c in obs.columns:
        if c.startswith("state_") and c.endswith("_score"):
            out[c] = obs[c].astype(float).to_numpy()
    return out


def run(donor: str, replicates: pd.DataFrame | None = None) -> pathlib.Path:
    """``replicates``: cell_id + one ``seed:<k>`` label column per re-run of the typer at another
    RANDOM_STATE (``scanvi_replicates.py``); None = the incumbent single run, stability unjudged."""
    t0 = time.time()
    calls = read_sidecar(donor)
    cores = pins.identity_cores()
    rec = gate.ArmRecord(
        arm=ARM, donor_id=str(donor), config_id="primary",
        identity_inputs={L: len(g) for L, g in cores.items()},
        n_cells=int(len(calls)), n_called=int(calls["lineage_native"].isin(pins.LINEAGES).sum()),
        n_ambiguous=0, n_unavailable=0, seconds=time.time() - t0, env=sys.executable, stochastic=True,
        pinned={"sidecar": str(pins.sidecar_path(donor)), "sidecar_sha256": pins.sha256_file(pins.sidecar_path(donor)),
                "base": str(pins.base_path(donor)), "panel": pins.panel_snapshot()},
        params={"source": "analysis/xenium/scripts/hierarchical_phenotype.py --output slim (run 2026-09-03)",
                "replicates": ("none: one run at random_state 0; stability is unjudged" if replicates is None else
                               f"{len(replicates.columns) - 1} re-runs of the same typer at RANDOM_STATE 1.. (scanvi_replicates.py); "
                               "primary sidecar untouched")},
        notes="incumbent per-donor scANVI typer, read from its sidecar; nothing re-run")
    return gate.write_calls(calls, rec, replicates=replicates)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="+", default=["6436", "6442"])
    ap.add_argument("--judge", action="store_true")
    a = ap.parse_args(argv)
    for d in a.donors:
        p = run(d)
        print(f"{ARM} {d}: wrote {p}")
        if a.judge:
            gate.judge(ARM, d)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
