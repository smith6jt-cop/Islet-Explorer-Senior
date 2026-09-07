#!/home/smith6jt/miniconda3/envs/scvi-env/bin/python
"""Seed replicates for arm X1 ``scanvi-per-donor`` so its stability criterion can be judged.

The incumbent typer (``analysis/xenium/scripts/hierarchical_phenotype.py``) hard-codes ``RANDOM_STATE = 0``
and writes its slim sidecar beside the base.  This runner imports that module unchanged, sets its
``RANDOM_STATE`` to the replicate seed, and redirects ONLY the sidecar write to
``arms/scanvi-per-donor/donor_id=<d>/replicates/<serial>_hierarchical__seed<k>.h5ad`` (the same slim
columns), so the primary sidecar at ``/data/xenium/processed`` is never touched and the Xenium lane's
code is not edited.  Every other setting is the typer's own (identity cores, seeds, kNN, state panels).
Afterwards ``replicate_labels.parquet`` (cell_id + ``seed:<k>`` = ``celltype_broad``) is written into the
arm directory through ``nullgate_typer.run(..., replicates=)`` and the gate re-judges.

    cd core/phenocycler && scripts/xenium_typing/scanvi_replicates.py --donors 6436 6442 --seeds 1 2 3 4 --judge
"""
from __future__ import annotations

import argparse
import importlib
import pathlib
import sys
import time

import numpy as np
import pandas as pd

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pins  # noqa: E402

XA = pins.XENIUM_ANALYSIS
ARM = pins.ARM_SCANVI


def replicate_dir(donor: str) -> pathlib.Path:
    return pins.arm_dir(ARM, donor) / "replicates"


def replicate_path(donor: str, seed: int) -> pathlib.Path:
    return replicate_dir(donor) / f"{pins.serial(donor)}_hierarchical__seed{seed}.h5ad"


def run_seed(donor: str, seed: int, *, dl_workers: int = 0) -> pathlib.Path:
    if str(XA) not in sys.path:
        sys.path.insert(0, str(XA))
    hp = importlib.import_module("scripts.hierarchical_phenotype")
    mod = importlib.import_module("scripts.tissue_markers.pancreas")
    out = replicate_path(donor, seed)
    out.parent.mkdir(parents=True, exist_ok=True)
    if out.exists():
        print(f"[replicate] {donor} seed {seed}: exists, skipping {out}", flush=True)
        return out

    def write_output(adata, sample, mode, input_path=None):
        import anndata as ad
        keep = [c for c in hp._SLIM_OBS if c in adata.obs.columns]
        keep += [c for c in adata.obs.columns if c.startswith("leiden_sub_") or c.startswith("broadscore_") or c.endswith("_score")]
        slim = ad.AnnData(obs=adata.obs[keep].copy())
        slim.obs_names = adata.obs_names
        for k in ("spatial", "X_scvi_broad"):
            if k in adata.obsm:
                slim.obsm[k] = np.asarray(adata.obsm[k])
        slim.uns["replicate"] = {"arm": ARM, "donor_id": donor, "seed": int(seed), "random_state": int(hp.RANDOM_STATE)}
        slim.write_h5ad(out)
        print(f"[replicate] {donor} seed {seed}: wrote {out} {slim.shape}", flush=True)
        return out

    hp.RANDOM_STATE = int(seed)
    hp.DL_WORKERS = int(dl_workers)
    if not hasattr(hp, "LOAD_SPARSE"):
        hp.LOAD_SPARSE = True
    hp.write_output = write_output
    use_gpu = hp.setup_torch_tf32()
    rsc = hp.try_import_rsc()
    t0 = time.time()
    hp.process_sample(pins.serial(donor), mod, pins.SUBPANEL_DIR, use_gpu, rsc, "slim", False, input_path=pins.base_path(donor))
    print(f"[replicate] {donor} seed {seed}: {(time.time() - t0) / 60:.1f} min", flush=True)
    return out


def replicate_table(donor: str, seeds: list, cell_ids: np.ndarray) -> pd.DataFrame:
    import anndata as ad
    rep = pd.DataFrame({"cell_id": cell_ids.astype(str)})
    for s in seeds:
        a = ad.read_h5ad(replicate_path(donor, s), backed="r")
        lab = pd.Series(a.obs["celltype_broad"].astype(str).to_numpy(), index=np.asarray(a.obs_names).astype(str))
        col = lab.reindex(rep["cell_id"])
        if col.isna().any():
            raise RuntimeError(f"{donor} seed {s}: {int(col.isna().sum())} primary cells absent from the replicate sidecar")
        rep[f"seed:{s}"] = col.to_numpy()
    return rep


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--donors", nargs="+", default=["6436", "6442"])
    ap.add_argument("--seeds", nargs="+", type=int, default=[1, 2, 3, 4])
    ap.add_argument("--dl-workers", type=int, default=0)
    ap.add_argument("--judge", action="store_true")
    ap.add_argument("--no-run", action="store_true", help="only assemble replicate_labels from existing sidecars")
    a = ap.parse_args(argv)
    import gate
    import nullgate_typer
    for d in a.donors:
        if not a.no_run:
            for s in a.seeds:
                run_seed(d, s, dl_workers=a.dl_workers)
        calls = nullgate_typer.read_sidecar(d)
        rep = replicate_table(d, a.seeds, calls["cell_id"].to_numpy())
        p = nullgate_typer.run(d, replicates=rep)
        print(f"[replicate] {d}: record rewritten with replicates {list(rep.columns[1:])} -> {p}", flush=True)
        if a.judge:
            gate.judge(ARM, d)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
