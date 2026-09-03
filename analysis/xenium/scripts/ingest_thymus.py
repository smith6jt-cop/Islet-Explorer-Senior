"""Ingest the two thymus Xenium bundles into data/raw/{sample}.zarr.

Mirrors the logic of notebooks/00_ingest_xenium.ipynb cell 1 but skips
nbconvert overhead. Runs both samples sequentially.
"""
import shutil
import os
from pathlib import Path

import spatialdata as sd
from spatialdata_io import xenium

# Two roots, deliberately separate.  ROOT is where the CODE and the tracked tables live;
# DATA_ROOT is where the multi-GB h5ads live and is NEVER inside the repo.  Collapsing them
# -- as the single HiPerGator ROOT did -- would put data/processed/ under analysis/xenium/,
# and the repo-root .gitignore pattern "/data/" is ANCHORED, so it would not be ignored.
ROOT = Path(__file__).resolve().parents[1]                     # analysis/xenium
REPO = Path(__file__).resolve().parents[3]                     # monorepo root
DATA_ROOT = Path(os.environ.get("XENIUM_DATA_ROOT", "/data/xenium"))
RAW_DIR = DATA_ROOT / "raw"
RAW_DIR.mkdir(parents=True, exist_ok=True)

BUNDLES = {
    "0041435": DATA_ROOT / "output-XETG00298__0041435__TMA172_THYMUS__20250606__213055",
    "0041534": DATA_ROOT / "output-XETG00298__0041534__HDL073__20250606__213055",
}

XENIUM_KWARGS = dict(
    cells_boundaries=True,
    nucleus_boundaries=True,
    cells_labels=True,
    nucleus_labels=True,
    transcripts=True,
    morphology_mip=False,
    morphology_focus=False,
    aligned_images=False,
    cells_table=True,
    n_jobs=4,
)


def ingest_one(sample_id, bundle_dir):
    zarr_path = RAW_DIR / f"{sample_id}.zarr"
    print(f"\n=== {sample_id} ===", flush=True)
    print(f"bundle: {bundle_dir}", flush=True)
    print(f"zarr:   {zarr_path}", flush=True)
    if not bundle_dir.exists():
        raise FileNotFoundError(f"Bundle directory not found: {bundle_dir}")
    if zarr_path.exists():
        try:
            sdata = sd.read_zarr(zarr_path)
            print(f"zarr already valid ({list(sdata.tables.keys())}); skipping re-import.",
                  flush=True)
            return sdata
        except Exception as err:
            print(f"existing zarr unreadable ({err!r}); recreating ...", flush=True)
            shutil.rmtree(zarr_path)
    print("running spatialdata_io.xenium(...) ...", flush=True)
    sdata = xenium(bundle_dir, **XENIUM_KWARGS)
    sdata.write(zarr_path)
    print(f"wrote {zarr_path}", flush=True)
    adata = sdata.tables["table"]
    print(f"table shape: {adata.shape}", flush=True)
    return sdata


if __name__ == "__main__":
    for sid, bdir in BUNDLES.items():
        ingest_one(sid, bdir)
    print("\nDONE", flush=True)
