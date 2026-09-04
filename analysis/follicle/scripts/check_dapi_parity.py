"""Verify the QuPath export and the Python inference path agree to the byte.

The DAPI pipeline implements the same pixel transform twice — once in Groovy
(`export_DAPI_training_pairs.groovy`, which writes the training PNGs) and once
in Python (`DAPI_unet_utils`, which normalises slides at inference). If they
ever drift, the model trains on one intensity scale and predicts on another,
and nothing else in the pipeline would notice.

This script re-derives each exported PNG from the raw OME-TIFF using the Python
half of the contract and asserts byte-equality. It needs no QuPath and no GPU,
so it can run on any machine that has the images and `data_dapi/`.

  python scripts/check_dapi_parity.py --data-dir data_dapi --images-dir Images

Two levels of check per donor:

  * `--use-sidecar` (default) applies the lo/hi the EXPORT recorded. This
    isolates `dapi_to_uint8` — it answers "does Python's normalisation
    arithmetic match Groovy's?".
  * `--recompute` additionally re-derives otsu/lo/hi from the raw pixels. This
    exercises the whole chain (`hist256` -> `otsu3_from_hist256` ->
    `percentile_from_hist256`) and is the check that actually matters for
    inference, where no sidecar exists.

WHY THE PIXEL CHECK IS TOLERANCE-BASED (and the parameter check is not)
-----------------------------------------------------------------------
The two sides reach the 2 µm/px grid by different routes. QuPath resamples with
its own filter at a possibly non-integer downsample — for the 0.5078 µm/px
Phenocycler slides the target downsample is 3.94, so it reads a pyramid level
and rescales — while Python reads the nearest suitable level and resizes with
PIL. Both are correct; they are not bit-identical. Raw values can therefore
differ by a few counts, which is invisible after normalisation (a ~50 000-count
window maps to 256 levels, so ~200 raw counts per output level) but not
guaranteed to be zero.

So:
  * `otsu_bin` / `lo` / `hi` are compared EXACTLY. They are integer bin indices
    and bin low edges, robust to resampling noise, and a difference there means
    the arithmetic genuinely drifted.
  * The uint8 images are compared with a tolerance (`--max-abs-diff`,
    `--max-mismatch-frac`). A handful of ±1 pixels is resampling; a large
    fraction, or a difference bigger than a grey level or two, is drift.

The transform arithmetic itself is pinned exactly, without any resampling in the
way, by the imperative Groovy transcriptions in tests/test_DAPI_unet_utils.py.

Exit status is non-zero if any donor exceeds tolerance, so it can gate a run.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image

import DAPI_unet_utils as D
import HE_unet_utils as U


def find_source(images_dir: Path, base: str, source_name: str) -> Path | None:
    """Locate the OME-TIFF a sidecar came from.

    `source_name` is QuPath's image name, which may carry a
    " - resolution #1" suffix and a .qptiff extension, so match on the base name
    first and fall back to a prefix search.
    """
    stem = source_name.split(" ")[0] if source_name else base
    for cand in (stem, base):
        for ext in (".ome.tiff", ".ome.tif", ".qptiff", ".tiff", ".tif"):
            p = images_dir / f"{cand}{ext}"
            if p.exists():
                return p
        direct = images_dir / cand
        if direct.exists() and direct.is_file():
            return direct
    matches = sorted(p for p in images_dir.rglob(f"{base}*")
                     if p.is_file() and p.suffix.lower() in
                     (".tif", ".tiff", ".qptiff"))
    return matches[0] if matches else None


def _diff_stats(a: np.ndarray, b: np.ndarray) -> dict:
    d = np.abs(a.astype(np.int16) - b.astype(np.int16))
    nz = int((d > 0).sum())
    return {"mismatched_px": nz,
            "mismatch_frac": float(nz / d.size),
            "max_abs_diff": int(d.max()) if d.size else 0,
            "p999_abs_diff": int(np.percentile(d, 99.9)) if d.size else 0}


def check_donor(sample: U.DonorSample, images_dir: Path, recompute: bool,
                max_abs_diff: int, max_mismatch_frac: float) -> dict:
    meta = json.loads(Path(sample.meta_path).read_text())
    base = meta.get("image", sample.name)
    src = find_source(images_dir, base, meta.get("source_name", ""))
    if src is None:
        return {"name": sample.name, "status": "source_not_found", "detail": base}

    target_px = float(meta.get("target_px_size_um", D.TARGET_UM_PER_PX))
    raw, eff_px = D.read_dapi_channel(src, target_px)
    exported = np.array(Image.open(sample.image_path).convert("L"), dtype=np.uint8)

    if raw.shape != exported.shape:
        return {"name": sample.name, "status": "shape_mismatch",
                "detail": f"read {raw.shape} vs exported {exported.shape}",
                "source": str(src)}

    result: dict = {"name": sample.name, "source": str(src),
                    "shape": list(raw.shape), "eff_px_size_um": eff_px}

    t = D.sidecar_transform(meta)
    sidecar_stats = _diff_stats(D.dapi_to_uint8(raw, int(t["lo"]), int(t["hi"])),
                                exported)
    result["sidecar"] = sidecar_stats

    # Parameters are compared EXACTLY — they are integer bin indices / bin low
    # edges, so resampling noise does not move them and a difference is drift.
    params_match = True
    if recompute:
        p = D.dapi_transform_params(raw)
        keys = ("otsu_bin", "lo", "hi")
        result["params"] = {"recomputed": {k: int(p[k]) for k in keys},
                            "sidecar": {k: int(t[k]) for k in keys}}
        params_match = all(int(p[k]) == int(t[k]) for k in keys)
        result["recomputed"] = _diff_stats(
            D.dapi_to_uint8(raw, p["lo"], p["hi"]), exported)
    result["params_match"] = params_match

    def within(s: dict) -> bool:
        return (s["max_abs_diff"] <= max_abs_diff
                and s["mismatch_frac"] <= max_mismatch_frac)

    pixels_ok = within(sidecar_stats) and (
        not recompute or within(result["recomputed"]))
    result["status"] = "ok" if (pixels_ok and params_match) else "MISMATCH"
    return result


def check_mask_region_agreement(sample: U.DonorSample, data_dir: Path) -> dict | None:
    """Confirm `masks == (regions == Follicle)`, the exporter's own invariant."""
    reg_path = data_dir / "regions" / sample.image_path.name
    if not reg_path.exists():
        return None
    reg = np.array(Image.open(reg_path).convert("L"), dtype=np.uint8)
    msk = np.array(Image.open(sample.mask_path).convert("L"), dtype=np.uint8)
    if reg.shape != msk.shape:
        return {"name": sample.name, "status": "shape_mismatch"}
    agree = bool(np.array_equal(msk > 127, reg == D.FOLLICLE_INDEX))
    return {"name": sample.name, "status": "ok" if agree else "MISMATCH",
            "disagreeing_px": int(((msk > 127) != (reg == D.FOLLICLE_INDEX)).sum())}


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", default=str(D.DATA_DIR))
    p.add_argument("--images-dir", default=str(U.PROJECT / "Images"))
    p.add_argument("--recompute", action="store_true", default=True,
                   help="also re-derive otsu/lo/hi from raw pixels (default on)")
    p.add_argument("--use-sidecar", dest="recompute", action="store_false",
                   help="only check dapi_to_uint8 against the recorded lo/hi")
    p.add_argument("--max-abs-diff", type=int, default=2,
                   help="largest tolerated per-pixel grey-level difference "
                        "(default 2; resampling noise, not drift)")
    p.add_argument("--max-mismatch-frac", type=float, default=0.02,
                   help="largest tolerated fraction of differing pixels "
                        "(default 0.02)")
    p.add_argument("--strict", action="store_true",
                   help="require byte-exact images (tolerances 0). Correct for "
                        "synthetic data, where both sides read the same raster "
                        "with no resampling in between.")
    args = p.parse_args()
    if args.strict:
        args.max_abs_diff = 0
        args.max_mismatch_frac = 0.0

    data_dir = Path(args.data_dir)
    images_dir = Path(args.images_dir)
    samples = U.discover_samples(data_dir)
    if not samples:
        print(f"No exported donors under {data_dir}. Run "
              "export_DAPI_training_pairs.groovy first.")
        return 2
    if not images_dir.exists():
        print(f"Images dir {images_dir} does not exist.")
        return 2

    failures = 0
    tol = ("byte-exact" if args.strict else
           f"<= {args.max_abs_diff} grey levels on <= "
           f"{args.max_mismatch_frac:.1%} of pixels")
    print(f"Checking {len(samples)} donors against {images_dir}")
    print(f"Image tolerance: {tol}; otsu/lo/hi compared exactly.\n")
    for s in samples:
        r = check_donor(s, images_dir, args.recompute,
                        args.max_abs_diff, args.max_mismatch_frac)
        if r["status"] == "source_not_found":
            failures += 1
            print(f"  NO SOURCE {r['name']} ({r['detail']})")
            continue
        if r["status"] == "shape_mismatch":
            failures += 1
            print(f"  SHAPE     {r['name']}: {r['detail']}")
            continue
        sc = r["sidecar"]
        detail = (f"{sc['mismatched_px']} px differ "
                  f"({sc['mismatch_frac']:.4%}, max {sc['max_abs_diff']})")
        if r["status"] == "ok":
            print(f"  OK       {r['name']}  {r['shape'][1]}x{r['shape'][0]} "
                  f"@ {r['eff_px_size_um']} µm/px  [{detail}]")
        else:
            failures += 1
            print(f"  MISMATCH {r['name']}  [{detail}]")
            if not r["params_match"]:
                print(f"    transform parameters DIFFER: {r.get('params')}")
                print("    ^ this is real drift, not resampling — reconcile "
                      "the Groovy and Python transforms.")
        a = check_mask_region_agreement(s, data_dir)
        if a and a["status"] != "ok":
            failures += 1
            print(f"  MASK!=REGION {a['name']}: {a.get('disagreeing_px')} px "
                  "disagree — masks/ should be derived from regions/")

    print()
    if failures:
        print(f"{failures} donor(s) outside tolerance. If the transform "
              "PARAMETERS differ, the Groovy export and DAPI_unet_utils.py have "
              "drifted; reconcile them before training. If only pixels differ, "
              "widen --max-abs-diff only after confirming the differences are "
              "±1-2 grey levels scattered across the image, not systematic.")
        return 1
    print(f"All {len(samples)} donors within tolerance.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
