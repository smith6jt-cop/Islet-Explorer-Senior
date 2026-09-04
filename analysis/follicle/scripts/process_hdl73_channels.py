#!/usr/bin/env python3
"""
Signal isolation for remaining HDL73 spleen channels.

Processes all unprocessed channels from HDL73_SPL_Registered/ by subtracting
autofluorescence using matched blank pairs, then saves results to
HDL73_SPL_Processed/ImageJ/.

Usage:
    conda run -n KINTSUGI python scripts/process_hdl73_channels.py [--force]

Flags:
    --force    Re-process channels that already exist in ImageJ/
"""

import argparse
import shutil
import sys
from pathlib import Path

import numpy as np
import tifffile

from kintsugi.signal import (
    analyze_for_subtraction,
    compute_subtraction_quality,
    subtract_autofluorescence,
)

# ── Paths ──────────────────────────────────────────────────────────────────────
CODE_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(CODE_ROOT))
from paths import DATA_ROOT, OUTPUT_ROOT  # noqa: E402

REGISTERED_DIR = DATA_ROOT / "FromHipergator" / "HDL73_SPL_Registered"
OUTPUT_DIR = OUTPUT_ROOT / "FromHipergator" / "HDL73_SPL_Processed" / "ImageJ"
PARAMS_DIR = OUTPUT_ROOT / "FromHipergator" / "HDL73_SPL_Processed" / "Processing_parameters"

# ── Channel → blank-position mapping ──────────────────────────────────────────
# Derived from channelnames.txt cycle layout (positions a/b/c per cycle)
POSITION_A = [
    "CD20", "CD31", "CD34", "CD35", "Lyve1", "PanCK", "SMActin",
]
POSITION_B = [
    "CD8", "CD15", "CD21", "CD44", "CD45RO", "CD5", "CollagenIV",
    "ECAD", "FoxP3", "Ki67", "Podoplanin",
]
POSITION_C = [
    "CD3e", "CD4", "CD11c", "CD107a", "CD163", "CD1c", "CD45",
    "CD68", "HLADR", "Vimentin",
]

CHANNEL_TO_POSITION: dict[str, str] = {}
for ch in POSITION_A:
    CHANNEL_TO_POSITION[ch] = "a"
for ch in POSITION_B:
    CHANNEL_TO_POSITION[ch] = "b"
for ch in POSITION_C:
    CHANNEL_TO_POSITION[ch] = "c"

ALL_SIGNAL_CHANNELS = sorted(CHANNEL_TO_POSITION.keys())


def load_image(path: Path) -> np.ndarray:
    """Load a single-channel TIF as uint16."""
    img = tifffile.imread(str(path))
    return img.astype(np.uint16)


def compute_blank_average(pos: str) -> np.ndarray:
    """Load blank pair for a position and return their average as uint16."""
    blank1 = load_image(REGISTERED_DIR / f"Blank1{pos}.tif")
    blank13 = load_image(REGISTERED_DIR / f"Blank13{pos}.tif")
    avg = ((blank1.astype(np.float32) + blank13.astype(np.float32)) / 2.0)
    return avg.astype(np.uint16)


def parse_param_file(marker: str) -> dict | None:
    """Parse an existing parameter file and return subtraction params.

    Returns dict with blank_clip_factor and blank_scale_factor, or None
    if no param file exists.
    """
    param_path = PARAMS_DIR / f"{marker}_param.txt"
    if not param_path.exists():
        return None

    text = param_path.read_text()
    params = {}

    for line in text.splitlines():
        if line.startswith("blank_clip_factor:"):
            try:
                params["blank_clip_factor"] = int(line.split(":", 1)[1].strip())
            except ValueError:
                pass
        elif line.startswith("background_scale_factor:"):
            try:
                params["blank_scale_factor"] = float(line.split(":", 1)[1].strip())
            except ValueError:
                pass

    if "blank_clip_factor" in params and "blank_scale_factor" in params:
        return params
    return None


def process_channel(
    marker: str,
    blank_avg: np.ndarray,
    use_auto: bool = False,
) -> dict:
    """Process a single channel and return result info."""
    signal_path = REGISTERED_DIR / f"{marker}.tif"
    if not signal_path.exists():
        # Try HLA-DR variant naming
        if marker == "HLADR":
            signal_path = REGISTERED_DIR / "HLADR.tif"
            if not signal_path.exists():
                signal_path = REGISTERED_DIR / "HLA-DR.tif"
        if not signal_path.exists():
            return {"marker": marker, "status": "MISSING", "error": f"No file found"}

    signal = load_image(signal_path)
    result_info = {
        "marker": marker,
        "shape": signal.shape,
        "input_p1": int(np.percentile(signal, 1)),
        "input_p99": int(np.percentile(signal, 99)),
    }

    # Get parameters
    params = None if use_auto else parse_param_file(marker)

    if params is not None:
        result_info["param_source"] = "file"
        result_info["blank_clip_factor"] = params["blank_clip_factor"]
        result_info["blank_scale_factor"] = params["blank_scale_factor"]
    else:
        # Auto-analyze
        suggested = analyze_for_subtraction(
            signal, blank_avg,
            tissue_type="spleen",
            marker_name=marker,
        )
        params = {
            "blank_clip_factor": suggested["blank_clip_factor"],
            "blank_scale_factor": suggested["blank_scale_factor"],
        }
        result_info["param_source"] = "auto"
        result_info["confidence"] = suggested.get("confidence", None)
        result_info["blank_clip_factor"] = params["blank_clip_factor"]
        result_info["blank_scale_factor"] = params["blank_scale_factor"]

    # Run subtraction
    subtracted = subtract_autofluorescence(
        signal, blank_avg,
        blank_clip_factor=params["blank_clip_factor"],
        blank_scale_factor=params["blank_scale_factor"],
    )

    # Quality metrics
    quality = compute_subtraction_quality(signal, subtracted, blank_avg)
    result_info["quality_score"] = quality.get("quality_score", None)
    result_info["signal_preservation"] = quality.get("signal_preservation", None)
    result_info["af_removal"] = quality.get("af_removal", None)
    result_info["snr_improvement"] = quality.get("snr_improvement", None)

    result_info["output_p1"] = int(np.percentile(subtracted, 1))
    result_info["output_p99"] = int(np.percentile(subtracted, 99))

    # Save
    output_path = OUTPUT_DIR / f"{marker}.tif"
    tifffile.imwrite(str(output_path), subtracted)
    result_info["status"] = "OK"

    return result_info


def main():
    parser = argparse.ArgumentParser(description="HDL73 signal isolation")
    parser.add_argument("--force", action="store_true",
                        help="Re-process already existing channels")
    args = parser.parse_args()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Determine which channels to skip
    existing = {p.stem for p in OUTPUT_DIR.glob("*.tif")}
    if not args.force:
        to_process = [ch for ch in ALL_SIGNAL_CHANNELS if ch not in existing]
        skipped = [ch for ch in ALL_SIGNAL_CHANNELS if ch in existing]
    else:
        to_process = list(ALL_SIGNAL_CHANNELS)
        skipped = []

    if skipped:
        print(f"Skipping {len(skipped)} already-processed channels: {', '.join(skipped)}")
    print(f"Processing {len(to_process)} channels")
    print()

    # Pre-compute blank averages (one per position)
    print("Computing blank averages...")
    blank_avgs = {}
    for pos in ("a", "b", "c"):
        blank_avgs[pos] = compute_blank_average(pos)
        print(f"  Position {pos}: shape={blank_avgs[pos].shape}, "
              f"median={int(np.median(blank_avgs[pos]))}")
    print()

    # Process each channel
    results = []
    for i, marker in enumerate(to_process, 1):
        pos = CHANNEL_TO_POSITION[marker]
        print(f"[{i}/{len(to_process)}] {marker} (blank position {pos})...", end=" ", flush=True)
        try:
            info = process_channel(marker, blank_avgs[pos])
            results.append(info)
            if info["status"] == "OK":
                src = info["param_source"]
                qs = info.get("quality_score")
                qs_str = f"{qs:.3f}" if qs is not None else "N/A"
                print(f"OK  [src={src}, clip={info['blank_clip_factor']}, "
                      f"scale={info['blank_scale_factor']:.1f}, quality={qs_str}]")
            else:
                print(f"SKIPPED: {info.get('error', 'unknown')}")
        except Exception as e:
            print(f"ERROR: {e}")
            results.append({"marker": marker, "status": "ERROR", "error": str(e)})

    # Copy DAPI
    dapi_src = REGISTERED_DIR / "DAPI.tif"
    dapi_dst = OUTPUT_DIR / "DAPI.tif"
    if dapi_dst.exists() and not args.force:
        print("\nDAPI: already exists, skipping")
    elif dapi_src.exists():
        shutil.copy2(str(dapi_src), str(dapi_dst))
        print("\nDAPI: copied (no subtraction needed)")
    else:
        print("\nDAPI: WARNING - source file not found!")

    # Summary table
    print("\n" + "=" * 100)
    print(f"{'Marker':<14} {'Status':<8} {'Source':<6} {'Clip':>6} {'Scale':>6} "
          f"{'Quality':>8} {'SigPres':>8} {'AF_Rem':>8} {'SNR_Imp':>8} "
          f"{'In_p99':>8} {'Out_p99':>8}")
    print("-" * 100)
    for r in sorted(results, key=lambda x: x["marker"]):
        if r["status"] == "OK":
            qs = r.get("quality_score")
            sp = r.get("signal_preservation")
            af = r.get("af_removal")
            snr = r.get("snr_improvement")
            print(f"{r['marker']:<14} {r['status']:<8} {r.get('param_source',''):<6} "
                  f"{r.get('blank_clip_factor',''):>6} {r.get('blank_scale_factor',''):>6.1f} "
                  f"{qs:>8.3f} {sp:>8.3f} {af:>8.3f} {snr:>+8.2f} "
                  f"{r.get('input_p99',''):>8} {r.get('output_p99',''):>8}")
        else:
            print(f"{r['marker']:<14} {r['status']:<8} {r.get('error','')}")
    print("=" * 100)

    # Final file count
    final_count = len(list(OUTPUT_DIR.glob("*.tif")))
    print(f"\nTotal files in ImageJ/: {final_count}")

    # Check for any concerning quality scores
    low_quality = [r for r in results
                   if r["status"] == "OK"
                   and r.get("quality_score") is not None
                   and r["quality_score"] < 0.5]
    if low_quality:
        print("\nWARNING - Low quality scores (<0.5):")
        for r in low_quality:
            print(f"  {r['marker']}: quality={r['quality_score']:.3f}, "
                  f"signal_preservation={r.get('signal_preservation', 'N/A')}")

    low_signal = [r for r in results
                  if r["status"] == "OK"
                  and r.get("signal_preservation") is not None
                  and r["signal_preservation"] < 0.3]
    if low_signal:
        print("\nWARNING - Low signal preservation (<0.3):")
        for r in low_signal:
            print(f"  {r['marker']}: signal_preservation={r['signal_preservation']:.3f}")


if __name__ == "__main__":
    main()
