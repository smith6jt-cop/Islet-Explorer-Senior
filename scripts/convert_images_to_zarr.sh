#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Batch-convert Islet Explorer OME-TIFF images to OME-Zarr (Phase 2.1).
#
# The app currently serves pyramidal OME-TIFFs via nginx with HTTP Range
# headers. For the scaling work we want the image tiles to live in
# S3-compatible object storage (MinIO on-prem, or cloud later), which means
# switching to OME-Zarr so Viv / Avivator can stream chunks directly.
#
# This script:
#   1. Enumerates every .ome.tiff / .ome.tif under a source directory.
#   2. Skips images whose .zarr store already exists in the destination.
#   3. Delegates each conversion to the single-image wrapper
#      `scripts/convert_to_omezarr.sh`, which prefers Dockerised
#      bioformats2raw and falls back to the bundled JAR.
#
# Usage:
#   scripts/convert_images_to_zarr.sh                         # defaults
#   scripts/convert_images_to_zarr.sh app/shiny_app/www/local_images data/zarr_images
#
# Environment:
#   B2R_TILE_WIDTH   (default 1024) -- tile width passed to bioformats2raw
#   B2R_TILE_HEIGHT  (default 1024)
#   B2R_MAX_WORKERS  (default 4)
#   B2R_COMPRESSION  (default blosc)
#
# Notes on CORS: after upload to MinIO, the bucket must advertise CORS so
# Avivator/Viv can fetch chunks from the browser. See docs/scaling.md (if
# present) or the Phase 2.2 section of the scaling brief.
# -----------------------------------------------------------------------------
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

INPUT_DIR="${1:-$REPO_ROOT/app/shiny_app/www/local_images}"
OUTPUT_DIR="${2:-$REPO_ROOT/data/zarr_images}"

TILE_WIDTH="${B2R_TILE_WIDTH:-1024}"
TILE_HEIGHT="${B2R_TILE_HEIGHT:-1024}"
MAX_WORKERS="${B2R_MAX_WORKERS:-4}"
COMPRESSION="${B2R_COMPRESSION:-blosc}"

CONVERTER="$SCRIPT_DIR/convert_to_omezarr.sh"

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "[ERROR] Input directory does not exist: $INPUT_DIR" >&2
  exit 2
fi

if [[ ! -x "$CONVERTER" ]]; then
  echo "[ERROR] Single-image converter not found or not executable: $CONVERTER" >&2
  echo "        Expected scripts/convert_to_omezarr.sh to be present." >&2
  exit 3
fi

mkdir -p "$OUTPUT_DIR"

shopt -s nullglob nocaseglob
tiffs=("$INPUT_DIR"/*.ome.tif "$INPUT_DIR"/*.ome.tiff)
shopt -u nullglob nocaseglob

if (( ${#tiffs[@]} == 0 )); then
  echo "[INFO] No .ome.tif/.ome.tiff files found under $INPUT_DIR"
  exit 0
fi

echo "[INFO] Converting ${#tiffs[@]} image(s)"
echo "       source : $INPUT_DIR"
echo "       target : $OUTPUT_DIR"
echo "       tile   : ${TILE_WIDTH}x${TILE_HEIGHT}"
echo "       workers: $MAX_WORKERS  compression: $COMPRESSION"
echo

converted=0
skipped=0
failed=0

for tiff in "${tiffs[@]}"; do
  base="$(basename "$tiff")"
  # Strip a single trailing .ome.tif[f] regardless of case.
  stem="$(echo "$base" | sed -E 's/\.ome\.tiff?$//i')"
  out_zarr="$OUTPUT_DIR/${stem}.zarr"

  if [[ -d "$out_zarr" && -n "$(ls -A "$out_zarr" 2>/dev/null)" ]]; then
    echo "[SKIP] $base -> $out_zarr already populated"
    skipped=$((skipped + 1))
    continue
  fi

  echo "[CONVERT] $base -> $out_zarr"
  if "$CONVERTER" "$tiff" "$out_zarr" \
        --tile_width "$TILE_WIDTH" \
        --tile_height "$TILE_HEIGHT" \
        --max_workers "$MAX_WORKERS" \
        --compression "$COMPRESSION"; then
    converted=$((converted + 1))
  else
    echo "[FAIL] $base" >&2
    failed=$((failed + 1))
  fi
done

echo
echo "[DONE] converted=$converted skipped=$skipped failed=$failed"

if (( failed > 0 )); then
  exit 1
fi
