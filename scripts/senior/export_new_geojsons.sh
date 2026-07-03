#!/usr/bin/env bash
# Export per-cell GeoJSON (REDSEA input) for any donor image that doesn't have one yet.
#
# Each image is exported headlessly via QuPath's `script` runner + the project's
# export_cells_geojson.groovy -> data/redsea_scratch/geojson/cells__<tag>.geojson,
# where <tag> is the QuPath image name with [^A-Za-z0-9._-] -> '_' (exactly what
# redsea_full.py:resolve_paths reconstructs from the cells parquet's `image` column,
# e.g. "6450panc_Scan1.er.qptiff - resolution #1").
#
# By default it exports ONLY images missing a GeoJSON — i.e. the 7 new batch-2 donors
# (the 15 batch-1 GeoJSONs already exist and are reused). Pass --all to re-export every
# image (force).
#
# ⚠️ CLOSE THE QuPath GUI FIRST. The GUI holds the JVM heap and can lock project.qpproj;
#    each export needs a fresh headless JVM. This script refuses to run if a GUI QuPath
#    process is detected (override with --ignore-gui if you know it's safe).
#
#   bash scripts/senior/export_new_geojsons.sh              # export the missing (7 new)
#   bash scripts/senior/export_new_geojsons.sh --all        # re-export all 22
#   bash scripts/senior/export_new_geojsons.sh --ignore-gui # skip the GUI-running guard
set -euo pipefail

QUPATH="$HOME/QuPath/bin/QuPath"
PROJECT="$HOME/IO60panc2nd/project.qpproj"
GROOVY="$HOME/IO60panc2nd/scripts/export_helpers/export_cells_geojson.groovy"
IMAGES="$HOME/IO60panc2nd/Images"
GEOJSON="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../data" && pwd)/redsea_scratch/geojson"

FORCE=0; IGNORE_GUI=0
for arg in "$@"; do
  case "$arg" in
    --all) FORCE=1 ;;
    --ignore-gui) IGNORE_GUI=1 ;;
    *) echo "unknown arg: $arg"; exit 2 ;;
  esac
done

for f in "$QUPATH" "$PROJECT" "$GROOVY"; do
  [ -e "$f" ] || { echo "ERROR: missing $f"; exit 1; }
done

# --- GUI guard: any QuPath java process whose cmdline lacks ' script ' is the GUI.
# (The headless runner is exec'd as java with 'script' in argv, so match on the jar/app
# name case-insensitively and drop the ' script ' lines + this script's own line.) ---
if [ "$IGNORE_GUI" = 0 ]; then
  gui=$(pgrep -af -i qupath 2>/dev/null | grep -iv -e ' script ' -e 'export_new_geojsons' -e 'pgrep' || true)
  if [ -n "$gui" ]; then
    echo "ERROR: a QuPath GUI appears to be running:"; echo "$gui" | sed 's/^/  /'
    echo "Close it first (it holds the heap and can lock the project), then re-run."
    echo "Override with --ignore-gui if you are sure it is safe."
    exit 1
  fi
fi

mkdir -p "$GEOJSON"
sanitize() { printf '%s' "$1" | sed 's/[^A-Za-z0-9._-]/_/g'; }

echo "== export_new_geojsons ($([ "$FORCE" = 1 ] && echo 'ALL/force' || echo 'missing only')) =="
echo "  images : $IMAGES"
echo "  out    : $GEOJSON"
n=0; skipped=0; failed=0
for q in "$IMAGES"/*.qptiff; do
  base=$(basename "$q")                 # e.g. 6450panc_Scan1.er.qptiff
  img="$base - resolution #1"           # QuPath project entry name (matches the cells parquet)
  out="$GEOJSON/cells__$(sanitize "$img").geojson"
  if [ -e "$out" ] && [ "$FORCE" = 0 ]; then
    echo "skip (exists): $base"; skipped=$((skipped+1)); continue
  fi
  n=$((n+1))
  echo "[export $n] $img"
  if "$QUPATH" script --project "$PROJECT" --image "$img" "$GROOVY"; then
    ls -lh "$out" 2>/dev/null | awk '{print "         ->", $5, $9}'
  else
    echo "         FAILED: $img"; failed=$((failed+1))
  fi
done
echo "== done: exported=$n skipped=$skipped failed=$failed =="
echo "GeoJSONs in $GEOJSON (total: $(ls -1 "$GEOJSON"/cells__*.geojson 2>/dev/null | wc -l))"
[ "$failed" = 0 ] || { echo "some exports failed — re-run for the failures"; exit 1; }
echo "Next: python scripts/senior/redsea_full.py --all   (recomputes cells_redsea for all 22)"
