#!/usr/bin/env bash
# Archive the 15-donor cohort before the 22-donor rebuild, and reclaim disk.
#
#   - MOVE  load-bearing / derived 15-donor outputs -> data/archive_15donor_<date>/
#           (move = free/instant on the same filesystem; data/ is gitignored)
#   - DELETE pure regeneratable scratch (redsea_scratch, ~44 GB) and the deprecated
#           pre-REDSEA restore family + stray leftovers
#   - cells/ and cells_redsea/ are LEFT IN PLACE by default (huge + reproducible; the
#           rerun regenerates them). Pass --archive-cells to snapshot them too.
#
# SAFETY: dry-run by default (prints what it would do). Pass --apply to execute.
# Run this only when you are ready to rebuild — moving the cohort files makes the app's
# app_data/ symlinks dangle until the 22-donor files are rebuilt.
#
#   bash scripts/senior/archive_15donor.sh            # dry run (default)
#   bash scripts/senior/archive_15donor.sh --apply    # execute
#   bash scripts/senior/archive_15donor.sh --apply --archive-cells
set -euo pipefail

DATA="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../data" && pwd)"
STAMP="2026-07-02"
ARCH="$DATA/archive_15donor_$STAMP"
APPLY=0
ARCHIVE_CELLS=0
for arg in "$@"; do
  case "$arg" in
    --apply) APPLY=1 ;;
    --archive-cells) ARCHIVE_CELLS=1 ;;
    *) echo "unknown arg: $arg"; exit 2 ;;
  esac
done

LOG="$ARCH/archive_manifest.txt"
say() { echo "$@"; [ "$APPLY" = 1 ] && echo "$@" >> "$LOG" || true; }
run() { if [ "$APPLY" = 1 ]; then "$@"; else echo "DRY: $*"; fi; }

echo "== archive_15donor ($([ "$APPLY" = 1 ] && echo APPLY || echo DRY-RUN)) =="
echo "data dir: $DATA"
df -h "$DATA" | sed 's/^/  /'
[ "$APPLY" = 1 ] && { mkdir -p "$ARCH"; : > "$LOG"; }

# ---- 1. MOVE: derived / cohort 15-donor outputs (provenance snapshot) ----
MOVE=(
  singlecell_protein.h5ad
  islets_core_senior.h5ad islets_peri_senior.h5ad islet_explorer_senior.h5ad
  restore_gated_redsea restore_redsea restore_thresholds_redsea.csv
  phenotype/broad phenotype/broad_lineage_composition.png
  phenotype/celltype_marker_dotplot.png phenotype/celltype_marker_heatmap.png
  phenotype/qupath_class phenotype/qupath_import
  neighborhood_metrics.csv islet_targets.parquet islet_composition.parquet parquet
  app_data/cells app_data/json app_data/islet_spatial_lookup.csv
  redsea_full_validation.png redsea_cohort_validation.png redsea_approx_validation.png
  redsea_residual_investigation.png redsea_doublet_investigation.png
  restore_reference_fix_validation.png immune_false_positive_investigation.png
  dpt_gap_diagnosis.png marker_audit_distributions.pdf marker_audit_cluster_means.csv
  marker_audit_marker_audit.csv cluster_labels_curated.csv cluster_labels_curated_with_rationale.csv
)
[ "$ARCHIVE_CELLS" = 1 ] && MOVE+=(cells cells_redsea)

say "-- MOVE -> $ARCH --"
for item in "${MOVE[@]}"; do
  src="$DATA/$item"
  if [ -e "$src" ]; then
    dst="$ARCH/$item"; run mkdir -p "$(dirname "$dst")"; run mv "$src" "$dst"; say "moved: $item"
  else
    say "skip (absent): $item"
  fi
done

# ---- 2. DELETE: pure scratch + deprecated pre-REDSEA family + strays ----
DEL=(
  redsea_scratch
  restore restore_gated restore_thresholds.csv restore_redsea_thresholds.csv
  restore_raw6539_thresholds.csv
)
say "-- DELETE (regeneratable scratch / deprecated) --"
for item in "${DEL[@]}"; do
  tgt="$DATA/$item"
  if [ -e "$tgt" ]; then run rm -rf "$tgt"; say "deleted: $item"; else say "skip (absent): $item"; fi
done
# stray QuPath "Root object" export leftovers in the drill-down dir (if cells not moved)
if [ "$ARCHIVE_CELLS" = 0 ]; then
  while IFS= read -r f; do run rm -f "$f"; say "deleted stray: ${f#$DATA/}"; done \
    < <(find "$DATA/app_data/cells" -name '*Root object*'.csv 2>/dev/null || true)
fi

echo "== done ($([ "$APPLY" = 1 ] && echo applied || echo dry-run)) =="
df -h "$DATA" | sed 's/^/  /'
[ "$APPLY" = 1 ] && echo "manifest: $LOG" || echo "(re-run with --apply to execute)"
