#!/usr/bin/env bash
#
# convert_viewer_ometiff.sh — batch-convert Akoya PhenoCycler .qptiff slides to
# pyramidal, tiled OME-TIFF for the Islet Explorer "Viewer" tab (Avivator/Viv).
#
# Pipeline (no Docker; local zip-app launchers on Java 21):
#   bioformats2raw  -s 0 --use-existing-resolutions  <qptiff> <zarr>   (Baseline series only)
#   raw2ometiff     --compression LZW               <zarr>  <out.ome.tiff>
#
# Output: $OUT_DIR/<case_id>.ome.tiff  (LZW = lossless; matches source MFI).
# The <case_id> basename is what the Viewer dropdown lists and what
# find_image_for_case() substring-matches for the drill-down "jump to Viewer" path.
#
# Design notes:
#   * STORAGE — the intermediate Zarr is thousands of small chunk files (random I/O);
#     it MUST live on SSD/NVMe. On a spinning disk it saturates the drive and stalls
#     the whole box (dirty-page backlog → every `sync` hangs; this happened on the
#     18 TB /data HDD, 2026-07-02). So TMPBASE defaults to NVMe and the script
#     HARD-ERRORS if TMPBASE is on a rotational disk. The final OME-TIFF is one big
#     sequential write; OUT_DIR defaults to NVMe too and only WARNS if it's an HDD.
#     A flock guarantees ONE conversion at a time — concurrent writers on one spindle
#     are what compounded the stall.
#   * Idempotent: a donor whose final .ome.tiff already exists is skipped, so the
#     batch is resumable after an interrupt. Partial output is written to a
#     .partial path and atomically mv'd on success.
#   * --use-existing-resolutions reuses the qptiff's own 5-level pyramid (fast,
#     exact). If a slide's reader can't expose sub-resolutions, it retries with
#     bioformats2raw regenerating the pyramid.
#
# Usage:
#   scripts/senior/convert_viewer_ometiff.sh 6539              # one donor
#   scripts/senior/convert_viewer_ometiff.sh 6539 6374 115     # several
#   scripts/senior/convert_viewer_ometiff.sh all15             # all 15 app donors
#
# Env overrides: OUT_DIR, TMPBASE (both default to NVMe; TMPBASE MUST be SSD/NVMe),
#                IMAGES_DIR, TILE, COMPRESSION, B2R_MAX_WORKERS, MIN_FREE_GB, LOCKFILE.
set -uo pipefail

# ---- config ---------------------------------------------------------------
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
IMAGES_DIR="${IMAGES_DIR:-/home/smith6jt/IO60panc2nd/Images}"
# Default OUT_DIR/TMPBASE to NVMe — NEVER the /data HDD (see "STORAGE" note above). Override per run.
OUT_DIR="${OUT_DIR:-$HOME/islet_ome_tiff}"       # final OME-TIFFs (persistent; put on NVMe)
TMPBASE="${TMPBASE:-/tmp/b2r_scratch}"           # intermediate Zarr — MUST be SSD/NVMe (/tmp is NVMe here)
TILE="${TILE:-1024}"
COMPRESSION="${COMPRESSION:-LZW}"
MIN_FREE_GB="${MIN_FREE_GB:-300}"                # abort a slide if OUT_DIR free < this
# Python with generate-tiff-offsets (Viv IFD index -> <case_id>.offsets.json). Optional.
OFFSETS_PYTHON="${OFFSETS_PYTHON:-/home/smith6jt/miniconda3/envs/xenium_analysis/bin/python}"

B2R="$(ls -d "$REPO"/scripts/tools/bioformats2raw-*/bin/bioformats2raw 2>/dev/null | sort -V | tail -1)"
R2O="$(ls -d "$REPO"/scripts/tools/raw2ometiff-*/bin/raw2ometiff 2>/dev/null | sort -V | tail -1)"
B2R_MAX_WORKERS="${B2R_MAX_WORKERS:-}"

# ---- case_id -> qptiff filename ----
# Batch 1 (15): shared 59-ch panel (SST@35, CD38@16). Batch 2 (7): different panel
# (SST@1; 6450/6505/6534/6623 are 59-ch with b-Catenin1@16, GCG@46; 6523/6566/6591 are 58-ch,
# GCG@45, no CD38/b-Catenin1). Per-image channel names still come from each OME-TIFF's OME-XML.
declare -A MAP=(
  [115]="HDL115pancLN_Scan1.er.qptiff"
  [6374]="6374_Scan1.er.qptiff"
  [6380]="6380_Scan1.er.qptiff"
  [6436]="6436_Scan1.er.qptiff"
  [6442]="6442_Scan1.er.qptiff"
  [6457]="6457_Scan1.er.qptiff"
  [6461]="6461_Scan1.er.qptiff"
  [6476]="6476_Scan2.er.qptiff"
  [6482]="6482_Scan1.er.qptiff"
  [6493]="6493_Scan1.er.qptiff"
  [6521]="6521_Scan1.er.qptiff"
  [6539]="6539_Scan1.er.qptiff"
  [6551]="6551_Scan1.er.qptiff"
  [6579]="6579_Scan2.er.qptiff"
  [6593]="6593_Scan1.er.qptiff"
  [6450]="6450panc_Scan1.er.qptiff"
  [6505]="6505panc_Scan1.er.qptiff"
  [6523]="6523panc_Scan1.er.qptiff"
  [6534]="6534pacn_Scan1.er.qptiff"
  [6566]="6566panc_Scan1.er.qptiff"
  [6591]="6591panc_Scan1.er.qptiff"
  [6623]="6623panc_Scan1.er.qptiff"
)
ALL15=(115 6374 6380 6436 6442 6457 6461 6476 6482 6493 6521 6539 6551 6579 6593)
ALL22=("${ALL15[@]}" 6450 6505 6523 6534 6566 6591 6623)

usage() { echo "Usage: $0 <case_id...> | all15 | all22" >&2; exit 1; }
[[ $# -ge 1 ]] || usage
case "$1" in
  all22) TARGETS=("${ALL22[@]}") ;;
  all15) TARGETS=("${ALL15[@]}") ;;
  *)     TARGETS=("$@") ;;
esac

[[ -x "$B2R" ]] || { echo "ERROR: bioformats2raw launcher not found under scripts/tools/ — run scripts/get_bioformats2raw.sh or unzip the release." >&2; exit 1; }
[[ -x "$R2O" ]] || { echo "ERROR: raw2ometiff launcher not found under scripts/tools/." >&2; exit 1; }
mkdir -p "$OUT_DIR" "$TMPBASE"

# ---- storage safety: keep heavy random I/O off spinning disks --------------
is_rotational() {  # true if the path's backing block device is an HDD (ROTA=1)
  local src; src=$(df --output=source "$1" 2>/dev/null | tail -1)
  [[ -n "$src" ]] && [[ "$(lsblk -no ROTA "$src" 2>/dev/null | head -1 | tr -d ' ')" == "1" ]]
}
if is_rotational "$TMPBASE"; then
  echo "ERROR: TMPBASE ($TMPBASE) is on a ROTATIONAL disk (HDD). bioformats2raw writes thousands of" >&2
  echo "       small Zarr chunk files (random I/O) — on an HDD that saturates the disk and stalls the" >&2
  echo "       whole machine (dirty-page backlog → every sync hangs). Point TMPBASE at SSD/NVMe:" >&2
  echo "         TMPBASE=/tmp/b2r_scratch $0 $*" >&2
  exit 1
fi
if is_rotational "$OUT_DIR"; then
  echo "WARNING: OUT_DIR ($OUT_DIR) is on a ROTATIONAL disk (HDD). The final OME-TIFF is a large" >&2
  echo "         sequential write (a CMR HDD handles it) and serving is read-only, but you asked to" >&2
  echo "         keep output off the HDD. Set OUT_DIR to NVMe to avoid it entirely." >&2
fi

# ---- single-instance lock: never run two conversions on one spindle --------
LOCKFILE="${LOCKFILE:-/tmp/convert_viewer_ometiff.lock}"
exec 9>"$LOCKFILE"
if ! flock -n 9; then
  echo "ERROR: another conversion holds the lock ($LOCKFILE). Refusing to run concurrently —" >&2
  echo "       parallel writers on one disk are what compounded the 2026-07-02 stall." >&2
  exit 1
fi

log() { echo "[$(date '+%F %T')] $*"; }
free_gb() { df -PBG "$1" | awk 'NR==2{gsub("G","",$4); print $4}'; }

WORKERS_ARG=()
[[ -n "$B2R_MAX_WORKERS" ]] && WORKERS_ARG=(--max-workers "$B2R_MAX_WORKERS")

log "converter: $B2R"
log "converter: $R2O"
log "output:    $OUT_DIR   (free $(free_gb "$OUT_DIR")GB)"
log "targets:   ${TARGETS[*]}"

rc_all=0
for cid in "${TARGETS[@]}"; do
  fn="${MAP[$cid]:-}"
  if [[ -z "$fn" ]]; then log "SKIP $cid: not in app-donor map"; rc_all=1; continue; fi
  in="$IMAGES_DIR/$fn"
  out="$OUT_DIR/$cid.ome.tiff"
  partial="$OUT_DIR/.$cid.ome.tiff.partial"

  if [[ ! -f "$in" ]]; then log "SKIP $cid: input missing ($in)"; rc_all=1; continue; fi
  if [[ -f "$out" ]]; then log "SKIP $cid: already converted ($(du -h "$out" | cut -f1))"; continue; fi

  free=$(free_gb "$OUT_DIR")
  if [[ "${free:-0}" -lt "$MIN_FREE_GB" ]]; then
    log "ABORT $cid: only ${free}GB free on $OUT_DIR (< ${MIN_FREE_GB}GB)"; rc_all=1; break
  fi

  insz=$(du -h "$in" | cut -f1)
  log "=== $cid: $fn ($insz) ==="
  zarr="$(mktemp -d "$TMPBASE/${cid}.XXXXXX")/raw"
  rm -f "$partial"
  t0=$SECONDS

  ok=1
  log "$cid: bioformats2raw (series 0, use-existing-resolutions) ..."
  if ! "$B2R" -s 0 --use-existing-resolutions -w "$TILE" -h "$TILE" "${WORKERS_ARG[@]}" "$in" "$zarr"; then
    log "$cid: use-existing-resolutions failed; retrying with pyramid regeneration ..."
    rm -rf "$zarr"
    if ! "$B2R" -s 0 -w "$TILE" -h "$TILE" "${WORKERS_ARG[@]}" "$in" "$zarr"; then
      log "$cid: bioformats2raw FAILED"; ok=0
    fi
  fi

  if [[ $ok -eq 1 ]]; then
    log "$cid: raw2ometiff (--compression $COMPRESSION) ..."
    if ! "$R2O" --compression "$COMPRESSION" "$zarr" "$partial"; then
      log "$cid: raw2ometiff FAILED"; ok=0
    fi
  fi

  rm -rf "$(dirname "$zarr")"

  if [[ $ok -eq 1 && -f "$partial" ]]; then
    mv -f "$partial" "$out"
    log "$cid: DONE -> $out ($(du -h "$out" | cut -f1), $((SECONDS-t0))s)"
    # Viv IFD-index sidecar -> <case_id>.offsets.json (removes Avivator's "no IFD index" banner)
    if [[ -x "$OFFSETS_PYTHON" ]] && "$OFFSETS_PYTHON" -c "import generate_tiff_offsets" 2>/dev/null; then
      if "$OFFSETS_PYTHON" -c "from generate_tiff_offsets import generate_tiff_offsets as g; g('$out')" 2>/dev/null; then
        log "$cid: offsets.json written"
      else
        log "$cid: offsets generation failed (non-fatal — Avivator falls back to IFD scan)"
      fi
    fi
  else
    rm -f "$partial"
    log "$cid: FAILED after $((SECONDS-t0))s"; rc_all=1
  fi
done

log "batch complete (exit $rc_all)"
exit $rc_all
