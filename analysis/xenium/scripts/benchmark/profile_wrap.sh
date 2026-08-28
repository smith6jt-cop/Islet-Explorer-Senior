#!/usr/bin/env bash
# =============================================================================
# profile_wrap.sh — uniform wall-time + peak-RSS profiler for ONE benchmark
# method run, across Python and R. Writes a metrics JSON the evaluator reads.
#
# The CALLER is responsible for the environment (conda activate + any
# LD_PRELOAD); this wrapper only measures. Each method should be run as its own
# process so peak RSS is clean.
#
# Usage:
#   profile_wrap.sh <sample> <variant> <method> -- <command ...>
# Optional env:
#   BENCH_TIMEOUT=21600   # seconds; on exceed -> status="capped" (no silent cap)
#
# Output: data/processed/benchmark/<sample>/results/<method>_<variant>_metrics.json
# =============================================================================
set -uo pipefail

ROOT="${XENIUM_ANALYSIS_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
sample="${1:?sample}"; variant="${2:?variant}"; method="${3:?method}"; shift 3
[ "${1:-}" = "--" ] && shift

RESDIR="${ROOT}/data/processed/benchmark/${sample}/results"
mkdir -p "${RESDIR}"
TIMEFILE="$(mktemp)"
METRICS="${RESDIR}/${method}_${variant}_metrics.json"

# Optional hard timeout for the scalability sweep.
RUNNER=(/usr/bin/time -v -o "${TIMEFILE}")
if [[ -n "${BENCH_TIMEOUT:-}" ]]; then
  RUNNER+=(timeout --signal=TERM "${BENCH_TIMEOUT}")
fi

echo "[profile] ${method} ${sample}/${variant} :: $*"
"${RUNNER[@]}" "$@"
EXIT=$?

# ---- parse /usr/bin/time -v ----
RSS_KB="$(awk -F': ' '/Maximum resident set size/{print $2}' "${TIMEFILE}" 2>/dev/null)"
WALL_RAW="$(awk -F'): ' '/Elapsed \(wall clock\)/{print $2}' "${TIMEFILE}" 2>/dev/null)"
WALL_S="$(awk -v t="${WALL_RAW:-}" 'BEGIN{n=split(t,a,":");
  if(n==3) printf "%.2f", a[1]*3600+a[2]*60+a[3];
  else if(n==2) printf "%.2f", a[1]*60+a[2];
  else if(n==1 && t!="") printf "%.2f", a[1];
  else printf "null"}')"
RSS_MB="$(awk -v k="${RSS_KB:-}" 'BEGIN{if(k=="") printf "null"; else printf "%.1f", k/1024}')"

# ---- status ----
if [[ "${EXIT}" -eq 124 || "${EXIT}" -eq 137 ]]; then
  STATUS="capped"            # timeout (TERM=124 / KILL=137)
elif [[ "${EXIT}" -eq 0 ]]; then
  STATUS="ok"
else
  STATUS="failed"
fi

# ---- n_domains (domain methods) from labels; misty -> null ----
LABELS="${RESDIR}/${method}_${variant}_labels.csv"
if [[ -f "${LABELS}" ]]; then
  NDOM="$(tail -n +2 "${LABELS}" | cut -d, -f2 | sort -u | wc -l)"
else
  NDOM="null"
fi

# ---- n_cells from cropinfo (best-effort) ----
CROP="${ROOT}/data/processed/benchmark/${sample}/${sample}_bench_${variant}_cropinfo.json"
NCELLS="$(grep -o '"n_cells": *[0-9]*' "${CROP}" 2>/dev/null | head -1 | grep -o '[0-9]*' || true)"
[[ -z "${NCELLS}" ]] && NCELLS="null"

cat > "${METRICS}" <<EOF
{
  "sample": "${sample}",
  "variant": "${variant}",
  "method": "${method}",
  "n_cells": ${NCELLS},
  "wall_s": ${WALL_S},
  "peak_rss_mb": ${RSS_MB},
  "n_domains": ${NDOM},
  "exit_code": ${EXIT},
  "status": "${STATUS}"
}
EOF
rm -f "${TIMEFILE}"
echo "[profile] ${method} ${sample}/${variant} -> status=${STATUS} wall=${WALL_S}s rss=${RSS_MB}MB ndom=${NDOM}"
exit "${EXIT}"
