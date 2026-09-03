#!/bin/bash
# The measurement the production sizing depends on, taken on the hardware that will do the
# work, before the hardware is committed to it.
#
# Everything downstream is sized from one number: how KODAMA's cost at full donor scale
# splits between the Monte-Carlo loop (which shards, and which cannot use more than about
# sqrt(m)/2 cores in one process) and the fixed kNN + dissimilarity (which does not shard,
# and does use its cores). The workstation ladder measured that split only up to 50,000
# cells, and an attempt to measure it here at n=20,000 was contended by other work and is
# an upper bound at best -- so it is measured properly, once, on an idle node.
#
# Three M points at fixed n and fixed cores give the line: T(M) = fixed + M * per_iteration.
# The smallest cohort donor (6450, 255,149 cells) and the largest (6476, 1,109,281) bracket
# the range; the array runs the ladder on whichever donors are named.
#
# LANDMARKS IS MEASURED, NOT ASSUMED. The workstation ladder on the 5,101-gene matrix showed
# cost rising steeply with landmarks (n=5,000: 500/M5 1583 s, 1000/M10 5490 s, 2000/M20
# 17246 s), but an attempt to repeat that on 50 PCs came out the other way -- n=20,000, M=5,
# 8 cores: landmarks 500 787.8 s vs landmarks 2000 451.9 s -- on a contended workstation
# where the two runs did not have the machine to themselves. The direction is genuinely
# unknown on the production input, and landmarks is the one knob whose scientific return is
# provably capped (KODAMA derives `neighbors` as min(landmarks, .., 500), so past 500 it buys
# a better per-iteration partition but NOT a wider neighbour graph). So it is the knob to
# measure, and the knob to trade against if anything has to give -- never the cell count,
# and never M.
#
#   sbatch --array=0-11 --dependency=afterok:$PREP scripts/xenium_celltype/hpg/slurm/05_calibrate.sh
#SBATCH --job-name=xct_cal
#SBATCH --account=maigan
#SBATCH --qos=maigan
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64gb
#SBATCH --time=24:00:00
#SBATCH --exclusive
#SBATCH --array=0-11
#SBATCH --output=logs/xct_cal_%A_%a.out
#SBATCH --error=logs/xct_cal_%A_%a.err
set -uo pipefail
source "$(dirname "$(dirname "${BASH_SOURCE[0]}")")/env.sh"
cd "${XCT_ROOT}"; mkdir -p logs "${XCT_RUN}/calibration"
xct_py                      # the manifest read below needs a python

# --exclusive is deliberate: a timing taken next to someone else's job is not a timing.
DONORS=(${XCT_CAL_DONORS:-6450 6476})
MS=(1 2 4)
LANDS=(${XCT_CAL_LANDMARKS:-500 2000})
nM=${#MS[@]}; nL=${#LANDS[@]}
d="${DONORS[$((SLURM_ARRAY_TASK_ID / (nM * nL)))]}"
LAND="${LANDS[$(( (SLURM_ARRAY_TASK_ID / nM) % nL ))]}"
M="${MS[$((SLURM_ARRAY_TASK_ID % nM))]}"
IN="${XCT_RUN}/donor_id=${d}/roi=panc/input"
OUT="${XCT_RUN}/calibration/donor_id=${d}/landmarks=${LAND}"
mkdir -p "${OUT}"
echo "[cal] donor=${d} M=${M} landmarks=${LAND} cores=${SLURM_CPUS_PER_TASK}"

xct_R
/usr/bin/time -v Rscript "${XCT_HPG}/kodama_shard.R" \
  --input "${IN}" --output "${OUT}" --arm nonspatial \
  --m "${M}" --shard 0 --base-seed 0 --landmarks "${LAND}" \
  --n-cores "${SLURM_CPUS_PER_TASK:-4}" --force 2>&1 \
  | tee "${OUT}/cal_M${M}.log"
mv "${OUT}/res_nonspatial_shard00.json" "${OUT}/cal_M${M}.json" 2>/dev/null
rm -f "${OUT}/res_nonspatial_shard00.rds"

echo "[cal] fit the line with:  python ${XCT_HPG}/fit_calibration.py --run ${XCT_RUN}"
