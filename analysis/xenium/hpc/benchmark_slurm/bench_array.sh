#!/bin/bash
# Parallel methods step: one array task per (sample x method) for a SINGLE
# variant. Run prep first and eval after (see the sbatch chain below).
#
#   # 2 samples x 5 methods = 10 tasks (--array=0-9)
#   PREP=$(sbatch --parsable --export=ALL,STEPS=prep,VARIANTS=crop50k scripts/benchmark/slurm/bench.sh)
#   METH=$(sbatch --parsable --dependency=afterok:$PREP --export=ALL,VARIANTS=crop50k scripts/benchmark/slurm/bench_array.sh)
#   sbatch --dependency=afterok:$METH --export=ALL,STEPS=eval,VARIANTS=crop50k scripts/benchmark/slurm/bench.sh
#SBATCH --job-name=spbench_arr
#SBATCH --account=maigan
#SBATCH --qos=maigan-b
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=24:00:00
#SBATCH --array=0-9
#SBATCH --output=logs/benchmark/slurm_arr_%A_%a.out
#SBATCH --error=logs/benchmark/slurm_arr_%A_%a.err
set -uo pipefail
cd "${XENIUM_ANALYSIS_ROOT:?set XENIUM_ANALYSIS_ROOT on the cluster}"
mkdir -p logs/benchmark
export NCORES="${SLURM_CPUS_PER_TASK:-16}"

SAMPLES_ARR=(${SAMPLES:-0041323 0041435})
METHODS_ARR=(banksy utag baseline kodama misty)
VARIANT="${VARIANTS:-crop50k}"
n=${#METHODS_ARR[@]}
s="${SAMPLES_ARR[$((SLURM_ARRAY_TASK_ID / n))]}"
m="${METHODS_ARR[$((SLURM_ARRAY_TASK_ID % n))]}"
echo "array task ${SLURM_ARRAY_TASK_ID}: sample=${s} method=${m} variant=${VARIANT}"

STEPS="${m}" SAMPLES="${s}" VARIANTS="${VARIANT}" \
  bash scripts/benchmark/run_benchmark.sh
