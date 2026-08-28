#!/bin/bash
# All-in-one benchmark job (prep -> methods -> eval) on a single CPU node.
# Recommended for the headline crop50k. KODAMA is the slow method; give time.
#   sbatch scripts/benchmark/slurm/bench.sh
#   SAMPLES="0041435" VARIANTS="crop50k" sbatch --export=ALL scripts/benchmark/slurm/bench.sh
#SBATCH --job-name=spatial_bench
#SBATCH --account=maigan
#SBATCH --qos=maigan-b
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128gb
#SBATCH --time=24:00:00
#SBATCH --output=logs/benchmark/slurm_bench_%j.out
#SBATCH --error=logs/benchmark/slurm_bench_%j.err
set -uo pipefail
cd "${XENIUM_ANALYSIS_ROOT:?set XENIUM_ANALYSIS_ROOT on the cluster}"
mkdir -p logs/benchmark
export NCORES="${SLURM_CPUS_PER_TASK:-16}"
SAMPLES="${SAMPLES:-0041323 0041435}" \
VARIANTS="${VARIANTS:-crop50k}" \
STEPS="${STEPS:-prep methods eval}" \
GENE_PANEL="${GENE_PANEL:-embedding_panel}" \
  bash scripts/benchmark/run_benchmark.sh
