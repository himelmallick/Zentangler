#!/bin/bash
#SBATCH --job-name=Zentangler_Paired
#SBATCH --partition=scu-cpu
#SBATCH --output=Zentangler_Paired_%A_%a.out
#SBATCH --error=Zentangler_Paired_%A_%a.err
#SBATCH --time=150:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

source /athena/mallicklab/scratch/naa4050/miniconda3/etc/profile.d/conda.sh
conda activate r-env

PACKAGE_DIR=${PACKAGE_DIR:-/home/fs01/naa4050/Zentangler}
GRID_FILE=${GRID_FILE:-/home/fs01/naa4050/Zentangler/intersim_zentangler_grid.csv}
SIM_DIR=${SIM_DIR:-/home/fs01/naa4050/Zentangler/simmba_cache}
OUT_DIR=${OUT_DIR:-/home/fs01/naa4050/Zentangler/zentangler_paired_results}
TASK_OFFSET=${TASK_OFFSET:-0}
TASK_SCRIPT=${TASK_SCRIPT:-${PACKAGE_DIR}/run_intersim_hpc_task.R}

which Rscript
Rscript --version

echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "TASK_OFFSET=${TASK_OFFSET}"
echo "Effective grid row=$((SLURM_ARRAY_TASK_ID + TASK_OFFSET))"
echo "GRID_FILE=${GRID_FILE}"
echo "SIM_DIR=${SIM_DIR}"
echo "OUT_DIR=${OUT_DIR}"
echo "PACKAGE_DIR=${PACKAGE_DIR}"
echo "TASK_SCRIPT=${TASK_SCRIPT}"

mkdir -p "${OUT_DIR}"

Rscript "${TASK_SCRIPT}" \
  --grid "${GRID_FILE}" \
  --task-id "${SLURM_ARRAY_TASK_ID}" \
  --task-offset "${TASK_OFFSET}" \
  --out-dir "${OUT_DIR}" \
  --package-dir "${PACKAGE_DIR}" \
  --sim-dir "${SIM_DIR}"
