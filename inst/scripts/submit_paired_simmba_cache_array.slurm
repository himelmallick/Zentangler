#!/bin/bash
#SBATCH --job-name=Zentangler_SimCache
#SBATCH --partition=scu-cpu
#SBATCH --output=Zentangler_SimCache_%A_%a.out
#SBATCH --error=Zentangler_SimCache_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1

source /athena/mallicklab/scratch/naa4050/miniconda3/etc/profile.d/conda.sh
conda activate r-env

PACKAGE_DIR=${PACKAGE_DIR:-/home/fs01/naa4050/Zentangler}
GRID_FILE=${GRID_FILE:-/home/fs01/naa4050/Zentangler/intersim_zentangler_grid.csv}
SIM_DIR=${SIM_DIR:-/home/fs01/naa4050/Zentangler/simmba_cache}
TASK_OFFSET=${TASK_OFFSET:-0}
TASK_SCRIPT=${TASK_SCRIPT:-${PACKAGE_DIR}/make_paired_simmba_cache_hpc_task.R}
OVERWRITE=${OVERWRITE:-FALSE}

which Rscript
Rscript --version

echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "TASK_OFFSET=${TASK_OFFSET}"
echo "Effective grid row=$((SLURM_ARRAY_TASK_ID + TASK_OFFSET))"
echo "GRID_FILE=${GRID_FILE}"
echo "SIM_DIR=${SIM_DIR}"
echo "PACKAGE_DIR=${PACKAGE_DIR}"
echo "TASK_SCRIPT=${TASK_SCRIPT}"

mkdir -p "${SIM_DIR}"

Rscript "${TASK_SCRIPT}" \
  --grid "${GRID_FILE}" \
  --task-id "${SLURM_ARRAY_TASK_ID}" \
  --task-offset "${TASK_OFFSET}" \
  --out-dir "${SIM_DIR}" \
  --package-dir "${PACKAGE_DIR}" \
  --overwrite "${OVERWRITE}"
