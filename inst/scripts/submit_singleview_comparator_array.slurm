#!/bin/bash
#SBATCH --job-name=Zentangler_SV_Compare
#SBATCH --partition=scu-cpu
#SBATCH --output=Zentangler_SV_Compare_%A_%a.out
#SBATCH --error=Zentangler_SV_Compare_%A_%a.err
#SBATCH --time=150:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

source /athena/mallicklab/scratch/naa4050/miniconda3/etc/profile.d/conda.sh
conda activate r-env

PACKAGE_DIR=${PACKAGE_DIR:-/home/fs01/naa4050/Zentangler}
GRID_FILE=${GRID_FILE:-/home/fs01/naa4050/Zentangler/intersim_zentangler_grid.csv}
OUT_DIR=${OUT_DIR:-/home/fs01/naa4050/Zentangler/singleview_comparator_results}
TASK_OFFSET=${TASK_OFFSET:-0}
METHODS=${METHODS:-hima,hima2,multimedia}
HIMA_PENALTY=${HIMA_PENALTY:-DBlasso}
HIMA2_PENALTY=${HIMA2_PENALTY:-DBlasso}
MULTIMEDIA_MODEL=${MULTIMEDIA_MODEL:-lm}
MULTIMEDIA_BOOTSTRAP=${MULTIMEDIA_BOOTSTRAP:-100}
MULTIMEDIA_RF_TREES=${MULTIMEDIA_RF_TREES:-500}

which Rscript
Rscript --version

echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "TASK_OFFSET=${TASK_OFFSET}"
echo "Effective grid row=$((SLURM_ARRAY_TASK_ID + TASK_OFFSET))"
echo "GRID_FILE=${GRID_FILE}"
echo "OUT_DIR=${OUT_DIR}"
echo "PACKAGE_DIR=${PACKAGE_DIR}"
echo "METHODS=${METHODS}"

mkdir -p "${OUT_DIR}"

Rscript "${PACKAGE_DIR}/inst/scripts/run_singleview_comparator_hpc_task.R" \
  --grid "${GRID_FILE}" \
  --task-id "${SLURM_ARRAY_TASK_ID}" \
  --task-offset "${TASK_OFFSET}" \
  --out-dir "${OUT_DIR}" \
  --package-dir "${PACKAGE_DIR}" \
  --methods "${METHODS}" \
  --hima-penalty "${HIMA_PENALTY}" \
  --hima2-penalty "${HIMA2_PENALTY}" \
  --multimedia-model "${MULTIMEDIA_MODEL}" \
  --multimedia-bootstrap "${MULTIMEDIA_BOOTSTRAP}" \
  --multimedia-rf-trees "${MULTIMEDIA_RF_TREES}"
