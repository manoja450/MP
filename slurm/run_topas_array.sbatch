#!/bin/bash
# Run one TOPAS job per array task. Usage:
#   sbatch --array=1-N slurm/run_topas_array.sbatch /path/to/jobs.txt
# Adjust partition/time/memory and the module lines for YOUR cluster.
#SBATCH --job-name=topas
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=12:00:00
#SBATCH --output=slurm-%A_%a.out
##SBATCH --partition=YOUR_PARTITION
##SBATCH --account=YOUR_ACCOUNT

set -euo pipefail
JOBS_FILE="${1:?give the jobs.txt file as the first argument}"
RUN_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$JOBS_FILE")
[ -n "$RUN_FILE" ] || { echo "no line ${SLURM_ARRAY_TASK_ID} in $JOBS_FILE"; exit 1; }

# --- environment: edit for your cluster ---------------------------------
# module load gcc cmake
# export PATH=$HOME/Applications/TOPAS/OpenTOPAS-install/bin:$PATH
# export TOPAS_G4_DATA_DIR=$HOME/Applications/GEANT4/G4DATA
# ------------------------------------------------------------------------

JOB_DIR=$(dirname "$RUN_FILE")
cd "$JOB_DIR"
if [ -f done.flag ]; then echo "already done: $JOB_DIR"; exit 0; fi
echo "host=$(hostname) start=$(date -Is) run=$RUN_FILE"
/usr/bin/time -v topas run.txt > topas.log 2>&1
touch done.flag
echo "end=$(date -Is)"
