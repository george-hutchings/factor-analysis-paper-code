#!/bin/bash
#SBATCH --job-name=compress_jan18
#SBATCH --partition=general
#SBATCH --array=1-4
#SBATCH --output=logs/compress_%A_%a.out
#SBATCH --error=logs/compress_%A_%a.err
#SBATCH --time=3-00:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4

# Array job: tasks 1-4 correspond to the 4 checkpoints
# %A = master job ID, %a = array task ID

module load george-bundle/0.1-foss-2023a-R-4.3.2

cd /home/uu85g9/factor-analysis/noms-analysis/mri/clinical_analysis

# Create logs directory if needed
mkdir -p logs

echo "Running array task $SLURM_ARRAY_TASK_ID"
Rscript compress_checkpoints.R
