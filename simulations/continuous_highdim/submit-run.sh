#!/bin/bash
#SBATCH --job-name=sim-cont-hd
#SBATCH --constraint="skylake"
#SBATCH --partition=short
#SBATCH --array=0-1499
#SBATCH --output=logs/%x-%A_%4a.out
#SBATCH --error=logs/%x-%A_%4a.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --mem=8GB

# BMRC. 3 methods x 500 realisations = 1500 tasks, one per (method, realisation);
# run.R maps the index as slurm_id %% 3 for the method and slurm_id / 3 for the
# realisation. Fits go to /well (see save_dir in run.R).
# Create logs/ before the first submission: SLURM opens the log files before the
# job script runs.

module load R-bundle-CRAN/2023.12-foss-2023a

# Single-threaded on purpose: 1-CPU jobs schedule far more readily than an 8-CPU
# request across a 1500-task array. D = 1956 makes the M-step BLAS-bound, so pin
# the threads explicitly rather than letting BLAS oversubscribe the one core.
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

Rscript run.R
