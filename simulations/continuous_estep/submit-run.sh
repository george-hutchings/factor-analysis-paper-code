#!/bin/bash
#SBATCH --job-name=sim-estep
#SBATCH --constraint="skylake"
#SBATCH --partition=short
#SBATCH --array=0-2499
#SBATCH --output=logs/%x-%A_%4a.out
#SBATCH --error=logs/%x-%A_%4a.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=10:00:00
#SBATCH --mem=2GB

# --time is 10h against a measured worst case far below it: M = 400, the top arm, was timed
# at about 0.005 min per iteration, so four rungs stalled at maxiter = 2000 is roughly 40
# minutes. The headroom is because a stalled rung is the failure mode this study exists to
# measure, cost is not monotone in M, and a task killed on the limit has to be found and
# requeued.

module load R-bundle-CRAN/2023.12-foss-2023a
Rscript run.R
