#!/bin/bash
#SBATCH --job-name=analyse-ss-final
#SBATCH --constraint="skylake"
#SBATCH --partition=short
#SBATCH --output=logs/%x-%A.out
#SBATCH --error=logs/%x-%A.err
#SBATCH --time=01:00:00
#SBATCH --mem=16GB

# Run once, after the array. analyse.R walks the n<N>/ subdirectories itself and writes
# all_results.rds and elementwise_mse.rds, which plot_elementwise.R then reads, so the order
# matters. plot_elementwise.R is cheap and reads only what analyse.R wrote, so it can also be
# rerun on its own without resubmitting the analysis.
#
# 16GB rather than 8: this holds 3000 runs' worth of loading matrices in memory at once.

module load R-bundle-CRAN/2023.12-foss-2023a

BASE=/well/nichols-nvs/users/peo100/factor-analysis/paper-code/simulations/continuous_samplesize/

Rscript analyse.R          "$BASE"
Rscript plot_elementwise.R "$BASE"
