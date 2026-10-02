#!/bin/bash
#SBATCH --job-name=analyse-cont-hd
#SBATCH --constraint="skylake"
#SBATCH --partition=short
#SBATCH --output=logs/%x-%A.out
#SBATCH --error=logs/%x-%A.err
#SBATCH --time=01:00:00
#SBATCH --mem=64GB

# 64GB, not the 8GB the low-dimensional scenarios use. analyse.R holds every
# realisation in memory at once: evaluate_run keeps four 1956 x 20 list-columns
# per row (lambda_hat, lambda_raw, w_prob, lambda_true), about 1.25 MB a row,
# so 1500 rows is roughly 2 GB before rbindlist's copy and the elementwise
# Reduce/Map over 500 matrices per method.

module load R-bundle-CRAN/2023.12-foss-2023a
Rscript ~/factor-analysis/simulations/analyse.R /well/nichols-nvs/users/peo100/factor-analysis/paper-code/simulations/continuous_highdim/
