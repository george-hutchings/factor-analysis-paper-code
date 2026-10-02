#!/bin/bash
#SBATCH --job-name=ana-estep
#SBATCH --partition=short
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem=16GB
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err

module load R-bundle-CRAN/2023.12-foss-2023a
Rscript analyse.R /well/nichols-nvs/users/peo100/factor-analysis/paper-code/simulations/continuous_estep/
