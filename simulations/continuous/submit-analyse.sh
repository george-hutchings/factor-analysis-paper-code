#!/bin/bash
#SBATCH --job-name=analyse-cont
#SBATCH --constraint="skylake"
#SBATCH --partition=short
#SBATCH --output=logs/%x-%A.out
#SBATCH --error=logs/%x-%A.err
#SBATCH --time=00:30:00
#SBATCH --mem=4GB

module load R-bundle-CRAN/2023.12-foss-2023a
Rscript ~/factor-analysis/simulations/analyse.R /well/nichols-nvs/users/peo100/factor-analysis/simulations/continuous/  # existing paper results in simulations_new/continuous/
