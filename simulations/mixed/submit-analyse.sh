#!/bin/bash
#SBATCH --job-name=analyse-mixed
#SBATCH --partition=general
#SBATCH --output=logs/%x-%A.out
#SBATCH --error=logs/%x-%A.err
#SBATCH --time=00:30:00
#SBATCH --mem=4G

module load george-bundle/0.1-foss-2023a-R-4.3.2
Rscript ~/factor-analysis/simulations/analyse.R /data/users/uu85g9/factor-analysis/simulations/mixed/  # existing paper results in simulations_new/mixed/
