#!/bin/bash
#SBATCH --job-name=sim-bin
#SBATCH --constraint="skylake"
#SBATCH --partition=short
#SBATCH --account=nichols-nvs.prj.low
#SBATCH --array=0-1999
#SBATCH --output=logs/%x-%A_%4a.out
#SBATCH --error=logs/%x-%A_%4a.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --mem=8GB

module load R-bundle-CRAN/2023.12-foss-2023a
Rscript run.R
