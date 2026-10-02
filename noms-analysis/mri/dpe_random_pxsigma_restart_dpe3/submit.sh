#!/bin/bash
#SBATCH --job-name=fa-mri-jan18-restart-dpe3
#SBATCH --partition=general
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --time=56-00:00:00

module load george-bundle/0.1-foss-2023a-R-4.3.2
mkdir -p logs
Rscript run.R
