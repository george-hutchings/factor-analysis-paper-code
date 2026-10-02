#!/bin/bash
#SBATCH --job-name=fa-clin-habib-v2-px
#SBATCH --partition=general,legacy
#SBATCH --output=logs/%x-%A_%4a.out
#SBATCH --error=logs/%x-%A_%4a.err
#SBATCH --array=1-10
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1

module load george-bundle/0.1-foss-2023a-R-4.3.2

mkdir -p logs
Rscript run.R
