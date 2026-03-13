#!/bin/bash
#SBATCH --job-name=sim-mixed
#SBATCH --partition=general
#SBATCH --array=0-2999
#SBATCH --output=logs/%x-%A_%4a.out
#SBATCH --error=logs/%x-%A_%4a.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=02:00:00
#SBATCH --mem=2G

module load george-bundle/0.1-foss-2023a-R-4.3.2
Rscript run.R
