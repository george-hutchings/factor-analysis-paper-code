#!/bin/bash
#SBATCH --job-name=cox_surv
#SBATCH --partition=general
#SBATCH --output=logs/cox_surv_%j.out
#SBATCH --error=logs/cox_surv_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G

module load george-bundle/0.1-foss-2023a-R-4.3.2

cd /home/uu85g9/factor-analysis/noms-analysis/mri/clinical_analysis
Rscript analyse_cox_survival.R
