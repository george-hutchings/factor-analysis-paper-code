#!/bin/bash
#SBATCH --job-name=edss_corr
#SBATCH --partition=general
#SBATCH --output=logs/edss_corr_%j.out
#SBATCH --error=logs/edss_corr_%j.err
#SBATCH --time=01:00:00
#SBATCH --mem=32G

module load george-bundle/0.1-foss-2023a-R-4.3.2

cd /home/uu85g9/factor-analysis/noms-analysis/mri/clinical_analysis
Rscript analyse_edss_correlation.R
