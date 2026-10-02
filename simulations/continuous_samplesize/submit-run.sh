#!/bin/bash
#SBATCH --job-name=sim-ss-final
#SBATCH --constraint="skylake"
#SBATCH --partition=short
#SBATCH --array=0-2999
#SBATCH --output=logs/%x-%A_%4a.out
#SBATCH --error=logs/%x-%A_%4a.err
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --time=03:00:00
#SBATCH --mem=4GB

# BMRC. The paper's sample-size study, at full scale:
#   6 sample sizes x 500 realisations = 3000 tasks, about 300 core-hours.
# run.R decodes the index as
#   N = c(50, 250, 500, 1000, 2500, 5000)[id %% 6 + 1], r = id %/% 6 + 1
#
#   sbatch submit-run.sh
#   sbatch --dependency=afterany:<run> submit-analyse.sh
#
# SLURM opens the log files before the job runs, so logs/ must exist. It is kept in the
# repository by logs/.gitkeep.
#
# --time is 3h against a worst cell of about 18 minutes (N = 5000, 400 iterations).
# The headroom is for the shared-filesystem stalls that can kill a task at an I/O boundary
# during array ramp-up, not for the compute.
#
# A fit that fails is caught by run.R, and analyse.R reports the number of failures.

module load R-bundle-CRAN/2023.12-foss-2023a

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

Rscript run.R
