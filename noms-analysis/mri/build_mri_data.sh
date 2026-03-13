#!/bin/bash
# build_mri_data.sh
#
# Recreates processed_fa_scans.rds from individual T2 lesion maps.
#
# NO.MS has two trial arms (D2301, D2309). Each subject has a T2 lesion map
# in MNI space. This script selects baseline visits (ses-V1x), merges them
# into 4D volumes with fslmerge, applies the brain mask, and stacks both
# arms into a subjects x voxels R matrix. D2301 rows come first.
#
# Verified 2026-03-12: fslmerge output is identical to original 4D files
# (fslmaths -sub + fslstats -R = 0.0 0.0 for both arms), and the resulting
# .rds is identical() to the version used in all MRI factor analysis runs.
#
# Requirements: FSL, R with RNifti
#   On HPC: module load george-bundle/0.1-foss-2023a-R-4.3.2
#
# Usage: bash build_mri_data.sh

module load george-bundle/0.1-foss-2023a-R-4.3.2

IND_DIR="/data/ms/processed/mri/MS_Share/4George/T2lesion_ind"
MASK="/data/ms/processed/mri/MS_Share/4George/final_mask.nii.gz"
OUTPUT_DIR="/data/users/uu85g9/factor_analysis/consolodation_may-24"

# Merge baseline (ses-V1x) scans into 4D volumes, sorted alphabetically
echo "Merging D2301 baseline scans..."
fslmerge -t "$OUTPUT_DIR/FTY2301_baseline.nii.gz" \
  $(ls "$IND_DIR"/sub-CFTY720D2301*ses-V1x*.nii.gz | sort)

echo "Merging D2309 baseline scans..."
fslmerge -t "$OUTPUT_DIR/FTY2309_baseline.nii.gz" \
  $(ls "$IND_DIR"/sub-CFTY720D2309*ses-V1x*.nii.gz | sort)

# Apply mask, stack into subjects x voxels matrix, save as .rds
echo "Building R matrix..."

Rscript -e "
library(RNifti)
mask <- readNifti('$MASK') > 0
scans1 <- readNifti('$OUTPUT_DIR/FTY2301_baseline.nii.gz')
scans2 <- readNifti('$OUTPUT_DIR/FTY2309_baseline.nii.gz')
Y <- rbind(
  matrix(scans1[mask], ncol = sum(mask), byrow = TRUE),
  matrix(scans2[mask], ncol = sum(mask), byrow = TRUE)
)
cat('Final matrix:', nrow(Y), 'subjects x', ncol(Y), 'voxels\n')
saveRDS(Y, '$OUTPUT_DIR/processed_fa_scans.rds')
"

echo "Done: $OUTPUT_DIR/processed_fa_scans.rds"
