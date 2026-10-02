#!/bin/bash
# make_pca_supp_figure.sh
#
# Reproducibly builds the PCA loading maps for the first several principal
# components of the NO.MS T2 lesion data (Figures 5, S8 and S9). The whole
# pipeline is scripted and deterministic.
#
#   1. make_pca_loading_maps.R : PCA via the (n x n) gram matrix, first NPC
#      loading vectors mapped into the mask -> pca_loadings.nii.gz (4D).
#   2. fsleyes (headless) renders each PC as its OWN axial-lightbox figure over
#      the MNI152 1mm brain-only template (pca1.png .. pcaNPC.png) - one per PC.
#   3. ImageMagick autocrops the black border off each figure.
#
# The loadings are shown UNthresholded on purpose: the figure's point is that
# PCA components are dense and non-localised, unlike the sparse MoCos in Fig 4.
#
# Requirements (HPC, compute/interactive node - NOT the login node):
#   module load george-bundle/0.1-foss-2023a-R-4.3.2   (R + RNifti)
#   FSL (fsleyes, pngappend) is already on the PATH at /apps/fsl - no module needed
#
# Usage: bash make_pca_supp_figure.sh

set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"

module load george-bundle/0.1-foss-2023a-R-4.3.2
# FSL (fsleyes, pngappend) is already on the PATH at /apps/fsl - no module load needed

# ---- parameters (the only two knobs are NPC and DOSCALE) ----
export DATA="/data/users/uu85g9/factor_analysis-oldold/consolodation_may-24/processed_fa_scans.rds"
export MASK="/data/ms/processed/mri/MS_Share/4George/final_mask.nii.gz"
export OUTDIR="/data/users/uu85g9/factor-analysis/noms/mri/pca_supp"
export NPC=3
export DOSCALE=TRUE        # repo PCA convention (scale columns); set FALSE for raw
LB_NCOLS=6                 # columns of axial slices per PC (the paper uses 6)
ZLO=0.32                   # z-range as a fraction of the axial extent; trims the very
ZHI=0.69                   #   bottom (inferior) and top (vertex) slices of the brain
SS=0.031                   # slice spacing (fraction); SMALLER = more slices.
                           #   ~ (ZHI-ZLO)/11 gives about 12 slices -> two rows of 6

mkdir -p "$OUTDIR"

# ---- 1. PCA loadings -> 4D NIfTI (deterministic) ----
Rscript "$HERE/make_pca_loading_maps.R"

# ---- 2. render each PC as its OWN axial-lightbox figure over the MNI template ----
#   House style of Fig 5: a 2 x 6 grid of axial slices on a
#   grayscale MNI underlay, red-yellow (+) / blue (-) overlay, loadings shown
#   UNthresholded (dense, non-localised - that's the point). One figure per PC.
#   Tuning: SS sets the slice COUNT (smaller = more slices); ZLO/ZHI set the
#   RANGE and trim the very bottom/top of the brain. Units are fractions (0-1).
NIFTI="$OUTDIR/pca_loadings.nii.gz"
TEMPLATE="${FSLDIR}/data/standard/MNI152_T1_1mm_brain.nii.gz"   # skull-stripped 1mm underlay
for ((k = 0; k < NPC; k++)); do
  fsleyes render --scene lightbox --zaxis 2 --hideCursor \
    --ncols "$LB_NCOLS" --zrange "$ZLO" "$ZHI" --sliceSpacing "$SS" \
    --size 1800 700 --outfile "$OUTDIR/pca$((k + 1)).png" \
    "$TEMPLATE" \
    "$NIFTI" --volume "$k" \
    --cmap red-yellow --negativeCmap blue-lightblue --useNegativeCmap
done

# ---- 3. autocrop the black border off each figure (ImageMagick, if available) ----
CROP=""
command -v magick  >/dev/null 2>&1 && CROP="magick"
command -v convert >/dev/null 2>&1 && CROP="convert"
if [ -n "$CROP" ]; then
  for ((k = 1; k <= NPC; k++)); do
    "$CROP" "$OUTDIR/pca$k.png" -fuzz 2% -trim +repage "$OUTDIR/pca$k.png"
  done
else
  echo "NOTE: ImageMagick (convert/magick) not found - skipping autocrop; try 'module load ImageMagick'."
fi

echo "Done. Three separate figures, one per PC:"
for ((k = 1; k <= NPC; k++)); do echo "  $OUTDIR/pca$k.png"; done
echo "Variance explained per PC: $OUTDIR/pca_variance_explained.txt"
echo "scp pca1.png pca2.png pca3.png into the paper repo figures/ for the supplement."
