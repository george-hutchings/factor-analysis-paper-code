# === EDIT THESE ===
base_path <- "/data/users/uu85g9/factor-analysis/noms/mri"  # previous: .../mri/runs_jan18_v0_0.1
condition <- "dpe_random_pxsigma_restart_dpe3"  # paper model (previously: dpe_pca_varimax_rotate_pxsigma)
iter <- -1  # -1 = res.rds, otherwise intermediate checkpoint
dpe_stage <- -1  # -1 = final stage, 1-4 = specific DPE stage (restart_dpe3 has stages 1-2 only)

# Load result
if (iter == -1) {
  res <- readRDS(paste0(base_path, "/", condition, "/res.rds"))
} else {
  prefix <- ""
  if (dpe_stage > 0) prefix <- paste0("dpe", dpe_stage, "_")
  res <- readRDS(paste0(base_path, "/", condition, "/intermediate_results/", prefix, "iter", iter, ".rds"))
}

# Extract DPE stage if needed (dpe_* conditions return a list of results)
dpe_stage_used <- dpe_stage
if (is.list(res) && is.list(res[[1]]) && "lambda" %in% names(res[[1]])) {
  if (dpe_stage == -1) dpe_stage_used <- length(res)
  res <- res[[dpe_stage_used]]
}

# Lambda and thresholded version
lambda <- res$lambda
lambda_thresholded <- lambda
lambda_thresholded[res$w < 0.5] <- 0

# Active factor count
K <- sum(colSums(res$w > 0.5) > 0)
cat("Condition:", condition, "\n")
cat("Active factors (K):", K, "\n")

# MRI mask
library(RNifti)
mask_nifti <- readNifti("/data/ms/processed/mri/MS_Share/4George/final_mask.nii.gz")
mask <- as.logical(mask_nifti > 0)
n_voxels <- sum(mask)

# Convert lambda matrix to 4D NIfTI array
lambda_to_nifti <- function(mat) {
  if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
  arr <- array(0, c(dim(mask_nifti), ncol(mat)))
  for (i in seq_len(ncol(mat))) {
    vec <- rep(0, length(mask))
    vec[mask] <- mat[1:n_voxels, i]
    arr[, , , i] <- vec
  }
  arr
}

# Save lambda NIfTIs
out_dir <- paste0(base_path, "/", condition)
suffix <- ""
if (iter != -1) suffix <- paste0(suffix, "_iter", iter)
if (dpe_stage != -1) suffix <- paste0(suffix, "_dpe", dpe_stage)
writeNifti(lambda_to_nifti(lambda), paste0(out_dir, "/lambda_", condition, suffix, ".nii.gz"), template = mask_nifti)
writeNifti(lambda_to_nifti(lambda_thresholded), paste0(out_dir, "/lambda_thresholded_", condition, suffix, ".nii.gz"), template = mask_nifti)

# Sort thresholded lambda by sparsity (number of non-zero voxels)
non_zeros <- colSums(lambda_thresholded != 0)
order_indices <- order(non_zeros, decreasing = TRUE)
lambda_thresholded_sorted <- lambda_thresholded[, order_indices, drop = FALSE]
writeNifti(lambda_to_nifti(lambda_thresholded_sorted), paste0(out_dir, "/lambda_thresholded_sorted_", condition, suffix, ".nii.gz"), template = mask_nifti)

cat("Saved to:", out_dir, "\n")

# Create factor atlas (winner-takes-all parcellation)
# Each voxel labeled with the factor having highest |λ|, 0 if none active
abs_lambda <- abs(lambda_thresholded[1:n_voxels, ])
has_active <- rowSums(abs_lambda > 0) > 0
max_factor <- rep(0L, n_voxels)
max_factor[has_active] <- max.col(abs_lambda[has_active, , drop = FALSE], ties.method = "first")

atlas_vec <- rep(0L, length(mask))
atlas_vec[mask] <- max_factor
factor_atlas <- array(atlas_vec, dim(mask_nifti))
writeNifti(factor_atlas, paste0(out_dir, "/factor_atlas_", condition, suffix, ".nii.gz"),
           template = mask_nifti)
cat("Factor atlas saved (labels 1-", max(max_factor), ", 0=inactive)\n", sep = "")

# Create overlap count image (number of active factors per voxel)
overlap_count <- rowSums(abs_lambda > 0)
overlap_vec <- rep(0L, length(mask))
overlap_vec[mask] <- overlap_count
factor_overlap <- array(overlap_vec, dim(mask_nifti))
writeNifti(factor_overlap, paste0(out_dir, "/factor_overlap_", condition, suffix, ".nii.gz"),
           template = mask_nifti)
cat("Factor overlap saved (0-", max(overlap_count), " factors per voxel)\n", sep = "")
