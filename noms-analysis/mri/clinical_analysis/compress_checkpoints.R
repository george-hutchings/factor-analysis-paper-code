# Compress specific checkpoints from runs_jan18_v0_0.1 experiments
# Outputs Eeta (expected latent scores) for each subject
#
# Uses mri-clinical_processed.rds so subjects match clinical-only file for downstream joins

library(Rcpp)
library(RcppArmadillo)

# Source C++ E-step (uses correct eta sampling: L*eps where LL'=Sigma)
sourceCpp('~/factor-analysis/src/mc_e_step.cpp')

# === Compression function with rank-based Z init ===
compress_data <- function(Y, lambda, alpha, nburn = 1000, nsamples = 250) {
  N <- nrow(Y)
  D <- ncol(Y)

  # Identify binary columns and remap Y to ordinal categories 1,2,...
  binary <- rep(FALSE, D)
  for (d in seq_len(D)) {
    uniq_Y <- sort(unique(Y[, d]))  # ignores NA
    if (length(uniq_Y) == 2) {
      binary[d] <- TRUE
    }
    Y[, d] <- match(Y[, d], uniq_Y)  # remap to 1,2,... (ignores NA)
  }
  binary_rcpp <- as.integer(binary)

  # Group indices: which rows belong to each category, per column
  group_indices <- lapply(seq_len(D), function(d) split(seq_len(N), Y[, d]))

  # Rank-based Z initialization
  rank_with_na <- function(x) {
    tmp <- rank(x, na.last = 'keep')
    tmp[is.na(tmp)] <- sample(tmp[!is.na(tmp)], sum(is.na(tmp)), replace = TRUE)
    tmp <- rank(tmp, ties.method = "average")
    return(tmp)
  }
  tmp <- apply(Y, 2, rank_with_na) / (N + 1)
  Z <- qnorm(tmp)
  Z <- scale(Z)

  EZdTZd <- colSums(Z^2)

  # Burn-in
  cat(sprintf('Burn-in (%d iters, rank-based Z init)...\n', nburn))
  t_start <- Sys.time()
  mce_res <- mce_step_rcpp(
    Y = Y, lambda = lambda, Z0 = Z, alpha = alpha,
    group_indices = group_indices, EZdTZd0 = EZdTZd,
    binary = binary_rcpp, nsteps = nburn
  )
  Z <- mce_res$Z
  EZdTZd <- mce_res$EZdTZd
  cat(sprintf('Burn-in done: %.1f min\n', as.numeric(difftime(Sys.time(), t_start, units = 'mins'))))

  # Sampling to compute E[eta]
  cat(sprintf('Sampling (%d iters)...\n', nsamples))
  t_start <- Sys.time()
  mce_res <- mce_step_rcpp(
    Y = Y, lambda = lambda, Z0 = Z, alpha = alpha,
    group_indices = group_indices, EZdTZd0 = EZdTZd,
    binary = binary_rcpp, nsteps = nsamples
  )
  Eeta <- mce_res$Eeta
  cat(sprintf('Sampling done: %.1f min\n', as.numeric(difftime(Sys.time(), t_start, units = 'mins'))))

  return(Eeta)
}

# === Configuration ===
base_path <- "/data/users/uu85g9/factor-analysis/noms/mri"  # existing paper results stored under runs_jan18_v0_0.1/
dir.create(file.path(base_path, "clinical_analysis"), showWarnings = FALSE)
output_dir <- file.path(base_path, "clinical_analysis/compressed_data")
dir.create(output_dir, showWarnings = FALSE)

# Load MRI-clinical data and extract MRI-only columns
# Using mri-clinical_processed.rds so rows match mri-clinical_clinical-only_processed.rds (for USUBJID joins)
tmp <- readRDS("/data/users/uu85g9/factor-analysis-old/mri-clinical/mri-clinical_processed.rds")
Y <- as.matrix(tmp[, grepl("^X\\d+", colnames(tmp))])
cat(sprintf("Loaded Y: %d subjects x %d voxels\n", nrow(Y), ncol(Y)))

# Define checkpoints to compress (indexed 1-4 for array job)
checkpoints <- list(
  list(
    name = "restart_dpe3_dpe1_iter20",
    condition = "dpe_random_pxsigma_restart_dpe3",
    dpe_stage = 1,
    iter = 20
  ),
  list(
    name = "restart_dpe3_dpe1_iter70",
    condition = "dpe_random_pxsigma_restart_dpe3",
    dpe_stage = 1,
    iter = 70
  ),
  list(
    name = "random_pxsigma_dpe2_iter70",
    condition = "dpe_random_pxsigma",
    dpe_stage = 2,
    iter = 70
  ),
  list(
    name = "pca_nopx_dpe3_iter90",
    condition = "dpe_pca_nopx",
    dpe_stage = 3,
    iter = 90
  )
)

# === Get checkpoint from SLURM array task ID (1-4) ===
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
cp <- checkpoints[[task_id]]

cat(sprintf("Array task %d: %s\n", task_id, cp$name))

# Load checkpoint
checkpoint_file <- file.path(
  base_path, cp$condition, "intermediate_results",
  sprintf("dpe%d_iter%d.rds", cp$dpe_stage, cp$iter)
)
checkpoint <- readRDS(checkpoint_file)

# Extract parameters
lambda <- checkpoint$lambda
alpha <- checkpoint$alpha

cat(sprintf("Lambda dims: %d x %d\n", nrow(lambda), ncol(lambda)))
cat(sprintf("Y dims: %d x %d\n", nrow(Y), ncol(Y)))
cat(sprintf("Active factors (w > 0.5): %d\n", sum(colSums(checkpoint$w > 0.5) > 0)))

# Compress with rank-based Z init (longer burn-in since not using checkpoint Z)
Eeta <- compress_data(Y, lambda, alpha, nburn = 1000, nsamples = 250)

# Save Eeta
output_file <- file.path(output_dir, sprintf("Eeta_%s.rds", cp$name))
saveRDS(Eeta, output_file)
cat(sprintf("Saved: %s\n", output_file))

# Save metadata
meta <- list(
  checkpoint_file = checkpoint_file,
  condition = cp$condition,
  dpe_stage = cp$dpe_stage,
  iter = cp$iter,
  n_subjects = nrow(Eeta),
  n_factors = ncol(Eeta),
  active_factors = sum(colSums(checkpoint$w > 0.5) > 0),
  w = checkpoint$w
)
saveRDS(meta, file.path(output_dir, sprintf("meta_%s.rds", cp$name)))

cat("\nCompression complete!\n")
