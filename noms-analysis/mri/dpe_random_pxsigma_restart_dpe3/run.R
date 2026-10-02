# Restart from DPE2 iter70 checkpoint, skip directly to DPE3
# Original run: dpe_random_pxsigma with v0s = [0.1, 0.01, 0.001, 0.0001]
# This run: starts from DPE3 (v0=0.001) using DPE2 iter70 as initialization

source("/home/uu85g9/factor-analysis/mcem_fa_algorithm.R")

# Load data
Y = readRDS("/data/users/uu85g9/factor_analysis-oldold/consolodation_may-24/processed_fa_scans.rds")

# Output directory for this restart experiment
save_dir = "/data/users/uu85g9/factor-analysis/noms/mri/dpe_random_pxsigma_restart_dpe3/"
dir.create(save_dir, showWarnings=FALSE)
dir.create(paste0(save_dir, "intermediate_results/"), showWarnings=FALSE)

# Load checkpoint from DPE2 iteration 70
checkpoint_path = "/data/users/uu85g9/factor-analysis/noms/mri/dpe_random_pxsigma/intermediate_results/dpe2_iter70.rds"
checkpoint = readRDS(checkpoint_path)

# Extract parameters from checkpoint
lambda_init = checkpoint$lambda
Z_init = checkpoint$Z
alpha_init = checkpoint$alpha
thetas_init = checkpoint$thetas

# Apply varimax to lambda (matching what dpe_func does between stages)
D = nrow(lambda_init)
K = ncol(lambda_init)
lambda_init = matrix(varimax(lambda_init + 1e-6)$loadings, D, K)

# DPE schedule starting from DPE3: v0 = 0.001 -> 0.0001
v0s = c(0.001, 0.0001)

# Run DPE from stage 3 onward
res = dpe_func(
  v0s = v0s,
  Y = Y,
  K = 100,
  v1 = 10,
  do_stick_breaking = TRUE,
  tol = 0.01,
  conv_param = TRUE,
  maxiter = 2000,
  save_dir = paste0(save_dir, "intermediate_results/"),
  verbose = 0,
  lambda = lambda_init,
  Z = Z_init,
  alpha = alpha_init,
  thetas = thetas_init,
  px_rotate = FALSE,
  varimax_every = -1,
  px_sigma = TRUE,
  nburn = 100,
  n_mcmc = 50,
  save_every = 10
)

saveRDS(res, paste0(save_dir, "res.rds"))

# Also save source checkpoint info for provenance
provenance = list(
  source_checkpoint = checkpoint_path,
  source_run = "dpe_random_pxsigma",
  restart_from = "dpe2_iter70",
  original_v0s = c(0.1, 0.01, 0.001, 0.0001),
  restart_v0s = v0s
)
saveRDS(provenance, paste0(save_dir, "provenance.rds"))
