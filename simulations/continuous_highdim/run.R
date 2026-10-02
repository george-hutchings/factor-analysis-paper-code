### Simulation Scenario 1 at high dimension (Rockova & George 2016 dimensions).
### SLURM array job: each task runs one (method, realisation) pair.
### n = 100, G = 1956, K_true = 5 overlapping factors, fit with K = 20.
###
### factanal is omitted: with p > n the sample correlation matrix is singular and
### stats::factanal stops before it iterates (measured rcond 2e-20 at p = 200,
### 8e-22 at p = 1956, both n = 100).

# ---- CONFIG (BMRC cluster) ----
# SIM_OUT overrides save_dir, for a local test run against a scratch directory.
algo_path = "~/factor-analysis/mcem_fa_algorithm.R"
save_dir  = Sys.getenv(
  "SIM_OUT",
  "/well/nichols-nvs/users/peo100/factor-analysis/paper-code/simulations/continuous_highdim/")
# ----------------

slurm_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

# saveRDS is the last line of the script, so a missing save_dir would only bite
# after the fit has already run. Create it up front instead.
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

method_names = c(
  'our method',
  'rockova 2016',
  'pca K known'
)
n_methods = length(method_names)

slurm_id_modn = slurm_id %% n_methods + 1
r = floor(slurm_id / n_methods) + 1

# ---- True loading matrix: R&G 2016 synthetic example ----
# G = 1956, K_true = 5, each factor 500 nonzeros, adjacent factors overlap by 136
# (offset 135 -> inclusive overlap 500 - 365 + 1 = 136), loadings = 1, residual = I.
G = 1956
dd = 5
blk = 500
offset = 135
lambda = matrix(0, G, dd)
end = 1
for (i in 1:dd) {
  start = end - (i != 1) * offset
  end = start + blk - 1
  lambda[start:end, i] = 1
}
true_K = ncol(lambda)
D = nrow(lambda)

# ---- Generate data ----
set.seed(12345 + r)
N = 100
eta = matrix(rnorm(N * true_K), N, true_K)
Z = tcrossprod(eta, lambda) + matrix(rnorm(N * D), N, D)
Y = Z  # continuous: observed = latent

# ---- Method selection ----
K = 20
v0s = NULL

source("~/factor-analysis/simulations/methods/pca_factanal_nmf.R")

method_name = method_names[slurm_id_modn]

if (slurm_id_modn == 1) {
  # Our method. Ladder chosen empirically at this dimension. Two settings differ
  # from the low-dimensional scenarios:
  #
  #   px_sigma and px_rotate are ON. Without them the fitted loadings settle near
  #   half their true magnitude here; with them they reach it.
  #
  #   tol = 0.045, not the low-dimensional 0.016. conv_param compares
  #   max|abs(lambda) - abs(lambda_old)| over all 1956 x 20 loadings, against
  #   12 x 10 low-dimensional. Measured, that statistic stays in 0.041-0.050 here,
  #   so 0.016 is unreachable and every rung would run to maxiter.
  source(algo_path)
  v0s = c(0.1, 0.075, 0.05)
} else if (slurm_id_modn == 2) {
  # Rockova & George (2016)
  source("~/factor-analysis/simulations/methods/veronika2016_wrapper.R")
  # varimax_all=TRUE here, unlike the low-dimensional scenario. FACTOR_ROTATE's
  # varimax argument is a period, so TRUE rotates after every EM iteration.
  # Without it the lambda0=20 stage did not terminate here, and with it support
  # recovery was better (FP+FN 37 against 302 over 10 realisations).
  dpe_func = function(v0s, ...) { list(veronika2016(..., varimax_all = TRUE)) }
} else if (slurm_id_modn == 3) {
  # PCA + varimax (K known)
  K = ncol(lambda)
  dpe_func = function(v0s, ...) {
    lambda_corr = pca_varimax(...)
    list(list(w = lambda_corr * NA, lambda = lambda_corr * NA, lambda_corr = lambda_corr))
  }
}

# ---- Run ----
print(paste("Method:", method_name, "| Realisation:", r))
set.seed(12345)
tmp = Sys.time()

if (slurm_id_modn == 1) {
  # Our method: explicit hyperparameters. maxiter caps a rung that never dips
  # below tol, bounding the single-threaded walltime. DPE warm-starts each rung
  # from the previous one, so a capped rung still passes its state on.
  results = dpe_func(v0s, Y = Y, K = K, do_stick_breaking = TRUE,
                     px_sigma = TRUE, px_rotate = TRUE, varimax_every = -1,
                     conv_param = TRUE, nburn = 100, n_mcmc = 100, tol = 0.045,
                     maxiter = 100)
} else {
  results = dpe_func(v0s, Y = Y, K = K, do_stick_breaking = TRUE)
}

time_taken = difftime(Sys.time(), tmp, units = 'mins')
print(paste('Time taken:', time_taken, 'mins'))

# ---- Evaluate ----
source("~/factor-analysis/simulations/shared/evaluate_highdim.R")
results_df = evaluate_run(results, lambda, method_name, r, time_taken)
saveRDS(results_df, paste0(save_dir, slurm_id, '.rds'))
