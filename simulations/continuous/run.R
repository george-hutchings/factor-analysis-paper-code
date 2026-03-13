### Simulation Scenario 1: Continuous Data
### SLURM array job: each task runs one (method, realisation) pair
# Note: Legacy simulation data and methods are available in older git commits under simulations/comparison/

# ---- CONFIG (BMRC cluster) ----
algo_path = "~/factor-analysis/mcem_fa_algorithm.R"
save_dir  = "/well/nichols-nvs/users/peo100/factor-analysis/simulations/continuous/"  # existing paper results stored in simulations_new/continuous/
# ----------------

slurm_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

method_names = c(
  'our method',
  'our method long v0',
  'rockova 2016',
  'factanal K known',
  'pca K known'
)
n_methods = length(method_names)

slurm_id_modn = slurm_id %% n_methods + 1
r = floor(slurm_id / n_methods) + 1

# ---- True loading matrix ----
lambda = matrix(0, 12, 3)
lambda[1:6, 1] = 1
lambda[7:10, 2] = 1
lambda[11:12, 3] = 1
true_K = ncol(lambda)
D = nrow(lambda)

# ---- Generate data ----
set.seed(12345 + r)
N = 200
eta = matrix(rnorm(N * true_K), N, true_K)
Z = tcrossprod(eta, lambda) + matrix(rnorm(N * D), N, D)
Y = Z  # continuous: observed = latent

# ---- Method selection ----
K = 10
v0s = NULL

# Load baselines (needed for methods 3-5)
source("~/factor-analysis/simulations/methods/pca_factanal_nmf.R")

method_name = method_names[slurm_id_modn]

if (slurm_id_modn == 1) {
  # Our method - short DPE
  source(algo_path)
  v0s = 10^-seq(1, 4, length.out = 4)
} else if (slurm_id_modn == 2) {
  # Our method - long DPE
  source(algo_path)
  v0s = c(0.02, 0.018, 0.015, 0.01, 0.005, 0.001, 0.0005, 0.0001)
} else if (slurm_id_modn == 3) {
  # Rockova & George (2016)
  source("~/factor-analysis/simulations/methods/veronika2016_wrapper.R")
  dpe_func = function(v0s, ...) { list(veronika2016(...)) }
} else if (slurm_id_modn == 4) {
  # factanal (K known)
  K = ncol(lambda)
  dpe_func = function(v0s, ...) {
    lambda_corr = my_factanal(...)
    list(list(w = lambda_corr * NA, lambda = lambda_corr * NA, lambda_corr = lambda_corr))
  }
} else if (slurm_id_modn == 5) {
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

if (slurm_id_modn <= 2) {
  # Our method: explicit hyperparameters
  results = dpe_func(v0s, Y = Y, K = K, do_stick_breaking = TRUE,
                     px_sigma = FALSE, px_rotate = FALSE, varimax_every = -1,
                     conv_param = TRUE, nburn = 100, n_mcmc = 100, tol=0.016)
} else {
  results = dpe_func(v0s, Y = Y, K = K, do_stick_breaking = TRUE)
}

time_taken = difftime(Sys.time(), tmp, units = 'mins')
print(paste('Time taken:', time_taken, 'mins'))

# ---- Evaluate ----
source("~/factor-analysis/simulations/shared/evaluate.R")
results_df = evaluate_run(results, lambda, method_name, r, time_taken)
saveRDS(results_df, paste0(save_dir, slurm_id, '.rds'))
