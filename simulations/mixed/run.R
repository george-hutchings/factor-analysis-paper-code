### Simulation Scenario 3: Mixed Type and Missing Data
### SLURM array job: each task runs one (condition, realisation) pair
### NOTE: Runs on different cluster from scenarios 1-2

# ---- CONFIG (different cluster) ----
# Note: Legacy simulation data and methods are available in older git commits under simulations/comparison/
algo_path = "~/factor-analysis/mcem_fa_algorithm.R"
save_dir  = "/data/users/uu85g9/factor-analysis/simulations/mixed/"  # existing paper results stored in simulations_new/mixed/
missingness_data = "/data/users/uu85g9/factor-analysis/simulations/mixed/missingness_analysis.rds"
# ----------------

# Load missingness models and ECDFs
list2env(readRDS(missingness_data), .GlobalEnv)

# Load data generation functions
source("~/factor-analysis/simulations/mixed/missingness-predictor/simulate_missing_data_functions.R")

slurm_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

method_names = c(
  'full data',
  'missing data',
  'complete cases',
  'full data long v0',
  'missing data long v0',
  'complete cases long v0'
)
n_methods = length(method_names)

slurm_id_modn = slurm_id %% n_methods + 1
r = floor(slurm_id / n_methods) + 1

# ---- True loading matrix (8 vars, 3 factors, partially overlapping) ----
lambda = matrix(0, 8, 3)
lambda[1:4, 1] = 1
lambda[4:6, 2] = 1
lambda[7:8, 3] = 1

# ---- Generate data ----
N = 200
set.seed(1234 + r)
Y_complete = simulate_complete_data(ecdf_list, lambda, N, ordinal_cols, discrete_cols)
missing_mask = simulate_missingness(Y_complete, model_list)
Y = Y_complete = data.matrix(Y_complete)

# ---- Condition selection ----
source(algo_path)

method_name = method_names[slurm_id_modn]

data_condition = ((slurm_id_modn - 1) %% 3) + 1  # 1=full, 2=missing, 3=complete cases

if (data_condition == 2) {
  Y[missing_mask] = NA
} else if (data_condition == 3) {
  Y[missing_mask] = NA
  Y = Y[complete.cases(Y), ]
}

if (slurm_id_modn <= 3) {
  v0s = 10^-seq(1, 4, length.out = 4)
} else {
  v0s = c(0.02, 0.018, 0.015, 0.01, 0.005, 0.001, 0.0005, 0.0001)
}

# ---- Run ----
K = 10
print(paste("Method:", method_name, "| Realisation:", r, "| N:", nrow(Y)))
set.seed(12345)
tmp = Sys.time()

results = dpe_func(v0s, Y = Y, K = K, do_stick_breaking = TRUE,
                   px_sigma = FALSE, px_rotate = FALSE, varimax_every = -1,
                   conv_param = TRUE, nburn = 100, n_mcmc = 100, tol=0.016, verbose = -1)

time_taken = difftime(Sys.time(), tmp, units = 'mins')
print(paste('Time taken:', time_taken, 'mins'))

# ---- Evaluate ----
source("~/factor-analysis/simulations/shared/evaluate.R")
results_df = evaluate_run(results, lambda, method_name, r, time_taken)
saveRDS(results_df, paste0(save_dir, slurm_id, '.rds'))
