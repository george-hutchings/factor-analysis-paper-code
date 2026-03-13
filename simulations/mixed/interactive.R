### Interactive mixed/missing sim - set slurm_id and source()
### Requires missingness_analysis.rds
slurm_id = 0

algo_path = "~/factor-analysis/mcem_fa_algorithm.R"
missingness_data = "/data/users/uu85g9/factor-analysis/simulations/mixed/missingness_analysis.rds"

list2env(readRDS(missingness_data), .GlobalEnv)
source("~/factor-analysis/simulations/mixed/missingness-predictor/simulate_missing_data_functions.R")

method_names = c('full data', 'missing data', 'complete cases',
                 'full data long v0', 'missing data long v0', 'complete cases long v0')
n_methods = length(method_names)
slurm_id_modn = slurm_id %% n_methods + 1
r = floor(slurm_id / n_methods) + 1

lambda = matrix(0, 8, 3)
lambda[1:4, 1] = 1
lambda[4:6, 2] = 1
lambda[7:8, 3] = 1

N = 200
set.seed(1234 + r)
Y_complete = simulate_complete_data(ecdf_list, lambda, N, ordinal_cols, discrete_cols)
missing_mask = simulate_missingness(Y_complete, model_list)
Y = data.matrix(Y_complete)

source(algo_path)
data_condition = ((slurm_id_modn - 1) %% 3) + 1
if (data_condition == 2) {
  Y[missing_mask] = NA
}
if (data_condition == 3) {
  Y[missing_mask] = NA
  Y = Y[complete.cases(Y), ]
}

if (slurm_id_modn <= 3) {
  v0s = 10^-seq(1, 4, length.out = 4)
} else {
  v0s = c(0.02, 0.018, 0.015, 0.01, 0.005, 0.001, 0.0005, 0.0001)
}

K = 10
set.seed(12345)
results = dpe_func(v0s, Y = Y, K = K, do_stick_breaking = TRUE, px_sigma = FALSE,
                   px_rotate = FALSE, varimax_every = -1, conv_param = TRUE,
                   nburn = 100, n_mcmc = 100, tol=0.016, verbose=0)

source("~/factor-analysis/simulations/shared/evaluate.R")
print_evaluation(results, lambda, method_names[slurm_id_modn], r, N = nrow(Y))
