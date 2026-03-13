### Interactive continuous sim - set slurm_id and source()
slurm_id = 1

algo_path = "~/factor-analysis/mcem_fa_algorithm.R"
method_names = c('our method', 'our method long v0', 'rockova 2016', 'factanal K known', 'pca K known')
n_methods = length(method_names)
slurm_id_modn = slurm_id %% n_methods + 1
r = floor(slurm_id / n_methods) + 1

lambda = matrix(0, 12, 3)
lambda[1:6, 1] = 1
lambda[7:10, 2] = 1
lambda[11:12, 3] = 1
D = nrow(lambda)

set.seed(12345 + r)
N = 200
eta = matrix(rnorm(N * ncol(lambda)), N, ncol(lambda))
Z = tcrossprod(eta, lambda) + matrix(rnorm(N * D), N, D)
Y = Z

K = 10
v0s = NULL
source("~/factor-analysis/simulations/methods/pca_factanal_nmf.R")

if (slurm_id_modn == 1) {
  source(algo_path)
  v0s = 10^-seq(1, 4, length.out = 4)
} else if (slurm_id_modn == 2) {
  source(algo_path)
  v0s = c(0.02, 0.018, 0.015, 0.01, 0.005, 0.001, 0.0005, 0.0001)
} else if (slurm_id_modn == 3) {
  source("~/factor-analysis/simulations/methods/veronika2016_wrapper.R")
  dpe_func = function(v0s, ...) list(veronika2016(...))
} else if (slurm_id_modn == 4) {
  K = ncol(lambda)
  dpe_func = function(v0s, ...) list(list(w = NA, lambda = NA, lambda_corr = my_factanal(...)))
} else if (slurm_id_modn == 5) {
  K = ncol(lambda)
  dpe_func = function(v0s, ...) list(list(w = NA, lambda = NA, lambda_corr = pca_varimax(...)))
}

set.seed(12345)
if (slurm_id_modn <= 2) {
  results = dpe_func(v0s, Y = Y, K = K, do_stick_breaking = TRUE, px_sigma = FALSE,
                     px_rotate = FALSE, varimax_every = -1, conv_param = TRUE,
                     nburn = 100, n_mcmc = 100, tol=0.016, verbose=0)
} else {
  results = dpe_func(v0s, Y = Y, K = K, do_stick_breaking = TRUE)
}

source("~/factor-analysis/simulations/shared/evaluate.R")
print_evaluation(results, lambda, method_names[slurm_id_modn], r)
