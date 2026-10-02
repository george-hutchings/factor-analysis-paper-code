### E-step robustness check: Scenario 1 with the Monte Carlo E-step sample size varied,
### and nothing else varied.
###
###   M (n_mcmc)      50, 100, 200, 300, 400
###   everything else exactly continuous/run.R, the paper's Scenario 1
###   realisations    500
###
### M = 100 is the paper's setting. The data seed (12345 + r), the fit seed (12345), the
### ladder, the tolerance and every other argument are unchanged.
###
### nburn is held at 100 in every arm rather than tied to M. Burn-in equilibrates the Gibbs
### chain at the start of a DPE stage, M sets the Monte Carlo noise in the M-step averages,
### and tying them would confound the two. It is also close to moot: dpe_func resets nburn
### to 100 after the first stage whatever is passed, so only stage 1 ever sees NBURN.
###
### Usage: sbatch submit-run.sh

# ---- CONFIG (BMRC cluster) ----
algo_path = "~/factor-analysis/mcem_fa_algorithm.R"
base_dir = Sys.getenv(
  "SIM_OUT",
  "/well/nichols-nvs/users/peo100/factor-analysis/paper-code/simulations/continuous_estep/")

# Monte Carlo E-step samples, the only factor. The array index maps to an arm by
# `id %% length(MCS)` and to a realisation by `id %/% length(MCS) + 1`.
MCS    = c(50, 100, 200, 300, 400)
NBURN  = 100                          # held fixed, see header
n_reps = 500

# Scenario 1's settings, copied from continuous/run.R and not to be varied here
V0S = 10^-seq(1, 4, length.out = 4)
TOL = 0.016
N   = 200
K   = 10
# ----------------

n_arms = length(MCS)

slurm_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(slurm_id)) stop("SLURM_ARRAY_TASK_ID is not set")
if (slurm_id >= n_arms * n_reps) {
  stop("SLURM_ARRAY_TASK_ID = ", slurm_id,
       " is past the end of the study (", n_arms * n_reps, " tasks)")
}

mc_n = MCS[slurm_id %% n_arms + 1]
r    = floor(slurm_id / n_arms) + 1

dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

# ---- True loading matrix (identical to continuous/run.R, the paper's Scenario 1) ----
lambda = matrix(0, 12, 3)
lambda[1:6, 1] = 1
lambda[7:10, 2] = 1
lambda[11:12, 3] = 1
true_K = ncol(lambda)
D = nrow(lambda)

# ---- Generate data ----
# Same seed as continuous/run.R, so realisation r is the identical dataset and the M = 100
# arm is comparable to the paper's run fit by fit.
set.seed(12345 + r)
eta = matrix(rnorm(N * true_K), N, true_K)
Z = tcrossprod(eta, lambda) + matrix(rnorm(N * D), N, D)
Y = Z  # continuous: observed = latent

# ---- Run ----
source(algo_path)

print(paste("M:", mc_n, "| Realisation:", r, "| nburn:", NBURN))
set.seed(12345)
tmp = Sys.time()

# Argument for argument continuous/run.R, with n_mcmc the only difference. maxiter and
# verbose are deliberately not passed, so the defaults that produced the paper's numbers apply.
results = tryCatch(
  dpe_func(V0S, Y = Y, K = K, do_stick_breaking = TRUE,
           px_sigma = FALSE, px_rotate = FALSE, varimax_every = -1,
           conv_param = TRUE, nburn = NBURN, n_mcmc = mc_n, tol = TOL),
  error = function(e) e)

time_taken = difftime(Sys.time(), tmp, units = 'mins')
failed = inherits(results, "error")
print(paste('Time taken:', time_taken, 'mins', if (failed) '(FAILED)' else ''))

# ---- Evaluate ----
source("~/factor-analysis/simulations/shared/evaluate.R")

method_name = paste0("M = ", mc_n)

if (failed) {
  print(paste("FAILED:", conditionMessage(results)))
  results_df = data.frame(method = method_name, realisation = r,
                          time_taken = time_taken, iterations = NA_integer_)
  results_df$failed = TRUE
  results_df$error_message = conditionMessage(results)
} else {
  # The paper's own columns, in the paper's own column ordering convention
  results_df = evaluate_run(results, lambda, method_name, r, time_taken)
  results_df$failed = FALSE
  results_df$error_message = NA_character_
}

results_df$n_mcmc   = mc_n
results_df$nburn    = NBURN
results_df$schedule = paste(V0S, collapse = "-")
results_df$tol      = TOL
results_df$N        = N

saveRDS(results_df, paste0(base_dir, slurm_id, '.rds'))
