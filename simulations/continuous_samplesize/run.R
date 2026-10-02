### Sample-size study (supplementary Table S2 and Fig S4): scenario 1 at six sample sizes.
###
### Only N is varied. Everything else is the configuration adopted for this study:
###
###   ladder      0.05, 0.03, 0.015, 0.005   (4 rungs)
###   budget      400 total, 100 per rung
###   tol         1e-6
###   n_mcmc      100
###   likelihood  rank, with px_sigma on
###   N           50, 250, 500, 1000, 2500, 5000
###   reps        500
###
### One task is one (N, realisation) pair: N = Ns[id %% 6 + 1], realisation = id %/% 6 + 1.
###
### Usage: sbatch submit-run.sh.

# ---- CONFIG (BMRC cluster) ----
algo_rank = "~/factor-analysis/mcem_fa_algorithm.R"
base_dir = Sys.getenv(
  "SIM_OUT",
  "/well/nichols-nvs/users/peo100/factor-analysis/paper-code/simulations/continuous_samplesize/")

Ns     = c(50, 250, 500, 1000, 2500, 5000)
n_reps = 500

SCHEDULE = c(0.05, 0.03, 0.015, 0.005)
BUDGET   = 400    # total iterations, split evenly over the rungs
N_MCMC   = 100

arm = list(name = 'rank', algo = algo_rank)
# ----------------

n_cells = length(Ns)   # 6

slurm_id = as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(slurm_id)) stop("SLURM_ARRAY_TASK_ID is not set")
if (slurm_id >= n_cells * n_reps) {
  stop("SLURM_ARRAY_TASK_ID = ", slurm_id,
       " is past the end of the study (", n_cells * n_reps, " tasks)")
}

cell = slurm_id %% n_cells
r    = floor(slurm_id / n_cells) + 1
N    = Ns[cell + 1]

v0s     = SCHEDULE
maxiter = BUDGET / length(v0s)
mc_n    = N_MCMC

save_dir = paste0(base_dir, "n", N, "/")
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

# ---- True loading matrix (identical to continuous/run.R, the paper's Scenario 1) ----
lambda = matrix(0, 12, 3)
lambda[1:6, 1] = 1
lambda[7:10, 2] = 1
lambda[11:12, 3] = 1
true_K = ncol(lambda)
D = nrow(lambda)

# ---- Generate data ----
# The data seed depends only on the realisation, as in the paper's continuous scenario, and the
# fit seed below is a constant, so each (N, realisation) pair is reproducible on its own.
set.seed(12345 + r)
eta = matrix(rnorm(N * true_K), N, true_K)
Z = tcrossprod(eta, lambda) + matrix(rnorm(N * D), N, D)
Y = Z  # continuous: observed = latent

# ---- Run ----
K = 10
source(arm$algo)

print(paste("Arm:", arm$name, "| N:", N, "| Realisation:", r,
            "| rungs:", length(v0s), "| maxiter/rung:", maxiter,
            "| total:", maxiter * length(v0s), "| n_mcmc:", mc_n))
set.seed(12345)
tmp = Sys.time()

results = tryCatch(
  dpe_func(v0s, Y = Y, K = K, do_stick_breaking = TRUE,
           px_sigma = TRUE, px_rotate = FALSE, varimax_every = -1,
           conv_param = TRUE, nburn = mc_n, n_mcmc = mc_n,
           tol = 1e-6, maxiter = maxiter, verbose = -1),
  error = function(e) e)

time_taken = difftime(Sys.time(), tmp, units = 'mins')
failed = inherits(results, "error")
print(paste('Time taken:', time_taken, 'mins', if (failed) '(FAILED)' else ''))

# ---- Evaluate ----
# Only the final rung is kept, since the paper reports the fitted model, which is the last rung.
source("~/factor-analysis/simulations/shared/evaluate.R")

if (failed) {
  print(paste("FAILED:", conditionMessage(results)))
  results_df = data.frame(method = arm$name, realisation = r, time_taken = time_taken,
                          iterations = NA_integer_)
  results_df$failed = TRUE
  results_df$error_message = conditionMessage(results)
} else {
  # evaluate_run() gives the columns the tables report - structure_recovered, K_correct,
  # lambda_cor_mse, false_positives, false_negatives - computed in colSums(w) order.
  results_df = evaluate_run(results, lambda, arm$name, r, time_taken)
  results_df$failed = FALSE
  results_df$error_message = NA_character_
}

results_df$N        = N
results_df$arm      = arm$name
results_df$schedule = paste(SCHEDULE, collapse = "-")
results_df$maxiter  = maxiter
results_df$n_mcmc   = mc_n

saveRDS(results_df, paste0(save_dir, slurm_id, '.rds'))
