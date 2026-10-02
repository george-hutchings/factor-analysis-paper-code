# Consolidate Habib V2 PX clinical results
# Prints thresholded lambda for each of 10 runs

variables <- c("EDSS", "T25FWM", "HPT9M", "PASAT", "SDMT", "VOLT2", "NBV2", "NUMGDT1", "RELAPSE")
base_path <- '/data/users/uu85g9/factor-analysis/paper-code/noms/clinical'

count_factors <- function(lambda) sum(colSums(abs(lambda) > 0) > 0)

# Threshold lambda using w
threshold_lambda <- function(res, variables) {
  n_dpe <- length(res)
  lam <- res[[n_dpe]]$lambda
  w <- res[[n_dpe]]$w
  rownames(lam) <- variables
  lam[w < 0.5] <- 0
  # Drop all-zero columns
  active <- colSums(abs(lam) > 0) > 0
  lam <- lam[, active, drop = FALSE]
  # Make columns positive (flip sign so colSums > 0; leave unchanged if colSum==0)
  cs <- colSums(lam)
  cs[cs == 0] <- 1
  lam <- t(t(lam) * sign(cs))
  # Normalise loadings to correlations with the latent Gaussian Z_d
  lam <- lam / sqrt(1 + rowSums(lam^2))
  # Reorder rows into the interpretive grouping used in Table 4
  lam <- lam[c("EDSS", "T25FWM", "HPT9M", "SDMT", "PASAT",
               "NUMGDT1", "VOLT2", "NBV2", "RELAPSE"), , drop = FALSE]
  # Order columns by the row of each column's dominant loading, which
  # yields the staircase and places brain damage before relapse
  lam[, order(apply(abs(lam), 2, which.max)), drop = FALSE]
}

# Summary table header
cat(sprintf("%-6s %8s %10s %12s\n", "Run", "Factors", "Tot Iter", "Q Function"))
cat(strrep('-', 40), '\n')

all_factors <- c()

for (r in 1:10) {
  file_path <- paste0(base_path, '/results_', r, '.rds')
  if (!file.exists(file_path)) {
    cat(sprintf("%-6d %8s\n", r, "MISSING"))
    next
  }
  res <- readRDS(file_path)
  lam_t <- threshold_lambda(res, variables)
  nf <- count_factors(lam_t)
  all_factors <- c(all_factors, nf)
  total_iter <- sum(sapply(res, function(x) x$iterations))
  logprior <- tail(res[[length(res)]]$q_function, 1)
  cat(sprintf("%-6d %8d %10d %12.1f\n", r, nf, total_iter, logprior))
}

cat(strrep('-', 40), '\n')
cat(sprintf("Mean factors: %.2f\n\n", mean(all_factors)))

# Print thresholded lambda for each run
for (r in 1:10) {
  file_path <- paste0(base_path, '/results_', r, '.rds')
  if (!file.exists(file_path)) next
  res <- readRDS(file_path)
  lam_t <- threshold_lambda(res, variables)
  cat(sprintf("\n===== Run %d (%d factors) =====\n", r, count_factors(lam_t)))
  print(round(lam_t, 3))
}
