### Paper table for the sample-size study (supplementary Table S2).
###
### Usage: Rscript analyse.R <base_save_dir>
###
### Produces one table, rows are N, columns are those of simulations/analyse.R:
###
###   struct_rec  structure_recovered   the whole support recovered exactly
###   K_correct   the right number of factors
###   false_pos   false_positives       element counts, colSums(w) column order
###   false_neg   false_negatives
###   lambda_mse  lambda_cor_mse        MSE of the correlation-scale loadings, same order
###   mean_time   n
###   iters       the mean total number of iterations over the four rungs (400 if none stopped early)
###
### Standard errors over realisations in brackets, as there.

library(dplyr)
library(data.table)
options(width = 1000, tibble.width = Inf)

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript analyse.R <base_save_dir>")
base_dir = args[1]

n_dirs = list.files(base_dir, pattern = "^n[0-9]+$", full.names = TRUE)
if (length(n_dirs) == 0) stop("no n<N>/ subdirectories found in ", base_dir)
n_dirs = n_dirs[order(as.numeric(sub("^.*/n", "", n_dirs)))]

df_list = list()
for (d in n_dirs) {
  rds = list.files(d, pattern = "^[0-9]+\\.rds$", full.names = TRUE)
  cat("Found", length(rds), "result files in", d, "\n")
  if (length(rds) == 0) next
  df_list[[d]] = rbindlist(lapply(rds, readRDS), use.names = TRUE, fill = TRUE)
}
all_df = rbindlist(df_list, use.names = TRUE, fill = TRUE)
if (is.null(all_df$failed)) all_df$failed = FALSE
all_df$failed[is.na(all_df$failed)] = FALSE

ARMS = "rank"
LABELS = c(rank = "Rank likelihood")
cat("\n", nrow(all_df), "runs total\n\n")

# ---- Failures ----
# A fit that fails is caught by run.R. The count is reported rather than hidden, and the tables
# carry an n column so no row can be read without its sample size.
fail_df = all_df %>%
  group_by(arm, N) %>%
  summarise(n_runs = n(), n_failed = sum(failed),
            fail_rate = sprintf("%.3f", sum(failed) / n()), .groups = "drop") %>%
  arrange(arm, N)
cat("Failures:\n")
print(as.data.frame(fail_df))
if (sum(all_df$failed)) {
  cat("\nDistinct errors:\n"); print(table(all_df$error_message[all_df$failed]))
}

ok = all_df[!all_df$failed, ]
if (nrow(ok) == 0) stop("every run failed")

se  = function(x) sd(x) / sqrt(length(x))
fmt = function(x, d = 3) sprintf(paste0("%.", d, "f (%.", d + 1, "f)"), mean(x), se(x))

# ---- The paper table ----
paper_table = function(s) {
  s %>%
    group_by(N) %>%
    summarise(
      struct_rec = fmt(structure_recovered),
      K_correct  = fmt(K_correct),
      false_pos  = fmt(false_positives),
      false_neg  = fmt(false_negatives),
      lambda_mse = sprintf("%.6f (%.7f)", mean(lambda_cor_mse), se(lambda_cor_mse)),
      mean_time  = sprintf("%.2f mins", mean(as.numeric(time_taken))),
      n          = n(),
      iters      = sprintf("%.0f", mean(iterations)),
      .groups = "drop"
    ) %>%
    arrange(N)
}

tables = list()
for (a in ARMS) {
  s = ok[ok$arm == a, ]
  if (!nrow(s)) { cat("\n\nNo surviving runs for arm", a, "\n"); next }
  tables[[a]] = paper_table(s)
  cat("\n\n=== ", LABELS[[a]], " ===\n\n", sep = "")
  print(as.data.frame(tables[[a]]), row.names = FALSE)
}
cat("\nTruth: K = 3, 12 variables in factors of size 6 / 4 / 2, raw loadings 1.0.\n")
cat("Ladder 0.05, 0.03, 0.015, 0.005 at 100 iterations per rung. Standard errors in brackets.\n")
cat("Columns and their definitions match simulations/analyse.R.\n")

# ---- Per-element MSE, in the same shape simulations/analyse.R writes ----
# So plot_elementwise.R can draw the grid later without any of this being rerun. Built from
# lambda_hat, i.e. the colSums(w) column order, to stay consistent with the table above.
elementwise = list()
for (a in ARMS) for (Nv in sort(unique(ok$N))) {
  rows = ok[ok$arm == a & ok$N == Nv, ]
  if (!nrow(rows)) next
  elementwise[[paste0(a, " N=", Nv)]] = list(
    mse  = Reduce(`+`, Map(function(x, y) (x - y)^2, rows$lambda_hat, rows$lambda_true)) / nrow(rows),
    true = rows$lambda_true[[1]],
    n    = nrow(rows))
}
saveRDS(elementwise, file.path(base_dir, "elementwise_mse.rds"))

saveRDS(ok,     file.path(base_dir, "all_results.rds"))
saveRDS(all_df, file.path(base_dir, "all_runs.rds"))
for (a in names(tables)) {
  write.csv(tables[[a]], file.path(base_dir, paste0("table_", a, ".csv")), row.names = FALSE)
}
write.csv(fail_df, file.path(base_dir, "failures.csv"), row.names = FALSE)

cat("\n\nWrote, in", base_dir, ":\n")
cat("  table_rank.csv         paper table, rank likelihood\n")
cat("  all_results.rds        surviving runs, one row each\n")
cat("  all_runs.rds           every run attempted, failures included\n")
cat("  elementwise_mse.rds    per-element MSE per N, for plot_elementwise.R\n")
cat("  failures.csv           failure count and rate per cell\n")
