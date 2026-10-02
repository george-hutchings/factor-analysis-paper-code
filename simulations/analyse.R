library(dplyr)
library(data.table)
options(width = 1000, tibble.width = Inf)

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript analyse.R <save_dir>")
save_dir = args[1]

# Only read numeric-named RDS files (SLURM job outputs), skip others like missingness_analysis.rds
rds_files = list.files(path = save_dir, pattern = "^[0-9]+\\.rds$", full.names = TRUE)
cat("Found", length(rds_files), "result files in", save_dir, "\n")

df_list = lapply(rds_files, readRDS)
big_df = rbindlist(df_list, use.names = TRUE, fill = TRUE)

summary_df = big_df %>%
  group_by(method) %>%
  summarise(
    struct_rec = sprintf("%.3f (%.4f)", mean(structure_recovered), sd(structure_recovered) / sqrt(n())),
    K_correct  = sprintf("%.3f (%.4f)", mean(K_correct), sd(K_correct) / sqrt(n())),
    false_pos  = sprintf("%.3f (%.4f)", mean(false_positives), sd(false_positives) / sqrt(n())),
    false_neg  = sprintf("%.3f (%.4f)", mean(false_negatives), sd(false_negatives) / sqrt(n())),
    lambda_mse = sprintf("%.6f (%.7f)", mean(lambda_cor_mse), sd(lambda_cor_mse) / sqrt(n())),
    mean_time  = sprintf("%.2f mins", mean(as.numeric(time_taken))),
    n = n()
  )
print(summary_df)

# Per-element MSE of the estimated loadings, averaged over realisations, one
# matrix per method. Plotted by plot_elementwise_mse.R.
methods = unique(big_df$method)
elementwise = lapply(methods, function(m) {
  rows = big_df[big_df$method == m, ]
  list(
    mse  = Reduce(`+`, Map(function(a, b) (a - b)^2, rows$lambda_hat, rows$lambda_true)) / nrow(rows),
    true = rows$lambda_true[[1]],
    n    = nrow(rows)
  )
})
names(elementwise) = methods
saveRDS(elementwise, paste0(save_dir, "elementwise_mse.rds"))

# Every realisation's loading matrices in one file, so further summaries never
# need the fits rerun. Copy this off the cluster alongside the summary.
saveRDS(big_df, paste0(save_dir, "all_results.rds"))
cat("Wrote", paste0(save_dir, "elementwise_mse.rds"), "and",
    paste0(save_dir, "all_results.rds"), "\n")
