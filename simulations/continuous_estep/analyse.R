### Summarise the E-step robustness check, one row per number of E-step samples M.
###
### Usage: Rscript analyse.R <save_dir>

library(dplyr)
library(data.table)
options(width = 1000, tibble.width = Inf)

args = commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: Rscript analyse.R <save_dir>")
save_dir = args[1]

rds_files = list.files(path = save_dir, pattern = "^[0-9]+\\.rds$", full.names = TRUE)
cat("Found", length(rds_files), "result files in", save_dir, "\n")
if (!length(rds_files)) stop("nothing to analyse")

big_df = rbindlist(lapply(rds_files, readRDS), use.names = TRUE, fill = TRUE)

n_failed = sum(big_df$failed %in% TRUE)
if (n_failed) {
  cat("\n", n_failed, "runs FAILED:\n")
  print(big_df[big_df$failed %in% TRUE, .N, by = .(method, error_message)])
}
ok = big_df[!(big_df$failed %in% TRUE), ]

# ---- The paper's columns, one row per arm ----
summary_df = ok %>%
  group_by(n_mcmc, method) %>%
  summarise(
    struct_rec = sprintf("%.3f (%.4f)", mean(structure_recovered), sd(structure_recovered) / sqrt(n())),
    K_correct  = sprintf("%.3f (%.4f)", mean(K_correct), sd(K_correct) / sqrt(n())),
    false_pos  = sprintf("%.3f (%.4f)", mean(false_positives), sd(false_positives) / sqrt(n())),
    false_neg  = sprintf("%.3f (%.4f)", mean(false_negatives), sd(false_negatives) / sqrt(n())),
    lambda_mse = sprintf("%.6f (%.7f)", mean(lambda_cor_mse), sd(lambda_cor_mse) / sqrt(n())),
    mean_time  = sprintf("%.2f mins", mean(as.numeric(time_taken))),
    n = n(),
    .groups = "drop"
  ) %>%
  arrange(n_mcmc)
cat("\n")
print(summary_df)

saveRDS(big_df, paste0(save_dir, "all_results.rds"))
cat("\nWrote", paste0(save_dir, "all_results.rds"), "\n")
