# EDSS and Functional Subscore Correlation Analysis
# Partial correlations controlling for baseline covariates

library(data.table)
library(dplyr)
library(tidyr)
library(ppcor)  # For partial correlation

set.seed(1234)

# === Configuration ===
base_path <- "/data/users/uu85g9/factor-analysis/noms/mri/clinical_analysis"
compressed_dir <- file.path(base_path, "compressed_data")
output_dir <- file.path(base_path, "results")
dir.create(output_dir, showWarnings = FALSE)

# Clinical data paths
CLINICAL_PATH <- "/data/ms/processed/clinical/projects/longitudinal"
SUBJECTS_PATH <- "/data/users/uu85g9/factor-analysis/noms/mri/inputs/subjects.rds"

# Variables
edss_vars <- c("EDSS", "VISFNC", "BRNFNC", "BOWFNC", "CLRFNC", "CRBFNC", "PYRFNC", "SENFNC")
baseline_vars <- c("AGE", "ACTVTRT", "SEX", "DURFS", "RELPST1Y", "NUMGDT1", "NBV2")
analysis_months <- c(-1, 6, 12, 24)  # Baseline + follow-up timepoints

# === Load Data ===

# Subjects with MRI
subjects_df <- readRDS(SUBJECTS_PATH) %>% dplyr::select(USUBJID)

# Longitudinal clinical data
clinical_data <- fread(file.path(CLINICAL_PATH, "longitudinal-table-CORE.csv")) %>%
  filter(STUDY %in% c("CFTY720D2301", "CFTY720D2309")) %>%
  mutate(SEX = ifelse(SEX == "Male", 1, 0)) %>%
  filter(USUBJID %in% subjects_df$USUBJID)

# Print sample sizes per month
for (m in analysis_months) {
  month_data <- clinical_data %>%
    filter(MONTH == m) %>%
    dplyr::select(USUBJID, all_of(edss_vars), all_of(baseline_vars)) %>%
    na.omit()
  cat(sprintf("Month %d: %d subjects with complete data\n", m, nrow(month_data)))
}

# === Correlation Function ===
perform_partial_correlations <- function(df, target_vars, factor_vars, control_vars, method = "spearman") {
  results_list <- list()

  for (target in target_vars) {
    for (factor_name in factor_vars) {
      test_data <- df[, c(target, factor_name, control_vars), drop = FALSE] %>% na.omit()

      test <- tryCatch(
        pcor.test(test_data[[target]], test_data[[factor_name]],
                  test_data[, control_vars, drop = FALSE], method = method),
        error = function(e) {
          cat(sprintf("  Warning: Failed for %s vs %s (n=%d): %s\n", 
                      target, factor_name, nrow(test_data), e$message))
          NULL
        }
      )

      if (!is.null(test)) {
        results_list[[length(results_list) + 1]] <- data.frame(
          clinical_var = target,
          factor = factor_name,
          correlation = test$estimate,
          p_value = test$p.value,
          n = nrow(test_data)
        )
      }
    }
  }

  results <- bind_rows(results_list)
  if (nrow(results) > 0) {
    results$p_value_bh <- p.adjust(results$p_value, method = "BH")
  }
  return(results)
}

# === Process Each Compressed Dataset x Month ===

# Find all compressed Eeta files
eeta_files <- list.files(compressed_dir, pattern = "^Eeta_.*\\.rds$", full.names = TRUE)

all_results <- list()

for (eeta_file in eeta_files) {
  dataset_name <- gsub("Eeta_|\\.rds", "", basename(eeta_file))
  cat(sprintf("\n=== Processing: %s ===\n", dataset_name))

  # Load Eeta and metadata
  Eeta <- readRDS(eeta_file)
  meta_file <- file.path(compressed_dir, sprintf("meta_%s.rds", dataset_name))
  meta <- if (file.exists(meta_file)) readRDS(meta_file) else NULL

  # Create factor names
  colnames(Eeta) <- paste0("V", seq_len(ncol(Eeta)))

  # Combine with subject IDs
  eeta_df <- cbind(subjects_df, as.data.frame(Eeta))

  # Identify non-redundant factors (if metadata available)
  if (!is.null(meta)) {
    active_cols <- which(colSums(meta$w > 0.5) > 0)
    factor_vars <- paste0("V", active_cols)
    cat(sprintf("Active factors: %d\n", length(factor_vars)))
  } else {
    factor_vars <- colnames(Eeta)
  }

  # Loop over timepoints
  for (month in analysis_months) {
    cat(sprintf("\n--- Month %d ---\n", month))

    # Get clinical data for this month
    month_clinical <- clinical_data %>%
      filter(MONTH == month) %>%
      dplyr::select(USUBJID, all_of(edss_vars), all_of(baseline_vars)) %>%
      na.omit()

    # Merge with Eeta
    merged_df <- inner_join(eeta_df, month_clinical, by = "USUBJID") %>%
      as.data.frame()

    cat(sprintf("Merged subjects: %d\n", nrow(merged_df)))

    # Run correlations (Spearman for ordinal EDSS data)
    results <- perform_partial_correlations(
      merged_df, edss_vars, factor_vars, baseline_vars, method = "spearman"
    )

    results$dataset <- dataset_name
    results$month <- month
    all_results[[paste(dataset_name, month, sep = "_m")]] <- results

    # Print significant results
    sig_results <- results %>% filter(p_value_bh < 0.05)
    if (nrow(sig_results) > 0) {
      cat(sprintf("Significant correlations (BH p < 0.05): %d\n", nrow(sig_results)))
    } else {
      cat("No significant correlations found.\n")
    }
  }
}

# === Save Combined Results ===
combined_results <- bind_rows(all_results)
saveRDS(combined_results, file.path(output_dir, "edss_correlation_results.rds"))

# Summary table by dataset and month
summary_table <- combined_results %>%
  filter(p_value_bh < 0.05) %>%
  group_by(dataset, month, clinical_var) %>%
  summarise(n_sig_factors = n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = clinical_var, values_from = n_sig_factors, values_fill = 0)

cat("\n=== Summary: Number of Significant Factors per Clinical Variable ===\n")
print(summary_table)

# Top significant correlations across all
cat("\n=== Top 20 Significant Correlations (all datasets/months) ===\n")
top_results <- combined_results %>%
  filter(p_value_bh < 0.05) %>%
  arrange(p_value_bh) %>%
  head(20)
print(top_results)

write.csv(summary_table, file.path(output_dir, "edss_correlation_summary.csv"), row.names = FALSE)
write.csv(combined_results, file.path(output_dir, "edss_correlation_all.csv"), row.names = FALSE)
cat(sprintf("\nResults saved to: %s\n", output_dir))
