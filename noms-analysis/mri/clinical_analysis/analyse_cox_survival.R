# Cox Proportional Hazards Analysis
# Predicting 3-month confirmed disability worsening (M3CDW)

library(data.table)
library(dplyr)
library(tidyr)
library(survival)
library(RNifti)

set.seed(1234)

# === Configuration ===
base_path <- "/data/users/uu85g9/factor-analysis/noms/mri/clinical_analysis"
compressed_dir <- file.path(base_path, "compressed_data")
output_dir <- file.path(base_path, "results")
dir.create(output_dir, showWarnings = FALSE)

# Data paths
LONGITUDINAL_PATH <- "/data/ms/processed/mri/MS_Share/4George/longitudinal.csv"
SUBJECTS_PATH <- "/data/users/uu85g9/factor-analysis/noms/mri/inputs/subjects.rds"
MASK_FILE <- "/data/ms/processed/mri/MS_Share/4George/final_mask.nii.gz"

# Baseline covariates to include in Cox model
baseline_vars <- c("AGE", "ACTVTRT", "SEX", "DURFS", "RELPST1Y", "EDSS", "NUMGDT1", "NBV2")

# === Load Survival Data ===

# Subjects with MRI
subjects_df <- readRDS(SUBJECTS_PATH) %>% dplyr::select(USUBJID)

# Longitudinal data
df_all <- read.csv(LONGITUDINAL_PATH)

# Events: M3CDW = 1
m3cdw <- df_all %>%
  filter(M3CDW == 1) %>%
  mutate(status = 1)
m3cdw <- m3cdw[!duplicated(m3cdw$USUBJID), ]

# Censored: last observation for non-events
cns <- df_all %>%
  filter(!USUBJID %in% m3cdw$USUBJID) %>%
  group_by(USUBJID) %>%
  filter(DAY == max(DAY)) %>%
  mutate(status = 0)

# Combine
surv_df <- bind_rows(m3cdw, cns) %>%
  dplyr::select(USUBJID, STUDY, MONTH, DAY, YEARS, M3CDW, status)

# Baseline covariates
base <- df_all %>%
  filter(DAY == 1) %>%
  dplyr::select(USUBJID, STUDY, ACTVTRT, AGE, SEX, DURFS, RELPST1Y, EDSS, T25FWM, HPT9M, PASAT, NUMGDT1, VOLT2, NBV2)

surv_df <- left_join(surv_df, base, by = c("USUBJID", "STUDY"))
surv_df <- surv_df[complete.cases(surv_df), ]

# Filter survival data to subjects with MRI factor scores
surv_df <- surv_df %>% filter(USUBJID %in% subjects_df$USUBJID)

cat(sprintf("Subjects with complete survival data: %d\n", nrow(surv_df)))
cat(sprintf("Events (M3CDW): %d (%.1f%%)\n", sum(surv_df$status), 100 * mean(surv_df$status)))

# === VOLT2 Comparison Model (baseline) ===
# Run Cox model with just VOLT2 + baseline covariates as comparison
# Uses same subject set as factor models for valid comparison
cat("\n=== VOLT2 Comparison Model ===\n")
volt2_formula <- as.formula(paste("Surv(DAY, status) ~ VOLT2 +", paste(baseline_vars, collapse = " + ")))
fit_volt2 <- coxph(volt2_formula, data = surv_df)
cat("VOLT2 model summary:\n")
print(summary(fit_volt2))

# Extract VOLT2 results for comparison
volt2_summary <- summary(fit_volt2)
volt2_results <- data.frame(
  variable = rownames(volt2_summary$coefficients),
  coef = volt2_summary$coefficients[, "coef"],
  exp_coef = volt2_summary$coefficients[, "exp(coef)"],
  se = volt2_summary$coefficients[, "se(coef)"],
  z = volt2_summary$coefficients[, "z"],
  p_value = volt2_summary$coefficients[, "Pr(>|z|)"],
  dataset = "VOLT2_only"
)
saveRDS(list(fit = fit_volt2, results = volt2_results), file.path(output_dir, "cox_volt2_comparison.rds"))

# === Load Mask for NIfTI Output ===
mask_nifti <- readNifti(MASK_FILE)
mask <- as.logical(mask_nifti > 0)

# === Cox Analysis Function ===
run_cox_analysis <- function(eeta_df, surv_df, factor_vars, baseline_vars) {
  # Merge
  my_df <- inner_join(eeta_df, surv_df, by = "USUBJID") %>%
    as.data.frame()
  my_df <- my_df[complete.cases(my_df[, c(factor_vars, baseline_vars, "DAY", "status")]), ]

  cat(sprintf("Cox analysis on %d subjects\n", nrow(my_df)))

  # Build formula
  formula_str <- paste("Surv(DAY, status) ~", paste(c(factor_vars, baseline_vars), collapse = " + "))
  formula_obj <- as.formula(formula_str)

  # Fit model
  fit <- coxph(formula_obj, data = my_df)

  # Extract results
  summary_fit <- summary(fit)
  results <- data.frame(
    variable = rownames(summary_fit$coefficients),
    coef = summary_fit$coefficients[, "coef"],
    exp_coef = summary_fit$coefficients[, "exp(coef)"],
    se = summary_fit$coefficients[, "se(coef)"],
    z = summary_fit$coefficients[, "z"],
    p_value = summary_fit$coefficients[, "Pr(>|z|)"]
  )

  return(list(fit = fit, results = results, n = nrow(my_df)))
}

# === Process Each Compressed Dataset ===

eeta_files <- list.files(compressed_dir, pattern = "^Eeta_.*\\.rds$", full.names = TRUE)

all_results <- list()

for (eeta_file in eeta_files) {
  dataset_name <- gsub("Eeta_|\\.rds", "", basename(eeta_file))
  cat(sprintf("\n=== Processing: %s ===\n", dataset_name))

  # Load Eeta and metadata
  Eeta <- readRDS(eeta_file)
  meta_file <- file.path(compressed_dir, sprintf("meta_%s.rds", dataset_name))
  meta <- if (file.exists(meta_file)) readRDS(meta_file) else NULL

  # Load lambda from original checkpoint
  lambda <- if (!is.null(meta) && file.exists(meta$checkpoint_file)) readRDS(meta$checkpoint_file)$lambda else NULL

  # Create factor names
  colnames(Eeta) <- paste0("V", seq_len(ncol(Eeta)))

  # Combine with subject IDs
  eeta_df <- cbind(subjects_df, as.data.frame(Eeta))

  # Identify active factors
  if (!is.null(meta)) {
    active_cols <- which(colSums(meta$w > 0.5) > 0)
    factor_vars <- paste0("V", active_cols)
    cat(sprintf("Active factors: %d\n", length(factor_vars)))
  } else {
    factor_vars <- colnames(Eeta)
  }

  # Run Cox analysis
  cox_result <- run_cox_analysis(eeta_df, surv_df, factor_vars, baseline_vars)

  # Filter to significant latent factors
  sig_factors <- cox_result$results %>%
    filter(p_value < 0.05, grepl("^V", variable)) %>%
    arrange(p_value) %>%
    mutate(factor_num = as.numeric(gsub("V", "", variable)))

  cat(sprintf("Significant factors (p < 0.05): %d\n", nrow(sig_factors)))

  if (nrow(sig_factors) > 0) {
    if (!is.null(meta)) {
      n_voxels <- sum(mask)
      w <- meta$w
      sig_factor_nums <- sig_factors$factor_num

      # Compute average non-zero loading per significant factor
      if (!is.null(lambda)) {
        avg_loadings <- sapply(sig_factor_nums, function(k) {
          active <- w[1:n_voxels, k] > 0.5
          if (any(active)) mean(lambda[1:n_voxels, k][active]) else NA
        })
        sig_factors$avg_loading <- avg_loadings
      }
    }

    print(sig_factors %>% dplyr::select(variable, coef, exp_coef, p_value,
                                        any_of("avg_loading")))

    if (!is.null(meta)) {
      # Winner-takes-all atlas (3D) — kept from original
      w_sig <- w[, sig_factor_nums, drop = FALSE]
      has_active <- rowSums(w_sig[1:n_voxels, , drop = FALSE] > 0.5) > 0
      max_factor <- rep(0L, n_voxels)
      if (any(has_active)) {
        max_factor[has_active] <- sig_factor_nums[max.col(w_sig[has_active, , drop = FALSE], ties.method = "first")]
      }
      atlas_vec <- rep(0L, length(mask))
      atlas_vec[mask] <- max_factor
      factor_atlas <- array(atlas_vec, dim(mask_nifti))

      output_atlas <- file.path(output_dir, sprintf("coxph_sig_factors_%s.nii.gz", dataset_name))
      writeNifti(factor_atlas, output_atlas, template = mask_nifti)
      cat(sprintf("Atlas NIfTI saved: %s\n", output_atlas))

      # 3D NIfTI per significant factor (thresholded loading columns)
      if (!is.null(lambda)) {
        loadings_dir <- file.path(output_dir, sprintf("coxph_sig_loadings_%s", dataset_name))
        dir.create(loadings_dir, showWarnings = FALSE)

        for (i in seq_along(sig_factor_nums)) {
          k <- sig_factor_nums[i]
          active <- w[1:n_voxels, k] > 0.5
          vol_vec <- rep(0, length(mask))
          vol_vec[mask] <- ifelse(active, lambda[1:n_voxels, k], 0)
          vol_3d <- array(vol_vec, dim(mask_nifti))

          output_nifti <- file.path(loadings_dir, sprintf("moco_V%d.nii.gz", k))
          writeNifti(vol_3d, output_nifti, template = mask_nifti)
        }

        cat(sprintf("Loadings NIfTIs saved (%d files): %s/\n", length(sig_factor_nums), loadings_dir))
      }
    }
  }

  # Add avg_loading to cox_result for all factor variables
  if (!is.null(lambda) && nrow(sig_factors) > 0 && "avg_loading" %in% names(sig_factors)) {
    cox_result$results$avg_loading <- NA
    for (j in seq_len(nrow(sig_factors))) {
      row_idx <- which(cox_result$results$variable == sig_factors$variable[j])
      cox_result$results$avg_loading[row_idx] <- sig_factors$avg_loading[j]
    }
  } else {
    cox_result$results$avg_loading <- NA
  }

  cox_result$results$dataset <- dataset_name
  all_results[[dataset_name]] <- cox_result$results
}

# === Save Combined Results ===
combined_results <- bind_rows(all_results)
saveRDS(combined_results, file.path(output_dir, "cox_survival_results.rds"))

# Summary of significant factors per dataset
summary_cox <- combined_results %>%
  filter(p_value < 0.05, grepl("^V", variable)) %>%
  group_by(dataset) %>%
  summarise(
    n_sig_factors = n(),
    n_protective = sum(coef < 0),
    n_risk = sum(coef > 0),
    .groups = "drop"
  )

cat("\n=== Summary: Significant Factors in Cox Models ===\n")
print(summary_cox)

# Compare with VOLT2 model
cat("\n=== VOLT2 Comparison ===\n")
volt2_row <- volt2_results %>% filter(variable == "VOLT2")
cat(sprintf("VOLT2: HR = %.3f (95%% CI: %.3f-%.3f), p = %.4f\n",
            volt2_row$exp_coef,
            exp(volt2_row$coef - 1.96 * volt2_row$se),
            exp(volt2_row$coef + 1.96 * volt2_row$se),
            volt2_row$p_value))

write.csv(summary_cox, file.path(output_dir, "cox_survival_summary.csv"), row.names = FALSE)
write.csv(combined_results, file.path(output_dir, "cox_survival_all.csv"), row.names = FALSE)
cat(sprintf("\nResults saved to: %s\n", output_dir))
