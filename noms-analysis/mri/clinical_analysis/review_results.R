library(dplyr)

# Methods
# -------
# EDSS analysis: separate Spearman partial correlation for each
# (factor, clinical variable) pair, controlling for AGE, ACTVTRT, SEX, DURFS,
# RELPST1Y, NUMGDT1, NBV2. Clinical variables: EDSS + 7 functional subscores
# (VISFNC, BRNFNC, BOWFNC, CLRFNC, CRBFNC, PYRFNC, SENFNC).
# Tested at months -1, 6, 12, 24. BH correction applied to all tests within
# each dataset-month (e.g. 8 clinical vars x K factors corrected together).
#
# Cox analysis: Cox proportional hazards model predicting 3-month confirmed
# disability worsening (M3CDW). All active factors entered jointly alongside
# baseline covariates (AGE, ACTVTRT, SEX, DURFS, RELPST1Y, EDSS, NUMGDT1,
# NBV2). Raw p < 0.05 threshold (single model, no multiple testing correction).
# VOLT2-only model fitted on the same subjects as a baseline comparison.
#
# Brain maps: winner-takes-all NIfTI showing voxels where significant Cox
# factors have w > 0.5. View in FSLeyes.

base_path <- "/data/users/uu85g9/factor-analysis/noms/mri/clinical_analysis"
compressed_dir <- file.path(base_path, "compressed_data")
results_dir <- file.path(base_path, "results")

edss <- readRDS(file.path(results_dir, "edss_correlation_results.rds"))
cox <- readRDS(file.path(results_dir, "cox_survival_results.rds"))
volt2 <- readRDS(file.path(results_dir, "cox_volt2_comparison.rds"))

meta_files <- list.files(compressed_dir, pattern = "^meta_.*\\.rds$", full.names = TRUE)
datasets <- gsub("meta_|\\.rds", "", basename(meta_files))

hr_fmt <- function(coef, se) sprintf("%.2f (%.2f-%.2f)", exp(coef), exp(coef - 1.96*se), exp(coef + 1.96*se))

for (ds in datasets) {
  meta <- readRDS(file.path(compressed_dir, sprintf("meta_%s.rds", ds)))
  cat("\n", strrep("=", 60), "\n ", ds, "\n", strrep("=", 60), "\n")
  cat(sprintf("%s | DPE %d iter %d | %d active factors\n",
              meta$condition, meta$dpe_stage, meta$iter, meta$active_factors))

  # EDSS
  sig_edss <- edss %>% filter(dataset == ds, p_value_bh < 0.05)
  cat(sprintf("\nEDSS correlations (BH<0.05): %d\n", nrow(sig_edss)))
  print(sig_edss %>% arrange(p_value_bh) %>% head(10) %>%
          mutate(p_value_bh = sprintf("%.4f", p_value_bh), correlation = round(correlation, 3)) %>%
          dplyr::select(month, clinical_var, factor, correlation, p_value_bh) %>% as.data.frame(),
        row.names = FALSE)

  # Cox
  sig_cox <- cox %>% filter(dataset == ds, p_value < 0.05, grepl("^V", variable))
  cat(sprintf("\nCox significant factors: %d\n", nrow(sig_cox)))
  print(sig_cox %>% arrange(p_value) %>%
          mutate(HR = hr_fmt(coef, se), p_value = sprintf("%.4f", p_value),
                 direction = ifelse(coef < 0, "protective", "risk"),
                 lesion_assoc = ifelse(avg_loading > 0, "pro-lesion", "anti-lesion"),
                 avg_loading = sprintf("%.4f", avg_loading)) %>%
          dplyr::select(variable, HR, p_value, direction, lesion_assoc, avg_loading) %>% as.data.frame(),
        row.names = FALSE)
  cat(sprintf("Brain map: %s/coxph_sig_factors_%s.nii.gz\n", results_dir, ds))
}

# VOLT2 comparison
v <- volt2$results %>% filter(variable == "VOLT2")
cat("\n", strrep("=", 60), "\n VOLT2 baseline\n", strrep("=", 60), "\n")
cat(sprintf("HR = %s, p = %.4f\n", hr_fmt(v$coef, v$se), v$p_value))

# Summary table
cat("\n", strrep("=", 60), "\n Comparison\n", strrep("=", 60), "\n")
summary_df <- data.frame(
  dataset = datasets,
  active_K = sapply(datasets, function(ds) readRDS(file.path(compressed_dir, sprintf("meta_%s.rds", ds)))$active_factors),
  edss_sig = sapply(datasets, function(ds) sum(edss$dataset == ds & edss$p_value_bh < 0.05)),
  cox_sig = sapply(datasets, function(ds) sum(cox$dataset == ds & cox$p_value < 0.05 & grepl("^V", cox$variable)))
)
print(summary_df, row.names = FALSE)
