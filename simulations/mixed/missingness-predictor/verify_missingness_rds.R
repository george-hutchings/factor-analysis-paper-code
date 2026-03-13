## Verify that the newly generated missingness_analysis.rds is consistent
## with the old one from noms-analysis/clinical/simulations.
##
## Differences expected:
##   - Old model was trained on 90% of data (train/test split)
##   - New model is trained on 100% of data
##   - Different seed for createDataPartition
## So: ECDFs, levels, metadata should be IDENTICAL.
##     Model coefficients should be SIMILAR but not exact.
##
## Usage (on Novartis cluster, after running missingness-clinical.R):
##   Rscript verify_missingness_rds.R

old_path = "/data/users/uu85g9/factor-analysis/noms/clinical/missingness_analysis.rds"
new_path = "/data/users/uu85g9/factor-analysis/simulations/mixed/missingness_analysis.rds"

old = readRDS(old_path)
new = readRDS(new_path)

# --- 1. Same top-level structure ---
cat("Same list names:", identical(names(old), names(new)), "\n")

# --- 2. Same metadata (deterministic, should be identical) ---
cat("discrete_cols identical:", identical(old$discrete_cols, new$discrete_cols), "\n")
cat("ordinal_cols identical:", identical(old$ordinal_cols, new$ordinal_cols), "\n")
cat("levels_list identical:", identical(old$levels_list, new$levels_list), "\n")

# --- 3. ECDFs (deterministic from data, should give same values) ---
# identical() is too strict for ecdf objects (different closures/environments)
ecdf_match = sapply(names(old$ecdf_list), function(v) {
  grid = knots(old$ecdf_list[[v]])
  all(old$ecdf_list[[v]](grid) == new$ecdf_list[[v]](grid))
})
cat("\nECDFs match per variable:\n")
print(ecdf_match)

# --- 4. Same variables have missingness models ---
old_has_model = !sapply(old$model_list, is.null)
new_has_model = !sapply(new$model_list, is.null)
cat("\nSame variables modelled:", identical(old_has_model, new_has_model), "\n")

# --- 5. Compare model coefficients (similar but not exact) ---
cat("\nModel coefficient comparison (old=90% train, new=100% train):\n")
for (v in names(old$model_list)) {
  if (!is.null(old$model_list[[v]])) {
    cat("\n", v, ":\n")
    print(rbind(old = coef(old$model_list[[v]]), new = coef(new$model_list[[v]])))
  }
}
