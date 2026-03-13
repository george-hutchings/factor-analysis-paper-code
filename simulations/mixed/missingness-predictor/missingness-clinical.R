## Fit missingness models and marginal ECDFs from NO.MS clinical data.
## Saves an RDS used by simulate_missing_data_functions.R to generate
## realistic simulated datasets with matching missingness patterns.
library(pROC)
library(caret)

save_path = "/data/users/uu85g9/factor-analysis/simulations/mixed/missingness_analysis.rds"
set.seed(1234)

# Load real clinical data
whole_data = read.csv("/data/ms/processed/clinical/projects/fahmm/interim_tables/03_baseline_results.csv")
habib_cols = c("RELAPSE", "EDSS", "T25FWM", "HPT9M", "PASAT", "NUMGDT1", "VOLT2", "NBV2")
Y = as.data.frame(whole_data[habib_cols])
missing_mask = is.na(Y)

# Column types
discrete_cols = c("RELAPSE")
ordinal_cols = c("EDSS", "NUMGDT1")
for (col in discrete_cols) Y[[col]] = as.factor(Y[[col]])
for (col in ordinal_cols) Y[[col]] = as.ordered(Y[[col]])

# Storage
model_list = ecdf_list = levels_list = roc_values = vector("list", ncol(Y))
names(model_list) = names(ecdf_list) = names(levels_list) = names(roc_values) = names(Y)

for (i in seq_len(ncol(Y))) {
  ecdf_list[[i]] = ecdf(Y[[i]])
  if (!is.numeric(Y[[i]])) levels_list[[i]] = levels(Y[[i]])

  if (sum(missing_mask[, i]) == 0) next

  # P(variable i missing | other variables)
  complete_idx = complete.cases(Y[, -i])
  target = missing_mask[complete_idx, i]
  predictors = as.data.frame(Y[complete_idx, -i])
  predictors$target = target

  # Diagnostic: 90/10 split stratified on NUMGDT1 (most levels of any ordinal predictor)
  # to ensure rare factor levels appear in both train and test sets
  train_idx = createDataPartition(Y$NUMGDT1, p = 0.9, list = FALSE)
  diag_model = glm(target ~ ., data = predictors[train_idx, ], family = binomial)
  preds = predict(diag_model, newdata = predictors[-train_idx, ], type = "response")
  roc_values[[i]] = roc(predictors$target[-train_idx], preds)

  # Final model: train on all data
  model_list[[i]] = glm(target ~ ., data = predictors, family = binomial)
}

# Diagnostics
cat("N:", nrow(Y), "\n\nMissing values per variable:\n")
for (v in names(Y)) cat("  ", v, ": ", sum(missing_mask[, v]), "\n", sep = "")
cat("\nAUC per variable:\n")
for (v in names(roc_values)) {
  if (!is.null(roc_values[[v]])) cat("  ", v, ": ", round(auc(roc_values[[v]]), 3), "\n", sep = "")
}
cat("\nModel summaries:\n")
for (v in names(model_list)) {
  if (!is.null(model_list[[v]])) { cat("\n---", v, "---\n"); print(summary(model_list[[v]])) }
}

saveRDS(list(
  model_list = model_list, ecdf_list = ecdf_list,
  levels_list = levels_list, discrete_cols = discrete_cols, ordinal_cols = ordinal_cols
), save_path)
cat("\nSaved to", save_path, "\n")
