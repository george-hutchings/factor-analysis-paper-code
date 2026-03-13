# load data
# SET THIS: must match save_path in missingness-clinical.R
missingness_rds_path = "/data/users/uu85g9/factor-analysis/simulations/mixed/missingness_analysis.rds"
# unpack list
list2env(readRDS(missingness_rds_path), .GlobalEnv)

# generation of data

simulate_complete_data = function(ecdf_list, lambda, N, ordinal_cols, discrete_cols){
  types = ifelse(names(ecdf_list) %in% c(discrete_cols, ordinal_cols), 1, 5)
  
  D = nrow(lambda)
  stopifnot(D==length(ecdf_list))
  
  true_K = ncol(lambda)
  
  eta = matrix(rnorm(N*true_K), N, true_K)
  Z = tcrossprod(eta, lambda) + matrix(rnorm(N*D), N, D)
  
  Y <- matrix(NA, nrow = N, ncol = D)
  for (d in seq_len(D)) {
    # convert from normal to uniform
    mu = 0
    sd = sqrt(sum(lambda[d, ]**2) + 1)
    rand_uniform = pnorm(Z[,d], mean=mu, sd=sd)
    Y[,d] = quantile(ecdf_list[[d]], rand_uniform, type = types[d])
    
  }
  colnames(Y) = names(ecdf_list)
  
  
  
  # Convert columns to appropriate types
  convert_columns <- function(df, cols, type) {
    for (col in cols) {
      df[[col]] <- type(df[[col]])
    }
    return(df)
  }
  # Convert columns ensuring they have same name levels as in proportion_summary
  Y <- as.data.frame(Y)
  Y <- convert_columns(Y, discrete_cols, as.factor)
  Y <- convert_columns(Y, ordinal_cols, as.ordered)
  
  
  # Function to set levels for discrete and ordinal columns
  set_levels <- function(df, cols, levels_list) {
    for (col in cols) {
      if (col %in% names(levels_list)) {
        levels(df[[col]]) <- levels_list[[col]]
      }
    }
    return(df)
  }
  Y <- set_levels(Y, c(discrete_cols, ordinal_cols), levels_list)
  
  return(Y)
}

simulate_missingness = function(Y, model_list){
  D = ncol(Y)
  N = nrow(Y)
  ## apply missingness model (glm object) to get missingness probabilities (if it is NULL do nothing)
  missingness_prob = matrix(0, N, D)
  for(d in 1:D){
    if(!is.null(model_list[[d]])){
      missingness_prob[,d] = predict(model_list[[d]], newdata=data.frame(Y), type='response')
    }
  }
  # Apply missingness model to generate missing data
  missing_mask = matrix(runif(N*D) < missingness_prob, N, D)
}

## generate data

# lambda = matrix(0, 8, 3)
# lambda[1:3, 1] = 1
# lambda[4:7, 2] = 1
# lambda[7:8, 3] = 1
# 
# N=10000
# Y_complete = simulate_complete_data(ecdf_list, lambda, N, ordinal_cols, discrete_cols)
# 
# missing_mask = simulate_missingness(Y_complete, model_list)
# Y_complete = as.matrix(Y_complete)

