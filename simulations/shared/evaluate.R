evaluate_run = function(results, lambda, method_name, r, time_taken) {
  final = results[[length(results)]]
  w = final$w > 0.5

  sorted_lambda_corr = final$lambda_corr
  sorted_lambda_corr[!w] = 0
  # Sign-correct: flip columns so the largest loading (by magnitude) is positive
  max_loading = apply(abs(sorted_lambda_corr), 2, which.max)
  signs = sign(sorted_lambda_corr[cbind(max_loading, seq_len(ncol(sorted_lambda_corr)))])
  sorted_lambda_corr = t(signs * t(sorted_lambda_corr))

  col_order = order(colSums(w), decreasing = TRUE)
  sorted_lambda_corr = sorted_lambda_corr[, col_order]
  w = w[, col_order]

  expanded_lambda = matrix(0, nrow(w), ncol(w))
  expanded_lambda[, seq_len(ncol(lambda))] = (lambda != 0)

  expanded_lambda_cor = expanded_lambda / sqrt(1 + diag(lambda %*% t(lambda)))
  lambda_cor_mse = mean((expanded_lambda_cor - sorted_lambda_corr)^2)

  true_k = sum(colSums(lambda > 0) > 0)
  estimated_k = sum(colSums(w) > 0)

  iterations = tryCatch(
    sum(sapply(results, function(x) x$iterations)),
    error = function(e) NA_integer_
  )

  data.frame(
    method = method_name,
    structure_recovered = all(expanded_lambda == w),
    K_correct = (true_k == estimated_k),
    lambda_cor_mse = lambda_cor_mse,
    realisation = r,
    time_taken = time_taken,
    iterations = iterations
  )
}

print_evaluation = function(results, lambda, method_name, r, N = NULL) {
  final = results[[length(results)]]
  w = final$w
  lc = final$lambda_corr
  if (!is.matrix(w)) { w = matrix(w, ncol = 1) }
  if (!is.matrix(lc)) { lc = matrix(lc, ncol = 1) }
  w = w > 0.5

  lc[!w] = 0
  # Sign-correct: flip columns so the largest loading (by magnitude) is positive
  max_loading = apply(abs(lc), 2, which.max)
  signs = sign(lc[cbind(max_loading, seq_len(ncol(lc)))])
  # Multiply each column by its sign (transposes needed since * is row-wise)
  lc = t(signs * t(lc))
  ord = order(colSums(w), decreasing = TRUE)
  lc = lc[, ord, drop = FALSE]
  w = w[, ord, drop = FALSE]

  K_est = ncol(w)
  K_max = max(K_est, ncol(lambda))
  expanded = matrix(0, nrow(lambda), K_max)
  expanded[, seq_len(ncol(lambda))] = (lambda != 0)
  w_full = matrix(FALSE, nrow(lambda), K_max)
  w_full[, seq_len(K_est)] = w
  lc_full = matrix(0, nrow(lambda), K_max)
  lc_full[, seq_len(K_est)] = lc

  header = paste0("Method: ", method_name, " | r: ", r)
  if (!is.null(N)) header = paste0(header, " | N: ", N)
  cat(header, "\n\nTRUE:\n")
  print(expanded)
  cat("\nESTIMATED (thresholded):\n")
  print(round(lc_full, 3))
  cat("\nK:", sum(colSums(w) > 0), "/", ncol(lambda),
      "| struct:", all(expanded == w_full), "\n")
}
