evaluate_run = function(results, lambda, method_name, r, time_taken) {
  final = results[[length(results)]]
  w = final$w > 0.5

  sorted_lambda_corr = final$lambda_corr
  sorted_lambda_corr[!w] = 0
  # Sign-correct: flip columns so the largest loading (by magnitude) is positive
  max_loading = apply(abs(sorted_lambda_corr), 2, which.max)
  signs = sign(sorted_lambda_corr[cbind(max_loading, seq_len(ncol(sorted_lambda_corr)))])
  # A column thresholded entirely to zero gives sign(0) = 0, which would then multiply
  # lambda_raw's copy of it by zero and destroy exactly the sub-threshold values that
  # field exists to preserve. Such a column has no orientation to fix, so leave it alone.
  signs[signs == 0] = 1
  sorted_lambda_corr = t(signs * t(sorted_lambda_corr))

  # Same signs and column order, but before thresholding, so the run can be
  # re-analysed at a different threshold without refitting
  lambda_raw = t(signs * t(final$lambda_corr))
  w_prob = final$w

  col_order = order(colSums(w), decreasing = TRUE)
  sorted_lambda_corr = sorted_lambda_corr[, col_order]
  w = w[, col_order]
  lambda_raw = lambda_raw[, col_order]
  w_prob = w_prob[, col_order]

  expanded_lambda = matrix(0, nrow(w), ncol(w))
  expanded_lambda[, seq_len(ncol(lambda))] = (lambda != 0)

  expanded_lambda_cor = expanded_lambda / sqrt(1 + diag(lambda %*% t(lambda)))
  lambda_cor_mse = mean((expanded_lambda_cor - sorted_lambda_corr)^2)

  true_k = sum(colSums(lambda > 0) > 0)
  estimated_k = sum(colSums(w) > 0)

  # Support recovery errors (NA for the K-known methods, which have no support)
  false_positives = sum(w & !expanded_lambda)
  false_negatives = sum(!w & expanded_lambda)

  iterations = tryCatch(
    sum(sapply(results, function(x) x$iterations)),
    error = function(e) NA_integer_
  )

  data.frame(
    method = method_name,
    structure_recovered = all(expanded_lambda == w),
    K_correct = (true_k == estimated_k),
    lambda_cor_mse = lambda_cor_mse,
    false_positives = false_positives,
    false_negatives = false_negatives,
    realisation = r,
    time_taken = time_taken,
    iterations = iterations,
    # Kept so the per-element MSE can be aggregated across realisations, and so
    # any other loading-level summary can be made without rerunning the fits
    lambda_hat = I(list(sorted_lambda_corr)),
    lambda_raw = I(list(lambda_raw)),
    w_prob = I(list(w_prob)),
    lambda_true = I(list(expanded_lambda_cor))
  )
}
