### Evaluation for the high-dimensional (Rockova & George 2016) continuous scenario.
###
### Identical to shared/evaluate.R's evaluate_run(), EXCEPT the estimated factors
### are matched to the true factors by summed support overlap (align_to_truth)
### rather than by colSums size. R&G's five true factors are all 500 wide (equal
### size, 136 overlap), so they cannot be ordered by nonzero count as the standard
### evaluate_run() does; they are told apart by *where* they load, not how much.

# Match estimated columns to true factors by how many of each column's selected
# variables (w > 0.5) fall in each true factor's rows. crossprod(w_thr, T) is
# t(w_thr) %*% T: for 0/1 matrices, entry (i, j) = number of rows active in both,
# i.e. the summed-row overlap of estimated column i with true factor j. Greedy
# one-to-one assignment (largest overlap first) is optimal for this design because
# each factor overlaps its own true region by ~500 but any other region by <= 136.
align_to_truth = function(w_thr, true_support) {
  Kfit = ncol(w_thr)
  Ktrue = ncol(true_support)
  overlap = crossprod(w_thr, true_support)        # Kfit x Ktrue: shared active vars
  matched = integer(Ktrue)
  avail = rep(TRUE, Kfit)
  for (s in seq_len(Ktrue)) {                     # take the best (est, true) pair each round
    o = overlap
    o[!avail, ] = -1L
    o[, matched != 0] = -1L
    ij = which(o == max(o), arr.ind = TRUE)[1, ]  # best remaining pair
    matched[ij[2]] = ij[1]
    avail[ij[1]] = FALSE
  }
  c(matched, which(avail))                        # true-matched columns first, rest after
}

evaluate_run = function(results, lambda, method_name, r, time_taken) {
  final = results[[length(results)]]
  w = final$w > 0.5

  sorted_lambda_corr = final$lambda_corr
  sorted_lambda_corr[!w] = 0
  # Sign-correct: flip columns so the largest loading (by magnitude) is positive
  max_loading = apply(abs(sorted_lambda_corr), 2, which.max)
  signs = sign(sorted_lambda_corr[cbind(max_loading, seq_len(ncol(sorted_lambda_corr)))])
  # See shared/evaluate.R: sign(0) on a fully-thresholded column would zero lambda_raw.
  signs[signs == 0] = 1
  sorted_lambda_corr = t(signs * t(sorted_lambda_corr))

  # Same signs and column order, but before thresholding, so the run can be
  # re-analysed at a different threshold without refitting
  lambda_raw = t(signs * t(final$lambda_corr))
  w_prob = final$w

  # --- the only change from shared/evaluate.R: match by support overlap, not size ---
  # The K-known methods (PCA) have no support, so w is all NA. Match those on the
  # loading magnitude in each true region instead, which is the continuous analogue
  # of the support overlap and keeps their MSE comparable.
  w_for_align = if (all(is.na(w))) abs(sorted_lambda_corr) else w
  col_order = align_to_truth(w_for_align, (lambda != 0))
  sorted_lambda_corr = sorted_lambda_corr[, col_order]
  w = w[, col_order]
  lambda_raw = lambda_raw[, col_order]
  w_prob = w_prob[, col_order]

  expanded_lambda = matrix(0, nrow(w), ncol(w))
  expanded_lambda[, seq_len(ncol(lambda))] = (lambda != 0)

  # rowSums(lambda^2) == diag(lambda %*% t(lambda)) but avoids a D x D matrix (D ~ 2000)
  expanded_lambda_cor = expanded_lambda / sqrt(1 + rowSums(lambda^2))
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
