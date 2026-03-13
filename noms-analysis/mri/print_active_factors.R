base_path <- "/data/users/uu85g9/factor-analysis/noms/mri"  # previous: same base but scanned 10 conditions under runs_jan18_v0_0.1

conditions <- c("dpe_random_pxsigma", "dpe_random_pxsigma_restart_dpe3")  # paper model conditions only

for (cond in conditions) {
  res_path <- paste0(base_path, "/", cond, "/res.rds")

  if (file.exists(res_path)) {
    res <- readRDS(res_path)
    if (is.list(res) && is.list(res[[1]]) && "lambda" %in% names(res[[1]])) {
      res <- res[[length(res)]]
    }
    src <- "res.rds"
  } else {
    files <- list.files(paste0(base_path, "/", cond, "/intermediate_results"), pattern = "iter.*\\.rds$", full.names = TRUE)
    if (length(files) == 0) { cat(cond, ": no results\n"); next }
    
    dpe_nums <- as.numeric(gsub(".*dpe([0-9]+)_iter.*", "\\1", files))
    iter_nums <- as.numeric(gsub(".*iter([0-9]+)\\.rds$", "\\1", files))
    files <- files[order(dpe_nums, iter_nums)]
    res <- readRDS(files[length(files)])
    src <- basename(files[length(files)])
  }

  K <- sum(colSums(res$w > 0.5) > 0)
  lambda_max <- max(abs(res$lambda))
  sigma_info <- if (!is.null(res$sigma)) paste0(", mean(σ²)=", round(mean(res$sigma), 4)) else ""
  cat(cond, "(", src, "): K =", K, ", max|Λ| =", round(lambda_max, 4), sigma_info, "\n")
}
