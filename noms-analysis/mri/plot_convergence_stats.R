# Plot convergence statistics for factor analysis runs
# Usage:
#   plot_convergence("/path/to/condition")           # single path
#   plot_convergence("/path/to/condition", save=FALSE)  # display only
#   plot_all_conditions()                            # all jan18 conditions
library(ggplot2)
library(gridExtra)

# Helper to get Ereconstruction from checkpoints if missing
get_Erecon_from_checkpoints <- function(cond_path) {
  ckpt_dir <- paste0(cond_path, "/intermediate_results")
  if (!dir.exists(ckpt_dir)) return(NULL)
  files <- list.files(ckpt_dir, pattern = "iter.*\\.rds$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  # Sort by DPE stage then iteration to get truly latest
  dpe_nums <- as.numeric(gsub(".*dpe([0-9]+)_iter.*", "\\1", files))
  iter_nums <- as.numeric(gsub(".*iter([0-9]+)\\.rds$", "\\1", files))
  files <- files[order(dpe_nums, iter_nums)]
  ckpt <- readRDS(files[length(files)])
  ckpt$Ereconstruction
}

#' Plot convergence statistics for a single run
#' @param cond_path Path to the condition directory (containing res.rds or intermediate_results/)
#' @param save If TRUE, save PNG to cond_path. If FALSE, just display.
#' @param label Optional label for plot title (defaults to basename of path)
plot_convergence <- function(cond_path, save = TRUE, label = NULL) {
  if (is.null(label)) label <- basename(cond_path)
  res_path <- paste0(cond_path, "/res.rds")

  # Load most recent result (res.rds or latest checkpoint)
  if (file.exists(res_path)) {
    res <- readRDS(res_path)
    # Handle DPE list structure - take final stage
    if (is.list(res) && is.list(res[[1]]) && "lambda" %in% names(res[[1]])) {
      dpe_stage <- length(res)
      res <- res[[dpe_stage]]
      suffix <- paste0("_dpe", dpe_stage)
    } else {
      suffix <- ""
    }
    src <- "res.rds"
  } else {
    files <- list.files(paste0(cond_path, "/intermediate_results"),
                        pattern = "iter.*\\.rds$", full.names = TRUE)
    if (length(files) == 0) {
      cat(label, ": no results found\n")
      return(invisible(NULL))
    }
    # Sort by DPE stage then iteration
    dpe_nums <- as.numeric(gsub(".*dpe([0-9]+)_iter.*", "\\1", files))
    iter_nums <- as.numeric(gsub(".*iter([0-9]+)\\.rds$", "\\1", files))
    files <- files[order(dpe_nums, iter_nums)]
    res <- readRDS(files[length(files)])
    src <- basename(files[length(files)])
    suffix <- paste0("_", gsub("\\.rds$", "", src))
  }

  # Check for convergence metrics
  if (is.null(res$parameter_diff) || is.null(res$logpriors)) {
    cat(label, "(", src, "): missing convergence metrics\n")
    return(invisible(NULL))
  }

  n_iter <- length(res$parameter_diff)

  # Get Ereconstruction (fallback to checkpoints if missing)
  Erecon <- res$Ereconstruction
  if (is.null(Erecon)) Erecon <- get_Erecon_from_checkpoints(cond_path)
  if (!is.null(Erecon)) Erecon <- Erecon[1:min(n_iter, length(Erecon))]
  has_Erecon <- !is.null(Erecon) && length(Erecon) > 0

  # Build data frame
  df <- data.frame(
    iter = seq_len(n_iter),
    param_diff = res$parameter_diff,
    logprior = res$logpriors
  )
  if (has_Erecon) df$Ereconstruction <- Erecon

  # Create plots
  p1 <- ggplot(df, aes(iter, param_diff)) +
    geom_line() +
    scale_y_log10() +
    labs(title = "Parameter diff", x = "Iteration",
         y = expression(paste("max ", group("|", group("|", Lambda^{OLD}, "|") - group("|", Lambda, "|"), "|")))) +
    theme_minimal()

  p2 <- ggplot(df, aes(iter, logprior)) +
    geom_line() +
    labs(title = "Log prior (objective)", x = "Iteration", y = "Log prior") +
    theme_minimal()

  if (has_Erecon) {
    p3 <- ggplot(df, aes(iter, Ereconstruction)) +
      geom_line() +
      labs(title = "E[reconstruction]", x = "Iteration", y = "Ereconstruction") +
      theme_minimal()
    p <- grid.arrange(p1, p2, p3, ncol = 3, top = paste0(label, " (", src, ")"))
    if (save) {
      out_file <- paste0(cond_path, "/convergence_", label, suffix, ".png")
      ggsave(out_file, p, width = 15, height = 4)
      cat(label, "(", src, "): saved", out_file, "\n")
    }
  } else {
    p <- grid.arrange(p1, p2, ncol = 2, top = paste0(label, " (", src, ")"))
    if (save) {
      out_file <- paste0(cond_path, "/convergence_", label, suffix, ".png")
      ggsave(out_file, p, width = 10, height = 4)
      cat(label, "(", src, "): saved", out_file, "\n")
    }
  }

  invisible(p)
}

#' Plot convergence for paper model conditions
#' (previously pointed at .../mri/runs_jan18_v0_0.1 with 10 experimental conditions)
plot_all_conditions <- function() {
  base_path <- "/data/users/uu85g9/factor-analysis/noms/mri"  # previous: .../mri/runs_jan18_v0_0.1
  conditions <- c("dpe_random_pxsigma", "dpe_random_pxsigma_restart_dpe3")  # paper model conditions only

  for (cond in conditions) {
    plot_convergence(paste0(base_path, "/", cond), save = TRUE, label = cond)
  }
  cat("\nDone!\n")
}
