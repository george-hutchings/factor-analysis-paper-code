### Loading paths along the DPE schedule, simulation Scenario 1: every element of
### the loading matrix traced across the spike variances, showing the large ones
### stabilise while the negligible ones are shrunk to zero.
###
###   Rscript simulations/dpe_path.R
###     -> simulations/dpe_maintext_pca_r1.{rds,pdf}
###
### Fits on first call, caches to .rds, then plots. Nothing on disk could produce
### this: continuous/run.R never passes save_dir, and shared/evaluate.R:2
### keeps only the final rung. Seeds match the paper (data 12345+r, fit 12345),
### so the realisation is the same dataset the published metrics were run on.

library(ggplot2)

## Paths resolve from the script's own location, so it runs from anywhere:
## the algorithm is one level up and the output is written next to this script.
HERE <- dirname(normalizePath(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1])))
ALGO <- normalizePath(file.path(HERE, "..", "mcem_fa_algorithm.R"))
OUT  <- HERE
f <- function(...) file.path(OUT, paste0(...))

## The paper's Scenario 1 on the main-text schedule, as in continuous/run.R.
## PCA init (lambda NULL) is what the paper uses.
N <- 200; K <- 10; R <- 1
V0S <- 10^-seq(1, 4, length.out = 4)

TRUTH <- 1/sqrt(2)                    # see "truth" below
RED <- "#d7191c"; BLUE <- "#2c7bb6"   # validated: deutan dE 23.1, normal 32.1
INK <- "grey20"; MUTED <- "grey45"; HAIR <- "grey88"

L <- matrix(0, 12, 3); L[1:6,1] <- 1; L[7:10,2] <- 1; L[11:12,3] <- 1  # 6/4/2 blocks
stem <- sprintf("dpe_maintext_pca_r%d", R)

## Must be source()d at top level, never inside a function: the algorithm finds
## its .cpp via dirname(sys.frame(1)$ofile), which is NULL inside a call frame.
if (!file.exists(f(stem, ".rds"))) source(ALGO)

fit <- function(stem) {
  set.seed(12345 + R)
  Y <- tcrossprod(matrix(rnorm(N*3), N, 3), L) + matrix(rnorm(N*12), N, 12)
  ck <- f(".ck_", stem, "/"); dir.create(ck, showWarnings = FALSE)
  set.seed(12345)
  ## save_dir makes mcem_algorithm write its iter0 checkpoint, which is the
  ## "init" point on the plot; save_every huge suppresses all the others.
  res <- dpe_func(V0S, Y = Y, K = K, do_stick_breaking = TRUE, px_sigma = FALSE,
                  px_rotate = FALSE, varimax_every = -1, conv_param = TRUE,
                  nburn = 100, n_mcmc = 100, tol = 0.016, maxiter = 2000,
                  verbose = -1, save_dir = ck, save_every = 1e6)
  i0 <- readRDS(file.path(ck, "dpe1_iter0.rds"))
  saveRDS(list(stages = lapply(res, function(f) f[c("v0","w","lambda_corr")]),
               init = i0[c("w","lambda_corr")]), f(stem, ".rds"))
  unlink(ck, recursive = TRUE)
}

## Two silent traps. dpe_func varimax-rotates between EVERY rung regardless
## of varimax_every=-1, so column k is a different factor at each rung and its
## sign is free. Match every rung to the final one, whose own signs follow
## evaluate.R:8-14 (largest entry positive); display order is evaluate.R:21.
sgn <- function(M) { s <- sign(M[cbind(apply(abs(M),2,which.max), 1:ncol(M))])
                     t(ifelse(s == 0, 1, s) * t(M)) }

draw <- function(stem) {
  o <- readRDS(f(stem, ".rds")); st <- o$stages; fin <- st[[length(st)]]
  ref <- sgn(fin$lambda_corr); disp <- order(colSums(fin$w > 0.5), decreasing = TRUE)
  match_cols <- function(X) {                 # greedy |inner product|, best pair first
    M <- abs(crossprod(X, ref)); ord <- integer(K); uX <- uR <- logical(K)
    for (i in 1:K) { m <- M; m[uX,] <- -Inf; m[,uR] <- -Inf; j <- which.max(m)
      p <- (j-1) %% K + 1; q <- (j-1) %/% K + 1; ord[q] <- p; uX[p] <- uR[q] <- TRUE }
    list(o = ord, s = vapply(1:K, function(q)
      if (sum(X[,ord[q]] * ref[,q]) < 0) -1 else 1, 0))
  }

  xs <- -log10(vapply(st, function(s) s$v0, 0))
  xin <- min(xs) - .18*diff(range(xs))        # init has no v0, sits left
  d <- do.call(rbind, Map(function(s, x) { m <- match_cols(s$lambda_corr)
    data.frame(x = x, id = seq_len(12*K),
               value = as.vector(t(m$s * t(s$lambda_corr[, m$o]))[, disp]),
               incl  = as.vector(s$w[, m$o][, disp] > 0.5)) },
    c(list(o$init), st), c(xin, xs)))
  d <- d[order(d$id, d$x), ]

  ## Each segment takes the colour of the rung it arrives at; mapping colour
  ## onto geom_line would split the line into groups rather than recolour it.
  ## Alpha on segments only, so 120 traces can overlap without the zero band
  ## going solid; points stay opaque and carry the per-rung colour.
  sg <- do.call(rbind, lapply(split(d, d$id), function(z)
    data.frame(x = head(z$x,-1), xend = tail(z$x,-1), y = head(z$value,-1),
               yend = tail(z$value,-1), incl = tail(z$incl,-1))))

  ## Truth: every variable loads on one factor with raw loading 1, so
  ## rowSums(lambda^2) = 1 and eq. (6) sends the whole support to 1/sqrt(2);
  ## the other 108 entries are truly 0. Dashed because these are thresholds, not
  ## grid, drawn as segments so they stop at the last rung, and keyed rather
  ## than captioned so dashed is explained where red and blue are.
  tl <- data.frame(x = xin, xend = max(xs), y = c(0, TRUTH), yend = c(0, TRUTH),
                   key = "true value")

  p <- ggplot() +
    geom_segment(data = tl, aes(x, y, xend = xend, yend = yend, linetype = key),
                 colour = MUTED, linewidth = .35) +
    scale_linetype_manual(values = c(`true value` = "22"), name = NULL) +
    geom_vline(xintercept = (xin + min(xs))/2, colour = HAIR, linewidth = .35) +
    geom_segment(data = sg, aes(x, y, xend = xend, yend = yend, colour = incl),
                 linewidth = .35, alpha = .55, show.legend = FALSE) +
    ## a white surface ring, not a border, separates overlapping markers
    geom_point(data = d, aes(x, value, fill = incl), shape = 21, colour = "white",
               stroke = .25, size = 1.5) +
    scale_colour_manual(values = c(`TRUE` = RED, `FALSE` = BLUE), guide = "none") +
    scale_fill_manual(values = c(`TRUE` = RED, `FALSE` = BLUE), name = "assigned to",
                      breaks = c(TRUE, FALSE), labels = c("slab (w > 0.5)", "spike")) +
    scale_x_continuous(breaks = c(xin, xs), expand = expansion(mult = c(.05, .05)),
                       labels = parse(text = c("init", sprintf("10^-%d", xs)))) +
    scale_y_continuous(breaks = sort(c(seq(-0.75, 0.5, by = .25), TRUTH)),
                       labels = function(b) ifelse(abs(b - TRUTH) < 1e-8, "0.71",
                                                   formatC(b, format = "g"))) +
    guides(fill = guide_legend(order = 1, override.aes = list(size = 2.6)),
           linetype = guide_legend(order = 2, override.aes = list(colour = MUTED))) +
    ## "schedule", not "ladder": the manuscript's own word for the v0 sequence.
    ## No in-plot title and no equation number anywhere: the caption describes
    ## the plot, and the y-axis "normalised loading" names the scale that the
    ## Methods (eq. 6) fixes as the paper-wide reporting convention.
    labs(x = expression("spike variance"~v[0]~"(decreasing"%->%")"),
         y = "normalised loading") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          panel.grid.major.x = element_line(colour = HAIR, linewidth = .3),
          ## no baseline rule: it would trail past the last rung. The rung
          ## gridlines already mark the positions.
          axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text = element_text(colour = MUTED, size = 8),
          axis.title = element_text(colour = INK, size = 9),
          axis.title.x = element_text(margin = margin(t = 7)),
          axis.title.y = element_text(margin = margin(r = 7)),
          legend.position = "bottom", legend.margin = margin(t = 2),
          legend.key.height = unit(9, "pt"),
          legend.text = element_text(colour = INK, size = 8),
          legend.title = element_text(colour = INK, size = 8),
          plot.margin = margin(10, 12, 8, 10))

  ggsave(f(stem, ".pdf"), p, width = 6.4, height = 4.3)
}

if (!file.exists(f(stem, ".rds"))) fit(stem)
draw(stem)
