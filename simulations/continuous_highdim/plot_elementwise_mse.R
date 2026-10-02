### Per-element MSE of the estimated loadings, high-dimensional scenario.
###   Rscript plot_elementwise_mse.R <elementwise_mse.rds> <out.pdf> [methods]
###     methods: optional comma-separated list of internal method names, in the
###     panel order wanted. If omitted, all methods are shown in stored order.
###
### Same colour scale, labels, shared legend and patchwork layout as
### simulations/plot_elementwise_mse.R. Two things change, both forced by
### D = 1956:
###   geom_raster rather than geom_tile. 1956 x 20 = 39120 outlined rectangles a
###   panel is unworkable in a vector PDF, and at this row count the grey cell
###   borders cover the cells entirely.
###   The true support is outlined as one rectangle per contiguous block instead
###   of per cell, which is the same information as a 500-row-tall outline
###   rather than 500 stacked ones.

library(ggplot2)
library(reshape2)
library(patchwork)

# Display labels matching the manuscript tables; unlisted methods keep their key.
LABELS = c(
  "our method"       = "MIDFA",
  "rockova 2016"     = "Rockova & George (2016)",
  "pca K known"      = "PCA (K known)",
  "factanal K known" = "factanal (K known)"
)

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript plot_elementwise_mse.R <elementwise_mse.rds> <out.pdf> [methods]")
res = readRDS(args[1])

if (length(args) >= 3) {
  want = trimws(strsplit(args[3], ",")[[1]])
  missing = setdiff(want, names(res))
  if (length(missing)) stop("methods not found in RDS: ", paste(missing, collapse = ", "))
  res = res[want]
}

D     = nrow(res[[1]]$mse)
K_max = max(sapply(res, function(x) ncol(x$mse)))

# The K-known methods estimate fewer columns than the truncation level, so pad
# them to K_max with NA. Those cells are drawn blank rather than as zero error,
# which would wrongly credit them with finding the empty factors.
pad = function(mat) {
  out = matrix(NA_real_, D, K_max)
  out[, seq_len(ncol(mat))] = mat
  out
}

# One colour scale across every panel, so the panels are directly comparable.
# Breaks are evenly spaced on the sqrt scale, so they do not bunch up at the top.
mse_max = max(sapply(res, function(x) max(x$mse)))
mse_breaks = signif(seq(0, sqrt(mse_max), length.out = 5)^2, 2)
mse_lim = mse_max * 1.04  # slight headroom so the top break label is not clipped

# Contiguous runs of true-support rows within each column, as rectangles.
support_blocks = function(true) {
  out = NULL
  for (k in seq_len(ncol(true))) {
    v = true[, k] != 0
    if (all(is.na(v))) next          # padded column of a K-known method
    rl = rle(v)
    hi = cumsum(rl$lengths)
    lo = hi - rl$lengths + 1
    keep = which(rl$values)
    if (length(keep))
      out = rbind(out, data.frame(Factor = k,
                                  ymin = lo[keep] - 0.5,
                                  ymax = hi[keep] + 0.5))
  }
  out
}

make_panel = function(x, title) {
  df = melt(pad(x$mse), varnames = c("Variable", "Factor"), value.name = "MSE")
  blocks = support_blocks(pad(x$true))

  p = ggplot(df, aes(x = Factor, y = Variable, fill = MSE)) +
    geom_raster() +
    scale_fill_gradient(low = "white", high = "firebrick", trans = "sqrt",
                        limits = c(0, mse_lim), breaks = mse_breaks,
                        na.value = "grey92",
                        guide = guide_colourbar(barwidth = 10, barheight = 0.6)) +
    scale_y_reverse(expand = c(0, 0), breaks = c(1, seq(500, D, by = 500), D)) +
    scale_x_continuous(expand = c(0, 0),
                       breaks = if (K_max <= 10) seq_len(K_max) else seq(5, K_max, by = 5),
                       position = "top") +
    ggtitle(title) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          panel.grid = element_blank(),
          axis.text = element_text(size = 7),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = 10))

  # Outline the blocks that are nonzero in the true loading matrix
  if (!is.null(blocks))
    p = p + geom_rect(data = blocks, inherit.aes = FALSE,
                      aes(xmin = Factor - 0.5, xmax = Factor + 0.5,
                          ymin = ymin, ymax = ymax),
                      fill = NA, colour = "black", linewidth = 0.3)
  p
}

titles = ifelse(names(res) %in% names(LABELS), LABELS[names(res)], names(res))
panels = Map(make_panel, res, titles)
fig = wrap_plots(panels, nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(args[2], fig, width = 2.6 * length(panels) + 1, height = 6)
cat("Figure saved to", args[2], "\n")
cat("Shared colour scale: 0 to", signif(mse_max, 3), "(sqrt spacing)\n")
for (i in seq_along(res))
  cat(sprintf("  %-24s %d realisations, K_est %d, mean MSE %.6f\n",
              names(res)[i], res[[i]]$n, ncol(res[[i]]$mse), mean(res[[i]]$mse)))
