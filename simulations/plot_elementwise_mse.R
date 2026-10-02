### Per-element MSE of the estimated loading matrix, one panel per method.
### Usage: Rscript plot_elementwise_mse.R <elementwise_mse.rds> <out.pdf> [methods]
###   methods: optional comma-separated list of internal method names, in the
###   panel order wanted. If omitted, all methods are shown in stored order.
### Commands for the paper, where <save_dir> is the output directory of the study given to analyse.R:
###   Fig S1: Rscript plot_elementwise_mse.R <save_dir>/elementwise_mse.rds mse_continuous.pdf 'our method,factanal K known,pca K known,rockova 2016'
###   Fig S5: Rscript plot_elementwise_mse.R <save_dir>/elementwise_mse.rds mse_binary.pdf 'our method,nmf K known,li 2023 mirt'
###   Fig S7: Rscript plot_elementwise_mse.R <save_dir>/elementwise_mse.rds mse_mixed.pdf 'full data,missing data,complete cases'
###
### The single shared linear-in-sqrt colour scale below is right for these scenarios, where the
### panels are competing methods at one sample size and their MSEs are within an order of
### magnitude of each other. It is wrong when panels span sample sizes: over N = 50 to 5000 the
### MSE falls by a factor of a thousand, and a difference-based scale leaves every panel but the
### smallest N blank white. See continuous_samplesize/plot_elementwise.R, which uses a
### log10 scale for its cross-N figure.

library(ggplot2)
library(reshape2)
library(patchwork)

# Display labels matching the manuscript tables; unlisted methods keep their key.
LABELS = c(
  "our method"     = "MIDFA",
  "factanal K known" = "factanal (K known)",
  "pca K known"    = "PCA (K known)",
  "rockova 2016"   = "Rockova & George (2016)",
  "nmf K known"    = "NMF (K known)",
  "li 2023 mirt"   = "Li et al. (2023)",
  "full data"      = "Full data (benchmark)",
  "missing data"   = "Structured missingness",
  "complete cases" = "Complete case analysis"
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

make_panel = function(x, title) {
  df = melt(pad(x$mse), varnames = c("Variable", "Factor"), value.name = "MSE")
  df$true_nonzero = melt(pad(x$true) != 0)$value
  df$Variable = factor(df$Variable, levels = rev(seq_len(D)))
  df$Factor   = factor(df$Factor, levels = seq_len(K_max))

  ggplot(df, aes(x = Factor, y = Variable, fill = MSE)) +
    geom_tile(colour = "grey85") +
    # Outline the cells that are nonzero in the true loading matrix
    geom_tile(data = subset(df, true_nonzero %in% TRUE),
              colour = "black", linewidth = 0.5, fill = NA) +
    scale_fill_gradient(low = "white", high = "firebrick", trans = "sqrt",
                        limits = c(0, mse_lim), breaks = mse_breaks,
                        na.value = "grey92",
                        guide = guide_colourbar(barwidth = 10, barheight = 0.6)) +
    ggtitle(title) +
    theme_minimal() +
    theme(axis.title = element_blank(),
          panel.grid = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = 10))
}

titles = ifelse(names(res) %in% names(LABELS), LABELS[names(res)], names(res))
panels = Map(make_panel, res, titles)
fig = wrap_plots(panels, nrow = 1) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(args[2], fig, width = 2.6 * length(panels) + 1, height = 4.5)
cat("Figure saved to", args[2], "\n")
cat("Shared colour scale: 0 to", signif(mse_max, 3), "(sqrt spacing)\n")
