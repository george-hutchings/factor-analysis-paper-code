### Per-element loading MSE at each sample size (supplementary Fig S4).
###
### Usage: Rscript plot_elementwise.R <base_save_dir> [out_dir]
###
### Writes
###   elementwise_mse_rank.pdf   all six N in a row, on one log10 colour scale
###
### White to firebrick, true support outlined in black. One shared log10 scale across all six
### panels, so the figure reads "how does the error shrink with N?". Log rather than a linear
### or sqrt scale because the mean MSE falls by a factor of several hundred from N = 50 to
### N = 5000: a difference-based scale spends everything on N = 50 and leaves the last three
### panels blank. The floor is the smallest nonzero cell present, capped at
### five decades below the maximum so that one near-zero cell cannot stretch the ramp until
### everything else bunches at the top.
###
### Built from lambda_hat in the colSums(w) column order, matching the table.

library(ggplot2)
library(reshape2)
library(patchwork)
library(scales)

args = commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript plot_elementwise.R <base_save_dir> [out_dir]")
base_dir = args[1]
out_dir  = if (length(args) >= 2) args[2] else base_dir
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ARMS = "rank"

f = file.path(base_dir, "elementwise_mse.rds")
if (!file.exists(f)) stop("no elementwise_mse.rds in ", base_dir, "; run analyse.R first")
ew = readRDS(f)

# analyse.R keys these as "<arm> N=<N>". Parse rather than assume an ordering.
cells = lapply(names(ew), function(k) {
  p = regmatches(k, regexec("^(\\S+) N=([0-9]+)$", k))[[1]]
  if (length(p) != 3) stop("unexpected key in elementwise_mse.rds: ", k)
  list(arm = p[2], N = as.numeric(p[3]), x = ew[[k]])
})

D     = nrow(cells[[1]]$x$mse)
K_max = max(sapply(cells, function(c) ncol(c$x$mse)))

# A cell built from a handful of survivors is a statement about the survivors, so a panel where
# fewer than half the runs survived is flagged.
n_full = max(sapply(cells, function(c) c$x$n))

panel = function(cell, fill_scale, title, title_size = 9) {
  df = melt(cell$x$mse, varnames = c("Variable", "Factor"), value.name = "MSE")
  df$true_nonzero = melt(cell$x$true != 0)$value
  df$Variable = factor(df$Variable, levels = rev(seq_len(D)))
  df$Factor   = factor(df$Factor,   levels = seq_len(K_max))

  # A blank space rather than NULL on the unflagged panels. The subtitle occupies vertical
  # space, so omitting it entirely would make a flagged panel shorter than its neighbours.
  sub = if (cell$x$n < 0.5 * n_full) paste0("n = ", cell$x$n, " of ", n_full,
                                            ", survivors only") else " "

  ggplot(df, aes(x = Factor, y = Variable, fill = MSE)) +
    geom_tile(colour = "grey85") +
    # Outline the cells that are nonzero in the true loading matrix
    geom_tile(data = subset(df, true_nonzero %in% TRUE),
              colour = "black", linewidth = 0.5, fill = NA) +
    fill_scale +
    labs(title = title, subtitle = sub) +
    theme_minimal(base_size = 9) +
    theme(axis.title = element_blank(), panel.grid = element_blank(),
          legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = title_size),
          plot.subtitle = element_text(hjust = 0.5, size = title_size - 1.5,
                                       colour = "grey30"))
}

# ---- Every N in a row, on one log10 scale ----
# The scale is computed over every panel, so a colour means the same MSE in all six.
#
# log10 rather than a sqrt scale. A sqrt scale is difference-based, which is right when
# everything in the figure is within about 5x. Here the figure spans N = 50 to N = 5000 and the
# mean MSE falls by a factor of several hundred, so a difference-based scale spends its whole
# range on the first panel and leaves the last three blank white. A ratio-based scale gives a
# cell 100x below the maximum a definite position.
hi = max(sapply(cells, function(c) max(c$x$mse)))
# The floor has to sit below the smallest value actually present, or the panels at large N are
# squished to white and the figure stops showing the thing it exists to show. Take the true
# minimum, capped at five decades so that a single near-zero cell cannot stretch the ramp until
# everything else bunches at the top.
lo = max(min(unlist(lapply(cells, function(c) c$x$mse[c$x$mse > 0]))), hi / 1e5)
sc = scale_fill_gradient(
  low = "white", high = "firebrick", trans = "log10",
  limits = c(lo, hi), breaks = 10^seq(ceiling(log10(lo)), floor(log10(hi))),
  labels = scales::label_scientific(digits = 1),
  oob = scales::squish, na.value = "grey92",
  guide = guide_colourbar(barwidth = 14, barheight = 0.6))

for (a in ARMS) {
  here = Filter(function(c) c$arm == a, cells)
  here = here[order(sapply(here, function(c) c$N))]
  if (!length(here)) next

  # The titles carry only N: the likelihood belongs in the caption, not the image.
  fig = wrap_plots(Map(function(c) panel(c, sc, title = paste0("N = ", c$N), title_size = 9),
                       here), nrow = 1) +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  out = file.path(out_dir, paste0("elementwise_mse_", a, ".pdf"))
  ggsave(out, fig, width = 2.3 * length(here) + 1, height = 4.2)
  cat("Wrote", out, "\n")
}
cat("All panels on the same log10 scale,", signif(lo, 3), "to", signif(hi, 3), "(",
    sprintf("%.1f", log10(hi / lo)), "decades ).\n")
if (lo > min(unlist(lapply(cells, function(c) c$x$mse[c$x$mse > 0])))) {
  cat("Note: the five-decade cap is binding, so the smallest cells are squished to white.\n")
}
