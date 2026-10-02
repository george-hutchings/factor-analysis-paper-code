### True loading matrix for the high-dimensional scenario (R&G 2016 dimensions).
###   Rscript plot_true_loading.R <out.pdf>
###
### The low-dimensional scenarios draw this with pheatmap (see
### simulations/create_true_loading_matrix.R), which lays down one
### grid rectangle per cell and labels every row. At G = 1956 that is neither
### legible nor a sensible PDF, so this uses geom_raster with the rows left
### unlabelled and tick marks every 500 variables, which is what makes the block
### structure and the 136-variable overlaps readable.

library(ggplot2)
library(reshape2)

args = commandArgs(trailingOnly = TRUE)
out = if (length(args) >= 1) args[1] else "true_loading_highdim.pdf"

G = 1956
dd = 5
blk = 500
offset = 135
lambda = matrix(0, G, dd)
end = 1
for (i in 1:dd) {
  start = end - (i != 1) * offset
  end = start + blk - 1
  lambda[start:end, i] = 1
}

df = melt(lambda, varnames = c("Variable", "Factor"), value.name = "Loading")

fig = ggplot(df, aes(x = Factor, y = Variable, fill = Loading)) +
  geom_raster() +
  scale_fill_gradient(low = "white", high = "#2166AC",
                      breaks = c(0, 1), limits = c(0, 1),
                      guide = guide_colourbar(barwidth = 8, barheight = 0.6)) +
  # Rows count downwards, so variable 1 is at the top, as in the other heatmaps.
  scale_y_reverse(expand = c(0, 0), breaks = c(1, seq(500, G, by = 500), G)) +
  scale_x_continuous(expand = c(0, 0), breaks = seq_len(dd), position = "top") +
  labs(x = "Factor", y = "Variable") +
  theme_minimal() +
  theme(panel.grid = element_blank(),
        legend.position = "bottom",
        axis.text = element_text(size = 8))

ggsave(out, fig, width = 3.4, height = 4.2)
cat("Figure saved to", out, "\n")
cat(sprintf("G = %d, K_true = %d, %d nonzeros per factor, adjacent overlap %d\n",
            G, dd, blk, offset + 1))
