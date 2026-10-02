### Marginal distributions of the NO.MS variables behind Scenario 3 (supplement Fig S6, fig:margins_mixed).
###
### Usage: Rscript plot_margins.R [csv] [out_dir]
###
### Writes margins_mixed.pdf, one histogram per variable, and prints the N and
### missing fraction of each for the caption.
###
### Scenario 3 generates data by pushing the latent Gaussians through the inverse
### empirical CDF of each real variable (simulate_missing_data_functions.R:26), so
### these observed margins are the F_d the simulation draws through.
###
### Needs the clinical extract, so it runs on Novartis.

library(ggplot2)

# Columns as missingness-clinical.R reads them. RELAPSE is binary and EDSS and
# NUMGDT1 are ordinal there; all eight are drawn as histograms here.
CSV  = "/data/ms/processed/mri/MS_Share/4George/03_baseline_results.csv"
COLS = c("RELAPSE", "EDSS", "T25FWM", "HPT9M", "PASAT", "NUMGDT1", "VOLT2", "NBV2")

# Panel titles follow the manuscript: Table 4's names, and the acronyms it defines in
# the real-data section, which are EDSS, SDMT and PASAT only (the walk and peg tests
# are always written out there, so they are here too). Units added: both timed tests
# are recorded in seconds and both volumes in mm^3 (SIENAX normalised brain volume
# averages about 1.6e6 mm^3 in MS, and T2 lesion volume a few thousand).
LABELS = c(RELAPSE = "Relapse status",
           EDSS    = "EDSS",
           T25FWM  = "Timed 25-Foot Walk (s)",
           HPT9M   = "9-Hole Peg Test (s)",
           PASAT   = "PASAT",
           NUMGDT1 = "Number of Gd-enhancing\nT1 lesions",
           VOLT2   = "T2 lesion volume (mm³)",
           NBV2    = "Normalised brain volume\n(mm³)")

args    = commandArgs(trailingOnly = TRUE)
csv     = if (length(args) >= 1) args[1] else CSV
out_dir = if (length(args) >= 2) args[2] else "."

Y = read.csv(csv)[COLS]

## Aim for 30 bins, but round the width up to a whole multiple of the spacing the
## variable's values actually take. Binning finer than that spacing puts empty bins
## between the legal values (EDSS moves in half steps and combs at 30 bins), and
## binning at a fractional multiple makes bins swallow alternately 2 then 1 value,
## which is what made PASAT jagged. ggplot2's own advice for discrete x is
## stat_count(); this is that, generalised so one code path covers all 8 variables.
binned = do.call(rbind, lapply(COLS, function(v) {
  x = Y[[v]][!is.na(Y[[v]])]
  u = sort(unique(x))
  grain = if (length(u) > 1) min(diff(u)) else 1
  w = grain * max(1, ceiling(diff(range(x)) / 30 / grain))
  # floor(x/w + 0.5), not round(), which is half-to-even and so breaks ties two ways
  tab = table(floor(x / w + 0.5) * w)
  data.frame(variable = v, value = as.numeric(names(tab)), n = as.numeric(tab), w = w)
}))
binned$variable = factor(binned$variable, levels = COLS)

p = ggplot(binned, aes(value, n, width = w)) +
  geom_col() +
  facet_wrap(~ variable, scales = "free", ncol = 4,
             labeller = labeller(variable = LABELS)) +
  # Keep only whole-number ticks. Every variable here is measured on an integer or
  # half-integer scale, and without this the binary relapse panel is labelled
  # -0.5, 0.0, 0.5, 1.0, 1.5, three of which it cannot take.
  scale_x_continuous(
    breaks = function(l) {
      b = scales::extended_breaks()(l)
      b[b == round(b)]
    },
    # Abbreviate the volumes. Written in full, 1800000 is both hard to read and wide
    # enough that the last brain-volume tick runs off the page.
    labels = function(b) {
      ifelse(is.na(b), "",
      ifelse(abs(b) >= 1e6, paste0(signif(b / 1e6, 3), "M"),
      ifelse(abs(b) >= 1e4, paste0(signif(b / 1e3, 3), "k"),
             format(b, trim = TRUE, scientific = FALSE))))
    }) +
  labs(x = NULL, y = "Count") +
  theme_minimal(base_size = 9) +
  theme(
    # Right margin only, so the last x tick label is not cut off by the page edge.
    plot.margin     = margin(0, 6, 0, 0),
    panel.spacing.x = unit(4, "pt"),
    panel.spacing.y = unit(14, "pt"),
    strip.text      = element_text(margin = margin(1, 0, 2, 0))
  )

out = file.path(out_dir, "margins_mixed.pdf")
ggsave(out, p, width = 7.5, height = 4.2)
cat("Wrote", out, "\n\n")

for (v in COLS) {
  x = Y[[v]][!is.na(Y[[v]])]
  cat(sprintf("%-9s N = %5d  %.4g to %.4g  median %-9.4g %.1f%% missing\n",
              v, length(x), min(x), max(x), median(x), 100 * mean(is.na(Y[[v]]))))
}
