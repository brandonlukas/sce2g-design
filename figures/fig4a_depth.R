# Fig 4a — scE2G-Multiome on the CRISPR benchmark across ATAC fragment counts: AUPRC (left), precision
# (centre) and recall (right) at the model threshold, one point per subsampled dataset (varying UMI and
# cell counts). Three ggplots joined with patchwork.
#
# Everything below is measured from the published vector panel:
#   panels : 39.5 x 39.5 mm each, 0.25 pt axes, 2.56 pt ticks; x 5-9 by 1 (limits 4.75-8.84 with the
#            default 5 % expansion), y 0-0.7 (recall 0-0.75) with breaks every 0.2;
#   points : #6B2063 discs, 4.3 pt diameter, 70 % opacity, no outline — the paper's 196 points per panel
#            are digitised in figures/data/fig4a_values.csv;
#   lines  : 0.5 pt, 2 on / 2 off (lty "33" = 2.25 pt, the closest R allows: dash units are 0.75 pt for lwd < 1); black vertical at 2 million fragments; horizontal references for ABC
#            (#008F99) and distance to TSS (#6A768B) on the full-depth data.
#
#   Rscript figures/fig4a_depth.R   # -> figures/fig4a_depth.{pdf,png}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
source("R/theme_sce2g.R")

d <- read.csv("figures/data/fig4a_values.csv")
ref <- tibble(metric = c("auprc", "precision", "recall"),
              title  = c("AUPRC", "Precision at model threshold", "Recall at model threshold"),
              abc    = c(0.482, 0.419, 0.700), dist = c(0.436, 0.315, 0.704), top = c(0.7, 0.7, 0.75))
panel <- function(m) {
  r <- ref[ref$metric == m, ]
  ggplot(d[d$metric == m, ], aes(log10_fragments, value)) +
    geom_hline(yintercept = r$dist, colour = "#6A768B", linewidth = pt_mm(.5), linetype = "33") +
    geom_hline(yintercept = r$abc,  colour = "#008F99", linewidth = pt_mm(.5), linetype = "33") +
    geom_vline(xintercept = log10(2e6), colour = "black", linewidth = pt_mm(.5), linetype = "33") +
    geom_point(colour = "#6B2063", alpha = 0.7, size = 4.3 / .pt / 0.75, stroke = 0, shape = 16) +   # pch 16 draws 0.75 x fontsize
    scale_x_continuous(breaks = 5:9, limits = c(4.75, 8.84)) +
    scale_y_continuous(breaks = seq(0, 0.6, 0.2), labels = c("0", "0.2", "0.4", "0.6"), limits = c(0, r$top)) +
    labs(x = expression(log[10] * "(number of unique ATAC fragments)"), y = r$title) +
    theme_sce2g()
}
fig <- (panel("auprc") | panel("precision") | panel("recall")) +
  plot_layout(widths = unit(rep(39.5, 3), "mm"), heights = unit(39.5, "mm"))
save_panel(fig, "figures/fig4a_depth.pdf", width_mm = 152, height_mm = 52)
save_panel(fig, "figures/fig4a_depth.png", width_mm = 152, height_mm = 52)
cat("wrote figures/fig4a_depth.{pdf,png}\n")
