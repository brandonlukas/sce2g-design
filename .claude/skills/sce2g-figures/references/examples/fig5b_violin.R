# Fig 5b — enrichment of eQTLs in scE2G-Multiome predictions when the prediction cell type is unmatched,
# in the same lineage or the same cell type as the eQTL biosample. Violins with a ggdist-style point
# interval and stacked comparison brackets.
#
# Measured from the published vector panel (Fig. 5 is printed at 95 %; sizes here are the 100 % design
# values, i.e. measured / 0.95):
#   panel   : 34.7 x 50.1 mm; y 0-15.5 with the default 5 % expansion, breaks every 5; x slots 1-3 with the
#             default 0.6 discrete padding, ticks but no tick labels (labels are drawn as three text rows);
#   violins : fills #8B98AC / #9A3983 / #260032, no outline; outlines digitised to figures/data/fig5b_violins.csv;
#   interval: black; 95 % line 0.77 pt, 66 % line 1.79 pt, median dot 3.05 pt (values digitised below);
#   brackets: 0.25 pt black with 1.1 pt down-ticks; stars 8.5 pt below, P values 4.5 pt above (plotmath superscripts for the exponents).
#
#   Rscript fig5b_violin.R   # -> fig5b_violin.{pdf,png}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
source(Sys.getenv("SCE2G_THEME", "../../assets/theme_sce2g.R"))
S <- sce2g$type$size_pt
txt <- function(size = S$tick_label) size / .pt
u <- 1 / 7.911                                            # data units per printed pt (for label offsets)

violins <- read.csv("data/fig5b_violins.csv")
cols <- c("#8B98AC", "#9A3983", "#260032")                # unmatched, same lineage, same cell type
pi <- tibble(cat = 1:3, lo95 = c(3.05, 4.34, 4.50), hi95 = c(7.39, 11.35, 11.99),
             lo66 = c(4.01, 5.69, 5.66), hi66 = c(5.85, 8.47, 9.14), med = c(4.93, 6.76, 7.60))
br <- tibble(x = c(1.064, 1.064, 2.064), xend = c(1.936, 2.936, 2.936), y = c(14.50, 15.71, 14.50))
lab <- function(x, y, label, size = S$tick_label, ...) annotate("text", x = x, y = y, label = label, size = txt(size), family = sce2g_font, vjust = 0, ...)

p <- ggplot() +
  geom_polygon(data = violins, aes(x, y, group = cat, fill = factor(cat))) +
  scale_fill_manual(values = cols, guide = "none") +
  geom_linerange(data = pi, aes(x = cat, ymin = lo95, ymax = hi95), linewidth = pt_mm(0.77)) +
  geom_linerange(data = pi, aes(x = cat, ymin = lo66, ymax = hi66), linewidth = pt_mm(1.79)) +
  geom_point(data = pi, aes(cat, med), size = 3.05 / .pt / 0.75, stroke = 0, shape = 16) +
  geom_segment(data = br, aes(x = x, xend = xend, y = y, yend = y), linewidth = pt_mm(sce2g$stroke_pt$bracket)) +
  geom_segment(data = br, aes(x = x, xend = x, y = y, yend = y - 0.136), linewidth = pt_mm(sce2g$stroke_pt$bracket)) +
  geom_segment(data = br, aes(x = xend, xend = xend, y = y, yend = y - 0.136), linewidth = pt_mm(sce2g$stroke_pt$bracket)) +
  # stars below each bracket, P values above (positions measured; hjust = 0 at the run start as printed)
  lab(1.416, 14.19, "***", 8.5, hjust = 0) + lab(1.916, 15.40, "***", 8.5, hjust = 0) + lab(2.472, 14.19, "*", 8.5, hjust = 0) +
  lab(1.504, 16.27, "italic(P)*' = 2.45 × '*10^'–23'", S$annotation, hjust = 0, parse = TRUE) +
  lab(1.004, 15.11, "italic(P)*' = 4.39 × '*10^'–41'", S$annotation, hjust = 0, parse = TRUE) +
  lab(2.143, 15.12, "italic(P)*' = 0.0324'", S$annotation, hjust = 0, parse = TRUE) +
  # category labels: up to three rows under the axis, 5.9 pt apart
  lab(1, -1.913, "Unmatched") + lab(1, -2.656, "italic(n)*' = 2,259'", parse = TRUE) +
  lab(2, -1.913, "Same") + lab(2, -2.656, "lineage") + lab(2, -3.383, "italic(n)*' = 204'", parse = TRUE) +
  lab(3, -1.913, "Same") + lab(3, -2.656, "cell type") + lab(3, -3.383, "italic(n)*' = 97'", parse = TRUE) +
  lab(2, -4.68, "Prediction cell type or eQTL biosample category", S$axis_title) +
  scale_x_continuous(breaks = 1:3) + scale_y_continuous(breaks = seq(0, 15, 5)) +
  coord_cartesian(xlim = c(0.4, 3.6), ylim = c(-0.775, 16.27), expand = FALSE, clip = "off") +
  labs(x = NULL, y = "Enrichment\neQTLs versus common variants") +
  theme_sce2g() + theme(axis.text.x = element_blank(), plot.margin = margin(6, 4, 38, 4))

fig <- p + plot_layout(widths = unit(34.7, "mm"), heights = unit(50.1, "mm"))   # panel size measured from the paper
save_panel(fig, "fig5b_violin.pdf", width_mm = 62, height_mm = 70)
save_panel(fig, "fig5b_violin.png", width_mm = 62, height_mm = 70)
cat("wrote fig5b_violin.{pdf,png}\n")
