# Fig 3a — forward feature selection for scE2G-Multiome: delta AUPRC from adding each feature.
#
# Values digitised from the published vector panel (figures/data/fig3a_values.csv). Anatomy:
#   * rows in selection order, "None" at the top; y title "Feature added"; features retained in the
#     final model are underlined (plotmath underline() in the axis labels);
#   * horizontal bars (#445369, 90 % of the 9.85 pt row pitch) = delta AUPRC, with a 95 % CI error
#     bar (0.32 pt, 2 pt caps); 2.3 pt black dot = cumulative AUPRC after adding the feature;
#   * grey dashed vertical at the final model's AUPRC (0.675);
#   * 7.5 pt significance stars rotated 90 deg beside the CI, bold, black for positive and red for negative
#     deltas (#E82829); exact bootstrap P in a right-hand column headed P_bootstrap;
#   * x axis 0-1 in quarters ("0", "0.25", ...), zero offset ~9 pt from the y axis; "Delta AUPRC".
#
#   Rscript fig3a_features.R   # -> fig3a_features.{pdf,png}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(readr); library(patchwork) })
source(Sys.getenv("SCE2G_THEME", "../../assets/theme_sce2g.R"))
S <- sce2g$type$size_pt; N <- sce2g$color$neutral
txt <- function(size = S$tick_label) size / .pt
BAR <- "#445369"; FINAL <- 0.675

d <- read_csv("data/fig3a_values.csv", col_types = cols(feature = col_character(), p = col_character(), stars = col_character(), star_side = col_character())) %>%
  mutate(row = n() + 1 - row_number(),                                     # top row = highest y
         negative = !is.na(delta) & delta < 0,
         star_x = ifelse(star_side == "L", pmin(lo, 0) - 0.045, pmax(hi, 0) + 0.055),   # ~0.05 clear of the CI, as in the paper
         star_col = ifelse(negative, N$significance_negative, "black"))
# label widths in mm, measured on the same renderer used for output (plotmath's underline() cannot
# be made as fine as the paper's 0.36 pt rule, so labels and rules are drawn by hand)
label_width_mm <- function(x, size) {
  f <- tempfile(fileext = ".png"); ragg::agg_png(f, width = 10, height = 10, units = "mm", res = 300); on.exit({ dev.off(); unlink(f) })
  vapply(x, function(s) grid::convertWidth(grid::stringWidth(s), "mm", valueOnly = TRUE), 1, USE.NAMES = FALSE) *
    size / 12   # stringWidth uses the default 12 pt fontsize (gpar fontfamily is picked up by the device default)
}
d <- d %>% mutate(w_mm = label_width_mm(feature, S$tick_label))
LAB_W <- 40   # label column width, mm (1 data unit = 1 mm in that panel)
p_lab <- ggplot(d) +
  geom_text(aes(x = -1.2, y = row, label = feature), hjust = 1, size = txt(), family = sce2g_font) +
  geom_segment(data = filter(d, in_model), aes(x = -1.2 - w_mm, xend = -1.2, y = row - 0.30, yend = row - 0.30), linewidth = pt_mm(.36)) +
  annotate("text", x = -LAB_W + 2.5, y = (nrow(d) + 1) / 2, label = "Feature added", angle = 90, size = txt(S$axis_title), family = sce2g_font) +
  coord_cartesian(xlim = c(-LAB_W, 0), ylim = c(0.4, nrow(d) + 0.6), expand = FALSE, clip = "off") +
  theme_void() + theme(plot.margin = margin(4, 0, 4, 0))

p <- ggplot(d, aes(y = row)) +
  geom_vline(xintercept = FINAL, colour = "#C6C6C6", linetype = "33", linewidth = pt_mm(.5)) +   # paper 2/2 pt; R dash unit is 0.75 pt
  geom_col(aes(x = delta), fill = BAR, width = 0.9, orientation = "y", na.rm = TRUE) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0.2, linewidth = pt_mm(.32), orientation = "y", na.rm = TRUE) +
  geom_point(aes(x = cum_auprc), size = 2.28 / .pt / 0.75, shape = 16, stroke = 0, colour = "black") +   # 2.3 pt disc, measured
  geom_text(data = filter(d, !is.na(stars), stars != ""), aes(x = star_x, label = stars, colour = star_col),
            angle = 90, fontface = "bold", size = 7.5 / .pt, family = sce2g_font) +
  geom_text(data = filter(d, !is.na(p)), aes(x = 0.78, label = p), hjust = 0, size = txt(S$axis_title), family = sce2g_font) +
  annotate("text", x = 0.77, y = max(d$row), label = "italic(P)[bootstrap]", parse = TRUE, hjust = 0, size = txt(S$axis_title), family = sce2g_font) +
  scale_colour_identity() +
  scale_x_continuous(breaks = seq(0, 1, .25), labels = c("0", "0.25", "0.50", "0.75", "1.00")) +
  scale_y_continuous(breaks = d$row, labels = NULL) +
  coord_cartesian(xlim = c(-0.11, 1.05), ylim = c(0.4, nrow(d) + 0.6), expand = FALSE, clip = "off") +   # 0.6 row above/below, measured
  labs(x = "Delta AUPRC", y = NULL) +
  theme_sce2g() +
  theme(axis.ticks.y = element_line(), plot.margin = margin(4, 6, 4, 0))
fig <- (p_lab | p) + plot_layout(widths = unit(c(LAB_W, 1), c("mm", "null")))

save_panel(fig, "fig3a_features.pdf", width_mm = 78, height_mm = 53)
save_panel(fig, "fig3a_features.png", width_mm = 78, height_mm = 53)
cat("wrote fig3a_features.{pdf,png}\n")
