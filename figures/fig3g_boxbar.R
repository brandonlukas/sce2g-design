# Fig 3g — TFRC in K562: (left) expression rank by promoter ATAC state, (right) quantile of the
# Kendall correlation at single-cell vs meta-cell resolution. Two ggplots joined with patchwork.
#
# Everything below is measured from the published vector panel:
#   left  : boxes 0.8 of the slot, fill #D3D3D3, outline 0.64 pt, median 1 pt, whiskers dashed
#           (1.9 pt on / 3.2 pt off with ROUND caps, so each dash reads ~2.5 / 2.5; lty "34" at 0.75 pt units),
#           staples half the box width; panel 17.4 x 35.6 mm; y 0-8,000 with
#           a 5 % top margin; P value (5.5 pt) right-aligned at the panel's top right; two label rows
#           under the axis ("Closed / Open" then counts, 5.5 pt) with 6.5 pt row headers
#           "ATAC state" / "Cell number" right-aligned at the y axis;
#   right : panel 24.1 x 34.5 mm, half a slot of space either side of the bars; bars 0.66 of the slot, #641E61 (single-cell) and #00783D (meta-cell) with a 0.25 pt black
#           outline, and a white outlined "remainder" up to exactly 100 %; y 90-100 % with NO
#           expansion; count labels 4.5 pt (above the bar when there is no room inside); 6.5 pt
#           category labels; no x ticks.
#
#   Rscript figures/fig3g_boxbar.R   # -> figures/fig3g_boxbar.{pdf,png}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
source("R/theme_sce2g.R")
S <- sce2g$type$size_pt; N <- sce2g$color$neutral; MV <- sce2g$color$model_variant
txt <- function(size = S$tick_label) size / .pt

# ---- left: box plot from digitised quartiles ------------------------------------------------------
bx <- tibble(state = factor(c("Closed", "Open"), c("Closed", "Open")), n = c(6225, 1596),
             ymin = c(113, 113), lower = c(1797, 2643), middle = c(3742, 4492), upper = c(5800, 6123), ymax = c(7815, 7815))
y_per_pt <- 2000 / 22.93                       # data units per pt, from the tick spacing
row1 <- -(733.1 - 722.3) * y_per_pt; row2 <- -(739.6 - 722.3) * y_per_pt   # label rows below the axis
p_box <- ggplot(bx, aes(x = state)) +
  geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax), stat = "identity",
               fill = N$box_fill, colour = "black", width = 0.8, linewidth = pt_mm(.639), median.linewidth = pt_mm(1),
               whisker.linetype = "34", staplewidth = 0.5, lineend = "round") +   # dash units are 0.75 pt (R clamps lwd < 1); round caps add one line width per dash as in the paper
  annotate("text", x = 2.59, y = 8130, label = "italic(P) == 8.88 %*% 10^-17", parse = TRUE, hjust = 1, vjust = 0,
           size = txt(), family = sce2g_font) +
  annotate("text", x = 1:2, y = row1, label = levels(bx$state), size = txt(), family = sce2g_font) +
  annotate("text", x = 1:2, y = row2, label = scales::comma(bx$n), size = txt(), family = sce2g_font) +
  annotate("text", x = 0.41, y = c(row1, row2), label = c("ATAC state", "Cell number"), hjust = 1,
           size = txt(S$axis_title), family = sce2g_font) +
  scale_y_continuous(breaks = seq(0, 8000, 2000), labels = scales::comma) +
  coord_cartesian(xlim = c(0.41, 2.59), ylim = c(-200, 8600), expand = FALSE, clip = "off") +
  labs(x = NULL, y = expression("Rank of " * italic(TFRC) * " gene expression")) +
  theme_sce2g() + theme(axis.text.x = element_blank(), plot.margin = margin(4, 4, 18, 8))

# ---- right: quantile bars with white remainders -------------------------------------------------
qb <- tibble(res = factor(c("Single-cell", "Meta-cell"), c("Single-cell", "Meta-cell")), k = c(9225, 8560), n = 9267,
             q = 100 * k / n, fill = c("#641E61", MV$multiome_meta_cell_bar), i = 1:2)
p_bar <- ggplot(qb, aes(x = res)) +
  geom_rect(aes(xmin = i - 0.33, xmax = i + 0.33, ymin = 90, ymax = q, fill = fill), colour = "black", linewidth = pt_mm(.25)) +   # from the 90 % floor
  geom_rect(aes(xmin = i - 0.33, xmax = i + 0.33, ymin = q, ymax = 100), fill = "white", colour = "black", linewidth = pt_mm(.25)) +
  annotate("text", x = 1, y = 100.25, label = "9,225 of 9,267", vjust = 0, size = txt(S$annotation), family = sce2g_font) +
  annotate("text", x = 2, y = qb$q[2] + 0.55, label = "8,560 of\n9,267", vjust = 0, lineheight = 0.85,
           size = txt(S$annotation), family = sce2g_font) +
  scale_fill_identity() +
  scale_y_continuous(breaks = seq(90, 100, 2), labels = function(v) paste0(v, "%")) +
  coord_cartesian(xlim = c(0.5, 2.5), ylim = c(90, 100), expand = FALSE, clip = "off") +   # x gap measured from the paper
  labs(x = NULL, y = "Quantile of Kendall correlation") +
  theme_sce2g() + theme(axis.ticks.x = element_blank(), axis.text.x = element_text(size = S$axis_title, margin = margin(t = 5)),
                        axis.title.y = element_text(margin = margin(r = 5)), plot.margin = margin(8, 4, 18, 4))

fig <- (p_box | p_bar) + plot_layout(widths = unit(c(17.4, 24.1), "mm"), heights = unit(35, "mm"))   # panel sizes measured from the paper
save_panel(fig, "figures/fig3g_boxbar.pdf", width_mm = 70, height_mm = 48)
save_panel(fig, "figures/fig3g_boxbar.png", width_mm = 70, height_mm = 48)
cat("wrote figures/fig3g_boxbar.{pdf,png}\n")
