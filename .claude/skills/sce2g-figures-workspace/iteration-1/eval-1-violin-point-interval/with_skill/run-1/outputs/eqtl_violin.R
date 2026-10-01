# eQTL enrichment by prediction cell type category, in the style of scE2G Fig. 5b:
# violins (no outline) + ggdist-style point-interval, Wilcoxon brackets, category/n label rows.
#   Rscript eqtl_violin.R   # -> eqtl_violin.{pdf,png}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
source("R/theme_sce2g.R")
S <- sce2g$type$size_pt
txt <- function(size = S$tick_label) size / .pt

d <- read.csv("eqtl_enrichment.csv") %>%
  mutate(category = factor(category, c("unmatched", "same_lineage", "same_cell_type")))
cols <- c(unmatched = "#8B98AC", same_lineage = "#9A3983", same_cell_type = "#260032")
names_rows <- list(unmatched = "Unmatched", same_lineage = c("Same", "lineage"), same_cell_type = c("Same", "cell type"))
n <- table(d$category)

# --- stats: two-sided Wilcoxon rank-sum, plotmath labels with the exponent as a superscript ---
wt <- function(a, b) wilcox.test(d$enrichment[d$category == a], d$enrichment[d$category == b])$p.value
p_lab <- function(p) {
  if (p >= 1e-3) return(sprintf("italic(P)*' = %s'", signif(p, 3)))
  e <- floor(log10(p)); sprintf("italic(P)*' = %.2f × '*10^'–%d'", p / 10^e, -e)
}
stars <- function(p) if (p < 1e-3) "***" else if (p < 1e-2) "**" else if (p < 0.05) "*" else "NS"
cmp <- tibble(x1 = c(1, 2, 1), x2 = c(2, 3, 3), row = c(1, 1, 2)) %>%
  rowwise() %>% mutate(p = wt(levels(d$category)[x1], levels(d$category)[x2])) %>% ungroup() %>%
  group_by(row) %>% mutate(hj = if (n() == 1) 0.5 else ifelse(x1 == min(x1), 1, 0)) %>% ungroup()

# --- vertical layout in printed points, converted to data units (panel is 50.1 mm = 142 pt tall) ---
top <- max(d$enrichment)
y_max_axis <- ceiling(top / 5) * 5                        # tick every 5
# data units per printed pt: the panel (142 pt) holds 5 % padding below 0, the violins, and a 28.8 pt bracket stack
u <- (top + 0.05 * y_max_axis) / (142 - 8 - 11.5 - 4.8 - 4.5)
pt <- function(x) x * u
# bracket rows: first row 8 pt above the tallest violin, second 11.5 pt higher (paper uses 9.6, which lets the superscripts touch the bar above); P text baseline 4.8 pt above each bar
y_row <- top + pt(8) + (cmp$row - 1) * pt(11.5)
y_lim <- c(-0.05 * y_max_axis, max(y_row) + pt(4.8) + pt(4.5))

lab <- function(x, y, label, size = S$tick_label, ...) annotate("text", x = x, y = y, label = label, size = txt(size), family = sce2g_font, vjust = 0, ...)
row_y <- function(i) y_lim[1] - pt(9) - (i - 1) * pt(5.9)  # label rows under the axis, 5.9 pt apart

p <- ggplot(d, aes(category, enrichment)) +
  geom_violin(aes(fill = category), colour = NA, width = 0.9, scale = "width") +
  scale_fill_manual(values = cols, guide = "none") +
  stat_pointinterval_sce2g() +
  # brackets: 0.25 pt, 1.1 pt down-ticks; stars 8.5 pt below the bar, P 4.5 pt above
  lapply(seq_len(nrow(cmp)), function(i) list(
    annotate("segment", x = cmp$x1[i] + .064, xend = cmp$x2[i] - .064, y = y_row[i], yend = y_row[i], linewidth = pt_mm(sce2g$stroke_pt$bracket)),
    annotate("segment", x = c(cmp$x1[i] + .064, cmp$x2[i] - .064), xend = c(cmp$x1[i] + .064, cmp$x2[i] - .064),
             y = y_row[i], yend = y_row[i] - pt(1.1), linewidth = pt_mm(sce2g$stroke_pt$bracket)),
    lab((cmp$x1[i] + cmp$x2[i]) / 2, y_row[i] - pt(2.5), stars(cmp$p[i]), S$stars, fontface = "bold"),
    # P values: centred when a bracket has its row to itself; brackets sharing a row align their P outward from the shared category
    lab(c(cmp$x2[i] - .064, (cmp$x1[i] + cmp$x2[i]) / 2, cmp$x1[i] + .064)[3 - cmp$hj[i] * 2], y_row[i] + pt(4.8), p_lab(cmp$p[i]), S$annotation, parse = TRUE, hjust = cmp$hj[i]))) +
  # category labels as rows under the axis, italic n on the last line
  lapply(1:3, function(i) { r <- names_rows[[i]]; c(
    lapply(seq_along(r), function(j) lab(i, row_y(j), r[j])),
    list(lab(i, row_y(length(r) + 1), sprintf("italic(n)*' = %s'", format(n[[i]], big.mark = ",")), parse = TRUE))) }) +
  lab(2, row_y(4) - pt(2), "Prediction cell type or eQTL biosample category", S$axis_title) +
  scale_x_discrete(labels = NULL) + scale_y_continuous(breaks = seq(0, y_max_axis, 5)) +
  coord_cartesian(xlim = c(0.4, 3.6), ylim = y_lim, expand = FALSE, clip = "off") +
  labs(x = NULL, y = "Enrichment\neQTLs versus common variants") +
  theme_sce2g() + theme(axis.text.x = element_blank(), plot.margin = margin(6, 4, 38, 4))

fig <- p + plot_layout(widths = unit(34.7, "mm"), heights = unit(50.1, "mm"))   # Fig. 5b panel size
save_panel(fig, "eqtl_violin.pdf", width_mm = 62, height_mm = 70)
save_panel(fig, "eqtl_violin.png", width_mm = 62, height_mm = 70)
print(cmp)
