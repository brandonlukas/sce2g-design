# eQTL enrichment by prediction cell type / eQTL biosample category, in the style of scE2G Fig. 5b:
# violins filled by category, black ggdist-style point interval (median dot, 66 % / 95 % lines), pairwise
# Wilcoxon brackets with stars below and the exact P above, category labels + italic n under the axis.
#   Rscript eqtl_violin.R   # -> eqtl_violin.{svg,pdf,png}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
source("R/theme_sce2g.R")
S <- sce2g$type$size_pt
txt <- function(size = S$tick_label) size / .pt

d <- read.csv("eqtl_enrichment.csv") %>%
  mutate(category = factor(category, c("unmatched", "same_lineage", "same_cell_type")))
cols <- c(unmatched = "#8B98AC", same_lineage = "#9A3983", same_cell_type = "#260032")
n <- table(d$category)

# ---- stats: unpaired two-sided Wilcoxon rank-sum, one bracket per pair ----
wp <- function(a, b) wilcox.test(d$enrichment[d$category == a], d$enrichment[d$category == b])$p.value
fmt_p <- function(p) {
  if (p >= 0.001) return(sprintf("italic(P)*' = %s'", signif(p, 3)))
  e <- floor(log10(p)); sprintf("italic(P)*' = %.2f × '*10^'–%d'", p / 10^e, -e)
}
stars <- function(p) if (p < 1e-3) "***" else if (p < 1e-2) "**" else if (p < 0.05) "*" else "NS"
tests <- tibble(x1 = c(1, 2, 1), x2 = c(2, 3, 3)) %>%
  rowwise() %>% mutate(p = wp(levels(d$category)[x1], levels(d$category)[x2])) %>% ungroup()

# ---- geometry: y spans ylim over the 50.1 mm panel, so u = data units per printed pt ----
ylim <- c(-0.775, 18.4); u <- diff(ylim) / (50.1 / 25.4 * 72)
y1 <- max(d$enrichment) + 8 * u                          # first bracket row clears the tallest violin
y2 <- y1 + 14 * u; y3 <- y2 + 14 * u                      # 14 pt rows clear the exponent of the P below
tip <- 1.1 * u; line <- 5.4 * u
lab <- function(x, y, label, size = S$tick_label, ...) annotate("text", x = x, y = y, label = label, size = txt(size), family = sce2g_font, vjust = 0, ...)
rows <- -(1.913 + 0.743 * 0:2)                            # label rows under the axis, 5.9 pt pitch

p <- ggplot(d, aes(category, enrichment)) +
  geom_violin(aes(fill = category), colour = NA, width = 0.9) +
  scale_fill_manual(values = cols, guide = "none") +
  stat_pointinterval_sce2g() +
  # brackets: one row each (two exponent-length P values do not fit side by side on a 34.7 mm panel)
  sig_bracket(1, 2, y1, stars(tests$p[1]), p = fmt_p(tests$p[1]), tip = tip, line = line) +
  sig_bracket(2, 3, y2, stars(tests$p[2]), p = fmt_p(tests$p[2]), tip = tip, line = line) +
  sig_bracket(1, 3, y3, stars(tests$p[3]), p = fmt_p(tests$p[3]), tip = tip, line = line) +
  # category labels: label lines first, italic n on the very next row (Fig. 5b)
  lab(1, rows[1], "Unmatched") + lab(1, rows[2], sprintf("italic(n)*' = %s'", format(n[[1]], big.mark = ",")), parse = TRUE) +
  lab(2, rows[1], "Same") + lab(2, rows[2], "lineage") + lab(2, rows[3], sprintf("italic(n)*' = %s'", n[[2]]), parse = TRUE) +
  lab(3, rows[1], "Same") + lab(3, rows[2], "cell type") + lab(3, rows[3], sprintf("italic(n)*' = %s'", n[[3]]), parse = TRUE) +
  lab(2, -4.68, "Prediction cell type or eQTL biosample category", S$axis_title) +
  scale_y_continuous(breaks = seq(0, 15, 5)) +
  coord_cartesian(xlim = c(0.4, 3.6), ylim = ylim, expand = FALSE, clip = "off") +
  labs(x = NULL, y = "Enrichment\neQTLs versus common variants") +
  theme_sce2g() + theme(axis.text.x = element_blank(), plot.margin = margin(6, 4, 38, 4))

fig <- p + plot_layout(widths = unit(34.7, "mm"), heights = unit(50.1, "mm"))
save_panel(fig, "eqtl_violin", width_mm = 62, height_mm = 70)
print(tests)
cat("wrote eqtl_violin.{svg,pdf,png}\n")
