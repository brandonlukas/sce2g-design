# Fig. 5b-style violin: eQTL enrichment by biosample match, point-interval, Wilcoxon brackets.
setwd(dirname(normalizePath(sys.frames()[[1]]$ofile %||% "eqtl_violin.R", mustWork = FALSE)))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("R/theme_sce2g.R")
txt <- function(pt = sce2g$type$size_pt$tick_label) pt / .pt

lv   <- c("unmatched", "same_lineage", "same_cell_type")
cols <- setNames(c("#8B98AC", "#9A3983", "#260032"), lv)
d <- read.csv("eqtl_enrichment.csv"); d$category <- factor(d$category, lv)
n  <- table(d$category)

# two-sided Mann-Whitney (Wilcoxon rank-sum) per pair
pw <- function(a, b) wilcox.test(d$enrichment[d$category == a], d$enrichment[d$category == b])$p.value
plab <- function(p) { e <- floor(log10(p)); sprintf("italic(P)*' = %.2f × '*10^'–%d'", p / 10^e, -e) }
stars <- function(p) if (p < 1e-3) "***" else if (p < 1e-2) "**" else if (p < .05) "*" else "NS"
p12 <- pw(lv[1], lv[2]); p23 <- pw(lv[2], lv[3]); p13 <- pw(lv[1], lv[3])

# y in data units: panel 50.1 mm = 142 pt over 17.2 units -> 8.26 pt/unit; text line ~0.6 units, tick 1.1 pt ~0.13
ln <- 0.6; tp <- 0.13; y1 <- 12.3; y2 <- y1 + 15 / 8.26   # stacked brackets 15 pt apart (lower P exponent clears the upper bar)
lab_rows <- -0.9 - c(0.75, 1.55, 2.35)                       # three text rows under the axis, 6.6 pt pitch
cat_lab <- list(c("Unmatched", ""), c("Same", "lineage"), c("Same", "cell type"))

p <- ggplot(d, aes(category, enrichment)) +
  geom_violin(aes(fill = category), colour = NA, width = 0.9) +
  scale_fill_manual(values = cols, guide = "none") +
  stat_pointinterval_sce2g() +
  # shared-row brackets: each P set outward from the shared middle category so they do not collide
  sig_bracket(1, 2, y1, stars(p12), tip = tp, line = ln) +
  sig_bracket(2, 3, y1, stars(p23), tip = tp, line = ln) +
  annotate("text", x = c(1.93, 2.07), y = y1 + 1.5 * ln, label = c(plab(p12), plab(p23)), hjust = c(1, 0), vjust = 0,
           parse = TRUE, size = txt(sce2g$type$size_pt$annotation), family = sce2g_font) +
  sig_bracket(1, 3, y2, stars(p13), p = plab(p13), tip = tp, line = ln) +
  annotate("text", x = 1:3, y = lab_rows[1], label = sapply(cat_lab, `[`, 1), size = txt(), family = sce2g_font, vjust = 0) +
  annotate("text", x = 1:3, y = lab_rows[2], label = sapply(cat_lab, `[`, 2), size = txt(), family = sce2g_font, vjust = 0) +
  annotate("text", x = 1:3, y = lab_rows[3], label = sprintf("italic(n)*' = %s'", prettyNum(as.integer(n), big.mark = ",")),
           parse = TRUE, size = txt(), family = sce2g_font, vjust = 0) +
  scale_y_continuous(breaks = seq(0, 15, 5)) +
  coord_cartesian(xlim = c(0.4, 3.6), ylim = c(-0.9, 16.3), expand = FALSE, clip = "off") +
  labs(x = NULL, y = "Enrichment\neQTLs versus common variants") +
  theme_sce2g() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(6, 4, 30, 4))

fig <- p + plot_layout(widths = unit(34.7, "mm"), heights = unit(50.1, "mm"))
save_panel(fig, "eqtl_violin.pdf", width_mm = 55, height_mm = 75)
save_panel(fig, "eqtl_violin.png", width_mm = 55, height_mm = 75)
cat(sprintf("n: %s\nP unmatched-lineage %.3g, lineage-cell type %.3g, unmatched-cell type %.3g\n",
            paste(names(n), n, collapse = ", "), p12, p23, p13))
