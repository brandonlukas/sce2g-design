# Fig 5a — similarity of scE2G-Multiome predictions between 46 cell types: Jaccard similarity of the
# thresholded links (upper left, red) and Pearson correlation of scores (lower right, blue), with the
# anti-diagonal in navy, a dataset side bar and two colour bars.
#
# The whole panel is drawn in the paper's PDF point coordinates (x right, y negated) so every element
# sits exactly where it was measured:
#   cells     : 3.227 pt squares from (103.7, 65.7); values digitised from the cell colours against the
#               paper's own colour-bar rasters (figures/data/fig5a_cells.csv, 0-0.6 for both metrics);
#               columns run in the reverse order of the rows, so the diagonal is the anti-diagonal (#001636);
#   side bar  : 2.9 pt wide at x = 254.8, one 3.227 pt block per row (figures/data/fig5a_sidebar.csv);
#   dataset key: four 7.8 pt squares at x = 264.7 (cell line #142D4A, BMMC #42566B, PBMC #6A768B, islets #8B98AC);
#   colour bars: 7.8 x 38.3 pt, the paper's 160-row ramps (figures/data/fig5a_ramps.csv), 0 at the bottom,
#               0.6 at the top, 1.5 pt white ticks on both edges every 0.1 (0.238 pt);
#   text      : every label as extracted (figures/data/fig5a_text.csv): row labels 2.85 pt right-aligned,
#               column labels rotated 60 deg, superscripts 1.9 pt, legend keys 4.28 pt, titles 5.22 pt
#               (the paper sets one column label, "BMMC erythroblasts", at 3.32 pt; regularised to 2.85).
#
#   Rscript fig5a_heatmap.R   # -> fig5a_heatmap.{pdf,png}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
source(Sys.getenv("SCE2G_THEME", "../../assets/theme_sce2g.R"))

sz <- 3.227; x0 <- 103.7; y0 <- 65.7; n <- 46
cells <- read.csv("data/fig5a_cells.csv")
ramps <- read.csv("data/fig5a_ramps.csv")           # row_from_top 0 = 0.6, 159 = 0
sidebar <- read.csv("data/fig5a_sidebar.csv")
text <- read.csv("data/fig5a_text.csv")
ramp_hex <- function(metric, v) {                            # value -> the paper's sampled ramp colour
  idx <- pmin(pmax(round((1 - v / 0.6) * 159), 0), 159)
  ramps$hex[match(paste(ifelse(metric == "jaccard", "red", "blue"), idx), paste(ramps$bar, ramps$row_from_top))]
}
cells <- cells %>% mutate(x = x0 + (col + 0.5) * sz, y = -(y0 + (row + 0.5) * sz),
                          fill = ifelse(metric == "diag", "#001636", ramp_hex(metric, value)))
sidebar <- sidebar %>% mutate(x = 254.81 + 1.45, y = -(y0 + (row + 0.5) * sz))
keys <- tibble(fill = c("#142D4A", "#42566B", "#6A768B", "#8B98AC"),        # Dataset legend: 7.8 pt squares
               x = 264.74 + 3.91, y = -(c(72.86, 80.68, 88.49, 96.32) + 3.91))
bar <- function(name, top) {                                 # 160 thin tiles + white edge ticks
  r <- ramps[ramps$bar == name, ]; h <- 38.25 / 160
  ticks <- top + 38.25 - seq(0, 0.5, 0.1) / 0.6 * 38.25
  list(geom_tile(data = r, aes(x = 264.84 + 3.91, y = -(top + (row_from_top + 0.5) * h), fill = hex), width = 7.82, height = h),
       annotate("segment", x = 264.92, xend = 266.44, y = -ticks, yend = -ticks, colour = "white", linewidth = pt_mm(.238)),
       annotate("segment", x = 271.02, xend = 272.55, y = -ticks, yend = -ticks, colour = "white", linewidth = pt_mm(.238)))
}
p <- ggplot() +
  geom_tile(data = cells, aes(x, y, fill = fill, colour = fill), width = sz, height = sz, linewidth = pt_mm(.2)) +   # hairline in the cell colour hides viewer seams
  geom_tile(data = sidebar, aes(x, y, fill = fill), width = 2.9, height = sz) +
  geom_tile(data = keys, aes(x, y, fill = fill), width = 7.81, height = 7.81) +
  bar("red", 120.13) + bar("blue", 174.74) +
  geom_text(data = text, aes(x, -y, label = text, size = size / .pt, angle = angle), hjust = 0, vjust = 0, family = sce2g_font) +
  scale_fill_identity() + scale_colour_identity() + scale_size_identity() +
  coord_cartesian(xlim = c(58, 292), ylim = c(-262, -62), expand = FALSE, clip = "off") +
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

pt2mm <- 25.4 / 72
save_panel(p, "fig5a_heatmap.pdf", width_mm = 234 * pt2mm, height_mm = 200 * pt2mm)
save_panel(p, "fig5a_heatmap.png", width_mm = 234 * pt2mm, height_mm = 200 * pt2mm)
cat("wrote fig5a_heatmap.{pdf,png}\n")
