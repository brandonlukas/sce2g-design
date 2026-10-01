# Fig 4 legend — the shared horizontal key under panels a/b: "Model", then the ABC (teal dashed) and
# distance-to-TSS (grey dashed) reference lines and the scE2G-Multiome point.
#
# Drawn in the paper's own PDF point coordinates (x right, y down) so every offset is the measured one:
# text baseline y = 356.47, keys centred 1.7 pt above it; 6.5 pt labels with 4.5 pt superscripts raised
# 3 pt; 11.9 pt dashed keys (0.5 pt, 2 on / 2 off -> lty "33", see fig4a_depth.R); 4.3 pt point at 70 %.
# The canvas is sized so 1 data unit = 1 pt (theme_void, no margins).
#
#   Rscript fig4_legend.R   # -> fig4_legend.{pdf,png}
suppressPackageStartupMessages(library(ggplot2))
source(Sys.getenv("SCE2G_THEME", "../../assets/theme_sce2g.R"))
S <- sce2g$type$size_pt; txt <- function(size) size / .pt
lab <- function(x, text, size = S$axis_title, y = 356.47) annotate("text", x = x, y = -y, label = text, hjust = 0, vjust = 0, size = txt(size), family = sce2g_font)
key <- function(x, colour) annotate("segment", x = x, xend = x + 11.92, y = -354.74, yend = -354.74, colour = colour, linewidth = pt_mm(.5), linetype = "33")
xlim <- c(170, 412); ylim <- c(-362, -348)   # y negated: PDF y grows downward (coord_cartesian sorts its limits)

p <- ggplot() +
  lab(174.9, "Model") +
  key(202.31, "#008F99") + lab(214.45, "ABC") + lab(227.62, "A=scATAC, C=power law", S$annotation, 353.47) +
  key(291.46, "#6A768B") + lab(303.99, "Distance to TSS") +
  annotate("point", x = 367.89, y = -354.74, colour = "#6B2063", alpha = 0.7, size = 4.3 / .pt / 0.75, stroke = 0, shape = 16) +
  lab(373.06, "scE2G") + lab(392.15, "Multiome", S$annotation, 353.47) +
  coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE, clip = "off") +
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

pt2mm <- 25.4 / 72
save_panel(p, "fig4_legend.pdf", width_mm = diff(xlim) * pt2mm, height_mm = diff(ylim) * pt2mm)
save_panel(p, "fig4_legend.png", width_mm = diff(xlim) * pt2mm, height_mm = diff(ylim) * pt2mm)
cat("wrote fig4_legend.{pdf,png}\n")
