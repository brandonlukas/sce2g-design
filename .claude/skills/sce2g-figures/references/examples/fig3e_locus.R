# Fig 3e — scE2G-Multiome and scE2G-ATAC (single-cell) predictions for TFRC in K562, the CRISPR-
# validated enhancer, normalised ATAC and the tested CRISPR elements, as a locus plot.
#
# Everything is measured from the published vector panel and digitised into figures/data/fig3e_locus.json
# (page-pt coordinates; 100 kb = 183.0 pt, so x converts to bp at 546.4 bp/pt from chr3:196,070,000):
#   rows (15.6 pt pitch): scale bar (0.5 pt, round caps, 0.4 x 2.6 pt end ticks) -> genes (11.7 pt exon blocks, 0.36 pt; #C9C9C9 intron line with one
#   strand chevron; italic 6.5 pt names) -> scE2G-Multiome arcs -> scE2G-ATAC arcs -> CRISPR-positive arc
#   -> ATAC track (9 pt, #006EAE) -> tested elements (#DC6464, 5.9 pt tall).
#   arcs : arc_paths() quadratics of constant 14.8 pt height from a 0.8 pt TSS tick to a 0.8 pt element foot,
#          colour and width scaling with score (measured hex / 0.12-0.61 pt kept as data);
#   text : 5.5 pt track names right-aligned in a 47 pt label column, superscripts via plotmath.
#
#   Rscript fig3e_locus.R   # -> fig3e_locus.{pdf,png}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(jsonlite) })
source(Sys.getenv("SCE2G_THEME", "../../assets/theme_sce2g.R"))
S <- sce2g$type$size_pt; G <- sce2g$color$genomic
txt <- function(size = S$tick_label) size / .pt
g <- fromJSON("data/fig3e_locus.json", simplifyVector = TRUE)
LAB_X <- 302.0                                                   # right edge of the label column (page pt)
ypt <- function(y) -y                                            # page y grows downward

# ---- arcs: the theme's constant-height quadratic arcs (arc_paths), one row at a time -------------
arcs <- g$arcs %>% mutate(x = x_to, xend = x_from, score = sw) %>% group_by(row) %>%
  group_modify(~ arc_paths(.x, height = .x$height[1]) %>% mutate(y = y + ypt(.y$row))) %>% ungroup()
feet <- g$feet %>% mutate(xmin = x, xmax = x + w, ymax = ypt(row - 0.786), ymin = ypt(row))   # `row` is the foot's bottom edge; the arc lands on its top
tss  <- g$arcs %>% group_by(row) %>% summarise(hex = hex[which.max(sw)], sw = max(sw), .groups = "drop")   # darkest arc marks the TSS
atac <- as.data.frame(g$atac); names(atac) <- c("x", "h")
ex   <- g$exons; ex$xmax <- ex$x + ex$w; ex$ymin <- ypt(ex$y_bottom); ex$ymax <- ypt(ex$y_bottom - ex$h)
chev <- data.frame(x = c(320.75, 462.31), y = c(466.93, 465.89), s = c(1, -1))     # "<" and ">"
els  <- g$crispr_elements
lw   <- function(pt) pt / .pt / 0.75                                                 # measured pt -> ggplot linewidth

p <- ggplot() +
  # scale bar
  annotate("segment", x = 312.14, xend = 495.15, y = ypt(456.85), yend = ypt(456.85), linewidth = pt_mm(.5), lineend = "round") +
  annotate("rect", xmin = c(311.94, 494.95), xmax = c(312.34, 495.35), ymin = ypt(458.15), ymax = ypt(455.55), fill = "black") +   # 0.4 x 2.6 pt end ticks
  annotate("text", x = 497.14, y = ypt(458.74), label = "100 kb", hjust = 0, vjust = 0, size = txt(S$axis_title), family = sce2g_font) +
  # genes: intron line, chevrons, exon blocks, italic names
  geom_segment(data = g$introns, aes(x = x, xend = x + len, y = ypt(y), yend = ypt(y)), colour = "#C9C9C9", linewidth = pt_mm(.4)) +
  geom_segment(data = chev, aes(x = x, xend = x + s * 2.711, y = ypt(y), yend = ypt(y) + 1.356), colour = "#C9C9C9", linewidth = pt_mm(.43)) +
  geom_segment(data = chev, aes(x = x, xend = x + s * 2.711, y = ypt(y), yend = ypt(y) - 1.356), colour = "#C9C9C9", linewidth = pt_mm(.43)) +
  geom_rect(data = ex, aes(xmin = x, xmax = xmax, ymin = ymin, ymax = ymax), fill = "black", colour = "black", linewidth = pt_mm(.363)) +
  annotate("text", x = c(330.35, 401.31), y = ypt(c(469.04, 467.17)), label = c("TFRC", "LINC00885"), fontface = "italic",
           hjust = 0, vjust = 0, size = txt(S$axis_title), family = sce2g_font) +
  # arcs with TSS ticks and element feet
  geom_path(data = arcs %>% mutate(lw = lw(sw)), aes(x, y, group = id, colour = hex, linewidth = lw), lineend = "round") +
  geom_rect(data = feet, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = hex)) +
  geom_segment(data = tss %>% mutate(lw = lw(sw)), aes(x = 327.22, xend = 327.22, y = ypt(row), yend = ypt(row + 0.786), colour = hex, linewidth = lw)) +
  # ATAC track: filled profile on a 0.21 pt baseline
  geom_ribbon(data = atac, aes(x = x, ymin = ypt(533.13), ymax = ypt(533.13) + h), fill = "#006EAE") +
  annotate("segment", x = 304.93, xend = 538.78, y = ypt(533.13), yend = ypt(533.13), colour = "#0072B2", linewidth = pt_mm(.213)) +
  # tested CRISPR elements
  geom_rect(data = els, aes(xmin = x, xmax = x + w, ymin = ypt(y_top + h), ymax = ypt(y_top)), fill = G$crispr_tested_element) +
  # track labels, right-aligned in the label column
  annotate("text", x = LAB_X, y = ypt(c(468.17, 491.85, 507.84, 521.71, 530.81, 538.60, 546.40)),
           label = c("Genes", "(single-cell)", "(single-cell)", "CRISPR positives", "ATAC-seq", "Tested CRISPR", "elements"),
           hjust = 1, vjust = 0, size = txt(), family = sce2g_font) +
  annotate("text", x = LAB_X, y = ypt(c(485.26, 501.24)), label = c("scE2G^Multiome", "scE2G^ATAC"), parse = TRUE,
           hjust = 1, vjust = 0, size = txt(), family = sce2g_font) +
  scale_colour_identity() + scale_fill_identity() + scale_linewidth_identity() +
  coord_cartesian(xlim = c(250, 545), ylim = ypt(c(551, 452)), expand = FALSE, clip = "off") +
  theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

pt2mm <- 25.4 / 72
save_panel(p, "fig3e_locus.pdf", width_mm = 295 * pt2mm, height_mm = 99 * pt2mm)
save_panel(p, "fig3e_locus.png", width_mm = 295 * pt2mm, height_mm = 99 * pt2mm)
cat("wrote fig3e_locus.{pdf,png}\n")
