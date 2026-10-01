# scE2G figure design system for ggplot2.
# Source this file; it reads ../tokens.json relative to itself.
# Provides: sce2g (tokens), theme_sce2g(), theme_embedding(), scale_*_predictor(), scale_*_sce2g_seq(),
#           scale_*_anchored(), scale_*_fraction(),
#           stat_pointinterval_sce2g(), geom_arc_sce2g(), arrow_grob_sce2g(), annotate_arrow_sce2g(), save_panel().

suppressPackageStartupMessages({ library(ggplot2); library(jsonlite) })

# repo root = first ancestor of getwd() (or SCE2G_ROOT) containing tokens.json
.sce2g_root <- local({
  d <- normalizePath(Sys.getenv("SCE2G_ROOT", getwd()))
  while (!file.exists(file.path(d, "tokens.json")) && dirname(d) != d) d <- dirname(d)
  if (!file.exists(file.path(d, "tokens.json"))) stop("tokens.json not found; run from the repo or set SCE2G_ROOT")
  d
})
sce2g <- fromJSON(file.path(.sce2g_root, "tokens.json"), simplifyVector = TRUE)

# ---- helpers ---------------------------------------------------------------
pt_mm <- function(pt) pt / .pt / 0.75   # R lwd 1 = 1/96 in = 0.75 pt, so ggplot "mm" linewidth draws at 0.75 mm; verified on mutool SVG of our PDF
sce2g_font <- "Helvetica"

# named vectors of predictor colours: names are index ("1".."19") and method name
sce2g_predictor_cols <- vapply(sce2g$color$predictor, `[[`, "", "hex")
sce2g_predictor_names <- vapply(sce2g$color$predictor, `[[`, "", "name")
names(sce2g_predictor_cols) <- names(sce2g$color$predictor)

# ---- theme -----------------------------------------------------------------
theme_sce2g <- function(base_size = 6.5, base_family = sce2g_font) {
  s <- sce2g$type$size_pt; w <- sce2g$stroke_pt
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      line = element_line(linewidth = pt_mm(w$axis), colour = "black"),
      text = element_text(colour = "black"),
      axis.line = element_line(linewidth = pt_mm(w$axis)),
      axis.ticks = element_line(linewidth = pt_mm(w$tick)),
      axis.ticks.length = unit(sce2g$mark_pt$tick_length, "pt"),   # 2.4 pt, measured on Fig. 2b/2c
      axis.text = element_text(size = s$tick_label, colour = "black"),
      axis.title = element_text(size = s$axis_title),
      axis.title.y = element_text(margin = margin(r = 4)),
      axis.title.x = element_text(margin = margin(t = 4)),
      plot.title = element_text(size = s$panel_title, hjust = 0.5, margin = margin(b = 4)),
      plot.subtitle = element_text(size = s$annotation, hjust = 0.5),
      plot.tag = element_text(size = s$panel_letter, face = "bold"),
      plot.tag.position = c(0, 1),
      legend.text = element_text(size = s$legend),
      legend.title = element_text(size = s$axis_title),
      legend.key.size = unit(sce2g$mark_pt$legend_square, "pt"),
      legend.key.spacing.y = unit(1, "pt"),
      legend.background = element_blank(),
      legend.margin = margin(0, 0, 0, 0),
      strip.background = element_blank(),
      strip.text = element_text(size = s$panel_title),
      panel.spacing = unit(sce2g$layout$gutter_h_pt[1], "pt"),
      plot.margin = margin(4, 6, 4, 4),
      plot.background = element_blank()
    )
}

# Embedding add-on (UMAP / PCA / t-SNE): arrowed, corner-cropped axis guides, no ticks or
# text, square aspect. A list of layers ADDED on top of theme_sce2g(), not a replacement.
# ---- the paper's arrow (svgkit.arrow_poly in R): a stem running into a swept head whose back edge is
# notched toward the tip. Head proportions measured on the Fig. 3f axes and the Fig. 5h variant arrow:
# width 0.82 x length, notch 0.24 x length. grid::arrow() cannot draw the notch, hence a polygon.
arrow_grob_sce2g <- function(x0, y0, dx, dy, head_len = 3, stem = 0.5, colour = "black") {
  # x0, y0: grid units of the stem start; dx, dy: vector to the TIP in pt
  L <- head_len; W <- 0.82 * L; N <- 0.24 * L
  len <- sqrt(dx^2 + dy^2); ux <- dx / len; uy <- dy / len; px <- -uy; py <- ux
  base <- c(dx - L * ux, dy - L * uy); notch <- c(dx - (L - N) * ux, dy - (L - N) * uy)
  head <- rbind(c(dx, dy), base + W / 2 * c(px, py), notch, base - W / 2 * c(px, py))
  pt <- function(v) grid::unit(v, "pt")
  grid::gList(
    grid::segmentsGrob(x0, y0, x0 + pt(notch[1]), y0 + pt(notch[2]),
                       gp = grid::gpar(col = colour, lwd = stem / 0.75, lineend = "butt")),   # lwd 1 = 0.75 pt
    grid::polygonGrob(x0 + pt(head[, 1]), y0 + pt(head[, 2]), gp = grid::gpar(fill = colour, col = NA)))
}
# Arrow of `length_pt` whose tip sits at data point (x, y), pointing at `angle` degrees (90 = up).
annotate_arrow_sce2g <- function(x, y, length_pt, angle = 90, head_len = 3, stem = 0.5, colour = "black") {
  a <- angle * pi / 180; o <- grid::unit(0.5, "npc")
  g <- arrow_grob_sce2g(o - grid::unit(length_pt * cos(a), "pt"), o - grid::unit(length_pt * sin(a), "pt"),
                        length_pt * cos(a), length_pt * sin(a), head_len, stem, colour)
  annotation_custom(grid::gTree(children = g), xmin = x, xmax = x, ymin = y, ymax = y)
}

theme_embedding <- function(arrow_pt = 15) {
  # Measured on Fig. 3f: 15 pt arrows (~14.5 % of a 106 pt panel) from the panel's bottom-left corner, 0.5 pt
  # stems, notched heads 3.0 long x 2.4 wide, 5.5 pt labels starting at the corner ~3 pt clear of the lines.
  o <- grid::unit(0, "npc")
  axes <- grid::gTree(children = grid::gList(arrow_grob_sce2g(o, o, arrow_pt, 0), arrow_grob_sce2g(o, o, 0, arrow_pt)))
  list(
    annotation_custom(axes),
    theme(axis.ticks = element_blank(), axis.text = element_blank(), axis.line = element_blank(), aspect.ratio = 1,
          axis.title = element_text(hjust = 0, size = sce2g$type$size_pt$tick_label),
          axis.title.x = element_text(margin = margin(t = 3)), axis.title.y = element_text(margin = margin(r = 3))),
    labs(x = "UMAP1", y = "UMAP2"))
}

# ---- scales ----------------------------------------------------------------
# Map predictor index (character "1".."19") or name to the fixed colour.
.predictor_lookup <- function(x) {
  x <- as.character(x)
  out <- sce2g_predictor_cols[x]
  miss <- is.na(out)
  if (any(miss)) out[miss] <- sce2g_predictor_cols[match(x[miss], sce2g_predictor_names)]
  unname(out)
}
scale_colour_predictor <- function(...) scale_colour_manual(values = c(sce2g_predictor_cols, setNames(sce2g_predictor_cols, sce2g_predictor_names)), ...)
scale_fill_predictor   <- function(...) scale_fill_manual(values = c(sce2g_predictor_cols, setNames(sce2g_predictor_cols, sce2g_predictor_names)), ...)
scale_color_predictor  <- scale_colour_predictor

scale_fill_sce2g_seq   <- function(palette = "red", ...) scale_fill_gradientn(colours = sce2g$color$sequential[[palette]], ...)
scale_colour_sce2g_seq <- function(palette = "red", ...) scale_colour_gradientn(colours = sce2g$color$sequential[[palette]], ...)

# ---- axis expansion rules (DESIGN_SYSTEM.md §1.5) --------------------------
# Anchored marks (bars/areas): zero gap at the baseline, 5 % headroom at the far end.
scale_y_anchored <- function(...) scale_y_continuous(expand = expansion(mult = c(0, .05)), ...)
scale_x_anchored <- function(...) scale_x_continuous(expand = expansion(mult = c(0, .05)), ...)
# Unanchored marks on a 0-1 metric: quarter breaks, default 5 % padding both ends (never clip at 0).
scale_y_fraction <- function(full = FALSE, limits = if (full) c(0, 1) else NULL, ...) scale_y_continuous(breaks = seq(0, 1, .25), limits = limits, ...)
scale_x_fraction <- function(full = FALSE, limits = if (full) c(0, 1) else NULL, ...) scale_x_continuous(breaks = seq(0, 1, .25), limits = limits, ...)

# ---- point interval (ggdist-style: median dot, 66% thick, 95% thin) --------
# Widths measured on Fig. 5b (printed at 95 %): 0.77 pt / 1.79 pt lines, 3.05 pt dot.
stat_pointinterval_sce2g <- function(colour = "black", ...) list(
  stat_summary(fun.min = function(z) quantile(z, .025), fun.max = function(z) quantile(z, .975),
               geom = "linerange", colour = colour, linewidth = pt_mm(0.77), ...),
  stat_summary(fun.min = function(z) quantile(z, .17), fun.max = function(z) quantile(z, .83),
               geom = "linerange", colour = colour, linewidth = pt_mm(1.79), ...),
  stat_summary(fun = median, geom = "point", colour = colour, size = 3.05 / .pt / 0.75, stroke = 0, shape = 16, ...))

# ---- locus arcs: data frame with x, xend, score -> quadratic arcs of one constant height ----
# The paper's arcs (Fig. 1b-d, 3e, 5f) all peak at the same height in a row regardless of span, so a short
# element-to-TSS link is a narrow full-height loop, not a small bump. Measured on Fig. 3e: 14.8 pt peaks
# (0.95 of the 15.6 pt row pitch); the quadratic Bezier's control point sits at 2 x height.
arc_paths <- function(df, height = 1, n = 60) {
  do.call(rbind, lapply(seq_len(nrow(df)), function(i) {
    t <- seq(0, 1, length.out = n); x0 <- df$x[i]; x1 <- df$xend[i]; xc <- (x0 + x1) / 2
    data.frame(id = i, x = (1 - t)^2 * x0 + 2 * (1 - t) * t * xc + t^2 * x1, y = 4 * height * (1 - t) * t,
               score = df$score[i], df[i, setdiff(names(df), c("x", "xend", "score")), drop = FALSE])
  }))
}
geom_arc_sce2g <- function(df, height = 1, palette = "arc_multiome") list(
  geom_path(data = arc_paths(df, height), aes(x, y, group = id, colour = score, linewidth = score, alpha = score),
            lineend = "round"),
  scale_colour_gradientn(colours = sce2g$color$genomic[[palette]], limits = c(0, 1), guide = "none"),
  scale_linewidth(range = pt_mm(c(sce2g$stroke_pt$arc_min, sce2g$stroke_pt$arc_max)), limits = c(0, 1), guide = "none"),
  scale_alpha(range = c(.4, 1), limits = c(0, 1), guide = "none"))

# ---- significance bracket -------------------------------------------------
# y, tip and line are in data units: tip = bracket tick length, line = one text line height.
# inset pulls both bracket ends in (x units) so adjacent brackets sharing a category do not touch.
# Bold stars sit directly on the bracket; the P value is bottom-aligned 1.5 lines above it,
# so labels on brackets that share a y share a baseline regardless of superscripts.
sig_bracket <- function(x1, x2, y, label, p = NULL, tip = 0.25, line = 1, star_size = sce2g$type$size_pt$stars, inset = 0.06) {
  s <- sce2g$type$size_pt; xm <- (x1 + x2) / 2; x1 <- x1 + inset; x2 <- x2 - inset
  c(list(annotate("segment", x = x1, xend = x2, y = y, yend = y, linewidth = pt_mm(sce2g$stroke_pt$bracket)),
         annotate("segment", x = c(x1, x2), xend = c(x1, x2), y = y, yend = y - tip, linewidth = pt_mm(sce2g$stroke_pt$bracket)),
         annotate("text", x = xm, y = y - tip * .2, label = label, vjust = 0, size = star_size / .pt, family = sce2g_font, fontface = "bold")),
    if (!is.null(p)) list(annotate("text", x = xm, y = y + 1.5 * line, label = p, vjust = 0, parse = TRUE,
                                   size = s$annotation / .pt, family = sce2g_font)))
}

# ---- output ----------------------------------------------------------------
save_panel <- function(p, file, width_mm = sce2g$layout$full_width_mm, height_mm = 60) {
  ext <- tools::file_ext(file)
  if (ext == "pdf") {
    ggsave(file, p, width = width_mm, height = height_mm, units = "mm", device = cairo_pdf)
  } else {
    # cairo, not ragg: ragg silently drops glyphs under ~8 px (1.9 pt superscripts at 300 dpi)
    ggsave(file, p, width = width_mm, height = height_mm, units = "mm", dpi = 300, device = png, type = "cairo", bg = "white")
  }
  invisible(file)
}
