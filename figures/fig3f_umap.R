# Fig 3f — K562 multiome UMAP triptych: meta-cell membership, TFRC expression, TFRC-promoter ATAC state.
#
# Layout measured from the published panel (a 319 x 122 pt raster plus vector legends/axes):
#   * three equal panels, no axes except a short arrowed L in the first (15 pt arrows, 0.5 pt, the
#     paper's notched 2.4 x 3 pt heads) with "UMAP1"/"UMAP2" (5.5 pt) starting at the corner;
#   * legends INSIDE the second and third panels, lower right: a 4.1 x 17.6 pt vertical bar with
#     white ticks, the title rotated 90 deg to the LEFT of the bar (5.5 pt), tick labels to its
#     right (4.5 pt); the binary ATAC state is still shown as a grey-to-red gradient bar labelled
#     Closed/Open;
#   * points ~1 pt, white-outlined only in the categorical panel; meta-cell centroids as 2.6 pt black discs with a
#     0.8 pt dot in the meta-cell colour at the centre.
# Cells are simulated to echo the paper's shape (dense body, upper arm, ring with a hole, tail).
#
#   Rscript figures/fig3f_umap.R   # -> figures/fig3f_umap.{pdf,png}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
source("R/theme_sce2g.R")
set.seed(11)
S <- sce2g$type$size_pt; N <- sce2g$color$neutral; SEQ <- sce2g$color$sequential
txt <- function(size = S$tick_label) size / .pt

# ---- simulated embedding ----------------------------------------------------------------------
seg <- function(n, from, to, sd) { t <- runif(n); tibble(x = from[1] + t * (to[1] - from[1]) + rnorm(n, 0, sd), y = from[2] + t * (to[2] - from[2]) + rnorm(n, 0, sd)) }
blob <- function(n, cx, cy, R, shape) {                       # uniform fill of a deformed disc -> crisp outline
  th <- runif(3 * n, 0, 2 * pi); rr <- sqrt(runif(3 * n)) * R * shape(th)
  tibble(x = cx + rr * cos(th), y = cy + rr * sin(th)) %>% slice_head(n = n)
}
ang <- runif(1500, 0, 2 * pi); r <- 2.2 + rnorm(1500, 0, .22)
cells <- bind_rows(
  blob(5200, -3.0, 0.2, 2.9, function(t) 1 + .22 * sin(2 * t) + .12 * cos(3 * t) - .1 * cos(t)),   # body
  seg(800, c(-0.6, 2.4), c(1.3, 5.3), .42),                                  # upper arm
  tibble(x = 3.2 + r * cos(ang), y = -0.4 + r * sin(ang)),                    # ring
  seg(420, c(-2.3, -2.3), c(-0.4, -5.2), .38),                                # tail (attached to the body)
  seg(280, c(5.3, 0.3), c(6.8, 1.6), .32),                                    # right spur
  seg(350, c(-0.2, 0.3), c(1.3, -0.2), .35)) %>%                              # bridge body -> ring
  filter(!((x - 3.2)^2 + (y + 0.4)^2 < 1.75^2)) %>%                           # keep the hole clear
  # expression: moderate mean in the body with per-cell noise (salt-and-pepper, mostly light purple);
  # open chromatin at the promoter in a scattered minority, enriched where expression is high
  mutate(expr = pmin(3, pmax(0, rnorm(n(), 1.5 * exp(-((x + 3)^2 + (y - 0.2)^2) / 16), 0.75))),
         open = runif(n()) < plogis(1.5 * expr - 2.4))
km <- kmeans(cells[, c("x", "y")], centers = 110, iter.max = 30, nstart = 2)
cells$meta <- factor(km$cluster)
centroids <- as_tibble(km$centers) %>% mutate(meta = factor(seq_len(n())))
meta_cols <- hsv(runif(110), runif(110, .45, .7), runif(110, .65, .9))

# ---- shared pieces ------------------------------------------------------------------------------
# categorical panel: 1.1 pt filled circles with a hairline white outline; continuous panels: no outline
pt <- 1.1 / .pt; cell_pt <- list(shape = 21, colour = "white", stroke = 0.08, size = pt)
base <- function(d) ggplot(d, aes(x, y)) + coord_equal(clip = "off")
legend_inside <- theme(legend.position = "inside", legend.position.inside = c(0.97, 0.05), legend.justification = c(1, 0),
                       legend.title.position = "left",
                       legend.title = element_text(angle = 90, size = S$tick_label, hjust = 0.5, margin = margin(r = 1.5)),
                       legend.text = element_text(size = S$annotation, margin = margin(l = 1.5)),
                       legend.key.width = unit(4.1, "pt"), legend.key.height = unit(17.6 / 5, "pt"),
                       legend.ticks.length = unit(0.8, "pt"), legend.margin = margin(0, 0, 0, 0))
no_axes <- list(theme(axis.line = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(),
                      axis.ticks = element_blank(), axis.text = element_blank()), labs(x = NULL, y = NULL))

p1 <- base(cells) + do.call(geom_point, c(list(mapping = aes(fill = meta)), cell_pt)) +
  geom_point(data = centroids, size = 2.6 / .pt / 0.75, colour = "black", shape = 16, stroke = 0) +              # centroid: black disc ...
  geom_point(data = centroids, aes(colour = meta), size = 1.0 / .pt / 0.75, shape = 16, stroke = 0) +           # ... with a dot in its colour
  scale_fill_manual(values = meta_cols, guide = "none") + scale_colour_manual(values = meta_cols, guide = "none") +
  theme_sce2g() + theme_embedding()
p2 <- base(cells) + geom_point(aes(colour = expr), size = pt, shape = 16) +                 # continuous: no outline
  scale_colour_gradientn(colours = SEQ$umap_expression, limits = c(0, 3), breaks = 0:3, name = "ln(TP10k+1)") +
  theme_sce2g() + theme(aspect.ratio = 1) + no_axes + legend_inside +
  guides(colour = guide_colourbar(ticks.colour = "white", ticks.linewidth = pt_mm(.25), frame.colour = NA))
p3 <- base(cells) + geom_point(aes(colour = as.numeric(open)), size = pt, shape = 16) +
  scale_colour_gradientn(colours = SEQ$umap_binary, limits = c(0, 1), breaks = c(0, 1), labels = c("Closed", "Open"), name = "ATAC state") +
  theme_sce2g() + theme(aspect.ratio = 1) + no_axes + legend_inside +
  guides(colour = guide_colourbar(ticks = FALSE, frame.colour = NA))

fig <- (p1 | p2 | p3) & theme(plot.margin = margin(2, 2, 2, 2))
save_panel(fig, "figures/fig3f_umap.pdf", width_mm = 113, height_mm = 43)
save_panel(fig, "figures/fig3f_umap.png", width_mm = 113, height_mm = 43)
cat("wrote figures/fig3f_umap.{pdf,png}\n")
