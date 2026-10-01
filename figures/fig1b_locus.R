# Fig 1b — scE2G-Multiome predictions for EBF1 across BMMC cell types, as a locus plot.
#
# Layout (top to bottom), measured from the published panel (single column, 88 mm wide):
#   scale bar -> gene models (three rows) -> five cell-type rows, each arcs over a normalised
#   ATAC track -> 49-row cell-type x element heatmap with a grey band on the B-lineage rows ->
#   horizontal score colour bar. Arc colour AND width scale with the scE2G score; heatmap cells
#   use the same red ramp. Data here are simulated with the paper's structure (B-lineage rows carry
#   many distal enhancers, everything else mostly promoter-proximal); swap in real predictions by
#   replacing `elements`, `scores` and `atac`.
#
#   Rscript figures/fig1b_locus.R   # -> figures/fig1b_locus.{pdf,png}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(patchwork) })
source("R/theme_sce2g.R")
set.seed(7)
S <- sce2g$type$size_pt; G <- sce2g$color$genomic; N <- sce2g$color$neutral
txt <- function(size = S$tick_label) size / .pt

# ---- genomic frame -----------------------------------------------------------------------------
# The published track spans 215.5 pt for chr1:158,330,000-159,250,000 (hg38); we keep that mapping so
# gene/element positions measured in pt convert directly to coordinates.
x0_pt <- 97.4; bp_per_pt <- 800e3 / (289.2 - 103.8)
gpos <- function(x_pt) 158330000 + (x_pt - x0_pt) * bp_per_pt
region <- c(gpos(97.4), gpos(312.9))

# label side follows Gviz: to the left of the gene unless that runs off the region (then right);
# the focal gene gets the larger label
genes <- tribble(
  ~name,        ~row, ~start_pt, ~end_pt, ~strand, ~side,   ~size,
  "LINC02227",   1,    99.7,      116.2,   "-",     "right", 3.3,
  "EBF1",        1,    182.2,     276.1,   "-",     "left",  5.1,
  "LINC01932",   1,    283.0,     288.0,   "+",     "right", 3.3,
  "LINC02202",   2,    275.7,     279.9,   "+",     "left",  3.3,
  "RNF145",      3,    289.0,     300.2,   "+",     "left",  3.3
) %>% mutate(start = gpos(start_pt), end = gpos(end_pt))
# exon (start, end) in pt, measured from the PDF; narrow ones widened to a 0.9 pt minimum so they print
exon_pt <- list(LINC02227 = c(99.7, 100.4, 104.6, 105.3, 106.2, 106.9, 115.6, 116.2),
                EBF1 = c(182.2, 183.2, 184.7, 185.4, 185.7, 186.8, 190.0, 190.7, 200.8, 201.4, 205.2, 205.8, 211.4, 212.0,
                         215.3, 215.9, 269.4, 270.0, 272.0, 272.6, 274.5, 276.1),
                LINC01932 = c(283.0, 283.6, 287.0, 288.0), LINC02202 = c(275.7, 276.3, 278.3, 278.9, 279.3, 279.9),
                RNF145 = c(289.0, 290.8, 291.5, 292.3, 292.7, 294.0, 294.5, 295.2, 297.5, 298.1, 299.5, 300.2))
exons <- bind_rows(lapply(names(exon_pt), function(g) { v <- exon_pt[[g]]; tibble(name = g, s_pt = v[c(TRUE, FALSE)], e_pt = v[c(FALSE, TRUE)]) })) %>%
  left_join(genes %>% select(name, row), by = "name") %>%
  mutate(e_pt = pmax(e_pt, s_pt + 0.9), x = gpos((s_pt + e_pt) / 2), w = (e_pt - s_pt) * bp_per_pt)
# strand chevrons along the intron line, every 18 pt, only where there is intron to draw on
chev <- genes %>% filter(end_pt - start_pt > 12) %>% rowwise() %>%
  reframe(name, row, strand, x_pt = seq(start_pt + 5, end_pt - 5, by = 12)) %>%
  anti_join(exons %>% transmute(name, s_pt, e_pt) %>% rowwise() %>% reframe(name, x_pt = seq(floor(s_pt) - 2, ceiling(e_pt) + 2)) %>% distinct(),
            by = c("name", "x_pt")) %>%   # drop chevrons that would sit on an exon
  mutate(x = gpos(x_pt))
tss <- gpos(276.1)   # EBF1 promoter (minus strand: TSS at the right end)

# ---- cell types (heatmap rows, top to bottom) --------------------------------------------------
# sup: a = BMMC source, b = PBMC source (as in the paper's row labels)
rows <- tribble(
  ~name, ~sup, ~lineage,
  "All BMMCs", "", "agg", "All PBMCs", "", "agg", "All pancreatic islets", "", "agg",
  "Proerythroblast", "", "ery", "Erythroblast", "", "ery", "Normoblast", "", "ery",
  "K562 (ref. 71)", "", "line", "K562 (this study)", "", "line", "K562 (ref. 39)", "", "line",
  "WTC11", "", "line", "HCT116", "", "line", "Jurkat", "", "line",
  "Mesenchymal cell (islets)", "", "islet", "Islet endothelial cells", "", "islet", "Islet duct cell", "", "islet",
  "Islet acinar cell", "", "islet", "Islet delta cell", "", "islet", "Islet beta cell", "", "islet",
  "Islet alpha cell", "", "islet", "Islet polypeptide cell", "", "islet", "Islet immune cell", "", "islet",
  "HSC", "", "prog", "GMP", "", "prog", "MEP", "", "prog", "pDC", "a", "myeloid",
  "Naive CD20+ B cell", "a", "B", "B1 B cell", "a", "B", "B cell", "b", "B", "Transitional B cell", "b", "B",
  "Lymphoid progenitor", "", "B", "Plasma cell", "", "B", "GM12878", "", "B",
  "CD14+ monocyte (1)", "b", "myeloid", "CD14+ monocyte (2)", "b", "myeloid", "CD14+ monocyte", "a", "myeloid",
  "CD16+ monocyte", "b", "myeloid", "CD16+ monocyte", "a", "myeloid", "cDC", "b", "myeloid", "cDC2", "a", "myeloid",
  "CD8+ CM T cell", "b", "T", "CD8+ T cell", "a", "T", "CD4+ naive T cell", "a", "T", "CD4+ activated T cell", "a", "T",
  "CD4+ memory T cell", "b", "T", "CD4+ naive T cell", "b", "T", "CD8+ naive T cell", "a", "T",
  "NK cell", "b", "T", "NK cell", "a", "T", "Innate lymphoid cell", "", "T"
) %>% mutate(id = row_number())
band_rows <- range(which(rows$lineage == "B"))   # the shaded B-lineage block

# plotmath label: "CD14+ monocyte (1)" with sup "b" -> paste("CD14"^"+", " monocyte (1)")^"b"
pm_label <- function(name, sup) {
  base <- if (grepl("\\+", name)) {
    p <- strsplit(name, "+", fixed = TRUE)[[1]]
    sprintf('paste("%s"^"+", "%s")', p[1], p[2])
  } else sprintf('"%s"', name)
  if (nzchar(sup)) sprintf("%s^\"%s\"", base, sup) else base
}
row_labels <- parse(text = mapply(pm_label, rows$name, rows$sup))

# ---- elements and scores (simulated) ------------------------------------------------------------
el_pt <- c(seq(101, 118, 2.5), seq(128, 265, 3.2), seq(266.5, 279, 1.6), seq(283, 296, 2.6), 303, 306)
elements <- tibble(el = seq_along(el_pt), x = gpos(el_pt),
                   prox = abs(el_pt - 276.1) < 6, distal_B = el_pt > 125 & el_pt < 268,
                   left = el_pt < 120, right = el_pt > 282)
lineage_p <- c(agg = .32, ery = .03, line = .02, islet = .02, prog = .08, myeloid = .05, B = .55, T = .04)
scores <- crossing(rows %>% select(id, lineage), elements) %>%
  mutate(lineage = ifelse(id == 3, "islet", lineage),
         p = case_when(prox ~ .95, distal_B ~ lineage_p[lineage], left ~ .05, right ~ .08, TRUE ~ .02),
         p = ifelse(id %in% c(13, 22, 23), pmax(p, .22 * distal_B), p),          # mesenchymal, HSC, GMP: scattered hits
         hit = runif(n()) < p,
         score = ifelse(hit, pmin(.95, .18 + rbeta(n(), 2, 5) * .8 + prox * .3), NA_real_)) %>%
  filter(hit)

# ---- five displayed cell-type rows: arcs + ATAC ------------------------------------------------
tracks <- tibble(ct = c("All BMMCs", "HSCs", "Lymphoid progenitors", "B1 B cells", "CD8+ T cells"),
                 lab = c('"All BMMCs"', '"HSCs"', '"progenitors"', '"B1 B cells"', '"CD8"^"+"~"T cells"'),
                 heat_id = c(1, 22, 30, 27, 41), k = 1:5)
arcs <- tracks %>% inner_join(scores %>% select(id, x, score), by = c("heat_id" = "id")) %>%
  filter(!(ct == "CD8+ T cells" & abs(x - tss) > 6e4 & runif(n()) < .5)) %>%
  transmute(k, x = x, xend = tss, score)
grid_bp <- seq(region[1], region[2], by = 2500)
atac <- crossing(tracks %>% select(k), tibble(x = grid_bp)) %>%
  left_join(arcs %>% group_by(k) %>% summarise(pk = list(x), .groups = "drop"), by = "k") %>%
  rowwise() %>%
  mutate(y = 0.04 + abs(rnorm(1, 0, .03)) +
           sum(exp(-((x - unlist(pk))^2) / (2 * 4000^2))) * runif(1, .25, .6) +
           exp(-((x - tss)^2) / (2 * 3500^2)) * .95 +                              # promoter peak in every row
           exp(-((x - gpos(298))^2) / (2 * 3500^2)) * runif(1, .3, .8)) %>%       # RNF145 promoter
  ungroup() %>% mutate(y = pmin(1, y))

# vertical layout of the track block (arbitrary units): each row = arcs (1.0) + ATAC (0.3) + gap
row_h <- 1.55; arc_h <- 1.0; atac_h <- 0.3; foot_h <- 0.055   # 0.786 pt at ~14 pt per unit
y_row <- function(k) 6.2 - (k - 1) * row_h   # row baseline (ATAC baseline)
y_genes <- c(8.75, 8.42, 8.09); y_scale <- 9.25
lbl_x <- region[1] - 12e3
p_tracks <- ggplot() +
  # scale bar
  annotate("segment", x = gpos(103.8), xend = gpos(289.2), y = y_scale, yend = y_scale, linewidth = pt_mm(.25)) +
  annotate("segment", x = c(gpos(103.8), gpos(289.2)), xend = c(gpos(103.8), gpos(289.2)), y = y_scale - .13, yend = y_scale + .13, linewidth = pt_mm(.25)) +
  annotate("text", x = gpos(291), y = y_scale, label = "800 kb", hjust = 0, size = txt(S$annotation), family = sce2g_font) +
  # genes: intron line, exon blocks, italic names
  geom_segment(data = genes, aes(x = start, xend = end, y = y_genes[row], yend = y_genes[row]), colour = "#B7BECB", linewidth = pt_mm(.25)) +
  geom_text(data = chev, aes(x, y_genes[row], label = ifelse(strand == "-", "<", ">")), colour = "#B7BECB", size = txt(4), family = sce2g_font, vjust = .45) +
  geom_tile(data = exons, aes(x = x, y = y_genes[row], width = w, height = .26), fill = "black") +
  geom_text(data = genes, aes(x = ifelse(side == "left", start - 1.6e3 * size, end + 1.6e3 * size), y = y_genes[row], label = name, hjust = ifelse(side == "left", 1, 0),
                              size = size / .pt), fontface = "italic", family = sce2g_font, vjust = .4) +
  scale_size_identity(guide = "none") +
  annotate("text", x = lbl_x, y = mean(y_genes), label = "Genes", hjust = 1, size = txt(S$annotation), family = sce2g_font) +
  # per cell type: arcs, then the ATAC ribbon on its baseline
  geom_path(data = arc_paths(arcs %>% select(x, xend, score, k), arc_h) %>% mutate(y = y + y_row(k) + atac_h + .05),
            aes(x, y, group = interaction(k, id), colour = score, linewidth = score, alpha = score), lineend = "round") +
  geom_tile(data = arcs, aes(x, y_row(k) + atac_h + .05 - foot_h / 2, width = 3000, height = foot_h, fill = score, alpha = score)) +   # 0.8 pt element feet
  scale_fill_gradientn(colours = G$arc_multiome, limits = c(0, 1), guide = "none") +
  geom_ribbon(data = atac %>% mutate(ymin = y_row(k), ymax = y_row(k) + y * atac_h), aes(x = x, ymin = ymin, ymax = ymax, group = k),
              fill = G$atac_track, colour = NA) +
  geom_text(data = tracks, aes(x = lbl_x, y = y_row(k) + atac_h + .45 - (k == 3) * .17, label = lab), hjust = 1, parse = TRUE,
            size = txt(S$annotation), family = sce2g_font) +
  annotate("text", x = lbl_x, y = y_row(3) + atac_h + .45 + .21, label = "Lymphoid", hjust = 1, size = txt(S$annotation), family = sce2g_font) +
  annotate("text", x = lbl_x, y = y_row(1) + .12, label = "ATAC-seq", hjust = 1, size = txt(3.3), family = sce2g_font) +
  scale_colour_gradientn(colours = G$arc_multiome, limits = c(0, 1), guide = "none") +
  scale_linewidth(range = pt_mm(c(sce2g$stroke_pt$arc_min, sce2g$stroke_pt$arc_max)), limits = c(0, 1), guide = "none") +
  scale_alpha(range = c(.35, 1), limits = c(0, 1), guide = "none") +
  coord_cartesian(xlim = region, ylim = c(y_row(5) - .05, y_scale + .25), expand = FALSE, clip = "off") +
  theme_void(base_family = sce2g_font) + theme(plot.margin = margin(2, 8, 1, 2))

# ---- heatmap ------------------------------------------------------------------------------------
el_w <- 3000
p_heat <- ggplot() +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = band_rows[2] + .5, ymax = band_rows[1] - .5, fill = G$row_highlight) +
  geom_tile(data = scores, aes(x, id, fill = score), width = el_w, height = .85) +
  scale_fill_gradientn(colours = sce2g$color$sequential$red[-1], limits = c(.15, .95), breaks = c(.2, .4, .6, .8),
                       name = expression(scE2G^Multiome ~ score), na.value = "white") +
  scale_y_reverse(breaks = rows$id, labels = row_labels, expand = c(0, 0)) +
  coord_cartesian(xlim = region, expand = FALSE, clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_sce2g() +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
        axis.text.y = element_text(size = S$heatmap_row_label_min, margin = margin(r = 2)),
        legend.position = "bottom", legend.direction = "horizontal", legend.title.position = "left",
        legend.title = element_text(size = S$legend, vjust = 0), legend.text = element_text(size = S$annotation),
        legend.key.height = unit(6, "pt"), legend.key.width = unit(10, "pt"), legend.text.position = "top",   # bar = 5 x key.width
        legend.margin = margin(4, 0, 0, 0),
        plot.margin = margin(1, 8, 2, 2)) +
  guides(fill = guide_colourbar(ticks.colour = "white", ticks.linewidth = pt_mm(.3), frame.colour = NA, draw.ulim = FALSE, draw.llim = FALSE))

# ---- assemble ------------------------------------------------------------------------------------
fig <- (p_tracks / p_heat) + plot_layout(heights = c(1.0, 2.05))
dir.create("figures", showWarnings = FALSE)
save_panel(fig, "figures/fig1b_locus.pdf", width_mm = sce2g$layout$single_column_mm, height_mm = 150)
save_panel(fig, "figures/fig1b_locus.png", width_mm = sce2g$layout$single_column_mm, height_mm = 150)
cat("wrote figures/fig1b_locus.{pdf,png}\n")
