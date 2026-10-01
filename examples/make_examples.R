# Generates example panels in the scE2G style from simulated data.
# Run: Rscript examples/make_examples.R   (from the repo root)
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr); library(patchwork); library(ggbeeswarm); library(ggseqlogo) })
source("R/theme_sce2g.R")  # run from the repo root
set.seed(1)
theme_set(theme_sce2g())
S <- sce2g$type$size_pt; N <- sce2g$color$neutral; M <- sce2g$mark_pt
txt <- function(size = S$tick_label) list(size = size / .pt, family = sce2g_font)

# ---- a. precision-recall curves --------------------------------------------
pr_curve <- function(id, top, decay) {
  r <- seq(0, 1, .01); p <- top * (1 - r^decay) + .05
  data.frame(id = id, recall = r, precision = pmin(1, pmax(.04, p + cumsum(rnorm(101, 0, .006)))))
}
pars <- tribble(~id, ~top, ~decay, "1", .90, 4, "2", .85, 3.5, "5", .78, 3, "3", .75, 2.5, "7", .70, 1.2, "9", .45, .8, "12", .30, .6)
pr <- bind_rows(Map(pr_curve, pars$id, pars$top, pars$decay))
thr <- pr %>% group_by(id) %>% slice(which.min(abs(precision - 0.5 - as.numeric(id) / 100))) %>% ungroup()
ends <- pr %>% filter(recall == .78) %>% arrange(desc(precision)) %>% mutate(ly = seq(.66, .18, length.out = n()))
pts <- data.frame(id = c("16", "17"), recall = c(.55, .5), precision = c(.48, .5))
p_a <- ggplot(pr, aes(recall, precision, colour = id)) +
  geom_hline(yintercept = .05, linetype = "44", colour = N$chance_line, linewidth = pt_mm(.5)) +
  geom_vline(xintercept = .7, linetype = "44", linewidth = pt_mm(.5)) +
  geom_path(linewidth = pt_mm(.42)) +
  geom_point(data = thr, size = M$threshold_dot_diameter / .pt / 0.75, shape = 16, stroke = 0) +
  geom_point(data = pts, size = M$threshold_dot_diameter / .pt / 0.75, shape = 16, stroke = 0) +
  geom_segment(data = ends, aes(x = recall, xend = 1.08, y = precision, yend = ly), linetype = "62", linewidth = pt_mm(.42)) +
  geom_text(data = ends, aes(x = 1.1, y = ly, label = id), hjust = 0, size = txt()$size, family = sce2g_font) +
  scale_colour_predictor(guide = "none") +
  scale_x_fraction(limits = c(0, 1.15)) + scale_y_fraction(full = TRUE) +   # PR axes span 0-1 (x extended for labels)
  labs(x = expression(Recall ~ (italic(n) == 461)), y = "Precision", title = "CRISPR benchmark (training)") +
  coord_cartesian(clip = "off")

# ---- b. grouped bars with CI, index tick labels and group table ------------
ids <- c("1", "2", "3", "5", "6", "15", "4", "7", "8", "9", "12")
grp <- tibble(dist = c("All", "[0,20)", "[20,100)", "[100,1000)"), pos_n = c(461, 195, 189, 77), neg_n = c(9225, 224, 1188, 7813), gi = 1:4,
              scale = c(1, 1.2, .9, .4))
bars <- expand_grid(id = ids, gi = 1:4) %>% left_join(grp, by = "gi") %>%
  mutate(auprc = pmax(.03, (12 - match(id, ids)) / 16 * scale + rnorm(n(), 0, .03)), ci = runif(n(), .02, .06))
all_order <- bars %>% filter(gi == 1) %>% arrange(desc(auprc)) %>% pull(id)      # rank in "All", kept in every bin
bars <- bars %>% mutate(pos = match(id, all_order) + (gi - 1) * (length(ids) + 1.7))
grp <- grp %>% mutate(xc = (gi - 1) * (length(ids) + 1.7) + (length(ids) + 1) / 2)
tbl_y <- c(-0.09, -0.145, -0.20)   # three annotation rows in data units below the axis (panel is clipped at 0)
p_b <- ggplot(bars, aes(pos, auprc, fill = id)) +
  geom_col(width = 1) +   # grouped bars touch within a group
  geom_errorbar(aes(ymin = auprc - ci, ymax = auprc + ci), width = M$error_cap_fraction_of_bar * .9, linewidth = pt_mm(.3)) +
  annotate("segment", x = -Inf, xend = Inf, y = 0, yend = 0, linewidth = pt_mm(sce2g$stroke_pt$axis)) +  # explicit baseline
  scale_fill_predictor(guide = "none") +
  geom_text(aes(y = -0.015, label = id), size = S$floor / .pt, family = sce2g_font, vjust = 1) +
  scale_x_continuous(expand = expansion(add = c(1, .5))) +
  scale_y_continuous(breaks = seq(0, .75, .25)) +
  annotate("text", x = grp$xc, y = tbl_y[1], label = grp$dist, size = txt()$size, family = sce2g_font) +
  annotate("text", x = grp$xc, y = tbl_y[2], label = grp$pos_n, size = txt()$size, family = sce2g_font) +
  annotate("text", x = grp$xc, y = tbl_y[3], label = scales::comma(grp$neg_n), size = txt()$size, family = sce2g_font) +
  annotate("text", x = -0.5, y = tbl_y, label = c("Distance to TSS (kb):", "CRISPR positives:", "CRISPR negatives:"), hjust = 1,
           size = txt(S$axis_title)$size, family = sce2g_font) +
  labs(x = NULL, y = "Performance (AUPRC)", title = "Prediction performance versus distance to TSS\n(CRISPR, training)") +
  # anchored bars: baseline at 0 with no gap (coord pins the range because the table annotations live below 0)
  coord_cartesian(clip = "off", ylim = c(0, .95), expand = FALSE) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = margin(4, 6, 30, 44))

# ---- c. violin + point interval + brackets ---------------------------------
cat_cols <- unlist(sce2g$color$category_match[c("unmatched", "same_lineage", "same_cell_type")])
vio <- bind_rows(data.frame(g = "Unmatched", v = rgamma(2259, 20, 4)), data.frame(g = "Same lineage", v = rgamma(204, 14, 2.1)),
                 data.frame(g = "Same cell type", v = rgamma(97, 16, 2.1))) %>% mutate(g = factor(g, unique(g)))
names(cat_cols) <- levels(vio$g)
p_c <- ggplot(vio, aes(g, v, fill = g)) + geom_violin(colour = NA, width = .8, scale = "width") +
  stat_pointinterval_sce2g() + scale_fill_manual(values = cat_cols, guide = "none") +
  sig_bracket(1, 2, 14, "***", "italic(P) == 4.39 %*% 10^-41", tip = .3, line = 1.3) + sig_bracket(2, 3, 14, "*", "italic(P) == 0.0324", tip = .3, line = 1.3) +
  sig_bracket(1, 3, 17.6, "***", "italic(P) == 2.45 %*% 10^-23", tip = .3, line = 1.3) +
  annotate("text", x = 1:3, y = -4.4, label = c("italic(n) == \"2,259\"", "italic(n) == 204", "italic(n) == 97"), parse = TRUE, size = txt()$size, family = sce2g_font) +
  scale_y_continuous(breaks = seq(0, 15, 5)) +
  labs(x = "Prediction cell type or eQTL biosample category", y = "Enrichment\n(eQTLs versus common variants)") +
  coord_cartesian(clip = "off", ylim = c(0, 20.5)) + theme(axis.title.x = element_text(margin = margin(t = 14)))

# ---- d. point-range + grey beeswarm ----------------------------------------
pr_d <- tibble(cell = c("Normoblasts", "Erythroblasts", "Proerythroblasts", "MEPs"), est = c(13.5, 12, 11.5, 11),
               lo = est - c(2, 1.5, 1.5, 1.5), hi = est + c(2.6, 1.6, 2, 1.7), sig = c("***", "*", "*", ""),
               col = c(rep(sce2g$color$category_match$same_cell_type, 3), sce2g$color$category_match$same_lineage))
other <- tibble(cell = "Other biosamples", est = rnorm(38, 5.5, 1))
lv <- c("Other biosamples", rev(pr_d$cell))
p_d <- ggplot() +
  geom_quasirandom(data = other, aes(x = est, y = factor(cell, lv)), colour = sce2g$color$category_match$other_points,
                   size = M$beeswarm_dot_diameter / .pt / 0.75, shape = 16, stroke = 0, width = .3, orientation = "y") +
  geom_errorbar(data = pr_d, aes(y = factor(cell, lv), xmin = lo, xmax = hi, colour = col), width = .25, linewidth = pt_mm(.6), orientation = "y") +
  geom_point(data = pr_d, aes(x = est, y = factor(cell, lv), colour = col), size = M$point_range_dot_diameter / .pt / 0.75, shape = 16, stroke = 0) +
  geom_text(data = pr_d, aes(x = hi + .6, y = factor(cell, lv), label = sig), hjust = 0, size = txt(7)$size, family = sce2g_font) +
  scale_colour_identity() + scale_x_continuous(limits = c(0, 18), breaks = seq(0, 15, 5)) +
  labs(x = "Enrichment\nGWAS variants versus common variants\nn = 1,109 MCV variants", y = NULL) +
  scale_y_discrete()

# ---- e. scatter with reference lines; lines and points share one legend -----
sc <- tibble(x = runif(180, 4.9, 8.5)) %>% mutate(y = .62 / (1 + exp(-(x - 5.6) * 2.2)) + rnorm(n(), 0, .025), model = "1")
refs <- tibble(model = c("5", "16"), y = c(.48, .43))
model_cols <- sce2g_predictor_cols[c("5", "16", "1")]
model_labs <- c(expression(ABC^"A=scATAC, C=power law"), "Distance to TSS", expression(scE2G^Multiome))
p_e <- ggplot(sc, aes(x, y)) +
  geom_hline(data = refs, aes(yintercept = y, colour = model), linetype = "44", linewidth = pt_mm(.5)) +
  geom_vline(xintercept = log10(2e6), linetype = "44", linewidth = pt_mm(.5)) +
  geom_point(aes(colour = model), alpha = M$scatter_alpha, size = M$scatter_dot_diameter / .pt / 0.75, shape = 16, stroke = 0) +
  scale_colour_manual(values = model_cols, breaks = names(model_cols), labels = model_labs, name = "Model") +
  # one legend for both geoms: dashed line + no point for baselines, point + no line for scE2G
  guides(colour = guide_legend(override.aes = list(linetype = c("44", "44", "blank"), shape = c(NA, NA, 16), alpha = 1,
                                                   size = M$threshold_dot_diameter / .pt / 0.75))) +
  scale_x_continuous(breaks = 5:9) + scale_y_continuous(breaks = seq(0, .6, .2)) +   # unanchored: default 5 % padding, no clipping at 0
  labs(x = expression(log[10] * "(number of unique ATAC fragments)"), y = "AUPRC") +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.key.width = unit(10, "pt"),
        legend.key.spacing.x = unit(6, "pt"), legend.text = element_text(margin = margin(l = 2)))

# ---- f. UMAP (categorical / continuous / binary) --------------------------
um <- bind_rows(lapply(1:6, function(k) { a <- k / 6 * 2 * pi; tibble(cl = k, x = 3 * cos(a) + rnorm(500, 0, .9), y = 3 * sin(a) + rnorm(500, 0, .9)) })) %>%
  mutate(expr = pmax(0, 3 * exp(-((x - 2)^2 + (y - 1)^2) / 8) + rnorm(n(), 0, .3)), open = runif(n()) < plogis(expr - 1.5))
cl_cols <- c("#E07AA0", "#8FBF5E", "#6FB0D8", "#F0A050", "#B48FD8", "#5BC0B0")
base_um <- function(d) ggplot(d, aes(x, y)) + theme_embedding() +
  theme(legend.position = "right", legend.key.width = unit(4, "pt"), legend.key.height = unit(6, "pt"))
p_f1 <- base_um(um) + geom_point(aes(colour = factor(cl)), size = .25, shape = 16) + scale_colour_manual(values = cl_cols, guide = "none")
p_f2 <- base_um(um) + geom_point(aes(colour = expr), size = .25, shape = 16) + scale_colour_sce2g_seq("umap_expression", name = "ln(TP10k+1)")
p_f3 <- base_um(um) + geom_point(aes(colour = open), size = .25, shape = 16) +
  scale_colour_manual(values = c(`FALSE` = N$box_fill, `TRUE` = sce2g$color$sequential$umap_binary[2]), labels = c("Closed", "Open"), name = "ATAC state")

# ---- g. box plot (base-R style: dashed whiskers with staples, no outliers) ------
bx <- bind_rows(tibble(g = "Closed", v = rbeta(6225, 2, 2.6) * 8000), tibble(g = "Open", v = rbeta(1596, 2.4, 2) * 8000))
row_lbl <- -480; row_n <- -1000   # tick-label row and cell-number row, data units below a panel pinned at 0
p_g <- ggplot(bx, aes(g, v)) +
  geom_boxplot(fill = N$box_fill, width = .8, linewidth = pt_mm(.6), outliers = FALSE, staplewidth = .5,
               whisker.linetype = "34", lineend = "round", colour = "black", median.linewidth = pt_mm(1)) +
  annotate("text", x = 1.5, y = 8300, label = "italic(P) == 8.88 %*% 10^-17", parse = TRUE, size = txt()$size, family = sce2g_font) +
  annotate("text", x = 1:2, y = row_n, label = c("6,225", "1,596"), size = txt()$size, family = sce2g_font) +
  annotate("text", x = 0.45, y = c(row_lbl, row_n), label = c("ATAC state", "Cell number"), hjust = 1, size = txt(S$axis_title)$size, family = sce2g_font) +
  scale_y_continuous(labels = scales::comma, breaks = seq(0, 8000, 2000)) +
  coord_cartesian(ylim = c(0, 8600), expand = FALSE, clip = "off") +
  labs(x = NULL, y = expression("Rank of " * italic(TFRC) * " gene expression")) +
  theme(plot.margin = margin(4, 6, 12, 40))

# ---- h. symmetric split heatmap --------------------------------------------
ct <- paste0("Cell type ", sprintf("%02d", 1:18)); grp <- rep(1:3, each = 6)
sim <- outer(grp, grp, function(a, b) ifelse(a == b, .45, .12)) + matrix(runif(324, 0, .12), 18); sim <- (sim + t(sim)) / 2
hm <- expand_grid(r = 1:18, c = 1:18) %>% mutate(jac = sim[cbind(r, c)], pear = pmin(.6, jac * 1.1 + .02), half = ifelse(r < c, "Jaccard", ifelse(r > c, "Pearson", "diag")))
blue_f <- scales::col_numeric(sce2g$color$sequential$blue, c(0, .6)); red_f <- scales::col_numeric(sce2g$color$sequential$red, c(0, .6))
hm <- hm %>% mutate(col = case_when(half == "Jaccard" ~ red_f(jac), half == "Pearson" ~ blue_f(pear), TRUE ~ "black"))
p_h_main <- ggplot(hm, aes(c, factor(r, 18:1), fill = col)) + geom_tile() + scale_fill_identity() +
  scale_x_reverse(breaks = 1:18, labels = ct, expand = c(0, 0)) + scale_y_discrete(labels = rev(ct), expand = c(0, 0)) +
  coord_equal() + labs(x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = S$heatmap_row_label_min), axis.text.y = element_text(size = S$heatmap_row_label_min),
        axis.line = element_blank(), axis.ticks = element_blank())
colourbar <- function(pal, title) ggplot(data.frame(v = seq(0, .6, length.out = 60)), aes(1, v, fill = v)) + geom_tile(width = 1) +
  scale_fill_gradientn(colours = sce2g$color$sequential[[pal]], guide = "none") + scale_y_continuous(breaks = seq(0, .6, .1), expand = c(0, 0), position = "right") +
  scale_x_continuous(expand = c(0, 0)) + labs(title = title, x = NULL, y = NULL) +
  theme(axis.line = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.ticks.length = unit(1, "pt"),
        plot.title = element_text(hjust = 0, size = S$tick_label), plot.margin = margin(2, 2, 2, 2))
p_h <- wrap_elements(full = p_h_main + ((colourbar("red", "Jaccard\nsimilarity") / colourbar("blue", "Pearson\ncorrelation") / plot_spacer()) + plot_layout(heights = c(1, 1, .6))) +
  plot_layout(widths = c(10, .5)))

# ---- i. motif logo with ref/alt sequence block --------------------------------
pwm <- matrix(runif(4 * 12, 0, .05), 4, dimnames = list(c("A", "C", "G", "T"), NULL))
for (i in seq_len(12)) pwm[c("C", "A", "T", "G", "T", "G", "A", "C", "A", "T", "G", "T")[i], i] <- c(.9, .7, .6, .9, .6, .4, .5, .9, .7, .6, .9, .5)[i]
pwm <- sweep(pwm, 2, colSums(pwm), "/")
xsc <- scale_x_continuous(breaks = 1:12, limits = c(0.4, 12.6), expand = c(0, 0))
p_logo <- ggseqlogo(pwm, method = "bits", font = "helvetica_regular") +
  scale_fill_manual(values = unlist(sce2g$color$motif), guide = "none") +
  theme_sce2g() + xsc + scale_y_continuous(breaks = seq(0, 2, .5), limits = c(0, 2), expand = c(0, 0)) +
  labs(x = NULL, y = "Bits", title = "TP73 motif") + theme(plot.margin = margin(4, 24, 0, 4))
ref <- strsplit("ACATGTGACATG", "")[[1]]; alt <- ref; alt[10] <- "G"; seq_col <- "#7F8FA6"
seq_df <- bind_rows(tibble(x = 1:12, y = 2, base = ref, col = ifelse(1:12 == 10, "black", seq_col)),
                    tibble(x = 1:12, y = 1.45, base = alt, col = ifelse(1:12 == 10, sce2g$color$genomic$crispr_positive_arc, seq_col)))
p_seq <- ggplot(seq_df, aes(x, y)) +
  geom_text(aes(label = base, colour = col), family = "Courier", size = txt(7)$size) + scale_colour_identity() +
  annotate("text", x = 0.3, y = c(2, 1.45), label = c("Ref:", "Alt:"), hjust = 1, size = txt()$size, family = sce2g_font) +
  annotate("text", x = 12.7, y = c(2, 1.45), label = c("(9.981)", "(13.231)"), hjust = 0, size = txt(S$annotation)$size, family = sce2g_font) +
  annotate("segment", x = 10, xend = 10, y = 0.8, yend = 1.17, arrow = arrow(length = unit(1.5, "pt"), type = "closed"), linewidth = pt_mm(.4)) +
  annotate("text", x = 10, y = 0.45, label = "rs7696969", size = txt()$size, family = sce2g_font) +
  scale_x_continuous(breaks = 1:12) + coord_cartesian(xlim = c(0.4, 12.6), ylim = c(0.2, 2.3), expand = FALSE, clip = "off") +
  theme_void() + theme(plot.margin = margin(0, 24, 2, 4))
p_i <- wrap_elements(full = (p_logo / p_seq) + plot_layout(heights = unit(c(1, 30), c("null", "pt"))))  # sequence block fixed at 30 pt

# ---- j. simplified locus plot ---------------------------------------------
G <- sce2g$color$genomic
genes <- tibble(name = c("EBF1", "LINC02202"), start = c(400, 620), end = c(560, 660), strand = c("-", "+"),
                exons = list(c(400, 430, 470, 520, 552), c(620, 645)))
exon_df <- genes %>% rowwise() %>% reframe(name, x = exons, w = 6)
cell_types <- c("All BMMCs", "HSCs", "B1 B cells")
links <- bind_rows(lapply(seq_along(cell_types), function(k) tibble(ct = cell_types[k], x = sort(runif(c(14, 5, 22)[k], 40, 740)), xend = 555, score = runif(length(x), .15, 1))))
atac <- expand_grid(ct = cell_types, x = seq(0, 800, 2)) %>%
  mutate(y = pmin(1, abs(rnorm(n(), .05, .05)) + ifelse(abs(x - 555) < 4, .9, 0) + ifelse(abs(x - 150) < 3, .5, 0) + ifelse(abs(x - 700) < 3, .4, 0)))
row_h <- 1.6; arc_h <- 1
locus_rows <- lapply(seq_along(cell_types), function(k) {
  off <- -(k - 1) * row_h
  l <- filter(links, ct == cell_types[k]); a <- filter(atac, ct == cell_types[k])
  list(geom_path(data = arc_paths(l, arc_h) %>% mutate(y = y + off + .3), aes(x, y, group = id, colour = score, linewidth = score, alpha = score), lineend = "round"),
       geom_ribbon(data = a %>% mutate(ymax = y * .28 + off, ymin = off), aes(x, ymin = ymin, ymax = ymax), fill = G$atac_track, colour = NA),
       annotate("segment", x = 0, xend = 800, y = off, yend = off, linewidth = pt_mm(.2), colour = G$atac_track),
       annotate("text", x = -10, y = off + .7, label = cell_types[k], hjust = 1, size = txt(6)$size, family = sce2g_font),
       if (k == 1) annotate("text", x = -10, y = off + .12, label = "ATAC-seq", hjust = 1, size = txt(S$annotation)$size, family = sce2g_font))
})
p_j <- ggplot() +
  annotate("segment", x = 550, xend = 800, y = 2.55, yend = 2.55, linewidth = pt_mm(.5)) + annotate("text", x = 810, y = 2.55, label = "250 kb", hjust = 0, size = txt()$size, family = sce2g_font) +
  geom_segment(data = genes, aes(x = start, xend = end, y = 2.05, yend = 2.05), colour = N$hairline, linewidth = pt_mm(.5)) +
  geom_tile(data = exon_df, aes(x + 3, 2.05, width = w, height = .22), fill = "black") +
  geom_text(data = genes, aes(x = ifelse(strand == "-", start - 8, end + 8), y = 2.05, label = name, hjust = ifelse(strand == "-", 1, 0)), fontface = "italic", size = txt(6)$size, family = sce2g_font) +
  annotate("text", x = -10, y = 2.05, label = "Genes", hjust = 1, size = txt(6)$size, family = sce2g_font) +
  geom_vline(xintercept = 150, linetype = "33", colour = G$variant_line, linewidth = pt_mm(.36)) +
  annotate("text", x = 150, y = 2.6, label = "rs0000001", colour = G$gwas_credible_set, size = txt(S$annotation)$size, family = sce2g_font) +
  annotate("tile", x = 150, y = 2.35, width = 30, height = .12, fill = G$gwas_credible_set) +
  locus_rows +
  scale_colour_gradientn(colours = G$arc_multiome, limits = c(0, 1), name = expression(scE2G^Multiome ~ score)) +
  scale_linewidth(range = pt_mm(c(sce2g$stroke_pt$arc_min, sce2g$stroke_pt$arc_max)), limits = c(0, 1), guide = "none") +
  scale_alpha(range = c(.4, 1), limits = c(0, 1), guide = "none") +
  coord_cartesian(xlim = c(-120, 860), ylim = c(-(length(cell_types) - 1) * row_h - .1, 2.7), expand = FALSE, clip = "off") +
  theme_void(base_family = sce2g_font) +
  theme(legend.position = "bottom", legend.direction = "horizontal", legend.key.height = unit(4, "pt"), legend.key.width = unit(14, "pt"),
        legend.title = element_text(size = S$axis_title, vjust = 1), legend.text = element_text(size = S$tick_label), plot.tag = element_text(size = S$panel_letter, face = "bold"))

# ---- assemble --------------------------------------------------------------
dir.create(file.path(.sce2g_root, "examples", "out"), showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(.sce2g_root, "examples", "out", f)
fig1 <- (p_a | p_b) + plot_layout(widths = c(1, 2.2))
fig2 <- (p_c | p_d | p_e) + plot_layout(widths = c(1.2, 1.1, 1))
fig3 <- (p_f1 | p_f2 | p_f3)
fig4 <- (p_h | p_g | p_i) + plot_layout(widths = c(1.5, .7, 1.1))
full <- (fig1 / fig2 / fig3 / p_j / fig4) + plot_layout(heights = c(1.1, 1.1, .9, 1, 1.4)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(size = S$panel_letter, face = "bold", family = sce2g_font))
save_panel(full, out("sce2g_examples.pdf"), height_mm = 250)
save_panel(full, out("sce2g_examples.png"), height_mm = 250)
for (nm in c("p_a", "p_b", "p_c", "p_d", "p_e", "p_g", "p_h", "p_i", "p_j")) save_panel(get(nm) + labs(tag = NULL), out(paste0(nm, ".pdf")), width_mm = if (nm %in% c("p_b", "p_j")) 120 else if (nm == "p_e") 90 else 60, height_mm = 50)
cat("wrote", out(""), "\n")
