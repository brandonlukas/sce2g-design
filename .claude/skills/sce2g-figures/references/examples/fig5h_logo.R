# Fig 5h — TP73 motif (JASPAR MA0861.1) as an information-content logo, with the reference and alternate
# alleles of rs7696969 written beneath in monospace and the variant base marked by an arrow.
#
# Measured from the published vector panel (Fig. 5 is printed at 95 %; values here are / 0.95):
#   logo   : ggseqlogo letter polygons; panel 37.2 x 17.7 mm; y 0-2 bits with no expansion, breaks 0.5;
#            18 positions, letters 0.95 of a slot; A #00904A, C #1A5F9D, G #F1AD30, T #D92641; per-position
#            letter heights digitised from the polygon boxes into figures/data/fig5h_logo.csv;
#   text   : 4.5 pt tick labels, 6.5 pt "Bits" and centred title; sequences 8 pt Courier, #8B98AC, 6.1 pt
#            advance (not aligned to the logo slots), variant base black (ref) / #9A1D20 (alt); "Ref:"/"Alt:"
#            6 pt; FIMO scores 4 pt; rsID 6.2 pt under a 5 pt arrow with the paper's notched 3 pt head.
#
#   Rscript fig5h_logo.R   # -> fig5h_logo.{pdf,png}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork); library(ggseqlogo) })
source(Sys.getenv("SCE2G_THEME", "../../assets/theme_sce2g.R"))
S <- sce2g$type$size_pt; k <- 1 / 0.95
txt <- function(size) size / .pt

logo <- read.csv("data/fig5h_logo.csv")
m <- t(as.matrix(logo[, c("A", "C", "G", "T")])); colnames(m) <- logo$position
cols <- make_col_scheme(chars = c("A", "C", "G", "T"), cols = c("#00904A", "#1A5F9D", "#F1AD30", "#D92641"))

# sequences: 20 letters at a 5.8 pt printed advance starting 7.7 pt left of slot 1 (x in slot units, y in bits)
seqs <- c(ref = "ATGCAAGTCTTTACTGGTCT", alt = "ATGCAAGTCTTGACTGGTCT")
sx <- ((401.7 + 5.8 * 0:19) - 409.39) / 5.57 + 1
seq_df <- bind_rows(lapply(names(seqs), function(n) tibble(row = n, letter = strsplit(seqs[[n]], "")[[1]], x = sx, i = 0:19))) %>%
  mutate(y = ifelse(row == "ref", -0.616, -0.877),
         colour = case_when(i != 11 ~ "#8B98AC", row == "ref" ~ "black", TRUE ~ "#9A1D20"))
lab <- function(x, y, label, size, ...) annotate("text", x = x, y = y, label = label, size = txt(size), family = sce2g_font, vjust = 0, ...)

p <- ggplot() +
  geom_logo(m, method = "custom", seq_type = "dna", col_scheme = cols) +
  geom_text(data = seq_df, aes(x, y, label = letter, colour = colour), size = txt(8), family = "Courier", hjust = 0.5, vjust = 0) +
  scale_colour_identity() +
  lab(-3.83, -0.616, "Ref:", 6, hjust = 0) + lab(-3.83, -0.877, "Alt:", 6, hjust = 0) +
  lab(20.05, -0.575, "(9.981)", 4, hjust = 0) + lab(20.05, -0.836, "(13.231)", 4, hjust = 0) +
  annotate_arrow_sce2g(11.03, -0.93, length_pt = 4.74 * k, head_len = 2.84 * k, stem = 0.475 * k) +   # the paper's notched head
  lab(11.03, -1.35, "rs7696969", 6.2) +
  scale_x_continuous(breaks = 1:18) +
  scale_y_continuous(breaks = seq(0, 2, 0.5), labels = c("0", "0.5", "1.0", "1.5", "2.0")) +
  coord_cartesian(xlim = c(0.55, 18.51), ylim = c(0, 2), expand = FALSE, clip = "off") +
  labs(x = NULL, y = "Bits", title = "TP73 motif") +
  theme_sce2g() + theme(axis.text = element_text(size = S$annotation), axis.ticks.x = element_line(), plot.margin = margin(4, 30, 42, 4))
fig <- p + plot_layout(widths = unit(100.065 * k, "pt"), heights = unit(47.598 * k, "pt"))
save_panel(fig, "fig5h_logo.pdf", width_mm = 58, height_mm = 38)
save_panel(fig, "fig5h_logo.png", width_mm = 58, height_mm = 38)
cat("wrote fig5h_logo.{pdf,png}\n")
