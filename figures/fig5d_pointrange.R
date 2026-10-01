# Fig 5d — enrichment of fine-mapped random-glucose GWAS variants in scE2G-Multiome predictions per cell
# type: eight pancreatic cell types as a dot-with-CI plot (islet cell types in dark purple) over a grey
# beeswarm of all other biosamples. Two ggplots stacked with patchwork, sharing one y axis line.
#
# Measured from the published vector panel (Fig. 5 is printed at 95 %; values here are / 0.95):
#   panels : 42.8 mm wide; rows 18.8 mm (8 rows, 6.5 pt each, 0.6-row padding) over 12.4 mm (swarm);
#            x 0-21.7 with the default 5 % expansion, breaks every 5; no x axis on the upper panel;
#   dots   : 6.2 pt, #260032 (islet) / #9A3983 (other), CI line 0.42 pt with 1.3 pt end caps in the same
#            colour; 8.5 pt star 3.8 pt right of the CI for the two significant rows;
#   swarm  : 29 points, 2.9 pt, #A7A9AB, positions digitised (figures/data/fig5d_swarm.csv). The layout is a
#            beeswarm "swarm" packing (every point touches a neighbour at ~1 diameter, 7 sit on the centre
#            line), not quasirandom; ggbeeswarm::geom_beeswarm(method = "swarm", cex = 5.9) at this panel size
#            gives an equivalent but not identical packing, so the digitised positions are used;
#   x title: three lines, "Enrichment" / comparison / italic n, 6.5 pt, 1.2 line spacing.
#
#   Rscript figures/fig5d_pointrange.R   # -> figures/fig5d_pointrange.{pdf,png}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
source("R/theme_sce2g.R")
S <- sce2g$type$size_pt
txt <- function(size = S$tick_label) size / .pt
k <- 1 / 0.95                                             # print -> design scale

d <- read.csv("figures/data/fig5d_values.csv") %>% mutate(cell_type = factor(cell_type, rev(cell_type)), colour = ifelse(group == "islet", "#260032", "#9A3983"))
swarm <- read.csv("figures/data/fig5d_swarm.csv")
xlim <- c(-1.086, 22.82); cap <- 0.1                      # cap half-height in row units (1.3 pt of a 6.5 pt row)
common <- list(scale_x_continuous(breaks = seq(0, 20, 5)), scale_colour_identity(), theme_sce2g())

p_top <- ggplot(d, aes(y = cell_type, colour = colour)) +
  geom_segment(aes(x = lo, xend = hi, yend = cell_type), linewidth = pt_mm(0.401 * k)) +
  geom_segment(aes(x = lo, xend = lo, y = as.numeric(cell_type) - cap, yend = as.numeric(cell_type) + cap), linewidth = pt_mm(0.401 * k)) +
  geom_segment(aes(x = hi, xend = hi, y = as.numeric(cell_type) - cap, yend = as.numeric(cell_type) + cap), linewidth = pt_mm(0.401 * k)) +
  geom_point(aes(x = enrichment), size = 5.92 * k / .pt / 0.75, stroke = 0, shape = 16) +
  geom_text(data = filter(d, group == "islet"), aes(x = hi + 0.75, y = as.numeric(cell_type) - 0.7, label = "*"),
            colour = "black", size = txt(8.5), family = sce2g_font, hjust = 0, vjust = 0) +
  common + coord_cartesian(xlim = xlim, ylim = c(0.4, 8.6), expand = FALSE, clip = "off") + labs(x = NULL, y = NULL) +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), plot.margin = margin(4, 8, 0, 4))

# lower panel: y in printed pt around the swarm centre line; the three-line x title is drawn here (clip off)
h <- 33.31 / 2; ty <- -(h + c(14.0, 21.4, 28.8))   # line tops (baselines 18.7/26.1/33.5 pt below the axis, minus the 4.7 pt cap height)
p_bot <- ggplot(swarm, aes(enrichment, y_offset_pt)) +
  geom_point(size = 2.721 * k / .pt / 0.75, stroke = 0, shape = 16, colour = "#A7A9AB") +
  annotate("text", x = mean(xlim), y = ty[1], label = "Enrichment", size = txt(S$axis_title), family = sce2g_font, vjust = 1) +
  annotate("text", x = mean(xlim), y = ty[2], label = "GWAS variants versus common variants", size = txt(S$axis_title), family = sce2g_font, vjust = 1) +
  annotate("text", x = mean(xlim), y = ty[3], label = "italic(n)*' = 210 glucose variants'", parse = TRUE, size = txt(S$axis_title), family = sce2g_font, vjust = 1) +
  scale_y_continuous(breaks = 0, labels = "Other biosamples") +
  common + coord_cartesian(xlim = xlim, ylim = c(-h, h), expand = FALSE, clip = "off") + labs(x = NULL, y = NULL) +
  theme(axis.ticks.y = element_line(), plot.margin = margin(0, 8, 40, 4))

fig <- (p_top / p_bot) + plot_layout(widths = unit(115.315 * k, "pt"), heights = unit(c(50.677, 33.307) * k, "pt"))
save_panel(fig, "figures/fig5d_pointrange.pdf", width_mm = 69, height_mm = 48)
save_panel(fig, "figures/fig5d_pointrange.png", width_mm = 69, height_mm = 48)
cat("wrote figures/fig5d_pointrange.{pdf,png}\n")
