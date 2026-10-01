# Fig 5g — PoPS scores for the genes within 1 Mb of rs7696969 (lymphocyte count): horizontal ranked
# bars, highest at the top, the two prioritised genes in the dark slate and the rest de-emphasised.
#
# Measured from the published vector panel (Fig. 5 is printed at 95 %; values here are / 0.95):
#   panel : 36.6 x 15.5 mm; x 0-1 with the default 5 % padding on BOTH ends (the bars start 0.05 units
#           right of the y axis), breaks every 0.25 labelled 0 / 0.25 / 0.50 / 0.75 / 1.00;
#   rows  : default 0.6-row padding, bars 0.9 of the row; fills #6A768B (IL15, INPP4B) / #B8C0CC (others);
#   text  : 5.5 pt tick labels, upright gene names; 6.5 pt titles "Gene" / "PoPS".
#
#   Rscript fig5g_pops.R   # -> fig5g_pops.{pdf,png}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(patchwork) })
source(Sys.getenv("SCE2G_THEME", "../../assets/theme_sce2g.R"))

d <- tibble(gene = c("IL15", "INPP4B", "GAB1", "USP38"), pops = c(0.8525, 0.5431, 0.0909, 0.0319),   # digitised bar lengths
            highlight = c(TRUE, TRUE, FALSE, FALSE)) %>% mutate(gene = factor(gene, rev(gene)))
p <- ggplot(d, aes(pops, gene, fill = highlight)) +
  geom_col(width = sce2g$mark_pt$bar_width_fraction) +
  scale_fill_manual(values = c(`TRUE` = "#6A768B", `FALSE` = "#B8C0CC"), guide = "none") +
  scale_x_continuous(breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.50", "0.75", "1.00"), limits = c(0, 1)) +
  labs(x = "PoPS", y = "Gene") +
  theme_sce2g()
fig <- p + plot_layout(widths = unit(98.493 / 0.95, "pt"), heights = unit(41.635 / 0.95, "pt"))
save_panel(fig, "fig5g_pops.pdf", width_mm = 50, height_mm = 27)
save_panel(fig, "fig5g_pops.png", width_mm = 50, height_mm = 27)
cat("wrote fig5g_pops.{pdf,png}\n")
