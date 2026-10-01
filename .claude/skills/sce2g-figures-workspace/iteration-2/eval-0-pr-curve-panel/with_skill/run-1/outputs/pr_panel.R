# PR-curve panel in the scE2G Fig. 2b style, from pr_curves.csv (predictor_id, recall, precision).
# Adapted from the reference implementation figures/fig2b_pr.R in the sce2g-design repo.
#   Rscript pr_panel.R   # -> pr_panel.{pdf,png}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
source("R/theme_sce2g.R")
S <- sce2g$type$size_pt; N <- sce2g$color$neutral
txt <- function(size = S$tick_label) size / .pt
chance <- 0.06; op_recall <- 0.7

pr <- read.csv("pr_curves.csv") %>% mutate(id = factor(predictor_id, names(sce2g_predictor_cols))) %>% arrange(id, recall)

# labels right of the operating-recall crossing, ordered by precision there; dashed leader from anchor to text
lab <- pr %>% filter(recall >= op_recall) %>% group_by(id) %>% slice(1) %>% ungroup() %>%
  arrange(desc(precision)) %>%
  mutate(x = recall, y = precision, lx = .93, ly = c(.60, .555, .51, .40), ex = lx - .012, ey = ly)

p <- ggplot(pr, aes(recall, precision, colour = id)) +
  geom_hline(yintercept = chance, linetype = "33", colour = N$chance_line, linewidth = pt_mm(.5)) +   # "33" = 2.25/2.25 pt (paper 2/2)
  geom_vline(xintercept = op_recall, linetype = "33", linewidth = pt_mm(.5)) +
  geom_path(linewidth = pt_mm(.42)) +
  geom_segment(data = lab, aes(x = x, y = y, xend = ex, yend = ey, colour = id), linetype = "31", linewidth = pt_mm(.42)) +
  geom_text(data = lab, aes(lx, ly, label = id, colour = id), hjust = 0, size = txt(S$annotation), family = sce2g_font) +
  scale_colour_predictor(guide = "none") +
  scale_x_fraction(labels = c("0", "0.25", "0.50", "0.75", "1.00")) +
  scale_y_fraction(labels = c("0", "0.25", "0.50", "0.75", "1.00")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(x = "Recall", y = "Precision", title = "CRISPR benchmark (training)") +
  theme_sce2g() + theme(aspect.ratio = 1)

save_panel(p, "pr_panel.pdf", width_mm = 62, height_mm = 60)
save_panel(p, "pr_panel.png", width_mm = 62, height_mm = 60)
cat("wrote pr_panel.{pdf,png}\n")
