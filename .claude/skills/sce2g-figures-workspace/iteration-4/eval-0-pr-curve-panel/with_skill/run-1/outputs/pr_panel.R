# Precision-recall panel, scE2G Fig. 2b style: four predictors (paper ids 1, 2, 5, 15) from pr_curves.csv.
# Positives: the smallest recall step is 1/400, so n = 400 (derived from the curve, not supplied).
#   Rscript pr_panel.R   # -> pr_panel.{svg,pdf,png}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
source("R/theme_sce2g.R")
S <- sce2g$type$size_pt; N <- sce2g$color$neutral
txt <- function(size = S$tick_label) size / .pt

pr <- read.csv("pr_curves.csv") %>% mutate(id = factor(predictor_id, c(1, 2, 5, 15))) %>% arrange(id, recall)
n_pos <- round(1 / min(diff(sort(unique(pr$recall)))))   # 400
chance <- 0.06; op_recall <- 0.7

# labels right of the 0.7 crossing, stacked by precision there; dashed leader from the crossing to the label
lab <- pr %>% filter(recall >= op_recall) %>% group_by(id) %>% slice(1) %>% ungroup() %>%
  arrange(desc(precision)) %>%
  mutate(lx = .93, ly = c(.62, .575, .53, .38), ex = lx - .012, ey = ly)

p <- ggplot(pr, aes(recall, precision, colour = id)) +
  geom_hline(yintercept = chance, linetype = "33", colour = N$chance_line, linewidth = pt_mm(.5)) +
  geom_vline(xintercept = op_recall, linetype = "33", linewidth = pt_mm(.5)) +
  geom_path(linewidth = pt_mm(.42)) +
  geom_segment(data = lab, aes(recall, precision, xend = ex, yend = ey, colour = id), linetype = "31", linewidth = pt_mm(.42)) +
  geom_text(data = lab, aes(lx, ly, label = id, colour = id), hjust = 0, size = txt(S$annotation), family = sce2g_font) +
  scale_colour_predictor(guide = "none") +
  scale_x_fraction(labels = c("0", "0.25", "0.50", "0.75", "1.00")) +
  scale_y_fraction(labels = c("0", "0.25", "0.50", "0.75", "1.00")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(x = bquote("Recall (" * italic(n) * " = " * .(n_pos) * ")"), y = "Precision", title = "CRISPR benchmark (training)") +
  theme_sce2g() + theme(aspect.ratio = 1)

save_panel(p, "pr_panel", width_mm = 62, height_mm = 60)
cat("n_pos =", n_pos, "; wrote pr_panel.{svg,pdf,png}\n")
