# Fig 2b — precision-recall curves on the CRISPR training benchmark (K562, n = 461 positives).
#
# Anatomy, as in the paper: square 0-1 axes with quarter breaks; one curve per predictor in its
# fixed palette colour (0.42 pt); binary predictors (16-18) as 3.7 pt discs; black dashed line at the
# operating recall of 0.7 and grey dashed line at chance precision; curves labelled by index number
# in their own colour (4.5 pt), tied to the curve by a dashed leader in the same colour. Labels are placed
# by hand in three groups (right of the crossing at recall 0.7; mid-plot for the weaker curves;
# bottom-left for the four weakest), exactly as the published panel does.
#
# Curves here are computed from SIMULATED scores (positives ~ N(mu, 1) vs negatives ~ N(0, 1), plus
# a few confidently wrong negatives for the weaker methods) so they carry real PR jaggedness. Swap
# `scores()` for real predictions to plot data.
#
#   Rscript fig2b_pr.R   # -> fig2b_pr.{pdf,png}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })
source(Sys.getenv("SCE2G_THEME", "../../assets/theme_sce2g.R"))
set.seed(2)
S <- sce2g$type$size_pt; N <- sce2g$color$neutral; M <- sce2g$mark_pt
txt <- function(size = S$tick_label) size / .pt
n_pos <- 461; n_neg <- 9225; chance <- n_pos / (n_pos + n_neg)

# ---- simulated predictors: separation (mu) and number of confidently wrong negatives ----------
# Positives ~ N(mu, 1), negatives ~ N(0, 1). mu is calibrated so each curve crosses the paper's
# precision at its label anchor (recall 0.7 strong / 0.25 mid / 0.15 weak). `bad` confidently wrong
# negatives shape the low-recall end: for the strong methods they are spread through the top of the
# ranking (N(mu + 1.2, 0.8)) and pull precision to ~0.9 by recall 0.1 as in the paper; for the weak
# methods they sit at the very top (N(mu + 3.5, 1)) and give the characteristic rise from zero.
pars <- tribble(~id, ~mu,  ~bad, ~spread,
                "1",  2.45, 12, TRUE,  "2",  2.38, 12, TRUE,  "3",  2.30, 14, TRUE,  "5",  2.25, 12, TRUE,
                "6",  2.15, 14, TRUE,  "15", 1.96, 16, TRUE,
                "4",  1.44, 12, FALSE, "7",  1.32, 3, FALSE,  "11", 1.16, 4, FALSE,  "9",  1.11, 3, FALSE,  "8",  0.97, 3, FALSE,
                "10", 0.72, 14, FALSE, "12", 0.70, 16, FALSE, "13", 0.68, 18, FALSE, "14", 0.66, 20, FALSE)
pr_curve <- function(mu, bad, spread) {
  s <- c(rnorm(n_pos, mu), rnorm(n_neg), if (spread) rnorm(bad, mu + 1.2, 0.8) else rnorm(bad, mu + 3.5))
  y <- c(rep(1, n_pos), rep(0, n_neg + bad))
  o <- order(-s); tp <- cumsum(y[o]); k <- seq_along(o)
  tibble(recall = tp / n_pos, precision = tp / k) %>% distinct(recall, .keep_all = TRUE) %>% filter(recall > 0)
}
pr <- pars %>% rowwise() %>% reframe(id, pr_curve(mu, bad, spread)) %>% mutate(id = factor(id, pars$id))
at_recall <- function(id_, r) pr %>% filter(id == id_, recall >= r) %>% slice(1) %>% transmute(x = recall, y = precision)

# binary predictors: single operating points
pts <- tibble(id = c("16", "17", "18"), x = c(.555, .52, .83), y = c(.468, .488, .215),
              lx = c(.575, .49, .80), ly = c(.44, .49, .19), hj = c(0, 1, 1))

# ---- label plan: (id, anchor recall, label x, label y, hjust) in the paper's three groups --------
lab <- tribble(~id,  ~r,   ~lx,  ~ly,   ~hj,
               "1",  .70, .93, .600, 0,   "2",  .70, .93, .555, 0,   "3",  .70, .93, .510, 0,
               "5",  .70, .93, .465, 0,   "6",  .70, .93, .420, 0,   "15", .70, .93, .375, 0,
               "4",  .25, .57, .360, 0,   "7",  .25, .57, .320, 0,   "11", .25, .57, .280, 0,
               "9",  .26, .57, .245, 0,   "8",  .27, .57, .205, 0,
               "10", .15, .01, .120, 0,   "12", .15, .01, .080, 0,   "13", .15, .01, .040, 0,   "14", .15, .01, .000, 0) %>%
  rowwise() %>% mutate(a = list(at_recall(id, r))) %>% unnest(a) %>%
  mutate(id = factor(id, pars$id),
         # leader ends just short of the text: left edge for right-side/mid groups, right edge for bottom-left
         ex = ifelse(lx > .5, lx - .012, lx + .045), ey = ly)

p <- ggplot(pr, aes(recall, precision, colour = id)) +
  # reference lines: 0.5 pt, 2 pt on / 2 pt off in the paper; R dash units are 0.75 pt (lwd < 1 is clamped), so "33" = 2.25
  geom_hline(yintercept = chance, linetype = "33", colour = N$chance_line, linewidth = pt_mm(.5)) +
  geom_vline(xintercept = .7, linetype = "33", linewidth = pt_mm(.5)) +
  geom_path(linewidth = pt_mm(.42)) +
  geom_segment(data = lab, aes(x = x, y = y, xend = ex, yend = ey, colour = id), linetype = "31", linewidth = pt_mm(.42)) +   # leaders: paper 2.52 on / 0.84 off; "31" = 2.25 / 0.75
  geom_text(data = lab, aes(lx, ly, label = id, colour = id, hjust = hj), size = txt(S$annotation), family = sce2g_font) +   # 4.5 pt, measured
  geom_point(data = pts, aes(x, y, colour = id), size = M$threshold_dot_diameter / .pt / 0.75, shape = 16, stroke = 0) +   # 3.7 pt disc (3.3 fill + 0.42 outline in the paper)
  geom_text(data = pts, aes(lx, ly, label = id, hjust = hj), colour = "black", size = txt(S$annotation), family = sce2g_font) +
  scale_colour_predictor(guide = "none") +
  scale_x_fraction(labels = c("0", "0.25", "0.50", "0.75", "1.00")) +
  scale_y_fraction(labels = c("0", "0.25", "0.50", "0.75", "1.00")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(x = expression("Recall (" * italic(n) * " = 461)"), y = "Precision", title = "CRISPR benchmark (training)") +
  theme_sce2g() + theme(aspect.ratio = 1)

save_panel(p, "fig2b_pr.pdf", width_mm = 62, height_mm = 60)
save_panel(p, "fig2b_pr.png", width_mm = 62, height_mm = 60)
cat("wrote fig2b_pr.{pdf,png}\n")
