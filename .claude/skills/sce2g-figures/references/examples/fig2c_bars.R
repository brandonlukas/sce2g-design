# Fig 2c — AUPRC on the CRISPR training benchmark by enhancer-TSS distance bin.
#
# Values, 95 % CIs, bar order and group sizes were digitised from the published vector panel
# (figures/data/fig2c_values.csv). Layout facts measured there and reproduced here:
#   * bars are 3.4 pt wide and TOUCH within a group; groups are separated by 1.7 bar widths;
#   * the order (1 2 3 5 6 15 4 7 11 8 10 9 12 13 14) is the "All" ranking and is kept in every bin;
#   * index numbers sit between the bar baseline and the x-axis line, which runs 5.7 pt below zero
#     and carries one 2.5 pt tick per group at the group centre;
#   * a three-row table below the axis: bin label, CRISPR positives, CRISPR negatives, with
#     right-aligned row headers (6.5 pt) ending under the left of the first group;
#   * y axis 0-0.75 in quarters, labelled "0", "0.25", "0.50", "0.75"; error bar caps 0.9 pt.
#
#   Rscript fig2c_bars.R   # -> fig2c_bars.{pdf,png}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(readr) })
source(Sys.getenv("SCE2G_THEME", "../../assets/theme_sce2g.R"))
S <- sce2g$type$size_pt; N <- sce2g$color$neutral
txt <- function(size = S$tick_label) size / .pt

d <- read_csv("data/fig2c_values.csv", col_types = cols(id = col_character(), group = col_character()))
order_ids <- d %>% filter(group == "All") %>% arrange(desc(auprc)) %>% pull(id)   # ranking in the All bin
bins <- c("All", "[0,20)", "[20,100)", "[100,1000)")
slot <- 1; gap <- 1.7                                                                 # bar width 1; group gap 1.7 bars
d <- d %>% mutate(gi = match(group, bins), rank = match(id, order_ids),
                  pos = (gi - 1) * (length(order_ids) + gap) + rank - 0.5)
grp <- d %>% group_by(group, gi, positives, negatives) %>% summarise(xc = mean(pos), .groups = "drop")

# vertical anchors in data units, from the measured layout (85.9 pt per 1.0 AUPRC)
y_axis <- -5.7 / 85.9; y_idx <- -2.9 / 85.9
y_rows <- -c(15.0, 22.8, 30.6) / 85.9
x_hdr  <- (131.6 - 117.1) / 3.4          # header right edge, in bar widths right of the first bar's left edge

p <- ggplot(d, aes(pos, auprc, fill = id)) +
  geom_col(width = slot) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.9 / 3.4, linewidth = pt_mm(.25)) +
  geom_text(aes(y = y_idx, label = id), size = 2.8 / .pt, family = sce2g_font) +          # paper: 2.5 pt (see note)
  annotate("segment", x = -0.6, xend = max(d$pos) + 0.6, y = y_axis, yend = y_axis, linewidth = pt_mm(.25)) +
  annotate("segment", x = grp$xc, xend = grp$xc, y = y_axis, yend = y_axis - 2.5 / 85.9, linewidth = pt_mm(.25)) +   # one tick per group
  # three-row table
  annotate("text", x = grp$xc, y = y_rows[1], label = grp$group, size = txt(), family = sce2g_font) +
  annotate("text", x = grp$xc, y = y_rows[2], label = grp$positives, size = txt(), family = sce2g_font) +
  annotate("text", x = grp$xc, y = y_rows[3], label = scales::comma(grp$negatives), size = txt(), family = sce2g_font) +
  annotate("text", x = x_hdr, y = y_rows, label = c("Distance to TSS (kb):", "CRISPR positives:", "CRISPR negatives:"),
           hjust = 1, size = txt(S$axis_title), family = sce2g_font) +
  scale_fill_predictor(guide = "none") +
  scale_x_continuous(expand = expansion(add = c(0.6, 0.6))) +
  scale_y_continuous(breaks = seq(0, .75, .25), labels = c("0", "0.25", "0.50", "0.75")) +
  coord_cartesian(ylim = c(y_axis, 0.885), expand = FALSE, clip = "off") +
  labs(x = NULL, y = "Performance (AUPRC)", title = "Prediction performance versus distance to TSS\n(CRISPR, training)") +
  theme_sce2g() +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(),
        plot.margin = margin(4, 6, 34, 40))   # room for the row headers left of the y axis

save_panel(p, "fig2c_bars.pdf", width_mm = 99, height_mm = 46)
save_panel(p, "fig2c_bars.png", width_mm = 99, height_mm = 46)
cat("wrote fig2c_bars.{pdf,png}\n")
