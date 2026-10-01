# Component recipes

Each recipe: what the paper's panel looks like, the measured sizes, a working
ggplot2 (or svgkit) core, and the full reference implementation bundled in
`references/examples/<name>.R|py` (self-contained and runnable — see its
README — with digitised data in `references/examples/data/`).
All snippets assume `source("R/theme_sce2g.R")` and
`S <- sce2g$type$size_pt; txt <- function(pt = S$tick_label) pt / .pt`.

Contents: 1 PR curve · 2 grouped bars · 3 ranked horizontal bars · 4 box plot + filled bars ·
5 scatter small multiples + shared legend · 6 violin with point-interval · 7 point-range + beeswarm ·
8 symmetric heatmap · 9 locus plot · 10 motif logo · 11 UMAP triptych · 12 schematic (svgkit) ·
13 matplotlib fallback · 14 finishing a figure

## 1. Precision–recall / enrichment–recall curves (Fig. 2b) — `fig2b_pr.R`

Square 40 × 40 mm panel, 0–1 axes in quarters, one 0.42 pt curve per predictor
in its palette colour, binary predictors as 3.7 pt discs, black dashed vertical
at the operating recall (0.7), grey dashed horizontal at chance, curves
labelled at their ends with the index number in 4.5 pt with a dashed leader.

```r
p <- ggplot(pr, aes(recall, precision, colour = id)) +
  geom_hline(yintercept = chance, linetype = "33", colour = N$chance_line, linewidth = pt_mm(.5)) +
  geom_vline(xintercept = .7, linetype = "33", linewidth = pt_mm(.5)) +
  geom_path(linewidth = pt_mm(.42)) +
  geom_segment(data = lab, aes(x, y, xend = ex, yend = ey, colour = id), linetype = "31", linewidth = pt_mm(.42)) +
  geom_text(data = lab, aes(lx, ly, label = id, colour = id, hjust = hj), size = txt(S$annotation), family = sce2g_font) +
  geom_point(data = pts, aes(x, y, colour = id), size = 3.74 / .pt / 0.75, shape = 16, stroke = 0) +
  scale_colour_predictor(guide = "none") +
  scale_x_fraction(labels = c("0", "0.25", "0.50", "0.75", "1.00")) + scale_y_fraction(labels = c("0", "0.25", "0.50", "0.75", "1.00")) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(x = expression("Recall (" * italic(n) * " = 461)"), y = "Precision", title = "CRISPR benchmark (training)") +
  theme_sce2g() + theme(aspect.ratio = 1)
save_panel(p, "figures/pr_panel", width_mm = 62, height_mm = 60)   # writes .svg, .pdf and 900 dpi .png
```
Label placement is by hand in three groups (strong curves right of the 0.7
crossing, mid curves in the open area, weakest bottom-left); `lab` holds
(id, anchor recall, label x/y, hjust). The recall title always carries the
number of positives, "Recall (*n* = 461)": if *n* is not supplied, derive it
from the curve (the smallest non-zero recall step is 1/*n*) and say so, or
ask — do not drop it.

## 2. Grouped performance bars (Fig. 2c) — `fig2c_bars.R`, data `fig2c_values.csv`

Bars touch within a group (width = slot), groups 1.7 slots apart, order fixed
by the reference group's ranking; 0.25 pt error bars with 0.9 pt caps; index
numbers in a strip between the bars and an axis line 5.7 pt below zero, no x
ticks; a three-row table under the axis (group label, positives, negatives)
with right-aligned 6.5 pt row headers; y 0–0.75 in quarters.

```r
d <- d %>% mutate(pos = (gi - 1) * (n_ids + 1.7) + rank - 0.5)
p <- ggplot(d, aes(pos, auprc, fill = id)) +
  geom_col(width = 1) +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.9 / 3.4, linewidth = pt_mm(.25)) +
  geom_text(aes(y = -2.9 / 85.9, label = id), size = 2.8 / .pt, family = sce2g_font) +     # paper 2.5 pt; 2.8 for legibility
  annotate("segment", x = -0.6, xend = max(d$pos) + 0.6, y = -5.7 / 85.9, yend = -5.7 / 85.9, linewidth = pt_mm(.25)) +
  annotate("text", x = grp$xc, y = y_rows[1], label = grp$group, size = txt(), family = sce2g_font) +
  annotate("text", x = x_hdr, y = y_rows, label = c("Distance to TSS (kb):", "CRISPR positives:", "CRISPR negatives:"), hjust = 1, size = txt(S$axis_title), family = sce2g_font) +
  scale_fill_predictor(guide = "none") +
  scale_y_continuous(breaks = seq(0, .75, .25), labels = c("0", "0.25", "0.50", "0.75")) +
  coord_cartesian(ylim = c(-5.7 / 85.9, 0.885), expand = FALSE, clip = "off") +
  theme_sce2g() + theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank(), plot.margin = margin(4, 6, 34, 40))
```

## 3. Ranked horizontal bars (Fig. 3a forward selection; Fig. 5g PoPS) — `fig3a_features.R`, `fig5g_pops.R`

Categories on y, highest at top, bars 0.9 of the row; fill `#445369`
(emphasis) or `#6A768B` with `#B8C0CC` for de-emphasised rows; x keeps the
default 5 % padding on both ends (bars start 0.05 units right of the axis).
Fig. 3a adds a 0.32 pt CI with 2 pt caps, a 2.3 pt black dot for the
cumulative metric, 7.5 pt stars rotated 90° beside the CI (red `#E82829` for
negative deltas), and a right-hand *P* column headed *P*~bootstrap~.

```r
p <- ggplot(d, aes(pops, gene, fill = highlight)) +
  geom_col(width = sce2g$mark_pt$bar_width_fraction) +
  scale_fill_manual(values = c(`TRUE` = "#6A768B", `FALSE` = "#B8C0CC"), guide = "none") +
  scale_x_continuous(breaks = seq(0, 1, .25), labels = c("0", "0.25", "0.50", "0.75", "1.00"), limits = c(0, 1)) +
  labs(x = "PoPS", y = "Gene") + theme_sce2g()
fig <- p + plot_layout(widths = unit(36.6, "mm"), heights = unit(15.5, "mm"))
```

## 4. Box plot + filled/stacked bars (Fig. 3g) — `fig3g_boxbar.R`

Left: boxes 0.8 of the slot, fill `#D3D3D3`, **black** 0.64 pt outline
(ggplot's default is grey20), 1 pt median, dashed whiskers `"34"` with round
caps, staples half the box width, no outliers; two label rows under the axis
with 6.5 pt right-aligned row headers; the *P* value right-aligned at the top.
Right: bars 0.66 of the slot with a 0.25 pt outline, a white outlined
"remainder" up to exactly 100 %, y 90–100 % with **no** expansion, counts in
4.5 pt. Panels 17.4 × 35.6 and 24.1 × 34.5 mm.

```r
geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax), stat = "identity",
             fill = N$box_fill, colour = "black", width = 0.8, linewidth = pt_mm(.639), median.linewidth = pt_mm(1),
             whisker.linetype = "34", staplewidth = 0.5, lineend = "round")
# filled bar from the 90 % floor plus its white remainder
geom_rect(aes(xmin = i - .33, xmax = i + .33, ymin = 90, ymax = q, fill = fill), colour = "black", linewidth = pt_mm(.25)) +
geom_rect(aes(xmin = i - .33, xmax = i + .33, ymin = q, ymax = 100), fill = "white", colour = "black", linewidth = pt_mm(.25)) +
coord_cartesian(xlim = c(0.5, 2.5), ylim = c(90, 100), expand = FALSE, clip = "off")
```

## 5. Scatter small multiples + shared legend (Fig. 4a) — `fig4a_depth.R`, `fig4_legend.R`

Three 39.5 mm square panels; 4.3 pt discs at alpha 0.7 in `#6B2063`; 0.5 pt
2/2-dashed references — black vertical at the data threshold, ABC `#008F99`
and distance `#6A768B` horizontals; y from 0 with default padding; x breaks at
integers with a log10 title. The shared "Model" key under the row shows each
reference as a 12 pt dashed line and the model as a point.

```r
panel <- function(m) ggplot(d[d$metric == m, ], aes(log10_fragments, value)) +
  geom_hline(yintercept = r$dist, colour = "#6A768B", linewidth = pt_mm(.5), linetype = "33") +
  geom_hline(yintercept = r$abc,  colour = "#008F99", linewidth = pt_mm(.5), linetype = "33") +
  geom_vline(xintercept = log10(2e6), linewidth = pt_mm(.5), linetype = "33") +
  geom_point(colour = "#6B2063", alpha = 0.7, size = 4.3 / .pt / 0.75, stroke = 0, shape = 16) +
  scale_x_continuous(breaks = 5:9, limits = c(4.75, 8.84)) +
  scale_y_continuous(breaks = seq(0, .6, .2), labels = c("0", "0.2", "0.4", "0.6"), limits = c(0, r$top)) +
  labs(x = expression(log[10] * "(number of unique ATAC fragments)"), y = r$title) + theme_sce2g()
fig <- (panel("auprc") | panel("precision") | panel("recall")) + plot_layout(widths = unit(rep(39.5, 3), "mm"), heights = unit(39.5, "mm"))
```

## 6. Violin with point-interval and brackets (Fig. 5b) — `fig5b_violin.R`

Panel 34.7 × 50.1 mm; violins filled in the category colour with no outline
(`#8B98AC` unmatched, `#9A3983` same lineage, `#260032` same cell type);
inner black point-interval: 3.05 pt median dot, 1.79 pt line for the central
66 %, 0.77 pt for the central 95 % (`stat_pointinterval_sce2g()` on raw
values); brackets 0.25 pt with 1.1 pt down-ticks, 8.5 pt stars just below the
bar, *P* in 4.5 pt above with a plotmath exponent; category labels as rows
under the axis, filled from the top: the label's line(s) first, then italic
*n* on the **very next** row (so "Unmatched / *n* = 2,259" occupies rows 1–2
while "Same / lineage / *n* = 204" occupies rows 1–3 — the *n* rows are not
aligned across categories, exactly as in Fig. 5b). The paper stacks brackets
9.6 pt apart; with Helvetica the exponent of the lower *P* touches the upper
bar, so space the two rows 12–15 pt apart (whatever clears the exponent).
Keep the paper's two-row layout — adjacent pairs (1–2, 2–3) on the lower row,
the outer pair (1–3) above — rather than one row per bracket; the two *P*
labels on the lower row fit when each is aligned *away* from the shared
middle category (`hjust = 1` at the left bracket's right end, `hjust = 0` at
the right bracket's left end).

```r
p <- ggplot(d, aes(category, value)) +
  geom_violin(aes(fill = category), colour = NA, width = 0.9) + scale_fill_manual(values = cols, guide = "none") +
  stat_pointinterval_sce2g() +
  sig_bracket(1, 2, 14.5, "***", p = "italic(P)*' = 4.39 × '*10^'–41'") +
  annotate("text", x = 1, y = -1.9, label = "Unmatched", size = txt(), family = sce2g_font, vjust = 0) +
  annotate("text", x = 1, y = -2.65, label = "italic(n)*' = 2,259'", parse = TRUE, size = txt(), family = sce2g_font, vjust = 0) +
  scale_y_continuous(breaks = seq(0, 15, 5)) +
  coord_cartesian(xlim = c(0.4, 3.6), ylim = c(-0.78, 16.3), expand = FALSE, clip = "off") +
  labs(x = NULL, y = "Enrichment\neQTLs versus common variants") +
  theme_sce2g() + theme(axis.text.x = element_blank(), plot.margin = margin(6, 4, 38, 4))
```

## 7. Point-range with a beeswarm reference row (Fig. 5d) — `fig5d_pointrange.R`

Two stacked panels sharing one y axis line: rows at a 6.5 pt pitch with 6.2 pt
discs and 0.42 pt CI lines with 1.3 pt caps in the category colour
(`#260032` significant, `#9A3983` others), 8.5 pt stars right of the CI; below,
a 12.4 mm panel with the other biosamples as a grey `#A7A9AB` 2.9 pt
**swarm** (touching packing — `ggbeeswarm::geom_beeswarm(method = "swarm")`,
tune `cex` so dots touch), the x axis only under the swarm, and a three-line
x title with italic *n*.

```r
p_top <- ggplot(d, aes(y = cell_type, colour = colour)) +
  geom_segment(aes(x = lo, xend = hi, yend = cell_type), linewidth = pt_mm(0.42)) +
  geom_point(aes(x = enrichment), size = 6.2 / .pt / 0.75, stroke = 0, shape = 16) +
  scale_colour_identity() + theme_sce2g() +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank())
p_bot <- ggplot(other, aes(enrichment, 0)) +
  ggbeeswarm::geom_beeswarm(method = "swarm", cex = 5.9, size = 2.9 / .pt / 0.75, stroke = 0, shape = 16, colour = "#A7A9AB") +
  scale_y_continuous(breaks = 0, labels = "Other biosamples") + theme_sce2g()
fig <- (p_top / p_bot) + plot_layout(widths = unit(42.8, "mm"), heights = unit(c(18.8, 12.4), "mm"))
```

## 8. Symmetric split heatmap (Fig. 5a) — `fig5a_heatmap.R`, data `fig5a_*.csv`

3.23 pt square cells, no visible borders (give tiles a hairline outline in
their own colour to hide viewer seams); two metrics split along the divider
(red Jaccard upper left, blue Pearson lower right), divider cells navy
`#001636`; a 2.9 pt dataset side bar, two 7.8 × 38 pt colour bars stacked at
the right with 1.5 pt white edge ticks; 2.85 pt row labels (1.9 pt
superscripts), columns rotated 60°. Draw in point coordinates with
`theme_void()`; map values to colour through the sampled ramps
(`scale_fill_sce2g_seq("red")` / `"blue"` for new data).

```r
ggplot() +
  geom_tile(data = cells, aes(x, y, fill = fill, colour = fill), width = sz, height = sz, linewidth = pt_mm(.2)) +
  geom_text(data = labels, aes(x, -y, label = text, size = size / .pt, angle = angle), hjust = 0, vjust = 0, family = sce2g_font) +
  scale_fill_identity() + scale_colour_identity() + scale_size_identity() +
  coord_cartesian(xlim = c(58, 292), ylim = c(-262, -62), expand = FALSE, clip = "off") + theme_void()
```

## 9. Locus plot (Fig. 1b, 3e) — `fig1b_locus.R`, `fig3e_locus.R`

Top to bottom: 0.5 pt scale bar with round caps and 0.4 × 2.6 pt end ticks →
gene models (black exon blocks, `#C9C9C9` intron line with a strand chevron,
italic 6.5 pt names) → per-track arc rows on a 15.6 pt pitch → ATAC track
(`#006EAE` fill, no y axis) → tested elements (`#DC6464`). Arcs are
**constant-height quadratics** (0.95 of the row pitch) from a 0.8 pt TSS tick
to a 0.8 pt element foot, colour and width ∝ score on the purple ramp
(`arc_paths()` / `geom_arc_sce2g()`); track names 5.5 pt right-aligned in a
label column. Heatmap variant: 4.2 pt row labels minimum, `#DEDFE1` band on
rows of interest, red ramp, horizontal colour bar.

```r
geom_arc_sce2g(arcs %>% select(x, xend, score), height = 0.95 * pitch)   # x = element, xend = TSS
geom_tile(data = arcs, aes(x, baseline - foot/2, width = el_w, height = foot, fill = score))  # element feet
geom_ribbon(data = atac, aes(x, ymin = base, ymax = base + h), fill = G$atac_track)
```

## 10. Motif logo (Fig. 5h) — `fig5h_logo.R`

ggseqlogo letters (A `#00904A`, C `#1A5F9D`, G `#F1AD30`, T `#D92641`) on a
37.2 × 17.7 mm panel, bits 0–2 with no expansion, 4.5 pt tick labels; the
reference/alternate sequences beneath in 8 pt Courier `#8B98AC` with the
variant base black / `#9A1D20`, 4 pt scores, and the rsID under a notched
`annotate_arrow_sce2g()`.

```r
ggplot() + ggseqlogo::geom_logo(m, method = "custom", seq_type = "dna", col_scheme = cols) +
  geom_text(data = seq_df, aes(x, y, label = letter, colour = colour), size = txt(8), family = "Courier") +
  annotate_arrow_sce2g(11.03, -0.93, length_pt = 5, head_len = 3, stem = 0.5) +
  scale_y_continuous(breaks = seq(0, 2, .5), labels = c("0", "0.5", "1.0", "1.5", "2.0")) +
  coord_cartesian(xlim = c(0.55, 18.51), ylim = c(0, 2), expand = FALSE, clip = "off") + labs(y = "Bits", title = "TP73 motif") +
  theme_sce2g() + theme(axis.text = element_text(size = S$annotation))
```
ggseqlogo installs its own x scale and coord; yours replace them (two notices are expected).

## 11. UMAP triptych (Fig. 3f) — `fig3f_umap.R`

Three equal square panels; ~1 pt points (white hairline outline only in the
categorical panel); axes replaced by a 15 pt L of notched arrows at the
bottom-left of the first panel with "UMAP1"/"UMAP2" in 5.5 pt
(`theme_sce2g() + theme_embedding()`); legends inside the panel at the lower
right: a 4 × 18 pt vertical bar with white ticks, title rotated 90° to its
left. Expression `#D3D3D3 → #0000FF`, binary state `#D3D3D3 → #FF0000`.

## 12. Schematic (Fig. 1a, 2a, 3b) — `fig1a_schematic.py`, `fig2a_table.py`, `fig3b_strategies.py`

```python
from svgkit import Canvas
c = Canvas(ox=100, oy=50, w=350, h=240)                     # page-point frame; emits mm SVG
c.add(c.text(155, 72, ["(1) Prepare single-cell", "multiome data"], 6.5, "middle"))
c.add(c.rect(230, 60, 330, 120, "#DEDFE1"))                 # container for a data table / track / feature box
c.add(c.arrow_poly(250, 310, 150, "#54AB58", stem=1.85))    # gene arrow (horizontal, filled)
# connectors LAST, tips exactly on the target edge (x = 230 here), never inside a container:
c.add(*c.connector([(200, 90), (230, 90)]))                  # straight, 0.93 pt, notched head
c.add(*c.connector([(280, 120), (280, 140), (350, 140)]))    # one 90-degree elbow
# cartoon embedding: on WHITE, never inside a grey container (Fig. 1a, 3b)
for x, y in cells: c.add(c.circle(x, y, 3.19, "#1F84FF"))  # 3.2 pt cells; a meta-cell = grey blob drawn first:
c.add(c.path([("M", 120, 100), ("C", 118, 90, 140, 88, 142, 98), ("C", 144, 108, 122, 110, 120, 100)], stroke="none", fill="#D3D3D3"))
c.save("figures/schematic.svg")                              # + .pdf/.png via Inkscape
```
Keep to three stroke weights (0.25 / 0.9 / 1.85 pt), grey containers only
around data-like content (tables, tracks, feature boxes — the UMAP cell cloud
sits on white with a small notched-arrow axis pair), one magenta `#9C4787`
model box, Courier for data snippets, and the notched arrow everywhere.
Draw order matters because the SVG paints in `add()` order: containers, then
their contents, then connectors last. Use `c.connector()` for every
connector — it makes stem and head one seamless mark with the tip on the
target edge; a hand-built line + separate head shows a seam and tends to
overshoot into the box (Fig. 1a's connectors stop at the container edge).

## 13. matplotlib fallback

```python
import sys; sys.path.insert(0, "<skill>/assets"); from matplotlib_rc import apply, PRED, TOK; apply()
fig, ax = plt.subplots(figsize=(40 / 25.4, 40 / 25.4))
ax.plot(recall, precision, color=PRED["1"], lw=TOK["stroke_pt"]["curve"])
ax.axvline(0.7, color="k", lw=0.5, ls=(0, (4, 4)))           # 2/2 pt at 0.5 pt
ax.set_xticks([0, .25, .5, .75, 1]); ax.set_xticklabels(["0", "0.25", "0.50", "0.75", "1.00"])
fig.savefig("fig.pdf", bbox_inches="tight")
```

## 14. Finishing a figure

- Panel letters 9.3 pt bold outside the plot, top-left, lowercase.
- Panels in a row share height; identical axes share one title.
- One legend per figure, numbered entries, in whitespace.
- Run `scripts/measure_svg.py` on the PDF and read the PNG at 300 %.
