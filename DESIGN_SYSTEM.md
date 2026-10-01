# scE2G Figure Design System

Rules for building panels that sit visually next to the main figures of
Sheth, Qiu et al., *Nature Genetics* 2026 ("Mapping enhancer–gene regulatory
interactions from single-cell data", `refs/sce2g-2026.pdf`).

Every value below was measured from the PDF's vector content (fill/stroke
colours, stroke widths, font names and point sizes) or sampled from 300-dpi
renders where a panel was rasterised (Fig. 1b–d heatmaps, Fig. 3f UMAPs).
Machine-readable copies of the tokens live in `tokens.json`; every component
has a reference implementation under `figures/` that reproduces the paper's
panel from digitised data (`figures/data/`) and was verified against the paper
by re-measuring the rendered PDF.

**Print scale.** The journal prints some figures below design size: Fig. 1 at
93 %, Fig. 5 at 95 %, Figs. 2–4 at 100 %. All values here and in `tokens.json`
are *design* (100 %) values: a 5.22 pt tick label measured on Fig. 5 is the
5.5 pt token. Build at design size and let the layout scale the whole figure.

---

## 1. Foundations

### 1.1 Canvas and layout

| Token | Value | Notes |
|---|---|---|
| Page text width | 522 pt = 184 mm | Full-width figure (Fig. 2, 3, 4, 5) spans this. |
| Single-column width | 249 pt = 88 mm | Half-width panel columns (Fig. 1b vs 1c/d). |
| Panel letter | Graphik Black 9.3 pt, lowercase **a b c**, black | Sits outside the plot, top-left, baseline aligned with the top of the panel's tallest element. |
| Gutter between panels | 24–36 pt horizontal, 28–40 pt vertical | Panels are aligned on a shared left edge per column and share y-axis baselines per row. |
| Panel title | 6.5 pt regular, centered above the plot | Used for the *comparison being shown* ("CRISPR benchmark (held-out)"), not for restating the axis. |
| Aspect ratios | PR curves ~1:1; bar plots 2:1 to 3:1 wide; locus plots 3:1 to 5:1 wide; violin/box 1:1.4 tall | Keep a row of panels equal in height. |
| Background | White. No panel fill, no border, no grid lines. | Only the conceptual-schematic boxes use a filled background (`#DEDFE1`). |

Axes are **open** (left and bottom spines only, no top/right), ticks point
outward, and the plot area has no frame — equivalent to `theme_classic()` in
ggplot or `spines[['top','right']].set_visible(False)` in matplotlib.

**Panel sizes** (plot area, excluding axis text; design scale). Panels in one
row share a height; use these as the starting sizes for each component:

| Component | Panel (w × h, mm) | Example |
|---|---|---|
| PR curve, square | 40 × 40 | Fig. 2b |
| Scatter small multiple | 39.5 × 39.5 | Fig. 4a |
| Grouped bars | 66 × 30 | Fig. 2c |
| Box plot / filled bars pair | 17.4 × 35.6 / 24.1 × 34.5 | Fig. 3g |
| Violin | 34.7 × 50.1 | Fig. 5b |
| Point-range + swarm | 42.8 × (18.8 + 12.4) | Fig. 5d |
| Ranked bars, 4 rows | 36.6 × 15.5 | Fig. 5g |
| Symmetric heatmap, 46 rows | 52 × 52 (3.23 pt cells) | Fig. 5a |
| Motif logo | 37.2 × 17.7 | Fig. 5h |
| UMAP triptych | 3 × 37 × 37 | Fig. 3f |

Set them explicitly (`patchwork::plot_layout(widths = unit(w, "mm"),
heights = unit(h, "mm"))`) rather than letting axis text decide the panel.

### 1.2 Colour

#### 1.2.1 Predictor palette (Figs. 2, 4, ED Figs. 4–6)

Every method has one fixed colour and one fixed index number that never
changes between panels or figures. Methods developed in this paper are
saturated jewel tones; published methods are a desaturated earth ramp;
baselines are greys and navies. This three-tier logic is the core of the
paper's look.

| # | Predictor | Hex | Tier | Contrast vs white |
|---|---|---|---|---|
| 1 | scE2G^Multiome | `#792374` | this study (primary) | 9.1 |
| 2 | scE2G^ATAC | `#006479` | this study | 6.8 |
| 3 | ARC-E2G score | `#B778B3` | this study | 3.3 |
| 4 | Kendall correlation | `#D3A9CE` | this study | 2.0 |
| 5 | ABC (A = scATAC, C = power law) | `#0096A0` | other study (teal) | 3.6 |
| 6 | STARE (gABC, power law) | `#96CED3` | other study (teal) | 1.7 |
| 7 | SnapATAC | `#5E5948` | other study (earth) | 7.0 |
| 8 | Signac | `#916953` | other study (earth) | 4.8 |
| 9 | Cicero | `#888364` | other study (earth) | 3.8 |
| 10 | FigR | `#BC9678` | other study (earth) | 2.7 |
| 11 | SCENT | `#A5A083` | other study (earth) | 2.6 |
| 12 | SCENIC+ | `#C5C1A5` | other study (earth) | 1.8 |
| 13 | DIRECT-NET | `#E0DCCA` | other study (earth) | 1.4 |
| 14 | ArchR | `#DCBCA1` | other study (earth) | 1.8 |
| 15 | Distance to TSS | `#404040` | baseline | 10.4 |
| 16 | Closest TSS | `#6E788D` | baseline | 4.4 |
| 17 | Closest gene | `#435369` | baseline | 7.8 |
| 18 | Within 100 kb of TSS | `#1C2A43` | baseline | 14.4 |
| 19 | In ABC element and distance to TSS | `#96A0B3` | baseline | 2.6 |

Rules:

- Bars are ordered by performance in the **reference group** ("All"), and that
  order is kept in every other group so the same colour sits in the same slot
  across bins (Fig. 2c, 2d); the index number is printed under each bar
  (2.5 pt in the paper; see §3 for the legibility caveat).
- Curves are labelled at their right end with the index number in the curve's
  own colour, connected by a short dashed leader when crowded (Fig. 2b).
- When a panel shows only a subset, keep the original colours and numbers;
  never renumber.

#### 1.2.2 Two-model and category palettes

| Use | Hex | Where |
|---|---|---|
| scE2G^Multiome single-cell | `#792374` | Fig. 3d, 3g |
| scE2G^Multiome meta-cell | `#1E871B` (bar `#00783D`) | Fig. 3d, 3g |
| Kendall single-cell K562 | `#D3A9CE` | Fig. 3c |
| Kendall single-cell K562+PBMC | `#DC6464` | Fig. 3c |
| Kendall meta-cell | `#A8C7A7` | Fig. 3c |
| Retrained ABC+Kendall (meta-cell) | `#228742` | Fig. 3c |
| Same cell type | `#260032` | Fig. 5b–e |
| Same lineage | `#9A3983` | Fig. 5b–e |
| Unmatched / other biosamples | `#8B98AC` (points `#A7A9AB`) | Fig. 5b–e |
| Dataset: cell line / BMMC / PBMC / islets | `#142D4A` `#42566B` `#6A768B` `#8B98AC` | Fig. 5a side bar and key |
| Heatmap divider cells | `#001636` | Fig. 5a anti-diagonal |
| Depth-series scatter point / ABC reference / distance reference | `#6B2063` / `#008F99` / `#6A768B` | Fig. 4 |
| Ranked bars: emphasis / secondary / muted | `#445369` / `#6A768B` / `#B8C0CC` | Fig. 3a, 5g |

Ordinal categories (unmatched → lineage → cell type) use a **single-hue
lightness ramp** of the primary purple, darkest = most specific.

#### 1.2.3 Genomic-track palette (locus plots, Figs. 1b–d, 3e, 5f, ED Fig. 3, 10)

| Element | Hex | Stroke |
|---|---|---|
| Gene models, scale bar | `#000000` (introns `#CED0D2`) | 0.23 pt |
| ATAC-seq signal track | `#006CAD` fill, no outline | — |
| scE2G^Multiome arcs | `#260032` → `#C29CBC` mapped to score | 0.1–0.6 pt, width ∝ score |
| scE2G^ATAC arcs | `#006479` → `#96CED3` | same |
| CRISPR-positive arcs | `#9B241C` | 0.6 pt |
| Tested CRISPR elements | `#DC6464` | bar |
| GWAS credible-set blocks, variant lollipops | `#00909A` | 0.35 pt |
| Highlight variant (dashed vertical) | `#003A4E`, dashed | 0.36 pt |
| Row highlight band | `#DEDFE1` | full-width behind heatmap rows |
| Locus heatmap (score 0.2 → 0.9) | `#DFA3A3` → `#D99296` → `#B05C58` → `#9B241C` | cells have no border |

#### 1.2.4 Continuous scales

| Scale | Stops | Used for |
|---|---|---|
| Sequential red | `#FFFFFF` → `#F2C6C8` → `#D99296` → `#B05C58` → `#9B241C` | scE2G score, Jaccard similarity |
| Sequential blue | `#FFFFFF` → `#B8D5DF` → `#63A8D7` → `#1A5F9D` → `#001636` | Pearson correlation, enrichment (ED Fig. 9d) |
| Sequential purple (Multiome) | `#F3EAF2` → `#B778B3` → `#792374` → `#3E0F3B` | Multiome robustness heatmaps (ED Fig. 7g) |
| Sequential teal (ATAC) | `#E6F2F3` → `#96CED3` → `#0096A0` → `#003A4E` | ATAC robustness heatmaps (ED Fig. 7i) |
| UMAP expression | `lightgrey` `#D3D3D3` → `blue` `#0000FF`, interpolated in CIELAB (Seurat `FeaturePlot` default; mid-tone `#9C76EC`) | Fig. 3f centre |
| UMAP binary state | `lightgrey` `#D3D3D3` → `red` `#FF0000`, CIELAB (mid-tone `#F9896E`) | Fig. 3f right |

Do not use rainbow or viridis; the paper's only rainbow-like scale (ED Fig. 2d
density) is a legacy default and should not be imitated.

#### 1.2.5 Schematic palette (Fig. 1a, 2a, 3b)

| Role | Hex |
|---|---|
| Container / panel background | `#DEDFE1` |
| Secondary container, non-selected cells | `#B9C0CB` |
| Cell nuclei, muted glyphs | `#6A768B` |
| Selected cell cluster | `#7DBCC3` |
| Accessibility peak light / dark | `#85C0E1` / `#006CAD` |
| Gene / transcription arrow | `#54AB58` |
| Model box (scE2G) | `#9C4787` |
| Emphasis (rare) | `#EC9146` |
| Cell-type dots in cartoon UMAP | `#1F84FF` `#F40014` `#00BF2B` |

#### 1.2.6 Neutrals

`#000000` text and axes · `#404040` strong grey · `#6E6E6E` mid grey ·
`#AEAEAE` dashed baseline · `#C9C9C9` / `#CED0D2` hairlines · `#D3D3D3`
box-plot fill · `#DEDFE1` bands · `#E4E4E4` heatmap grid · `#F3F3F3` table
fill.

### 1.3 Typography

The typeset paper uses **Graphik Naturel** (Nature's house sans) for all
figure text and **Harding Text** for captions. Graphik is not freely
licensed; author-supplied figures should use **Helvetica** or **Arial**,
which Nature will substitute 1:1. Set Helvetica at the **nominal** size: its
lowercase runs 10–15 % narrower than Graphik's (measured on Fig. 2a), so
measured line breaks fit without shrinking. Line pitch is 1.2 × the nominal
size (6.6 pt at 5.5 pt) and must not be scaled with any glyph adjustment.

| Role | Size (pt, at final print size) | Weight | Colour |
|---|---|---|---|
| Panel letter | 9.3 | Black / Bold (the only bold text in a figure) | `#000000` |
| Panel title / legend title ("Predictors") | 6.5–8 | Regular (never bold) | black |
| Axis title | 6.5 | Regular | black |
| Tick label, legend entry | 5.5 | Regular | black |
| Secondary annotation (n, P value, dashed-line label) | 4.5–5.5 | Regular; *P*, *n* italic | black or series colour |
| Locus track name | 5.1–6.0 | Regular | black |
| Heatmap row label (many rows) | 4.2 minimum | Regular | black |
| Bar index number | 2.5 in paper — **use 4.5 minimum** | Regular | black |
| Significance stars | 8.5 | Bold | black, `#E82829` for negative direction |
| Sequence / data snippets | Courier 8 (`#8B98AC`, variant base black / `#9A1D20`) | Regular | grey |
| Superscripts (CD4^+, 10^−17, ref. numbers) | 0.7 × the base size, raised ~0.45 em | Regular | as base |
| Caption | Harding 7.0 | Bold lead-in, regular body | black |

Conventions: sentence case everywhere; gene names italic; superscripts for
model variants (scE2G^Multiome, scE2G^ATAC); *n* and *P* italic; units in
parentheses; multi-line axis titles put the metric first and the definition on
the following line(s) in the same size ("Recall / (Fraction of credible sets
linked to target gene) / 159 credible sets from 9 traits"). Numbers use
thin-space or comma thousands separators ("9,225"), and en dashes for ranges
and pairs ("enhancer–gene", "[20,100)").

### 1.4 Strokes and marks

| Element | Value |
|---|---|
| Axis line, tick | 0.25 pt black; ticks 2.4 pt long, outward |
| Reference line (threshold, random baseline) | 0.5 pt, dash 2 pt on / 2 pt off (ggplot `linetype = "33"`, 2.25 pt: R dash units are 0.75 pt for hairlines); black for a chosen threshold, `#AEAEAE` for chance |
| Label leader (curve index numbers) | 0.42 pt in the curve's colour, dash 2.5 pt on / 0.85 pt off (ggplot `linetype = "31"`, 2.25 / 0.75 pt) |
| Data curve (PR, enrichment–recall) | 0.42 pt; highlighted curves up to 0.6–0.8 pt |
| Error bar | 0.25 pt black (0.32 pt in Fig. 3a), cap width ≈ 40 % of bar width |
| Bar outline | none (fill only); stacked bars 0.5 pt black outline |
| Threshold marker on curve | filled circle ⌀ 3.7 pt (3.3 fill + 0.42 outline), series colour |
| Scatter point (many) | ⌀ 4.3 pt, alpha 0.7, no outline (Fig. 4) |
| Point-range dot | ⌀ 6.2 pt filled, 0.42 pt CI line with 1.3 pt caps (Fig. 5d) |
| Violin point-interval | ⌀ 3.05 pt black dot, 1.79 pt (66 %) and 0.77 pt (95 %) lines |
| Cumulative-metric dot beside bars | ⌀ 2.3 pt black (Fig. 3a) |
| Beeswarm point | ⌀ 2.9 pt `#A7A9AB`, swarm-packed (touching) |
| Box plot | 0.64 pt black outline, 1 pt median, whiskers 1.9 pt on / 3.2 pt off with round caps |
| Violin outline | none; fill only |
| Locus arc | 0.1–0.6 pt, alpha 0.4–1, both ∝ score |
| Gene model exon | filled black block, 3.7 pt tall |
| Hairline table rule | 0.35 pt |
| Arrow | one shape everywhere: stem into a swept head with a notched back edge (§2.1); 0.5 pt stem with a 3 pt head for annotations and UMAP axes |

### 1.5 Axis ranges and expansion

The paper follows one rule, even where it looks inconsistent at first glance
(Fig. 4 leaves space below AUPRC = 0; Fig. 3d does not; Fig. 5g, a horizontal bar chart, keeps the default 5 % on both ends of x):

| Mark type | Baseline end | Far end | Example |
|---|---|---|---|
| **Anchored at zero** — bars, columns, stacked bars, areas, histograms, ATAC tracks | **0 expansion**: the first bar sits on the axis line | 5 % headroom (`expansion(mult = c(0, .05))`) | Fig. 2c, 2f, 3d (not 5g, measured) |
| **Unanchored** — points, curves, violins, box plots, point-ranges | default 5 % padding on **both** ends, even for quantities bounded at 0 (AUPRC, recall) | default 5 % | Fig. 2b, 4, 5b–e |
| **Tiles** — heatmaps, value grids, logos | 0 expansion both axes | 0 | Fig. 5a, ED Fig. 7 |
| Categorical axis of a bar chart | default 0.6-unit padding | default | Fig. 3d, 5g |

Additional rules:

- Never clip data with hard limits; if a point sits at 0 or 1 it must remain
  fully visible, which is exactly why unanchored marks keep the 5 % gap.
- Fraction scales (precision, recall, AUPRC, PoPS) use breaks 0, 0.25, 0.50,
  0.75, 1.00. PR-curve axes always span the full 0–1 range (Fig. 2b); for
  other fractions the far limit follows the data (Fig. 4 stops at 0.7).
- Log axes label the transform in the title ("log~10~(…)") and break at
  integers; expansion is the default.
- A truncated baseline (Fig. 3g, 90–100 %) is allowed only for a
  filled/stacked bar whose remainder is the visible complement, with every
  tick label shown and the truncation stated in the caption.
- Do not use `coord_cartesian(expand = FALSE)` on point or curve panels; do
  not use `expand = c(0, 0)` on both ends of a bar chart (bars would touch the
  top of the panel).

In ggplot: `scale_y_anchored()` / `scale_x_anchored()` for bars,
`scale_y_fraction()` / `scale_x_fraction()` for 0–1 metrics on unanchored
marks, and plain `scale_*_continuous()` (defaults) for everything else.

### 1.6 Spacing and density

- Tick spacing: 3–5 major ticks per axis (0, 0.25, 0.50, 0.75, 1.00 for
  fractions). No minor ticks. Ticks point outward and are **2.4 pt** long at
  0.25 pt (measured on Fig. 2b/2c); grouped bar charts carry one 2.5 pt tick
  under the axis line at each group centre instead of per-bar ticks.
- Axis title sits 4–6 pt from tick labels; tick labels sit 2 pt from ticks.
- Bars in grouped charts touch (width = slot); groups separated by ~1.7 slots.
  Single-series bar charts use ~90 % of the slot.
- Legends are placed in whitespace inside the figure grid, never overlapping
  data; a legend shared across a figure appears once (Fig. 2, right of e/f).
- A row of related panels shares one x-axis title when the axes are identical
  (Fig. 4, ED Fig. 8).

---

## 2. Components

Each component lists anatomy, the spec measured from the paper, and the
constraints that keep new panels consistent.

### 2.1 Conceptual overview / schematic (Fig. 1a, 2a, 3b)

Anatomy: numbered stage headings "(1) …" in 6.5 pt across the top; boxed
stages on a light grey ground; black 0.9 pt arrows with 90° elbows; a single
magenta box for the model; simplified data tables in Courier.

- Grey containers `#DEDFE1`, no border, 6 pt inner padding, square corners —
  around data-like content only (fragment files, count matrices, signal
  tracks, feature boxes). Cartoon embeddings (the UMAP cell clouds of Fig. 1a
  and 3b) sit directly on white with a small notched-arrow axis pair.
- Signal tracks drawn as smooth peaks: light `#85C0E1` unselected, dark
  `#006CAD` selected; genes as green `#54AB58` arrows; contact as a black
  0.9 pt curve.
- Text inside containers 5.1–6 pt regular; math set in the same sans.
- Cartoon UMAP dots use pure RGB (`#1F84FF`, `#F40014`, `#00BF2B`) only when
  the colours carry no data meaning; grey `#B9C0CB` blobs indicate grouping.
- Maximum three visual weights: hairline (0.25), connector (0.9), emphasis
  (1.85 pt).
- **One arrow shape everywhere.** Gene arrows, TSS glyphs, the promoter-class
  arrow and the Fig. 2a table glyphs all use the same swept, notched head,
  proportional to the stem: head length 4 × stem, head width 3.25 × stem,
  back-edge notch 1.2 × stem. Stems: 1.85 pt (gene, green), 3.7 pt (emphasis),
  0.92 pt (TSS glyphs, grey), 2.0 pt (table glyphs, grey), 0.93 pt (black
  elbow connectors, whose 3.7 × 3.1 pt head follows the same ratios).
  Implemented as `svgkit.Canvas.arrow_poly(stem=…)` and, in R, `arrow_grob_sce2g()` /
  `annotate_arrow_sce2g()` (head width 0.82 × length, notch 0.24 × length).
- Cartoon UMAP schematic (Fig. 3b): 3.2 pt dots in the cartoon RGB colours,
  meta-cells as `#D3D3D3` blobs behind the dots, strategies separated by
  0.5 pt `#BEBEBE` 2/2-dashed verticals, two-line 6.5 pt headings centred over
  each block, and a ~20 pt UMAP axis pair with the notched arrowheads.
- Reference implementations: `figures/fig1a_schematic.py` rebuilds Fig. 1a as
  print-size SVG from coordinates measured in the PDF (1 unit = 1 mm,
  123.5 × 83.6 mm), with reusable helpers for peaks, gene arrows, elbow
  connectors, containers and mono tables; `figures/fig3b_strategies.py` draws
  Fig. 3b from digitised dot and blob geometry (`figures/data/fig3b_geometry.json`).

### 2.2 Locus plot (Fig. 1b–d, 3e, 5f; ED Fig. 3, 10)

Anatomy, top to bottom: scale bar with kb label → gene models (black exons,
grey introns, arrowheads for strand; gene names italic 6 pt) → optional
credible-set blocks and variant lollipops in teal → per-cell-type row pairs
(arcs above, ATAC track below) → optional cell-type × element heatmap →
score colour bar.

- Each cell-type block: arcs 22 pt tall, ATAC track 6 pt tall, 4 pt gap.
  Track name left-aligned in 6 pt with "ATAC-seq" in 5.1 pt under the first
  row only.
- Arcs: colour and width scale with score (light thin → dark thick). Arcs
  anchor at element midpoint and gene TSS; small tick at element foot. Every
  arc in a row peaks at the **same height** (~0.95 of the row pitch) whatever
  its span, so short links read as narrow full-height loops (Fig. 1b, 3e);
  `arc_paths()` / `geom_arc_sce2g()` draw them as quadratic curves.
- ATAC tracks are normalised per cell type, no y-axis, filled `#006CAD`.
- Variant of interest: dashed vertical line `#003A4E` spanning all rows,
  rsID label in teal at the top.
- Heatmap rows: one per cell type, labelled 4.2 pt (2.7 pt in Fig. 1b is below
  legibility — avoid). Rows of interest carry a `#DEDFE1` background band.
  Cells: red ramp, no border. Colour bar horizontal, 0.2–0.8 ticks, labelled
  "scE2G^Multiome score".
- Row order: fixed global cell-type order (cell lines → islets → BMMC/PBMC by
  lineage), identical across every locus plot.
- In the compact single-gene variant (Fig. 3e) the arcs run from a 0.8 pt TSS
  tick to a 0.8 pt element foot at 0.12–0.61 pt, colour and width ∝ score;
  the CRISPR-positive arc is `#9B241C` at 0.61 pt. Rows sit on
  a 15.6 pt pitch; the scale bar is 0.5 pt with round caps and 0.4 × 2.6 pt end
  ticks; tested elements are 5.9 pt tall `#DC6464` blocks under a 9 pt ATAC
  track.
- Reference implementations: `figures/fig1b_locus.R` rebuilds Fig. 1b (EBF1)
  in ggplot2 from the measured layout, with the 49 cell-type rows, superscript
  source markers, score-scaled arcs and the shared red ramp;
  `figures/fig3e_locus.R` draws Fig. 3e (TFRC) from the digitised geometry in
  `figures/data/fig3e_locus.json`, including the ATAC profile.

### 2.3 Precision–recall and enrichment–recall curves (Fig. 2b, 2e, 2g, 3c)

- Square plot area; x = Recall, y = Precision (or Enrichment). Axis range
  0–1 with ticks at 0.25 steps; recall axis title carries *n* ("Recall
  (*n* = 461)").
- One curve per predictor, 0.42 pt, palette colour; point at the model's
  suggested threshold as a ⌀ 3.7 pt filled circle; baseline-only predictors
  appear as points without curves.
- Reference lines: black dashed vertical at the chosen operating recall
  (0.7), grey dashed horizontal at random-classifier precision, both 0.5 pt
  with a 2/2 pt dash.
- Curves labelled at their end with the index number in matching colour;
  dashed leaders resolve crowding. No boxed legend inside the plot.
- For binary predictors show error bars (95 % CI) through the point: vertical
  for precision, horizontal for recall (Fig. 2h, i).
- Label placement is by hand in three groups: strong curves labelled to the
  right of their crossing with the recall-0.7 line, mid-strength curves in the
  open area below the crossing, the weakest stacked bottom-left; each leader
  is a dashed 0.42 pt segment in the curve's colour ending just short of the
  number. Reference implementation: `figures/fig2b_pr.R`.

### 2.4 Bar plots

**Grouped performance bars (Fig. 2c, 2d, 3c right; ED Fig. 4)**

- Bars ordered by the "All" group's ranking, kept fixed across groups; fill =
  palette colour, no outline; bars **touch** within a group (3.4 pt each in
  Fig. 2c) and groups are separated by ~1.7 bar widths.
- Mean ± 95 % CI error bar centred on each bar, 0.25 pt with 0.9 pt caps.
- Index numbers sit in a strip between the bar baseline and the x-axis line,
  which is drawn ~5.7 pt below zero; no x ticks.
- Group label and sample sizes in a three-row table below the axis, values
  centred under each group, right-aligned row headers ("Distance to TSS (kb):
  / CRISPR positives: / CRISPR negatives:") in 6.5 pt ending under the left
  edge of the first group.
- Y from 0; ticks 0, 0.25, 0.50, 0.75 labelled "0", "0.25", "0.50", "0.75".
- Reference implementation: `figures/fig2c_bars.R`, plotting values digitised
  from the paper (`figures/data/fig2c_values.csv`).

**Horizontal ranked bars (Fig. 3a, 5g)**

- Categories on y, highest at top; fill `#445369` (Fig. 3a) or `#6A768B`
  with `#B8C0CC` for de-emphasised rows (Fig. 5g); bars 0.9 of the row.
- The x axis keeps the default 5 % padding on both ends (bars start 0.05 units
  right of the y axis; Fig. 5g), unlike vertical bar charts.
- Significance stars right of the bar; exact *P* values as a right-hand
  text column headed *P*~bootstrap~.
- Features in the final model are underlined in the category label.
- Forward-selection variant (Fig. 3a): bar = ΔAUPRC with 95 % CI, black dot =
  cumulative AUPRC, grey dashed line at the final model's AUPRC, stars rotated
  90° beside the CI (red for negative deltas), *P*~bootstrap~ column at
  x ≈ 0.78; dot ⌀ 2.3 pt, stars 7.5 pt. Reference implementation:
  `figures/fig3a_features.R` with digitised values in
  `figures/data/fig3a_values.csv`; Fig. 5g in `figures/fig5g_pops.R`.

**Paired bars (Fig. 3d)** — two colours only, `#792374` and `#1E871B`,
touching, no gap inside a pair; legend as square swatches below the plot.

**Quantile / filled bars (Fig. 3g right)** — bars 0.66 of the slot with a
0.25 pt black outline; the "remainder" up to 100 % is a white outlined box;
y axis 90–100 % with **no** expansion (the 100 % line is the top of the
remainder); counts in 4.5 pt inside the remainder, or just above the 100 %
line when there is no room; 6.5 pt category labels; no x ticks. Single-cell
bar `#641E61`, meta-cell `#00783D`. Reference implementation:
`figures/fig3g_boxbar.R` (box plot + bars joined with patchwork).

### 2.5 Point-range (dot plot with CI) (Fig. 5d, e)

- Categories on y sorted by estimate, 6.5 pt row pitch; ⌀ 6.2 pt filled dot
  at the estimate, 0.42 pt CI line with 1.3 pt end caps, both in the category
  colour (`#260032` significant/dark, `#9A3983` others). Measured on Fig. 5d.
- Non-selected biosamples collapsed into one row as a grey beeswarm
  (`#A7A9AB`, ⌀ 2.9 pt) in a second, 12.4 mm panel below the 18.8 mm row
  panel; the y axis line runs through both, the x axis only under the swarm.
- Stars right of the CI for adjusted significance; *n* in the x-axis title.

### 2.6 Scatter plots (Fig. 4; ED Fig. 2c, 7, 8, 9b)

- Points ⌀ 4.3 pt, alpha 0.7, no outline, so overlap reads as density;
  single colour (`#6B2063`) when one model is shown.
- Horizontal dashed reference lines (0.5 pt, 2/2 dash) for competing models
  in *their* colours (ABC `#008F99`, Distance `#6A768B`); black dashed
  vertical for the recommended data threshold. The shared "Model" legend
  beneath the row shows each as a 12 pt dashed key or a point
  (`figures/fig4_legend.R`).
- Panels 39.5 mm square; y from 0 with the default 5 % padding both ends
  (0–0.7 for AUPRC/precision, 0–0.75 for recall); x breaks at integers.
- Reference implementation: `figures/fig4a_depth.R` (all 588 points
  digitised in `figures/data/fig4a_values.csv`).
- Small multiples share a y-axis title per row and an x-axis title per
  column; equal plot sizes; ~24 pt gutters.
- Log axes labelled explicitly as "log~10~(…)"; ticks at integers.

### 2.7 UMAP (Fig. 3f; ED Fig. 9)

- No axes: a short L-shaped arrow pair in the lower-left (15 pt arms, 0.5 pt,
  the paper's notched heads 3 × 2.4 pt, see §2.1), labelled "UMAP1" / "UMAP2"
  in 5.5 pt starting at the corner with ~3 pt clearance. In ggplot add
  `theme_embedding()` on top of `theme_sce2g()`: it blanks the axis lines,
  draws the arrow pair from the panel's bottom-left corner with
  `arrow_grob_sce2g()`, and forces a square aspect. Only the first panel of a
  multi-panel strip carries the axes. For a single arrow at a data position
  (Fig. 5h) use `annotate_arrow_sce2g()`.
- Points ⌀ ≈ 1 pt. Categorical panels draw them as filled circles with a
  hairline white outline (so neighbouring hues stay separable); continuous
  panels draw plain filled points with no outline. Three encodings side by
  side at equal size: categorical (many-hue; meta-cell centroids as 2.6 pt
  black discs with a 1 pt dot of the meta-cell colour at the centre),
  continuous expression (grey → purple), binary state (grey → red).
- Legends sit **inside** the panel at the lower right: a 4 × 18 pt vertical
  bar with white ticks, the title rotated 90° to the *left* of the bar
  (5.5 pt), tick labels to its right (4.5 pt). A binary state is still drawn
  as a gradient bar labelled with its two states.
- Reference implementation: `figures/fig3f_umap.R`.

### 2.8 Box plot (Fig. 3g left; ED Fig. 5g, h)

- Fill `#D3D3D3` (or category colour when comparing methods), 0.64 pt black
  outline, 1 pt median, whiskers 1.5 × IQR drawn **dashed** (1.9 pt on /
  3.2 pt off with round caps, so each dash reads ~2.5/2.5; ggplot
  `whisker.linetype = "34", lineend = "round"`) with staples half the box
  width; no outlier points. Set `colour = "black"`: ggplot's default box
  outline is grey20.
- Box width 80 % of the slot. Two label rows beneath the axis (category, then
  *n*, 5.5 pt) with 6.5 pt row headers ("ATAC state / Cell number")
  right-aligned at the y axis.
- *P* value (5.5 pt) right-aligned at the panel's top-right corner, baseline
  just above the tallest whisker. Reference: `figures/fig3g_boxbar.R`.

### 2.9 Violin with point-interval (Fig. 5b, c; ED Fig. 9a)

- Violin fill = category colour, no outline; inner black point-interval:
  ⌀ 3.05 pt median dot, thick 1.79 pt line for the central 66 %, thin 0.77 pt
  line for the central 95 % (ggdist "point interval"; measured on Fig. 5b).
- Comparison brackets: 0.25 pt black with 1.1 pt down-ticks at the ends,
  8.5 pt stars centred just below the bar, exact *P* value in 4.5 pt above it;
  stacked brackets step 10 pt apart. Panel 34.7 x 50.1 mm.
- Category labels two lines ("Same / lineage") with *n* on a third line.

### 2.10 Beeswarm (Fig. 5d, e "Other biosamples")

- Used only for the reference distribution; grey `#A7A9AB`, ⌀ 2.9 pt,
  no outline, symmetric swarm about the category line.
- The layout is the beeswarm **swarm** packing (every point touches a
  neighbour at one diameter; points sit on the centre line whenever free),
  not quasirandom. `ggbeeswarm::geom_beeswarm(method = "swarm")` reproduces
  it; `cex` must be tuned to the panel size (≈ 6 for a 12 mm panel) so that
  dots touch rather than overlap.

### 2.11 Heatmaps

**Symmetric correlation heatmap (Fig. 5a)** — 3.23 pt square cells, no visible
borders; columns run in the reverse order of the rows, so the two metrics split
along the anti-diagonal (red Jaccard upper left, blue Pearson lower right) and
the divider cells are navy `#001636`; dataset annotation as a 2.9 pt side bar;
two 7.8 x 38 pt colour bars (0-0.6) stacked to the right with 1.5 pt white edge
ticks; row/column labels 2.85 pt (superscripts 1.9 pt), columns rotated 60°;
legend keys 4.3 pt, legend titles 5.2 pt. Built in `figures/fig5a_heatmap.R`
from the digitised cells.

**Value grid (ED Fig. 7c, d, g, i)** — single-hue ramp, value printed inside
each cell in 4 pt (white on dark cells, black on light), axis titles
"Fragments per cell" / "# cells".

**Clustered enrichment heatmap (ED Fig. 9d)** — blue ramp capped at a stated
maximum, significant cells marked with a small dot, marginal bar plots
aligned to rows/columns, cell borders `#E4E4E4` 0.29 pt.

### 2.12 Motif logo (Fig. 5h)

- Bits on y (0–2, no expansion), positions on x, letters A `#00904A`, C `#1A5F9D`
  (measured; not ggseqlogo's default blue), G `#F1AD30`, T `#D92641`. Panel
  37.2 x 17.7 mm, 4.5 pt tick labels.
- Reference and alternate sequences in 8 pt Courier `#8B98AC` beneath at a
  fixed 6.1 pt advance (not aligned to the logo slots); the variant base black
  (ref) and `#9A1D20` (alt); 4 pt FIMO scores in parentheses to the right;
  rsID (6.2 pt) under a 3 pt closed upward arrow. Built in `figures/fig5h_logo.R`.

### 2.13 Tables inside figures (Fig. 2a)

- 0.5 pt `#202020` rules, `#F3F3F3` fill across the whole table, 5.5 pt text,
  6.5 pt title above; mini-glyphs in the first column are drawn in the mid
  grey (`#6E6E6E`) so they read as icons, not data.
- Text setting: line pitch 1.2 × size (6.6 pt at 5.5 pt); text cells
  left-aligned 5 pt in from their column rule; a multi-line block is centred
  vertically in its row; first-column labels are centred over their glyph;
  headers centred in their column.
- Reference implementation: `figures/fig2a_table.py` (uses the shared
  `figures/svgkit.py` emitter).

### 2.14 Legends

- Line series: 6 pt colour dash (Fig. 2) or a 12 pt dashed line for
  reference lines (Fig. 4); point series: a disc at the data size; bar
  series: 6 × 6 pt square.
- Entries prefixed with the index number "(1) scE2G^Multiome".
- Legend titles and group headings ("Predictors", "Developed in this study")
  are regular weight, 6.5 pt; the only bold text in a figure is the panel
  letter.
- Single shared legend per figure where possible.
- Reference lines and data points that belong to the same "Model" legend
  share one key row (Fig. 4): map both geoms to the same colour aesthetic,
  then set `guide_legend(override.aes = list(linetype = c("33", "33", "blank"),
  shape = c(NA, NA, 16)))` so each entry shows only its own glyph. See panel e
  in `examples/make_examples.R`.

### 2.15 Statistical annotation

- `*` *P* < 0.05, `***` *P* < 0.001, "n.s." otherwise. Stars are **bold,
  8.5 pt**, centred just below the bracket bar (baseline 2.5 pt under it);
  the exact *P* (4.5 pt, scientific notation "× 10^−17^" with a 0.7× superscript)
  sits above the bar so neighbouring labels share a baseline. Brackets are
  0.25 pt with 1.1 pt down-ticks; stacked brackets step 10 pt.
- In ggplot write the exponent as a plotmath superscript
  (`italic(P)*' = 2.45 × '*10^'–23'`) rather than a separately placed run:
  Helvetica sets the mantissa wider than Graphik and the runs collide.
- Test and correction named in the caption, never in the panel.
- Black stars = positive/expected direction; red `#E82829` stars = negative
  direction (Fig. 3a).

---

## 3. Accessibility

### 3.1 Contrast

Measured WCAG contrast ratios against white (see §1.2.1). Only 1, 2, 7, 8,
15, 17, 18 and the dark category colours reach 4.5:1; use those for
**coloured text**. Colours below 3:1 (4, 6, 10–14, 19) must never carry text
or thin lines alone — use them for filled areas ≥ 3 pt wide, and print the
index number in black beneath.

### 3.2 Colour-vision deficiency

Simulated deuteranopia shows the two headline colours, scE2G^Multiome purple
(`#792374`) and scE2G^ATAC teal (`#006479`), collapse to nearly the same
dark blue-grey. ABC `#0096A0` also merges with ARC-E2G `#B778B3` and Closest
TSS `#6E788D`. The paper mitigates this with **redundant encoding** — every
curve, bar and point carries its index number — and this is mandatory in new
panels. Additional rules:

- Never distinguish two series by hue alone; add number, label, line style
  or position.
- Prefer lightness contrast between series that must be compared (dark
  purple vs light teal works; dark purple vs dark teal does not).
- Sequential scales are single-hue lightness ramps (safe); the red/blue split
  heatmap is safe because the two halves are spatially separated.

### 3.3 Legibility

- Nature minimum is 5 pt at final size. The paper's 2.5 pt bar indices and
  2.7 pt heatmap row labels fall below this; do not repeat that. Use 4.5 pt
  as the absolute floor and 5.5 pt as the default.
- Minimum stroke 0.25 pt; hairlines under 0.2 pt disappear in print.
- Dashed lines: reference lines 2 pt on / 2 pt off at 0.5 pt; leaders 2.5 pt
  on / 0.85 pt off at 0.42 pt. Never use ggplot's default "dashed"/"22" at
  these weights, which prints as dots.

---

## 4. Information density

The figures are dense but disciplined. Targets:

| Panel type | Series / items before splitting |
|---|---|
| PR curves | ≤ 19 curves only with end-labels; ≤ 6 for readability |
| Grouped bars | ≤ 4 groups × 14 bars per row |
| Locus plot | ≤ 6 cell-type rows of arcs + ≤ 45 heatmap rows |
| Small multiples | 3 × 2 grid at full width, each ≥ 120 pt wide |
| Violin/box | ≤ 4 categories per panel |
| Heatmap | ≤ 45 × 45 labelled cells |

Strategies used in the paper to hold density without clutter: shared
legends; numbers instead of names on marks; sample sizes moved into axis
titles; explanatory sub-axis tables under bar groups; suppression of top/right
spines and grid; light-grey de-emphasis for reference distributions.

---

## 5. Panel checklist

1. Fixed colour and index for every predictor; index printed on the mark.
2. Open axes, 0.25 pt, outward ticks, no grid, white background.
3. Text ≥ 4.5 pt (5.5 default), Helvetica/Arial, sentence case, italic
   *n*/*P*/gene names.
4. Error bars 95 % CI, described in caption.
5. Reference lines dashed: black = chosen threshold, grey = chance,
   coloured = competing model.
6. Legend once per figure, in whitespace, numbered entries.
7. Panel letter 9 pt bold, outside plot, top-left.
8. Passes the CVD check: no two compared series differ only by hue.

---

## 6. Implementation

### ggplot2

`R/theme_sce2g.R` implements the system (reads `tokens.json`):
`theme_sce2g()`, `theme_embedding()` (UMAP/PCA add-on: notched corner arrows,
square aspect), `scale_colour_predictor()` / `scale_fill_predictor()` (keyed by
index or method name), `scale_fill_sce2g_seq("red"|"blue"|"purple"|"teal")`,
`scale_*_anchored()` / `scale_*_fraction()`, `stat_pointinterval_sce2g()`,
`sig_bracket()`, `arc_paths()` / `geom_arc_sce2g()` for locus arcs,
`arrow_grob_sce2g()` / `annotate_arrow_sce2g()` for the paper's arrow, and
`save_panel()` (mm sizing; PDF via cairo_pdf, PNG via cairo — ragg silently
drops glyphs under ~8 px, i.e. the 1.9 pt superscripts at 300 dpi).
Deliverables are SVG (editable), PDF and a 900 dpi PNG per panel — `save_panel(p, "figures/name", w, h)` writes all three.
`examples/make_examples.R` builds one panel of every type into `examples/out/`;
`figures/*.R` are the verified reference implementations.

```r
source("R/theme_sce2g.R"); theme_set(theme_sce2g())
ggplot(df, aes(recall, precision, colour = predictor_index)) +
  geom_path(linewidth = pt_mm(0.42)) + scale_colour_predictor(guide = "none")
```

**Unit conversions that bite.** ggplot's `linewidth` is nominally mm but draws
at 0.75 mm per unit (R's lwd is 1/96 in): convert points with
`pt_mm(pt) = pt / .pt / 0.75`. Point diameters follow the same factor: a
`d` pt disc is `geom_point(size = d / .pt / 0.75, shape = 16, stroke = 0)`
(shape 19, the default, adds a 0.42 pt outline). Text is the exception:
`size = pt / .pt` is exact. Dash patterns (`lty` digits) are in units of
max(lwd, 1) = 0.75 pt for the hairlines used here, so a 2/2 pt dash is `"33"`
(2.25 pt) and 2.5/0.85 is `"31"`. Superscripts via plotmath are 0.7×, which
matches the paper.

**Composition.** Legends, keys and dense labelled panels (Fig. 4 legend, Fig.
5a) are easiest to reproduce as a `theme_void()` ggplot in absolute point
coordinates (x right, y *negated* — `coord_cartesian` silently sorts reversed
limits), one data unit per point, with the canvas sized to match. Text placed
with `vjust = 0` sits on its baseline for strings without descenders; anchor
multi-line blocks by their tops (`vjust = 1`) so descenders do not shift lines.

**Verifying.** Convert both PDFs to SVG (`mutool draw -o out.svg in.pdf`) and
compare stroke widths, dash arrays, glyph scales and circle diameters directly;
`figures/data/*.csv` hold the digitised values the reference panels draw.

### matplotlib

```python
import json, matplotlib as mpl
tok = json.load(open("tokens.json"))
mpl.rcParams.update({
  "font.family": "Helvetica", "font.size": 6.5,
  "axes.titlesize": 6.5, "axes.labelsize": 6.5,
  "xtick.labelsize": 5.5, "ytick.labelsize": 5.5, "legend.fontsize": 5.5,
  "axes.linewidth": 0.25, "xtick.major.width": 0.25, "ytick.major.width": 0.25,
  "xtick.direction": "out", "ytick.direction": "out",
  "axes.spines.top": False, "axes.spines.right": False,
  "lines.linewidth": 0.42, "axes.grid": False, "figure.dpi": 300,
  "pdf.fonttype": 42, "svg.fonttype": "none",
})
PRED = tok["color"]["predictor"]   # {"1": {"name": ..., "hex": ...}, ...}
```
