# Verifying a panel against the paper

## Why measure

On screen, a 0.19 pt axis and a 0.25 pt axis look identical; in print the
panel reads lighter than its neighbours. The three systematic errors that
produced every early discrepancy (line widths 25 % thin, discs 25 % small,
dashes 50 % long) were invisible until the rendered PDF was measured. So:
render, convert, measure, then look.

## The tool

```bash
python scripts/measure_svg.py panel.pdf
python scripts/measure_svg.py paper.pdf --page 7 --bbox 363 615 567 740
```

Requires `mutool` (MuPDF). It converts the page to SVG and reports, for the
whole page or a bounding box in PDF points (origin top-left, y down):

- `strokes` — (width pt, colour) with counts. Axes should be 0.25, curves 0.42,
  reference lines 0.5, box outlines 0.64.
- `dashes` — (dasharray, width). 2.25,2.25 at 0.5 is the reference dash.
- `text_sizes` — glyph scale in pt with counts. 5.5 / 6.5 / 4.5 for a data panel.
- `discs` — (diameter, fill) for circular paths. Compare with the intended
  diameter; cairo circles are two Béziers, the paper's are four — both handled.
- `fills` — distinct fill colours (a heatmap will have hundreds).

## Print scale

The journal printed Fig. 1 at 93 % and Fig. 5 at 95 % of design size, Figs.
2–4 at 100 %. A 5.22 pt label measured on a Fig. 5 panel is therefore the
5.5 pt token, a 0.24 pt axis is 0.25, a 0.73 pt interval line is 0.77. Always
divide a Fig. 1 or Fig. 5 measurement by its print scale before comparing with
tokens or building at design size.

## Digitising a reference panel

When a panel should carry the paper's data (bars, points, violins, heatmap
cells), digitise from the vector SVG rather than eyeballing:

1. `mutool draw -o page.svg paper.pdf <page>`; every mark is a `<path>` with a
   `transform="matrix(1,0,0,-1,x,y)"` (origin) or absolute coordinates under
   `matrix(1,0,0,-1,0,H)` with H = page height (790.866 pt for this paper),
   and text is `<use data-text="…" transform="matrix(size,0,0,-size,x,y)">`.
2. Find the axes: tick paths give the pt ↔ data mapping (two ticks suffice).
3. Filter marks by fill/stroke and region, convert to data units, write a CSV
   under `figures/data/` and read it from the script. Keep the raw page
   coordinates in the CSV when a panel is drawn in point space.
4. Colour-encoded values (heatmaps): decode the base64 colour-bar raster,
   sample its column, and map each cell colour to the nearest ramp position.

Gotchas met in practice: a circle path's origin is its left edge, not its
centre; rectangle origins may be the bottom edge; glyph runs are split at
every superscript; a single `C` command can carry several Bézier segments.

## SVG versus PDF

svglite (the `.svg` output) keeps text as text and strokes in pt, but it sizes
dash patterns at 1 pt per `lty` unit, so a `"33"` reference line reads
3/3 pt in the SVG and 2.25/2.25 pt in the PDF. The PDF is the print
reference; measure that one.

## Visual checks that measurement misses

- Nothing clipped at the canvas edge (fixed-size layouts overflow silently).
- Superscripts and exponents do not collide with the mantissa (use plotmath).
- Arc ends land on their element feet; bars start on the axis line.
- Category labels do not touch; if the paper's size is too tight for Helvetica,
  widen the panel to the paper's measured width before shrinking text.
