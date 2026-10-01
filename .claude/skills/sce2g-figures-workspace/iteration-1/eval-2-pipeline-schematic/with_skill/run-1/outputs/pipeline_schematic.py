"""Pipeline schematic in the scE2G Fig. 1a style (no data): three numbered stages, left to right.

    python pipeline_schematic.py   # -> pipeline_schematic.{svg,pdf,png}, 120 x 45 mm

Coordinates are page points on a 340 x 128 pt frame; svgkit converts to mm. Colours from tokens.json.
"""
import json
import math
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from svgkit import Canvas  # noqa: E402

TOK = json.loads((HERE / "tokens.json").read_text())
SCH, NEU = TOK["color"]["schematic"], TOK["color"]["neutral"]
BOX, PEAK_L, PEAK_D, GENE, MODEL = SCH["container"], SCH["peak_light"], SCH["peak_dark"], SCH["gene_arrow"], SCH["model_box"]
DOTS = SCH["cartoon_dots"]
FS_HEAD, FS_LABEL, FS_TICK = 6.5, 5.1, 5.5

W, H = 340.0, 128.0                     # 120 x 45 mm
BOXES = [(4, 96), (124, 216), (244, 336)]   # x0, x1 of the three containers
BY0, BY1 = 26.0, 122.0
c = Canvas(0, 0, W, H)


def rot(item, deg, x, y):
    return f'<g transform="rotate({deg} {c.X(x):.3f} {c.Y(y):.3f})">{item}</g>'


def peak(x0, x1, y_base, h, fill):
    xm, w = (x0 + x1) / 2, x1 - x0
    return c.path([("M", x0, y_base), ("C", x0 + 0.40 * w, y_base, xm - 0.17 * w, y_base - h, xm, y_base - h),
                   ("C", xm + 0.17 * w, y_base - h, x1 - 0.40 * w, y_base, x1, y_base), ("Z",)], stroke="none", sw=0, fill=fill)


def axis_arrow(x0, tip, y):     # UMAP axis: 0.5 pt stem, 3.0 x 2.4 pt head (Fig. 3b)
    return c.line(x0, y, tip - 2.3, y, sw=0.5) + c.arrow_poly(tip - 3.0, tip, y, c.ink, stem=0.5, head_len=3.0, head_w=2.44, notch=0.71)


# ---- headings and containers ------------------------------------------------------------------
for (x0, x1), lines in zip(BOXES, (("(1) Single-cell", "multiome data"), ("(2) Candidate", "enhancer–gene pairs"),
                                   ("(3) Model", "prediction"))):
    c.add(c.text((x0 + x1) / 2, 10.0, lines, FS_HEAD, "middle"))
    c.add(c.rect(x0, BY0, x1, BY1, BOX))

# ---- stage 1: cartoon UMAP ---------------------------------------------------------------------
rng = random.Random(7)
for colour, (cx, cy, rx, ry) in zip(DOTS, ((36, 56, 12, 9), (72, 50, 9, 13), (62, 90, 13, 8))):
    for _ in range(13):
        a, r = rng.uniform(0, 6.2832), rng.uniform(0, 1) ** 0.5
        c.add(c.circle(cx + rx * r * math.cos(a), cy + ry * r * math.sin(a), 3.19, colour))
ax, ay = 14.0, 114.0
c.add(axis_arrow(ax, ax + 20.5, ay))
c.add(rot(axis_arrow(ax, ax + 20.5, ay), -90, ax, ay))
c.add(c.text(ax, ay + 6.5, "UMAP1", FS_TICK))
c.add(rot(c.text(ax - 2.0, ay, "UMAP2", FS_TICK), -90, ax - 2.0, ay))

# ---- stage 2: peak track with gene arrow -------------------------------------------------------
yb = 80.0
c.add(c.text(130, 40.0, "Pseudobulk ATAC", FS_LABEL))
c.add(peak(132, 148, yb, 7, PEAK_L), peak(153, 177, yb, 16, PEAK_D), peak(180, 194, yb, 10, PEAK_L))
c.add(c.line(132, yb, 200, yb, sw=0.25))
c.add(c.rect(132, yb + 1.9, 148, yb + 4.5, PEAK_L), c.rect(153, yb + 1.9, 177, yb + 4.5, PEAK_D), c.rect(180, yb + 1.9, 194, yb + 4.5, PEAK_L))
c.add(c.arrow_poly(193, 211, yb + 0.3, GENE, stem=1.85))
c.add(c.path([("M", 165, yb + 6), ("C", 165, yb + 20, 197, yb + 20, 197, yb + 6)], sw=0.5, dash=f"{c.S(2):.3f} {c.S(2):.3f}"))
c.add(c.text(165, yb + 27, "Enhancer", FS_LABEL, "middle"), c.text(200, yb + 27, "Gene", FS_LABEL, "middle"))

# ---- stage 3: model box and predictions --------------------------------------------------------
mx, my = 290.0, 54.0
c.add(c.rect(mx - 25, my - 12, mx + 25, my + 12, MODEL))
c.add(c.text(mx, my + 2.3, "scE2G", FS_HEAD, "middle", fill="#FFFFFF"))
c.add(c.line(mx, my + 12, mx, my + 24, sw=0.93))
c.add(rot(c.arrow_poly(mx - 2.0, mx + 2.7, my + 24, c.ink, stem=0.93), 90, mx, my + 24))
c.add(c.text(mx, my + 39, ["Genome-wide", "enhancer–gene", "predictions"], FS_LABEL, "middle"))

# ---- connectors between stages -----------------------------------------------------------------
ym = (BY0 + BY1) / 2
for (_, xa), (xb, _) in zip(BOXES[:-1], BOXES[1:]):
    c.add(c.line(xa + 2.5, ym, xb - 2.0, ym, sw=0.93))
    c.add(c.arrow_poly(xb - 3.0, xb + 0.8, ym, c.ink, stem=0.93))

c.save(HERE / "pipeline_schematic.svg")
