"""Pipeline schematic in the scE2G Fig. 1a style: (1) single-cell multiome data (cartoon UMAP on white),
(2) candidate enhancer-gene pairs (peak track + gene arrow in a grey container), (3) model prediction
(magenta model box). Coordinates in PDF points; svgkit emits a 120 x 45 mm SVG (+ PDF/PNG via Inkscape).

    python pipeline_schematic.py
"""
import json
import random
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from svgkit import Canvas  # noqa: E402

TOK = json.loads((HERE / "tokens.json").read_text())
SCH, FS = TOK["color"]["schematic"], TOK["type"]["size_pt"]
BOX, PEAK_L, PEAK_D, GENE, MODEL = SCH["container"], SCH["peak_light"], SCH["peak_dark"], SCH["gene_arrow"], SCH["model_box"]
DOTS = SCH["cartoon_dots"]
FS_HEAD, FS_LABEL, FS_NOTE = FS["panel_title"], FS["tick_label"], FS["annotation"]
MM = 25.4 / 72

c = Canvas(0, 0, 120 / MM, 45 / MM)          # 340.2 x 127.6 pt
CX = (58, 170, 285)                          # stage centres


def heading(x, lines):
    return c.text(x, 12, lines, FS_HEAD, "middle")


def peak(x0, x1, yb, h, fill):
    xm, w = (x0 + x1) / 2, x1 - x0
    return c.path([("M", x0, yb), ("C", x0 + .4 * w, yb, xm - .17 * w, yb - h, xm, yb - h),
                   ("C", xm + .17 * w, yb - h, x1 - .4 * w, yb, x1, yb), ("Z",)], stroke="none", fill=fill)


def axis_arrow(x0, y0, x1, y1):              # Fig. 3b UMAP axis: 0.5 pt stem, 3.0 x 2.4 pt notched head
    return c.connector([(x0, y0), (x1, y1)], sw=0.5, head_len=2.985, head_w=2.44, notch=0.71)


# ---- (1) single-cell multiome data: cartoon UMAP on white -------------------------------------
c.add(heading(CX[0], ["(1) Single-cell", "multiome data"]))
rng = random.Random(7)
for (cx, cy, rx, ry), col in zip(((48, 52, 11, 8), (80, 62, 9, 10), (62, 88, 12, 7)), DOTS):
    for _ in range(16):
        a, r = rng.uniform(0, 6.2832), rng.uniform(0, 1) ** .5
        import math
        c.add(c.circle(cx + rx * r * math.cos(a), cy + ry * r * math.sin(a), 3.19, col))
c.add(*axis_arrow(26, 108, 26, 87.5), *axis_arrow(26, 107.5, 46.5, 107.5))
c.add(c.text(26, 115, "UMAP1", FS_LABEL))
c.add(f'<text font-family="Helvetica, Arial, sans-serif" font-size="{c.S(FS_LABEL):.3f}" '
      f'transform="rotate(-90 {c.X(22.5):.3f} {c.Y(106.5):.3f})" x="{c.X(22.5):.3f}" y="{c.Y(106.5):.3f}">UMAP2</text>')

# ---- (2) candidate enhancer-gene pairs: peak track in a grey container ------------------------
c.add(heading(CX[1], ["(2) Candidate", "enhancer–gene pairs"]))
X0, X1, YB = 118, 222, 84
c.add(c.rect(X0, 36, X1, 110, BOX))
c.add(peak(126, 144, YB, 7, PEAK_L), peak(148, 176, YB, 16, PEAK_D), peak(180, 194, YB, 10, PEAK_L))
c.add(c.line(126, YB, 197, YB, sw=0.25))
c.add(c.rect(126, YB + 1.9, 144, YB + 4.5, PEAK_L), c.rect(148, YB + 1.9, 176, YB + 4.5, PEAK_D),
      c.rect(180, YB + 1.9, 194, YB + 4.5, PEAK_L))
c.add(c.arrow_poly(199, 216, YB + .3, GENE, stem=1.85))
c.add(c.text(162, YB + 10.5, "Enhancer", FS_NOTE, "middle"), c.text(207, YB + 10.5, "Gene", FS_NOTE, "middle"))
c.add(c.path([("M", 162, 66), ("C", 166, 54, 200, 54, 206, 81)], sw=0.25, dash="0.66 0.66"))  # 1.87/1.87 pt
c.add(c.text(170, 104.5, "ATAC signal, pseudobulk", FS_NOTE, "middle"))

# ---- (3) model prediction: the magenta model box ----------------------------------------------
c.add(heading(CX[2], ["(3) Model", "prediction"]))
c.add(c.rect(252, 58, 318, 88, MODEL))
c.add(c.text(CX[2], 75.5, "scE2G model", FS_HEAD, "middle", fill="#FFFFFF"))
c.add(c.text(CX[2], 100, ["Genome-wide", "enhancer–gene predictions"], FS_LABEL, "middle"))

# ---- connectors last, tips on the target edges -------------------------------------------------
c.add(*c.connector([(103, 73), (X0, 73)]), *c.connector([(X1, 73), (252, 73)]))

c.save(HERE / "pipeline_schematic.svg")
