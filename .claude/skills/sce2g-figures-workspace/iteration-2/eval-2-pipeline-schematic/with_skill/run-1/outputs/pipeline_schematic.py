"""Conceptual pipeline schematic in the scE2G Fig. 1a style (no data), 120 x 45 mm.

Three numbered stages left to right: (1) single-cell multiome data (cartoon UMAP on white with a
notched-arrow axis pair), (2) candidate enhancer-gene pairs (peak track + green gene arrow in a
grey container), (3) model prediction (the magenta model box). Coordinates in PDF points; svgkit
converts to mm. Colours and sizes from tokens.json.

    python pipeline_schematic.py      # -> pipeline_schematic.{svg,pdf,png}
"""
import json
import random
from pathlib import Path

from svgkit import Canvas

HERE = Path(__file__).resolve().parent
TOK = json.loads((HERE / "tokens.json").read_text())
SCH, SZ, SW = TOK["color"]["schematic"], TOK["type"]["size_pt"], TOK["stroke_pt"]
BOX, PEAK_L, PEAK_D, GENE, MODEL = SCH["container"], SCH["peak_light"], SCH["peak_dark"], SCH["gene_arrow"], SCH["model_box"]
DOTS = SCH["cartoon_dots"]
FS_HEAD, FS_LABEL, FS_NOTE = SZ["panel_title"], SZ["tick_label"], SZ["annotation"]   # 6.5 / 5.5 / 4.5
HAIR, CONN, EMPH = SW["axis"], 0.93, SW["schematic_emphasis"]                       # the three weights

W, H = 120 / Canvas.S(1), 45 / Canvas.S(1)          # 340.2 x 127.6 pt
c = Canvas(0, 0, W, H)


def arrow(x0, y0, x1, y1, stem, head_len=None, head_w=None, notch=None):
    """The paper's notched arrow in any direction (head 4 x stem long, 3.25 x wide, notch 1.2 x)."""
    head_len, head_w, notch = head_len or 4 * stem, head_w or 3.25 * stem, notch or 1.2 * stem
    dx, dy = x1 - x0, y1 - y0
    L = (dx * dx + dy * dy) ** 0.5
    ux, uy, px, py = dx / L, dy / L, -dy / L, dx / L
    bx, by = x1 - head_len * ux, y1 - head_len * uy
    nx, ny = x1 - (head_len - notch) * ux, y1 - (head_len - notch) * uy
    head = [(x1, y1), (bx + head_w / 2 * px, by + head_w / 2 * py), (nx, ny), (bx - head_w / 2 * px, by - head_w / 2 * py)]
    return [c.line(x0, y0, nx, ny, c.ink, stem), c.poly(head, c.ink)]


def peak(x0, x1, y, h, fill):
    xm, w = (x0 + x1) / 2, x1 - x0
    return c.path([("M", x0, y), ("C", x0 + .4 * w, y, xm - .17 * w, y - h, xm, y - h),
                   ("C", xm + .17 * w, y - h, x1 - .4 * w, y, x1, y), ("Z",)], stroke="none", fill=fill)


# ---- headings ----------------------------------------------------------------------------------
for cx, lines in ((52, ("(1) Single-cell", "multiome data")), (180, ("(2) Candidate", "enhancer–gene pairs")),
                  (300, ("(3) Model", "prediction"))):
    c.add(c.text(cx, 12, lines, FS_HEAD, "middle"))

# ---- stage 1: cartoon UMAP on white, small notched axis pair -----------------------------------
rng = random.Random(3)
for (cx, cy, rx, ry), col in (((44, 52, 13, 8), DOTS[0]), ((76, 66, 11, 10), DOTS[1]), ((54, 86, 15, 7), DOTS[2])):
    for _ in range(14):
        t, r = rng.uniform(0, 6.2832), rng.uniform(0, 1) ** .5
        c.add(c.circle(cx + rx * r * __import__("math").cos(t), cy + ry * r * __import__("math").sin(t), 3.19, col))
c.add(*arrow(20, 102, 20, 82, 0.5, 2.985, 2.44, 0.71), *arrow(20, 102, 40, 102, 0.5, 2.985, 2.44, 0.71))
c.add(c.text(20, 109, "UMAP1", FS_LABEL))
c.add(c.text(0, 0, "UMAP2", FS_LABEL).replace("<text ", f'<text transform="rotate(-90 {c.X(16.5):.3f} {c.Y(101):.3f})" ')
      .replace('x="0.000" y="0.000"', f'x="{c.X(16.5):.3f}" y="{c.Y(101):.3f}"'))

# ---- connectors --------------------------------------------------------------------------------
c.add(c.arrow_poly(104, 128, 68, c.ink, stem=CONN), c.arrow_poly(236, 260, 68, c.ink, stem=CONN))

# ---- stage 2: peak track + gene arrow in a grey container -------------------------------------
c.add(c.rect(132, 30, 232, 106, BOX))
yb = 78
c.add(peak(140, 158, yb, 7, PEAK_L), peak(162, 188, yb, 17, PEAK_D), peak(191, 204, yb, 10, PEAK_L))
c.add(c.line(140, yb, 222, yb, sw=HAIR))
for x0, x1, col in ((140, 158, PEAK_L), (162, 188, PEAK_D), (191, 204, PEAK_L)):
    c.add(c.rect(x0, yb + 1.9, x1, yb + 4.5, col))
c.add(c.arrow_poly(208, 226, yb + .3, GENE, stem=EMPH))
c.add(c.path([("M", 175, yb - 18.5), ("C", 183, yb - 30, 204, yb - 28, 208, yb - 2.5)], sw=0.5))   # enhancer apex -> TSS
c.add(c.text(175, 42, "ATAC signal", FS_LABEL, "middle"))
c.add(c.text(175, 91, "Enhancer", FS_NOTE, "middle"), c.text(217, 91, "Gene", FS_NOTE, "middle"))
c.add(c.text(182, 101, "Peak–gene pairs within 5 Mb", FS_NOTE, "middle"))

# ---- stage 3: the model box and its output -----------------------------------------------------
c.add(c.rect(274, 40, 326, 64, MODEL))
c.add(c.text(300, 54.3, "Model", FS_HEAD, "middle", fill="#FFFFFF"))
c.add(*arrow(300, 64, 300, 80, CONN))
c.add(c.text(300, 90, ["Enhancer–gene", "prediction score"], FS_LABEL, "middle"))

c.save(HERE / "pipeline_schematic.svg")
