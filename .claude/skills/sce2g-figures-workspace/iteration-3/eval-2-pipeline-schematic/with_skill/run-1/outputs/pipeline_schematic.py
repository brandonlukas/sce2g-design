"""Pipeline schematic in the scE2G Fig. 1a style: (1) single-cell multiome data (cartoon UMAP on
white) -> (2) candidate enhancer-gene pairs (peak track + gene arrow in a grey container) ->
(3) model prediction (magenta model box, predictions container). 120 x 45 mm, no data.

    python pipeline_schematic.py      # -> pipeline_schematic.{svg,pdf,png}
"""
import json, math, random, sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from svgkit import Canvas  # noqa: E402

TOK = json.loads((HERE / "tokens.json").read_text())
SCH, FS = TOK["color"]["schematic"], TOK["type"]["size_pt"]
BOX, PEAK_L, PEAK_D, GENE, MODEL = SCH["container"], SCH["peak_light"], SCH["peak_dark"], SCH["gene_arrow"], SCH["model_box"]
DOTS, PRED = SCH["cartoon_dots"], TOK["color"]["predictor"]["1"]["hex"]
FS_HEAD, FS_LABEL, FS_NOTE, FS_AXIS = 6.5, 5.1, 4.5, 5.5

c = Canvas(ox=0, oy=0, w=340, h=128)          # 120 x 45 mm in page points


def arrow(x0, y0, x1, y1, stem=0.5, head_len=2.985, head_w=2.44, notch=0.71):
    """Notched arrow in any direction (fig3b_strategies.py); defaults are the UMAP-axis size."""
    dx, dy = x1 - x0, y1 - y0
    L = (dx * dx + dy * dy) ** 0.5
    ux, uy = dx / L, dy / L
    px, py = -uy, ux
    bx, by = x1 - head_len * ux, y1 - head_len * uy
    nx, ny = x1 - (head_len - notch) * ux, y1 - (head_len - notch) * uy
    head = [(x1, y1), (bx + head_w / 2 * px, by + head_w / 2 * py), (nx, ny), (bx - head_w / 2 * px, by - head_w / 2 * py)]
    return [c.line(x0, y0, nx, ny, c.ink, stem), c.poly(head, c.ink)]


def connector(x0, x1, y, sw=0.93):
    """Horizontal 0.93 pt connector with the paper's swept head (stem-scaled, tip overlaps the box ~1 pt)."""
    return [c.line(x0, y, x1 - 2.0, y, c.ink, sw), c.arrow_poly(x1 - 2.0, x1 + 2.9 * sw, y, c.ink, stem=sw)]


def peak(x0, x1, y, h, fill):
    xm, w = (x0 + x1) / 2, x1 - x0
    return c.path([("M", x0, y), ("C", x0 + 0.40 * w, y, xm - 0.17 * w, y - h, xm, y - h),
                   ("C", xm + 0.17 * w, y - h, x1 - 0.40 * w, y, x1, y), ("Z",)], stroke="none", fill=fill)


def track(y, peaks, gene):
    """Peaks (x0, x1, h, fill) on a hairline baseline with element bars beneath, then a gene arrow."""
    out = [peak(x0, x1, y, h, f) for x0, x1, h, f in peaks]
    out.append(c.line(peaks[0][0], y, gene[0] + 2, y, c.ink, 0.25))
    out += [c.rect(x0, y + 1.9, x1, y + 4.5, f) for x0, x1, _, f in peaks]
    out.append(c.arrow_poly(gene[0], gene[1], y + 0.3, GENE, stem=1.85))
    return out


# ---- (1) single-cell multiome data: cartoon UMAP on white ---------------------------------------
c.add(c.text(50, 13, ["(1) Single-cell", "multiome data"], FS_HEAD, "middle"))
rng = random.Random(3)
for (cx, cy, rx, ry), col in zip(((42, 52, 13, 9), (78, 66, 10, 12), (56, 88, 14, 8)), DOTS):
    for _ in range(14):
        a, r = rng.uniform(0, 6.2832), rng.uniform(0, 1) ** 0.5
        c.add(c.circle(cx + rx * r * math.cos(a), cy + ry * r * math.sin(a), 3.19, col))
c.add(*arrow(16, 108, 16, 87.5), *arrow(16, 107.55, 36.4, 107.55))          # UMAP axis pair, 20.5 pt arms
c.add(c.text(15.85, 114.9, "UMAP1", FS_AXIS))
c.add(f'<text font-family="Helvetica, Arial, sans-serif" font-size="{c.S(FS_AXIS):.3f}" '
      f'transform="rotate(-90 {c.X(12.4):.3f} {c.Y(106.8):.3f})" x="{c.X(12.4):.3f}" y="{c.Y(106.8):.3f}">UMAP2</text>')

# ---- (2) candidate enhancer-gene pairs: peak track in a grey container --------------------------
c.add(c.text(170, 13, ["(2) Candidate", "enhancer–gene pairs"], FS_HEAD, "middle"))
c.add(c.rect(118, 30, 222, 118, BOX))
c.add(c.text(124, 41, "Pseudobulk ATAC peaks × genes", FS_LABEL))
c.add(c.path([("M", 156, 92), ("C", 156, 92 - 1.35 * 24, 200, 92 - 1.35 * 24, 200, 92)], c.ink, 0.5,
             dash=f"{c.S(2):.3f} {c.S(2):.3f}"))                                  # the candidate pair
c.add(*track(92, [(126, 142, 6, PEAK_L), (146, 166, 14, PEAK_D), (170, 184, 9, PEAK_L)], (198, 216)))
c.add(c.text(156, 105, "Enhancer", FS_NOTE, "middle"), c.text(206, 105, "Gene", FS_NOTE, "middle"))
c.add(c.text(178, 66, "Pair?", FS_NOTE, "middle"))

# ---- (3) model prediction: magenta model box, predictions container -----------------------------
c.add(c.text(290, 13, "(3) Model prediction", FS_HEAD, "middle"))
c.add(c.rect(262, 36, 318, 58, MODEL))
c.add(c.text(290, 49.3, "scE2G", FS_HEAD, "middle", fill="#FFFFFF"))
c.add(*arrow(290, 58, 290, 71, stem=0.93, head_len=3.72, head_w=3.02, notch=1.12))   # 0.93 pt connector, head 4 x stem
c.add(c.rect(240, 70, 340, 118, BOX))
c.add(c.text(290, 80, "Enhancer–gene predictions", FS_LABEL, "middle"))
c.add(*track(104, [(248, 262, 5, PEAK_L), (266, 282, 9, PEAK_D), (288, 300, 6, PEAK_L)], (312, 330)))
c.add(c.arc(274, 314, 104, 15, PRED, 0.9), c.arc(255, 314, 104, 9, PRED, 0.9))

# ---- stage connectors --------------------------------------------------------------------------
c.add(*connector(101, 119, 47), *connector(223, 263, 47))

c.save(HERE / "pipeline_schematic.svg")
