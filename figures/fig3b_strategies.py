"""Fig 3b — three strategies for computing the Kendall correlation: (1) single cells of one cell type,
(2) meta-cells of one cell type, (3) single cells of several cell types. A cartoon UMAP schematic drawn
as print-size SVG.

Geometry measured from page 7 of refs/sce2g-2026.pdf (points, page coordinates): every cell dot
(3.19 pt, cartoon RGB from tokens.json schematic.cartoon_dots) and the four grey meta-cell blobs are
digitised into figures/data/fig3b_geometry.json; 0.5 pt #BEBEBE 2/2-dashed separators at x = 336.0 and
411.7; two-line 6.5 pt headings centred over each strategy; a ~20 pt UMAP axis pair with the paper's
notched arrowheads (0.5 pt stems, 3.0 x 2.4 pt heads) and 5.5 pt "UMAP1"/"UMAP2" labels.

    python figures/fig3b_strategies.py      # -> figures/fig3b_strategies.{svg,pdf,png}
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from svgkit import Canvas  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TOK = json.loads((ROOT / "tokens.json").read_text())
GEO = json.loads((ROOT / "figures/data/fig3b_geometry.json").read_text())
DOT = dict(zip(("blue", "red", "green"), TOK["color"]["schematic"]["cartoon_dots"]))
BLOB, SEP = TOK["color"]["neutral"]["box_fill"], "#BEBEBE"
FS_HEAD, FS_AXIS = TOK["type"]["size_pt"]["axis_title"], TOK["type"]["size_pt"]["tick_label"]

c = Canvas(ox=268, oy=164, w=236, h=106)


def arrow(x0, y0, x1, y1, stem=0.5, head_len=2.985, head_w=2.44, notch=0.71):
    """Stem from (x0, y0) to a notched head whose tip is (x1, y1): the paper's arrow (svgkit.arrow_poly
    generalised to any direction; head width 0.82 x length, notch 0.24 x length as in R's arrow_grob_sce2g)."""
    dx, dy = x1 - x0, y1 - y0
    L = (dx * dx + dy * dy) ** 0.5
    ux, uy = dx / L, dy / L
    px, py = -uy, ux
    bx, by = x1 - head_len * ux, y1 - head_len * uy            # head base centre
    nx, ny = x1 - (head_len - notch) * ux, y1 - (head_len - notch) * uy
    head = [(x1, y1), (bx + head_w / 2 * px, by + head_w / 2 * py), (nx, ny), (bx - head_w / 2 * px, by - head_w / 2 * py)]
    return [c.line(x0, y0, nx, ny, c.ink, stem), c.poly(head, c.ink)]


# ---- headings and separators ----------------------------------------------------------------
for cx, lines in ((302, ("(1) Single-cell", "one cell type")), (374, ("(2) Meta-cell", "one cell type")),
                  (456, ("(3) Single-cell", "several cell types"))):
    c.add(c.text(cx, 173.42, lines, FS_HEAD, "middle"))
for x in (336.0, 411.71):
    c.add(c.line(x, 164.0, x, 267.6, SEP, 0.5, dash=f"{c.S(2):.3f} {c.S(2):.3f}"))

# ---- meta-cell blobs (behind the dots), then the cells --------------------------------------
for b in GEO["blobs"]:
    c.add(c.path([tuple(seg) for seg in b], stroke="none", sw=0, fill=BLOB))
for d in GEO["dots"]:
    c.add(c.circle(d["x"], d["y"], 3.19, DOT[d["c"]]))

# ---- UMAP axes: corner at (278.9, 256), arrows 20.5 pt, labels starting at the corner ----------
c.add(*arrow(278.9, 256.0, 278.9, 235.19))
c.add(*arrow(278.9, 255.55, 299.27, 255.55))
c.add(c.text(278.75, 262.89, "UMAP1", FS_AXIS))
c.add(f'<text font-family="{c.text(0, 0, "").split(chr(34))[1]}" font-size="{c.S(FS_AXIS):.3f}" '
      f'transform="rotate(-90 {c.X(275.33):.3f} {c.Y(254.84):.3f})" x="{c.X(275.33):.3f}" y="{c.Y(254.84):.3f}">UMAP2</text>')

out = c.save(ROOT / "figures/fig3b_strategies.svg")
print(out, f"({c.S(c.w):.1f} × {c.S(c.h):.1f} mm)")
