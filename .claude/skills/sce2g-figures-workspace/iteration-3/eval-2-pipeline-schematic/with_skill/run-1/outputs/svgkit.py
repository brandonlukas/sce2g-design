"""Minimal SVG emitter for measured schematic panels (1 user unit = 1 mm).

Coordinates are given in PDF points relative to a panel origin (so numbers measured from the
published PDF can be used verbatim); the canvas converts to millimetres on output. Colours and
type sizes are the caller's business (take them from tokens.json).

    c = Canvas(ox=88, oy=70, w=254, h=96)     # panel frame in page points
    c.add(c.text(100, 80, "hello", 5.5)); c.save("out.svg")
"""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

MM = 25.4 / 72
FONT = "Helvetica, Arial, sans-serif"
MONO = "Courier, 'Courier New', monospace"
# Glyph-size multiplier applied to text (line pitch is NOT scaled). Default 1.0: Helvetica at the
# nominal size is already ~10-15 % narrower than the journal's Graphik Naturel, so measured line
# breaks fit. Pass a smaller value per canvas only where a measured box would otherwise overflow.
TEXT_SCALE = 1.0


class Canvas:
    def __init__(self, ox: float, oy: float, w: float, h: float, ink: str = "#000000", text_scale: float = TEXT_SCALE):
        self.ox, self.oy, self.w, self.h, self.ink, self.text_scale = ox, oy, w, h, ink, text_scale
        self.items: list[str] = []

    # ---- coordinate helpers (pt -> mm) ---------------------------------------------------------
    def X(self, x): return (x - self.ox) * MM
    def Y(self, y): return (y - self.oy) * MM
    @staticmethod
    def S(v): return v * MM

    def add(self, *items):
        for it in items:
            self.items.extend(it if isinstance(it, list) else [it])

    # ---- primitives ----------------------------------------------------------------------------
    def text(self, x, y, s, size=5.5, anchor="start", fill=None, font=FONT, lead=1.2, style=""):
        lines = s if isinstance(s, (list, tuple)) else [s]
        step = self.S(size) * lead          # baseline step from the NOMINAL size (1.2 x 5.5 pt = 6.6 pt, as measured)
        size *= self.text_scale             # optional glyph shrink; line pitch is untouched
        spans = "".join(f'<tspan x="{self.X(x):.3f}" y="{self.Y(y) + i * step:.3f}">{ln}</tspan>'
                        for i, ln in enumerate(lines))
        st = f' style="{style}"' if style else ""
        return (f'<text xml:space="preserve" font-family="{font}" font-size="{self.S(size):.3f}" '
                f'text-anchor="{anchor}" fill="{fill or self.ink}"{st}>{spans}</text>')

    def rect(self, x0, y0, x1, y1, fill, stroke="none", sw=0.0):
        return (f'<rect x="{self.X(x0):.3f}" y="{self.Y(y0):.3f}" width="{self.S(x1 - x0):.3f}" '
                f'height="{self.S(y1 - y0):.3f}" fill="{fill}" stroke="{stroke}" stroke-width="{self.S(sw):.3f}"/>')

    def poly(self, pts, fill, stroke="none", sw=0.0):
        p = " ".join(f"{self.X(x):.3f},{self.Y(y):.3f}" for x, y in pts)
        return f'<polygon points="{p}" fill="{fill}" stroke="{stroke}" stroke-width="{self.S(sw):.3f}"/>'

    def path(self, d_pt, stroke=None, sw=0.25, fill="none", dash=None, cap="butt", marker_end=None, marker_start=None):
        """d_pt: list of (cmd, *coords) in pt, e.g. ('M',x,y), ('L',x,y), ('C',x1,y1,x2,y2,x,y), ('Z',)."""
        parts = []
        for seg in d_pt:
            cmd, nums = seg[0], seg[1:]
            parts.append(cmd + " ".join(f"{self.X(v):.3f}" if i % 2 == 0 else f"{self.Y(v):.3f}" for i, v in enumerate(nums)))
        extra = (f' stroke-dasharray="{dash}"' if dash else "") + \
                (f' marker-end="url(#{marker_end})"' if marker_end else "") + \
                (f' marker-start="url(#{marker_start})"' if marker_start else "")
        return (f'<path d="{" ".join(parts)}" fill="{fill}" stroke="{stroke or self.ink}" stroke-width="{self.S(sw):.3f}" '
                f'stroke-linecap="{cap}" stroke-linejoin="miter"{extra}/>')

    def line(self, x0, y0, x1, y1, stroke=None, sw=0.25, **kw):
        return self.path([("M", x0, y0), ("L", x1, y1)], stroke, sw, **kw)

    def polyline(self, pts, stroke=None, sw=0.93, **kw):
        return self.path([("M", *pts[0])] + [("L", *p) for p in pts[1:]], stroke, sw, **kw)

    def circle(self, cx, cy, d, fill, stroke="none", sw=0.0):
        return (f'<circle cx="{self.X(cx):.3f}" cy="{self.Y(cy):.3f}" r="{self.S(d / 2):.3f}" fill="{fill}" '
                f'stroke="{stroke}" stroke-width="{self.S(sw):.3f}"/>')

    def star(self, cx, cy, r, fill, stroke="none", sw=0.0, n=5, inner=0.45):
        import math
        pts = []
        for k in range(2 * n):
            rr = r if k % 2 == 0 else r * inner
            a = -math.pi / 2 + k * math.pi / n
            pts.append((cx + rr * math.cos(a), cy + rr * math.sin(a)))
        return self.poly(pts, fill, stroke, sw)

    def arc(self, x0, x1, y, h, stroke=None, sw=0.5):
        """Semi-elliptical arc from (x0, y) to (x1, y) rising h above."""
        return self.path([("M", x0, y), ("C", x0, y - 1.35 * h, x1, y - 1.35 * h, x1, y)], stroke, sw)

    def gene_arrow(self, x0, x1, y, colour, sw=1.85, head=(7.0, 5.7)):
        """Thick line with a filled triangular head whose base is at x1 (head extends past it)."""
        hw, hh = head
        return self.path([("M", x0, y), ("L", x1, y)], colour, sw) + \
            self.poly([(x1 - 0.2, y - hh / 2), (x1 + hw - 0.2, y), (x1 - 0.2, y + hh / 2)], colour)

    def arrow_poly(self, x0, tip, y, colour, stem=2.0, head_len=None, head_w=None, notch=None):
        """The paper's arrow, used for table glyphs, gene arrows and TSS glyphs alike: a stem running
        into a swept head whose back edge is notched toward the tip. Measured across Fig. 1a/2a the
        head scales with the stem — length 4 x stem, width 3.25 x stem, notch 1.2 x stem — so passing
        only `stem` reproduces every instance (1.85 pt gene arrows, 3.7 pt promoter arrow, 0.92 pt
        TSS glyphs). tip < x0 points left; notch=0 gives a plain triangle."""
        head_len = 4.0 * stem if head_len is None else head_len
        head_w = 3.25 * stem if head_w is None else head_w
        notch = 1.2 * stem if notch is None else notch
        d = 1 if tip >= x0 else -1
        xb = tip - d * head_len
        xn = xb + d * notch * (1 - stem / head_w)   # stem edges meet the notch edges: one seamless polygon
        return self.poly([(x0, y - stem / 2), (xn, y - stem / 2), (xb, y - head_w / 2), (tip, y),
                          (xb, y + head_w / 2), (xn, y + stem / 2), (x0, y + stem / 2)], colour)

    # ---- output --------------------------------------------------------------------------------
    def defs(self):
        return ("<defs>"
                f'<marker id="ah" viewBox="0 0 4 4" refX="3.6" refY="2" markerWidth="1.3" markerHeight="1.3" '
                f'markerUnits="userSpaceOnUse" orient="auto"><path d="M0,0 L4,2 L0,4 z" fill="{self.ink}"/></marker>'
                "</defs>")

    def svg(self, background="#FFFFFF"):
        w, h = self.S(self.w), self.S(self.h)
        bg = f'<rect width="{w:.3f}" height="{h:.3f}" fill="{background}"/>' if background else ""
        return "\n".join([f'<svg xmlns="http://www.w3.org/2000/svg" width="{w:.2f}mm" height="{h:.2f}mm" viewBox="0 0 {w:.3f} {h:.3f}">',
                          self.defs(), bg, *self.items, "</svg>"])

    def save(self, path: str | Path, png_dpi: int = 900) -> Path:
        out = Path(path); out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(self.svg())
        ink = shutil.which("inkscape")
        if ink:
            for kind, extra in (("pdf", []), ("png", [f"--export-dpi={png_dpi}", "--export-background=white"])):
                subprocess.run([ink, str(out), f"--export-type={kind}", *extra, f"--export-filename={out.with_suffix('.' + kind)}"],
                               check=True, capture_output=True)
        print(f"{out} ({self.S(self.w):.1f} × {self.S(self.h):.1f} mm){'' if ink else '  [inkscape missing: no pdf/png]'}")
        return out
