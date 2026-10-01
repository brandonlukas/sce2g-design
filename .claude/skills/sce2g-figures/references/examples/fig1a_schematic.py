"""Fig 1a — the scE2G overview schematic, written directly as print-size SVG (1 user unit = 1 mm).

Three stages, left to right: (1) PREPARE single-cell multiome data (a cartoon embedding with
one selected cluster funnelling into an ATAC fragment file and an RNA count matrix), (2) ANNOTATE
peak–gene pairs (four feature boxes, each a small peak-track cartoon, plus the promoter-class
question), and (3) APPLY scE2G (the model box, its training recipe, and the CRISPR perturbation
data it was trained on).

Every coordinate below was measured from the vector content of the published PDF (page 3 of
refs/sce2g-2026.pdf) and is written in PDF points relative to the panel's top-left corner; the
emit helpers convert to millimetres. Colours and type sizes come from tokens.json.

    python figures/fig1a_schematic.py      # -> figures/fig1a_schematic.{svg,pdf,png}
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from svgkit import Canvas  # noqa: E402  (shared swept-arrow geometry)

ROOT = Path(__file__).resolve().parent
TOK = json.loads((ROOT / "tokens.json").read_text())
SCH, GEN, NEU = TOK["color"]["schematic"], TOK["color"]["genomic"], TOK["color"]["neutral"]

INK = NEU["text"]
BOX, BOX2, MUTED = SCH["container"], SCH["container_secondary"], SCH["glyph_muted"]
PEAK_L, PEAK_D, GENE, MODEL = SCH["peak_light"], SCH["peak_dark"], SCH["gene_arrow"], SCH["model_box"]
CLUSTER, DOTS, CONE = SCH["selected_cluster"], SCH["cluster_dots"], SCH["cone"]
PROT_D, PROT_L = SCH["protein_dark"], SCH["protein_light"]

FONT = "Helvetica, Arial, sans-serif"
MONO = "Courier, 'Courier New', monospace"
# type sizes measured in the typeset panel (pt): stage title, box label, annotation, mono table
FS_TITLE, FS_LABEL, FS_NOTE, FS_MONO, FS_SMALL = 6.0, 5.1, 4.2, 3.65, 4.65
# Glyph-size multiplier (line pitch is not scaled). Graphik Naturel is wider than Helvetica for
# lowercase but not for the caps-heavy "CRISPRi-" line, which overflows its box at 1.0; 0.97 fits.
TEXT_SCALE = 0.97

# ---- panel frame (PDF points, page coordinates from top-left); everything is relative to OX, OY
OX, OY = 125.0, 58.0
W_PT, H_PT = 350.0, 237.0
MM = 25.4 / 72


KIT = Canvas(OX, OY, W_PT, H_PT)   # same frame; used for the shared arrow polygon


def X(x):  # pt (page) -> mm (panel)
    return (x - OX) * MM


def Y(y):
    return (y - OY) * MM


def S(v):  # length pt -> mm
    return v * MM


# ---- emit helpers ----------------------------------------------------------------------------
def text(x, y, s, size=FS_NOTE, anchor="start", fill=INK, font=FONT, lead=1.2):
    """One or more lines (list) baseline-anchored at (x, y) pt; later lines step down by lead*size."""
    lines = s if isinstance(s, (list, tuple)) else [s]
    step = S(size) * lead        # baseline pitch from the nominal size
    size = size * TEXT_SCALE
    spans = "".join(
        f'<tspan x="{X(x):.3f}" y="{Y(y) + i * step:.3f}">{ln}</tspan>' for i, ln in enumerate(lines)
    )
    return (f'<text xml:space="preserve" font-family="{font}" font-size="{S(size):.3f}" text-anchor="{anchor}" '
            f'fill="{fill}">{spans}</text>')


def rect(x0, y0, x1, y1, fill, stroke="none", sw=0.0, opacity=None):
    op = f' fill-opacity="{opacity}"' if opacity is not None else ""
    return (
        f'<rect x="{X(x0):.3f}" y="{Y(y0):.3f}" width="{S(x1 - x0):.3f}" height="{S(y1 - y0):.3f}" '
        f'fill="{fill}" stroke="{stroke}" stroke-width="{S(sw):.3f}"{op}/>'
    )


def poly(pts, fill, stroke="none", sw=0.0):
    p = " ".join(f"{X(x):.3f},{Y(y):.3f}" for x, y in pts)
    return f'<polygon points="{p}" fill="{fill}" stroke="{stroke}" stroke-width="{S(sw):.3f}"/>'


def path(d_pt, stroke=INK, sw=0.23, fill="none", dash=None, cap="butt"):
    """d_pt: list of (cmd, *coords) in pt, e.g. ('M',x,y), ('L',x,y), ('C',x1,y1,x2,y2,x,y)."""
    parts = []
    for seg in d_pt:
        cmd, nums = seg[0], seg[1:]
        parts.append(cmd + " ".join(f"{X(v):.3f}" if i % 2 == 0 else f"{Y(v):.3f}" for i, v in enumerate(nums)))
    da = f' stroke-dasharray="{" ".join(f"{S(float(v)):.3f}" for v in dash.split())}"' if dash else ""   # dash given in pt
    return (
        f'<path d="{" ".join(parts)}" fill="{fill}" stroke="{stroke}" stroke-width="{S(sw):.3f}" '
        f'stroke-linecap="{cap}" stroke-linejoin="miter"{da}/>'
    )


def line(x0, y0, x1, y1, stroke=INK, sw=0.23, arrow=False, dash=None):
    return path([("M", x0, y0), ("L", x1, y1)], stroke, sw, dash=dash) if not arrow else path(
        [("M", x0, y0), ("L", x1, y1)], stroke, sw, dash=dash
    ).replace("/>", ' marker-end="url(#ah)"/>')


def polyline(pts, stroke=INK, sw=0.93, arrow=False):
    """Elbow connector. arrow=True adds the shared swept head (svgkit.arrow_poly, sized from the
    stroke) on the final, horizontal segment: measured tip = line end + 2.7 pt for a 0.93 pt line."""
    d = [("M", *pts[0])] + [("L", *p) for p in pts[1:]]
    out = path(d, stroke, sw)
    if arrow:
        (xa, ya), (xb, yb) = pts[-2], pts[-1]
        direction = 1 if xb >= xa else -1
        out += KIT.arrow_poly(xb - direction * 2.0, xb + direction * 2.9 * sw, yb, stroke, stem=sw)
    return out


def ellipse(cx, cy, rx, ry, fill, rot=0.0):
    return (
        f'<ellipse cx="{X(cx):.3f}" cy="{Y(cy):.3f}" rx="{S(rx):.3f}" ry="{S(ry):.3f}" fill="{fill}" '
        f'transform="rotate({rot} {X(cx):.3f} {Y(cy):.3f})"/>'
    )


def dot(cx, cy, fill, d=2.4):
    return f'<circle cx="{X(cx):.3f}" cy="{Y(cy):.3f}" r="{S(d / 2):.3f}" fill="{fill}"/>'


def peak(x0, x1, y_base, h, fill):
    """A smooth accessibility peak: a bell from (x0, y_base) to (x1, y_base) with apex h above."""
    xm, w = (x0 + x1) / 2, x1 - x0
    # concave flanks (first control far along the baseline, second close to the apex) as in the paper
    return path(
        [("M", x0, y_base), ("C", x0 + 0.40 * w, y_base, xm - 0.17 * w, y_base - h, xm, y_base - h),
         ("C", xm + 0.17 * w, y_base - h, x1 - 0.40 * w, y_base, x1, y_base), ("Z",)],
        stroke="none", fill=fill,
    )


def gene_arrow(x0, tip, y, colour=GENE, stem=1.85):
    """The paper's gene arrow: the shared swept arrow (svgkit.arrow_poly); head proportions follow
    from the stem weight (1.85 pt here, 3.7 pt for the promoter-class arrow, 0.92 pt for TSS glyphs)."""
    return KIT.arrow_poly(x0, tip, y, colour, stem=stem)


# ---- stage 1: prepare single-cell multiome data ------------------------------------------------
def stage1():
    out = [text(166.4, 65.6, ["(1) Prepare single-cell", "multiome data"], FS_TITLE, "middle")]
    # embedding axes (L-shape) and three clusters; the selected cluster is teal
    grey_dots = [(146.9, 86.3), (159.4, 87.5), (153.2, 88.7), (185.2, 88.7), (178.8, 89.2), (143.4, 92.6),
                 (161.6, 93.0), (149.1, 93.2), (154.8, 94.4), (176.4, 95.6), (182.7, 95.6), (189.6, 97.2),
                 (159.4, 98.4), (144.6, 99.0), (151.0, 99.0), (156.0, 102.4), (177.0, 102.4), (183.9, 102.4),
                 (149.1, 104.3), (189.6, 105.5), (183.9, 109.9)]
    teal_dots = [(163.5, 114.7), (154.8, 118.2), (170.8, 119.3), (161.1, 119.7), (175.8, 123.2),
                 (166.2, 124.4), (158.2, 125.5), (180.0, 125.5), (172.6, 129.5), (165.0, 130.7)]
    # cone from the selected cluster down to the fragment-file box (drawn first, under the blob)
    out.append(poly([(152, 122), (182, 122), (195.6, 158.5), (134.5, 158.5)], CONE))
    out.append(path([("M", 134.5, 76.9), ("L", 134.5, 138.1), ("L", 195.6, 138.1)], INK, 0.23))
    out.append(ellipse(152.7, 94.0, 13.8, 13.2, BOX2, rot=-15))
    out.append(ellipse(183.1, 98.3, 10.2, 16.3, BOX2, rot=12))
    out.append(ellipse(166.9, 122.3, 16.9, 11.2, CLUSTER, rot=-8))
    out += [dot(x, y, MUTED) for x, y in grey_dots]
    out += [dot(x, y, DOTS) for x, y in teal_dots]
    # data boxes
    out.append(rect(128.3, 158.5, 203.8, 186.3, BOX))
    out.append(text(166, 166.7, "ATAC-seq fragment file", FS_TITLE, "middle"))
    out.append(text(132.1, 172.6, ["Chr. 1 10011   10132  CGTAGCGG", "Chr. 1 10012   10151  ATGTCATC",
                                    "               ..."], FS_MONO, font=MONO))
    out.append(rect(128.3, 189.3, 203.8, 221.9, BOX))
    out.append(text(166, 195.8, "RNA count matrix", FS_TITLE, "middle"))
    out.append(text(130.1, 201.5, ["       CGTAGCGG   ATGTCATC     ...", "SAMD11         0          1    ...",
                                    "NOC2L          2          0    ...", "KLHL17         1          0    ...",
                                    "               ..."], FS_MONO, font=MONO))
    return out


def connectors():
    """Elbow connectors from the data boxes into the feature boxes. Drawn LAST so the arrowheads,
    whose tips overlap the box edge by ~1 pt as in the paper, sit on top of the grey boxes."""
    out = []
    # elbow connectors: two trunks, four arrowheads into the feature boxes
    out.append(polyline([(203.8, 166.6), (212, 166.6), (212, 146.2), (229.7, 146.2)], arrow=True))   # ATAC -> chromatin
    out.append(polyline([(203.8, 178.0), (220, 178.0)]))                                              # ATAC -> outer trunk
    out.append(polyline([(220, 193.5), (220, 102.4), (229.7, 102.4)], arrow=True))                    # outer trunk -> ABC
    out.append(polyline([(220, 193.5), (229.7, 193.5)], arrow=True))                                  # outer trunk -> gene density
    out.append(polyline([(203.8, 211.5), (212, 211.5), (212, 240.1), (229.7, 240.1)], arrow=True))    # RNA (+ATAC) -> kendall
    out.append(polyline([(212, 178.0), (212, 211.5)]))                                                # ATAC joins the kendall trunk
    return out


# ---- stage 2: annotate peak–gene pairs ---------------------------------------------------------
def track(y, promoter_box=False, activity_box=False, band=None):
    """The shared peak-track cartoon: four peaks (one dark), element bars beneath, gene arrow.
    band = (x0, y0, x1, y1) draws a grey annotation band under the element bars."""
    out = []
    if band:
        out.append(rect(*band, BOX2))
    if activity_box:
        out.append(rect(281.5, y - 16.6, 304.8, y + 5.3, BOX2))
    if promoter_box:
        out.append(rect(323.3, y - 11.5, 334.5, y + 4.0, BOX2))
    out += [peak(258.7, 278.2, y, 6.5, PEAK_L), peak(282.4, 304.2, y, 14.6, PEAK_D),
            peak(305.5, 318.5, y, 11.6, PEAK_L), peak(323.9, 334.0, y, 10.1, PEAK_L)]
    out.append(line(258.7, y, 332.1, y, INK, 0.23))
    out += [rect(258.8, y + 1.9, 277.9, y + 4.5, PEAK_L), rect(282.4, y + 1.9, 304.2, y + 4.5, PEAK_D),
            rect(305.5, y + 1.9, 318.5, y + 4.5, PEAK_L)]
    out.append(gene_arrow(328.9, 347.6, y + 0.3))
    return out


def stage2():
    out = [text(291.1, 65.6, "(2) Annotate peak–gene pairs", FS_TITLE, "middle")]
    # ABC score
    out.append(rect(231.3, 75.6, 350.7, 127.6, BOX))
    out.append(text(234.1, 83.9, "ABC score", FS_LABEL))
    out += track(103.0, activity_box=True)
    out.append(text(293.1, 85.0, "Activity", FS_NOTE, "middle"))
    out.append(text(257.2, 96.1, ["Pseudobulk", "ATAC"], FS_NOTE, "end"))
    # contact-vs-distance curve: a slow rise through the marked point (293.1, 119.7), steepening into a plateau
    out.append(path([("M", 258.8, 123.4), ("C", 275.0, 123.4, 286.0, 121.8, 293.1, 119.7),
                     ("C", 300.0, 117.4, 306.0, 113.4, 310.5, 111.2), ("L", 328.9, 111.2)], INK, 0.93))
    out.append(line(293.1, 103.0, 293.1, 119.7, INK, 0.23, dash="1.87 1.87"))
    out.append(dot(293.1, 119.7, INK, 1.8))
    out.append(text(257.2, 115.0, ["Contact ≈", "1/distance"], FS_NOTE, "end"))
    out.append(text(273.1, 118.2, "Contact", FS_NOTE))
    out.append(text(304.5, 120.5, "ABC score =", FS_NOTE))
    out.append(text(336.4, 118.2, "A × C", FS_NOTE, "middle"))
    out.append(line(328.8, 119.2, 343.9, 119.2, INK, 0.23))
    out.append(text(336.4, 123.7, "Σ(A × C)", FS_NOTE, "middle"))
    # Chromatin accessibility
    out.append(rect(231.3, 130.7, 350.7, 170.4, BOX))
    out.append(text(234.1, 137.8, ["Chromatin", "accessibility"], FS_LABEL))
    out += track(154.8, promoter_box=True, band=(266.1, 154.9, 320.0, 160.6))
    out.append(text(328.6, 136.9, ["ATAC signal", "at promoter"], FS_NOTE, "middle", lead=1.12))
    out.append(line(268.5, 162.4, 317.6, 162.4, INK, 0.46))      # "within 5 kb" span, heads both ways
    out.append(KIT.arrow_poly(270.0, 266.1, 162.4, INK, stem=0.46, head_len=2.8, head_w=2.2, notch=0.9))
    out.append(KIT.arrow_poly(316.1, 320.0, 162.4, INK, stem=0.46, head_len=2.8, head_w=2.2, notch=0.9))
    out.append(text(293.1, 167.7, "Number of peaks within 5 kb", FS_NOTE, "middle"))
    # Gene density
    out.append(rect(231.3, 173.7, 350.7, 211.2, BOX))
    out.append(text(234.1, 181.3, ["Gene", "density"], FS_LABEL))
    out += track(191.0, band=(293.5, 191.7, 329.0, 196.8))
    # two small TSS glyphs, one per strand: the same swept arrow at a quarter of the size
    out.append(gene_arrow(312.1, 318.6, 194.2, MUTED, stem=0.92))
    out.append(gene_arrow(327.1, 320.6, 194.2, MUTED, stem=0.92))
    out.append(path([("M", 293.6, 197.2), ("L", 293.6, 199.6), ("M", 328.8, 197.2), ("L", 328.8, 199.6),
                     ("M", 293.6, 198.4), ("L", 328.8, 198.4)], INK, 0.46))
    out.append(text(309.4, 204.0, ["Number of peaks, number of TSSs", "between peak and gene"], FS_NOTE, "middle", lead=1.12))
    # Kendall correlation
    out.append(rect(231.2, 214.5, 350.6, 264.1, BOX))
    out.append(text(234.1, 222.3, ["Kendall", "correlation"], FS_LABEL))
    out.append(text(304.3, 222.3, "Peak accessibility × gene expression", FS_NOTE, "middle"))
    for label, yb, pk, box, sq, n_sq in (
            ("Cell 1", 238.6, (280.1, 297.0, 9.2), (279.3, 228.5, 297.9, 240.4), (311.3, 229.4, 323.7, 234.9), 2),
            ("Cell 2", 256.8, (282.6, 293.4, 6.1), (281.7, 249.9, 294.0, 258.5), (311.3, 247.6, 323.7, 253.0), 1)):
        out.append(text(258.5, yb + 0.5, label, FS_NOTE))
        out.append(rect(*box, BOX2))
        out.append(peak(pk[0], pk[1], yb, pk[2], PEAK_D))
        out.append(line(272.2, yb, 315.8, yb, INK, 0.23))
        out.append(rect(*sq, BOX2))
        sx0, sy0, sx1, sy1 = sq
        ym = (sy0 + sy1) / 2
        if n_sq == 2:   # two transcripts: a short wave above a longer one
            out.append(path([("M", sx0 + 5.0, ym - 0.6), ("C", sx0 + 7.0, ym - 3.0, sx0 + 8.5, ym + 0.4, sx1 - 1.0, ym - 2.4)], GENE, 0.35))
            out.append(path([("M", sx0 + 1.2, ym + 1.6), ("C", sx0 + 4.0, ym - 1.0, sx0 + 6.5, ym + 3.2, sx0 + 9.5, ym + 0.6),
                             ("C", sx0 + 10.5, ym - 0.2, sx1 - 2.0, ym - 0.6, sx1 - 1.0, ym + 0.2)], GENE, 0.35))
        else:
            out.append(path([("M", sx0 + 3.0, ym - 0.2), ("C", sx0 + 4.5, ym + 2.6, sx0 + 6.5, ym + 1.0, sx0 + 8.0, ym - 0.4),
                             ("C", sx0 + 9.0, ym - 1.2, sx1 - 3.0, ym - 1.4, sx1 - 2.4, ym - 0.6)], GENE, 0.35))
        # dashed link: leaves the highlight box's top edge ~2/3 across, apex ~5 pt up, lands on the
        # expression box's top-left corner (apex ~4 pt up: clear of the label above)
        bx0, by0, bx1 = box[0], box[1], box[2]
        xa = bx0 + 0.66 * (bx1 - bx0)
        out.append(path([("M", xa, by0 - 0.2), ("C", xa + 4.0, by0 - 5.4, sx0 - 3.0, by0 - 5.6, sx0 + 1.0, sy0 - 0.2)],
                        INK, 0.23, dash="1.87 1.87"))
        out.append(gene_arrow(312.7, 331.4, yb + 0.1))
    # Promoter class
    out.append(rect(244.5, 266.7, 336.3, 290.8, BOX))
    out.append(text(290.2, 273.6, ["Promoter class", "ubiquitously", "expressed?"], FS_LABEL, "end"))
    out.append(gene_arrow(294.1, 328.6, 278.4, stem=3.7))
    # bracket gathering the stage, and the arrow into the model
    out.append(path([("M", 346.6, 72.8), ("L", 358.8, 72.8), ("L", 358.8, 293.1), ("L", 346.6, 293.1)], INK, 0.93))
    out.append(polyline([(358.8, 105.5), (392.9, 105.5)], arrow=True))
    return out


# ---- stage 3: apply scE2G ------------------------------------------------------------------------
def stage3():
    out = [text(420.6, 65.6, ["(3) Apply scE2G to make", "genome-wide enhancer–", "gene predictions"], FS_TITLE, "middle")]
    out.append(rect(395.5, 93.2, 444.5, 117.7, MODEL))
    out.append(text(420, 107.3, "scE2G", FS_TITLE, "middle", fill="#FFFFFF"))
    out.append(poly([(395.5, 117.7), (444.5, 117.7), (473.0, 134.4), (367.0, 134.4)], BOX))
    out.append(rect(367.4, 134.4, 473.0, 181.3, BOX))
    out.append(text(420.0, 127.9, "Model training", FS_TITLE, "middle"))
    out.append(text(375.0, 141.8, "(1)", FS_TITLE))
    out.append(text(385.7, 141.0, ["Compute features on K562", "multiome data"], FS_TITLE))
    out.append(text(375.0, 159.0, "(2)", FS_TITLE))
    out.append(text(385.7, 158.3, ["Train logistic regression", "classifier to predict CRISPRi-", "validated enhancers"], FS_TITLE))
    # CRISPR perturbation data
    out.append(rect(367.4, 186.5, 473.0, 249.5, BOX))
    out.append(text(420.2, 194.5, "CRISPR perturbation data", FS_TITLE, "middle"))
    # gRNA arcs first (the bean is drawn over their origins): leave the bean at ~30 deg, swing outward,
    # arrive steeply on the black bar; measured feet at 380.6 / 393.7 / 430.8 / 443.9
    for x_end in (380.6, 393.7, 430.8, 443.9):
        xs = 414.4 if x_end < 412 else 410.2
        dx = x_end - xs
        out.append(path([("M", xs, 206.5), ("C", xs + 0.58 * dx, 209.5, x_end - 0.15 * dx, 221.0, x_end, 233.3)],
                        INK, 0.92, dash="1.85 1.85"))
        out.append(rect(x_end - 2.4, 233.5, x_end + 2.4, 234.7, INK))                      # the targeted DHS
        # hairpin: thin tail left of the foot, loop just right of it (bbox ~8 x 6 pt)
        out.append(path([("M", x_end - 4.0, 231.6), ("C", x_end - 2.6, 231.7, x_end - 1.4, 231.5, x_end + 0.6, 230.8),
                         ("L", x_end + 1.6, 229.2), ("C", x_end + 1.9, 225.3, x_end + 5.0, 225.5, x_end + 4.6, 228.5),
                         ("C", x_end + 4.3, 230.5, x_end + 2.8, 231.6, x_end + 1.6, 231.6)], PROT_D, 0.28))
    # dCas9-KRAB: a light domain tucked behind the shoulder of a dark bean (bbox 403.6-420.3 x 200.1-211.1)
    out.append(ellipse(421.5, 201.3, 4.6, 3.6, PROT_L, rot=-25))
    out.append(path([("M", 404.0, 207.0), ("C", 403.0, 201.5, 409.5, 198.6, 414.5, 200.4),
                     ("C", 417.5, 201.5, 419.4, 204.2, 420.3, 207.2), ("C", 421.0, 210.3, 415.5, 212.0, 411.0, 211.0),
                     ("C", 407.0, 210.2, 404.6, 209.4, 404.0, 207.0), ("Z",)], stroke="none", fill=PROT_D))
    out.append(text(369.0, 204.9, "dCas9-KRAB", FS_LABEL))
    out.append(text(430.9, 204.6, ["466 positives", "9,876 negatives"], FS_LABEL))
    out.append(text(412.6, 222.3, ["gRNAs to", "DHS elements"], FS_SMALL, "middle"))
    out.append(line(376.5, 237.7, 449.8, 237.7, INK, 0.23))
    out.append(gene_arrow(406.0, 419.7, 237.8))
    for xq in (380.4, 393.1, 430.9, 443.1):
        out.append(text(xq, 245.1, "?", FS_TITLE, "middle"))
    return out


def build() -> str:
    w_mm, h_mm = S(W_PT), S(H_PT)
    defs = ""   # arrowheads are drawn as polygons (svgkit.arrow_poly); no SVG markers needed
    return "\n".join([
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w_mm:.2f}mm" height="{h_mm:.2f}mm" viewBox="0 0 {w_mm:.3f} {h_mm:.3f}">',
        defs,
        f'<rect width="{w_mm:.3f}" height="{h_mm:.3f}" fill="#FFFFFF"/>',
        *stage1(), *stage2(), *stage3(), *connectors(),
        "</svg>",
    ])


def main() -> None:
    out = ROOT / "fig1a_schematic.svg"
    out.write_text(build())
    ink = shutil.which("inkscape")
    if ink:
        for kind, extra in (("pdf", []), ("png", ["--export-dpi=600", "--export-background=white"])):
            subprocess.run([ink, str(out), f"--export-type={kind}", *extra, f"--export-filename={out.with_suffix('.' + kind)}"],
                           check=True, capture_output=True)
    print(f"{out} ({S(W_PT):.1f} × {S(H_PT):.1f} mm){'' if ink else '  [inkscape missing: no pdf/png]'}")


if __name__ == "__main__":
    main()
