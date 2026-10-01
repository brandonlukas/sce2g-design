"""Fig 2a — "Overview of benchmarking data": a three-row table (CRISPR / eQTL / GWAS) with a
mini-schematic of each perturbation type, drawn as print-size SVG.

Geometry measured from page 5 of refs/sce2g-2026.pdf (points, page coordinates): outer frame
90–339.7 × 82.9–162.8, header rule at y = 98, row rules at 119.4 and 141.1, column rules at
x = 207.7 and 272.9. Table fill and hairline colours come from tokens.json (neutral.table_fill,
neutral.strong-ish #202020 rule at 0.5 pt); glyphs use neutral.mid grey.

    python figures/fig2a_table.py      # -> figures/fig2a_table.{svg,pdf,png}
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from svgkit import Canvas  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
TOK = json.loads((ROOT / "tokens.json").read_text())
NEU = TOK["color"]["neutral"]
FILL, RULE, GREY = NEU["table_fill"], "#202020", NEU["mid"]
FS_TITLE, FS_BODY = 6.5, 5.5

c = Canvas(ox=86, oy=70, w=258, h=96)

# ---- frame and rules -------------------------------------------------------------------------
X0, X1, Y0, Y1 = 90.0, 339.7, 82.9, 162.8
COLS, ROWS = (207.7, 272.9), (98.0, 119.4, 141.1)
c.add(c.rect(X0, Y0, X1, Y1, FILL, RULE, 0.5))
c.add(*[c.line(X0, y, X1, y, RULE, 0.5) for y in ROWS])
c.add(*[c.line(x, Y0, x, Y1, RULE, 0.5) for x in COLS])

# ---- title and header (all baselines/centres measured; text cells are left-aligned ~5 pt in from
# their column rule, first-column labels are centred over their glyph, line pitch 1.2 x size) ------
c.add(c.text((X0 + X1) / 2, 77.2, "Overview of benchmarking data", FS_TITLE, "middle"))
c.add(c.text(148.6, 92.8, "Genetic perturbation data", FS_BODY, "middle"))
c.add(c.text(240.2, 91.9, "Advantages", FS_BODY, "middle"))
c.add(c.text(306.4, 92.8, "Limitations", FS_BODY, "middle"))


# ---- the perturbation glyph: baseline, element marker, arc to a grey gene arrow -----------------
def glyph(y, marker):
    xe = 154.1   # element position on the baseline
    out = [c.line(138.0, y, 187.3, y, "#000000", 0.5),
           c.arc(xe, 170.5, y - 1.0, 6.0, GREY, 0.5),                       # lands on the stem's top-left corner
           c.arrow_poly(170.5, 179.7, y, GREY, stem=2.0, head_len=7.2, head_w=5.8, notch=1.8)]  # table variant: slightly slimmer head                              # measured swept arrow
    if marker == "circle":                                                # markers straddle the baseline
        out.append(c.circle(xe, y, 3.5, GREY))
    else:
        out.append(c.star(xe, y + 0.2, 2.7, GREY, "#000000", 0.21))
    return out


# CRISPR: element as a filled circle, hit by an inhibition bar ("CRISPRi or deletion")
c.add(c.text(94.1, 109.5, "CRISPR", FS_BODY))
c.add(*glyph(113.9, "circle"))
c.add(c.line(146.8, 103.4, 151.2, 111.1, "#000000", 0.5))                    # the perturbation stem
c.add(c.poly([(149.2, 111.0), (153.0, 109.4), (153.4, 110.4), (149.6, 112.0)], "#000000"))  # its bar
c.add(c.text(134.6, 104.0, "CRISPRi", FS_BODY, "middle"))
c.add(c.text(134.6, 110.3, "or deletion", FS_BODY, "middle"))
c.add(c.text(212.6, 106.3, ["Clear causation;", "inhibit entire element"], FS_BODY))
c.add(c.text(279.5, 110.1, "Mostly K562", FS_BODY))

# eQTL: fine-mapped eVariant (star) linked to its eGene
c.add(c.text(94.6, 132.0, "eQTL", FS_BODY))
c.add(*glyph(135.2, "star"))
c.add(c.text(139.4, 126.9, "Fine-mapped", FS_BODY, "middle"))
c.add(c.text(139.4, 133.1, "eVariant", FS_BODY, "middle"))
c.add(c.text(178.8, 133.0, "eGene", FS_BODY))
c.add(c.text(212.6, 127.7, ["Available in multiple", "cell types"], FS_BODY))
c.add(c.text(279.3, 128.2, ["Uncertainty in causal", "variants"], FS_BODY))

# GWAS: fine-mapped variant (star) linked to a known gene (independent coding variant)
c.add(c.text(94.2, 153.1, "GWAS", FS_BODY))
c.add(*glyph(153.7, "star"))
c.add(c.text(135.5, 146.2, "Fine-mapped", FS_BODY, "middle"))
c.add(c.text(135.5, 152.4, "variant", FS_BODY, "middle"))
c.add(c.text(170.2, 149.5, "Known gene*", FS_BODY))
c.add(c.text(124.7, 160.7, "*Independent coding variant", FS_BODY))
c.add(c.text(212.6, 152.6, "Application use case", FS_BODY))
c.add(c.text(279.3, 147.0, ["Uncertainty in causal", "variants cell types,", "and genes"], FS_BODY))

if __name__ == "__main__":
    c.save(ROOT / "figures" / "fig2a_table.svg")
