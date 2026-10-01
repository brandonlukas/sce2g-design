#!/usr/bin/env python3
"""Measure the vector content of a PDF page (or an SVG): stroke widths, dash arrays, glyph sizes and
disc diameters, optionally inside a bounding box in PDF points (origin top-left, y down).

    python measure_svg.py panel.pdf
    python measure_svg.py paper.pdf --page 7 --bbox 363 615 567 740

Needs `mutool` (MuPDF) for PDF input. Output units are PDF points.
"""
import argparse, collections, math, re, subprocess, sys, tempfile
from pathlib import Path

def to_svg(path, page):
    if path.suffix.lower() == ".svg":
        return path.read_text()
    out = Path(tempfile.mkdtemp()) / "page.svg"
    subprocess.run(["mutool", "draw", "-o", str(out), str(path), str(page)], check=True, capture_output=True)
    return out.read_text()

def measure(svg, bbox=None):
    H = None
    m = re.search(r'<svg[^>]*height="([\d.]+)"', svg)
    if m: H = float(m.group(1))
    def inbox(x, y): return bbox is None or (bbox[0] <= x <= bbox[2] and bbox[1] <= y <= bbox[3])
    strokes, dashes, fills, discs, texts = (collections.Counter() for _ in range(5))
    for pm in re.finditer(r"<path[^>]*>", svg):
        p = pm.group(0); d = re.search(r'd="([^"]*)"', p).group(1)
        nums = [float(v) for v in re.findall(r"-?[\d.]+", d)]
        if not nums: continue
        t = re.search(r'transform="matrix\(([^)]*)\)"', p)
        if t:
            a = [float(v) for v in t.group(1).split(",")]
            page_abs = a[4] == 0 and H and abs(a[5] - H) < 0.01
            x, y = (nums[0], H - nums[1]) if page_abs else (a[4], a[5])
        else:
            x, y = nums[0], nums[1]
        if not inbox(x, y): continue
        g = lambda k: (re.search(k + r'="([^"]*)"', p) or [None, None])[1]
        sw, st, f, da = g("stroke-width"), g(" stroke"), g(" fill"), g("stroke-dasharray")
        if sw and st and st != "none": strokes[(round(float(sw), 2), st)] += 1
        if da: dashes[(da, round(float(sw), 2) if sw else None)] += 1
        if f and f != "none": fills[f] += 1
        if "C" in d and 8 <= len(nums) <= 40:
            xs, ys = nums[0::2], nums[1::2]; w = max(xs) - min(xs); h = max(ys) - min(ys)
            if len(nums) == 14 and not t: h *= 0.75          # cairo: 2 Beziers, control points overshoot
            if w > 0 and abs(w - h) / w < 0.15: discs[(round(w, 2), f or "default")] += 1
    for tm in re.finditer(r'<use data-text="[^"]*"[^>]*transform="matrix\(([^)]*)\)"', svg):
        a = [float(v) for v in tm.group(1).split(",")]
        if inbox(a[4], a[5]): texts[round(math.hypot(a[0], a[1]), 2)] += 1
    return {"strokes": strokes.most_common(15), "dashes": dashes.most_common(8), "text_sizes": sorted(texts.items()),
            "discs": discs.most_common(10), "fills": len(fills), "top_fills": fills.most_common(8)}

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("file"); ap.add_argument("--page", type=int, default=1)
    ap.add_argument("--bbox", type=float, nargs=4, metavar=("X0", "Y0", "X1", "Y1"))
    a = ap.parse_args()
    for k, v in measure(to_svg(Path(a.file), a.page), a.bbox).items(): print(k, v)
