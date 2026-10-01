"""Grade each run's outputs against the eval's assertions by measuring the PDF/SVG. Writes grading.json per run."""
import json, re, sys, glob, subprocess
from pathlib import Path
sys.path.insert(0, '/Users/brandonlukas/code/sce2g-design/.claude/skills/sce2g-figures/scripts')
from measure_svg import measure, to_svg
IT = Path(sys.argv[1])

def near(v, t, tol): return abs(v - t) <= tol
def svg_png_check(out, stem, width_mm):
    svg = out / f'{stem}.svg'; png = out / f'{stem}.png'
    dpi = 0
    if png.exists():
        try:
            from PIL import Image; w = Image.open(png).size[0]; dpi = w / (width_mm / 25.4)
        except Exception as e: dpi = -1
    return svg.exists() and dpi >= 800, f'svg={svg.exists()} png_dpi≈{dpi:.0f} (needs ≥800; assumed width {width_mm} mm)'
def n_glyph_rows(svg):
    ys = [round(float(m_.group(1)), 1) for m_ in re.finditer(r'<use data-text="n"[^>]*transform="matrix\([^,]+,[^,]+,[^,]+,[^,]+,[^,]+,([\d.]+)\)"', svg)]
    return ys
def stroke_widths(m): return {w: c for (w, col), c in m['strokes']}
def stroke_cols(m, w, tol=0.03): return {col.lower() for (ww, col), c in m['strokes'] if near(ww, w, tol)}
def texts(m): return [s for s, c in m['text_sizes'] if c > 0]

def grade_pr(out, m, svg):
    ex = out / 'pr_panel.pdf'; png = out / 'pr_panel.png'
    r = []
    r.append(('pr_panel.pdf and pr_panel.png both exist', ex.exists() and png.exists(), f'pdf={ex.exists()} png={png.exists()}'))
    ax = any(near(w, .25, .03) and col.lower() in ('#000000', 'black') for (w, col), c in m['strokes'])
    r.append(('Axis lines and ticks are 0.25 pt black in the PDF', ax, f"strokes={m['strokes'][:6]}"))
    pal = {'#792374', '#006479', '#0096a0', '#404040'}
    curves = {col.lower() for (w, col), c in m['strokes'] if near(w, .42, .04)}
    r.append(("Curves are 0.42 pt strokes in the fixed predictor colours (#792374, #006479, #0096A0, #404040)", pal <= curves, f'0.42 pt colours={sorted(curves)}'))
    dash_ok = any(near(w or 0, .5, .05) and 1.8 <= float(str(da).replace(',', ' ').split()[0]) <= 2.4 for (da, w), c in m['dashes'])
    r.append(('Reference lines are 0.5 pt with a 2 to 2.25 pt dash', dash_ok, f"dashes={m['dashes']}"))
    ts = texts(m); t_ok = bool(ts) and any(near(s, 5.5, .1) for s in ts) and any(near(s, 6.5, .1) for s in ts) and min(ts) >= 4.4
    r.append(('Text sizes are 5.5 pt tick labels and 6.5 pt titles, nothing below 4.5 pt', t_ok, f'sizes={ts}'))
    # index labels in curve colours: text fills — count <use> elements whose enclosing fill is a palette colour
    coloured = len(re.findall(r'<use data-text="[^"]*"[^>]*fill="(#792374|#006479|#0096a0|#404040)"', svg, re.I))
    r.append(("Curve index numbers are drawn at the curve ends in the curve's own colour", coloured >= 3, f'coloured glyphs={coloured}'))
    ok, ev = svg_png_check(out, 'pr_panel', 62); r.append(('An SVG is saved alongside and the PNG is at least 800 dpi', ok, ev))
    has_n = bool(re.search(r'data-text="n"', svg)) and bool(re.search(r'data-text="="', svg)) and bool(re.search(r'data-text="\("', svg))
    r.append(('The recall axis title carries the number of positives, "Recall (n = …)"', has_n, f'n/=/( glyphs present: {has_n}'))
    return r

def grade_violin(out, m, svg):
    ex = out / 'eqtl_violin.pdf'; png = out / 'eqtl_violin.png'; r = []
    r.append(('eqtl_violin.pdf and eqtl_violin.png both exist', ex.exists() and png.exists(), f'pdf={ex.exists()} png={png.exists()}'))
    fills = {f.lower() for f, c in m['top_fills']}
    allf = {f.lower() for f in re.findall(r' fill="(#[0-9a-fA-F]{6})"', svg)}
    want = {'#8b98ac', '#9a3983', '#260032'}
    outlined = any(col.lower() in want for (w, col), c in m['strokes'] if w > 0.3)
    r.append(('Violins are filled #8B98AC / #9A3983 / #260032 with no outline stroke', want <= allf and not outlined, f'fills∩={sorted(want & allf)} outlined={outlined}'))
    ws = stroke_widths(m); thick = any(near(w, 1.79, .25) for w in ws); thin = any(near(w, .77, .12) for w in ws)
    r.append(('Point-interval lines are two weights, about 0.77 pt and 1.79 pt, in black', thick and thin, f'widths={sorted(ws)}'))
    br = any(near(w, .25, .03) for w in ws); stars = any(7.5 <= s <= 9 for s in texts(m)); pv = any(near(s, 4.5, .15) for s in texts(m))
    r.append(('Brackets are 0.25 pt lines with stars below and P values above', br and stars and pv, f'0.25={br} star-size={stars} 4.5pt={pv} sizes={texts(m)}'))
    ts = texts(m); t_ok = bool(ts) and any(near(s, 5.5, .1) for s in ts) and any(near(s, 6.5, .1) for s in ts) and min(ts) >= 3.0
    body=[x for x in ts if x >= 4.0]   # plotmath superscripts (0.7 x 4.5 = 3.15 pt) are allowed
    r.append(('Axis lines are 0.25 pt and text sizes include 5.5 and 6.5 pt with nothing below 4.5 pt (superscripts excepted)', br and t_ok and (not body or min(body) >= 4.4), f'sizes={ts}'))
    italic_n = bool(re.search(r'data-text="n"', svg)) and bool(re.search(r'data-text="="', svg))
    r.append(('Category labels with italic n appear under the axis', italic_n, f'n-glyph and = present: {italic_n}'))
    ok, ev = svg_png_check(out, 'eqtl_violin', 62); r.append(('An SVG is saved alongside and the PNG is at least 800 dpi', ok, ev))
    ys = n_glyph_rows(svg); rows = len(set(ys))
    r.append(('Each n line sits directly under its own label (n rows are not all on one shared baseline)', len(ys) >= 3 and rows >= 2, f'n-glyph baselines={sorted(set(ys))}'))
    return r

def grade_schematic(out, m, svg):
    r = []; f = out / 'pipeline_schematic.svg'
    r.append(('pipeline_schematic.svg exists (pdf/png optional)', f.exists(), f'svg={f.exists()}'))
    s = f.read_text() if f.exists() else ''
    cont = bool(re.search(r'fill="#dedfe1"', s, re.I)); box = bool(re.search(r'#9c4787', s, re.I))
    r.append(('Containers use #DEDFE1 without a border and the model box is #9C4787', cont and box, f'container={cont} model box={box}'))
    gene = bool(re.search(r'#54ab58', s, re.I)); markers = 'marker-end=' in s or 'marker-start=' in s
    polys = len(re.findall(r'<polygon', s)) + len(re.findall(r'<path[^>]*Z"[^>]*fill="#', s))
    r.append(('The gene arrow is #54AB58 and arrowheads are filled polygons (the notched shape), not marker triangles', gene and polys >= 2 and not markers, f'gene={gene} polygons={polys} markers={markers}'))
    widths = {round(float(w), 2) for w in re.findall(r'stroke-width="([\d.]+)"', s) if float(w) > 0}
    mm = 25.4 / 72; fam = {0.25, 0.9, 0.93, 1.85, 0.23, 0.5}
    pt_widths = {round(w / mm, 2) for w in widths} if widths and max(widths) < 1 else widths
    clusters = []
    for w in sorted(pt_widths):
        if not clusters or w - clusters[-1][-1] > 0.05: clusters.append([w])
        else: clusters[-1].append(w)
    ok = 0 < len(clusters) <= 3 and all(any(abs(w - t) <= 0.08 for t in fam) for w in pt_widths)
    r.append(('At most three distinct stroke widths are used, from the 0.25 / 0.9 / 1.85 pt family', ok, f'widths(mm or pt)={sorted(widths)} -> pt={sorted(pt_widths)}'))
    heads = bool(re.search(r'\(1\)', s)) and bool(re.search(r'\(2\)', s)) and bool(re.search(r'\(3\)', s))
    sizes = [float(v) for v in re.findall(r'font-size="([\d.]+)', s)]
    six = any(abs(v - 6.5) <= .2 or abs(v - 6.5 * mm) <= .05 for v in sizes)
    helv = bool(re.search(r'Helvetica|Arial', s))
    r.append(("Stage headings are numbered '(1)', '(2)', '(3)' at about 6.5 pt in Helvetica/Arial", heads and six and helv, f'numbered={heads} 6.5pt={six} helv={helv} sizes={sorted(set(sizes))[:8]}'))
    # cartoon cells (pure-RGB dots) must not sit inside a grey container
    rects = [tuple(map(float, m_.groups())) for m_ in re.finditer(r'<rect x="([\d.]+)" y="([\d.]+)" width="([\d.]+)" height="([\d.]+)" fill="#dedfe1"', s, re.I)]
    cells = [tuple(map(float, m_.groups())) for m_ in re.finditer(r'<circle cx="([\d.]+)" cy="([\d.]+)" r="[\d.]+" fill="(?:#1f84ff|#f40014|#00bf2b)"', s, re.I)]
    inside = sum(1 for cx, cy in cells for x, y, w, h in rects if x <= cx <= x + w and y <= cy <= y + h)
    r.append(('The cartoon embedding (cell dots) sits on white, not inside a grey container', bool(cells) and inside == 0, f'cells={len(cells)} inside grey rects={inside}'))
    ok, ev = svg_png_check(out, 'pipeline_schematic', 120); r.append(('The PNG is at least 800 dpi', ok, ev))
    # connector heads: black polygons must not intrude >1 pt into a grey container, and any container they
    # touch must be drawn before them (document order) so nothing paints over the head
    mm = 25.4 / 72; bad = []; heads = 0
    rect_pos = [(m_.start(), tuple(map(float, m_.groups()))) for m_ in re.finditer(r'<rect x="([\d.]+)" y="([\d.]+)" width="([\d.]+)" height="([\d.]+)" fill="#dedfe1"', s, re.I)]
    for m_ in re.finditer(r'<polygon points="([^"]*)" fill="#000000"', s, re.I):
        pts = [tuple(map(float, q.split(','))) for q in m_.group(1).split()]
        if len(pts) not in (4, 7): continue
        heads += 1; xs = [q[0] for q in pts]; ys = [q[1] for q in pts]
        for pos, (x, y, w, h) in rect_pos:
            ox = min(max(xs), x + w) - max(min(xs), x); oy = min(max(ys), y + h) - max(min(ys), y)
            if ox > 0 and oy > 0:
                if min(ox, oy) > 1.0 * mm: bad.append(f'intrudes {min(ox, oy) / mm:.1f} pt')
                if pos > m_.start(): bad.append('container drawn over head')
    r.append(('Connector heads end on the container edge (no intrusion > 1 pt) and are drawn above the containers', heads > 0 and not bad, f'heads={heads} issues={bad[:4]}'))
    return r

for ev in sorted(IT.glob('eval-*')):
    meta = json.load(open(ev / 'eval_metadata.json'))
    for run in ('with_skill', 'without_skill'):
        out = ev / run / 'run-1' / 'outputs'
        if not out.exists(): continue
        pdfs = [p for p in out.glob('*.pdf')]
        m = {'strokes': [], 'dashes': [], 'text_sizes': [], 'discs': [], 'fills': 0, 'top_fills': []}; svg = ''
        if pdfs:
            try: svg = to_svg(pdfs[0], 1); m = measure(svg)
            except Exception as e: print('measure failed', pdfs[0], e)
        g = {0: grade_pr, 1: grade_violin, 2: grade_schematic}[meta['eval_id']](out, m, svg)
        res = {'expectations': [{'text': t, 'passed': bool(p), 'evidence': str(e)[:300]} for t, p, e in g]}
        res['pass_rate'] = sum(x['passed'] for x in res['expectations']) / len(res['expectations'])
        n=len(res['expectations']); k=sum(x['passed'] for x in res['expectations']); res['summary']={'pass_rate':k/n,'passed':k,'failed':n-k,'total':n}
        json.dump(res, open(ev / run / 'run-1' / 'grading.json', 'w'), indent=2)
        print(ev.name, run, f"{res['pass_rate']:.2f}", [x['passed'] for x in res['expectations']])
