"""Conceptual pipeline schematic (no data). 120 x 45 mm; SVG + PDF + PNG."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({"font.family": "Helvetica", "svg.fonttype": "none", "pdf.fonttype": 42})
MM = 1 / 25.4
fig = plt.figure(figsize=(120 * MM, 45 * MM))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 120); ax.set_ylim(0, 45); ax.axis("off")

GREY, DARK, ACCENT, PEAK = "#bfbfbf", "#333333", "#d95f02", "#1b9e77"
panels = [(2, "Single-cell multiome"), (42, "Candidate E-G pairs"), (82, "Model prediction")]
for i, (x0, title) in enumerate(panels):
    ax.add_patch(FancyBboxPatch((x0, 4), 36, 34, boxstyle="round,pad=0,rounding_size=2",
                                fc="white", ec=GREY, lw=0.8))
    ax.add_patch(plt.Circle((x0 + 4, 34), 2.6, fc=DARK, ec="none"))
    ax.text(x0 + 4, 34, str(i + 1), ha="center", va="center", color="white", fontsize=7, weight="bold")
    ax.text(x0 + 8, 34, title, ha="left", va="center", fontsize=6.5, color=DARK)
    if i < 2:  # connector arrow to next panel
        ax.add_patch(FancyArrowPatch((x0 + 36.5, 21), (x0 + 41.5, 21), arrowstyle="-|>",
                                     mutation_scale=8, color=DARK, lw=1))

# (1) cartoon UMAP: three blobs of cells
rng = np.random.default_rng(0)
for cx, cy, col in [(13, 21, "#7570b3"), (25, 24, "#1b9e77"), (22, 13, "#e7298a")]:
    pts = np.clip(rng.normal([cx, cy], [2.4, 1.8], (60, 2)), [7, 9], [36, 30])
    ax.scatter(pts[:, 0], pts[:, 1], s=3, c=col, ec="none", alpha=0.9)
ax.annotate("", (12, 8), (5, 8), arrowprops=dict(arrowstyle="-|>", lw=0.6, color=GREY, mutation_scale=6))
ax.annotate("", (5, 15), (5, 8), arrowprops=dict(arrowstyle="-|>", lw=0.6, color=GREY, mutation_scale=6))
ax.text(5, 6.5, "UMAP", fontsize=4.5, color=GREY, va="top")
ax.text(34, 7, "RNA + ATAC", fontsize=5, color=DARK, ha="right")

# (2) signal peak (enhancer) + gene arrow, joined by an arc
xs = np.linspace(46, 60, 100); ys = 12 + 12 * np.exp(-((xs - 51) / 2.2) ** 2)
ax.fill_between(xs, 12, ys, color=PEAK, alpha=0.85, lw=0)
ax.plot([45, 75], [12, 12], color=DARK, lw=0.8)                     # genomic axis
ax.plot([66, 66], [12, 15.5], color=DARK, lw=0.8)                   # TSS
ax.annotate("", (72, 15.5), (66, 15.5), arrowprops=dict(arrowstyle="-|>", lw=0.8, color=DARK, mutation_scale=7))
ax.add_patch(FancyBboxPatch((66, 10.8), 8, 2.4, boxstyle="square,pad=0", fc=DARK, ec="none"))
ax.text(51, 9.5, "enhancer", fontsize=5, ha="center", va="top", color=PEAK)
ax.text(70, 9.5, "gene", fontsize=5, ha="center", va="top", color=DARK)
ax.annotate("", (67, 17), (51, 24.5), arrowprops=dict(arrowstyle="-", connectionstyle="arc3,rad=-0.5",
                                                     ls="--", lw=0.8, color=ACCENT))
ax.text(59, 27.5, "E-G pair", fontsize=5, ha="center", color=ACCENT)

# (3) highlighted model box with inputs -> score
ax.add_patch(FancyBboxPatch((93, 15), 16, 12, boxstyle="round,pad=0,rounding_size=1.5",
                            fc="#fde9d9", ec=ACCENT, lw=1.2))
ax.text(101, 21, "scE2G\nmodel", ha="center", va="center", fontsize=6, color=ACCENT, weight="bold")
for y, lab in [(24.5, "ATAC"), (21, "Kendall"), (17.5, "distance")]:
    ax.text(91, y, lab, ha="right", va="center", fontsize=4.5, color=DARK)
    ax.annotate("", (93, y), (91.5, y), arrowprops=dict(arrowstyle="-", lw=0.6, color=DARK))
ax.annotate("", (113, 21), (109.5, 21), arrowprops=dict(arrowstyle="-|>", lw=0.8, color=DARK, mutation_scale=7))
ax.text(113.5, 21, "score", ha="left", va="center", fontsize=5, color=DARK)

for ext in ("svg", "pdf", "png"):
    fig.savefig(f"pipeline_schematic.{ext}", dpi=300)
