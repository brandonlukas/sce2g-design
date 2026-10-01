# /// script
# dependencies = ["pandas", "scipy", "matplotlib"]
# ///
"""Violin + median point-interval plot of eQTL enrichment by biosample match category
(scE2G Fig. 5b style). Run: uv run eqtl_violin.py"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

HERE = Path(__file__).parent
df = pd.read_csv(HERE.parent.parent.parent.parent.parent / "sce2g-figures/evals/files/eqtl_enrichment.csv")

order = ["unmatched", "same_lineage", "same_cell_type"]
labels = {"unmatched": "Unmatched", "same_lineage": "Same\nlineage", "same_cell_type": "Same\ncell type"}
colors = {"unmatched": "#B0B0B0", "same_lineage": "#7FB2D8", "same_cell_type": "#2C6FAD"}
groups = [df.loc[df.category == c, "enrichment"].to_numpy() for c in order]

plt.rcParams.update({"font.size": 7, "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "pdf.fonttype": 42, "axes.linewidth": 0.5})
fig, ax = plt.subplots(figsize=(2.3, 2.4))

v = ax.violinplot(groups, positions=range(3), widths=0.8, showextrema=False)
for body, c in zip(v["bodies"], order):
    body.set_facecolor(colors[c]); body.set_edgecolor("none"); body.set_alpha(0.7)

# median point with IQR interval
for i, g in enumerate(groups):
    q1, med, q3 = pd.Series(g).quantile([0.25, 0.5, 0.75])
    ax.plot([i, i], [q1, q3], color="black", lw=1)
    ax.plot(i, med, "o", color="black", ms=3.5, zorder=3)

# Wilcoxon rank-sum (Mann-Whitney U) brackets
ymax = df.enrichment.max(); step = (ymax - df.enrichment.min()) * 0.09
for k, (a, b) in enumerate([(0, 1), (1, 2), (0, 2)]):
    p = mannwhitneyu(groups[a], groups[b], alternative="two-sided").pvalue
    y = ymax + step * (1 + 1.25 * k)
    ax.plot([a, a, b, b], [y - step * 0.2, y, y, y - step * 0.2], color="black", lw=0.6)
    ax.text((a + b) / 2, y + step * 0.1, f"P = {p:.1e}" if p < 1e-3 else f"P = {p:.3f}",
            ha="center", va="bottom", fontsize=6)

ax.set_xticks(range(3))
ax.set_xticklabels([f"{labels[c]}\n(n = {len(g)})" for c, g in zip(order, groups)])
ax.set_ylabel("eQTL enrichment")
ax.set_xlabel("Biosample match")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.tick_params(width=0.5, length=2)
fig.tight_layout()
fig.savefig(HERE / "eqtl_violin.pdf")
fig.savefig(HERE / "eqtl_violin.png", dpi=300)
