import csv
from collections import defaultdict
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

CSV = "/Users/brandonlukas/code/sce2g-design/.claude/skills/sce2g-figures/evals/files/pr_curves.csv"
OUT = "/Users/brandonlukas/code/sce2g-design/.claude/skills/sce2g-figures-workspace/iteration-1/eval-0-pr-curve-panel/without_skill/outputs/pr_panel"
OPERATING_RECALL, CHANCE = 0.7, 0.06

curves = defaultdict(list)
with open(CSV) as f:
    for row in csv.DictReader(f):
        curves[int(row["predictor_id"])].append((float(row["recall"]), float(row["precision"])))

plt.rcParams.update({"font.family": "Helvetica", "font.size": 7, "pdf.fonttype": 42,
                     "axes.linewidth": 0.5, "xtick.major.width": 0.5, "ytick.major.width": 0.5})
fig, ax = plt.subplots(figsize=(2.2, 2.2))
colors = {1: "#D62728", 2: "#1F77B4", 5: "#2CA02C", 15: "#7F7F7F"}
for pid in sorted(curves):
    pts = sorted(curves[pid])
    r, p = zip(*pts)
    ax.plot(r, p, lw=1, color=colors.get(pid, "black"))
# labels at right end, spread apart so they don't overlap (curves converge near chance)
ends = sorted((sorted(curves[pid])[-1][1], pid) for pid in curves)
GAP, y = 0.06, -1
for p_end, pid in ends:
    y = max(p_end, y + GAP)
    ax.annotate(str(pid), (1.0, y), xytext=(3, 0), textcoords="offset points", va="center", ha="left",
                fontsize=6, color=colors.get(pid, "black"))
ax.axvline(OPERATING_RECALL, ls="--", lw=0.5, color="black")
ax.axhline(CHANCE, ls="--", lw=0.5, color="grey")
ax.set(xlim=(0, 1.05), ylim=(0, 1.02), xlabel="Recall", ylabel="Precision", title="CRISPR benchmark (training)")
ax.set_aspect("equal")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(OUT + ".pdf")
fig.savefig(OUT + ".png", dpi=300)
