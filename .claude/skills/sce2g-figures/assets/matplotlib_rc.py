"""matplotlib fallback for the scE2G style. `from matplotlib_rc import apply, TOK, PRED; apply()`.
Sizes in pt at final print size; the open-axes / outward-tick / no-grid conventions of theme_sce2g()."""
import json
from pathlib import Path
import matplotlib as mpl

TOK = json.loads((Path(__file__).with_name("tokens.json")).read_text())
PRED = {k: v["hex"] for k, v in TOK["color"]["predictor"].items()}      # index -> hex
S, W = TOK["type"]["size_pt"], TOK["stroke_pt"]

def apply():
    mpl.rcParams.update({
        "font.family": "Helvetica", "font.size": S["axis_title"],
        "axes.titlesize": S["panel_title"], "axes.labelsize": S["axis_title"],
        "xtick.labelsize": S["tick_label"], "ytick.labelsize": S["tick_label"], "legend.fontsize": S["legend"],
        "axes.linewidth": W["axis"], "xtick.major.width": W["tick"], "ytick.major.width": W["tick"],
        "xtick.major.size": TOK["mark_pt"]["tick_length"], "ytick.major.size": TOK["mark_pt"]["tick_length"],
        "xtick.direction": "out", "ytick.direction": "out",
        "axes.spines.top": False, "axes.spines.right": False, "axes.grid": False,
        "lines.linewidth": W["curve"], "lines.markersize": TOK["mark_pt"]["threshold_dot_diameter"],
        "legend.frameon": False, "figure.dpi": 300, "savefig.dpi": 300,
        "pdf.fonttype": 42, "svg.fonttype": "none",
    })

def dashed(ax_line, kind="reference"):
    """Apply the paper's dash pattern (2/2 pt reference, 2.52/0.84 leader) to a Line2D."""
    on, off = TOK["mark_pt"]["dash"][kind]
    ax_line.set_linestyle((0, (on / ax_line.get_linewidth(), off / ax_line.get_linewidth())))
