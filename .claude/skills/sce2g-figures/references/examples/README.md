# Reference implementations

Full worked examples behind the snippets in `../recipes.md`, each producing the
panel described there. Self-contained: `svgkit.py` and `tokens.json` are copied
in for the Python scripts; `data/` holds the digitised paper data every script
reads.

Run an R one directly:
    SCE2G_ROOT=../../assets Rscript fig2b_pr.R      # -> fig2b_pr.pdf/.png here

Run a Python one directly:
    python3 fig1a_schematic.py                       # -> fig1a_schematic.svg/.pdf/.png here

`SCE2G_ROOT` only matters for the R scripts (it tells `theme_sce2g.R` where to
find `tokens.json`); the Python scripts find their copy next to themselves with
no environment variable needed. Outputs land in this directory — move or
delete them after inspecting.
