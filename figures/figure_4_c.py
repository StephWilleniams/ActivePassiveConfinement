"""Figure 4(c): local wall-normal diffusivity D_y(y_eff) at the four kick-parameter corners; faint raw per-bin curves under a LOESS trend.

Renders a simple version of this panel from figure_4_c.csv alone (no processing).
Run: python3 figure_4_c.py
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = "figure_4_c"
df = pd.read_csv(HERE / f"{STEM}.csv")

CORNER_COLORS = {"LR LS": "#4C72B0", "LR HS": "#55A868",
                 "HR LS": "#C44E52", "HR HS": "#8172B2"}

fig, ax = plt.subplots(figsize=(5.8, 4.2))
for corner, color in CORNER_COLORS.items():
    raw = df[(df.corner == corner) & (df.kind == "raw")]
    ax.plot(raw.x, raw.y, lw=1.0, color=color, alpha=0.30)
    tr = df[(df.corner == corner) & (df.kind == "trend")]
    ax.plot(tr.x, tr.y, lw=2.0, color=color, label=corner)
# wall furniture: wall line at y_eff = -0.5, unvisited strip shaded
ax.axvline(-0.5, color="#374151", lw=2.5)
ax.axvspan(-0.5, 0.0, color="#9CA3AF", alpha=0.22, lw=0)
ax.set_xlim(-0.55, 2.5)
ax.set_xlabel(r"$y_{\mathrm{eff}}$")
ax.set_ylabel(r"$\langle D_y(y_{\mathrm{eff}})\rangle$")
ax.legend(fontsize=8, ncol=2)

fig.tight_layout()
fig.savefig(HERE / f"{STEM}.png", dpi=200)
print(f"saved {STEM}.png")
