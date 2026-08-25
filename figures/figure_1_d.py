"""Figure 1(d): stationary distributions of both species, no-kick baseline (k = 0).

Renders a simple version of this panel from figure_1_d.csv alone (no processing).
Run: python3 figure_1_d.py
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = "figure_1_d"
df = pd.read_csv(HERE / f"{STEM}.csv")

COLORS = {"passive": "tab:orange", "swimmer": "tab:blue"}
BIN_WIDTH = 0.03           # histogram bin width h from the caption

fig, ax = plt.subplots(figsize=(5.5, 4.0))
for sp in ("passive", "swimmer"):
    h = df[df.series == f"{sp}_hist"]
    ax.bar(h.x, h.y, width=BIN_WIDTH, color=COLORS[sp], alpha=0.35,
           edgecolor="none", label=f"{sp} (histogram)")
    k = df[df.series == f"{sp}_kde"]
    ax.plot(k.x, k.y, color=COLORS[sp], lw=1.8, label=f"{sp} (KDE)")
for edge in (-1.5, 1.5):   # inner edge of the kick band, both walls
    ax.axvline(edge, ls="--", lw=1.2, color="0.45", zorder=0)
ax.set_xlabel("$y$")
ax.set_ylabel("$p(y)$")
ax.set_xlim(-3.0, 3.0)
ax.set_ylim(bottom=0)
ax.legend(fontsize=8, frameon=False)

fig.tight_layout()
fig.savefig(HERE / f"{STEM}.png", dpi=200)
print(f"saved {STEM}.png")
