"""Figure 3(c): passive accumulation E^p vs kick rate k, coloured by kick-strength tercile.

Renders a simple version of this panel from figure_3_c.csv alone (no processing).
Run: python3 figure_3_c.py
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = "figure_3_c"
df = pd.read_csv(HERE / f"{STEM}.csv")

COLORS = ["#4C72B0", "#DD8452", "#55A868"]

fig, ax = plt.subplots(figsize=(5.6, 4.0))
for (i, label), g in df.groupby(["tercile_index", "tercile_label"]):
    ax.scatter(g.x, g.E_p, s=16, alpha=0.75, color=COLORS[i], linewidths=0,
               label=label)
ax.set_xlabel(r"$k$ (s$^{-1}$)")
ax.set_ylabel(r"$E^{\,p}$")
ax.grid(True, color="0.85", lw=0.5, alpha=0.6)
ax.set_axisbelow(True)
ax.legend(title=r"$\tau$ $\in$", fontsize=8, ncol=3, loc="lower center",
          bbox_to_anchor=(0.5, 1.0), frameon=False)

fig.tight_layout()
fig.savefig(HERE / f"{STEM}.png", dpi=200)
print(f"saved {STEM}.png")
