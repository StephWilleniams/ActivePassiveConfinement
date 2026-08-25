"""SM Figure S1(d): numerical-limiter engagement rates vs integration time step dt; the dashed line marks the production step dt = 6.25e-4 s.

Renders a simple version of this panel from figure_s1_d.csv alone (no processing).
Run: python3 figure_s1_d.py
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = "figure_s1_d"
df = pd.read_csv(HERE / f"{STEM}.csv")

STYLE = {"translational cap (passive)": ("tab:red", "o"),
         "translational cap (swimmer)": ("tab:orange", "s"),
         "force cap (swimmer)": ("tab:blue", "^"),
         "y-clamp (swimmer)": ("tab:green", "v")}

fig, ax = plt.subplots(figsize=(5.2, 4.0))
for series, g in df.groupby("series", sort=False):
    color, marker = STYLE[series]
    ax.plot(g.x, g.y, marker=marker, ms=5, lw=1.4, color=color, label=series)
ax.axvline(6.25e-4, color="black", ls="--", lw=1.5)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"integration time step, $dt$ (s)")
ax.set_ylabel("fraction of particle-steps")
ax.legend(fontsize=8, frameon=False, loc="lower right")

fig.tight_layout()
fig.savefig(HERE / f"{STEM}.png", dpi=200)
print(f"saved {STEM}.png")
