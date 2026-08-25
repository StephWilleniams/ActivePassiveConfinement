"""Figure 4(d): LRLS passive stationary distribution -- measured KDE vs the MTM and jump-KM reconstructions.

Renders a simple version of this panel from figure_4_d.csv alone (no processing).
Run: python3 figure_4_d.py
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = "figure_4_d"
df = pd.read_csv(HERE / f"{STEM}.csv")

STYLE = {"KDE": ("black", "-", 3.0), "MTM full": ("darkorange", "--", 1.4),
         "MTM moments": ("teal", "--", 1.4), "Jump KM corr.": ("olive", "--", 1.4),
         "Jump KM flag": ("firebrick", "--", 1.4)}

fig, ax = plt.subplots(figsize=(5.8, 4.2))
for series, (color, ls, lw) in STYLE.items():
    g = df[df.series == series]
    ax.plot(g.x, g.y, color=color, linestyle=ls, lw=lw, label=series)
ax.axvline(-0.5, color="#374151", lw=2.5)
ax.axvspan(-0.5, 0.0, color="#9CA3AF", alpha=0.22, lw=0)
ax.set_xlim(-0.55, 2.5)
ax.set_ylim(bottom=0)
ax.set_xlabel(r"$y_{\mathrm{eff}}$")
ax.set_ylabel(r"$p(y_{\mathrm{eff}})$")
ax.set_title("LRLS")
ax.legend(fontsize=8)

fig.tight_layout()
fig.savefig(HERE / f"{STEM}.png", dpi=200)
print(f"saved {STEM}.png")
