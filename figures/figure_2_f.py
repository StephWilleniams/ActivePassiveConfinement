"""Figure 2(f): mean kick duration vs kick strength tau, for three rotational diffusivities.

Renders a simple version of this panel from figure_2_f.csv alone (no processing).
Run: python3 figure_2_f.py
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = "figure_2_f"
df = pd.read_csv(HERE / f"{STEM}.csv")

COLORS = {"Dr=0.00": "#2a78d6", "Dr=0.05": "#008300", "Dr=0.10": "#eb6834"}

fig, ax = plt.subplots(figsize=(5.0, 4.0))
for series, g in df.groupby("series", sort=False):
    if series == "Poisson":
        ax.plot(g.x, g.y, "k--", lw=1.8, label="Poisson")
    else:
        label = rf"$D_r = {series.split('=')[1]}$"
        ax.errorbar(g.x, g.y, yerr=g.yerr, fmt="o-", ms=3.5, lw=2.0, capsize=2,
                    color=COLORS[series], label=label)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"$\tau$")
ax.set_ylabel(r"$T_{\mathrm{kick}}$ (s)")
ax.grid(True, which="both", color="0.85", lw=0.5, alpha=0.6)
ax.set_axisbelow(True)
ax.legend(fontsize=8)

fig.tight_layout()
fig.savefig(HERE / f"{STEM}.png", dpi=200)
print(f"saved {STEM}.png")
