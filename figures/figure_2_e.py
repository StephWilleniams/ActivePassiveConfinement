"""Figure 2(e): mean pre-kick occupancy time vs kick rate k, for three rotational diffusivities, with the theoretical Poisson trigger time.

Renders a simple version of this panel from figure_2_e.csv alone (no processing).
Run: python3 figure_2_e.py
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = "figure_2_e"
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
ax.set_xlabel(r"$k$ (s$^{-1}$)")
ax.set_ylabel(r"$T_{\mathrm{pre\text{-}kick}}$ (s)")
ax.grid(True, which="both", color="0.85", lw=0.5, alpha=0.6)
ax.set_axisbelow(True)
ax.legend(fontsize=8)

fig.tight_layout()
fig.savefig(HERE / f"{STEM}.png", dpi=200)
print(f"saved {STEM}.png")
