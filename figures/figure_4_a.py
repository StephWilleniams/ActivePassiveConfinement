"""Figure 4(a): Sobol first-order (S1) and total-order (ST) sensitivity indices of passive and active accumulation, with 95% bootstrap CIs and the Monte Carlo noise-floor band hatched onto each total-order bar.

Renders a simple version of this panel from figure_4_a.csv alone (no processing).
Run: python3 figure_4_a.py
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = "figure_4_a"
df = pd.read_csv(HERE / f"{STEM}.csv")

import numpy as np

BAR_COLORS = {("passive", "S1"): "#4C72B0", ("passive", "ST"): "#DD8452",
              ("swimmer", "S1"): "#55A868", ("swimmer", "ST"): "#8172B2"}
BAR_LABELS = {("passive", "S1"): "S1P", ("passive", "ST"): "STP",
              ("swimmer", "S1"): "S1A", ("swimmer", "ST"): "STA"}

params = list(dict.fromkeys(df["param"]))
x = np.arange(len(params))
width = 0.18
half_gap = (width + 0.08) / 2
offsets = {("passive", "S1"): -half_gap - width, ("passive", "ST"): -half_gap,
           ("swimmer", "S1"): half_gap, ("swimmer", "ST"): half_gap + width}

fig, ax = plt.subplots(figsize=(5.6, 4.2))
for species in ("passive", "swimmer"):
    sub = df[df["species"] == species].set_index("param").loc[params]
    for order in ("S1", "ST"):
        ax.bar(x + offsets[(species, order)], sub[order].values, width,
               yerr=sub[f"{order}_conf"].values, capsize=4,
               color=BAR_COLORS[(species, order)], label=BAR_LABELS[(species, order)])
    # noise floor, hatched from S1 up to S1 + F on the total-order bar
    band = sub.iloc[0]
    if np.isfinite(band["F_corrected"]):
        x_st = x + offsets[(species, "ST")]
        ax.bar(x_st, band["F_corrected"], width, bottom=sub["S1"].values,
               fill=False, hatch="////", edgecolor="black", linewidth=1.0, zorder=4)
        for xi, s1i in zip(x_st, sub["S1"].values):
            for edge in ("F_ci_lo", "F_ci_hi"):
                ax.plot([xi - width / 2, xi + width / 2], [s1i + band[edge]] * 2,
                        color="black", lw=0.8, linestyle=":", zorder=5)
ax.axhline(0, color="gray", lw=0.8, zorder=0)
ax.set_xticks(x)
ax.set_xticklabels(["$k$", r"$\tau$"])
ax.set_ylabel(r"$s_{\rm sobol}$")
ax.set_ylim(0, 0.82)
ax.legend(ncol=2, fontsize=8)

fig.tight_layout()
fig.savefig(HERE / f"{STEM}.png", dpi=200)
print(f"saved {STEM}.png")
