"""SM Figure S1(c): near-wall passive diffusivity D_y(y_eff) at four integration time steps (multiplier m); the production step is m = 160.

Renders a simple version of this panel from figure_s1_c.csv alone (no processing).
Run: python3 figure_s1_c.py
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = "figure_s1_c"
df = pd.read_csv(HERE / f"{STEM}.csv")

STYLE = {40: dict(color="#9ecae1", lw=1.5, ls="-"),
         80: dict(color="#4292c6", lw=1.5, ls="-"),
         160: dict(color="#08306b", lw=2.4, ls="-"),
         320: dict(color="#f16913", lw=1.6, ls="--")}

fig, ax = plt.subplots(figsize=(5.2, 4.0))
for mult, g in df.groupby("multiplier"):
    style = STYLE.get(int(mult), dict(color="grey", lw=1.4, ls="-"))
    label = rf"$m={mult}$" + (" (results)" if mult == 160 else "")
    ax.plot(g.x, g.y, label=label, **style)
ax.set_xlabel(r"distance from accessible edge, $y_{\mathrm{eff}}$")
ax.set_ylabel(r"$D_y$")
ax.set_yscale("log")
ax.legend(fontsize=8, frameon=False)

fig.tight_layout()
fig.savefig(HERE / f"{STEM}.png", dpi=200)
print(f"saved {STEM}.png")
