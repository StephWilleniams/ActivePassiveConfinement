"""Figure 3(b): heatmap of passive boundary accumulation E^p over the (k, tau) kick-parameter plane (25 x 25 grid over 4096 Saltelli runs).

Renders a simple version of this panel from figure_3_b.csv alone (no processing).
Run: python3 figure_3_b.py
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = "figure_3_b"
df = pd.read_csv(HERE / f"{STEM}.csv")

piv = df.pivot(index="tau", columns="k", values="E")
k = piv.columns.values
tau = piv.index.values
dk, dt = k[1] - k[0], tau[1] - tau[0]

fig, ax = plt.subplots(figsize=(5.6, 4.2))
im = ax.imshow(piv.values, origin="lower", aspect="auto", cmap="viridis",
               extent=[k[0] - dk / 2, k[-1] + dk / 2,
                       tau[0] - dt / 2, tau[-1] + dt / 2])
fig.colorbar(im, ax=ax, label=r"$E^{\,p} = \phi^{p}/\phi_{\rm unif}$")
ax.set_xlabel(r"$k$ (s$^{-1}$)")
ax.set_ylabel(r"$\tau$")

fig.tight_layout()
fig.savefig(HERE / f"{STEM}.png", dpi=200)
print(f"saved {STEM}.png")
