"""Figure 2(g): polar rose of post-scatter exit angle for five kick strengths (k = 45, D_r = 0.05). theta_out = 0 points along the boundary (south), pi/2 directly away from it (east).

Renders a simple version of this panel from figure_2_g.csv alone (no processing).
Run: python3 figure_2_g.py
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = "figure_2_g"
df = pd.read_csv(HERE / f"{STEM}.csv")

import numpy as np

taus = sorted(df.tau.unique())
colors = plt.cm.viridis(np.linspace(0, 1, len(taus)))
fig = plt.figure(figsize=(5.5, 5.0))
ax = fig.add_subplot(projection="polar")
ax.set_theta_zero_location("S")
ax.set_theta_direction(1)
ax.set_thetamin(0)
ax.set_thetamax(180)
width = np.pi / 180        # 180 bins over [0, pi]
for tau, c in zip(taus, colors):
    g = df[df.tau == tau]
    ax.bar(g.x, g.y, width=width, color=(*c[:3], 0.85), edgecolor=c,
           linewidth=0.6, label=rf"$\tau = {tau:.1f}$")
ax.legend(loc="center left", bbox_to_anchor=(1.05, 0.5), fontsize=8)

fig.tight_layout()
fig.savefig(HERE / f"{STEM}.png", dpi=200)
print(f"saved {STEM}.png")
