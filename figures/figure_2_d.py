"""Figure 2(d): scattering time and exit-angle histograms, HRHS (k = 45, tau = 13.5).

Renders a simple version of this panel from figure_2_d.csv alone (no processing).
Run: python3 figure_2_d.py
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = "figure_2_d"
df = pd.read_csv(HERE / f"{STEM}.csv")

COLOR_TIME = "#2a78d6"
COLOR_ANGLE = "#eb6834"

t = df[df.series == "time_hist"]
a = df[df.series == "angle_hist"]
fig, (ax_t, ax_a) = plt.subplots(1, 2, figsize=(8.5, 3.6))
ax_t.bar(t.x, t.y, width=t.x.iloc[1] - t.x.iloc[0], color=COLOR_TIME, alpha=0.55)
ax_t.set_xlabel(r"$T_{\mathrm{total}}$ (s)")
ax_t.set_ylabel(r"$p(T_{\mathrm{total}})$")
ax_a.bar(a.x, a.y, width=a.x.iloc[1] - a.x.iloc[0], color=COLOR_ANGLE, alpha=0.55)
ax_a.set_xlabel(r"$\theta_{\mathrm{out}}$ (rad)")
ax_a.set_ylabel(r"$p(\theta_{\mathrm{out}})$")
ax_a.set_xlim(0, 3.14159)

fig.tight_layout()
fig.savefig(HERE / f"{STEM}.png", dpi=200)
print(f"saved {STEM}.png")
