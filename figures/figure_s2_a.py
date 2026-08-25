"""SM Figure S2(a): the WCA interaction potential; the dotted line marks the cutoff r = 2^(1/6) sigma.

Renders a simple version of this panel from figure_s2_a.csv alone (no processing).
Run: python3 figure_s2_a.py
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = "figure_s2_a"
df = pd.read_csv(HERE / f"{STEM}.csv")

fig, ax = plt.subplots(figsize=(5.2, 4.0))
ax.plot(df.x, df.y, color="tab:blue", lw=2.0)
ax.axvline(2 ** (1 / 6), color="grey", ls=":", lw=1.5)
ax.axhline(0.0, color="black", lw=0.8)
ax.set_xlabel(r"$r/\sigma$")
ax.set_ylabel(r"$U_{\mathrm{WCA}}/\epsilon$")
ax.set_ylim(-0.2, 4.0)
ax.set_xlim(df.x.min(), df.x.max())

fig.tight_layout()
fig.savefig(HERE / f"{STEM}.png", dpi=200)
print(f"saved {STEM}.png")
