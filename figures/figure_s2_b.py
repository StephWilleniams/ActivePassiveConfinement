"""SM Figure S2(b): the corresponding WCA repulsive force magnitude; the dotted line marks the cutoff r = 2^(1/6) sigma.

Renders a simple version of this panel from figure_s2_b.csv alone (no processing).
Run: python3 figure_s2_b.py
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = "figure_s2_b"
df = pd.read_csv(HERE / f"{STEM}.csv")

fig, ax = plt.subplots(figsize=(5.2, 4.0))
ax.plot(df.x, df.y, color="tab:red", lw=2.0)
ax.axvline(2 ** (1 / 6), color="grey", ls=":", lw=1.5)
ax.set_xlabel(r"$r/\sigma$")
ax.set_ylabel(r"$F_{\mathrm{WCA}}\,\sigma/\epsilon$")
ax.set_yscale("log")
ax.set_ylim(1e-1, 1e3)
ax.set_xlim(df.x.min(), df.x.max())

fig.tight_layout()
fig.savefig(HERE / f"{STEM}.png", dpi=200)
print(f"saved {STEM}.png")
