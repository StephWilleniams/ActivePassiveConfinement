"""SM Figure S1(b): near-wall passive density vs |y| (HRHS trajectory). Dotted line: WCA standoff (accessible edge, 2.4388); dashed line: hard wall y = 3.

Renders a simple version of this panel from figure_s1_b.csv alone (no processing).
Run: python3 figure_s1_b.py
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = "figure_s1_b"
df = pd.read_csv(HERE / f"{STEM}.csv")

fig, ax = plt.subplots(figsize=(5.2, 4.0))
ax.plot(df.x, df.y, color="tab:blue", lw=2.0)
ax.axvline(2.4388, color="grey", ls=":", lw=1.5)
ax.axvline(3.0, color="black", ls="--", lw=1.5)
ax.set_xlabel(r"$|y|$")
ax.set_ylabel("passive density")
ax.set_xlim(2.2, 3.02)
ax.set_yscale("log")

fig.tight_layout()
fig.savefig(HERE / f"{STEM}.png", dpi=200)
print(f"saved {STEM}.png")
