"""SM Figure S1(a): histogram of the per-frame closest passive-passive approach r_min/sigma (HRHS trajectory). Dashed line: WCA cutoff 2^(1/6); dotted line: sigma.

Renders a simple version of this panel from figure_s1_a.csv alone (no processing).
Run: python3 figure_s1_a.py
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = "figure_s1_a"
df = pd.read_csv(HERE / f"{STEM}.csv")

fig, ax = plt.subplots(figsize=(5.2, 4.0))
ax.bar(df.x, df.y, width=df.x.iloc[1] - df.x.iloc[0], color="tab:blue", alpha=0.85)
ax.axvline(2 ** (1 / 6), color="black", ls="--", lw=1.5)
ax.axvline(1.0, color="grey", ls=":", lw=1.5)
ax.set_xlabel(r"closest passive pair approach, $r_{\min}/\sigma$")
ax.set_ylabel("frames")
ax.set_yscale("log")

fig.tight_layout()
fig.savefig(HERE / f"{STEM}.png", dpi=200)
print(f"saved {STEM}.png")
