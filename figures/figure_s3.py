"""SM Figure S3: deterministic escape kinematics -- exit angle swept by an initially wall-parallel swimmer escaping the boundary band, vs kick strength tau. Red markers denote the chosen bounds tau in [0.05, 13.5].

Renders a simple version of this panel from figure_s3.csv alone (no processing).
Run: python3 figure_s3.py
"""
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STEM = "figure_s3"
df = pd.read_csv(HERE / f"{STEM}.csv")

curve = df[df.series == "curve"]
bound = df[df.series == "bound"]

fig, ax = plt.subplots(figsize=(5.6, 4.0))
ax.plot(curve.x, curve.y, color="tab:blue", lw=2.0)
ax.axhline(90.0, color="black", lw=1.2, ls="--")
for _, row in bound.iterrows():
    ax.plot([row.x], [row.y], "o", color="tab:red", ms=8, zorder=5)
    ax.axvline(row.x, color="grey", ls=":", lw=1.2)
ax.set_xscale("log")
ax.set_xlabel(r"kick reorientation rate $\tau$  (rad s$^{-1}$)")
ax.set_ylabel(r"swept angle $\theta_{\mathrm{exit}}$  (deg)")
ax.set_ylim(0, 100)

fig.tight_layout()
fig.savefig(HERE / f"{STEM}.png", dpi=200)
print(f"saved {STEM}.png")
