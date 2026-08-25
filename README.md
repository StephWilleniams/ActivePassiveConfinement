# Near-boundary active kicks control passive distributions in confined active-passive suspensions

Minimal release package for the paper: the numerical simulation, its configuration,
and the plotted data behind every figure panel.

## Contents

- `main.cpp` — the N-swimmer / N-passive Brownian-dynamics simulation
  (three-segment self-propelled swimmers, single-circle passives, WCA interactions,
  stochastic boundary kicks near the top/bottom walls).
- `config.yml` — the simulation configuration (HRHS kick parameters shown;
  `kick_rate` / `kick_torque` are the swept pair).
- `runscript.sh` — builds and runs the simulator: `./runscript.sh [multiplier]`
  (default multiplier = 160, the time step `dt = 0.1/160 = 6.25e-4 s`).
  Requires `g++` (C++17) and yaml-cpp (`brew install yaml-cpp` /
  `apt install libyaml-cpp-dev`).
- `figures/` — one `figure_<n>_<panel>.csv` + `figure_<n>_<panel>.py` pair per
  panel of the manuscript (main Figs. 1–4) and the SM (Figs. S1–S4).
  Each CSV holds the final plotted values (x, y per series — no
  processing needed); each script renders a simple version of its panel from its
  sibling CSV alone and writes `figure_<n>_<panel>.png`.

  ```bash
  cd figures
  python3 figure_3_b.py        # or any other panel
  ```

  Requires Python 3 with `pandas`, `numpy`, `matplotlib`.

## Panel index

| Pair | Panel |
|---|---|
| `figure_1_d`, `figure_1_e`  | Stationary distributions of both species: no-kick baseline and strongest kick regime (Fig. 1 a–c are schematics with no data) |
| `figure_2_a`–`figure_2_d`   | Scattering-time and exit-angle histograms at the four kick-parameter extremes (baseline, LRHS, HRLS, HRHS) |
| `figure_2_e`, `figure_2_f`  | Mean pre-kick occupancy time vs $k$; mean kick duration vs $\tau$ |
| `figure_2_g`                | Polar plot of exit angle across five kick strengths |
| `figure_3_a`, `figure_3_b`  | Accumulation heatmaps $E^a$, $E^p$ over the $(k, \tau)$ plane |
| `figure_3_c`, `figure_3_d`  | Passive-accumulation marginals vs $k$ and $\tau$, coloured by tercile of the other parameter |
| `figure_4_a`                | Sobol sensitivity indices with bootstrap CIs and noise-floor bands |
| `figure_4_b`, `figure_4_c`  | Wall-normal drift $V_y(y_{\rm eff})$ and diffusivity $D_y(y_{\rm eff})$ at the four corners (raw + LOESS trend) |
| `figure_4_d`, `figure_4_e`  | LRLS / HRHS stationary distributions: measured KDE vs MTM and jump-KM reconstructions |
| `figure_s1_a`–`figure_s1_d` | Numerical integrity: closest pair approach, near-wall density, $D_y$ time-step convergence, limiter engagement |
| `figure_s2_a`, `figure_s2_b`| WCA potential and force (closed form) |
| `figure_s3`                 | Kick-strength bounds from deterministic escape kinematics (closed form) |
| `figure_s4_a`, `figure_s4_b`| Fig. 4(d, e) as CDFs |

## Simulator output format

CSV, no header, one row per particle per recorded step:
`time, type(S=swimmer|P=passive), x, y, theta, interacting_flag`.
