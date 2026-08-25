#!/bin/bash

# Build and run the active-passive suspension simulator.
# Usage: ./runscript.sh [multiplier]
#   multiplier scales base_steps up / base_dt down together (fixed total
#   simulated time, finer internal resolution). Production runs use 160.
set -e

MULTIPLIER="${1:-160}"

COMPILER="g++"
STANDARD="-std=c++17"
# -ffast-math is intentionally omitted: it perturbs the last decimal digit of
# the output for a negligible (~3%) speed gain on this branchy/scalar hot path.
OPTIMISATION="-O3 -march=native -ftree-vectorize -funroll-loops"

# Requires yaml-cpp (e.g. `brew install yaml-cpp` on macOS,
# `apt install libyaml-cpp-dev` on Debian/Ubuntu).
if command -v brew >/dev/null 2>&1; then
  INCLUDES="-I$(brew --prefix yaml-cpp)/include"
  LIBS="-L$(brew --prefix yaml-cpp)/lib -lyaml-cpp"
else
  INCLUDES=""
  LIBS="-lyaml-cpp"
fi

mkdir -p outputs

echo "--> Compiling main.cpp..."
$COMPILER $STANDARD $OPTIMISATION $INCLUDES main.cpp -o simulate_active_matter $LIBS

echo "--> Running (multiplier=$MULTIPLIER)..."
./simulate_active_matter "$MULTIPLIER" 1 config.yml "outputs/output_${MULTIPLIER}.txt"

echo "--> Done. Trajectory written to outputs/output_${MULTIPLIER}.txt"
