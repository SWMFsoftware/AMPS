#!/usr/bin/env bash
# Example shieldSim command lines.  Run from the build directory after compiling.

set -euo pipefail

# Single run with normal-incidence beam.
./shieldSim --source-mode=beam --shield=G4_Al:2 --events=50000

# Isotropic source over the upstream plane.
./shieldSim --source-mode=isotropic --shield=G4_Al:2 --events=50000

# Synthetic tabulated spectrum example.
./shieldSim --source-mode=isotropic --spectrum=examples/sep_spectrum.dat \
            --shield=G4_Al:2 --events=100000

# Dose-vs-thickness sweep.
./shieldSim --sweep --source-mode=isotropic --sweep-material=G4_Al \
            --sweep-tmin=0.5 --sweep-tmax=30 --sweep-n=15 \
            --sweep-log --events=20000
