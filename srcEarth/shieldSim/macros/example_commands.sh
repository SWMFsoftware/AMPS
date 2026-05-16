#!/usr/bin/env bash
# Example shieldSim command lines.  Run from the build directory after compiling.

set -euo pipefail

# Show material catalogs.
./shieldSim --list-materials
./shieldSim --list-target-materials

# Single run with normal-incidence beam.
./shieldSim --source-mode=beam --shield=Al:2 --events=50000

# Isotropic source over the upstream plane with high-precision neutron transport.
./shieldSim --physics-list=FTFP_BERT_HP --source-mode=isotropic \
            --shield=Al:2 --events=50000

# Synthetic tabulated spectrum example.
./shieldSim --source-mode=isotropic --spectrum=examples/sep_spectrum.dat \
            --shield=Al:2 --events=100000

# Dose-vs-thickness sweep over HDPE.
./shieldSim --physics-list=Shielding --sweep --source-mode=isotropic --sweep-material=HDPE \
            --sweep-tmin=0.5 --sweep-tmax=30 --sweep-n=15 \
            --sweep-log --events=20000

# Tissue and silicon scoring behind aluminum.
./shieldSim --source-mode=isotropic --shield=Al:2 --target=BFO:50,Si:1 --events=50000

# Electronics detector-material scoring stack.
./shieldSim --shield=Al:2 --target=Si:1,SiO2:0.01,GaAs:1,InGaAs:1,Ge:1 --events=50000
