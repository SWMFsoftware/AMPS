#!/bin/sh
set -eu

# This runner never fabricates native or external evidence.  Projects with a
# linked AMPS/SWMF campaign pass its exact command through SRCSEP_NATIVE_GATE;
# dependency-light distributions report BLOCKED and still return success so the
# routine source suite can distinguish unavailable infrastructure from failure.
if [ -z "${SRCSEP_NATIVE_GATE:-}" ]; then
  echo "SRCSEP_SUITE_RESULT=SKIP"
  echo "SKIP WP59-WP60: set SRCSEP_NATIVE_GATE to the linked AMPS matrix/restart command"
  echo "SKIP WP63: real SWMF replay and held-out observational manifests are not bundled"
  echo "SKIP WP64: multi-node scaling requires a frozen hardware environment"
  exit 0
fi

# The command is explicitly operator supplied.  eval is intentional here so a
# complete MPI launcher plus arguments can be provided; production CI should
# set this variable from a reviewed, immutable workflow rather than user data.
eval "$SRCSEP_NATIVE_GATE"
