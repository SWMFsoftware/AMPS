#!/bin/bash
# Example SLURM job for one production case of the I3 campaign qfit_20261123 (adapt partition, time limit
# and launcher to the machine). The campaign used one exclusive node with 128 ranks:
#
#   for c in 1 2 3 4 5 6; do
#     sbatch --partition=<partition> --time=7-00:00:00 --ntasks=128 --job-name=i3qfit-case$c \
#            --export=ALL,CASE=$c,AMPS=/path/to/build/amps,OUT=/path/to/campaign_output \
#            srcHartley2/i3/campaigns/qfit_20261123/scripts/production.sbatch
#   done
#
# Optional: CHEMISTRY=state_resolved (default sa1982_vh), MPIRUN (default mpiexec).
# Wall time in the campaign: 1.7-3.1 days (convergence after 70000-100000 iterations).
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mem=0
set -euo pipefail
: "${CASE:?CASE must be set}" "${AMPS:?AMPS must be set}" "${OUT:?OUT must be set}"
SCRIPTS=${SCRIPTS:-$(cd "$(dirname "$(scontrol show job "$SLURM_JOB_ID" | awk -F= '/Command=/{print $2}')")" && pwd)}
exec "$SCRIPTS/run_case.sh" --amps "$AMPS" --case "$CASE" --out "$OUT" --np "$SLURM_NTASKS" \
     --chemistry "${CHEMISTRY:-sa1982_vh}"
