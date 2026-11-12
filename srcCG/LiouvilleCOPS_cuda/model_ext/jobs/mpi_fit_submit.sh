#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/mpi_fit_submit.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/mpi_fit_submit.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: f2fc180dd46e2522d7455779dc2082feea9e4997e104f0e21ca71af94332fe55
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: submit helper. Writes and submits the sbatch line of jobs/mpi_fit.sbatch for NODES 1 | 2 | auto (sbatch
#   --nodes=1-2), exports the VAR=value settings and logs the line in model_ext/logs/submissions.txt (2026-11-10).
# RUN: mpi_fit_submit.sh NODES SPEC [VAR=value ...] [-- <sbatch options>]   (MPI tests 123629-123649)
# ==== end of header ====
# Submit jobs/mpi_fit.sbatch for one fit (2026-11-10) and log the sbatch line in model_ext/logs/submissions.txt.
#   jobs/mpi_fit_submit.sh NODES SPEC [VAR=value ...] [-- <more sbatch options>]
#     NODES  1 (8 ranks) | 2 (16 ranks) | auto (sbatch --nodes=1-2: 2 nodes if both are free at the start, else 1)
#     SPEC   the fit spec (see jobs/mpi_fit.sbatch)
#     VAR=value  exported to the job (MODEL_EXT=..., ACTIVITY=..., TRANSPORT=..., spec variables); no commas
#     after --: sbatch options, e.g. --time=04:00:00 --nodelist=gpu001,gpu002 --dependency=afterany:<job> --job-name=x
#   DESC="text" in the environment: the description written into submissions.txt.
# Prints the job id. Resume after a wall-time stop (exit 75): run the same line again (any NODES).
set -euo pipefail
ROOT=/data001/highresolu_vt/cops_output
[ $# -ge 2 ] || { sed -n '2,10p' "$0"; exit 2; }
NODES=$1; SPEC=$(readlink -f "$2"); shift 2
[ -f "$SPEC" ] || { echo "no spec $SPEC" >&2; exit 2; }
case $NODES in 1) NOPT=--nodes=1;; 2) NOPT=--nodes=2;; auto) NOPT=--nodes=1-2;; *) echo "NODES = 1, 2 or auto" >&2; exit 2;; esac
EXP="ALL,NODES=$NODES,SPEC=$SPEC"; SB=()
while [ $# -gt 0 ]; do
    case $1 in
        --) shift; SB+=("$@"); break;;
        *=*) case $1 in *,*) echo "no commas in exported values: $1" >&2; exit 2;; esac; EXP+=",$1";;
        *) echo "not VAR=value: $1" >&2; exit 2;;
    esac; shift
done
CMD=(sbatch --parsable "$NOPT" "${SB[@]}" "--export=$EXP" "$ROOT/jobs/mpi_fit.sbatch")
J=$("${CMD[@]}")
echo "$(date -u +%FT%TZ) $J $ROOT/jobs/mpi_fit.sbatch NODES=$NODES ${DESC:-}: ${CMD[*]:2}" >> "$ROOT/model_ext/logs/submissions.txt"
echo "$J"
