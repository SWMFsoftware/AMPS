#!/bin/bash
#======================================================================================
# amps_gpu_select.sh -- start a GPU amps on one of several GPUs, chosen by the MPI
# local rank.
#
# The engine's device pool uses the default device of its process, so with several MPI
# ranks per node every rank would use GPU 0.  Use this script in place of the amps
# executable, e.g.
#
#   GEOSPACE_AMPS=$GEOSPACE_RUNS/builds/gpu_double/amps GEOSPACE_GPUS=0,1 \
#     srcGeospace/test/tools/run_srcearth_tests.sh --run-root <dir> \
#       --amps srcGeospace/scripts/amps_gpu_select.sh --np 2 --nt 16
#
#   GEOSPACE_AMPS   the GPU amps executable (required)
#   GEOSPACE_GPUS   comma-separated device ids (default: 0)
#
# Local rank r (OMPI_COMM_WORLD_LOCAL_RANK, MPI_LOCALRANKID or SLURM_LOCALID; 0 for a
# run without mpirun) runs with CUDA_VISIBLE_DEVICES set to entry (r mod n) of the
# list, so each process sees exactly one GPU.  mpirun passes the environment of the
# test runner to the ranks.
#======================================================================================
set -euo pipefail
: "${GEOSPACE_AMPS:?set GEOSPACE_AMPS to the GPU amps executable}"
IFS=, read -r -a gpus <<< "${GEOSPACE_GPUS:-0}"
local_rank=${OMPI_COMM_WORLD_LOCAL_RANK:-${MPI_LOCALRANKID:-${SLURM_LOCALID:-0}}}
export CUDA_VISIBLE_DEVICES=${gpus[$(( local_rank % ${#gpus[@]} ))]}
exec "$GEOSPACE_AMPS" "$@"
