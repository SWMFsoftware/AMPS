#!/bin/bash
# One LiouvilleCOPS fit over several GPU nodes (SCALING_PLAN.md P3/P0; template, not yet run on Voyager).
#
# One MPI rank per GPU. Rank 0 runs the unchanged fitter and the first block of rows on its own GPU;
# ranks >= 1 serve the other blocks (cops_cuda/transport.py). With --gpus-per-task=1 every rank sees
# exactly one device (device 0 of the rank). The node count and the GPUs per node are NOT in this file
# (user decision 2026-10-22: nothing hardwired); give them at submission, matching the partition:
#
#   sbatch --nodes=N --ntasks-per-node=<GPUs per node> tools/cluster/fit_mpi.slurm \
#       <base run> <record store> <out dir> [fit_cops arguments ...]
#
# (`sinfo -h -o %G` shows the partition's GPUs per node.) The lines marked SITE must be set for the cluster
# (partition, account, modules) once the allocation exists. The campaign runner
# tools/cluster/cops_pipeline.py writes complete scripts of its own (submit).
#SBATCH --job-name=cops_fit
#SBATCH --gpus-per-task=1
#SBATCH --gpu-bind=closest             # each rank on the socket of its GPU (SCALING_PLAN 3.4)
##SBATCH --cpus-per-task=SITE          # cores per node / GPUs per node
#SBATCH --mem=0                        # all node memory: pinned host buffers for streamed records
#SBATCH --time=02:00:00
##SBATCH --partition=SITE
##SBATCH --account=SITE
#SBATCH --output=%x-%j.out

set -euo pipefail
if [ -z "${SLURM_NTASKS_PER_NODE:-}" ]; then
    echo "fit_mpi.slurm: submit with --nodes=N --ntasks-per-node=<GPUs per node> (one rank per GPU)" >&2
    exit 2
fi
BASE=$1; STORE=$2; OUT=$3; shift 3
APP=${COPS_CUDA_APP:-$SLURM_SUBMIT_DIR}          # the LiouvilleCOPS_cuda checkout
# SITE: module load cuda/12.x openmpi/4.x (the MPI mpi4py was built against)
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1} PYTHONUNBUFFERED=1
FITARGS=${FITARGS:-"--profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6"}

echo "fit_mpi.slurm: $SLURM_JOB_NUM_NODES node(s) x $SLURM_NTASKS_PER_NODE rank(s) = $SLURM_NTASKS ranks"
nvidia-smi topo -m > "$OUT.topology.txt" 2>&1 || true
cd "$APP"
srun --mpi=pmix --kill-on-bad-exit=1 .venv/bin/python -W ignore tools/gpu_fit.py \
    --transport mpi --precision mixed --record-store "$STORE" --headroom-fraction 0.12 \
    -- "$BASE" $FITARGS --out "$OUT" "$@"
