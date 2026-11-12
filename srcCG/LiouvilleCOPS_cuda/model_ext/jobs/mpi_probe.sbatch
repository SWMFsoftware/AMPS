#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/mpi_probe.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/mpi_probe.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: e30b83c63a5a7d77f23dae6866176511b57f30c1b669dd178754383e2537b803
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (test). MPI probe across the GPU nodes (2026-11-10; no fit): srun --mpi=pmix with mpi4py on every node, rank
#   -> GPU map, ping-pong, Bcast and the transport's request/reply round (mpi_multinode/mpi_probe.py); per node
#   nvidia-smi, IB state, UCX devices, free memory.
# RUN: sbatch --nodes=2 --nodelist=gpu001,gpu002 mpi_probe.sbatch   (job 123623 -> mpi_multinode/probe-123623/)
# ==== end of header ====
# MPI probe across the GPU nodes (test, 2026-11-10; no fit, no physics). Checks that srun --mpi=pmix starts mpi4py
# ranks (Open MPI 4.1.7rc1 of /usr/mpi/gcc, the MPI mpi4py 4.0.3 is built with) on every node, which GPU each rank
# takes (cops_cuda/transport.py: local rank modulo the visible devices), and measures point-to-point, Bcast and the
# transport's request/reply round with host buffers. Per node: nvidia-smi -L, IB state, UCX devices, free memory.
#   sbatch --nodes=2 --nodelist=gpu001,gpu002 jobs/mpi_probe.sbatch      (one node: --nodes=1)
# Products: $COPS_OUTPUT_ROOT/mpi_multinode/probe-<job>/ (rank JSON, bench.json, per-node files, logs).
#SBATCH --job-name=cops_mpi_probe
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=8
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=00:15:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --exclude=gpu003
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_mpi_probe-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
O=$COPS_OUTPUT_ROOT/mpi_multinode/probe-$SLURM_JOB_ID; mkdir -p "$O"
cp -p "$0" "$O/mpi_probe.sbatch"; cp -p "$COPS_OUTPUT_ROOT/mpi_multinode/mpi_probe.py" "$O/"
PY="$APP/.venv/bin/python -W ignore"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID: $SLURM_JOB_NUM_NODES node(s) $SLURM_JOB_NODELIST, $SLURM_NTASKS tasks"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$O/git.txt" 2>&1
{ ompi_info --version; ompi_info | grep -E 'MCA (pml|btl|mtl|coll: hcoll|pmix|ess: slurm|plm: slurm)'; srun --mpi=list; sinfo --version; } > "$O/mpi_stack.txt" 2>&1
# per node
srun --nodes="$SLURM_JOB_NUM_NODES" --ntasks-per-node=1 --overlap bash -c '
h=$(hostname -s); f='"$O"'/node_$h.txt
{ echo "== $h"; nvidia-smi -L; nvidia-smi topo -m; echo "== ibstat"; ibstat 2>&1 | grep -E "CA |State|Rate|Link layer|Port [0-9]"
  echo "== ibv_devinfo -l"; ibv_devinfo -l 2>&1; echo "== ucx_info -v"; ucx_info -v 2>&1
  echo "== ucx devices"; ucx_info -d 2>&1 | grep -E "^#\s+(Transport|Device)"; echo "== free -g"; free -g
  echo "== CUDA_VISIBLE_DEVICES in a 1-task step: ${CUDA_VISIBLE_DEVICES:-unset}"; } > "$f" 2>&1
nvidia-smi -q > '"$O"'/nvidia_smi_q_$h.txt 2>&1'
# ranks: default MCA selection, with the PML selection log (one rank per node prints it)
export OMPI_MCA_pml_base_verbose=10
srun --mpi=pmix --kill-on-bad-exit=1 --nodes="$SLURM_JOB_NUM_NODES" --ntasks-per-node=8 \
    $PY "$O/mpi_probe.py" "$O/default" --devices 0,1,2,3,4,5,6,7 > "$O/default.log" 2>&1; rc=$?
unset OMPI_MCA_pml_base_verbose
echo "srun default: exit $rc"; grep '^probe:' "$O/default.log"; grep -i -E 'select|pml' "$O/default.log" | sort | uniq -c | head -20
# the same with UCX forced (fails loudly if UCX cannot be used)
OMPI_MCA_pml=ucx srun --mpi=pmix --kill-on-bad-exit=1 --nodes="$SLURM_JOB_NUM_NODES" --ntasks-per-node=8 \
    $PY "$O/mpi_probe.py" "$O/ucx" --devices 0,1,2,3,4,5,6,7 > "$O/ucx.log" 2>&1; rc2=$?
echo "srun pml=ucx: exit $rc2"; grep -E '^probe: (pingpong|bcast|request)' "$O/ucx.log"
# mpirun (Open MPI's own launcher inside the allocation) as the alternative
mpirun -np "$SLURM_NTASKS" --map-by ppr:8:node --bind-to none -x PATH -x LD_LIBRARY_PATH -x CUPY_CACHE_DIR \
    -x OMP_NUM_THREADS -x OPENBLAS_NUM_THREADS -x MKL_NUM_THREADS -x PYTHONUNBUFFERED \
    $PY "$O/mpi_probe.py" "$O/mpirun" --devices 0,1,2,3,4,5,6,7 > "$O/mpirun.log" 2>&1; rc3=$?
echo "mpirun: exit $rc3"; grep -E '^probe: (rank|pingpong|bcast|request)' "$O/mpirun.log" | head -40
grep -i -E 'warn|error' "$O"/*.log | sort | uniq -c | head -20
echo "$(date -u +%FT%TZ) done: default $rc, ucx $rc2, mpirun $rc3"
exit $rc
