#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/probe_gpu.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/probe_gpu.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: d098b2e2c87b2eb4f1b5a450b1f4bdc7b08835564bd365212dc2e14662bc34ba
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: site probe (handoff §3): GPU model, memory, MIG, driver, host CPU/RAM, file systems. Read only. Job 120809.
# RUN: sbatch probe_gpu.sbatch
# ==== end of header ====
#SBATCH --job-name=cops_probe
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=00:20:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/%x-%j.out
# Site probe (handoff §3): GPU model/memory/MIG/driver, host CPU/RAM, file-system visibility.
# Read-only: writes only into $D. No builds, no GPU computation.
set -u
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=/data001/highresolu_vt/cops_output/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
D=/data001/highresolu_vt/cops_output/cluster_site
cd "$D" || exit 1
cp -p "$0" "$D/probe_gpu.sbatch.used" 2>/dev/null
echo "$(date -u +%FT%TZ) start on $(hostname)"
nvidia-smi > gpu_nvidia_smi.txt 2>&1
nvidia-smi -q > gpu_query.txt 2>&1
nvidia-smi topo -m > gpu_topology.txt 2>&1
nvidia-smi --query-gpu=index,name,uuid,compute_cap,memory.total,memory.free,memory.used,mig.mode.current,ecc.mode.current,pcie.link.gen.max,pcie.link.width.max,driver_version --format=csv > gpu_memory.csv 2>&1
{ hostname; lscpu; echo; numactl -H; echo; free -g; free -b; echo; ulimit -a; echo; nproc; echo; cat /proc/self/cgroup;
  echo; grep -E 'MemTotal|MemAvailable|HugePages_Total' /proc/meminfo; echo "SLURM_MEM_PER_NODE=${SLURM_MEM_PER_NODE:-unset}";
  echo "SLURM_CPUS_ON_NODE=${SLURM_CPUS_ON_NODE:-unset} CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}";
  echo; df -hT ${TMPDIR:-/tmp}; } > compute_node.txt 2>&1
scontrol show job "$SLURM_JOB_ID" > job_probe.txt 2>&1
{ for p in /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python \
           /data001/highresolu_vt/cops_cluster_transfer_20261026/SHA256SUMS \
           /home/vtenishev/SPICE/cspice/lib/cspice.a /home/vtenishev/SPICE/Kernels/ROSETTA/esa_2022 \
           /scratch/rosetta; do
    if [ -e "$p" ]; then echo "visible  $p"; else echo "MISSING  $p"; fi; done
  echo; df -hT /data001/highresolu_vt /home/vtenishev /scratch 2>&1; mount | grep -E ' /data001/highresolu_vt | /home | /scratch '
  echo; ls -l /usr/local | grep -i cuda; readlink -f /usr/local/cuda; ls /opt/nvidia/hpc_sdk/Linux_x86_64/25.9/cuda 2>&1
  echo; which nvcc gcc mpicc 2>&1; gcc --version | head -1; nvcc --version | tail -1
  echo; module avail 2>&1 | grep -i -E 'cuda|gcc|mpi|nvhpc'; } > compute_paths.txt 2>&1
# CuPy sees the GPUs? (device query only; no kernels compiled, so no ~/.cupy cache written)
/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python -c "
import cupy as cp
print('cupy', cp.__version__, 'runtime', cp.cuda.runtime.runtimeGetVersion(), 'driver', cp.cuda.runtime.driverGetVersion())
for i in range(cp.cuda.runtime.getDeviceCount()):
    with cp.cuda.Device(i):
        p = cp.cuda.runtime.getDeviceProperties(i); f, t = cp.cuda.runtime.memGetInfo()
        print(i, p['name'].decode(), 'cc', p['major'], p['minor'], 'free', f, 'total', t)
" > cupy_devices.txt 2>&1
curl -sI --max-time 10 https://pypi.org/simple/ | head -1 > internet_compute.txt 2>&1
echo "$(date -u +%FT%TZ) done"
