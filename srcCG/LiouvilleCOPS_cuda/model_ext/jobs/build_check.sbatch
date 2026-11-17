#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/build_check.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/build_check.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 8da5191447dc5833af05297268c42f7345f10c8b59a196e2fcee3e00abb99a70
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: site setup (handoff §5): builds liouville_cops and native/liouville_cops_cuda on a GPU node (-march=native),
#   then cuda_doctor, pytest and the native self test. Jobs 120810, 120814.
# RUN: sbatch build_check.sbatch
# ==== end of header ====
#SBATCH --job-name=cops_build
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/%x-%j.out
# Handoff §5: build the two executables on a GPU node (-march=native), then cuda_doctor, pytest, native selftest.
# Each step runs even if an earlier check fails; exit codes are summarized at the end.
set -u
export AMPS=/data001/highresolu_vt/AMPS
export APP=$AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
CSPICE=/home/vtenishev/SPICE/cspice
P=$COPS_OUTPUT_ROOT/build/job-$SLURM_JOB_ID; mkdir -p "$P" || exit 1
cp -p "$0" "$P/build_check.sbatch.used"
stamp() { echo "$(date -u +%FT%TZ) $*"; }
stamp "start on $(hostname), products in $P"

# provenance (CLAUDE.md §3.8)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; g++ --version | head -1; grep -m1 CSPICE $CSPICE/doc/version.txt; mpicc --showme:version 2>&1; } > "$P/toolchain.txt" 2>&1
"$APP/.venv/bin/python" -m pip list > "$P/pip_list.txt" 2>&1 || ~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1

cd "$APP" || exit 1
stamp "make spice_angular_velocity"
make -C ../LiouvilleCOPS CSPICE=$CSPICE spice_angular_velocity > "$P/make_spice_angular_velocity.log" 2>&1; rc_sav=$?
stamp "make liouville_cops_cuda (sm_80)"
make -C native CSPICE=$CSPICE CUDA=$CUDA_HOME ARCH=sm_80 > "$P/make_liouville_cops_cuda.log" 2>&1; rc_nat=$?
ls -l ../LiouvilleCOPS/spice_angular_velocity native/liouville_cops_cuda > "$P/executables.txt" 2>&1
sha256sum ../LiouvilleCOPS/spice_angular_velocity native/liouville_cops_cuda >> "$P/executables.txt" 2>&1

stamp "cuda_doctor"
.venv/bin/python tools/cuda_doctor.py --json "$P/cuda_doctor.json" > "$P/cuda_doctor.log" 2>&1; rc_doc=$?
cp -p "$P/cuda_doctor.json" "$COPS_OUTPUT_ROOT/cluster_site/cuda_doctor_job$SLURM_JOB_ID.json" 2>/dev/null

stamp "pytest"
.venv/bin/python -m pytest -q -rs tests > "$P/pytest.log" 2>&1; rc_pyt=$?

stamp "native selftest"
( cd ../LiouvilleCOPS && env COPS_CUDA_DEVICE=0 ../LiouvilleCOPS_cuda/native/liouville_cops_cuda selftest.in ) > "$P/selftest.log" 2>&1; rc_st=$?

stamp "summary: make_sav=$rc_sav make_native=$rc_nat cuda_doctor=$rc_doc pytest=$rc_pyt selftest=$rc_st" | tee "$P/summary.txt"
tail -n 3 "$P/pytest.log"
git -C "$AMPS" status --short > "$P/git_after.txt" 2>&1
[ $((rc_sav + rc_nat + rc_doc + rc_pyt + rc_st)) -eq 0 ]
