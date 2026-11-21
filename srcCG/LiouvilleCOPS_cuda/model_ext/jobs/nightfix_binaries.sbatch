#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/nightfix_20261118/jobs/nightfix_build.sbatch (on Voyager: /data001/highresolu_vt/cops_output/tests/nightfix_20261118/jobs/nightfix_build.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 80bedbdd0594351ff98d8576dd49cacacb064b23955d1986822224933e9744e8
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (2026-11-18, H2O night-side fix). Rebuild of liouville_cops, spice_angular_velocity, native/liouville_cops_cuda
#   on gpu002 (old binaries moved aside), GetQ probe, selftests, A6 gate, cuda_doctor, LiouvilleCOPS_cuda pytest.
# RUN: sbatch tests/nightfix_20261118/jobs/nightfix_build.sbatch (job 136406).
# ==== end of header ====
#SBATCH --job-name=nightfix_build
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodelist=gpu002
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/tests/nightfix_20261118/jobs/%x-%j.out
# H2O night-side fix (user decision 2026-11-18: "(b) -- use 0.02; x10 was not intended"; lead 2026-11-18): rebuild every
# binary that compiles LiouvilleCOPS/src/Source.cpp, on gpu002 (-march=native on the GPU nodes' CPU, AMD EPYC 7513 =
# the earlier build_check.sbatch jobs 120810/120814), then the checks:
#   1 provenance (CLAUDE.md rule 8)
#   2 old binaries moved aside (never deleted) into cops_output/tests/backups_before_20261118_nightfix/binaries/moved_<UTC>
#   3 make ../LiouvilleCOPS liouville_cops spice_angular_velocity (CSPICE=); make native (CSPICE= CUDA= ARCH=sm_80)
#   4 GetQ probe against the rebuilt LiouvilleCOPS objects: H2O qmin/qmax = 0.02
#   5 C++ selftests (CPU liouville_cops, native liouville_cops_cuda)
#   6 LiouvilleCOPS gate tools/run_a6_validation.py (unittest suite + selftest), then --integration
#   7 cuda_doctor, LiouvilleCOPS_cuda pytest (GPU)
# Every step runs even if an earlier one fails; return codes in summary.txt. Products in $P.
set -u
export AMPS=/data001/highresolu_vt/AMPS
export APP=$AMPS/srcCG/LiouvilleCOPS_cuda REF=$AMPS/srcCG/LiouvilleCOPS
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
CSPICE=/home/vtenishev/SPICE/cspice
W=$COPS_OUTPUT_ROOT/tests/nightfix_20261118
P=$W/build/job-$SLURM_JOB_ID; mkdir -p "$P" || exit 1
UTC=$(date -u +%Y%m%dT%H%M%SZ)
MOVED=$COPS_OUTPUT_ROOT/tests/backups_before_20261118_nightfix/binaries/moved_$UTC
cp -p "$0" "$P/nightfix_build.sbatch.used"
stamp() { echo "$(date -u +%FT%TZ) $*"; }
stamp "start on $(hostname), products in $P"

# 1 provenance
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
git -C "$AMPS" diff -- srcCG/LiouvilleCOPS/src srcCG/LiouvilleCOPS/tools/fit_cops.py srcCG/LiouvilleCOPS_cuda/cops_cuda > "$P/git_diff_code.txt" 2>&1
module list > "$P/module_list.txt" 2>&1
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
lscpu > "$P/lscpu.txt" 2>&1
{ nvcc --version; gcc --version | head -1; g++ --version | head -1; grep -m1 CSPICE $CSPICE/doc/version.txt; mpicc --showme:version 2>&1; } > "$P/toolchain.txt" 2>&1
"$APP/.venv/bin/python" -m pip list > "$P/pip_list.txt" 2>&1 || ~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
sha256sum "$REF/src/Source.cpp" "$REF/tools/fit_cops.py" "$APP/cops_cuda/evaluator.py" > "$P/sources_sha256.txt"

# 2 old binaries aside (copies were also made beforehand in backups_before_20261118_nightfix/binaries/)
mkdir -p "$MOVED/LiouvilleCOPS" "$MOVED/native" || exit 1
for f in liouville_cops spice_angular_velocity; do [ -e "$REF/$f" ] && mv "$REF/$f" "$MOVED/LiouvilleCOPS/"; done
[ -d "$REF/build" ] && mv "$REF/build" "$MOVED/LiouvilleCOPS/build"
[ -e "$APP/native/liouville_cops_cuda" ] && mv "$APP/native/liouville_cops_cuda" "$MOVED/native/"
[ -d "$APP/native/build" ] && mv "$APP/native/build" "$MOVED/native/build"
( cd "$MOVED" && find . -type f | sort | xargs sha256sum ) > "$P/moved_aside_sha256.txt" 2>&1
echo "$MOVED" > "$P/moved_aside_dir.txt"

# 3 builds
cd "$APP" || exit 1
stamp "make liouville_cops spice_angular_velocity"
make -C ../LiouvilleCOPS CSPICE=$CSPICE liouville_cops spice_angular_velocity > "$P/make_liouville_cops.log" 2>&1; rc_cpu=$?
stamp "make liouville_cops_cuda (sm_80)"
make -C native CSPICE=$CSPICE CUDA=$CUDA_HOME ARCH=sm_80 > "$P/make_liouville_cops_cuda.log" 2>&1; rc_nat=$?
ls -l ../LiouvilleCOPS/liouville_cops ../LiouvilleCOPS/spice_angular_velocity native/liouville_cops_cuda > "$P/executables.txt" 2>&1
sha256sum ../LiouvilleCOPS/liouville_cops ../LiouvilleCOPS/spice_angular_velocity native/liouville_cops_cuda >> "$P/executables.txt" 2>&1

# 4 GetQ probe (the CPU objects; the native objects build/ref/Source.o are compiled from the same Source.cpp)
stamp "GetQ probe"
OBJS=$(ls ../LiouvilleCOPS/build/*.o | grep -v '/main\.o$')
g++ -O2 -std=c++17 -fopenmp -I$CSPICE/include -I../LiouvilleCOPS/src -o "$P/getq_probe" "$W/build/getq_probe.cpp" $OBJS \
    $CSPICE/lib/cspice.a -lm > "$P/getq_probe_build.log" 2>&1 && "$P/getq_probe" > "$P/getq_probe.txt" 2>&1; rc_probe=$?
OBJS_N=$(ls native/build/ref/*.o | grep -v '/main\.o$')
g++ -O2 -std=c++17 -fopenmp -I$CSPICE/include -I../LiouvilleCOPS/src -o "$P/getq_probe_native_objects" "$W/build/getq_probe.cpp" \
    $OBJS_N native/build/access_shim.o native/build/gpu_liouville.o $CSPICE/lib/cspice.a -L$CUDA_HOME/lib64 -lcudart_static \
    -ldl -lrt -lpthread -lm > "$P/getq_probe_native_build.log" 2>&1 && "$P/getq_probe_native_objects" > "$P/getq_probe_native_objects.txt" 2>&1; rc_probe_n=$?
cat "$P/getq_probe.txt" "$P/getq_probe_native_objects.txt" 2>/dev/null

# 5 selftests
stamp "CPU selftest"
( cd ../LiouvilleCOPS && env OMP_NUM_THREADS=16 ./liouville_cops selftest.in ) > "$P/selftest_cpu.log" 2>&1; rc_st_cpu=$?
stamp "native selftest"
( cd ../LiouvilleCOPS && env COPS_CUDA_DEVICE=0 OMP_NUM_THREADS=16 ../LiouvilleCOPS_cuda/native/liouville_cops_cuda selftest.in ) > "$P/selftest_native.log" 2>&1; rc_st_nat=$?

# 6 LiouvilleCOPS gate (unittest discover + selftest), then the integration gate
stamp "LiouvilleCOPS A6 gate"
( cd ../LiouvilleCOPS && "$APP/.venv/bin/python" tools/run_a6_validation.py --threads 16 --report "$P/a6_validation.json" ) > "$P/a6_validation.log" 2>&1; rc_a6=$?
stamp "LiouvilleCOPS A6 gate --integration"
( cd ../LiouvilleCOPS && timeout 3600 "$APP/.venv/bin/python" tools/run_a6_validation.py --integration --threads 16 --report "$P/a6_validation_integration.json" ) > "$P/a6_validation_integration.log" 2>&1; rc_a6i=$?

# 7 GPU application checks
stamp "cuda_doctor"
.venv/bin/python tools/cuda_doctor.py --json "$P/cuda_doctor.json" > "$P/cuda_doctor.log" 2>&1; rc_doc=$?
stamp "pytest LiouvilleCOPS_cuda/tests"
.venv/bin/python -m pytest -q -rs tests > "$P/pytest.log" 2>&1; rc_pyt=$?

stamp "summary: make_cpu=$rc_cpu make_native=$rc_nat getq_probe=$rc_probe getq_probe_native=$rc_probe_n selftest_cpu=$rc_st_cpu selftest_native=$rc_st_nat a6=$rc_a6 a6_integration=$rc_a6i cuda_doctor=$rc_doc pytest=$rc_pyt" | tee "$P/summary.txt"
tail -n 3 "$P/pytest.log"
git -C "$AMPS" status --short > "$P/git_after.txt" 2>&1
exit 0
