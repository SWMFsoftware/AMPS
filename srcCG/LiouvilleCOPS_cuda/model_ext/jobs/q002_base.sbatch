#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/q002_base_20261118/jobs/q002_base.sbatch (on Voyager: /data001/highresolu_vt/cops_output/tests/q002_base_20261118/jobs/q002_base.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 4088a1985b1619fc44f735aff428d1567feaf6670ffc92ab6c2f3b7283cd5b40
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (2026-11-18, Phase 0). Sample re-runs A-D of the lategrid base with the rebuilt C++, compare-sample, the
#   full new-law table (WriteKernels off), assembly gated on the comparison (job 136506: compare crashed in v1, so the
#   assembly was skipped; the full table ran).
# ==== end of header ====
#SBATCH --job-name=q002_base
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodelist=gpu002
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=64
#SBATCH --mem=200G
#SBATCH --time=03:30:00
#SBATCH --output=/data001/highresolu_vt/cops_output/tests/q002_base_20261118/logs/%x-%j.out
# Phase 0 of the H2O night-side fix (lead 2026-11-18 02:15 CST; user "(ii); refit D2 first, then all stage-2 chains"):
#   1 the REBUILT C++ (native/liouville_cops_cuda, job 136406: new Source.cpp) re-runs sample rows of
#     runs/shap8_012m_lategrid_dz1.5/base with the base's own input (boundary.in); only OutputDirectory and
#     FirstPoint / LastPoint / PointStep change: A 0..1541 step 67 (24 rows over the descent), B 1079-1080 (10:33:10,
#     10:33:10.282), C 1549-1550 (10:39:26, .282), D 927 (10:24:38, last NG before the gap)
#   2 nightfix_base_q002.py compare-sample: kernel records bitwise vs the base; C++ new-law pNG/pRG vs fit_cops' own
#     new-law model from the EXISTING kernels (fit_cops.main's consistency check, FP64)
#   3 the full new-law table: the same input with WriteKernels = off (all 1551 rows; no kernel files written)
#   4 only if 2 passed: nightfix_base_q002.py assemble -> runs/shap8_012m_lategrid_dz1.5/base_q002 (hard links +
#     the full run's table and manifest); refuses on any geometry/data column difference
set -u
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
OLD=$R/base NEW_BASE=$R/base_q002
W=$COPS_OUTPUT_ROOT/tests/q002_base_20261118
CODE=${CODE:?CODE (frozen nightfix_base_q002 snapshot) required}
BIN=$APP/native/liouville_cops_cuda
PY=$APP/.venv/bin/python
P=$W/provenance/job-$SLURM_JOB_ID; mkdir -p "$P" "$W/inputs" || exit 1
stamp() { echo "$(date -u +%FT%TZ) $*"; }
cp -p "$0" "$P/q002_base.sbatch.used"; cp -p "$CODE/nightfix_base_q002.py" "$P/"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1; module list > "$P/module_list.txt" 2>&1; lscpu > "$P/lscpu.txt" 2>&1
"$PY" -m pip list > "$P/pip_list.txt" 2>&1
sha256sum "$BIN" "$AMPS/srcCG/LiouvilleCOPS/src/Source.cpp" "$AMPS/srcCG/LiouvilleCOPS/tools/fit_cops.py" \
    "$APP/cops_cuda/evaluator.py" "$OLD/boundary.in" "$OLD/run_manifest.json" "$OLD/COPS_model.dat" > "$P/sha256.txt"
stamp "job $SLURM_JOB_ID on $(hostname)"

mkin() {  # $1 name, $2 first, $3 last, $4 step, $5 writekernels
  sed -e "s#^OutputDirectory *=.*#OutputDirectory = $W/cpp_$1#" -e "s#^FirstPoint *=.*#FirstPoint = $2#" \
      -e "s#^LastPoint *=.*#LastPoint  = $3#" -e "s#^PointStep *=.*#PointStep  = $4#" \
      -e "s#^WriteKernels *=.*#WriteKernels = $5#" "$OLD/boundary.in" > "$W/inputs/$1.in"
  diff "$OLD/boundary.in" "$W/inputs/$1.in" > "$W/inputs/$1.in.diff"
}
runcpp() {  # $1 name
  mkdir "$W/cpp_$1" || return 1
  stamp "C++ $1"
  ( cd "$W/cpp_$1" && env COPS_CUDA_DEVICE=0 "$BIN" "$W/inputs/$1.in" ) > "$W/cpp_$1/build.log" 2>&1
  local rc=$?; stamp "C++ $1 rc=$rc"; tail -1 "$W/cpp_$1/build.log"; return $rc
}
mkin sample_A 0 1541 67 on; mkin sample_B 1079 1080 1 on; mkin sample_C 1549 1550 1 on; mkin sample_D 927 927 1 on
mkin full_q002 0 1550 1 off
rc_s=0
for n in sample_A sample_B sample_C sample_D; do runcpp $n || rc_s=1; done
stamp "compare-sample"
cd "$APP" || exit 1
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$PY" -W ignore "$CODE/nightfix_base_q002.py" compare-sample \
    --old-base "$OLD" --new "$W/cpp_sample_A" --new "$W/cpp_sample_B" --new "$W/cpp_sample_C" --new "$W/cpp_sample_D" \
    --out "$W/sample_report.json"; rc_cmp=$?
stamp "compare-sample rc=$rc_cmp; full new-law table"
runcpp full_q002; rc_full=$?
rc_asm=skipped
if [ $rc_s -eq 0 ] && [ $rc_cmp -eq 0 ] && [ $rc_full -eq 0 ]; then
  stamp "assemble $NEW_BASE"
  "$PY" -W ignore "$CODE/nightfix_base_q002.py" assemble --old-base "$OLD" --cpp-run "$W/cpp_full_q002" \
      --cpp-input "$W/inputs/full_q002.in" --sample-report "$W/sample_report.json" --out "$NEW_BASE"; rc_asm=$?
fi
stamp "summary: samples=$rc_s compare=$rc_cmp full=$rc_full assemble=$rc_asm" | tee "$P/summary.txt"
exit 0
