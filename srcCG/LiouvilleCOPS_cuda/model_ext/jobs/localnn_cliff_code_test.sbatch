#!/bin/bash
# ==== version-controlled copy: header added 2026-11-22, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_cliff_code_test.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_cliff_code_test.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 3291cc432d59f082f40713d16525f9605798c128bc8b5b4f11938cc05aebc68f
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: code test of --hold-mult BIN=V; no fit
# RUN: sbatch --export=ALL,TOOLS=<frozen cliff folder> <this file>
# ==== end of header ====
#SBATCH --job-name=cops_cliff_code
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_cliff_code-%j.out
# Code test of localnn_fit.py --hold-mult BIN=V (cliff test, 2026-11-22): evaluations only, the optimizer is a stub.
# The slope-test suite against the previous frozen code (default path identical), then test_localnn_cliff.py. Env: TOOLS.
set -euo pipefail
R=/data001/highresolu_vt/cops_output
A=/data001/highresolu_vt/AMPS
export AMPS=$A APP=$A/srcCG/LiouvilleCOPS_cuda COPS_OUTPUT_ROOT=$R
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export LOCALNN_NONNEG_CODE=$R/model_ext_frozen/nonneg_20261118T072729Z
TOOLS=${TOOLS:?TOOLS is required}
OLD=$R/model_ext_frozen/localnn_slopetest_20261122T033334Z
N=$R/runs/shap8_012m_lategrid_dz1.5
L=$N/nonneg/localnn
REF=$L/stage2r_q002_f10_pm0.01_tploff_20261118T162407Z/c1/fit
PY=$APP/.venv/bin/python
O=$R/tests/localnn_cliff_code_$SLURM_JOB_ID
mkdir "$O" "$O/provenance"
P=$O/provenance
cp -p "$0" "$P/localnn_cliff_code_test.sbatch.used"
cp -p "$TOOLS/code/SHA256SUMS" "$P/CODE_SHA256SUMS"
{ git -C $A rev-parse HEAD; git -C $A status --short; } > "$P/git.txt"
{ echo "host $(hostname) job $SLURM_JOB_ID $(date -Is)"; echo "CPU job: no GPU, nvcc or CSPICE used"; gcc --version | head -1; $PY --version; } > "$P/toolchain.txt"
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "no module command" >> "$P/module_list.txt"
$PY -c "import importlib.metadata as m; print('\n'.join(sorted(d.metadata['Name']+' '+d.version for d in m.distributions())))" > "$P/pip_list.txt"
(cd "$TOOLS/code" && sha256sum -c --quiet SHA256SUMS)
cd "$APP"
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
COM=(--design "$L/grid_150m/design_h2oqmin0.02_d2q002" --base "$N/base_q002" --d2 "$N/clip_refit/D2_q002_warm/seg01"
     --point-table "$N/late_grid_points.txt" --grid "$L/grid_150m" --templates off --start-file "$REF/start"
     --smooth-lambda 8.249495752351283e-05 --smooth-rule-file "$L/lambda_rule_anchored_20261114.json"
     --mult-penalty-f 0.01 --max-nfev 4000 --checkpoint-every 250)
heavy "$PY" -W ignore "$TOOLS/code/test_localnn_cliff.py" --code "$TOOLS/code" --old-code "$OLD" --reference "$REF/seg02" --out "$O/cliff" -- "${COM[@]}"
heavy "$PY" -W ignore "$TOOLS/code/test_localnn_slopetest.py" --code "$TOOLS/code" --old-code "$OLD" --reference "$REF/seg02" --out "$O/slopetest" -- "${COM[@]}"
echo "code test finished rc 0"
