#!/bin/bash
# ==== version-controlled copy: header added 2026-11-21, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_slopetest_code_test.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_slopetest_code_test.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 403074d0e7ea9dcc16930939bbf1fdc67ff1d356bf06d1c7acfba1b1a49159ff
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: code test of localnn_fit.py --slope-test; no fit
# RUN: sbatch --export=ALL,CODE=<frozen code> jobs/localnn_slopetest_code_test.sbatch
# ==== end of header ====
#SBATCH --job-name=cops_slopetest_code
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_slopetest_code-%j.out
# Test of localnn_fit.py --slope-test (decision 3c, 2026-11-21): evaluations only, the optimizer is a stub and no fit
# runs. Then the direct facet-sum check of the lowered-bound test state. CODE = the frozen code to test.
set -euo pipefail
R=/data001/highresolu_vt/cops_output
A=/data001/highresolu_vt/AMPS
export AMPS=$A APP=$A/srcCG/LiouvilleCOPS_cuda COPS_OUTPUT_ROOT=$R
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 PYTHONUNBUFFERED=1 PYTHONDONTWRITEBYTECODE=1
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export LOCALNN_NONNEG_CODE=$R/model_ext_frozen/nonneg_20261118T072729Z
CODE=${CODE:?CODE (frozen code) is required}
OLD=$R/model_ext_frozen/localnn_q002_20261118T162327Z
N=$R/runs/shap8_012m_lategrid_dz1.5
L=$N/nonneg/localnn
REF=$L/stage2r_q002_f10_pm0.01_tploff_20261118T162407Z/c1/fit
DES=$L/grid_150m/design_h2oqmin0.02_d2q002
D2=$N/clip_refit/D2_q002_warm/seg01
PT=$N/late_grid_points.txt
PY=$APP/.venv/bin/python
O=$R/tests/localnn_slopetest_$SLURM_JOB_ID
mkdir "$O"
P=$O/provenance
mkdir "$P"
cp -p "$0" "$P/localnn_slopetest_code_test.sbatch.used"
cp -p "$CODE/SHA256SUMS" "$P/CODE_SHA256SUMS"
{ git -C $A rev-parse HEAD; git -C $A status --short; } > "$P/git.txt"
{ echo "host $(hostname) job $SLURM_JOB_ID $(date -Is)"; echo "CPU job: no GPU, nvcc or CSPICE used"; gcc --version | head -1; $PY --version; } > "$P/toolchain.txt"
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "no module command" >> "$P/module_list.txt"
$PY -c "import importlib.metadata as m; print('\n'.join(sorted(d.metadata['Name']+' '+d.version for d in m.distributions())))" > "$P/pip_list.txt"
cd "$APP"
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
heavy "$PY" -W ignore "$CODE/test_localnn_slopetest.py" --code "$CODE" --old-code "$OLD" --reference "$REF/seg02" --out "$O/cases" -- \
  --design "$DES" --base "$N/base_q002" --d2 "$D2" --point-table "$PT" --grid "$L/grid_150m" --templates off \
  --start-file "$REF/start" --smooth-lambda 8.249495752351283e-05 --smooth-rule-file "$L/lambda_rule_anchored_20261114.json" \
  --mult-penalty-f 0.01 --max-nfev 4000 --checkpoint-every 250
echo "direct facet-sum check of the lowered-bound test state (m = 0.7)"
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 heavy "$PY" -W ignore "$CODE/localnn_direct.py" --seg "$O/cases/free_min_TEST/seg01" \
  --base "$N/base_q002" --d2 "$D2" --point-table "$PT" --design "$DES" \
  --contrib-code $R/model_ext_frozen/nonneg_contrib_20261118T072729Z --procs "$SLURM_CPUS_PER_TASK"
echo "test finished rc 0"
