#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_q002_preflight.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_q002_preflight.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: d0e3bea2effde4b73e546883f3306ab3cfd3f2229f8369eaab8c1dfaf8d405a0
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (CPU): the preflight of the q002 option C chains: the refitted D2 vs the 0.2-law D2 and base_q002
#   vs base (preflight_compare.py; stops when a rebuild is needed), the start files start_q002_templates_{on,off} of
#   the new D2 on the 0.02 design, test_localnn (8: reproduction of the new D2 by the 0.02 design).
# RUN 2026-11-18: planned (model_ext/logs/plan_q002_stage2r_C.md).
# ==== end of header ====
# Option C stage-2r chains of the new H2O law (q_min 0.02 q_max; lead / user 2026-11-18 ~02:40 CST): PREFLIGHT, after the
# D2 refit on base_q002 and before the five chains (CPU, ~15 min):
#   1. tests/localnn_q002/preflight_compare.py: the refitted D2 vs the 0.2-law D2 (data selection, references, ft templates,
#      enhancement scale, TS vs the design's), base_q002 vs base (kernels, elements, observation columns of COPS_model.dat)
#      -> preflight.json with rebuild_design / rebuild_late (the job stops with exit 3 when a rebuild is needed: a
#      decision; nothing is rebuilt silently)
#   2. the start files of the stage-1 start rule from the NEW D2: localnn_fit.py start --start d2 --templates on|off on the
#      0.02 design -> GRID_DIR/start_q002_templates_{on,off} (costs at the start, the lambda record)
#   3. test_localnn.py on the 0.02 design with the NEW D2 and base_q002: 3a / 3b design vs direct sums, 5 Jacobian (the
#      multi-step check), 8 reproduction of the NEW D2 (data cost and tables = the new-law fitted model, the decisive
#      check of the design); the c9 reproduction (0.2-law state) is expected to fail
# Settings: D2_DIR (required, the refitted D2 final), BASE_DIR (required, base_q002), MODEL_EXT (frozen localnn snapshot),
#   DESIGN_DIR (grid_150m/design_h2oqmin0.02), LATE (late_kernels_150m_h2oqmin0.02.npz), CONTRIB, OUT (preflight folder).
#SBATCH --job-name=cops_localnn_q002_pre
#SBATCH --partition=normal128cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=01:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_localnn_q002_pre-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 PYTHONUNBUFFERED=1
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; L=$R/nonneg/localnn; G=$L/grid_150m
export MODEL_EXT=$(readlink -f "${MODEL_EXT:?MODEL_EXT}")
D2N=$(readlink -f "${D2_DIR:?D2_DIR}"); BQ=$(readlink -f "${BASE_DIR:?BASE_DIR}")
DESIGN_DIR=$(readlink -f "${DESIGN_DIR:-$G/design_h2oqmin0.02}"); LATE=$(readlink -f "${LATE:-$L/late_kernels_150m_h2oqmin0.02.npz}")
CONTRIB=$(readlink -f "${CONTRIB:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_contrib_20261118T072729Z}")
C9=$R/nonneg/softplus/chain_5cycles_smm1_20261113T041215Z/c9
OUT=${OUT:-$COPS_OUTPUT_ROOT/tests/localnn_q002/preflight_$SLURM_JOB_ID}
PY=$APP/.venv/bin/python
stamp() { echo "$(date -u +%FT%TZ) $*"; }
FREE=$(df -P -B1G "$L" | awk 'NR==2 {print $4}'); [ "$FREE" -ge 20 ] || { stamp "disk guard: $FREE GB free < 20 GB"; exit 4; }
mkdir -p "$OUT" || exit 1; P=$OUT/provenance; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): D2_DIR=$D2N BASE_DIR=$BQ DESIGN_DIR=$DESIGN_DIR LATE=$LATE MODEL_EXT=$MODEL_EXT ($FREE GB free)"
cp -p "$0" "$P/localnn_q002_preflight.sbatch.used"; cp -p "$MODEL_EXT"/localnn_*.py "$MODEL_EXT"/test_localnn.py "$COPS_OUTPUT_ROOT/tests/localnn_q002/preflight_compare.py" "$P/"
( cd "$P" && sha256sum *.py localnn_q002_preflight.sbatch.used > code_sha256.txt )
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
echo "CPU job on $(hostname): no GPU used" > "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$PY" > "$P/pip_list.txt" 2>&1 || "$PY" -m pip list > "$P/pip_list.txt" 2>&1
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "export MODEL_EXT=$MODEL_EXT"; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; "$@"; }
# 1. compare
run "$PY" -W ignore "$COPS_OUTPUT_ROOT/tests/localnn_q002/preflight_compare.py" --old-d2 "$R/clip_refit/D2/final" --new-d2 "$D2N" \
  --old-base "$R/base" --new-base "$BQ" --facets-from "$C9/fit/seg01/nonneg_facets.npz" --design "$DESIGN_DIR" --late "$LATE" \
  --out "$OUT/preflight.json" > "$OUT/compare.log" 2>&1 || { stamp "compare failed"; tail -5 "$OUT/compare.log"; exit 1; }
cat "$OUT/compare.log" | cut -c1-600
NEED=$("$PY" -c "import json; j=json.load(open('$OUT/preflight.json')); print(int(j['rebuild_design']), int(j['rebuild_late']))")
if [ "$NEED" != "0 0" ]; then stamp "a rebuild is needed (design, late = $NEED): stop, see preflight.json (a decision)"; exit 3; fi
# 2. the start files of the new D2
COM=(--design "$DESIGN_DIR" --base "$BQ" --d2 "$D2N" --point-table "$R/late_grid_points.txt")
for t in on off; do
  [ -e "$G/start_q002_templates_$t" ] && { stamp "start_q002_templates_$t exists"; continue; }
  run "$PY" -W ignore "$MODEL_EXT/localnn_fit.py" start "${COM[@]}" --grid "$G" --templates $t --start d2 \
    --out "$G/start_q002_templates_$t" > "$OUT/start_$t.log" 2>&1 || { stamp "start $t failed"; tail -5 "$OUT/start_$t.log"; exit 1; }
  grep -o '"costs_at_start": {[^}]*}' "$OUT/start_$t.log" | head -1
done
# 3. tests on the 0.02 design with the new D2
run "$PY" -W ignore "$MODEL_EXT/test_localnn.py" --grid "$G" --design "$DESIGN_DIR" --base "$BQ" --d2 "$D2N" \
  --point-table "$R/late_grid_points.txt" --c9-ft "$C9/ft/final" --c9-fit "$C9/fit/final" --contrib-code "$CONTRIB" \
  --start-file "$G/start_q002_templates_off" --out "$OUT/test_localnn.json" > "$OUT/test_localnn.log" 2>&1
echo "test_localnn rc $? (the c9 reproduction is a 0.2-law state: expected to fail)"
grep -E "^PASS|^FAIL" "$OUT/test_localnn.log" | cut -c1-220
stamp "done"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_localnn_q002_pre-$SLURM_JOB_ID.out" "$P/job.log" 2>/dev/null
exit 0
