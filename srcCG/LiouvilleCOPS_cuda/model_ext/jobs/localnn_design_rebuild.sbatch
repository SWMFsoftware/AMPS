#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_design_rebuild.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_design_rebuild.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 89a97c6e1d6342c3d4cfaa7176375b00041c1f1364e3a5a9d857fd3e733a4ad9
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (CPU): the design matrices and late kernels of option C for the current production law into new
#   folders (design_<TAG>, late_kernels_150m_<TAG>.npz; the old ones kept), then test_localnn / test_localnn_stage2
#   / test_localnn_multpen on them; disk guard.
# RUN 2026-11-18: 136507 (TAG h2oqmin0.02).
# ==== end of header ====
# Option C (lead 2026-11-18 02:15 CST, the H2O night-law fix q_min = 0.02 q_max): rebuild the design matrices and the
# late kernels with the contribution code of the current law into NEW folders (the 0.2-law ones are kept):
#   1. localnn_design.py -> GRID_DIR/design_h2oqmin0.02 (same grid, slope bins, seen / TS file and test states as the
#      0.2-law build of job 125004; the contribution code CONTRIB)
#   2. localnn_late.py   -> nonneg/localnn/late_kernels_150m_h2oqmin0.02.npz (late_only facets of NONNEG_CODE)
#   3. tests: test_localnn.py (the new design; 3a / 3b = design vs direct facet sums with the new law), test_localnn_stage2.py
#      (the new design and late kernels), test_localnn_multpen.py (the new design)
# Settings: MODEL_EXT (frozen localnn snapshot), CONTRIB (nonneg_contrib_20261118T072729Z), NONNEG_CODE
#   (nonneg_20261118T072729Z), TAG (h2oqmin0.02), STEPS (default "design late tests").
#SBATCH --job-name=cops_localnn_design
#SBATCH --partition=normal128cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=160G
#SBATCH --time=02:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_localnn_design-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; L=$R/nonneg/localnn; G=$L/grid_150m
export MODEL_EXT=$(readlink -f "${MODEL_EXT:?MODEL_EXT}")
CONTRIB=$(readlink -f "${CONTRIB:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_contrib_20261118T072729Z}")
NONNEG_CODE=$(readlink -f "${NONNEG_CODE:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_20261118T072729Z}")
TAG=${TAG:-h2oqmin0.02}; STEPS=${STEPS:-design late tests}
DOUT=$G/design_$TAG; LOUT=$L/late_kernels_150m_$TAG.npz
C9=$R/nonneg/softplus/chain_5cycles_smm1_20261113T041215Z/c9
PY=$APP/.venv/bin/python
stamp() { echo "$(date -u +%FT%TZ) $*"; }
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
FREE=$(df -P -B1G "$L" | awk 'NR==2 {print $4}'); [ "$FREE" -ge 20 ] || { stamp "disk guard: $FREE GB free < 20 GB"; exit 4; }
P=$L/provenance-$SLURM_JOB_ID-design_$TAG; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): MODEL_EXT=$MODEL_EXT CONTRIB=$CONTRIB NONNEG_CODE=$NONNEG_CODE TAG=$TAG STEPS='$STEPS' ($FREE GB free)"
cp -p "$0" "$P/localnn_design_rebuild.sbatch.used"; cp -p "$MODEL_EXT"/localnn_*.py "$MODEL_EXT"/test_localnn*.py "$P/"
mkdir "$P/contrib_code" && cp -p "$CONTRIB"/*.py "$P/contrib_code/"
( cd "$P" && sha256sum *.py contrib_code/*.py localnn_design_rebuild.sbatch.used > code_sha256.txt )
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded)" >> "$P/module_list.txt"
echo "CPU job on $(hostname): no GPU used" > "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$PY" > "$P/pip_list.txt" 2>&1 || "$PY" -m pip list > "$P/pip_list.txt" 2>&1
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "export MODEL_EXT=$MODEL_EXT"; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; "$@"; }
rc=0
for s in $STEPS; do
  case $s in
    design)
      [ -e "$DOUT" ] && { stamp "$DOUT exists: not rebuilt"; continue; }
      run heavy "$PY" -W ignore "$MODEL_EXT/localnn_design.py" --base "$R/base" --grid "$G" \
        --slope-bins "$R/slope_mult/bundle/slope_bins.npz" --facets-from "$C9/fit/seg01/nonneg_facets.npz" \
        --d2 "$R/clip_refit/D2/seg01" --c9-ft "$C9/ft/seg01" --c9-fit "$C9/fit/seg01" --contrib-code "$CONTRIB" \
        --out "$DOUT" --procs "${SLURM_CPUS_PER_TASK:-64}" > "$P/design.log" 2>&1; r=$?;;
    late)
      [ -e "$LOUT" ] && { stamp "$LOUT exists: not rebuilt"; continue; }
      run heavy "$PY" -W ignore "$MODEL_EXT/localnn_late.py" --base "$R/base" --d2 "$R/clip_refit/D2/final" \
        --point-table "$R/late_grid_points.txt" --contrib-code "$CONTRIB" --nonneg-code "$NONNEG_CODE" --out "$LOUT" \
        --procs "${SLURM_CPUS_PER_TASK:-64}" > "$P/late.log" 2>&1; r=$?;;
    tests)
      export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
      T1=$P/tests; mkdir -p "$T1"
      run "$PY" -W ignore "$MODEL_EXT/test_localnn.py" --grid "$G" --design "$DOUT" --base "$R/base" \
        --d2 "$R/clip_refit/D2/final" --point-table "$R/late_grid_points.txt" --c9-ft "$C9/ft/final" --c9-fit "$C9/fit/final" \
        --contrib-code "$CONTRIB" --start-file "$G/start_templates_off" --out "$T1/test_localnn.json" > "$T1/test_localnn.log" 2>&1
      echo "test_localnn rc $?"
      run "$PY" -W ignore "$MODEL_EXT/test_localnn_stage2.py" --grid "$G" --design "$DOUT" --base "$R/base" \
        --d2 "$R/clip_refit/D2/final" --point-table "$R/late_grid_points.txt" --late "$LOUT" --contrib-code "$CONTRIB" \
        --state "$L/stage1_anch_f10_pm0.01_cont_20261115T020243Z/fit/seg01" --rule "$L/lambda_rule_anchored_20261114.json" \
        --work "$T1/work_stage2" --out "$T1/test_localnn_stage2.json" > "$T1/test_localnn_stage2.log" 2>&1
      echo "test_localnn_stage2 rc $?"
      run "$PY" -W ignore "$MODEL_EXT/test_localnn_multpen.py" --grid "$G" --design "$DOUT" --base "$R/base" \
        --d2 "$R/clip_refit/D2/final" --point-table "$R/late_grid_points.txt" \
        --frozen "$COPS_OUTPUT_ROOT/model_ext_frozen/localnn_anch_20261114T090437Z" \
        --state "$L/stage1_anch_f10_20261114T094259Z/fit/final/localnn_state.npz" --out "$T1/test_localnn_multpen.json" \
        > "$T1/test_localnn_multpen.log" 2>&1
      echo "test_localnn_multpen rc $?"
      grep -h -E "^PASS|^FAIL|tests," "$T1"/*.log | cut -c1-220
      export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1; r=0;;
    *) stamp "unknown step $s"; r=2;;
  esac
  stamp "step $s rc $r"; tail -n 2 "$P/$s.log" 2>/dev/null | cut -c1-300; [ $r -eq 0 ] || rc=$r
done
du -sh "$DOUT" "$LOUT" 2>/dev/null
stamp "done rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_localnn_design-$SLURM_JOB_ID.out" "$P/job.log" 2>/dev/null
exit $rc
