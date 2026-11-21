#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_fit.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_fit.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: ae493f51a27da71e0016b0e034c22d4652e1e1244af5ca743544e10e448eb6fc
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (CPU, normal128cores, option C): one fit (TEMPLATES on|off, SMOOTH_FACTOR or SMOOTH_LAMBDA [+
#   SMOOTH_RULE_FILE], MULT_PENALTY_F (P_m, 2026-11-14), START, NFEV), wall-time stop / exit 75 / RUNROOT resume, then
#   localnn_direct.py; run root runs/.../nonneg/localnn/<RUN>_<UTC>.
# RUN: timing test 124995; start-rule stage 1 125014-125016 (cancelled); pilot 125071; anchored stage 1 125107, 125119,
#   125113 (node failure) -> 125177; P_m check 125360.
# 2026-11-18: 2026-11-18: DESIGN_DIR and CONTRIB defaults of the current law (design_h2oqmin0.02,
#   nonneg_contrib_20261118T072729Z).
# RUN 2026-11-18: no job since the change (no C fit started).
# ==== end of header ====
# Option C (local non-negative basis, user 2026-11-14): ONE fit on a CPU node (normal128cores), no GPU.
#   A_e = m_b(e) (sum_j c_j phi_j(e) [+ A_s sum_k ft_ck]),  c_j >= 0, m_[0,40) = 1, six multipliers >= 1, kappa >= 0,
#   grain A, B >= 0; D2's configuration otherwise (lategrid base, W10x5 from 10:33:10 UTC, grain at the time stamp,
#   crossing weight 10, original T law held, lambda 0 for the fitter's Tikhonov term); smoothness lambda_C x_c^T Lg x_c.
#   Model and cost: model_ext/localnn_problem.py; driver localnn_fit.py fit; after a complete fit localnn_direct.py
#   re-evaluates the final per-facet activity by direct facet sums at every row and checks the products.
# Inputs: the build of jobs/localnn_build.sbatch ($GRID_DIR with design/, start_templates_{off,on}/).
# Settings (sbatch --export=ALL,VAR=value,...):
#   TEMPLATES     on | off (required): D2's ft_c1..ft_c6 held at 1 additively inside the activity, or not at all
#   SMOOTH_FACTOR f (lambda_C = f * data cost at the start / x_c0^T Lg x_c0)   } exactly one of the two (required)
#   SMOOTH_LAMBDA lambda_C itself                                              }
#   SMOOTH_RULE_FILE  with SMOOTH_LAMBDA: the JSON of the rule that gave it (localnn_fit.py --smooth-rule-file; checked and
#                 recorded; 2026-11-14 the solution-anchored rule lambda_C = f C_ref / R_pilot)
#   NFEV          max_nfev (required; total over segments)
#   MULT_PENALTY_F  f_m of the multiplier penalty P_m = f_m C_ref sum_b (ln m_b)^2 (user 2026-11-14 ~17:30 CST;
#                 localnn_fit.py --mult-penalty-f; unset = off, the earlier fits unchanged)
#   RUN           run name (default localnn_t<TEMPLATES>_f<SMOOTH_FACTOR> or _l<SMOOTH_LAMBDA>); the run root is
#                 $R/nonneg/localnn/<RUN>_<UTC> (new) or RUNROOT (an existing run root: continue it after exit 75)
#   GRID_DIR      default $R/nonneg/localnn/grid_150m; START: default $GRID_DIR/start_templates_<TEMPLATES>
#   BLAS_THREADS  OpenBLAS threads of the dense SVD of the TRF step (default 16: fastest of 16/32/64/128 for the
#                 8516 x 2754 SVD on node032, 5.7 s; recorded)
#   WALL          seconds after which the segment stops at its best state and exits 75 (default: time limit - 20 min)
#   CHECK         1: analytic vs finite-difference Jacobian check before segment 1 (default 1)
#   MODEL_EXT     the code (production: a frozen snapshot), CONTRIB: the frozen contribution code
# Exit 75 = wall-time stop; resubmit with RUNROOT=<the run root> and the same settings to continue.
#SBATCH --job-name=cops_localnn_fit
#SBATCH --partition=normal128cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_localnn_fit-%j.out
set -uo pipefail
# 2026-11-18 (H2O night law q_min = 0.02 q_max): the products of the current law by default (DESIGN_DIR, CONTRIB, LATE,
# NONNEG_CODE); the loads refuse a design / late kernels of another law (localnn_problem.law_check)
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
BLAS_THREADS=${BLAS_THREADS:-16}
export OMP_NUM_THREADS=$BLAS_THREADS OPENBLAS_NUM_THREADS=$BLAS_THREADS MKL_NUM_THREADS=$BLAS_THREADS PYTHONUNBUFFERED=1
UNSET=$(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*|COPS_NONNEG_[A-Z_]*)' | sort | paste -sd' ' -)
[ -n "$UNSET" ] && unset $UNSET
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
BASE=$R/base; D2=$R/clip_refit/D2/final; PT=$R/late_grid_points.txt
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
CONTRIB=${CONTRIB:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_contrib_20261118T072729Z}
DESIGN_DIR=${DESIGN_DIR:-$GRID_DIR/design_h2oqmin0.02}
TEMPLATES=${TEMPLATES:-}; SMOOTH_FACTOR=${SMOOTH_FACTOR:-}; SMOOTH_LAMBDA=${SMOOTH_LAMBDA:-}; NFEV=${NFEV:-}
CHECK=${CHECK:-1}; RUNROOT=${RUNROOT:-}; SMOOTH_RULE_FILE=${SMOOTH_RULE_FILE:-}; MULT_PENALTY_F=${MULT_PENALTY_F:-}
GRID_DIR=${GRID_DIR:-$R/nonneg/localnn/grid_150m}
case $TEMPLATES in on|off) ;; *) echo "TEMPLATES must be on or off" >&2; exit 2;; esac
[ -n "$NFEV" ] || { echo "NFEV is required" >&2; exit 2; }
if [ -n "$SMOOTH_FACTOR" ] && [ -z "$SMOOTH_LAMBDA" ]; then SM=(--smooth-factor "$SMOOTH_FACTOR"); TAG=f$SMOOTH_FACTOR
elif [ -z "$SMOOTH_FACTOR" ] && [ -n "$SMOOTH_LAMBDA" ]; then SM=(--smooth-lambda "$SMOOTH_LAMBDA"); TAG=l$SMOOTH_LAMBDA
else echo "give exactly one of SMOOTH_FACTOR and SMOOTH_LAMBDA" >&2; exit 2; fi
if [ -n "$SMOOTH_RULE_FILE" ]; then [ -n "$SMOOTH_LAMBDA" ] || { echo "SMOOTH_RULE_FILE needs SMOOTH_LAMBDA" >&2; exit 2; }
  SM+=(--smooth-rule-file "$(readlink -f "$SMOOTH_RULE_FILE")"); fi
[ -n "$MULT_PENALTY_F" ] && SM+=(--mult-penalty-f "$MULT_PENALTY_F")
START=${START:-$GRID_DIR/start_templates_$TEMPLATES}
RUN=${RUN:-localnn_t${TEMPLATES}_${TAG}}
if [ -z "$RUNROOT" ]; then RUNROOT=$R/nonneg/localnn/${RUN}_$(date -u +%Y%m%dT%H%M%SZ); mkdir -p "$RUNROOT" || exit 1
else [ -d "$RUNROOT" ] || { echo "RUNROOT $RUNROOT does not exist" >&2; exit 2; }; fi
LIMIT=$(squeue -h -j "$SLURM_JOB_ID" -o %l 2>/dev/null)
secs() { local t=$1 d=0; [[ $t == *-* ]] && { d=${t%%-*}; t=${t#*-}; }; IFS=: read -r a b c <<< "$t"; [ -z "$c" ] && { c=$b; b=$a; a=0; }; echo $(( d*86400 + 10#$a*3600 + 10#$b*60 + 10#$c )); }
WALL=${WALL:-$(( $(secs "${LIMIT:-08:00:00}") - 1200 ))}
PY=$APP/.venv/bin/python
stamp() { echo "$(date -u +%FT%TZ) $*"; }
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
P=$RUNROOT/provenance-$SLURM_JOB_ID
mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): RUNROOT=$RUNROOT TEMPLATES=$TEMPLATES ${SM[*]} NFEV=$NFEV BLAS_THREADS=$BLAS_THREADS WALL=$WALL MODEL_EXT=$MODEL_EXT"

# provenance (CLAUDE.md rule 8)
cp -p "$0" "$P/localnn_fit.sbatch.used"
for f in localnn_mesh.py localnn_design.py localnn_problem.py localnn_fit.py localnn_direct.py test_localnn.py positive_problem.py; do
  cp -p "$MODEL_EXT/$f" "$P/" || exit 1
done
[ -f "$MODEL_EXT/SHA256SUMS" ] && cp -p "$MODEL_EXT/SHA256SUMS" "$P/MODEL_EXT_SHA256SUMS"
( cd "$P" && sha256sum *.py localnn_fit.sbatch.used > code_sha256.txt )
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded)" >> "$P/module_list.txt"
{ command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -q; } > "$P/nvidia_smi_q.txt" 2>&1 || echo "CPU job on $(hostname): no GPU used (nvidia-smi not available)" >> "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; echo "CSPICE: not used"; echo "nvcc: not used"; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$PY" > "$P/pip_list.txt" 2>&1 || "$PY" -m pip list > "$P/pip_list.txt" 2>&1
env | grep -E '^(SLURM_JOB_ID|SLURM_JOB_NODELIST|SLURM_CPUS_PER_TASK|MODEL_EXT|OMP_NUM_THREADS|OPENBLAS_NUM_THREADS)=' > "$P/env.txt"
echo "TEMPLATES=$TEMPLATES SMOOTH_FACTOR=$SMOOTH_FACTOR SMOOTH_LAMBDA=$SMOOTH_LAMBDA SMOOTH_RULE_FILE=$SMOOTH_RULE_FILE MULT_PENALTY_F=$MULT_PENALTY_F NFEV=$NFEV GRID_DIR=$GRID_DIR START=$START RUNROOT=$RUNROOT WALL=$WALL CHECK=$CHECK" >> "$P/env.txt"
{ echo '#!/bin/bash'; grep '^export ' "$0" | sed "s/\$BLAS_THREADS/$BLAS_THREADS/g"; echo "export MODEL_EXT=$MODEL_EXT"; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; "$@"; }
free -g > "$P/free_start.txt" 2>&1
JC=(); [ "$CHECK" = 1 ] && ! ls -d "$RUNROOT"/fit/seg[0-9][0-9] >/dev/null 2>&1 && JC=(--check-jacobian)
stamp "fit -> $RUNROOT/fit"
run heavy "$PY" -W ignore "$MODEL_EXT/localnn_fit.py" fit --design "$DESIGN_DIR" --base "$BASE" --d2 "$D2" \
  --point-table "$PT" --grid "$GRID_DIR" --templates "$TEMPLATES" --slope-mult-min 1 --start-file "$START" "${SM[@]}" \
  --max-nfev "$NFEV" --wall-time "$WALL" "${JC[@]}" --out "$RUNROOT/fit" > "$P/fit.log" 2>&1; rc=$?
tail -n 6 "$P/fit.log" | cut -c1-600
if [ $rc -eq 0 ] && [ -e "$RUNROOT/fit/final" ] && [ ! -e "$RUNROOT/fit/direct_check_$(readlink "$RUNROOT/fit/final")" ]; then
  stamp "direct facet-sum check of the final state"
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 run heavy "$PY" -W ignore "$MODEL_EXT/localnn_direct.py" \
    --seg "$RUNROOT/fit/$(readlink "$RUNROOT/fit/final")" --base "$BASE" --d2 "$D2" --point-table "$PT" \
    --design "$DESIGN_DIR" --contrib-code "$CONTRIB" --procs "${SLURM_CPUS_PER_TASK:-16}" > "$P/direct.log" 2>&1; rc=$?
  tail -n 2 "$P/direct.log" | cut -c1-800
fi
free -g > "$P/free_end.txt" 2>&1
stamp "done rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_localnn_fit-$SLURM_JOB_ID.out" "$P/job.log" 2>/dev/null
exit $rc
