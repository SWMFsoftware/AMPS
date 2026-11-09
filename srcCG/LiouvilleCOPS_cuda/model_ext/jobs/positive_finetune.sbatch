#!/bin/bash
# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/positive_finetune.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/positive_finetune.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 0dd26efd173abf346f00a581eec6377075c727cd0b3b164ae8005442b14f3bc0
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. Per-element fine tuning (positive_finetune.py; user decisions 2026-11-07/08); FIT, OUT, FACETS (required: late_only |
#   all_seen), WINDOW_MIN, MAXNFEV (500), MODEL_EXT, LSMR_MAXITER (50), VALIDATE (1), CHECKPOINT_MIN; env COPS_POS_* from
#   the fit's manifest; exit 0 complete / 75 wall time (resubmit the same command) / 1-4 failure; OUT/final + top-level products.
# RUN: sbatch --time=02:00:00 --nodelist=gpu002 --exclude=gpu003 --export=ALL,FIT=<c k fit>/final,OUT=<c k>/ft,
#   FACETS=late_only,WINDOW_MIN=5,MAXNFEV=500,MODEL_EXT=<snapshot> positive_finetune.sbatch
# 2026-11-09 update: ratio 1e-5 (user 2026-11-08), CO2 flux floor (COPS_POS_FLOOR_CO2), start override --start-set (kappa = 1).
# ==== end of header ====
# Per-element fine tuning of a positive-activity 12M fit (model_ext/positive_finetune.py; user decisions 2026-11-07):
# every facet seen from a point the objective uses gets its own exponent increment delta (A = exp(s_fit + held +
# delta)), lower bound A >= FLOOR_RATIO * A_orig, no upper bound; the fit's full objective (data, W10x5 + grain in the
# late window, grain A, B refitted, crossing penalty, lambda/2 D^2 with the fit's lambda on the total ln A);
# temperatures, a, c, kappa held. Matrix-free Gauss-Newton (scipy least_squares TRF + LSMR, LinearOperator J) on all
# GPUs of the node, mixed precision (FP32 records, FP64 arithmetic), one bind of the fit's base (~150-160 GB host
# memory on 7-8 GPUs, ~3 min).
#
# Variables (sbatch --export=ALL,VAR=...):
#   FIT           the tuned fit's segment directory or its <fit>/final link (required; positive_fit.py products,
#                 restart status complete)
#   OUT           output directory (required; new, or an earlier run of THIS command to be continued, see below)
#   FACETS        the parameter facets (REQUIRED, no default; user 2026-11-08): late_only = facets with a structural
#                 record (Kn, KPhi or KV != 0) at an objective point inside the last WINDOW_MIN minutes (ending at the
#                 last selected data point) and none at any objective point before the window start (the old
#                 late-window fine tuning's set; earlier pressures unchanged exactly); all_seen = every facet seen
#                 from a point the objective uses (the 2026-11-07 procedure)
#   WINDOW_MIN    window length [min] for FACETS=late_only (required with it; user 2026-11-08: 5)
#   MAXNFEV       optimizer residual evaluations over all segments (default 500: user decision 2026-11-08)
#   LSMR_MAXITER  LSMR iterations per outer iteration (default 50; with x_scale 1 the method choices of the tool,
#                 recorded by the user 2026-11-08)
#   FLOOR_RATIO   r of the lower bound A >= r A_orig (default: the fit's COPS_POS_FLOOR_RATIO; a different value is
#                 refused by the driver: user 2026-11-07)
#   VALIDATE      1 (default): checks 1-6 in segment 1 first (cost at delta = 0 == the fit's, adjoint, J v / J^T u vs
#                 finite differences, bounds, timings), no fit if one fails; 0: skip them
#   VALIDATE_ONLY 1: the checks only (OUT/validate, exit 0 when they pass)
#   RESUME        finetune_state.npz / finetune_checkpoint.npz of another fine tuning of the same fit: start of segment 1
#   CHECKPOINT_MIN rolling OUT/seg<NN>/finetune_checkpoint.npz at most every N minutes while improving (default 30;
#                 0: none)
#   MODEL_EXT     the model_ext directory the job copies and runs (default $COPS_OUTPUT_ROOT/model_ext); production
#                 chains use a frozen snapshot (model_ext_frozen/<name>, never edited), so later edits cannot reach a
#                 pending job. The job copies positive_finetune.py, positive_problem.py, positive_activity.py,
#                 grain_crossing.py, test_positive_finetune.py and kernels/ into OUT/provenance-<job>/ at its start and
#                 runs from that copy (recorded in provenance-<job>/segment.txt and commands.sh).
# The COPS_POS_* environment is the fit's (FIT/fit_manifest.json model_extension.positive.environment); the driver
# refuses any difference.
# COPS_POS_FLOOR_CO2 (2026-11-08: the CO2 flux floor of the ratio-1e-5 fits) is exported only when the fit recorded
# it (else unset); the driver compares it like the other COPS_POS_* (positive_finetune.ENV_KEYS).
#
# Layout and segments (the positive_fit.sbatch convention): every job runs one segment OUT/seg<NN>/ (the driver's
# products). Wall time: the driver stops at its best evaluated state 15 min (+5 min margin) before the job's end, writes
# the segment's products (finetune_result.txt "stop wall_time", "complete no") and the job exits 75; resubmitting the
# SAME command (same FIT, OUT, MAXNFEV) continues: segment NN+1 starts from seg<NN>/finetune_state.npz with
# MAXNFEV - (evaluations of all earlier segments). A segment directory without finetune_state.npz (killed, failed) is
# moved aside as seg<NN>.incomplete_<UTC> and that segment runs again. When a segment finishes (stop tolerance or
# max_nfev) the job checks its finetune_state.npz (held_next, grain finite, complete) and links OUT/final -> seg<NN>
# and OUT/{finetune_state.npz, finetune_result.txt, COPS_model.dat, late_window_model.dat, pressure_comparison.txt,
# fine_tuning.json, validation.json} -> final/... (what the next cycle and the figure helpers read), then exits 0.
# Exit codes: 0 complete (final linked); 75 wall-time stop (resubmit to continue); anything else is a failure
# (2 usage, 3 inconsistent state or a bound violated, 4 product check failed, 1 driver error / failed validation).
# A completed OUT is not rerun (exit 0). Time limit: on the command line (sbatch --time=HH:MM:SS); default 1 h.
# Nodes: gpu001 / gpu002 only (user 2026-11-08: never gpu003; --exclude=gpu003 above, --nodelist on the command line).
# Example (the lambda 1e-2 fit with that fit's ratio 0.02):
#   F=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5/positive/C30_n4000_exp_cliff_w10x5_grain_cross10_lam1e-2/seg01
#   sbatch --time=05:00:00 --export=ALL,FIT=$F,OUT=<dir>,FACETS=late_only,WINDOW_MIN=5,MAXNFEV=500 positive_finetune.sbatch
#SBATCH --job-name=cops12m_posft
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --exclude=gpu003
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_posft-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
X=$MODEL_EXT; PY="$APP/.venv/bin/python -W ignore"
MODFILES=(positive_finetune.py positive_problem.py positive_activity.py grain_crossing.py test_positive_finetune.py)
PRODUCTS=(finetune_state.npz finetune_result.txt COPS_model.dat late_window_model.dat pressure_comparison.txt
          fine_tuning.json validation.json)
[ -n "${FIT:-}" ] || { echo "set FIT (the tuned fit's segment directory)" >&2; exit 2; }
[ -n "${OUT:-}" ] || { echo "set OUT (the output directory)" >&2; exit 2; }
FIT_GIVEN=$FIT; FIT=$(readlink -f "$FIT"); OUT=$(readlink -m "$OUT")
VALIDATE=${VALIDATE:-1}; VALIDATE_ONLY=${VALIDATE_ONLY:-0}; LSMR_MAXITER=${LSMR_MAXITER:-50}
RESUME=${RESUME:-}; CHECKPOINT_MIN=${CHECKPOINT_MIN:-30}; MAXNFEV=${MAXNFEV:-500}
FACETS=${FACETS:-}; WINDOW_MIN=${WINDOW_MIN:-}
case $FACETS in
    late_only) [ -n "$WINDOW_MIN" ] || { echo "FACETS=late_only needs WINDOW_MIN (e.g. 5)" >&2; exit 2; }
               FACETOPTS=(--facets late_only --window-min "$WINDOW_MIN");;
    all_seen) [ -z "$WINDOW_MIN" ] || { echo "WINDOW_MIN applies to FACETS=late_only only" >&2; exit 2; }
              FACETOPTS=(--facets all_seen);;
    *) echo "set FACETS=late_only (with WINDOW_MIN) or FACETS=all_seen (required since 2026-11-08; no default)" >&2; exit 2;;
esac
for f in "$FIT/fit_restart.npz" "$FIT/fit_manifest.json" "$FIT/activity_exp_facets.npz" "$X/kernels/finetune.cu" \
         "$X/kernels/pressure.cu" "$X/kernels/activity.cu"; do [ -e "$f" ] || { echo "missing: $f" >&2; exit 2; }; done
for f in "${MODFILES[@]}"; do [ -e "$X/$f" ] || { echo "missing: $X/$f" >&2; exit 2; }; done
# ---- the fit's GPU-patch environment (the driver checks it against the manifest again). COPS_POS_FLOOR_CO2
# (2026-11-08, the CO2 flux floor of the ratio-1e-5 fits) is exported only when the fit recorded it, else unset (the
# earlier fits and the 0.02 / 0.002 fits have none)
unset COPS_POS_FLOOR_CO2
ENVLINES=$($PY - "$FIT" <<'PY'
import json, sys
m = json.load(open(sys.argv[1] + "/fit_manifest.json"))
env = m["model_extension"]["positive"]["environment"]
for k in ("COPS_POS_ENABLE", "COPS_POS_NIGHT_RATIO", "COPS_POS_ACTIVITY", "COPS_POS_FLOOR_RATIO", "COPS_POS_FLOOR_CO2"):
    if k == "COPS_POS_FLOOR_CO2" and env.get(k) is None:
        continue
    v = str(env[k])
    assert v.replace(".", "", 1).replace("e-", "", 1).isalnum(), v
    print("export %s=%s" % (k, v))
PY
) || { echo "cannot read the COPS_POS_* environment of $FIT" >&2; exit 2; }
eval "$ENVLINES"
FLOOR_RATIO=${FLOOR_RATIO:-$COPS_POS_FLOOR_RATIO}
# ---- segment state of OUT: last seg<NN>, its stop, the evaluations used, its fit
segstate() {   # <OUT> <NN>: "stop nfev_total fit_restart_sha256" of seg<NN>/finetune_state.npz, or "nostate 0 -"
    $PY - "$1" "$2" <<'PY'
import json, os, sys
import numpy as np
p = os.path.join(sys.argv[1], "seg%02d" % int(sys.argv[2]), "finetune_state.npz")
if not os.path.isfile(p):
    print("nostate 0 -"); sys.exit()
with np.load(p, allow_pickle=False) as z:
    meta = json.loads(str(z["meta_json"]))
print(meta.get("stop", "unknown"), int(meta.get("nfev_total", meta.get("nfev", 0))), meta.get("fit_restart_sha256", "-"))
PY
}
FITSHA=$(sha256sum "$FIT/fit_restart.npz" | cut -c1-64)
LAST=0; LSTOP=none; LUSED=0; MOVED=""
if [ -e "$OUT" ]; then
    [ -d "$OUT" ] || { echo "OUT is not a directory: $OUT" >&2; exit 2; }
    for e in "$OUT"/*; do
        [ -e "$e" ] || [ -L "$e" ] || continue
        b=$(basename "$e")
        case $b in seg[0-9][0-9]|seg[0-9][0-9].incomplete_*|provenance-*|validate|final) ;;
            finetune_state.npz|finetune_result.txt|COPS_model.dat|late_window_model.dat|pressure_comparison.txt|fine_tuning.json|validation.json) ;;
            *) echo "OUT $OUT holds $b: not a directory of this job script (choose a new OUT)" >&2; exit 2;; esac
    done
    LAST=$(ls -d "$OUT"/seg[0-9][0-9] 2>/dev/null | sed 's/.*seg//' | sort -n | tail -1); LAST=$((10#${LAST:-0}))
fi
if [ "$LAST" -ge 1 ]; then
    read -r LSTOP LUSED LSHA < <(segstate "$OUT" "$LAST")
    if [ "$LSTOP" = nostate ]; then                      # killed or failed: move aside, run that segment again
        MOVED="$OUT/seg$(printf %02d "$LAST").incomplete_$(date -u +%Y%m%dT%H%M%SZ)"
        mv "$OUT/seg$(printf %02d "$LAST")" "$MOVED"
        LAST=$((LAST - 1)); LSTOP=none; LUSED=0
        [ "$LAST" -ge 1 ] && read -r LSTOP LUSED LSHA < <(segstate "$OUT" "$LAST")
    fi
    if [ "$LAST" -ge 1 ] && [ "$LSHA" != "$FITSHA" ]; then
        echo "$OUT/seg$(printf %02d "$LAST") tunes the fit restart $LSHA, FIT has $FITSHA: not the same fine tuning" >&2; exit 3
    fi
fi
final_link() {   # <NN>: check the segment's state, link OUT/final and the products
    local seg=seg$(printf %02d "$1")
    $PY - "$OUT/$seg/finetune_state.npz" <<'PY' || { echo "product check failed: $OUT/$seg/finetune_state.npz" >&2; return 4; }
import json, sys
import numpy as np
with np.load(sys.argv[1], allow_pickle=False) as z:
    meta = json.loads(str(z["meta_json"]))
    h, d, hp = z["held_next"], z["delta"], z["held_prev"]
    A, B = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
    ok = (meta.get("complete") is True and meta.get("stop") in ("tolerance", "max_nfev") and h.shape == d.shape == hp.shape
          and np.all(np.isfinite(h)) and np.array_equal(h, hp + d) and np.isfinite(A) and np.isfinite(B)
          and A >= 0 and B >= 0)
print("finetune_state.npz: stop %s, complete %s, %d evaluations, %d facets, grain A %.6g uPa B %.6g uPa/100s: %s"
      % (meta.get("stop"), meta.get("complete"), meta.get("nfev_total", -1), h.size, A, B, "OK" if ok else "FAILED"))
sys.exit(0 if ok else 1)
PY
    if [ -L "$OUT/final" ]; then
        [ "$(readlink "$OUT/final")" = "$seg" ] || { echo "$OUT/final -> $(readlink "$OUT/final"), expected $seg" >&2; return 3; }
    else ln -s "$seg" "$OUT/final"; fi
    for f in "${PRODUCTS[@]}"; do
        [ -e "$OUT/$seg/$f" ] || continue
        if [ -L "$OUT/$f" ]; then [ "$(readlink "$OUT/$f")" = "final/$f" ] || { echo "$OUT/$f is not final/$f" >&2; return 3; }
        else ln -s "final/$f" "$OUT/$f" || return 3; fi
    done
    echo "$(date -u +%FT%TZ) fine tuning complete: $OUT/final -> $seg"
}
if [ "$VALIDATE_ONLY" = 1 ]; then
    SEG=0; O=$OUT/validate; [ -e "$O" ] && { echo "exists: $O" >&2; exit 2; }
else
    case $LSTOP in
        none) SEG=1; LEFT=$MAXNFEV;;
        wall_time) SEG=$((LAST + 1)); LEFT=$((MAXNFEV - LUSED));;
        tolerance|max_nfev) echo "$(date -u +%FT%TZ) $OUT: segment $LAST finished ($LSTOP, $LUSED evaluations); nothing to do"
                            final_link "$LAST"; exit $?;;
        *) echo "$OUT/seg$(printf %02d "$LAST")/finetune_state.npz has stop '$LSTOP'; not continued" >&2; exit 3;;
    esac
    if [ "$LEFT" -le 0 ]; then
        echo "$(date -u +%FT%TZ) $OUT: $LUSED of $MAXNFEV evaluations used and the last segment stopped at the wall time;" \
             "nothing is left to run: the budget is used up but the state is marked incomplete (decide: raise MAXNFEV)" >&2
        exit 3
    fi
    O=$OUT/seg$(printf %02d "$SEG"); [ -e "$O" ] && { echo "exists: $O" >&2; exit 2; }
fi
mkdir -p "$OUT"; P=$OUT/provenance-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
# ---- provenance; the fine tuning runs from this copy of model_ext
cp -p "$0" "$P/positive_finetune.sbatch"
for f in "${MODFILES[@]}"; do cp -p "$X/$f" "$P/$f" || exit 1; done
cp -rp "$X/kernels" "$P/kernels" || exit 1
( cd "$P" && sha256sum "${MODFILES[@]}" kernels/* ) > "$P/model_ext.sha256" 2>&1
sha256sum "$FIT/fit_restart.npz" "$FIT/activity_exp_facets.npz" "$FIT/fit_manifest.json" > "$P/fit.sha256" 2>&1
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z0-9_]*|FIT|OUT|MAXNFEV|FACETS|WINDOW_MIN|LSMR_MAXITER|FLOOR_RATIO|VALIDATE[A-Z_]*|RESUME|CHECKPOINT_MIN|MODEL_EXT|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH|SLURM_JOB_GPUS|CUDA_VISIBLE_DEVICES)=' | sort > "$P/environment.txt"
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); AMPS $(git -C "$AMPS" rev-parse HEAD)"
  grep '^export ' "$0" | grep -v MODEL_EXT; echo "export MODEL_EXT=$X"
  echo "# the tuned fit's COPS_POS_* (fit_manifest.json model_extension.positive.environment)"
  echo "unset COPS_POS_FLOOR_CO2"; echo "$ENVLINES"; echo "cd $APP"; } > "$P/commands.sh"
rec() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; }
free -g > "$P/free_at_start.txt"; df -h "$OUT" > "$P/df_at_start.txt" 2>&1
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then NG=$(echo "$CUDA_VISIBLE_DEVICES" | tr , "\n" | wc -l); else NG=$(nvidia-smi -L | wc -l); fi
DEV=$(seq -s, 0 $((NG - 1)))
OPTS=()
if [ "$VALIDATE_ONLY" = 1 ]; then OPTS=(--validate-only)
else
    OPTS=(--max-nfev "$LEFT" --segment "$SEG" --nfev-before "$LUSED")
    [ "$SEG" = 1 ] && [ "$VALIDATE" = 1 ] && OPTS+=(--validate)
    if [ "$SEG" -gt 1 ]; then OPTS+=(--resume-from "$OUT/seg$(printf %02d "$LAST")/finetune_state.npz")
    elif [ -n "$RESUME" ]; then OPTS+=(--resume-from "$(readlink -f "$RESUME")"); fi
fi
{ echo "fine tuning of $FIT (given as $FIT_GIVEN; restart sha256 $FITSHA) -> $O"
  echo "model_ext copied and run: $X (MODEL_EXT) -> $P"
  echo "env: $(echo $ENVLINES | tr '\n' ' '); floor ratio $FLOOR_RATIO; LSMR_MAXITER $LSMR_MAXITER; FACETS $FACETS WINDOW_MIN ${WINDOW_MIN:--}"
  if [ "$VALIDATE_ONLY" = 1 ]; then echo "validation only"
  else echo "segment $SEG; previous segment $LAST ($LSTOP, $LUSED evaluations used); evaluations for this segment $LEFT of $MAXNFEV"; fi
  [ -n "$MOVED" ] && echo "moved aside: $MOVED"; } | tee "$P/segment.txt"
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): $NG GPUs ($DEV)"
# ---- GPU and host memory samples during the run
nvidia-smi --query-gpu=timestamp,index,memory.used,memory.total,utilization.gpu --format=csv -l 20 > "$P/nvidia_smi_samples.csv" 2>&1 &
SMI=$!
( while true; do echo "$(date -u +%FT%TZ) $(free -g | awk '/Mem:/ {print "used_GB", $3, "free_GB", $4, "available_GB", $7}')"; sleep 30; done ) > "$P/free_samples.txt" 2>&1 &
FREE=$!
CMD=($PY "$P/positive_finetune.py" --fit "$FIT" --out "$O" --floor-ratio "$FLOOR_RATIO" --devices "$DEV"
     --lsmr-maxiter "$LSMR_MAXITER" --wall-time slurm --checkpoint-minutes "$CHECKPOINT_MIN" "${FACETOPTS[@]}" "${OPTS[@]}")
rec "${CMD[@]}"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/finetune.log" 2>&1; rc=$?
kill $SMI $FREE 2>/dev/null
echo "$(date -u +%FT%TZ) fine tuning: exit $rc"
grep -E '^positive_finetune: (bound|parameters|FACETS|objective|the fit|check|timing|validation|start|final|products|checkpoint|stopped|complete|ERROR)' "$P/finetune.log" | cut -c1-320
tail -5 "$P/finetune.log" | cut -c1-250
free -g > "$P/free_at_end.txt"; df -h "$OUT" > "$P/df_at_end.txt" 2>&1; du -sh "$O" >> "$P/df_at_end.txt" 2>&1
sacct -j "$SLURM_JOB_ID" --format=JobID,State,Elapsed,MaxRSS,ExitCode > "$P/sacct.txt" 2>&1 || true
if [ $rc -eq 0 ] && [ "$VALIDATE_ONLY" != 1 ]; then
    read -r st _ _ < <(segstate "$OUT" "$SEG")
    if [ "$st" = tolerance ] || [ "$st" = max_nfev ]; then final_link "$SEG" || rc=$?
    else echo "segment $SEG exited 0 with stop '$st'" >&2; rc=3; fi
fi
[ $rc -eq 75 ] && echo "stopped before the wall-time limit: resubmit the same command (same FIT, OUT, MAXNFEV) to continue with segment $((SEG + 1))"
exit $rc
