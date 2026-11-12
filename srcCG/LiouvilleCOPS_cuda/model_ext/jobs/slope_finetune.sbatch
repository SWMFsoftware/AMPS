#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/slope_finetune.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/slope_finetune.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 68373490e131d3f46da6a800748aaa19bb98f6a0d09459a05a0dc18aa323cc53
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. Late-facet fine tuning of a slope-multiplier / bump fit (model_ext/slope_finetune.py; user decisions
#   2026-11-10): FACETS=late_only WINDOW_MIN=5, delta >= ln(F/s) (the floor as the bound), grain refitted, MAXNFEV 500,
#   L_THREADS, VALIDATE (checks 1-8); segments, exit 75 = resubmit.
# RUN: sbatch
#   --export=ALL,FIT=<c>/fit/final,OUT=<c>/ft,FACETS=late_only,WINDOW_MIN=5,MAXNFEV=500,L_THREADS=32,MODEL_EXT=<snapshot>
#   slope_finetune.sbatch   (3a chain 124064, 124073, 124082; 124015 failed validation; tests inside
#   slope_bumps_test.sbatch)
# ==== end of header ====
# Late-facet fine tuning of a slope-multiplier / bump + slope-multiplier fit (model_ext/slope_finetune.py; user decisions
# 2026-11-10: "add fine tuning -- fitting cycles as in other fits"): FACETS=late_only WINDOW_MIN=5 (the facets seen
# exclusively in the last 5 min; earlier pressures bitwise unchanged), A_e -> max(F_e, m_b L_e exp(h_e + delta_e)),
# delta_e >= ln(F_e / s_e) (the floor as the bound; late facets at the floor in the fit are not tuned), no cap, grain A, B
# refitted, lambda 0, the fit's full objective (W10x5, grain, crossing), positive_finetune's TRF + LSMR (x_scale 1,
# LSMR_MAXITER 50), MAXNFEV 500. G = dp/dA of the parameter facets is computed once on the GPU; every evaluation is a
# small host product (exact: pressures are linear in A at the held T and kappa).
#
# Variables (sbatch --export=ALL,VAR=...):
#   FIT          the tuned fit's segment directory or its <fit>/final link (required; slope_mult_fit.py products, complete)
#   OUT          output directory (required; new, or an earlier run of THIS command to be continued)
#   FACETS       late_only (required); WINDOW_MIN (required, user: 5)
#   MAXNFEV      optimizer evaluations over all segments (default 500, user 2026-11-08)
#   LSMR_MAXITER default 50
#   L_THREADS    chunked bitwise L thread pool of the fitter process (the GPU forward's L) and of the 12.4M-facet host
#                fields (default unset: 1 thread; the 3a production chain uses 32)
#   VALIDATE     1 (default): checks 1-8 in segment 1 first (start cost = the fit's to 1e-8, adjoint, J v vs FD, host model
#                vs a direct GPU evaluation, strict support, single facets, G vs a GPU FD, zero sensitivity), no fit if
#                one fails; 0: skip; VALIDATE_ONLY=1: the checks only (OUT/validate)
#   MODEL_EXT    the model_ext directory the job copies and runs (a frozen snapshot for chains)
# Layout and segments (positive_finetune.sbatch's convention): OUT/seg<NN>; a wall-time stop exits 75 and the SAME command
# continues with seg<NN+1> from seg<NN>/finetune_state.npz (MAXNFEV - the evaluations used); a segment without state is
# moved aside as seg<NN>.incomplete_<UTC>; on completion OUT/final -> seg<NN> and OUT/<product> -> final/<product>.
# Exit: 0 complete; 75 wall time; 2 usage; 3 inconsistent state; 4 product check failed; 1 driver error / validation.
#SBATCH --job-name=cops12m_slopeft
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
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_slopeft-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
UNSET=$(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*)' | sort | paste -sd' ' -)
[ -n "$UNSET" ] && unset $UNSET
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
X=$MODEL_EXT; PY="$APP/.venv/bin/python -W ignore"
MODFILES=(slope_finetune.py slope_mult.py slope_mult_problem.py positive_finetune.py positive_problem.py grain_crossing.py
          clip_problem.py test_slope_mult.py)
PRODUCTS=(finetune_state.npz finetune_result.txt COPS_model.dat late_window_model.dat pressure_comparison.txt
          fine_tuning.json validation.json)
[ -n "${FIT:-}" ] || { echo "set FIT" >&2; exit 2; }
[ -n "${OUT:-}" ] || { echo "set OUT" >&2; exit 2; }
FIT=$(readlink -f "$FIT"); OUT=$(readlink -m "$OUT")
VALIDATE=${VALIDATE:-1}; VALIDATE_ONLY=${VALIDATE_ONLY:-0}; LSMR_MAXITER=${LSMR_MAXITER:-50}; MAXNFEV=${MAXNFEV:-500}
FACETS=${FACETS:-}; WINDOW_MIN=${WINDOW_MIN:-}; L_THREADS=${L_THREADS:-}
[ "$FACETS" = late_only ] || { echo "set FACETS=late_only (required)" >&2; exit 2; }
[ -n "$WINDOW_MIN" ] || { echo "FACETS=late_only needs WINDOW_MIN (user: 5)" >&2; exit 2; }
for f in "$FIT/fit_restart.npz" "$FIT/fit_manifest.json" "$FIT/slope_mult_facets.npz"; do
    [ -e "$f" ] || { echo "missing: $f" >&2; exit 2; }; done
for f in "${MODFILES[@]}"; do [ -e "$X/$f" ] || { echo "missing: $X/$f" >&2; exit 2; }; done
segstate() {   # <OUT> <NN>: "stop nfev_total fit_restart_sha256" of seg<NN>/finetune_state.npz, or "nostate 0 -"
    $PY - "$1" "$2" <<'PY'
import json, os, sys
import numpy as np
p = os.path.join(sys.argv[1], "seg%02d" % int(sys.argv[2]), "finetune_state.npz")
if not os.path.isfile(p):
    print("nostate 0 -"); sys.exit()
with np.load(p, allow_pickle=False) as z:
    meta = json.loads(str(z["meta_json"]))
print(meta.get("stop", "unknown"), int(meta.get("nfev_total", 0)), (meta.get("fit") or {}).get("fit_restart_sha256", "-"))
PY
}
FITSHA=$(sha256sum "$FIT/fit_restart.npz" | cut -c1-64)
LAST=0; LSTOP=none; LUSED=0; MOVED=""
if [ -e "$OUT" ]; then
    for e in "$OUT"/*; do
        [ -e "$e" ] || [ -L "$e" ] || continue
        b=$(basename "$e")
        case $b in seg[0-9][0-9]|seg[0-9][0-9].incomplete_*|provenance-*|validate|final) ;;
            finetune_state.npz|finetune_result.txt|COPS_model.dat|late_window_model.dat|pressure_comparison.txt|fine_tuning.json|validation.json) ;;
            *) echo "OUT $OUT holds $b: not a directory of this job script" >&2; exit 2;; esac
    done
    LAST=$(ls -d "$OUT"/seg[0-9][0-9] 2>/dev/null | sed 's/.*seg//' | sort -n | tail -1); LAST=$((10#${LAST:-0}))
fi
if [ "$LAST" -ge 1 ]; then
    read -r LSTOP LUSED LSHA < <(segstate "$OUT" "$LAST")
    if [ "$LSTOP" = nostate ]; then
        MOVED="$OUT/seg$(printf %02d "$LAST").incomplete_$(date -u +%Y%m%dT%H%M%SZ)"
        mv "$OUT/seg$(printf %02d "$LAST")" "$MOVED"
        LAST=$((LAST - 1)); LSTOP=none; LUSED=0
        [ "$LAST" -ge 1 ] && read -r LSTOP LUSED LSHA < <(segstate "$OUT" "$LAST")
    fi
    if [ "$LAST" -ge 1 ] && [ "$LSHA" != "$FITSHA" ]; then
        echo "$OUT/seg$(printf %02d "$LAST") tunes the fit restart $LSHA, FIT has $FITSHA" >&2; exit 3
    fi
fi
final_link() {
    local seg=seg$(printf %02d "$1")
    $PY - "$OUT/$seg/finetune_state.npz" <<'PY' || { echo "product check failed: $OUT/$seg/finetune_state.npz" >&2; return 4; }
import json, sys
import numpy as np
with np.load(sys.argv[1], allow_pickle=False) as z:
    meta = json.loads(str(z["meta_json"]))
    h, d, hp = z["held_next"], z["delta"], z["held_prev"]
    A, B = float(z["grain_A_uPa"]), float(z["grain_B_uPa_per_100s"])
    ok = (meta.get("complete") is True and meta.get("stop") in ("tolerance", "max_nfev") and h.shape == d.shape == hp.shape
          and np.all(np.isfinite(h)) and np.array_equal(h, hp + d) and np.isfinite(A) and np.isfinite(B) and A >= 0 and B >= 0)
print("finetune_state.npz: stop %s, complete %s, %d evaluations, %d parameters, grain A %.6g uPa B %.6g: %s" % (
    meta.get("stop"), meta.get("complete"), meta.get("nfev_total", -1), meta.get("parameters", -1), A, B,
    "OK" if ok else "FAILED"))
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
RES=()
if [ "$VALIDATE_ONLY" = 1 ]; then
    SEG=0; O=$OUT/validate; [ -e "$O" ] && { echo "exists: $O" >&2; exit 2; }
    LEFT=0
else
    case $LSTOP in
        none) SEG=1; LEFT=$MAXNFEV;;
        wall_time) SEG=$((LAST + 1)); LEFT=$((MAXNFEV - LUSED)); RES=(--resume-from "$OUT/seg$(printf %02d "$LAST")/finetune_state.npz");;
        tolerance|max_nfev) echo "$(date -u +%FT%TZ) $OUT: segment $LAST finished ($LSTOP); nothing to do"; final_link "$LAST"; exit $?;;
        *) echo "$OUT/seg$(printf %02d "$LAST") has stop '$LSTOP'; not continued" >&2; exit 3;;
    esac
    [ "$LEFT" -le 0 ] && { echo "$OUT: $LUSED of $MAXNFEV evaluations used and incomplete: raise MAXNFEV" >&2; exit 3; }
    O=$OUT/seg$(printf %02d "$SEG"); [ -e "$O" ] && { echo "exists: $O" >&2; exit 2; }
fi
mkdir -p "$OUT"; P=$OUT/provenance-$SLURM_JOB_ID; mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/slope_finetune.sbatch"
for f in "${MODFILES[@]}"; do cp -p "$X/$f" "$P/$f" || exit 1; done
[ -f "$X/SHA256SUMS" ] && cp -p "$X/SHA256SUMS" "$P/"
(cd "$P" && sha256sum "${MODFILES[@]}") > "$P/run_copy.sha256" 2>&1
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1; nvidia-smi -L > "$P/nvidia_smi_L.txt" 2>&1
{ nvcc --version; gcc --version | head -1; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z_]*|FIT|OUT|FACETS|WINDOW_MIN|MAXNFEV|LSMR_MAXITER|L_THREADS|VALIDATE|VALIDATE_ONLY|MODEL_EXT|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PATH|LD_LIBRARY_PATH|SLURM_JOB_NODELIST)=' | sort > "$P/environment.txt"
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); AMPS $(git -C "$AMPS" rev-parse HEAD)"
  grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
{ echo "fine tuning of $FIT -> $OUT segment $SEG (previous $LAST: $LSTOP, $LUSED evaluations); evaluations $LEFT of $MAXNFEV"
  echo "FACETS $FACETS WINDOW_MIN $WINDOW_MIN; L_THREADS ${L_THREADS:-unset}; MODEL_EXT $X; node $(hostname)"
  [ -n "$MOVED" ] && echo "moved aside: $MOVED"; } | tee "$P/segment.txt"
free -g > "$P/free_at_start.txt"
NG=$(nvidia-smi -L | wc -l); DEV=$(seq -s, 0 $((NG - 1)))
VOPT=(); if [ "$VALIDATE_ONLY" = 1 ]; then VOPT=(--validate-only); elif [ "$SEG" = 1 ] && [ "$VALIDATE" = 1 ]; then VOPT=(--validate); fi
TOPT=(); [ -n "$L_THREADS" ] && TOPT=(--l-threads "$L_THREADS")
CMD=($PY "$P/slope_finetune.py" --fit "$FIT" --out "$O" --facets late_only --window-min "$WINDOW_MIN" --devices "$DEV"
     --lsmr-maxiter "$LSMR_MAXITER" --wall-time slurm --segment "$((SEG > 0 ? SEG : 1))" --nfev-before "$LUSED"
     "${VOPT[@]}" "${TOPT[@]}" "${RES[@]}")
[ "$VALIDATE_ONLY" = 1 ] || CMD+=(--max-nfev "$LEFT")
printf '%q ' bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" >> "$P/commands.sh"; echo >> "$P/commands.sh"
( while true; do echo "$(date -u +%T) $(free -m | awk '/Mem/{print $3, $7}')"; sleep 20; done ) > "$P/memory_used_avail_MB.txt" 2>&1 &
MEMLOG=$!
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${CMD[@]}" > "$P/finetune.log" 2>&1; rc=$?
kill $MEMLOG 2>/dev/null
echo "$(date -u +%FT%TZ) fine tuning segment $SEG: exit $rc"
grep -E '^slope_finetune|^slope_mult: (GPU patch installed .*role driver|L = Y a)|Traceback|Error' "$P/finetune.log" | cut -c1-500
tail -4 "$P/finetune.log" | cut -c1-300
free -g > "$P/free_at_end.txt"
echo "peak used host memory: $(sort -k2 -n "$P/memory_used_avail_MB.txt" | tail -1) (MB used, available)"
if [ $rc -eq 0 ] && [ "$VALIDATE_ONLY" != 1 ]; then final_link "$SEG" || rc=$?; fi
[ $rc -eq 75 ] && echo "stopped before the wall-time limit: resubmit the same command to continue"
exit $rc
