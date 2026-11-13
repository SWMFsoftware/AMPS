#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/chain_onejob.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/chain_onejob.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 76f39cad56f81647bc0996188cd17ca9aaa18e4c1cb7d931a3645acbb1ee1fc8
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: submit helper. Submits one-job chains (jobs/chain_onejob.sbatch; user decision 2026-11-10): dry run by default,
#   --submit logs every job in model_ext/logs/submissions.txt; attempts with afterok:<after>,afternotok:A1:.. (an
#   upstream failure cancels every attempt); --kind nonneg --model softplus|squared, --l-threads, --start-file
#   (2026-11-11).
# RUN: chain_onejob.sh --kind nonneg --model softplus --chain-root <dir> --model-ext <frozen> --chain-code <frozen>
#   --cycles 3 --nfev 4000 --check 1 --ft-maxnfev 500 --ft-facets late_only --ft-window-min 5 --l-threads 32 --node
#   gpu001 --time 08:00:00 --attempts 3 --submit   (124162-124164; squared with --start-file on gpu002: 124165-124167)
# 2026-11-13: --slope-mult-min V (softplus, value 1); production 124664 (5 cycles), 124772 (extended to 9).
# ==== end of header ====
# Submit ONE-JOB fit / fine-tuning chains (jobs/chain_onejob.sbatch, model_ext/chain_onejob.py; user decision
# 2026-11-10 "yes, implement (b)"): the whole chain c1 fit -> c1 ft -> ... -> c<N> ft of one lambda in one job that
# binds the base's records once. The counterpart of jobs/positive_cycles.sh (one job per stage); same products, same
# layout <run>/positive_r<R>/lam<L>/c<k>/{fit,ft}, same stage settings.
#
# Default: DRY RUN (prints the sbatch lines, job ids as <placeholders>, and the chain's current stage states).
# --submit submits and logs each job in model_ext/logs/submissions.txt.
#
#   jobs/chain_onejob.sh --lambda 1e-4 --ratio 1e-5 --ft-maxnfev 500 --ft-facets late_only --ft-window-min 5 \
#       [--source R.npz] [--cycles 3] [--nfev 4000] [--check 1] [--c1-start-set "kappa=1"] [--model-ext DIR] \
#       [--chain-code DIR] [--chain-root DIR] [--node gpu002] [--gpus 8] [--time 08:00:00] [--attempts 3] \
#       [--sbatch "OPTS"] [--ft-export "FT_LSMR_MAXITER=50,FT_VALIDATE=1,FT_CHECKPOINT_MIN=30,FT_FLOOR_RATIO=r"] \
#       [--min-minutes "FIT FT"] [--after JOBID] [--submit]
#
# Attempts (as positive_cycles.sh): attempt a > 1 depends on afternotok of ALL earlier attempts (it runs only when every
# earlier one ended non-zero: 75 = wall time, the chain continues where it stopped; a real failure stops the chain again
# with its exit code, e.g. 3 for a fit segment without restart file, the CHAIN=1 rule); every job
# --kill-on-invalid-dep=yes. --after JOBID: EVERY attempt also depends on afterok:JOBID (attempt a > 1:
# "afterok:JOBID,afternotok:A1:...": if JOBID fails, all attempts are cancelled, none starts).
# --source default: positive_cycles.sh's (the first scan's fit of this lambda, last segment, status complete).
# --chain-root default: <run>/positive_r<R>/lam<L> (production); tests: a directory under $COPS_OUTPUT_ROOT/tests.
# --time: the job's limit (default 08:00:00 = the gpu partition's maximum; a longer chain resumes in the next attempt).
# Node: gpu003 is always excluded (user 2026-11-08: no fine tuning on gpu003); --node adds --nodelist.
# --kind slope | nonneg (2026-11-11, v6): the slope-multiplier families (KIND of jobs/chain_onejob.sbatch):
#   jobs/chain_onejob.sh --kind nonneg --model softplus|squared --chain-root DIR --model-ext <frozen> \
#       --ft-maxnfev 500 --ft-facets late_only --ft-window-min 5 --l-threads 32 [--nfev 4000] [--cycles 3] \
#       [--source <cycle-1 --start-restart, default D2 final>] [--start-file S (squared, nonneg_start.py)] \
#       [--slope-start X] [--slope-bundle B] [--ft-templates-from M] [--checkparams a:b] --node gpu001 --time T
#   --kind slope --model bumps: the 3a model (slope_mult_fit / slope_finetune).
#   --slope-mult-min 1 (nonneg softplus only; user 2026-11-12 refit: multipliers >= 1): lower bound of the free slope multipliers, every
#   cycle's fit (SLOPE_MULT_MIN); absent: unchanged.
#   --lambda / --ratio / --c1-start-set are positive-only; --chain-root is required for these kinds.
set -uo pipefail
ROOT_OUT=/data001/highresolu_vt/cops_output
J=$ROOT_OUT/jobs; RUN=$ROOT_OUT/runs/shap8_012m_lategrid_dz1.5; LOG=$ROOT_OUT/model_ext/logs/submissions.txt
PY="/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python -W ignore"
LAMBDA=""; RATIO=""; SOURCE=""; CYCLES=3; NFEV=4000; CHECK=1; C1_START_SET=""; MODEL_EXT=""; CHAIN_CODE=""
CHAIN_ROOT=""; NODE=""; GPUS=""; TIME="08:00:00"; ATTEMPTS=3; OPTS_EXTRA=""; FT_MAXNFEV=""; FT_FACETS=""
FT_WINDOW_MIN=""; FT_EXPORT=""; MINMIN=""; AFTER=""; SUBMIT=0
KIND=positive; SMODEL=""; SSTART=""; START_FILE=""; LT=""; LMODE=chunked; CHECKPARAMS=""; SBUNDLE=""; FTT=""; BUMPS=""
SMM=""
die() { echo "chain_onejob.sh: $*" >&2; exit 2; }
while [ $# -gt 0 ]; do
    case $1 in
        --lambda) LAMBDA=$2; shift 2;;          --ratio) RATIO=$2; shift 2;;
        --source) SOURCE=$2; shift 2;;          --cycles) CYCLES=$2; shift 2;;
        --nfev) NFEV=$2; shift 2;;              --check) CHECK=$2; shift 2;;
        --c1-start-set) C1_START_SET=$2; shift 2;;
        --model-ext) MODEL_EXT=$2; shift 2;;    --chain-code) CHAIN_CODE=$2; shift 2;;
        --chain-root) CHAIN_ROOT=$2; shift 2;;  --node) NODE=$2; shift 2;;
        --gpus) GPUS=$2; shift 2;;              --time) TIME=$2; shift 2;;
        --attempts) ATTEMPTS=$2; shift 2;;      --sbatch) OPTS_EXTRA=$2; shift 2;;
        --ft-maxnfev) FT_MAXNFEV=$2; shift 2;;  --ft-facets) FT_FACETS=$2; shift 2;;
        --ft-window-min) FT_WINDOW_MIN=$2; shift 2;;
        --ft-export) FT_EXPORT=$2; shift 2;;    --min-minutes) MINMIN=$2; shift 2;;
        --after) AFTER=$2; shift 2;;            --submit) SUBMIT=1; shift;;
        --kind) KIND=$2; shift 2;;              --model) SMODEL=$2; shift 2;;
        --slope-start) SSTART=$2; shift 2;;     --start-file) START_FILE=$2; shift 2;;
        --l-threads) LT=$2; shift 2;;           --l-mode) LMODE=$2; shift 2;;
        --checkparams) CHECKPARAMS=$2; shift 2;; --slope-bundle) SBUNDLE=$2; shift 2;;
        --ft-templates-from) FTT=$2; shift 2;;  --bumps) BUMPS=$2; shift 2;;
        --slope-mult-min) SMM=$2; shift 2;;
        -h|--help) sed -n '2,39p' "$0"; exit 0;;
        *) die "unknown argument $1 (--help)";;
    esac
done
case $KIND in positive|slope|nonneg) ;; *) die "--kind positive|slope|nonneg";; esac
if [ "$KIND" = positive ]; then
[ -n "$LAMBDA" ] || die "--lambda is required (verbatim, e.g. 1e-4)"
case $RATIO in 0.02|0.002|1e-5) ;; *) die "--ratio 0.02|0.002|1e-5 is required (the user's decisions)";; esac
else
[ -z "$LAMBDA$RATIO$C1_START_SET" ] || die "--lambda / --ratio / --c1-start-set are positive-kind options"
case $KIND/$SMODEL in slope/bumps|nonneg/softplus|nonneg/squared) ;;
    *) die "--kind $KIND --model $SMODEL: slope/bumps, nonneg/softplus or nonneg/squared";; esac
[ -n "$CHAIN_ROOT" ] || die "--kind $KIND: --chain-root is required"
[ "$FT_FACETS" = late_only ] || die "--kind $KIND: --ft-facets late_only"
case $LT in ""|*[!0-9]*) [ -z "$LT" ] || die "--l-threads N";; esac
case $LMODE in chunked|openblas) ;; *) die "--l-mode chunked|openblas";; esac
if [ "$SMODEL" = squared ] && [ -z "$START_FILE" ] && [ -z "$SSTART" ]; then die "--model squared: --start-file (nonneg_start.py)"; fi
if [ -n "$START_FILE" ] && [ ! -f "$START_FILE" ]; then
    [ "$SUBMIT" = 1 ] && die "--start-file $START_FILE does not exist (run jobs/nonneg_start.sbatch first)"
    echo "# NOTE: --start-file $START_FILE does not exist yet (jobs/nonneg_start.sbatch writes it): --submit refuses"
fi
fi
case $CYCLES in [1-9]) ;; *) die "--cycles must be 1..9";; esac
case $ATTEMPTS in [1-9]) ;; *) die "--attempts must be 1..9";; esac
case $CHECK in 0|1) ;; *) die "--check 0|1";; esac
[ -n "$FT_MAXNFEV" ] || die "--ft-maxnfev is required (the fine tuning's MAXNFEV; no default here)"
case $FT_FACETS in
    late_only) [ -n "$FT_WINDOW_MIN" ] || die "--ft-facets late_only needs --ft-window-min";;
    all_seen) [ -z "$FT_WINDOW_MIN" ] || die "--ft-window-min applies to late_only only";;
    *) die "--ft-facets late_only|all_seen is required";;
esac
case $TIME in [0-9]*:[0-9][0-9]:[0-9][0-9]|[0-9]*-[0-9]*:[0-9][0-9]:[0-9][0-9]) ;; *) die "--time HH:MM:SS";; esac
case $AFTER in ""|*[!0-9]*) [ -z "$AFTER" ] || die "--after takes ONE job id (Slurm cannot AND it with the attempts' afternotok when it is an OR list)";; esac
if [ -n "$C1_START_SET" ]; then
    C1_START_SET=${C1_START_SET//,/:}
    case $C1_START_SET in *[!A-Za-z0-9_=.:+-]*) die "--c1-start-set NAME=VALUE[:NAME=VALUE...]";; esac
fi
case $FT_EXPORT in "") ;; *[!A-Za-z0-9_=.,+-]*) die "--ft-export VAR=value,VAR=value";; esac
for kv in ${FT_EXPORT//,/ }; do
    case ${kv%%=*} in FT_LSMR_MAXITER|FT_VALIDATE|FT_CHECKPOINT_MIN|FT_FLOOR_RATIO) ;;
        *) die "--ft-export: only FT_LSMR_MAXITER, FT_VALIDATE, FT_CHECKPOINT_MIN, FT_FLOOR_RATIO";; esac
done
[ -z "$MODEL_EXT" ] || MODEL_EXT=$(readlink -f "$MODEL_EXT")
case $KIND in positive) NEED="positive_fit.py positive_finetune.py";; slope) NEED="slope_mult_fit.py slope_finetune.py slope_mult.py";;
    nonneg) NEED="nonneg_fit.py nonneg_finetune.py nonneg_mult.py";; esac
for f in $NEED; do
    if [ ! -f "${MODEL_EXT:-$ROOT_OUT/model_ext}/$f" ]; then
        [ "$SUBMIT" = 1 ] && die "MODEL_EXT ${MODEL_EXT:-$ROOT_OUT/model_ext} has no $f"
        echo "# NOTE: MODEL_EXT ${MODEL_EXT:-$ROOT_OUT/model_ext} has no $f (yet)"
    fi
done
CHAIN_CODE=$(readlink -f "${CHAIN_CODE:-$ROOT_OUT/model_ext}")
[ -f "$CHAIN_CODE/chain_onejob.py" ] || die "--chain-code $CHAIN_CODE has no chain_onejob.py"
if [ "$KIND" = positive ]; then CHAIN_ROOT=$(readlink -m "${CHAIN_ROOT:-$RUN/positive_r$RATIO/lam$LAMBDA}")
else CHAIN_ROOT=$(readlink -m "$CHAIN_ROOT"); fi
case $CHAIN_ROOT in "$RUN"/*|"$ROOT_OUT"/tests/*) ;; *) die "--chain-root under $RUN or $ROOT_OUT/tests";; esac
for f in chain_onejob.sbatch; do [ -f "$J/$f" ] || die "missing $J/$f"; done
if [ -z "$SOURCE" ] && [ "$KIND" != positive ]; then SOURCE=$RUN/clip_refit/D2/final/fit_restart.npz; fi   # D2 final
if [ -z "$SOURCE" ]; then                         # positive_cycles.sh's default
    FIRST=$RUN/positive/C30_n4000_exp_cliff_w10x5_grain_cross10_lam$LAMBDA
    SEGLAST=$(ls -d "$FIRST"/seg[0-9][0-9] 2>/dev/null | tail -1)
    [ -n "$SEGLAST" ] || die "no first-scan fit $FIRST/seg<NN> for lambda $LAMBDA (give --source)"
    SOURCE=$SEGLAST/fit_restart.npz
fi
case $SOURCE in /*) ;; *) SOURCE=$(readlink -m "$SOURCE");; esac
[ -f "$SOURCE" ] || die "missing cycle-1 source $SOURCE"
read -r SSTATUS SLAMBDA SCOST < <($PY - "$SOURCE" <<'PY'
import json, sys
import numpy as np
with np.load(sys.argv[1], allow_pickle=False) as z:
    m = json.loads(str(z["meta_json"]))
print(m.get("status"), (m.get("settings") or {}).get("lambda"), (m.get("optimizer") or {}).get("cost"))
PY
)
echo "# cycle-1 source $SOURCE: status $SSTATUS, lambda $SLAMBDA, final cost $SCOST"
[ "$SSTATUS" = complete ] || die "the cycle-1 source has status $SSTATUS, not complete"
if [ "$KIND" = positive ]; then
$PY -c "import sys; sys.exit(0 if float(sys.argv[1]) == float(sys.argv[2]) else 1)" "$SLAMBDA" "$LAMBDA" ||
    die "the cycle-1 source has lambda $SLAMBDA, not $LAMBDA"
VARS="ALL,LAMBDA=$LAMBDA,RATIO=$RATIO,SOURCE=$SOURCE,CYCLES=$CYCLES,NFEV=$NFEV,CHECK=$CHECK,FT_MAXNFEV=$FT_MAXNFEV"
VARS+=",FT_FACETS=$FT_FACETS,CHAIN_CODE=$CHAIN_CODE,CHAIN_ROOT=$CHAIN_ROOT"
else
VARS="ALL,KIND=$KIND,SLOPE_MODEL=$SMODEL,SOURCE=$SOURCE,CYCLES=$CYCLES,NFEV=$NFEV,CHECK=$CHECK,FT_MAXNFEV=$FT_MAXNFEV"
VARS+=",FT_FACETS=$FT_FACETS,CHAIN_CODE=$CHAIN_CODE,CHAIN_ROOT=$CHAIN_ROOT"
[ -n "$SSTART" ] && VARS+=",SLOPE_START=$SSTART"
[ -n "$START_FILE" ] && VARS+=",START_FILE=$(readlink -f "$START_FILE")"
[ -n "$LT" ] && VARS+=",L_THREADS=$LT,L_MODE=$LMODE"
[ -n "$CHECKPARAMS" ] && VARS+=",CHECKPARAMS=$CHECKPARAMS"
[ -n "$SBUNDLE" ] && VARS+=",SLOPE_BUNDLE=$(readlink -f "$SBUNDLE")"
[ -n "$FTT" ] && VARS+=",FT_TEMPLATES_FROM=$(readlink -f "$FTT")"
[ -n "$BUMPS" ] && VARS+=",BUMPS=$(readlink -f "$BUMPS")"
if [ -n "$SMM" ]; then
    [ "$KIND/$SMODEL" = nonneg/softplus ] || die "--slope-mult-min: nonneg softplus only (user 2026-11-12)"
    case $SMM in 1|1.0) ;; *) die "--slope-mult-min $SMM: the user's decision is 1";; esac
    VARS+=",SLOPE_MULT_MIN=$SMM"
fi
fi
[ -n "$FT_WINDOW_MIN" ] && VARS+=",FT_WINDOW_MIN=$FT_WINDOW_MIN"
[ -n "$MODEL_EXT" ] && VARS+=",MODEL_EXT=$MODEL_EXT"
[ -n "$C1_START_SET" ] && VARS+=",C1_START_SET=$C1_START_SET"
[ -n "$FT_EXPORT" ] && VARS+=",$FT_EXPORT"
if [ -n "$MINMIN" ]; then read -r MF MT <<< "$MINMIN"; VARS+=",MIN_FIT_MINUTES=$MF,MIN_FT_MINUTES=${MT:-$MF}"; fi
OPTS=(--kill-on-invalid-dep=yes --exclude=gpu003 --time="$TIME")
[ -n "$NODE" ] && OPTS+=(--nodelist="$NODE")
[ -n "$GPUS" ] && OPTS+=(--gpus-per-node="$GPUS")
if [ -n "$OPTS_EXTRA" ]; then read -r -a extra <<< "$OPTS_EXTRA"; OPTS+=("${extra[@]}"); fi
echo "# one-job $KIND chain ${SMODEL:+model $SMODEL }${LAMBDA:+for lambda $LAMBDA, ratio $RATIO}${LT:+ L_THREADS $LT ($LMODE)}: $CHAIN_ROOT/c<k>/{fit,ft}; $CYCLES cycles; fit NFEV $NFEV CHECK $CHECK; fine tuning MAXNFEV $FT_MAXNFEV FACETS $FT_FACETS ${FT_WINDOW_MIN:+WINDOW_MIN $FT_WINDOW_MIN}; MODEL_EXT ${MODEL_EXT:-(live model_ext)}; chain code $CHAIN_CODE; $ATTEMPTS attempts; $([ "$SUBMIT" = 1 ] && echo SUBMIT || echo 'DRY RUN')"
[ -n "$C1_START_SET" ] && echo "# cycle-1 fit START_SET=$C1_START_SET (cycles 2+ unchanged)"
echo "# current stage states (chain_onejob.py --plan-only):"
if [ "$KIND" = positive ]; then KARGS=(--lambda "$LAMBDA" --ratio "$RATIO")
else KARGS=(--kind "$KIND" --slope-model "$SMODEL"); fi
$PY "$CHAIN_CODE/chain_onejob.py" --plan-only --chain-root "$CHAIN_ROOT" "${KARGS[@]}" \
    --cycles "$CYCLES" --source "$SOURCE" --ft-maxnfev "$FT_MAXNFEV" --ft-facets "$FT_FACETS" \
    ${FT_WINDOW_MIN:+--ft-window-min "$FT_WINDOW_MIN"} --model-ext "${MODEL_EXT:-$ROOT_OUT/model_ext}" | sed 's/^/#   /'
ids=()
for a in $(seq 1 "$ATTEMPTS"); do
    # attempt 1: afterok:AFTER; attempt a > 1: afterok:AFTER AND afternotok of every earlier attempt (comma = AND).
    # Without the upstream condition on a > 1, an attempt 1 cancelled by --kill-on-invalid-dep (upstream failed) would
    # count as "not ok" and start attempt 2 (seen by the lead on Slurm 22.05.9, 2026-11-10).
    if [ "$a" = 1 ]; then d=${AFTER:+afterok:$AFTER}
    else d="${AFTER:+afterok:$AFTER,}afternotok:$(IFS=:; echo "${ids[*]}")"; fi
    if [ "$KIND" = positive ]; then JN="chain_l${LAMBDA}_r${RATIO}_$a"; else JN="chain_${SMODEL}_$a"; fi
    args=("${OPTS[@]}" ${d:+--dependency=$d} --job-name="$JN" --export="$VARS" "$J/chain_onejob.sbatch")
    if [ "$SUBMIT" = 1 ]; then
        JOBID=$(sbatch --parsable "${args[@]}") || die "sbatch failed for attempt $a"
        echo "attempt $a: job $JOBID"
        echo "$(date -u +%FT%TZ) $JOBID $J/chain_onejob.sbatch ${NODE:-any} chain_onejob.sh $KIND ${SMODEL:-lambda $LAMBDA ratio $RATIO} attempt $a: $(printf '%s ' "${args[@]}" | cut -c1-900)" >> "$LOG"
    else
        JOBID="<attempt$a>"
        printf '%s\n  sbatch' "# attempt $a"
        for x in "${args[@]}"; do case $x in *[?\<\>*]*) printf " '%s'" "$x";; *) printf ' %s' "$x";; esac; done; echo
    fi
    ids+=("$JOBID")
done
[ "$SUBMIT" = 1 ] || echo "# dry run: nothing submitted (add --submit)"
