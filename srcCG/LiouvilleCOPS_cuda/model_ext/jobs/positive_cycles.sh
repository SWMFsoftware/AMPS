#!/bin/bash
# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/positive_cycles.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/positive_cycles.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: e22a7af52f01ba58def3e4e2433f4e0beb3f08fa4d86db11feeed3f033e6aab5
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job chain. Fit / fine-tuning cycles for one lambda (user decisions 2026-11-07/08): c1 fit -> c1 ft -> c2 fit -> c2 ft -> c3 fit
#   -> c3 ft with sbatch dependencies (attempts with afternotok, next stage afterok:A1?afterok:A2..., --kill-on-invalid-dep);
#   the next fit holds the fine tuning's held_next and takes its grain. Dry run by default.
# RUN: positive_cycles.sh --lambda 1e-3 --ratio 0.002 --model-ext <c1 snapshot> --ft-model-ext <ft snapshot> --fit-node gpu002
#   --ft-node gpu002 --start-at c1ft --ft-maxnfev 500 --ft-time 02:00:00 --ft-export "FACETS=late_only,WINDOW_MIN=5" --submit
# 2026-11-09 update: ratio 1e-5 (user 2026-11-08), CO2 flux floor (COPS_POS_FLOOR_CO2), start override --start-set (kappa = 1).
# ==== end of header ====
# Fit / fine-tuning cycles of the positive-activity 12M fit for ONE lambda (user decisions 2026-11-07: "repeat all 6 fits
# starting from the best restart of each of the corresponding fit", "3 cycles of fitting and fine tuning for each of the
# fit"; in cycles 2 and 3 the fine-tuning change is HELD as a fixed per-element term in the exponent, accumulated over
# cycles, and D measures the TOTAL ln A from s0; grain A, B of the next fit start from the previous fine tuning).
#
#   c1 fit -> c1 ft -> c2 fit -> c2 ft -> c3 fit -> c3 ft     (sbatch --dependency, --kill-on-invalid-dep=yes)
#
# Default: DRY RUN (prints every sbatch command, job ids as <placeholders>). --submit submits and logs each job in
# model_ext/logs/submissions.txt.
#
#   jobs/positive_cycles.sh --lambda 1e-2 [--ratio 0.002] [--model-ext DIR] [--ft-model-ext DIR]
#       [--fit-node gpu001] [--ft-node gpu002] [--fit-gpus 8] [--ft-gpus G] [--fit-sbatch "OPTS"] [--ft-sbatch "OPTS"]
#       [--fit-attempts 3] [--ft-attempts 4] [--nfev 4000] [--check 1] [--cycles 3] [--source R.npz]
#       [--start-at c1fit|c1ft|c2fit|c2ft|c3fit|c3ft] [--after JOBID]
#       --ft-maxnfev N --ft-time HH:MM:SS [--ft-export "VAR=v,VAR=v"] [--c1-start-set "NAME=VALUE"] [--submit]
#
# Products (run = runs/shap8_012m_lategrid_dz1.5, L = --lambda verbatim, R = --ratio):
#   <run>/positive_r<R>/lam<L>/c<k>/fit/seg<NN>   fit of cycle k (jobs/positive_fit.sbatch CYCLE=k); <fit>/final -> seg<NN>
#                                                 when the fit is complete (written by positive_fit.sbatch)
#   <run>/positive_r<R>/lam<L>/c<k>/ft/seg<NN>    fine tuning of cycle k (jobs/positive_finetune.sbatch FIT=<fit>/final
#                                                 OUT=c<k>/ft); on completion it links c<k>/ft/final -> seg<NN> and the
#                                                 products c<k>/ft/finetune_state.npz (held_next = held_prev + delta,
#                                                 grain_A_uPa, grain_B_uPa_per_100s), finetune_result.txt, ...
# Cycle 1 fit: START_RESTART = --source (default: the first scan's fit of this lambda, runs/.../positive/
#   C30_n4000_exp_cliff_w10x5_grain_cross10_lam<L>, its last segment's fit_restart.npz, which must have status complete;
#   every first-scan fit has one segment, seg01), no held term, grain from that restart.
# Cycle k > 1 fit: START_RESTART = c<k-1>/fit/final/fit_restart.npz (a, c, kappa; prior s0 unchanged),
#   HELD = GRAIN_FROM = c<k-1>/ft/finetune_state.npz (the top-level link written on completion: held_next; grain).
#
# Wall time and retries (fits and fine tunings alike; both exit 0 only when finished, 75 at the wall-time limit, and
# resubmitting the SAME command continues with the next segment): every stage is a group of attempts (--fit-attempts,
# default 3; --ft-attempts, default 4). Attempt a > 1 depends on afternotok of ALL earlier attempts of the stage (it
# runs only when every earlier attempt ended non-zero; when one succeeded Slurm cancels it, no node is allocated); the
# next stage depends on afterok of ANY attempt ("afterok:A1?afterok:A2?..."). All jobs: --kill-on-invalid-dep=yes, so
# when every attempt of a stage failed, the rest of the chain is cancelled. Fit jobs run with CHAIN=1: an attempt that
# finds a segment directory without restart file (a real failure, not 75) exits 3 instead of rerunning it. A fine-tuning
# attempt after a real failure (exit 1-4) is handled by positive_finetune.sbatch itself (no CHAIN mode there).
# (Dependency pattern tested by the lead with dummy jobs 122455-122466 on this Slurm 22.05.9.)
#
# --start-at continues a chain (earlier stages done or submitted): e.g. the c1 fits were submitted by hand ->
#   --start-at c1ft --after <c1 fit job id>    (afterok; if that job stops with 75, resubmit it and use the new id), or,
#   once the c1 fit is complete (<fit>/final exists): --start-at c1ft (no dependency).
# Node / GPU choice: --fit-node / --ft-node (sbatch --nodelist), --fit-gpus / --ft-gpus (--gpus-per-node), any other
# sbatch options with --fit-sbatch / --ft-sbatch (both job scripts request --exclusive, 64 CPUs, --mem=0).
# Fine-tuning jobs always get --exclude=gpu003 and --time=<--ft-time> (required; the script's 1 h default is too short).
# --model-ext: MODEL_EXT of the fits (a frozen snapshot, model_ext_frozen/<name>; default: positive_fit.sbatch's
# default, the live model_ext). --ft-model-ext: MODEL_EXT of the fine tunings (default: --model-ext); it must contain
# the fine tuning's files (positive_finetune.py, test_positive_finetune.py, positive_problem.py, positive_activity.py,
# grain_crossing.py, kernels/{finetune,pressure,activity}.cu). The fine tuning gets FIT, OUT, MAXNFEV (--ft-maxnfev,
# required: the fine tuning's own setting, no default here), MODEL_EXT and any --ft-export variables (LSMR_MAXITER,
# VALIDATE, CHECKPOINT_MIN, FLOOR_RATIO).
# --c1-start-set "NAME=VALUE[:NAME=VALUE...]" (2026-11-08; user "kappa = 1 where stuck", e.g. "kappa=1"): passed as
# START_SET to the cycle-1 fit ONLY (all its attempts; positive_fit.sbatch -> positive_fit.py --start-set after
# --start-restart = the cycle-1 source); refused unless --start-at c1fit; commas are turned into colons (sbatch --export
# splits commas). Cycles 2 and 3 start from the previous fit's restart as usual. Absent: nothing changes.
set -uo pipefail
ROOT_OUT=/data001/highresolu_vt/cops_output
J=$ROOT_OUT/jobs; RUN=$ROOT_OUT/runs/shap8_012m_lategrid_dz1.5; LOG=$ROOT_OUT/model_ext/logs/submissions.txt
PY="/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python -W ignore"
LAMBDA=""; RATIO=0.002; SOURCE=""; CYCLES=3; START_AT=c1fit; AFTER=""; FIT_ATTEMPTS=3; FT_ATTEMPTS=4; NFEV=4000; CHECK=1
FIT_NODE=""; FT_NODE=""; FIT_GPUS=""; FT_GPUS=""; FIT_OPTS=""; FT_OPTS=""; MODEL_EXT=""; FT_MODEL_EXT=""; SUBMIT=0
FT_MAXNFEV=""; FT_TIME=""; FT_EXPORT=""; C1_START_SET=""
FT_FILES=(positive_finetune.py test_positive_finetune.py positive_problem.py positive_activity.py grain_crossing.py
          kernels/finetune.cu kernels/pressure.cu kernels/activity.cu)
die() { echo "positive_cycles.sh: $*" >&2; exit 2; }
while [ $# -gt 0 ]; do
    case $1 in
        --lambda) LAMBDA=$2; shift 2;;               --ratio) RATIO=$2; shift 2;;
        --source) SOURCE=$2; shift 2;;               --cycles) CYCLES=$2; shift 2;;
        --start-at) START_AT=$2; shift 2;;           --after) AFTER=$2; shift 2;;
        --fit-attempts|--attempts) FIT_ATTEMPTS=$2; shift 2;;
        --ft-attempts) FT_ATTEMPTS=$2; shift 2;;     --nfev) NFEV=$2; shift 2;;
        --check) CHECK=$2; shift 2;;
        --model-ext) MODEL_EXT=$2; shift 2;;         --ft-model-ext) FT_MODEL_EXT=$2; shift 2;;
        --fit-node) FIT_NODE=$2; shift 2;;           --ft-node) FT_NODE=$2; shift 2;;
        --fit-gpus) FIT_GPUS=$2; shift 2;;           --ft-gpus) FT_GPUS=$2; shift 2;;
        --fit-sbatch) FIT_OPTS=$2; shift 2;;         --ft-sbatch) FT_OPTS=$2; shift 2;;
        --ft-maxnfev) FT_MAXNFEV=$2; shift 2;;       --ft-time) FT_TIME=$2; shift 2;;
        --ft-export) FT_EXPORT=$2; shift 2;;
        --c1-start-set) C1_START_SET=$2; shift 2;;
        --submit) SUBMIT=1; shift;;
        -h|--help) sed -n '2,61p' "$0"; exit 0;;
        *) die "unknown argument $1 (--help)";;
    esac
done
[ -n "$LAMBDA" ] || die "--lambda is required (verbatim, e.g. 1e-2)"
case $RATIO in 0.02|0.002|1e-5) ;; *) die "--ratio $RATIO: the user decided 0.02, 0.002 and 1e-5";; esac
case $CYCLES in [1-9]) ;; *) die "--cycles must be 1..9";; esac
case $FIT_ATTEMPTS in [1-9]) ;; *) die "--fit-attempts must be 1..9";; esac
case $FT_ATTEMPTS in [1-9]) ;; *) die "--ft-attempts must be 1..9";; esac
case $START_AT in c[1-9]fit|c[1-9]ft) ;; *) die "--start-at c<k>fit or c<k>ft";; esac
[ "${START_AT:1:1}" -le "$CYCLES" ] || die "--start-at $START_AT is beyond --cycles $CYCLES"
if [ -n "$C1_START_SET" ]; then
    [ "$START_AT" = c1fit ] || die "--c1-start-set applies to the cycle-1 fit only: it needs --start-at c1fit"
    C1_START_SET=${C1_START_SET//,/:}
    case $C1_START_SET in *[!A-Za-z0-9_=.:+-]*|"") die "--c1-start-set NAME=VALUE[:NAME=VALUE...]";; esac
fi
case $FT_TIME in "") ;; [0-9]*:[0-9][0-9]:[0-9][0-9]|[0-9]*-[0-9]*:[0-9][0-9]:[0-9][0-9]) ;; *) die "--ft-time HH:MM:SS";; esac
case $FT_EXPORT in *ALL*|*FIT=*|*OUT=*|*MAXNFEV=*|*MODEL_EXT=*)
    die "--ft-export must not set ALL, FIT, OUT, MAXNFEV or MODEL_EXT";; esac
[ -z "$MODEL_EXT" ] || MODEL_EXT=$(readlink -f "$MODEL_EXT")
[ -z "$MODEL_EXT" ] || [ -f "$MODEL_EXT/positive_fit.py" ] || die "--model-ext $MODEL_EXT has no positive_fit.py"
FT_MODEL_EXT=${FT_MODEL_EXT:-$MODEL_EXT}
[ -z "$FT_MODEL_EXT" ] || FT_MODEL_EXT=$(readlink -f "$FT_MODEL_EXT")
for f in positive_fit.sbatch positive_finetune.sbatch; do [ -f "$J/$f" ] || die "missing $J/$f"; done
for f in "${FT_FILES[@]}"; do
    if [ ! -e "${FT_MODEL_EXT:-$ROOT_OUT/model_ext}/$f" ]; then
        [ "$SUBMIT" = 1 ] && die "the fine tuning's MODEL_EXT ${FT_MODEL_EXT:-$ROOT_OUT/model_ext} lacks $f (give --ft-model-ext)"
        echo "# NOTE: the fine tuning's MODEL_EXT ${FT_MODEL_EXT:-$ROOT_OUT/model_ext} lacks $f: printed as <ft-model-ext>"
        FT_MODEL_EXT="<ft-model-ext>"; break
    fi
done
if [ -z "$FT_MAXNFEV" ]; then
    [ "$SUBMIT" = 1 ] && die "--ft-maxnfev is required (positive_finetune.sbatch MAXNFEV; no default here)"
    echo "# NOTE: --ft-maxnfev not given: printed as <ft-maxnfev>"
fi
if [ -z "$FT_TIME" ]; then
    [ "$SUBMIT" = 1 ] && die "--ft-time is required (the fine tuning's time limit; its script default 1 h is too short)"
    echo "# NOTE: --ft-time not given: printed as <ft-time>"
fi
L=$RUN/positive_r$RATIO/lam$LAMBDA
# ---- cycle-1 source: the first scan's last segment (status complete)
if [ -z "$SOURCE" ]; then
    FIRST=$RUN/positive/C30_n4000_exp_cliff_w10x5_grain_cross10_lam$LAMBDA
    SEGLAST=$(ls -d "$FIRST"/seg[0-9][0-9] 2>/dev/null | tail -1)
    [ -n "$SEGLAST" ] || die "no first-scan fit $FIRST/seg<NN> for lambda $LAMBDA (give --source)"
    SOURCE=$SEGLAST/fit_restart.npz
fi
[ -f "$SOURCE" ] || [ "$START_AT" != c1fit ] || die "missing cycle-1 source $SOURCE"
if [ -f "$SOURCE" ] && [ "$START_AT" = c1fit ]; then
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
    $PY -c "import sys; sys.exit(0 if float(sys.argv[1]) == float(sys.argv[2]) else 1)" "$SLAMBDA" "$LAMBDA" ||
        die "the cycle-1 source has lambda $SLAMBDA, not $LAMBDA"
fi
# ---- what --start-at assumes is done (no --after): the previous stage's product must exist
if [ -z "$AFTER" ] && [ "$START_AT" != c1fit ]; then
    k0=${START_AT:1:1}
    case $START_AT in
        c*ft) need=$L/c$k0/fit/final/fit_restart.npz;;
        c*fit) need=$L/c$((k0 - 1))/ft/finetune_state.npz;;
    esac
    if [ -e "$need" ]; then echo "# --start-at $START_AT: found $need -> $(readlink -f "$need")"
    elif [ "$SUBMIT" = 1 ]; then die "--start-at $START_AT without --after, but $need does not exist"
    else echo "# NOTE: --start-at $START_AT: $need does not exist yet"; fi
fi
# ---- submission (or its dry-run print)
submit() {   # $1 = label, rest = sbatch arguments; sets JOBID
    local label=$1; shift
    if [ "$SUBMIT" = 1 ]; then
        JOBID=$(sbatch --parsable "$@") || die "sbatch failed for $label"
        echo "$label: job $JOBID"
        local node="any"; for a in "$@"; do case $a in --nodelist=*) node=${a#--nodelist=};; esac; done
        echo "$(date -u +%FT%TZ) $JOBID ${!#} $node positive_cycles.sh lambda $LAMBDA ratio $RATIO $label: $(printf '%s ' "$@" | cut -c1-700)" >> "$LOG"
    else
        JOBID="<$label>"
        printf '%s\n  sbatch' "# $label"
        for a in "$@"; do case $a in *[?\<\>*]*) printf " '%s'" "$a";; *) printf ' %s' "$a";; esac; done; echo
    fi
}
opts() {   # $1 node, $2 gpus, $3 extra -> OPTS array
    OPTS=(--kill-on-invalid-dep=yes)
    [ -n "$1" ] && OPTS+=(--nodelist="$1")
    [ -n "$2" ] && OPTS+=(--gpus-per-node="$2")
    if [ -n "$3" ]; then read -r -a extra <<< "$3"; OPTS+=("${extra[@]}"); fi
}
stage() {   # $1 label prefix, $2 attempts, $3 export, $4 script; uses DEP; sets DEP to "afterok:A1?afterok:A2?..."
    local ids=() okdep="" a d
    for a in $(seq 1 "$2"); do
        if [ "$a" = 1 ]; then d=$DEP; else d="afternotok:$(IFS=:; echo "${ids[*]}")"; fi
        submit "$1.$a" "${OPTS[@]}" ${d:+--dependency=$d} --job-name="pos_l${LAMBDA}_${1}$a" --export="$3" "$4"
        ids+=("$JOBID"); okdep+="${okdep:+?}afterok:$JOBID"
    done
    DEP=$okdep
}
DEP=""; [ -n "$AFTER" ] && DEP="afterok:$AFTER"
started=0
echo "# chain for lambda $LAMBDA, ratio $RATIO: $L/c<k>/{fit,ft}; $CYCLES cycles from $START_AT; attempts: fit $FIT_ATTEMPTS, fine tuning $FT_ATTEMPTS; fit MODEL_EXT ${MODEL_EXT:-(live model_ext)}; fine-tuning MODEL_EXT ${FT_MODEL_EXT:-(live model_ext)}; $([ "$SUBMIT" = 1 ] && echo SUBMIT || echo 'DRY RUN')"
[ -n "$C1_START_SET" ] && echo "# cycle-1 fit START_SET=$C1_START_SET (positive_fit.py --start-set; cycles 2+ unchanged)"
for k in $(seq 1 "$CYCLES"); do
    for st in fit ft; do
        [ "$started" = 1 ] || [ "c$k$st" = "$START_AT" ] || continue
        started=1
        if [ "$st" = fit ]; then
            VARS="ALL,LAMBDA=$LAMBDA,RATIO=$RATIO,CYCLE=$k,NFEV=$NFEV,CHECK=$CHECK,CHAIN=1"
            [ -n "$MODEL_EXT" ] && VARS+=",MODEL_EXT=$MODEL_EXT"
            if [ "$k" = 1 ]; then
                VARS+=",START_RESTART=$SOURCE"
                [ -n "$C1_START_SET" ] && VARS+=",START_SET=$C1_START_SET"
            else
                VARS+=",START_RESTART=$L/c$((k - 1))/fit/final/fit_restart.npz"
                VARS+=",HELD=$L/c$((k - 1))/ft/finetune_state.npz,GRAIN_FROM=$L/c$((k - 1))/ft/finetune_state.npz"
            fi
            opts "$FIT_NODE" "$FIT_GPUS" "$FIT_OPTS"
            stage "c${k}fit" "$FIT_ATTEMPTS" "$VARS" "$J/positive_fit.sbatch"
        else
            VARS="ALL,FIT=$L/c$k/fit/final,OUT=$L/c$k/ft,MAXNFEV=${FT_MAXNFEV:-<ft-maxnfev>}"
            [ -n "$FT_MODEL_EXT" ] && VARS+=",MODEL_EXT=$FT_MODEL_EXT"
            [ -n "$FT_EXPORT" ] && VARS+=",$FT_EXPORT"
            opts "$FT_NODE" "$FT_GPUS" "$FT_OPTS"
            OPTS+=(--exclude=gpu003 --time="${FT_TIME:-<ft-time>}")
            stage "c${k}ft" "$FT_ATTEMPTS" "$VARS" "$J/positive_finetune.sbatch"
        fi
    done
done
[ "$SUBMIT" = 1 ] || echo "# dry run: nothing submitted (add --submit)"
