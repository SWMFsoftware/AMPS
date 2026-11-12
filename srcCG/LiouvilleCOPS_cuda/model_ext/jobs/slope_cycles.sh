#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/slope_cycles.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/slope_cycles.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: f9319b9795869a8df784453afa972b6f9c87c2ca24330c60363c0a403affae20
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job chain. Fit / fine-tuning cycles of the bump + slope-multiplier model (option 3a; user decisions 2026-11-10
#   "implement option 3a with slope multipliers", "add fine tuning -- fitting cycles as in other fits", "execute the fit
#   using the optimal number of threads"): 3 cycles with sbatch dependencies, from 2026-11-10 21:50 UTC a guard job per
#   stage (attempt 1 afterok:G, attempt a > 1 afterok:G,afternotok:A1..); dry run by default.
# RUN: slope_cycles.sh --model-ext <v8> --chain <chain dir> --start-at c1ft --after 124063 --l-threads 32 --cycles 3
#   --nfev 4000 --ft-maxnfev 500 --node gpu002 --submit   (124064-124086; the c1 fit 124008 was submitted by hand from
#   its dry run)
# ==== end of header ====
# Fit / fine-tuning cycles of the bump + slope-multiplier model (option 3a; user decisions 2026-11-10: "implement option
# 3a with slope multipliers", "add fine tuning -- fitting cycles as in other fits", "execute the fit using the optimal
# number of threads"): 3 cycles, the positive cycles' settings (jobs/positive_cycles.sh):
#
#   c1 fit -> c1 ft -> c2 fit -> c2 ft -> c3 fit -> c3 ft     (sbatch --dependency, --kill-on-invalid-dep=yes)
#
#   c1 fit:  jobs/slope_mult_fit.sbatch BUMPS=<bump bundle> START=nnls START_RESTART=<D2 final restart> (bump weights =
#            the NNLS start, multipliers 1, kappa 1, D2's grain), NFEV 4000, CHECK 1, L_THREADS (chunked, bitwise)
#   ck ft:   jobs/slope_finetune.sbatch FIT=<ck>/fit/final FACETS=late_only WINDOW_MIN=5 MAXNFEV=500 (A_e -> max(F, m L
#            exp(h + delta)), delta >= ln(F/s), grain refitted), L_THREADS
#   ck+1 fit: START=chain START_RESTART=<ck>/fit/final/fit_restart.npz (weights, multipliers, kappa by name),
#            HELD=GRAIN_FROM=<ck>/ft/finetune_state.npz (held_next = accumulated delta; grain A, B of the fine tuning)
# Products: <chain>/c<k>/fit/seg<NN> (+ final), <chain>/c<k>/ft/seg<NN> (+ final, finetune_state.npz, ...).
# Default: DRY RUN (prints every sbatch command, job ids as <placeholders>); --submit submits and logs each job in
# model_ext/logs/submissions.txt.
#
#   jobs/slope_cycles.sh --model-ext DIR [--chain DIR] [--bumps NPZ] [--source R.npz] [--cycles 3] [--nfev 4000]
#       [--l-threads 32] [--ft-maxnfev 500] [--fit-time 06:00:00] [--ft-time 02:00:00] [--node gpu002]
#       [--fit-attempts 3] [--ft-attempts 4] [--start-at c1fit|c1ft|...] [--after JOBID[,JOBID...]]
#       [--guard-partition normal128cores] [--submit]
#   --after: a gate for the first submitted stage. One id: that job is the gate (e.g. a GPU test job: the stage runs only
#   if it exits 0). A list: a guard job waits for afterok of ANY of them, checks the completion marker of the stage
#   before --start-at (c<k>ft: c<k>/fit; c<k>fit, k > 1: c<k-1>/ft; c1fit: none) and becomes the gate (e.g. the
#   attempts of the preceding stage submitted elsewhere).
#
# Retries and guards (2026-11-10, after the flaw 124014 -> 124015: an attempt cancelled by --kill-on-invalid-dep counts
# as "not ok", so an attempt a > 1 that depended only on afternotok of the earlier attempts started although the
# upstream had failed; layout tested with dummy CPU jobs 124035-124058, checks/chain_deps_20261110T214653Z, 24 of 24):
#   every stage is a group of attempts behind ONE gate job G (Slurm 22.05 cannot mix "?" and "," in a dependency):
#     attempt 1: afterok:G;  attempt a > 1: afterok:G,afternotok:A1:...:A(a-1)   (no G for a chain started without
#     --after: attempt 1 has no dependency, attempt a > 1 afternotok only)
#   after the attempts, a guard job (jobs/slope_chain_guard.sbatch, CPU partition, seconds) depends on
#   afterok:A1?...?afterok:An and exits 0 only if the stage's completion marker is present (fit: <c>/fit/final with a
#   complete fit_restart.npz; ft: <c>/ft/final and finetune_state.npz passing final_link's product check); it is the
#   next stage's G. Every job has --kill-on-invalid-dep=yes: when every attempt of a stage failed, its guard is
#   cancelled and so is every later job (none starts). The last stage has a guard as well (a record that the chain
#   completed). Fit jobs run with CHAIN=1 (a segment without restart file stops the chain: exit 3). Both job scripts exit
#   75 at the wall-time limit and continue with the next segment when resubmitted (the next attempt does that).
set -uo pipefail
ROOT_OUT=/data001/highresolu_vt/cops_output
J=$ROOT_OUT/jobs; RUN=$ROOT_OUT/runs/shap8_012m_lategrid_dz1.5; LOG=$ROOT_OUT/model_ext/logs/submissions.txt
MODEL_EXT=""; CHAINDIR=""; BUMPS=$RUN/slope_bumps/bundle/bumps32.npz; SOURCE=$RUN/clip_refit/D2/final/fit_restart.npz
CYCLES=3; NFEV=4000; LT=32; FT_MAXNFEV=500; FIT_TIME=06:00:00; FT_TIME=02:00:00; NODE=gpu002
FIT_ATTEMPTS=3; FT_ATTEMPTS=4; START_AT=c1fit; AFTER=""; SUBMIT=0; GUARD_PARTITION=normal128cores
die() { echo "slope_cycles.sh: $*" >&2; exit 2; }
while [ $# -gt 0 ]; do
    case $1 in
        --model-ext) MODEL_EXT=$2; shift 2;;   --chain) CHAINDIR=$2; shift 2;;   --bumps) BUMPS=$2; shift 2;;
        --source) SOURCE=$2; shift 2;;         --cycles) CYCLES=$2; shift 2;;    --nfev) NFEV=$2; shift 2;;
        --l-threads) LT=$2; shift 2;;          --ft-maxnfev) FT_MAXNFEV=$2; shift 2;;
        --fit-time) FIT_TIME=$2; shift 2;;     --ft-time) FT_TIME=$2; shift 2;;  --node) NODE=$2; shift 2;;
        --fit-attempts) FIT_ATTEMPTS=$2; shift 2;; --ft-attempts) FT_ATTEMPTS=$2; shift 2;;
        --start-at) START_AT=$2; shift 2;;     --after) AFTER=$2; shift 2;;      --submit) SUBMIT=1; shift;;
        --guard-partition) GUARD_PARTITION=$2; shift 2;;
        -h|--help) sed -n '2,39p' "$0"; exit 0;;
        *) die "unknown argument $1 (--help)";;
    esac
done
[ -n "$MODEL_EXT" ] || die "--model-ext <frozen snapshot> is required"
MODEL_EXT=$(readlink -f "$MODEL_EXT")
for f in slope_mult_fit.py slope_finetune.py slope_mult.py slope_mult_problem.py positive_finetune.py; do
    [ -f "$MODEL_EXT/$f" ] || die "--model-ext $MODEL_EXT lacks $f"; done
[ -f "$BUMPS" ] || die "missing bump bundle $BUMPS"
[ -f "$SOURCE" ] || [ "$START_AT" != c1fit ] || die "missing cycle-1 source $SOURCE"
case $CYCLES in [1-9]) ;; *) die "--cycles 1..9";; esac
case $START_AT in c[1-9]fit|c[1-9]ft) ;; *) die "--start-at c<k>fit or c<k>ft";; esac
CHAINDIR=$(readlink -m "${CHAINDIR:-$RUN/slope_bumps/chain_$(basename "$MODEL_EXT")}")   # absolute: the jobs cd
BUMPS=$(readlink -f "$BUMPS"); [ -f "$SOURCE" ] && SOURCE=$(readlink -f "$SOURCE")
[ -f "$J/slope_chain_guard.sbatch" ] || die "missing $J/slope_chain_guard.sbatch"
submit() {   # $1 label, $2 where (node / partition, for the log), rest sbatch arguments; sets JOBID
    local label=$1 where=$2; shift 2
    if [ "$SUBMIT" = 1 ]; then
        JOBID=$(sbatch --parsable "$@") || die "sbatch failed for $label"
        echo "$label: job $JOBID"
        echo "$(date -u +%FT%TZ) $JOBID ${!#} $where slope_cycles.sh $label (bump + slope-multiplier chain $CHAINDIR): $(printf '%s ' "$@" | cut -c1-900)" >> "$LOG"
    else
        JOBID="<$label>"
        printf '%s\n  sbatch' "# $label"
        for a in "$@"; do case $a in *[?\<\>*]*) printf " '%s'" "$a";; *) printf ' %s' "$a";; esac; done; echo
    fi
}
guard() {   # $1 label, $2 any-of dependency ("afterok:A1?afterok:A2?..."), $3 KIND, $4 TARGET; sets GATE
    submit "$1" "$GUARD_PARTITION" --kill-on-invalid-dep=yes --partition="$GUARD_PARTITION" --dependency="$2" \
        --job-name="bumps_$1" --export="ALL,KIND=$3,TARGET=$4" "$J/slope_chain_guard.sbatch"
    GATE=$JOBID
}
stage() {   # $1 label, $2 attempts, $3 export, $4 time, $5 script, $6 guard KIND, $7 guard TARGET; uses and sets GATE
    local ids=() okdep="" a d
    for a in $(seq 1 "$2"); do
        if [ "$a" = 1 ]; then d=${GATE:+afterok:$GATE}
        else d="${GATE:+afterok:$GATE,}afternotok:$(IFS=:; echo "${ids[*]}")"; fi
        submit "$1.$a" "$NODE" --kill-on-invalid-dep=yes --nodelist="$NODE" --exclude=gpu003 --time="$4" \
            ${d:+"--dependency=$d"} --job-name="bumps_${1}$a" --export="$3" "$5"
        ids+=("$JOBID"); okdep+="${okdep:+?}afterok:$JOBID"
    done
    guard "$1.guard" "$okdep" "$6" "$7"
}
GATE=""
echo "# bump + slope-multiplier chain: $CHAINDIR/c<k>/{fit,ft}; $CYCLES cycles from $START_AT; MODEL_EXT $MODEL_EXT;" \
     "L_THREADS $LT (chunked, bitwise); node $NODE; guards on $GUARD_PARTITION;" \
     "$([ "$SUBMIT" = 1 ] && echo SUBMIT || echo 'DRY RUN')"
if [ -n "$AFTER" ]; then
    n=0; anydep=""
    for j in ${AFTER//,/ }; do
        case $j in *[!0-9]*|'') die "--after: job ids (digits), comma separated: $AFTER";; esac
        anydep+="${anydep:+?}afterok:$j"; n=$((n + 1)); GATE=$j
    done
    if [ "$n" -gt 1 ]; then          # the guard checks the marker of the stage before START_AT (none before c1fit)
        k0=${START_AT:1:1}
        case $START_AT in
            c1fit) guard "after.guard" "$anydep" gate "";;
            c*fit) guard "after.guard" "$anydep" ft "$CHAINDIR/c$((k0 - 1))/ft";;
            c*ft) guard "after.guard" "$anydep" fit "$CHAINDIR/c$k0/fit";;
        esac
    fi
fi
started=0
for k in $(seq 1 "$CYCLES"); do
    for st in fit ft; do
        [ "$started" = 1 ] || [ "c$k$st" = "$START_AT" ] || continue
        started=1
        if [ "$st" = fit ]; then
            V="ALL,BUMPS=$BUMPS,FITDIR=$CHAINDIR/c$k/fit,NFEV=$NFEV,CHECK=1,CHAIN=1,L_THREADS=$LT,MODEL_EXT=$MODEL_EXT"
            if [ "$k" = 1 ]; then V+=",START=nnls,START_RESTART=$SOURCE"
            else
                p=$CHAINDIR/c$((k - 1))
                V+=",START=chain,START_RESTART=$p/fit/final/fit_restart.npz,HELD=$p/ft/finetune_state.npz"
                V+=",GRAIN_FROM=$p/ft/finetune_state.npz"
            fi
            stage "c${k}fit" "$FIT_ATTEMPTS" "$V" "$FIT_TIME" "$J/slope_mult_fit.sbatch" fit "$CHAINDIR/c$k/fit"
        else
            V="ALL,FIT=$CHAINDIR/c$k/fit/final,OUT=$CHAINDIR/c$k/ft,FACETS=late_only,WINDOW_MIN=5"
            V+=",MAXNFEV=$FT_MAXNFEV,L_THREADS=$LT,MODEL_EXT=$MODEL_EXT"
            stage "c${k}ft" "$FT_ATTEMPTS" "$V" "$FT_TIME" "$J/slope_finetune.sbatch" ft "$CHAINDIR/c$k/ft"
        fi
    done
done
[ "$SUBMIT" = 1 ] || echo "# dry run: nothing submitted (add --submit)"
