#!/bin/bash
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/nonneg_alt_stage1.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/nonneg_alt_stage1.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 77057755e84350ca8c3a928234adb7cfbbc4785cdf296b635d4a5dfc6fd3e29b
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: launcher. Stage-1 fits of options A, B, D (user 2026-11-14): one 8-GPU job per fit, 3 attempts
#   (afternotok), guard, summary job; run roots nonneg/<model>/stage1_<setting>_<UTC>/ with README.md,
#   HOW-TO-RERUN.md.
# RUN: nonneg_alt_stage1.sh --model-ext <snapshot> --utc <UTC> --runs bg:1e-3,...,sos:3,ls:0.1 --spread-cref C
#   [--submit] (2026-11-14: B x3 + A jobs 125032-125051; D see the README).
# 2026-11-14: run copy updated (summary snapshot nonneg_summary_20261114T090853Z)
# 2026-11-14: run copy updated (summary snapshot nonneg_summary_20261114T090853Z; sos start text from the start file, --note; template sign text)
# ==== end of header ====
# Stage 1 of the alternative non-negative activity models (options A, B, D; user 2026-11-14 "implement the
# recommendations"): single cycle-1 fits, one 8-GPU job each, through jobs/nonneg_fit.sbatch.
# Common scheme: lategrid base, D2's setup, ft_c1..ft_c6 held at 1, slope-bin multipliers free >= 1 (start 1),
# max_nfev 4000, the full Jacobian check in segment 1, 32 host threads.
#
#   nonneg_alt_stage1.sh --model-ext <frozen snapshot> --utc <UTC> --runs SPEC[,SPEC...] [--submit]
# SPEC (one fit each, submitted in the given order):
#   bg:<V>     option B, --model background --bg-min V, start d2 (A_bg = V, its lower bound)
#   sos:<r>    option A, --model sos --sos-rank r, start gram (--sos-start-file, default the Gram start of job 124994)
#   ls:<f>     option D, --model logspread --spread-f f --spread-cref <C_ref>, start uniform (needs --spread-cref)
# Run root: runs/shap8_012m_lategrid_dz1.5/nonneg/<model>/stage1_<setting>_<UTC>/ (setting bg<V>, r<r>, f<f>) with the
# fit in <root>/fit (fit/final when complete), README.md and HOW-TO-RERUN.md.
# Jobs: --attempts N (default 3) per fit: attempt 1 has no dependency, attempt a > 1 afternotok of ALL earlier attempts
# (it runs only after a time-limit stop or failure of every earlier one; nonneg_fit.sbatch continues from the last
# segment), every job --kill-on-invalid-dep=yes; then a guard job (jobs/positive_chain_guard.sbatch KIND=fit, CPU,
# seconds) with afterok:A1?...?afterok:AN that exits 0 only when <root>/fit/final is complete (the record that the
# fit finished). No dependency between fits: they start in submission order as GPU nodes free (gpu001, gpu002; the
# fit sbatch excludes gpu003). --node NODE pins them.
# Options: --sos-start-file NPZ, --spread-cref C, --time T (default 08:00:00, the gpu maximum), --attempts N,
#   --node NODE, --guard-partition P (normal128cores).
set -uo pipefail
COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; J=$COPS_OUTPUT_ROOT/jobs
SUBLOG=$COPS_OUTPUT_ROOT/model_ext/logs/submissions.txt
AMPS=/data001/highresolu_vt/AMPS
MX=""; UTC=""; RUNS=""; SOSFILE=$R/nonneg/bundle/sos_start_r3.npz; CREF=""; TIME=08:00:00; ATTEMPTS=3; NODE=""
SX=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_summary_20261114T090853Z
GP=normal128cores; SUBMIT=0; NOTE=""
while [ $# -gt 0 ]; do
    case $1 in
        --model-ext) MX=$(readlink -f "$2"); shift 2;;
        --utc) UTC=$2; shift 2;;
        --runs) RUNS=$2; shift 2;;
        --sos-start-file) SOSFILE=$(readlink -f "$2"); shift 2;;
        --spread-cref) CREF=$2; shift 2;;
        --summary-ext) SX=$(readlink -f "$2"); shift 2;;
        --note) NOTE=$2; shift 2;;
        --time) TIME=$2; shift 2;;
        --attempts) ATTEMPTS=$2; shift 2;;
        --node) NODE=$2; shift 2;;
        --guard-partition) GP=$2; shift 2;;
        --submit) SUBMIT=1; shift;;
        *) echo "unknown option $1" >&2; exit 2;;
    esac
done
die() { echo "nonneg_alt_stage1: $*" >&2; exit 2; }
[ -n "$MX" ] && [ -f "$MX/nonneg_alt.py" ] && [ -f "$MX/SHA256SUMS" ] || die "--model-ext <frozen alternative-model snapshot>"
case $UTC in [0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]T[0-9][0-9][0-9][0-9][0-9][0-9]Z) ;; *) die "--utc YYYYMMDDTHHMMSSZ";; esac
[ -n "$RUNS" ] || die "--runs SPEC[,SPEC...]"
case $ATTEMPTS in [1-9]) ;; *) die "--attempts 1..9";; esac
[ -f "$J/positive_chain_guard.sbatch" ] || die "missing $J/positive_chain_guard.sbatch"
[ -f "$SX/nonneg_alt_summary.py" ] || die "--summary-ext: no nonneg_alt_summary.py in $SX"
NOPT=(); [ -n "$NODE" ] && NOPT=(--nodelist="$NODE")
COMMON="SLOPE_MULT_MIN=1,NFEV=4000,CHECK=1,L_THREADS=32,MODEL_EXT=$MX"
GIT=$(git -C "$AMPS" rev-parse HEAD 2>/dev/null)
echo "# stage 1: $RUNS; MODEL_EXT $MX; time $TIME; $ATTEMPTS attempts + guard ($GP); ${NODE:-any gpu node but gpu003}; $([ "$SUBMIT" = 1 ] && echo SUBMIT || echo 'DRY RUN')"
IFS=, read -r -a SPECS <<< "$RUNS"
for spec in "${SPECS[@]}"; do
    kind=${spec%%:*}; val=${spec#*:}
    [ -n "$val" ] && [ "$val" != "$spec" ] || die "SPEC $spec: kind:value"
    case $kind in
        bg) MODEL=background; SET=bg$val; EXP="MODEL=background,BG_MIN=$val,START=d2"
            DESC="option B (background): A = m_b (A_bg A_orig + s ln(1 + exp(L/s))) exp(h), A_bg >= $val (start $val), start d2"
            START_TXT="D2's harmonics a00..a24, kappa and grain by name (as the softplus), multipliers 1, A_bg = $val (its lower bound)";;
        sos) MODEL=sos; SET=r$val; EXP="MODEL=sos,SOS_RANK=$val,START=gram,START_FILE=$SOSFILE"
            [ -f "$SOSFILE" ] || die "missing sos start file $SOSFILE"
            read -r SREF SJOB < <(/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python -c "
import json, sys, numpy as np, os
m = json.loads(str(np.load(sys.argv[1])['meta_json']))
print(m.get('refine', 'none'), os.path.basename(os.path.dirname(m['tool']['path'])))" "$SOSFILE")
            DESC="option A (sum of squares): A = m_b (||L^T y||^2 / A_s + A_s sum T) exp(h), L 25 x $val lower-trapezoidal, start gram ($( [ "$SREF" = lm ] && echo "LM-refined rank-$val LSQ" || echo "Gram/PSD"))"
            if [ "$SREF" = lm ]; then
                START_TXT="L from $SOSFILE (nonneg_start.py --model sos --rank $val --refine lm, $SJOB: the rank-$val least-squares fit of D2's activity by Levenberg-Marquardt, from the Gram/PSD start), kappa and grain by name from D2, multipliers 1"
            else
                START_TXT="L from $SOSFILE (nonneg_start.py --model sos --rank $val, $SJOB: area-weighted LSQ Gram matrix of D2's activity, PSD projection, top-$val eigenpairs), kappa and grain by name from D2, multipliers 1"
            fi;;
        ls) MODEL=logspread; SET=f$val; [ -n "$CREF" ] || die "ls:<f> needs --spread-cref"
            EXP="MODEL=logspread,SPREAD_F=$val,SPREAD_CREF=$CREF,START=uniform"
            DESC="option D (log-normal spread): A = m_b (A_orig exp(u) + A_s sum T) exp(h), u = Y a / A_s, penalty f C_ref Var_w(u), f = $val, C_ref = $CREF, start uniform"
            START_TXT="uniform: u = ln s (a00 only, the other harmonics 0), s = the least-squares fit in ln p of the model pressures to D2's modelled pressures over the fit's used NG and RG rows (recorded in the manifest, model_check.alternative.uniform_start); kappa and grain by name from D2, multipliers 1";;
        *) die "SPEC kind $kind: bg, sos or ls";;
    esac
    ROOT=$R/nonneg/$MODEL/stage1_${SET}_$UTC; FD=$ROOT/fit
    EXPORT="ALL,$EXP,$COMMON,FITDIR=$FD"
    if [ "$SUBMIT" != 1 ]; then
        echo "## $spec -> $ROOT"
        for a in $(seq 1 "$ATTEMPTS"); do
            d=""; [ "$a" -gt 1 ] && d="--dependency=afternotok:<attempts 1..$((a - 1))> "
            echo "sbatch --parsable ${NOPT[*]} --kill-on-invalid-dep=yes ${d}--time=$TIME --job-name=nonneg_s1_${SET} --export=$EXPORT $J/nonneg_fit.sbatch"
        done
        echo "sbatch --parsable --kill-on-invalid-dep=yes --partition=$GP --dependency=afterok:<A1>?...?afterok:<A$ATTEMPTS> --job-name=nonneg_s1_${SET}_guard --export=ALL,KIND=fit,TARGET=$FD,NFEV=4000 $J/positive_chain_guard.sbatch"
        echo "sbatch --parsable --kill-on-invalid-dep=yes --dependency=afterok:<guard> --export=ALL,FIT=$FD/final,OUT=$ROOT/summary.txt,MODEL_EXT=$SX $J/nonneg_alt_summary.sbatch"
        continue
    fi
    [ -e "$ROOT" ] && die "exists: $ROOT (a new stage-1 run needs a new --utc; resume an existing one by resubmitting its HOW-TO-RERUN.md command)"
    mkdir -p "$ROOT"
    ids=(); okdep=""; CMDS=()
    for a in $(seq 1 "$ATTEMPTS"); do
        args=(--parsable "${NOPT[@]}" --kill-on-invalid-dep=yes --time="$TIME" --job-name="nonneg_s1_${SET}" --export="$EXPORT")
        [ "$a" -gt 1 ] && args+=(--dependency="afternotok:$(IFS=:; echo "${ids[*]}")")
        id=$(sbatch "${args[@]}" "$J/nonneg_fit.sbatch") || die "sbatch failed: $spec attempt $a"
        ids+=("$id"); okdep+="${okdep:+?}afterok:$id"
        CMDS+=("sbatch ${args[*]} $J/nonneg_fit.sbatch")
        echo "$(date -u +%FT%TZ) $id jobs/nonneg_fit.sbatch ${NODE:-gpu (any but gpu003)} stage 1 $spec attempt $a/$ATTEMPTS (nonneg_alt_stage1.sh): $DESC; $EXP,$COMMON; FITDIR=$FD; time $TIME${a:+}$([ "$a" -gt 1 ] && echo "; afternotok:$(IFS=:; echo "${ids[*]:0:$((a - 1))}")")" >> "$SUBLOG"
    done
    gargs=(--parsable --kill-on-invalid-dep=yes --partition="$GP" --dependency="$okdep" --job-name="nonneg_s1_${SET}_guard"
           --export="ALL,KIND=fit,TARGET=$FD,NFEV=4000")
    gid=$(sbatch "${gargs[@]}" "$J/positive_chain_guard.sbatch") || die "sbatch failed: $spec guard"
    echo "$(date -u +%FT%TZ) $gid jobs/positive_chain_guard.sbatch $GP stage 1 $spec guard (KIND=fit TARGET=$FD NFEV=4000): $okdep" >> "$SUBLOG"
    sargs=(--parsable --kill-on-invalid-dep=yes --dependency="afterok:$gid" --export="ALL,FIT=$FD/final,OUT=$ROOT/summary.txt,MODEL_EXT=$SX")
    sid=$(sbatch "${sargs[@]}" "$J/nonneg_alt_summary.sbatch") || die "sbatch failed: $spec summary"
    echo "$(date -u +%FT%TZ) $sid jobs/nonneg_alt_summary.sbatch normal128cores stage 1 $spec summary (FIT=$FD/final OUT=$ROOT/summary.txt MODEL_EXT=$SX): afterok:$gid" >> "$SUBLOG"
    now=$(date -u +%FT%TZ)
    cat > "$ROOT/README.md" <<EOF
# Stage 1, $DESC

Created $now by \`jobs/nonneg_alt_stage1.sh\` (fit attempts ${ids[*]}, guard $gid, summary $sid).

- **Why** (user 2026-11-14): "prepare suggestions for modification of the minimization polynomial: the current approach
  shuts down almost all the surface, and I need a polynomial that is non-negative and would allow fit with more seen
  surface contributing"; "implement the recommendations" (2026-11-14 ~02:35 CST, relayed by the lead).
- **Model and options:** \`srcCG/LiouvilleCOPS_cuda/model_ext/README.md\`, section "Alternative non-negative activity
  models (options A-D)"; this run: \`$EXP\`.
- **Common scheme:** lategrid base (gauge +1.5 m), D2's setup (W10x5 from 10:33:10 UTC, grain at the time stamp,
  NG = RG crossing weight 10, original temperature law held: Tmax 182.09922, T75 136.11149, Tmin 133 K; kappa free;
  lambda 0; --tbounds 30 400), D2's ft_c1..ft_c6 held at 1, slope-bin multipliers [0,40) held at 1 and six free >= 1
  (start 1, --slope-mult-min 1), max_nfev 4000.
- **Templates (lead decision 1, 2026-11-14; check \`checks/ft_template_sign_20261114/\`, job 125029):** D2's
  ft_c1..ft_c6 are each > 0 on the same 265 facets and 0 on the other 12,364,025; negative on none. A_s sum T_k has
  min 0 and max 7.894e23 (3.55e5 A_s), 350.9 m^2 in total; all 265 facets are seen; latitude 0.5 to 69.0 deg,
  longitude -11.6 to 166.7 deg. Since the sum is non-negative everywhere it is added as is: for A and D outside the
  non-negative map (times m_b exp(h)); for B inside the softplus argument, as in the softplus model.
- **Start:** $START_TXT.${NOTE:+
- **Note:** $NOTE}
- **Code:** frozen snapshot \`$MX\` (SHA256SUMS inside); AMPS $GIT.
- **Products:** \`fit/seg<NN>/\`, \`fit/final\` -> the last segment when complete (fit_result.txt, COPS_model.dat,
  late_window_model.dat, nonneg_facets.npz, nonneg_bins.txt, NONNEG_MODEL.txt, fit_manifest.json, fit_restart.npz,
  evaluations.txt), \`fit/provenance-<job>/\` (commands.sh with every argument resolved, git, nvidia-smi, pip list).
- **Summary:** \`summary.txt\` + \`summary.txt.json\` (job $sid after the guard: costs incl. the penalties, termination,
  rms ln NG/RG overall / before the gap / late window, parameters at bounds, Q_H2O at 10:33:10, active seen area,
  seen area holding 50 % / 90 % of the seen production, NG/RG shares from the active area; definitions in
  \`$SX/nonneg_alt_summary.py\`).
EOF
    {
        echo "# How to rerun: stage 1, $spec"
        echo
        echo "$DESC. Model: \`srcCG/LiouvilleCOPS_cuda/model_ext/README.md\`, section \"Alternative non-negative"
        echo "activity models (options A-D)\". Run root: \`$ROOT\`."
        echo
        echo "## Inputs"
        echo "- base \`$R/base\`, points \`$R/late_grid_points.txt\`, slope bins \`$R/slope_mult/bundle/slope_bins.npz\`;"
        echo "- D2 \`$R/clip_refit/D2/final\` (fit_restart.npz: kappa, grain$([ "$MODEL" = background ] && echo ", harmonics"); fit_manifest.json: the ft templates);"
        [ "$MODEL" = sos ] && echo "- sos start \`$SOSFILE\` (refine $SREF; its provenance and commands.sh next to it, $SJOB);"
        echo "- code: \`$MX\` (frozen; SHA256SUMS); AMPS $GIT."
        echo
        echo "## Commands (submitted $now)"
        echo '```bash'
        echo "export COPS_OUTPUT_ROOT=$COPS_OUTPUT_ROOT"
        for k in "${!CMDS[@]}"; do
            echo "# attempt $((k + 1)) = job ${ids[$k]}"
            echo "${CMDS[$k]}"
        done
        echo "# guard = job $gid"
        echo "sbatch ${gargs[*]} $J/positive_chain_guard.sbatch"
        echo '```'
        echo "To rerun from scratch, use a new FITDIR (or a new --utc with \`jobs/nonneg_alt_stage1.sh --model-ext $MX --utc <UTC> --runs $spec$([ "$MODEL" = logspread ] && echo " --spread-cref $CREF") --submit\`)."
        echo "The same nonneg_fit.sbatch command on the existing FITDIR continues the fit from its last segment (after a"
        echo "time-limit stop). Each segment's \`fit/provenance-<job>/commands.sh\` reruns that segment alone with every"
        echo "argument resolved."
        echo
        echo "## Summary"
        echo '```bash'
        echo "# job $sid (after the guard)"
        echo "sbatch ${sargs[*]} $J/nonneg_alt_summary.sbatch"
        echo '```'
    } > "$ROOT/HOW-TO-RERUN.md"
    echo "$spec: attempts ${ids[*]}, guard $gid, summary $sid -> $ROOT"
done
