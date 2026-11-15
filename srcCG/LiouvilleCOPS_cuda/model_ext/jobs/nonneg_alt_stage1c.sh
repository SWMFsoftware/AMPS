#!/bin/bash
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/nonneg_alt_stage1c.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/nonneg_alt_stage1c.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 8a4420c3d44a5b58cfc31442d9ac271f4b2431e4af76e1920a454afca3e93467
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: launcher. Continuations of stage-1 fits (user 2026-11-14 ~17:30 CST): --start restart from each source fit's
#   final state (optimizer scale rebuilt), +4000 evaluations, settings read from the source, --spread-u0 for D;
#   3 attempts, guard, summary; run roots nonneg/<model>/stage1c_<setting>_<UTC>/.
# RUN: nonneg_alt_stage1c.sh --model-ext <snapshot> --utc <UTC> --from ROOT[,...] [--spread-u0 U] [--submit].
# ==== end of header ====
# Continuations of stage-1 fits of options A, B, D (user 2026-11-14 ~17:30 CST, after the convergence diagnosis
# checks/stage1_diag_20261114/): a NEW fit from each source fit's final state (--start restart: every parameter by name,
# no held term), so that scipy's x_scale='jac' scale is rebuilt from the Jacobian there; NFEV more evaluations (4000).
# The model settings are read from the source fit (model, bg_min, sos_rank, spread_f, spread_cref); --spread-u0 U adds
# the fixed-level spread penalty to logspread continuations (user 2026-11-14: P = f C_ref [Var_w(u) + (u_bar - u0)^2]).
#
#   nonneg_alt_stage1c.sh --model-ext <frozen snapshot> --utc <UTC> --from ROOT[,ROOT...] [--spread-u0 U] [--submit]
# ROOT: a stage-1 run root runs/.../nonneg/<model>/stage1_<setting>_<UTC0> with a complete fit/final. The fits are
# submitted in the given order; new run roots runs/.../nonneg/<model>/stage1c_<setting>_<UTC>/ (fit/, README.md,
# HOW-TO-RERUN.md, summary.txt after the guard). Jobs: --attempts N (3) per fit (attempt a > 1 afternotok of all
# earlier), a guard (positive_chain_guard.sbatch KIND=fit, CPU) and the summary job (nonneg_alt_summary.sbatch), every
# job --kill-on-invalid-dep=yes; no dependency between fits (submission order; any gpu node but gpu003 unless --node).
# Options: --nfev N (4000), --time T (08:00:00), --attempts N, --node NODE, --summary-ext DIR, --guard-partition P.
set -uo pipefail
COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; J=$COPS_OUTPUT_ROOT/jobs
SUBLOG=$COPS_OUTPUT_ROOT/model_ext/logs/submissions.txt
AMPS=/data001/highresolu_vt/AMPS
PY=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python
MX=""; UTC=""; FROM=""; U0=""; NFEV=4000; TIME=08:00:00; ATTEMPTS=3; NODE=""; GP=normal128cores; SUBMIT=0
SX=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_summary_20261115T020314Z
while [ $# -gt 0 ]; do
    case $1 in
        --model-ext) MX=$(readlink -f "$2"); shift 2;;
        --utc) UTC=$2; shift 2;;
        --from) FROM=$2; shift 2;;
        --spread-u0) U0=$2; shift 2;;
        --nfev) NFEV=$2; shift 2;;
        --time) TIME=$2; shift 2;;
        --attempts) ATTEMPTS=$2; shift 2;;
        --node) NODE=$2; shift 2;;
        --summary-ext) SX=$(readlink -f "$2"); shift 2;;
        --guard-partition) GP=$2; shift 2;;
        --submit) SUBMIT=1; shift;;
        *) echo "unknown option $1" >&2; exit 2;;
    esac
done
die() { echo "nonneg_alt_stage1c: $*" >&2; exit 2; }
[ -n "$MX" ] && [ -f "$MX/nonneg_alt.py" ] && [ -f "$MX/SHA256SUMS" ] || die "--model-ext <frozen snapshot>"
grep -q '"restart"' "$MX/nonneg_problem.py" || die "$MX has no --start restart (an older snapshot)"
case $UTC in [0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]T[0-9][0-9][0-9][0-9][0-9][0-9]Z) ;; *) die "--utc YYYYMMDDTHHMMSSZ";; esac
[ -n "$FROM" ] || die "--from ROOT[,ROOT...]"
case $ATTEMPTS in [1-9]) ;; *) die "--attempts 1..9";; esac
[ -f "$SX/nonneg_alt_summary.py" ] || die "--summary-ext: no nonneg_alt_summary.py in $SX"
NOPT=(); [ -n "$NODE" ] && NOPT=(--nodelist="$NODE")
GIT=$(git -C "$AMPS" rev-parse HEAD 2>/dev/null)
echo "# stage 1 continuations: $FROM; MODEL_EXT $MX; NFEV $NFEV; time $TIME; $ATTEMPTS attempts + guard + summary; ${NODE:-any gpu node but gpu003}${U0:+; logspread spread u0 $U0}; $([ "$SUBMIT" = 1 ] && echo SUBMIT || echo 'DRY RUN')"
IFS=, read -r -a SRCS <<< "$FROM"
for src in "${SRCS[@]}"; do
    src=$(readlink -f "$src"); [ -d "$src" ] || die "no run root $src"
    RS=$(readlink -f "$src/fit/final/fit_restart.npz"); [ -f "$RS" ] || die "$src has no complete fit/final"
    read -r MODEL BG SR SF SC SU0 MSG NF COST < <($PY - "$src/fit/final" <<'PYEOF'
import json, sys
m = json.load(open(sys.argv[1] + "/fit_manifest.json")); c = m["model_extension"]["nonneg"]["config"]; o = m["optimizer_results"][0]
g = lambda k: "-" if c.get(k) is None else repr(c[k])
print(c["model"], g("bg_min"), g("sos_rank"), g("spread_f"), g("spread_cref"), g("spread_u0"),
      str(o.get("message")).replace(" ", "_")[:60], o.get("nfev"), repr(float(o["data_cost"])))
PYEOF
)
    base=$(basename "$src"); SET=${base#stage1_}; SET=${SET%_*}
    case $MODEL in
        background) EXP="MODEL=background,BG_MIN=$BG";;
        sos) EXP="MODEL=sos,SOS_RANK=$SR";;
        logspread) EXP="MODEL=logspread,SPREAD_F=$SF,SPREAD_CREF=$SC"; [ -n "$U0" ] && EXP+=",SPREAD_U0=$U0";;
        *) die "$src: model $MODEL is not an alternative model";;
    esac
    EXP+=",START=restart,START_RESTART=$RS,SLOPE_MULT_MIN=1,NFEV=$NFEV,CHECK=1,L_THREADS=32,MODEL_EXT=$MX"
    ROOT=$(dirname "$src")/stage1c_${SET}_$UTC; FD=$ROOT/fit
    EXPORT="ALL,$EXP,FITDIR=$FD"
    DESC="continuation of $(basename "$(dirname "$src")")/$base ($MODEL $SET; source data cost $COST, $NF evaluations, ${MSG//_/ })$( [ "$MODEL" = logspread ] && [ -n "$U0" ] && echo "; fixed-level spread penalty u0 = $U0")"
    if [ "$SUBMIT" != 1 ]; then
        echo "## $src -> $ROOT"
        echo "sbatch --parsable ${NOPT[*]} --kill-on-invalid-dep=yes --time=$TIME --job-name=nonneg_s1c_$SET --export=$EXPORT $J/nonneg_fit.sbatch   # attempt 1 of $ATTEMPTS; $DESC"
        continue
    fi
    [ -e "$ROOT" ] && die "exists: $ROOT"
    mkdir -p "$ROOT"
    ids=(); okdep=""; CMDS=()
    for a in $(seq 1 "$ATTEMPTS"); do
        args=(--parsable "${NOPT[@]}" --kill-on-invalid-dep=yes --time="$TIME" --job-name="nonneg_s1c_$SET" --export="$EXPORT")
        [ "$a" -gt 1 ] && args+=(--dependency="afternotok:$(IFS=:; echo "${ids[*]}")")
        id=$(sbatch "${args[@]}" "$J/nonneg_fit.sbatch") || die "sbatch failed: $src attempt $a"
        ids+=("$id"); okdep+="${okdep:+?}afterok:$id"; CMDS+=("sbatch ${args[*]} $J/nonneg_fit.sbatch")
        echo "$(date -u +%FT%TZ) $id jobs/nonneg_fit.sbatch ${NODE:-gpu (any but gpu003)} stage-1 continuation attempt $a/$ATTEMPTS (nonneg_alt_stage1c.sh): $DESC; $EXP; FITDIR=$FD; time $TIME$([ "$a" -gt 1 ] && echo "; afternotok:$(IFS=:; echo "${ids[*]:0:$((a - 1))}")")" >> "$SUBLOG"
    done
    gargs=(--parsable --kill-on-invalid-dep=yes --partition="$GP" --dependency="$okdep" --job-name="nonneg_s1c_${SET}_guard"
           --export="ALL,KIND=fit,TARGET=$FD,NFEV=$NFEV")
    gid=$(sbatch "${gargs[@]}" "$J/positive_chain_guard.sbatch") || die "sbatch failed: guard"
    echo "$(date -u +%FT%TZ) $gid jobs/positive_chain_guard.sbatch $GP continuation guard (KIND=fit TARGET=$FD NFEV=$NFEV): $okdep" >> "$SUBLOG"
    sargs=(--parsable --kill-on-invalid-dep=yes --dependency="afterok:$gid" --export="ALL,FIT=$FD/final,OUT=$ROOT/summary.txt,MODEL_EXT=$SX")
    sid=$(sbatch "${sargs[@]}" "$J/nonneg_alt_summary.sbatch") || die "sbatch failed: summary"
    echo "$(date -u +%FT%TZ) $sid jobs/nonneg_alt_summary.sbatch normal128cores continuation summary (FIT=$FD/final OUT=$ROOT/summary.txt MODEL_EXT=$SX): afterok:$gid" >> "$SUBLOG"
    now=$(date -u +%FT%TZ)
    cat > "$ROOT/README.md" <<EOF
# Stage 1 continuation: $MODEL $SET

Created $now by \`jobs/nonneg_alt_stage1c.sh\` (fit attempts ${ids[*]}, guard $gid, summary $sid).

- **Continuation of** \`$src\`: a new fit from that fit's final state (\`--start restart\`, every parameter by name
  in exact optimizer coordinates from \`$RS\`, no held term), with $NFEV more evaluations.
- **Why** (user 2026-11-14 ~17:30 CST, after the convergence diagnosis \`checks/stage1_diag_20261114/\`): B and A
  stopped at max_nfev far from a stationary point (the Gauss-Newton model valid only for small steps), and D f = 0.1
  stopped on ftol in a curved valley where scipy's cumulative x_scale='jac' scale had compressed the multiplier
  directions by more than 1000x. A new least_squares call rebuilds the scale from the Jacobian at the start point.
- **Source fit:** data cost $COST, $NF evaluations, termination "${MSG//_/ }".
- **Model settings (from the source):** \`$EXP\`.
$( [ "$MODEL" = logspread ] && [ -n "$U0" ] && printf '%s\n' "- **Spread penalty (user 2026-11-14):** about the fixed level u0 = $U0 (the stage-1 uniform start level ln s), P = f C_ref sum_seen w (u - u0)^2 = f C_ref [Var_w(u) + (u_bar - u0)^2]; recorded apart (model_check.alternative.spread: Var_w, (u_bar - u0)^2, both penalty parts). The source fit used the mean-centred penalty f C_ref Var_w(u), so the start cost of this fit includes the level term.")
- **Common scheme:** as stage 1 (lategrid base, D2's setup, ft_c1..ft_c6 held at 1, multipliers >= 1, --tbounds 30 400,
  lambda 0, crossing weight 10).
- **Code:** frozen snapshot \`$MX\`; AMPS $GIT.
- **Products:** \`fit/seg<NN>/\`, \`fit/final\`, \`fit/provenance-<job>/commands.sh\`; \`summary.txt\` (+ .json) after the
  guard (snapshot \`$SX\`). The start check in seg01's fit_manifest.json compares the start with the source's final
  cost (equal data cost expected).
EOF
    {
        echo "# How to rerun: continuation of $src"
        echo
        echo "$DESC. Model: \`srcCG/LiouvilleCOPS_cuda/model_ext/README.md\`, section \"Alternative non-negative activity"
        echo "models (options A-D)\", continuation and fixed-level penalty paragraphs. Run root: \`$ROOT\`."
        echo
        echo "## Inputs"
        echo "- the source fit \`$RS\` (sha256 $(sha256sum "$RS" | cut -d' ' -f1));"
        echo "- base, points, slope bins, D2 templates as stage 1 (\`$src/HOW-TO-RERUN.md\`);"
        echo "- code: \`$MX\` (frozen; SHA256SUMS); AMPS $GIT."
        echo
        echo "## Commands (submitted $now)"
        echo '```bash'
        echo "export COPS_OUTPUT_ROOT=$COPS_OUTPUT_ROOT"
        for k in "${!CMDS[@]}"; do echo "# attempt $((k + 1)) = job ${ids[$k]}"; echo "${CMDS[$k]}"; done
        echo "# guard = job $gid"
        echo "sbatch ${gargs[*]} $J/positive_chain_guard.sbatch"
        echo "# summary = job $sid"
        echo "sbatch ${sargs[*]} $J/nonneg_alt_summary.sbatch"
        echo '```'
        echo "All at once: \`jobs/nonneg_alt_stage1c.sh --model-ext $MX --utc <new UTC> --from $src$( [ "$MODEL" = logspread ] && [ -n "$U0" ] && echo " --spread-u0 $U0") --submit\`."
        echo "Each segment's \`fit/provenance-<job>/commands.sh\` reruns that segment with every argument resolved."
    } > "$ROOT/HOW-TO-RERUN.md"
    echo "$src -> $ROOT: attempts ${ids[*]}, guard $gid, summary $sid"
done
