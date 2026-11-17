#!/bin/bash
# ==== version-controlled copy: header added 2026-11-16, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/nonneg_stage2r_chains.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/nonneg_stage2r_chains.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 9c4672b4da00136de1d30466cff3d803fd8c2b82800f419970779b78b370dea8
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: launcher. Stage-2 reruns with the regularized late-facet fine tuning (user 2026-11-16): softplus reference
#   (start D2), background bg-min 0.01, logspread f 0.1 with u0, sos rank 3 (the last three with cycle 1 adopted from
#   stage1c); 5 cycles, --slope-mult-min 1, --ft-delta-* 10 / 4.2688066575e-03 / 4.605170185988092; roots
#   nonneg/<model>/stage2r_<setting>_<UTC>/ with README.md / HOW-TO-RERUN.md; dry run by default.
# RUN: jobs/nonneg_stage2r_chains.sh --utc 20261116T060629Z --after 135366 --submit (jobs 135367-135378).
# 2026-11-16: chain code chain_onejob_20261116T160712Z (v7.2); resubmitted --utc 20261116T160802Z --after 135862
#   (jobs 135863-135874; the first submission 135367-135378 was cancelled with its failed gate 135365).
# ==== end of header ====
# Stage-2 RERUNS with the REGULARIZED late-facet fine tuning (user decisions 2026-11-16 ~00:10 CST, via the lead):
# ONE one-job chain (jobs/chain_onejob.sh -> jobs/chain_onejob.sbatch -> chain_onejob.py, kind nonneg) per model,
# 5 fit / fine-tuning cycles each, in the order a, b, c, d:
#   a softplus   the reference: the setting of nonneg/softplus/chain_5cycles_smm1_20261113T041215Z (bounded softplus,
#                cycle-1 fit from D2) with the regularized fine tuning -- the fair reference for the option comparison
#   b background --bg-min 0.01                       (c1 fit adopted from nonneg/background/stage1c_bg1e-2_20261115T023450Z)
#   c logspread  --spread-f 0.1 --spread-cref 4.2688066575e-03 --spread-u0 -1.864126788612443
#                                                    (c1 fit adopted from nonneg/logspread/stage1c_f0.1_20261115T023450Z)
#   d sos        --sos-rank 3                        (c1 fit adopted from nonneg/sos/stage1c_r3_20261115T023450Z)
# all: --slope-mult-min 1, D2's ft templates (held at 1); fit NFEV 4000 per cycle, CHECK 1; fine tuning late_only,
# WINDOW_MIN 5, MAXNFEV 500 with --delta-penalty-f 10 --delta-cref 4.2688066575e-03 --delta-cap 4.605170185988092
# (FT_DELTA_*); L_THREADS 32 (chunked); MODEL_EXT model_ext_frozen/nonneg_20261116T031115Z.
# Chain roots: runs/shap8_012m_lategrid_dz1.5/nonneg/<model>/stage2r_<setting>_<UTC>/ with README.md and HOW-TO-RERUN.md.
# Nodes: any free GPU node but gpu003 (chain_onejob.sh always excludes it), 8 GPUs, 08:00:00, 3 attempts per chain
# (attempt a > 1: afternotok of every earlier attempt; --kill-on-invalid-dep=yes); --after JOBID (afterok, e.g. the
# bitwise split-vs-one-job test of the regularized fine tuning) gates every attempt.
#   jobs/nonneg_stage2r_chains.sh [--utc <UTC>] [--after JOBID] [--models "softplus background logspread sos"] [--submit]
# DRY RUN by default (prints the chain_onejob.sh sbatch lines; nothing is created). --submit creates the roots, writes
# README.md / HOW-TO-RERUN.md, submits (chain_onejob.sh logs every job in model_ext/logs/submissions.txt) and appends the
# job ids to HOW-TO-RERUN.md.
set -uo pipefail
ROOT_OUT=/data001/highresolu_vt/cops_output
J=$ROOT_OUT/jobs; RUN=$ROOT_OUT/runs/shap8_012m_lategrid_dz1.5; NN=$RUN/nonneg
SNAP=$ROOT_OUT/model_ext_frozen/nonneg_20261116T031115Z
CC=$ROOT_OUT/model_ext_frozen/chain_onejob_20261116T160712Z
DELTA=(--ft-delta-penalty-f 10 --ft-delta-cref 4.2688066575e-03 --ft-delta-cap 4.605170185988092)
PY="/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python -W ignore"
UTC=$(date -u +%Y%m%dT%H%M%SZ); AFTER=""; MODELS="softplus background logspread sos"; SUBMIT=0
die() { echo "nonneg_stage2r_chains.sh: $*" >&2; exit 2; }
while [ $# -gt 0 ]; do
    case $1 in
        --utc) UTC=$2; shift 2;; --after) AFTER=$2; shift 2;; --models) MODELS=$2; shift 2;;
        --submit) SUBMIT=1; shift;; -h|--help) sed -n '2,25p' "$0"; exit 0;; *) die "unknown argument $1";;
    esac
done
case $UTC in [0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]T[0-9][0-9][0-9][0-9][0-9][0-9]Z) ;; *) die "--utc YYYYmmddTHHMMSSZ";; esac
case $AFTER in ""|*[!0-9]*) [ -z "$AFTER" ] || die "--after takes ONE job id";; esac
for d in "$SNAP" "$CC"; do [ -d "$d" ] || die "missing $d"; done
AMPSC=$(git -C /data001/highresolu_vt/AMPS rev-parse HEAD)
for m in $MODELS; do
    case $m in
        softplus)   SET=smm1; SRC=""; TOK=()
                    LABEL="softplus reference (A = m_b s ln(1 + exp(L/s)), s = 1e-5 A_orig; the setting of nonneg/softplus/chain_5cycles_smm1_20261113T041215Z: cycle-1 fit from D2's final state)";;
        background) SET=bg1e-2; SRC=$NN/background/stage1c_bg1e-2_20261115T023450Z; TOK=(--bg-min 0.01)
                    LABEL="option B (background floor A_bg >= bg_min = 0.01)";;
        logspread)  SET=f0.1; SRC=$NN/logspread/stage1c_f0.1_20261115T023450Z
                    TOK=(--spread-f 0.1 --spread-cref 4.2688066575e-03 --spread-u0 -1.864126788612443)
                    LABEL="option D (log form with the spread penalty f = 0.1, C_ref = 4.2688066575e-03, fixed level u0 = -1.864126788612443)";;
        sos)        SET=r3; SRC=$NN/sos/stage1c_r3_20261115T023450Z; TOK=(--sos-rank 3)
                    LABEL="option A (sum of squares, rank 3)";;
        *) die "--models: softplus, background, logspread, sos";;
    esac
    L=$NN/$m/stage2r_${SET}_$UTC
    if [ -n "$SRC" ]; then C1=$SRC/fit; C1OPT=(--c1-from "$C1")
        [ -f "$C1/final/fit_restart.npz" ] || die "missing the stage-1c fit $C1/final/fit_restart.npz"
    else C1=$RUN/clip_refit/D2/final; C1OPT=(); fi
    [ -e "$L" ] && die "exists: $L (choose a new --utc)"
    ARGS=(--kind nonneg --model "$m" "${TOK[@]}" --slope-mult-min 1 "${C1OPT[@]}" --chain-root "$L" --model-ext "$SNAP"
          --chain-code "$CC" --cycles 5 --nfev 4000 --check 1 --ft-maxnfev 500 --ft-facets late_only --ft-window-min 5
          --l-threads 32 "${DELTA[@]}" --time 08:00:00 --attempts 3 ${AFTER:+--after "$AFTER"})
    echo "# ===== $m: $L"
    if [ "$SUBMIT" != 1 ]; then $J/chain_onejob.sh "${ARGS[@]}" | grep -E '^  sbatch|^# (one-job|cycle-1|attempt)'; continue; fi
    if [ -n "$SRC" ]; then RS=$C1/final/fit_restart.npz; else RS=$C1/fit_restart.npz; fi
    read -r SCOST SNFEV SMSG < <($PY - "$RS" <<'PYEOF'
import json, sys
import numpy as np
with np.load(sys.argv[1], allow_pickle=False) as z:
    m = json.loads(str(z["meta_json"]))
o = m.get("optimizer") or {}
print(o.get("cost"), sum(int(s.get("nfev", 0)) for s in m.get("segments", [])), o.get("message"))
PYEOF
)
    SSHA=$(sha256sum "$RS" | cut -c1-64)
    mkdir -p "$L"
    if [ -n "$SRC" ]; then
        C1TXT="- **Cycle 1's fit (adopted, not rerun):** \`$C1\` (final -> \`$(readlink -f "$C1/final")\`; fit_restart.npz sha256 $SSHA;
  optimizer cost $SCOST, $SNFEV evaluations, \"$SMSG\"). \`c1/fit/final\` links to it and \`c1/fit/ADOPTED_C1_FIT.txt\` records the
  checks (status complete, no held term, model tokens = this chain's)."
    else
        C1TXT="- **Cycle 1's fit (run here):** from D2's final state \`$C1/fit_restart.npz\` (sha256 $SSHA; cost $SCOST), --start d2,
  NFEV 4000, CHECK 1 -- as c1 of \`nonneg/softplus/chain_5cycles_smm1_20261113T041215Z\` (same fit code; only nonneg_finetune.py
  changed in the snapshot)."
    fi
    cat > "$L/README.md" <<EOF
# Stage 2 rerun (regularized fine tuning): $m $SET, 5 fit / fine-tuning cycles (one job per chain)

Created $UTC by \`jobs/nonneg_stage2r_chains.sh\` (agent 3), \`jobs/chain_onejob.sh\` -> \`jobs/chain_onejob.sbatch\` -> \`chain_onejob.py\` (kind nonneg). Job ids: HOW-TO-RERUN.md.

- **Model:** $LABEL; slope-bin multipliers >= 1 (\`--slope-mult-min 1\`), D2's ft_c1..ft_c6 templates held at 1; D2's scheme otherwise
  (lategrid base, --tbounds 30 400, lambda 0, crossing weight 10). Model description: \`srcCG/LiouvilleCOPS_cuda/model_ext/README.md\`,
  "Alternative non-negative activity models (options A-D)".
- **Decision** (user 2026-11-16 ~00:10 CST, via the lead): rerun stage 2 with a REGULARIZED late-facet fine tuning, templates on:
  --delta-penalty-f 10 --delta-cref 4.2688066575e-03 --delta-cap 4.605170185988092 (P = f C_ref sum w (h + delta)^2 / sum w over
  the tuned facets, cap |h + delta| <= X); four chains: the softplus reference (first), B bg-min 1e-2, D logspread f 0.1 with u0,
  A sos rank 3; 5 cycles, --slope-mult-min 1, fit max_nfev 4000, fine tuning late_only, WINDOW_MIN 5, MAXNFEV 500. The
  unregularized stage-2 chains (nonneg/<model>/stage2_*_20261115T184949Z) are kept as the unregularized references.
$C1TXT
- **Cycles:** c1 ft (late_only, WINDOW_MIN 5, MAXNFEV 500, regularized) -> c2 fit (--start chain: every parameter of c1's fit, grain A/B of c1 ft,
  held term = c1 ft's accumulated term; NFEV 4000, CHECK 1) -> c2 ft -> ... -> c5 ft. L_THREADS 32 (chunked).
- **Code:** model \`$SNAP\` (frozen, SHA256SUMS); chain driver \`$CC\` (frozen; v7.2: the fine-tuning delta options; all attempts gated afterok:${AFTER:-none} = the split-vs-one-job GPU test of the regularized fine tuning);
  AMPS $AMPSC.
- **Products:** \`c<k>/fit/seg<NN>\`, \`c<k>/fit/final\`, \`c<k>/ft/seg<NN>\`, \`c<k>/ft/final\` + product links; per stage
  \`c<k>/{fit,ft}/provenance-<job>-c<k>{fit,ft}/\` (commands.sh reruns that stage); per job \`provenance-<job>/\` (job script,
  code copy that ran, chain.log, chain_summary.json, git / module / nvidia-smi / toolchain / pip records);
  \`chain_onejob_events.jsonl\`. A job that reaches the 8-h limit stops at its best state (exit 75) and the next attempt continues.
- Reaching max_nfev is not convergence; compare fits by final cost, termination and active bounds.
EOF
    {
        echo "# How to rerun: stage-2 rerun (regularized fine tuning) of $m $SET (5 cycles, one job)"
        echo
        if [ -n "$SRC" ]; then echo "Run root: \`$L\`. Inputs: the adopted cycle-1 fit \`$C1/final/fit_restart.npz\` (sha256 $SSHA; its own HOW-TO-RERUN:"
            echo "\`$SRC/HOW-TO-RERUN.md\`); base, points, slope bins and D2 templates as stage 1; code \`$SNAP\` and \`$CC\`; AMPS $AMPSC."
        else echo "Run root: \`$L\`. Inputs: D2's final state \`$C1/fit_restart.npz\` (sha256 $SSHA), base, points, slope bins and D2"
            echo "templates as \`nonneg/softplus/chain_5cycles_smm1_20261113T041215Z\`; code \`$SNAP\` and \`$CC\`; AMPS $AMPSC."; fi
        echo
        echo "## Commands (submitted $UTC)"
        echo '```bash'
        echo "export COPS_OUTPUT_ROOT=$ROOT_OUT"
        printf '%q ' "$J/chain_onejob.sh" "${ARGS[@]}" --submit; echo
        echo '```'
        echo "The sbatch lines chain_onejob.sh submitted (job ids below):"
        echo '```'
    } > "$L/HOW-TO-RERUN.md"
    OUT=$($J/chain_onejob.sh "${ARGS[@]}" --submit) || { echo "$OUT"; die "chain_onejob.sh failed for $m"; }
    echo "$OUT" | grep -v '^#   '
    $J/chain_onejob.sh "${ARGS[@]}" 2>/dev/null | grep -E '^  sbatch' | sed 's/^  //' >> "$L/HOW-TO-RERUN.md"
    {
        echo '```'
        echo "$OUT" | grep -E '^attempt [0-9]+: job [0-9]+$' | sed 's/^/- /'
        echo
        echo "To continue after all attempts ended with exit 75 (wall time), rerun the same chain_onejob.sh command: completed stages are"
        echo "checked and skipped, the current stage continues with its next segment. Each stage's"
        echo "\`c<k>/{fit,ft}/provenance-<job>-c<k>{fit,ft}/commands.sh\` reruns that stage with every argument resolved; the job's"
        echo "\`provenance-<job>/commands.sh\` reruns the driver."
    } >> "$L/HOW-TO-RERUN.md"
done
[ "$SUBMIT" = 1 ] || echo "# dry run: nothing created or submitted (add --submit)"
