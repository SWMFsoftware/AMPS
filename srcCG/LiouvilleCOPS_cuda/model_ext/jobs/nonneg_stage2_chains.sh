#!/bin/bash
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/nonneg_stage2_chains.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/nonneg_stage2_chains.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 5ba3b8e7e00f2791b0e078fe55e95a3be5c792c6bd51cb4a3aa1454953ae39fb
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: launcher. Stage 2 of the alternative non-negative models (user 2026-11-15): one one-job chain per model
#   (background bg-min 0.01, logspread f 0.1 with the fixed-level u0, sos rank 3), 5 fit / fine-tuning cycles,
#   --slope-mult-min 1, cycle 1's fit adopted from the stage-1c fit (--c1-from); creates the chain roots
#   nonneg/<model>/stage2_<setting>_<UTC>/ with README.md / HOW-TO-RERUN.md; dry run by default.
# RUN: jobs/nonneg_stage2_chains.sh --utc 20261115T184949Z --after 131988 --submit (jobs 131989-131997).
# ==== end of header ====
# Stage 2 of the alternative non-negative activity models (user decision 2026-11-15, via the lead): ONE one-job chain
# (jobs/chain_onejob.sh -> jobs/chain_onejob.sbatch -> chain_onejob.py, kind nonneg) per model, 5 fit / fine-tuning
# cycles each, in the order B, D, A:
#   B background --bg-min 0.01                                     (from nonneg/background/stage1c_bg1e-2_20261115T023450Z)
#   D logspread  --spread-f 0.1 --spread-cref 4.2688066575e-03 --spread-u0 -1.864126788612443  (fixed-level u0;
#                                                                     from nonneg/logspread/stage1c_f0.1_20261115T023450Z)
#   A sos        --sos-rank 3                                      (from nonneg/sos/stage1c_r3_20261115T023450Z)
# all with --slope-mult-min 1 and D2's ft templates (held at 1); fit NFEV 4000 per cycle, CHECK 1; fine tuning late_only,
# WINDOW_MIN 5, MAXNFEV 500; L_THREADS 32 (chunked). Cycle 1's fit is the existing complete stage-1 continuation fit,
# ADOPTED (--c1-from: c1/fit/final -> its final segment, c1/fit/ADOPTED_C1_FIT.txt): each chain runs c1 ft, then
# c2..c5 fit (--start chain) + ft.
# Chain roots: runs/shap8_012m_lategrid_dz1.5/nonneg/<model>/stage2_<setting>_<UTC>/ with README.md and HOW-TO-RERUN.md.
# Nodes: any free GPU node but gpu003 (chain_onejob.sh always excludes it; no --nodelist), 8 GPUs, 08:00:00, 3 attempts
# per chain (attempt a > 1: afternotok of every earlier attempt; --kill-on-invalid-dep=yes); --after JOBID (afterok,
# e.g. the bitwise split-vs-one-job test) gates every attempt.
#   jobs/nonneg_stage2_chains.sh [--utc <UTC>] [--after JOBID] [--models "background logspread sos"] [--submit]
# DRY RUN by default (prints the chain_onejob.sh sbatch lines; nothing is created). --submit creates the roots, writes
# README.md / HOW-TO-RERUN.md, submits (chain_onejob.sh logs every job in model_ext/logs/submissions.txt) and appends the
# job ids to HOW-TO-RERUN.md.
set -uo pipefail
ROOT_OUT=/data001/highresolu_vt/cops_output
J=$ROOT_OUT/jobs; RUN=$ROOT_OUT/runs/shap8_012m_lategrid_dz1.5; NN=$RUN/nonneg
SNAP=$ROOT_OUT/model_ext_frozen/nonneg_20261115T020314Z
CC=$ROOT_OUT/model_ext_frozen/chain_onejob_20261115T184608Z
PY="/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python -W ignore"
UTC=$(date -u +%Y%m%dT%H%M%SZ); AFTER=""; MODELS="background logspread sos"; SUBMIT=0
die() { echo "nonneg_stage2_chains.sh: $*" >&2; exit 2; }
while [ $# -gt 0 ]; do
    case $1 in
        --utc) UTC=$2; shift 2;; --after) AFTER=$2; shift 2;; --models) MODELS=$2; shift 2;;
        --submit) SUBMIT=1; shift;; -h|--help) sed -n '2,24p' "$0"; exit 0;; *) die "unknown argument $1";;
    esac
done
case $UTC in [0-9][0-9][0-9][0-9][0-9][0-9][0-9][0-9]T[0-9][0-9][0-9][0-9][0-9][0-9]Z) ;; *) die "--utc YYYYmmddTHHMMSSZ";; esac
case $AFTER in ""|*[!0-9]*) [ -z "$AFTER" ] || die "--after takes ONE job id";; esac
for d in "$SNAP" "$CC"; do [ -d "$d" ] || die "missing $d"; done
AMPSC=$(git -C /data001/highresolu_vt/AMPS rev-parse HEAD)
for m in $MODELS; do
    case $m in
        background) SET=bg1e-2; SRC=$NN/background/stage1c_bg1e-2_20261115T023450Z; TOK=(--bg-min 0.01)
                    LABEL="option B (background floor A_bg >= bg_min = 0.01)";;
        logspread)  SET=f0.1; SRC=$NN/logspread/stage1c_f0.1_20261115T023450Z
                    TOK=(--spread-f 0.1 --spread-cref 4.2688066575e-03 --spread-u0 -1.864126788612443)
                    LABEL="option D (log form with the spread penalty f = 0.1, C_ref = 4.2688066575e-03, fixed level u0 = -1.864126788612443)";;
        sos)        SET=r3; SRC=$NN/sos/stage1c_r3_20261115T023450Z; TOK=(--sos-rank 3)
                    LABEL="option A (sum of squares, rank 3)";;
        *) die "--models: background, logspread, sos";;
    esac
    C1=$SRC/fit; L=$NN/$m/stage2_${SET}_$UTC
    [ -f "$C1/final/fit_restart.npz" ] || die "missing the stage-1c fit $C1/final/fit_restart.npz"
    [ -e "$L" ] && die "exists: $L (choose a new --utc)"
    ARGS=(--kind nonneg --model "$m" "${TOK[@]}" --slope-mult-min 1 --c1-from "$C1" --chain-root "$L" --model-ext "$SNAP"
          --chain-code "$CC" --cycles 5 --nfev 4000 --check 1 --ft-maxnfev 500 --ft-facets late_only --ft-window-min 5
          --l-threads 32 --time 08:00:00 --attempts 3 ${AFTER:+--after "$AFTER"})
    echo "# ===== $m: $L"
    if [ "$SUBMIT" != 1 ]; then $J/chain_onejob.sh "${ARGS[@]}" | grep -E '^  sbatch|^# (one-job|cycle-1|attempt)'; continue; fi
    read -r SCOST SNFEV SMSG < <($PY - "$C1/final/fit_restart.npz" <<'PYEOF'
import json, sys
import numpy as np
with np.load(sys.argv[1], allow_pickle=False) as z:
    m = json.loads(str(z["meta_json"]))
o = m.get("optimizer") or {}
print(o.get("cost"), sum(int(s.get("nfev", 0)) for s in m.get("segments", [])), o.get("message"))
PYEOF
)
    SSHA=$(sha256sum "$C1/final/fit_restart.npz" | cut -c1-64)
    mkdir -p "$L"
    cat > "$L/README.md" <<EOF
# Stage 2: $m $SET, 5 fit / fine-tuning cycles (one job per chain)

Created $UTC by \`jobs/nonneg_stage2_chains.sh\` (agent 3), \`jobs/chain_onejob.sh\` -> \`jobs/chain_onejob.sbatch\` -> \`chain_onejob.py\` (kind nonneg). Job ids: HOW-TO-RERUN.md.

- **Model:** $LABEL; slope-bin multipliers >= 1 (\`--slope-mult-min 1\`), D2's ft_c1..ft_c6 templates held at 1; D2's scheme otherwise
  (lategrid base, --tbounds 30 400, lambda 0, crossing weight 10). Model description: \`srcCG/LiouvilleCOPS_cuda/model_ext/README.md\`,
  "Alternative non-negative activity models (options A-D)".
- **Decision** (user 2026-11-15, via the lead): stage 2 = 5 fit / fine-tuning cycles each for B (bg-min 1e-2), D (logspread f 0.1
  with the fixed-level u0) and A (sos rank 3), all with --slope-mult-min 1 and D2's ft templates; fit max_nfev 4000 per cycle;
  fine tuning late_only, WINDOW_MIN 5, MAXNFEV 500. Cycle 1's fit = the existing complete stage-1 continuation fit.
- **Cycle 1's fit (adopted, not rerun):** \`$C1\` (final -> \`$(readlink -f "$C1/final")\`; fit_restart.npz sha256 $SSHA;
  optimizer cost $SCOST, $SNFEV evaluations, "$SMSG"). \`c1/fit/final\` links to it and \`c1/fit/ADOPTED_C1_FIT.txt\` records the
  checks (status complete, no held term, model tokens = this chain's).
- **Cycles:** c1 ft (late_only, WINDOW_MIN 5, MAXNFEV 500) -> c2 fit (--start chain: every parameter of c1's fit, grain A/B of c1 ft,
  held term = c1 ft's accumulated term; NFEV 4000, CHECK 1) -> c2 ft -> ... -> c5 ft. L_THREADS 32 (chunked).
- **Code:** model \`$SNAP\` (frozen, SHA256SUMS); chain driver \`$CC\` (frozen; v7, GPU-tested bitwise against the split job scripts);
  AMPS $AMPSC.
- **Products:** \`c<k>/fit/seg<NN>\`, \`c<k>/fit/final\`, \`c<k>/ft/seg<NN>\`, \`c<k>/ft/final\` + product links; per stage
  \`c<k>/{fit,ft}/provenance-<job>-c<k>{fit,ft}/\` (commands.sh reruns that stage); per job \`provenance-<job>/\` (job script,
  code copy that ran, chain.log, chain_summary.json, git / module / nvidia-smi / toolchain / pip records);
  \`chain_onejob_events.jsonl\`. A job that reaches the 8-h limit stops at its best state (exit 75) and the next attempt continues.
- Reaching max_nfev is not convergence; compare fits by final cost, termination and active bounds.
EOF
    {
        echo "# How to rerun: stage 2 of $m $SET (5 cycles, one job)"
        echo
        echo "Run root: \`$L\`. Inputs: the adopted cycle-1 fit \`$C1/final/fit_restart.npz\` (sha256 $SSHA; its own HOW-TO-RERUN:"
        echo "\`$SRC/HOW-TO-RERUN.md\`); base, points, slope bins and D2 templates as stage 1; code \`$SNAP\` and \`$CC\`; AMPS $AMPSC."
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
