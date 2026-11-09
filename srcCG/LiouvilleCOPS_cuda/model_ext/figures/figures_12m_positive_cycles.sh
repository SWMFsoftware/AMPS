#!/bin/bash
# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_12m_positive_cycles.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_12m_positive_cycles.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 0e6bf83f1fde1d30fdd0dff37f7333ff4833b106e79598ea78158db83a39bfed
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure script. The positive-activity cycles against the COPS observations: final state of every lambda, per
#   lambda first-scan fit / cycle-1 fit / final; per lambda every stage (_stages); cycles_table.txt; HOW-TO-RERUN.md (every stage's commands.sh, AMPS
#   commits, code commits); code/ (frozen snapshots, job and figure scripts). 2026-11-08.
# RUN: figures_12m_positive_cycles.sh   (login node, no GPU)
# 2026-11-09 update: ratio 1e-5 (user 2026-11-08), CO2 flux floor (COPS_POS_FLOOR_CO2), start override --start-set (kappa = 1).
# ==== end of header ====
# 12M positive-activity fit / fine-tuning cycles against the COPS observations (user decisions 2026-11-07/08: night,
# shadow and floor ratio 0.002; "repeat all 6 fits starting from the best restart of each of the corresponding fit";
# "3 cycles of fitting and fine tuning"; fine tuning per facet, "only on the faces ... seen in the last 5 min only"
# (FACETS=late_only, WINDOW_MIN=5); the change held in the next fit; the fit's own lambda; no cap; MAXNFEV 500).
# Writes, into figures_cluster/12M/model_vs_observations_positive_cycles/:
#   12M_positive_cycles_final_*        the final state (cycle-3 fine tuning) of every lambda
#   12M_positive_cycles_lam<L>_*       per lambda: first-scan fit (ratio 0.02), cycle-1 fit (ratio 0.002), final state
#   12M_positive_cycles_lam<L>_stages_* per lambda: every stage (cycle 1-3 fit and fine tuning)
#   cycles_table.txt                   every stage (positive_cycles_table.py)
#   HOW-TO-RERUN.md, code/             every command with all arguments, AMPS commits, code commits, the code that ran
# The fits' and fine tunings' 401 model-only grid rows are removed (positive_table.py). A lambda whose chain is not
# complete is left out of the final figure and gets no per-lambda figure.
# RATIO (default 0.002; 2026-11-09 also 1e-5: user 2026-11-08 'decrease the floor to 1.0E-5 and also use that for nightside
# and shadows', 'use the same floor for CO2', kappa start 1 for lambda 1e-3/1e-2/1e-1): products positive_r<RATIO>,
# figures into model_vs_observations_positive_cycles (0.002) or model_vs_observations_positive_cycles_r<RATIO>; for
# 1e-5 the per-lambda figure compares ratio 0.002 after 3 cycles, the ratio-1e-5 cycle-1 fit and the ratio-1e-5 final state.
# Login node, no GPU:   [RATIO=1e-5] figures_12m_positive_cycles.sh
set -u
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; J=$COPS_OUTPUT_ROOT/jobs
RATIO=${RATIO:-0.002}; case $RATIO in 0.002|1e-5) ;; *) echo "RATIO 0.002 or 1e-5" >&2; exit 2;; esac
RUN=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; SCAN=$RUN/positive; CY=$RUN/positive_r$RATIO; CY2=$RUN/positive_r0.002
N=C30_n4000_exp_cliff_w10x5_grain_cross10_lam
O=$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_positive_cycles; [ "$RATIO" = 0.002 ] || O=${O}_r$RATIO
TB=$O/tables; C=$O/code
mkdir -p "$TB" "$C"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
PY=".venv/bin/python -W ignore"
LAMS="0 1e-5 1e-4 1e-3 1e-2 1e-1"
done_l=""
for L in $LAMS; do [ -f "$CY/lam$L/c3/ft/finetune_result.txt" ] && done_l="$done_l $L"; done
echo "complete chains:$done_l"
STY=(dotted dashed dashdot solid solid dashed)
FINAL=(); k=0
for L in $done_l; do
    if [ "$RATIO" = 0.002 ]; then
        $PY $S/positive_table.py "$SCAN/$N$L/seg01" "$TB/lam${L}_scan.dat" || exit 1
        CMP=(--fit "first scan, ratio 0.02=$TB/lam${L}_scan.dat:dotted")
    else
        $PY $S/positive_table.py "$CY2/lam$L/c3/ft" "$TB/lam${L}_r0.002final.dat" || exit 1
        CMP=(--fit "ratio 0.002 after 3 cycles=$TB/lam${L}_r0.002final.dat:dotted")
    fi
    $PY $S/positive_table.py "$CY/lam$L/c1/fit/final" "$TB/lam${L}_c1fit.dat" || exit 1
    $PY $S/positive_table.py "$CY/lam$L/c3/ft" "$TB/lam${L}_final.dat" || exit 1
    LAB="lambda $L"; [ "$L" = 0 ] && LAB="lambda 0 (fits not converged)"
    FINAL+=(--fit "$LAB, after 3 cycles=$TB/lam${L}_final.dat:${STY[$k]}"); k=$((k + 1))
    $PY $S/model_obs_figures.py --bundle $RUN/bundle/restart_templates.npz --overwrite \
        --out "$O/12M_positive_cycles_lam$L" --ref "$TB/lam${L}_final.dat" \
        --title "SHAP8 12M, positive activity, lambda $L, ratio $RATIO: $( [ $RATIO = 0.002 ] && echo 'ratio 0.02 fit' || echo 'ratio 0.002 final'), cycle-1 fit, after 3 fit + late-facet fine-tuning cycles" \
        "${CMP[@]}" \
        --fit "cycle-1 fit, ratio $RATIO=$TB/lam${L}_c1fit.dat:dashed" \
        --fit "after 3 cycles (final)=$TB/lam${L}_final.dat:solid" || exit 3
    # every stage of the chain (added 2026-11-08 on request "prepare comparison plots for them"): fits dashed,
    # fine tunings solid, cycle by colour order
    STG=()
    for c in 1 2 3; do
        $PY $S/positive_table.py "$CY/lam$L/c$c/fit/final" "$TB/lam${L}_c${c}fit.dat" > /dev/null || exit 1
        $PY $S/positive_table.py "$CY/lam$L/c$c/ft" "$TB/lam${L}_c${c}ft.dat" > /dev/null || exit 1
        STG+=(--fit "cycle $c fit=$TB/lam${L}_c${c}fit.dat:dashed" --fit "cycle $c fine tuning=$TB/lam${L}_c${c}ft.dat:solid")
    done
    $PY $S/model_obs_figures.py --bundle $RUN/bundle/restart_templates.npz --overwrite \
        --out "$O/12M_positive_cycles_lam${L}_stages" --ref "$TB/lam${L}_final.dat" \
        --title "SHAP8 12M, positive activity, lambda $L, ratio $RATIO: every fit and late-facet fine tuning of the 3 cycles" \
        "${STG[@]}" || exit 3
done
REF=$TB/lam1e-2_final.dat; [ -f "$REF" ] || REF=$(ls $TB/lam*_final.dat | head -1)
$PY $S/model_obs_figures.py --bundle $RUN/bundle/restart_templates.npz --overwrite \
    --out "$O/12M_positive_cycles_final" --ref "$REF" \
    --title 'SHAP8 12M, positive activity exp(Y a + c T$_{\rm cliff}$), ratio '"$RATIO"', after 3 cycles of fit + late-facet fine tuning' \
    "${FINAL[@]}" || exit 3
$PY $S/positive_cycles_table.py "$CY" "$SCAN" "$O/cycles_table.txt" > /dev/null || exit 4
# code: frozen snapshots, job scripts as run (provenance copies), figure scripts
SNAPS="c1_20261108T000327Z ft_20261108T032615Z ft_20261108T034023Z"; [ "$RATIO" = 1e-5 ] && SNAPS="r1e-5k_20261109T023540Z"
for snap in $SNAPS; do
    mkdir -p "$C/model_ext_frozen/$snap"; cp -rp "$COPS_OUTPUT_ROOT/model_ext_frozen/$snap/." "$C/model_ext_frozen/$snap/"
done
chmod -R u+w "$C/model_ext_frozen"
mkdir -p "$C/jobs" "$C/figures"
cp -p $(ls -d $CY/lam1e-3/c1/fit/provenance-* | tail -1)/positive_fit.sbatch "$C/jobs/positive_fit.sbatch"
cp -p $(ls -d $CY/lam1e-3/c1/ft/provenance-* | tail -1)/positive_finetune.sbatch "$C/jobs/positive_finetune.sbatch"
cp -p $J/positive_cycles.sh "$C/jobs/"
cp -p $S/figures_12m_positive_cycles.sh $S/positive_table.py $S/positive_cycles_table.py $S/model_obs_figures.py \
    $S/rerun_commands.py "$C/figures/"
find "$C" -name __pycache__ -prune -exec rm -rf {} + 2>/dev/null
(cd "$C" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
# HOW-TO-RERUN.md
prov() { ls -d "$1"/provenance-* 2>/dev/null | sort -t- -k2 -n | tail -1; }   # the last job of a stage (the one that completed)
CM=$($PY - "$C" <<'PY'
import hashlib, sys
from pathlib import Path
sys.path.insert(0, "/data001/highresolu_vt/cops_output/figures_cluster/scripts")
from rerun_commands import commit_of
for p in sorted(Path(sys.argv[1]).rglob("*")):
    if p.is_file() and p.suffix in (".py", ".sbatch", ".sh", ".cu"):
        h = hashlib.sha256(p.read_bytes()).hexdigest()
        c = commit_of(p.name, h)
        print("| `%s` | `%s` | %s |" % (p.relative_to(sys.argv[1]), h[:16], "`%s`" % c if c else "**not in git**"))
PY
)
if [ "$RATIO" = 0.002 ]; then
WHAT="- Night/shadow ratio and H2O floor ratio 0.002 (before: 0.02); CO2 no floor.
- Cycle 1 fit: started from the first scan's fit of the same lambda (ratio 0.02), \`--start-restart\`.
- Code: \`model_ext_frozen/c1_20261108T000327Z\` for every fit, \`ft_20261108T032615Z\` for the fine tunings,
  \`ft_20261108T034023Z\` for the lambda-0 fine tunings (the same code plus the zero-sensitivity handling)."
A2="Cycle-1 fits (2026-11-08T00:17Z), then the chains from the cycle-1 fine tuning (2026-11-08T03:33Z; lambda 0 again
at 04:12Z with the zero-sensitivity snapshot):

\`\`\`bash
cd $J
F=$COPS_OUTPUT_ROOT/model_ext_frozen/c1_20261108T000327Z; P=$SCAN
for L in 1e-3 1e-2 1e-4 1e-1 1e-5 0; do
  sbatch --nodelist=gpu001 --export=ALL,MODEL_EXT=\$F,LAMBDA=\$L,RATIO=0.002,CYCLE=1,START_RESTART=\$P/$N\$L/seg01/fit_restart.npz positive_fit.sbatch
done
FT=$COPS_OUTPUT_ROOT/model_ext_frozen/ft_20261108T032615Z     # lambda 0: ft_20261108T034023Z
./positive_cycles.sh --lambda <L> --ratio 0.002 --model-ext \$F --ft-model-ext \$FT --fit-node <gpu001|gpu002> --ft-node <same> \\
    --start-at c1ft --ft-maxnfev 500 --ft-time 02:00:00 --ft-export \"FACETS=late_only,WINDOW_MIN=5\" --submit
\`\`\`"
else
WHAT="- Night/shadow ratio (H2O and CO2), H2O floor and CO2 floor all 1e-5 (user 2026-11-08: 'decrease the floor to 1.0E-5 and
  also use that for nightside and shadows', 'use the same floor for CO2'): per record F_CO2 = max(kappa Q_CO2 ACO2,
  r qmax_C ACO2) (ACO2 = the run's CO2 map), switch COPS_POS_FLOOR_CO2=1.
- Cycle 1 fit: started from the ratio-0.002 cycle-1 fit of the same lambda (\`--start-restart\`); for lambda 1e-3, 1e-2,
  1e-1 kappa started at 1 (\`--start-set kappa=1\`; user 2026-11-08 'kappa = 1 where stuck': below kappa ~1e-5 every
  CO2 record is at the floor and d/dkappa = 0).
- Code: \`model_ext_frozen/r1e-5k_20261109T023540Z\` for every fit and fine tuning."
A2="All six chains, 2026-11-09T02:48Z (exact lines, \`model_ext/logs/agent_kappastart_six_lines.txt\`):

\`\`\`bash
cd $J
$(cat $COPS_OUTPUT_ROOT/model_ext/logs/agent_kappastart_six_lines.txt)
\`\`\`"
fi
{
cat <<MD
# How to rerun: 12M positive-activity fit / fine-tuning cycles, ratio $RATIO

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_12m_positive_cycles.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never delete
earlier attempts (move them aside as \`*.incomplete_<UTC>\`).

**What was run** (user decisions 2026-11-07 and 2026-11-08):
$WHAT
- Activity A = exp(Y a + c T_cliff + held), D on the
  total ln A from the prior s0, the original temperature law held, kappa free, grain + crossing weight 10, W10x5 for
  the observations from 10:33:10, gauge +1.5 m along +Z, max_nfev 4000.
- Fine tuning: per facet, A_e = exp(s_e + delta_e), only the facets seen exclusively in the last 5 min
  (window 10:34:26-10:39:26 UTC; 271 facets, 355.6 m2 on this base); delta >= ln($RATIO A_orig / A) (floor), no cap;
  the fit's own lambda; grain A, B refitted; MAXNFEV 500. Earlier pressures are unchanged exactly (strict support).
- Cycles 2 and 3: the fit starts from the previous fit's restart, holds the accumulated fine-tuning change
  (\`--held-exponent\`), takes the grain from the fine tuning.
- Base: \`$RUN/base\` (1150 observations + 401 model-only 1-s grid points; grid rows removed in the figures).
- Complete chains:$done_l. Stage-by-stage numbers: \`cycles_table.txt\`.

Code: the "commit" column names the git commit (in \`srcCG/LiouvilleCOPS_cuda/model_ext\`) whose version-controlled
copy equals the file that ran ("not in git": only in \`code/\` here); \`code/\` has the frozen snapshots the jobs ran,
with \`SHA256SUMS\`. The commands below are complete (every argument as run).

## A0. The figures only (login node, no GPU)

\`\`\`bash
RATIO=$RATIO $S/figures_12m_positive_cycles.sh
\`\`\`

## A1. The base and the first scan

The base (job 122291) and the six first-scan fits: \`$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_positive_lambda_scan/HOW-TO-RERUN.md\`.

## A2. Chain submission (as run)

$A2

## A3. Every stage, commands as run

MD
for L in $done_l; do
    echo "### lambda $L"; echo
    for c in 1 2 3; do
        for st in fit ft; do
            d=$CY/lam$L/c$c/$st; p=$(prov "$d"); [ -n "$p" ] || continue
            echo "**cycle $c $( [ $st = fit ] && echo fit || echo 'fine tuning') ($(basename $p)):** AMPS \`$(head -1 $p/git.txt)\`"
            echo; echo '```bash'; cat "$p/commands.sh"; echo '```'; echo
        done
    done
done
cat <<MD
## Code that ran (\`code/\`, with \`SHA256SUMS\`) and the git commit with an identical copy

| file | SHA256 (first 16) | commit |
|---|---|---|
$CM
MD
} > "$O/HOW-TO-RERUN.md"
echo "wrote $O/HOW-TO-RERUN.md"
