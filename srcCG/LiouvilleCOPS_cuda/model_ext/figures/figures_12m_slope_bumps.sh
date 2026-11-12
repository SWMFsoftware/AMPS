#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_12m_slope_bumps.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_12m_slope_bumps.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 91dbddaf9022581a0a83642c2abb94bde9390a1f727ac0f333449fdcb9930173
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure script. The bump + slope-multiplier chain (option 3a, 3 cycles, jobs 124008-124086; user decisions
#   2026-11-10) against the COPS observations next to the harmonic slope-multiplier fit, D2 and the best positive fit;
#   per-stage table, c1 / c3 multiplier and bump tables, summary.txt, HOW-TO-RERUN.md, code/. 2026-11-10.
# RUN: figures_12m_slope_bumps.sh   (login node, no GPU) -> figures_cluster/12M/model_vs_observations_slope_bumps/
# ==== end of header ====
# The bump + slope-multiplier chain (option 3a; user decisions 2026-11-10: "implement option 3a with slope multipliers",
# "add fine tuning -- fitting cycles as in other fits", "execute the fit using the optimal number of threads", "Fix the
# check step") against the COPS observations, next to the slope-multiplier fit with the harmonic basis (srcCG start), the
# clipped-model refit D2 and the best positive-activity fit (ratio 1e-5, lambda 1e-4, after 3 cycles); lategrid base
# (gauge +1.5 m, W10x5). Writes figures_cluster/12M/model_vs_observations_slope_bumps/: 12M_slope_bumps_*, the per-stage
# cycle table, the c1/c3 multiplier and bump tables (copies), summary.txt, HOW-TO-RERUN.md, code/.
#     Login node, no GPU:   figures_12m_slope_bumps.sh
set -u
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; J=$COPS_OUTPUT_ROOT/jobs
RUN=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
CH=$RUN/slope_bumps/chain_slopemul_20261110T201212Z
BF=$CH/c3/ft/final                                   # final state of the chain (c3 fine tuning)
SM=$RUN/slope_mult/fit/final; D2=$RUN/clip_refit/D2/final; PB=$RUN/positive_r1e-5/lam1e-4/c3/ft
V6=$COPS_OUTPUT_ROOT/model_ext_frozen/slopemul_20261110T201212Z    # code of the c1 fit (124008)
V8=$COPS_OUTPUT_ROOT/model_ext_frozen/slopemul_20261110T221713Z    # code of every later stage
O=$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_slope_bumps; TB=$O/tables; C=$O/code
mkdir -p "$TB" "$C"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
PY=".venv/bin/python -W ignore"
$PY $S/positive_table.py "$BF" "$TB/slope_bumps_c3ft.dat" || exit 1
$PY $S/positive_table.py "$SM" "$TB/slope_mult.dat" || exit 1
$PY $S/positive_table.py "$D2" "$TB/clip_D2.dat" || exit 1
$PY $S/positive_table.py "$PB" "$TB/positive_r1e-5_lam1e-4.dat" || exit 1
$PY $S/model_obs_figures.py --bundle $RUN/bundle/restart_templates.npz --overwrite --out "$O/12M_slope_bumps" \
    --ref "$TB/slope_bumps_c3ft.dat" \
    --title "SHAP8 12M (gauge +1.5 m, W10x5): 32 Gaussian bumps x slope multipliers, floor 1e-5 A_orig, after 3 fit / fine-tuning cycles" \
    --fit "bumps x slope multipliers, 3 cycles=$TB/slope_bumps_c3ft.dat:solid" \
    --fit "slope multipliers x harmonics (srcCG start)=$TB/slope_mult.dat:dashed" \
    --fit "clipped refit D2 (cliff + 6 templates)=$TB/clip_D2.dat:dotted" \
    --fit "positive, ratio 1e-5, lambda 1e-4, 3 cycles=$TB/positive_r1e-5_lam1e-4.dat:dashdot" || exit 3
for c in c1 c3; do
    cp -p "$CH/$c/fit/final/slope_mult_bins.txt" "$O/${c}_fit_slope_mult_bins.txt"
    cp -p "$CH/$c/fit/final/slope_bumps_table.txt" "$O/${c}_fit_slope_bumps_table.txt"
done
cp -p "$CH/c1/fit/final/SLOPE_MULTIPLIER_MODEL.txt" "$O/"
# per-stage table of the chain
{
echo "# stage  job  data_cost  rms_ln_NG  rms_ln_RG  stop  evaluations   (chain $CH)"
for c in c1 c2 c3; do
    f=$CH/$c/fit/final/fit_result.txt; j=$(basename $(ls -d $CH/$c/fit/provenance-* | tail -1) | sed 's/provenance-//')
    printf "%-6s %-7s %s %s %s %s %s\n" "${c}fit" "$j" $(awk '/^final_data_cost/{print $2} /^rms_ln_residual_NG/{print $2} /^rms_ln_residual_RG/{print $2}' $f) \
        "$(grep -oE 'ftol|maximum number' $f | head -1 | sed 's/maximum number/max_nfev/')" "$(grep -oE 'nfev [0-9]+' $f | head -1 | cut -d' ' -f2)"
    f=$CH/$c/ft/final/finetune_result.txt; j=$(basename $(ls -d $CH/$c/ft/provenance-* | tail -1) | sed 's/provenance-//')
    printf "%-6s %-7s %s %s %s %s %s\n" "${c}ft" "$j" $(awk '/^cost_data/{print $5} /^rms_ln_NG/{print $5} /^rms_ln_RG/{print $5}' $f) \
        "$(awk '/^stop /{print $2}' $f)" "$(awk '/^nfev_total/{print $2}' $f)"
done
} > "$O/chain_stages.txt"
{
echo "# model  data_cost  rms_ln_NG  rms_ln_RG  stop"
printf "%-26s %s %s %s %s\n" "slope_bumps_3cycles" $(awk '/^cost_data/{print $5} /^rms_ln_NG/{print $5} /^rms_ln_RG/{print $5}' $BF/finetune_result.txt) "$(awk '/^stop /{print $2}' $BF/finetune_result.txt)"
for pair in "slope_mult:$SM/fit_result.txt" "clip_D2:$D2/fit_result.txt"; do n=${pair%%:*}; f=${pair#*:}
    printf "%-26s %s %s %s %s\n" "$n" $(awk '/^final_data_cost/{print $2} /^rms_ln_residual_NG/{print $2} /^rms_ln_residual_RG/{print $2}' $f) \
        "$(grep -oE 'ftol|maximum number' $f | head -1 | sed 's/maximum number/max_nfev/')"
done
printf "%-26s %s %s %s %s\n" "positive_r1e-5_lam1e-4" $(awk '/^cost_data/{print $5} /^rms_ln_NG/{print $5} /^rms_ln_RG/{print $5}' $PB/finetune_result.txt) "$(awk '/^stop /{print $2}' $PB/finetune_result.txt)"
} > "$O/summary.txt"
# code: both snapshots, the job scripts as run, the chain and guard scripts, the bump builder, the figure scripts
for z in "$V6" "$V8"; do mkdir -p "$C/model_ext_frozen/$(basename $z)"; cp -rp "$z/." "$C/model_ext_frozen/$(basename $z)/"; done
chmod -R u+w "$C/model_ext_frozen"; find "$C" -name __pycache__ -prune -exec rm -rf {} + 2>/dev/null
mkdir -p "$C/jobs" "$C/figures" "$C/bumps"
cp -p "$CH/c1/fit/provenance-124008/slope_mult_fit.sbatch" "$C/jobs/"
cp -p "$(ls -d $CH/c1/ft/provenance-* | tail -1)/slope_finetune.sbatch" "$C/jobs/"
cp -p $J/slope_cycles.sh $J/slope_chain_guard.sbatch $J/slope_bumps_test.sbatch "$C/jobs/"
cp -p $COPS_OUTPUT_ROOT/checks/slopebumps_20261110/slope_bumps.py "$C/bumps/"
cp -p $S/figures_12m_slope_bumps.sh $S/positive_table.py $S/model_obs_figures.py $S/rerun_commands.py "$C/figures/"
(cd "$C" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
CM=$($PY - "$C" <<'PY'
import hashlib, sys
from pathlib import Path
sys.path.insert(0, "/data001/highresolu_vt/cops_output/figures_cluster/scripts")
from rerun_commands import commit_of
for p in sorted(Path(sys.argv[1]).rglob("*")):
    if p.is_file() and p.suffix in (".py", ".sbatch", ".sh", ".cu"):
        h = hashlib.sha256(p.read_bytes()).hexdigest(); c = commit_of(p.name, h)
        print("| `%s` | `%s` | %s |" % (p.relative_to(sys.argv[1]), h[:16], "`%s`" % c if c else "**not in git**"))
PY
)
{
cat <<MD
# How to rerun: 12M bump + slope-multiplier chain (option 3a), 3 fit / fine-tuning cycles

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_12m_slope_bumps.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never delete
earlier attempts (move them aside as \`*.incomplete_<UTC>\`). AMPS commit of every stage: \`$(head -1 $CH/c1/fit/provenance-124008/git.txt)\`.

**Model** (user decisions 2026-11-10): H2O activity A = max(1e-5 A_orig, m_b L), L = A_scale sum_k w_k phi_k with 32
Gaussian bumps phi_k = exp(-|x - c_k|^2 / (2 s_k^2)) (26 centres on the seen surface, 6 on the unseen, s_k = the
nearest-centre distance, w_k >= 0), one multiplier m_b >= 0 per effective-gravitational-slope bin [0,40) (held at 1),
[40,50), [50,60), [60,70), [70,80), [80,90), [90,180]; lategrid base (gauge +1.5 m), W10x5 from 10:33:10, original
temperature law held, kappa free, grain + crossing 10, lambda 0, max_nfev 4000, L step with 32 chunked threads (bitwise).
Start of cycle 1: w = the NNLS fit to the original map, multipliers 1, kappa 1, grain from D2. Fine tuning: per-facet
factors on the facets seen only in the last 5 min (late_only), down to the floor, no cap, grain refitted, MAXNFEV 500;
cycle k+1 starts from cycle k's fit restart with the fine tuning's held term and grain. Stage table: \`chain_stages.txt\`;
multipliers and bump weights of the c1 and c3 fits: \`c1_fit_*\`, \`c3_fit_*\`; model text: \`SLOPE_MULTIPLIER_MODEL.txt\`.
Compared: the slope-multiplier fit with the harmonic basis (\`$SM\`; \`model_vs_observations_slope_mult/HOW-TO-RERUN.md\`),
D2 (\`$D2\`) and the positive fit ratio 1e-5, lambda 1e-4 after 3 cycles (\`$PB\`).

Code: the c1 fit ran snapshot \`$(basename $V6)\` (job 124008); every later stage ran \`$(basename $V8)\`, which differs
from it only in the fine tuning (worker request names, the pressure-sized steps of validation checks 6 and 7; user
2026-11-10 "Fix the check step") and the tests; the fit code is identical.

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_12m_slope_bumps.sh
\`\`\`

## A1. The bump bundle (CPU job 123981)

\`\`\`bash
$(sacct -j 123981 -o SubmitLine%2000 -n -X -P 2>/dev/null | head -1)
\`\`\`
Output: \`$RUN/slope_bumps/bundle/bumps32.npz\` (+ \`.json\`: centres, widths, coverage, NNLS start).

## A2. The chain as submitted

\`\`\`bash
# cycle-1 fit (124008; submitted by hand from the dry run of the line below with --model-ext $(basename $V6) and without --start-at)
$J/slope_cycles.sh --model-ext $V6 --l-threads 32 --cycles 3 --nfev 4000 --ft-maxnfev 500 --node gpu002   # dry run; its c1fit lines were submitted
# GPU test of the fine tuning (124063), then the rest of the chain behind it (guard layout):
sbatch --nodelist=gpu002 --export=ALL,STEP=ft,MODEL_EXT=$V8,TAG=$(basename $V8 | sed 's/slopemul_//') $J/slope_bumps_test.sbatch
$J/slope_cycles.sh --model-ext $V8 --chain $CH --start-at c1ft --after <test job> --l-threads 32 --cycles 3 --nfev 4000 --ft-maxnfev 500 --node gpu002 --submit
\`\`\`

## A3. Commands of every stage, as run (inside its Slurm job)
MD
for s in c1/fit c1/ft c2/fit c2/ft c3/fit c3/ft; do
    p=$(ls -d $CH/$s/provenance-* | tail -1)
    printf '\n### %s (%s)\n\n```bash\n%s\n```\n' "$s" "$(basename $p)" "$(cat $p/commands.sh)"
done
cat <<MD

## Code that ran (\`code/\`, with \`SHA256SUMS\`) and the git commit with an identical copy

| file | SHA256 (first 16) | commit |
|---|---|---|
$CM
MD
} > "$O/HOW-TO-RERUN.md"
echo "wrote $O"; cat "$O/chain_stages.txt"; cat "$O/summary.txt"
