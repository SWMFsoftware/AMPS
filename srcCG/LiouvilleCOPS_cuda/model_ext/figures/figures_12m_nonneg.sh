#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_12m_nonneg.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_12m_nonneg.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: a44d0e2ff778dfad4e58c0c53a696b758d64a82fdf7da55bb4610e42bf71ed94
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure script. The softplus and squared non-negative models (user decisions 2026-11-10/11; 3 cycles each in one
#   job, 124162 / 124165) against the COPS observations next to D2 and the best positive fit; since 2026-11-11 19:26 UTC
#   also every stage of each chain (fits dashed, fine tunings solid), cycle 1 vs final, and the best fit of every family
#   (user: "prepare comparison plots based on the recent fittings"); chain_stages.txt, summary.txt, bin tables,
#   HOW-TO-RERUN.md, code/.
# RUN: figures_12m_nonneg.sh   (login node, no GPU) -> figures_cluster/12M/model_vs_observations_nonneg/
# ==== end of header ====
# The two non-negative-by-construction activity models (user decisions 2026-11-10/11: "implement option 1 and 2 and
# execute them concurrently ... use 6 cliff templates as multipliers; add fine tuning cycles at the end. submit fitting +
# tuning in a single job and use multiple thread parallelization on the CPU") against the COPS observations, next to
# the clipped refit D2 and the best positive-activity fit; lategrid base (gauge +1.5 m, W10x5). Each model: 3 fit /
# late-facet fine-tuning cycles in ONE job (model_ext/chain_onejob.py, kind nonneg).
#   softplus: A = m_b * s ln(1 + exp(L/s)), L = Y a + A_s sum ft_c1..6 (held 1), s = 1e-5 A_orig (clip where A_orig = 0)
#   squared:  A = m_b * (A_s (Y b)^2 + A_s sum ft_c1..6)
#   m_b: 7 slope bins, [0,40) held at 1, six free >= 0; no floor; start from D2.
# Writes figures_cluster/12M/model_vs_observations_nonneg/: 12M_nonneg_*, chain_stages.txt, summary.txt, the c3 bin
# tables, HOW-TO-RERUN.md, code/.        Login node, no GPU:   figures_12m_nonneg.sh
set -u
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; J=$COPS_OUTPUT_ROOT/jobs
RUN=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; F12=$COPS_OUTPUT_ROOT/figures_cluster/12M
TAG=chain_onejob_nonneg_20261111T044920Z
SP=$RUN/nonneg/softplus/$TAG; SQ=$RUN/nonneg/squared/$TAG
D2=$RUN/clip_refit/D2/final; PB=$RUN/positive_r1e-5/lam1e-4/c3/ft
NX=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_20261111T044920Z
CX=$COPS_OUTPUT_ROOT/model_ext_frozen/chain_onejob_20261111T060742Z
O=$F12/model_vs_observations_nonneg; TB=$O/tables; C=$O/code
mkdir -p "$TB" "$C"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
PY=".venv/bin/python -W ignore"
$PY $S/positive_table.py "$SP/c3/ft/final" "$TB/softplus_c3ft.dat" || exit 1
$PY $S/positive_table.py "$SQ/c3/ft/final" "$TB/squared_c3ft.dat" || exit 1
$PY $S/positive_table.py "$D2" "$TB/clip_D2.dat" || exit 1
$PY $S/positive_table.py "$PB" "$TB/positive_r1e-5_lam1e-4_c3.dat" || exit 1
cost() { awk '/^cost_data/{printf "%.2e", $5}' "$1/finetune_result.txt"; }
$PY $S/model_obs_figures.py --bundle $RUN/bundle/restart_templates.npz --overwrite --out "$O/12M_nonneg" \
    --ref "$TB/clip_D2.dat" \
    --title "SHAP8 12M (gauge +1.5 m, W10x5): non-negative activity models x slope multipliers, 3 cycles each, vs D2" \
    --fit "softplus x slope multipliers, 3 cycles, $(cost $SP/c3/ft/final)=$TB/softplus_c3ft.dat:solid" \
    --fit "squared harmonics x slope multipliers, 3 cycles, $(cost $SQ/c3/ft/final)=$TB/squared_c3ft.dat:dashed" \
    --fit "clipped refit D2, 4.89e-03=$TB/clip_D2.dat:dotted" \
    --fit "positive, ratio 1e-5, lambda 1e-4, 3 cycles, 8.54e-03=$TB/positive_r1e-5_lam1e-4_c3.dat:dashdot" || exit 3
# every stage of each chain (user 2026-11-11 "prepare comparison plots based on the recent fittings"): fits dashed,
# fine tunings solid; and cycle 1 vs final against D2
for m in softplus squared; do
    CH=$RUN/nonneg/$m/$TAG; STG=()
    for c in 1 2 3; do
        $PY $S/positive_table.py "$CH/c$c/fit/final" "$TB/${m}_c${c}fit.dat" || exit 1
        $PY $S/positive_table.py "$CH/c$c/ft/final" "$TB/${m}_c${c}ft.dat" || exit 1
        STG+=(--fit "cycle $c fit, $(awk '/^final_data_cost/{printf "%.3e", $2}' $CH/c$c/fit/final/fit_result.txt)=$TB/${m}_c${c}fit.dat:dashed"
              --fit "cycle $c fine tuning, $(cost $CH/c$c/ft/final)=$TB/${m}_c${c}ft.dat:solid")
    done
    $PY $S/model_obs_figures.py --bundle $RUN/bundle/restart_templates.npz --overwrite --out "$O/12M_nonneg_${m}_stages" \
        --ref "$TB/${m}_c3ft.dat" \
        --title "SHAP8 12M (gauge +1.5 m, W10x5): $m x slope multipliers, every stage of the 3 fit / fine-tuning cycles" \
        "${STG[@]}" || exit 3
    $PY $S/model_obs_figures.py --bundle $RUN/bundle/restart_templates.npz --overwrite --out "$O/12M_nonneg_${m}_c1_vs_final" \
        --ref "$TB/${m}_c3ft.dat" \
        --title "SHAP8 12M (gauge +1.5 m, W10x5): $m x slope multipliers, cycle-1 fit vs after 3 cycles, and D2" \
        --fit "cycle-1 fit, $(awk '/^final_data_cost/{printf "%.3e", $2}' $CH/c1/fit/final/fit_result.txt)=$TB/${m}_c1fit.dat:dashed" \
        --fit "after 3 cycles, $(cost $CH/c3/ft/final)=$TB/${m}_c3ft.dat:solid" \
        --fit "clipped refit D2, 4.886e-03=$TB/clip_D2.dat:dotted" || exit 3
done
# overview: every activity-model family fitted on the lategrid base 2026-11-08 ... 11-11, best result of each
SM=$RUN/slope_mult/fit/final; BB=$RUN/slope_bumps/chain_slopemul_20261110T201212Z/c3/ft/final
$PY $S/positive_table.py "$SM" "$TB/slope_mult_harmonics.dat" || exit 1
$PY $S/positive_table.py "$BB" "$TB/slope_bumps_c3ft.dat" || exit 1
$PY $S/model_obs_figures.py --bundle $RUN/bundle/restart_templates.npz --overwrite --out "$O/12M_all_families" \
    --ref "$TB/clip_D2.dat" \
    --title "SHAP8 12M (gauge +1.5 m, W10x5): best fit of every activity-model family, 2026-11-08 ... 11-11" \
    --fit "softplus x slope multipliers, 3 cycles, $(cost $SP/c3/ft/final)=$TB/softplus_c3ft.dat:solid" \
    --fit "clipped refit D2, 4.89e-03=$TB/clip_D2.dat:dotted" \
    --fit "squared harmonics x slope multipliers, 3 cycles, $(cost $SQ/c3/ft/final)=$TB/squared_c3ft.dat:dashed" \
    --fit "positive, ratio 1e-5, lambda 1e-4, 3 cycles, 8.54e-03=$TB/positive_r1e-5_lam1e-4_c3.dat:dashdot" \
    --fit "slope multipliers x harmonics, floor 1e-5, 8.72e-03=$TB/slope_mult_harmonics.dat:dotted" \
    --fit "bumps x slope multipliers, 3 cycles, 1.41e-02=$TB/slope_bumps_c3ft.dat:dashdot" || exit 3
for m in softplus squared; do
    CH=$RUN/nonneg/$m/$TAG
    cp -p "$CH/c3/fit/final/nonneg_bins.txt" "$O/${m}_c3_fit_nonneg_bins.txt"
    cp -p "$CH/c3/fit/final/NONNEG_MODEL.txt" "$O/${m}_NONNEG_MODEL.txt"
done
{
echo "# stage  data_cost  rms_ln_NG  rms_ln_RG  stop/evaluations   (one job per model; chains under $RUN/nonneg/<model>/$TAG)"
for m in softplus squared; do
    CH=$RUN/nonneg/$m/$TAG; echo "## $m (job $(ls -d $CH/provenance-* | head -1 | sed 's/.*provenance-//'))"
    for c in c1 c2 c3; do
        f=$CH/$c/fit/final/fit_result.txt
        printf "%-6s %s %s %s %s\n" "${c}fit" $(awk '/^final_data_cost/{print $2} /^rms_ln_residual_NG/{print $2} /^rms_ln_residual_RG/{print $2}' $f) \
            "$(grep -oE 'ftol|maximum number' $f | head -1 | sed 's/maximum number/max_nfev/')/$(grep -oE 'nfev [0-9]+' $f | head -1 | cut -d' ' -f2)"
        g=$CH/$c/ft/final/finetune_result.txt
        printf "%-6s %s %s %s %s\n" "${c}ft" $(awk '/^cost_data/{print $5} /^rms_ln_NG/{print $5} /^rms_ln_RG/{print $5}' $g) \
            "$(awk '/^stop /{print $2}' $g)/$(awk '/^nfev_total/{print $2}' $g)"
    done
done
} > "$O/chain_stages.txt"
{
echo "# model  data_cost  rms_ln_NG  rms_ln_RG  stop   kappa(last fit)"
for m in softplus squared; do CH=$RUN/nonneg/$m/$TAG; g=$CH/c3/ft/final/finetune_result.txt
    printf "%-28s %s %s %s %s %s\n" "${m}_3cycles" $(awk '/^cost_data/{print $5} /^rms_ln_NG/{print $5} /^rms_ln_RG/{print $5}' $g) "$(awk '/^stop /{print $2}' $g)" "$(awk '/^kappa /{print $3; exit}' $CH/c3/fit/final/fit_result.txt)"; done
f=$D2/fit_result.txt
printf "%-28s %s %s %s %s %s\n" "clip_refit_D2" $(awk '/^final_data_cost/{print $2} /^rms_ln_residual_NG/{print $2} /^rms_ln_residual_RG/{print $2}' $f) max_nfev "$(awk '/^kappa /{print $3; exit}' $f)"
printf "%-28s %s %s %s %s %s\n" "positive_r1e-5_lam1e-4_c3" $(awk '/^cost_data/{print $5} /^rms_ln_NG/{print $5} /^rms_ln_RG/{print $5}' $PB/finetune_result.txt) "$(awk '/^stop /{print $2}' $PB/finetune_result.txt)" "$(awk '/^kappa /{print $3; exit}' $RUN/positive_r1e-5/lam1e-4/c3/fit/final/fit_result.txt)"
} > "$O/summary.txt"
# code: the two frozen snapshots, the job scripts, the figure scripts
for z in "$NX" "$CX"; do mkdir -p "$C/model_ext_frozen/$(basename $z)"; cp -rp "$z/." "$C/model_ext_frozen/$(basename $z)/"; done
chmod -R u+w "$C/model_ext_frozen"; find "$C" -name __pycache__ -prune -exec rm -rf {} + 2>/dev/null
mkdir -p "$C/jobs" "$C/figures"
cp -p "$SP/provenance-124162/chain_onejob.sbatch" "$C/jobs/"
cp -p $J/chain_onejob.sh $J/nonneg_start.sbatch $J/nonneg_fit.sbatch $J/nonneg_finetune.sbatch "$C/jobs/"
cp -p $S/figures_12m_nonneg.sh $S/positive_table.py $S/model_obs_figures.py $S/rerun_commands.py "$C/figures/"
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
# How to rerun: 12M softplus and squared-harmonics models x slope multipliers (3 cycles each, one job per model)

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_12m_nonneg.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never delete
earlier attempts (move them aside as \`*.incomplete_<UTC>\`). AMPS commit: \`$(head -1 $SP/provenance-124162/git.txt)\`.

**Models** (user decisions 2026-11-10/11): H2O activity, 7 effective-gravitational-slope bins with multipliers m_b >= 0
([0,40) held at 1); D2's six ft templates ft_c1..ft_c6 held at 1; no cliff_ge_60 template; no floor.
softplus: A = m_b s ln(1 + exp(L/s)), L = Y a (25 signed degree-<=4 harmonics) + A_s sum ft, s = 1e-5 A_orig per facet,
max(0, L) where A_orig = 0; squared: A = m_b (A_s (Y b)^2 + A_s sum ft). Start of cycle 1 from D2 (softplus: D2's
harmonic coefficients; squared: b by area-weighted least squares to D2's activity, CPU job 124138), multipliers 1, kappa
and grain from D2. Fit: lategrid base, W10x5 from 10:33:10, original T law held, kappa free, grain + crossing 10,
lambda 0, max_nfev 4000; fine tuning per facet exp(delta) on the facets seen only in the last 5 min, no lower bound,
MAXNFEV 500, grain refitted; 3 cycles; one job per model (records bound once), L step with 32 chunked threads.
Stage table: \`chain_stages.txt\`; per-bin multipliers, zero area, production: \`<model>_c3_fit_nonneg_bins.txt\`.

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_12m_nonneg.sh
\`\`\`

## A1. The squared start (CPU job 124138)

\`\`\`bash
sbatch --export=ALL,MODEL_EXT=$NX $J/nonneg_start.sbatch
\`\`\`

## A2. The two chains as submitted (one job each; attempts 2 and 3 continue after a time-limit stop)

\`\`\`bash
cd $J
./chain_onejob.sh --kind nonneg --model softplus --chain-root $SP --model-ext $NX --chain-code $CX --cycles 3 --nfev 4000 --check 1 --ft-maxnfev 500 --ft-facets late_only --ft-window-min 5 --l-threads 32 --node gpu001 --time 08:00:00 --attempts 3 --submit
./chain_onejob.sh --kind nonneg --model squared --chain-root $SQ --model-ext $NX --chain-code $CX --start-file $RUN/nonneg/bundle/squared_start.npz --cycles 3 --nfev 4000 --check 1 --ft-maxnfev 500 --ft-facets late_only --ft-window-min 5 --l-threads 32 --node gpu002 --time 08:00:00 --attempts 3 --submit
\`\`\`
Each stage also has its own \`commands.sh\` (runs that stage alone): \`<chain>/c<k>/{fit,ft}/provenance-<job>-<stage>/commands.sh\`; the job's records: \`<chain>/provenance-<job>/\`.

## A3. Commands of every stage, as run
MD
for m in softplus squared; do CH=$RUN/nonneg/$m/$TAG
  for s in c1/fit c1/ft c2/fit c2/ft c3/fit c3/ft; do p=$(ls -d $CH/$s/provenance-* | tail -1)
    printf '\n### %s %s (%s)\n\n```bash\n%s\n```\n' "$m" "$s" "$(basename $p)" "$(cat $p/commands.sh)"
  done
done
cat <<MD

## Code that ran (\`code/\`, with \`SHA256SUMS\`) and the git commit with an identical copy

| file | SHA256 (first 16) | commit |
|---|---|---|
$CM
MD
} > "$O/HOW-TO-RERUN.md"
echo "wrote $O"; cat "$O/chain_stages.txt"; cat "$O/summary.txt"
