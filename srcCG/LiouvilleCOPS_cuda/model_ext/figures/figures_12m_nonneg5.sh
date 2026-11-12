#!/bin/bash
# ==== version-controlled copy: header added 2026-11-12, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_12m_nonneg5.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_12m_nonneg5.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: db20975c16bf36a736465ef612d3d02d7c30268bb537aa8d399c96650a5f4d31
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure script. The 5-cycle softplus and squared chains (124310, 124313; user 2026-11-11 "repeat softplus and squared runs adding 2
#   extra cycles", 2026-11-12 "prepare comparison plots") against the COPS observations next to D2 and the best
#   positive fit; per model the cycle-1 fit and the state after each cycle; chain_stages.txt, summary.txt, bin
#   tables, MAPS.txt (where the maps are), HOW-TO-RERUN.md (chain lines, map commands, every stage), code/.
# RUN: figures_12m_nonneg5.sh   (login node, no GPU) -> figures_cluster/12M/model_vs_observations_nonneg_5cycles/
# ==== end of header ====
# The 5-cycle softplus and squared-harmonics chains (user 2026-11-11: "repeat softplus and squared runs adding 2 extra
# cycles of fitting and fine tuning; save the activity and source rate maps for visualization"; answers: fresh 5-cycle
# runs, maps of the final state) against the COPS observations (user 2026-11-12: "prepare comparison plots"); lategrid
# base (gauge +1.5 m, W10x5). Same models and settings as figures_12m_nonneg.sh (3 cycles), cycles 1-3 bitwise equal.
# Writes figures_cluster/12M/model_vs_observations_nonneg_5cycles/:
#   12M_nonneg5_final_*        softplus and squared after 5 cycles vs D2 and the best positive fit
#   12M_nonneg5_<model>_cycles_* the cycle-1 fit and the state after each of the 5 cycles (fine tunings)
#   chain_stages.txt, summary.txt, the cycle-5 bin tables, MAPS.txt (where the activity / source-rate maps are),
#   HOW-TO-RERUN.md, code/.                Login node, no GPU:   figures_12m_nonneg5.sh
set -u
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; J=$COPS_OUTPUT_ROOT/jobs
RUN=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; F12=$COPS_OUTPUT_ROOT/figures_cluster/12M
TAG=chain_5cycles_20261112T012553Z; NC=5
declare -A JOB=([softplus]=124310 [squared]=124313) MAPJ=([softplus]=124324 [squared]=124325)
SP=$RUN/nonneg/softplus/$TAG; SQ=$RUN/nonneg/squared/$TAG
D2=$RUN/clip_refit/D2/final; PB=$RUN/positive_r1e-5/lam1e-4/c3/ft
NX=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_20261111T044920Z
CX=$COPS_OUTPUT_ROOT/model_ext_frozen/chain_onejob_20261111T060742Z
MX=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_maps_20261112T021623Z
O=$F12/model_vs_observations_nonneg_5cycles; TB=$O/tables; C=$O/code
mkdir -p "$TB" "$C"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
PY=".venv/bin/python -W ignore"
cost() { awk '/^cost_data/{printf "%.4e", $5}' "$1/finetune_result.txt"; }
fcost() { awk '/^final_data_cost/{printf "%.4e", $2}' "$1/fit_result.txt"; }
$PY $S/positive_table.py "$D2" "$TB/clip_D2.dat" || exit 1
$PY $S/positive_table.py "$PB" "$TB/positive_r1e-5_lam1e-4_c3.dat" || exit 1
for m in softplus squared; do CH=$RUN/nonneg/$m/$TAG
    for c in $(seq 1 $NC); do
        $PY $S/positive_table.py "$CH/c$c/fit/final" "$TB/${m}_c${c}fit.dat" || exit 1
        $PY $S/positive_table.py "$CH/c$c/ft/final" "$TB/${m}_c${c}ft.dat" || exit 1
    done
done
$PY $S/model_obs_figures.py --bundle $RUN/bundle/restart_templates.npz --overwrite --out "$O/12M_nonneg5_final" \
    --ref "$TB/clip_D2.dat" \
    --title "SHAP8 12M (gauge +1.5 m, W10x5): softplus and squared harmonics x slope multipliers after 5 cycles, vs D2" \
    --fit "softplus x slope multipliers, 5 cycles, $(cost $SP/c5/ft/final)=$TB/softplus_c5ft.dat:solid" \
    --fit "squared harmonics x slope multipliers, 5 cycles, $(cost $SQ/c5/ft/final)=$TB/squared_c5ft.dat:dashed" \
    --fit "clipped refit D2, 4.8864e-03=$TB/clip_D2.dat:dotted" \
    --fit "positive, ratio 1e-5, lambda 1e-4, 3 cycles, 8.5367e-03=$TB/positive_r1e-5_lam1e-4_c3.dat:dashdot" || exit 3
# per model: the cycle-1 fit and the state after each of the 5 cycles (6 curves: the figure palette has 8 colours, so
# not all 10 stages fit in one plot; every stage's numbers are in chain_stages.txt)
for m in softplus squared; do CH=$RUN/nonneg/$m/$TAG
    STY=(dashdot dashed dashdot dashed solid)
    CYC=(--fit "cycle-1 fit, $(fcost $CH/c1/fit/final)=$TB/${m}_c1fit.dat:dotted")
    for c in $(seq 1 $NC); do
        CYC+=(--fit "after cycle $c, $(cost $CH/c$c/ft/final)=$TB/${m}_c${c}ft.dat:${STY[$((c - 1))]}")
    done
    $PY $S/model_obs_figures.py --bundle $RUN/bundle/restart_templates.npz --overwrite --out "$O/12M_nonneg5_${m}_cycles" \
        --ref "$TB/${m}_c5ft.dat" \
        --title "SHAP8 12M (gauge +1.5 m, W10x5): $m x slope multipliers, cycle-1 fit and the state after each of 5 cycles" \
        "${CYC[@]}" || exit 3
    cp -p "$CH/c5/fit/final/nonneg_bins.txt" "$O/${m}_c5_fit_nonneg_bins.txt"
    cp -p "$CH/c5/fit/final/NONNEG_MODEL.txt" "$O/${m}_NONNEG_MODEL.txt"
done
{
echo "# stage  data_cost  rms_ln_NG  rms_ln_RG  stop/evaluations   (one job per model; chains under $RUN/nonneg/<model>/$TAG)"
for m in softplus squared; do CH=$RUN/nonneg/$m/$TAG; echo "## $m (job ${JOB[$m]})"
    for c in $(seq 1 $NC); do
        f=$CH/c$c/fit/final/fit_result.txt
        printf "%-6s %s %s %s %s\n" "c${c}fit" $(awk '/^final_data_cost/{print $2} /^rms_ln_residual_NG/{print $2} /^rms_ln_residual_RG/{print $2}' $f) \
            "$(grep -oE 'ftol|maximum number' $f | head -1 | sed 's/maximum number/max_nfev/')/$(grep -oE 'nfev [0-9]+' $f | head -1 | cut -d' ' -f2)"
        g=$CH/c$c/ft/final/finetune_result.txt
        printf "%-6s %s %s %s %s\n" "c${c}ft" $(awk '/^cost_data/{print $5} /^rms_ln_NG/{print $5} /^rms_ln_RG/{print $5}' $g) \
            "$(awk '/^stop /{print $2}' $g)/$(awk '/^nfev_total/{print $2}' $g)"
    done
done
} > "$O/chain_stages.txt"
{
echo "# model  data_cost  rms_ln_NG  rms_ln_RG  stop  kappa(last fit)"
for m in softplus squared; do CH=$RUN/nonneg/$m/$TAG; g=$CH/c5/ft/final/finetune_result.txt
    printf "%-28s %s %s %s %s %s\n" "${m}_5cycles" $(awk '/^cost_data/{print $5} /^rms_ln_NG/{print $5} /^rms_ln_RG/{print $5}' $g) "$(awk '/^stop /{print $2}' $g)" "$(awk '/^kappa /{print $3; exit}' $CH/c5/fit/final/fit_result.txt)"; done
f=$D2/fit_result.txt
printf "%-28s %s %s %s %s %s\n" "clip_refit_D2" $(awk '/^final_data_cost/{print $2} /^rms_ln_residual_NG/{print $2} /^rms_ln_residual_RG/{print $2}' $f) max_nfev "$(awk '/^kappa /{print $3; exit}' $f)"
printf "%-28s %s %s %s %s %s\n" "positive_r1e-5_lam1e-4_c3" $(awk '/^cost_data/{print $5} /^rms_ln_NG/{print $5} /^rms_ln_RG/{print $5}' $PB/finetune_result.txt) "$(awk '/^stop /{print $2}' $PB/finetune_result.txt)" "$(awk '/^kappa /{print $3; exit}' $RUN/positive_r1e-5/lam1e-4/c3/fit/final/fit_result.txt)"
} > "$O/summary.txt"
{
echo "Activity and instantaneous source-rate maps of the final state (after cycle 5), written by model_ext/nonneg_maps.py"
echo "(frozen code $MX; job script $J/nonneg_maps.sbatch). Each folder has README.txt (fields, units, law) and maps_manifest.json."
for m in softplus squared; do echo "  $m: $RUN/nonneg/$m/$TAG/maps   (job ${MAPJ[$m]})"; done
echo "ParaView: <model>_flux.pvd (13 times) + <model>_static.vtu; Tecplot 360: <model>_maps.dat."
} > "$O/MAPS.txt"
for z in "$NX" "$CX" "$MX"; do mkdir -p "$C/model_ext_frozen/$(basename $z)"; cp -rp "$z/." "$C/model_ext_frozen/$(basename $z)/"; done
chmod -R u+w "$C/model_ext_frozen"; find "$C" -name __pycache__ -prune -exec rm -rf {} + 2>/dev/null
mkdir -p "$C/jobs" "$C/figures"
cp -p "$SP/provenance-${JOB[softplus]}/chain_onejob.sbatch" "$C/jobs/"
cp -p $J/chain_onejob.sh $J/nonneg_start.sbatch $J/nonneg_fit.sbatch $J/nonneg_finetune.sbatch $J/nonneg_maps.sbatch "$C/jobs/"
cp -p $S/figures_12m_nonneg5.sh $S/positive_table.py $S/model_obs_figures.py $S/rerun_commands.py "$C/figures/"
(cd "$C" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
CM=$($PY - "$C" <<'PY'
import hashlib, sys
from pathlib import Path
sys.path.insert(0, "/data001/highresolu_vt/cops_output/figures_cluster/scripts")
from rerun_commands import commit_of
for p in sorted(Path(sys.argv[1]).rglob("*")):
    if p.is_file() and p.suffix in (".py", ".sbatch", ".sh", ".cu", ".cpp"):
        h = hashlib.sha256(p.read_bytes()).hexdigest(); c = commit_of(p.name, h)
        print("| `%s` | `%s` | %s |" % (p.relative_to(sys.argv[1]), h[:16], "`%s`" % c if c else "**not in git**"))
PY
)
{
cat <<MD
# How to rerun: 12M softplus and squared-harmonics models x slope multipliers, 5 cycles each (one job per model)

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_12m_nonneg5.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never delete
earlier attempts (move them aside as \`*.incomplete_<UTC>\`). AMPS commit: \`$(head -1 $SP/provenance-${JOB[softplus]}/git.txt)\`.

**Models and settings:** as in \`model_vs_observations_nonneg/HOW-TO-RERUN.md\` (3-cycle runs) and the README of the
version-controlled \`srcCG/LiouvilleCOPS_cuda/model_ext\` (section "Softplus and squared-harmonics models"), with
\`--cycles 5\`. Cycles 1-3 of these runs are bitwise equal to the 3-cycle runs (checked: every array of the fit
restarts, facet files and fine-tuning states). Maps of the final state: \`MAPS.txt\`.

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_12m_nonneg5.sh
\`\`\`

## A1. Inputs: lategrid base, D2, slope bins, squared start

See "Running them from scratch" in \`srcCG/LiouvilleCOPS_cuda/model_ext/README.md\` (steps 0-5); the squared start used
here is \`$RUN/nonneg/bundle/squared_start.npz\` (CPU job 124138).

## A2. The two chains as submitted (one job each, concurrently)

\`\`\`bash
cd $J
./chain_onejob.sh --kind nonneg --model softplus --chain-root $SP --model-ext $NX --chain-code $CX --cycles 5 --nfev 4000 --check 1 --ft-maxnfev 500 --ft-facets late_only --ft-window-min 5 --l-threads 32 --node gpu001 --time 08:00:00 --attempts 3 --submit
./chain_onejob.sh --kind nonneg --model squared --chain-root $SQ --model-ext $NX --chain-code $CX --start-file $RUN/nonneg/bundle/squared_start.npz --cycles 5 --nfev 4000 --check 1 --ft-maxnfev 500 --ft-facets late_only --ft-window-min 5 --l-threads 32 --node gpu002 --time 08:00:00 --attempts 3 --submit
\`\`\`

## A3. The maps (CPU job per model after its chain; \`$J/nonneg_maps.sbatch\`), commands as run

### softplus (job ${MAPJ[softplus]})

\`\`\`bash
$(cat $SP/maps/provenance-${MAPJ[softplus]}/commands.sh)
\`\`\`

### squared (job ${MAPJ[squared]})

\`\`\`bash
$(cat $SQ/maps/provenance-${MAPJ[squared]}/commands.sh)
\`\`\`

## A4. Commands of every stage, as run
MD
for m in softplus squared; do CH=$RUN/nonneg/$m/$TAG
  for c in $(seq 1 $NC); do for st in fit ft; do p=$(ls -d $CH/c$c/$st/provenance-* | tail -1)
    printf '\n### %s c%s/%s (%s)\n\n```bash\n%s\n```\n' "$m" "$c" "$st" "$(basename $p)" "$(cat $p/commands.sh)"
  done; done
done
cat <<MD

## Code that ran (\`code/\`, with \`SHA256SUMS\`) and the git commit with an identical copy

| file | SHA256 (first 16) | commit |
|---|---|---|
$CM
MD
} > "$O/HOW-TO-RERUN.md"
echo "wrote $O"; cat "$O/chain_stages.txt"; cat "$O/summary.txt"
