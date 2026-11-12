#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_12m_slope_mult.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_12m_slope_mult.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: f254ffb6059e666b9b1a6bd718c0650ac5339b09d8be0621ecda46fc6910c077
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure script. The slope-multiplier fit (job 123786; user decisions 2026-11-10: bins 0-40 ... 80-90 + 90-180,
#   floor 1e-5 A_orig, srcCG start) against the COPS observations next to D2 and the best positive fit;
#   slope_mult_bins.txt, summary.txt, HOW-TO-RERUN.md, code/. 2026-11-10.
# RUN: figures_12m_slope_mult.sh   (login node, no GPU) -> figures_cluster/12M/model_vs_observations_slope_mult/
# ==== end of header ====
# The slope-multiplier fit against the COPS observations (user 2026-11-10: "define multiple cliff templates (0-40, 40-50,
# 50-60, 60-70, 70-80, and 80-90) ... activity multiplier ... multiplied by activity from spherical harmonics ...; each of
# these multipliers will be fit independently; ... clip the activity not by 0, but by 1E-5"; answers: floor 1e-5 x A_orig,
# 7th bin 90-180, D2 setup without templates, multipliers >= 0, start srcCG), next to the clipped-model refit D2 and the
# best positive-activity fit (ratio 1e-5, lambda 1e-4, after 3 cycles), all on the lategrid base (gauge +1.5 m, W10x5).
# Writes figures_cluster/12M/model_vs_observations_slope_mult/: 12M_slope_mult_*, slope_mult_bins.txt (copy),
# summary.txt, HOW-TO-RERUN.md, code/.     Login node, no GPU:   figures_12m_slope_mult.sh
set -u
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; J=$COPS_OUTPUT_ROOT/jobs
RUN=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; SM=$RUN/slope_mult/fit/final; D2=$RUN/clip_refit/D2/final
PB=$RUN/positive_r1e-5/lam1e-4/c3/ft
O=$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_slope_mult; TB=$O/tables; C=$O/code
mkdir -p "$TB" "$C"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
PY=".venv/bin/python -W ignore"
$PY $S/positive_table.py "$SM" "$TB/slope_mult.dat" || exit 1
$PY $S/positive_table.py "$D2" "$TB/clip_D2.dat" || exit 1
$PY $S/positive_table.py "$PB" "$TB/positive_r1e-5_lam1e-4.dat" || exit 1
$PY $S/model_obs_figures.py --bundle $RUN/bundle/restart_templates.npz --overwrite --out "$O/12M_slope_mult" \
    --ref "$TB/slope_mult.dat" \
    --title "SHAP8 12M (gauge +1.5 m, W10x5): slope-binned activity multipliers, floor 1e-5 A_orig, vs clipped refit D2 and best positive fit" \
    --fit "slope multipliers, floor 1e-5 (srcCG start)=$TB/slope_mult.dat:solid" \
    --fit "clipped refit D2 (cliff + 6 templates)=$TB/clip_D2.dat:dashed" \
    --fit "positive, ratio 1e-5, lambda 1e-4, 3 cycles=$TB/positive_r1e-5_lam1e-4.dat:dotted" || exit 3
cp -p "$SM/slope_mult_bins.txt" "$SM/SLOPE_MULTIPLIER_MODEL.txt" "$O/"
{
echo "# model  total_cost  data_cost  rms_ln_NG  rms_ln_RG  stop"
for pair in "slope_mult:$SM/fit_result.txt" "clip_D2:$D2/fit_result.txt"; do n=${pair%%:*}; f=${pair#*:}
    printf "%-26s %s %s %s %s %s\n" "$n" $(awk '/^final_cost/{print $2} /^final_data_cost/{print $2} /^rms_ln_residual_NG/{print $2} /^rms_ln_residual_RG/{print $2}' $f) \
        "$(grep -oE 'ftol|maximum number' $f | head -1 | sed 's/maximum number/max_nfev/')"
done
printf "%-26s %s %s %s %s %s\n" "positive_r1e-5_lam1e-4" $(awk '/^cost_total/{print $5} /^cost_data/{print $5} /^rms_ln_NG/{print $5} /^rms_ln_RG/{print $5}' $PB/finetune_result.txt) "$(awk '/^stop /{print $2}' $PB/finetune_result.txt)"
} > "$O/summary.txt"
# code: snapshot, job script as run, figure scripts
P=$(ls -d $RUN/slope_mult/fit/provenance-* | tail -1)
mkdir -p "$C/model_ext_frozen/slopemul_20261110T094013Z" "$C/jobs" "$C/figures"
cp -rp "$COPS_OUTPUT_ROOT/model_ext_frozen/slopemul_20261110T094013Z/." "$C/model_ext_frozen/slopemul_20261110T094013Z/"
chmod -R u+w "$C/model_ext_frozen"; find "$C" -name __pycache__ -prune -exec rm -rf {} + 2>/dev/null
cp -p "$P/slope_mult_fit.sbatch" "$C/jobs/"
cp -p $S/figures_12m_slope_mult.sh $S/positive_table.py $S/model_obs_figures.py $S/rerun_commands.py "$C/figures/"
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
# How to rerun: 12M slope-multiplier fit (job $(basename $P | sed 's/provenance-//'))

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_12m_slope_mult.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never delete
earlier attempts (move them aside as \`*.incomplete_<UTC>\`).

**Model** (user decisions 2026-11-10): H2O activity A = max(1e-5 A_orig, m_b L), L = the 25 signed degree-4 harmonics,
one multiplier m_b >= 0 per effective-gravitational-slope bin [0,40) (held at 1), [40,50), [50,60), [60,70), [70,80), [80,90),
[90,180]; no cliff or fine-tuning templates; lategrid base (gauge +1.5 m), W10x5 from 10:33:10, original temperature law
held, kappa free, grain + crossing 10, lambda 0, srcCG night ratios, max_nfev 4000; start = the original srcCG harmonic
coefficients, multipliers 1, kappa 1, grain from D2. Details: \`SLOPE_MULTIPLIER_MODEL.txt\`, per bin: \`slope_mult_bins.txt\`.
Compared: clipped refit D2 (\`$D2\`) and the positive fit ratio 1e-5, lambda 1e-4 after 3 cycles (\`$PB\`; its commands:
\`model_vs_observations_positive_cycles_r1e-5/HOW-TO-RERUN.md\`, D2's: \`model_vs_observations_clip_vs_positive/HOW-TO-RERUN.md\`).

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_12m_slope_mult.sh
\`\`\`

## A1. The fit (one GPU node; inside a Slurm job), commands as run ($(basename $P); AMPS \`$(head -1 $P/git.txt)\`)

\`\`\`bash
$(cat $P/commands.sh)
\`\`\`

Submission: \`sbatch --export=ALL,MODEL_EXT=$COPS_OUTPUT_ROOT/model_ext_frozen/slopemul_20261110T094013Z,START=srccg $J/slope_mult_fit.sbatch\`.
The slope bins: \`$RUN/slope_mult/bundle/slope_bins.npz\` (model_ext/slope_bins.py, CPU job 123624).

## Code that ran (\`code/\`, with \`SHA256SUMS\`) and the git commit with an identical copy

| file | SHA256 (first 16) | commit |
|---|---|---|
$CM
MD
} > "$O/HOW-TO-RERUN.md"
echo "wrote $O"; cat "$O/summary.txt"
