#!/bin/bash
# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_12m_clip_vs_positive.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_12m_clip_vs_positive.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 9e76f2cf350ffdac534a7379295b0e18944a76a7004438610a62424077d583f1
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure script. The last fitting round (ratio-1e-5 positive chains, after 3 cycles) vs the diagnostic refit of the
#   clipped original-T model on the same base (D2), and the diagnostic refits D0/D1/D2 against the COPS observations;
#   summary.txt, HOW-TO-RERUN.md, code/. 2026-11-09.
# RUN: figures_12m_clip_vs_positive.sh   (login node, no GPU)
# ==== end of header ====
# The last fitting round against the COPS observations (user 2026-11-09: "prepare comparison plot from the last fitting
# round"): the ratio-1e-5 positive-activity chains (final state after 3 cycles) next to the diagnostic refit of the
# earlier clipped model on the same base (user 2026-11-08: "refit the earlier clipped model (original temperature, with its
# templates) on the new base, with W10x5 and the +1.5 m offset"; control D0 "yes"). Writes into
# figures_cluster/12M/model_vs_observations_clip_vs_positive/:
#   12M_clip_vs_positive_*   D2 (clipped, lategrid base, W10x5) vs ratio-1e-5 positive finals for lambda 1e-4, 1e-3, 1e-2
#   12M_clip_diagnostic_*    D0 (c6 continued on its own base), D1 (lategrid base, +1.5 m), D2 (+ W10x5)
#   summary.txt, HOW-TO-RERUN.md, code/
# All tables have the 1150 observation rows (grid rows removed by positive_table.py). D0 ran on the old base (gauge at
# the trajectory point): its rows are the same 1150 observation times and data; on the altitude axis (the lategrid
# base's, gauge +1.5 m) it is drawn up to 1.5 m off its own altitude, the time axis is exact.
# Login node, no GPU:   figures_12m_clip_vs_positive.sh
set -u
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; J=$COPS_OUTPUT_ROOT/jobs
RUN=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; E=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier
CL=$RUN/clip_refit; CY=$RUN/positive_r1e-5
O=$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_clip_vs_positive; TB=$O/tables; C=$O/code
mkdir -p "$TB" "$C"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
PY=".venv/bin/python -W ignore"
for v in D1 D2; do $PY $S/positive_table.py "$CL/$v/final" "$TB/clip_$v.dat" || exit 1; done
cp -p "$E/clip_refit/D0/final/COPS_model.dat" "$TB/clip_D0.dat"
for L in 1e-4 1e-3 1e-2; do $PY $S/positive_table.py "$CY/lam$L/c3/ft" "$TB/pos_r1e-5_lam${L}_final.dat" || exit 1; done
B=$RUN/bundle/restart_templates.npz
$PY $S/model_obs_figures.py --bundle $B --overwrite --out "$O/12M_clip_vs_positive" --ref "$TB/clip_D2.dat" \
    --title "SHAP8 12M, same base (gauge +1.5 m, W10x5): clipped original-T model refit vs positive activity (ratio 1e-5, after 3 cycles)" \
    --fit "clipped model refit (D2), zero activity allowed=$TB/clip_D2.dat:solid" \
    --fit "positive, lambda 1e-4=$TB/pos_r1e-5_lam1e-4_final.dat:dashed" \
    --fit "positive, lambda 1e-3=$TB/pos_r1e-5_lam1e-3_final.dat:dashdot" \
    --fit "positive, lambda 1e-2=$TB/pos_r1e-5_lam1e-2_final.dat:dotted" || exit 3
$PY $S/model_obs_figures.py --bundle $B --overwrite --out "$O/12M_clip_diagnostic" --ref "$TB/clip_D2.dat" \
    --title "SHAP8 12M, clipped original-T cycle-6 model refit: old base (D0), new base +1.5 m (D1), + W10x5 (D2)" \
    --fit "D0, old base (control)=$TB/clip_D0.dat:dotted" \
    --fit "D1, new base, gauge +1.5 m=$TB/clip_D1.dat:dashed" \
    --fit "D2, new base, gauge +1.5 m, W10x5=$TB/clip_D2.dat:solid" || exit 3
# summary: costs and rms of every plotted model (fit_result.txt / finetune_result.txt)
{
echo "# model  total_cost  data_cost  rms_ln_NG  rms_ln_RG  stop   (fitter units; lambda 0 for D0-D2)"
for v in D0 D1 D2; do d=$CL/$v/final; [ $v = D0 ] && d=$E/clip_refit/D0/final
    printf "%-26s %s %s %s %s %s\n" "clip_$v" $(awk '/^final_cost/{print $2} /^final_data_cost/{print $2} /^rms_ln_residual_NG/{print $2} /^rms_ln_residual_RG/{print $2}' $d/fit_result.txt) \
        "$(grep -oE 'ftol|maximum number' $d/fit_result.txt | head -1 | sed 's/maximum number/max_nfev/')"
done
for L in 1e-4 1e-3 1e-2; do f=$CY/lam$L/c3/ft/finetune_result.txt
    printf "%-26s %s %s %s %s %s\n" "positive_r1e-5_lam$L" $(awk '/^cost_total/{print $5} /^cost_data/{print $5} /^rms_ln_NG/{print $5} /^rms_ln_RG/{print $5}' $f) \
        "$(awk '/^stop /{print $2}' $f)"
done
} > "$O/summary.txt"
# code: snapshots, job scripts as run, figure scripts
for snap in clip_20261108T234013Z r1e-5k_20261109T023540Z; do
    mkdir -p "$C/model_ext_frozen/$snap"; cp -rp "$COPS_OUTPUT_ROOT/model_ext_frozen/$snap/." "$C/model_ext_frozen/$snap/"
done
chmod -R u+w "$C/model_ext_frozen"; mkdir -p "$C/jobs" "$C/figures"
cp -p $(ls -d $CL/D2/provenance-* | tail -1)/clip_refit.sbatch "$C/jobs/clip_refit.sbatch"
cp -p $S/figures_12m_clip_vs_positive.sh $S/positive_table.py $S/model_obs_figures.py $S/rerun_commands.py "$C/figures/"
find "$C" -name __pycache__ -prune -exec rm -rf {} + 2>/dev/null
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
# How to rerun: last fitting round vs the clipped-model refit (12M)

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_12m_clip_vs_positive.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never delete
earlier attempts (move them aside as \`*.incomplete_<UTC>\`).

**What is compared** (numbers: \`summary.txt\`):
- D0, D1, D2: the original-temperature cycle-6 clipped model (A = max(0, Y a + c_cliff T_cliff + sum ft_ck), ft_c1..ft_c6 held
  at 1, Tmax/T75/Tmin 182.09922/136.11149/133 K held, kappa free, grain + crossing 10, lambda 0, srcCG night ratios, no floor)
  refit from c6's final state with up to 4000 evaluations: D0 on its own base \`$E/base\` (control), D1 on
  \`$RUN/base\` (gauge +1.5 m along +Z, model at the observation time stamps), D2 = D1 + W10x5 for the observations from
  10:33:10 UTC.
- Positive activity, ratio 1e-5 (night/shadow both species, H2O and CO2 floor), after 3 cycles of fit + late-facet fine
  tuning, lambda 1e-4, 1e-3, 1e-2 (same base and W10x5 as D2): every command in
  \`$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_positive_cycles_r1e-5/HOW-TO-RERUN.md\`.

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_12m_clip_vs_positive.sh
\`\`\`

## A1. The diagnostic refits D0, D1, D2 (each on one GPU node; inside a Slurm job), commands as run

MD
for v in D0 D1 D2; do d=$CL/$v; [ $v = D0 ] && d=$E/clip_refit/D0
    p=$(ls -d $d/provenance-* | tail -1)
    echo "**$v ($(basename $p)):** AMPS \`$(head -1 $p/git.txt)\`"; echo; echo '```bash'; cat "$p/commands.sh"; echo '```'; echo
done
cat <<MD
Submission (2026-11-09T00:02Z and 00:17Z): \`cd $J; F=$COPS_OUTPUT_ROOT/model_ext_frozen/clip_20261108T234013Z;
sbatch --nodelist=gpu001 --export=ALL,VARIANT=D1,MODEL_EXT=\$F clip_refit.sbatch; sbatch --nodelist=gpu002
--export=ALL,VARIANT=D2,MODEL_EXT=\$F clip_refit.sbatch; sbatch --nodelist=gpu001 --dependency=afterany:<D1 job>
--export=ALL,VARIANT=D0,MODEL_EXT=\$F clip_refit.sbatch\`.

## Code that ran (\`code/\`, with \`SHA256SUMS\`) and the git commit with an identical copy

| file | SHA256 (first 16) | commit |
|---|---|---|
$CM
MD
} > "$O/HOW-TO-RERUN.md"
echo "wrote $O"; cat "$O/summary.txt"
