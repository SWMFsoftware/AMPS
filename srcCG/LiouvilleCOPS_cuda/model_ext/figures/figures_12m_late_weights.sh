#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_12m_late_weights.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_12m_late_weights.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 4e1acabbe88b065ac7be8133d38a20bcfd45e5c4b223124e39abda05c44c2d13
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison figures (model vs COPS observations; .eps + .png 300 dpi) + HOW-TO-RERUN.md for: 12M late-weights fit vs cycle 3 (also writes code/ with the code that ran).
# RUN: on the login node, no GPU: <this> [<job id>]
# OUT: $COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_late_weights/
# ==== end of header ====
# 12M fit from cycle 3 with the last 5 points weighted x10 (NG and RG; job 122236, jobs/grain_late_weights.sbatch,
# model_ext/grain_fit_weighted.py + late_weights.py; user 2026-11-06), against the COPS observations, with cycle 3;
# writes HOW-TO-RERUN.md (with the AMPS commit of each fit) and code/ (a copy of the code that ran; user 2026-11-06).
# Login node, no GPU:   figures_12m_late_weights.sh [<job id>]
set -u
JOB=${1:-122236}
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; J=$COPS_OUTPUT_ROOT/jobs; X=$COPS_OUTPUT_ROOT/model_ext
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; CY=$R/grain/cycles
N=C30_n4000_grain_cross10_c3_last5_w10_constant_ngrg; LW=$R/grain/late_weights/$N
O=$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_late_weights; TB=$O/tables; C=$O/code
mkdir -p "$TB" "$C"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
F3=$CY/c3/C30_n4000_grain_cross10_ftc3; FT3=$CY/c3/finetune/C30_n4000_grain_cross10_ftc3_m5_i5_cpu
FW=$LW/fit; FTW=$LW/finetune/${N}_m5_i5_cpu
.venv/bin/python $S/finetune_table.py "$FT3" "$F3" "$TB/c3_finetuned.dat" || exit 1
.venv/bin/python $S/finetune_table.py "$FTW" "$FW" "$TB/last5w10_finetuned.dat" || exit 1
T='SHAP8 12M, cliff $\geq$60$^\circ$ template, free T, dust grain in NG + crossing constraint'
.venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R/bundle/restart_templates.npz --overwrite \
    --out "$O/12M_last5_w10" --ref "$FW" --title "$T: last 5 points weighted x10 (NG and RG)" \
    --fit "cycle 3 fit=$F3:dotted" --fit "cycle 3, fine-tuned=$TB/c3_finetuned.dat:dashdot" \
    --fit "last 5 points x10, fit=$FW:dashed" --fit "last 5 points x10, fine-tuned=$TB/last5w10_finetuned.dat:solid" \
    && echo "wrote $O/12M_last5_w10_*"
# code copy: the job script and model_ext modules as they ran (provenance copies), the figure scripts
P=$LW/provenance-$JOB; P3=$(ls -d $CY/c3/provenance-* | head -1)
mkdir -p "$C/late_weights_fit-$JOB" "$C/cycle3_fit" "$C/figures"
cp -p "$P"/*.sbatch "$P"/*.py "$P"/model_ext.sha256 "$P"/git.txt "$C/late_weights_fit-$JOB/"
cp -p "$P3"/*.sbatch "$P3"/*.py "$P3"/model_ext.sha256 "$P3"/git.txt "$C/cycle3_fit/"
cp -p $S/figures_12m_late_weights.sh $S/model_obs_figures.py $S/finetune_table.py "$C/figures/"
(cd "$C" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
commit() { head -1 "$1/git.txt"; }; dirty() { tail -n +2 "$1/git.txt" | grep -v '^??' | tr '\n' ' '; }
cat > "$O/HOW-TO-RERUN.md" <<MD
# How to rerun: 12M fit from cycle 3 with the last 5 points weighted x10

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_12m_late_weights.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never
delete earlier attempts (move them aside as \`*.incomplete_<UTC>\`).

## Code version

| Plotted case | Job | AMPS commit | Tracked local changes |
|---|---|---|---|
| last 5 points x10, fit + fine tuning | $JOB | \`$(commit $P)\` | $(dirty $P) |
| cycle 3, fit + fine tuning | $(basename $P3 | sed 's/provenance-//') | \`$(commit $P3)\` | $(dirty $P3) |

Untracked code that ran (\`$X\`, outside git) is copied in \`code/\` (job scripts and model_ext modules from each
job's provenance, with their \`model_ext.sha256\` and \`git.txt\`; the figure scripts; \`code/SHA256SUMS\`).

## What was done

\`jobs/grain_late_weights.sbatch\` (job $JOB, gpu002, 3 GPUs; provenance \`$P/\`): \`model_ext/grain_fit_weighted.py\`
= the cycle-3 C30 command (grain + crossing 10, mixed, --tbounds 30 400, degree 4/4, lambda 0, cliff_ge_60deg +
ft_c1..ft_c3 held at 1, \`--restart-from\` cycle 3's restart with its fine-tuned grain, --template-start 0,
max_nfev 4000) with \`--weight-last 5 --weight-factor 10 --weight-profile constant --weight-gauges ng,rg\`: the cost
terms of the last 5 used NG and RG points (10:39:10-10:39:26 UTC, 21-7 m) x10 (residuals x sqrt 10). Reported costs
and the restart file are unweighted; the weighted objective is in \`fit/fit_manifest.json\` ("late_weights").
Then \`grain_finetune.py\` (5 min, 5 intervals, cpu). Earlier attempts: 122226 (device list; test dir
\`late_weights_test-122226.incomplete_*\`), 122229 (weighted NG only: sbatch split WGAUGES at the comma;
\`${N%ngrg}ng.incomplete_*\`).

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_12m_late_weights.sh $JOB
\`\`\`

## A1. The fit and fine tuning (from cycle 3's restart files)

Move \`$LW\` aside first (as \`$LW.incomplete_<UTC>\`), then (both gauges are the default; to change them,
\`export WGAUGES=ng\` in the shell, since sbatch splits --export values at commas):

\`\`\`bash
sbatch --gpus-per-node=3 --cpus-per-task=16 --mem=80G --time=06:00:00 \\
    --export=ALL,WLAST=5,WFACTOR=10,WPROFILE=constant $J/grain_late_weights.sbatch
\`\`\`

(the script's defaults are 6 GPUs, 32 CPUs, 150 GB on gpu002; 3 GPUs were used because only 3 were free.)

## B. Cycle 3 and everything before it (from scratch)

See \`$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_ft_cycles/HOW-TO-RERUN.md\`, then A1, A0.
MD
echo "wrote $O/HOW-TO-RERUN.md"
