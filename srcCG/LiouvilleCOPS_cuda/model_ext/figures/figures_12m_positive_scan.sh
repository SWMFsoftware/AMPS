#!/bin/bash
# ==== version-controlled copy: header added 2026-11-07, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_12m_positive_scan.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_12m_positive_scan.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 6356078a906446cb5b3d053d446b85128e1bb8bb2e29232de051fbad719e496d
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure script. The six 12M positive-activity lambda fits against the COPS observations, scan table and
#   L-curve, HOW-TO-RERUN.md (base preparation and every fit, all arguments; rerun_commands.py) and code/, into
#   figures_cluster/12M/model_vs_observations_positive_lambda_scan/. 2026-11-07.
# RUN: figures_12m_positive_scan.sh   (login node, no GPU)
# ==== end of header ====
# 12M positive-activity lambda scan against the COPS observations (user decisions 2026-11-07: A = exp(Y a + c T_cliff),
# c >= 0, night/shadow ratio 0.02 for H2O and CO2, H2O floor 0.02 q_max A_orig, start and prior = emission-weighted
# projection of ln A_orig, area-weighted log Tikhonov, lambda 0, 1e-5 ... 1e-1 "the user picks from the L-curve";
# original temperature law held, kappa free, grain, crossing weight 10, W10x5 for the observations from 10:33:10,
# gauge offset +1.5 m along +Z, max_nfev 4000). The fits' 401 model-only grid rows are removed (positive_table.py).
# Writes the figures, the scan table and L-curve (analysis/positive_lambda_scan_12m), HOW-TO-RERUN.md with every
# command and all arguments (base preparation, each fit; rerun_commands.py), the AMPS commit and the commits of the
# code, and code/ (the code that ran).
# Login node, no GPU:   figures_12m_positive_scan.sh
set -u
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; J=$COPS_OUTPUT_ROOT/jobs
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; PD=$R/positive; N=C30_n4000_exp_cliff_w10x5_grain_cross10_lam
AN=$COPS_OUTPUT_ROOT/analysis/positive_lambda_scan_12m
O=$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_positive_lambda_scan; TB=$O/tables; C=$O/code
mkdir -p "$TB" "$C"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
# lambda, job, line style, label
CASES=("0 122307 dotted" "1e-5 122311 dashed" "1e-4 122320 dashdot" "1e-3 122310 solid" "1e-2 122319 solid" "1e-1 122312 dashed")
FITS=()
for c in "${CASES[@]}"; do
    set -- $c; L=$1; JOB=$2; ST=$3
    .venv/bin/python $S/positive_table.py "$PD/$N$L/seg01" "$TB/lam$L.dat" || exit 1
    LAB="lambda $L"; [ "$L" = 0 ] && LAB="lambda 0 (max_nfev 4000, not converged)"
    FITS+=(--fit "$LAB=$TB/lam$L.dat:$ST")
done
T='SHAP8 12M, positive activity exp(Y a + c T$_{\rm cliff}$), night 0.02, floor 0.02, original T held, W10x5 from 10:33:10, gauge +1.5 m Z'
.venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R/bundle/restart_templates.npz --overwrite \
    --out "$O/12M_positive_lambda_scan" --ref "$TB/lam1e-2.dat" --title "$T" "${FITS[@]}" \
    && echo "wrote $O/12M_positive_lambda_scan_*" || exit 3
cp -p "$AN/scan_table.txt" "$AN/lcurve.eps" "$AN/lcurve.png" "$O/"
# code copy: the job scripts and modules as they ran (provenance copies), the figure and analysis scripts
for c in "${CASES[@]}"; do
    set -- $c; d=$PD/$N$1/provenance-$2; k=lam$1-$2; mkdir -p "$C/$k"
    cp -p "$d"/*.sbatch "$d"/*.py "$d"/model_ext.sha256 "$d"/git.txt "$d"/commands.sh "$C/$k/"; cp -rp "$d/kernels" "$C/$k/"
done
B=$R/provenance/prepare-122291; mkdir -p "$C/base-prepare-122291"
cp -p "$B"/*.sbatch "$B"/*.py "$B"/git.txt "$B"/commands.sh "$B"/executables.sha256 "$C/base-prepare-122291/"
mkdir -p "$C/figures"
cp -p $S/figures_12m_positive_scan.sh $S/positive_table.py $S/model_obs_figures.py $S/rerun_commands.py "$AN/scan_summary.py" "$C/figures/"
(cd "$C" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
RC=".venv/bin/python $S/rerun_commands.py"
{
cat <<MD
# How to rerun: 12M positive-activity lambda scan

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_12m_positive_scan.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never
delete earlier attempts (move them aside as \`*.incomplete_<UTC>\`).

What was fitted (user decisions 2026-11-07): H2O activity A = exp(sum_j a_j Y_j + c T_cliff) per facet (25 real
spherical harmonics of degree <= 4, template cliff_ge_60deg, c >= 0; no fine-tuning templates), night/shadow ratio
0.02 for H2O and CO2, H2O flux floor F = max(q A, 0.02 q_max A_orig) (CO2 none); start and prior = the
emission-weighted (w = area A_orig) LS projection of ln A_orig; Tikhonov D^2 = sum area (s - s0)^2 / sum area; original
temperature law held (182.09922 / 136.11149 / 133 K); kappa free; grain A start 6.96024299190473 uPa, B = 0;
crossing weight 10; observations from 10:33:10 UTC compared with W10x5 (mean of five 10-s running averages ending at
T, T-2, ..., T-8 s, T = L2 STOP); gauge offset +1.5 m along +Z of the spacecraft (base
\`$R/base\`, 1551 points = 1150 observations + 401 model-only 1-s grid points, removed from the figures by
\`positive_table.py\`); max_nfev 4000. Termination: lambda 0 reached max_nfev (not converged); the others \`ftol\`.
Table and L-curve: \`scan_table.txt\`, \`lcurve.png\` (copies of \`$AN\`, see its README.md).

The commands below are complete (every argument as run), so they do not depend on the current job scripts. Code:
the "commit" column names the git commit (in \`srcCG/LiouvilleCOPS_cuda/model_ext\`) whose version-controlled copy
is identical to the code that ran ("not in git": only in \`code/\` here); \`code/\` has the files as they ran, with
\`SHA256SUMS\`. Job script for the same fits: \`sbatch --export=ALL,LAMBDA=<lambda> $J/positive_fit.sbatch\`, not needed
for the commands below.

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_12m_positive_scan.sh
\`\`\`

The scan table and L-curve (login node): see \`$AN/README.md\`.

## A1. The fitting base (job 122291, one GPU; inside a Slurm job)

AMPS commit \`$(head -1 "$B/git.txt")\`. Commands as run (\`$B/commands.sh\`; every argument resolved):

\`\`\`bash
$(cat "$B/commands.sh")
\`\`\`

Consistency of this base with the unchanged GPU forward model (job 122304, \`$J/base_consistency_check.sbatch\`,
\`$R/checks/consistency-122304.provenance/commands.sh\`): NG 2.02e-9, RG 3.70e-9 (tolerance 1e-6).

## A2. The six fits (each on one GPU node; inside a Slurm job)

MD
for c in "${CASES[@]}"; do
    set -- $c
    $RC --fit "$PD/$N$1/seg01" --provenance "$PD/$N$1/provenance-$2" --label "lambda $1 (job $2)"
done
} > "$O/HOW-TO-RERUN.md"
echo "wrote $O/HOW-TO-RERUN.md"
