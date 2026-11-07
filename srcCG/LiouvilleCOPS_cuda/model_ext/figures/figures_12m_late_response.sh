#!/bin/bash
# ==== version-controlled copy: header added 2026-11-07, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_12m_late_response.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_12m_late_response.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: c1d7e25c951c058901e72ddd3eb2ba6ad54ce4ac0368272357d897658fa98894
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure script. Figures of a late_response job + HOW-TO-RERUN.md (the job's commands.sh, AMPS commit, git
#   commit of each code file via rerun_commands.commit_of) + code/ into
#   figures_cluster/12M/model_vs_observations_late_response/. User 2026-11-07.
# RUN: figures_12m_late_response.sh <late_response job id>   (login node, no GPU)
# ==== end of header ====
# Figures of the late-window forward check (job late_response.sbatch; user 2026-11-07: "1: both" time-response windows,
# "3: make the sensitivity check" of the gauge offset, "for the forward run use two 12M best fits with and without
# extra weight of the last point observations"), against the COPS observations, with HOW-TO-RERUN.md (every command
# with all arguments, the AMPS commit, the commits of the code) and code/ (the code that ran).
# Login node, no GPU:   figures_12m_late_response.sh <job id>
set -u
JOB=${1:?usage: figures_12m_late_response.sh <late_response job id>}
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; R=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier
D=$R/late_response/$JOB; P=$D/provenance; E=$D/eval
O=$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_late_response; C=$O/code
[ -e "$E/summary.txt" ] || { echo "no evaluation in $E" >&2; exit 2; }
mkdir -p "$O" "$C/figures" "$C/job-$JOB"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
.venv/bin/python -W ignore $S/late_response_figures.py --eval "$E" --out "$O/12M_late_response" --overwrite \
    --pair c3,c3_last5_w10_constant \
    --label c3="cycle 3 (no extra weight)" --label c3_last5_w10_constant="cycle 3 + last 5 x10" \
    --label c3_last6_w1to20_linear="cycle 3 + last 6 ramp 1-20" --label torig_c6="original T cycle 6" \
    --label torig_c6_last6_w1to20_linear="original T c6 + last 6 ramp 1-20" || exit 3
cp -p "$E/summary.txt" "$E/validation.txt" "$O/"
cp -p "$P"/*.sbatch "$P"/*.py "$P/model_ext.sha256" "$P/git.txt" "$P/commands.sh" "$C/job-$JOB/"
cp -p "$S/figures_12m_late_response.sh" "$S/late_response_figures.py" "$C/figures/"
(cd "$C" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
# git commit holding an identical version-controlled copy (header block removed), else "not in git"
CM=$(.venv/bin/python - "$C" <<'PY'
import hashlib, sys
from pathlib import Path
sys.path.insert(0, "/data001/highresolu_vt/cops_output/figures_cluster/scripts")
from rerun_commands import commit_of
for p in sorted(Path(sys.argv[1]).rglob("*")):
    if p.is_file() and p.suffix in (".py", ".sbatch", ".sh"):
        c = commit_of(p.name, hashlib.sha256(p.read_bytes()).hexdigest())
        print("| `%s` | `%s` | %s |" % (p.relative_to(sys.argv[1]), hashlib.sha256(p.read_bytes()).hexdigest()[:16], "`%s`" % c if c else "**not in git**"))
PY
)
{
cat <<MD
# How to rerun: late-window forward check (time response, gauge offset)

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_12m_late_response.sh $JOB\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never delete
earlier attempts (move them aside as \`*.incomplete_<UTC>\`).

What it is: no fit. The finished 12M fits keep their per-facet activity, temperature law, CO2 scale and grain; the
geometry is recomputed on a 1-s grid of the late window (10:32:46.282-10:39:26.282 UTC, plus the 95 original
observation times) for gauge offsets along +Z of the spacecraft (sensitivity check: the COPS position relative to the
trajectory point is not known), and the model is averaged as the gauges measure (ROSINA COPS instrument modes
RO-ROS-MAN-1019 sect. 4.2): W10 = 10-s running average ending at the L2 STOP time, W10x5 = mean of five W10 ending
2 s apart (DPU 5-readout average); inst = the model at the time stamp (what the fits used); W4 (L2 interval) in the
tables only. Products: \`$D\` (bases \`base_dz*\`, evaluation \`eval/\`, provenance).
Validation (\`validation.txt\`): the dz0 base at the original times reproduces each fit's and each fine tuning's own
model values.

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_12m_late_response.sh $JOB
\`\`\`

## A1. The job (one GPU node, 3 GPUs; inside a Slurm job). AMPS commit: \`$(head -1 "$P/git.txt")\`; tracked changes: $(tail -n +2 "$P/git.txt" | grep -v '^??' | tr -s ' ' | tr '\n' ',' | sed 's/,$//')

Commands as run (\`$P/commands.sh\`; every argument resolved):

\`\`\`bash
$(cat "$P/commands.sh")
\`\`\`

The submission: \`sbatch $COPS_OUTPUT_ROOT/jobs/late_response.sbatch\` (DZS="0 1.5 3.0", 3 GPUs, 48 CPUs, 150 GB).

## Code that ran (\`code/\`, with \`SHA256SUMS\`) and the git commit with an identical copy in \`srcCG/LiouvilleCOPS_cuda/model_ext\`

| file | SHA256 (first 16) | commit |
|---|---|---|
$CM

The GPU forward model \`native/liouville_cops_cuda\` is built from the tracked sources of the AMPS commit above
(its SHA256 is in \`code/job-$JOB/model_ext.sha256\`).
MD
} > "$O/HOW-TO-RERUN.md"
echo "wrote $O/HOW-TO-RERUN.md"
