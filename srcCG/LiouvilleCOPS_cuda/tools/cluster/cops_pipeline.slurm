#!/bin/bash
# Group "@GROUP@" of a tools/cluster/cops_pipeline.py campaign. Template: `cops_pipeline.py submit` fills the
# @...@ fields (nodes and GPUs from --nodes/--gpus-per-node or the partition's GRES; nothing is hardwired)
# and writes the result to <restart dir>/slurm/@GROUP@.sbatch. Lines marked SITE are for the partition,
# account and modules of the cluster (submit --slurm-partition/--slurm-account/--slurm-module fill them).
#SBATCH --job-name=@JOB_NAME@
#SBATCH --nodes=@NODES@
@TASK_LINES@
#SBATCH --mem=0                        # all node memory: pinned host buffers for streamed records
#SBATCH --time=@TIME@
#SBATCH --output=@OUTPUT@
@SITE_LINES@

set -uo pipefail
@MODULE_LINES@
export PYTHONUNBUFFERED=1
cd @APP@
nvidia-smi topo -m 2>/dev/null | head -20 || true

# The group's stages; the runner resumes at the first incomplete stage and exits 0 (done),
# 75 (stopped before the wall-time deadline, resumable) or 1/2 (failed, refused).
@RUN@
rc=$?

@CHAIN@
exit $rc
