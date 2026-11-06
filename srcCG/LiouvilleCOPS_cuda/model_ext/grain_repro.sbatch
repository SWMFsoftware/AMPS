#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/grain_repro.sbatch (on Voyager: /data001/highresolu_vt/cops_output/model_ext/grain_repro.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: f4938cd08c8eb3bd8fabe526d0808e281de0af3ed2cc3004eb1668a2aa5f000e
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: reruns a finished grain fit with exactly its recorded command (fit_manifest model_extension.command, only
#   --out changed), then fine-tunes it. Jobs 121862, 121863 (bitwise identical to the originals).
# RUN: sbatch --export=ALL,COPS_RUN=<run>,SRC_FIT=<fit dir name> grain_repro.sbatch
# ==== end of header ====
# Reproduce one finished grain fit with exactly its recorded command (user request 2026-11-03: "start a new fit that
# reproduces that fit for further minor tweaks"), then fine-tune it as grain_finetune.sbatch does.
#   sbatch --export=ALL,COPS_RUN=<run>,SRC_FIT=<fit dir name> grain_repro.sbatch
# The command is read from <run>/grain/fits/<SRC_FIT>/fit_manifest.json["model_extension"]["command"]; only --out
# changes: <run>/grain/repro/<SRC_FIT>_repro_<jobid>/ (fine tuning: .../finetune/<that name>_m5_i5_cpu/).
# Same node type, GPU count, precision, code (model_ext frozen, its SHA256 is checked against the source fit) and
# start; expect agreement at rounding level, not bitwise (say "bitwise identical" only after checking).
#SBATCH --job-name=cops_grain_repro
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=03:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/model_ext/logs/cops_grain_repro-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=${COPS_RUN:?set COPS_RUN}; SRC=${SRC_FIT:?set SRC_FIT}
X=$COPS_OUTPUT_ROOT/model_ext; PY="$APP/.venv/bin/python -W ignore"
S=$R/grain/fits/$SRC; N=${SRC}_repro_$SLURM_JOB_ID
O=$R/grain/repro/$N; F=$R/grain/repro/finetune/${N}_m5_i5_cpu; P=$R/grain/repro/provenance-$SLURM_JOB_ID
[ -f "$S/fit_manifest.json" ] || { echo "grain_repro: $S/fit_manifest.json missing" >&2; exit 2; }
[ -e "$O" ] && { echo "grain_repro: $O exists" >&2; exit 2; }
mkdir -p "$P" "$(dirname "$F")"
cd "$APP" || exit 1
cp -p "$0" "$P/grain_repro.sbatch"; cp -p "$X/grain_crossing.py" "$X/grain_fit.py" "$X/grain_finetune.py" "$P/"
sha256sum "$X"/*.py > "$P/model_ext.sha256"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_|CUPY_|COPS_|OMP_|OPENBLAS_|MKL_|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH)=' | sort > "$P/environment.txt"
# the source fit's command, with only --out replaced; the grain code must be the one that made the source fit
CMD=$($PY - "$S/fit_manifest.json" "$O" <<'EOF'
import json, shlex, sys
m = json.load(open(sys.argv[1]))["model_extension"]
argv = shlex.split(m["command"]) if isinstance(m["command"], str) else list(m["command"])
argv[argv.index("--out") + 1] = sys.argv[2]
print(shlex.join(argv))
EOF
) || exit 2
WANT=$($PY -c "import json,sys; print(json.load(open(sys.argv[1]))['model_extension']['module_sha256'])" "$S/fit_manifest.json")
HAVE=$(sha256sum "$X/grain_crossing.py" | cut -d' ' -f1)
[ "$WANT" = "$HAVE" ] || { echo "grain_repro: grain_crossing.py differs from the one of $S ($HAVE vs $WANT)" >&2; exit 2; }
echo "$CMD" > "$P/command.txt"
echo "$(date -u +%FT%TZ) fit $N: start"; echo "  $CMD"
$PY $CMD > "$P/fit.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fit $N: exit $rc"
[ $rc -eq 0 ] || exit $rc
echo "$(date -u +%FT%TZ) fine tuning: start"
export COPS_CUDA_RECORD_STORE=$R/store
DEV=$(nvidia-smi --query-gpu=index --format=csv,noheader | paste -sd, -)
$PY "$X/grain_finetune.py" "$R/base" --restart-from "$O/fit_restart.npz" --minutes 5 --intervals 5 \
    --method cpu --max-nfev 2000 --devices "$DEV" --precision mixed --out "$F" > "$P/finetune.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) fine tuning: exit $rc"
exit $rc
