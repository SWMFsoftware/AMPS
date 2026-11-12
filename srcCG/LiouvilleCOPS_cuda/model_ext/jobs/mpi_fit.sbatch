#!/bin/bash
# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/mpi_fit.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/mpi_fit.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: d07cc807e52d083d5198541fb6d61f774703044fcd29926c5a9f61c19b8c54d5
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. ONE fit split across one or two GPU nodes (user decisions 2026-11-10 "one fit split across both nodes",
#   "make it also runnable on one node if only one node is available"): generic launcher of a model_ext driver over
#   tools/gpu_fit.py --transport mpi (srun --mpi=pmix, one rank per GPU, every rank through model_ext/mpi_support.py);
#   SPEC file, NODES=1|2|auto, ACTIVITY device|host, resume across node counts, per-rank provenance.
# RUN: jobs/mpi_fit_submit.sh 1|2|auto <spec> MODEL_EXT=<snapshot> ...   This version (commands.sh exports, 2026-11-10
#   10:27 UTC) has not been submitted; the MPI tests ran the earlier versions (.before_nfevenv_*: 123629-123632, 123638;
#   .before_commandsexport_*: 123645-123649).
# ==== end of header ====
# ONE fit split across one or two GPU nodes (user decisions 2026-11-10: "one fit split across both nodes" and "make it
# also runnable on one node if only one node is available"). Generic launcher for a model_ext driver command (a
# driver that wraps the unchanged tools/gpu_fit.py in-process: clip_refit.py, grain_fit*.py, ...). No physics here.
#
# How the split works (tracked code, unchanged): tools/gpu_fit.py --transport mpi (quick_guide.md 4e2,
# cops_cuda/transport.py): one MPI rank per GPU; rank 0 runs the fitter and the first block of rows on its own GPU,
# ranks >= 1 serve their blocks; the rows are split so that the kernel records are balanced. EVERY rank runs the
# same driver command (through model_ext/mpi_support.py: per-rank provenance, node memory, optional instrumentation),
# so whatever the driver installs at start-up exists in every process; gpu_fit.main sends ranks >= 1 into
# transport.mpi_worker_main. Launch: srun --mpi=pmix (Open MPI 4.1.7rc1 of /usr/mpi/gcc = mpi4py's MPI; UCX over
# the HDR InfiniBand, selected by default; checked by jobs/mpi_probe.sbatch, job 123623).
#
# Submission (jobs/mpi_fit_submit.sh writes the sbatch line and logs it in model_ext/logs/submissions.txt):
#   NODES=2     sbatch --nodes=2   --export=ALL,SPEC=<spec>,MODEL_EXT=<snapshot> mpi_fit.sbatch   (16 ranks)
#   NODES=1     sbatch --nodes=1   --export=ALL,SPEC=<spec>,MODEL_EXT=<snapshot> mpi_fit.sbatch   (8 ranks)
#   NODES=auto  sbatch --nodes=1-2 --export=ALL,SPEC=<spec>,MODEL_EXT=<snapshot> mpi_fit.sbatch   (Slurm gives 2 nodes
#               if both are free at the start, else 1; ranks = 8 x SLURM_JOB_NUM_NODES)
#
# Settings (--export=ALL,VAR=value,...; values must not contain commas):
#   SPEC        (required) bash file sourced by this script; defines NAME (fit name), D (fit directory: seg<NN>,
#               provenance-<job>, segments.tsv, final), NFEV (evaluations over all segments), DRIVER (file name in
#               MODEL_EXT), optional REQUIRED=(input files that must exist), optional ACTIVITY_DEFAULT, and a function
#               fit_args that sets the array ARGS = the driver's arguments (driver options, gpu_fit options, '--',
#               fit_cops arguments with --out "$OUT") from SEG, OUT, PREV_RESTART (empty in segment 1) and DEVICES;
#               optional SPEC_EXPORTS=(names of environment variables the spec exports; written into commands.sh).
#               This launcher inserts before the first '--': --transport, --activity, --max-nfev LEFT, --segment SEG,
#               --wall-time (and --restart-from PREV_RESTART from segment 2 on); ARGS must not contain them.
#   MODEL_EXT   model_ext directory (default $COPS_OUTPUT_ROOT/model_ext); production: a frozen snapshot. Its *.py are
#               copied into provenance-<job>/model_ext/ and the fit runs from that copy (SHA256SUMS checked if present)
#   NODES       1 | 2 | auto (default auto): checked against SLURM_JOB_NUM_NODES
#   TRANSPORT   mpi (default) | local (one node only: one process with spawned workers, the one-node path of all
#               earlier fits; for comparisons)
#   ACTIVITY    device | host (default: the spec's ACTIVITY_DEFAULT, else gpu_fit's default for the transport: mpi ->
#               device (user decision Q3), local -> host). host: exact host field broadcast every evaluation (96 MB at
#               12M), bitwise the local fit; device: every GPU forms the activity of its own facets, the same bytes
#               for any rank count. A resume must keep the activity of the earlier segments (ALLOW_ACTIVITY_CHANGE=1
#               overrides; recorded).
#   INSTRUMENT  1 (default) | 0: model_ext/mpi_support.py times every evaluation, counts the bytes per message and
#               records the optimizer's cost history with SHA-256 of x, residuals and Jacobian (observation only)
#   WALL        gpu_fit --wall-time (default slurm = the job's end); WALL_MARGIN, SAVE_SECONDS: gpu_fit --wall-margin,
#               --save-seconds (default: gpu_fit's 300 s / 600 s)
#   CHAIN       1: a segment directory without restart file stops the job (exit 3) instead of being moved aside
#   DRYRUN      1: print the segment and the driver command, run nothing
#   MPI_SUPPORT the per-rank wrapper (default: MODEL_EXT/mpi_support.py if the snapshot has it, else
#               $COPS_OUTPUT_ROOT/model_ext/mpi_support.py); copied into provenance and run from there
#
# Segments / resume across node counts: the fit stops at its best evaluated state before the job's end (gpu_fit wall
# clock), writes an interrupted fit_restart.npz and exits 75. Resubmitting the same SPEC (with ANY node count: NODES=1,
# 2 or auto) continues with segment k+1 (gpu_fit --restart-from seg<k>/fit_restart.npz), --max-nfev = NFEV - the
# evaluations of the earlier segments. The evaluations of a state are identical for 8 and 16 ranks with the same
# ACTIVITY (row blocks are independent; device activity is partition-invariant), so the node count may change between
# segments; a segmented fit is not identical to one continuous run (new least_squares start). $D/segments.tsv records
# every segment (job, nodes, node list, ranks, transport, activity, exit). A segment directory without restart file is
# moved aside as seg<k>.incomplete_<UTC> (never deleted). A completed fit is not rerun (exit 0).
#
# Provenance: $D/provenance-<job>/ : this script, the spec, model_ext copy (+ SHA-256), commands.sh (every command
# with all arguments resolved, runnable as is), git.txt, module_list.txt, toolchain.txt, pip_list.txt, environment.txt,
# segment.txt, per node nodes/<host>.txt (nvidia-smi -L, IB, free) and nvidia_smi_q_<host>.txt, per rank
# ranks/rank<NNN>.json (host, rank, local rank, GPU index + UUID + PCI bus, CUDA_VISIBLE_DEVICES, pid, oom_score_adj,
# MPI library) + rank<NNN>.end.json (peak RSS), ranks/mem_<host>.tsv (node memory every 10 s), ranks/task-<NN>.log
# (stdout+stderr of every rank; fit.log -> ranks/task-00.log), comm_rank000.json / optimizer_history.json
# (INSTRUMENT=1), free/df at start and end. Every process runs with oom_score_adj 1000.
#SBATCH --job-name=cops_mpi_fit
#SBATCH --nodes=1-2
#SBATCH --ntasks-per-node=8
#SBATCH --cpus-per-task=8
#SBATCH --gpus-per-node=8
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=08:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --exclude=gpu003
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_mpi_fit-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
PY="$APP/.venv/bin/python -W ignore"
SELF=$(readlink -f "$0")
die() { echo "mpi_fit: $*" >&2; exit 2; }
# ---- settings
SPEC=${SPEC:-}; [ -n "$SPEC" ] && [ -f "$SPEC" ] || die "set SPEC=<spec file> (got '${SPEC}')"
SPEC=$(readlink -f "$SPEC")
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
NODES=${NODES:-auto}; TRANSPORT=${TRANSPORT:-mpi}; INSTRUMENT=${INSTRUMENT:-1}; CHAIN=${CHAIN:-0}
WALL=${WALL:-slurm}; WALL_MARGIN=${WALL_MARGIN:-}; SAVE_SECONDS=${SAVE_SECONDS:-}; ALLOW_ACTIVITY_CHANGE=${ALLOW_ACTIVITY_CHANGE:-0}
NN=${SLURM_JOB_NUM_NODES:?not in a Slurm job}
case $NODES in
    1|2) [ "$NN" = "$NODES" ] || die "NODES=$NODES but the job has $NN node(s): submit with --nodes=$NODES";;
    auto) [ "$NN" = 1 ] || [ "$NN" = 2 ] || die "NODES=auto expects 1 or 2 nodes, the job has $NN";;
    *) die "NODES must be 1, 2 or auto";;
esac
case $TRANSPORT in mpi) ;; local) [ "$NN" = 1 ] || die "TRANSPORT=local runs on one node (the job has $NN)";; *) die "TRANSPORT=mpi|local";; esac
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then G=$(echo "$CUDA_VISIBLE_DEVICES" | tr , "\n" | wc -l); else G=$(nvidia-smi -L | wc -l); fi
RPN=${SLURM_NTASKS_PER_NODE:-8}; RPN=${RPN%%(*}
[ "$TRANSPORT" = local ] || [ "$RPN" = "$G" ] || die "ranks per node ($RPN) must equal the GPUs per node ($G)"
NR=$((NN * G)); [ "$TRANSPORT" = local ] && NR=1
DEVICES=$(seq -s, 0 $((G - 1)))
# ---- the fit spec
NAME=; D=; DRIVER=; REQUIRED=(); ACTIVITY_DEFAULT=; SPEC_EXPORTS=()    # NFEV (and the spec's own variables) come from --export
# shellcheck disable=SC1090
source "$SPEC" || die "cannot source $SPEC"
[ -n "$NAME" ] && [ -n "$D" ] && [ -n "$NFEV" ] && [ -n "$DRIVER" ] || die "$SPEC must set NAME, D, NFEV, DRIVER"
declare -F fit_args > /dev/null || die "$SPEC must define the function fit_args"
if [ -z "${ACTIVITY:-}" ]; then
    if [ -n "$ACTIVITY_DEFAULT" ]; then ACTIVITY=$ACTIVITY_DEFAULT; elif [ "$TRANSPORT" = mpi ]; then ACTIVITY=device; else ACTIVITY=host; fi
fi
case $ACTIVITY in host|device) ;; *) die "ACTIVITY=host|device";; esac
[ -f "$MODEL_EXT/$DRIVER" ] || die "driver not found: $MODEL_EXT/$DRIVER"
if [ -z "${MPI_SUPPORT:-}" ]; then
    if [ -f "$MODEL_EXT/mpi_support.py" ]; then MPI_SUPPORT=$MODEL_EXT/mpi_support.py; else MPI_SUPPORT=$COPS_OUTPUT_ROOT/model_ext/mpi_support.py; fi
fi
MPI_SUPPORT=$(readlink -f "$MPI_SUPPORT"); [ -f "$MPI_SUPPORT" ] || die "missing $MPI_SUPPORT"
for f in "${REQUIRED[@]}"; do [ -e "$f" ] || die "missing: $f"; done
mkdir -p "$D"
exec 9> "$D/.mpi_fit.lock"; flock -n 9 || die "another job holds $D/.mpi_fit.lock (the same fit is running)"
# ---- segment state: the last seg<NN>, its restart's status, the evaluations used so far
state() {  # $1 = directory, $2 = segment (optional: the last one)
$PY - "$@" <<'PY'
import glob, json, os, sys
import numpy as np
d = sys.argv[1]
segs = sorted(int(os.path.basename(p)[3:]) for p in glob.glob(os.path.join(d, "seg[0-9][0-9]")))
k = int(sys.argv[2]) if len(sys.argv) > 2 else (segs[-1] if segs else 0)
if k == 0:
    print(0, "none", 0); sys.exit()
rs = os.path.join(d, "seg%02d" % k, "fit_restart.npz")
if not os.path.isfile(rs):
    print(k, "norestart", 0); sys.exit()
with np.load(rs, allow_pickle=False) as z:
    meta = json.loads(str(z["meta_json"]))
used = sum(int(s.get("nfev", 0)) for s in meta.get("segments", []))
print(k, meta.get("status", "unknown"), used)
PY
}
read -r LAST LSTATE LUSED < <(state "$D")
[ -n "${LAST:-}" ] || die "cannot read the segment state of $D"
MOVED=""
if [ "$LSTATE" = norestart ]; then                     # a killed or failed attempt
    [ "$CHAIN" = 1 ] && { echo "$D/seg$(printf %02d "$LAST") has no restart file (a failed attempt): chain stopped (CHAIN=1)" >&2; exit 3; }
    MOVED="$D/seg$(printf %02d "$LAST").incomplete_$(date -u +%Y%m%dT%H%M%SZ)"   # move aside, run that segment again
    mv "$D/seg$(printf %02d "$LAST")" "$MOVED"
    LAST=$((LAST - 1))
    if [ "$LAST" -ge 1 ]; then read -r LAST LSTATE LUSED < <(state "$D" "$LAST"); else LSTATE=none; LUSED=0; fi
fi
final_link() {
    local seg=seg$(printf %02d "$1")
    if [ -L "$D/final" ]; then
        [ "$(readlink "$D/final")" = "$seg" ] || { echo "$D/final -> $(readlink "$D/final"), expected $seg" >&2; return 3; }
    else ln -s "$seg" "$D/final"; fi
    echo "$(date -u +%FT%TZ) $NAME complete: $D/final -> $seg"
}
case $LSTATE in
    none) SEG=1; LEFT=$NFEV;;
    interrupted) SEG=$((LAST + 1)); LEFT=$((NFEV - LUSED));;
    complete) echo "$(date -u +%FT%TZ) $NAME: segment $LAST is complete ($LUSED evaluations); nothing to do"; final_link "$LAST"; exit $?;;
    *) echo "$D/seg$(printf %02d "$LAST")/fit_restart.npz has status '$LSTATE'; not continued" >&2; exit 3;;
esac
if [ "$LEFT" -le 0 ]; then echo "$(date -u +%FT%TZ) $NAME: $LUSED of $NFEV evaluations used; nothing to do"; final_link "$LAST"; exit $?; fi
PREV_RESTART=""; [ "$SEG" -gt 1 ] && PREV_RESTART=$D/seg$(printf %02d "$LAST")/fit_restart.npz
# the activity of the earlier segments (segments.tsv column 9)
if [ "$SEG" -gt 1 ] && [ -f "$D/segments.tsv" ]; then
    PREV_ACT=$(awk -F'\t' -v s="$LAST" '$1 == s && $13 != "" {a = $9} END {print a}' "$D/segments.tsv")
    if [ -n "$PREV_ACT" ] && [ "$PREV_ACT" != "$ACTIVITY" ] && [ "$ALLOW_ACTIVITY_CHANGE" != 1 ]; then
        die "segment $LAST ran with ACTIVITY=$PREV_ACT, this job has ACTIVITY=$ACTIVITY (ALLOW_ACTIVITY_CHANGE=1 to override)"
    fi
fi
OUT=$D/seg$(printf %02d "$SEG")
[ -e "$OUT" ] && die "exists: $OUT"
# ---- the command
ARGS=()
fit_args || die "fit_args failed"
split=-1; for i in "${!ARGS[@]}"; do [ "${ARGS[$i]}" = "--" ] && { split=$i; break; }; done
[ "$split" -ge 0 ] || die "ARGS has no '--'"
for t in "${ARGS[@]:0:$split}"; do
    case $t in --transport|--transport=*|--activity|--activity=*|--max-nfev|--max-nfev=*|--segment|--segment=*|\
               --wall-time|--wall-time=*|--deadline|--deadline=*|--wall-margin|--save-seconds|--restart-from|--restart-from=*)
        die "the spec's ARGS must not contain $t (the launcher sets it)";; esac
done
LAUNCH=(--transport "$TRANSPORT" --activity "$ACTIVITY" --max-nfev "$LEFT" --segment "$SEG" --wall-time "$WALL")
[ -n "$WALL_MARGIN" ] && LAUNCH+=(--wall-margin "$WALL_MARGIN")
[ -n "$SAVE_SECONDS" ] && LAUNCH+=(--save-seconds "$SAVE_SECONDS")
[ "$SEG" -gt 1 ] && LAUNCH+=(--restart-from "$PREV_RESTART")
ARGS=("${ARGS[@]:0:$split}" "${LAUNCH[@]}" "${ARGS[@]:$split}")
if [ "${DRYRUN:-0}" = 1 ]; then      # print what would run, change nothing (the lock file and $D exist)
    echo "DRYRUN: $NAME segment $SEG ($LEFT of $NFEV evaluations; last $LAST $LSTATE $LUSED), $NN node(s), $NR rank(s), transport $TRANSPORT, activity $ACTIVITY"
    echo -n "DRYRUN: "; printf '%q ' "$MODEL_EXT/$DRIVER" "${ARGS[@]}"; echo; exit 0
fi
# ---- provenance; the fit runs from the provenance copy of model_ext
P=$D/provenance-$SLURM_JOB_ID; X=$P/model_ext
mkdir -p "$X" "$P/ranks" "$P/nodes"; cd "$APP" || exit 1
cp -p "$SELF" "$P/mpi_fit.sbatch"; cp -p "$SPEC" "$P/spec.sh"
cp -p "$MODEL_EXT"/*.py "$X/"
[ -f "$MODEL_EXT/SHA256SUMS" ] && { cp -p "$MODEL_EXT/SHA256SUMS" "$X/"; (cd "$X" && sha256sum -c --quiet SHA256SUMS) > "$P/model_ext_check.txt" 2>&1 \
    || die "the copy of $MODEL_EXT does not match its SHA256SUMS (see $P/model_ext_check.txt)"; }
[ "$MPI_SUPPORT" = "$MODEL_EXT/mpi_support.py" ] || cp -p "$MPI_SUPPORT" "$X/mpi_support.py"
(cd "$X" && sha256sum *.py) > "$P/run_copy.sha256"
sha256sum "$MODEL_EXT"/*.py "$MPI_SUPPORT" "$SPEC" > "$P/model_ext.sha256" 2>&1
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; mpirun --version | head -1
  ompi_info | grep -E 'MCA pml: ucx|MCA pmix'; ucx_info -v | head -1; srun --version; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
free -g > "$P/free_at_start.txt"; df -h "$D" > "$P/df_at_start.txt" 2>&1
NODE_LIST=$(scontrol show hostnames "$SLURM_JOB_NODELIST" | paste -sd, -)
# per node (one task per node, sharing the job's resources with the later steps)
srun --overlap --nodes="$NN" --ntasks="$NN" --ntasks-per-node=1 --cpus-per-task=1 bash -c '
h=$(hostname -s); { echo "== $h $(date -u +%FT%TZ)"; nvidia-smi -L; echo "== ibstat"; ibstat 2>&1 | grep -E "CA |State|Rate|Link layer"
  echo "== free -g"; free -g; echo "== CUDA_VISIBLE_DEVICES ${CUDA_VISIBLE_DEVICES:-unset}"; uptime; } > "$1/nodes/$h.txt" 2>&1
nvidia-smi -q > "$1/nvidia_smi_q_$h.txt" 2>&1' _ "$P"
export COPS_MPI_RANK_DIR=$P/ranks COPS_MPI_INSTRUMENT=$INSTRUMENT
env | grep -E '^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z_]*|SPEC|NODES|TRANSPORT|ACTIVITY|INSTRUMENT|WALL[A-Z_]*|SAVE_SECONDS|CHAIN|MODEL_EXT|MPI_SUPPORT|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH|SLURM_JOB_[A-Z_]*|SLURM_NTASKS[A-Z_]*|OMPI_[A-Z_]*|UCX_[A-Z_]*)=' | sort > "$P/environment.txt"
OOM=(bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _)
RUN=($PY "$X/mpi_support.py" "$X/$DRIVER" "${ARGS[@]}")
if [ "$TRANSPORT" = mpi ]; then
    LAUNCHER=(srun --mpi=pmix --kill-on-bad-exit=1 --nodes="$NN" --ntasks="$NR" --ntasks-per-node="$G" --cpus-per-task="${SLURM_CPUS_PER_TASK:-8}"
              --distribution=block --output="$P/ranks/task-%2t.log")
    ln -s ranks/task-00.log "$P/fit.log"
else
    LAUNCHER=()
fi
# commands.sh: every command of this job with all arguments resolved, runnable as is (inside an allocation of the same
# shape for srun)
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $NODE_LIST ($NN node(s), $NR rank(s), transport $TRANSPORT, activity $ACTIVITY), $(date -u +%FT%TZ); AMPS $(git -C "$AMPS" rev-parse HEAD)"
  grep '^export ' "$SELF" | grep -v '^export COPS_MPI_RANK_DIR='; echo "export COPS_MPI_RANK_DIR=$P/ranks COPS_MPI_INSTRUMENT=$INSTRUMENT"
  for v in "${SPEC_EXPORTS[@]}"; do printf 'export %s=%q\n' "$v" "${!v-}"; done; echo "cd $APP"
  printf '%q ' "${LAUNCHER[@]}" "${OOM[@]}" "${RUN[@]}"; [ "$TRANSPORT" = local ] && printf '> %q 2>&1' "$P/fit.log"; echo; } > "$P/commands.sh"
{ echo "fit $NAME, directory $D, spec $SPEC"
  echo "model_ext copied and run: $MODEL_EXT (MODEL_EXT) -> $X; wrapper $MPI_SUPPORT"
  echo "job $SLURM_JOB_ID: $NN node(s) $NODE_LIST, $NR rank(s), $G GPUs per node, transport $TRANSPORT, activity $ACTIVITY, NODES=$NODES"
  echo "segment $SEG; previous segment $LAST ($LSTATE, $LUSED evaluations used); evaluations for this segment $LEFT of $NFEV"
  [ -n "$PREV_RESTART" ] && echo "resume from $PREV_RESTART"
  [ -n "$MOVED" ] && echo "moved aside: $MOVED"; } | tee "$P/segment.txt"
[ -f "$D/segments.tsv" ] || printf 'seg\tjob\tstart_utc\tend_utc\tnodes\tnode_list\tranks\ttransport\tactivity\tmodel_ext\tevaluations_left\tstatus\texit\n' > "$D/segments.tsv"
T0=$(date -u +%FT%TZ)
echo "$T0 job $SLURM_JOB_ID: fit $NAME segment $SEG ($LEFT evaluations) on $NODE_LIST, $NR rank(s)"
if [ "$TRANSPORT" = mpi ]; then
    "${LAUNCHER[@]}" "${OOM[@]}" "${RUN[@]}"; rc=$?
else
    "${OOM[@]}" "${RUN[@]}" > "$P/fit.log" 2>&1; rc=$?
fi
# rank 0's own exit status (srun reports the largest task status; a task killed after rank 0's exit 75 is not a failure)
R0=$($PY -c 'import json,sys; e=json.load(open(sys.argv[1])).get("exit"); print(e if isinstance(e,int) else "")' "$P/ranks/rank000.end.json" 2>/dev/null)
echo "$(date -u +%FT%TZ) fit segment $SEG: launcher exit $rc, rank 0 exit ${R0:-unknown}"
[ -n "$R0" ] && rc=$R0
grep -E '^(clip_refit|positive_fit|consistency|initial cost|lambda |gpu_model: bound|gpu_model: basis|gpu_fit: |grain_fit: |fit data|selection profile|mpi_support: |Function evaluations)' "$P/fit.log" | cut -c1-400
tail -4 "$P/fit.log" | cut -c1-250
if [ "$TRANSPORT" = mpi ]; then         # tracebacks or errors on any rank
    grep -l -E 'Traceback|Error|MPI_ABORT|out of memory' "$P"/ranks/task-*.log 2>/dev/null | head -5 | while read -r f; do echo "== $f"; grep -m5 -E 'Traceback|Error|MPI_ABORT|out of memory' "$f" | cut -c1-300; done
fi
free -g > "$P/free_at_end.txt"; df -h "$D" > "$P/df_at_end.txt" 2>&1; du -sh "$OUT" >> "$P/df_at_end.txt" 2>&1
ST=failed
if [ -f "$OUT/fit_restart.npz" ]; then read -r _ ST _ < <(state "$D" "$SEG"); fi
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$SEG" "$SLURM_JOB_ID" "$T0" "$(date -u +%FT%TZ)" "$NN" "$NODE_LIST" \
    "$NR" "$TRANSPORT" "$ACTIVITY" "$MODEL_EXT" "$LEFT" "$ST" "$rc" >> "$D/segments.tsv"
if [ "$rc" -eq 0 ] && [ "$ST" = complete ]; then final_link "$SEG" || rc=$?; fi
[ "$rc" -eq 75 ] && echo "stopped before the wall-time limit: resubmit the same SPEC (NODES=1, 2 or auto) to continue with segment $((SEG + 1))"
exit "$rc"
