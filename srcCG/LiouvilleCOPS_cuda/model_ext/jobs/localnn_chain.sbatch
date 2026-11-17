#!/bin/bash
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_chain.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_chain.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 87be390c9e61576bfeae766625fc8afefb6221ef36c672de9c88a875f1aa2bff
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (CPU or 1 GPU, option C stage 2): ONE 5-cycle fit / fine-tuning chain (start from a state + grain
#   + held, fit with checkpoints and the wall-time guard (exit 75), direct check, CPU fine tuning (late_only,
#   WINDOW_MIN 5, MAXNFEV 500, validation), links); chain_events.jsonl (device per stage); HANDOVER file; DEVICE=gpu
#   (CUDA environment, takeover with CANCEL_JOBS; 2026-11-17 TAKEOVER_AFTER_CHECKPOINT=1: the takeover cancels the
#   CPU attempt only seconds after the active fit's next checkpoint write, or when the CPU attempts have ended);
#   TEMPLATES, C1_START_EXTRA, FT_EXTRA, NONNEG_CODE; disk guard DISK_MIN_GB (20) before every stage, exit 4.
# RUN 2026-11-15: stage2r f10 / f1 pm1e-4 GPU attempts 136092-136093, 136095-136096 (gpu002).
# ==== end of header ====
# Option C stage 2 (user / lead 2026-11-15): ONE fit / fine-tuning chain of the local non-negative basis, cycles 1..CYCLES,
# the layout of the nonneg chains: <CHAIN>/c<k>/fit (localnn_fit.py: segNN, final; checkpoints) and <CHAIN>/c<k>/ft
# (localnn_ft.py: the frozen nonneg_finetune.run on CPU: seg01, final + links); <CHAIN>/chain_events.jsonl (one line per
# stage: stage, rc, device, host, job, UTC); <CHAIN>/provenance-<job>/ per attempt.
#   cycle 1 fit: from the state C1_FROM (a fit's localnn_state.npz: the stage-1 end point), or C1_FIT = an existing fit
#                directory linked as c1/fit (it must have final);
#   cycle k > 1 fit: from c<k-1>/fit/final's state, grain A, B from c<k-1>/ft/final, held h = c<k-1>/ft/final held_next;
#   fine tuning of every cycle: late_only facets, WINDOW_MIN 5, no lower bound, grain refitted, lambda 0, MAXNFEV 500,
#                validation checks 1-8, the softplus chains' optimizer settings (lsmr 50, ftol/xtol/gtol 1e-10).
# The job runs the stages in order and skips complete ones; a fit that reaches the wall-time guard exits 75 (the next
# attempt continues from its segment or checkpoint). Submit attempts with jobs/localnn_chain_submit.sh (afternotok chain).
# Handover to another device: a file <CHAIN>/HANDOVER containing a device name (e.g. "gpu"): a job of another DEVICE
# stops before its next stage with exit 75 (event "handover"); the job of that DEVICE continues the chain (a fit in
# progress continues from its checkpoint / segment).
# Settings: CHAIN (required), LAMBDA (required, exact repr) + RULE (the rule file), FM (multiplier penalty f_m),
#   C1_FROM or C1_FIT, CYCLES (5), NFEV (4000), FT_MAXNFEV (500), WINDOW_MIN (5), CHECKPOINT (250), LATE (late kernels),
#   DEVICE (cpu; gpu = the GPU solver path when available), MODEL_EXT (frozen snapshot), BLAS_THREADS (16), GRID_DIR,
#   CHECK (1: localnn_fit --check-jacobian at the first segment of every cycle's fit), GPU_SVD_MODE (qr; DEVICE=gpu).
#   DEVICE=gpu: submit with the GPU options on the sbatch line (localnn_chain_submit.sh --device gpu --sbatch-extra
#   "--partition=gpu --account=general --gpus-per-node=1 ..."); the fine tuning always runs on the CPU (seconds).
#SBATCH --job-name=cops_localnn_chain
#SBATCH --partition=normal128cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_localnn_chain-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
BLAS_THREADS=${BLAS_THREADS:-16}
export OMP_NUM_THREADS=$BLAS_THREADS OPENBLAS_NUM_THREADS=$BLAS_THREADS MKL_NUM_THREADS=$BLAS_THREADS PYTHONUNBUFFERED=1
# numpy without AVX-512 kernels on every node type (the stages reproduce A bitwise across attempts: np.exp differs
# between the AVX2 and AVX-512 loops; as jobs/nonneg_contrib.sbatch)
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
if [ "${DEVICE:-cpu}" = gpu ]; then   # the GPU solver path (localnn_gpu_solver.py: CuPy FP64, the CUDA 12.9 toolkit)
  export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=/data001/highresolu_vt/cops_output/.cupy
  export PATH=/usr/local/cuda-12.9/bin:$PATH LD_LIBRARY_PATH=/usr/local/cuda-12.9/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
fi
UNSET=$(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*|COPS_NONNEG_[A-Z_]*)' | sort | paste -sd' ' -)
[ -n "$UNSET" ] && unset $UNSET
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
BASE=$R/base; D2=$R/clip_refit/D2/final; PT=$R/late_grid_points.txt
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
GRID_DIR=${GRID_DIR:-$R/nonneg/localnn/grid_150m}
CHAIN=$(readlink -f "${CHAIN:?CHAIN is required}")
LAMBDA=${LAMBDA:?LAMBDA is required}; RULE=$(readlink -f "${RULE:?RULE is required}")
FM=${FM:?FM is required}; CYCLES=${CYCLES:-5}; NFEV=${NFEV:-4000}; FT_MAXNFEV=${FT_MAXNFEV:-500}
WINDOW_MIN=${WINDOW_MIN:-5}; CHECKPOINT=${CHECKPOINT:-250}; DEVICE=${DEVICE:-cpu}; CHECK=${CHECK:-1}
LATE=$(readlink -f "${LATE:-$R/nonneg/localnn/late_kernels_150m.npz}")
C1_FROM=${C1_FROM:-}; C1_FIT=${C1_FIT:-}
# 2026-11-16 (user / lead: the stage-2 reruns): TEMPLATES (on | off; D2's templates in the model), C1_START_EXTRA (extra
# arguments of the cycle-1 start, e.g. --state-drop-templates), FT_EXTRA (extra arguments of every fine tuning, e.g. the
# regularization --delta-penalty-f 10 --delta-cref 4.2688066575e-03 --delta-cap 4.605170185988092), NONNEG_CODE (the
# frozen nonneg_finetune snapshot of localnn_ft.py: LOCALNN_NONNEG_CODE)
TEMPLATES=${TEMPLATES:-on}; C1_START_EXTRA=${C1_START_EXTRA:-}; FT_EXTRA=${FT_EXTRA:-}
[ -n "${NONNEG_CODE:-}" ] && export LOCALNN_NONNEG_CODE=$(readlink -f "$NONNEG_CODE")
CONTRIB=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_contrib_20261112T210012Z
PY=$APP/.venv/bin/python
stamp() { echo "$(date -u +%FT%TZ) $*"; }
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
LIMIT=$(squeue -h -j "$SLURM_JOB_ID" -o %l 2>/dev/null)
secs() { local t=$1 d=0; [[ $t == *-* ]] && { d=${t%%-*}; t=${t#*-}; }; IFS=: read -r a b c <<< "$t"; [ -z "$c" ] && { c=$b; b=$a; a=0; }; echo $(( d*86400 + 10#$a*3600 + 10#$b*60 + 10#$c )); }
DEADLINE=$(( $(date +%s) + $(secs "${LIMIT:-12:00:00}") - 1200 ))
mkdir -p "$CHAIN" || exit 1
P=$CHAIN/provenance-$SLURM_JOB_ID; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): CHAIN=$CHAIN LAMBDA=$LAMBDA FM=$FM CYCLES=$CYCLES NFEV=$NFEV DEVICE=$DEVICE MODEL_EXT=$MODEL_EXT"
cp -p "$0" "$P/localnn_chain.sbatch.used"; cp -p "$MODEL_EXT"/*.py "$P/" 2>/dev/null
[ -f "$MODEL_EXT/SHA256SUMS" ] && cp -p "$MODEL_EXT/SHA256SUMS" "$P/MODEL_EXT_SHA256SUMS"
( cd "$P" && sha256sum *.py localnn_chain.sbatch.used > code_sha256.txt )
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded)" >> "$P/module_list.txt"
{ command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -q; } > "$P/nvidia_smi_q.txt" 2>&1 || echo "CPU job on $(hostname): no GPU used" >> "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$PY" > "$P/pip_list.txt" 2>&1 || "$PY" -m pip list > "$P/pip_list.txt" 2>&1
env | grep -E '^(SLURM_JOB_ID|SLURM_JOB_NODELIST|SLURM_CPUS_PER_TASK|MODEL_EXT|OMP_NUM_THREADS|OPENBLAS_NUM_THREADS|CUDA_VISIBLE_DEVICES)=' > "$P/env.txt"
echo "CHAIN=$CHAIN LAMBDA=$LAMBDA RULE=$RULE FM=$FM C1_FROM=$C1_FROM C1_FIT=$C1_FIT CYCLES=$CYCLES NFEV=$NFEV FT_MAXNFEV=$FT_MAXNFEV WINDOW_MIN=$WINDOW_MIN CHECKPOINT=$CHECKPOINT LATE=$LATE DEVICE=$DEVICE GRID_DIR=$GRID_DIR CHECK=$CHECK TEMPLATES=$TEMPLATES C1_START_EXTRA=$C1_START_EXTRA FT_EXTRA=$FT_EXTRA LOCALNN_NONNEG_CODE=${LOCALNN_NONNEG_CODE:-}" >> "$P/env.txt"
{ echo '#!/bin/bash'; grep '^export ' "$0" | sed "s/\$BLAS_THREADS/$BLAS_THREADS/g"; echo "export MODEL_EXT=$MODEL_EXT"; echo "cd $APP"
  echo "heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec \"\$@\"' heavy \"\$@\"; }"; } > "$P/commands.sh"
[ -n "${LOCALNN_NONNEG_CODE:-}" ] && echo "export LOCALNN_NONNEG_CODE=$LOCALNN_NONNEG_CODE" >> "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; "$@"; }
event() { printf '{"stage": "%s", "rc": %s, "device": "%s", "host": "%s", "job": "%s", "utc": "%s", "note": "%s"}\n' \
  "$1" "$2" "$DEVICE" "$(hostname)" "$SLURM_JOB_ID" "$(date -u +%FT%TZ)" "${3:-}" >> "$CHAIN/chain_events.jsonl"; }
COM=(--design "$GRID_DIR/design" --base "$BASE" --d2 "$D2" --point-table "$PT")
DEVOPT=(); [ "$DEVICE" != cpu ] && DEVOPT=(--device "$DEVICE" --gpu-svd-mode "${GPU_SVD_MODE:-qr}")
# disk guard (2026-11-17, after /data001/highresolu_vt filled at 2026-11-17 03:10Z): before every stage the free space of
# the chain's file system must be >= DISK_MIN_GB (20); otherwise the job stops cleanly with exit 4 (event "disk_guard";
# a fit in progress continues from its checkpoint when resubmitted)
DISK_MIN_GB=${DISK_MIN_GB:-20}
diskguard() {
  local free; free=$(df -P -B1G "$CHAIN" | awk 'NR==2 {print $4}')
  if [ -z "$free" ] || [ "$free" -lt "$DISK_MIN_GB" ]; then
    stamp "disk guard before $1: ${free:-?} GB free < $DISK_MIN_GB GB: exit 4"; event disk_guard 4 "before $1: ${free:-?} GB free"; exit 4; fi; }
handover() {   # exit 75 before the next stage when <CHAIN>/HANDOVER names another device
  if [ -s "$CHAIN/HANDOVER" ] && [ "$(head -1 "$CHAIN/HANDOVER" | tr -d ' \n')" != "$DEVICE" ]; then
    stamp "handover to $(head -1 "$CHAIN/HANDOVER") before $1: exit 75"; event handover 75 "before $1"; exit 75; fi; }
# device takeover (2026-11-15, lead "option 3"): a job of DEVICE != cpu submitted with CANCEL_JOBS (the chain's CPU
# attempts and its old post job, all submitted by this session) writes HANDOVER = DEVICE, cancels them and waits until
# they have left the queue, then continues the chain (a fit resumes from its checkpoint / segment); skipped when the
# chain is already complete
if [ "$DEVICE" != cpu ] && [ -n "${CANCEL_JOBS:-}" ] && ! grep -q '"stage": "chain_complete", "rc": 0' "$CHAIN/chain_events.jsonl" 2>/dev/null; then
  echo "$DEVICE" > "$CHAIN/HANDOVER"
  alive() { for jb in "$@"; do squeue -h -o %i -j "$jb" 2>/dev/null; done | paste -sd' ' -; }   # one id per call
  LIVE=$(alive $(echo "$CANCEL_JOBS" | tr ":," "  "))
  # 2026-11-17 (lead): TAKEOVER_AFTER_CHECKPOINT=1: do not cancel a running CPU fit at an arbitrary moment: wait until
  # the active fit (the first cycle without fit/final) writes its next checkpoint.npz (atomic os.replace; the cancel
  # then comes seconds after the write, nothing in progress) or until the CPU attempts have left the queue (e.g. exit 75
  # at the HANDOVER check); at most TAKEOVER_WAIT_MIN (90) min, then cancel anyway
  if [ -n "$LIVE" ] && [ "${TAKEOVER_AFTER_CHECKPOINT:-0}" = 1 ]; then
    AF=; for kk in $(seq 1 "$CYCLES"); do [ -e "$CHAIN/c$kk/fit/final" ] || { AF=$CHAIN/c$kk/fit; break; }; done
    m0=$(stat -c %Y "$AF/checkpoint.npz" 2>/dev/null || echo 0)
    stamp "takeover: waiting for the next checkpoint of $AF (now $(date -u -d @$m0 +%FT%TZ 2>/dev/null)) or the end of $LIVE"
    for i in $(seq 1 $(( ${TAKEOVER_WAIT_MIN:-90} * 6 ))); do
      head -c 1 "$AF/checkpoint.npz" >/dev/null 2>&1   # open(): NFS close-to-open revalidation (acregmin 600)
      m1=$(stat -c %Y "$AF/checkpoint.npz" 2>/dev/null || echo 0)
      { [ "$m1" != "$m0" ] && [ "$m1" != 0 ]; } && { sleep 5; break; }
      [ -z "$(alive $LIVE)" ] && break
      sleep 10
    done
    event takeover_wait 0 "checkpoint $(date -u -d @$(stat -c %Y "$AF/checkpoint.npz" 2>/dev/null || echo 0) +%FT%TZ)"
    LIVE=$(alive $LIVE)
  fi
  if [ -n "$LIVE" ]; then
    scancel $LIVE
    for jb in $LIVE; do echo "$(date -u +%FT%TZ) $jb scancel (takeover of chain $CHAIN by $DEVICE job $SLURM_JOB_ID)" >> "$COPS_OUTPUT_ROOT/model_ext/logs/submissions.txt"; done
    for i in $(seq 1 90); do [ -z "$(alive $LIVE)" ] && break; sleep 10; done
  fi
  event takeover 0 "cancelled: ${LIVE:-none}"; stamp "takeover by $DEVICE: cancelled ${LIVE:-none}"
fi
if [ -n "$C1_FIT" ] && [ ! -e "$CHAIN/c1/fit" ]; then mkdir -p "$CHAIN/c1"; ln -s "$(readlink -f "$C1_FIT")" "$CHAIN/c1/fit"; event c1fit_linked 0 "$C1_FIT"; fi
for k in $(seq 1 "$CYCLES"); do
  FIT=$CHAIN/c$k/fit; FT=$CHAIN/c$k/ft
  if [ ! -e "$FIT/final" ]; then
    handover "c${k}fit"; diskguard "c${k}fit"
    mkdir -p "$FIT"
    HELD=();
    if [ "$k" -gt 1 ]; then
      PF=$CHAIN/c$((k-1)); PFT=$(readlink -f "$PF/ft/final")
      HELD=(--held "$PFT/finetune_state.npz" --late-kernels "$LATE")
      SRC=(--state-file "$(readlink -f "$PF/fit/final")/localnn_state.npz" --grain-from "$PFT/finetune_state.npz")
    else
      [ -n "$C1_FROM" ] || { stamp "c1: neither C1_FROM nor C1_FIT"; exit 2; }
      SRC=(--state-file "$(readlink -f "$C1_FROM")")
      # shellcheck disable=SC2206
      [ -n "$C1_START_EXTRA" ] && SRC+=($C1_START_EXTRA)
    fi
    if [ ! -e "$FIT/start" ]; then
      stamp "c$k start"
      run heavy "$PY" -W ignore "$MODEL_EXT/localnn_fit.py" start "${COM[@]}" --grid "$GRID_DIR" --templates "$TEMPLATES" --start state \
        "${SRC[@]}" "${HELD[@]}" --out "$FIT/start" > "$P/c${k}_start.log" 2>&1; rc=$?
      event "c${k}start" $rc; [ $rc -eq 0 ] || { tail -5 "$P/c${k}_start.log"; exit $rc; }
    fi
    WALL=$(( DEADLINE - $(date +%s) ))
    [ $WALL -gt 1800 ] || { stamp "not enough time for c$k fit ($WALL s): exit 75"; event "c${k}fit" 75 "no time"; exit 75; }
    stamp "c$k fit (wall $WALL s)"
    JC=(); [ "$CHECK" = 1 ] && ! ls -d "$FIT"/seg[0-9][0-9] >/dev/null 2>&1 && [ ! -e "$FIT/checkpoint.npz" ] && JC=(--check-jacobian)
    run heavy "$PY" -W ignore "$MODEL_EXT/localnn_fit.py" fit "${COM[@]}" --grid "$GRID_DIR" --templates "$TEMPLATES" --slope-mult-min 1 \
      --start-file "$FIT/start" --smooth-lambda "$LAMBDA" --smooth-rule-file "$RULE" --mult-penalty-f "$FM" "${HELD[@]}" \
      --max-nfev "$NFEV" --wall-time "$WALL" --checkpoint-every "$CHECKPOINT" "${DEVOPT[@]}" "${JC[@]}" --out "$FIT" \
      >> "$P/c${k}_fit.log" 2>&1; rc=$?
    event "c${k}fit" $rc; tail -n 3 "$P/c${k}_fit.log" | cut -c1-300
    [ $rc -eq 75 ] && exit 75
    [ $rc -eq 0 ] || exit $rc
  fi
  # the direct facet-sum check of the fit's final segment (2026-11-15: outside the fit block, so that an attempt
  # stopped between the fit and its check, e.g. by a device takeover, runs it in the next attempt)
  FSEG=$(basename "$(readlink -f "$FIT/final")")
  if [ ! -e "$FIT/direct_check_$FSEG" ]; then
    diskguard "c${k}direct"
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 run heavy "$PY" -W ignore "$MODEL_EXT/localnn_direct.py" --seg "$(readlink -f "$FIT/final")" \
      --base "$BASE" --d2 "$D2" --point-table "$PT" --design "$GRID_DIR/design" --contrib-code "$CONTRIB" \
      --procs "${SLURM_CPUS_PER_TASK:-16}" > "$P/c${k}_direct.log" 2>&1; rc=$?
    event "c${k}direct" $rc; tail -n 1 "$P/c${k}_direct.log" | cut -c1-300; [ $rc -eq 0 ] || exit $rc
  fi
  if [ ! -e "$FT/final" ]; then
    handover "c${k}ft"; diskguard "c${k}ft"
    mkdir -p "$FT"
    if [ -e "$FT/seg01" ]; then mv "$FT/seg01" "$FT/seg01.incomplete_$(date -u +%Y%m%dT%H%M%SZ)"; fi
    stamp "c$k fine tuning"
    run heavy "$PY" -W ignore "$MODEL_EXT/localnn_ft.py" run --fit "$(readlink -f "$FIT/final")" --late "$LATE" \
      --design "$GRID_DIR/design" --base "$BASE" --d2 "$D2" --out "$FT/seg01" --facets late_only --window-min "$WINDOW_MIN" \
      --max-nfev "$FT_MAXNFEV" --lsmr-maxiter 50 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --x-scale-grain 1.0 --validate \
      --seed 20261110 --save-seconds 900 ${FT_EXTRA} > "$P/c${k}_ft.log" 2>&1; rc=$?
    event "c${k}ft" $rc; tail -n 2 "$P/c${k}_ft.log" | cut -c1-300
    [ $rc -eq 0 ] || exit $rc
    ln -s seg01 "$FT/final"
  fi
  for f in finetune_state.npz finetune_result.txt fine_tuning.json COPS_model.dat late_window_model.dat pressure_comparison.txt validation.json localnn_ft.json; do
    [ -e "$FT/final/$f" ] && [ ! -e "$FT/$f" ] && ln -s "final/$f" "$FT/$f"
  done
done
event chain_complete 0
stamp "chain complete"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_localnn_chain-$SLURM_JOB_ID.out" "$P/job.log" 2>/dev/null
exit 0
