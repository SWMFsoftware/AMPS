#!/bin/bash
# ==== version-controlled copy: header added 2026-11-17, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_chain_post.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_chain_post.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: d6d0e94536b66bca376e8558f9cbdee7df531e8c196db96035d7b6b6c553fd86
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (CPU): the post-fit tools of a finished chain: degeneracy, report (last fit), tools with FT_DIR
#   (TOOLS_STEPS, default export contrib maps; 2026-11-17 the reruns without maps, ~15-20 GB each), chain docs;
#   POST_STEPS; disk guard (>= 20 GB per step, maps 45 GB; else the step is skipped, rc 4); exits 0 when the chain
#   is not complete.
# 2026-11-17: 2026-11-17: NPY_DISABLE_CPU_FEATURES exported for all steps (the bitwise A check of the tools).
# RUN 2026-11-17: 136123, 136097.
# ==== end of header ====
# Option C stage 2 (lead 2026-11-15): the post-fit tools of ONE finished localnn chain (jobs/localnn_chain.sbatch), CPU:
#   1. localnn_degeneracy.sbatch on <CHAIN>/c<CYCLES>       (m_b vs c_j diagnostics of the last fit)
#   2. localnn_report.sbatch     on <CHAIN>/c<CYCLES>       (summary numbers of the last fit, its README / HOW-TO-RERUN)
#   3. localnn_tools.sbatch      on <CHAIN>/c<CYCLES> with FT_DIR=<CHAIN>/c<CYCLES>/ft/final (maps and contributions of
#      the chain's final state = the last fine tuning: tools_state_ft, contributions_ft, maps_ft)
#   4. localnn_chain_docs.py     (CHAIN/README.md, HOW-TO-RERUN.md)
# The three job scripts run as plain bash inside this job (their provenance folders carry this job's id). The GPU
# cross-check of the final state (localnn_gpu_eval.sbatch with FT_DIR) is a separate GPU job.
# Exits 0 without doing anything when the chain is not complete (no chain_complete event).
# 2026-11-17: TOOLS_STEPS (default "export contrib maps"; the maps are ~15-20 GB per chain, "export contrib" ~2.2 GB) and a
#   disk guard (each step needs >= DISK_MIN_GB = 20 GB free, the maps step 20 GB + 25 GB; otherwise the step is skipped,
#   rc 4).
# Settings: CHAIN (required), CYCLES (5), MODEL_EXT (the frozen localnn snapshot), POST_STEPS (default "degeneracy report
#   tools docs"; 2026-11-15: renamed from STEPS, which leaked into localnn_tools.sbatch's own STEPS, job 132812).
#SBATCH --job-name=cops_localnn_post
#SBATCH --partition=normal128cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=320G
#SBATCH --time=05:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_localnn_post-%j.out
set -uo pipefail
# numpy without AVX-512 kernels, as the chain driver (2026-11-17: on an AVX-512 node (node028, EPYC 9554) np.exp differs
# in the last bit at 38 of the 271 held facets and the tools' bitwise reproduction of A failed, post 136094)
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
J=$COPS_OUTPUT_ROOT/jobs
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
CHAIN=$(readlink -f "${CHAIN:?CHAIN is required}"); CYCLES=${CYCLES:-5}
POST_STEPS=${POST_STEPS:-degeneracy report tools docs}
unset STEPS
TOOLS_STEPS=${TOOLS_STEPS:-export contrib maps}; DISK_MIN_GB=${DISK_MIN_GB:-20}
freegb() { df -P -B1G "$CHAIN" | awk 'NR==2 {print $4}'; }
PY=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda/.venv/bin/python
stamp() { echo "$(date -u +%FT%TZ) $*"; }
stamp "job $SLURM_JOB_ID on $(hostname): CHAIN=$CHAIN CYCLES=$CYCLES POST_STEPS='$POST_STEPS' MODEL_EXT=$MODEL_EXT"
if ! grep -q '"stage": "chain_complete", "rc": 0' "$CHAIN/chain_events.jsonl" 2>/dev/null; then
  stamp "chain not complete: nothing done"; "$PY" "$MODEL_EXT/localnn_chain_docs.py" --chain "$CHAIN"; exit 0; fi
RR=$CHAIN/c$CYCLES
rc=0
for s in $POST_STEPS; do
  need=$DISK_MIN_GB; [ "$s" = tools ] && case " $TOOLS_STEPS " in *" maps "*) need=$(( DISK_MIN_GB + 25 ));; esac
  fg=$(freegb); if [ -z "$fg" ] || [ "$fg" -lt "$need" ]; then stamp "step $s skipped: ${fg:-?} GB free < $need GB (disk guard)"; rc=4; continue; fi
  case $s in
    degeneracy) RUNROOT=$RR bash "$J/localnn_degeneracy.sbatch"; r=$?;;
    report)     RUNROOT=$RR bash "$J/localnn_report.sbatch"; r=$?;;
    tools)      RUNROOT=$RR FT_DIR=$RR/ft/final STEPS="$TOOLS_STEPS" bash "$J/localnn_tools.sbatch"; r=$?;;
    docs)       "$PY" "$MODEL_EXT/localnn_chain_docs.py" --chain "$CHAIN"; r=$?;;
    *) stamp "unknown step $s"; r=2;;
  esac
  stamp "step $s rc $r"; [ $r -eq 0 ] || rc=$r
done
printf '{"stage": "post", "rc": %s, "device": "cpu", "host": "%s", "job": "%s", "utc": "%s", "note": "%s"}\n' \
  "$rc" "$(hostname)" "$SLURM_JOB_ID" "$(date -u +%FT%TZ)" "$POST_STEPS" >> "$CHAIN/chain_events.jsonl"
"$PY" "$MODEL_EXT/localnn_chain_docs.py" --chain "$CHAIN" >/dev/null 2>&1
mkdir -p "$CHAIN/provenance-$SLURM_JOB_ID-post" && cp -p "$0" "$CHAIN/provenance-$SLURM_JOB_ID-post/localnn_chain_post.sbatch.used"
cp -p "$J/logs/cops_localnn_post-$SLURM_JOB_ID.out" "$CHAIN/provenance-$SLURM_JOB_ID-post/job.log" 2>/dev/null
stamp "done rc $rc"
exit $rc
