#!/bin/bash
# ==== version-controlled copy: header added 2026-11-07, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/late_response.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/late_response.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 82d1c0abfe02d62452506799b7614b9b87bb8f78fee55397246171039bb2450c
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job. Late-window forward check (user 2026-11-07): 1-s time grid, one base per gauge offset along +Z
#   (default 0 1.5 3.0 m) built in parallel on 3 GPUs with the unchanged native/liouville_cops_cuda, then
#   late_response_eval.py for cycle 3, cycle 3 + last 5 x10, cycle 3 + last 6 ramp, original-T cycle 6 and its ramp
#   fit (fit and fine-tuned). First run: job 122289.
# RUN: sbatch late_response.sbatch   [--export=ALL,DZS="0 1.5 3.0",THREADS=15]
#   products runs/shap8_012m_flat_earlier/late_response/<job>/ (late_grid.h, base_dz*, eval/, provenance/commands.sh).
# ==== end of header ====
# Forward-only check of the late window (user 2026-11-07): the COPS time response and the gauge position.
#   "1: both" (compare the model with the measurement's 10-s running average and with the 10 s + DPU 5-readout average,
#   ROSINA COPS instrument modes RO-ROS-MAN-1019 sect. 4.2); "3: make the sensitivity check" (gauge offset along +Z of
#   the spacecraft: the COPS position relative to the trajectory point (centre of gravity) is not known; 0, 1.5, 3 m);
#   "for the forward run use two 12M best fits with and without extra weight of the last point observations".
# 1. model_ext/late_time_grid.py: observation file with a 1-s grid 10:32:46.282 .. 10:39:26.282 (the last L2 STOP) plus
#    the 95 original observations from 10:33:10 (reference values unchanged; grid points are gaps).
# 2. One late-window base per offset with the unchanged GPU forward model native/liouville_cops_cuda (the fitting base's
#    boundary.in with only ObservationFile, OutputDirectory, LastPoint, GaugeOffsetSC = 0 0 DZ and Threads changed),
#    one GPU each, in parallel.
# 3. model_ext/late_response_eval.py: the fits' and fine tunings' per-facet activity, temperature law, kappa and grain
#    on each base; windows inst / W4 / W10 / W10x5; validation of the dz0 base at the original times against each
#    fit's own values (must agree to 1e-4).
# Models: cycle 3 (no extra weight; lowest-cost 12M fit), cycle 3 + last 5 x10 (122236), cycle 3 + last 6 ramp 1->20
# (122248), original-T cycle 6 (122221), original-T cycle 6 + ramp (122250); each as fit and fine-tuned.
#   sbatch late_response.sbatch          (DZS="0 1.5 3.0" by default; products runs/<12M>/late_response/<job>/)
#SBATCH --job-name=cops12m_lresp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=3
#SBATCH --cpus-per-task=48
#SBATCH --mem=150G
#SBATCH --time=03:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_lresp-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
DZS=${DZS:-0 1.5 3.0}; THREADS=${THREADS:-15}
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; X=$COPS_OUTPUT_ROOT/model_ext; PY="$APP/.venv/bin/python -W ignore"
BIN=$APP/native/liouville_cops_cuda; G=$R/grain; LW=$G/late_weights
D=$R/late_response/$SLURM_JOB_ID; P=$D/provenance; [ -e "$D" ] && { echo "exists: $D" >&2; exit 2; }
mkdir -p "$P"; cd "$AMPS/srcCG/LiouvilleCOPS" || exit 1     # the base was built from here (rebuild.json cwd)
cp -p "$0" "$P/late_response.sbatch"; cp -p "$X/late_time_grid.py" "$X/late_response_eval.py" "$P/"
sha256sum "$X"/*.py "$BIN" > "$P/model_ext.sha256"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z_]*|DZS|THREADS|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH|SLURM_JOB_GPUS|CUDA_VISIBLE_DEVICES)=' | sort > "$P/environment.txt"
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); AMPS $(git -C "$AMPS" rev-parse HEAD)"; grep '^export ' "$0"; echo "cd $AMPS/srcCG/LiouvilleCOPS"; } > "$P/commands.sh"
rec() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; }
free -g > "$P/free_at_start.txt"
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then NGPU=$(echo "$CUDA_VISIBLE_DEVICES" | tr , "\n" | wc -l); else NGPU=$(nvidia-smi -L | wc -l); fi
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): $NGPU GPUs visible (CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-unset}); offsets $DZS"
# 1. time grid
rec $PY "$X/late_time_grid.py" --observations "$AMPS/srcCG/RosinaMeasurements.h" --hk "$COPS_OUTPUT_ROOT/data/COPS_HK_20160930.txt" \
    --start 2016-09-30T10:32:46.282 --stop 2016-09-30T10:39:26.282 --step 1 --keep-from 2016-09-30T10:33:10 \
    --out "$D/late_grid.h" --table "$D/late_grid_points.txt"
$PY "$X/late_time_grid.py" --observations "$AMPS/srcCG/RosinaMeasurements.h" --hk "$COPS_OUTPUT_ROOT/data/COPS_HK_20160930.txt" \
    --start 2016-09-30T10:32:46.282 --stop 2016-09-30T10:39:26.282 --step 1 --keep-from 2016-09-30T10:33:10 \
    --out "$D/late_grid.h" --table "$D/late_grid_points.txt" || exit 3
NP=$(grep -vc '^#' "$D/late_grid_points.txt")
# 2. bases, one GPU each
i=0; pids=()
for dz in $DZS; do
    B=$D/base_dz$dz; mkdir -p "$B"
    sed -e "s|^ObservationFile *=.*|ObservationFile = $D/late_grid.h|" \
        -e "s|^OutputDirectory *=.*|OutputDirectory = $B|" \
        -e "s|^FirstPoint *=.*|FirstPoint = 0|" -e "s|^LastPoint *=.*|LastPoint  = $((NP - 1))|" \
        -e "s|^GaugeOffsetSC *=.*|GaugeOffsetSC  = 0 0 $dz            # late_response.sbatch: sensitivity offset along +Z (s/c frame) [m]|" \
        -e "s|^Threads *=.*|Threads = $THREADS|" "$R/base/boundary.in" > "$B/boundary.in"
    diff "$R/base/boundary.in" "$B/boundary.in" > "$B/boundary.in.diff"
    [ "$(grep -c '^[<>]' "$B/boundary.in.diff")" -eq 10 ] || { echo "unexpected boundary.in changes:"; cat "$B/boundary.in.diff"; exit 4; }
    echo "env COPS_CUDA_DEVICE=$i bash -c 'echo 1000 > /proc/self/oom_score_adj; exec \"\$@\"' _ $BIN $B/boundary.in" >> "$P/commands.sh"
    COPS_CUDA_DEVICE=$i bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "$BIN" "$B/boundary.in" > "$B/build.log" 2>&1 &
    pids+=($!); i=$((i + 1))
    echo "$(date -u +%FT%TZ) base dz=$dz: started on device $((i - 1)) (pid ${pids[-1]})"
done
rc=0; for p in "${pids[@]}"; do wait "$p" || rc=1; done
for dz in $DZS; do echo "--- dz=$dz"; tail -3 "$D/base_dz$dz/build.log" | cut -c1-200; done
[ $rc -eq 0 ] || { echo "$(date -u +%FT%TZ) a base build failed" >&2; exit 5; }
echo "$(date -u +%FT%TZ) bases done"; du -sh "$D"/base_dz* | sed 's/^/  /'
# 3. evaluation (the label dz0 must be the zero offset: it is validated against the fits)
M=(--model "c3=$G/cycles/c3/C30_n4000_grain_cross10_ftc3,$G/cycles/c3/finetune/C30_n4000_grain_cross10_ftc3_m5_i5_cpu")
for n in c3_last5_w10_constant_ngrg c3_last6_w1to20_linear_ngrg; do
    M+=(--model "${n%%_ngrg}=$LW/C30_n4000_grain_cross10_$n/fit,$LW/C30_n4000_grain_cross10_$n/finetune/C30_n4000_grain_cross10_${n}_m5_i5_cpu")
done
M+=(--model "torig_c6=$G/cycles_Torig/c6/C30_n4000_Torig_grain_cross10_ftc6,$G/cycles_Torig/c6/finetune/C30_n4000_Torig_grain_cross10_ftc6_m5_i5_cpu")
M+=(--model "torig_c6_last6_w1to20_linear=$LW/C30_n4000_grain_cross10_torigc6_last6_w1to20_linear_ngrg/fit,$LW/C30_n4000_grain_cross10_torigc6_last6_w1to20_linear_ngrg/finetune/C30_n4000_grain_cross10_torigc6_last6_w1to20_linear_ngrg_m5_i5_cpu")
EV=(--points "$D/late_grid_points.txt")
for dz in $DZS; do l=dz$dz; [ "$dz" = 0 ] && l=dz0; EV+=(--base "$l=$D/base_dz$dz"); done
rec $PY "$X/late_response_eval.py" "${EV[@]}" "${M[@]}" --out "$D/eval"
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ $PY "$X/late_response_eval.py" "${EV[@]}" "${M[@]}" --out "$D/eval" > "$P/eval.log" 2>&1; rc=$?
echo "$(date -u +%FT%TZ) evaluation: exit $rc"; tail -45 "$P/eval.log" | cut -c1-220
exit $rc
