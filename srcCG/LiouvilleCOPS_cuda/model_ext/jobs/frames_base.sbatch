#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/frames_base.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/frames_base.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 8318598df235b534c38a9f53a69cdc1c7a369d639721aa54d54e4515b5e52a53
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Fitting-base kernels at the movie frame times: the lategrid base's boundary.in with ObservationFile /
#   OutputDirectory / LastPoint / Threads (FirstPoint, constant T_RG after the HK table: FIRST, TRG_CONST) changed, the
#   base's executable (sha256 checked), one GPU.
# RUN: sbatch --nodelist=gpu002 --export=ALL (FB, OBSF[, FIRST, TRG_CONST])
# ==== end of header ====
# Fitting-base kernels at the movie frame times (2026-11-14, user: "colour = contribution per unit surface area, exactly
# as the fitted model computes it"): the 12M lategrid base's boundary.in with only ObservationFile (the frame times),
# OutputDirectory, LastPoint and Threads changed (checked), run with the same executable native/liouville_cops_cuda
# (sha256 of the base) on one GPU -> <run>/frames_base_pan1959 (kernels/point=NNNN.bin = frame NNNN, COPS_model.dat).
# Environment: OBSF (observation file), FB (output directory); FIRST (FirstPoint, default 0) and TRG_CONST (2026-11-14:
# TemperatureRGMode = constant with TemperatureRG = TRG_CONST [K], for frames after the end of the RG temperature table;
# T_RG enters only the RG pressure conversion, not the kernels) - default: unchanged.
#SBATCH --job-name=cops_frames_base
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=96G
#SBATCH --time=04:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_frames_base-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
RUN=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; BASE=$RUN/base; FB=${FB:?}; OBSF=${OBSF:?}
THREADS=${SLURM_CPUS_PER_TASK:-48}; FIRST=${FIRST:-0}; TRG_CONST=${TRG_CONST:-}
BIN=$APP/native/liouville_cops_cuda; BIN_SHA=34a7551fc5d1381d3c2d5e59f7ab56a2cbbdb8c0222811f126bc72039c0dc0bd
now() { date -u +%FT%TZ; }
[ "$(sha256sum < "$BIN" | cut -d' ' -f1)" = "$BIN_SHA" ] || { echo "liouville_cops_cuda differs from the executable of the 12M base" >&2; exit 2; }
[ -e "$FB" ] && { echo "$FB exists" >&2; exit 2; }
mkdir -p "$FB"; P=$FB.provenance-$SLURM_JOB_ID; mkdir "$P" || exit 1
cp -p "$0" "$P/frames_base.sbatch.used"; cp -p "$OBSF" "$P/"
sha256sum "$BIN" "$OBSF" > "$P/executables.sha256"
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (CUDA from $CUDA_HOME)" >> "$P/module_list.txt"
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
NP=$(grep -c '^"2016' "$OBSF")
SED=(sed -e "s|^ObservationFile *=.*|ObservationFile = $OBSF|" -e "s|^OutputDirectory *=.*|OutputDirectory = $FB|"
     -e "s|^LastPoint *=.*|LastPoint  = $((NP - 1))|" -e "s|^Threads *=.*|Threads = $THREADS|" -e "s|^FirstPoint *=.*|FirstPoint = $FIRST|")
[ -n "$TRG_CONST" ] && SED+=(-e "s|^TemperatureRGMode *=.*|TemperatureRGMode = constant      # frames_base.sbatch TRG_CONST: after the RG temperature table|"
                             -e "s|^TemperatureRG  *=.*|TemperatureRG  = $TRG_CONST            # frames_base.sbatch TRG_CONST: the table's last measured value, held|"
                             -e "s|^\(TemperatureRGTable.*\)|# frames_base.sbatch TRG_CONST (not read in constant mode): \1|")
SED+=("$BASE/boundary.in")
"${SED[@]}" > "$FB/boundary.in"
diff "$BASE/boundary.in" "$FB/boundary.in" > "$FB/boundary.in.diff"
KEYS='^(ObservationFile|OutputDirectory|LastPoint|Threads|FirstPoint|TemperatureRGMode|TemperatureRG) *=|TemperatureRGTable'
diff <(grep -vE "$KEYS" "$BASE/boundary.in") <(grep -vE "$KEYS" "$FB/boundary.in") > /dev/null || { echo "boundary.in differs outside the intended keys" >&2; exit 4; }
for want in "ObservationFile = $OBSF" "OutputDirectory = $FB" "LastPoint  = $((NP - 1))" "Threads = $THREADS" "FirstPoint = $FIRST"; do
    [ "$(grep -cxF "$want" "$FB/boundary.in")" -eq 1 ] || { echo "boundary.in: missing '$want'" >&2; exit 4; }
done
echo "$(now) job $SLURM_JOB_ID on $(hostname): $NP points; boundary.in changes:"; cat "$FB/boundary.in.diff"
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $AMPS/srcCG/LiouvilleCOPS"; echo "COPS_CUDA_DEVICE=0 $BIN $FB/boundary.in > $FB/build.log 2>&1"; } > "$P/commands.sh"
cd "$AMPS/srcCG/LiouvilleCOPS" || exit 1
COPS_CUDA_DEVICE=0 bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "$BIN" "$FB/boundary.in" > "$FB/build.log" 2>&1; rc=$?
echo "$(now) build exit $rc"; tail -4 "$FB/build.log" | cut -c1-220; du -sh "$FB"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_frames_base-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
