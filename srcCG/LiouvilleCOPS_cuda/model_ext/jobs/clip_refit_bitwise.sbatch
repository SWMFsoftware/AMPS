#!/bin/bash
# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/clip_refit_bitwise.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/clip_refit_bitwise.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 03fba160b5d77aad2c1e183d08ac81370c124f1f24a10f9ef044e9527b61a8a0
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job (test). On c6's own base: grain_fit_fixed.py --restart-from c6 vs clip_refit.py --variant D1
#   --start-restart c6, compared bitwise (job 122854). 2026-11-08.
# ==== end of header ====
# GPU test (2026-11-08): with the W10x5 option off, model_ext/clip_refit.py adds nothing to the grain_fit_fixed.py path,
# bitwise. On c6's OWN base (runs/shap8_012m_flat_earlier/base, 1150 rows, no grid rows; record store, mixed
# precision, 8 GPUs, c6's command) two short fits from c6's final state, on the same node:
#   A  grain_fit_fixed.py ... --restart-from <c6 fit_restart.npz>   (gpu_fit's own restart mechanism, as c6 itself)
#   B  clip_refit.py --variant D1 --start-restart <c6 fit_restart.npz> --check-full-jacobian ...  (the D1 production
#      path: start restart by name, start check, full Jacobian check, records, restart marker)
# then model_ext/clip_compare.py fits A B (x, a, T, kappa, templates, activity_linear, pred_ng/rg, optimizer
# cost/nfev/status/message/optimality, COPS_model.dat, grain_term.txt, grain_components.dat, ActivityFit.txt,
# pressure_residuals.dat, fit_result.txt, manifest costs): exit 0 iff all identical. B's start check also gives c6's
# final cost reproduced on its own base (same problem) and the Jacobian check at c6's final state.
#   sbatch --nodelist=gpu00N --export=ALL,MODEL_EXT=<snapshot>[,NFEV=5] clip_refit_bitwise.sbatch
#SBATCH --job-name=cops12m_clipbit
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=8
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --time=01:00:00
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --exclude=gpu003
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops12m_clipbit-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:/usr/mpi/gcc/openmpi-4.1.7rc1/bin:$PATH
export LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
POSUNSET=$(env | grep -o '^COPS_POS_[A-Z_]*' | sort | paste -sd' ' -)
[ -n "$POSUNSET" ] && unset $POSUNSET
NFEV=${NFEV:-5}
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
E=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier; X=$MODEL_EXT; PY="$APP/.venv/bin/python -W ignore"
C6R=$E/grain/cycles_Torig/c6/C30_n4000_Torig_grain_cross10_ftc6/fit_restart.npz
TA=(--activity-template "$E/bundle/restart_templates.npz" cliff_ge_60deg)
FIX=(--fix Tmax=182.09922 --fix T75=136.11149 --fix Tmin=133.0)
for j in 1 2 3 4 5 6; do
    if [ $j -le 3 ]; then TD=$E/grain/cycles/c$j; else TD=$E/grain/cycles_Torig/c$j; fi
    TA+=(--activity-template "$TD/template_ft_c$j/ft_c$j.templates.npz" ft_c$j); FIX+=(--fix c_ft_c$j=1)
done
D=$E/clip_refit_tests/bitwise-$SLURM_JOB_ID; P=$D/provenance-$SLURM_JOB_ID
FILES=(clip_refit.py clip_problem.py clip_compare.py grain_fit_fixed.py fixed_parameters.py grain_fit.py grain_crossing.py
       positive_problem.py late_response_eval.py test_clip_refit.py test_positive_fit.py)
[ -e "$D" ] && { echo "exists: $D" >&2; exit 2; }
mkdir -p "$P"; cd "$APP" || exit 1
cp -p "$0" "$P/clip_refit_bitwise.sbatch"
for f in "${FILES[@]}"; do cp -p "$X/$f" "$P/"; done
sha256sum "$X"/*.py "$C6R" > "$P/model_ext.sha256" 2>&1
(cd "$P" && sha256sum "${FILES[@]}") > "$P/run_copy.sha256" 2>&1
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
module list > "$P/module_list.txt" 2>&1 || true
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; grep -m1 CSPICE /home/vtenishev/SPICE/cspice/doc/version.txt; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
env | grep -E '^(CUDA_[A-Z_]*|CUPY_[A-Z_]*|COPS_[A-Z_]*|NFEV|MODEL_EXT|OMP_[A-Z_]*|OPENBLAS_[A-Z_]*|MKL_[A-Z_]*|PYTHONUNBUFFERED|LD_LIBRARY_PATH|PATH|SLURM_JOB_GPUS|CUDA_VISIBLE_DEVICES)=' | sort > "$P/environment.txt"
{ echo '#!/bin/bash'; echo "# job $SLURM_JOB_ID on $(hostname), $(date -u +%FT%TZ); AMPS $(git -C "$AMPS" rev-parse HEAD)"
  grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
rec() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; }
free -g > "$P/free_at_start.txt"; df -h "$E" > "$P/df_at_start.txt" 2>&1
echo "model_ext copied and run: $X (MODEL_EXT); NFEV $NFEV${POSUNSET:+; unset $POSUNSET}" | tee "$P/segment.txt"
if [ -n "${CUDA_VISIBLE_DEVICES:-}" ]; then NG=$(echo "$CUDA_VISIBLE_DEVICES" | tr , "\n" | wc -l); else NG=$(nvidia-smi -L | wc -l); fi
DEV=$(seq -s, 0 $((NG - 1)))
COMMON=(--crossing-weight 10 --check-jacobian --devices "$DEV" --precision mixed --record-store "$E/store"
        --max-nfev "$NFEV" --wall-time slurm)
FITARGS=(--profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400
         --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6)
CA=($PY "$P/grain_fit_fixed.py" "${FIX[@]}" "${COMMON[@]}" --restart-from "$C6R" --template-start 0 -- "$E/base"
    "${FITARGS[@]}" --out "$D/A_grain_fit_fixed" "${TA[@]}")
CB=($PY "$P/clip_refit.py" --variant D1 --start-restart "$C6R" --check-full-jacobian "${FIX[@]}" "${COMMON[@]}"
    --template-start 0 -- "$E/base" "${FITARGS[@]}" --out "$D/B_clip_refit_D1" "${TA[@]}")
CC=($PY "$P/clip_compare.py" fits "$D/A_grain_fit_fixed" "$D/B_clip_refit_D1" --json "$D/compare.json")
rc=0
for k in A B; do
    if [ $k = A ]; then C=("${CA[@]}"); else C=("${CB[@]}"); fi
    rec bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${C[@]}"
    echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname): fit $k ($NG GPUs)"
    bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "${C[@]}" > "$P/fit_$k.log" 2>&1; r=$?
    echo "$(date -u +%FT%TZ) fit $k: exit $r"
    grep -E '^(clip_refit|consistency|initial cost|lambda |gpu_model: bound|gpu_fit: |grain_fit: |restart mapping|fit data)' "$P/fit_$k.log" | cut -c1-400
    [ $r -eq 0 ] || rc=$r
done
rec "${CC[@]}"
"${CC[@]}" 2>&1 | tee "$P/compare.log"; r=${PIPESTATUS[0]}
echo "$(date -u +%FT%TZ) compare: exit $r"; [ $r -eq 0 ] || [ $rc -ne 0 ] || rc=$r
df -h "$E" > "$P/df_at_end.txt" 2>&1; du -sh "$D" >> "$P/df_at_end.txt" 2>&1
exit $rc
