#!/bin/bash
# ==== version-controlled copy: header added 2026-11-16, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_20261114/ftreg/ftreg_runs.sh (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_20261114/ftreg/ftreg_runs.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 8a174816cac13fd1c75d29c3b01cb7e3ba655df1168a0b4979c8f1a1519bb371
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: the fine-tuning regularization test runs on CPU: c5 fine tuning of C f1 / 0.01, C f10 / 0.01 (localnn_ft.py
#   run), B and D (run-nonneg, plus an unregularized control each) with f_delta 0.1, 1, 10 (C_ref 4.2688066575e-03)
#   and the cap ln 100, MAXNFEV 500 (GROUP=C | BD; NONNEG_CODE = agent 1's snapshot nonneg_20261116T031115Z).
# RUN 2026-11-16: 134504 (C), 134505 (BD).
# ==== end of header ====
# fine-tuning regularization tests (lead / user 2026-11-15 ~21:10 CST): the c5 fine tuning re-run on CPU from the c5 fit
# state with P_delta = f_delta C_ref sum w (h+delta)^2 / sum w (f_delta 0.1, 1, 10; C_ref 4.2688066575e-03) and the cap
# |h+delta| <= ln 100, MAXNFEV 500, otherwise the chains' settings; no chain is rerun or changed.
#   GROUP=C : the two finished option C chains (localnn_ft.py run; LocalnnFitContext)
#   GROUP=BD: agent 3's B (background bg1e-2) and D (logspread f0.1) chains (localnn_ft.py run-nonneg: the real
#             NonnegFitContext of NONNEG_CODE; CPU engine with the fit-state pressures of late_kernels_BD.npz), plus
#             one UNREGULARIZED CPU control each (to compare with their GPU c5 fine tuning)
# Settings: GROUP, MODEL_EXT (frozen localnn snapshot), NONNEG_CODE (agent 1's frozen snapshot with the new options)
set -uo pipefail
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16 MKL_NUM_THREADS=16 PYTHONUNBUFFERED=1
W=/data001/highresolu_vt/cops_output/tests/localnn_20261114/ftreg
APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda; PY=$APP/.venv/bin/python; cd $APP
M=$(readlink -f "${MODEL_EXT:?MODEL_EXT}"); export LOCALNN_NONNEG_CODE=$(readlink -f "${NONNEG_CODE:?NONNEG_CODE}")
R=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5; N=$R/nonneg; L=$N/localnn
P=$W/provenance-$SLURM_JOB_ID-$GROUP; mkdir -p $P; cp -p "$0" $P/ftreg_runs.sh.used
cp -p $M/localnn_ft.py $M/localnn_late.py $M/localnn_problem.py $P/; mkdir $P/nonneg_code && cp -p $LOCALNN_NONNEG_CODE/*.py $P/nonneg_code/
( cd $P && sha256sum *.py nonneg_code/*.py > code_sha256.txt )
{ git -C /data001/highresolu_vt/AMPS rev-parse HEAD; git -C /data001/highresolu_vt/AMPS status --short; } > $P/git.txt 2>&1
echo "CPU job on $(hostname): no GPU used" > $P/nvidia_smi_q.txt; { gcc --version | head -1; lscpu | grep 'Model name'; } > $P/toolchain.txt
~/.local/bin/uv pip list --python $PY > $P/pip_list.txt 2>&1
echo "MODEL_EXT=$M LOCALNN_NONNEG_CODE=$LOCALNN_NONNEG_CODE GROUP=$GROUP" > $P/env.txt
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "export LOCALNN_NONNEG_CODE=$LOCALNN_NONNEG_CODE"; echo "cd $APP"; } > $P/commands.sh
COMMON=(--facets late_only --window-min 5 --max-nfev 500 --lsmr-maxiter 50 --ftol 1e-10 --xtol 1e-10 --gtol 1e-10
        --x-scale-grain 1.0 --validate --seed 20261110 --save-seconds 900)
CAP=(--delta-cap 4.605170185988092); CREF=(--delta-cref 4.2688066575e-03)
run() { echo "$(date -u +%FT%TZ) $1"; shift; printf '%q ' "$@" >> $P/commands.sh; echo >> $P/commands.sh
        bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ "$@" > $P/$(basename "${@: -1}").log 2>&1; echo "  rc $?"; }
rc=0
if [ "$GROUP" = C ]; then
  for ch in stage2_f1_pm0.01 stage2_f10_pm0.01; do
    C=$(ls -d $L/${ch}_2*/ | head -1); FIT=$(readlink -f ${C}c5/fit/final)
    for f in 0.1 1 10; do
      O=$W/C_${ch#stage2_}_fd${f}_cap
      run "C $ch f_delta $f" $PY -W ignore $M/localnn_ft.py run --fit $FIT --late $L/late_kernels_150m.npz --design $L/grid_150m/design \
        --base $R/base --d2 $R/clip_refit/D2/final "${COMMON[@]}" --delta-penalty-f $f "${CREF[@]}" "${CAP[@]}" --out $O
    done
  done
else
  for spec in "B background/stage2_bg1e-2_20261115T184949Z B_c5" "D logspread/stage2_f0.1_20261115T184949Z D_c5"; do
    set -- $spec; G=$1; FIT=$N/$2/c5/fit/seg01; ST=$3
    run "$G control (unregularized, CPU)" $PY -W ignore $M/localnn_ft.py run-nonneg --fit $FIT --late $W/late_kernels_BD.npz \
      --state-name $ST "${COMMON[@]}" --out $W/${G}_control_cpu
    for f in 0.1 1 10; do
      run "$G f_delta $f" $PY -W ignore $M/localnn_ft.py run-nonneg --fit $FIT --late $W/late_kernels_BD.npz --state-name $ST \
        "${COMMON[@]}" --delta-penalty-f $f "${CREF[@]}" "${CAP[@]}" --out $W/${G}_fd${f}_cap
    done
  done
fi
echo "$(date -u +%FT%TZ) done"
