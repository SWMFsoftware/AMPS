#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/cops_movies_prep.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/cops_movies_prep.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 224a55686c3e6631b5aaef41f70f02317a97dc0626c0f968da0a8db654b3aa26
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (GPU, 1 GPU): movie preparation with provenance: build the SPICE helper; cops_movies.py prep, fov,
#   sc, plan, validate, validate_fov, validate_sc (STEPS); 2026-11-13 also copy (FROM_WORK; with contact.json when
#   present), plan with CAMERA=fixed or CAMERA=detach (CAMERA_ALT / DETACH_ALT, FOV_DEG), validate_fixed (MASK_RG,
#   MASK_NG, SAMPLES, TEST_FRAMES), touchdown (TRACK_BASE, IMPACT_JSON), contact (HGA_CK, HGA_HOLD_UTC, SA_HOLD_UTC),
#   prep with EXTEND_TO_CONTACT / EXTEND_TO_UTC, sc with HGA_CK; MOVIES_PY selects the cops_movies.py to run.
# RUN: STEPS=... WORK=... GEOM=... sbatch --export=ALL cops_movies_prep.sbatch (OUTPUT_VISUALIZATION.md 3.7). Jobs:
#   124615-124634, 124708, 124724, 124725, 124786, 124787, 124863, 124864, 124894 (detach-200 m), 124909 (detach-300 m)
# ==== end of header ====
# COPS movies, step 1 (user request 2026-11-12, one-minute movies of the descent as seen from RG / NG):
# build the SPICE geometry helper (figures_cluster/scripts/cops_movies_geom.cpp, CSPICE), then
# figures_cluster/scripts/cops_movies.py prep (frame times uniform in log altitude, SPICE geometry, shield table, checks
# against the base), plan (camera framing, sub-spacecraft points; GPU), validate (shadow map vs the maps geometry at its
# 13 times; cos / shadow / acceptance weights / visibility vs the kernel records of base rows; test frames; GPU).
# Inputs from the environment (defaults = the softplus 5-cycle final state; values without commas in --export=):
#   RUN, CHAIN, GEOM (geometry folder of the maps job), WORK (default $CHAIN/movies_work), STEPS (prep,plan,validate;
#   2026-11-13 also fov = IK FOVs of both gauges at every frame / row time -> $WORK/fov.npz, and validate_fov = the
#   ik_fov masks against an analytic instrument-frame test, weight ranges, kernel-record shares, masked test frames;
#   sc = spacecraft geometry for the shadow (bus + solar-array wings, SA orientation held after SA_HOLD_UTC) -> sc.npz,
#   validate_sc = shadow checks; movies only, user decisions 2026-11-13)
# 2026-11-13 (fixed-camera movies and the touchdown flag; new steps, the defaults unchanged):
#   copy = copy prep.npz, fov.npz, sc.npz and their checks / inputs from FROM_WORK (a finished gauge-camera work folder;
#   sha256 listed and compared) instead of recomputing them; plan with CAMERA=fixed (CAMERA_ALT, default 100 m; FOV_DEG,
#   default 90) = cops_movies.py plan --camera fixed; validate_fixed (MASK_RG, MASK_NG default ik_fov; SAMPLES default
#   2000); touchdown (TRACK_BASE, default runs/shap8_012m_flat_earlier/base = the ground-track figures' base;
#   IMPACT_JSON: their _impact.json) = the touchdown point, its visibility per frame and the SPICE check;
#   MOVIES_PY = the cops_movies.py to run (default figures_cluster/scripts/cops_movies.py).
# 2026-11-13 (detach-200 m set; new steps, defaults unchanged): contact = first contact of the spacecraft model with the
#   surface (HGA_CK, HGA_HOLD_UTC, SA_HOLD_UTC) -> $WORK/contact.json, run before prep; prep with EXTEND_TO_CONTACT=1
#   (frames after the last observation to that contact; or EXTEND_TO_UTC); sc with HGA_CK (HGA and solar-array frames in
#   sc.npz); plan with CAMERA=detach (DETACH_ALT, default 200 m; FOV_DEG); validate_fixed with TEST_FRAMES.
# Provenance in $WORK/provenance-<job>/. Shared-node rule: heavy processes with oom_score_adj 1000.
#SBATCH --job-name=cops_movies_prep
#SBATCH --partition=gpu
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=24G
#SBATCH --time=01:30:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_movies_prep-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export CUDA_PATH=/usr/local/cuda-12.9 CUDA_HOME=/usr/local/cuda-12.9 CUPY_CACHE_DIR=$COPS_OUTPUT_ROOT/.cupy
export PATH=$CUDA_HOME/bin:$PATH LD_LIBRARY_PATH=$CUDA_HOME/lib64${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
R=${RUN:-$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5}; RUN=$R; BASE=$R/base
CH=${CHAIN:-$R/nonneg/softplus/chain_5cycles_20261112T012553Z}
GEOM=${GEOM:-$R/maps_geometry/job_124324/geometry}
WORK=${WORK:-$CH/movies_work}; STEPS=${STEPS:-prep,plan,validate}
CSPICE=/home/vtenishev/SPICE/cspice; S=$COPS_OUTPUT_ROOT/figures_cluster/scripts
MPY=${MOVIES_PY:-$S/cops_movies.py}; FROM_WORK=${FROM_WORK:-}; CAMERA=${CAMERA:-}; CAMERA_ALT=${CAMERA_ALT:-100}; FOV_DEG=${FOV_DEG:-90}
MASK_RG=${MASK_RG:-ik_fov}; MASK_NG=${MASK_NG:-ik_fov}; SAMPLES=${SAMPLES:-2000}
TRACK_BASE=${TRACK_BASE:-$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier/base}; IMPACT_JSON=${IMPACT_JSON:-}
HGA_CK=${HGA_CK:-}; HGA_HOLD_UTC=${HGA_HOLD_UTC:-}; EXTEND_TO_CONTACT=${EXTEND_TO_CONTACT:-}; EXTEND_TO_UTC=${EXTEND_TO_UTC:-}
DETACH_ALT=${DETACH_ALT:-200}; TEST_FRAMES=${TEST_FRAMES:-}; SA_HOLD_UTC=${SA_HOLD_UTC:-}
PY="$APP/.venv/bin/python -W ignore"
stamp() { echo "$(date -u +%FT%TZ) $*"; }
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
mkdir -p "$WORK" || exit 1
P=$WORK/provenance-$SLURM_JOB_ID; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): STEPS=$STEPS WORK=$WORK"
for v in RUN BASE CH GEOM WORK STEPS MPY FROM_WORK CAMERA CAMERA_ALT FOV_DEG MASK_RG MASK_NG SAMPLES TRACK_BASE IMPACT_JSON HGA_CK HGA_HOLD_UTC EXTEND_TO_CONTACT EXTEND_TO_UTC DETACH_ALT TEST_FRAMES SA_HOLD_UTC; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/cops_movies_prep.sbatch.used"
mkdir "$P/code" && cp -p "$MPY" "$S/cops_movies_geom.cpp" "$APP/cops_cuda/figures.py" "$APP/cops_cuda/kernels/render.cu" "$P/code/" || exit 1
(cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded; CUDA from $CUDA_HOME)" >> "$P/module_list.txt"
nvidia-smi -q > "$P/nvidia_smi_q.txt" 2>&1
{ nvcc --version; gcc --version | head -1; g++ --version | head -1; grep -m1 CSPICE $CSPICE/doc/version.txt; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; "$@"; }
free -g > "$P/free_start.txt" 2>&1
nfail=0
helper() {      # the geometry helper built in this job (once), for steps after prep (2026-11-13)
    if [ ! -x "$P/cops_movies_geom_job" ]; then
        grep "^Kernel = " "$BASE/boundary.in" | sed 's/^Kernel = //; s/ *#.*//' > "$P/kernels_job.txt"
        cmp -s "$P/kernels_job.txt" "$WORK/kernels.txt" || { stamp "kernel list differs from $WORK/kernels.txt"; return 1; }
        run g++ -O2 -std=c++17 -I"$CSPICE/include" "$S/cops_movies_geom.cpp" "$CSPICE/lib/cspice.a" -lm -o "$P/cops_movies_geom_job" || return 1
        sha256sum "$P/cops_movies_geom_job" > "$P/geom_exe_job.sha256"
    fi
}
if [[ ",$STEPS," == *,copy,* ]]; then
    stamp "copy the gauge-camera products from $FROM_WORK"
    [ -n "$FROM_WORK" ] && [ "$FROM_WORK" != "$WORK" ] || { stamp "copy needs FROM_WORK (another work folder)"; exit 1; }
    for f in prep.npz prep_checks.json fov.npz fov_checks.json sc.npz sc_checks.json kernels.txt frame_utc.txt rows_utc.txt \
             maps_utc.txt frame_geom.txt rows_geom.txt maps_geom.txt shield_table.txt; do
        [ -e "$WORK/$f" ] && { stamp "$WORK/$f exists"; exit 1; }
        cp -p "$FROM_WORK/$f" "$WORK/$f" || exit 1
        cmp -s "$FROM_WORK/$f" "$WORK/$f" || { stamp "copy of $f differs"; exit 1; }
    done
    (cd "$FROM_WORK" && sha256sum prep.npz prep_checks.json fov.npz fov_checks.json sc.npz sc_checks.json kernels.txt \
        frame_utc.txt rows_utc.txt maps_utc.txt frame_geom.txt rows_geom.txt maps_geom.txt shield_table.txt) > "$P/copied_from_$(basename "$FROM_WORK").sha256"
    echo "$FROM_WORK" > "$P/copied_from.txt"
    if [ -e "$FROM_WORK/contact.json" ]; then      # 2026-11-13: the first-contact record of a detach set (end of the frames)
        cp -p "$FROM_WORK/contact.json" "$WORK/contact.json" && cmp -s "$FROM_WORK/contact.json" "$WORK/contact.json" || exit 1
        (cd "$FROM_WORK" && sha256sum contact.json) >> "$P/copied_from_$(basename "$FROM_WORK").sha256"
    fi
fi
if [[ ",$STEPS," == *,contact,* ]]; then
    stamp "contact (first contact of the spacecraft model with the surface)"
    [ -n "$HGA_CK" ] || { stamp "contact needs HGA_CK"; exit 1; }
    [ -e "$WORK/kernels.txt" ] || grep "^Kernel = " "$BASE/boundary.in" | sed 's/^Kernel = //; s/ *#.*//' > "$WORK/kernels.txt"
    helper || exit 1
    run heavy $PY "$MPY" contact --work "$WORK" --geometry "$GEOM" --kernels "$WORK/kernels.txt" --geom-exe "$P/cops_movies_geom_job" \
        --hga-ck "$HGA_CK" ${HGA_HOLD_UTC:+--hga-hold-utc "$HGA_HOLD_UTC"} ${SA_HOLD_UTC:+--sa-hold-utc "$SA_HOLD_UTC"} > "$P/contact.log" 2>&1 \
        || { stamp "contact FAILED"; tail -n 20 "$P/contact.log"; exit 1; }
    tail -n 3 "$P/contact.log"
fi
if [ -n "$EXTEND_TO_CONTACT" ]; then
    EXTEND_TO_UTC=$($PY -c "import json,sys; print(json.load(open(sys.argv[1]))['contact_utc'])" "$WORK/contact.json") || exit 1
    stamp "frames extended to the first contact $EXTEND_TO_UTC"
fi
if [[ ",$STEPS," == *,prep,* ]]; then
    stamp "build the geometry helper"
    grep "^Kernel = " "$BASE/boundary.in" | sed 's/^Kernel = //; s/ *#.*//' > "$WORK/kernels.txt"
    run g++ -O2 -std=c++17 -I"$CSPICE/include" "$S/cops_movies_geom.cpp" "$CSPICE/lib/cspice.a" -lm -o "$WORK/cops_movies_geom" || exit 1
    sha256sum "$WORK/cops_movies_geom" > "$P/geom_exe.sha256"
    stamp "prep"
    run $PY "$MPY" prep --work "$WORK" --base "$BASE" --geometry "$GEOM" --kernels "$WORK/kernels.txt" \
        --geom-exe "$WORK/cops_movies_geom" --gauge-offset 0,0,1.5 ${EXTEND_TO_UTC:+--extend-to-utc "$EXTEND_TO_UTC"} > "$P/prep.log" 2>&1 || { stamp "prep FAILED"; cat "$P/prep.log"; exit 1; }
    tail -n 5 "$P/prep.log"
fi
if [[ ",$STEPS," == *,fov,* ]]; then
    stamp "fov"
    if [[ ",$STEPS," != *,prep,* ]]; then           # the helper with the fov modes (built in prep otherwise)
        grep "^Kernel = " "$BASE/boundary.in" | sed 's/^Kernel = //; s/ *#.*//' > "$P/kernels.txt"
        cmp -s "$P/kernels.txt" "$WORK/kernels.txt" || { stamp "kernel list differs from $WORK/kernels.txt"; exit 1; }
        run g++ -O2 -std=c++17 -I"$CSPICE/include" "$S/cops_movies_geom.cpp" "$CSPICE/lib/cspice.a" -lm -o "$P/cops_movies_geom" || exit 1
        GEXE=$P/cops_movies_geom
    else
        GEXE=$WORK/cops_movies_geom
    fi
    sha256sum "$GEXE" > "$P/geom_exe_fov.sha256"
    run $PY "$MPY" fov --work "$WORK" --kernels "$WORK/kernels.txt" --geom-exe "$GEXE" > "$P/fov.log" 2>&1 \
        || { stamp "fov FAILED"; cat "$P/fov.log"; exit 1; }
    tail -n 5 "$P/fov.log"
fi
if [[ ",$STEPS," == *,sc,* ]]; then
    stamp "sc (spacecraft geometry for the shadow)"
    grep "^Kernel = " "$BASE/boundary.in" | sed 's/^Kernel = //; s/ *#.*//' > "$P/kernels_sc.txt"
    cmp -s "$P/kernels_sc.txt" "$WORK/kernels.txt" || { stamp "kernel list differs from $WORK/kernels.txt"; exit 1; }
    run g++ -O2 -std=c++17 -I"$CSPICE/include" "$S/cops_movies_geom.cpp" "$CSPICE/lib/cspice.a" -lm -o "$P/cops_movies_geom_sc" || exit 1
    sha256sum "$P/cops_movies_geom_sc" > "$P/geom_exe_sc.sha256"
    run $PY "$MPY" sc --work "$WORK" --kernels "$WORK/kernels.txt" --geom-exe "$P/cops_movies_geom_sc" \
        ${SA_HOLD_UTC:+--sa-hold-utc "$SA_HOLD_UTC"} ${HGA_CK:+--hga-ck "$HGA_CK"} ${HGA_HOLD_UTC:+--hga-hold-utc "$HGA_HOLD_UTC"} > "$P/sc.log" 2>&1 || { stamp "sc FAILED"; cat "$P/sc.log"; exit 1; }
    tail -n 3 "$P/sc.log"
fi
if [[ ",$STEPS," == *,plan,* ]]; then
    stamp "plan${CAMERA:+ (camera $CAMERA)}"
    if [ "$CAMERA" = detach ]; then
        run heavy $PY "$MPY" plan --work "$WORK" --geometry "$GEOM" --device 0 --camera detach --detach-altitude "$DETACH_ALT" \
            --fov-deg "$FOV_DEG" > "$P/plan.log" 2>&1 || { stamp "plan FAILED"; nfail=$((nfail+1)); }
    elif [ "$CAMERA" = fixed ]; then
        helper || exit 1
        run heavy $PY "$MPY" plan --work "$WORK" --geometry "$GEOM" --device 0 --camera fixed --camera-altitude "$CAMERA_ALT" \
            --fov-deg "$FOV_DEG" --base "$BASE" --kernels "$WORK/kernels.txt" --geom-exe "$P/cops_movies_geom_job" \
            --gauge-offset 0,0,1.5 > "$P/plan.log" 2>&1 || { stamp "plan FAILED"; nfail=$((nfail+1)); }
    else
        run heavy $PY "$MPY" plan --work "$WORK" --geometry "$GEOM" --device 0 > "$P/plan.log" 2>&1 || { stamp "plan FAILED"; nfail=$((nfail+1)); }
    fi
    tail -n 4 "$P/plan.log"
fi
if [[ ",$STEPS," == *,validate_fixed,* ]]; then
    stamp "validate_fixed"
    run heavy $PY "$MPY" validate_fixed --work "$WORK" --geometry "$GEOM" --device 0 --mask-rg "$MASK_RG" --mask-ng "$MASK_NG" \
        --samples "$SAMPLES" ${TEST_FRAMES:+--test-frames "$TEST_FRAMES"} > "$P/validate_fixed.log" 2>&1 || { stamp "validate_fixed FAILED"; nfail=$((nfail+1)); }
    tail -n 30 "$P/validate_fixed.log"
fi
if [[ ",$STEPS," == *,touchdown,* ]]; then
    stamp "touchdown"
    helper || exit 1
    run heavy $PY "$MPY" touchdown --work "$WORK" --geometry "$GEOM" --device 0 --track-base "$TRACK_BASE" \
        ${IMPACT_JSON:+--impact-json "$IMPACT_JSON"} --kernels "$WORK/kernels.txt" --geom-exe "$P/cops_movies_geom_job" \
        --gauge-offset 0,0,1.5 > "$P/touchdown.log" 2>&1 || { stamp "touchdown FAILED"; nfail=$((nfail+1)); }
    tail -n 12 "$P/touchdown.log"
fi
if [[ ",$STEPS," == *,validate,* ]]; then
    stamp "validate"
    run heavy $PY "$MPY" validate --work "$WORK" --geometry "$GEOM" --base "$BASE" --device 0 > "$P/validate.log" 2>&1 || { stamp "validate FAILED"; nfail=$((nfail+1)); }
    tail -n 40 "$P/validate.log"
fi
if [[ ",$STEPS," == *,validate_fov,* ]]; then
    stamp "validate_fov"
    run heavy $PY "$MPY" validate_fov --work "$WORK" --geometry "$GEOM" --base "$BASE" --device 0 > "$P/validate_fov.log" 2>&1 || { stamp "validate_fov FAILED"; nfail=$((nfail+1)); }
    tail -n 40 "$P/validate_fov.log"
fi
if [[ ",$STEPS," == *,validate_sc,* ]]; then
    stamp "validate_sc"
    run heavy $PY "$MPY" validate_sc --work "$WORK" --geometry "$GEOM" --device 0 > "$P/validate_sc.log" 2>&1 || { stamp "validate_sc FAILED"; nfail=$((nfail+1)); }
    tail -n 20 "$P/validate_sc.log"
fi
free -g > "$P/free_end.txt" 2>&1
stamp "done: $nfail failures"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_movies_prep-$SLURM_JOB_ID.out" "$P/job.log"
[ $nfail -eq 0 ]
