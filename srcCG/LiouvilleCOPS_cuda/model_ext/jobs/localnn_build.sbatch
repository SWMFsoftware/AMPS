#!/bin/bash
# ==== version-controlled copy: header added 2026-11-14, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_build.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_build.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: cc58c869dbf950dec35a4c6b52ddf608a8142546e8bfa866a79d94ade80bda6c
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (CPU, normal128cores, option C): mesh, design, verify, start (templates off / on), tests
#   -> runs/shap8_012m_lategrid_dz1.5/nonneg/localnn/grid_<SPACING>m; MODEL_EXT = the frozen snapshot.
# RUN: 125004 (grid_150m, frozen model_ext_frozen/localnn_20261114T081912Z).
# ==== end of header ====
# Option C of the non-negative activity models (user 2026-11-14, relayed by the lead): LOCAL non-negative basis
#   A_e = m_b(e) (sum_j c_j phi_j(e) [+ A_s sum ft_ck]), c_j >= 0, phi_j = hat functions of a coarse triangulation.
# This job builds everything the fits need, on a CPU node (normal128cores), no GPU:
#   mesh    localnn_mesh.py: coarse mesh (geodesic Voronoi / Lloyd on the 12M facet graph, dual triangles), hat
#           functions at every facet, partition-of-unity check, stiffness Lg, mass M             -> $GRID_DIR
#   design  localnn_design.py: design matrices G^(b)[row, node] (NG, RG, 7 slope bins), template and CO2 columns, from
#           the base's 4.74e9 kernel records with the frozen nonneg_contrib arithmetic; direct facet sums of D2, c9 and
#           two test activities                                                                     -> $GRID_DIR/design
#   verify  localnn_fit.py verify: D2's and the c9 final state's data costs and tables reproduced  -> $GRID_DIR/verify
#   start   localnn_fit.py start --templates off|on: NNLS projection of D2, costs, lambda_C of the proposed rule
#                                                                                   -> $GRID_DIR/start_templates_{off,on}
#   tests   test_localnn.py (partition of unity, stiffness, design vs direct, Jacobian, bounds, NNLS, reproduction)
#                                                                                                -> $GRID_DIR/tests.json
# Settings (sbatch --export=ALL,VAR=value,...):
#   MODEL_EXT  the code (production: a frozen snapshot model_ext_frozen/localnn_<UTC>; default $COPS_OUTPUT_ROOT/model_ext)
#   SPACING    node spacing h [m] (default 150; user/lead 2026-11-14)
#   GRID_DIR   default $R/nonneg/localnn/grid_${SPACING}m
#   STEPS      default "mesh design verify start tests" (a step whose product exists is skipped)
#   LLOYD      Lloyd iterations (default 6); GEOM: validated nonneg_maps geometry (default maps_geometry/job_124775)
#   CONTRIB    the frozen contribution code (default model_ext_frozen/nonneg_contrib_20261112T210012Z)
#   PROCS      default SLURM_CPUS_PER_TASK
#SBATCH --job-name=cops_localnn_build
#SBATCH --partition=normal128cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=400G
#SBATCH --time=02:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_localnn_build-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
UNSET=$(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*|COPS_NONNEG_[A-Z_]*)' | sort | paste -sd' ' -)
[ -n "$UNSET" ] && unset $UNSET
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
BASE=$R/base; D2=$R/clip_refit/D2/final; PT=$R/late_grid_points.txt
C9=$R/nonneg/softplus/chain_5cycles_smm1_20261113T041215Z/c9
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
SPACING=${SPACING:-150}; LLOYD=${LLOYD:-6}
GRID_DIR=${GRID_DIR:-$R/nonneg/localnn/grid_${SPACING}m}
STEPS=${STEPS:-mesh design verify start tests}
GEOM=${GEOM:-$R/maps_geometry/job_124775/geometry}
CONTRIB=${CONTRIB:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_contrib_20261112T210012Z}
PROCS=${PROCS:-${SLURM_CPUS_PER_TASK:-8}}
PY=$APP/.venv/bin/python
stamp() { echo "$(date -u +%FT%TZ) $*"; }
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
J=$R/nonneg/localnn/build_jobs/job_$SLURM_JOB_ID; P=$J/provenance-$SLURM_JOB_ID
mkdir -p "$R/nonneg/localnn/build_jobs" || exit 1
mkdir "$J" "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): GRID_DIR=$GRID_DIR SPACING=$SPACING STEPS='$STEPS' MODEL_EXT=$MODEL_EXT PROCS=$PROCS"

# provenance (CLAUDE.md rule 8)
cp -p "$0" "$P/localnn_build.sbatch.used"
for f in localnn_mesh.py localnn_design.py localnn_problem.py localnn_fit.py localnn_direct.py test_localnn.py positive_problem.py; do
  cp -p "$MODEL_EXT/$f" "$P/" || exit 1
done
[ -f "$MODEL_EXT/SHA256SUMS" ] && cp -p "$MODEL_EXT/SHA256SUMS" "$P/MODEL_EXT_SHA256SUMS"
( cd "$P" && sha256sum *.py localnn_build.sbatch.used > code_sha256.txt )
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded)" >> "$P/module_list.txt"
{ command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -q; } > "$P/nvidia_smi_q.txt" 2>&1 || echo "CPU job on $(hostname): no GPU used (nvidia-smi not available)" >> "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; echo "CSPICE: not used"; echo "nvcc: not used"; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$PY" > "$P/pip_list.txt" 2>&1 || "$PY" -m pip list > "$P/pip_list.txt" 2>&1
env | grep -E '^(SLURM_JOB_ID|SLURM_JOB_NODELIST|SLURM_CPUS_PER_TASK|MODEL_EXT|OMP_NUM_THREADS|OPENBLAS_NUM_THREADS)=' > "$P/env.txt"
echo "SPACING=$SPACING LLOYD=$LLOYD GRID_DIR=$GRID_DIR STEPS=$STEPS GEOM=$GEOM CONTRIB=$CONTRIB PROCS=$PROCS" >> "$P/env.txt"
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "export MODEL_EXT=$MODEL_EXT"; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; "$@"; }
free -g > "$P/free_start.txt" 2>&1; df -h "$COPS_OUTPUT_ROOT" > "$P/df_start.txt" 2>&1
COM=(--design "$GRID_DIR/design" --base "$BASE" --d2 "$D2" --point-table "$PT")
rc=0
for step in $STEPS; do
  case $step in
    mesh)
      [ -e "$GRID_DIR/grid.npz" ] && { stamp "mesh: $GRID_DIR exists, skipped"; continue; }
      stamp "mesh -> $GRID_DIR"
      run heavy "$PY" -W ignore "$MODEL_EXT/localnn_mesh.py" --base "$BASE" --geometry "$GEOM" --spacing "$SPACING" \
        --lloyd "$LLOYD" --out "$GRID_DIR" --seen-from "$C9/fit/final/nonneg_facets.npz" > "$J/mesh.log" 2>&1; rc=$?
      tail -n 14 "$J/mesh.log" | cut -c1-600;;
    design)
      [ -e "$GRID_DIR/design" ] && { stamp "design exists, skipped"; continue; }
      stamp "design -> $GRID_DIR/design"
      run heavy "$PY" -W ignore "$MODEL_EXT/localnn_design.py" --base "$BASE" --grid "$GRID_DIR" \
        --slope-bins "$R/slope_mult/bundle/slope_bins.npz" --facets-from "$C9/fit/final/nonneg_facets.npz" --d2 "$D2" \
        --c9-ft "$C9/ft/final" --c9-fit "$C9/fit/final" --contrib-code "$CONTRIB" --out "$GRID_DIR/design" \
        --procs "$PROCS" > "$J/design.log" 2>&1; rc=$?
      tail -n 5 "$J/design.log" | cut -c1-800;;
    verify)
      [ -e "$GRID_DIR/verify" ] && { stamp "verify exists, skipped"; continue; }
      stamp "verify -> $GRID_DIR/verify"
      OPENBLAS_NUM_THREADS=16 run heavy "$PY" -W ignore "$MODEL_EXT/localnn_fit.py" verify "${COM[@]}" \
        --c9-ft "$C9/ft/final" --c9-fit "$C9/fit/final" --out "$GRID_DIR/verify" > "$J/verify.log" 2>&1; rc=$?
      tail -n 4 "$J/verify.log" | cut -c1-1500;;
    start)
      for t in off on; do
        [ -e "$GRID_DIR/start_templates_$t" ] && { stamp "start $t exists, skipped"; continue; }
        stamp "start --templates $t -> $GRID_DIR/start_templates_$t"
        OPENBLAS_NUM_THREADS=16 run heavy "$PY" -W ignore "$MODEL_EXT/localnn_fit.py" start "${COM[@]}" --grid "$GRID_DIR" \
          --templates $t --out "$GRID_DIR/start_templates_$t" > "$J/start_$t.log" 2>&1; rc=$?
        tail -n 3 "$J/start_$t.log" | cut -c1-1500
        [ $rc -eq 0 ] || break
      done;;
    tests)
      [ -e "$GRID_DIR/tests.json" ] && { stamp "tests exist, skipped"; continue; }
      stamp "tests -> $GRID_DIR/tests.json"
      OPENBLAS_NUM_THREADS=16 run heavy "$PY" -W ignore "$MODEL_EXT/test_localnn.py" --grid "$GRID_DIR" "${COM[@]}" \
        --c9-ft "$C9/ft/final" --c9-fit "$C9/fit/final" --contrib-code "$CONTRIB" \
        --start-file "$GRID_DIR/start_templates_off" --out "$GRID_DIR/tests.json" > "$J/tests.log" 2>&1; rc=$?
      grep -E '^(PASS|FAIL)' "$J/tests.log" | cut -c1-300; tail -n 1 "$J/tests.log";;
    *) stamp "unknown step $step"; rc=2;;
  esac
  [ $rc -eq 0 ] || { stamp "step $step FAILED (rc $rc)"; break; }
done
free -g > "$P/free_end.txt" 2>&1; df -h "$COPS_OUTPUT_ROOT" > "$P/df_end.txt" 2>&1
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_localnn_build-$SLURM_JOB_ID.out" "$P/job.log" 2>/dev/null
if [ -d "$GRID_DIR" ]; then
  cp -rp "$P" "$GRID_DIR/provenance-$SLURM_JOB_ID"; cp -p "$J"/*.log "$GRID_DIR/provenance-$SLURM_JOB_ID/" 2>/dev/null
fi
stamp "done rc $rc"
exit $rc
