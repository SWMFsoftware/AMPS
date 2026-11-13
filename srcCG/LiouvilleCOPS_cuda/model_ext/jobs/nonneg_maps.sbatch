#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/nonneg_maps.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/nonneg_maps.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: a6cbbf085a8defb212ed7e42a0a178cf22045d18a4acebee0c2cc150a0f24fa5
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (CPU, normal128cores). Activity and source-rate maps of one final state (FT and FIT dirs): tests,
#   geometry build + validation against the kernel records, maps (.vtu + Tecplot), verify and ParaView read-back;
#   products <chain>/<MAPS_NAME> (default maps; 2026-11-13 option MAPS_NAME, e.g. maps_c9, so a later state of the same
#   chain does not touch the earlier maps), per-job geometry and build in runs/<run>/maps_geometry/job_<id>/;
#   EXPECT_CYCLES guard; MODEL_EXT = model_ext_frozen/nonneg_maps_20261112T021623Z.
# RUN: test 124322 (3-cycle states); production 124324 (softplus), 124325 (squared), 124667 (softplus, multipliers >=
#   1, c5); 124775 (the same chain, c9: MAPS_NAME=maps_c9 EXPECT_CYCLES=9)
# ==== end of header ====
# Activity and instantaneous source-rate maps of the non-negative H2O activity models (softplus, squared harmonics x
# slope multipliers) for visualization: model_ext/nonneg_maps.py. CPU node (normal128cores), no GPU.
# User request 2026-11-11: "save the activity and source rate maps for visualization" (final state only; H2O and CO2
# flux per facet at 00:00 .. 10:00, 10:33:10, 10:39:26 UTC; VTK .vtu + Tecplot ASCII).
#
# Steps (each product is new; nothing is overwritten):
#   1. tests: model_ext/test_nonneg_maps.py
#   2. build nonneg_maps_geom: model_ext/nonneg_maps_geom.cpp + nonneg_maps_shadow.cpp + the UNCHANGED reference
#      sources ../LiouvilleCOPS/src/*.cpp except main.cpp, the reference flags; BVH.cpp and the shadow loop with
#      -ffp-contract=off (the base builder ran them on the GPU with --fmad=false). Objects in $J/build (not in the tree)
#   3. times file, all-facet geometry (cos SZA, cast shadow) at the 13 map times + 17 validation rows of the base,
#      validate-geometry (mesh vs elements.bin and the BDF; Sun, cos, shadow vs the kernel records; Q vs COPS_model.dat)
#   4. maps of each chain in CHAINS: its last cycle c<N> (ft/final + fit/final; the chain's last event must be c<N>ft
#      with rc 0) -> <chain>/maps (written as maps.incomplete_<UTC>, renamed on success)
#   5. ParaView read-back of each maps dir (/usr/bin/pvpython model_ext/nonneg_maps_pvcheck.py)
# Environment: MODEL_EXT (code dir; default $COPS_OUTPUT_ROOT/model_ext; production: a frozen snapshot), CHAINS
#   (space-separated chain roots), EXPECT_CYCLES (optional: required last cycle number), PROCS (default
#   SLURM_CPUS_PER_TASK), MAPS_NAME (output folder name inside each chain, default maps; 2026-11-13: maps_c9 for the
#   extended chain, so that the c5 maps stay). Job products: $R/maps_geometry/job_<jobid>/ (provenance, build, times,
#   geometry, logs).
#SBATCH --job-name=cops_nonneg_maps
#SBATCH --partition=normal128cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=192G
#SBATCH --time=04:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_nonneg_maps-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
UNSET=$(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*|COPS_NONNEG_[A-Z_]*)' | sort | paste -sd' ' -)
[ -n "$UNSET" ] && unset $UNSET
CSPICE=/home/vtenishev/SPICE/cspice
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
BASE=$R/base
REF=$AMPS/srcCG/LiouvilleCOPS/src
export MODEL_EXT=$(readlink -f "${MODEL_EXT:-$COPS_OUTPUT_ROOT/model_ext}")
PROCS=${PROCS:-${SLURM_CPUS_PER_TASK:-8}}
CHAINS=${CHAINS:-}
EXPECT_CYCLES=${EXPECT_CYCLES:-}
MAPS_NAME=${MAPS_NAME:-maps}
case $MAPS_NAME in */*|"") echo "MAPS_NAME must be a plain folder name"; exit 2;; esac
PY=$APP/.venv/bin/python
PVPY=/usr/bin/pvpython
CXXF="-O3 -march=native -std=c++17 -fopenmp -Wall -Wextra -Wno-unused-parameter"
# maps step only: numpy without its AVX-512 kernels. The fits and fine tunings ran on gpu001/gpu002 (AMD EPYC 7513,
# AVX2); on the EPYC 9554 (AVX-512) CPU nodes numpy 1.24 dispatches other exp/log1p kernels, which change the last bit
# at 136 (softplus) / 29 (squared) of 12.4M facets (max rel 4.4e-16; diag job of 2026-11-12, maps_geometry/job_124318/
# diag_cpu_dispatch.txt), so the bitwise check "A_fit reproduced from nonneg_facets.npz" needs the AVX2 path. The maps
# themselves use the stored A_final (no transcendental function), so their values do not depend on this setting.
NPY_AVX2="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
stamp() { echo "$(date -u +%FT%TZ) $*"; }
# heavy processes on a shared node: the OOM killer takes them first (CLAUDE.md rule 7)
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
J=$R/maps_geometry/job_$SLURM_JOB_ID
P=$J/provenance-$SLURM_JOB_ID
mkdir -p "$R/maps_geometry" || exit 1
mkdir "$J" "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): MODEL_EXT=$MODEL_EXT CHAINS=[$CHAINS] EXPECT_CYCLES=[$EXPECT_CYCLES] MAPS_NAME=$MAPS_NAME PROCS=$PROCS J=$J"

# ---------------------------------------------------------------- provenance (CLAUDE.md rule 8)
cp -p "$0" "$P/nonneg_maps.sbatch.used"
for f in nonneg_maps.py nonneg_maps_geom.cpp nonneg_maps_shadow.cpp nonneg_maps_pvcheck.py test_nonneg_maps.py nonneg_maps_diag_dispatch.py; do
  cp -p "$MODEL_EXT/$f" "$P/" || exit 1
done
[ -f "$MODEL_EXT/SHA256SUMS" ] && cp -p "$MODEL_EXT/SHA256SUMS" "$P/MODEL_EXT_SHA256SUMS"
( cd "$P" && sha256sum nonneg_maps.py nonneg_maps_geom.cpp nonneg_maps_shadow.cpp nonneg_maps_pvcheck.py test_nonneg_maps.py \
  nonneg_maps_diag_dispatch.py nonneg_maps.sbatch.used > code_sha256.txt )
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded)" >> "$P/module_list.txt"
{ command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -q; } > "$P/nvidia_smi_q.txt" 2>&1 || echo "CPU job on $(hostname): no GPU used (nvidia-smi not available)" >> "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; g++ --version | head -1; grep -m1 CSPICE "$CSPICE/doc/version.txt" 2>/dev/null || echo "CSPICE: $CSPICE";
  command -v nvcc >/dev/null 2>&1 && nvcc --version | tail -2 || echo "nvcc: not used"; lscpu | grep -E 'Model name|^CPU\(s\)';
  "$PVPY" --version 2>&1 | head -1; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$PY" > "$P/pip_list.txt" 2>&1 || "$PY" -m pip list > "$P/pip_list.txt" 2>&1
env | grep -E '^(SLURM_JOB_ID|SLURM_JOB_NODELIST|SLURM_CPUS_PER_TASK|MODEL_EXT|CHAINS|EXPECT_CYCLES|MAPS_NAME|PROCS|OMP_NUM_THREADS)=' > "$P/env.txt"
echo "MAPS_NAME=$MAPS_NAME" >> "$P/env.txt"
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "export MODEL_EXT=$MODEL_EXT"; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; "$@"; }
free -g > "$P/free_start.txt" 2>&1; df -h "$COPS_OUTPUT_ROOT" > "$P/df_start.txt" 2>&1

# ---------------------------------------------------------------- 1. tests
stamp "1. tests"
run "$PY" -W ignore "$MODEL_EXT/test_nonneg_maps.py" --scratch "$J/tests" > "$J/tests.log" 2>&1; rc=$?
tail -n 3 "$J/tests.log"
[ $rc -eq 0 ] || { stamp "tests FAILED (rc $rc); see $J/tests.log"; exit 1; }

# ---------------------------------------------------------------- 2. build
stamp "2. build nonneg_maps_geom in $J/build"
BD=$J/build; mkdir "$BD" || exit 1
for f in "$REF"/*.cpp; do
  b=$(basename "$f" .cpp); [ "$b" = main ] && continue
  x=""; [ "$b" = BVH ] && x="-ffp-contract=off"
  run g++ $CXXF $x -I"$CSPICE/include" -I"$REF" -c "$f" -o "$BD/ref_$b.o" >> "$BD/build.log" 2>&1 || { stamp "compile $f FAILED"; exit 1; }
done
run g++ $CXXF -I"$CSPICE/include" -I"$REF" -c "$MODEL_EXT/nonneg_maps_geom.cpp" -o "$BD/nonneg_maps_geom.o" >> "$BD/build.log" 2>&1 || exit 1
run g++ $CXXF -ffp-contract=off -I"$CSPICE/include" -I"$REF" -c "$MODEL_EXT/nonneg_maps_shadow.cpp" -o "$BD/nonneg_maps_shadow.o" >> "$BD/build.log" 2>&1 || exit 1
run g++ -o "$BD/nonneg_maps_geom" "$BD"/*.o "$CSPICE/lib/cspice.a" -lm -fopenmp >> "$BD/build.log" 2>&1 || exit 1
( cd "$BD" && sha256sum nonneg_maps_geom ref_*.o nonneg_maps_*.o > sha256.txt; sha256sum "$REF"/*.cpp "$REF"/*.h > reference_sources_sha256.txt )
stamp "built $(sha256sum "$BD/nonneg_maps_geom" | cut -c1-16)"

# ---------------------------------------------------------------- 3. geometry + validation
stamp "3. times, geometry, validate-geometry"
run "$PY" -W ignore "$MODEL_EXT/nonneg_maps.py" times --base "$BASE" --out "$J/times.txt" || exit 1
run heavy "$BD/nonneg_maps_geom" "$BASE/boundary.in" "$J/times.txt" "$J/geometry" "$PROCS" > "$J/geometry.log" 2>&1; rc=$?
grep -E "nonneg_maps_geom|Shape model|Surface:" "$J/geometry.log"
[ $rc -eq 0 ] || { stamp "nonneg_maps_geom FAILED (rc $rc)"; exit 1; }
run heavy "$PY" -W ignore "$MODEL_EXT/nonneg_maps.py" validate-geometry --geometry "$J/geometry" --base "$BASE" > "$J/validate_geometry.log" 2>&1; rc=$?
cat "$J/validate_geometry.log"
[ $rc -eq 0 ] || { stamp "validate-geometry FAILED (rc $rc): no maps"; exit 1; }

# ---------------------------------------------------------------- 4. maps + 5. ParaView read-back
nfail=0
for C in $CHAINS; do
  C=${C%/}
  MODEL=$(basename "$(dirname "$C")")
  N=$(ls -d "$C"/c[0-9]* 2>/dev/null | sed 's#.*/c##' | sort -n | tail -1)
  FT=$C/c$N/ft/final; FIT=$C/c$N/fit/final; OUT=$C/$MAPS_NAME
  stamp "4. maps $MODEL: $C (last cycle c$N)"
  if [ -z "$N" ] || [ ! -e "$FT/finetune_state.npz" ] || [ ! -e "$FIT/nonneg_facets.npz" ]; then stamp "  no final state in $C"; nfail=$((nfail+1)); continue; fi
  if [ -n "$EXPECT_CYCLES" ] && [ "$N" != "$EXPECT_CYCLES" ]; then stamp "  last cycle c$N, expected c$EXPECT_CYCLES"; nfail=$((nfail+1)); continue; fi
  last=$(tail -n 1 "$C/chain_onejob_events.jsonl")
  if ! echo "$last" | grep -q "\"stage\": \"c${N}ft\"" || ! echo "$last" | grep -q '"rc": 0'; then
    stamp "  the chain's last event is not c${N}ft with rc 0: $(echo "$last" | cut -c1-200)"; nfail=$((nfail+1)); continue; fi
  grep -q '^complete yes$' "$FT/finetune_result.txt" || { stamp "  $FT not complete"; nfail=$((nfail+1)); continue; }
  [ -e "$OUT" ] && { stamp "  $OUT exists: not overwritten"; nfail=$((nfail+1)); continue; }
  run heavy env NPY_DISABLE_CPU_FEATURES="$NPY_AVX2" "$PY" -W ignore "$MODEL_EXT/nonneg_maps.py" maps --ft-dir "$FT" --fit-dir "$FIT" --base "$BASE" \
    --geometry "$J/geometry" --out "$OUT" --expect-model "$MODEL" --procs "$PROCS" > "$J/maps_${MODEL}.log" 2>&1; rc=$?
  tail -n 25 "$J/maps_${MODEL}.log"
  if [ $rc -ne 0 ]; then stamp "  maps $MODEL FAILED (rc $rc)"; nfail=$((nfail+1)); continue; fi
  cp -rp "$P" "$OUT/provenance-$SLURM_JOB_ID"
  cp -p "$J/maps_${MODEL}.log" "$J/validate_geometry.log" "$J/geometry.log" "$J/times.txt" "$OUT/provenance-$SLURM_JOB_ID/"
  cp -p "$J/geometry/geometry_validation.json" "$J/geometry/geometry_validation.txt" "$OUT/provenance-$SLURM_JOB_ID/"
  stamp "5. ParaView read-back of $OUT"
  if [ -x "$PVPY" ]; then
    run "$PVPY" "$MODEL_EXT/nonneg_maps_pvcheck.py" "$OUT" > "$J/pvcheck_${MODEL}.txt" 2>&1; rc=$?
    grep -v "^PASS" "$J/pvcheck_${MODEL}.txt" | tail -n 5; cp -p "$J/pvcheck_${MODEL}.txt" "$OUT/provenance-$SLURM_JOB_ID/"
    [ $rc -eq 0 ] || { stamp "  pvcheck $MODEL FAILED"; nfail=$((nfail+1)); }
  else
    stamp "  $PVPY not available on $(hostname): ParaView read-back skipped"
  fi
  du -sh "$OUT"; ls -l "$OUT" | head -40
done
free -g > "$P/free_end.txt" 2>&1; df -h "$COPS_OUTPUT_ROOT" > "$P/df_end.txt" 2>&1
stamp "done: $nfail failures"
for C in $CHAINS; do C=${C%/}; [ -d "$C/$MAPS_NAME/provenance-$SLURM_JOB_ID" ] && cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_nonneg_maps-$SLURM_JOB_ID.out" "$C/$MAPS_NAME/provenance-$SLURM_JOB_ID/job.log"; done
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_nonneg_maps-$SLURM_JOB_ID.out" "$P/job.log"
[ $nfail -eq 0 ]
