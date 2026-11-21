#!/bin/bash
# ==== version-controlled copy: header added 2026-11-20, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/finish_full_set.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/finish_full_set.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 0548789da52165f84e694dd117e33ca795f4b73d9f29ce59bd99bf81c7ec6ea7
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: last CPU steps of one full set: visibility captions, regions, model_data captions, facts, CAPTIONS.md
#   (2026-11-17); 2026-11-20 (full sets of the new-law fits): BASE and REF_DIR overrides (defaults unchanged)
# RUN: see model_ext/README.md; new-law sets: jobs/full_set_q002_finish.sbatch
# ==== end of header ====
# Last CPU steps of ONE full figure set after its jobs (full_set_submit.sh): visibility captions (copied; set-1 figures
# checked pixel-identical to the c9 set's), regions/ (regions_figures.py), model_data captions (model_data_captions.py),
# checks/facts.json + facts.txt (full_set_facts.py), CAPTIONS.md (collect_captions.py). 2026-11-17. Run inside a CPU job
# (srun on normal128cores) from anywhere. Environment: KEY, FT, FIT, MAPS, CONTRIB, PHASE, MODEL_NAME, JOB_SET4 (job of
# pub_maps_figures set1/set4/pressures), CHAIN + NCYC (chains; unset for D2), CAVEATS ('|'-separated, optional).
set -uo pipefail
: "${KEY:?}" "${FT:?}" "${FIT:?}" "${MAPS:?}" "${CONTRIB:?}" "${PHASE:?}" "${MODEL_NAME:?}" "${JOB_SET4:?}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda; PY="$APP/.venv/bin/python -W ignore"
S=/data001/highresolu_vt/cops_output/figures_cluster/scripts; F=/data001/highresolu_vt/cops_output/figures_cluster/12M
R=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5; O=$F/full_set_$KEY; C9=$F/full_set_softplus_smm1_9cycles
cd "$APP" || exit 1
rc=0
# 1. visibility captions (fit-independent set 1; as the c9 set: copies + check)
for c in 12M_groundtrack_visibility_3d_caption.txt 12M_groundtrack_visibility_latlon_caption.txt; do
  [ -e "$O/visibility/$c" ] || cp -p "$C9/visibility/$c" "$O/visibility/" || rc=1
done
$PY - "$O/visibility" "$C9/visibility" > "$O/visibility/CAPTIONS_COPIED.txt" <<'PYEOF' || rc=1
import hashlib, sys
from pathlib import Path
import numpy as np
from PIL import Image
o, c = Path(sys.argv[1]), Path(sys.argv[2])
print("12M_groundtrack_visibility_3d_caption.txt and 12M_groundtrack_visibility_latlon_caption.txt are copies (cp -p) of the")
print("captions in ../../full_set_softplus_smm1_9cycles/visibility/ (themselves copies from full_set_softplus_5cycles, see")
print("there). The set-1 figures do not depend on the fit (observation-only base and its bundle). Checks of this set's")
print("regenerated set-1 products against the c9 set's:")
for f in ("12M_groundtrack_visibility_3d.png", "12M_groundtrack_visibility_latlon.png", "12M_templates_visibility.png"):
    a, b = np.asarray(Image.open(o / f)), np.asarray(Image.open(c / f))
    print("  %s pixel-identical: %s" % (f, a.shape == b.shape and bool(np.array_equal(a, b))))
for f in ("12M_templates_visibility_caption.txt", "12M_groundtrack_visibility_impact.json"):
    h = [hashlib.sha256((d / f).read_bytes()).hexdigest() for d in (o, c)]
    print("  %s byte-identical: %s" % (f, h[0] == h[1]))
PYEOF
cat "$O/visibility/CAPTIONS_COPIED.txt"
# 2. regions
BASE=${BASE:-$R/base}; REF_DIR=${REF_DIR:-$R/clip_refit/D2/final}   # 2026-11-20: overrides (defaults unchanged)
[ -e "$O/regions" ] || $PY "$S/regions_figures.py" --data "$PHASE" --base "$BASE" \
   --bundle /data001/highresolu_vt/cops_output/runs/shap8_012m_flat_earlier/bundle/restart_templates.npz --out-dir "$O/regions" \
   --model-name "$MODEL_NAME" > "$O/checks/regions_figures.log" 2>&1 || { echo "regions FAILED"; tail -5 "$O/checks/regions_figures.log"; rc=1; }
# 3. model_data captions
ls "$O/model_data/"*_caption.txt >/dev/null 2>&1 || $PY "$S/model_data_captions.py" --model-data "$O/model_data" --tag "$KEY" \
   --ft-dir "$FT" --ref-dir "$REF_DIR" --model-text "$MODEL_NAME" --job-set4 "$JOB_SET4" --job-pressures "$JOB_SET4" \
   > "$O/checks/model_data_captions.log" 2>&1 || { echo "model_data captions FAILED"; tail -5 "$O/checks/model_data_captions.log"; rc=1; }
# 4. facts
CAV=(); IFS='|' read -ra CS <<< "${CAVEATS:-}"; for c in "${CS[@]}"; do [ -n "$c" ] && CAV+=(--caveat "$c"); done
CHO=(); [ -n "${CHAIN:-}" ] && CHO=(--chain "$CHAIN" --cycles "${NCYC:-5}")
[ -e "$O/checks/facts.json" ] || $PY "$S/full_set_facts.py" --set "$O" --ft "$FT" --fit "$FIT" --maps "$MAPS" --contrib "$CONTRIB" \
   --phase "$PHASE" --base "$BASE" "${CHO[@]}" "${CAV[@]}" > "$O/checks/facts.txt" 2> "$O/checks/facts.log" || { echo "facts FAILED"; tail -5 "$O/checks/facts.log"; rc=1; }
cat "$O/checks/facts.txt"
# 5. captions
[ -e "$O/CAPTIONS.md" ] || $PY "$S/collect_captions.py" "$O" --order visibility,model_data,activity,production_rate,source_areas,cliffs,contribution_density,regions \
   --title "Figure captions: 12M $MODEL_NAME (complete figure set)" || rc=1
echo "finish $KEY rc $rc"
exit $rc
