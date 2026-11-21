#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/tests/localnn_q002/d2check.sh (on Voyager: /data001/highresolu_vt/cops_output/tests/localnn_q002/d2check.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 56ec33dfe5e43950957fae2be5497e7435e524d50d421760ae499e033391167b
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: the decisive preflight check of the 0.02 design against the NEW D2: the design build with --d2 =
#   D2_q002_warm/final into a test folder (matrices bitwise = design_h2oqmin0.02; the stored direct sums are the new
#   D2's), then localnn_fit verify (the new D2's data cost and tables).
# RUN 2026-11-18: 136558: D2 data cost rel 6.7e-16, COPS_model 0 printed differences, late window <= 3.5e-10.
# ==== end of header ====
# the decisive preflight check, done properly: localnn_fit verification (test_localnn 8) compares the D2 tables with
# the design's STORED direct sums of the D2 given at the design build (the 0.2-law D2 for design_h2oqmin0.02), so the
# preflight's test 8 paired old-D2 sums with the new D2's tables. Here: the design build again with --d2 = the NEW D2
# into a TEST folder (not the runs tree), its matrices compared bitwise with design_h2oqmin0.02, then the verification
# (test_localnn 8) with that design and the new D2 / base_q002.
set -uo pipefail
export NPY_DISABLE_CPU_FEATURES="AVX512F AVX512CD AVX512_KNL AVX512_KNM AVX512_SKX AVX512_CLX AVX512_CNL AVX512_ICL"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 PYTHONUNBUFFERED=1
APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda; PY=$APP/.venv/bin/python; cd $APP
M=/data001/highresolu_vt/cops_output/model_ext_frozen/localnn_s2law_20261118T081952Z
CC=/data001/highresolu_vt/cops_output/model_ext_frozen/nonneg_contrib_20261118T072729Z
R=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5; L=$R/nonneg/localnn; G=$L/grid_150m
C9=$R/nonneg/softplus/chain_5cycles_smm1_20261113T041215Z/c9; D2N=$R/clip_refit/D2_q002_warm/final; BQ=$R/base_q002
W=/data001/highresolu_vt/cops_output/tests/localnn_q002/d2check_$SLURM_JOB_ID; mkdir -p $W
echo "$(date -u +%FT%TZ) job $SLURM_JOB_ID on $(hostname)"; cp -p "$0" $W/
bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' _ $PY -W ignore $M/localnn_design.py --base $BQ --grid $G \
  --slope-bins $R/slope_mult/bundle/slope_bins.npz --facets-from $C9/fit/seg01/nonneg_facets.npz --d2 $D2N \
  --c9-ft $C9/ft/seg01 --c9-fit $C9/fit/seg01 --contrib-code $CC --out $W/design_d2q002 --procs 64 > $W/design.log 2>&1
echo "design rc $?"; tail -2 $W/design.log | cut -c1-300
$PY - <<PYEOF
import numpy as np
a=np.load("$W/design_d2q002/design.npz"); b=np.load("$G/design_h2oqmin0.02/design.npz")
for k in ("G_NG","G_RG","T_NG","T_RG","C_NG","C_RG","iPoint"):
    print(k, "bitwise equal to design_h2oqmin0.02:", bool(np.array_equal(a[k], b[k])))
print("sum_D2 new vs build-time (old D2) max rel", float(np.max(np.abs(a["sum_D2"]/np.where(b["sum_D2"]!=0,b["sum_D2"],1)-1))))
PYEOF
export OMP_NUM_THREADS=16 OPENBLAS_NUM_THREADS=16
$PY -W ignore $M/localnn_fit.py verify --design $W/design_d2q002 --base $BQ --d2 $D2N --point-table $R/late_grid_points.txt \
  --c9-ft $C9/ft/final --c9-fit $C9/fit/final --out $W/verify > $W/verify.log 2>&1
echo "verify rc $?"; grep -i "D2\b\|reproduction\|data cost\|COPS_model\|late_window" $W/verify.log | tail -8 | cut -c1-400
echo "$(date -u +%FT%TZ) done"
