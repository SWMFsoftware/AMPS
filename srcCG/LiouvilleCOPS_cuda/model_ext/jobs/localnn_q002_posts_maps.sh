#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/localnn_q002_posts_maps.sh (on Voyager: /data001/highresolu_vt/cops_output/jobs/localnn_q002_posts_maps.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 6364bb3f6590b4da31a2d8a5fc11d1db0b2e1aa0e05b685042e235881420818b
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: login-node helper: the five q002 option C posts WITH maps (localnn_chain_post.sbatch; afterany on each
#   chain's CPU and GPU attempts), with GEOM = maps_geometry/job_136691/geometry (validated against base_q002) and
#   CONTRIB = nonneg_contrib_20261118T192126Z (the --shape-manifest fix of nonneg_maps.py); replaces the held posts
#   (scancel) and records submissions.txt and chain_events.jsonl.
# RUN 2026-11-18: 2026-11-18 19:25Z: posts 136694-136698 replacing 136616, 136619, 136601, 136606, 136611.
# ==== end of header ====
# Agent C 2026-11-18 (lead 13:27 CST): the five q002 option C posts WITH maps, replacing the held posts 136601 136606 136611
# 136616 136619 (same dependencies; the maps step needs a geometry validated against base_q002 and the contribution code
# nonneg_contrib_20261118T192126Z). Run on the login node only AFTER job 136691 completed with its geometry validation passed.
# localnn_tools.sbatch / localnn_report.sbatch read GEOM and CONTRIB from the environment (defaults unchanged); the post
# job runs them as plain bash with its environment, so no job-script line is needed.
set -euo pipefail
ROOT=/data001/highresolu_vt/cops_output; R=$ROOT/runs/shap8_012m_lategrid_dz1.5; L=$R/nonneg/localnn; U=20261118T162407Z
F=$ROOT/model_ext_frozen/localnn_q002_20261118T162327Z
B=$(readlink -f $R/base_q002); D=$(readlink -f $R/clip_refit/D2_q002_warm/final); G=$L/grid_150m/design_h2oqmin0.02_d2q002
GEOM=$R/maps_geometry/job_136691/geometry; CONTRIB=$ROOT/model_ext_frozen/nonneg_contrib_20261118T192126Z
for x in "f1_pm1e-4_tplon on 136589:136590:136591:136592:136612:136613 136616" \
         "f10_pm1e-4_tplon on 136593:136594:136595:136596:136614:136615 136619" \
         "f1_pm0.01_tplon on 136597:136598:136599:136600 136601" \
         "f10_pm0.01_tplon on 136602:136603:136604:136605 136606" \
         "f10_pm0.01_tploff off 136607:136608:136609:136610 136611"; do
  set -- $x; C=$L/stage2r_q002_$1_$U; NM=lnn_r_q002_${1%_tpl*}_$2_post
  PO=(--parsable --partition=normal128cores --exclude="node002,node[019-024]" --job-name=$NM --dependency=afterany:$3
      --export="ALL,CHAIN=$C,CYCLES=5,MODEL_EXT=$F,BASE_DIR=$B,D2_DIR=$D,DESIGN_DIR=$G,GEOM=$GEOM,CONTRIB=$CONTRIB,TOOLS_STEPS=export contrib maps")
  J=$(sbatch "${PO[@]}" $ROOT/jobs/localnn_chain_post.sbatch)
  echo "$1: new post $J (replaces held $4)"
  echo "$(date -u +%FT%TZ) $J jobs/localnn_chain_post.sbatch normal128cores localnn q002 chain r_q002_$1 post WITH maps (replaces held $4; lead 13:27 CST: geometry job_136691 validated against base_q002, contrib code nonneg_contrib_20261118T192126Z): ${PO[*]:1}" >> $ROOT/model_ext/logs/submissions.txt
  printf '{"stage": "post_submitted", "rc": 0, "device": "cpu", "host": "%s", "job": "%s", "utc": "%s", "note": "replaces held post %s (scancelled): afterany:%s, TOOLS_STEPS export contrib maps, GEOM %s, CONTRIB %s, DESIGN_DIR %s"}\n' \
    "$(hostname)" "$J" "$(date -u +%FT%TZ)" "$4" "$3" "$GEOM" "$CONTRIB" "$G" >> $C/chain_events.jsonl
  scancel "$4" && echo "$(date -u +%FT%TZ) (scancel) $4: held post of r_q002_$1 replaced by $J" >> $ROOT/model_ext/logs/submissions.txt
done
