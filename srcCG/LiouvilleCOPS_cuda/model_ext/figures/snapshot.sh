#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/snapshot.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/snapshot.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 729e9375192f1a26f1b3d22a7bb0e06834201f6bdc3e9badcbce2c2b8199ffff
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: copies a campaign's fit results, every restart file, bundles, job and plot scripts and the code version into
#   figures_cluster/<RES>/snapshot_<UTC>/ (+ snapshot_latest, SHA256SUMS).
# RUN: snapshot.sh <RES, e.g. 12M> <runner restart dir>
# ==== end of header ====
# Snapshot of a campaign next to its plots (user request 2026-10-27): the fit results and settings the plots show and
# every restart file needed to rerun or continue a fit, copied (not linked) so they survive changes in the run dir.
#     snapshot.sh <resolution label, e.g. 12M> <runner restart dir>
# Output: $COPS_OUTPUT_ROOT/figures_cluster/<label>/snapshot_<UTC>/ (+ symlink snapshot_latest), with SHA256SUMS.
# The run directory itself is not changed; nothing is deleted. Copies take ~5 GB per resolution (fit_restart.npz
# ~400 MB per 12M fit). Run on the login node or in a CPU job; it only copies files.
set -eu
RES=${1:?resolution label}; R=${2:?runner restart dir}
COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
D=$COPS_OUTPUT_ROOT/figures_cluster/$RES
U=$(date -u +%Y%m%dT%H%M%SZ); S=$D/snapshot_$U
mkdir -p "$S"
cpf() { local dst=$1; shift; mkdir -p "$dst"; for f in "$@"; do [ -e "$f" ] && cp -p "$f" "$dst/"; done; return 0; }
FITFILES="fit_restart.npz fit_result.txt fit_manifest.json fit_report.json cops_fit.in ActivityFit.txt
          ActivityEnhancementFit.txt selection_profile.json pressure_residuals.dat fit_data_used.txt"
n=0
# free-temperature (runner, "earlier") fits: <run>/fits/<M>_n2000/seg01
for d in "$R"/fits/*/seg01; do
    [ -f "$d/fit_result.txt" ] || continue
    m=$(basename "$(dirname "$d")"); cpf "$S/fits_free_temperatures/$m" $(for f in $FITFILES; do echo "$d/$f"; done); n=$((n+1))
done
# monotonic penalty + ordered temperatures: <run>/c16/fits/<M>_n2000_mono10_ordered
for d in "$R"/c16/fits/*_mono10_ordered; do
    [ -f "$d/fit_result.txt" ] || continue
    cpf "$S/fits_monotonic_ordered_temperatures/$(basename "$d")" $(for f in $FITFILES; do echo "$d/$f"; done); n=$((n+1))
done
# fine tuning: state and prepared arrays are its restart files
for d in "$R"/finetune/*; do
    [ -f "$d/fine_tuning.json" ] || continue
    cpf "$S/finetune/$(basename "$d")" "$d"/fine_tuning.json "$d"/fine_tuning_state.npz "$d"/fine_tuning_prepared.npz \
        "$d"/cops_fine_tuning.in "$d"/ActivityFit.txt "$d"/ActivityEnhancementFit.txt "$d"/pressure_comparison.txt
done
# dust grain + crossing (model_ext, user request 2026-10-27): <run>/grain/fits/<fit>, <run>/grain/finetune/<fit>_m5_i5_cpu
for d in "$R"/grain/fits/*_grain_cross10; do
    [ -f "$d/fit_result.txt" ] || continue
    cpf "$S/fits_grain/$(basename "$d")" $(for f in $FITFILES grain_term.txt grain_term.json grain_components.dat COPS_model.dat; do echo "$d/$f"; done); n=$((n+1))
done
for d in "$R"/grain/finetune/*_m5_i5_cpu; do
    [ -f "$d/fine_tuning.json" ] || continue
    cpf "$S/finetune_grain/$(basename "$d")" "$d"/fine_tuning.json "$d"/fine_tuning_state.npz "$d"/fine_tuning_prepared.npz \
        "$d"/cops_fine_tuning.in "$d"/ActivityFit.txt "$d"/ActivityEnhancementFit.txt "$d"/pressure_comparison.txt "$d"/grain_finetune.json
done
[ -d "$R/grain" ] && cpf "$S/model_ext" "$COPS_OUTPUT_ROOT"/model_ext/*.py "$COPS_OUTPUT_ROOT"/model_ext/*.sbatch
for p in "$R"/grain/provenance-*; do [ -d "$p" ] && cpf "$S/model_ext/$(basename "$p")" "$p"/*; done
# activity templates
cpf "$S/templates" "$R"/bundle/restart_templates.npz "$R"/bundle/restart_templates.json "$R"/c16/bundle/slope_templates.npz
# how the fits were run and plotted
cpf "$S/jobs" "$R"/pipeline.json "$R"/slurm/*.sbatch "$COPS_OUTPUT_ROOT"/jobs/c16.sbatch "$COPS_OUTPUT_ROOT"/jobs/submit_campaign.sh
cpf "$S/plot_scripts" "$COPS_OUTPUT_ROOT"/figures_cluster/scripts/*.sh "$COPS_OUTPUT_ROOT"/figures_cluster/scripts/*.py
{ echo "snapshot $U of $R"; git -C "$APP" rev-parse HEAD; git -C "$APP" status --short; } > "$S/code_version.txt" 2>&1
[ -f "$D/STATUS.md" ] && cp -p "$D/STATUS.md" "$S/STATUS_at_snapshot.md"
cat > "$S/README.txt" <<EOF
SHAP8 $RES snapshot $U of the run directory
  $R
Copies (not links) of what the plots in ../model_vs_observations show, and every restart file:
  fits_free_temperatures/<fit>/fit_restart.npz                  restart of the free-temperature fits
  fits_monotonic_ordered_temperatures/<fit>/fit_restart.npz     restart of the monotonic + ordered-T fits
  finetune/<fit>/fine_tuning_state.npz, fine_tuning_prepared.npz  fine-tuning restart
  fits_grain/<fit>/fit_restart.npz, grain_term.*            restart of the fits with the dust grain + crossing
  finetune_grain/<fit>/...                                     fine tuning with the grain (grain_finetune.json)
  model_ext/  the grain code as run (+ per-job provenance); plots in ../model_vs_observations_grain
  templates/  activity template bundles;  jobs/  job scripts and pipeline.json;  plot_scripts/
  code_version.txt  git commit of the application;  STATUS_at_snapshot.md  ../STATUS.md at this time
The base, record store and gravity cache are not copied (hundreds of GB); a rerun needs them from the run directory.
Model codes in the directory names: see ../model_vs_observations/NOMENCLATURE.txt.
SHA256SUMS lists every file.
EOF
(cd "$S" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
ln -sfn "snapshot_$U" "$D/snapshot_latest"
echo "snapshot $S: $n fits, $(du -sh "$S" | cut -f1)"
