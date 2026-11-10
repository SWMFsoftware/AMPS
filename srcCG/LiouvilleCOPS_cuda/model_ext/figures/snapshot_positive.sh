#!/bin/bash
# ==== version-controlled copy: header added 2026-11-09, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/snapshot_positive.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/snapshot_positive.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 0c3153e8b2544e254008aee4e6dbef6792d3135a08b5b65a76b9ad2ba69162ef
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: snapshot of the 12M positive-activity work and the clipped refits next to the plots (restart files hard-linked,
#   everything else copied; SHA256SUMS). User 2026-11-09. First run: CPU job 123566 -> snapshot_positive_20261109T233216Z.
# RUN: snapshot_positive.sh   (or as a small CPU job: it hashes ~25 GB)
# ==== end of header ====
# Snapshot of the 12M work since the last campaign snapshot (2026-10-28) next to the plots (user 2026-11-09: "update the
# 12M STATUS.md summary and make the snapshot"; answer "Hard links"): the first positive-activity lambda scan, the
# ratio-0.002 and ratio-1e-5 fit / fine-tuning cycles (runs/shap8_012m_lategrid_dz1.5) and the clipped-model refits
# D0 (runs/shap8_012m_flat_earlier/clip_refit), D1, D2.
# Large restart files (fit_restart.npz, finetune_state.npz) are HARD LINKS (same filesystem; they are never modified in
# place, so the link keeps the data even if the run directory is moved or deleted; LINKS.txt lists them); every other
# file is a copy. The bases, record stores and per-facet text products are not included (hundreds of GB; a rerun needs
# the base from the run directory). Output: figures_cluster/12M/snapshot_positive_<UTC>/ (+ link
# snapshot_positive_latest); the campaign snapshots (snapshot_latest) are not touched. Nothing is deleted.
#     sbatch snapshot_positive.sbatch   (or directly: snapshot_positive.sh; it only links/copies and hashes ~25 GB)
set -eu
COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
RUN=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; E=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier
D=$COPS_OUTPUT_ROOT/figures_cluster/12M; U=$(date -u +%Y%m%dT%H%M%SZ); S=$D/snapshot_positive_$U
[ -e "$S" ] && { echo "exists: $S" >&2; exit 2; }
mkdir -p "$S"; : > "$S/LINKS.txt"
cpf() { local dst=$1; shift; mkdir -p "$dst"; for f in "$@"; do [ -e "$f" ] && cp -pL "$f" "$dst/"; done; return 0; }
lnk() { local dst=$1 f=$2; mkdir -p "$dst"; [ -e "$f" ] || return 0
        local real; real=$(readlink -f "$f")
        ln "$real" "$dst/$(basename "$f")" || { echo "hard link failed: $real" >&2; exit 3; }
        echo "$dst/$(basename "$f") -> $real" >> "$S/LINKS.txt"; }
FITSMALL="fit_result.txt fit_manifest.json fit_report.json ACTIVITY_PARAMETERIZATION.txt CLIP_REFIT.txt grain_term.txt
          grain_term.json grain_components.dat late_window_model.dat COPS_model.dat pressure_residuals.dat
          selection_profile.json fit_data_used.txt cops_fit.in ActivityFit.txt"
FTSMALL="finetune_result.txt fine_tuning.json validation.json COPS_model.dat late_window_model.dat pressure_comparison.txt"
PROV="commands.sh git.txt segment.txt environment.txt model_ext.sha256 run_copy.sha256"
nf=0; nt=0
fit() {  # $1 = destination, $2 = fit directory (seg or final), $3 = stage directory holding provenance-*
    local dst=$1 d; d=$(readlink -f "$2")
    [ -f "$d/fit_result.txt" ] || return 0
    cpf "$dst" $(for f in $FITSMALL; do echo "$d/$f"; done); lnk "$dst" "$d/fit_restart.npz"
    for p in "$3"/provenance-*; do [ -d "$p" ] && cpf "$dst/$(basename "$p")" $(for f in $PROV; do echo "$p/$f"; done); done
    echo "$(basename "$d")" > "$dst/SEGMENT.txt"; nf=$((nf + 1)); }
ft() {
    local dst=$1 d=$2
    [ -f "$d/finetune_result.txt" ] || return 0
    cpf "$dst" $(for f in $FTSMALL; do echo "$d/$f"; done); lnk "$dst" "$d/finetune_state.npz"
    for p in "$d"/provenance-*; do [ -d "$p" ] && cpf "$dst/$(basename "$p")" $(for f in $PROV; do echo "$p/$f"; done); done
    nt=$((nt + 1)); }
# first positive-activity scan (ratio 0.02)
for d in "$RUN"/positive/C30_n4000_exp_cliff_w10x5_grain_cross10_lam*; do
    L=${d##*_lam}; fit "$S/positive_scan_r0.02/lam$L" "$d/seg01" "$d"
done
# cycles, ratio 0.002 and 1e-5
for r in 0.002 1e-5; do
    for ld in "$RUN"/positive_r$r/lam*; do
        L=$(basename "$ld")
        for c in 1 2 3; do
            fit "$S/positive_cycles_r$r/$L/c$c/fit" "$ld/c$c/fit/final" "$ld/c$c/fit"
            ft "$S/positive_cycles_r$r/$L/c$c/ft" "$ld/c$c/ft"
        done
    done
done
# clipped-model refits
fit "$S/clip_refit/D0" "$E/clip_refit/D0/final" "$E/clip_refit/D0"
for v in D1 D2; do fit "$S/clip_refit/$v" "$RUN/clip_refit/$v/final" "$RUN/clip_refit/$v"; done
# base definition (not the base), frozen code, job and plot scripts
cpf "$S/base_definition" "$RUN/late_grid.h" "$RUN/late_grid_points.txt" "$RUN/base/boundary.in" "$RUN/base/run_manifest.json" \
    "$RUN/bundle/LINKED_FROM.txt"
for z in "$COPS_OUTPUT_ROOT"/model_ext_frozen/{c1_20261108T000327Z,ft_20261108T032615Z,ft_20261108T034023Z,clip_20261108T234013Z,r1e-5_20261109T010228Z,r1e-5k_20261109T023540Z}; do
    mkdir -p "$S/model_ext_frozen/$(basename "$z")"; cp -rpL "$z/." "$S/model_ext_frozen/$(basename "$z")/"
done
chmod -R u+w "$S/model_ext_frozen"; find "$S/model_ext_frozen" -name __pycache__ -prune -exec rm -rf {} + 2>/dev/null || true
J=$COPS_OUTPUT_ROOT/jobs
cpf "$S/jobs" $J/positive_fit.sbatch $J/positive_finetune.sbatch $J/positive_cycles.sh $J/positive_held_test.sbatch \
    $J/positive_activity_check.sbatch $J/prepare_12m_lategrid.sbatch $J/base_consistency_check.sbatch $J/clip_refit.sbatch \
    $J/clip_refit_bitwise.sbatch
cpf "$S/plot_scripts" "$COPS_OUTPUT_ROOT"/figures_cluster/scripts/*.sh "$COPS_OUTPUT_ROOT"/figures_cluster/scripts/*.py
{ echo "snapshot $U"; git -C "$APP" rev-parse HEAD; git -C "$APP" status --short; } > "$S/code_version.txt" 2>&1
cp -p "$D/STATUS.md" "$S/STATUS_at_snapshot.md"
cat > "$S/README.txt" <<EOF
SHAP8 12M snapshot $U: the work since the campaign snapshot of 2026-10-28 (positive activity, cycles, clipped refits)
  runs: $RUN, $E/clip_refit
  positive_scan_r0.02/lam<L>/                      first lambda scan (fit_restart.npz + results + provenance)
  positive_cycles_r0.002/lam<L>/c<k>/{fit,ft}/     ratio 0.002 cycles (fit_restart.npz; ft: finetune_state.npz = held_next, grain)
  positive_cycles_r1e-5/lam<L>/c<k>/{fit,ft}/      ratio 1e-5 + CO2 floor cycles
  clip_refit/{D0,D1,D2}/                           clipped original-T c6 model refits
  base_definition/  late grid and boundary.in of the lategrid base (the base itself is not included)
  model_ext_frozen/ the frozen code every job ran;  jobs/, plot_scripts/;  code_version.txt;  STATUS_at_snapshot.md
fit_restart.npz and finetune_state.npz are HARD LINKS to the run directories' files (user 2026-11-09), listed in
LINKS.txt; everything else is a copy. SHA256SUMS lists every file (hard-linked ones included).
Plots: ../model_vs_observations_positive_lambda_scan, _positive_cycles, _positive_cycles_r1e-5, _clip_vs_positive.
EOF
(cd "$S" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
ln -sfn "snapshot_positive_$U" "$D/snapshot_positive_latest"
echo "snapshot $S: $nf fits, $nt fine tunings, $(wc -l < "$S/LINKS.txt") hard links, $(du -sh --apparent-size "$S" | cut -f1) apparent / $(du -sh "$S" | cut -f1) on disk (hard links counted once per tree)"
