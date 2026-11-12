# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/mpi_multinode/specs/clip_D2.sh (on Voyager: /data001/highresolu_vt/cops_output/mpi_multinode/specs/clip_D2.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: c1fd7dbfdc192c0fca71fd11c29e4c9069394d71ad954076934642c4b5f0694e
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: spec for jobs/mpi_fit.sbatch (test, 2026-11-10): clip_refit.py --variant D2 exactly as job 122866 (lategrid
#   base, W10x5, c6 start) except the output directory and the launcher's options (--max-nfev, --segment, --wall-time,
#   --transport, --activity); CASE, TESTSET, NFEV, CHECK.
# USED BY: sourced by jobs/mpi_fit.sbatch (MPI tests 123629-123632, 123638, 123645, 123647, 123649).
# ==== end of header ====
# Spec for jobs/mpi_fit.sbatch: clip_refit.py --variant D2 (the clipped original-T cycle-6 model with W10x5 on the
# lategrid base), exactly the command of job 122866 (runs/shap8_012m_lategrid_dz1.5/clip_refit/D2/provenance-122866/
# commands.sh) except: the output directory, --max-nfev/--segment/--wall-time/--transport/--activity (set by the
# launcher). Start: c6's final state (clip_refit.py --start-restart), as jobs/clip_refit.sbatch segment 1.
# Test of the MPI transport across nodes (2026-11-10), not a production fit.
#   CASE     test case name (required): products in runs/shap8_012m_lategrid_dz1.5/clip_refit_tests/mpi_<TESTSET>/<CASE>
#   TESTSET  test set name (default 20261110)
#   NFEV     evaluations over all segments (default 5)
#   CHECK    1 (default): in segment 1 grain_fit --check-jacobian and clip_refit --check-full-jacobian (as the job)
[ -n "${CASE:-}" ] || { echo "clip_D2.sh: set CASE" >&2; return 1; }
if env | grep -q '^COPS_POS_'; then echo "clip_D2.sh: COPS_POS_* set; this model runs the unchanged GPU physics" >&2; return 1; fi
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; E=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier
C6=$E/grain/cycles_Torig/c6/C30_n4000_Torig_grain_cross10_ftc6
TESTSET=${TESTSET:-20261110}
NAME=clip_D2_mpitest_${CASE}
D=$R/clip_refit_tests/mpi_$TESTSET/$CASE
NFEV=${NFEV:-5}
DRIVER=clip_refit.py
TA=(--activity-template "$E/bundle/restart_templates.npz" cliff_ge_60deg)
FIX=(--fix Tmax=182.09922 --fix T75=136.11149 --fix Tmin=133.0)
for j in 1 2 3 4 5 6; do
    if [ $j -le 3 ]; then TD=$E/grain/cycles/c$j; else TD=$E/grain/cycles_Torig/c$j; fi
    TA+=(--activity-template "$TD/template_ft_c$j/ft_c$j.templates.npz" ft_c$j); FIX+=(--fix c_ft_c$j=1)
done
REQUIRED=("$R/base/COPS_model.dat" "$R/base/elements.bin" "$R/base/kernels" "$R/late_grid_points.txt" "$C6/fit_restart.npz"
          $(printf '%s\n' "${TA[@]}" | grep '\.npz$'))
fit_args() {
    local start=() chk=() gchk=()
    [ "$SEG" = 1 ] && start=(--start-restart "$C6/fit_restart.npz")
    [ "$SEG" = 1 ] && [ "${CHECK:-1}" = 1 ] && { chk=(--check-full-jacobian); gchk=(--check-jacobian); }
    ARGS=(--variant D2 --point-table "$R/late_grid_points.txt" "${start[@]}" "${chk[@]}" "${FIX[@]}"
          --crossing-weight 10 "${gchk[@]}" --devices "$DEVICES" --precision mixed --template-start 0 -- "$R/base"
          --profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400
          --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6 --out "$OUT" "${TA[@]}")
}
