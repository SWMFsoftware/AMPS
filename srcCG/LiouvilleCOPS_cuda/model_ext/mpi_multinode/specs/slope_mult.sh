# ==== version-controlled copy: header added 2026-11-11, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/mpi_multinode/specs/slope_mult.sh (on Voyager: /data001/highresolu_vt/cops_output/mpi_multinode/specs/slope_mult.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: ddbd7b0d7fb85587c9c7dc2de781e6fefec4b7a97f32891b9086190f41aed9b4
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: spec for jobs/mpi_fit.sbatch: the slope-multiplier fit (slope_mult_fit.py) with the command of
#   jobs/slope_mult_fit.sbatch over the MPI transport, ACTIVITY_DEFAULT=host; TESTNAME (empty: the production directory
#   slope_mult/mpi_fit), NFEV, CHECK, CHECKPARAMS, START srccg | d2, START_RESTART, BUNDLE (2026-11-10).
# USED BY: sourced by jobs/mpi_fit.sbatch (tests 123646, 123648, bitwise = the one-node 123635); the production line was
#   not used (user 2026-11-10: run on one node now).
# ==== end of header ====
# Spec for jobs/mpi_fit.sbatch: the slope-multiplier H2O activity fit (model_ext/slope_mult_fit.py of the slope agent;
# user decisions 2026-11-10), the command of jobs/slope_mult_fit.sbatch with the MPI transport. Nothing of the model
# is set here beyond that job's command: lategrid base, D2's late-window treatment (W10x5), held temperatures
# Tmax 182.09922 / T75 136.11149 / Tmin 133 K and c_slope_00_40 = 1, crossing weight 10, lambda 0, --tbounds 30 400,
# ftol/xtol/gtol 1e-10, --precision mixed; the start is the driver's --start (START below).
# The driver installs its GPU patch in EVERY process from its own command line and needs --activity host under MPI
# (the activity is a host field): ACTIVITY_DEFAULT=host; the launcher passes --activity host.
#   TESTNAME  empty: the production fit runs/shap8_012m_lategrid_dz1.5/slope_mult/mpi_fit/seg<NN>
#             set:   test runs  runs/shap8_012m_lategrid_dz1.5/slope_mult_tests/mpi_<TESTNAME>/seg<NN>
#   NFEV      evaluations over all segments (default 4000)
#   CHECK     1 (default): segment 1 runs grain_fit --check-jacobian and --check-full-jacobian (as the one-node job)
#   CHECKPARAMS  ':'-separated parameters of the full Jacobian check (default: the driver's)
#   START     the driver's --start (required by snapshot slopemul_20261110T094013Z and later; given on EVERY segment):
#             srccg (default; user decision 2026-11-10: original srcCG harmonic coefficients, multipliers 1, kappa 1,
#             grain from D2's restart) | d2
#   START_RESTART  the --start-restart of segment 1 (default D2's final fit_restart.npz: the grain for srccg, the whole
#             start for d2); BUNDLE the slope-bin bundle
# The inherited COPS_POS_* / COPS_SLOPEMUL_* variables are unset (as jobs/slope_mult_fit.sbatch does).
_unset=$(env | grep -o -E '^(COPS_POS_[A-Z_]*|COPS_SLOPEMUL_[A-Z_]*)' | sort | paste -sd' ' -)
[ -n "$_unset" ] && { echo "slope_mult.sh: unset (inherited): $_unset"; unset $_unset; }
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5
BUNDLE=$(readlink -f "${BUNDLE:-$R/slope_mult/bundle/slope_bins.npz}")
START_RESTART=$(readlink -f "${START_RESTART:-$R/clip_refit/D2/final/fit_restart.npz}")
if [ -n "${TESTNAME:-}" ]; then NAME=slope_mult_mpitest_$TESTNAME; D=$R/slope_mult_tests/mpi_$TESTNAME
else NAME=slope_mult_mpi_fit; D=$R/slope_mult/mpi_fit; fi
NFEV=${NFEV:-4000}
START=${START:-srccg}
case $START in srccg|d2) ;; *) echo "slope_mult.sh: START=srccg|d2" >&2; return 1;; esac
DRIVER=slope_mult_fit.py
ACTIVITY_DEFAULT=host
REQUIRED=("$R/base/COPS_model.dat" "$R/base/elements.bin" "$R/base/kernels" "$R/late_grid_points.txt" "$START_RESTART" "$BUNDLE")
fit_args() {
    local start=() chk=() gchk=() chkp=() dev=()
    [ "$TRANSPORT" = local ] && dev=(--devices "$DEVICES")    # MPI: no --devices (the slope agent's READY, section 3)
    [ "$SEG" = 1 ] && start=(--start-restart "$START_RESTART")
    if [ "$SEG" = 1 ] && [ "${CHECK:-1}" = 1 ]; then
        chk=(--check-full-jacobian); gchk=(--check-jacobian)
        [ -n "${CHECKPARAMS:-}" ] && chkp=(--check-parameters "${CHECKPARAMS//:/,}")
    fi
    ARGS=(--point-table "$R/late_grid_points.txt" --slope-bundle "$BUNDLE" --start "$START" "${start[@]}" "${chk[@]}" "${chkp[@]}"
          --fix Tmax=182.09922 --fix T75=136.11149 --fix Tmin=133.0 --fix c_slope_00_40=1 --crossing-weight 10 "${gchk[@]}"
          "${dev[@]}" --precision mixed --template-start 0 -- "$R/base"
          --profile all_available --degree 4 --enhancement-degree 4 --lambda 0 --jobs 1 --tbounds 30 400
          --ftol 1e-10 --xtol 1e-10 --gtol 1e-10 --check-tolerance 1e-6 --out "$OUT")
}
