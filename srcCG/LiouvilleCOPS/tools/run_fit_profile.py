#!/usr/bin/env python3
"""Run one or all named LiouvilleCOPS selection profiles safely.

The runner is deliberately thin: ``fit_cops.py`` remains the only optimizer
and ``cops_selection.py`` remains the only profile interpreter.  This script
validates profile names, assigns distinct output directories, prints the exact
commands, and runs profiles sequentially.  Sequential profile execution avoids
loading the multi-gigabyte kernel set once per profile at the same time.

``--jobs`` is passed to the fitter's multiprocessing pool.  That pool runs
independent lambda values; a single-lambda fit therefore has only one optimizer
task even when ``--jobs 14`` is recorded.  BLAS/OpenMP threads default to one
per worker to prevent a 14-process scan from oversubscribing the host.
"""

import argparse
import math
import os
from pathlib import Path
import shlex
import subprocess
import sys

from cops_paths import output_path
from cops_selection import (SelectionError, available_profiles, load_profile,
                            profile_directory)


APP_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PROCESSES = 14
DEFAULT_PROFILE_ORDER = ("legacy_base_mask", "all_available", "transition_cut_100m")


def build_parser():
    names = available_profiles(profile_directory(APP_ROOT))
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", default="output_fit_full_default_base",
                        help="static LiouvilleCOPS fitting base (default: %(default)s)")
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument("--profile", action="append", choices=names,
                        help="profile to run; repeat for several profiles")
    choice.add_argument("--all-profiles", action="store_true",
                        help="run every checked-in profile in documented comparison order")
    parser.add_argument("--out", help="output directory; valid only with one --profile")
    parser.add_argument("--out-prefix",
                        help="prefix for multi-profile outputs (default: <run-name>_fit_profile)")
    parser.add_argument("--jobs", type=int, default=DEFAULT_PROCESSES,
                        help="parallel lambda-fit processes passed to fit_cops.py (default: %(default)s)")
    parser.add_argument("--threads-per-process", type=int, default=1,
                        help="OMP/BLAS threads in each fit process (default: %(default)s)")
    parser.add_argument("--lambda", dest="lambdas", type=float, nargs="+", default=[0.0],
                        help="one or more Tikhonov weights (default: literal fitter default 0)")
    parser.add_argument("--degree", type=int, default=4,
                        help="maximum fitted legacy activity degree (default: %(default)s)")
    parser.add_argument("--enhancement-degree", type=int, choices=range(4, 9), default=4,
                        help="maximum residual H2O degree; 4 is the zero-length control (default: %(default)s)")
    parser.add_argument("--fix-legacy-activity", action="store_true",
                        help="fix all 25 legacy activity coefficients (required for residual fits)")
    parser.add_argument("--enhancement-scale", type=float,
                        help="assert the required area-weighted legacy-linear RMS scale")
    parser.add_argument("--fix-temperature", action="store_true")
    parser.add_argument("--ordered-temperature", action="store_true")
    parser.add_argument("--fix-kappa", action="store_true")
    parser.add_argument("--tbounds", type=float, nargs=2, default=(30.0, 400.0),
                        metavar=("MIN_K", "MAX_K"))
    parser.add_argument("--max-nfev", type=int, default=1000)
    parser.add_argument("--ftol", type=float, default=1e-10)
    parser.add_argument("--xtol", type=float, default=1e-10)
    parser.add_argument("--gtol", type=float, default=1e-10)
    parser.add_argument("--check-jacobian", action="store_true")
    parser.add_argument("--check-tolerance", type=float, default=1e-6)
    parser.add_argument("--data", help="replacement pressure table passed to fit_cops.py")
    parser.add_argument("--hk-table", help="housekeeping table used to identify saturated NG")
    parser.add_argument("--allow-legacy-base", action="store_true")
    parser.add_argument("--data-only", action="store_true",
                        help="write audited selections and plots without optimizing")
    parser.add_argument("--dry-run", action="store_true",
                        help="validate arguments and print commands without running them")
    return parser


def output_names(args, profiles):
    """Return one non-overlapping fitter output argument per profile."""
    if args.out and len(profiles) != 1:
        raise ValueError("--out requires exactly one --profile")
    if args.out and args.out_prefix:
        raise ValueError("--out and --out-prefix cannot be combined")
    if args.out:
        return {profiles[0]: args.out}
    run_name = Path(args.run.rstrip("/")).name
    prefix = args.out_prefix or (run_name + "_fit_profile")
    return {name: prefix + "_" + name for name in profiles}


def fit_command(args, profile, destination):
    """Build the complete argv for one profile without invoking a shell."""
    command = [sys.executable, str(APP_ROOT / "tools" / "fit_cops.py"),
               args.run, "--profile", profile, "--out", destination,
               "--degree", str(args.degree),
               "--enhancement-degree", str(args.enhancement_degree), "--lambda"]
    command.extend("%.17g" % value for value in args.lambdas)
    command.extend(["--jobs", str(args.jobs), "--tbounds",
                    "%.17g" % args.tbounds[0], "%.17g" % args.tbounds[1],
                    "--max-nfev", str(args.max_nfev),
                    "--ftol", "%.17g" % args.ftol,
                    "--xtol", "%.17g" % args.xtol,
                    "--gtol", "%.17g" % args.gtol,
                    "--check-tolerance", "%.17g" % args.check_tolerance])
    for enabled, option in (
            (args.fix_temperature, "--fix-temperature"),
            (args.fix_legacy_activity, "--fix-legacy-activity"),
            (args.ordered_temperature, "--ordered-temperature"),
            (args.fix_kappa, "--fix-kappa"),
            (args.check_jacobian, "--check-jacobian"),
            (args.allow_legacy_base, "--allow-legacy-base"),
            (args.data_only, "--data-only")):
        if enabled:
            command.append(option)
    for value, option in ((args.data, "--data"), (args.hk_table, "--hk-table")):
        if value:
            command.extend([option, value])
    if args.enhancement_scale is not None:
        command.extend(["--enhancement-scale", "%.17g" % args.enhancement_scale])
    return command


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.threads_per_process < 1:
        parser.error("--threads-per-process must be positive")
    if any(value < 0 for value in args.lambdas):
        parser.error("--lambda values must be nonnegative")
    if args.ordered_temperature and args.fix_temperature:
        parser.error("--ordered-temperature has no effect with --fix-temperature")
    if not 0 <= args.degree <= 4:
        parser.error("--degree must be in 0..4")
    if args.enhancement_scale is not None and (
            not math.isfinite(args.enhancement_scale) or args.enhancement_scale <= 0):
        parser.error("--enhancement-scale must be finite and positive")
    if args.enhancement_degree > 4:
        missing = [option for enabled, option in (
            (args.fix_legacy_activity, "--fix-legacy-activity"),
            (args.fix_temperature, "--fix-temperature"),
            (args.fix_kappa, "--fix-kappa")) if not enabled]
        if missing:
            parser.error("--enhancement-degree 5..8 currently requires " + " ".join(missing))
    if args.fix_legacy_activity and not (args.fix_temperature and args.fix_kappa):
        parser.error("--fix-legacy-activity currently also requires --fix-temperature --fix-kappa")

    if args.all_profiles:
        available = available_profiles(profile_directory(APP_ROOT))
        profiles = ([name for name in DEFAULT_PROFILE_ORDER if name in available]
                    + [name for name in available if name not in DEFAULT_PROFILE_ORDER])
    else:
        profiles = args.profile
    try:
        for name in profiles:
            load_profile(name, profile_directory(APP_ROOT))
        destinations = output_names(args, profiles)
    except (SelectionError, ValueError) as exc:
        parser.error(str(exc))

    # Refuse the batch before its first expensive cache audit if any result
    # destination already exists.  fit_cops.py independently repeats this
    # protection immediately before creating each directory.
    resolved = {}
    for name, destination in destinations.items():
        try:
            path = output_path(destination)
        except ValueError as exc:
            parser.error(str(exc))
        if path.exists() and not args.dry_run:
            parser.error("output already exists for %s: %s" % (name, path))
        resolved[name] = path
    if len({str(path.resolve()) for path in resolved.values()}) != len(resolved):
        parser.error("profiles resolve to duplicate output directories")

    environment = os.environ.copy()
    for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        environment[variable] = str(args.threads_per_process)

    independent_fits = len(set(args.lambdas))
    if args.data_only:
        print("Parallelism: data-only mode runs no optimizer; --jobs=%d is recorded but creates no worker pool." %
              args.jobs, flush=True)
    else:
        effective = min(args.jobs, independent_fits)
        print("Parallelism: %d distinct lambda value(s), so at most %d of the requested %d fit process(es) can run." %
              (independent_fits, effective, args.jobs), flush=True)

    for index, name in enumerate(profiles, 1):
        command = fit_command(args, name, destinations[name])
        print("[%d/%d] %s" % (index, len(profiles), name), flush=True)
        print("env OMP_NUM_THREADS=%d OPENBLAS_NUM_THREADS=%d MKL_NUM_THREADS=%d %s" %
              (args.threads_per_process, args.threads_per_process,
               args.threads_per_process, shlex.join(command)), flush=True)
        if not args.dry_run:
            completed = subprocess.run(command, cwd=str(APP_ROOT), env=environment)
            if completed.returncode:
                return completed.returncode
    return 0


if __name__ == "__main__":
    sys.exit(main())
