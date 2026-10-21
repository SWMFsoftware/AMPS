#!/usr/bin/env python3
"""Run the unchanged ``fit_cops`` fitter with GPU pressure/Jacobian evaluation (CUDA env).

Usage::

    .venv/bin/python tools/gpu_fit.py [--device N | --devices 0,1] [--max-nfev N]
        [--template-start C] [--cpu-comparison] [--mode tree|sequential] [--restart-from R]
        [--no-save-restart] [--no-figures] -- <fit_cops.py arguments>

Run settings decided by the user (2026-10-20):

* ``--max-nfev`` is a driver parameter: **2000 for production fits (default)** and
  **1000 for comparisons with the CPU** (``--cpu-comparison``). It is passed to the fitter, so
  ``--max-nfev`` must not also appear after ``--``.
* ``--template-start C``, **default 1.0**: the optimizer starts every nonnegative
  template amplitude at C instead of 0. The lower bound stays 0. With one clip,
  ``max(0, baseline + c·A_scale·T)``, a template on facets the free baseline drives
  negative has a zero derivative at c = 0 and cannot leave it (STATUS.md C07, N30).
  Only the start of ``scipy.optimize.least_squares`` changes. ``P.x0`` stays the
  fitter's prior (templates 0), as do the regularization, D, the "initial" column of
  ``fit_result.txt`` and ``parameters.initial_scaled`` in the manifest. The actual start
  vector is recorded in ``compute_backend.fit_start``.
* ``--cpu-comparison``: ``--max-nfev 1000`` and ``--template-start 0``, i.e. exactly the
  CPU fitter's command.

Everything after ``--`` goes verbatim to ``fit_cops.main()``: the run, the
profile, ``--out``, degrees, templates, tolerances, and so on. The fitter is the
reference module, imported from ``../LiouvilleCOPS/tools`` and unmodified on disk:
the same argument parsing, cache and consistency audits, selection,
``Problem``, ``scipy.optimize.least_squares`` (TRF), writers and manifest.

What differs, inside this process only:

* ``fit_cops.Model`` is ``cops_cuda.gpu_model.gpu_model_class()``, a subclass of
  the reference ``Model`` that reads kernel headers only and evaluates
  pressures, Jacobian columns and ``seen_faces`` on the GPU (FP64,
  approved budgets in VALIDATION.md);
* with ``--devices 0,1`` one fit runs on both GPUs (``cops_cuda.multi_gpu``: one
  optimizer, row-partitioned persistent workers; bitwise equal to one GPU);
* the fit manifest gains a ``compute_backend`` section, written on every
  manifest save, which states the CUDA backend, the device and dataset
  identities, the binding events and the call counts and times. The
  manifest's ``command`` field lists the fit_cops arguments (argparse reads
  ``sys.argv``); ``compute_backend.driver_command`` is the actual process
  command.

Refused: several ``--lambda`` values. The fitter runs an L-curve in forked
processes, and a forked child cannot use the parent's CUDA context. Run one
λ per process instead.

At the end of a fit the driver also writes, into ``<out>``:

* ``fit_restart.npz``: the restart file (ROADMAP C12; ``--no-save-restart`` skips it);
* ``figures/pressures_{time,altitude}.{eps,png}``: the publication pressure figures
  (``cops_cuda/figures.py``, user request 2026-10-20). They show the base model, the fit and
  the COPS data, with the late-visibility window shaded. ``--no-figures`` skips them, and a
  figure failure is recorded in ``compute_backend.figures`` without failing the fit.

The process-failure handling of ``fit_cops``'s ``__main__`` block is reproduced:
after the output directory exists, a failure is recorded as
``states.process = "failed"`` in the manifest.
"""

import argparse
import os
from pathlib import Path

import numpy as np
import sys
import time

APP_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(APP_ROOT))


def lambda_count(fit_args):
    """Number of values given to ``--lambda`` (nargs='+'), 1 if absent."""
    count = 1
    for k, token in enumerate(fit_args):
        if token == "--lambda":
            values = 0
            for value in fit_args[k + 1:]:
                if value.startswith("--"):
                    break
                values += 1
            count = values
        elif token.startswith("--lambda="):
            count = 1
    return count


def install_template_start(fit_cops, value, record, restart_state=None, captured=None):
    """Start the optimizer with template amplitudes at ``value``; the prior is unchanged.

    ``fit_cops.fit_worker`` calls ``fit_one(lam)`` through the module global, so
    replacing that attribute changes only the start passed to ``least_squares``.
    Parameter layout (``fit_cops.Problem``): a (na), harmonics (nh), templates (nt), T, κ.

    With ``restart_state`` (``cops_cuda.restart.Restart``), every parameter the restart defines
    (matched by name; see ``restart.start_vector``) then overrides that start; parameters
    the restart lacks keep the template start or the fitter's start. ``captured`` (a
    list) receives each result of ``fit_one`` (for the restart file written at the end).
    """
    original = fit_cops.fit_one

    def fit_one(lam, x_start=None):
        P = fit_cops.P
        if x_start is None and value != 0.0 and getattr(P, "nt", 0):
            x_start = P.x0.copy()
            first = P.na + P.nh
            names = list(P.names[first:first + P.nt])
            if names != ["c_" + n for n in P.template_names]:
                raise RuntimeError("unexpected parameter layout %s" % names)
            x_start[first:first + P.nt] = value
            if not (np.all(x_start >= P.lo) and np.all(x_start <= P.hi)):
                raise RuntimeError("template start outside the bounds")
            record.update({"applied": True, "template_parameters": names,
                           "x_start_scaled": [float(v) for v in x_start],
                           "prior_x0_unchanged": True})
        if x_start is None and restart_state is not None:
            x_start = P.x0.copy()
        if restart_state is not None:
            from cops_cuda import restart as rs
            xr, report = rs.start_vector(P, restart_state)
            for k, name in enumerate(P.names):
                if name in report["from_restart"]:
                    x_start[k] = xr[k]
            record.update({"applied": True, "restart": report,
                           "x_start_scaled": [float(v) for v in x_start], "prior_x0_unchanged": True})
        result = original(lam, x_start)
        if captured is not None:
            captured.append(result)
        return result

    fit_cops.fit_one = fit_one
    return original


def restart_template_args(state, fit_args, info):
    """Templates of a restart, as --activity-template arguments, unless the user gave some.

    Each source bundle must still be the file the restart recorded (same SHA-256).
    """
    from cops_cuda import provenance
    if any(a == "--activity-template" for a in fit_args) or not state.template_names:
        info["templates_from_restart"] = False
        return fit_args
    extra = []
    for source in state.meta["templates"]["sources"]:
        bundle = source["bundle"]
        path = bundle["path"]
        if not Path(path).is_file():
            raise SystemExit("restart template bundle is missing: %s" % path)
        if provenance.file_identity(path)["sha256"] != bundle["sha256"]:
            raise SystemExit("restart template bundle changed since the restart was written: %s" % path)
        extra += ["--activity-template", path, source["name"]]
    info["templates_from_restart"] = [s["name"] for s in state.meta["templates"]["sources"]]
    return fit_args + extra


def check_restart_mesh(state, fit_args, fit_cops):
    """The restart must belong to the base run's mesh (otherwise: restart_tool.py translate)."""
    from cops_cuda import reference
    run = next(a for a in fit_args if not a.startswith("-"))
    run = Path(reference.load("cops_paths").product_input(run))
    elements = fit_cops.read_elements(str(run / "elements.bin"))
    identity = reference.load("activity_templates").mesh_identity(elements["centroid"], elements["area"])
    if identity != state.mesh_identity():
        raise SystemExit("restart mesh %s differs from the base run's mesh %s; translate it with "
                         "tools/restart_tool.py translate" % (state.mesh_identity()[:16], identity[:16]))


def headroom_now():
    from cops_cuda import memory
    return memory.headroom_fraction_default()


def write_figures(fit_cops, fit_args, info):
    """Publication pressure figures of the finished fit (cops_cuda/figures.py).

    Base model, this fit and the COPS data with the fit's selection; the late-visibility window
    comes from the first ``--activity-template`` bundle, else the c05_write_bundle default.
    A failure here is recorded in the manifest and never fails the fit.
    """
    from cops_cuda import figures
    try:
        out = Path(fit_cops.ACTIVE_MANIFEST_PATH).parent
        base_dir = Path(fit_cops.ACTIVE_MANIFEST["inputs"]["base_elements"]["path"]).parent
        bundle = None
        for k, token in enumerate(fit_args[:-1]):
            if token == "--activity-template":
                bundle = fit_args[k + 1]
                break
        base = figures.read_table(base_dir)
        fit = figures.read_table(out)
        window = figures.late_window(base, bundle)
        info["paths"] = figures.pressure_figures(base, [(out.name, fit)], fit, window,
                                                 out / "figures" / "pressures", overwrite=False)
        info["late_window"] = {k: window[k] for k in ("start", "stop", "source", "altitude_min_m",
                                                       "altitude_max_m")}
        print("gpu_fit: figures %s" % ", ".join(info["paths"]), flush=True)
    except Exception as exc:                      # figures never fail a fit
        info["error"] = "%s: %s" % (exc.__class__.__name__, exc)
        print("gpu_fit: figures not written (%s)" % info["error"], flush=True)


def main():
    argv = sys.argv[1:]
    if "--" not in argv:
        raise SystemExit("usage: gpu_fit.py [--device N | --devices 0,1] [--mode tree|sequential] "
                         "-- <fit_cops args>")
    split = argv.index("--")
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--devices", default=None,
                        help="comma-separated devices for one fit split over several GPUs, e.g. 0,1")
    parser.add_argument("--mode", choices=("tree", "sequential"), default="tree",
                        help="tree: production; sequential: bitwise CPU order (slow, verification)")
    parser.add_argument("--records-per-block", type=int, default=4096)
    parser.add_argument("--max-nfev", type=int, default=None,
                        help="optimizer evaluation limit (default 2000; 1000 with --cpu-comparison)")
    parser.add_argument("--template-start", type=float, default=None,
                        help="optimizer start of template amplitudes (default 1.0; 0 with --cpu-comparison)")
    parser.add_argument("--cpu-comparison", action="store_true",
                        help="exactly the CPU fitter's command: --max-nfev 1000, --template-start 0")
    parser.add_argument("--restart-from", default=None,
                        help="start from a restart file (fit_restart.npz); parameters matched by name")
    parser.add_argument("--no-save-restart", action="store_true",
                        help="do not write <out>/fit_restart.npz at the end of the fit")
    parser.add_argument("--headroom-fraction", type=float, default=None,
                        help="device memory never planned, fraction of total (default 0.20; 0.05-0.5). "
                             "Needed below 0.20 only for very large bases (SHAP5 3M: 0.15)")
    parser.add_argument("--precision", choices=("fp64", "mixed"), default="fp64",
                        help="fp64 (default, the validated path) or mixed (experiment 2026-10-21): kernel "
                             "records stored in FP32, all arithmetic FP64 (meshes too large for FP64 records)")
    parser.add_argument("--reference-validation", action="store_true",
                        help="validate the base with the reference's per-record Python loop (identical "
                             "result; default: cops_cuda.base_validation, NumPy workers, ~1.8 h less on 12M)")
    parser.add_argument("--no-figures", action="store_true",
                        help="do not write the publication figures <out>/figures/pressures_{time,altitude}.{eps,png}")
    args = parser.parse_args(argv[:split])
    fit_args = argv[split + 1:]
    if any(a == "--max-nfev" or a.startswith("--max-nfev=") for a in fit_args):
        raise SystemExit("give --max-nfev to gpu_fit.py (before --), not to the fitter")
    if args.cpu_comparison and (args.max_nfev not in (None, 1000) or
                                args.template_start not in (None, 0.0)):
        raise SystemExit("--cpu-comparison fixes --max-nfev 1000 and --template-start 0")
    max_nfev = args.max_nfev if args.max_nfev is not None else (1000 if args.cpu_comparison else 2000)
    template_start = (args.template_start if args.template_start is not None
                      else (0.0 if args.cpu_comparison else 1.0))
    if max_nfev < 1 or not template_start >= 0.0:
        raise SystemExit("--max-nfev must be >= 1 and --template-start >= 0 (amplitudes are nonnegative)")
    fit_args = fit_args + ["--max-nfev", str(max_nfev)]
    if lambda_count(fit_args) != 1:
        raise SystemExit("gpu_fit.py runs exactly one --lambda value per process "
                         "(the fitter forks for an L-curve; CUDA contexts do not survive fork)")

    if args.headroom_fraction is not None:
        from cops_cuda import memory
        os.environ[memory.HEADROOM_ENVIRONMENT] = repr(float(args.headroom_fraction))
        memory.headroom_fraction_default()            # validates the range
    from cops_cuda import precision as cuda_precision
    os.environ[cuda_precision.RECORD_PRECISION_ENVIRONMENT] = "fp32" if args.precision == "mixed" else "fp64"
    from cops_cuda import gpu_model, manifest, provenance, reference
    from cops_cuda import restart as rs
    reference.bind()
    fit_cops = reference.load("fit_cops")
    from cops_cuda import base_validation
    validation_info = ({"policy": "reference"} if args.reference_validation
                       else base_validation.install(fit_cops))
    restart_state, restart_info = None, {"save": not args.no_save_restart}
    figures_info = {"write": not args.no_figures}
    if args.restart_from:
        restart_state = rs.read(args.restart_from)
        restart_info["loaded"] = provenance.file_identity(args.restart_from, "fit_restart")
        fit_args = restart_template_args(restart_state, fit_args, restart_info)
        check_restart_mesh(restart_state, fit_args, fit_cops)
    devices = [int(d) for d in args.devices.split(",")] if args.devices else [args.device]
    gpu_model.REGISTRY.configure(devices[0], args.mode, args.records_per_block, devices=devices)
    fit_cops.Model = gpu_model.gpu_model_class()

    started = provenance.utc_now()
    software = [manifest.backend_manifest("cuda_fp64_v1", d) for d in devices]
    fit_start = {"template_start": template_start, "applied": False, "max_nfev": max_nfev,
                 "cpu_comparison": bool(args.cpu_comparison),
                 "policy": "user 2026-10-20: max_nfev 2000 for production fits, 1000 for CPU "
                           "comparison; template amplitudes start at 1.0 (bound 0)"}
    original_save = fit_cops.save_active_manifest

    def save_with_backend():
        if fit_cops.ACTIVE_MANIFEST is not None:
            fit_cops.ACTIVE_MANIFEST["compute_backend"] = {
                "kind": "cuda_fp64_v1 via LiouvilleCOPS_cuda gpu_model (fit_cops.Model replaced "
                        "in-process; fit_cops.py unmodified)" +
                        (" on devices %s (multi_gpu)" % devices if len(devices) > 1 else ""),
                "driver": provenance.file_identity(Path(__file__), "tool"),
                "driver_command": [sys.executable, str(Path(__file__).resolve())] + argv,
                "fit_cops": provenance.file_identity(Path(fit_cops.__file__), "reference"),
                "started_utc": started,
                "validation_policy": ("VALIDATION.md (pressure <= 512 eps M, Jacobian <= 1024 eps M_J)"
                                      if args.precision == "fp64" else
                                      "EXPERIMENT mixed precision: FP32 records, FP64 arithmetic; "
                                      "tests/test_evaluator_mixed.py (<= 1e-6 M); not VALIDATION.md"),
                "precision": {"requested": args.precision,
                              "record_storage": "fp32" if args.precision == "mixed" else "fp64",
                              "policy": (cuda_precision.MIXED_POLICY_NAME if args.precision == "mixed"
                                         else cuda_precision.POLICY_NAME)},
                "software": software,
                "registry": gpu_model.REGISTRY.stats(),
                "fit_start": fit_start,
                "restart": restart_info,
                "figures": figures_info,
                "base_validation": validation_info,
                "memory_headroom_fraction": headroom_now(),
            }
        original_save()

    fit_cops.save_active_manifest = save_with_backend
    captured = []
    install_template_start(fit_cops, template_start, fit_start, restart_state, captured)
    sys.argv = [str(Path(fit_cops.__file__))] + fit_args
    t0 = time.perf_counter()
    try:
        fit_cops.main()
        if captured and not args.no_save_restart and fit_cops.ACTIVE_MANIFEST is not None:
            out = Path(fit_cops.ACTIVE_MANIFEST_PATH).parent
            elements = fit_cops.read_elements(fit_cops.ACTIVE_MANIFEST["inputs"]["base_elements"]["path"])
            state = rs.from_fit(fit_cops.P, captured[-1], out, fit_cops.ACTIVE_MANIFEST, elements)
            if restart_state is not None:
                state.meta["history"] = list(restart_state.meta.get("history", [])) + state.meta["history"]
                state.meta["started_from"] = restart_info["loaded"]
            restart_info["saved"] = rs.write(out / "fit_restart.npz", state)
            print("gpu_fit: restart file %s" % (out / "fit_restart.npz"), flush=True)
        if not args.no_figures and fit_cops.ACTIVE_MANIFEST is not None:
            write_figures(fit_cops, fit_args, figures_info)
        # A last save so the final call counts and times are recorded.
        save_with_backend()
    except BaseException as exc:
        state = (fit_cops.ACTIVE_MANIFEST or {}).get("states", {})
        if fit_cops.ACTIVE_MANIFEST is not None and state.get("process") == "running":
            state["process"] = "failed"
            fit_cops.ACTIVE_MANIFEST["updated_utc"] = fit_cops.utc_now()
            fit_cops.ACTIVE_MANIFEST["error"] = str(exc) or exc.__class__.__name__
            save_with_backend()
        raise
    finally:
        stats = gpu_model.REGISTRY.stats()
        print("gpu_fit: %.1f s wall; GPU calls %s; seconds %s" %
              (time.perf_counter() - t0, stats["calls"],
               {k: round(v, 1) for k, v in stats["seconds"].items()}), flush=True)
        gpu_model.REGISTRY.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
