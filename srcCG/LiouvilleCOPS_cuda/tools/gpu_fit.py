#!/usr/bin/env python3
"""Run the unchanged ``fit_cops`` fitter with GPU pressure/Jacobian evaluation (CUDA env).

Usage::

    .venv/bin/python tools/gpu_fit.py [--device N] [--mode tree|sequential] -- <fit_cops.py arguments>

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
* the fit manifest gains a ``compute_backend`` section, written on every
  manifest save, which states the CUDA backend, the device and dataset
  identities, the binding events and the call counts and times. The
  manifest's ``command`` field lists the fit_cops arguments (argparse reads
  ``sys.argv``); ``compute_backend.driver_command`` is the actual process
  command.

Refused: several ``--lambda`` values. The fitter runs an L-curve in forked
processes, and a forked child cannot use the parent's CUDA context. Run one
λ per process instead.

The process-failure handling of ``fit_cops``'s ``__main__`` block is reproduced:
after the output directory exists, a failure is recorded as
``states.process = "failed"`` in the manifest.
"""

import argparse
from pathlib import Path
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


def main():
    argv = sys.argv[1:]
    if "--" not in argv:
        raise SystemExit("usage: gpu_fit.py [--device N] [--mode tree|sequential] -- <fit_cops args>")
    split = argv.index("--")
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--device", type=int, default=0)
    parser.add_argument("--mode", choices=("tree", "sequential"), default="tree",
                        help="tree: production; sequential: bitwise CPU order (slow, verification)")
    parser.add_argument("--records-per-block", type=int, default=4096)
    args = parser.parse_args(argv[:split])
    fit_args = argv[split + 1:]
    if lambda_count(fit_args) != 1:
        raise SystemExit("gpu_fit.py runs exactly one --lambda value per process "
                         "(the fitter forks for an L-curve; CUDA contexts do not survive fork)")

    from cops_cuda import gpu_model, manifest, provenance, reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    gpu_model.REGISTRY.configure(args.device, args.mode, args.records_per_block)
    fit_cops.Model = gpu_model.gpu_model_class()

    started = provenance.utc_now()
    software = manifest.backend_manifest("cuda_fp64_v1", args.device)
    original_save = fit_cops.save_active_manifest

    def save_with_backend():
        if fit_cops.ACTIVE_MANIFEST is not None:
            fit_cops.ACTIVE_MANIFEST["compute_backend"] = {
                "kind": "cuda_fp64_v1 via LiouvilleCOPS_cuda gpu_model (fit_cops.Model replaced "
                        "in-process; fit_cops.py unmodified)",
                "driver": provenance.file_identity(Path(__file__), "tool"),
                "driver_command": [sys.executable, str(Path(__file__).resolve())] + argv,
                "fit_cops": provenance.file_identity(Path(fit_cops.__file__), "reference"),
                "started_utc": started,
                "validation_policy": "VALIDATION.md (pressure <= 512 eps M, Jacobian <= 1024 eps M_J)",
                "software": software,
                "registry": gpu_model.REGISTRY.stats(),
            }
        original_save()

    fit_cops.save_active_manifest = save_with_backend
    sys.argv = [str(Path(fit_cops.__file__))] + fit_args
    t0 = time.perf_counter()
    try:
        fit_cops.main()
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
