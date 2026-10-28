#!/usr/bin/env python3
"""gpu_fit.py with the dust-grain term in NG and the NG = RG crossing constraint (model_ext/grain_crossing.py).

    .venv/bin/python <this> [--crossing-weight MU] [--check-jacobian] <gpu_fit.py options> -- <fit_cops args>

Everything except --crossing-weight (default 10, user 2026-10-27) and --check-jacobian goes to the unchanged
tools/gpu_fit.py. The fitter (fit_cops.py) and the git tree are not modified: this driver loads the reference
fitter, installs the extension (Problem subclass, least_squares and write_fit wrappers) and runs gpu_fit.main().
--check-jacobian: before fitting, compare the analytic Jacobian columns of the two grain parameters and the
crossing rows with central finite differences at the start (printed, recorded in the manifest).
The extension's provenance is written by the running process into fit_manifest.json["model_extension"].
"""

import hashlib
import sys
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(APP), str(APP / "tools"), str(HERE)]


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    argv = sys.argv[1:]
    weight, check = 10.0, False
    if "--crossing-weight" in argv[:argv.index("--")]:
        i = argv.index("--crossing-weight")
        weight = float(argv[i + 1])
        del argv[i:i + 2]
    if "--check-jacobian" in argv[:argv.index("--")]:
        argv.remove("--check-jacobian")
        check = True

    from cops_cuda import reference
    reference.bind()
    fit_cops = reference.load("fit_cops")
    import grain_crossing
    record = {"driver": str(Path(__file__).resolve()), "driver_sha256": sha256(__file__),
              "module": str(Path(grain_crossing.__file__).resolve()),
              "module_sha256": sha256(grain_crossing.__file__), "command": sys.argv}
    grain_crossing.install(fit_cops, weight, record)

    original_save = fit_cops.save_active_manifest

    def save():
        if fit_cops.ACTIVE_MANIFEST is not None:
            fit_cops.ACTIVE_MANIFEST["model_extension"] = record
        original_save()
    fit_cops.save_active_manifest = save

    if check:
        original_fit_one = None

        def checked(lam, x_start=None):
            P = fit_cops.P
            x = np.array(P.x0, float)
            Jd = P.data_jacobian(x)
            out = {}
            for k in (len(x) - 2, len(x) - 1):
                h = max(abs(x[k]) * 1e-6, 1e-6)
                xp, xm = x.copy(), x.copy()
                xp[k] += h
                xm[k] = max(x[k] - h, 0.0)
                fd = (P.data_residuals(xp) - P.data_residuals(xm)) / (xp[k] - xm[k])
                scale = max(np.abs(fd).max(), 1e-300)
                out[P.names[k]] = float(np.abs(fd - Jd[:, k]).max() / scale)
            print("grain_fit: Jacobian check (max |analytic - central FD| / max|FD|): %s" % out, flush=True)
            record["jacobian_check"] = out
            return original_fit_one(lam, x_start)

        import gpu_fit
        original_install = gpu_fit.install_template_start

        def install_template_start(fc, *a, **k):
            nonlocal original_fit_one
            original_install(fc, *a, **k)
            original_fit_one = fc.fit_one
            fc.fit_one = checked
        gpu_fit.install_template_start = install_template_start

    import gpu_fit
    sys.argv = [str(APP / "tools" / "gpu_fit.py")] + argv
    return gpu_fit.main()


if __name__ == "__main__":
    sys.exit(main())
