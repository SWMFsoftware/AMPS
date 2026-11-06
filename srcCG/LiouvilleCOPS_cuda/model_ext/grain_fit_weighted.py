#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/grain_fit_weighted.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/grain_fit_weighted.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: cf0020b9c4fc0a97ead8fc00ec251519f9aa6317ee3ff0c40b52350f425a0fdd
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: driver = grain_fit.py with larger weights for the last N data points (late_weights.py) and, optionally,
#   held parameters (--fix, as grain_fit_fixed.py). Without --weight-last (or N = 0) it is grain_fit_fixed.py.
# USED FOR: the 12M late-weights fit (job 122236, jobs/grain_late_weights.sbatch: last 5 points x10, NG and RG).
# RUN: .venv/bin/python <this> --weight-last 5 --weight-factor 10 --weight-profile constant|linear
#        --weight-gauges ng,rg [--fix NAME=VALUE ...] <grain_fit.py arguments>
# OUT: a normal fit directory; costs/restart unweighted; the weighted objective in
#   fit_manifest.json["model_extension"]["late_weights"].
# ==== end of header ====
"""grain_fit.py with larger weights for the last N data points (model_ext/late_weights.py; user 2026-11-06) and,
optionally, named parameters held at given values (model_ext/fixed_parameters.py, as grain_fit_fixed.py).

    .venv/bin/python <this> --weight-last 5 --weight-factor 10 [--weight-profile constant|linear]
        [--weight-gauges ng,rg] [--fix NAME=VALUE ...] <grain_fit.py arguments>

* --weight-last N: the last N used points (time order) of each selected gauge get weight W (default 0: off);
* --weight-factor F: W of the last point (W = 1 is the fitter's own weight; default 10);
* --weight-profile constant (W = F for all N, default) or linear (W from 1 to F over the N points);
* --weight-gauges: ng, rg or ng,rg (default ng,rg).
W multiplies the point's cost term (residual x sqrt(W)). Reported costs, rms and the restart file are unweighted;
the weighted objective is in fit_manifest.json["model_extension"]["late_weights"].
The unchanged grain_fit.py runs; this driver wraps grain_crossing.install: the weights are installed first (inside
the crossing wrapper, so the crossing penalty sees unweighted residuals), then the grain problem, then the held
parameters. fit_manifest.json["model_extension"]["command"] is this driver's full command.
"""

import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def pop_option(argv, name, convert=str, default=None):
    head = argv[:argv.index("--")]
    if name in head:
        i = argv.index(name)
        value = convert(argv[i + 1])
        del argv[i:i + 2]
        return value
    return default


def main():
    command = [str(Path(__file__).resolve())] + sys.argv[1:]
    argv = sys.argv[1:]
    if "--" not in argv:
        sys.exit("grain_fit_weighted: need '--' before the fit_cops arguments")
    fixed = {}
    while "--fix" in argv[:argv.index("--")]:
        i = argv.index("--fix")
        name, value = argv[i + 1].split("=")
        fixed[name] = float(value)
        del argv[i:i + 2]
    cfg = {"last": pop_option(argv, "--weight-last", int, 0),
           "factor": pop_option(argv, "--weight-factor", float, 10.0),
           "profile": pop_option(argv, "--weight-profile", str, "constant"),
           "gauges": [g.strip().upper() for g in pop_option(argv, "--weight-gauges", str, "ng,rg").split(",")]}
    import fixed_parameters
    import grain_crossing
    import grain_fit
    import late_weights
    original = grain_crossing.install

    def install(fit_cops, weight, record):
        lw = {"module": str(Path(late_weights.__file__).resolve()), "module_sha256": sha256(late_weights.__file__)}
        if cfg["last"] > 0:
            late_weights.install(fit_cops, cfg, lw)
        else:
            lw["off"] = True
        result = original(fit_cops, weight, record)
        record["late_weights"] = lw
        record["grain_fit_command"] = record.get("command")
        record["command"] = command
        record["driver_weighted"] = {"driver": str(Path(__file__).resolve()), "driver_sha256": sha256(__file__)}
        if fixed:
            sub = {"module": str(Path(fixed_parameters.__file__).resolve()),
                   "module_sha256": sha256(fixed_parameters.__file__)}
            record["fixed_parameters"] = sub
            fixed_parameters.install(fit_cops, fixed, sub)
        return result
    grain_crossing.install = install
    sys.argv = [str(HERE / "grain_fit.py")] + argv
    return grain_fit.main()


if __name__ == "__main__":
    sys.exit(main())
