#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/grain_fit_fixed.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/grain_fit_fixed.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: c9bdd755be9f976747a90203b5ed574c2fa41dc21866c8841a9dfbae0916abbb
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: driver = grain_fit.py (grain + NG = RG crossing) with named parameters held (fixed_parameters.py).
# USED FOR: the 12M cycles 1-6 (jobs/grain_ft_cycle*.sbatch), 20M option 2 (grain_ft20.sbatch), 44M cases.
# RUN (from $APP = srcCG/LiouvilleCOPS_cuda): .venv/bin/python <this> --fix c_ft_c1=1 [--fix Tmin=130 ...]
#        --crossing-weight 10 --devices 0,...,7 --precision mixed --record-store <run>/store --max-nfev 4000
#        --restart-from <restart.npz> --template-start 0 -- <run>/base <fit_cops args> --out <new dir>
# OUT: a normal fit directory (fit_restart.npz, fit_manifest.json with model_extension.fixed_parameters, ...).
# ==== end of header ====
"""grain_fit.py with named parameters held at given values (model_ext/fixed_parameters.py; user 2026-11-04).

    .venv/bin/python <this> --fix T75=133 --fix Tmin=133 --fix c_<template>=1 <grain_fit.py arguments>

The unchanged grain_fit.py runs; this driver only wraps grain_crossing.install so that, after the grain problem
is installed, fixed_parameters.install subclasses it. The provenance goes into the same
fit_manifest.json["model_extension"] record, under "fixed_parameters" (driver, module, SHA256, values);
its "command" is this driver's full command (so grain_repro.sbatch reruns it), "grain_fit_command" the inner one.
"""

import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    command = [str(Path(__file__).resolve())] + sys.argv[1:]   # the full command, for grain_repro.sbatch
    argv, fixed = sys.argv[1:], {}
    while "--fix" in argv[:argv.index("--")]:
        i = argv.index("--fix")
        name, value = argv[i + 1].split("=")
        fixed[name] = float(value)
        del argv[i:i + 2]
    import fixed_parameters
    import grain_crossing
    import grain_fit
    original = grain_crossing.install

    def install(fit_cops, weight, record):
        result = original(fit_cops, weight, record)
        sub = {"driver": str(Path(__file__).resolve()), "driver_sha256": sha256(__file__),
               "module": str(Path(fixed_parameters.__file__).resolve()),
               "module_sha256": sha256(fixed_parameters.__file__)}
        record["fixed_parameters"] = sub
        record["grain_fit_command"] = record.get("command")
        record["command"] = command
        fixed_parameters.install(fit_cops, fixed, sub)
        return result
    grain_crossing.install = install
    sys.argv = [str(HERE / "grain_fit.py")] + argv
    return grain_fit.main()


if __name__ == "__main__":
    sys.exit(main())
