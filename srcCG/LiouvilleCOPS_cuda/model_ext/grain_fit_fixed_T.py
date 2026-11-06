#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/grain_fit_fixed_T.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/grain_fit_fixed_T.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 46bd99b3db6bc7bc11c3783da0a2723333759d72d19077552e79cf07a710c31a
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: driver = grain_fit.py with some temperatures held (fixed_temperature.py; older form of grain_fit_fixed.py).
# USED FOR: the 12M cliff >=60 fit with T75 = Tmin = 133 K (job 121932, model_ext/grain_fixT.sbatch).
# RUN: .venv/bin/python <this> --fix-T T75=133 --fix-T Tmin=133 <grain_fit.py arguments>
# ==== end of header ====
"""grain_fit.py with some analytical temperatures fixed (model_ext/fixed_temperature.py; user 2026-11-04).

    .venv/bin/python <this> --fix-T T75=133 --fix-T Tmin=133 <grain_fit.py arguments>

The unchanged grain_fit.py runs; this driver only wraps grain_crossing.install so that, after the grain problem
is installed, fixed_temperature.install subclasses it. The provenance goes into the same
fit_manifest.json["model_extension"] record, under "fixed_temperature" (driver, module, SHA256, values);
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
    while "--fix-T" in argv[:argv.index("--")]:
        i = argv.index("--fix-T")
        name, value = argv[i + 1].split("=")
        fixed[name] = float(value)
        del argv[i:i + 2]
    import fixed_temperature
    import grain_crossing
    import grain_fit
    original = grain_crossing.install

    def install(fit_cops, weight, record):
        result = original(fit_cops, weight, record)
        sub = {"driver": str(Path(__file__).resolve()), "driver_sha256": sha256(__file__),
               "module": str(Path(fixed_temperature.__file__).resolve()),
               "module_sha256": sha256(fixed_temperature.__file__)}
        record["fixed_temperature"] = sub
        record["grain_fit_command"] = record.get("command")
        record["command"] = command
        fixed_temperature.install(fit_cops, fixed, sub)
        return result
    grain_crossing.install = install
    sys.argv = [str(HERE / "grain_fit.py")] + argv
    return grain_fit.main()


if __name__ == "__main__":
    sys.exit(main())
