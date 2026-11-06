# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/fixed_temperature.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/fixed_temperature.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: e0e204c061c2c6bd1c6a7b5e6ada31dd1435dc7675fd0d4798bccb58a44cb3c6
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: module (older, temperatures only). Holds some of Tmax, T75, Tmin at given values; superseded by
#   fixed_parameters.py, kept because the 12M T75 = Tmin = 133 K fit (job 121932) ran with it.
# USED BY: grain_fit_fixed_T.py (--fix-T NAME=VALUE), model_ext/grain_fixT.sbatch.
# ==== end of header ====
"""Fix some of the analytical temperatures (Tmax, T75, Tmin) at given values; the others stay free.

User request 2026-11-04: "run the fit with Tmin and T75 fixed at 133 K" (cliff >=60 deg template, free
temperatures, dust grain + crossing, default start; Tmax free within --tbounds).

Local, untracked extension like grain_crossing.py: fit_cops.py, gpu_fit.py and grain_crossing.py are not edited.
``install`` subclasses the Problem class installed at that moment (the GrainProblem of grain_crossing.py) and
removes the fixed temperatures from the optimizer vector (names, x0, lo, hi, Jacobian columns). ``unpack`` puts
them back at their fixed values, so every evaluation, writer, restart file (its ``T`` array) and the fine tuning
see Tmax (fitted), T75 and Tmin (fixed). The grain start A (the observed jump minus the surface model at the
start) is recomputed with the fixed temperatures, as grain_crossing.py defines it.
Only the plain (not ordered, not --fix-temperature) temperature parameterization is supported.
"""

import numpy as np

TNAMES = ("Tmax", "T75", "Tmin")


def install(fit_cops, fixed, record):
    """fixed: {name: value [K]} with names from TNAMES. record (dict) receives the provenance."""
    bad = [k for k in fixed if k not in TNAMES]
    if bad or not fixed or len(fixed) == len(TNAMES):
        raise ValueError("fix one or two of %s (use --fix-temperature for all three): %s" % (TNAMES, fixed))
    base = fit_cops.Problem
    record.update(fixed=dict(fixed), free=[k for k in TNAMES if k not in fixed])

    class FixedTemperatureProblem(base):
        def __init__(self, M, el, coef, args):
            if args.fix_temperature or args.ordered_temperature:
                raise SystemExit("fixed_temperature: needs free, not ordered temperatures")
            lo, hi = args.tbounds
            for k, v in fixed.items():
                if not lo <= v <= hi:
                    raise SystemExit("fixed_temperature: %s = %g K outside --tbounds %g %g" % (k, v, lo, hi))
            self._full = None
            super().__init__(M, el, coef, args)
            names = list(self.names)
            k0 = names.index("Tmax")
            if names[k0:k0 + 3] != list(TNAMES):
                raise SystemExit("fixed_temperature: unexpected parameter layout %s" % names[k0:k0 + 3])
            full = np.array(self.x0, float)
            for i, n in enumerate(TNAMES):
                if n in fixed:
                    full[k0 + i] = fixed[n]
            drop = [k0 + i for i, n in enumerate(TNAMES) if n in fixed]
            self._keep = np.array([i for i in range(len(full)) if i not in drop])
            self._full = full
            self.names = [names[i] for i in self._keep]
            self.x0, self.lo, self.hi = full[self._keep], np.asarray(self.lo)[self._keep], np.asarray(self.hi)[self._keep]
            # grain start with the fixed temperatures (grain_crossing.py: A = observed jump at t_end, clipped at 0)
            if hasattr(self, "et_end"):
                a, b, T, kappa = self.unpack(self.x0)
                pNG, _ = M.evaluate(a, T, kappa, b=b if self.nb else None,
                                    enhancement_basis=self.enhancement_basis, enhancement_scale=self.enhancement_scale)
                i_end = int(np.nonzero(self.useNG & (self.et_fit == self.et_end))[0][0])
                jump = max(float(M.refNG[M.rows][i_end] - pNG[i_end]), 0.0)
                A_UNIT = getattr(__import__("grain_crossing"), "A_UNIT")
                self.x0[-2] = jump / A_UNIT
                record.update(grain_start_with_fixed_T={"grain_A_uPa": jump / A_UNIT,
                                                        "observed_ng_at_t_end_Pa": float(M.refNG[M.rows][i_end]),
                                                        "surface_model_ng_at_t_end_Pa": float(pNG[i_end])})
            record.update(start_T=[float(v) for v in T] if hasattr(self, "et_end") else None,
                          names=list(self.names))

        def expand(self, x):
            if self._full is None:          # during the base constructor
                return x
            xf = self._full.copy()
            xf[self._keep] = x
            return xf

        def unpack(self, x):
            return super().unpack(self.expand(np.asarray(x, float)))

        def data_jacobian(self, x):
            J = super().data_jacobian(x)    # the base builds all three temperature columns
            return J if self._full is None else J[:, self._keep]

    fit_cops.Problem = FixedTemperatureProblem
    return FixedTemperatureProblem
