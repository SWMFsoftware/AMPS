# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/fixed_parameters.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/fixed_parameters.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 70d414491d8a197dfc2ddea5aa9d683aa694ee6d0d61a8347e2bfe0b5a128606
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: module. Holds any named fit parameters at given values (temperatures, template amplitudes, ...): they are
#   removed from the optimizer vector and put back by unpack(); the regularization term, its Jacobian and D are
#   evaluated on the full vector (fix of 2026-11-06). Generalizes fixed_temperature.py.
# USED BY: grain_fit_fixed.py and grain_fit_weighted.py (--fix NAME=VALUE). Not run directly.
# USED FOR: 20M fit from the 12M fine-tuned state (grain_ft20.sbatch), the 44M cases (grain_44m_case.sbatch), all 12M
#   fit/fine-tuning cycles (jobs/grain_ft_cycle*.sbatch: correction templates c_ft_c<k> = 1, Tmin/Tmax/T75 held), the
#   late-weights fit (jobs/grain_late_weights.sbatch).
# CALL: fixed_parameters.install(fit_cops, {"Tmin": 130.0, "c_ft_c1": 1.0}, record) after grain_crossing.install.
# NOTE: needs free (not ordered, not --fix-temperature) temperatures; a held value must be inside its bounds.
# ==== end of header ====
"""Hold named fit parameters at given values; the others stay free (generalization of fixed_temperature.py).

User request 2026-11-04: "run option 2 on 20M" (a 20M fit starting from the 12M fine-tuned state: the 12M
fine-tuning activity correction carried to 20M as a fixed template, amplitude held at 1; T75 = Tmin = 133 K held as in
the 12M fit).

Local, untracked extension: fit_cops.py, gpu_fit.py, grain_crossing.py are not edited. ``install`` subclasses the
Problem class installed at that moment (the GrainProblem of grain_crossing.py) and removes the named parameters
from the optimizer vector (names, x0, lo, hi, Jacobian columns); ``unpack`` puts them back at their values, so every
evaluation, writer, restart file and the fine tuning see them. Values are in optimizer units: K for Tmax/T75/Tmin,
the amplitude for a template ``c_<name>`` (its activity is amplitude * A_scale * column), kappa as is.
Temperature names need the plain (not ordered, not --fix-temperature) parameterization. The grain start A is
recomputed with the held values, as grain_crossing.py defines it (relevant for the default start only).
"""

import numpy as np

TNAMES = ("Tmax", "T75", "Tmin")


def install(fit_cops, fixed, record):
    """fixed: {parameter name: value}. record (dict) receives the provenance."""
    if not fixed:
        raise ValueError("no parameter to hold")
    base = fit_cops.Problem
    record.update(fixed=dict(fixed))

    class FixedParameterProblem(base):
        def __init__(self, M, el, coef, args):
            if any(k in TNAMES for k in fixed) and (args.fix_temperature or args.ordered_temperature):
                raise SystemExit("fixed_parameters: temperatures need free, not ordered temperatures")
            lo_t, hi_t = args.tbounds
            for k, v in fixed.items():
                if k in TNAMES and not lo_t <= v <= hi_t:
                    raise SystemExit("fixed_parameters: %s = %g K outside --tbounds %g %g" % (k, v, lo_t, hi_t))
            self._fp_full = None
            super().__init__(M, el, coef, args)
            names = list(self.names)
            missing = [k for k in fixed if k not in names]
            if missing:
                raise SystemExit("fixed_parameters: %s not among the parameters %s" % (missing, names))
            if len(fixed) == len(names):
                raise SystemExit("fixed_parameters: every parameter would be held")
            full = np.array(self.x0, float)
            self._fp_prior = full.copy()    # the base prior x0 (full vector), for the Tikhonov term and D
            lo, hi = np.asarray(self.lo, float), np.asarray(self.hi, float)
            drop = []
            for k, v in fixed.items():
                i = names.index(k)
                if not lo[i] <= v <= hi[i]:
                    raise SystemExit("fixed_parameters: %s = %g outside its bounds [%g, %g]" % (k, v, lo[i], hi[i]))
                full[i] = v
                drop.append(i)
            self._fp_keep = np.array([i for i in range(len(full)) if i not in drop])
            self._fp_full = full
            self.names = [names[i] for i in self._fp_keep]
            self.x0, self.lo, self.hi = full[self._fp_keep], lo[self._fp_keep], hi[self._fp_keep]
            if hasattr(self, "et_end"):     # grain start A with the held values (grain_crossing.py definition)
                a, b, T, kappa = self.unpack(self.x0)
                pNG, _ = M.evaluate(a, T, kappa, b=b if self.nb else None,
                                    enhancement_basis=self.enhancement_basis, enhancement_scale=self.enhancement_scale)
                i_end = int(np.nonzero(self.useNG & (self.et_fit == self.et_end))[0][0])
                jump = max(float(M.refNG[M.rows][i_end] - pNG[i_end]), 0.0)
                A_UNIT = getattr(__import__("grain_crossing"), "A_UNIT")
                self.x0[-2] = jump / A_UNIT
                record.update(default_start_grain_A_uPa=jump / A_UNIT, default_start_T=[float(v) for v in T])
            record.update(names=list(self.names))

        def _fp_expand(self, x):
            if self._fp_full is None:       # during the base constructor
                return x
            xf = self._fp_full.copy()
            xf[self._fp_keep] = x
            return xf

        def unpack(self, x):
            return super().unpack(self._fp_expand(np.asarray(x, float)))

        # Fix 2026-11-06 (cycle 5 of the Tmin130 cycles, job 122118, crashed in reg_jacobian: "could not broadcast
        # (31,31) into (31,30)"): reg_slice indexes the FULL parameter vector, so the Tikhonov term, its Jacobian
        # and the diagnostic D are evaluated on the full vector (held values filled in, prior = the base x0) and the
        # Jacobian columns of the held parameters dropped. Before the fix, held parameters inside reg_slice shifted
        # it: D was misaligned in every fit with held templates; with lambda = 0 the cost, the Jacobian and the
        # optimizer path were unaffected (the term is sqrt(0) * ...).
        def regularization(self, x):
            if self._fp_full is None:
                return super().regularization(x)
            xf = self._fp_expand(np.asarray(x, float))
            return np.sqrt(self.lam) * (self.Lt @ (xf[self.reg_slice] - self._fp_prior[self.reg_slice]))

        def reg_jacobian(self):
            if self._fp_full is None:
                return super().reg_jacobian()
            R = np.zeros((self.Lt.shape[0], len(self._fp_full)))
            R[:, self.reg_slice] = np.sqrt(self.lam) * self.Lt
            return R[:, self._fp_keep]

        def deviation(self, x):
            if self._fp_full is None:
                return super().deviation(x)
            xf = self._fp_expand(np.asarray(x, float))
            return np.linalg.norm(self.Lt @ (xf[self.reg_slice] - self._fp_prior[self.reg_slice]))

        def data_jacobian(self, x):
            J = super().data_jacobian(x)    # the base builds every column
            return J if self._fp_full is None else J[:, self._fp_keep]

        def activity_linear_and_derivative(self, x, parameter_index):
            # fit_cops --check-jacobian only; the index is in the reduced vector, mapped to the full one
            k = parameter_index if self._fp_full is None else int(self._fp_keep[parameter_index])
            a, b, _, _ = self.unpack(x)
            total = self.M.Y @ a
            if self.nb:
                total = total + self.enhancement_scale * (self.enhancement_basis @ b)
            if k < self.na:
                return total, self.ascale * self.M.Y[:, self.fit_j[k]]
            if k < self.na + self.nb:
                return total, self.enhancement_scale * self.enhancement_basis[:, k - self.na]
            return total, None

    fit_cops.Problem = FixedParameterProblem
    return FixedParameterProblem
