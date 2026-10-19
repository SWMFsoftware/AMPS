"""Pressure/Jacobian backend contract (ROADMAP C02 item 1; Python 3.8 compatible).

The CPU reference ``fit_cops.Model.evaluate`` defines what a backend computes.
This module makes that interface explicit, so the CPU and CUDA
implementations are interchangeable behind it and compared through it:

    backend.evaluate(ParameterState, ColumnRequest) -> EvaluationResult

* ``ParameterState`` is the physical parameter vector: the 25 legacy
  activity coefficients ``a``, the analytical temperatures
  ``T = (Tmax, T75, Tmin)``, the CO2 scale ``kappa`` and the enhancement
  coefficients ``b``. Mapping the optimizer vector ``x`` to it (activity
  scaling, ordered temperatures, fixed parameters) stays in the CPU
  ``Problem.unpack``.
* ``ColumnRequest`` names the analytic-Jacobian columns wanted. Only those
  are computed (ROADMAP C06 item 3).
* ``EvaluationResult`` holds FP64 ``pNG`` and ``pRG`` for every loaded point,
  in the backend's row order, and the requested columns as
  ``(dpNG/dq, dpRG/dq)`` pairs.

Exact existing order
--------------------
The keys and their order are those ``Model.evaluate`` inserts into its ``J``
dictionary:

    ("a", j) for j in basis_index, in the caller's order
    ("b", k) for k in enhancement_index, in the caller's order
    ("T", 0), ("T", 1), ("T", 2)    = d/dTmax, d/dT75, d/dTmin
    ("kappa", 0)

``fit_cops.Problem.data_jacobian`` turns these into optimizer columns. It
applies the activity scale, the ordered-temperature chain rule, the log and
the per-gauge weights. That code is not duplicated here: ``ModelFacade``
lets the unmodified ``Problem`` run on any backend. Residuals, weights,
regularization, SciPy TRF and every diagnostic therefore remain the
reference code, as ROADMAP §3.3 and C06 item 6 require.

Joint families and identities (ROADMAP rev. 2, C02 extensions)
---------------------------------------------------------------
The primary families B29/N30/C30/NC31 (free baseline, temperatures and κ,
plus nonnegative template amplitudes) use this contract unchanged. Columns
come from ``ColumnRequest.from_problem``: a(fit_j), templates, T, κ, 31 for
NC31.
* ``problem_contract(problem)`` identifies everything a backend must match:
  the parameter layout, the template names and the basis digest, the
  enhancement scale, and the ORIGINAL prior (run coefficients, x0,
  regularization operator), kept distinct from any initialization.
* ``install`` refuses a backend whose bound basis digest or scale differs from
  the problem's: a stale or reordered template set.
* ``ParameterState.to_dict``/``from_dict`` and ``ColumnRequest.as_dict``/``from_dict``
  round-trip exactly (doubles as ``float.hex``). Fixtures and restart records
  can therefore name states and requests without loss.

Static versus per-call data
---------------------------
The mesh, the transport kernels, the model settings and the per-facet
enhancement basis are bound when a backend is constructed. Only the parameter
values change between calls. ``Model.evaluate`` accepts the enhancement basis
on every call; the facade checks that the basis passed is the bound one, so
the two cannot diverge.
"""

import abc
import collections
import hashlib

import numpy as np

from . import provenance, reference
from .errors import BackendError, CacheFormatError, IdentityError, UnsupportedConfiguration


CONTRACT_VERSION = 1
NBASIS = 25
TEMPERATURE_KEYS = (("T", 0), ("T", 1), ("T", 2))
KAPPA_KEY = ("kappa", 0)


# -------------------------------------------------------------------- settings
class ModelSettings:
    """Run-level model choices that fix the evaluator's formula.

    Values are read from the ``COPS_model.dat`` header with the reference's
    own parsers (``fit_cops.header_value`` and ``fit_cops.sensitivities``).
    The constants (Boltzmann constant, masses, the 75-degree factor, the
    pressure floor) are taken from the reference module, not copied.
    """

    FIELDS = ("sza_mode", "heliocentric_model", "temperature_ng", "rg_model",
              "sensitivities", "analytical_T0")

    def __init__(self, sza_mode, heliocentric_model, temperature_ng, rg_model,
                 sensitivities, analytical_T0):
        self.sza_mode = str(sza_mode)
        self.heliocentric_model = str(heliocentric_model)
        self.temperature_ng = float(temperature_ng)
        self.rg_model = str(rg_model)
        self.sensitivities = tuple(float(v) for v in sensitivities)
        self.analytical_T0 = tuple(float(v) for v in analytical_T0)
        if len(self.sensitivities) != 2 or len(self.analytical_T0) != 3:
            raise ValueError("sensitivities need (H2O, CO2) and T0 needs (Tmax, T75, Tmin)")

    @classmethod
    def from_header(cls, header):
        """Parse exactly as ``fit_cops.Model.__init__`` does, without loading kernels."""
        fit_cops = reference.load("fit_cops")
        try:
            # The reference refuses these modes; so does every backend.
            if fit_cops.header_value(header, "TemperatureModel") != "analytical":
                raise UnsupportedConfiguration("the run must use TemperatureModel=analytical")
            if fit_cops.header_value(header, "FlightTimeRotation") != "none":
                raise UnsupportedConfiguration(
                    "the run must use FlightTimeRotation=none (version 1 kernels)")
            return cls(fit_cops.header_value(header, "SZAMode"),
                       fit_cops.header_value(header, "HeliocentricModel"),
                       fit_cops.header_value(header, "TemperatureNG"),
                       fit_cops.header_value(header, "RGPressureModel"),
                       fit_cops.sensitivities(header),
                       [fit_cops.header_value(header, "Analytical" + k)
                        for k in fit_cops.TNAMES])
        except SystemExit as exc:
            # The reference parser exits on a missing key (CLI behaviour);
            # inside a library that becomes an ordinary refusal.
            raise CacheFormatError("model table header: %s" % exc) from None

    @classmethod
    def from_model(cls, model):
        """Settings of an already constructed ``fit_cops.Model``."""
        return cls(model.SZAMode, model.helio, model.TNG, model.RGModel, model.S, model.T0)

    def as_dict(self):
        return {name: (list(value) if isinstance(value, tuple) else value)
                for name, value in ((f, getattr(self, f)) for f in self.FIELDS)}

    def __eq__(self, other):
        return isinstance(other, ModelSettings) and self.as_dict() == other.as_dict()

    def __repr__(self):
        return "ModelSettings(%s)" % self.as_dict()


def reference_constants():
    """Physical and numerical constants of the reference evaluator, as used."""
    fit_cops = reference.load("fit_cops")
    return {"KB": float(fit_cops.KB), "AMU": float(fit_cops.AMU), "AU": float(fit_cops.AU),
            "MASS": {k: float(v) for k, v in fit_cops.MASS.items()},
            "KT": float(fit_cops.KT), "COS75": float(fit_cops.COS75),
            "FLOOR": float(fit_cops.FLOOR), "NBASIS": int(fit_cops.NBASIS)}


# ------------------------------------------------------------------- parameters
def _float64_vector(value, length, label):
    array = np.array(value, dtype=np.float64, copy=True).reshape(-1)
    if length is not None and array.shape != (length,):
        raise ValueError("%s must have %d values; got %d" % (label, length, array.size))
    if not np.isfinite(array).all():
        raise ValueError("%s contains non-finite values" % label)
    array.setflags(write=False)
    return array


class ParameterState:
    """Physical parameters of one evaluation (immutable FP64 copies)."""

    def __init__(self, a, T, kappa, b=None):
        self.a = _float64_vector(a, NBASIS, "a")
        self.T = _float64_vector(T, 3, "T")
        self.kappa = float(kappa)
        if not np.isfinite(self.kappa):
            raise ValueError("kappa is not finite")
        self.b = None if b is None else _float64_vector(b, None, "b")

    @classmethod
    def from_problem(cls, problem, x):
        """Map an optimizer vector with the reference ``Problem.unpack``."""
        a, b, T, kappa = problem.unpack(np.asarray(x, dtype=np.float64))
        return cls(a, T, kappa, b if problem.nb else None)

    def to_dict(self):
        """Exact, JSON-safe form: every double as ``float.hex`` (bitwise round trip)."""
        def hexes(values):
            return [float(v).hex() for v in values]
        return {"a": hexes(self.a), "T": hexes(self.T), "kappa": float(self.kappa).hex(),
                "b": None if self.b is None else hexes(self.b)}

    @classmethod
    def from_dict(cls, value):
        def floats(values):
            return [float.fromhex(v) for v in values]
        return cls(floats(value["a"]), floats(value["T"]), float.fromhex(value["kappa"]),
                   None if value.get("b") is None else floats(value["b"]))

    def digest(self):
        """Exact content key, for the C06 cache of repeated evaluations."""
        digest = hashlib.sha256()
        for part in (self.a, self.T, np.array([self.kappa]),
                     np.empty(0) if self.b is None else self.b):
            digest.update(np.ascontiguousarray(part, dtype="<f8").tobytes())
            digest.update(b"|")
        digest.update(b"b" if self.b is not None else b"-")
        return digest.hexdigest()


class ColumnRequest:
    """The analytic-Jacobian columns one evaluation must return."""

    def __init__(self, basis_index=(), enhancement_index=(), temperature=False, kappa=False):
        self.basis_index = tuple(int(j) for j in basis_index)
        self.enhancement_index = tuple(int(k) for k in enhancement_index)
        self.temperature = bool(temperature)
        self.kappa = bool(kappa)
        if any(j < 0 or j >= NBASIS for j in self.basis_index):
            raise ValueError("legacy basis index outside [0, %d)" % NBASIS)
        if len(set(self.basis_index)) != len(self.basis_index):
            raise ValueError("duplicate legacy basis index")
        if any(k < 0 for k in self.enhancement_index) or (
                len(set(self.enhancement_index)) != len(self.enhancement_index)):
            raise ValueError("enhancement indices must be unique and non-negative")

    @classmethod
    def from_problem(cls, problem):
        """Exactly the columns ``Problem.data_jacobian`` reads."""
        return cls(problem.fit_j, range(problem.nb),
                   temperature=not problem.args.fix_temperature,
                   kappa=not problem.args.fix_kappa)

    def keys(self):
        """Column keys in the reference order (see module docstring)."""
        keys = [("a", j) for j in self.basis_index]
        keys += [("b", k) for k in self.enhancement_index]
        if self.temperature:
            keys += list(TEMPERATURE_KEYS)
        if self.kappa:
            keys.append(KAPPA_KEY)
        return keys

    @property
    def empty(self):
        return not self.keys()

    def as_dict(self):
        return {"basis_index": list(self.basis_index),
                "enhancement_index": list(self.enhancement_index),
                "temperature": self.temperature, "kappa": self.kappa}

    @classmethod
    def from_dict(cls, value):
        return cls(value["basis_index"], value["enhancement_index"],
                   value["temperature"], value["kappa"])

    def __eq__(self, other):
        return isinstance(other, ColumnRequest) and self.as_dict() == other.as_dict()


PRESSURES_ONLY = ColumnRequest()


class EvaluationResult:
    """Pressures and requested derivative columns, validated on construction."""

    def __init__(self, pNG, pRG, columns, request):
        self.pNG = np.asarray(pNG, dtype=np.float64)
        self.pRG = np.asarray(pRG, dtype=np.float64)
        self.request = request
        npt = self.pNG.shape[0] if self.pNG.ndim == 1 else -1
        if npt < 0 or self.pRG.shape != (npt,):
            raise BackendError("pressures must be two vectors of equal length")
        self.columns = collections.OrderedDict()
        for key in request.keys():
            if key not in columns:
                raise BackendError("backend did not return requested column %s" % (key,))
            dNG, dRG = columns[key]
            dNG, dRG = np.asarray(dNG, dtype=np.float64), np.asarray(dRG, dtype=np.float64)
            if dNG.shape != (npt,) or dRG.shape != (npt,):
                raise BackendError("column %s has shapes %s/%s; expected (%d,)" %
                                   (key, dNG.shape, dRG.shape, npt))
            self.columns[key] = (dNG, dRG)

    @property
    def npoints(self):
        return self.pNG.shape[0]

    def jacobian_dict(self):
        """The ``J`` dictionary of ``Model.evaluate`` restricted to the request."""
        return collections.OrderedDict(self.columns)

    def matrix(self, gauge):
        """``(npoints, ncolumns)`` derivative matrix of one gauge, in key order."""
        index = {"NG": 0, "RG": 1}[gauge]
        if not self.columns:
            return np.empty((self.npoints, 0))
        return np.column_stack([pair[index] for pair in self.columns.values()])


# ------------------------------------------------------------ problem contract
def array_digest(array):
    return provenance.array_sha256(np.ascontiguousarray(np.asarray(array, dtype=np.float64)))


def problem_contract(problem):
    """Identity of everything in a ``fit_cops.Problem`` that a backend must match.

    The parameter layout (names, bounds, start, scale), the bound per-facet
    template basis and its names, the enhancement scale, and the ORIGINAL prior
    (run coefficients, ``x0`` and the regularization operator ``Lt``). The prior
    is recorded separately from any later initialization: a warm start must
    never redefine it (ROADMAP §3.3). Row selection is identified by the
    model's table rows.
    """
    names = list(problem.names)
    templates = list(getattr(problem, "template_names", []))
    contract = {
        "names": names, "lo": [float(v).hex() for v in problem.lo],
        "hi": [float(v).hex() for v in problem.hi],
        "ascale": float(problem.ascale).hex(),
        "fit_j": [int(j) for j in problem.fit_j], "nh": int(problem.nh), "nt": int(problem.nt),
        "template_names": templates,
        "template_sources": [str(s) for s in getattr(problem, "template_sources", [])],
        "enhancement_scale": float(problem.enhancement_scale).hex(),
        "enhancement_basis_sha256": (array_digest(problem.enhancement_basis)
                                     if problem.nb else None),
        "prior": {"coefficients_sha256": array_digest(problem.coef),
                  "x0_sha256": array_digest(problem.x0),
                  "regularization_operator_sha256": array_digest(problem.Lt),
                  "reg_slice": [problem.reg_slice.start, problem.reg_slice.stop]},
        "rows_sha256": provenance.canonical_sha256([int(r) for r in problem.M.rows]),
        "temperature_free": not problem.args.fix_temperature,
        "ordered_temperature": bool(problem.ordered), "kappa_free": not problem.args.fix_kappa,
    }
    contract["sha256"] = provenance.canonical_sha256(contract)
    return contract


# --------------------------------------------------------------------- backend
class PressureBackend(abc.ABC):
    """Abstract evaluator bound to one dataset (rows, mesh, settings, basis).

    Subclasses implement ``capabilities``, ``_evaluate`` and ``identity``.
    The public ``evaluate`` validates input and output for every backend
    alike, so a subclass cannot skip the checks.
    """

    name = "abstract"

    def __init__(self, rows, settings, nbasis_enhancement):
        self.rows = tuple(int(r) for r in rows)
        self.settings = settings
        self.nbasis_enhancement = int(nbasis_enhancement)

    @abc.abstractmethod
    def capabilities(self):
        """Dict with the supported ``sza_modes``, ``heliocentric_models``,
        ``rg_models``, ``column_kinds`` (subset of a, b, T, kappa) and
        ``enhancement`` (bool)."""

    @abc.abstractmethod
    def _evaluate(self, state, request):
        """Return an ``EvaluationResult`` (inputs already validated)."""

    @abc.abstractmethod
    def identity(self):
        """Manifest section: backend name, contract version, precision, data."""

    def close(self):
        """Release resources; the base backend holds none."""

    def unsupported(self, request=None):
        """Reasons this backend cannot evaluate its bound configuration."""
        caps, reasons = self.capabilities(), []
        settings = self.settings
        for field, allowed in (("sza_mode", "sza_modes"),
                               ("heliocentric_model", "heliocentric_models"),
                               ("rg_model", "rg_models")):
            value = getattr(settings, field)
            if value not in caps[allowed]:
                reasons.append("%s=%s not in %s" % (field, value, sorted(caps[allowed])))
        if self.nbasis_enhancement and not caps["enhancement"]:
            reasons.append("enhancement basis not supported")
        if request is not None:
            kinds = set(key[0] for key in request.keys())
            missing = kinds - set(caps["column_kinds"])
            if missing:
                reasons.append("derivative kinds %s not supported" % sorted(missing))
        return reasons

    def check_supported(self, request=None):
        reasons = self.unsupported(request)
        if reasons:
            raise UnsupportedConfiguration("%s backend: %s" % (self.name, "; ".join(reasons)))

    def evaluate(self, state, request=PRESSURES_ONLY):
        if not isinstance(state, ParameterState):
            raise TypeError("evaluate needs a ParameterState")
        nb = 0 if state.b is None else state.b.size
        if nb != self.nbasis_enhancement:
            raise ValueError("state has %d enhancement coefficients; backend is bound to %d"
                             % (nb, self.nbasis_enhancement))
        if any(k >= self.nbasis_enhancement for k in request.enhancement_index):
            raise ValueError("enhancement column outside the bound basis")
        self.check_supported(request)
        result = self._evaluate(state, request)
        if result.npoints != len(self.rows):
            raise BackendError("backend returned %d points for %d rows" %
                               (result.npoints, len(self.rows)))
        if list(result.columns) != request.keys():
            raise BackendError("backend column order differs from the contract")
        return result


class ReferenceBackend(PressureBackend):
    """The CPU reference ``fit_cops.Model`` behind the contract.

    This is the oracle for C06: the CUDA backend is compared with it on the
    same ``ParameterState``s and ``ColumnRequest``s. It calls
    ``Model.evaluate`` unchanged. Unrequested columns, which the reference
    always computes, are dropped.
    """

    name = "cpu_reference_fit_cops"

    def __init__(self, model, enhancement_basis=None, enhancement_scale=1.0):
        nb = 0 if enhancement_basis is None else int(np.asarray(enhancement_basis).shape[1])
        super().__init__(model.rows, ModelSettings.from_model(model), nb)
        self.model = model
        self.enhancement_basis = enhancement_basis
        self.enhancement_scale = float(enhancement_scale)
        self.enhancement_basis_sha256 = (None if enhancement_basis is None
                                         else array_digest(enhancement_basis))

    @classmethod
    def from_problem(cls, problem):
        """Bind exactly the basis and scale of a ``fit_cops.Problem``."""
        return cls(problem.M, problem.enhancement_basis if problem.nb else None,
                   problem.enhancement_scale)

    def capabilities(self):
        # Every branch of Model: any SZA mode (continuous or binned) and
        # heliocentric model, all three RG models, and all derivative kinds.
        return {"sza_modes": {self.settings.sza_mode},
                "heliocentric_models": {self.settings.heliocentric_model},
                "rg_models": {"equilibrium_sphere", "calibration_srccg", "dsmc_srccg"},
                "column_kinds": {"a", "b", "T", "kappa"}, "enhancement": True}

    def _evaluate(self, state, request):
        kwargs = {}
        if state.b is not None:
            kwargs = {"b": state.b, "enhancement_basis": self.enhancement_basis,
                      "enhancement_scale": self.enhancement_scale}
        if request.empty:
            pNG, pRG = self.model.evaluate(state.a, state.T, state.kappa, **kwargs)
            return EvaluationResult(pNG, pRG, {}, request)
        pNG, pRG, J = self.model.evaluate(
            state.a, state.T, state.kappa, jac=True, basis_index=list(request.basis_index),
            enhancement_index=list(request.enhancement_index), **kwargs)
        return EvaluationResult(pNG, pRG, J, request)

    def identity(self):
        return {"backend": self.name, "contract_version": CONTRACT_VERSION,
                "precision": "float64 (NumPy/SciPy reference)",
                "settings": self.settings.as_dict(), "constants": reference_constants(),
                "rows": len(self.rows), "rows_sha256": provenance.canonical_sha256(list(self.rows)),
                "enhancement_columns": self.nbasis_enhancement,
                "enhancement_scale": self.enhancement_scale,
                "enhancement_basis_sha256": self.enhancement_basis_sha256}


# ---------------------------------------------------------------------- facade
class ModelFacade:
    """Stand-in for ``fit_cops.Model`` inside an unmodified ``fit_cops.Problem``.

    ``Problem`` reads selection data from its model (``rows``, ``useNG``,
    ``refNG``, ``T0``, ``Y``, ...) and calls ``evaluate``. The facade
    forwards the data attributes to the CPU model and routes ``evaluate``
    through a backend. Install it with ``install``.

    One deliberate restriction: with ``jac=True`` the facade returns the
    temperature and kappa columns only if the bound ``ColumnRequest`` asks
    for them. ``Problem.data_jacobian`` reads exactly those, as
    ``ColumnRequest.from_problem`` encodes, and computing unused columns is
    the cost C06 item 3 removes.
    """

    def __init__(self, model, backend, request, enhancement_basis=None,
                 enhancement_scale=1.0):
        if tuple(int(r) for r in model.rows) != backend.rows:
            raise IdentityError("backend rows differ from the model rows (row order matters)")
        if ModelSettings.from_model(model) != backend.settings:
            raise IdentityError("backend model settings differ from the CPU model")
        self._model = model
        self._backend = backend
        self._request = request
        self._basis = enhancement_basis
        self._basis_digest = (None if enhancement_basis is None else
                              provenance.array_sha256(enhancement_basis))
        self._scale = float(enhancement_scale)
        self.calls = {"pressures": 0, "jacobian": 0}

    def __getattr__(self, name):
        # Only reached for attributes the facade does not define itself.
        # Private names are never forwarded (avoids recursion before __init__
        # has set _model, e.g. during copy or pickling).
        if name.startswith("_"):
            raise AttributeError(name)
        return getattr(self._model, name)

    def _check_basis(self, b, basis, scale):
        if b is None:
            return
        if basis is not self._basis:
            # A different object is accepted only if its bytes are identical.
            if basis is None or provenance.array_sha256(basis) != self._basis_digest:
                raise IdentityError("evaluate called with an enhancement basis that is not "
                                    "the one bound to the backend")
        if float(scale) != self._scale:
            raise IdentityError("enhancement scale %r differs from the bound %r" %
                                (scale, self._scale))

    def evaluate(self, a, T, kappa, b=None, enhancement_basis=None, enhancement_scale=1.0,
                 jac=False, basis_index=None, enhancement_index=None):
        self._check_basis(b, enhancement_basis, enhancement_scale)
        state = ParameterState(a, T, kappa, b)
        if not jac:
            self.calls["pressures"] += 1
            result = self._backend.evaluate(state, PRESSURES_ONLY)
            return result.pNG, result.pRG
        self.calls["jacobian"] += 1
        request = ColumnRequest(basis_index or (), enhancement_index or (),
                                temperature=self._request.temperature,
                                kappa=self._request.kappa)
        result = self._backend.evaluate(state, request)
        return result.pNG, result.pRG, result.jacobian_dict()


def install(problem, backend):
    """Route an existing ``fit_cops.Problem``'s evaluations through ``backend``.

    Returns the facade. The CPU model stays reachable as ``facade._model``,
    and ``problem.M = facade._model`` undoes the installation.
    """
    if backend.nbasis_enhancement != problem.nb:
        raise IdentityError("backend binds %d enhancement columns; problem has %d" %
                            (backend.nbasis_enhancement, problem.nb))
    # Every backend states the basis it holds (a CUDA backend: the uploaded
    # array's digest). A stale or reordered template set is refused here, before
    # any evaluation, not discovered later as a wrong pressure.
    bound = backend.identity().get("enhancement_basis_sha256")
    expected = array_digest(problem.enhancement_basis) if problem.nb else None
    if bound != expected:
        raise IdentityError("backend enhancement basis %s does not match the problem's %s "
                            "(stale or reordered templates)" % (bound, expected))
    scale = backend.identity().get("enhancement_scale")
    if problem.nb and scale is not None and float(scale) != float(problem.enhancement_scale):
        raise IdentityError("backend enhancement scale %r differs from the problem's %r" %
                            (scale, problem.enhancement_scale))
    facade = ModelFacade(problem.M, backend, ColumnRequest.from_problem(problem),
                         problem.enhancement_basis if problem.nb else None,
                         problem.enhancement_scale)
    problem.M = facade
    return facade
