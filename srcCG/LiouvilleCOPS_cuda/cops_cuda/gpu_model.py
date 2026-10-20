"""GPU-backed ``fit_cops.Model`` for the unchanged fitter (ROADMAP C06 to C07, step 1).

Purpose
-------
``tools/gpu_fit.py`` runs the *unchanged* ``fit_cops.main()``: argument
parsing, cache audits, selection profiles, ``Problem``, SciPy TRF, every
output writer and the fit manifest. The only difference is that, inside that
process, the module attribute ``fit_cops.Model`` is the class returned by
``gpu_model_class()``. No file of ``../LiouvilleCOPS`` changes.

The class is a subclass of the header-only model (``host_model.py``), which
is itself a subclass of the reference ``fit_cops.Model``:

* ``__init__`` is the reference's, so the parsing and refusals are unchanged;
* ``load_kernels`` reads headers only, so no kernel record is ever held on
  the host;
* ``evaluate`` has the reference signature and return format (``pNG, pRG`` or
  ``pNG, pRG, J`` with ``J[("a", j)]``, ``J[("b", k)]``, ``J[("T", k)]``,
  ``J[("kappa", 0)]``) and runs on a ``CudaBackend``. With ``jac=True`` it
  returns the temperature and κ columns, as the reference always does;
* ``seen_faces`` runs on the device with the reference predicate
  (``kernels/masks.cu: seen_positive``).

Backends
--------
A backend is bound lazily on the first ``evaluate``/``seen_faces`` and
cached per (run, rows, device, mode). The fitter's two models (``MF`` over all
rows and ``M`` over the fit rows) therefore share one device dataset when
their rows are equal, as in the campaign profile. A backend holds at most one
enhancement basis:

* an evaluation with ``b is None`` ignores the basis, as the reference does,
  and uses any backend for those rows;
* an evaluation with ``b`` needs the bound basis to be that array (same object,
  or identical bytes by digest) with the same scale. Otherwise the backend is
  rebound once (dataset reload, about 70 s for the 788k base).

All binding and evaluation statistics are kept in ``REGISTRY.stats()`` for
the fit manifest's ``compute_backend`` section, written by the driver.
"""

import time

import numpy as np

from . import backend as bk
from . import host_model, masks
from .errors import IdentityError

_CLASS = None


class _Registry:
    """Process-wide backends, budgets and statistics."""

    def __init__(self):
        self.device_id = 0
        self.devices = [0]
        self.mode = "tree"
        self.records_per_block = 4096
        self.backends = {}            # (run, rows, device, mode) -> (backend, basis object)
        self.budgets = {}
        self.events = []
        self.calls = {"pressures": 0, "jacobian": 0, "seen_faces": 0}
        self.seconds = {"pressures": 0.0, "jacobian": 0.0, "bind": 0.0}

    def configure(self, device_id=0, mode="tree", records_per_block=4096, devices=None):
        """``devices`` with more than one entry selects ``multi_gpu.MultiGpuBackend``."""
        if self.backends:
            raise RuntimeError("configure before the first binding")
        self.devices = [int(d) for d in (devices if devices else [device_id])]
        self.device_id, self.mode = self.devices[0], str(mode)
        self.records_per_block = int(records_per_block)

    def budget(self):
        from .memory import DeviceBudget
        if self.device_id not in self.budgets:
            self.budgets[self.device_id] = DeviceBudget(self.device_id)
        return self.budgets[self.device_id]

    def backend_for(self, model, basis, scale):
        """A backend for ``model``'s rows; ``basis`` None means "any basis"."""
        from . import evaluator
        key = (str(model.run), tuple(int(r) for r in model.rows), tuple(self.devices), self.mode)
        entry = self.backends.get(key)
        if entry is not None:
            backend, bound = entry
            if basis is None or self._same_basis(backend, bound, basis, scale):
                return backend
            backend.release()
            del self.backends[key]
            self.events.append({"event": "rebind_for_new_basis", "rows": len(key[1])})
        t0 = time.perf_counter()
        scale = scale if basis is not None else 1.0
        if len(self.devices) > 1:
            from .multi_gpu import MultiGpuBackend
            backend = MultiGpuBackend(model._spec, model, basis, scale, devices=self.devices,
                                      mode=self.mode, records_per_block=self.records_per_block)
        else:
            backend = evaluator.CudaBackend(model, model.run, basis, scale,
                                            device_id=self.device_id, mode=self.mode,
                                            budget=self.budget(),
                                            records_per_block=self.records_per_block)
        seconds = time.perf_counter() - t0
        self.seconds["bind"] += seconds
        self.backends[key] = (backend, basis)
        self.events.append({"event": "bind", "rows": len(key[1]), "devices": list(self.devices),
                            "records": (int(backend.dataset.nrecords) if hasattr(backend, "dataset")
                                        else int(sum(backend.block_records))),
                            "enhancement_columns": backend.nbasis_enhancement,
                            "enhancement_basis_sha256": backend.enhancement_basis_sha256,
                            "seconds": seconds})
        print("gpu_model: bound %s in %.1f s" % (backend.summary(), seconds), flush=True)
        return backend

    @staticmethod
    def _same_basis(backend, bound, basis, scale):
        if backend.enhancement_basis is None:
            return False
        if float(scale) != backend.enhancement_scale:
            return False
        if basis is bound:
            return True
        return bk.array_digest(basis) == backend.enhancement_basis_sha256

    def stats(self):
        identities = [backend.identity() for backend, _ in self.backends.values()
                      if not backend.released]
        return {"device_id": self.device_id, "devices": list(self.devices), "mode": self.mode,
                "records_per_block": self.records_per_block,
                "calls": dict(self.calls), "seconds": dict(self.seconds),
                "events": list(self.events), "backends": identities}

    def close(self):
        for backend, _ in self.backends.values():
            backend.release()
        self.backends.clear()
        for budget in self.budgets.values():
            budget.close()
        self.budgets.clear()


REGISTRY = _Registry()


def gpu_model_class():
    """The GPU-backed subclass of the reference ``fit_cops.Model`` (created once)."""
    global _CLASS
    if _CLASS is not None:
        return _CLASS
    base = host_model._header_model_class()

    class GpuModel(base):
        """``fit_cops.Model`` whose kernels live on a GPU (see module docstring)."""

        def __init__(self, run, header, cols, rows, elements, coef):
            # Frozen copies of the constructor inputs: a multi-GPU worker rebuilds
            # its model from exactly these with the reference parser.
            import copy
            self._spec = {"run": str(run), "header": list(header), "cols": dict(cols),
                          "rows": copy.deepcopy(rows), "coef": np.array(coef, dtype=float)}
            super().__init__(run, header, cols, rows, elements, coef)

        def evaluate(self, a, T, kappa, b=None, enhancement_basis=None,
                     enhancement_scale=1.0, jac=False, basis_index=None,
                     enhancement_index=None):
            # The reference's own argument refusal (Model.evaluate).
            if b is not None and enhancement_basis is None:
                raise ValueError("enhancement coefficients need a facet basis")
            basis = enhancement_basis if b is not None else None
            backend = REGISTRY.backend_for(self, basis, enhancement_scale)
            if b is not None and len(np.asarray(b)) != backend.nbasis_enhancement:
                raise IdentityError("b has %d entries; the bound basis has %d columns" %
                                    (len(np.asarray(b)), backend.nbasis_enhancement))
            state = bk.ParameterState(a, T, kappa, b)
            t0 = time.perf_counter()
            if not jac:
                result = backend.evaluate(state, bk.PRESSURES_ONLY)
                REGISTRY.calls["pressures"] += 1
                REGISTRY.seconds["pressures"] += time.perf_counter() - t0
                return result.pNG, result.pRG
            request = bk.ColumnRequest(basis_index or (), enhancement_index or (),
                                       temperature=True, kappa=True)
            result = backend.evaluate(state, request)
            REGISTRY.calls["jacobian"] += 1
            REGISTRY.seconds["jacobian"] += time.perf_counter() - t0
            J = {key: (np.asarray(ng), np.asarray(rg)) for key, (ng, rg) in result.columns.items()}
            return result.pNG, result.pRG, J

        def seen_faces(self):
            """Facets contributing to a loaded point (reference predicate, on the device)."""
            backend = REGISTRY.backend_for(self, None, 1.0)
            REGISTRY.calls["seen_faces"] += 1
            if hasattr(backend, "seen_faces"):              # multi-GPU: union over workers
                return backend.seen_faces()
            return masks.seen_positive(backend.dataset)

    _CLASS = GpuModel
    return _CLASS
