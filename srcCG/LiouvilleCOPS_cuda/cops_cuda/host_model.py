"""Host-light stand-in for ``fit_cops.Model``: kernel headers only (ROADMAP C06/C09).

Why
---
``fit_cops.Model.__init__`` reads every kernel record of every table row into
host memory: about 15 arrays of 8 bytes per record, i.e. tens of GB for the
245M records of the 788k base. That is too much to sit beside the CPU campaign.
On the CUDA path the records live only on the device (``DeviceDataset``), and
``fit_cops.Problem`` needs only the model's per-point data (``rows``,
``useNG``, ``refNG``, ``T0``, ``Y``, ...) and its ``evaluate``. Through the
C02 facade, ``evaluate`` is routed to the GPU.

What
----
``header_model(...)`` returns an instance of a subclass of the reference
``fit_cops.Model``. The CPU code itself is unchanged and nothing here
monkeypatches it. Only ``load_kernels`` is overridden:

* each kernel file is opened for its header only (``formats.read_kernel_header``,
  with the same version, iPoint, count and size checks as the loader);
* ``VdotB`` and ``rh`` are the header's geometry words 1 and 8, exactly the
  values ``fit_cops.read_kernel`` returns (``{"VdotB": h[1], "rh": h[8]}``);
* ``indptr`` is the exact record CSR offsets; record arrays are empty;
* then the *reference* ``Model.prepare`` runs on those empty record arrays, so
  ``bm``, ``cNG`` and ``cRG`` come from the reference's own expressions and
  are not copied here.

Parsing in ``Model.__init__`` (header checks, use/ref columns, ``T0``, ``S``,
``SZAMode``, ...) runs unchanged, because ``__init__`` is inherited.

Every record-dependent method (``evaluate``, ``csr``, ``seen_faces``) refuses
with ``UnsupportedConfiguration``, and ``PhiCO2`` is ``None``. A header model is
therefore usable only with a GPU backend installed (``backend.install``). If
anything tried to evaluate it on the CPU, the error would be loud, never a
silent zero pressure.
"""

from pathlib import Path

import numpy as np

from . import formats, reference
from .errors import UnsupportedConfiguration

_CLASS = None


def _header_model_class():
    """The subclass of the reference ``Model``, created once per process."""
    global _CLASS
    if _CLASS is not None:
        return _CLASS
    fit_cops = reference.load("fit_cops")

    class HeaderModel(fit_cops.Model):
        """``fit_cops.Model`` whose kernels are represented by headers only."""

        records_loaded = False

        def load_kernels(self, rowsel):
            rowsel = np.asarray(rowsel, dtype=int)
            kernel_dir = Path(self.run) / "kernels"
            counts = np.empty(len(rowsel), dtype=np.int64)
            vb = np.empty(len(rowsel), dtype=np.float64)
            rh = np.empty(len(rowsel), dtype=np.float64)
            for k, irow in enumerate(rowsel):
                point = int(self.ipoint[irow])
                path = kernel_dir / formats.kernel_file_name(point)
                _, count, geometry = formats.read_kernel_header(path, point, self.nfaces)
                counts[k] = count
                vb[k], rh[k] = geometry[1], geometry[8]      # fit_cops.read_kernel: h[1], h[8]
            self.npt = len(rowsel)
            self.rows = rowsel
            self.indptr = np.concatenate([[0], np.cumsum(counts)]).astype(np.int64)
            self.VdotB, self.rh = vb, rh
            empty_f, empty_i = np.empty(0, dtype=np.float64), np.empty(0, dtype=np.int32)
            self.pt, self.face = empty_i, empty_i
            self.Kn, self.KPhi, self.KV, self.cos = empty_f, empty_f, empty_f, empty_f
            self.night = np.empty(0, dtype=bool)
            self.Q = {"H2O": empty_f, "CO2": empty_f}
            self.VB = empty_f
            print("kernel headers: %d points, %d element contributions (records stay on disk/device)"
                  % (self.npt, int(self.indptr[-1])))
            # The reference prepare(): bm, cNG, cRG from its own expressions;
            # the per-record factors are empty. PhiCO2 would be zeros here
            # (no records), so it is withdrawn rather than left misleading.
            self.prepare()
            self.PhiCO2 = None

        def _refuse(self, *args, **kwargs):
            raise UnsupportedConfiguration(
                "a header-only model holds no kernel records; install a GPU backend "
                "(backend.install) before evaluating")

        evaluate = _refuse
        csr = _refuse
        seen_faces = _refuse

    _CLASS = HeaderModel
    return _CLASS


def header_model(run, header, cols, rows, elements, coef):
    """``fit_cops.Model(run, header, cols, rows, elements, coef)`` without kernel records.

    Use like the reference: ``M = header_model(...); M.load_kernels(M.sel)``.
    """
    return _header_model_class()(str(run), header, cols, rows, elements, coef)


def is_header_model(model):
    return isinstance(model, _header_model_class())
