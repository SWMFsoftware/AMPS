"""Fit restart files: one self-contained file per fitted (or edited) state (ROADMAP C12).

A restart file is a compressed NPZ, ``fit_restart.npz``, written by ``tools/gpu_fit.py`` at the
end of a fit and read by ``gpu_fit.py --restart-from`` and ``tools/restart_tool.py``.
``../LiouvilleCOPS`` is unchanged: the fitter's own products are written as before, and the
restart file is an additional product of the CUDA driver.

Contents
--------
The **physical state** (independent of the optimizer's scaling):

* ``a``: the 25 legacy H2O coefficients;
* ``harmonic_values``: the degree 5-8 residual amplitudes (by mode name);
* ``template_values``: the template amplitudes (by template name);
* ``T`` (Tmax, T75, Tmin) and ``kappa``.

The **data that make the state self-contained**:

* ``template_columns``: per-facet template fields, in the restart's template order;
* ``activity_linear``: the per-facet H2O linear activity ``Y a + A_scale (B b)`` before the clip;
* ``centroid``, ``area``: the mesh geometry, with its identity (``activity_templates.mesh_identity``);
* ``coef_prior``: the run's 3x25 activity coefficients;
* per observation row: ``point_ids``, ``utc``, ``altitude_m``, ``ref_ng``, ``ref_rg`` [Pa],
  ``use_ng``, ``use_rg``, ``saturated``, ``pred_ng``, ``pred_rg`` [Pa]. The predictions belong
  to this state.

The **metadata** (``meta_json``): schema, the base-run identities (``elements.bin``, table,
manifest), the selection profile, the fit settings (degrees, λ, bounds, fixes), the parameter
names, ``x``/``x0``/``lo``/``hi`` of the producing problem, the enhancement scale, the template
sources (bundle identities and definitions), the optimizer record, and the history of edits.

Mapping to a fit
----------------
``start_vector(problem, restart)`` builds an optimizer start for any ``fit_cops.Problem``
**by parameter name**:

* legacy ``aNN`` -> ``a[N] / ascale``;
* harmonics -> by mode name;
* ``c_<name>`` -> by template name;
* temperatures in either parameterization (free or ordered);
* ``kappa``.

Names absent from the restart keep the problem's start. Values are clipped to the bounds, and
every decision is reported. The problem's prior (``P.x0``) is never changed.
"""

import datetime as dt
import json
from pathlib import Path

import numpy as np

from . import provenance, reference

SCHEMA = "liouville_cops_cuda.fit_restart"
SCHEMA_VERSION = 1
ARRAY_KEYS = ("a", "coef_prior", "T", "kappa", "harmonic_values", "template_values", "template_columns",
              "activity_linear", "centroid", "area", "point_ids", "utc", "altitude_m", "ref_ng", "ref_rg",
              "use_ng", "use_rg", "saturated", "pred_ng", "pred_rg", "x", "x0", "lo", "hi")


class Restart:
    """In-memory restart state (see the module docstring)."""

    def __init__(self, meta, arrays):
        self.meta = meta
        self.arrays = arrays

    # convenient views
    def __getattr__(self, name):
        if name in ("meta", "arrays"):
            raise AttributeError(name)
        if name in self.arrays:
            return self.arrays[name]
        raise AttributeError(name)

    @property
    def template_names(self):
        return list(self.meta["templates"]["names"])

    @property
    def harmonic_modes(self):
        return [tuple(m) for m in self.meta["harmonics"]["modes"]]

    @property
    def enhancement_scale(self):
        return float(self.meta["enhancement_scale"])

    def mesh_identity(self):
        return self.meta["mesh"]["identity"]

    def copy(self):
        return Restart(json.loads(json.dumps(self.meta)),
                       {k: np.array(v, copy=True) for k, v in self.arrays.items()})


# ----------------------------------------------------------------------------- I/O
def write(path, restart):
    """Write ``restart`` to a new file (never replaces); returns its file identity."""
    path = Path(path)
    if path.exists():
        raise FileExistsError("restart file exists: %s" % path)
    meta = dict(restart.meta)
    meta["schema"], meta["schema_version"] = SCHEMA, SCHEMA_VERSION
    meta["array_sha256"] = {k: provenance.array_sha256(np.ascontiguousarray(v))
                            for k, v in restart.arrays.items() if v.dtype.kind in "fiub"}
    arrays = {k: np.asarray(v) for k, v in restart.arrays.items()}
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        np.savez_compressed(stream, meta_json=np.array(json.dumps(meta, sort_keys=True)), **arrays)
    return provenance.file_identity(path, "fit_restart")


def read(path):
    """Read and validate a restart file."""
    with np.load(str(path), allow_pickle=False) as z:
        meta = json.loads(str(z["meta_json"]))
        arrays = {k: np.array(z[k]) for k in z.files if k != "meta_json"}
    if meta.get("schema") != SCHEMA or meta.get("schema_version") != SCHEMA_VERSION:
        raise ValueError("%s is not a %s v%d restart file" % (path, SCHEMA, SCHEMA_VERSION))
    missing = [k for k in ARRAY_KEYS if k not in arrays]
    if missing:
        raise ValueError("restart file lacks %s" % missing)
    for k, digest in meta.get("array_sha256", {}).items():
        if provenance.array_sha256(np.ascontiguousarray(arrays[k])) != digest:
            raise ValueError("restart array %s is corrupt (digest mismatch)" % k)
    r = Restart(meta, arrays)
    nf = len(r.area)
    if r.centroid.shape != (nf, 3) or r.activity_linear.shape != (nf,):
        raise ValueError("restart mesh arrays are inconsistent")
    if r.template_columns.shape != (nf, len(r.template_names)):
        raise ValueError("restart template columns do not match the template names")
    templates = reference.load("activity_templates")
    if templates.mesh_identity(r.centroid, r.area) != r.mesh_identity():
        raise ValueError("restart mesh identity does not match its geometry")
    return r


# ---------------------------------------------------------------------- state maps
def harmonic_basis(centroid, modes):
    """The degree 5..L orthonormal basis columns for ``modes`` (activity_basis order)."""
    if not modes:
        return np.empty((len(centroid), 0))
    fit_cops = reference.load("fit_cops")
    lmax = max(l for l, _ in modes)
    full_modes = list(fit_cops.enhancement_modes(fit_cops.ENHANCEMENT_MIN_DEGREE, lmax))
    full = fit_cops.evaluate_enhancement_basis(centroid, fit_cops.ENHANCEMENT_MIN_DEGREE, lmax)
    index = [full_modes.index(tuple(m)) for m in modes]
    return np.ascontiguousarray(full[:, index])


def enhancement_basis(restart, elements_y=None):
    """Harmonic columns then template columns, as fit_cops.Problem orders them."""
    return np.ascontiguousarray(np.column_stack((harmonic_basis(restart.centroid, restart.harmonic_modes),
                                                 restart.template_columns)))


def b_vector(restart):
    return np.concatenate((np.asarray(restart.harmonic_values, float),
                           np.asarray(restart.template_values, float)))


def activity_linear(restart, Y):
    """``Y a + A_scale (B b)``: the H2O linear activity of the state (before the clip)."""
    act = Y @ np.asarray(restart.a, float)
    b = b_vector(restart)
    if b.size:
        act = act + restart.enhancement_scale * (enhancement_basis(restart) @ b)
    return act


def mode_name(mode):
    return "b_l%02d_m%+03d" % tuple(mode)


def start_vector(problem, restart):
    """Optimizer start for ``problem`` from ``restart`` (by name); returns (x, report)."""
    P = problem
    x = np.array(P.x0, dtype=float, copy=True)
    report = {"from_restart": [], "kept_problem_start": [], "clipped": []}
    # Exact optimizer coordinates when the restart came from the same parameterization
    # (same name and, for legacy coefficients, the same ascale): no x*ascale/ascale rounding.
    exact = {}
    names_r = restart.meta.get("names") or []
    if len(names_r) == len(restart.x):
        same_scale = float(restart.meta.get("ascale", np.nan)) == float(P.ascale)
        for k, name in enumerate(names_r):
            if not (name.startswith("a") and name[1:].isdigit()) or same_scale:
                exact[name] = float(restart.x[k])
    harmonics = {mode_name(m): float(v) for m, v in zip(restart.harmonic_modes, restart.harmonic_values)}
    templates = {"c_" + n: float(v) for n, v in zip(restart.template_names, restart.template_values)}
    Tmax, T75, Tmin = (float(v) for v in restart.T)
    temperature = {"Tmax": Tmax, "T75": T75, "Tmin": Tmin, "dT75": T75 - Tmin, "dTmax": Tmax - T75}
    for k, name in enumerate(P.names):
        value = None
        if name in exact:
            value = exact[name]
        elif name.startswith("a") and name[1:].isdigit():
            value = float(restart.a[int(name[1:])]) / P.ascale
        elif name in harmonics:
            value = harmonics[name]
        elif name in templates:
            value = templates[name]
        elif name in temperature:
            value = temperature[name]
        elif name == "kappa":
            value = float(np.asarray(restart.kappa).reshape(-1)[0])
        if value is None:
            report["kept_problem_start"].append(name)
            continue
        clipped = min(max(value, P.lo[k]), P.hi[k])
        if clipped != value:
            report["clipped"].append({"name": name, "restart": value, "used": clipped})
        x[k] = clipped
        report["from_restart"].append(name)
    # Legacy coefficients the problem keeps fixed come from the run prior, not the restart.
    fixed = [j for j in range(25) if j not in P.fit_j]
    differing = [j for j in fixed if float(restart.a[j]) != float(P.coef[0][j])]
    report["fixed_legacy_coefficients_differing_from_restart"] = differing
    return x, report


# ------------------------------------------------------------------- from a fit
def _read_prediction_table(path):
    names, rows = None, []
    for line in open(path):
        if line.startswith("#iPoint"):
            names = line[1:].split()
        elif not line.startswith("#") and line.strip():
            rows.append(line.split())
    col = {n: i for i, n in enumerate(names)}
    get = lambda n, f=float: np.array([f(r[col[n]]) for r in rows])
    return {"point_ids": get("iPoint", int), "utc": np.array([r[col["UTC"]] for r in rows]),
            "altitude_m": get("altitude"), "ref_ng": get("refNG"), "ref_rg": get("refRG"),
            "use_ng": get("useNG", int), "use_rg": get("useRG", int),
            "pred_ng": get("pNG"), "pred_rg": get("pRG")}


def from_fit(problem, result, fit_dir, manifest, elements):
    """Restart of a finished fit: ``problem`` (fit_cops.P), ``result`` (fit_one's dict)."""
    P = problem
    x = np.asarray(result["x"], dtype=float)
    a, b, T, kappa = P.unpack(x)
    nh = P.nh
    table = _read_prediction_table(Path(fit_dir) / "COPS_model.dat")
    templates = reference.load("activity_templates")
    centroid = np.asarray(elements["centroid"], dtype=float)
    area = np.asarray(elements["area"], dtype=float)
    columns = (np.asarray(P.enhancement_basis[:, nh:], dtype=float) if P.nt
               else np.empty((len(area), 0)))
    act = elements["Y"] @ a
    if P.nb:
        act = act + P.enhancement_scale * (P.enhancement_basis @ b)
    saturated = np.asarray(getattr(P, "saturated", np.zeros(len(table["point_ids"]), bool)), dtype=bool)
    if saturated.shape != table["point_ids"].shape:
        saturated = np.zeros(len(table["point_ids"]), bool)
    meta = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "kind": "fit",
        "mesh": {"identity": templates.mesh_identity(centroid, area), "faces": int(len(area))},
        "base": {k: manifest.get("inputs", {}).get(k) for k in ("base_elements", "base_model_table", "base_input")},
        "base_run_manifest": manifest.get("parent"),
        "selection": manifest.get("selection"),
        "settings": {"degree": int(P.args.degree), "enhancement_degree": int(P.enhancement_degree),
                     "lambda": float(result.get("lam", 0.0)), "fix_legacy_activity": bool(P.args.fix_legacy_activity),
                     "fix_temperature": bool(P.args.fix_temperature), "fix_kappa": bool(P.args.fix_kappa),
                     "ordered_temperature": bool(P.args.ordered_temperature),
                     "tbounds": [float(v) for v in P.args.tbounds]},
        "names": list(P.names), "ascale": float(P.ascale), "enhancement_scale": float(P.enhancement_scale),
        "harmonics": {"modes": [list(m) for m in P.enhancement_modes]},
        "templates": {"names": list(P.template_names), "sources": list(P.template_sources)},
        "optimizer": {k: (v if isinstance(v, (int, float, str, bool)) else str(v))
                      for k, v in result.items() if k in ("message", "status", "nfev", "optimality", "cost",
                                                          "success", "lam")},
        "fit_directory": str(Path(fit_dir).resolve()),
        "fit_manifest": provenance.file_identity(Path(fit_dir) / "fit_manifest.json"),
        "history": [{"utc": dt.datetime.now(dt.timezone.utc).isoformat(), "action": "fit",
                     "fit_directory": str(Path(fit_dir).resolve())}],
    }
    arrays = {"a": np.asarray(a, float), "coef_prior": np.asarray(P.coef, float), "T": np.asarray(T, float),
              "kappa": np.array([float(kappa)]), "harmonic_values": np.asarray(b[:nh], float),
              "template_values": np.asarray(b[nh:], float), "template_columns": columns,
              "activity_linear": np.asarray(act, float), "centroid": centroid, "area": area,
              "saturated": saturated, "x": x, "x0": np.asarray(P.x0, float),
              "lo": np.asarray(P.lo, float), "hi": np.asarray(P.hi, float), **table}
    return Restart(meta, arrays)


# --------------------------------------------------------------- GPU evaluation
TOOL_SCRATCH_PER_FACE = 32    # bytes per device facet of the restart/fine-tuning per-call arrays
RECORDS_PER_BLOCK = 4096      # the evaluator's production block size (tests use small blocks)


def record_factor_args(backend, records=None):
    """(kernel-name suffix, per-record factor arguments) of ``kernels/restart.cu`` for a backend.

    FP64 records: the evaluator's stored factors (wnH, wvH, wpH, xlit, lit). FP32 records (C14
    experiment, ``backend.mixed``): the records and per-point constants from which the ``_m``
    kernels recompute the same factors in FP64, as the mixed evaluator does. ``records``: one
    chunk of ``CudaBackend.record_chunks()`` (required when the backend streams records).
    """
    if records is None and getattr(backend.dataset, "streamed", False):
        from .errors import UnsupportedConfiguration
        raise UnsupportedConfiguration("internal: a streamed backend's records are passed chunk by chunk "
                                       "(CudaBackend.record_chunks)")
    if getattr(backend, "mixed", False):
        d = backend.dataset.arrays
        r = records or {n: d["kernel." + n] for n in ("kn", "kphi", "kv", "cos")}
        return "_m", (r["kn"], r["kphi"], r["kv"], r["cos"],
                      backend.arrays["pt"], np.int64(backend.dataset.npoints),
                      np.float64(backend.constants.COS75))
    a = backend.arrays
    return "", (a["wnH"], a["wvH"], a["wpH"], a["xlit"], a["lit"])


class StateEvaluator:
    """Pressures of restart states on a base run (one GPU backend per basis)."""

    def __init__(self, run, device_id=0, row_indices=None):
        """``row_indices``: optional contiguous block of table rows to load (default: all).

        Per-row results (pressures, point changes) then refer to these rows, in table order;
        ``point_ids`` are their ids. ``rows``/``cols`` stay the full table.
        """
        from . import host_model
        self.run = Path(run).resolve()
        self.device_id = int(device_id)
        fit_cops = reference.load("fit_cops")
        header, names, rows = fit_cops.read_model_file(str(self.run / "COPS_model.dat"))
        cols = {name: i for i, name in enumerate(names)}
        self.elements = fit_cops.read_elements(str(self.run / "elements.bin"))
        coef = fit_cops.read_activity_coefficients(str(reference.bind()), header)
        self.coef = coef
        self.header = header
        self.model = host_model.header_model(self.run, header, cols, rows, self.elements, coef)
        self.row_indices = (np.arange(len(rows)) if row_indices is None
                            else np.asarray(row_indices, dtype=int))
        self.model.load_kernels(self.row_indices)
        self.point_ids = np.asarray(self.model.ipoint, dtype=int)[self.row_indices]
        self.rows = rows
        self.cols = cols
        templates = reference.load("activity_templates")
        self.mesh_identity = templates.mesh_identity(self.elements["centroid"], self.elements["area"])
        self._backend, self._key = None, None

    def enhancement_scale(self):
        """A_scale of this base: sqrt(sum(area*(Y a0)^2)/sum(area)), as fit_cops.Problem."""
        legacy = self.elements["Y"] @ self.coef[0]
        n0 = float(np.sum(self.elements["area"] * legacy ** 2))
        return float(np.sqrt(n0 / np.sum(self.elements["area"])))

    def pressures(self, restart, basis=None, b=None, scale=None):
        """pNG, pRG (all table rows) of ``restart``'s state; optional basis/b/scale override."""
        from . import backend as bk, evaluator
        if restart.mesh_identity() != self.mesh_identity:
            raise ValueError("restart mesh differs from the base run's mesh (translate it first)")
        basis = enhancement_basis(restart) if basis is None else basis
        b = b_vector(restart) if b is None else np.asarray(b, float)
        scale = restart.enhancement_scale if scale is None else float(scale)
        has_b = b.size > 0
        key = (provenance.array_sha256(np.ascontiguousarray(basis)) if has_b else None, scale if has_b else None)
        if self._backend is None or key != self._key:
            if self._backend is not None:
                self._backend.release()
            self._backend = evaluator.CudaBackend(self.model, self.run, basis if has_b else None,
                                                  scale if has_b else 1.0, device_id=self.device_id,
                                                  records_per_block=RECORDS_PER_BLOCK,
                                                  scratch_per_face=TOOL_SCRATCH_PER_FACE)
            self._key = key
        state = bk.ParameterState(restart.a, restart.T, float(np.asarray(restart.kappa).reshape(-1)[0]),
                                  b if has_b else None)
        result = self._backend.evaluate(state)
        return np.asarray(result.pNG), np.asarray(result.pRG)

    def facet_point_shares(self, restart, delta_af, flagged_points):
        """Per facet: max over flagged points of |dpNG|/pNG and |dpRG|/|pRG| for a change dAF.

        ``delta_af`` is the per-facet change of the clipped H2O activity
        (``max(act', 0) - max(act, 0)``). The state's own pressures are the reference.
        Returns (share_ng, share_rg, pNG, pRG) with the shares as float64 per facet
        (``kernels/restart.cu``).
        """
        import cupy
        from . import cuda_sources
        from .memory import scratch_reservation
        pNG, pRG = self.pressures(restart)          # binds the backend with the restart's basis
        be = self._backend
        ds = be.dataset
        flags = np.isin(self.point_ids, np.asarray(list(flagged_points), dtype=int)).astype(np.uint8)
        c = be.constants
        Tmax, T75, Tmin = (float(v) for v in restart.T)
        nf = ds.nfaces_device
        with cupy.cuda.Device(self.device_id), scratch_reservation(
                be.budget, "facet point shares", {"dAF": (nf, np.float64), "max_ng": (nf, np.uint64),
                                                  "max_rg": (nf, np.uint64)}):
            d_daf = cupy.asarray(ds.gather_faces(np.asarray(delta_af, dtype=np.float64)))
            d_flag = cupy.asarray(flags)
            d_png, d_prg = cupy.asarray(pNG), cupy.asarray(pRG)
            d_crg = cupy.asarray(np.asarray(c.cRG["H2O"], dtype=np.float64))
            max_ng = cupy.zeros(nf, dtype=cupy.uint64)
            max_rg = cupy.zeros(nf, dtype=cupy.uint64)
            fn = None
            for g0, g1, table, records in be.record_chunks():     # atomicMax: any order
                suffix, factors = record_factor_args(be, records)
                fn = fn or cuda_sources.function("restart", "facet_point_shares" + suffix, self.device_id)
                fn((g1 - g0,), (256,),
                   (table, records["face"]) + factors + (
                    d_daf, d_flag, d_png, d_prg, d_crg, np.float64(c.cNG["H2O"] * c.bm["H2O"]),
                    np.float64(c.bm["H2O"]), np.float64(Tmax), np.float64(T75), np.float64(Tmin),
                    np.float64(c.KT), max_ng, max_rg))
            share_ng = ds.scatter_faces(max_ng.get().view(np.float64), 0.0)
            share_rg = ds.scatter_faces(max_rg.get().view(np.float64), 0.0)
            del d_daf, d_flag, d_png, d_prg, d_crg, max_ng, max_rg
        return share_ng, share_rg, pNG, pRG

    def point_changes(self, restart, delta_af):
        """Total change of every point's (pNG, pRG) for a per-facet change of the clipped activity.

        Exact: the pressures are linear in each facet's clipped H2O activity
        (``kernels/restart.cu: point_change_partials``). Requires the backend bound for
        ``restart`` (call ``pressures`` or ``facet_point_shares`` first).
        """
        import cupy
        from . import cuda_sources
        from .memory import scratch_reservation
        be = self._backend
        ds, c = be.dataset, be.constants
        Tmax, T75, Tmin = (float(v) for v in restart.T)
        nblocks = int(be._first_block[-1])
        with cupy.cuda.Device(self.device_id), scratch_reservation(
                be.budget, "point changes", {"dAF": (ds.nfaces_device, np.float64),
                                             "partial": (2 * nblocks, np.float64)}):
            d_daf = cupy.asarray(ds.gather_faces(np.asarray(delta_af, dtype=np.float64)))
            d_crg = cupy.asarray(np.asarray(c.cRG["H2O"], dtype=np.float64))
            partial = cupy.empty(2 * nblocks, dtype=np.float64)
            fn = None
            for g0, g1, table, records in be.record_chunks():     # partials by global block
                suffix, factors = record_factor_args(be, records)
                fn = fn or cuda_sources.function("restart", "point_change_partials" + suffix, self.device_id)
                fn((g1 - g0,), (256,),
                   (table, records["face"]) + factors + (
                    d_daf, d_crg, np.float64(c.cNG["H2O"] * c.bm["H2O"]), np.float64(c.bm["H2O"]),
                    np.float64(Tmax), np.float64(T75), np.float64(Tmin), np.float64(c.KT),
                    partial[2 * g0:2 * g1]))
            host = partial.get().reshape(nblocks, 2)
            del d_daf, d_crg, partial
        first = np.asarray(be._first_block[:-1], dtype=np.int64)
        return np.add.reduceat(host[:, 0], first), np.add.reduceat(host[:, 1], first)

    def table_arrays(self):
        """Observation rows of the base table (for restarts translated to this base)."""
        c, rows = self.cols, self.rows
        # The base table stores the reference pressures already converted to Pa (main.cpp: RefToPa).
        get = lambda n, f=float: np.array([f(r[c[n]]) for r in rows])
        return {"point_ids": get("iPoint", int), "utc": np.array([r[c["UTC"]] for r in rows]),
                "altitude_m": get("altitude"), "ref_ng": get("refNG"), "ref_rg": get("refRG"),
                "use_ng": get("useNG", int), "use_rg": get("useRG", int)}

    def release(self):
        if self._backend is not None:
            self._backend.release()
            self._backend = None
