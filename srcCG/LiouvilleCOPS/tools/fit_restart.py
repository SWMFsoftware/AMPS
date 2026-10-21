"""CPU implementation of the CUDA C12 restart wire format (schema version 1).

Ported from LiouvilleCOPS_cuda/cops_cuda/restart.py, SHA256
f44007885c48d2b84cf74a24c89185518cbc1a520071f33e83e3fa8592b1d926.
The sibling is read-only and is NOT a runtime dependency. Retain the original
schema string and array hashes so its reader/editor accepts these files.

This is a physical-state warm start, not a serialized SciPy trust-region state.
Geometry, template fields, fitted coefficients, temperatures, kappa, predictions,
selection, parameter bounds and provenance are carried in a pickle-free NPZ.
A warm start never replaces the new problem\'s prior or selection policy.
"""

import datetime as dt
import json
from pathlib import Path

import numpy as np

import hashlib
from cops_io import file_identity


def array_sha256(array):
    """C12 wire digest: dtype string, shape tuple, then contiguous C-order bytes."""
    array = np.ascontiguousarray(array)
    digest = hashlib.sha256()
    digest.update(array.dtype.str.encode("ascii"))
    digest.update(repr(tuple(array.shape)).encode("ascii"))
    digest.update(array.tobytes(order="C"))
    return digest.hexdigest()


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
    meta["array_sha256"] = {k: array_sha256(np.ascontiguousarray(v))
                            for k, v in restart.arrays.items() if v.dtype.kind in "fiub"}
    arrays = {k: np.asarray(v) for k, v in restart.arrays.items()}
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("xb") as stream:
        np.savez_compressed(stream, meta_json=np.array(json.dumps(meta, sort_keys=True)), **arrays)
    return file_identity(path, "fit_restart")


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
        if array_sha256(np.ascontiguousarray(arrays[k])) != digest:
            raise ValueError("restart array %s is corrupt (digest mismatch)" % k)
    r = Restart(meta, arrays)
    nf = len(r.area)
    if r.centroid.shape != (nf, 3) or r.activity_linear.shape != (nf,):
        raise ValueError("restart mesh arrays are inconsistent")
    if r.template_columns.shape != (nf, len(r.template_names)):
        raise ValueError("restart template columns do not match the template names")
    # Malformed physical state must fail before SciPy sees it. Infinite
    # optimizer bounds are legitimate; infinite coefficients/temperatures are not.
    expected = {"a": (25,), "coef_prior": (3, 25), "T": (3,), "kappa": (1,),
                "harmonic_values": (len(r.harmonic_modes),),
                "template_values": (len(r.template_names),)}
    for key, shape in expected.items():
        if arrays[key].shape != shape or not np.isfinite(arrays[key]).all():
            raise ValueError("invalid restart physical array: %s" % key)
    for key in ("centroid", "area", "activity_linear", "template_columns", "x", "x0"):
        if not np.isfinite(arrays[key]).all():
            raise ValueError("nonfinite restart array: %s" % key)
    if not np.isfinite(r.enhancement_scale) or r.enhancement_scale <= 0:
        raise ValueError("invalid restart enhancement scale")
    if len(set(r.template_names)) != len(r.template_names):
        raise ValueError("duplicate restart template names")
    import activity_templates as templates
    if templates.mesh_identity(r.centroid, r.area) != r.mesh_identity():
        raise ValueError("restart mesh identity does not match its geometry")
    return r


def restore_template_arguments(restart, explicit):
    """Mirror CUDA's default bundle restoration; explicit CLI templates win.

    Embedded columns describe the old state, but do not bypass current cliff
    provenance/mesh validation. The normal fitter reloads these source bundles.
    """
    if explicit or not restart.template_names:
        return explicit, False
    sources = restart.meta["templates"]["sources"]
    if [s["name"] for s in sources] != restart.template_names:
        raise ValueError("restart lacks ordered template sources; supply --activity-template")
    arguments = []
    for source in sources:
        bundle = source["bundle"]
        if file_identity(bundle["path"])["sha256"] != bundle["sha256"]:
            raise ValueError("restart template bundle changed: %s" % bundle["path"])
        arguments.append([bundle["path"], source["name"]])
    return arguments, list(restart.template_names)


def check_mesh(restart, elements):
    """Reject direct restarts on different geometry or facet ordering."""
    from activity_templates import mesh_identity
    if mesh_identity(elements["centroid"], elements["area"]) != restart.mesh_identity():
        raise ValueError("restart mesh differs from the fitting base; translate it first "
                         "with the CUDA restart_tool.py translate workflow")


# ---------------------------------------------------------------------- state maps
def harmonic_basis(centroid, modes):
    """The degree 5..L orthonormal basis columns for ``modes`` (activity_basis order)."""
    if not modes:
        return np.empty((len(centroid), 0))
    import fit_cops
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
        # Enhancement amplitudes are dimensionless relative to A_scale. If the
        # new run prior changes that scale, preserve the physical residual,
        # not the old optimizer coordinate. At identical scales use the exact
        # saved coordinate for bitwise same-problem continuation.
        if name in harmonics or name in templates:
            value *= restart.enhancement_scale / P.enhancement_scale
        clipped = min(max(value, P.lo[k]), P.hi[k])
        if clipped != value:
            report["clipped"].append({"name": name, "restart": value, "used": clipped})
        x[k] = clipped
        report["from_restart"].append(name)
    # Legacy coefficients the problem keeps fixed come from the run prior, not the restart.
    fixed = [j for j in range(25) if j not in P.fit_j]
    differing = [j for j in fixed if float(restart.a[j]) != float(P.coef[0][j])]
    report["fixed_legacy_coefficients_differing_from_restart"] = differing
    report["enhancement_scale_ratio"] = restart.enhancement_scale / P.enhancement_scale
    report["fixed_temperature"] = bool(P.args.fix_temperature)
    report["fixed_kappa"] = bool(P.args.fix_kappa)
    return x, report


# ------------------------------------------------------------------- from a fit
def _read_prediction_table(path):
    names, rows = None, []
    with open(path) as stream:
        for line in stream:
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


def from_fit(problem, result, fit_dir, manifest, elements, manifest_path=None):
    """Restart of a finished fit: ``problem`` (fit_cops.P), ``result`` (fit_one's dict)."""
    P = problem
    x = np.asarray(result["x"], dtype=float)
    a, b, T, kappa = P.unpack(x)
    nh = P.nh
    table = _read_prediction_table(Path(fit_dir) / "COPS_model.dat")
    import activity_templates as templates
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
        "fit_manifest": file_identity(manifest_path or (Path(fit_dir) / "fit_manifest.json")),
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
