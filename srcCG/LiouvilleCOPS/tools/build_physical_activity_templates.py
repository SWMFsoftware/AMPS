#!/usr/bin/env python3
"""Build physically declared localized H2O activity templates.

Two different hypotheses are kept separate:

``late_ng_not_rg``
    A binary geometry template.  A facet is one when at least one static
    kernel in the explicitly supplied late interval has positive NG response
    while that same observation geometry has zero RG aperture response.
    Exact kernel zeros define "not seen"; no fitted pressure or arbitrary
    sensitivity-ratio threshold enters this mask.  Integrated per-unit-
    activity pressure sensitivities are diagnostics, not region thresholds.

``cliff_ge_<angle>deg``
    A binary hypothesis template based on effective gravitational slope.  The
    acceleration is homogeneous-polyhedron self-gravity plus centrifugal
    acceleration from the same SPICE body frame/rate used by LiouvilleCOPS.
    The activity amplitude is fitted later; this tool does not assert or
    invent a production-versus-slope law.

The output bundle is tied to the exact ``elements.bin`` facet order.  It can be
loaded by ``fit_cops.py --activity-template BUNDLE NAME`` only with that mesh.
The external ``polyhedral_gravity`` package is ESA/TU Munich's validated
implementation of the Tsoulis homogeneous-polyhedron formulation.
"""

import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import re
import subprocess
import sys
import time

import numpy as np

from activity_templates import write_bundle
from cops_io import file_identity, validate_fit_base
from cops_paths import output_path, product_input
from cops_selection import parse_utc
from fit_cops import (AMU, KB, MASS, SQRTPI, COS75, KT, header_value,
                      read_elements, read_kernel, sensitivities)
from gravity_cache import load_cache, make_definition, write_cache
from gravity_topology import require_outward_closed_body, ORIENTATION_POLICY


ROOT = Path(__file__).resolve().parents[1]


def nastran_real(text):
    """Parse standard and NASTRAN's exponent-without-E real syntax."""
    value = text.strip().replace("D", "E").replace("d", "E")
    if not value:
        raise ValueError("empty NASTRAN real field")
    for index in range(1, len(value)):
        if value[index] in "+-" and value[index - 1] not in "Ee":
            value = value[:index] + "E" + value[index:]
            break
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("nonfinite NASTRAN real")
    return result


def read_nastran_surface(path):
    """Read GRID/GRID* and CTRIA3 cards with the solver's field conventions."""
    nodes, triangles = {}, []
    with Path(path).open("r", errors="strict") as stream:
        lines = iter(stream)
        for line_number, line in enumerate(lines, 1):
            if line.startswith("GRID*"):
                try:
                    continuation = next(lines)
                except StopIteration:
                    raise ValueError("missing GRID* continuation at end of %s" % path)
                node = int(line[8:24])
                xyz = (nastran_real(line[40:56]), nastran_real(line[56:72]),
                       nastran_real(continuation[8:24]))
                if node in nodes:
                    raise ValueError("duplicate GRID %d in %s" % (node, path))
                nodes[node] = xyz
            elif line.startswith("GRID") and not line.startswith("GRID*"):
                node = int(line[8:16])
                xyz = (nastran_real(line[24:32]), nastran_real(line[32:40]),
                       nastran_real(line[40:48]))
                if node in nodes:
                    raise ValueError("duplicate GRID %d in %s" % (node, path))
                nodes[node] = xyz
            elif line.startswith("CTRIA3"):
                fields = line.split()
                if len(fields) < 6:
                    # Eight-digit IDs occupy the whole NASTRAN small field:
                    # two adjacent IDs need not have intervening whitespace.
                    # Match the C++ reader's fallback without changing winding.
                    if len(line.rstrip("\r\n")) < 48:
                        raise ValueError("invalid CTRIA3 at %s:%d" % (path, line_number))
                    values = [int(line[start:start + 8]) for start in range(8, 48, 8)]
                    triangles.append(tuple(values[2:5]))
                else:
                    triangles.append(tuple(int(value) for value in fields[3:6]))
    if not triangles:
        raise ValueError("gravity model has no CTRIA3 surface: %s" % path)
    used = sorted({node for triangle in triangles for node in triangle})
    missing = [node for node in used if node not in nodes]
    if missing:
        raise ValueError("CTRIA3 references undefined GRID %d" % missing[0])
    index = {node: local for local, node in enumerate(used)}
    vertices = np.asarray([nodes[node] for node in used], dtype=float)
    faces = np.asarray([[index[node] for node in triangle]
                        for triangle in triangles], dtype=np.int64)
    return vertices, faces


def utc_midpoint(start, stop):
    epoch = dt.datetime(2016, 1, 1)
    return (epoch + dt.timedelta(seconds=0.5 * (parse_utc(start) + parse_utc(stop)))).strftime(
        "%Y-%m-%dT%H:%M:%S.%f").rstrip("0").rstrip(".")


def spice_omega(executable, input_file, utc_values):
    """Evaluate the exact forward-model SPICE angular velocity and summarize it."""
    command = [str(executable), str(input_file)] + list(utc_values)
    completed = subprocess.run(command, cwd=str(ROOT), text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if completed.returncode:
        raise ValueError("SPICE angular-velocity helper failed:\n" + completed.stdout)
    records = []
    for line in completed.stdout.splitlines():
        if line.startswith("OMEGA "):
            fields = line.split()
            records.append([float(value) for value in fields[2:5]])
    if len(records) != len(utc_values):
        raise ValueError("SPICE helper returned %d OMEGA rows for %d epochs" %
                         (len(records), len(utc_values)))
    values = np.asarray(records)
    median = np.median(values, axis=0)
    deviation = np.linalg.norm(values - median, axis=1)
    return median, float(np.max(deviation)), command


def point_columns(names):
    return {name: index for index, name in enumerate(names)}


def analytical_temperature(cos_sza, shadow, temperatures):
    """Match fit_cops.py's analytical gas-temperature law per contribution."""
    tmax, t75, tmin = temperatures
    c = np.where((cos_sza <= 0.0) | shadow, -1.0,
                 np.clip(cos_sza, -1.0, 1.0))
    lit = c >= COS75
    xlit = np.where(lit, 1.0 / np.where(lit, c, 1.0) - 1.0, 0.0)
    illuminated = tmax + KT * (tmax - t75) * xlit
    return np.where(lit & (illuminated > tmin), illuminated, tmin)


def late_gauge_template(run, header, names, rows, nfaces, start, stop):
    """Return the exact visibility mask and diagnostic integrated sensitivities."""
    columns = point_columns(names)
    t0, t1 = parse_utc(start), parse_utc(stop)
    selected = [row for row in rows
                if t0 <= parse_utc(row[columns["UTC"]]) <= t1 and
                not int(row[columns["inside"]])]
    if not selected:
        raise ValueError("late interval selects no outside-nucleus observation rows")

    model = header_value(header, "HeliocentricModel")
    exponent = 4.2143229 if model == "shap5" else 4.8
    temperatures = np.array([float(header_value(header, "Analytical" + name))
                             for name in ("Tmax", "T75", "Tmin")])
    tng = float(header_value(header, "TemperatureNG"))
    sensitivity_h2o = sensitivities(header)[0]
    rg_model = header_value(header, "RGPressureModel")
    sza_mode = header_value(header, "SZAMode")
    if sza_mode not in ("continuous", "srccg_table"):
        raise ValueError("unsupported SZAMode %s" % sza_mode)
    mass = MASS["H2O"]
    bm = math.sqrt(mass / (2.0 * KB)) * SQRTPI / 2.0
    cng = KB * tng / sensitivity_h2o
    ng_seen = np.zeros(nfaces, dtype=bool)
    rg_seen = np.zeros(nfaces, dtype=bool)
    ng_not_rg = np.zeros(nfaces, dtype=bool)
    ng_response = np.zeros(nfaces)
    rg_response = np.zeros(nfaces)

    for row in selected:
        point = int(row[columns["iPoint"]])
        saved_point, kernel_header, records = read_kernel(
            os.path.join(run, "kernels", "point=%04d.bin" % point))
        if saved_point != point:
            raise ValueError("kernel point mismatch for %d" % point)
        face = records["face"]
        kn = records["Kn"]
        kphi = records["KPhi"]
        kv = records["KV"]
        point_ng = kn > 0.0
        point_rg = (kphi != 0.0) | (kv != 0.0)
        ng_seen[face[point_ng]] = True
        rg_seen[face[point_rg]] = True
        # The union is evaluated point by point: the same surface facet can be
        # NG-only at one geometry and enter the RG aperture after the nucleus
        # rotates.  Requiring it never to enter RG during the whole interval
        # would answer a different, overly restrictive question.
        ng_not_rg[face[point_ng & ~point_rg]] = True

        rh_au = kernel_header["rh"] / 149597870700.0
        qmin, qmax = 0.2 * 600.0 / rh_au ** exponent, 600.0 / rh_au ** exponent
        cosine = np.clip(records["cos"], -1.0, 1.0)
        night = (cosine <= 0.0) | (records["shadow"] != 0)
        if sza_mode == "continuous":
            production_cosine = cosine
        else:
            # Match Model.__init__ and cSourceModel exactly: the historical
            # srcCG table mode truncates illuminated SZA to an integer degree
            # and caps it at 89 degrees.  This changes only the saved response
            # diagnostic, never the exact-zero visibility mask above.
            angle = np.minimum(np.floor(np.degrees(np.arccos(cosine))), 89.0)
            production_cosine = np.cos(np.radians(angle))
        q = np.where(night, qmin,
                     qmin + (qmax - qmin) * production_cosine)
        temperature = analytical_temperature(cosine, records["shadow"] != 0,
                                             temperatures)
        inv_sqrt_t = 1.0 / np.sqrt(temperature)
        trg = float(row[columns["T_RG"]])
        if rg_model == "equilibrium_sphere":
            crg = math.sqrt(2.0 * math.pi * mass * KB * trg) / sensitivity_h2o
        elif rg_model == "calibration_srccg":
            crg = math.sqrt(math.pi * mass * KB * trg / 2.0) / sensitivity_h2o
        elif rg_model == "dsmc_srccg":
            ktrg = KB * trg
            crg = ktrg * math.pi * math.sqrt(mass / (2.0 * ktrg)) / 2.0 / sensitivity_h2o
        else:
            raise ValueError("unsupported RGPressureModel %s" % rg_model)
        dng = cng * bm * q * kn * inv_sqrt_t
        drg = crg * q * (kphi + bm * kernel_header["VdotB"] * kv * inv_sqrt_t)
        np.add.at(ng_response, face, np.maximum(dng, 0.0))
        np.add.at(rg_response, face, np.maximum(drg, 0.0))

    interval_exclusive = ng_seen & ~rg_seen
    return ng_not_rg.astype(float), ng_response, rg_response, {
        "interval_utc": [start, stop],
        "selected_point_count": len(selected),
        "selected_point_ids": [int(row[columns["iPoint"]]) for row in selected],
        "definition": "union_over_points_of_Kn_gt_0_and_same_point_KPhi_eq_0_and_KV_eq_0",
        "ng_seen_face_count": int(np.count_nonzero(ng_seen)),
        "rg_seen_face_count": int(np.count_nonzero(rg_seen)),
        "ng_not_rg_union_face_count": int(np.count_nonzero(ng_not_rg)),
        "never_rg_over_interval_face_count": int(np.count_nonzero(interval_exclusive)),
    }


def effective_slope(centroid, outward_normal, gravity_vertices, gravity_faces,
                    density, omega, offset, offset_tolerance, chunk, parallel):
    """Evaluate self-gravity plus centrifugal acceleration and surface slope."""
    try:
        import polyhedral_gravity as pg
    except ImportError:
        raise ValueError("polyhedral_gravity is required; install the pinned version in "
                         "requirements-gravity.txt") from None
    # HEAL corrupts metre-unit kilometre-scale meshes (HEAL-corruption.md).
    # Validate topology independently and preserve the supplied face order.
    topology = require_outward_closed_body(gravity_vertices, gravity_faces)
    polyhedron = pg.Polyhedron(
        polyhedral_source=(gravity_vertices, gravity_faces), density=density,
        normal_orientation=pg.NormalOrientation.OUTWARDS,
        integrity_check=pg.PolyhedronIntegrity.DISABLE,
        metric_unit=pg.MetricUnit.METER)
    evaluator = pg.GravityEvaluable(polyhedron=polyhedron)
    centrifugal = -np.cross(omega, np.cross(omega, centroid))

    def evaluate_acceleration(points):
        result = evaluator(computation_points=points.tolist(), parallel=parallel)
        return np.asarray([record[1] for record in result])

    def slopes(indices, gravity_values):
        effective_values = gravity_values + centrifugal[indices]
        magnitude_values = np.linalg.norm(effective_values, axis=1)
        cosine_values = (np.sum(-outward_normal[indices] * effective_values, axis=1) /
                         magnitude_values)
        return np.degrees(np.arccos(np.clip(cosine_values, -1.0, 1.0)))

    # Surface-limit formulas can lose precision exactly on an edge/face.  The
    # declared outward offset avoids the singular representation; this audit
    # proves on a deterministic surface-wide sample that halving the offset
    # does not materially change the derived slope classification.
    sample = np.unique(np.linspace(0, len(centroid) - 1,
                                   min(512, len(centroid)), dtype=int))
    gravity_offset = evaluate_acceleration(
        centroid[sample] + offset * outward_normal[sample])
    gravity_half = evaluate_acceleration(
        centroid[sample] + 0.5 * offset * outward_normal[sample])
    offset_difference = np.abs(slopes(sample, gravity_offset) -
                               slopes(sample, gravity_half))
    offset_max = float(np.max(offset_difference))
    if not np.isfinite(offset_max) or offset_max > offset_tolerance:
        raise ValueError(
            "gravity slope changes by %.6g deg when the evaluation offset is halved; "
            "this exceeds --gravity-offset-tolerance %.6g deg" %
            (offset_max, offset_tolerance))

    points = centroid + offset * outward_normal
    gravity = np.empty_like(points)
    started = time.monotonic()
    for begin in range(0, len(points), chunk):
        end = min(begin + chunk, len(points))
        gravity[begin:end] = evaluate_acceleration(points[begin:end])
        elapsed = time.monotonic() - started
        print("gravity: %d/%d facets (%.1f%%, %.1f s)" %
              (end, len(points), 100.0 * end / len(points), elapsed), flush=True)
    effective = gravity + centrifugal
    magnitude = np.linalg.norm(effective, axis=1)
    if not np.all(np.isfinite(effective)) or np.any(magnitude <= 0.0):
        raise ValueError("nonfinite/zero effective gravity at a target facet")
    inward = -outward_normal
    cosine = np.sum(inward * effective, axis=1) / magnitude
    slope = np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0)))
    return slope, gravity, centrifugal, {
        "integrity_check": "DISABLE after topology audit",
        "orientation_policy": ORIENTATION_POLICY,
        "topology": topology,
        "package_version": pg.__version__,
        "package_binary": file_identity(Path(pg.__file__)),
        "package_commit": getattr(pg, "__commit__", None),
        "parallelization": getattr(pg, "__parallelization__", None),
        "formulation": "Tsoulis homogeneous polyhedron via ESA polyhedral-gravity",
        "offset_halving_sample_count": int(len(sample)),
        "offset_halving_slope_difference_deg_max": offset_max,
        "offset_halving_slope_difference_deg_p95": float(np.percentile(
            offset_difference, 95.0)),
        "offset_tolerance_deg": offset_tolerance,
        "evaluation_seconds": time.monotonic() - started,
    }


def validate_gravity_resolution(centroid, outward_normal, omega, offset,
                                density, reference_slope, vertices, faces,
                                tolerance, parallel):
    """Compare a second gravity polyhedron on a deterministic target sample."""
    try:
        import polyhedral_gravity as pg
    except ImportError:
        raise ValueError("polyhedral_gravity is required") from None
    topology = require_outward_closed_body(vertices, faces, "validation gravity mesh")
    polyhedron = pg.Polyhedron(
        polyhedral_source=(vertices, faces), density=density,
        normal_orientation=pg.NormalOrientation.OUTWARDS,
        integrity_check=pg.PolyhedronIntegrity.DISABLE,
        metric_unit=pg.MetricUnit.METER)
    evaluator = pg.GravityEvaluable(polyhedron=polyhedron)
    sample = np.unique(np.linspace(0, len(centroid) - 1,
                                   min(512, len(centroid)), dtype=int))
    points = centroid[sample] + offset * outward_normal[sample]
    result = evaluator(computation_points=points.tolist(), parallel=parallel)
    gravity = np.asarray([record[1] for record in result])
    centrifugal = -np.cross(omega, np.cross(omega, centroid[sample]))
    effective = gravity + centrifugal
    magnitude = np.linalg.norm(effective, axis=1)
    if not np.all(np.isfinite(effective)) or np.any(magnitude <= 0):
        raise ValueError("nonfinite/zero effective gravity in resolution audit")
    slope = np.degrees(np.arccos(np.clip(
        np.sum(-outward_normal[sample] * effective, axis=1) / magnitude,
        -1.0, 1.0)))
    difference = np.abs(slope - reference_slope[sample])
    p95 = float(np.percentile(difference, 95.0))
    if p95 > tolerance:
        raise ValueError(
            "gravity-model resolution changes the sampled slope by %.6g deg at p95; "
            "this exceeds --gravity-resolution-tolerance %.6g deg" %
            (p95, tolerance))
    return {
        "integrity_check": "DISABLE after topology audit",
        "orientation_policy": ORIENTATION_POLICY,
        "topology": topology,
        "sample_count": int(len(sample)),
        "slope_difference_deg_median": float(np.median(difference)),
        "slope_difference_deg_p95": p95,
        "slope_difference_deg_max": float(np.max(difference)),
        "p95_tolerance_deg": tolerance,
    }


def area_fraction(mask, area):
    return float(np.sum(area[np.asarray(mask, dtype=bool)]) / np.sum(area))


def write_map(path, centroid, area, fields):
    """Write compact area-weighted longitude/latitude diagnostic maps."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    radius = np.linalg.norm(centroid, axis=1)
    longitude = np.degrees(np.arctan2(centroid[:, 1], centroid[:, 0]))
    latitude = np.degrees(np.arcsin(np.clip(centroid[:, 2] / radius, -1.0, 1.0)))
    lon_edges = np.linspace(-180.0, 180.0, 241)
    lat_edges = np.linspace(-90.0, 90.0, 121)
    area_grid = np.histogram2d(latitude, longitude, bins=(lat_edges, lon_edges),
                               weights=area)[0]
    fig, axes = plt.subplots(len(fields), 1, figsize=(14, 3.4 * len(fields)),
                             squeeze=False, sharex=True, sharey=True)
    for axis, (name, values, cmap, limits) in zip(axes[:, 0], fields):
        weighted = np.histogram2d(latitude, longitude,
                                  bins=(lat_edges, lon_edges),
                                  weights=area * values)[0]
        grid = np.divide(weighted, area_grid, out=np.full_like(weighted, np.nan),
                         where=area_grid > 0.0)
        image = axis.pcolormesh(lon_edges, lat_edges, grid, shading="auto",
                                cmap=cmap, vmin=limits[0], vmax=limits[1])
        fig.colorbar(image, ax=axis, label=name)
        axis.set_ylabel("body-frame latitude [deg]")
        axis.set_title(name)
    axes[-1, 0].set_xlabel("body-frame longitude [deg]")
    fig.tight_layout()
    fig.savefig(path, dpi=140)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run", help="manifested static LiouvilleCOPS fitting base")
    parser.add_argument("--out", required=True, help="new .npz bundle under COPS_OUTPUT_ROOT or absolute path")
    parser.add_argument("--late-start", required=True, help="UTC start of the declared late NG interval")
    parser.add_argument("--late-stop", required=True, help="UTC stop of the declared late NG interval")
    parser.add_argument("--gravity-shape", required=True,
                        help="closed NASTRAN CTRIA3 surface defining the homogeneous gravity body")
    parser.add_argument("--gravity-cache",
                        help="reusable gravity .npz under COPS_OUTPUT_ROOT or absolute path "
                             "(default: <out stem>.gravity.npz)")
    parser.add_argument("--gravity-validation-shape",
                        help="optional second closed surface for a sampled gravity-slope resolution check")
    parser.add_argument("--density", required=True, type=float,
                        help="homogeneous bulk density [kg/m^3]")
    parser.add_argument("--cliff-threshold", type=float, action="append", default=None,
                        help="effective-slope threshold [deg]; repeatable (default: 60)")
    parser.add_argument("--gravity-offset", type=float, default=0.1,
                        help="outward evaluation offset from each facet centroid [m] (default: 0.1)")
    parser.add_argument("--gravity-offset-tolerance", type=float, default=0.25,
                        help="maximum sampled slope change when offset is halved [deg] (default: 0.25)")
    parser.add_argument("--gravity-resolution-tolerance", type=float, default=5.0,
                        help="maximum p95 slope change for --gravity-validation-shape [deg] (default: 5)")
    parser.add_argument("--gravity-chunk", type=int, default=2048,
                        help="polyhedral-gravity evaluation batch size (default: 2048)")
    parser.add_argument("--serial-gravity", action="store_true",
                        help="disable the package's parallel evaluator")
    parser.add_argument("--omega-executable", default="spice_angular_velocity",
                        help="helper built by 'make spice_angular_velocity'")
    args = parser.parse_args()
    if parse_utc(args.late_stop) < parse_utc(args.late_start):
        parser.error("--late-stop must not precede --late-start")
    if not math.isfinite(args.density) or args.density <= 0.0:
        parser.error("--density must be finite and positive")
    if not math.isfinite(args.gravity_offset) or args.gravity_offset <= 0.0:
        parser.error("--gravity-offset must be finite and positive")
    if (not math.isfinite(args.gravity_offset_tolerance) or
            args.gravity_offset_tolerance <= 0.0):
        parser.error("--gravity-offset-tolerance must be finite and positive")
    if (not math.isfinite(args.gravity_resolution_tolerance) or
            args.gravity_resolution_tolerance <= 0.0):
        parser.error("--gravity-resolution-tolerance must be finite and positive")
    if args.gravity_chunk <= 0:
        parser.error("--gravity-chunk must be positive")
    thresholds = sorted(set(args.cliff_threshold or [60.0]))
    if any(not math.isfinite(value) or not 0.0 < value < 90.0
           for value in thresholds):
        parser.error("every --cliff-threshold must be finite and in (0,90)")

    # Resolve and preflight every final template product before the expensive
    # kernel/gravity work.  An existing gravity cache is intentionally allowed:
    # it is loaded only after its physical inputs and exact facet ordering pass
    # the checks in gravity_cache.py.  Template products are never overwritten.
    destination = output_path(args.out)
    if destination.suffix != ".npz":
        raise ValueError("--out must end in .npz")
    json_path = destination.with_suffix(".json")
    map_path = destination.with_suffix(".png")
    for product in (destination, json_path, map_path):
        if product.exists():
            raise ValueError("template product already exists: %s" % product)
    if args.gravity_cache:
        gravity_cache_path = output_path(args.gravity_cache)
    else:
        gravity_cache_path = destination.with_suffix(".gravity.npz")
    if gravity_cache_path.suffix != ".npz":
        raise ValueError("--gravity-cache must end in .npz")
    if gravity_cache_path in (destination, json_path, map_path):
        raise ValueError("--gravity-cache must differ from the template products")

    run = Path(product_input(args.run)).resolve()
    audit = validate_fit_base(run)
    if audit["legacy_manifest"]:
        raise ValueError("template construction requires a manifested fitting base")
    header, names, rows = (audit[key] for key in ("headers", "names", "rows"))
    elements = read_elements(run / "elements.bin")
    nfaces = len(elements["area"])
    late, ng_response, rg_response, late_meta = late_gauge_template(
        str(run), header, names, rows, nfaces, args.late_start, args.late_stop)
    print("late NG-not-RG: %d facets, %.6f of surface area" %
          (np.count_nonzero(late), area_fraction(late, elements["area"])), flush=True)

    gravity_shape = Path(args.gravity_shape).expanduser().resolve()
    gravity_vertices, gravity_faces = read_nastran_surface(gravity_shape)
    # Also audit on cache reuse: geometry provenance must remain inspectable
    # even when no package evaluator is constructed in this invocation.
    gravity_topology = require_outward_closed_body(gravity_vertices, gravity_faces)
    gravity_shape_identity = file_identity(gravity_shape)
    print("gravity surface: %d vertices, %d triangles" %
          (len(gravity_vertices), len(gravity_faces)), flush=True)
    omega_executable = Path(args.omega_executable)
    if not omega_executable.is_absolute():
        omega_executable = (ROOT / omega_executable).resolve()
    input_snapshot = run / "input.snapshot.in"
    omega_epochs = [args.late_start, utc_midpoint(args.late_start, args.late_stop),
                    args.late_stop]
    omega, omega_variation, omega_command = spice_omega(
        omega_executable, input_snapshot, omega_epochs)
    print("SPICE omega [rad/s]: %s; max sampled deviation %.3e" %
          (" ".join("%.17e" % value for value in omega), omega_variation), flush=True)

    gravity_definition = make_definition(
        elements["centroid"], elements["normal"], elements["area"],
        gravity_shape_identity, len(gravity_vertices), len(gravity_faces),
        args.density, args.gravity_offset)
    if gravity_cache_path.exists():
        cached = load_cache(
            gravity_cache_path, elements["centroid"], elements["normal"],
            elements["area"], gravity_definition, omega)
        gravity = cached["self_gravity"]
        centrifugal = cached["centrifugal"]
        effective = cached["effective"]
        slope = cached["slope"]
        gravity_meta = cached["metadata"].get("calculation", {})
        if (gravity_meta.get("orientation_policy") != ORIENTATION_POLICY or
                gravity_meta.get("topology") != gravity_topology):
            raise ValueError("gravity cache lacks matching independently audited topology")
        offset_max = gravity_meta.get("offset_halving_slope_difference_deg_max")
        if offset_max is None:
            raise ValueError("gravity cache lacks the required offset-halving audit")
        if float(offset_max) > args.gravity_offset_tolerance:
            raise ValueError(
                "cached gravity slope offset audit %.6g deg exceeds the requested "
                "--gravity-offset-tolerance %.6g deg" %
                (float(offset_max), args.gravity_offset_tolerance))
        cache_status = "reused"
        print("reused gravity cache %s (%s reference omega)" %
              (gravity_cache_path,
               "same" if cached["reference_omega_reused"] else "recomputed"),
              flush=True)
    else:
        slope, gravity, centrifugal, gravity_meta = effective_slope(
            elements["centroid"], elements["normal"], gravity_vertices,
            gravity_faces, args.density, omega, args.gravity_offset,
            args.gravity_offset_tolerance,
            args.gravity_chunk, not args.serial_gravity)
        effective = gravity + centrifugal
        gravity_meta = dict(gravity_meta,
                            created_utc=dt.datetime.now(dt.timezone.utc).isoformat())
        cached = write_cache(
            gravity_cache_path, elements["centroid"], elements["normal"],
            elements["area"], gravity_definition, gravity, omega, gravity_meta)
        # Use the serialized representation returned by the cache writer so
        # this first build and every future reuse follow the identical path.
        gravity = cached["self_gravity"]
        centrifugal = cached["centrifugal"]
        effective = cached["effective"]
        slope = cached["slope"]
        cache_status = "created"
        print("wrote reusable gravity cache %s" % gravity_cache_path, flush=True)
    validation_identity = None
    if args.gravity_validation_shape:
        validation_path = Path(args.gravity_validation_shape).expanduser().resolve()
        validation_vertices, validation_faces = read_nastran_surface(validation_path)
        resolution = validate_gravity_resolution(
            elements["centroid"], elements["normal"], omega,
            args.gravity_offset, args.density, slope, validation_vertices,
            validation_faces, args.gravity_resolution_tolerance,
            not args.serial_gravity)
        validation_identity = file_identity(validation_path)
        gravity_meta["resolution_validation"] = {
            "shape": validation_identity,
            "vertices": int(len(validation_vertices)),
            "triangles": int(len(validation_faces)),
            **resolution,
        }
        print("gravity resolution validation: p95 %.4f deg, max %.4f deg" %
              (resolution["slope_difference_deg_p95"],
               resolution["slope_difference_deg_max"]), flush=True)
    template_names = ["late_ng_not_rg"]
    columns = [late]
    template_meta = {
        "late_ng_not_rg": dict(late_meta,
            area_fraction=area_fraction(late, elements["area"]),
            interpretation="geometry-defined candidate region; fitted amplitude is not proof of a physical source"),
    }
    for threshold in thresholds:
        label = ("cliff_ge_%gdeg" % threshold).replace(".", "p")
        mask = (slope >= threshold).astype(float)
        template_names.append(label)
        columns.append(mask)
        template_meta[label] = {
            "definition": "effective_gravitational_slope_deg >= %.17g" % threshold,
            "threshold_deg": threshold,
            "face_count": int(np.count_nonzero(mask)),
            "area_fraction": area_fraction(mask, elements["area"]),
            # These facets can change classification under tiny numerical
            # perturbations. Expose indices in the model's original face order
            # instead of silently hiding threshold sensitivity in area totals.
            "threshold_margin_deg": 1e-6,
            "faces_within_threshold_margin_zero_based": np.flatnonzero(
                np.abs(slope - threshold) <= 1e-6).tolist(),
            "interpretation": "cliff/steep-terrain hypothesis; amplitude fitted nonnegative; no production law asserted",
        }
        print("%s: %d facets, %.6f of surface area" %
              (label, np.count_nonzero(mask), template_meta[label]["area_fraction"]),
              flush=True)

    base_manifest = run / "run_manifest.json"
    metadata = {
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "base_run": str(run),
        "base_run_manifest": file_identity(base_manifest),
        "elements": file_identity(run / "elements.bin"),
        "templates": template_meta,
        "late_sensitivity_diagnostics": {
            "ng_integrated_min_max": [float(np.min(ng_response)), float(np.max(ng_response))],
            "rg_integrated_min_max": [float(np.min(rg_response)), float(np.max(rg_response))],
            "not_used_to_threshold_template": True,
        },
        "gravity": {
            "shape": gravity_shape_identity,
            "shape_vertices": int(len(gravity_vertices)),
            "shape_triangles": int(len(gravity_faces)),
            "density_kg_m3": args.density,
            "homogeneous_density_assumption": True,
            "evaluation_offset_m": args.gravity_offset,
            "acceleration_definition": "self_gravity - omega_cross_omega_cross_r",
            "omega_body_frame_rad_s": omega.tolist(),
            "omega_sample_epochs_utc": omega_epochs,
            "omega_max_sample_deviation_rad_s": omega_variation,
            "omega_command": [str(value) for value in omega_command],
            "slope_definition": "angle_between_inward_facet_normal_and_effective_gravity",
            "slope_deg_min_median_max": [float(np.min(slope)), float(np.median(slope)),
                                          float(np.max(slope))],
            "gravity_m_s2_min_max": [float(np.min(np.linalg.norm(gravity, axis=1))),
                                      float(np.max(np.linalg.norm(gravity, axis=1)))],
            "centrifugal_m_s2_min_max": [float(np.min(np.linalg.norm(centrifugal, axis=1))),
                                          float(np.max(np.linalg.norm(centrifugal, axis=1)))],
            "effective_gravity_m_s2_min_max": [
                float(np.min(np.linalg.norm(effective, axis=1))),
                float(np.max(np.linalg.norm(effective, axis=1)))],
            "cache": {
                "status": cache_status,
                "file": file_identity(gravity_cache_path),
                "definition": gravity_definition,
                "reference_omega_reused": bool(cached["reference_omega_reused"]),
            },
            **gravity_meta,
        },
        "scientific_scope": {
            "species": "H2O",
            "templates_are_predictors_not_established_causal_laws": True,
            "late_ng_transition_requires_profile_sensitivity_test": True,
        },
        "software": {"script": file_identity(__file__)},
    }
    payload = write_bundle(destination, template_names,
                           np.column_stack(columns), elements["centroid"],
                           elements["area"], metadata)
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    write_map(map_path, elements["centroid"], elements["area"], [
        ("late NG-visible / RG-not-visible mask", late, "viridis", (0.0, 1.0)),
        ("effective gravitational slope [deg]", slope, "magma", (0.0, 90.0)),
    ] + [(name, values, "viridis", (0.0, 1.0))
         for name, values in zip(template_names[1:], columns[1:])])
    print("wrote %s\nwrote %s\nwrote %s" % (destination, json_path, map_path))


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError) as exc:
        sys.exit("Error: %s" % exc)
