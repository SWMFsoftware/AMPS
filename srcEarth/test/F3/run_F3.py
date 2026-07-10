#!/usr/bin/env python3
"""
F3 — dipole cutoff-filtered flux validation test.

Run from the directory containing the AMPS executable:

    python srcEarth/test/F3/run_F3.py -np 4 -nt 16

F3 links the validated dipole cutoff calculation to the density/flux folding
path.  It runs gridless density/spectrum at explicit points on the 9000 km shell
and compares the output fluxes with a Størmer hard-cutoff power-law reference:

    F_ch = 4*pi*T_open*int_{max(E1,Ecut)}^{E2} J0*(E/E0)^(-gamma) dE

where Ecut is obtained by converting the vertical Størmer rigidity to proton
kinetic energy, and T_open is the analytic open-sky fraction outside the inner
absorbing sphere.  The comparison is approximate because AMPS evaluates a full
directional access map with penumbra, while the reference collapses access to a
vertical step function.
"""

from __future__ import print_function
import argparse
import csv
import json
import math
import os
import re
import shutil
import subprocess
import sys
import textwrap
from pathlib import Path

TEST_ID = "F3"
TEST_NAME = "Dipole cutoff-filtered flux"

RE_KM = 6371.2
R_INNER_RE = 1.01
R_INNER_KM = R_INNER_RE * RE_KM
ALT_KM_DEFAULT = 9000.0
STORMER_R0_GV = 14.9
MP_MEV = 938.2720813
C_LIGHT = 299792458.0
J0 = 1.0
E0 = 10.0
GAMMA = 3.5
EMIN = 1.0
EMAX = 1000.0
DEFAULT_LATS = [-70.0, -60.0, -45.0, -30.0, 0.0, 30.0, 45.0, 60.0, 70.0]
DEFAULT_LONS = [0.0, 90.0]
ENERGY_BINS = [
    ("BIN01_1_3", 1.0, 3.0),
    ("BIN02_3_10", 3.0, 10.0),
    ("BIN03_10_30", 10.0, 30.0),
    ("BIN04_30_100", 30.0, 100.0),
    ("BIN05_100_300", 100.0, 300.0),
    ("BIN06_300_1000", 300.0, 1000.0),
    ("FULL", 1.0, 1000.0),
]


def parse_float_list(text, name):
    out = []
    for part in text.split(','):
        part = part.strip()
        if not part:
            continue
        out.append(float(part))
    if not out:
        raise SystemExit("%s must contain at least one numeric value" % name)
    return out


def normalize_negative_list_options(argv, option_names):
    """Return argv with negative comma-list values attached to their options.

    argparse normally accepts ``--opt value`` for options with one argument, but
    a value such as ``-70,-60,0`` can be mistaken for another option because it
    begins with ``-``.  Users can always avoid that with ``--opt=-70,-60,0``.
    This helper also supports the more natural form used in the test README:

        --lats -70,-60,-45,-30,0,30,45,60,70

    It only rewrites values that look like negative numeric comma lists, so
    missing arguments and real options such as ``--skip-run`` are still reported
    by argparse in the usual way.
    """
    if argv is None:
        argv = sys.argv[1:]
    argv = list(argv)
    out = []
    i = 0
    negative_numeric_list = re.compile(r"^-([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][+-]?[0-9]+)?(,.*)?$")
    option_names = set(option_names)

    while i < len(argv):
        token = argv[i]
        if token in option_names and i + 1 < len(argv):
            value = argv[i + 1]
            if negative_numeric_list.match(value):
                out.append(token + "=" + value)
                i += 2
                continue
        out.append(token)
        i += 1

    return out


def point_label(lon, lat):
    def fmt(x):
        if abs(x - round(x)) < 1.0e-10:
            return str(int(round(x))).replace('-', 'm')
        return ("%.6g" % x).replace('-', 'm').replace('.', 'p')
    return "lon%s_lat%s" % (fmt(lon), fmt(lat))


def build_points(alt_km, lons, lats):
    r = RE_KM + alt_km
    points = []
    for lon in lons:
        lon_rad = math.radians(lon)
        for lat in lats:
            lat_rad = math.radians(lat)
            c = math.cos(lat_rad)
            x = r * c * math.cos(lon_rad)
            y = r * c * math.sin(lon_rad)
            z = r * math.sin(lat_rad)
            points.append({
                "label": point_label(lon, lat),
                "lon_deg": lon,
                "lat_deg": lat,
                "alt_km": alt_km,
                "x_km": x,
                "y_km": y,
                "z_km": z,
            })
    return points


def stormer_vertical_cutoff_gv(lat_deg, alt_km):
    r_re = (RE_KM + alt_km) / RE_KM
    return STORMER_R0_GV * (math.cos(math.radians(lat_deg)) ** 4) / (r_re * r_re)


def energy_from_rigidity_mev(R_GV):
    pc_mev = 1000.0 * abs(R_GV)
    return math.sqrt(pc_mev * pc_mev + MP_MEV * MP_MEV) - MP_MEV


def rigidity_from_energy_gv(E_MeV):
    pc_mev = math.sqrt(max(0.0, E_MeV * (E_MeV + 2.0 * MP_MEV)))
    return pc_mev / 1000.0


def open_sky_fraction(alt_km, r_inner_km):
    r = RE_KM + alt_km
    if r <= r_inner_km:
        return 0.0
    s = r_inner_km / r
    return 0.5 * (1.0 + math.sqrt(max(0.0, 1.0 - s * s)))


def j_boundary(E_MeV):
    return J0 * (E_MeV / E0) ** (-GAMMA)


def speed_from_energy(E_MeV):
    gamma = 1.0 + E_MeV / MP_MEV
    beta2 = max(0.0, 1.0 - 1.0 / (gamma * gamma))
    return C_LIGHT * math.sqrt(beta2)


def flux_power_law(E1, E2):
    if E2 <= E1:
        return 0.0
    if abs(GAMMA - 1.0) < 1.0e-14:
        integ = J0 * (E0 ** GAMMA) * math.log(E2 / E1)
    else:
        integ = J0 * (E0 ** GAMMA) / (GAMMA - 1.0) * (E1 ** (1.0 - GAMMA) - E2 ** (1.0 - GAMMA))
    return 4.0 * math.pi * integ


def density_high_resolution(E1, E2, n=200000):
    if E2 <= E1:
        return 0.0
    log1 = math.log(E1)
    log2 = math.log(E2)
    total = 0.0
    prev_E = None
    prev_g = None
    for i in range(n):
        a = float(i) / float(n - 1)
        E = math.exp(log1 + a * (log2 - log1))
        g = j_boundary(E) / speed_from_energy(E)
        if prev_E is not None:
            total += 0.5 * (prev_g + g) * (E - prev_E)
        prev_E = E
        prev_g = g
    return 4.0 * math.pi * total


def analytic_reference_for_point(point):
    rc = stormer_vertical_cutoff_gv(point["lat_deg"], point["alt_km"])
    ecut = energy_from_rigidity_mev(rc)
    topen = open_sky_fraction(point["alt_km"], R_INNER_KM)
    ref = {
        "Rc_vert_GV": rc,
        "Ecut_MeV": ecut,
        "T_open": topen,
        "T_at_Emax": topen if rigidity_from_energy_gv(EMAX) >= rc else 0.0,
        "density_total_m3": topen * density_high_resolution(max(EMIN, ecut), EMAX, n=40000),
        "flux_total_m2s1": topen * flux_power_law(max(EMIN, ecut), EMAX),
    }
    for name, e1, e2 in ENERGY_BINS:
        ref["flux_%s_m2s1" % name] = topen * flux_power_law(max(e1, ecut), e2)
    return ref


def compute_reference_rows(points):
    rows = []
    for p in points:
        ref = analytic_reference_for_point(p)
        common = {
            "point": p["label"],
            "lon_deg": p["lon_deg"],
            "lat_deg": p["lat_deg"],
            "alt_km": p["alt_km"],
            "reference_type": "stormer_hard_cutoff_open_sky",
        }
        quantities = [
            ("Rc_vert_GV", ref["Rc_vert_GV"], "GV", "vertical Størmer cutoff rigidity"),
            ("Ecut_MeV", ref["Ecut_MeV"], "MeV", "proton kinetic energy equivalent to Rc_vert"),
            ("T_open", ref["T_open"], "1", "analytic straight-line open-sky fraction"),
            ("T_at_Emax", ref["T_at_Emax"], "1", "hard-cutoff transmission at SPEC_EMAX"),
            ("density_total_m3", ref["density_total_m3"], "m^-3", "T_open*4*pi*int J(E)/v(E)dE above Ecut"),
            ("flux_total_m2s1", ref["flux_total_m2s1"], "m^-2 s^-1", "T_open*4*pi*int J(E)dE above Ecut"),
        ]
        for name, _, _ in ENERGY_BINS:
            quantities.append(("flux_%s_m2s1" % name, ref["flux_%s_m2s1" % name], "m^-2 s^-1",
                               "T_open*4*pi*int over channel above Ecut"))
        for q, value, units, note in quantities:
            row = dict(common)
            row.update({
                "quantity": q,
                "expected_value": value,
                "units": units,
                "notes": note,
            })
            rows.append(row)
    return rows


def write_reference_csv(path, points):
    rows = compute_reference_rows(points)
    fields = ["point", "lon_deg", "lat_deg", "alt_km", "quantity", "expected_value",
              "units", "reference_type", "notes"]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            out = dict(row)
            for key in ("lon_deg", "lat_deg", "alt_km", "expected_value"):
                out[key] = "%.17e" % float(out[key])
            w.writerow(out)


def reference_map_from_rows(points):
    refs = {}
    for p in points:
        ref = analytic_reference_for_point(p)
        refs[p["label"]] = ref
    return refs


def _norm_key(k):
    return k.lower().replace(" ", "_").replace("-", "_").replace("^", "").replace("/", "_")


def _parse_variables(line):
    return [x.strip() for x in re.findall(r'"([^"]+)"', line)]


def _skip_line(line):
    s = line.strip()
    if not s:
        return True
    u = s.upper()
    return s.startswith("#") or s.startswith("!") or u.startswith("TITLE") or u.startswith("AUXDATA")


def read_table_file(path):
    records = []
    variables = []
    in_zone = False
    with path.open("r", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            u = line.upper()
            if u.startswith("VARIABLES"):
                variables = _parse_variables(line)
                continue
            if u.startswith("ZONE"):
                in_zone = True
                continue
            if _skip_line(line) or not in_zone:
                continue
            parts = line.split()
            if len(parts) != len(variables):
                continue
            try:
                records.append(dict((k, float(v)) for k, v in zip(variables, parts)))
            except ValueError:
                pass
    return variables, records


def read_spectrum_file(path):
    zones = []
    variables = []
    current = None
    with path.open("r", errors="replace") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            u = line.upper()
            if u.startswith("VARIABLES"):
                variables = _parse_variables(line)
                continue
            if u.startswith("ZONE"):
                if current is not None:
                    zones.append(current)
                current = {k: [] for k in variables}
                continue
            if _skip_line(line) or current is None:
                continue
            parts = line.split()
            if len(parts) != len(variables):
                continue
            try:
                for k, v in zip(variables, parts):
                    current[k].append(float(v))
            except ValueError:
                pass
    if current is not None:
        zones.append(current)
    return zones


def find_col(record, candidates):
    norm = {}
    for k in record.keys():
        norm[_norm_key(k)] = k
    for c in candidates:
        key = _norm_key(c)
        if key in norm:
            return norm[key]
    for c in candidates:
        key = _norm_key(c)
        for nk, orig in norm.items():
            if key in nk:
                return orig
    return None


def trapz(xs, ys):
    if len(xs) < 2:
        return 0.0
    s = 0.0
    for i in range(len(xs) - 1):
        s += 0.5 * (ys[i] + ys[i + 1]) * (xs[i + 1] - xs[i])
    return s


def interp(xs, ys, x):
    if x <= xs[0]:
        return ys[0]
    if x >= xs[-1]:
        return ys[-1]
    lo = 0
    hi = len(xs) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if xs[mid] <= x:
            lo = mid
        else:
            hi = mid
    t = (x - xs[lo]) / (xs[hi] - xs[lo])
    return ys[lo] + t * (ys[hi] - ys[lo])


def interval_grid(E, E1, E2):
    lo = max(E1, E[0])
    hi = min(E2, E[-1])
    if lo >= hi:
        return []
    xs = [lo]
    for x in E:
        if x > lo and x < hi:
            xs.append(x)
    xs.append(hi)
    return xs


def flux_from_spectrum(E, T, E1, E2):
    xs = interval_grid(E, E1, E2)
    if len(xs) < 2:
        return 0.0
    ys = [interp(E, T, x) * j_boundary(x) for x in xs]
    return 4.0 * math.pi * trapz(xs, ys)


def density_from_spectrum(E, T, E1, E2):
    xs = interval_grid(E, E1, E2)
    if len(xs) < 2:
        return 0.0
    ys = [interp(E, T, x) * j_boundary(x) / speed_from_energy(x) for x in xs]
    return 4.0 * math.pi * trapz(xs, ys)


def max_rel_boundary_spectrum(E, Jb):
    m = 0.0
    for e, j in zip(E, Jb):
        ref = j_boundary(e)
        m = max(m, abs(j - ref) / max(abs(ref), 1.0e-300))
    return m


def max_abs_jlocal_residual(T, Jb, Jl):
    m = 0.0
    for t, jb, jl in zip(T, Jb, Jl):
        denom = max(abs(t * jb), 1.0e-300)
        m = max(m, abs(jl - t * jb) / denom)
    return m


def replace_key(text, key, value):
    return re.sub(r'(?m)^(\s*' + re.escape(key) + r'\s+)\S+', r'\g<1>' + str(value), text)


def points_block(points):
    lines = ["  OUTPUT_MODE POINTS", "  POINTS_BEGIN"]
    for p in points:
        lines.append("    ! %s lon=%g lat=%g alt=%g km" % (p["label"], p["lon_deg"], p["lat_deg"], p["alt_km"]))
        lines.append("    POINT  %.12e  %.12e  %.12e" % (p["x_km"], p["y_km"], p["z_km"]))
    lines.append("  POINTS_END")
    return "\n".join(lines)


def render_input(template_path, output_path, points, args):
    text = template_path.read_text()
    for key, value in [
        ("MODE3D_THREADS", args.nt),
        ("GRIDLESS_MPI_SCHEDULER", args.scheduler),
        ("GRIDLESS_MPI_DYNAMIC_CHUNK", args.dynamic_chunk),
        ("DS_TRANSMISSION_SCAN_N", args.scan_n),
        ("DS_NINTERVALS", args.nintervals),
        ("DS_MAX_TRAJ_TIME", args.max_trace_time),
        ("MAX_TRACE_TIME", args.max_trace_time),
        ("DT_TRACE", args.dt_trace),
    ]:
        text = replace_key(text, key, value)
    block = points_block(points)
    text = re.sub(r'(?ms)^#OUTPUT_DOMAIN\s+.*?^\s*POINTS_END\s*$', "#OUTPUT_DOMAIN\n" + block, text)
    text += ("\n! F3 harness settings\n"
             "! F3_ALT_KM %.17g\n"
             "! F3_LONS %s\n"
             "! F3_LATS %s\n"
             "! F3_NT %d\n"
             "! F3_GRIDLESS_MPI_SCHEDULER %s\n"
             "! F3_GRIDLESS_MPI_DYNAMIC_CHUNK %d\n"
             "! F3_DS_TRANSMISSION_SCAN_N %d\n"
             "! F3_DS_NINTERVALS %d\n") % (
                 args.alt_km,
                 ",".join("%.12g" % x for x in args.lons_values),
                 ",".join("%.12g" % x for x in args.lats_values),
                 args.nt, args.scheduler, args.dynamic_chunk, args.scan_n, args.nintervals)
    output_path.write_text(text)


def read_outputs(workdir, points):
    density_path = workdir / "gridless_points_density.dat"
    flux_path = workdir / "gridless_points_flux.dat"
    spectrum_path = workdir / "gridless_points_spectrum.dat"
    missing = [str(p) for p in (density_path, flux_path, spectrum_path) if not p.exists()]
    if missing:
        raise RuntimeError("Missing required AMPS output files: " + ", ".join(missing))
    _, density_rows = read_table_file(density_path)
    _, flux_rows = read_table_file(flux_path)
    spectra = read_spectrum_file(spectrum_path)
    if not (len(density_rows) == len(flux_rows) == len(spectra) == len(points)):
        raise RuntimeError("Output size mismatch: density=%d flux=%d spectra=%d expected=%d" %
                           (len(density_rows), len(flux_rows), len(spectra), len(points)))
    out = []
    for i, (p, drow, frow, spec) in enumerate(zip(points, density_rows, flux_rows, spectra)):
        n_col = find_col(drow, ["N_m^-3", "N_m3", "density", "n"])
        rc_eff_col = find_col(drow, ["Rc_effective_GV", "Rc_eff_GV", "Rc_effective"])
        rc_low_col = find_col(drow, ["Rc_lower_GV", "Rc_lower"])
        rc_up_col = find_col(drow, ["Rc_upper_GV", "Rc_upper"])
        t_high_col = find_col(drow, ["T_high", "THigh"])
        f_tot_col = find_col(frow, ["F_tot_m2s1", "F_total", "F_tot"])
        if n_col is None or f_tot_col is None:
            raise RuntimeError("Could not identify density or total-flux column in point %d" % i)
        for key in ("E_MeV", "T", "J_boundary_perMeV", "J_local_perMeV"):
            if key not in spec:
                raise RuntimeError("Spectrum output is missing column %s" % key)
        flux = {"F_tot_m2s1": frow[f_tot_col]}
        for name, e1, e2 in ENERGY_BINS:
            col = find_col(frow, ["F_%s_m2s1" % name, name])
            if col is not None:
                flux["F_%s_m2s1" % name] = frow[col]
        out.append({
            "point": p,
            "density": drow[n_col],
            "Rc_effective_GV": drow[rc_eff_col] if rc_eff_col is not None else float("nan"),
            "Rc_lower_GV": drow[rc_low_col] if rc_low_col is not None else float("nan"),
            "Rc_upper_GV": drow[rc_up_col] if rc_up_col is not None else float("nan"),
            "T_high": drow[t_high_col] if t_high_col is not None else float("nan"),
            "flux": flux,
            "E_MeV": spec["E_MeV"],
            "T": spec["T"],
            "J_boundary_perMeV": spec["J_boundary_perMeV"],
            "J_local_perMeV": spec["J_local_perMeV"],
        })
    return out


def add_check(rows, check, point, quantity, value, expected_value, rel_tol=None, abs_tol=None,
              units="", check_type="error_metric", note=""):
    abs_error = abs(value - expected_value)
    if expected_value == 0.0:
        rel_error = abs_error
    else:
        rel_error = abs_error / max(abs(expected_value), 1.0e-300)
    ok = True
    if rel_tol is not None:
        ok = ok and (rel_error <= rel_tol)
    if abs_tol is not None:
        ok = ok and (abs_error <= abs_tol)
    rows.append({
        "check": check,
        "point": point,
        "quantity": quantity,
        "check_type": check_type,
        "passed": ok,
        "value": value,
        "expected_value": expected_value,
        "abs_error": abs_error,
        "rel_error": rel_error,
        "rel_tol": rel_tol,
        "abs_tol": abs_tol,
        "units": units,
        "note": note,
    })
    return ok


def mean(values):
    return sum(values) / float(len(values)) if values else float("nan")


def analyze(workdir, points, reference_csv, args):
    refs = reference_map_from_rows(points)
    outputs = read_outputs(workdir, points)
    rows = []
    passed = True

    for data in outputs:
        p = data["point"]
        label = p["label"]
        ref = refs[label]
        E = data["E_MeV"]
        T = data["T"]
        Jb = data["J_boundary_perMeV"]
        Jl = data["J_local_perMeV"]
        if not (len(E) == len(T) == len(Jb) == len(Jl)) or len(E) < 2:
            raise RuntimeError("Malformed spectrum zone for %s" % label)

        passed = add_check(rows, "boundary_power_law", label, "max_rel_J_boundary", max_rel_boundary_spectrum(E, Jb),
                           0.0, abs_tol=args.differential_tol, units="1", check_type="setup_identity",
                           note="spectrum file must preserve the imposed power law") and passed
        passed = add_check(rows, "local_spectrum_closure", label, "max_rel_Jlocal_minus_TJb",
                           max_abs_jlocal_residual(T, Jb, Jl), 0.0, abs_tol=args.consistency_tol,
                           units="1", check_type="exact_internal_identity",
                           note="J_local(E) must equal T(E)*J_boundary(E)") and passed

        dens_spec = density_from_spectrum(E, T, EMIN, EMAX)
        passed = add_check(rows, "density_file_vs_spectrum", label, "density_total_relative_residual",
                           abs(data["density"] - dens_spec) / max(abs(dens_spec), 1.0e-300), 0.0,
                           abs_tol=args.consistency_tol, units="1", check_type="exact_internal_identity",
                           note="density file must match integration of saved T(E)J(E)/v(E)") and passed

        flux_spec_total = flux_from_spectrum(E, T, EMIN, EMAX)
        passed = add_check(rows, "flux_file_vs_spectrum", label, "flux_total_relative_residual",
                           abs(data["flux"]["F_tot_m2s1"] - flux_spec_total) / max(abs(flux_spec_total), 1.0e-300), 0.0,
                           abs_tol=args.consistency_tol, units="1", check_type="exact_internal_identity",
                           note="total flux file must match integration of saved T(E)J(E)") and passed

        # Analytical hard-cutoff comparison.  This is intentionally approximate;
        # the directional solver has penumbra and non-vertical access.
        analytic_tol = args.analytic_flux_tol
        for q_file, q_ref, units in [
            ("F_tot_m2s1", "flux_total_m2s1", "m^-2 s^-1"),
        ]:
            expected = ref[q_ref]
            value = data["flux"][q_file]
            if expected == 0.0:
                passed = add_check(rows, "stormer_hard_cutoff_flux", label, q_ref, value, expected,
                                   abs_tol=args.blocked_flux_abs_tol, units=units,
                                   check_type="approx_analytic_reference",
                                   note="vertical Størmer hard-cutoff flux; zero means Ecut exceeds channel upper edge") and passed
            else:
                passed = add_check(rows, "stormer_hard_cutoff_flux", label, q_ref, value, expected,
                                   rel_tol=analytic_tol, units=units,
                                   check_type="approx_analytic_reference",
                                   note="vertical Størmer hard-cutoff flux with analytic open-sky factor") and passed

        for name, e1, e2 in ENERGY_BINS:
            out_key = "F_%s_m2s1" % name
            ref_key = "flux_%s_m2s1" % name
            if out_key not in data["flux"]:
                passed = add_check(rows, "flux_channel_present", label, out_key, 0.0, 1.0,
                                   abs_tol=0.0, units="1", check_type="setup_identity",
                                   note="missing requested flux-channel column") and passed
                continue
            # Compare file-vs-spectrum for every channel tightly.
            f_spec = flux_from_spectrum(E, T, e1, e2)
            passed = add_check(rows, "flux_channel_file_vs_spectrum", label,
                               "%s_relative_residual" % out_key,
                               abs(data["flux"][out_key] - f_spec) / max(abs(f_spec), 1.0e-300),
                               0.0, abs_tol=args.consistency_tol, units="1",
                               check_type="exact_internal_identity",
                               note="channel flux file must match integration of saved T(E)J(E)") and passed
            # Compare analytical channel only when it is not nearly zero.  Tiny
            # high-cutoff channels are better covered by the blocked-total check.
            expected = ref[ref_key]
            if expected > args.channel_min_reference:
                passed = add_check(rows, "stormer_hard_cutoff_channel_flux", label, ref_key,
                                   data["flux"][out_key], expected, rel_tol=args.analytic_channel_tol,
                                   units="m^-2 s^-1", check_type="approx_analytic_reference",
                                   note="channel flux compared with vertical Størmer hard cutoff") and passed

        passed = add_check(rows, "stormer_cutoff_diagnostic", label, "Rc_effective_minus_Rc_vert_GV",
                           data["Rc_effective_GV"], ref["Rc_vert_GV"], rel_tol=args.rc_tol,
                           units="GV", check_type="approx_analytic_reference",
                           note="effective cutoff from T(E) should be close to vertical Størmer cutoff away from penumbra") and passed
        passed = add_check(rows, "high_energy_transmission", label, "T_high", data["T_high"],
                           ref["T_at_Emax"], abs_tol=args.t_high_tol, units="1",
                           check_type="approx_analytic_reference",
                           note="T at SPEC_EMAX compared with hard-cutoff open-sky value") and passed

    # Symmetry and trend checks use actual AMPS output and are exact for a centered
    # aligned dipole up to directional quadrature/noise.
    by_lat = {}
    by_lon_lat = {}
    for data in outputs:
        lat = data["point"]["lat_deg"]
        lon = data["point"]["lon_deg"]
        by_lat.setdefault(lat, []).append(data)
        by_lon_lat[(lon, lat)] = data

    for lat, vals in sorted(by_lat.items()):
        flux_vals = [v["flux"]["F_tot_m2s1"] for v in vals]
        spread = (max(flux_vals) - min(flux_vals)) / max(abs(mean(flux_vals)), 1.0e-300)
        passed = add_check(rows, "dipole_longitude_symmetry", "lat_%g" % lat,
                           "flux_total_lon_spread", spread, 0.0, abs_tol=args.longitude_tol,
                           units="1", check_type="symmetry_identity",
                           note="centered aligned dipole should be longitude independent") and passed

    lons = sorted(set(p["lon_deg"] for p in points))
    lats = sorted(set(p["lat_deg"] for p in points))
    for lon in lons:
        for lat in lats:
            if lat <= 0.0:
                continue
            a = by_lon_lat.get((lon, lat))
            b = by_lon_lat.get((lon, -lat))
            if a is None or b is None:
                continue
            fa = a["flux"]["F_tot_m2s1"]
            fb = b["flux"]["F_tot_m2s1"]
            rel = abs(fa - fb) / max(abs(0.5 * (fa + fb)), 1.0e-300)
            passed = add_check(rows, "dipole_north_south_symmetry", "lon_%g_lat_%g" % (lon, lat),
                               "flux_total_NS_relative_difference", rel, 0.0,
                               abs_tol=args.ns_tol, units="1", check_type="symmetry_identity",
                               note="centered aligned dipole should be symmetric in north/south latitude") and passed

    abs_lat_means = []
    for alat in sorted(set(abs(x) for x in lats)):
        vals = []
        for lat in (alat, -alat):
            vals.extend([v["flux"]["F_tot_m2s1"] for v in by_lat.get(lat, [])])
        if vals:
            abs_lat_means.append((alat, mean(vals)))
    # Flux should generally increase toward high absolute latitude because the
    # Størmer cutoff decreases.  Allow tiny numerical non-monotonicity.
    prev = None
    for alat, fmean in abs_lat_means:
        if prev is not None:
            alat_prev, fprev = prev
            residual = max(0.0, fprev - fmean) / max(abs(fprev), abs(fmean), 1.0e-300)
            passed = add_check(rows, "latitude_access_trend", "abs_lat_%g" % alat,
                               "flux_total_monotonic_residual", residual, 0.0,
                               abs_tol=args.trend_tol, units="1", check_type="physical_trend",
                               note="mean flux should not decrease as |latitude| increases") and passed
        prev = (alat, fmean)

    summary_csv = workdir / "F3_summary.csv"
    fields = ["check", "point", "quantity", "check_type", "passed", "value",
              "expected_value", "abs_error", "rel_error", "rel_tol", "abs_tol",
              "units", "note"]
    with summary_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for row in rows:
            w.writerow(row)

    result = {
        "test_id": TEST_ID,
        "test_name": TEST_NAME,
        "passed": passed,
        "n_points": len(points),
        "n_checks": len(rows),
        "n_failed": sum(1 for r in rows if not r["passed"]),
        "workdir": str(workdir),
        "summary_csv": str(summary_csv),
        "reference_csv": str(reference_csv),
        "settings": {
            "alt_km": args.alt_km,
            "lons": args.lons_values,
            "lats": args.lats_values,
            "np": args.np,
            "nt": args.nt,
            "scheduler": args.scheduler,
            "dynamic_chunk": args.dynamic_chunk,
            "scan_n": args.scan_n,
            "nintervals": args.nintervals,
        },
        "tolerances": {
            "analytic_flux_tol": args.analytic_flux_tol,
            "analytic_channel_tol": args.analytic_channel_tol,
            "blocked_flux_abs_tol": args.blocked_flux_abs_tol,
            "rc_tol": args.rc_tol,
            "t_high_tol": args.t_high_tol,
            "longitude_tol": args.longitude_tol,
            "ns_tol": args.ns_tol,
            "trend_tol": args.trend_tol,
            "consistency_tol": args.consistency_tol,
            "differential_tol": args.differential_tol,
        },
    }
    result_json = workdir / "F3_result.json"
    with result_json.open("w") as f:
        json.dump(result, f, indent=2, sort_keys=True)
        f.write("\n")

    print("F3 %s" % ("PASS" if passed else "FAIL"))
    print("Summary:", summary_csv)
    print("Result :", result_json)
    if not passed:
        print("Failed checks:")
        for r in rows:
            if not r["passed"]:
                print("  - {check} [{point}] {quantity}: value={value:.6e}, expected={expected_value:.6e}, rel={rel_error:.3e}, abs={abs_error:.3e}".format(**r))
    return passed


def parse_args(argv=None):
    p = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""
        F3: dipole cutoff-filtered flux.

        The test renders points on the 9000 km shell, runs gridless
        density/spectrum with FIELD_MODEL=DIPOLE and DS_TRANSMISSION_MODE=SCAN,
        and compares fluxes with a Størmer hard-cutoff power-law reference.
        """),
        epilog=textwrap.dedent("""
        Examples:
          python srcEarth/test/F3/run_F3.py
          python srcEarth/test/F3/run_F3.py -np 4 -nt 16 --scan-n 100
          python srcEarth/test/F3/run_F3.py --lons 0,90,180,270 --lats -70,-60,-45,-30,0,30,45,60,70
          python srcEarth/test/F3/run_F3.py --analytic-flux-tol 0.75 --rc-tol 0.5
          python srcEarth/test/F3/run_F3.py --dry-run --workdir test_output/F3_dryrun
          python srcEarth/test/F3/run_F3.py --skip-run --workdir test_output/F3_gridless
        """))
    p.add_argument("-np", type=int, default=4, help="number of MPI ranks; default: 4")
    p.add_argument("-nt", type=int, default=16, help="threads per rank; default: 16")
    p.add_argument("--amps", default="./amps", help="AMPS executable path; default: ./amps")
    p.add_argument("--mpirun", default="mpirun", help="MPI launcher; default: mpirun")
    p.add_argument("--workdir", default="test_output/F3_gridless", help="run/output directory")
    p.add_argument("--scheduler", choices=["DYNAMIC", "BLOCK_CYCLIC", "STATIC"], default="DYNAMIC")
    p.add_argument("--dynamic-chunk", type=int, default=0)
    p.add_argument("--scan-n", type=int, default=100, help="DS_TRANSMISSION_SCAN_N value; default: 100")
    p.add_argument("--nintervals", type=int, default=120, help="DS_NINTERVALS value; default: 120")
    p.add_argument("--alt-km", type=float, default=ALT_KM_DEFAULT, help="shell altitude in km; default: 9000")
    p.add_argument("--lons", default=",".join(str(x) for x in DEFAULT_LONS), help="comma-separated longitudes in degrees")
    p.add_argument("--lats", default=",".join(str(x) for x in DEFAULT_LATS), help="comma-separated latitudes in degrees")
    p.add_argument("--dt-trace", type=float, default=0.1, help="DT_TRACE value")
    p.add_argument("--max-trace-time", type=float, default=900.0, help="MAX_TRACE_TIME / DS_MAX_TRAJ_TIME value")
    p.add_argument("--analytic-flux-tol", type=float, default=0.75,
                   help="relative tolerance for approximate total flux vs Størmer hard-cutoff reference")
    p.add_argument("--analytic-channel-tol", type=float, default=0.75,
                   help="relative tolerance for nonzero channel fluxes vs Størmer hard-cutoff reference")
    p.add_argument("--channel-min-reference", type=float, default=1.0e-8,
                   help="skip approximate channel comparison when analytical channel flux is below this value")
    p.add_argument("--blocked-flux-abs-tol", type=float, default=1.0e-8,
                   help="absolute tolerance when analytical total flux is zero")
    p.add_argument("--rc-tol", type=float, default=0.75,
                   help="relative tolerance for Rc_effective vs vertical Størmer Rc")
    p.add_argument("--t-high-tol", type=float, default=0.35,
                   help="absolute tolerance for T_high vs hard-cutoff high-energy transmission")
    p.add_argument("--longitude-tol", type=float, default=5.0e-2,
                   help="allowed relative longitude spread in total flux at each latitude")
    p.add_argument("--ns-tol", type=float, default=5.0e-2,
                   help="allowed north/south relative difference in total flux")
    p.add_argument("--trend-tol", type=float, default=5.0e-2,
                   help="allowed monotonicity residual in the flux-vs-|latitude| trend")
    p.add_argument("--consistency-tol", type=float, default=2.0e-5,
                   help="tolerance for exact file-vs-spectrum consistency checks; allows Tecplot ASCII precision")
    p.add_argument("--differential-tol", type=float, default=2.0e-5,
                   help="tolerance for imposed power-law values in the spectrum file; allows Tecplot ASCII precision")
    p.add_argument("--skip-run", action="store_true", help="analyze existing files in --workdir")
    p.add_argument("--keep", action="store_true", help="do not delete an existing workdir")
    p.add_argument("--dry-run", action="store_true", help="render input and print the AMPS command without executing it")
    normalized_argv = normalize_negative_list_options(argv, ["--lats", "--lons"])
    return p.parse_args(normalized_argv)


def main(argv=None):
    args = parse_args(argv)
    if args.np < 1:
        raise SystemExit("-np must be >= 1")
    if args.nt < 1:
        raise SystemExit("-nt must be >= 1")
    if args.dynamic_chunk < 0:
        raise SystemExit("--dynamic-chunk must be >= 0")
    if args.scan_n < 2:
        raise SystemExit("--scan-n must be >= 2")
    if args.nintervals < 2:
        raise SystemExit("--nintervals must be >= 2")
    if args.alt_km <= 0.0:
        raise SystemExit("--alt-km must be positive")
    args.lons_values = parse_float_list(args.lons, "--lons")
    args.lats_values = parse_float_list(args.lats, "--lats")
    for lat in args.lats_values:
        if lat < -90.0 or lat > 90.0:
            raise SystemExit("latitude out of range: %g" % lat)

    points = build_points(args.alt_km, args.lons_values, args.lats_values)
    launch_dir = Path.cwd().resolve()
    script_dir = Path(__file__).resolve().parent
    template = script_dir / "AMPS_PARAM_F3_gridless.in"
    reference_repo_csv = script_dir / "reference_F3_dipole_cutoff.csv"
    workdir = (launch_dir / args.workdir).resolve()

    if not args.skip_run:
        if workdir.exists() and not args.keep:
            shutil.rmtree(str(workdir))
        workdir.mkdir(parents=True, exist_ok=True)
        render_input(template, workdir / "AMPS_PARAM_F3.in", points, args)
    else:
        if not workdir.exists():
            raise SystemExit("--skip-run requested but workdir does not exist: %s" % workdir)
    workdir.mkdir(parents=True, exist_ok=True)

    reference_used_csv = workdir / "reference_F3_dipole_cutoff_used.csv"
    write_reference_csv(reference_used_csv, points)
    if reference_repo_csv.exists() and not args.skip_run:
        shutil.copy2(str(reference_repo_csv), str(workdir / "reference_F3_dipole_cutoff_repository.csv"))

    if not args.skip_run:
        amps_path = Path(args.amps)
        if not amps_path.is_absolute():
            amps_path = (launch_dir / amps_path).resolve()
        cmd = [
            args.mpirun, "-np", str(args.np), str(amps_path),
            "-mode", "gridless",
            "-i", "AMPS_PARAM_F3.in",
            "-mode3d-parallel", "THREADS",
            "-mode3d-threads", str(args.nt),
            "-mode3d-mpi-scheduler", args.scheduler,
            "-mode3d-mpi-dynamic-chunk", str(args.dynamic_chunk),
            "-density-transmission", "SCAN",
        ]
        print("Running:", " ".join(cmd))
        if args.dry_run:
            print("F3 dry run complete. Generated input under %s" % workdir)
            print("Reference:", reference_used_csv)
            return 0
        with (workdir / "F3_amps.log").open("w") as log:
            log.write("# Command: %s\n" % " ".join(cmd))
            log.flush()
            rc = subprocess.call(cmd, cwd=str(workdir), stdout=log, stderr=subprocess.STDOUT)
        if rc != 0:
            print("AMPS failed with exit code %d; see %s" % (rc, workdir / "F3_amps.log"))
            return rc

    ok = analyze(workdir, points, reference_used_csv, args)
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
