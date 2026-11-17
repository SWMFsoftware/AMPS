#!/usr/bin/env python3
"""Independent check of the drag coefficients of the AMPS sphere drag-coefficient test.

The test (MakefileTest/Table, "Sphere-drag-coefficient") runs one case of input/test/sphere-drag-coefficient.in per
gas-surface interaction model and writes PT/plots/drag-coefficient-test.dat.  This script computes the drag
coefficient of the same free-molecular flow about the same triangulated body from the definitions of the models,
by numerical integration of the mean momentum of the re-emitted molecules over the incident velocity distribution
on every face, and compares it with the AMPS value:

    PASS|FAIL <case> Cd=<AMPS> reference=<this script> sigma=<AMPS statistical sigma> deviation=<n> sigma (limit <k>)

The free-molecular force on a face of a convex body (no face shadows another, a re-emitted molecule does not return)
with the outward normal n and the area A, in a stream of number density n_inf, velocity U and temperature T, is

    F = A m n_inf  int_{v.n<0} |v.n| f(v) (v - <u>(v)) d^3v ,

with f the Maxwellian of the stream and <u>(v) the mean velocity that the model re-emits for the incident velocity v.
The drag is F.U/|U| summed over the faces, and Cd = 2 F_D/(rho U^2 A_proj), with A_proj the projected area that AMPS
uses (Orbiter::CalculateProjectionArea, a pixel map; the exact projected area of the faces is reported too).  On a face
the stream is a drifting Maxwellian with the normal drift W = U cos(theta) toward the wall and the tangential drift
U sin(theta) (theta: the angle between -n and U); the normal momentum flux (pressure P) and the tangential one
(shear S, along the tangential drift) depend only on theta, so they are tabulated in theta and interpolated.
The models (src/models/surface/Surface.h), with a = sqrt(2kTw/m), w = -v.n > 0, v_t the tangential velocity:

    adsorption  <u> = 0
    specular    <u> = v + 2 w n                              (P = 2 P_in, S = 0)
    diffuse     <u> = sqrt(pi)/2 a n
    maxwell     alpha diffuse + (1-alpha) specular
    cll         <u.n> = a E|w_m + sqrt(alpha_n/2) Z|, w_m = sqrt(1-alpha_n) w/a, Z a complex standard normal variable
                (the mean of the Rice distribution, closed form with Bessel functions); <u_t> = (1-sigma_t) v_t
    nocilla     <u.n> = c_r <x>, <x> = int x^2 e^{-(x-s)^2} / int x e^{-(x-s)^2}, s = U_r cos(theta_r)/c_r,
                c_r = sqrt(2kT_r/m) (closed form with erf); <u_t> = U_r sin(theta_r) v_t/|v_t|;
                with accommodation coefficients (Cook & Hoffbauer 1997, Eqs. 1, 2): <u.n> = (1-sigma') w +
                sigma' sqrt(pi)/2 a, <u_t> = (1-sigma) v_t
    schamberg   <u> = V_r (M_n(theta_r) n + M_t(theta_r) v_t/|v_t|), V_r^2 = (1-alpha) |v|^2 + alpha 4kTw/m,
                sin(theta_r) = (|v_t|/|v|)^nu; M_n, M_t: the mean direction of the beam ~ cos(pi phi/(2 phi0)) per unit
                solid angle within phi0 of the axis and above the surface (the part below the surface is excluded,
                as in AMPS), by quadrature over phi with the azimuth integrated analytically

The incident-velocity integrals use Gauss-Legendre quadrature in w on [0, max(W,0)+8c], and, where the direction of
v_t enters, polar coordinates in the tangent plane (Gauss-Legendre in |v_t|, trapezoidal in the azimuth); doubling the
nodes changes the result by less than 1e-9.  The flow, the model and its parameters are read from the runtime input
(the case overlaid on the common block, as AMPS::Input::ResolveCase does); the values that AMPS echoes in its result
file must agree with them, so the check fails if AMPS ran another input.  The mass and the Boltzmann constant are those
that AMPS reports (18 amu, 1.3806503e-23 J/K).

Usage:
  sphere_drag_check.py --in sphere-drag-coefficient.in --conf sphere-drag-coefficient.conf --mesh sphere.bdf \
      --result CASE=PATH [--result CASE=PATH ...] [--sigmas 4] [--report FILE]
The exit status is 0 when every case passes.
"""

import argparse
import math
import re
import sys

import numpy as np

SQRT_PI = math.sqrt(math.pi)


# ------------------------------------------------------------------------------------------------- runtime input
class Block:
    def __init__(self, kind="", name=""):
        self.kind, self.name, self.keys, self.blocks = kind, name, {}, []

    def find(self, kind, name=""):
        for b in self.blocks:
            if b.kind == kind and b.name == name:
                return b
        return None


def parse_input(path):
    """The subset of the AMPS runtime-input grammar (src/INPUT_FILE.md) that the test file uses: '!' comments,
    '\\' continuation, '#begin KIND[: NAME]', '#end KIND' and 'KEY = VALUE' (no #include)."""
    root = Block()
    stack = [root]
    pending = ""
    with open(path) as f:
        for raw in f:
            line = raw.split("!", 1)[0].rstrip()
            if line.endswith("\\"):
                pending += line[:-1] + " "
                continue
            line = (pending + line).strip()
            pending = ""
            if not line:
                continue
            m = re.match(r"#begin\s+([a-z][a-z0-9_]*)\s*(?::\s*([A-Za-z0-9_-]+))?$", line, re.I)
            if m:
                b = Block(m.group(1).lower(), m.group(2) or "")
                stack[-1].blocks.append(b)
                stack.append(b)
                continue
            m = re.match(r"#end\s+([a-z][a-z0-9_]*)$", line, re.I)
            if m:
                if len(stack) < 2 or stack[-1].kind != m.group(1).lower():
                    raise ValueError("%s: unbalanced '%s'" % (path, line))
                stack.pop()
                continue
            if line.startswith("#"):
                raise ValueError("%s: directive not supported by this checker: '%s'" % (path, line))
            m = re.match(r"([A-Za-z][A-Za-z0-9_.-]*)\s*=\s*(.+)$", line)
            if not m:
                raise ValueError("%s: cannot parse '%s'" % (path, line))
            stack[-1].keys[m.group(1)] = m.group(2).strip()
    if len(stack) != 1:
        raise ValueError("%s: a block is not closed" % path)
    return root


def overlay(base, over):
    """ResolveCase: the case replaces the keys of the common block and overlays the matching child blocks."""
    res = Block(base.kind, base.name)
    res.keys = dict(base.keys)
    res.blocks = [b for b in base.blocks if b.kind != "case"]
    res.keys.update(over.keys)
    for ob in over.blocks:
        target = res.find(ob.kind, ob.name)
        if target is None:
            res.blocks.append(ob)
        else:
            res.blocks[res.blocks.index(target)] = overlay(target, ob)
    return res


def resolve_case(path, case):
    root = parse_input(path)
    app = root.find("application", "orbiter")
    if app is None:
        raise ValueError("%s: no 'application: orbiter' block" % path)
    sel = app.find("case", case)
    if sel is None:
        raise ValueError("%s: no case '%s'" % (path, case))
    return overlay(app, sel)


def case_parameters(path, case, species):
    app = resolve_case(path, case)
    up = app.find("upstream")
    gsi = app.find("physics", "gas_surface_interaction")
    sp = gsi.find("species", species)
    p = {
        "velocity": np.array([float(x) for x in up.keys["velocity_m_s"].replace(",", " ").split()]),
        "temperature": float(up.keys["temperature_k"]),
        "density": float(up.find("species", species).keys["number_density_m-3"]),
        "model": gsi.keys["model"],
        "wall_temperature": float(gsi.keys["wall_temperature_k"]) if "wall_temperature_k" in gsi.keys else None,
        "par": {k: float(v) for k, v in (sp.keys.items() if sp else [])},
    }
    for k, v in gsi.keys.items():
        m = re.match(r"species\.%s\.(.+)$" % re.escape(species), k)
        if m:
            p["par"][m.group(1)] = float(v)
    return p


def read_result(path):
    res = {}
    with open(path) as f:
        for line in f:
            if line.startswith("#") or "=" not in line:
                continue
            k, v = line.rstrip("\n").split("=", 1)
            res[k.strip()] = v.strip()
    return res


# ------------------------------------------------------------------------------------------------- the body
def read_conf_scale(path):
    with open(path) as f:
        for line in f:
            m = re.match(r"\s*SurfaceMeshScalingFactor\s*=\s*([-+0-9.eE]+)", line.split("!", 1)[0])
            if m:
                return float(m.group(1))
    raise ValueError("%s: no SurfaceMeshScalingFactor" % path)


def read_nastran(path, scale):
    """GRID* (long format, continuation '*') and CTRIA3 records of a NASTRAN file."""
    nodes, tris = {}, []
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        l = lines[i]
        if l.startswith("GRID*"):
            nid = int(l[8:24])
            x1, x2 = float(l[40:56]), float(l[56:72])
            x3 = float(lines[i + 1][8:24])
            nodes[nid] = (x1, x2, x3)
            i += 2
            continue
        if l.startswith("CTRIA3"):
            f8 = [l[k:k + 8] for k in range(0, 48, 8)]
            tris.append((int(f8[3]), int(f8[4]), int(f8[5])))
        i += 1
    ids = sorted(nodes)
    index = {n: k for k, n in enumerate(ids)}
    xyz = np.array([nodes[n] for n in ids]) * abs(scale)
    tri = np.array([[index[a], index[b], index[c]] for a, b, c in tris])
    return xyz, tri


def faces(xyz, tri):
    """Areas and outward unit normals (oriented away from the centroid of the convex body)."""
    p0, p1, p2 = xyz[tri[:, 0]], xyz[tri[:, 1]], xyz[tri[:, 2]]
    c = np.cross(p1 - p0, p2 - p0)
    area = 0.5 * np.linalg.norm(c, axis=1)
    n = c / (2.0 * area[:, None])
    centre = xyz.mean(axis=0)
    out = np.einsum("ij,ij->i", n, (p0 + p1 + p2) / 3.0 - centre)
    n[out < 0] *= -1.0
    return area, n


# ------------------------------------------------------------------------------------------------- quadrature
GL64 = np.polynomial.legendre.leggauss(64)
GL48 = np.polynomial.legendre.leggauss(48)


def gl(a, b, rule=GL64):
    x, w = rule
    return 0.5 * (b - a) * x + 0.5 * (b + a), 0.5 * (b - a) * w


def normal_nodes(W, c):
    """Nodes and weights of the incident normal speed w > 0 with the weight w g(w), g the 1D Maxwellian of mean W."""
    w, ww = gl(0.0, max(W, 0.0) + 8.0 * c)
    g = np.exp(-((w - W) / c) ** 2) / (SQRT_PI * c)
    return w, ww * w * g          # number flux per unit density at the nodes


def tangent_polar_nodes(Ut, c, nr=96, npsi=128):
    """Nodes (vt, vb, rho) and weights of the 2D Maxwellian in the tangent plane (mean (Ut,0), variance c^2/2 per
    component) in polar coordinates about the origin; the weights sum to 1."""
    rmax = Ut + 8.0 * c
    r, wr = gl(0.0, rmax, np.polynomial.legendre.leggauss(nr))
    psi = (np.arange(npsi) + 0.5) * 2.0 * math.pi / npsi
    R, PSI = np.meshgrid(r, psi, indexing="ij")
    vt, vb = R * np.cos(PSI), R * np.sin(PSI)
    dens = np.exp(-((vt - Ut) ** 2 + vb ** 2) / c ** 2) / (math.pi * c ** 2)
    wgt = dens * R * wr[:, None] * (2.0 * math.pi / npsi)
    return vt.ravel(), vb.ravel(), R.ravel(), wgt.ravel()


def rice_mean(nu, sig):
    """E|nu + sig (X + iY)|, X, Y standard normal (the mean of the Rice distribution), for an array nu >= 0:
    sig sqrt(pi/2) L_{1/2}(-nu^2/(2 sig^2)) = sig sqrt(pi/2) ((1+2y) I0e(y) + 2y I1e(y)), y = nu^2/(4 sig^2)
    (I0e, I1e: exponentially scaled modified Bessel functions, scipy).  Without scipy: quadrature in polar coordinates
    about the origin when nu < 6 sig (the distribution is a wide ring), about nu otherwise (a narrow peak)."""
    nu = np.asarray(nu, dtype=float)
    try:
        from scipy.special import i0e, i1e
        y = nu * nu / (4.0 * sig * sig)
        return sig * math.sqrt(0.5 * math.pi) * ((1 + 2 * y) * i0e(y) + 2 * y * i1e(y))
    except ImportError:
        pass
    out = np.empty_like(nu)
    nphi = 256
    phi = (np.arange(nphi) + 0.5) * 2.0 * math.pi / nphi
    for k, v in enumerate(nu):
        if v < 6.0 * sig:
            r, wr = gl(0.0, v + 10.0 * sig, np.polynomial.legendre.leggauss(128))
            R, PH = np.meshgrid(r, phi, indexing="ij")
            dens = np.exp(-((R * np.cos(PH) - v) ** 2 + (R * np.sin(PH)) ** 2) / (2 * sig ** 2)) / (2 * math.pi * sig ** 2)
            out[k] = np.sum(R * dens * R * wr[:, None]) * 2.0 * math.pi / nphi
        else:
            r, wr = gl(0.0, 10.0, np.polynomial.legendre.leggauss(128))
            R, PH = np.meshgrid(r, phi, indexing="ij")
            mod = np.sqrt(v * v + (sig * R) ** 2 + 2 * v * sig * R * np.cos(PH))
            out[k] = np.sum(mod * R * np.exp(-0.5 * R * R) * wr[:, None]) / nphi
    return out


def schamberg_beam_mean(theta_r, phi0):
    """(M_n, M_t): the mean unit vector of the beam ~ cos(pi phi/(2 phi0)) per solid angle, phi <= phi0, about the axis
    at theta_r from the normal, without the directions below the surface; components along n and along the tangential
    direction of the axis.  For a direction d = cos(phi) a + sin(phi) (cos(psi) b1 + sin(psi) b2), a = (cos tr, sin tr),
    b1 = (-sin tr, cos tr): d.n = cos(phi) cos(tr) - sin(phi) sin(tr) cos(psi) > 0 is the arc cos(psi) < C."""
    if phi0 <= 0.0:
        return math.cos(theta_r), math.sin(theta_r)
    ct, st = math.cos(theta_r), math.sin(theta_r)
    # split the phi integral where the truncation starts (phi* = pi/2 - theta_r)
    edges = [0.0, phi0]
    ps = 0.5 * math.pi - theta_r
    if 0.0 < ps < phi0:
        edges = [0.0, ps, phi0]
    N = Mn = Mt = 0.0
    for a, b in zip(edges[:-1], edges[1:]):
        phi, wphi = gl(a, b)
        w = np.cos(0.5 * math.pi * phi / phi0) * np.sin(phi) * wphi
        cph, sph = np.cos(phi), np.sin(phi)
        with np.errstate(divide="ignore", invalid="ignore"):
            C = np.where(sph * st > 0, cph * ct / (sph * st), np.inf)
        L = np.where(C >= 1, 2 * math.pi, np.where(C <= -1, 0.0, 2 * math.pi - 2 * np.arccos(np.clip(C, -1, 1))))
        J = np.where((C >= 1) | (C <= -1), 0.0, -2 * np.sqrt(np.clip(1 - C * C, 0, 1)))
        N += np.sum(w * L)
        Mn += np.sum(w * (cph * ct * L - sph * st * J))
        Mt += np.sum(w * (cph * st * L + sph * ct * J))
    return Mn / N, Mt / N


# ------------------------------------------------------------------------------------------------- the models
def plate(cos_th, U, T, mass, kB, Tw, model, par, beam_table=None):
    """Pressure P (momentum flux into the wall) and shear S (along the tangential drift) per unit area and per unit
    number density, divided by the mass."""
    c = math.sqrt(2 * kB * T / mass)
    sin_th = math.sqrt(max(0.0, 1.0 - cos_th * cos_th))
    W, Ut = U * cos_th, U * sin_th
    w, fw = normal_nodes(W, c)                 # flux per unit density
    Ndot = fw.sum()
    P_in = np.sum(fw * w)
    S_in = Ndot * Ut
    aw = math.sqrt(2 * kB * Tw / mass) if Tw else None

    if model == "adsorption":
        return P_in, S_in
    if model == "specular":
        return 2 * P_in, 0.0
    if model == "diffuse":
        return P_in + Ndot * SQRT_PI / 2 * aw, S_in
    if model == "maxwell":
        al = par["accommodation"]
        return P_in + al * Ndot * SQRT_PI / 2 * aw + (1 - al) * P_in, al * S_in
    if model in ("cll", "quasi_specular_cll"):
        an, st = par["normal_energy_accommodation"], par["tangential_momentum_accommodation"]
        nu = math.sqrt(1 - an) * w / aw
        un = aw * (rice_mean(nu, math.sqrt(an / 2)) if an > 0 else nu)
        return P_in + np.sum(fw * un), st * S_in
    if model == "nocilla" and "normal_momentum_accommodation" in par:
        # accommodation coefficients: the mean re-emitted velocity is (1-sigma')|v.n| + sigma' sqrt(pi)/2 a along n and
        # (1-sigma) v_t (Cook & Hoffbauer 1997, Eqs. 1 and 2); the energy coefficient does not change the mean
        sn, st = par["normal_momentum_accommodation"], par["tangential_momentum_accommodation"]
        return P_in + (1 - sn) * P_in + sn * Ndot * SQRT_PI / 2 * aw, st * S_in
    if model == "nocilla":
        Ur, thr, Tr = par["drift_speed_m_s"], math.radians(par["drift_angle_deg"]), par["temperature_k"]
        cr = math.sqrt(2 * kB * Tr / mass)
        s = Ur * math.cos(thr) / cr
        E = 1 + math.erf(s)
        mean_x = (s * math.exp(-s * s) + (0.5 + s * s) * SQRT_PI * E) / (math.exp(-s * s) + s * SQRT_PI * E)
        vt, vb, rho, wt = tangent_polar_nodes(Ut, c)
        ecos = np.sum(wt * vt / rho)           # <cos psi> of the incident tangential direction
        return P_in + Ndot * cr * mean_x, S_in - Ndot * Ur * math.sin(thr) * ecos
    if model == "schamberg":
        al, nu, phi0 = par["energy_accommodation"], par["nu"], math.radians(par["beam_half_angle_deg"])
        vt, vb, rho, wt = tangent_polar_nodes(Ut, c)
        # all combinations of the normal nodes and the tangential nodes
        V2 = w[:, None] ** 2 + rho[None, :] ** 2
        Vr = np.sqrt((1 - al) * V2 + (al * 4 * kB * Tw / mass if al > 0 else 0.0))
        sin_i = np.minimum(rho[None, :] / np.sqrt(V2), 1.0)
        tr = np.arcsin(np.minimum(sin_i ** nu if math.isfinite(nu) else np.where(sin_i < 1, 0.0, 1.0), 1.0))
        tgrid, Mn_tab, Mt_tab = beam_table
        Mn, Mt = np.interp(tr, tgrid, Mn_tab), np.interp(tr, tgrid, Mt_tab)
        W2 = fw[:, None] * wt[None, :]
        un = np.sum(W2 * Vr * Mn)
        ut = np.sum(W2 * Vr * Mt * (vt / rho)[None, :])
        return P_in + un, S_in - ut
    raise ValueError("unknown model '%s'" % model)


def drag_coefficient(p, mass, kB, area, n, A_proj_amps):
    U = p["velocity"]
    speed = np.linalg.norm(U)
    Uh = U / speed
    rho = p["density"] * mass
    model, par, Tw = p["model"], p["par"], p["wall_temperature"]
    beam = None
    if model == "schamberg":
        tgrid = np.linspace(0.0, 0.5 * math.pi, 2001)
        tab = np.array([schamberg_beam_mean(t, math.radians(par["beam_half_angle_deg"])) for t in tgrid])
        beam = (tgrid, tab[:, 0], tab[:, 1])
    # the table of P and S in theta (Gauss-Legendre nodes in theta; also used for the exact sphere)
    th_nodes, th_w = gl(0.0, math.pi, np.polynomial.legendre.leggauss(160))
    PS = np.array([plate(math.cos(t), speed, p["temperature"], mass, kB, Tw, model, par, beam) for t in th_nodes])
    PS *= mass * p["density"]
    # the exact sphere: Cd = 4/(rho U^2) int_0^pi (P cos + S sin) sin dtheta
    cd_sphere = 4.0 / (rho * speed ** 2) * np.sum(th_w * (PS[:, 0] * np.cos(th_nodes) + PS[:, 1] * np.sin(th_nodes)) * np.sin(th_nodes))
    # the faces: theta of each face, P and S interpolated in theta (cubic, on the Gauss nodes)
    cos_f = np.clip(-n @ Uh, -1.0, 1.0)
    th_f = np.arccos(cos_f)
    from numpy.polynomial import legendre as L
    x = th_nodes / (0.5 * math.pi) - 1.0
    cP = L.legfit(x, PS[:, 0], 100)
    cS = L.legfit(x, PS[:, 1], 100)
    xf = th_f / (0.5 * math.pi) - 1.0
    Pf, Sf = L.legval(xf, cP), L.legval(xf, cS)
    F_D = np.sum(area * (Pf * cos_f + Sf * np.sin(th_f)))
    A_exact = np.sum(area * np.maximum(cos_f, 0.0))
    return {
        "cd_faces_amps_area": 2.0 * F_D / (rho * speed ** 2 * A_proj_amps),
        "cd_faces_exact_area": 2.0 * F_D / (rho * speed ** 2 * A_exact),
        "cd_sphere": cd_sphere,
        "A_exact": A_exact,
        "fit_residual": max(np.max(np.abs(L.legval(x, cP) - PS[:, 0])) / np.max(np.abs(PS[:, 0])),
                            np.max(np.abs(L.legval(x, cS) - PS[:, 1])) / max(np.max(np.abs(PS[:, 1])), 1e-300)),
    }


# ------------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--in", dest="inp", required=True)
    ap.add_argument("--conf", required=True)
    ap.add_argument("--mesh", required=True)
    ap.add_argument("--result", action="append", required=True, help="CASE=PATH of drag-coefficient-test.dat")
    ap.add_argument("--sigmas", type=float, default=4.0)
    ap.add_argument("--report")
    args = ap.parse_args()

    out = []
    failed = 0
    xyz, tri = read_nastran(args.mesh, read_conf_scale(args.conf))
    area, n = faces(xyz, tri)
    out.append("INFO body: %d faces, area %.6e m^2" % (len(area), area.sum()))

    for item in args.result:
        case, path = item.split("=", 1)
        try:
            r = read_result(path)
        except OSError as e:
            out.append("FAIL %s result missing: %s" % (case, e))
            failed += 1
            continue
        species = [k[len("model["):-1] for k in r if k.startswith("model[")]
        if len(species) != 1:
            out.append("FAIL %s: one species expected in %s" % (case, path))
            failed += 1
            continue
        sp = species[0]
        p = case_parameters(args.inp, case, sp)
        mass, kB = float(r["mass_kg[%s]" % sp]), float(r["boltzmann_constant"])

        # the input that AMPS used must be the input of this case
        echo = [
            (r["case"] == case, "case"),
            (r["model[%s]" % sp] == p["model"], "model"),
            (np.allclose([float(x) for x in r["upstream_velocity_m_s"].split()], p["velocity"], rtol=1e-8, atol=1e-8), "velocity"),
            (math.isclose(float(r["upstream_temperature_k"]), p["temperature"], rel_tol=1e-8), "temperature"),
            (math.isclose(float(r["upstream_number_density_m-3[%s]" % sp]), p["density"], rel_tol=1e-8), "density"),
        ]
        if p["wall_temperature"] is not None:
            echo.append((math.isclose(float(r["wall_temperature_k"]), p["wall_temperature"], rel_tol=1e-8), "wall temperature"))
        # the parameters of the model, as AMPS installed them ("species.SPECIES.KEY=VALUE", the units of the input)
        amps_par = {k.split(".", 2)[2]: float(v) for k, v in r.items() if k.startswith("species.%s." % sp)}
        echo.append((set(amps_par) == set(p["par"]), "parameter names"))
        for k in p["par"]:
            if k in amps_par:
                echo.append((math.isclose(amps_par[k], p["par"][k], rel_tol=1e-8, abs_tol=1e-12), "parameter " + k))
        bad = [name for ok, name in echo if not ok]
        if bad:
            out.append("FAIL %s: AMPS ran another input (%s differ)" % (case, ", ".join(bad)))
            failed += 1
            continue

        cd, sigma, A = float(r["Cd"]), float(r["sigma_Cd"]), float(r["projection_area_m2"])
        ref = drag_coefficient(p, mass, kB, area, n, A)
        dev = (cd - ref["cd_faces_amps_area"]) / sigma if sigma > 0 else float("inf")
        ok = (sigma > 0) and (sigma < 0.02 * cd) and abs(dev) <= args.sigmas
        failed += 0 if ok else 1
        out.append("%s %s Cd=%.6f reference=%.6f sigma=%.6f deviation=%+.2f sigma (limit %g) hits=%s" % (
            "PASS" if ok else "FAIL", case, cd, ref["cd_faces_amps_area"], sigma, dev, args.sigmas, r.get("hits", "?")))
        out.append("INFO %s exact-sphere Cd=%.6f, faces with their exact projected area Cd=%.6f; projected area: AMPS %.6e, "
                   "faces %.6e m^2 (ratio %.6f); fit residual %.1e" % (
                       case, ref["cd_sphere"], ref["cd_faces_exact_area"], A, ref["A_exact"], A / ref["A_exact"],
                       ref["fit_residual"]))

    text = "\n".join(out) + "\n"
    sys.stdout.write(text)
    if args.report:
        with open(args.report, "w") as f:
            f.write(text)
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
