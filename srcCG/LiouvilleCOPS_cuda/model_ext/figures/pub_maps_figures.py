#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/pub_maps_figures.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/pub_maps_figures.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 74d100307a9e508af0aab99ec6809c699c51ffc4415d99df3df22002419adaea
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure script (full set 3.4: activity / production-rate template and ground-track figures with fixed ranges)
#   and the helper module of the version-2 scripts (fields, tracks, captions); 2026-11-13: MODEL_NAME, MULT_NOTE and
#   the CO2 sentence of the captions data-driven.
# RUN: jobs/pub_maps_figures.sbatch (GPU); imported by pub_maps_v2_figures.py and pub_source_area_figures.py
# ==== end of header ====
"""Publication template and ground-track figures of ONE fitted state with the surface coloured by its H2O activity
or by its production rate (H2O, CO2, total) instead of the late-visibility (NG-only / RG-only) facets.

User request 2026-11-12: "generate the full set of plots that you have done before (templates, gravity, ground
tracks, model/data comparison, etc.); i need three sets of the template and groundtrack plots: 1) as now, 2) surface
activity map instead of visibility maps, 3) production rate [m^-2 s^-1] (H2O, CO2, total) instead of visibility
maps". User answers (same day): the fit = the softplus model after 5 cycles (final state); production-rate times
10:33:10 and 10:39:26 UTC on 2016-09-30; colour scale LOGARITHMIC (log10 over the non-zero range, zero in a separate
grey colour).

Set 1 (visibility) is made by the UNCHANGED tools/pub_figures.py (templates, groundtrack). This script makes sets 2
and 3 with the same projections, panel sizes, impact box, ground track and impact marker (the same
cops_cuda.figures.groundtrack_data of the same base and bundle), and writes the captions of the set-1 ground-track
figures (pub_figures.py writes a caption for the templates figure only).

CODE REUSE. cops_cuda.figures (AMPS commit bad75ec; figures.py last changed in 47fccd5, sha256 e40a66ef...) is
imported unchanged: style, save, read_table, late_window, bundle_metadata, groundtrack_data, gravity_slope,
face_colors, lonlat, OrthoView, render_ortho, render_lonlat, _minute_marks, _scale_bar, write_caption, the colour
constants. Two functions hard-wire the visibility colouring and are COPIED here (marked "copy of ..."), with only the
row / surface colour content changed:
  values_templates_figure     <- cops_cuda.figures.templates_figure (47fccd5): the 'late' row is replaced by the
                                 value rows; header text, legend surface entries, value colour bars and caption.
  values_groundtrack_figures  <- cops_cuda.figures.groundtrack_figures (47fccd5): the surface colours
                                 (face_colors(len(F), masks) -> the value colours); the legend's NG-only / RG-only
                                 entries -> the zero (and below-range) entries; a second colour bar for the surface
                                 values (the altitude colour bar moved below it); the track table is not rewritten
                                 (it is the set-1 one).
  The zoom boxes (templates: late track + NG-only/RG-only facets; 3-D: the same masks) are computed exactly as in the
  originals, from the same bundle masks, so the three sets show the same areas.

FIELDS (one final state; the law of the maps job, model_ext nonneg_maps.py, imported from the frozen copy that made
the maps and sha256-checked against their maps_manifest.json):
  A       = <ft-dir>/finetune_state.npz A_final (the model's H2O activity factor)
  F_H2O   = q_H2O(cos SZA, night/cast shadow, r_h) A,  F_CO2 = kappa q_CO2 A_CO2,  total = F_H2O + F_CO2  [m^-2 s^-1]
            nonneg_maps.source_rates with the maps job's geometry (cos SZA and cast shadow of every facet at t11 =
            10:33:10 and t12 = 10:39:26 UTC, r_h), elements.bin ACO2, fit_restart.npz kappa, the base's
            HeliocentricModel / SZAMode; law_guard checks the reference fitter's expressions.
  VERIFIED bitwise against the maps' .vtu files for ALL facets (A_H2O, F_H2O_per_m2_s, F_CO2_per_m2_s, cos_SZA,
  illumination) and a printed sample of facets; Q_s = sum F area against maps_manifest.json.

TRACK. --track-base (default --base): the base whose table gives the ground track, late window and impact point (user
decision 2026-11-12: the observation-only base runs/shap8_012m_flat_earlier/base, as in the earlier late-visibility
figures, for all three sets); its elements.bin must equal the fit base's (checked: centroid, normal, area, ACO2).

COLOUR. Values > 0: log10(value) on a linear scale (matplotlib YlOrRd, its near-white first 12 % left out), drawn
without relief shading so that the colour matches the colour bar (the colours are divided by the renderer's own
relief intensity, which render_lonlat / render_ortho multiply back; intensity formulas copied from them).
Value = 0 exactly: grey with relief shading (cops_cuda.figures.SURFACE_COLOR, as the background of the set-1 maps).
--range KEY=LO,HI for A, F_H2O, F_CO2 and total (mode fixed; user decision 2026-11-12: A 13,23, F_H2O 12,22, F_CO2
-28,-18, total 12,22): non-zero values below the range green (relief-shaded), above it the top colour (arrow).
Without --range (mode nonzero, the first request): the colour range of a quantity is [min, max] of its non-zero values (over
all facets; production rates over both times). --range-mode span (a proposal, only with --preview-dir): top = the
decade above the --top-percentile of the non-zero values, bottom = top - --span-decades; non-zero values below the
range get a separate colour, values above it the top colour.

Run from <AMPS>/srcCG/LiouvilleCOPS_cuda with .venv/bin/python on a GPU node (Slurm job; jobs/pub_maps_figures.sbatch):

    pub_maps_figures.py --base <run>/base --bundle <run>/bundle/restart_templates.npz \
        --ft-dir <chain>/c5/ft/final --fit-dir <chain>/c5/fit/final --maps-dir <chain>/maps \
        --geometry <run>/maps_geometry/job_<J>/geometry --maps-code <frozen dir>/nonneg_maps.py \
        --out-dir <figure dir> [--visibility-prefix <figure dir>/visibility/12M_groundtrack_visibility] \
        [--preview-dir <figure dir>/colour_range_proposal --span-decades 10 --top-percentile 99.99] [--device 0]

Writes (new files only; the work goes to <out-dir>/pub_maps.incomplete_<UTC>, moved into place at the end):
  activity/12M_templates_activity, activity/12M_groundtrack_activity_{3d,latlon}          (.eps, .png, _caption.txt)
  production_rate/12M_templates_production_rate_<hhmmss>UTC,
  production_rate/12M_groundtrack_production_rate_<hhmmss>UTC_{3d,latlon}                  (.eps, .png, _caption.txt)
  checks/fields_verification.json, checks/value_statistics.txt, checks/pub_maps_manifest.json
  <visibility-prefix>_{3d,latlon}_caption.txt                    (captions of the set-1 ground-track figures)
  <preview-dir>/...  the template figures with the span colour range (a proposal for the user)
"""

import argparse
import datetime as dt
import hashlib
import importlib.util
import json
import os
import socket
import sys
import time
from pathlib import Path

import numpy as np

APP = Path(os.environ.get("APP", "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda"))
if str(APP) not in sys.path:
    sys.path.insert(0, str(APP))
from cops_cuda import figures as FIG  # noqa: E402

# Name of the fitted state in the captions and the bound note of the slope multipliers; set by the calling figure script
# (pub_maps_v2_figures.py --model-name, and its data-driven multiplier note, 2026-11-13 for the refit with multipliers
# >= 1); the defaults keep the captions of the softplus 5-cycle fit of 2026-11-12 unchanged.
MODEL_NAME = "the softplus activity model, final state after 5 fit and fine-tuning cycles"
MULT_NOTE = ""

SCHEMA = "pub_maps_figures_v1"
TIMES = (("t11", "2016-09-30T10:33:10"), ("t12", "2016-09-30T10:39:26"))
CMAP_NAME, CMAP_FROM = "YlOrRd", 0.12
BELOW_COLOR = np.array([0x4D, 0xAF, 0x4A]) / 255.0       # span mode only: 0 < value < bottom of the range (green:
                                                            # apart from grey, YlOrRd and the cliff blues)
ZERO_PATCH = FIG.SURFACE_COLOR * 0.85                       # legend swatch of the relief-shaded zero grey
REL_TOL = 1e-14
T0 = time.time()


def say(text):
    print("pub_maps_figures [%6.0f s]: %s" % (time.time() - T0, text), flush=True)


def sha256_file(path, bufsize=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(bufsize)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def fmt_e(x, d=2):
    return ("%%.%de" % d) % x


def pct(x):
    return "%.1f %%" % (100.0 * x)


def thousands(n):
    return "{:,}".format(int(n))


# ============================================================================================ fields
def load_fields(args):
    """A, F_H2O / F_CO2 / total at both times, seen, bins, and the checks against the maps' files."""
    maps = Path(args.maps_dir)
    man = json.load(open(maps / "maps_manifest.json"))
    code_sha = sha256_file(args.maps_code)
    if code_sha != man["code_sha256"]["nonneg_maps.py"]:
        raise SystemExit("%s sha256 %s is not the maps' nonneg_maps.py (%s)" % (
            args.maps_code, code_sha, man["code_sha256"]["nonneg_maps.py"]))
    for key, d in (("ft_dir_resolved", args.ft_dir), ("fit_dir_resolved", args.fit_dir)):
        if Path(man["state"][key]) != Path(d).resolve():
            raise SystemExit("maps %s were made from %s, not %s" % (maps, man["state"][key], Path(d).resolve()))
    if Path(man["base"]) != Path(args.base).resolve():
        raise SystemExit("maps made on base %s, not %s" % (man["base"], args.base))
    if Path(man["geometry"]) != Path(args.geometry).resolve():
        raise SystemExit("maps made with geometry %s, not %s" % (man["geometry"], args.geometry))
    nm = load_module(args.maps_code, "nonneg_maps_frozen")
    fc, header, _, _ = nm.base_table(args.base)
    law = nm.law_guard(fc)
    helio, sza = fc.header_value(header, "HeliocentricModel"), fc.header_value(header, "SZAMode")
    el = nm.read_elements(Path(args.base) / "elements.bin")
    nf = el["nfaces"]
    G = nm.read_geometry(args.geometry)
    if not (np.array_equal(G["facegeom"][:, 6], el["area"]) and np.array_equal(G["facegeom"][:, 0:3], el["centroid"])):
        raise SystemExit("geometry mesh is not the base's (area / centroid)")
    del G["facegeom"], G["nodes"], G["faces"]
    with np.load(Path(args.ft_dir) / "finetune_state.npz", allow_pickle=False) as z:
        A = np.asarray(z["A_final"], np.float64)
    with np.load(Path(args.fit_dir) / "fit_restart.npz", allow_pickle=False) as z:
        kappa = float(np.asarray(z["kappa"]).reshape(-1)[0])
        tv = np.asarray(z["template_values"], np.float64)
    with np.load(Path(args.fit_dir) / "nonneg_facets.npz", allow_pickle=False) as z:
        seen = np.asarray(z["seen"]).astype(bool)
        bins = np.asarray(z["bin"]).astype(np.intp)
        meta = json.loads(str(z["meta_json"]))
    edges = meta["model"]["bin_edges_deg"]
    mult = tv[:len(edges) - 1]
    if len(A) != nf or len(seen) != nf:
        raise SystemExit("facet counts differ: A %d, seen %d, base %d" % (len(A), len(seen), nf))
    chk = {"maps_code_sha256": code_sha, "law_guard": law["checked"], "heliocentric_model": helio, "sza_mode": sza,
           "kappa": kappa, "kappa_manifest": man["state_checks"]["kappa"]}
    chk["kappa_equal_manifest"] = kappa == man["state_checks"]["kappa"]
    say("verify A against %s (all facets)" % (maps / "softplus_static.vtu"))
    rd, fd, _ = nm.read_vtu(maps / ("%s_static.vtu" % man["model"]))
    chk["A_equal_vtu_all_facets"] = bool(np.array_equal(rd["A_H2O"], A))
    chk["seen_equal_vtu"] = bool(np.array_equal(rd["seen"].astype(bool), seen))
    chk["kappa_ACO2_equal_vtu"] = bool(np.array_equal(rd["CO2_activity_kappa_ACO2"], kappa * el["ACO2"]))
    chk["kappa_equal_vtu_fielddata"] = fd.get("kappa", [None])[0] == kappa
    rng = np.random.default_rng(20261112)
    sample = np.unique(np.concatenate((rng.choice(nf, 8, replace=False), rng.choice(np.nonzero(A > 0)[0], 8, replace=False),
                                       [int(np.argmax(A))])))
    del rd
    fields = {"A": A, "seen": seen, "bins": bins, "mult": mult, "edges": edges, "kappa": kappa, "area": el["area"],
              "ACO2": el["ACO2"], "normal": el["normal"], "times": {}}
    chk["times"] = {}
    for lab, utc in TIMES:
        T = G["times"][lab]
        if T["utc"] != utc:
            raise SystemExit("geometry %s is %s, expected %s" % (lab, T["utc"], utc))
        cos, sh = nm.geometry_time(args.geometry, lab, nf)
        Fh, Fc, night, qf = nm.source_rates(cos, sh, T["rh_m"], helio, sza, A, el["ACO2"], kappa, fc.AU)
        ill = nm.illumination(cos, sh)
        tot = Fh + Fc
        stamp = utc.replace("-", "").replace(":", "")
        vtu = maps / ("%s_flux_%s.vtu" % (man["model"], stamp))
        say("verify F_H2O, F_CO2, cos SZA, illumination at %s against %s (all facets)" % (utc, vtu.name))
        rd, fd, _ = nm.read_vtu(vtu)
        c = {"vtu": str(vtu), "F_H2O_equal_vtu_all_facets": bool(np.array_equal(rd["F_H2O_per_m2_s"], Fh)),
             "F_CO2_equal_vtu_all_facets": bool(np.array_equal(rd["F_CO2_per_m2_s"], Fc)),
             "F_H2O_max_rel_diff_vtu": max_rel(rd["F_H2O_per_m2_s"], Fh),
             "F_CO2_max_rel_diff_vtu": max_rel(rd["F_CO2_per_m2_s"], Fc),
             "cos_SZA_equal_vtu": bool(np.array_equal(rd["cos_SZA"], cos)),
             "illumination_equal_vtu": bool(np.array_equal(rd["illumination"], ill)),
             "TimeValue_vtu": fd["TimeValue"][0], "rh_m_vtu": fd["rh_m"][0], "rh_m": T["rh_m"]}
        c["sample"] = [{"facet": int(e), "A": float(A[e]), "cos_SZA": float(cos[e]), "illumination": int(ill[e]),
                        "F_H2O": float(Fh[e]), "F_H2O_vtu": float(rd["F_H2O_per_m2_s"][e]), "F_CO2": float(Fc[e]),
                        "F_CO2_vtu": float(rd["F_CO2_per_m2_s"][e]), "total": float(tot[e])} for e in sample]
        mt = [t for t in man["times"] if t["label"] == lab][0]
        Qh, Qc = float(np.sum(Fh * el["area"])), float(np.sum(Fc * el["area"]))
        c.update({"Q_H2O_s-1": Qh, "Q_CO2_s-1": Qc, "Q_H2O_manifest": mt["Q_H2O_s-1"], "Q_CO2_manifest": mt["Q_CO2_s-1"],
                  "Q_equal_manifest": Qh == mt["Q_H2O_s-1"] and Qc == mt["Q_CO2_s-1"],
                  "Q_rel_diff_manifest": max(abs(Qh / mt["Q_H2O_s-1"] - 1), abs(Qc / mt["Q_CO2_s-1"] - 1))})
        chk["times"][lab] = c
        del rd
        cnt = {name: int(np.sum(ill == k)) for k, name, _ in nm.ILLUMINATION}
        ar = {name: float(np.sum(el["area"][ill == k])) for k, name, _ in nm.ILLUMINATION}
        fields["times"][lab] = {"utc": utc, "F_H2O": Fh, "F_CO2": Fc, "total": tot, "illumination": ill,
                                "rh_AU": qf["rhA"], "q": {"H2O": qf["H2O"], "CO2": qf["CO2"]}, "p": qf["p"],
                                "Q_H2O": Qh, "Q_CO2": Qc, "facets": cnt, "area": ar,
                                "sun": [float(v) for v in T["sun"]]}
        say("%s: Q_H2O %.6e (manifest %.6e), Q_CO2 %.6e; F_H2O %s, F_CO2 %s bitwise = vtu" % (
            utc, Qh, mt["Q_H2O_s-1"], Qc, c["F_H2O_equal_vtu_all_facets"], c["F_CO2_equal_vtu_all_facets"]))
    # A, cos, illumination: bitwise (stored / read values); F: bitwise expected (elementwise IEEE arithmetic), a
    # last-bit difference (another CPU's libm pow for q) is accepted up to REL_TOL and reported; Q likewise
    chk["bitwise_all"] = bool(all(v["F_H2O_equal_vtu_all_facets"] and v["F_CO2_equal_vtu_all_facets"] and
                                  v["Q_equal_manifest"] for v in chk["times"].values()))
    chk["passed"] = bool(chk["kappa_equal_manifest"] and chk["A_equal_vtu_all_facets"] and chk["seen_equal_vtu"] and
                         chk["kappa_ACO2_equal_vtu"] and all(
                             v["F_H2O_max_rel_diff_vtu"] <= REL_TOL and v["F_CO2_max_rel_diff_vtu"] <= REL_TOL and
                             v["cos_SZA_equal_vtu"] and v["illumination_equal_vtu"] and
                             v["Q_rel_diff_manifest"] <= REL_TOL for v in chk["times"].values()))
    fields["checks"] = chk
    fields["maps_manifest"] = man
    return fields


def max_rel(a, b):
    """max |a - b| / max(|a|, |b|) over the facets (0 where both are 0)."""
    d = np.abs(a - b)
    m = np.maximum(np.abs(a), np.abs(b))
    with np.errstate(divide="ignore", invalid="ignore"):
        r = np.where(m > 0, d / m, 0.0)
    return float(r.max())


def value_stats(v, seen, area):
    nz = v > 0
    out = {"n": int(v.size), "zero": int(np.sum(~nz)), "zero_frac": float(np.mean(~nz)),
           "zero_area_frac": float(np.sum(area[~nz]) / np.sum(area)),
           "zero_frac_seen": float(np.mean(~nz[seen])),
           "zero_area_frac_seen": float(np.sum(area[~nz & seen]) / np.sum(area[seen])),
           "nonzero": int(np.sum(nz)), "nonzero_unseen": int(np.sum(nz & ~seen))}
    if nz.any():
        x = v[nz]
        out.update({"min_nonzero": float(x.min()), "max": float(x.max()),
                    "percentiles_nonzero": {str(p): float(np.percentile(x, p)) for p in
                                            (0.1, 1, 5, 10, 25, 50, 75, 90, 99, 99.9, 99.99)}})
        lg = np.log10(x)
        e = np.arange(np.floor(lg.min() / 10) * 10, np.ceil(lg.max() / 10) * 10 + 10, 10)
        h, e = np.histogram(lg, bins=e)
        out["log10_histogram_10dec"] = [[float(a), int(b)] for a, b in zip(e[:-1], h) if b]
    return out


def colour_range(values_list, mode, span_decades, top_percentile):
    x = np.concatenate([v[v > 0] for v in values_list])
    if mode == "nonzero":
        return float(np.log10(x.min())), float(np.log10(x.max()))
    hi = float(np.ceil(np.log10(np.percentile(x, top_percentile))))
    return hi - span_decades, hi


# ============================================================================================ colours
def value_cmap(plt):
    from matplotlib.colors import ListedColormap
    return ListedColormap(plt.get_cmap(CMAP_NAME)(np.linspace(CMAP_FROM, 1.0, 256)), name=CMAP_NAME + "_from012")


def value_rgb(values, lo, hi, cmap, below=False):
    """Per-facet RGB: grey for value 0, cmap(log10 value) otherwise; returns (rgb, coloured mask)."""
    rgb = np.tile(FIG.SURFACE_COLOR, (len(values), 1))
    nz = values > 0
    lg = np.log10(values[nz])
    rgb[nz] = cmap(np.clip((lg - lo) / (hi - lo), 0.0, 1.0))[:, :3]
    if below:
        b = np.zeros(len(values), bool)
        b[nz] = lg < lo
        rgb[b] = BELOW_COLOR
    return rgb, nz


def lonlat_intensity(V, F, N):
    """The relief intensity of cops_cuda.figures.render_lonlat at shade_weight 1 (its light, copied)."""
    c = V[F].mean(axis=1)
    up = c / np.linalg.norm(c, axis=1)[:, None]
    east = np.cross([0.0, 0.0, 1.0], up)
    east /= np.maximum(np.linalg.norm(east, axis=1), 1e-12)[:, None]
    north = np.cross(up, east)
    L = np.cos(np.radians(40)) * (north - east) / np.sqrt(2) + np.sin(np.radians(40)) * up
    lam = np.clip(np.einsum("ij,ij->i", N, L), 0.0, 1.0)
    return 1.0 - (0.70 - 0.70 * lam)


def ortho_intensity(N, view):
    """The relief intensity of cops_cuda.figures.render_ortho at shade_weight 1 with its default light (copied);
    returns (intensity, front-facing mask)."""
    L = 0.75 * view.d - 0.35 * view.right + 0.45 * view.up
    L = L / np.linalg.norm(L)
    lam = np.clip(N @ L, 0.0, 1.0)
    facing = N @ view.d
    return np.where(facing < 0, 0.25, 1.0 - (0.70 - 0.70 * lam)), facing >= 0


def flat(rgb, coloured, intensity, front=None):
    """Colours divided by the renderer's relief intensity on the coloured facets, so that they come out unshaded."""
    out = rgb.copy()
    m = coloured if front is None else (coloured & front)
    out[m] = rgb[m] / intensity[m][:, None]
    return out


class ValueLayer:
    """One quantity to colour the surface with."""

    def __init__(self, key, title, values, lo, hi, cb_label, cmap, below_mode):
        self.key, self.title, self.values, self.lo, self.hi, self.cb_label = key, title, values, lo, hi, cb_label
        self.cmap = cmap
        nz = values > 0
        lg = np.log10(values[nz])
        self.n_below = int(np.sum(lg < lo)) if below_mode else 0
        self.n_above = int(np.sum(lg > hi))
        self.below = below_mode and self.n_below > 0
        self.rgb, self.coloured = value_rgb(values, lo, hi, cmap, below=self.below)
        if self.below:                      # the below-range class keeps the relief shading (it is a category)
            b = np.zeros(len(values), bool)
            b[nz] = lg < lo
            self.coloured = self.coloured & ~b

    def mappable(self):
        from matplotlib.cm import ScalarMappable
        from matplotlib.colors import Normalize
        cm = self.cmap.copy()
        cm.set_under(BELOW_COLOR)
        cm.set_over(self.cmap(1.0))
        extend = {(False, False): "neither", (True, False): "min", (False, True): "max", (True, True): "both"}[
            (self.below, self.n_above > 0)]
        return ScalarMappable(norm=Normalize(self.lo, self.hi), cmap=cm), extend


# ============================================================================================ figures
def zoom_box(gt, late):
    """The impact-area box of cops_cuda.figures.templates_figure (copied): late track + NG-only/RG-only facets."""
    V, F, masks = gt["vertices"], gt["faces"], gt["masks"] or {}
    glon, glat = FIG.lonlat(gt["ground"])
    lons, lats = [glon[late]], [glat[late]]
    for name in ("late_ng_not_rg", "late_rg_not_ng"):
        if name in masks and masks[name].any():
            lo, la = FIG.lonlat(V[F[masks[name]]].mean(axis=1))
            lons.append(lo)
            lats.append(la)
    lons, lats = np.concatenate(lons), np.concatenate(lats)
    pad = 0.10 * max(np.ptp(lons), np.ptp(lats), 1.0)
    return (lons.min() - pad, lons.max() + pad), (lats.min() - pad, lats.max() + pad)


_RENDER_CACHE = {}


def values_templates_figure(gt, prefix, layers_in, header, legend_surface, caption, I_ll, gravity_cache=None,
                            device_id=0, overwrite=False, title=None, pixels=900):
    """Copy of cops_cuda.figures.templates_figure (AMPS 47fccd5) with the late-visibility row replaced by the value
    rows ``layers_in`` (ValueLayer list); header text, legend surface entries, value colour bars and caption are the
    caller's. Rows after the value rows: effective gravitational slope and cliff templates, unchanged. The slope and
    cliff images are cached between figures (identical inputs give identical images)."""
    plt = FIG.style()
    import matplotlib.patheffects as pe
    from matplotlib.collections import LineCollection
    from matplotlib.colors import Normalize
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    V, F, N = gt["vertices"], gt["faces"], gt["normals"]
    masks, window, impact, table = gt["masks"] or {}, gt["window"], gt["impact"], gt["table"]
    meta = FIG.bundle_metadata(gt["bundle"]) if gt["bundle"] else {}
    slope, cache_path = FIG.gravity_slope(gt["bundle"], gravity_cache, len(F)) if gt["bundle"] else (None, None)
    vrows = {L.key: L for L in layers_in}
    rows = [L.key for L in layers_in] + (["slope", "cliff"] if slope is not None else [])
    late = np.isin(table["iPoint"], window["points"])
    glon, glat = FIG.lonlat(gt["ground"])
    step = np.linalg.norm(np.diff(gt["ground"], axis=0), axis=1)
    ok_seg = (np.abs(np.diff(glon)) <= 180.0) & (step <= max(150.0, 5.0 * np.median(step)))
    zl, za = zoom_box(gt, late)
    coslat = np.cos(np.radians(0.5 * (za[0] + za[1])))

    cmap = plt.get_cmap("magma")
    snorm = Normalize(0.0, 90.0)
    layers = {}
    base = FIG.face_colors(len(F), None)
    for L in layers_in:                                         # changed: value rows instead of the 'late' row
        layers[L.key] = (flat(L.rgb, L.coloured, I_ll), 1.0)
    if slope is not None:
        layers["slope"] = (cmap(snorm(np.clip(slope, 0, 90)))[:, :3], 0.0)
        cc = base.copy()
        c45 = masks.get("cliff_ge_45deg", slope >= 45.0)
        c60 = masks.get("cliff_ge_60deg", slope >= 60.0)
        cc[c45] = FIG.CLIFF45_COLOR
        cc[c60] = FIG.CLIFF60_COLOR
        layers["cliff"] = (cc, 1.0)

    nrow = len(rows)
    row_h = 2.05
    height = row_h * nrow + 1.55
    fig = plt.figure(figsize=(FIG.DOUBLE_COLUMN_IN, height))
    gs = fig.add_gridspec(nrow, 3, width_ratios=[2.0, 1.0, 0.05], wspace=0.10, hspace=0.30,
                          left=0.07, right=0.93, top=1.0 - 0.10 / nrow, bottom=1.35 / height)
    halo = [pe.Stroke(linewidth=1.5, foreground="white"), pe.Normal()]
    halo_thick = [pe.Stroke(linewidth=2.9, foreground="white"), pe.Normal()]
    letters = iter("abcdefghijkl")
    row_titles = {"slope": "effective gravitational slope", "cliff": "cliff templates"}
    row_titles.update({L.key: L.title for L in layers_in})
    for r, row in enumerate(rows):
        colors, shade_w = layers[row]
        for c, (lon_r, lat_r) in enumerate((((-180, 180), (-90, 90)), (zl, za))):
            ax = fig.add_subplot(gs[r, c])
            width = 2 * pixels if c == 0 else pixels
            key = (row, c, width, tuple(lon_r), tuple(lat_r))
            if row in ("slope", "cliff") and key in _RENDER_CACHE:
                image = _RENDER_CACHE[key]
            else:
                image, _ = FIG.render_lonlat(V, F, N, lon_r, lat_r, width, colors, device_id, shade_weight=shade_w)
                if row in ("slope", "cliff"):
                    _RENDER_CACHE[key] = image
            ax.imshow(image, extent=(lon_r[0], lon_r[1], lat_r[0], lat_r[1]), interpolation="lanczos",
                      aspect="equal" if c == 0 else 1.0 / coslat)
            pts = np.column_stack((glon, glat))
            seg = np.stack((pts[:-1], pts[1:]), axis=1)
            keep = gt["outside"][:-1] & gt["outside"][1:] & ok_seg
            ax.add_collection(LineCollection(seg[keep], colors="k", linewidths=0.8, zorder=5, path_effects=halo))
            lk = keep & late[:-1] & late[1:]
            ax.add_collection(LineCollection(seg[lk], colors="k", linewidths=2.0, zorder=6, path_effects=halo_thick))
            ax.plot([impact["lon_deg"]], [impact["lat_deg"]], marker="*", ms=8 if c == 0 else 11,
                    mfc="#F0E442", mec="k", mew=0.6, zorder=9)
            if c == 0:
                ax.plot([glon[0]], [glat[0]], marker="s", ms=3.6, mfc="white", mec="k", mew=0.6, zorder=9)
                ax.plot([zl[0], zl[1], zl[1], zl[0], zl[0]], [za[0], za[0], za[1], za[1], za[0]], color="k",
                        lw=0.6, zorder=7, path_effects=halo)
                ax.set_xticks(range(-180, 181, 60))
                ax.set_yticks(range(-90, 91, 45))
                ax.set_ylabel("latitude [°]")
                ax.set_title(row_titles[row], fontsize=8, loc="left", pad=3)
            else:
                FIG._minute_marks(ax, gt, late, lambda p: FIG.lonlat(p), label=False)
                if r == 0:
                    ax.annotate("impact %s UTC" % impact["utc"].split("T")[1],
                                xy=(impact["lon_deg"], impact["lat_deg"]), xytext=(18, -48),
                                textcoords="offset points", fontsize=6.0, ha="center",
                                arrowprops={"arrowstyle": "-|>", "lw": 0.6, "color": "k", "shrinkB": 5},
                                bbox={"facecolor": "white", "edgecolor": "0.3", "lw": 0.4,
                                      "boxstyle": "round,pad=0.2"}, zorder=10)
            ax.set_xlim(*lon_r)
            ax.set_ylim(*lat_r)
            if r == nrow - 1:
                ax.set_xlabel("longitude [°]")
            ax.text(0.015, 0.96, "(%s)" % next(letters), transform=ax.transAxes, va="top", fontsize=7.5,
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.8})
        cax = fig.add_subplot(gs[r, 2])
        if row == "slope":
            from matplotlib.cm import ScalarMappable
            cb = fig.colorbar(ScalarMappable(norm=snorm, cmap=cmap), cax=cax, extend="max")
            cb.set_label("slope [°] (> 90: overhang)", fontsize=6.5)
            cb.set_ticks([0, 15, 30, 45, 60, 75, 90])
        elif row in vrows:                                      # changed: colour bar of a value row
            sm, extend = vrows[row].mappable()
            cb = fig.colorbar(sm, cax=cax, extend=extend)
            cb.set_label(vrows[row].cb_label, fontsize=6.5)
        else:
            cax.axis("off")
    handles = list(legend_surface)                              # changed: zero / below-range instead of NG/RG-only
    if slope is not None:
        for name, col, thr in (("cliff_ge_45deg", FIG.CLIFF45_COLOR, "45–60°"),
                               ("cliff_ge_60deg", FIG.CLIFF60_COLOR, "≥ 60°")):
            frac = meta.get("templates", {}).get(name, {}).get("area_fraction", np.nan)
            handles.append(Patch(facecolor=col, label="slope %s (%s: %.1f%% of area)" % (thr, name, 100 * frac)))
        handles[-2].set_label(handles[-2].get_label().replace("of area)", "of area, incl. ≥ 60°)"))
    handles += [Line2D([], [], color="k", lw=0.8, label="ground track (sub-spacecraft intercept)"),
                Line2D([], [], color="k", lw=2.0, label="ground track in the late-visibility window"),
                Line2D([], [], marker="o", ls="none", ms=3, mfc="white", mec="k", mew=0.5,
                       label="1-min marks from %s UTC" % window["start"].split("T")[1]),
                Line2D([], [], marker="s", ls="none", ms=3.6, mfc="white", mec="k", mew=0.6,
                       label="track start %s UTC, %.1f km" % (table["time"][0].strftime("%H:%M:%S"),
                                                             table["altitude"][0] / 1000)),
                Line2D([], [], marker="*", ls="none", ms=8, mfc="#F0E442", mec="k", mew=0.6,
                       label="Rosetta impact, lon %.2f°, lat %.2f°" % (impact["lon_deg"], impact["lat_deg"]))]
    fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(0.07, 0.0), ncol=2, fontsize=6.3,
               frameon=False, handlelength=2.0, columnspacing=1.2)
    fig.text(0.07, 1.0 + 0.02 / nrow, header, fontsize=6.8, va="bottom", ha="left", linespacing=1.35)
    if title:
        fig.suptitle(title, fontsize=8, y=1.0 + 0.06 / nrow)
    paths = FIG.save(fig, prefix, overwrite)
    plt.close(fig)
    paths.append(FIG.write_caption(prefix, caption, overwrite))
    return paths


def values_groundtrack_figures(gt, prefix, layer, legend_surface, captions, I_ll, device_id=0, overwrite=False,
                               title=None, pixels=800):
    """Copy of cops_cuda.figures.groundtrack_figures (AMPS 47fccd5): ``<prefix>_3d`` and ``<prefix>_latlon`` with
    the surface coloured by ``layer`` (ValueLayer) instead of the NG-only / RG-only facets. Changed: surface colours;
    legend surface entries; a second colour bar (surface values) above the altitude colour bar; captions written;
    the track table is not rewritten (it is the visibility set's)."""
    plt = FIG.style()
    from matplotlib.colors import LogNorm
    from matplotlib.lines import Line2D
    table, window, impact = gt["table"], gt["window"], gt["impact"]
    V, F, N = gt["vertices"], gt["faces"], gt["normals"]
    alt = np.maximum(table["altitude"], 0.1)
    norm = LogNorm(vmin=max(alt.min(), 0.1), vmax=alt.max())
    cmap = plt.get_cmap("viridis")
    late = np.zeros(len(alt), bool)
    late[[int(np.nonzero(table["iPoint"] == p)[0][0]) for p in window["points"]]] = True
    ip = np.asarray(impact["point_m"])
    ground = gt["ground"]
    step = np.linalg.norm(np.diff(ground, axis=0), axis=1)
    no_jump = step <= max(150.0, 5.0 * np.median(step))
    written = []
    legend_extra = [Line2D([], [], marker="*", ls="none", ms=8, mfc="#F0E442", mec="k", mew=0.6,
                           label="Rosetta impact (lon %.2f°, lat %.2f°)" % (impact["lon_deg"], impact["lat_deg"])),
                    Line2D([], [], marker="s", ls="none", ms=4, mfc="white", mec="k", mew=0.6,
                           label="track start, %s UTC, %.1f km" % (table["time"][0].strftime("%H:%M:%S"),
                                                                    table["altitude"][0] / 1000.0))]

    def region_legend(anchor, ncol, hidden):
        # copy of cops_cuda.figures._region_legend (47fccd5); changed: surface entries instead of NG-only / RG-only
        handles = list(legend_surface)
        handles.append(Line2D([], [], color="#21918C", lw=2.6, label="ground track, late-visibility window"))
        handles.append(Line2D([], [], color="#21918C", lw=1.1, label="ground track (sub-spacecraft intercept)"))
        if hidden:
            handles.append(Line2D([], [], color="0.45", lw=0.6, ls=(0, (1, 1.5)), label="ground track, hidden"))
        handles.append(Line2D([], [], marker="o", ls="none", ms=3, mfc="white", mec="k", mew=0.5,
                              label="1-min marks from %s" % "window start"))
        handles.extend(legend_extra)
        fig.legend(handles=handles, loc="lower left", bbox_to_anchor=anchor, fontsize=6.3, ncol=ncol, frameon=False,
                   handlelength=2.2, columnspacing=1.0)

    # ---------------- 3D orthographic views
    approach = np.column_stack((table["x"], table["y"], table["z"]))[late][0] - ip
    view = FIG.OrthoView(approach)
    I3, front = ortho_intensity(N, view)
    colors = flat(layer.rgb, layer.coloured, I3, front)        # changed: value colours
    body_half = 1.04 * np.max(np.abs(np.column_stack(view.project(V - V.mean(axis=0))[:2])))
    zoom_half = 1.15 * max(np.max(np.linalg.norm(ground[late] - ip, axis=1)), 1.0)
    if gt["masks"] is not None:
        for name in ("late_ng_not_rg", "late_rg_not_ng"):
            if name in gt["masks"] and gt["masks"][name].any():
                cen = V[F[gt["masks"][name]]].mean(axis=1)
                zoom_half = max(zoom_half, 1.08 * np.max(np.abs(np.column_stack(view.project(cen - ip)[:2]))))
    close_half = 1.3 * max(np.max(np.abs(np.column_stack(view.project(ground[late] - ip)[:2]))), 20.0)
    panels = ((V.mean(axis=0), body_half), (ip, zoom_half), (ip, close_half))
    fig = plt.figure(figsize=(FIG.DOUBLE_COLUMN_IN, 3.35))
    width = 0.315
    axs = [fig.add_axes([0.005 + k * (width + 0.0175), 0.27, width, width * FIG.DOUBLE_COLUMN_IN / 3.35])
           for k in range(3)]
    lc = None
    for k, (ax, (center, half)) in enumerate(zip(axs, panels)):
        image, to_image, visible = FIG.render_ortho(V, F, N, view, center, half, pixels, colors, device_id)
        ax.imshow(image, interpolation="lanczos", origin="upper", extent=(0, pixels, pixels, 0))
        gx, gy, _ = to_image(ground)
        tol = max(3.0, 4.0 * 2.0 * half / pixels)
        vis = visible(ground, tol) & gt["outside"]
        hidden = ~vis & gt["outside"]
        FIG._colored_track(ax, gx, gy, alt, norm, cmap, 0.6, mask=hidden, seg_mask=no_jump, zorder=4,
                           solid_color="0.45", dotted=True)
        FIG._colored_track(ax, gx, gy, alt, norm, cmap, 2.6, mask=vis & late, seg_mask=no_jump, zorder=5)
        lc = FIG._colored_track(ax, gx, gy, alt, norm, cmap, 1.1, mask=vis, seg_mask=no_jump, zorder=6)
        if vis[0]:
            ax.plot([gx[0]], [gy[0]], marker="s", ms=4, mfc="white", mec="k", mew=0.6, zorder=9)
        ix, iy, _ = to_image(ip[None, :])
        ax.plot([ix[0]], [iy[0]], marker="*", ms=9 if k < 2 else 11, mfc="#F0E442", mec="k", mew=0.6, zorder=9)
        if k < 2:
            s = panels[k + 1][1] / half * pixels / 2
            ax.plot([ix[0] - s, ix[0] + s, ix[0] + s, ix[0] - s, ix[0] - s],
                    [iy[0] - s, iy[0] - s, iy[0] + s, iy[0] + s, iy[0] - s], color="k", lw=0.6, zorder=7)
        if k == 2:
            FIG._minute_marks(ax, gt, late, lambda p: to_image(p)[:2], label=True)
            ax.annotate("Rosetta impact\n%s UTC" % impact["utc"].split("T")[1], xy=(ix[0], iy[0]),
                        xytext=(-45, -55), textcoords="offset points", fontsize=6.3, ha="center",
                        arrowprops={"arrowstyle": "-|>", "lw": 0.7, "color": "k", "shrinkB": 6},
                        bbox={"facecolor": "white", "edgecolor": "0.3", "lw": 0.5, "boxstyle": "round,pad=0.25"},
                        zorder=10)
        else:
            FIG._minute_marks(ax, gt, late, lambda p: to_image(p)[:2], label=False)
        ax.set_xlim(0, pixels)
        ax.set_ylim(pixels, 0)
        ax.set_xticks([])
        ax.set_yticks([])
        FIG._scale_bar(ax, half, pixels)
        ax.text(0.03, 0.97, "(%s)" % "abc"[k], transform=ax.transAxes, va="top", fontsize=8,
                bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
    region_legend((0.005, 0.0), 2, True)
    sm, extend = layer.mappable()                               # changed: surface colour bar, altitude bar below it
    cb = fig.colorbar(sm, cax=fig.add_axes([0.70, 0.165, 0.29, 0.035]), orientation="horizontal", extend=extend)
    cb.set_label(layer.cb_label, fontsize=7)
    cb = fig.colorbar(lc, cax=fig.add_axes([0.70, 0.015, 0.29, 0.035]), orientation="horizontal")
    cb.set_label("spacecraft altitude [m]", fontsize=7)
    if title:
        fig.suptitle(title, fontsize=8, y=1.0)
    written += FIG.save(fig, str(prefix) + "_3d", overwrite)
    plt.close(fig)
    written.append(FIG.write_caption(str(prefix) + "_3d", captions["3d"], overwrite))

    # ---------------- longitude/latitude map
    colors = flat(layer.rgb, layer.coloured, I_ll)              # changed: value colours
    glon, glat = FIG.lonlat(ground)
    no_seam = (np.abs(np.diff(glon)) <= 180.0) & no_jump
    lons, lats = [glon[late]], [glat[late]]
    if gt["masks"] is not None:
        for name in ("late_ng_not_rg", "late_rg_not_ng"):
            if name in gt["masks"] and gt["masks"][name].any():
                lo, la = FIG.lonlat(V[F[gt["masks"][name]]].mean(axis=1))
                lons.append(lo)
                lats.append(la)
    lons, lats = np.concatenate(lons), np.concatenate(lats)
    pad = 0.10 * max(np.ptp(lons), np.ptp(lats), 1.0)
    zl = (lons.min() - pad, lons.max() + pad)
    za = (lats.min() - pad, lats.max() + pad)
    coslat = np.cos(np.radians(0.5 * (za[0] + za[1])))
    fig = plt.figure(figsize=(FIG.DOUBLE_COLUMN_IN, 7.0))
    ax = fig.add_axes([0.08, 0.60, 0.90, 0.38])
    image, _ = FIG.render_lonlat(V, F, N, (-180, 180), (-90, 90), 2 * pixels, colors, device_id)
    ax.imshow(image, extent=(-180, 180, -90, 90), interpolation="lanczos", aspect="equal")
    FIG._colored_track(ax, glon, glat, alt, norm, cmap, 2.6, mask=late & gt["outside"], seg_mask=no_seam, zorder=5)
    lc = FIG._colored_track(ax, glon, glat, alt, norm, cmap, 1.0, mask=gt["outside"], seg_mask=no_seam, zorder=6)
    ax.plot([impact["lon_deg"]], [impact["lat_deg"]], marker="*", ms=9, mfc="#F0E442", mec="k", mew=0.6, zorder=9)
    ax.plot([glon[0]], [glat[0]], marker="s", ms=4, mfc="white", mec="k", mew=0.6, zorder=9)
    ax.plot([zl[0], zl[1], zl[1], zl[0], zl[0]], [za[0], za[0], za[1], za[1], za[0]], color="k", lw=0.6, zorder=7)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xticks(range(-180, 181, 60))
    ax.set_yticks(range(-90, 91, 30))
    ax.set_xlabel("longitude [°]")
    ax.set_ylabel("latitude [°]")
    ax.text(0.01, 0.97, "(a)", transform=ax.transAxes, va="top", fontsize=8,
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
    axz = fig.add_axes([0.08, 0.13, 0.90, 0.40])
    image, _ = FIG.render_lonlat(V, F, N, zl, za, 2 * pixels, colors, device_id)
    axz.imshow(image, extent=(zl[0], zl[1], za[0], za[1]), interpolation="lanczos", aspect=1.0 / coslat)
    FIG._colored_track(axz, glon, glat, alt, norm, cmap, 2.6, mask=late, seg_mask=no_seam, zorder=5)
    FIG._colored_track(axz, glon, glat, alt, norm, cmap, 1.1, mask=gt["outside"], seg_mask=no_seam, zorder=6)
    FIG._minute_marks(axz, gt, late, lambda p: FIG.lonlat(p), label=False)
    axz.plot([impact["lon_deg"]], [impact["lat_deg"]], marker="*", ms=11, mfc="#F0E442", mec="k", mew=0.6, zorder=9)
    axz.annotate("Rosetta impact\n%s UTC\nlon %.2f°, lat %.2f°" % (
        impact["utc"].split("T")[1], impact["lon_deg"], impact["lat_deg"]),
        xy=(impact["lon_deg"], impact["lat_deg"]), xytext=(-75, 45), textcoords="offset points",
        fontsize=6.3, ha="center", arrowprops={"arrowstyle": "-|>", "lw": 0.7, "color": "k", "shrinkB": 6},
        bbox={"facecolor": "white", "edgecolor": "0.3", "lw": 0.5, "boxstyle": "round,pad=0.25"}, zorder=10)
    axz.set_xlim(*zl)
    axz.set_ylim(*za)
    axz.set_xlabel("longitude [°]")
    axz.set_ylabel("latitude [°]")
    axz.text(0.01, 0.97, "(b)", transform=axz.transAxes, va="top", fontsize=8,
             bbox={"facecolor": "white", "edgecolor": "none", "pad": 1})
    region_legend((0.08, 0.0), 2, False)
    sm, extend = layer.mappable()                               # changed: surface colour bar, altitude bar below it
    cb = fig.colorbar(sm, cax=fig.add_axes([0.70, 0.055, 0.28, 0.014]), orientation="horizontal", extend=extend)
    cb.set_label(layer.cb_label, fontsize=7)
    cb = fig.colorbar(lc, cax=fig.add_axes([0.70, -0.015, 0.28, 0.014]), orientation="horizontal")
    cb.set_label("spacecraft altitude [m]", fontsize=7)
    if title:
        fig.suptitle(title, fontsize=8, y=1.0)
    written += FIG.save(fig, str(prefix) + "_latlon", overwrite)
    plt.close(fig)
    written.append(FIG.write_caption(str(prefix) + "_latlon", captions["latlon"], overwrite))
    return written


# ============================================================================================ captions
def gauge_offset(base):
    """GaugeOffsetSC [m] of a base run (boundary.in, else input.snapshot.in); None if not found."""
    for name in ("boundary.in", "input.snapshot.in"):
        path = Path(base) / name
        if path.exists():
            for line in path.read_text().splitlines():
                key, _, value = line.partition("=")
                if key.strip() == "GaugeOffsetSC":
                    return [float(v) for v in value.split("#")[0].split()]
    return None


def fmt_vec(v):
    return " ".join("%g" % x for x in v) if v is not None else "?"


def track_text(ctx):
    gt, imp, tab = ctx["gt"], ctx["gt"]["impact"], ctx["gt"]["table"]
    obs = (tab["refNG"] > 0) | (tab["refRG"] > 0)

    def kind(i):
        return "a COPS observation row" if obs[i] else "a model-only point of the 1-s late grid (no data)"
    w = gt["window"]
    off = ctx["track_offset"]
    if off is not None and not any(off):
        where = "the spacecraft reference point (no gauge offset)"
    else:
        where = "the COPS gauges, offset %s m from the spacecraft reference point (spacecraft frame)" % fmt_vec(off)
    text = (
        "Line: ground track of Rosetta on 2016-09-30, the sub-spacecraft intercept (first surface hit on the line from "
        "the spacecraft toward the body-frame origin), from %s UTC at %.1f km altitude (square) to impact; thick within "
        "the late-visibility window %s–%s UTC (%d rows of the table), with 1-minute marks from its start. Star: impact "
        "point, the straight-line extrapolation of the last two rows of the table (%s; last row %.2f m above the "
        "surface, %.2f m/s) to the shape model: %s UTC, longitude %.3f°, latitude %.3f° "
        "(planetocentric, body-fixed frame of the SPICE kernels). Track and impact point come from the table of the "
        "base %s (%d rows, %d of them observations; positions and altitudes of %s), the same in the late-visibility, "
        "activity and production-rate versions of the figures." % (
            tab["time"][0].strftime("%H:%M:%S"), tab["altitude"][0] / 1000, w["start"].split("T")[1],
            w["stop"].split("T")[1], len(w["points"]),
            ("the COPS observations of %s and %s UTC" % (tab["utc"][-2].split("T")[1], tab["utc"][-1].split("T")[1])
             if obs[-2] and obs[-1] else "%s UTC, %s, and %s UTC, %s" % (
                 tab["utc"][-2].split("T")[1], kind(-2), tab["utc"][-1].split("T")[1], kind(-1))),
            tab["altitude"][-1], imp["speed_m_s"], imp["utc"].split("T")[1],
            imp["lon_deg"], imp["lat_deg"], ctx["track_base"], len(tab["utc"]), int(obs.sum()), where))
    if ctx["track_base"] != ctx["fit_base"]:
        text += (" The fit itself was made on the base %s (COPS gauges at GaugeOffsetSC = %s m, i.e. 1.5 m from the "
                 "spacecraft reference point toward the comet, plus model-only 1-s grid rows in the late window); "
                 "the shape model and its facets are the same (elements.bin identical)." % (
                     ctx["fit_base"], fmt_vec(ctx["fit_offset"])))
    return text


def mesh_text(ctx):
    return ("on the %s shape model (%s facets), longitude/latitude, outermost surface where overhangs put several "
            "surface layers at one longitude/latitude. Left: global; right: the impact area (box in the left panels; "
            "the same box as in the late-visibility version of this figure)" % (
                Path(ctx["gt"]["shape"]).name, thousands(len(ctx["gt"]["faces"]))))


def model_text(ctx):
    f = ctx["fields"]
    nb = len(f["mult"])
    mult = ", ".join("[%g,%g%s %s" % (a, b, "]" if k == nb - 1 else ")", "%.3g" % m)
                     for k, (a, b, m) in enumerate(zip(f["edges"][:-1], f["edges"][1:], f["mult"])))
    return ("Model: %s (its last fine "
            "tuning applied), fitted to the COPS NG and RG data of 2016-09-30 on this model base (COPS gauges 1.5 m "
            "toward the comet). Its H2O activity factor is A = m_b s ln(1 + exp(L/s)) exp(h): L = the "
            "spherical-harmonic expansion plus six held correction templates, s = 1e-5 times the facet's original "
            "activity, m_b the multiplier of the facet's effective-gravitational-slope bin (%s; bin edges in degrees%s), "
            "h the accumulated per-facet fine-tuning term. The emitted H2O flux is F_H2O = q_H2O A; the code gives "
            "units only to F (m⁻² s⁻¹), not separately to q and A." % (MODEL_NAME, mult, MULT_NOTE))


def seen_text(ctx):
    f, s = ctx["fields"], ctx["stats"]["A"]
    seen = f["seen"]
    return ("%s of the %s facets (%s of the area) contribute to at least one fit point ('seen'); the activity of the "
            "other %s facets is not constrained by the data and is shown as the model has it." % (
                thousands(seen.sum()), thousands(len(seen)), pct(np.sum(f["area"][seen]) / np.sum(f["area"])),
                thousands((~seen).sum())))


def scale_text(ctx, mode, what, layers):
    rng = "; ".join("%s: 10^%g to 10^%g" % (L.title_short, round(L.lo, 1), round(L.hi, 1)) for L in layers)
    if mode == "nonzero":
        t = ("Colour: log10 of the value on a linear colour scale from the smallest to the largest non-zero value of "
             "each quantity (%s%s). " % (rng, "; the same range at both times" if what == "production" else ""))
    elif mode == "fixed":
        t = ("Colour (user decision 2026-11-12): log10 of the value on a linear colour scale over a fixed range per "
             "quantity (%s; the same at both times); non-zero values below the range are green (relief-shaded), "
             "values above it get the top colour (arrow on the colour bar). " % rng)
    else:
        t = ("Colour (proposal, not the user's choice): log10 of the value on a linear colour scale over %g decades "
             "below the decade above the %g-th percentile of the non-zero values (%s); non-zero values below the "
             "range are green (relief-shaded), values above it get the top colour (arrows). " % (
                 ctx["args"].span_decades, ctx["args"].top_percentile, rng))
    return t + ("Facets whose value is exactly 0 are grey, relief-shaded like the background of the other rows; "
                "coloured facets are drawn without relief shading so that the colour matches the colour bar.")


def slope_text(ctx):
    g = ctx["meta"].get("gravity", {})
    return ("Effective gravitational slope: the angle between the inward facet normal and the effective gravity "
            "(self-gravity of a homogeneous polyhedron, density %s kg/m³, evaluated %s m above each facet centroid, "
            "minus the centrifugal term ω×(ω×r) of the SPICE rotation at %s). Values above 90° are overhangs. Cliff "
            "templates: slope ≥ 45° (light) and ≥ 60° (dark). Gravity cache: %s (genus-4 mesh admitted as an "
            "experiment)." % (g.get("density_kg_m3", "?"), g.get("evaluation_offset_m", "?"),
                              ", ".join(g.get("omega_sample_epochs_utc", [])) or "the late interval",
                              ctx["cache_path"]))


def activity_numbers(ctx):
    s = ctx["stats"]["A"]
    A = ctx["fields"]["A"]
    f = ctx["fields"]
    small = (A > 0) & (A < 1e10)
    tiny_bin = f["mult"][f["bins"]] < 1e-20
    return ("A = 0 on %s of all facets (%s of the area) and on %s of the seen facets (%s of their area); A > 0 on %s "
            "facets, %s of them unseen. The non-zero values span %s to %s: %s facets have 0 < A < 1e10 (%s of them in "
            "the slope bins whose multiplier is below 1e-20), %s have A >= 1e10." % (
                pct(s["zero_frac"]), pct(s["zero_area_frac"]), pct(s["zero_frac_seen"]), pct(s["zero_area_frac_seen"]),
                thousands(s["nonzero"]), thousands(s["nonzero_unseen"]), fmt_e(s["min_nonzero"]), fmt_e(s["max"]),
                thousands(small.sum()), thousands((small & tiny_bin).sum()), thousands((A >= 1e10).sum())))


def time_text(ctx, lab):
    t = ctx["fields"]["times"][lab]
    qh, qc = t["q"]["H2O"], t["q"]["CO2"]
    n = t["facets"]
    return ("Time %s UTC (r_h = %.6f AU). Instantaneous emitted flux per facet, as the model computes it: F_H2O = "
            "q_H2O A and F_CO2 = κ q_CO2 A_CO2, total = F_H2O + F_CO2, all in m⁻² s⁻¹. q_s = q_min,s on the night side "
            "(cos SZA ≤ 0) and in cast shadow (Sun ray from the facet blocked by the nucleus), else q_min,s + (q_max,s "
            "− q_min,s) cos SZA, with q_max,H2O = 600 / r_h^%.7g = %.4f and q_min,H2O = 0.2 q_max,H2O, q_max,CO2 = "
            "600 / r_h² = %.4f and q_min,CO2 = 0.1 q_max,CO2 (r_h in AU); facets at this time: %s day, %s night, %s "
            "in cast shadow. A_CO2: the CO2 activity map of the base (max(0, the spherical-harmonic CO2 coefficients "
            "of srcCG)); κ = %.4g: the fitted CO2 scale. Total production Q_H2O = Σ F_H2O area = %.4e s⁻¹, Q_CO2 = "
            "%.4e s⁻¹ (Q_CO2 / Q_H2O = %.1e)" % (
                t["utc"].split("T")[1], t["rh_AU"], t["p"], qh[1], qc[1], thousands(n["day"]), thousands(n["night"]),
                thousands(n["cast_shadow"]), ctx["fields"]["kappa"], t["Q_H2O"], t["Q_CO2"],
                t["Q_CO2"] / t["Q_H2O"]) + (
            # data-driven (2026-11-13): the 'negligible' statement only where it holds
            ": the CO2 production of this fit is negligible; it is shown as it is (maximum F_CO2 %.2e m⁻² s⁻¹), and the "
            "total equals F_H2O to machine precision wherever F_H2O > 0 (it differs from F_H2O only where A = 0 and "
            "A_CO2 > 0, where it equals F_CO2)." % t["F_CO2"].max() if t["Q_CO2"] / t["Q_H2O"] < 1e-12 else
            "; maximum F_CO2 %.2e m⁻² s⁻¹." % t["F_CO2"].max()))


def production_numbers(ctx, lab):
    out = []
    for key, name in (("F_H2O", "F_H2O"), ("F_CO2", "F_CO2"), ("total", "total")):
        s = ctx["stats"]["%s_%s" % (key, lab)]
        out.append("%s = 0 on %s of all facets (%s of the area), non-zero values %s to %s" % (
            name, pct(s["zero_frac"]), pct(s["zero_area_frac"]), fmt_e(s["min_nonzero"]), fmt_e(s["max"])))
    f = ctx["fields"]
    co2_only = int(np.sum((f["A"] == 0) & (f["times"][lab]["F_CO2"] > 0)))
    return "; ".join(out) + (". The total equals F_CO2 alone (A = 0, A_CO2 > 0) on %s facets (%s of all facets)." % (
        thousands(co2_only), pct(co2_only / len(f["A"]))))


def gt_caption_3d(ctx, surface):
    return ("End-of-mission ground track of Rosetta and its impact point on the %s shape model (%s facets), "
            "orthographic views rendered on the GPU (z-buffer). (a) The whole nucleus seen down the final approach "
            "(view direction: from the start of the late-visibility window toward the impact point); (b) the impact "
            "region (box in a); (c) close-up of the late descent (box in b), with labelled 1-minute marks (UTC, "
            "altitude). Scale bars in each panel. Surface: %s Track: coloured by altitude (log colour bar 'spacecraft "
            "altitude'), thick within the late-visibility window, dotted grey where hidden behind the nucleus. %s" % (
                Path(ctx["gt"]["shape"]).name, thousands(len(ctx["gt"]["faces"])), surface, track_text(ctx)))


def gt_caption_latlon(ctx, surface):
    return ("End-of-mission ground track of Rosetta and its impact point on the %s shape model (%s facets), "
            "longitude/latitude (equirectangular), outermost surface where overhangs put several surface layers at "
            "one longitude/latitude. (a) Global map; (b) the impact area (box in a). Surface: %s Track: coloured by "
            "altitude (log colour bar 'spacecraft altitude'), thick within the late-visibility window. %s" % (
                Path(ctx["gt"]["shape"]).name, thousands(len(ctx["gt"]["faces"])), surface, track_text(ctx)))


def visibility_surface_text(ctx):
    m = ctx["gt"]["masks"]
    meta = ctx["meta"]
    n = meta.get("templates", {}).get("late_ng_not_rg", {}).get("selected_point_count")
    return ("NG-only facets (orange, %s) and RG-only facets (blue, %s) of the late-visibility templates of the "
            "template bundle %s: facets seen by the nude gauge (Kn > 0) and not by the ram gauge (KΦ = KV = 0) from "
            "at least one of the %s COPS observation points of %s–%s UTC, resp. the reverse; 'seen' = nonzero "
            "collisionless transport kernel. %s Other facets grey with relief shading." % (
                thousands(m["late_ng_not_rg"].sum()), thousands(m["late_rg_not_ng"].sum()), ctx["args"].bundle, n,
                meta.get("templates", {}).get("late_ng_not_rg", {}).get("interval_utc", ["?", "?"])[0].split("T")[1],
                meta.get("templates", {}).get("late_ng_not_rg", {}).get("interval_utc", ["?", "?"])[1].split("T")[1],
                ("The masks were computed on the bundle's own base %s, the base of the track." % meta.get("base_run")
                 if Path(str(meta.get("base_run"))).resolve() == Path(ctx["track_base"]).resolve() else
                 "These masks were computed on the bundle's base %s and were not recomputed for the track base."
                 % meta.get("base_run"))))


# ============================================================================================ main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--base", required=True, help="the fit's base run (fields, maps checks)")
    ap.add_argument("--track-base", default=None,
                    help="base run whose table gives the ground track, late window and impact point (default --base); "
                         "its elements.bin must equal the fit base's")
    ap.add_argument("--range", action="append", default=[],
                    help="KEY=LO,HI (log10) fixed colour range of A, F_H2O, F_CO2, total (all four; mode 'fixed')")
    ap.add_argument("--bundle", required=True)
    ap.add_argument("--ft-dir", required=True)
    ap.add_argument("--fit-dir", required=True)
    ap.add_argument("--maps-dir", required=True)
    ap.add_argument("--geometry", required=True)
    ap.add_argument("--maps-code", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--visibility-prefix", default=None)
    ap.add_argument("--preview-dir", default=None)
    ap.add_argument("--span-decades", type=float, default=10.0)
    ap.add_argument("--top-percentile", type=float, default=99.99)
    ap.add_argument("--gravity-cache", default=None)
    ap.add_argument("--device", type=int, default=0)
    ap.add_argument("--sets", default="activity,production",
                    help="comma list of activity, production (sets 2 and 3); a smoke test uses 'activity' alone")
    args = ap.parse_args()
    sets = set(args.sets.split(","))
    args.track_base = args.track_base or args.base
    fixed = {}
    for item in args.range:
        key, _, val = item.partition("=")
        lo, hi = (float(v) for v in val.split(","))
        fixed[key] = (lo, hi)
    if fixed and set(fixed) != {"A", "F_H2O", "F_CO2", "total"}:
        raise SystemExit("--range: give A, F_H2O, F_CO2 and total")
    main_mode = "fixed" if fixed else "nonzero"
    if not sets <= {"activity", "production"}:
        raise SystemExit("--sets: activity and/or production")
    out = Path(args.out_dir)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    stage = out / ("pub_maps.incomplete_%s" % stamp)
    targets = [out / "checks"] + [out / d for k, d in (("activity", "activity"), ("production", "production_rate"))
                                  if k in sets]
    if args.preview_dir:
        targets.append(Path(args.preview_dir))
    vis_caps = []
    if args.visibility_prefix:
        vis_caps = [Path(args.visibility_prefix + "_3d_caption.txt"), Path(args.visibility_prefix + "_latlon_caption.txt")]
    for t in targets + vis_caps:
        if t.exists():
            raise SystemExit("%s exists (never overwritten; move it aside first)" % t)
    stage.mkdir(parents=True)
    say("stage %s" % stage)

    # ground track, impact, masks, window: exactly as tools/pub_figures.py templates / groundtrack
    base_tab = FIG.read_table(args.track_base)
    window = FIG.late_window(base_tab, args.bundle, None, None)
    gt = FIG.groundtrack_data(args.track_base, None, args.bundle, window, args.device)
    say("ground track: %d rows, impact %s UTC lon %.3f lat %.3f; window %s .. %s (%d rows)" % (
        len(gt["table"]["utc"]), gt["impact"]["utc"], gt["impact"]["lon_deg"], gt["impact"]["lat_deg"],
        window["start"], window["stop"], len(window["points"])))
    fields = load_fields(args)
    if not np.array_equal(fields["normal"], gt["normals"]):
        raise SystemExit("elements.bin normals differ between the figure mesh and the fields")
    if Path(args.track_base).resolve() != Path(args.base).resolve():
        nm = load_module(args.maps_code, "nonneg_maps_frozen_tb")
        e2 = nm.read_elements(Path(args.track_base) / "elements.bin")
        same = {k: bool(np.array_equal(e2[k], fields[{"normal": "normal", "area": "area"}.get(k, k)]))
                for k in ("normal", "area")}
        e1 = nm.read_elements(Path(args.base) / "elements.bin")
        same["centroid"] = bool(np.array_equal(e1["centroid"], e2["centroid"]))
        same["ACO2"] = bool(np.array_equal(e1["ACO2"], e2["ACO2"]))
        del e1, e2
        fields["checks"]["track_base_elements_equal_fit_base"] = same
        if not all(same.values()):
            raise SystemExit("track base %s has another mesh than the fit base %s: %s" % (args.track_base, args.base, same))
        say("track base %s: elements.bin fields equal to the fit base's (%s)" % (args.track_base, same))
    if not fields["checks"]["passed"]:
        (stage / "checks").mkdir()
        json.dump(fields["checks"], open(stage / "checks" / "fields_verification.json", "w"), indent=1)
        raise SystemExit("field checks FAILED (see %s/checks)" % stage)
    say("field checks passed")
    seen, area = fields["seen"], fields["area"]
    stats = {"A": value_stats(fields["A"], seen, area)}
    for lab, _ in TIMES:
        for key in ("F_H2O", "F_CO2", "total"):
            stats["%s_%s" % (key, lab)] = value_stats(fields["times"][lab][key], seen, area)
    meta = FIG.bundle_metadata(args.bundle)
    _, cache_path = FIG.gravity_slope(args.bundle, args.gravity_cache, len(gt["faces"]))
    ctx = {"gt": gt, "fields": fields, "stats": stats, "meta": meta, "cache_path": cache_path, "args": args,
           "track_base": str(Path(args.track_base).resolve()), "fit_base": str(Path(args.base).resolve()),
           "track_offset": gauge_offset(args.track_base), "fit_offset": gauge_offset(args.base)}
    plt = FIG.style()
    cmap = value_cmap(plt)
    say("relief intensity of the lon/lat renderer")
    I_ll = lonlat_intensity(gt["vertices"], gt["faces"], gt["normals"])

    def ranges(mode):
        if mode == "fixed":
            return dict(fixed)
        r = {"A": colour_range([fields["A"]], mode, args.span_decades, args.top_percentile)}
        for key in ("F_H2O", "F_CO2", "total"):
            r[key] = colour_range([fields["times"][lab][key] for lab, _ in TIMES], mode, args.span_decades,
                                  args.top_percentile)
        return r

    def layer(key, title, short, values, rng, cb_label, mode):
        L = ValueLayer(key, title, values, rng[0], rng[1], cb_label, cmap, below_mode=(mode in ("span", "fixed")))
        L.title_short = short
        return L

    def surface_entries(layers, zero_label):
        from matplotlib.patches import Patch
        h = [Patch(facecolor=ZERO_PATCH, edgecolor="none", label=zero_label)]
        if any(L.below for L in layers):
            h.append(Patch(facecolor=BELOW_COLOR, edgecolor="none",
                           label="non-zero, below the colour range"))   # thresholds: colour bar and caption
        return h

    def make_set(mode, dest, only_templates):
        rg = ranges(mode)
        written = []
        tag = "_span%g" % args.span_decades if mode == "span" else ""
        # ---- set 2: activity
        LA = layer("act", "H₂O activity factor A (softplus model, final state)", "A", fields["A"], rg["A"],
                   "log₁₀ A", mode)
        hdr = ("Surface H₂O activity factor A of the softplus activity model, final state after 5 fit and fine-tuning "
               "cycles\n(H₂O flux F_H₂O = q_H₂O A; the code gives units only to F, in m⁻² s⁻¹).\nColour: log₁₀ A %s; "
               "grey: A = 0 (%s of all facets,\n%s of the facets seen from the fit points). Rows (c)–(f) as in the "
               "late-visibility figure." % (
                   "over its non-zero values, $10^{%.1f}$ to $10^{%.1f}$" % (LA.lo, LA.hi) if mode == "nonzero" else
                   "over the %srange $10^{%g}$ to $10^{%g}$; green: 0 < A < $10^{%g}$" % (
                       "proposed " if mode == "span" else "", LA.lo, LA.hi, LA.lo),
                   pct(stats["A"]["zero_frac"]), pct(stats["A"]["zero_frac_seen"])))
        cap = ("Surface H2O activity of the softplus activity model with the effective gravitational slope and the "
               "cliff templates, %s. (a, b) H2O activity factor A. %s %s %s %s (c, d) Effective gravitational slope "
               "and (e, f) cliff templates, as in the late-visibility version of this figure. %s %s" % (
                   mesh_text(ctx), model_text(ctx), activity_numbers(ctx), seen_text(ctx),
                   scale_text(ctx, mode, "activity", [LA]), slope_text(ctx), track_text(ctx)))
        zero_lab = "A = 0 (grey, relief-shaded)"
        if "activity" not in sets:
            only_templates_act = None
        else:
            only_templates_act = only_templates
        if only_templates_act is not None:
            written += values_templates_figure(gt, dest / "activity" / ("12M_templates_activity%s" % tag), [LA], hdr,
                                               surface_entries([LA], zero_lab), cap, I_ll, args.gravity_cache,
                                               args.device)
            say("wrote %s" % written[-1])
        if only_templates_act is False:
            surf = ("the H2O activity factor A of the softplus activity model, final state (colour bar 'log10 A'). "
                    "%s %s %s %s" % (model_text(ctx), activity_numbers(ctx), seen_text(ctx),
                                     scale_text(ctx, mode, "activity", [LA])))
            written += values_groundtrack_figures(
                gt, dest / "activity" / ("12M_groundtrack_activity%s" % tag), LA, surface_entries([LA], zero_lab),
                {"3d": gt_caption_3d(ctx, surf), "latlon": gt_caption_latlon(ctx, surf)}, I_ll, args.device)
            say("wrote %s" % written[-1])
        # ---- set 3: production rate at each time
        for lab, utc in (TIMES if "production" in sets else ()):
            hhmmss = utc.split("T")[1].replace(":", "")
            hm = utc.split("T")[1]
            t = fields["times"][lab]
            LH = layer("h2o", "H₂O production rate, %s UTC" % hm, "F_H2O", t["F_H2O"], rg["F_H2O"],
                       "log₁₀ F_H₂O [m⁻² s⁻¹]", mode)
            LC = layer("co2", "CO₂ production rate, %s UTC" % hm, "F_CO2", t["F_CO2"], rg["F_CO2"],
                       "log₁₀ F_CO₂ [m⁻² s⁻¹]", mode)
            LT = layer("tot", "total production rate (H₂O + CO₂), %s UTC" % hm, "total", t["total"], rg["total"],
                       "log₁₀ F_total [m⁻² s⁻¹]", mode)
            hdr = ("Production rate (emitted molecular flux) of the softplus activity model, final state after 5 fit "
                   "and\nfine-tuning cycles, 2016-09-30 %s UTC (r_h = %.4f AU); night side and cast shadow as in the "
                   "model.\nColour: log₁₀ of the flux %s;\ngrey: zero. CO₂ shown as it is (Q_CO₂ / Q_H₂O = %.0e). "
                   "Rows (g)–(j) as in the late-visibility figure." % (
                       hm, t["rh_AU"], "over the non-zero values of each quantity at both times" if mode == "nonzero"
                       else "over the %scolour-bar ranges (same at both times); green: non-zero below the range" % (
                           "proposed " if mode == "span" else ""), t["Q_CO2"] / t["Q_H2O"]))
            cap = ("Production rate (instantaneous emitted molecular flux per facet) of the softplus activity model at "
                   "%s UTC on 2016-09-30, with the effective gravitational slope and the cliff templates, %s. (a, b) "
                   "H2O, (c, d) CO2, (e, f) total production rate. %s %s %s %s %s (g, h) Effective gravitational slope and "
                   "(i, j) cliff templates, as in the late-visibility version of this figure. %s %s" % (
                       hm, mesh_text(ctx), time_text(ctx, lab), model_text(ctx), production_numbers(ctx, lab),
                       seen_text(ctx), scale_text(ctx, mode, "production", [LH, LC, LT]), slope_text(ctx),
                       track_text(ctx)))
            zero_lab = "zero flux (grey, relief-shaded)"
            written += values_templates_figure(
                gt, dest / "production_rate" / ("12M_templates_production_rate_%sUTC%s" % (hhmmss, tag)),
                [LH, LC, LT], hdr, surface_entries([LH, LC, LT], zero_lab), cap, I_ll, args.gravity_cache,
                args.device)
            say("wrote %s" % written[-1])
            if not only_templates:
                surf = ("the total production rate F_H2O + F_CO2 [m⁻² s⁻¹] of the softplus activity model at %s UTC "
                        "(colour bar). %s %s %s %s %s" % (
                            hm, time_text(ctx, lab), model_text(ctx), production_numbers(ctx, lab), seen_text(ctx),
                            scale_text(ctx, mode, "production", [LT])))
                written += values_groundtrack_figures(
                    gt, dest / "production_rate" / ("12M_groundtrack_production_rate_%sUTC%s" % (hhmmss, tag)), LT,
                    surface_entries([LT], zero_lab),
                    {"3d": gt_caption_3d(ctx, surf), "latlon": gt_caption_latlon(ctx, surf)}, I_ll, args.device)
                say("wrote %s" % written[-1])
        return written, rg

    if vis_caps:
        surf = visibility_surface_text(ctx)
        (stage / "visibility").mkdir()
        FIG.write_caption(stage / "visibility" / (Path(args.visibility_prefix).name + "_3d"), gt_caption_3d(ctx, surf))
        FIG.write_caption(stage / "visibility" / (Path(args.visibility_prefix).name + "_latlon"),
                          gt_caption_latlon(ctx, surf))
    written, rg_main = make_set(main_mode, stage, False)
    rg_prev = None
    if args.preview_dir:
        pv = stage / "preview"
        w2, rg_prev = make_set("span", pv, True)
        written += w2
    # ---- checks
    (stage / "checks").mkdir()
    json.dump(fields["checks"], open(stage / "checks" / "fields_verification.json", "w"), indent=1)
    with open(stage / "checks" / "value_statistics.txt", "w") as f:
        f.write("# pub_maps_figures.py (%s): per-facet value statistics of the softplus final state; seen = the fit's "
                "nonneg_facets.npz seen; area fractions over all facets and over the seen facets\n" % SCHEMA)
        for k, s in stats.items():
            f.write("\n[%s]\n%s\n" % (k, json.dumps(s, indent=1)))
        f.write("\n[colour ranges, log10, mode %s]\n%s\n" % (main_mode, json.dumps(rg_main, indent=1)))
        if rg_prev:
            f.write("\n[colour ranges, log10, mode span (preview)]\n%s\n" % json.dumps(rg_prev, indent=1))
    manifest = {"schema": SCHEMA, "created_utc": stamp, "host": socket.gethostname(),
                "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "command": sys.argv,
                "script_sha256": sha256_file(__file__), "figures_py_sha256": sha256_file(FIG.__file__),
                "figures_py": str(Path(FIG.__file__).resolve()),
                "inputs_sha256": {str(p): sha256_file(p) for p in (
                    Path(args.ft_dir) / "finetune_state.npz", Path(args.fit_dir) / "fit_restart.npz",
                    Path(args.fit_dir) / "nonneg_facets.npz", Path(args.bundle), Path(args.maps_dir) / "maps_manifest.json",
                    Path(args.base) / "COPS_model.dat", Path(args.track_base) / "COPS_model.dat")},
                "impact": gt["impact"], "window": {k: gt["window"][k] for k in ("start", "stop", "source",
                                                                              "altitude_min_m", "altitude_max_m")},
                "window_rows": len(gt["window"]["points"]), "colour_ranges_log10": rg_main, "colour_mode": main_mode,
                "track_base": ctx["track_base"], "fit_base": ctx["fit_base"],
                "colour_ranges_log10_preview": rg_prev, "colormap": "%s from %g" % (CMAP_NAME, CMAP_FROM),
                "statistics": stats, "field_checks": fields["checks"],
                "files": [str(Path(p).relative_to(stage)) for p in written], "wall_seconds": time.time() - T0}
    json.dump(manifest, open(stage / "checks" / "pub_maps_manifest.json", "w"), indent=1, default=float)
    # ---- move into place
    for sub, dest in (("activity", out / "activity"), ("production_rate", out / "production_rate"),
                      ("checks", out / "checks")):
        if (stage / sub).exists():
            (stage / sub).rename(dest)
    if args.preview_dir:
        Path(args.preview_dir).parent.mkdir(parents=True, exist_ok=True)
        (stage / "preview").rename(Path(args.preview_dir))
    for p in vis_caps:
        (stage / "visibility" / p.name).rename(p)
    if vis_caps:
        (stage / "visibility").rmdir()
    stage.rmdir()
    say("done: %d files" % (len(written) + len(vis_caps) + 3))
    return 0


if __name__ == "__main__":
    sys.exit(main())
