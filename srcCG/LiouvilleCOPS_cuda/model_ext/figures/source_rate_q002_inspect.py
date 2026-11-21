#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/source_rate_q002_inspect.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/source_rate_q002_inspect.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 0b903f2f2b20eb4da9725107b319e3594c30c94a0999c317f93bd1d2aa63ba1f
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: inspection figure of the H2O source rate after the night-side fix: new-law F_H2O | old-law F_H2O of the
#   counterpart (common log classes) | new q/q_max, per map time, with numbers (Q_H2O total / seen, night + shadow
#   share, active area) and checks (2026-11-18)
# RUN: GPU job jobs/source_rate_q002.sbatch (STEPS inspect); --numbers-only on CPU
# ==== end of header ====
"""Inspection of the H2O source rate after the H2O night-side fix (2026-11-18): the production rate of a state fitted with
the NEW law (q_min = 0.02 q_max, srcCG's value) next to the production rate of its OLD-law counterpart (q_min = 0.2 q_max,
the law of every fit and figure before 2026-11-18, not intended) and next to the new illumination factor q/q_max.

User request (2026-11-18, relayed by the lead): "prepare the surface source rate maps for completed runs for inspection of
the source rate". Lead: per completed run the source-rate figure set (pub_maps_v2_figures.py) plus, per time, an
inspection figure: new-law production rate side by side with the old-law production rate and with the new q/q_max.

    source_rate_q002_inspect.py --time HH:MM:SS --out-dir <run folder> \
        --new-base B --new-ft D --new-fit D --new-maps D --new-geometry G --new-maps-code nonneg_maps.py \
        --old-base B --old-ft D --old-fit D --old-maps D --old-geometry G --old-maps-code nonneg_maps.py \
        --track-base B --bundle restart_templates.npz --new-name "..." --old-name "..." [--low-cut-area 0.02]
        [--template-relief 0.7] [--classes 15,20] [--device 0]
GPU job (the lon/lat rasterizer is cops_cuda.figures.rasterize, a CUDA kernel). Run from <AMPS>/srcCG/LiouvilleCOPS_cuda.

FIGURE (lon/lat, global and the impact area, the layout of illumination_maps_figures.py's side-by-side figure)
  (a, b) F_H2O of the new-law state, (c, d) F_H2O of the old-law state: colour only on the facets seen during the fit
         (seen F = 0 light green, unseen grey, relief-shaded), ONE set of log10 classes for both rows per window
         (pub_maps_v2_figures.Scale on the union of the two states' displayed seen non-zero values, area cut
         --low-cut-area at the bottom, top kept), so that the night / shadow drop and the refit's changes are visible
  (e, f) the new law's illumination factor q/q_max on the whole surface (illumination_maps_figures' factor classes:
         0.02 + 0.98 cos SZA on day facets in 0.1 classes, 0.02 on night and cast-shadow facets)
NUMBERS (<out>/inspection/checks/inspect_<hhmmss>.{json,txt}) per state: Q_H2O total and of the seen facets, the share of
  Q_H2O from night and cast-shadow facets (whole, seen), their area shares, the active area (A > 0; whole, seen), the
  law (q_min / q_max from the maps' own q). CHECKS: both states' fields reproduce their maps' .vtu (pub_maps_figures
  .load_fields), the same mesh, cos SZA and illumination equal between the two geometries (bitwise), the factor law
  (illumination_maps_figures.factor_checks: F = q A bitwise, q_min / q_max = the maps' H2O_QMIN_RATIO).
CODE REUSE: pub_maps_figures.py, pub_maps_v2_figures.py, illumination_maps_figures.py imported unchanged; the figure
  layout is a copy (with changes: rows) of illumination_maps_figures.cmd_figures' side-by-side figure.
Refuses existing products.
"""
import argparse
import json
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import illumination_maps_figures as IM  # noqa: E402

V2, M, FIG = IM.V2, IM.M, IM.FIG
T0 = time.time()
SCHEMA = "source_rate_q002_inspect_v1"


def say(t):
    print("source_rate_q002_inspect [%6.0f s]: %s" % (time.time() - T0, t), flush=True)


def state_numbers(F, A, seen, area, ill, q):
    """Numbers of one state at one time (F = F_H2O per facet, ill: 0 day, 1 night, 2 cast shadow)."""
    P = F * area
    dark = ill != 0
    Q, Qs = float(P.sum()), float(P[seen].sum())
    act = A > 0
    return {"Q_H2O_s-1": Q, "Q_H2O_seen_s-1": Qs, "seen_share_of_Q": Qs / Q if Q > 0 else float("nan"),
            "dark_share_of_Q": float(P[dark].sum() / Q) if Q > 0 else float("nan"),
            "night_share_of_Q": float(P[ill == 1].sum() / Q) if Q > 0 else float("nan"),
            "cast_shadow_share_of_Q": float(P[ill == 2].sum() / Q) if Q > 0 else float("nan"),
            "dark_share_of_seen_Q": float(P[dark & seen].sum() / Qs) if Qs > 0 else float("nan"),
            "dark_area_fraction": float(area[dark].sum() / area.sum()),
            "dark_area_fraction_seen": float(area[dark & seen].sum() / area[seen].sum()),
            "active_area_m2": float(area[act].sum()), "active_area_fraction": float(area[act].sum() / area.sum()),
            "active_area_fraction_seen": float(area[act & seen].sum() / area[seen].sum()),
            "active_facets": int(act.sum()), "seen_area_m2": float(area[seen].sum()), "total_area_m2": float(area.sum()),
            "q_min_H2O": float(q[0]), "q_max_H2O": float(q[1]), "q_ratio_H2O": float(q[0] / q[1])}


def old_fields_from_maps(a, lab, utc):
    """The OLD-law state's fields as its maps wrote them (2026-11-18): the old maps code's law guard refuses to run since
    the reference source code (fit_cops) has the new law, so the old fields are read from the old maps' own .vtu files
    (the files its full-set figures were drawn from) and checked: A and seen against the state's files (bitwise), the
    static and flux files against the manifest (sizes, sha256 not recomputed), F_H2O = q A with the manifest's q_H2O
    (q_min on night / cast-shadow facets, q_min + (q_max - q_min) min(cos SZA, 1) on day facets), Q_H2O against the
    manifest."""
    maps = Path(a.old_maps)
    man = json.load(open(maps / "maps_manifest.json"))
    code_sha = M.sha256_file(a.old_maps_code)
    chk = {"source": "old maps' .vtu files (law guard of the old maps code refuses the new reference law)",
           "maps_code_sha256_equal_manifest": code_sha == man["code_sha256"]["nonneg_maps.py"]}
    for key, d in (("ft_dir_resolved", a.old_ft), ("fit_dir_resolved", a.old_fit)):
        chk[key + "_equal"] = Path(man["state"][key]) == Path(d).resolve()
    chk["base_equal"] = Path(man["base"]) == Path(a.old_base).resolve()
    chk["geometry_equal"] = Path(man["geometry"]) == Path(a.old_geometry).resolve()
    nm = M.load_module(a.old_maps_code, "nonneg_maps_old_reader")
    el = nm.read_elements(Path(a.old_base) / "elements.bin")
    with np.load(Path(a.old_ft) / "finetune_state.npz", allow_pickle=False) as z:
        A = np.asarray(z["A_final"], np.float64)
    with np.load(Path(a.old_fit) / "nonneg_facets.npz", allow_pickle=False) as z:
        seen = np.asarray(z["seen"]).astype(bool)
    rd, _, _ = nm.read_vtu(maps / ("%s_static.vtu" % man["model"]))
    chk["A_equal_vtu_all_facets"] = bool(np.array_equal(rd["A_H2O"], A))
    chk["seen_equal_vtu"] = bool(np.array_equal(rd["seen"].astype(bool), seen))
    del rd
    mt = [t for t in man["times"] if t["label"] == lab][0]
    if mt["utc"] != utc:
        raise SystemExit("old maps time %s is %s, not %s" % (lab, mt["utc"], utc))
    rd, fd, _ = nm.read_vtu(maps / mt["vtu"])
    F = np.asarray(rd["F_H2O_per_m2_s"], np.float64)
    cos = np.asarray(rd["cos_SZA"], np.float64)
    ill = np.asarray(rd["illumination"]).astype(np.int64)
    del rd
    qmin, qmax = (float(v) for v in mt["q_H2O_min_max"])
    q = np.where(ill != 0, qmin, qmin + (qmax - qmin) * np.minimum(cos, 1.0))
    chk["F_equal_q_times_A_bitwise"] = bool(np.array_equal(F, q * A))
    chk["F_max_rel_diff_q_times_A"] = M.max_rel(F, q * A)
    Q = float(np.sum(F * el["area"]))
    chk["Q_H2O_rel_diff_manifest"] = abs(Q / mt["Q_H2O_s-1"] - 1.0)
    chk["passed"] = bool(chk["maps_code_sha256_equal_manifest"] and chk["ft_dir_resolved_equal"] and
                         chk["fit_dir_resolved_equal"] and chk["base_equal"] and chk["geometry_equal"] and
                         chk["A_equal_vtu_all_facets"] and chk["seen_equal_vtu"] and
                         chk["F_max_rel_diff_q_times_A"] <= M.REL_TOL and chk["Q_H2O_rel_diff_manifest"] <= M.REL_TOL)
    return {"A": A, "seen": seen, "area": el["area"], "normal": el["normal"], "cos": cos,
            "times": {lab: {"F_H2O": F, "illumination": ill, "q": {"H2O": (qmin, qmax)}, "Q_H2O": Q}},
            "checks": dict(chk, times={lab: {"vtu": str(maps / mt["vtu"]), "Q_H2O_manifest": mt["Q_H2O_s-1"]}})}


def numbers_text(r):
    L = ["# source_rate_q002_inspect.py (%s), %s UTC; Q in s^-1, areas in m^2; dark = night (cos SZA <= 0) + cast shadow"
         % (SCHEMA, r["utc"]), ""]
    keys = [("q_ratio_H2O", "q_min / q_max (H2O)", "%.4g"), ("Q_H2O_s-1", "Q_H2O total", "%.4e"),
            ("Q_H2O_seen_s-1", "Q_H2O of the seen facets", "%.4e"), ("seen_share_of_Q", "seen share of Q_H2O", "%.4f"),
            ("dark_share_of_Q", "night + cast-shadow share of Q_H2O", "%.4f"),
            ("night_share_of_Q", "  night share", "%.4f"), ("cast_shadow_share_of_Q", "  cast-shadow share", "%.4f"),
            ("dark_share_of_seen_Q", "night + cast-shadow share of the seen Q_H2O", "%.4f"),
            ("dark_area_fraction", "night + cast-shadow area fraction (whole)", "%.4f"),
            ("dark_area_fraction_seen", "night + cast-shadow area fraction (seen)", "%.4f"),
            ("active_area_m2", "active area A > 0", "%.4e"), ("active_area_fraction", "active area fraction (whole)", "%.4f"),
            ("active_area_fraction_seen", "active area fraction (seen)", "%.4f"), ("active_facets", "active facets", "%d")]
    L.append("%-48s %-16s %-16s" % ("quantity", "new law", "old law"))
    for k, lab, f in keys:
        L.append("%-48s %-16s %-16s" % (lab, f % r["new"][k], f % r["old"][k]))
    L.append("%-48s %.4f" % ("Q_H2O new / old", r["new"]["Q_H2O_s-1"] / r["old"]["Q_H2O_s-1"]))
    L.append("%-48s %.4f" % ("seen Q_H2O new / old", r["new"]["Q_H2O_seen_s-1"] / r["old"]["Q_H2O_seen_s-1"]))
    L.append("")
    c = r["checks"]
    short = {k: v for k, v in c.items() if not isinstance(v, dict)}
    short["fields_new_passed"] = c["fields_new"]["passed"]
    short["fields_old_passed"] = c["fields_old"]["passed"]
    short["old_F_equal_q_times_A_bitwise"] = c["fields_old"].get("F_equal_q_times_A_bitwise")
    short["new_F_equal_vtu_bitwise"] = c["vtu_new"]["F_H2O_equal_vtu_all_facets"]
    short["factor_new_passed"] = c["factor_new"]["passed"]
    if "factor_classes" in c:
        short["factor_classes_passed"] = c["factor_classes"]["passed"]
    L.append("checks (all in the .json): %s" % json.dumps(short, default=str))
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for side in ("new", "old"):
        for x in ("base", "ft", "fit", "maps", "geometry", "maps-code"):
            ap.add_argument("--%s-%s" % (side, x), required=True)
        ap.add_argument("--%s-name" % side, required=True)
    for x in ("--track-base", "--bundle", "--time", "--out-dir"):
        ap.add_argument(x, required=True)
    ap.add_argument("--low-cut-area", type=float, default=0.02)
    ap.add_argument("--template-relief", type=float, default=0.7)
    ap.add_argument("--classes", default="15,20")
    ap.add_argument("--device", type=int, default=0)
    ap.add_argument("--numbers-only", action="store_true",
                    help="CPU: fields, checks and numbers only (inspection/checks/numbers_<hhmmss>.{json,txt}), no figure")
    a = ap.parse_args()
    hit = [t for t in IM.MAP_TIMES if t[0] == a.time]
    if len(hit) != 1:
        raise SystemExit("--time must be one of %s" % [t[0] for t in IM.MAP_TIMES])
    hhmmss, lab, utc = hit[0]
    tag = hhmmss.replace(":", "")
    out = Path(a.out_dir) / "inspection"
    prefix = out / ("12M_inspect_h2o_production_rate_new_vs_old_law_%sUTC" % tag)
    nstem = "numbers" if a.numbers_only else "inspect"
    finals = [out / "checks" / ("%s_%s.%s" % (nstem, tag, s)) for s in ("json", "txt")]
    if not a.numbers_only:
        finals += [prefix.with_name(prefix.name + s) for s in (".png", ".eps", "_caption.txt")]
    for p in finals:
        if p.exists():
            raise SystemExit("%s exists (never overwritten; move it aside first)" % p)
    (out / "checks").mkdir(parents=True, exist_ok=True)
    V2.NMIN, V2.NMAX = (int(v) for v in a.classes.split(","))

    # ---- fields (pub_maps_figures.load_fields for each state, this time only)
    M.TIMES = ((lab, utc),)
    st = {}
    ns = SimpleNamespace(base=a.new_base, ft_dir=a.new_ft, fit_dir=a.new_fit, maps_dir=a.new_maps,
                         geometry=a.new_geometry, maps_code=a.new_maps_code)
    st["new"] = M.load_fields(ns)
    st["old"] = old_fields_from_maps(a, lab, utc)
    for side in ("new", "old"):
        if not st[side]["checks"]["passed"]:
            raise SystemExit("%s state: field checks failed: %s" % (side, st[side]["checks"]))
        say("%s state fields loaded (checks passed; Q_H2O %.4e)" % (side, st[side]["times"][lab]["Q_H2O"]))
    n_, o_ = st["new"], st["old"]
    checks = {"fields_new": {k: v for k, v in n_["checks"].items() if k != "times"},
              "fields_old": {k: v for k, v in o_["checks"].items() if k != "times"},
              "vtu_new": n_["checks"]["times"][lab], "vtu_old": o_["checks"]["times"][lab]}
    for k in ("area", "normal"):
        if not np.array_equal(n_[k], o_[k]):
            raise SystemExit("the two states' meshes differ (%s)" % k)
    nm_new = M.load_module(a.new_maps_code, "nonneg_maps_new_inspect")
    nm_old = M.load_module(a.old_maps_code, "nonneg_maps_old_inspect")
    IM.set_qratio(nm_new)
    nf = len(n_["A"])
    cos, sh = nm_new.geometry_time(a.new_geometry, lab, nf)
    cos_o, sh_o = nm_old.geometry_time(a.old_geometry, lab, nf)
    checks["old_vtu_cos_equal_old_geometry"] = bool(np.array_equal(o_["cos"], cos_o))
    Tn, To = n_["times"][lab], o_["times"][lab]
    ill = Tn["illumination"]
    checks["geometry_cos_equal_bitwise"] = bool(np.array_equal(cos, cos_o))
    checks["geometry_shadow_equal"] = bool(np.array_equal(sh, sh_o))
    checks["illumination_equal"] = bool(np.array_equal(ill, To["illumination"]))
    checks["seen_equal"] = bool(np.array_equal(n_["seen"], o_["seen"]))
    if not (checks["geometry_cos_equal_bitwise"] and checks["illumination_equal"]):
        raise SystemExit("the two geometries differ at %s: %s" % (utc, checks))
    Fn, Fo, An, Ao = Tn["F_H2O"], To["F_H2O"], n_["A"], o_["A"]
    fchk, f = IM.factor_checks(cos, sh, ill, Fn, An, Tn["q"]["H2O"], nm_new)
    checks["factor_new"] = fchk
    if not fchk["passed"]:
        raise SystemExit("factor checks failed: %s" % fchk)
    r_old = To["q"]["H2O"][0] / To["q"]["H2O"][1]
    r_new = Tn["q"]["H2O"][0] / Tn["q"]["H2O"][1]
    checks["q_ratio_new_maps"], checks["q_ratio_old_maps"] = r_new, r_old
    if abs(r_new - IM.QRATIO) > 1e-12 * IM.QRATIO:
        raise SystemExit("new maps' q_min/q_max %r is not the code's %r" % (r_new, IM.QRATIO))
    say("checks: %s" % {k: v for k, v in checks.items() if not isinstance(v, dict)})
    seen_n, seen_o, area = n_["seen"], o_["seen"], n_["area"]
    res = {"schema": SCHEMA, "utc": utc, "new_name": a.new_name, "old_name": a.old_name,
           "new": state_numbers(Fn, An, seen_n, area, ill, Tn["q"]["H2O"]),
           "old": state_numbers(Fo, Ao, seen_o, area, To["illumination"], To["q"]["H2O"]), "checks": checks,
           "inputs": {k: str(v) for k, v in vars(a).items()}}
    say("numbers: new Q %.4e (seen %.4e, dark share %.4f), old Q %.4e (seen %.4e, dark share %.4f)" % (
        res["new"]["Q_H2O_s-1"], res["new"]["Q_H2O_seen_s-1"], res["new"]["dark_share_of_Q"], res["old"]["Q_H2O_s-1"],
        res["old"]["Q_H2O_seen_s-1"], res["old"]["dark_share_of_Q"]))
    if a.numbers_only:
        res["mode"] = "numbers only (no figure)"
        res["wall_seconds"] = time.time() - T0
        with open(out / "checks" / ("numbers_%s.json" % tag), "x") as fh:
            json.dump(res, fh, indent=1, default=str)
        with open(out / "checks" / ("numbers_%s.txt" % tag), "x") as fh:
            fh.write(numbers_text(res))
        say("numbers only: wrote %s" % (out / "checks" / ("numbers_%s.*" % tag)))
        return 0

    # ---- track (GPU: cops_cuda.figures.groundtrack_data) and images
    base_tab = FIG.read_table(a.track_base)
    window = FIG.late_window(base_tab, a.bundle, None, None)
    gt = FIG.groundtrack_data(a.track_base, None, a.bundle, window, a.device)
    if not np.array_equal(n_["normal"], gt["normals"]):
        raise SystemExit("normals differ between the track mesh and the fields")

    # ---- images (copy of illumination_maps_figures.cmd_figures, side-by-side part, with changes)
    plt = FIG.style()
    V, F, N = gt["vertices"], gt["faces"], gt["normals"]
    I_ll = M.lonlat_intensity(V, F, N)
    table, impact = gt["table"], gt["impact"]
    late = np.isin(table["iPoint"], window["points"])
    glon, glat = FIG.lonlat(gt["ground"])
    step = np.linalg.norm(np.diff(gt["ground"], axis=0), axis=1)
    ok_seg = (np.abs(np.diff(glon)) <= 180.0) & (step <= max(150.0, 5.0 * np.median(step)))
    zl, za = M.zoom_box(gt, late)
    coslat = np.cos(np.radians(0.5 * (za[0] + za[1])))
    wins = {"global": ((-180, 180), (-90, 90), V2.GLOBAL_W), "zoom": (zl, za, V2.ZOOM_W)}
    QF = [Q for Q in IM.make_quantities() if Q.key == "illumination_factor"][0]
    codes, v = QF.codes(ill, cos, f)
    qc = QF.check(codes, v, ill, cos, sh, f)
    checks["factor_classes"] = qc
    if not qc["passed"]:
        raise SystemExit("factor class check failed: %s" % qc)
    del v
    imgs, scales, own, tables = {}, {}, {}, {}
    for wname, (lon_r, lat_r, width) in wins.items():
        face, ss = V2.lonlat_faces(V, F, lon_r, lat_r, width, a.device)
        disp = V2.displayed(face)
        sel_n = disp[seen_n[disp] & (Fn[disp] > 0)]
        sel_o = disp[seen_o[disp] & (Fo[disp] > 0)]
        vals = np.concatenate((Fn[sel_n], Fo[sel_o]))
        ars = np.concatenate((area[sel_n], area[sel_o]))
        sc = V2.Scale("h2o_common_%s/%s" % (tag, wname), vals, ars, a.low_cut_area)
        scales[wname] = sc
        own[wname] = {s: V2.Scale("h2o_%s_%s/%s" % (s, tag, wname), F_[sel], area[sel], a.low_cut_area)
                      for s, F_, sel in (("new", Fn, sel_n), ("old", Fo, sel_o))}
        for s, F_, sn in (("new", Fn, seen_n), ("old", Fo, seen_o)):
            tab, _ = V2.class_table("h2o %s law %s, lon/lat %s (common scale)" % (s, tag, wname), F_, sn, area, disp, sc)
            tables["%s_%s" % (s, wname)] = tab
            imgs[(s, wname)] = V2.shade(face, V2.facet_rgb(F_, sn, sc, I_ll, a.template_relief), ss)
        imgs[("f", wname)] = V2.shade(face, IM.surface_rgb(QF, codes, seen_n, I_ll, "whole", a.template_relief), ss)
        say("window %s: %d displayed facets; common scale %s" % (wname, disp.size, sc.text()))
    res["scales"] = {w: {"common": {"edges_log10": [float(e) for e in s.edges], "step": s.step,
                                    "cut_value": s.cut_value, "vmin": s.vmin, "vmax": s.vmax},
                         "own_new": {"edges_log10": [float(e) for e in own[w]["new"].edges], "vmin": own[w]["new"].vmin,
                                     "vmax": own[w]["new"].vmax},
                         "own_old": {"edges_log10": [float(e) for e in own[w]["old"].edges], "vmin": own[w]["old"].vmin,
                                     "vmax": own[w]["old"].vmax}} for w, s in scales.items()}
    res["class_tables"] = tables

    def draw_map(ax, img, wname, c, ylabel, annotate):
        """copy of illumination_maps_figures / pub_maps_v2_figures draw_map"""
        lon_r, lat_r, _ = wins[wname]
        ax.imshow(img, extent=(lon_r[0], lon_r[1], lat_r[0], lat_r[1]), interpolation="lanczos",
                  aspect="equal" if c == 0 else 1.0 / coslat)
        V2.latlon_track_templates(ax, gt, late, glon, glat, ok_seg, impact, c, zl, za, annotate=annotate)
        V2.map_axes_setup(ax, c, lon_r, lat_r, xlabel=False, ylabel=ylabel)

    row_h, hsp = 2.05, 0.30
    rows_in = 3 * row_h + 0.05
    height = (rows_in + 3 * hsp * rows_in / 4 + 1.35) / 0.995
    fig = plt.figure(figsize=(FIG.DOUBLE_COLUMN_IN, height))
    gs = fig.add_gridspec(4, 3, width_ratios=[2.0, 1.0, 0.05], height_ratios=[row_h, row_h, 0.05, row_h],
                          wspace=0.10, hspace=hsp, left=0.07, right=0.93, top=0.995, bottom=1.35 / height)
    letters = iter("abcdef")
    axes = {}
    titles = {"new": "H₂O production rate, %s UTC, NEW law: q$_\\mathrm{min}$ / q$_\\mathrm{max}$ %g (seen facets)" % (
                  hhmmss, r_new),
              "old": "H₂O production rate, %s UTC, OLD law: q$_\\mathrm{min}$ / q$_\\mathrm{max}$ %g (seen facets)" % (
                  hhmmss, r_old),
              "f": "illumination factor %s of the new law, whole surface" % IM.QSYM}
    for r, row in ((0, "new"), (1, "old"), (3, "f")):
        for c, wname in enumerate(("global", "zoom")):
            ax = fig.add_subplot(gs[r, c])
            draw_map(ax, imgs[(row, wname)], wname, c, ylabel=(c == 0), annotate=(row == "new"))
            if c == 0:
                ax.set_title(titles[row], fontsize=7.5, loc="left", pad=3)
            if r == 3:
                ax.set_xlabel("longitude [°]")
            ax.text(0.015, 0.96, "(%s)" % next(letters), transform=ax.transAxes, va="top", fontsize=7.5,
                    bbox={"facecolor": "white", "edgecolor": "none", "pad": 0.8})
            axes[(row, wname)] = ax
        cax = fig.add_subplot(gs[r, 2])
        if row != "f":
            cax.axis("off")
        else:
            IM.fixed_colorbar(fig, cax, QF, "%s (triangle: %s)" % (IM.QSYM, IM.qr()["r"]), narrow=True,
                              orientation="vertical", fontsize=6.5)
    cbl = "H₂O production rate [m⁻² s⁻¹] (seen facets; log classes common to a-d)"
    for wname, ml in (("global", 7), ("zoom", 5)):
        V2.discrete_colorbar(fig, V2.under_axes(fig, axes[("old", wname)], 0.33, 0.11), scales[wname],
                             cbl if wname == "global" else "H₂O production rate (impact area; c, d as a, b)",
                             max_labels=ml)
    fig.legend(handles=V2.surface_handles("F_H2O") + IM.surface_handles_q(QF, "whole") + V2.track_handles_templates(gt),
               loc="lower left", bbox_to_anchor=(0.07, 0.0), ncol=2, fontsize=6.3, frameon=False, handlelength=2.0,
               columnspacing=1.2)

    def sc_txt(s, where):
        return ("%s: the seen non-zero values of both states drawn in this window range from %.3e to %.3e m^-2 s^-1; "
                "%s." % (where, s.vmin, s.vmax, s.text().replace("$10^{", "10^").replace("}$", "")) + s.below_text())

    nn, oo = res["new"], res["old"]
    mesh = "on the %s shape model (%s facets), longitude/latitude" % (Path(gt["shape"]).name, M.thousands(nf))
    cap = (
        "Inspection of the H2O source rate after the H2O night-side fix, %s UTC on 2016-09-30, %s. Left: global; right: "
        "the impact area (box in the left panels). (a, b) H2O production rate F_H2O = q_H2O A (instantaneous emitted H2O "
        "molecular flux per facet, m^-2 s^-1) of %s, computed with the NEW law q_min = %g q_max (srcCG's value; the H2O "
        "night-side fix of 2026-11-18): q_H2O = q_min on night (cos SZA <= 0) and cast-shadow facets, q_min + (q_max - "
        "q_min) cos SZA on day facets, q_max = 600 / r_h^%.7g = %.4f. (c, d) The same quantity for %s, fitted and "
        "drawn with the OLD law q_min = %g q_max (the law of every fit and figure before 2026-11-18, not intended; its "
        "full-set figures used it). a-d: colour only on the facets seen during the fit (seen F_H2O = 0 light green, unseen "
        "grey, relief-shaded); ONE set of log10 classes for a-d in each window, from the union of the two states' "
        "displayed seen non-zero values (pub_maps_v2_figures.Scale, area cut %g at the bottom, the top kept), so the "
        "colours compare directly (the old full-set figure used its own classes). %s %s (e, f) The new law's "
        "illumination factor q_H2O/q_max on the whole surface: %s on day facets in %d classes of 0.1, %s on night and "
        "cast-shadow facets (the triangle at the low end of the colour bar); the old law had 0.2 + 0.8 cos SZA and 0.2. "
        "Numbers at this time (checks/inspect_%s.txt): Q_H2O new %.3e s^-1 (seen facets %.3e, %.1f %%; night + cast "
        "shadow %.1f %% of Q_H2O, %.1f %% of the seen Q_H2O), old %.3e s^-1 (seen %.3e, %.1f %%; night + cast shadow "
        "%.1f %% of Q_H2O); new / old Q_H2O = %.3f; night + cast-shadow area %.1f %% of the surface (%.1f %% of the seen "
        "area); active area (A > 0) new %.2f %% of the surface (%.2f %% of the seen area), old %.2f %% (%.2f %%). Checks: "
        "both states' F_H2O reproduce their maps' .vtu; the two geometries give bitwise equal cos SZA and illumination "
        "classes; F_H2O = q A with q_min / q_max = %g (the maps code's H2O_QMIN_RATIO). %s" % (
            hhmmss, mesh, a.new_name, r_new, Tn["p"], Tn["q"]["H2O"][1], a.old_name, r_old, a.low_cut_area,
            sc_txt(scales["global"], "(a, c)"), sc_txt(scales["zoom"], "(b, d)"), IM.qr()["day"], IM.qr()["n_f"],
            IM.qr()["r"], tag, nn["Q_H2O_s-1"], nn["Q_H2O_seen_s-1"], 100 * nn["seen_share_of_Q"],
            100 * nn["dark_share_of_Q"], 100 * nn["dark_share_of_seen_Q"], oo["Q_H2O_s-1"], oo["Q_H2O_seen_s-1"],
            100 * oo["seen_share_of_Q"], 100 * oo["dark_share_of_Q"], nn["Q_H2O_s-1"] / oo["Q_H2O_s-1"],
            100 * nn["dark_area_fraction"], 100 * nn["dark_area_fraction_seen"], 100 * nn["active_area_fraction"],
            100 * nn["active_area_fraction_seen"], 100 * oo["active_area_fraction"], 100 * oo["active_area_fraction_seen"],
            r_new, M.track_text({"gt": gt, "track_base": str(Path(a.track_base).resolve()),
                                 "track_offset": M.gauge_offset(a.track_base),
                                 "fit_base": str(Path(a.new_base).resolve()),
                                 "fit_offset": M.gauge_offset(a.new_base)})))
    written = FIG.save(fig, prefix)
    plt.close(fig)
    with open(prefix.with_name(prefix.name + "_caption.txt"), "x") as fh:
        fh.write(cap + "\n")
    res["figure"] = written
    res["caption"] = cap
    res["wall_seconds"] = time.time() - T0
    with open(out / "checks" / ("inspect_%s.json" % tag), "x") as fh:
        json.dump(res, fh, indent=1, default=str)
    with open(out / "checks" / ("inspect_%s.txt" % tag), "x") as fh:
        fh.write(numbers_text(res))
    say("wrote %s" % written)
    return 0


if __name__ == "__main__":
    sys.exit(main())
