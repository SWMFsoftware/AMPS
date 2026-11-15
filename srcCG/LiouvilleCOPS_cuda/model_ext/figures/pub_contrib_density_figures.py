#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/pub_contrib_density_figures.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/pub_contrib_density_figures.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 0ba98e630638ff06982afdfbf1446c45390ee9ae39726368f7b3b8e8f0f0fe89
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Contribution-density maps (lon/lat, global + close-up fitted to the 90 % source areas near the impact, log
#   colour per map spanning at least 3 decades, smoothed 50/90 % source-area contours, close-up drawn-signal share) of
#   the snapshots and phase means.
# RUN: jobs/pub_contrib_density.sbatch (gpu002, 1 GPU, short)
# ==== end of header ====
"""Maps of the contribution per unit surface area of the fitted model (2026-11-14, user: "make all plots for the current
output"; lead: "maps of the contribution per unit surface area exactly as the fit computes it ... NG density per m² incl.
shield side, RG flux per m²; no spacecraft blocking; snapshots and phase means; global + impact-area close-up; NG and RG
side by side and separately; log colour, range per map; 50/90 % source-area contours").

Data: contrib_phase_data.py (<data>/snap_*.npz, phase_means.npz, manifest.json): per facet e and base row, with the fit's
kernels and arithmetic: NG n_e / A_e = sum_s F_s,e (sqrt(pi) beta_s,e / 2) K_n,e / A_e [m^-3 per m^2], RG Phi_e / A_e =
sum_s F_s,e [K_Phi,e + (V.b_RG)(sqrt(pi) beta_s,e / 2) K_V,e] / A_e [m^-2 s^-1 per m^2]; phase means = time means of these
over the phase's observation rows (late window: W10x5 of the 1-s grid rows, as the fit).
Maps: longitude / latitude (outermost surface), global and the impact area of the other figures (pub_maps_figures.zoom_box),
ground track and impact as in the template figures; colour = log10 of the per-area value on a continuous scale (viridis)
from a per-map bottom (the per-area level above which the facets carry 99 % of the map's surface signal, sum_e value x
area) to a per-map top (the 99.99th percentile of the values > 0); below the bottom: grey relief; above the top: the top
colour. Contours: the 50 % (solid) and 90 % (dashed) source areas of the map (the smallest facet sets, in decreasing
order of the facet contribution to the modelled surface pressure, that supply 50 / 90 % of it).
Figures (each .eps, .png 300 dpi, _caption.txt) in <out-dir>/contribution_density/:
  12M_contrib_density_<HHMMSS>_{NG_RG,NG,RG}   snapshots (base rows)
  12M_contrib_density_{pre_gap,late}_{NG_RG,NG,RG}, 12M_contrib_density_gap_RG   phase means
Run (GPU for the lon/lat face images, short): pub_contrib_density_figures.py --base --track-base --bundle --data
    --maps-code --out-dir [--model-name]
"""

import argparse
import datetime as dt
import json
import os
import socket
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import pub_maps_figures as M  # noqa: E402
import pub_maps_v2_figures as V2  # noqa: E402

FIG = M.FIG
T0 = time.time()
SMOOTH_PX = 5.0                    # output pixels: Gaussian smoothing of the source-area masks before the contour
SMOOTH_LEVEL = 0.3                 # contour level of the smoothed mask
MIN_DECADES = 3.0                  # the colour range spans at least this many decades below the top
UNIT = {"NG": "m⁻³ per m² of surface", "RG": "m⁻² s⁻¹ per m²"}
QNAME = {"NG": "NG: density per m² of surface [m⁻³ m⁻²]", "RG": "RG: flux into RG per m² [m⁻² s⁻¹ m⁻²]"}


def say(text):
    print("pub_contrib_density [%6.0f s]: %s" % (time.time() - T0, text), flush=True)


def source_set(c, p):
    """facets of the P source area of the per-facet contributions c (dense)."""
    o = np.argsort(-c)
    cs = np.cumsum(c[o])
    j = int(np.searchsorted(cs, p * cs[-1]))
    m = np.zeros(len(c), bool)
    m[o[:j + 1]] = True
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for a in ("--base", "--track-base", "--bundle", "--data", "--maps-code", "--out-dir"):
        ap.add_argument(a, required=True)
    ap.add_argument("--device", type=int, default=0)
    ap.add_argument("--model-name", default=M.MODEL_NAME)
    args = ap.parse_args()
    out = Path(args.out_dir)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    if (out / "contribution_density").exists():
        raise SystemExit("%s/contribution_density exists (move it aside first)" % out)
    stage = out / ("contribution_density.incomplete_%s" % stamp)
    stage.mkdir(parents=True)
    plt = FIG.style()
    import matplotlib
    from matplotlib.colors import LogNorm
    from matplotlib.cm import ScalarMappable
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    data = Path(args.data)
    man = json.load(open(data / "manifest.json"))
    base_tab = FIG.read_table(args.track_base)
    window = FIG.late_window(base_tab, args.bundle, None, None)
    gt = FIG.groundtrack_data(args.track_base, None, args.bundle, window, args.device)
    nm = M.load_module(args.maps_code, "nonneg_maps_frozen_cd")
    e1 = nm.read_elements(Path(args.base) / "elements.bin")
    if not np.array_equal(e1["normal"], gt["normals"]):
        raise SystemExit("the fit base mesh differs from the track mesh")
    area = e1["area"]
    nf = len(area)
    del e1
    V, F, N = gt["vertices"], gt["faces"], gt["normals"]
    I_ll = M.lonlat_intensity(V, F, N)
    table, impact = gt["table"], gt["impact"]
    late = np.isin(table["iPoint"], window["points"])
    glon, glat = FIG.lonlat(gt["ground"])
    step = np.linalg.norm(np.diff(gt["ground"], axis=0), axis=1)
    ok_seg = (np.abs(np.diff(glon)) <= 180.0) & (step <= max(150.0, 5.0 * np.median(step)))
    zl0, za0 = M.zoom_box(gt, late)
    face0, _ = V2.lonlat_faces(V, F, zl0, za0, V2.ZOOM_W, args.device)
    shown0 = np.zeros(len(F), bool)
    shown0[V2.displayed(face0)] = True
    del face0
    grey = V2.NOT_SEEN[None, :] * I_ll[:, None]
    cmap = matplotlib.colormaps["viridis"]
    # maps
    maps = []
    PM = np.load(data / "phase_means.npz")
    for r, s in sorted(man["snapshots"].items(), key=lambda t: t[1]["utc"]):
        z = np.load(data / s["file"])
        d = {}
        for g, k in (("NG", "ng"), ("RG", "rg")):
            v = np.zeros(nf)
            c = np.zeros(nf)
            v[z["face"]] = z[k]
            c[z["face"]] = z[k + "_Pa"]
            d[g] = (v, c)
        hh = s["utc"].split("T")[1]
        maps.append({"key": hh.replace(":", "").split(".")[0], "data": d, "gauges": ("NG", "RG"),
                     "title": "%s UTC (gauge %.1f m)" % (hh[:12], s["altitude_m"]),
                     "what": "at the base row of %s UTC (gauge altitude %.2f m; the model at that instant, no time "
                             "response)" % (hh, s["altitude_m"])})
    cl = man["classes"]
    maps.append({"key": "pre_gap", "data": {"NG": (PM["NG_pre"], PM["NG_pre_Pa"]), "RG": (PM["RG_pre"], PM["RG_pre_Pa"])},
                 "gauges": ("NG", "RG"), "title": "mean before the NG gap",
                 "what": "time mean over the observation rows before the gap in the NG data (to %s UTC; %d NG rows, %d RG "
                         "rows)" % (cl["last_NG_row_before_gap"].split("T")[1], cl["NG_pre_rows"], cl["RG_pre_rows"])})
    maps.append({"key": "gap", "data": {"RG": (PM["RG_gap"], PM["RG_gap_Pa"])}, "gauges": ("RG",),
                 "title": "mean in the NG gap (RG only)",
                 "what": "time mean over the %d RG observation rows of the gap in the NG data (after %s, before 10:33:10 "
                         "UTC; no NG data)" % (cl["RG_gap_rows"], cl["last_NG_row_before_gap"].split("T")[1])})
    maps.append({"key": "late", "data": {"NG": (PM["NG_late"], PM["NG_late_Pa"]), "RG": (PM["RG_late"], PM["RG_late_Pa"])},
                 "gauges": ("NG", "RG"), "title": "mean of the late window",
                 "what": "time mean over the %d late observations (10:33:10-10:39:26 UTC), each the W10x5 response "
                         "average of the 1-s grid rows as in the fit" % cl["late_observations"]})
    # close-up (as pub_source_area_figures.py): the 90 % source areas of all maps drawn in the standard impact area, with the
    # impact point, padded by 15 % (at least 1 deg half-size); one box for all maps
    u90 = np.zeros(nf, bool)
    for mp_ in maps:
        for g in mp_["gauges"]:
            u90 |= source_set(mp_["data"][g][1], 0.9)
    near = u90 & shown0
    ulon, ulat = FIG.lonlat(V[F[near]].mean(axis=1))
    lo = np.append(ulon, impact["lon_deg"])
    la = np.append(ulat, impact["lat_deg"])
    hl, ha = max(0.575 * np.ptp(lo), 1.0), max(0.575 * np.ptp(la), 1.0)
    cx, cy = 0.5 * (lo.min() + lo.max()), 0.5 * (la.min() + la.max())
    zl, za = (cx - hl, cx + hl), (cy - ha, cy + ha)
    n_far = int(u90.sum() - near.sum())
    coslat = np.cos(np.radians(0.5 * (za[0] + za[1])))
    wins = {"global": ((-180, 180), (-90, 90), V2.GLOBAL_W), "zoom": (zl, za, V2.ZOOM_W)}
    faces = {w: V2.lonlat_faces(V, F, lr, ar, width, args.device) for w, (lr, ar, width) in wins.items()}
    say("close-up lon %.2f..%.2f lat %.2f..%.2f; face images rendered" % (zl[0], zl[1], za[0], za[1]))
    shown = {}
    for w in wins:
        shown[w] = np.zeros(nf, bool)
        shown[w][V2.displayed(faces[w][0])] = True
    from scipy.ndimage import gaussian_filter
    mapdesc = ("Longitude/latitude maps of the %s shape model (%s facets), outermost surface where overhangs put several "
               "surface layers at one longitude/latitude: global map (top) and a close-up around the impact (bottom; "
               "lon %.1f-%.1f deg, lat %.1f-%.1f deg, box in the global map; fitted, as in the source-area figures, to the "
               "90 %% source areas of all these maps that lie in the impact area of the other figures and to the impact "
               "point; %s facets of those source areas lie farther away or under other surface layers and appear only in "
               "the global maps)." % (Path(gt["shape"]).name, M.thousands(nf), zl[0], zl[1], za[0], za[1],
                                       M.thousands(n_far)))
    captions, written, ranges = [], [], {}

    def save(fig, rel, caption):
        prefix = stage / rel
        written.extend(FIG.save(fig, prefix))
        plt.close(fig)
        written.append(FIG.write_caption(prefix, caption))
        captions.append((str(Path("contribution_density") / rel), caption))

    def layout(gauges):
        ncol = len(gauges)
        W = FIG.DOUBLE_COLUMN_IN if ncol == 2 else FIG.DOUBLE_COLUMN_IN * 0.62
        fig = plt.figure(figsize=(W, 6.6))
        axes = {}
        for j, g in enumerate(gauges):
            x0 = 0.09 + j * 0.47 if ncol == 2 else 0.14
            wd = 0.42 if ncol == 2 else 0.82
            axes[(g, "global")] = fig.add_axes([x0, 0.70, wd, 0.26])
            axes[(g, "zoom")] = fig.add_axes([x0, 0.29, wd, 0.35])
            axes[(g, "cbar")] = fig.add_axes([x0, 0.215, wd, 0.016])
        return fig, axes

    trk = ("Black line: ground track of Rosetta (sub-spacecraft intercept), thick in the late-visibility window %s-%s UTC "
           "with 1-minute marks; square: track start; star: Rosetta impact (%s UTC, lon %.3f°, lat %.3f°)." % (
               window["start"].split("T")[1], window["stop"].split("T")[1], impact["utc"].split("T")[1],
               impact["lon_deg"], impact["lat_deg"]))
    defs = ("Colour: the contribution per unit surface area as the fitted model (%s) computes it, with the fit's kernels at "
            "the gauge position (base rows of the fitting base) and its arithmetic (fitted H2O activity, CO2 scale, "
            "temperature law, night and cast-shadow rule; H2O + CO2): NG = density at the nude gauge F (√π β/2) K_n / A "
            "[m⁻³ per m²], the fit's NG acceptance (front factor 1 − l·b_NG and the shield fraction behind) inside K_n, so "
            "the back side through the shield is included; RG = flux into the ram gauge F (K_Φ + (V·b_RG)(√π β/2) K_V) / "
            "A [m⁻² s⁻¹ per m²] over its whole front half-space (the fit's acceptance); no spacecraft blocking (as the "
            "fit); facets the gauge does not see contribute 0. Log scale from the per-map bottom (the level above which "
            "the facets carry 99 %% of the map's surface signal, or 3 decades below the top if that is lower) to the "
            "per-map top (the 99.99th percentile of the values > 0); grey: relief-shaded surface below the bottom or without contribution; above the top: the top "
            "colour. White contours: the 50 %% (solid) and 90 %% (dashed) source areas of the map (the smallest facet "
            "sets, in decreasing order of the facet contribution to the modelled surface pressure, that supply 50 / 90 %% "
            "of it; drawn as the %.1f level of the set's pixel mask after a Gaussian smoothing of %.0f pixels, so the "
            "outline encloses clusters of facets rather than single facets). Facets on overhangs and under other "
            "surface layers are not drawn in the longitude/latitude maps; the share of each map's signal from the facets "
            "drawn in the close-up is given below. Data: %s (contrib_phase_data.py)." % (args.model_name, SMOOTH_LEVEL, SMOOTH_PX, data))
    for mp_ in maps:
        rng = {}
        rgbs = {}
        sets = {}
        for g in mp_["gauges"]:
            v, c = mp_["data"][g]
            sig = v * area
            o = np.argsort(-v)
            cs = np.cumsum(sig[o])
            j = int(np.searchsorted(cs, 0.99 * cs[-1]))
            bottom = float(v[o][min(j, nf - 1)])
            top = float(np.percentile(v[v > 0], 99.99))
            if not top > bottom:
                top = float(v.max())
            bottom = min(bottom, top / 10 ** MIN_DECADES)
            norm = LogNorm(vmin=bottom, vmax=top)
            rgb = grey.copy()
            m = v > bottom
            rgb[m] = cmap(norm(np.clip(v[m], bottom, top)))[:, :3]
            rgbs[g] = (rgb, norm)
            sets[g] = (source_set(c, 0.5), source_set(c, 0.9))
            rng[g] = {"bottom": bottom, "top": top, "decades": float(np.log10(top / bottom)), "unit": UNIT[g],
                      "facets_coloured": int(m.sum()), "P50_facets": int(sets[g][0].sum()),
                      "P90_facets": int(sets[g][1].sum()), "P50_area_m2": float(area[sets[g][0]].sum()),
                      "P90_area_m2": float(area[sets[g][1]].sum()), "signal_total": float(sig.sum()),
                      "signal_share_drawn_global": float(sig[shown["global"]].sum() / sig.sum()),
                      "signal_share_drawn_zoom": float(sig[shown["zoom"]].sum() / sig.sum())}
        ranges[mp_["key"]] = rng
        variants = [mp_["gauges"]] + ([(g,) for g in mp_["gauges"]] if len(mp_["gauges"]) == 2 else [])
        for gauges in variants:
            fig, axes = layout(gauges)
            for j, g in enumerate(gauges):
                rgb, norm = rgbs[g]
                for w, cc in (("global", 0), ("zoom", 1)):
                    lr, ar, _ = wins[w]
                    face, ss = faces[w]
                    ax = axes[(g, w)]
                    ax.imshow(V2.shade(face, rgb, ss), extent=(lr[0], lr[1], ar[0], ar[1]), interpolation="lanczos",
                              aspect="equal" if cc == 0 else 1.0 / coslat)
                    for k, ls in ((0, "-"), (1, "--")):
                        msk = np.where(face >= 0, sets[g][k][np.maximum(face, 0)], False).astype(float)
                        Hs, Ws = msk.shape
                        msk = msk.reshape(Hs // ss, ss, Ws // ss, ss).mean(axis=(1, 3))
                        msk = gaussian_filter(msk, SMOOTH_PX)    # outline of the facet set (not of single facets)
                        xs = np.linspace(lr[0], lr[1], msk.shape[1])
                        ys = np.linspace(ar[1], ar[0], msk.shape[0])
                        if msk.max() > SMOOTH_LEVEL > msk.min():
                            ax.contour(xs, ys, msk, levels=[SMOOTH_LEVEL], colors="white", linewidths=0.55 if cc else 0.4,
                                       linestyles=ls, zorder=4)
                    V2.latlon_track_templates(ax, gt, late, glon, glat, ok_seg, impact, cc, zl, za, annotate=False)
                    V2.map_axes_setup(ax, cc, lr, ar, xlabel=True, ylabel=(j == 0))
                axes[(g, "global")].set_title("%s, %s" % (g, mp_["title"]), fontsize=8, loc="left", pad=3)
                cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=axes[(g, "cbar")], orientation="horizontal")
                cb.ax.tick_params(labelsize=6.0)
                cb.set_label(QNAME[g], fontsize=6.5)
            hs = [Line2D([], [], color="0.25", lw=0.8, ls="-", label="50 % source area"),
                  Line2D([], [], color="0.25", lw=0.8, ls="--", label="90 % source area"),
                  Patch(facecolor=V2.NOT_SEEN * 0.85, edgecolor="none", label="below the colour range / none")]
            fig.legend(handles=hs + V2.track_handles_templates(gt), loc="lower left", bbox_to_anchor=(0.02, 0.0),
                       ncol=3 if len(gauges) == 2 else 2, fontsize=6.0, frameon=False, handlelength=1.8, columnspacing=1.0)
            name = "12M_contrib_density_%s_%s" % (mp_["key"], "_".join(gauges))
            cap = ("Contribution per unit surface area of the fitted model to %s, %s. %s %s Ranges: %s. Source areas: "
                   "%s. Share of the map's signal from the facets drawn in the close-up (facets outside it, hidden under other "
                   "surface layers such as overhangs, or smaller than a pixel are not drawn): %s. %s" % (
                       " and ".join({"NG": "the gas density at NG", "RG": "the flux into RG"}[g] for g in gauges),
                       mp_["what"], defs, mapdesc,
                       "; ".join("%s %.3g-%.3g %s (%.1f decades)" % (g, rng[g]["bottom"], rng[g]["top"], UNIT[g],
                                                                       rng[g]["decades"]) for g in gauges),
                       "; ".join("%s P50 %s facets (%.3g m²), P90 %s facets (%.3g m²)" % (
                           g, M.thousands(rng[g]["P50_facets"]), rng[g]["P50_area_m2"], M.thousands(rng[g]["P90_facets"]),
                           rng[g]["P90_area_m2"]) for g in gauges),
                       "; ".join("%s %.0f %%" % (g, 100 * rng[g]["signal_share_drawn_zoom"]) for g in gauges), trk))
            save(fig, name, cap)
        say("map %s written" % mp_["key"])
    json.dump({"created_utc": stamp, "host": socket.gethostname(), "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
               "command": sys.argv, "script_sha256": M.sha256_file(__file__), "data": str(data), "ranges": ranges,
               "captions": [{"figure": n, "caption": c} for n, c in captions],
               "files": [str(Path(p).relative_to(stage)) for p in written]},
              open(stage / "contribution_density_manifest.json", "w"), indent=1, ensure_ascii=False)
    stage.rename(out / "contribution_density")
    say("done: %d files" % (len(written) + 1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
