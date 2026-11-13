#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/pub_source_area_figures.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/pub_source_area_figures.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 5c4ca96228e3eab70fce008cdbdfcfc3a348dffdebe4d9cdc30bbf44f4c8460b
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure script: P% source-area union and occurrence maps of the 95 post-gap observations (NG, RG, NG_RG);
#   2026-11-13: --model-name, the contribution job and verification from contrib_manifest.json.
# RUN: jobs/pub_source_areas.sbatch STEPS=sources (GPU)
# ==== end of header ====
"""P% source areas of the post-gap COPS observations (softplus model, final state), NG and RG: union maps and
source-area occurrence maps.

User request 2026-11-12 (second round): per gauge (NG, RG) and P = 50, 60, 70, 80, 90 %: the union of the per-observation
highest-contribution facet sets and the per-facet count of the post-gap observations whose set contains the facet;
"put NG and RG plots side by side in the same figure as well as save them separately". Names (proposed by the maps
agent, confirmed by the lead): "P% source area" S_P(k, g) = the smallest set of facets, taken in decreasing order of
their contribution, that supplies P % of observation k's modelled surface signal at gauge g (after the footprint
'source area' of Schmid 1994); "source-area occurrence" N_P,g(e) = the number of the 95 post-gap observations
(10:33:10-10:39:26 UTC) whose S_P contains facet e. Data: <chain>/contributions/source_areas_{NG,RG}.npz
(model_ext nonneg_contrib.py, job 124549; README.txt there).

FIGURES (lon/lat, outermost surface; global map + the impact area of the other figures; ground track, late-window track,
1-min marks and the Rosetta impact star exactly as in the template figures; grey relief elsewhere):
  12M_source_area_union_{NG,RG,NG_RG}          nested unions U_P = {N_P > 0}, P = 50 (darkest) .. 90 (lightest);
                                               the S_P of one observation are nested (S_50 c S_60 c ...), so are U_P
  12M_source_area_occurrence_P<P>_{NG,RG,NG_RG}  N_P in classes 1, 2-4, 5-9, 10-19, 20-39, 40-59, 60-79, 80-94, 95
Every figure: .eps, .png (300 dpi), _caption.txt. checks/source_area_tables.txt: facets and area per class (counted from
the occurrence arrays; the union counts are checked against the files' union arrays and the nesting is checked).

Run from <AMPS>/srcCG/LiouvilleCOPS_cuda with .venv/bin/python on a GPU node (jobs/pub_source_areas.sbatch):
    pub_source_area_figures.py --base <fit base> --track-base <observation-only base> --bundle <bundle .npz>
        --contrib-dir <chain>/contributions --maps-code <frozen nonneg_maps.py> --out-dir <v2 folder>
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
SCHEMA = "pub_source_area_v1"
OCC_EDGES = (1, 2, 5, 10, 20, 40, 60, 80, 95, 96)
OCC_LABELS = ("1", "2–4", "5–9", "10–19", "20–39", "40–59", "60–79", "80–94", "95")
T0 = time.time()


def say(text):
    print("pub_source_areas [%6.0f s]: %s" % (time.time() - T0, text), flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for a in ("--base", "--track-base", "--bundle", "--contrib-dir", "--maps-code", "--out-dir"):
        ap.add_argument(a, required=True)
    ap.add_argument("--device", type=int, default=0)
    ap.add_argument("--model-name", default=M.MODEL_NAME,
                    help="name of the fitted state in the captions (default: %(default)s; 2026-11-13: set per fit)")
    args = ap.parse_args()
    out = Path(args.out_dir)
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    if (out / "source_areas").exists():
        raise SystemExit("%s/source_areas exists (move it aside first)" % out)
    stage = out / ("source_areas.incomplete_%s" % stamp)
    stage.mkdir(parents=True)
    plt = FIG.style()
    import matplotlib
    from matplotlib.colors import BoundaryNorm, ListedColormap
    from matplotlib.cm import ScalarMappable
    from matplotlib.patches import Patch

    base_tab = FIG.read_table(args.track_base)
    window = FIG.late_window(base_tab, args.bundle, None, None)
    gt = FIG.groundtrack_data(args.track_base, None, args.bundle, window, args.device)
    say("track: %d rows, impact %s UTC" % (len(gt["table"]["utc"]), gt["impact"]["utc"]))
    nm = M.load_module(args.maps_code, "nonneg_maps_frozen_sa")
    e1 = nm.read_elements(Path(args.base) / "elements.bin")
    if not np.array_equal(e1["normal"], gt["normals"]):
        raise SystemExit("the fit base mesh differs from the track mesh")
    area = e1["area"]
    del e1
    cdir = Path(args.contrib_dir)
    # job id and verification of the contribution data from their manifest (data-driven, 2026-11-13)
    cman = json.load(open(cdir / "contrib_manifest.json"))
    cjob = cman.get("slurm_job_id") or "?"
    ver = cman["verification"]
    vmax = max(ver["NG final"]["max_rel"], ver["RG final"]["max_rel"])
    say("contributions: job %s, per-row sums vs the model tables max rel %.3g" % (cjob, vmax))
    nneg = cman["negative_record_contributions"]
    if any(int(v) for v in nneg.values()):
        raise SystemExit("negative record contributions in %s: %s (the caption states all >= 0)" % (cdir, nneg))
    data, checks = {}, {}
    for g in ("NG", "RG"):
        with np.load(cdir / ("source_areas_%s.npz" % g), allow_pickle=False) as z:
            occ = np.asarray(z["occurrence"])
            uni = np.asarray(z["union"])
            pl = [int(p) for p in z["p_levels"]]
            grain_share = np.asarray(z["grain_Pa"]) / (np.asarray(z["model_surface_Pa"]) + np.asarray(z["grain_Pa"]))
            kk, aa = np.asarray(z["k"]), np.asarray(z["area_m2"])
            nobs = len(z["iPoint"])
        if occ.shape != (len(pl), len(area)):
            raise SystemExit("occurrence of %s has shape %s" % (g, occ.shape))
        c = {"union_equals_file": bool(np.array_equal(occ > 0, uni)),
             "nested": bool(all(np.all(~uni[i] | uni[i + 1]) for i in range(len(pl) - 1))),
             "max_occurrence": [int(v) for v in occ.max(axis=1)], "observations": nobs}
        if not (c["union_equals_file"] and c["nested"]):
            raise SystemExit("source-area check failed for %s: %s" % (g, c))
        checks[g] = c
        data[g] = {"occ": occ, "uni": uni, "pl": pl, "grain_share": grain_share, "k": kk, "a": aa, "nobs": nobs}
        say("%s: union facets per P %s; checks %s" % (g, [int(u.sum()) for u in uni], c))
    pl = data["NG"]["pl"]

    V, F, N = gt["vertices"], gt["faces"], gt["normals"]
    I_ll = M.lonlat_intensity(V, F, N)
    table, impact = gt["table"], gt["impact"]
    late = np.isin(table["iPoint"], window["points"])
    glon, glat = FIG.lonlat(gt["ground"])
    step = np.linalg.norm(np.diff(gt["ground"], axis=0), axis=1)
    ok_seg = (np.abs(np.diff(glon)) <= 180.0) & (step <= max(150.0, 5.0 * np.median(step)))
    zl0, za0 = M.zoom_box(gt, late)
    # close-up of the source areas: the P90 unions of both gauges inside the standard impact box, with the impact point,
    # padded by 15 % (at least 1.5 deg half-size); farther parts of the unions are visible in the global map
    # (only the union facets DRAWN in the impact-area image count: facets hidden under overhangs in the outermost-
    # surface map would stretch the box)
    u90 = data["NG"]["uni"][-1] | data["RG"]["uni"][-1]
    face0, _ = V2.lonlat_faces(V, F, zl0, za0, V2.ZOOM_W, args.device)
    shown = np.zeros(len(F), bool)
    shown[V2.displayed(face0)] = True
    del face0
    near = u90 & shown
    ulon, ulat = FIG.lonlat(V[F[near]].mean(axis=1))
    lo = np.append(ulon, impact["lon_deg"])
    la = np.append(ulat, impact["lat_deg"])
    hl, ha = max(0.575 * np.ptp(lo), 1.0), max(0.575 * np.ptp(la), 1.0)
    cx, cy = 0.5 * (lo.min() + lo.max()), 0.5 * (la.min() + la.max())
    zl, za = (cx - hl, cx + hl), (cy - ha, cy + ha)
    n_far = int(u90.sum() - near.sum())
    say("close-up box lon %.2f..%.2f lat %.2f..%.2f; P90 union facets not drawn in the impact area: %d" % (
        zl[0], zl[1], za[0], za[1], n_far))
    coslat = np.cos(np.radians(0.5 * (za[0] + za[1])))
    wins = {"global": ((-180, 180), (-90, 90), V2.GLOBAL_W), "zoom": (zl, za, V2.ZOOM_W)}
    faces = {}
    for w, (lr, ar, width) in wins.items():
        faces[w] = V2.lonlat_faces(V, F, lr, ar, width, args.device)
    say("face images rendered")
    grey = V2.NOT_SEEN[None, :] * I_ll[:, None]
    # nested union colours: P = 90 lightest ... 50 darkest (Oranges)
    ucols = matplotlib.colormaps["Oranges"](np.linspace(0.30, 0.95, len(pl)))[::-1][:, :3]   # index 0 = P50 (darkest)
    ocols = matplotlib.colormaps["YlOrRd"](np.linspace(0.12, 1.0, len(OCC_LABELS)))[:, :3]
    tables_txt, captions, written = [], [], []

    def save(fig, rel, caption):
        prefix = stage / rel
        written.extend(FIG.save(fig, prefix))
        plt.close(fig)
        written.append(FIG.write_caption(prefix, caption))
        captions.append((str(Path("source_areas") / rel), caption))

    def union_rgb(g):
        rgb = grey.copy()
        uni = data[g]["uni"]
        for i in range(len(pl) - 1, -1, -1):          # paint P90 first, then the nested smaller unions on top
            rgb[uni[i]] = ucols[i]
        return rgb

    def occ_rgb(g, i):
        rgb = grey.copy()
        n = data[g]["occ"][i]
        m = n > 0
        cls = np.searchsorted(np.asarray(OCC_EDGES), n[m], side="right") - 1
        rgb[m] = ocols[cls]
        return rgb, n

    def draw(ax, rgb, w, c, ylabel=True):
        lr, ar, _ = wins[w]
        face, ss = faces[w]
        ax.imshow(V2.shade(face, rgb, ss), extent=(lr[0], lr[1], ar[0], ar[1]), interpolation="lanczos",
                  aspect="equal" if c == 0 else 1.0 / coslat)
        V2.latlon_track_templates(ax, gt, late, glon, glat, ok_seg, impact, c, zl, za, annotate=False)
        V2.map_axes_setup(ax, c, lr, ar, xlabel=True, ylabel=ylabel)
        ax.text(0.015, 0.96, "", transform=ax.transAxes)

    def layout(gauges):
        """Figure with one column per gauge: global map above, impact area below."""
        ncol = len(gauges)
        W = FIG.DOUBLE_COLUMN_IN if ncol == 2 else FIG.DOUBLE_COLUMN_IN * 0.62
        H = 6.3
        fig = plt.figure(figsize=(W, H))
        axes = {}
        for j, g in enumerate(gauges):
            x0 = 0.09 + j * 0.47 if ncol == 2 else 0.14
            wd = 0.42 if ncol == 2 else 0.82
            axes[(g, "global")] = fig.add_axes([x0, 0.69, wd, 0.27])
            axes[(g, "zoom")] = fig.add_axes([x0, 0.25, wd, 0.36])
        return fig, axes

    track_h = V2.track_handles_templates(gt)
    defs = ("P%% source area S_P(k, g) of a post-gap observation k at gauge g: the smallest set of facets, taken in "
            "decreasing order of their contribution c_k,g,e to the observation's modelled surface pressure (the fitted "
            "forward model of %s, term by "
            "term; for the 95 observations of 10:33:10-10:39:26 UTC averaged over the gauges' 10-s and 5-readout "
            "response like the fit), that supplies P %% of sum_e c_k,g,e (after the footprint 'source area' of Schmid "
            "1994). Source-area occurrence N_P,g(e): the number of the 95 post-gap observations whose S_P contains "
            "facet e (0-95). All contributions are >= 0; the facet sums reproduce the model tables to <= %.2g "
            "(print precision); the NG grain term is not a facet contribution and is excluded (its share of the model "
            "NG pressure is %.2f-%.2f, median %.2f, in these observations). Data: %s (job %s)." % (
                args.model_name, vmax, data["NG"]["grain_share"].min(), data["NG"]["grain_share"].max(),
                float(np.median(data["NG"]["grain_share"])), cdir, cjob))
    mapdesc = ("Longitude/latitude maps of the %s shape model (%s facets), outermost surface where overhangs put several "
               "surface layers at one longitude/latitude: global map (top) and a close-up of the source areas around "
               "the impact (bottom; box in the global map: lon %.1f-%.1f deg, lat %.1f-%.1f deg, fitted to the union "
               "facets drawn in the impact area; %s facets of the P90 unions lie outside the impact area of the other "
               "figures or are hidden under other surface layers in the lon/lat projection). Grey: "
               "relief-shaded surface. " % (Path(gt["shape"]).name, M.thousands(len(F)), zl[0], zl[1], za[0], za[1],
                                            M.thousands(n_far)))
    trk = ("Black line: ground track of Rosetta (sub-spacecraft intercept), thick in the late-visibility window "
           "%s-%s UTC with 1-minute marks; square: track start; star: Rosetta impact (%s UTC, lon %.3f°, lat %.3f°); "
           "track and impact from the observation rows (%s)." % (
               window["start"].split("T")[1], window["stop"].split("T")[1], impact["utc"].split("T")[1],
               impact["lon_deg"], impact["lat_deg"], args.track_base))

    def stats_text(g):
        u = data[g]["uni"]
        return "; ".join("U_%d: %s facets, %.3g m^2" % (p, M.thousands(u[i].sum()), area[u[i]].sum())
                         for i, p in enumerate(pl))

    # ---- union maps
    for gauges in (("NG",), ("RG",), ("NG", "RG")):
        fig, axes = layout(gauges)
        for j, g in enumerate(gauges):
            rgb = union_rgb(g)
            for w, c in (("global", 0), ("zoom", 1)):
                draw(axes[(g, w)], rgb, w, c, ylabel=(j == 0))
            axes[(g, "global")].set_title("%s: union of the P%% source areas" % g, fontsize=8, loc="left", pad=3)
        hs = [Patch(facecolor=ucols[i], edgecolor="none", label="in a %d%% source area" % p) for i, p in enumerate(pl)]
        fig.legend(handles=hs + [Patch(facecolor=V2.NOT_SEEN * 0.85, edgecolor="none", label="in none")] + track_h,
                   loc="lower left", bbox_to_anchor=(0.02, 0.0), ncol=3 if len(gauges) == 2 else 2, fontsize=6.0,
                   frameon=False, handlelength=1.8, columnspacing=1.0)
        name = "12M_source_area_union_%s" % "_".join(gauges)
        cap = ("Union of the P%% source areas of the 95 post-gap COPS observations, %s: U_P = the facets that belong to "
               "S_P of at least one observation (N_P > 0), for P = 50 (darkest) to 90 %% (lightest); the unions are "
               "nested (S_50 within S_60 ... within S_90 for every observation), so each facet is coloured by the "
               "smallest P whose union contains it. %s%s %s. %s" % (
                   " and ".join(("%s (%s)" % (g, "left" if j == 0 and len(gauges) == 2 else "right" if j else ""))
                                .replace(" ()", "") for j, g in enumerate(gauges)),
                   mapdesc, defs, "; ".join("%s: %s" % (g, stats_text(g)) for g in gauges), trk))
        save(fig, name, cap)
    say("union maps written")

    # ---- occurrence maps per P
    cm = ListedColormap(ocols)
    norm = BoundaryNorm(np.arange(len(OCC_LABELS) + 1), cm.N)
    for i, p in enumerate(pl):
        rgbs = {g: occ_rgb(g, i) for g in ("NG", "RG")}
        for g in ("NG", "RG"):
            n = rgbs[g][1]
            m = n > 0
            cls = np.searchsorted(np.asarray(OCC_EDGES), n[m], side="right") - 1
            L = ["[%s, P = %d %%] facets with N > 0: %d (= union: %s), max N %d; area %.6e m^2" % (
                g, p, m.sum(), bool(m.sum() == data[g]["uni"][i].sum()), n.max(), area[m].sum())]
            for j, lab in enumerate(OCC_LABELS):
                s = cls == j
                L.append("  N %-6s colour %s facets %8d area %.6e m^2" % (lab, V2.hexcol(ocols[j]), s.sum(),
                                                                          area[m][s].sum()))
            tables_txt.append("\n".join(L))
        for gauges in (("NG",), ("RG",), ("NG", "RG")):
            fig, axes = layout(gauges)
            for j, g in enumerate(gauges):
                for w, c in (("global", 0), ("zoom", 1)):
                    draw(axes[(g, w)], rgbs[g][0], w, c, ylabel=(j == 0))
                axes[(g, "global")].set_title("%s: occurrence in the %d%% source area" % (g, p), fontsize=8,
                                              loc="left", pad=3)
            pos = axes[(gauges[-1], "zoom")].get_position()
            cax = fig.add_axes([0.55 if len(gauges) == 2 else 0.14, 0.135, 0.40 if len(gauges) == 2 else 0.8, 0.018])
            cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cm), cax=cax, orientation="horizontal",
                              ticks=np.arange(len(OCC_LABELS)) + 0.5)
            cb.ax.set_xticklabels(OCC_LABELS)
            cb.ax.tick_params(labelsize=6.0, length=0)
            cb.set_label("N_%d: post-gap observations (of 95) whose %d%% source area contains the facet" % (p, p),
                         fontsize=6.3)
            fig.legend(handles=[Patch(facecolor=V2.NOT_SEEN * 0.85, edgecolor="none", label="N = 0")] + track_h,
                       loc="lower left", bbox_to_anchor=(0.02, 0.0), ncol=2 if len(gauges) == 1 else 3, fontsize=6.0,
                       frameon=False, handlelength=1.8, columnspacing=1.0)
            name = "12M_source_area_occurrence_P%d_%s" % (p, "_".join(gauges))
            cap = ("Source-area occurrence N_%d of the facets, %s: the number of the 95 post-gap COPS observations "
                   "whose %d%% source area contains the facet, in the classes %s (colour bar); grey: N = 0. %s%s "
                   "Facets with N_%d > 0: %s. %s" % (
                       p, " and ".join(gauges), p, ", ".join(OCC_LABELS), mapdesc, defs, p,
                       "; ".join("%s %s (max N %d)" % (g, M.thousands((data[g]["occ"][i] > 0).sum()),
                                                       data[g]["occ"][i].max()) for g in gauges), trk))
            save(fig, name, cap)
        say("occurrence maps P = %d written" % p)

    (stage / "checks").mkdir()
    with open(stage / "checks" / "source_area_tables.txt", "w") as f:
        f.write("# pub_source_area_figures.py (%s): facets and area per occurrence class (counted from the "
                "source_areas_{NG,RG}.npz occurrence arrays); union = N > 0, checked against the files' union arrays; "
                "nesting U_50 c U_60 c ... checked: %s\n\n" % (SCHEMA, json.dumps(checks)))
        f.write("\n\n".join(tables_txt) + "\n")
    with open(stage / "checks" / "source_area_captions.json", "w") as f:
        json.dump([{"figure": n, "caption": c} for n, c in captions], f, indent=1, ensure_ascii=False)
    json.dump({"schema": SCHEMA, "created_utc": stamp, "host": socket.gethostname(),
               "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "command": sys.argv,
               "script_sha256": M.sha256_file(__file__), "contrib_dir": str(cdir), "checks": checks,
               "files": [str(Path(p).relative_to(stage)) for p in written]},
              open(stage / "checks" / "source_area_manifest.json", "w"), indent=1)
    # move into place: figures to out/source_areas/, checks to out/checks/
    (out / "checks").mkdir(exist_ok=True)
    for p in sorted((stage / "checks").iterdir()):
        dest = out / "checks" / p.name
        if dest.exists():
            raise SystemExit("%s exists" % dest)
        p.rename(dest)
    (stage / "checks").rmdir()
    stage.rename(out / "source_areas")
    say("done: %d files" % (len(written) + 3))
    return 0


if __name__ == "__main__":
    sys.exit(main())
