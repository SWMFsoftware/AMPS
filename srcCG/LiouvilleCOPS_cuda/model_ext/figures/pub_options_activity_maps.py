#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/pub_options_activity_maps.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/pub_options_activity_maps.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 9ce74fcc5ad2d389f3a1e000a27bb720a7108731f27438906107cf5fed8d8d5f
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: activity maps (log A on the seen surface, global + impact close-up, one log scale, v2 colours) of the best
#   stage-1 setting of each non-negative option and the bounded softplus c9 (2026-11-14)
# RUN: inside jobs/options_activity_maps.sbatch (1 GPU): pub_options_activity_maps.py --folder <comparison folder>
#   --base --track-base --bundle --maps-code --c9-ft --c9-fit
# ==== end of header ====
"""Activity maps of the best stage-1 setting of each non-negative option and the bounded softplus c9, side by side
(2026-11-14, lead: "log A on the seen surface, global + impact close-up; same style as the activity v2 maps ... to show how
much surface contributes").

A per facet: the fit's nonneg_facets.npz A (options; the fit's activity incl. the multipliers and the held fine-tuning term)
and, for the softplus c9 final state, nonneg_maps.state_activity of c9/ft/final (the activity of the maps and the
contribution data). Seen: the fit's seen flag (nonneg_facets.npz seen; c9: the state's seen). Colour: log10 A on the seen
facets with A above the active threshold (1 % of the median of the srcCG activity A_orig over the seen facets, the
summaries' definition), one continuous scale for all panels up to the largest 99.9th percentile of the seen A; seen facets
below the threshold light grey; facets never seen by a fitted observation darker grey (relief). Longitude/latitude,
outermost surface; global map and the impact area of the other figures; ground track and impact as the template figures.
Writes <folder>/activity_maps/12M_options_activity_best.{eps,png} + _caption.txt and activity_maps_manifest.json.
Run (GPU for the lon/lat face images, short): pub_options_activity_maps.py --folder <comparison folder> --base <fit base>
    --track-base <observation-only base> --bundle <bundle> --maps-code <frozen nonneg_maps.py> --c9-ft --c9-fit
"""
import argparse
import csv
import json
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import pub_maps_figures as M  # noqa: E402
import pub_maps_v2_figures as V2  # noqa: E402

FIG = M.FIG
NAMES = {"A": "A: sum of squares", "B": "B: background", "C": "C: local basis", "D": "D: log spread"}
T0 = time.time()


def say(t):
    print("options_activity [%5.0f s]: %s" % (time.time() - T0, t), flush=True)


def main():
    ap = argparse.ArgumentParser()
    for k in ("--folder", "--base", "--track-base", "--bundle", "--maps-code", "--c9-ft", "--c9-fit"):
        ap.add_argument(k, required=True)
    ap.add_argument("--device", type=int, default=0)
    ap.add_argument("--procs", type=int, default=8)
    a = ap.parse_args()
    O = Path(a.folder)
    out = O / "activity_maps"
    out.mkdir(exist_ok=False)
    plt = FIG.style()
    import matplotlib
    from matplotlib.colors import LogNorm
    from matplotlib.cm import ScalarMappable
    from matplotlib.patches import Patch
    runs = json.load(open(O / "runs.json"))
    best = json.load(open(O / "best.json"))
    rows = {(r["option"], r["setting"]): r for r in csv.DictReader(open(O / "summary.csv"))}
    nm = M.load_module(a.maps_code, "nonneg_maps_frozen_oa")
    el = nm.read_elements(Path(a.base) / "elements.bin")
    area = el["area"]
    nf = len(area)
    models = []
    for opt, b in sorted(best.items()):
        r = [x for x in runs["runs"] if x["option"] == opt and x["setting"] == b["setting"]][0]
        with np.load(Path(r["root"]) / "fit" / "final" / "nonneg_facets.npz", allow_pickle=False) as z:
            A, seen = np.asarray(z["A"], np.float64), np.asarray(z["seen"]).astype(bool)
        models.append({"label": "%s, %s (data cost %.4e)" % (NAMES[opt], b["setting"], b["data_cost"]), "A": A, "seen": seen,
                       "row": rows[(opt, b["setting"])], "src": str(Path(r["root"]) / "fit/final/nonneg_facets.npz")})
    S = nm.state_activity(Path(a.c9_ft), Path(a.c9_fit), el, a.procs)
    if not S["checks"]["passed"]:
        raise SystemExit("c9 state checks failed")
    models.append({"label": "bounded softplus c9 (data cost %.4e)" % float(rows[("ref", "9 cycles")]["data_cost"]),
                   "A": np.asarray(S["A"], np.float64), "seen": np.asarray(S["seen"]).astype(bool),
                   "row": rows[("ref", "9 cycles")], "src": "%s (nonneg_maps.state_activity)" % a.c9_ft})
    say("activities read: %s" % [m["label"] for m in models])
    med = float(np.median(el["AH2O"][models[-1]["seen"]]))
    thr = 0.01 * med
    top = max(float(np.percentile(m["A"][m["seen"] & (m["A"] > thr)], 99.9)) for m in models if (m["seen"] & (m["A"] > thr)).any())
    bottom = thr
    norm = LogNorm(vmin=bottom, vmax=top)
    cmap = matplotlib.colormaps[V2.CMAP_NAME]
    base_tab = FIG.read_table(a.track_base)
    window = FIG.late_window(base_tab, a.bundle, None, None)
    gt = FIG.groundtrack_data(a.track_base, None, a.bundle, window, a.device)
    V, F, N = gt["vertices"], gt["faces"], gt["normals"]
    if not np.array_equal(el["normal"], N):
        raise SystemExit("the fit base mesh differs from the track mesh")
    I_ll = M.lonlat_intensity(V, F, N)
    table, impact = gt["table"], gt["impact"]
    late = np.isin(table["iPoint"], window["points"])
    glon, glat = FIG.lonlat(gt["ground"])
    step = np.linalg.norm(np.diff(gt["ground"], axis=0), axis=1)
    ok_seg = (np.abs(np.diff(glon)) <= 180.0) & (step <= max(150.0, 5.0 * np.median(step)))
    zl, za = M.zoom_box(gt, late)
    coslat = np.cos(np.radians(0.5 * (za[0] + za[1])))
    wins = {"global": ((-180, 180), (-90, 90), V2.GLOBAL_W), "zoom": (zl, za, V2.ZOOM_W)}
    faces = {w: V2.lonlat_faces(V, F, lr, ar, width, a.device) for w, (lr, ar, width) in wins.items()}
    say("face images rendered")
    nm_ = len(models)
    BM, PITCH, RH, TM = 1.45, 1.65, 1.40, 0.35                # inches: bottom margin, row pitch, row height, top margin
    H = BM + (nm_ - 1) * PITCH + RH + TM
    fig = plt.figure(figsize=(FIG.DOUBLE_COLUMN_IN, H))
    stats = []
    for i, m in enumerate(models):
        A, seen = m["A"], m["seen"]
        act = seen & (A > thr)
        rgb = V2.NOT_SEEN[None, :] * I_ll[:, None]                       # as pub_maps_v2_figures.facet_rgb
        low = seen & ~act
        rgb[low] = V2.ZERO_SEEN[None, :] * I_ll[low, None]
        rgb[act] = cmap(norm(np.clip(A[act], bottom, top)))[:, :3]
        rgb = np.clip(rgb, 0.0, 1.0)
        st = {"label": m["label"], "active_seen_area_frac_recomputed": float(area[act].sum() / area[seen].sum()),
              "active_seen_area_frac_summary": m["row"].get("active_seen_area_frac"), "source": m["src"]}
        stats.append(st)
        y0 = (BM + (nm_ - 1 - i) * PITCH) / H
        for j, (w, cc) in enumerate((("global", 0), ("zoom", 1))):
            lr, ar, _ = wins[w]
            face, ss = faces[w]
            ax = fig.add_axes([0.07 + j * 0.48, y0, 0.43 if j == 0 else 0.40, RH / H])
            ax.imshow(V2.shade(face, rgb, ss), extent=(lr[0], lr[1], ar[0], ar[1]), interpolation="lanczos",
                      aspect="equal" if cc == 0 else 1.0 / coslat)
            V2.latlon_track_templates(ax, gt, late, glon, glat, ok_seg, impact, cc, zl, za, annotate=False)
            V2.map_axes_setup(ax, cc, lr, ar, xlabel=(i == nm_ - 1), ylabel=(j == 0))
            ax.tick_params(labelsize=5.5)
            if j == 0:
                ax.set_title("%s; active seen area %.1f %%" % (m["label"], 100 * st["active_seen_area_frac_recomputed"]),
                             fontsize=6.5, loc="left", pad=2)
    cax = fig.add_axes([0.2, 0.68 / H, 0.6, 0.10 / H])
    cb = fig.colorbar(ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.ax.tick_params(labelsize=6)
    cb.set_label("H2O activity A (facets above 1 % of the srcCG median; same scale in every panel)", fontsize=6.5)
    fig.legend(handles=[Patch(facecolor=V2.ZERO_SEEN * 0.92, edgecolor="none", label="seen, A below 1 % of the srcCG median"),
                        Patch(facecolor=V2.NOT_SEEN * 0.85, edgecolor="none", label="not seen by a fitted observation")],
               loc="lower center", bbox_to_anchor=(0.5, 0.02 / H), ncol=2, fontsize=6.0, frameon=False)
    pre = out / "12M_options_activity_best"
    written = FIG.save(fig, pre)
    plt.close(fig)
    cap = ("Fitted H2O activity A of the best stage-1 setting (lowest data cost) of each non-negative activity option and of "
           "the bounded softplus c9 (rows, top to bottom: %s), on the seen surface: longitude/latitude maps of the SHAP8 12M "
           "shape model (outermost surface where overhangs put several layers at one longitude/latitude), global (left) and "
           "the impact area of the other figures (right; lon %.1f-%.1f deg, lat %.1f-%.1f deg, box in the global maps). "
           "Colour (%s, logarithmic, continuous): A on the seen facets with A above the active threshold %.3g (1 %% of the "
           "median of the srcCG activity map A_orig over the seen facets, the definition of the summaries' active seen area), "
           "one scale for all panels up to %.3g (the largest 99.9th percentile of the seen A above the threshold; larger "
           "values drawn in the top colour); light green: seen facets below the threshold (incl. A = 0); grey: facets that no "
           "fitted observation sees (their activity does not enter the fit); relief shading as the activity v2 maps. A per "
           "facet: each fit's fit/final/nonneg_facets.npz A (the full activity of the option incl. the slope-bin multipliers "
           "and D2's fine-tuning templates held at 1, as written by the fit; definitions in the file's meta_json); "
           "c9: the state's activity (nonneg_maps.state_activity of c9/ft/final, as the maps and contribution data). Active seen area per "
           "row (from these arrays; the summaries' values in brackets): %s. Black line: ground track, thick in the "
           "late-visibility window; star: Rosetta impact." % (
               "; ".join(m["label"] for m in models), zl[0], zl[1], za[0], za[1], V2.CMAP_NAME, thr, top,
               "; ".join("%s %.1f %% (%s)" % (s["label"].split(" (")[0], 100 * s["active_seen_area_frac_recomputed"],
                                              ("%.1f %%" % (100 * float(s["active_seen_area_frac_summary"])))
                                              if s["active_seen_area_frac_summary"] not in (None, "", "-") else "-")
                         for s in stats)))
    written.append(FIG.write_caption(pre, cap))
    json.dump({"threshold": thr, "srcCG_median_seen": med, "top": top, "stats": stats, "files": written,
               "command": sys.argv}, open(out / "activity_maps_manifest.json", "x"), indent=1)
    say("done: %s" % written)
    return 0


if __name__ == "__main__":
    sys.exit(main())
