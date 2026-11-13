#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/source_area_compare.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/source_area_compare.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 04c938a1f10037ee8ee5cfdf47c641d0e487b377df26766c77a25dd5f36b06d1
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: analysis. Two contribution folders: P% source-area unions (facets, area, area-weighted location, distance
#   bands from the last gauge position, Jaccard overlap, centre shift) and slope-bin pressure shares (2026-11-13:
#   bounded vs unbounded softplus).
# RUN: jobs/source_area_compare.sbatch (CPU)
# ==== end of header ====
"""Where the P% source areas and the slope-bin pressure shares of two fitted states lie, and how they differ (2026-11-13,
lead's report items for the softplus refit with slope multipliers >= 1: "where the P50 source areas lie (lon/lat, facet
counts) and whether they moved vs the unbounded fit; cliff/slope-bin pressure shares; overhang (90-180) share near the
surface").

Inputs: two contribution folders (model_ext nonneg_contrib.py: source_areas_{NG,RG}.npz, slope_bins.npz,
contrib_manifest.json), the fit base (elements.bin: facet centroids and areas; COPS_model.dat: altitudes), the frozen
nonneg_maps.py (elements reader). Read only. Writes <out>.txt and <out>.json (refuses existing files).

Definitions (all area-weighted with the facet areas):
  union U_P,g = facets in the P% source area of at least one of the 95 post-gap observations (the files' union arrays);
  location: area-weighted mean of the unit centroid vectors -> lon/lat; lon/lat 5-95 % ranges of the union facets
  (lon unwrapped around the mean); overlap of the two fits' unions: area of U_a & U_b over the area of U_a | U_b (Jaccard)
  and the fraction of each union's area inside the other; centre distance on the surface between the two means [m] (chord);
  per observation: facet count k and area of S_P (median, range);
  slope-bin shares of the modelled surface pressure per base row (slope_bins.npz fraction_NG / fraction_RG): mean, min, max
  over the observation rows before the late window, the 95 late observations, and the late observations below 20 m and
  at the last one (5.5 m).

Run (CPU job; ~4 GB of memory for the 12M mesh):
  .venv/bin/python -W ignore source_area_compare.py --a <contrib dir> --a-label "..." --b <contrib dir> --b-label "..."
      --base <fit base> --maps-code <frozen nonneg_maps.py> --out <prefix>
"""

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np

APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
sys.path.insert(0, str(APP))
from cops_cuda import figures as FIG  # noqa: E402


def load_module(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def wpct(v, w, q):
    o = np.argsort(v)
    c = np.cumsum(w[o]) / w.sum()
    return [float(v[o][min(np.searchsorted(c, qq), len(v) - 1)]) for qq in q]


BANDS_M = (0.0, 20.0, 100.0, 500.0, 2000.0, np.inf)


def location(mask, cent, area, ref):
    u = cent[mask] / np.linalg.norm(cent[mask], axis=1)[:, None]
    w = area[mask]
    m = (u * w[:, None]).sum(axis=0)
    m /= np.linalg.norm(m)
    lon0, lat0 = np.degrees(np.arctan2(m[1], m[0])), np.degrees(np.arcsin(m[2]))
    lon, lat = FIG.lonlat(cent[mask])
    dl = (lon - lon0 + 180.0) % 360.0 - 180.0
    r = np.linalg.norm(cent[mask], axis=1)
    rm = float((r * w).sum() / w.sum())
    dist = np.linalg.norm(cent[mask] - ref[None, :], axis=1)
    bands = {"%g-%g m" % (a, b): float(w[(dist >= a) & (dist < b)].sum() / w.sum()) for a, b in zip(BANDS_M[:-1], BANDS_M[1:])}
    return {"facets": int(mask.sum()), "area_m2": float(w.sum()), "lon_mean_deg": float(lon0), "lat_mean_deg": float(lat0),
            "area_fraction_by_distance_from_ref": bands, "distance_median_m": wpct(dist, w, (0.5,))[0],
            "lon_5_95_deg": [float(lon0 + v) for v in wpct(dl, w, (0.05, 0.95))],
            "lat_5_95_deg": wpct(lat, w, (0.05, 0.95)), "radius_mean_m": rm, "unit_mean": m.tolist()}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    for a in ("--a", "--a-label", "--b", "--b-label", "--base", "--maps-code", "--out"):
        ap.add_argument(a, required=True)
    args = ap.parse_args()
    out = Path(args.out)
    for suf in (".txt", ".json"):
        if Path(str(out) + suf).exists():
            raise SystemExit("%s%s exists" % (out, suf))
    nm = load_module(args.maps_code, "nonneg_maps_frozen_cmp")
    el = nm.read_elements(Path(args.base) / "elements.bin")
    cent, area = el["centroid"], el["area"]
    del el
    tab = FIG.read_table(args.base)
    row = {int(p): i for i, p in enumerate(tab["iPoint"])}
    # reference point: the gauge position of the last observation (smallest altitude of the observation rows)
    obs = (tab["refNG"] > 0) | (tab["refRG"] > 0)
    iref = np.nonzero(obs)[0][np.argmin(tab["altitude"][obs])]
    ref = np.array([tab["x"][iref], tab["y"][iref], tab["z"][iref]], np.float64)
    if np.linalg.norm(ref) < 100.0:                       # km in the table?
        raise SystemExit("unexpected position units in %s" % args.base)
    sets = {"a": (Path(args.a), args.a_label), "b": (Path(args.b), args.b_label)}
    res = {"inputs": {k: {"dir": str(d), "label": l} for k, (d, l) in sets.items()}, "base": args.base,
           "reference_point": {"iPoint": int(tab["iPoint"][iref]), "utc": tab["utc"][iref],
                               "altitude_m": float(tab["altitude"][iref]), "xyz_m": ref.tolist()},
           "source_areas": {}, "slope_bins": {}}
    L = ["# source_area_compare.py: %s (a) vs %s (b)" % (args.a_label, args.b_label),
         "# a: %s\n# b: %s" % (args.a, args.b),
         "# distances from the reference point = the gauge position of the last observation (iPoint %d, %s UTC, %.2f m)"
         % (tab["iPoint"][iref], tab["utc"][iref], tab["altitude"][iref]), ""]
    for g in ("NG", "RG"):
        z = {k: np.load(d / ("source_areas_%s.npz" % g), allow_pickle=False) for k, (d, _) in sets.items()}
        pl = [int(p) for p in z["a"]["p_levels"]]
        if pl != [int(p) for p in z["b"]["p_levels"]]:
            raise SystemExit("different P levels")
        res["source_areas"][g] = {}
        for i, p in enumerate(pl):
            U = {k: np.asarray(z[k]["union"][i]).astype(bool) for k in z}
            loc = {k: location(U[k], cent, area, ref) for k in U}
            inter = area[U["a"] & U["b"]].sum()
            union = area[U["a"] | U["b"]].sum()
            ma, mb = np.array(loc["a"]["unit_mean"]), np.array(loc["b"]["unit_mean"])
            rr = 0.5 * (loc["a"]["radius_mean_m"] + loc["b"]["radius_mean_m"])
            ent = {"a": loc["a"], "b": loc["b"], "jaccard_area": float(inter / union),
                   "a_area_in_b": float(inter / area[U["a"]].sum()), "b_area_in_a": float(inter / area[U["b"]].sum()),
                   "facets_common": int((U["a"] & U["b"]).sum()),
                   "centre_shift_m": float(rr * np.linalg.norm(ma - mb)),
                   "per_obs": {k: {"k_median": float(np.median(z[k]["k"][:, i])), "k_min": int(z[k]["k"][:, i].min()),
                                   "k_max": int(z[k]["k"][:, i].max()),
                                   "area_median_m2": float(np.median(z[k]["area_m2"][:, i])),
                                   "area_min_m2": float(z[k]["area_m2"][:, i].min()),
                                   "area_max_m2": float(z[k]["area_m2"][:, i].max())} for k in z}}
            res["source_areas"][g][str(p)] = ent
            L.append("%s P%d union: a %d facets %.4g m^2 at lon %.2f lat %.2f (5-95 %%: lon %.2f..%.2f, lat %.2f..%.2f); "
                     "b %d facets %.4g m^2 at lon %.2f lat %.2f (lon %.2f..%.2f, lat %.2f..%.2f); common facets %d, "
                     "Jaccard (area) %.3f, a in b %.3f, b in a %.3f, centre shift %.1f m; per observation S_%d: a k median "
                     "%g (%d-%d), area median %.4g m^2; b k median %g (%d-%d), area median %.4g m^2\n"
                     "    area by distance from the reference point: a %s (median %.0f m); b %s (median %.0f m)" % (
                         g, p, loc["a"]["facets"], loc["a"]["area_m2"], loc["a"]["lon_mean_deg"], loc["a"]["lat_mean_deg"],
                         *loc["a"]["lon_5_95_deg"], *loc["a"]["lat_5_95_deg"], loc["b"]["facets"], loc["b"]["area_m2"],
                         loc["b"]["lon_mean_deg"], loc["b"]["lat_mean_deg"], *loc["b"]["lon_5_95_deg"],
                         *loc["b"]["lat_5_95_deg"], ent["facets_common"], ent["jaccard_area"], ent["a_area_in_b"],
                         ent["b_area_in_a"], ent["centre_shift_m"], p,
                         ent["per_obs"]["a"]["k_median"], ent["per_obs"]["a"]["k_min"], ent["per_obs"]["a"]["k_max"],
                         ent["per_obs"]["a"]["area_median_m2"], ent["per_obs"]["b"]["k_median"],
                         ent["per_obs"]["b"]["k_min"], ent["per_obs"]["b"]["k_max"], ent["per_obs"]["b"]["area_median_m2"],
                         ", ".join("%s %.3f" % kv for kv in loc["a"]["area_fraction_by_distance_from_ref"].items()),
                         loc["a"]["distance_median_m"],
                         ", ".join("%s %.3f" % kv for kv in loc["b"]["area_fraction_by_distance_from_ref"].items()),
                         loc["b"]["distance_median_m"]))
        L.append("")
    window = FIG.late_window(tab)
    for k, (d, lab) in sets.items():
        z = np.load(d / "slope_bins.npz", allow_pickle=False)
        names = [str(s) for s in z["bin_names"]]
        kind = np.asarray(z["kind"]).astype(str)
        ip = np.asarray(z["iPoint"])
        alt = np.array([tab["altitude"][row[int(p)]] for p in ip])
        late = kind == "late_obs"
        early = kind == "obs"
        sel = {"before the late window (obs rows)": early, "late window (95 obs)": late,
               "late, altitude < 20 m": late & (alt < 20.0),
               "last observation (%.1f m)" % alt[late].min(): late & (alt == alt[late].min())}
        ent = {"multiplier": dict(zip(names, map(float, z["multiplier"]))), "shares": {}}
        L.append("slope-bin shares of the modelled surface pressure, %s (%s); multipliers %s" % (
            lab, k, ", ".join("%s %.4g" % (n, m) for n, m in ent["multiplier"].items())))
        for g in ("NG", "RG"):
            fr = np.asarray(z["fraction_%s" % g])
            use = np.asarray(z["use_%s" % g.lower()]).astype(bool)
            ent["shares"][g] = {}
            for sl, m in sel.items():
                m = m & use if "late" not in sl and "last" not in sl else m
                if not m.any():
                    continue
                v = fr[m]
                ent["shares"][g][sl] = {"rows": int(m.sum()), "mean": dict(zip(names, v.mean(axis=0).tolist())),
                                        "min": dict(zip(names, v.min(axis=0).tolist())),
                                        "max": dict(zip(names, v.max(axis=0).tolist())),
                                        "ge60_mean": float(v[:, 3:].sum(axis=1).mean()),
                                        "ge60_range": [float(v[:, 3:].sum(axis=1).min()), float(v[:, 3:].sum(axis=1).max())]}
                L.append("  %s %-34s rows %3d: mean %s; >= 60 deg mean %.3f (%.3f-%.3f)" % (
                    g, sl, m.sum(), ", ".join("%s %.3f" % (n.replace("slope_", ""), x) for n, x in zip(names, v.mean(axis=0))),
                    v[:, 3:].sum(axis=1).mean(), v[:, 3:].sum(axis=1).min(), v[:, 3:].sum(axis=1).max()))
        res["slope_bins"][k] = ent
        L.append("")
    Path(str(out) + ".txt").write_text("\n".join(L) + "\n")
    json.dump(res, open(str(out) + ".json", "x"), indent=1)
    print("\n".join(L))
    return 0


if __name__ == "__main__":
    sys.exit(main())
