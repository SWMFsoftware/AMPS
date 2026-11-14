#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/ng_source_regions.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/ng_source_regions.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2e330881e41c86b9d123fafcb75bd844f8a97575bfc68117ff1ed11b15649a36
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Lobes (neck plane), touchdown lobe, NG shares by lobe / distance / azimuth, dominant regions, skyline and
#   relief around the touchdown point, map figure (.png/.eps 300 dpi).
# RUN: srun on normal128cores (CPU), after ng_source_map.py
# ==== end of header ====
"""Regions of the modelled NG surface signal (from ng_source_map.py) for the panorama camera path (lead / user 2026-11-13
~20:05 CST): the two lobes, the lobe of the touchdown point, shares by lobe / distance / azimuth, the dominant regions, the
relief (skyline, highest terrain, cliffs) around the touchdown point, and a map figure (.png + .eps, 300 dpi).

Lobes: the plane of the neck = the plane of minimum surface area per unit length (facet area in a 20 m slab / 20 m) between
the two area-weighted 2-means centres of the facet centroids, its normal searched within 40 deg of the centre-centre axis
(subsample of every 10th facet for the search; the classification uses all facets). The lobe with the larger area is called
the big lobe (body), the other the small lobe (head). Local frame at the touchdown point: up = area-weighted mean facet
normal within 30 m, east = Z x up (Z = the body-frame spin axis), north = up x east; azimuth from north through east.
Shares: late = (sum over the 95 late observations of the per-observation shares) / 95 (mean share of an observation's NG
surface signal); all = the same over the 95 + 839 NG observations of the descent; Pa variants: share of the summed Pa.
Usage: ng_source_regions.py --geometry G --src <ng_source_map>.npz --touchdown <work>/touchdown.npz --work <work> --out P
"""

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

DEG = np.pi / 180


def kmeans2(c, w, it=40):
    m = (c * w[:, None]).sum(0) / w.sum()
    cov = np.cov((c - m).T, aweights=w)
    e = np.linalg.eigh(cov)[1][:, -1]
    s = (c - m) @ e
    lab = s > 0
    for _ in range(it):
        mu = [(c[lab == j] * w[lab == j, None]).sum(0) / w[lab == j].sum() for j in (0, 1)]
        new = np.linalg.norm(c - mu[1], axis=1) < np.linalg.norm(c - mu[0], axis=1)
        if np.array_equal(new, lab):
            break
        lab = new
    return np.array(mu), lab


def neck_plane(c, a, mu):
    u = mu[1] - mu[0]
    u /= np.linalg.norm(u)
    t1 = np.cross(u, [0.0, 0.0, 1.0]) if abs(u[2]) < 0.9 else np.cross(u, [1.0, 0.0, 0.0])
    t1 /= np.linalg.norm(t1)
    t2 = np.cross(u, t1)
    best = None
    for th in np.arange(0, 41, 2.5):
        for ph in (np.arange(0, 360, 15) if th > 0 else [0.0]):
            n = np.cos(th * DEG) * u + np.sin(th * DEG) * (np.cos(ph * DEG) * t1 + np.sin(ph * DEG) * t2)
            s = c @ n
            s0, s1 = mu[0] @ n, mu[1] @ n
            lo, hi = s0 + 0.15 * (s1 - s0), s0 + 0.85 * (s1 - s0)
            edges = np.arange(lo - 10, hi + 10 + 5, 5.0)
            h, _ = np.histogram(s, edges, weights=a)
            sl = np.convolve(h, np.ones(4), "valid") / 20.0          # 20 m slab, per metre
            ctr = edges[:-3][:len(sl)] + 10.0
            ok = (ctr >= lo) & (ctr <= hi)
            j = np.argmin(np.where(ok, sl, np.inf))
            if best is None or sl[j] < best[0]:
                best = (float(sl[j]), n.copy(), float(ctr[j]), float(th), float(ph))
    return best


def lonlat(p):
    p = np.atleast_2d(p)
    r = np.linalg.norm(p, axis=1)
    return np.degrees(np.arctan2(p[:, 1], p[:, 0])), np.degrees(np.arcsin(p[:, 2] / r))


def main():
    ap = argparse.ArgumentParser()
    for k in ("--geometry", "--src", "--touchdown", "--work", "--out"):
        ap.add_argument(k, required=True)
    a = ap.parse_args()
    t0 = time.time()
    out = Path(a.out)
    if Path(str(out) + ".json").exists():
        raise SystemExit("%s.json exists" % out)
    g = Path(a.geometry)
    fg = np.memmap(g / "mesh_facegeom.f8", dtype="<f8", mode="r").reshape(-1, 8)
    cent = np.ascontiguousarray(fg[:, 0:3])
    norm = np.ascontiguousarray(fg[:, 3:6])
    area = np.ascontiguousarray(fg[:, 6])
    del fg
    nf = len(area)
    S = np.load(a.src)
    sh = {"late": S["late_share"] / 95.0, "all": S["all_share"] / 934.0,
          "late_Pa": S["late_Pa"] / S["late_Pa"].sum(), "all_Pa": S["all_Pa"] / S["all_Pa"].sum()}
    for k, v in sh.items():
        if len(v) != nf:
            raise SystemExit("facet count differs")
    TD = np.load(a.touchdown)
    ptd = np.asarray(TD["point"], float)
    W = Path(a.work)
    P = np.load(W / "prep.npz")
    PL = np.load(W / "plan.npz")
    R = {"definition": __doc__, "inputs": vars(a), "touchdown_point_m": ptd.tolist(),
         "touchdown_lonlat_deg": [float(v[0]) for v in lonlat(ptd)]}
    # lobes
    sub = slice(None, None, 10)
    mu, lab_km = kmeans2(cent[sub], area[sub])
    best = neck_plane(cent[sub], area[sub], mu)
    _, n_neck, s_neck, th, ph = best
    side = cent @ n_neck > s_neck                       # True on the mu[1] side
    A1 = area[side].sum()
    A0 = area[~side].sum()
    big_is_1 = A1 > A0
    lobe = np.where(side == big_is_1, 1, 0).astype(np.uint8)  # 1 = big lobe (body), 0 = small lobe (head)
    np.save(str(out) + "_lobe.npy", lobe)
    lab_full_km = np.linalg.norm(cent - mu[1], axis=1) < np.linalg.norm(cent - mu[0], axis=1)
    agree = float(area[(lab_full_km == side)].sum() / area.sum())
    td_side = bool(ptd @ n_neck > s_neck)
    td_lobe = 1 if td_side == big_is_1 else 0
    names = {0: "small lobe (head)", 1: "big lobe (body)"}
    R["lobes"] = {"neck_plane_normal": n_neck.tolist(), "neck_plane_offset_m": s_neck,
                  "neck_min_area_per_length_m": best[0], "neck_normal_tilt_from_centre_axis_deg": th,
                  "kmeans_centres_m": mu.tolist(), "area_km2": {names[0]: float(area[lobe == 0].sum() / 1e6),
                                                                 names[1]: float(area[lobe == 1].sum() / 1e6)},
                  "plane_vs_2means_area_agreement": agree,
                  "touchdown_lobe": names[td_lobe],
                  "touchdown_distance_from_neck_plane_m": float(abs(ptd @ n_neck - s_neck))}
    other = 1 - td_lobe
    # local frame
    d = np.linalg.norm(cent - ptd, axis=1)
    near = d < 30.0
    up = (norm[near] * area[near, None]).sum(0)
    up /= np.linalg.norm(up)
    east = np.cross([0.0, 0.0, 1.0], up)
    east /= np.linalg.norm(east)
    north = np.cross(up, east)
    rad = ptd / np.linalg.norm(ptd)
    R["local_frame"] = {"up": up.tolist(), "east": east.tolist(), "north": north.tolist(),
                        "up_vs_radial_deg": float(np.degrees(np.arccos(np.clip(up @ rad, -1, 1)))),
                        "up_vs_touchdown_facet_normal_deg": float(np.degrees(np.arccos(np.clip(up @ TD["normal"], -1, 1))))}
    v = cent - ptd
    h = v @ up
    xe, xn = v @ east, v @ north
    dh = np.hypot(xe, xn)
    az = (np.degrees(np.arctan2(xe, xn)) + 360) % 360

    def dir_of(p, p0=ptd):
        w = p - p0
        dd = np.linalg.norm(w)
        return {"distance_m": float(dd), "azimuth_deg": float((np.degrees(np.arctan2(w @ east, w @ north)) + 360) % 360),
                "elevation_deg": float(np.degrees(np.arcsin((w @ up) / dd))), "height_above_touchdown_m": float(w @ up),
                "horizontal_distance_m": float(np.hypot(w @ east, w @ north))}
    # shares
    dist_edges = [0, 10, 30, 100, 300, 1000, 3000, 1e9]
    R["shares"] = {}
    for k, s in sh.items():
        e = {"by_lobe": {names[j]: float(s[lobe == j].sum()) for j in (0, 1)}}
        e["by_distance_from_touchdown_m"] = {"%g-%g" % (dist_edges[i], dist_edges[i + 1]): float(
            s[(d >= dist_edges[i]) & (d < dist_edges[i + 1])].sum()) for i in range(len(dist_edges) - 1)}
        e["by_azimuth_30deg_beyond_30m"] = {"%d-%d" % (i, i + 30): float(s[(az >= i) & (az < i + 30) & (d >= 30)].sum())
                                           for i in range(0, 360, 30)}
        e["above_touchdown_height_beyond_30m"] = float(s[(h > 0) & (d >= 30)].sum())
        e["share_weighted_centroid"] = dir_of((cent * s[:, None]).sum(0) / s.sum())
        so = s * (lobe == other)
        e["other_lobe_share_weighted_centroid"] = dir_of((cent * so[:, None]).sum(0) / so.sum())
        lo_, la_ = lonlat((cent * so[:, None]).sum(0) / so.sum())
        e["other_lobe_share_weighted_centroid"]["lonlat_deg"] = [float(lo_[0]), float(la_[0])]
        # top 100 m voxels (all facets) and top 200 m voxels of the other lobe
        for vname, vsz, mask in (("top_voxels_100m", 100.0, np.ones(nf, bool)),
                                 ("other_lobe_top_voxels_200m", 200.0, lobe == other)):
            key = np.floor(cent[mask] / vsz).astype(np.int64)
            kk = (key[:, 0] + 100) * 1000000 + (key[:, 1] + 100) * 1000 + (key[:, 2] + 100)
            uk, inv = np.unique(kk, return_inverse=True)
            ssum = np.bincount(inv, weights=s[mask])
            order = np.argsort(-ssum)[:8]
            lst = []
            for o in order:
                m = inv == o
                ww = s[mask][m]
                pc = (cent[mask][m] * ww[:, None]).sum(0) / ww.sum()
                lo_, la_ = lonlat(pc)
                dd = dir_of(pc)
                dd.update({"share": float(ssum[o]), "centre_m": pc.tolist(), "lonlat_deg": [float(lo_[0]), float(la_[0])],
                           "lobe": names[int(round(float(np.average(lobe[mask][m], weights=ww))))]})
                lst.append(dd)
            e[vname] = lst
        R["shares"][k] = e
    # relief around the touchdown point
    p0h = 2.0
    far = dh > 3.0
    el = np.degrees(np.arctan2(h - p0h, dh))
    azb = np.floor(az / 2).astype(int)
    sky = np.full(180, -90.0)
    skyd = np.zeros(180)
    for b in range(180):
        m = far & (azb == b)
        if m.any():
            j = np.argmax(np.where(m, el, -np.inf))
            sky[b] = el[j]
            skyd[b] = dh[j]
    np.savez(str(out) + "_skyline.npz", az_deg=np.arange(180) * 2 + 1.0, elev_deg=sky, dist_m=skyd)
    R["skyline_from_touchdown_2m"] = {"max_elev_deg": float(sky.max()), "az_of_max_deg": float(np.argmax(sky) * 2 + 1),
                                      "dist_of_max_m": float(skyd[np.argmax(sky)]),
                                      "by_30deg": {"%d-%d" % (i, i + 30): {"max_elev_deg": float(sky[i // 2:(i + 30) // 2].max()),
                                                   "at_dist_m": float(skyd[i // 2 + np.argmax(sky[i // 2:(i + 30) // 2])])}
                                                   for i in range(0, 360, 30)}}
    relief = {}
    for rmax in (500.0, 1500.0):
        m = dh < rmax
        sec = {}
        for i in range(0, 360, 30):
            mm = m & (az >= i) & (az < i + 30)
            j = np.flatnonzero(mm)[np.argmax(h[mm])]
            lo_, la_ = lonlat(cent[j])
            sec["%d-%d" % (i, i + 30)] = {"max_height_m": float(h[j]), "horizontal_distance_m": float(dh[j]),
                                          "lonlat_deg": [float(lo_[0]), float(la_[0])], "lobe": names[int(lobe[j])]}
        relief["within_%dm" % rmax] = sec
    slope = np.degrees(np.arccos(np.clip(norm @ up, -1, 1)))
    cl = {}
    for i in range(0, 360, 30):
        mm = (dh < 500) & (dh > 20) & (az >= i) & (az < i + 30)
        cl["%d-%d" % (i, i + 30)] = {"area_slope_gt_60deg_m2": float(area[mm & (slope > 60)].sum()),
                                     "area_m2": float(area[mm].sum())}
    relief["cliffs_20_500m_slope_vs_local_up"] = cl
    R["relief"] = relief
    # the gauges at the last frame and the detach camera, in the local frame
    k_end = len(P["utc"]) - 1
    kd = int(PL["k_detach"])
    R["gauges_last_frame"] = {"utc": str(P["utc"][k_end]), "position": dir_of(P["x"][k_end])}
    for gname in ("b_rg", "b_ng"):
        b = P[gname][k_end]
        R["gauges_last_frame"][gname] = {"azimuth_deg": float((np.degrees(np.arctan2(b @ east, b @ north)) + 360) % 360),
                                         "elevation_deg": float(np.degrees(np.arcsin(b @ up)))}
    R["detach_camera"] = {"utc": str(PL["detach_utc"]), "position": dir_of(PL["cam_x_frames"][kd]),
                          "axis_b_rg": {"azimuth_deg": float((np.degrees(np.arctan2(PL["cam_b_rg_frames"][kd] @ east,
                                                                                    PL["cam_b_rg_frames"][kd] @ north)) + 360) % 360),
                                        "elevation_deg": float(np.degrees(np.arcsin(PL["cam_b_rg_frames"][kd] @ up)))}}
    np.savez(str(out) + "_frame.npz", up=up, east=east, north=north, ptd=ptd, n_neck=n_neck, s_neck=s_neck)
    # figure
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    vd = up - (up @ n_neck) * n_neck                     # view along the touchdown up, both lobes side by side
    vd /= np.linalg.norm(vd)
    ex = n_neck - (n_neck @ vd) * vd
    ex /= np.linalg.norm(ex)
    if (mu[1] - mu[0]) @ ex * (1 if big_is_1 else -1) < 0:
        ex = -ex                                          # big lobe to the right
    ey = np.cross(vd, ex)
    X, Y, Z = cent @ ex, cent @ ey, cent @ vd
    res = 3.0
    ix = np.floor((X - X.min()) / res).astype(np.int64)
    iy = np.floor((Y - Y.min()) / res).astype(np.int64)
    nx, ny = ix.max() + 1, iy.max() + 1
    pix = iy * nx + ix
    o = np.lexsort((-Z, pix))
    first = np.ones(len(o), bool)
    first[1:] = pix[o][1:] != pix[o][:-1]
    fr = o[first]
    fr = fr[norm[fr] @ vd > 0]
    shade = np.full(nx * ny, np.nan)
    shade[pix[fr]] = 0.25 + 0.75 * np.clip(norm[fr] @ vd, 0, 1)
    extent = [X.min() / 1e3, (X.min() + nx * res) / 1e3, Y.min() / 1e3, (Y.min() + ny * res) / 1e3]
    fig, axs = plt.subplots(1, 3, figsize=(16, 5.4), gridspec_kw={"width_ratios": [1, 1, 1.15]})
    for ax, k, ttl in ((axs[0], "late", "(a) 95 late observations (W10x5)"), (axs[1], "all", "(b) all 934 NG observations")):
        dens = np.full(nx * ny, np.nan)
        val = sh[k][fr] / area[fr] * 1e6                  # share per km^2
        dens[pix[fr]] = np.where(val > 0, val, np.nan)
        ax.imshow(shade.reshape(ny, nx), cmap="gray", vmin=0, vmax=1.2, origin="lower", extent=extent, interpolation="nearest")
        vmax = np.nanpercentile(dens, 99.9)
        im = ax.imshow(dens.reshape(ny, nx), cmap="inferno", norm=LogNorm(vmin=vmax * 1e-4, vmax=vmax), origin="lower",
                       extent=extent, interpolation="nearest", alpha=0.9)
        ax.plot([(ptd @ ex) / 1e3], [(ptd @ ey) / 1e3], marker="*", ms=14, mfc="cyan", mec="k", ls="none", label="touchdown")
        c_o = R["shares"][k]["other_lobe_share_weighted_centroid"]
        po = (cent * (sh[k] * (lobe == other))[:, None]).sum(0) / (sh[k] * (lobe == other)).sum()
        ax.plot([(po @ ex) / 1e3], [(po @ ey) / 1e3], marker="D", ms=8, mfc="lime", mec="k", ls="none",
                label="other-lobe NG centroid (%.0f m)" % c_o["distance_m"])
        ax.set_title(ttl + "\n%s share %.3f, %s %.3f" % (names[td_lobe], R["shares"][k]["by_lobe"][names[td_lobe]], names[other],
                                                       R["shares"][k]["by_lobe"][names[other]]), fontsize=10)
        ax.set_xlabel("km (along the neck normal)")
        ax.set_ylabel("km")
        ax.legend(loc="lower left", fontsize=8)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.02, label="mean share of the NG surface signal per km²")
    ax = axs[2]
    azc = np.arange(180) * 2 + 1.0
    ax.plot(azc, sky, color="0.3", lw=1.2, label="skyline elevation from 2 m above touchdown")
    ax.set_xlabel("azimuth from touchdown (deg, N through E)")
    ax.set_ylabel("skyline elevation (deg)")
    ax.set_xlim(0, 360)
    ax.set_xticks(range(0, 361, 45))
    ax2 = ax.twinx()
    for k, col in (("late", "C3"), ("all", "C0")):
        hs = np.bincount(np.floor(az[d >= 30] / 10).astype(int), weights=sh[k][d >= 30], minlength=36)
        ax2.step(np.arange(36) * 10, hs, where="post", color=col, label="NG share per 10° beyond 30 m: %s" % k)
    for gname, col in (("b_ng", "C2"), ("b_rg", "C1")):
        ax.axvline(R["gauges_last_frame"][gname]["azimuth_deg"], color=col, ls="--", lw=1,
                   label="%s boresight azimuth at contact" % gname[2:].upper())
    ax2.set_ylabel("share")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7, loc="upper right")
    ax.set_title("(c) relief and NG source direction around the touchdown point", fontsize=10)
    fig.tight_layout()
    for ext in ("png", "eps"):
        fig.savefig(str(out) + "_map." + ext, dpi=300)
    R["figure"] = str(out) + "_map.png"
    R["seconds"] = time.time() - t0
    json.dump(R, open(str(out) + ".json", "x"), indent=1)
    print(json.dumps({k: R[k] for k in ("lobes", "local_frame")}, indent=1))
    for k in sh:
        print(k, json.dumps(R["shares"][k]["by_lobe"]), json.dumps(R["shares"][k]["by_distance_from_touchdown_m"]))
    print("done in %.0f s" % R["seconds"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
