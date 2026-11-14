#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/pan_candidates.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/pan_candidates.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 5e5daf0df624400ef25491eb6d8ba3ec782e7cd3d30894f7e1aa0914f4d1b7b3
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Candidate final camera poses of the panorama movie: grid around the touchdown point, framing, visible shares
#   of the big-lobe and relief NG sources (depth map), clearance, quick-looks.
# RUN: jobs/pan_candidates.sbatch (gpu002, 1 GPU)
# ==== end of header ====
"""Candidate final camera poses of the panorama movie (2026-11-14, user: "in the last frame the camera shows the fields of
view of both gauges as well as a panoramic view on the mountains and the other comet's lobe that contributes most to the
NG data"). GPU job (cops_movies Renderer).

Grid of camera positions around the touchdown point p (local frame of ng_source_regions.py: up, east, north):
c = p + Dh (cos az north + sin az east) + Hc up. Per position (--framing sum, version 1): forward = the normalized sum of
the unit directions to p, to the share-weighted mean of the body's NG sources (all 934 NG observations; facets carrying
90 % of the body's share) and (weight 0.5) to the near-east relief of the late NG signal (head facets beyond 30 m carrying
90 % of the head's late share); image up = local up projected (level horizon); the horizontal FOV that contains the
touchdown point, COPS, the body sources' mean direction and the relief (15 % margin), limited to [90, 120] deg.
--framing targets (version 2, 2026-11-14): yaw, pitch (image up = local up projected) and the focal length such that the
touchdown point lands at image row --row-td (fraction of H from the top), the body sources' mean direction at --row-body,
and their mean column at the image centre (least squares); the FOV limited to [--fov-min, --fov-max] (then yaw / pitch
re-fitted). Also the visible fraction of the body's late share (95 late observations). The camera view (face buffer,
cops_movies fixed_geometry, no supersampling) gives a depth map; a source facet counts as visible when its centroid
projects into the image, its facet faces the camera and its distance is within 0.2 % + 2 m of the depth map. Metrics:
visible fraction of the body's all-observation NG share, of the head relief's late share, the touchdown point visible
(exact ray), the clearance (nearest surface distance), the pixels per metre at the touchdown point.
Writes <out>.json, <out>.tsv and per candidate pose files <out>_pose_<id>.json (pan_to json of cops_movies.py plan),
and quick-look images of the best candidates (grey relief, sources tinted, touchdown marked)."""

import argparse
import importlib.util
import json
import math
import sys
import time
from pathlib import Path

import numpy as np


def frame_targets(a_td, b, up, east, north, W, H, opt):
    """yaw, pitch (level horizon: image up = local up projected) and focal length F that put the direction a_td at row
    opt.row_td * H, b at row opt.row_body * H and their mean column at W / 2; F limited to the FOV range."""
    from scipy.optimize import least_squares

    def fwd(yaw, pitch):
        return (math.cos(pitch) * (math.cos(yaw) * north + math.sin(yaw) * east) + math.sin(pitch) * up)

    def pix(f, F, v):
        u = up - np.dot(up, f) * f
        u /= np.linalg.norm(u)
        r = np.cross(f, u)
        z = v @ f
        return W / 2.0 + F * (v @ r) / z, H / 2.0 - F * (v @ u) / z

    Fmin = (W / 2.0) / math.tan(math.radians(opt.fov_max) / 2.0)
    Fmax = (W / 2.0) / math.tan(math.radians(opt.fov_min) / 2.0)
    m = a_td + b
    m /= np.linalg.norm(m)
    y0 = math.atan2(m @ east, m @ north)
    p0 = float(np.clip(math.asin(np.clip(m @ up, -1, 1)), -1.45, 1.45))

    def res(x, F=None):
        f = fwd(x[0], x[1])
        FF = x[2] if F is None else F
        if (a_td @ f) < 0.05 or (b @ f) < 0.05:
            return np.array([1e4, 1e4, 1e4])
        xa, ya = pix(f, FF, a_td)
        xb, yb = pix(f, FF, b)
        return np.array([(ya - opt.row_td * H) / H, (yb - opt.row_body * H) / H, (0.5 * (xa + xb) - W / 2.0) / W])
    r1 = least_squares(res, [y0, p0, 0.5 * (Fmin + Fmax)], bounds=([-10, -1.5, Fmin], [10, 1.5, Fmax]))
    return fwd(r1.x[0], r1.x[1]), float(r1.x[2])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--movies-py", required=True)
    for k in ("--work", "--geometry", "--src", "--regions", "--out"):
        ap.add_argument(k, required=True)
    ap.add_argument("--device", type=int, default=0)
    ap.add_argument("--az", default="0:360:15")
    ap.add_argument("--dh", default="30,60,120")
    ap.add_argument("--hc", default="30,60,100,150,220")
    ap.add_argument("--quicklook", type=int, default=30)
    ap.add_argument("--framing", default="sum", choices=("sum", "targets"))
    ap.add_argument("--row-td", type=float, default=0.75)
    ap.add_argument("--row-body", type=float, default=0.28)
    ap.add_argument("--fov-min", type=float, default=90.0)
    ap.add_argument("--fov-max", type=float, default=120.0)
    a = ap.parse_args()
    t0 = time.time()
    sys.path.insert(0, str(Path(a.movies_py).parent))
    spec = importlib.util.spec_from_file_location("cops_movies_pan", a.movies_py)
    CM = importlib.util.module_from_spec(spec)
    sys.argv = [a.movies_py]
    spec.loader.exec_module(CM)
    import cupy as cp
    from PIL import Image, ImageDraw
    work = Path(a.work)
    P = CM.load_prep(work)
    rnd = CM.Renderer(a.geometry, P, a.device)
    S = np.load(a.src)
    allsh, latesh = S["all_share"] / 934.0, S["late_share"] / 95.0
    lobe = np.load(a.regions + "_lobe.npy")
    Fr = np.load(a.regions + "_frame.npz")
    up, east, north, p = Fr["up"], Fr["east"], Fr["north"], Fr["ptd"]
    cent = cp.asnumpy(rnd.cent)
    norm = cp.asnumpy(rnd.norm)
    d_p = np.linalg.norm(cent - p, axis=1)
    td_lobe = int(lobe[np.argmin(d_p)])

    def top(mask, w, frac=0.9):
        i = np.flatnonzero(mask & (w > 0))
        o = i[np.argsort(-w[i])]
        c = np.cumsum(w[o])
        return o[:int(np.searchsorted(c, frac * c[-1])) + 1]
    B = top(lobe != td_lobe, allsh)
    BL = top(lobe != td_lobe, latesh)
    wBL = latesh[BL]
    Cl = top((lobe == td_lobe) & (d_p >= 30.0), latesh)
    wB, wC = allsh[B], latesh[Cl]
    xg = np.asarray(P["x"][-1], float)
    say = CM.say
    say("body sources: %d facets (90 %% of the body's all-observation share %.4f); relief: %d facets" % (
        len(B), allsh[lobe != td_lobe].sum(), len(Cl)))
    az = np.arange(*[float(v) for v in a.az.split(":")])
    dhs = [float(v) for v in a.dh.split(",")]
    hcs = [float(v) for v in a.hc.split(",")]
    cands = []
    for A in az:
        for dh in dhs:
            for hc in hcs:
                c = p + dh * (math.cos(math.radians(A)) * north + math.sin(math.radians(A)) * east) + hc * up
                cands.append({"az_deg": float(A), "dh_m": dh, "h_m": hc, "c": c})
    nodes_h = cp.asnumpy(rnd.nodes)
    faces_h = cp.asnumpy(rnd.faces).reshape(-1, 3)
    import cops_cuda.figures as FIG
    _, clear, _ = FIG.ground_track(nodes_h, faces_h, np.array([q["c"] for q in cands]))
    vis_p, _, _, _ = rnd.points_visible(np.array([q["c"] for q in cands]), np.tile(p, (len(cands), 1)),
                                        np.full(len(cands), -1))
    W, H = CM.W, CM.H
    for i, q in enumerate(cands):
        c = q["c"]
        q["clearance_m"] = float(clear[i])
        q["touchdown_ray_visible"] = bool(vis_p[i])
        dB = cent[B] - c
        dB /= np.linalg.norm(dB, axis=1)[:, None]
        mB = CM.unit((dB * wB[:, None]).sum(0))
        dC = cent[Cl] - c
        dC /= np.linalg.norm(dC, axis=1)[:, None]
        mC = CM.unit((dC * wC[:, None]).sum(0))
        f = CM.unit(CM.unit(p - c) + mB + 0.5 * mC)
        if a.framing == "targets":
            f, Ft = frame_targets(CM.unit(p - c), mB, up, east, north, W, H, a)
            q["F_targets"] = Ft
        u = CM.unit(up - np.dot(up, f) * f)
        r = np.cross(f, u)
        need = 0.0
        for v in (CM.unit(p - c), CM.unit(xg - c), mB, mC):
            z = v @ f
            if z <= 0.05:
                need = 1e9
                break
            need = max(need, abs(v @ r) / z, abs(v @ u) / z * W / H)
        hfov = 2 * math.degrees(math.atan(need * 1.15)) if need < 1e8 else 180.0
        q["hfov_required_deg"] = hfov
        q["hfov_deg"] = float(min(max(hfov, 90.0), 120.0))
        if a.framing == "targets":
            q["hfov_deg"] = 2 * math.degrees(math.atan((W / 2.0) / q["F_targets"]))
            q["hfov_required_deg"] = max(hfov, q["hfov_deg"]) if need < 1e8 else 180.0
        q["forward"], q["up"] = f, u
    ok = [q for q in cands if q["clearance_m"] >= 10.0 and q["touchdown_ray_visible"] and q["hfov_required_deg"] <= 120.0]
    say("%d candidates, %d with clearance >= 10 m, the touchdown point visible and FOV <= 120 deg (%.0f s)" % (
        len(cands), len(ok), time.time() - t0))
    for j, q in enumerate(ok):
        F = (W / 2.0) / math.tan(math.radians(q["hfov_deg"]) / 2.0)
        vw = rnd.fixed_geometry(q["c"], q["forward"], q["up"], "perspective", F, (0.0, 0.0), ss=1)
        dmap = np.full(W * H, np.inf)
        dmap[cp.asnumpy(vw["idx"])] = np.linalg.norm(cp.asnumpy(vw["pts"]) - q["c"], axis=1)
        f, u = q["forward"], q["up"]
        r = np.cross(f, u)

        def visible(ix):
            v = cent[ix] - q["c"]
            z = v @ f
            dist = np.linalg.norm(v, axis=1)
            okz = z > 0.05
            X = np.where(okz, (v @ r) / np.where(okz, z, 1), 0)
            Y = np.where(okz, (v @ u) / np.where(okz, z, 1), 0)
            px = W / 2.0 + F * X
            py = H / 2.0 - F * Y
            inside = okz & (px >= 0) & (px < W) & (py >= 0) & (py < H)
            pi = np.where(inside, np.floor(py).astype(np.int64) * W + np.floor(px).astype(np.int64), 0)
            facing = np.einsum("ij,ij->i", norm[ix], -v) > 0
            return inside & facing & (dist <= dmap[pi] * 1.002 + 2.0)
        vb = visible(B)
        vc = visible(Cl)
        q["body_late_share_visible_frac"] = float(wBL[visible(BL)].sum() / wBL.sum())
        q["body_all_share_visible_frac"] = float(wB[vb].sum() / wB.sum())
        q["relief_late_share_visible_frac"] = float(wC[vc].sum() / wC.sum())
        q["px_per_m_at_touchdown"] = float(F / np.linalg.norm(p - q["c"]))
        q["distance_to_touchdown_m"] = float(np.linalg.norm(p - q["c"]))
        q["F"] = F
        q["sky_frac"] = float(1 - vw["idx"].size / (W * H))
        if (j + 1) % 25 == 0:
            say("evaluated %d / %d (%.0f s)" % (j + 1, len(ok), time.time() - t0))
        del vw
    for q in ok:   # score: both source sets visible, close to the spacecraft
        q["score"] = (q["body_all_share_visible_frac"] * q["relief_late_share_visible_frac"]
                      * min(1.0, q["px_per_m_at_touchdown"] / 10.0) ** 0.5)
    ok.sort(key=lambda q: -q["score"])
    out = Path(a.out)
    keys = ["az_deg", "dh_m", "h_m", "clearance_m", "distance_to_touchdown_m", "hfov_required_deg", "hfov_deg",
            "body_all_share_visible_frac", "body_late_share_visible_frac", "relief_late_share_visible_frac", "px_per_m_at_touchdown", "sky_frac", "score"]
    with open(str(out) + ".tsv", "x") as fo:
        fo.write("id\t" + "\t".join(keys) + "\n")
        for i, q in enumerate(ok):
            q["id"] = "g%03d" % i
            fo.write(q["id"] + "\t" + "\t".join("%.4g" % q[k] for k in keys) + "\n")
    sun = P["sun"][-1]
    cosz = cp.asnumpy(rnd.norm @ cp.asarray(CM.unit(sun)))
    for i, q in enumerate(ok[:a.quicklook]):
        F = q["F"]
        vw = rnd.fixed_geometry(q["c"], q["forward"], q["up"], "perspective", F, (0.0, 0.0), ss=1)
        img = np.zeros((H * W, 3))
        idx = cp.asnumpy(vw["idx"])
        fid = cp.asnumpy(vw["fid"])
        br = 0.15 + 0.85 * np.clip(cosz[fid], 0, 1)
        img[idx] = br[:, None] * 0.8
        tB = np.zeros(len(cent), bool)
        tB[B] = True
        tC = np.zeros(len(cent), bool)
        tC[Cl] = True
        img[idx[tB[fid]]] = np.array([1.0, 0.45, 0.1]) * (0.4 + 0.6 * br[tB[fid]][:, None])
        img[idx[tC[fid]]] = np.array([0.2, 0.9, 1.0]) * (0.4 + 0.6 * br[tC[fid]][:, None])
        im = Image.fromarray((np.clip(img, 0, 1) * 255).astype(np.uint8).reshape(H, W, 3))
        dr = ImageDraw.Draw(im)
        f, u = q["forward"], q["up"]
        r = np.cross(f, u)
        for pt, col in ((p, (255, 255, 0)), (xg, (255, 0, 255))):
            v = pt - q["c"]
            z = v @ f
            if z > 0:
                x_, y_ = W / 2 + F * (v @ r) / z, H / 2 - F * (v @ u) / z
                dr.ellipse((x_ - 10, y_ - 10, x_ + 10, y_ + 10), outline=col, width=3)
        dr.rectangle((0, 0, W, 40), fill=(0, 0, 0))
        dr.text((10, 8), "%s az %.0f dh %.0f m h %.0f m fov %.0f: body vis %.2f relief vis %.2f %.1f px/m d %.0f m" % (
            q["id"], q["az_deg"], q["dh_m"], q["h_m"], q["hfov_deg"], q["body_all_share_visible_frac"],
            q["relief_late_share_visible_frac"], q["px_per_m_at_touchdown"], q["distance_to_touchdown_m"]),
            font=CM.font(24), fill=(255, 255, 255))
        im.resize((W // 2, H // 2)).save(str(out) + "_ql_%s.png" % q["id"])
        del vw
    for q in ok:
        q["camera_position_m"] = q.pop("c").tolist()
        q["forward"] = q["forward"].tolist()
        q["up"] = q["up"].tolist()
    json.dump({"definition": __doc__, "inputs": vars(a), "pivot_m": p.tolist(), "local_up": up.tolist(),
               "body_source_facets": int(len(B)), "relief_facets": int(len(Cl)), "candidates_total": len(cands),
               "candidates_ok": ok, "seconds": time.time() - t0}, open(str(out) + ".json", "x"), indent=1)
    say("done in %.0f s; best: %s" % (time.time() - t0, json.dumps({k: ok[0][k] for k in keys + ["id"]})))
    return 0


if __name__ == "__main__":
    sys.exit(main())
