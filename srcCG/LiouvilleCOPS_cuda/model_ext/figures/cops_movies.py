#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/cops_movies.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/cops_movies.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: c2cf00fb8a00494bb8a25e8688239e0d817312f0a42b9a260d8426940d175e15
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: COPS movies (prep, plan, render, compose, encode, checks); 2026-11-14: plan --camera detach --pan-to / --pan-
#   start-frame / --pan-end-frame (camera moves after the detach to a final pose: smootherstep, arc around the pivot,
#   quaternion slerp, focal length geometric; held after the end frame, labelled) and compose --no-sub-mark; defaults
#   unchanged (tests/movies_detach300_pan_20261113, validation/regression of the panorama set).
# RUN: jobs/cops_movies_{prep,render,encode,qc}.sbatch (Slurm, gpu002 / normal128cores); see 12M/movies_*/HOW-TO-
#   RERUN.md
# ==== end of header ====
"""One-minute movies of Rosetta's final descent (2016-09-30) as seen from the COPS gauges, with the softplus activity
model after 5 fit / fine-tuning cycles (final state).

User request 2026-11-12: "create 1 min movies: 1) the nucleus as seen from RG, 2) the nucleus as seen from NG, and 3) a
movie that has 3 panels that shows what is seen from RG and NG, and the third panel is comparison of the model with COPS
RG and NG with time and vertical lines on the plot showing what corresponds to the current time; make the resolution
appropriate for presentation projectors and for placing online". User answers (2026-11-12):
  * whole descent 00:00:04 -> 10:39:26 UTC (16.7 km -> 5.5 m), frame times UNIFORM IN LOG ALTITUDE (equal movie time per
    factor of altitude); UTC and altitude on every frame;
  * ONE camera-like view for both gauges: a camera at the COPS position looking along the RG axis (spacecraft +Z),
    framed on the nucleus, perspective, a fisheye (azimuthal equidistant) only where the nucleus is too wide for a
    perspective view (said on the frame); sunlight with cast shadows; movie 1 tints the surface by the RG acceptance
    weight -l.b_RG, movie 2 by the NG acceptance weight (1 - l.b_NG in front of the NG, the shield fraction behind it),
    with a colour bar of the weight; movie 3: both views side by side + NG and RG time series (COPS data and model) with
    a vertical line at the frame's altitude;
  * H.264 MP4 1920x1080 30 fps (CRF 18) + 1280x720 (CRF 26, faststart); ffmpeg = imageio-ffmpeg installed in
    $COPS_OUTPUT_ROOT/tools/imageio_ffmpeg (not in the .venv).

GEOMETRY (as the fitted model; LiouvilleCOPS src/Geometry.cpp, src/Gauges.cpp, src/Liouville.h)
  gauge position, Sun direction (LT+S), NG / RG boresights per frame time: cops_movies_geom.cpp (Geometry.cpp's SPICE
  calls; GaugeOffsetSC 0 0 1.5 m of the fit's base), NG boresight = ROS_ROSINA_COPS_G0 (-226331) = -Y_sc, RG = _G1
  (-226332) = +Z_sc. l = unit vector from a facet (centroid) to the gauge (molecule direction).
  RG weight (RGAcceptance = hemisphere; flux into the aperture): w_RG = max(0, -l.b_RG).
  NG weight (NGAcceptance = srccg_shield, NGHeadOnWeight = on): mu = l.b_NG; mu < 0 (in front): w_NG = 1 - mu (1..2);
  mu >= 0 (behind): w_NG = SF(mu) if l.b_RG < 0 else 0, SF = the shield-fraction table of Gauges.cpp BuildShieldTable
  (reproduced bit for bit by cops_movies_geom.cpp shield).
  illumination: cos SZA = n.s (n = elements.bin / maps-geometry facet normal, s = Sun unit vector); cast shadow by a GPU
  shadow map (orthographic depth buffer along s, SHADOW_RES^2 texels over the nucleus; a facet is lit when it is the
  front-most surface at its centroid's texel or within a slope-scaled depth tolerance); validated against the exact
  BVH shadows of the maps geometry (validate). brightness = AMBIENT + DIFFUSE * max(0, cos SZA) * lit.
  camera: position = gauge, axis = b_RG, image up = b_NG projected, right = axis x up; perspective (X, Y) = (x/z, y/z)
  while the nucleus fits within PERSPECTIVE_MAX_DEG of the axis, then azimuthal equidistant fisheye (X, Y) = theta
  (cos phi, sin phi) for the rest of the movie; pixel = (W/2 + F (X - Xc), H/2 - F (Y - Yc)): the principal point
  (Xc, Yc) is the centre of the nucleus' bounding box in that frame (the camera axis stays the RG axis; a cropped wide
  sensor), F from the box with margin MARGIN, smoothed (min over +-15 frames, then mean over +-15 frames, never larger
  than the frame's own fit) (plan); z-buffer rasterization of the 12M facets (cops_cuda kernels/render.cu raster_faces,
  unchanged), front-facing facets only, SS x SS supersampling; empty pixels enclosed by the drawn nucleus are ray-cast
  exactly (render.cu ray_first_hit) and take the first facet hit.
  tint per supersampled pixel from its exact ray direction d (l = -d): w_RG = max(0, d.b_RG), NG mu = -d.b_NG; viridis
  of w / w_max (RG 0..1, NG 0..2), grey where w = 0, times the facet brightness.

Subcommands (run from <AMPS>/srcCG/LiouvilleCOPS_cuda with .venv/bin/python; jobs/cops_movies_*.sbatch):
  prep      frame times, SPICE geometry of the frames and of the base rows, shield table, geometry checks   (CPU)
  plan      camera framing per frame, sub-spacecraft intercept points                                       (GPU)
  validate  shadow map vs the maps geometry (13 times), cos / shadow / weights / visibility vs the kernel records of
            base rows, test frames                                                                         (GPU)
  render    clean view frames (RG tint, NG tint) of the frames it claims (dynamic claims: split-independent) (GPU)
  compose   movie frames 1, 2, 3 (overlays, colour bars, time series)                                       (CPU)
  compare   two frame folders pixel by pixel (split independence)                                           (CPU)
  encode    MP4 (1920x1080 CRF 18; 1280x720 CRF 26 faststart)                                               (CPU)
Added 2026-11-13 (new options; the defaults reproduce the first movies bit for bit):
  fov          IK FOVs (getfov -226331 NG, -226332 RG) and their body-frame bound vectors at every frame and base-row
               time -> fov.npz (cops_movies_geom fov / fovinfo)                                             (CPU)
  validate_fov the ik_fov masks against an analytic instrument-frame test, weight ranges inside the FOVs, kernel-record
               shares inside the FOVs, masked test frames                                                   (GPU)
  sc           spacecraft geometry for the shadow (bus box of ROS_SPACECRAFT_C1..C8, the two solar-array wings
               ROS_SA+-Y_C1..C4, attitude, SA orientation held after SA_HOLD_UTC) -> sc.npz              (CPU)
  validate_sc  shadow checks: bus-shadow centre vs the analytic line-mesh intersection, disc-sampling convergence,
               umbra / penumbra widths on a synthetic plane vs 2 alpha L                                    (GPU)
  diagnose     classify the dark (ambient-only) pixels of rendered frames by facet: cos <= 0, cast shadow (shadow map
               and exact shadow ray), false shadow; annotated stills                                        (GPU)
  frametable   frames.tsv from frames/meta (with the mask statistics of masked frames)                    (CPU)
  render options: --mask-rg / --mask-ng none|fit|hemisphere|ik_fov|cone:<deg> (observable part; tint inside, sunlit grey
  outside, boundary line), --tint-rg / --tint-ng weight|uniform, --sc-shadow (spacecraft shadow: direct sunlight x the
  unblocked fraction of the solar disc, SC_SAMPLES points), --no-boundary.
Added 2026-11-13, second set (new options; the defaults reproduce every earlier movie bit for bit):
  plan --camera fixed --camera-altitude 100 --fov-deg 90 --base --kernels --geom-exe (user: "freeze the location of the
               camera at about 100 m altitude in the frame moving with the nucleus"): the camera = the COPS position
               (ROSETTA + GaugeOffsetSC, cops_movies_geom geom as prep) and the RG / NG boresights at the time the COPS
               altitude is --camera-altitude (time interpolated in log altitude over the observation rows, as prep's frame
               times), frozen in 67P/C-G_CK for every frame; axis = that RG boresight, image up = that NG boresight
               projected; perspective with --fov-deg horizontal field of view, principal point at the centre (CPU+GPU)
  render (with a fixed-camera plan.npz; --mask-rg / --mask-ng ik_fov or cone:<deg> required): the camera view is the same
               in every frame (face buffer once per worker, its sha256 in every frame's meta); per frame, from the TRUE
               gauge position and attitude of the frame time: a surface point P is tinted for a gauge when (1) d = P -
               gauge is inside the gauge's mask (IK FOV, Geometry.cpp InsideFOV on the frame's body-frame bounds), (2) its
               facet faces the gauge (n.(gauge - P) > 0), (3) the segment P -> gauge crosses no facet (line of sight:
               Renderer.gauge_los, a perspective depth map from the gauge over the window of the selected directions +
               exact segment tests against the front-most facets of the 3 x 3 texels around P and against the facets cut
               by the map's near plane); weight from l = unit(gauge - P) as the fitted model; sunlit relief (sunlight and
               cast shadows of the frame time, spacecraft shadow) everywhere else; white line = the FOV boundary on the
               visible surface (condition 1 only)
  validate_fixed  the fixed camera and its view; per test frame and gauge: condition (1) against an independent
               instrument-frame rectangle test (getfov half-angles, the frame's rotation), (2), and (3) against exact rays
               (ray_first_hit over all facets) at random and boundary samples; test frames and line-of-sight diagnostic
               images                                                                                       (GPU)
  touchdown    the touchdown (Rosetta impact) point = cops_cuda/figures.impact_point on this mesh with the ground-track
               figures' table (--track-base; tools/pub_figures.py groundtrack's code path), compared with their
               _impact.json (--impact-json); per frame its pixel, the projected surface normal and its visibility from the
               camera (in front, in the image, facet facing the camera, exact first-hit ray camera -> point); check: the
               SPICE trajectory (RORB, ROSETTA and the COPS gauge) intersected with the mesh -> touchdown.npz / .json (GPU)
  compose --touchdown <touchdown.npz>: a generic chequered finish flag (no logo) on a pole planted at the touchdown point
               (pole along the projected surface normal when it is >= FLAG_NORMAL_MIN_DEG from the line of sight, else
               image up), waving (travelling sine wave, phase = frame / FPS in all three movies), anti-aliased, fixed size
               in output pixels, drawn when the point is visible from the camera; legend "touchdown (impact) <UTC> UTC"
2026-11-14 (panorama set; new options, the defaults unchanged):
  plan --camera detach --pan-to <pose.json> [--pan-start-frame K] [--pan-end-frame K]: after the detach the camera moves
               to the final pose of the json (camera_position_m, forward, up, fov_horizontal_deg, pivot_m): smootherstep in
               the frame index, position on an arc around the pivot (direction slerp, distance geometric), orientation by
               quaternion slerp, focal length geometric; checks (clearance, pivot and COPS in the image, steps) in
               camera_checks.json "pan". render / validate_fixed / touchdown / frametable use the pose of each frame.
  compose --no-sub-mark: the sub-spacecraft point (+) is not drawn (user: "remove the white cross")
"""

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import socket
import struct
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

APP = Path(os.environ.get("APP", "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda"))
if str(APP) not in sys.path:
    sys.path.insert(0, str(APP))
from cops_cuda import figures as FIG  # noqa: E402

SCHEMA = "cops_movies_v1"
W, H, SS = 1920, 1080, 2
NFRAMES, FPS = 1800, 30
SHADOW_RES = 24576
AMBIENT, DIFFUSE = 0.16, 0.84
PERSPECTIVE_MAX_DEG = 60.0
MARGIN = 0.92
ZNEAR = 0.05
HOLE_MAX_RAYS = 200000      # enclosed empty pixels (supersampled) ray-cast exactly per frame; above this: not filled, logged
GREY_ZERO = np.array([0.55, 0.55, 0.55])
# observable-part masks (2026-11-13, user: "colored should be only that part of the nucleus that can be observed
# (separately for NG and RG -- use appropriate); use actual aperture of RG"); defaults 'none' = the first movies unchanged
MASK_LINE = np.array([1.0, 1.0, 1.0])     # aperture boundary line (white, alpha MASK_LINE_ALPHA)
MASK_LINE_ALPHA = 0.85
# spacecraft shadow (2026-11-13, user decisions: bus + two solar-array wings only, the last solar-array orientation held
# after the SA C-kernel ends, the Sun as a disc; movies only, not in the fit); default off
R_SUN_M = 6.957e8                         # IAU 2015 nominal solar radius [m]
SA_HOLD_UTC = "2016-09-30T10:39:19.851"   # end of ROS_SA_2016_V0042.BC coverage (ckbrief), orientation held after it
SC_SAMPLES = 64
# detach-200 m movies (2026-11-13, user decisions; new options, defaults unchanged): spacecraft model and FOV cones
HGA_DIAMETER = 2.2             # m: "A 2.2m diameter HGA" (PDS INSTHOST.CAT; ESA); drawn as a flat disc (no depth documented)
SA_PANELS, SA_PANEL_LEN = 5, 2.736   # five panels per wing, substrate 2.25 x 2.736 m (PDS INSTHOST.CAT)
SC_NEAR = 1.0                  # m: spacecraft geometry closer than this to the camera is not drawn (near clipping)
SC_TESS = 0.25                 # m: tessellation of the drawn spacecraft model (perspective depth interpolation)
SC_COL = {"bus": (0.17, 0.17, 0.18), "bus_back": (0.17, 0.17, 0.18), "cells": (0.12, 0.16, 0.40),
          "sa_back": (0.21, 0.21, 0.22), "gap": (0.06, 0.06, 0.07), "hga": (0.15, 0.16, 0.16)}
SC_LINE_COL = (0.72, 0.72, 0.75)     # V-yoke arms and the HGA strut (drawn as lines: their diameters are not documented)
CONE_N = 24                    # generator lines per FOV cone
CONE_COL = {"RG": (1.0, 0.62, 0.11), "NG": (0.31, 0.76, 0.97)}
NO_SUB_MARK = False                 # 2026-11-14 (user): compose --no-sub-mark drops the sub-spacecraft "+" (default: drawn)
CONE_ALPHA = 0.55
CONTACT_TOL = 0.01             # m: the model "touches" the surface when a sample comes this close (unsigned distance)
# fixed camera (2026-11-13, user decisions: the COPS position and the RG / NG boresights when the COPS altitude is 100 m,
# frozen in 67P/C-G_CK; perspective, 90 deg horizontal field of view; the gauges' visible surface from the true gauge
# position: IK FOV, facing the gauge, line of sight); default: the camera at the gauge (the earlier movies)
LOS_RES = 8192                 # gauge line-of-sight depth map: texels on the longer side of its window
LOS_EXTRA_MAX = 200000         # facets cut by the map's near plane, tested exactly for every query (more: refused)
LOS_REFINE_MAX = 400000        # exact boundary rays per gauge and frame (more: the rest keeps the map result; logged)
# touchdown flag (2026-11-13, user: "make the touchdown location a Formula 1 waving finish flag in the movies"): a generic
# chequered flag (no logo), FLAG_NX x FLAG_NY squares on a pole at the touchdown point; sizes in output pixels at
# 1920 x 1080 and text scale 1
FLAG_POLE, FLAG_W, FLAG_H = 64.0, 54.0, 36.0
FLAG_NX, FLAG_NY = 6, 4
FLAG_WAVE_HZ, FLAG_WAVELENGTH, FLAG_AMP = 1.6, 0.9, 0.16  # travelling wave: Hz, wavelength [flag widths], amplitude [heights]
FLAG_NORMAL_MIN_DEG = 20.0     # the pole follows the projected surface normal when it is >= this angle from the line of sight
FLAG_SS = 4                    # supersampling of the flag drawing (anti-aliasing)
LOS_KERNEL = r'''
__device__ __forceinline__ bool seg_hit(const double* __restrict__ V, const long long* __restrict__ Fc, long long f,
                                        double ox, double oy, double oz, double dx, double dy, double dz) {
    // Moller-Trumbore, FP64: does the segment o -> o + d (parameter 0 < t < 1) cross facet f?
    const double* A = V + 3 * Fc[3 * f];
    const double* B = V + 3 * Fc[3 * f + 1];
    const double* C = V + 3 * Fc[3 * f + 2];
    const double e1x = B[0] - A[0], e1y = B[1] - A[1], e1z = B[2] - A[2];
    const double e2x = C[0] - A[0], e2y = C[1] - A[1], e2z = C[2] - A[2];
    const double hx = dy * e2z - dz * e2y, hy = dz * e2x - dx * e2z, hz = dx * e2y - dy * e2x;
    const double det = e1x * hx + e1y * hy + e1z * hz;
    if (det == 0.0) return false;
    const double inv = 1.0 / det;
    const double sx = ox - A[0], sy = oy - A[1], sz = oz - A[2];
    const double u = (sx * hx + sy * hy + sz * hz) * inv;
    if (u < 0.0 || u > 1.0) return false;
    const double qx = sy * e1z - sz * e1y, qy = sz * e1x - sx * e1z, qz = sx * e1y - sy * e1x;
    const double v = (dx * qx + dy * qy + dz * qz) * inv;
    if (v < 0.0 || u + v > 1.0) return false;
    const double t = (e2x * qx + e2y * qy + e2z * qz) * inv;
    return t > 0.0 && t < 1.0;
}
// line of sight from surface points to the gauge: blocked[r] = 1 when the segment orig[r] -> gauge crosses one of the
// front-most facets of the 3 x 3 depth-map texels around (qi[r], qj[r]) or one of the `extra` facets; own[r] skipped
extern "C" __global__ void los_blocked(const double* __restrict__ V, const long long* __restrict__ Fc,
    const unsigned long long* __restrict__ zbuf, const int nx, const int ny,
    const long long* __restrict__ qi, const long long* __restrict__ qj, const double* __restrict__ orig,
    const long long* __restrict__ own, const double gx, const double gy, const double gz,
    const long long* __restrict__ extra, const long long nextra, const long long n,
    unsigned char* __restrict__ blocked, unsigned char* __restrict__ via_map)
{
    const long long r = blockIdx.x * (long long)blockDim.x + threadIdx.x;
    if (r >= n) return;
    const double ox = orig[3 * r], oy = orig[3 * r + 1], oz = orig[3 * r + 2];
    const double dx = gx - ox, dy = gy - oy, dz = gz - oz;
    const long long me = own[r];
    long long seen[9];
    int ns = 0;
    unsigned char b = 0, m = 0;
    for (int dj = -1; dj <= 1 && !b; ++dj) {
        for (int di = -1; di <= 1 && !b; ++di) {
            const long long i = qi[r] + di, j = qj[r] + dj;
            if (i < 0 || j < 0 || i >= nx || j >= ny) continue;
            const unsigned long long key = zbuf[j * (long long)nx + i];
            if (key == 0xffffffffffffffffULL) continue;
            const long long f = (long long)(key & 0xffffffffULL);
            if (f == me) continue;
            bool dup = false;
            for (int s = 0; s < ns; ++s) dup |= (seen[s] == f);
            if (dup) continue;
            seen[ns++] = f;
            if (seg_hit(V, Fc, f, ox, oy, oz, dx, dy, dz)) { b = 1; m = 1; }
        }
    }
    for (long long s = 0; s < nextra && !b; ++s) {
        const long long f = extra[s];
        if (f != me && seg_hit(V, Fc, f, ox, oy, oz, dx, dy, dz)) b = 1;
    }
    blocked[r] = b;
    via_map[r] = m;
}
'''
CMAP = "viridis"
DAY0 = dt.datetime(2016, 9, 30)
KREC = np.dtype([("face", "<i4"), ("Kn", "<f8"), ("KPhi", "<f8"), ("KV", "<f8"), ("cos", "<f8"), ("shadow", "<i4")])
T0 = time.time()


def say(text):
    print("cops_movies [%7.1f s]: %s" % (time.time() - T0, text), flush=True)


def sha256_file(path, bufsize=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(bufsize)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def unit(v):
    v = np.asarray(v, float)
    return v / np.linalg.norm(v)


# ============================================================================================ inputs
def read_geom_lines(path):
    out = {"utc": [], "rows": []}
    for line in open(path):
        t = line.split()
        out["utc"].append(t[0])
        out["rows"].append([float(v) for v in t[1:]])
    a = np.array(out["rows"])
    return {"utc": out["utc"], "et": a[:, 0], "x": a[:, 1:4], "v": a[:, 4:7], "sun": a[:, 7:10], "rh": a[:, 10],
            "b_ng": a[:, 11:14], "b_rg": a[:, 14:17]}


def load_mesh(gdir):
    gdir = Path(gdir)
    nodes = np.fromfile(gdir / "mesh_nodes.f8", "<f8").reshape(-1, 3)
    faces = np.fromfile(gdir / "mesh_faces.i4", "<i4").reshape(-1, 3)
    fg = np.memmap(gdir / "mesh_facegeom.f8", dtype="<f8", mode="r").reshape(-1, 8)
    cent = np.ascontiguousarray(fg[:, 0:3])
    norm = np.ascontiguousarray(fg[:, 3:6])
    area = np.ascontiguousarray(fg[:, 6])
    del fg
    return nodes, faces, cent, norm, area


def read_kernel(path):
    with open(path, "rb") as f:
        raw = f.read()
    version, ip = struct.unpack_from("<ii", raw, 0)
    h = struct.unpack_from("<10d", raw, 8)
    (nk,) = struct.unpack_from("<i", raw, 88)
    rec = np.frombuffer(raw, dtype=KREC, offset=92, count=nk)
    return ip, {"et": h[0], "sc": np.array(h[2:5]), "sun": np.array(h[5:8]), "rh": h[8], "alt": h[9]}, rec


def kernel_header(path):
    with open(path, "rb") as f:
        raw = f.read(92)
    h = struct.unpack_from("<10d", raw, 8)
    return {"et": h[0], "sc": np.array(h[2:5]), "sun": np.array(h[5:8]), "rh": h[8], "alt": h[9]}


# ============================================================================================ prep
def parse_mask(spec):
    """'none' | 'fit' | 'hemisphere' | 'ik_fov' | 'cone:<half-angle deg>' -> (mode, half_angle_deg or None)."""
    if spec.startswith("cone:"):
        return "cone", float(spec.split(":", 1)[1])
    if spec not in ("none", "fit", "hemisphere", "ik_fov"):
        raise SystemExit("unknown mask %s (none, fit, hemisphere, ik_fov, cone:<deg>)" % spec)
    return spec, None


def fov_planes(bounds):
    """normals m_i = B_i x B_(i+1) of the FOV pyramid's side planes (Geometry.cpp InsideFOV), float64 (n, 3)."""
    b = np.asarray(bounds, np.float64)
    return np.cross(b, np.roll(b, -1, axis=0))


def mask_text(gauge, spec, fovinfo=None):
    mode, ang = parse_mask(spec)
    if mode == "ik_fov":
        h = fovinfo["half_angles_deg"][gauge] if fovinfo else None
        return ("the %s IK field of view (%s%s)" % (gauge, "±%g° × ±%g° rectangle" % tuple(h) if h else "SPICE getfov",
                                                   ", SPICE getfov of %s" % fovinfo["id"][gauge] if fovinfo else ""))
    if mode == "hemisphere":
        return "the front half-space of the %s boresight" % gauge
    if mode == "cone":
        return "a cone of %g° half-angle around the %s boresight" % (ang, gauge)
    if mode == "fit":
        return ("the fit's %s acceptance (%s)" % (gauge, "front half-space" if gauge == "RG" else
                                                   "front, or behind through the shield with l·b_RG < 0"))
    return "everywhere"


def cmd_fov(a):
    """IK FOVs of both gauges (getfov, instrument frame) and their boresight + bound vectors in the body frame at every
    frame time and base row time -> <work>/fov.npz; boresights checked against prep.npz (same SPICE calls)."""
    work = Path(a.work)
    P = load_prep(work)
    exe = a.geom_exe

    def run(args, stdin=None):
        r = subprocess.run([exe] + args, input=stdin, capture_output=True, text=True, check=True)
        return r.stdout

    info = {}
    for line in run(["fovinfo", a.kernels]).split("\n"):
        if not line.strip():
            continue
        t = line.split()
        n = int(t[3])
        v = np.array([float(x) for x in t[4:]], np.float64)
        g = "NG" if int(t[0]) == -226331 else "RG"
        info[g] = {"id": int(t[0]), "shape": t[1], "frame": t[2], "n": n, "boresight": v[:3],
                   "bounds": v[3:3 + 3 * n].reshape(n, 3)}
    for g, f in info.items():
        if f["shape"] != "RECTANGLE":
            raise SystemExit("%s FOV shape %s (RECTANGLE expected)" % (g, f["shape"]))
        # half-angles of the rectangle from its corners: angles of a corner projected on the two edge-normal planes
        bs = f["boresight"]
        e1 = f["bounds"][0] - f["bounds"][1]
        e1 -= e1.dot(bs) * bs
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(bs, e1)
        c = f["bounds"][0]
        f["half_angles_deg"] = sorted([float(np.degrees(np.arctan2(abs(c.dot(e)), c.dot(bs)))) for e in (e1, e2)])
        f["corner_angle_deg"] = float(np.degrees(np.arccos(np.clip(f["bounds"].dot(bs), -1, 1).min())))
        f["axes"] = np.array([e1, e2])
        say("%s FOV (%d, %s): %s, %d bounds, half-angles %s deg, corner %.4f deg from the boresight" % (
            g, f["id"], f["frame"], f["shape"], f["n"], ["%.6f" % v for v in f["half_angles_deg"]], f["corner_angle_deg"]))
    out = {}
    for lab, utcs in (("frame", [str(u) for u in P["utc"]]), ("rows", [str(u) for u in P["rows_utc"]])):
        txt = run(["fov", a.kernels], "\n".join(utcs) + "\n")
        L = [ln.split() for ln in txt.split("\n") if ln.strip()]
        if len(L) != len(utcs):
            raise SystemExit("fov: %d lines for %d times" % (len(L), len(utcs)))
        for g_i, g in enumerate(("NG", "RG")):
            n = info[g]["n"]
            bs, bd, rot = [], [], []
            for t in L:
                # utc et | NG: n bs(3) bounds(3n) rot(9) | RG: ...
                off = 2
                for k in range(g_i):
                    off += 1 + 3 + 3 * int(t[off]) + 9
                if int(t[off]) != n:
                    raise SystemExit("bound count changed")
                v = np.array([float(x) for x in t[off + 1: off + 1 + 3 + 3 * n + 9]])
                bs.append(v[:3]); bd.append(v[3:3 + 3 * n].reshape(n, 3)); rot.append(v[3 + 3 * n:].reshape(3, 3))
            out["%s_bs_%s" % (lab, g.lower())] = np.array(bs)
            out["%s_bounds_%s" % (lab, g.lower())] = np.array(bd)
            out["%s_rot_%s" % (lab, g.lower())] = np.array(rot)
        if lab == "frame":
            out["frame_utc"] = np.array(utcs)
    chk = {}
    for lab, pre in (("frame", ""), ("rows", "rows_")):
        for g in ("ng", "rg"):
            chk["%s_boresight_vs_prep_max_abs_%s" % (lab, g)] = float(np.abs(out["%s_bs_%s" % (lab, g)] -
                                                                             P["%sb_%s" % (pre, g)]).max())
    say("boresights vs prep.npz: %s" % json.dumps(chk))
    if max(chk.values()) > 1e-14:
        raise SystemExit("FOV boresights differ from prep.npz")
    np.savez(work / "fov.npz", **out, **{"info_%s_%s" % (g, k): np.asarray(v) for g, f in info.items()
                                         for k, v in f.items()})
    json.dump({"info": {g: {k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in f.items()}
                        for g, f in info.items()}, "checks": chk}, open(work / "fov_checks.json", "w"), indent=1)
    say("fov done")
    return 0


def load_fov(work):
    p = Path(work) / "fov.npz"
    if not p.exists():
        return None
    z = np.load(p, allow_pickle=False)
    F = {k: z[k] for k in z.files}
    info = {"id": {}, "half_angles_deg": {}, "corner_angle_deg": {}, "boresight_cos_corner": {}}
    for g in ("NG", "RG"):
        info["id"][g] = int(F["info_%s_id" % g])
        info["half_angles_deg"][g] = [float(v) for v in F["info_%s_half_angles_deg" % g]]
        info["corner_angle_deg"][g] = float(F["info_%s_corner_angle_deg" % g])
    F["_info"] = info
    return F


def mask_design(gauge, spec, FV):
    """text and the analytic acceptance-weight range of directions inside the mask (for the captions)."""
    mode, ang = parse_mask(spec)
    info = FV["_info"] if FV is not None else None
    out = {"spec": spec, "text": mask_text(gauge, spec, info)}
    out["label"] = ("the %s IK field of view, ±%g° × ±%g°" % ((gauge,) + tuple(info["half_angles_deg"][gauge]))
                    if mode == "ik_fov" else out["text"])
    if mode in ("ik_fov", "cone"):
        cmin = (math.cos(math.radians(info["corner_angle_deg"][gauge])) if mode == "ik_fov"
                else math.cos(math.radians(ang)))
        if gauge == "RG":
            out["weight_range"] = [max(cmin, 0.0), 1.0]
        elif cmin >= 0.0:                     # inside a front cone / rectangle of NG: front only, w = 1 + d.b_NG
            out["weight_range"] = [1.0 + cmin, 2.0]
        out["max_angle_deg"] = info["corner_angle_deg"][gauge] if mode == "ik_fov" else ang
    elif mode == "hemisphere":
        out["weight_range"] = [0.0, 1.0] if gauge == "RG" else [1.0, 2.0]
    return out


def disc_samples(n):
    """n points uniformly over the unit disc (Fibonacci / sunflower), float64 (n, 2)."""
    i = np.arange(n) + 0.5
    r = np.sqrt(i / n)
    th = np.pi * (3.0 - np.sqrt(5.0)) * i
    return np.column_stack([r * np.cos(th), r * np.sin(th)])


def cmd_sc(a):
    """spacecraft geometry for the shadow (movies only): bus corners and solar-array wing corners per frame time, the
    attitude, checks against SPICE -> <work>/sc.npz, sc_checks.json."""
    work = Path(a.work)
    P = load_prep(work)
    utcs = [str(u) for u in P["utc"]]
    r = subprocess.run([a.geom_exe, "scgeom", a.kernels, a.sa_hold_utc], input="\n".join(utcs) + "\n",
                       capture_output=True, text=True, check=True)
    L = [ln.split() for ln in r.stdout.split("\n") if ln.strip()]
    if len(L) != len(utcs):
        raise SystemExit("scgeom: %d lines for %d times" % (len(L), len(utcs)))
    V = np.array([[float(x) for x in t[1:]] for t in L])
    et, et_sa, held = V[:, 0], V[:, 1], V[:, 2].astype(bool)
    o = 3
    xros = V[:, o:o + 3]; o += 3
    R = V[:, o:o + 9].reshape(-1, 3, 3); o += 9
    bus = V[:, o:o + 24].reshape(-1, 8, 3); o += 24
    wing = V[:, o:o + 24].reshape(-1, 8, 3); o += 24
    bus_b = V[:, o:o + 24].reshape(-1, 8, 3); o += 24
    wing_b = V[:, o:o + 24].reshape(-1, 8, 3); o += 24
    sun_sc = V[:, o:o + 3]; o += 3
    rh = V[:, o]
    # checks: corners rotated into the body frame vs SPICE directly; gauge = ROSETTA + R (0,0,1.5) vs prep; Sun
    # the rotation check (the direct values are relative to -226000 in the body frame)
    tb = np.einsum("kij,knj->kni", R, bus)
    tw = np.einsum("kij,knj->kni", R, wing)
    ok = ~held
    # placement: ROS_SPACECRAFT origin (-226000, launch-vehicle interface point) = ROSETTA (-226, centre of gravity)
    # + R origin_offset (option; no loaded kernel gives the offset; default 0 = the lead's specification 2026-11-13)
    off = np.array([float(v) for v in a.origin_offset.split(",")])
    gauge = xros + np.einsum("kij,j->ki", R, np.array([0.0, 0.0, 1.5]))     # the fit's gauge: ROSETTA (CoG) + R (0,0,1.5)
    x_cog = xros.copy()
    xros = xros + np.einsum("kij,j->ki", R, off)
    sun_from_prep = np.einsum("kji,kj->ki", R, np.array([unit(v) for v in P["sun"]]))
    bmin, bmax = bus.min(axis=1), bus.max(axis=1)
    chk = {"origin_offset_sc_frame_m": off.tolist(),
           "note_origin": "ROS_STRUCT_V8.BSP gives the structure relative to ROS_SPACECRAFT (-226000, the launch-vehicle "
                          "interface point); ROSETTA (-226) is the centre of gravity (ROS_V38.TF); no loaded kernel gives "
                          "the offset; -226000 placed at ROSETTA + R origin_offset",
           "bus_corners_vs_spice_max_m": float(np.abs(tb - bus_b).max()),
           "wing_corners_vs_spice_max_m (frames not held)": float(np.abs(tw[ok] - wing_b[ok]).max()) if ok.any() else None,
           "gauge_position_vs_prep_max_m": float(np.abs(gauge - P["x"]).max()),
           "gauge_in_sc_frame_relative_to_-226000_m": (np.array([0.0, 0.0, 1.5]) - off).tolist(),
           "sun_in_sc_frame_vs_prep_max_abs": float(np.abs(sun_from_prep - sun_sc).max()),
           "frames_with_sa_held": int(held.sum()), "first_held_frame": int(np.argmax(held)) if held.any() else None,
           "sa_hold_utc": a.sa_hold_utc,
           "bus_box_sc_frame_m": {"min": bmin[0].tolist(), "max": bmax[0].tolist(),
                                  "corners_are_box_corners": bool(np.all([np.all(np.isclose(np.sort(np.abs(
                                      bus[k] - 0.5 * (bmin[k] + bmax[k])), axis=0), 0.5 * (bmax[k] - bmin[k])[None, :]))
                                      for k in (0, len(bus) - 1)]))},
           "wing_size_m": {g: [float(np.linalg.norm(wing[0, j + 1] - wing[0, j])) for j in (i0, i0 + 1)]
                           for g, i0 in (("+Y", 0), ("-Y", 4))},
           "sun_sc_frame_last": sun_sc[-1].tolist(), "rh_m_first_last": [float(rh[0]), float(rh[-1])],
           "solar_angular_radius_deg_first_last": [math.degrees(R_SUN_M / rh[0]), math.degrees(R_SUN_M / rh[-1])]}
    say("sc checks: %s" % json.dumps(chk))
    if chk["bus_corners_vs_spice_max_m"] > 1e-3 or (chk["wing_corners_vs_spice_max_m (frames not held)"] or 0) > 1e-3 \
            or chk["gauge_position_vs_prep_max_m"] > 1e-3:
        raise SystemExit("spacecraft geometry checks failed")
    extra = {}
    if getattr(a, "hga_ck", None):            # detach-200 m set (2026-11-13): HGA dish and solar-array frames
        extra, chk_h = sc_hga_geometry(a, utcs)
        chk.update(chk_h)
        say("HGA checks: %s" % json.dumps(chk_h))
    np.savez(work / "sc.npz", et=et, et_sa=et_sa, held=held, x_rosetta=xros, x_rosetta_cog=x_cog, R_sc2body=R,
             bus_sc=bus, wing_sc=wing, sun_sc=sun_sc, rh_m=rh, origin_offset=off, **extra)
    json.dump(chk, open(work / "sc_checks.json", "w"), indent=1)
    say("sc done")
    return 0


def hgageom_rows(a, utcs):
    """cops_movies_geom hgageom for the UTC strings: dict of arrays (HGA dish centre, ROS_HGA rotation, gimbal; solar-array
    frame rotations and gimbals), all in ROS_SPACECRAFT relative to -226000."""
    r = subprocess.run([a.geom_exe, "hgageom", a.kernels, a.hga_ck, a.hga_hold_utc, a.sa_hold_utc],
                       input="\n".join(utcs) + "\n", capture_output=True, text=True, check=True)
    L = [ln.split() for ln in r.stdout.split("\n") if ln.strip()]
    if len(L) != len(utcs):
        raise SystemExit("hgageom: %d lines for %d times" % (len(L), len(utcs)))
    V = np.array([[float(x) for x in t[1:]] for t in L])
    o = 4
    out = {"et_hga": V[:, 1], "hga_held": V[:, 2].astype(bool)}
    out["hga_c_sc"] = V[:, o:o + 3]; o += 3
    out["hga_R"] = V[:, o:o + 9].reshape(-1, 3, 3); o += 9
    out["hga_gim_sc"] = V[:, o:o + 3]; o += 3
    out["sa_R_py"] = V[:, o:o + 9].reshape(-1, 3, 3); o += 9
    out["sa_R_my"] = V[:, o:o + 9].reshape(-1, 3, 3); o += 9
    out["sa_gim_py"] = V[:, o:o + 3]; o += 3
    out["sa_gim_my"] = V[:, o:o + 3]
    return out


def sc_hga_geometry(a, utcs):
    """sc step with --hga-ck: the HGA and solar-array frame geometry of every frame time -> extra sc.npz keys + checks."""
    G = hgageom_rows(a, utcs)
    d = np.linalg.norm(G["hga_c_sc"] - G["hga_gim_sc"], axis=1)
    bs = G["hga_R"][:, :, 2]
    chk = {"hga_ck": str(Path(a.hga_ck).resolve()), "hga_ck_sha256": sha256_file(a.hga_ck), "hga_hold_utc": a.hga_hold_utc,
           "hga_frames_held": int(G["hga_held"].sum()), "hga_diameter_m": HGA_DIAMETER,
           "hga_gimbal_to_dish_centre_m_min_max": [float(d.min()), float(d.max())],
           "hga_boresight_last_sc_frame": bs[-1].tolist(),
           "hga_rotation_orthonormal_max_err": float(np.abs(np.einsum("kij,kil->kjl", G["hga_R"], G["hga_R"]) -
                                                            np.eye(3)[None]).max()),
           "sa_gimbals_sc_m": [G["sa_gim_py"][0].tolist(), G["sa_gim_my"][0].tolist()]}
    return G, chk


def sc_model(scz, k, tess=SC_TESS):
    """Spacecraft model of frame k in ROS_SPACECRAFT relative to -226000 (detach-200 m movies, 2026-11-13): the bus box of
    ROS_SPACECRAFT_C1..C8 (outward normals, drawn one-sided), the two solar-array wings (each the rectangle of its SPICE
    corners ROS_SA+-Y_C1..C4, split into SA_PANELS panels of SA_PANEL_LEN along its length with equal gaps; cell side =
    +Z of ROS_SA+-Y; two-sided), the HGA dish (a flat disc of HGA_DIAMETER around ROS_HGA, normal = its boresight +Z;
    two-sided); lines: the V-yoke arms (SA gimbal -> the two inner wing corners) and the HGA strut (gimbal -> dish
    centre). Returns dict V (n,3), T (m,3), part (m,), col_front, col_back (m,3), two_sided (m,), normal (m,3), lines."""
    Vs, Ts, part, cf, cb, two, nrm = [], [], [], [], [], [], []

    def add(Vq, Tq, pid, c_front, c_back, two_sided):
        base = sum(len(v) for v in Vs)
        Vs.append(Vq)
        Ts.append(Tq + base)
        e1, e2 = Vq[Tq[:, 1]] - Vq[Tq[:, 0]], Vq[Tq[:, 2]] - Vq[Tq[:, 0]]
        n = np.cross(e1, e2)
        nrm.append(n / np.linalg.norm(n, axis=1)[:, None])
        part.append(np.full(len(Tq), pid)); cf.append(np.tile(c_front, (len(Tq), 1))); cb.append(np.tile(c_back, (len(Tq), 1)))
        two.append(np.full(len(Tq), two_sided))

    def grid(p0, e1, e2, step):
        n1 = max(1, int(math.ceil(np.linalg.norm(e1) / step)))
        n2 = max(1, int(math.ceil(np.linalg.norm(e2) / step)))
        s_, t_ = np.meshgrid(np.linspace(0, 1, n1 + 1), np.linspace(0, 1, n2 + 1), indexing="ij")
        Vq = p0[None, :] + s_.reshape(-1, 1) * e1[None, :] + t_.reshape(-1, 1) * e2[None, :]
        idx = np.arange((n1 + 1) * (n2 + 1)).reshape(n1 + 1, n2 + 1)
        a_, b_, c_, d_ = idx[:-1, :-1].ravel(), idx[1:, :-1].ravel(), idx[1:, 1:].ravel(), idx[:-1, 1:].ravel()
        return Vq, np.concatenate([np.stack([a_, b_, c_], 1), np.stack([a_, c_, d_], 1)])

    bus = np.asarray(scz["bus_sc"][k], float)
    lo, hi = bus.min(axis=0), bus.max(axis=0)
    L = hi - lo
    for ax in range(3):
        b_, c_ = (ax + 1) % 3, (ax + 2) % 3
        eb, ec = np.zeros(3), np.zeros(3)
        eb[b_], ec[c_] = L[b_], L[c_]
        p_hi, p_lo = lo.copy(), lo.copy()
        p_hi[ax] = hi[ax]
        Vq, Tq = grid(p_hi, eb, ec, tess)                  # + face: eb x ec along +ax
        add(Vq, Tq, 0, SC_COL["bus"], SC_COL["bus_back"], False)
        Vq, Tq = grid(p_lo, ec, eb, tess)                  # - face
        add(Vq, Tq, 0, SC_COL["bus"], SC_COL["bus_back"], False)
    wing = np.asarray(scz["wing_sc"][k], float)
    lines = []
    for wi, (i0, gkey, rkey) in enumerate(((0, "sa_gim_py", "sa_R_py"), (4, "sa_gim_my", "sa_R_my"))):
        w = wing[i0:i0 + 4]
        g = np.asarray(scz[gkey][k], float)
        zc = np.asarray(scz[rkey][k], float)[:, 2]           # cell side
        # w0 -> w1 is the long edge (14.36 m), w1 -> w2 the width (2.25 m); the inner end is the one nearer the gimbal
        i_in = 0 if np.linalg.norm(w[0] - g) < np.linalg.norm(w[1] - g) else 1
        p_inner, p_outer = w[i_in], w[1 - i_in]
        p_inner2 = w[3] if i_in == 0 else w[2]           # the other corner of the inner edge
        el = p_outer - p_inner
        Lw = float(np.linalg.norm(el))
        ul = el / Lw
        gap = (Lw - SA_PANELS * SA_PANEL_LEN) / (SA_PANELS - 1)
        for pnl in range(SA_PANELS):
            s0 = pnl * (SA_PANEL_LEN + gap)
            p0 = p_inner + ul * s0
            e1, e2 = ul * SA_PANEL_LEN, (p_inner2 - p_inner)
            if np.dot(np.cross(e1, e2), zc) < 0:
                p0, e1, e2 = p0 + e1, -e1, e2                    # flip so that the triangle normal is the cell side
            Vq, Tq = grid(p0, e1, e2, tess)
            add(Vq, Tq, 1 + wi, SC_COL["cells"], SC_COL["sa_back"], True)
        lines += [(g, p_inner, SC_LINE_COL), (g, p_inner2, SC_LINE_COL)]
    c_h = np.asarray(scz["hga_c_sc"][k], float)
    Rh = np.asarray(scz["hga_R"][k], float)
    xh, yh, zh = Rh[:, 0], Rh[:, 1], Rh[:, 2]
    rr = np.linspace(0.0, 0.5 * HGA_DIAMETER, max(2, int(math.ceil(0.5 * HGA_DIAMETER / tess))) + 1)
    na = 48
    ang = np.linspace(0, 2 * np.pi, na, endpoint=False)
    Vq = [c_h[None, :]]
    for r_ in rr[1:]:
        Vq.append(c_h[None, :] + r_ * (np.cos(ang)[:, None] * xh[None, :] + np.sin(ang)[:, None] * yh[None, :]))
    Vq = np.vstack(Vq)
    Tq = []
    for j in range(na):
        Tq.append([0, 1 + j, 1 + (j + 1) % na])
    for ri in range(1, len(rr) - 1):
        o0, o1 = 1 + (ri - 1) * na, 1 + ri * na
        for j in range(na):
            j1 = (j + 1) % na
            Tq += [[o0 + j, o1 + j, o1 + j1], [o0 + j, o1 + j1, o0 + j1]]
    Tq = np.array(Tq)
    # orient: normal along the boresight +Z
    e1, e2 = Vq[Tq[0, 1]] - Vq[Tq[0, 0]], Vq[Tq[0, 2]] - Vq[Tq[0, 0]]
    if np.dot(np.cross(e1, e2), zh) < 0:
        Tq = Tq[:, [0, 2, 1]]
    add(Vq, Tq, 3, SC_COL["hga"], SC_COL["hga"], True)
    lines.append((np.asarray(scz["hga_gim_sc"][k], float), c_h, SC_LINE_COL))
    return {"V": np.vstack(Vs), "T": np.vstack(Ts).astype(np.int64), "part": np.concatenate(part),
            "col_front": np.vstack(cf), "col_back": np.vstack(cb), "two_sided": np.concatenate(two),
            "normal": np.vstack(nrm), "lines": lines}


def sc_model_body(scz, k, tess=SC_TESS):
    """sc_model in the body frame: x_rosetta + R_sc2body p."""
    M = sc_model(scz, k, tess)
    R = np.asarray(scz["R_sc2body"][k], float)
    x = np.asarray(scz["x_rosetta"][k], float)
    M["V"] = x[None, :] + M["V"] @ R.T
    M["normal"] = M["normal"] @ R.T
    M["lines"] = [(x + R @ p0, x + R @ p1, c) for p0, p1, c in M["lines"]]
    return M


class ScShadow:
    """Blocked fraction of the solar disc by the spacecraft (bus box + two flat solar-array wings) seen from surface
    points; geometry in the spacecraft frame (bus = axis-aligned box of the 8 corner points; each wing = 2 triangles of
    its 4 corners); the Sun a uniform disc of angular radius R_SUN / r_h sampled with n points (sunflower)."""

    def __init__(self, cp, scz, k, sun_body, nsamp=SC_SAMPLES):
        self.cp = cp
        self.R = np.asarray(scz["R_sc2body"][k], float)
        self.x = np.asarray(scz["x_rosetta"][k], float)
        bus = np.asarray(scz["bus_sc"][k], float)
        self.bmin, self.bmax = bus.min(axis=0), bus.max(axis=0)
        w = np.asarray(scz["wing_sc"][k], float)
        self.tris = [(w[i0], w[i0 + 1], w[i0 + 2]) for i0 in (0, 4)] + [(w[i0], w[i0 + 2], w[i0 + 3]) for i0 in (0, 4)]
        allp = np.vstack([bus, w])
        self.c_sc = 0.5 * (allp.min(axis=0) + allp.max(axis=0))
        self.rb = float(np.linalg.norm(allp - self.c_sc, axis=1).max())
        self.hga = None                        # detach-200 m set (2026-11-13): the HGA dish, a flat 2.2 m disc
        if "hga_c_sc" in getattr(scz, "files", ()):
            c_h = np.asarray(scz["hga_c_sc"][k], float)
            self.hga = (c_h, unit(np.asarray(scz["hga_R"][k], float)[:, 2]), 0.5 * HGA_DIAMETER)
            self.rb_hga = max(self.rb, float(np.linalg.norm(c_h - self.c_sc)) + 0.5 * HGA_DIAMETER)
        self.s_body = unit(sun_body)
        self.s_sc = self.R.T @ self.s_body
        self.alpha = R_SUN_M / float(scz["rh_m"][k])
        e1 = unit(np.cross(self.s_sc, [1.0, 0, 0] if abs(self.s_sc[0]) < 0.9 else [0, 1.0, 0]))
        e2 = np.cross(self.s_sc, e1)
        self.dirs = np.array([unit(self.s_sc + self.alpha * (u * e1 + v * e2)) for u, v in disc_samples(nsamp)])

    def blocked(self, Pb, parts=("bus", "wings")):
        """Pb: (n, 3) cupy surface points (body frame) -> blocked fraction (n,) of the solar disc."""
        cp = self.cp
        o = (Pb - cp.asarray(self.x)[None, :]) @ cp.asarray(self.R)          # = R^T (P - x): spacecraft frame
        # candidates: the central ray passes the bounding sphere (with the disc's angular radius)
        wv = cp.asarray(self.c_sc)[None, :] - o
        tc = wv @ cp.asarray(self.s_sc)
        d2 = cp.sum(wv * wv, axis=1) - tc * tc
        rb = self.rb_hga if ("hga" in parts and self.hga is not None) else self.rb
        cand = cp.nonzero((tc > 0) & (d2 < (rb + 2.0 * self.alpha * cp.abs(tc) + 0.01) ** 2))[0]
        out = cp.zeros(len(Pb))
        if len(cand) == 0:
            return out, 0
        oc = o[cand]
        acc = cp.zeros(len(cand))
        bmin, bmax = cp.asarray(self.bmin), cp.asarray(self.bmax)
        for dvec in self.dirs:
            hit = cp.zeros(len(cand), dtype=bool)
            if "bus" in parts:
                with np.errstate(divide="ignore"):
                    inv = 1.0 / dvec
                t1 = (bmin[None, :] - oc) * cp.asarray(inv)[None, :]
                t2 = (bmax[None, :] - oc) * cp.asarray(inv)[None, :]
                tmin = cp.max(cp.minimum(t1, t2), axis=1)
                tmax = cp.min(cp.maximum(t1, t2), axis=1)
                hit |= (tmax >= cp.maximum(tmin, 0.0)) & (tmax > 0.0)
            if "wings" in parts:
                for v0, v1, v2 in self.tris:
                    e1, e2 = v1 - v0, v2 - v0
                    h = np.cross(dvec, e2)
                    det = float(np.dot(e1, h))
                    if abs(det) < 1e-14:
                        continue
                    sv = oc - cp.asarray(v0)[None, :]
                    u = (sv @ cp.asarray(h)) / det
                    q = cp.cross(sv, cp.asarray(e1)[None, :])
                    vv = (q @ cp.asarray(dvec)) / det
                    t = (q @ cp.asarray(e2)) / det
                    hit |= (u >= 0) & (u <= 1) & (vv >= 0) & (u + vv <= 1) & (t > 0)
            if "hga" in parts and self.hga is not None:
                c_h, n_h, r_h = self.hga
                dn = float(np.dot(dvec, n_h))
                if abs(dn) > 1e-14:
                    t = ((cp.asarray(c_h)[None, :] - oc) @ cp.asarray(n_h)) / dn
                    pnt = oc + t[:, None] * cp.asarray(dvec)[None, :] - cp.asarray(c_h)[None, :]
                    hit |= (t > 0) & (cp.sum(pnt * pnt, axis=1) <= r_h * r_h)
            acc += hit
        out[cand] = acc / len(self.dirs)
        return out, int(len(cand))


def sc_rows(a, utcs):
    """scgeom + hgageom for UTC strings -> dict of arrays as sc.npz (bus_sc, wing_sc, R_sc2body, x_rosetta (offset 0),
    HGA and solar-array keys)."""
    r = subprocess.run([a.geom_exe, "scgeom", a.kernels, a.sa_hold_utc], input="\n".join(utcs) + "\n",
                       capture_output=True, text=True, check=True)
    V = np.array([[float(x) for x in l.split()[1:]] for l in r.stdout.split("\n") if l.strip()])
    o = 3
    out = {"et": V[:, 0], "x_rosetta": V[:, o:o + 3]}
    o += 3
    out["R_sc2body"] = V[:, o:o + 9].reshape(-1, 3, 3); o += 9
    out["bus_sc"] = V[:, o:o + 24].reshape(-1, 8, 3); o += 24
    out["wing_sc"] = V[:, o:o + 24].reshape(-1, 8, 3)
    out.update(hgageom_rows(a, utcs))
    return out


def cmd_contact(a):
    """First contact of the spacecraft model (sc_model: bus box, the two solar-array wings as their five panels, the HGA
    disc; solar arrays and HGA held after their C-kernels end) with the 12M surface on the SPICE trajectory (ROSETTA,
    RORB; structure origin -226000 at ROSETTA): contact at time t when a sample point of the model (5 cm grid) lies on or
    below the surface (signed distance to the nearest of the 48 facets with the nearest centroids, sign by the facet
    normal) or a mesh vertex lies inside the bus box or has crossed a panel or the HGA disc (sign of its distance to the
    plate changed since --t-start, inside the plate outline). Coarse scan (0.1 s, 25 cm model), fine scan (5 ms, 5 cm)
    from 0.2 s before the coarse distance first drops below 0.4 m, bisection to 0.1 ms -> <work>/contact.json.
    The distance used is UNSIGNED (the sign from the nearest facet's normal is not reliable next to cliffs and
    overhangs): contact = the first time a model sample comes within CONTACT_TOL of the surface or a mesh vertex crosses
    a plate / enters the bus; the model is outside the nucleus before (ray parity check of its closest point at
    --t-start); fine-scan steps of 5 ms (< 1 cm of motion at 0.9 m/s) so no sample can pass through unseen."""
    from scipy.spatial import cKDTree
    work = Path(a.work)
    work.mkdir(parents=True, exist_ok=True)
    if (work / "contact.json").exists():
        raise SystemExit("%s/contact.json exists" % work)
    nodes, faces, cent, norm, area = load_mesh(a.geometry)
    tree = cKDTree(cent)
    say("contact: mesh %d facets, KD-tree built" % len(faces))
    e0 = None

    def utc_of(et):
        return (DAY0 + dt.timedelta(seconds=float(et - OFF))).strftime("%Y-%m-%dT%H:%M:%S.%f")

    # et <-> UTC offset (no leap second in the interval): from one SPICE call
    r = subprocess.run([a.geom_exe, "geom", a.kernels, "0", "0", "0"], input=a.t_start + "\n", capture_output=True,
                       text=True, check=True)
    et_s = float(r.stdout.split()[1])
    t_s = dt.datetime.strptime(a.t_start[:19], "%Y-%m-%dT%H:%M:%S")
    OFF = et_s - (t_s - DAY0).total_seconds() - (float("0" + a.t_start[19:]) if len(a.t_start) > 19 else 0.0)

    def signed_min(G, k, tess):
        M = sc_model(G, k, tess)
        R, x = G["R_sc2body"][k], G["x_rosetta"][k]
        Pb = x[None, :] + M["V"] @ R.T
        _, idx = tree.query(Pb, k=48)
        tri = nodes[faces[idx]]
        pc, dist = FIG._closest_on_triangles(Pb, tri[:, :, 0], tri[:, :, 1], tri[:, :, 2])
        j = np.argmin(dist, axis=1)
        rows = np.arange(len(Pb))
        sd = dist[rows, j]                    # unsigned (see the docstring)
        i = int(np.argmin(sd))
        vpart = np.repeat(M["part"], 1)
        # part of the closest vertex: map vertex -> a triangle's part
        vp = np.zeros(len(M["V"]), int)
        vp[M["T"][:, 0]] = M["part"]
        vp[M["T"][:, 1]] = M["part"]
        vp[M["T"][:, 2]] = M["part"]
        return float(sd[i]), int(vp[i]), Pb[i]

    plates_ref = {}

    def vertex_test(G, k, ref=False):
        """mesh vertices inside the bus box or crossing a panel / the HGA disc (sign change since the reference)."""
        R, x = G["R_sc2body"][k], G["x_rosetta"][k]
        near = tree.query_ball_point(x, 40.0)
        if not near:
            return False, None
        vid = np.unique(faces[np.asarray(near)].ravel())
        p = (nodes[vid] - x[None, :]) @ R               # spacecraft frame
        bus = G["bus_sc"][k]
        lo, hi = bus.min(axis=0), bus.max(axis=0)
        inside_bus = np.all((p >= lo) & (p <= hi), axis=1)
        hit = bool(inside_bus.any())
        what = "mesh vertex inside the bus" if hit else None
        M = sc_model(G, k, 10.0)                      # coarse: only the plate outlines are needed
        plates = []
        wing = G["wing_sc"][k]
        for wi, i0 in enumerate((0, 4)):
            w = wing[i0:i0 + 4]
            plates.append(("wing%+d" % (1 - 2 * wi), w[0], w[1] - w[0], w[3] - w[0]))
        c_h, Rh = G["hga_c_sc"][k], G["hga_R"][k]
        for name, p0, e1, e2 in plates:
            n = np.cross(e1, e2)
            n /= np.linalg.norm(n)
            q = p - p0[None, :]
            u = (q @ e1) / np.dot(e1, e1)
            v = (q @ e2) / np.dot(e2, e2)
            wv = q @ n
            inr = (u >= 0) & (u <= 1) & (v >= 0) & (v <= 1)
            key = name
            if ref:
                plates_ref[key] = (vid, np.sign(wv))
            else:
                rv, rs = plates_ref[key]
                pos = dict(zip(rv.tolist(), rs.tolist()))
                s_ref = np.array([pos.get(int(i), 0.0) for i in vid])
                cross = inr & (s_ref != 0) & (np.sign(wv) != s_ref)
                if cross.any():
                    hit, what = True, "mesh vertex crossed %s" % name
        q = p - c_h[None, :]
        wv = q @ Rh[:, 2]
        rad = np.linalg.norm(q - wv[:, None] * Rh[:, 2][None, :], axis=1)
        if ref:
            plates_ref["hga"] = (vid, np.sign(wv))
        else:
            rv, rs = plates_ref["hga"]
            pos = dict(zip(rv.tolist(), rs.tolist()))
            s_ref = np.array([pos.get(int(i), 0.0) for i in vid])
            cross = (rad <= 0.5 * HGA_DIAMETER) & (s_ref != 0) & (np.sign(wv) != s_ref)
            if cross.any():
                hit, what = True, "mesh vertex crossed the HGA disc"
        return hit, what

    G0 = sc_rows(a, [a.t_start])
    sd0, p0_, pt0 = signed_min(G0, 0, 0.05)
    vertex_test(G0, 0, ref=True)
    say("contact: at %s the model is %.3f m above the surface" % (a.t_start, sd0))
    if sd0 <= CONTACT_TOL:
        raise SystemExit("contact: the model already touches the surface at --t-start")
    # the closest model point at --t-start is outside the nucleus: ray parity over all facets (3 directions)
    A_, B_, C_ = nodes[faces[:, 0]], nodes[faces[:, 1]], nodes[faces[:, 2]]
    E1, E2 = B_ - A_, C_ - A_

    def crossings(o, d):
        h = np.cross(d[None, :], E2)
        det = np.einsum("ij,ij->i", E1, h)
        ok = np.abs(det) > 1e-14
        inv = np.where(ok, 1.0 / np.where(ok, det, 1.0), 0.0)
        s_ = o[None, :] - A_
        u_ = np.einsum("ij,ij->i", s_, h) * inv
        q_ = np.cross(s_, E1)
        v_ = (q_ @ d) * inv
        t_ = np.einsum("ij,ij->i", E2, q_) * inv
        return int((ok & (u_ >= 0) & (v_ >= 0) & (u_ + v_ <= 1) & (t_ > 1e-9)).sum())
    par = [crossings(np.asarray(pt0, float), unit(dd)) for dd in ([0, 0, 1.0], [0.3, -0.5, 0.81], [-0.6, 0.2, 0.77])]
    del A_, B_, C_, E1, E2
    if any(c % 2 for c in par):
        raise SystemExit("contact: the closest model point at --t-start is inside the nucleus (crossings %s)" % par)
    say("contact: closest model point at --t-start outside the nucleus (ray crossings %s)" % par)
    # coarse scan
    t = et_s
    tc = None
    while t < et_s + a.max_seconds:
        t += 0.1
        G = sc_rows(a, [utc_of(t)])
        sd, _, _ = signed_min(G, 0, SC_TESS)
        if sd < 0.4:
            tc = t
            break
    if tc is None:
        raise SystemExit("contact: no approach within %g s" % a.max_seconds)
    # fine scan
    t = tc - 0.2
    t_hit = None
    while t < tc + 2.0:
        G = sc_rows(a, [utc_of(t)])
        sd, part, pnt = signed_min(G, 0, 0.05)
        vh, what = vertex_test(G, 0)
        if sd <= CONTACT_TOL or vh:
            t_hit = t
            break
        t += 0.005
    if t_hit is None:
        raise SystemExit("contact: no contact found")
    t_lo, t_hi = t_hit - 0.005, t_hit
    for _ in range(14):
        tm = 0.5 * (t_lo + t_hi)
        G = sc_rows(a, [utc_of(tm)])
        sd, part, pnt = signed_min(G, 0, 0.05)
        vh, what = vertex_test(G, 0)
        if sd <= CONTACT_TOL or vh:
            t_hi = tm
        else:
            t_lo = tm
    G = sc_rows(a, [utc_of(t_hi)])
    sd_hi, part_hi, pnt_hi = signed_min(G, 0, 0.05)
    vh_hi, what_hi = vertex_test(G, 0)
    G = sc_rows(a, [utc_of(t_lo)])
    sd_lo, part_lo, pnt_lo = signed_min(G, 0, 0.05)
    names = {0: "bus", 1: "+Y solar-array wing", 2: "-Y solar-array wing", 3: "HGA dish"}
    res = {"definition": cmd_contact.__doc__, "t_start": a.t_start, "contact_utc": utc_of(t_hi), "contact_et": t_hi,
           "last_free_utc": utc_of(t_lo), "bisection_resolution_s": t_hi - t_lo,
           "min_distance_at_contact_m": sd_hi, "min_distance_before_m": sd_lo, "contact_tolerance_m": CONTACT_TOL,
           "closest_point_at_t_start_ray_crossings": par, "distance_at_t_start_m": sd0,
           "part_at_contact": names.get(part_hi), "point_at_contact_body_m": pnt_hi.tolist(),
           "vertex_test_at_contact": what_hi, "model_sample_spacing_m": 0.05,
           "sa_hold_utc": a.sa_hold_utc, "hga_hold_utc": a.hga_hold_utc, "hga_ck": a.hga_ck,
           "structure_origin": "ROS_SPACECRAFT (-226000) at ROSETTA (-226), offset 0"}
    json.dump(res, open(work / "contact.json", "x"), indent=1, default=float)
    say("contact: %s" % json.dumps({k: v for k, v in res.items() if k != "definition"}))
    return 0


def cmd_prep(a):
    out = Path(a.work)
    out.mkdir(parents=True, exist_ok=True)
    tab = FIG.read_table(a.base)
    alt, et = tab["altitude"], tab["et"]
    if not np.all(np.diff(alt) < 0):
        raise SystemExit("altitude of %s is not strictly decreasing" % a.base)
    obs = (tab["refNG"] > 0) | (tab["refRG"] > 0)
    k_end = int(np.nonzero(obs)[0][-1])
    secs = np.array([(t - DAY0).total_seconds() for t in tab["time"]])
    off = et - secs
    if np.ptp(off) > 1e-5:
        raise SystemExit("et - UTC not constant over the table (%g s)" % np.ptp(off))
    off = float(np.median(off))
    alt_f = np.geomspace(alt[0], alt[k_end], NFRAMES)
    et_f = np.interp(-np.log(alt_f), -np.log(alt[:k_end + 1]), et[:k_end + 1])
    et_f[0], et_f[-1] = et[0], et[k_end]
    utc_f = [(DAY0 + dt.timedelta(seconds=float(e - off))).strftime("%Y-%m-%dT%H:%M:%S.%f") for e in et_f]
    say("frames: %d, altitude %.2f m (%s) -> %.3f m (%s), uniform in log altitude" % (
        NFRAMES, alt_f[0], utc_f[0], alt_f[-1], utc_f[-1]))
    ext = None
    if getattr(a, "extend_to_utc", None):     # detach-200 m set (2026-11-13): frames after the last observation
        tt = dt.datetime.strptime(a.extend_to_utc[:19], "%Y-%m-%dT%H:%M:%S")
        et_end = (tt - DAY0).total_seconds() + (float("0" + a.extend_to_utc[19:]) if len(a.extend_to_utc) > 19 else 0.0) + off
        step = float(et_f[-1] - et_f[-2])
        n_ext = int(math.floor((et_end - et_f[-1]) / step + 1e-9))
        et_x = et_f[-1] + step * np.arange(1, n_ext + 1)
        if et_end - (et_x[-1] if n_ext else et_f[-1]) > 1e-6:
            et_x = np.append(et_x, et_end)
        ext = {"n_obs_frames": len(et_f), "step_s": step, "end_utc": a.extend_to_utc, "frames": len(et_x)}
        et_f = np.concatenate([et_f, et_x])
        utc_f = utc_f + [(DAY0 + dt.timedelta(seconds=float(e - off))).strftime("%Y-%m-%dT%H:%M:%S.%f") for e in et_x]
        say("frames extended: %d frames uniform in time (step %.4f s, the last log-altitude step) after the last "
            "observation, to %s: %d frames in all" % (len(et_x), step, utc_f[-1], len(et_f)))
    (out / "frame_utc.txt").write_text("\n".join(utc_f) + "\n")
    (out / "rows_utc.txt").write_text("\n".join(tab["utc"]) + "\n")
    gtimes = [l.split() for l in open(Path(a.geometry) / "times.txt") if not l.startswith("#") and l.split()]
    gt_map = [(t[0], t[1], np.array([float(v) for v in t[3:6]])) for t in gtimes if t[0].startswith("t")]
    (out / "maps_utc.txt").write_text("\n".join(t[1] for t in gt_map) + "\n")
    off_sc = [float(v) for v in a.gauge_offset.split(",")]
    for name in ("frame", "rows", "maps"):
        with open(out / ("%s_utc.txt" % name)) as fin, open(out / ("%s_geom.txt" % name), "w") as fout:
            subprocess.run([a.geom_exe, "geom", a.kernels] + ["%g" % v for v in off_sc], stdin=fin, stdout=fout,
                           check=True)
    with open(out / "shield_table.txt", "w") as f:
        subprocess.run([a.geom_exe, "shield", "10.3", "33.0", "20.5", "2001", "20000"], stdout=f, check=True)
    fr, rw, mp = (read_geom_lines(out / ("%s_geom.txt" % n)) for n in ("frame", "rows", "maps"))
    # ---- checks of the SPICE geometry against the base (positions: COPS_model.dat; Sun, rh, et: kernel headers)
    chk = {"rows": len(rw["utc"])}
    xt = np.column_stack((tab["x"], tab["y"], tab["z"]))
    chk["gauge_position_vs_COPS_model_dat_max_m"] = float(np.abs(rw["x"] - xt).max())
    chk["et_vs_COPS_model_dat_max_s"] = float(np.abs(rw["et"] - et).max())
    kh = [kernel_header(Path(a.base) / "kernels" / ("point=%04d.bin" % int(p))) for p in tab["iPoint"]]
    # |difference of the unit vectors| ~ angle (acos of a dot product near 1 is limited to ~2e-8 by rounding)
    chk["sun_vs_kernel_headers_max_angle_rad"] = float(max(
        np.linalg.norm(unit(h["sun"]) - rw["sun"][i]) for i, h in enumerate(kh)))
    chk["rh_vs_kernel_headers_max_rel"] = float(max(abs(h["rh"] / rw["rh"][i] - 1) for i, h in enumerate(kh)))
    chk["gauge_position_vs_kernel_headers_max_m"] = float(max(np.abs(h["sc"] - rw["x"][i]).max() for i, h in enumerate(kh)))
    chk["sun_vs_maps_geometry_max_angle_rad"] = float(max(
        np.linalg.norm(unit(s) - mp["sun"][i]) for i, (_, _, s) in enumerate(gt_map)))
    if ext is not None:                       # COPS altitude of the added frames: nearest distance to the shape model
        nodes_m, faces_m = load_mesh(a.geometry)[:2]
        n0 = ext["n_obs_frames"]
        _, d_x, _ = FIG.ground_track(nodes_m, faces_m, np.vstack([fr["x"][n0 - 1:n0], fr["x"][n0:]]))
        ext["last_observation_nearest_distance_m"] = float(d_x[0])
        ext["last_observation_table_altitude_m"] = float(alt_f[-1])
        alt_f = np.concatenate([alt_f, d_x[1:]])
        ext["altitude_m_first_last_added"] = [float(d_x[1]), float(d_x[-1])]
        chk["extension"] = ext
        del nodes_m, faces_m
    sf = np.loadtxt(out / "shield_table.txt")
    chk["shield_table"] = {"bins": len(sf), "SF(0)": float(sf[0]), "SF(0.5)": float(sf[1000]), "SF(1)": float(sf[-1])}
    say("geometry checks: %s" % json.dumps(chk))
    # COPS_model.dat prints positions with 7 significant digits (<= 0.01 m at 15 km); the kernel headers hold the
    # base builder's full-precision values
    if (chk["gauge_position_vs_COPS_model_dat_max_m"] > 0.02 or chk["gauge_position_vs_kernel_headers_max_m"] > 1e-6
            or chk["sun_vs_kernel_headers_max_angle_rad"] > 1e-12 or chk["et_vs_COPS_model_dat_max_s"] > 1e-5):
        raise SystemExit("SPICE geometry does not reproduce the base")
    np.savez(out / "prep.npz", utc=np.array(fr["utc"]), et=fr["et"], x=fr["x"], v=fr["v"], sun=fr["sun"], rh=fr["rh"],
             b_ng=fr["b_ng"], b_rg=fr["b_rg"], alt=alt_f, shield=sf, k_end=k_end,
             maps_labels=np.array([t[0] for t in gt_map]), maps_utc=np.array([t[1] for t in gt_map]),
             maps_sun=mp["sun"], rows_utc=np.array(rw["utc"]), rows_et=rw["et"], rows_x=rw["x"], rows_sun=rw["sun"],
             rows_b_ng=rw["b_ng"], rows_b_rg=rw["b_rg"], rows_ipoint=tab["iPoint"], rows_alt=alt,
             **({} if ext is None else {"n_obs_frames": ext["n_obs_frames"], "extend_step_s": ext["step_s"],
                                        "extend_end_utc": np.array(ext["end_utc"])}))
    json.dump(chk, open(out / "prep_checks.json", "w"), indent=1)
    say("prep done: %s" % out)
    return 0


# ============================================================================================ GPU renderer
class Renderer:
    def __init__(self, gdir, prep, device):
        import cupy as cp
        from cops_cuda import cuda_sources
        self.cp = cp
        self.dev = cp.cuda.Device(int(device))
        self.dev.use()
        self.gdir = Path(gdir)
        nodes, faces, cent, norm, area = load_mesh(gdir)
        self.nf, self.nn = len(faces), len(nodes)
        self.nodes = cp.asarray(nodes)
        self.faces = cp.asarray(faces.astype(np.int64).ravel())
        self.cent = cp.asarray(cent)
        self.norm = cp.asarray(norm)
        self.area_h = area
        self.R = float(np.linalg.norm(nodes, axis=1).max()) * 1.001
        self.idx = cp.arange(self.nf, dtype=cp.uint64)
        self.raster = cuda_sources.function("render", "raster_faces", int(device))
        self.rayhit = cuda_sources.function("render", "ray_first_hit", int(device))
        self.draw_all = cp.ones(self.nf, dtype=cp.uint8)
        sf = prep["shield"]
        self.sf_x = cp.asarray(np.linspace(0.0, 1.0, len(sf)))
        self.sf_y = cp.asarray(sf)
        import matplotlib
        self.lut = cp.asarray(matplotlib.colormaps[CMAP](np.linspace(0, 1, 256))[:, :3])
        del nodes, faces, cent, norm

    def _raster(self, u, v, w, draw, width, height, max_box):
        cp = self.cp
        zbuf = cp.full(width * height, np.iinfo(np.uint64).max, dtype=cp.uint64)
        threads = 256
        blocks = int(min((self.nf + threads - 1) // threads, 65535 * 4))
        self.raster((blocks,), (threads,), (u, v, w, self.faces, draw, np.int64(self.nf), np.int32(width),
                                            np.int32(height), np.int32(max_box), zbuf))
        return zbuf

    def lit(self, sun, res=SHADOW_RES, ktol=1.0, tanmax=10.0):
        """cos SZA per facet and the lit flag (cos > 0 and not in cast shadow) from the shadow map."""
        cp = self.cp
        s = unit(sun)
        a = np.array([1.0, 0, 0]) if abs(s[0]) < 0.9 else np.array([0, 1.0, 0])
        e1 = unit(np.cross(s, a))
        e2 = np.cross(s, e1)
        R = self.R
        texel = 2 * R / res
        cosz = self.norm @ cp.asarray(s)
        u = (self.nodes @ cp.asarray(e1) + R) / texel
        v = (R - self.nodes @ cp.asarray(e2)) / texel
        w = R - self.nodes @ cp.asarray(s)
        zbuf = self._raster(u, v, w, self.draw_all, res, res, res)
        del u, v, w
        cu = cp.floor((self.cent @ cp.asarray(e1) + R) / texel).astype(cp.int64)
        cv = cp.floor((R - self.cent @ cp.asarray(e2)) / texel).astype(cp.int64)
        cu = cp.clip(cu, 0, res - 1)
        cv = cp.clip(cv, 0, res - 1)
        key = zbuf[cv * res + cu]
        del zbuf
        empty = key == cp.uint64(np.iinfo(np.uint64).max)
        fid = key & cp.uint64(0xffffffff)
        zf = (key >> cp.uint64(32)).astype(cp.uint32).view(cp.float32).astype(cp.float64)
        dc = R - self.cent @ cp.asarray(s)
        c = cp.clip(cosz, 1e-6, 1.0)
        tan = cp.minimum(cp.sqrt(1 - c * c) / c, tanmax)
        tol = texel * (ktol + tan)
        lit = (cosz > 0) & (empty | (fid == (self.idx & cp.uint64(0xffffffff))) | (dc <= zf + tol))
        return cosz, lit

    def weights(self, x, b_ng, b_rg):
        cp = self.cp
        l = cp.asarray(np.asarray(x, float))[None, :] - self.cent
        l /= cp.linalg.norm(l, axis=1)[:, None]
        lrg = l @ cp.asarray(unit(b_rg))
        mu = l @ cp.asarray(unit(b_ng))
        w_rg = cp.maximum(0.0, -lrg)
        sf = cp.interp(cp.clip(mu, 0.0, 1.0), self.sf_x, self.sf_y)
        w_ng = cp.where(mu < 0.0, 1.0 - mu, cp.where(lrg < 0.0, sf, 0.0))
        return w_rg, w_ng, l

    def camera(self, x, b_rg, b_ng):
        f = unit(b_rg)
        up = unit(np.asarray(b_ng) - np.dot(b_ng, f) * f)
        right = np.cross(f, up)
        return f, up, right

    def extent(self, x, b_rg, b_ng):
        """angular extent of the nucleus about the camera axis: max theta, and the perspective / fisheye scales that fit."""
        cp = self.cp
        f, up, right = self.camera(x, b_rg, b_ng)
        q = self.nodes - cp.asarray(np.asarray(x, float))[None, :]
        z = q @ cp.asarray(f)
        xx = q @ cp.asarray(right)
        yy = q @ cp.asarray(up)
        rho = cp.sqrt(xx * xx + yy * yy)
        th = cp.arctan2(rho, z)
        thmax = float(th.max())
        out = {"theta_max_deg": math.degrees(thmax)}
        # bounding box of the nucleus in the normalized image coordinates X, Y of each projection (perspective: x/z,
        # y/z; equidistant fisheye: theta cos(phi), theta sin(phi)); the frame is the box centred (principal-point
        # shift: the camera still looks along the RG axis) and scaled to fit with MARGIN
        if float(z.min()) > ZNEAR:
            X, Y = xx / z, yy / z
            out["box_persp"] = [float(X.min()), float(X.max()), float(Y.min()), float(Y.max())]
        else:
            out["box_persp"] = None
        r = cp.where(rho > 0, th / cp.maximum(rho, 1e-300), 0.0)
        X, Y = r * xx, r * yy
        out["box_fish"] = [float(X.min()), float(X.max()), float(Y.min()), float(Y.max())]
        return out

    def project(self, x, b_rg, b_ng, proj, F, center, ss=SS):
        cp = self.cp
        Xc, Yc = center
        f, up, right = self.camera(x, b_rg, b_ng)
        q = self.nodes - cp.asarray(np.asarray(x, float))[None, :]
        z = q @ cp.asarray(f)
        xx = q @ cp.asarray(right)
        yy = q @ cp.asarray(up)
        dist = cp.sqrt(z * z + xx * xx + yy * yy)
        Fs = F * ss
        cx, cy = W * ss / 2.0, H * ss / 2.0
        if proj == "perspective":
            ok_v = z > ZNEAR
            zz = cp.where(ok_v, z, 1.0)
            u = cx + Fs * (xx / zz - Xc)
            v = cy - Fs * (yy / zz - Yc)
        else:
            rho = cp.sqrt(xx * xx + yy * yy)
            th = cp.arctan2(rho, z)
            ok_v = th < math.radians(170.0)
            sc = cp.where(rho > 0, Fs * th / cp.maximum(rho, 1e-300), 0.0)
            u = cx + sc * xx - Fs * Xc
            v = cy - sc * yy + Fs * Yc
        return u, v, dist, ok_v, (f, up, right)

    def view(self, x, b_rg, b_ng, sun, proj, F, center, ktol=1.0, ss=SS):
        """face index per supersampled pixel, per-facet cos, lit, per-facet weights (validation), camera."""
        cp = self.cp
        cosz, lit = self.lit(sun, ktol=ktol)
        u, v, dist, ok_v, basis = self.project(x, b_rg, b_ng, proj, F, center, ss)
        fa = self.faces.reshape(-1, 3)
        okf = ok_v[fa[:, 0]] & ok_v[fa[:, 1]] & ok_v[fa[:, 2]]
        facing = cp.sum(self.norm * (cp.asarray(np.asarray(x, float))[None, :] - self.cent), axis=1) > 0
        draw = (okf & facing).astype(cp.uint8)
        Wp, Hp = W * ss, H * ss
        zbuf = self._raster(u, v, dist, draw, Wp, Hp, Wp // 3)
        del u, v, dist
        empty = zbuf == cp.uint64(np.iinfo(np.uint64).max)
        face = cp.where(empty, -1, (zbuf & cp.uint64(0xffffffff)).astype(cp.int64)).reshape(Hp, Wp)
        del zbuf
        holes = self._fill_holes(face, x, proj, F, center, basis, ss)
        w_rg, w_ng, _ = self.weights(x, b_ng, b_rg)
        return {"face": face, "cos": cosz, "lit": lit, "w_rg": w_rg, "w_ng": w_ng, "drawn": draw, "proj": proj, "F": F,
                "center": center, "basis": basis, "b_ng": unit(b_ng), "b_rg": unit(b_rg), "holes": holes}

    def _pixel_dirs(self, idx, Ws, proj, F, center, basis, ss):
        """world ray directions of the supersampled pixels with flat indices idx (same pixel model as pixel_weights)."""
        cp = self.cp
        f, up, right = basis
        Fs = F * ss
        Xc, Yc = center
        Hs = H * ss
        ii = (idx % Ws).astype(cp.float64) + 0.5
        jj = (idx // Ws).astype(cp.float64) + 0.5
        X = (ii - Ws / 2.0) / Fs + Xc
        Y = (Hs / 2.0 - jj) / Fs + Yc
        if proj == "perspective":
            n = cp.sqrt(X * X + Y * Y + 1.0)
            dr, du, df = X / n, Y / n, 1.0 / n
        else:
            th = cp.sqrt(X * X + Y * Y)
            st = cp.where(th > 0, cp.sin(th) / cp.maximum(th, 1e-300), 1.0)
            dr, du, df = X * st, Y * st, cp.cos(th)
        return (dr[:, None] * cp.asarray(right)[None, :] + du[:, None] * cp.asarray(up)[None, :] +
                df[:, None] * cp.asarray(f)[None, :])

    def _fill_holes(self, face, x, proj, F, center, basis, ss):
        """Empty pixels enclosed by the drawn nucleus (not connected to the frame border) are ray-cast exactly
        (kernels/render.cu ray_first_hit, every facet, FP64) and take the first facet hit. The rasterizer draws
        front-facing facets with straight projected edges, which can leave slivers at grazing incidence close to the
        surface. Pixels whose ray hits nothing stay background. Returns counts (supersampled pixels)."""
        cp = self.cp
        from cupyx.scipy import ndimage as cndi
        empty = face < 0
        out = {"enclosed_px": 0, "filled_px": 0, "filled_backfacing_px": 0, "no_hit_px": 0, "skipped": False}
        if not bool(empty.any()):
            return out
        lab, n = cndi.label(empty)
        edge = cp.unique(cp.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]]))
        hole = empty & ~cp.isin(lab, edge)
        del lab
        idx = cp.nonzero(hole.ravel())[0]
        out["enclosed_px"] = int(idx.size)
        if idx.size == 0:
            return out
        if idx.size > HOLE_MAX_RAYS:
            out["skipped"] = True
            return out
        Hp, Wp = face.shape
        dirs = self._pixel_dirs(idx, Wp, proj, F, center, basis, ss)
        orig = cp.asarray(np.asarray(x, float))
        key = cp.empty(idx.size, dtype=cp.uint64)
        for r0 in range(0, int(idx.size), 512):
            r1 = min(r0 + 512, int(idx.size))
            best = cp.full(r1 - r0, np.iinfo(np.uint64).max, dtype=cp.uint64)
            o = cp.ascontiguousarray(cp.broadcast_to(orig[None, :], (r1 - r0, 3)))
            self.rayhit((256, r1 - r0), (256,), (self.nodes, self.faces, np.int64(self.nf), o,
                                                 cp.ascontiguousarray(dirs[r0:r1]), best))
            key[r0:r1] = best
        hit = key != cp.uint64(np.iinfo(np.uint64).max)
        fh = (key & cp.uint64(0xffffffff)).astype(cp.int64)
        ih = idx[hit]
        face[ih // Wp, ih % Wp] = fh[hit]
        back = cp.sum(self.norm[fh[hit]] * dirs[hit], axis=1) > 0      # facet normal along the ray: seen from behind
        out["filled_px"] = int(hit.sum())
        out["filled_backfacing_px"] = int(back.sum())
        out["no_hit_px"] = int(idx.size - hit.sum())
        return out

    def pixel_weights(self, vw, ss=SS):
        """acceptance weights per supersampled pixel from the pixel's ray direction d (l = -d): exact for the surface
        point seen in the pixel; w_RG = max(0, d.b_RG), w_NG = 1 + d.b_NG (in front of NG, i.e. mu = -d.b_NG < 0), else
        SF(mu) when d.b_RG > 0, else 0."""
        cp = self.cp
        f, up, right = vw["basis"]
        Fs = vw["F"] * ss
        Xc, Yc = vw["center"]
        Hs, Ws = H * ss, W * ss
        ii = cp.arange(Ws, dtype=cp.float64)[None, :] + 0.5
        jj = cp.arange(Hs, dtype=cp.float64)[:, None] + 0.5
        X = (ii - Ws / 2.0) / Fs + Xc
        Y = (Hs / 2.0 - jj) / Fs + Yc
        if vw["proj"] == "perspective":
            n = cp.sqrt(X * X + Y * Y + 1.0)
            dr, du, df = X / n, Y / n, 1.0 / n
        else:
            th = cp.sqrt(X * X + Y * Y)
            st = cp.where(th > 0, cp.sin(th) / cp.maximum(th, 1e-300), 1.0)
            dr, du, df = X * st, Y * st, cp.cos(th)
        bng = vw["b_ng"]
        brg = vw["b_rg"]
        d_bng = dr * float(bng @ right) + du * float(bng @ up) + df * float(bng @ f)
        d_brg = dr * float(brg @ right) + du * float(brg @ up) + df * float(brg @ f)
        w_rg = cp.maximum(0.0, d_brg)
        mu = -d_bng
        sf = cp.interp(cp.clip(mu, 0.0, 1.0).ravel(), self.sf_x, self.sf_y).reshape(mu.shape)
        w_ng = cp.where(mu < 0.0, 1.0 - mu, cp.where(d_brg > 0.0, sf, 0.0))
        return w_rg, w_ng

    def pixel_dirs_world(self, vw, ss=SS):
        """unit ray direction d (gauge -> surface, world frame) of every supersampled pixel (3 arrays Hs x Ws)."""
        cp = self.cp
        f, up, right = vw["basis"]
        Fs = vw["F"] * ss
        Xc, Yc = vw["center"]
        Hs, Ws = H * ss, W * ss
        ii = cp.arange(Ws, dtype=cp.float64)[None, :] + 0.5
        jj = cp.arange(Hs, dtype=cp.float64)[:, None] + 0.5
        X = (ii - Ws / 2.0) / Fs + Xc
        Y = (Hs / 2.0 - jj) / Fs + Yc
        if vw["proj"] == "perspective":
            n = cp.sqrt(X * X + Y * Y + 1.0)
            dr, du, df = X / n, Y / n, 1.0 / n
        else:
            th = cp.sqrt(X * X + Y * Y)
            st = cp.where(th > 0, cp.sin(th) / cp.maximum(th, 1e-300), 1.0)
            dr, du, df = X * st, Y * st, cp.cos(th)
        return tuple(dr * float(right[i]) + du * float(up[i]) + df * float(f[i]) for i in range(3))

    def pixel_points(self, vw, dirs, x):
        """surface point of every nucleus pixel: the pixel ray intersected with its facet's plane (exact for the
        drawn facet). Returns (flat pixel indices, points (n, 3))."""
        cp = self.cp
        face = vw["face"].ravel()
        idx = cp.nonzero(face >= 0)[0]
        f = face[idx]
        d = cp.stack([dirs[i].ravel()[idx] for i in range(3)], axis=1)
        xo = cp.asarray(np.asarray(x, float))
        n = self.norm[f]
        t = cp.sum((self.cent[f] - xo[None, :]) * n, axis=1) / cp.sum(d * n, axis=1)
        return idx, xo[None, :] + d * t[:, None]

    def mask_dirs(self, d, gauge, spec, b_rg, b_ng, bounds=None):
        """observability mask for directions d = (dx, dy, dz) (gauge -> surface; cupy arrays of any shape):
        fit (the fitted model's acceptance > 0), hemisphere, ik_fov (Geometry.cpp InsideFOV with the body-frame bounds
        of this time), cone:<deg>."""
        cp = self.cp
        mode, ang = parse_mask(spec)
        dot = lambda v: d[0] * float(v[0]) + d[1] * float(v[1]) + d[2] * float(v[2])   # noqa: E731
        bg = unit(b_rg if gauge == "RG" else b_ng)
        if mode == "fit":
            d_brg = dot(unit(b_rg))
            if gauge == "RG":
                return d_brg > 0.0
            mu = -dot(unit(b_ng))
            sf = cp.interp(cp.clip(mu, 0.0, 1.0).ravel(), self.sf_x, self.sf_y).reshape(mu.shape)
            return (mu < 0.0) | ((d_brg > 0.0) & (sf > 0.0))
        db = dot(bg)
        if mode == "hemisphere":
            return db > 0.0
        if mode == "cone":
            return db >= math.cos(math.radians(ang))
        if mode == "ik_fov":
            if bounds is None:
                raise SystemExit("ik_fov needs fov.npz (cops_movies.py fov)")
            m = db > 0.0
            for mi in fov_planes(bounds):
                m &= (dot(mi) * float(np.dot(mi, bg))) >= 0.0          # InsideFOV: outside if Dot(m,d)*Dot(m,b) < 0
            return m
        raise SystemExit("mask %s" % spec)

    def colorize_masked(self, vw, which, pw, mask, tint="weight", uniform_color=None, boundary=True, ss=SS,
                        sunfac=None):
        """as colorize, with an observability mask (supersampled bool): inside the mask the tint (the acceptance weight
        on the viridis scale, or one colour), outside sunlit grey; the mask boundary drawn as a thin line."""
        cp = self.cp
        bright = AMBIENT + DIFFUSE * cp.maximum(vw["cos"], 0.0) * vw["lit"]
        w = pw[0] if which == "rg" else pw[1]
        wmax = 1.0 if which == "rg" else 2.0
        if tint == "weight":
            k = cp.clip((w / wmax * 255.0 + 0.5).astype(cp.int64), 0, 255)
            col = self.lut[k]
        else:
            col = cp.broadcast_to(cp.asarray(uniform_color, dtype=cp.float64)[None, None, :], mask.shape + (3,))
        col = cp.where(mask[..., None], col, cp.asarray(GREY_ZERO)[None, None, :])
        face = vw["face"]
        m = face >= 0
        b = cp.zeros(face.shape, dtype=cp.float64)
        if sunfac is None:
            b[m] = bright[face[m]]
        else:
            b[m] = AMBIENT + (bright[face[m]] - AMBIENT) * sunfac[m]
        img = (col * b[..., None]).reshape(H, ss, W, ss, 3).mean(axis=(1, 3))
        if boundary:
            md = mask.reshape(H, ss, W, ss).mean(axis=(1, 3))
            mb = md >= 0.5
            edge = (md > 0.0) & (md < 1.0)
            edge[:, 1:] |= mb[:, 1:] != mb[:, :-1]
            edge[1:, :] |= mb[1:, :] != mb[:-1, :]
            e2 = edge.copy()                              # one-pixel dilation: a 2-3 px line at 1920 x 1080
            e2[:, 1:] |= edge[:, :-1]
            e2[1:, :] |= edge[:-1, :]
            img = cp.where(e2[..., None], img * (1.0 - MASK_LINE_ALPHA) + MASK_LINE_ALPHA * cp.asarray(MASK_LINE), img)
        return cp.asnumpy(cp.clip(img * 255.0 + 0.5, 0, 255).astype(cp.uint8))

    def mask_stats(self, vw, x, which_mask, w, wmax_scale=None):
        """visible-surface area fraction inside the mask (drawn facets, centroid direction), nucleus-pixel fraction, and
        the weight range on nucleus pixels inside the mask."""
        cp = self.cp
        face = vw["face"]
        nuc = face >= 0
        out = {"nucleus_pixel_frac": float(cp.asnumpy((which_mask & nuc).sum() / cp.maximum(nuc.sum(), 1)))}
        sel = which_mask & nuc
        if bool(sel.any()):
            out["weight_min_on_nucleus_in_mask"] = float(cp.asnumpy(w[sel].min()))
            out["weight_max_on_nucleus_in_mask"] = float(cp.asnumpy(w[sel].max()))
        return out

    def colorize(self, vw, which, ss=SS, pw=None, sunfac=None):
        cp = self.cp
        bright = AMBIENT + DIFFUSE * cp.maximum(vw["cos"], 0.0) * vw["lit"]
        w_rg, w_ng = pw if pw is not None else self.pixel_weights(vw, ss)
        w = w_rg if which == "rg" else w_ng
        wmax = 1.0 if which == "rg" else 2.0
        k = cp.clip((w / wmax * 255.0 + 0.5).astype(cp.int64), 0, 255)
        tint = self.lut[k]
        tint = cp.where((w > 0)[..., None], tint, cp.asarray(GREY_ZERO)[None, None, :])
        face = vw["face"]
        m = face >= 0
        b = cp.zeros(face.shape, dtype=cp.float64)
        if sunfac is None:
            b[m] = bright[face[m]]
        else:                                  # spacecraft shadow: direct sunlight scaled per pixel (2026-11-13)
            b[m] = AMBIENT + (bright[face[m]] - AMBIENT) * sunfac[m]
        img = tint * b[..., None]
        img = img.reshape(H, ss, W, ss, 3).mean(axis=(1, 3))
        return cp.asnumpy(cp.clip(img * 255.0 + 0.5, 0, 255).astype(cp.uint8))

    # ---------------------------------------------------------------- fixed camera, line of sight (2026-11-13)
    def maxedge(self):
        """MaxEdge per facet (mesh_facegeom column 7), loaded on first use."""
        if getattr(self, "_maxedge", None) is None:
            fg = np.memmap(self.gdir / "mesh_facegeom.f8", dtype="<f8", mode="r").reshape(-1, 8)
            self._maxedge = self.cp.asarray(np.ascontiguousarray(fg[:, 7]))
            del fg
        return self._maxedge

    def fixed_geometry(self, x, b_rg, b_ng, proj, F, center, ss=SS):
        """camera view of a camera at x looking along b_rg (image up = b_ng projected), the same in every frame of the
        fixed-camera movies: the face buffer exactly as view() draws it (front-facing facets, exact hole fill), the
        nucleus pixels (flat supersampled indices), their surface points (pixel ray x facet plane), facets and normals,
        and the sha256 of the face buffer."""
        cp = self.cp
        u, v, dist, ok_v, basis = self.project(x, b_rg, b_ng, proj, F, center, ss)
        fa = self.faces.reshape(-1, 3)
        okf = ok_v[fa[:, 0]] & ok_v[fa[:, 1]] & ok_v[fa[:, 2]]
        facing = cp.sum(self.norm * (cp.asarray(np.asarray(x, float))[None, :] - self.cent), axis=1) > 0
        draw = (okf & facing).astype(cp.uint8)
        Wp, Hp = W * ss, H * ss
        zbuf = self._raster(u, v, dist, draw, Wp, Hp, Wp // 3)
        del u, v, dist
        empty = zbuf == cp.uint64(np.iinfo(np.uint64).max)
        face = cp.where(empty, -1, (zbuf & cp.uint64(0xffffffff)).astype(cp.int64)).reshape(Hp, Wp)
        del zbuf
        holes = self._fill_holes(face, x, proj, F, center, basis, ss)
        vw = {"face": face, "proj": proj, "F": F, "center": center, "basis": basis, "holes": holes,
              "b_ng": unit(b_ng), "b_rg": unit(b_rg)}
        dirs = self.pixel_dirs_world(vw, ss)
        idx, pts = self.pixel_points(vw, dirs, x)
        del dirs
        fid = face.ravel()[idx]
        vw.update(idx=idx, pts=pts, fid=fid, nrm=self.norm[fid],
                  face_sha256=hashlib.sha256(cp.asnumpy(face).tobytes()).hexdigest())
        return vw

    def point_weights(self, x, b_ng, b_rg, pts):
        """acceptance weights of surface points pts (n, 3) for the gauge at x, as weights() (l = unit vector from the
        point to the gauge): w_RG = max(0, -l.b_RG); mu = l.b_NG: w_NG = 1 - mu if mu < 0, else SF(mu) if l.b_RG < 0,
        else 0."""
        cp = self.cp
        l = cp.asarray(np.asarray(x, float))[None, :] - pts
        l /= cp.linalg.norm(l, axis=1)[:, None]
        lrg = l @ cp.asarray(unit(b_rg))
        mu = l @ cp.asarray(unit(b_ng))
        w_rg = cp.maximum(0.0, -lrg)
        sf = cp.interp(cp.clip(mu, 0.0, 1.0), self.sf_x, self.sf_y)
        w_ng = cp.where(mu < 0.0, 1.0 - mu, cp.where(lrg < 0.0, sf, 0.0))
        return w_rg, w_ng, l

    def _los_kernel(self):
        if getattr(self, "_losk", None) is None:
            self._losk = self.cp.RawKernel(LOS_KERNEL, "los_blocked")
        return self._losk

    def gauge_los(self, x, pts, fid, sel, axis, res=LOS_RES):
        """line of sight from the gauge at x to the surface points pts (n, 3) on facets fid, for the points with sel:
        visible = the segment from the point (+ 1e-6 MaxEdge along its facet normal) to the gauge crosses no facet.
        A perspective depth map from x around `axis` (z along axis, X = x/z, Y = y/z) over the bounding box of the
        selected directions (res texels on its longer side, 2 texels margin; every facet with its three vertices at
        z > ZNEAR, nearest distance wins) gives the front-most facet per texel; the segment is tested exactly
        (Moller-Trumbore, FP64; LOS_KERNEL) against the front-most facets of the point's texel and its 8 neighbours, and
        against every facet cut by the map's near plane that can reach the window (centroid within MaxEdge /
        sin(half the angle between the window's edge and 90 deg) + MaxEdge of the gauge); the point's own facet is
        skipped. Returns (visible (n,) bool, statistics)."""
        cp = self.cp
        n = int(pts.shape[0])
        vis = cp.zeros(n, dtype=bool)
        q = cp.nonzero(sel)[0]
        st = {"queries": int(q.size)}
        if q.size == 0:
            return vis, st
        xg = cp.asarray(np.asarray(x, float))
        f = unit(axis)
        a0 = np.array([1.0, 0, 0]) if abs(f[0]) < 0.9 else np.array([0, 1.0, 0])
        e1 = unit(np.cross(f, a0))
        e2 = np.cross(f, e1)
        d = pts[q] - xg[None, :]
        z = d @ cp.asarray(f)
        if float(z.min()) <= ZNEAR:
            raise SystemExit("line of sight: a selected point is not in front of the gauge's map plane (z %.3g m)"
                             % float(z.min()))
        X = (d @ cp.asarray(e1)) / z
        Y = (d @ cp.asarray(e2)) / z
        del d, z
        x0, x1, y0, y1 = float(X.min()), float(X.max()), float(Y.min()), float(Y.max())
        span = max(x1 - x0, y1 - y0, 1e-9)
        texel = span / (res - 4)
        X0, Y1 = x0 - 2 * texel, y1 + 2 * texel
        nx = int(math.ceil((x1 - x0) / texel)) + 4
        ny = int(math.ceil((y1 - y0) / texel)) + 4
        dn = self.nodes - xg[None, :]
        zn = dn @ cp.asarray(f)
        okn = zn > ZNEAR
        zs = cp.where(okn, zn, 1.0)
        u = ((dn @ cp.asarray(e1)) / zs - X0) / texel
        v = (Y1 - (dn @ cp.asarray(e2)) / zs) / texel
        w = cp.sqrt(cp.sum(dn * dn, axis=1))
        del dn, zs, zn
        fa = self.faces.reshape(-1, 3)
        okf = okn[fa[:, 0]] & okn[fa[:, 1]] & okn[fa[:, 2]]
        zbuf = self._raster(u, v, w, okf.astype(cp.uint8), nx, ny, max(nx, ny) + 1)
        del u, v, w
        # facets cut by the near plane (not in the map) that can reach the window: angular radius asin(MaxEdge / r)
        # around their centroid must span from the window's edge (angle thw from the axis) to ~90 deg
        thw = math.atan(math.hypot(max(abs(X0), abs(X0 + nx * texel)), max(abs(Y1), abs(Y1 - ny * texel))))
        half = max(0.5 * (0.5 * math.pi - thw) - 0.01, 1e-3)
        me = self.maxedge()
        rc = cp.sqrt(cp.sum((self.cent - xg[None, :]) ** 2, axis=1))
        extra = cp.nonzero(~okf & (rc <= me / math.sin(half) + me))[0].astype(cp.int64)
        del rc, okf, okn
        if int(extra.size) > LOS_EXTRA_MAX:
            raise SystemExit("line of sight: %d near-plane facets (more than LOS_EXTRA_MAX)" % int(extra.size))
        qi = cp.floor((X - X0) / texel).astype(cp.int64)
        qj = cp.floor((Y1 - Y) / texel).astype(cp.int64)
        fq = fid[q].astype(cp.int64)
        o = cp.ascontiguousarray(pts[q] + self.norm[fq] * (1e-6 * me[fq])[:, None])
        blocked = cp.zeros(int(q.size), dtype=cp.uint8)
        via = cp.zeros(int(q.size), dtype=cp.uint8)
        threads = 128
        self._los_kernel()(((int(q.size) + threads - 1) // threads,), (threads,),
                           (self.nodes, self.faces, zbuf, np.int32(nx), np.int32(ny), qi, qj, o, fq,
                            np.float64(x[0]), np.float64(x[1]), np.float64(x[2]),
                            extra if int(extra.size) else cp.zeros(1, dtype=cp.int64), np.int64(extra.size),
                            np.int64(q.size), blocked, via))
        vis[q] = blocked == 0
        st.update({"map_px": [nx, ny], "texel_tan": texel, "window_tan": [x0, x1, y0, y1],
                   "window_max_angle_deg": math.degrees(thw), "near_plane_facets_tested": int(extra.size),
                   "blocked": int((blocked != 0).sum()), "blocked_by_map_facets": int((via != 0).sum())})
        del zbuf, qi, qj, o, blocked, via
        return vis, st

    def first_hits_faces(self, orig, dirs, sub):
        """first hit of rays orig (n, 3) + t dirs (n, 3, unit; cupy) among the facets sub (cupy int64 indices; render.cu
        ray_first_hit on that facet list, FP64). Returns (face index or -1, t or inf) as cupy arrays."""
        cp = self.cp
        n = int(orig.shape[0])
        fsub = cp.ascontiguousarray(self.faces.reshape(-1, 3)[sub].ravel())
        ns = int(sub.size)
        key = cp.full(n, np.iinfo(np.uint64).max, dtype=cp.uint64)
        if ns > 0:
            nbx = max(1, min(256, (ns + 255) // 256))
            o = cp.ascontiguousarray(orig)
            d = cp.ascontiguousarray(dirs)
            for r0 in range(0, n, 60000):
                r1 = min(r0 + 60000, n)
                self.rayhit((nbx, r1 - r0), (256,), (self.nodes, fsub, np.int64(ns), o[r0:r1], d[r0:r1], key[r0:r1]))
        hit = key != cp.uint64(np.iinfo(np.uint64).max)
        loc = (key & cp.uint64(0xffffffff)).astype(cp.int64)
        hf = cp.where(hit, sub[cp.where(hit, loc, 0)] if ns > 0 else -1, -1)
        th = cp.where(hit, (key >> cp.uint64(32)).astype(cp.uint32).view(cp.float32).astype(cp.float64), cp.inf)
        return hf, th

    def refine_los(self, x, vw0, fov, seen, max_rays=LOS_REFINE_MAX):
        """exact line of sight at the edges of the unseen parts: every seen camera pixel with an 8-neighbour inside the
        FOV but not seen (hidden, or facet turned away from the gauge) gets the exact segment test (first hit over the
        facets whose centroid lies in the bounding box of the segments and the gauge, widened by the largest MaxEdge:
        every facet that can cross a segment); pixels found hidden become unseen and their neighbours are tested in
        the next iteration, until no pixel changes. The depth-map test can only miss occluders (a blocked verdict is an
        exact segment-facet intersection), so this removes its misses along the edges. Returns (seen, statistics)."""
        cp = self.cp
        idx, pts, fid, nrm = vw0["idx"], vw0["pts"], vw0["fid"], vw0["nrm"]
        Hs, Ws = H * SS, W * SS
        me = self.maxedge()
        mmax = float(cp.asnumpy(me.max()))
        xg = cp.asarray(np.asarray(x, float))
        seen = seen.copy()
        tested = cp.zeros(int(idx.size), dtype=bool)
        st = {"exact_rays": 0, "hidden_by_exact_rays": 0, "iterations": 0, "truncated": False, "facets_tested_max": 0}
        while True:
            dark = cp.zeros(Hs * Ws, dtype=bool)
            dark[idx] = fov & ~seen
            dark = dark.reshape(Hs, Ws)
            nb = dark.copy()
            nb[1:, :] |= dark[:-1, :]
            nb[:-1, :] |= dark[1:, :]
            nb[:, 1:] |= dark[:, :-1]
            nb[:, :-1] |= dark[:, 1:]
            nb[1:, 1:] |= dark[:-1, :-1]
            nb[:-1, :-1] |= dark[1:, 1:]
            nb[1:, :-1] |= dark[:-1, 1:]
            nb[:-1, 1:] |= dark[1:, :-1]
            cand = cp.nonzero(nb.ravel()[idx] & seen & ~tested)[0]
            del dark, nb
            if cand.size == 0:
                break
            room = max_rays - st["exact_rays"]
            if cand.size > room:
                cand = cand[:room]
                st["truncated"] = True
            if cand.size == 0:
                break
            fc = fid[cand]
            o = pts[cand] + nrm[cand] * (1e-6 * me[fc])[:, None]
            dv = xg[None, :] - o
            dd = cp.sqrt(cp.sum(dv * dv, axis=1))
            lo = cp.minimum(o.min(axis=0), xg) - mmax
            hi = cp.maximum(o.max(axis=0), xg) + mmax
            sub = cp.nonzero(cp.all((self.cent >= lo[None, :]) & (self.cent <= hi[None, :]), axis=1))[0]
            hf, th = self.first_hits_faces(o, dv / dd[:, None], sub)
            blk = (hf >= 0) & (th < dd) & (hf != fc)
            seen[cand[blk]] = False
            tested[cand] = True
            st["exact_rays"] += int(cand.size)
            st["hidden_by_exact_rays"] += int(blk.sum())
            st["iterations"] += 1
            st["facets_tested_max"] = max(st["facets_tested_max"], int(sub.size))
            if not bool(blk.any()) or st["truncated"]:
                break
        return seen, st

    def colorize_fixed(self, vw, which, w, mask, fovm, tint="weight", uniform_color=None, boundary=True, ss=SS,
                       sunfac=None):
        """fixed camera: the tint (acceptance weight w per supersampled pixel on the viridis scale, or one colour) where
        mask (the gauge sees the surface point), sunlit grey elsewhere; the FOV boundary on the visible surface (fovm,
        nucleus pixels only; not the limb) as a thin line."""
        cp = self.cp
        bright = AMBIENT + DIFFUSE * cp.maximum(vw["cos"], 0.0) * vw["lit"]
        wmax = 1.0 if which == "rg" else 2.0
        if tint == "weight":
            k = cp.clip((w / wmax * 255.0 + 0.5).astype(cp.int64), 0, 255)
            col = self.lut[k]
        else:
            col = cp.broadcast_to(cp.asarray(uniform_color, dtype=cp.float64)[None, None, :], mask.shape + (3,))
        col = cp.where(mask[..., None], col, cp.asarray(GREY_ZERO)[None, None, :])
        face = vw["face"]
        m = face >= 0
        b = cp.zeros(face.shape, dtype=cp.float64)
        if sunfac is None:
            b[m] = bright[face[m]]
        else:
            b[m] = AMBIENT + (bright[face[m]] - AMBIENT) * sunfac[m]
        img = (col * b[..., None]).reshape(H, ss, W, ss, 3).mean(axis=(1, 3))
        if boundary:
            md = fovm.reshape(H, ss, W, ss).mean(axis=(1, 3))
            nd = m.reshape(H, ss, W, ss).mean(axis=(1, 3))
            full = nd >= 1.0
            mb = md >= 0.5
            edge = (md > 0.0) & (md < nd)
            edge[:, 1:] |= (mb[:, 1:] != mb[:, :-1]) & full[:, 1:] & full[:, :-1]
            edge[1:, :] |= (mb[1:, :] != mb[:-1, :]) & full[1:, :] & full[:-1, :]
            e2 = edge.copy()
            e2[:, 1:] |= edge[:, :-1]
            e2[1:, :] |= edge[:-1, :]
            img = cp.where(e2[..., None], img * (1.0 - MASK_LINE_ALPHA) + MASK_LINE_ALPHA * cp.asarray(MASK_LINE), img)
        return cp.asnumpy(cp.clip(img * 255.0 + 0.5, 0, 255).astype(cp.uint8))

    # ---------------------------------------------------------------- detach-200 m movies (2026-11-13)
    def raster_model(self, cam, basis, F, center, V, T, draw, near=SC_NEAR, ss=SS):
        """perspective z-buffer of a small triangle model (body-frame vertices V (n, 3), triangles T (m, 3), draw (m,))
        for the camera at cam (basis f, up, right): per supersampled pixel the triangle index (-1: none) and its distance
        (linear in screen space; the model is tessellated to SC_TESS). Triangles with a vertex closer than `near` to the
        camera or behind its near plane are not drawn (near clipping)."""
        cp = self.cp
        f, up, right = basis
        q = cp.asarray(np.asarray(V, float) - np.asarray(cam, float)[None, :])
        z = q @ cp.asarray(f)
        xx = q @ cp.asarray(right)
        yy = q @ cp.asarray(up)
        dist = cp.sqrt(cp.sum(q * q, axis=1))
        Fs = F * ss
        Ws_, Hs_ = W * ss, H * ss
        zz = cp.where(z > ZNEAR, z, 1.0)
        u = cp.ascontiguousarray(Ws_ / 2.0 + Fs * (xx / zz - center[0]))
        v = cp.ascontiguousarray(Hs_ / 2.0 - Fs * (yy / zz - center[1]))
        okv = (z > ZNEAR) & (dist >= near)
        Tg = cp.asarray(np.asarray(T, np.int64))
        d = cp.asarray(np.asarray(draw, bool)) & okv[Tg[:, 0]] & okv[Tg[:, 1]] & okv[Tg[:, 2]]
        zbuf = cp.full(Ws_ * Hs_, np.iinfo(np.uint64).max, dtype=cp.uint64)
        nt = int(len(T))
        blocks = int(max(1, min((nt + 255) // 256, 65535 * 4)))
        self.raster((blocks,), (256,), (u, v, cp.ascontiguousarray(dist), cp.ascontiguousarray(Tg.ravel()),
                                         d.astype(cp.uint8), np.int64(nt), np.int32(Ws_), np.int32(Hs_),
                                         np.int32(max(Ws_, Hs_) + 1), zbuf))
        empty = zbuf == cp.uint64(np.iinfo(np.uint64).max)
        tri = cp.where(empty, -1, (zbuf & cp.uint64(0xffffffff)).astype(cp.int64))
        dep = cp.where(empty, cp.inf, (zbuf >> cp.uint64(32)).astype(cp.uint32).view(cp.float32).astype(cp.float64))
        return tri.reshape(Hs_, Ws_), dep.reshape(Hs_, Ws_), int(cp.asnumpy(d.sum()))

    def paint_lines(self, segs, cam, basis, F, center, depth, ss=SS, half=1):
        """depth-tested lines in the supersampled image: segs = [(p0, p1, rgb, alpha)] in the body frame (perspective
        camera); a sample point is drawn where it is nearer than the scene (depth (Hs, Ws), inf = sky) by 1e-4 of its
        distance, as a (2 half + 1)^2 supersampled-pixel dot (half = 1: 1.5 output pixels wide). Returns (alpha, rgb)."""
        cp = self.cp
        f, up, right = (np.asarray(v, float) for v in basis)
        cam = np.asarray(cam, float)
        Hs_, Ws_ = H * ss, W * ss
        A = cp.zeros(Hs_ * Ws_)
        RGB = cp.zeros((Hs_ * Ws_, 3))
        Fs = F * ss
        dflat = depth.ravel()
        for p0, p1, rgb, al in segs:
            p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
            seg = p1 - p0
            Ls = float(np.linalg.norm(seg))
            if Ls <= 0:
                continue
            tt = np.clip(np.dot(cam - p0, seg) / (Ls * Ls), 0.0, 1.0)
            dmin = max(float(np.linalg.norm(p0 + tt * seg - cam)), 0.05)
            n = int(min(400000, 2 + math.ceil(2.0 * Ls / dmin * Fs)))
            t = cp.linspace(0.0, 1.0, n)
            pts = cp.asarray(p0)[None, :] + t[:, None] * cp.asarray(seg)[None, :]
            q = pts - cp.asarray(cam)[None, :]
            z = q @ cp.asarray(f)
            ok = z > ZNEAR
            zz = cp.where(ok, z, 1.0)
            px = Ws_ / 2.0 + Fs * ((q @ cp.asarray(right)) / zz - center[0])
            py = Hs_ / 2.0 - Fs * ((q @ cp.asarray(up)) / zz - center[1])
            dist = cp.sqrt(cp.sum(q * q, axis=1))
            ii = cp.floor(px).astype(cp.int64)
            jj = cp.floor(py).astype(cp.int64)
            for dj in range(-half, half + 1):
                for di in range(-half, half + 1):
                    i2, j2 = ii + di, jj + dj
                    m = ok & (i2 >= 0) & (i2 < Ws_) & (j2 >= 0) & (j2 < Hs_)
                    idx = (j2 * Ws_ + i2)[m]
                    vis = dist[m] < dflat[idx] * (1.0 - 1e-4)
                    idx = idx[vis]
                    A[idx] = cp.maximum(A[idx], al)
                    RGB[idx] = cp.asarray(rgb)[None, :]
        return A.reshape(Hs_, Ws_), RGB.reshape(Hs_, Ws_, 3)

    def colorize_detach(self, vw, which, w, mask, fovm, sc_front, sc_img, lines, tint="weight", uniform_color=None,
                        boundary=True, ss=SS, sunfac=None):
        """detached frames: colorize_fixed's nucleus (tint where the gauge sees the surface, sunlit grey elsewhere), the
        spacecraft model where it is in front (sc_front, sc_img), the depth-tested lines (cone generators, yoke, strut)
        blended in, then the FOV boundary on the visible nucleus (not where the spacecraft covers it)."""
        cp = self.cp
        bright = AMBIENT + DIFFUSE * cp.maximum(vw["cos"], 0.0) * vw["lit"]
        wmax = 1.0 if which == "rg" else 2.0
        if tint == "weight":
            k = cp.clip((w / wmax * 255.0 + 0.5).astype(cp.int64), 0, 255)
            col = self.lut[k]
        else:
            col = cp.broadcast_to(cp.asarray(uniform_color, dtype=cp.float64)[None, None, :], mask.shape + (3,))
        col = cp.where(mask[..., None], col, cp.asarray(GREY_ZERO)[None, None, :])
        face = vw["face"]
        m = face >= 0
        b = cp.zeros(face.shape, dtype=cp.float64)
        if sunfac is None:
            b[m] = bright[face[m]]
        else:
            b[m] = AMBIENT + (bright[face[m]] - AMBIENT) * sunfac[m]
        img = col * b[..., None]
        if sc_front is not None:
            img = cp.where(sc_front[..., None], sc_img, img)
        if lines is not None:
            Al, RGBl = lines
            img = img * (1.0 - Al[..., None]) + Al[..., None] * RGBl
        img = img.reshape(H, ss, W, ss, 3).mean(axis=(1, 3))
        if boundary:
            mm = m if sc_front is None else (m & ~sc_front)
            fm = fovm if sc_front is None else (fovm & ~sc_front)
            md = fm.reshape(H, ss, W, ss).mean(axis=(1, 3))
            nd = mm.reshape(H, ss, W, ss).mean(axis=(1, 3))
            full = nd >= 1.0
            mb = md >= 0.5
            edge = (md > 0.0) & (md < nd)
            edge[:, 1:] |= (mb[:, 1:] != mb[:, :-1]) & full[:, 1:] & full[:, :-1]
            edge[1:, :] |= (mb[1:, :] != mb[:-1, :]) & full[1:, :] & full[:-1, :]
            e2 = edge.copy()
            e2[:, 1:] |= edge[:, :-1]
            e2[1:, :] |= edge[:-1, :]
            img = cp.where(e2[..., None], img * (1.0 - MASK_LINE_ALPHA) + MASK_LINE_ALPHA * cp.asarray(MASK_LINE), img)
        return cp.asnumpy(cp.clip(img * 255.0 + 0.5, 0, 255).astype(cp.uint8))

    def first_hits(self, orig, dirs):
        """first facet hit of rays orig (n, 3) + t dirs (n, 3, unit) over all facets (render.cu ray_first_hit, FP64;
        t stored as float32 in the key). Returns (face or -1, t or inf) as numpy arrays."""
        cp = self.cp
        orig = np.ascontiguousarray(np.asarray(orig, float).reshape(-1, 3))
        dirs = np.ascontiguousarray(np.asarray(dirs, float).reshape(-1, 3))
        key = np.empty(len(orig), np.uint64)
        for r0 in range(0, len(orig), 512):
            r1 = min(r0 + 512, len(orig))
            best = cp.full(r1 - r0, np.iinfo(np.uint64).max, dtype=cp.uint64)
            self.rayhit((256, r1 - r0), (256,), (self.nodes, self.faces, np.int64(self.nf), cp.asarray(orig[r0:r1]),
                                                 cp.asarray(dirs[r0:r1]), best))
            key[r0:r1] = cp.asnumpy(best)
        hit = key != np.iinfo(np.uint64).max
        hf = np.where(hit, (key & np.uint64(0xffffffff)).astype(np.int64), -1)
        th = np.where(hit, (key >> np.uint64(32)).astype(np.uint32).view(np.float32).astype(float), np.inf)
        return hf, th

    def points_visible(self, cams, pts, own):
        """exact visibility of points pts (n, 3) from cams (n, 3): the first facet hit on the ray cam -> point is the
        point's own facet (own >= 0) or lies at the point or beyond (t >= dist (1 - 1e-5) - 1 mm), or nothing is hit.
        Returns (visible, first-hit face, first-hit distance, point distance)."""
        cams = np.asarray(cams, float).reshape(-1, 3)
        pts = np.asarray(pts, float).reshape(-1, 3)
        dv = pts - cams
        dist = np.linalg.norm(dv, axis=1)
        hf, th = self.first_hits(cams, dv / dist[:, None])
        own = np.asarray(own).reshape(-1)
        vis = (hf < 0) | ((hf == own) & (own >= 0)) | (th >= dist * (1.0 - 1e-5) - 1e-3)
        return vis, hf, th, dist


def load_prep(work):
    z = np.load(Path(work) / "prep.npz", allow_pickle=False)
    return {k: z[k] for k in z.files}


# ============================================================================================ plan
def frame_from_box(box):
    """scale F [px per unit of X, Y] and centre (Xc, Yc) of a frame that fits the box with MARGIN."""
    x0, x1, y0, y1 = box
    F = MARGIN * min(W / max(x1 - x0, 1e-12), H / max(y1 - y0, 1e-12))
    return F, (0.5 * (x0 + x1), 0.5 * (y0 + y1))


def choose_projection(ext):
    return ("perspective" if (ext["theta_max_deg"] <= PERSPECTIVE_MAX_DEG and ext["box_persp"] is not None)
            else "fisheye")


def to_pixel(v_cam, proj, F, center):
    """pixel (x right, y down) of a direction given in camera coordinates (right, up, forward)."""
    xx, yy, zz = v_cam
    Xc, Yc = center
    if proj == "perspective":
        if zz <= 0:
            return None
        X, Y = xx / zz, yy / zz
    else:
        rho = math.hypot(xx, yy)
        th = math.atan2(rho, zz)
        X, Y = (th * xx / rho, th * yy / rho) if rho > 0 else (0.0, 0.0)
    return [W / 2 + F * (X - Xc), H / 2 - F * (Y - Yc)]


def to_pixel_np(v_cam, proj, F, center):
    """vectorized to_pixel for an (n, 3) array of camera-frame directions; NaN where not projectable."""
    xx, yy, zz = v_cam[:, 0], v_cam[:, 1], v_cam[:, 2]
    Xc, Yc = center
    if proj == "perspective":
        ok = zz > 0
        zs = np.where(ok, zz, 1.0)
        X, Y = np.where(ok, xx / zs, np.nan), np.where(ok, yy / zs, np.nan)
    else:
        rho = np.hypot(xx, yy)
        th = np.arctan2(rho, zz)
        r = np.where(rho > 0, th / np.maximum(rho, 1e-300), 0.0)
        X, Y = r * xx, r * yy
    return np.column_stack((W / 2 + F * (X - Xc), H / 2 - F * (Y - Yc)))


def cmd_plan(a):
    import cupy as cp
    work = Path(a.work)
    P = load_prep(work)
    rnd = Renderer(a.geometry, P, a.device)
    if getattr(a, "camera", "gauge") == "fixed":          # 2026-11-13: the fixed camera (new option)
        return plan_fixed(a, work, P, rnd)
    NF = len(P["utc"])                      # number of frames (1800; more with prep --extend-to-utc)
    ext = []
    for k in range(NF):
        ext.append(rnd.extent(P["x"][k], P["b_rg"][k], P["b_ng"][k]))
        if k % 300 == 0:
            say("plan: frame %d theta_max %.1f deg" % (k, ext[-1]["theta_max_deg"]))
    th = np.array([e["theta_max_deg"] for e in ext])
    over = np.nonzero([choose_projection(e) == "fisheye" for e in ext])[0]
    k_switch = int(over[0]) if len(over) else NF
    proj = np.array(["perspective" if k < k_switch else "fisheye" for k in range(NF)])
    box = np.array([e["box_persp"] if proj[k] == "perspective" else e["box_fish"] for k, e in enumerate(ext)])
    Freq = np.empty(NF)
    cen = np.empty((NF, 2))
    for k in range(NF):
        Freq[k], cen[k] = frame_from_box(box[k])
    hw = 15
    F = np.empty(NF)
    C = np.empty((NF, 2))
    for s0, s1 in ((0, k_switch), (k_switch, NF)):
        if s1 <= s0:
            continue
        n = s1 - s0
        r = Freq[s0:s1]
        mn = np.array([r[max(0, i - hw):i + hw + 1].min() for i in range(n)])
        F[s0:s1] = [mn[max(0, i - hw):i + hw + 1].mean() for i in range(n)]
        # the frame centre follows the nucleus in every frame (the camera axis stays the RG axis; the attitude changes
        # between frames move the RG axis on the frame, not the nucleus). F: erosion (min over +-hw) then mean over
        # +-hw, which never exceeds the F that fits the frame's own box.
        C[s0:s1] = cen[s0:s1]
    nclip = 0
    for k in range(NF):                    # the smoothed frame must still contain the nucleus
        x0, x1, y0, y1 = box[k]
        hx = max(x1 - C[k, 0], C[k, 0] - x0)
        hy = max(y1 - C[k, 1], C[k, 1] - y0)
        Ffit = MARGIN * min((W / 2) / hx, (H / 2) / hy)
        if F[k] > Ffit * (1 + 1e-12):
            nclip += 1
        F[k] = min(F[k], Ffit)
    say("plan: frames where the smoothed scale had to be reduced to fit the nucleus: %d" % nclip)
    sub, sub_f = FIG.intercept_track(cp.asnumpy(rnd.nodes), cp.asnumpy(rnd.faces).reshape(-1, 3), P["x"], a.device)
    if getattr(a, "camera", "gauge") == "detach":         # 2026-11-13: detach-200 m set (new option)
        return plan_detach_save(a, work, P, proj, F, C, th, box, Freq, cen, k_switch, sub, sub_f, rnd)
    np.savez(work / "plan.npz", proj=proj, F=F, center=C, theta_max_deg=th, box=box, F_req=Freq, center_req=cen,
             k_switch=k_switch, sub=sub, sub_face=sub_f)
    say("plan: perspective for frames 0..%d (theta_max <= %.0f deg), fisheye from frame %d (%s UTC, %.1f m); "
        "sub-spacecraft hits %d / %d" % (k_switch - 1, PERSPECTIVE_MAX_DEG, k_switch,
                                         P["utc"][min(k_switch, NF - 1)], P["alt"][min(k_switch, NF - 1)],
                                         int((sub_f >= 0).sum()), NF))
    return 0


def plan_detach_save(a, work, P, proj, F, C, th, box, Freq, cen, k_switch, sub, sub_f, rnd=None):
    """plan of the detach-200 m movies (2026-11-13, user decisions): frames before the detach frame k_d = the first frame
    whose COPS altitude is <= --detach-altitude use the gauge camera exactly as the gauge-view plan (camera at COPS
    along the RG axis; the same framing: perspective, fisheye from k_switch); from k_d on the camera stays at the COPS
    position of frame k_d with the RG / NG boresights of frame k_d (frozen in 67P/C-G_CK) and a perspective projection
    with --fov-deg horizontal field of view, principal point at the image centre."""
    NF = len(P["utc"])
    alt = P["alt"]
    k_d = int(np.argmax(alt <= a.detach_altitude))
    if not alt[k_d] <= a.detach_altitude or k_d == 0:
        raise SystemExit("detach: no frame at or below %g m" % a.detach_altitude)
    Fd = (W / 2.0) / math.tan(math.radians(a.fov_deg) / 2.0)
    proj = np.array([str(v) for v in proj])
    proj[k_d:] = "perspective"
    F = np.asarray(F, float).copy()
    F[k_d:] = Fd
    C = np.asarray(C, float).copy()
    C[k_d:] = 0.0
    th = np.asarray(th, float).copy()
    th[k_d:] = math.degrees(math.atan(math.hypot(W / 2.0, H / 2.0) / Fd))
    cx = np.asarray(P["x"], float).copy()
    cbr = np.asarray(P["b_rg"], float).copy()
    cbn = np.asarray(P["b_ng"], float).copy()
    cx[k_d:] = P["x"][k_d]
    cbr[k_d:] = P["b_rg"][k_d]
    cbn[k_d:] = P["b_ng"][k_d]
    chk = {"detach_frame": k_d, "detach_utc": str(P["utc"][k_d]), "detach_altitude_m": float(alt[k_d]),
           "previous_frame_altitude_m": float(alt[k_d - 1]), "detach_altitude_target_m": a.detach_altitude,
           "camera_position_m": P["x"][k_d].tolist(), "axis_b_rg": unit(P["b_rg"][k_d]).tolist(),
           "b_ng": unit(P["b_ng"][k_d]).tolist(), "fov_horizontal_deg": a.fov_deg, "F_px": Fd,
           "fov_vertical_deg": 2 * math.degrees(math.atan((H / 2.0) / Fd)), "frames": NF,
           "frames_attached": k_d, "frames_detached": NF - k_d, "gauge_plan_fisheye_from_frame": int(k_switch),
           "gauge_plan_at_detach": {"projection": "fisheye" if k_d >= k_switch else "perspective"}}
    say("detach camera: %s" % json.dumps(chk))
    extra = {}
    if getattr(a, "pan_to", ""):               # 2026-11-14 (user): the camera moves after the detach (new option)
        extra = plan_pan(a, P, rnd, k_d, cx, cbr, cbn, F, th, chk)
    np.savez(work / "plan.npz", proj=proj, F=F, center=C, theta_max_deg=th, box=box, F_req=Freq, center_req=cen,
             k_switch=k_switch, sub=sub, sub_face=sub_f, camera_mode=np.array("detach"), k_detach=k_d,
             detach_utc=np.array(str(P["utc"][k_d])), detach_altitude_m=float(alt[k_d]), cam_fov_deg=float(a.fov_deg),
             cam_x_frames=cx, cam_b_rg_frames=cbr, cam_b_ng_frames=cbn, cam_checks=np.array(json.dumps(chk)), **extra)
    json.dump(chk, open(work / "camera_checks.json", "w"), indent=1)
    say("plan (detach) done: attached frames 0..%d, detached %d..%d" % (k_d - 1, k_d, NF - 1))
    return 0


def quat_from_basis(right, up, f):
    """unit quaternion (w, x, y, z) of the rotation whose columns are (right, up, -f) (camera -> comet frame; right =
    f x up, so (right, up, -f) is right-handed)."""
    m = np.column_stack((right, up, -np.asarray(f, float)))
    if np.linalg.det(m) < 0.999999:
        raise SystemExit("quat_from_basis: not a rotation")
    tr = m[0, 0] + m[1, 1] + m[2, 2]
    if tr > 0:
        S_ = math.sqrt(tr + 1.0) * 2
        q = np.array([0.25 * S_, (m[2, 1] - m[1, 2]) / S_, (m[0, 2] - m[2, 0]) / S_, (m[1, 0] - m[0, 1]) / S_])
    elif m[0, 0] > m[1, 1] and m[0, 0] > m[2, 2]:
        S_ = math.sqrt(1.0 + m[0, 0] - m[1, 1] - m[2, 2]) * 2
        q = np.array([(m[2, 1] - m[1, 2]) / S_, 0.25 * S_, (m[0, 1] + m[1, 0]) / S_, (m[0, 2] + m[2, 0]) / S_])
    elif m[1, 1] > m[2, 2]:
        S_ = math.sqrt(1.0 + m[1, 1] - m[0, 0] - m[2, 2]) * 2
        q = np.array([(m[0, 2] - m[2, 0]) / S_, (m[0, 1] + m[1, 0]) / S_, 0.25 * S_, (m[1, 2] + m[2, 1]) / S_])
    else:
        S_ = math.sqrt(1.0 + m[2, 2] - m[0, 0] - m[1, 1]) * 2
        q = np.array([(m[1, 0] - m[0, 1]) / S_, (m[0, 2] + m[2, 0]) / S_, (m[1, 2] + m[2, 1]) / S_, 0.25 * S_])
    return q / np.linalg.norm(q)


def basis_from_quat(q):
    w, x, y, z = q / np.linalg.norm(q)
    m = np.array([[1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
                  [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
                  [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)]])
    return m[:, 0], m[:, 1], -m[:, 2]


def slerp(a, b, s):
    """spherical linear interpolation of unit vectors / quaternions a -> b (b's sign chosen for the short arc by the
    caller for quaternions)."""
    d = float(np.clip(np.dot(a, b), -1.0, 1.0))
    om = math.acos(d)
    if om < 1e-12:
        return np.asarray(a, float).copy()
    return (math.sin((1 - s) * om) * a + math.sin(s * om) * b) / math.sin(om)


def smootherstep(t):
    t = min(max(t, 0.0), 1.0)
    return t * t * t * (t * (6 * t - 15) + 10)


def plan_pan(a, P, rnd, k_d, cx, cbr, cbn, F, th, chk):
    """camera path after the detach (2026-11-14, user: "after detaching the camera from the spacecraft, slowly move the
    camera and change its orientation such that in the last frame ..."): --pan-to <json> gives the final pose in the
    comet frame (camera_position_m, forward, up, fov_horizontal_deg, pivot_m). From the detach pose (frame k_d) to the
    final pose between frames --pan-start-frame (default k_d) and --pan-end-frame (default the last frame):
    s = smootherstep (6t^5 - 15t^4 + 10t^3, t linear in the frame index: zero speed and acceleration at both ends);
    position on an arc around the pivot (direction from the pivot by slerp, distance geometric in s); orientation by
    quaternion slerp (short arc) of the camera bases; focal length geometric in s (perspective, principal point at the
    image centre). Frames after the end keep the final pose. Fills cx, cbr (= forward), cbn (= up), F, th in place for
    k >= k_d; checks per frame: clearance (nearest surface distance), the pivot and the gauge in the image, steps."""
    import cupy as cp
    spec = json.load(open(a.pan_to))
    NF = len(P["utc"])
    ks = k_d if a.pan_start_frame < 0 else int(a.pan_start_frame)
    ke = NF - 1 if a.pan_end_frame < 0 else int(a.pan_end_frame)
    if not (k_d <= ks < ke <= NF - 1):
        raise SystemExit("pan: need detach frame <= start < end <= last frame (%d %d %d %d)" % (k_d, ks, ke, NF - 1))
    piv = np.asarray(spec["pivot_m"], float)
    x0 = np.asarray(cx[k_d], float)
    f0 = unit(cbr[k_d])
    u0 = unit(np.asarray(cbn[k_d], float) - np.dot(cbn[k_d], f0) * f0)
    r0_ = np.cross(f0, u0)
    x1 = np.asarray(spec["camera_position_m"], float)
    f1 = unit(spec["forward"])
    u1 = unit(np.asarray(spec["up"], float) - np.dot(spec["up"], f1) * f1)
    r1_ = np.cross(f1, u1)
    q0, q1 = quat_from_basis(r0_, u0, f0), quat_from_basis(r1_, u1, f1)
    if np.dot(q0, q1) < 0:
        q1 = -q1
    d0, d1 = x0 - piv, x1 - piv
    R0, R1 = np.linalg.norm(d0), np.linalg.norm(d1)
    n0, n1 = d0 / R0, d1 / R1
    F0 = float(F[k_d])
    F1 = (W / 2.0) / math.tan(math.radians(float(spec["fov_horizontal_deg"])) / 2.0)
    sv = np.zeros(NF)
    for k in range(k_d, NF):
        s = smootherstep((k - ks) / float(ke - ks))
        sv[k] = s
        rr = R0 ** (1 - s) * R1 ** s
        cx[k] = piv + rr * slerp(n0, n1, s)
        r_, u_, f_ = basis_from_quat(slerp(q0, q1, s))
        cbr[k] = f_
        cbn[k] = u_
        F[k] = F0 ** (1 - s) * F1 ** s
        th[k] = math.degrees(math.atan(math.hypot(W / 2.0, H / 2.0) / F[k]))
    # checks
    ksel = np.arange(k_d, NF)
    nodes_h = cp.asnumpy(rnd.nodes)
    faces_h = cp.asnumpy(rnd.faces).reshape(-1, 3)
    _, clear, _ = FIG.ground_track(nodes_h, faces_h, cx[ksel])
    vis_piv, _, _, dpiv = rnd.points_visible(cx[ksel], np.tile(piv, (len(ksel), 1)), np.full(len(ksel), -1))
    pix_piv, pix_g, ang_step, pos_step = [], [], [0.0], [0.0]
    for j, k in enumerate(ksel):
        f_, u_ = unit(cbr[k]), unit(cbn[k])
        r_ = np.cross(f_, u_)
        for pt, lst in ((piv, pix_piv), (np.asarray(P["x"][k], float), pix_g)):
            v = pt - cx[k]
            lst.append(to_pixel((float(v @ r_), float(v @ u_), float(v @ f_)), "perspective", float(F[k]), (0.0, 0.0)))
        if j:
            ang_step.append(math.degrees(math.acos(np.clip(np.dot(f_, unit(cbr[k - 1])), -1, 1))))
            pos_step.append(float(np.linalg.norm(cx[k] - cx[k - 1])))

    def margin(px):
        if px is None:
            return -1e9
        return float(min(px[0], W - px[0], px[1], H - px[1]))
    mp_ = np.array([margin(v) for v in pix_piv])
    far_g = np.linalg.norm(np.asarray(P["x"])[ksel] - cx[ksel], axis=1) > 1.0   # the gauge is the camera at the detach
    mg_ = np.array([margin(v) if fg else np.inf for v, fg in zip(pix_g, far_g)])
    pc = {"spec_file": str(Path(a.pan_to).resolve()), "spec": spec, "start_frame": ks, "end_frame": ke,
          "start_utc": str(P["utc"][ks]), "end_utc": str(P["utc"][ke]), "easing": "smootherstep in the frame index",
          "position": "arc around the pivot: slerp of the direction, distance geometric", "orientation": "quaternion slerp",
          "focal_length": "geometric", "start": {"x": x0.tolist(), "forward": f0.tolist(), "up": u0.tolist(),
                                                 "fov_horizontal_deg": float(chk["fov_horizontal_deg"]),
                                                 "pivot_distance_m": float(R0)},
          "end": {"x": x1.tolist(), "forward": f1.tolist(), "up": u1.tolist(),
                  "fov_horizontal_deg": float(spec["fov_horizontal_deg"]), "pivot_distance_m": float(R1)},
          "rotation_angle_deg": math.degrees(2 * math.acos(min(1.0, abs(float(np.dot(q0, q1)))))),
          "pivot_arc_angle_deg": math.degrees(math.acos(np.clip(np.dot(n0, n1), -1, 1))),
          "mid_frame": int(ks + np.argmin(np.abs(sv[ks:ke + 1] - 0.5))),
          "min_clearance_m": float(clear.min()), "min_clearance_frame": int(ksel[np.argmin(clear)]),
          "pivot_visible_frames": int(vis_piv.sum()), "frames": int(len(ksel)),
          "pivot_not_visible_frames": [int(v) for v in ksel[~vis_piv]][:50],
          "pivot_min_margin_px": float(mp_.min()), "pivot_min_margin_frame": int(ksel[np.argmin(mp_)]),
          "gauge_min_margin_px": float(mg_.min()), "gauge_min_margin_frame": int(ksel[np.argmin(mg_)]),
          "gauge_margin_frames": "frames with COPS more than 1 m from the camera: %d" % int(far_g.sum()),
          "max_axis_step_deg_per_frame": float(max(ang_step)), "max_position_step_m_per_frame": float(max(pos_step)),
          "final_pivot_pixel": pix_piv[-1], "final_gauge_pixel": pix_g[-1]}
    chk["pan"] = pc
    say("pan: %s" % json.dumps({k: v for k, v in pc.items() if k not in ("spec",)}))
    return {"pan_s": sv, "pan_start_frame": ks, "pan_end_frame": ke, "pan_spec": np.array(json.dumps(spec)),
            "pan_clearance_m": np.concatenate([np.full(k_d, np.nan), clear])}


def plan_fixed(a, work, P, rnd):
    """plan of the fixed-camera movies (2026-11-13): the camera = the COPS position and the RG / NG boresights at the
    time t_c at which the COPS altitude of the base table is --camera-altitude (t_c interpolated in -log(altitude) over
    the observation rows, exactly as prep's frame times), from SPICE (cops_movies_geom geom with --gauge-offset, as
    prep), frozen in the body frame; perspective with --fov-deg horizontal field of view, principal point at the image
    centre, the same in every frame. Sub-spacecraft points of the frames as the gauge-camera plan (true positions)."""
    import cupy as cp
    for opt in ("base", "kernels", "geom_exe"):
        if not getattr(a, opt, None):
            raise SystemExit("--camera fixed needs --%s" % opt.replace("_", "-"))
    tab = FIG.read_table(a.base)
    alt, et = tab["altitude"], tab["et"]
    k_end = int(P["k_end"])
    h = float(a.camera_altitude)
    if not (alt[k_end] < h < alt[0]):
        raise SystemExit("--camera-altitude %g m outside the descent (%g .. %g m)" % (h, alt[k_end], alt[0]))
    et_c = float(np.interp(-math.log(h), -np.log(alt[:k_end + 1]), et[:k_end + 1]))
    secs = np.array([(t - DAY0).total_seconds() for t in tab["time"]])
    off = float(np.median(et - secs))
    utc_c = (DAY0 + dt.timedelta(seconds=float(et_c - off))).strftime("%Y-%m-%dT%H:%M:%S.%f")
    off_sc = [float(v) for v in a.gauge_offset.split(",")]
    r = subprocess.run([a.geom_exe, "geom", a.kernels] + ["%g" % v for v in off_sc], input=utc_c + "\n",
                       capture_output=True, text=True, check=True)
    t = r.stdout.split()
    v = np.array([float(s) for s in t[1:]])
    cam = {"utc": t[0], "et": float(v[0]), "x": v[1:4], "sun": v[7:10], "b_ng": v[11:14], "b_rg": v[14:17]}
    i1 = int(np.searchsorted(et, et_c))
    i0 = i1 - 1
    fr = (et_c - et[i0]) / (et[i1] - et[i0])
    xt = np.column_stack((tab["x"], tab["y"], tab["z"]))
    x_tab = xt[i0] + fr * (xt[i1] - xt[i0])
    nodes_h = cp.asnumpy(rnd.nodes)
    faces_h = cp.asnumpy(rnd.faces).reshape(-1, 3)
    _, nd, nf_ = FIG.ground_track(nodes_h, faces_h, cam["x"][None, :])
    f, up, right = rnd.camera(cam["x"], cam["b_rg"], cam["b_ng"])
    F = (W / 2.0) / math.tan(math.radians(a.fov_deg) / 2.0)
    hf, th = rnd.first_hits(cam["x"][None, :], f[None, :])
    axis_point = None
    if hf[0] >= 0:
        n0, c0 = cp.asnumpy(rnd.norm[int(hf[0])]), cp.asnumpy(rnd.cent[int(hf[0])])
        ta = float(np.dot(c0 - cam["x"], n0) / np.dot(f, n0))
        axis_point = (cam["x"] + ta * f).tolist()
    sub, sub_f = FIG.intercept_track(nodes_h, faces_h, P["x"], a.device)
    kf = int(np.argmin(np.abs(P["et"] - et_c)))
    chk = {"camera_altitude_target_m": h, "camera_utc": cam["utc"], "camera_et": cam["et"],
           "camera_et_interpolated": et_c, "et_spice_minus_interpolated_s": cam["et"] - et_c,
           "rows_bracketing": [str(tab["utc"][i0]), str(tab["utc"][i1])],
           "camera_position_m": cam["x"].tolist(),
           "camera_vs_table_position_linear_interp_m": float(np.linalg.norm(cam["x"] - x_tab)),
           "camera_nearest_surface_distance_m": float(nd[0]), "camera_nearest_face": int(nf_[0]),
           "gauge_offset_sc_m": off_sc, "axis_b_rg": unit(cam["b_rg"]).tolist(), "b_ng": unit(cam["b_ng"]).tolist(),
           "image_up_angle_from_b_ng_deg": math.degrees(math.acos(np.clip(np.dot(up, unit(cam["b_ng"])), -1, 1))),
           "fov_horizontal_deg": a.fov_deg, "F_px": F,
           "fov_vertical_deg": 2 * math.degrees(math.atan((H / 2.0) / F)),
           "fov_diagonal_half_deg": math.degrees(math.atan(math.hypot(W / 2.0, H / 2.0) / F)),
           "axis_surface_point_m": axis_point, "axis_surface_distance_m": float(th[0]) if hf[0] >= 0 else None,
           "axis_face": int(hf[0]),
           "nadir_angle_of_axis_deg": math.degrees(math.acos(np.clip(np.dot(f, -unit(cam["x"])), -1, 1))),
           "nearest_frame": kf, "nearest_frame_utc": str(P["utc"][kf]), "nearest_frame_altitude_m": float(P["alt"][kf]),
           "frames_before_camera_time": int((P["et"] < cam["et"]).sum())}
    say("fixed camera: %s" % json.dumps(chk))
    if abs(chk["et_spice_minus_interpolated_s"]) > 1e-5 or chk["camera_vs_table_position_linear_interp_m"] > 0.1:
        raise SystemExit("fixed camera: SPICE does not reproduce the base table at the camera time")
    hd = chk["fov_diagonal_half_deg"]
    NF = len(P["utc"])
    np.savez(work / "plan.npz", proj=np.array(["perspective"] * NF), F=np.full(NF, F),
             center=np.zeros((NF, 2)), theta_max_deg=np.full(NF, hd), k_switch=NF, sub=sub,
             sub_face=sub_f, camera_mode=np.array("fixed"), cam_utc=np.array(cam["utc"]), cam_et=cam["et"],
             cam_x=cam["x"], cam_b_rg=cam["b_rg"], cam_b_ng=cam["b_ng"], cam_f=f, cam_up=up, cam_right=right,
             cam_altitude_target_m=h, cam_fov_deg=float(a.fov_deg), cam_checks=np.array(json.dumps(chk)))
    json.dump(chk, open(work / "camera_checks.json", "w"), indent=1)
    say("plan (fixed camera) done: sub-spacecraft hits %d / %d" % (int((sub_f >= 0).sum()), NF))
    return 0


def is_fixed(Q):
    return "camera_mode" in Q.files and str(Q["camera_mode"]) == "fixed"


def is_detach(Q):
    return "camera_mode" in Q.files and str(Q["camera_mode"]) == "detach"


def frame_camera(P, Q, k):
    """camera position and (b_rg, b_ng) defining its basis in frame k: the gauge of the frame, or the fixed camera."""
    if is_detach(Q):
        return (np.asarray(Q["cam_x_frames"][k], float), np.asarray(Q["cam_b_rg_frames"][k], float),
                np.asarray(Q["cam_b_ng_frames"][k], float))
    if is_fixed(Q):
        return np.asarray(Q["cam_x"], float), np.asarray(Q["cam_b_rg"], float), np.asarray(Q["cam_b_ng"], float)
    return np.asarray(P["x"][k], float), np.asarray(P["b_rg"][k], float), np.asarray(P["b_ng"][k], float)


# ============================================================================================ validate
def cmd_validate(a):
    import cupy as cp
    from PIL import Image
    work = Path(a.work)
    P = load_prep(work)
    Q = np.load(work / "plan.npz", allow_pickle=False)
    rnd = Renderer(a.geometry, P, a.device)
    area = rnd.area_h
    res = {"shadow_map_vs_maps_geometry": [], "kernel_records": []}
    nf = rnd.nf
    # ---- V1: shadow map vs the exact (BVH) cast shadows of the maps geometry at its 13 times
    for i, lab in enumerate(P["maps_labels"]):
        cos_x = np.fromfile(Path(a.geometry) / ("%s.cos.f8" % lab), "<f8")
        sh_x = np.fromfile(Path(a.geometry) / ("%s.shadow.u1" % lab), np.uint8) != 0
        for resn, ktol in ((16384, 1.0), (16384, 2.0), (16384, 0.5), (24576, 1.0)):
            t0 = time.time()
            cosz, lit = rnd.lit(P["maps_sun"][i], res=resn, ktol=ktol)
            dt_s = time.time() - t0
            cosz, lit = cp.asnumpy(cosz), cp.asnumpy(lit)
            day = cos_x > 0
            ours_sh = day & ~lit
            ex_sh = day & sh_x
            disagree = day & (ours_sh != ex_sh)
            res["shadow_map_vs_maps_geometry"].append({
                "label": str(lab), "utc": str(P["maps_utc"][i]), "res": resn, "ktol": ktol, "seconds": dt_s,
                "cos_max_abs_diff": float(np.abs(cosz - cos_x).max()), "day_facets": int(day.sum()),
                "exact_shadowed": int(ex_sh.sum()), "ours_shadowed": int(ours_sh.sum()),
                "false_lit": int((ex_sh & ~ours_sh).sum()), "false_shadow": int((ours_sh & ~ex_sh).sum()),
                "disagree_frac_of_day_facets": float(disagree.sum() / day.sum()),
                "disagree_frac_of_day_area": float(area[disagree].sum() / area[day].sum())})
            r = res["shadow_map_vs_maps_geometry"][-1]
            say("V1 %s %s res %d ktol %g: day %d, shadowed exact %d ours %d, disagree %.3f %% facets %.3f %% area (%.2f s)" % (
                lab, r["utc"], resn, ktol, r["day_facets"], r["exact_shadowed"], r["ours_shadowed"],
                100 * r["disagree_frac_of_day_facets"], 100 * r["disagree_frac_of_day_area"], dt_s))
    # ---- V2: kernel records of base rows: cos, shadow, weights, visibility in the rendered view
    rows = [int(v) for v in a.rows.split(",")]
    rnd_cent_h = cp.asnumpy(rnd.cent)
    ip_all = list(P["rows_ipoint"].astype(int))
    for ip in rows:
        i = ip_all.index(ip)
        _, hdr, rec = read_kernel(Path(a.base) / "kernels" / ("point=%04d.bin" % ip))
        x, sun, b_ng, b_rg = P["rows_x"][i], P["rows_sun"][i], P["rows_b_ng"][i], P["rows_b_rg"][i]
        ext = rnd.extent(x, b_rg, b_ng)
        proj = choose_projection(ext)
        F, center = frame_from_box(ext["box_persp"] if proj == "perspective" else ext["box_fish"])
        vw = rnd.view(x, b_rg, b_ng, sun, proj, F, center)
        face_img = cp.asnumpy(vw["face"])
        # in-frame test of the record facets: centroid direction -> pixel
        fcam = np.column_stack([(rnd_cent_h[rec["face"]] - x) @ v for v in (vw["basis"][2], vw["basis"][1], vw["basis"][0])])
        drawn = np.zeros(nf, bool)
        drawn[np.unique(face_img[face_img >= 0])] = True
        cosz, lit = cp.asnumpy(vw["cos"]), cp.asnumpy(vw["lit"])
        w_rg, w_ng = cp.asnumpy(vw["w_rg"]), cp.asnumpy(vw["w_ng"])
        f = rec["face"].astype(np.int64)
        rc = {"iPoint": ip, "utc": str(P["rows_utc"][i]), "altitude_m": float(P["rows_alt"][i]), "records": int(len(rec)),
              "projection": proj, "theta_max_deg": ext["theta_max_deg"]}
        rc["cos_max_abs_diff"] = float(np.abs(cosz[f] - rec["cos"]).max())
        day = rec["cos"] > 0
        rsh = day & (rec["shadow"] != 0)
        osh = day & ~lit[f]
        rc["shadow_agree_frac_day_records"] = float(np.mean(rsh[day] == osh[day])) if day.any() else None
        rc["shadow_agree_frac_Kn_weighted"] = (float(rec["Kn"][day][rsh[day] == osh[day]].sum() / rec["Kn"][day].sum())
                                               if day.any() else None)
        m = rec["KV"] > 0
        rc["RG_weight_KPhi_over_KV_vs_ours"] = {
            "n": int(m.sum()), "median_abs": float(np.median(np.abs(rec["KPhi"][m] / rec["KV"][m] - w_rg[f][m]))),
            "p99_abs": float(np.percentile(np.abs(rec["KPhi"][m] / rec["KV"][m] - w_rg[f][m]), 99))}
        rel = np.abs(rec["Kn"][m] / rec["KV"][m] / np.maximum(w_ng[f][m], 1e-300) - 1)
        rc["NG_weight_Kn_over_KV_vs_ours"] = {"n": int(m.sum()), "median_rel": float(np.median(rel)),
                                              "p99_rel": float(np.percentile(rel, 99))}
        rc["records_Kn_pos_with_our_w_ng_pos_frac"] = float(np.mean(w_ng[f][rec["Kn"] > 0] > 0))
        rc["records_KV_pos_with_our_w_rg_pos_frac"] = float(np.mean(w_rg[f][rec["KV"] > 0] > 0))
        pix = to_pixel_np(fcam, proj, F, center)
        inside = np.isfinite(pix[:, 0]) & (pix[:, 0] >= 0) & (pix[:, 0] < W) & (pix[:, 1] >= 0) & (pix[:, 1] < H)
        for key in ("Kn", "KPhi"):
            mm = rec[key] > 0
            rc["%s_in_frame_frac" % key] = float(rec[key][mm][inside[mm]].sum() / rec[key][mm].sum())
            rc["%s_on_drawn_facets_frac" % key] = float(rec[key][mm][drawn[f][mm]].sum() / rec[key][mm].sum())
            rc["%s_records_drawn_frac" % key] = float(drawn[f][mm].mean())
        has_rec = np.zeros(nf, bool)
        has_rec[f] = True
        rc["drawn_facets"] = int(drawn.sum())
        rc["drawn_facets_with_record_frac"] = float(has_rec[drawn].mean())
        res["kernel_records"].append(rc)
        say("V2 row %d %s (%.1f m, %s): %s" % (ip, rc["utc"], rc["altitude_m"], proj, json.dumps(
            {k: v for k, v in rc.items() if k not in ("iPoint", "utc", "altitude_m", "projection")})))
        del vw
    # ---- test frames (clean views) for the visual check and the split check
    td = work / "test_frames"
    td.mkdir(exist_ok=True)
    for k in [int(v) for v in a.test_frames.split(",")]:
        vw = rnd.view(P["x"][k], P["b_rg"][k], P["b_ng"][k], P["sun"][k], str(Q["proj"][k]), float(Q["F"][k]),
                      tuple(Q["center"][k]))
        pw = rnd.pixel_weights(vw)
        for which in ("rg", "ng"):
            Image.fromarray(rnd.colorize(vw, which, pw=pw)).save(td / ("f%05d_%s.png" % (k, which)), compress_level=1)
        say("test frame %d written; enclosed empty pixels (supersampled): %s" % (k, json.dumps(vw["holes"])))
        res.setdefault("test_frame_holes", {})[str(k)] = vw["holes"]
        del vw
    json.dump(res, open(work / ("validation_%s.json" % socket.gethostname()), "w"), indent=1, default=float)
    say("validate done")
    return 0


# ============================================================================================ render
def cmd_render(a):
    from PIL import Image
    work = Path(a.work)
    P = load_prep(work)
    Q = np.load(work / "plan.npz", allow_pickle=False)
    rnd = Renderer(a.geometry, P, a.device)
    specs = {"RG": a.mask_rg, "NG": a.mask_ng}
    masked = any(v != "none" for v in specs.values())
    for v in specs.values():
        parse_mask(v)
    FV = load_fov(work) if any(parse_mask(v)[0] == "ik_fov" for v in specs.values()) else None
    if masked and any(parse_mask(v)[0] == "ik_fov" for v in specs.values()) and FV is None:
        raise SystemExit("ik_fov needs %s/fov.npz (cops_movies.py fov)" % work)
    tints = {"RG": a.tint_rg, "NG": a.tint_ng}
    SCZ = None
    if a.sc_shadow:
        if not (work / "sc.npz").exists():
            raise SystemExit("--sc-shadow needs %s/sc.npz (cops_movies.py sc)" % work)
        SCZ = np.load(work / "sc.npz", allow_pickle=False)
    ucol = np.array([int(a.uniform_color[i:i + 2], 16) / 255.0 for i in (1, 3, 5)])
    fd = Path(a.frames_dir) if a.frames_dir else work / "frames"
    for sub in ("rg_clean", "ng_clean", "meta", "claims"):
        (fd / sub).mkdir(parents=True, exist_ok=True)
    order = [int(v) for v in a.only.split(",")] if a.only else list(range(len(P["utc"])))
    nw = max(1, a.nworkers)
    order = order[a.worker::nw] + [k for i, k in enumerate(order) if i % nw != a.worker]  # own stride, then help
    done, t_sum = 0, 0.0
    host = socket.gethostname()
    if is_fixed(Q):                            # 2026-11-13: the fixed camera (plan --camera fixed; new option)
        return render_fixed(a, P, Q, rnd, specs, FV, SCZ, tints, ucol, fd, order)
    if is_detach(Q):                           # 2026-11-13: the detach-200 m set (plan --camera detach; new option)
        return render_detach(a, P, Q, rnd, specs, FV, SCZ, tints, ucol, fd, order)
    for k in order:
        if (fd / "meta" / ("f%05d.json" % k)).exists():
            continue
        try:
            fdesc = os.open(fd / "claims" / ("f%05d" % k), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            continue
        os.write(fdesc, ("%s %s worker %d device %s\n" % (host, os.environ.get("SLURM_JOB_ID"), a.worker, a.device)).encode())
        os.close(fdesc)
        t0 = time.time()
        proj, F, center = str(Q["proj"][k]), float(Q["F"][k]), tuple(float(v) for v in Q["center"][k])
        vw = rnd.view(P["x"][k], P["b_rg"][k], P["b_ng"][k], P["sun"][k], proj, F, center)
        t1 = time.time()
        pw = rnd.pixel_weights(vw)
        mstats = {}
        sunfac, scst = None, None
        dirs = None
        if SCZ is not None:                    # spacecraft shadow (movies only, 2026-11-13)
            dirs = rnd.pixel_dirs_world(vw)
            sh = ScShadow(rnd.cp, SCZ, k, P["sun"][k], a.sc_samples)
            pidx, pts = rnd.pixel_points(vw, dirs, P["x"][k])
            bl, ncand = sh.blocked(pts)
            sunfac = rnd.cp.ones(vw["face"].size)
            sunfac[pidx] = 1.0 - bl
            sunfac = sunfac.reshape(vw["face"].shape)
            nb = rnd.cp.asnumpy(bl)
            scst = {"samples": int(a.sc_samples), "candidate_pixels_ss": ncand, "umbra_pixels_ss": int((nb >= 1.0).sum()),
                    "penumbra_pixels_ss": int(((nb > 0) & (nb < 1.0)).sum()), "sa_held": bool(SCZ["held"][k]),
                    "solar_angular_radius_deg": math.degrees(sh.alpha)}
            if (nb > 0).any():
                w_ = nb[nb > 0]
                scst["shadow_centroid_body_m"] = (rnd.cp.asnumpy(pts)[nb > 0] * w_[:, None]).sum(axis=0).tolist()
                scst["shadow_centroid_body_m"] = [v / w_.sum() for v in scst["shadow_centroid_body_m"]]
            del pts, bl
        if not masked:                         # the first movies (2026-11-12): unchanged path
            for which in ("rg", "ng"):
                Image.fromarray(rnd.colorize(vw, which, pw=pw, sunfac=sunfac)).save(
                    fd / ("%s_clean" % which) / ("f%05d.png" % k), compress_level=1)
        else:
            dirs = dirs if dirs is not None else rnd.pixel_dirs_world(vw)
            # facet directions of the drawn facets (visible-surface area fraction inside the mask)
            cp_ = rnd.cp
            fdraw = cp_.unique(vw["face"][vw["face"] >= 0])
            dc = rnd.cent[fdraw] - cp_.asarray(np.asarray(P["x"][k], float))[None, :]
            dc /= cp_.linalg.norm(dc, axis=1)[:, None]
            adraw = cp_.asarray(rnd.area_h)[fdraw]
            for which, g in (("rg", "RG"), ("ng", "NG")):
                bnd = FV["frame_bounds_%s" % which][k] if FV is not None else None
                if specs[g] == "none":
                    img = rnd.colorize(vw, which, pw=pw, sunfac=sunfac)
                    mstats[g] = {"mask": "none"}
                else:
                    mk = rnd.mask_dirs(dirs, g, specs[g], P["b_rg"][k], P["b_ng"][k], bnd)
                    img = rnd.colorize_masked(vw, which, pw, mk, tints[g], ucol, not a.no_boundary, sunfac=sunfac)
                    mf = rnd.mask_dirs((dc[:, 0], dc[:, 1], dc[:, 2]), g, specs[g], P["b_rg"][k], P["b_ng"][k], bnd)
                    st = rnd.mask_stats(vw, P["x"][k], mk, pw[0] if which == "rg" else pw[1])
                    st["visible_area_frac"] = float(cp_.asnumpy(adraw[mf].sum() / adraw.sum()))
                    st["visible_area_m2_in_mask"] = float(cp_.asnumpy(adraw[mf].sum()))
                    st["mask"], st["tint"] = specs[g], tints[g]
                    mstats[g] = st
                    del mk, mf
                Image.fromarray(img).save(fd / ("%s_clean" % which) / ("f%05d.png" % k), compress_level=1)
            del dirs, dc
        del pw
        cp = rnd.cp
        drawn = int((vw["face"] >= 0).sum())
        meta = {"k": k, "utc": str(P["utc"][k]), "et": float(P["et"][k]), "altitude_m": float(P["alt"][k]),
                "projection": proj, "F": F, "theta_max_deg": float(Q["theta_max_deg"][k]),
                "sub_point": [float(v) for v in Q["sub"][k]], "lit_frac_day": float(rnd.cp.asnumpy(
                    (vw["lit"].sum() / cp.maximum((vw["cos"] > 0).sum(), 1)))),
                "drawn_pixels": drawn, "host": host, "job": os.environ.get("SLURM_JOB_ID"), "worker": a.worker,
                "device": a.device, "t_view_s": t1 - t0, "t_total_s": time.time() - t0, "holes": vw["holes"]}
        if scst is not None:
            meta["sc_shadow"] = scst
        if masked:
            meta["masks"] = mstats
            meta["mask_design"] = {g: mask_design(g, specs[g], FV) for g in ("RG", "NG")}
        f, up, right = rnd.camera(P["x"][k], P["b_rg"][k], P["b_ng"][k])
        meta["center"] = list(center)
        meta["rg_axis_pixel"] = to_pixel((0.0, 0.0, 1.0), proj, F, center)
        bng = unit(P["b_ng"][k])
        meta["ng_axis_pixel"] = to_pixel((float(bng @ right), float(bng @ up), float(bng @ f)), proj, F, center)
        sp = np.asarray(Q["sub"][k]) - P["x"][k]
        if np.all(np.isfinite(sp)):
            zz, xx, yy = sp @ f, sp @ right, sp @ up
            px, py = to_pixel((xx, yy, zz), proj, F, center)
            dsub = float(np.linalg.norm(sp))
            th = math.atan2(math.hypot(xx, yy), zz)
            meta["sub_pixel"] = [px, py]
            meta["sub_distance_m"] = dsub
            meta["sub_theta_deg"] = math.degrees(th)
            # tangential scale at the sub-spacecraft point (perspective: d cos(theta) / F; equidistant fisheye:
            # d sin(theta) / (F theta), -> d / F on the axis)
            meta["m_per_px_at_sub"] = (dsub * math.cos(th) / F if proj == "perspective" else
                                       (dsub * math.sin(th) / (F * th) if th > 1e-9 else dsub / F))
        json.dump(meta, open(fd / "meta" / ("f%05d.json" % k), "w"))
        del vw
        done += 1
        t_sum += time.time() - t0
        if done % 25 == 0:
            say("worker %d (%s dev %s): %d frames, mean %.2f s/frame" % (a.worker, host, a.device, done, t_sum / done))
    say("worker %d (%s device %s) finished: %d frames, %.1f s, mean %.2f s/frame" % (
        a.worker, host, a.device, done, t_sum, t_sum / max(done, 1)))
    return 0


def fixed_check_masks(specs):
    for g, sp in specs.items():
        mode, ang = parse_mask(sp)
        if mode not in ("ik_fov", "cone") or (mode == "cone" and ang >= 80.0):
            raise SystemExit("fixed camera: --mask-%s must be ik_fov or cone:<deg < 80> (the gauge's view is a cone from "
                             "the moving gauge), got %s" % (g.lower(), sp))


def fixed_camera_meta(Q, vw0, axis_dist):
    f, up, right = vw0["basis"]
    return {"mode": "fixed", "utc": str(Q["cam_utc"]), "et": float(Q["cam_et"]),
            "altitude_target_m": float(Q["cam_altitude_target_m"]), "fov_deg": float(Q["cam_fov_deg"]),
            "x": [float(v) for v in Q["cam_x"]], "f": [float(v) for v in f], "up": [float(v) for v in up],
            "right": [float(v) for v in right], "axis_surface_distance_m": axis_dist}


def fixed_frame(rnd, P, FV, SCZ, specs, k, vw0, cen_d, fdraw, adraw, sc_samples, refine=True, sc_parts=("bus", "wings")):
    """per-frame products of the fixed camera: lighting, per gauge the weight, FOV mask, seen mask (FOV, facing,
    line of sight: depth map + exact refinement along the edges of the unseen parts), statistics (the area statistics
    of the camera-visible facets use the depth-map test at their centroids); spacecraft shadow factor."""
    cp = rnd.cp
    Hs, Ws = H * SS, W * SS
    idx, pts, fid, nrm = vw0["idx"], vw0["pts"], vw0["fid"], vw0["nrm"]
    cosz, lit = rnd.lit(P["sun"][k])
    vw = {"face": vw0["face"], "cos": cosz, "lit": lit}
    xg = np.asarray(P["x"][k], float)
    brg, bng = P["b_rg"][k], P["b_ng"][k]
    w_rg, w_ng, l = rnd.point_weights(xg, bng, brg, pts)
    sunfac, scst = None, None
    if SCZ is not None:
        sh = ScShadow(cp, SCZ, k, P["sun"][k], sc_samples)
        bl, ncand = sh.blocked(pts, sc_parts)
        sunfac = cp.ones(Hs * Ws)
        sunfac[idx] = 1.0 - bl
        sunfac = sunfac.reshape(Hs, Ws)
        nb = cp.asnumpy(bl)
        scst = {"samples": int(sc_samples), "candidate_pixels_ss": ncand, "umbra_pixels_ss": int((nb >= 1.0).sum()),
                "penumbra_pixels_ss": int(((nb > 0) & (nb < 1.0)).sum()), "sa_held": bool(SCZ["held"][k]),
                "solar_angular_radius_deg": math.degrees(sh.alpha)}
        if (nb > 0).any():
            w_ = nb[nb > 0]
            c_ = (cp.asnumpy(pts)[nb > 0] * w_[:, None]).sum(axis=0) / w_.sum()
            scst["shadow_centroid_body_m"] = c_.tolist()
        del bl
    facing = cp.sum(nrm * l, axis=1) > 0.0                    # n . (gauge - P) > 0 (l: unit vector P -> gauge)
    dvec = (-l[:, 0], -l[:, 1], -l[:, 2])                      # gauge -> P
    lc = cp.asarray(xg)[None, :] - cen_d
    lc /= cp.linalg.norm(lc, axis=1)[:, None]
    facing_c = cp.sum(rnd.norm[fdraw] * lc, axis=1) > 0.0
    dvec_c = (-lc[:, 0], -lc[:, 1], -lc[:, 2])
    allp = cp.concatenate([pts, cen_d])
    allf = cp.concatenate([fid, fdraw])
    out = {"vw": vw, "sunfac": sunfac, "sc": scst, "g": {}}
    npx = int(idx.size)
    for which, g in (("rg", "RG"), ("ng", "NG")):
        bnd = FV["frame_bounds_%s" % which][k] if FV is not None else None
        fov = rnd.mask_dirs(dvec, g, specs[g], brg, bng, bnd)
        fov_c = rnd.mask_dirs(dvec_c, g, specs[g], brg, bng, bnd)
        sel, sel_c = fov & facing, fov_c & facing_c
        vis, lst = rnd.gauge_los(xg, allp, allf, cp.concatenate([sel, sel_c]), brg if g == "RG" else bng)
        seen, seen_c = sel & vis[:npx], sel_c & vis[npx:]
        seen_map = seen
        if refine:                             # exact rays along the edges of the unseen parts
            seen, lst["refine"] = rnd.refine_los(xg, vw0, fov, seen)
        wv = w_rg if g == "RG" else w_ng
        st = {"mask": specs[g], "line_of_sight": True, "nucleus_pixel_frac": float(cp.asnumpy(seen.sum())) / max(npx, 1),
              "fov_pixel_frac": float(cp.asnumpy(fov.sum())) / max(npx, 1),
              "fov_facing_pixel_frac": float(cp.asnumpy(sel.sum())) / max(npx, 1),
              "visible_area_frac": float(cp.asnumpy(adraw[seen_c].sum() / adraw.sum())),
              "visible_area_m2_in_mask": float(cp.asnumpy(adraw[seen_c].sum())),
              "fov_area_m2": float(cp.asnumpy(adraw[fov_c].sum())), "fov_facing_area_m2": float(cp.asnumpy(adraw[sel_c].sum())),
              "los": lst}
        if bool(seen.any()):
            st["weight_min_on_nucleus_in_mask"] = float(cp.asnumpy(wv[seen].min()))
            st["weight_max_on_nucleus_in_mask"] = float(cp.asnumpy(wv[seen].max()))
        out["g"][g] = {"fov": fov, "seen": seen, "seen_map": seen_map, "w": wv, "stats": st}
        del vis, fov_c, sel_c, seen_c
    return out


def _full(rnd, idx, vals, dtype):
    cp = rnd.cp
    a = cp.zeros(H * SS * W * SS, dtype=dtype)
    a[idx] = vals
    return a.reshape(H * SS, W * SS)


def render_fixed(a, P, Q, rnd, specs, FV, SCZ, tints, ucol, fd, order):
    """render of the fixed-camera movies (2026-11-13): the camera view once (fixed_geometry), then per frame the
    lighting of the frame time and per gauge the surface it sees (fixed_frame); clean frames rg_clean / ng_clean and
    the frame metadata (camera, face-buffer sha256, statistics, sub-spacecraft point and gauge in the view)."""
    from PIL import Image
    cp = rnd.cp
    fixed_check_masks(specs)
    xc = np.asarray(Q["cam_x"], float)
    proj, F, center = "perspective", float(Q["F"][0]), (0.0, 0.0)
    t0 = time.time()
    vw0 = rnd.fixed_geometry(xc, Q["cam_b_rg"], Q["cam_b_ng"], proj, F, center)
    f_, up_, right_ = vw0["basis"]
    fdraw = cp.unique(vw0["fid"])
    adraw = cp.asarray(rnd.area_h)[fdraw]
    cen_d = rnd.cent[fdraw]
    hf, th = rnd.first_hits(xc[None, :], np.asarray(f_)[None, :])
    axis_dist = float(th[0]) if hf[0] >= 0 else None
    camd = fixed_camera_meta(Q, vw0, axis_dist)
    say("fixed camera view: %d nucleus pixels (supersampled), %d facets, holes %s, face buffer sha256 %s (%.1f s)" % (
        int(vw0["idx"].size), int(fdraw.size), json.dumps(vw0["holes"]), vw0["face_sha256"][:16], time.time() - t0))
    host = socket.gethostname()
    done, t_sum = 0, 0.0
    for k in order:
        if (fd / "meta" / ("f%05d.json" % k)).exists():
            continue
        try:
            fdesc = os.open(fd / "claims" / ("f%05d" % k), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            continue
        os.write(fdesc, ("%s %s worker %d device %s\n" % (host, os.environ.get("SLURM_JOB_ID"), a.worker, a.device)).encode())
        os.close(fdesc)
        t0 = time.time()
        R = fixed_frame(rnd, P, FV, SCZ, specs, k, vw0, cen_d, fdraw, adraw, a.sc_samples)
        t1 = time.time()
        mstats = {}
        for which, g in (("rg", "RG"), ("ng", "NG")):
            G = R["g"][g]
            img = rnd.colorize_fixed(R["vw"], which, _full(rnd, vw0["idx"], G["w"], cp.float64),
                                     _full(rnd, vw0["idx"], G["seen"], bool), _full(rnd, vw0["idx"], G["fov"], bool),
                                     tints[g], ucol, not a.no_boundary, sunfac=R["sunfac"])
            Image.fromarray(img).save(fd / ("%s_clean" % which) / ("f%05d.png" % k), compress_level=1)
            st = G["stats"]
            st["tint"] = tints[g]
            mstats[g] = st
        vw = R["vw"]
        meta = {"k": k, "utc": str(P["utc"][k]), "et": float(P["et"][k]), "altitude_m": float(P["alt"][k]),
                "projection": proj, "F": F, "theta_max_deg": float(Q["theta_max_deg"][k]),
                "sub_point": [float(v) for v in Q["sub"][k]], "lit_frac_day": float(cp.asnumpy(
                    (vw["lit"].sum() / cp.maximum((vw["cos"] > 0).sum(), 1)))),
                "drawn_pixels": int(vw0["idx"].size), "host": host, "job": os.environ.get("SLURM_JOB_ID"),
                "worker": a.worker, "device": a.device, "t_view_s": t1 - t0, "holes": vw0["holes"],
                "camera": camd, "face_sha256": vw0["face_sha256"], "masks": mstats,
                "mask_design": {g: dict(mask_design(g, specs[g], FV), line_of_sight=True) for g in ("RG", "NG")}}
        if R["sc"] is not None:
            meta["sc_shadow"] = R["sc"]
        meta["center"] = list(center)
        if axis_dist is not None:
            meta["m_per_px_at_center"] = axis_dist / F
        # the sub-spacecraft point (+) and the gauge in the fixed view (exact visibility from the camera)
        sp = np.asarray(Q["sub"][k], float)
        xg = np.asarray(P["x"][k], float)
        pts = [sp if np.all(np.isfinite(sp)) else xc + f_, xg]
        vis, _, _, dist = rnd.points_visible(np.vstack([xc, xc]), np.vstack(pts), [int(Q["sub_face"][k]), -1])
        for name, p_, ok_ in (("sub", pts[0], np.all(np.isfinite(sp))), ("gauge", pts[1], True)):
            v_ = p_ - xc
            pix = to_pixel((float(v_ @ right_), float(v_ @ up_), float(v_ @ f_)), proj, F, center) if ok_ else None
            meta["%s_pixel" % name] = pix
            meta["%s_visible" % name] = bool(ok_ and pix is not None and vis[0 if name == "sub" else 1])
            meta["%s_distance_from_camera_m" % name] = float(dist[0 if name == "sub" else 1])
        meta["t_total_s"] = time.time() - t0
        json.dump(meta, open(fd / "meta" / ("f%05d.json" % k), "w"))
        del R, vw
        done += 1
        t_sum += time.time() - t0
        if done % 25 == 0:
            say("worker %d (%s dev %s): %d frames, mean %.2f s/frame" % (a.worker, host, a.device, done, t_sum / done))
    say("worker %d (%s device %s) finished: %d frames, %.1f s, mean %.2f s/frame (fixed camera)" % (
        a.worker, host, a.device, done, t_sum, t_sum / max(done, 1)))
    return 0


def cone_segments(rnd, x_g, axis, half_deg, rgb, alt):
    """generator lines of a FOV cone (CONE_N, translucent) from the gauge point x_g around axis: each to its first hit on
    the nucleus (exact ray, all facets); generators without a hit to the longest hit length (2 x the COPS altitude if
    none hits)."""
    b = unit(axis)
    a0 = np.array([1.0, 0, 0]) if abs(b[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = unit(np.cross(b, a0))
    e2 = np.cross(b, e1)
    ang = np.radians(half_deg)
    phi = np.linspace(0, 2 * np.pi, CONE_N, endpoint=False)
    dirs = np.cos(ang) * b[None, :] + np.sin(ang) * (np.cos(phi)[:, None] * e1[None, :] + np.sin(phi)[:, None] * e2[None, :])
    hf, th = rnd.first_hits(np.tile(np.asarray(x_g, float), (CONE_N, 1)), dirs)
    hit = hf >= 0
    Lmax = float(th[hit].max()) if hit.any() else 2.0 * max(float(alt), 1.0)
    L = np.where(hit, th, Lmax)
    segs = [(np.asarray(x_g, float), np.asarray(x_g, float) + L[i] * dirs[i], rgb, CONE_ALPHA) for i in range(CONE_N)]
    return segs, {"generators": CONE_N, "hitting": int(hit.sum()), "length_m_min_max": [float(L.min()), float(L.max())]}


def render_detach(a, P, Q, rnd, specs, FV, SCZ, tints, ucol, fd, order):
    """render of the detach-200 m movies (2026-11-13): frames before the detach frame as the gauge-view movies (camera
    at COPS along the RG axis; tint by the pixel direction inside the gauge's cone; spacecraft shadow incl. the HGA);
    from the detach frame on the frozen perspective camera (fixed_geometry once): per frame the surface each gauge sees
    from its true position (cone, facing, line of sight, as the fixed camera), the spacecraft model (sc_model: bus, wings,
    HGA; Lambert, Sun of the frame time; near clipping SC_NEAR), its shadow (bus, wings, HGA), the depth-tested
    translucent FOV cone of the movie's gauge and the yoke / strut lines."""
    from PIL import Image
    cp = rnd.cp
    if SCZ is None or "hga_c_sc" not in SCZ.files:
        raise SystemExit("detach render needs --sc-shadow and an sc.npz with the HGA (cops_movies.py sc --hga-ck)")
    for g, sp in specs.items():
        if parse_mask(sp)[0] != "cone":
            raise SystemExit("detach render: --mask-%s must be cone:<deg> (user decision: circular FOVs)" % g.lower())
    fixed_check_masks(specs)
    parts = ("bus", "wings", "hga")
    work = Path(a.work)
    td = load_touchdown(work / "touchdown.npz") if (work / "touchdown.npz").exists() else None
    k_d = int(Q["k_detach"])
    n_obs = int(P["n_obs_frames"]) if "n_obs_frames" in P else len(P["utc"])
    xc = np.asarray(Q["cam_x_frames"][k_d], float)
    Fd = float(Q["F"][k_d])
    cbr_c, cbn_c = Q["cam_b_rg_frames"][k_d], Q["cam_b_ng_frames"][k_d]
    pan = "pan_s" in Q.files                   # 2026-11-14: the moving camera after the detach (plan --pan-to)
    vkey = None
    camd_base = {"mode": "detach", "detach_frame": k_d, "detach_utc": str(Q["detach_utc"]),
                 "detach_altitude_m": float(Q["detach_altitude_m"]), "fov_deg": float(Q["cam_fov_deg"])}
    host = socket.gethostname()
    vw0 = None
    done, t_sum = 0, 0.0
    for k in order:
        if (fd / "meta" / ("f%05d.json" % k)).exists():
            continue
        try:
            fdesc = os.open(fd / "claims" / ("f%05d" % k), os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            continue
        os.write(fdesc, ("%s %s worker %d device %s\n" % (host, os.environ.get("SLURM_JOB_ID"), a.worker, a.device)).encode())
        os.close(fdesc)
        t0 = time.time()
        if k < k_d:
            meta = _detach_attached_frame(a, P, Q, rnd, specs, FV, SCZ, tints, ucol, fd, k, parts, host, t0)
            meta["camera"] = dict(camd_base, attached=True)
        else:
            if pan:                            # the frame's own camera pose; the view is recomputed when it changes
                xk, brk, bnk = frame_camera(P, Q, k)
                key = (xk.tobytes(), np.asarray(brk).tobytes(), np.asarray(bnk).tobytes(), float(Q["F"][k]))
                if key != vkey:
                    vw0 = None
                    xc, cbr_c, cbn_c, Fd, vkey = xk, brk, bnk, float(Q["F"][k]), key
            if vw0 is None:
                vw0 = rnd.fixed_geometry(xc, cbr_c, cbn_c, "perspective", Fd, (0.0, 0.0))
                fdraw = cp.unique(vw0["fid"])
                adraw = cp.asarray(rnd.area_h)[fdraw]
                cen_d = rnd.cent[fdraw]
                hf0, th0 = rnd.first_hits(xc[None, :], np.asarray(vw0["basis"][0])[None, :])
                axis_dist = float(th0[0]) if hf0[0] >= 0 else None
                nuc_depth = cp.full(H * SS * W * SS, cp.inf)
                nuc_depth[vw0["idx"]] = cp.sqrt(cp.sum((vw0["pts"] - cp.asarray(xc)[None, :]) ** 2, axis=1))
                nuc_depth = nuc_depth.reshape(H * SS, W * SS)
                say("detached camera view: %d nucleus pixels, %d facets, face buffer %s" % (
                    int(vw0["idx"].size), int(fdraw.size), vw0["face_sha256"][:16]))
            f_, up_, right_ = vw0["basis"]
            R = fixed_frame(rnd, P, FV, SCZ, specs, k, vw0, cen_d, fdraw, adraw, a.sc_samples, sc_parts=parts)
            t1 = time.time()
            # spacecraft model
            M = sc_model_body(SCZ, k)
            cen_t = M["V"][M["T"]].mean(axis=1)
            facing = np.einsum("ij,ij->i", M["normal"], xc[None, :] - cen_t) > 0
            draw = M["two_sided"] | facing
            tri, dep, ndrawn = rnd.raster_model(xc, vw0["basis"], Fd, (0.0, 0.0), M["V"], M["T"], draw)
            s_b = unit(P["sun"][k])
            nv = np.where(facing[:, None], M["normal"], -M["normal"])
            colt = np.where(facing[:, None], M["col_front"], M["col_back"])
            bt = AMBIENT + DIFFUSE * np.maximum(0.0, nv @ s_b)
            rgb_t = cp.asarray(colt * bt[:, None])
            sc_front = (tri >= 0) & (dep < nuc_depth)
            sc_img = rgb_t[cp.clip(tri, 0, None)]
            depth = cp.where(sc_front, dep, nuc_depth)
            segs = [(p0, p1, c, 1.0) for p0, p1, c in M["lines"]]
            xg = np.asarray(P["x"][k], float)
            cone_st = {}
            mstats = {}
            for which, g in (("rg", "RG"), ("ng", "NG")):
                cs, cone_st[g] = cone_segments(rnd, xg, P["b_rg"][k] if g == "RG" else P["b_ng"][k],
                                               parse_mask(specs[g])[1], CONE_COL[g], P["alt"][k])
                lines = rnd.paint_lines(segs + cs, xc, vw0["basis"], Fd, (0.0, 0.0), depth)
                G = R["g"][g]
                img = rnd.colorize_detach(R["vw"], which, _full(rnd, vw0["idx"], G["w"], cp.float64),
                                          _full(rnd, vw0["idx"], G["seen"], bool), _full(rnd, vw0["idx"], G["fov"], bool),
                                          sc_front, sc_img, lines, tints[g], ucol, not a.no_boundary, sunfac=R["sunfac"])
                Image.fromarray(img).save(fd / ("%s_clean" % which) / ("f%05d.png" % k), compress_level=1)
                st = G["stats"]
                st["tint"] = tints[g]
                mstats[g] = st
                del lines
            vw = R["vw"]
            meta = {"k": k, "utc": str(P["utc"][k]), "et": float(P["et"][k]), "altitude_m": float(P["alt"][k]),
                    "projection": "perspective", "F": Fd, "theta_max_deg": float(Q["theta_max_deg"][k]),
                    "sub_point": [float(v) for v in Q["sub"][k]], "lit_frac_day": float(cp.asnumpy(
                        (vw["lit"].sum() / cp.maximum((vw["cos"] > 0).sum(), 1)))),
                    "drawn_pixels": int(vw0["idx"].size), "host": host, "job": os.environ.get("SLURM_JOB_ID"),
                    "worker": a.worker, "device": a.device, "t_view_s": t1 - t0, "holes": vw0["holes"],
                    "face_sha256": vw0["face_sha256"], "masks": mstats,
                    "mask_design": {g: dict(mask_design(g, specs[g], FV), line_of_sight=True) for g in ("RG", "NG")},
                    "camera": dict(camd_base, attached=False, x=xc.tolist(), f=[float(v) for v in f_],
                                   up=[float(v) for v in up_], right=[float(v) for v in right_],
                                   axis_surface_distance_m=axis_dist, **({} if not pan else {
                                       "pan": True, "pan_s": float(Q["pan_s"][k]),
                                       "fov_deg": 2 * math.degrees(math.atan((W / 2.0) / Fd))})),
                    "sc_model": {"triangles": int(len(M["T"])), "drawn": ndrawn,
                                 "pixels_ss": int(cp.asnumpy(sc_front.sum())), "hga_held": bool(SCZ["hga_held"][k]),
                                 "sa_held": bool(SCZ["held"][k])},
                    "cones": cone_st, "center": [0.0, 0.0]}
            if R["sc"] is not None:
                meta["sc_shadow"] = R["sc"]
            if axis_dist is not None:
                meta["m_per_px_at_center"] = axis_dist / Fd
            sp = np.asarray(Q["sub"][k], float)
            okp = bool(np.all(np.isfinite(sp)))
            pts_ = [sp if okp else xc + f_, xg]
            vis, _, _, dist = rnd.points_visible(np.vstack([xc, xc]), np.vstack(pts_), [int(Q["sub_face"][k]), -1])
            for name, p_, ok_ in (("sub", pts_[0], okp), ("gauge", pts_[1], True)):
                v_ = p_ - xc
                pix = to_pixel((float(v_ @ right_), float(v_ @ up_), float(v_ @ f_)), "perspective", Fd, (0.0, 0.0)) if ok_ else None
                meta["%s_pixel" % name] = pix
                vis_ = bool(ok_ and pix is not None and vis[0 if name == "sub" else 1])
                if name == "sub" and vis_ and 0 <= pix[0] < W and 0 <= pix[1] < H:   # hidden by the spacecraft?
                    ii, jj = int(pix[0] * SS), int(pix[1] * SS)
                    vis_ = not bool(cp.asnumpy(sc_front[jj, ii]))
                meta["%s_visible" % name] = vis_ if name == "sub" else bool(pix is not None)   # the gauge point: always marked
                meta["%s_distance_from_camera_m" % name] = float(dist[0 if name == "sub" else 1])
            if td is not None and bool(td["visible"][k]):
                bx, by = td["base_px"][k]
                if 0 <= bx < W and 0 <= by < H:
                    ii, jj = int(bx * SS), int(by * SS)
                    meta["sc_hides_touchdown"] = bool(cp.asnumpy(sc_front[jj, ii] & (dep[jj, ii] < float(td["distance_m"][k]))))
            del R, vw, tri, dep, sc_front, sc_img, depth
        meta["after_last_obs"] = bool(k >= n_obs)
        meta["t_total_s"] = time.time() - t0
        json.dump(meta, open(fd / "meta" / ("f%05d.json" % k), "w"))
        done += 1
        t_sum += time.time() - t0
        if done % 25 == 0:
            say("worker %d (%s dev %s): %d frames, mean %.2f s/frame" % (a.worker, host, a.device, done, t_sum / done))
    say("worker %d (%s device %s) finished: %d frames, %.1f s, mean %.2f s/frame (detach)" % (
        a.worker, host, a.device, done, t_sum, t_sum / max(done, 1)))
    return 0


def _detach_attached_frame(a, P, Q, rnd, specs, FV, SCZ, tints, ucol, fd, k, parts, host, t0):
    """an attached frame of the detach-200 m movies: the gauge-view render of cmd_render (masked path) with the
    spacecraft shadow incl. the HGA (parts); returns the frame metadata."""
    from PIL import Image
    cp = rnd.cp
    proj, F, center = str(Q["proj"][k]), float(Q["F"][k]), tuple(float(v) for v in Q["center"][k])
    vw = rnd.view(P["x"][k], P["b_rg"][k], P["b_ng"][k], P["sun"][k], proj, F, center)
    t1 = time.time()
    pw = rnd.pixel_weights(vw)
    mstats = {}
    dirs = rnd.pixel_dirs_world(vw)
    sh = ScShadow(cp, SCZ, k, P["sun"][k], a.sc_samples)
    pidx, pts = rnd.pixel_points(vw, dirs, P["x"][k])
    bl, ncand = sh.blocked(pts, parts)
    sunfac = cp.ones(vw["face"].size)
    sunfac[pidx] = 1.0 - bl
    sunfac = sunfac.reshape(vw["face"].shape)
    nb = cp.asnumpy(bl)
    scst = {"samples": int(a.sc_samples), "parts": list(parts), "candidate_pixels_ss": ncand,
            "umbra_pixels_ss": int((nb >= 1.0).sum()), "penumbra_pixels_ss": int(((nb > 0) & (nb < 1.0)).sum()),
            "sa_held": bool(SCZ["held"][k]), "hga_held": bool(SCZ["hga_held"][k]),
            "solar_angular_radius_deg": math.degrees(sh.alpha)}
    del pts, bl
    fdraw = cp.unique(vw["face"][vw["face"] >= 0])
    dc = rnd.cent[fdraw] - cp.asarray(np.asarray(P["x"][k], float))[None, :]
    dc /= cp.linalg.norm(dc, axis=1)[:, None]
    adraw = cp.asarray(rnd.area_h)[fdraw]
    for which, g in (("rg", "RG"), ("ng", "NG")):
        bnd = FV["frame_bounds_%s" % which][k] if FV is not None else None
        mk = rnd.mask_dirs(dirs, g, specs[g], P["b_rg"][k], P["b_ng"][k], bnd)
        img = rnd.colorize_masked(vw, which, pw, mk, tints[g], ucol, not a.no_boundary, sunfac=sunfac)
        mf = rnd.mask_dirs((dc[:, 0], dc[:, 1], dc[:, 2]), g, specs[g], P["b_rg"][k], P["b_ng"][k], bnd)
        st = rnd.mask_stats(vw, P["x"][k], mk, pw[0] if which == "rg" else pw[1])
        st["visible_area_frac"] = float(cp.asnumpy(adraw[mf].sum() / adraw.sum()))
        st["visible_area_m2_in_mask"] = float(cp.asnumpy(adraw[mf].sum()))
        st["mask"], st["tint"] = specs[g], tints[g]
        mstats[g] = st
        Image.fromarray(img).save(fd / ("%s_clean" % which) / ("f%05d.png" % k), compress_level=1)
        del mk, mf
    meta = {"k": k, "utc": str(P["utc"][k]), "et": float(P["et"][k]), "altitude_m": float(P["alt"][k]),
            "projection": proj, "F": F, "theta_max_deg": float(Q["theta_max_deg"][k]),
            "sub_point": [float(v) for v in Q["sub"][k]], "lit_frac_day": float(cp.asnumpy(
                (vw["lit"].sum() / cp.maximum((vw["cos"] > 0).sum(), 1)))),
            "drawn_pixels": int((vw["face"] >= 0).sum()), "host": host, "job": os.environ.get("SLURM_JOB_ID"),
            "worker": a.worker, "device": a.device, "t_view_s": t1 - t0, "holes": vw["holes"], "sc_shadow": scst,
            "masks": mstats, "mask_design": {g: mask_design(g, specs[g], FV) for g in ("RG", "NG")}}
    f, up, right = rnd.camera(P["x"][k], P["b_rg"][k], P["b_ng"][k])
    meta["center"] = list(center)
    meta["rg_axis_pixel"] = to_pixel((0.0, 0.0, 1.0), proj, F, center)
    bng = unit(P["b_ng"][k])
    meta["ng_axis_pixel"] = to_pixel((float(bng @ right), float(bng @ up), float(bng @ f)), proj, F, center)
    sp = np.asarray(Q["sub"][k]) - P["x"][k]
    if np.all(np.isfinite(sp)):
        zz, xx, yy = sp @ f, sp @ right, sp @ up
        px, py = to_pixel((xx, yy, zz), proj, F, center)
        dsub = float(np.linalg.norm(sp))
        th = math.atan2(math.hypot(xx, yy), zz)
        meta["sub_pixel"] = [px, py]
        meta["sub_distance_m"] = dsub
        meta["sub_theta_deg"] = math.degrees(th)
        meta["m_per_px_at_sub"] = (dsub * math.cos(th) / F if proj == "perspective" else
                                   (dsub * math.sin(th) / (F * th) if th > 1e-9 else dsub / F))
    del vw, pw, dirs, dc
    return meta


# ============================================================================================ compose
def font(size):
    import matplotlib
    from PIL import ImageFont
    return ImageFont.truetype(str(Path(matplotlib.get_data_path()) / "fonts" / "ttf" / "DejaVuSans.ttf"), size)


def colorbar_img(which, height, width=26):
    import matplotlib
    lut = matplotlib.colormaps[CMAP](np.linspace(1, 0, height))[:, :3]
    arr = (np.repeat(lut[:, None, :], width, axis=1) * 255).astype(np.uint8)
    return arr


def nice_length(target):
    p = 10 ** math.floor(math.log10(target))
    for m in (5, 2, 1):
        if m * p <= target:
            return m * p
    return p


def fixed_texts(meta, which, gname, md):
    """title, subtitle and short title of the fixed-camera frames (2026-11-13)."""
    title = ("RG: where RG sees the surface (tint), seen from a fixed camera" if which == "rg" else
             "NG: where NG sees the surface (tint), seen from the same fixed camera")
    what = "RG acceptance weight -l·b_RG" if which == "rg" else "NG acceptance weight 1 - l·b_NG"
    wr = md.get("weight_range")
    lab = md["label"].replace("the %s IK field of view, " % gname, "")
    tint = meta["masks"][gname].get("tint", "weight")
    sub1 = ("tint: %s%s where the surface is inside the %s IK FOV %s (white line), faces %s and is in its line of "
            "sight; sunlit relief elsewhere" % (what if tint == "weight" else "one colour",
                                                (" (%.2f–%.3g)" % tuple(wr)) if (wr and tint == "weight") else "",
                                                gname, lab, gname))
    short = "%s (tint: %s weight where %s sees the surface: IK FOV, facing, line of sight)" % (gname, gname, gname)
    return title, sub1, short


def detach_texts(meta, which, gname, md):
    """title, subtitle and short title of the detached frames of the detach-200 m movies (2026-11-13)."""
    cam = meta["camera"]
    title = ("RG: where RG sees the surface (tint); camera detached at %s" % fmt_alt(cam["detach_altitude_m"])
             if which == "rg" else "NG: where NG sees the surface (tint); the same detached camera")
    what = "RG acceptance weight -l·b_RG" if which == "rg" else "NG acceptance weight 1 - l·b_NG"
    wr = md.get("weight_range")
    ang = md.get("max_angle_deg")
    sub1 = ("tint: %s%s where the surface is inside its %g° cone (white line on the surface, translucent lines "
            "from COPS), faces %s and is in its line of sight" % (
                what, (" (%.2f–%.3g)" % tuple(wr)) if wr else "", ang, gname))
    short = "%s (tint: %s weight in its %g° cone, facing, line of sight)" % (gname, gname, ang)
    return title, sub1, short


def flag_layer(k, t, up, static=False):
    """RGBA image of the chequered finish flag on its pole (2026-11-13), drawn FLAG_SS times larger and reduced
    (anti-aliasing); returns (image, (x, y) of the pole base in it). up: unit vector of the pole in image coordinates
    (x right, y down); the cloth hangs from the pole top to the right of the pole (seen along up); travelling wave of
    phase 2 pi (s / FLAG_WAVELENGTH - FLAG_WAVE_HZ k / FPS) (s = 0 .. 1 from the pole), amplitude FLAG_AMP (flag heights)
    times s; static: no wave (legend icon)."""
    from PIL import Image, ImageDraw
    ux, uy = float(up[0]), float(up[1])
    rx, ry = -uy, ux
    L, fw, fh = FLAG_POLE * t, FLAG_W * t, FLAG_H * t
    tsec = k / float(FPS)

    def phase(s):
        return 2.0 * math.pi * (s / FLAG_WAVELENGTH - FLAG_WAVE_HZ * tsec)

    def pt(s, v):
        ph = phase(s)
        dsp = 0.0 if static else FLAG_AMP * fh * s * math.sin(ph)
        along = s * fw * (1.0 if static else 1.0 - 0.05 * s * (1.0 - math.cos(ph)))
        down = v * fh - dsp
        return (ux * L + rx * along - ux * down, uy * L + ry * along - uy * down)

    nsub = 4
    polys = []
    for i in range(FLAG_NX):
        for j in range(FLAG_NY):
            base = 238.0 if (i + j) % 2 == 0 else 22.0
            for m in range(nsub):
                s0, s1 = (i + m / nsub) / FLAG_NX, (i + (m + 1) / nsub) / FLAG_NX
                v0, v1 = j / FLAG_NY, (j + 1) / FLAG_NY
                shade = 1.0 if static else 0.86 + 0.14 * math.cos(phase(0.5 * (s0 + s1)) + 0.6)
                c = int(round(min(255.0, base * shade if base > 100 else base + 18.0 * (1.0 - shade))))
                polys.append(([pt(s0, v0), pt(s1, v0), pt(s1, v1), pt(s0, v1)], (c, c, c, 255)))
    ns = 48
    outline = ([pt(s / ns, 0.0) for s in range(ns + 1)] + [pt(1.0, v / 8.0) for v in range(1, 9)] +
               [pt(1.0 - s / ns, 1.0) for s in range(1, ns + 1)] + [pt(0.0, 1.0 - v / 8.0) for v in range(1, 9)])
    top = (ux * (L + 4.0 * t), uy * (L + 4.0 * t))
    allp = [p for q, _ in polys for p in q] + [(0.0, 0.0), top]
    mg = 6.0 * t
    x0 = math.floor(min(p[0] for p in allp) - mg)
    y0 = math.floor(min(p[1] for p in allp) - mg)
    x1 = math.ceil(max(p[0] for p in allp) + mg)
    y1 = math.ceil(max(p[1] for p in allp) + mg)
    S = FLAG_SS
    im = Image.new("RGBA", ((x1 - x0) * S, (y1 - y0) * S), (0, 0, 0, 0))
    d = ImageDraw.Draw(im)

    def tr(p):
        return ((p[0] - x0) * S, (p[1] - y0) * S)

    d.line([tr((0.0, 0.0)), tr(top)], fill=(15, 15, 15, 255), width=max(1, int(round(4.4 * t * S))))
    for q, col in polys:
        d.polygon([tr(p) for p in q], fill=col)
    d.line([tr(p) for p in outline] + [tr(outline[0])], fill=(15, 15, 15, 255), width=max(1, int(round(1.7 * t * S))),
           joint="curve")
    d.line([tr((0.0, 0.0)), tr(top)], fill=(215, 215, 215, 255), width=max(1, int(round(2.0 * t * S))))
    r_b = 3.4 * t * S
    bx, by = tr((0.0, 0.0))
    d.ellipse((bx - r_b, by - r_b, bx + r_b, by + r_b), fill=(15, 15, 15, 255))
    r_b2 = 1.5 * t * S
    d.ellipse((bx - r_b2, by - r_b2, bx + r_b2, by + r_b2), fill=(245, 245, 245, 255))
    tx, ty = tr(top)
    r_t = 2.6 * t * S
    d.ellipse((tx - r_t, ty - r_t, tx + r_t, ty + r_t), fill=(15, 15, 15, 255))
    im = im.convert("RGBa").resize((x1 - x0, y1 - y0), Image.LANCZOS).convert("RGBA")
    return im, (-x0, -y0)


def draw_flag(img, base_xy, up, k, t):
    """the waving flag with its pole base at base_xy (output pixels of img, sub-pixel)."""
    lay, (bx, by) = flag_layer(k, t, up)
    # integer placement; the base point sits within half a pixel of base_xy (sub-pixel offset applied by a shift of
    # the drawing would need a re-draw: the pole base marker covers > 3 px)
    ox, oy = int(round(base_xy[0] - bx)), int(round(base_xy[1] - by))
    img.paste(lay.convert("RGB"), (ox, oy), lay.split()[3])
    return img


def load_touchdown(path):
    z = np.load(path, allow_pickle=False)
    return {key: z[key] for key in z.files}


def td_show(td, k, meta):
    """draw the flag in frame k: depth-tested (visible from the camera, not behind the spacecraft model), or with
    compose --flag-on-top (2026-11-13) whenever its point projects into the image."""
    if td.get("_on_top"):
        return bool(td["in_front"][k]) and bool(td["in_image"][k])
    return bool(td["visible"][k]) and not meta.get("sc_hides_touchdown")


def td_legend(img, d, td, k, t, pad, full=True, sc_hidden=False):
    """legend of the touchdown flag: a small static flag icon and "touchdown (impact) <UTC> UTC" (right-aligned, above
    the scale-bar line)."""
    Wi, Hi = img.size
    fn = font(int(round(20 * t)))
    vis = (bool(td["visible"][k]) and not sc_hidden) or bool(td.get("_on_top"))
    text = "touchdown (impact) %s UTC%s" % (str(td["utc"])[11:], "" if vis else (
        ": not visible from the camera" if full else ": not visible"))
    icon, (ibx, iby) = flag_layer(0, 0.42 * t, (0.0, -1.0), static=True)
    tw = d.textlength(text, font=fn)
    y_base = Hi - pad - int(round(70 * t))
    x_text = Wi - pad - tw
    x_icon = int(round(x_text - 8 * t - icon.size[0]))
    y_icon = int(round(y_base - icon.size[1] + 4 * t))
    d.rectangle((x_icon - 4 * t, y_icon - 3 * t, Wi - pad + 4 * t, y_base + 8 * t), fill=(0, 0, 0))
    img.paste(icon.convert("RGB"), (x_icon, y_icon), icon.split()[3])
    d.text((x_text, y_base), text, font=fn, fill=(255, 255, 255), anchor="ls")


def overlay(img, meta, which, scale=1.0, full=True, text_scale=None, td=None, k=None):
    """Text, colour bar, scale bar, sub-spacecraft mark on a view frame (PIL image, modified in place). scale: image size
    relative to the rendered W x H (geometry); text_scale: size of the text and the colour bar (default scale).
    td, k (2026-11-13, compose --touchdown): the touchdown legend of frame k."""
    from PIL import Image, ImageDraw
    d = ImageDraw.Draw(img)
    Wi, Hi = img.size
    s = scale
    t = scale if text_scale is None else text_scale
    f1, f2 = font(int(round(30 * t))), font(int(round(21 * t)))
    pad = int(14 * t)
    lw = max(1, int(round(2 * t)))
    gname = "RG" if which == "rg" else "NG"
    title = ("RG view: camera at COPS along the RG axis" if which == "rg" else
             "NG view: the same camera, tinted by the NG acceptance")
    sub1 = ("tint: RG acceptance weight -l·b_RG (flux into the RG aperture)" if which == "rg" else
            "tint: NG acceptance weight, 1 - l·b_NG in front of NG, shield fraction behind")
    md = meta.get("mask_design", {}).get(gname)
    if md and md["spec"] != "none":           # observable-part mask (2026-11-13)
        wr = md.get("weight_range")
        what = ("RG acceptance weight -l·b_RG" if which == "rg" else "NG acceptance weight 1 - l·b_NG")
        tint = meta["masks"][gname].get("tint", "weight")
        sub1 = ("tint inside %s (outlined): %s%s; grey outside" % (
            (md["label"].replace("field of view (", "field of view, ")[:-1] if "field of view (" in md["label"]
             else md["label"]), what if tint == "weight" else "one colour",
            (" (%.2f–%.3g there)" % tuple(wr)) if (wr and tint == "weight") else ""))
    when = "%s UTC   COPS altitude %s" % (meta["utc"][:19].replace("T", " "), fmt_alt(meta["altitude_m"]))
    short = ("%s view (tint: %s acceptance weight)" % (gname, gname) if not (md and md["spec"] != "none") else
             "%s view (tint: %s weight inside the %s %s, outlined)" % (
                 gname, gname, gname, "IK FOV" if md["spec"] == "ik_fov" else md["spec"]))
    cam = meta.get("camera")
    fixed = cam is not None and cam.get("mode") == "fixed"
    if fixed:                                  # 2026-11-13: the fixed camera
        title, sub1, short = fixed_texts(meta, which, gname, md)
    detached = cam is not None and cam.get("mode") == "detach" and not cam.get("attached")
    attached_d = cam is not None and cam.get("mode") == "detach" and bool(cam.get("attached"))
    if detached:                               # 2026-11-13: the detach-200 m set
        title, sub1, short = detach_texts(meta, which, gname, md)
    if attached_d and md and parse_mask(md["spec"])[0] == "cone":
        short = "%s view (tint: %s weight inside its %g° cone, outlined)" % (gname, gname, parse_mask(md["spec"])[1])
    if meta.get("after_last_obs"):
        when = "%s UTC   COPS altitude %s (after the last COPS observation; SPICE trajectory)" % (
            meta["utc"][:21].replace("T", " "), fmt_alt(meta["altitude_m"]))
    lines = [(title, f1), (sub1, f2), (when, f1)] if full else [(short, f2), (when, f2)]
    y = pad
    for text, fnt in lines:
        bb = d.textbbox((pad, y), text, font=fnt)
        d.rectangle((bb[0] - 4 * t, bb[1] - 3 * t, bb[2] + 4 * t, bb[3] + 3 * t), fill=(0, 0, 0))
        d.text((pad, y), text, font=fnt, fill=(255, 255, 255))
        y = bb[3] + int(8 * t)
    # projection note
    if meta["projection"] == "perspective":
        Xc = meta["center"][0]
        hw = (W / 2) / meta["F"]
        hf = math.degrees(math.atan(Xc + hw) - math.atan(Xc - hw))
        note = "perspective view, horizontal field of view %.1f°" % hf
    else:
        note = ("fisheye view (azimuthal equidistant): the nucleus reaches %.0f° from the RG axis, too wide for a "
                "perspective view" % meta["theta_max_deg"])
        if not full:
            note = "fisheye view (azimuthal equidistant), nucleus up to %.0f° from the RG axis" % meta["theta_max_deg"]
    if fixed:
        note = ("fixed camera (comet frame): COPS position and RG axis of %s UTC (COPS altitude %s); perspective, "
                "horizontal field of view %.0f°" % (cam["utc"][11:19], fmt_alt(cam["altitude_target_m"]), hf)
                if full else "fixed camera: COPS at %s UTC (%s), perspective, %.0f°" % (
                    cam["utc"][11:19], fmt_alt(cam["altitude_target_m"]), hf))
    if detached:
        note = ("camera detached at %s UTC (COPS %s): position and axis frozen in the comet frame; perspective, "
                "horizontal field of view %.0f°" % (cam["detach_utc"][11:19], fmt_alt(cam["detach_altitude_m"]), hf)
                if full else "camera detached at %s UTC (%s), perspective, %.0f°" % (
                    cam["detach_utc"][11:19], fmt_alt(cam["detach_altitude_m"]), hf))
        if cam.get("pan"):                     # 2026-11-14: the moving camera (plan --pan-to)
            note = ("camera detached at %s UTC (COPS %s), then moving to a panoramic view; perspective, horizontal "
                    "field of view %.0f°" % (cam["detach_utc"][11:19], fmt_alt(cam["detach_altitude_m"]), hf)
                    if full else "camera detached at %s UTC (%s), moving; perspective, %.0f°" % (
                        cam["detach_utc"][11:19], fmt_alt(cam["detach_altitude_m"]), hf))
            if cam.get("pan_s", 0.0) >= 1.0:   # the final pose reached (held until the end)
                note = ("camera detached at %s UTC (COPS %s), now held at its final panoramic view; perspective, "
                        "horizontal field of view %.0f°" % (cam["detach_utc"][11:19], fmt_alt(cam["detach_altitude_m"]), hf)
                        if full else "camera detached at %s UTC (%s), final view held; perspective, %.0f°" % (
                            cam["detach_utc"][11:19], fmt_alt(cam["detach_altitude_m"]), hf))
    if attached_d:
        note = note + ("; camera at COPS until %s, then detached" % fmt_alt(cam["detach_altitude_m"]) if full else
                       " (camera at COPS until %s)" % fmt_alt(cam["detach_altitude_m"]))
    fn = font(int(round(20 * t)))
    bb = d.textbbox((pad, Hi - pad), note, font=fn, anchor="ls")
    d.rectangle((bb[0] - 4 * t, bb[1] - 3 * t, bb[2] + 4 * t, bb[3] + 3 * t), fill=(0, 0, 0))
    d.text((pad, Hi - pad), note, font=fn, fill=(255, 255, 255), anchor="ls")
    # scale bar at the sub-spacecraft point, and the point itself (+)
    if "m_per_px_at_sub" in meta:
        mpp = meta["m_per_px_at_sub"] / s
        L = nice_length(0.18 * Wi * mpp)
        npx = L / mpp
        x0, y0 = pad, bb[1] - int(22 * t)
        d.rectangle((x0, y0 - 4 * t, x0 + npx, y0), fill=(255, 255, 255))
        lab = ("%g km" % (L / 1000)) if L >= 1000 else ("%g m" % L)
        d.text((x0 + npx + 8 * t, y0 + 2 * t), lab + (" at the sub-spacecraft point" if NO_SUB_MARK else
                                                      " at the sub-spacecraft point (+)"), font=fn, fill=(255, 255, 255),
               anchor="ls")
        px, py = meta["sub_pixel"]
        px, py = px * s, py * s
        if 0 <= px < Wi and 0 <= py < Hi and not NO_SUB_MARK:
            r = 9 * t
            d.line((px - r, py, px + r, py), fill=(255, 255, 255), width=lw)
            d.line((px, py - r, px, py + r), fill=(255, 255, 255), width=lw)
    if fixed or detached:                      # scale bar at the image centre; the sub-spacecraft point; the gauge
        if "m_per_px_at_center" in meta:
            mpp = meta["m_per_px_at_center"] / s
            L = nice_length(0.18 * Wi * mpp)
            npx = L / mpp
            x0, y0 = pad, bb[1] - int(22 * t)
            d.rectangle((x0, y0 - 4 * t, x0 + npx, y0), fill=(255, 255, 255))
            lab = ("%g km" % (L / 1000)) if L >= 1000 else ("%g m" % L)
            d.text((x0 + npx + 8 * t, y0 + 2 * t), lab + (((" at the image centre; ring: COPS" if NO_SUB_MARK else
                                                            " at the image centre; +: sub-spacecraft point; ring: COPS"))
                                                           if full else " at the image centre"),
                   font=fn, fill=(255, 255, 255), anchor="ls")
        if meta.get("sub_visible") and meta.get("sub_pixel") and not NO_SUB_MARK:
            px, py = meta["sub_pixel"][0] * s, meta["sub_pixel"][1] * s
            if 0 <= px < Wi and 0 <= py < Hi:
                r = 9 * t
                for col, wd in (((0, 0, 0), lw + 2), ((255, 255, 255), lw)):
                    d.line((px - r, py, px + r, py), fill=col, width=wd)
                    d.line((px, py - r, px, py + r), fill=col, width=wd)
        if meta.get("gauge_visible") and meta.get("gauge_pixel"):
            px, py = meta["gauge_pixel"][0] * s, meta["gauge_pixel"][1] * s
            if 0 <= px < Wi and 0 <= py < Hi:
                r = 8 * t
                d.ellipse((px - r, py - r, px + r, py + r), outline=(0, 0, 0), width=lw + 2)
                d.ellipse((px - r, py - r, px + r, py + r), outline=(255, 120, 200), width=lw)
                d.text((px - r - 4 * t, py), "COPS" if fixed else "COPS (fit)", font=fn, fill=(255, 120, 200), anchor="rm",
                       stroke_width=max(1, int(round(2 * t))), stroke_fill=(0, 0, 0))
    # colour bar (right); the block is right-aligned inside the frame
    ch = int(Hi * 0.46)
    cw = int(26 * t)
    wmax = 1.0 if which == "rg" else 2.0
    ticks = np.linspace(0, 1, 5) if which == "rg" else np.linspace(0, 2, 5)
    tw = max(d.textlength("%g" % v, font=fn) for v in ticks)
    title_w = d.textlength("%s weight" % gname, font=fn)
    block_w = max(cw + 9 * t + tw, title_w)
    cx0 = int(Wi - pad - block_w)
    cy0 = int(Hi * 0.25)
    d.rectangle((cx0 - int(8 * t), cy0 - int(40 * t), Wi - pad + int(6 * t), cy0 + ch + int(16 * t)), fill=(0, 0, 0))
    img.paste(Image.fromarray(colorbar_img(which, ch, cw)), (cx0, cy0))
    for v in ticks:
        yy = cy0 + ch * (1 - v / wmax)
        d.line((cx0 + cw, yy, cx0 + cw + 6 * t, yy), fill=(255, 255, 255), width=max(1, int(t)))
        d.text((cx0 + cw + 9 * t, yy), "%g" % v, font=fn, fill=(255, 255, 255), anchor="lm")
    d.text((cx0, cy0 - int(12 * t)), "%s weight" % gname, font=fn, fill=(255, 255, 255), anchor="ls")
    if td is not None:                         # 2026-11-13: touchdown legend
        td_legend(img, d, td, k, t, pad, full, bool(meta.get("sc_hides_touchdown")))
    return img


def fmt_alt(a):
    return ("%.2f km" % (a / 1000)) if a >= 1000 else ("%.1f m" % a if a < 100 else "%.0f m" % a)


TS_LABEL_DEFAULT = "COPS data and the\nsoftplus activity model\n(final state, 5 cycles)"


def timeseries_background(table_path, alt_frames, width, height, label=TS_LABEL_DEFAULT):
    """Background of the NG / RG time-series panel and its pixel mapping for log10(altitude)."""
    plt = FIG.style()
    import matplotlib
    matplotlib.rcParams.update({"font.size": 15, "axes.labelsize": 15, "xtick.labelsize": 13, "ytick.labelsize": 13,
                                "legend.fontsize": 13})
    tab = FIG.read_table(table_path)
    alt = tab["altitude"]
    o = np.argsort(-alt)
    fig = plt.figure(figsize=(width / 100, height / 100), dpi=100)
    axs = [fig.add_axes([0.075, 0.555, 0.80, 0.40]), fig.add_axes([0.075, 0.11, 0.80, 0.40])]
    # late-visibility window (CLAUDE.md rule 11: shaded, labelled with its UTC and altitude range); the observation rows
    # of the table in figures.DEFAULT_LATE_WINDOW (= the template bundle's late_ng_not_rg interval)
    win = FIG.late_window(tab)
    for ax in axs:
        ax.axvspan(win["altitude_min_m"], win["altitude_max_m"], facecolor=FIG.LATE_FILL, edgecolor=FIG.LATE_EDGE,
                   lw=0.8, zorder=0)
    import matplotlib.transforms as mtransforms
    tr = mtransforms.blended_transform_factory(axs[0].transData, axs[0].transAxes)
    axs[0].text(math.sqrt(win["altitude_min_m"] * win["altitude_max_m"]), 0.05,
                "late visibility window: %s–%s UTC, %s–%s m" % (win["start"].split("T")[1], win["stop"].split("T")[1],
                                                               FIG._hm(win["altitude_max_m"]),
                                                               FIG._hm(win["altitude_min_m"])),
                transform=tr, ha="center", va="bottom", fontsize=13, color="0.15", zorder=5)
    for ax, g in zip(axs, ("NG", "RG")):
        obs, use = tab["ref" + g], tab["use" + g]
        ok = obs > 0
        ax.plot(alt[ok], obs[ok], "o", color="k", ms=2.6, mew=0, label="COPS %s" % g, zorder=3)
        ax.plot(alt[o], tab["p" + g][o], "-", color="#0072B2" if g == "NG" else "#D55E00", lw=1.8,
                label="model %s" % g, zorder=4)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(alt_frames[0] * 1.15, alt_frames[-1] * 0.85)
        ax.set_ylabel("%s pressure [Pa]" % g)
        ax.legend(loc="upper left", frameon=False, ncol=2)
    axs[0].tick_params(labelbottom=False)
    axs[1].set_xlabel("COPS altitude above the shape model [m] (descent →)")
    fig.text(0.885, 0.93, label, fontsize=13, va="top")
    fig.canvas.draw()
    bg = np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy()
    xmap = []
    for ax in axs:
        x0 = ax.transData.transform((alt_frames[0], 1.0))[0]
        x1 = ax.transData.transform((alt_frames[-1], 1.0))[0]
        bbox = ax.get_window_extent()
        xmap.append({"x0": x0, "x1": x1, "la0": math.log10(alt_frames[0]), "la1": math.log10(alt_frames[-1]),
                     "y_top": height - bbox.y1, "y_bot": height - bbox.y0})
    model = {g: (np.log10(alt[o]), np.log10(np.maximum(tab["p" + g][o], 1e-30))) for g in ("NG", "RG")}
    ylims = [ax.get_ylim() for ax in axs]
    plt.close(fig)
    return bg, xmap, model, ylims


def compose_one(args):
    k, work, out, ts = args[:4]
    td = args[4] if len(args) > 4 else None    # 2026-11-13: the touchdown flag (compose --touchdown)
    from PIL import Image, ImageDraw
    fd = Path(work) / "frames"
    meta = json.load(open(fd / "meta" / ("f%05d.json" % k)))
    for i, which in enumerate(("rg", "ng")):
        p = Path(out) / ("movie%d" % (i + 1)) / ("f%05d.png" % k)
        if not p.exists():
            img = Image.open(fd / ("%s_clean" % which) / ("f%05d.png" % k)).convert("RGB")
            if td is None:
                overlay(img, meta, which, 1.0, True).save(p, compress_level=1)
            else:
                if td_show(td, k, meta):
                    draw_flag(img, td["base_px"][k], td["up_px"][k], k, 1.0)
                overlay(img, meta, which, 1.0, True, td=td, k=k).save(p, compress_level=1)
    p3 = Path(out) / "movie3" / ("f%05d.png" % k)
    if not p3.exists():
        bg, xmap, model, ylims, axes_pix = ts
        canvas = Image.new("RGB", (W, H), (0, 0, 0))
        for i, which in enumerate(("rg", "ng")):
            img = Image.open(fd / ("%s_clean" % which) / ("f%05d.png" % k)).convert("RGB").resize((W // 2, H // 2),
                                                                                                  Image.LANCZOS)
            if td is None:
                canvas.paste(overlay(img, meta, which, 0.5, False, text_scale=0.75), (i * W // 2, 0))
            else:
                if td_show(td, k, meta):
                    draw_flag(img, 0.5 * np.asarray(td["base_px"][k]), td["up_px"][k], k, 0.75)
                canvas.paste(overlay(img, meta, which, 0.5, False, text_scale=0.75, td=td, k=k), (i * W // 2, 0))
        ts_img = Image.fromarray(bg.copy())
        d = ImageDraw.Draw(ts_img)
        la = math.log10(meta["altitude_m"])
        after = bool(meta.get("after_last_obs"))      # 2026-11-13: frames after the last COPS observation
        if after:
            la = xmap[0]["la1"]                       # the line stays at the end of the time series (10:39:26 UTC)
        for j, (m, g) in enumerate(zip(xmap, ("NG", "RG"))):
            x = m["x0"] + (la - m["la0"]) / (m["la1"] - m["la0"]) * (m["x1"] - m["x0"])
            d.line((x, m["y_top"], x, m["y_bot"]), fill=(200, 0, 0), width=3)
        tlab = "%s UTC, %s" % (meta["utc"][11:19], fmt_alt(meta["altitude_m"]))
        if after:
            tlab = "%s UTC, COPS %s: after the last observation (time series ends at 10:39:26 UTC)" % (
                meta["utc"][11:21], fmt_alt(meta["altitude_m"]))
        x = xmap[0]["x0"] + (la - xmap[0]["la0"]) / (xmap[0]["la1"] - xmap[0]["la0"]) * (xmap[0]["x1"] - xmap[0]["x0"])
        fn = font(18)
        bb = d.textbbox((x, xmap[0]["y_top"] - 4), tlab, font=fn, anchor="mb")
        shift = max(0, bb[2] - (W - 6)) - min(0, bb[0] - 6)
        d.text((x - shift, xmap[0]["y_top"] - 4), tlab, font=fn, fill=(200, 0, 0), anchor="mb")
        canvas.paste(ts_img, (0, H // 2))
        dd = ImageDraw.Draw(canvas)
        dd.line((W // 2, 0, W // 2, H // 2), fill=(255, 255, 255), width=2)
        canvas.save(p3, compress_level=1)
    return k


def cmd_compose(a):
    import multiprocessing as mp
    out = Path(a.out)
    for i in (1, 2, 3):
        (out / ("movie%d" % i)).mkdir(parents=True, exist_ok=True)
    P = load_prep(a.work)
    alt_ts = P["alt"][:int(P["n_obs_frames"])] if "n_obs_frames" in P else P["alt"]   # time series ends at the last observation
    bg, xmap, model, ylims = timeseries_background(a.table, alt_ts, W, H // 2,
                                                   a.ts_label.replace("\\n", "\n") if a.ts_label else TS_LABEL_DEFAULT)
    ts = (bg, xmap, model, ylims, None)
    ks = [int(v) for v in a.frames.split(",")] if a.frames else list(range(len(P["utc"])))
    missing = [k for k in ks if not (Path(a.work) / "frames" / "meta" / ("f%05d.json" % k)).exists()]
    if missing:
        raise SystemExit("%d frames not rendered (first %d)" % (len(missing), missing[0]))
    td = None
    if getattr(a, "touchdown", ""):            # 2026-11-13: the touchdown flag (new option)
        td = load_touchdown(a.touchdown)
        say("touchdown flag from %s: %s UTC, visible in %d of %d frames" % (a.touchdown, str(td["utc"]),
                                                                          int(np.sum(td["visible"])), len(P["utc"])))
        if getattr(a, "flag_on_top", False):   # 2026-11-13 (user): drawn on top, no depth test
            td["_on_top"] = True
            say("flag drawn on top (no depth test) where its point projects into the image: %d frames" % int(
                np.sum(td["in_front"] & td["in_image"])))
    if getattr(a, "no_sub_mark", False):      # 2026-11-14 (user): "remove the white cross" (new option)
        global NO_SUB_MARK
        NO_SUB_MARK = True
        say("sub-spacecraft mark (+) not drawn")
    jobs = [(k, a.work, str(out), ts) if td is None else (k, a.work, str(out), ts, td) for k in ks]
    with mp.get_context("fork").Pool(a.procs) as pool:
        for i, _ in enumerate(pool.imap_unordered(compose_one, jobs, chunksize=4)):
            if (i + 1) % 200 == 0:
                say("composed %d / %d" % (i + 1, len(ks)))
    say("compose done: %d frames" % len(ks))
    return 0


# ============================================================================================ encode
def cmd_encode(a):
    out = Path(a.out)
    fdir = Path(a.frames_dir)
    res = {}
    for i, name in ((1, "movie1_RG_view"), (2, "movie2_NG_view"), (3, "movie3_combined")):
        src = fdir / ("movie%d" % i) / "f%05d.png"
        for tag, extra in (("1920x1080", ["-crf", "18"]), ("1280x720_web", ["-vf", "scale=1280:720:flags=lanczos", "-crf", "26"])):
            dst = out / ("%s_%s.mp4" % (name, tag))
            if dst.exists():
                raise SystemExit("%s exists" % dst)
            cmd = [a.ffmpeg, "-hide_banner", "-loglevel", "error", "-framerate", str(FPS), "-i", str(src), "-c:v", "libx264",
                   "-preset", "slow"] + extra + ["-pix_fmt", "yuv420p", "-movflags", "+faststart", "-r", str(FPS), str(dst)]
            t0 = time.time()
            subprocess.run(cmd, check=True)
            res[dst.name] = {"command": cmd, "bytes": dst.stat().st_size, "seconds": time.time() - t0,
                             "sha256": sha256_file(dst)}
            say("encoded %s: %.1f MB in %.0f s" % (dst.name, dst.stat().st_size / 1e6, time.time() - t0))
    json.dump(res, open(out / "encode_manifest.json", "x"), indent=1)
    return 0


def cmd_frametable(a):
    """frames.tsv: one row per frame from frames/meta (geometry, framing, hole fill, GPU, time; masked frames also the
    visible-surface area fraction and nucleus-pixel fraction inside each gauge's mask and the weight range there)."""
    fd = Path(a.frames_dir) if a.frames_dir else Path(a.work) / "frames"
    out = Path(a.out)
    if out.exists():
        raise SystemExit("%s exists" % out)
    M = [json.load(open(fd / "meta" / ("f%05d.json" % k))) for k in range(len(load_prep(a.work)["utc"]))]
    masked = "masks" in M[0]
    fixed = "camera" in M[0]                   # 2026-11-13: fixed-camera frames (extra columns, camera check)
    if fixed:
        MD = M
        if M[0]["camera"].get("mode") == "detach":         # detach-200 m set: the detached frames only
            MD = [m for m in M if not m["camera"].get("attached")]
        if MD[0]["camera"].get("pan"):         # 2026-11-14: moving camera: per frame the plan's pose; same pose, same view
            return frametable_pan(a, fd, out, M, MD)
        ref = (MD[0]["camera"], MD[0]["face_sha256"])
        bad = [m["k"] for m in MD if (m.get("camera"), m.get("face_sha256")) != ref]
        chk = {"frames": len(MD), "first_frame": MD[0]["k"], "frames_with_a_different_camera_or_face_buffer": bad, "camera": ref[0],
               "face_buffer_sha256": ref[1], "hosts_devices": sorted({"%s:%s" % (m["host"], m["device"]) for m in M}),
               "jobs": sorted({str(m["job"]) for m in M})}
        json.dump(chk, open(str(out) + ".camera_check.json", "x"), indent=1)
        say("fixed camera check: %d frames, %d differ (camera or face buffer); %d host:device pairs" % (
            len(M), len(bad), len(chk["hosts_devices"])))
        if bad:
            raise SystemExit("fixed camera: the camera or the face buffer differs in %d frames" % len(bad))
    cols = ["frame", "movie_s", "utc", "altitude_m", "projection", "theta_max_deg", "F_px_per_unit", "center_X",
            "center_Y", "m_per_px_at_sub", "sub_distance_m", "lit_frac_day", "enclosed_px", "filled_px", "no_hit_px"]
    if masked:
        for g in ("RG", "NG"):
            cols += ["%s_mask" % g, "%s_visible_area_frac_in_mask" % g, "%s_visible_area_m2_in_mask" % g,
                     "%s_nucleus_pixel_frac_in_mask" % g, "%s_w_min_in_mask" % g, "%s_w_max_in_mask" % g]
    if fixed:
        for g in ("RG", "NG"):
            cols += ["%s_fov_pixel_frac" % g, "%s_fov_facing_pixel_frac" % g, "%s_fov_area_m2" % g,
                     "%s_fov_facing_area_m2" % g, "%s_los_blocked_queries" % g, "%s_los_near_plane_facets" % g,
                     "%s_los_map_px" % g]
        cols += ["sub_visible", "gauge_visible", "gauge_px_x", "gauge_px_y", "gauge_distance_from_camera_m"]
    cols += ["host", "device", "t_total_s"]
    with open(out, "x") as f:
        f.write("# one row per movie frame (%s/meta); movie time = frame / %d s; masked frames: visible-surface area "
                "fraction inside the gauge's mask = drawn facets whose centroid direction is inside, area-weighted\n"
                % (fd, FPS))
        if fixed:
            f.write("# fixed camera: 'in mask' = the gauge sees the surface (inside its mask (IK FOV) from the true gauge "
                    "position, facet facing the gauge, line of sight); pixel fractions of the camera's nucleus pixels "
                    "(supersampled); areas of the camera-visible facets by their centroid\n")
        f.write("\t".join(cols) + "\n")
        for m in M:
            h = m["holes"]
            row = [m["k"], "%.3f" % (m["k"] / FPS), m["utc"], "%.4f" % m["altitude_m"], m["projection"],
                   "%.3f" % m["theta_max_deg"], "%.3f" % m["F"], "%.6f" % m["center"][0], "%.6f" % m["center"][1],
                   "%.5g" % m.get("m_per_px_at_sub", float("nan")), "%.4f" % m.get("sub_distance_m", float("nan")),
                   "%.4f" % m["lit_frac_day"], h["enclosed_px"], h["filled_px"], h["no_hit_px"]]
            if masked:
                for g in ("RG", "NG"):
                    st = m["masks"][g]
                    row += [st.get("mask"), "%.5f" % st.get("visible_area_frac", float("nan")),
                            "%.6g" % st.get("visible_area_m2_in_mask", float("nan")),
                            "%.5f" % st.get("nucleus_pixel_frac", float("nan")),
                            "%.5f" % st.get("weight_min_on_nucleus_in_mask", float("nan")),
                            "%.5f" % st.get("weight_max_on_nucleus_in_mask", float("nan"))]
            if fixed:
                for g in ("RG", "NG"):
                    st = m["masks"][g]
                    ls = st.get("los", {})
                    nan = float("nan")
                    row += ["%.5f" % st.get("fov_pixel_frac", nan), "%.5f" % st.get("fov_facing_pixel_frac", nan),
                            "%.6g" % st.get("fov_area_m2", nan), "%.6g" % st.get("fov_facing_area_m2", nan),
                            ls.get("blocked", 0), ls.get("near_plane_facets_tested", 0),
                            "x".join(str(v) for v in ls.get("map_px", []))]
                gp = m.get("gauge_pixel") or [float("nan"), float("nan")]
                row += [int(bool(m.get("sub_visible"))), int(bool(m.get("gauge_visible"))), "%.2f" % gp[0], "%.2f" % gp[1],
                        "%.4f" % m.get("gauge_distance_from_camera_m", float("nan"))]
            row += [m["host"], m["device"], "%.3f" % m["t_total_s"]]
            f.write("\t".join(str(v) for v in row) + "\n")
    say("wrote %s" % out)
    return 0


def frametable_pan(a, fd, out, M, MD):
    """frametable of the moving-camera (pan) frames (2026-11-14): checks that each detached frame used the plan's pose
    of its frame (position, forward, up, focal length) and that frames with the same pose have the same face buffer;
    then writes the table as cmd_frametable (camera columns added)."""
    Q = np.load(Path(a.work) / "plan.npz", allow_pickle=False)
    bad, groups = [], {}
    for m in MD:
        k = m["k"]
        f = unit(Q["cam_b_rg_frames"][k])
        up = unit(np.asarray(Q["cam_b_ng_frames"][k], float) - np.dot(Q["cam_b_ng_frames"][k], f) * f)
        c = m["camera"]
        if (np.max(np.abs(np.asarray(c["x"]) - Q["cam_x_frames"][k])) > 0 or np.max(np.abs(np.asarray(c["f"]) - f)) > 1e-12
                or np.max(np.abs(np.asarray(c["up"]) - up)) > 1e-12 or m["F"] != float(Q["F"][k])):
            bad.append(k)
        groups.setdefault(json.dumps([c["x"], c["f"], c["up"], m["F"]]), set()).add(m["face_sha256"])
    fb = [g for g, v in groups.items() if len(v) > 1]
    chk = {"frames": len(MD), "first_frame": MD[0]["k"], "frames_not_matching_the_plan_pose": bad,
           "distinct_poses": len(groups), "poses_with_more_than_one_face_buffer": len(fb),
           "pan_start_frame": int(Q["pan_start_frame"]), "pan_end_frame": int(Q["pan_end_frame"]),
           "hosts_devices": sorted({"%s:%s" % (m["host"], m["device"]) for m in M}),
           "jobs": sorted({str(m["job"]) for m in M})}
    json.dump(chk, open(str(out) + ".camera_check.json", "x"), indent=1)
    say("pan camera check: %d detached frames, %d not matching the plan, %d poses, %d with more than one face buffer" % (
        len(MD), len(bad), len(groups), len(fb)))
    if bad or fb:
        raise SystemExit("pan camera: %d frames differ from the plan, %d poses with several face buffers" % (len(bad), len(fb)))
    cols = ["frame", "movie_s", "utc", "altitude_m", "projection", "F_px_per_unit", "fov_horizontal_deg", "pan_s",
            "cam_x", "cam_y", "cam_z", "lit_frac_day", "enclosed_px", "filled_px", "no_hit_px"]
    for g in ("RG", "NG"):
        cols += ["%s_mask" % g, "%s_visible_area_frac_in_mask" % g, "%s_visible_area_m2_in_mask" % g,
                 "%s_nucleus_pixel_frac_in_mask" % g, "%s_fov_area_m2" % g, "%s_fov_facing_area_m2" % g]
    cols += ["gauge_px_x", "gauge_px_y", "gauge_distance_from_camera_m", "host", "device", "t_total_s"]
    nan = float("nan")
    with open(out, "x") as f:
        f.write("# one row per movie frame (%s/meta); movie time = frame / %d s; detach camera moving after the detach "
                "(plan --pan-to); attached frames: camera at COPS\n" % (fd, FPS))
        f.write("\t".join(cols) + "\n")
        for m in M:
            h = m["holes"]
            c = m.get("camera", {})
            x = c.get("x", [nan] * 3) if not c.get("attached") else list(Q["cam_x_frames"][m["k"]])
            row = [m["k"], "%.3f" % (m["k"] / FPS), m["utc"], "%.4f" % m["altitude_m"], m["projection"], "%.3f" % m["F"],
                   "%.3f" % (2 * math.degrees(math.atan((W / 2.0) / m["F"]))) if m["projection"] == "perspective" else "nan",
                   "%.6f" % c.get("pan_s", 0.0), "%.4f" % x[0], "%.4f" % x[1], "%.4f" % x[2], "%.4f" % m["lit_frac_day"],
                   h["enclosed_px"], h["filled_px"], h["no_hit_px"]]
            for g in ("RG", "NG"):
                st = m["masks"][g]
                row += [st.get("mask"), "%.5f" % st.get("visible_area_frac", nan),
                        "%.6g" % st.get("visible_area_m2_in_mask", nan), "%.5f" % st.get("nucleus_pixel_frac", nan),
                        "%.6g" % st.get("fov_area_m2", nan), "%.6g" % st.get("fov_facing_area_m2", nan)]
            gp = m.get("gauge_pixel") or [nan, nan]
            row += ["%.2f" % gp[0], "%.2f" % gp[1], "%.4f" % m.get("gauge_distance_from_camera_m", nan), m["host"],
                    m["device"], "%.3f" % m["t_total_s"]]
            f.write("\t".join(str(v) for v in row) + "\n")
    say("wrote %s" % out)
    return 0


def cmd_validate_fov(a):
    """IK-FOV masks: (1) the mask (Geometry.cpp InsideFOV on the body-frame bounds) against an independent analytic
    test in the instrument frame (rectangle half-angles) at sample directions, at several frame times; (2) the
    acceptance-weight range inside each FOV (corners); (3) kernel records of base rows: share of the kernel sums (Kn for
    NG, KPhi for RG; facet-centroid directions) inside each FOV; (4) masked test frames."""
    import cupy as cp
    from PIL import Image
    work = Path(a.work)
    P = load_prep(work)
    Q = np.load(work / "plan.npz", allow_pickle=False)
    FV = load_fov(work)
    if FV is None:
        raise SystemExit("run cops_movies.py fov first")
    rnd = Renderer(a.geometry, P, a.device)
    info = json.load(open(work / "fov_checks.json"))
    res = {"fov_checks": info["checks"], "fov_info": {g: {k: info["info"][g][k] for k in
                                                          ("id", "shape", "frame", "half_angles_deg",
                                                           "corner_angle_deg")} for g in ("NG", "RG")},
           "analytic_test": [], "weight_range": {}, "kernel_records": [], "test_frames": {}}
    rng = np.random.default_rng(20261113)
    for k in [int(v) for v in a.test_frames.split(",")]:
        for g in ("RG", "NG"):
            gi = info["info"][g]
            bsI = np.array(gi["boresight"])
            e1, e2 = (np.array(v) for v in gi["axes"])
            h1, h2 = gi["half_angles_deg"]
            # half-angle of e1 / e2 from the corner (sorted list: assign by recomputing)
            c0 = np.array(gi["bounds"][0])
            ha = {0: np.degrees(np.arctan2(abs(c0.dot(e1)), c0.dot(bsI))),
                  1: np.degrees(np.arctan2(abs(c0.dot(e2)), c0.dot(bsI)))}
            al = np.concatenate([np.arange(-89.5, 90.0, 0.5), [ha[0] - 1e-6, ha[0] + 1e-6, -ha[0] + 1e-6, -ha[0] - 1e-6]])
            be = np.concatenate([np.arange(-89.5, 90.0, 0.5), [ha[1] - 1e-6, ha[1] + 1e-6, -ha[1] + 1e-6, -ha[1] - 1e-6]])
            A_, B_ = np.meshgrid(np.radians(al), np.radians(be), indexing="ij")
            dI = (np.tan(A_)[..., None] * e1 + np.tan(B_)[..., None] * e2 + bsI).reshape(-1, 3)
            rnd_dirs = rng.normal(size=(200000, 3))
            dI = np.vstack([dI, -dI[:5000], rnd_dirs])
            dI /= np.linalg.norm(dI, axis=1)[:, None]
            fwd = dI.dot(bsI)
            aa = np.degrees(np.arctan2(dI.dot(e1), fwd))
            bb = np.degrees(np.arctan2(dI.dot(e2), fwd))
            expect = (fwd > 0) & (np.abs(aa) <= ha[0]) & (np.abs(bb) <= ha[1])
            R = FV["frame_rot_%s" % g.lower()][k]
            dB = dI @ R.T
            d = tuple(cp.asarray(dB[:, i]) for i in range(3))
            got = cp.asnumpy(rnd.mask_dirs(d, g, "ik_fov", P["b_rg"][k], P["b_ng"][k],
                                           FV["frame_bounds_%s" % g.lower()][k]))
            near = (np.abs(np.abs(aa) - ha[0]) < 1e-9) | (np.abs(np.abs(bb) - ha[1]) < 1e-9)
            mis = (got != expect) & ~near
            res["analytic_test"].append({"frame": k, "gauge": g, "directions": int(len(dI)), "inside": int(expect.sum()),
                                         "mismatches": int(mis.sum()), "mismatches_within_1e-9_deg_of_an_edge":
                                         int(((got != expect) & near).sum()), "half_angles_deg": [ha[0], ha[1]]})
            say("ik_fov test frame %d %s: %d directions (%d inside), mismatches %d (edge %d)" % (
                k, g, len(dI), expect.sum(), mis.sum(), ((got != expect) & near).sum()))
    for g in ("RG", "NG"):
        md = mask_design(g, "ik_fov", FV)
        res["weight_range"][g] = md
        say("%s inside the IK FOV: corner %.4f deg from the boresight, weight range %s" % (
            g, md["max_angle_deg"], md.get("weight_range")))
    # kernel records of base rows: share of Kn (NG) / KPhi (RG) inside each FOV (facet-centroid direction)
    base = Path(a.base)
    ip_all = list(P["rows_ipoint"].astype(int))
    cent = cp.asnumpy(rnd.cent)
    for ip in [int(v) for v in a.rows.split(",")]:
        i = ip_all.index(ip)
        kp, h, rec = read_kernel(base / "kernels" / ("point=%04d.bin" % ip))
        x = np.asarray(h["sc"], float)
        dd = cent[rec["face"]] - x[None, :]
        dd /= np.linalg.norm(dd, axis=1)[:, None]
        d = tuple(cp.asarray(dd[:, j]) for j in range(3))
        row = {"iPoint": ip, "utc": str(P["rows_utc"][i]), "altitude_m": float(P["rows_alt"][i]), "records": int(len(rec))}
        for g, key in (("RG", "KPhi"), ("NG", "Kn")):
            m = cp.asnumpy(rnd.mask_dirs(d, g, "ik_fov", P["rows_b_rg"][i], P["rows_b_ng"][i],
                                         FV["rows_bounds_%s" % g.lower()][i]))
            v = np.asarray(rec[key], float)
            row["%s_%s_share_inside_ik_fov" % (g, key)] = float(v[m].sum() / max(v.sum(), 1e-300))
            row["%s_records_inside" % g] = int(m.sum())
        res["kernel_records"].append(row)
        say("records row %d (%.1f m): %s" % (ip, row["altitude_m"], json.dumps(
            {kk: vv for kk, vv in row.items() if "share" in kk})))
    # masked test frames
    td = work / "test_frames_masked"
    td.mkdir(exist_ok=True)
    for k in [int(v) for v in a.test_frames.split(",")]:
        vw = rnd.view(P["x"][k], P["b_rg"][k], P["b_ng"][k], P["sun"][k], str(Q["proj"][k]), float(Q["F"][k]),
                      tuple(Q["center"][k]))
        pw = rnd.pixel_weights(vw)
        dirs = rnd.pixel_dirs_world(vw)
        for which, g in (("rg", "RG"), ("ng", "NG")):
            mk = rnd.mask_dirs(dirs, g, "ik_fov", P["b_rg"][k], P["b_ng"][k], FV["frame_bounds_%s" % which][k])
            img = rnd.colorize_masked(vw, which, pw, mk)
            Image.fromarray(img).save(td / ("f%05d_%s_ikfov.png" % (k, which)), compress_level=1)
            st = rnd.mask_stats(vw, P["x"][k], mk, pw[0] if which == "rg" else pw[1])
            res["test_frames"]["%d_%s" % (k, g)] = st
            say("test frame %d %s: %s" % (k, g, json.dumps(st)))
        del vw, pw, dirs
    json.dump(res, open(work / ("validation_fov_%s.json" % socket.gethostname()), "w"), indent=1, default=float)
    bad = sum(t["mismatches"] for t in res["analytic_test"])
    say("validate_fov done: %d mismatches" % bad)
    return 0 if bad == 0 else 1


def cmd_diagnose(a):
    """Dark-surface diagnosis of rendered frames (2026-11-13, user: "starting from about 20 m altitude, the movies have
    some artifact on the surface that I thought is the spacecraft shadow. what is it actually?"). Per frame: the face
    buffer as rendered; every nucleus pixel at the ambient level (brightness = AMBIENT: no direct sunlight) is classified
    by its facet: 'cos <= 0' (facet turned away from the Sun), 'cast shadow' (cos > 0, the shadow map says shadow and the
    EXACT shadow ray says shadow), 'false shadow' (shadow map shadow, exact ray lit); and 'false lit' (shadow map lit,
    exact ray shadowed) among the lit pixels; holes = empty pixels enclosed by the nucleus. Exact shadow ray as
    nonneg_maps_shadow.cpp / Source.cpp: x0 = centroid + 1e-6 MaxEdge normal, toward the Sun, nucleus only (here a
    brute-force first-hit over all facets, render.cu ray_first_hit, FP64; a hit on the facet itself is ignored).
    Writes <out>/artifact_check_f<k>.png (annotated), artifact_check.json and artifact_check.txt."""
    import cupy as cp
    from PIL import Image, ImageDraw
    work = Path(a.work)
    P = load_prep(work)
    Q = np.load(work / "plan.npz", allow_pickle=False)
    rnd = Renderer(a.geometry, P, a.device)
    fg = np.memmap(Path(a.geometry) / "mesh_facegeom.f8", dtype="<f8", mode="r").reshape(-1, 8)
    maxedge = cp.asarray(np.ascontiguousarray(fg[:, 7]))
    del fg
    out = Path(a.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    res = {"definition": cmd_diagnose.__doc__, "frames": {}}
    T = []
    boxes = {}
    for spec in (a.box or []):
        k, x0, x1, y0, y1 = (int(v) for v in spec.split(","))
        boxes.setdefault(k, []).append((x0, x1, y0, y1))
    for k in [int(v) for v in a.frames.split(",")]:
        x = np.asarray(P["x"][k], float)
        sun = unit(P["sun"][k])
        vw = rnd.view(x, P["b_rg"][k], P["b_ng"][k], sun, str(Q["proj"][k]), float(Q["F"][k]), tuple(Q["center"][k]))
        face = vw["face"]
        cosz, lit = vw["cos"], vw["lit"]
        nuc = face >= 0
        fpix = face[nuc]
        drawn, npix = cp.unique(fpix, return_counts=True)
        # exact shadow ray for the drawn facets with cos > 0
        cand = drawn[cosz[drawn] > 0.0]
        o = rnd.cent[cand] + rnd.norm[cand] * (1e-6 * maxedge[cand])[:, None]
        dsun = cp.broadcast_to(cp.asarray(sun)[None, :], o.shape)
        key = cp.empty(len(cand), dtype=cp.uint64)
        for r0 in range(0, len(cand), 512):
            r1 = min(r0 + 512, len(cand))
            best = cp.full(r1 - r0, np.iinfo(np.uint64).max, dtype=cp.uint64)
            rnd.rayhit((256, r1 - r0), (256,), (rnd.nodes, rnd.faces, np.int64(rnd.nf), cp.ascontiguousarray(o[r0:r1]),
                                                cp.ascontiguousarray(dsun[r0:r1]), best))
            key[r0:r1] = best
        hit = key != cp.uint64(np.iinfo(np.uint64).max)
        hf = (key & cp.uint64(0xffffffff)).astype(cp.int64)
        self_hit = hit & (hf == cand)
        exact_sh = cp.zeros(rnd.nf, dtype=bool)
        exact_sh[cand] = hit & ~self_hit
        occ = cp.full(rnd.nf, -1, dtype=cp.int64)              # the facet that blocks the Sun ray, and the distance
        occ[cand] = cp.where(hit & ~self_hit, hf, -1)
        occ_t = cp.zeros(rnd.nf)
        occ_t[cand] = ((key >> cp.uint64(32)).astype(cp.uint32).view(cp.float32)).astype(cp.float64)
        # per-pixel classes
        cpx = cosz[face.clip(0)]
        lpx = lit[face.clip(0)]
        epx = exact_sh[face.clip(0)]
        cls = cp.full(face.shape, -1, dtype=cp.int8)          # -1 sky, 0 directly lit (map), 1 cos<=0, 2 cast, 3 false shadow
        cls[nuc] = 0
        cls[nuc & (cpx <= 0.0)] = 1
        cls[nuc & (cpx > 0.0) & ~lpx & epx] = 2
        cls[nuc & (cpx > 0.0) & ~lpx & ~epx] = 3
        false_lit = nuc & (cpx > 0.0) & lpx & epx
        names = {1: "cos <= 0", 2: "cast shadow (map and exact)", 3: "false shadow (map shadow, exact lit)"}
        ntot = int(nuc.sum())
        fr = {"nucleus_pixels_ss": ntot, "holes": vw["holes"], "altitude_m": float(P["alt"][k]), "utc": str(P["utc"][k]),
              "drawn_facets": int(len(drawn)), "facet_pixels_median_ss": float(cp.asnumpy(cp.median(npix))),
              "facet_pixels_max_ss": int(cp.asnumpy(npix.max())),
              "exact_rays": int(len(cand)), "exact_self_hits_ignored": int(self_hit.sum()),
              "false_lit_pixel_frac": float(cp.asnumpy(false_lit.sum() / max(ntot, 1)))}
        for c, n in names.items():
            fr["pixel_frac_" + n] = float(cp.asnumpy((cls == c).sum() / max(ntot, 1)))
        # facets under the dark pixels: by class, with cos, size in pixels (output pixels = ss^2 supersampled)
        dark = nuc & (cls > 0)
        fd, nd = cp.unique(face[dark], return_counts=True)
        fdh, ndh = cp.asnumpy(fd), cp.asnumpy(nd) / float(SS * SS)
        cd = cp.asnumpy(cosz[fd])
        ld = cp.asnumpy(lit[fd])
        ed = cp.asnumpy(exact_sh[fd])
        o_ = np.argsort(-ndh)
        fr["largest_dark_facets"] = [{"facet": int(fdh[i]), "pixels": float(ndh[i]), "cos_sza": float(cd[i]),
                                      "map_lit": bool(ld[i]), "exact_shadow": bool(ed[i]),
                                      "class": names[1] if cd[i] <= 0 else names[2] if ed[i] else names[3]}
                                     for i in o_[:15]]
        fr["boxes"] = []
        for (x0, x1, y0, y1) in boxes.get(k, []):
            sub = cls[y0 * SS:y1 * SS, x0 * SS:x1 * SS]
            fsub = face[y0 * SS:y1 * SS, x0 * SS:x1 * SS]
            dm = sub > 0
            fb, nb = cp.unique(fsub[dm], return_counts=True)
            fbh = cp.asnumpy(fb)
            fr["boxes"].append({"box_px": [x0, x1, y0, y1], "dark_pixel_frac": float(cp.asnumpy(dm.mean())),
                                "class_counts": {names[c]: int(cp.asnumpy((sub == c).sum())) for c in names},
                                "facets": [{"facet": int(f), "pixels": float(n) / SS ** 2,
                                            "cos_sza": float(cp.asnumpy(cosz[f])), "map_lit": bool(cp.asnumpy(lit[f])),
                                            "exact_shadow": bool(cp.asnumpy(exact_sh[f])),
                                            "occluder_facet": int(cp.asnumpy(occ[f])),
                                            "occluder_distance_m": float(cp.asnumpy(occ_t[f])),
                                            "occluder_height_above_facet_m": float(cp.asnumpy(
                                                cp.dot(rnd.cent[occ[f]] - rnd.cent[f], rnd.norm[f]))) if int(cp.asnumpy(occ[f])) >= 0 else None,
                                            "facet_distance_from_camera_m": float(np.linalg.norm(cp.asnumpy(rnd.cent[f]) - x)),
                                            "facet_area_m2": float(rnd.area_h[int(f)])}
                                           for f, n in sorted(zip(fbh, cp.asnumpy(nb)), key=lambda t: -t[1])[:20]]})
        res["frames"][str(k)] = fr
        T.append("frame %d (%s UTC, %.1f m): dark pixels: %s; false lit %.4f; holes %s; drawn facets %d, median facet "
                 "size %.1f output px (max %.0f)" % (
                     k, fr["utc"][11:19], fr["altitude_m"],
                     ", ".join("%s %.4f" % (n, fr["pixel_frac_" + n]) for n in names.values()),
                     fr["false_lit_pixel_frac"], json.dumps(vw["holes"]), fr["drawn_facets"],
                     fr["facet_pixels_median_ss"] / SS ** 2, fr["facet_pixels_max_ss"] / SS ** 2))
        for b in fr["boxes"]:
            T.append("  box x %d-%d y %d-%d: dark fraction %.3f, %s; largest facets: %s" % (
                *b["box_px"], b["dark_pixel_frac"], json.dumps(b["class_counts"]),
                "; ".join("%d (%.0f px, cos %.3f, map %s, exact %s)" % (f["facet"], f["pixels"], f["cos_sza"],
                                                                        "lit" if f["map_lit"] else "shadow",
                                                                        "shadow" if f["exact_shadow"] else "lit")
                          for f in b["facets"][:6])))
        # annotated still: the rendered RG frame, dark classes outlined in colour
        src = Path(a.frames_dir) / "rg_clean" / ("f%05d.png" % k)
        img = np.asarray(Image.open(src).convert("RGB")).astype(float) / 255.0
        cl_ds = cp.asnumpy(cls.reshape(H, SS, W, SS).max(axis=(1, 3)))
        colors = {1: (0.2, 0.5, 1.0), 2: (0.7, 0.3, 1.0), 3: (1.0, 0.1, 0.1)}
        ov = img.copy()
        for c, col in colors.items():
            m = cl_ds == c
            ov[m] = 0.45 * img[m] + 0.55 * np.array(col)
        im = Image.fromarray((np.clip(ov, 0, 1) * 255 + 0.5).astype(np.uint8))
        d = ImageDraw.Draw(im)
        fnt = font(22)
        lines = ["frame %d, %s UTC, %.1f m: nucleus pixels at the ambient level, by facet:" % (k, fr["utc"][11:19],
                                                                                                fr["altitude_m"]),
                 "blue: cos SZA <= 0 (facet turned away from the Sun)",
                 "violet: cast shadow (shadow map and exact shadow ray agree)",
                 "red: shadow-map shadow, exact ray lit (false shadow)"]
        y = 12
        for t in lines:
            bb = d.textbbox((12, y), t, font=fnt)
            d.rectangle((bb[0] - 4, bb[1] - 3, bb[2] + 4, bb[3] + 3), fill=(0, 0, 0))
            d.text((12, y), t, font=fnt, fill=(255, 255, 255))
            y = bb[3] + 8
        for (x0, x1, y0, y1) in boxes.get(k, []):
            d.rectangle((x0, y0, x1, y1), outline=(255, 255, 0), width=2)
        im.save(out / ("artifact_check_f%05d.png" % k))
        say(T[-1 - len(fr["boxes"])])
        del vw, face, cls
    json.dump(res, open(out / "artifact_check.json", "x"), indent=1, default=float)
    (out / "artifact_check.txt").write_text("\n".join(T) + "\n")
    print("\n".join(T))
    return 0


def cmd_validate_sc(a):
    """spacecraft shadow checks (movies only): (1) geometry vs SPICE (sc_checks.json); (2) per test frame the bus shadow
    on the mesh (facet centroids within 60 m, bus only, 1024 disc samples) vs the analytic intersection of the line from
    the bus centre along -Sun with the mesh (ray_first_hit); (3) disc-sampling convergence (16, 64, 256 vs 1024 samples);
    (4) synthetic plane perpendicular to the Sun at distance L below the bus: geometric (point-Sun) shadow width W, umbra
    and penumbra widths along a line vs the analytic W - 2 alpha L, 2 alpha L (umbra gone beyond L = W / (2 alpha))."""
    import cupy as cp
    work = Path(a.work)
    P = load_prep(work)
    SCZ = np.load(work / "sc.npz", allow_pickle=False)
    rnd = Renderer(a.geometry, P, a.device)
    res = {"geometry_checks": json.load(open(work / "sc_checks.json")), "frames": [], "synthetic": []}
    for k in [int(v) for v in a.test_frames.split(",")]:
        sun = unit(P["sun"][k])
        sh = {n: ScShadow(cp, SCZ, k, sun, n) for n in (16, 64, 256, 1024)}
        s0 = sh[1024]
        cb = s0.x + s0.R @ (0.5 * (s0.bmin + s0.bmax))
        best = cp.full(1, np.iinfo(np.uint64).max, dtype=cp.uint64)
        rnd.rayhit((256, 1), (256,), (rnd.nodes, rnd.faces, np.int64(rnd.nf), cp.asarray(cb[None, :]),
                                      cp.asarray(-sun[None, :]), best))
        kk = int(cp.asnumpy(best)[0])
        if kk == np.iinfo(np.uint64).max:
            res["frames"].append({"frame": k, "note": "the line from the bus along -Sun misses the nucleus"})
            continue
        f0 = kk & 0xffffffff
        c0, n0 = cp.asnumpy(rnd.cent[f0]), cp.asnumpy(rnd.norm[f0])
        t0 = float(np.dot(c0 - cb, n0) / np.dot(-sun, n0))
        S0 = cb - sun * t0
        near = cp.nonzero(cp.sum((rnd.cent - cp.asarray(S0)[None, :]) ** 2, axis=1) < 60.0 ** 2)[0]
        pts = rnd.cent[near]
        fb = {n: cp.asnumpy(sh[n].blocked(pts, ("bus",))[0]) for n in sh}
        fa = {n: cp.asnumpy(sh[n].blocked(pts)[0]) for n in sh}
        area = rnd.area_h[cp.asnumpy(near)]
        w = fb[1024] * area
        ent = {"frame": k, "utc": str(P["utc"][k]), "altitude_m": float(P["alt"][k]), "sa_held": bool(SCZ["held"][k]),
               "solar_angular_radius_deg": math.degrees(s0.alpha), "bus_to_surface_along_sun_m": t0,
               "analytic_bus_shadow_centre_m": S0.tolist(), "facets_within_60m": int(len(near))}
        if w.sum() > 0:
            Sb = (cp.asnumpy(pts) * w[:, None]).sum(axis=0) / w.sum()
            ent["bus_shadow_centroid_m"] = Sb.tolist()
            ent["centroid_vs_analytic_m"] = float(np.linalg.norm(Sb - S0))
            ent["bus_shadow_area_m2"] = float(w.sum())
            ent["bus_umbra_facets"] = int((fb[1024] >= 1.0).sum())
        ent["convergence_vs_1024"] = {str(n): {"max_abs": float(np.abs(fa[n] - fa[1024]).max()),
                                               "mean_abs_on_shadowed": float(np.abs(fa[n] - fa[1024])[fa[1024] > 0].mean())
                                               if (fa[1024] > 0).any() else 0.0} for n in (16, 64, 256)}
        ent["shadowed_facets_bus_and_wings"] = int((fa[1024] > 0).sum())
        res["frames"].append(ent)
        say("sc frame %d: %s" % (k, json.dumps({kk_: vv for kk_, vv in ent.items() if kk_ not in (
            "analytic_bus_shadow_centre_m", "bus_shadow_centroid_m")})))
    # synthetic plane perpendicular to the Sun below the bus (frame of the last test frame)
    k = int(a.test_frames.split(",")[-1])
    sun = unit(P["sun"][k])
    s1 = ScShadow(cp, SCZ, k, sun, 1)
    sD = ScShadow(cp, SCZ, k, sun, 1024)
    cb = sD.x + sD.R @ (0.5 * (sD.bmin + sD.bmax))
    e = unit(sD.R @ np.array([1.0, 0, 0]) - np.dot(sD.R @ np.array([1.0, 0, 0]), sun) * sun)
    for L in (20.0, 100.0, 500.0, 800.0, 850.0, 900.0, 2000.0):
        half = 4.0 + 3.0 * sD.alpha * L
        xs = np.linspace(-half, half, 8001)
        pts = cb[None, :] - sun[None, :] * L + xs[:, None] * e[None, :]
        b1 = cp.asnumpy(s1.blocked(cp.asarray(pts), ("bus",))[0])
        bD = cp.asnumpy(sD.blocked(cp.asarray(pts), ("bus",))[0])
        dx = xs[1] - xs[0]
        Wg = (b1 > 0.5).sum() * dx
        umb = (bD >= 1.0 - 1e-12).sum() * dx
        pen = ((bD > 0) & (bD < 1.0 - 1e-12)).sum() * dx / 2.0
        ent = {"L_m": L, "geometric_width_m": Wg, "umbra_width_m": umb, "penumbra_width_each_side_m": pen,
               "analytic_umbra_m": max(Wg - 2 * sD.alpha * L, 0.0), "analytic_penumbra_m": 2 * sD.alpha * L,
               "umbra_length_analytic_m": Wg / (2 * sD.alpha), "max_blocked": float(bD.max())}
        res["synthetic"].append(ent)
        say("synthetic L %.0f m: %s" % (L, json.dumps(ent)))
    json.dump(res, open(work / ("validation_sc_%s.json" % socket.gethostname()), "w"), indent=1, default=float)
    say("validate_sc done")
    return 0


def _cmp_one(args):
    k, a, b = args
    from PIL import Image
    out = [k]
    for which in ("rg", "ng"):
        x = np.asarray(Image.open(Path(a) / ("%s_clean" % which) / ("f%05d.png" % k)))
        y = np.asarray(Image.open(Path(b) / ("%s_clean" % which) / ("f%05d.png" % k)))
        out.append(int((x != y).any(axis=2).sum()))
    ma = json.load(open(Path(a) / "meta" / ("f%05d.json" % k)))
    mb = json.load(open(Path(b) / "meta" / ("f%05d.json" % k)))
    skip = {"host", "job", "worker", "device", "t_view_s", "t_total_s"}
    out.append(sorted(key for key in set(ma) | set(mb) if key not in skip and ma.get(key) != mb.get(key)))
    out.append((ma["host"], ma["device"], mb["host"], mb["device"]))
    return out


def cmd_compare(a):
    import multiprocessing as mp
    ks = [int(v) for v in a.frames.split(",")] if a.frames else list(range(len(load_prep(a.work)["utc"])))
    with mp.get_context("fork").Pool(a.procs) as pool:
        res = pool.map(_cmp_one, [(k, a.a, a.b) for k in ks], chunksize=8)
    bad = [r for r in res if r[1] or r[2] or r[3]]
    pairs = sorted({r[4] for r in res})
    rep = {"a": a.a, "b": a.b, "frames": len(res), "frames_differing": len(bad),
           "pixels_differing_rg": int(sum(r[1] for r in res)), "pixels_differing_ng": int(sum(r[2] for r in res)),
           "meta_keys_differing": sorted({key for r in res for key in r[3]}), "device_pairs": [list(p) for p in pairs],
           "first_differing": bad[:20]}
    if a.report:
        json.dump(rep, open(a.report, "x"), indent=1)
    say("compare: %d frames, %d differ (pixels rg %d, ng %d; meta keys %s); %d (host, device) pairs" % (
        rep["frames"], rep["frames_differing"], rep["pixels_differing_rg"], rep["pixels_differing_ng"],
        rep["meta_keys_differing"], len(pairs)))
    return 0 if not bad else 1


def cmd_touchdown(a):
    """Touchdown (Rosetta impact) point for the movie flag (2026-11-13, user: "mark on the nucleus surface the location
    of the touchdown"; "make the touchdown location a Formula 1 waving finish flag"): cops_cuda/figures.impact_point
    (the code path of tools/pub_figures.py groundtrack: straight-line extrapolation of the last two observation rows of
    the ground-track table --track-base to the shape model, Moller-Trumbore on the facets near the last row) on this
    mesh, compared with the ground-track figures' _impact.json (--impact-json; refused if the point or facet differs);
    per frame: its pixel, the pixel of a short step along the facet normal (pole direction: the projected normal when
    the normal is >= FLAG_NORMAL_MIN_DEG from the line of sight, else image up), and its visibility from the camera:
    in front of the camera, inside the image, facet facing the camera, and the exact first-hit ray camera -> point
    (ray_first_hit over all facets) reaches the point. Check (--kernels, --geom-exe): the SPICE trajectory (RORB:
    ROSETTA, as the ground-track table; and the COPS gauge = ROSETTA + GaugeOffsetSC) intersected with the mesh, its
    distance from the touchdown point. Writes <work>/touchdown.npz and touchdown.json."""
    import cupy as cp
    work = Path(a.work)
    for name in ("touchdown.npz", "touchdown.json"):
        if (work / name).exists():
            raise SystemExit("%s exists" % (work / name))
    P = load_prep(work)
    Q = np.load(work / "plan.npz", allow_pickle=False)
    rnd = Renderer(a.geometry, P, a.device)
    nodes_h = cp.asnumpy(rnd.nodes)
    faces_h = cp.asnumpy(rnd.faces).reshape(-1, 3)
    ttab = FIG.read_table(a.track_base)
    imp = FIG.impact_point(nodes_h, faces_h, ttab)
    ptd = np.asarray(imp["point_m"], float)
    ftd = int(imp["face"])
    ntd = cp.asnumpy(rnd.norm[ftd])
    A_, B_, C_ = nodes_h[faces_h[ftd]]
    nrm_plane = np.cross(B_ - A_, C_ - A_)
    res = {"definition": cmd_touchdown.__doc__, "track_base": str(Path(a.track_base).resolve()), "impact": imp,
           "facet_normal": ntd.tolist(), "point_off_facet_plane_m": float(np.dot(ptd - A_, nrm_plane / np.linalg.norm(nrm_plane))),
           "mesh": str(Path(a.geometry).resolve()), "camera_mode": "fixed" if is_fixed(Q) else "gauge"}
    if a.impact_json:
        ref = json.load(open(a.impact_json))
        cmpd = {"file": str(Path(a.impact_json).resolve()), "ref_shape": ref.get("shape"), "ref_base_run": ref.get("base_run"),
                "point_diff_m": float(np.linalg.norm(ptd - np.asarray(ref["point_m"], float))),
                "same_face": int(ref["face"]) == ftd, "same_utc": ref["utc"] == imp["utc"],
                "lon_diff_deg": imp["lon_deg"] - ref["lon_deg"], "lat_diff_deg": imp["lat_deg"] - ref["lat_deg"]}
        res["vs_ground_track_figures"] = cmpd
        say("touchdown vs %s: %s" % (a.impact_json, json.dumps(cmpd)))
        if cmpd["point_diff_m"] > 1e-6 or not cmpd["same_face"] or not cmpd["same_utc"]:
            raise SystemExit("the touchdown point differs from the ground-track figures' impact point")
    say("touchdown: %s UTC, lon %.3f, lat %.3f, facet %d, point %s" % (imp["utc"], imp["lon_deg"], imp["lat_deg"], ftd,
                                                                        ptd.tolist()))
    N = len(P["utc"])
    base = np.full((N, 2), np.nan)
    tip = np.full((N, 2), np.nan)
    upv = np.tile(np.array([0.0, -1.0]), (N, 1))
    ndef = np.zeros(N, bool)
    front = np.zeros(N, bool)
    inimg = np.zeros(N, bool)
    facing = np.zeros(N, bool)
    dist = np.zeros(N)
    cams = np.zeros((N, 3))
    smin = math.sin(math.radians(FLAG_NORMAL_MIN_DEG))
    for k in range(N):
        xc, brg, bng = frame_camera(P, Q, k)
        cams[k] = xc
        f, up_, right = rnd.camera(xc, brg, bng)
        proj, F, center = str(Q["proj"][k]), float(Q["F"][k]), tuple(float(v) for v in Q["center"][k])
        v = ptd - xc
        vc = (float(v @ right), float(v @ up_), float(v @ f))
        dist[k] = float(np.linalg.norm(v))
        front[k] = vc[2] > ZNEAR if proj == "perspective" else math.atan2(math.hypot(vc[0], vc[1]), vc[2]) < math.radians(170.0)
        facing[k] = float(np.dot(ntd, xc - ptd)) > 0.0
        if not front[k]:
            continue
        px = to_pixel(vc, proj, F, center)
        if px is None:
            front[k] = False
            continue
        base[k] = px
        inimg[k] = 0 <= px[0] < W and 0 <= px[1] < H
        vt = ptd + ntd * (0.02 * dist[k]) - xc
        pt_ = to_pixel((float(vt @ right), float(vt @ up_), float(vt @ f)), proj, F, center)
        sinang = float(np.linalg.norm(np.cross(ntd, v / dist[k])))
        if pt_ is not None:
            tip[k] = pt_
            dd = np.asarray(pt_) - np.asarray(px)
            if sinang >= smin and np.hypot(*dd) > 1e-9:
                upv[k] = dd / np.hypot(*dd)
                ndef[k] = True
    if is_fixed(Q):
        v1, hf1, th1, d1 = rnd.points_visible(cams[:1], ptd[None, :], [ftd])
        rayvis = np.repeat(v1, N)
        hf, th = np.repeat(hf1, N), np.repeat(th1, N)
    else:
        rayvis, hf, th, _ = rnd.points_visible(cams, np.tile(ptd, (N, 1)), np.full(N, ftd))
    visible = front & inimg & facing & rayvis
    res["frames"] = {"visible": int(visible.sum()), "first_visible": int(np.argmax(visible)) if visible.any() else None,
                     "last_visible": int(N - 1 - np.argmax(visible[::-1])) if visible.any() else None,
                     "in_front": int(front.sum()), "in_image": int(inimg.sum()), "facing": int(facing.sum()),
                     "ray_reaches_point": int(rayvis.sum()), "pole_along_projected_normal": int((ndef & visible).sum()),
                     "pole_image_up": int((~ndef & visible).sum())}
    if visible.any():
        ks = np.nonzero(visible)[0]
        res["frames"]["visible_runs"] = [[int(r[0]), int(r[-1])] for r in np.split(ks, np.nonzero(np.diff(ks) > 1)[0] + 1)]
    say("touchdown visibility: %s" % json.dumps(res["frames"]))
    # check: the SPICE trajectory intersected with the mesh
    if a.kernels and a.geom_exe:
        et_t = ttab["et"]
        secs = np.array([(t - DAY0).total_seconds() for t in ttab["time"]])
        off = float(np.median(et_t - secs))
        e0 = float(et_t[-1]) - 4.0
        ets = e0 + 0.05 * np.arange(int(round(40.0 / 0.05)) + 1)
        utcs = [(DAY0 + dt.timedelta(seconds=float(e - off))).strftime("%Y-%m-%dT%H:%M:%S.%f") for e in ets]
        res["spice_check"] = {"times": "%s .. %s UTC, every 0.05 s" % (utcs[0], utcs[-1]), "tracks": []}
        for lab, offs in (("ROSETTA (RORB; the ground-track table's positions)", "0,0,0"),
                          ("COPS gauge (ROSETTA + GaugeOffsetSC %s m)" % a.gauge_offset.replace(",", " "), a.gauge_offset)):
            r = subprocess.run([a.geom_exe, "geom", a.kernels] + ["%g" % float(v) for v in offs.split(",")],
                               input="\n".join(utcs) + "\n", capture_output=True, text=True, check=True)
            Lw = [ln.split() for ln in r.stdout.split("\n") if ln.strip()]
            X = np.array([[float(t) for t in ln[2:5]] for ln in Lw])
            E = np.array([float(ln[1]) for ln in Lw])
            seg = X[1:] - X[:-1]
            sl = np.linalg.norm(seg, axis=1)
            hfs, ths = rnd.first_hits(X[:-1], seg / sl[:, None])
            okh = (hfs >= 0) & (ths <= sl * (1 + 1e-6))
            tr = {"track": lab, "points": len(X)}
            if lab.startswith("ROSETTA"):
                kt = int(np.argmin(np.abs(E - et_t[-1])))
                xt_last = np.array([ttab["x"][-1], ttab["y"][-1], ttab["z"][-1]])
                tr["vs_table_last_row_m"] = float(np.linalg.norm(X[kt] - xt_last))
            if okh.any():
                i = int(np.argmax(okh))
                fh = int(hfs[i])
                n_, c_ = cp.asnumpy(rnd.norm[fh]), cp.asnumpy(rnd.cent[fh])
                dvec = seg[i] / sl[i]
                tt = float(np.dot(c_ - X[i], n_) / np.dot(dvec, n_))
                Xs = X[i] + tt * dvec
                when = E[i] + tt / sl[i] * (E[i + 1] - E[i])
                utc_s = (DAY0 + dt.timedelta(seconds=float(when - off)))
                tr.update({"surface_point_m": Xs.tolist(), "face": fh,
                           "utc": utc_s.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3],
                           "lon_deg": float(np.degrees(np.arctan2(Xs[1], Xs[0]))),
                           "lat_deg": float(np.degrees(np.arcsin(Xs[2] / np.linalg.norm(Xs)))),
                           "distance_from_touchdown_point_m": float(np.linalg.norm(Xs - ptd)),
                           "time_minus_touchdown_utc_s": float(when - (float(et_t[-1]) + imp["seconds_after_last_row"])),
                           "speed_m_s": float(sl[i] / (E[i + 1] - E[i]))})
            else:
                tr["note"] = "no crossing of the surface in this interval"
            res["spice_check"]["tracks"].append(tr)
            say("SPICE check %s: %s" % (lab, json.dumps(tr)))
    np.savez(work / "touchdown.npz", utc=np.array(imp["utc"]), point=ptd, face=ftd, normal=ntd,
             lon_deg=imp["lon_deg"], lat_deg=imp["lat_deg"], base_px=base, tip_px=tip, up_px=upv, normal_defined=ndef,
             visible=visible, in_front=front, in_image=inimg, facing=facing, ray_visible=rayvis, ray_first_face=hf,
             ray_first_t=th, distance_m=dist, camera=cams)
    json.dump(res, open(work / "touchdown.json", "x"), indent=1, default=float)
    say("touchdown done")
    return 0


def cmd_validate_fixed(a):
    """Fixed-camera checks before the render (2026-11-13): (1) the camera (plan.npz: SPICE vs the base table, its altitude
    above the mesh) and its view, computed twice (identical face buffer); (2) per test frame and gauge: the FOV condition
    of every camera nucleus pixel against an independent instrument-frame rectangle test (getfov boresight and
    half-angles in the instrument frame, the frame's rotation from fov.npz), as validate_fov; the line of sight (the
    render's result = depth map + exact edge refinement, and the depth map alone) against exact rays (ray_first_hit
    over ALL facets from the point + 1e-6 MaxEdge along the normal toward the gauge; blocked when a facet other than the
    point's is hit before the gauge) at --samples random points and --samples points on the boundary of the seen
    region (in FOV and facing; of the depth-map result and of the refined result); (3) test frames (as render) and
    diagnostic images (grey: outside the FOV; blue: inside, facet facing away from the gauge; red: inside, facing,
    no line of sight; tint: seen; yellow squares: samples where the render's result and the exact ray disagree)."""
    import cupy as cp
    from PIL import Image, ImageDraw
    work = Path(a.work)
    P = load_prep(work)
    Q = np.load(work / "plan.npz", allow_pickle=False)
    if not (is_fixed(Q) or is_detach(Q)):
        raise SystemExit("validate_fixed: plan.npz is not a fixed-camera or detach plan")
    FV = load_fov(work)
    if FV is None:
        raise SystemExit("validate_fixed needs fov.npz")
    SCZ = np.load(work / "sc.npz", allow_pickle=False) if (work / "sc.npz").exists() else None
    rnd = Renderer(a.geometry, P, a.device)
    info = json.load(open(work / "fov_checks.json"))
    specs = {"RG": a.mask_rg, "NG": a.mask_ng}
    fixed_check_masks(specs)
    if is_detach(Q):                           # detach-200 m set: the frozen camera of the detached frames
        k_d = int(Q["k_detach"])
        xc = np.asarray(Q["cam_x_frames"][k_d], float)
        F = float(Q["F"][k_d])
        cbr, cbn = Q["cam_b_rg_frames"][k_d], Q["cam_b_ng_frames"][k_d]
        parts = ("bus", "wings", "hga")
    else:
        xc = np.asarray(Q["cam_x"], float)
        F = float(Q["F"][0])
        cbr, cbn = Q["cam_b_rg"], Q["cam_b_ng"]
        parts = ("bus", "wings")
    t0 = time.time()
    vw0 = rnd.fixed_geometry(xc, cbr, cbn, "perspective", F, (0.0, 0.0))
    t_view = time.time() - t0
    vw1 = rnd.fixed_geometry(xc, cbr, cbn, "perspective", F, (0.0, 0.0))
    fdraw = cp.unique(vw0["fid"])
    adraw = cp.asarray(rnd.area_h)[fdraw]
    cen_d = rnd.cent[fdraw]
    res = {"definition": cmd_validate_fixed.__doc__, "camera_checks": json.loads(str(Q["cam_checks"])),
           "view": {"nucleus_pixels_ss": int(vw0["idx"].size), "nucleus_pixel_frac": int(vw0["idx"].size) / (W * H * SS * SS),
                    "facets": int(fdraw.size), "facets_area_m2": float(cp.asnumpy(adraw.sum())), "holes": vw0["holes"],
                    "face_sha256": vw0["face_sha256"], "repeat_face_sha256_equal": vw0["face_sha256"] == vw1["face_sha256"],
                    "seconds": t_view}, "frames": []}
    del vw1
    say("view: %s" % json.dumps(res["view"]))
    me = rnd.maxedge()
    rng = np.random.default_rng(20261113)
    od = work / "test_frames_fixed"
    od.mkdir(exist_ok=True)
    idx, pts, fid, nrm = vw0["idx"], vw0["pts"], vw0["fid"], vw0["nrm"]
    Hs, Ws = H * SS, W * SS
    pan = is_detach(Q) and "pan_s" in Q.files  # 2026-11-14: moving camera: each test frame with its own pose
    vkey = (xc.tobytes(), np.asarray(cbr).tobytes(), np.asarray(cbn).tobytes(), F)
    for k in [int(v) for v in a.test_frames.split(",")]:
        if pan and k >= int(Q["k_detach"]):
            xk, brk, bnk = frame_camera(P, Q, k)
            key = (xk.tobytes(), np.asarray(brk).tobytes(), np.asarray(bnk).tobytes(), float(Q["F"][k]))
            if key != vkey:
                vw0 = None
                vw0 = rnd.fixed_geometry(xk, brk, bnk, "perspective", float(Q["F"][k]), (0.0, 0.0))
                xc, vkey = xk, key
                fdraw = cp.unique(vw0["fid"])
                adraw = cp.asarray(rnd.area_h)[fdraw]
                cen_d = rnd.cent[fdraw]
                idx, pts, fid, nrm = vw0["idx"], vw0["pts"], vw0["fid"], vw0["nrm"]
        t0 = time.time()
        R = fixed_frame(rnd, P, FV, SCZ, specs, k, vw0, cen_d, fdraw, adraw, a.sc_samples, sc_parts=parts)
        t_fr = time.time() - t0
        xg = np.asarray(P["x"][k], float)
        l = cp.asarray(xg)[None, :] - pts
        l /= cp.linalg.norm(l, axis=1)[:, None]
        facing = cp.sum(nrm * l, axis=1) > 0.0
        dB = cp.asnumpy(-l)
        ent = {"frame": k, "utc": str(P["utc"][k]), "altitude_m": float(P["alt"][k]), "seconds_frame": t_fr,
               "gauge_distance_from_camera_m": float(np.linalg.norm(xg - xc))}
        for which, g in (("rg", "RG"), ("ng", "NG")):
            G = R["g"][g]
            gi = info["info"][g]
            bsI = np.array(gi["boresight"])
            e1, e2 = (np.array(v) for v in gi["axes"])
            c0 = np.array(gi["bounds"][0])
            ha = (np.degrees(np.arctan2(abs(c0.dot(e1)), c0.dot(bsI))), np.degrees(np.arctan2(abs(c0.dot(e2)), c0.dot(bsI))))
            Rm = FV["frame_rot_%s" % which][k]
            dI = dB @ Rm
            fwd = dI @ bsI
            aa = np.degrees(np.arctan2(dI @ e1, fwd))
            bb = np.degrees(np.arctan2(dI @ e2, fwd))
            expect = (fwd > 0) & (np.abs(aa) <= ha[0]) & (np.abs(bb) <= ha[1])
            near = (np.abs(np.abs(aa) - ha[0]) < 1e-9) | (np.abs(np.abs(bb) - ha[1]) < 1e-9)
            if parse_mask(specs[g])[0] == "cone":   # circular cone: the angle from the boresight (prep), in degrees
                bP = unit(P["b_rg"][k] if g == "RG" else P["b_ng"][k])
                angd = np.degrees(np.arccos(np.clip(dB @ bP, -1.0, 1.0)))
                expect = angd <= parse_mask(specs[g])[1]
                near = np.abs(angd - parse_mask(specs[g])[1]) < 1e-6
            got = cp.asnumpy(G["fov"])
            sel = G["fov"] & facing
            seen = G["seen"]
            selh, seenh, seenmh = cp.asnumpy(sel), cp.asnumpy(seen), cp.asnumpy(G["seen_map"])
            st = {"fov_analytic_mismatches": int(((got != expect) & ~near).sum()),
                  "fov_mismatches_within_1e-9_deg_of_an_edge": int(((got != expect) & near).sum()),
                  "fov_pixels": int(got.sum()), "fov_facing_pixels": int(selh.sum()), "seen_pixels": int(seenh.sum()),
                  "seen_pixels_map_only": int(seenmh.sum()),
                  "half_angles_deg": [float(ha[0]), float(ha[1])], "los": G["stats"]["los"]}
            # line of sight vs exact rays: random samples, and samples on the boundary of the seen region within sel
            # (of the depth-map result and of the refined result); agreement of both results with the exact rays
            full_sel = _full(rnd, idx, sel, bool)
            full_seen = _full(rnd, idx, seen, bool)
            pos = cp.full(Hs * Ws, -1, dtype=cp.int64)
            pos[idx] = cp.arange(int(idx.size))
            pos = pos.reshape(Hs, Ws)

            def boundary_of(fs):
                bnd = cp.zeros((Hs, Ws), dtype=bool)
                for sl_a, sl_b in (((slice(None), slice(1, None)), (slice(None), slice(None, -1))),
                                   ((slice(1, None), slice(None)), (slice(None, -1), slice(None)))):
                    dif = full_sel[sl_a] & full_sel[sl_b] & (fs[sl_a] != fs[sl_b])
                    bnd[sl_a] |= dif
                    bnd[sl_b] |= dif
                return cp.asnumpy(pos[bnd])

            cand_r = np.nonzero(selh)[0]
            cand_bm = boundary_of(_full(rnd, idx, G["seen_map"], bool))
            cand_b = boundary_of(full_seen)
            samp = {}
            for lab, cand in (("random", cand_r), ("boundary_of_map_result", cand_bm), ("boundary", cand_b)):
                if len(cand) == 0:
                    samp[lab] = {"n": 0}
                    continue
                pick = rng.choice(cand, size=min(a.samples, len(cand)), replace=False)
                pk = cp.asarray(pick)
                fq = fid[pk]
                o = cp.asnumpy(pts[pk] + nrm[pk] * (1e-6 * me[fq])[:, None])
                dv = xg[None, :] - o
                dd = np.linalg.norm(dv, axis=1)
                hf, th = rnd.first_hits(o, dv / dd[:, None])
                blocked = (hf >= 0) & (th < dd) & (hf != cp.asnumpy(fq))
                mapvis = seenh[pick]
                maponly = seenmh[pick]
                samp[lab] = {"n": int(len(pick)), "candidates": int(len(cand)),
                             "agree": int((mapvis == ~blocked).sum()),
                             "map_visible_exact_blocked": int((mapvis & blocked).sum()),
                             "map_blocked_exact_visible": int((~mapvis & ~blocked).sum()),
                             "exact_blocked": int(blocked.sum()),
                             "map_only_agree": int((maponly == ~blocked).sum()),
                             "map_only_visible_exact_blocked": int((maponly & blocked).sum()),
                             "map_only_blocked_exact_visible": int((~maponly & ~blocked).sum())}
                samp[lab]["agree_frac"] = samp[lab]["agree"] / max(samp[lab]["n"], 1)
                samp[lab]["map_only_agree_frac"] = samp[lab]["map_only_agree"] / max(samp[lab]["n"], 1)
                samp[lab]["_bad_idx"] = pick[mapvis != ~blocked]
            st["los_vs_exact_rays"] = {kk: {k2: v2 for k2, v2 in vv.items() if k2 != "_bad_idx"} for kk, vv in samp.items()}
            ent[g] = st
            say("frame %d %s: %s" % (k, g, json.dumps({kk: vv for kk, vv in st.items() if kk != "los"})))
            # test frame (as render) and the diagnostic image
            img = rnd.colorize_fixed(R["vw"], which, _full(rnd, idx, G["w"], cp.float64), full_seen,
                                     _full(rnd, idx, G["fov"], bool), "weight", None, True, sunfac=R["sunfac"])
            Image.fromarray(img).save(od / ("f%05d_%s.png" % (k, which)), compress_level=1)
            cls = cp.zeros(int(idx.size), dtype=cp.int8)
            cls[G["fov"] & ~facing] = 1
            cls[sel & ~seen] = 2
            cls[seen] = 3
            vwc = R["vw"]
            bright = AMBIENT + DIFFUSE * cp.maximum(vwc["cos"], 0.0) * vwc["lit"]
            bpx = bright[fid]
            wmax = 1.0 if which == "rg" else 2.0
            kk_ = cp.clip((G["w"] / wmax * 255.0 + 0.5).astype(cp.int64), 0, 255)
            col = cp.broadcast_to(cp.asarray(GREY_ZERO)[None, :], (int(idx.size), 3)).copy()
            col[cls == 1] = cp.asarray([0.25, 0.45, 1.0])
            col[cls == 2] = cp.asarray([1.0, 0.15, 0.15])
            col[cls == 3] = rnd.lut[kk_[cls == 3]]
            imgd = cp.zeros((Hs * Ws, 3))
            imgd[idx] = col * bpx[:, None]
            imgd = imgd.reshape(H, SS, W, SS, 3).mean(axis=(1, 3))
            im = Image.fromarray(cp.asnumpy(cp.clip(imgd * 255.0 + 0.5, 0, 255).astype(cp.uint8)))
            dr = ImageDraw.Draw(im)
            for lab in samp:
                for p_ in samp[lab].get("_bad_idx", []):
                    fl = int(cp.asnumpy(idx[int(p_)]))
                    yy, xx = (fl // Ws) / SS, (fl % Ws) / SS
                    dr.rectangle((xx - 3, yy - 3, xx + 3, yy + 3), outline=(255, 230, 0), width=2)
            fnt = font(22)
            lines = ["frame %d, %s UTC, COPS %.1f m: %s line of sight from the true gauge position" % (
                         k, str(P["utc"][k])[11:19], float(P["alt"][k]), g),
                     "grey: outside the %s IK FOV; blue: inside, facet turned away from %s; red: inside, facing, hidden"
                     % (g, g), "tint: seen (%s weight); yellow squares: samples where the map test and the exact ray "
                               "disagree" % g]
            y = 12
            for t_ in lines:
                bbx = dr.textbbox((12, y), t_, font=fnt)
                dr.rectangle((bbx[0] - 4, bbx[1] - 3, bbx[2] + 4, bbx[3] + 3), fill=(0, 0, 0))
                dr.text((12, y), t_, font=fnt, fill=(255, 255, 255))
                y = bbx[3] + 8
            im.save(od / ("f%05d_%s_los_diagnostic.png" % (k, which)))
            del full_sel, full_seen, pos, cls, col, imgd
        res["frames"].append(ent)
        del R
    json.dump(res, open(work / ("validation_fixed_%s.json" % socket.gethostname()), "w"), indent=1, default=float)
    bad = sum(e[g]["fov_analytic_mismatches"] for e in res["frames"] for g in ("RG", "NG"))
    say("validate_fixed done: %d FOV mismatches; view repeat identical %s" % (bad, res["view"]["repeat_face_sha256_equal"]))
    return 0 if (bad == 0 and res["view"]["repeat_face_sha256_equal"]) else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = argparse.ArgumentParser(add_help=False)
    c.add_argument("--work", required=True)
    c.add_argument("--geometry", default=None)
    c.add_argument("--base", default=None)
    c.add_argument("--device", type=int, default=0)
    p = sub.add_parser("prep", parents=[c])
    p.add_argument("--kernels", required=True)
    p.add_argument("--geom-exe", required=True)
    p.add_argument("--gauge-offset", default="0,0,1.5")
    p.add_argument("--extend-to-utc", default=None, help="2026-11-13: add frames uniform in time (the last log-altitude "
                   "step) after the last observation up to this UTC (default: none)")
    p = sub.add_parser("plan", parents=[c])
    p.add_argument("--camera", default="gauge", choices=("gauge", "fixed", "detach"),
                   help="gauge (default: the camera at the gauge, the earlier movies), fixed or detach (2026-11-13)")
    p.add_argument("--detach-altitude", type=float, default=200.0, help="detach camera: the COPS altitude [m] at which "
                   "the camera leaves the gauge and stays")
    p.add_argument("--camera-altitude", type=float, default=100.0, help="fixed camera: the COPS altitude [m] of its time")
    p.add_argument("--fov-deg", type=float, default=90.0, help="fixed camera: horizontal field of view [deg]")
    p.add_argument("--pan-to", default="", help="2026-11-14: detach camera: move to this final pose after the detach "
                   "(json: camera_position_m, forward, up, fov_horizontal_deg, pivot_m; default: frozen camera)")
    p.add_argument("--pan-start-frame", type=int, default=-1, help="pan: first frame of the motion (default the detach frame)")
    p.add_argument("--pan-end-frame", type=int, default=-1, help="pan: frame at which the final pose is reached (default last)")
    p.add_argument("--kernels", default=None)
    p.add_argument("--geom-exe", default=None)
    p.add_argument("--gauge-offset", default="0,0,1.5")
    p = sub.add_parser("validate", parents=[c])
    p.add_argument("--rows", default="0,300,600,900,1079,1300,1450,1549")
    p.add_argument("--test-frames", default="0,30,600,687,688,1200,1500,1700,1799")
    p = sub.add_parser("render", parents=[c])
    p.add_argument("--worker", type=int, default=0)
    p.add_argument("--nworkers", type=int, default=1)
    p.add_argument("--frames-dir", default="", help="output folder (default <work>/frames); a check re-render goes elsewhere")
    p.add_argument("--only", default="", help="comma list of frames (default all)")
    for g in ("rg", "ng"):
        p.add_argument("--mask-%s" % g, default="none", help="observable part of %s: none (default: the first movies), "
                       "fit, hemisphere, ik_fov, cone:<half-angle deg>" % g.upper())
        p.add_argument("--tint-%s" % g, default="weight", choices=("weight", "uniform"))
    p.add_argument("--uniform-color", default="#E69F00", help="tint 'uniform': its colour")
    p.add_argument("--no-boundary", action="store_true", help="masked frames: do not draw the mask boundary line")
    p.add_argument("--sc-shadow", action="store_true", help="draw the spacecraft shadow (bus + wings; needs sc.npz)")
    p.add_argument("--sc-samples", type=int, default=SC_SAMPLES, help="points on the solar disc")
    p = sub.add_parser("fov", parents=[c])
    p.add_argument("--kernels", required=True)
    p.add_argument("--geom-exe", required=True)
    p = sub.add_parser("sc", parents=[c])
    p.add_argument("--kernels", required=True)
    p.add_argument("--geom-exe", required=True)
    p.add_argument("--sa-hold-utc", default=SA_HOLD_UTC)
    p.add_argument("--origin-offset", default="0,0,0", help="position of ROS_SPACECRAFT's origin (-226000) relative to "
                   "ROSETTA (-226) in ROS_SPACECRAFT [m]; no kernel gives it (default 0)")
    p.add_argument("--hga-ck", default=None, help="2026-11-13: HGA C-kernel -> HGA and solar-array frames in sc.npz")
    p.add_argument("--hga-hold-utc", default=SA_HOLD_UTC, help="HGA orientation held after this UTC (end of the C-kernel)")
    p = sub.add_parser("contact", parents=[c])
    p.add_argument("--kernels", required=True)
    p.add_argument("--geom-exe", required=True)
    p.add_argument("--hga-ck", required=True)
    p.add_argument("--hga-hold-utc", default=SA_HOLD_UTC)
    p.add_argument("--sa-hold-utc", default=SA_HOLD_UTC)
    p.add_argument("--t-start", default="2016-09-30T10:39:26")
    p.add_argument("--max-seconds", type=float, default=15.0)
    p = sub.add_parser("validate_sc", parents=[c])
    p.add_argument("--test-frames", default="1200,1500,1600,1700,1799")
    p = sub.add_parser("validate_fov", parents=[c])
    p.add_argument("--rows", default="0,300,600,900,1079,1300,1450,1549")
    p.add_argument("--test-frames", default="0,300,687,688,904,1198,1500,1799")
    p = sub.add_parser("diagnose", parents=[c])
    p.add_argument("--frames", default="1500,1600,1700,1799")
    p.add_argument("--frames-dir", required=True, help="the rendered clean frames (rg_clean/) to annotate")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--box", action="append", help="frame,x0,x1,y0,y1 (output pixels): a region to list facet by facet")
    p = sub.add_parser("frametable", parents=[c])
    p.add_argument("--frames-dir", default="")
    p.add_argument("--out", required=True)
    p = sub.add_parser("compose", parents=[c])
    p.add_argument("--out", required=True)
    p.add_argument("--table", required=True)
    p.add_argument("--ts-label", default="", help="text right of the movie-3 time series ('\\n' = new line; default: the "
                   "first movies' label)")
    p.add_argument("--procs", type=int, default=16)
    p.add_argument("--frames", default="")
    p.add_argument("--touchdown", default="", help="touchdown.npz (cops_movies.py touchdown): draw the touchdown flag "
                   "(2026-11-13; default: none)")
    p.add_argument("--no-sub-mark", action="store_true", help="2026-11-14: do not draw the sub-spacecraft point (+) "
                   "(default: drawn)")
    p.add_argument("--flag-on-top", action="store_true", help="2026-11-13: draw the flag on top (no depth test against "
                   "the nucleus or the spacecraft) wherever its point projects into the image")
    p = sub.add_parser("compare", parents=[c])
    p.add_argument("--a", required=True)
    p.add_argument("--b", required=True)
    p.add_argument("--frames", default="")
    p.add_argument("--procs", type=int, default=16)
    p.add_argument("--report", default="")
    p = sub.add_parser("touchdown", parents=[c])
    p.add_argument("--track-base", required=True, help="the ground-track figures' base (tools/pub_figures.py groundtrack --base)")
    p.add_argument("--impact-json", default="", help="their <out>_impact.json (the point must agree)")
    p.add_argument("--kernels", default=None)
    p.add_argument("--geom-exe", default=None)
    p.add_argument("--gauge-offset", default="0,0,1.5")
    p = sub.add_parser("validate_fixed", parents=[c])
    p.add_argument("--test-frames", default="0,300,687,688,904,1198,1500,1600,1700,1799")
    p.add_argument("--mask-rg", default="ik_fov")
    p.add_argument("--mask-ng", default="ik_fov")
    p.add_argument("--samples", type=int, default=2000, help="exact-ray samples per class (random, boundary)")
    p.add_argument("--sc-samples", type=int, default=SC_SAMPLES)
    p = sub.add_parser("encode", parents=[c])
    p.add_argument("--out", required=True)
    p.add_argument("--frames-dir", required=True)
    p.add_argument("--ffmpeg", required=True)
    a = ap.parse_args()
    return {"prep": cmd_prep, "plan": cmd_plan, "validate": cmd_validate, "render": cmd_render, "compose": cmd_compose,
            "compare": cmd_compare, "encode": cmd_encode, "fov": cmd_fov, "validate_fov": cmd_validate_fov,
            "frametable": cmd_frametable, "diagnose": cmd_diagnose, "sc": cmd_sc, "validate_sc": cmd_validate_sc,
            "touchdown": cmd_touchdown, "validate_fixed": cmd_validate_fixed, "contact": cmd_contact}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
