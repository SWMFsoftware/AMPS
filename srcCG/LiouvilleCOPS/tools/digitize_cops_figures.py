#!/usr/bin/env python3
"""Digitize the published end-of-mission (2016-09-30) COPS pressure figures.

DATA SOURCES (image files in <COPS_DATA_ROOT>/sources/)
  1. EGU2019_Tenishev_slide02.jpg
     Origin   V. Tenishev, M. Rubin, C.-Y. Tzou, X. Shi, M. R. Combi, K. Altwegg, T. Gombosi, Y. Shou, Z. Huang,
              G. Toth, H. Sierks, M. Hofstadter, "End-of-Mission ROSINA/COPS observation of the innermost coma of
              comet 67P/Churyumov-Gerasimenko", EGU General Assembly 2019, slide 2 ("End of mission COPS
              measurements ... profiles 'flip' at altitude of about 360 m"). Extracted from ~/Downloads/2019-EGU.pptx
              (ppt/media/image2.jpg, 790x703 px).
     Content  Tecplot, log-log. x: "Altitude" [m] (definition not given), 10^1.24 .. ~10^4.3; y: "Pressure" (unit not
              given; the values are the mbar values of the COPS archive). NG: blue dashed ("Nude Gauge COPS
              Measuremetns"), RG: red dashed ("Ram Gauge COPS Measurements").
  2. EGU2019_Tenishev_slide16_data.jpg
     Origin   the same presentation, slide 16 "Data/Model comparison", left figure "Observed RG and NG pressure"
              (ppt/media/image26.jpg, 789x703 px).
     Content  Tecplot, log-log, x: "Altitude" [m] 10^0.65 .. ~10^4.25, y: "Pressure" (mbar values), NG blue dashed
              ("Nude Gauge Pressure (data)"), RG red dashed ("Ram Gauge Pressure (data)").
              The NG/RG crossing is at ~700 m (slide 2: ~360 m).
     Slide 2 vs. slide 16  The two slides show different images (different files, md5 8bffb68d.. and eba25ad9..),
              not the same plot. EGU2019_slide02_on_slide16_axes.png (--check) shows slide 2 resampled onto the axes
              of slide 16: above ~4 km the curves coincide pixel by pixel, including the RG and NG structure at
              10-16 km (this also checks both axis calibrations); below ~3 km they differ:
                altitude [m]        20     100    300    700    1000   2000
                NG slide 16/slide 2 1.7    1.9    1.7    1.3    1.25   1.18
                RG slide 16/slide 2  -     1.08   1.10   1.10   1.12   1.15
              Compared with the archive L2 data vs. the altitude above M004, slide 16 agrees (RG within 1-3% at all
              altitudes, NG above 700 m; below 300 m it follows the ESA blog curve). Slide 2 at 500-2000 m has the
              NG and RG of the archive records at 370-520 m larger M004 altitude (e.g. slide 2 at 1000 m = archive
              10:09 UTC, 1386 m); below ~300 m the slide 2 NG is one straight segment from 2.3e-8 mbar (~300 m) to
              8e-8 mbar (18 m). The altitude used in slide 2 is not stated in the presentation.
  3. ESA_blog_2016-12-15_COPS_EoM.jpg
     Origin   ESA Rosetta blog, "Rosetta's last words: science descending to a comet", 2016-12-15,
              https://blogs.esa.int/rosetta/2016/12/15/rosettas-last-words-science-descending-to-a-comet/
              image https://blogs.esa.int/rosetta/files/2016/12/COPS_EoM.jpg (3216x2461 px), credit K. Altwegg.
              Caption: "ROSINA's Comet Pressure Sensor (COPS) measured the pressure of gas in the vicinity of the
              spacecraft during its final hours ... Small breaks in the black line are due to internal
              adjustments of the COPS instrument made during the descent to cope with the increasing pressure."
     Content  Origin plot. x: "UTC, 30. Sept. 2016" 00:00 .. 12:00; left y: "pressure (mbar)", log, 10^-10 .. ~2e-7
              (values above the frame are clipped); "COPS pressure": black dots (the nude gauge; it matches the NG
              of the archive); right y (red): "distance to comet center (km)" (red line, not digitized).
     Note     The figure was made from the data of 2016 (preliminary processing). It contains NG values from
              10:24:42 to 10:33:14, where the archived L2 V2.0 data have none, and it is ~230 s earlier than the
              archived NG (06:00-10:18; see tools/collect_cops_sources.py).

AXIS CALIBRATION
  Linear in log10(pressure), log10(altitude) or UTC between two major tick marks; the tick pixel positions were
  measured on the images (dark tick marks next to the axis lines):
    slide 2   x: 10^2 at 247.5 px, 10^4 at 628.5 px;   y: 1e-7 at 88.5 px, 1e-9 at 510.5 px
    slide 16  x: 10^1 at 161.5 px, 10^4 at 657.5 px;   y: 1e-7 at 132.5 px, 1e-10 at 558.5 px
    blog      x: 00:00 at 576 px, 12:00 at 2769 px;     y: 1e-7 at 446 px, 1e-10 at 2052 px

METHOD
  EGU figures: pixels of the curve colour (blue: B>150, R,G<120; red: R>120, G,B<90) inside the frame, legend
    masked; for each pixel column the median row of the curve pixels. The dashes leave gaps.
  ESA blog: black pixels (max channel < 110, grey), legend and the inward ticks of the top axis masked;
    centroids of the isolated dots, and for merged dots (a dot is ~12 px = 4 min wide; the 1-min data overlap)
    the centre line: median row of each pixel column and, after 10:18 (nearly vertical curve), median column of
    each pixel row.
  Resolution in pressure: 0.005 (slide 2), 0.007 (slide 16), 0.002 (blog) decades per pixel. The blog times are
    uncertain by ~ +-2 min (dot width) where the curve is steep.

OUTPUT (<AMPS>/cops_output/data/sources/, or --out-dir under COPS_OUTPUT_ROOT)
  digitized_EGU2019_Tenishev_slide02_{NG,RG}.csv, digitized_EGU2019_Tenishev_slide16_data_{NG,RG}.csv:
      altitude_m, pressure_mbar
  digitized_ESA_blog_2016-12-15_COPS_EoM_NG.csv: utc, pressure_mbar
  digitized_check.png (--check): the digitized points over the images
  EGU2019_slide02_on_slide16_axes.png (--check): slide 2 resampled onto the axes of slide 16 with the two
      calibrations above (slide 16: NG blue, RG red; slide 2: NG cyan, RG orange)

usage: digitize_cops_figures.py [--check] [--dir IMAGE_DIR] [--out-dir OUTPUT_DIR]
"""
import argparse
import datetime as dt
import os

from cops_paths import output_path, sources_root

import numpy as np
import matplotlib

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--check", action="store_true")
ap.add_argument("--dir", default=str(sources_root()),
                help="source-image directory (default: <COPS_DATA_ROOT>/sources; override with COPS_SOURCES_ROOT)")
ap.add_argument("--out-dir", default="data/sources", help="CSV/check-plot directory under cops_output (or absolute path outside LiouvilleCOPS); --dir remains image input only")
args = ap.parse_args()
args.out_dir = str(output_path(args.out_dir))
os.makedirs(args.out_dir, exist_ok=True)
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import ndimage

# pixel -> data calibration: (pixel of decade a, log10 a), (pixel of decade b, log10 b)
EGU = {
    "EGU2019_Tenishev_slide02.jpg": {
        "x": ((247.5, 2.0), (628.5, 4.0)), "y": ((88.5, -7.0), (510.5, -9.0)),
        "frame": (104, 700, 0, 624),                     # x0, x1, y0, y1 (interior; no right frame line)
        "legend": (110, 225, 500, 565),                   # x0, x1, y0, y1 (masked)
    },
    "EGU2019_Tenishev_slide16_data.jpg": {
        "x": ((161.5, 1.0), (657.5, 4.0)), "y": ((132.5, -7.0), (558.5, -10.0)),
        "frame": (104, 700, 0, 624),
        "legend": (160, 280, 425, 480),
    },
}
BLOG = {"file": "ESA_blog_2016-12-15_COPS_EoM.jpg",
        "x": ((576.0, 0.0), (2769.0, 12.0)),               # UTC hours
        "y": ((446.0, -7.0), (2052.0, -10.0)),
        "frame": (582, 2600, 291, 2046),                  # data end at ~10:40 (x=2520); right-axis ticks beyond
        "legend": (1150, 2150, 330, 640),
        "top_ticks": (2480, 2560, 325)}


def lin(cal, p):
    (p0, v0), (p1, v1) = cal
    return v0 + (np.asarray(p, float) - p0) * (v1 - v0) / (p1 - p0)


def digitize_dashed(fname, cal):
    im = plt.imread(os.path.join(args.dir, fname)).astype(float)
    R, G, B = im[..., 0], im[..., 1], im[..., 2]
    masks = {"NG": (B > 150) & (R < 120) & (G < 120),          # blue
             "RG": (R > 120) & (G < 90) & (B < 90)}            # red
    x0, x1, y0, y1 = cal["frame"]
    lx0, lx1, ly0, ly1 = cal["legend"]
    out = {}
    for gauge, m in masks.items():
        m = m.copy()
        m[:, :x0] = m[:, x1:] = False
        m[:y0, :] = m[y1:, :] = False
        m[ly0:ly1, lx0:lx1] = False
        xs, ys = [], []
        for x in range(x0, x1):
            rows = np.nonzero(m[:, x])[0]
            if len(rows):
                xs.append(x)
                ys.append(np.median(rows))
        alt = 10 ** lin(cal["x"], xs)
        p = 10 ** lin(cal["y"], ys)
        out[gauge] = (np.array(xs), np.array(ys), alt, p)
    return out


def digitize_blog():
    im = plt.imread(os.path.join(args.dir, BLOG["file"])).astype(float)
    mx, mn = im.max(axis=2), im.min(axis=2)
    m = (mx < 110) & (mx - mn < 45)                          # black (not the red distance curve)
    x0, x1, y0, y1 = BLOG["frame"]
    lx0, lx1, ly0, ly1 = BLOG["legend"]
    m[:, :x0] = m[:, x1:] = False
    m[:y0, :] = m[y1:, :] = False
    m[ly0:ly1, lx0:lx1] = False
    tx0, tx1, ty1 = BLOG["top_ticks"]                        # inward ticks of the top axis (data only at the end)
    m[:ty1, :tx0] = m[:ty1, tx1:] = False
    lab, n = ndimage.label(m)
    cy, cx = np.array(ndimage.center_of_mass(m, lab, range(1, n + 1))).T
    size = ndimage.sum(m, lab, range(1, n + 1))
    keep = (size >= 20) & (size <= 400)                      # one dot ~ 12x12 px; larger regions = merged dots
    big = size > 400
    # merged dots (a dot is ~12 px = 4 min wide; 1-min data overlap): centre line of the merged region,
    # median row per pixel column, and in the steep part at the end also median column per pixel row
    xs, ys = list(cx[keep]), list(cy[keep])
    xsteep = lin(((0.0, BLOG["x"][0][0]), (12.0, BLOG["x"][1][0])), 10.3)    # pixel column of 10:18 UTC
    for k in np.nonzero(big)[0]:
        sl = ndimage.find_objects(lab == k + 1)[0]
        sub = (lab[sl] == k + 1)
        for j in range(0, sub.shape[1], 3):
            r = np.nonzero(sub[:, j])[0]
            if len(r):
                xs.append(sl[1].start + j)
                ys.append(sl[0].start + np.median(r))
        for i in range(0, sub.shape[0], 3):
            c = np.nonzero(sub[i, :])[0]
            if len(c) and sl[1].start + np.median(c) > xsteep:
                xs.append(sl[1].start + np.median(c))
                ys.append(sl[0].start + i)
    o = np.argsort(xs)
    xs, ys = np.array(xs)[o], np.array(ys)[o]
    hours = lin(BLOG["x"], xs)
    p = 10 ** lin(BLOG["y"], ys)
    return xs, ys, hours, p


results = {}
for fname, cal in EGU.items():
    d = digitize_dashed(fname, cal)
    base = os.path.splitext(fname)[0]
    for gauge, (xs, ys, alt, p) in d.items():
        out = os.path.join(args.out_dir, "digitized_%s_%s.csv" % (base, gauge))
        with open(out, "w") as f:
            f.write("# digitized from %s (%s, dashed line; median curve row per pixel column) by tools/digitize_cops_figures.py\n" % (fname, gauge))
            f.write("# altitude as plotted in the figure (definition not given in the figure)\n")
            f.write("altitude_m,pressure_mbar\n")
            for a, v in zip(alt, p):
                f.write("%.4g,%.5e\n" % (a, v))
        print("%s: %d points" % (out, len(alt)))
        results[(fname, gauge)] = (xs, ys)

xs, ys, hours, p = digitize_blog()
out = os.path.join(args.out_dir, "digitized_ESA_blog_2016-12-15_COPS_EoM_NG.csv")
t0 = dt.datetime(2016, 9, 30)
with open(out, "w") as f:
    f.write("# digitized from %s (COPS pressure, black dots; nude gauge) by tools/digitize_cops_figures.py\n" % BLOG["file"])
    f.write("utc,pressure_mbar\n")
    for h, v in zip(hours, p):
        f.write("%s,%.5e\n" % ((t0 + dt.timedelta(hours=float(h))).strftime("%Y-%m-%dT%H:%M:%S"), v))
print("%s: %d points" % (out, len(hours)))

if args.check:
    fig, ax = plt.subplots(1, 3, figsize=(20, 7))
    for a, (fname, cal) in zip(ax[:2], EGU.items()):
        a.imshow(plt.imread(os.path.join(args.dir, fname)))
        for gauge, col in (("NG", "c"), ("RG", "y")):
            x, y = results[(fname, gauge)]
            a.plot(x, y, ".", ms=2, color=col)
        a.set_title(fname, fontsize=9)
    ax[2].imshow(plt.imread(os.path.join(args.dir, BLOG["file"])))
    ax[2].plot(xs, ys, ".", ms=2, color="c")
    for a in ax:
        a.axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(args.out_dir, "digitized_check.png"), dpi=90)
    print("check plot: %s" % os.path.join(args.out_dir, "digitized_check.png"))

    # slide 2 resampled onto the axes of slide 16 (nearest pixel, through the two calibrations)
    (f2, c2), (f16, c16) = EGU.items()
    i2 = plt.imread(os.path.join(args.dir, f2)).astype(float)
    i16 = plt.imread(os.path.join(args.dir, f16)).astype(float)
    H, W = i16.shape[:2]
    py, px = np.mgrid[0:H, 0:W]
    inv = lambda cal, v: lin(((cal[0][1], cal[0][0]), (cal[1][1], cal[1][0])), v)     # data -> pixel
    qx = np.rint(inv(c2["x"], lin(c16["x"], px))).astype(int)
    qy = np.rint(inv(c2["y"], lin(c16["y"], py))).astype(int)
    x0, x1, y0, y1 = c2["frame"]
    ok = (qx >= x0) & (qx < x1) & (qy >= y0) & (qy < y1) & (px >= c16["frame"][0]) & (py < c16["frame"][3])
    src = np.full_like(i16, 255.0)
    src[ok] = i2[qy[ok], qx[ok]]
    blue = lambda im: (im[..., 2] > 150) & (im[..., 0] < 120) & (im[..., 1] < 120)
    red = lambda im: (im[..., 0] > 120) & (im[..., 1] < 90) & (im[..., 2] < 90)
    out = 0.35 * i16 + 0.65 * 255.0
    out[blue(i16)] = (0, 0, 200)
    out[red(i16)] = (200, 0, 0)
    out[blue(src)] = (0, 190, 255)
    out[red(src)] = (255, 150, 0)
    fig = plt.figure(figsize=(W / 90, H / 90))
    plt.imshow(out.astype(np.uint8))
    plt.axis("off")
    plt.title("%s (NG blue, RG red) with %s resampled onto its axes (NG cyan, RG orange)" % (f16, f2), fontsize=7)
    fig.savefig(os.path.join(args.out_dir, "EGU2019_slide02_on_slide16_axes.png"), dpi=90, bbox_inches="tight")
    print("check plot: %s" % os.path.join(args.out_dir, "EGU2019_slide02_on_slide16_axes.png"))
