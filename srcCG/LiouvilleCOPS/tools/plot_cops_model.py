#!/usr/bin/env python3
"""Interactive comparison of the COPS reference data with LiouvilleCOPS model runs.

Opens a matplotlib window with four panels: NG and RG pressure together (NG solid lines/dots, RG dashed
lines/squares), NG pressure, RG pressure and the NG/RG ratio. The reference
data are drawn as dots; points not selected by the run chosen with ``--ref``
(useNG/useRG=0: data gaps, excluded intervals, gauge inside the nucleus) are
drawn hollow. When displaying a Python fit, choose that fit with ``--ref`` so
the markers show the fit's actual selection rather than the base input's
selection. Each model run is a line; the check boxes on the right switch
the runs (and the reference data) on and off. The 'All' and 'None' buttons switch all runs.

Schema-v2 Python fit tables are accepted only when their bytes and associated
fit_report.json match fit_manifest.json and their report/profile/selection/
basis identities agree. Their labels show profile, enhancement degree, lambda,
and a short report identity. Unmanifested historical tables are explicitly
labelled legacy.


DATA SOURCES
  Reference data   columns refNG, refRG of COPS_model.dat [Pa]: the COPS NG and RG pressures of the srcCG table
                   srcCG/RosinaMeasurements.h (1150 points, 2016-09-30T00:00:04 .. 10:39:26, given in mbar and
                   converted by LiouvilleCOPS with ReferenceDataUnit = mbar). The table is identical to the ROSINA COPS
                   L2 BG housekeeping pressures of the PSA/PDS3 data set RO-C-ROSINA-2-EXT3-V2.0 (N2-equivalent
                   pressures, 20 C; time = end of each housekeeping interval). 0 in the table = no valid value
                   (e.g. NG 10:24:42-10:33:14). useNG/useRG = 0 marks the points not selected in the run chosen by
                   --ref (no data, its fit selection or ExcludeNG/ExcludeRG/FitInterval, gauge inside the shape).
  Model            columns pNG, pRG [Pa] (and the species columns) of each LiouvilleCOPS run, or of a fit written by
                   tools/fit_cops.py (species columns -1).
  Geometry         UTC and altitude columns: SPICE geometry of the run (altitude = distance to the nearest point of
                   the run's shape model: SHAP5, M001 or M004).

usage:
  plot_cops_model.py [run directories or COPS_model.dat files ...] [--x time|altitude] [--ref RUN] [--save FILE] [--theme light|dark]

Without arguments, output*/COPS_model.dat under <AMPS>/cops_output is loaded;
the legacy application directory is searched only if no new runs exist.
Relative --save paths are written under cops_output (or COPS_OUTPUT_ROOT).
The reference data, useNG/useRG marker masks, and reference geometry are taken from the most recently modified run
that has the useNG/useRG columns, or from --ref RUN. With --x altitude the reference points are placed at the
altitudes of that run (shape model dependent).
"""
import argparse
import glob
import json
import os
import sys
from datetime import datetime
from pathlib import Path

from cops_paths import APP_ROOT, output_root, prepare_output, product_input
from cops_io import (FIT_SCHEMA, RUN_SCHEMA, ProductError, load_manifest,
                     manifest_summary, sha256_file)

import numpy as np
import matplotlib

if "--save" in sys.argv:
    matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.widgets import CheckButtons, Button

#reference data-viz palette: 8 categorical slots (fixed order); more runs reuse the slots with other line styles
PALETTE = {
    "light": ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"],
    "dark":  ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"],
}
CHROME = {
    "light": {"surface": "#fcfcfb", "ink": "#0b0b0b", "ink2": "#52514e", "muted": "#898781", "grid": "#e1e0d9", "axis": "#c3c2b7"},
    "dark":  {"surface": "#1a1a19", "ink": "#f0efec", "ink2": "#c3c2b7", "muted": "#898781", "grid": "#2c2c2a", "axis": "#383835"},
}
LINESTYLES = ["-", "--", ":", "-."]


class Run:
    def __init__(self, fname):
        self.fname = fname
        self.name = os.path.basename(os.path.dirname(os.path.abspath(fname)))
        self.header = []
        names = None
        rows = []

        with open(fname) as f:
            for line in f:
                if line.startswith("#iPoint"):
                    names = line[1:].split()
                elif line.startswith("#"):
                    self.header.append(line[1:].strip())
                elif line.strip():
                    rows.append(line.split())

        if names is None:
            raise ValueError("%s: no column header line (#iPoint ...)" % fname)
        if not rows:
            raise ValueError("%s: no data" % fname)

        # New products declare a table schema.  Unknown schemas are rejected
        # rather than plotted under legacy assumptions.  Unmanifested legacy
        # tables remain readable and are visibly labelled below.
        schemas = [item.split("=", 1)[1] for line in self.header for item in line.split()
                   if item.startswith("TableSchema=")]
        supported = {"liouville_cops_forward_v1", "liouville_cops_fit_prediction_v2"}
        if schemas and schemas[0].rstrip(";") not in supported:
            raise ValueError("unsupported TableSchema=%s" % schemas[0].rstrip(";"))
        try:
            self.kind, self.states = manifest_summary(fname)
        except ProductError as exc:
            raise ValueError("invalid product manifest: %s" % exc) from None

        def header_value(key):
            values = [item.split("=", 1)[1].rstrip(";")
                      for line in self.header for item in line.split()
                      if item.startswith(key + "=")]
            if len(set(values)) > 1:
                raise ValueError("conflicting %s values in table header" % key)
            return values[0] if values else None

        self.table_schema = schemas[0].rstrip(";") if schemas else None
        self.profile = header_value("FitProfile")
        self.selection_identity = header_value("SelectionSHA256")
        self.report_identity = header_value("ReportIdentity")
        self.basis_id = header_value("ActivityBasisID")
        self.enhancement_degree = header_value("EnhancementDegree")
        self.fit_lambda = header_value("Lambda")
        self.regularization = header_value("Regularization")
        self.fit_manifest = None
        self.fit_record = None
        self.run_manifest = None
        self.activity_identity = None
        self.provenance_label = "legacy-unmanifested"
        if self.table_schema == "liouville_cops_fit_prediction_v2":
            self._load_fit_identity()
        elif self.kind != "legacy-unmanifested":
            self._load_forward_identity()

        self.col = {n: i for i, n in enumerate(names)}
        self.t = np.array([datetime.strptime(r[self.col["UTC"]], "%Y-%m-%dT%H:%M:%S") for r in rows])

        def num(name):
            if name not in self.col:
                return None
            return np.array([float(r[self.col[name]]) for r in rows])

        self.ipoint = num("iPoint").astype(int)
        self.alt = num("altitude")
        self.inside = num("inside")
        self.pNG, self.pRG = num("pNG"), num("pRG")
        self.refNG, self.refRG = num("refNG"), num("refRG")
        self.useNG, self.useRG = num("useNG"), num("useRG")

        #the gauge inside the nucleus has no model value
        bad = self.inside > 0
        self.pNG = np.where(bad | (self.pNG <= 0), np.nan, self.pNG)
        self.pRG = np.where(bad | (self.pRG <= 0), np.nan, self.pRG)

        self.mtime = os.path.getmtime(fname)
        self.has_flags = self.useNG is not None

    def _load_fit_identity(self):
        """Require a matching manifest/report for every schema-v2 fit table."""
        table = Path(self.fname).resolve()
        manifest_path = table.parent / "fit_manifest.json"
        if not manifest_path.is_file() and table.parent.name.startswith("lambda="):
            manifest_path = table.parent.parent / "fit_manifest.json"
        if not manifest_path.is_file():
            raise ValueError("fit prediction has no fit_manifest.json")
        manifest = load_manifest(manifest_path, FIT_SCHEMA)
        outputs = manifest.get("outputs", {})
        if table.parent == manifest_path.parent:
            declared = outputs
            result_directory = "."
        else:
            result_directory = table.parent.name
            declared = outputs.get("lambda_directories", {}).get(result_directory, {})
        identity = declared.get("COPS_model.dat")
        if identity is None:
            raise ValueError("fit manifest does not declare this COPS_model.dat")
        if table.stat().st_size != identity.get("size_bytes") or \
                sha256_file(table) != identity.get("sha256"):
            raise ValueError("COPS_model.dat does not match its fit manifest identity")

        records = [record for record in manifest.get("optimizer_results", [])
                   if record.get("result_directory", ".") == result_directory]
        if len(records) != 1:
            raise ValueError("fit manifest has no unique optimizer result for %s" %
                             result_directory)
        record = records[0]
        manifest_report = record.get("report_identity")
        if not self.report_identity or self.report_identity != manifest_report:
            raise ValueError("table and manifest report identities do not agree")

        report_path = table.parent / "fit_report.json"
        report_identity_record = declared.get("fit_report.json")
        if report_identity_record is None or not report_path.is_file():
            raise ValueError("fit report is absent from the declared outputs")
        if (report_path.stat().st_size != report_identity_record.get("size_bytes") or
                sha256_file(report_path) != report_identity_record.get("sha256")):
            raise ValueError("fit_report.json does not match its manifest identity")
        try:
            report = json.loads(report_path.read_text())
        except (OSError, ValueError) as exc:
            raise ValueError("cannot read fit_report.json: %s" % exc) from None
        if report.get("report_identity") != self.report_identity:
            raise ValueError("table and fit_report.json identities do not agree")

        enhancement = manifest.get("activity_enhancement", {})
        expected_basis = (enhancement.get("basis_id") if enhancement.get("enabled")
                          else "legacy_srcCG_25")
        if not self.basis_id or self.basis_id != expected_basis:
            raise ValueError("table and manifest activity basis IDs do not agree")
        expected_profile = manifest.get("selection", {}).get("profile") or \
            "legacy_cli_selection"
        if self.profile != expected_profile:
            raise ValueError("table and manifest selection profiles do not agree")
        if self.selection_identity != manifest.get("selection", {}).get("sha256"):
            raise ValueError("table and manifest selection identities do not agree")

        self.fit_manifest = manifest
        self.fit_record = record
        self.provenance_label = "fit manifest/report verified"

    def _load_forward_identity(self):
        """Verify and expose activity identity for a manifested C++ table."""
        table = Path(self.fname).resolve()
        manifest_path = table.parent / "run_manifest.json"
        if not manifest_path.is_file():
            raise ValueError("manifested forward table has no run_manifest.json")
        manifest = load_manifest(manifest_path, RUN_SCHEMA)
        identity = manifest.get("outputs", {}).get("model_table")
        if identity is None:
            raise ValueError("run manifest does not declare its COPS_model.dat")
        if table.stat().st_size != identity.get("size_bytes") or \
                sha256_file(table) != identity.get("sha256"):
            raise ValueError("COPS_model.dat does not match its run manifest identity")
        enhancement = manifest.get("activity_enhancement", {})
        self.run_manifest = manifest
        if enhancement.get("enabled"):
            self.basis_id = enhancement.get("basis_id")
            self.enhancement_degree = str(enhancement.get("maximum_degree"))
            self.fit_lambda = str(enhancement.get("lambda"))
            coefficient = enhancement.get("coefficient_file", {})
            self.activity_identity = coefficient.get("sha256")
            self.provenance_label = "forward manifest verified; enhancement enabled"
        else:
            self.basis_id = "legacy_srcCG_25"
            self.enhancement_degree = "4"
            self.activity_identity = enhancement.get("legacy_activity_identity")
            self.provenance_label = "forward manifest verified; enhancement disabled"

    def description(self):
        keys = ("ShapeModel=", "FlightTimeRotation=", "SurfaceBC=", "TemperatureModel=", "NGAcceptance=",
                "RGPressureModel=", "TemperatureRGMode=", "ReferenceDataUnit=", "NGHeadOnWeight=")
        out = []
        for h in self.header:
            for item in h.split():
                for k in keys:
                    if item.startswith(k):
                        v = item[len(k):]
                        if k == "ShapeModel=":
                            v = os.path.basename(v)
                        out.append(k + v)
        process = self.states.get("process", "unknown")
        optimizer = self.states.get("optimizer")
        status = "%s process=%s" % (self.kind, process)
        if optimizer and optimizer != "not_applicable":
            status += " optimizer=" + optimizer
        if self.fit_manifest is not None:
            out.extend(["profile=" + self.profile,
                        "basis=" + self.basis_id,
                        "L=" + self.enhancement_degree,
                        "lambda=" + self.fit_lambda,
                        "regularization=" + self.regularization,
                        "report=" + self.report_identity[:12]])
        elif self.run_manifest is not None:
            activity_short = str(self.activity_identity).split(":")[-1][:12]
            out.extend(["basis=" + str(self.basis_id),
                        "L=" + str(self.enhancement_degree),
                        "activity=" + activity_short])
            if self.fit_lambda is not None:
                out.append("lambda=" + self.fit_lambda)
        elif self.kind == "legacy-unmanifested":
            out.append("LEGACY/UNMANIFESTED")
        return status + " " + " ".join(out)

    def display_name(self):
        """Short checkbox label that prevents silent fit-profile mixing."""
        if self.fit_manifest is None:
            if self.run_manifest is not None:
                activity = str(self.activity_identity).split(":")[-1][:8]
                return "%s [forward,L%s,activity=%s]" % (
                    self.name, self.enhancement_degree, activity)
            suffix = " [legacy]" if self.kind == "legacy-unmanifested" else ""
            return self.name + suffix
        return "%s [%s,L%s,lambda=%s]" % (
            self.name, self.profile, self.enhancement_degree, self.fit_lambda)


def find_runs(args):
    if args:
        files = []
        for a in args:
            a = product_input(a)
            files += [os.path.join(a, "COPS_model.dat")] if os.path.isdir(a) else [a]
    else:
        files = sorted(glob.glob(os.path.join(str(output_root()), "output*", "COPS_model.dat")))
        if not files:
            files = sorted(glob.glob(os.path.join(str(APP_ROOT), "output*", "COPS_model.dat")))
            if files:
                print("No new runs found; loading legacy application output directories")

    runs = []
    for f in files:
        try:
            runs.append(Run(f))
        except (OSError, ValueError, KeyError) as e:
            print("skipped %s: %s" % (f, e))
    return runs


def main(default_x="time"):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("runs", nargs="*", help="run directories or COPS_model.dat files (default: all output*/)")
    p.add_argument("--x", choices=["time", "altitude"], default=default_x, help="horizontal axis (default: %(default)s)")
    p.add_argument("--save", default=None, help="save under cops_output (or absolute path outside LiouvilleCOPS) instead of opening a window")
    p.add_argument("--theme", choices=["light", "dark"], default="light")
    p.add_argument("--ref", default=None,
                   help="run that provides reference data, useNG/useRG marker masks, and altitudes "
                        "(default: newest run with useNG/useRG)")
    args = p.parse_args()

    runs = find_runs(args.runs)
    if not runs:
        sys.exit("No model runs found")

    # Directory basenames are not unique in arbitrary comparisons (lambda scan
    # children are a common example).  Use the full scientific label as the UI
    # key and disambiguate exact duplicates deterministically.
    used_labels = {}
    for run in runs:
        base_label = run.display_name()
        used_labels[base_label] = used_labels.get(base_label, 0) + 1
        run.plot_label = (base_label if used_labels[base_label] == 1 else
                          "%s #%d" % (base_label, used_labels[base_label]))

    flagged = [r for r in runs if r.has_flags]
    if args.ref:
        sel = [r for r in runs if r.name == args.ref or r.name == "output_" + args.ref]
        if not sel:
            sys.exit("--ref %s: no such run among the loaded runs" % args.ref)
        ref = sel[0]
    else:
        ref = max(flagged, key=lambda r: r.mtime) if flagged else max(runs, key=lambda r: r.mtime)

    print("Reference data and selection masks from: %s" % ref.fname)
    print("Runs:")
    for i, r in enumerate(runs):
        print("  %2d %-24s %s" % (i, r.name, r.description()))

    c = CHROME[args.theme]
    pal = PALETTE[args.theme]

    plt.rcParams["font.family"] = "sans-serif"
    fig = plt.figure(figsize=(14, 12))
    fig.patch.set_facecolor(c["surface"])

    #panels on the left (both pressures together, NG, RG, NG/RG), controls on the right
    axBoth = fig.add_axes([0.07, 0.755, 0.66, 0.205])
    axNG = fig.add_axes([0.07, 0.525, 0.66, 0.205], sharex=axBoth)
    axRG = fig.add_axes([0.07, 0.295, 0.66, 0.205], sharex=axBoth)
    axR = fig.add_axes([0.07, 0.065, 0.66, 0.205], sharex=axBoth)

    def X(r):
        return r.t if args.x == "time" else r.alt

    artists = {}

    #reference data: used (filled) and not used (hollow) points
    def ref_points(ax, y, use, marker="o"):
        x = X(ref)
        ok = y > 0
        used = ok & (use > 0) if use is not None else ok
        unused = ok & ~used
        ms = 2.6 if marker == "o" else 2.4
        a1 = ax.plot(x[used], y[used], marker, ls="none", ms=ms, color=c["ink"], label="reference", zorder=4)[0]
        a2 = ax.plot(x[unused], y[unused], marker, ls="none", ms=ms+0.4, mfc="none", mec=c["muted"], mew=0.7, zorder=4)[0]
        return [a1, a2]

    artists["reference"] = ref_points(axNG, ref.refNG, ref.useNG) + ref_points(axRG, ref.refRG, ref.useRG)

    #combined panel: NG = solid line / dots, RG = dashed line / squares
    artists["reference"] += ref_points(axBoth, ref.refNG, ref.useNG, "o") + ref_points(axBoth, ref.refRG, ref.useRG, "s")

    rr = np.where((ref.refNG > 0) & (ref.refRG > 0), ref.refNG / np.where(ref.refRG > 0, ref.refRG, 1.0), np.nan)
    useR = ((ref.useNG > 0) & (ref.useRG > 0)) if ref.has_flags else np.isfinite(rr)
    artists["reference"] += [axR.plot(X(ref)[useR], rr[useR], "o", ms=2.6, color=c["ink"], zorder=4)[0],
                             axR.plot(X(ref)[~useR & np.isfinite(rr)], rr[~useR & np.isfinite(rr)], "o", ms=3.0,
                                      mfc="none", mec=c["muted"], mew=0.7, zorder=4)[0]]

    colors = {}
    for i, r in enumerate(runs):
        col = pal[i % len(pal)]
        ls = LINESTYLES[(i // len(pal)) % len(LINESTYLES)]
        colors[r.plot_label] = (col, ls)
        kw = dict(color=col, ls=ls, lw=1.3, zorder=3)
        artists[r.plot_label] = [axNG.plot(X(r), r.pNG, **kw)[0],
                                 axRG.plot(X(r), r.pRG, **kw)[0],
                                 axR.plot(X(r), r.pNG / r.pRG, **kw)[0],
                                 axBoth.plot(X(r), r.pNG, color=col, ls="-", lw=1.4, zorder=3)[0],
                                 axBoth.plot(X(r), r.pRG, color=col, ls="--", lw=1.2, zorder=3)[0]]

    from matplotlib.lines import Line2D
    keyNG = Line2D([], [], color=c["ink2"], ls="-", marker="o", ms=3.5, lw=1.4, label="NG: model solid, data dots")
    keyRG = Line2D([], [], color=c["ink2"], ls="--", marker="s", ms=3.2, lw=1.2, label="RG: model dashed, data squares")
    leg = axBoth.legend(handles=[keyNG, keyRG], loc="best", fontsize=8, frameon=True, facecolor=c["surface"], edgecolor=c["grid"], framealpha=0.9)
    for t in leg.get_texts():
        t.set_color(c["ink2"])

    for ax, ylab in ((axBoth, "NG and RG pressure [Pa]"), (axNG, "NG pressure [Pa]"), (axRG, "RG pressure [Pa]"), (axR, "NG / RG")):
        ax.set_yscale("log")
        ax.set_ylabel(ylab, color=c["ink2"])
        ax.set_facecolor(c["surface"])
        ax.grid(True, color=c["grid"], lw=0.6)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            ax.spines[s].set_color(c["axis"])
        ax.tick_params(colors=c["muted"], labelcolor=c["ink2"], labelsize=8)

    axR.axhline(1.0, color=c["axis"], lw=0.8)

    if args.x == "time":
        axR.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
        axR.set_xlabel("UTC, %s" % ref.t[0].date(), color=c["ink2"])
    else:
        axR.set_xscale("log")
        axR.set_xlabel("altitude [m]", color=c["ink2"])
    for ax in (axBoth, axNG, axRG):
        plt.setp(ax.get_xticklabels(), visible=False)

    ref_context = ("profile=%s, report=%s" %
                   (ref.profile, ref.report_identity[:12])
                   if ref.fit_manifest is not None else ref.provenance_label)
    axBoth.set_title("COPS reference data and LiouvilleCOPS runs "
                     "(hollow: not selected by %s; %s)" % (ref.name, ref_context),
                     loc="left", fontsize=11, color=c["ink"])

    #controls: one check box per run, colored like its line
    labels = ["reference"] + [r.plot_label for r in runs]
    h = min(0.86, 0.035 * len(labels))
    axCheck = fig.add_axes([0.75, 0.96 - h, 0.24, h])
    axCheck.set_facecolor(c["surface"])
    for s in axCheck.spines.values():
        s.set_visible(False)
    check = CheckButtons(axCheck, labels, [True] * len(labels))

    #identity: a line sample (color and line style of the run) between the box and the label
    for i, lab in enumerate(labels):
        txt = check.labels[i]
        txt.set_color(c["ink2"])
        txt.set_fontsize(8)
        y = txt.get_position()[1]
        txt.set_x(0.30)
        check.rectangles[i].set_edgecolor(c["axis"])
        check.rectangles[i].set_facecolor(c["surface"])
        for ln in check.lines[i]:
            ln.set_color(c["ink"])
        if lab == "reference":
            axCheck.plot([0.20], [y], "o", ms=4, color=c["ink"], transform=axCheck.transAxes)
        else:
            col, ls = colors[lab]
            axCheck.plot([0.13, 0.27], [y, y], color=col, ls=ls, lw=1.6, transform=axCheck.transAxes)
    axCheck.set_xlim(0, 1), axCheck.set_ylim(0, 1)

    def toggle(label):
        for a in artists[label]:
            a.set_visible(not a.get_visible())
        fig.canvas.draw_idle()

    check.on_clicked(toggle)

    def set_all(state):
        for i, lab in enumerate(labels):
            if lab == "reference":
                continue
            if check.get_status()[i] != state:
                check.set_active(i)   #calls toggle()

    axAll = fig.add_axes([0.75, 0.90 - h - 0.02, 0.07, 0.035])
    axNone = fig.add_axes([0.83, 0.90 - h - 0.02, 0.07, 0.035])
    bAll, bNone = Button(axAll, "All"), Button(axNone, "None")
    bAll.on_clicked(lambda e: set_all(True))
    bNone.on_clicked(lambda e: set_all(False))

    if args.save:
        saved = prepare_output(args.save)
        identities = ",".join(
            r.report_identity or
            ("forward:" + str(r.activity_identity)
             if r.run_manifest is not None else "legacy-unmanifested")
            for r in runs)
        fig.savefig(saved, dpi=130, facecolor=c["surface"], metadata={
            "Title": "LiouvilleCOPS pressure comparison",
            "Description": "report identities=" + identities,
        })
        print("Figure written: %s" % saved)
    else:
        plt.show()


if __name__ == "__main__":
    main()
