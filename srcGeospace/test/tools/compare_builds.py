#!/usr/bin/env python3
"""Cross-build comparison report with plots (e.g. GPU vs CPU, float vs double).

Usage:
    compare_builds.py A_ROOT B_ROOT -o OUT_DIR [--label-a CPU-double] [--label-b GPU-double]
                      [--max-plots 80]

A_ROOT and B_ROOT are run roots written by `run_srcearth_tests.sh --run-root ... --amps ...`
for the same test list with two executables; A is the reference.  The report goes to
OUT_DIR (outside the source tree):

  REPORT.md       test status and run time of both builds, output-file comparison per
                  test, per-variable difference statistics, links to the plots
  files.csv       compare_outputs.py classification of every output file
  variables.csv   variable_stats.py statistics of every differing Tecplot variable
  plots/*.png     for each Tecplot file with differing values:
                  * shell grids (lon_deg/lat_deg columns, one zone per altitude): maps of
                    the relative difference (B-A)/max(|A|,|B|) of the main variable
                    (Rc_GV if present) per zone, and B against A;
                  * other files: B against A for up to four differing variables.

The individual test runners still make their own comparison plots (model vs reference
data or analytic solution) inside each run root; this report adds the comparison of
the two builds with each other.  Neither root is modified.
"""

import argparse
import csv
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import compare_outputs   # noqa: E402
import variable_stats    # noqa: E402

VARS = re.compile(r'VARIABLES\s*=\s*(.*)', re.IGNORECASE)
ZONE = re.compile(r'^ZONE\b(?:.*?T\s*=\s*"([^"]*)")?', re.IGNORECASE)
COORDS = ("lon_deg", "lat_deg", "alt_km", "x", "y", "z", "E_MeV", "energy_MeV")


def read_zones(path):
    """(variable names, [(zone title, rows)]) of a Tecplot ASCII POINT file."""
    names, zones = None, []
    with open(path, errors="replace") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            m = VARS.match(s)
            if m:
                names = [v.strip().strip('"') for v in re.split(r'[,\s]+', m.group(1)) if v.strip()]
                continue
            z = ZONE.match(s)
            if z:
                zones.append((z.group(1) or "zone %d" % (len(zones) + 1), []))
                continue
            if s[0].isalpha():
                continue
            try:
                vals = [float(t) for t in s.replace(",", " ").split()]
            except ValueError:
                continue
            if not zones:
                zones.append(("", []))
            zones[-1][1].append(vals)
    return names, zones


def load_report(root):
    path = os.path.join(root, "test_runner_report.json")
    if not os.path.exists(path):
        return {}
    with open(path) as f:
        return {r["command"]: r for r in json.load(f)}


def short(cmd, width=70):
    cmd = cmd.replace("python3 ", "").replace("./", "")
    cmd = re.sub(r"src(?:Earth|Geospace)/test/", "", cmd)
    return cmd if len(cmd) <= width else cmd[: width - 3] + "..."


def plot_file(plt, a_path, b_path, test, name, la, lb, out_png):
    na, za = read_zones(a_path)
    nb, zb = read_zones(b_path)
    if not na or na != nb or len(za) != len(zb):
        return False
    if any(len(x[1]) != len(y[1]) for x, y in zip(za, zb)):
        return False
    differing = []
    for j, var in enumerate(na):
        if var in COORDS:
            continue
        if any(ra[j] != rb[j] for (_, A), (_, B) in zip(za, zb) for ra, rb in zip(A, B)
               if len(ra) == len(na) and len(rb) == len(na)):
            differing.append(j)
    if not differing:
        return False
    is_map = "lon_deg" in na and "lat_deg" in na
    if is_map:
        main = na.index("Rc_GV") if "Rc_GV" in na and na.index("Rc_GV") in differing else differing[0]
        ilon, ilat = na.index("lon_deg"), na.index("lat_deg")
        nz = len(za)
        fig, axes = plt.subplots(1, nz + 1, figsize=(4.2 * (nz + 1), 3.6), squeeze=False)
        allx, ally = [], []
        for k, ((title, A), (_, B)) in enumerate(zip(za, zb)):
            lon = [r[ilon] for r in A]
            lat = [r[ilat] for r in A]
            rel = []
            for ra, rb in zip(A, B):
                s = max(abs(ra[main]), abs(rb[main]))
                rel.append((rb[main] - ra[main]) / s if s > 0 else 0.0)
                allx.append(ra[main])
                ally.append(rb[main])
            lim = max(1e-12, max(abs(v) for v in rel))
            sc = axes[0][k].scatter(lon, lat, c=rel, cmap="coolwarm", vmin=-lim, vmax=lim, s=28)
            axes[0][k].set_title(f"{title}\nmax |rel diff| {lim:.2e}", fontsize=8)
            axes[0][k].set_xlabel("lon (deg)")
            axes[0][k].set_ylabel("lat (deg)")
            fig.colorbar(sc, ax=axes[0][k], fraction=0.046)
        ax = axes[0][nz]
        ax.scatter(allx, ally, s=8)
        lo, hi = min(allx + ally), max(allx + ally)
        ax.plot([lo, hi], [lo, hi], "k--", lw=0.8)
        ax.set_xlabel(f"{na[main]} {la}")
        ax.set_ylabel(f"{na[main]} {lb}")
        ax.set_title("B vs A", fontsize=8)
        fig.suptitle(f"{test}/{name}: {na[main]}, relative difference ({lb} - {la})", fontsize=9)
    else:
        sel = differing[:4]
        fig, axes = plt.subplots(1, len(sel), figsize=(4.0 * len(sel), 3.6), squeeze=False)
        for k, j in enumerate(sel):
            xs = [ra[j] for (_, A), (_, B) in zip(za, zb) for ra, rb in zip(A, B)]
            ys = [rb[j] for (_, A), (_, B) in zip(za, zb) for ra, rb in zip(A, B)]
            ax = axes[0][k]
            ax.scatter(xs, ys, s=8)
            lo, hi = min(xs + ys), max(xs + ys)
            ax.plot([lo, hi], [lo, hi], "k--", lw=0.8)
            pos = [v for v in xs + ys if v > 0]
            if len(pos) == len(xs + ys) and pos and max(pos) / min(pos) > 100:
                ax.set_xscale("log")
                ax.set_yscale("log")
            ax.set_xlabel(f"{na[j]} {la}", fontsize=8)
            ax.set_ylabel(f"{na[j]} {lb}", fontsize=8)
        fig.suptitle(f"{test}/{name}: {lb} vs {la}", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a_root")
    ap.add_argument("b_root")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--label-a", default="A")
    ap.add_argument("--label-b", default="B")
    ap.add_argument("--max-plots", type=int, default=80)
    args = ap.parse_args()
    la, lb = args.label_a, args.label_b
    os.makedirs(os.path.join(args.out, "plots"), exist_ok=True)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    md = [f"# Build comparison: {lb} vs {la}", "",
          f"* A ({la}): `{os.path.abspath(args.a_root)}`",
          f"* B ({lb}): `{os.path.abspath(args.b_root)}`", ""]

    # 1. Test status and run time.
    ra, rb = load_report(args.a_root), load_report(args.b_root)
    if ra or rb:
        md += ["## 1. Test status and run time", "",
               f"| test | {la} | {lb} | {la} s | {lb} s | {lb}/{la} |", "|---|---|---|---|---|---|"]
        ta = tb = 0.0
        for cmd in list(ra) + [c for c in rb if c not in ra]:
            a, b = ra.get(cmd), rb.get(cmd)
            sa = f"{a['actual']}{'' if a['matched_reference'] else '!'}" if a else "-"
            sb = f"{b['actual']}{'' if b['matched_reference'] else '!'}" if b else "-"
            ea = a["elapsed_s"] if a else None
            eb = b["elapsed_s"] if b else None
            ta += ea or 0.0
            tb += eb or 0.0
            ratio = f"{eb / ea:.2f}" if ea and eb else "-"
            fa = f"{ea:.1f}" if ea is not None else "-"
            fb = f"{eb:.1f}" if eb is not None else "-"
            md.append(f"| `{short(cmd)}` | {sa} | {sb} | {fa} | {fb} | {ratio} |")
        md += [f"| **total** | | | {ta:.0f} | {tb:.0f} | {tb / ta:.2f} |" if ta else "", "",
               "`!` marks a status that differs from the expectation in the list.", ""]

    # 2. Output files.
    rows = compare_outputs.compare_roots(args.a_root, args.b_root)
    with open(os.path.join(args.out, "files.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["file", "status", "max_abs", "max_rel", "note"])
        w.writerows(rows)
    per_test = {}
    for r in rows:
        per_test.setdefault(r[0].split(os.sep)[0], {}).setdefault(r[1], 0)
        per_test[r[0].split(os.sep)[0]][r[1]] += 1
    statuses = ["IDENTICAL", "IDENTICAL_AFTER_PATHS", "NONDETERMINISTIC", "NUMERIC", "DIFFERENT", "MISSING"]
    md += ["## 2. Output files (`files.csv`)", "",
           "| test output | " + " | ".join(statuses) + " |", "|---" * (len(statuses) + 1) + "|"]
    for t in sorted(per_test):
        md.append(f"| {t} | " + " | ".join(str(per_test[t].get(s, "")) for s in statuses) + " |")
    md.append("")

    # 3. Variables.
    vrows = variable_stats.collect(args.a_root, args.b_root)
    with open(os.path.join(args.out, "variables.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["test", "file", "variable", "values", "identical_fraction",
                    "median_rel_of_differing", "p90_rel", "p99_rel", "max_rel"])
        w.writerows(vrows)
    md += ["## 3. Differing Tecplot variables (`variables.csv`)", "",
           "Relative difference |a-b|/max(|a|,|b|) over the values that differ; values are",
           "compared as printed.", "",
           "| test | file | variable | values | identical | median | p90 | max |", "|---|---|---|---|---|---|---|---|"]
    for t, fn, var, n, ident, med, p90, p99, mx in vrows:
        md.append(f"| {t} | {fn} | {var} | {n} | {ident:.3f} | {med:.1e} | {p90:.1e} | {mx:.1e} |")
    md.append("")

    # 4. Plots.
    md += ["## 4. Plots", ""]
    seen, nplots = set(), 0
    for t, fn, *_ in vrows:
        if (t, fn) in seen or nplots >= args.max_plots:
            continue
        seen.add((t, fn))
        for dirpath, _, files in os.walk(os.path.join(args.a_root, "test_output", t)):
            if fn not in files:
                continue
            a_path = os.path.join(dirpath, fn)
            relp = os.path.relpath(a_path, os.path.join(args.a_root, "test_output"))
            b_path = os.path.join(args.b_root, "test_output", relp)
            if not os.path.exists(b_path) or nplots >= args.max_plots:
                continue
            png = "plots/" + relp.replace(os.sep, "__").replace(".dat", ".png")
            if plot_file(plt, a_path, b_path, t, relp, la, lb, os.path.join(args.out, png)):
                md += [f"### {relp}", "", f"![{relp}]({png})", ""]
                nplots += 1
    if nplots >= args.max_plots:
        md.append(f"(plot limit {args.max_plots} reached; use --max-plots)")
    with open(os.path.join(args.out, "REPORT.md"), "w") as f:
        f.write("\n".join(md) + "\n")
    print(f"report: {os.path.join(args.out, 'REPORT.md')} ({nplots} plots)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
