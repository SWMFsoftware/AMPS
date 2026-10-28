#!/usr/bin/env python3
"""Number of surface elements seen by the COPS nude gauge (NG) only or by the ram gauge (RG) only.

Builds and runs surface_contributions.cpp (this directory) with the unchanged LiouvilleCOPS sources
(../../LiouvilleCOPS/src) for the observation times in a UTC and/or altitude interval, and reports the number
(and area) of the surface elements that are

  * NG only / RG only at least once in the interval,
  * NG only / RG only all the time, strict: at every sample of the interval,
  * NG only / RG only all the time, loose:  at every sample at which either gauge sees the element,

plus the elements seen by NG, RG, both and either at least once.

"Only" (user decision 2026-10-27): the contribution of an element is the gauge pressure response to a unit
uniform H2O flux emitted with the temperature --temperature; NG only means H_NG > 0 and |H_RG| <= eps H_NG,
RG only means |H_RG| > 0 and H_NG <= eps |H_RG| (eps = --epsilon, default 1e-3). The gauge acceptance,
pressure conversion and transport settings are those of the LiouvilleCOPS input template (--template, default
LiouvilleCOPS/fit_base.in, the settings of the fits), with the paths replaced by the command-line values.

Example (12M, the late window by UTC; run it inside a Slurm job, it compiles with -march=native):
  python3 surface_contributions.py --mesh <...>/cg-spg_mspcd-shap8-v0.4cmb-12Mfacets.bdf \\
      --kernel-dir /home/vtenishev/SPICE/Kernels/ROSETTA/esa_2022 --cspice /home/vtenishev/SPICE/cspice \\
      --utc-start 2016-09-30T10:34:26 --utc-end 2016-09-30T10:39:26 --temperature 150 --out <dir>
or by altitude:  --alt-min 7 --alt-max 300   (metres above the closest surface point)

Products in --out: input.in (the generated LiouvilleCOPS input), build/ (binary, compile log), run.log,
samples.dat (per-sample counts), counts.json (interval counts and areas), results.txt (the table printed at the
end), facets.npz with --write-facets (per-element sample counts and class flags), command.json.
"""

import argparse
import json
import os
import re
import shlex
import socket
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

HERE = Path(__file__).resolve().parent
LC = (HERE / ".." / ".." / "LiouvilleCOPS").resolve()
DEFAULT_OBS = (HERE / ".." / ".." / "RosinaMeasurements.h").resolve()

# keys of the template that this tool sets itself (the rest — gauge, transport, SPICE frames — is kept)
REPLACED = {"ShapeModel", "ShapeModelFormat", "KernelDirectory", "ObservationFile", "OutputDirectory", "WriteKernels",
            "FirstPoint", "LastPoint", "PointStep", "Threads", "ExcludeNG", "ExcludeRG", "FitInterval", "Mode",
            "AllowExistingOutputDirectory", "GeometryOnly", "TemperatureRGTable"}

LABELS = [("seen_by_NG_at_least_once", "seen by NG, at least once"),
          ("seen_by_RG_at_least_once", "seen by RG, at least once"),
          ("seen_by_both_at_least_once", "seen by both, at least once"),
          ("seen_by_either_at_least_once", "seen by either, at least once"),
          ("NG_only_at_least_once", "NG only, at least once"),
          ("RG_only_at_least_once", "RG only, at least once"),
          ("NG_only_all_the_time_strict", "NG only, all the time (strict: every sample)"),
          ("RG_only_all_the_time_strict", "RG only, all the time (strict: every sample)"),
          ("NG_only_all_the_time_loose", "NG only, all the time (loose: whenever seen)"),
          ("RG_only_all_the_time_loose", "RG only, all the time (loose: whenever seen)")]


def parse():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog=__doc__[__doc__.index("\n") + 1:])
    g = p.add_argument_group("locations")
    g.add_argument("--mesh", required=True, help="shape model (NASTRAN .bdf CTRIA3 surface, or SPICE DSK with --mesh-format dsk)")
    g.add_argument("--mesh-format", choices=["nastran", "dsk"], help="default: dsk for .bds/.dsk, else nastran")
    g.add_argument("--kernel-dir", required=True, help="SPICE kernel directory (KernelDirectory)")
    g.add_argument("--kernels", nargs="+", help="SPICE kernels, relative to --kernel-dir or absolute (default: the template's list)")
    g.add_argument("--cspice", required=True, help="CSPICE toolkit root (include/, lib/cspice.a)")
    g.add_argument("--observations", default=str(DEFAULT_OBS), help="observation table (default %(default)s)")
    g.add_argument("--rg-temperature-table", help="COPS HK table for TemperatureRGMode=hk_table (default: the template's, "
                   "relative paths under $COPS_OUTPUT_ROOT)")
    g.add_argument("--template", default=str(LC / "fit_base.in"), help="LiouvilleCOPS input with the gauge/transport settings (default %(default)s)")
    g.add_argument("--set", action="append", default=[], metavar="KEY=VALUE", help="override a template key (repeatable)")
    g.add_argument("--liouville-src", default=str(LC / "src"), help="LiouvilleCOPS sources (default %(default)s)")
    s = p.add_argument_group("interval and definition")
    s.add_argument("--utc-start", help="first UTC (e.g. 2016-09-30T10:34:26)")
    s.add_argument("--utc-end", help="last UTC")
    s.add_argument("--alt-min", type=float, help="lowest altitude above the surface [m]")
    s.add_argument("--alt-max", type=float, help="highest altitude above the surface [m]")
    s.add_argument("--stride", type=int, default=1, help="use every n-th observation time in the interval (default 1)")
    s.add_argument("--epsilon", type=float, default=1e-3, help="'only' threshold on the contribution ratio (default %(default)g)")
    s.add_argument("--temperature", type=float, required=True, help="gas emission temperature [K] of the unit H2O flux")
    r = p.add_argument_group("run")
    r.add_argument("--out", required=True, help="output directory (new, or empty)")
    r.add_argument("--threads", type=int, default=os.cpu_count(), help="OpenMP threads (default: all cores, %(default)s)")
    r.add_argument("--cxx", default="g++")
    r.add_argument("--cxxflags", default="-O3 -march=native -std=c++17 -fopenmp")
    r.add_argument("--write-facets", action="store_true", help="also write facets.npz (per-element counts and flags)")
    a = p.parse_args()
    if not any(v is not None for v in (a.utc_start, a.utc_end, a.alt_min, a.alt_max)):
        p.error("give a UTC interval (--utc-start/--utc-end) and/or an altitude interval (--alt-min/--alt-max)")
    if not 0 < a.epsilon < 1:
        p.error("--epsilon must be in (0, 1)")
    return a


def read_template(path):
    """(key, value) pairs of a LiouvilleCOPS input, comments stripped, order kept."""
    pairs = []
    for line in Path(path).read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if "=" in line:
            k, v = (x.strip() for x in line.split("=", 1))
            pairs.append((k, v))
    return pairs


def write_input(a, out):
    pairs = read_template(a.template)
    keys = {k for k, _ in pairs}
    sets = dict(x.split("=", 1) for x in a.set)
    rg_table = a.rg_temperature_table or dict(pairs).get("TemperatureRGTable")
    if rg_table and not os.path.isabs(rg_table):
        root = os.environ.get("COPS_OUTPUT_ROOT")
        if not root:
            sys.exit("error: the template's TemperatureRGTable is relative (%s); set COPS_OUTPUT_ROOT or give "
                     "--rg-temperature-table" % rg_table)
        rg_table = str(Path(root) / rg_table)
    mode = sets.get("TemperatureRGMode", dict(pairs).get("TemperatureRGMode", "constant"))
    if mode == "hk_table" and not (rg_table and Path(rg_table).is_file()):
        sys.exit("error: TemperatureRGMode=hk_table but the table is not found: %s" % rg_table)
    fmt = a.mesh_format or ("dsk" if a.mesh.lower().endswith((".bds", ".dsk")) else "nastran")
    lines = ["# generated by %s from %s" % (Path(__file__).name, a.template)]
    for k, v in pairs:
        if k in REPLACED or k in sets or (k == "Kernel" and a.kernels):
            continue
        lines.append("%s = %s" % (k, v))
    for k in a.kernels or []:
        lines.append("Kernel = %s" % k)
    lines += ["KernelDirectory = %s" % a.kernel_dir, "ShapeModel = %s" % os.path.abspath(a.mesh),
              "ShapeModelFormat = %s" % fmt, "ObservationFile = %s" % os.path.abspath(a.observations),
              "Threads = %d" % a.threads]
    if rg_table:
        lines.append("TemperatureRGTable = %s" % rg_table)
    lines += ["%s = %s" % kv for kv in sets.items()]
    for k in ("Kernel",):
        if k not in keys and not a.kernels:
            sys.exit("error: no SPICE kernels (template has none; give --kernels)")
    path = out / "input.in"
    path.write_text("\n".join(lines) + "\n")
    return path


def compile_tool(a, out):
    build = out / "build"
    build.mkdir()
    src = Path(a.liouville_src)
    sources = sorted(p for p in src.glob("*.cpp") if p.name != "main.cpp") + [HERE / "surface_contributions.cpp"]
    flags = shlex.split(a.cxxflags) + ["-I%s" % src, "-I%s" % Path(a.cspice, "include")]
    log = open(build / "compile.log", "w")

    def cc(p):
        cmd = [a.cxx] + flags + ["-c", str(p), "-o", str(build / (p.stem + ".o"))]
        r = subprocess.run(cmd, capture_output=True, text=True)
        return cmd, r
    print("compiling %d sources on %s ..." % (len(sources), socket.gethostname()), flush=True)
    with ThreadPoolExecutor(max_workers=min(len(sources), os.cpu_count() or 1)) as ex:
        results = list(ex.map(cc, sources))
    for cmd, r in results:
        log.write(" ".join(cmd) + "\n" + r.stdout + r.stderr)
        if r.returncode:
            log.close()
            sys.exit("error: compilation failed (%s); see %s\n%s" % (cmd[-3], build / "compile.log", r.stderr[-2000:]))
    exe = build / "surface_contributions"
    cmd = [a.cxx] + [str(build / (p.stem + ".o")) for p in sources] + ["-o", str(exe), str(Path(a.cspice, "lib", "cspice.a")),
                                                                      "-lm", "-fopenmp"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    log.write(" ".join(cmd) + "\n" + r.stdout + r.stderr)
    log.write(subprocess.run([a.cxx, "--version"], capture_output=True, text=True).stdout)
    log.close()
    if r.returncode:
        sys.exit("error: link failed; see %s\n%s" % (build / "compile.log", r.stderr[-2000:]))
    return exe


def main():
    a = parse()
    out = Path(a.out).resolve()
    if out.exists() and any(out.iterdir()):
        sys.exit("error: %s exists and is not empty (products are never overwritten)" % out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "command.json").write_text(json.dumps({"argv": sys.argv, "cwd": os.getcwd(), "host": socket.gethostname(),
                                                  "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                                                  "args": vars(a)}, indent=1))
    inp = write_input(a, out)
    exe = compile_tool(a, out)
    opt = lambda v: "-" if v is None else str(v)
    cmd = [str(exe), str(inp), str(out), repr(a.epsilon), repr(a.temperature), str(a.stride), opt(a.utc_start),
           opt(a.utc_end), opt(a.alt_min), opt(a.alt_max), "1" if a.write_facets else "0"]
    env = dict(os.environ, OMP_NUM_THREADS=str(a.threads))
    print("running: %s" % " ".join(cmd), flush=True)
    t0 = time.time()
    with open(out / "run.log", "w") as log:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
        for line in proc.stdout:
            log.write(line)
            log.flush()
            if line.startswith(("point ", "samples:", "Done", "Error")):
                print(line, end="", flush=True)
        rc = proc.wait()
    if rc:
        sys.exit("error: surface_contributions exited %d; see %s" % (rc, out / "run.log"))

    res = json.loads((out / "counts.json").read_text())
    if a.write_facets:
        import numpy as np
        raw = (out / "facets.bin").read_bytes()
        nf, ns = np.frombuffer(raw, np.int32, 2)
        u = np.frombuffer(raw, np.uint32, 5 * nf, offset=8).reshape(5, nf)
        flags = np.frombuffer(raw, np.uint8, nf, offset=8 + 20 * nf)
        np.savez_compressed(out / "facets.npz", n_samples=ns, n_seen_ng=u[0], n_seen_rg=u[1], n_seen_any=u[2],
                            n_ng_only=u[3], n_rg_only=u[4], flags=flags,
                            flag_bits=np.array(["NG_only_at_least_once", "RG_only_at_least_once", "NG_only_strict",
                                                "RG_only_strict", "NG_only_loose", "RG_only_loose", "seen_by_NG",
                                                "seen_by_RG"]))
        (out / "facets.bin").unlink()

    interval = []
    if a.utc_start or a.utc_end:
        interval.append("UTC %s .. %s" % (a.utc_start or "start", a.utc_end or "end"))
    if a.alt_min is not None or a.alt_max is not None:
        interval.append("altitude %s .. %s m" % (opt(a.alt_min), opt(a.alt_max)))
    txt = ["Surface elements seen by NG only / RG only",
           "mesh %s (%d elements, %.6e m^2)" % (a.mesh, res["n_facets"], res["total_area_m2"]),
           "interval %s; samples %d (%s .. %s; %d in the interval, stride %d, %d skipped: gauge inside the nucleus)" % (
               "; ".join(interval), res["n_samples"], res["first_sample_utc"], res["last_sample_utc"],
               res["n_samples_in_interval"], res["stride"], res["n_samples_gauge_inside_nucleus"]),
           "only: contribution ratio <= epsilon %g (pressure response to a unit uniform H2O flux at %g K)" % (
               res["epsilon"], res["temperature_K"]),
           "",
           "%-48s %12s %14s %9s" % ("class", "elements", "area [m^2]", "area [%]")]
    for k, label in LABELS:
        txt.append("%-48s %12d %14.6e %9.4f" % (label, res["counts"][k], res["areas_m2"][k],
                                                 100 * res["areas_m2"][k] / res["total_area_m2"]))
    txt.append("")
    txt.append("per-sample counts: %s; run time %.1f s" % (out / "samples.dat", time.time() - t0))
    (out / "results.txt").write_text("\n".join(txt) + "\n")
    print("\n".join(txt))


if __name__ == "__main__":
    main()
