#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/nightfix_base_q002.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/nightfix_base_q002.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: f66841b7c5dd10036a0a821820658ef69340a78de9f5907aaa209b83276739bc
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: tool (2026-11-18, Phase 0 of the H2O night-side fix, route ii). compare-sample: C++ re-runs of sample rows of
#   the lategrid base vs the base (kernel files / fields bitwise, geometry/data columns as text, C++ new-law pNG/pRG vs
#   fit_cops' own new-law model from the existing kernels); assemble: runs/shap8_012m_lategrid_dz1.5/base_q002 from hard
#   links of the base's kernels / elements.bin and the full new-law C++ run's table and manifest.
#   v1 (frozen q002_base_20261118T081402Z) crashed (reloaded only the used rows); v2 (..082517Z) fixed that; v3 (..082717Z,
#   this version): NGBelowDeckFraction is flux-weighted (Liouville.cpp:159/173/347), reported, not gated.
# USED BY: jobs/q002_base.sbatch (136506), jobs/q002_compare.sbatch (136519, v2), jobs/q002_compare_v3.sbatch (136521),
#   jobs/q002_assemble.sbatch.
# RUN: see tests/q002_base_20261118/ (jobs/, README of runs/shap8_012m_lategrid_dz1.5/base_q002)
# 2026-11-18: v4 (frozen q002_base_20261118T083938Z): assemble for a WriteKernels = off run (no elements.bin there: the sample runs' byte-identical elements.bin is the evidence; no run_manifest.json in the base, the C++ manifest kept as cpp_run_manifest.json); job 136522 had failed on the missing elements.bin
# ==== end of header ====
"""Phase 0 of the H2O night-side fix (lead 2026-11-18, user decision "(ii); refit D2 first, then all stage-2 chains"):
a fitting base for the new law (q_min/q_max = 0.02) WITHOUT rebuilding the 172 GB kernel set.

The kernel records (face, Kn, KPhi, KV, cos, shadow) are geometry: they do not depend on the production law. Only the
model columns of COPS_model.dat (written by the C++ at base-build time with the old law) do, and fit_cops.main compares
its model at the base's own parameters with them before every fit (tolerance 1e-6). This tool:

  compare-sample  NON-CIRCULAR check on rows re-run by the REBUILT C++ (native/liouville_cops_cuda, new Source.cpp, the
                  base's own input with only FirstPoint / LastPoint / PointStep / OutputDirectory changed):
                  (a) every re-run kernel file vs the existing base's: sha256 of the whole file and, field by field,
                      face, Kn, KPhi, KV, cos, shadow and the header (bitwise); elements.bin sha256;
                      the geometry / data columns of the re-run table rows (UTC, et, x, y, z, altitude, inside, rh_AU,
                      VdotB_RG, nVisible, refNG, refRG, NGBelowDeckFraction, T_RG, useNG, useRG) = the base's, as text
                  (b) the re-run C++ table's pNG / pRG (new law) vs fit_cops' OWN model (fit_cops.Model.evaluate with
                      the run's coefficients, T0 and kappa = 1, i.e. exactly fit_cops.main's consistency check, FP64
                      records) computed from the EXISTING base's kernels of the same rows; also vs the old table
  assemble        the new base folder: kernels/ = HARD LINKS of the existing base's kernel files, elements.bin = a hard
                  link of the existing one (only if the C++ run's elements.bin is byte-identical), COPS_model.dat,
                  run_manifest.json, input.snapshot.in, boundary.in, build.log = the full C++ run's own files (copied
                  unchanged); refuses unless every geometry / data column of the full run's table equals the existing
                  table as text for all rows and the compare-sample report passed. Writes q002_assembly.json and
                  README.txt. Nothing in the existing base is written (hard links add a name; content and mode
                  untouched).
usage (cwd $APP, .venv/bin/python):
  nightfix_base_q002.py compare-sample --old-base BASE --new RUN_DIR [--new RUN_DIR ...] --out REPORT.json
  nightfix_base_q002.py assemble --old-base BASE --cpp-run FULL_RUN_DIR --sample-report REPORT.json --out NEW_BASE
"""
import argparse
import datetime as dt
import hashlib
import json
import os
import shutil
import socket
import sys
import time
from pathlib import Path

import numpy as np

APP = Path(os.environ.get("APP", "/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda"))
if str(APP) not in sys.path:
    sys.path.insert(0, str(APP))
GEOMETRY_COLUMNS = ("UTC", "et", "x", "y", "z", "altitude", "inside", "rh_AU", "VdotB_RG", "nVisible", "refNG", "refRG",
                    "T_RG", "useNG", "useRG")
# v3 (2026-11-18): NGBelowDeckFraction is NOT geometry: Liouville.cpp weights it by the emitted flux (st.Flux), so it
# changes with the law (v2 listed it and flagged rows 1474, 1541, 1549, 1550); reported, not gated
LAW_DEPENDENT_REPORTED = ("NGBelowDeckFraction",)
SAMPLE_GATE = 1e-8          # (b): max rel of the C++ new-law table vs fit_cops' new-law model (the check needs 1e-6)
T0 = time.time()


def say(t):
    print("nightfix_base_q002 [%6.0f s]: %s" % (time.time() - T0, t), flush=True)


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def table(path):
    """(header lines, column names, rows keyed by iPoint as token lists)."""
    hdr, names, rows = [], None, {}
    for line in open(path):
        if line.startswith("#iPoint"):
            names = line[1:].split()
        elif line.startswith("#"):
            hdr.append(line.rstrip("\n"))
        elif line.strip():
            t = line.split()
            rows[int(t[0])] = t
    return hdr, names, rows


def cmd_compare(a):
    from cops_cuda import reference
    fc = reference.load("fit_cops")
    old = Path(a.old_base)
    _, onames, orows = table(old / "COPS_model.dat")
    oc = {n: i for i, n in enumerate(onames)}
    rep = {"schema": "nightfix_base_q002_sample_v1", "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "host": socket.gethostname(), "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "command": sys.argv,
           "old_base": str(old.resolve()), "runs": [], "rows": {}}
    old_el = sha256_file(old / "elements.bin")
    sample = []
    for d in a.new:
        d = Path(d)
        _, nnames, nrows = table(d / "COPS_model.dat")
        nc = {n: i for i, n in enumerate(nnames)}
        if nnames != onames:
            raise SystemExit("%s: table columns differ" % d)
        el = sha256_file(d / "elements.bin")
        rep["runs"].append({"dir": str(d.resolve()), "rows": sorted(nrows), "elements_sha256": el,
                            "elements_equal_old": el == old_el, "manifest_sha256": sha256_file(d / "run_manifest.json")})
        for ip, t in sorted(nrows.items()):
            o = orows[ip]
            r = {"run": str(d), "UTC": t[nc["UTC"]]}
            r["geometry_columns_equal_text"] = all(t[nc[k]] == o[oc[k]] for k in GEOMETRY_COLUMNS)
            r["geometry_columns_differing"] = [k for k in GEOMETRY_COLUMNS if t[nc[k]] != o[oc[k]]]
            r["law_dependent_columns_new_old"] = {k: [t[nc[k]], o[oc[k]]] for k in LAW_DEPENDENT_REPORTED}
            kn, ko = d / "kernels" / ("point=%04d.bin" % ip), old / "kernels" / ("point=%04d.bin" % ip)
            sn, so = sha256_file(kn), sha256_file(ko)
            r["kernel_file_sha256_equal"] = sn == so
            pn, hn, recn = fc.read_kernel(str(kn))
            po, ho, reco = fc.read_kernel(str(ko))
            r["kernel_records"] = int(len(reco))
            r["kernel_fields_bitwise"] = {k: bool(len(recn) == len(reco) and
                                                  np.asarray(recn[k]).tobytes() == np.asarray(reco[k]).tobytes())
                                          for k in ("face", "Kn", "KPhi", "KV", "cos", "shadow")}
            r["kernel_header_equal"] = bool(pn == po and all(np.asarray(hn[k]).tobytes() == np.asarray(ho[k]).tobytes()
                                                             for k in ho))
            r["cpp_new_pNG"], r["cpp_new_pRG"] = float(t[nc["pNG"]]), float(t[nc["pRG"]])
            r["cpp_old_pNG"], r["cpp_old_pRG"] = float(o[oc["pNG"]]), float(o[oc["pRG"]])
            rep["rows"][str(ip)] = r
            sample.append(ip)
        say("%s: %d rows compared" % (d.name, len(nrows)))
    # (b) fit_cops' own model (fit_cops.main's consistency check) from the EXISTING kernels of the sample rows
    header, names, rows = fc.read_model_file(str(old / "COPS_model.dat"))
    cols = {n: i for i, n in enumerate(names)}
    sub = [r for r in rows if int(r[cols["iPoint"]]) in set(sample)]
    el = fc.read_elements(str(old / "elements.bin"))
    srcdir = str(Path(fc.__file__).resolve().parents[1])
    coef = fc.read_activity_coefficients(srcdir, header)
    # the constructor loads the kernels of EVERY given row (fit_cops.Model.__init__: load_kernels(arange(len(rows)))),
    # as fit_cops.main before its consistency check; 2026-11-18 fix: v1 then reloaded MF.sel (the 21 USED rows only:
    # the 8 late-grid rows of the sample have use 0) and crashed in job 136506 (IndexError)
    MF = fc.Model(str(old), header, cols, sub, el, coef)
    if MF.npt != len(sub):
        raise SystemExit("model loaded %d rows, sample has %d" % (MF.npt, len(sub)))
    pNG0, pRG0 = MF.evaluate(coef[0], MF.T0, 1.0)
    ip_sub = [int(r[cols["iPoint"]]) for r in sub]
    dn, dr, don, dor = [], [], [], []
    for k, ip in enumerate(ip_sub):
        r = rep["rows"][str(ip)]
        r["fit_cops_new_law_pNG"], r["fit_cops_new_law_pRG"] = float(pNG0[k]), float(pRG0[k])
        r["rel_cpp_new_vs_fit_cops_new"] = [abs(pNG0[k] / r["cpp_new_pNG"] - 1), abs(pRG0[k] / r["cpp_new_pRG"] - 1)]
        r["rel_cpp_old_vs_fit_cops_new"] = [abs(pNG0[k] / r["cpp_old_pNG"] - 1), abs(pRG0[k] / r["cpp_old_pRG"] - 1)]
        dn.append(r["rel_cpp_new_vs_fit_cops_new"][0]); dr.append(r["rel_cpp_new_vs_fit_cops_new"][1])
        don.append(r["rel_cpp_old_vs_fit_cops_new"][0]); dor.append(r["rel_cpp_old_vs_fit_cops_new"][1])
    rows_ok = [rep["rows"][str(ip)] for ip in sample]
    rep["summary"] = {
        "rows": len(sample), "records": int(sum(r["kernel_records"] for r in rows_ok)),
        "kernel_files_sha256_equal": int(sum(r["kernel_file_sha256_equal"] for r in rows_ok)),
        "kernel_fields_bitwise_all": bool(all(all(r["kernel_fields_bitwise"].values()) for r in rows_ok)),
        "kernel_headers_equal_all": bool(all(r["kernel_header_equal"] for r in rows_ok)),
        "geometry_columns_equal_all": bool(all(r["geometry_columns_equal_text"] for r in rows_ok)),
        "elements_equal_old_all_runs": bool(all(x["elements_equal_old"] for x in rep["runs"])),
        "max_rel_cpp_new_vs_fit_cops_new": {"NG": max(dn), "RG": max(dr)},
        "max_rel_cpp_old_vs_fit_cops_new": {"NG": max(don), "RG": max(dor)},
        "gate_b": SAMPLE_GATE, "model_records": int(len(MF.face)), "T0": [float(v) for v in MF.T0]}
    s = rep["summary"]
    s["passed"] = bool(s["kernel_fields_bitwise_all"] and s["kernel_headers_equal_all"] and s["geometry_columns_equal_all"]
                       and s["elements_equal_old_all_runs"] and max(max(dn), max(dr)) <= SAMPLE_GATE)
    rep["wall_seconds"] = time.time() - T0
    json.dump(rep, open(a.out, "x"), indent=1, default=float)
    say("summary: %s" % json.dumps(s))
    return 0 if s["passed"] else 1


def cmd_assemble(a):
    old, run, out = Path(a.old_base).resolve(), Path(a.cpp_run).resolve(), Path(a.out)
    if out.exists():
        raise SystemExit("%s exists (never overwritten)" % out)
    srep = json.load(open(a.sample_report))
    if not srep["summary"]["passed"]:
        raise SystemExit("the sample report did not pass")
    _, onames, orows = table(old / "COPS_model.dat")
    _, nnames, nrows = table(run / "COPS_model.dat")
    if onames != nnames or sorted(orows) != sorted(nrows):
        raise SystemExit("full run table: columns or rows differ from the base")
    oc = {n: i for i, n in enumerate(onames)}
    bad = [(ip, k) for ip in orows for k in GEOMETRY_COLUMNS if orows[ip][oc[k]] != nrows[ip][oc[k]]]
    if bad:
        raise SystemExit("full run table: %d geometry/data cells differ from the base (first %s)" % (len(bad), bad[:5]))
    say("full run table: %d rows, every geometry / data column equal to the base as text" % len(orows))
    mold = json.load(open(old / "run_manifest.json"))
    mrun = json.load(open(run / "run_manifest.json"))
    if mrun.get("states", {}).get("process") != "complete":
        raise SystemExit("the C++ run's manifest does not mark the process complete")
    sh_old_el = sha256_file(old / "elements.bin")
    if (run / "elements.bin").exists():
        sh_new_el = sha256_file(run / "elements.bin")
        if sh_old_el != sh_new_el:
            raise SystemExit("elements.bin of the C++ run differs from the base's")
        el_evidence = "the full C++ run's own elements.bin is byte-identical"
    else:
        # v4 (2026-11-18): a WriteKernels = off run writes no elements.bin (main.cpp:460-463 writes it with the kernels);
        # the evidence is then the sample re-runs (same binary, same input but the point range): byte-identical
        if not srep["summary"]["elements_equal_old_all_runs"] or "elements" in mrun.get("outputs", {}):
            raise SystemExit("no elements.bin in the C++ run and no sample evidence for it")
        sh_new_el = None
        el_evidence = ("the full run (WriteKernels off) wrote none; the 4 sample re-runs' elements.bin are byte-identical "
                       "to the base's (sample report)")
    stage = out.parent / ("%s.incomplete_%s" % (out.name, dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")))
    (stage / "kernels").mkdir(parents=True)
    names = sorted(p.name for p in (old / "kernels").iterdir())
    if len(names) != len(orows) or set(names) != {"point=%04d.bin" % ip for ip in orows}:
        raise SystemExit("kernel files of the base do not match its table")
    for n in names:
        os.link(old / "kernels" / n, stage / "kernels" / n)
    os.link(old / "elements.bin", stage / "elements.bin")
    for n in ("COPS_model.dat", "input.snapshot.in", "build.log"):
        shutil.copy2(run / n, stage / n)
    # v4: the C++ run's manifest records kernels version 0 / an empty set hash (it wrote none); as run_manifest.json it
    # would make cops_io.validate_fit_base refuse ("kernel set does not match"). It is kept UNCHANGED under another name;
    # the base has no run_manifest.json (validate_fit_base's legacy path: every file is still validated structurally,
    # table vs kernel headers, and the kernel set hash is computed; compare it with old_manifest_kernels.set_sha256)
    shutil.copy2(run / "run_manifest.json", stage / "cpp_run_manifest.json")
    shutil.copy2(Path(a.cpp_input), stage / "boundary.in")
    inode = {n: (os.stat(old / "kernels" / n).st_ino == os.stat(stage / "kernels" / n).st_ino) for n in names}
    rec = {"schema": "nightfix_base_q002_v1", "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
           "host": socket.gethostname(), "slurm_job_id": os.environ.get("SLURM_JOB_ID"), "command": sys.argv,
           "law": "H2O q_min = 0.02 q_max (user decision 2026-11-18); CO2 q_min = 0.1 q_max",
           "old_base": str(old), "cpp_run": str(run), "cpp_input": str(Path(a.cpp_input).resolve()),
           "hard_links": {"kernels": len(names), "kernels_same_inode_as_old": int(sum(inode.values())),
                          "elements.bin_same_inode_as_old": os.stat(old / "elements.bin").st_ino ==
                          os.stat(stage / "elements.bin").st_ino},
           "old_manifest_kernels": mold["outputs"].get("kernels"),
           "files_from_cpp_run_sha256": {n: sha256_file(stage / n) for n in ("COPS_model.dat", "cpp_run_manifest.json",
                                                                              "input.snapshot.in", "build.log",
                                                                              "boundary.in")},
           "run_manifest": "none on purpose (v4): the C++ run's manifest (cpp_run_manifest.json, unchanged) lists kernels "
                           "version 0 with an empty set hash; validate_fit_base then uses its legacy path",
           "cpp_run_manifest_outputs": mrun.get("outputs"), "elements_sha256_old": sh_old_el,
           "elements_sha256_cpp_run": sh_new_el, "elements_evidence": el_evidence, "old_COPS_model_sha256": sha256_file(old / "COPS_model.dat"),
           "geometry_columns_equal_all_rows": True, "sample_report": str(Path(a.sample_report).resolve()),
           "sample_summary": srep["summary"]}
    if rec["hard_links"]["kernels_same_inode_as_old"] != len(names) or not rec["hard_links"]["elements.bin_same_inode_as_old"]:
        raise SystemExit("hard links failed: %s" % rec["hard_links"])
    json.dump(rec, open(stage / "q002_assembly.json", "x"), indent=1, default=float)
    stage.rename(out)
    say("assembled %s" % out)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compare-sample")
    c.add_argument("--old-base", required=True)
    c.add_argument("--new", action="append", required=True)
    c.add_argument("--out", required=True)
    s = sub.add_parser("assemble")
    for k in ("old-base", "cpp-run", "cpp-input", "sample-report", "out"):
        s.add_argument("--" + k, required=True)
    a = ap.parse_args(argv)
    return cmd_compare(a) if a.cmd == "compare-sample" else cmd_assemble(a)


if __name__ == "__main__":
    sys.exit(main())
