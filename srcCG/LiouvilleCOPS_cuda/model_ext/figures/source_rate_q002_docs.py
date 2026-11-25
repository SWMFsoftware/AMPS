#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/source_rate_q002_docs.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/source_rate_q002_docs.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 592f4b23b809d667ef77b51f5e532f2b58296369413ab3fae6fd15bed1d08a23
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: README.txt, CAPTIONS.md, HOW-TO-RERUN.md and code/ of one new-law source-rate folder (2026-11-18)
# RUN: login node after the figure job(s) of the folder
# CHANGED: 2026-11-23 full set of the all-held candidate: request, date, node, submit day and command as options
# ==== end of header ====
"""README.txt, CAPTIONS.md, HOW-TO-RERUN.md and code/ of one source-rate folder of the new H2O law (2026-11-18;
figures_cluster/12M/source_rate_q002/<run>/, made by jobs/source_rate_q002.sbatch). Login node, light.

    source_rate_q002_docs.py --folder <run folder> --title "<run title>" [--note "<line>" ...]

Reads the job's provenance-<job>/inputs.txt and commands.sh, the inspection numbers (inspection/checks/inspect_*.json)
and the maps manifest of the new-law state. Refuses existing outputs.
"""
import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

S = Path(__file__).resolve().parent
APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
JOB = Path("/data001/highresolu_vt/cops_output/jobs/source_rate_q002.sbatch")
USER = "prepare the surface source rate maps for completed runs for inspection of the source rate"


def inputs(p):
    d = {}
    for line in (p / "inputs.txt").read_text().splitlines():
        k, _, v = line.partition("=")
        d[k] = v
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folder", required=True)
    ap.add_argument("--title", required=True)
    ap.add_argument("--note", action="append", default=[])
    # 2026-11-23 (source-rate set of the all-held candidate): request, its date, submit node and day, the submit command;
    # the defaults are the 2026-11-18 texts (unchanged)
    ap.add_argument("--request", default=USER)
    ap.add_argument("--request-date", default="2026-11-18, relayed by the lead")
    ap.add_argument("--node", default="gpu003")
    ap.add_argument("--submitted", default="2026-11-18 on gpu003, one job at a time, lead's conditions")
    ap.add_argument("--submit-command", default=None)
    a = ap.parse_args()
    O = Path(a.folder).resolve()
    for n in ("README.txt", "CAPTIONS.md", "HOW-TO-RERUN.md", "code"):
        if (O / n).exists():
            raise SystemExit("%s exists" % (O / n))
    provs = sorted(O.glob("provenance-*"))
    if not provs:
        raise SystemExit("no provenance-<job> folder")
    P = provs[-1]
    inp = inputs(P)
    man = json.load(open(Path(inp["MAPS"]) / "maps_manifest.json"))
    insp = {}
    for t in ("103310", "103926"):
        f = O / "inspection" / "checks" / ("inspect_%s.json" % t)
        if f.exists():
            insp[t] = json.load(open(f))
    L = [a.title, "=" * len(a.title), "",
         "User request (%s): \"%s\". Figures for inspection of the H2O source rate of a run "
         "fitted with the NEW H2O law q_min = 0.02 q_max (srcCG's value; the night-side fix of 2026-11-18, FINDINGS A18): "
         "q_H2O = q_min on night (cos SZA <= 0) and cast-shadow facets, q_min + (q_max - q_min) cos SZA on day facets. "
         "Every fit and figure before 2026-11-18 used q_min = 0.2 q_max (not intended); the old-law figures of the "
         "counterpart are in its full set." % (a.request_date, a.request), ""]
    L += ["STATE (new law)", "  maps      %s (job %s, %s)" % (inp["MAPS"], man.get("slurm_job_id"), man.get("created_utc")),
          "  fit / ft  %s / %s" % (inp["FIT"], inp["FT"]), "  base      %s" % inp["BASE"],
          "  geometry  %s" % inp["GEOM"], "  maps code %s/nonneg_maps.py" % inp["MX"],
          "OLD-LAW COUNTERPART (inspection figure only)", "  maps      %s" % inp.get("OLD_MAPS", ""),
          "  fit / ft  %s / %s" % (inp.get("OLD_FIT", ""), inp.get("OLD_FT", "")), ""]
    L += ["CONTENTS",
          "  activity/          H2O activity factor A: templates figure (+ panels/), lon/lat and 3-D ground tracks "
          "(pub_maps_v2_figures.py; seen-only colouring, area cut 0.02, flat colours)",
          "  production_rate/   H2O production rate F_H2O = q_H2O A at 10:33:10 and 10:39:26 UTC, the same layouts "
          "(relief 0.7); captions state the law (q_min = 0.02 q_max)",
          "  inspection/        per time: new-law F_H2O | old-law F_H2O (one set of classes for both) | new q/q_max "
          "(whole surface); inspection/checks/inspect_<hhmmss>.{json,txt}: numbers and checks",
          "  checks/            the v2 tools' manifests and class tables",
          "  provenance-<job>/  job script, inputs, commands.sh, logs, code copy, git, nvidia-smi, pip list", ""]
    if insp:
        L.append("NUMBERS (inspection/checks/inspect_<hhmmss>.txt; Q in s^-1; dark = night + cast shadow)")
        L.append("  %-40s %-26s %-26s" % ("", "10:33:10 UTC new (old)", "10:39:26 UTC new (old)"))
        rows = [("Q_H2O total", "Q_H2O_s-1", "%.3e"), ("Q_H2O of the seen facets", "Q_H2O_seen_s-1", "%.3e"),
                ("seen share of Q_H2O", "seen_share_of_Q", "%.3f"), ("dark share of Q_H2O", "dark_share_of_Q", "%.3f"),
                ("dark share of the seen Q_H2O", "dark_share_of_seen_Q", "%.3f"),
                ("dark area fraction (whole)", "dark_area_fraction", "%.3f"),
                ("active area A > 0 [m^2]", "active_area_m2", "%.3e"),
                ("active area fraction (whole)", "active_area_fraction", "%.4f"),
                ("active area fraction (seen)", "active_area_fraction_seen", "%.4f"),
                ("q_min / q_max", "q_ratio_H2O", "%g")]
        for lab, k, f in rows:
            cells = []
            for t in ("103310", "103926"):
                r = insp.get(t)
                cells.append(("%s (%s)" % (f % r["new"][k], f % r["old"][k])) if r else "-")
            L.append("  %-40s %-26s %-26s" % (lab, cells[0], cells[1]))
        for t, r in insp.items():
            L.append("  %s UTC: Q_H2O new / old = %.3f; seen Q_H2O new / old = %.3f" % (
                r["utc"].split("T")[1], r["new"]["Q_H2O_s-1"] / r["old"]["Q_H2O_s-1"],
                r["new"]["Q_H2O_seen_s-1"] / r["old"]["Q_H2O_seen_s-1"]))
        L.append("")
        c0 = list(insp.values())[0]["checks"]
        L += ["CHECKS",
              "  both states' fields reproduce their maps' .vtu and manifests (pub_maps_figures.load_fields: passed); "
              "the two geometries give bitwise equal cos SZA (%s) and illumination classes (%s); seen masks equal: %s; "
              "F_H2O = q A with q_min / q_max = %g (new maps) and %g (old maps)" % (
                  c0["geometry_cos_equal_bitwise"], c0["illumination_equal"], c0["seen_equal"], c0["q_ratio_new_maps"],
                  c0["q_ratio_old_maps"]), ""]
    if a.note:
        L.append("NOTES")
        L += ["  " + n for n in a.note]
        L.append("")
    L.append("JOBS (provenance-<job>/: job script, inputs, commands.sh, logs, code copy with SHA256SUMS)")
    for q in provs:
        env_ = (q / "env.txt").read_text() if (q / "env.txt").exists() else ""
        node = env_.split("SLURM_JOB_NODELIST=")[-1].split()[0] if "SLURM_JOB_NODELIST=" in env_ else "?"
        st = inputs(q).get("STEPS", "")
        tail = (q / "job.log").read_text().strip().splitlines()[-1] if (q / "job.log").exists() else ""
        L.append("  %s on %s, STEPS '%s': %s" % (q.name.split("-")[1], node, st, tail))
    L.append("Rerun: HOW-TO-RERUN.md.")
    (O / "README.txt").open("x").write("\n".join(L) + "\n")
    subprocess.run([sys.executable, str(S / "collect_captions.py"), str(O), "--title", "Figure captions: " + a.title,
                    "--order", "activity,production_rate,inspection"], check=True)
    H = ["# How to rerun: %s" % a.title, "",
         "The tools refuse existing products: move this folder aside (`mv X X.superseded_<UTC>`), never delete it.", "",
         "## Submit (GPU; %s)" % a.submitted, "", "```bash",
         "cd /data001/highresolu_vt/cops_output/jobs"]
    for k in ("KEY", "BASE", "FT", "FIT", "MAPS", "GEOM", "MX", "OLD_BASE", "OLD_FT", "OLD_FIT", "OLD_MAPS", "OLD_GEOM",
              "OLD_MX", "MODEL_NAME", "MODEL_TEXT", "OLD_NAME"):
        if inp.get(k):
            H.append("export %s=%s" % (k, json.dumps(inp[k])))
    H += ["sbatch --nodelist=%s --export=ALL source_rate_q002.sbatch" % a.node]
    if a.submit_command:
        H += ["# (as submitted: %s)" % a.submit_command]
    H += [
          "# then on the login node:",
          "%s/.venv/bin/python %s --folder %s --title %s%s%s" % (APP, Path(__file__).resolve(), O, json.dumps(a.title),
                                                               "".join(" --note %s" % json.dumps(n) for n in a.note),
                                                               "".join(" --%s %s" % (k.replace("_", "-"), json.dumps(v)) for k, v in (
                                                                   ("request", a.request), ("request_date", a.request_date),
                                                                   ("node", a.node), ("submitted", a.submitted),
                                                                   ("submit_command", a.submit_command)) if v is not None
                                                                   and a.request_date != "2026-11-18, relayed by the lead")),
          "```", ""]
    for q in provs:
        H += ["## The exact commands as run (%s/commands.sh)" % q.name, "", "```bash", (q / "commands.sh").read_text().rstrip(),
              "```", ""]
    (O / "HOW-TO-RERUN.md").open("x").write("\n".join(H))
    code = O / "code"
    code.mkdir()
    for f in (S / "pub_maps_v2_figures.py", S / "pub_maps_figures.py", S / "illumination_maps_figures.py",
              S / "source_rate_q002_inspect.py", Path(__file__).resolve(), S / "collect_captions.py", JOB):
        shutil.copy2(f, code / f.name)
    with open(code / "SHA256SUMS", "w") as fh:
        for f in sorted(code.iterdir()):
            if f.name != "SHA256SUMS":
                fh.write("%s  %s\n" % (hashlib.sha256(f.read_bytes()).hexdigest(), f.name))
    print("wrote README.txt, CAPTIONS.md, HOW-TO-RERUN.md, code/ in %s" % O)
    return 0


if __name__ == "__main__":
    sys.exit(main())
