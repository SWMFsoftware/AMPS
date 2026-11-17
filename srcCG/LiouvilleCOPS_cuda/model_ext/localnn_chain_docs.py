#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-15, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/model_ext/localnn_chain_docs.py (on Voyager: /data001/highresolu_vt/cops_output/model_ext/localnn_chain_docs.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: efdf760e95e6cfd6ba123d6653effdd49107a94807a3730ef26df781593246da
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: README.md / HOW-TO-RERUN.md of one stage-2 chain from its files (submission.json, chain_events.jsonl, each
#   cycle's fit meta and fine_tuning.json); earlier versions kept as *.before_<UTC>; the README text names the
#   templates setting, the nonneg code and the fine-tuning regularization (2026-11-16); HOW-TO-RERUN names the
#   continuation (--chain) and the GPU move.
# RUN 2026-11-15: localnn_chain_submit.sh, the chain post jobs.
# ==== end of header ====
"""README.md / HOW-TO-RERUN.md of ONE option C stage-2 chain (jobs/localnn_chain.sbatch layout), from its files only.

    localnn_chain_docs.py --chain CHAIN

Reads CHAIN/submission.json (written by jobs/localnn_chain_submit.sh: settings, the submit command, the sbatch lines and
job ids), CHAIN/chain_events.jsonl (one line per stage: stage, rc, device, host, job, UTC), every cycle's
c<k>/fit/final/localnn_state.npz meta (costs, evaluations, termination) and c<k>/ft/final/fine_tuning.json (fine tuning:
costs, evaluations, stop) and writes CHAIN/README.md and CHAIN/HOW-TO-RERUN.md (earlier versions are kept as
*.before_<UTC>). It changes nothing else. Light (login node is fine): only small npz members and json are read.
"""

import argparse
import json
import time
from pathlib import Path

import numpy as np

UTC = lambda: time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())  # noqa: E731


def fit_info(fit):
    f = fit / "final" / "localnn_state.npz"
    if not f.exists():
        segs = sorted(p.name for p in fit.glob("seg[0-9][0-9]")) if fit.exists() else []
        ck = fit / "checkpoint.npz"
        r = {"done": False, "segments": segs}
        if ck.exists():
            with np.load(ck, allow_pickle=False) as z:
                r["checkpoint_nfev_total"] = int(z["nfev_total"])
        return r
    with np.load(f, allow_pickle=False) as z:
        m = json.loads(str(z["meta_json"]))
    c = m.get("costs", {})
    segs = sorted(p.name for p in fit.glob("seg[0-9][0-9]"))
    return {"done": True, "final": str((fit / "final").resolve()), "segments": segs, "nfev_total": m.get("nfev_total"),
            "message": m.get("message"), "objective": c.get("objective"), "data": c.get("data_cost"),
            "crossing": c.get("crossing_cost"), "penalty": c.get("penalty"), "mult_penalty": c.get("mult_penalty"),
            "held": (m.get("held") or {}).get("nonzero_facets"), "resumed": bool(m.get("resumed_from_checkpoint")),
            "device": m.get("device", "cpu")}


def ft_info(ft):
    f = ft / "final" / "fine_tuning.json"
    if not f.exists():
        return {"done": False}
    j = json.load(open(f))
    o = j.get("optimizer", {})
    ci, cf = j.get("costs_initial", {}), j.get("costs_final", {})
    pick = lambda d: d.get("total", d.get("cost_total")) if isinstance(d, dict) else d  # noqa: E731
    return {"done": True, "final": str((ft / "final").resolve()), "nfev": o.get("nfev"),
            "termination": o.get("termination"), "stop": j.get("stop"), "initial": pick(ci), "final_cost": pick(cf),
            "host_vs_gpu": j.get("final_host_vs_gpu_max_rel")}


def fmt(v, f="%.10e"):
    return "-" if v is None else (f % v if isinstance(v, (int, float)) else str(v))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--chain", required=True)
    args = ap.parse_args(argv)
    ch = Path(args.chain).resolve()
    sub = json.load(open(ch / "submission.json")) if (ch / "submission.json").exists() else {}
    st = sub.get("settings", {})
    ev = []
    if (ch / "chain_events.jsonl").exists():
        for line in open(ch / "chain_events.jsonl"):
            try:
                ev.append(json.loads(line))
            except ValueError:
                pass
    ncyc = int(st.get("CYCLES", 5))
    rows = []
    for k in range(1, ncyc + 1):
        rows.append((k, fit_info(ch / ("c%d" % k) / "fit"), ft_info(ch / ("c%d" % k) / "ft")))
    complete = any(e.get("stage") == "chain_complete" and str(e.get("rc")) == "0" for e in ev)
    L = ["# Option C stage 2 chain: %s" % ch.name, "",
         "Local non-negative basis (option C), A_e = m_b(e) (sum_j c_j phi_j(e) + A_s TS_e) exp(h_e); %d cycles of "
         "fit (localnn_fit.py, max_nfev %s, smoothness lambda_C %s, multiplier penalty f_m %s, templates %s, "
         "slope-mult-min 1) + CPU fine tuning (localnn_ft.py: the frozen nonneg_finetune %s, late_only facets, "
         "WINDOW_MIN %s, grain refitted, lambda 0, MAXNFEV %s; %s); the accumulated "
         "h = held_next of cycle k is held in the fit of cycle k+1." % (
             ncyc, st.get("NFEV"), st.get("LAMBDA"), st.get("FM"), st.get("TEMPLATES", "on"),
             st.get("NONNEG_CODE", "model_ext_frozen/nonneg_20261113T035423Z"), st.get("WINDOW_MIN"),
             st.get("FT_MAXNFEV"), ("regularization " + st["FT_EXTRA"]) if st.get("FT_EXTRA") else
             "no lower bound on delta, no regularization"), "",
         "Status (%s): %s." % (time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime()),
                               "chain complete" if complete else "in progress"), ""]
    if sub.get("c1_note"):
        L += ["Cycle 1: " + sub["c1_note"], ""]
    L += ["| cycle | fit evals | fit termination | objective | data | P_m | held facets | ft evals | ft stop | "
          "ft final cost |", "|---|---|---|---|---|---|---|---|---|---|"]
    for k, fi, ti in rows:
        if fi["done"]:
            fcol = [fmt(fi["nfev_total"], "%d"), (fi["message"] or "-").replace("|", "/"), fmt(fi["objective"]),
                    fmt(fi["data"]), fmt(fi["mult_penalty"], "%.4e"), fmt(fi["held"], "%d")]
        else:
            fcol = ["(running: %s, checkpoint at %s)" % (",".join(fi["segments"]) or "-",
                                                        fi.get("checkpoint_nfev_total", "-"))] + ["-"] * 5
        tcol = ([fmt(ti["nfev"], "%d"), "%s (%s)" % (ti["stop"], ti["termination"]), fmt(ti["final_cost"])]
                if ti["done"] else ["-", "-", "-"])
        L.append("| c%d | " % k + " | ".join(fcol + tcol) + " |")
    L += ["", "Reaching max_nfev is not convergence. Costs: objective = data + crossing + smoothness + P_m (fit); "
          "the fine tuning's cost is data + crossing at lambda 0 (no smoothness, no P_m).", "",
          "## Stages (chain_events.jsonl: device, host, job)", "", "| stage | rc | device | host | job | UTC | note |",
          "|---|---|---|---|---|---|---|"]
    for e in ev:
        L.append("| %s | %s | %s | %s | %s | %s | %s |" % tuple(str(e.get(k, "")) for k in
                                                              ("stage", "rc", "device", "host", "job", "utc", "note")))
    L += ["", "## Layout", "",
          "* `c<k>/fit`: localnn_fit.py (seg01, seg02, ... = wall-time segments; `final` -> the last; "
          "`checkpoint.npz` while running, moved into the segment as checkpoint_last.npz; resumed checkpoints kept as "
          "checkpoint_resumed_segNN_<UTC>.npz); `direct_check_segNN` = the direct facet-sum check.",
          "* `c<k>/ft`: localnn_ft.py run (seg01, `final` -> seg01; links to its finetune_state.npz, "
          "finetune_result.txt, fine_tuning.json, late_window_model.dat, pressure_comparison.txt, validation.json, "
          "localnn_ft.json).",
          "* `provenance-<job>/`: per attempt: code copies + sha256, the job script used, git, modules, nvidia-smi, "
          "toolchain, pip list, `commands.sh` (every command with all arguments), the stage logs.",
          "* After completion: `c%d/analysis` (localnn_report.py), `c%d/analysis_degeneracy`, `c%d/tools_state_ft`, "
          "`c%d/contributions_ft`, `c%d/maps_ft` (the final state = the last fine tuning), `c%d/gpu_eval_ft` (GPU "
          "cross-check, separate job)." % ((ncyc,) * 6), ""]
    if sub:
        L += ["Submission: see HOW-TO-RERUN.md (jobs %s)." % ", ".join(str(j) for j in sub.get("jobs", []))]
    H = ["# How to rerun: %s" % ch.name, ""]
    if sub:
        H += ["Settings: " + ", ".join("%s=%s" % kv for kv in sorted(st.items())), "",
              "Code: MODEL_EXT `%s`; AMPS %s." % (st.get("MODEL_EXT"), sub.get("amps_commit")), "",
              "## Submit command (%s)" % sub.get("submitted_utc"), "", "```bash", sub.get("command", ""), "```", "",
              "The sbatch lines it ran (job ids):", "", "```"]
        for s_, j in zip(sub.get("sbatch_lines", []), sub.get("jobs", [])):
            H.append("# job %s" % j)
            H.append(s_)
        H += ["```", ""]
    H += ["To continue after all attempts ended with exit 75 (wall time), submit the same command again with "
          "`--chain <this root>` added: completed stages are skipped, a running fit continues from its segment / "
          "checkpoint. To move the remaining cycles to a GPU: `jobs/localnn_chain_to_gpu.sh --chain <this root> "
          "--model-ext <frozen GPU snapshot> --submit`. Each attempt's `provenance-<job>/commands.sh` lists every "
          "command with all arguments resolved (start, fit, direct check, fine tuning)."]
    for name, text in (("README.md", L), ("HOW-TO-RERUN.md", H)):
        p = ch / name
        if p.exists():
            p.rename(ch / ("%s.before_%s" % (name, UTC())))
        p.write_text("\n".join(text) + "\n")
    print("localnn_chain_docs: wrote %s/README.md, HOW-TO-RERUN.md (%s)" % (ch, "complete" if complete else
                                                                           "in progress"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
