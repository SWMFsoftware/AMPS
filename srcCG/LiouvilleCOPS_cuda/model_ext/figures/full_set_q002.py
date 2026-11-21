#!/usr/bin/env python3
# ==== version-controlled copy: header added 2026-11-20, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/full_set_q002.py (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/full_set_q002.py)
# SHA256 of the run copy = SHA256 of this file without these header lines: 3619181a95904ef96dfdf96153d919411ed944f960ac2540d0ce3aad291ac374
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: full figure sets of the ten new-law (H2O q_min = 0.02 q_max) states: spec, publication figures and copies of
#   the source-rate products (prepare), submission (full_set_submit.sh), finish and docs (2026-11-20)
# RUN: full_set_q002.py spec | prepare (srun CPU) | submit | finish (jobs/full_set_q002_finish.sbatch) | docs
# ==== end of header ====
"""Full figure sets of the ten NEW-LAW states (H2O q_min = 0.02 q_max; 2026-11-20, lead; user: "prepare the full set of
plots for the last round of fits with the fixed surface source distribution"). The sets are made exactly like the
2026-11-17 full sets (full_set_submit.sh, finish_full_set.sh, full_set_docs.py; the job scripts with the 2026-11-20
BASE / D2REF / MCODE / MODEL_TEXT / MEM_CAP / SKIP overrides), on base_q002 with the new-law snapshots
(nonneg_maps / nonneg_contrib_20261118T192126Z) and each state's own maps / contributions.

    full_set_q002.py spec                       write data/full_sets_q002_20261120.json (refuses an existing file)
    full_set_q002.py prepare [--keys k1,k2]     CPU (srun): publication NG/RG figure folder publication_<key>/ and the
                                                products copied from source_rate_q002/<run>/ (activity/, production_rate/,
                                                inspection/, their checks/ and provenance; COPIED_FROM_source_rate_q002.txt)
    full_set_q002.py submit  [--keys ...] [--dry-run]   full_set_submit.sh per set (gpu002, MEM_CAP 24G, SKIP activity,h2o)
    full_set_q002.py finish  --key k --job-set4 J        finish_full_set.sh (run inside an srun CPU allocation)
    full_set_q002.py docs    --key k                     full_set_docs.py --spec <spec> --key k
Inputs of every state from its maps manifest (base, geometry, ft / fit folders, maps code). Keys: q002_D2, q002_softplus,
q002_A, q002_B, q002_D, q002_C_f1_pm1em4, q002_C_f10_pm1em4, q002_C_f1_pm0p01, q002_C_f10_pm0p01, q002_C_f10_pm0p01_tploff.
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path("/data001/highresolu_vt/cops_output")
R = ROOT / "runs/shap8_012m_lategrid_dz1.5"
NN = R / "nonneg"
F12 = ROOT / "figures_cluster/12M"
S = ROOT / "figures_cluster/scripts"
J = ROOT / "jobs"
APP = Path("/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda")
PY = str(APP / ".venv/bin/python")
SPEC = S / "data/full_sets_q002_20261120.json"
RUNS = S / "data/source_rate_q002_runs_20261118.json"
FZ = ROOT / "model_ext_frozen"
MX = FZ / "nonneg_maps_20261118T192126Z"
CX = FZ / "nonneg_contrib_20261118T192126Z"
LNN = FZ / "localnn_q002_20261118T162327Z/localnn_tools.py"
D2T = FZ / "d2_tools_20261117T225752Z/d2_tools.py"
BASE = R / "base_q002"
D2REF = R / "clip_refit/D2_q002_warm/final"
LAW = ("new H2O law q_min = 0.02 q_max (srcCG's value; the night-side fix of 2026-11-18); the earlier full sets used "
       "q_min = 0.2 q_max")
INTRO = ("Lead 2026-11-20 14:45 CST (user: \"prepare the full set of plots for the last round of fits with the fixed surface "
         "source distribution\"): a FULL figure set of a NEW-LAW state (H2O q_min = 0.02 q_max, the night-side fix of "
         "2026-11-18; base_q002, new-law snapshots nonneg_maps / nonneg_contrib_20261118T192126Z), exactly like the "
         "2026-11-17 full sets (../full_set_stage2_*/, ../full_set_D2/; same tools, styles and user decisions), no movies. "
         "activity/, production_rate/ and inspection/ are copies (cp -p, sha256 in COPIED_FROM_source_rate_q002.txt) of "
         "../source_rate_q002/<run>/ (the same tool, settings, state, maps, geometry and code); every other folder was made "
         "for this set.")

# key: (source-rate run, model text key, root (relative to NN, or D2), kind, title phrase, model name phrase, label, caveats,
#       old-law full set)
RUNLIST = [
    ("q002_D2", "D2_q002_warm", "clipped", "D2", "d2", "the clipped activity model D2 refitted",
     "the clipped activity model D2 refitted with the new law (D2_q002_warm: warm start from D2, one segment of 4000 "
     "evaluations, max_nfev, not converged)", "clipped D2 refitted", "last fit max_nfev (4000): not converged|no slope "
     "multipliers (clipped linear activity, no fine tuning)|kappa at its lower bound (CO2 negligible)", "full_set_D2"),
    ("q002_softplus", "softplus_q002_reg", "softplus", "softplus/stage2r_q002_smm1_20261118T193956Z", "nonneg",
     "the softplus activity model", "the softplus activity model (slope multipliers >= 1) refitted with the new law, "
     "regularized stage-2 chain, final state after 5 fit and fine-tuning cycles", "softplus model",
     "last fit ftol (591)|kappa at its lower bound (CO2 negligible)", "full_set_stage2_softplus_reg"),
    ("q002_A", "A_q002_reg", "sos", "sos/stage2r_q002_r3_20261118T193956Z", "nonneg",
     "the option-A sum-of-squares model (rank 3)", "the option-A sum-of-squares model (rank 3, slope multipliers >= 1) "
     "refitted with the new law, regularized stage-2 chain, final state after 5 fit and fine-tuning cycles",
     "sum-of-squares model (A)", "last fit max_nfev (4000): not converged|kappa at its lower bound (CO2 negligible)",
     "full_set_stage2_A_reg"),
    ("q002_B", "B_q002_reg", "background", "background/stage2r_q002_bg1e-2_20261118T193956Z", "nonneg",
     "the option-B background model (A_bg >= 1e-2)", "the option-B background model (A_bg >= 1e-2, slope multipliers "
     ">= 1) refitted with the new law, regularized stage-2 chain, final state after 5 fit and fine-tuning cycles",
     "background model (B)", "last fit max_nfev (4000): not converged|kappa and A_bg at their lower bounds",
     "full_set_stage2_B_reg"),
    ("q002_D", "D_q002_reg", "logspread", "logspread/stage2r_q002_f0.1_20261118T193956Z", "nonneg",
     "the option-D log-spread model (f 0.1)", "the option-D log-spread model (f 0.1, slope multipliers >= 1) refitted "
     "with the new law, regularized stage-2 chain, final state after 5 fit and fine-tuning cycles", "log-spread model (D)",
     "last fit ftol (31)|kappa at its lower bound (CO2 negligible)", "full_set_stage2_D_reg"),
    ("q002_C_f1_pm1em4", "C_f1_pm1em4_q002_reg", "localnn_tplon", "localnn/stage2r_q002_f1_pm1e-4_tplon_20261118T162407Z",
     "localnn", "the option-C local non-negative basis (f 1, f_m 1e-4, templates on)", "the option-C local non-negative "
     "basis model (f 1, multiplier penalty f_m 1e-4, templates on, slope multipliers >= 1) refitted with the new law, "
     "regularized stage-2 chain, final state after 5 fit and fine-tuning cycles", "local basis C (f 1, f_m 1e-4)",
     "last fit max_nfev (4000): not converged|slope multipliers not determined at f_m 1e-4 (70-80 / 80-90 deg 2.5e4 / "
     "8.1e4)|kappa at its lower bound", "full_set_stage2_C_f1_pm1em4_reg"),
    ("q002_C_f10_pm1em4", "C_f10_pm1em4_q002_reg", "localnn_tplon",
     "localnn/stage2r_q002_f10_pm1e-4_tplon_20261118T162407Z", "localnn",
     "the option-C local non-negative basis (f 10, f_m 1e-4, templates on)", "the option-C local non-negative basis model "
     "(f 10, multiplier penalty f_m 1e-4, templates on, slope multipliers >= 1) refitted with the new law, regularized "
     "stage-2 chain, final state after 5 fit and fine-tuning cycles", "local basis C (f 10, f_m 1e-4)",
     "last fit ftol (812, segment 2)|slope multipliers not determined at f_m 1e-4 (70-80 / 80-90 deg 5.5e4 / 1.8e5)|kappa "
     "at its lower bound", None),
    ("q002_C_f1_pm0p01", "C_f1_pm0p01_q002_reg", "localnn_tplon", "localnn/stage2r_q002_f1_pm0.01_tplon_20261118T162407Z",
     "localnn", "the option-C local non-negative basis (f 1, f_m 0.01, templates on)", "the option-C local non-negative "
     "basis model (f 1, multiplier penalty f_m 0.01, templates on, slope multipliers >= 1) refitted with the new law, "
     "regularized stage-2 chain, final state after 5 fit and fine-tuning cycles", "local basis C (f 1, f_m 0.01)",
     "last fit ftol (7)|kappa at its lower bound", None),
    ("q002_C_f10_pm0p01", "C_f10_pm0p01_q002_reg", "localnn_tplon",
     "localnn/stage2r_q002_f10_pm0.01_tplon_20261118T162407Z", "localnn",
     "the option-C local non-negative basis (f 10, f_m 0.01, templates on)", "the option-C local non-negative basis model "
     "(f 10, multiplier penalty f_m 0.01, templates on, slope multipliers >= 1) refitted with the new law, regularized "
     "stage-2 chain, final state after 5 fit and fine-tuning cycles", "local basis C (f 10, f_m 0.01)",
     "last fit ftol (6)|kappa at its lower bound", None),
    ("q002_C_f10_pm0p01_tploff", "C_f10_pm0p01_tploff_q002_reg", "localnn_tploff",
     "localnn/stage2r_q002_f10_pm0.01_tploff_20261118T162407Z", "localnn",
     "the option-C local non-negative basis (f 10, f_m 0.01, templates off)", "the option-C local non-negative basis model "
     "(f 10, multiplier penalty f_m 0.01, templates off, slope multipliers >= 1) refitted with the new law, regularized "
     "stage-2 chain, final state after 5 fit and fine-tuning cycles", "local basis C (f 10, f_m 0.01, templates off)",
     "last fit ftol (3)|kappa at its lower bound", None),
]


def maps_state(maps):
    m = json.load(open(Path(maps) / "maps_manifest.json"))
    if m["code_sha256"]["nonneg_maps.py"] != hashlib.sha256((MX / "nonneg_maps.py").read_bytes()).hexdigest():
        raise SystemExit("%s: maps code is not %s" % (maps, MX))
    if Path(m["base"]) != BASE:
        raise SystemExit("%s: base %s, not base_q002" % (maps, m["base"]))
    return m


def build_spec():
    texts = json.load(open(RUNS))["model_texts"]
    sets = {}
    for key, srk, mt, root, kind, title, name, label, cav, old in RUNLIST:
        if kind == "d2":
            chain = R / "clip_refit/D2_q002_warm"
            ft = fit = chain / "tools_state"
            maps, contrib = chain / "maps", chain / "contributions"
            xc, phase = chain / "v2_total_crosscheck", chain / "contrib_phase"
            adapter, chain_env, ncyc = D2T, str(chain), None
            state = "%s (tools: %s)" % (chain / "final", ft)
        elif kind == "localnn":
            chain = NN / root / "c5"
            ft, fit = chain / "tools_state_ft", chain / "tools_state_ft" / "fit_state"
            maps, contrib = chain / "maps_ft", chain / "contributions_ft"
            xc, phase = chain / "v2_total_crosscheck_ft", chain / "contrib_phase_ft"
            adapter, chain_env, ncyc = LNN, str(chain), 5
            state = "%s (tools: %s)" % (chain / "ft/final", ft)
        else:
            chain = NN / root
            ft, fit = chain / "c5/ft/final", chain / "c5/fit/final"
            maps, contrib = chain / "maps_c5", chain / "contributions_c5"
            xc, phase = chain / "v2_total_crosscheck_c5", chain / "contrib_phase_c5"
            adapter, chain_env, ncyc = None, str(chain), 5
            state = str(ft)
        m = maps_state(maps)
        geom = m["geometry"]
        pub = "publication_%s" % key
        mname = "%s, %s" % (name, LAW)
        env = {"KEY": key, "CHAIN": chain_env, "FT": str(ft), "FIT": str(fit), "MAPS": str(maps), "GEOM": geom,
               "CONTRIB": str(contrib), "XC": str(xc), "PHASE": str(phase), "MX": str(MX), "CX": str(CX),
               "ADAPTER": str(adapter) if adapter else "", "PUB": str(F12 / pub), "PUBSTEM": "12M_%s_NG_RG_altitude" % key,
               "PUBTABLE": "%s_final.dat" % key, "MODEL_NAME": mname,
               "MODEL_LABEL": "%s, new law q_min 0.02 q_max, final state" % label,
               "MODEL_SHORT": "%s, new law" % label,
               "TITLE_PRESSURES": "SHAP8 12M: %s, new H2O law (q_min 0.02 q_max) vs COPS" % label,
               "MODEL_TEXT": texts[mt] or "", "BASE": str(BASE), "D2REF": str(D2REF), "MCODE": str(MX / "nonneg_maps.py"),
               "MEM_CAP": "24G", "SKIP": "activity,h2o", "NODE": "gpu002"}
        if ncyc:
            env["NCYC"] = str(ncyc)
        fin = {"KEY": key, "FT": str(ft), "FIT": str(fit), "MAPS": str(maps), "CONTRIB": str(contrib), "PHASE": str(phase),
               "MODEL_NAME": mname, "BASE": str(BASE), "REF_DIR": str(D2REF), "CAVEATS": cav}
        if kind == "localnn":
            fin["CHAIN"], fin["NCYC"] = str(chain / ".."), "5"
        elif kind == "nonneg":
            fin["CHAIN"], fin["NCYC"] = str(chain), "5"
        notes = ["Law: %s." % LAW,
                 "State inputs from its maps manifest: maps %s (job %s), geometry %s (content identical to "
                 "maps_geometry/job_136691/geometry, checked: all cos / shadow / mesh files)." % (maps, m.get("slurm_job_id"), geom),
                 "activity/, production_rate/, inspection/ and their checks/ files copied from ../source_rate_q002/%s/ "
                 "(COPIED_FROM_source_rate_q002.txt); the jobs that made them: the copied provenance-<job>/ folders." % srk,
                 "Model-data reference (model_data/, regions): the refitted D2 with the same law (%s)." % D2REF,
                 "Caveats: %s." % cav.replace("|", "; ")]
        if old:
            notes.append("The same state with the old law (q_min = 0.2 q_max): ../%s/." % old)
        if kind == "d2":
            notes.append("contributions: d2_tools.sbatch STEPS='contrib verify' job 137493 (2026-11-20; verify: NG model vs "
                         "D2's per-row model 5.1e-12, RG 2.2e-15, data cost 5.6707116e-3 = fit_result to its print "
                         "precision).")
        sets[key] = {
            "title": "Complete figure set (12M) of %s, new H2O law (q_min = 0.02 q_max), final state" % title,
            "model_name": mname, "state": state, "chain": chain_env if kind != "d2" else None, "ft": str(ft),
            "fit": str(fit), "maps": str(maps), "contrib": str(contrib), "phase": str(phase), "xc": str(xc), "geom": geom,
            "srkey": srk, "jobs_state": ["maps job %s" % m.get("slurm_job_id")] + (["contributions 137493 (d2_tools.sbatch)"]
                                                                                    if kind == "d2" else []),
            "frozen": {"maps / contributions": "model_ext_frozen/nonneg_maps_20261118T192126Z, nonneg_contrib_20261118T192126Z",
                       **({"adapter": str(adapter)} if adapter else {})},
            "notes": notes, "comparisons": ["vs the old-law full set of the same model (../%s/)" % old] if old else [],
            "intro": INTRO, "reference": "the refitted D2 with the same law (reference)", "pub": pub,
            "written": "2026-11-20 (CST)",
            "extra_dirs": ["inspection/        copy of ../source_rate_q002/%s/inspection/ (new-law vs old-law production "
                           "rate, q/q_max)" % srk if (F12 / "source_rate_q002" / srk / "inspection").exists() else
                           "inspection/        not made (no old-law maps of this setting)"],
            "env": env, "finish_env": fin,
            "submit_command": "full_set_q002.py submit --keys %s (the environment: env of this spec entry; bash "
                              "full_set_submit.sh)" % key,
            "finish_command": "srun -p normal128cores -A general -N1 -n1 -c 8 --mem=96G -t 00:45:00 <full_set_q002.py "
                              "finish --key %s --job-set4 <set1/set4 job>>" % key,
            "state_commands": (["STEPS='contrib verify' MODEL_EXT=%s D2DIR=%s D2_BASE=%s GEOM=%s CONTRIB=%s sbatch "
                                "--export=ALL d2_tools.sbatch" % (D2T.parent, chain, BASE, geom, CX)] if kind == "d2"
                               else ["# maps / contributions: the chain's post (its README.md / HOW-TO-RERUN.md)"])}
    return sets


def sh(cmd, env=None, cwd=None):
    print("$ " + " ".join(cmd), flush=True)
    return subprocess.run(cmd, env=env, cwd=cwd, check=True)


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def prepare(spec, keys):
    for k in keys:
        sp = spec["sets"][k]
        e = sp["env"]
        pub = Path(e["PUB"])
        if not pub.exists():
            pub.mkdir()
            ad = sp["env"]["ADAPTER"]                  # the state's own COPS_model.dat (as the 2026-11-17 publication figures)
            ftdir = str(D2REF) if ad.endswith("d2_tools.py") else (str(Path(e["CHAIN"]) / "ft" / "final") if ad else e["FT"])
            sh([PY, "-W", "ignore", str(S / "positive_table.py"), ftdir, str(pub / e["PUBTABLE"])], cwd=str(APP))
            sh([PY, "-W", "ignore", str(S / "pub_ng_rg_altitude.py"), "--table", str(pub / e["PUBTABLE"]), "--out",
                str(pub / e["PUBSTEM"])], cwd=str(APP), env=dict(os.environ, MPLBACKEND="Agg"))
            (pub / "code").mkdir()
            for f in ("positive_table.py", "pub_ng_rg_altitude.py", "full_set_q002.py"):
                shutil.copy2(S / f, pub / "code" / f)
            subprocess.run(["bash", "-c", "cd '%s' && sha256sum * > SHA256SUMS" % (pub / "code")], check=True)
            (pub / "HOW-TO-RERUN.md").write_text(
                "# How to rerun: publication figure, COPS and %s NG / RG vs altitude\n\nWritten 2026-11-20 (CST) for "
                "../full_set_%s/ (lead 2026-11-20: full sets of the new-law states). Same figure style and scripts as "
                "../publication_stage2_softplus_reg/ (unchanged). Law: %s. State: %s.\n\n```bash\ncd %s\n%s -W ignore %s "
                "%s %s\nMPLBACKEND=Agg %s -W ignore %s --table %s --out %s\n```\n\nCode: code/ (SHA256SUMS).\n" % (
                    sp["model_name"].split(",")[0], k, LAW, sp["state"], APP, PY, S / "positive_table.py", ftdir,
                    pub / e["PUBTABLE"], PY, S / "pub_ng_rg_altitude.py", pub / e["PUBTABLE"], pub / e["PUBSTEM"]))
        O = F12 / ("full_set_%s" % k)
        src = F12 / "source_rate_q002" / sp["srkey"]
        marker = O / "COPIED_FROM_source_rate_q002.txt"
        if marker.exists():
            continue
        O.mkdir(exist_ok=True)
        (O / "checks").mkdir(exist_ok=True)
        lines = ["# copied (cp -p, not regenerated) on the login / CPU node by full_set_q002.py prepare from %s" % src,
                 "# the same tool (pub_maps_v2_figures.py), settings (LOWCUT 0.02, classes 15-20, relief 0.7 rate / 0 "
                 "activity), state, maps, geometry, maps code and caption texts as a full-set job would use",
                 "# sha256 destination = source (checked) | path"]
        items = [d for d in ("activity", "production_rate", "inspection") if (src / d).exists()]
        provs = [p.name for p in sorted(src.glob("provenance-*")) if p.name != "provenance-136740"]
        for d in items + provs:
            if (O / d).exists():
                raise SystemExit("%s exists" % (O / d))
            shutil.copytree(src / d, O / d, copy_function=shutil.copy2)
        for f in sorted((src / "checks").iterdir()):
            if (O / "checks" / f.name).exists():
                raise SystemExit("%s exists" % (O / "checks" / f.name))
            shutil.copy2(f, O / "checks" / f.name)
        n = 0
        for d in items + provs + ["checks"]:
            for f in sorted((O / d).rglob("*")):
                if f.is_file():
                    rel = f.relative_to(O)
                    if d == "checks" and not (src / rel).exists():
                        continue
                    a, b = sha(f), sha(src / rel)
                    if a != b:
                        raise SystemExit("copy differs: %s" % rel)
                    lines.append("%s  %s" % (a, rel))
                    n += 1
        lines.append("# %d files; source provenance folders copied: %s (provenance-136740 of D2_q002_warm not copied: its "
                     "inspection was superseded)" % (n, ", ".join(provs)))
        marker.write_text("\n".join(lines) + "\n")
        print("%s: copied %d files from %s" % (k, n, src))


def submit(spec, keys, dry):
    for k in keys:
        e = dict(os.environ)
        e.update({x: y for x, y in spec["sets"][k]["env"].items()})
        if k == "q002_D2":
            e["CONTRIB_JOB"] = ""
        cmd = ["bash", str(J / "full_set_submit.sh")]
        print("# %s: %s" % (k, " ".join("%s=%s" % (x, json.dumps(y)) for x, y in spec["sets"][k]["env"].items())))
        if not dry:
            sh(cmd, env=e, cwd=str(J))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=("spec", "prepare", "submit", "finish", "docs"))
    ap.add_argument("--keys", default=None)
    ap.add_argument("--key", default=None)
    ap.add_argument("--job-set4", default=None)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()
    if a.cmd == "spec":
        if SPEC.exists():
            raise SystemExit("%s exists" % SPEC)
        sets = build_spec()
        SPEC.write_text(json.dumps({"_about": "Full figure sets of the ten new-law states (lead 2026-11-20 14:45 CST); "
                                              "full_set_q002.py; full_set_docs.py --spec this --key <key>",
                                    "sets": sets}, indent=1) + "\n")
        print("wrote %s (%d sets)" % (SPEC, len(sets)))
        return 0
    spec = json.load(open(SPEC))
    keys = a.keys.split(",") if a.keys else list(spec["sets"])
    if a.cmd == "prepare":
        prepare(spec, keys)
    elif a.cmd == "submit":
        submit(spec, keys, a.dry_run)
    elif a.cmd == "finish":
        e = dict(os.environ)
        e.update(spec["sets"][a.key]["finish_env"])
        e["JOB_SET4"] = a.job_set4
        sh(["bash", str(J / "finish_full_set.sh")], env=e)
    elif a.cmd == "docs":
        sh([PY, str(S / "full_set_docs.py"), "--spec", str(SPEC), "--key", a.key])
    return 0


if __name__ == "__main__":
    sys.exit(main())
