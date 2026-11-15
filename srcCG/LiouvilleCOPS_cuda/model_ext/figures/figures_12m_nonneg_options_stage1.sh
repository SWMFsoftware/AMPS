#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_12m_nonneg_options_stage1.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_12m_nonneg_options_stage1.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: aacb334a1ef371df62176cc17346252acbb4228530ef5a1ed7ccf9f192f999ad
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: stage-1 comparison of the non-negative activity options (A sum of squares, B background, C local basis, D log
#   spread) with the COPS data, the bounded softplus c9 and D2: tables, per-option and best-of-each model/data figures,
#   captions (2026-11-14); v2 options DIAG, RESID, ACT_FROM; final options RUNS (candidate list file), TEXTS (caption
#   texts), CONV_TERM (convergence from the terminations) (2026-11-15); defaults unchanged
# RUN: OUT=<new folder> [RUNS=<json> TEXTS=<json> CONV_TERM=1 | DIAG=<json>] [RESID=1] [ACT_FROM=<folder>] srun -p
#   normal128cores -A general -N1 -n1 -c 8 --mem=64G -t 00:50:00 bash figures_12m_nonneg_options_stage1.sh (CPU); then
#   jobs/options_activity_maps.sbatch unless ACT_FROM
# ==== end of header ====
# Comparison of the stage-1 fits of the four non-negative activity options with the COPS data and the bounded softplus c9
# (user 2026-11-14: "make comparison plot of the results with the data and the current softplus fit"; lead: per option all
# settings vs the data and c9; a summary figure with the best (lowest data cost) setting of each option + c9 + D2; a table
# per fit from each run's summary; activity maps of the best settings and c9).
#   options: B background (A = m_b [A_bg A_orig + softplus]), A sum of squares rank 3, D log form with spread penalty,
#            C local non-negative basis (CPU fits; COPS_model.dat from the design matrices)
# Writes $OUT (default figures_cluster/12M/model_vs_observations_nonneg_options_stage1; must not exist):
#   tables/*.dat (positive_table.py: observation rows of each fit's COPS_model.dat), runs.json, summary.txt / .csv / best.json
#   (options_table.py), 12M_options_<B|A|D|C>_* and 12M_options_best_* (model_obs_figures.py), captions (options_captions.py),
#   CAPTIONS.md (collect_captions.py), code/. The activity maps: jobs/options_activity_maps.sbatch (GPU, short) afterwards.
# Only runs whose summary exists are included (rerun into a new OUT when the remaining fits finish). CPU (srun / login).
# Options (2026-11-14, v2; unset = the first version's set):
#   DIAG=<json>     convergence verdicts + diag files (options_table.py --diag; e.g. scripts/data/stage1_convergence_20261114.json)
#   RESID=1         also 12M_options_best_residuals (options_residual_figure.py: ln(model/obs) vs altitude, NG and RG)
#   ACT_FROM=<dir>  copy <dir>/activity_maps (and its provenance) instead of a new GPU job, after checking that the best
#                   settings and their inputs are unchanged (same best.json and roots, inputs older than the maps)
# Options (2026-11-15, final; unset = as before):
#   RUNS=<json>     the candidate runs {"runs": [[option, name, setting, root relative to nonneg/, kind], ...]} instead of
#                   the built-in list (e.g. scripts/data/stage1_final_runs_20261115.json)
#   TEXTS=<json>    caption texts (options_captions.py --texts; e.g. scripts/data/stage1_final_texts_20261115.json)
#   CONV_TERM=1     convergence column from each fit's termination (options_table.py --conv-from-termination)
set -u
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; RUN=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; NN=$RUN/nonneg
F12=$COPS_OUTPUT_ROOT/figures_cluster/12M
O=${OUT:-$F12/model_vs_observations_nonneg_options_stage1}; TB=$O/tables
C9=$NN/softplus/chain_5cycles_smm1_20261113T041215Z/c9/ft/final; D2=$RUN/clip_refit/D2/final
[ -e "$O" ] && { echo "$O exists: use OUT=<new folder>"; exit 2; }
mkdir -p "$TB" "$O/code" || exit 1
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
PY=".venv/bin/python -W ignore"
cp -p "$S/figures_12m_nonneg_options_stage1.sh" "$S/options_table.py" "$S/options_captions.py" "$S/positive_table.py" "$S/model_obs_figures.py" "$S/collect_captions.py" "$O/code/" || exit 1   # 2026-11-14: absolute path (a relative $0 failed after the cd)
{ git -C /data001/highresolu_vt/AMPS rev-parse HEAD; git -C /data001/highresolu_vt/AMPS status --short; } > "$O/code/git.txt" 2>&1
# the runs (only those with a summary / report)
$PY - "$NN" "$C9" "$D2" > "$O/runs.json" <<'PYEOF' || exit 1
import json, os, sys
from pathlib import Path
nn, c9, d2 = Path(sys.argv[1]), sys.argv[2], sys.argv[3]
cand = [("B", "background", "bg 1e-3", "background/stage1_bg1e-3_20261114T084711Z", "summary"),
        ("B", "background", "bg 1e-2", "background/stage1_bg1e-2_20261114T084711Z", "summary"),
        ("B", "background", "bg 1e-1", "background/stage1_bg1e-1_20261114T084711Z", "summary"),
        ("A", "sum of squares", "rank 3", "sos/stage1_r3_20261114T091354Z", "summary"),
        ("D", "log spread", "f 0.1", "logspread/stage1_f0.1_20261114T090251Z", "summary"),
        ("D", "log spread", "f 1", "logspread/stage1_f1_20261114T090251Z", "summary"),
        ("D", "log spread", "f 10", "logspread/stage1_f10_20261114T090251Z", "summary"),
        ("C", "local basis", "f 0.1", "localnn/stage1_anch_f0.1_20261114T094256Z", "localnn"),
        ("C", "local basis", "f 1", "localnn/stage1_anch_f1_20261114T222136Z", "localnn"),   # rerun (the 094258Z attempt is .incomplete)
        ("C", "local basis", "f 10", "localnn/stage1_anch_f10_20261114T094259Z", "localnn")]
if os.environ.get("RUNS"):                              # 2026-11-15: candidate list from a file
    cand = [tuple(x) for x in json.load(open(os.environ["RUNS"]))["runs"]]
runs = []
for o, n, s, rel, k in cand:
    root = nn / rel
    need = root / ("summary.txt.json" if k == "summary" else "analysis/report.json")
    if need.exists() and (root / "fit/final/COPS_model.dat").exists():
        runs.append({"option": o, "name": n, "setting": s, "root": str(root), "kind": k, "tag": "%s_%s" % (o, s.replace(" ", "").replace(".", "p"))})
refs = [{"option": "ref", "name": "softplus c9", "setting": "9 cycles", "root": c9, "result": c9 + "/finetune_result.txt",
         "kind": "finetune", "report_key": "c9", "tag": "ref_softplus_c9",
         "multipliers": {"40_50": 33.51204, "50_60": 1.0, "60_70": 14.12697, "70_80": 1.0, "80_90": 17.03273, "90_180": 5.943917}},
        {"option": "ref", "name": "clipped D2", "setting": "D2", "root": d2, "result": d2 + "/fit_result.txt", "kind": "fit",
         "report_key": "D2", "tag": "ref_clip_D2"}]
json.dump({"runs": runs, "refs": refs}, sys.stdout, indent=1)
PYEOF
for t in $($PY -c "import json,sys; d=json.load(open('$O/runs.json')); print(' '.join(r['tag']+'|'+r['root'] for r in d['runs']))"); do
    $PY $S/positive_table.py "${t#*|}/fit/final" "$TB/${t%%|*}.dat" || exit 3
done
$PY $S/positive_table.py "$C9" "$TB/ref_softplus_c9.dat" || exit 3
$PY $S/positive_table.py "$D2" "$TB/ref_clip_D2.dat" || exit 3
DIAG=${DIAG:-}; RESID=${RESID:-0}; ACT_FROM=${ACT_FROM:-}; RUNS=${RUNS:-}; TEXTS=${TEXTS:-}; CONV_TERM=${CONV_TERM:-0}
for f in "$RUNS" "$TEXTS"; do [ -n "$f" ] && { cp -p "$f" "$O/code/" || exit 4; }; done
[ "$RESID" = "1" ] && { cp -p "$S/options_residual_figure.py" "$O/code/" || exit 4; }
if [ "$CONV_TERM" = "1" ]; then
    $PY $S/options_table.py --spec "$O/runs.json" --out-dir "$O" --conv-from-termination || exit 4
elif [ -n "$DIAG" ]; then
    cp -p "$DIAG" "$S/options_residual_figure.py" "$O/code/" || exit 4
    $PY $S/options_table.py --spec "$O/runs.json" --out-dir "$O" --diag "$DIAG" || exit 4
else
    $PY $S/options_table.py --spec "$O/runs.json" --out-dir "$O" || exit 4
fi
# figures: per option all settings + c9; the best of each option + c9 + D2
$PY - "$O" > "$O/figure_plan.sh" <<'PYEOF' || exit 5
import json, os, shlex, sys
from pathlib import Path
O = Path(sys.argv[1])
d = json.load(open(O / "runs.json"))
best = json.load(open(O / "best.json"))
rows = {}
import csv
for r in csv.DictReader(open(O / "summary.csv")):
    rows[(r["option"], r["setting"])] = r
styles = ["dashed", "dashdot", "dotted", "dashed", "dashdot", "dotted"]
lines = []
c9 = rows[("ref", "9 cycles")]
ref_fit = '--fit %s' % shlex.quote("softplus c9 (reference), %.4e=%s:solid" % (float(c9["data_cost"]), O / "tables/ref_softplus_c9.dat"))
names = {"A": "sum of squares", "B": "background", "C": "local non-negative basis", "D": "log form with spread penalty"}
for opt in sorted({r["option"] for r in d["runs"]}):
    fits = [r for r in d["runs"] if r["option"] == opt]
    args = [ref_fit]
    for i, r in enumerate(fits):
        dc = float(rows[(opt, r["setting"])]["data_cost"])
        args.append('--fit %s' % shlex.quote("%s %s, %.4e=%s:%s" % (opt, r["setting"], dc, O / "tables" / (r["tag"] + ".dat"), styles[i])))
    lines.append('$PY $S/model_obs_figures.py --bundle $RUN/bundle/restart_templates.npz --out %s --ref %s --title %s %s || exit 6' % (
        shlex.quote(str(O / ("12M_options_%s" % opt))), shlex.quote(str(O / "tables/ref_softplus_c9.dat")),
        shlex.quote("SHAP8 12M, option %s (%s), stage 1, vs the data and the bounded softplus c9" % (opt, names[opt])),
        " ".join(args)))
args = [ref_fit, '--fit %s' % shlex.quote("clipped D2 (start), %.4e=%s:dotted" % (float(rows[("ref", "D2")]["data_cost"]),
                                                                                 O / "tables/ref_clip_D2.dat"))]
for i, (opt, b) in enumerate(sorted(best.items())):
    tag = [r["tag"] for r in d["runs"] if r["option"] == opt and r["setting"] == b["setting"]][0]
    args.append('--fit %s' % shlex.quote("%s %s %s, %.4e=%s:%s" % (opt, b["model"], b["setting"], b["data_cost"],
                                                                   O / "tables" / (tag + ".dat"), ["dashed", "dashdot", "dashed", "dashdot"][i % 4])))
lines.append('$PY $S/model_obs_figures.py --bundle $RUN/bundle/restart_templates.npz --out %s --ref %s --title %s %s || exit 6' % (
    shlex.quote(str(O / "12M_options_best")), shlex.quote(str(O / "tables/ref_softplus_c9.dat")),
    shlex.quote("SHAP8 12M: best stage-1 fit of each non-negative option vs the data, the bounded softplus c9 and D2"),
    " ".join(args)))
if os.environ.get("RESID", "0") == "1":                 # v2: residuals of the same fits as the best figure
    rargs = []
    for x in args:
        lab = shlex.split(x)[1].split("=", 1)[0]
        key = {"softplus c9 (reference)": "ref|9 cycles", "clipped D2 (start)": "ref|D2"}.get(lab.split(",")[0])
        if key is None:
            o_ = lab.split()[0]
            key = "%s|%s" % (o_, best[o_]["setting"])
        rargs += [x, "--key %s" % shlex.quote("%s=%s" % (lab, key))]
    lines.append('$PY $S/options_residual_figure.py --bundle $RUN/bundle/restart_templates.npz --out %s --summary %s --title %s %s || exit 6' % (
        shlex.quote(str(O / "12M_options_best_residuals")), shlex.quote(str(O / "summary.csv")),
        shlex.quote("SHAP8 12M: log residuals of the best stage-1 fit of each option, the bounded softplus c9 and D2"),
        " ".join(rargs)))
print("\n".join(lines))
PYEOF
. "$O/figure_plan.sh"
if [ -n "$TEXTS" ]; then
    $PY $S/options_captions.py --out-dir "$O" --texts "$TEXTS" || exit 7
else
    $PY $S/options_captions.py --out-dir "$O" || exit 7
fi
if [ -n "$ACT_FROM" ]; then                            # v2: the activity maps of an earlier folder, inputs unchanged
    $PY - "$ACT_FROM" "$O" "$C9" <<'PYEOF' || exit 9
import hashlib, json, os, shutil, sys
from pathlib import Path
src, dst, c9 = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
b0, b1 = json.load(open(src / "best.json")), json.load(open(dst / "best.json"))
if {k: v["setting"] for k, v in b0.items()} != {k: v["setting"] for k, v in b1.items()}:
    sys.exit("best settings differ: %s vs %s" % (b0, b1))
r0 = {(r["option"], r["setting"]): r["root"] for r in json.load(open(src / "runs.json"))["runs"]}
r1 = {(r["option"], r["setting"]): r["root"] for r in json.load(open(dst / "runs.json"))["runs"]}
man = src / "activity_maps" / "activity_maps_manifest.json"
t_maps = man.stat().st_mtime
inputs = [Path(r1[(k, v["setting"])]) / "fit/final/nonneg_facets.npz" for k, v in b1.items()]
inputs += [c9 / "finetune_state.npz", c9 / "finetune_result.txt", c9.parent.parent / "fit/final/nonneg_facets.npz"]
lines = ["activity_maps/ copied from %s (figures_12m_nonneg_options_stage1.sh ACT_FROM; lead 2026-11-14: the inputs did not "
         "change). Checks:" % (src / "activity_maps"),
         "  best settings equal: %s" % {k: v["setting"] for k, v in b1.items()}]
for k, v in b1.items():
    if r0.get((k, v["setting"])) != r1[(k, v["setting"])]:
        sys.exit("run root of %s %s differs" % (k, v["setting"]))
lines.append("  run roots of the best settings equal")
for f in inputs:
    if not f.exists():
        sys.exit("missing input %s" % f)
    if f.stat().st_mtime >= t_maps:
        sys.exit("%s is newer than the maps (%s)" % (f, man))
    lines.append("  input older than the maps manifest: %s (%s)" % (f, __import__("time").strftime("%Y-%m-%dT%H:%M:%S", __import__("time").gmtime(f.stat().st_mtime))))
shutil.copytree(src / "activity_maps", dst / "activity_maps", copy_function=shutil.copy2)
for p in sorted((src / "activity_maps").iterdir()):
    h = [hashlib.sha256(open(q, "rb").read()).hexdigest() for q in (p, dst / "activity_maps" / p.name)]
    if h[0] != h[1]:
        sys.exit("copy differs: %s" % p.name)
    lines.append("  %s  %s (identical copy)" % (h[1], p.name))
for pv in sorted((src / "provenance").glob("activity_maps-*")):
    shutil.copytree(pv, dst / "provenance" / pv.name, copy_function=shutil.copy2)
    lines.append("  provenance copied: provenance/%s" % pv.name)
(dst / "activity_maps" / "COPIED_FROM.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
PYEOF
fi
$PY $S/collect_captions.py "$O" --title "Figure captions: stage-1 fits of the non-negative activity options vs the COPS data and the bounded softplus c9" --order "" || exit 8
(cd "$O/code" && sha256sum * > SHA256SUMS)
echo "done: $O"
