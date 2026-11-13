#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_12m_nonneg5_smm1.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_12m_nonneg5_smm1.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 6a614e4adcf4e0114b1219985cefaa025a5316d2e51fbabcbb3612da232e5846
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: figure wrapper. The softplus refit with slope multipliers >= 1 (chain 124664) vs the unbounded softplus and
#   D2; the bounded chain's cycles; chain_stages.txt / summary.txt with termination, nfev and parameters at a bound;
#   HOW-TO-RERUN.md, code/ (2026-11-13).
# RUN: login node: figures_cluster/scripts/figures_12m_nonneg5_smm1.sh ->
#   figures_cluster/12M/model_vs_observations_softplus_smm1/
# ==== end of header ====
# The softplus refit with slope multipliers >= 1 (user 2026-11-12: refit the softplus model with the six free
# slope-bin multipliers 40-50 ... 90-180 bounded below by 1, start 1, no upper bound, 0-40 held at 1; everything else as
# chain_5cycles_20261112T012553Z: start D2, 5 fit / fine-tuning cycles, fine tuning late_only 5 min; nonneg_fit.py
# --slope-mult-min 1, frozen code model_ext_frozen/nonneg_20261113T035423Z; chain job 124664) against the COPS observations
# (user 2026-11-12 CST: "yes, make maps, contributions and comparison plots"). Copy of figures_12m_nonneg5.sh, adapted.
# Writes figures_cluster/12M/model_vs_observations_softplus_smm1/:
#   12M_softplus_smm1_final_*    the bounded softplus after 5 cycles vs the unbounded softplus after 5 cycles vs D2
#   12M_softplus_smm1_cycles_*   the bounded chain: the cycle-1 fit and the state after each of the 5 cycles
#   chain_stages.txt (both chains: every stage's cost, rms, termination / evaluations, active bounds),
#   summary.txt (final states; multipliers and which sit at their lower bound; A = 0 fractions), the cycle-5 bin tables
#   and NONNEG_MODEL.txt of both runs, MAPS.txt, HOW-TO-RERUN.md, code/.     Login node, no GPU:  figures_12m_nonneg5_smm1.sh
set -u
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; J=$COPS_OUTPUT_ROOT/jobs
RUN=$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5; F12=$COPS_OUTPUT_ROOT/figures_cluster/12M
NC=5
# the bounded refit (B) and the unbounded softplus run (U)
TAGB=chain_5cycles_smm1_20261113T041215Z; JOBB=${JOBB:-124664}; MAPJB=124667; CONJB=124668
TAGU=chain_5cycles_20261112T012553Z;      JOBU=124310;            MAPJU=124324; CONJU=124549
SB=$RUN/nonneg/softplus/$TAGB; SU=$RUN/nonneg/softplus/$TAGU
D2=$RUN/clip_refit/D2/final
NXB=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_20261113T035423Z; CXB=$COPS_OUTPUT_ROOT/model_ext_frozen/chain_onejob_20261113T041215Z
NXU=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_20261111T044920Z; CXU=$COPS_OUTPUT_ROOT/model_ext_frozen/chain_onejob_20261111T060742Z
MX=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_maps_20261112T021623Z
CTX=$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_contrib_20261112T210012Z
O=${OUT:-$F12/model_vs_observations_softplus_smm1}; TB=$O/tables; C=$O/code
[ -e "$O/12M_softplus_smm1_final_time.png" ] && { echo "$O has figures: move it aside first"; exit 2; }
for f in $SB/c$NC/ft/final/finetune_result.txt $SU/c$NC/ft/final/finetune_result.txt $D2/fit_result.txt; do
    [ -e "$f" ] || { echo "missing $f"; exit 2; }; done
mkdir -p "$TB" "$C"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
PY=".venv/bin/python -W ignore"
cost() { awk '/^cost_data/{printf "%.4e", $5}' "$1/finetune_result.txt"; }
fcost() { awk '/^final_data_cost/{printf "%.4e", $2}' "$1/fit_result.txt"; }
$PY $S/positive_table.py "$D2" "$TB/clip_D2.dat" || exit 1
$PY $S/positive_table.py "$SU/c$NC/ft/final" "$TB/softplus_unbounded_c${NC}ft.dat" || exit 1
for c in $(seq 1 $NC); do
    $PY $S/positive_table.py "$SB/c$c/fit/final" "$TB/softplus_smm1_c${c}fit.dat" || exit 1
    $PY $S/positive_table.py "$SB/c$c/ft/final" "$TB/softplus_smm1_c${c}ft.dat" || exit 1
done
$PY $S/model_obs_figures.py --bundle $RUN/bundle/restart_templates.npz --overwrite --out "$O/12M_softplus_smm1_final" \
    --ref "$TB/clip_D2.dat" \
    --title "SHAP8 12M (gauge +1.5 m, W10x5): softplus x slope multipliers after 5 cycles, multipliers ≥ 1 vs ≥ 0, and D2" \
    --fit "softplus, multipliers ≥ 1, 5 cycles, $(cost $SB/c$NC/ft/final)=$TB/softplus_smm1_c${NC}ft.dat:solid" \
    --fit "softplus, multipliers ≥ 0, 5 cycles, $(cost $SU/c$NC/ft/final)=$TB/softplus_unbounded_c${NC}ft.dat:dashed" \
    --fit "clipped refit D2, $(fcost $D2)=$TB/clip_D2.dat:dotted" || exit 3
STY=(dashdot dashed dashdot dashed solid)
CYC=(--fit "cycle-1 fit, $(fcost $SB/c1/fit/final)=$TB/softplus_smm1_c1fit.dat:dotted")
for c in $(seq 1 $NC); do
    CYC+=(--fit "after cycle $c, $(cost $SB/c$c/ft/final)=$TB/softplus_smm1_c${c}ft.dat:${STY[$((c - 1))]}")
done
$PY $S/model_obs_figures.py --bundle $RUN/bundle/restart_templates.npz --overwrite --out "$O/12M_softplus_smm1_cycles" \
    --ref "$TB/softplus_smm1_c${NC}ft.dat" \
    --title "SHAP8 12M (gauge +1.5 m, W10x5): softplus x slope multipliers ≥ 1, cycle-1 fit and the state after each of 5 cycles" \
    "${CYC[@]}" || exit 3
for m in smm1 unbounded; do CH=$([ $m = smm1 ] && echo $SB || echo $SU)
    cp -p "$CH/c$NC/fit/final/nonneg_bins.txt" "$O/softplus_${m}_c${NC}_fit_nonneg_bins.txt"
    cp -p "$CH/c$NC/fit/final/NONNEG_MODEL.txt" "$O/softplus_${m}_NONNEG_MODEL.txt"
done
# stage table and summary (cost, rms, termination / evaluations, active bounds, multipliers): python, read only
$PY - "$O" "$SB" "$SU" "$D2" "$NC" "$JOBB" "$JOBU" <<'PY' || exit 4
import json, re, sys
from pathlib import Path
import numpy as np
out, SB, SU, D2, NC, JOBB, JOBU = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]), int(sys.argv[5]), sys.argv[6], sys.argv[7]

def fit_info(d):
    t = (d / "fit_result.txt").read_text()
    g = lambda k: re.search(r"^%s\s+(\S+)" % k, t, re.M).group(1)
    opt = re.search(r"^optimizer: (.*)$", t, re.M).group(1)
    st = re.search(r"status (\d+), nfev (\d+)", opt)
    term = "ftol" if "ftol" in opt else "max_nfev" if "maximum number" in opt else "gtol" if "gtol" in opt else \
        "xtol" if "xtol" in opt else opt
    with np.load(d / "fit_restart.npz", allow_pickle=False) as z:
        names = json.loads(str(z["meta_json"]))["names"]
        x, lo, hi = (np.asarray(z[k], np.float64) for k in ("x", "lo", "hi"))
    tol = 1e-9 * np.maximum(1.0, np.abs(lo))
    atlo = [n for n, xi, l, tt in zip(names, x, lo, tol) if np.isfinite(l) and xi - l <= tt]
    athi = [n for n, xi, h in zip(names, x, hi) if np.isfinite(h) and h - xi <= 1e-9 * max(1.0, abs(h))]
    mult = {n[2:]: (xi, l) for n, xi, l in zip(names, x, lo) if n.startswith("c_slope_")}
    return {"cost": float(g("final_data_cost")), "rms_ng": float(g("rms_ln_residual_NG")),
            "rms_rg": float(g("rms_ln_residual_RG")), "term": term, "nfev": int(st.group(2)) if st else -1,
            "at_lower": atlo, "at_upper": athi, "mult": mult,
            "kappa": float(x[names.index("kappa")]), "grain_B": float(x[names.index("grain_B_uPa_per_100s")])}

def ft_info(d):
    t = (d / "finetune_result.txt").read_text()
    fin = lambda k: re.search(r"^%s\s+initial\s+\S+\s+final\s+(\S+)" % k, t, re.M).group(1)
    one = lambda k: re.search(r"^%s\s+(.*)$" % k, t, re.M).group(1).strip()
    return {"cost": float(fin("cost_data")), "rms_ng": float(fin("rms_ln_NG")), "rms_rg": float(fin("rms_ln_RG")),
            "term": one("termination"), "stop": one("stop"), "nfev": int(one("nfev_total")),
            "facets_at_lower": one("facets_at_lower_bound").split()[-1], "active": one("active_bounds")}

def zero_seen(d):
    for line in (d / "nonneg_bins.txt").read_text().splitlines():
        if line.startswith("all "):
            f = line.split()
            return {"zero_facets_frac": float(f[7]), "zero_area_frac": float(f[8]), "seen_zero_area_frac": float(f[9])}

L = ["# stage  data_cost  rms_ln_NG  rms_ln_RG  termination/evaluations  active bounds (fit: parameters at a bound, "
     "|x - bound| <= 1e-9 max(1, |bound|); fine tuning: its active_bounds line)"]
for lab, CH, job in (("softplus, multipliers >= 1 (bounded refit)", SB, JOBB), ("softplus, multipliers >= 0", SU, JOBU)):
    L.append("## %s: %s (job %s)" % (lab, CH, job))
    for c in range(1, NC + 1):
        f = fit_info(CH / ("c%d" % c) / "fit" / "final")
        L.append("c%dfit  %.6e %.4f %.4f %s/%d  at lower bound: %s%s" % (
            c, f["cost"], f["rms_ng"], f["rms_rg"], f["term"], f["nfev"], ", ".join(f["at_lower"]) or "none",
            ("; at upper bound: " + ", ".join(f["at_upper"])) if f["at_upper"] else ""))
        t = ft_info(CH / ("c%d" % c) / "ft" / "final")
        L.append("c%dft   %.10e %.6f %.6f %s/%d  facets at lower bound %s; %s" % (
            c, t["cost"], t["rms_ng"], t["rms_rg"], t["term"], t["nfev"], t["facets_at_lower"], t["active"]))
(out / "chain_stages.txt").write_text("\n".join(L) + "\n")

S = ["# final states (c%d/ft/final) and D2; multipliers from the last fit (c%d/fit/final; the fine tuning keeps them)" % (NC, NC),
     "# model  data_cost  rms_ln_NG  rms_ln_RG  termination  kappa"]
fin = {}
for lab, CH in (("softplus_multipliers_ge1_5cycles", SB), ("softplus_multipliers_ge0_5cycles", SU)):
    t, f = ft_info(CH / ("c%d" % NC) / "ft" / "final"), fit_info(CH / ("c%d" % NC) / "fit" / "final")
    fin[lab] = (t, f)
    S.append("%-36s %.10e %.6f %.6f %s %+.6e" % (lab, t["cost"], t["rms_ng"], t["rms_rg"], t["term"], f["kappa"]))
d2 = (D2 / "fit_result.txt").read_text()
g = lambda k: re.search(r"^%s\s+(\S+)" % k, d2, re.M).group(1)
S.append("%-36s %.6e %.4f %.4f %s %s" % ("clip_refit_D2", float(g("final_data_cost")), float(g("rms_ln_residual_NG")),
                                          float(g("rms_ln_residual_RG")), "max_nfev",
                                          re.search(r"^kappa\s+\S+\s+(\S+)", d2, re.M).group(1)))
cb, cu = fin["softplus_multipliers_ge1_5cycles"][0]["cost"], fin["softplus_multipliers_ge0_5cycles"][0]["cost"]
S.append("")
S.append("cost ratio bounded / unbounded = %.4f (difference %+.4e); bounded / D2 = %.4f" % (
    cb / cu, cb - cu, cb / float(g("final_data_cost"))))
S.append("")
S.append("# slope-bin multipliers (value, lower bound; '*' = at the lower bound within 1e-9 max(1, |bound|))")
names = list(fin["softplus_multipliers_ge1_5cycles"][1]["mult"])
S.append("%-14s %24s %24s" % ("bin", "multipliers >= 1", "multipliers >= 0"))
for n in names:
    row = []
    for lab in ("softplus_multipliers_ge1_5cycles", "softplus_multipliers_ge0_5cycles"):
        x, lo = fin[lab][1]["mult"][n]
        row.append("%.6e (>= %g)%s" % (x, lo, "*" if x - lo <= 1e-9 * max(1.0, abs(lo)) else " "))
    S.append("%-14s %24s %24s" % (n, row[0], row[1]))
S.append("slope_00_40    held at 1 in both")
S.append("")
S.append("# A = 0 (no H2O activity), from the cycle-%d fit's nonneg_bins.txt (the fine tuning changes only facets with A > 0)" % NC)
for lab, CH in (("multipliers >= 1", SB), ("multipliers >= 0", SU)):
    z = zero_seen(CH / ("c%d" % NC) / "fit" / "final")
    S.append("%-18s A = 0 on %.4f of the facets, %.4f of the area, %.4f of the seen area" % (
        lab, z["zero_facets_frac"], z["zero_area_frac"], z["seen_zero_area_frac"]))
(out / "summary.txt").write_text("\n".join(S) + "\n")
print("\n".join(L)); print("\n".join(S))
PY
{
echo "Activity and instantaneous source-rate maps of the final state (after cycle $NC), written by model_ext/nonneg_maps.py"
echo "(frozen code $MX; job script $J/nonneg_maps.sbatch). Each folder has README.txt (fields, units, law) and maps_manifest.json."
echo "  softplus, multipliers >= 1: $SB/maps   (job $MAPJB; geometry $RUN/maps_geometry/job_$MAPJB/geometry)"
echo "  softplus, multipliers >= 0: $SU/maps   (job $MAPJU)"
echo "Facet contributions, P% source areas and slope-bin fractions (model_ext/nonneg_contrib.py, frozen $CTX):"
echo "  multipliers >= 1: $SB/contributions (job $CONJB); multipliers >= 0: $SU/contributions (job $CONJU)"
echo "ParaView: <model>_flux.pvd (13 times) + <model>_static.vtu; Tecplot 360: <model>_maps.dat."
} > "$O/MAPS.txt"
for z in "$NXB" "$CXB" "$MX" "$CTX"; do mkdir -p "$C/model_ext_frozen/$(basename $z)"; cp -rp "$z/." "$C/model_ext_frozen/$(basename $z)/"; done
chmod -R u+w "$C/model_ext_frozen"; find "$C" -name __pycache__ -prune -exec rm -rf {} + 2>/dev/null
mkdir -p "$C/jobs" "$C/figures"
cp -p "$SB/provenance-$JOBB/chain_onejob.sbatch" "$C/jobs/"
cp -p $J/chain_onejob.sh $J/nonneg_maps.sbatch $J/nonneg_contrib.sbatch "$C/jobs/"
cp -p $S/figures_12m_nonneg5_smm1.sh $S/positive_table.py $S/model_obs_figures.py $S/rerun_commands.py "$C/figures/"
(cd "$C" && find . -type f ! -name SHA256SUMS -print0 | sort -z | xargs -0 sha256sum > SHA256SUMS)
CM=$($PY - "$C" <<'PY'
import hashlib, sys
from pathlib import Path
sys.path.insert(0, "/data001/highresolu_vt/cops_output/figures_cluster/scripts")
from rerun_commands import commit_of
for p in sorted(Path(sys.argv[1]).rglob("*")):
    if p.is_file() and p.suffix in (".py", ".sbatch", ".sh", ".cu", ".cpp"):
        h = hashlib.sha256(p.read_bytes()).hexdigest(); c = commit_of(p.name, h)
        print("| `%s` | `%s` | %s |" % (p.relative_to(sys.argv[1]), h[:16], "`%s`" % c if c else "**not in git**"))
PY
)
SUBL=/data001/highresolu_vt/cops_output/model_ext/logs/submissions.txt
{
cat <<MD
# How to rerun: 12M softplus x slope multipliers, refit with the multipliers bounded below by 1 (5 cycles)

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_12m_nonneg5_smm1.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never delete
earlier attempts (move them aside as \`*.incomplete_<UTC>\`). AMPS commit of the chain job $JOBB:
\`$(head -1 $SB/provenance-$JOBB/git.txt)\` (tracked changes: $(tail -n +2 $SB/provenance-$JOBB/git.txt | tr '\n' ' ')).

**Model and settings (user 2026-11-12):** the softplus model of \`chain_5cycles_20261112T012553Z\` (A = m_b s ln(1 +
exp(L/s)) exp(h)) with the six free slope-bin multipliers (40-50 ... 90-180 deg) bounded below by 1 (start 1, no upper
bound; 0-40 deg held at 1): \`nonneg_fit.py --slope-mult-min 1\` (frozen code \`$NXB\`, chain driver \`$CXB\`).
Everything else as the unbounded run: start D2, 5 fit / fine-tuning cycles, fine tuning late_only, 5-min window,
max_nfev 4000 (fit) / 500 (fine tuning). Maps and contributions of the final state: \`MAPS.txt\`.

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_12m_nonneg5_smm1.sh
\`\`\`

## A1. The chain as submitted (one job, 8 GPUs on gpu002; retries 124665/124666 only after a failure)

From \`$SUBL\`:

\`\`\`
$(grep -E "^[^ ]+ (124664|124665|124666) " $SUBL)
\`\`\`
i.e. \`cd $J && ./chain_onejob.sh --kind nonneg --model softplus --chain-root $SB --model-ext $NXB --chain-code $CXB --cycles 5 --nfev 4000 --check 1 --ft-maxnfev 500 --ft-facets late_only --ft-window-min 5 --l-threads 32 --slope-mult-min 1 --node gpu002 --time 08:00:00 --attempts 3 --submit\`
(the wrapper's \`--submit\` writes these sbatch lines). Inside the job (provenance-$JOBB/commands.sh):

\`\`\`bash
$(cat $SB/provenance-$JOBB/commands.sh)
\`\`\`

## A2. Maps and contributions of the final state (CPU jobs)

\`\`\`
$(grep -E "^[^ ]+ ($MAPJB|$CONJB) " $SUBL)
\`\`\`

### maps (job $MAPJB), commands as run

\`\`\`bash
$(cat $SB/maps/provenance-$MAPJB/commands.sh 2>/dev/null || echo "(see $SB/maps/provenance-$MAPJB)")
\`\`\`

### contributions (job $CONJB), commands as run

\`\`\`bash
$(cat $SB/contributions/provenance-$CONJB/commands.sh 2>/dev/null || echo "(see $SB/contributions/provenance-$CONJB)")
\`\`\`

## A3. Commands of every stage, as run
MD
for c in $(seq 1 $NC); do for st in fit ft; do p=$(ls -d $SB/c$c/$st/provenance-* 2>/dev/null | tail -1)
    [ -n "$p" ] && printf '\n### c%s/%s (%s)\n\n```bash\n%s\n```\n' "$c" "$st" "$(basename $p)" "$(cat $p/commands.sh)"
done; done
cat <<MD

## A4. The unbounded run and D2 (comparison curves)

\`$SU\` (job $JOBU): see \`../model_vs_observations_nonneg_5cycles/HOW-TO-RERUN.md\`; D2: \`$D2\`.

## Code that ran (\`code/\`, with \`SHA256SUMS\`) and the git commit with an identical copy

| file | SHA256 (first 16) | commit |
|---|---|---|
$CM
MD
} > "$O/HOW-TO-RERUN.md"
echo "wrote $O"
