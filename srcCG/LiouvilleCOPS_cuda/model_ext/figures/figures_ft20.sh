#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_ft20.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_ft20.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 88acd053664824cd1c54c40c61f065e1b605f2ada25224eac84952a9f2c7bf30
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison figures (model vs COPS observations; .eps + .png 300 dpi) + HOW-TO-RERUN.md for: 20M fit from the 12M fine-tuned state.
# RUN: on the login node, no GPU: <this> 
# OUT: $COPS_OUTPUT_ROOT/figures_cluster/20M/model_vs_observations_grain_from_12M_finetuned_state/
# ==== end of header ====
# Figures of the 20M grain fit started from the 12M fine-tuned state (model_ext/grain_ft20.sbatch; user request
# 2026-11-04 "run option 2 on 20M", integral-conserving correction template held at amplitude 1) and its 20M fine
# tuning, against the COPS observations and the translated 12M fit (model_obs_figures.py, as figures_grain.sh), and
# HOW-TO-RERUN.md. Login node:   figures_ft20.sh <12M run> <20M run> <job id of the fit>
set -u
R12=${1:?}; R20=${2:?}; JOB=${3:?}
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; X=$COPS_OUTPUT_ROOT/model_ext
N=C30_n2000_grain_cross10_T75_133K_Tmin_133K
FD=$R20/grain/ft12m/${N}_ft12m; FT=$R20/grain/ft12m/finetune/${N}_ft12m_m5_i5_cpu
XLD=$R20/grain/translated/${N}_from_12M; XFT=$R20/grain/translated/finetune/${N}_from_12M_m5_i5_cpu
SF=$R12/grain/fixed_T/$N; SFT=$R12/grain/fixed_T/finetune/${N}_m5_i5_cpu
O=$COPS_OUTPUT_ROOT/figures_cluster/20M/model_vs_observations_grain_from_12M_finetuned_state; TB=$O/tables
mkdir -p "$TB"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
.venv/bin/python $S/finetune_table.py "$FT" "$FD" "$TB/$(basename $FT).dat" || exit 1
.venv/bin/python $S/finetune_table.py "$SFT" "$SF" "$TB/12M_$(basename $SFT).dat" || exit 1
# the translated 12M fit on 20M (no correction): 'before' of its 20M fine tuning, geometry of the 20M base
.venv/bin/python - "$R20/base/COPS_model.dat" "$XFT/pressure_comparison.txt" "$TB/20M_translated_12M_fit.dat" <<'PY' || exit 1
import sys
from pathlib import Path
def table(path, header):
    names, rows = None, {}
    for line in Path(path).read_text().splitlines():
        if line.startswith(header):
            names = line[1:].split(); continue
        if line.startswith("#") or not line.strip() or names is None:
            continue
        v = line.split(); rows[int(v[0])] = dict(zip(names, v))
    return rows
base, pc = table(sys.argv[1], "#iPoint"), table(sys.argv[2], "# iPoint")
lines = ["# the translated 12M fit on 20M ('before' of %s); geometry from %s" % (sys.argv[2], sys.argv[1]),
         "#iPoint UTC et altitude pNG pRG refNG refRG useNG useRG"]
for i in sorted(pc):
    b, t = base[i], pc[i]
    lines.append(" ".join([str(i), b["UTC"], b["et"], b["altitude"], t["beforeNG"], t["beforeRG"], t["refNG"], t["refRG"],
                           t["useNG"], t["useRG"]]))
Path(sys.argv[3]).write_text("\n".join(lines) + "\n")
PY
P=".venv/bin/python -W ignore $S/model_obs_figures.py --bundle $R20/bundle/restart_templates.npz --overwrite"
L='cliff $\geq$60$^\circ$ template, T75 and Tmin fixed at 133 K, dust grain'
$P --out "$O/20M_fit_from_12M_finetuned_state" --ref "$FD" \
    --title "SHAP8 20M, $L: fit started from the 12M fine-tuned state (12M correction held fixed)" \
    --fit "12M fit on 20M (translated, no correction)=$TB/20M_translated_12M_fit.dat:dotted" \
    --fit "20M fit with the 12M correction held fixed=$FD:dashed" \
    --fit "20M fit with the 12M correction, fine-tuned on 20M=$TB/$(basename $FT).dat:solid" \
    --fit "12M fit, fine-tuned on 12M (reference)=$TB/12M_$(basename $SFT).dat:dashdot" \
    && echo "wrote 20M_fit_from_12M_finetuned_state_{time,altitude,altitude_NG,altitude_RG,altitude_NG_RG}"
TPL=$(dirname "$FD")/template_ft12m
cat > "$O/HOW-TO-RERUN.md" <<MD
# How to rerun: the 20M grain fit started from the 12M fine-tuned state, and its 20M fine tuning

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_ft20.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never delete
earlier attempts (move them aside as \`*.incomplete_<UTC>\`).

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_ft20.sh $R12 $R20 $JOB
\`\`\`

## A1. From the restart files: template, fit and fine tuning (one 8-GPU job, about 65 min)

The job (provenance of this result: \`$R20/grain/ft12m/provenance-$JOB/\`) runs \`model_ext/ft_template.py\` (the
12M fine-tuning correction carried to 20M with its integral conserved, template \`ft12m_late\`), then
\`model_ext/grain_fit_fixed.py\` (T75 = Tmin = 133 K and c_ft12m_late = 1 held; start from the translated restart
\`$XLD/fit_restart.npz\`; --template-start 0), then \`grain_finetune.py\`. It refuses an existing output: move
\`$FD\`, \`$TPL\` and \`$FT\` aside first.

\`\`\`bash
sbatch --export=ALL,COPS_RUN=$R20,SRC12=$SF,FT12=$SFT,XL=$XLD $X/grain_ft20.sbatch
\`\`\`

The translated restart itself: \`$COPS_OUTPUT_ROOT/figures_cluster/20M/model_vs_observations_grain_translated_finetune/HOW-TO-RERUN.md\` §A1.

## B. The 12M fit and its fine tuning, from the restart files or from scratch

\`$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_fixed_T/HOW-TO-RERUN.md\`, then the
translation (above), A1 and A0.
MD
echo "wrote $O/HOW-TO-RERUN.md"
