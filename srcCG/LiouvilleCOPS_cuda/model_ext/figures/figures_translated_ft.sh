#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/figures_cluster/scripts/figures_translated_ft.sh (on Voyager: /data001/highresolu_vt/cops_output/figures_cluster/scripts/figures_translated_ft.sh)
# SHA256 of the run copy = SHA256 of this file without these header lines: 8e80f10e2f331e9b742d2c66075d62c0c0de03c611a56b67dadcadf6a8ac457d
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: comparison figures (model vs COPS observations; .eps + .png 300 dpi) + HOW-TO-RERUN.md for: a grain fit translated to another mesh and fine-tuned.
# RUN: on the login node, no GPU: <this> <source run> <source fit name> <target run> <target label>
# OUT: $COPS_OUTPUT_ROOT/figures_cluster/20M/model_vs_observations_grain_translated_finetune/
# ==== end of header ====
# Figures of a grain fit translated to another mesh and fine-tuned there (model_ext/grain_translate_finetune.sbatch;
# user request 2026-11-04: the 12M fit with T75 = Tmin = 133 K fixed -> fine tuning on the 20M surface model),
# against the COPS observations (model_obs_figures.py, as figures_grain.sh), and HOW-TO-RERUN.md. Login node:
#     figures_translated_ft.sh <runner dir of the source fit> <source fit name> <target runner dir> <target label>
set -u
RS=${1:?}; SRC=${2:?}; RT=${3:?}; RES=${4:?}
export COPS_OUTPUT_ROOT=${COPS_OUTPUT_ROOT:-/data001/highresolu_vt/cops_output}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts; X=$COPS_OUTPUT_ROOT/model_ext
SF=$RS/grain/fixed_T/$SRC; SFT=$RS/grain/fixed_T/finetune/${SRC}_m5_i5_cpu
N=${SRC}_from_12M; TD=$RT/grain/translated/$N; TFT=$RT/grain/translated/finetune/${N}_m5_i5_cpu
O=$COPS_OUTPUT_ROOT/figures_cluster/$RES/model_vs_observations_grain_translated_finetune; TB=$O/tables
mkdir -p "$TB"
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
.venv/bin/python $S/finetune_table.py "$SFT" "$SF" "$TB/12M_$(basename $SFT).dat" || exit 1
# target tables: 'before' (the translated state, with the grain) and 'after' of the fine tuning, geometry of the target base
.venv/bin/python - "$RT/base/COPS_model.dat" "$TFT/pressure_comparison.txt" "$TB" <<'PY' || exit 1
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
base, pc, out = table(sys.argv[1], "#iPoint"), table(sys.argv[2], "# iPoint"), Path(sys.argv[3])
for kind in ("before", "after"):
    lines = ["# %s fine tuning (%s); geometry from %s" % (kind, sys.argv[2], sys.argv[1]),
             "#iPoint UTC et altitude pNG pRG refNG refRG useNG useRG"]
    for i in sorted(pc):
        b, t = base[i], pc[i]
        lines.append(" ".join([str(i), b["UTC"], b["et"], b["altitude"], t[kind + "NG"], t[kind + "RG"], t["refNG"],
                               t["refRG"], t["useNG"], t["useRG"]]))
    (out / ("20M_translated_%s_finetuning.dat" % kind)).write_text("\n".join(lines) + "\n")
PY
P=".venv/bin/python -W ignore $S/model_obs_figures.py --bundle $RT/bundle/restart_templates.npz --overwrite"
L='cliff $\geq$60$^\circ$ template, T75 and Tmin fixed at 133 K, dust grain'
$P --out "$O/12M_fit_finetuned_on_20M" --ref "$SF" --title "$L: 12M fit, fine-tuned on 12M and on SHAP8 20M" \
    --fit "12M fit on 20M (translated)=$TB/20M_translated_before_finetuning.dat:dashed" \
    --fit "12M fit on 20M, fine-tuned on 20M=$TB/20M_translated_after_finetuning.dat:solid" \
    --fit "12M fit, fine-tuned on 12M=$TB/12M_$(basename $SFT).dat:dotted" \
    && echo "wrote 12M_fit_finetuned_on_20M_{time,altitude,altitude_NG,altitude_RG,altitude_NG_RG}"
J=$(ls -d $RT/grain/translated/provenance-* | tail -1)
cat > "$O/HOW-TO-RERUN.md" <<MD
# How to rerun: the 12M grain fit (T75 = Tmin = 133 K) translated to SHAP8 $RES and fine-tuned there

Written $(date -u '+%Y-%m-%d %H:%M UTC') by \`figures_cluster/scripts/figures_translated_ft.sh\`.
Rules (CLAUDE.md): GPU work only inside Slurm jobs; never write into an existing product directory; never
delete earlier attempts (move them aside as \`*.incomplete_<UTC>\`).

## A0. The figures only (login node, no GPU)

\`\`\`bash
$S/figures_translated_ft.sh $RS $SRC $RT $RES
\`\`\`

## A1. Translation and fine tuning from the 12M fit's restart file (one 8-GPU job, about 20 min)

The job (provenance of this result: \`$J/\`) runs \`model_ext/translate_grain.py\` (restart_tool.py translate, the
cliff template from the $RES bundle, $RES predictions with the source grain, the source fit's data selection)
and then \`model_ext/grain_finetune.py\` (5 min, 5 intervals, cpu, mixed, max_nfev 2000). It refuses an existing
output: move \`$TD\` and \`$TFT\` aside first.

\`\`\`bash
sbatch --export=ALL,SRC_FIT=$SF,COPS_RUN=$RT,NAME=$N $X/grain_translate_finetune.sbatch
\`\`\`

## B. The 12M fit itself, from its restart file or from scratch

See \`$COPS_OUTPUT_ROOT/figures_cluster/12M/model_vs_observations_grain_fixed_T/HOW-TO-RERUN.md\`, then A1, A0.
The $RES base, store and bundle: \`$RT/pipeline.json\` (runner) and
\`$COPS_OUTPUT_ROOT/figures_cluster/$RES/model_vs_observations/HOW-TO-RERUN.md\` §B1.
MD
echo "wrote $O/HOW-TO-RERUN.md"
