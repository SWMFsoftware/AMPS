#!/bin/bash
# ==== version-controlled copy: header added 2026-11-18, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/illumination_numbers.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/illumination_numbers.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 58b9304712cf251e0bffb40da8f14128b1632cf502c1727ea3846739a5e80d9a
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: Slurm job (CPU): numbers and cross-case check of the illumination maps (illumination_maps_figures.py numbers;
#   2026-11-18)
# RUN: sbatch illumination_numbers.sbatch (defaults: stage-2 softplus reg state, the six full-set cases)
# ==== end of header ====
# Numbers and cross-case check of the illumination maps at the production-rate map times 10:33:10 and 10:39:26 UTC
# (figures_cluster/scripts/illumination_maps_figures.py numbers; 2026-11-18, user: "check the production rate plots --
# i dont see the effect of the shadow or night side decrease of the productions -- is it there? generate illumination
# maps ..."). CPU. The reference state's fields (pub_maps_figures.load_fields, bitwise against its maps' .vtu), the
# illumination factor q/q_max and its checks, areas / mean factor / night and shadow shares of the non-zero production
# area and of the production (whole and seen surface), class changes between the times; cos_SZA and illumination of the
# .vtu of every case against the reference, the geometry folders, the seen masks.
# Environment (defaults = the stage-2 softplus reg state and the six full-set cases): RUN, CHAIN, NCYC, FT, FIT, MAPS,
# GEOM, MX, OUT, REF_KEY, CASES ("KEY=MAPS_DIR KEY=MAPS_DIR ..."), MODEL_NAME (set it in the submitting shell: commas).
# Products: $OUT/checks/illumination_{numbers,crosscase}.{json,txt}; provenance in $OUT/provenance-<job>/.
# Shared node: heavy processes with oom_score_adj 1000 (CLAUDE.md rule 7).
#SBATCH --job-name=cops_illum_num
#SBATCH --partition=normal128cores
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=00:45:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_illum_num-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export COPS_DATA_ROOT=$T/rosetta_data_end_of_mission_cops
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
R=${RUN:-$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5}; BASE=$R/base
CH=${CHAIN:-$R/nonneg/softplus/stage2r_smm1_20261116T160802Z}; NCYC=${NCYC:-5}
FT=${FT:-$CH/c$NCYC/ft/final}; FIT=${FIT:-$CH/c$NCYC/fit/final}; MAPS=${MAPS:-$CH/maps_c$NCYC}
GEOM=${GEOM:-$R/maps_geometry/job_136152/geometry}; MX=${MX:-$COPS_OUTPUT_ROOT/model_ext_frozen/nonneg_maps_20261112T021623Z}
O=${OUT:-$COPS_OUTPUT_ROOT/figures_cluster/12M/illumination_maps_103310_103926}
REF_KEY=${REF_KEY:-stage2_softplus_reg}
N=$R/nonneg
CASES=${CASES:-stage2_softplus_reg=$MAPS stage2_A_reg=$N/sos/stage2r_r3_20261116T160802Z/maps_c5 stage2_B_reg=$N/background/stage2r_bg1e-2_20261116T160802Z/maps_c5 stage2_C_f1_pm1em4_reg=$N/localnn/stage2r_f1_pm1e-4_tplon_20261116T055923Z/c5/maps_ft stage2_D_reg=$N/logspread/stage2r_f0.1_20261116T160802Z/maps_c5 D2=$R/clip_refit/D2/maps}
MODEL_NAME=${MODEL_NAME:-the softplus activity model with slope multipliers >= 1, regularized stage-2 chain, final state after 5 fit and fine-tuning cycles}
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts
PY="$APP/.venv/bin/python -W ignore"
stamp() { echo "$(date -u +%FT%TZ) $*"; }
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
CARGS=(); for c in $CASES; do CARGS+=(--case "$c"); [ -e "${c#*=}/maps_manifest.json" ] || { stamp "missing: ${c#*=}/maps_manifest.json"; exit 2; }; done
for f in "$BASE/COPS_model.dat" "$BASE/elements.bin" "$FT/finetune_state.npz" "$FIT/fit_restart.npz" \
         "$FIT/nonneg_facets.npz" "$MAPS/maps_manifest.json" "$GEOM/times.txt" "$MX/nonneg_maps.py" \
         "$S/illumination_maps_figures.py" "$S/pub_maps_v2_figures.py" "$S/pub_maps_figures.py"; do
    [ -e "$f" ] || { stamp "missing: $f"; exit 2; }; done
mkdir -p "$O" || exit 1
P=$O/provenance-$SLURM_JOB_ID; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): illumination numbers, OUT=$O REF_KEY=$REF_KEY"
RUN=$R
for v in RUN BASE CH NCYC FT FIT MAPS GEOM MX O REF_KEY CASES MODEL_NAME; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/illumination_numbers.sbatch.used"
mkdir "$P/code" && cp -p "$S/illumination_maps_figures.py" "$S/pub_maps_v2_figures.py" "$S/pub_maps_figures.py" \
    "$APP/cops_cuda/figures.py" "$MX/nonneg_maps.py" "$P/code/" || exit 1
(cd "$P/code" && sha256sum * > SHA256SUMS)
# read-only git queries for provenance (--no-optional-locks: no index refresh is written)
{ git -C "$AMPS" rev-parse HEAD; git --no-optional-locks -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded)" >> "$P/module_list.txt"
echo "CPU job: no GPU (nvidia-smi not run)" > "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1 || $PY -m pip list > "$P/pip_list.txt" 2>&1
env | grep -E '^(SLURM_JOB_ID|SLURM_JOB_NODELIST|SLURM_CPUS_PER_TASK|SLURM_MEM_PER_NODE|OUT|RUN|CHAIN|NCYC|FT|FIT|MAPS|GEOM|MX|REF_KEY|CASES|OMP_NUM_THREADS|MPLBACKEND)=' | sort > "$P/env.txt"
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; "$@"; }
free -g > "$P/free_start.txt" 2>&1; df -h "$COPS_OUTPUT_ROOT" > "$P/df_start.txt" 2>&1
run heavy $PY "$S/illumination_maps_figures.py" numbers --out-dir "$O" --base "$BASE" --ft-dir "$FT" --fit-dir "$FIT" \
    --maps-dir "$MAPS" --geometry "$GEOM" --maps-code "$MX/nonneg_maps.py" --model-name "$MODEL_NAME" \
    --ref-key "$REF_KEY" "${CARGS[@]}" > "$P/illumination_numbers.log" 2>&1; rc=$?
tail -n 30 "$P/illumination_numbers.log"
free -g > "$P/free_end.txt" 2>&1; df -h "$COPS_OUTPUT_ROOT" > "$P/df_end.txt" 2>&1
stamp "done: rc $rc"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_illum_num-$SLURM_JOB_ID.out" "$P/job.log"
exit $rc
