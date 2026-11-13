#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/cops_movies_encode.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/cops_movies_encode.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: d2a2152aa2d0eaddbe0d5adc5db4c04b88b56f37be22248884941b61965c4696
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (CPU, normal128cores): compare, frametable, compose (overlays, colour bars, the movie-3 time
#   series; TS_LABEL; 2026-11-13 TOUCHDOWN = touchdown.npz: the waving touchdown flag), encode (H.264 1920x1080 CRF 18,
#   1280x720 CRF 26 faststart); MOVIES_PY selects the cops_movies.py to run.
# RUN: STEPS=frametable,compose,encode OUT=... TABLE=... [TS_LABEL=... TOUCHDOWN=...] sbatch --export=ALL
#   cops_movies_encode.sbatch (OUTPUT_VISUALIZATION.md 3.7). Jobs: 124631-124634 (first movies), 124732, 124734 (c5),
#   124866 (c9 gauge view), 124868 (c9 fixed camera)
# ==== end of header ====
# COPS movies, step 3 (CPU): compare (split check of two rendered frame folders), compose (movie frames 1, 2, 3 with
# overlays and the NG / RG time series), encode (H.264 MP4, 1920x1080 CRF 18 and 1280x720 CRF 26, faststart).
# figures_cluster/scripts/cops_movies.py; ffmpeg = imageio-ffmpeg's static binary in $COPS_OUTPUT_ROOT/tools/imageio_ffmpeg.
# Environment (values with commas: set them in the submitting shell and use --export=ALL):
#   WORK (default $CHAIN/movies_work), STEPS (default compose,encode; also compare), COMPOSED (composed frames, default
#   $WORK/composed), FRAMES (compose only these frames, comma list; for sample stills), OUT (MP4 folder, default
#   figures_cluster/12M/movies_softplus_5cycles), TABLE (time-series table), CMP_A, CMP_B (frame folders to compare),
#   PROCS (default 32); STEPS may also hold frametable (2026-11-13: frames.tsv from frames/meta -> FRAMETABLE, default
#   $OUT/validation/frames.tsv); TS_LABEL (text right of the movie-3 time series, '\n' = new line; default: the first
#   movies' "COPS data and the softplus activity model (final state, 5 cycles)"); TOUCHDOWN (2026-11-13: touchdown.npz of
#   cops_movies.py touchdown -> compose --touchdown: the touchdown flag; default none); MOVIES_PY (the cops_movies.py to
#   run, default figures_cluster/scripts/cops_movies.py).
#SBATCH --job-name=cops_movies_encode
#SBATCH --partition=normal128cores
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_movies_encode-%j.out
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1 MPLBACKEND=Agg
R=${RUN:-$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5}
CH=${CHAIN:-$R/nonneg/softplus/chain_5cycles_20261112T012553Z}
WORK=${WORK:-$CH/movies_work}; STEPS=${STEPS:-compose,encode}
COMPOSED=${COMPOSED:-$WORK/composed}; FRAMES=${FRAMES:-}; PROCS=${PROCS:-32}
OUT=${OUT:-$COPS_OUTPUT_ROOT/figures_cluster/12M/movies_softplus_5cycles}
TABLE=${TABLE:-$COPS_OUTPUT_ROOT/figures_cluster/12M/full_set_softplus_5cycles/model_data/tables/softplus_5cycles_final.dat}
CMP_A=${CMP_A:-$WORK/frames}; CMP_B=${CMP_B:-}; FRAMETABLE=${FRAMETABLE:-$OUT/validation/frames.tsv}; TS_LABEL=${TS_LABEL:-}
FFMPEG=$COPS_OUTPUT_ROOT/tools/imageio_ffmpeg/imageio_ffmpeg/binaries/ffmpeg-linux-x86_64-v9.0.2
S=$COPS_OUTPUT_ROOT/figures_cluster/scripts
MPY=${MOVIES_PY:-$S/cops_movies.py}; TOUCHDOWN=${TOUCHDOWN:-}
PY="$APP/.venv/bin/python -W ignore"
stamp() { echo "$(date -u +%FT%TZ) $*"; }
heavy() { bash -c 'echo 1000 > /proc/self/oom_score_adj; exec "$@"' heavy "$@"; }
P=$WORK/provenance-$SLURM_JOB_ID; mkdir "$P" || exit 1
cd "$APP" || exit 1
stamp "job $SLURM_JOB_ID on $(hostname): STEPS=$STEPS"
for v in WORK STEPS COMPOSED FRAMES PROCS OUT TABLE CMP_A CMP_B FFMPEG TS_LABEL TOUCHDOWN MPY; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
cp -p "$0" "$P/cops_movies_encode.sbatch.used"
mkdir "$P/code" && cp -p "$MPY" "$APP/cops_cuda/figures.py" "$P/code/" && (cd "$P/code" && sha256sum * > SHA256SUMS)
{ git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$P/git.txt" 2>&1
{ type module >/dev/null 2>&1 && module list; } > "$P/module_list.txt" 2>&1 || echo "module command not available (no modules loaded)" >> "$P/module_list.txt"
echo "CPU job: no GPU (nvidia-smi not applicable)" > "$P/nvidia_smi_q.txt"
{ gcc --version | head -1; "$FFMPEG" -hide_banner -version | head -3; sha256sum "$FFMPEG"; lscpu | grep -E 'Model name|^CPU\(s\)'; } > "$P/toolchain.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$P/pip_list.txt" 2>&1
{ echo '#!/bin/bash'; grep '^export ' "$0"; echo "cd $APP"; } > "$P/commands.sh"
run() { printf '%q ' "$@" >> "$P/commands.sh"; echo >> "$P/commands.sh"; "$@"; }
free -g > "$P/free_start.txt" 2>&1
nfail=0
if [[ ",$STEPS," == *,compare,* ]]; then
    stamp "compare $CMP_A $CMP_B"
    run heavy $PY "$MPY" compare --work "$WORK" --a "$CMP_A" --b "$CMP_B" --procs "$PROCS" \
        --report "$P/compare_report.json" > "$P/compare.log" 2>&1 || { stamp "compare: frames differ or FAILED"; nfail=$((nfail+1)); }
    tail -n 3 "$P/compare.log"
fi
if [[ ",$STEPS," == *,frametable,* ]]; then
    stamp "frametable -> $FRAMETABLE"
    mkdir -p "$(dirname "$FRAMETABLE")"
    run $PY "$MPY" frametable --work "$WORK" --out "$FRAMETABLE" > "$P/frametable.log" 2>&1 || { stamp "frametable FAILED"; nfail=$((nfail+1)); }
fi
if [[ ",$STEPS," == *,compose,* ]]; then
    stamp "compose -> $COMPOSED"
    mkdir -p "$COMPOSED" || exit 1
    run heavy $PY "$MPY" compose --work "$WORK" --out "$COMPOSED" --table "$TABLE" --procs "$PROCS" \
        ${FRAMES:+--frames "$FRAMES"} ${TS_LABEL:+--ts-label "$TS_LABEL"} ${TOUCHDOWN:+--touchdown "$TOUCHDOWN"} > "$P/compose.log" 2>&1 || { stamp "compose FAILED"; nfail=$((nfail+1)); }
    tail -n 3 "$P/compose.log"
fi
if [[ ",$STEPS," == *,encode,* ]] && [ $nfail -eq 0 ]; then
    stamp "encode -> $OUT"
    mkdir -p "$OUT" || exit 1
    run heavy $PY "$MPY" encode --work "$WORK" --out "$OUT" --frames-dir "$COMPOSED" --ffmpeg "$FFMPEG" \
        > "$P/encode.log" 2>&1 || { stamp "encode FAILED"; nfail=$((nfail+1)); }
    tail -n 8 "$P/encode.log"
fi
free -g > "$P/free_end.txt" 2>&1
stamp "done: $nfail failures"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_movies_encode-$SLURM_JOB_ID.out" "$P/job.log"
[ $nfail -eq 0 ]
