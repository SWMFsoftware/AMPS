#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/cops_movies_qc.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/cops_movies_qc.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: dcddd5aadf37c9dc9c38efc3583d53982d7c5e0c5372793a787a457b85ae0451
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job script (CPU): PSNR / SSIM of the MP4s against the composed frames and decoded stills.
# RUN: sbatch --export=ALL cops_movies_qc.sbatch
# ==== end of header ====
# COPS movies, step 4 (CPU): quality check of the encoded MP4s against the composed PNG frames (ffmpeg psnr / ssim filters;
# the 1280x720 version against the frames scaled the same way, lanczos) and decoded stills of a few frames for a visual
# check of the compression. Environment: WORK, COMPOSED (default $WORK/composed), OUT (MP4 folder), STILLS (frame list).
#SBATCH --job-name=cops_movies_qc
#SBATCH --partition=normal128cores
#SBATCH --account=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/logs/cops_movies_qc-%j.out
set -uo pipefail
export COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
R=${RUN:-$COPS_OUTPUT_ROOT/runs/shap8_012m_lategrid_dz1.5}
CH=${CHAIN:-$R/nonneg/softplus/chain_5cycles_20261112T012553Z}
WORK=${WORK:-$CH/movies_work}; COMPOSED=${COMPOSED:-$WORK/composed}
OUT=${OUT:-$COPS_OUTPUT_ROOT/figures_cluster/12M/movies_softplus_5cycles}
STILLS=${STILLS:-0 688 1200 1799}
FFMPEG=$COPS_OUTPUT_ROOT/tools/imageio_ffmpeg/imageio_ffmpeg/binaries/ffmpeg-linux-x86_64-v9.0.2
stamp() { echo "$(date -u +%FT%TZ) $*"; }
P=$WORK/provenance-$SLURM_JOB_ID; mkdir "$P" || exit 1
cp -p "$0" "$P/cops_movies_qc.sbatch.used"
for v in WORK COMPOSED OUT STILLS FFMPEG; do printf '%s=%s\n' "$v" "${!v}"; done > "$P/inputs.txt"
{ "$FFMPEG" -hide_banner -version | head -1; sha256sum "$FFMPEG"; lscpu | grep 'Model name'; } > "$P/toolchain.txt" 2>&1
stamp "job $SLURM_JOB_ID on $(hostname)"
for i in 1 2 3; do
    case $i in 1) n=movie1_RG_view;; 2) n=movie2_NG_view;; 3) n=movie3_combined;; esac
    for tag in 1920x1080 1280x720_web; do
        mp4=$OUT/${n}_${tag}.mp4
        if [ $tag = 1920x1080 ]; then ref="[1:v]format=yuv420p[r]"; else ref="[1:v]scale=1280:720:flags=lanczos,format=yuv420p[r]"; fi
        res=$("$FFMPEG" -hide_banner -nostats -i "$mp4" -framerate 30 -i "$COMPOSED/movie$i/f%05d.png" \
              -lavfi "[0:v]format=yuv420p[m];$ref;[m]split[m1][m2];[r]split[r1][r2];[m1][r1]psnr;[m2][r2]ssim" -f null - 2>&1 \
              | grep -E "PSNR|SSIM" | sed 's/.*\] //')
        echo "$n $tag: $res" | tee -a "$P/quality.txt"
        for k in $STILLS; do
            "$FFMPEG" -hide_banner -loglevel error -i "$mp4" -vf "select=eq(n\,$k)" -vframes 1 "$P/decoded_${n}_${tag}_f$(printf %05d $k).png"
        done
    done
done
stamp "done"
cp -p "$COPS_OUTPUT_ROOT/jobs/logs/cops_movies_qc-$SLURM_JOB_ID.out" "$P/job.log"
