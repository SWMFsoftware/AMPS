#!/bin/bash
# ==== version-controlled copy: header added 2026-11-13, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/analysis/overhang_illumination_12m/run.sbatch (on Voyager: /data001/highresolu_vt/cops_output/analysis/overhang_illumination_12m/run.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 2cbb6e2ef328a9eb788cf899bbb3ab13bb01405fd964bf825b411f491d719e4b
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: job of analysis/overhang_illumination.py (CPU); job 124578 (out_124578/summary.txt).
# ==== end of header ====
#SBATCH --job-name=cops_overhang_ill
#SBATCH --partition=normal128cores
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=/data001/highresolu_vt/cops_output/analysis/overhang_illumination_12m/run-%j.out
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
D=/data001/highresolu_vt/cops_output/analysis/overhang_illumination_12m
R=/data001/highresolu_vt/cops_output/runs/shap8_012m_lategrid_dz1.5
cd /data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda || exit 1
{ git -C /data001/highresolu_vt/AMPS rev-parse HEAD; sha256sum $D/overhang_illumination.py $0; } > $D/provenance-$SLURM_JOB_ID.txt 2>&1
.venv/bin/python -W ignore $D/overhang_illumination.py $R/nonneg/softplus/chain_5cycles_20261112T012553Z $R/maps_geometry/job_124324/geometry $D/out_$SLURM_JOB_ID
