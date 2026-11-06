#!/bin/bash
# ==== version-controlled copy: header added 2026-11-06, not part of the code that ran ====
# Run copy: $COPS_OUTPUT_ROOT/jobs/hash_check.sbatch (on Voyager: /data001/highresolu_vt/cops_output/jobs/hash_check.sbatch)
# SHA256 of the run copy = SHA256 of this file without these header lines: 16aeeb62c4bea17069516d1c7dfac72e1ffa85ad606f1eb54ed3608c4eb5fc92
# Cluster paths, partition and account are as run (Voyager); change them to run elsewhere (README.md).
# WHAT: site check (handoff §9): SHA256 of the 12M base kernels and GPU-product arrays against the passion hashes.
#   CPU only, read only. Job 120867.
# RUN: sbatch hash_check.sbatch
# ==== end of header ====
#SBATCH --job-name=cops_hashes
#SBATCH --account=general
#SBATCH --partition=normal128cores
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=16 --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/data001/highresolu_vt/cops_output/jobs/cops_hashes-%j.out
# handoff §9: base kernel files and GPU-product arrays of the 12M earlier-set run against the passion hashes.
# CPU only, read only on the products; results in $COPS_OUTPUT_ROOT/checks/hash-<jobid>/.
set -uo pipefail
export AMPS=/data001/highresolu_vt/AMPS APP=/data001/highresolu_vt/AMPS/srcCG/LiouvilleCOPS_cuda
export T=/data001/highresolu_vt/cops_cluster_transfer_20261026 COPS_OUTPUT_ROOT=/data001/highresolu_vt/cops_output
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONUNBUFFERED=1
R=$COPS_OUTPUT_ROOT/runs/shap8_012m_flat_earlier
E=$T/expected_hashes/012m
O=$COPS_OUTPUT_ROOT/checks/hash-$SLURM_JOB_ID; mkdir -p "$O"
cp -p "$0" "$O/hash_check.sbatch"
{ hostname; date -u; git -C "$AMPS" rev-parse HEAD; git -C "$AMPS" status --short; } > "$O/provenance.txt" 2>&1
~/.local/bin/uv pip list --python "$APP/.venv/bin/python" > "$O/pip_list.txt" 2>&1

# 1. base elements.bin + 1,150 kernel files (paths in the list are relative to the base directory)
t0=$(date +%s)
( cd "$R/base" && awk '{print $2}' "$E/base_location_independent.sha256" \
    | xargs -P 16 -n 8 sha256sum ) | sort -k2 > "$O/base_location_independent.cluster.sha256"
sort -k2 "$E/base_location_independent.sha256" > "$O/base_location_independent.passion.sha256"
join -1 2 -2 2 "$O/base_location_independent.passion.sha256" "$O/base_location_independent.cluster.sha256" \
    | awk '$2!=$3{print $1}' > "$O/base_mismatches.txt"
echo "base: expected $(wc -l < "$E/base_location_independent.sha256") hashed $(wc -l < "$O/base_location_independent.cluster.sha256")" \
     "mismatches $(wc -l < "$O/base_mismatches.txt") ($(( $(date +%s) - t0 )) s)" | tee "$O/summary.txt"

# 2. arrays of the gravity field, cache, bundle and visibility templates (runner names -> package names)
"$APP/.venv/bin/python" -I "$T/scripts/012m_mono/npz_array_hashes.py" \
    "$R"/gravity/cuda_c04/equivalence_*/gpu_field.npz "$R/cache/admitted_gravity_cache.npz" \
    "$R/bundle/restart_templates.npz" "$R/visibility/visibility_templates.npz" \
  | sed -e 's/admitted_gravity_cache\.npz:/admitted_gravity_cache_shap8_012m_mixed.npz:/' \
        -e 's/restart_templates\.npz:/restart_templates_shap8_012m.npz:/' \
        -e 's/visibility_templates\.npz:/visibility_templates_shap8_012m.npz:/' \
  > "$O/gpu_products_arrays.cluster.sha256"
"$APP/.venv/bin/python" -I - "$E/gpu_products_arrays.sha256" "$O/gpu_products_arrays.cluster.sha256" <<'EOF' | tee -a "$O/summary.txt"
import sys
def load(p):
    d = {}
    for line in open(p):
        h, rest = line.rstrip("\n").split("  ", 1)
        d[rest] = h
    return d
ref, got = load(sys.argv[1]), load(sys.argv[2])
for key in sorted(set(ref) | set(got)):
    state = "missing" if key not in got else "extra" if key not in ref else "identical" if ref[key] == got[key] else "DIFFERS"
    print("%-10s %s" % (state, key))
EOF
date -u >> "$O/provenance.txt"
