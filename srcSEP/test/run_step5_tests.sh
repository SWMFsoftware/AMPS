#!/bin/sh
set -eu

# Step 5 is a source-scope change.  The checks deliberately inspect the files
# that enter the srcSEP production library instead of searching documentation,
# migration notes, or this test itself, where removed names must remain
# recorded.  That boundary makes SCOPE01 a useful regression guard rather than
# a brittle repository-wide word filter.
src_root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)

production_files="
$src_root/sep.h
$src_root/sep.dfn
$src_root/sep.cpp
$src_root/mover_state.cpp
$src_root/focused_transport_dmumu.cpp
$src_root/focused_transport_mfp.cpp
$src_root/parker_mover.cpp
$src_root/transport_common.cpp
$src_root/main.cpp
$src_root/main_lib.cpp
$src_root/mesh.cpp
$src_root/inner_boundary_injection.cpp
$src_root/sampling.cpp
$src_root/sampling_output.cpp
$src_root/field_line.cpp
$src_root/solar_wind.cpp
$src_root/QLT1.h
$src_root/QLT1.cpp
$src_root/makefile
"

forbidden_symbols='ParticleMover_Parker3D_MeanFreePath|ParticleMover_HE_2019_AJL|ParticleMover__He_2019_AJL|ParticleMover_Kartavykh_2016_AJ|ParticleMover_BOROVIKOV_2019_ARXIV|ParticleMover_default|Relativistic::Boris|_SEP_MOVER_BOROVIKOV_2019_ARXIV_|_SEP_MOVER_HE_2019_AJL_|_SEP_MOVER_KARTAVYKH_2016_AJ_|_SEP_MOVER_DRIFT_|ParticleTrajectoryCalculation|GetDriftVelocity|InitDriftVelData|b_times_grad_absB_offset|CurlB_offset|b_b_Curl_B_offset|tempParticleMovingListTable|FirstCellParticleTable|Sample3D|sample3d'

if grep -En "$forbidden_symbols" $production_files >/dev/null; then
  echo "FAIL SCOPE01: transferred full-3D symbols remain in production sources" >&2
  grep -En "$forbidden_symbols" $production_files >&2 || true
  exit 1
fi

for removed_file in drift.cpp output.cpp sample3d.h sample3d_init.cpp sample3d_sampling.cpp sample3d_output.cpp; do
  if test -e "$src_root/$removed_file"; then
    echo "FAIL SCOPE01: transferred source still exists: $removed_file" >&2
    exit 1
  fi
done
echo "PASS SCOPE01: transferred mover and sample3d symbols are absent"

# Production movers may convert field-line coordinates to Cartesian vectors in
# order to evaluate a three-dimensional background.  They must not advance a
# particle by writing Cartesian position, invoke a Boris pusher, diffuse it
# across field lines, or attach it to an AMR mesh-cell particle list.
advance_files="$src_root/sep.h $src_root/mover_state.cpp $src_root/focused_transport_dmumu.cpp $src_root/focused_transport_mfp.cpp $src_root/parker_mover.cpp $src_root/transport_common.cpp $src_root/production_mover_runtime.cpp"
cartesian_advance='ParticleBuffer::SetX|PB::SetX|Relativistic::Boris|calculatePerpendicularDiffusion|_PIC_PARTICLE_LIST_ATTACHING_NODE_'
if grep -En "$cartesian_advance" $advance_files >/dev/null; then
  echo "FAIL SCOPE02: Cartesian particle-advance mechanism remains" >&2
  grep -En "$cartesian_advance" $advance_files >&2 || true
  exit 1
fi
grep -q '_PIC_PARTICLE_LIST_ATTACHING_FL_SEGMENT_' "$src_root/production_mover_runtime.cpp" || {
  echo "FAIL SCOPE02: production adapter does not enforce field-line attachment" >&2
  exit 1
}
echo "PASS SCOPE02: production advance is field-line-only"

# Three-dimensional *geometry* is intentionally retained.  These checks guard
# the public flux-tube abstraction, Cartesian field-line embedding, imported
# vector magnetic field, and field-line observer sampling that replaced the
# removed spatial-volume sample3d module.
grep -q 'namespace FluxTubeGeometry' "$src_root/sep.h" || {
  echo "FAIL SCOPE03: flux-tube geometry interface is missing" >&2
  exit 1
}
grep -q 'GetCartesian' "$src_root/field_line.cpp" || {
  echo "FAIL SCOPE03: three-dimensional field-line embedding is missing" >&2
  exit 1
}
grep -q 'GetMagneticField' "$src_root/flux_tube_geometry.cpp" || {
  echo "FAIL SCOPE03: vector magnetic-field access is missing" >&2
  exit 1
}
grep -q 'InitSingleFieldLineSampling' "$src_root/sampling.cpp" || {
  echo "FAIL SCOPE03: field-line observer sampling is missing" >&2
  exit 1
}
echo "PASS SCOPE03: 3-D field-line geometry and field-line sampling are retained"

# Re-run the exact production-registry check: removing implementations must not
# broaden or damage the three canonical public choices introduced in Step 4.
"$src_root/test/run_step4_tests.sh"
