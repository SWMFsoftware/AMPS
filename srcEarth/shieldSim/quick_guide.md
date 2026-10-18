# shieldSim CRaTER parallel quick guide

This guide runs the clearly bounded `foreground_only_axial` CRaTER benchmark.
It is a diagnostic comparison with published observations, not a complete
CRaTER instrument model or Stage-8 validation.

Run every command from the `shieldSim` source directory. Keep generated files
under `test_output/`; do not redirect test or simulation files to `/tmp`.

## 1. Load Geant4 and locate the external reference package

```bash
source /data/vtenishe/G4/install/bin/geant4.sh
export Geant4_DIR=/data/vtenishe/G4/install/lib/cmake/Geant4
export SHIELDSIM_REFERENCE_DATA_DIR=/data/vtenishe/shieldSim_reference_data

test -f "${Geant4_DIR}/Geant4Config.cmake"
test -f "${SHIELDSIM_REFERENCE_DATA_DIR}/reference_crater/manifest.sha256"
```

The distributed example case points to this same external CRaTER directory.
If the package is relocated, change `source.path` in
`examples/stage8_crater_foreground_axial.json` and set the CMake path below to
the identical location.

## 2. Configure, compile, and run the Stage-8 gates

```bash
cmake -S . -B test_output/build \
  -DCMAKE_BUILD_TYPE=Release \
  -DBUILD_TESTING=ON \
  -DGeant4_DIR="${Geant4_DIR}" \
  -DSHIELDSIM_CRATER_REFERENCE_DIR="${SHIELDSIM_REFERENCE_DATA_DIR}/reference_crater"

cmake --build test_output/build --parallel 12 \
  --target shieldSimCraterBenchmark stage8_crater_core_tests

cd test_output/build
ctest -R '^stage8_(reference_inputs|crater_core|crater_transport)$' \
  --output-on-failure
cd ../..
```

## 3. Run with 12 local processes

```bash
python3 scripts/run_crater_parallel.py \
  --benchmark test_output/build/shieldSimCraterBenchmark \
  --case examples/stage8_crater_foreground_axial.json \
  --workers 12
```

The `histories` field in the case is the global count. The launcher divides it
by quotient/remainder, gives every shard a distinct deterministic seed, keeps
one Geant4 log and partial result per shard, and reports live shard status once
per minute. It refuses a failed or missing shard. The merger checks scientific
identity and sums raw importance-weight moments before recomputing rates and
Monte Carlo standard errors; normalized spectra are never averaged.

The final files are:

```text
test_output/crater_foreground_axial/
├── crater_bon2020_2009259_2010065__foreground_only_axial_result.json
├── crater_bon2020_2009259_2010065__foreground_only_axial_spectrum.csv
├── crater_bon2020_2009259_2010065__foreground_only_axial_comparison.png
└── crater_bon2020_2009259_2010065__foreground_only_axial_comparison.svg
```

Before using a result, require:

```text
result_type = complete_parallel_merge
parallel.merge_status = complete_all_shards_present
```

Also inspect `accuracy.status`; completing the requested histories does not
guarantee that every spectrum bin meets the requested 10% statistical target.

## 4. Optional `mpiexec` launch

The executable does not link MPI. Open MPI and PMI launchers provide rank and
size through environment variables that the executable converts into the same
shard definition:

```bash
mpiexec -n 12 test_output/build/shieldSimCraterBenchmark \
  examples/stage8_crater_foreground_axial.json

python3 scripts/merge_crater_shards.py \
  --output test_output/crater_foreground_axial/crater_bon2020_2009259_2010065__foreground_only_axial_result.json \
  test_output/crater_foreground_axial/crater_bon2020_2009259_2010065__foreground_only_axial__shard_*_partial_result.json

python3 scripts/plot_crater_foreground_benchmark.py \
  test_output/crater_foreground_axial/crater_bon2020_2009259_2010065__foreground_only_axial_result.json \
  test_output/crater_foreground_axial/crater_bon2020_2009259_2010065__foreground_only_axial_comparison
```

Do not use an older, non-shard-aware executable with `mpiexec`: each rank would
repeat the full run and could overwrite another rank's output. For detailed
build, schema, output, troubleshooting, and scientific-boundary information,
see [`RUN.md`](RUN.md) and
[`development_and_validatoin.md`](development_and_validatoin.md).

## Stage-9 verification-interface quick start

Stage 9 is available for verification workflows only; anomaly and preliminary-
mission submissions remain disabled until the prior scientific gates and
independent regression review are complete.

```bash
python3 scripts/stage9_interface.py preview \
  --request examples/stage9_ui_request.json

python3 scripts/stage9_ui_server.py \
  --host 127.0.0.1 --port 8765 \
  --output-root test_output/stage9_ui_jobs \
  --executable test_output/build/shieldSim
```

Open `http://127.0.0.1:8765/`, preview the request, compile a uniquely named
bundle, and run it. Use the UI `Cancel`/`Resume` controls or the equivalent CLI
commands documented in [`RUN.md`](RUN.md#10-stage-9-verification-interface).
Two shields × three absorbers means six transport cases regardless of how many
quantities are selected; quantity-series count is reported separately.

## Stage-10 analytic and mesh shape examples

The bounded non-slab providers support explicitly transformed analytic boxes,
spheres, cylinders, Boolean material regions, and closed outward-oriented
convex triangle meshes. They are verification-only and accept only
`purpose: "verification_test"`; science-facing slab workflows are unchanged.

```bash
test_output/build/shieldSim \
  --validate-case=examples/stage10_analytic_csg_case.json

for case_file in examples/stage10_shapes/*.json; do
  test_output/build/shieldSim --validate-case="$case_file" || exit 1
done

mkdir -p test_output/stage10_example
cd test_output/stage10_example
../build/shieldSim --case=../../examples/stage10_analytic_csg_case.json
cd ../build
ctest --output-on-failure -R '^stage10_(analytic_csg|tessellated_mesh)$'
cd ../..
```

The test writes numerical evidence, PNG/SVG comparison figures, and all run
artifacts under `test_output/build/`. The six gallery inputs are analytic box,
cylinder, sphere and meshed box, tetrahedron, octahedron. Exactly one sensitive
region is allowed; multiple non-overlapping passive material regions are
supported. Analytic entry surfaces use `local_negative_z_face` or
`entire_exterior`; meshes use explicit coplanar `tagged_facets` or
`entire_exterior`. `tessellated_mesh.v1` does not import or repair STL/CAD and
rejects concave/open/non-manifold meshes. See
[`RUN.md`](RUN.md#11-stage-10-analytic-and-tessellated-mesh-absorbers) for the
complete schemas, run commands, geometry/mass checks, hashing, limitations,
and test boundaries.

In every gallery JSON, select shielding in
`geometry.shield_stack.layers`, the absorber material in the sensitive
`geometry.absorbers[0].geometry.regions` entry, and calculations in
`quantities`. The committed cases request all five specified functions for
silicon. DDD/`n_eq` require `damage_response`; LET requires `let_response`.
Unsupported material/function combinations fail validation explicitly.
The separate `examples/stage10_analytic_csg_case.json` Boolean-enclosure
example follows the same input-selection contract.

For ready-to-edit combinations rather than one common default, use
[`examples/configuration_cases/`](examples/configuration_cases/):

```bash
for case_file in examples/configuration_cases/*.json; do
  test_output/build/shieldSim --validate-case="$case_file" || exit 1
done

mkdir -p test_output/configuration_cases/01_si_box_proton_all_quantities
cd test_output/configuration_cases/01_si_box_proton_all_quantities
../../build/shieldSim \
  --case=../../../examples/configuration_cases/01_si_box_proton_all_quantities.json
cd ../../..
```

The directory README maps every file to its incident species, absorber
material/shape, shield size and ordered layers, and quantities. Proton cases
can request all five names, but their neutron-only `n_eq` is defined as zero;
use either neutron damage example to exercise a nonzero `n_eq` calculation.
