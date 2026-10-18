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
