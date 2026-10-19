# shieldSim parallel quick guide

The first part of this guide runs the clearly bounded `foreground_only_axial`
CRaTER benchmark. It is a diagnostic comparison with published observations,
not a complete CRaTER instrument model or Stage-8 validation. A later section
records the exact eight-process procedure used for the high-statistics
BON2020-proton configuration example.

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

## Stage-11 calculated AMPS-spectrum input

The runnable case selects the spectrum filename and component in JSON:

```json
"spectrum": {
  "format": "spectrum_exchange_v1",
  "path": "../tests/data/stage11_spectrum_exchange.v1.json",
  "component_id": "local_proton_points"
}
```

Build, validate, run, and execute the two focused tests:

```bash
cmake -S . -B test_output/build -DBUILD_TESTING=ON \
  -DGeant4_DIR=/path/to/Geant4/lib/cmake/Geant4
cmake --build test_output/build -j 8 \
  --target stage11_spectrum_exchange_tests shieldSim

test_output/build/shieldSim \
  --validate-case=examples/stage11_spectrum_exchange_case.json
test_output/build/shieldSim \
  --case=examples/stage11_spectrum_exchange_case.json

cd test_output/build
ctest -R '^stage11_spectrum_exchange(_transport)?$' --output-on-failure
cd ../..
```

The format is strict JSON, uses local `J_local_perMeV` rather than an already
angular-integrated AMPS column, and currently transports one explicit proton
or alpha component per case. See
[`SPECTUM_INPUT_FORMAT.md`](SPECTUM_INPUT_FORMAT.md) for the complete producer
schema, units, interpolation, quality/provenance rules, and failure boundaries.

## High-statistics BON2020 proton example: exact eight-process procedure

This section records the command sequence used for the completed
`01_si_box_bon2020_proton_all_quantities_directional` test. The input is the
prepared 400-node proton component of the checksum-pinned CRaTER BON2020 GCR
spectrum. The simulation used 100,000 independent histories, divided into
eight workers of 12,500 histories each. The completed run is under:

```text
test_output/configuration_cases/
└── 01_si_box_bon2020_proton_all_quantities_directional/
    └── parallel_100000/
        ├── shard_0000/ ... shard_0007/
        ├── parallel_run_manifest.json
        ├── parallel_completion.json
        └── pooled/
```

The generic `shieldSim` executable does **not** implement the rank/shard
contract provided by `shieldSimCraterBenchmark`. Consequently, this case was
not run as `mpiexec -n 8 shieldSim ...`: doing that would make every rank repeat
the entire input and would allow ranks to overwrite the same files. Instead,
the procedure below constructs eight independent input files, launches eight
ordinary `shieldSim` processes concurrently, and pools only review plots and
summary statistics after all workers succeed. Each shard result JSON remains a
canonical result; the pooled files are explicitly visualization products and
are not represented as one canonical merged transport result.

### 1. Load the runtime and name the exact paths

Run this block from the `shieldSim` source directory:

```bash
source /data/vtenishe/G4/install/bin/geant4.sh

export SHIELDSIM_EXE="$PWD/test_output/build/shieldSim"
export SHIELDSIM_CASE_ROOT="$PWD/test_output/configuration_cases/01_si_box_bon2020_proton_all_quantities_directional"
export SHIELDSIM_BASE_CASE="${SHIELDSIM_CASE_ROOT}/01_si_box_bon2020_proton_all_quantities_directional.json"
export SHIELDSIM_RUN_ROOT="${SHIELDSIM_CASE_ROOT}/parallel_100000"

test -x "${SHIELDSIM_EXE}"
test -f "${SHIELDSIM_BASE_CASE}"
```

`source` supplies the Geant4 shared-library and physics-data paths required at
runtime. The four exported names prevent path changes between preparation,
execution, checking, and plotting. The two `test` commands fail immediately if
the executable or input case is missing.

The executable had already been configured and compiled in
`test_output/build`, so no CMake command was part of this particular execution.
For a clean machine, perform the configure/build commands in Section 2 first
and include the `shieldSim` target. Likewise, the BON2020 source adapter had
already produced `SHIELDSIM_BASE_CASE`; its spectrum and adapter provenance are
recorded beside that case in
`bon2020_crater_solmin_proton_10MeV_50GeV.metadata.v1.json` and
`source_adapter_manifest.json`. Those preparation products were inputs to this
run, not files downloaded or regenerated while the eight workers were active.

`SHIELDSIM_RUN_ROOT` above is the location of the completed run. Do not execute
the preparation block below over that nonempty directory. For a new repeat,
set it to a new name such as `parallel_100000_repeat_01`; never mix new shards
with an earlier run.

### 2. Create the eight independent shard inputs

The following is the preparation command used for this run. It changes only
the case ID, history controls, random seed, output prefix, and paths that must
remain valid from the deeper shard directories. It does not alter the source,
geometry, material definitions, requested quantities, physics list, cuts, or
step limit.

```bash
PYTHONDONTWRITEBYTECODE=1 python3 - <<'PY'
import copy
import json
import os
from pathlib import Path

source_path = Path(os.environ["SHIELDSIM_BASE_CASE"]).resolve()
run_root = Path(os.environ["SHIELDSIM_RUN_ROOT"]).resolve()
repository = Path.cwd().resolve()
worker_count = 8
total_histories = 100_000
histories_per_worker = total_histories // worker_count
base_seed = 1_020_200

if total_histories % worker_count:
    raise SystemExit("this recorded run requires an equal integral shard size")
if run_root.exists() and any(run_root.iterdir()):
    raise SystemExit(f"refusing to mix with nonempty run directory: {run_root}")

base_case = json.loads(source_path.read_text(encoding="utf-8"))
run_root.mkdir(parents=True, exist_ok=True)

for index in range(worker_count):
    shard = copy.deepcopy(base_case)
    tag = f"shard_{index:04d}_of_{worker_count:04d}"
    shard["case_id"] = f"configuration_example_01_si_box_bon2020_proton_parallel_{tag}"
    shard["accuracy"]["sampling"]["min_histories"] = histories_per_worker
    shard["accuracy"]["sampling"]["max_histories"] = histories_per_worker
    shard["accuracy"]["sampling"]["batch_size"] = histories_per_worker
    shard["transport"]["random_seed"] = base_seed + index
    shard["output"]["prefix"] = f"bon2020_proton_parallel_{tag}"

    # The base case uses paths relative to its original directory. Absolute
    # paths preserve the identical libraries when the case moves into a shard.
    shard["material_library"]["path"] = str(
        repository / "data/material_library.v1.json"
    )
    shard["damage_response"]["path"] = str(
        repository / "data/niel/damage_response_library.v1.json"
    )

    shard_directory = run_root / f"shard_{index:04d}"
    shard_directory.mkdir(parents=True, exist_ok=False)
    (shard_directory / "case.json").write_text(
        json.dumps(shard, indent=2) + "\n", encoding="utf-8"
    )

manifest = {
    "schema_version": "shieldSim.parallel_visualization_run.v1",
    "worker_count": worker_count,
    "histories_per_worker": histories_per_worker,
    "total_histories": total_histories,
    "seed_strategy": "base_seed_1020200_plus_worker_index",
    "merge_rule": (
        "event-weighted pooling of differential per-primary spectrum bins; "
        "workers have equal histories"
    ),
    "source_case": str(source_path),
    "purpose": (
        "higher-statistics review visualization; eight independent serial "
        "shieldSim processes because shieldSim has no MPI shard contract"
    ),
}
(run_root / "parallel_run_manifest.json").write_text(
    json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
)
print(run_root / "parallel_run_manifest.json")
PY
```

Setting `min_histories`, `max_histories`, and `batch_size` to the same value
forces each process to run exactly 12,500 histories; it disables an
accuracy-triggered early stop that would otherwise give unequal shard sizes.
Seeds `1020200` through `1020207` select distinct deterministic Geant4 random
streams. Unique output prefixes and working directories prevent file-name
collisions. `PYTHONDONTWRITEBYTECODE=1` prevents creation of `__pycache__` in
the source tree.

### 3. Launch and wait for all eight processes

This is the parallel-launch command used. `subprocess.Popen` starts every
worker before any `wait` call, so the processes execute concurrently rather
than sequentially. Standard output and errors are kept in a separate
`run.log` for each shard.

```bash
PYTHONDONTWRITEBYTECODE=1 python3 - <<'PY'
import json
import os
import subprocess
import time
from pathlib import Path

executable = Path(os.environ["SHIELDSIM_EXE"]).resolve()
run_root = Path(os.environ["SHIELDSIM_RUN_ROOT"]).resolve()
workers = sorted(run_root.glob("shard_[0-9][0-9][0-9][0-9]"))
if len(workers) != 8:
    raise SystemExit(f"expected 8 shard directories, found {len(workers)}")

started = time.perf_counter()
processes = []
logs = []
for index, shard_directory in enumerate(workers):
    log = (shard_directory / "run.log").open("wb")
    logs.append(log)
    process = subprocess.Popen(
        [str(executable), "--case=case.json"],
        cwd=shard_directory,
        stdout=log,
        stderr=subprocess.STDOUT,
    )
    processes.append((index, process))

exit_codes = {str(index): process.wait() for index, process in processes}
for log in logs:
    log.close()

completion = {
    "worker_exit_codes": exit_codes,
    "wall_seconds": time.perf_counter() - started,
    "complete": all(code == 0 for code in exit_codes.values()),
}
(run_root / "parallel_completion.json").write_text(
    json.dumps(completion, indent=2) + "\n", encoding="utf-8"
)
print(json.dumps(completion, indent=2))
if not completion["complete"]:
    raise SystemExit("at least one shieldSim worker failed; do not pool outputs")
PY
```

For the recorded run all eight exit codes were zero and the elapsed wall time
was 86.104 seconds on the test host. This time is hardware- and cache-dependent
and is evidence for this run, not a guaranteed performance requirement.

### 4. Check the shard boundary before post-processing

The integrity check used after execution was:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 - <<'PY'
import json
import os
from pathlib import Path

run_root = Path(os.environ["SHIELDSIM_RUN_ROOT"])
documents = []
for shard_directory in sorted(run_root.glob("shard_[0-9][0-9][0-9][0-9]")):
    result_paths = list(shard_directory.glob("*_result.json"))
    if len(result_paths) != 1:
        raise SystemExit(f"expected one result in {shard_directory}")
    document = json.loads(result_paths[0].read_text(encoding="utf-8"))
    if document["events"] != 12_500:
        raise SystemExit(f"wrong history count in {result_paths[0]}")
    documents.append(document)

if len(documents) != 8:
    raise SystemExit(f"expected 8 result documents, found {len(documents)}")
if len({item["random_seed"] for item in documents}) != 8:
    raise SystemExit("random seeds are not distinct")
if len({item["source_hash"] for item in documents}) != 1:
    raise SystemExit("source hashes differ")
if len({item["geometry_hash"] for item in documents}) != 1:
    raise SystemExit("geometry hashes differ")
if len({item["physics_list"] for item in documents}) != 1:
    raise SystemExit("physics-list identities differ")
print("8 shards, 100000 histories: source, geometry, and physics identities match")
PY

if rg -n -i \
  'fatal exception|segmentation fault|terminate called|(^|[^a-z])error([^a-z]|$)' \
  "${SHIELDSIM_RUN_ROOT}"/shard_*/run.log; then
  echo "fatal/error text found in a worker log" >&2
  exit 1
fi
```

The JSON check prevents accidental pooling of different sources, geometries,
or physics settings and confirms both the total event count and independent
seeds. The `rg` command searches all worker logs for fatal/error signatures;
no match is the successful outcome.

### 5. Pool review data and render the figures

After the checks, the spectrum and LET text products were pooled by history
count. For each transport-spectrum bin, the per-primary columns were combined
as

```text
pooled bin = sum(worker histories * worker bin) / sum(worker histories).
```

The LET standard errors were not averaged. Each worker's event sum of squares
was reconstructed from its mean, standard error, and event count; those raw
moments were added and the standard error was recomputed for 100,000 events.
The same event-moment rule, including numerator/denominator covariance, was
used for the LET means and `H100/10`. The pooler first required identical bin
edges, source hash, geometry hash, physics list, and distinct seeds. Its inputs,
SHA-256 hashes, and pooling rule are recorded in
`pooled/pooling_manifest.json`. This was visualization/statistical
post-processing only: it did not combine eight result JSON files into a
purported canonical result.

The exact plot commands used on the resulting pooled text files were:

```bash
export SHIELDSIM_POOLED="${SHIELDSIM_RUN_ROOT}/pooled"

python3 scripts/plot_shield_spectra.py \
  "${SHIELDSIM_POOLED}/bon2020_proton_100000_pooled_spectra.dat" \
  --save --format pdf --no-show
python3 scripts/plot_shield_spectra.py \
  "${SHIELDSIM_POOLED}/bon2020_proton_100000_pooled_spectra.dat" \
  --save --format png --dpi 240 --no-show
python3 scripts/plot_shield_spectra.py \
  "${SHIELDSIM_POOLED}/bon2020_proton_100000_pooled_spectra.dat" \
  --save --format svg --no-show

PYTHONDONTWRITEBYTECODE=1 python3 scripts/plot_let_spectrum.py \
  "${SHIELDSIM_POOLED}/bon2020_proton_100000_pooled_let_spectrum.dat" \
  --summary-json "${SHIELDSIM_POOLED}/pooled_quantity_summary.json" \
  --output-prefix "${SHIELDSIM_POOLED}/bon2020_proton_100000_pooled_let_spectrum" \
  --formats pdf png svg --dpi 240
```

The three `plot_shield_spectra.py` calls create vector PDF, high-resolution
PNG, and vector SVG versions of the input/transmitted energy spectra.
`--no-show` selects batch operation without an X display. The LET command
creates the same three formats in one call; `--summary-json` adds the pooled
track-mean and dose-mean reference lines and the second panel exposes the
per-bin relative standard error.

The recorded pooled scalar results met the requested 10% relative standard
uncertainty for TID, track-mean LET, dose-mean LET, and `H100/10`. DDD and
`n_eq` remained `null`/`not_assessable`, not zero, because the broad source and
generated secondaries produced unsupported-species or out-of-table crossings.
No clamping, extrapolation, or silent omission was used. Sparse secondary-alpha
bins may still appear as separated segments; only 66 alpha shield-exit
crossings occurred in 100,000 histories, and drawing through empty bins would
misrepresent the Monte Carlo evidence.
