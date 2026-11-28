# Enceladus forward-model tests

Run the complete application-independent suite from `srcEnceladus`:

```bash
python3 test/run_tests.py
```

The compiler can be supplied as an argv vector for a site wrapper or alternate
toolchain:

```bash
python3 test/run_tests.py --cxx clang++
```

The same receipt's 45 registered results (`U00`–`U10` and the `R_*`
references) are imported by the repository-wide procedure, which uses its
configured C++ compiler. The local receipt additionally lists one
`S_<module>.<Class>.<method>` row per standalone test (described below); the
global owner validates those rows and turns any `S_` failure into a suite
error, but does not import them as registered identities:

```bash
cd /home/vtenishe/Mars2/AMPS
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test APPLICATION-ENCELADUS-FORWARD-001
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test APPLICATION-ENCELADUS-NATIVE-STARTUP-001
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test APPLICATION-ENCELADUS-E3-001
```

The first owner imports the 45 dependency-light checks below. The native owner
performs a clean configured AMPS build and imports six startup gates and six case/MPI/mesh groups. The
E3 owner retains all 13 long-run identities but defaults to explicit `SKIP`;
enable the optional `enceladus-e3` machine capability only on a host approved
for the 10,000-step MPI campaign. The global catalog and manifest are the
authorities for exact IDs and current campaign membership.

An operator may instead request every expensive, portability, distributed,
and qualification owner with one explicit global key; this also activates E3:

```bash
python3 src/test/run_tests.py --full-validation --jobs 8 \
  --config test/qualification-machine/passion.json
```

Run one case, with its own JSON/SVG receipt, by its stable ID:

```bash
python3 test/run_case.py U06_source_law
python3 test/run_case.py U08_validation
```

`python3 test/run_case.py --help` lists all exporter, compile, parser,
chemistry, source-law, numerics, and negative-validation cases.

The exporter owner also checks automatic per-fit source snapshots, exact
nonuniform and zero strengths for all four geometries, exclusion of nuisance
backgrounds, deferred apparent-CO2 absolute exports, preservation of fitted
kinematics, and independent per-observation background recovery. The new
`U03c_catalog_source_geometry` runs 98-, 41- and 99-source absolute/shape
files through the actual dependency-light C++ reader, then corrupts source
membership and requires rejection. These are synthetic contract fixtures;
real catalog coordinates and source-rate inference require campaign evidence.
The same case requires exact preservation of a finite subnormal binary64
strength and rejection of overflow and underflow to zero; acceptance alone
cannot hide an accidental zeroing of a tiny fitted coefficient.
`U06_source_law` additionally compares negative shifted-Maxwellian moments
and quantiles against independent integration through shift -8. Its original
analytic-range references and tolerances are retained. No fitting result is
used as a numerical expected value.

Select the application through the repository workflow (never by editing
`.amps.conf`), confirm no Mars2 AMPS build/run is active, remove only the
generated `build` directory, regenerate, and compile before validating the
pre-MPI startup gates:

```bash
cd /home/vtenishe/Mars2/AMPS
rm -rf -- build
./Config.pl -application=enceladus
./ampsConfig.pl -input enceladus.input -no-compile
env MAKEFLAGS="-j16" make amps
```

```bash
cd /home/vtenishe/Mars2/AMPS/srcEnceladus
python3 test/run_native_startup.py \
  --executable /home/vtenishe/Mars2/AMPS/amps
```

This checks CLI help, validation of the complete checked-in E3 deck, explicit rejection of a
deliberately corrupted inward copy of the projected stripe receipt, pre-MPI
rejection of an inaccessible configured plasma-cache directory using the
production projection, and final generated pointer-mover and loss-callback
macros. These tests do not initialize MPI, allocate a mesh, or claim an
evolving physics run.

The global runner executes every test of the standalone plume model
(`standalone/tests/test_<module>.py`) and reports each test as its own row
`S_<module>.<Class>.<method>`, for example
`S_physics.PhysicsTests.test_inverse_square_and_horizon`. Remove `S_` and put
back `test_` to rerun one test from `srcEnceladus`:
`python3 -m unittest standalone.tests.test_physics.PhysicsTests.test_inverse_square_and_horizon`.
The four exporter/publisher modules (`test_export_forward_sources`,
`test_fit_source_export`, `test_prepare_forward_cases`,
`test_export_forward_observations`) run together in the caller's environment;
the registered case `U00_exporter` passes when none of their tests fails, and
their per-test rows follow it. Every other module runs in its own interpreter
and reads the Cassini validation data and, for the `S_spice_ground_tracks.*`
tests, CSPICE and its kernels:

```bash
python3 srcEnceladus/test/run_tests.py --output /tmp/enceladus-tests \
  --validation-data-root /absolute/path/to/cassini_validation_data \
  --spice-path /absolute/path/to/cspice --spice-kernels /absolute/path/to/Kernels
```

Each option defaults to its environment variable (`ENCELADUS_VALIDATION_DATA_ROOT`,
`ENCELADUS_SPICE_TOOLKIT_ROOT`, `ENCELADUS_SPICE_KERNEL_ROOT`). Without them the
tests that need these inputs are skipped; each skipped test is reported as
`SKIP` with its reason. `test/standalone_test_receipt.py` records the outcomes
(its docstring gives the rules): a unittest error is `FAIL`, any failing
subtest fails its test, a test that never ran because a `setUpClass` or
`setUpModule` fixture was skipped or failed inherits that outcome, a failing
fixture that explains no test (such as `tearDownClass`) is its own `FAIL` row,
and a module that cannot be imported or has no tests is one `FAIL` row
`S_<module>`. The per-test JSON of each interpreter is kept in
`<output>/standalone_tests/`. The AMPS nightly passes the machine file's
`enceladus_validation_data`, `spice_path` and `spice_kernels`.

It then exports a complete 8-jet/7-stripe fixture, builds the dependency-free C++ core, and runs
focused source-parser, runtime-parser, Burger-chemistry, Tenishev flux-law,
canonical ID/geometry, radial stripe projection, parent-loss-only planning,
Courant/weight, conditional flyby-mesh-refinement, strict inactive-control, and dimensionless
velocity-shift population-control tests. The source-law case compares the production K3
moment and inverse CDF with direct composite-Simpson quadrature; it does not
reuse the production recurrence. The chemistry case verifies that Burger's
resonant channel contains both H2O and H2O+. It then compares named outputs with
`reference/forward_model_reference.json` using independent equations in the
Python runner. Production code never writes the expected values.

Results are `PASS`, `FAIL`, `SKIP`, or an explicit compilation/runner failure in
`test_output/srcEnceladus/test_receipt.json`. The adjacent
`reference_comparison.svg` plots relative errors and is generated with the
Python standard library, so visualization does not depend on matplotlib.

The fixture records a 252.100 km minimum and maximum radius for every exported
stripe node and the same explicit reference radius. The exporter unit test
independently recomputes every vertex norm; the native negative case proves
that a one-percent inward corruption is still rejected.

The exporter tests also validate both E3 and E5 H2O+CO2 startup manifests:
each species has all 15 equal strengths, the sums reproduce the selected H2O
and apparent-CO2 totals exactly, and the CO2 provenance retains the factor-2--5
fragmentation bracket and its derived physical-rate interval.

Run the bounded one-rank E3 integration and conservation gate in a new output
directory (the runner refuses to overwrite prior evidence):

```bash
python3 test/run_e3_simulation.py \
  --executable /home/vtenishe/Mars2/AMPS/amps \
  --launcher mpiexec \
  --output test_output/srcEnceladus/e3_first_run_qualified
```

The runner requires normal process completion, the exact E3/no-plasma receipt,
a positive global timestep, 800 completed steps and injection calls, exactly
80,000 injected particles of each primary, zero CO2 photolytic losses, no
plasma-output file, native finalization, and exact per-species conservation:

```text
injected = photolytic loss + surface loss + outer escape + resident
```

`--launcher direct` is available only for one-rank MPI singleton diagnostics.
Some restricted containers deny the local sockets needed by both launcher
forms; that is an MPI-environment failure before model startup, not permission
to skip the native run. Use an authorized host MPI environment and retain the
failed environmental receipt separately.

The 2026-11-01 qualified E3 receipt passes all nine gates, but this does not
qualify four-rank ownership, plasma/coupled chemistry, numerical convergence,
or observational accuracy. Those integration gates remain in
`development_and_validation_plan.md`.

For the refined 10,000-step run, the runner additionally requires
`flyby_mesh_refinement=on` in the native receipt and verifies
`cell_max_edge_m <= 20000` for every in-domain row in every Tecplot zone. This
is a direct check of the allocated mesh at the observation locations; it does
not infer resolution from the requested configuration.

## End-to-end Table tests: E3 INMS and UVIS 2010

`MakefileTest/Table` holds four end-to-end tests of the production
executable, documented in [`end_to_end/README.md`](end_to_end/README.md):
`Enceladus--E3-INMS` and `Enceladus--UVIS-2010` (nightly, 5000 steps) and
`Long--Enceladus--E3-INMS` and `Long--Enceladus--UVIS-2010` (long tier,
10000 steps, 8 ranks). Each loads the `tenishev2014_hybrid_8` source fitted
to its observation alone from `data/input/Enceladus/fitted_plume`, compares
the model with the observation in `data/input/Enceladus/observation_files`,
gates on completion, rank count, zone count, particle ledger and an exact
regression of its products (`output/test_Enceladus/<NAME>--*.ref_np8.gz`),
and plots the model-observation comparison without gating on it.

## Repeatable INMS and UVIS comparison run

The checked-in runtime deck requests the currently usable INMS flybys and UVIS
occultations through repeatable named `observation` blocks. E7 supplies a
calibrated apparent-CO2 density comparison; E13 remains model-only; E21
retains H2O/CO2 detector counts as display-only native context because no
validated scalar density response exists. Additional UVIS events without
qualified reductions/rays are deliberately not treated as comparison data.
To add a new data set, add another block; to omit one from a run, comment its complete
`#begin observation` through `#end observation` range. The parser contract,
units, required columns, ray geometry, and fail-closed validation rules are
documented in `../RUNTIME_INPUT.md`.

Regenerate the repository-released E7/E21 tables and their checksummed
provenance manifest from the validation archive. The archive location is a
required input (`--validation-data-root` or `ENCELADUS_VALIDATION_DATA_ROOT`);
the exporter refuses to run without one rather than assuming a machine path:

```bash
python3 standalone/export_forward_observations.py \
  --validation-data-root /absolute/path/to/cassini_validation_data
```

By default this rewrites `standalone/observations/`. To compare a rebuild with
the checked-in release first, add `--output-dir /new/empty/directory` and
compare the table SHA-256 values in the two `release_manifest.json` files.

Run the current E3 source distribution for ten 1,000-step sampling windows in
a fresh output directory:

```bash
python3 test/run_e3_simulation.py \
  --executable /home/vtenishe/Mars2/AMPS/amps \
  --input /home/vtenishe/Mars2/AMPS/input/enceladus.in \
  --output test_output/srcEnceladus/e3_10000_extended_observations \
  --launcher direct --timeout-s 10800 \
  --expected-iterations 10000 --expected-sample-length 1000 \
  --expect-volume-fields off
```

Each configured request produces a ten-zone Tecplot file named
`PT/plots/observation.<request>.dat`. Generate the final-window comparison
PNGs and the machine-readable summary with:

```bash
python3 test/plot_observation_comparisons.py \
  test_output/srcEnceladus/e3_10000_extended_observations/PT/plots/observation.*.dat \
  --output-dir test_output/srcEnceladus/e3_10000_extended_observations/PT/plots/comparisons \
  --summary test_output/srcEnceladus/e3_10000_extended_observations/PT/plots/comparisons/observation_comparison_summary.json \
  --expected-zones 10
```

The plotter overlays observations only when the input block declares a
response-qualified `comparison_species`. A request with
`comparison_species=none` is either model-only or, when the three
`context_*` keys are present, adds a third native-observable panel. Context is
never scored against model density and is labelled explicitly in the plot and
JSON receipt.
UVIS panels compare the simulated line-of-sight H2O column with the optical-
depth shape using one fitted nonnegative scale plus a signed baseline, and
separately report the AMPS peak/published-column ratio. A plotter `PASS` means
that all products have the requested zone count, stable row sequence, valid
schema, and a generated PNG; it is not a claim of observational agreement or
physical convergence.


## Case runner and runtime mesh regression

The case campaign uses explicitly synthetic geometry and nonuniform source
rates; it never substitutes them for the real fitted-data catalog. Four
stdlib-only publisher tests run inside U00. U09 drives the actual C++ resolver
through nine test methods, including relocation, hash changes, source symlink
escape, missing/unknown cases, absolute/shape totals and removal of the original
fitting data. U10 checks independent point/block geometry, shell cavities,
segment-to-box distances, polar refinement and compiled AMR limits.

| Global test ID | What is checked / why |
|---|---|
| ENCELADUS-U09-CASE-RUNTIME | Actual CLI/catalog parser preserves all 28 nonuniform rate vectors and rejects wrong or unavailable inputs; no silent source fallback. |
| ENCELADUS-U10-RUNTIME-MESH | Independent box/sphere/shell/cone and exact-distance fixtures verify refinement intersection and compiled-depth rejection. |
| ENCELADUS-N06-CASE-LOCAL | Native list/describe/local validation resolves all 28 measurement–basis cases, including both UVIS ray schemas, without mesh allocation. |
| ENCELADUS-N07-CASE-MPI | One/four ranks agree on content; rank-local mismatch and missing data fail collectively without hanging and without a particle buffer. |
| ENCELADUS-N08-CASE-TRANSPORT | Nine short native runs cover every source basis, H2O and CO2 plus four-rank ownership; each closes an exact 40-particle global injection ledger. |
| ENCELADUS-N09-CASE-UVIS | Both native UVIS ray operators advance and output all three fixture rays in each of two sampling windows with positive domain paths. |
| ENCELADUS-N10-RUNTIME-MESH | A noncubic domain with box/sphere/shell refinements produces actual observation cells within 80 km and the independently calculated global timestep. |
| ENCELADUS-N11-CASE-ERRORS | Native local validation rejects a domain cutting the surface, an unattainable cell target and a misspelled mesh key before allocation. |

Native groups execute through the existing serialized native owner:

```bash
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test APPLICATION-ENCELADUS-NATIVE-STARTUP-001
```

This owner creates its synthetic bundle in retained artifacts, performs the
required root/idle/clean/configure/hooks/top-level `-j16` build, and then runs
the unchanged six legacy gates plus the six groups above. It temporarily holds
and restores existing application test output during source generation.
`tools.mpiexec` and `execution.mpi_rank_flag` come from the machine profile;
missing MPI is reported as SKIP for MPI-dependent groups, never PASS. Native
runs use fresh subdirectories and retain individual command logs, manifest,
source snapshots and the application receipt. Results are imported into the
single global JSON/JUnit reports with stable IDs.

For an already built executable whose configured root contains a generated
fixture bundle:

```bash
python3 -B srcEnceladus/test/run_case_native_tests.py \
  --executable ./amps --fixture-deck /path/to/fixture/cases.in \
  --output /new/output/directory --launcher-json '["mpiexec"]'
```

The native script does not build/configure or alter fitted data. Its synthetic
sources and four-step windows check integration and bookkeeping; they do not
measure Cassini agreement, chemistry convergence or steady plume convergence.
New tests must be registered in the global manifest/case registry and
`src/test/TEST_CATALOG.md` in the same change; agents must regenerate the
catalog and its manual snapshot rather than leave undocumented case IDs.
