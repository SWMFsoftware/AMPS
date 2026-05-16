# shieldSim Test Suite

This directory contains automated and semi-automated test scripts for the
`shieldSim` Geant4 shielding application.

The tests are organized in layers. Lower-numbered layers test basic software
functionality, while higher-numbered layers test geometry, source sampling,
scoring, physics behavior, numerical convergence, and regression stability.

## Test Philosophy

The purpose of the test suite is to separate several different questions:

1. Does the code build?
2. Does the command-line interface work?
3. Are invalid inputs rejected cleanly?
4. Are the geometry, source, and scoring algorithms behaving correctly?
5. Are the physics results consistent with simple reference problems?
6. Are Monte Carlo results statistically converged?
7. Do future code changes preserve validated behavior?

Layer-1 tests should be quick and should be run frequently. Layer-2 tests are
also intended to be lightweight, but they run short Geant4 simulations and use
diagnostic output to check geometry/source/scoring behavior. Higher-layer
physics and regression tests may require many events and should be run before
major releases or science production runs.

---

## Directory Layout

```text
tests/
├── README.md
├── run_layer1_tests.sh
├── run_layer2_tests.sh
├── run_layer3_tests.sh              # planned
├── run_regression_tests.sh          # planned
├── data/
│   ├── mono_100MeV_proton.dat
│   └── mono_100MeV_alpha.dat
├── expected/
│   └── README.md
└── logs/
    └── .gitignore
```

### `run_layer1_tests.sh`

Build and smoke tests.

This script verifies that the code builds, that basic CLI commands work, and
that invalid input is rejected cleanly.

### `run_layer2_tests.sh`

Geometry, source, and scoring tests.

This script verifies beam-source sampling, isotropic-source sampling,
near-vacuum transmission through the shield rear face, target ordering, and
basic material alias behavior.

### `data/`

Small deterministic input spectra for tests.  These files are intentionally
simple and should be version controlled.

### `expected/`

Reference output summaries or metadata for future regression tests.  Large raw
Geant4 output files should generally not be committed unless they are small and
essential.

### `logs/`

Generated test logs.  This directory should not be version controlled except for
its `.gitignore` file.

---

## Diagnostic CLI Options Used by Tests

Layer-2 tests require several diagnostic options implemented in the main code:

```bash
--random-seed=<integer>
--output-prefix=<name>
--dump-source-samples=<file>
--dump-exit-particles=<file>
--diagnostic-max-rows=<n>
```

### `--random-seed=<integer>`

Sets the CLHEP/Geant4 random seed before the run.  This makes source-sampling
and short Monte Carlo tests reproducible.

### `--output-prefix=<name>`

Changes standard output file names.  With the default prefix `shieldSim`, files
keep their historical names, such as:

```text
shieldSim_spectra.dat
shieldSim_quantities.dat
shieldSim_let_spectrum.dat
shieldSim_dose_sweep.dat
shieldSimOutput*.csv
```

With `--output-prefix=caseA`, they become:

```text
caseA_spectra.dat
caseA_quantities.dat
caseA_let_spectrum.dat
caseA_dose_sweep.dat
caseAOutput*.csv
```

### `--dump-source-samples=<file>`

Writes one row per generated primary:

```text
row species E_MeV x_mm y_mm z_mm ux uy uz
```

Positions are global coordinates in millimeters.  The vector `(ux,uy,uz)` is the
unit momentum direction requested by `PrimaryGeneratorAction`.

This file is used to verify that:

- beam mode generates `x = y = 0` and direction `+z`,
- isotropic mode keeps positions on the upstream plane,
- isotropic mode samples the inward-hemisphere crossing distribution
  `p(mu) = 2 mu`, where `mu = uz`.

### `--dump-exit-particles=<file>`

Writes one row for each particle accepted as crossing the downstream shield
face:

```text
row species E_MeV xg_mm yg_mm zg_mm xl_mm yl_mm zl_mm uxl uyl uzl
```

The `g` coordinates are global coordinates.  The `l` coordinates and directions
are in the local coordinate frame of the shield.  This file is used to verify
that shield rear-face scoring is based on the correct local-coordinate crossing
condition.

### `--diagnostic-max-rows=<n>`

Limits the number of rows written to each diagnostic file.  This prevents
accidentally generating huge diagnostic dumps.  Use a value less than or equal
to zero to remove the explicit cap.

---

## Layer-1 Tests

Layer-1 tests are software-level build and CLI smoke tests.

They check:

- clean CMake configuration,
- clean build,
- executable discovery,
- `--help`,
- `--list-materials`,
- `--list-target-materials`,
- `--list-detector-materials`,
- `--list-quantities`,
- rejection of invalid shielding material,
- rejection of invalid physics list,
- rejection of invalid source mode,
- rejection of negative target thickness,
- rejection of invalid computed quantity,
- rejection of malformed shield specification,
- rejection of negative event count.

Run from the top-level package directory:

```bash
chmod +x tests/run_layer1_tests.sh
tests/run_layer1_tests.sh
```

To see script options:

```bash
tests/run_layer1_tests.sh --help
```

Optional environment variables:

```bash
BUILD_DIR=build_test JOBS=8 tests/run_layer1_tests.sh
```

Available variables:

```text
BUILD_DIR   Build directory used by the test script.
            Default: build_layer1_tests

EXE_NAME    Executable name.
            Default: shieldSim

LOG_DIR     Directory where test logs are written.
            Default: layer1_test_logs

JOBS        Number of parallel build jobs.
            Default: nproc, or 2 if nproc is unavailable.
```

Expected result:

```text
All Layer-1 tests passed.
```

---

## Layer-2 Tests

Layer-2 tests verify geometry, source sampling, and scoring logic using short
Geant4 runs and diagnostic files.

The implemented Layer-2 script checks:

1. The diagnostic options are present in `--help`.
2. Beam mode generates 100 MeV protons at `x = y = 0` with direction `+z`.
3. Isotropic mode samples the source plane and direction distribution correctly:
   - all directions point inward, `0 < mu <= 1`,
   - `<mu> ≈ 2/3`,
   - `<mu^2> ≈ 1/2`,
   - source positions stay within the finite upstream plane.
4. A near-vacuum shield transmits most 100 MeV protons through the downstream
   shield face.
5. Exit diagnostics use shield-local coordinates and positive downstream local
   direction cosine.
6. Detector/target material order is preserved in output metadata.
7. `Al` and `G4_Al` both run successfully as shield material specifications.

Run from the top-level package directory:

```bash
chmod +x tests/run_layer2_tests.sh
tests/run_layer2_tests.sh
```

To see script options:

```bash
tests/run_layer2_tests.sh --help
```

Optional environment variables:

```bash
BUILD_DIR=build_l2 JOBS=8 tests/run_layer2_tests.sh
```

Available variables:

```text
BUILD_DIR   Build directory used by the test script.
            Default: build_layer2_tests

EXE_NAME    Executable name.
            Default: shieldSim

LOG_DIR     Directory where test logs are written.
            Default: layer2_test_logs

RUN_DIR     Directory where each Layer-2 test case writes its run outputs.
            Default: layer2_test_runs

JOBS        Number of parallel build jobs.
            Default: nproc, or 2 if nproc is unavailable.

SKIP_BUILD  Set to 1 to skip CMake configure/build.
            Default: 0

EXE_PATH    Explicit executable path when using SKIP_BUILD=1.
```

Example using an already-built executable:

```bash
SKIP_BUILD=1 EXE_PATH=build/shieldSim tests/run_layer2_tests.sh
```

Expected result:

```text
All Layer-2 tests passed.
```

Layer-2 test outputs are written under:

```text
layer2_test_runs/
```

Logs are written under:

```text
layer2_test_logs/
```

These generated directories should not be committed.

---

## Layer-3 Tests

Layer-3 tests should verify physics behavior.  They are planned but not yet
implemented in this package.

Recommended tests:

1. Proton stopping power in thin Al, water, and Si targets.
2. Alpha stopping power.
3. CSDA range behavior as a function of shield thickness.
4. TID unit normalization.
5. Source-intensity scaling.
6. LET spectrum behavior.
7. DDD and `n_eq` proxy sanity checks.
8. Physics-list comparison:
   - `FTFP_BERT`
   - `FTFP_BERT_HP`
   - `Shielding`
   - `QGSP_BIC_HP`

These tests may require many more events than Layer-1 or Layer-2 tests.

---

## Regression Tests

Regression tests should be created after a set of reference results is accepted.

Recommended regression cases:

```text
regression_001_beam_100MeV_p_Al2mm_Si1mm
regression_002_isotropic_SEP_Al2mm_BFO50mm
regression_003_GCR_HDPE10mm_Si1mm
regression_004_physicslist_comparison_Al20mm
regression_005_material_catalog_all_targets
```

For each regression case, compare quantities such as:

```text
total number of events
transmitted proton integral
transmitted alpha integral
transmitted neutron integral
TID
DDD
n_eq
H100/10
LET-spectrum integral
```

Because this is a Monte Carlo code, comparisons should use statistical
tolerances.  A useful starting point is a three-sigma tolerance for integrated
quantities.

---

## What Should Be Committed

Commit:

```text
tests/*.sh
tests/README.md
tests/data/*.dat
tests/expected/*.md
tests/expected/*.json
tests/expected/*.txt
tests/logs/.gitignore
```

Do not commit generated build directories, logs, test-run directories, or large
raw output files.

Recommended `.gitignore` entries:

```text
build_layer1_tests/
build_layer2_tests/
build_layer3_tests/
layer1_test_logs/
layer2_test_logs/
layer3_test_logs/
layer2_test_runs/
layer3_test_runs/
tests/logs/
```

---

## Notes on Monte Carlo Testing

Monte Carlo tests are not exactly reproducible unless the random seed and all
Geant4 settings are fixed.  For reliable regression testing, use:

```bash
--random-seed=<integer>
--output-prefix=<name>
```

Useful future diagnostic or convenience options include:

```bash
--mono-proton=<MeV>
--mono-alpha=<MeV>
--disable-protons
--disable-alphas
--set-production-cut=<mm>
--max-step=<mm>
```

The current package already implements source and exit diagnostic dumps, which
are enough for the Layer-2 source/geometry/scoring tests.
