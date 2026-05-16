# shieldSim Test Suite

This directory contains automated and semi-automated test scripts for the
`shieldSim` Geant4 shielding application.

The tests are organized in layers. Lower-numbered layers test basic software
functionality, while higher-numbered layers test physics behavior, numerical
convergence, and regression stability.

## Test Philosophy

The purpose of the test suite is to separate several different questions:

1. Does the code build?
2. Does the command-line interface work?
3. Are invalid inputs rejected cleanly?
4. Are the geometry, source, and scoring algorithms behaving correctly?
5. Are the physics results consistent with simple reference problems?
6. Are Monte Carlo results statistically converged?
7. Do future code changes preserve validated behavior?

Not every test is expected to be fast. Layer-1 tests should be quick and should
be run frequently. Higher-layer physics and regression tests may require many
events and should be run before major releases or science production runs.

---

## Directory Layout

Recommended structure:

```text
tests/
├── README.md
├── run_layer1_tests.sh
├── run_layer2_tests.sh
├── run_layer3_tests.sh
├── run_regression_tests.sh
├── data/
│   ├── mono_100MeV_proton.dat
│   ├── flat_test_spectrum.dat
│   └── zero_flux.dat
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

Examples include vacuum transmission, shield-exit scoring, monoenergetic beam
transport, beam versus isotropic source checks, and detector ordering tests.

### `run_layer3_tests.sh`

Physics verification tests.

Examples include stopping-power comparisons, dose normalization checks,
LET-spectrum tests, secondary neutron production, and physics-list comparison.

### `run_regression_tests.sh`

Longer reference-output tests.

These should compare important integrated quantities against previously saved
reference results, using tolerances appropriate for Monte Carlo statistics.

### `data/`

Small input spectra and other fixed test data.

These files should be version controlled if they are small and deterministic.

### `expected/`

Reference output summaries or metadata for regression tests.

Large raw Geant4 output files should generally not be committed unless they are
small and essential.

### `logs/`

Generated test logs.

This directory should not be version controlled, except for a `.gitignore` file.

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

Layer-2 tests should verify geometry, source sampling, and scoring logic.

Recommended tests:

1. Vacuum or near-vacuum transmission.
2. Thin-shield monoenergetic proton transmission.
3. Shield rear-face scoring.
4. Beam source direction.
5. Isotropic source angular distribution.
6. Isotropic source position distribution.
7. Detector/target ordering.
8. Material alias equivalence, for example `Al` versus `G4_Al`.

These tests are intended to catch implementation errors such as:

- incorrect coordinate transformations,
- lost particles at boundaries,
- incorrect source direction sampling,
- incorrect detector indexing,
- target-volume scoring mistakes.

---

## Layer-3 Tests

Layer-3 tests should verify physics behavior.

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

These tests may require many more events than Layer-1 tests.

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
tolerances. A useful starting point is a three-sigma tolerance for integrated
quantities.

---

## What Should Be Committed

Commit:

```text
tests/*.sh
tests/README.md
tests/data/*.dat
tests/expected/*.json
tests/expected/*.txt
```

Do not commit generated build directories, logs, or large raw output files.

Recommended `.gitignore` entries:

```text
build_layer1_tests/
build_layer2_tests/
build_layer3_tests/
layer1_test_logs/
layer2_test_logs/
layer3_test_logs/
tests/logs/
```

---

## Notes on Monte Carlo Testing

Monte Carlo tests are not exactly reproducible unless the random seed and all
Geant4 settings are fixed. For reliable regression testing, the application
should support:

```bash
--random-seed=<integer>
--output-prefix=<name>
```

Useful future diagnostic options include:

```bash
--mono-proton=<MeV>
--mono-alpha=<MeV>
--disable-protons
--disable-alphas
--dump-source-samples=<file>
--dump-exit-particles=<file>
--set-production-cut=<mm>
--max-step=<mm>
```

These options make it much easier to diagnose differences between runs.
