# SWCME validation

The validation suite is a standalone C++ executable that calls production
SWCME interfaces. It uses the repository's existing Make-based build approach;
no second build system or external test dependency is required.

## Directory layout

```text
test/
  core/       Common registry, runner, reporting, CFG01, and CFG02
  1d/         Tests specific to the production 1-D model
  3d/         Tests specific to the production 3-D model
  reference/  Reviewed reference data for future comparison tests
  profiles/   Inputs describing future validation sampling profiles
  output/     Generated executable and results (ignored by Git)
```

Every validation is linked into the single `output/test_swcme` executable.

## Build

From `srcSEP/swcme`:

```sh
make -C test clean all
```

From `srcSEP/swcme/test`, the equivalent command is `make clean all`.

## Command-line interface

```sh
./output/test_swcme                 # run all tests; same as --all
./output/test_swcme --all           # run all tests in registry order
./output/test_swcme --list          # list tests without executing them
./output/test_swcme --test CFG01    # run exactly CFG01
./output/test_swcme --test CFG02    # run exactly CFG02
./output/test_swcme --help          # usage, options, and examples
```

Options are mutually exclusive. Unknown test IDs, unknown options, a missing
value after `--test`, and combined modes return exit code 2. A completed test
run returns 0 when every mandatory test passes and 1 when any test fails.
Subcheck SKIPs identify unavailable production paths and do not by themselves
fail an otherwise successful test.

The registry in `core/test_swcme.cpp` is the only source for `--list`, `--all`,
and test lookup, so displayed and executed tests cannot silently diverge.

## Validation classifications

- `COMMON`: contracts shared by both models, with both public paths exercised.
- `1D`: behavior specific to `swcme1d`.
- `3D`: behavior specific to `swcme3d`.
- `1D<->3D`: direct equivalence or consistency between the two implementations.

The current registry uses `COMMON` and `1D`; future comparison tests may use
`1D<->3D` when their primary purpose is cross-model equivalence.

## CFG01 status

The intended CFG01 definition is **Configuration rejection and physical-range
validation**, classified `COMMON`. Production SWCME currently has no unified
configuration-validation API and does not return rejection diagnostics. The
1-D model silently clamps some inputs, while the 3-D model follows different
policies. Inventing a new rejection API is outside the present validation task.

CFG01 therefore temporarily retains the pre-existing physical-constants
baseline and exact 1-D/3-D constant-consistency checks. It explicitly reports
the unavailable configuration-rejection path as SKIP. The unit conversions
formerly mixed into CFG01 moved to CFG02, so no existing coverage was removed.
A future production configuration contract should allow CFG01 to be completed
and its temporary constants checks to move to a dedicated constants test.

## CFG02: unit and dimensional consistency

CFG02 is a `COMMON` test of the complete unit contract currently reachable
through production SWCME. Its purpose is to detect scale errors that can look
numerically plausible but substantially change CME kinematics, density,
magnetic field, Alfvén speed, and shock Mach number.

CFG02 covers:

- forward velocity conversion for 0, 1, 321, 400, 1500, and 3000 km/s;
- forward magnetic-field conversion for 0.1 through 1000 nT;
- forward number-density conversion for 0.01 through 10000 cm^-3;
- AU and nominal-solar-radius distance conversion, including AU/Rs derived from
  the adopted constants rather than an independently rounded ratio;
- the DBM inverse-length conversion from km^-1 to m^-1;
- the 3-D degree-derived cone default at 40, 90, and 180 degrees;
- finite extreme values intended to detect overflow, underflow, integer
  conversion, or hidden scale assumptions rather than event realism;
- explicit 1-D-versus-3-D comparisons for independently implemented paths;
- the closest available equivalent external/SI configuration checks; and
- an Alfvén-speed smoke test using production-converted `B` and `n`, the shared
  production SI Alfvén helper, and a separate scalar reference calculation.

The equivalent configuration uses 400 km/s, 5 nT, 5 cm^-3, 1 AU, 24 hours,
and 40 degrees. SWCME has no unified direct-SI `Params` alternative: radius and
time already enter evaluators in SI, while velocity, field, and density are
heliophysics-unit fields. CFG02 tests the closest real paths and marks the
missing direct-SI configuration comparison as SKIP rather than constructing a
fictitious API.

Likewise, production accepts time in seconds and does not expose an hours input
conversion. The 1-D API accepts angles in radians, and only the 3-D default
contains a degree-to-radian conversion. The 1-D Tecplot writers provide inverse
velocity, AU, and solar-radius output paths; CFG02 parses temporary writer output
to test those round trips at the writer's documented decimal precision. Other
inverse paths are not exposed and are printed as SKIP; test-side arithmetic is
never presented as a production round trip.

CFG02 prints every executed value, units, production result, independent
reference, factor, absolute and relative error, tolerance, and status. It also
prints structured counters and maxima for conversion error, round-trip error,
cross-model difference, and Alfvén-speed error.

### Tolerances

Exact SI decimal-prefix conversions use exact comparison when the selected
binary values permit it. Short floating-point calculation chains use
`64 * std::numeric_limits<double>::epsilon()`. AU and solar radius checks test
the documented adopted SWCME constants, not a competing astronomical
convention. Physical-constant tolerances remain documented in CFG01.

Analytical reference expressions must remain independent of the production
routine under test. A tolerance or reference value must never be relaxed merely
to turn a production mismatch into PASS. A CFG02 failure should be investigated
as a unit-contract or production-policy issue and reported before physics is
changed.

## Current SWCME unit contract

External `Params` values use:

- velocity in km/s;
- magnetic-field magnitude in nT;
- number density in cm^-3;
- launch radius in nominal solar radii;
- drag coefficient in km^-1;
- sheath/ejecta dimensions in AU normalized at 1 AU;
- cone angle in radians, despite the 3-D default being initialized from 40
  degrees; and
- temperature in kelvin.

Model evaluators and prepared states use SI meters, seconds, m/s, tesla, m^-3,
and radians. SWCME adopts `1 AU = 149597870700 m`, nominal solar radius
`6.957e8 m`, CODATA 2022 proton mass `1.67262192595e-27 kg`, vacuum
permeability `1.25663706127e-6 N/A^2`, exact Boltzmann constant
`1.380649e-23 J/K`, and the model solar-rotation convention `2.86533e-6 rad/s`.

The 1-D and 3-D modules still contain separate hard-coded decimal conversion
factors for km/s, cm^-3, nT, and km^-1 even though physical constants are now
centralized. CFG02 exercises and compares both paths. Consolidating these
factors is recommended technical debt, but this test does not automatically
change them because doing so could conceal the behavior it is meant to audit.
