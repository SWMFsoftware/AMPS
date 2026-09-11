# SWCME validation

The validation suite is a standalone C++ executable that calls production
SWCME interfaces. It uses the repository's existing Make-based build approach;
no second build system or external test dependency is required.

## Directory layout

```text
test/
  core/       Common registry, runner, reporting, CFG01, CFG02, and DEN01
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
./output/test_swcme --test DEN01    # run exactly DEN01
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

## DEN01: Leblanc density normalization at the reference distance

DEN01 is a `COMMON` physics/implementation test of the defining normalization
property of the production Leblanc, Dulk & Bougeret (1998) density profile. It
does not merely compare the two models with each other: an independent scalar
reference evaluates

```text
n_L(r) = A r^-2 + B r^-4 + C r^-6,
A = 3.3e5, B = 4.1e6, C = 8.0e7,
```

where `r` is heliocentric radius in nominal solar radii and the result is in
`cm^-3`. The test calculates `r_ref_Rs = AU_m / Rs_m` from the adopted constants
rather than maintaining a rounded AU-to-solar-radius constant. For each model,
the intended normalization is `n(r) = S n_L(r)`, with
`S = n_ref / n_L(1 AU)` and a fixed production reference radius of one AU.

The positive reference-density fixtures are 1, 5, 8, 20, and 100 `cm^-3`.
Each fixture constructs a fresh production model and step cache, evaluates the
public 1-D radial or 3-D Cartesian field interface at one AU, and checks

```text
abs(n_model(1 AU) / n_ref - 1) < 1e-12.
```

The acceptance criterion is strict and must not be weakened to turn a mismatch
into PASS. Production evaluators return density in `m^-3`; the test reports in
`cm^-3` for direct comparison with the configured normalization while retaining
the actual SI evaluation path.

The 3-D coverage evaluates `(+AU,0,0)`, `(0,+AU,0)`, `(0,0,+AU)`, and a
normalized non-axis-aligned `(1,2,3)` direction to verify that the ambient
profile is radial. For every reference density, DEN01 also compares 1-D against
3-D density and inferred scale factors at roundoff level. The same reference
location is supplied both as the production AU and as `(AU_m/Rs_m)*Rs_m` to
check AU/solar-radius coordinate equivalence without repeating CFG02.

Both production implementations expose their prepared SI coefficients `C2`,
`C4`, and `C6`. DEN01 independently recovers the nominal `A`, `B`, and `C` from
those caches and verifies the coefficients and inferred normalization factor at
`64 * double epsilon`. That roundoff tolerance covers only the short
floating-point operation chain; it is not a relaxed density-physics tolerance.
The model APIs do not expose a separately named scale factor, so inference from
`C2` is the closest production diagnostic.

The 1-D and 3-D modules currently duplicate the Leblanc coefficient literals
and normalization arithmetic in their respective `prepare_step()` functions.
Both use the shared adopted AU (`149597870700 m`) and nominal solar radius
(`6.957e8 m`), and both cache scaled coefficients in a value-type `StepState`.
No coefficient or normalization factor is static, global, or mutable shared
state. DEN01 nevertheless executes an A=5, B=20, C=1 `cm^-3` construction and
re-evaluation sequence independently for both models. Existing model A and B
must retain their original results after later instances are constructed, and
their cached coefficients must remain unchanged. This guards against future
cross-instance contamination that would be blocking in parallel AMPS use.

The production APIs currently fix density normalization at one AU, so tests of
a configurable 0.5- or 2-AU normalization radius are reported as `SKIP` rather
than supported by a fictitious test-only API. Values at 0.5 and 2 AU are checked
only for finite positivity as optional diagnostics. Detailed verification of
the radial `r^-2 + r^-4 + r^-6` behavior belongs to DEN02, not DEN01.

DEN01 prints the reference radius, each nominal analytical term, unscaled
density, expected and inferred scale, production density, absolute error,
normalization residual, tolerance, direction coordinates, and every individual
PASS/FAIL. Its structured text metrics include fixture/pass/fail/skip counts,
maximum normalization residual, maximum 1-D/3-D difference, coefficient
mutation status, and cross-instance contamination status.

A DEN01 failure may indicate incorrect radial units or coefficient powers, an
incorrect normalization factor, a `cm^-3`/`m^-3` conversion error, inconsistent
1-D and 3-D implementations, mutable global/static normalization state, or
accidentally modified nominal coefficients. The independent analytical helper
must remain separate from the production routines under test; references and
tolerances must never be changed solely to make a production failure pass.

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

## PAR01-PAR03: corrected 3-D Parker magnetic field

The production 3-D Parker field now uses an explicit solar-rotation axis and a
local spherical basis at every evaluation point.  The azimuthal direction is

```text
e_phi = (Omega_hat x e_r) / |Omega_hat x e_r|
```

and the local pitch is proportional to
`|Omega_hat x e_r| = sin(theta_local)`.  This corrects the previous
implementation, which used `(Omega_hat x e_r) x e_r` (a meridional direction)
and one global `sin_theta` for the entire 3-D domain.

`Params::solar_rotation_axis` specifies the global solar-rotation axis and is
normalized when a `StepState` is prepared.  A zero or non-finite axis is
rejected because it cannot define a Parker azimuthal direction.  The existing
`Params::sin_theta` field is retained only for backward-compatible
normalization of `B1AU_nT`: it specifies the reference `sin(colatitude)` at
which `B1AU_nT` is interpreted as total field magnitude at 1 AU.  It no longer
sets the local 3-D Parker winding.  The step cache stores `k_AU = Omega*AU/Vsw`;
each point multiplies this by its own `r_AU*sin(theta_local)`.

The pole is handled analytically.  When the radial direction is parallel to
the rotation axis, `sin(theta_local)=0`, so `B_phi=0` and the field is purely
radial.  The implementation does not manufacture an arbitrary azimuthal unit
vector at the coordinate singularity.

The following deterministic 3-D tests were added:

- `PAR01` — equatorial Parker-vector orientation.  With the rotation axis along
  +Z and the point at +X, it requires the outward radial component to lie along
  +X, the Parker azimuthal component to lie along -Y for the current polarity
  convention, and the meridional Z component to vanish.  It also verifies the
  requested total-field normalization at the reference latitude.
- `PAR02` — arbitrary latitude and arbitrary rotation axis.  It compares the
  production vector with an independent analytical Parker reference away from
  the equator and repeats the check after rotating the solar axis, preventing a
  hard-coded +Z implementation from passing accidentally.
- `PAR03` — polar-limit regularity.  It verifies the exact north and south
  rotation poles and a near-pole point, requiring the transverse field to tend
  continuously to zero without NaN, Inf, or an arbitrary transverse direction.

The analytical reference in `3d/test_parker.cpp` is intentionally independent
of the production Parker helper.  Component comparisons use roundoff-level
absolute tolerances scaled to the expected Tesla magnitude.  These tests do not
cover numerical solenoidality or field-line tangent/path-length validation;
those remain separate planned tests (`PAR04` and later) so the individual
validation requirements remain diagnostically focused.
