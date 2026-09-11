# SWCME validation

The validation suite is a standalone C++ executable that calls production
SWCME interfaces. It uses the repository's existing Make-based build approach;
no second build system or external test dependency is required.

## Directory layout

```text
test/
  core/       Common registry, runner, reporting, CFG/DEN/KIN/SHK tests
  1d/         Tests specific to the production 1-D model
  3d/         Tests specific to the production 3-D model
  reference/  Reviewed reference/event data used by comparison tests
  profiles/   SMOKE/ROUTINE/FULL/EVENT test selections
  python/     Regression tests for the campaign manager
  output/     Generated executables and campaign artifacts (ignored by Git)
  run_tests.py  Python campaign manager
  event_config.example.json  executable EVENT/sweep example
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
./output/test_swcme --test PST01    # prepared-state immutability
./output/test_swcme --test PST02    # prepared-state ownership rejection
./output/test_swcme --test PST03    # configuration-state ownership rejection
./output/test_swcme --test PST06    # prepared-record integrity rejection
./output/test_swcme --test PST04    # concurrent prepared-state evaluation
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

## PST01: prepared-state immutability

`PST01` verifies that a successfully prepared 1-D or 3-D state is independent
of every later attempt to modify its model configuration.  The model remains
configurable during setup, but its first successful `prepare_step()` freezes
the configuration.  Failed preparation leaves the setup phase unlocked so an
invalid model can be corrected and retried.

The 1-D fixture evaluates density, velocity, Parker-field components, field
magnitude, and divergence at three fixed radii and serializes all 18 values.
It then attempts every retained legacy mutation path:

- `SetParams`, `SetCME`, and both kinematics setters;
- `SetAmbient`, including the historically defective `sin_theta` case;
- region and shock-acceleration mode setters;
- geometry, smoothing, and sheath/ejecta setters;
- copy assignment.

Each attempt must throw `std::logic_error` with the immutable-lifecycle
diagnostic before changing the configuration digest or model identity.  The
same state is reevaluated after every attempt, and its serialized result must
be bitwise identical to the baseline.

Compile-time API inspection separately requires the legacy 1-D
`MutableParams()` raw-reference escape hatch to be unavailable.  Merely checking
that accessor at call time would be insufficient because a reference obtained
before preparation could be retained and used after the model freezes.

The 3-D fixture uses compile-time inspection to confirm that no mutable Params
accessor exists, then verifies that copy assignment is rejected after
preparation.  Density, all velocity and magnetic-field components, and
divergence at three Cartesian points remain bitwise identical.  Both fixtures
exercise `reconfigured(params)`: it must create a new owner and configuration
snapshot, permit independent preparation and different physics, and leave the
original state unchanged.

Run the gate directly with:

```sh
./output/test_swcme --test PST01
```

`PST01` is the first test in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it
through their `@ALL` expansion.

## PST02: cross-model prepared-state rejection

`PST02` verifies that each prepared state belongs to the exact model instance
that created it.  This is intentionally stronger than configuration equality:
two separately constructed models with identical parameters receive different
process-local identities, and neither accepts the other's `StepState`.

The fixture prepares a state with Model A and attempts to consume it with Model
B using equal configurations and with Model C using a different solar-wind
configuration.  It covers:

- checked 1-D background, full-field, shock-source, and Tecplot-writer paths;
- checked 3-D background, magnetic-field, divergence, shock, acceleration, and
  Tecplot bundle paths;
- legacy 3-D geometry, mesh, and connectivity wrappers; and
- 1-D/3-D AMPS-facing background, point source, surface source, and cobpoint
  source adapters.

Every checked call must return `STATE_MODEL_MISMATCH` with
`has_model_identities=true`, the receiver's identity in
`expected_model_identity`, and Model A's identity in
`supplied_model_identity`.  Numerical arrays and source/background objects are
initialized with sentinels and must remain unchanged.  Pre-existing 1-D and
3-D output files contain fixed byte strings and must remain byte-for-byte
identical, proving rejection occurs before `fopen("w")`.  Legacy wrappers must
throw an exception whose diagnostic retains `STATE_MODEL_MISMATCH`; they may
not turn ownership misuse into a valid radius, mesh, or disconnected result.

Run the gate directly with:

```sh
./output/test_swcme --test PST02
```

`PST02` is included in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it through
their `@ALL` expansion.  A default-constructed state has owner identity zero
and is rejected by the same contract.

## PST03: cross-configuration state rejection

`PST03` verifies that prepared-state provenance protects the complete resolved
configuration, not only the `gamma_ad` value that originally exposed the
defect.  `prepare_step()` records an allocation-free deterministic digest in
the state.  Every state-consuming API compares the receiver's current digest
before it evaluates physics or modifies caller-owned output.

The test prepares one reviewed 3-D baseline, then constructs receivers that
differ by exactly one field representing each validation-plan family:

- adiabatic index (`gamma_ad`);
- geometry (`axis_ratio_y`, including inactive-field coverage);
- Parker orientation (`solar_rotation_axis`);
- kinematics mode;
- thermal-closure input (`T_K`);
- Parker field/polarity normalization input (`B1AU_nT`); and
- region mode.

Parker radial polarity and the proton-only pressure closure are compile-time
resolved conventions rather than independently mutable `Params` fields.  They
are therefore included as explicit versioned digest tags; `B1AU_nT` and `T_K`
exercise their runtime normalization/closure inputs.  A multi-field receiver
checks that no parameter-specific rejection branch exists.

For every foreign receiver, the checked evaluator must return
`STATE_MODEL_MISMATCH` before writing its sentinel outputs.  PST02 ownership has
intentional precedence, while `has_configuration_digests=true` and the
receiver/prepared digest pair prove the additional configuration mismatch.
An explicitly populated default-equivalent receiver must have the same digest
as the baseline yet remain rejected as a foreign model.

Because PST01 rejects supported post-prepare model mutation, PST03 deliberately
alters the digest tag only on a copied negative-test state.  The same-owner
consumer must return `STATE_CONFIGURATION_MISMATCH`, carry current and supplied
digests, and leave density/velocity sentinels unchanged.  A reviewed golden
digest plus an independently constructed equal configuration guards
reproducibility across runs and compiler rebuilds.

Run the gate directly with:

```sh
./output/test_swcme --test PST03
```

`PST03` follows `PST02` in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it
through their `@ALL` expansion.

## PST06: prepared-state record integrity

`PST06` verifies that model ownership and configuration equality cannot be
bypassed by modifying a cached `StepState` value.  Each successful preparation
stores a private seal produced by explicit field serialization.  The seal
covers canonical common solar-wind/kinematic state, region boundaries,
acceleration configuration, 1-D shock/RH data, 3-D frame and geometry caches,
and every top-level legacy compatibility mirror.

Compile-time assertions require `integrity_digest()` to return by value rather
than expose a mutable reference.  They also require both state types to remain
copy- and move-constructible.  Runtime fixtures exercise copy construction,
move construction, copy assignment, and move assignment; valid copies retain
the original owner, seal, and evaluability.

The negative matrix copies a valid state and changes exactly one record field
at a time.  It includes every top-level 1-D and 3-D mirror plus representative
members of every nested canonical record.  Each checked evaluator must return
`STALE_PREPARED_STATE`, set `has_state_integrity`, report the prepared and
recomputed seals, and leave all numerical output sentinels unchanged.  A legacy
1-D wrapper must throw with the same status name.  Separate owner/configuration
tag tests remain in PST02/PST03 because those higher-precedence diagnostics are
intentional.

Run the gate directly with:

```sh
./output/test_swcme --test PST06
```

`PST06` follows PST03 in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it
through their `@ALL` expansion.

## PST04: concurrent prepared-state evaluation

`PST04` verifies the production usage pattern required by AMPS: configuration
and `prepare_step()` finish on one setup thread, after which multiple workers
read the same model and the same immutable prepared state.  It does not grant
permission to call model mutators, prepare a new state concurrently, share
caller-owned destination arrays, or invoke file writers from multiple threads.

One 1-D and one 3-D state are each prepared once.  The fixture then exercises
nine operation families:

- scalar 1-D and 3-D AMPS background queries;
- direct 1-D and 3-D full-field batch evaluators;
- 1-D shock-source conversion;
- 3-D directional shock and source evaluation;
- 3-D Parker-line connectivity with complete cobpoint shock records; and
- simultaneous, intentionally different 1-D/3-D failures that retain their
  own status code, context, sample index, offending-value flag, and outputs.

Every public result is serialized field by field into fixed-width words.  The
serializer includes every `ModelStatus` field, all background/source/shock
numbers and flags, every connectivity root, and all batch output elements.  It
does not compare raw structure memory, so unspecified padding cannot cause a
false race diagnosis.  Each concurrent result must be bitwise identical to a
serial oracle; no roundoff allowance is currently required.

The workload runs eight repetitions at 1, 2, 4, and 8 threads under
forward-interleaved, reverse-interleaved, and operation-grouped job orders.  A
one-shot start gate creates overlap, while each output slot has exactly one
writer.  The test framework is used only after joining the workers, preventing
the harness itself from introducing a race.  Finally, both prepared-state
integrity seals are recomputed to prove the stress run did not mutate them.

Run the normal concurrency gate with:

```sh
./output/test_swcme --test PST04
```

When the compiler and runtime support ThreadSanitizer, rebuild the entire
executable with instrumentation so both the test and production implementation
are observed:

```sh
make clean
make CXXFLAGS="-O1 -g -std=c++17 -Wall -Wextra -Wpedantic -fsanitize=thread -fno-omit-frame-pointer"
TSAN_OPTIONS=halt_on_error=1 ./output/test_swcme --test PST04
```

A supported ThreadSanitizer run must report no race.  Some container kernels
cannot initialize ThreadSanitizer's shadow memory; that platform limitation is
not a passing sanitizer result and should be recorded as unavailable.  `PST04`
follows `PST06` in `SMOKE`; `ROUTINE`, `FULL`, and `EVENT` include it through
their `@ALL` expansion.

## Python campaign manager and reproducible run artifacts

`run_tests.py` is the manager for multi-test development gates and validation
campaigns.  It does not contain independent SWCME physics: each C++ validation
is still executed through the single `output/test_swcme` registry executable.
Builds also produce `output/sep_reference`, a public-interface consumer used for
reference histories and sweep probes.

Typical commands from `swcme/test` are:

```sh
python3 run_tests.py --list
python3 run_tests.py --profile SMOKE
python3 run_tests.py --profile ROUTINE
python3 run_tests.py --all
python3 run_tests.py --test SEP03
python3 run_tests.py --profile EVENT --event-config event_config.example.json
python3 python/test_run_tests.py
```

The named profiles are version-controlled text files in `profiles/`:

- `SMOKE` is the short development gate;
- `ROUTINE` is the broad deterministic gate and excludes `MSH05` and `CON05`;
- `FULL` expands to every test in the C++ registry;
- `EVENT` first runs FULL and then executes the supplied event-analysis JSON.

`--no-build` is available when the executables are already current;
`--build-only` compiles the test and reference tools without running a campaign.
`--reference-export` forces the default SEP history export and
`--no-reference-export` disables it. FULL and EVENT export it by default.

Each campaign creates an isolated output directory containing:

```text
manifest.json
summary.json
summary.csv
logs/<TEST_ID>.log
sep_reference.csv              # when reference export is enabled
sweeps/*.csv + per-case logs   # when configured
convergence/*.json             # when configured
comparisons/*.json             # when configured
plots/*                        # optional
```

The manifest records the random seed, selected tests, git metadata, compiler
path/version/flags, host and Python details, visible MPI/OpenMP environment,
source-tree SHA-256, complete resolved SWCME+SEP model configuration and hash,
and the event JSON/hash.  Use `--resolved-config FILE` when a driver has already
written the exact resolved parameter block; otherwise the manager obtains the
default block from `sep_reference --print-manifest`.

Campaign exit codes are part of the automation contract:

```text
0  all required tests/analyses passed
1  validation, reference export, or EVENT analysis failed
2  command-line/configuration error
3  build failure
```

### EVENT JSON

`event_config.example.json` is executable and demonstrates a parameter sweep.
The event object may contain `sweeps`, `convergence`, `comparisons`, and `plots`.
A sweep defines a Cartesian product of parameter arrays and a command template.
Any command token can use a sweep field plus `{seed}`, `{root}`, `{test_dir}` or
`{output_dir}`.  `metric_regex` captures a scalar from stdout; optional
`metric_min`/`metric_max` turn the captured quantity into an acceptance gate.

A convergence entry either supplies explicit positive `x` and `error` arrays or
references a previous sweep through `sweep`, `x_parameter`, and `error_metric`.
The manager fits the slope of `log(error)` versus `log(x)` and checks optional
`min_order`/`max_order` limits.

A comparison entry names model/reference CSV files, optional key/key tolerance,
and one or more columns with absolute and/or relative tolerances.  A plot entry
selects a CSV, one x column and one or more y columns.  Plotting is optional;
when Matplotlib is unavailable a non-required plot is reported as SKIP rather
than invalidating the physics campaign.

The Python manager itself is regression-tested in `python/test_run_tests.py`.
Those tests cover profile expansion against the live C++ registry, exact
second-order convergence fitting, Cartesian sweeps/metric capture, CSV
comparison and EVENT aggregation, manifest/report creation, and deterministic
CLI configuration-error handling.

## Validation classifications

- `COMMON`: contracts shared by both models, with both public paths exercised.
- `1D`: behavior specific to `swcme1d`.
- `3D`: behavior specific to `swcme3d`.
- `1D<->3D`: direct equivalence or consistency between the two implementations.

The current registry uses all four classifications.  In particular, the
`1D3D01`-`1D3D03` and `SEP03` gates use direct dimensional-equivalence fixtures
where appropriate.

## CFG01: centralized configuration rejection and physical-range validation

`CFG01` is a `COMMON` test of the production configuration boundary.  Both
model classes expose a side-effect-free `validate()` method returning
`swcme::config::ValidationResult`.  Every issue contains the public field name,
a stable error code, and the violated rule.  `prepare_step()` invokes the same
validator and throws `std::invalid_argument` before any physics if the result is
invalid, so callers cannot bypass validation accidentally.

The deterministic fixtures cover zero solar-wind speed, negative density,
negative DBM Gamma, invalid Parker normalization latitude, NaN temperature,
negative smoothing widths, malformed DATA_DRIVEN tables, zero CME direction,
zero solar-rotation axis, invalid finite-SSE half widths, non-positive
ellipsoid axis ratios, and negative solar-rotation rate.  A multiple-error
fixture confirms that one validation call reports all independent bad fields
instead of stopping after the first.  For representative invalid cases CFG01
also calls `prepare_step()` and requires rejection.

Validation is intentionally separate from unit conversion.  The common unit
layer converts `0 km/s -> 0 m/s` exactly; CFG01 separately rejects `V_sw<=0`
because the Parker/DBM baseline requires a positive ambient wind.  This is a
permanent regression guard against the former 1-D `max(1,V_sw*1000)` behavior.
The compact shared-constant baseline remains in CFG01 to guard the immutable
AU, nominal solar radius, proton mass, vacuum permeability, and Boltzmann
constant used by both interfaces.

## CFG02: centralized unit-conversion and dimensional-consistency test

`CFG02` validates `swcme_units.hpp`, the single production source for unit
conversion.  Forward and round-trip checks cover km/s, nT, cm^-3, km^-1, AU,
nominal solar radii, hours, and degrees.  Exact decimal/defined conversions are
checked exactly where possible; short floating-point chains use
`64 * std::numeric_limits<double>::epsilon()`.  No physics tolerance is used.

The test then constructs valid 1-D and 3-D models and confirms that their
prepared ambient speed, one-AU density, and one-AU magnetic-field magnitude
agree to roundoff.  An independent SI calculation of
`V_A=B/sqrt(mu0*rho)` is compared with the shared production Alfvén-speed
helper for both model paths.  Zero wind speed is deliberately *not* passed to
`prepare_step()` in CFG02: it is a valid conversion input but an invalid model
configuration and is therefore covered by CFG01.

The model wrappers now use the same conversion helpers when building their SI
state, so CFG02 no longer audits duplicated `1e3`, `1e6`, `1e-9`, or `1e-3`
conversion literals in 1-D versus 3-D.  A failure should be investigated as a
unit-contract or wrapper-integration defect; reference values or tolerances
must never be changed merely to obtain PASS.

## DEF01-DEF04: canonical defaults, scope, and resolved metadata

`swcme_defaults.hpp` owns every dimensionality-independent default used by the
1-D and 3-D public parameter structures.  These tests are release guards: a
change to a science default must be deliberate, documented, and changed in one
place rather than drifting independently between the two interfaces.

- `DEF01` compares all shared default ambient, kinematic, source, region,
  smoothing, sheath, and ejecta values between default-constructed 1-D and 3-D
  `Params`.  It also checks that the DBM reference radius and 1-AU density come
  from the canonical default constants.
- `DEF02` verifies the declared science convention: `SHOCK_ONLY + SOURCE` maps
  to `CONTROLLED_SEP_PRE_SHOCK`, the default 3-D geometry is finite `SSE`, the
  Parker normalization is total positive |B| at the equatorial 1-AU reference,
  radial polarity is outward, and the baseline frame name is stable.
- `DEF03` verifies observer-local scope gating.  A pre-shock observer is in the
  controlled Parker scope; after the modeled front reaches that observer the
  status is explicitly out of scope.  A late-time observer outside a finite SSE
  angular cap remains in scope because no front exists on that ray.
- `DEF04` verifies deterministic complete resolved-configuration serialization.
  Common fields, model-specific geometry, data-table counts, compatibility
  parameters, derived scope, and explicit event overrides must all appear.  The
  serialization is repeated byte-for-byte to guard deterministic campaign
  manifests/hashes.

Tests that exercise the optional `FULL_ICME/RESOLVED_COMPRESSION` path set both
options explicitly; they no longer inherit that behavior from defaults.  This is
important because a verification fixture should declare when it is testing an
optional phenomenological model rather than the controlled SEP baseline.

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

The Leblanc coefficient literals and normalization arithmetic now live only in
production `swcme_solarwind.hpp`.  Both dimensional wrappers store the same
`StepState::common.solar_wind` cache and merely mirror `C2/C4/C6` for backward
source compatibility.  DEN01 still independently reconstructs the published
coefficients from those mirrors and executes an A=5, B=20, C=1 `cm^-3`
construction/re-evaluation sequence. Existing model A and B must retain their
original results after later instances are constructed, guarding against any
future cross-instance contamination in parallel AMPS use.

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

## CFG01: centralized configuration rejection and physical-range validation

`CFG01` now exercises the production validation API rather than reporting the
configuration layer as a skip.  Both model classes expose a side-effect-free
`validate()` method returning `swcme::config::ValidationResult`.  Every issue
contains the public field name, an error code, and the violated rule.
`prepare_step()` calls this same validator and throws `std::invalid_argument`
before any physics if the result is invalid.

The deterministic fixtures cover zero solar-wind speed, negative density,
negative DBM Gamma, invalid Parker normalization latitude, NaN values, negative
smoothing widths, malformed DATA_DRIVEN tables, zero CME direction, zero solar
rotation axis, invalid finite-SSE half widths, non-positive ellipsoid axis
ratios, and negative solar rotation rate.  A multiple-error fixture confirms
that one validation call reports all independent bad fields rather than only
the first.  For representative invalid cases CFG01 also calls `prepare_step()`
and requires rejection, proving that no caller can bypass the validation layer
by ignoring `validate()`.

Validation is intentionally distinct from conversion.  For example, the unit
layer maps `0 km/s -> 0 m/s` exactly; CFG01 separately rejects `V_sw<=0` because
the baseline Parker/DBM model requires positive ambient wind speed.  This is a
permanent regression guard against the former 1-D `max(1,V_sw*1000)` behavior.

## CFG02: centralized unit-conversion and dimensional-consistency test

`CFG02` validates the production helpers in `swcme_units.hpp` for forward and
round-trip conversions of velocity, magnetic field, number density, inverse
length, AU, solar radii, hours, and angles.  Exact decimal/defined conversions
are checked at exact or roundoff precision; no physics tolerance is used.

The test then constructs valid 1-D and 3-D models and confirms that their
prepared solar-wind speed, 1-AU density, and 1-AU magnetic-field magnitude agree
to roundoff.  Finally, an independent SI calculation of
`V_A=B/sqrt(mu0*rho)` is compared with the shared production Alfven-speed helper
for both model paths.  The previous CFG02 failure at zero wind speed is removed
by design: zero is a valid unit-conversion input and an invalid model
configuration, so it is tested in the correct layer rather than by forcing a
model setup with inadmissible parameters.

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

Unit conversion is now centralized in production `swcme_units.hpp`; the 1-D
and 3-D `prepare_step()` paths no longer maintain independent km/s, cm^-3, nT,
or km^-1 conversion factors.  Conversion and physical validation are separate:
zero values remain zero under conversion, while `swcme_config.hpp` decides
whether a model parameter is admissible.  CFG02 verifies the common conversion
helpers directly and then checks that both dimensional wrappers produce the
same prepared SI state for valid physical inputs.

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

## GEO01-GEO08: corrected finite shock geometry, normals, and flank speed

The 3-D shock geometry has been reworked so that the production model no longer
uses the former `ConeSSE` approximation

```text
R(theta) = R_apex cos(theta)^m
```

with a radial normal and an artificial clamped radius beyond the configured
half width.  `ShockShape::SSE` now denotes a true finite self-similar-expansion
spherical cap.  `ShockShape::ConeSSE` is retained only as a source-compatible
enum alias and has the same corrected SSE semantics; `flank_slowdown_m` remains
in `Params` only for source/input compatibility and is ignored by the SSE
geometry.

For an apex distance `R_apex` and angular half width `lambda`, the generating
sphere has center distance and radius

```text
c = R_apex / (1 + sin(lambda))
a = c sin(lambda).
```

A heliocentric ray separated from the CME axis by `alpha` intersects the
outward front at

```text
R(alpha) = c cos(alpha) + sqrt(a^2 - c^2 sin(alpha)^2),
```

provided `alpha <= lambda`.  At `alpha=lambda` the ray is tangent to the
generating sphere.  Directions beyond the half width return `exists=false`;
the public geometry API also returns zero radius/normal sentinels so a caller
cannot mistake an absent surface for a physical flank.

The exact outward normal is the normalized level-set gradient

```text
n_hat = (R e_r - c e_CME) / a.
```

All supported shapes are treated as self-similar.  For a fixed surface
direction,

```text
dR/dt = V_apex R/R_apex,
V_sh,n = V_apex (R/R_apex) (e_r dot n_hat).
```

This replaces the former independent `cos(theta)^m` flank-speed factor and
also corrects the ellipsoid, whose flanks previously inherited the full apex
speed.  For a sphere the formula reduces exactly to `V_sh,n=V_apex`; at the
mathematical SSE tangent boundary the normal projection tends to zero.

`Model::shape_radius_normal()` and `Model::diagnose_direction()` now return a
boolean surface-existence flag.  Existing callers that ignore the return value
remain source-compatible, but finite-width-aware code should always test it.
The Cartesian plasma/field evaluators do so internally: outside the SSE angular
support they return the undisturbed ambient solar wind/Parker field and do not
construct a sheath, ejecta, shock normal, or magnetic amplification from a
fabricated flank.

The following deterministic tests validate the corrected geometry:

- `GEO01` — Sun-centered spherical reference: radius is direction independent
  and the outward normal equals the radial unit vector.
- `GEO02` — true SSE apex: the directional radius equals the configured apex
  distance and the normal equals the CME propagation direction.
- `GEO03` — SSE tangent flank: the boundary point at `alpha=lambda` remains a
  finite valid surface point; the ray is tangent and the normal radial
  projection tends to zero.
- `GEO04` — finite-width enforcement: any direction beyond the half width
  returns `exists=false`, zero geometry sentinels, `rc=1`, and `V_sh,n=0` from
  the diagnostic API.
- `GEO05` — SSE level-set residual: returned surface points satisfy
  `|x-c e_CME|=a` to near machine precision over multiple angles/azimuths.
- `GEO06` — analytical normal validation: SSE normals are compared with an
  independent finite-difference gradient of the dimensionless spherical level
  set; the ellipsoid normal is independently checked against its level-set
  gradient.
- `GEO07` — normal-speed validation: the reported `V_sh,n` for sphere, SSE, and
  ellipsoid is compared with the normal projection of centered finite-
  difference surface motion at neighboring times.
- `GEO08` — rotational covariance: rotating the CME axis, solar axis, and query
  direction together must rotate the normal while leaving radius, normal
  speed, and scalar physical compression unchanged.

The analytical/reference calculations in `3d/test_geometry.cpp` intentionally
do not call the production geometry helper.  They are independent checks of
surface radius, level-set membership, normal orientation, and motion.  Test
tolerances are roundoff- or finite-difference-level tolerances appropriate to
each calculation and must not be loosened simply to make a production result
pass.

The geometry tests remain diagnostically focused on shape, normals, finite
angular support, and normal speed.  Shock existence/compression/downstream
physics is covered independently by `SHK01`-`SHK12` below, while the now-
corrected triangulation/topology is qualified independently by `MSH01`-`MSH05`.

## MSH01-MSH05: shock-surface mesh topology, quality, and area sampling

The production mesh is no longer a rectangular theta-phi array.  It stores one
apex, unique periodic rings, and either one finite-SSE boundary ring or one rear
pole for a closed Sphere/Ellipsoid.  The seam is represented by wrapped
connectivity, so there is no second copy of the `phi=0` vertex at `phi=2*pi`.
Triangle metrics reject repeated indices, scale-aware degenerate area, and
inward winding as `INVALID_MESH`; validation must never make a bad mesh pass by
filtering cells after construction.

- `MSH01` — **nondegeneracy**.  Builds minimum, typical, and fine Sphere/SSE
  meshes, including narrow and broad SSE caps.  Every triangle must have finite
  positive area and `compute_triangle_metrics()` must accept the complete mesh.
  The test prints minimum/median/maximum area for each fixture.
- `MSH02` — **outward cell orientation**.  For Sphere, rotated Ellipsoid, and
  rotated SSE meshes, each triangle cross-product normal is compared with the
  analytical surface normal evaluated in the centroid direction.  Every dot
  product must be positive; the worst alignment is reported.
- `MSH03` — **surface-area convergence**.  Refines `nTheta` and `nPhi` by factors
  of two and sums physical triangle area.  The Sphere reference is `4*pi*R^2`.
  For true SSE with apex radius `R` and half width `lambda`, the independent
  translated-sphere reference is `2*pi*R^2*sin(lambda)^2/(1+sin(lambda))`.
  Error must decrease monotonically, the final observed order must exceed 1.5,
  and the fine-grid error must be below 0.5%.
- `MSH04` — **unique apex / periodic seam**.  Inspects connectivity rather than
  only coordinates.  A finite SSE cap must have `1+nTheta*nPhi` vertices,
  `nPhi*(2*nTheta-1)` triangles, no duplicate coordinate pairs, exactly
  `nPhi` boundary edges, manifold edge incidence, apex valence `nPhi`, disk
  Euler characteristic one, and explicit last-to-first ring adjacency.
- `MSH05` — **area-weighted source-patch sampling**.  Verifies the production
  cumulative-area table is strictly increasing and ends exactly at one.  Three
  fixed RNG seeds each draw 200,000 cells through `sample_triangle_by_area()`.
  Aggregated counts are compared with exact physical area probabilities using a
  chi-square statistic and independent incomplete-gamma p-value calculation;
  every deterministic fixture requires `p > 1e-3`.

Run this mesh gate directly with:

```sh
./output/test_swcme --test MSH01
./output/test_swcme --test MSH02
./output/test_swcme --test MSH03
./output/test_swcme --test MSH04
./output/test_swcme --test MSH05
```

The random-number generator intentionally remains outside the production mesh
class.  Production code provides only the area CDF and deterministic mapping
from a caller-supplied uniform variate to a triangle, so AMPS controls seeds and
parallel RNG policy without being able to accidentally revert to uniform-by-cell
sampling.


## DIV01-DIV03: velocity-divergence treatment

These tests qualify the production divergence operators used by SEP adiabatic
energy change.  They enforce the rule that a radial identity may only be used
for a radial velocity field and that the general 3-D operator must demonstrate
its numerical order on an independent manufactured solution.

- `DIV01` — **analytical constant radial wind**.  Configures `SHOCK_ONLY` in
  both 1-D and 3-D, samples several radii and several unrelated Cartesian
  directions, and requires the production result to equal `2*Vsw/r` to
  roundoff.  Directional dependence is a failure.  This test also guarantees
  that the baseline calculation does not reintroduce finite-difference noise.
- `DIV02` — **manufactured nonconstant radial flow**.  Uses
  `Vr=V0(1+a x+b x^2)`, `x=r/r0`, for which the symbolic spherical divergence is
  known exactly.  `swcme::divergence::radial_terms()` must agree at better than
  `1e-10` relative error.  The same test then differentiates representative
  production shock/sheath/LE/TE velocity samples independently and confirms the
  analytical `RadialVelocityState::d_velocity_dr_s_inv` follows the actual
  transport profile.
- `DIV03` — **general 3-D Cartesian divergence**.  A cubic manufactured vector
  field is evaluated at four successively halved steps.  Because centered
  differences are not exact for the cubic terms, the observed error must show
  second-order convergence.  The test additionally forces the production 3-D
  Cartesian operator on an off-axis constant radial wind and verifies
  convergence toward the exact `2*Vsw/r` result.

Run the divergence gates directly with:

```sh
./output/test_swcme --test DIV01
./output/test_swcme --test DIV02
./output/test_swcme --test DIV03
```

A passing implementation must not obtain `FULL_ICME` divergence by projecting
`V` onto `e_r` and differentiating only along the ray.  That formula drops
non-radial Rankine-Hugoniot velocity and angular gradients of a finite shock
surface.  Conversely, `SHOCK_ONLY` should not be made noisier by forcing the
Cartesian finite-difference operator when the exact analytical result is known.



## KIN01-KIN08: shared CME/shock-apex kinematics

`swcme_kinematics.hpp` is the production source of apex radius and speed for
both `swcme1d` and `swcme3d`.  The kinematics tests are classified `COMMON` and
exercise both the common SI solver and the two public dimensional wrappers.
The independent DBM reference in `core/test_kinematics.cpp` evaluates the
closed-form equations directly and does not call the production DBM helper.

The production DBM solves

```text
DeltaV0 = V0 - Vsw
DeltaV(t) = DeltaV0 / (1 + Gamma |DeltaV0| t)
R(t) = R0 + Vsw t
       + sign(DeltaV0) log(1 + Gamma |DeltaV0| t) / Gamma.
```

This sign-aware form is required for both fast and slow CMEs.  `Gamma=0` uses
the exact ballistic solution.  For very small `x=Gamma*|DeltaV0|*t`, the
production code evaluates the logarithmic distance through a short
`log1p(x)/x` series; `KIN04` explicitly straddles that numerical branch
boundary to guard against a discontinuity.

The data-driven mode uses monotone PCHIP radius interpolation.  The input time
table must be strictly increasing and radius must be nondecreasing.  The
interpolant passes through each knot exactly and its derivative is returned as
the apex speed.  Out-of-range queries return `OUTSIDE_TIME` by default;
ballistic endpoint continuation is available only when requested explicitly.

The tests are:

- `KIN01` — fast-CME DBM versus the independent closed form at multiple times,
  including exact 1-D/3-D wrapper equivalence;
- `KIN02` — slow-CME sign-aware branch; verifies acceleration toward `Vsw`
  without clipping or overshoot and checks both wrappers;
- `KIN03` — exact `Gamma=0` ballistic radius/speed in the common, 1-D, and 3-D
  paths;
- `KIN04` — small-`Gamma` analytical accuracy and continuity across the
  series/direct-log numerical branch;
- `KIN05` — long-time stability, monotonic radius, finite state, and monotonic
  approach of speed toward the ambient wind for fast and slow CMEs;
- `KIN06` — data-driven PCHIP exactness at every height-time knot, derivative
  consistency, and 1-D/3-D use of the same trajectory;
- `KIN07` — dense-grid monotonicity, nonnegative propagation speed, and absence
  of cubic overshoot between data knots; and
- `KIN08` — explicit out-of-time status, optional ballistic continuation, and
  rejection of duplicate-time or decreasing-radius tables.

The DBM analytical comparisons use `1e-12` relative accuracy where specified
by the validation plan.  PCHIP knot radii are required to be exact to
`1e-13` relative.  Tolerances must not be relaxed merely to obtain PASS.

The 1-D and 3-D public `Params` structures expose the same kinematic mode plus
`data_time_s`, `data_radius_Rs`, and `data_extrapolation`.  The 1-D convenience
method `SetDataDrivenKinematics()` sets these fields and switches the mode to
`DataDriven`.  The wrappers convert radii from solar radii to SI only when
constructing the common configuration.

A data-driven query outside the measurement interval causes `prepare_step()`
to throw a clear runtime error unless a continuation policy was selected.  This
is intentional: silent cubic extrapolation is a modeling assumption and is
not allowed to masquerade as measured/constrained kinematics.

## SHK01-SHK14: fast-shock existence, Rankine-Hugoniot validation, and shock-surface ownership

The `SHK` group validates the shared production shock solver in
`swcme_shock.hpp` and its integration into both the 1-D and 3-D models.  SHK01-
SHK12 exercise the shared jump physics and are classified `COMMON`; SHK13-
SHK14 are 3-D integration guards for the rule that local shock physics belongs
to the shock surface rather than to an arbitrary Cartesian query point.

The solver first evaluates the upstream fast-mode speed along the shock normal.
A geometric front is not automatically a shock.  The physical condition is

```text
U1n = Vsh,n - V1.n > c_fast,
M_fast = U1n/c_fast > 1.
```

If this condition is not met, the result is `has_shock=false`, compression is
exactly one, and downstream equals upstream.  No sheath compression floor is
allowed to override this decision.

For a physical fast shock, the solver works in a frame moving with the normal
shock speed.  For a trial compression `r=rho2/rho1`, mass conservation fixes
`u2n=u1n/r`.  Tangential electric-field and tangential momentum continuity form
a 2x2 linear system for `B2t` and `u2t`; normal momentum determines `p2`.  A
bracketed solve of total-energy-flux conservation selects the compressive
fast-shock branch.  The trivial `r=1` solution is divided out of the scalar
residual so it cannot be mistaken for the physical shock root.

The individual tests are:

- `SHK01` — fast-shock existence and no-shock threshold.  Includes permanent
  1-D and 3-D regression checks proving that `sheath_comp_floor` cannot create
  a shock when `Vsh=Vsw`.
- `SHK02` — acute `theta_Bn` calculation and invariance under `B -> -B`.
- `SHK03` — strictly parallel limit against the independent gas-dynamic normal
  shock solution.
- `SHK04` — strictly perpendicular ideal-MHD benchmark against an independent
  test-side scalar reduction.
- `SHK05` — oblique benchmark grid and continuity of the selected compression
  branch as shock speed/obliquity vary.
- `SHK06` — mass-flux conservation.  It also samples the production 3-D field
  immediately downstream and verifies that density, velocity, and magnetic
  field approach the exact RH downstream state instead of the ambient wind.
- `SHK07` — continuity of the normal magnetic-field component.
- `SHK08` — tangential ideal-MHD electric-field continuity.
- `SHK09` — vector momentum-flux conservation.
- `SHK10` — total ideal-MHD energy-flux conservation.
- `SHK11` — physical admissibility: positive downstream pressure/density,
  compressive branch, strong-shock bound for gamma=5/3, and entropy increase.
- `SHK12` — near-Mach-one conditioning and smooth approach of compression to
  unity without an empirical floor.
- `SHK13` — shock-state independence from arbitrary query radius.  The test
  calls the legacy scalar wrapper with query radii ranging from well inside to
  far outside the front and requires identical compression, normal speed, and
  `theta_Bn`.  It also proves the fixture is sensitive by confirming that the
  ambient density at those radii differs materially from the density at the
  physical shock surface.
- `SHK14` — canonical-state consistency across `diagnose_direction()` and shock
  mesh nodes.  Every reported radius, normal, compression, and normal speed is
  compared with a fresh `shock_state_direction()` query.  This prevents future
  diagnostic/mesh code from reintroducing a separate shock-strength calculation.

The independent parallel/perpendicular references do not call the production
nonlinear shock solver.  Conservation tests recompute the conserved fluxes from
the returned primitive states.  This prevents a test from passing merely by
reusing the same internal algebra that produced the result.

Typical direct use is:

```sh
./output/test_swcme --test SHK01
./output/test_swcme --test SHK04
./output/test_swcme --test SHK10
./output/test_swcme --test SHK13
./output/test_swcme --test SHK14
./output/test_swcme --all
```

A shock test failure must not be repaired by increasing the compression floor
or loosening conservation tolerances.  Diagnose the frame transformation,
normal direction, upstream state, root bracket, and downstream reconstruction
first.  A super-fast state for which the nonlinear solver cannot identify an
admissible root is reported with `solver_converged=false` and must not be used
for SEP source physics.


## CON01-CON08: observer-shock magnetic connectivity and cobpoint tracking

`CON01`-`CON08` validate the production 3-D connectivity API implemented by
`swcme3d::Model::observer_connectivity()`.  The connectivity solver does not
maintain a second copy of the shock model: candidate points are tested against
`shape_radius_normal()` and final cobpoints obtain their local plasma/shock
state from `shock_state_direction()`.

The observer Parker line is analytical.  For the same Parker field used by the
production field evaluator,

```text
Delta phi(r) = -Omega (r-r_obs) / V_sw.
```

The observer radial direction is rotated around the configured solar axis by
this amount; colatitude remains constant.  The production parameter
`solar_rotation_rate_rad_s` is shared by field evaluation and connectivity, so
there cannot be a hidden Omega mismatch between the two modules.  The default
value remains the SWCME solar-rotation convention; zero rotation is allowed to
exercise the exact radial reference case.

The intersection residual is

```text
h(r) = r - R_shock[u_Parker(r)].
```

The solver scans from a configurable inner radius to the observer.  It refines
sign-changing roots with bisection, searches local minima of `abs(h)` so tangent
roots that merely touch zero are not missed, and explicitly refines transitions
of the finite-SSE surface-existence flag.  All accepted roots must satisfy the
configured surface-residual tolerance.  Roots are retained in increasing
radius; the default selected cobpoint is the outermost root, corresponding to
the first surface encountered while tracing inward from the observer.

Each `ConnectivityRoot` stores:

- cobpoint radius and Cartesian position;
- signed surface residual;
- analytical Parker arc length from cobpoint to observer; and
- the complete production `LocalShockState`, including `has_shock`, normal,
  normal speed, `theta_Bn`, fast Mach number, compression, and full upstream /
  downstream MHD states.

A geometrical connection and a physical fast shock remain separate concepts.
The connectivity state can therefore be geometrically connected while the
embedded `LocalShockState::has_shock` is false; an SEP source must check the
latter before injection.

The tests are:

- `CON01` — zero-solar-rotation radial limit.  A radial observer line is
  intersected with sphere and SSE geometries and compared with exact analytic
  radius/position/path-length references.
- `CON02` — nonzero-rotation Parker line intersecting a Sun-centered sphere.
  The sphere fixes the root radius exactly while an independent Rodrigues
  rotation verifies the longitude/sign of the production Parker mapping.
- `CON03` — no connection to a finite-width SSE shock.  A field line wholly
  outside the configured cap must return `connected=false` and no fabricated
  flank root.
- `CON04` — exact and near-tangent SSE connection.  The exact half-width ray is
  retained as a valid tangent root, a slightly interior ray connects, and a
  slightly exterior ray remains disconnected.
- `CON05` — multiple intersections and deterministic root selection.  A tightly
  wound Parker line through a strongly non-spherical ellipsoid generates
  multiple physical roots; all are returned in radial order and the outermost
  selection is stable under scan refinement.
- `CON06` — time-continuous cobpoint history.  A ballistic finite SSE front
  evolves from disconnected to connected; the history must show one physical
  onset and continuous/monotone cobpoint motion thereafter without one-step
  classification flicker.
- `CON07` — cobpoint-to-`ShockState` consistency.  Shock quantities stored in a
  connectivity root must be identical to a direct production
  `shock_state_direction()` query at that cobpoint direction.
- `CON08` — analytical Parker path length.  The path length carried by a
  cobpoint is compared with an independent closed-form arc-length reference
  and, away from the pole, must exceed simple radial separation.

The default connectivity search tolerance is much tighter than the validation
campaign's `1e-8 AU` root-position target.  Tests deliberately repeat selected
cases with different radial scan densities so correctness cannot depend on a
particular subdivision.  Tolerances must not be relaxed merely to mask missed
roots or unstable tangent classification.

Typical direct use is:

```sh
./output/test_swcme --test CON01
./output/test_swcme --test CON04
./output/test_swcme --test CON06
./output/test_swcme --test CON08
./output/test_swcme --all
```


## 1D3D01-1D3D03: common-core dimensional-equivalence tests

These tests qualify the common-core refactor rather than a new physical model.
They deliberately configure the 3-D model as a Sun-centered sphere propagating
along +X with the solar rotation axis along +Z and compare it with an equatorial
1-D ray.  In this limit the geometry is exactly equivalent, so any difference is
a software duplication/regression rather than a legitimate dimensional effect.

- `1D3D01` compares the canonical `StepState::common` solar-wind cache and apex
  kinematics field-by-field, then evaluates both public model interfaces at
  1 AU.  Density, radial velocity, `Br`, and `Bphi`/Cartesian azimuthal field
  must agree to roundoff; transverse 3-D velocity and meridional field must
  vanish in the chosen symmetry plane.
- `1D3D02` compares the complete shock state at the spherical apex.  The 1-D
  `StepState::shock_jump` is compared against the 3-D
  `shock_state_direction()` result for shock existence, solver status, radius,
  normal speed, `theta_Bn`, fast-mode speed/Mach number, compression, upstream
  density, full upstream/downstream velocity and magnetic-field vectors,
  downstream density, and downstream pressure.

The tolerances are roundoff-level because both interfaces are now expected to
call the same production core.  A failure must be diagnosed as a wrapper input
mismatch, a reintroduced duplicate equation, or a geometry reduction error; the
tolerance must not be relaxed to hide the difference.

Typical use:

```sh
./output/test_swcme --test 1D3D01
./output/test_swcme --test 1D3D02
./output/test_swcme --all
```

`1D3D03` now validates the shared SOURCE acceleration record.  In the exact
spherical/+X reduction it requires 1-D and 3-D source position, normal, shock
speed, compression, `theta_Bn`, Mach number, upstream density/|B|, DSA
phase-space slope, relative area weight, and deterministic serialized record to
be identical before transport is allowed to diverge.

### Shock-state query-location policy

A local shock state is defined by `(time, shock-surface direction)` and not by
the radius of a background-field sample.  The production sequence is:

1. intersect the selected geometry along the requested direction;
2. evaluate upstream Leblanc/Parker plasma at that surface point;
3. compute the self-similar normal shock speed at that surface;
4. solve the shared ideal-MHD jump;
5. let field/mesh/connectivity consumers use that immutable result.

This ordering is physically important because density and magnetic field vary
with heliocentric radius.  Sampling them at an arbitrary query point would make
Mach number and compression depend on the observer's diagnostic location rather
than on the shock.  `SHK13` and `SHK14` are permanent regression tests for this
contract.

## REG01-REG05: sheath/ejecta region model and mode validation

`REG01`-`REG05` qualify the shared `swcme_regions.hpp` contract used by the 1-D
and 3-D field evaluators.  The exact RH discontinuity remains available in the
shock diagnostic API, while transport-facing FULL_ICME fields use the common C1
shock layer selected by RESOLVED_COMPRESSION.  SHOCK_ONLY/SOURCE keeps the
transport background analytical on both sides of the mathematical source.

- `REG01` — **SHOCK_ONLY upstream-field identity**.  Samples points ahead of
  and geometrically behind the expanding front in both 1-D and 3-D.  Density,
  velocity, and Parker magnetic field must be identical to the analytical
  upstream state.  The test deliberately changes FULL_ICME sheath/ejecta
  parameters to extreme valid values and proves they have no effect in
  SHOCK_ONLY mode.
- `REG02` — **FULL_ICME resolved-shock inner RH boundary**.  Evaluates the
  inner endpoint `R_sh-w_sh/2` of the RESOLVED_COMPRESSION C1 shock layer and
  verifies that density, velocity, and magnetic field equal the complete
  production Rankine-Hugoniot downstream state.  This preserves the exact RH
  boundary while allowing the transport field to resolve the jump numerically.
- `REG03` — **magnetic-ejecta density and velocity factors**.  Evaluates the
  middle of the ejecta, away from LE/TE blends, for factors below, equal to,
  and above unity.  `f_ME=0.5` must give `0.5*n_up` and
  `V_ME_factor=0.8` must give `0.8*V_sw`.  Negative factors must fail
  centralized configuration validation rather than be clipped.
- `REG04` — **self-similar local layer nesting**.  Samples a finite SSE cap from
  apex to near-flank and verifies
  `R_LE/R_sh = 1-f_sheath` and
  `R_TE/R_sh = 1-f_sheath-f_ejecta` at every direction.  Layer ordering must
  never invert.  A public configuration whose thickness fractions sum to one
  or more is rejected.
- `REG05` — **continuity/smoothness at artificial region transitions**.  Checks
  the common symmetric LE/TE smoothstep weights and samples both public model
  interfaces around every transition endpoint.  Equivalent spherical +X
  1-D/3-D fixtures must agree to roundoff.  One-sided finite-difference
  derivatives converge across the C1 artificial interfaces.  The resolved
  shock itself is covered separately by `ACC03`/`ACC04` because its smoothing
  is controlled by the acceleration representation rather than by LE/TE region
  phenomenology.

Run them directly with:

```sh
./output/test_swcme --test REG01
./output/test_swcme --test REG02
./output/test_swcme --test REG03
./output/test_swcme --test REG04
./output/test_swcme --test REG05
./output/test_swcme --all
```

The region tests must not be made green by reintroducing a compression floor,
clipping sub-unity ejecta factors, or sorting invalid boundary radii at runtime.
Such behavior defeats the physical/configuration contracts the tests are meant
to protect.


## ACC01-ACC05: single shock-acceleration representation

These tests qualify `swcme_acceleration.hpp` and the shock-smoothing portion of
`swcme_regions.hpp`.  Their purpose is to prevent one SEP population from seeing
a pre-imposed DSA source and a second resolved compression accelerator at the
same shock.

- `ACC01` — **SOURCE explicit source / SHOCK_ONLY flow**.  Requires a physical
  fast shock to produce `source_enabled=true` and
  `resolved_compression_enabled=false`, checks the phase-space DSA slope
  `q=3r/(r-1)`, and samples both sides of the mathematical front to prove the
  transport background remains the exact Parker/Leblanc solar wind.
- `ACC02` — **mutual-exclusion validation**.  Rejects `SOURCE+FULL_ICME`,
  `RESOLVED_COMPRESSION+SHOCK_ONLY`, zero shock width in resolved mode, and a
  negative relative source weight.  These are configuration errors rather than
  runtime conventions.
- `ACC03` — **resolved-compression shock profile**.  Verifies a finite total
  width, exact upstream outer endpoint, exact RH downstream inner endpoint,
  an intermediate midpoint, and zero-slope/C1 matching at both ends of the
  common numerical shock layer.
- `ACC04` — **1-D/3-D resolved-profile identity**.  Samples five normalized
  positions across an equivalent spherical/+X shock and requires density,
  radial velocity, Br/Bx, and Bphi/By to agree to roundoff.  This prevents
  dimension-dependent shock smoothing from masquerading as a transport effect.
- `ACC05` — **resolved mode disables prescribed DSA source**.  Requires
  `source_enabled=false`, `resolved_compression_enabled=true`, no active DSA
  slope, zero source weight, and an explicit `NA` slope in deterministic audit
  serialization.

`1D3D03` complements these tests by comparing the complete SOURCE acceleration
record from equivalent 1-D and 3-D configurations and requiring byte-identical
serialization.

Run the acceleration gates directly with:

```sh
./output/test_swcme --test ACC01
./output/test_swcme --test ACC02
./output/test_swcme --test ACC03
./output/test_swcme --test ACC04
./output/test_swcme --test ACC05
./output/test_swcme --test 1D3D03
```

A passing test must not be obtained by merely masking `div(V)` after constructing
an RH velocity jump in SOURCE mode.  SOURCE removes the resolved shock from the
transport background by using SHOCK_ONLY; RESOLVED_COMPRESSION owns the single
finite-width compression profile.

## SEP01-SEP06: SWCME-to-SEP/AMPS integration contract

These tests qualify `swcme_sep_source.hpp` and `swcme_sep_interface.hpp`.  The
integration layer is a consumer of the production SWCME model, not a second
shock/connectivity implementation.

- `SEP01` — **background adapter identity**.  A 3-D AMPS-facing single-point
  query must reproduce the direct checked production density, velocity,
  magnetic field, `|B|`, and `div(V)` result in SI.
- `SEP02` — **spectrum/unit convention**.  Checks differential-intensity unit
  conversion round trips, relativistic proton rigidity, DSA slope-to-intensity
  indices, unity at `E_ref`, and physical `J(E_ref)` when explicit reference
  normalization is selected.
- `SEP03` — **1-D/3-D complete source-record identity**.  Equivalent spherical
  +X configurations must emit byte-identical deterministic `SEPSourceState`
  records and matching adapter background values.  This is the final interface
  guard beyond lower-level `1D3D03`.
- `SEP04` — **finite-SSE source-surface weighting**.  Builds source records from
  the corrected triangular mesh, requires physical patch areas, normalizes
  active area fractions to unity, and verifies that total relative patch
  weight equals the configured uniform source weight.
- `SEP05` — **cobpoint source reuse**.  Obtains an observer source through the
  production Parker connectivity solver and verifies that source position and
  shock quantities correspond to the selected production cobpoint rather than
  a separately solved surface.
- `SEP06` — **resolved-compression/source exclusion and manifest stability**.
  `RESOLVED_COMPRESSION` must expose no prescribed source spectrum,
  `relative_intensity_shape()` must return `SOURCE_INACTIVE`, and the resolved
  SWCME+SEP manifest must be deterministic.

`SEPSourceState` makes source and connectivity state explicit through
`active`, `connection_evaluated`, and `connected`.  Its geometry and shock
fields use SI units.  Default `RELATIVE_ONLY` normalization remains
dimensionless; a dimensional spectrum exists only when
`REFERENCE_DIFFERENTIAL_INTENSITY` supplies a positive SI `J(E_ref)`.

Run the integration gates directly with:

```sh
./output/test_swcme --test SEP01
./output/test_swcme --test SEP02
./output/test_swcme --test SEP03
./output/test_swcme --test SEP04
./output/test_swcme --test SEP05
./output/test_swcme --test SEP06
```

## ERR01-ERR05: explicit numerical-status propagation

These tests enforce the rule that the physics layer may not convert a failed
numerical state into a plausible physical value.  They qualify the shared
`swcme_status.hpp` status contract, the checked 1-D/3-D evaluators, the explicit
RH solver outcome, and output-data validation.

- `ERR01` — **1-D outside-domain query**.  Requests a radius below the supported
  `1.05 R_sun` Parker/Leblanc boundary and requires
  `OUTSIDE_MODEL_DOMAIN` with `sample_index=0`.  Sentinel outputs at the failed
  sample must remain unchanged.  The source-compatible void evaluator must
  throw rather than clipping the radius to the domain boundary.
- `ERR02` — **3-D non-finite Cartesian input**.  Inserts a NaN into the second
  point of a two-sample batch and requires `NONFINITE_INPUT` with
  `sample_index=1`.  The failed sample must not be rewritten as zero/ambient
  plasma, and the legacy wrapper must throw.
- `ERR03` — **Rankine-Hugoniot outcome classification**.  A degenerate normal
  and a non-finite upstream magnetic component must be `INVALID_INPUT`; a
  well-formed sub-fast front must be `NO_SHOCK`; a well-conditioned fast shock
  must be `SOLVED`.  The test prevents a bad primitive state from masquerading
  as an unmagnetized/no-shock solution.
- `ERR04` — **no arbitrary +X normalization fallback**.  A zero shock direction
  must return `DEGENERATE_VECTOR`, leave `surface_exists=false`, and make the
  old bool shock wrapper throw rather than constructing a +X surface.
- `ERR05` — **writer rejects corrupt physics**.  A shock mesh containing NaN is
  passed to the checked surface writer.  The writer must return
  `NONFINITE_RESULT` before creating the requested file; substituting zero for
  the bad coordinate is forbidden.

Run the status gates directly with:

```sh
./output/test_swcme --test ERR01
./output/test_swcme --test ERR02
./output/test_swcme --test ERR03
./output/test_swcme --test ERR04
./output/test_swcme --test ERR05
./output/test_swcme --all
```

`NO_SURFACE` and RH `NO_SHOCK` are deliberately not treated as numerical
failures: they represent valid physical/geometrical outcomes.  The SEP layer
similarly treats `NO_CONNECTION` and `SOURCE_INACTIVE` as explicit expected
absence states.  Conversely, a solver, normalization, domain, or
non-finite-value failure must remain visible to the caller and must never be
made green by restoring `finite_or`, radius clipping, or ambient/zero
substitution.
