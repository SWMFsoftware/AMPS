# Explicit physical sodium case preparation

This producer prepares input for the existing experimental native sodium owner.
Preparation is software/source preflight; it does not execute particles, admit a
model error budget, change source permissions, or qualify I21/M3.

## Contract recorded before implementation

Purpose and entry points: `prepare_physical_case.prepare()` authenticates a real
state profile using the maintained `profile_projection.project()` and writes the
exact scalar/block grammar consumed by `sodium_cohort/native.cpp::Parse`.
The production profile, composer, application parser and native owner remain
authoritative at execution. This helper does not reproduce their physics.

Inputs: an explicit profile JSON, its eight-column energy crosswalk CSV, an
application `.in` file, and a case-control JSON. Every native scientific control
must be supplied; no tolerances, particle counts, seeds, rates, state corrections,
timesteps, model budgets, or reference distances are selected here. Numbers are
decimal strings (unsigned integers may also be JSON integers). Native spatial
coordinates/velocities are LSO metres and metres/second; intervals are physical
seconds. Model/ledger bounds retain the native SI units. Geometry epoch and
physical clock are supplied by normal startup, never invented during preparation.

Oracle and procedure: compare emitted assignments to explicit inputs and the
native owner schema; independently verify SHA-256, CSV state/source bindings,
exact binary64 step/horizon equality and ordered profile projection. Administrative
fixtures test only parser/provenance/transaction semantics, never observations.
Tests cover valid complete preparation, missing/null/unknown controls, corrupt
hashes/crosswalks, original-source mutation, path and line injection, output
non-overwrite, include graphs and CLI status. No acceptance threshold is changed.

Artifacts: `native.in`, byte-preserved `profile.native.in`, `case.json`,
`run_manifest.json`, and `result.json` in a previously absent output directory.
The manifest hashes original inputs, complete application include documents,
ordered source files, the generated files and preparer dependencies. Originals
remain in place. A refusal publishes no directory or runnable input.

Status: PASS means this preparation contract ran successfully; SKIPPED/77 means
required explicit configuration/state is absent; ERROR/2 means malformed input,
source corruption, an unsafe path, or infrastructure failure. There is no
observational agreement or formal predecessor claim. A later native refusal is
preserved as a separate result.

## Commands and mandatory controls

Write a new, deliberately incomplete template:

```sh
/usr/bin/python3 -B srcMoon/test/sodium_cohort/prepare_physical_case.py \
  --write-template /tmp/my-explicit-sodium-controls.json
```

This command writes all **64** mandatory control names with `null` values and
returns 77, because the template is not runnable. Fill every value from the
reviewed case specification. Decimal physical/numerical values must be JSON
**strings**; JSON floating-point numbers are refused because their original
precision cannot be recovered. Counts/indices/seeds accept canonical uint64
strings or integer JSON tokens. The top-level schema is exactly:

```json
{
  "schema": "srcMoon-native-sodium-case-controls-v1",
  "cohort_id": "an-explicit-cohort-identity",
  "controls": {
    "physical_step_s": null
  }
}
```

The excerpt above is intentionally incomplete. `--write-template` is the full
list. Its controls comprise all three state-specific inertial/acceleration/
relativistic error limits; mass/momentum/energy/threshold/rotation/number/rate
and planner absolute/relative bounds; explicit root/exposure limits; immutable
seed/genealogy capacity; compiled indices, particle/step counts and represented
population; LSO position/velocity; solar endpoint, frame, origin and source
reference distance/identity; and hash-bound transport-policy evidence. There
is no assumed energy correction, source rate or physical model budget.

Once every control is explicitly supplied:

```sh
/usr/bin/python3 -B srcMoon/test/sodium_cohort/prepare_physical_case.py \
  --profile /absolute/path/to/profile.json \
  --crosswalk /absolute/path/to/na_source_energy_crosswalk.csv \
  --application-input /absolute/path/to/application.in \
  --case-json /tmp/my-explicit-sodium-controls.json \
  --output-dir /absolute/path/to/a-new-prepared-case
```

Relative policy-evidence paths belong to the case JSON, profile file paths to
the profile JSON, and application includes to their assigning input document.
Every referenced file remains hash-bound at its original location. The output
parent directory must exist; the case directory must not exist. No input file,
source packet, application configuration, generated build, or database lock is
modified. Unrepresentable newline/comment/quote injection is refused.

The application must explicitly select either new force-cut-surface order,
`tracing_clock_mode = common-time`, `simulation_environment_mode = frozen`, and
`m2_environment = off` or `disabled`. Homogeneous cases additionally require
explicitly zero source velocity. Step count times the **binary64** step must
equal the binary64 horizon exactly; the preparer never truncates a final step.
The native parser and linked physics admission checks still run afterward.

`reference_weight_mode` must be explicitly `production` or `explicit`.
Production mode forbids extra weight assignments. Explicit mode additionally
requires `parent_reference_weight` and `ion_reference_weight` as positive decimal
strings representable in binary64. These are normal block reference weights,
not represented atom counts; the native owner assigns them before cohort
creation. The preparer preserves any declared unequal ratio and never derives
either reference weight from the represented number or the other species.

The eight-column CSV must be the profile's pinned atomic crosswalk. All fields
must agree with its pinned metadata; the selected three states must be available
and explicitly describe the same 23Na nucleus, neutral/+1/-1 charges, and common
zero ground reference. The ion's exact formation source identity is checked.
The archived threshold/rate tables are authenticated and the source threshold is
compared in the same binary64 representation as the production registry. This
neither replaces the retained D04 threshold with the physical ionization energy
nor assigns a rate-bias bound to that distinction. Source-effective photon
kinematics remain the explicitly named experimental convention.

Run administrative verification independently of an application build:

```sh
/usr/bin/python3 -B srcMoon/test/sodium_cohort/test_prepare_physical_case.py
```

These controls use labelled fictional parser/provenance records only. A source
authentication check with the actual retained profile is a separate scoped
preflight; neither kind of check is a physical cohort or observational test.

One control independently compiles the maintained `amps_input.h` parser and
checks the emitted raw `Assignment.value` values. The native owner's `Text`
reader does not remove quotes: state IDs, solar endpoint, photon policy and the
full energy-reference token (including its inline `#` record) are emitted bare.
File paths alone are double-quoted for the owner's `Path` reader. This verifies
the actual native grammar without configuring AMPS or publishing particles.
