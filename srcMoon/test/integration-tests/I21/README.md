# I21 — Sodium photoionization survival and Na-to-Na+ conversion

The complete specification is [reference/acceptance.json](reference/acceptance.json).
Original acceptance: Ensemble lifetime relative error <=2%; inverse-square rate error <=0.1%; Na-to-Na+ species budget closes under the reviewed physical convention.

The [approved2026-11-18 addendum](../../M3_QUALIFICATION_ADDENDUM_20261118.md)
allocates relative0.002 numerical bias and fixes endpoint-survival lifetime
as the homogeneous primary. The owner snapshots that approved protocol;
actual source/state providers, native execution and formal predecessors remain
required. Numerical/statistical approval is no longer a missing permission.

This owner checks formal capabilities and predecessors. It also accepts an
explicit homogeneous production cohort under the [execution contract](EXECUTION_CONTRACT.md).
The cohort uses the existing source-bound Na/Na⁺ producer and independent
timing, inverse-square rate and ledger audits. Six launches cover matched
8/16/32 steps at 1 and 8 MPI ranks with one OpenMP thread per rank. The full
distance, spectral, shadow, transport and weight families remain required
for I21 qualification.

From the repository root (an existing dedicated campaign directory):

```sh
python3 -B srcMoon/test/integration-tests/I21/test.py \
  --amps "$PWD/amps" --input srcMoon/examples/m3_reaction_rates.in \
  --mpi-ranks 8 --campaign-mode regression \
  --output-dir test_output/srcMoon/M3-formal-NEW
```

The runner creates `<output>/I21/` with its real result/manifest,
stdout/stderr and an exact input snapshot when available. The default command
evaluates prerequisites only. Exit codes follow the parent result:
`0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`.

The owner now retains `protocol_preparation.json`, binding the concrete
[physical model](../../M3_PHYSICAL_MODEL_SPECIFICATION.md) and
[qualification protocol](../../M3_QUALIFICATION_PROTOCOL.md). Missing inputs
remain SKIPPED and corrupt identities are ERROR.

### Run the pinned development cohort

For the fixed-point lunar shadow case, select a sphere, a frozen environment,
zero initial velocity, the reviewed zero-force producer and Earth occultation
off. The auditor independently projects the lunar centre onto the finite
particle-to-Sun segment using LSO coordinates and the lunar radius from the
authenticated producer. Each `lunar_shadow-XX.json` records the closest
distance in metres, radius, intersection class, native visibility and source
identities. Every original parent must survive through the full common
horizon. The six launch-specific event ledgers must be empty. A fully dark
fixture tests the zero-hazard invariant only. Positive-rate lifetimes,
branching precision and the M1/M2 predecessors have their own gates.

First prepare the source-bound profile, application, case controls and
isolated zero-acceleration producer as in [the commissioning workflow](../../../RUN.md).
Retain that producer checkout and its generated ABI. Supply the approved
reference input receipt and the SHA-256 from its original owner manifest. For
the commissioning case this receipt provides the previously declared geometry
and visibility. The preparer derives the 1-AU-scaled expectation independently
from D04 and reads no event or survival counts. `--branch-claims 6` records the
six selected development launches. The full qualification family has its own
predeclared claim count and power design.

```sh
/usr/bin/python3 -B srcMoon/test/integration-tests/I21/prepare_cohort.py \
  --profile "$NA_PROFILE" --crosswalk "$NA_CROSSWALK" \
  --application-input "$NA_APPLICATION" --controls "$NA_CONTROLS" \
  --reference "$NA_REFERENCE_RECEIPT" --reference-sha256 "$NA_REFERENCE_SHA256" \
  --producer "$NA_PRODUCER_MANIFEST" --branch-claims 6 \
  --output-dir "$NEW_PLAN_DIR"

PLAN="$NEW_PLAN_DIR/cohort_plan.json"
PLAN_SHA256=$(sha256sum "$PLAN" | cut -d ' ' -f 1)
/usr/bin/python3 -B srcMoon/test/run_all_test.py --only I21 --skip-build \
  --campaign-mode regression --m3-input "$NA_APPLICATION" \
  --i21-cohort-plan "$PLAN" --i21-cohort-plan-sha256 "$PLAN_SHA256" \
  --test-runtime-seconds 1800 --output-dir "$NEW_RUN_DIR"
```

`--test-id I22` (default `I21`) freezes the same matrix for the I22 atomic
photoionization pilot (approved 2026-11-26). The parent must be He, Ne or Ar,
taken from the authenticated reference receipt, on the illuminated zero-force
x-axis line with the SPICE Sun. The plan then runs through
`run_all_test.py --only I22 --i22-cohort-plan PLAN --i22-cohort-plan-sha256 SHA`.
See [the I22 README](../I22/README.md).

The runner authenticates the reference geometry, D04, the native input graph,
and the producer's recorded source checkout. The configured wall allowance
covers prerequisite handling, all launches and the independent audit. The
`I21/homogeneous/` directory contains per-launch native records, the launch
journal, per-check status, and `independent-audit/numerical_comparison.*` in
CSV, PNG and EPS. The development path reports endpoint survival as diagnostic
and retains the numerical gate verbatim. An evaluated supporting FAIL or
execution ERROR is surfaced by the aggregate. Formal I21 remains SKIPPED.
With `--campaign-mode qualification`, the candidate cohort additionally
requires the same-campaign I13–I20 and U14 PASS receipts before launch and
applies the frozen statistics scorer. One candidate cohort leaves the overall
I21 gate open.

## Varying-hazard reference

Run the standalone prescribed radial case with
`/usr/bin/python3 -B srcMoon/test/test_m3_trajectory_exposure.py --output-dir test_output/srcMoon/i21-radial-NEW`.
The runner writes `radial_refinement.csv`, PNG/EPS plots and the numerical
check records. See the [varying-exposure contract](VARYING_EXPOSURE_CONTRACT.md)
for the closed-form clock and fixed controls.

To run the linked straight-line case, use the same frozen preparation and
aggregate commands above. Supply a dedicated controls JSON declaring
`case_kind = trajectory-integrated-hazard`, explicit x-only LSO velocity,
zero-force source-bound producer, and the planned exposure tolerance.
The plan builder uses the same authenticated fixed-position receipt for the
initial geometry and preserves the moving launch inputs at all six levels.
Use `--campaign-mode regression` for this numerical development case.
The output under `I21/trajectory/` retains actual native receipts, the event
time/position comparison, the moving-path survival diagnostic and the
independent inverse-square and budget checks. The constant cohort's result
remains in its original campaign.

### Select the initial optical-clock precision

For the illuminated fixed-position refinement, add
`"optical_clock_precision": "compensated-extended-log"` to the `controls`
object in your copy of the case JSON **before preparing the plan**. If you
prepare a single native input directly, the equivalent assignment is
`optical_clock_precision = compensated-extended-log` inside its existing
`sodium_qualification` block. The six-launch preparer writes the assignment
and freezes the selection. Each `native-output/input_receipt.json` reports
the mode used by the compiled binary. Existing input files that omit the
selection use `binary64-log`.

The production controller computes the first-hit optical depth from the same
open uniform for either choice. Extended mode evaluates the logarithm before
splitting its result into the existing two binary64 words. A continuing parent
retains both words during movement and MPI migration. Each new daughter gets
its own keyed draw under the selected convention. The mode is part of the
opaque particle configuration fingerprint. The experimental extended mode
requires a binary `long double` with 54–106 significand bits.

Inspect the timing maxima in `independent-audit/numerical.json`, the
per-parent differences in `numerical_comparison.csv` and the matching
PNG/EPS figures. The independent comparison uses exact native event times
and an 80-digit logarithm of each original parent's keyed uniform. Timing,
rate, chronology and signed ledger criteria are scored individually.
Production-time error contributions and the preserved historical comparison
are documented in the [numerical audit guide](../../sodium_cohort/NUMERICAL_README.md).

### Reproduce the lunar shadow fixture

The checked-in [sphere input](reference/shadow_sphere.in) and
[cohort controls](reference/shadow_controls.json) specify the compact mesh,
quiet solar endpoint, fixed spherical occluder, seed, source point, 65,536 s
horizon and 128 original parents. Obtain the Na/Na⁺ producer, ground-state
profile and crosswalk using the source-qualification build instructions
above. Generate the source-bound input and an independent native geometry
receipt first. From the repository root run

```bash
case_root=srcMoon/test/integration-tests/I21/reference
NA_APPLICATION="$PWD/$case_root/shadow_sphere.in"
NA_CONTROLS="$PWD/$case_root/shadow_controls.json"
work=$(mktemp -d)
export NA_APPLICATION NA_CONTROLS
/usr/bin/python3 -B srcMoon/test/sodium_cohort/prepare_physical_case.py \
  --profile "$NA_PROFILE" --crosswalk "$NA_CROSSWALK" \
  --application-input "$NA_APPLICATION" --case-json "$NA_CONTROLS" \
  --output-dir "$work/reference-input"
producer_root=$(/usr/bin/python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["source_root"])' "$NA_PRODUCER_MANIFEST")
/usr/bin/python3 -B "$producer_root/srcMoon/test/run_m3_sodium_cohort.py" \
  --mode run --mpi-ranks 1 --producer-manifest "$NA_PRODUCER_MANIFEST" \
  --input "$work/reference-input/native.in" --output-dir "$work/reference"
NA_REFERENCE_RECEIPT="$work/reference/native-output/input_receipt.json"
NA_REFERENCE_SHA256=$(sha256sum "$NA_REFERENCE_RECEIPT" | cut -d ' ' -f 1)
NEW_PLAN_DIR="$work/plan"
NEW_RUN_DIR="$work/execution"
```

Use these paths in the plan-preparation and aggregate commands above. The
geometry-only comparison uses the Sun direction actually returned by SPICE.
The reference-segment test also checks the Moon's radius from the authenticated
build. No observation files are needed for this analytical shadow case.

## Decision-1 cases (approved 2026-11-26)

The user approved decision 1 of the M3 completion follow-up and its case
design on 2026-11-26. The [decision-1 contract](DECISION1_CONTRACT.md) records
the rules (the constant-cohort FAIL stays visible; the existing varying-rate
pair satisfies item 2), the cases, the production changes and the gates.
Formal I21 remains SKIPPED.

| File | Role |
| --- | --- |
| `reference/decision1_approval_20261126.json` | The approvals and case parameters. |
| `prepare_decision1_cases.py` | Generates `reference/decision1/` (case controls and the straight, active and gravity decks for both movers) from `reference/shadow_sphere.in` and `reference/shadow_controls.json`; `--check` requires byte equality. |
| `compare_distance_endpoints.py` | Scores the native 1 AU / 2 AU rate ratio (4 within 0.1%) and the active/quiet D04 ratio across the four prescribed-Sun cases. |
| `../../test_m3_i21_decision1.py` | U16 controls: the shadow and Kepler oracles against direct integration and an RK4 orbit, the generator and the gravity recipe. |

Each case uses the documented chain: `sodium_cohort/prepare_physical_case.py`
with the case JSON and deck, a one-rank reference run with
`run_m3_sodium_cohort.py --mode run`, `prepare_cohort.py` for each mover
(`--trajectory-oracle moon-shadow-crossing-zero-force-straight-line-v1` for
`moon-shadow`; the Kepler oracle is selected automatically for
`lunar-gravity`), then `run_all_test.py --only I21 --i21-cohort-plan ...` as in
the commands above. Use a Na/Na⁺ producer built with
`--reference-problem homogeneous-zero-acceleration` for the straight-line and
shadow cases, and with `--reference-problem lunar-point-mass-gravity` for the
gravity case ([builder](../../source_qualification_build/README.md)). Select
the `_active` decks for the active endpoint; the normal runtime refuses a deck
whose solar endpoint differs from the cohort.

The independent audit writes `independent-audit/numerical.json`. For the shadow
case its `trajectory_oracle` block records the entry and exit times and the
depth at entry and at the horizon. For the gravity case it records the
production lunar GM, semi-major axis and period, and `position_refinement`
lists the maximum event-position error per level, which must strictly
decrease.
