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
