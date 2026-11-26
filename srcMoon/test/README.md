# `srcMoon` subsystem verification tests

## November 19 current-core preservation and supervision

Master integrates `60a156ddcc3deeaf48e2981b0beb14738245457e` with the preserved local source (eleven incoming commits). The completed corrective CORE registry has **844 PASS, 240 SKIP, 1 FAIL and 1 ERROR**. Exact historical/current-upstream PASS preservation audits pass; the known nonpasses remain visible. Moon M0 is **10/10 PASS**, with a zero-byte reviewed nightly diff, and **1533 targeted native checks PASS**. The U-series has 15 PASS/11 SKIPPED and zero FAIL/ERROR.
See [CORE_INTEGRATION.md](../CORE_INTEGRATION.md) for exact source bases, the
retained original infrastructure ERROR, native manifests and complete report.
A preservation PASS does not make the registry all green or qualify M3.

`CampaignBudget` observes descendants of its calling process. Do not run two
independent native budgets in threads sharing that PID: either budget can
adopt/wait for/terminate the other's child. Use one native owner per OS process,
or separate OS worker processes with independent work/configuration trees.
Include setup/build/run/cleanup/reporting in each original deadline; preserve
failed attempts and never infer individual PASS from incomplete suite logs.
The corrected campaign changes only this orchestration and build parallelism,
using the same registered tests, cases, scientific criteria and authentic data.

## M3 completion proposal — November 17

The [November 24 follow-up](M3_COMPLETION_FOLLOWUP_20261124.md) sets out the
recommended order, physical model choices and qualification packages after
the first linked I21 sodium matrix. It records the remaining decisions and
source evidence for I22–I24 and the same-campaign M1/M2 predecessors.

The [I21 trajectory case](integration-tests/I21/VARYING_EXPOSURE_CONTRACT.md)
uses the existing I21 owner. The [I23 helium pilot](integration-tests/I23/README.md)
is opt-in and retains its supporting rate and event-admission results.
The global extended adapter selects 40 U/I cases including I21–I24. The
nightly `app-moon` family continues to select U01–U26.

[M3_COMPLETION_PROPOSAL.md](M3_COMPLETION_PROPOSAL.md) specifies the recommended
remaining R3.1–R3.4 work: conservative reaction packets, actual force-aware mover
wiring, photon/electron/charge-exchange/molecular contracts, existing I21–I24
runners, independent oracles, statistical design and separate one-hour/cluster
tiers. It records the available D04/D05 data, missing kinetic/driver evidence
and decisions requiring review. The user approved implementation and the
additional first-/second-order local particle tracers.
[M3_IMPLEMENTATION_CONTRACT.md](M3_IMPLEMENTATION_CONTRACT.md) records that scope.
Production services now have local/native controls; unresolved physical
providers and formal predecessor gates still keep M3 open. Dated initial
rate-only results below retain their original scope.

## M3 registered support and force tracing — November 17

The existing U15/U17 owners now invoke the maintained photon statistics,
reaction primitives, runtime configuration and total-capture rate controls.
Their actual `supporting_checks` records, counts, runner and log paths survive
the global aggregation. A component PASS is separate from the full-ID SKIPPED
status: molecular kinetics, complete source/state mapping and unchanged
predecessors are still required. Child FAIL/ERROR is never hidden by SKIPPED.

The extended `--phase m3-development --isolated-checkout` also executes the
actual PIC packet adapter on one/eight ranks, both tracing orders with the
actual Moon force and the combined constant-history controller, using one
fresh all-species producer. These fixtures test
unequal local weights/time steps, all-or-nothing publication, allocation
failure, local time and convergence. Their analytical elastic/constant-field
cases are not lunar measurements or qualified CEX/daughter laws. Without
`--isolated-checkout`, the derived local-weight/time build is explicitly
SKIPPED. Commands and artifact descriptions are in [RUN.md](../RUN.md).

The generic core tracer has its own normal nightly registration
`unit_test_force_cut_surface_mover`; it does not borrow a lunar U/I ID. See
[its README](../../test/srcForceCutSurfaceMoverTest/README.md) for clean build,
one/eight-rank fixed/adaptive runs, negative mode controls and PNG/CSV/EPS.
The six-body reflection regression additionally covers both new orders and
the legacy default. A native core PASS does not qualify Moon environmental
time dependence or missing reaction physics.

I21–I24 now have individual registered owners and detailed contracts. They
currently implement predecessor/capability gating and evidence preservation;
their full physical execution and scoring remain unavailable. They produce
honest SKIPPED / NOT QUALIFIED receipts instead of placeholder success.
No formal M3 gate or acceptance criterion changes through this registration.

Supporting results must declare `supporting_checks` together with
`supporting_counts`. Every child has a unique nonempty `check`, a valid
`status`, and explicit `runner` and `log` paths. Counts must be exactly the
four integer fields `PASS`, `FAIL`, `ERROR`, `SKIPPED`, matching those rows.
An outer PASS cannot hide a non-PASS child; a formally SKIPPED owner may retain
successful/unavailable components but cannot hide a child FAIL/ERROR. An
inconsistent receipt becomes ERROR and preserves the original child diagnostics.
Child CLI exits remain0/1/2/77 for PASS/FAIL/ERROR/SKIPPED. The shared supporting
deadline covers descendants and later groups; timeouts/resource failure are
ERROR, and no later group receives another full allowance. The dedicated
capacity-refusal control accepts only its intended exit1 plus its specific
actionable diagnostic, never a signal/crash carrying matching text. The
[protocol contract](M3_SUPPORTING_PROTOCOL_CONTRACT.md) and registered U16
orchestration controls exercise these requirements without science fixtures.

## Initial M3 rate service and research — November 17

The user approved early M3 implementation/local verification while formal
M1/M2 remain open. U16 now compiles the maintained production
`species/ElectronImpact` service, checks all four native neutral atom fits
against independent high-precision values, and runs the actual runtime input
parser/install controls. It remains registered through the existing literal
U16 entry in `run_all_test.py`; no permanent ID was added. A PASS verifies
rates/configuration only, not product trajectories, conservative reactions
or measured plasma forcing. The optional explicit global
`test_runner.py --phase m3-development --mpi-ranks 8` runs the unchanged M0
prerequisite and an actual eight-species callback fixture, then restores normal
Moon configuration in the selected checkout. All builds/children share a
3600-second supervisor budget. Use a new output directory and an isolated
checkout to preserve a separately configured active application.

The complete contract is [M3_INITIAL_IMPLEMENTATION_CONTRACT.md](M3_INITIAL_IMPLEMENTATION_CONTRACT.md).
Charge-exchange and molecular-kinetics research proposals do not authorize
new channels or supply missing daughter distributions. Formal I21-I24 remain
SKIPPED / NOT QUALIFIED. See STATUS.md section 5.15 and RUN.md for evidence,
commands and the distinction between scalar rates and particle reactions.

The final native M3 development campaign completes in 16.3 minutes: U16 PASS,
four formal chemistry SKIPPED, three supporting PASS, no FAIL/ERROR. The
complete M0 baseline passes with zero-byte nightly diff. U16 includes 170 pure
kernel controls, 17 native runtime controls and 31 orchestration controls;
the eight-rank callback separately verifies 32 rows and actual masses/charges.
The all-species fixture exposed a shared He+ `mas`/`mass` input typo; its
intended 4.0-AMU value is now actually installed, and the full 40.8-minute
CORE audit preserves all 843 prior PASS cases (known nonpasses unchanged).
Detailed evidence and the final diagnostic-log-only extension are recorded
in STATUS.md section 5.15 and the work-package completion report. A passing
rate service still does not establish any I21-I24 product/ensemble gate.

## November 16 remaining-M2 implementation

The existing global `--phase m2-development` campaign now includes completion
protocol checks, compiled input-control checks and bounded real I15/I20
exporter fixtures. They retain separate supporting outcomes; formal I14-I20
still require same-campaign I13 and each individual science/data gate. These
fixtures introduce no new permanent U/I IDs and contain no observations.

The final November 16 isolated eight-rank campaign completed in 37.4 minutes:
58 supporting PASS, two SKIPPED, no FAIL/ERROR; the formal registry retains
six local PASS and seven M2 SKIPPED. M0's ten linked gates and zero-byte
nightly diff pass. All 36 short numerical cases pass. A separate 22.64-second
native source-only check validates the requested complete I19 initial-state
export after canonical-build restoration; the global normalizer now requests
that export as well. Its `candidate_initial_state` is explicitly unqualified:
the normalizer's READY is not initial-state approval or formal I19/I20 PASS.
Commands, source/build/data identities, earlier attempts and qualification
limits are in `../STATUS.md` section 5.14 and
`../../test_output/srcMoon/m2-completion-20261116T195511Z/completion_report.json`.
Historical data-package status notes below retain their dates; the current
read-only audit records D01/D03 selected-input readiness and the unresolved
D02 temporal operator, without granting observational validation.

U09 extends its original exact point-ray criteria with independent full/half
disc radiometry, correlated masks and production-BVH penumbra controls. U11
checks the same energy solver's integrated W/m2 boundary and the production
model-depth interpolation kernel. Original acceptance criteria are unchanged.
U10 also checks the approved explicit Diviner grey-facet inverse, production
prescribed-temperature caches, impact callbacks and Na desorption across two
SPICE epochs, including real missing-map behavior. Native I20 support checks
surface Na/Ar residence laws and retention independently of checkpoint writes;
subsurface rows cannot claim surface kinetics. Its scorer binds actual cohort
indices/elements/seeds to frozen sites and reports terrain-only geometric PSR
separately from Earth-policy visibility and thermal retention.

Individual new infrastructure commands are documented in the I15/I16/I17/I19/I20
READMEs and `M2_FORMAL_CONTRACTS.md`; every new check must also remain registered
in the development driver so the global entry cannot silently omit it.
The corresponding individual I-series runners also execute their metadata/
orchestration controls when selected by `run_all_test.py`; a control PASS does
not change a predecessor-blocked scientific SKIPPED. Their failure writes an
ERROR result and remains visible in that runner's aggregate log.

Formal dispatch uses a reviewed `m2_formal_contract_file` in the CLI `.in`.
I19 revalidates all retained native cases and independently re-scores them;
its aggregate binds the real receipt family instead of inventing one native
`receipt.txt`. New long execution requires `--m2-cluster-execute`, an approved
execution contract, measured feasibility, all predecessors and an explicit
wall/allocation budget. The development phase never opts in. Merely selecting
qualification does not authorize a missing contract or unavailable data.

Retained Na/Ar binaries now preserve their actual generated chemistry/layout
and RNG headers alongside their build receipt. Formal execution checks those
immutable snapshots, executable bytes and current native source identity;
an unrelated later Na-only restoration cannot invalidate that evidence.
Development/probe reuse additionally checks the live generated tree, because
new probe objects must use the matching ABI. Missing or changed snapshots
fail closed; this is not permission to reuse a stale executable.

The native source identity also binds root `ampsConfigLib.pm`, the delegated
preprocessor used by `ampsConfig.pl`. Changing its macro-guard processing can
change generated compilation while leaving the launcher and input identical.
The retained-build protocol includes a mutation check for that library in both
strict live-tree and self-contained retained modes. This is provenance
verification, not a physical or numerical qualification gate.

Current-core integration has a separate regression-preservation audit, described
in [CORE_INTEGRATION.md](../CORE_INTEGRATION.md). The existing global core
registry remains the authority for suite/case membership. The helper
`compare_core_regressions.py` reads its aggregate reports and `selection.json`,
protects every actual baseline PASS, and refuses missing cases or a later
failed retry that attempts to borrow an earlier PASS. Repeat `--baseline` and
`--candidate` in explicit attempt/retry order, and choose a new output directory.
It emits `result.json` plus `comparison.csv`, with runner/log locations. Its
three protocol checks execute automatically; protocol fixtures are never
observations. A preservation PASS does not mean the registry is all green or
that a Moon milestone is qualified. Existing FAIL/ERROR/SKIP remain visible.

For relocated Enceladus regression inputs, the optional
`prepare_core_regression_inputs.py` helper runs the existing production exporter
and changes only three absolute host-path comment lines. It requires both
complete non-path payloads to remain identical, backs up their original bytes,
and records hashes/commands. Missing or changed scientific content is ERROR.
It performs no calibration, numerical-reference replacement or observation
processing. The missing native fixture provenance is retained separately in
[`src/test/NATIVE_FIXTURE_PROVENANCE.md`](../../src/test/NATIVE_FIXTURE_PROVENANCE.md).
These helpers support existing registry tests; neither introduces a new U/I ID.

The approved M2 development-order exception has a separate global mode:

```sh
python3 -B srcMoon/test/test_runner.py --phase m2-development --mpi-ranks 8 \
  --output-dir test_output/srcMoon/M2-development-NEW
```

It runs U08-U13 and supporting production I14/I18/I19 checks plus read-only
data QA, preserving formal I14-I20 as I13-blocked SKIPPED. The dedicated clean
Na/Ar build has chemistry OFF and individual weights ON; the canonical
Moon/SPICE configuration is restored afterward. A retained binary may be
supplied through `--m2-amps` only with its matching current generated tree.
See [RUN.md](../RUN.md#m2-supporting-verification-before-m1-closure) and
[M2_IMPLEMENTATION_PLAN.md](M2_IMPLEMENTATION_PLAN.md) for commands, physics,
oracles, resource guards and output paths. Supporting checks have separate
counts and cannot open formal gates. FAIL/ERROR in either namespace makes the
global command fail; an actual SKIPPED remains visible (and fails the command
when `--strict-skips` is requested).

The data-QA child automatically runs ten D02 native-metadata/error protocol
checks through `data_processing/audit_m2_packages.py`. They have no new U/I
identifier and supply infrastructure evidence only. The child preserves
`d02_protocol_checks.log` and `d02_native_temporal_audit.json`; missing exact
averaging metadata remains SKIPPED, and altered frozen bytes are ERROR.
See [D02_ILLUMINATION_RESEARCH.md](D02_ILLUMINATION_RESEARCH.md) for the native
observable and remaining authoritative-data gaps. Spatial QA v2 distinguishes
a completed transformation from the unresolved temporal operator.

The user approved the v3 numerical design on November 13. The runners now
consume [M2_I19_I20_PROTOCOL.json](M2_I19_I20_PROTOCOL.json); the reviewed
proposal remains as design history. See [M2_NUMERICAL_PROTOCOL.md](M2_NUMERICAL_PROTOCOL.md)
for executable commands, physics, analytical oracles, input controls and exact
status semantics. The M2 development phase automatically runs the protocol
infrastructure checks, U12 imposed-warming kernel comparison, original
positive-coverage I19 fixture, 36 short production cases, and deterministic
source-normalization reference. Each grid stage retains CSV/PNG/EPS diagnostics.
All phases, required clean builds and plotting share one monotonic deadline;
new sessions and orphan MPI children remain supervised. The dedicated M2 build
explicitly selects the existing AMPS Mersenne Twister option; M0/M1 RNG defaults
and the reviewed nightly reference remain unchanged.

`--phase m2-cluster` exports a 288-case plan by default. Execution requires an
explicit `--cluster-execute`, same-campaign linked predecessors, READY data,
frozen native sites and physical normalization, justified statistical/RNG
assumptions, and a measured allocation forecast. The default 300 CPU-core-hour
allowance covers the whole selected bundle, never each seed. Missing evidence
is SKIPPED; corrupt evidence or exhausted resource/deadline is ERROR. Short
supporting PASS cannot close formal I19/I20, M1, M2 or observational validation.

New tests must remain registered in the global runner and the frozen test
registry. The new I19 parser/builder unit checks are invoked automatically by
`integration-tests/I19/test.py` and the development driver, under existing I19,
not new permanent IDs. They produce a unittest exit status; the child runner
then writes its standard `result.json`. Every new supporting check must emit
`status` as PASS/FAIL/SKIPPED/ERROR, its `check` name, runner/log paths for
executed failures, and an explicit supporting scope. Never insert supporting
PASS into the formal `results` namespace or fabricate a predecessor result.

For a one-hour M1 launch budget after a clean Moon/SPICE build, the existing
global runner can launch four isolated meshes concurrently:

```sh
python3 -B srcMoon/test/run_all_test.py --skip-build \
  --campaign-mode qualification --only I11 --only I12 --only I13 \
  --m1-input srcMoon/examples/m1_memory200.in --mpi-ranks 8 \
  --i11-parallel-levels 4 --test-runtime-seconds 3600 \
  --output-dir test_output/srcMoon/m1-parallel-NEW
```

This uses 32 hardware threads; `--mpi-ranks` applies to each mesh. I11 first
executes its linked sampling-window and event-budget/Na-chemistry prerequisites
automatically, so new storage
regressions remain registered through its existing global I11 entry. See the
[I11 README](integration-tests/I11/README.md) and
[sampling contract](integration-tests/I11/SAMPLING_CONTRACT.md) for exact
oracles, guards, artifacts and scope. The default serial schedule is retained.

The event prerequisite is specified in
[BUDGET_CONTRACT.md](integration-tests/I11/BUDGET_CONTRACT.md). Its actual
sunlit/shadow survival check protects the production chemistry dispatcher,
which a local lifetime test cannot verify. `test_diagnostics.py` is also
automatically run by I11, and final stages generate an independent per-metric
statistics CSV plus PNG/EPS plots. These are supporting checks under the
existing I11 ID, never extra permanent IDs or observational validation.

For the bounded coarsest-grid diagnostic pilot:

```sh
python3 -B srcMoon/test/integration-tests/I11/run_duration_pilot.py \
  --amps ./amps --input srcMoon/examples/m1_memory200.in --mpi-ranks 8 \
  --wall-seconds 3600 --memory-limit-gib 200 \
  --output-dir test_output/srcMoon/I11-duration-pilot-NEW
```

A cooperative wall stop writes a PARTIAL diagnostic with closing ledger and
receipt; it never reports formal I11 PASS or changes the 5000-step formal cap.
I12/I13 reuse verified finest-mesh initialization receipts when available and
retain separate scores. A timeout produces ERROR for incomplete children;
I13 remains SKIPPED until I11 passes and its candidate is explicitly reviewed.

The 2026-11-11 corrected-sampler campaign on `passion` completed this
four-by-eight-rank schedule in 37.14 minutes at 64.35 GiB peak aggregate RSS,
retaining about 1.28 GB. Its formal result was I11 FAIL, I12 PASS, I13 SKIPPED,
zero ERROR. All I11 stages reached the cap without statistical qualification;
the runtime measurement is operational evidence, not numerical qualification.
See [`STATUS.md` section 5.4](../STATUS.md#54-resumed-verification--2026-11-11)
for exact results, manifests, commands and remaining decisions.

The authoritative campaign procedure is
[`TESTING_PROCEDURE.md`](TESTING_PROCEDURE.md), and the current formal
milestone/evidence ledger is [`../STATUS.md`](../STATUS.md).  The single
top-level entry point is:

```bash
python3 -B srcMoon/test/test_runner.py --phase audit
python3 -B srcMoon/test/test_runner.py --phase implemented \
  --campaign-mode qualification --mpi-ranks 8 \
  --output-dir test_output/srcMoon/qualification-YYYYMMDD
```

The first command checks the U01-U26/I01-I24 runner, README, and reference
contract structure without making a physics claim.  The second delegates the
roadmap-ordered scientific work to the established phase runners. I21-I24
currently execute prerequisite/capability gating and preserve evidence; their
physical execution/scoring is not implemented. I25-I46
remain frozen roadmap IDs without executable owners and are reported as
explicit gaps by `test_runner.py --phase release`; they are never represented
by placeholder PASS results.

The 200 GB allocation check reuses the registered I12 runner with
`--only I12 --m1-input srcMoon/examples/m1_memory200.in --mpi-ranks 8`.
See [the profile contract and direct command](../examples/README.md#200-gb-memory-budget-mesh-2026-10-14).
The approved memory profile can now pass I12 wiring, but M1 still requires
I11 convergence and reviewed I13 baseline evidence. Qualification defaults
to 480/240/120/60 km; explicitly select `--m1-input srcMoon/examples/m1_sphere.in`
only to run the legacy memory-intensive FULL sequence.

For a presentation-oriented but scientifically labelled LOLA product, use
[`lola_surface_demo/run.py`](lola_surface_demo/README.md). It runs the verified
production surface kernels, compares every mesh vertex with an independent
native D01 decode, and writes CSV/PNG/EPS plus CEA/Tecplot surfaces. It is
deliberately outside `run_all_test.py`: its `execution_status=PASS` is kept
separate from formal I14 `execution=NOT_RUN`, `status=SKIPPED`, and
`validation_status=NOT_VALIDATED`; it cannot open any roadmap gate.

For production illumination, topographic-shadow, Earth-occultation, and
depth-resolved thermal model tables plus PNG/EPS figures, use
[`environment_demo/run.py`](environment_demo/README.md). It requires local
U08/U09/U11 PASS and exports every logical surface cell and every subsurface
layer through a thin adapter linked to the production kernels. It remains
outside `run_all_test.py`; its demonstration PASS is recorded separately from
formal I15/I17/I20 `SKIPPED / NOT VALIDATED` status. The default zero-duration
result is the explicit input initial condition, while any transient or
periodic-spinup run is selected visibly on the command line and preserved in
the receipt.

For the bounded Tenishev et al. (2013) Case-1 neutral-Na model, use
[`sodium_exosphere_demo/run.py`](sodium_exosphere_demo/README.md). It executes
the ordinary MPI production source, mover, surface-impact, photochemistry and
sampling paths after U03/U04/U05/U12/U14/U22 pass. It is intentionally outside
`run_all_test.py`: successful artifact generation is `EXPERIMENTAL`, while
convergence and the unavailable qualified 1991 observation comparison remain
separately NOT RUN/SKIPPED.

I17's runner now includes its launcher/plot-reader infrastructure checks
automatically and enforces the approved 2026-10-14 three-step thermal numerical
error budget before residual computation. See the
[I17 contract](integration-tests/I17/README.md) for the 10%-of-observational-sigma
criterion and the `dt_div_4/` location of the finest-step experimental figures.
Infrastructure success is not I17 validation.

This directory implements the local/stand-alone `U01`-`U26` layer and the
currently implemented linked M1 `I11`-`I13` and registered M2 `I14`-`I20`
layers defined by
[`../development_and_validation.md`](../development_and_validation.md).  It does
not replace the linked `I01`-`I46` campaign.  A local-kernel pass is evidence
for that kernel only; it is not evidence that the callback is reached in a
production run, that the numerical solution converges, or that observations
agree with the model.

For milestone M0, every literal U01-U26 ID remains registered and is executed
with the four-state semantics below.  Local gates used directly by the linked
I01-I10 baseline must pass; a U-test for production physics scheduled in a
later milestone remains honestly `SKIPPED` until that capability exists.  M0
does not convert such a skip into a pass and does not pull M2-M5 physics into
Phase 0.  All linked I01-I10 tests themselves must pass before M0 is complete.

The normal M0 science configuration remains neutral Na only.  I08 will use a
separate committed `Na,NA_PLUS` configuration solely to exercise the real ion
mover in analytic uniform fields; this does not advertise Na+ as enabled in
the normal production input.

The production input now explicitly selects
`PIC::Mover::UniformWeight_UniformTimeStep_SecondOrder`. This is required
because AMPS's global default is the no-force trajectory/boundary mover; merely
installing `Moon::TotalParticleAcceleration` does not cause that default mover
to call it. I05-I10 therefore test the same force-aware mover selected for
ordinary Moon particles, rather than a test-only integration path.

## Layout

Each frozen test ID has its own directory:

```text
UXX_descriptive_name/
  README.md
  test.py
  reference/acceptance.json
```

`acceptance.json` is the machine-readable test contract and reference
solution.  It records the purpose, production entry points, configuration,
inputs, independent oracle, units, frame/time convention, stochastic controls,
procedure, metrics, thresholds, expected result and failure modes, required
data and hashes, predecessor gates, artifacts, and status semantics.

The common C++ probe links against the generated production library
`build/libAMPS.a`.  It calls compiled `srcMoon`/exosphere functions directly;
it does not copy the production formulas into a test-only implementation.  Its
expected values come from analytical invariants or explicit control points.
Source guards are deliberately labelled as such and never reported as linked
runtime verification.

### Runner and probe design

`common/moon_testlib.py` is orchestration only. It discovers frozen IDs, reads
each `acceptance.json`, compiles the thin C++ probe when required, preserves
stdout/stderr and the compile command, assigns the four-state result, and
writes JSON artifacts. It does not contain lunar physics.

`common/production_kernel_probe.cpp` links the generated `libAMPS.a`, so calls
resolve to the same configured production functions used by `srcMoon`. The
probe contains only independent expected values and invariant calculations:
closed-form point gravity, inverse-square ratios, independently sourced table
control points, and topology/orientation checks. A probe PASS means only that
named kernel and fixture passed; it is not a full executable, convergence,
conservation, or validation result.

U05 additionally verifies the absolute one-AU Na radiation-pressure curve
against a reproducible two-resolution digitization of Combi et al. (1997),
Figure 7. The historical digitization that produced the embedded array is
unavailable and is not claimed as provenance. The runner verifies every
reference SHA-256 before comparing 14 unobscured curve points with the compiled
production kernel using a predeclared `0.5 cm s^-2` graphical uncertainty.

U08 has an additional two-part status:

1. the local production parser/decoder/geometry/writer probe runs whenever the
   exact raw IMG/LBL pair is available;
2. the enclosing U08 result checks D01 package qualification independently.

Consequently the kernel can report PASS while U08 reports
`SKIPPED / NOT VALIDATED`. This is deliberate and prevents raw directory
presence from being treated as qualified validation evidence.

## Running

### Repository-wide procedure

The global routine owns stable identities for every implemented test in this
directory. U01-U26 and I11-I20 are imported through an exclusive clean
Moon/SPICE regression build; I01-I10 are separately visible as explicit SKIPs
until the historical multi-build M0 orchestrator is converted to the current
mandatory native-build preflight.

List both global owners without running them:

```sh
cd /home/vtenishe/Mars2/AMPS
python3 src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --component srcMoon --list
```

To execute the regression owner, provide site paths in an overlay and enable
the optional capability:

```json
{
  "schema_version": 1,
  "paths": {
    "spice_path": "/site/SPICE/cspice",
    "spice_kernels": "/site/SPICE/Kernels"
  },
  "capabilities": {
    "moon-spice": {"available": true, "policy": "optional"}
  }
}
```

```sh
python3 src/test/run_tests.py \
  --test APPLICATION-MOON-REGRESSION-001 \
  --config /path/to/moon-machine.json
```

The repository-wide full selector requests the local qualification mode rather
than regression mode and selects every other registered core/application
owner as well:

```sh
python3 src/test/run_tests.py --full-validation --jobs 8 \
  --config test/qualification-machine/passion.json
```

This explicit request does not supply SPICE paths or missing M2/observational
packages. Such prerequisites and the still-preflight-blocked I01-I10 owner
remain visible as typed failures or skips in the combined report.

The global adapter performs the exact root/process check, removes only root
`build`, configures Moon/SPICE, and runs `MAKEFLAGS=-j16 make amps` before
calling this directory's aggregate with `--skip-build`. Local SKIPPED results
remain global SKIP; they are never upgraded from a zero convenience exit.

### Complete campaign

`run_all_test.py` is the phase engine for all 36 currently registered
`srcMoon` U-series and linked tests: U01-U26 followed by I11-I20. The global
entry point is `test_runner.py`, which also orders the separate I01-I10 M0
engine and exposes the later roadmap gaps. By default `run_all_test.py` runs
the bounded `regression` campaign with eight MPI ranks: I11 and I13
produce explicit `SKIPPED` non-execution records, while I12 runs the actual
production initialization and wiring with `srcMoon/examples/m1_regression.in`.
The default regression campaign builds only the canonical Na executable,
because I13 and hence I18-I20 are deliberately not launchable. A qualification
selection containing both I13 and any of I18-I20 first creates and preserves a
dedicated Na+Ar M2 executable and its hashed build receipt. It then removes the
generated `build/` tree again, restores the canonical Na-only Moon/SPICE
configuration, runs `make -j`, and launches every registered test
through that test's own `test.py` runner. Configuration and compilation are
both recorded in `<output-dir>/build.log` and as separate steps in
`run_all_summary.json`; the M2 build has its own log and JSON receipt:

```sh
python3 srcMoon/test/run_all_test.py
```

Before each linked I-test, the script prints and flushes a `RUNNING` line with
the test ID, individual runner and log path. Completed PASS/FAIL/SKIPPED/ERROR
lines are flushed too, including when redirected to a file. `RUNNING` is only
progress text and is excluded from result counts. For concurrent I11, inspect
`I11/campaign_resources.csv` for live aggregate resources and each
`I11/<resolution>km/stdout.log` for particle-step progress.

The default build sequence is equivalent to:

```sh
rm -rf build
./Config.pl -application=moon -amps-test=off \
  -spice-path=/home/vtenishe/SPICE/cspice \
  -spice-kernels=/home/vtenishe/SPICE/Kernels
make -j "OPT3=-O3 -g"
```

The command-line optimization override retains debug symbols while restoring
production optimization; it does not edit the generated Makefile or
Makefile.conf. The resulting canonical `./amps` executable and the bounded
`srcMoon/examples/m1_regression.in` is passed to I11--I13 and
`srcMoon/examples/lola_surface.in` is passed to I14--I20. In default
regression mode. Request the bounded qualification campaign explicitly:

```sh
python3 srcMoon/test/run_all_test.py \
  --campaign-mode qualification --mpi-ranks 8
```

Qualification selects the reviewed `srcMoon/examples/m1_memory200.in` receipt.
I11 uses 500 particle steps per sample and at most 5000 particle steps per
resolution, with particle restart output disabled at build time. Ten complete
samples allow the unchanged nine-batch minimum to be evaluated once when raw
blocks are already independent, but do not guarantee equilibrium or acceptable
confidence widths. A case that reaches the cap without every criterion reports
`FAIL`; it does not authorize I13 or observational scoring. After every
completed resolution, cumulative raw diagnostic CSV/PNG/EPS/JSON artifacts are
written with an explicit unqualified/non-observational label and are excluded
from formal scoring. Full volumetric and legacy visualization text writers are
suppressed only inside linked I11 because they are not scored; normal
production sampling and the linked scoring callback remain active.

The separate physical-duration pilot is intentionally not registered in
`run_all_test.py`, because running it as I11 would confuse pilot completion
with a formal four-resolution result. It is launched explicitly with
`integration-tests/I11/run_duration_pilot.py`; see that directory's README for
the command, resource guards, artifact contract, and status semantics.
Override `--amps` or `--m1-input` only when the
alternate executable/input identity is intentional and will be captured in the
linked manifests. These M1 meshes can be substantially more expensive than
the U-series kernel checks.

Relocated installations can be selected explicitly with `--spice-path` and
`--spice-kernels`. These arguments are ignored only when `--skip-build` is
used, because that option deliberately reuses the caller's existing generated
tree.

The default output directory has the form
`test_output/srcMoon/run-all-YYYYMMDDTHHMMSSZ`. An explicit `--output-dir`
must be new or empty; the runner refuses to overwrite an earlier campaign.

Use `--skip-build` only when the existing generated tree is known to match the
source and SPICE configuration. If a qualification selection includes I13 and
I18, I19, or I20 with `--skip-build`, `--m2-amps` must identify a separately
built Na+Ar executable; the canonical `--amps` remains Na-only. `--list`
prints the explicit campaign
registry, `--only UXX` or `--only IXX` may be repeated for a focused runner
check, `--stream` mirrors individual logs to the terminal, and
`--case-id-prefix` changes only the live progress label used by a parent
orchestrator (never the local result IDs). `--strict-skips`
makes any
SKIPPED test fail the aggregate command. Without `--strict-skips`, SKIPPED
remains a recorded incomplete capability but does not change the campaign exit
status code. The machine-readable `campaign_status` remains `SKIPPED` whenever
one or more prerequisites/capabilities are skipped and no FAIL/ERROR takes
precedence; a zero developer-workflow exit never relabels incomplete evidence
as PASS.

The final terminal summary reports the number of tests executed, passed,
failed, skipped, and ended with an error. For every FAIL or ERROR it prints the
test ID, the exact individual runner, and that runner's combined stdout/stderr
log. The same information is written to
`<output-dir>/run_all_summary.json`; individual logs are under
`<output-dir>/runner_logs/UXX.log` or `IXX.log`. A normally completed linked
runner keeps `receipt.txt`, `stdout.log`, `stderr.log`, `run_manifest.json`, an
input snapshot, and `result.json` under `<output-dir>/IXX/`. An infrastructure
ERROR may retain only the subset created before failure. The aggregate JSON
records the artifact directory and manifest path without flattening or
overwriting the individual evidence. Its `artifacts` map indexes only files
that actually exist, including partial evidence from an ERROR run. A reported
PASS, FAIL, or SKIPPED missing any mandatory linked artifact is reclassified
as ERROR.

The aggregate exit values are `0` when no test failed or errored, `1` for an
evaluated FAIL (or a strict-skip rejection), and `2` for ERROR. A Moon/SPICE
configuration failure or clean-build failure is a campaign infrastructure
ERROR: no individual test is counted as executed, the failed step is named in
`run_all_summary.json`, and diagnostics are preserved in
`<output-dir>/build.log`.

### Adding a test to the aggregate campaign

Adding a `UXX_*` or linked I-series directory is not sufficient. Every new permanent
test must also be added explicitly to the `TESTS` registry in
`run_all_test.py` as its
stable ID and relative `test.py` path. The runner performs a preflight
comparison between that registry and all `U[0-9][0-9]_*` directories and exits
with ERROR if either side is missing. Permanent IDs must first be added to the
roadmap/test registry consistently, as required by `AGENTS.md`.

Each U-series `test.py` runner must accept these arguments:

```text
--build-dir <generated-AMPS-build>
--output-dir <campaign-output-root>
```

It must write exactly `<output-dir>/<UXX>/result.json`. That JSON object must
contain `id` equal to its registered `UXX` and `status` equal to one of `PASS`,
`FAIL`, `SKIPPED`, or `ERROR`. It should also contain a descriptive `name` and,
for SKIPPED or ERROR, a precise `reason`. The process exit must agree with the
JSON status: `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. Missing, malformed,
stale, ID-mismatched, status-invalid, or exit-inconsistent output is classified
as ERROR by `run_all_test.py`; the individual runner log is retained for
diagnosis.

The aggregate removes an existing per-test `result.json` immediately before
launching that test, so a crashed runner cannot be credited with a stale
result. Tests must keep all other artifacts beneath their own
`<output-dir>/<UXX>/` directory and must never write into validation-data
`raw/` directories.

The registered I11-I20 runners instead accept `--amps`, `--input`,
`--output-dir`, `--campaign-mode`, and `--mpi-ranks`, and write
`<output-dir>/<IXX>/result.json` using the same
status/exit-code protocol. They also preserve the linked artifacts required by
`AGENTS.md`. Successful process launch or a wiring invariant never substitutes
for science evidence. I11 now runs and scores the complete four-resolution
campaign. I13 remains `SKIPPED` until I11 passes and both its sampled-science
CSV and deterministic mesh/species/source/configuration JSON candidates are
reviewed and frozen; the runner never installs either reference automatically.

### Lower-level runners

List contracts:

```sh
python3 srcMoon/test/run_tests.py --list
```

Run all currently specified U tests:

```sh
python3 srcMoon/test/run_tests.py --all \
  --output-dir test_output/srcMoon/unit
```

Run selected tests or an individual directory runner:

```sh
python3 srcMoon/test/run_tests.py --test U03 --test U12
python3 srcMoon/test/U03_lunar_gravity/test.py
```

Use `--strict-skips` when automation must reject an incomplete campaign.  By
default, skipped capabilities are recorded but do not make the runner fail.
Exit values are `0` for PASS, `1` for FAIL, `2` for ERROR, and `77` for an
individual SKIPPED test.

The optional `--rebuild` flag first removes `build/` and then runs `make -j`, as
required by this checkout's build discipline.  It is not CMake-based.  The
canonical application regression remains:

```sh
rm -rf build
make -j test_Moon
ls -l test_Moon.diff
```

The required regression result is a zero-byte `test_Moon.diff`.  Preserve its
build log under `test_output/srcMoon/`.

The reviewed four-rank reference was updated on 2026-10-11 for the production
force-aware mover after two complete runs produced identical bytes. Its
SHA-256 is
`abc1d9c91f61a1c44e06e98bde8739b2de1650b7a6caef783bb97081f2eb80b4`.
The prior no-force reference hash was
`9386ea2d353db638ce84bf1f202d55cb7276443f800938c6365807be8a8deef3`;
that physical configuration is intentionally retired rather than preserved as
a legacy option. The campaign runner only compares output and never rewrites
the reference. Any later nonzero diff therefore remains a reviewed FAIL, not
permission for an automatic baseline update.

To run only the currently executable local checks:

```sh
python3 srcMoon/test/run_tests.py \
  --test U02 --test U03 --test U04 --test U05 --test U07 \
  --test U08 --test U12 --test U14 --test U22 \
  --output-dir test_output/srcMoon/kernel-tests
```

U08 returns exit code 77 while D01 qualification files are absent even when
its nested kernel result is PASS. Inspect `U08/result.json` rather than
interpreting that exit as a numerical failure.

U04 additionally requires the SPICE-enabled Moon configuration and exact
frozen kernel bytes. Configure and clean-build before running it:

```sh
./Config.pl -application=moon -amps-test=off \
  -spice-path=/home/vtenishe/SPICE/cspice \
  -spice-kernels=/home/vtenishe/SPICE/Kernels
rm -rf build
make -j
python3 srcMoon/test/U04_rotating_frame/test.py \
  --output-dir test_output/srcMoon/U04
```

The runner verifies every SHA-256 declared by the U04 contract before calling
CSPICE. It writes `kernel_hashes.json`, component-wise `term_vectors.json`, and
the complete `result.json`. Missing kernels mean `SKIPPED`; bytes that disagree
with the contract mean `ERROR`.

### Linked M0 baseline

The linked `I01`-`I10` contracts and individual runners are under
`integration-tests/`. The campaign entry point performs separate clean builds
for the normal SPICE baseline, the no-SPICE two-body fixture, and the dedicated
no-SPICE `Na,NA_PLUS` fixture. It also runs the complete U-series registry,
reruns U06 in the ion fixture, runs the canonical `test_Moon` regression, and
restores the normal SPICE build at the end:

```sh
python3 srcMoon/test/run_tests.py --phase baseline \
  --output-dir test_output/srcMoon/M0-baseline
```

Equivalently, call `run_integration_tests.py` directly and optionally provide
relocated `--spice-path` and `--spice-kernels` roots. Every output directory
must be new or empty; scientifically meaningful earlier runs are never
overwritten. The final `baseline_summary.json` reports separate local-
prerequisite, I01-I10, and combined M0 gate states.

The M0 orchestrator passes `OPT3=-O3 -g` to `make -j`, including the
canonical nightly target, to retain debug symbols and the optimization used
for the diagnostic executable. Exact commands and compiler output are kept
in `build_logs/`. Its nightly run lives under the campaign's `nightly/`
directory, preserving earlier reviewed root `run_test_Moon/` evidence. The
generated root `Makefile` remains a local configuration artifact.

Individual linked runners require an executable built for their declared
configuration. For example:

```sh
python3 srcMoon/test/integration-tests/I02/test.py \
  --amps ./amps --output-dir test_output/srcMoon/I02
```

Each I-test writes `stdout.log`, `stderr.log`, `receipt.txt`, `result.json`, and
`run_manifest.json`. I01 additionally writes an identical-control repeat and a
single-control-change manifest. I09 freezes seed 731993 and classifies a finite
domain exit as physical lunar escape only when the Moon-relative inertial
specific energy is positive; LSO velocity is transformed as
`v_inertial = v_LSO + Omega × r` before applying that criterion.

U01 runs its inside-sphere guard in a separate process. Acceptance requires
the exact production fatal diagnostic and a zero-byte, incomplete receipt.
The shell return code is only recorded: the legacy AMPS exit helper invokes
`MPI_Abort` with error code zero, which OpenMPI may expose as process status
zero even though the guard terminated execution.

### M1 numerical gates

The reviewed FULL spherical input for M1 is
`srcMoon/examples/m1_sphere.in`. I11 and I13 fail preflight if given QUICK or
CUSTOM because neither may substitute for the frozen convergence campaign.
Focused linked mesh and sphere-inventory runs
use:

```sh
python3 srcMoon/test/integration-tests/I11/test.py \
  --amps ./amps --input srcMoon/examples/m1_sphere.in \
  --output-dir test_output/srcMoon/I11
python3 srcMoon/test/integration-tests/I12/test.py \
  --amps ./amps --input srcMoon/examples/m1_sphere.in \
  --output-dir test_output/srcMoon/I12
python3 srcMoon/test/integration-tests/I13/test.py \
  --amps ./amps --input srcMoon/examples/m1_sphere.in \
  --output-dir test_output/srcMoon/I13
```

For an approximately 20-minute I12 production-wiring run, use the explicitly
coarser QUICK receipt:

```sh
python3 srcMoon/test/integration-tests/I12/test.py \
  --amps ./amps --input srcMoon/examples/m1_quick.in \
  --output-dir test_output/srcMoon/I12-quick
```

The committed QUICK targets are three times the FULL values. The first
factor-two candidate reached 1,074,788 blocks and was stopped after 28.7
minutes at 223.7 GiB RSS in an eight-rank I12 run. Factor three crosses the
next coarser dyadic AMR level. Its cold eight-rank run completed in 1,409.9
seconds (23.50 minutes), with 191.9 GiB peak RSS; the result was correctly
`SKIPPED / NOT QUALIFIED` after every wiring criterion passed. QUICK does not
establish M1 convergence or a FULL I12 result.

They are also included in the default aggregate campaign. Run only the linked
M1 bounded-regression set while reusing a verified SPICE build with:

```sh
python3 srcMoon/test/run_all_test.py --skip-build \
  --only I11 --only I12 --only I13 \
  --campaign-mode regression --mpi-ranks 8 \
  --amps ./amps --m1-input srcMoon/examples/m1_regression.in \
  --output-dir test_output/srcMoon/M1-linked
```

This command must yield two `SKIPPED` non-executions plus the evaluated I12
status; it is not M1 completion evidence. Use `--campaign-mode qualification`
with `m1_sphere.in` for the formal gates.

I11 and I12 are complete executable gates. I13 evaluates the analytical sphere area
even before a baseline exists, but reports `SKIPPED` until a passing I11
campaign supplies equilibrated profiles, return flux, column density, and
brightness and both the science CSV and deterministic configuration JSON are
reviewed. The runner never generates or updates either reference automatically.
I11 preserves separate derived inputs,
logs, manifests, production block/criterion metrics, and resource histories for all four
levels. Unequal-epoch, unequal-seed, or unequal-layout runs are not numerical-
convergence evidence.

## Current implementation status

| ID | Scope | Runner state | Meaning |
|---|---|---|---|
| U01 | LOS geometry | implemented linked executable fixture | Production ray limits, record layout, and the inside-origin guard are evaluated after full initialization. |
| U02 | target isolation/wiring | implemented source guard | Configuration evidence only. |
| U03 | lunar gravity | implemented linked probe | Point-mass kernel at a closed-form point. |
| U04 | rotating frame | implemented linked probe | Exact J2000/LSO/MOON_ME_DE421 names, kernel hashes, differential gravity, fictitious terms, transform derivative, and orthogonality are checked at the frozen epoch. |
| U05 | Na radiation pressure/shadow | implemented linked probe | Earth-umbra logic, inverse-square invariant, and absolute Combi-curve agreement to frozen graphical uncertainty; full mover gating remains I07. |
| U06 | Lorentz force | implemented for dedicated ion fixture | Normal Na-only builds honestly SKIP; the M0 runner rebuilds with Na/Na+ and evaluates the production force callback. |
| U07 | refinement selector | implemented linked probe | Legacy no-input callbacks and the input-driven radial, SZA, terminator, and body-fixed polar production kernels are evaluated. |
| U08 | LOLA geometry | implemented linked probe; PASS | Production parser/decoder/icosphere/writers and the strict receipt are exercised; the runner parses D01 readiness metadata and independently re-hashes the immutable IMG/LBL. Linked I14 remains separate. |
| U09 | terrain shadow | implemented linked probe; PASS | Production BVH ray/facet horizon service passes flat, ridge, backface, self-hit, and blocking fixtures; I15 owns D02 comparison. |
| U10 | Diviner temperature | implemented linked probe; PASS | Strict processed-D03 schema, node reproduction, periodic axes, independent trilinear fixture, native RDR UTC/flag-mask logic, radiance-coordinate calibration inverse, Paige C3--C9 substitution, and grey bolometric operator pass; I16/I17 own observations. |
| U11 | thermal inertia | implemented linked probe; PASS | Depth-resolved nonlinear finite-volume solver passes forcing, convergence, continuity, and exact enthalpy-form energy audit; I17 owns Diviner histories. |
| U12 | surface interaction | implemented linked probe; PASS | Na/Ar sticking and residence kernels plus the selectable linear-speed/kinetic-energy accommodation identity are checked; linked stochastic I18 remains separate. |
| U13 | cold trapping | implemented linked probe; PASS | Continuous finite residence, warming release, inventory conservation, classification, and checkpoint round trip pass; I20 owns D02/D03 linkage. |
| U14 | Na photoionization | implemented linked probe | Legacy constant lifetime only. |
| U15 | multispecies photoionization | full gate SKIPPED; supporting services implemented | Complete D04 registry, frozen 12-family photon statistics, conservative primitives and runtime controls execute; physical product mapping/kinetics and unchanged predecessors remain unqualified. |
| U16 | electron impact | implemented; rate/control PASS | Native D05 Maxwellian/Voronov coefficients, runtime installation, actual callback and registered orchestration controls pass; this does not qualify I23 particle kinetics. |
| U17 | charge exchange | full gate SKIPPED; total-rate services implemented | He/Na source-identified rates and independent source/domain checks execute; differential/state-resolved kinetics and physical bath/event providers remain missing. |
| U18 | plasma driver | SKIPPED | D06 reports ERROR and no qualified driver is wired. |
| U19 | helium source | SKIPPED | Source file exists but He is disabled and reservoir physics is unresolved. |
| U20 | neon source | SKIPPED | Source file exists but Ne is disabled. |
| U21 | argon source | local uniform kernel verified; full ID SKIPPED | The species-specific continuous source conserves its 3.8e20 atoms/s reference-sphere area partition and rejects unsupported maps. D10 is regional/QA-incomplete, Ar-linked stochastic injection and transient/K-weighted modes are absent, so this is not source validation. |
| U22 | Na sources | implemented linked probe and configuration guard | M0 built-in impact-on/`MySource`-off invariant, impact normalization, PSD density kernel, and exact Tenishev-2013 Case-1 parser/runtime scalars; no linked prescribed-reservoir or night-to-day migration claim. |
| U23 | meteoroid driver | SKIPPED | No production forcing path. |
| U24 | H2O/OH chemistry | SKIPPED | Source-bound OH/H₂O residence and finite-time release run as a supporting check. Molecular surface and chemistry dispatch remain open. [Scope and r14 lock](U24_water_chemistry/README.md) |
| U25 | data provenance | implemented 2026-11-23; aggregate SKIPPED | Evaluates the 12 datasets declared by U15/U16/U18/U21/U23/U26 (`--dataset`/`--for-test` select one): QA readiness, frozen observation splits, hash manifests, raw inventory, provenance hashes, and a rerun of the package's own processing on a private copy. D04/D05 PASS; D08 has an unfrozen hold-out; the other nine are not QA-ready. Data-package QA only, never validation. See `U25_data_provenance/README.md`. |
| U26 | observation geometry | SKIPPED | The M1 observation-grid settings are parsed strictly, but the complete two-limb, body-fixed pole/shell, finite-beam, and exposure-time operators are not yet implemented or qualified; D14 lacks the radial WHAM response and numerical full raster. |

The status is capability-sensitive.  Downloaded files, directory presence, or
a source/preflight check never changes a missing linked capability into PASS.

## Results and reproducibility

Results go under `test_output/srcMoon/` and never under the validation-data
repository. Each test writes `result.json`; `run_tests.py` writes
`summary.json`, while the individual-runner campaign writes
`run_all_summary.json` and one log per registered runner. Probe compilation
details are retained in the result. Linked I-series runs require the larger
manifest/artifact set specified in the roadmap; I11-I20 preserve that set when
launched by either their individual runners or `run_all_test.py`.

### M2 linked gates

I14-I20 are registered after I13. Each runner checks the I13 `result.json` in
the same campaign directory before launching AMPS. If that predecessor is
absent or not `PASS`, it writes `application.in`, empty `stderr.log`, an
explanatory `stdout.log`, `run_manifest.json` with `execution=not-run`, and
`result.json=SKIPPED`. This prevents the passing U09-U13 kernel layer or an old
I13 artifact from being promoted to current linked evidence.

The bounded M2 preflight with reused builds is:

```bash
python3 srcMoon/test/run_all_test.py --skip-build \
  --m2-amps /path/to/amps-moon-m2-na-ar \
  --only I14 --only I15 --only I16 --only I17 \
  --only I18 --only I19 --only I20 --mpi-ranks 8
```

Formal execution uses `--campaign-mode qualification` and includes I11-I13 in
the same command. I14 requires D01; I15 requires D01+D02; I16/I17 require the
appropriate D03 products; I20 requires D02+D03. Directory presence is not
readiness: package README, provenance, QA, hashes, native labels, and the exact
processed product are checked. I18 uses the input-selected Tenishev et al.
(2013) Na law and Kegerreis et al. (2017) Ar law through a dedicated Na+Ar
build, but remains predecessor-gated until the same FULL campaign produces
I13 PASS. Missing physics or evidence is `SKIPPED`, never `PASS`.

When an I14-I20 science criterion is actually evaluated, its runner must
preserve the row-level values as `comparison.csv` and render those exact rows
as both `comparison.png` and `comparison.eps`. The strict renderer is
`common/write_comparison_artifacts.py`; it rejects non-finite values, mixed
units within one comparison ID, negative uncertainty, and undeclared selection
labels. Multiple named comparisons are rendered as separate panels without
changing their row order. It never interpolates or changes the CSV. A
`SKIPPED`/not-run test must not
create placeholder comparison files, because a graphic without evaluated
evidence could be mistaken for validation.

See [`TEST_GAP_PROPOSALS.md`](TEST_GAP_PROPOSALS.md) for missing coverage.  No
new permanent IDs are assigned there.

### M2 reference-data inventory (audited 2026-10-13)

D02 currently contains immutable north/south PDS3 IMG/LBL pairs for
`AVGVISIB_65N/65S_240M_201608` and `LPSR_65N/65S_240M_201608`. The native
labels define 6420-by-6420 signed 16-bit little-endian grids, 240 m/pixel polar
stereographic projection, planetocentric latitude, positive-east longitude,
the DE421 mean-Earth/polar-axis body-fixed frame, and the exact DN scale and
offset. AVGVISIB covers the 2009-07-13 through 2013-07-18 calculation interval
and maps average solar-disc visibility continuously toward zero; LPSR is the
binary permanent-shadow classification. D02's immutable raw package remains
separately labelled `PARTIAL`; the committed `prepare_d02_reference.py` now
reproducibly creates campaign-owned `illumination_reference.csv` and
`d02_processed_qa.json` files on a predeclared approximately 100 km lattice
without changing `raw/`. This makes the derived reference inspectable but
does not supply the missing production extended-Sun temporal operator required
for I15.

D03 currently contains all 18 version-1.0 global cumulative-product latitude
bands (`global_cumul_avg_cyl_*_002`), both version-1.0 polar resource products
(`dlre_prp_north/south`), and all 144 ten-minute ZIP/XML RDR products for each
frozen 2011-06-15 training and 2011-12-10 hold-out date. The Diviner team
identifies these as RDR processing version 3.1; the individual PDS4 product
labels retain product version 1.0. The GCP labels provide 0.5-degree latitude/longitude
by 0.25-hour local-time bins, calibrated channel-3 through channel-9 and
bolometric brightness temperatures in kelvin, positive-east/planetocentric
coordinates, and `-9998/-9999` invalid/unknown constants. PRP supplies annual
average 2-cm temperature, annual maximum surface temperature, and ice-depth
fields on polar triangles. RDR labels provide UTC/JD, geometry, channel,
detector, radiance/brightness-temperature fields, and native quality/status
metadata. The committed GCP processor creates a hashed
`temperature_table.csv`, hold-out `benchmark_points.csv`, and a processing-QA
record; the PRP processor likewise creates a hashed polar reference. The I17
processor selects native footprints/flags, aggregates radiance before
temperature, applies the published Paige C3-C9 sensitivity/T_bol algorithm,
and derives uncertainty from calibration plus within-bin dispersion. GCP/PRP
campaign records and the canonical `diviner/processed/i17` RDR product qualify
only their deterministic transformations. They do not make D03 READY: the
canonical `diviner/` package currently has no
`README.md`, `provenance.json`, or `qa_report.json`. The linked runner checks
those three package records before reading a processed D03 table and therefore
reports I16/I17/I20 as `SKIPPED / NOT VALIDATED`. I17 now has a local
multi-channel/grey-bolometric production operator, reproducible RDR processor,
input-controlled finite-disc eclipse forcing, observation-free periodic
full-column spin-up, and a linked June-only training comparison path.
Canonical processed RDR evidence is written only by an explicit
data-maintenance command under `diviner/processed/i17`; the linked runner
verifies and stages it but never rewrites the validation-data package. Before
threshold freeze, staging exposes only the 2,039 June rows to the executable.
A completed run writes `comparison.csv`, `comparison.png`, `comparison.eps`,
`channel_diagnostics.csv`, and a review-required
`june_threshold_candidate.json`, then reports `SKIPPED / NOT VALIDATED` rather
than PASS. I17 still requires D03 package qualification, a same-campaign FULL
I13 PASS, investigator review/freeze of the June threshold, and a later
untouched December execution. Raw RDR presence or an experimental plot alone
never establishes I17.

The canonical 2026-10-13 I17 processing receipt is `READY` for the deterministic
RDR transformation itself: 77,075,901 records from 87 archives were examined,
25,848,333 passed the frozen native cuts, and 6,353 complete C3--C9 bins were
published (2,039 June training and 4,314 December hold-out).  The history,
contract, and processor SHA-256 values are respectively
`870577c82b430b098d816d210b528c714ec6cfa3a20ae9be3f6482a3639ae182`,
`88fc5196fee5863691937a505873d420990b98056db8a1e8167255cd736189b0`, and
`5a5d29d0f6dba88cd49c899e541a980eb84ce24717481b0cad471f8012459a65`.
This receipt does not override the missing D03 package qualification or linked
I17 capability gates and therefore is not observational validation evidence.
# M2 effective compile-configuration integration contract (2026-11-17)

The bounded I19 build has two authoritative compile-input controls for
photochemistry: `PhotolyticReactions=off` and the application block's explicit
`define _PIC_PHOTOLYTIC_REACTIONS_MODE_ ..._OFF_`. Both must select OFF only in
the derived M2 application; canonical M0/M1 remains ON. The upstream guarded
insertion order makes an unmatched late application definition effective, so
finding an earlier textual OFF definition cannot establish the compiled mode.

Before compiling and again after compiling, the builder must preprocess the
normal `pic.h` then `Moon.h` production include chain and require effective
photochemistry OFF, individual weights ON, and Mersenne Twister. This is a
compiler configuration check, not a physics or observational PASS. Preserve
the probe, complete macro output, command, hashes and executable binding. A
contradictory late definition must be rejected before the native build. Bind
the builder recipe itself in source identity and require this effective-mode
proof in both live-tree and retained-binary verification; legacy receipts
containing only a textual OFF match are insufficient. Retain the generated
`Exosphere.dfn` alongside the existing configuration/layout/RNG snapshots,
because that application header can override the global chemistry macro.

Protocol fixtures use temporary minimal headers and the real C++ preprocessor
to test late contradictory definitions and missing proof. They are explicitly
infrastructure fixtures, never physical observations. Existing I19 physics,
numerical horizons, seeds, inputs, and thresholds are unchanged. The first
upstream integration campaign remains ERROR; a corrected run is new evidence.

The independent M3 services are included in the existing U15 supporting
children (declared piecewise exposure, the optical-clock bridge, common-time
clock and qualification-protocol controls). U16 also discovers the protocol
controls through its orchestration suite. These child PASS results do not
complete the frozen U15 physics scope. The extended `m3-development` campaign
also runs the actual common-time source/mover fixture against the same fresh
native producer, before ordinary Moon restoration. All children consume the
same remaining wall allowance, preserve runner/log/result paths and propagate
FAIL/ERROR; no new permanent test ID was assigned.

Read [`M3_PHYSICAL_MODEL_SPECIFICATION.md`](M3_PHYSICAL_MODEL_SPECIFICATION.md)
and [`M3_QUALIFICATION_PROTOCOL.md`](M3_QUALIFICATION_PROTOCOL.md) for concrete
model inputs, missing distributions, unchanged criteria, proposed cohort
planning and exact preparation commands. `m3_protocol.py` emits a structural
preparation receipt and honest SKIPPED formal status, never a physical run.
The source/data roadmap hashes and controlling differences remain explicit.

The extended common-time checks distinguish manual production-hook integration
from full normal initialization. After a successful ordinary Moon restore the
campaign binds source-before/after, actual successful Config/make log, archive,
retained executable and generated headers in `normal_build_evidence.json`.
`run_m3_transport_clock.py --normal-startup-only --build-evidence ...` then
runs both actual `amps_init`/two-`amps_time_step` paths at eight ranks without
another application build. It chains the existing event observer to preserve
production budgets. Actual furnished SPICE dependencies, sampling/source order,
common endpoints and frozen environmental epoch are retained. This consumes
the same original campaign deadline; it adds no physical qualification claim.

The source-backed [chemical database/Hartley2 research report](M3_CHEMISTRY_DATABASE_RESEARCH.md)
records available water/OH photorates, activity-dependent comet tables, electron
and capture candidates, exact reuse boundaries and source interpretation gaps.
It is a research audit, not a new U/I test or a qualified physical provider.
Acquisition manifests preserve failed/challenged downloads as well as usable
native numbers. Canonical data, frozen rate selections and test gates remain
unchanged. The existing comet steady-flux executor and neutral-only ionization
sinks are not the Moon physical-time conservative packet implementation.

[`M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.md`](M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.md)
and its [JSON register](M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.json) give nine
concrete choices approved on2026-11-18, source limitations and independent
oracles. P01–P09 are decision labels, not test IDs. The proposal itself added
no runner; the approved implementation below adds supporting children to
existing owners. Reduced photon/three-body/CEX models and the ideal ion
collector require explicit selections and declared experimental scope;
they cannot qualify unsupported lunar kinetics or material recycling.
The existing frozen physical/protocol JSON and U/I criteria remain unchanged.

## Approved P01–P09 component implementation

The global U15 supporting list now includes `reduced-daughter-kinematics`,
`ideal-ion-collector` and `oxygen-state-source-controls`; the U17 supporting list includes `sodium-partial-capture`
and `helium-silica-conditional`. Each invokes its actual maintained runner,
preserves structured PASS/FAIL/ERROR/SKIPPED and its log/result/artifact paths,
and shares the existing bounded supervisor. Their experimental/laboratory
scope does not change the full U15/U17 or formal I21–I24 status. Surface
subchecks are named supporting components, not a newly assigned lunar
material-validation ID. The gap registry keeps that physical coverage open.

Individual runner/source commands and independent oracles are in
[reduced_kinematics](reduced_kinematics/README.md),
[sodium_partial_capture](sodium_partial_capture/README.md),
[ion_collector](ion_collector/README.md),
[helium_silica](helium_silica/README.md) and
[physical_state_sources](physical_state_sources/README.md). The separate
[O I source procedure](state_cascades/README.md) preserves native levels,
transition multipoles, missing A values and five explicit cascade levels.
Its offline controls check source processing and corruption handling; their
PASS is not a physical lifetime, branching law or cascade-transport PASS.
Source tables must exist and match their provenance; no fresh-machine PASS
is inferred from the existence of a research directory. Missing dependencies
or source data are SKIPPED; changed data are ERROR. The mass-source reader
re-extracts the derived inventory from authentic raw bytes and rejects edited
derived masses/thresholds/constants. It does not infer source uncertainties or
isotope shifts. Use an interpreter with Matplotlib for requested PNG/EPS.

The native packet-backend runner additionally invokes the actual normal Moon
sphere dispatcher and PIC collector adapter at1/8MPI. Its explicit nominal
state and identity geometry are analytical controls, not lunar thermochemistry
or a material benchmark. It checks actual local weight/correction, charge/
nuclear/energy transfer, distinct physical/geometry epochs, install refusal,
neutral byte/ledger preservation and deletion ownership. Ordinary M0/nightly
regression separately protects the original neutral path. Static orbit-OFF
terrain geometry is not an actual SPICE/LSO epoch qualification.

The [shared chemical database tests](../../test/chemical_database/README.md)
exercise the actual shared C++17 implementation, Python API, malformed inputs,
immutable source/profile/export bindings, permissions and installed/relocated
SDK. Run `python3 -B test/chemical_database/run_all_tests.py --help` for the
stand-alone protocol/data options. The registered shared-core owner is
`CORE-CHEMICAL-DATABASE-001` with case `CHEMICAL-DATABASE-PROTOCOL-001`;
this adds no U/I identifier. The nine existing Moon scalar/adapter suites use
`chemical_database_build.py` to compile the actual common production sources.
The global Moon runner continues to execute their existing owners. Actual
native M0/nightly and full CORE preservation evidence is recorded in STATUS.md
section 5.30. Catalogue/API integrity PASS does not install missing physical
providers or qualify I21–I24. Existing criteria, reviewed references, scientific
status meanings and permission/data gates remain binding. The broader
[design proposal](../../doc/CHEMISTRY_DATABASE_PROPOSAL.md) retains unimplemented
application-profile migration work, including Hartley2.

The additional [r4 source-input update](../../utility/chemical_database/PHYSICAL_INPUT_STATUS_20261119_R4.md)
is tested by the same automatically discovering shared runner. New molecular,
surface, quenching, metadata-amendment and recipe-dependency controls retain
administrative versus actual-source scopes. The shared test README documents
explicit native-file selectors and the new lock. No permanent lunar ID,
scientific criterion, profile default or formal I21–I24 predecessor is changed.

The [targeted physical-input assessment](M3_PHYSICAL_INPUT_FOLLOWUP_20261120.md)
records available source primitives, eleven-channel molecular threshold
derivations and the remaining explicit model choices. New NIST source-view
checks belong to the automatically discovering shared database runner, with
native selectors documented in its README. The optional Na fit-source audit
retains104PASS/2FAIL printed-coefficient comparisons; it is separate from
source-reader software tests and the frozen Voronov production-rate choice.

The [November 24 archive search](M3_ARCHIVE_SOURCE_SEARCH_20261124.md)
documents same-reaction sodium scattering measurements, a verified numerical
helium charge-exchange archive, and new cross-section diagnostics for the
I22 and I23 source reviews. It includes native hashes, reproduction commands
and the outstanding physics evidence required for transport.

The [shared r6 release](../../utility/chemical_database/NATIVE_M3_ARCHIVE_SOURCES.md)
adds permission-cleared numerical He⁺–He charge-transfer data and the live
IAEA corroboration. It is selected by its own database lock. The shared
database checks cover exact source rows and API behavior. The Moon I24
sodium angle and physical-error gates remain open.

The [r7 literature follow-up](../../utility/chemical_database/NATIVE_M3_LITERATURE_SOURCES.md)
adds source-specific water onsets and ion signals, neutral H/He integral
ionization fits, and measured and calculated He⁺–He capture percentages.
It records the still-missing joint water fragment data, I23 uncertainty
envelope and sodium angular data. The source tables are accessible through
the shared database's explicit r7 lock.

The [r8 follow-up](../../utility/chemical_database/NATIVE_M3_PARTIAL_SOURCES.md)
adds He, Ne, Na and Ar theoretical subshell photon cross sections, a
He photoionization calculations with ground and excited residuals, and gas-phase He⁺ excitation
ratios at 200 eV. Its shared API checks cover source fidelity. The
atomic broadband daughter mapping, complete water fragment kinetics,
I23 recoil and I24 sodium angular-data gates remain open.

The [r10 laboratory source follow-up](../../utility/chemical_database/NATIVE_M3_LAB_BENCHMARKS_SOURCES.md)
adds atmospheric-pressure OH loss on fused silica and gas-phase
He⁺ ratios at 200 eV. These laboratory observables retain their
conditions. Lunar OH sticking/accommodation, excited fractions and
recoil bounds, broadband photon weighting and sodium capture data
still need source-qualified input for M3.

### Nightly coverage audit — 2026-11-23

The nightly `moon` entry in `utility/nightly/runners.json` invokes
`srcMoon/test/run_tests.py --all`, whose scope is exactly **U01–U26**. Its
26 report rows do not count every supporting method inside each owner. The
shared chemical-database tests run independently through the nightly `core`
entry and `CORE-CHEMICAL-DATABASE-001`; see their
[nightly command and scope](../../test/chemical_database/README.md).

The local `run_all_test.py --list` instead lists **40 cases**: U01–U26 plus
I11–I24. I01–I10 have a separate M0 multi-build runner. The nightly Moon entry
does not invoke that aggregate or its linked I11–I24 cases. A further registry
gap remains: the global `APPLICATION-MOON-REGRESSION-001` adapter in `src/test`
still exposes only 36 cases, U01–U26/I11–I20, in the extended/application-smoke
campaigns. I21–I24 are present locally but absent from that adapter. These
findings do not grant any linked qualification or change the nightly selection;
STATUS.md section 5.42 records the archive evidence, fresh CORE execution and
remaining source/release/linked coverage gaps.

## Source-bound M3 profiles, trajectories and conditional cohorts

The [I21 owner](integration-tests/I21/README.md) accepts an explicit plan and
its SHA-256 through this aggregate. It runs the same native Na owner from
the producer checkout. The existing registry still uses literal I21, with
its homogeneous development results and formal status reported separately.
Its candidate qualification mode enforces the existing U14 and I13–I20 gates.
The [He/Ne/Ar state-profile preparation](../../utility/chemical_database/state_profiles/ATOMIC_FAMILY_CONTRACT.md)
and [U16 fit slope](U16_electron_impact/README.md) provide inputs to subsequent
I22/I23 provider integration.

The U15 supporting namespace also runs source-state profiles, moving-parent
photon composition, accepted-trajectory exposure, exact compensated clocks,
approved qualification/RNG controls and independent sodium event-ledger checks.
Each uses the same production implementation or an independent verifier of its
actual output; none replaces missing physical data. Their detailed contracts
are in [physical_reaction_events](physical_reaction_events/CONTRACT.md),
[trajectory_exposure](trajectory_exposure/README.md),
[qualification](qualification/CONTRACT.md),
[semantic_cohort](semantic_cohort/CONTRACT.md) and
[sodium_cohort](sodium_cohort/README.md). Global registration retains child
status, runner, log, result, comparison CSV/PNG/EPS and source identity. Add
future supporting runners to `m3_supporting_services.RUNNERS` and the relevant
existing owner's tuple in `common/moon_testlib.py`; add protocol controls rather
than assigning a permanent U/I identifier informally.

`m3_qualification_executor.py` now has an actual conditional I21 execution path.
The frozen plan must bind its complete source profile, a named installed
producer, generated build identity, every root/included `.in` path+SHA, the
explicit three-state mapping and at least three matched common-step refinements
at BOTH1/8MPI. Only after actual source and same-campaign predecessor admission
does it call the normal native owner. It independently audits every original
parent, complete per-rank event set, actual committed child, signed reservoirs,
source-clock timing, endpoint survival and exact branch confidence. The original
one-hour resource clock covers launches and all auditing/plots; it never resets
between levels. A larger experiment is deferred to a separately declared cluster
package with the same criteria. Actual partial launches survive in
`execution_journal.json` even if a later native process fails.

The default partial profile still leaves Na+ unavailable. The opt-in atomic
profile now supplies source-bound ground-ion evidence and an explicit
source-effective photon crosswalk; physical execution still requires reviewed
model-error controls and the full native admission checks. I22–I24
still lack complete source-consistent product providers and remain disabled.
Readiness or one successful cohort is not full I21, I21–I24 or M3 closure.
The approved numerical0.2/0.3percent budgets, frozen original protocol, exact
branch intervals, data gates and full qualification matrix are unchanged.
Use [the numerical contract](sodium_cohort/NUMERICAL_CONTRACT.md) for the
homogeneous reference assumptions and [ledger instructions](sodium_cohort/LEDGER_README.md)
for exact native hexadecimal input and independent reconstruction.


### Supporting M3 environment-stage and admission controls

`coefficient-stage-controls` checks exact zero-acceleration macro admission and
actual1/8MPI administrative consensus before physical startup. It is registered
with the global supporting runner and U15. Heavy normal-SPICE stage checks use
`run_m3_coefficient_stage.py` with an authentic successful normal producer;
the developer orchestration supplies that same proof after restoration.
See [coefficient_stage/README.md](coefficient_stage/README.md). None of these
software/source checks closes formal M3 or supplies missing Na+ state data.

U15's expanded supporting collection has the named600s execution allowance
`M3_PHOTON_SUPPORT_WALL_SECONDS`. It covers additional source/profile/ledger/
statistics checks; every enclosing native campaign still retains its original
one-hour deadline. This changes no scientific acceptance threshold. A previous
300s timeout is preserved as ERROR; a corrective rerun is separate evidence.


### November 19 source, reservoir and cascade additions

U15 now includes `atomic-state-profile`, `source-build-protocol`,
`physical-case-preparation`, `oxygen-cascade` and `reservoir-products`. U17 also
includes `reservoir-products`. These names live in `m3_supporting_services.RUNNERS`
and the corresponding owner tuples in `common/moon_testlib.py`, which the global
runner calls. Each child uses the normal `--output-dir NEW` protocol, preserves
method/metric results plus logs/source hashes, and returns 0/1/2/77 for
PASS/FAIL/ERROR/SKIPPED. The global supporting summary lists every runner/log and
retains child counts. Whole U/I status remains separate from supporting results.

For example, run the new case-preparation controls without Config or MPI:

```sh
/usr/bin/python3 -B srcMoon/test/test_m3_physical_case_preparation.py \
  --output-dir test_output/srcMoon/NEW-case-preparation-controls
```

See [sodium_cohort/PREPARATION.md](sodium_cohort/PREPARATION.md) for explicit
scientific inputs and [source_qualification_build/README.md](source_qualification_build/README.md)
for a disjoint clean build. Source-preflight is a real MPI inspection of source
states and compiled mass/charge, but launches no particle cohort. Do not use its
PASS as I21, or reservoir/cascade numerical checks as I22–I24 qualification.
The shared CORE test discovers reservoir-product checks independently; overlapping
local/CORE counts must not be added as unique science tests. Existing predecessors,
original-parent confidence rules and frozen numerical allocations still apply.

### Approved P02 and physical commissioning audits

`p02-independent-controls` and `commissioning-audit-controls` are additional U15
supporting groups in the existing registry. They run 30 and 20 maintained
software/source controls respectively, report all four normal statuses, and
retain per-method logs and source hashes. These controls do not fabricate a
successful physical owner. Actual event audits are explicit postprocessing of
the completed native owner, documented in
[sodium_cohort/P02_README.md](sodium_cohort/P02_README.md) and
[sodium_cohort/COMMISSIONING_AUDIT_README.md](sodium_cohort/COMMISSIONING_AUDIT_README.md).

The approved mechanical budget and controlled source inputs are documented in
[M3_P02_NUMERICAL_BUDGET.md](M3_P02_NUMERICAL_BUDGET.md). A small commissioning
cohort may verify event wiring and ledgers, but cannot replace the frozen
statistical design, both MPI layouts, three refinements or formal predecessors.
The original native timestep and rounded-affine failures are retained in the
November 19 approved-P02 package. Bounds and source rates are unchanged by their
fixes. The actual rounded-matrix regression uses an independent high-precision
inverse; it is frame arithmetic, not observational geometry validation.

Shared r11 physical input and application-specific admission decisions are
documented in the [I22](integration-tests/I22/README.md),
[I23](integration-tests/I23/README.md) and [I24](integration-tests/I24/README.md)
runner descriptions. The database software suite checks the archive contracts
and reports its results separately from the formal linked M3 gates.

For I21–I24 experimental comparisons, see the
[survey of published model choices](../../utility/chemical_database/NATIVE_M3_MODEL_PRACTICE_SOURCES.md).
It identifies which lunar and comet observations constrain neutral survival
or bulk pickup, the starting velocities other models use, and the evidence
still needed for daughter transport and OH surface interactions.

The [lunar and sibling-AMPS reference comparison](M3_LUNAR_MODEL_REFERENCE_CROSSWALK.md)
identifies the reviewed I04 output, the native and analytical I21 event
records, and the fragment observables that a matched Europa AMPS run can
provide. It also documents what independent lunar models assume about
ionization, water surface coverage and ion return.

The r14 shared database provides named reference cases for bulk lunar Na⁺
pickup and for the 4.0 keV/u ground-product charge-exchange event. The shared
database suite checks their physical units, frames, model bindings and source
provenance. The I24 runner retains responsibility for linked reaction budgets
and verification. [I24 contract](integration-tests/I24/README.md#selectable-database-reference-for-modeled-recoil)

Select the sodium I24 cohort with `--only I24 --i24-sodium-lock <r14-lock>`.
Add `--i24-pic-producer <normal_build_evidence.json>` to execute both PIC
tracing orders at one and eight MPI ranks in the same development run.
The aggregate retains each supporting outcome, launch log and native
artifact directory. [Commands and result locations](integration-tests/I24/README.md#running-both-the-mpi-cohort-and-the-linked-pic-movers)

Add `--i24-pic-refinement` to that command to include second-order checks
at three time-step and three reference-weight levels, keeping the proton
beam, parent population and final time fixed. [Levels and outputs](integration-tests/I24/README.md#time-step-and-particle-weight-refinement)

### M3 family runner

Use the [frozen-family campaign](qualification/README.md#frozen-chemistry-family)
to replay several source-bound I21 cohorts and both I24 sodium comparisons
under one wall, memory and artifact budget. It uses the existing registered
science tests and formal readiness owners. The current I22 and I23
readiness results therefore appear alongside the sodium development results.
The family plan also accepts the two I23 electron-bath kinds,
`i23-helium-ensemble` (profile pin) and `i23-helium-pic` (profile and
Na/He/He⁺ producer pins), and the two I24 total-capture kinds,
`i24-total-capture-ensemble` (case pin) and `i24-total-capture-pic` (case
and Na/Na⁺/H producer pins).

### I21 decision-1 cases

The I21 owner runs the decision-1 plans (approved 2026-11-26) through its
existing `--i21-cohort-plan` path: 1:8 and 8:1 weights, a prescribed static
Sun at 1 and 2 AU, a Moon-shadow crossing and a lunar-gravity path, each with
both movers. `test_m3_i21_decision1.py` runs in the U16 orchestration battery.
[Contract](integration-tests/I21/DECISION1_CONTRACT.md) ·
[Commands](integration-tests/I21/README.md#decision-1-cases-approved-2026-11-26)

### Python interpreter check for the I23/I24 supporting runners

The six I23/I24 supporting runners check their third-party packages through
`common/python_environment.py` after input authentication and before any
native work. A missing package is SKIPPED (naming the interpreter and
`/usr/bin/python3`); a corrupt input remains ERROR. Run manifests and the
orchestrator summaries record the interpreter. `test_m3_python_environment.py`
runs in the U16 orchestration battery: with matplotlib hidden, each runner must
report SKIPPED with no native command, and a corrupt-input stub must remain
ERROR. [Machine Python policy](../RUN.md#python-environment-for-this-machine)

### I24 total proton–sodium capture

The I24 owner accepts `--i24-total-capture` in regression mode. It runs the
user-approved (2026-11-25) capture H⁺ + Na(3s) → Na⁺ + H(n) at 1–5 keV/u as
a supporting MPI ensemble and, with `--i24-pic-producer` (and optionally
`--i24-pic-refinement`), a linked PIC matrix. No r14 lock is needed for this
case; `--i24-sodium-lock` still selects the 4.0 keV/u H(1s) fixture, and both
may be given together. Formal I24 remains SKIPPED. The generator and oracle
controls (`test_m3_i24_total_capture.py`) run inside the registered U16
orchestration battery. [Contract](integration-tests/I24/TOTAL_CAPTURE_CONTRACT.md) ·
[Source review and decisions](integration-tests/I24/SOURCE_REVIEW_20261125.md)

### I23 prescribed electron-bath coupling

The I23 owner accepts `--i23-helium-profile`, `--i23-pic-producer` and
`--i23-pic-refinement` in regression mode. They run the approved controlled
10 eV bath as two supporting children, an MPI provider ensemble and a linked
PIC matrix, while formal I23 remains SKIPPED. The admission and protocol
controls (`test_m3_i23_electron_bath.py`) run inside the registered U16
orchestration battery. The He/He⁺ `.in` key controls run in
`test_m3_runtime_configuration.py`. [Contract](integration-tests/I23/ELECTRON_BATH_COUPLING_CONTRACT.md) ·
[Approval record](integration-tests/I23/HELIUM_ADMISSION_20261125.md)
