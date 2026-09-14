# WP31–WP41 implementation and evidence contract

This document records the production contracts added for WP31 through WP41,
their physical and numerical meaning, and the evidence that can be produced
from this source-only handoff. All internal dimensional quantities use SI.
Passing `make test-wp31-wp41-unit` is analytical-core and source-integration
evidence; it is not a linked AMPS, MPI, SWMF, observational, or hardware-scaling
result.

## WP31 — turbulence splitting and stiffness control

Why: advection previously selected a CFL subcycle while reflection and cascade
used the full particle step. A fast local exchange could therefore exceed its
accuracy timescale even when advection was stable.

What: `Turbulence::PlanAdvance` estimates advection, source, reflection, and
cascade limits for every segment. It selects one finite positive stage size,
rejects steps below `minimumSubstepS` or above `maximumSubsteps`, identifies the
limiting operator, and publishes per-operator limits and work counts.

How: the declared scheme is `OperatorSplitting::LieFirstOrder`. Additive source
application and exponential two-branch reflection are individually exact over
their stage; the nonlinear cascade uses positivity-preserving backward Euler.
Noncommuting reflection/cascade and conservative advection use the common stage
count. The configuration exposes CFL, operator accuracy safety, maximum source
fraction, maximum cascade fraction, minimum stage, and maximum stages. These
values are serialized by turbulence checkpoint schema 2 and run-configuration
schema 2. The focused manufactured cascade case measures order 0.97, consistent
with the declared first-order method.

## WP32 — typed turbulence exceptional paths

Why: a floor, skipped contribution, or nonfinite-to-zero conversion can hide a
model-domain failure and break a global energy account.

What: `UpdateDisposition` classifies applied, physical-zero, regularized,
corrected, rejected, and fatal paths. `OperatorEvent` records the operator,
requested/applied joules, and reason. Diagnostics count physical zeros,
regularizations, corrections, rejections, operator stages, and deterministic
work. `EnergyLedger::rejectedSourceJ` retains the signed unapplied request.

How: negative source requests use one conservative positivity correction and
book both applied energy and the signed rejection. Advection roundoff correction
emits an event. Nonfinite rates and invalid density/volume reject the plan. The
former hidden `max(abs(mu),1e-12)` resonance floor was removed: `mu=0` now
returns an explicit out-of-domain status.

## WP33 — particle merge and split invariants

Why: population control can bias weight, energy, momentum, pitch moments, and
rare tails, and fixed production arguments made the effective policy difficult
to audit.

What: `sep_population_control.*` defines SI particle observations, population
moments, invariant reports, stable lineage identities, exact clone splitting,
and compatible-species merging. `RunConfiguration::populationControl` owns the
three phase-space bin counts, population thresholds, selection mode, invariant
tolerance, and lineage schema.

How: splitting copies phase-space state and partitions statistical weight; the
last child receives the floating residual. Tagged lineage hashing excludes rank
and thread. Merging preserves weight, charge, energy, parallel momentum, and
the first pitch moment exactly for one species/charge state. The second pitch
moment is measured explicitly because one representative `mu` cannot preserve
both first and second moments for an arbitrary distribution. `main_lib.cpp`
passes the frozen bin counts and thresholds to the actual PIC merge/split calls.
The native adapter campaign must observe PIC particles before and after those
calls to validate host-algorithm bias.

## WP34 — native production adapter harness contract

Why: a pure kernel can pass while a wrapper selects another provider, drops an
event, bypasses the turbulence driver, or flushes a coupling queue twice.

What: `sep_evidence.*` defines `NativeHarnessRequest`, `NativeObservation`, and
`ValidateNativeObservation`. A valid observation must show entry into the
registered mover, coefficient adapter, and turbulence driver; exactly one
queue-flush owner; and exact binary, run fingerprint, source generation, thread,
and rank identity.

How: the validator refuses promotion from source-integration to native evidence.
The source gate additionally checks both standalone and coupled drivers call
`PICAdapter::Advance` and that only its transaction phase calls
`FlushWaveContributions`. Completing the native physics fixture still requires
the enclosing AMPS tree and generated PIC types; absent that executable, the
native gate is BLOCKED rather than skipped or passed.

## WP35 — generated configuration compatibility matrix

Why: mover, coefficient source, turbulence ownership, and coupling can be valid
separately but invalid together.

What: `sep_configuration_matrix.*` enumerates all 90 combinations of three
movers, three coefficient sources, five turbulence sources, and two coupling
policies. Every row is supported, conditional, or unsupported and has a stable
diagnostic plus an alternative where applicable.

How: CLI parsing and frozen-run validation call the same `Preflight` function
before model initialization. Streaming coupling with immutable waves,
self-consistent coefficients without local waves, and SWMF coefficients without
an SWMF generation fail closed. SWMF-to-local handoff is conditional on its
explicit phase. `make print-configuration-matrix` renders the authoritative
Markdown table directly from the registry; the focused test proves it contains
exactly 90 data rows.

## WP36 — global conservation restart and remap ledger

Why: subsystem ledgers can close while energy or particles disappear at source,
escape, coupling, remap, checkpoint, or output seams.

What: `sep_system_ledger.*` defines signed number, charge, energy, and parallel
momentum transactions with seam and provenance. It closes final state against
initial state plus every transaction. `CampaignCheckpoint` serializes field-line
generation, stable particle IDs, lineage generations, RNG counters, residual
hazards, turbulence hash, output checksum, and the global ledger.

How: canonical checkpoint schema 1 ends with a checksum over all preceding
bytes. Truncation or corruption is rejected before state publication. The
focused campaign closes injection and escape exactly, round-trips all identity
state, and proves a mutated checkpoint fails. A linked campaign must populate
the same transaction type around actual PIC mutation seams.

## WP37 — multi-seed statistical validation

Why: bitwise repeatability for one seed does not measure Monte Carlo bias,
variance, confidence coverage, or test power.

What: `sep_validation_tools.*` creates versioned, observable-domain-separated
seed panels and uses stable online mean/variance accumulation. `CompareMean`
records estimate, standard error, z score, confidence level, and sample count.

How: fixed-seed regression remains separate. Ensemble gates require at least
two samples, nonzero standard error, an explicit confidence level, and a
predeclared maximum absolute z score. The focused case verifies deterministic
panels and domain separation. Scheduled production panels remain a distinct
native evidence artifact.

## WP38 — property generation and failure injection

Why: named cases alone do not cover IEEE endpoints, degenerate geometry,
malformed schema, and partial adapter failures.

What: bounded case generation always includes zero, signed zero, pitch
endpoints and adjacent values, finite extremes, infinities, and NaN before
pseudo-random values. `CheckProperty` returns the first exact seed/index/value
counterexample. `FaultInjector` deterministically fires at a configured hit for
allocation, file, checksum, snapshot, coupler, or communication seams.

How: the focused case proves stable counterexample reproduction and third-hit
fault behavior under ASan/UBSan. New repaired branches should add their property
callback to this runner; stable failures should become named regression cases.

## WP39 — observation forward model and external gate

Why: simulation differential intensity is not an instrument count until energy
response, aperture, cadence, species, detector effects, background, and
uncertainty are applied.

What: `sep_observation_forward_model.*` folds an SI differential-intensity grid
through channel energy overlap, dimensionless response, geometric factor
`m² sr`, and cadence. It applies a nonparalyzable dead-time model, saturation
flag, background subtraction, and combined Poisson/background/model variance.

How: a two-bin synthetic spectrum yields analytically known counts and
uncertainty. Observational manifests now require a structured forward-operator
record. Real held-out event status still requires checksum-verified spacecraft
inputs and the existing six metric families; source-only synthetic results
cannot satisfy it.

## WP40 — performance and scaling gates

Why: more accurate hazards, reductions, ledgers, and spectral operations add
cost, but physics tolerances must never be relaxed to meet a runtime target.

What: `WorkCounters` covers mover stages, rejected steps, hazard/provider calls,
queue peak, reduction bytes, spectral updates, output/checkpoint bytes, and peak
resident bytes. The turbulence core/adapter exports face, cell, spectral,
limiter, rejection, and stage counts.

How: deterministic work must match exactly on every host. Wall time is enforced
only when sample and baseline share an explicit hardware/compiler environment
identity; otherwise it is recorded but does not create a false regression.
Native strong/weak scaling remains BLOCKED without the linked MPI/OpenMP build
and a reviewed hardware baseline.

## WP41 — evidence-level documentation governance

Why: a source grep, pure analytical test, native run, SWMF replay, and held-out
observation support different scientific claims.

What: the evidence vocabulary is `analytical-core`, `source-integration`,
`native-amps`, `swmf-replay`, and `observational-validation`. Claims state the
required and observed level plus PASS, FAIL, BLOCKED, or INCOMPLETE. A PASS below
its required level or without a reproducible command/artifact is invalid.

How: README tables use this vocabulary. The generated configuration matrix and
source checks consume production registries/call sites. The focused gate rejects
an overstated native claim and checks exactly one coupling-queue flush owner.

## Commands and evidence boundary

```sh
make test-wp31-wp41-unit
make print-configuration-matrix
make test-sanitizer
make test-native-amps-validation SEP_EXECUTABLE=/path/to/amps
make test-swmf-validation SWMF_MANIFEST=/evidence/swmf.json
make test-observational-validation \
  OBSERVATIONAL_MANIFESTS="/evidence/event-1.json /evidence/event-2.json"
```

The first two commands are available from this archive. Native AMPS,
decomposition/scaling, SWMF replay, and observational validation require their
named external dependencies and remain separate release gates.
