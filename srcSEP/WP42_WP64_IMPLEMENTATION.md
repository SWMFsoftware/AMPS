# WP42--WP64 implementation notes

## Scope and evidence boundary

This directory contains dependency-light prototypes for the work packages in
numerical order. B05 deliberately does not relabel those prototypes as active
production behavior merely because they compile. The source tree remains
field-line-only and retains exactly the three public movers `parker`,
`fte-dmumu`, and `fte-mfp`.

The dependency-light gate compiles the production-independent implementations
with strict C++11 warnings, AddressSanitizer, and UndefinedBehaviorSanitizer.
It is analytical-core and source-integration evidence.  It does not claim that
native AMPS, MPI/OpenMP decomposition, real SWMF replay, held-out spacecraft
data, or multi-node performance ran in the source-only environment.

The machine-readable authority is
[`WP42_WP64_DISPOSITION.json`](WP42_WP64_DISPOSITION.json). The current
classification is intentionally conservative:

| Status | Work packages | Meaning |
|---|---|---|
| Accepted supporting data contract | `ParkerMeasure`; per-event `WaveContribution` branch and momentum fields | These small contracts repair current header compatibility and event bookkeeping. They do not promote an entire work package. |
| Experimental component API | WP42--WP58, WP61--WP62 | The four extension implementations compile under strict warnings/sanitizers and have controlled tests, but are absent from `MAINLIBOBJ` and have no production selector/call path. |
| External qualification gate | WP59--WP60, WP63--WP64 | Validators exist, but linked AMPS, decomposition/restart, real SWMF/observational, or frozen-hardware evidence must be supplied externally. |

No WP42--WP64 package is currently classified as a production-active
end-to-end feature. This avoids changing stable physics merely to satisfy a
stale overlay. Promotion requires the per-package gate recorded in the
disposition manifest and a separate review of the resulting native evidence.

## WP42 -- persistent production turbulence ownership

**Why.** Reconstructing turbulence from PIC display arrays every timestep loses
pending sources, operator phase, restart metadata, and accumulated ledgers. It
also permits two apparent owners of the same wave energy.

**Prototype.** `TurbulenceRuntimeStore` can own one `TurbulenceLineRecord` per
integer field-line ID. The record contains the complete common-core `State`,
geometry generation, owner rank, and configuration fingerprint. The production
adapter does not yet use this owner; it still imports/exports host arrays.
The prototype refreshes geometry/background coefficients without importing
wave energy and its checkpoints
hex-encode complete core checkpoints and are deserialized into staged objects
before live state is replaced.

**Contract.** A geometry-generation change is rejected while unresolved source
transactions remain unless the caller explicitly selects apply-before-remap.
The component `Serialize`/`Deserialize` API is not an application restart seam
until a configured AMPS run owns it across checkpoint/restart.

## WP43 -- atomic particle--wave exchange

**Why.** Worker-side mutation and a second legacy coupling call can double count
energy, depend on thread order, or lose work when a particle exits a boundary.

**Prototype and accepted supporting fix.** `ApplyCouplingTransactions`
validates identity, generation, cell, branch, bin, uniqueness, and available
wave energy before committing a staged prototype store; invalid batches roll
back completely. It is not called by the active production turbulence driver.
The accepted current-mover fix is narrower: each completed MFP event now emits
its actual selected branch and its immediate pre/post wave-frame momenta.
Deterministic intervals use branch zero and zero event momenta. Records are
deposited only after the complete event/no-event update succeeds.

**Contract.** PIC particle pointers, worker IDs, rank IDs, and queue addresses
are never part of transaction identity or ordering. The deprecated void flush
wrapper is retained only for source compatibility and is not a production owner.

## WP44 -- distinct velocity derivatives

**Why.** `div(U)`, `bb:grad(U)`, and `d(U dot b)/ds` are different on a curved,
expanding field line. Reusing one scalar for all three changes focusing and
adiabatic terms.

**Prototype.** `VelocityGradientInput` accepts a full Cartesian velocity
gradient, field unit vector, and curvature. `ComputeVelocityDerivatives`
calculates divergence, field-aligned strain, parallel derivative, and curvature
contribution separately. The current PIC adapter has not migrated to the full
tensor because the required provider gradient is not yet available. Density-
history experiments can report the continuity residual
`D ln(rho)/Dt + div(U)` rather than silently asserting exact continuity.

**Limitation.** The current host exposes a field-line chord reconstruction at
this seam. A coupled provider that supplies the complete 3-D gradient should
populate the same input directly and validate the residual in the native gate.

## WP45 -- validated shock state and upstream flux

**Why.** Downstream density times shock swept volume overcounts processed
particles, while replacing a slow shock with a configured minimum speed creates
unphysical source particles. Silent compression/exponent clipping masks invalid
shock states.

**Prototype.** `ShockState` records the normal, upstream/downstream normal
velocities and densities, compression, Mach numbers, obliquity, provider, and
provenance. Validation distinguishes no-compressive-shock from inconsistent
Rankine--Hugoniot data. `ProcessedUpstreamParticles` uses
`n_up A max(V_sh,n-U_up,n,0) dt`. The active analytic/SWMF injection branches
have not all migrated and are therefore not claimed fixed by B05. Their
replacement belongs to the dedicated production shock/source work package.
The prototype DSA contract requires finite `r>1` and uses `q=3r/(r-1)` without
hidden caps or floors.

## WP46 -- versioned source-event identity

**Why.** Floating timestamps, buffer handles, and shared RNG consumption make
injection change under restart, timestep subdivision, or decomposition.

**Prototype.** `SourceEventKey` contains schema, campaign, source, integer
event, field line, species, and ordinal. `StableSourceIdentity` hashes exactly
that tuple. Random purposes use separate deterministic streams, so adding a
gyrophase draw cannot move the energy sequence. The component scheduler
serializes campaign seed and next event, but the active injection/restart path
has not adopted it and must not be described as timestep-invariant yet.

## WP47 -- Parker measure and geometric drift

**Why.** A walker representing density per arc length does not obey the same
SDE as one representing density per physical volume on a variable-area tube.

**What and how.** `ParkerMeasure` is explicit. For per-arc-length walkers the
Itô drift is `U + d(kappa)/ds`; for per-volume walkers it additionally contains
`kappa d ln(A)/ds`. `WalkerToPhysicalDensity` applies the one-and-only area
conversion. Invalid area or coefficient inputs fail with typed status.

## WP48 -- signed dynamic resonance and finite-bin overlap

**Why.** The magnetostatic `k~Omega/(v|mu|)` shortcut is singular near 90
degrees and erases charge, magnetic polarity, wave propagation, harmonic, and
finite spectral-bin information.

**What and how.** `SolveDynamicResonance` retains signed charge, signed field,
relativistic gyrofrequency, wave branch, harmonic, and Alfvénic phase speed. It
returns resolved/no-root/singular/out-of-band metadata. A resolved root is
partitioned linearly between adjacent logarithmic bin centers so coupling does
not jump discontinuously at bin midpoints.

## WP49 -- controlled 90-degree scattering closure

**Why.** A hard floor or arbitrary interpolation at `mu=0` changes parallel
transport without naming the added physics.

**What and how.** `NinetyDegreeClosure` is a named, provenance-bearing Gaussian
resonance-broadening term. The implementation returns both `D_mumu` and its
analytic derivative, is C1 through zero, preserves `D(+/-1)=0`, and marks when
the closure contributes. Zero amplitude exactly reproduces the pure slab
baseline.

## WP50 -- wave action and geometric conservation

**Why.** Energy density, cell-integrated energy, and wave action transform
differently when volume and intrinsic frequency change.

**What and how.** `WaveInvariant` declares the authority. Conversion to energy
is centralized, and `ApplyGeometricConservation` preserves either integrated
energy or wave action exactly while recording background work as the difference
in physical energy. Background work remains separate from dissipation.

## WP51 -- conservative spectral cascade and heat partition

**Why.** A non-telescoping cascade can create energy internally, and energy
leaving the high-k boundary must not disappear as an unnamed sink.

**What and how.** `AdvanceConservativeCascade` consumes N cell energies and
N+1 interface fluxes for each wave branch. The finite-volume difference
telescopes exactly. High-k outflow is partitioned between explicit electron and
ion heat ledgers, and the returned closure residual measures all sources and
sinks. Negative trial energies fail rather than being silently clipped.

## WP52 -- versioned C1 analytic profiles

**Why.** Anonymous piecewise fits with implicit extrapolation introduce jumps
in derivative-driven transport and are impossible to identify in evidence.

**Prototype.** `AnalyticProfile` stores ID, semantic version, units,
provenance, ordered Hermite knots, analytic slopes, and a named extrapolation
policy. Validation rejects incomplete or non-finite profiles. Existing
production analytic fits have not been migrated; each requires reviewed knots,
units, provenance, and an explicit domain policy before promotion.

## WP53 -- Brownian-tree adaptive SDE integration

**Why.** Drawing new noise after step rejection changes the stochastic path and
biases adaptive solutions.

**What and how.** `BrownianAddress` names campaign, particle, operator, physical
interval, depth, and node. Deterministic bridge increments satisfy parent equals
left plus right. `AdvanceAdaptiveMultiplicativeSde` compares a full step with
two half steps on that same path, recursively subdivides only within configured
absolute/relative tolerances, and commits only after the whole interval passes.
It returns accepted/rejected counts, maximum depth, and normalized error.

## WP54 -- second-order monotone advection

**Why.** First-order upwind is robust but excessively diffusive for transported
turbulence gradients.

**What and how.** `AdvectPeriodic` retains first-order upwind as an explicit
baseline and adds MUSCL reconstruction with named minmod or monotonized-central
limiting and SSPRK2 time integration. It reports face work, limiter activity,
and conservation residual. The refinement test uses translated cell averages
and requires improved error with measured order above the declared threshold.

## WP55 -- conservative moving-tube remap

**Why.** Point interpolation does not conserve cell-integrated wave energy when
the field-line mesh or cross-sectional volume changes.

**What and how.** `RemapConservative` operates in monotonically increasing
physical-volume coordinates and integrates exact old/new overlaps. Piecewise
constant and limited piecewise-linear reconstructions are named policies.
Uncovered measure is reported as boundary exchange or initialized integral;
the numerical residual must close the old and new integrals.

## WP56 -- typed mover failure and commit semantics

**Why.** Numerical failure must not leave particle state advanced while its
coupling record is absent, or delete a particle before completed work is
published.

**What and how.** `MoverDisposition`, `ParticleFate`, `FailurePolicy`, and
`FailureContext` classify completion, boundary exit, rejected step, quarantine,
configuration error, and fatal internal error. `ParticleTransaction` stages a
complete state and commits atomically; rollback restores the initial value.
Context contains stable particle/background/coefficient identities and replay
coordinates rather than transient pointers.

## WP57 -- lineage-aware estimators

**Why.** Split particles are correlated descendants, not independent samples;
treating every child independently understates uncertainty.

**What and how.** Samples are first accumulated by root lineage and sampling
window, then combined into bin estimates, variance, covariance, effective
sample size, and independent-root count. Snapshot-density and crossing-flux
estimators are distinct. Statistical weight, relativistic particle speed, and
Jacobian are particle-resolved rather than replaced by a bin-center speed.

## WP58 -- multidimensional instrument response

**Why.** An energy-only efficiency cannot represent species contamination,
look direction, dead time, background, saturation, or count statistics.

**What and how.** `InstrumentResponseMatrix` is versioned and flattened over
channel, species, direction, and true energy. `FoldInstrumentCounts` applies
the response, exposure, nonparalyzable/paralyzable dead time, and Poisson log
likelihood. Saturated observations are treated as censored thresholds rather
than exact counts. Calibration validity and checksum fields are mandatory.

## WP59 -- linked native configuration traces

**Why.** Source scans prove reachability but cannot prove which callbacks a
linked executable actually entered.

**What and how.** `NativeMatrixTrace` records configuration/executable identity,
compiler flags, selected mover/provider/owners, and observed callback entries.
Validation refuses evidence below `native-amps` and requires the coupling entry
when the row requests coupling. The native runner is fail-closed and requires a
site-supplied `SRCSEP_NATIVE_GATE` command.

## WP60 -- decomposition and restart invariance

**Why.** Reproducibility claims require physical identity/order and restart
continuity, not merely similar plots.

**What and how.** `DecompositionSignature` compares configuration fingerprint,
particle identity hash, transaction order, restart hash, event/transaction
counts, and signed ledgers across rank/thread layouts. Discrete fields are exact;
floating ledgers use only the explicitly supplied relative tolerance.

## WP61 -- manufactured-equation refinement

**Why.** One favorable coarse/fine pair can hide pre-asymptotic behavior or an
omitted PDE term.

**What and how.** `ManufacturedCase` records all active terms plus at least
three resolutions and positive errors. `FitRefinementOrder` performs a
least-squares log(error)-versus-log(resolution) fit and reports order,
intercept, and R-squared. Cases without an active-term declaration fail.

## WP62 -- predeclared statistical power and distribution tests

**Why.** Choosing ensemble size, significance, or observable after viewing the
result invalidates rare-event and stochastic comparisons.

**What and how.** `PowerSpecification` freezes observable, standardized effect,
alpha, power, multiplicity, and seed-panel version. `PlanNormalMeanPower`
applies Bonferroni adjustment and computes the required sample count.
`KolmogorovSmirnovStatistic` compares the complete empirical distribution with
precomputed expected CDF values instead of checking moments alone.

## WP63 -- external cross-model and observational protocol

**Why.** Synthetic fixtures or undocumented preprocessing cannot support real
SWMF or observational claims.

**What and how.** `ExternalEventManifest` requires event role, background,
geometry, shock, and calibration checksums; preprocessing; exclusions;
configuration prior; and metrics. Synthetic records are rejected for coupled
or observational evidence, and a calibration event cannot masquerade as held-
out validation. Normalized multi-metric distance retains onset, peak, fluence,
anisotropy, spectrum, and profile-shape uncertainties.

## WP64 -- scaling and resource regression gate

**Why.** Wall time alone is noisy and incomparable across machines; performance
regressions can also appear first in work, communication, memory, or queue size.

**What and how.** `ScalingSample` records frozen workload/environment identity,
ranks, threads, item count, operations, communicated bytes, peak resident bytes,
peak queue records, and wall time. `CompareScaling` reports parallel efficiency
and per-item resource metrics, requires matching environment fingerprints for
timing, and enforces declared work and memory ratios.

## Files and commands

The implementation is in:

- `util/sep_runtime_contracts.*` for WP42--WP46;
- `util/sep_physics_extensions.*` for WP47--WP52;
- `util/sep_numerical_extensions.*` for WP53--WP56;
- `util/sep_validation_extensions.*` for WP57--WP64;
- `src/models/sep_common/sep_transport_common.h` and
  `util/sep_focused_transport_core.h` plus the MFP core/adapter for the two
  accepted supporting data-contract corrections;
- `test/test_wp42_wp64.cpp` and `test/run_wp42_wp64_tests.sh` for the bounded
  implementation gate.

The four extension `.cpp` files are deliberately absent from `MAINLIBOBJ`.
Run `make test-wp42-wp64-experimental`. Run native/external/scaling evidence separately
with `make test-wp59-wp64-native SRCSEP_NATIVE_GATE='reviewed command ...'`.
