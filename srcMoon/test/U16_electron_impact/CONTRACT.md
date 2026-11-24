# U16 — atomic Maxwellian Voronov coefficient contract

This contract precedes production changes. It extends existing U16 without
assigning a new ID or changing its relative 1e-12 numerical gate. Linked I23
retains its independent 2% coefficient/3% ensemble-lifetime requirements.

## Inventory and scope

`Moon::ElectronImpactIonizationRate` currently uses a fixed cross section and
400 km/s speed. D05 is archived under the canonical photochemistry package,
not the older illustrative `electron_impact/raw` roadmap path. Native
`cfit.dat` has 406 fixed-width atomic-stage rows; `cfit.txt` supplies positions,
units and native temperature bounds; `cfit.f` supplies the independent fit
implementation. Both processed full and neutral-selected CSVs retain bounds.
Package README, provenance, QA and processing summary have been inspected and
raw/processed hashes must be verified before every data-dependent U16 run.

The source roadmap controls the U16 production-kernel/API/numerical contract;
the data roadmap controls immutable acquisition/provenance. Their differing
hashes and common D05 requirements are recorded with each run rather than
silently reconciling documents. No plasma observation or measured electron
forcing is introduced. The user approved an explicit Maxwellian electron
closure with configured electron density and temperature. Neutral He, Na, Ne,
Ar direct ionization only is covered; ions, molecules, product kinematics and
charge exchange are outside this kernel.

## Production interface and configuration

`srcMoon/species/ElectronImpact.h/.cpp` owns the shared production formula.
`MaxwellianVoronovTable::LoadNativeTable` reads all 406 native rows; it rejects
malformed widths/separators/numbers, duplicates, missing atomic stages,
incorrect labels/stages and invalid physical coefficient/bound values. It
selects neutral He(Z=N=2), Ne(10), Na(11), Ar(18), retaining native bounds.
The pure parser does not guess a data root or implement a second JSON
provenance system: the application/data preflight must hash-bind its selected
file, while U16 performs this binding independently before invoking it.

Temperature is explicitly eV, electron density m^-3, coefficient m^3 s^-1,
rate s^-1. The optional Kelvin converter uses exact SI k_B/e; it is a unit
conversion and adds no electron distribution assumption. Species names are
exact case-sensitive `He`, `Na`, `Ne`, `Ar`. No masses or plasma moments are
invented. Native Tmin is eV; Tmax is keV and is multiplied by exactly 1000.
The four archived neutral fits have inclusive domain [1,20000] eV.

For U=dE/Te the coefficient is
A*(1+P*sqrt(U))/(X+U)*U^K*exp(-U), followed by cm^3→m^3 factor1e-6.
The production evaluator uses the logarithm of this expression so a small
finite rate is not lost through avoidable intermediate overflow/underflow.
It does not impose the legacy Fortran U>80 early-zero cutoff; only natural
floating-point underflow produces zero. The archived neutral native domains
have U<=24.6, so a separately labeled altered-byte arithmetic fixture, never
physical evidence, exercises the deep-underflow path. Domain bounds still
apply; no extrapolation is permitted.

The reduced I23 reservoir preparation also needs the dimensionless derivative
`d ln k / d ln Te` of this same fit. The table returns it for the same four
species and inclusive native temperature domain. It supplies the analytic
slope of the fitted rate. Using it to infer a collision-conditioned mean
requires an explicitly declared unshifted Maxwellian, source/state threshold
crosswalk, and independent model error allowance, as described in
`M3_RESERVOIR_PRODUCTS_CONTRACT.md`. No daughter event or bath is activated by
the slope API. U16 compares the native slope with a 70-digit symmetric
log-temperature difference of the independent direct fit at the existing
temperature grid. The arithmetic check uses the existing `1e-12` scale with
unit scale at a zero crossing. Its outcome is separate from the frozen
coefficient/rate criterion and the I23 ensemble gate.

## Procedure, independent oracle and controls

The thin standalone C++ adapter calls this exact production kernel, compiled
with g++ in the unique test output directory; it contains no ionization law.
The implementation is maintained under `src/models/chemical_database`, with
the historical Moon namespace preserved as an alias. Both shared `DataIO` and
`ElectronImpact` CPP/headers and the common build recipe are source-bound in
the receipt; the empty compatibility CPP is never treated as the kernel.
A separately implemented Python Decimal high-precision evaluator reads native
fixed columns and cross-checks all five coefficient fields against the
archived Fortran DATA blocks and processed neutral CSV. Closed-form density
scaling/zero-density identities supplement numerical coefficient checks.
The grid is declared before evaluation: native endpoints plus logarithmically
spaced internal temperatures over each supported neutral domain.

The original unmodified Fortran is compiled/executed where a Fortran compiler
is available. Its single-precision arithmetic and rounded 8.617385e-5 eV/K
constant are recorded separately; the exact-SI/double production evaluator is
not incorrectly required to agree to 1e-12 with rounded single precision.
This cross-check neither replaces nor relaxes the 1e-12 Decimal formula gate.
Missing optional Fortran capability is explicitly SKIPPED in its supporting
record, not fabricated success. A missing required C++ compiler/data/preflight
means U16 SKIPPED; malformed/hash-mismatched data, build/runtime/parser failure
is ERROR; an evaluated formula discrepancy is FAIL.

Tests exercise uninitialized service, exact neutral mapping, missing/duplicate/
malformed/incorrect-stage records, transactional reload, explicit out-of-domain
and nonfinite/negative input refusal, inclusive endpoint acceptance, zero
and scaled density, Kelvin conversion and arithmetic underflow. Failed calls
leave outputs and an earlier valid table unchanged; the pure layer reports
status/error rather than aborting. No stochastic controls, frame or epoch
enter these coefficients.

## Artifacts and predecessor gates

Outputs belong under a unique source-tree test_output/srcMoon directory, never
canonical raw data. Preserve command/compiler identity and build stdout/stderr,
production/header/adapter/script hashes, both roadmap hashes, raw/metadata/QA/
provenance/processed hashes, native table selection/units/domains, per-point
rate_vs_temperature.csv, Fortran comparison, control-case logs and result.json.
The existing U25 units prerequisite remains declared. Production wiring,
neutral-to-ion allocation, lifetime/event budgets, source drivers and formal
I23 gates remain separate; local U16 PASS cannot close M3 or validate measured
plasma variability. Legacy mode remains unchanged by runtime integration.
