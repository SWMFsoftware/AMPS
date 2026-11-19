# Source-effective Na photon / physical-state crosswalk

Preimplementation contract, 2026-11-19. User-approved research establishes a
distinct frozen D04 effective photon threshold and a physical Na ground-state
ionization threshold. Keep both. No source rate, uncertainty, numerical bound,
canonical file, old installed release or legacy configuration is altered.

## Inventory and entry points

`PhysicalReactionEvents::SodiumPhotoionization::Create` currently requires exact
native representation agreement between D04 threshold and formation-state gap.
`ReducedKinematics::PhotonTwoBody` has the corresponding strict default.
`MoonSodiumCohort::Prepared::Create` authenticates a five-column channel/state
CSV, source profile, numerical controls and planner policy. Actual normal
startup/provider/controller/particle transaction already exist. Native cohort
receipts and independent exact-rational/Decimal ledger/numerical verifiers own
the event evidence. The actual source profile remains separate original-research
input; a native source loader grants no managed-database storage permission.

## Additive explicit model

Keep the current `source-equals-state` behavior as the default. Add an explicit
`source-effective-photon` energy policy with a pinned crosswalk. Its eight-column
CSV is:

```text
channel_id,parent_state_id,ion_state_id,electron_state_id,energy_reference_identity,photon_energy_policy,source_threshold_eV,physical_threshold_source_identity
```

Require exactly one complete Na source channel and exact mapped state IDs/common
energy reference. The explicit native eV threshold must equal the actual source
registry value, and the physical-threshold evidence must equal the resolved ion
formation-energy source identity. The caller cannot supply an arbitrary energy
correction. Derive `delta = I_physical - I_source` from authenticated states and
source. Use the frozen representative photon `I_source + mean_excess`; subtract
the physical formation gap in the existing compensated full energy budget.
Thus the rest-frame available translation is `mean_excess - delta` before COM
recoil. A negative/unrepresentable result refuses, never clips.

The same moving-parent photon transform, affine LSO velocity transformation,
actual reference/physical masses, full photon/electron reservoirs, unused local
interval, semantic draws and actual committed-particle readbacks remain required.
The alternative convention changes source-to-state bookkeeping, not the frozen
hazard. Source/state/model errors remain separate from numerical tolerances.

## Verification contract

- ID/name: extend existing U14/U15 supporting event/crosswalk coverage; no new ID.
- Purpose: verify the explicit source-effective/physical-threshold mapping and
  keep strict admission and legacy behavior unchanged.
- Production entries: the kernels and normal owner listed above.
- Configuration: default strict policy; optional authenticated eight-column
  source-effective crosswalk; all numerical/model controls remain explicit.
- Inputs: real D04 quiet/active source endpoints; authentic original research
  Na profile; labelled mathematical fixtures only for local mechanics guards.
- Oracle: independently evaluate source photon energy, physical formation Q,
  two-body COM, Doppler and signed affine reservoirs using Decimal/Fraction.
- Units/frame/time: J, kg, m/s, seconds; one declared inertial frame and actual
  event/geometry epochs; original source eV converted by exact SI definition.
- Randomness: deterministic supplied open quantiles for kernels; existing
  parent-keyed RNG/seed and MPI-independent lineage for linked runs.
- Procedure: strict-default guards first; explicit-policy happy/failure cases;
  quiet/active rest/moving/rotated fixtures; actual native source admission;
  linked cohorts only after source/model-domain admission and fresh producer.
- Metrics/criteria: retain all existing arithmetic/statistical/conservation
  bounds; compare the physical gap and source threshold separately. No relaxed
  threshold guard, invented source interval or new physical-error default.
- Expected PASS: each evaluated kernel/binding/ledger criterion met; this is
  not automatically I21 qualification or observational validation.
- Failure modes: missing crosswalk/evidence; incorrect native threshold or
  mapped source identity; invented correction; negative energy; improper frame;
  stale profile/producer; missing required model-error limits; budget mismatch.
- Required data/hashes: canonical D04 README/provenance/QA and source products;
  original isotope/CODATA/Na threshold/ground-state evidence and prepared profile
  hashes; all recorded with each actual run. Raw data are immutable.
- Predecessors: same-campaign M1/M2 and all formal I21 gates unchanged; approved
  early M3 development allows local/linked supporting verification only.
- Artifacts: exact input/config/source/build identities, compiler and runtime
  logs, JSON receipts, comparison CSV/PNG/EPS, actual event/ledger records and
  retained failures under the dated M3 work package.
- Status: unavailable physical evidence/capability -> SKIPPED; corruption/build/
  parser/runtime -> ERROR; evaluated mismatch -> FAIL; evaluated criterion met
  -> PASS. No local/source PASS replaces a missing formal predecessor.

The source/data roadmap hashes and their different controlling roles are in
the work-package manifest. Scope-specific mechanical error caps require the
review explicitly required by P02; independent implementation proceeds while
that numerical input decision is pending.
