# `srcMoon` species-specific production modules

This directory contains lunar physics that belongs to one simulated species
and should not be embedded in the generic Exosphere library. Generic terrain,
thermal, inventory, mover, and chemistry services remain in their existing
modules; moving them here would incorrectly make shared physics species-local.

The application-independent scalar data services now live in
[`src/models/chemical_database`](../../src/models/chemical_database/README.md).
`Photochemistry.h`, `ElectronImpact.h` and `ChargeExchange.h` preserve the
Moon-facing names through namespace/type aliases; their compatibility `.cpp`
files contain no second copy of the physics. The normal build compiles the
shared implementations once with unique `ChemDb` object names. This preserves
the existing source/domain and rate conventions. Application-owned reaction
execution, surface laws and particle-time/weight accounting remain here.
The installed archive has its own [API and examples](/data/vtenishe/CHEMICAL_DATABASE/README.md);
its catalogue-only profile does not install a physical reaction provider.

`SodiumCase1.{h,cpp}` is the narrowly scoped production adapter for the exact
Tenishev et al. (2013) Case-1 demonstration. It supplies input-validated MIV
rate/temperature and photoionization lifetime to the generic production
kernels, initializes and maintains the prescribed uniform surface reservoir,
and verifies the generated PSD/TD/source switches before a run. It does not
duplicate the PSD distribution, force law, impact handler, chemistry, or
samplers. Inputs without `sodium_exosphere_model=tenishev-2013-case1` do not
install any of its callbacks.

`Argon.cpp` implements the optional continuous radiogenic `40Ar` source. Its
pure parameter validator and uniform reference-sphere area partition are used
by U21. Its production callbacks use the normal AMPS user-defined-source
dispatcher, source accounting, MPI ownership, body/LSO state transform, local
surface-temperature callback, and resolved-LOLA injection adapter. It does not
duplicate a test-only injection law.

The canonical `input/moon.input` keeps `RadiogenicArgon` off because the frozen
M0/M1 executable contains Na only. An experimental Ar run requires all of:

1. `SpeciesList=Na,Ar` at configuration time;
2. `SourceProcessCode=RadiogenicArgon` set to `on` while historical `MySource`
   remains off;
3. a clean Config.pl build with SPICE;
4. `argon_source_enabled=true` in a copied runtime `.in` receipt.

For I18-I20, `srcMoon/test/build_m2_qualification.py` performs the first three
steps in a uniquely named temporary Config.pl application, hashes the inputs
and executable, verifies that generated photolytic chemistry is OFF for the
approved I19 budget, and removes every application-specific generated file.
The aggregate runner then restores the canonical Na-only build. Its linked
runner derives the fourth setting in the preserved per-test `application.in`;
neither source input file is edited in place.

The M2 build also requires individual particle weights ON to represent a final
fractional finite-release packet. The shared Exosphere finite-inventory
implementation owns depletion and departure sampling; Argon retains its
ordinary physical rate and Maxwell-flux source law. During resolved moving
source injection, early spherical MPI ownership is deferred until the final
terrain position at the sampled departure epoch. An intermediate spherical
point can have a different owner; discarding it early would lose valid atoms.
The eight-rank supporting I18/I19 checks exercise this production path, while
formal qualification and temporal/weight convergence remain open.

The runtime parser then accepts only a continuous uniform-global source and
fails if Ar is absent from the compiled species table. Potassium-weighted
geography is deliberately unavailable while D10 is regional and QA-incomplete.
See [`../SURFACE_INTERACTION_LITERATURE.md`](../SURFACE_INTERACTION_LITERATURE.md)
for the rate, residence energy, sticking interpretation, citations, and data
limitations.

`ReactionExecution.h/.cpp` supplies PIC-independent conservative reaction
primitives for the opt-in M3 adapter: exact integer material stoichiometry and
state identities; local physical intervals and compensated optical-depth
continuation; combined constant hazards; immutable represented-number product
plans; explicit signed bath/model-mass/momentum/energy ledgers; and a
reserve/stage/atomic-commit/rollback interface. It does not choose a lunar
recoil/angular distribution or enable any reaction. Transport planning
requires caller-qualified classical kinematics, an explicit event frame/speed
domain and supplied bath exchange; rate-only topology plans cannot mutate
particles. Closure bounds are explicit per-physical-event numerical inputs,
not default empirical corrections. Radiation energy/momentum enters through
an explicit bath ledger rather than a fictitious massive photon species.
The actual PIC adapter, both local-time tracers and the constant-history
controller have separate scoped native checks (828/89/73 PASS in the retained
November17 package). Those analytical mechanics do not qualify physical recoil
providers, variable-hazard exposure or an evolving asynchronous lunar
environment. See
[`../test/reaction_execution/README.md`](../test/reaction_execution/README.md)
for the supporting production-service runner and scope limits.

`ElectronImpact` loads the native D05 Voronov coefficients for He/Na/Ne/Ar
and evaluates an explicitly prescribed Maxwellian bath within the source
temperature range. U16 verifies coefficients/configuration; no measured
plasma distribution or physical ionization executor is inferred.

`Photochemistry` loads the immutable D04 v1 complete channel inventory,
separate quiet/active rates and source state labels. U15 supporting checks
sample all12 parent/solar families with frozen99% statistical bounds. Mean
energies and total rates do not determine daughter velocities.

`ChargeExchange` implements the source-identified He total-capture fit and
Na source-knot interpolation, with explicit units and strict energy bounds.
U17 supporting checks use independently evaluated publication numbers;
physical differential/state-resolved capture and collision-biased bath
sampling remain unavailable. The normal library object is uniquely named
`ChemDbChargeExchange.o` because archive member names share a flat namespace
with the upstream core's unrelated `ChargeExchange.o`. `MoonChargeExchange.o`
now contains only the compatibility translation unit.

The CLI examples and restrictions are in [the Moon README](../README.md).
Source/data contracts and individual fast runners are in
[the rate README](../test/common/m3_rates/README.md),
[the PIC backend README](../test/reaction_execution/NATIVE_BACKEND.md) and
[the controller README](../test/reaction_execution/CONTROLLER.md). Formal
I21-I24 and their M1/M2 predecessors remain open.

The 2026-11-18 independent M3 package adds `PiecewiseHazard` and
`ReactionExposure` to the normal application library. Declared affine rate
histories use relative physical seconds and dimensionless exposure, retain
real support gaps/discontinuities, and locate first crossings without RNG.
The bridge shares the production compensated optical clock and instantaneous
channel selector. It refuses lossy pending-depth projection, zero-rate event
CDFs and unrepresentable rates; it enables no physical provider. The normal
controller remains constant-history-only. See
[`exposure controls`](../test/reaction_execution/exposure_bridge/README.md)
and [`model specification`](../test/M3_PHYSICAL_MODEL_SPECIFICATION.md).

The approved P01–P09 package adds four production component services:

- `ReducedKinematics` resolves additive isotope/electron reference masses,
  preserves separate source formation/excitation terms through cancellation,
  evaluates explicit model-error limits, and implements parent-rest photon
  two-body recoil plus mass-weighted microcanonical three-neutral breakup.
  Directions are normalized after admission, tiny open polar draws retain
  their transverse component, and missing states/thresholds/controls refuse.
- `SodiumPartialCapture` authenticates the twelve D119990 ground-product
  partial-channel nodes. Immutable receipts bind density, source node and
  actual constants. The optional forward-COM reference model uses source-
  derived Q and the shared two-body solver; it requires an explicit target-
  rest frame. Native initial Na state metadata remain unavailable:3s is the
  approved reference cohort, not measured state evidence. This is not the
  all-state Morgan rate or physical lunar CEX activation.
- `IonCollector` records the explicitly selected ideal boundary's incident
  material/charge/momentum/kinetic/formation/excitation transfers. Normal
  sphere/static/moving terrain dispatch reaches its actual PIC weight/epoch
  adapter before neutral sticking. Installation requires a matching source-
  state resolver; the CLI switch alone cannot provide it. No transfer is
  silently interpreted as heat, neutralization or implantation/release.
- `HeliumSilicaChargeSurvival` authenticates seven Tellez2012 TableI fits and
  evaluates the source's explicitly supplied reciprocal-velocity coordinate.
  Mechanism gaps and reported uncertainties are retained. Conversion from
  arbitrary incident/outgoing velocity vectors remains Unavailable because
  the fitted coordinate convention is unresolved. Conditional scattered-ion
  charge survival is not an incident retention probability.

These services are linked in the normal Moon library and tested through
maintained thin runners. Full moving-parent photon providers, trajectory-
dependent chemistry dispatch, complete molecular/cascade data and physical
ion-material laws remain separate unfinished requirements. The generic event
planner's `qualified` flag is never set merely because a model was approved.
See [the approved decisions](../test/M3_DAUGHTER_AND_ION_SURFACE_DECISIONS.md),
[implementation contract](../test/M3_APPROVED_MODELS_IMPLEMENTATION_CONTRACT.md)
and the individual [recoil](../test/reduced_kinematics/README.md),
[Na partial capture](../test/sodium_partial_capture/README.md),
[collector](../test/ion_collector/README.md) and
[He/silica](../test/helium_silica/README.md) test descriptions.

P05 additionally has a reproducible [O I source tool](../test/state_cascades/README.md)
and an offline global supporting runner. The acquired NIST ASD export retains
614 levels, 910 lines, all 56 missing A values and twelve selected/descendant
multipole components across five explicit levels. These are stored source
records, not a certificate that every physical decay exit is included.
Isotope/state applicability, omitted-transition bounds, alignment and quenching
remain unresolved; physical lifetimes and normalized branches stay null and
no cascade is installed in the particle controller. State-source QA and
transport qualification are reported separately.
