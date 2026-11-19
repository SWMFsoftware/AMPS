# `srcMoon` test-gap proposals

These are coverage gaps found during the source audit.  They do not allocate
new permanent test IDs.  Each item should be reconciled with the test registry
before implementation.

## November 13 approved implementation checkpoint

The user approved v3. Existing U12/I19/I20 contracts and source roadmap are
extended consistently; no new permanent ID is assigned. Active numerical
controls, passive cohorts and the whole-command/cluster execution infrastructure
are described in [M2_NUMERICAL_PROTOCOL.md](M2_NUMERICAL_PROTOCOL.md).
The review history below retains the original gaps at their recorded dates;
statements there that the design is unapproved describe those earlier dates.
The current implementation must be verified by real native execution. Formal
convergence, measured cluster feasibility, native D02 temporal/operator evidence,
independent thermal/PSR comparisons and M1/I13 closure remain open.

## Legitimate extensions of existing IDs

### November 12 follow-up: D02 operator and formal I19/I20 protocol

Primary-source research is recorded in D02_ILLUMINATION_RESEARCH.md. The native
AVGVISIB observable is the fraction of sampled times with any solar-disc
visibility, not mean visible area or flux. Neither general PGDA descriptions,
older six-hour calculations nor native date fields establish the complete
201608 averaging contract. Exact epochs/sampling and an independent
time-resolved I15 benchmark remain missing. No synthetic epoch labels can be
derived from the averaged maps.

The spatial processor's unsupported exact-product timing claims are removed
from QA, with all 366 selected CSV rows byte-identical. Overall QA remains
PARTIAL while spatial processing is READY. Ten automatically registered
metadata protocol checks protect this distinction; they extend the existing
M2 data-audit entry and introduce no scientific test ID or validation PASS.

M2_I19_I20_PROTOCOL_PROPOSAL.md/JSON propose synchronized physical-time
accounting, independent dt/weight refinement, fixed independent replicates,
simultaneous equivalence budgets and censored reversible-retention diagnostics.
These remain unapproved, unconsumed proposals. Full source/mover coverage,
physical-time diagnostics, formal convergence and I20 qualification remain
open; the short supporting ledger cannot close them. Review the protocol and
update the existing I19/I20 contracts/roadmap consistently before implementation
or formal execution. No new permanent ID is proposed.

The subsequent v2 review identifies the Ar timestep copy (0.1 times Na after
initial messages), the supporting fixture's hard-coded seed, and missing
physical-time/cohort diagnostics. It defines two separate physical cases,
288 potential full-duration runs, first-release versus re-adsorption metrics,
an explicitly conditional statistical approximation, and bounded capacity
planning. These proposals are not implementation or approval. The complete
clock/seed/tag/source-path contracts must be verified before formal execution;
current 20-step checks do not establish month-long convergence or precision.

### November 12 execution-budget follow-up

V3 records the user-required one-hour aggregate development campaign and
deferred cluster qualification of order a few hundred CPU-core hours. The
current runner resets its child timeout for each linked test and leaves
build/QA/plot/cleanup outside that guard. Whole-command deadline propagation
and safe process-group termination remain infrastructure gaps, not completed
features. The proposed short 24/48/96-step, two-seed cases provide scoped
production/conservation diagnostics only; long convergence and full seed
qualification remain separately SKIPPED until executed. A candidate 288-run
cluster matrix must be costed against an aggregate allocation before launch,
not granted an unbounded per-job allowance. No acceptance gate is relaxed.

### November 12 M2 supporting coverage

The approved development-order exception allows implementation/supporting
checks while I11/I13 remain unresolved; formal I14-I20 retain their gates.
No new permanent IDs or relaxed thresholds were introduced.

- U09 now reproduces and fixes shared-vertex misses on a closed icosphere.
  Supporting I14 verifies 513/513 LOLA radial rays and collisions on all seven
  eligible ranks in an eight-rank run. The actual selected force mover now
  enters the moving-surface callback, and source position/flight share a
  departure epoch. Mesh/time/collision convergence remains a formal gate.
- U11 adds the analytical uniform radiative-equilibrium invariant at three
  step sizes. This does not replace I17 spin-up or transient convergence.
- U12 checks integrated competing finite-release hazards, including stiff and
  semigroup limits. Supporting I19 closes Na/Ar gas and evolving finite stocks
  over 20 ordinary steps, with actual adsorption/release/escape and weighted
  packets. This closes that scoped production/accounting gap, while formal
  I19, other source/mover/chemistry layouts and temporal/weight convergence
  remain open. Population-independent sputtering is explicitly unsupported
  with finite stocks until its depletion coupling is specified.
- D03 selected-input package qualification is now reproducible and READY;
  D02 remains PARTIAL. Actual D02 integration epochs and the extended-disc
  temporal comparison are missing. I16/I17/I20 operator, hold-out and
  predecessor requirements remain; data readiness is never validation.

The older audit bullets below describe initial gaps. This dated entry records
the precise subset addressed; it does not declare their complete scopes closed.

- **U01 / I01:** expose a thin adapter around the same production ray/sphere
  intersection and LOS integration routines used by sampling; add tangent,
  miss, inside-origin, and reversed-direction analytical fixtures.
- **U03 / I05:** add two-body trajectory energy/angular-momentum invariants and
  timestep convergence after the local acceleration gate.
- **I06 (U04 local prerequisite implemented):** U04 now freezes the
  J2000/LSO/MOON_ME_DE421 convention, epoch, exact DE421 kernel hashes, and
  independently checks the production differential/fictitious kernels and
  frame derivative. I06 must still propagate a particle through the actual
  linked AMPS mover and compare the complete trajectory; U04 is not a
  substitute for that wiring test.
- **I07 (U05 local prerequisite implemented):** U05 now qualifies 14
  unobscured points of the absolute one-AU Na radiation-pressure curve against
  a reproducible 150/300 dpi digitization of Combi et al. (1997), Figure 7,
  with a predeclared `0.5 cm s^-2` graphical uncertainty. The original
  machine-readable digitization remains unavailable. I07 must still propagate
  particles through the actual linked mover and verify force application and
  suppression through lunar and terrestrial shadows; U05 is not a substitute
  for that production-wiring test.
- **U06 / I08:** enable a declared ion species, mass, charge, E/B fixture, and
  compare the production mover with the analytic gyro-orbit.
- **I11 statistical-support audit (2026-11-11):** `I11Sampling::Analyze`
  assigns zero relative CI width when both the five-batch mean and its
  empirical variance are zero. A deterministic structural zero and a sparse
  Monte Carlo bin with no sampled events need different evidence. Extend I11
  with boundary-case statistics fixtures and an explicitly reviewed support
  policy before relying on this shortcut for qualification; do not invent
  event counts or a probability law for weighted/correlated samples. No
  threshold or current run is changed by this proposal. The corrected live
  480/240 km runs each contain 1175 metrics and have no all-zero nine-window
  series; their observed first failure is independence, not this potential
  edge case. No new permanent ID is assigned.
- **I13 reference review (I11/I12 implementation complete):** I11 now runs the
  approved 480/240/120/60 km campaign (explicit legacy FULL remains
  40/20/10/5 km), and I12 verifies the linked realized mesh. I13 cannot PASS
  until a passing I11 run produces its candidate science
  CSV and deterministic configuration JSON, and an investigator reviews and
  freezes them as `I13_archived_baseline_metrics.csv` and
  `I13_archived_configuration.json`. Do not automate that scientific approval
  step or treat candidate generation as a baseline PASS.
- **U07 / I12 QUICK timing qualification:** the first FULL I12 attempt spent
  106.5 minutes building the AMR tree and then exhausted `/data` by writing a
  739 GB optional `mesh.dat`; it is correctly classified `ERROR`, not a mesh
  verification result. Full-volume Tecplot output is now input-gated and off
  in committed receipts. A factor-two QUICK attempt reached 1,074,788 blocks
  and was stopped after 28.7 minutes at 223.7 GiB RSS on eight ranks. The
  revised QUICK receipt uses three times the FULL targets to cross the next
  coarser dyadic level. Its cold eight-rank wiring run completed in 1,409.9
  seconds (23.50 minutes), with 191.9 GiB peak RSS and a 6,128,421,120-byte
  center-measure cache. It was correctly `SKIPPED / NOT QUALIFIED`; even a
  successful QUICK wiring check cannot replace FULL I12.
- **U08 / I14:** convert a declared D01 LOLA product into the production
  triangulated boundary using a committed converter; verify area, normals,
  watertightness, intersections, and longitude/latitude convention. The
  2026-10-11 linked smoke run also exposed a production-wiring gap: the
  first-order no-force mover's initial internal-boundary switch rejects the
  registered Nastran descriptor at iteration 1. Repair and exercise that exact
  mover path before I14; a local geometry-kernel PASS is not a substitute.
- **U09 / I15:** wire terrain-horizon calculations to the production
  illumination service; independently ray-cast selected facets and preserve
  real D02 gaps.
- **U10-U11 / I16-I17:** local processed-Diviner interpolation, the published
  C3--C9 bolometric operator, grey-emission conversion, and the depth-resolved
  dynamic thermal kernel are implemented and locally verified.  Remaining
  coverage is the linked observation-free periodic spin-up, restart continuity,
  production phase/amplitude and time-step convergence, D03 qualification,
  June-only threshold freeze, and one-pass December hold-out score.  Until
  those gates exist I17 remains `SKIPPED / NOT VALIDATED`.
- **U12 / I18-I19:** add seeded Bernoulli and Maxwellian/accommodation tests;
  verify adsorption/desorption bookkeeping, non-negativity, conservation,
  residence-time laws, and checkpoint continuity.  The present surface-density
  configuration overwrites the inventory with a prescribed constant.
- **U13 / I20:** define the PSR mask, cold-trapping criterion, release law, and
  coupling to production illumination/thermal services.
- **U14-U15 / I21-I22:** generate species-specific rates from qualified D04
  inputs and verify shadow gating, branching, daughter creation, and survival
  statistics without using AMPS output as the oracle.
  The November 12 supporting I11 fixture now checks the actual configured
  neutral-Na processor against exponential sunlit survival and exact shadow
  survival. It exposed and corrected a generic H/O-only processor bypass of
  the Moon Na lifetime callback. This closes that specific dispatch gap;
  it does not qualify the legacy coefficient against D04, daughter kinetics,
  multispecies branching, or the complete I21/I22 scope. The neutral-Na sphere
  event ledger closes actual gas events, but complete I19 surface-reservoir,
  evolving-geometry and multispecies conservation remain open.
- **U16 / I23:** replace or qualify the fixed cross-section/speed constants
  against D05 and test rate integration over a declared electron distribution.
- **U17 / I24:** declare reactants/products and authoritative charge-exchange
  cross sections, then test production stoichiometry and momentum/weight
  bookkeeping.
- **U18 / I25-I26:** repair and qualify D06 before testing interpolation,
  frames, gaps, and production E/B/plasma driver wiring.
- **U19-U21 / I27-I28:** declare the He reservoir, Ne accommodation law, and Ar
  geography/transient equations and source references before enabling their
  production dispatch.
- **U22 / I29:** detect duplicate source registration, test each enabled Na
  generator statistically, verify surface depletion policy, and close the
  injected-particle/source-rate budget.
- **U23 / I30:** wire qualified D11 stream/sporadic forcing to the production
  impact source and test timing, units, and mass-to-vapor conversion.
- **U24 / I31:** add H2O/OH species, reactions, surface migration, branching,
  and conserved H/O accounting before observational use.
- **U25 / I41-I46:** add a package selector that reads each required README,
  `provenance.json`, and `qa_report.json`, verifies every referenced hash, and
  emits an immutable run manifest with build, executable, layout, seeds, and
  external-data identities.
- **U26 / I32-I40:** implement independent geometry fixtures for Kaguya TVIS,
  LADEE NMS, LAMP, LACE, PACE, and ground-based Na operators before any
  observational score.  Instrument selection, calibration, and hold-out sets
  must be frozen in advance.  For I33, preserve native per-pointing UTC and
  J2000 RA/Dec (or the documented ECLIPJ2000 transform), implement physical
  exposure averaging, and do not select a WHAM Gaussian/top-hat beam law until
  an authoritative radial response is available.  D14 currently supplies only
  four exact brightest-beam positions plus a circular diameter and exposure,
  not the numerical raster, per-beam timestamps, spectra, or beam response.

## Cross-cutting gaps

### Asynchronous local-time sampling — 2026-11-18

The new force/event tracers preserve physical intervals and represented number,
but this alone does not qualify normal density/LOS/brightness/surface-flux
sampling with different particle endpoints. The canonical Na-only input uses
SpeciesGlobalTimeStep and SpeciesGlobalParticleWeight. Its snapshot regression
is a different scope. `PIC::Sampling` samples before injection/move; `ProcessCell`
adds number and number/volume once, and mesh output divides by sample count.
Moon's direct LOS accumulator/output likewise uses number/volume and iteration
counts, not stored particle epochs or residence exposure. This is a concrete
missing observation/sampling capability for asynchronous science, not a failure
of the completed mechanical number/continuation controls.

Declare common-time transient snapshots versus stationary residence-time
estimation before implementing either. Extend the existing I11-I13 numerical
and appropriate observation-operator contracts with an independent steady-beam
number-density/residence oracle and controlled common-time snapshots as
appropriate; retain the existing I21-I24 event/exposure gates separately. Those
extensions need a consistent reviewed registry/contract, not a new permanent
ID assigned here. No unimplemented estimator or asynchronous density product
is called verified/validated. Data/kinetic and formal predecessor gaps remain.

### M3 chemistry event/product audit — 2026-11-17

The consolidated [M3 completion proposal](M3_COMPLETION_PROPOSAL.md) adds a
source-backed implementation/test sequence, upstream API reuse limits, explicit
physics/data decisions and development/cluster execution design. It proposes
extensions under the existing U14–U17/U24 and I21–I24 meanings, with no new IDs
or changed criteria. U24's M3/M4 overlap, I23's measured-forcing claim and the
roadmap's nonexistent existing U17 channel require reviewed clarification;
the proposal does not resolve them by silently changing this registry.

The initial M3 audit identified uncovered production seams in the existing
generic electron and photolytic executors: independently sampled product
birth/parent-loss packets, unequal species dt/weight transfer, a shadowed
electron product-loop cell index, and collapsed O(3P/1D/1S) identities. A
correct coefficient table cannot verify event-level nucleus/charge/momentum
closure or physical daughter transport. Molecular archived mean energies do
not supply joint distributions; the experimental OH source also needs
authoritative wavelength/threshold-tail interpretation.

Extend the EXISTING I21-I24 contracts with actual conservative transfer,
same-cell insertion and unequal-weight/time-step cases when those production
capabilities are implemented. Extend U15/U24 with complete channel/state
identity and required-data refusal. Preserve every original numeric criterion;
these are legitimate existing-test extensions, not new permanent IDs. The
current U16 verifies coefficients and rate/configuration wiring only. No
generic executor is activated or fixed by silently substituting a test-only
law. Details and proposed physics choices are in `M3_INITIAL_IMPLEMENTATION_CONTRACT.md`
and the two `M3_*_PROPOSAL.md` research documents.

- Add deterministic random-seed control and statistical confidence intervals
  to every stochastic source, surface, chemistry, and mover test.
- Add MPI-rank/OpenMP-thread invariants and statistical equivalence rather than
  assuming bitwise identity.
- Exercise production checkpoint/restart for particles, surface inventory,
  time-dependent thermal fields, sampling accumulators, and RNG state.
- Record compiler, MPI implementation, compile macros, generated-source hash,
  source revision, executable hash, and runtime layout for linked runs.
- Audit diagnostic labels and output units; several source comments and legacy
  names refer to Mercury and cannot be treated as lunar scientific metadata.

### Approved common-time follow-up — 2026-11-18

The user approved common physical time as the initial qualification direction.
`MoonTransportClock` now supplies guarded equal-step chronology and actual
source/mover Begin/Commit wiring, with explicit frozen external epoch and
advancing transport births. Local controls and real one/eight-rank hook checks
are separate from full eight-rank normal application startup. These extend
existing U15/U16 supporting coverage and do not assign a permanent ID. The
asynchronous estimator gap above remains; common chronology does not establish
numerical convergence, a source model or an observation operator.

Declared affine exposure integration and the optical-clock bridge have local
analytical controls. A physical trajectory-dependent provider/controller is
still missing. A compensated pending depth not exactly representable by the
current scalar exposure interface is explicitly refused. The full physical
model and qualification specifications preserve unavailable distributions and
source-interpretation gaps rather than filling them with chosen kinematics.
See `M3_PHYSICAL_MODEL_SPECIFICATION.md` and `M3_QUALIFICATION_PROTOCOL.md`.

### Approved P01–P09 component follow-up — 2026-11-18

The user approved explicit reduced photon/neutral-three-body models, a bounded
Na partial-capture reference experiment, an ideal verification collector and
the He/silica conditional laboratory scope. Their supporting runners execute
under the existing global component namespace. This adds no permanent ID or
lunar material-validation gate and changes no frozen acceptance threshold.
U15/U17 remain incomplete without full physical products/predecessors; literal
U24 still includes actual water/OH surface migration and budgets.

Ion-boundary component coverage now distinguishes strict refusal from an
explicit mathematical collector. A complete lunar ion-material test still
needs independently sourced joint reflection/neutralization/implantation,
target secondary yields and later-release laws. Existing neutral Na/Ar U12
does not establish those laws. Assigning any new permanent material-test ID
requires a consistent roadmap/registry review; none is assigned here.

The D119990 native label specifies outgoing H(1s), Na+(2p6), target frame and
units, but not the initial Na electronic state or a full uncertainty/source
crosswalk. Na3s is the approved reference input, not an inferred measurement.
The He/silica fitted inverse-velocity projection also remains unresolved.
Neither gap can be closed by the scalar law or a conservative sampler.
Complete state/threshold/isotope/energy error bounds, physical photon/bath
providers, cascades and actual trajectory-dependent chemistry remain required.

The acquired O I ASD packet now retains the complete stored level/line exports,
selected O(1D)/O(1S) exits and 3P fine-structure descendants. This resolves
source acquisition/transcription coverage, not the physical-completeness gap:
all stored records are not all possible exits. Applicable isotope/state
mapping, alignment/quenching and quantitative omitted-process bounds remain
missing. The processor therefore retains missing A values and leaves physical
lifetimes/normalized branches null. Fifteen offline source controls are named
U15 supporting children; no full cascade PASS or permanent new ID is assigned.

### Source-bound completion follow-up — 2026-11-18/19

The earlier implementation gaps above now have supporting production services:
`PhysicalStateProfile`, `PhysicalReactionEvents`, `TrajectoryExposure`, normal
explicit-ET geometry/controller wiring, exact compensated native clocks and
full semantic daughter RNG records. The source-bound experimental Na composer
retains moving-parent photon and electron reservoirs; the actual conservative
backend validates initial PIC child readbacks. A named conditional I21 cohort
owner and independent ledger/first-hit/statistical auditors are installed.
These extend existing U14–U17 supporting coverage without a new permanent ID.

Their implementation does not close the physical-data gaps. Current profiles
retain fourteen unavailable state records, including Na+; D04 and evaluated
Na ground thresholds differ beyond native representation. No source/Q interval
is inferred from printed digits or widened through a numerical tolerance.
Molecular initial/product populations and Q, collision-biased electron joint
kinetics, complete CEX families/differential states, applicable cascades and
physical lunar ion-surface outcomes remain unresolved. Unsupported channels
stay disabled and full I21–I24 remain SKIPPED/NOT QUALIFIED. Seven available
reference states and scalar rates are not a complete reaction family.

The normal full-force cohort also requires the actual normal environment-stage
refresh before its direct production walker; source-only/zero-force admission
must not disguise that separate wiring requirement. Its extraction and native
proof are tracked explicitly in STATUS.md. The full qualification matrix,
same-campaign predecessors,0.2/0.3percent numerical allocations, exact branch
confidence and independent cohort membership remain unchanged. No supporting
fixture, regression output or downloaded table is promoted to validation.
