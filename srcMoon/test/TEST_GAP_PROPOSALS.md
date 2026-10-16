# `srcMoon` test-gap proposals

These are coverage gaps found during the source audit.  They do not allocate
new permanent test IDs.  Each item should be reconciled with the test registry
before implementation.

## Legitimate extensions of existing IDs

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
- **I13 reference review (I11/I12 implementation complete):** I11 now runs the
  full 40/20/10/5 km production campaign and I12 verifies the linked realized
  mesh. I13 cannot PASS until a passing I11 run produces its candidate science
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
