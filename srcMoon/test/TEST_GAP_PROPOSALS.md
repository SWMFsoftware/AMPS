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
- **U07 / I11-I13:** test actual AMR selection and convergence after the
  intended angular/refinement policy is declared.  The current surface callback
  erases its computed subsolar angle and returns a constant.
- **U08 / I14:** convert a declared D01 LOLA product into the production
  triangulated boundary using a committed converter; verify area, normals,
  watertightness, intersections, and longitude/latitude convention.
- **U09 / I15:** wire terrain-horizon calculations to the production
  illumination service; independently ray-cast selected facets and preserve
  real D02 gaps.
- **U10-U11 / I16-I17:** implement and verify simple, Diviner-backed, and
  dynamic thermal states separately, including units, local time, interpolation
  rules, energy balance, phase lag, and restart state.
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
  must be frozen in advance.

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
