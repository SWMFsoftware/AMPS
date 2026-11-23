# Standalone core-model integration: stages 0–4 test contract

This contract implements the bounded scope of
[`INTEGRATION_ROADMAP_STANDALON.md`](../INTEGRATION_ROADMAP_STANDALON.md).
The first physical authority is the existing reduced SSE front/ambient
provider. No particle source, new shock model, acceleration prescription or
downstream-volume closure is introduced by these stages.

## Stage 0: baseline and capability inventory

Capture HEAD, the binary working-tree patch, required untracked source hashes,
input/asset hashes and build settings before implementation. Execute the five
registered prerequisite suites named in the roadmap. Preserve their individual
PASS/FAIL/SKIP/ERROR results. Existing full-volume failures do not become
reduced-front passes. A skipped native row is not native qualification.

The reduced core exposes ambient HCI SI state, branch-preserving Cartesian
derivatives, front trajectory/geometry and immediate one-sided jump states. It
does not expose a physical downstream/sheath/ejecta volume. SWCME and the
neutral coronal-CME bundle retain their maintained owners and capability limits;
they are not silently substituted for missing reduced-model state.

## Planned stable integration cases

The SINT prefix is reserved here for the new stage-1–4 integration cases. Each
case must exercise production adapter/kernel code and emit a stable typed
result. Model-relative expectations test adapter fidelity only; independent
analytical references or separately implemented numerical quadrature must
grade physical/numerical claims.

| ID | Purpose and production scope | Independent/reference criterion |
| --- | --- | --- |
| SINT01 | Core ambient conversion, SI/HCI metadata and generation | Native-neutral values equal direct core values at identical positions/epochs; changed identity rejects |
| SINT02 | Electron/proton/mass density and thermal meanings | Explicit charge-neutral H/He composition, with nonzero alpha abundance, detects electron/proton conflation |
| SINT03 | Complete magnetic and velocity derivatives | Tensor contraction, finite values, branch identity and provider derivative comparison; no density-history gradient proxy |
| SINT04 | Immutable interval ownership and failure atomicity | Invalid time/domain/capability leaves outputs and committed provider unchanged; retained view cannot be rebound by reconfiguration |
| SINT05 | Full-interval trajectory integration | Exact translation/radial/rotating manufactured flow; fixed supplied endpoints and physical dt |
| SINT06 | Actual core Parker material trajectory and state | Independent radial/time and azimuth quadrature; continuity and induction convergence within the declared exterior domain |
| SINT07 | Target-epoch magnetic tracing and geometry initialization | Radial/curved analytic references and actual core B-tangent convergence; either magnetic polarity |
| SINT07R | Advected curved Parker tangent refinement | Fixed 600-s flow interval, 9/17 vertices; tangent residual <2e-6 and ratio >3 |
| SINT08 | Native geometry owner and lifecycle | Actual PIC container, endpoints, empty insertion, whole-chord trimming, particle relocation/loss, exposure and ownership tables |
| SINT09 | Native late-failure transaction | No partial native geometry/physical-state publication after a failed callback or MPI readiness |
| SINT10 | Supported SSE intersections | Independent polar-axis position/speed and off-axis generating-sphere roots; direct core front query grades state fidelity |
| SINT11 | Unsupported rear roots and finite support | Reject rear/cone-exterior roots without inventing shocked state or dropping a numerical failure |
| SINT12 | Endpoint, multiple-hit and tangent ownership | Canonical segment coordinates, duplicate coalescing and explicit tangency classification |
| SINT13 | Query epoch and material-coordinate mapping | Same material worldline at explicit physical epochs; round-trip trajectory reference; committed provider unchanged |
| SINT14 | Forward Parker tube metric and diffusion drift | Physical-number measure yields drift d(kappa)/ds+kappa*dln(A)/ds; exact conditional increments and independent moments |
| SINT15 | Focused deterministic coefficients | Analytical focusing and expansion/strain rates; isotropic angular quadrature recovers Parker cooling |
| SINT16 | Common scattering/transport limits | Equivalent coefficients and independently predicted moments; preserve the distinction between reset-event and pitch-diffusion processes |

| SINT13N | Actual native explicit-time SSE connectivity | 12000-s polar core front, independently computed relative cut speed; immutable metadata/geometry |
| SINT14N | Actual Parker native dispatch | 300-s unchanged step, no bulk displacement, exact cooling and complete Cartesian state commit |
| SINT15N | Actual focused-diffusion native dispatch | Manufactured zero-scattering path against independently integrated gyrotropic ODE |
| SINT16N | Actual focused-event native dispatch | Conditioned no-event path against independent ODE, residual clock and normal native commit |

## Acceptance and evidence

- Use fixed manufactured inputs, documented SI units, seeds and sample sizes.
  Freeze analytical error or statistical-confidence criteria before observing
  outcomes. Numerical accuracy controls do not modify AMPS' physical timestep.
- Register cases through the shared registry and the srcSEP orchestrator with
  per-case result readers, logs and JSON/JUnit evidence. Native and pure scopes
  remain distinguishable. A dependency-adapter pass is not a real MPI pass.
- Run focused tests before dependent changes; preserve failing receipts and
  isolated negative controls. Do not replace expected results with production
  output or change physical gates to obtain a pass.
- Rebuild configured native objects after class/API changes. Record source,
  generated configuration and executable hashes. Compare the retained baseline
  and final source so unrelated local changes are preserved.
- Stage completion requires the roadmap's actual production path and required
  tests. A bounded kernel PASS alone cannot complete stages 2 or 4. Missing
  capabilities and any unexecuted required case remain explicit.

Stage 5 source choices and downstream-boundary physics are separate. Tests may
manufacture a source or initial distribution for verification, but must label
that scope and cannot make it a production SEP injection prescription.
# Stages 5–6 source and absorbing-boundary gates

These additional gates are fixed before observing implementation results.
Manufactured distances and source strengths in fixtures are verification
inputs; they are not production defaults or calibrated SEP parameters.

- SINT17: reduced-source options preserve disabled operation, require explicit
  positive release distance when enabled, and reject pass-through/gross modes
  without mutating the supplied configuration.
- SINT18: local-compression DSA number spectrum agrees with the independent
  analytic integral of p^(2-3r/(r-1)); no phase-space/number Jacobian ambiguity.
- SINT19: release is an intersection of the same material line with the
  sphere's positive normal-distance offset, using the actual birth epoch;
  field-line and front-support residuals have absolute SI bounds.
- SINT20: normal incident proton flux times tube/surface area and explicit
  net release fraction gives the dimensional source rate. Electron density
  and whole-front area substitutions must be detected.
- SINT21: focused escaping pitch distribution has positive relative normal
  flux and agrees with the independent integral of the linear positive part.
- SINT22: conditional Brownian first passage retains positive-endpoint hidden
  hits with probability exp(-d0*d1/(kappa*dt)). For 65,536 fixed-key samples,
  counts must lie within six binomial standard deviations of that independent
  value. This is an absorbing-boundary kernel gate, not native qualification.
- SINT23: unconditional Brownian first-hit probability and early-hit CDF agree
  with erfc(d/sqrt(4*kappa*t)), also within six binomial standard deviations
  for 65,536 fixed-key samples. Deterministic and invalid inputs consume no
  random draws; invalid arguments return failure. A second independent ensemble
  grades a linearly moving front at 1 m/s with d0=2 m, kappa=1 m^2/s and T=1 s.
  Its independent first-passage CDF is
  `0.5*erfc((2-t)/sqrt(4*t))+0.5*exp(2)*erfc((2+t)/sqrt(4*t))`;
  total and early (t=0.25 s) counts retain the same N=65,536 and six-sigma rule.
- SINT24: a between-event MFP shock hit stops at the independent ballistic
  first-hit time 0.1 s; no branch, pitch, or next-clock variate is drawn after
  absorption. Event timing and deterministic interpolation remain separate.
- SINT26: exact rotation-axis background evaluation agrees with independent
  axial-dipole/radial-wind relations to off-axis core queries at equal radius.
  The existing strict scalar bound is unchanged. The off-axis evaluation
  retains the complete Parker-winding quadrature; the optimization may skip
  only an exactly identity rotation on the Cartesian axis, without skipping
  magnetic-null, finite-value or coverage checks.
- SINT27: the core's exactly parallel, normal-incidence RH bracket shortcut
  agrees with independent gas-dynamic compression, pressure and velocity
  identities over both gamma values, several sonic Mach numbers and both
  magnetic signs, including B=0. Transverse-flow controls use the ordinary
  scan. Every case retains the production residual bound and entropy/branch
  checks; nonzero transverse B and near-identity weak roots retain the scan.
- SINT28: 4096 fixed-key non-dyadic internal MFP durations with one early
  scatter must report the exact closed endpoint, with no epsilon. Key 20261122
  fixes the cases. Absorbing controls at fractions zero, one and interior
  values must retain the closed interval and consume no later random draws.
  These internal durations do not change AMPS physical dt.
- SINT29: 64 forward/reverse constant-flow trajectories at heliocentric
  coordinates exercise short final numerical intervals. Their independent
  endpoints must remain within the unchanged 1-mm position tolerance and
  128-attempt work limit. This positive stress case does not reproduce or
  resolve the separate full-application iteration-84 trajectory failure.
- SINT30: exact trajectory inputs retained from the identical one-/four-rank
  iteration-84 failure must complete with the original 8.503654-mm position
  tolerance, 37.5-s maximum numerical step and 100000-attempt work limit.
  An independent fixed-step midpoint oracle accumulates in long double on
  8192/16384/32768-point grids; both reference refinements must differ by at
  most tolerance/64 before grading the production endpoint. The retained old
  controller fails this case. No physical AMPS timestep is changed.
- SINT31: both polar hemispheres, all signed-zero azimuthal inputs, three radii
  and three epochs compare the ambient against uncached canonical PFSS with
  explicitly nonaxisymmetric harmonics. Exact equality is required. A separate
  magnetic-null control must still reject; caching cannot bypass admission.
- SINT32: a manufactured quadratic front has exact linear-chord error h^2/4.
  Cold and hinted production refiners must accept the same 128 intervals of
  0.0625 s at the unchanged 1-mm bound and finish at 8 s. The warm controller
  must use at most 200 checks (149 observed), versus the fixed 1024 cold
  checks. A failed 64-attempt refinement must publish no interval or history.
- SINT25N/SINT26N/SINT27N: actual net-source births, explicit remaining duration,
  unchanged 300-s physical steps, and joined source/mover/geometry/energy
  accounting through Parker, Dmumu and MFP respectively. Incident number and
  birth angular laws use independent polar-axis formulas, with relative number
  bound 2e-10 and pitch bound 1e-10. Ledger residual bounds are 1e-10 relative;
  these are floating-point accounting limits, not observational accuracy.
  A repeated exact boundary query must return identical coordinates/duration.
  The full advected-polyline oracle still grades the polar specialization to
  1 mm. The retained full-range binary search additionally requires exactly
  equal root coordinates after the search-hint optimization. Step two changes
  the immutable material-view owner after motion and
  prefix insertion; old arc-label cache entries must not survive that change.
  The actual output writer must identify the reduced SSE/Parker authority,
  report its core heliocentric apex at the requested epoch, and serialize
  current/previous magnetic fields in the advertised nanotesla units. Native
  stored fields remain in tesla; only text serialization converts their units.
  The one-rank SINT26N fixture additionally relocates one fully transported
  record into the outgoing tail after sampling, solely to grade geometry
  retirement and its immutable source birth-energy accounting. It must remove
  exactly one computational particle. That controlled relocation is not a
  physical trajectory and does not qualify a source spatial distribution.
- SINT28N: explicitly labelled manufactured monoenergetic source.
- SINT29N: zero efficiency yields zero source particles without changing dt.
- SINT30N: an inadmissible finite-reference geometry rejects before allocating
  particles or committing source counters.
- SINT31N: a physically valid material line outside the finite SSE cap remains
  disconnected and receives zero shock particles.
- SINT32N: a net source with an incompatible 1–2 eV momentum support rejects
  before allocating particles or committing source counters. Silently reducing
  its specified net rate by excluding slow proposals must fail this gate.

Native fixtures perform direct mover dispatch, join all destination lists and
execute the existing field-line MPI exchange before closing a birth window.
Particles can stream across rank-owned segments; joining only the original
owner's destinations leaves outbound temporary lists unjoined, while joining
without migration strands particles outside the next owned-segment sweep.
They do not impersonate PIC::TimeStep. The final gate must also execute the
linked complete application with the same explicit source/return configuration
and require completed clock, sampling, geometry and signed ledgers. Missing
output, MPI startup failures and unimplemented modes remain errors.
