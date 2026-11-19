# Source-explicit heavy products and reservoir moments — preimplementation contract

2026-11-19; source basis `8fbbd5e`. This extends the user's approved
parent-velocity/reservoir scope; it does not claim a full collision-conditioned
joint daughter law, silently choose states or activate any application provider.

## Inventory and production boundary

`src/models/chemical_database/ReducedModels.h` supplies transactional explicit-SI
`InheritParentVelocity` and `TwoBodyProducts`. These mathematical services select
no source, distribution, species or profile. `srcMoon/species/ReactionExecution`
adds integer nuclear/charge topology, actual particle intervals/weights and
plan/reserve/commit. `ReducedKinematics` resolves explicit states/formation Q;
`PhysicalReactionEvents` currently binds a source-specific photon provider.
Those application classes must remain outside the new shared module dependency.
Existing Voronov total rates do not by themselves provide individual electron
kinematics. The shared API must remain C++17/PIC/MPI/RNG independent.

New `AMPS::ChemicalDatabase::ReservoirProducts` consumes explicit source/state
profiles, kinetic and physical masses, integer isotope nuclei/charge, state
energies on one reference, actual parent velocity/position/time/represented N,
and explicit numerical/model error budgets. Typed `Runtime::Status` reports
software evaluation only. Missing source/Q/crosswalk/error proof is unavailable;
invalid contradictions are errors/domain refusals. Every failure preserves the
output object. No installation, default reaction set or physical activation.

## Proposed admitted reduced models

* Electron impact: `e+A -> A+ +2e`, one incoming prescribed electron and TWO
  outgoing electrons represented by aggregate number/momentum/kinetic/state
  energy/charge moments. The ion inherits the actual parent velocity. Incoming
  and outgoing moments and their signed difference are retained separately;
  no individual outgoing energy/angle PDF is inferred.
* For a T-independent cross section and unshifted isotropic Maxwellian around
  the resting target, differentiation of `k(T)*T^(3/2)` gives collision-conditioned
  mean incoming energy `kBT*(3/2+d ln k/d ln T)`. The caller supplies the exact
  source/fit derivative, domain and threshold; the helper does not duplicate
  Voronov formulas or call a fit-derived moment a measurement. Unsupported drift,
  missing derivative, outside temperature domain or mean below threshold refuses.
* Forward partial CEX: `A+B+ -> A+ +B`, source/state-specific partial only. Use
  actual incoming velocities and source-derived Q in the shared conservative
  two-body service, with heavy-product COM direction parallel to its incoming
  COM direction. The removed prescribed plasma ion and produced fast neutral
  are distinct reservoir terms. A bare electron sink is insufficient. Zero
  incoming relative velocity/direction, impossible energies or absent partial
  source/state/error profile refuses. No all-state rate/branch renormalization.

## Quantitative limited-observable bounds

For the electron model, with parent-rest conditional incident electron mean
energy Ein and zero-ion-recoil outgoing aggregate budget Eout=Ein+Q, conservation
and Cauchy/Minkowski give an upper RMS omitted ion impulse
`sqrt(2*me*Ein)+sqrt(4*me*Eout)`. Divide by ion kinetic mass for an RMS speed bound.
A corresponding mean absolute lab-frame kinetic-energy discrepancy is bounded
by `M*|v_parent|*dv_rms + M*dv_rms^2/2`. This is an ensemble moment bound, not an
all-events support bound or a fitted confidence interval. Formation/mass
conventions and source assumptions remain explicit. No downstream orbit/LOS
error is inferred from this local bound.

For forward CEX, no measured angle is invented. With explicit source-supported
maximum COM deflection theta, the heavy-product speed discrepancy relative to
forward is `2*p_relative/m_heavy*sin(theta/2)`. Its lab kinetic discrepancy bound
follows the same triangle identity. A full pi bound is permissible only as an
explicit mathematical all-angle bound; it is not a measured narrow-angle claim.
Source angle support/error evidence and caller observable limits are required.
Physical-versus-kinetic mass error is reported separately, never hidden in a
numerical conservation tolerance. Caller limits are not relaxed on failure.

## Verification contract (support under existing U16/U17/I23/I24 meanings)

Purpose: first-principles software/moment/conservation checks, not observations.
Independent oracles: integer isotope/charge identities; analytic Maxwellian
constant-sigma moment2kBT; source-supplied power-law derivative identities;
independent vector/moment/energy sums and forward two-body expressions; closed
RMS and angle bounds. Administrative nominal-mass fixtures are NOT measured
Na/CEX or lunar profiles. Units: SI kg,m,s,J,C; input thermal energy J. One named
inertial frame/time scale per event; actual physical interval and represented
number are copied, without product-dt resampling or RNG. Every event records
separate incoming/outgoing bath number/P/E/charge/nuclei, net bath exchange,
state Q, numerical residuals, source/profile/crosswalk identities and bound scope.

Procedure: compile a thin C++ probe against the actual shared production header;
check successful EI/CEX cases and exact topology/weight/time preservation; rotate
and boost analytical cases; vary masses/thresholds/temperature; check missing
proof/drift/domain/impossible-energy/zero-direction/overflow and transactional
output guards. Python independently evaluates expected moment algebra; no
production output is frozen as physical truth. Native source/build/executable
hashes, commands, stdout/stderr and `result.json` retained in the source-tree
campaign output. Missing compiler/source is SKIPPED, numerical/parser/build
fault ERROR, evaluated analytic criterion miss FAIL, met criterion PASS. No new
permanent U/I ID. Formal I23/I24 and earlier predecessor gates remain unchanged.

## Source protocol and remaining physical gaps

The bounded source partial currently documented for supporting Na CEX is IAEA
D119990 at its native4–25keV/u nodes. The initial-Na-state source crosswalk and
angular discrepancy evidence remain explicit prerequisites; the shared helper
contains neither Na constants nor a selected node. For Maxwellian EI, an analytic
source-rate derivative defines a MODEL moment, not an independent electron
measurement. The full outgoing joint PDF is unnecessary for the approved reduced
heavy-ion observable scope; a source-consistent profile, state/Q mapping and
quantified local model error remain necessary. No mean electron spectrum, fixed
recoil energy or default isotropy is manufactured.

Source roadmap SHA256:712e1152bf50af6d9ee1dad8f8a585809ca8b92cf42ac6fc7a009f6c0253ee17.
Data roadmap SHA256:d10eca1c539f9cf82940236d9e1aa36af665d43e82ba475925950d582fbf49c4.
Source amendments control approved development order; unchanged scientific,
source interpretation and formal predecessor rules control qualification.
