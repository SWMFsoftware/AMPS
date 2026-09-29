# Stage 7: transactional stand-alone shock provider

Stage 7 combines the immutable coronal background, ellipsoid tessellation,
and oblique-MHD jump kernel into a provider-neutral surface snapshot.  The
implementation is deliberately independent of AMPS and MPI.  An application
adapter may distribute the immutable patches, but rank ownership is never
part of a patch's physical identity.

## Patch physics and classification

For each outward-oriented patch normal **n**, the provider computes the
positive upstream inflow in the front frame,

`u_1n^in = V_sh,n - u_1 dot n`,

and the fast Mach number `M_f=u_1n^in/c_f`.  `M_f<=1` is a geometrically valid
candidate front, not a shock.  Such a patch retains geometry and diagnostics
but contributes exactly zero physical source.  A fast patch is passed to the
Stage-6 Rankine--Hugoniot solver; a failed jump invalidates the complete
candidate generation rather than leaving a patch with a partial downstream
state.

The immutable patch records distinguish `geometric`, `fast`,
`supercritical`, `sourceEligibleBeforeClearance`, `sourceActive`, and
`sourceTerminated`.  Termination is history state: a later fast evaluation
cannot infer that the source should reactivate.  `LocateFirstFastCrossing`
uses a bracketed bisection of `M_f(t)-1`, so activation time is independent of
the outer provider cadence.

## One-sided interfaces and transition clearance

`PatchInterfaceEvidence` binds each child patch to its background generation,
interface identity, selected policy, and topology/kinematic/balance results.
Production accepts only evidence that passed either the bounded-approximation
or stationary tangential-discontinuity policy.  Relabeling a diagnostic
kinematic state without its evidence therefore fails before publication.

`SplitShockPatch` replaces an interface-crossing parent by deterministic
one-sided children.  The positive fractions must sum to one.  Area, incident
number rate, and incident kinetic-energy rate use the same fractions, giving
exact closure without rank-dependent renormalization.

Transition-sheet masking is applied after forming the counterfactual eligible
set.  Candidate and excluded area, number rate, and kinetic-energy rate remain
separate.  A zero counterfactual denominator is represented by the shared
typed `InapplicableZeroDenominator` state.  Exceeding any independent cap
rejects preparation; surviving patches are never rescaled.

## Transaction and lifetime contract

`TransactionalShockProvider::Prepare` builds a private candidate, checks every
patch, gate, jump, provenance record, and integral, then publishes one
`shared_ptr<const ShockSurfaceSnapshot>`.  Generation advances only at this
commit point.  Any failure leaves the prior generation untouched, and owning
handles returned by `PreparedSurface()` remain valid across later successful
updates.

The initial gates are intentionally separate from source physics:

- `None` permits a physically valid all-sub-fast, zero-source initialization.
- `AnyFastPatch` requires nonzero fast area.
- `MinimumFastAreaFraction` requires the configured geometric area fraction.

The recommended stand-alone baseline is `None`.

## Verification

`make test-stage7` runs the 117 earlier tests plus `SHK3D05--14` and
`SNAP3D09` (128 cumulative tests).  The new records cover mixed fast/sub-fast
surfaces, root-located activation, all initial gates, conservative one-sided
splitting, persistent termination, owning snapshot lifetime, cadence
independence, transition budgets, interface provenance, and rollback after a
failed update.
