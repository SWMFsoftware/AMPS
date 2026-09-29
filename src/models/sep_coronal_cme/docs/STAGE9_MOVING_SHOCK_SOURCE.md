# Stage 9: moving-shock source and deterministic particle allocation

Stage 9 converts a prepared Stage-7 shock snapshot into a physically
normalized upstream first-passage source.  Physical number and energy rates
are computed before selecting a Monte Carlo sample count.  Increasing samples
therefore changes estimator variance, not the physical source strength.

## Conditional momentum laws

`MomentumSpectrum` represents normalized `dN/dp`.  A fixed comparison law is
specified directly as `(dN/dp) proportional to p^a`.  For the local
compression DSA option,

`q=3r/(r-1)`, `f(p) proportional to p^-q`,

and spherical momentum volume gives

`dN/dp = C 4*pi p^2 f(p) proportional to p^(2-q)`.

The implementation normalizes analytically, including the logarithmic
exponent `a=-1`, and supplies exact CDF/inverse-CDF functions.  This is a
conditional net-first-passage spectrum at the finite upstream reference
surface; production rejects gross-shock-emission and unimplemented
accelerator-derived free-escape interpretations.

## Reference surface and transport-depth diagnostic

Every reference point must satisfy

`x_ref-x_sh = L_ref n`, with `L_ref>0`,

within the registered SI tolerance.  Stable patch IDs must be unique, the
offset must be upstream, and tangential leakage/folds fail before allocation.
`L_ref` is fixed while mesh and placement thickness converge.

The dimensionless diagnostic is evaluated by positive, ordered trapezoidal
quadrature,

`P_a(p)=integral_0^Lref [u_1n^in(d)/kappa_nn(d,p)] dd`.

Both `u_1n^in` and `kappa_nn` must be finite and strictly positive.  Constant
coefficients recover `u_1n^in L_ref/kappa_nn`.  The value is diagnostic: it is
never treated as a return/escape probability and never changes `L_ref`.

## Focused outward-flux sampling

For pitch cosine `mu`, the implementation forms the signed front-relative
normal speed from detector-frame advection and field-aligned motion, retaining
magnetic polarity explicitly.  With isotropic proposal density `a(mu)=1/2`,

`Z=integral_-1^1 a(mu)[w_n(mu)]_+ dmu`.

The tabulated conditional CDF is constructed only from the positive part.
`Z=0` gives typed `NoFocusedEscape`; a positive value below registered
absolute numerical resolution is `UnresolvedPositiveFluxNumerics`, never a
physical no-escape result.  Tangent fields are evaluated without division by
`b dot n`: positive advection can admit all pitch angles, while nonpositive
advection can have zero support.  A sampled production particle is therefore
never born with nonpositive outward relative speed.

The Parker no-through-flow sensitivity uses the conormal

`kappa n/(n dot kappa n)`.

It checks `n dot kappa n` before division and does not substitute geometric
normal reflection.  The shock-adjacent absorbing verification branch uses

`P_esc(delta) = expm1(u delta/kappa)/expm1(u L_ref/kappa)`,

which correctly tends to zero as placement distance tends to zero and remains
unavailable to production intent.

## Allocation, frames, and source support

Largest-remainder allocation assigns an exact integer sample total using
physical rates and stable input order.  Counter-based random draws are keyed
by campaign, stream, generation, tick, species, patch, and sample; adding a
new pitch or diagnostic draw cannot perturb a momentum stream.

Every compiled species retains stable semantic ID, compiled slot, chemical
symbol, mass, charge, and integer nucleon count.  All compiled entries are
validated; a source-enabled neutral is rejected.  Number-flux and nonthermal
energy budgets are independent.  Closed-field diagnose-only and
transition-clearance patches preserve diagnostics but return exactly zero
eligible rate.

The radial envelope is either a hard cutoff or the `C2` smooth quintic
`1-(10t^3-15t^4+6t^5)` between full-strength and zero radii.  It alters source
rate only, not front geometry or shock classification.  Four-momentum is
Lorentz-boosted from the declared source frame before insertion, with
`E^2-c^2|p|^2` invariant.

## Immutable release and loss ledgers

`InjectionCommitRegistry` prevents duplicate generation/tick commits.
`CohortLedger` retains the schema-stable terms separately:

- `CandidateReferenceRelease` and `NoFocusedEscapeExcluded`;
- `CommittedFirstPassageRelease` and `NetFirstPassageRelease`;
- `ImmediateShockAdjacentReturn` and `DelayedFrontReturn`;
- time-dependent `SurvivingUpstreamInventory`;
- verification-only `GrossShockAdjacentEmission` and
  `HitEscapeSurfaceBeforeShock`;
- independent `TransitionSheetContact`.

The cohort key contains species, stable patch lineage, provider generation,
and birth tick.  Finite-horizon reduction reports
`1-DelayedFrontReturn/CommittedFirstPassageRelease`, its dimensional numerator
and denominator, surviving inventory, available follow-up, right censoring,
and a typed zero denominator.  It is not labeled an asymptotic probability and
returned particles are not re-emitted.

Runtime loss caps use represented physical number and immutable birth kinetic
energy.  Event-time energy/four-momentum remain separate dimensional
diagnostics.  A nonremoving boundary rejects active front-loss caps.

## Verification

`make test-stage9` runs 173 cumulative tests.  Stage 9 adds `SRC3D01--19` and
`LOS3D02`: spectrum/CDF normalization, exact allocation closure, keyed-stream
independence, envelope smoothness, species and budget guards, Lorentz
invariance, transition exclusion, time integration, source semantics,
duplicate prevention, reference-surface orientation, focused support,
conormal/absorbing sensitivities, transport-depth quadrature, finite-horizon
accounting, and represented-number/birth-energy caps.
