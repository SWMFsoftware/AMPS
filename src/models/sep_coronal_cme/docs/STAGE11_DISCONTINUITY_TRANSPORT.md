# Stage 11: independent discontinuity transport capabilities

The public API is `sep_coronal_cme/discontinuity_transport.h`. The supplied
families have explicit geometric/physical domains; passing their independent
manufactured gates does not qualify an arbitrary CME sheath or a finite
PFSS/SCS transition sheet. The application moving-source baseline remains
available without either optional provider.

## 11A: finite exterior HCS

`FiniteHcsSheet::Create` prepares an owning immutable, stationary planar
force-free rotational sheet. With normal `n`, tangent `t`, signed distance `d`,
and independently supplied physical half-thickness `a`, its field is

`B = B0 [tanh(d/a) t + sech(d/a) (n cross t)]`.

Density, pressure, and physical outward/inward wave energies are constant;
fluid velocity and electric field are zero in this stationary family.
`B.n=0`, `div B=0`, `|B|=B0`, and `J cross B=0`. Thus the plasma and magnetic
pressure are in equilibrium and no artificial weak-field null is introduced.
The center has explicitly invalid signed-sector wave labels, with finite
sentinels; physical outward/inward energies remain defined and conserved.

`Advance` integrates full relativistic momentum using a reversible
half-drift/Boris-rotation/half-drift map. It bounds both gyro angle and travel
relative to the physical thickness. A signed-distance crossing is bisected
on the actual orbit map, the step is split at the root, and the remainder is
recomputed from the crossing state. Magnetic rotation conserves momentum
magnitude and energy; pitch angle is derived from the local full field rather
than passed through an ideal sign reversal. Negative duration provides the
reverse manufactured experiment. An explicit substep budget fails instead of
silently enlarging the permitted step.

`fieldGridSpacingM` is zero for the analytic field. Positive spacing (no larger
than `a`) samples its direction on an explicit normal grid and renormalizes
field magnitude. This still has zero normal field, solenoidal geometry and
constant magnetic pressure, and supplies an independent field-resolution
convergence control. It is not an AMPS interpolation callback.

Ideal HCS, finite HCS, separatrix, and PFSS/SCS transition identities are
distinct enum values. A transition/separatrix passed to the finite-HCS factory
fails. The Stage-3 signed-gauge PFSS/SCS join and transition-clearance rule
are unchanged. No finite sheet is glued to that overlap by this implementation.
Moving/curved finite sheets, guiding-center approximations and wave-force
feedback require a separately constructed and qualified provider.

## 11B: finite downstream planar sheath

`PlanarShockSheath::Create` solves the existing regular oblique fast RH branch
for a constant upstream state and a uniformly moving planar front. Its finite
rear boundary has `L(t)=L0+Ldot(t-epoch)` and speed `Vshock-Ldot`; the entire
requested time interval must have positive thickness. Every downstream query
returns the same solved RH state. Queries behind the rear boundary or outside
time coverage fail; the exact front requires an explicit side.

This is a global planar extension, not independently assembled local ellipsoid
patches. Normal magnetic field is constant throughout each volume and
continuous at the front. Opposite tangential side-face fluxes cancel. The
finite-volume audit checks front mass/normal flux/energy continuity and the
moving-volume inventory balance. The rear's mass and energy throughput are
accounted explicitly, so an expanding control volume does not acquire
unreported mass or energy. The independent audit gates factory publication.
Uniformity makes spatial/time refinement exact to floating-point roundoff;
this is the stated analytic accuracy of this family.

Waves are passive test fields: each directional wave's normal energy flux in
the shock frame transmits completely, with zero reflection and downstream
energy density scaled by the ratio of upstream/downstream characteristic
speeds `|u_n-Vshock +/- sector*B_n/sqrt(mu0 rho)|`. The sector is supplied by
the background; the shock-normal sign of `B` cannot define the physical
outward/inward labels. A singular characteristic fails.
This supplies a specified transmission policy and a separate wave-flux audit;
it does not claim nonlinear turbulence transmission, wave backreaction, or
oblique-wave reflection from a general CME.

The factory checks the upstream and downstream wave energy against the
declared `maximumPassiveWavePressureFraction` of thermal plus normal ram
pressure. A state above that bound fails rather than dropping significant
wave backreaction from the RH equations.

The zero-potential mathematical crossing leaves inertial four-momentum
continuous. `Cross` explicitly transforms it into the outgoing local plasma
frame and records local pitch angle and shock-frame energy. It invents no
instantaneous DSA kick; subsequent scattering/transport owns physical energy
changes. A particle/surface/generation ledger makes repeated callbacks
idempotent and can be restored for restart. There is one committed crossing
per registered surface generation and particle lineage.

`Advance` supplies the actual unscattered full-vector particle step, using
one-sided uniform `B` and ideal-MHD `E=-u cross B`. Its relativistic Boris map
bounds gyro angle and displacement through the minimum covered sheath width.
The shock root is bisected on the incoming map; the outgoing field is used
only for the remainder after the crossing transformation. It never blends the
shock into a grid cell. Negative duration supplies the reverse experiment.
The ledger is copied for the step and published only with a successful final
particle state; an exit through the finite rear or a failed coverage query
cannot leave a partially committed event. The clock is anchored to the
requested interval, avoiding coverage failures from repeated-step roundoff.
Scattering is a separate operator and must use the recorded local plasma
state and transmitted waves. The host remains responsible for applying the
coincident-event ordering before combining magnetic and shock operators.

`LocatePlaneCrossing` event-locates a moving plane on a space-time segment
independently of mesh/cadence. Endpoint ownership avoids duplicate roots.
Coincident events process unqualified transition/separatrix exclusions first,
qualified magnetic operators next, and shock transformations last; identity
breaks ties deterministically.

## Independent selection and integration boundary

`ValidateDiscontinuityCapabilities` requires only the provider actually
requested. An HCS handle does not enable a sheath, and vice versa. Requesting
a finite composite transition still fails. These public owning providers are
usable by a host's full-vector mover; schema-5/application selectors are not
silently repurposed to activate them. Existing application input grammar lacks
the geometric/state records needed to select these supplied families.

A production curved CME extension must independently construct its global
field/flow, pass its global conservation audit, and supply a compatible host
adapter. The present passing planar gate is not permission to overwrite the
analytic PFSS/SCS/Parker background with per-patch downstream values.

## Verification

```bash
make -j8 test-stage11
```

The cumulative gate passed 202/202 tests. New canonical records:

- `HCS3D04`: profile field/pressure/density/wave invariants, independently
  differenced force-free condition, signed-flux cancellation and thickness guard.
- `HCS3D05`: located crossing, drift/time convergence with particle step and
  field-grid spacing, forward/reverse energy and orbit convergence, endpoint ownership.
- `HCS3D06`: thin-sheet one-sided limits, separate ideal operator and rejected
  reuse of transition/separatrix identities.
- `SHEATH3D01`: moving-front cadence independence, one commitment per lineage/
  generation, frame/mass-shell conservation, reversal, ledger restore and
  coincident ordering, actual split forward/reverse full-vector orbits, and
  atomic ledger rollback on a failed rear-boundary step.
- `SHEATH3D02`: finite global divergence, mass, magnetic-flux, energy and wave
  audits at several spatial/time resolutions, inventory accounting, one-sided
  state and coverage/subfast rejection, both background-sector wave labels,
  and rejection above the declared passive-wave bound.
- `SHEATH3D03`: independent capability gates and rejected composite promotion.

These are shared-model verification tests. Full AMPS/MPI mover integration and
an event-specific observational qualification remain separate evidence.
