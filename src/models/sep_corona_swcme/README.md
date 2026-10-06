# Corona–SWCME shared analytical background

This sibling model owns the one continuous background event used first by
`srcSEP3D` and, only after G1 qualifies, by `srcSEP`.  It is dependency-light:
shared code has no AMPS, PIC, MPI, source, or particle type.  `model.md` is the
physics/design authority, `BG3D0_CONTRACT.md` freezes the selected G1 profile,
and this file documents the implemented surface of the current stage.

## Status by G1 stage

| Stage | Implemented | Independently validated |
|---|---|---|
| BG3D-1 input, launch and outer trajectory | Yes | Yes: `CMBGU01`–`02` |
| BG3D-2 ambient plasma/IMF | Yes | Yes: `CMBGU03` |
| BG3D-3 front and local shocks | Yes | Yes: `CMBGU04` plus maintained RH/ellipsoid suites |
| BG3D-3 shared contact/interface | Fixed-fraction reference only | No: reopened because it is not the BG3D-4 material contact authority |
| BG3D-4 spatial sheath/compression | Diagnostic map implemented; ejecta-driven per-ray production closure specified only | No: the selected production equations, finite ambient ray inventory, common contact authority, computed shock and acceptance suite are not implemented |
| BG3D-5 ejecta/remaining regions | No | No |
| BG3D-6 material handoff | No | No |
| BG3D-7/8 native storage, MPI and restart | No | No |
| BG3D-9 native low-corona-to-1-AU campaign | No | No |

G1 is therefore open. In particular, the analytical trajectory reaching 1 AU
is not a native propagation receipt, and no planned owner/ghost or temporal-slot
behavior is described here as implemented.

## Implemented and specified physics through active BG3D-4

`cme_event.h` and `cme_event.cpp` resolve the background-only event envelope.
The resolver accepts a closed `key=value` grammar, obtains bytes through an
injected reader, verifies lower-case SHA-256 checksums, and returns an immutable
configuration only after all validation succeeds.  Paths are provenance; byte
checksums and normalized physical assignments own identity.  Three assets are
mandatory even though their detailed physics is implemented later:

- `history`: a strict CSV containing time and value/rate for ellipsoid centre,
  radial semiaxis, and both lateral semiaxes;
- `ambient`: the BG3D-2 PFSS/Parker state asset;
- `ejecta_reference`: the future BG3D-5 reference map/vector-potential asset.

The history uses cubic Hermite interpolation. Positivity and outward apex
motion are checked using polynomial extrema between knots, not snapshots.
Exact triaxial heliocentric radial extrema establish the selected
`attached-then-detached` policy. The first valid plasma radius is separate from
the physical solar surface and from the initial front apex. Coverage must reach
at least 1.05 AU, while the frozen example's actual end state passes 1 AU.

At the handoff start, the complete ellipsoid supplies the crossing state. The
sign-aware constant-wind DBM law has no separately entered outer radius or
speed. All centre/axis dimensions continue from one reference surface. During
the finite transition, a quintic weight blends one position construction and
the code analytically differentiates it, including the weight-derivative terms.
Nested scans bound the non-polynomial transition acceleration. Queries outside
the declared interval fail rather than extrapolate.

The implemented stages do **not** yet provide ejecta material state, native
mesh publication, or G1 qualification. The vector-potential part of the asset
in `examples/bg3d1` remains a typed BG3D-5 input and cannot yet be presented as
spatial ejecta-field evidence.

### Units, frame and signs

Public values are SI: seconds, metres, kg m^-3, pascals, m s^-1, tesla and
webers. Vectors use the configured inertial HCI Cartesian frame. The event
basis has radial, first-lateral and second-lateral unit vectors; ellipsoid
polar/azimuth parameters are radians. Front normals point from the CME toward
upstream ambient. The positive shock-frame inflows used for admission are
`w1=Vn-U1.n` and `w2=Vn-U2.n`; an admissible forward shock has `w1>c_fast` and
`w2>0`. Magnetic sector is the ambient radial polarity, not the sign of
`B.n`. A sub-fast surface element remains geometric but admits no sheath mass.

### BG3D-2 ambient state

`ambient_state.h` and `ambient_state.cpp` consume the typed BG3D-1 ambient
record and compose maintained `sep_coronal_cme` kernels. From the physical Sun
through the PFSS source surface, the magnetic field is the checksummed harmonic
solution. Explicit two-direction RK4 tracing distinguishes open tubes from
closed loops. Closed plasma is rotating isothermal hydrostatic plasma. Open
tubes use the one transonic isothermal Parker speed and conserve mass per unit
magnetic flux. Outside the source surface, the same provider supplies the
retarded, rotating Parker IMF, angular-momentum velocity, and exact radial mass
flux through the declared coverage beyond 1 AU.

Each sample returns SI `rho`, electron density, pressure, electron/proton
temperatures, `U`, `B`, region and signed sector. A closed loop has sector zero;
open/exterior state has sector -1 or +1. A derivative sample additionally owns
event fingerprint, epoch, nonzero background generation, Cartesian `grad(B)`
and `grad(U)`, logarithmic magnetic gradient, and each central/one-sided stencil
choice. Stencils never cross a region or sector branch. The PFSS source-surface
query uses an explicit one-sided interior coordinate to avoid floating-point
shell overshoot; it does not alter or floor a magnetic value.

Invalid radial/time coverage, generation zero, a magnetic null, unresolved
topology, or a derivative stencil that cannot stay on one physical branch is a
typed failure. No Parker field, minimum density, or minimum magnetic magnitude
is substituted. The selected ambient is prescribed and isothermal: it is not a
global MHD solution or an observationally calibrated solar-wind forecast.

### BG3D-3 surfaces and local shocks

`surface_shock.h` and `surface_shock.cpp` bind the same immutable event and
ambient identities into one transactional epoch. They tessellate the complete
moving ellipsoid, evaluate every patch's analytical normal speed, sample its
upstream ambient state, and invoke the maintained canonical MHD jump solver.
The resulting records retain geometric, sub-fast, fast and supercritical
distinctions; local Mach number, obliquity and compression; upstream and
downstream states; Rankine--Hugoniot residuals; stable physical IDs; epoch and
background generation. A failed preparation leaves the last committed epoch
unchanged.

The legacy contact/ejecta reference is a rear-aligned nested ellipsoid derived
from the front. The configured `contact_apex_fraction` is retained for input
compatibility but is interpreted as the fraction of the front's rear-to-apex
span: all semiaxes are scaled by that fraction and the rear support point is
shared. It is a complete nested reference shape, not an independently reset
apex or radius, but it is **not qualified as the shared material contact**.
BG3D-4 currently uses the oldest material cohort instead. The resulting
cross-stage mismatch is retained as a failing physical qualification rather
than hidden by `contact_apex_fraction=0.99`.

This model is background-only. It explicitly disables particle-source
eligibility on every diagnostic shock patch and leaves all number and energy
source measures exactly zero. The compatibility default in the maintained
shock-provider input remains source-capable, so existing particle users retain
their baseline behavior.

### BG3D-4 shock-fed sheath

#### Selected production replacement: per-ray Lagrangian piston

The reviewed production direction is now
[`docs/BG3D4_PISTON_CLOSURE.md`](docs/BG3D4_PISTON_CLOSURE.md).  It is a
partially implemented capability: the independent contact authority,
checksummed coarse ray quadrature, conservative planar gas/MHD core and
spherical gas-dynamic geometry, curved magnetic energy/flux identities and a
frozen projected-ambient source path are implemented and independently tested.
Moving Parker-wind boundary transport, production ray initialization,
computed 3-D shock and assembly are not yet qualified.  One differentiable
ejecta/contact reconstruction is the piston and cross-stage contact authority.
For each supported, non-grazing HCI ray, a finite ambient column with physical
area `A_q=DeltaOmega_q*r^2` evolves in Lagrangian mass coordinates:

```text
partial_t r = u,
partial_t v = partial_m(A_q u),
partial_t u = -A_q partial_m(P+Q)
              - |B_t|^2/(mu0 rho r) - GM_sun/r^2 + f_amb,
partial_t e = -(p+Q) partial_t v + h_amb/rho,
partial_t [B_t/(rho r)] = S_b.
```

Here `v=1/rho`, `P=p+|B_t|^2/(2mu0)`, and `Q` is compression-only shock
dissipation.  The fixed ambient-maintaining sources make the projected
undisturbed ambient a steady solution and are applied to shocked plasma as an
explicit approximation.  The contact motion is prescribed; the compression
region and Level B shock are computed.  The current prescribed front remains
the Level A authority and becomes a Level B validation target, never a second
production shock.

The review identified pre-implementation issues that remain open:

- fixed-ray Lagrangian dynamics evolve `u*q_hat`; by user decision the
  published state adds unevolved ambient `U_t,a`, while omitted cross-ray
  transport is reported through `epsilon_t` and full-normal contact slip;
- `B_r/|B|` is not the physical `|B.n|/|B|` on a nonspherical contact or shock;
- the reduced radial compression speed differs from the canonical oblique-MHD
  fast speed when passive `B_r` is nonzero;
- ray solid-angle weights, contact-history handoff and all ambient magnetic/
  mechanical source work require explicit fingerprinted contracts; and
- any 3-D divergence-cleaning projection changes magnetic energy and interface
  flux and therefore needs a separate ledger.

The strict `pfss-parker-piston-sheath-v1` event selector requires separate
checksummed `assets.piston_contact_*` and `assets.piston_rays_*` records.  The example
`examples/bg3d4_piston/contact.asset` supplies a fixed-HCI triaxial ellipsoid
with a global 1200-s quintic rate ramp.  Its apex accelerates no faster than
`1875 m/s2`, reaches the prescribed `1.2e6 m/s` coronal speed, crosses the
`20 R_sun` handoff radius at `11528.3 s`, transitions without a position or
velocity reset to self-similar DBM motion, and reaches `1.81242e11 m` by the
end of the fixture.  `history.csv` remains a different, synthetic Level-A
shock surface; comparing the two is not event validation.

`examples/bg3d4_piston/rays.csv` is an explicit rotated 7-by-24 equal-solid-
angle quadrature.  Its 168 HCI unit vectors, physical weights and reciprocal
neighbor graph are parsed strictly and participate in the event fingerprint.
The common rigid `0.037 rad` rotation about HCI `+y` preserves the quadrature
and avoids putting qualification rays exactly on the analytical PFSS current-
sheet null; no magnetic floor is introduced.  The weights close to `4*pi`
steradians to parser precision.  At launch, four rays are supported
(`0.299199 sr`) and 164 are typed unsupported (`12.2672 sr`).  This asset
qualifies identity, parsing, complete unsupported-area accounting and the
first multi-ray smoke only; it is not sufficient angular resolution for
ray-count convergence or a cross-ray `epsilon_t` qualification.

`PistonContactModel::EvaluateRay` selects the outermost positive ellipsoid
intersection.  It computes the physical normal and `R_dot`, `R_ddot` by
analytic implicit differentiation of `g_c(R q_hat,t)=0`; it never differences
large heliocentric positions in production.  Rays below ambient support,
without an intersection, or with `q_hat.n_c < mu_min` receive explicit typed
dispositions.  `CheckStartupCompatibility` independently samples the ambient
on a supplied ray set and applies the frozen `|u_a|/c_f,a` gate.

`PlanarPistonSolver` is currently the reusable planar/spherical tube numerical
core (the historical class name is retained while its interface is stabilized).  It owns
immutable cell masses, staggered nodes, a moving material piston, an ambient-
comoving coordinate representation, compression-only quadratic VNR pressure
(`C2=2`, selected `C1=0`), RK4 time integration and a piston/outer work ledger.
The comoving coordinate makes a uniformly translating ambient an exact steady
discrete solution instead of repeatedly subtracting large translated node
positions to recover small widths.  Candidate evolution is private until all
cell volumes, densities, pressures and internal energies remain positive and
finite.  Shock-zone admission from `Q/P` is separate from the later availability
of two-sided plateau states.

The planar reference uses an independently evaluated exact piston shock.  Its
finest relative errors are `1.38e-4` in position, `2.28e-4` in compression and
`1.36e-3` in jump-derived speed.  A fixed-mesh CFL study reduces the normalized
energy/work residual from `6.87e-7` to a `1.3e-9` numerical floor.  The
accelerating-piston detector uses the required vanishing threshold
`epsilon_sh=2.6/N` for quadratic VNR because a shock is born at zero amplitude;
a fixed nonzero threshold has no continuum formation time.  The 6400-cell
formation result is within `0.33%` in time and `0.16%` in position.  The
linear-pulse density/velocity and propagation-centroid errors converge to
`5.64e-4` and `4.16e-4`.

The same core now evolves both signed tangent-basis components of the
transverse field, preserves their planar/spherical frozen-flux invariants and includes
magnetic pressure.  The cold perpendicular `gamma=2` piston and an independent
finite-beta canonical RH comparison pass at three resolutions; the finest
shock-position errors are `9.92e-5` and `1.10e-4`, respectively.  Exact
spherical-sector face areas and shell volumes are also implemented.  Starting
from an independently integrated finite-time Taylor similarity state avoids
misrepresenting the singular strong-shock solution as an ambient-only startup.
At 150/300/600 cells, spherical shock-radius error decreases
`6.93e-3, 3.43e-3, 2.04e-3`, piston-pressure error decreases
`3.16e-2, 7.87e-3, 2.03e-3`, and the finest shock-speed error is `1.90e-2`.
This manufactured initialized state is a reference/restart interface only;
production still starts from compatible ambient material at the contact.
For a dynamic spherical magnetic compression, normalized energy residuals
decrease `7.55e-6, 3.79e-6, 1.90e-6`; both signed
`B_t/(rho*r)` components remain invariant to less than `2.5e-16` relative.

`PistonAmbientProjection` freezes one maintained 3-D ambient epoch on a fixed
HCI ray and evaluates the declared `f_amb`, `h_amb` and two-component `S_b`
formulas with branch-preserving radial stencils.  A deterministic source table
avoids rerunning PFSS topology tracing at every RK stage.  It never
interpolates across a topology/sector or PFSS/Parker branch change; those rare
intervals fall back to a one-sided authority evaluation.  The solver applies
body acceleration, volume heating and induction at every RK stage and records
their work separately from piston and outer-boundary work.

`CSWC0625` uses the actual nonuniform closed-corona PFSS/plasma profile over a
spherical tube.  At 64/128/256 cells, the density, pressure, velocity and
magnetic equilibrium errors decrease by approximately four per refinement;
their finest values are `6.36e-7`, `1.06e-6`, `6.09e-7`, and `6.38e-7`.
The fixture reports `max |f_amb|/|g|=0.559104`; its closed-ray heating is
exactly zero.  A separate uniform manufactured source exercises nonzero
heating and both induction components, closing the independent source-work
ledger at `2.12e-15`.  The same test advects the actual Parker profile between
independently integrated ambient material boundaries.  At 64/128/256 cells
its finest density, pressure, velocity and full-field errors are `3.43e-9`,
`4.75e-9`, `6.37e-9` and `9.59e-10`; its energy residual decreases
`4.88e-7,2.44e-7,1.22e-7`, and `B_r r^2` closes to `4.02e-16`.

`PlanarPistonSolver::AppendAmbient` transactionally extends the immutable
material domain with complete ambient mass, momentum, thermal and magnetic
state.  It records transported total energy and rejects a mismatched shared
face without modifying committed state.  `CSWC0626` verifies this through the
actual Parker projection; the finest post-append density, pressure, velocity
and magnetic errors are `8.57e-9`, `9.12e-9`, `2.48e-8` and `2.40e-9`, with
exact reported mass closure and decreasing energy error.

`PistonSheathModel` is the first contact-driven collection of Level-B tubes.
The checksummed generic contact is the only piston authority; Level-A
`history.csv` is deliberately absent from its evolution.  It initializes
finite ambient material, advances every supported ray privately, maintains an
undisturbed outer buffer by transactional append, commits only when all ray
solves succeed, detects the numerical shock, and queries typed ejecta-side,
compression/sheath and ambient regions.  `CSWC0627` evolves four supported
rays to `1800 s` and exercises these interfaces.  Its sampled ray has
`R_c=2.13043e9 m`, `R_sh=2.59383e9 m`, compression `4.34277`, and 42 ambient
buffer cells.  It is quasi-parallel (`theta_Bn=0.067924 rad`), so the
quasi-perpendicular canonical-RH tolerance is explicitly not applicable; its
compression/pressure/velocity discrepancies are still reported as model
error.  This smoke PASS therefore does not qualify the production-RH subset,
angular convergence, `epsilon_t`, 3-D assembly or BG3D-4.

The detailed specification records SI units and HCI signs, unknowns, equations,
initial/boundary conditions, shock detection, region assembly, strict inputs,
fourteen independent acceptance cases and limitations.  Contact, ray and VNR
finite-volume controls are now typed and fingerprinted.  The first production
mass-refinement attempt remains failed because disjoint shock-zone count and
jump/contact quantities do not converge; input identity alone does not qualify
the closure.  Until the remaining convergence/validity/assembly tests pass,
`BG3D-4` remains open.

#### Retained relaxation diagnostic

The implemented `rh-relaxing-material-map-v1` is an **experimental prescribed
map**, not a fluid evolution solver or a qualified production closure. A parcel
label is `(theta,phi,tau)`: two complete front coordinates and its unique
shock-crossing time. Define the current RH/front velocity deficit
`D(q,t)=U2(q,t)-d_t X(q,t)`, age `a=t-tau`, and

```text
L(a)  = kappa*a + (1-kappa)*T*[1-exp(-a/T)],
L'(a) = kappa + (1-kappa)*exp(-a/T),
x(q,tau,t) = X(q,t) + L(t-tau)*D(q,t).
```

The frozen example uses `kappa=0.25` and `T=25 s`. Consequently `x=X` and
`U=U2` as `t -> tau+`; shock compression is applied once. The constants avoid
the singular old-cohort volume of complete drift arrest, but they are empirical
and do not follow from momentum balance, a piston solution, or an ejecta model.
The derivative bases with respect to `(theta,phi,tau)` at birth and at the
requested epoch define `F_rel=A*A0^-1` and `J=det(F_rel)`. The returned state is

```text
rho = rho2/J,   p = p2 J^(-gamma),
B = F B2/J,     U = partial_t x,     J > minimum_jacobian.
```

This is the Cauchy ideal-induction construction with zero added heating. It
kinematically preserves mass, flux and adiabatic pressure; it does not by
itself prove momentum balance, energy balance, or observational accuracy.
Fourth-order centred
derivatives evaluate the smooth birth-deficit field, with explicit one-sided
stencils at physical support boundaries. If no stencil stays on the same fast
branch, if the crossing-volume basis is singular, or if `J` falls below the
configured minimum, the query fails. No field, density, pressure or Jacobian
floor is substituted.

Inventory quadrature uses the actual shock area at each birth epoch:

```text
dM = rho1 w1 dA dtau = rho2 w2 dA dtau.
```

Every patch/time cell has a unique lineage. Photospheric and outer-front exits
are found from signed geometric event functions with bounded bisection to
`root_tolerance_s`; exited material is ledgered at that time and is never
extrapolated into a later map fold. The closed angular domain has no unrecorded
lateral edge. Fast/sub-fast support is recorded separately. A sub-fast
candidate is not relabeled as a solved compression wave: candidate preparation
fails transactionally, while `QueryCommitted` continues to return the stored
material cells from the last complete epoch without re-running current RH
admission. No compression-region state is claimed.

The current startup is a declared zero-volume limiting fixture. It is useful
for admission tests but is not a justified finite production inventory or
prehistory. Running the same empirical relaxation for longer would not repair
that deficiency. A replacement evolution law will require a separately
reviewed closed mathematical construction with unknowns, equations, compatible
finite initial data, contact and shock boundary conditions, and a determinate
solution procedure.

The example's `contact_apex_fraction=0.99` describes only the unqualified
fixed-fraction BG3D-3/ejecta reference. It is not used to define the BG3D-4
oldest-cohort material surface. The cross-authority regression currently
measures a maximum relative position mismatch of approximately `1e-2`, so the
shared contact/interface contract is explicitly open.

Current validity is deliberately limited to admitted fast patches and their
first physical solar/outer exit. Sub-fast compression-layer plasma, the ejecta
interior, post-event recovery, a qualified contact, full-run material handoff,
and native publication remain missing. A local RH state or retained parcel
inventory must not be used to fill those regions with quiet ambient.

### Independent fixtures and production residual diagnostics

`test/bg3d4_reference_fixtures.h` contains two manufactured solutions that do
not call the production map for their expected values:

- a finite planar slab with an initial material inventory, material piston,
  constant-speed Mach-3 shock, exact ideal-MHD jump, nonzero uniform normal
  magnetic field, and independently balanced moving-control-volume energy;
- a finite spherical shell with `X=lambda(t)*a`, `lambda=1+alpha*t`,
  `rho=rho0/lambda^3`, `p=p0/lambda^(3 gamma)`, and a uniform Cartesian
  `B=B0/lambda^2`. The field is divergence-free and threads a declared
  transverse-field contact; it is not a closed flux rope. With linear lambda,
  acceleration, pressure gradient, curl(B), gravity and manufactured body
  forcing are zero. Pressure and Maxwell boundary work exactly account for the
  internal and magnetic energy changes.

`sheath_diagnostics.{h,cpp}` samples the actual public production map. Smooth
fourth-order material-coordinate/time stencils independently form inertia,
pressure-gradient, Lorentz and solar-gravity force densities, their residual,
and the conservative ideal-MHD total-energy residual. Volume integration uses
the production map's curved physical cell volumes and reports absolute force,
volume-weighted local P99, signed/absolute energy residual, resolved work terms,
and signed/absolute residual-force work. The current distributed result is
stable under the three tested time steps near an integrated force ratio of
`0.646`; it is therefore evidence of persistent approximation/model error, not
a demonstrated discretization error. Conversely, the integrated independently
differentiated contact-leakage estimates do not converge and have comparable
derivative uncertainty. That is a numerical verification failure: it neither
proves a resolved physical leak nor permits the identically zero self-reported
flux to qualify the interface. These diagnostics do not qualify the map.

## Inputs and commented example

`event.conf` is a closed grammar; unknown or duplicate keys fail. The selected
profile's three legacy or five piston assets are acquired by the host and
SHA-256 checked before parsing.
The principal input groups are:

| Group | Meaning and constraints |
|---|---|
| `run.*`, `domain.*` | profile, inertial frame, time interval, physical Sun, first plasma radius and >1-AU support |
| `geometry.*` | fixed direction/tilt and attachment policy; complete center/three-axis history comes from the CSV |
| `plasma.*` | one composition, species temperatures, alpha abundance, electron-mass convention and `gamma>1` |
| `ambient.*` asset | PFSS harmonics, source/reference radii, density, rotation, trace and Parker-table controls |
| `regions.*` | exact sheath/ejecta closure identifiers; unsupported names reject |
| regional asset | admission start, contact span fraction, ejecta reference state/fluxes, minimum J and force/heat/work caps |
| `handoff.*` | transition endpoints, ambient wind and signed DBM drag; outer state is derived from the crossing |
| `numerics.*` | root-time and maximum prescribed-acceleration gates; these are not physical floors |
| `assets.piston_*` | independent contact and explicit HCI ray quadrature; both checksums join the event fingerprint |
| `bg3d4.*` | finite ambient cells/buffer, source and trajectory resolution, append policy, disturbance/shock thresholds, CFL, VNR coefficients and well-balanced-source selector |

Commented example fragments are in
`examples/bg3d1/event.conf`, `history.csv`, `ambient.asset`, and
`ejecta.asset`. The supplied event begins in the low corona, transitions from
the checksummed Hermite history through 100–200 s, and its analytical apex is
`1.59878e11 m` at 200000 s. That number verifies trajectory coverage only.

The piston example in `examples/bg3d4_piston/event.conf` supplies the accepted
`bg3d4.*` controls with comments.  They populate `PistonNumericsInput` and join
the immutable fingerprint; `PistonSheathModel::Create` has no programmatic
override.  Changing resolution or viscosity therefore resolves a distinct
event.  The controls are not aliased to legacy relaxation or particle scalars.
Validity/assembly tolerances not yet implemented remain absent rather than
being accepted and ignored.

## Build and test

From the AMPS root:

```text
make -C src/models/sep_corona_swcme -j16 test
```

This is a standalone shared-library build, not a native AMPS rebuild. It runs:

- `CMBGU01`: strict grammar, checksums, asset roles, closure capabilities,
  path-independent identity, and missing/corrupt/unsupported negatives;
- `CMBGU02`: continuous history, exact extrema, attachment/detachment,
  no-extrapolation, actual 1-AU coverage, phase selection, and an independent
  numerical derivative of transition position versus reported velocity;
- `CMBGU03`: open/closed/exterior plasma, signed sectors, EOS, Parker mass,
  momentum, winding and angular-momentum identities, typed null/coverage
  failures, a separate five-point tensor oracle, divergence, and three sphere
  resolutions for signed magnetic flux;
- `CMBGU04`: complete front and fixed-fraction reference geometry, stable
  identities, analytical versus finite-difference normal motion, three-level area convergence,
  sub-fast classification, both-sided fast-shock states, Mach/compression/
  obliquity and RH residuals, transaction failure integrity, and zero particle
  eligibility/rates/measures;
- `CMBGU05`: exact finite-inventory planar shock/piston and homologous curved
  fixtures; exact RH boundary limit and positive-J Cauchy state; three-level
  curved-volume, admission-time and material-identity checks; explicit
  cross-stage contact mismatch; independently differentiated absolute leakage;
  transactional sub-fast rejection with committed-material readback; and
  distributed production momentum/energy/work diagnostics. Passing this test
  validates the diagnostics and manufactured solutions, not BG3D-4 physics;
- `CSWC0620`: independent Level-B contact asset/selector, separation from the
  Level-A shock, checksum/fingerprint mutation, global C2 startup ramp,
  attachment/detachment, smooth contact handoff, 1-AU coverage, outermost ray
  intersection, analytic implicit derivatives versus independent history
  differences, unsupported-flank typing and ambient startup compatibility;
- `CSWC0621`: conservative staggered planar Lagrangian gas piston, exact
  translating-ambient balance, independently evaluated piston-shock position/
  compression/speed, mass and CFL refinement, energy/piston-work closure and
  transactional failure;
- `CSWC0622`: uniformly accelerating piston shock-formation time/height with a
  vanishing quadratic-VNR detector, plus a small-amplitude sub-fast pulse's
  density/velocity invariant and propagation-centroid convergence;
- `CSWC0623`: cold perpendicular and finite-beta MHD piston references with
  canonical RH comparison, signed transverse frozen flux and energy ledgers;
- `CSWC0624`: independently integrated Taylor spherical-piston reference,
  exact curved metrics and a separate spherical magnetic-work refinement;
- `CSWC0625`: actual PFSS/closed-corona and Parker projected-ambient balance,
  plus independent nonzero source-work and passive-radial-flux checks;
- `CSWC0626`: transactional Parker ambient append, immutable material mass and
  complete appended-energy transport;
- `CSWC0627`: contact-driven multi-ray evolution, computed shock, maintained
  buffer, regional query and all-ray transactional commit.  The sampled
  production ray is outside the quasi-perpendicular RH acceptance subset and
  is reported as `NOT-APPLICABLE`, not counted as that physical gate;
- `ARCHCSWC01`: source dependency scan, external public-header compile, and
archive-symbol audit excluding AMPS, PIC, and MPI dependencies.

`CSWC0620`--`CSWC0627` are registered and passing for their explicitly bounded
scopes. `CSWC0628`--`CSWC0633` remain reserved requirements covering
mass/time/ray convergence, validity diagnostics, 3-D assembly, complete
transactional failure, sanitizers/regressions and comparison with the Level A
shock.  The unexercised quasi-perpendicular production-RH subgate also remains
open.  Component PASS results must not be promoted to BG3D-4 qualification.

Focused BG3D-4 execution is:

```text
make -C src/models/sep_corona_swcme -j16 build/test_bg3d4
(cd src/models/sep_corona_swcme && build/test_bg3d4)
make -C src/models/sep_corona_swcme -j16 build/test_bg3d4_piston_contact
(cd src/models/sep_corona_swcme && build/test_bg3d4_piston_contact)
make -C src/models/sep_corona_swcme -j16 build/test_bg3d4_planar_piston \
  build/test_bg3d4_planar_formation build/test_bg3d4_planar_mhd \
  build/test_bg3d4_spherical_piston build/test_bg3d4_ambient_balance \
  build/test_bg3d4_append build/test_bg3d4_piston_sheath
(cd src/models/sep_corona_swcme && build/test_bg3d4_planar_piston && \
  build/test_bg3d4_planar_formation && build/test_bg3d4_planar_mhd && \
  build/test_bg3d4_spherical_piston && build/test_bg3d4_ambient_balance && \
  build/test_bg3d4_append && build/test_bg3d4_piston_sheath)
```

The maintained coronal regression, which includes the exact moving-planar
sheath and weak/oblique/sub-fast RH references, is:

```text
make -C src/models/sep_coronal_cme -j16 test
```

Expected outputs are named `CMBGU*-EVIDENCE`, `*-REFERENCE`,
`*-CONTACT-AUTHORITY`, `*-CONTACT-INTEGRAL`, `*-SUBFAST-TRANSACTION`,
`*-DIFFERENTIAL`, and `*-DISTRIBUTED*`. They report dimensional inventory,
force, leakage, energy and work values separately from normalized numerical
residuals. Qualification uses the frozen budgets in `BG3D0_CONTRACT.md` and the
pre-implementation BG3D-4 criteria in `CODEX_CME_PLAN.md`; a diagnostic test
passes when it correctly exposes an unqualified physical construction. Required
campaign skips and failed physical gates remain open rather than becoming PASS.

## Ownership, epochs and native boundary

`EventConfiguration` is immutable and owns the physics fingerprint. One
`AmbientModel` uses it for all radii. Surface/shock and sheath preparation
build private candidates and publish an immutable epoch only after every patch,
state and ledger passes; zero generation is reserved and failure preserves the
previous handle. Background preparation has no source plan, random stream,
particle allocation or feedback prerequisite.

The shared library contains no AMPS, PIC or MPI ABI. BG3D-7 must add a thin
srcSEP3D adapter that stages current/previous native slots privately, performs
owner write/readback and actual halo exchange, joins collective readiness, and
only then publishes the generation. BG3D-8 must validate restart and
pre/inside/post-handoff epochs. None of that native storage or synchronization
is implemented or validated at the present stage.

The existing `ELL3D01`–`ELL3D10` tests independently retain ellipsoid implicit
geometry, normals, surface velocity, area convergence, solar clipping, stable
patch identity, equivalent center/apex construction, analytic histories,
knot continuity, basis orientation, and piston nesting. They remain shared
geometry evidence rather than being duplicated under new IDs.

Before a later native AMPS rebuild, follow the repository-root `AGENTS.md`
four-step clean/regenerate/hooks/`-j16` workflow exactly. The command above does
not remove or regenerate the native root `build/` directory.
