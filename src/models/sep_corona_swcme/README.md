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
| BG3D-4 shock-fed spatial sheath | Partial | No: the empirical map and component identities are diagnostic only; a compatible finite initial state and common contact authority are absent |
| BG3D-5 ejecta/remaining regions | No | No |
| BG3D-6 material handoff | No | No |
| BG3D-7/8 native storage, MPI and restart | No | No |
| BG3D-9 native low-corona-to-1-AU campaign | No | No |

G1 is therefore open. In particular, the analytical trajectory reaching 1 AU
is not a native propagation receipt, and no planned owner/ghost or temporal-slot
behavior is described here as implemented.

## Implemented physics: BG3D-1 through BG3D-4

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

`event.conf` is a closed grammar; unknown or duplicate keys fail. The three
referenced assets are acquired by the host and SHA-256 checked before parsing.
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

Commented example fragments are in
`examples/bg3d1/event.conf`, `history.csv`, `ambient.asset`, and
`ejecta.asset`. The supplied event begins in the low corona, transitions from
the checksummed Hermite history through 100–200 s, and its analytical apex is
`1.59878e11 m` at 200000 s. That number verifies trajectory coverage only.

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
- `ARCHCSWC01`: source dependency scan, external public-header compile, and
archive-symbol audit excluding AMPS, PIC, and MPI dependencies.

Focused BG3D-4 execution is:

```text
make -C src/models/sep_corona_swcme -j16 build/test_bg3d4
(cd src/models/sep_corona_swcme && build/test_bg3d4)
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
