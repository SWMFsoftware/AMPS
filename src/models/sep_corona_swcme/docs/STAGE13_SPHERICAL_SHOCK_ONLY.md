# First spherical SHOCK_ONLY core candidate

The user-approved Stage 13 work starts with the maintained PFSS/open-field
and exterior Parker ambient and the exact Sun-centred spherical SWCME3D front.
The first physical particle spheres are **1.05 solar radii and 5 AU**. Provider
coverage is a separate explicit input: retained straddling material, derivatives
and numerical stages may require actual guards beyond those physical spheres.
No fixed 8-AU or 10-AU guard is installed as a production default.

The shared core provides the original kinematic material map and a resolved
thermal-only wind over a spatial band of open tubes. Either can drive a
spherical event's ambient evaluation, material map, and local shock state.
Native SEP transport needs an observationally qualified wind realization and
the corresponding moving-line lifecycle.

## Resolved thermal-only wind on an open tube

`OpenTubeIsoPolyWind::Prepare(event->Ambient(), request)` uses the actual
axisymmetric PFSS field from the immutable event. A request supplies the
source-surface colatitude in radians, a positive mass loading in kg/(s Wb),
the reference mass density in kg/m³, separate thermal-transition radii and
polytropic indices, and explicit trace and wind accuracy. The event supplies
the species reference temperatures and composition. The controlled example
used by `S13SC04` is in `test/test_spherical_shock_only.cpp`. Its numbers are
manufactured. An observational realization requires separately evidenced
thermal and mass-loading inputs over its declared radial support.

Below the source surface, the solver traces the meridional field with
`dtheta/dr = B_theta/(r B_r)` and uses cross section perpendicular to `B`.
The speed in the projected equation is along the tube. The radial coordinate
parametrizes the monotone curve, so every derivative acquires the same
`ds/dr` factor. Rotation enters through the centrifugal potential. Above the
source surface, the azimuthal velocity satisfies `r U_phi = constant`, fixed
by corotation at the handoff. The equivalent radial momentum equation uses
radial speed, radial flux area and `dPhi/dr = GM/r² - U_phi²/r`. Both
representations carry the same normal mass flux, normal magnetic flux, and
tangential ideal electric field at the handoff. The exported velocity and
magnetic field are inertial Cartesian vectors in SI units. The tube's
colatitude and any requested inertial azimuth determine the sampled position.

The shared `IsoPolyTubeWind` solves the sonic numerator and denominator on
the same pressure law. It starts both directional integrations from the
regular sonic slope and verifies strict subsonic and supersonic endpoints.
The isothermal limit additionally requires a unique global nozzle minimum.
For a species with a thermal transition, the solver blends log temperature
with the full quintic derivative in the momentum equation. It returns the
net species heat rates implied by that closure in W/m³. Those rates require
an independent energy budget before field use. The tube export uses `U_r`
for advection of radial derivatives. In the nonradial corona this is
`u dr/ds`, so the heat carries the same material derivative that the steady
speed and density imply. The steady closure derivative and the adiabatic
derivative used in the shock characteristic are distinct.

`Event::SurfaceOnTube` evaluates the same canonical spherical radius and
normal speed at the query epoch. It accepts the resolved tube from that event's
ambient owner and applies the existing oblique fast-mode and RH gates to its
actual upstream primitive. Missing support and unresolved branches return an
explicit status. Sub-fast geometry has no jump. This narrow API supplies
along-tube derivatives. A continuous family of tubes and its transverse
derivatives are required for the native material map. A nonaxisymmetric
magnetogram requires a conservative longitude map when winding depends on
the source tube. Radial turning points, unresolved magnetic sectors, and
out-of-range epochs are outside the supported chart.

## Spatial family and the spherical event

Set `SphericalShockOnly::Definition::resolvedWind` to an
`OpenTubeIsoPolyFamilyInput` with a source-surface colatitude band and one
identified population. Supply the population's mass loading in kg/(s Wb),
reference mass density, species thermal bands, and numerical trace/solve
controls. The angle bounds are material surfaces measured at the PFSS source
surface. Construction uses a nonpolar, axisymmetric, monotone open band of
one magnetic polarity. The manufactured example is `ThermalFamily()` in
`test/spherical_shock_only_fixture.h`.

Inside the corona, each HCI point seeds the meridional trace in both radial
directions. The outward branch locates its source label, and the inward branch
completes the declared guard coverage. Outside the source surface the
colatitude is already the source label because the meridional component is
zero. The transonic solution uses that tube's area, rotation and thermal law.
This implementation solves at each query. It makes the solution independent of
an angular interpolation or cache and carries a measurable cost per query.
The bounds and population affect the fingerprint used by the event and
detached interval. A failed tube query leaves the field unavailable at that
point. The manufactured population supplies software verification inputs.

The Cartesian derivative stencil resolves neighboring tubes at a common epoch
and uses the same branch for velocity, field, density and total pressure. A
one-sided stencil handles the PFSS–Parker interface. A failed transonic solve
at a derivative point propagates its status. The event reports the spatial
continuity and induction residuals, the vector force/work needed to sustain
the background, and total and per-species implied heat. The momentum solve is
field-aligned. The transverse force diagnostic remains a qualification input.
`Evaluate`, `Move`, `AdvanceFixed`, `Surface` and `Connect` all use the same
resolved owner when configured. `SurfaceOnTube` remains a separate single-tube
diagnostic for an event without a spatial family.

`S13SC11` verifies neighboring tubes against an analytical dipole flux label,
an independent Cartesian stencil, signed mass/field budgets, the sourced
heat sum and one-sided branch behavior. `S13SC12` exercises full 300 s steps
through the source surface and at 1 AU. It verifies the Cauchy magnetic update,
mass Jacobian, inverse map, inertial longitude advance and oblique fast Mach
number from the upstream state. The linked `S13AD01` also checks the spatial
profile identity, local e/p/alpha temperatures and gradient translation.

```sh
make -C src/models/sep_coronal_cme -j8 test-tube-isopoly
make -C src/models/sep_corona_swcme -j8 test-spherical-shock-only
```

## Kinematic material API and physical ownership

`sep_corona_swcme/spherical_shock_only.h` declares an immutable
`SphericalShockOnly::Event`. `Definition` explicitly contains:

- the complete `AmbientDefinition`, including composition/EOS, PFSS harmonic
  bytes, source surface, reference density and radius, thermal anchors,
  rotation and actual spatial/time support;
- canonical SWCME ballistic, DBM or monotone data-driven front kinematics,
  with no extrapolation;
- physical inner/outer particle radii, independent of provider guards; and
- an explicit input identity.

`Create` seals the complete typed ambient and event identity. Changing a
harmonic, species temperature, profile convention, front history or support
cannot silently reuse a stale caller label. It rejects the spherical-density
shortcut for this nonradial open-field candidate: below the source surface the
actual canonical density follows mass loading per signed magnetic flux, not a
radial inverse-square prescription applied to a nonradial tube.

`Evaluate` retains the actual canonical open PFSS/Parker state and both full
Cartesian derivative tensors. The density and pressure gradients share the
same branch-fixed finite-difference stencil as the existing vector gradients;
closed loops and unresolved HCS/null states reject. Species temperatures,
charge neutrality, alpha abundance and the electron-mass convention remain
owned by the canonical EOS. Vectors and tensors use inertial HCI Cartesian
components and SI units.

`PhysicalSupport` distinguishes inner guard, physical particle domain and
outer guard. `FrontSide` also distinguishes upstream, exactly on the front and
`BehindFrontReferenceOnly`. The latter is an undisturbed reference used for
geometry/diagnostics, not a downstream ambient volume. A consumer must retain
the prescribed first-return absorbing shock rule; this core does not invent
pass-through transport, a sheath, ejecta or a downstream heating partition.

## Material map and exact physical step

`Move` integrates the actual canonical HCI velocity and its variational
system together:

```
dX/dt = U(X,t)
dF/dt = grad(U)(X,t) F,  F(begin) = I,  J = det(F).
```

A controlled RK4/step-doubling controller subdivides the mathematical interval
on an endpoint-preserving bisection mesh. Position and deformation tolerances,
maximum numerical time/arc step and work count are explicit caller inputs.
Each smooth panel receives a duration-proportional share of the interval
budget. The resolved thermal flow reserves half the error budget for panels
straddling the PFSS–Parker interface. Its continuous velocity has a jump in
the one-sided gradient there. The solver refines those panels against an
absolute increment estimated from stage rates, since the smooth RK4 error
estimate loses its order at the interface. Every stage and endpoint query real open-field/guard
support; missing support and exhausted numerical work reject. Reverse-time
queries implement the inverse flow, not reverse particle stepping.

`AdvanceFixed` requires `end == begin + fixedDt`, positive unchanged
`fixedDt` and a complete physical interval. It integrates using that explicit
duration instead of inferring a shortened physical step from a host endpoint.
RK stages are numerical quadrature, not particle sampling timesteps. No method
reads or changes PIC clocks, particle buffers, sampling lengths or global dt.

Independent tests use `rho_end J = rho_begin` and
`B_end = F B_begin/J`, including a nonradial coronal parcel and an actual
source-surface crossing. These are verification of the declared kinematic
map/induction law; they do not prove momentum/energy closure.

## Explicit balance diagnostics and remaining physical limit

The canonical background is stationary in its rotating magnetic frame. For
`Vrot = Omega cross X`, Eulerian derivatives obey

```
scalar_t = -grad(scalar) dot Vrot
vector_t = Omega cross vector - grad(vector) Vrot.
```

The core reports the actual continuity and ideal-induction residuals and the
sustaining force required by its prescribed flow:

```
f_required = rho DU/Dt + grad(p) - rho g - (curl B) cross B / mu0.
```

The field-aligned component divided by density is explicitly available in
m/s^2, as is mechanical power `U dot f_required` in W/m^3. The internal-energy
source needed to retain the specified isothermal state is

```
Q_required = (Dp/Dt + gamma p div(U)) / (gamma - 1).
```

When the declared mass continuity holds this reduces to `p div(U)`. These
quantities diagnose a missing sustaining force/heat authority; they are not
silently applied as a newly approved physical heating model.

In particular, reusing a radial Parker speed magnitude along a nonradial PFSS
tube need not satisfy the field-aligned nozzle momentum equation with that
tube's actual expansion. The maintained manufactured nonradial test exposes a
finite required field-aligned force and compares it against independent
material/time and pressure differencing. It must not be described as an
observational prediction or a force-balanced coronal wind. Simultaneous coronal
density/measured-wind qualification still needs the user's actual assets and,
if this candidate fails them, a selected flux-tube wind/heating closure.
No such observational data or closure is fabricated by these tests.

## Spherical front and line connections

`swcme_spherical_geometry.hpp` contains the exact canonical sphere relation
previously in SWCME3D: radius equals the front radius, outward normal equals
the normalized heliocentric ray and normal speed equals the canonical front
speed. The legacy SWCME3D sphere branch delegates to it, preserving its
normalization/radius arithmetic. The new event therefore acquires no irrelevant
legacy Leblanc, sheath or source parameters merely to use sphere geometry.

`Surface` evaluates the actual same-epoch ambient at the physical front. Fast
Mach uses `V_front,n - U dot n` and the local oblique fast speed. A sub-fast
geometric front has no jump; a fast point must pass the unchanged canonical
RH compression/entropy/residual gates. A physical launch at 1.05 Rs can remain
sub-fast. Solver nonconvergence propagates and never becomes a compression
floor or gas-dynamic fallback.

`Connect` intersects actual supplied same-epoch polyline segments with the
same spherical front, retaining both crossings and tangencies. Long-double
quadratic coefficients start from the exact input doubles. A negative
discriminant remains a miss even when close to zero; a spatial residual
tolerance does not promote it to a tangent. Shared endpoints use the next
segment at fraction zero, matching native field-line ownership. Surface
residuals are graded in metres. The supplied polyline represents the existing
field-line geometry; the method does not add solar-wind velocity to relative
particle displacement or remap a line to another epoch implicitly.

## Verification and use

The explicitly manufactured, shared test fixture lives in
`test/spherical_shock_only_fixture.h`. It inherits the Stage 14 dipole and
isothermal control anchors and is labelled test input only. Its 8/10-AU guard
controls do not change the physical 5-AU sphere and are not scientific defaults.
The original ten `S13SC01`--`S13SC10` core cases cover admission/identity,
Parker/EOS invariants, full material deformation and inverse maps, nonradial
mass/induction/momentum/heat, source-surface continuity, local shock physics,
all sphere roots, unchanged fixed dt/guard semantics, existing front laws and
both magnetic polarities. No old case is removed or reclassified.

```sh
make -C src/models/sep_corona_swcme -j8 test-spherical-shock-only
src/models/sep_corona_swcme/build/test_spherical_shock_only --test S13SC07
python3 -B src/models/sep_corona_swcme/test/check_architecture.py
```

The focused executable is also a dependency/execution of this model's local
`make test`; the srcSEP validation owner runs the ten IDs in its maintained
routine. Native material-profile admission, source/scattering calibration,
physical escape/loss callbacks and observation comparison require their own
application/data gates. The exact scope here remains the shared spherical
background/material/front candidate.

### Independent matched 1-D/3-D spherical verification

The added test-only `test/spherical_limit_reference.h` resolves exactly the
same controlled isothermal PFSS-polar/Parker laws, without calling production
wind, EOS, PFSS, derivative, trajectory, shock or critical-selection routines.
The polar ray of the declared flux-balanced dipole is open; it does not require
or introduce a magnetic monopole. Analytic dipole flux-tube area below the
source surface and radial flux above it determine the same mass loading.
A long-double scalar Parker invariant solve selects the unique transonic
radial branch, independently of any sampled tube critical-point algorithm.

The reference predicts full local Cartesian tensors, including transverse
tube expansion and rotation. Independent scalar travel-time quadrature and
radius inversion determine material position; flux-tube area, radial travel
time sensitivity and integrated angular rotation determine the full F and J.
The regular parallel fast-shock limit has an independent analytic compression,
pressure, velocity and signed normal field, with its fast-branch applicability
checked explicitly. Both hemispheres and both dipole polarities are tested.

These checks extend the existing S13SC02, S13SC03 and S13SC06 records, preserve
every original assertion and add no new scientific profile. The Stage 12
adiabatic, constant-speed profile remains separate and unchanged. Agreement
here verifies the declared isothermal kinematic candidate; it does not solve
its outstanding force/heating closure or turn it into observed-wind evidence.
The existing Cartesian PFSS pole uses its documented finite angular limit;
normalized vector/tensor comparisons include that numerical error rather than
assuming an exact zero transverse component. Added limits are fixed in the
independent reference header before execution and scale tiny SI quantities
with their physical vector/tensor magnitude.
