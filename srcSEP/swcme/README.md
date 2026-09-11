# SWCME semi-analytical solar-wind / CME-shock model

SWCME provides lightweight 1-D and 3-D solar-wind/CME-shock backgrounds for
energetic-particle transport studies.  The 3-D model combines a Parker magnetic
field, Leblanc density profile, analytical CME/shock kinematics, configurable
shock geometry, and phenomenological sheath/ejecta fields.

## Current 3-D shock geometries

`swcme3d::ShockShape` currently provides:

- `Sphere` — a Sun-centered spherical verification geometry;
- `Ellipsoid` — a Sun-centered ellipsoid whose three axes scale self-similarly
  with the shock apex distance; and
- `SSE` — the recommended finite self-similar-expansion spherical cap.

`ConeSSE` is retained as an enum alias for source compatibility but now has the
same corrected semantics as `SSE`.  It no longer means the old
`R=R_apex*cos(theta)^m` cosine-cap model.  The legacy `flank_slowdown_m`
parameter is retained in `Params` only so existing source/input code continues
to compile; it is ignored for the corrected SSE geometry.

### Finite SSE geometry

For apex distance `R_apex` and half width `lambda`, the generating sphere is

```text
c = R_apex / (1 + sin(lambda))
a = c sin(lambda),
```

where `c` is the heliocentric distance of the sphere center along the CME axis
and `a` is its radius.  A ray separated by `alpha` from the CME axis intersects
the outward cap at

```text
R(alpha) = c cos(alpha) + sqrt(a^2 - c^2 sin(alpha)^2)
```

for `alpha <= lambda`.  The ray is tangent at `alpha=lambda`; no surface exists
for larger angular separation.  The exact outward normal is

```text
n_hat = (R e_r - c e_CME) / a.
```

`Model::shape_radius_normal()` returns a boolean surface-existence flag.  When
it returns `false`, the radius and normal outputs are zero and are not physical
values.  `Model::diagnose_direction()` follows the same convention and returns
`rc=1`, `V_sh,n=0` for an absent finite surface.  The Cartesian evaluators use
this flag internally and return the undisturbed ambient wind/Parker field
outside the SSE angular support.

### Self-similar flank speed

All supported shapes use the same geometric normal-speed relation.  If the
shape scales with apex distance so that `R(u,t)=f(u) R_apex(t)`, then

```text
dR/dt   = V_apex R/R_apex
V_sh,n  = V_apex (R/R_apex) (e_r dot n_hat).
```

This replaces the legacy independent cosine flank-speed factor and corrects the
ellipsoid flank speed, which previously used the full apex speed away from the
apex.


## Shared CME/shock-apex kinematics

The 1-D and 3-D models now obtain the shock-apex radius and speed from the same
`swcme_kinematics.hpp` implementation.  This removes the former divergence in
which the 1-D slow-CME branch clipped `V0-Vsw` to zero while the 3-D branch used
the fast-CME formula for negative speed differences and divided by `Gamma` at
`Gamma=0`.

Three kinematic modes are supported:

- `swcme::kinematics::Mode::Ballistic` uses
  `R=R0+V0*t`, `V=V0` exactly;
- `Mode::DBM` uses the sign-aware constant-background drag-based solution; and
- `Mode::DataDriven` uses monotone PCHIP interpolation of a supplied
  height-time table and returns the derivative of that same interpolant as the
  apex speed.

For DBM, define `DeltaV0=V0-Vsw` and `a=abs(DeltaV0)`.  The common solution is

```text
DeltaV(t) = DeltaV0 / (1 + Gamma a t)
R(t)      = R0 + Vsw t
            + sign(DeltaV0) log(1 + Gamma a t) / Gamma.
```

The absolute value in the drag denominator is essential: fast CMEs decelerate
toward the ambient wind and slow CMEs accelerate toward it.  `Gamma=0` is an
explicit ballistic branch, so there is no division by zero and no physics-level
`finite_or()` fallback.  A short series for `log1p(x)/x` is used at very small
`x=Gamma*a*t` so the nonzero-drag solution remains continuous with the exact
ballistic limit.

### Data-driven kinematics

`Params::data_time_s` contains strictly increasing times in seconds and
`Params::data_radius_Rs` contains nondecreasing apex radii in nominal solar
radii.  PCHIP was selected because it passes exactly through the supplied
height-time knots while preserving monotonicity and avoiding cubic overshoot.
The local derivative is used as `V_sh`, so the reported speed is kinematically
consistent with the radius curve.

The default extrapolation policy is
`swcme::kinematics::ExtrapolationPolicy::OutsideTime`: a query before the first
or after the last knot is rejected explicitly.  `Ballistic` continuation may be
selected when an explicit endpoint continuation assumption is desired; it uses
the endpoint PCHIP derivative and does not silently extrapolate the cubic.

Both dimensional wrappers convert their public parameters to the same SI
`swcme::kinematics::Config`.  Therefore identical kinematic inputs are required
to produce identical apex radius and speed to roundoff.  The recommended
default DBM reference radius is now `20 R_s`, reflecting the intended use of
this simple drag model in the drag-dominated heliosphere rather than at
`~1.05 R_s`.  Event-specific calculations may still set another radius when
there is a documented physical justification; observationally constrained
`DataDriven` mode is preferred when height-time measurements are available.

The deterministic kinematics validation block is `KIN01`-`KIN08`; see
`test/README.md` for individual purposes and acceptance criteria.

## Observer-to-shock magnetic connectivity and cobpoint tracking

The 3-D model now provides an analytical Parker-field-line connectivity solver
through `Model::observer_connectivity()`.  The solver is designed for the
controlled upstream-Parker experiment used by the SEP study: it traces the
observer's nominal Parker line inward and intersects that line with the same
production shock geometry used by `shape_radius_normal()` and the same local
shock physics used by `shock_state_direction()`.

For the production Parker field

```text
B_phi/B_r = -Omega r sin(theta) / V_sw,
```

a field-line tangent satisfies

```text
r sin(theta) dphi/dr = B_phi/B_r,
```

so `dphi/dr=-Omega/V_sw`.  The exact observer-anchored field line is therefore
constructed by rotating the observer radial direction about the configured
solar-rotation axis by

```text
Delta phi = -Omega (r-r_obs) / V_sw.
```

`Params::solar_rotation_rate_rad_s` now explicitly carries the rotation rate
used by both the Parker magnetic field and connectivity mapping.  Its default
is the existing SWCME solar-rotation convention.  Setting it to zero provides
the exact radial-field limit used by `CON01`; the normal production default is
unchanged.

`ConnectivityState` retains every geometrical field-line/shock intersection in
increasing radial order.  Each `ConnectivityRoot` contains the Cartesian
cobpoint, shock radius residual, analytical Parker path length to the observer,
and the complete `LocalShockState`.  The default selected cobpoint is the
outermost root, i.e. the first shock surface encountered when tracing inward
from the observer.  Retaining all roots makes this choice explicit and keeps
the infrastructure usable for future non-convex geometries.

The root search is deliberately robust to connection boundaries.  It combines
radial scanning, bisection of sign-changing roots, local minimization of the
surface residual to detect tangent roots that do not change sign, and explicit
refinement of finite-SSE surface-validity transitions.  No artificial SSE
flank is introduced when the Parker line remains outside the configured cap.

The Parker path length returned with a cobpoint is analytical.  With

```text
k = Omega sin(theta) / V_sw,
```

SWCME integrates

```text
ds/dr = sqrt(1 + (k r)^2)
```

in closed form.  The exact zero-rotation/polar limit is the radial distance.
This length, rather than radial separation, is the quantity intended for
field-aligned SEP transport timing.

`Model::observer_connectivity_history()` evaluates a stationary observer at a
requested set of times.  Every time step is solved independently from the
production kinematics, Parker line, geometry, and shock state; the history
contains no hidden hysteresis.  Consequently connection onset/loss and
cobpoint motion can be interpreted as model physics rather than state retained
by the tracker.

The deterministic connectivity validation block is `CON01`-`CON08`; see
`test/README.md` for the individual fixtures and acceptance checks.

## Centralized configuration validation and unit handling

SWCME now separates **unit conversion** from **physical admissibility**.
`swcme_units.hpp` is the single production source for conversions between the
public heliophysics units and SI, including km/s, nT, cm^-3, km^-1, AU, solar
radii, hours, and degrees.  Conversion functions are deliberately pure: for
example, `0 km/s` converts to `0 m/s`; it is not silently replaced by a
positive speed.

`swcme_config.hpp` provides the common validation contract used by both the
1-D and 3-D models.  `Model::validate()` is side-effect free and returns every
invalid field with a structured code and requirement.  `prepare_step()` calls
the same validator before any basis normalization, unit conversion, kinematic
evaluation, or field/shock calculation.  Invalid input therefore fails once at
setup rather than being clipped into a plausible-looking state.

Common checks include positive solar-wind speed, reference density and
temperature; non-negative magnetic field and DBM drag coefficient; `gamma>1`;
`sin(theta)` in `[0,1]`; valid region/smoothing parameters; and well-formed
DATA_DRIVEN tables.  The 3-D interface additionally validates non-zero finite
CME and solar-rotation axes, positive ellipsoid axis ratios, non-negative solar
rotation rate, and SSE half width in `(0,pi/2]`.

The production models now use the centralized conversion helpers when preparing
their SI state.  This fixes the prior 1-D `V_sw` path that used
`max(1,V_sw*1000)` and therefore changed `0 km/s` into `1 m/s`.  Zero wind speed
is now converted exactly and rejected by the configuration validator because a
positive wind speed is required by the Parker/DBM baseline.

See `CONFIGURATION_UNITS_FIX_NOTES.md` and the `CFG01`/`CFG02` sections in
`test/README.md` for the full contract and validation coverage.

## Validation

Build the validation executable from `test/`:

```sh
cd test
make
./output/test_swcme --list
./output/test_swcme --all
```

The finite shock geometry is covered by `GEO01`-`GEO08` in
`test/3d/test_geometry.cpp`.  See `test/README.md` for detailed test purposes,
reference calculations, and acceptance criteria.

## Ideal-MHD shock existence and downstream state

The production model now separates the existence of a geometric CME front from
the existence of a physical fast shock.  A local surface point is classified as
a shock only when its normal speed exceeds the upstream normal flow by more
than the local oblique fast-mode speed.  If that criterion is not met, the
model returns `has_shock=false`, `compression=1`, and an unchanged downstream
state.  The legacy `sheath_comp_floor` parameter no longer changes the physical
shock compression and therefore cannot manufacture a shock.

The shared `swcme_shock.hpp` solver evaluates the ideal-MHD Rankine-Hugoniot
conditions in the shock frame.  For a trial density compression it enforces
mass conservation, tangential momentum conservation, and tangential electric
field continuity; normal momentum gives the downstream pressure and a
bracketed scalar solve enforces total-energy-flux conservation.  Accepted
solutions must be compressive, have positive downstream pressure, increase the
entropy proxy, and satisfy the stored conservation residual tolerances.

`swcme3d::Model::shock_state_direction()` is the preferred 3-D API for local
shock diagnostics.  It returns the shock-surface radius and normal, normal shock
speed, fast Mach number, `theta_Bn`, density compression, complete upstream and
downstream primitive states, and conservation residuals.  Upstream quantities
are always evaluated at the actual shock surface rather than at an arbitrary
query radius.

The 3-D Cartesian field evaluators now use the exact MHD downstream state in
the limit immediately behind the shock.  The shock itself is treated as a
physical discontinuity: points at or ahead of the surface return the upstream
state, while the limit just behind the surface returns the RH downstream state.
The interior sheath relaxation remains phenomenological and is a separate model
component.

The 1-D model uses the same shared ideal-MHD jump solver.  Its radial direction
is the shock normal and the Parker azimuthal field is tangential to that normal.
This removes the former 1-D compression-floor shock and makes the downstream
normal speed satisfy the same physical jump conditions as 3-D.

The shock validation suite `SHK01`-`SHK12` covers shock/no-shock classification,
obliquity, independent parallel and perpendicular limiting solutions, oblique
branch continuity, mass flux, normal magnetic field, tangential electric field,
momentum and energy fluxes, entropy/admissibility, and the weak-shock limit.

## Remaining remediation items

The current model still requires separate remediation of the 1-D ejecta
density/velocity-factor bugs, shock-mesh apex/seam topology, centralized
configuration/unit validation.  Magnetic connectivity/cobpoint tracking is now
implemented and validated by CON01-CON08; the remaining items retain independent
validation gates.  CME/shock-apex kinematics are now
shared and validated by KIN01-KIN08.
