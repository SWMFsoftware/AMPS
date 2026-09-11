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


## Shared common physics core for 1-D and 3-D

SWCME now prepares dimensionality-independent solar-wind and apex-kinematic
physics through two shared production components:

- `swcme_solarwind.hpp` owns the Leblanc density coefficients and normalization,
  Parker radial-field normalization, Parker scalar components, the Cartesian
  Parker vector for an arbitrary solar axis, and the current proton thermal
  pressure closure; and
- `swcme_core.hpp` converts one common public-unit configuration to SI, prepares
  the shared solar-wind cache, and evaluates the shared ballistic/DBM/data-driven
  apex kinematics.

Both `swcme1d::StepState` and `swcme3d::StepState` retain legacy mirror fields
(`C2/C4/C6`, `Br1AU_T`, `k_AU`, `V_sw`, apex radius/speed) so existing source
continues to compile.  Those fields are no longer independently calculated; they
are copied from `StepState::common`.  This makes the common prepared state the
authoritative source while preserving the current public interfaces.

The dimensional wrappers now differ only where geometry genuinely differs.  The
1-D model supplies a fixed ray `sin(theta)` to the common Parker-component
function.  The 3-D model derives local latitude from the solar rotation axis and
radial direction and asks the same common Parker implementation to construct the
Cartesian field.  Both models call the same common Leblanc density evaluator and
the same proton-pressure closure before entering the already shared ideal-MHD
shock solver.

The new `1D3D01` and `1D3D02` validation tests are permanent guards against a
return of duplicated physics.  `1D3D01` compares the canonical common cache and
the public upstream density/velocity/Parker field for an exactly equivalent
equatorial geometry.  `1D3D02` compares the complete MHD shock state in the
spherical +X limit, including shock existence, Mach number, compression,
upstream/downstream vectors, density, and pressure.  Both must agree to
roundoff-level tolerances.

This refactor intentionally does **not** unify the phenomenological sheath/ejecta
region shaping; that behavior is still dimension-specific and is covered by the
separate region-model remediation item.  Likewise, `1D3D03` remains reserved for
the future SEP-source contract, which has not yet been implemented.

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

### Surface-owned shock state and query-point invariance

The shock state is now explicitly **owned by the physical shock surface**.  The
Cartesian field evaluators, connectivity solver, directional diagnostics, and
shock-mesh builder all obtain local shock properties through the same
`shock_state_direction()` path.  In particular, ambient density and magnetic
field used to form the Mach number are sampled at `R_shock * u`, not at the
radius of a point where a caller happens to request `n`, `V`, or `B`.  This
prevents the same physical shock from acquiring different compression ratios
when queried from different upstream/downstream locations.

`local_oblique_rc()` remains only as a source-compatible scalar wrapper for
older callers.  Its historical `r_eval_m`, `Rdir_m`, and `n_hat` inputs no
longer control shock physics; the wrapper recomputes the canonical surface state
from the direction and returns scalar projections of that state.  New production
code should call `shock_state_direction()` directly.

The mesh builder and `diagnose_direction()` were also converted to consume the
canonical `LocalShockState` directly, eliminating the last internal secondary
shock-strength path.  Validation tests `SHK13` and `SHK14` permanently guard
query-radius invariance and mesh/diagnostic consistency.

The canonical 3-D shock API retains the exact MHD upstream/downstream state at
the mathematical shock surface.  Transport-facing FULL_ICME fields may represent
that discontinuity with the finite C1 shock layer selected by
`RESOLVED_COMPRESSION`; the inner edge of that numerical layer is pinned to the
exact RH downstream state.  SOURCE mode instead uses SHOCK_ONLY and therefore
contains no resolved RH compression in the transport background.  The interior
sheath relaxation remains phenomenological and separate from the exact shock
diagnostic state.

The 1-D model uses the same shared ideal-MHD jump solver.  Its radial direction
is the shock normal and the Parker azimuthal field is tangential to that normal.
This removes the former 1-D compression-floor shock and makes the downstream
normal speed satisfy the same physical jump conditions as 3-D.

The shock validation suite `SHK01`-`SHK14` covers shock/no-shock classification,
obliquity, independent parallel and perpendicular limiting solutions, oblique
branch continuity, mass flux, normal magnetic field, tangential electric field,
momentum and energy fluxes, entropy/admissibility, the weak-shock limit, query-radius invariance, and canonical-state consistency across diagnostics/mesh outputs.

## Remaining remediation items

The dimensionality-independent constants, units, configuration validation,
Leblanc/Parker ambient state, apex kinematics, and ideal-MHD shock solver are now
shared between 1-D and 3-D.  Remaining work is intentionally focused on infrastructure or physics outside
the now-shared region and acceleration models: shock-mesh apex/seam topology,
explicit numerical-status propagation in the remaining 3-D field evaluators,
the full AMPS-facing SEP source adapter/units contract, velocity-divergence
cleanup, and the higher-level Python validation campaign runner.

## Repaired sheath/ejecta region model and SHOCK_ONLY/FULL_ICME modes

The phenomenological downstream-region model is now shared through
`swcme_regions.hpp`.  This closes the former 1-D/3-D divergence in layer
geometry, ejecta factors, and artificial leading/trailing-edge smoothing.

Two explicit modes are available through `Params::region_mode`:

- `swcme::regions::Mode::ShockOnly` returns the undisturbed analytical
  Parker/Leblanc solar-wind state everywhere.  Shock geometry, connectivity,
  and local MHD shock/source diagnostics remain available through their
  dedicated APIs, but no sheath/ejecta plasma modification is applied to the
  transport-facing background.  This is the recommended controlled baseline
  for the SEP connectivity/perpendicular-diffusion study.
- `swcme::regions::Mode::FullICME` adds the optional phenomenological sheath and
  magnetic-ejecta profile behind the geometric CME/shock surface.

The exact physical shock remains an explicit surface in the diagnostic API and
its downstream state is the ideal-MHD Rankine-Hugoniot solution.  For transport,
FULL_ICME is paired with `RESOLVED_COMPRESSION`: a symmetric finite-width C1
layer is centered on the mathematical shock, reaches the analytical upstream
state at its outer edge, and reaches the exact RH downstream state at its inner
edge.  The sheath then relaxes smoothly toward a leading-edge target.  The
empirical `sheath_comp_floor` no longer participates in shock or region physics
and is retained only for source compatibility.

### Self-similar local layer geometry

The public sheath/ejecta thickness inputs are specified as AU at a 1-AU shock.
They are now interpreted as dimensionless self-similar fractions.  For a local
shock-surface radius `R_sh(u)`,

```text
R_LE = (1 - f_sheath) R_sh
R_TE = (1 - f_sheath - f_ejecta) R_sh.
```

This is applied to the **local** SSE/ellipsoid radius, not the apex radius.
Consequently the shock, leading-edge, and trailing-edge surfaces remain nested
and geometrically similar from apex to flank.  Configurations with
`f_sheath + f_ejecta >= 1` are rejected instead of repaired by runtime clipping.

The LE and TE smoothing inputs are likewise local self-similar fractions.  Each
configured width is the total width of a symmetric C1 smoothstep transition
centered on the nominal boundary and is capped so neighboring finite layers
cannot overlap.

### Sheath and magnetic-ejecta targets

For a physical fast shock, the sheath starts at the complete MHD downstream
state and relaxes toward the local Parker/upstream state at the leading edge.
If a geometric CME front exists but the fast shock has decayed, the model does
not fabricate a sheath compression; that portion of the profile remains
ambient while the optional ejecta may still be represented.

Magnetic-ejecta factors are now honored exactly:

```text
n_ME = f_ME * n_up
V_ME = V_ME_factor * V_sw.
```

Values below unity therefore correctly produce density depletion and slower
bulk ejecta.  Negative factors are rejected by centralized configuration
validation rather than clipped in the evaluator.  The baseline ejecta magnetic
field remains Parker; a flux-rope/ejecta-field model is intentionally outside
the controlled one-year scope.

Artificial LE and TE interfaces are C1 smooth.  The numerical shock layer is
also C1, but only when `RESOLVED_COMPRESSION` is selected.  Its total width is
`edge_smooth_shock_AU_at1AU * R_sh(local)` and is capped so it cannot consume the
entire sheath.  In SOURCE mode that shock-layer width is forced to zero and the
SHOCK_ONLY background remains analytical on both sides of the mathematical
source surface.

Validation tests `REG01`-`REG05` cover SHOCK_ONLY identity, the exact RH state at
the inner edge of the resolved shock layer, sub-unity ejecta factors, local
self-similar surface nesting, and continuity/smoothness of the modeled region
transitions.

See `REGION_MODEL_FIX_NOTES.md` and `test/README.md` for detailed equations,
configuration conventions, and validation procedures.


## Single shock-acceleration representation

`swcme_acceleration.hpp` defines one authoritative acceleration switch shared by
1-D and 3-D:

```cpp
swcme::acceleration::Mode::Source
swcme::acceleration::Mode::ResolvedCompression
```

The model deliberately does **not** expose independent `use_dsa_source` and
`use_compression_acceleration` booleans.  A single enum prevents a run from
enabling both representations of first-order shock acceleration for the same
particle population.

### SOURCE mode

`SOURCE` is the recommended baseline for the connectivity/perpendicular-
diffusion experiment.  Centralized validation currently requires

```text
acceleration = SOURCE
regions      = SHOCK_ONLY
```

so the transport-facing velocity remains the analytical solar wind through the
mathematical shock surface.  A physical fast shock instead produces a
`ShockAccelerationState` with `source_enabled=true`, local shock position and
normal, normal shock speed, compression, `theta_Bn`, fast Mach number, upstream
density and magnetic-field magnitude, and the diagnostic DSA phase-space slope

```text
q = 3 r_c / (r_c - 1).
```

The baseline `relative_source_weight_per_area` is dimensionless and is intended
only for controlled relative weighting.  Absolute particle/spectral units are
reserved for the later AMPS-facing source adapter so this layer does not invent
an injection-rate convention prematurely.

### RESOLVED_COMPRESSION mode

Centralized validation currently requires

```text
acceleration = RESOLVED_COMPRESSION
regions      = FULL_ICME
edge_smooth_shock_AU_at1AU > 0
```

and no prescribed DSA source is exposed.  The exact RH discontinuity remains in
`shock_state_direction()` / the 1-D cached jump for diagnostics, while the
transport-facing primitive fields use one common symmetric C1 shock profile.
For a total local width `w_sh`, the profile spans

```text
R_sh + w_sh/2   upstream endpoint
R_sh            midpoint
R_sh - w_sh/2   exact RH downstream endpoint.
```

The smoothstep derivative vanishes at both endpoints.  The phenomenological
sheath begins at the inner endpoint, and its own profile also has zero slope
there, so there is no second numerical kink at the resolved-shock/sheath join.
The same `swcme_regions.hpp` formulas are called by 1-D and 3-D, including at
finite-SSE/ellipsoidal flanks where `w_sh` scales with the **local** shock radius.

In `RESOLVED_COMPRESSION`, `ShockAccelerationState::source_enabled` is false,
`resolved_compression_enabled` is true, the active DSA slope is intentionally
reported as unavailable, and the source weight is zero.  This prevents a
consumer from silently applying a pre-imposed DSA source on top of the resolved
`div(V)` accelerator.

The deterministic helper `swcme::acceleration::serialize_csv()` exists for
regression/audit comparisons.  It is not the final science-data format.  Tests
`ACC01`-`ACC05` and `1D3D03` verify the no-double-counting gate, C1 resolved shock
profile, RH endpoints, 1-D/3-D smoothing identity, disabled source in resolved
mode, and byte-identical SOURCE records in the spherical radial limit.

See `SHOCK_ACCELERATION_FIX_NOTES.md` and `test/README.md` for the validation
fixtures and implementation details.

## Explicit numerical status and error propagation

`swcme_status.hpp` provides the shared `swcme::ModelStatus` contract used by
numerically sensitive 1-D/3-D evaluation paths.  The production physics layer
no longer converts an invalid intermediate state into zero, ambient plasma, a
clipped radius, or an arbitrary +X direction merely to keep a calculation
running.

New code should prefer the checked APIs:

```cpp
swcme1d::Model::evaluate_radii_fast_checked(...)
swcme1d::Model::evaluate_radii_with_B_div_checked(...)

swcme3d::Model::shock_state_direction_checked(...)
swcme3d::Model::evaluate_cartesian_fast_checked(...)
swcme3d::Model::evaluate_cartesian_with_B_checked(...)
swcme3d::Model::evaluate_cartesian_with_B_div_checked(...)
swcme3d::Model::compute_divV_radial_checked(...)
```

They return `ModelStatus` with an explicit code, context, and failed sample
index for batch queries.  The existing void/bool science wrappers remain for
source compatibility but convert numerical failures to exceptions instead of
silently repairing the result.  `NO_SURFACE` remains an expected geometrical
outcome for directions outside a finite front and is distinguishable from an
actual numerical failure.

The analytical Leblanc/Parker lower boundary, `1.05 R_sun`, is now a real model
domain: science queries below it return `OUTSIDE_MODEL_DOMAIN`.  The low-level
solar-wind equations no longer clip the caller's radius.  Vector norms use
`std::hypot` and degenerate/non-finite vectors are rejected; the previous +X
normalization fallback has been removed.

The shared ideal-MHD jump result also carries an explicit `SolveStatus`
(`NO_SHOCK`, `SOLVED`, `INVALID_INPUT`, `NO_PHYSICAL_BRACKET`,
`INVALID_ACCEPTED_STATE`, or `CONSERVATION_FAILURE`).  A failed super-fast RH
solve is propagated as `SHOCK_SOLVER_FAILURE`; it is not converted to a
compression-one ambient state.

Mesh and checked Tecplot output paths reject non-finite physics instead of
writing a sanitized surrogate.  `ERR01`-`ERR05` protect outside-domain behavior,
non-finite batch input, explicit RH outcomes, degenerate-vector rejection, and
writer rejection of corrupt mesh data.  See `NUMERICAL_STATUS_FIX_NOTES.md` and
`test/README.md` for the complete status semantics and validation fixtures.
