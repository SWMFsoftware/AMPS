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

The current shock correction deliberately does **not** repair the independent
DBM slow-CME/zero-drag issues, the 1-D ejecta density/velocity-factor bugs, the
shock-mesh apex/seam topology, centralized configuration/unit validation, or
magnetic connectivity/cobpoint tracking.  Those remain separate remediation
items with their own validation gates.
