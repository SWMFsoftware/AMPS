# SWCME finite shock geometry correction

This update replaces the former `ConeSSE` cosine-cap approximation with a true
finite self-similar-expansion (SSE) spherical-cap geometry and corrects the
normal speed for every supported self-similar shape.

## Physics change

For a CME/shock apex distance `R_apex` and angular half width `lambda`, the SSE
front is represented by the outward portion of a sphere centered on the CME
axis.  The center distance `c` and sphere radius `a` are

```text
c = R_apex / (1 + sin(lambda))
a = c sin(lambda).
```

The outward heliocentric ray at angular separation `alpha` from the CME axis
intersects the front at

```text
R = c cos(alpha) + sqrt(a^2 - c^2 sin(alpha)^2)
```

for `alpha <= lambda`.  The discriminant vanishes at `alpha=lambda`, where the
ray is tangent to the front.  No surface exists beyond that angular extent.

The exact local outward normal is

```text
n_hat = (R e_r - c e_CME)/a.
```

Because the front is self-similar, a point at fixed angular coordinates has
radial velocity

```text
dR/dt = V_apex R/R_apex
```

and the shock-normal speed is

```text
V_sh,n = V_apex (R/R_apex) (e_r dot n_hat).
```

The same self-similar speed expression is now used for the sphere and
ellipsoid, fixing the previous ellipsoid flank-speed error.

## API behavior

`ShockShape::SSE` is the corrected geometry.  `ShockShape::ConeSSE` remains as
an enum alias for source compatibility but now has the same SSE semantics.  The
legacy parameter `flank_slowdown_m` is retained in `Params` so older source/input
code still compiles, but it is intentionally ignored for SSE because flank
radius and speed are now determined by the geometry itself.

`Model::shape_radius_normal()` returns `bool`:

- `true`: a physical surface intersects the query ray;
- `false`: no finite SSE surface exists in that direction.

On `false`, radius and normal outputs are zero.  `diagnose_direction()` follows
the same convention and returns `rc=1` and `V_sh,n=0`.  Cartesian field/plasma
evaluators explicitly return the undisturbed ambient state outside the SSE
angular support.

## Validation

Tests `GEO01` through `GEO08` were added under `test/3d/test_geometry.cpp` to
cover the spherical reference, SSE apex/tangent/finite extent, level-set
residual, normals, finite-difference normal speed, ellipsoid speed, and
rotational covariance.

## Scope not changed by this update

This update intentionally does not modify the existing compression-ratio proxy,
Rankine-Hugoniot physics, immediate downstream sheath speed, DBM sign/zero-drag
handling, or shock-mesh apex topology.  Those are independent remediation items
and should be changed only with their own validation tests.
