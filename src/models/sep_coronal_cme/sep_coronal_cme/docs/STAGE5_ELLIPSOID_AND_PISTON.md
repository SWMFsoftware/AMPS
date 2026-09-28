# Stage 5: candidate-front ellipsoid and independent piston

Stage 5 implements the fixed-orientation schema-5 geometry in
`ellipsoid_geometry.{h,cpp}`.  It deliberately does not add the future
deflection/rotation capability: the center remains on one declared radial
direction and the right-handed principal basis is fixed while center distance
and all three semiaxes evolve independently.

## Principal basis and parameterization

`BuildRadialPrincipalBasis` constructs `e_r`, local east, and local north from
the configured latitude/longitude, then applies the declared right-handed
tilt about `e_r`.  It verifies `e_r cross e_1=e_2`; a handedness failure is a
configuration error, not an axis relabeling.

The surface uses

`F=(x-c)^T Q (x-c)-1`,
`c=d_c e_r`,
`Q=R diag(a_r^-2,a_1^-2,a_2^-2) R^T`.

Exactly one radial form is supplied to a factory.  `FromCenter` consumes
`d_c`; `FromApex` consumes `r_apex` and derives the complete center history
componentwise from `d_c=r_apex-a_r`, including first and second derivatives.
Both produce the same immutable `FixedOrientationEllipsoid`.

## Independent kinematics

`EvaluateSmoothKinematics` implements the Section-8.4 quintic rate transition.
It continues the initial rate to the component's own transition time, blends
to the final rate with zero endpoint acceleration, integrates that rate
analytically, and then continues at the final rate.  Each of `d_c`, `a_r`,
`a_1`, and `a_2` can use independent transition parameters.

`EvaluateCubicHermiteHistory` is the tabulated `C1` alternative.  Every knot
contains value and rate, time ordering is strict, derivatives are analytic,
and extrapolation is rejected.

## Normal, speed, curvature, and quadrature

The normal is `Q(x-c)/|Q(x-c)|`.  The exact implicit normal speed is

`V_n=-F_t/|grad F|`,
`F_t=-2 c_dot dot Q(x-c)+(x-c)^T Q_dot (x-c)`.

For the fixed basis, each diagonal component of `Q_dot` is
`-2 a_dot/a^3`.  This retains separate translation and expansion and prevents
the common error of using apex speed everywhere on the flank.  Mean and
Gaussian curvature are evaluated analytically from `Q`; the sphere limit is
`H=1/a` and `K=1/a^2`.

`Tessellate` uses midpoint surface quadrature in ellipsoidal polar/azimuthal
coordinates.  The exact parameter Jacobian supplies patch area.  IDs are the
global parameter-cell indices, so solar clipping leaves gaps instead of
renumbering patches and MPI decomposition cannot change physical identity.
Only patch centers on or above the configured solar sphere are published.
`PatchSpatialIndex` provides an immutable ID query; an AMPS adapter can replace
its linear scan with a distributed acceleration structure without redefining
the surface.

## Piston nesting

The piston is a second, independently configured ellipsoid.  No hidden
standoff is inferred.  `CheckPistonNesting` first requires every sampled piston
point to lie inside the front.  It then intersects the front with the same
heliocentric ray through that piston point and checks the requested minimum
radial separation.  Comparing equal ellipsoid parameters would be wrong for
a displaced center because those points need not lie on one solar ray.

## Validation

`make test-stage5` runs every earlier gate plus `ELL3D01--10`.  The tests cover
sphere/triaxial normals, exact normal speed, convergent sphere area, solar-dome
clipping, stable IDs/spatial lookup, center/apex equivalence, independent
smooth histories, tabulated `C1` knots, basis/Q-dot/curvature regressions, and
valid/invalid independently specified pistons.  `ELL3D11` remains reserved for
the future dynamic-attitude and center-deflection capability.
