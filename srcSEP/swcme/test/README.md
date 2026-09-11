# SWCME validation

The validation suite is a standalone C++ executable that calls production
SWCME interfaces. It uses the repository's existing Make-based build approach;
no second build system or external test dependency is required.

## Directory layout

```text
test/
  core/       Common registry, runner, reporting, CFG/DEN/KIN/SHK tests
  1d/         Tests specific to the production 1-D model
  3d/         Tests specific to the production 3-D model
  reference/  Reviewed reference data for future comparison tests
  profiles/   Inputs describing future validation sampling profiles
  output/     Generated executable and results (ignored by Git)
```

Every validation is linked into the single `output/test_swcme` executable.

## Build

From `srcSEP/swcme`:

```sh
make -C test clean all
```

From `srcSEP/swcme/test`, the equivalent command is `make clean all`.

## Command-line interface

```sh
./output/test_swcme                 # run all tests; same as --all
./output/test_swcme --all           # run all tests in registry order
./output/test_swcme --list          # list tests without executing them
./output/test_swcme --test CFG01    # run exactly CFG01
./output/test_swcme --test CFG02    # run exactly CFG02
./output/test_swcme --test DEN01    # run exactly DEN01
./output/test_swcme --help          # usage, options, and examples
```

Options are mutually exclusive. Unknown test IDs, unknown options, a missing
value after `--test`, and combined modes return exit code 2. A completed test
run returns 0 when every mandatory test passes and 1 when any test fails.
Subcheck SKIPs identify unavailable production paths and do not by themselves
fail an otherwise successful test.

The registry in `core/test_swcme.cpp` is the only source for `--list`, `--all`,
and test lookup, so displayed and executed tests cannot silently diverge.

## Validation classifications

- `COMMON`: contracts shared by both models, with both public paths exercised.
- `1D`: behavior specific to `swcme1d`.
- `3D`: behavior specific to `swcme3d`.
- `1D<->3D`: direct equivalence or consistency between the two implementations.

The current registry uses `COMMON` and `1D`; future comparison tests may use
`1D<->3D` when their primary purpose is cross-model equivalence.

## CFG01 status

The intended CFG01 definition is **Configuration rejection and physical-range
validation**, classified `COMMON`. Production SWCME currently has no unified
configuration-validation API and does not return rejection diagnostics. The
1-D model silently clamps some inputs, while the 3-D model follows different
policies. Inventing a new rejection API is outside the present validation task.

CFG01 therefore temporarily retains the pre-existing physical-constants
baseline and exact 1-D/3-D constant-consistency checks. It explicitly reports
the unavailable configuration-rejection path as SKIP. The unit conversions
formerly mixed into CFG01 moved to CFG02, so no existing coverage was removed.
A future production configuration contract should allow CFG01 to be completed
and its temporary constants checks to move to a dedicated constants test.

## CFG02: unit and dimensional consistency

CFG02 is a `COMMON` test of the complete unit contract currently reachable
through production SWCME. Its purpose is to detect scale errors that can look
numerically plausible but substantially change CME kinematics, density,
magnetic field, Alfvén speed, and shock Mach number.

CFG02 covers:

- forward velocity conversion for 0, 1, 321, 400, 1500, and 3000 km/s;
- forward magnetic-field conversion for 0.1 through 1000 nT;
- forward number-density conversion for 0.01 through 10000 cm^-3;
- AU and nominal-solar-radius distance conversion, including AU/Rs derived from
  the adopted constants rather than an independently rounded ratio;
- the DBM inverse-length conversion from km^-1 to m^-1;
- the 3-D degree-derived cone default at 40, 90, and 180 degrees;
- finite extreme values intended to detect overflow, underflow, integer
  conversion, or hidden scale assumptions rather than event realism;
- explicit 1-D-versus-3-D comparisons for independently implemented paths;
- the closest available equivalent external/SI configuration checks; and
- an Alfvén-speed smoke test using production-converted `B` and `n`, the shared
  production SI Alfvén helper, and a separate scalar reference calculation.

The equivalent configuration uses 400 km/s, 5 nT, 5 cm^-3, 1 AU, 24 hours,
and 40 degrees. SWCME has no unified direct-SI `Params` alternative: radius and
time already enter evaluators in SI, while velocity, field, and density are
heliophysics-unit fields. CFG02 tests the closest real paths and marks the
missing direct-SI configuration comparison as SKIP rather than constructing a
fictitious API.

Likewise, production accepts time in seconds and does not expose an hours input
conversion. The 1-D API accepts angles in radians, and only the 3-D default
contains a degree-to-radian conversion. The 1-D Tecplot writers provide inverse
velocity, AU, and solar-radius output paths; CFG02 parses temporary writer output
to test those round trips at the writer's documented decimal precision. Other
inverse paths are not exposed and are printed as SKIP; test-side arithmetic is
never presented as a production round trip.

CFG02 prints every executed value, units, production result, independent
reference, factor, absolute and relative error, tolerance, and status. It also
prints structured counters and maxima for conversion error, round-trip error,
cross-model difference, and Alfvén-speed error.

### Tolerances

Exact SI decimal-prefix conversions use exact comparison when the selected
binary values permit it. Short floating-point calculation chains use
`64 * std::numeric_limits<double>::epsilon()`. AU and solar radius checks test
the documented adopted SWCME constants, not a competing astronomical
convention. Physical-constant tolerances remain documented in CFG01.

Analytical reference expressions must remain independent of the production
routine under test. A tolerance or reference value must never be relaxed merely
to turn a production mismatch into PASS. A CFG02 failure should be investigated
as a unit-contract or production-policy issue and reported before physics is
changed.

## DEN01: Leblanc density normalization at the reference distance

DEN01 is a `COMMON` physics/implementation test of the defining normalization
property of the production Leblanc, Dulk & Bougeret (1998) density profile. It
does not merely compare the two models with each other: an independent scalar
reference evaluates

```text
n_L(r) = A r^-2 + B r^-4 + C r^-6,
A = 3.3e5, B = 4.1e6, C = 8.0e7,
```

where `r` is heliocentric radius in nominal solar radii and the result is in
`cm^-3`. The test calculates `r_ref_Rs = AU_m / Rs_m` from the adopted constants
rather than maintaining a rounded AU-to-solar-radius constant. For each model,
the intended normalization is `n(r) = S n_L(r)`, with
`S = n_ref / n_L(1 AU)` and a fixed production reference radius of one AU.

The positive reference-density fixtures are 1, 5, 8, 20, and 100 `cm^-3`.
Each fixture constructs a fresh production model and step cache, evaluates the
public 1-D radial or 3-D Cartesian field interface at one AU, and checks

```text
abs(n_model(1 AU) / n_ref - 1) < 1e-12.
```

The acceptance criterion is strict and must not be weakened to turn a mismatch
into PASS. Production evaluators return density in `m^-3`; the test reports in
`cm^-3` for direct comparison with the configured normalization while retaining
the actual SI evaluation path.

The 3-D coverage evaluates `(+AU,0,0)`, `(0,+AU,0)`, `(0,0,+AU)`, and a
normalized non-axis-aligned `(1,2,3)` direction to verify that the ambient
profile is radial. For every reference density, DEN01 also compares 1-D against
3-D density and inferred scale factors at roundoff level. The same reference
location is supplied both as the production AU and as `(AU_m/Rs_m)*Rs_m` to
check AU/solar-radius coordinate equivalence without repeating CFG02.

Both production implementations expose their prepared SI coefficients `C2`,
`C4`, and `C6`. DEN01 independently recovers the nominal `A`, `B`, and `C` from
those caches and verifies the coefficients and inferred normalization factor at
`64 * double epsilon`. That roundoff tolerance covers only the short
floating-point operation chain; it is not a relaxed density-physics tolerance.
The model APIs do not expose a separately named scale factor, so inference from
`C2` is the closest production diagnostic.

The 1-D and 3-D modules currently duplicate the Leblanc coefficient literals
and normalization arithmetic in their respective `prepare_step()` functions.
Both use the shared adopted AU (`149597870700 m`) and nominal solar radius
(`6.957e8 m`), and both cache scaled coefficients in a value-type `StepState`.
No coefficient or normalization factor is static, global, or mutable shared
state. DEN01 nevertheless executes an A=5, B=20, C=1 `cm^-3` construction and
re-evaluation sequence independently for both models. Existing model A and B
must retain their original results after later instances are constructed, and
their cached coefficients must remain unchanged. This guards against future
cross-instance contamination that would be blocking in parallel AMPS use.

The production APIs currently fix density normalization at one AU, so tests of
a configurable 0.5- or 2-AU normalization radius are reported as `SKIP` rather
than supported by a fictitious test-only API. Values at 0.5 and 2 AU are checked
only for finite positivity as optional diagnostics. Detailed verification of
the radial `r^-2 + r^-4 + r^-6` behavior belongs to DEN02, not DEN01.

DEN01 prints the reference radius, each nominal analytical term, unscaled
density, expected and inferred scale, production density, absolute error,
normalization residual, tolerance, direction coordinates, and every individual
PASS/FAIL. Its structured text metrics include fixture/pass/fail/skip counts,
maximum normalization residual, maximum 1-D/3-D difference, coefficient
mutation status, and cross-instance contamination status.

A DEN01 failure may indicate incorrect radial units or coefficient powers, an
incorrect normalization factor, a `cm^-3`/`m^-3` conversion error, inconsistent
1-D and 3-D implementations, mutable global/static normalization state, or
accidentally modified nominal coefficients. The independent analytical helper
must remain separate from the production routines under test; references and
tolerances must never be changed solely to make a production failure pass.

## Current SWCME unit contract

External `Params` values use:

- velocity in km/s;
- magnetic-field magnitude in nT;
- number density in cm^-3;
- launch radius in nominal solar radii;
- drag coefficient in km^-1;
- sheath/ejecta dimensions in AU normalized at 1 AU;
- cone angle in radians, despite the 3-D default being initialized from 40
  degrees; and
- temperature in kelvin.

Model evaluators and prepared states use SI meters, seconds, m/s, tesla, m^-3,
and radians. SWCME adopts `1 AU = 149597870700 m`, nominal solar radius
`6.957e8 m`, CODATA 2022 proton mass `1.67262192595e-27 kg`, vacuum
permeability `1.25663706127e-6 N/A^2`, exact Boltzmann constant
`1.380649e-23 J/K`, and the model solar-rotation convention `2.86533e-6 rad/s`.

The 1-D and 3-D modules still contain separate hard-coded decimal conversion
factors for km/s, cm^-3, nT, and km^-1 even though physical constants are now
centralized. CFG02 exercises and compares both paths. Consolidating these
factors is recommended technical debt, but this test does not automatically
change them because doing so could conceal the behavior it is meant to audit.

## PAR01-PAR03: corrected 3-D Parker magnetic field

The production 3-D Parker field now uses an explicit solar-rotation axis and a
local spherical basis at every evaluation point.  The azimuthal direction is

```text
e_phi = (Omega_hat x e_r) / |Omega_hat x e_r|
```

and the local pitch is proportional to
`|Omega_hat x e_r| = sin(theta_local)`.  This corrects the previous
implementation, which used `(Omega_hat x e_r) x e_r` (a meridional direction)
and one global `sin_theta` for the entire 3-D domain.

`Params::solar_rotation_axis` specifies the global solar-rotation axis and is
normalized when a `StepState` is prepared.  A zero or non-finite axis is
rejected because it cannot define a Parker azimuthal direction.  The existing
`Params::sin_theta` field is retained only for backward-compatible
normalization of `B1AU_nT`: it specifies the reference `sin(colatitude)` at
which `B1AU_nT` is interpreted as total field magnitude at 1 AU.  It no longer
sets the local 3-D Parker winding.  The step cache stores `k_AU = Omega*AU/Vsw`;
each point multiplies this by its own `r_AU*sin(theta_local)`.

The pole is handled analytically.  When the radial direction is parallel to
the rotation axis, `sin(theta_local)=0`, so `B_phi=0` and the field is purely
radial.  The implementation does not manufacture an arbitrary azimuthal unit
vector at the coordinate singularity.

The following deterministic 3-D tests were added:

- `PAR01` — equatorial Parker-vector orientation.  With the rotation axis along
  +Z and the point at +X, it requires the outward radial component to lie along
  +X, the Parker azimuthal component to lie along -Y for the current polarity
  convention, and the meridional Z component to vanish.  It also verifies the
  requested total-field normalization at the reference latitude.
- `PAR02` — arbitrary latitude and arbitrary rotation axis.  It compares the
  production vector with an independent analytical Parker reference away from
  the equator and repeats the check after rotating the solar axis, preventing a
  hard-coded +Z implementation from passing accidentally.
- `PAR03` — polar-limit regularity.  It verifies the exact north and south
  rotation poles and a near-pole point, requiring the transverse field to tend
  continuously to zero without NaN, Inf, or an arbitrary transverse direction.

The analytical reference in `3d/test_parker.cpp` is intentionally independent
of the production Parker helper.  Component comparisons use roundoff-level
absolute tolerances scaled to the expected Tesla magnitude.  These tests do not
cover numerical solenoidality or field-line tangent/path-length validation;
those remain separate planned tests (`PAR04` and later) so the individual
validation requirements remain diagnostically focused.

## GEO01-GEO08: corrected finite shock geometry, normals, and flank speed

The 3-D shock geometry has been reworked so that the production model no longer
uses the former `ConeSSE` approximation

```text
R(theta) = R_apex cos(theta)^m
```

with a radial normal and an artificial clamped radius beyond the configured
half width.  `ShockShape::SSE` now denotes a true finite self-similar-expansion
spherical cap.  `ShockShape::ConeSSE` is retained only as a source-compatible
enum alias and has the same corrected SSE semantics; `flank_slowdown_m` remains
in `Params` only for source/input compatibility and is ignored by the SSE
geometry.

For an apex distance `R_apex` and angular half width `lambda`, the generating
sphere has center distance and radius

```text
c = R_apex / (1 + sin(lambda))
a = c sin(lambda).
```

A heliocentric ray separated from the CME axis by `alpha` intersects the
outward front at

```text
R(alpha) = c cos(alpha) + sqrt(a^2 - c^2 sin(alpha)^2),
```

provided `alpha <= lambda`.  At `alpha=lambda` the ray is tangent to the
generating sphere.  Directions beyond the half width return `exists=false`;
the public geometry API also returns zero radius/normal sentinels so a caller
cannot mistake an absent surface for a physical flank.

The exact outward normal is the normalized level-set gradient

```text
n_hat = (R e_r - c e_CME) / a.
```

All supported shapes are treated as self-similar.  For a fixed surface
direction,

```text
dR/dt = V_apex R/R_apex,
V_sh,n = V_apex (R/R_apex) (e_r dot n_hat).
```

This replaces the former independent `cos(theta)^m` flank-speed factor and
also corrects the ellipsoid, whose flanks previously inherited the full apex
speed.  For a sphere the formula reduces exactly to `V_sh,n=V_apex`; at the
mathematical SSE tangent boundary the normal projection tends to zero.

`Model::shape_radius_normal()` and `Model::diagnose_direction()` now return a
boolean surface-existence flag.  Existing callers that ignore the return value
remain source-compatible, but finite-width-aware code should always test it.
The Cartesian plasma/field evaluators do so internally: outside the SSE angular
support they return the undisturbed ambient solar wind/Parker field and do not
construct a sheath, ejecta, shock normal, or magnetic amplification from a
fabricated flank.

The following deterministic tests validate the corrected geometry:

- `GEO01` — Sun-centered spherical reference: radius is direction independent
  and the outward normal equals the radial unit vector.
- `GEO02` — true SSE apex: the directional radius equals the configured apex
  distance and the normal equals the CME propagation direction.
- `GEO03` — SSE tangent flank: the boundary point at `alpha=lambda` remains a
  finite valid surface point; the ray is tangent and the normal radial
  projection tends to zero.
- `GEO04` — finite-width enforcement: any direction beyond the half width
  returns `exists=false`, zero geometry sentinels, `rc=1`, and `V_sh,n=0` from
  the diagnostic API.
- `GEO05` — SSE level-set residual: returned surface points satisfy
  `|x-c e_CME|=a` to near machine precision over multiple angles/azimuths.
- `GEO06` — analytical normal validation: SSE normals are compared with an
  independent finite-difference gradient of the dimensionless spherical level
  set; the ellipsoid normal is independently checked against its level-set
  gradient.
- `GEO07` — normal-speed validation: the reported `V_sh,n` for sphere, SSE, and
  ellipsoid is compared with the normal projection of centered finite-
  difference surface motion at neighboring times.
- `GEO08` — rotational covariance: rotating the CME axis, solar axis, and query
  direction together must rotate the normal while leaving radius, normal
  speed, and scalar physical compression unchanged.

The analytical/reference calculations in `3d/test_geometry.cpp` intentionally
do not call the production geometry helper.  They are independent checks of
surface radius, level-set membership, normal orientation, and motion.  Test
tolerances are roundoff- or finite-difference-level tolerances appropriate to
each calculation and must not be loosened simply to make a production result
pass.

The geometry tests remain diagnostically focused on shape, normals, finite
angular support, and normal speed.  Shock existence/compression/downstream
physics is now covered independently by `SHK01`-`SHK12` below.  Shock-surface
mesh apex/seam topology remains a separate planned work package.



## KIN01-KIN08: shared CME/shock-apex kinematics

`swcme_kinematics.hpp` is the production source of apex radius and speed for
both `swcme1d` and `swcme3d`.  The kinematics tests are classified `COMMON` and
exercise both the common SI solver and the two public dimensional wrappers.
The independent DBM reference in `core/test_kinematics.cpp` evaluates the
closed-form equations directly and does not call the production DBM helper.

The production DBM solves

```text
DeltaV0 = V0 - Vsw
DeltaV(t) = DeltaV0 / (1 + Gamma |DeltaV0| t)
R(t) = R0 + Vsw t
       + sign(DeltaV0) log(1 + Gamma |DeltaV0| t) / Gamma.
```

This sign-aware form is required for both fast and slow CMEs.  `Gamma=0` uses
the exact ballistic solution.  For very small `x=Gamma*|DeltaV0|*t`, the
production code evaluates the logarithmic distance through a short
`log1p(x)/x` series; `KIN04` explicitly straddles that numerical branch
boundary to guard against a discontinuity.

The data-driven mode uses monotone PCHIP radius interpolation.  The input time
table must be strictly increasing and radius must be nondecreasing.  The
interpolant passes through each knot exactly and its derivative is returned as
the apex speed.  Out-of-range queries return `OUTSIDE_TIME` by default;
ballistic endpoint continuation is available only when requested explicitly.

The tests are:

- `KIN01` — fast-CME DBM versus the independent closed form at multiple times,
  including exact 1-D/3-D wrapper equivalence;
- `KIN02` — slow-CME sign-aware branch; verifies acceleration toward `Vsw`
  without clipping or overshoot and checks both wrappers;
- `KIN03` — exact `Gamma=0` ballistic radius/speed in the common, 1-D, and 3-D
  paths;
- `KIN04` — small-`Gamma` analytical accuracy and continuity across the
  series/direct-log numerical branch;
- `KIN05` — long-time stability, monotonic radius, finite state, and monotonic
  approach of speed toward the ambient wind for fast and slow CMEs;
- `KIN06` — data-driven PCHIP exactness at every height-time knot, derivative
  consistency, and 1-D/3-D use of the same trajectory;
- `KIN07` — dense-grid monotonicity, nonnegative propagation speed, and absence
  of cubic overshoot between data knots; and
- `KIN08` — explicit out-of-time status, optional ballistic continuation, and
  rejection of duplicate-time or decreasing-radius tables.

The DBM analytical comparisons use `1e-12` relative accuracy where specified
by the validation plan.  PCHIP knot radii are required to be exact to
`1e-13` relative.  Tolerances must not be relaxed merely to obtain PASS.

The 1-D and 3-D public `Params` structures expose the same kinematic mode plus
`data_time_s`, `data_radius_Rs`, and `data_extrapolation`.  The 1-D convenience
method `SetDataDrivenKinematics()` sets these fields and switches the mode to
`DataDriven`.  The wrappers convert radii from solar radii to SI only when
constructing the common configuration.

A data-driven query outside the measurement interval causes `prepare_step()`
to throw a clear runtime error unless a continuation policy was selected.  This
is intentional: silent cubic extrapolation is a modeling assumption and is
not allowed to masquerade as measured/constrained kinematics.

## SHK01-SHK12: fast-shock existence and ideal-MHD Rankine-Hugoniot validation

The `SHK` group validates the shared production shock solver in
`swcme_shock.hpp` and its integration into both the 1-D and 3-D models.  The
classification is `COMMON`: the physical jump calculation is shared even when
geometry-specific inputs such as the 3-D surface normal are different.

The solver first evaluates the upstream fast-mode speed along the shock normal.
A geometric front is not automatically a shock.  The physical condition is

```text
U1n = Vsh,n - V1.n > c_fast,
M_fast = U1n/c_fast > 1.
```

If this condition is not met, the result is `has_shock=false`, compression is
exactly one, and downstream equals upstream.  No sheath compression floor is
allowed to override this decision.

For a physical fast shock, the solver works in a frame moving with the normal
shock speed.  For a trial compression `r=rho2/rho1`, mass conservation fixes
`u2n=u1n/r`.  Tangential electric-field and tangential momentum continuity form
a 2x2 linear system for `B2t` and `u2t`; normal momentum determines `p2`.  A
bracketed solve of total-energy-flux conservation selects the compressive
fast-shock branch.  The trivial `r=1` solution is divided out of the scalar
residual so it cannot be mistaken for the physical shock root.

The individual tests are:

- `SHK01` — fast-shock existence and no-shock threshold.  Includes permanent
  1-D and 3-D regression checks proving that `sheath_comp_floor` cannot create
  a shock when `Vsh=Vsw`.
- `SHK02` — acute `theta_Bn` calculation and invariance under `B -> -B`.
- `SHK03` — strictly parallel limit against the independent gas-dynamic normal
  shock solution.
- `SHK04` — strictly perpendicular ideal-MHD benchmark against an independent
  test-side scalar reduction.
- `SHK05` — oblique benchmark grid and continuity of the selected compression
  branch as shock speed/obliquity vary.
- `SHK06` — mass-flux conservation.  It also samples the production 3-D field
  immediately downstream and verifies that density, velocity, and magnetic
  field approach the exact RH downstream state instead of the ambient wind.
- `SHK07` — continuity of the normal magnetic-field component.
- `SHK08` — tangential ideal-MHD electric-field continuity.
- `SHK09` — vector momentum-flux conservation.
- `SHK10` — total ideal-MHD energy-flux conservation.
- `SHK11` — physical admissibility: positive downstream pressure/density,
  compressive branch, strong-shock bound for gamma=5/3, and entropy increase.
- `SHK12` — near-Mach-one conditioning and smooth approach of compression to
  unity without an empirical floor.

The independent parallel/perpendicular references do not call the production
nonlinear shock solver.  Conservation tests recompute the conserved fluxes from
the returned primitive states.  This prevents a test from passing merely by
reusing the same internal algebra that produced the result.

Typical direct use is:

```sh
./output/test_swcme --test SHK01
./output/test_swcme --test SHK04
./output/test_swcme --test SHK10
./output/test_swcme --all
```

A shock test failure must not be repaired by increasing the compression floor
or loosening conservation tolerances.  Diagnose the frame transformation,
normal direction, upstream state, root bracket, and downstream reconstruction
first.  A super-fast state for which the nonlinear solver cannot identify an
admissible root is reported with `solver_converged=false` and must not be used
for SEP source physics.


## CON01-CON08: observer-shock magnetic connectivity and cobpoint tracking

`CON01`-`CON08` validate the production 3-D connectivity API implemented by
`swcme3d::Model::observer_connectivity()`.  The connectivity solver does not
maintain a second copy of the shock model: candidate points are tested against
`shape_radius_normal()` and final cobpoints obtain their local plasma/shock
state from `shock_state_direction()`.

The observer Parker line is analytical.  For the same Parker field used by the
production field evaluator,

```text
Delta phi(r) = -Omega (r-r_obs) / V_sw.
```

The observer radial direction is rotated around the configured solar axis by
this amount; colatitude remains constant.  The production parameter
`solar_rotation_rate_rad_s` is shared by field evaluation and connectivity, so
there cannot be a hidden Omega mismatch between the two modules.  The default
value remains the SWCME solar-rotation convention; zero rotation is allowed to
exercise the exact radial reference case.

The intersection residual is

```text
h(r) = r - R_shock[u_Parker(r)].
```

The solver scans from a configurable inner radius to the observer.  It refines
sign-changing roots with bisection, searches local minima of `abs(h)` so tangent
roots that merely touch zero are not missed, and explicitly refines transitions
of the finite-SSE surface-existence flag.  All accepted roots must satisfy the
configured surface-residual tolerance.  Roots are retained in increasing
radius; the default selected cobpoint is the outermost root, corresponding to
the first surface encountered while tracing inward from the observer.

Each `ConnectivityRoot` stores:

- cobpoint radius and Cartesian position;
- signed surface residual;
- analytical Parker arc length from cobpoint to observer; and
- the complete production `LocalShockState`, including `has_shock`, normal,
  normal speed, `theta_Bn`, fast Mach number, compression, and full upstream /
  downstream MHD states.

A geometrical connection and a physical fast shock remain separate concepts.
The connectivity state can therefore be geometrically connected while the
embedded `LocalShockState::has_shock` is false; an SEP source must check the
latter before injection.

The tests are:

- `CON01` — zero-solar-rotation radial limit.  A radial observer line is
  intersected with sphere and SSE geometries and compared with exact analytic
  radius/position/path-length references.
- `CON02` — nonzero-rotation Parker line intersecting a Sun-centered sphere.
  The sphere fixes the root radius exactly while an independent Rodrigues
  rotation verifies the longitude/sign of the production Parker mapping.
- `CON03` — no connection to a finite-width SSE shock.  A field line wholly
  outside the configured cap must return `connected=false` and no fabricated
  flank root.
- `CON04` — exact and near-tangent SSE connection.  The exact half-width ray is
  retained as a valid tangent root, a slightly interior ray connects, and a
  slightly exterior ray remains disconnected.
- `CON05` — multiple intersections and deterministic root selection.  A tightly
  wound Parker line through a strongly non-spherical ellipsoid generates
  multiple physical roots; all are returned in radial order and the outermost
  selection is stable under scan refinement.
- `CON06` — time-continuous cobpoint history.  A ballistic finite SSE front
  evolves from disconnected to connected; the history must show one physical
  onset and continuous/monotone cobpoint motion thereafter without one-step
  classification flicker.
- `CON07` — cobpoint-to-`ShockState` consistency.  Shock quantities stored in a
  connectivity root must be identical to a direct production
  `shock_state_direction()` query at that cobpoint direction.
- `CON08` — analytical Parker path length.  The path length carried by a
  cobpoint is compared with an independent closed-form arc-length reference
  and, away from the pole, must exceed simple radial separation.

The default connectivity search tolerance is much tighter than the validation
campaign's `1e-8 AU` root-position target.  Tests deliberately repeat selected
cases with different radial scan densities so correctness cannot depend on a
particular subdivision.  Tolerances must not be relaxed merely to mask missed
roots or unstable tangent classification.

Typical direct use is:

```sh
./output/test_swcme --test CON01
./output/test_swcme --test CON04
./output/test_swcme --test CON06
./output/test_swcme --test CON08
./output/test_swcme --all
```
