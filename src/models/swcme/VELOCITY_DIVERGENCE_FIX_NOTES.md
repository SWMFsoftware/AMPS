# Velocity-divergence correction and validation

## Problem corrected

Before this change the 3-D model evaluated velocity divergence everywhere with

```text
div(V) = (1/r^2) d(r^2 V_r)/dr
```

by sampling only along a heliocentric ray.  That relation is exact for a purely
radial vector field.  It is not the divergence of the general FULL_ICME field:
an oblique Rankine-Hugoniot jump can contain tangential velocity, and finite
SSE/ellipsoid fronts make the region state vary with angular position.  The
radial-only operator therefore omitted real contributions to the transport
compression term.

The 1-D implementation also used a finite-difference approximation even though
its radial geometry and smooth profiles permit an analytical derivative.  That
introduced avoidable step-size noise into the SEP adiabatic energy-change term.

## Production design

`swcme_divergence.hpp` now contains the dimension-independent divergence
mathematics.

For a radial velocity field it exposes the exact decomposition

```text
div(V) = 2 V_r/r + dV_r/dr
       = geometric term + profile-derivative term.
```

The separate terms are retained in `RadialTerms` so tests and diagnostics can
identify whether an error comes from spherical geometry or from the profile
derivative.

For a general vector field the same header provides a second-order Cartesian
operator

```text
div(V) = dVx/dx + dVy/dy + dVz/dz.
```

The operator uses centered differences whenever all stencil points are inside
the model domain.  If exactly one side crosses the explicit inner model
boundary, it switches to the corresponding second-order one-sided derivative.
It never clips a sample radius into the model domain.

## 1-D treatment

The 1-D velocity is purely radial by construction.  In SHOCK_ONLY,
`dV_r/dr=0`, giving the exact baseline result

```text
div(V) = 2 V_sw/r.
```

For FULL_ICME, `swcme_regions.hpp` now returns `RadialVelocityState`, containing
both `V_r` and its analytical radial derivative.  The derivative is taken from
the exact same smoothstep functions used to construct the production velocity:

- resolved shock transition;
- phenomenological sheath relaxation;
- sheath/ejecta leading transition;
- ejecta plateau;
- ejecta/ambient trailing transition.

This eliminates the possibility that the divergence differentiates a profile
that differs from the velocity actually returned to the transport solver.
The legacy `dr_frac` argument remains in the 1-D API only for source
compatibility and is ignored.

## 3-D treatment

`Model::compute_divV_checked()` is the canonical public entry point.

- SHOCK_ONLY: returns the exact `2 V_sw/r` result.  The output is independent of
  finite-difference step and Cartesian direction to roundoff.
- FULL_ICME: calls `compute_divV_cartesian_checked()` and evaluates the complete
  Cartesian Jacobian trace from the production velocity field.

`compute_divV_cartesian_checked()` is also public so a caller can explicitly run
convergence studies in a case where the canonical result is analytical.

The historical `compute_divV_radial_checked()` and `compute_divV_radial()`
symbols remain for source compatibility.  They now delegate to the canonical
mode-aware implementation; despite their old names, they no longer apply the
radial approximation to FULL_ICME.

## Numerical step

The Cartesian operator uses

```text
h = max(1000 m, dr_frac * r).
```

`dr_frac` is validated as finite and positive in the numerical 3-D branch.
The 1000 m floor prevents a vanishing stencil very near the inner model
boundary while being negligible at heliospheric scales.  The convergence test
uses explicit step refinement and does not assume the default step is optimal
for every future science configuration.

## Validation

Three deterministic gates were added.

### DIV01 — analytical constant radial wind

SHOCK_ONLY 1-D and 3-D results are evaluated at several radii and at several
unrelated Cartesian directions.  Every result must equal

```text
2 V_sw/r
```

to roundoff.  Any direction dependence indicates an incorrect 3-D baseline
operator.

### DIV02 — manufactured radial-flow divergence

A positive polynomial profile

```text
V_r = V0 (1 + a x + b x^2),  x=r/r0
```

has an independent symbolic divergence.  The production radial helper must
match it to better than `1e-10` relative error.  The test also compares the
analytical derivative stored by representative FULL_ICME region profiles with
an independent centered numerical derivative of the public 1-D velocity.

### DIV03 — general Cartesian divergence convergence

A cubic manufactured vector field with known analytic divergence is evaluated
with four successively halved Cartesian steps.  The observed truncation error
must converge at second order.  A second part evaluates the production 3-D
Cartesian operator on an off-axis constant radial solar wind and verifies that
it converges toward the exact `2 V_sw/r` value.

These tests distinguish analytical-physics correctness from numerical-operator
convergence and prevent either path from being made green by substituting the
other.
