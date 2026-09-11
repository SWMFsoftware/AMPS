# SWCME explicit numerical-status and error-propagation update

## Purpose

This update removes the remaining policy in which a failed numerical operation
could be converted into a physically plausible value merely so evaluation could
continue.  Examples of the old behavior included clipping a radius to the lower
solar-wind boundary, replacing non-finite field values by zero, replacing a
failed vector normalization by `+X`, and allowing output code to hide non-finite
physics values.  Such behavior is especially dangerous in SEP calculations
because a fabricated ambient state can look numerically smooth while changing
shock strength, connectivity, or `div(V)`.

The governing rule is now:

> A numerical/physics failure is data.  It must reach the caller as an explicit
> status (or, for a source-compatible legacy wrapper, as an exception).  It may
> not be converted into a different physical state.

Expected physical/geometrical outcomes remain distinct from failures.  In
particular, `NO_SURFACE` for a direction outside a finite SSE cap and
`NO_SHOCK` for a sub-fast front are not numerical errors.

## Shared `ModelStatus`

`swcme_status.hpp` introduces the common status object used by checked hot-path
APIs:

```cpp
swcme::ModelStatus
```

The status carries a `StatusCode`, the first failed batch-sample index when
applicable, a short context string, and optionally the offending scalar value.
The current codes are:

```text
OK
NO_SURFACE
NULL_POINTER
NONFINITE_INPUT
OUTSIDE_MODEL_DOMAIN
DEGENERATE_VECTOR
INVALID_NUMERICAL_STEP
GEOMETRY_FAILURE
SHOCK_SOLVER_FAILURE
NONFINITE_RESULT
INVALID_MESH
FILE_OPEN_FAILURE
FILE_WRITE_FAILURE   (reserved for writers that distinguish write failures)
```

`OK` is success.  `NO_SURFACE` is an expected finite-geometry result and can be
handled without treating it as an error.  All other non-OK codes are failures.

`ModelStatus::summary()` produces a deterministic diagnostic suitable for log
messages.  `swcme::throw_if_error()` is used by the old void-returning APIs so
existing source continues to compile while failures stop being ignored.

## Checked evaluator APIs

The following status-returning paths are now available:

### 1-D

```cpp
Model::evaluate_radii_fast_checked(...)
Model::evaluate_radii_with_B_div_checked(...)
```

### 3-D

```cpp
Model::shock_state_direction_checked(...)
Model::evaluate_cartesian_fast_checked(...)
Model::evaluate_cartesian_with_B_checked(...)
Model::evaluate_cartesian_with_B_div_checked(...)
Model::compute_divV_radial_checked(...)
```

The existing public wrappers retain their original signatures.  They call the
checked implementation and throw `std::runtime_error` if the returned status is
a failure.  This preserves source compatibility while eliminating the previous
silent-repair semantics.  AMPS-facing code should prefer the checked forms so a
failure can be propagated through its own error/status machinery without
exceptions.

Batch evaluators return the **first failed sample index**.  Values at the failed
sample are not overwritten with zero, ambient flow, or any other substitute.
Earlier successfully evaluated samples may already have been written; callers
that require transactional batches should evaluate into temporary storage.

## Explicit radial model domain

The analytical Leblanc/Parker background has a documented lower domain boundary
of

```text
r_min = 1.05 R_sun.
```

`swcme_solarwind.hpp` no longer implements this as
`max(r, r_min)`.  Low-level density/Parker equations evaluate the radius they
are given, while public 1-D/3-D science evaluators first verify the domain and
return `OUTSIDE_MODEL_DOMAIN` for `r < r_min`.

This distinction is important: a request at `0.9 R_sun` must not silently return
the plasma state at `1.05 R_sun`.

Prepared shock surfaces are checked against the same domain.  If kinematics
place the front below the supported analytical background, `prepare_step()`
raises an explicit error instead of sampling a fabricated boundary state.

## Vector normalization and floating-point robustness

The former 3-D normalization helper could replace a zero vector with `(1,0,0)`.
That meant a bad CME/shock direction could become a legitimate +X calculation.
The replacement uses `std::hypot` with an explicit magnitude threshold and
returns failure for a zero or non-finite vector.

The same principle is applied to the shared RH solver.  Every upstream primitive
component, shock normal, shock speed, and `gamma` is validated before derived
wave speeds are evaluated.  Vector norms use `std::hypot`; a NaN magnetic-field
component can no longer collapse to a zero magnetic norm and masquerade as an
unmagnetized physical state.

Roundoff clamps are retained only where a quantity is mathematically bounded.
Examples include an `acos` argument or a tiny negative SSE discriminant caused
by floating-point roundoff.  Physical Mach numbers, compression ratios, and
failed solver states are not clamped into acceptance.

## Rankine-Hugoniot solver outcome

`swcme_shock.hpp` now exposes `SolveStatus` in every `JumpResult`:

```text
NO_SHOCK
SOLVED
INVALID_INPUT
NO_PHYSICAL_BRACKET
INVALID_ACCEPTED_STATE
CONSERVATION_FAILURE
```

This removes the need to infer failure from `compression == 1` or from the old
`solver_converged` compatibility flag.  A sub-fast state is explicitly
`NO_SHOCK`; malformed/non-finite input is `INVALID_INPUT`; a super-fast state
whose nonlinear RH solve fails is a solver failure and is not converted to
ambient downstream plasma.

The 1-D prepared state and 3-D surface shock query propagate a failed physical
RH solve as `SHOCK_SOLVER_FAILURE` rather than silently selecting
`compression=1`.

## Field and region evaluation

Removed behaviors include:

- radius clipping to `1.05 R_sun` in science evaluators;
- `n=0` or `V=0` replacement after a non-finite regional calculation;
- `Br/Bphi/Bmag=0` replacement after magnetic-field failure;
- `div(V)=0` replacement after a failed finite difference;
- arbitrary +X replacement for a degenerate direction;
- `finite_or()` repair in 3-D mesh/output data.

`SHOCK_ONLY` and `FULL_ICME` physics are otherwise unchanged by this update.
Expected absence of a finite shock surface is still handled as `NO_SURFACE`, so
a point outside an SSE angular cap correctly receives ambient background rather
than being classified as a numerical error.

## Divergence boundary policy

The current `div(V)` implementation still uses the established radial finite-difference formulation; replacing it with the analytical upstream expression
is a separate planned remediation.

What changes here is failure handling.  A centered stencil that would cross the
lower Parker/Leblanc domain uses an explicit one-sided lower-domain stencil with
its inner sample exactly at `r_min`.  This is a documented numerical boundary
condition, not a physical-value fallback.  Invalid `dr_frac`, non-finite
coordinates, degenerate denominators, or failed nested field evaluations return
an explicit status.

## Mesh and Tecplot behavior

Mesh construction and triangle-metric calculation no longer pass non-finite
coordinates, normals, compression, normal speed, or derived triangle quantities
through `finite_or()`.  Invalid mesh/connectivity data raise an explicit error.

Checked 3-D Tecplot writer entry points validate mesh/metric data before writing
and reject non-finite physics.  A file-open failure is reported separately from
physics-data failure.  The `FILE_WRITE_FAILURE` status is reserved for a later
writer-I/O refactor that checks every `fprintf`/`fclose` result; this update does
not claim transactional output or detailed mid-file I/O diagnostics.

At the time of Fix 11, the surface-mesh **topology** was intentionally left for
a separate remediation.  Fix 12 now removes the apex/seam duplicate vertices
and zero-area cells at construction time.  The status rule established here
remains unchanged: malformed or degenerate connectivity is reported as
`INVALID_MESH` rather than filtered or repaired after the fact.

## Validation tests

Five deterministic regression tests were added:

- **ERR01 — 1-D outside-domain propagation.** A query below `1.05 R_sun`
  returns `OUTSIDE_MODEL_DOMAIN`, records the failed sample, leaves the failed
  output untouched, and makes the legacy wrapper throw instead of clipping.
- **ERR02 — 3-D non-finite Cartesian input.** A NaN coordinate returns
  `NONFINITE_INPUT` at the correct array index and is not replaced by zero
  output.
- **ERR03 — explicit RH solver outcomes.** Degenerate/non-finite shock inputs,
  a valid sub-fast state, and a valid fast shock must report different
  `SolveStatus` values (`INVALID_INPUT`, `NO_SHOCK`, and `SOLVED`).
- **ERR04 — degenerate direction rejection.** A zero direction returns
  `DEGENERATE_VECTOR`; no +X shock surface is manufactured.
- **ERR05 — output rejection of non-finite physics.** A NaN shock-mesh value is
  rejected before the checked surface writer creates a file.

The complete deterministic suite after this update contains **63 tests**.

## Compatibility and caller guidance

Code that supplied valid inputs sees the same physical results to normal
floating-point tolerances.  The intentional behavior change is for invalid or
unsupported states: a legacy void evaluator can now throw where it previously
returned plausible-looking repaired data.

Recommended AMPS integration pattern:

```cpp
const swcme::ModelStatus s = model.evaluate_cartesian_with_B_checked(...);
if (!s.ok()) {
    // propagate s.code, s.sample_index, and s.summary() to the AMPS driver/log
    // rather than continuing with a fabricated background state.
}
```

Do not catch an SWCME error merely to replace the result with ambient plasma.
Recovery, if scientifically justified, should occur at the application layer
and should remain visible in run metadata.
