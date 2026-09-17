# Observer-to-shock connectivity / cobpoint implementation notes

## Purpose

This update adds the missing SWCME capability needed to distinguish evolving
magnetic connection from perpendicular SEP transport.  A stationary observer's
nominal Parker field line is traced analytically and intersected with the
production time-dependent shock surface.  Every accepted intersection carries
the local production shock state and field-line distance to the observer.

## Physics and geometry

The production Parker field uses

`B_phi/B_r = -Omega r sin(theta)/V_sw`.

Because a magnetic-field-line tangent satisfies
`r sin(theta) dphi/dr = B_phi/B_r`, the Parker line has
`dphi/dr=-Omega/V_sw` at fixed colatitude.  SWCME therefore constructs the
observer-anchored line exactly by rotating its radial direction around the
solar-rotation axis by `-Omega (r-r_obs)/V_sw`.

The line/shock residual is `h(r)=r-R_shock[u(r)]`.  Finite SSE geometry can be
absent along part of the line; such samples are treated as no-surface states,
not as a clamped/fabricated radius.

## Numerical root discovery

The implementation combines four mechanisms:

1. adaptive radial scanning (with additional samples if Parker winding is
   strong);
2. bisection for sign-changing roots;
3. golden-section minimization of `|h|` around local minima to detect tangent
   roots without a sign change; and
4. bisection of finite-surface validity transitions so first/last SSE
   connection at the angular boundary is not missed.

Candidate roots are deduplicated, reevaluated through the production shock API,
and accepted only if their final surface residual satisfies the configured
metric tolerance.  All roots are preserved in increasing radius.  The default
selected root is the outermost one (first encountered from the observer).

## Path length

For `k=Omega sin(theta)/V_sw`, the exact Parker arc-length integral is

`F(r)=0.5 [ r sqrt(1+(k r)^2) + asinh(k r)/k ]`.

The cobpoint path length is `|F(r_obs)-F(r_cob)|`.  The `k=0` limit is handled
analytically as radial distance.

## API additions

`swcme3d.hpp` now exposes:

- `ConnectivityStatus`;
- `ConnectivityOptions`;
- `ConnectivityRoot`;
- `ConnectivityState`;
- `ConnectivityHistorySample`;
- `Model::parker_field_line_point()`;
- `Model::parker_field_line_length()`;
- `Model::observer_connectivity()`; and
- `Model::observer_connectivity_history()`.

`Params::solar_rotation_rate_rad_s` and the corresponding `StepState` cache are
also explicit so field evaluation and connectivity always use the same Omega.
The default is unchanged; zero rotation is supported as an exact validation
limit.

## Validation

`CON01`-`CON08` cover radial, spiral, disconnected, tangent, multiple-root,
time-history, ShockState-consistency, and analytical path-length cases.  The
full deterministic suite should still have only the pre-existing CFG02 failure
associated with the legacy 1-D zero-velocity unit/configuration clamp.

## Scope / remaining work

This connectivity model intentionally follows the nominal upstream analytical
Parker field.  It does not attempt to trace through a disturbed ICME magnetic
field after shock passage, and it does not yet implement perpendicular particle
transport.  Those limitations are deliberate for the controlled SEP study:
connectivity is computed from the corrected upstream Parker geometry, while
perpendicular diffusion remains a separate transport parameter in AMPS.
