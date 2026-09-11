# Centralized configuration validation and unit handling

This update completes the configuration/unit remediation identified by CFG01 and
CFG02.  The physical models still expose their established heliophysics-facing
`Params` fields, but all repeated conversion factors and common physical-range
checks now pass through shared infrastructure.

## New common components

### `swcme_units.hpp`

This header is the single production source for conversions between public
heliophysics units and internal SI quantities:

- km/s <-> m/s;
- cm^-3 <-> m^-3;
- nT <-> T;
- km^-1 <-> m^-1;
- AU <-> m;
- solar radii <-> m;
- hours <-> seconds; and
- degrees <-> radians.

A conversion function never changes the physical value for admissibility.  In
particular, `0 km/s` converts to exactly `0 m/s`.  Whether zero is a valid
solar-wind speed is a separate configuration question.  This separation fixes
the previous 1-D behavior that silently replaced `0 km/s` with `1 m/s`.

### `swcme_config.hpp`

This header defines `ValidationResult`, structured validation issues, and the
common rules applied to both dimensional interfaces.  Each invalid issue records
its field name, a stable error code, and the violated requirement.  A caller can
inspect `Model::validate()` without changing model state; `prepare_step()` calls
the same validator and rejects invalid input before entering any physics.

Common validation covers:

- finite, strictly positive ambient solar-wind speed;
- finite, strictly positive reference density and temperature;
- non-negative magnetic-field magnitude and drag coefficient;
- adiabatic index greater than one;
- Parker reference `sin(theta)` in `[0,1]`;
- positive kinematic reference radius and non-negative launch speed;
- non-negative region/smoothing thicknesses;
- physically allowed sheath/profile factors; and
- complete DATA_DRIVEN table validity (matching sizes, finite positive radii,
  strictly increasing times, nondecreasing radii).

The 3-D wrapper adds geometry rules for non-zero finite CME/solar axes,
non-negative rotation rate, positive ellipsoid axis ratios, and finite SSE half
width in `(0,pi/2]`.

## Model behavior

Both 1-D and 3-D `prepare_step()` now perform the same validation before unit
conversion.  Invalid configurations throw `std::invalid_argument` containing a
summary of all detected fields.  Unit conversion is then deterministic and
contains no physical clipping.  The common kinematics component retains its own
defensive SI validation as a second internal guard.

The change intentionally does not remove every legacy `finite_or()` in the 3-D
field/region evaluators; that broader numerical-status cleanup remains a
separate remediation task.  It does remove configuration-derived fallbacks in
`prepare_step()` where valid parameters have already been established.

## Validation

`CFG01` is now the intended configuration-rejection test rather than a constants
placeholder.  It exercises the public side-effect-free validation API and
confirms that `prepare_step()` rejects representative invalid configurations,
including zero wind speed, negative density/Gamma, invalid Parker latitude,
NaN values, malformed data-driven tables, zero 3-D axes, invalid SSE width,
invalid ellipsoid ratio, and negative solar rotation rate.  It also verifies
that multiple invalid fields are reported in one pass.

`CFG02` now tests `swcme_units.hpp` directly for forward and inverse conversions,
then verifies that valid 1-D and 3-D prepared states reach the same SI values.
It includes an independent Alfven-speed dimensional smoke test.

At the time of this update the complete deterministic suite contains 43 tests
and all 43 pass.
