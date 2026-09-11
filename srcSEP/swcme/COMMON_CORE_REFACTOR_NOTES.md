# SWCME common-core refactor notes

## Purpose

This update removes the remaining duplicated **ambient solar-wind preparation**
from the 1-D and 3-D SWCME models.  Constants, unit conversion, configuration
validation, apex kinematics, and the ideal-MHD shock solver were already shared
by earlier updates; however, the two dimensional wrappers still independently
normalized the Leblanc density and Parker field and independently assembled the
common SI/kinematic state.  That duplication was a long-term regression risk.

## New production components

### `swcme_solarwind.hpp`

This is the authoritative implementation of dimensionality-independent ambient
physics:

- immutable Leblanc coefficients `A`, `B`, and `C`;
- normalization of the Leblanc profile to the configured 1-AU density;
- fast SI density evaluation `C2/r^2 + C4/r^4 + C6/r^6`;
- Parker radial-field normalization from the configured total field at a
  documented reference latitude;
- Parker `Br` and `Bphi` components for any local `sin(theta)`;
- Cartesian Parker-vector construction for an arbitrary normalized solar axis;
- the current proton-only pressure closure `p=n k_B T` used by the shock model.

The 1-D wrapper supplies its fixed latitude to `parker_components()`.  The 3-D
wrapper supplies its solar axis and local radial direction to
`parker_field_cartesian()`.  Thus the scalar Parker physics is identical even
though the geometrical representation differs.

### `swcme_core.hpp`

`swcme::core::prepare()` accepts one common configuration in public units and:

1. converts ambient inputs to SI;
2. prepares the shared `swcme::solarwind::PreparedState`;
3. converts the kinematic inputs to SI;
4. evaluates the shared BALLISTIC/DBM/DATA_DRIVEN kinematics.

Both dimensional models store this result in `StepState::common`.

## Backward compatibility

Existing `StepState` members such as `C2/C4/C6`, `Br1AU_T`, `k_AU`, ambient
speed, apex radius, and apex speed remain populated.  They are now **mirrors**
of the canonical common state instead of independent calculations.  This avoids
breaking existing callers while establishing one authoritative production path.

The 1-D `StepState` also stores the complete shared `swcme::shock::JumpResult`.
This provides an auditable shock-state bridge for dimensional-equivalence tests
and for the future common SEP-source interface.

## Deliberately not unified in this update

The phenomenological sheath/ejecta region shaping remains different between the
1-D and 3-D wrappers.  It includes known behavior scheduled for the separate
region-model remediation work package, so moving it into the common core now
would merely centralize behavior that has not yet been corrected and validated.

Likewise, full geometry is intentionally not common: the 1-D radial ray and the
3-D sphere/ellipsoid/SSE surfaces are genuinely different representations.

## Validation

Two cross-dimensional tests were added:

- `1D3D01` verifies identical common solar-wind caches, apex kinematics, and
  public upstream density/velocity/Parker field in an exactly equivalent
  equatorial geometry.
- `1D3D02` verifies identical local MHD shock states at the spherical apex,
  including shock classification, Mach number, compression, and complete
  upstream/downstream primitive vectors.

Both comparisons use roundoff-level tolerances.  The complete deterministic
validation suite must remain green after the refactor.
