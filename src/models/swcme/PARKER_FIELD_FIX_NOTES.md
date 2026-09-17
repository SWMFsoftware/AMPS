# 3-D Parker magnetic-field correction

This update corrects the 3-D Parker field implementation without changing the
1-D Parker model.

## Production changes

- Added `Params::solar_rotation_axis[3]`, defaulting to global +Z.
- `prepare_step()` normalizes and caches the solar axis as
  `StepState::solar_axis_hat` and rejects a zero/non-finite axis.
- The Parker azimuthal basis is now
  `e_phi = (Omega_hat x e_r)/|Omega_hat x e_r|`.
- The local winding now uses
  `sin(theta_local)=|Omega_hat x e_r|` at every point.
- The cached `StepState::k_AU` is now the latitude-independent equatorial
  coefficient `Omega*AU/V_sw`; the local pitch is
  `k_AU*r_AU*sin(theta_local)`.
- Exact rotation-pole points use the analytical limit `B_phi=0` and do not
  construct an arbitrary transverse direction.
- `Params::sin_theta` is retained for source/API compatibility, but it now has
  one limited role: it specifies the reference colatitude used to normalize
  `B1AU_nT` as total field magnitude at 1 AU. It no longer controls the local
  3-D Parker winding. Its default is now 1.0 (equatorial reference), matching
  the 1-D default and the usual 1-AU ecliptic interpretation.
- The 3-D demonstration programs explicitly define their solar-rotation axis so
  their existing +Z CME/observer geometry remains equatorial.

## Validation added

Three production-facing tests were added under `test/3d/test_parker.cpp`:

- `PAR01`: equatorial vector orientation and 1-AU field normalization.
- `PAR02`: arbitrary latitude and arbitrary solar-axis orientation.
- `PAR03`: north/south pole and near-pole regularity.

All three tests pass against independent analytical references. A separate
finite-difference diagnostic over 175 upstream points produced a maximum
normalized divergence

`|div B|/(|B|/r) = 5.94e-10`,

consistent with a solenoidal Parker field at the tested differencing step.
Formal campaign test `PAR04` should still implement the planned multi-step,
>=1000-point convergence study.

## Full-suite status

The updated test executable builds cleanly with `-Wall -Wextra -Wpedantic` and
both 3-D demos compile. The full current suite reports 6 PASS / 1 FAIL. The only
failure is the pre-existing CFG02 1-D zero-velocity conversion issue
(`0 km/s -> 1 m/s`); the Parker changes do not introduce any additional test
failure.
