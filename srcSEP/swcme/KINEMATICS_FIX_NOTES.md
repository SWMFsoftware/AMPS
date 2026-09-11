# SWCME shared CME/shock-apex kinematics correction

## Scope

This update replaces the separate 1-D and 3-D CME/shock-apex kinematic formulas
with one common implementation in `swcme_kinematics.hpp`.  It addresses the
known slow-CME DBM sign error, the 3-D `Gamma=0` division-by-zero/fallback bug,
and the absence of a controlled observational height-time mode.

No changes in this update are intended to alter the already corrected Parker
field, finite SSE geometry, or ideal-MHD Rankine-Hugoniot shock solution.

## Problems corrected

### 1-D slow-CME clipping

The legacy 1-D code computed `u0=max(0,V0-Vsw)`.  Any CME initially slower than
the wind therefore lost its actual speed deficit and immediately propagated at
`Vsw`.  The common solver retains the signed `DeltaV0`.

### 3-D slow-CME sign error

The legacy 3-D code used `1+Gamma*u0*t` directly.  For negative `u0`, the
denominator decreases and the CME can move farther away from the solar-wind
speed.  The correct drag equation depends on `|DeltaV|` and therefore uses
`1+Gamma*|DeltaV0|*t`.

### 3-D zero-drag failure

The previous 3-D radius contained `log(...)/Gamma`.  With `Gamma=0`, invalid
arithmetic was later hidden by `finite_or(r_sh,r0)`, allowing the shock radius
to remain at its launch value.  `Gamma=0` is now an exact ballistic branch.

## Common sign-aware DBM

For

```text
DeltaV0 = V0 - Vsw
a       = |DeltaV0|,
```

the production solver uses

```text
DeltaV(t) = DeltaV0 / (1 + Gamma a t)
R(t)      = R0 + Vsw t
            + sign(DeltaV0) log(1 + Gamma a t) / Gamma.
```

Fast CMEs decelerate monotonically toward `Vsw`; slow CMEs accelerate
monotonically toward it.  The solution never crosses the background-wind speed.

For `Gamma=0`:

```text
R(t) = R0 + V0 t
V(t) = V0.
```

For very small `x=Gamma*a*t`, the logarithmic distance uses a Taylor series for
`log1p(x)/x`.  This avoids a fragile `0/0` limit while remaining continuous with
the direct `log1p` expression.

## Kinematic modes

The common API supports:

- `Mode::Ballistic` — exact constant-speed propagation;
- `Mode::DBM` — the sign-aware constant-background drag solution above; and
- `Mode::DataDriven` — monotone PCHIP interpolation of radius versus time.

Both `swcme1d::Params` and `swcme3d::Params` include `kinematics_mode`.  The
existing `r0_Rs`, `V0_sh_kms`, `V_sw_kms`, and `Gamma_kmInv` fields remain
source-compatible DBM inputs.

## Data-driven PCHIP

`data_time_s` is strictly increasing and `data_radius_Rs` is nondecreasing.
The common layer converts no units itself beyond its SI contract; the 1-D and
3-D wrappers convert radii to meters when building the common configuration.

The monotone cubic Hermite slopes use the standard weighted harmonic-mean
construction.  The interpolant passes exactly through every knot and cannot
overshoot a monotone radius interval.  Its derivative is returned as the apex
speed, so `R(t)` and `V(t)` are derived from one consistent trajectory.

Default out-of-time behavior is `OUTSIDE_TIME`.  Optional `Ballistic`
extrapolation extends the nearest endpoint using the PCHIP endpoint derivative.
There is no silent cubic extrapolation.

## DBM start radius

The default 1-D and 3-D `r0_Rs` value is now `20 R_s`.  The simple constant-wind
DBM is intended for the drag-dominated heliosphere rather than the low corona.
Callers can still set an event-specific reference radius explicitly.  For
observational event studies, `DataDriven` mode is preferred when suitable
height-time measurements are available.

## Validation

Eight deterministic tests were added:

- KIN01 fast-CME closed form;
- KIN02 slow-CME sign-aware branch;
- KIN03 zero-drag ballistic limit;
- KIN04 small-Gamma continuity;
- KIN05 long-time asymptotic stability;
- KIN06 PCHIP knot exactness;
- KIN07 PCHIP monotonicity/no overshoot;
- KIN08 explicit extrapolation policy and invalid-table rejection.

KIN01-KIN03 also exercise both dimensional wrappers and require 1-D/3-D apex
radius and speed to agree to roundoff.

## Remaining related work

The broader centralized configuration/status API is still incomplete.  The
common kinematics layer already rejects malformed tables and negative drag, but
other SWCME parameters still use legacy clamping rules.  In particular, the
known CFG02 `V_sw=0` 1-D conversion issue is intentionally outside this update.

## Validation status of this update

A clean validation build with `-Wall -Wextra -Wpedantic` registers 35 tests.
All eight new kinematics tests pass, as do the previously corrected Parker,
geometry, density, and shock-physics tests.  The full suite result is currently
34 PASS / 1 FAIL.  The sole failure is the pre-existing CFG02 1-D unit-path
behavior in which `V_sw=0 km/s` is still clamped to `1 m/s`; this is part of the
separate centralized configuration/unit remediation and was intentionally not
changed by the kinematics update.

Both 1-D and 3-D demonstration programs compile cleanly after the update.
A separate 100,000-case randomized DBM stress check spanning fast/slow CME
speeds, ambient speeds, nonnegative drag coefficients and long times produced
no invalid/nonfinite states and no cases in which `|V-Vsw|` increased with time.
