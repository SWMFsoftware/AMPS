# SWCME shock-acceleration representation and smoothing fix

## Purpose

This update implements remediation item 10: **define one shock-acceleration
representation and make shock smoothing identical in 1-D and 3-D**.

The underlying problem is physical, not cosmetic.  In a focused-transport
calculation, first-order shock acceleration can be represented either by a
prescribed DSA-consistent source or by resolving a compressive velocity gradient
and retaining the transport equation's `div(V)` momentum-change term.  Applying
both to the same particle population double-counts the same shock acceleration.
Earlier SWCME revisions had no production-level switch that prevented this.
They also allowed 1-D and 3-D to treat the shock edge differently.

## 1. One authoritative acceleration mode

A new common header, `swcme_acceleration.hpp`, defines

```cpp
swcme::acceleration::Mode::Source
swcme::acceleration::Mode::ResolvedCompression
```

There are no independent source/compression enable booleans.  The enum is the
single source of truth.

`Params::shock_acceleration_mode` is available in both dimensional interfaces.
Fix 10 originally retained `ResolvedCompression` as the compatibility default.
Fix 14 subsequently standardized the science baseline: the canonical default is
now `SOURCE` together with `SHOCK_ONLY`.  `ResolvedCompression` remains the
explicit `FULL_ICME_DIAGNOSTIC` alternative.  The mutual-exclusion physics and
smoothing implementation described in this note are unchanged.

## 2. Safe mode/region combinations

Centralized validation in `swcme_config.hpp` currently accepts only:

```text
SOURCE               + SHOCK_ONLY
RESOLVED_COMPRESSION + FULL_ICME
```

`RESOLVED_COMPRESSION` also requires
`edge_smooth_shock_AU_at1AU > 0`.

This strict pairing is deliberate.  SOURCE/SHOCK_ONLY guarantees that the
transport-facing velocity contains no RH jump that could create a second DSA
process through `div(V)`.  RESOLVED_COMPRESSION/FULL_ICME guarantees that the
selected compression mechanism actually exists in the transport background.
Unsupported mixed configurations fail during setup with
`INCOMPATIBLE_OPTIONS`; the evaluator never guesses how to repair them.

## 3. SOURCE representation

When a physical fast shock exists, the common acceleration record contains:

- physical/surface flags;
- time;
- shock position and normal;
- normal shock speed;
- density compression ratio;
- `theta_Bn`;
- fast Mach number;
- upstream density and `|B|`;
- DSA phase-space momentum slope

```text
q = 3 r_c / (r_c - 1);
```

- a dimensionless relative source weight per unit sampled area.

`source_enabled=true` and `resolved_compression_enabled=false`.

At the Fix-10 layer the relative weight is deliberately not assigned physical
injection-rate units.  Fix 15 now supplies the AMPS-facing
`swcme_sep_source.hpp` contract: its default `RELATIVE_ONLY` mode preserves this
dimensionless controlled source, while optional
`REFERENCE_DIFFERENTIAL_INTENSITY` requires an explicitly declared physical
`J(E_ref)` and performs the documented unit conversion.  No empirical injection
efficiency is inferred implicitly.

In SOURCE mode the validated region mode is SHOCK_ONLY, so the field evaluator
returns the analytical Parker/Leblanc wind on both sides of the mathematical
shock surface.  The shock still exists for geometry, connectivity, and source
bookkeeping.

## 4. RESOLVED_COMPRESSION representation

When RESOLVED_COMPRESSION is selected, the prescribed DSA source is disabled:

```text
source_enabled = false
resolved_compression_enabled = true
dsa_q_phase_space = unavailable
relative_source_weight_per_area = 0
```

The exact discontinuous Rankine-Hugoniot state is still computed by the common
shock solver and remains available from the shock diagnostic API.  The
transport-facing field, however, represents the jump with a finite C1 layer.

For local shock radius `R_sh` and total width `w_sh`, the transition spans

```text
R_sh + w_sh/2   exact analytical upstream endpoint
R_sh            transition midpoint
R_sh - w_sh/2   exact RH downstream endpoint
```

The blend is the shared cubic smoothstep

```text
S(q) = q^2 (3 - 2q),    0 <= q <= 1.
```

Thus `S'(0)=S'(1)=0`.  Density, velocity, and magnetic field use the same blend
weight and the same geometrical endpoints.

## 5. Consistent sheath join

The phenomenological sheath no longer begins its relaxation at the mathematical
shock center when a resolved layer is present.  Its progress coordinate begins
at the **inner edge** `R_sh-w_sh/2`, where the state is exactly the RH downstream
state.  The sheath's own smoothstep also has zero derivative at its start.
Consequently the resolved shock and sheath match in value and first derivative
without a second numerical kink.

The total shock width is self-similar:

```text
w_sh = edge_smooth_shock_AU_at1AU * R_sh(local)
```

CFG03 validates this request at setup and rejects values above 90% of the
sheath thickness; boundary construction no longer silently caps it.  Every
accepted value is therefore preserved exactly as a self-similar fraction.  The
finite SSE/ellipsoid flanks use the same dimensionless smoothing as the apex;
they do not inherit an apex-sized absolute width.

## 6. Shared acceleration state

`swcme::acceleration::ShockAccelerationState` is the current common acceleration
contract.  It is intentionally smaller than the future AMPS `SEPSourceState`.
Its role is to make the physical representation auditable now and to provide a
stable 1-D/3-D comparison target.

Public wrappers are:

```cpp
swcme1d::Model::shock_acceleration_state(step)

swcme3d::Model::shock_acceleration_state(step, direction, state)
```

Both wrappers call the same `swcme::acceleration::make_state()` decision logic.

`serialize_csv()` provides deterministic scientific-notation serialization for
regression comparisons.  It is an audit/test format, not the final event-data
or AMPS interchange format.

## 7. Validation tests

Five new acceleration tests were added.

### ACC01 — SOURCE source / SHOCK_ONLY flow

Checks that a physical shock enables the explicit source, disables resolved
compression, returns the correct DSA phase-space slope, and leaves transport
fields analytical on both sides of the source surface.

### ACC02 — mode mutual exclusion

Checks rejection of:

- SOURCE + FULL_ICME;
- RESOLVED_COMPRESSION + SHOCK_ONLY;
- zero resolved-shock smoothing width;
- negative relative source weight.

### ACC03 — resolved C1 shock profile

Checks the finite width, exact upstream outer endpoint, exact RH downstream
inner endpoint, interior midpoint, and derivative matching at both ends.

### ACC04 — 1-D/3-D smoothing identity

Uses an equivalent spherical +X/equatorial fixture and samples five normalized
positions through the resolved shock.  Density, radial velocity, radial magnetic
field, and Parker-azimuthal magnetic field must agree to roundoff.

### ACC05 — resolved mode has no prescribed DSA source

Checks source-disable/resolved-enable flags, unavailable active DSA slope, zero
source weight, and deterministic serialization.

### 1D3D03 — SOURCE record identity

Compares every common source quantity between equivalent 1-D and 3-D models and
requires byte-identical serialized records.

## 8. What this fix does not do

This remediation does **not** yet define the final AMPS-facing source spectrum,
absolute injection-rate units, energy-grid conventions, or source-particle
sampler.  Those belong to the later dedicated SWCME-to-SEP interface task.

It also does not complete the separate `div(V)` modernization.  The important
contract established here is that only RESOLVED_COMPRESSION is allowed to put a
shock compression profile into the transport flow; SOURCE/SHOCK_ONLY cannot
silently reintroduce one.
