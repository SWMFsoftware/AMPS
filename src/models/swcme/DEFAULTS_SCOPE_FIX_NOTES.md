# SWCME canonical defaults and model-scope fix

## Purpose

This update implements remediation item 14: **standardize defaults and
model-scope conventions**.

Before this change, the 1-D and 3-D public parameter structures had evolved
independently.  Even after their physics equations were moved into common code,
default-constructed models could still differ because of hidden parameter
choices (notably the 1-AU density, ejecta thickness, and deprecated sheath
compression floor).  The previous default also selected the optional
FULL_ICME/RESOLVED_COMPRESSION path even though the controlled SEP study is
intended to use an upstream Parker/Leblanc background and an explicit shock
source before shock arrival at the observer.

The correction has three goals:

1. make one header the authoritative source of all shared defaults;
2. make the intended science scope explicit and queryable at each observer;
3. provide deterministic resolved-configuration text for run/validation
   metadata so event-specific overrides are auditable.

## 1. Single shared default set

The new `swcme_defaults.hpp` owns the common baseline.  Both `swcme1d::Params`,
`swcme3d::Params`, and `swcme::core::CommonConfig` refer to these values rather
than repeating literals.

The canonical baseline is:

```text
SWCME_CONFIG_VERSION       2
frame                      HCI_like_inertial
V_sw                       400 km/s
n(1 AU)                    5 cm^-3
|B|(1 AU, reference)       5 nT
T_p                        1.2e5 K
gamma                      5/3
Parker reference sin(theta) 1
solar rotation rate        shared SWCME constant (2.86533e-6 rad/s)
kinematics                 DBM
DBM reference radius       20 R_s
V0                         1500 km/s
Gamma                      1e-7 km^-1
data extrapolation         OUTSIDE_TIME
regions                    SHOCK_ONLY
acceleration               SOURCE
relative source weight     1
```

Inactive FULL_ICME defaults are also common:

```text
sheath thickness @1 AU     0.10 AU
ejecta thickness @1 AU     0.20 AU
shock smoothing @1 AU      0.01 AU
LE smoothing @1 AU         0.02 AU
TE smoothing @1 AU         0.03 AU
sheath ramp power          2
V_sheath,LE / V_sw         1.10
n_ME / n_up                0.50
V_ME / V_sw                0.80
```

The deprecated `sheath_comp_floor` is still present for source compatibility,
but its common neutral default is now `1.0`.  It remains ignored by physical
shock and region calculations; compression comes only from the MHD RH solver.

The former 1-D-specific defaults `n1AU_cm3=6` and
ejecta thickness `0.25 AU` are removed.  The former different 1-D/3-D sheath
floor literals are also removed.

## 2. Default 3-D science geometry

`swcme3d::Params::shape` now defaults to `ShockShape::SSE`, not `Sphere`.
The default half width remains 40 degrees.

A Sun-centered sphere remains available and remains important for analytical
verification and exact 1-D/3-D reduction tests, but it is not the intended
science front because it has infinite angular extent.  Event calculations that
need another geometry must set that choice explicitly and the resolved
configuration manifest records it.

## 3. Parker normalization and polarity convention

The public `B1AU_nT` value is explicitly defined as a **positive total field
magnitude** at 1 AU at the documented reference latitude
`sin(theta_ref)=1`.  The common solar-wind preparation converts that value to
`Br(1 AU)` using the Parker pitch.

The adopted radial polarity is outward (`+1`).  The local 3-D field is still
computed from

```text
sin(theta_local) = |Omega_hat x e_r|,
e_phi            = (Omega_hat x e_r)/sin(theta_local),
```

so the legacy 3-D `Params::sin_theta` member is no longer a spatial physics
parameter.  It is retained for source compatibility only as the reference
latitude at which `B1AU_nT` is interpreted.  The resolved manifest records both
the normalization convention and this reference value.

## 4. Explicit model scopes

`swcme::defaults::ModelScope` defines:

```text
CONTROLLED_SEP_PRE_SHOCK
FULL_ICME_DIAGNOSTIC
INVALID
```

The scope is derived from the existing validated option pair:

```text
SHOCK_ONLY + SOURCE                    -> CONTROLLED_SEP_PRE_SHOCK
FULL_ICME + RESOLVED_COMPRESSION       -> FULL_ICME_DIAGNOSTIC
anything else                          -> INVALID
```

It is deliberately **not** another independent input.  This avoids a new failure
mode in which metadata says "controlled SEP" while the velocity field actually
contains a resolved shock compression.

The canonical default is `CONTROLLED_SEP_PRE_SHOCK`.

### Observer-local validity

Both model classes expose `observer_scope_status(...)`.

For the controlled SEP scope, the Parker/Leblanc background is declared in
scope only while the modeled shock has not reached that observer.  Equality is
counted as shock arrival; the stated scope is pre-shock, not "up to and
including the shock".

The 3-D implementation does not compare only with the apex radius.  It queries
the same production surface geometry used by shock diagnostics and
connectivity.  Consequently, an observer outside a finite SSE cap remains
pre-shock/in-scope even after the CME apex has passed that heliocentric radius.

`FULL_ICME_DIAGNOSTIC` remains mathematically evaluable across its modeled
regions and therefore reports itself within its declared diagnostic scope.  The
name is intentional: this does not promote the phenomenological sheath/ejecta
model to a validated global ICME background, and the ejecta magnetic field
remains Parker-like.

## 5. Resolved-configuration manifests

Both public namespaces now provide:

```cpp
swcme1d::resolved_configuration_manifest(const swcme1d::Params&)
swcme3d::resolved_configuration_manifest(const swcme3d::Params&)
```

The format is deterministic plain `key=value` text.  It includes:

- configuration/schema version;
- model dimensionality and frame;
- derived model scope;
- Parker normalization and polarity convention;
- all ambient and kinematic public inputs;
- data-driven table counts and every table entry;
- region and acceleration choices;
- all source, smoothing, sheath, and ejecta parameters;
- deprecated compatibility parameters;
- in 3-D, shock shape, axis ratios, half width, CME direction, solar axis and
  solar rotation rate.

Inactive fields are deliberately still serialized.  A run manifest must show
an event override even if the selected mode means that field does not affect the
current calculation.  This makes configuration reviews and later regression
hashes reproducible.

The Fix 15 campaign/AMPS layer now embeds this block directly into its larger
run manifest together with git hash, compiler flags, MPI/OpenMP information,
random seeds, and the SEP spectrum/normalization configuration.

## 6. Validation

Four deterministic regression tests were added:

### DEF01 — shared default equivalence

Compares every common 1-D/3-D default, including ambient, kinematic, source,
region, smoothing, sheath and ejecta values.  It also checks the common 1-AU
density and DBM reference radius against `swcme_defaults.hpp`.

### DEF02 — science-scope and convention guard

Checks:

- default `SHOCK_ONLY + SOURCE`;
- derived `CONTROLLED_SEP_PRE_SHOCK` scope;
- default 3-D finite SSE geometry;
- stable HCI-like inertial frame label;
- outward radial magnetic polarity;
- explicit total-|B| Parker normalization convention.

### DEF03 — observer-local scope gating

Exercises both dimensional models with a ballistic front.  It verifies:

- pre-shock observer -> in scope;
- shock-at/past observer -> out of controlled scope;
- observer outside a finite SSE cap -> still in scope even after apex passage.

### DEF04 — resolved-configuration metadata

Checks that common and 3-D-specific keys exist, explicit event overrides appear,
derived scope is written, data tables are represented, and repeated
serialization is byte-identical.

Tests that are specifically about the optional FULL_ICME path now set both
`region_mode=FullICME` and
`shock_acceleration_mode=ResolvedCompression` explicitly.  They no longer rely
on an unrelated global default to select the physics they intend to test.

## 7. Compatibility notes

This update intentionally changes **default behavior**:

- 1-D default density changes from `6` to `5 cm^-3`;
- 1-D default ejecta thickness changes from `0.25` to `0.20 AU at 1 AU`;
- 3-D default geometry changes from `Sphere` to finite `SSE`;
- both models default to `SHOCK_ONLY + SOURCE` rather than
  `FULL_ICME + RESOLVED_COMPRESSION`;
- the deprecated sheath-floor default is neutralized to `1.0` in both models.

These are deliberate science-convention changes.  Existing event studies that
require old defaults should set them explicitly; the resolved manifest will
then make that choice visible.  Production formulas, RH physics, connectivity,
mesh topology, and divergence treatment are otherwise unchanged by this fix.
