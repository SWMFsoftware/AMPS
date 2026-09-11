# Shock-surface-owned state correction

## Purpose

This update completes the remediation item **"Evaluate shock properties at the
shock surface, not at the query point."**  Earlier SWCME revisions could pass a
Cartesian sampling radius into local shock-strength calculations.  Since the
Leblanc density and Parker magnetic field vary with heliocentric distance, that
allowed one physical shock to acquire different Alfven/fast speeds, Mach
numbers, and compression ratios depending only on where a field was sampled.

The ideal-MHD shock solver introduced in the previous remediation already
centralized most of the corrected physics in
`Model::shock_state_direction()`.  This update removes the remaining secondary
internal paths and makes that function the canonical owner of local 3-D shock
state throughout production diagnostics and mesh generation.

## Canonical physical contract

For a model time `t` and unit direction `u`, local shock physics is evaluated in
this order:

1. `shape_radius_normal()` locates the physical surface radius `R_sh` and
   outward normal `n`.
2. Upstream density and Parker magnetic field are evaluated **at** `R_sh*u`.
3. The self-similar surface motion is projected onto the local normal to obtain
   `V_sh,n`.
4. The common ideal-MHD solver determines shock existence and, when present,
   the complete Rankine-Hugoniot downstream state.
5. The resulting `LocalShockState` is consumed by field evaluators,
   connectivity, diagnostics, and surface-mesh products.

An arbitrary background query radius is therefore not an input to physical
shock strength.

## Production-code changes

### `diagnose_direction()`

The directional convenience diagnostic now calls `shock_state_direction()`
once and projects the requested scalar/vector fields from the returned
`LocalShockState`.  It no longer performs geometry and shock-strength queries
through separate paths.

### `build_shock_mesh()`

Each mesh node is now populated from one canonical `LocalShockState`.  Position,
analytical normal, compression, and normal shock speed therefore come from the
same state used by Cartesian field and connectivity calculations.

### `local_oblique_rc()` compatibility wrapper

The legacy API is retained so external callers continue to compile, but its
historical `r_eval_m`, `Rdir_m`, and `n_hat` inputs do not control the solution.
It recomputes the state from the supplied direction with
`shock_state_direction()` and returns only the legacy scalar diagnostics.
New code should use `shock_state_direction()` directly.

## Validation tests

### SHK13 — Shock-state independence from arbitrary query radius

A nontrivial oblique spherical-shock fixture is used.  The legacy scalar API is
called with radii from `0.5 R_sh` through `5 R_sh`; compression, normal shock
speed, and `theta_Bn` must remain equal to the canonical surface state to
roundoff.  The test also confirms that ambient density at those radii differs
substantially from the surface density, proving that it would detect the old
query-radius leak.

The test additionally samples the Cartesian field at several radii on both
sides of the front and then re-queries the shock state, verifying that field
sampling is side-effect-free and cannot mutate local shock diagnostics.

### SHK14 — Canonical state across diagnostics and mesh

`diagnose_direction()` is compared field-by-field with
`shock_state_direction()`.  A generated shock mesh is then traversed node by
node; each stored radius, normal, compression, and normal shock speed must agree
with a fresh canonical state query in that node direction.

## Expected impact

No validated physical result should change relative to the immediately previous
release because `shock_state_direction()` was already used by the Cartesian
field evaluator and connectivity solver.  This update is an architectural
hardening step: it removes remaining alternate internal shock-strength paths and
adds regression coverage so query-radius-dependent shock physics cannot return.
