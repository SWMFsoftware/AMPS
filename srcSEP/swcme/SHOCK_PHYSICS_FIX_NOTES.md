# SWCME shock-existence / Rankine-Hugoniot correction

This update replaces the previous fast-Mach-to-hydrodynamic-compression proxy
with a shared ideal-MHD Rankine-Hugoniot solver used by both SWCME 1-D and 3-D.

## Physical changes

1. **Geometric front and physical shock are separate states.**  A finite CME
   surface is a shock only when its normal relative inflow is super-fast.  A
   sub-fast front returns `has_shock=false`, `compression=1`, and no downstream
   jump.  `sheath_comp_floor` no longer creates or strengthens a physical shock.
2. **Full conservative downstream primitive state.**  The solver returns
   downstream density, thermal pressure, vector velocity, and vector magnetic
   field while enforcing the ideal-MHD jump conditions.
3. **Correct immediate 3-D downstream state.**  The Cartesian evaluator now
   approaches the RH state immediately behind the shock instead of returning
   the ambient solar-wind speed at the shock boundary.
4. **Shock state is evaluated at the shock surface.**  The new 3-D
   `shock_state_direction()` API evaluates density and Parker B at `R_shock`;
   the result no longer depends on the radius of an arbitrary field query.
5. **1-D uses the same shock solver.**  The radial direction is the normal and
   the Parker azimuthal component is tangential.  The old compression-floor
   no-shock bug is removed.

## Numerical method

For a candidate compression `r`, mass conservation fixes the downstream normal
shock-frame speed. Tangential momentum and electric-field continuity form a
2x2 linear system. Normal momentum determines downstream pressure. The scalar
total-energy residual is solved on the physical compression interval using a
bracketed bisection search. The trivial `r=1` state is removed from the root
function by dividing the energy residual by `r-1`.

Accepted states must be compressive, have positive downstream pressure, satisfy
the strong-shock compression bound, increase the entropy proxy, and meet the
conservation residual requirements. A super-fast state that does not meet these
conditions reports `solver_converged=false` and is not used as a valid shock.

## New validation

`SHK01` through `SHK12` validate shock thresholding, theta_Bn, parallel and
perpendicular independent limits, the versioned 80-digit full-system oblique
benchmark, mass flux, normal B, tangential electric field, momentum, energy,
entropy/admissibility, and the weak-shock limit. `SHK05` compares every
downstream primitive variable and the evolutionary-fast branch against twelve
independently generated cases; `SHK06` additionally verifies the 3-D production
evaluator's immediate downstream n/V/B state.

The complete current suite has one pre-existing unrelated failure, `CFG02`, due
to the 1-D zero-velocity unit/preparation path converting 0 km/s to 1 m/s. The
shock tests themselves pass.

## Intentionally deferred work

This update does not repair DBM slow-CME/Gamma=0 behavior, centralized
configuration/unit handling, the 1-D ejecta factor bugs, shock-mesh topology,
velocity-divergence treatment, or magnetic connectivity/cobpoint tracking.
