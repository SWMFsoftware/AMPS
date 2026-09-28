# Stage 6: oblique ideal-MHD jump and critical-Mach classification

Stage 6 implements a pure, dependency-light Rankine--Hugoniot solver in
`mhd_jump_solver.{h,cpp}`.  Its inputs are one upstream primitive state, the
outward surface normal, normal front speed, and `gamma_ad`.  It never reads a
mesh or mutable global state.  Failure returns no downstream value, so a
sub-fast front, singular branch, or residual failure cannot publish a partial
shock state.

## Characteristic speeds and conventions

`EvaluateMhdCharacteristics` derives

`c_s^2=gamma_ad p/rho`, `v_A^2=B^2/(mu_0 rho)`,
`v_An^2=B_n^2/(mu_0 rho)`,

and the two magnetosonic roots

`c_f,s^2=0.5[(c_s^2+v_A^2) +/- sqrt((c_s^2+v_A^2)^2-4c_s^2 v_An^2)]`.

Obliquity uses `acos(|B dot n|/|B|)` in `[0,pi/2]`.  The reported upstream
inflow is `V_sh,n-u_1 dot n`; both `M_f` and total-field `M_A` use that same
positive speed.  Normal-Alfvén criticality remains a distinct convention.

## Full oblique jump solve

For a trial compression `X`, mass conservation fixes
`rho_2=X rho_1` and `v_2n=v_1n/X`.  The coupled tangential electric-field and
momentum equations are solved analytically for `B_2t` and `v_2t`; normal
momentum then fixes `p_2`.  The remaining total-energy jump is a scalar
function of `X`.

`SolveObliqueFastShock` scans the physical interval

`1<X<(gamma_ad+1)/(gamma_ad-1)`

for the nontrivial compressive root and refines it with bisection.  It does not
choose a root by proximity to an initial guess.  Every candidate must have
positive downstream pressure and density.  The accepted root must increase
entropy, have upstream super-fast/downstream sub-fast characteristic ordering,
remain below the EOS strong-shock limit, and satisfy independently recomputed
normalized residuals for mass, normal magnetic field, tangential electric
field, normal/tangential momentum, and total energy.

The exact parallel regular branch follows continuously from the oblique
system.  The tangential vector equations are basis-free, so scalar outputs do
not depend on an arbitrary transverse basis.  This baseline does not select a
switch-on branch from an angle cutoff or an Alfvén Mach number shortcut.

## Versioned critical-Mach table

`CriticalMachTable` binds version, checksum, `gamma_ad`, Mach convention, beta
grid, obliquity grid, and row-major values.  A production table must span the
complete physical obliquity interval `[0,pi/2]`; malformed values, incomplete
coverage, EOS mismatch, and convention mismatch are fatal.  Bilinear
interpolation occurs only inside the beta domain--there is no endpoint clamp
or extrapolation.

Beta-domain misses are typed.  For a normal-Alfvén table, exact `B_n=0` is a
separate inapplicable state and never floating-point infinity or false
subcriticality.  `ApplyCriticalityPolicy` implements the three specified
routes: diagnostic-only fast source, fatal preflight, or unrenormalized
budgeted exclusion.

`PreflightCriticalCoverage` accumulates candidate and unavailable area,
incident number rate, and incident kinetic-energy rate using the supplied
physical measures.  It reports each unavailable/candidate ratio, uses a typed
no-candidate-support state for empty history, and fails when any independently
configured budget is exceeded.  Neighboring patches are never renormalized.

## Validation

`make test-stage6` runs the complete Stage 0--5 baseline plus `RH3D01--11`.
The tests cover independent hydrodynamic compression, parallel and
perpendicular limits, a general oblique state, weak/strong limits, rejection
of negative-pressure/sub-fast/expansion inputs, the versioned cold
quasi-perpendicular critical value, EOS and convention guards, near-parallel
continuation, exact-normal-Alfvén inapplicability, all miss policies, and
closure of area/number/energy exclusion ledgers.
