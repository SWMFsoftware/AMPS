# Stage 2: wind, closed plasma, EOS, and interface kernels

Stage 2 contains dependency-free physical kernels. It does not construct the
Stage-3 PFSS--SCS--Parker composite and therefore accepts completed tube
geometry, magnetic magnitude, potentials, or one-sided states as data. This
separation lets every conservation law be verified without a mesh or AMPS.

## Flux-tube and transonic wind

Magnetic flux fixes tube area through `A |B| = Phi_B`. The scalar velocity is
the positive speed along the geometrically outward tube tangent in the rigidly
corotating frame; magnetic polarity is deliberately absent from that scalar.
The effective potential retained by both wind and closed plasma is

`Phi_eff = -GM_sun/r - |Omega cross x|^2/2`.

`SolveRadialIsothermalParker` is the exact spherical verification kernel. With
`y=(u/a)^2` and critical radius `r_c`, it solves

`y - ln(y) = 4 ln(r/r_c) + 4 r_c/r - 3`

by bracketing the subsonic root for `r<r_c` and the supersonic root for
`r>r_c`; at `r_c`, `y=1` exactly. Density follows one independent mass flux,
so no base density is silently introduced. `CheckWindInvariants` independently
checks `rho u A`, the appropriate isothermal/polytropic Bernoulli invariant,
and `p/rho^gamma_w`.

For a supplied variable-area tube, `FindCriticalCandidates` reports every
sampled point with positive
`a_c^2=(d Phi_eff/ds)/(d ln A/ds)` and marks the outer global accelerating
candidate. It never chooses the nearest root or hides the rejected candidates.
`SolvePolytropicTube` enforces `1<gamma_w<3/2`, validates the selected critical
point against both geometrical and thermodynamic values of `a_c^2`, fixes the
mass-flux and Bernoulli eigenvalues there, and solves the two algebraic roots at
every other tube location. It selects the subsonic root inward and supersonic
root outward, failing transactionally if either global branch is absent.
`InvertTargetSpeed` is a bracketed outer solve and publishes nothing on a
non-bracketed or failed tube solution. `ResolveMassLoading` implements
`eta_m=rho u_s/|B|`; the radial record implements
`eta_m=F_m/|B_r|`. A speed alone or zero mapped field is rejected.

These kernels expose the mathematical Stage-2 building blocks. Full composite
tube tracing, transition iteration, and D7 event assets remain owned by the
later stages named in `model/testing_validation.md`.

## Composition-aware EOS

`EvaluatePlasmaFromElectronDensity` requires every ion's abundance relative to
protons, charge number, mass, and temperature. Quasineutrality gives

`n_p = n_e / sum_j(Z_j f_j)`,

and mass density is `n_p sum_j(m_j f_j)` plus electron mass only when the
recorded policy selects it. Pressure is the sum of all ion and electron partial
pressures. Sound, Alfvén, and oblique fast-mode speeds use that same density and
pressure and the independent `gamma_ad`; `gamma_w` is never reused for shock
characteristics.

## Empirical kinematic wind

Positive channel data are represented in `z=ln(y/y_ref)`. For each interval,
`CertifiedPositiveProfile` constructs the unique quintic that reproduces both
endpoint values and their first and second physical derivatives. It recursively
isolates every real stationary point of that quintic and rejects adjacent-node
overshoot. Evaluation returns analytic value, first derivative, and second
derivative and refuses extrapolation.

`BlendTwoZoneWind` uses the exact endpoint-flat weight
`w=10 xi^3-15 xi^4+6 xi^5` in log density. The inner density is authoritative
below `r_a`; the outer density `eta_m |B|/u_out` is authoritative above `r_b`;
one `eta_m` then derives velocity everywhere. The returned log mismatch remains
visible rather than being fit away.

Velocity component and frame are independent discriminants. Radial inertial
speed is divided by a channel-local positive projection only after its guard
passes. Field-aligned inertial speed has `(Omega cross x) dot t_hat` subtracted
exactly once. Field-aligned corotating speed is used directly. The unsupported
radial/corotating pair fails explicitly. `CoverageCensus` reduces physical flux,
area, source, observer, and export measures; raw trace count is never a
coverage authority and event-nominal missing support is fatal.

## Closed-field plasma

The isothermal provider evaluates

`rho=rho0 exp(-DeltaPhi/c_T^2)`, `p=p0 exp(-DeltaPhi/c_T^2)`.

The polytropic provider evaluates

`h=h0-DeltaPhi`,
`rho=rho0(h/h0)^(1/(gamma_c-1))`, and
`p=p0(h/h0)^(gamma_c/(gamma_c-1))`.

It rejects every point with nonpositive enthalpy. `gamma_c` is its own input;
it is neither `gamma_w` nor `gamma_ad`. Two-footpoint normalizations must agree
within the declared relative tolerance and are never averaged.

## Interface balance

For a sharp interface, the normal points from side A to side B and
`w_n=u dot n-V_I,n`. The complete moving-control-surface traction is

`t = rho w_n u + (p+B^2/(2 mu0)) n - (B_n/mu0) B`.

`EvaluateSharpInterface` reports all three components of `Delta t`, its norm,
and signed mass-flux jump. An open--closed separatrix additionally applies the
pointwise `B_n` and interface-relative `w_n` gates. A generic open--open class
boundary does not inherit those separatrix assumptions. Diagnostic policy
keeps mass continuity gating but only reports force imbalance; bounded policy
requires uncertainty and at least three convergence levels; stationary TD is
restricted to solver-produced/imported sharp open--closed states.

For a finite-width representation, `EvaluateVolumeMomentumResidual` returns
the vector

`R = d(rho u)/dt + div(momentum flux) - rho g - f_declared`

in `N/m^3`. This object is intentionally distinct from a surface traction
jump. Tests prove that a missing declared force reappears as the predicted
residual rather than being canceled or hidden by a mean.

## Validation ownership

`WND3D01--20`, `CLS3D01--10`, and `STR3D01` are individually runnable and are
also part of the cumulative Stage-2 gate. They exercise analytic Parker flow,
mass/Bernoulli/entropy closure, critical topology, bracket failures,
composition, normalization invariance, quintic and two-zone construction,
frame/component dispatch, physical-measure coverage, both hydrostatic
closures, moving-interface gates, policy provenance, full vector traction,
and the finite-width volume residual.
