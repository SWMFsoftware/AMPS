# Stage 3: finite-SCS, conservative Parker mapping, and plasma coupling

Stage 3 adds the first complete dependency-light source-surface coupling
kernel.  The implementation is in `source_surface_coupling.{h,cpp}` and is
used through typed records: magnetic sector, interface side, coupling mode,
mapping validity, calibration role, and budget validity are never represented
as interpolated floating-point flags.  AMPS mesh storage and halo exchange
remain adapter responsibilities for Stage 8.

## Finite-shell SCS field

`FiniteShellScs` solves `B=-grad(Psi)` on `R_i<=r<=R_scs` using the same real,
orthonormal Condon--Shortley harmonics as `PfssHarmonics`.  For each retained
`l>=1` mode it evaluates

`Psi_lm=c_lm[r^l-R_scs^(2l+1)r^(-l-1)]Y_lm`,

with `c_lm` fixed by the requested unsigned normal field at `R_i`.  This makes
the outer field radial analytically.  Degree zero is not discarded: the
unsigned boundary has nonzero net flux, so the implementation uses
`d_00=h_00 R_i^2` and fixes the harmless potential gauge with
`Psi_00(R_scs)=0`.  The factory rejects an absent or nonpositive monopole,
duplicate modes, malformed real coefficients, and non-finite radii.

`ScsHarmonicAttenuation` implements the exact finite-shell factor

`G_l=(2l+1)x^(l+1)/[l+(l+1)x^(2l+1)]`, `x=R_scs/R_i`,

and `EvaluateScsSpectrum` reports inner/outer non-monopole power and the
absolute outer zonal fraction.  A zero-thickness shell is rejected; it cannot
be advertised as radializing the field.

## Sector and interface semantics

`RestoreSector` applies only the categorical sign after the continuous
unsigned field has been evaluated.  A query exactly on an ideal HCS requires
an explicit side.  Cross-sector or HCS-drift transport is rejected unless a
separately qualified finite-sheet provider exists.  The pure sign reversal
therefore cannot manufacture a small magnetic magnitude.

`EvaluateMagneticInterface` requires the conservative normal trace to close.
It then reports the physical tangential kink and surface current
`K_s=n cross (B_scs-B_pfss)/mu_0`.  Production requires a resolved transition;
the sharp and no-SCS alternatives are verification-only contracts checked by
`ValidateCouplingMode`.

The resolved transition uses the endpoint-flat quintic `chi` and evaluates

`B=(1-chi)B_pfss+chi B_scs+grad(chi) cross (A_scs-A_pfss)`.

The last term is intentionally present: directly blending magnetic vectors
would not be the curl of the blended potential and would generally violate
the solenoidal construction.  Signed endpoint potentials must be flux
balanced, share the fixed Mie gauge, and have compatible tangential sheet
traces before this operation is legal.

## Rotating-footpoint map and plasma continuation

`IntegrateLongitudeMap` advances the characteristic equations

`dPhi/dr=-K`, `dA_phi/dr=-(partial K/partial phi)A_phi`

with RK4.  It publishes `J_phi=1/A_phi` only after every step remains above
the configured positive floor; a fold is an error, not a value to clamp.
`MapParkerState` then uses the same source label and Jacobian for both magnetic
and mass flux:

`B_r=(R_scs/r)^2 B_r0 J_phi`,
`B_phi=-r sin(theta) K B_r`,
`rho=rho0(u_r0/u_r)(R_scs/r)^2 J_phi`.

Consequently the field-aligned rotating-frame speed and the mass flux cannot
come from a second, inconsistent Parker prescription.

## Plasma sheet, D9, D6, and consumer budgets

The empirical plasma-sheet factor is applied to exactly one normalization
authority.  Base-density normalization scales density and pressure together,
preserving fixed temperature; mass-per-flux normalization changes only the
loading target.  Invalid contrast, width, or thermodynamic rules fail before
a wind solve.

`EvaluateTransitionDiagnostics` builds the one-sided D9 scalars from a
conservative normal mortar trace.  It rejects weak fields or a failed normal
trace, distinguishes supported antipodality from a zero-jump inapplicable
state, and retains absolute/net crossing flux and incidence cosines.

`BuildOpenFluxLifecycle` enforces the two-pass D6 sequence: derive one scale
from the construction asset, rebuild Pass B, and qualify that rebuilt field
against a distinct checksum.  Construction and qualification data cannot be
the same asset.  `EvaluateBudgetRatio` stores the physical numerator and
denominator and gives a zero denominator a typed inapplicable state rather
than fabricating a passing zero.  Point observers likewise receive a
valid/rejected-clearance state and never a fictitious magnetic-flux fraction.

## Validation

The cumulative `make test-stage3` gate runs the Stage 0--2 baseline plus
`SCS3D01--09`, `HCS3D01--03`, `CPL3D01--09`, `CPL3D11--12`,
`PLS3D01--04`, `OFX3D01`, and `LOS3D01`.  These tests cover inner/outer shell
boundary conditions, the analytic attenuation law, sector discontinuities,
surface current, conservative mapping and inverse closure, fold rejection,
transition endpoints/cross term, signed-potential qualification, plasma-sheet
authority, two-pass calibration provenance, and typed clearance budgets.
