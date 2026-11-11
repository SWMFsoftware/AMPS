# Conservative analytical material transport

Date: 2026-11-10. Status: **local transport kernel implemented and independently
tested; a production lateral map and coupled regional closure remain open**.
The user selected analytical/semi-empirical cross-ray transport without an MHD
runtime. `material_transport.{h,cpp}` implements that choice's universal
material identities. It does not invent the missing velocity/map closure.

## State, frames and equations

Use reference material coordinates a in HCI metres and one differentiable
map X(a,t), also in Sun-centred inertial HCI metres:

```text
F_ij = partial X_i / partial a_j     (dimensionless)
J = det(F) > 0                     (dimensionless)
v = partial_t X at fixed a          (m/s)
L = grad_X v = Fdot * inverse(F)    (1/s)
D J/Dt = J * trace(L).
```

For a regular ideal material parcel with no mass source, constant composition,
constant gamma, no resistive magnetic transport and adiabatic thermodynamics:

```text
rho = rho0/J                           (kg/m^3)
B = F*B0/J                             (tesla)
p = p0*J^(-gamma)                       (Pa)
D rho/Dt = -rho*trace(L)                (kg/m^3/s)
D p/Dt = -gamma*p*trace(L)              (Pa/s)
D B/Dt = Fdot*B0/J - B*trace(L)         (T/s).
```

For an actual compatible map, a material face transforms as
`dS = J*inverse(F)^T*dS0`. Hence `B.dS = B0.dS0`, including its sign and all
orientations, and `rho*dV=rho0*dV0`. A divergence-free reference field remains
divergence-free under the same compatible map. These are full Cartesian
Cauchy identities; diagonal-only field scaling or passive radial flux cannot
replace them during lateral transport.

The public `MapJet` supplies X, v, acceleration, F and Fdot. These must come
from one map. A pointwise algebraic API cannot prove global integrability,
invertibility, coverage, reference-field divergence or event ownership.
Those require actual map/inverse/mesh tests before a production provider may
publish this state. No per-particle coordinate, weight, source or PIC/MPI type
enters the public contract.

Pressure here is explicitly ADIABATIC. Retaining an isothermal ambient with
nonzero compression requires its independently specified heating/work law;
these four tests do not establish unforced ambient thermal preservation.
The kernel does not solve momentum balance, choose shock states, supply
entropy generation at a shock, or qualify force/heating/work residuals.

## Numerical algorithm and rejection

Normalize F by its largest absolute entry before computing cofactors and the
determinant in long double. Check orientation and a determinant roundoff bound
based on the actual six signed triple-product magnitudes before constructing
the inverse. Compute the physical J with the scale restored. An unresolved or
negative determinant returns typed numerical failure. There is no map repair,
rank regularization or clipping.

Evaluate pressure with `exp(log(p0)-gamma*log(J))` to avoid an unnecessary
intermediate overflow in J^gamma. Density, pressure, field, energy, gradient
and material rates must remain finite and physically positive where relevant.
Unrepresentable output rejects. Zero B is valid and stays zero.

A resolved positive J much smaller than one is valid algebraically. The
frozen PRODUCTION minimum-J gate remains separate: a kernel unit test of
J=1e-36 does not qualify that deformation for a CME. No physical floor or
force/work tolerance was changed in this implementation.

## Independent references and commands

From `src/models/sep_corona_swcme`:

```sh
make -j16 build/test_material_transport
build/test_material_transport
```

This is a dependency-light shared-library test, not a native AMPS build.
The four globally registered identities have prefix `CSWC-FULL-`:

| Local case | Independent what/how/why |
| --- | --- |
| `CSWC0628-XRAY-MATERIAL-INVARIANTS` | Closed-form triangular affine deformation, direct deformed-face cross products, positive/negative flux orientations and exact mass/entropy identities. Prevents a diagonal-only or unsigned magnetic update. |
| `CSWC0628-XRAY-MATERIAL-RATES` | Differentiate a separate closed-form transported-field oracle in time and an independent triangular-inverse Eulerian velocity in space. Checks continuity, adiabatic/induction rates and all nine velocity-gradient components. |
| `CSWC0628-XRAY-MATERIAL-NONLINEAR` | Nonlinear lateral map X=(a_x,a_y*exp(k*a_x*t),a_z), its explicit inverse, exact transported field and spatially varying J; 24 physical samples check mass and divergence. |
| `CSWC0628-XRAY-MATERIAL-REJECTION` | Folded/unresolved maps, overflowed J and nonfinite rates reject. Distinguishes a resolved small J from the separate production gate and verifies zero B has no floor. |

For the nonlinear reference `k=0.2/(m s)`, `t=2 s`, a known initially uniform
field B0=(1,-2,0.5) T becomes

```text
Bx = exp(-k*x*t)
By = k*t*y*exp(-k*x*t)-2
Bz = 0.5*exp(-k*x*t).
```

Analytically `partial_x Bx+partial_y By+partial_z Bz=0`. Independent spatial
finite differences give maximum `1.24345e-10 T/m` on this manufactured metre-
scale fixture. This is retained numerical error, not a coronal field floor
or a substituted production divergence threshold. The optimized/debug and
ASan/UBSan/LeakSanitizer controls pass; receipts are under
`test_output/remaining-failures-repair-20261110/material-*`.

## Required next regional work

The kernel does not close the model by itself. Next steps must supply:

1. A complete analytical lateral velocity/material-map law, with normal contact
   motion from the selected corotating geometry and compatible shock-side
   states. Merely choosing a desired divergence is not a complete momentum or
   thermodynamic closure. It must match the actual 3-D ambient at initialization
   and undisturbed boundaries without invented radial mass sources.
2. A common curved material geometry and inverse, positive maps, disjoint
   retained/admitted volumes and reference inventory. Material flow across
   fixed observation rays is real transport, while reference labels and mass
   ownership remain consistent. Do not infer finite-cell boundaries solely
   from legacy ray centres/weights without an explicit checked construction.
3. EOS, shock entropy and explicit prescribed force/heating/work. Compare
   one-sided smooth residuals and interface jumps independently; numerical
   exactness of mass/flux does not make a map an MHD solution.
4. Original computed-shock branch, mass/time/ray refinement, contact leakage,
   minimum-J, signed flux, force/work, supported-area and execution gates.
   Preserve the failed radial fixture and distinct expanded identity.
5. Actual BG3D-5 matching and later one-/four-rank, zero-particle native G1
   publication/handoff/1-AU evidence, only after preceding gates pass.

The production lateral law is still unselected/unimplemented. This document
states that gap explicitly rather than filling it with a guessed ambient-only
transverse velocity, a fitted density multiplier or a new MHD runtime.


The next design investigation is in
[ANALYTICAL_LATERAL_TRANSPORT_DESIGN.md](ANALYTICAL_LATERAL_TRANSPORT_DESIGN.md).
It first traces the already prescribed ambient material flow, then tests a
candidate composed deformation for representability of actual RH/contact
states before choosing a production ansatz. No scalar embedding or guessed
transverse velocity is selected merely because Cauchy identities pass.
