# Conservative cross-ray extension of the regional CME model

Date: 2026-11-10. Authorization: the user resumed full-CME work and selected
retaining the three-dimensional ambient while developing conservative
cross-ray transport. This supersedes the fixed-ray/no-lateral-transport
restriction only for the bounded regional extension described here. Existing
legacy/reduced providers, particle behavior and G2/G3/G4 boundaries are
preserved.

Status: implementation in progress, **not a qualified regional provider**.
The user subsequently selected the analytical/semi-empirical evolution path,
preserving the prohibition on an MHD runtime; see
[CROSS_RAY_EVOLUTION_SCOPE.md](CROSS_RAY_EVOLUTION_SCOPE.md). HLL, static CT
and fluid-update kernels are conservation references, not a selected production
MHD time integrator. Production lateral transport will use a compatible
material map and independently qualified analytical closure.
Implemented bounded increments are the face-flux/transfer kernel, exact
static sector geometry, immutable shared face/edge topology, and static
face-flux constrained transport (CT). Those kernels plus cell-field reconstruction, atomic fluid/CT epoch coupling
and the user-selected corotating threaded contact have 25 independent cases.
The subsequent full Cauchy material-transport kernel adds four independent
cases; see [MATERIAL_TRANSPORT.md](MATERIAL_TRANSPORT.md).
There is no complete MHD integrator, ideal upwind edge EMF, moving mesh,
well-balanced ambient, piston coupling, or native application yet.

## Physical and numerical contract

The regional fluid is ideal, single-velocity MHD with the existing configured
composition and constant adiabatic gamma. All quantities use SI and Sun-centred
HCI Cartesian components. The ambient authority and its observed closed/open
transitions remain unchanged. No radial mass source is invented to compensate
for omitted lateral transport.

Conservative densities are

```text
U = (rho, rho*v, E, B),
E = p/(gamma-1) + rho*|v|^2/2 + |B|^2/(2*mu0).
```

For a unit shared-face normal n directed from left to right, define
`vn=v.n`, `Bn=B.n`, and `Pt=p+|B|^2/(2*mu0)`. Stationary-face fluxes are

```text
F_rho = rho*vn,
F_momentum = rho*v*vn + Pt*n - (Bn/mu0)*B,
F_energy = (E+Pt)*vn - (B.v)*Bn/mu0,
F_B = vn*B - Bn*v.
```

For a face with physical normal speed `wn`, Reynolds transport gives the ALE
flux `G=F-wn*U`. Thermodynamic states and total energy remain in HCI; merely
subtracting wn from v everywhere would omit the required energy/momentum
transformation. In particular the normal component of the ALE B flux can be
`-wn*Bn`, even though the stationary induction flux has zero normal component.
Magnetic flux on the actual moving face must ultimately use its own geometric
and constrained-transport update.

The initial Riemann implementation is the HLL two-wave flux. Wave estimates
use normal fluid velocity relative to the face and the full canonical MHD fast
speed, not the passive-radial reduced characteristic. Its central branch is

```text
G_HLL = (sR*G_L - sL*G_R + sL*sR*(U_R-U_L))/(sR-sL).
```

It is a diffusive initial discretization, not a claim of exact rotational,
contact, or compound-wave resolution. Invalid states, invalid geometry,
unresolved/nonfinite waves, inconsistent face-normal magnetic fields, and
nonpositive candidate pressure/density reject with typed status. There are no
density/pressure/field floors and no unconditional positivity claim outside
the validated timestep/state domain.

A single evaluated face flux is transferred with opposite signs to both
neighbors. Every Cartesian momentum component and the full total energy are
included; no radial-only pressure update or unledgered transverse work is
permitted. The first transfer API supports stationary fixed volumes only.
It explicitly rejects moving-face use until geometric conservation and face
motion are implemented and independently tested. This prevents an ALE flux
unit test from pretending to validate a moving mesh.

## Ordered implementation and qualification gates

| Gate | Required implementation and independent evidence | Current state |
| --- | --- | --- |
| XR0 | Flux/state conversion; exact constant-state and ALE flux checks; pairwise conservation, positivity rejection and advection convergence | Kernel implemented; five independent cases PASS. Regional coupling remains open. |
| XR1 | One common angular/radial mesh and shared-face quadrature, overlap/remap ownership, volume/area closure and geometric conservation | Static geometry and shared topology implemented; six independent cases PASS. Motion/remap, piston cuts and legacy-ray conversion remain open. |
| XR2 | Consistent face magnetic flux and edge electric fields; constrained transport, divergence preservation and magnetic-energy/interface budgets | Static CT, cell-field reconstruction and atomic fluid/CT energy coupling implemented and independently tested. Ideal-MHD edge EMF, actual flux/time quadrature and moving-face coupling remain open. |
| XR3 | Conservative multidirectional evolution with all transverse momentum/energy terms and compatible source/work treatment | Atomic static conservative update implemented; four independent cases PASS. Production analytical material transport and its force/heating/work closure remain open; no MHD runtime is selected. |
| XR4 | Actual 3-D ambient continuity and unforced-reference preservation, including topology/interface budgets and rotating-frame consistency | Not implemented |
| XR5 | One common piston/contact boundary and actual computed shock; original mass/time/ray convergence gates, residuals and valid supported-area budgets | Not implemented |
| XR6 | Region-safe publication and qualified shared assembly before the remaining BG3D-5--9/native steps | Not implemented |

These are prerequisites within the resumed BG3D-4 correction. Their labels do
not replace the existing BG3D or CSWC gate identities. The old per-ray fixture
and original failing receipts remain retained. No old test is relabeled PASS
from XR0 evidence, and no physical/convergence tolerance or execution ceiling
is changed here.

The expanded closure needs a distinct immutable/fingerprinted identity before
production selection; legacy `per-ray-lagrangian-piston-v1` cannot silently
mean a different transport model. A model selector, actual coupling, and input
update are not implemented by the initial kernel. All original initial states,
ray coverage, refinement controls and comparison thresholds remain the
acceptance references for subsequent integrated work unless the user explicitly
approves a scientifically justified physical-contract change.

## Magnetic and geometric completion requirements

Conserving the volume integral of cell-centred B across one face does not prove
`div B=0`. Shared normal magnetic flux, edge circulation, changing face areas,
and discrete geometric closure need their own authority and tests. No cleaning
or component clipping may substitute for that evidence. The later ALE geometry
must satisfy uniform-state preservation and the discrete geometric conservation
law, with equal-and-opposite flux across matching overlap faces.

The regional method evolves only the prescribed CME/background domain. It is
not an implicit change to generic PIC, a global solar-eruption MHD model, or
permission to add particle feedback. Native allocation/injection remains zero
for the required G1 campaigns.

## References and validation provenance

The numerical flux follows the integral two-wave construction of
[Harten, Lax and van Leer (1983), SIAM Review](https://epubs.siam.org/doi/10.1137/1025002).
The need to distinguish conservative cell updates from face-based magnetic
transport on static/moving meshes is described by
[Mocz, Vogelsberger and Hernquist (2014), sections 2.1--2.3](https://arxiv.org/html/1402.5963).
These references motivate the method; they do not qualify this implementation.

Independent tests must check hand-evaluated SI fluxes, analytic transported
profiles, coordinate rotations, both upwind branches, zero-field limits,
physical canonical waves, conservation and deliberate defects. Published
shock/wave references and multidimensional magnetic-loop tests are required
before a production coupling is accepted. Retain exact compiler/input/binary
identity, runtime and numerical versus physical residual budgets.

## Current reproducible bounded evidence (2026-11-10)

The model Makefile and global `CORE-MODEL-SEP-CORONA-SWCME-FULL-001`
owner register the same binaries and case markers. From this model directory:

```sh
make -j16 build/test_cross_ray_transport build/test_spherical_sector build/test_cross_ray_mesh build/test_cross_ray_ct
build/test_cross_ray_transport
build/test_spherical_sector
build/test_cross_ray_mesh
build/test_cross_ray_ct
```

These are dependency-light shared-library tests, not native AMPS builds.
Inputs are explicit SI radial/mu/longitude edges, a right-handed HCI basis,
a positive resource budget and an owning event fingerprint. No ray boundaries
are inferred. All 336 reference cells and 792 internal faces are checked.
An independent eight-node Gaussian integration of `A=(B cross X)/2` supplies
initial flux, and a prescribed analytic curl-E supplies the Faraday oracle.
The retained CT run reports maximum raw divergence `3.20004e-20 T/m` after
its manufactured/random-edge campaign. Arithmetic bounds explain floating-
point cancellation; they do not replace the separate regional physical gates.
Evidence: repository-local `test_output/remaining-failures-repair-20261110/`,
`cross-ray-sector-controls.json`, `mesh-controls-current.json` and
`ct-controls-current.json`, with corresponding build/run logs.

## Field reconstruction, fluid epochs and contact startup

`ReconstructCellField` fits a constant HCI field to actual CT face flux using
physical curved-face areas. It preserves the authoritative flux/divergence,
retains fit/conditioning diagnostics, and rejects deficient geometric rank.
An independently integrated affine divergence-free field gives relative L1
errors 0.0152182, 0.00391239 and 0.000984985 on 4/8/16 cells per axis.
This is a bounded consistent reconstruction, not a higher-order field map.

`StaticFluidState` couples every shared mass/vector-momentum/total-energy
transfer to one new CT epoch. It records signed boundary and prescribed
volume force/power ledgers. Thermal pressure uses the NEW reconstructed B
and the CONSERVED E; no magnetic-energy reset or density/pressure floor is
allowed. Failure rejects the whole private epoch. Its face inputs are actual
area-integrated, time-averaged fluxes; supplying those integrals and ideal
edge EMF remains the qualified integrator's responsibility. Manufactured
external-curl-E controls exercise explicit electromagnetic work and are not
an ideal-MHD evolution qualification.

The user's explicit startup/contact choices and the implemented differentiated
geometry are documented in [CROSS_RAY_CONTACT_STARTUP.md](CROSS_RAY_CONTACT_STARTUP.md).
The full ellipsoid initially corotates, tapers angular rate to zero over
T=1200 s, and retains the resulting 0.0985-degree final offset. The selected
continuous-field ideal material contact retains threaded magnetic flux.
Manufactured two-sided tests do not supply actual BG3D-5 ejecta matching.
No old contact asset or selector is silently changed by these additions.

Additional exact shared-library commands:

```sh
make -j16 build/test_cross_ray_reconstruction build/test_cross_ray_fluid build/test_cross_ray_contact
build/test_cross_ray_reconstruction
build/test_cross_ray_fluid
build/test_cross_ray_contact
```

Optimized and debug execution passes all first 21 cross-ray cases. ASan,
UBSan and LeakSanitizer pass all six transport/geometry/CT/reconstruction/
fluid binaries with every new translation unit instrumented; the existing
coronal dependency remains the production uninstrumented archive. Contact
rotation has its own subsequent optimized receipt and profile checks.
`model-owner-current/amps-core-tests.json` retains the complete 38-row model
owner plus 13-row harness snapshot: 50 PASS, no FAIL/SKIP, one ERROR, 540.165 s.
That ERROR is the unchanged 180 s refinement ceiling. The later four contact
cases extend the registry to 42 model rows and need their own current-owner
receipt. Do not relabel the 51-case historical snapshot as 55 cases.

A separate smooth actual-ambient diagnostic checks 96 spacetime stencils on
four rays, eight radii and three step sizes. Maximum residuals normalized by
rho*abs(v)/r, abs(B)/r and abs(v)*abs(B)/r are respectively 3.491e-6, 6.945e-11
and 2.362e-7. It includes time variation and transverse derivatives. These
are diagnostics, not frozen acceptance gates, interface tests or evidence
of unforced finite-volume preservation; trace/interpolation and derivative
convergence still need to be separated. Sources and raw rows are retained in
`ambient-continuity-diagnostic.cpp`, `ambient-continuity.csv` and
`ambient-continuity-summary.json` under the evidence root.


The user selected the analytical/semi-empirical path after reviewing both
methods. No MHD runtime is authorized. Production evolution must supply a
complete conservative lateral map and its qualified physical closure; no
ambient-only transverse velocity or arbitrary interpolation is silently
adopted here. The current material kernel supplies mass, full signed frozen
flux, adiabatic entropy and material rates for a supplied compatible map.
Its four tests and the current 1,067-case core campaign are recorded in the
full-CME ledger. The earlier 55-case focused snapshot is 54 PASS/1 ERROR;
it does not include these subsequent four cases.
