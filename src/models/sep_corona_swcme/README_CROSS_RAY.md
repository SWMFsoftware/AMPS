# Conservative analytical cross-ray development

Updated 2026-11-10. This README documents the new shared APIs, their numerical
contracts, independently tested behavior, and remaining production work.
[Main model README](README.md) · [Model-local tests](test/README.md) ·
[Global adapter](../../test/tests/models/cme/sep_corona_swcme_full/README.md).

## Selected scope and current result

The user selected the existing 3-D PFSS/Parker ambient, conservative lateral
material transport, an analytical/semi-empirical regional closure **without an
MHD runtime**, a field-threaded ideal material contact, and startup corotation
with a 1200-second taper to zero angular rate. The full-volume program resumes
at G1/BG3D-4. Existing reduced-front, generic PIC, particle and legacy SWCME
behavior retains its own authority and validation scope.

The implementation provides 29 new cross-ray/material/contact cases in eight
binaries, plus two independent contact/shell controls. Together with the new
coverage-governance case these add 32 global cases relative to the 1035-case
remaining-failure repair baseline. The complete core campaign
now has 1,067 cases: **825 PASS, 1 FAIL, 240 SKIP, 1 ERROR**. All 792 previous
PASS identities remain PASS; CMBGU05 changes FAIL to PASS; all 32 new cases
PASS. The sole FAIL is strict incomplete coverage; the ERROR is the original
180-second production-refinement ceiling. No local PASS closes BG3D-4/G1.

A production lateral map, common moving material geometry, shock/thermal/work
closure, integrated runtime selector, actual ejecta matching and native
full-CME publication are still open. There is no input-file switch enabling a
new cross-ray provider. The existing per-ray event/contact assets are unchanged.
The 1200-second choice is an explicit argument of the new contact API and its
qualification fixture; it was not inserted into the legacy asset.

## Code and public contracts

Public headers live in `include/sep_corona_swcme`; implementations in `src`;
independent tests in `test`. None of the new public contracts requires PIC,
MPI, an AMPS mesh, or an application type. Every vector is Cartesian in the
same Sun-centred inertial HCI frame. Length, time, density, pressure and field
use metres, seconds, kg/m^3, pascals and tesla; angles are radians and angular
weights are steradians. Matrix storage is row-major.

| Header / implementation | API and intended use | Validated scope |
| --- | --- | --- |
| `spherical_shell.h` | `SphericalShellVolumePerSrM3`, `SphericalShellCentroidM`; exact stable shell moments used by the existing radial solver/provider | Broad, compressed and one-ulp shells; independent Gaussian moments |
| `cross_ray_transport.{h,cpp}` | `ToConserved`, `ToPrimitive`, `FastSpeedMPerS`, `HllFlux`, `TransferStationaryFace`; conservation/characteristic reference operators | SI fluxes, orientation/rotation, ALE terms, opposite transfers and analytic advection |
| `spherical_sector.{h,cpp}` | `MeasureSphericalSector`; exact static volumes, centroids and curved-face measures in an explicit orthonormal basis | Octant, rotations, whole shell, thin geometry and rejected invalid frames |
| `cross_ray_mesh.{h,cpp}` | `StaticSphericalMesh::Create`, `SampleEdge`; one immutable shared face/edge topology | Seam/pole/origin ownership, face incidence, curved edge line elements and exact integer boundary-of-boundary cancellation |
| `cross_ray_ct.{h,cpp}` | `StaticMagneticFlux::FromFaceFlux`, `FromEdgePotential`, `Advance`; static constrained-transport reference | Independent Stokes/Faraday, raw divergence and rejected monopoles/stale data |
| `cross_ray_reconstruction.{h,cpp}` | `ReconstructCellField`; constant Cartesian field fit to authoritative CT flux | Rotation/scaling, known affine field refinement, rank rejection and zero field |
| `cross_ray_fluid.{h,cpp}` | `StaticFluidState::Create`, `Advance`; atomic static conservative-fluid/CT reference update | Total-energy recovery with new B, opposite vector transfers, boundary/source ledgers and whole-epoch rejection |
| `cross_ray_contact.{h,cpp}` | `CorotatingPistonContact`, `EvaluateIdealContactResiduals`; selected startup geometry and dimensional threaded-contact diagnostics | Actual four-ray startup, independent time/ray derivatives, Maxwell/electric/work jumps and typed exclusions |
| `material_transport.{h,cpp}` | `Material::PushForwardAdiabatic`; exact local Cauchy/adiabatic transport for a supplied regular map | Signed deformed-face flux, mass/entropy, full gradient/rates, nonlinear transport and rejected folds |

The HLL/static-fluid operators remain independently tested reference tools.
They do not autonomously select a production timestep, flux quadrature, edge
EMF, physical source law or domain boundary. Their existence does not authorize
an MHD time integrator for the selected analytical model.

## Mass, frozen flux and thermodynamics

For one actual compatible material map X(a,t), with reference a in metres:

```text
F = dX/da, J = det(F) > 0
v = partial_t X at fixed a
L = grad_X(v) = Fdot * inverse(F)
rho = rho0/J
B = F*B0/J
p = p0*J^(-gamma)                     [adiabatic subset]
D rho/Dt = -rho*trace(L)
D p/Dt = -gamma*p*trace(L)
D B/Dt = Fdot*B0/J - B*trace(L).
```

A material face transforms as `dS=J*inverse(F)^T*dS0`, giving signed
`B.dS=B0.dS0`; mass satisfies `rho*dV=rho0*dV0`. The implementation uses the
full tensor, including shear and lateral deformation. It does not substitute
a passive radial-field invariant for this relation.

`MapJet` supplies X, v, acceleration, F and Fdot from one prescribed map. The
pointwise kernel checks finite data, positive resolved J and representable
physical output. Global compatibility, inverse uniqueness, reference-field
divergence, event ownership, shock matching and force/work balance require
separate actual-map evidence. Passing the algebra does not select a velocity
closure. The adiabatic subset does not preserve an isothermal expanding ambient
without its explicitly derived heating/work. See
[the material equations and tests](docs/MATERIAL_TRANSPORT.md) and
[the next analytical design investigation](docs/ANALYTICAL_LATERAL_TRANSPORT_DESIGN.md).

For example, a static isotropic deformation F=2I has J=8, rho=rho0/8,
B=B0/4 and p=p0/32 for gamma=5/3. A minimal neutral-header consumer is:

```cpp
#include "sep_corona_swcme/material_transport.h"
#include "sep_coronal_cme/constants.h"
#include <iostream>

int main() {
  namespace XR = SEP::CoronaSwcme::CrossRay;
  namespace M = SEP::CoronaSwcme::Material;
  XR::Primitive reference{2.0, {0,0,0}, 3.0, {1,-2,0.5}};
  M::MapJet map;
  map.positionHciM = {2,0,0};
  map.deformation = {2,0,0, 0,2,0, 0,0,2};
  const XR::Thermodynamics gas{
      5.0/3.0, SEP::CoronalCME::Constants::kVacuumPermeabilityHPerM};
  const auto state = M::PushForwardAdiabatic(reference, map, gas);
  if (!state.ok()) {
    std::cerr << state.status.message << '\n';
    return 1;
  }
  std::cout << state.value.jacobian << ' '
            << state.value.primitive.densityKgM3 << '\n'; // 8 0.25
}
```

These are manufactured SI values illustrating the API. A CME map, reference
inventory and physical parameters still have to be supplied and qualified.

## Exact geometry and immutable ownership

`StaticMeshInput` takes explicit radial edges, mu edges (`mu=sin(latitude)`),
longitude edges, right-handed HCI basis, periodic-longitude choice and positive
cell budget. The factory also requires an owning event SHA-256, finite epoch
and positive generation. It does not infer finite-cell boundaries from legacy
ray centres/weights. Native distribution, moving faces, remapping, piston cuts
and application publication remain unimplemented in this static authority.

On a curved face, distinguish `integral n dA` from scalar `integral dA`.
Their magnitudes generally differ. A volume centroid is an exact first moment;
it need not lie inside a broad nonconvex sector and is not automatically a
valid plasma sampling location. Poles, the origin and the periodic seam use
exact coordinate/topology identifications, not a distance-based merge floor.

Each face has one canonical left-to-right orientation and each curved edge
one canonical node direction. All incident cells/faces refer to that same
object. `StaticSphericalMesh`, `StaticMagneticFlux`, `StaticFluidState` and
`CorotatingPistonContact` have private default constructors: only validated
factories can create them. Failed candidate publication leaves the committed
state unchanged. Fingerprints bind geometry/event and actual state; epoch and
generation remain explicit. CT interval start and generation must match the
committed authority exactly; approximate timestamp equality cannot accept stale
data. There is no mutable shared hash cache or hidden MPI synchronization.

## HLL, canonical characteristics and ALE reference flux

For n directed left to right, vn=v.n, Bn=B.n and Pt=p+B^2/(2*mu0):

```text
F_rho = rho*vn
F_momentum = rho*v*vn + Pt*n - Bn*B/mu0
F_energy = (E+Pt)*vn - (B.v)*Bn/mu0
F_B = vn*B - Bn*v
G_ALE = F - wn*U.
```

wn is the physical face normal speed in the same HCI frame. Substituting
vn-wn only in selected terms would omit energy/Maxwell work. The HLL central
branch uses `(sR*G_L-sL*G_R+sL*sR*(U_R-U_L))/(sR-sL)` and the two supersonic
branches select the appropriate one-sided flux. Signal estimates use full
canonical fast magnetosonic speed, not the legacy reduced radial speed:

```text
a^2=gamma*p/rho, vA^2=B^2/(mu0*rho), vAn^2=Bn^2/(mu0*rho)
vAt^2=abs(B-Bn*n)^2/(mu0*rho)
cf^2=((a^2+vA^2)+sqrt((a^2-vA^2)^2+4*a^2*vAt^2))/2.
```

The implementation evaluates the discriminant with nonnegative acoustic/
transverse terms. Left/right normal B must agree within explicit dot-product
roundoff; no projection enforces agreement. Invalid states, unresolved waves
or nonpositive updated thermal energy reject. The two-wave reference is not
an exact contact/rotational-wave solver or an unconditional positivity theorem.
More equations and provenance are in
[the conservative-reference contract](docs/CONSERVATIVE_CROSS_RAY_TRANSPORT.md).

## Magnetic and total-energy references

CT stores a single signed magnetic flux per shared face. For each interval:

```text
Phi_new = Phi_old - dt * sum_oriented_edges(integral E.dl).
```

One common edge EMF enters every incident face. Integer incidence preserves
cell divergence in exact arithmetic. Raw floating-point divergence and upward-
rounded arithmetic uncertainty are retained separately; no cleaning, field
projection, physical floor or residual subtraction is used. `FromEdgePotential`
requires supplied line integrals of an actual vector potential; it does not
silently project arbitrary measured B onto another field.

The cell reconstruction minimizes
`sum_faces[((S/A).B)-(Phi/A)]^2`, where S is the area vector and A the physical
curved area. It retains normal-fit residuals and a Cholesky-pivot ratio. The
latter is a numerical conditioning diagnostic, not a full spectral condition
number or a physical validity verdict. A nonuniform divergence-free field can
have a nonzero constant-fit residual.

Static-fluid updates take actual time-averaged AREA-INTEGRATED mass, all
Cartesian momentum and total-energy fluxes, explicit volume force/power and
common edge EMF. Internal transfers enter neighbors with opposite signs;
missing-neighbor faces enter a signed outward-boundary ledger. No volume mass
source is supported. Pressure uses the conservatively updated total energy and
the NEW CT field:

```text
E = p/(gamma-1) + rho*abs(v)^2/2 + abs(B)^2/(2*mu0).
```

There is no magnetic-energy reset to hold pressure fixed. Invalid density or
thermal energy rejects the whole epoch. Manufactured externally prescribed
curl-E/source-work tests verify this coupling contract; they are not an ideal
MHD edge-solver or a CME heating-law qualification. Moving-face pair transfer
returns `UnsupportedCapability` until actual moving geometry/GCL is supplied.

## Corotating startup and threaded contact

The original stationary HCI contact has nonzero ambient normal flow on the
four initial supported rays. The user-selected geometry rotates the WHOLE
ellipsoid, including its centre, about HCI z. For tau=(t-t0)/T in [0,1]:

```text
omega = Omega_ambient*(1 - (10*tau^3-15*tau^4+6*tau^5))
theta = Omega_ambient*T*(tau-2.5*tau^4+3*tau^5-tau^6).
```

T=1200 s is explicit. After the transition omega=0 and theta=Omega*T/2,
about 0.0985 degrees for this fixture. No orientation or velocity reset occurs.
Material surface-label velocity includes `omega cross X`. Geometry constrains
physical NORMAL motion; the label velocity is not an independent specification
of later ejecta/sheath plasma tangential velocity. Fixed-HCI ray derivatives
include orientation derivatives as well as axis/centre motion.

The selected continuous-field ideal material interface has zero relative normal
mass flux and continuous v, B and p; density/entropy may jump. The same
continuous branch is retained through Bn=0 without a floating-point branch
switch or division by Bn. The separate closed/draped tangential-discontinuity
branch is not implemented. Diagnostics retain signed leakage, full field/
velocity/pressure jumps, Maxwell traction, contact-frame tangential electric
jump and HCI relative total-energy-flux jump. Actual two-sided ejecta matching
remains BG3D-5 work. See [startup/interface details](docs/CROSS_RAY_CONTACT_STARTUP.md).

## Contact derivative and stable radial-shell repairs

CMBGU05 independently differentiates sampled positions using the equivalent
fourth-order paired form:

```text
D_h X = [8*(X(t+h)-X(t-h))-(X(t+2h)-X(t-2h))]/(12*h).
```

The independent grader computes `rho*(v-D_h X).n`; an API-reported zero leakage
is never its oracle. The original h values 0.1/0.05/0.025 s, physical flux
bound and convergence clause remain. The user explicitly approved replacing
an old expectation that some distributed level fail with the stronger demand
that EVERY level pass unchanged physical limits. Manufactured plane/sphere,
imposed signed slip, convergence and malformed-data controls validate the
actual grader. This passes bounded software validation of the diagnostic map;
its physical force/work discrepancy and cross-stage mismatch remain reported.

The radial solver/provider use exact factored shell volume
`(rR-rL)*(rR^2+rR*rL+rL^2)/3` per steradian and a positive-offset centroid.
This removes cancellation in differences of large powers and avoids placing a
one-ulp shell centre beyond a face. Existing masses, RK4/CFL, viscosity,
source samples and physical equations are unchanged. Invalid/nonrepresentable
geometry returns NaN and is rejected by the caller's existing typed paths.
Independent Gaussian polynomial moments grade the geometry. The older
cancellation-sensitive centroid oracle was corrected independently, restoring
CSWC0625's radial flux error to 4.01794e-16 against its unchanged 2e-13 bound.

## Build, run and retained evidence

See [the model-local test README](test/README.md) for exact quick binaries and
[the global adapter README](../../test/tests/models/cme/sep_corona_swcme_full/README.md)
for the complete 46-case owner. Standalone model-library builds are not native
AMPS builds. Actual native rebuilds retain the required root/idle/build-removal,
configuration-generation and top-level -j16 workflow.

The full core command and identity comparison are recorded in
[`src/test/STATUS.md`](../../test/STATUS.md). Evidence is repository-local under
`test_output/remaining-failures-repair-20261110/`:

- `core-current/amps-core-tests.json` and `.xml`: complete 1,067-case result;
- `core-current-baseline-comparison.json`: all 792 old passes preserved, all
  32 new cases pass, only CMBGU05 changes baseline status;
- `final-kernel-validation.json`, `cross-ray-*-results.json`,
  `contact-latest-profiles.json`, `material-profiles.json`: focused build/
  debug/sanitizer receipts with exact argv and scope;
- `bounded-comparisons.png` and `bounded-comparisons-data.json`: visualization
  of recorded contact flux and analytic-reference convergence errors;
- frozen original binaries, stencil/profiler/coarse-run measurements and all
  first failed/superseded diagnostic receipts remain alongside final results.

All new kernel translation units have debug and ASan/UBSan/LeakSanitizer
coverage; existing linked archives remain uninstrumented production dependencies
where stated by the receipt. Restricted ptrace environments cannot run
LeakSanitizer reliably, so retained all-enabled runs use the explicit local
execution exception. This model-specific evidence does not expand T12's
separate bounded array-sanitizer claim into whole-AMPS coverage.

The original full refinement still times out at 180 s. Its complete original
coarse diagnostic detects two shock zones; a faster algebraically stable
coarse diagnostic also detects two. Neither a larger diagnostic execution
budget nor a one-level zero exit is a refinement PASS. The full material/map,
RH/contact, force/heating/work, validity, assembly and later native gates remain
open. Production source/input selection must not be inferred from private
numerical experiments stored only in `test_output`.


## Private ambient-map investigations

Standalone material-flow diagnostics under the evidence root integrate only
the already prescribed ambient velocity and its variational gradient. Closed
corotation reproduces density/field to roundoff; Parker density/field errors
are of order 1e-12 in the recorded 1200-s probe. The actual open-PFSS 60-s
probe has a step-independent density discrepancy about 4.88e-8. An independent
dipolar flux invariant finds approximately 1.8e-6 pointwise relative density error from
the existing numerical field-line endpoint trace on those points.

A PRIVATE exact-dipolar-trace prototype reduces the density-reference error to
about 1.6e-10 and the analytic endpoint-density error to about 1e-13; its full
CMBGU05 regression passes unchanged gates. Neither that prototype nor an
ambient-flow production API has been imported into source. The generic PFSS
field evaluator is unchanged. Diagnostic receipts are
`ambient-material-flow-receipt.json`, `ambient-material-density-receipt.json`,
`dipolar-trace-prototype-receipt.json` and
`contact-with-dipolar-trace-prototype.json`. The first mislabeled open seed and
first failed diagnostic build are retained separately. These investigations
identify the next numerical prerequisite; they do not select the lateral CME
velocity law or qualify full-volume evolution.
