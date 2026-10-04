# BG3D-4: review and authoritative corrections for Codex

Date: 2026-10-04. This is a separate review/addendum to
`Pasted text(20261004-061828).txt`; it preserves that submitted document. It
does not inspect the current code or verify its reported leakage. The input
requirements are broadly sound after the corrections below. They are a
requirements/validation document, not yet a fully closed sheath algorithm.

## Assessment

Retain the corrected reference-frame explanation, material-contact condition,
independent planar fixture, propagation/expansion distinction, consistent
initial-inventory choices, topology selection, regular map/inverse tests,
shock-side RH matching, conservation, failure integrity and regression gates.
Do not reset inventory at a prescribed front's local fast-mode activation.
Do not require a global MHD solver. The active work remains background-only
G1 in srcSEP3D; preserve baseline particle code and read-only Mars1.

The following corrections take precedence over the corresponding expressions
or requirements in the submitted document.

## 1. Comparisons of normal velocities

The numerical example comparing downstream and contact velocities is planar
or a locally aligned nose example. In a curved layer, shock and contact normals
are generally different. Evaluate each boundary condition at its own physical
point with its own normal, and use a declared geometric correspondence for
cross-layer comparisons. Do not subtract independently projected normal
components as if they used the same basis.

The shock mass-jump identity is correct for the stated convention:

\[
w_1=V_{\mathrm{sh},n}-\mathbf U_1\cdot\mathbf n_{\mathrm{sh}}>0,
\qquad
\rho_1w_1=\rho_2w_2,
\qquad
w_2=V_{\mathrm{sh},n}-\mathbf U_2\cdot\mathbf n_{\mathrm{sh}}.
\]

For the uniform planar constant-speed reference, use the complete admissible
upstream/downstream states, initial coincident piston/shock and the appropriate
uniform-field/EOS assumptions. Disable gravity or other nonuniform forces in
that exact reference. It is not the only imaginable force-free flow, but it is
the only ballistic profile presently proposed for qualification here.

## 2. Non-singular contact-flux tolerances

Do not normalize contact leakage by a local `rho1*w1` without a well-defined
contact-to-shock mapping and a guaranteed nonzero denominator. Admission can
vanish or become very weak, and a contact may exist without a corresponding
currently active shock patch.

At each contact sample, define physical flux

\[
f_c=\rho(\mathbf U-\mathbf V_c)\cdot\mathbf n_c.
\]

Use predeclared dimensional and relative bounds such as

\[
|f_c|\leq \epsilon_{f,\mathrm{abs}}
             +\epsilon_{f,\mathrm{rel}}f_{\mathrm{ref}},
\]

where `epsilon_f_abs` has units kg m^-2 s^-1 and `f_ref` is a finite,
independent case scale with the same units. This is diagnostic scaling, not a
physical density/velocity floor. Record the signed flux and the absolute
crossing measure:

\[
M_{c,\mathrm{abs}}=\int_{t_0}^{t_1}\int_{S_c(t)}|f_c|\,dA\,dt.
\]

A signed integral alone can hide large opposing local leaks. Test both
pointwise residual and accumulated absolute leakage under refinement. Mass
budget tolerances also need an absolute kg term for cases with zero initial
inventory; never divide by zero reference mass. Choose scales/tolerances
before grading and do not change them to pass a failed fixture.

## 3. Initial inventory and the pre-shock diagnostic

All three submitted initializations remain valid when they satisfy their
declared compatibility conditions: zero-inventory coincident surfaces for a
manufactured case, an owning checksummed initial map/state, or reconstructed
prehistory. Zero thickness is not a general consequence of first local shock
admissibility. Initial reference states must cover a finite positive-density
sheath whenever one exists.

The submitted integral of ambient density over contact-swept volume is at
most a **displaced-ambient reference mass**. It is not automatically the mass
retained in a pre-shock compression layer. Lateral diversion, initial motion,
nonuniform wind and incoming/outgoing boundary fluxes matter. Avoid overlapping
swept-volume double counting. A retained pre-shock inventory needs an actual
compatible spatial state or prehistory and a control-volume/material ledger.

Do not impose that approximate diagnostic as a required physical inventory.
If the reference case deliberately omits a pre-shock compression disturbance,
label and bound that approximation; do not present it as undisturbed ambient
in a profile claiming full disturbed-plasma coverage.

## 4. Exact material identities require a defined reference map

The submitted identities `rho=rho0/J` and `B=F*B0/J` are correct only when F
is deformation relative to physical reference-volume coordinates and the
reference deformation is the identity. They are not automatically correct
when labels are admission time and angular coordinates.

For a consistent regular reference map `X0(a)` for an initialized material
volume or appropriately constructed admitted cohort, define

\[
A(t)=\frac{\partial\mathbf X(\mathbf a,t)}{\partial\mathbf a},
\quad A_0=\frac{\partial\mathbf X_0(\mathbf a)}{\partial\mathbf a},
\quad F_{\mathrm{rel}}=A(t)A_0^{-1},
\quad J_{\mathrm{rel}}=\det F_{\mathrm{rel}}.
\]

Then, for smooth ideal evolution with a fixed-composition, constant-gamma
adiabatic closure,

\[
\rho=\frac{\rho_0}{J_{\mathrm{rel}}},\qquad
\mathbf B=\frac{F_{\mathrm{rel}}\mathbf B_0}{J_{\mathrm{rel}}},\qquad
p=p_0J_{\mathrm{rel}}^{-\gamma}.
\]

Equivalently, density in arbitrary labels obeys
`rho(t)*det A(t)=rho0*det A0`. Record the physical mass and magnetic flux of
each admitted reference volume. A boundary point by itself has no invertible
three-dimensional reference deformation: construct the admitted volume and
its time/angle metric correctly before applying these identities.

Shock-born parcel/cohort reference states are their actual RH downstream
states, not upstream values. Initialize a spatially consistent magnetic-flux
reference, not unrelated B vectors at independent points. The Cauchy identity
preserves divergence constraints only when the initial field/map has the
required compatibility. Entropy is conserved after admission only for the
declared ideal adiabatic profile with no further shocks/heating/dissipation.
Reject other requested thermal profiles unless explicitly implemented.

Use regular coordinate charts, including at angular poles, and a documented
dimensionless relative-Jacobian/conditioning test. Positive local J does not
prove global injectivity. Verify no overlaps, duplicates or holes and a unique
supported inverse. A deliberately zero-volume initial sheath requires its own
limiting initialization; do not invert it or hide it with a density/J floor.

## 5. Curved-layer conservation: correct the mass measure

The unweighted normal column `integral rho d ell` is not generally the mass per
unit reference area of a finite curved sheath. Define a regular layer map,
the reference-surface velocity/parameterization and its physical volume metric
explicitly. For normal coordinates, write

\[
dV=h(\mathbf q,\ell,t)\,d\ell\,dA_{\mathrm{ref}},
\qquad
\Sigma(\mathbf q,t)=\int \rho h\,d\ell,
\qquad
M(t)=\int_{S_{\mathrm{ref}}(t)}\Sigma\,dA_{\mathrm{ref}}.
\]

Here h is the geometrical volume factor, not an empirical density adjustment.
For concentric spherical layers using the inner radius Rc as reference,
`h=(1+ell/Rc)^2`. This provides a simple exact negative control against using
an unweighted column. For a general non-normal map use its full metric/Jacobian.
Normal coordinates must remain regular and cover the intended layer; do not
assume every finite triaxial sheath admits one global normal chart.

Let W be the chosen reference-surface parameterization velocity, `D_W` the
derivative at fixed reference labels, and q_s the appropriately mapped,
depth-integrated mass flux tangent to the reference surface. A conservative
surface balance has the form

\[
D_W\Sigma+\Sigma\nabla_s\cdot\mathbf W+
\nabla_s\cdot\mathbf q_s=s_{\mathrm{sh}}+s_{\mathrm{other}},
\]

where the shock admission source is obtained by a defined area correspondence,
`s_sh=rho1*w1*(dA_sh/dA_ref)`, only where that correspondence is regular.
`s_other` contains only explicitly declared physical transfers. Flank/domain
boundary fluxes use velocity relative to the moving edge, with a consistent
outward conormal. Do not count both a boundary flux and a duplicate source term.
For a closed layer, do not invent a flank exit to enforce steady inventory.

The shorthand `q_s=Sigma*u_bar_rel` is acceptable only when `u_bar_rel` is
defined from that metric-consistent mapped mass flux. It is not automatically
the ordinary arithmetic average of tangential velocity across the layer.

For an irrotational diagnostic closure, a weighted surface problem is

\[
\nabla_s\cdot(\Sigma\nabla_s\psi)=
s_{\mathrm{sh}}+s_{\mathrm{other}}-D_W\Sigma-
\Sigma\nabla_s\cdot\mathbf W.
\]

Specify the gauge, edge boundary conditions, integral compatibility, and a
physical/evolution closure for Sigma. At zero Sigma the operator degenerates;
use a justified limiting/startup formulation or explicitly unsupported state,
not a positive numerical mass floor. Any numerical regularization must vanish
appropriately under refinement, preserve the physical budget and be documented.

If instead using the leading-order thin-layer approximation h=1, declare
thickness/curvature ordering and measure the omitted-geometry error. Do not
claim exact full-volume mass closure using the approximate column.

## 6. A velocity-closure menu is not a closed physical algorithm

Before marking a reduced model implemented, document its complete state,
equations, driving/boundary conditions, unknowns, thermal/magnetic closure,
initialization, applicability and numerical solution. Explain which observed
surfaces are inputs and which quantities are solved/constrained. Do not claim
standoff is obtained from mass alone. Use candidate compatibility checks if
shock/contact histories are prescribed together.

Do not introduce an MHD runtime. A kinematic construction remains an
intermediate diagnostic until all required physical residuals, interfaces,
coverage and validation gates pass within frozen independent bounds. Adding
labels, mass closure and a regular map does not prove momentum or energy.

For a smooth isotropic-pressure ideal background, report at least

\[
R_\rho=\partial_t\rho+\nabla\cdot(\rho\mathbf U),
\quad
\mathbf R_B=\partial_t\mathbf B-\nabla\times(\mathbf U\times\mathbf B),
\quad
\nabla\cdot\mathbf B,
\]

\[
\mathbf R_{\mathrm{mom}}=\rho D\mathbf U/Dt+\nabla p-
\frac{\nabla\times\mathbf B}{\mu_0}\times\mathbf B-\rho\mathbf g.
\]

For the ideal constant-gamma adiabatic branch, also check
`Dp/Dt+gamma*p*div U=0`, total-energy/control-volume balances and interface
traction/work. If additional sources or stresses are permitted by the chosen
model, include them explicitly in the equations and evidence. Physical model
residuals and numerical errors have separate scales/budgets. Evaluate smooth
residuals one-sided near interfaces and validate jumps through the proper
interface conservation conditions; do not differentiate through a shock as
if it were smooth.

## 7. Topology and contact data

The submitted ideal isotropic interface conditions are broadly correct:

- A selected tangential discontinuity has zero relative normal mass flux,
  Bn=0 on both sides, and balanced thermal-plus-magnetic pressure; tangential
  velocity/field can jump.
- A selected transverse-field ideal contact has no relative normal mass flux
  and continuous U, B and p, with density/entropy potentially discontinuous.

Do not equate every flux-rope boundary with a closed tangential surface.
Sun-connected geometry and footpoint/cut boundaries need their own coverage
and flux contracts. Mixed or vanishing-normal-field cases need typed handling,
not a raw `Bn!=0` floating-point switch. The background topology selector is
explicit and validated.

BG3D-4 publishes sheath-side contact position, normal, normal speed, plasma
state, mass flux, Bn and magnetic/thermal stress. For isotropic ideal stress,
retain the traction derived from

\[
T=(p+B^2/(2\mu_0))I-\mathbf B\mathbf B/\mu_0,
\]

not only a scalar pressure where a general magnetic interface is requested.
Use compatible manufactured ejecta-side data for BG3D-4 boundary tests. Full
two-sided matching to the actual ejecta remains a BG3D-5/G1 gate; sheath-side
output alone cannot certify it. Contact motion depends only on normal speed;
tangential parameterization choices must not change physical leakage.

## 8. Reconnection estimate is optional sensitivity, not a required gate

The formula `delta_rec ~ M_rec*v_A*delta_t` is a dimensional local-inflow
estimate, not a demonstrated global ejecta-erosion depth or mass-transfer law.
The submitted arithmetic is approximately correct for its assumed parameters,
but M_rec=0.01–0.1 is not a universal applicable interval for this fixture.

Keep impermeability as an explicit selected approximation. If reporting an
erosion sensitivity, define local reconnecting field, density, inflow length,
duration, affected surface, geometry and uncertainties. Label hypothetical
parameter choices. Do not add a reconnection-physics prerequisite or compulsory
configuration key to BG3D-4 based on this order-of-magnitude estimate.
Any future reconnection/transfer closure needs its own topology, induction,
flux and mass/energy accounting and separate authorization/scope.

## 9. Required implementation and validation instructions

1. Reproduce the actual leakage and contact/flow/map mismatch using the real
   current source. Preserve inputs and receipts. This review is not evidence
   that the reported failure has been reproduced.
2. Retain the exact planar constant-speed fixture and add an independent
   expanding-curved reference with a complete stated velocity/state solution.
   Do not generate expected values with the production map itself.
3. Choose and document a fully defined regional closure in the analytical
   model's scope. Enforce its material contact pointwise; do not tune fraction
   or tolerances, suppress crossings, or delete physically retained parcels.
4. Implement a compatible initialized material volume or declared limiting
   zero-volume fixture. Preserve inventory across local shock activity changes.
5. Verify shock-side RH states, both contact-flux measures, all admissible
   boundary fluxes and initial-plus-admission mass. Verify exact curved volume
   measures or explicitly qualified thin-layer truncation.
6. Verify reference-normalized material/field/thermal identities, independent
   smooth residuals, regular/global map coverage and inverse, and no duplicate
   material or holes. Refine admission time, surface/material resolution and
   timestep separately. Check zero-inventory and vanishing-admission limits.
7. Freeze physical discrepancy and numerical error thresholds before grading;
   do not mark a large-momentum-residual diagnostic as qualified full plasma.
8. Preserve immutable owning candidates: a rejected preparation leaves the
   committed epoch and all state identities unchanged. Distinguish this shared
   preparation gate from later actual native write/halo rollback/fail-stop.
9. Run appropriate ASan/UBSan and architecture checks and maintained coronal
   regressions; report exact commands/PASS/FAIL/SKIP/ERROR and environment
   blockers. Keep required skipped qualifications open; preserve baseline
   tests and report any pre-existing failures separately.
10. Add detailed equations/units/HCI/sign/coordinate/reference-state/coverage/
    approximation/initialization comments, input documentation and READMEs.
    Checksums and event fingerprints cover all added physical/numerical assets.
    Freeze actual schema names through the current parser rather than assuming
    the submitted example keys already exist. Reject unsupported requests.
11. Update CODEX_CME_PLAN.md with the current BG3D-4 implementation and
    validation status. Keep baseline particle functionality unchanged, Mars1
    read-only and generic PIC/legacy SWCME functional. Follow the authorized
    clean native rebuild workflow in Mars2. Do not commit or push.

## Ready-to-paste instruction

```text
Use BG3D4_CODEX_REVIEW_AND_CORRECTIONS.md as the authoritative addendum to the
submitted revised BG3D-4 requirements. Its corrections take precedence.

Reproduce the current leakage, then implement a complete documented analytical
or reduced regional sheath construction with material contact compatibility
and valid initial inventory. Correct curved-volume metrics, reference-map
identities, zero-admission tolerances and degenerate startup handling. Do not
add numerical mass/density floors, a reconnection prerequisite or an MHD solver.

Add the independent physics, interface, inverse/coverage, residual and
convergence tests specified here. A velocity interpolation or surface Poisson
problem alone is an intermediate diagnostic, not qualification of BG3D-4/G1.

Preserve baseline particle code and read-only Mars1. Add detailed comments and
READMEs, update CODEX_CME_PLAN.md with actual evidence, and continue the active
background task under the existing milestone authorization. No commit/push.
```

## Primary references and their limited role

- Dziuk and Elliott, *L2-estimates for the evolving surface finite element
  method*, transport theorem and surface conservation formulation:
  <https://wrap.warwick.ac.uk/id/eprint/54318/1/WRAP_Elliott_0672442-ma-020513-dziell13mcomp.pdf>.
  Supports the moving-area and surface-velocity terms; it does not specify
  the CME closure or the finite-thickness reduction used here.
- Wang and Xin, *Existence of Multi-dimensional Contact Discontinuities for
  the Ideal Compressible Magnetohydrodynamics*:
  <https://arxiv.org/abs/2112.08580>.
  Supports the ideal transverse-field contact and Lagrangian/Cauchy formulation;
  it does not establish the validity of the proposed CME kinematic construction.
- Siscoe and Odstrcil (2008), *Ways in which ICME sheaths differ from
  magnetosheaths*, <https://doi.org/10.1029/2008JA013142>.
  Supports distinguishing propagation and expansion; no universal steady
  bow-shock assumption follows.
- Lulic et al., *Formation of Coronal Shock Waves*,
  <https://arxiv.org/abs/1303.2786>.
  Supports distinguishing wave steepening/formation from an imposed local
  front-admissibility classification.
