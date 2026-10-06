# BG3D-4 replacement production closure: per-ray Lagrangian piston

Status: **selected production closure under qualification; contact/ray
identity, planar and spherical gas/MHD references, projected ambient sources,
moving Parker boundaries, transactional ambient append and contact-driven
multi-ray tubes are implemented and independently tested within their stated
subsets; convergence, validity diagnostics and 3-D assembly remain
unimplemented**.  This
closure is intended to replace
`rh-relaxing-material-map-v1` for production BG3D-4.  The old map remains a
Level A+ diagnostic and cannot satisfy the BG3D-4 production gate.

This document records the equations, units, boundary and initial data,
algorithm, acceptance tests and known limitations before implementation.  A
passing test of the existing diagnostic map is not evidence for this closure.

The first bounded implementation is `piston_contact.{h,cpp}` with
`examples/bg3d4_piston/contact.asset`.  It implements the separate immutable
contact identity, global quintic startup, fixed-HCI ellipsoid, analytic
outermost-ray intersection and time derivatives, typed support, smooth
self-similar DBM handoff, 1-AU coverage and ambient startup preflight.
`CSWC0620` independently tests those claims and the strictly parsed checksummed
ray asset.  It does **not** implement or validate production ambient sources,
the computed 3-D shock, transverse-pressure diagnostic, 3-D assembly, or
production srcSEP3D publication.

The subsequent bounded increments are `piston_solver.{h,cpp}` with
`CSWC0621`--`CSWC0624`.  They exercise the actual moving-piston, material-cell,
shock-zone and energy-ledger interfaces using planar gas, planar transverse
MHD and spherical gas columns.  The selected controls
are staggered immutable masses, exact planar area, RK4 at CFL 0.25 and
quadratic-only VNR (`C2=2`, `C1=0`).  Quadratic-only viscosity is intentional:
an always-active linear term produced leading-order damping of the required
sub-fast wave, while a hard activation sensor introduced a nonconvergent shock
threshold.  RK4 supplies linear acoustic stability without that dissipation.

For formation convergence, `epsilon_sh=2.6/N` tends to zero with mass-cell
size.  A fixed threshold is valid for a finite-resolution production shock but
cannot converge to the zero-amplitude analytical birth event.  Failed fixed and
`N^(-2/3)` detector sequences are retained in the plan.  `CSWC0624` separately
qualifies exact spherical-sector area/volume geometry against an independently
integrated Taylor strong-shock similarity solution.  It initializes the
manufactured solution at finite time through `CreateInitialized`; that
interface is for exact fixtures and restarts and is not permission to invent
production prehistory.  `CSWC0624` additionally qualifies dynamic spherical
magnetic energy/work convergence and both signed frozen-flux components.
`CSWC0625` qualifies branch-safe frozen ambient-source projection, a
nonuniform closed-corona well balance, independent nonzero heating/induction
ledgers and moving Parker material boundaries.  `CSWC0626` qualifies passive
radial flux and transactional outer ambient insertion with complete incoming
energy.  `PistonSheathModel` then drives every supported tube from the one
analytical contact, appends an undisturbed buffer, detects a computed shock,
and exposes typed radial region queries.  `CSWC0627` qualifies that bounded
interface/transaction path at 1800 s.  Its sampled shock is quasi-parallel, so
the required quasi-perpendicular production-RH subgate is explicitly not
applicable and remains open rather than being inferred from the smoke PASS.

## 1. Pre-implementation consistency review

The supplied piston construction closes the principal defect in the old map:
the ejecta contact drives a material plasma volume, while the sheath and shock
are outcomes rather than independently prescribed surfaces.  The following
points must be resolved in the implementation contract; none may be hidden by
a tolerance or by assuming a spherical surface.

1. **Selected radial reduction and its declared discrepancy.**  The
   Lagrangian map `x(m,t)=r(m,t) q_hat` with fixed `q_hat` and
   `partial_t r=u` has material velocity exactly `U=u q_hat`.  Publishing a
   nonzero carried `U_t,a` as part of that material velocity would move plasma
   between rays without a lateral mass, momentum, induction or energy flux.
   At the contact, the inconsistency is visible in the exact normal condition

   ```text
   (U - V_c).n_c = (u - dR_c/dt) mu_c + U_t,a.n_c = 0.
   ```

   The user selected radial-only Level B *dynamics*: the conservation system
   evolves `u q_hat`, while the provider adds `U_t,a` as a prescribed,
   unevolved component.  The simultaneous fixed-mass boundary equations
   `r(0,t)=R_c(t)` and `partial_t r(0,t)=u(0,t)` require
   `u=dR_c/dt`; they are not altered to cancel transverse slip.  Consequently
   the reduced radial contact is impermeable, but the published full vector is
   not claimed to be a three-dimensional material map.  Report both the full
   normal slip/mass flux and the delegated model-error diagnostic
   `epsilon_t=|grad_t P|/|partial_r P|`.  Moving rays and cross-ray fluxes are
   outside BG3D-4 and are reconsidered only if these diagnostics invalidate the
   supported production rays.
2. **Normal magnetic field and shock obliquity use physical normals.**  On a
   nonspherical contact, `B_r/|B|` is a radial-field fraction, not
   `|B.n_c|/|B|`.  Likewise a fitted shock surface generally has
   `n_sh != q_hat`; its obliquity is
   `acos(|B_1.n_sh|/|B_1|)`.  Both the radial proxy and the physical-normal
   diagnostic may be reported, but they must have different names.
3. **The reduced PDE and the canonical MHD characteristic are distinct.**  In
   the stated model `B_r` is passive and the radial momentum equation contains
   only transverse magnetic pressure and curvature stress.  Its compression
   speed is therefore

   ```text
   c_red^2 = (gamma_ad p + |B_t|^2/mu0)/rho.
   ```

   The canonical ideal-MHD fast speed for propagation along `q_hat` is

   ```text
   a^2 = gamma_ad p/rho,
   v_A^2 = |B|^2/(mu0 rho),
   v_Ar^2 = B_r^2/(mu0 rho),
   c_f^2 = 0.5[(a^2+v_A^2)
              + sqrt((a^2+v_A^2)^2-4 a^2 v_Ar^2)].
   ```

   Use `c_red` for the reduced solver's CFL and viscosity unless a different
   characteristic derivation is supplied.  Use `c_f` for physical shock
   diagnostics.  Their disagreement is part of the oblique-ray validity
   error, not numerical truncation error.
4. **A ray needs an absolute area measure.**  Define
   `A_q(r)=DeltaOmega_q r^2`, where `DeltaOmega_q` is the ray's solid-angle
   quadrature weight.  Then `dm=rho A_q dr` has units kg.  A proportionality
   alone is insufficient for mass, energy, piston force or refinement
   comparisons.  Ray weights, neighbors and unsupported solid angle belong to
   the checksummed ray asset.
5. **The handoff must transfer the contact driver.**  The current event
   history calls its ellipsoid a front/shock.  Level B instead requires a
   reconstructed ejecta/contact body through the entire corona-to-SWCME
   transition and outer propagation.  The existing prescribed shock becomes
   a Level A diagnostic and Level B validation target.  It must not continue
   to drive the sheath simultaneously with the piston-computed shock.
6. **Energy accounting includes every maintained source.**  `W_amb` must
   include work by `rho f_amb`, heat `h_amb`, and magnetic-energy change caused
   by `S_b`.  Artificial viscosity transfers resolved kinetic energy to
   internal energy; it is not an unledgered heat source.  Appending ambient
   cells also carries mass, momentum, internal and magnetic energy through the
   moving outer boundary.
7. **Divergence cleaning is a physical-state modification.**  A projection
   requires declared boundary conditions and must report changes in magnetic
   energy, contact/shock normal flux and per-ray conservation.  It cannot be
   used merely to force `delta_B` below its limit.

These are implementation prerequisites, not reasons to retain the empirical
relaxation map.  They are recorded now so the production code and independent
tests cannot silently adopt inconsistent spherical or radial assumptions.

## 2. Conventions and supported geometry

- Public units are SI.  Coordinates are Sun-centred heliocentric inertial
  Cartesian (HCI), with the Sun at the origin.
- A fixed ray is a unit direction `q_hat` in HCI and has position
  `x=r q_hat`.  Rays do not rotate with the Sun or CME.
- The reconstructed ejecta surface is the one contact authority
  `g_c(x,t)=0`, with `g_c<0` in the ejecta.  Its normal
  `n_c=grad(g_c)/|grad(g_c)|` points from ejecta into sheath.
- `gamma_ad=5/3` is the production ideal-gas response unless a differently
  fingerprinted closure is selected.
- Ambient quantities carry subscript `a` and are sampled from the one ambient
  authority at the parcel's current HCI position and physical epoch.

For a ray, `R_c(q_hat,t)` is the unique root of
`g_c(R_c q_hat,t)=0`.  The ray is supported only when:

```text
the contact is star-shaped on the ray (one root),
mu_c = q_hat.n_c >= mu_min > 0,
R_c >= R_in,
and the complete initial/buffer interval has ambient coverage.
```

Grazing, multiply intersecting and unsupported flank rays return the typed
status `unsupported-flank`.  Their physical contact area and solid angle are
budgeted; they are not filled with quiet ambient and are never SEP-source
eligible.

The analytical contact speed follows from the level set:

```text
dR_c/dt = -partial_t(g_c)/(q_hat.grad(g_c)),
V_c,n = (dR_c/dt) mu_c.
```

This derivative must come from the differentiable contact history, including
handoff weight-derivative terms.  Subtracting two `O(10^9 m)` positions is not
the production velocity algorithm.

## 3. Per-ray state and material coordinate

Ray `q` owns the physical area `A_q(r)=DeltaOmega_q r^2` and Lagrangian mass
coordinate

```text
dm = rho A_q(r) dr.
```

Cell/node unknowns are

```text
r(m,t), u(m,t), v(m,t)=1/rho, e(m,t),
b_t(m,t)=B_t/(rho r),
```

where `b_t` has two HCI components in a deterministic ray-tangent basis.  The
basis and its behavior at poles are part of the ray asset.  Derived state is

```text
p = (gamma_ad-1) rho e,
B_t = rho r b_t,
B_r = B_r,a(r0) r0^2/r^2,
P = p + |B_t|^2/(2 mu0),
T = p m_bar/(rho k_B).
```

`B_r r^2` is passive flux per solid angle.  `B_t/(rho r)` is the invariant of
the source-free radial ideal-induction reduction.  The evolved material
velocity is `u q_hat`; the published reduced state is
`U=u q_hat+U_t,a`.  `U_t,a` does not enter the mass, momentum, energy or
induction evolution.  Its omitted cross-ray transport and its mismatch with
the radial material map are physical model errors measured with `epsilon_t`
and the full-vector contact-slip diagnostic, not numerical errors.

## 4. Governing equations

At fixed material mass `m`, the quasi-one-dimensional spherical equations are

```text
partial_t r = u,
partial_t v = partial_m(A_q u),

partial_t u = -A_q partial_m(P+Q)
              - |B_t|^2/(mu0 rho r)
              - G M_sun/r^2 + f_amb(r),

partial_t e = -(p+Q) partial_t v + h_amb(r)/rho,

partial_t b_t = S_b(r).
```

The first two equations are equivalent to spherical mass conservation because
`partial_m r=v/A_q`.  The second magnetic term is the declared curvature/hoop
stress of the transverse field.  This reduction is most faithful on
quasi-perpendicular rays; it is not a general oblique one-dimensional MHD
system because tangential momentum is absent.

The production discretization is a conservative staggered Lagrangian
finite-volume method or a Lagrangian Godunov method.  Cell masses are immutable
after insertion.  Node positions and velocities advance with a second-order
or higher time integrator; cell volumes are computed from exact spherical
shell volumes over `DeltaOmega_q`, not `A Delta r` at an arbitrary face.
Time-step limits include the reduced fast signal speed, compression, gravity,
source variation and moving-piston boundary.  Negative volume, density,
internal energy, nonfinite state or node crossing rejects the private candidate
epoch; no density, pressure or Jacobian floor is allowed.

For a von Neumann--Richtmyer option, compression-only viscosity is

```text
Q = C2 rho (Delta u)^2 + C1 rho c_red |Delta u|,  Delta u < 0,
Q = 0,                                             Delta u >= 0.
```

`Delta u` is the outward-face velocity difference under a documented sign
convention.  `C1`, `C2`, shock-zone width and their refinement behavior are
inputs/evidence, not adjustable post-run fitting parameters.

## 5. Ambient-maintaining sources

The only permitted volume accelerations/heating are gravity, numerical shock
dissipation and time-independent ambient-maintaining sources.  Along each ray,
define

```text
f_amb = u_a du_a/dr
      + (1/rho_a) d[p_a+|B_t,a|^2/(2 mu0)]/dr
      + |B_t,a|^2/(mu0 rho_a r)
      + G M_sun/r^2,

h_amb = rho_a u_a [de_a/dr + p_a d(1/rho_a)/dr],

S_b = u_a d[B_t,a/(rho_a r)]/dr.
```

The derivatives use the ambient provider's same-branch derivatives or an
independently converged branch-preserving stencil.  These sources make the
projected undisturbed ambient an exact steady solution of the reduced
equations.  They represent coronal heating/momentum deposition implicit in the
ambient authority and the residual introduced by projecting a 3-D wind onto
independent radial tubes.  Applying them unchanged to shocked plasma is an
explicit approximation.

Implemented representation: `PistonAmbientProjection` fixes the ambient epoch
and a deterministic right-handed HCI tangent basis once per ray.  It retains
both signed transverse components and uses a fourth-order central radial
stencil, falling back to a fourth-order one-sided stencil when necessary to
remain on one region/sector branch.  `PistonAmbientSourceTable` samples that
authority once and linearly interpolates only between samples with identical
region and sector; an interval containing a branch change is evaluated
directly and one-sided instead of being smoothed.  The solver consumes the
well-conditioned net acceleration `-GM_sun/r^2+f_amb`, but the projection
retains the two terms separately for the reported force ratio.  RK4 stages
advance both `B_t/(rho r)` components and accumulate separate body-force,
volume-heating and magnetic-source work.  Table-spacing convergence remains a
production acceptance requirement.

Each ray reports `max |f_amb|/(G M_sun/r^2)`, `max |h_amb|`, magnetic source
work, and the integrated work/heat entering the energy ledger.  Sources vanish
to numerical tolerance when the selected ambient already solves this reduced
steady problem.

Artificial viscosity must be localized:

```text
epsilon_Q = integral_outside_shock |Q partial_t(v)| dm
            / integral |p partial_t(v)| dm
          <= epsilon_Q,max.
```

Expected-zero denominators use predeclared dimensional absolute tolerances;
they are never replaced with a numerical floor.

## 6. Initial and boundary data

At `t0`, the contact apex satisfies
`R_launch >= R_in + delta_launch`.  Each supported ray contains a finite
ambient material column from `R_c(t0)` to `R_c(t0)+L_buf`.  There is no initial
*disturbed sheath*, but there is finite ambient mass, magnetic flux and energy.
This distinction replaces the old zero-volume production startup.

Compatibility conditions are

```text
dR_c/dt(t0) = u_a(R_c(t0))              [conditional radial form],
|d2R_c/dt2(t0)| < infinity,
contact history is C2 or smoother,
cell state at t0 equals the ambient state and satisfies its EOS.
```

The radial compatibility equation is exact only for the radial-only Level B
choice.  A model retaining transverse material velocity requires the expanded
moving-ray/lateral-flux construction in Section 1.  Finite subsequent
acceleration is a physical piston driver and may launch a compression wave; it
is not an initial discontinuity.

The inner face is the impermeable contact.  Its prescribed geometric motion
sets the boundary node and supplies piston work.  Ghost/cell reconstruction
must enforce the selected normal-velocity condition without independently
over-prescribing pressure.

The outer boundary retains at least `N_buf` ambient cells ahead of

```text
r_d(t) = max{r: |P-P_a(r)|/P_a(r) > epsilon_d}.
```

New cells are appended with the current ambient state and checksummed ray
weight.  Their incoming mass, momentum and total energy are included in the
boundary-flux ledger.  Failure to maintain the buffer rejects the candidate;
it does not extrapolate the last cell.

## 7. Magnetic topology and validity

The transverse field evolves through `b_t`; the ray-parallel field obeys
`B_r r^2=constant` and exerts no force in this reduced radial model.  No field
component is floored.

The impermeable piston does not enforce a tangential discontinuity with
`B.n_c=0`: normal-field draping requires the omitted tangential dynamics.  The
implementation therefore reports both

```text
beta_r = |B_r|/|B|,
beta_n = |B.n_c|/|B|,
```

and never labels `beta_r` as the physical contact-normal fraction.  A request
requiring a closed or perfectly draped contact field is unsupported.

The model is most applicable when the physical shock is quasi-perpendicular,
the contact is non-grazing, and tangential pressure gradients are small.
Large obliquity error, `epsilon_t`, `beta_n`, or unsupported solid-angle
coverage removes the associated capability rather than silently returning a
radial approximation.

## 8. Shock formation and diagnostics

A shock zone is a contiguous compressive set with `Q/P > epsilon_sh`.  Select
the outermost zone causally connected to the piston.  Its reported radius is
the `Q`-weighted centroid.  Upstream and downstream states come from the first
settled cells outside the numerical shock zone, with the sampling offset and
settling criterion refined with resolution.

For a resolved density jump,

```text
V_sh = (rho2 u2-rho1 u1)/(rho2-rho1).
```

When `rho2-rho1` is below a frozen dimensional/relative jump threshold, this
formula is ill-conditioned and the feature is a linear compression/no-shock
state; it is not assigned an arbitrary speed or compression floor.

Fit `R_sh(q,t)` over neighboring supported shocked rays to obtain the physical
normal `n_sh` and area metric.  Diagnostics are

```text
X = rho2/rho1,
theta_Bn = acos(|B1.n_sh|/|B1|),
M_f,n = (V_sh,n-U1.n_sh)/c_f,1(n_sh).
```

The radial jump quantities may also be reported as reduced-model diagnostics.
Canonical oblique RH comparisons are graded only inside the preregistered
quasi-perpendicular validity range.

## 9. Regions and three-dimensional assembly

Along a supported ray:

```text
ejecta:              r < R_c,
sheath:              R_c <= r <= R_sh, when a shock exists,
compression region:  R_c <= r <= r_d,  when no shock exists,
ambient:              r > max(R_sh,r_d).
```

BG3D-4 publishes the sheath/compression side only.  `r<R_c` remains a typed
BG3D-5 request and may not be filled by sheath or ambient state.

Neighboring rays interpolate only within the same physical region, using
normalized coordinates such as

```text
xi_s = (r-R_c)/(R_sh-R_c),
xi_a = (r-R_sh)/(r_max-R_sh).
```

No stencil crosses a contact or shock.  Region identity, ray identity,
interpolation provenance and one-sided derivative status accompany every
sample.  The assembled field reports

```text
delta_B = integral |div B| dV / integral (|B|/ell) dV.
```

An optional cleaning projection is a separate fingerprinted algorithm with
its own boundary conditions and magnetic-energy/interface-flux ledger.

The sheath-side contact publication for BG3D-5 contains `R_c`, `n_c`,
`V_c,n`, `rho`, `U`, `p`, `T`, `B`, `P`, `beta_r`, `beta_n`, and

```text
F_c = surface_integral P n_c dA,
W_c = time_integral surface_integral P V_c,n dA dt,
```

with Maxwell traction included separately when the physical-normal field is
nonzero.  Level B source publication uses only the computed shock surface and
its computed upstream/RH diagnostics.  The observed/Level-A shock history is
comparison evidence, not a second source surface.

## 10. Configuration and identity

The closure/contact/ray identity and finite-volume numerical controls are now
accepted by the strict parser.  They are stored in `PistonNumericsInput` and
join the event fingerprint; callers cannot replace them programmatically or
derive them from particle/source settings.  The example uses:

```text
regions.sheath_model=per-ray-lagrangian-piston-v1
assets.piston_rays_file=examples/bg3d4_piston/rays.csv
assets.piston_rays_sha256=<SHA-256>
assets.piston_contact_file=examples/bg3d4_piston/contact.asset
assets.piston_contact_sha256=<SHA-256>

# Finite ambient inventory and source/trajectory resolution.
bg3d4.initial_cells=80
bg3d4.initial_buffer_m=4000000000
bg3d4.source_table_points=257
bg3d4.trajectory_maximum_step_s=1000

# Maintain an undisturbed material buffer beyond r_d.
bg3d4.minimum_buffer_cells=24
bg3d4.append_cells=64
bg3d4.disturbance_threshold=0.0001

# Reduced radial RK4/VNR shock capture.
bg3d4.cfl=0.2
bg3d4.artificial_viscosity=vnr
bg3d4.quadratic_viscosity=2
bg3d4.linear_viscosity=0
bg3d4.shock_detection_threshold=0.01
bg3d4.well_balanced_sources=on
```

The contact, rays, solver/discretization and source policy all participate in
the fingerprint.  Unknown or missing keys, unsupported viscosity/source
selectors and out-of-range controls are rejected before a model exists.
Changing a refinement level therefore creates a distinct immutable event
identity.  Physical tolerances for `epsilon_Q`, RH, `epsilon_t` and `delta_B`
still require strict input keys before production qualification.

## 11. Independent acceptance program

Expected values must be produced outside the production solver path.

The currently implemented shared regression is run from the AMPS root with

```text
make -C src/models/sep_corona_swcme -j16 test
```

As of the 2026-10-04 contact-driven increment it reports executable PASS 14
(`CMBGU01`--`05`, `CSWC0620`--`CSWC0627`, and `ARCHCSWC01`), FAIL 0, SKIP 0,
ERROR 0.  The production quasi-perpendicular RH applicability subgate is
separately recorded as NOT-APPLICABLE/SKIP 1 and therefore keeps BG3D-4 open.
`CSWC0620` qualifies the independent contact plus ray-asset identity;
`CSWC0621` the exact planar gas piston and energy ledger; `CSWC0622` analytical
shock formation plus the linear sub-fast response; `CSWC0623` cold and
finite-beta perpendicular MHD; `CSWC0624` spherical Taylor geometry plus
curved-MHD flux/energy; and `CSWC0625` projected ambient well balance and
source work and Parker transport.  `CSWC0626` qualifies appended ambient
transport, and `CSWC0627` the bounded contact-driven multi-ray/query path.  Its
`CMBGU05` output deliberately continues to report the old
contact mismatch, nonconverged independent leakage and persistent force/energy
residual; that diagnostic PASS is not a production-closure PASS.

No native AMPS command applies at this specification stage.  Once BG3D-7 is
reached, native srcSEP3D evidence must use the repository's mandatory clean
root-build/configuration/production-hook workflow and separate handoff-smoke
and low-corona-to-1-AU decks.  Documentation or a source-only shared test
cannot substitute for those later 1-/4-rank zero-particle receipts.

1. **Planar gas piston:** for constant piston speed `u_p` into gas at rest,

   ```text
   V_sh = ((gamma+1)/4) u_p
          + sqrt([((gamma+1)/4)u_p]^2+c0^2),
   X = V_sh/(V_sh-u_p).
   ```
2. **Cold perpendicular MHD piston:** compare with the `gamma=2` gas analogue,
   `V_sh=3u_p/4+sqrt((3u_p/4)^2+v_A^2)`.  Finite-beta cases use the maintained
   canonical RH solver as an independent oracle.
3. **Expanding spherical piston:** for `B=0`, `A proportional r^2` and constant
   piston speed, compare shock-radius ratio and piston pressure with a separate
   implementation of Taylor's self-similar ODEs.  Freeze the exact equations,
   branch and reference values before grading.
4. **Formation:** a uniformly accelerating planar gas piston must converge to
   `t_b=2c0/[(gamma+1)a]`, `x_b=c0 t_b`; the cold perpendicular analogue uses
   `t_b=2v_A/(3a)`.
5. **Linear response:** a small piston perturbation satisfies
   `rho'/rho0=u'/c0` and travels at the appropriate gas/reduced fast speed.
6. **Well balance:** a piston comoving with ambient leaves the projected
   ambient steady within frozen `epsilon_wb` density, velocity, pressure,
   magnetic and source-work bounds.
7. **Per-ray conservation:** fixed cell masses, `B_r r^2`, the `b_t` source
   equation, and total energy close against piston work, all ambient source
   terms, cell-append flux and outer flux.
8. **Detected-shock RH:** settled states agree with the canonical solver over
   the declared quasi-perpendicular range.  Other rays are reported as model
   discrepancy, not forced to PASS.
9. **Convergence:** at least three mass, time-step and ray-count refinements
   for `R_sh`, `X`, `M_f,n`, contact pressure, shock formation and assembled
   3-D state.  Numerical error is separated from reduced-model error.
10. **Validity:** report and grade `epsilon_t`, `beta_r`, `beta_n`, `mu_c`,
    `epsilon_Q`, `delta_B`, unsupported area/solid angle, source maxima and
    reduced-versus-canonical characteristic differences.
11. **Assembly:** no cross-interface interpolation; region-wise neighboring
    ray continuity; physical shock normals; complete unsupported-flank budget.
12. **Transactional failure:** any invalid private ray/candidate leaves the
    complete committed epoch, inventory and generation unchanged.
13. **Robustness:** ASan/UBSan, architecture checks, maintained coronal tests
    and baseline particle behavior remain unchanged.
14. **Event validation:** compare the computed Level B shock with the Level A
    observed/prescribed history in standoff, Mach number and compression.  This
    is model-error validation, not a fitted acceptance oracle.

Physical accuracy criteria and exact fixtures are frozen before tuning
resolution, viscosity or source treatment.  PASS/FAIL/SKIP/ERROR totals are
recorded separately.  A required skip leaves BG3D-4 open.

## 12. Known limitations

- Independent rays contain no tangential mass drainage, so propagation-led
  sheath standoff may be overestimated.
- Tangential velocity jumps, switch-on/off structure and normal-field draping
  are not evolved.  Oblique and quasi-parallel rays are approximations.
- Fixed HCI rays cannot evolve transverse material transport.  The selected
  passive ambient transverse component leaves an explicitly measured
  ambient/interface approximation; moving rays are an expanded closure.
- Ray assembly is not solenoidal by construction.  Cleaning, if selected,
  changes the state and requires its own budgets.
- Ambient-maintaining sources are applied unchanged to shocked material.
- The closure predicts a compression/shock from a prescribed ejecta piston;
  it does not predict eruption onset, reconstruct the ejecta interior or solve
  global three-dimensional MHD.
- Observed shock agreement is validation.  It does not authorize resetting or
  forcing the solved shock to the observed history.

The full three-dimensional evolved MHD sheath with an immersed moving boundary
would remove some of these limitations and is outside the authorized model.
