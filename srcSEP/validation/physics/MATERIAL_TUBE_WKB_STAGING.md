# Controlled four-spectrum WKB component and material-tube staging

This document distinguishes the implemented, prepare-only value component from
the proposed native transaction. It does not authorize production evolution or
choose a physical MHD domain cutoff.

## Physical and representation contract

The immutable wave snapshot stores four independent magnetic energies
`E_(sigma,h),j` [J] for each physical cell and fixed positive-k bin. `sigma=+1`
propagates along B, `sigma=-1` against B. Initial left/right energies are equal
within each branch by the user's initialization policy; subsequent equality is
not imposed. Exact edges come from the existing logarithmic bin centers and
declared log width. The collision view and transport reconstruction both use
constant `Ebin/(V*Delta k)` within a bin; neither uses a point sample at its
center as an integrated energy.

For ideal, small-amplitude nondispersive Alfvén waves, let
`b=B/|B|`, `theta=div(U)`, `S=b_i*(partial_j U_i)*b_j` and
`v_A=|B|/sqrt(mu0*rho)`. With positive parallel k magnitude,

```text
D_t ln(v_A) = S - theta/2
d ln(k)/dt  = -S - sigma*b.grad(v_A)
epsilon_total = 2*W_B,       P_wave = W_B
action_bin = integral_bin 2*g_B(k)/(k*v_A) dk.
```

These follow from ideal induction, source-free continuity, and wave action;
they are not inferred from the availability of positive B/rho/U samples.
For packet-integrated magnetic energy, compression gives
`d ln(E_B)/dt=-theta/2`. The Eulerian/log-k density equation is

```text
partial_t(w_B) + div[(U+sigma*v_A*b)*w_B]
  + partial_ln(k)[(-S-sigma*d_s v_A)*w_B] = -theta*w_B/2.
```

The physical total-wave flux also contains `U*P_wave`. Thus the background
compression receipt is `-2*Delta E_B,work`, whereas the background mechanical
receipt for a physical control volume subtracts the outward `U*P_wave` surface
term as well. Material side walls have zero relative advective flux but can
have nonzero pressure work. All caps and side walls must enter that ledger.
This distinction follows the wave-pressure/action treatment of
[Jacques (1977)](https://adsabs.harvard.edu/pdf/1977ApJ...215..942J) and the total
wave flux in [van der Holst et al. (2014), Eq.27](https://doi.org/10.1088/0004-637x/782/2/81).

## Selected first-case boundaries

Spatial incoming characteristics continuously receive the initial magnetic
energy density [J/m³ per bin] from the reused initializer. There are four
independent reservoir arrays at each end. Outgoing characteristics leave.
The signed speed is `(U+sigma*v_A*b-v_grid).n_out`; a permanently forced outer
E-minus is not the selected policy. Capture is allowed only from initialization
revision1. The reservoir must survive restart unchanged even after the interior
develops unequal polarizations.
The shared API validates a versioned reservoir-content identity covering exact
grid, all four densities, line/side, epoch and origin/configuration/provenance.
This catches accidental public scratch/wire mutation. A self-consistent checksum
does not prove historical origin; the future immutable outer owner and trusted
capture/checkpoint path must preserve the actual initialized reservoir.

The modeled k range stays fixed. Incoming power from outside it is zero.
Outgoing magnetic energy and total-wave action go to external accounts at
their crossing time. They are neither clipped nor called heat. Reflection and
cascade are disabled for this first qualification case. Their full-plan
closures remain open.

## Implemented prepare-only API

`util/sep_polarized_wkb_transport.h` supplies:

- `CaptureInitialReservoir` and `EvaluateBoundaryFlux`, with explicit signed
  geometry/velocity integrals and separate energy/action/pressure-work fluxes;
- `BeginExternalAccount`, with immutable origin/last snapshot identity;
- `PrepareSpectralWork`, an independently testable frozen/common-profile
  characteristic work/remap stage; and
- `PrepareInterval`, conservative donor-cell spatial transport followed by the
  source/ray stage, returning a canonical `PrepareEvolutionCandidate` plus
  ledgers and a staged external account.

No function publishes native arrays, advances a particle, or owns a second
physical wave state. `qualifiedForProduction` remains false. Fixed exact grid,
cell namespace/order, initializer provenance, and initial projection ledgers
are preserved by the shared candidate factory. Operation provenance belongs
to the returned operation ledger.

Inputs contain actual cap area **normals**, not an assumption that a material
cap stays perpendicular to B. They explicitly provide endpoint B/rho/U/grid
velocities; time-integrated plasma/grid/Alfvén volume transports; mass and
action transports; full velocity-gradient and Alfvén-gradient integrals;
true initial/final cell mass; and lateral swept volume. Supplied geometry must
satisfy discrete GCL, constant signed tube flux, and mass continuity. A midpoint
rho sample times V is not substituted for the true integrated cell mass.
Material mode additionally requires full grid velocity U and zero relative
plasma/mass convection. Sign reversals require a caller-supplied resolved cut.

The current source stage requires fixed B direction during the interval. Its
coefficients describe a frozen/common-profile history; arbitrary temporal
histories need converged subintervals. Within-bin energy is remapped by exact
overlap of the stated reconstruction. Outgoing k-edge work/action is evaluated
at the crossing fraction, so a departed wave receives no later on-grid work.
Donor CFL violations, unresolved nonzero shifts, under/overflow, incompatible
geometry, unknown policies, reflection, and cascade return explicit status.
No positivity clipping, extrapolated exterior wave power, or k floor is used.
The module checks ideal induction and source-free mass transport; it does not
solve background momentum or thermal evolution. The affine manufactured maps
are explicit kinematic inputs and may require prescribed background forces.
Their tests do not claim an autonomous unforced full-MHD solution or define a
production provider's missing source/perturbation law. The static radial
convergence case separately states its force-free B, constant gas pressure,
zero U, and no-gravity equilibrium.

The action residual from spatial/source splitting and `Ebin/Delta k`
projection is recorded as numerical discretization error. Pressure is frozen
at the donor stage; its temporal/splitting error requires independent
convergence. These residuals do not become heating or physical background work.

External nonnegative accounts use `ExactNonnegativeSum`: 34 unsigned64-bit
limbs, integer units `2^-1074` of the declared SI quantity. Every accepted
finite nonnegative binary64 increment is retained exactly; totals above the
largest finite binary64 reject before replacing the account. Nearest scalar,
signed diagnostic remainder, and conservative remainder-rounding bound are
derived views. Their two-float sum is not asserted exact for arbitrary streams.
Future native ownership/checkpoint must retain all limbs and account provenance.
This numerical representation bound is not a physical spectral cutoff.
The source-stage crossing integrals themselves are finite numerical results;
the exact accumulator guarantees lossless addition of each accepted binary64
cell increment, not exact real-arithmetic evaluation of the ray equation.
Ordinary component-ledger sums remain floating convergence diagnostics.

## Neutral kinematic component and proposed native transaction

The now-implemented neutral `util/sep_material_tube_kinematics.*` provides the
local Cauchy/Nanson push-forward and an explicit global-affine cell GCL stage
using actual reference-surface area/first moments, true material mass and
caller-provided time quadrature/budgets. Curved reference cells are supported;
arbitrary nonaffine native interpolation is not supplied. See
[MATERIAL_KINEMATIC_CORE.md](MATERIAL_KINEMATIC_CORE.md) for equations, numerical
bounds and independent MTG01–06 references. This component publishes no native
geometry and does not enable the transaction proposed below.

The repaired fixed analytic startup is not an ideal-MHD evolution provider.
For example, its configured Parker field outside the source sphere has
`Br=C/r²`, `Bphi=-C*Omega*(r-Rs)*sin(theta)/(U*r²)` and purely radial
constant-speed U. Direct substitution gives
`[curl(U cross B)]phi=C*Omega*Rs*sin(theta)/r³`. This is nonzero for a
finite source radius away from the rotation axis, whereas a held-fixed field
has zero time derivative. Identical initial previous/current B fixes a storage
history; it does not eliminate that induction residual. No azimuthal flow or
different winding law is silently added to repair it. A material qualification
must evolve the declared initial field with its actual compatible deformation
and mass map, or use another explicitly specified consistent provider. The
neutral Cauchy/Nanson tests do not supply that native evolution.

The user selected a **material magnetic tube**. Label a fluid element by a and
let `X(a,t)` be the full three-component map, `F=partial X/partial a`, `J=det F>0`.
For a smooth ideal interval,

```text
B1 = F*B0/J,          rho1 = rho0/J,
oriented_cap1 = J*F^(-T)*oriented_cap0,
vertex_velocity = U, particle_coordinate_convection = U-v_grid = 0.
```

These give exact affine manufactured geometries and signed frozen flux. They
also show why moving only along the old tangent, assigning `A=Phi/|B|` to an
arbitrary tilted cap, or imposing a volume repair is insufficient.

A future native `MaterialGeometryPlan` should be an immutable staged component
of the **same outer transaction** as background, particle, and four-wave state.
It should pin the old background/snapshot/geometry identities and datum-layout
generation; expose prepared current/previous vertices and their complete
B/U/rho/position state; store cap area vectors, segments/length/tangent,
true mass, volumes, material grid velocity, flux, quadrature provenance/error,
and all shared face/side transport integrals. It must not publish a second
live owner or infer missing thermodynamic/source equations.

The proposed interface boundary is explicit (these signatures are design,
not existing callable APIs):

```cpp
// Inputs are already pinned immutable provider/geometry realizations.
GeometryPrepareResult PrepareMaterialGeometry(
    const NativeGeometryImage& before,
    const IdealMaterialInterval& providerInterval,
    const BackgroundSnapshot& expectedBackground);

// Returns the complete explicit WKB interval AND target native background.
// No lazy flux initialization or host writes occur during projection.
WkbProjectionResult ProjectControlledWkbInterval(
    const MaterialGeometryPlan& geometry,
    const PolarizedSnapshot& waves,
    const InitialBoundaryReservoirs& initial,
    const ExternalSpectralAccount& external);

// These join the outer particle/wave/background transaction. They do not
// provide a separate public geometry-only commit path.
Status ValidateMaterialGeometryReady(const MaterialGeometryPlan& plan);
Status ConfirmCollectiveMaterialGeometryReady(const MaterialGeometryPlan& plan);
```

`IdealMaterialInterval` must provide the material map/deformation and physical
interval quadrature or return unavailable. `NativeGeometryImage` must pin raw
layout identity, full current/previous primitives, and independently decoded
geometry values. `MaterialGeometryPlan` owns prepared native buffers, exact
shared-face/side integrals, the expected old owner/background pointers and
revision tokens, and its immutable target image. Its WKB projection carries
the actual cap normals and fixed-direction capability explicitly; an arbitrary
curved/rotating interval cannot acquire that capability from a name or flag.
The outer transaction alone may consume this plan at commit after particle,
wave, source, and checkpoint readiness succeed.

Readiness must verify, before any writes:

1. Complete line/segment namespace, provider capabilities and time cuts,
   finite positive J/volume/mass, ideal induction and continuity, signed flux,
   independent geometric conservation, valid all-rank identity/configuration,
   and no active particle/background read phase.
2. Every particle consumes the intended immutable interval and relative
   convection. Moving geometry updates Cartesian position caches from their
   authoritative field-line coordinate; it does not advance that coordinate a
   second time with U. All three movers need this integration proof.
3. Four-wave candidate, retained initial reservoirs, exact external accounts,
   particle edits, physical background/pressure ledgers, lifecycle/event clocks,
   current/previous native primitives, and derived diagnostics are staged
   together. A failed cell, budget, allocation, rank, or readback rejects all.
4. Native buffer/configuration identities remain unchanged at final readiness.
   All allocations, collective agreement and recoverable readback precede an
   allocation-free commit. Old owner lifetime is pinned until completion.

Source injection requires a separate matching contract. An exact full-particle
drifting Maxwell–Jüttner flux through a moving material plane is not
automatically the flux of a gyroaveraged one-dimensional field-line model.
Under a general deformation the cap normal need not align with B, and finite
gyro-orbits cross that plane differently from a guiding center. The native
bridge therefore needs either an explicit aligned-cap domain or a derived
guiding-center boundary/source matching law. Projecting a sampled Cartesian
momentum onto B alone does not establish that law. The controlled WKB component
neither samples particles nor asserts source-bridge qualification.

The cap and lateral sweeps must be independently integrated and compared with
`V1-V0`, not defined as the discrepancy required to pass GCL. Lazy native tube
flux initialization must either be finalized before staging or included in the
same transaction. Existing begin/mid/end Simpson volume quadrature needs its
own geometry-resolution/error evidence; native moment readback alone does not
provide that evidence. Rotating B requires resolved b/G quadrature before
general curved material providers can use the current source stage.

The retained fixed-Parker/nonideal profiles do not automatically meet these
requirements. In particular, constant radial U with a Leblanc-type mixture
`n=a/r²+b/r⁴+c/r⁶` has nonzero source-free continuity residual. A provider's
forcing/source/perturbation contract is required before applying the ideal
work law there. No such contract is guessed here.

## Bounded verification and remaining qualification

Run `srcSEP/test/run_polarized_wkb_tests.sh`. It compiles the same shared
callbacks as strict C++11 with ASan/UBSan. `WKB01`–`WKB08` cover signed initial
reservoirs/polarity/grid motion, finite-band expansion/strain crossing,
independent compact-spectrum refinement, full radial spatial/ray
characteristics, controlled ALE GCL, fail-before-publication inputs, material
side-wall pressure cancellation, and exact multiscale account arithmetic.
Expected profiles use independent polynomial/Gaussian antiderivatives or
analytic ideal-MHD characteristics, not production remapping.

These are controlled synthetic gates. They do not qualify native material
geometry, event-rate changes during particle/wave exchange, fresh-process full
checkpoint continuation, reflection/cascade/shocks, observational comparison,
or the MHD physical domain. The default initial `deltaB/B=0.3`, high-k controlled
ambient `k*rho_p` near981, and unresolved high-k warm-dispersion points remain
study inputs/issues. No automatic truncation or accepted numerical/physical
cutoff follows from the plotted1/5/10% study reference contours.
