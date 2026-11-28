# Wave evolution contracts requiring review

This is a proposed contract for the new four-spectrum model, not an enabled
production implementation. Initial L/R fractions of 0.5 within each branch are
already authorized; subsequent equality is not imposed. Stored quantities are
magnetic joules per bin. The local branch signs are along/opposite B.
The user subsequently selected material magnetic tubes, fixed initial
magnetic J/m³ per-bin incoming spatial reservoirs, zero exterior k inflow,
and external outgoing spectral accounts. Reflection and cascade are OFF for
the first qualification case. Non-WKB coherence, nonlinear perpendicular
turbulence, shock matching, and provider-compatible material geometry remain
outside that first controlled component.

## Determined ideal-MHD work and ray equations

Let `theta=div U`, `S=b_i (dU_i/dx_j) b_j`, `v_g=U+sigma v_A b`, and
`xi=ln k`, with positive parallel wavenumber magnitude. Ideal induction and
source-free continuity give

```text
D_t ln(v_A) = S - theta/2
d xi/dt     = -S - sigma d_s(v_A)
d ln(k v_A)/dt along the ray = -theta/2.
```

For magnetic energy density per logarithmic wavenumber `w_B`, action density
per logarithmic wavenumber is `N=2 w_B/(k v_A)`. Under the parallel WKB ray/tube
contract, conservative action transport implies

```text
partial_t w_B + div(v_g w_B) + partial_xi[(d xi/dt) w_B]
    = -w_B theta/2.
```

Thus total wave energy is `epsilon=2 W_B`, wave pressure is `P_w=W_B`, and
the integrated total-wave compression term is `-P_w div U`. Magnetic energy
alone receives half that total-wave work. Physical spectral-edge losses must
be recorded; bins leaving the support cannot be clipped or reflected silently.
These equations do not assert that an arbitrary prescribed provider obeys the
ideal background equations or that an oblique full-vector ray remains parallel.

The primary [Jacques (1977), §§II–IV](https://adsabs.harvard.edu/pdf/1977ApJ...215..942J)
distinguishes intrinsic energy transport from pressure-work energy flux. In
particular, the physical total-wave flux is `v_g epsilon + U P_w`. Equivalently:

```text
partial_t epsilon + div(v_g epsilon) = -P_w div U,
partial_t epsilon + div(v_g epsilon + U P_w) = U dot grad(P_w).
```

Recommendation: report the physical energy flux and background mechanical
receipt consistently. If the implementation uses the first, convenient
compression-work form, record the `U P_w` surface term explicitly when closing
the mechanical-energy ledger. Adding both local work expressions would double
count. Particle exchange remains the already adopted `Delta E_B=-Delta K/2`.

## Selected first-case contracts and remaining dependencies

| Contract | Recommendation | Clarification still needed |
|---|---|---|
| Embedded line motion | **Selected:** material magnetic tube, full vertex grid velocity U, frozen labeled magnetic flux, discrete geometric conservation law, and particle coordinate convection U−v_grid. | Native provider-compatible deformation/geometry staging remains to be implemented and qualified. Current static curves plus projected velocity are insufficient. |
| Background compatibility | First qualify ideal analytic manufactured providers. Require an explicit provider contract for induction, continuity and geometry before applying WKB work/ray operators. **Selected 2026-11-28 for the first coupled particle/wave case:** the approved Stage 12 SWCME1D `SHOCK_ONLY` material profile (constant radial wind, inverse-square mass density, induction-compatible Parker field, explicit sustaining force and work; 0.30–1.20 AU; decision D10-1 in `../STAGE10_BACKEND_CAPABILITY_AUDIT.md`). The canonical radial background below is not used for that case. | Native field-line, source and restart admission still rejects the SWCME1D owner; it must be qualified for the wave case. The canonical radial background would still need physical source/forcing equations with their perturbation response before any wave-work use. |
| Spatial boundaries | **Selected:** continuously fixed initial magnetic J/m³ per bin in each incoming component, outgoing characteristics free. Initial reservoirs remain immutable across evolution. | Native ownership/checkpoint of initial reservoirs and time-resolved sign-reversal splitting remain implementation dependencies. |
| Spectral boundaries | **Selected:** fixed modeled edges, zero exterior incoming power; outgoing magnetic energy and total-wave action enter explicit external accounts at crossing time. | Native ledger ownership/checkpoint and convergence remain dependencies. No heating reinterpretation or support extension is selected. |
| Non-WKB reflection | **Selected first case:** OFF. A later extension must derive linearized MHD amplitudes/cross-covariances. | Full-plan reflection remains open; four initial energies do not specify cross-branch coherence. |
| Nonlinear cascade | **Selected first case:** OFF. A later extension must state perpendicular degrees of freedom and governing closure. | Full-plan cascade remains open; no perpendicular spectrum/correlation-length contract has been chosen. |
| Shocks | Treat a shock as an interface with derived wave transmission/reflection and wavenumber jumps, outside the smooth WKB operator. | Which shock-wave interaction model and matching conditions should govern the four spectra? Smooth ray transport through a discontinuity is insufficient. |

Incoming/outgoing means the sign of
`(U + sigma v_A b - v_grid) dot n_out`, separately at each boundary. A cached
outer E-minus condition is not generally incoming: both branches can be outgoing
in an Eulerian super-Alfvénic flow. A material grid changes that classification.
Initial helicity equality applies to the captured initial reservoirs, but does
not impose equal interior spectra after evolution.
Changing to a material grid also requires all three particle movers to use
the relative convection `U-v_grid` in their field-line coordinate. Retaining
the present `U dot tangent` coordinate advection as well as moving vertices
with U would count convection twice.

For constant composition and radial constant U, a profile
`n(r)=a r^-2 + b r^-4 + c r^-6` has a nonzero source-free continuity residual
`U dn/dr + 2 U n/r = -2 U b r^-5 - 4 U c r^-7`.
The canonical Leblanc-type background therefore cannot inherit a source-free
wave-work law solely from its positive density and a radial velocity tensor.
The missing input is its physical source/perturbation contract, not an arbitrary
numerical density floor or a guessed temperature.

The [Chandran et al. non-WKB conservation analysis](https://arxiv.org/abs/1509.01135)
also distinguishes stationary-background wave-action statements from general
time-dependent cases. A closed energy-only reflection rate should not be
declared exact without its additional statistical assumptions.

## Independent qualification sequence

1. **Ideal isotropic expansion:** `U=H x`, `B~exp(-2Ht)`,
   `rho~exp(-3Ht)`, `v_A~exp(-Ht/2)`, `k~exp(-Ht)`, packet volume
   `~exp(3Ht)`, magnetic packet energy `~exp(-3Ht/2)`. Verify conserved
   action `2 E_B/(k v_A)`, exact signed background work and spectral shifts.
2. **Incompressible extension:** `U=(a x,-a y,0)`, B along x,
   `B,v_A~exp(at)`, `k~exp(-at)`, constant volume and packet energy.
   This isolates missing wavenumber transport without a compression source.
3. **Static radial inhomogeneity:** solenoidal `B~r^-2`, `rho~r^-6`,
   `v_A=a r`. Rays have `r~exp(sigma a t)`, `k~exp(-sigma a t)`,
   packet volume `~exp(3 sigma a t)`, and constant packet energy/action.
4. **Arbitrary ALE motion of a static physical medium:** changes to the
   computational grid must not create physical compression work. Verify the
   geometric conservation law, relative characteristic flux and exact boundary
   accounting before coupling to a provider.
5. **Finite-support tests:** move compact spectra through known k edges and
   spatial faces. Reconcile outgoing action, magnetic/total energy and pressure
   work. Test both B/tangent polarities and nonuniform physical cell sizes.
6. **Coupled particles:** a shared wave update changes resident-particle rates.
   Advance competing integrated hazards consistently, or explicitly document
   and converge a frozen-wave operator split. Never clip/skip a finite-budget
   event silently. P/mu/k resolution, statistical weight, step size and energy
   ledgers require joint convergence tests.

Use one authoritative bin-energy state and its exact edges. Any within-bin
reconstruction and energy/action conversion is a numerical representation that
needs its own convergence and balance checks. Native begin/mid/end Simpson
volume quadrature likewise needs a separate geometry-resolution study; native
moment readback alone proves consistent use of that volume, not its accuracy.

The prepare-only implementation and exact proposed native staging boundary are
documented in [MATERIAL_TUBE_WKB_STAGING.md](MATERIAL_TUBE_WKB_STAGING.md).
Its current fixed-direction controlled intervals and manufactured tests must
not be presented as native material-tube evolution or full-plan completion.

The controlled ambient dispersion scan does not choose limits. The reused
default amplitude `deltaB/B=0.3` and outer support near `k rho_p=981` in this
controlled ambient case are study
inputs requiring assessment, not accepted MHD-domain values. The completed
nominal scan retains 29 unresolved high-k right-hand modes. Neither those
numerical failures nor the plotted 1/5/10% references authorize truncation or
production qualification.
