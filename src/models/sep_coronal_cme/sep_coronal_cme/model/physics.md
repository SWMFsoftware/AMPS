# Shared Stand-Alone Low-Coronal CME-Shock Model for 3-D and Field-Aligned SEP Transport

## Document status

**Specification revision:** 2026-09-28-r4  
**Application input schema:** 5  
**Field-line bundle schema:** 3

Schema 5 remains a **pre-freeze draft** until the existing review-3 release
gates in Sections 14, 16, and 17 pass and the review-4 documentation and
future-capability dispositions are traceable in those contracts. The changes
below therefore refine the proposed
schema-5 contract before implementation; they do not reinterpret a released
schema. Once the first schema-5 deck is accepted as a release artifact, any
subsequent incompatible key or enum change requires a new schema number.

This document is the canonical physics and software contract for the shared
AMPS `src/models/sep_coronal_cme` model. `srcSEP3D` links the model through a
thin AMPS-facing adapter. Field-aligned `srcSEP` does not link or invoke the
three-dimensional model; it consumes an immutable, versioned field-line bundle
through neutral `sep_common` records. Neither application owns or may fork the
background, CME-front, shock, turbulence, source, or field-line-reduction
physics specified here. The shared model permits a prescribed CME
candidate front to begin at an apex height of
`1.01--1.05 R_sun` and propagate through the corona and heliosphere without
calling SWMF/AWSoM, MAS, ENLIL, EUHFORIA, or another external MHD solver during
the run. The front becomes a shock only where and when the local fast-mode and
Rankine--Hugoniot admissibility tests pass. The same authoritative
three-dimensional background and candidate-front/shock state
can be reduced onto one or more open magnetic field lines, exported through a
versioned interface, and consumed by `srcSEP` for one-dimensional Parker or
focused field-aligned transport. The export is a controlled dimensional
reduction of the 3-D model, not a separately reconstructed Parker background.

This revision incorporates the completed physics reviews: separate wind and
shock thermodynamics; authoritative closed-field plasma; a finite-shell
Schatten current-sheet layer; a conservative longitude-dependent Parker map;
structured fast/slow wind; delayed patchwise shock activation; fully specified
ellipsoid kinematics; bounded critical-Mach and near-parallel jump handling;
magnetogram filtering; smooth/broken scattering and source closures; and
non-lossy connection status for field-line export. Review 2 additionally
requires quantitative current-sheet radialization, absolute open-flux
comparison, an observation-constrained kinematic-wind branch, local
separatrix stress diagnostics, a topology-safe PFSS/SCS transition policy,
force-based turbulence validity, policy-dependent critical-Mach coverage, one
solar-rotation authority, and explicit upstream-release semantics. Review 3
adds a finite physical release surface, a two-zone observation-constrained
wind, honest open/closed-interface policies, quantitative consumer-impact
budgets, joint background/shock calibration, independently identified
field-aligned wind products, and an executable modular-document authority.
Review 4 sharpens the physical interpretation of the finite release surface,
separates diagnostic Péclet depth from a return probability, and records the
additional physics needed before a momentum-dependent release surface or
post-release recycling could be enabled. It also makes explicit that
self-generated waves, gradient/curvature drift, time-dependent CME deflection
and body rotation, and a flare-associated source are future model branches,
not hidden capabilities of schema 5. Finally, it defines the observational
content required for a future continuous wind-plausibility qualification while
preserving the present positivity, continuity, D7, and momentum-residual
gates. These review-4 clarifications do not change application input schema 5
or field-line bundle schema 3 and do not claim that any of the future branches
is implemented.
Where added realism would require a global MHD or kinetic calculation, the
document defines a typed failure or an explicitly bounded approximation
instead of implying that the stand-alone model solves that missing physics.

### Review-correction disposition

Accepted review findings are requirements, while partially correct proposals
are incorporated only after their physics and scope are corrected. Every
disposition is traceable in the specification and release tests:

{{GENERATED:REVIEW_DISPOSITION}}

Additional corrections needed for a defensible production baseline are also
normative: the PFSS/SCS join is divergence preserving and spatially resolved;
PFSS and composite topology are distinct; repeated observers use edge-defined
energy channels and finite empty-bin output; and 1-D `srcSEP` inputs are exact
reductions of the prepared 3-D background rather than reconstructed models.

### Acronyms and model names

| Term | Meaning |
|---|---|
| AMPS / AMR / API | Adaptive Mesh Particle Simulator / adaptive mesh refinement / application programming interface |
| CME / SEP / GLE | coronal mass ejection / solar energetic particle / ground-level enhancement |
| MHD / EOS / IMF | magnetohydrodynamics / equation of state / interplanetary magnetic field |
| PFSS / SCS / HCS / CSSS | potential-field source-surface / Schatten current-sheet / heliospheric current sheet / current-sheet source-surface |
| DSA / QLT / SOQLT | diffusive shock acceleration / quasilinear theory / second-order quasilinear theory |
| WKB / WSA / RTV / CIR / CFL | Wentzel--Kramers--Brillouin / Wang--Sheeley--Arge / Rosner--Tucker--Vaiana / corotating interaction region / Courant--Friedrichs--Lewy |
| HCI / HEEQ | heliocentric inertial / heliocentric Earth equatorial coordinate frame |
| GONG / HMI / EUV | Global Oscillation Network Group / Helioseismic and Magnetic Imager / extreme ultraviolet |
| SWMF / AWSoM / MAS | Space Weather Modeling Framework / Alfvén Wave Solar atmosphere Model / Magnetohydrodynamics Around a Sphere |
| ENLIL / EUHFORIA | named heliospheric MHD model / European Heliospheric Forecasting Information Asset |

The model is intended to be:

1. physically constrained and internally self-consistent at the level claimed;
2. inexpensive enough for baseline SEP transport and acceleration studies;
3. deterministic and reproducible from a versioned input file;
4. compatible with the existing AMPS mesh, particle, sampling, restart, and
   population-control infrastructure;
5. structured so that an imported time-dependent MHD background can later
   replace the analytic background without replacing the particle-transport
   implementation; and
6. able to drive matched 3-D and 1-D calculations from the same background,
   shock, turbulence, source, species, coordinate, and provenance contracts.

The model is **not** a substitute for a global coronal MHD simulation. It does
not predict CME initiation, shock formation, or the global downstream sheath
from first principles. It prescribes a time-dependent candidate-front geometry
and kinematics, evaluates whether each patch is an admissible MHD fast shock in
the analytic upstream corona, and uses only source-eligible shock patches as
an SEP source. This distinction is part of the model contract and must be
preserved in code, output metadata, and publications.

---

## 1. Scientific objective and scope

The primary objective is to model SEP release and transport from a CME-driven
shock in the low corona to one or more observers. The computational domain may
begin immediately above the photosphere, while AMPS' registered internal
spherical boundary represents the solar surface at `r = R_sun`.

The background and particle mesh are defined on the physical side of the solar
sphere, `R_sun<r<=R_dom`, where `R_dom=domain.outer_radius_m` is the exact
spherical escape boundary enclosed by the Cartesian AMR box. The separately
configured **qualified source inner radius** is

\[
  R_{\rm in}=1.01\text{--}1.05\,R_\odot,
  \qquad R_\odot\equiv R_{\rm sun},
\]

and only points with `r>=R_in` may support the low-coronal shock source; it is
not a second absorbing boundary. The configurable outer boundary extends
to at least 1 AU for Sun-to-Earth transport. The default baseline should use
`R_in = 1.05 R_sun`; the
`1.01 R_sun` option is retained for controlled sensitivity experiments because
the density, temperature, magnetic topology, and characteristic wave speeds
have especially steep and observation-dependent gradients at that height.

The model contains six coupled but separately testable components:

1. a potential-field source-surface (PFSS) low-coronal magnetic field;
2. a finite-shell Schatten-type current-sheet (SCS) field and a discrete
   heliospheric-current-sheet sector map;
3. composition-consistent open-tube wind and closed-field hydrostatic plasma;
4. a solar-anchored, expanding ellipsoidal CME piston and candidate front;
5. Parker or focused SEP transport using locally admissible shock patches as a
   moving source; and
6. a conservative field-line reduction and exchange layer through which one
   or more selected open lines are produced by the shared model, exported by a
   `srcSEP3D` adapter, and used by `srcSEP`.

The production magnetic background has three separately identified analytic
authorities and a resolved PFSS/SCS join. PFSS is solved to a radial outer
boundary `R_b`; the composite field transitions from PFSS to a finite-shell
SCS representation beginning at `R_i`; and the SCS solution is made radial at
`R_scs`, where the conservative rotating-footpoint Parker extension begins.
Direct PFSS-to-Parker coupling is retained only as a named analytic
verification mode. `R_b`, `R_i`, and `R_scs` remain distinct inputs in the
production composite model. No numerical value of `R_scs`, including
`2.5 R_sun`, is accepted solely because it is a conventional coupling radius:
the realized shell must report its complete spectral attenuation, while its
absolute outer zonal power and the completed field's registered latitude-
flatness metrics must pass.
If winding is required to begin at `2.5 R_sun` but that radius fails the SCS
radialization gate, the existing radial Parker map is inapplicable. A future
general solenoidal pushforward may use a distinct winding radius only after
its nonradial mapping, Jacobian, wind coupling, and topology tests pass.
In the no-SCS analytic verification branch, `R_i` and `R_scs` are derived as
`R_b` rather than supplied as independent radii. None of these radii is a
universal physical constant; an event calculation must constrain them using coronal
holes, streamer/current-sheet morphology, and open-flux diagnostics [2--4,
34--36, 54, 55].

### 1.1 Radial domain of applicability

The computational outer radius and the radial range over which the prescribed
shock remains physically credible are different limits. The analytic
background and SEP transport may extend to an observer at 1 AU even after the
low-coronal shock source has been deactivated. The recommended baseline is:

| Radial interval | Model role | Permitted interpretation |
|---|---|---|
| `R_sun--R_in` | fully initialized precipitation/guard shell with no SEP injection | boundary transport only; outside source-validation claims |
| `R_in--R_i` | signed PFSS field, PFSS/composite topology diagnostics, open wind or closed hydrostatic plasma selected by composite topology, and ellipsoidal candidate front | primary low-coronal source domain |
| `R_i--R_scs` | resolved PFSS/SCS join followed by finite-shell SCS field, discrete sectors/HCS, continued open-tube wind, and candidate front | outer-coronal current-sheet domain |
| `R_scs--min(20 R_sun,R_zero)` when nonempty | conservative Parker extension, continued wind, and candidate front/shock | baseline early-heliospheric acceleration domain; an empty interval is reported rather than written with inverted limits |
| `20--min(30 R_sun,R_zero)` when nonempty | transition/sensitivity interval | shock results require a stated model-form uncertainty; an empty interval is valid and explicit |
| configured zero-source radius through `1 AU` | composite SCS/Parker background and particle transport after low-coronal injection; Parker only where `r>=R_scs` | baseline Sun-to-observer transport domain; the source stops at `20 R_sun` nominally or up to `30 R_sun` only in the stated sensitivity branch, independently of the magnetic-authority boundary |
| beyond `1 AU` | unsupported by the baseline validation | requires an explicitly extended and independently validated heliospheric model |

For the first production release, the domain should extend to `1 AU`, while
the active low-coronal shock source should terminate at a configurable apex
radius whose recommended value is `20 R_sun`; `30 R_sun` is a sensitivity
case, not a new default. A surface patch is locally source-inactive whenever it
fails its configured fast-mode or supercritical criterion; it may reactivate
only if its recorded Mach history again satisfies the criterion before the
global source envelope reaches exact zero. Particles injected earlier continue
to propagate to all observers.

The approximate `20 R_sun` source termination is not a numerical limit. Beyond the
upper corona, CME propagation becomes increasingly controlled by interaction
with the ambient solar wind and is commonly represented by a drag-based
heliospheric model [26--28]. Extending the shock itself to 1 AU therefore
requires a new branch containing at least a propagated CME driver, a separately
evolved shock surface and standoff distance, local shock weakening in the
structured wind, and validation against heliospheric-imager tracks and in-situ
shock arrival. A drag law alone supplies CME kinematics; it does not supply the
local compression, obliquity, criticality, or SEP-injection efficiency.

The field-line exchange API must not encode `1 AU` as a hard limit. It may
export a longer line for a separately qualified `srcSEP` study, including a
future 5-AU transport calculation, but the bundle must record that the present
composite PFSS--SCS--Parker wind and turbulence realization is validated only
to 1 AU. A claim
beyond 1 AU requires additional validation of solar-wind structure, stream
interactions, turbulence and scattering evolution, and outer-heliospheric
particle physics.

### 1.2 Matched 3-D and field-aligned use

The 1-D mode is intended for controlled, computationally inexpensive
field-aligned transport along one or more selected magnetic connections. It
has three principal uses:

1. production 1-D SEP calculations along observer-connected field lines;
2. direct Parker-versus-focused and 1-D-versus-3-D intercomparisons; and
3. ensembles over multiple shock-connected field lines without allocating the
   complete 3-D AMR particle domain.

Each exported line is an oriented curve `x_l(s,t)` with outward arc length
`s`; schema 5 freezes the background geometry at its declared trace time, but
the time argument remains explicit in the exchange model.
The coordinate orientation is geometric and always points away from the Sun;
magnetic polarity is stored independently. Therefore

\[
 \left.\frac{\partial\mathbf x_l}{\partial s}\right|_t=\hat{\mathbf t}_l,
 \qquad
 \hat{\mathbf b}=\sigma_{B,l}\hat{\mathbf t}_l,
 \qquad \sigma_{B,l}\in\{-1,+1\}.
\]

This convention prevents a polarity reversal from silently reversing the
meaning of the inner and outer particle boundaries. Only open lines that reach
the requested outer radius may be exported for Sun-to-observer transport.
Closed lines are reported diagnostically and rejected as `srcSEP` transport
inputs. The first production focused-transport model is sector confined: an
exported line may approach but may not cross an unresolved zero-thickness HCS
or an unverified PFSS/SCS field-direction discontinuity.

---

## 2. Definitions, frames, and notation

### 2.1 Coordinates

All kernels use SI units internally. Let

\[
  \mathbf{x}=r\,\hat{\mathbf e}_r
\]

be the position relative to the configured heliocentric origin. Spherical
coordinates are `(r, theta, phi)`, where `theta` is colatitude and `phi`
increases in the right-handed sense about the configured solar-rotation axis.
The coordinate frame, epoch, solar north, zero-longitude direction, and
rotation convention are mandatory metadata.

With rigid field-line rotation, the analytic background is naturally stationary
in a frame rotating at `Omega_F`. A latitude-dependent rotation law has no
single corotating frame and must instead be evaluated in an inertial or a
declared fixed-rate frame as a time-dependent mapping. Particle positions and
observer ephemerides may be stored in an inertial frame, but every conversion
must be explicit and tested. A vector must never be identified as HCI, HEEQ,
Carrington, or an application-local frame solely because it is Sun-centered.

The first production profile accepts one rigid **sidereal** `Omega_F` expressed
relative to the declared inertial heliocentric axes. A synodic input is legal
only when an ephemeris-owned orbital angular rate and sign convention convert
it to that sidereal value before any centrifugal, longitude-map, source, or
observer calculation. Latitude-dependent rotation is an analytic-verification
branch until a time-dependent 3-D background and closed-corona provider are
qualified; it is not a production shortcut.

`[solar_rotation]` is the sole configuration authority for `Omega_F`, its
convention, and any synodic-conversion asset. The closed-plasma potential,
open-tube frame, Parker winding, observer transformations, and field-line
export receive the same resolved sidereal value from the immutable prepared
configuration. Provider-local copies may appear in manifests for regression
checking, but they are derived values and cannot be entered independently.

### 2.2 Plasma variables

The upstream and downstream states use subscripts 1 and 2. The principal
quantities are

| Symbol | Meaning | SI unit |
|---|---|---:|
| `rho` | total mass density | kg m^-3 |
| `n_s` | number density of species `s` | m^-3 |
| `p` | total thermal pressure | Pa |
| `T_s` | species temperature | K |
| `u` | bulk plasma velocity | m s^-1 |
| `B` | magnetic field | T |
| `v_A` | Alfvén speed | m s^-1 |
| `a_w` | effective polytropic/isothermal wind sound speed | m s^-1 |
| `c_s` | thermodynamic adiabatic sound speed used by waves and shocks | m s^-1 |
| `c_f` | fast-magnetosonic speed | m s^-1 |
| `w_out`, `w_in` | wave energies propagating geometrically outward/inward | J m^-3 |
| `w_+`, `w_-` | derived wave energies parallel/antiparallel to signed `B` | J m^-3 |
| `lambda_parallel` | parallel particle mean free path | m |
| `gamma_w` | effective open-wind polytropic index | dimensionless |
| `gamma_ad` | thermodynamic ideal-MHD adiabatic index | dimensionless |
| `gamma_c` | optional closed-plasma polytropic index | dimensionless |
| `sigma_B` | magnetic polarity relative to the outward tube tangent | `-1` or `+1` |
| `R_b` | PFSS radial outer-boundary radius | m |
| `R_i` | PFSS/SCS interface radius | m |
| `R_scs` | radial SCS boundary and Parker-coupling radius | m |
| `R_w` | Parker-winding start radius; equals `R_scs` for the schema-5 radial map | m |

For a proton-electron baseline with optional alpha particles,

\[
  \rho = m_p n_p + m_\alpha n_\alpha + m_e n_e,
  \qquad n_e=n_p+2n_\alpha,
\]

and

\[
  p=k_B(n_pT_p+n_eT_e+n_\alpha T_\alpha).
\]

Electron mass may be neglected in `rho` only if that approximation is recorded
in the resolved configuration. It must not be neglected in charge neutrality.

The wind and shock indices are intentionally distinct:

\[
 a_w^2=\gamma_w\frac{p}{\rho},
 \qquad
 c_s^2=\gamma_{\rm ad}\frac{p}{\rho}.
\]

`a_w` controls the steady wind critical point and represents the effective
heating closure. `c_s` controls rapid compressive disturbances and the MHD
jump conditions. The downstream shock state is not constrained to lie on the
effective wind polytrope.

### 2.3 Shock frame

At a point on the moving shock, `n_hat` points from the downstream region into
the upstream region. If the surface normal speed is `V_sh,n`, the plasma
velocity in the local shock frame is

\[
  \mathbf{v}_i=\mathbf{u}_i-V_{\rm sh,n}\hat{\mathbf n}.
\]

Normal and tangential components are

\[
  v_{in}=\mathbf{v}_i\cdot\hat{\mathbf n},\qquad
  \mathbf{v}_{it}=\mathbf{v}_i-v_{in}\hat{\mathbf n},
\]

with analogous definitions for `B_n` and `B_t`.

### 2.4 Particle kinematics, species, and momentum frame

For compiled species mass `m` and signed charge `q_s`, kinetic energy, momentum,
speed, and rigidity use the exact special-relativistic relations

\[
 E_k(p,m)=\sqrt{p^2c^2+m^2c^4}-mc^2,
 \qquad
 p(E_k,m)=\frac{\sqrt{E_k(E_k+2mc^2)}}{c},
\]

\[
 v(p,m)=\frac{pc^2}{E_k+mc^2},
 \qquad
 \mathcal R=\frac{pc}{|q_s|}.
\]

Here the species charge is henceforth `q_s`; bare `p` in particle/source/
transport equations is momentum magnitude, whereas `p` in plasma and MHD
equations is thermal pressure. These conventional symbols never denote both
quantities in the same equation. The DSA spectral slope is written `q_DSA`,
not `q_s`.

Rigidity is in volts when SI momentum and coulombs are used. “Energy per
nucleon” means `E_k/A`, where integer nucleon count `A` is supplied by
validated species metadata. It is undefined for electrons, unresolved isotope
mixtures, and any compiled symbol without an unambiguous `A`; such a selection
is rejected rather than estimated from the floating-point mass.

Unless a record says otherwise, mover variables `p`, `mu`, the source spectrum,
and scattering coefficients are defined in the **local plasma rest frame**,
while position and time are coordinates in the declared transport frame. A
spectrum defined in another supported momentum frame is sampled there and
transformed into the local-plasma mover frame through an exact Lorentz boost
of each particle four-momentum; a Galilean energy shift is not used for
relativistic particles. Such a source frame is selectable only when the
transformed distribution remains representable by the mover's isotropic or
gyrotropic state. An exact boost of individual particles does not justify
discarding gyrophase dependence or redefining mover `p,mu` in the coordinate
frame. The source ledger stores sampled-frame and local-plasma four-momenta
and, as separately derived diagnostic views, shock-, detector-, or Cartesian
storage-frame four-momenta. Those diagnostic/storage views never overwrite
the phase-space coordinates advanced by the mover.

For the generic transformation below, let the sampled inertial frame move
with velocity `V` relative to a chosen local orthonormal destination frame,
and let `E'_tot=E'_k+mc^2` and
`gamma_V=(1-|V|^2/c^2)^(-1/2)`.

The calculation is performed in local orthonormal inertial tetrads after any
rotating-coordinate basis/velocity conversion; a non-inertial coordinate
velocity is never inserted directly into a Lorentz formula. For source
insertion the destination is the local plasma frame. If AMPS requires a full
Cartesian stored momentum, the adapter reconstructs the sampled gyrophase and
transforms that individual four-vector to the documented storage frame; this
is a representation conversion, not a change to the mover's `p,mu` frame.

The boost is

\[
 E_{\rm tot}=\gamma_V(E'_{\rm tot}+\mathbf V\cdot\mathbf p'),
\]

\[
 \mathbf p=\mathbf p'+
 \left[
  \frac{\gamma_V-1}{|\mathbf V|^2}(\mathbf V\cdot\mathbf p')+
  \frac{\gamma_V E'_{\rm tot}}{c^2}
 \right]\mathbf V,
\]

with the continuous identity limit at `V=0`. The inverse uses `-V`. The frame
velocity is evaluated at the sampled patch and epoch and is stored with the
source event.

---

## 4. PFSS coronal magnetic field

### 4.1 Governing equations

PFSS assumes that electric currents do not contribute to the modeled
large-scale coronal field:

\[
  \nabla\times\mathbf B=0,
  \qquad \nabla\cdot\mathbf B=0.
\]

Therefore

\[
  \mathbf B=-\nabla\Psi,
  \qquad \nabla^2\Psi=0.
\]

In the spherical shell `R_sun <= r <= R_b`, the scalar potential is expanded
as

\[
 \Psi(r,\theta,\phi)=
 \sum_{\ell=1}^{L_{\max}}\sum_{m=-\ell}^{\ell}
 \left(a_{\ell m}r^\ell+b_{\ell m}r^{-(\ell+1)}\right)
 Y_{\ell m}(\theta,\phi).
\]

The monopole term is excluded because the net magnetic flux through the solar
surface must vanish to numerical precision.

### 4.2 Inner magnetic boundary

At the photosphere,

\[
 B_r(R_\odot,\theta,\phi)=B_{r,0}(\theta,\phi).
\]

`B_r,0` may be supplied in either of two ways:

1. **Analytic-harmonic baseline:** input coefficients define a dipole,
   quadrupole, and optional higher harmonics. This mode is data-independent and
   is required for exact regression tests.
2. **Magnetogram mode:** a local GONG or HMI synoptic/synchronic map is
   transformed into spherical-harmonic coefficients. The map, calibration,
   pole treatment, flux-balance correction, epoch, and checksum become part of
   the run manifest.

The code may remove a measured net-flux offset only through a documented
flux-balance operation. It must report the removed monopole flux and the ratio
of that correction to the unsigned flux.

Magnetogram processing must also specify the native angular grid, cell-area
quadrature, polar treatment, and an admissible `L_max` for that grid. The
available spectral choices are `none` and a named heat-kernel filter

\[
 \widetilde g_{\ell m}=g_{\ell m}
 \exp\left[-\frac{\ell(\ell+1)}
                  {\ell_a(\ell_a+1)}\right].
\]

Filtering is optional for exact low-order analytic inputs but is part of the
convergence study for a magnetogram used as low as `1.01 R_sun`. Raw and
filtered coefficients, the transfer function, map/grid metadata, weighted
reconstruction errors, unsigned-flux change, neutral-line displacement, and
open/closed topology changes are recorded. Filtering is a declared smoothing
of the measured boundary; it must never be described as an exact
reconstruction of the raw map. The harmonic truncation, grid quadrature, and
direct-versus-cached PFSS solution are qualified by convergence tests [42].

### 4.3 PFSS radial/source-surface boundary

The PFSS field is constrained to be radial at `R_b`:

\[
 B_\theta(R_b)=B_\phi(R_b)=0.
\]

Equivalently, the potential is constant on the source surface. For every
non-monopole harmonic,

\[
 a_{\ell m}R_b^{\ell}
 +b_{\ell m}R_b^{-(\ell+1)}=0,
\]

so that

\[
 b_{\ell m}=-a_{\ell m}R_b^{2\ell+1}.
\]

If the raw lower radial field is expanded in the same normalized basis,

\[
 B_{r,0}(\theta,\phi)=
 \sum_{\ell m}g_{\ell m}Y_{\ell m}(\theta,\phi),
 \qquad
 g_{\ell m}=\int B_{r,0}Y_{\ell m}^{*}\,d\Omega,
\]

define the selected coefficients

\[
 g_{\ell m}^{\rm sel}=g_{\ell m}
 \quad\hbox{for no spectral filter},
 \qquad
 g_{\ell m}^{\rm sel}=\widetilde g_{\ell m}
 \quad\hbox{for the heat-kernel filter}.
\]

The PFSS solve uses `g_lm_sel`, so

\[
 a_{\ell m}=
 -\frac{g_{\ell m}^{\rm sel}}
 {\ell R_\odot^{\ell-1}
 +(\ell+1)R_b^{2\ell+1}R_\odot^{-(\ell+2)}},
 \qquad
 b_{\ell m}=-a_{\ell m}R_b^{2\ell+1}.
\]

For a real magnetic field, the complex coefficients obey the corresponding
conjugate-symmetry relation, for example
`g_sel_(l,-m)=(-1)^m conj(g_sel_(l,m))` under the standard
complex-harmonic convention.
The selected raw or filtered boundary representation therefore determines
`a_lm` and `b_lm` without a numerical elliptic solve. Only the unfiltered,
untruncated representation is described as reproducing the raw map exactly;
every filtered or truncated solve reports its boundary residual.

The implementation should use normalized associated Legendre functions and an
explicit normalization convention; the convention is part of the fingerprint.

### 4.4 Open and closed topology

Field-line integration uses

\[
  \frac{d\mathbf x}{ds}=\pm\frac{\mathbf B}{|\mathbf B|}.
\]

A PFSS line that reaches `R_b` is `pfss_open_to_Rb`; a PFSS line that returns to
the solar boundary before `R_b` is `pfss_closed_below_Rb`. Independently, a
line in the production composite field that reaches `R_i` from below is
`composite_open_to_Ri`; a line whose two directions return to the solar
boundary below `R_i` is `composite_closed_below_Ri`. These classifications can
differ when `R_i<R_b`, because the SCS construction deliberately replaces the
PFSS field above `R_i` and can open a line that the unused PFSS continuation
would close.

Every background query returns both stable classifications (or `separatrix`)
and the generation in which they were made. Composite topology controls
open-wind versus closed-plasma ownership, particle-corridor eligibility, and
shock-patch topology. A target-speed asset separately declares whether its
coronal-hole boundary, expansion factor, and calibration use PFSS topology at
`R_b` or composite topology at `R_i`; applying it outside that declared
topology is fatal. Classification must converge with trace tolerance and may
not be inferred from the AMR cell containing the point.

The active SEP corridor must be centered on an open line. Configuration
validation fails if a requested Sun-to-observer line closes below `R_i`.
Shock diagnostics may be calculated over both topologies, but source injection
from a closed patch is disabled by the recommended baseline policy.

### 4.5 Authoritative closed-field plasma

The transonic-wind provider owns plasma variables only on open tubes. A
separate closed-field provider owns `rho`, `p`, species temperatures,
composition, and velocity on closed tubes. A visualization-only value may not
be consumed as shock input.

For the first release, a static closed provider is defined in the same rigidly
corotating frame with the resolved angular velocity `Omega_F` owned by
`[solar_rotation]`. The closed provider consumes that immutable value read-only;
there is no `Omega_c` input or second physical authority. Its effective
potential is

\[
 \Phi_{\rm eff}(\mathbf x)=-\frac{GM_\odot}{r}
 -\frac12|\boldsymbol\Omega_F\times\mathbf x|^2,
 \qquad
 \Delta\Phi_{\rm eff}=\Phi_{\rm eff}(\mathbf x)
 -\Phi_{\rm eff}(\mathbf x_0).
\]

Here `x_0` is the topology-consistent reference footpoint or reference point
on the configured sphere `r=R_c0`; the base density, composition, and species
temperatures are defined at that same point. The provider may not mix a base
state at `R_c0` with a potential difference referenced to another radius.
Schema 5 names this distinct radius
`closed_field_plasma.reference_radius_m`.

For the first-release `isothermal-hydrostatic` closure, the composition-aware
EOS supplies `rho_0` and `p_0` at `R_c0`, and

\[
 c_T^2=\frac{p_0}{\rho_0},\qquad
 \rho(\mathbf x)=\rho_0
 \exp\left[-\frac{\Delta\Phi_{\rm eff}}{c_T^2}\right],\qquad
 p(\mathbf x)=p_0
 \exp\left[-\frac{\Delta\Phi_{\rm eff}}{c_T^2}\right].
\]

Species temperatures remain at their declared base values. This form avoids
an ambiguous mean molecular mass: for a quasi-neutral proton-electron plasma,
for example, `p = n_p k_B(T_p+T_e)` while `rho` is approximately `m_p n_p`.
The baseline takes `u'=0` in this corotating frame; in an inertial frame the
closed plasma corotates. The manifest reports the maximum centrifugal-to-
gravity acceleration ratio. A gravity-only option is an explicitly named
small-rotation approximation, not the default. A latitude-dependent rotation
law has no globally static closed-corona frame and is rejected by this provider
unless a separately validated time-dependent closed model is selected.

An optional `polytropic-hydrostatic` closure uses its own `gamma_c>1`, not
`gamma_w`:

\[
 h=\frac{\gamma_c p}{(\gamma_c-1)\rho},\qquad
 h(\mathbf x)=h_0-\Delta\Phi_{\rm eff},
\]

\[
 \rho(\mathbf x)=\rho_0
 \left(\frac{h(\mathbf x)}{h_0}\right)^{1/(\gamma_c-1)},\qquad
 p(\mathbf x)=p_0
 \left(\frac{h(\mathbf x)}{h_0}\right)^{\gamma_c/(\gamma_c-1)}.
\]

Each declared species temperature follows the same polytropic scaling,

\[
 T_s(\mathbf x)=T_{s0}
 \left[\frac{\rho(\mathbf x)}{\rho_0}\right]^{\gamma_c-1},
\]

and the composition-aware sum of species pressures must reproduce `p`. The two
solar footpoints of a closed loop must predict the same pressure at a common
reference point within tolerance; incompatible footpoint boundary data are
rejected rather than blended.

Every requested point must satisfy `h>0`; otherwise preparation fails. These
closures enforce field-aligned hydrostatic balance but do not solve loop
heating, thermal conduction, evaporation, siphon flow, reconnection, or
CME-driven opening. Resolving dynamic coronal opening requires a global MHD
model rather than either static closure [30]. An RTV-type loop option would
require pressure and heating in addition to loop length and is therefore a
separate future model [31].

### 4.6 Open--closed interface and closed-loop policy

The separatrix is topology aligned. A sharp representation is queried from a
declared side; interpolation, derivative stencils, and shock quadrature may
not combine its open and closed states. A shock element intersected by it is
split and each child uses its own upstream authority. A finite-width
representation instead owns a physical signed-distance coordinate and an
explicit layer thickness; it is not produced by blending primitive variables
through an arbitrary angular weight.

Let `n` point from catalog side A to catalog side B; for an open--closed
separatrix, side A is closed and side B is open. The stored surface basis
`(n,t_1,t_2)` is deterministic, right handed, and bound to a stable basis ID.
A geometrical moving surface defines a unique **normal** speed `V_I,n`; a tangential "interface
velocity" is only a surface-parametrization choice and must not enter a
physical jump condition. In one fixed declared coordinate frame define

\[
 w_{n,k}=\mathbf u_k\mathbin{\cdot}\mathbf n-V_{I,n},\qquad
 B_{n,k}=\mathbf B_k\mathbin{\cdot}\mathbf n .
\]

Every **open--closed separatrix** policy below retains the two hard
topology/kinematic gates

\[
 |B_{n,k}|\le\tau_{B,\rm abs}+\tau_{B,\rm rel}|\mathbf B_k|,
 \qquad
 |w_{n,k}|\le\tau_{w,\rm abs}+\tau_{w,M}c_{f,k},
\]

at every qualified point on each side. The second expression is relative to
the moving interface; using `u dot n` without subtracting `V_I,n` is
incorrect for a moving surface. The combined absolute-plus-relative form
remains defined at magnetic nulls and vanishing flow. A relative diagnostic is
typed `inapplicable-near-zero` when its denominator is unresolved, while its
absolute gate remains active.

For a **sharp** interface, the relevant momentum-flux traction is the full
conservative moving-control-surface flux, not only its normal component. For
side `k`,

\[
 \mathbf t_k=
 \rho_k w_{n,k}\mathbf u_k+
 \left(p_k+\frac{|\mathbf B_k|^2}{2\mu_0}\right)\mathbf n
 -\frac{B_{n,k}}{\mu_0}\mathbf B_k,
 \qquad
 \Delta\mathbf t=\mathbf t_B-\mathbf t_A.
\]

This expression is invariant to a tangential reparametrization of the
surface. A formula written with a full relative-velocity vector is legal only
in one explicitly fixed Galilean frame and together with a verified zero
mass-flux jump; it is not the primary interface contract. The provider also
reports and gates

\[
 \Delta j_m=\rho_Bw_{n,B}-\rho_Aw_{n,A}
 \quad [\mathrm{kg\,m^{-2}\,s^{-1}}],
\]

including the signed and absolute values, a dimensionless value with a
declared nonzero normalization floor, its physical-area distribution, the
configured quantile, and the pointwise maximum.

A stationary ideal-MHD tangential discontinuity in the interface frame
satisfies `B_n=0`, `w_n=0`, `Delta j_m=0`, and `Delta t=0` pointwise. Only
under those conditions does the normal component reduce to continuity of
`p+B_t^2/(2*mu_0)`. The implementation must not assume that the two magnetic
traces are identical and then test `[p]=0`; such an assumption would conceal
a tangential magnetic-stress error. It emits the signed normal component,
both signed tangential components in a deterministic surface basis, the
Euclidean norm, absolute-plus-relative normalized values, the area-weighted
distribution, a configured quantile, and the pointwise maximum. A global
signed or absolute mean is diagnostic only and can never pass a gate.

For a **finite-width** interface or empirical open-wind/plasma-sheet
transition, a surface jump is not the correct balance diagnostic. In the
declared frame the provider evaluates the volume momentum residual

\[
 \mathbf R_{\rm mom}=
 \frac{\partial(\rho\mathbf u)}{\partial t}+
 \boldsymbol\nabla\mathbin{\cdot}
 \left[\rho\mathbf u\mathbf u+
 \left(p+\frac{B^2}{2\mu_0}\right)\mathbf I-
 \frac{\mathbf B\mathbf B}{\mu_0}\right]
 -\rho\mathbf g-\mathbf f_{\rm declared}.
\]

`f_declared` contains every explicitly retained nonideal, wave-pressure, and
frame-force term, with signs and units recorded in the manifest; omitting a
term does not make its residual vanish. The layer reports every vector
component, `|R_mom|`, an absolute force-density value in `N/m^3`, a
normalization by the local sum of retained force magnitudes with a registered
floor, pointwise values, volume-weighted quantiles, and the maximum. It also
demonstrates stability across preregistered mesh and physical-thickness
ensembles. A sharp traction jump and a smooth volume residual are different
objects and are never compared as though they were the same metric.

The same representation rule applies to **open--open** structure. A genuinely
sharp boundary between two open states reports the two one-sided states,
`Delta j_m`, and the complete vector traction jump. A resolved fast/slow-wind
transition, streamer/plasma-sheet layer, or any other continuously varying
open-state blend reports the volume momentum residual and its declared force
inventory instead; it must not manufacture one-sided traces merely to obtain a
small jump. D8 covers both open--closed and open--open interfaces and records
the stable interface ID, interface class, state provenance, deterministic
surface basis, topology classes, and representation type before selecting the
applicable diagnostic. The separatrix-specific `B_n≈0` and `w_n≈0` hard gates
are **not** imposed on a generic open--open boundary: an open--open MHD class
boundary may carry normal magnetic flux and mass flux. Its sharp contract is
instead the signed mass-flux jump and complete vector traction jump; its smooth
contract is the resolved volume momentum residual. A future open--open
tangential-discontinuity subclass would need its own explicit catalog tag and
gates rather than inheriting the open--closed assumption implicitly.

`open_closed_interface.policy` selects one of three contracts:

1. **`diagnostic-kinematic`.** This is the default for the independently
   assembled analytic PFSS/SCS field, open wind, and hydrostatic closed
   plasma. The hard `B_n` and relative-normal-flow gates still apply, but the
   traction jump or volume residual is reported rather than forced to zero.
   The state is a kinematic background, not an MHD equilibrium. It is valid
   for analytic verification and declared sensitivities, but cannot be called
   a stationary tangential discontinuity or used as an event-grade equilibrium.
2. **`bounded-approximation`.** This is the event-grade stand-alone option
   when no self-consistent global MHD state is available. In addition to the
   hard gates, every sharp-interface traction residual or every smooth-layer
   volume residual must stay within preregistered absolute and relative
   pointwise bounds. The manifest records observational input uncertainty,
   interface-location or layer-thickness uncertainty, the mesh/thickness
   convergence sequence, the selected quantile, and the maximum. Passing
   these bounds establishes only a quantified approximation, not exact force
   balance.
3. **`stationary-td-equilibrium`.** This strict sharp-interface policy is
   available only to a solver-produced or imported state whose provenance
   declares the equilibrium construction. It requires the full vector
   traction jump, mass-flux jump, `B_n`, and interface-relative `w_n` to pass
   at every point.
   An independently normalized analytic open/closed pair cannot select this
   policy merely because a post hoc fit makes one scalar residual small.

The closed provider may use prescribed base normalization, a single global
pressure-scale sensitivity, or a loop/footpoint-resolved constrained asset.
Those operations explore sensitivity; they do not solve the multidimensional
force-balance problem. In particular, copying an open-side temperature to the
closed side cannot enforce pressure or traction continuity because density,
composition, magnetic stress, and tangential momentum also enter the balance.
Opposite local residuals can cancel in an area mean, and a fit that minimizes a
quantile can leave an unacceptable maximum. Therefore all policies preserve
pointwise records, quantiles, and maxima, and event gates are evaluated on the
pointwise bounds. The null/cusp mask and its area are reported; no quantity is
divided by `B^2/(2*mu_0)` inside that mask. A loop with incompatible
two-footpoint hydrostatic data is rejected rather than averaged.

A finite-width layer must reconstruct an EOS-consistent state, preserve or
re-solve the open-tube invariants, and pass registered resolution and
thickness studies. A universal angular blend of density and temperature is
not a physical default.

One-way outward WKB turbulence is invalid on a line with two solar
footpoints. Closed-field scattering must explicitly select a bidirectional
prescribed model, a direct mean-free-path closure, or a future two-footpoint
reflection/cascade model. No automatic WKB-to-prescribed substitution is
allowed. With parallel-only transport, particles remain on the loop and may
precipitate into either absorbing solar footpoint; cross-field escape requires
an explicitly enabled and interface-aware operator.

No such open/closed-separatrix transmission operator belongs to the first
production profile. Therefore a nonzero perpendicular coefficient is rejected
whenever a reachable stochastic domain intersects the separatrix. Analytic
verification may use perpendicular diffusion only after proving geometrically
that the reachable domain remains on one smooth topology side. A failed proof
is a configuration error, not an absorbing or reflecting fallback.

### 4.7 PFSS limitations and numerical qualification

PFSS neglects active-region currents, magnetic free energy, and dynamical
opening by the CME. Its purpose here is to provide a reproducible large-scale
coronal topology and flux-tube expansion, not to model CME initiation. The
limitation is strongest close to a sheared active region and must be discussed
when event-specific results are interpreted [2--4].

High-degree spherical-harmonic evaluation can also suffer from map-grid
aliasing, polar errors, and ringing. Direct harmonic evaluation remains the
authority; any cached spherical grid is only an optimization and must reproduce
the direct field, derivatives, topology, and flux within registered error.
Magnetogram qualification reports reconstruction error, flux changes, and
topology sensitivity versus `L_max`, grid/quadrature choice, and filter
strength at both `1.01` and `1.05 R_sun`.

---

## 5. Transonic solar wind on an open flux tube

### 5.1 Flux-tube area

Magnetic-flux conservation on the completed PFSS--SCS--Parker line gives

\[
  A(s)|B(s)|=\Phi_B=\text{constant},
\]

where `s` is distance along the open field line. The production wind solver
uses the complete composite field; a PFSS-only tube is a kernel-verification
case. It must not assume `A proportional r^2` below `R_scs`, and a
longitude-structured Parker extension must also include the mapping Jacobian
described in Section 6.

The scalar wind speed is defined without reference to magnetic polarity. Let
`t_hat=dx/ds` point geometrically away from the Sun and let

\[
 \mathbf u'=\mathbf u-\boldsymbol\Omega_F\times\mathbf x
             =u_s\hat{\mathbf t},\qquad u_s>0,
\]

where `u` is the inertial velocity and `u'` is the velocity in the selected
rigidly rotating frame. Because
`B=sigma_B |B| t_hat`, `u'` is parallel or antiparallel to signed `B` according
to sector but is always outward. Define

\[
 u_r=\mathbf u\cdot\hat{\mathbf e}_r
     =u_s\hat{\mathbf t}\cdot\hat{\mathbf e}_r.
\]

Thus `u_s` is the speed in the one-dimensional nozzle equations, while WSA-like
targets and heliospheric observations refer explicitly to inertial radial
speed `u_r`. In a tube-labeled differential-rotation option, `Omega_F` is
constant along one tube but may vary between source latitudes. That branch is
a quasi-steady ensemble, not one globally stationary rotating-frame MHD
solution, and its validity requires a reported evolution-timescale check.

### 5.2 Conservation laws

For a steady, inviscid, single-fluid baseline without wave pressure, mass
conservation is

\[
  \rho u_s A=\dot m_{\rm tube}=\text{constant}.
\]

In the rigidly rotating tube frame, define

\[
 \Phi_{\rm eff}(\mathbf x)=-\frac{GM_\odot}{r}
 -\frac12|\boldsymbol\Omega_F\times\mathbf x|^2.
\]

The field-aligned momentum equation is

\[
  u_s\frac{du_s}{ds}
  =-\frac{1}{\rho}\frac{dp}{ds}
   -\frac{d\Phi_{\rm eff}}{ds}.
\]

The Coriolis force does no work along `u'`. The centrifugal potential is kept;
dropping it is a named small-rotation verification approximation. This
one-dimensional equation enforces along-tube dynamics but not transverse force
balance or the magnetic torque needed to derive the adopted azimuthal flow.

The effective polytropic wind closure is

\[
  p=K_w\rho^{\gamma_w},
  \qquad a_w^2=\frac{\partial p}{\partial\rho}
             =\gamma_w\frac{p}{\rho},
\]

where `a_w` is the effective wind sound speed. Combining these equations gives the
nozzle form

\[
 \left(u_s-\frac{a_w^2}{u_s}\right)\frac{du_s}{ds}
 =a_w^2\frac{d\ln A}{ds}
  -\frac{d\Phi_{\rm eff}}{ds}.
\]

At a smooth critical point,

\[
 u_{s,c}=a_{w,c},
 \qquad
 a_{w,c}^2\left.\frac{d\ln A}{ds}\right|_c
 =\left.\frac{d\Phi_{\rm eff}}{ds}\right|_c.
\]

The physical solution passes continuously through the appropriate critical
point and becomes supersonic. Rapidly expanding coronal-hole tubes can admit
multiple mathematical critical points; the admissible branch must be selected
using the global transonic topology, not a nearest-root rule [1, 5].

For the polytropic model with `gamma_w>1`, the Bernoulli invariant is

\[
  \frac{u_s^2}{2}+\frac{a_w^2}{\gamma_w-1}
  +\Phi_{\rm eff}=\mathcal E_w,
\]

provided no explicit heating, wave work, or thermal conduction is included.

The isothermal branch is selected explicitly, not by an exact floating-point
comparison inside the solver. With constant `a_w`,

\[
 p=a_w^2\rho,
 \qquad
 \frac{u_s^2}{2}+a_w^2\ln\left(\frac{\rho}{\rho_*}\right)
 +\Phi_{\rm eff}=\mathcal E_{\rm iso}.
\]

`rho_*` is a fixed positive reference density used only to make the logarithm
dimensionless; changing it shifts the isothermal Bernoulli constant by a
constant and cannot change the
wind solution.

For this no-heating accelerating-wind closure, the production polytropic range
is restricted to `1 < gamma_w < 3/2`; the spherical Parker problem identifies
`3/2` as the critical value. The restriction is a declared baseline-model
condition, not a theorem for every possible non-spherical heated wind. Every
accepted input must still pass the full variable-area transonic existence and
branch test [1, 5, 29].

### 5.3 Boundary-value formulation

The transonic wind is an eigenvalue problem. The input may prescribe a base
temperature/composition and one independent density or mass-flux constraint,
but it must not independently prescribe an arbitrary base speed, base density,
and 1-AU mass flux if those values overdetermine the transonic solution.

The open-wind base is one explicit sphere `r=R_w0`, with
`R_sun <= R_w0 <= R_in`; schema 5 calls it
`solar_wind.base_reference_radius_m`. The same `R_w0` is used for open-wind
temperature, density/mass loading, composition, any open turbulence datum
explicitly declared base-owned, and the Bernoulli boundary condition. The
closed-plasma and non-base turbulence records have their own named reference
radii and are not silently forced onto `R_w0`. The global transonic solution is
integrated both inward and outward as needed, so
choosing `R_w0>R_sun` does not leave the guard shell below `R_w0` undefined.
Changing `R_w0` while retaining the same numerical base values changes the
physical model and therefore the fingerprint.

The solver should use the following procedure:

1. trace an open composite tube through PFSS and SCS to `R_scs` and, when
   required by a target-speed condition, through the outer Parker domain;
2. calculate a differentiable `A(s)` separately inside every smooth magnetic
   region and identify every typed interface;
3. locate every critical-point candidate;
4. determine the global transonic branch by regularity and outer-boundary
   consistency;
5. integrate inward and outward from the critical point using a regularized
   derivative and apply the interface conditions below; and
6. determine the density normalization from the selected base density or
   resolution-independent outer mass-flux measure.

At an idealized PFSS/SCS kink, `|B|` and hence `A` can jump even though `B_r`
is continuous. The smooth nozzle ODE is not integrated through that jump. Its
one-sided states preserve the mass loading per magnetic flux,

\[
 \eta_m=\frac{\rho u_s}{|B|}
       =\frac{\dot m_{\rm tube}}{\Phi_B},
\]

the selected wind entropy/isothermal closure, and the appropriate Bernoulli
invariant. Continuity of `B_r` then gives continuity of normal mass flux
`rho u_n=eta_m B_r/sigma_B`, where at this spherical interface
`u_n=u dot e_r=u_r`. The prescribed surface-current layer may carry a
normal/tangential stress residual; that residual is reported and is a model
limitation, not silently forced to zero. The production recommendation is a
finite, divergence-preserving transition whose width is resolved by the wind,
AMR, and particle step. The zero-thickness kink remains a magnetic/kernel
verification option and cannot be crossed by a transport mover without a
separately verified conservative interface operator.

Mass flux, Bernoulli invariant, positivity, and monotonic passage through the
critical point are mandatory diagnostics. Interface mass-loading, energy,
stress-residual, and transition-width diagnostics are mandatory as well.

### 5.4 Derived plasma quantities

For the selected composition and thermodynamic plasma EOS,

\[
  v_A=\frac{|\mathbf B|}{\sqrt{\mu_0\rho}},
  \qquad
  c_s=\sqrt{\frac{\gamma_{\rm ad}p}{\rho}}.
\]

For propagation at angle `theta_Bn` to the upstream magnetic field, the fast
speed is

\[
 c_f^2=\frac12\left[
 v_A^2+c_s^2+
 \sqrt{(v_A^2+c_s^2)^2
       -4v_A^2c_s^2\cos^2\theta_{Bn}}
 \right].
\]

These values must be computed from the same mass density, pressure, and field
used by the shock solver. Electron density cannot be substituted for mass
density in `v_A`.

For the first release, `gamma_ad=5/3` represents an isotropic,
nonrelativistic, ideal single-fluid plasma. A later multispecies EOS may use

\[
 c_s^2=\frac{\sum_s\gamma_{{\rm ad},s}p_s}{\rho},
\]

where `gamma_ad,s` is the rapid adiabatic response of species `s`; such a
model must replace the scalar-pressure shock closure consistently.

### 5.5 Uniform, target-speed, and empirical-kinematic open-tube closures

The analytic verification mode uses one declared base-temperature closure on
all open tubes. This mode is reproducible and useful for isolating geometry,
but it cannot by itself reproduce the full observed fast/slow-wind contrast.

A `tube-from-target-speed` polytropic option may instead define

\[
 u_{t,i}=F_{\rm WSA}(f_{s,i},\theta_{b,i};\mathbf c),
\]

where `f_s` is the tube expansion, `theta_b` is distance to the selected
open-field boundary, and the exact relation and coefficients come from a
versioned, checksummed asset [39, 40]. That asset declares the topology
definition, expansion-factor radii, coronal-hole-boundary algorithm, target
speed definition, calibration interval, and coefficient units. For a
finite-radius target, each tube solves

\[
 F_i(T_{0,i})=
 u_{r,i}(r_t;T_{0,i},\gamma_w,A_i)-u_{t,i}=0,
\]

where `u_r`, not `u_s`, is the inertial radial speed. The asymptotic option
instead solves

\[
 \lim_{r\to\infty}\left[u_r(r)-u_t\right]=0
\]

on the asymptotic branch and requires `target_speed_radius_m=0`; a finite-radius option requires a positive
radius equal to that encoded by the coefficient asset. Both use a bracketed
inversion of the complete transonic boundary-value solution.
A closed-form Bernoulli temperature may seed the solve but is not the
accepted answer because it omits the finite-radius enthalpy and gravity,
base eigenvelocity, and its dependence on temperature.

This inversion is a useful controlled sensitivity branch, but it is not
automatically event grade. A no-heating polytrope has too few physical degrees
of freedom to fit, independently, a tube's terminal speed, low-coronal density
and temperature, and heliospheric mass flux. In the spherical limit the same
temperature change that produces a fast target speed can generate orders-of-
magnitude changes in the critical-point density and mass flux; variable area
and rotation change the numbers but not the need for observational gates
[29, 53]. Every target-speed solution therefore publishes D7 and may be used
only as a schema-5 sensitivity branch; it cannot be labeled event-nominal or
used as the production-qualified background. A sensitivity campaign may attach
an observational D7 asset and must then report its realized density,
temperature, and mass-flux mismatches, but passing them does not promote the
closure. Combining `tube-from-target-speed` with one uniform base density is
forbidden, and replacing it with a uniform mass-loading constraint does not
remove the requirement to emit D7.

Mass-flux normalization is an independent, resolution-invariant closure. It
may prescribe the mass loading per unit magnetic flux,

\[
 \eta_{m,i}=\frac{\rho_i u_{s,i}}{|B_i|}
 \quad [\mathrm{kg\,s^{-1}\,Wb^{-1}}],
\]

or the radial mass-flux density at a declared radius,

\[
 F_{m,i}(r_m)=\rho_i(r_m)u_{r,i}(r_m)
 \quad [\mathrm{kg\,m^{-2}\,s^{-1}}].
\]

Alternatively, a declared base density fixes the normalization. An absolute
`kg/s` assigned independently to every numerical trace is forbidden because
it changes under tube refinement. A versioned empirical density--speed asset
may supply either `eta_m` or `F_m` and must declare its represented magnetic or
solid-angle measure. The same datum may not independently normalize both the
base and outer states. The structured solution is an ensemble of
one-dimensional flux tubes, not a transverse 3-D MHD equilibrium; it does not
form stream-interaction regions or CIR shocks.

The recommended first event branch is `empirical-kinematic`, but an
event-nominal instance is **not** a speed-only density inference. A speed at
one radius does not determine mass loading: without a co-located density and
magnetic field, infinitely many values of `eta_m` produce the same speed. The
expression `rho=eta_m|B|/u_s` is also singular as `u_s->0`. The deficiency is
physical underdetermination (and a zero-speed singularity), not a numerical
conditioning problem that can be repaired by choosing a different
interpolator.

The event-nominal closure therefore has two observational zones on every
qualified open tube:

1. an **inner density zone**, constrained by an electron-density or
   composition-resolved mass-density reconstruction; and
2. an **outer velocity zone**, constrained by an independently processed
   velocity profile.

These are independent channel records, not columns whose meaning is imposed by
one file-wide velocity or abscissa selector. A checksummed channel manifest identifies
the inner density, outer velocity, and every species-temperature record
separately. Each channel owns its variable, units, velocity component/reference
frame and consumer selector where applicable, spatial frame, abscissa, line/
tube IDs, support intervals, interpolation certificate, uncertainty/covariance
identity, and data-use role. Density may, for example, be radial while velocity
is supplied on oriented arc length, provided both are evaluated on the same
fingerprinted tube and genuinely overlap. No channel inherits another
channel's coordinate or support, and the provider never relabels a radial
record as field aligned, radius as arc length, or one species temperature as
another.

The zones overlap over the explicit radial interval `[r_a,r_b]`, where
`R_w0<=r_a<r_b` and both products have support. Every use of this radial blend
requires `r(s)` to be strictly outward monotone, with its registered derivative
margin satisfied throughout the overlap. A nonmonotone line is split at each
certified turning point into maximal monotone segments; no blend, interpolant,
source, observer, or export request may span a turning point, and repeated
values of `r` are never silently assigned the same branch. Let `rho_in(s)` be
the positive inner-zone mass density after the composition conversion below. Let
`u_s,out(s)` be the positive outer velocity converted into the outward-
geometric-tangent, corotating quantity according to its declared component and
reference-frame pair. Exactly one
mass-per-flux authority `eta_m>0` is fixed before constructing the overlap.
The outer-zone density implied by that authority is

\[
 \rho_{\rm out}(s)=\frac{\eta_m|\mathbf B(s)|}{u_{s,{\rm out}}(s)}.
\]

With `xi=(r-r_a)/(r_b-r_a)`, a positive reference density `rho_*` recorded in
the channel manifest, and the `C2` endpoint-flat blend

\[
 w(\xi)=10\xi^3-15\xi^4+6\xi^5,
\]

the authoritative density is

\[
 \ln\!\left(\frac{\rho(s)}{\rho_*}\right)=
 \begin{cases}
   \ln[\rho_{\rm in}(s)/\rho_*], & r\le r_a,\\
   [1-w(\xi)]\ln[\rho_{\rm in}(s)/\rho_*]
     +w(\xi)\ln[\rho_{\rm out}(s)/\rho_*], & r_a<r<r_b,\\
   \ln[\rho_{\rm out}(s)/\rho_*], & r\ge r_b,
 \end{cases}
 \qquad
 u_s(s)=\frac{\eta_m|\mathbf B(s)|}{\rho(s)}.
\]

Thus positivity is preserved, density and its first two radial derivatives
join smoothly when the one-sided input interpolants are `C2`, and the final
velocity is derived **everywhere** from the same continuity authority. It
equals the prescribed outer velocity beyond `r_b`; inside the overlap the
velocity adjusts continuously to the density product rather than violating
mass conservation. The overlap is a compatibility test, not a hidden fit.
The provider reports

\[
 d(s)=\Delta_{\ln\rho}(s)=
 \ln\!\left[\frac{\rho_{\rm in}(s)}{\rho_{\rm out}(s)}\right]
\]

pointwise, by preregistered weighted quantile, and by maximum. It also reports
an explicit covariance-aware statistic. On the preregistered overlap
evaluation grid, collect the residuals into `d` and transform the declared
joint density/velocity/composition/magnetic-field/mass-loading covariance into
log-residual space with the complete Jacobian, `C_d=J C_x J^T`. The provider
verifies that `C_d` is symmetric positive semidefinite within a registered
eigenvalue tolerance, records its numerical rank `r_C`, and evaluates

\[
 \chi_d^2=\mathbf d^{\mathsf T}C_d^+\mathbf d,
 \qquad \chi_{d,\rm red}^2=\chi_d^2/r_C,
\]

where `C_d^+` is the Moore--Penrose inverse formed with the same preregistered
rank tolerance. A component of `d` in the covariance null space beyond its
absolute numerical tolerance is an inconsistent asset, not a direction of
zero penalty; `r_C=0` is typed covariance-inapplicable and cannot qualify an
event member. No diagonal-only approximation or regularization may silently discard the
declared cross-channel correlations. An event-nominal member fails
when either the absolute-log mismatch or covariance-normalized mismatch
exceeds its preregistered bound; it never rescales one product to manufacture
agreement.

`eta_m` itself must come from exactly one versioned normalization record. A
record may provide `eta_m` directly; a co-located positive `rho`, `u`, and
mapped `B`; or a radial mass-flux density `F_m=rho u_r` together with the
mapped nonzero `B_r` at the same position and epoch. The latter two give

\[
 \eta_m=\frac{\rho u_s}{|\mathbf B|}
       =\frac{\rho u_r}{|B_r|}
       =\frac{F_m}{|B_r|}.
\]

A lone speed, or a mass-flux product without its co-located mapped field,
cannot define `eta_m` and is rejected. Each normalization record states its
represented magnetic-flux or solid-angle measure so that refinement cannot
change the physical loading.

Electron density is converted to mass density only through the declared
composition and charge state. With ion abundance `f_j=n_j/n_p` and charge
number `Z_j`, quasineutrality gives

\[
 n_p=\frac{n_e}{\sum_j Z_j f_j},\qquad
 \rho=n_p\sum_j m_j f_j+\epsilon_e m_e n_e,
\]

where `epsilon_e` is exactly the `plasma_eos.electron_mass_in_density`
choice. The asset records abundance uncertainties and their covariance with
the density inversion; an alpha abundance is not silently assumed from the
`proton-electron` branch. Pressure follows from the same species abundances
and the declared species-temperature products.

Velocity component, velocity reference frame, and profile abscissa are
independent discriminated fields on each channel record. The pair
`velocity_component=radial, velocity_reference_frame=inertial` means the
channel contains inertial `u_r`. Position, tangent, and velocity are first
transformed to the channel's declared spatial frame and trace epoch. With
`q=t_hat dot e_r`, the provider must then use

\[
 u_s=\frac{u_r}{\hat{\mathbf t}\cdot\hat{\mathbf e}_r}
\]

over its complete support with the **same channel's**
`t_hat dot e_r >= minimum_radial_projection > 0`. A tangent with zero or
negative radial projection, an inward-turning segment, or a sub-margin value
is a typed geometry/profile mismatch; the provider never divides by zero,
changes the tangent sign, borrows another channel's tolerance, or silently
reinterprets radial data as field-aligned speed.

The pair `field-aligned,corotating` contains the rotating-frame quantity

\[
 u_s=(\mathbf u_{\rm inertial}-\boldsymbol\Omega_F\times\mathbf x)
       \mathbin{\cdot}\hat{\mathbf t}
\]

directly and does not apply the radial projection. The pair
`field-aligned,inertial` instead stores
`u_inertial dot t_hat`; the provider subtracts
`(Omega_F cross x) dot t_hat` exactly once at the registered epoch. The
combination `radial,corotating` is noncanonical and unsupported by schema 5
(rigid heliocentric rotation has no radial component), so it is rejected
rather than given a second meaning. Both field-aligned
pairs set the channel-local radial-projection guard to its inactive zero value
and are legal only for separately supplied products bound to stable field-line
IDs, trace epoch, spatial frame, rotation authority, and the immutable
background fingerprint on which they were reconstructed. Changing a label is
not a frame conversion, and a radial file can never activate either branch.

Here field aligned always means the outward geometric tangent `t_hat`; magnetic
polarity remains a separate `b_hat`/sector datum. Each outer-velocity channel
also owns a selector over either stable mutually exclusive topology classes or
an explicit set of stable line/tube IDs. A broad topology selector may subtract
explicit stable IDs so an exceptional field-aligned subset does not overlap the
ordinary radial population. Before evaluation, every selector is expanded
against the frozen line/tube catalogue. Over the requested physical support,
the expanded domains must be pairwise disjoint and every required source,
observer, and export consumer must match exactly one channel. Multiple matches
return `AmbiguousKinematicRoute`; zero matches follow the declared required-
consumer coverage policy. File order, support order, component type, and
priority are not routing rules, and failure of the selected channel's support,
frame, projection, or interpolation certificate never falls back to another
channel. This permits ordinary tubes to use a radial product while separately
identified nonradial/separatrix-adjacent tubes use a field-aligned product
without requiring one contradictory global projection setting.

Likewise, `abscissa=heliocentric-radius` evaluates the asset at `r`, whereas
`abscissa=oriented-field-line-arclength` evaluates it at the bundle's
oriented `s`, whose origin is `s=0` at the registered inner endpoint of that
specific trace. Arc-length data require that origin, the trace epoch and frame,
the same line-ID catalogue, and the exact background fingerprint; a
radius-tabulated product is never reinterpreted as arc length.

The required interpolation is the exact
`quintic-hermite-c2-certified` construction. For every positive channel the
asset supplies, at every node, the value and the first and second derivatives
with respect to its own abscissa. For a physical value `y>0` and its recorded
positive reference `y_*`, the provider forms
`z=ln(y/y_*)`, `z'=y'/y`, and
`z''=y''/y-(y'/y)^2`, then interpolates `z` with the unique quintic Hermite
polynomial matching `z,z',z''` at both endpoints. It verifies bit-consistent
shared node derivatives and hence `C2` continuity, finds every real stationary
point of each polynomial segment, and certifies finite positive reconstruction
and no overshoot outside the endpoint range. Failure of any certificate
rejects the segment; derivatives are never invented online by a different
application. Extrapolation across a
data edge, topology event, magnetic interface, turning point, or line-ID change
is forbidden.

Every density, velocity, temperature, composition, and normalization record
declares instrument/reconstruction source, processing version, epoch,
coordinate frame, units, topology class, stable tube/line IDs, radial or
arc-length support, checksum, data-use role (`construction`, `qualification`,
or `withheld-validation`), uncertainty, and covariance/ensemble identity.
Construction and qualification roles cannot reuse the same observations.
Coverage is evaluated by physical measure, not by counting numerical traces:
the manifest reports covered open magnetic flux and open area, potential
source incident-number/energy flux, every finite observer-footprint exposure,
and every requested export interval. In event-nominal mode, any uncovered
required source support, observer footprint, or export interval is fatal.
Diagnostic masking is available only to a sensitivity run and leaves an
explicit uncovered-measure census. The census preserves every originally
requested stable line/tube ID. A rejected ID remains in the manifest with its
original identity, requested support, physical measure when defined, and all
typed rejection reasons; accepted IDs are never renumbered to hide a rejection.

This branch satisfies continuity by construction but does **not** claim to
solve the field-aligned momentum equation. It reports the signed acceleration
residual

\[
 {\cal R}_{\rm mom}=
 u_s\frac{du_s}{ds}+\frac{1}{\rho}\frac{dp}{ds}
 +\frac{d\Phi_{\rm eff}}{ds},
\]

its absolute value, and a dimensionless value normalized by the sum of the
absolute retained accelerations plus a registered nonzero reference floor.
Production requires the residual distribution and maximum to remain below
declared tolerances over every source patch, active mover cell, observer
support, and exported line. Passing this gate makes the approximation measured
and reproducible; it does not turn it into a momentum solution. A future
heated or wave-driven branch is the physical successor.

The existing positivity, continuity, D7 observational comparisons, and
pointwise momentum-residual gates do not by themselves exclude an unphysical
speed spike or excessive deceleration between comparison radii. A future
event-grade qualification may therefore add a **continuous wind-plausibility
envelope**, but that envelope is not an active schema-5 input and is not a
universal solar-wind speed cap. It must be a versioned, checksummed,
event-specific observational asset that states, independently for every
channel:

- whether it constrains inertial radial speed `u_r`, corotating
  field-aligned speed `u_s`, or (for the quasi-steady schema-5 wind) the
  field-aligned advective acceleration `u_s*du_s/ds`; a future time-dependent
  asset instead constrains the consistently projected material derivative
  `b_hat dot [partial_t u+(u dot grad)u]` and may not omit the time derivative;
- coordinate frame, epoch, topology/line selector, physical support,
  uncertainty or covariance, interpolation certificate, and data-use role;
- whether monotonic acceleration is observationally justified on that support
  or whether bounded deceleration is allowed; and
- the continuous acceptance functional, including how extrema between asset
  nodes and interpolation uncertainty are bounded.

Qualification evaluates the authoritative continuous wind representation,
including every certified interior extremum, rather than sampling only the
asset nodes. A limit on `u_r` is never applied to `u_s` by relabeling, and an
acceleration envelope is never replaced by a speed envelope. Construction and
qualification observations remain disjoint under the same data-role rules as
D7. Until such an asset and its tests exist, schema 5 reports the current D7
and residual diagnostics and must not claim continuous observational
plausibility between the registered comparison points.

When the target speed varies with source longitude, the wind and longitude
map in Section 6 are coupled. A qualified polytropic preparation iterates the transonic
tube solution, longitude map, mapping Jacobian, composite `B`, and tube area
until velocity, mass flux, momentum residual, and mapping all meet tolerance.
Failure to converge or formation of a mapping fold is fatal. A one-pass
ballistic mapping is not licensed merely because the wind profiles are
empirical: the kinematic branch uses the same full longitude Jacobian and must
also fail before a fold. Its measured momentum residual remains explicit.

### 5.6 More advanced wind option

A later wave-driven option may add wave pressure and heating to the momentum
and energy equations. That option would be closer to a ZEPHYR-class open
flux-tube model [6], but it must be a distinct named closure. A polytropic
wind must not be described as wave heated merely because turbulence is also
initialized for particle scattering. Heating below and above the sonic point
affects mass flux and terminal speed differently [53], so the future closure
must solve a stated energy equation rather than tune one effective temperature.

---

## 6. Solenoidal PFSS--current-sheet--Parker magnetic field

### 6.1 Magnetic regions and radii

The production field contains three analytic authorities plus a resolved join:

1. signed PFSS for `R_sun <= r < R_i`;
2. a divergence-free PFSS/SCS transition for
   `R_i <= r <= R_i+w_tr`;
3. a finite-shell Schatten-type current-sheet solution for
   `R_i+w_tr < r < R_scs`; and
4. a conservative rotating-footpoint Parker extension for `r >= R_scs`, with
   its value at equality taken from the single shared radial interface trace.

PFSS itself is solved to the radial boundary `R_b`. A direct **PFSS/SCS
interface** uses `R_i=R_b`. An overlap-minimized PFSS/SCS option may use
`R_i<R_b`: PFSS remains solved to `R_b`, but its field at `R_i` supplies the
SCS inner boundary. `R_i` is an interface radius, not a CSSS cusp unless the
distinct CSSS equations are implemented. The selected radii, open-flux change,
interface kink statistics, and topology hashes are part of the physics
fingerprint [34, 35, 51].

A separate no-SCS analytic verification branch couples the signed radial PFSS
field directly to the Parker map at `R_b`. In that branch `R_i` and `R_scs` are
derived as `R_b`, all SCS-only input fields are inactive, no SCS sector map is
constructed, and current-sheet transport is `not-applicable`. Because the
signed PFSS radial field vanishes on its neutral line, this branch is not the
event-grade heliospheric field and cannot be selected by a production deck.

The `direct-pfss-scs` choice with `R_i=R_b` is likewise a **zero-width sharp
verification interface** and requires `w_tr=0`; it cannot satisfy the overlap
condition for a positive transition. Production requires `R_i<R_b` and

\[
 0<w_{\rm tr}\le\min(R_b,R_{\rm scs})-R_i,
\]

with the complete transition resolved and qualified as specified below.

The finite-shell construction below is related to, but not identical to, the
original Schatten solution whose potential outer boundary is at infinity.
The distinction is required in output and publications.

### 6.2 Unsigned finite-shell SCS field

Define the mathematical unsigned inner target

\[
 q_0(\theta,\phi)=
 |B_r^{\rm PFSS}(R_i,\theta,\phi)|.
\]

Represent it with an independently resolved, nonnegative-constrained SCS
spectral fit `q_L`, and define an unsigned potential field

\[
 \widetilde{\mathbf B}=-\nabla\widetilde\Psi,
 \qquad \nabla^2\widetilde\Psi=0,
 \qquad R_i\le r\le R_{\rm scs},
\]

with boundary conditions

\[
 \widetilde B_r(R_i,\theta,\phi)=q_L(\theta,\phi),
\]

\[
 \widetilde B_\theta(R_{\rm scs})=
 \widetilde B_\phi(R_{\rm scs})=0.
\]

Using

\[
 \widetilde\Psi=
 \sum_{\ell=0}^{L_{\rm scs}}\sum_{m=-\ell}^{\ell}
 \left(c_{\ell m}r^\ell+d_{\ell m}r^{-\ell-1}\right)Y_{\ell m},
\]

the `ell=0` mode is mandatory because the unsigned inner boundary has nonzero
net flux. If `h_lm` are the coefficients of the accepted constrained
representation,

\[
 q_L(\theta,\phi)=\sum_{\ell m}h_{\ell m}Y_{\ell m},
\]

then, for every `ell>=1`, the finite radial outer boundary gives

\[
 d_{\ell m}=-c_{\ell m}R_{\rm scs}^{2\ell+1},
\]

\[
 c_{\ell m}=-\frac{h_{\ell m}}
 {\ell R_i^{\ell-1}+
  (\ell+1)R_{\rm scs}^{2\ell+1}R_i^{-\ell-2}}.
\]

For `ell=0`, the tangential outer-boundary condition is automatic and does
not determine the additive constant in the scalar potential. The
implementation fixes the harmless monopole gauge by imposing
`Psi_00(R_scs)=0`; only with that explicit gauge does
`d_00=-c_00 R_scs` hold. The inner unsigned flux then determines
`d_00=h_00 R_i^2`. Magnetic fields and every diagnostic are gauge invariant.

The shell thickness is a physical/numerical model choice, not a cosmetic
coupling parameter. For an unsigned inner-boundary harmonic of degree
`ell>=1`, the amplitude at the radial outer boundary relative to the monopole,
divided by the same ratio at `R_i`, is

\[
 G_\ell(x)=
 \frac{(2\ell+1)x^{\ell+1}}
      {\ell+(\ell+1)x^{2\ell+1}},
 \qquad x=\frac{R_{\rm scs}}{R_i}.
\]

Thus a geometrically thin shell can be almost inert. For example,
`R_i=2.3 R_sun`, `R_scs=2.5 R_sun` gives `G_2=0.9800`; it removes only about
two percent of the dominant quadrupolar component of a dipole-like
`|cos(theta)|` pattern. This number is illustrative rather than universal:
the realized attenuation depends on `R_scs/R_i` and the complete accepted
harmonic spectrum. Preparation reports every retained `G_ell`, the input and
output non-monopole powers, and the dominant-mode attenuation before any
Parker winding is applied.

The harmonic normalization and quadrature convention are compatible with PFSS,
but `L_scs` and the positivity-preserving fit controls are independent
fingerprinted inputs.
The absolute-value cusp generates power above the PFSS truncation and a naive
finite reconstruction can ring negative. Therefore an adaptive angular
verification must establish `q_L>=0`, bounded unsigned-flux error, converged
neutral-line displacement, no unintended nulls, and successful tracing of
every sampled unsigned line from `R_i` to `R_scs`. The fit must conserve the
unsigned monopole flux and reduce the one-sided normal-field mismatch under
angular refinement. A failing reconstruction is rejected; negative values are
never clipped cell by cell. The mathematical sharp-interface limit uses `q_0`
exactly, while a finite representation declares and reports its maximum and
flux-weighted `q_L-q_0` error.

Production additionally applies two independent radialization gates. First,
let

\[
 P_{\rm nm}^{\rm in}=\sum_{\ell\ge1,m}|h_{\ell m}|^2,
 \qquad
 P_{\rm nm}^{\rm out}=\sum_{\ell\ge1,m}|h_{\ell m}|^2G_\ell^2,
 \qquad
 P_{\rm zon}^{\rm out}=\sum_{\ell\ge1}|h_{\ell 0}|^2G_\ell^2,
\]

and define the **absolute outer zonal non-monopole fraction**

\[
 F_{\rm zon}^{\rm out}=
 \frac{P_{\rm zon}^{\rm out}}
 {|h_{00}|^2+P_{\rm zon}^{\rm out}}.
\]

The denominator must be positive. `F_zon_out` is the quantity constrained by
`maximum_outer_zonal_nonmonopole_power_fraction`; it measures structure in the
longitude mean and allows both an already-flat input and genuine longitudinal
structure to avoid a false latitude penalty. The all-mode `P_nm_out`, its
absolute fraction, individual `G_ell`, and attenuation ratio
`P_nm_out/P_nm_in` are reported diagnostics. The ratio is defined as zero when
`P_nm_in=0`; none of these all-mode quantities is itself an event-acceptance
gate.

Second, the completed exterior solution is sampled at every explicitly listed
diagnostic radius. Every listed radius must be at or outside `R_scs` and
inside the exact spherical domain. The input list must contain 10 and 20 solar
radii whenever each is both exterior to `R_scs` and inside the domain; the
parser never appends an un-fingerprinted radius.

For `X=r^2|B_r|`, longitude samples inside the common HCS/transition clearance
mask are removed. In latitude bin `j`, form the quadrature-weighted longitude
mean `Xbar_j` from the remaining samples and require the configured minimum
unmasked longitude fraction. With `W_j` equal to the represented unmasked
solid angle and `Xbar=sum_j W_j Xbar_j/sum_j W_j`, define

\[
 \epsilon_{\rm lat,rms}=
 \left[\frac{\sum_jW_j(\overline X_j/\overline X-1)^2}
 {\sum_jW_j}\right]^{1/2}.
\]

The percentile metric is the area-weighted `P95(Xbar_j)/P05(Xbar_j)`, not a
percentile over all longitude pixels. A nonpositive/nonfinite mean or P05,
insufficient longitude coverage in any accepted bin, or nonfinite metric is a
typed gate failure. Sector-resolved longitude structure is still reported as
a diagnostic but does not contaminate the latitude-flatness measure. The RMS
and percentile gates are observation/configuration inputs, not hidden constants.
The outer rotating-footpoint map can redistribute longitude through `J_phi`
but cannot be credited with SCS latitude flattening. Failure of either gate
invalidates the event configuration; moving `R_scs`, changing `R_i`, or
changing SCS resolution requires a new fingerprinted preparation.

### 6.3 Polarity restoration and the HCS

For a point not on the current sheet, trace `B_tilde` to its footpoint
`x_i` on `R_i` and define

\[
 \sigma_B(\mathbf x)=
 \operatorname{sgn}\left[B_r^{\rm PFSS}(R_i,\mathbf x_i)\right],
 \qquad
 \mathbf B(\mathbf x)=\sigma_B(\mathbf x)\widetilde{\mathbf B}(\mathbf x).
\]

`sigma_B` is a discrete sector identifier, not an interpolated scalar. The HCS
is the union of unsigned field lines launched from the PFSS neutral line.
Consequently `B_tilde` is tangent to the HCS. Distributionally,

\[
 \nabla\cdot(\sigma_B\widetilde{\mathbf B})
 =\sigma_B\nabla\cdot\widetilde{\mathbf B}
  +\widetilde{\mathbf B}\cdot\nabla\sigma_B=0.
\]

The two one-sided limits satisfy

\[
 \mathbf B^+=-\mathbf B^-,\qquad
 |\mathbf B^+|=|\mathbf B^-|,\qquad
 \mathbf B^\pm\cdot\hat{\mathbf n}_{\rm HCS}=0.
\]

These statements apply to the pure finite-SCS region and its conservative
outer continuation. Inside the PFSS/SCS overlap, the inherited surface traced
from the SCS neutral line is denoted `S_tr`, the **transition sector surface**.
It is not assumed to be tangent to the final composite field. There the
antiparallel/tangency relations hold only if a separately constructed
common-flux-surface or finite-thickness branch passes its topology gate. The
SCS construction sector and the sector obtained from the unique non-neutral
photospheric footpoint of a final composite-field line are stored separately;
agreement is a diagnostic, not a relabeling operation.

There is no single vector value on an ideal zero-thickness sheet. A query on
the sheet must supply a side or receive a typed discontinuity result. Sector
connectivity and an explicit field-line-generated sheet surface are stored;
trilinear interpolation of the sign is forbidden. Signed, unsigned, and
sector-resolved flux must close on nested spheres. Zero signed flux is recovered
only when the PFSS boundary is flux balanced and the numerical sector map
preserves both polarity fluxes.

### 6.4 PFSS/SCS interface

At `R_i`, the mathematical boundary and polarity restoration guarantee

\[
 [B_r]=0,
\]

but generally

\[
 [\mathbf B_t]\ne0,
 \qquad
 \mathbf K_s=\frac{1}{\mu_0}\hat{\mathbf e}_r\times
   (\mathbf B_{\rm SCS}-\mathbf B_{\rm PFSS}).
\]

The discrete PFSS/SCS coupling uses one conservative normal-flux mortar trace;
it may not evaluate unrelated PFSS and SCS approximations and accept their
difference. Its `q_L-q_0` projection error must converge below the registered
interface tolerance. The tangential jump is a modeled surface current and
field-line kink, not only a derivative discontinuity. Evaluation, derivatives,
shock quadrature, and output are one-sided and record the jump and kink angle.

For production transport, replace the ideal jump by a declared overlap shell
`r_a=R_i <= r <= r_b=R_i+w_tr`, with
`r_b<=min(R_b,R_scs)`. The blend uses globally single-valued vector
potentials of the **signed, flux-balanced** PFSS and SCS fields. It never
constructs `A_SCS` from the unsigned `B_tilde`, whose required monopole flux
precludes a global single-valued vector potential on a spherical shell.

The reproducible baseline uses the Mie representation

\[
 \mathbf B=\nabla\times\nabla\times(P\mathbf r)
            +\nabla\times(T\mathbf r),
 \qquad
 \mathbf A=\nabla\times(P\mathbf r)+T\mathbf r.
\]

Only `ell>=1` signed modes are admitted; the surface means of `P` and `T` are
zero on every radius, and both fields use the same harmonic normalization,
frame, pole convention, and radial integration constants. `P` is recovered
from `B_r` and `T` from `(curl B)_r` through surface Poisson problems. The
signed SCS discontinuity is represented on a sheet-conforming sector grid or
an equivalent weak discretization, not spectrally smoothed across its sign
jump. A piecewise potential must satisfy `n_S cross [A]=0` to numerical
tolerance; otherwise its curl contains an unrequested sheet-delta field and
preparation fails. Signed flux closure, one-sided `curl(A)=B`, the tangential
traces, and the fixed-gauge constraints are verified before blending. A common
fixed gauge is mandatory because the interior blend depends on the potential
difference even though each endpoint field is gauge invariant.

After matching the conservative normal-flux trace and the verified potential
traces, use a quintic blend

\[
 \mathbf A_{\rm tr}=(1-\chi)\mathbf A_{\rm PFSS}
                    +\chi\mathbf A_{\rm SCS},
 \qquad
 \mathbf B_{\rm tr}=\nabla\times\mathbf A_{\rm tr},
\]

where, with `z=(r-r_a)/(r_b-r_a)`,

\[
 \chi(r)=6z^5-15z^4+10z^3.
\]

Thus `chi(r_a)=0`, `chi(r_b)=1`, and its first and second derivatives vanish
at both ends. Taking the curl after blending makes `div(B_tr)=0` analytically; blending
the magnetic components directly is forbidden. The transition's distributed
current, magnetic-energy change, maximum direction rotation per gyroradius/
mean-free-path, and residual force imbalance are reported. The wind is solved
on the completed transitioned field, not blended afterward. This construction
is a controlled divergence-free regularization of a model interface, not a
claim that its current or force balance was derived from coronal MHD.

Explicitly,

\[
 \mathbf B_{\rm tr}=(1-\chi)\mathbf B_{\rm PFSS}
 +\chi\mathbf B_{\rm SCS}
 +\nabla\chi\times(\mathbf A_{\rm SCS}-\mathbf A_{\rm PFSS}).
\]

Although `B_SCS` is tangent to the inherited transition sector surface
`S_tr`, in general

\[
 \mathbf B_{\rm tr}\cdot\hat{\mathbf n}_{S_{\rm tr}}=
 (1-\chi)\mathbf B_{\rm PFSS}\cdot\hat{\mathbf n}_{S_{\rm tr}}
 +[\nabla\chi\times(\mathbf A_{\rm SCS}-\mathbf A_{\rm PFSS})]
  \cdot\hat{\mathbf n}_{S_{\rm tr}}\ne0.
\]

Therefore `curl(A_tr)` proves solenoidality but not HCS tangency, sector
connectivity, or the absence of in-shell nulls and X-lines.

#### Transition-sheet topology and production policies

Schema 5 supports `exclude-clearance` as the first-release production policy.
The region `|d_tr|<=d_clear` around `S_tr`, together with every connected
component having ambiguous composite footpoints, is excluded from source
injection, particle transport, observer support, interpolation/derivative
stencils, and production field-line export. A particle segment reaching the
clearance is event-located and removed into a species-, weight-, momentum-,
and energy-resolved `TransitionSheetExcludedLoss` ledger; it is never passed
through by the ideal-HCS sign relabeling. Shock diagnostics remain available,
but source rate in the excluded region is exactly zero and separately
accounted. Clearance convergence is mandatory.

A future `common-flux-surface-qualified` policy must construct a sheet from
the final composite field and demonstrate both one-sided tangencies, a unique
two-sector partition, no self-intersection or terminating/null line, endpoint
matching to the exterior SCS sheet, and converged geometry/enclosed flux.
Footpoint labels alone cannot satisfy this test. A future
`finite-thickness-qualified` policy belongs to Stage 11A and requires its own
continuous-field crossing/drift operator. Neither future selector silently
falls back to an exclusion or to an ideal sign flip.

For the inherited zero-thickness surface, evaluate both one-sided transition
traces `B_tr^+` and `B_tr^-`; an unsided field value is undefined. The common
normal-flux mortar trace `B_n^*` exists only after

\[
 \epsilon_{[B_n]}=
 \frac{\int_{S_{\rm tr}}|\mathbf B_{\rm tr}^+\cdot\hat{\mathbf n}
 -\mathbf B_{\rm tr}^-\cdot\hat{\mathbf n}|\,dS}
 {\Phi_{\rm open}(r_b)}
\]

passes its registered normal-trace tolerance. The implementation then uses the
conservative mortar value for the grid-independent crossing measures

\[
 \Phi_{\rm cross}^{\rm abs}=\int_{S_{\rm tr}}
 |B_n^*|\,dS,
 \qquad
 \Phi_{\rm cross}^{\rm net}=\int_{S_{\rm tr}}
 B_n^*\,dS,
\]

\[
 \epsilon_{\rm cross}=\frac{\Phi_{\rm cross}^{\rm abs}}
 {\Phi_{\rm open}(r_b)}.
\]

Both raw one-sided fluxes are retained for audit. Directional diagnostics are
valid only when both `|B_tr^+|` and `|B_tr^-|` are at least the configured
positive `transition_minimum_field_tesla`; otherwise the point is typed
`weak-field-direction-inapplicable` and belongs to the exclusion mask. At a
valid point define the dimensionless jump support

\[
 a_{\rm jump}=\frac{|\mathbf B_{\rm tr}^+-\mathbf B_{\rm tr}^-|}
 {|\mathbf B_{\rm tr}^+|+|\mathbf B_{\rm tr}^-|}\in[0,1],
\]

the signed one-sided incidence cosines

\[
 \iota_\pm=\frac{\mathbf B_{\rm tr}^\pm\cdot\hat{\mathbf n}}
 {|\mathbf B_{\rm tr}^\pm|}\in[-1,1],
\]

and, only where
`a_jump>=transition_jump_support_fraction`, the antipodality defect

\[
 \delta_{\rm anti}=\arccos[\operatorname{clamp}
 (-\hat{\mathbf B}_{\rm tr}^+\cdot\hat{\mathbf B}_{\rm tr}^-,-1,1)]
 \in[0,\pi].
\]

Outside that support, `delta_anti` is typed `inapplicable-zero-jump`, because
coincident fields at a zero-amplitude endpoint would otherwise produce the
meaningless value `pi`; `iota_+/-` remain valid when their field thresholds
pass. These flux measures, incidence cosines, supported non-antipodality, null
census, sector mismatch, and immutable exclusion-mask geometry form the
pre-run **D9-background** product.

The checksummed clearance-convergence asset preregisters at least three
strictly ordered physical clearances together with the transition widths,
angular/radial grids, field-line tolerances, minimum-field threshold, request
set, and refinement order used at each member. It also declares absolute and
relative convergence bounds for every D9-background scalar/measure and for
each separate shock/source, observer/export, and runtime-loss consumer ledger.
The production clearance is an identified member and passes only when its
successive refined changes and all support-coverage fractions meet those
bounds; the asset cannot be generated or reweighted after inspecting SEP
results. Shock/source overlap, rejected observer/export requests, and runtime
particle losses are separate generation-indexed consumer ledgers and are never
hashed into the immutable background identity. A traced-line count may
accompany those ledgers but cannot be a D9 acceptance metric because it depends
on seed density.

The clearance is also subject to preregistered **consumer-impact budgets**;
convergence of the background mask is necessary but is not an acceptance cap.
For surface generation `g`, registered interval `I`, and every `t in I`, let
`P_cf(s,g,t)` be the physical shock patches for species `s` that satisfy every
source gate (geometry, upstream authority, fast/jump/criticality, topology,
species, abundance, spectrum, and radial envelope) **except** the
`S_tr`-clearance gate. After exact patch splitting at the mask boundary, define

\[
 A_{\rm cf}(t)=\sum_{j\in P_{\rm cf}(t)}A_j(t),
 \qquad
 A_{\rm tr}(t)=\sum_{j\in P_{\rm cf}(t)}A_j(t) I_{{\rm tr},j}(t),
 \qquad f_A^{\rm tr}(t)=A_{\rm tr}(t)/A_{\rm cf}(t),
\]

where `I_tr,j` is one only for the clearance-excluded child. The same
counterfactual set supplies the number- and energy-rate denominators

\[
 \dot N_{\rm cf,s}(t)=\sum_j\int \dot q_{{\rm cf},s,j}(p,t)\,dp,
 \qquad
 \dot E_{\rm cf,s}(t)=\sum_j\int K_s(p)\dot q_{{\rm cf},s,j}(p,t)\,dp,
\]

and the excluded numerators insert `I_tr,j` under the sums. Thus
`f_N_tr=dot(N_tr,s)/dot(N_cf,s)` and
`f_E_tr=dot(E_tr,s)/dot(E_cf,s)` measure what the clearance alone removes;
their denominators are not the whole candidate dome and do not omit another
failed source gate merely to make the fractions small. For an interval, the
normative measures are ratios of time integrals,

\[
 F_A^{\rm tr}(I)=\frac{\int_I A_{\rm tr}(t)\,dt}
                           {\int_I A_{\rm cf}(t)\,dt},\quad
 F_N^{\rm tr}(I)=\frac{\int_I \dot N_{\rm tr,s}(t)\,dt}
                           {\int_I \dot N_{\rm cf,s}(t)\,dt},\quad
 F_E^{\rm tr}(I)=\frac{\int_I \dot E_{\rm tr,s}(t)\,dt}
                           {\int_I \dot E_{\rm cf,s}(t)\,dt}.
\]

The code never substitutes an integral or arithmetic mean of instantaneous
fractions. It stores both sides of every ratio, the instantaneous maxima, and
the ratio-of-integrals for every registered interval. `A_cf(t)=0` and
`integral_I A_cf dt=0` are typed `inapplicable-no-counterfactual-area`;
zero number- or energy-rate denominators are typed
`inapplicable-no-counterfactual-source`. None is converted to numerical zero
or counted as a passing ratio.

Every finite observer or exported-line footprint `F` is budgeted by its
unsigned represented magnetic flux,

\[
 \Phi_F=\int_F |\mathbf B\!\cdot\!\hat{\mathbf n}_F|\,dA,
 \qquad
 f_{\Phi,F}^{\rm tr}=
 \frac{\int_{F\cap C_{\rm tr}}
 |\mathbf B\!\cdot\!\hat{\mathbf n}_F|\,dA}{\Phi_F},
\]

using the same traced footprint and quadrature that define its physical
measure. A finite footprint with `Phi_F=0` is typed
`inapplicable-zero-footprint-flux`; its ratio is not evaluated and a required
event-grade consumer follows the preregistered missing-measure action rather
than receiving a fabricated zero. For a moving footprint the interval measure
is likewise the ratio of the time-integrated excluded unsigned flux to the
time-integrated total unsigned flux. A characteristic line or mathematical point observer has no finite
flux measure: it receives a typed `valid` or `rejected-clearance` result and
its flux fraction is `inapplicable-point-measure`, never a fictitious zero or
one. The checksummed consumer-budget asset supplies aggregate and per-species,
birth-energy-bin, stable patch-lineage, and birth-time-bin maxima. Exceeding
any area, counterfactual-rate, footprint-flux, or runtime-loss maximum marks an
event run `not-event-grade` with the exact numerator, denominator, stratum,
generation, interval, and violated bound in the manifest. Thresholds are
registered before SEP output is inspected and cannot be reweighted afterward.

For isotropic Parker transport, merely labelling the direction jump is not a
conservative crossing rule: one-sided density, advection, and diffusion tensors
require continuity of normal particle flux. The production baseline therefore
uses a finite, divergence-preserving transition resolved by the field, wind,
mesh, and particle step. A zero-thickness Parker interface operator is legal
only after a finite-volume/Fokker--Planck and stochastic-interface test proves
the same conservative transmission law.

A gyrotropic focused particle cannot in general be mapped through an arbitrary
instantaneous field rotation from `(p,mu)` alone because gyrophase has been
averaged out. It likewise requires the resolved transition or a separately
verified operator that samples the unresolved gyrophase while preserving
physical velocity and energy. Until one is qualified, either mover is rejected
when the interface treatment is unresolved or the kink exceeds its verified
tolerance. Overlap-minimized coupling can reduce the kink but does not prove
vector continuity.

### 6.5 Conservative rotating-footpoint mapping

The finite-shell field is radial at `R_scs`, so the Parker extension begins
continuously with `B_theta=B_phi=0`. Let `alpha` label longitude at `R_scs` and
let `Phi(r,theta,alpha)` be the rotating-frame longitude of that labeled tube:

\[
 \Phi(R_{\rm scs},\theta,\alpha)=\alpha,
 \qquad
 \frac{d\Phi}{dr}=-K,
\]

\[
 K=\frac{\Omega_F-u_\phi/(r\sin\theta)}{u_r}.
\]

For the baseline angular-momentum-conserving azimuthal velocity,

\[
 u_\phi=\Omega_F\frac{R_{\rm scs}^2}{r}\sin\theta,
 \qquad
 K=\frac{\Omega_F(1-R_{\rm scs}^2/r^2)}{u_r}.
\]

Define the forward and inverse longitude-map Jacobians

\[
 A_\phi=\frac{\partial\Phi}{\partial\alpha},
 \qquad
 J_\phi=\frac{\partial\alpha}{\partial\phi}=A_\phi^{-1}.
\]

For an Eulerian `K(r,theta,phi)`,

\[
 \frac{dA_\phi}{dr}=-(\partial_\phi K)A_\phi,
 \qquad
 \frac{dJ_\phi}{dr}=(\partial_\phi K)J_\phi,
 \qquad A_\phi(R_{\rm scs})=1.
\]

At an evaluation point, the periodic monotone equation
`Phi(r,theta,alpha)=phi` is inverted for `alpha`. The field is

\[
 B_r=\left(\frac{R_{\rm scs}}{r}\right)^2
 B_{r,\rm scs}(\theta,\alpha)J_\phi,
 \qquad B_\theta=0,
\]

\[
 B_\phi=-r\sin\theta\,K B_r.
\]

With the same `u_phi` and `K`, the one-sided rotating-frame velocity obeys

\[
 \mathbf u'=\mathbf u-\boldsymbol\Omega_F\times\mathbf x
            =\frac{u_r}{B_r}\mathbf B,
 \qquad
 u_s=\frac{u_r|\mathbf B|}{|B_r|},
\]

away from the ideal HCS. Thus the outer magnetic mapping and the Stage 2
field-aligned wind use the same velocity definition; they are not two
independent Parker prescriptions.

This follows from conservation of flux,

\[
 r^2B_r\,d\phi=R_{\rm scs}^2B_{r,\rm scs}\,d\alpha.
\]

The old scalar mapping is recovered only when the Eulerian
`partial_phi K=0` along every characteristic (equivalently, no longitudinal
dependence in the source-labeled solution), for which
`A_phi=J_phi=1`. Longitude mapping, sector/HCS advection, `B_phi`, and the
Jacobian are one inseparable operation. A latitude-only rotation rate may be
used in an inertial or fixed-rate frame; a longitude-dependent rotation law
requires the same general Jacobian treatment.

`A_phi=0` is a characteristic/stream crossing and `A_phi<0` is a folded,
multivalued ballistic map. Production requires `A_phi>A_phi_min` throughout
the physical domain. The reported dimensionless fold margin is defined
unambiguously as `m_fold=A_phi-A_phi_min`; production requires `m_fold>0` and
fails transactionally before a fold. A diagnostic run
may mark downstream values invalid only if every consumer rejects them. The
Jacobian makes the pre-fold magnetic map solenoidal; it does not model stream
interaction or CIR formation after a fold.

The often-requested case `R_w<R_scs` is deliberately **not** obtained by
starting the radial formulas above at a nonradial surface. The defined future
baseline is a **steady corotating spatial pushforward**, not an arbitrary
time-dependent twist. Use reference tube coordinates
`X=(varrho,theta_0,alpha)` with `varrho>=R_w` and define the
orientation-preserving map `x=F(X)` by

\[
 \mathbf F(R_w,\theta_0,\alpha)
   =R_w\hat{\mathbf e}_r(\theta_0,\alpha),
 \qquad
 \frac{\partial\mathbf F}{\partial\varrho}
   =\frac{\mathbf u'(\mathbf F)}
   {\mathbf u'(\mathbf F)\cdot\hat{\mathbf e}_r(\mathbf F)},
 \qquad \mathbf u'=\mathbf u-\boldsymbol\Omega_F\times\mathbf x.
\]

The denominator must remain strictly positive, `|F|=varrho`, and the two
tangential boundary derivatives are those of the identity sphere. The full
corotating velocity, frame, epoch, outer-boundary condition, and any smooth
start/ramp are immutable inputs. The steady ideal condition
`curl(u' cross B)=0` is enforced; the baseline chooses the field-aligned branch
`u' cross B=0` on every open characteristic. These equations and the boundary
trace determine `F`; a different diffeomorphism is a different future model.

With deformation gradient `D F`, determinant `J=det(D F)>0`, and reference
field `B_0(X)`, the contravariant Piola transform is

\[
 \mathbf B(\mathbf x)=\frac{D\mathbf F}{J}\,
 \mathbf B_0(\mathbf X),\qquad \mathbf X=\mathbf F^{-1}(\mathbf x).
\]

At `R_w`, the Piola field must reproduce the complete one-sided unwound vector
field and its conservative normal-flux mortar trace. The smooth-start profile
must also match the declared derivative order; a tangential or normal mismatch
would introduce an unrequested surface current or magnetic charge and is
fatal. A later genuinely time-dependent inflow formulation would instead need
surface labels plus a parcel launch time `tau` and impose identity only at each
parcel's crossing, `F_{tau,tau}(X_tau)=X_tau`; it cannot reuse this steady
selector or impose a fixed-boundary identity on a Lagrangian flow with nonzero
outflow.

The Piola transform preserves magnetic flux and `div(B)=0`; it is not a
complete model by itself. The provider must also push forward the HCS,
separatrix, PFSS/SCS interfaces and their one-sided normals, transform all
gradients and metric factors, and re-solve the wind/mass-per-flux problem on
the deformed nonradial tubes. `J<=0`, nonpositive radial projection, loss of
inverse uniqueness, intersection of distinct tubes, or a failed topology/event
map is a transactional fold failure. The future `CPL3D10` contract therefore
verifies the steady-spatial equations and boundary trace above and requires
manufactured Piola flux/divergence tests, interface and sector pushforward,
mass-per-flux closure, gradient and focusing-length convergence, 3-D/1-D
parity, fold/turning rejection, and reduction to the radial schema-5 formulas
when `R_w=R_scs`. A time-dependent launch-time provider requires a different
test ID. Until that provider and test record exist, `R_w` is derived as
`R_scs` and no input selector can request the generalized branch.

### 6.6 Plasma continuation and coupled structured-wind solve

Plasma is carried by the same source label `alpha` used by the magnetic field.
For a steady corotating flow,

\[
 r^2\rho u_r=
 R_{\rm scs}^2\rho_{\rm scs}u_{r,\rm scs}J_\phi,
\]

or

\[
 \rho(r,\theta,\phi)=\rho_{\rm scs}(\theta,\alpha)
 \frac{u_{r,\rm scs}}{u_r}
 \left(\frac{R_{\rm scs}}{r}\right)^2J_\phi.
\]

Thus the angular stream-tube area scales as `r^2 A_phi`, not simply `r^2`,
when the wind is longitude structured. The wind entropy/polytropic constant is
advected with the same label and pressure is reconstructed from the selected
closure. Magnetic and plasma mappings may not use different inversions.

Because `K` depends on `u_r`, every structured branch uses the coupled mapping
iteration specified in Section 5.5. The polytropic branch obtains `u_r` from
its transonic tube solution and must satisfy the solved field-aligned momentum
equation. The empirical-kinematic branch uses its selected channel's declared
velocity component/reference frame, converts to `u_r` only through the
qualified tube geometry, joins the inner density and outer implied density,
and derives the final `u_s` from its single `eta_m`. It reports/gates rather
than solves its momentum residual. Both branches must converge `div(B)=0` and
`div(rho u')=0` and conserve mass and magnetic flux per source-label interval.
No independent one-AU density law may overwrite either solution.

### 6.7 Empirical heliospheric plasma sheet

An optional tube-labeled plasma-sheet closure uses angular distance `d_N`
between the tube footpoint on `R_i` and the PFSS neutral line:

\[
 C_{\rm ps}(d_N)=1+(C_{\rm ps,0}-1)
 \exp[-(d_N/w_{\rm ps})^2].
\]

Require `C_ps,0>=1`, `w_ps>0`, and the first-release
`fixed-temperature` thermodynamic rule. With base-density normalization, apply
the factor to the authoritative tube base state before solving the wind:

\[
 \rho_0\rightarrow C_{\rm ps}\rho_0,
 \qquad p_0\rightarrow C_{\rm ps}p_0,
 \qquad T_0\ \text{unchanged}.
\]

This preserves `p/rho` and the EOS and changes the tube mass-flux
normalization. With `mass-per-magnetic-flux`,
`radial-flux-density-with-mapped-field`, or
`colocated-density-velocity-field` normalization, the same factor instead
modifies that target normalization and
the base density remains an eigen-solution output; applying both changes would
overdetermine the tube and is rejected. It is never a density-only cell
overwrite. Wind, characteristic speeds, turbulence, and shock quantities are
derived only after the modified tube solution is complete. The closure is
empirical: it does not enforce transverse pressure balance or predict radial
sheet-thickness evolution [37, 38]. Consequently it is diagnosed as a smooth
volume transition, never as a sharp traction jump. The provider evaluates the
Section 4.6 `R_mom` vector through the sheet using the complete declared force
inventory and reports pointwise values, physical-volume quantiles, and the
maximum. Analytic use is `diagnostic-kinematic`. Event-grade use of a nonzero
sheet selects `bounded-approximation` and preregisters density/width
uncertainty, angular/radial resolution and width ensembles, and absolute and
relative force-density bounds; passing those bounds does not make the
empirical sheet an MHD equilibrium.

### 6.8 One-sided numerical representation and particle interfaces

AMPS center/vertex storage carries magnetic region, both PFSS and composite topology, sector,
interface side, discontinuity identity, signed HCS distance, and mapping
validity. Evaluation determines region and side first, evaluates the
continuous one-sided field, and only then applies the discrete sector sign.
Halo exchange includes these categorical fields. No derivative is formed by
differencing through `R_i`, the HCS, or the open/closed separatrix.

Parallel field lines are tangent to the ideal HCS only in the pure SCS/Parker
regions and in a separately common-flux-surface-qualified transition; they are
also tangent to the open/closed separatrix. The inherited `S_tr` inside the
overlap is not assumed to be a flux surface. Under `exclude-clearance`, no
particle step may enter its excluded neighborhood: the segment is clipped at
the clearance boundary and recorded as `TransitionSheetExcludedLoss` rather
than passed through or relabeled. A cross-field/full-velocity step that crosses
an enabled ideal sector surface **outside** the transition is split exactly. Position,
physical velocity, momentum magnitude, `w_out`, and `w_in` remain continuous;
for a pure sign reversal,

\[
 \mu_2=-\mu_1,
 \qquad \mu_2\hat{\mathbf b}_2=\mu_1\hat{\mathbf b}_1.
\]

This is a coordinate relabeling for the pure-SCS ideal HCS, not for `S_tr` and not a finite-thickness current-sheet or
current-sheet-drift model [47, 48]. A run requesting antisymmetric drift or HCS
drift must select a separately verified finite-thickness/sheet-drift provider
or fail. A general PFSS/SCS kink uses the interface policy in Section 6.4, not
the HCS sign-flip rule. Shock patches cut by the HCS or separatrix are split;
opposite vectors are never averaged into an artificial weak field.

---

## 7. Turbulence and particle-scattering closure

### 7.1 Directional wave energy

Let `w_+` and `w_-` be energy densities of Alfvénic fluctuations propagating
parallel and antiparallel to the local signed `B`. Those labels are not the
same as physical outward/inward propagation in both magnetic sectors. For an
open tube with outward tangent `t_hat`, define

\[
 \sigma_B=\operatorname{sgn}(\hat{\mathbf b}\cdot\hat{\mathbf t}),
\]

\[
 w_{\rm out}=\frac{1+\sigma_B}{2}w_+
              +\frac{1-\sigma_B}{2}w_-,
 \qquad
 w_{\rm in}=\frac{1-\sigma_B}{2}w_+
             +\frac{1+\sigma_B}{2}w_-.
\]

The physical and field-relative cross helicities are

\[
 \sigma_c^{\rm out}=\frac{w_{\rm out}-w_{\rm in}}
                           {w_{\rm out}+w_{\rm in}},
 \qquad
 \sigma_c^{B}=\sigma_B\sigma_c^{\rm out}.
\]

Thus an outward wave is `w_+` in a positive-polarity sector and `w_-` in a
negative-polarity sector. At an ideal sign reversal, `w_out` and `w_in` remain
the physical labels while `w_+` and `w_-` exchange. The authoritative stored
quantities are `w_out`, `w_in`, and `sigma_B`; field-relative values are
derived.

For total wave energy `w=w_out+w_in`, ideal Alfvén-wave equipartition gives

\[
  \delta B^2=\mu_0 w,
  \qquad
  \delta v^2=\frac{w}{\rho}.
\]

Every output must state whether `delta B^2` denotes the total variance or one
direction/polarization.

### 7.2 Minimal prescribed model

The existing prescribed model can define

\[
 w(r)=w_{\rm ref}\left(\frac{r}{r_w}\right)^{-\alpha_w},
 \qquad
 w_{\rm out,in}=\frac{1\pm\sigma_c^{\rm out}}{2}w,
\]

and

\[
 L_c(r)=L_{c,\rm ref}\left(\frac{r}{r_w}\right)^{\alpha_L}.
\]

This is a controlled parameterization, not a self-consistent turbulence
solution. Here `r_w` is the turbulence reference radius, distinct from the
source-termination radius used in Section 10.6. Its amplitudes and exponents
must be input values with provenance, and every prescribed or WKB realization
must pass the single perturbative-validity contract in Section 7.5.

### 7.3 WKB outward-wave option

For a stationary, non-dissipative outward wave on a flux tube, wave-action
conservation gives

\[
 \mathcal S_{\rm out}=
 A\,w_{\rm out}\frac{(u_s+v_A)^2}{v_A}=\text{constant}.
\]

Thus

\[
 w_{\rm out}(s)=w_{{\rm out},0}
 \frac{A_0}{A(s)}
 \frac{v_A(s)}{v_{A,0}}
 \left[\frac{u_{s,0}+v_{A,0}}{u_s(s)+v_A(s)}\right]^2.
\]

Subscript `0` denotes the unique outward-tube reference point `s_0` at which
`r(s_0)=r_w`, using the same `turbulence.reference_radius_m` as the prescribed
branch; `A_0`, `u_s,0`, `v_A,0`, and `w_out,0` are all evaluated there from one
prepared snapshot. A tube with no unique valid crossing of that reference
sphere is rejected rather than assigned a nearby node.

An inward component cannot be invented by applying the same expression through
the Alfvén critical point. It requires an outer boundary, a measured cross
helicity, or an explicit reflection/cascade equation. Therefore the first WKB
implementation uses `w_in=0`, equivalently `sigma_c_out=1`. A nonzero inward
component requires a separately identified, versioned reflection/boundary
closure and is not selectable in the first production schema. A global
`w_-=0` rule is forbidden because it would describe inward physical propagation
in a negative-polarity sector. Reflection-driven turbulence is a later roadmap
stage [6, 7]. On closed loops there is no unique outward direction; the
explicit closed-field policy from Section 4.6 applies.

### 7.4 Mean-free-path closure

For a controlled transport baseline, the parallel mean free path may be
specified as

\[
 \lambda_\parallel(r,\mathcal R)=
 \lambda_0
 \left(\frac{r}{r_\lambda}\right)^{\alpha_r}
 \left(\frac{\mathcal R}{\mathcal R_0}\right)^{\alpha_R},
\]

where rigidity is

\[
 \mathcal R=\frac{pc}{|q_s|}.
\]

The corresponding isotropic parallel diffusion coefficient is

\[
 \kappa_\parallel=\frac{v\lambda_\parallel}{3}.
\]

`r_lambda` is the mean-free-path reference radius and is unrelated to the
zero-source radius in Section 10.6.

This closure must remain available even after wave-energy initialization is
implemented, because a wave amplitude alone does not determine particle
scattering without specifying the spectrum, resonance range, polarization, and
nonlinear broadening [8].

The optional `smooth-broken-power-law` closure can represent a change in
radial regime without a derivative discontinuity. With `x=r/r_b` and
`nu>0`,

\[
 \lambda_\parallel=\lambda_b x^{\alpha_{\rm in}}
 \left(\frac{1+x^\nu}{2}\right)^{
   (\alpha_{\rm out}-\alpha_{\rm in})/\nu}
 \left(\frac{\mathcal R}{\mathcal R_0}\right)^{\alpha_R}.
\]

It is normalized at `r_b` and approaches the requested inner and outer slopes.
It remains an empirical closure, not an intrinsically more physical law. A
single power law remains valid for controlled comparisons but cannot represent
a radial change of regime and may be inadequate for joint near-Sun/1-AU fits.
Mean-free-path inference is event-, energy-, turbulence-, and radius-dependent
[45, 46]. A turbulence-derived QLT/SOQLT/nonlinear option is a separate model
that must define spectrum geometry, dissipation range, resonance broadening,
and the treatment of the 90-degree resonance problem.

### 7.5 Validity of an uncoupled wave field

The baseline wind omits wave pressure, wave heating, and back reaction. It is
therefore valid only while the omitted wave force is demonstrably small for
the selected background and while any small-amplitude assumptions made by the
selected turbulence/scattering closure remain valid. With the total-wave
convention in Section 7.1,

\[
 \epsilon_{\delta B}=\frac{\sqrt{\mu_0(w_{\rm out}+w_{\rm in})}}{|\mathbf B|},
 \qquad
 p_w=\frac{w_{\rm out}+w_{\rm in}}{2},
\]

\[
 \epsilon_{w,p}=\frac{p_w}{p},
 \qquad
 \epsilon_{w,\Pi}=\frac{p_w}{p+\rho u_s^2},
 \qquad
 \epsilon_{w,B}=\frac{p_w}{|\mathbf B|^2/(2\mu_0)}
                =\epsilon_{\delta B}^2.
\]

`epsilon_w,p` and `epsilon_w,Pi` are mandatory diagnostics but neither is the
normative back-reaction error: a pressure magnitude does not determine its
gradient. The omitted field-aligned wave acceleration is

\[
 a_{w,s}=-\frac{1}{\rho}\frac{dp_w}{ds}.
\]

Define the retained-acceleration scale without cancellation,

\[
 a_{\rm retained}=
 \left|u_s\frac{du_s}{ds}\right|
 +\left|\frac{1}{\rho}\frac{dp}{ds}\right|
 +\left|\frac{d\Phi_{\rm eff}}{ds}\right|.
\]

Let `tau_w,abs` be the maximum absolute wave-acceleration uncertainty obtained
from a checksummed derivative/refinement qualification asset; it is not a
free denominator floor. The asset-derived value must not exceed the configured
`maximum_absolute_wave_acceleration_m_per_s2` cap. The pointwise production
condition is

\[
 |a_{w,s}|\le \tau_{w,\rm abs}
 +\epsilon_{w,F}^{\max}a_{\rm retained}.
\]

For reporting and the distribution gate define

\[
 \epsilon_{w,F}^{\rm excess}=
 \frac{\max(0,|a_{w,s}|-\tau_{w,\rm abs})}{a_{\rm retained}}
\]

when `a_retained>0`. It is zero when both the numerator and retained scale are
zero, and typed `+infinity/fail` when the numerator is positive but the scale
is zero. Preparation evaluates signed `a_w,s`, the absolute allowance,
`epsilon_w,F^excess`, all four amplitude ratios, and their distributions at
every point where wave energy is applicable, including source patches,
observer support, active mover cells, and exported-line nodes. Production of
an uncoupled background requires both the pointwise maximum condition and the
separately configured area/volume-weighted quantile bound to pass. The
qualification asset records the refinement sequence and estimator; changing
it changes the fingerprint, and a coarser derivative may not claim a smaller
uncertainty than its verified convergence envelope.

The action associated with `epsilon_deltaB` is closure dependent. It remains
a fatal validity gate for linear WKB and every QLT/small-amplitude coefficient
that assumes `delta B/B << 1`. For an empirical prescribed-amplitude branch
whose independently prescribed coefficient makes no small-amplitude claim, it
may be diagnostic-only, but it is still emitted. `epsilon_w,p` is likewise
never used as a universal fatal threshold in a cold supersonic wind; a large
value can coexist with a small force gradient, while a modest value with a
sharp gradient can fail `epsilon_w,F`. The identity
`epsilon_w,B=epsilon_deltaB^2` remains a required consistency check.

A closed authority selected as `direct-mean-free-path` carries an explicitly
inapplicable wave state and is not assigned artificial zero wave energy merely
to pass these ratios. A closed authority selected as `invalid` is legal only
when topology/reachability analysis proves that no source, active particle
cell, observer support, or exported line can enter it; otherwise preparation
fails. Exceeding an applicable closure-specific amplitude gate or the
force-residual gate is fatal for the uncoupled model. These gates do not prove
that wave dynamics are negligible in all omitted equations, but they prevent a
nominally test-particle turbulence field from exerting an unmeasured force
while the wind equations ignore it. A wave-driven wind must instead solve the
coupled momentum and energy equations under a different authority.

### 7.6 Post-reference self-generated waves and future closure tiers

The schema-5 release does **not** evolve waves generated by the released SEP
population. Its prescribed and WKB wave energies are prepared independently
of the particles, and its direct mean-free-path closures remain immutable
background coefficients. Consequently, the present model cannot capture the
streaming-instability feedback in which an outwardly anisotropic SEP
population amplifies resonant Alfvénic fluctuations, reduces its own parallel
mean free path, and changes subsequent escape and transport [24, 25, 59, 60]. This
omission applies after particles cross the finite reference surface as well as
inside the unresolved shock-to-reference region; the latter contribution is
already part of the calibrated first-passage source and must not be counted a
second time by a transport-side correction.

Development should proceed through two explicitly different future tiers:

1. The future **`foreshock-distance-proxy` empirical post-reference scattering
   sensitivity** prescribes
   `lambda_parallel=F_foreshock*lambda_parallel,ambient` with
   `0<F_foreshock<=1`. The factor is a positive,
   smooth function of
   declared upstream distance, event time, species/rigidity, source patch,
   local shock obliquity/Mach state, and background generation. Its support is
   entirely on the modeled post-reference upstream side; it never modifies the
   unresolved front-to-reference accelerator already represented by
   `eta*g(p)`. The proxy must
   be a versioned calibration asset with finite support, uncertainty, bounded
   coefficients, recovery of the registered ambient mean free path at its
   outer edge, and convergence tests at every support edge. It is a sensitivity
   model only: it does not evolve wave
   energy, cannot be called self-generated turbulence, and cannot use a wave
   amplitude to infer scattering without an explicit spectrum and resonance
   closure.
2. A **coupled particle--wave model** must evolve a signed-direction,
   wavenumber-resolved wave spectrum with declared advection/refraction,
   SEP-driven growth, damping, reflection, and cascade operators. A schematic
   balance such as

   \[
     \frac{\partial W^\pm(k)}{\partial t}
     +{\cal T}^\pm[W^\pm]
     =2\,[\Gamma_{\rm SEP}^\pm-\Gamma_{\rm damp}^\pm]W^\pm
      +{\cal C}^\pm[W^+,W^-]
   \]

   is not an executable closure until every operator, frame, sign convention,
   resonance condition, spectral boundary, and normalization is specified.
   The evolved spectrum must determine `D_mumu` or `lambda_parallel` through a
   validated scattering theory, consume the same represented SEP distribution
   used by the mover, and close particle/wave energy and momentum exchange.
   If its wave force or heating violates Section 7.5, the solar-wind momentum
   and energy equations must be coupled as well rather than retaining the
   polytropic background unchanged.

   A separately qualified coupled solver may publish an immutable offline
   `W^+`/`W^-` spectral asset for later consumption. That asset must bind units,
   frame, wavenumber convention, resonance mapping, interpolation/coverage,
   solver version, tolerances, conservation residuals, and the exact
   background/shock, source spectrum and normalization, species/weight,
   transport/return-policy, particle-distribution, coupling-iteration, and
   coupling-cadence fingerprints. If those particle/source bindings are
   absent, the product is only a prescribed external wave field and cannot be
   described as a replay of a self-consistent particle--wave solution.
   Consuming a fully bound asset is an offline coupled-physics
   approximation, not runtime particle--wave feedback.

Neither tier is currently implemented. In particular, no configuration may
advertise a `srcSEP`-calibrated self-excited-wave option unless an actual
solver or immutable calibration asset with the contracts above exists. The
empirical proxy and coupled solver require separate capability identifiers,
provenance, restart state, diagnostics, and validation campaigns; one may not
silently fall back to the other.

---

## 8. CME piston and candidate-front geometry

### 8.1 Why an ellipsoid is used

A global spherical shell would place a shock around the entire Sun and would
give identical nose and flank curvature. A solar-anchored ellipsoidal dome can
represent unequal radial and lateral expansion, a finite angular extent, and
different shock normals and obliquities over the surface. Ellipsoidal shock
reconstruction has been used to distinguish a bubble-like shock from the
underlying CME ejecta [9, 10].

### 8.2 Implicit surface

The baseline uses a fully specified radial-principal-axis convention. For
model-frame latitude `lambda`, longitude `phi`, and lateral-axis tilt `gamma`,
define the right-handed basis

\[
 \hat{\mathbf e}_r=(\cos\lambda\cos\phi,
                     \cos\lambda\sin\phi,\sin\lambda),
\]

\[
 \hat{\mathbf e}_E=(-\sin\phi,\cos\phi,0),
 \qquad
 \hat{\mathbf e}_N=(-\sin\lambda\cos\phi,
                    -\sin\lambda\sin\phi,\cos\lambda),
\]

so that `e_r cross e_E = e_N`. A positive tilt is a right-handed rotation
about `e_r`:

\[
 \hat{\mathbf e}_1=\cos\gamma\,\hat{\mathbf e}_E
                  +\sin\gamma\,\hat{\mathbf e}_N,
 \qquad
 \hat{\mathbf e}_2=-\sin\gamma\,\hat{\mathbf e}_E
                  +\cos\gamma\,\hat{\mathbf e}_N.
\]

The center, orientation, and quadratic form are

\[
 \mathbf c(t)=d_c(t)\hat{\mathbf e}_r,
 \qquad
 R=[\hat{\mathbf e}_r,\hat{\mathbf e}_1,\hat{\mathbf e}_2],
\]

\[
 D(t)=\operatorname{diag}(a_r^{-2},a_1^{-2},a_2^{-2}),
 \qquad Q(t)=R D(t)R^T,
 \qquad \dot R=0.
\]

This convention is compatible with radial-axis ellipsoid reconstructions such
as the Kwon/PyThea parameterization after an explicit adapter conversion
[9, 50]. The
baseline orientation is fixed. A future arbitrary-orientation history must use
a normalized quaternion or rotation matrix and angular velocity; longitude,
latitude, and a scalar tilt are not sufficient to describe an arbitrary
time-dependent orientation.

Two physically distinct extensions must remain separate. **Center-path
deflection** changes the direction of the ellipsoid center, so that

\[
 \mathbf c(t)=d_c(t)\hat{\mathbf e}_r(t),\qquad
 \dot{\mathbf c}=\dot d_c\hat{\mathbf e}_r
                 +d_c\dot{\hat{\mathbf e}}_r .
\]

The second term is a nonradial translational velocity and therefore changes
the surface normal speed even if the body does not rotate. **Body rotation**
changes the principal-axis attitude `R(t)` and hence `Q(t)` about the moving
center. Changing `Q` cannot stand in for a deflected center trajectory, and
changing the center direction cannot stand in for body rotation. Schema 5
implements neither time-dependent operation: it keeps `e_r`, `R`, longitude,
latitude, and tilt fixed while evolving only center distance and semiaxes.
Observed and modeled low-coronal CME deflection/rotation demonstrate why these
future histories must be fit independently rather than replaced by a generic
angle perturbation [64, 65].

A future history must also declare how attitude is transported along a
deflected path. In a radial-axis-locked branch, the first principal axis
remains `e_r(t)`; the two lateral axes require a nonsingular transported basis
plus a separately fitted spin/tilt about `e_r`, and `Q_dot` includes that
geometric reorientation. In a fully arbitrary-attitude branch, `R(t)` is an
independent `SO(3)` history and radial alignment is not assumed. Mixing these
conventions or recomputing lateral axes independently at each snapshot would
create artificial angular velocity and shock-normal speed.

The surface is

\[
 F(\mathbf x,t)=
 [\mathbf x-\mathbf c(t)]^TQ(t)
 [\mathbf x-\mathbf c(t)]-1=0.
\]

The unnormalized surface normal is

\[
 \nabla F=2Q(\mathbf x-\mathbf c),
\]

and

\[
 \hat{\mathbf n}=\frac{\nabla F}{|\nabla F|}.
\]

For a moving implicit surface, the exact normal speed is

\[
 V_{\rm sh,n}=-\frac{\partial F/\partial t}{|\nabla F|},
\]

where

\[
 \frac{\partial F}{\partial t}=
 -2\dot{\mathbf c}\cdot Q(\mathbf x-\mathbf c)
 +(\mathbf x-\mathbf c)^T\dot Q(\mathbf x-\mathbf c).
\]

This equation includes center translation, independent axis expansion, and,
when explicitly implemented, orientation change. The code evaluates `Q_dot`
analytically. Finite differencing tessellated surfaces is not the authoritative
normal-speed calculation. Translation, expansion, and rotation contributions
may be reported separately, but their sum from the implicit-surface equation
is authoritative.

### 8.3 Solar-anchored dome

The full mathematical ellipsoid may extend below the photosphere. The physical
candidate front used by the application is

\[
  \mathcal S_{\rm front}(t)=
  \{\mathbf x:F(\mathbf x,t)=0,\ r\ge R_\odot\}.
\]

The below-surface portion is clipped, producing a dome rooted near the eruption
site. The intersection with `r=R_sun` must be computed topologically rather
than by sampling a few points. No patch may be generated inside the registered
AMPS solar boundary.

Exactly one initial radial parameterization is accepted:

1. `center-and-principal-axes`: input `d_c`, `a_r`, `a_1`, and `a_2`; or
2. `apex-and-principal-axes`: input `r_apex`, `a_r`, `a_1`, and `a_2`, then
   derive `d_c=r_apex-a_r`.

Under the radial-axis convention,

\[
 r_{\rm apex}=d_c+a_r,
 \qquad v_{\rm apex}=\dot d_c+\dot a_r.
\]

These sums are not valid for an arbitrary orientation. Validation requires
`d_c>0`, positive semiaxes throughout the declared history, the derived apex
inside the configured initial range, and a topologically valid solar-sphere
intersection. Within one geometry record, an apex is never supplied
independently when center distance is the selected radial parameterization.
The same exclusivity applies separately inside an enabled piston record; it
does not couple or forbid independently specified front and piston surfaces.

### 8.4 Smooth prescribed kinematics

Near-Sun CME motion uses an acceleration-phase model rather than a
heliospheric drag law. Apply the following independent law to every
`q in {d_c,a_r,a_1,a_2}`. Let `q_ref` be the value at `t_ref`; let the rate
transition begin at `t_s,q`, last `T_q>0`, and change from `qdot_0` to
`qdot_1`. Define

\[
 S(\tau)=10\tau^3-15\tau^4+6\tau^5,
 \qquad
 I(\tau)=\frac52\tau^4-3\tau^5+\tau^6,
 \qquad \tau=\frac{t-t_{s,q}}{T_q}.
\]

Before the transition, `q=q_ref+qdot_0(t-t_ref)`. During
`0 <= tau <= 1`,

\[
 \dot q(t)=\dot q_0+(\dot q_1-\dot q_0)S(\tau),
\]

\[
 q(t)=q_{\rm start}+\dot q_0(t-t_{s,q})
      +(\dot q_1-\dot q_0)T_q I(\tau),
\]

where `q_start` is the continuous kinematic value at `t_s,q`; it is unrelated
to the species charge `q_s` defined in Section 2.4. Afterward, the default
`constant-final-rate` branch continues from the transition endpoint at
`qdot_1`. Acceleration is zero at both endpoints. Axis rates are expansion
rates; only `d_c_dot` is a center-translation speed. Unequal start times and
durations allow radial and lateral overexpansion without conflating them.

A separate `tabulated-snapshots` evolution may carry time, center distance,
three semiaxes, covariance, and provenance while retaining the fixed direction
and tilt from the parent record. It uses a declared shape-preserving `C1`
interpolant with analytic derivative and rejects temporal extrapolation. A
snapshot that changes orientation is invalid in schema 5. Raw two-frame
differencing is not used to determine shock status. A future rotating branch
requires normalized quaternion interpolation, angular velocity, analytic
`Q_dot`, and separate tests. A drag-based history may be joined only in the
separately validated heliospheric branch.

A future deflecting/rotating history must carry an independently smooth center
vector and an `SO(3)` attitude history. The center interpolation must provide
analytic `c_dot`; the attitude interpolation must remain orthonormal, preserve
unit quaternion norm when quaternions are used, provide a continuous angular
velocity, and produce analytic `Q_dot`. Componentwise interpolation of raw
quaternion components is not sufficient unless followed by a method whose
normalization and derivative are derived and tested. Event sensitivity ranges
must come from the reconstruction covariance or an explicitly identified
exploratory prior; a universal `+/-10 degree` deflection or rotation is not a
physical uncertainty model.

### 8.5 Piston surface

A second, nested ellipsoid may represent the CME piston/ejecta. It uses the
same typed geometry record but an independent instance and history. The
candidate shock front must remain outside the piston on every common
heliocentric ray; the minimum standoff and its location are reported.

For the first release, the piston is either disabled or an independently
specified ellipsoid. A piston inferred inward from the shock is diagnostic
only. A future `standoff-from-piston` shock authority must define its closure,
coefficients, units, and range of validity and must reject a simultaneously
independent shock geometry. No ambiguous derived-standoff mode is permitted.

---

## 9. Shock existence, classification, and jump conditions

### 9.1 Local fast-Mach number

At every candidate surface point,

\[
 u_{1n}^{\rm in}=V_{\rm sh,n}-\mathbf u_1\cdot\hat{\mathbf n}
\]

is the positive inflow speed toward the shock. The obliquity is

\[
 \cos\theta_{Bn}=
 \frac{|\mathbf B_1\cdot\hat{\mathbf n}|}{|\mathbf B_1|}.
\]

The fast Mach number is

\[
 M_f=\frac{u_{1n}^{\rm in}}{c_f}.
\]

The reported Alfvén Mach number uses the same positive normal inflow speed and
the Alfvén speed based on the total upstream field,

\[
 M_A=\frac{u_{1n}^{\rm in}}{v_A}.
\]

If the jump solver also uses the normal-field Alfvén Mach number, it is stored
under the distinct name
`M_An=u_1n_in/(abs(B_n)/sqrt(mu_0*rho_1))`. Every critical-Mach asset declares
which convention it uses; `M_A`, `M_An`, and `M_f` are never interchangeable.
At exact `B_n=0`, `M_An` is an unbounded mathematical limit represented by a
typed state rather than floating-point infinity. A normal-Alfvén criticality
asset must define and pass a verified asymptotic treatment there or the patch
is `criticality-inapplicable`. That state is diagnostic-only for a fast-only
source; a supercritical source routes it through the selected
`fail-preflight` or `exclude-source-budgeted` policy and the same ledgers as
any other table-domain miss. It is never silently encoded as subcritical.

A point is not a fast shock when `M_f <= 1`. The code must report it as a
sub-fast candidate front, not force `M_f`, compression, or injection to a
convenient value.

### 9.2 Ideal-MHD Rankine--Hugoniot conditions

For a stationary discontinuity in the local shock frame, the jump conditions
are

**Mass:**

\[
 [\rho v_n]=0.
\]

**Normal magnetic field:**

\[
 [B_n]=0.
\]

**Tangential electric field:**

\[
 [v_n\mathbf B_t-B_n\mathbf v_t]=0.
\]

**Normal momentum:**

\[
 \left[\rho v_n^2+p+\frac{B_t^2}{2\mu_0}\right]=0.
\]

**Tangential momentum:**

\[
 \left[\rho v_n\mathbf v_t
       -\frac{B_n\mathbf B_t}{\mu_0}\right]=0.
\]

**Total energy flux:**

\[
 \left[
 v_n\left(
 \frac12\rho v^2+\frac{\gamma_{\rm ad}}{\gamma_{\rm ad}-1}p+
 \frac{B^2}{\mu_0}
 \right)
 -\frac{B_n(\mathbf v\cdot\mathbf B)}{\mu_0}
 \right]=0.
\]

These relations define candidate downstream branches from the upstream state,
shock normal, and shock speed. They have been used with coronal
white-light and spectroscopic measurements to infer CME-shock plasma and field
changes [11, 12].

### 9.3 Numerical jump solver

The implementation should solve for the compression ratio

\[
 X=\frac{\rho_2}{\rho_1}
  =\frac{v_{1n}}{v_{2n}}
\]

and the remaining downstream components using a bracketed nonlinear solve. It
must not select a root solely because it is numerically closest to an initial
guess. A candidate root is accepted only when:

1. every state variable is finite;
2. `rho_2 > rho_1 > 0`;
3. `p_2 > 0`;
4. the entropy increases across the shock;
5. all normalized Rankine--Hugoniot residuals meet tolerance;
6. the characteristic ordering corresponds to a fast shock; and
7. the compression does not exceed the physically allowed strong-shock limit
   for the selected equation of state.

For an ideal gas, the entropy diagnostic can be evaluated, up to an additive
constant, as

\[
 s_2-s_1=c_v\ln\left[
 \frac{p_2}{p_1}
 \left(\frac{\rho_1}{\rho_2}\right)^{\gamma_{\rm ad}}
 \right]>0.
\]

The corresponding ideal-gas strong hydrodynamic limit is

\[
 X_{\max}=\frac{\gamma_{\rm ad}+1}{\gamma_{\rm ad}-1}.
\]

For `gamma_ad=5/3`, this limit is four. Four is a limit, not a universal
CME-shock compression. Neither the downstream state nor this jump calculation
uses the effective wind index `gamma_w`.

At exact parallel incidence, the tangential jump system is degenerate and a
switch-on fast branch may coexist with the parallel branch over a restricted
parameter range. The solver continues the full oblique solution toward
`theta_Bn=0`, applies characteristic/evolutionary ordering, and records the
selected branch. A near-parallel angle is a numerical conditioning tolerance,
not a physical branch-selection rule, and results must converge as it is
reduced. At exactly `B_t1=0`, the azimuth of a nonzero switch-on `B_t2` is
physically degenerate; scalar outputs must be invariant under rotation of the
arbitrary tangential basis [43, 44]. The requirement not to choose a switch-on
solution from `M_A` plus an angle cutoff alone is an implementation inference:
the full jump and characteristic conditions, not those two scalars, define the
branch.

### 9.4 Critical Mach number

An MHD fast shock can exist for `M_f > 1`, but efficient ion reflection and
acceleration may require a supercritical shock. The first critical Mach number

\[
 M_c^\chi=M_c^\chi(\beta_1,\theta_{Bn},\gamma_{\rm ad}),
 \qquad
 \beta_1=\frac{2\mu_0p_1}{|\mathbf B_1|^2}.
\]

The first critical Mach number depends on upstream beta, obliquity, EOS, and
Mach-number convention. For `gamma_ad=5/3`, the Edmiston--Kennel surface spans
values near one in appropriate quasi-parallel/high-beta regimes and approaches
approximately `2.76` in the cold, quasi-perpendicular limit [13]. No fixed
“one-to-two” shortcut is permitted: beta and obliquity are evaluated locally,
and it is not valid to assume that the entire low corona has `beta<0.1` [41].

The built-in table has a version and checksum and declares the `gamma_ad` and
Mach convention `chi` for which it was generated. The code evaluates the
matching numerator `M_chi` (`M_f`, `M_A`, or `M_An`) and the table returns
`M_c^chi`; conversion among conventions is forbidden. EOS or Mach-convention
mismatch is always fatal. Interpolation coordinates and endpoint policy are
fingerprinted; the code never silently clamps or extrapolates.

Table-domain handling depends on whether criticality changes the physical
source. Before publication, a complete-history preflight evaluates the
`(beta_1,theta_Bn)` envelope over all topology-valid front patches and all
event-refined times that can contribute to the source. A production table is
required to cover the complete physical obliquity interval `[0,pi/2]`; a table
that does not is an incomplete asset and is fatal. The policy below applies to
beta-domain misses and to a verified exact-normal-Alfven asymptotic that
returns inapplicable, not to a missing obliquity domain, EOS mismatch, failed
checksum, or malformed table.

If `source.surface_requirement=fast`, criticality is diagnostic and a
beta-domain miss may return typed `criticality-inapplicable`; its patch area,
incident number/kinetic-energy flux, and time interval are reported without
changing the fast source. If `surface_requirement=supercritical`, a miss either
fails transactional preparation or, under the explicit
`exclude-source-budgeted` policy, assigns unrenormalized zero source to that
patch.

The exclusion budgets are reproducible. At time `t`, let `E(t)` be the set of
otherwise-eligible, topology-valid fast patches before the critical-Mach
query, and let `X(t)` be the subset whose criticality is unavailable under the
budgeted policy. Define

\[
 A_E=\sum_{j\in E}A_j,\qquad A_X=\sum_{j\in X}A_j,
\]

\[
 N_E=\sum_{j\in E}\sum_s n_{s,1,j}u^{\rm in}_{1n,j}A_j,
 \qquad
 N_X=\sum_{j\in X}\sum_s n_{s,1,j}u^{\rm in}_{1n,j}A_j,
\]

\[
 K_E=\sum_{j\in E}\sum_s
 \tfrac12m_sn_{s,1,j}(u^{\rm in}_{1n,j})^3A_j,
 \qquad K_X=\sum_{j\in X}\sum_s
 \tfrac12m_sn_{s,1,j}(u^{\rm in}_{1n,j})^3A_j.
\]

The species sum covers every compiled charged source species using the same
upstream abundance authority as injection. For `Y` in `{A,N,K}`, the budget
measure is

\[
 b_Y=\max\left[
 \sup_{t:Y_E(t)>0}\frac{Y_X(t)}{Y_E(t)},
 \frac{\int Y_X(t)dt}{\int Y_E(t)dt}
 \right].
\]

A ratio over an empty instantaneous support is defined as zero; if the
complete-history denominator is zero, the integrated ratio is also zero and
the run reports `no-candidate-support`, not a successful criticality coverage
claim. Area, number-flux, and kinetic-energy-flux budgets are independent.
At every time and in the time-integrated ledger,
`candidate query measure = valid table-query measure + unavailable excluded
measure`; the later split into supercritical and subcritical is stored
separately. Neighboring source weights are never renormalized to replace the
excluded physical flux. A coverage miss beyond any budget remains fatal.

The model therefore distinguishes value from validity:

- `active_fast_shock`: an admissible fast shock exists;
- `active_supercritical_shock`: a typed `{value,validity,reason}` whose Boolean
  value is meaningful only when `M_c^chi` is valid; unavailable is never
  encoded as `false`; and
- `particle_source_active`: the configured injection criterion is satisfied.

The source may require only a fast shock for baseline studies or a
supercritical shock for a more restrictive injection model. The selection is
part of the physics fingerprint.

An optional weak-field guard never edits the background or clamps a Mach
number. With declared reference field `B_ref` and dimensionless threshold
`epsilon_B`, its exact predicate is

\[
 |\mathbf B_1|<\epsilon_B B_{\rm ref}.
\]

The `diagnostic` policy only reports the predicate; `exclude-source` assigns a
typed zero source and accounts for the excluded area/rate. Both inputs are
strictly positive and sensitivity-tested.

### 9.5 Delayed, local shock activation

The prescribed ellipsoid is a candidate front before local shock admissibility
has been established. For stable patch `i`, define

\[
 H_i(t)=V_{{\rm sh},n,i}
 -\mathbf u_1(\mathbf x_i,t)\cdot\hat{\mathbf n}_i-c_{f,i}.
\]

A patch is fast-active only on the open set `H_i>0` and when the fast-branch
Rankine--Hugoniot solve is admissible; `H_i=0` is the inactive event surface.
The registered root-time and residual tolerances bound event-location error but
do not shift the physical predicate or create an unreported Mach-number
hysteresis. Source
activity is a separate state that additionally applies the selected
criticality, topology, weak-field, and source-envelope criteria.

For every patch, store the first fast-active and first source-active times and
radii, plus a list of half-open active intervals. A patch can activate,
deactivate, and reactivate as the front and background evolve. The event finder
brackets each status change between provider states. It root-locates every
continuous margin that can change the state: `M_f-1`,
`M_chi-M_c^chi`, jump-solver admissibility margins, the weak-field margin, and
the taper endpoints. It splits exactly at analytic topology/interface
crossings, kinematic knots, source-support-radius crossings, and the global
termination latch. When a jump branch loses existence without a usable smooth
margin, adaptive bisection localizes the first typed state transition to the
same event-time tolerance. Output history cadence does not define a physics
event and may not shift activation, criticality, or injection.

The recommended baseline initial requirement is `none`, so a front at
`1.01--1.05 R_sun` may initialize with zero active shock area and zero
particles. Optional already-formed-shock studies may require
`any-fast-patch` or `minimum-fast-area-fraction`. Both gates refer to locally
admissible fast-shock patches, not merely source-active or geometric patches.
For region `R`,

\[
 f_{{\rm fast},R}(t_0)=
 \frac{\sum_i A_i\chi_{R,i}\chi_{{\rm fast},i}}
      {\sum_i A_i\chi_{R,i}}.
\]

The denominator includes every physical clipped-dome patch; missing upstream
authority is fatal and cannot be removed from it. `R` is either the whole
clipped dome or a precisely defined apex cone. On failure, diagnostics report
normal-speed deficits

\[
 \Delta V_{n,i}=\max[0,c_{f,i}+
   \mathbf u_1\cdot\hat{\mathbf n}_i-V_{{\rm sh},n,i}],
\]

their area-weighted quantiles, Mach statistics, and the failed locations. A
unique minimum pair of radial and lateral speeds is not reported because the
inverse problem is underdetermined unless the input defines a one-parameter
multiplier of a fixed surface-velocity field.

### 9.6 Interpretation of a front beginning at `1.01--1.05 R_sun`

Observed metric type-II onset heights are commonly above this range. The
`1.20--1.93 R_sun` interval, with mean and median near `1.43` and
`1.38 R_sun`, is the range reported for one 32-event sample [14]; it is not a
universal prior and must not be applied as an event-independent acceptance
window. For the 2012 May 17 GLE specifically, the published radio/EUV analysis
places shock formation near `1.38 R_sun`, with its stated density-model and
emission-lane uncertainty [21]. Consequently, a candidate front initialized at
`1.01 R_sun` represents either an extreme event, a controlled parameter study,
or a prescribed disturbance below the observable onset. The model must not
present that initialization height as a typical shock-formation height.
Type-II onset remains an observational proxy affected by radio visibility,
fundamental/harmonic identification, projection, and density-model selection,
not an exact universal definition of shock formation.

---

## 10. From the candidate front to an SEP source

### 10.1 Surface discretization

The ellipsoid is tessellated into patches with known physical area. The patch
quadrature must converge in total surface area and in integrals of normal,
Mach number, and source rate. Surface resolution is independent of AMR mesh
resolution.

Each patch publishes at least:

- unique physical identity and generation;
- position and area;
- outward normal and curvature;
- normal shock speed;
- one-sided upstream authority or a structured rejection reason;
- downstream state when an admissible jump exists;
- `theta_Bn`, `M_f`, `M_A`, plasma beta, and compression;
- candidate-front, fast, jump-admissible, supercritical, and source-active
  flags;
- shock-branch and near-parallel-degeneracy identifiers;
- first-fast/first-source time and radius plus current activity-interval ID;
- topology, magnetic sector, HCS/separatrix child identity, and source
  eligibility reason; and
- the exact configuration/background generation used to derive it.

Candidate ellipsoid elements below the solar boundary may appear only in the
pre-clipping geometry ledger; they are not physical surface patches. Every
physical patch satisfies `r>=R_sun`. Physical patches below the qualified
source radius, outside the active AMPS domain, or failing shock admissibility
carry diagnostics but receive zero injection weight.

Closed-field patches are diagnosed, but the recommended source policy is
`diagnose-only`. Enabling closed-field injection requires a validated
closed-loop mover, explicit bidirectional scattering, two absorbing solar
footpoints, mirroring/loss-cone tests, and trapped/precipitating ledgers. The
shock tessellation is geometrically split at every HCS or open/closed
separatrix intersection so that no patch averages opposite sectors or plasma
authorities.

Patches inside the resolved PFSS/SCS transition are tagged separately because
their upstream current and force balance are a model regularization. The
source-region policy is either `exclude` (retain shock diagnostics but inject
zero) or `convergence-qualified`. The latter is legal only after an ensemble
of transition widths/profiles demonstrates convergence of `M_f`, obliquity,
compression, injected rate, and observer products; its contribution remains a
separate ledger category. Independently, the Section 6.4 `S_tr` clearance mask
has unconditional precedence in schema 5: any intersecting patch receives
exactly zero source even when the broader transition region is
`convergence-qualified`. Masked and region-policy rates close separate
ledgers, are not redistributed, and are never silently combined with PFSS or
SCS patches.

### 10.2 Test-particle DSA spectrum

For a steady, planar, test-particle shock with isotropic diffusion, the
phase-space distribution has

\[
  f(p)\propto p^{-q_{\rm DSA}},
  \qquad q_{\rm DSA}=\frac{3X}{X-1}.
\]

The density per unit momentum is

\[
  \frac{dN}{dp}\propto p^2f(p)\propto p^{2-q_{\rm DSA}}.
\]

The implementation must preserve this Jacobian when sampling momentum. This
spectrum is a closure, not a prediction of the injection efficiency, maximum
energy, or finite shock-age rollover [15, 16].

In the first-release moving-source model, this DSA slope is used only as a
**phenomenological conditional momentum spectrum of the net first-passage
flux through the finite upstream reference surface**. The steady shock
distribution `f(p)` is not, by
itself, an upstream escape spectrum. The model does not resolve downstream
acceleration, residence, or a free-escape boundary, so it must not describe
`local-compression-dsa` as a derived total acceleration or escape solution.

### 10.3 Source normalization

A production source is defined for every compiled AMPS species by a stable
semantic ID plus a verified application-local compiled slot, never by assuming
that index zero is a proton. After AMPS publishes `PIC::nTotalSpecies`, the
adapter reads `GetChemSymbol`, mass, and signed charge for each slot and checks
them against exactly one `[species.ID]` record. Distinct compiled populations
may share a chemical symbol, but their stable IDs and compiled slots remain
different; slot numbers are never exchanged between independently built
applications. A source-enabled production run then requires exactly one
`[source_species.ID]` record with the same suffix for every `[species.ID]`.
No compiled entry may be silently skipped. The schema-5 SEP movers transport
charged particles, so a source-enabled run rejects every compiled neutral
species before allocation; `[species.ID].transport_role=initialization-only`
is legal only in a source-disabled initialization/verification workflow.
Missing slots, repeated slots, duplicate stable IDs, mismatched symbol/mass/
charge, or unsupported neutrals are fatal; the runtime file never changes the
compiled species table.

The source boundary is part of the physical model, not a particle-placement
detail. Schema 5 has three noninterchangeable branches:

1. `upstream-reference-first-passage` is the preferred event baseline. For
   each shock patch it constructs the kinematically moving upstream reference surface
   `Sigma_ref,j(t)` at the positive signed-front distance
   `d_sh=L_ref>0`. `L_ref` is a finite, observation/model-calibration-owned SI
   length, not a cell fraction. It must be smaller than the local normal reach
   of the front, must not produce an offset-surface fold or patch overlap, and
   must remain inside the active domain for the complete source interval. The
   supplied release rate is the one-way outward first-hitting rate through this
   surface; "net first passage" names the shock-to-surface survivor population,
   not an outward-minus-inward Eulerian flux evaluated after later recrossings.
   Unresolved acceleration and all return/recycling between the shock
   and `Sigma_ref` are already included in that calibrated rate.
2. `front-conormal-flux-no-through-flow-sensitivity` applies a Parker
   phase-space flux directly at the candidate front. Its stochastic boundary
   map reflects a returning diffusive trajectory in the diffusion-tensor
   conormal direction, not in the geometric normal. It is a named numerical
   sensitivity closure, requires `source.scientific_role=sensitivity`, and
   must never be described as physical shock reflection. Schema 5 rejects
   this branch for a focused mover because a
   gyro-averaged `(p,mu)` state does not contain the gyrophase needed for a
   general oblique specular boundary map.
3. `shock-adjacent-absorbing-verification` is the former positive-`delta`
   placement layer with an absorbing front. It exists only to demonstrate the
   source/transport degeneracy, requires
   `source.scientific_role=verification`, and is forbidden in event-nominal
   and observation-comparison products.

Schema 5 deliberately leaves the physical choice of `L_ref` to a documented
calibration or uncertainty ensemble. It provides geometric admissibility and
convergence requirements, but it does not contain an observational inversion
that determines `L_ref`, and it does not select that distance by holding a
Péclet number fixed. The production baseline uses one declared physical
distance for the applicable release-surface record and holds it fixed while
mesh and source-kernel thickness are refined. Any result sensitive to that
choice must propagate the calibrated range of `L_ref`; geometric convergence
alone does not validate the release distance.

Where a diffusion representation is valid, the schema-5 D10 contract requires
the dimensionless
**release-depth diagnostic** for stable species `a`, momentum `p`, front
surface coordinate `sigma`, and time `t`,

\[
 {\cal P}_a(p,\sigma,t)=
 \int_{0}^{L_{\rm ref}}
 \frac{u^{\rm in}_{1n}(d,\sigma,t)}
      {\kappa_{nn,a}(d,p,\sigma,t)}\,dd,
 \qquad
 \kappa_{nn}=\hat{\mathbf n}\mathbin{\cdot}
             \boldsymbol\kappa\mathbin{\cdot}\hat{\mathbf n},
\]

evaluated on the same registered upstream normal ray used by the source.
Normal depth `d` increases from the front into the upstream region and
`u_1n^in>0` denotes the front-frame relative flow magnitude toward the front.
The Parker branch uses the actual transport tensor. The focused branch may
report only an explicitly labeled diffusion-approximation proxy formed from
`kappa_parallel=v*lambda_parallel/3` and its independently declared
`kappa_perp`; that proxy does not alter the focused mover and is not a second
surface-only transport coefficient. The integral is applicable only where the
sign convention and finite `kappa_nn>0` are certified over the complete
interval. For constant coefficients,

\[
 {\cal P}=\frac{u^{\rm in}_{1n}L_{\rm ref}}{\kappa_{nn}}
          =\frac{L_{\rm ref}}{\ell_d},
 \qquad \ell_d=\frac{\kappa_{nn}}{u^{\rm in}_{1n}}.
\]

Thus `P<1` places the surface *inside* one conventional upstream
diffusion length, whereas `P>1` places it farther than one such length
in the constant planar limit. Neither statement is a survival or return
probability. Indeed, on an unbounded one-dimensional half-line with constant
coefficients and drift toward an absorbing front, the eventual front-hitting
probability is unity in the absence of an outer escape boundary or another
killing process, irrespective of this finite starting-depth label. Pitch-angle
memory, anisotropic and spatially varying diffusion,
front curvature and motion, finite source duration, outer escape, and the
chosen evaluation horizon all affect return statistics. A branch with no
declared spatial-diffusion interpretation reports a typed inapplicable state;
it never invents a scalar `kappa_nn` merely to emit this diagnostic. The
directly measured, finite-horizon
`DelayedFrontReturn` and `SurvivingUpstreamInventory(t)` ledgers remain the
authoritative transport diagnostics.

For a cohort evaluated at a declared horizon `T`, the represented-number
quantity

\[
 F_{{\rm no\ front\ return},c}(T)=
 1-\frac{N_{{\rm DelayedFrontReturn},c}(\le T)}
         {N_{{\rm committed},c}(\le T)}
\]

may be reported only as the **finite-horizon no-front-return fraction**. It is
distinct from the live upstream-inventory fraction because particles lost at
the Sun, transition clearance, or outer boundary have not returned to the
front but are no longer live. It is neither an infinite-time escape
probability nor a quantity derivable from `P` alone; the cohort definition,
birth interval, horizon, and all competing terminal ledgers must accompany it.
A zero committed denominator is inapplicable, not a numerical zero. When
`T` is an administrative end time, members born at different times have
different exposure ages: the displayed quantity is then a crude cumulative
ledger fraction, not a survival estimator. Output must retain birth time,
exposure age, right-censoring state, and competing-loss cause for every cohort.
An age-conditioned return curve at lag `tau` uses only cohorts with the
declared follow-up or a preregistered competing-risk/survival estimator and
publishes its risk set and confidence interval; it must not relabel the
administrative-horizon ratio as asymptotic survival.

A future fixed-`P` rule with one common finite positive target would generally make
`L_a=L_a(p,sigma,t)` because `kappa_nn` is species-, momentum-, position-, and
time-dependent. It would therefore define a family of release surfaces and a
placement measure `varphi_ref,a(x,p,t)`, invalidate the present
spatial/momentum factorization, and require species/momentum/time-resolved
reach, offset-cap, fold, self/other-surface-overlap, interface-clearance,
active-domain, generation, interpolation, and ledger tests. Such a rule is a
new source model, not a change of default for schema 5. A target that itself
varies by species, momentum, patch, or time is a still more general calibrated
source law and requires a separate discriminator and inference contract.

For the baseline, let `g_s,j(p,t)` be the normalized momentum density of the
net first-passage population in the declared source frame,

\[
 \int_{p_{\min,s}}^{p_{\max,s}}g_{s,j}(p,t)\,dp=1,
\]

and let `varphi_ref,j(x,t)` be a nonnegative **one-sided approximate surface
delta** normalized over physical volume,

\[
 \int \varphi_{{\rm ref},j}(\mathbf x,t)\,d^3x=1.
\]

Writing `xi_ref=d_sh-L_ref`, its support satisfies
`0<=xi_ref<=epsilon_ref`; hence every birth is on the outward side of the
first-passage surface and `d_sh>=L_ref>0`. Its first normal moment satisfies

\[
 0\le \int \xi_{\rm ref}\varphi_{{\rm ref},j}\,d^3x
 \le \epsilon_{\rm ref}\longrightarrow0.
\]

`epsilon_ref` may scale with local mesh spacing only to realize the one-sided
surface delta numerically. Refinement sends it and the first moment to zero
while **holding `L_ref` fixed**. The complete support remains separated from
the front for every ensemble member, so convergence never sends the source
toward the absorbing front.

For a focused source, a momentum-independent `h(mu)` is generally invalid.
Let `a_s,j(p,mu,t)` be a nonnegative seed angular density normalized over
`[-1,1]`. Let `n_ref` be the outward normal, let `V_ref,n` be the unique normal
speed of the moving reference surface, and define the canonical normal
surface-velocity representative
`V_ref=V_ref,n*n_ref`; no physical result depends on an arbitrary tangential
surface parametrization. For every proposed `(p,mu)` the implementation first
applies the complete supported source-frame transformation from Section 2.4,
including a selected parallel wave-frame boost, and obtains the physical
mover-resolved guiding-center particle velocity `v_particle^T` in the declared
transport frame. The authoritative relative outward speed is

\[
 w_{n,s,j}(p,\mu,t)=
 [\mathbf v_{{\rm particle},s}^{T}(p,\mu,t)-\mathbf V_{\rm ref}]
 \mathbin{\cdot}\hat{\mathbf n}_{\rm ref}.
\]

Thus neither relativistic particle velocity and plasma velocity nor a wave
velocity are added Galileanly. The frame transformation occurs **before** the
normal-flux test and before `Z` is evaluated.

The conditional first-passage angular density is the normal-flux-weighted
law

\[
 h_{s,j}(\mu\mid p,t)=
 \frac{a_{s,j}(p,\mu,t)[w_{n,s,j}(p,\mu,t)]_+}
 {Z_{s,j}(p,t)},\qquad
 Z_{s,j}=\int_{-1}^{1}a_{s,j}[w_{n,s,j}]_+\,d\mu,
\]

when `Z_s,j>0`. Direct evaluation of this integral is authoritative; the state
is `NoFocusedEscape` if and only if `Z_s,j=0`. Thus the actual joint law is
`H_s,j(p,mu,t)=g_s,j(p,t)h_s,j(mu|p,t)` and is not silently refactored as a
momentum-independent pitch-angle law. The seed `a=1/2` is called isotropic;
the **released flux** is not isotropic after multiplication by `[w_n]_+`.
Gyrophase is uniform on `[0,2 pi)` only after this conditional law is accepted.

The exact zero decision is a support statement, not a floating-point cutoff.
Define

\[
 {cal S}(p,t)=\{\mu\in[-1,1]:a(p,\mu,t)>0
                  \ \hbox{and}\ w_n(p,\mu,t)>0\}.
\]

Because the integrand is nonnegative, `Z=0` exactly when this set has zero
measure. The numerical evaluator partitions `[-1,1]` at seed-support edges,
pitch-law knots, and every isolated root of the **unclipped** `w_n=0`, certifies
the sign on each open subinterval, and integrates only positive intervals.
Root-bracket uncertainty is included in the quadrature error. For `Z>0`, both
configured absolute and relative error bounds must pass; neither bound is a
zero test. A certified positive-measure interval remains admissible even when
its `Z` is smaller than the absolute error target. Failure to isolate a root or
certify the integral is `UnresolvedPositiveFluxNumerics`, which is fatal before
allocation and cannot debit `NoFocusedEscapeExcluded`. The implementation
never tests `Z<=epsilon`, floors `Z`, clamps a negative estimate to zero, or
redistributes unresolved measure.

If support changes within a configured momentum or event-time interval, the
preflight recursively brackets and splits that interval at the support boundary
using the registered momentum/time root tolerances **before** rate integration
and ledgering. The exact tangent condition remains
`b_hat dot n_hat=0`. The configured minimum `abs(b_hat dot n_hat)` controls only
whether the derived `mu_c` diagnostic is well conditioned and emitted; near-
tangent states still use the same transformed direct `w_n` and `Z`, so crossing
that diagnostic threshold cannot change escape status or the sampled law.

The Parker source has no pitch coordinate and uses its net reference-surface
flux directly.

For interpretation and an analytic sampler test only, consider the
nonrelativistic upstream-plasma-frame limit
`v_particle^T=u_1+v_s(p)*mu*b_hat`. If
`beta_n=abs(b_hat dot n_ref)>0`, write
`sigma=sign(b_hat dot n_ref)` and
`Delta u_n=V_ref,n-u_1 dot n_ref`. Outward focused release then requires

\[
 \sigma\mu>\mu_c(p,j,t),\qquad
 \mu_c=\frac{\Delta u_n}
 {v_s(p)\beta_n}.
\]

The signed `sigma` is mandatory: using only
`abs(b_hat dot n_ref)` reverses the allowed hemisphere in one magnetic sector.
The three exact support regimes are: `mu_c<=-1`, for which all pitch support is
outward apart from a possible measure-zero endpoint; `-1<mu_c<1`, for which
only `sigma*mu>mu_c` is outward; and `mu_c>=1`, for which no pitch is outward.
When `beta_n=0`, `sigma` and `mu_c` are undefined and are not evaluated. In
this tangent-field limit the same nonrelativistic expression gives
`w_n=u_1 dot n_ref-V_ref,n`, so positive relative bulk advection admits the
complete seed distribution and nonpositive relative advection admits none.
The general transformed-velocity `Z` test, not this simplified case table,
decides every production state. A zero `Z` is a typed `NoFocusedEscape`, not a
zero hidden inside an otherwise valid patch.
That status is carried per stable species ID, momentum interval, patch,
provider generation, and event-split time interval. Production either fails
preflight or applies the explicitly selected budgeted exclusion; it never
renormalizes the rejected momentum measure onto the surviving support.

The flux-fraction branch defines the physical patch rate

\[
 \dot N_{s,j}=
 H_{{\rm act},j}H_{\rm term}\,
 \eta_{s,j}\,n_{s,1,j}u^{\rm in}_{1n,j}A_j,
 \qquad 0\le\eta_{s,j}\le1,
\]

Here `eta_s,j` is the empirical fraction of the incident upstream species flux
that appears as **net outward first passage through `Sigma_ref`**. It is not
the thermal-to-nonthermal injection efficiency at the shock, the total number
accelerated on both sides, a downstream DSA normalization, or the shock's
total acceleration efficiency. Likewise, `g_s,j(p,t)` is conditional on that
finite-surface release and integrates to unity over the modeled release
interval. Consequently, an event-nominal focused deck must have admissible
outward phase space everywhere that `dot N*g` is positive; otherwise it fails
preflight rather than redefining the declared net flux. Only a labeled
sensitivity may interpret the positive but inaccessible measure as a
candidate demand and debit it to `NoFocusedEscapeExcluded` without
renormalization. The reference distance, surface generation, release law, and
population meaning are part of the configuration fingerprint and every
physical-rate ledger.

The corresponding **number-rate density per unit momentum** is

\[
 S_{N,s}(\mathbf x,p,t)=
 \sum_j \dot N_{s,j}\,g_{s,j}(p,t)\,
 \varphi_{{\rm ref},j}(\mathbf x,t).
\]

It integrates to the physical rate but is not yet the source term for a
phase-space density. With the conventions used in Section 11,

\[
 n_s=4\pi\int_0^\infty p^2 f_s^{P}\,dp
\]

for the isotropic Parker distribution and

\[
 n_s=2\pi\int_0^\infty p^2dp\int_{-1}^{1}f_s^{F}\,d\mu
\]

for the gyrotropic focused distribution. Therefore the correctly normalized
phase-space source terms are

\[
 Q_{f,s}^{P}=\frac{S_{N,s}}{4\pi p^2},
 \qquad
 Q_{f,s}^{F}=\sum_j
 \frac{\dot N_{s,j}\varphi_{{\rm ref},j}(\mathbf x,t)
       H_{s,j}(p,\mu,t)}{2\pi p^2}.
\]

Their momentum and angular integrals recover
`sum_j dot(N_s,j)` exactly. The displayed factored formulas are written in the
frame in which `g` and `h` are prescribed. For a supported wave-frame focused
source, use primed `(p',mu')` in that formula and push the complete measure
forward to the **local-plasma mover frame** with the exact parallel
four-momentum transformation and momentum-space Jacobian. The resulting
mover-frame distribution remains gyrotropic but need not remain separable in
`p` and `mu`; code must not refactor it. Direct Monte Carlo sampling in the
source frame followed by the same wave-to-plasma boost is the equivalent
change of variables. Shock-frame and Cartesian storage-frame momenta are
computed afterward as independent per-particle diagnostic/representation
views and do not define this source term.

Thus patch area occurs exactly once; an area-normalized surface delta and an
additional `A_j` factor cannot be combined. In the preferred moving-source
mode, `varphi_ref,j` approaches a delta on the **fixed physical reference
surface**, rather than a delta on the shock. A convergence sequence changes
kernel width and mesh scale but not `L_ref`, the reference-surface identity,
or the physical rate.

The `physical-rate` branch instead supplies a per-species total rate from a
versioned asset or explicit scalar and distributes it over eligible patches
with nonnegative weights that sum to one. It does not also apply `eta`. In
either branch, `n_s,1` comes from a background species density or an explicit
abundance relative to a stable background species ID; it is never inferred from
total mass density. Density and temperature alone do not determine a
macroparticle birth rate.
`H_act,j` contains the local fast/jump/criticality/topology decision, while
`H_term` is the independently selected radial source envelope.

Because the baseline has no downstream sheath, define signed front distance
`d_sh>0` on the outward/upstream side of the oriented candidate front,
`d_sh=0` on it, and `d_sh<0` on the unmodeled CME-interior/downstream side.
The reference surface is the level set `d_sh=L_ref`; its numerical source
kernel has support wholly inside `d_sh>0` and remains separated from the front
under the complete mesh/refinement ensemble. A baseline particle that later
reaches `d_sh=0` after a resolved positive-distance flight is event-located,
removed, and recorded in the schema-stable `DelayedFrontReturn` ledger. This is a transport
loss **after** the calibrated first-passage release, not a correction to
`eta*g`. Candidate/fast/supercritical state at return is recorded because the
domain boundary exists over the complete candidate front even when a patch is
temporarily sub-fast. The particle is never propagated into a fictitious
upstream state on the other side and is never silently reinjected.

In particular, a returned particle must not simply be re-emitted with its
incident momentum or with a new draw from `g(p,t)`. Either choice invents an
accelerator while omitting downstream transmission, residence time,
return probability, energy gain, pitch-angle transformation, and the source
of the added particle energy; it can also double count shock recycling already
contained in the calibrated first-passage rate `eta*g`. A future return model
must therefore be one of two independently validated constructions:

1. an explicit downstream/sheath accelerator with resolved shock crossings
   and work; or
2. a normalized return-to-release renewal kernel conditioned on incoming
   species, momentum, pitch/full velocity as required, patch, time, and front
   state, with complementary absorption/transmission probabilities.

The renewal construction must specify its residence-time distribution,
momentum and angular transition law, frame transformations, support, and
termination. It must close represented number, momentum, kinetic energy, and
shock-work ledgers, preserve immutable cohort ancestry through every renewal,
and prove that its calibrated normalization does not already include the same
recycling in `eta*g`. Until one of these constructions is implemented and
validated, `DelayedFrontReturn` remains a terminal loss in the preferred
branch.

The conormal sensitivity branch instead imposes, for the Parker equation, the
moving-surface normal phase-space flux

\[
 J_n=\left[(\mathbf u-\mathbf V_{\rm sh})f
           -\boldsymbol\kappa\!\cdot\!\nabla f\right]
           \!\cdot\!\hat{\mathbf n}=J_{n,{\rm src}}(p,t)
\]

and prevents through-flow on stochastic returns. For anisotropic diffusion the
reflection direction is proportional to
`kappa dot n/(n dot kappa dot n)`; geometric-normal reflection is valid only
for isotropic `kappa`. This construction requires the explicitly checked
normal diffusivity

\[
 \kappa_{nn}=\hat{\mathbf n}\mathbin{\cdot}
 \boldsymbol\kappa\mathbin{\cdot}\hat{\mathbf n}>0
\]

at every boundary quadrature point, with a registered positive numerical
margin. If `kappa_nn` is zero or sub-margin, the conormal direction is
undefined. Schema 5 records `degenerate-normal-diffusion` and rejects that
patch/branch; it does not divide by a floor, substitute geometric reflection,
or silently reinterpret a purely advective boundary. A future advective
inflow/outflow branch must derive and test its own characteristic boundary
condition. The qualified boundary must reproduce the imposed conormal flux
without an auxiliary volume source, and number conservation at the boundary
must close independently of the outer-domain loss ledger. This sensitivity
branch has no placement distance and no `DelayedFrontReturn` debit.

Return detection uses the immutable front and reference-surface generations
that bracket the particle segment. After a provider generation change, only
the unadvanced remainder is tested against the new generation, and a stable
`(particle,event,generation,boundary-model)` key makes each event exactly once
under MPI repartitioning. Coincident absorbing events retain all reason bits
but debit one primary ledger in the deterministic order solar boundary,
transition clearance, source-boundary/front return, then outer escape.
Baseline source birth is committed only after the kernel is proven disjoint
from the front. In the verification branch a birth/front coincidence is
classified as `ImmediateShockAdjacentReturn`, not hidden by an epsilon shift.

Every branch publishes additive number and kinetic-energy ledgers, plus the
four-momentum sum in the ledger's declared inertial frame, indexed by stable
species ID, momentum bin, patch ID, event-split time interval, provider
generation, and boundary model. The minimum source-boundary entries are
`CandidateReferenceRelease`, `NoFocusedEscapeExcluded`,
`CommittedFirstPassageRelease`, `ImmediateShockAdjacentReturn`,
`DelayedFrontReturn`, and `NetFirstPassageRelease`. The exact source identities
are

\[
 L_{\rm candidate}=L_{\rm no\ escape}+L_{\rm committed},
 \qquad
 L_{\rm net\ first\ passage}=L_{\rm committed},
\]

with no redistribution. At any ledger time `t`, the separately named
`SurvivingUpstreamInventory(t)` closes committed represented weight against
live particles and every terminal physical loss accumulated by that time; it
obeys, for each additive represented measure,

\[
 L_{\rm surviving}(t)=L_{\rm committed}(\le t)
 -\sum_{\ell\in{\rm terminal}}L_\ell(\le t),
\]

and is never relabelled as a first-passage rate. `ImmediateShockAdjacentReturn` is
identically zero in the preferred and conormal branches. The absorbing
verification branch alone publishes `GrossShockAdjacentEmission`,
`HitEscapeSurfaceBeforeShock`, and immediate/delayed return ledgers over its
declared diagnostic horizon. Those verification-only quantities do not change
the definition of `NetFirstPassageRelease` in the preferred branch. Counts,
represented weights, momentum, and energy are never mixed in one scalar
fraction.

Runtime loss acceptance uses represented physical measures, never raw
macroparticle counts. A cohort is the immutable tuple
`(stable species ID, birth-energy bin, stable birth-patch lineage,
birth-time interval)`, where both edge sets and the patch-lineage map come from
the preregistered consumer-budget asset. For cohort `c`,

\[
 N_{{\rm birth},c}=\sum_{i\in c}w_i,
 \qquad
 E_{{\rm birth},c}=\sum_{i\in c}w_iK_i(t_{{\rm birth},i}),
\]

are the explicit denominators. A loss channel `l` accumulates

\[
 N_{l,c}=\sum_{i\in l,c}w_i,\qquad
 E^{\rm birth}_{l,c}=\sum_{i\in l,c}w_iK_i(t_{{\rm birth},i}),\qquad
 E^{\rm event}_{l,c}=\sum_{i\in l,c}w_iK_i(t_{l,i}).
\]

The bounded loss fractions are `N_l,c/N_birth,c` and
`E_birth_l,c/E_birth,c`. Both lie in `[0,1]` for an exactly-once terminal loss
ledger and measure the represented number and represented **birth energy** of
the cohort that was lost. The separate event-energy impact ratio
`E_event_l,c/E_birth,c` is finite and nonnegative but deliberately is not
bounded by one: adiabatic, electric-field, or other retained work can change
particle kinetic energy before the event. The manifest reports all three
dimensional numerators, their denominators, and the unbounded impact ratio, so
energy evolution cannot be hidden or mislabelled as a loss fraction.
Split/merge descendants inherit the cohort exactly, and ledgers sum their
physical weights; changing
macroparticle multiplicity cannot change a budget.

`TransitionSheetExcludedLoss` is capped for every cohort and registered
aggregate. `DelayedFrontReturn` is capped only for a branch that
actually removes a particle at the front. A no-through-flow or other
nonabsorbing branch sets every front-loss asset and maximum inactive and keeps
front-contact/recrossing as a nonloss diagnostic; it must not manufacture a
zero loss fraction. A retained absorbing verification/sensitivity branch must
activate the front-return budget. A zero birth denominator is typed
`inapplicable-empty-cohort`, while a nonzero loss numerator with a zero
denominator is an invariant failure. Any exceeded cap makes an event run
`not-event-grade`; the physical ledgers remain available and are never
renormalized into the surviving population.

A future accelerator-derived `free-escape-boundary` release closure must solve
and normalize an
actual upstream diffusion-advection problem. For the restricted steady planar
case with constant upstream `u_1`, diffusion coefficient `kappa_1(p)`, and a
free-escape surface a distance `L_esc` ahead of the shock, the escaping
number-flux spectrum is proportional to

\[
 4\pi p^2 f_0(p)\,
 \frac{u_1}{\exp[u_1L_{\rm esc}/\kappa_1(p)]-1}.
\]

That branch must specify its pitch-angle distribution, finite-age/maximum-
momentum behavior, and whether particles are emitted at the escape surface or
at the shock; it cannot multiply the shock spectrum by an unspecified escape
probability or also use the shock-adjacent placement kernel [52]. It is a
future schema capability, not an alternative interpretation of the empirical
finite-reference-surface baseline.

The test-particle approximation has two mandatory budget gates. First, either
normalization branch verifies per patch and over the whole surface that the
injected number rate does not exceed the incoming upstream number flux of that
species. Second,

\[
 \sum_s \dot N_{s,j}\langle E_{k,s}\rangle
 \le \epsilon_{\rm nt}\,A_j F_{{\rm conv},j},
 \qquad 0<\epsilon_{\rm nt}\le\epsilon_{\rm tp}<1,
\]

where the baseline convertible-flux ceiling is computed from the same
shock-frame Rankine--Hugoniot state,

\[
 F_{{\rm conv},j}=\max\left[0,
 \frac{\rho_1u_{1n}^{\rm in}}{2}
 \left(|\mathbf v_1|^2-|\mathbf v_2|^2\right)\right].
\]

This is the loss of bulk kinetic-energy flux, not the conserved total MHD
energy flux and not a claim that all converted energy is available to SEPs. A
versioned event-specific asset may prescribe a smaller available flux. The
release profile supplies the preregistered test-particle perturbativity limit
`epsilon_tp`; the requested efficiency `epsilon_nt` may not exceed it. Neither
quantity has an undocumented universal default, and both are fingerprinted
and sensitivity-tested. The mean kinetic energy is evaluated after the sampled
distribution has been transformed into that same shock frame; mixing a
plasma-frame spectrum energy with a shock-frame flux is forbidden. This is a
test-particle ceiling, not an acceleration-efficiency prediction. A deck that
exceeds either number or energy budget is rejected; it is not repaired by
changing statistical weight.

The physical particle rate and numerical samples per step are separate:

\[
 W_{s,j,n}=\frac{1}{N_{{\rm macro},s,j,n}}
 \int_{t_n}^{t_{n+1}}\dot N_{s,j}(t)\,dt.
\]

This is the required represented physical weight of each newly allocated
macroparticle in that species/patch allocation group. Let the positive
configured AMPS base weight installed for species `s` be `W_base,s`. The
created particle receives

\[
 C_{s,j,n}=W_{s,j,n}/W_{{\rm base},s}
\]

as its weight-correction factor, so the authoritative represented weight is
the product `W_base,s*C_s,j,n=W_s,j,n`. The base value is a numerical scale, not a
second source normalization. Split/merge may redistribute this product among
model particles but must preserve its sum. Output reports the product (and its
logarithm only when positive), not either factor alone. The sum of represented
physical particles must close to the declared physical rate on every step and
across every MPI decomposition.

`N_macro,s,j,n` is positive only for a group with a strictly positive
integrated physical rate. If the event-split integral is zero, no particle is
allocated and no zero-weight placeholder is created; `samples_per_step` is an
allocation budget for active groups, not an instruction to manufacture
particles during a zero-source interval.

The time integral is split at every activation/deactivation, criticality,
topology, source-radius, kinematic-knot, taper, and termination event from
Section 9.5. Within each smooth subinterval it is evaluated analytically when
possible or by a registered error-controlled quadrature. Replacing it by the
endpoint rate times the full step is legal only when the rate is proven
constant on that step.

### 10.4 Injection direction and position

The preferred branch places particles in the one-sided `varphi_ref,j`
approximate delta immediately outward of the exact finite reference surface,
with `d_sh>=L_ref`; it does not straddle that surface. The surface offset is fixed in meters throughout a
placement-convergence sequence. Only the numerical kernel thickness is an
explicit fraction of local mesh scale with an absolute cap, and its complete
support must lie outward of the reference surface and remain separated from
the front. The offset surface,
kernel support, normal, velocity, reach margin, and active-domain clearance
are preflighted at every provider knot; interpolation between knots may not
change surface topology or cross the front.

Momentum direction is selected explicitly. The schema-5 Parker mover has no
pitch coordinate and consumes the scalar net first-passage rate in the
upstream-plasma frame. A focused mover samples the joint
`H(p,mu,t)=g(p,t)h(mu|p,t)` of Section 10.3; an unrestricted isotropic angular
birth is legal only in the shock-adjacent absorbing verification branch. The
focused baseline may use the upstream-plasma frame or a named parallel wave
frame. The outward-wave frame moves at `+v_A*t_hat` relative to the local
plasma and requires `w_out>0`; the inward-wave frame moves at `-v_A*t_hat` and
requires `w_in>0`. The seed law is sampled in its declared frame, but every
candidate is transformed through the supported exact frame chain before the
transport-frame `w_n` test and positive-flux weighting are evaluated. Thus a
wave-frame branch includes the wave velocity before acceptance; it cannot test
with wave-frame `p,mu` while omitting motion of that frame.

The implementation applies only the corresponding parallel wave-to-plasma
Lorentz boost to obtain mover `p,mu` and rejects any state with an invalid MHD
wave speed. Because that boost is parallel to `b_hat`, a gyrotropic law remains
gyrotropic but generally does not remain separable in momentum and pitch. A
subsequent plasma-to-transport transformation is used to construct the
physical guiding-center velocity required by the boundary test, but never to
redefine the local-plasma gyrotropic mover coordinates. If a proposed
nonparallel source-frame transformation would require unresolved gyrophase to
decide crossing, that source branch is unsupported and rejected rather than
gyro-averaged implicitly. A reconstructed full per-particle four-vector may be
transformed separately for AMPS storage or diagnostics after the declared
gyrophase sampling. A general shock-frame source boost likewise has a
perpendicular component and is not selectable by the gyrotropic schema-5
mover; it requires a future full-velocity or explicitly gyro-averaged source
operator. The exact wave-to-plasma transformation in Section 2.4 is applied
before the particle is handed to the mover. No default angular law, magnetic
orientation, or implicit frame is inferred from the spectrum name.

The absorbing verification branch deliberately retains shock-adjacent
placement `delta`. In a steady one-dimensional drift--diffusion manufactured
problem, let `x` increase outward from the absorbing front at zero to the
escape surface at `L`, and define `u_1>0` as the magnitude of drift **toward**
the front, so the stochastic drift is `dx=-u_1 dt+sqrt(2*kappa)dW`. Then

\[
 P_{\rm esc}(\delta)=
 \frac{\exp(u_1\delta/\kappa)-1}
      {\exp(u_1L/\kappa)-1}\longrightarrow0
 \quad\text{as}\quad\delta\longrightarrow0.
\]

Its acceptance test is to reproduce this degeneracy. It is not allowed to
pass a production placement-convergence gate by holding `delta` fixed in cell
units or by renormalizing survivors.

### 10.5 Maximum energy

The baseline may accept a prescribed energy range for controlled comparisons.
A later age-limited model should determine `p_max` by comparing acceleration,
escape, and loss times. Using the shock-normal projection of each one-sided
diffusion tensor, a representative acceleration time is

\[
 t_{\rm acc}(p)=\frac{3}{u_1-u_2}
 \left(\frac{\kappa_1(p)}{u_1}+
       \frac{\kappa_2(p)}{u_2}\right),
\]

Here `u_1>u_2>0` are the magnitudes of the upstream and downstream
shock-frame normal flow speeds, not inertial bulk-speed vectors, and
`kappa_i=n_hat dot kappa_i dot n_hat` is the corresponding normal diffusion
coefficient. It equals the parallel coefficient only for the appropriate
parallel geometry. The estimate is invalid until both one-sided diffusion
tensors and their obliquity treatment are defined. Self-generated waves and
their back reaction on coronal shock acceleration are
a distinct advanced closure rather than an automatic consequence of this
timescale estimate [24, 25].

### 10.6 Radial source termination

The candidate-front geometry and diagnostic shock history may continue after
particle injection ends. A hard cutoff remains available for regression. The
optional quintic taper uses

\[
 H_{\rm term}=\begin{cases}
 1, & r_{\rm apex}\le r_1,\\
 1-S\!\left(\dfrac{r_{\rm apex}-r_1}{r_0-r_1}\right),
    & r_1<r_{\rm apex}<r_0,\\
 0, & r_{\rm apex}\ge r_0,
 \end{cases}
\]

where `S` is the quintic smoothstep from Section 8.4 and `r_1<r_0`. `r_0` is
the exact zero-source radius; it is not the end of the geometry history. The
taper is an empirical model-form regularization, not a physical repair for an
invalid shock. Local patch admissibility remains mandatory throughout.
Production candidate-front histories require a nondecreasing apex over the
source interval. In addition, reaching `r_0` sets a persistent
`source_envelope_terminated` latch, so interpolation noise or a later
diagnostic inward motion cannot restart injection.

Output and ledgers retain the untapered physical rate, `H_term`, realized rate,
and fraction of the total source emitted inside the taper. The envelope scales
represented physical weight/rate and never silently changes the configured
Monte Carlo sample count. Results must converge with particle step, provider
cadence, and taper width; the zero-width limit recovers the hard cutoff.

### 10.7 Optional future flare-associated attribution source

The schema-5 baseline contains only the CME-shock-associated source defined in
Sections 10.2--10.6. It can test whether that source and the selected transport
are sufficient for an event, but it cannot infer from agreement alone that a
flare-associated contribution was absent. Conversely, an early observed onset
that the shock-only model misses is not by itself proof of flare acceleration;
magnetic connectivity, release-surface calibration, scattering, instrumental
response, and background uncertainty must first be propagated.

A future flare-associated branch is therefore an optional **source-attribution
sensitivity**, not a correction silently added to the shock baseline. It must
have an independent space--time support, stable source identity, species and
composition, momentum and angular law, source frame, normalization,
uncertainty, provenance, and additive number/momentum/energy ledgers. Hard
X-ray, gamma-ray, type-III radio, EUV, or reconnection-flux observations may
constrain timing, location, or a prior, but none is automatically a calibrated
escaping-ion rate or spectrum. Any proxy-to-particle conversion must be an
explicit, versioned inference model with uncertainties and withheld tests.

When shock and flare-associated sources coexist, their immutable source labels
must survive splitting, merging, transport, and observer sampling so that the
mixed prediction can be decomposed without rerunning or assigning particles
after arrival. A joint fit must expose degeneracy between source amplitude,
spectrum, release time, connectivity, and transport coefficients; it may not
reuse the same observation independently to construct both sources and then
claim it as validation.

The future branch must also define the relationship between every impulsive
birth and the moving candidate front. If a legal birth precedes the first
published front generation, it carries an explicit no-front-yet state until a
front exists; it is not guessed to be upstream. A birth in the unmodeled downstream/CME
interior is invalid unless Stage 11B or another independently validated
downstream provider supplies that state. If an upstream impulsive particle is
later overtaken by the front, the upstream-only approximation may absorb it
into a distinct `ImpulsiveFrontEncounter` terminal interaction ledger while
retaining `ImpulsiveCoronalRelease` as its immutable origin. Transmission or
reacceleration is permitted only through a validated downstream-transfer or
return-renewal operator that preserves the particle's impulsive ancestry,
records residence and shock work, and keeps it out of the independently
calibrated `eta*g(p)` first-passage source. Thus the same physical particle is
never counted once as flare release and again as a new shock-source birth.
Until this operator and its validation campaign exist, schema 5 reports a
shock-only result and identifies flare attribution as unresolved model-form
uncertainty.

---

## 11. SEP transport equations

### 11.1 Parker transport

For an isotropic distribution `f(x,p,t)`, the Parker equation is written in
the conventional transport form [17]:

\[
 \frac{\partial f}{\partial t}
 +\mathbf u\cdot\nabla f
 =\nabla\cdot(\boldsymbol\kappa\cdot\nabla f)
 +\frac{1}{3}(\nabla\cdot\mathbf u)
 p\frac{\partial f}{\partial p}+Q.
\]

Here `Q` is the phase-space-density source `Q_f^P` from Section 10.3, not the
number-rate density `S_N`.

The time derivative is taken at fixed coordinates in the declared transport
frame. Let `v_grid` be that frame's grid velocity relative to the inertial
heliocentric frame: `v_grid=0` for an inertial run and
`v_grid=Omega_F cross x` for the supported rigidly corotating run. Everywhere
in Sections 11.1--11.4,

\[
 \mathbf u\equiv\mathbf u_{\rm inertial}-\mathbf v_{\rm grid}.
\]

The phase-space density is a scalar under the coordinate rotation, `p` and
`mu` remain defined in the local plasma frame, and the tensor and spatial
derivatives are rotated consistently. The placement density in `Q` follows
the coordinate-volume mapping, but its momentum variables are already in the
local-plasma mover frame established in Sections 2.4 and 10.3; no
plasma-to-transport momentum boost is inserted into the Parker/focused
equation. A source or observer may not use inertial positions with a
corotating time derivative.
The first production profile uses either the inertial frame or the single
rigid sidereal frame normalized in Section 2.1; latitude-dependent rotation is
verification-only.

The diffusion tensor is

\[
 \boldsymbol\kappa=
 \kappa_\parallel\hat{\mathbf b}\hat{\mathbf b}
 +\kappa_\perp(\mathbf I-\hat{\mathbf b}\hat{\mathbf b})
 +\boldsymbol\kappa_A,
\]

where the antisymmetric term represents drift when enabled.

Schema 5 sets `kappa_A=0` for both Parker and focused production transport.
A weak-scattering gradient/curvature drift requires a complete relativistic
coefficient, reduction factor, weak-field regularization, and explicit
HCS/separatrix/interface behavior; none is yet part of this baseline. Naming a
future coefficient asset does not make that operator executable.

A future drift capability must be derived as one frame-consistent transport
operator rather than added as an independent velocity correction. For the
isotropic Parker equation, one admissible formulation is an explicitly signed
antisymmetric diffusion tensor whose weak-scattering coefficient has magnitude

\[
 \kappa_A^{(0)}=\frac{p v}{3|q_s|\,|\mathbf B|},
\]

with the charge sign carried by the tensor convention and with any turbulence
reduction factor stated and validated. The associated drift advection is
obtained from the divergence of that same tensor (equivalently from the curl
of its signed coefficient times `b_hat`, with the selected sign convention);
the code may not specify an unrelated drift velocity and tensor.

For focused transport, the spatial gradient/curvature drift and the
corresponding momentum and pitch characteristics must be reduced from the same
underlying relativistic guiding-center kinetic equation in the declared
transport frame. In an inertial full-orbit description the kinetic-energy work is
`q_s E dot v_gc`; in the present local-plasma formulation, however, the
existing `dot(p)` already represents work associated with large-scale plasma
compression and strain after a frame reduction. Appending
`q_s E dot v_d` to that `dot(p)` without re-deriving the complete transformed
characteristics can count the same electric-field work twice. The future
operator must demonstrate the energy identity of its chosen variables and
frame and replace the affected characteristic terms coherently; an additive
`v_d`/`qE dot v_d` patch is forbidden.

Both formulations require explicit guiding-center validity and numerical
guards: nonzero resolved `|B|`; gyroradius and gyroperiod small relative to the
registered field/flow variation scales and times; differentiable coefficients
with gradient convergence; bounded reduction factors; and a time step that
resolves drift across the local cell and coefficient scales. Zero-thickness
HCS, open/closed separatrices, the PFSS/SCS interface, and the transition
clearance are not ordinary smooth-gradient locations. Crossing or sheet drift
there requires a separately derived finite-thickness or discontinuity
operator [47, 48]; otherwise the drift branch must exclude or terminate before
the interface according to a budgeted policy. Parker antisymmetric diffusion
and focused guiding-center characteristics are equation-specific reductions
of one coherent future provider, not independently adjustable physics. Each
reduction nevertheless requires its own conservation, polarity-reversal,
weak-field, manufactured-solution/orbit, and convergence tests.
Analytic Parker-spiral drift and full-orbit SEP studies provide independent
benchmarks for the smooth-field limit [57, 58].

A future runtime drift-significance product must retain both signed/net drift
and accumulated pathwise drift exposure; cancellation in the former cannot
hide large drift excursions in the latter. Its stochastic comparison scale
must state its dimension. For locally isotropic perpendicular diffusion, the
one-component standard deviation is
`sigma_perp,1=sqrt(2*integral(kappa_perp dt))`, whereas the two-dimensional
transverse RMS radius is
`rho_perp,rms=sqrt(4*integral(kappa_perp dt))`. The implementation must name
which denominator it uses and may not interchange them. An anisotropic
perpendicular tensor instead publishes the integrated covariance and its
principal-axis widths.

### 11.2 Focused transport

For gyrotropic `f(x,p,mu,t)`, a convenient characteristic form follows the
standard focused-transport construction [18, 19]:

\[
 \frac{\partial f}{\partial t}
 +(\mathbf u+\mu v\hat{\mathbf b})\cdot\nabla f
 +\dot p\frac{\partial f}{\partial p}
 +\dot\mu\frac{\partial f}{\partial\mu}
 =\frac{\partial}{\partial\mu}
 \left(D_{\mu\mu}\frac{\partial f}{\partial\mu}\right)
 +\nabla\cdot(\boldsymbol\kappa_\perp\nabla f)+Q.
\]

Here `Q=Q_f^F` under the gyrotropic normalization in Section 10.3.
For the baseline isotropic cross-field closure,

\[
 \boldsymbol\kappa_\perp=\kappa_\perp
 (\mathbf I-\hat{\mathbf b}\hat{\mathbf b});
\]

in a future anisotropic
closure the time-step controller uses the trace/eigenvalues of the declared
positive-semidefinite perpendicular tensor rather than treating it as a
scalar.

Neglecting relativistic inertial corrections to the solar-wind frame, the
large-scale characteristic rates are

\[
 \dot p=-p\left[
 \frac{1-\mu^2}{2}\nabla\cdot\mathbf u
 +\frac{3\mu^2-1}{2}
 \hat{\mathbf b}\hat{\mathbf b}:\nabla\mathbf u
 \right],
\]

\[
 \dot\mu=\frac{1-\mu^2}{2}\left[
 v\nabla\cdot\hat{\mathbf b}
 +\mu\left(\nabla\cdot\mathbf u
 -3\hat{\mathbf b}\hat{\mathbf b}:\nabla\mathbf u\right)
 \right].
\]

The background provider must therefore supply or consistently derive
`div(b_hat)`, focusing length, `div(u)`, and field-aligned strain. The Parker
and focused movers must use the same background snapshot and shock-source
generation when they are compared [17--19].

The input distinguishes two physically different focused collision operators.
`focused-pitch-angle-diffusion` realizes the Fokker--Planck term with the Itô
increment

\[
 d\mu=\left[\dot\mu+\frac{\partial D_{\mu\mu}}{\partial\mu}\right]dt
       +\sqrt{2D_{\mu\mu}}\,dW_t,
\]

using a verified no-flux treatment at `mu=+/-1`. Its mean free path is

\[
 \lambda_\parallel=\frac{3v}{8}
 \int_{-1}^{1}\frac{(1-\mu^2)^2}{D_{\mu\mu}(\mu)}\,d\mu.
\]

`focused-discrete-scattering` instead replaces the diffusion term by a
normalized, detailed-balance collision kernel and samples Poisson scattering
events. It is not silently identified with pitch-angle diffusion; an input
mean free path determines its event rate only for a named kernel with a tested
velocity-autocorrelation relation. The small-angle sequence must converge to
the selected `D_mumu` operator before the two modes are compared. For Parker
transport, the corresponding diffusion-limit closure is
`kappa_parallel=v*lambda_parallel/3`; `kappa_perp` and antisymmetric drift have
independent, explicitly selected models.

The first-release `isotropic` pitch-angle-diffusion shape is executable rather
than merely named:

\[
 D_{\mu\mu}^{\rm iso}=\frac{v}{2\lambda_\parallel}(1-\mu^2).
\]

Substitution in the integral above returns exactly `lambda_parallel`. The
first-release `isotropic-poisson` kernel redraws `mu` uniformly on `[-1,1]`
and gyrophase uniformly on `[0,2 pi)` at Poisson rate

\[
 \nu=\frac{v}{\lambda_\parallel},
 \qquad \langle\mathbf v(0)\cdot\mathbf v(t)\rangle=v^2e^{-\nu t},
\]

so its diffusion limit is `kappa_parallel=v*lambda_parallel/3`. Other shapes
or kernels require a versioned normalization asset. The first-release focused
equation contains no gradient/curvature-drift characteristic, so focused
transport requires `drift_model=none`; an antisymmetric Parker drift model
cannot be reused by changing an enum.

### 11.3 Particle and background time-step constraints

The stochastic transport equations do not acquire a conventional explicit
diffusion stability limit when they are integrated as stochastic differential
equations, but their coefficients must remain approximately local over one
step. Define `Delta` as the local active-cell scale and `L_C` as the minimum
resolved coefficient-variation length. A conservative accuracy controller is

\[
 \Delta t \le C_t\min\left(
 \frac{\Delta}{|\mathbf u|+v},
 \frac{\Delta^2}{2(\kappa_\parallel+2\kappa_\perp)},
 \frac{L_B}{v},
 \frac{L_C}{|\mathbf u|+v},
 \frac{L_{\rm sh}}{|V_{\rm sh,n}|}
 \right),
\]

where `L_B=|1/(b_hat dot grad ln|B|)|`, `L_sh` is the smallest shock-geometry
or source-weight variation scale, and `0<C_t<1` is a declared accuracy factor.
The diffusion term is an accuracy bound on the root-mean-square displacement,
not a claim of Eulerian CFL instability.

For pitch-angle diffusion, the adaptive substep must additionally resolve the
local scattering time. One useful local condition is

\[
 \Delta t \le C_\mu
 \frac{(1-|\mu|+\epsilon_\mu)^2}
      {2D_{\mu\mu}(\mu,p,\mathbf x)},
\]

with the endpoint regularization and reflecting-boundary scheme defined by the
selected focused mover. Here `0<C_mu<1` and `0<epsilon_mu<1`; `epsilon_mu` is
the dimensionless endpoint floor supplied by
`pitch_endpoint_regularization`. The code must use the mover's mathematically verified
endpoint treatment rather than relying on this bound to prevent `|mu|>1`.

The background/shock publication cadence must also satisfy

\[
 \Delta t_{\rm shock}
 \max_{i\in\mathcal S}
 \frac{|V_{{\rm sh},n,i}|}
      {\min(\Delta_i,L_{{\rm patch},i})}
 \le C_{\rm sh},
 \qquad 0<C_{\rm sh}<1,
\]

unless exact continuous-time surface evaluation is used. `L_patch` is the
local patch scale. Time-step convergence is part of the validation campaign;
the constants `C_t`, `C_mu`, and `C_sh` are numerical controls, not physical
parameters.

### 11.4 Field-aligned reduction for `srcSEP`

For an exported open line `l`, let `s` be outward arc length, `A_l(s,t)` the
represented flux-tube area, and `x_l(s,t)` its position in the inertial frame.
At fixed line coordinate define

\[
 \mathbf v_{\rm line}=\left.\frac{\partial\mathbf x_l}{\partial t}\right|_s,
 \qquad
 u_s=(\mathbf u_{\rm inertial}-\mathbf v_{\rm line})
       \cdot\hat{\mathbf t},
 \qquad
 \kappa_{ss}=\hat{\mathbf t}\cdot\boldsymbol\kappa\cdot\hat{\mathbf t}.
\]

For the rigidly rotating analytic line,
`v_line=Omega_F cross x`; in coordinates corotating at that same rate the line
is stationary and `u_s=u dot t_hat` using the transport-frame velocity already
defined in Section 11.1. The bundle stores `v_line`, the transport-frame ID,
and `u_s`, so neither application subtracts rotation twice. In the equation
below, the time derivative is taken at fixed `s`, each spatial derivative at
fixed `t`, and `Q_l` is the source transformed into its represented tube
measure.

When cross-field diffusion and drift are disabled, the 3-D isotropic Parker
equation reduces on that tube to

\[
 \frac{\partial f}{\partial t}
 +u_s\frac{\partial f}{\partial s}
 =\frac{1}{A_l}\frac{\partial}{\partial s}
   \left(A_l\kappa_{ss}\frac{\partial f}{\partial s}\right)
 +\frac{1}{3}(\nabla\cdot\mathbf u)
 p\frac{\partial f}{\partial p}+Q_l.
\]

For a characteristic line this is the infinitesimal-tube equation. For a
finite footprint it is a cross-section-integrated closure only when `f` and
the transported coefficients are uniform to the registered transverse error;
otherwise the footprint must be subdivided into a quadrature group and the
solutions combined with the derived flux weights. A single averaged line is
not claimed to reproduce arbitrary transverse structure merely because
perpendicular diffusion is zero.

The corresponding focused equation uses the same `dot(p)` and `dot(mu)` from
Section 11.2, evaluated from the exported `div(u)`,
`b_hat b_hat : grad(u)`, `div(b_hat)`, and magnetic polarity. `srcSEP` must not
recompute these derivatives from an unrelated Parker formula. This is
particularly important at `R_i` and `R_scs`: `R_i` may contain a typed vector
jump, while the field is continuous at `R_scs` but its one-sided derivatives
may differ. A focused line containing an unresolved interface or HCS crossing
is rejected rather than interpolated through it.

Pitch-angle cosine `mu` is measured relative to `b_hat`, not the outward
tangent. Consequently the field-aligned streaming characteristic in the
outward coordinate is

\[
 \frac{ds}{dt}=u_s+\sigma_{B,l}\mu v.
\]

The polarity-independent outward derivative and the signed focusing
quantities are related by

\[
 \hat{\mathbf b}\cdot\nabla\ln|B|
 =\sigma_{B,l}\frac{d\ln|B|}{ds},
 \qquad
 \nabla\cdot\hat{\mathbf b}
 =-\sigma_{B,l}\frac{d\ln|B|}{ds},
 \qquad
 L=\frac{1}{\nabla\cdot\hat{\mathbf b}}.
\]

The infinite-`L` uniform-field limit is represented by a validity/state flag,
not a floating-point infinity in the exchange table.

Magnetic flux is the conserved measure of a tube. For a finite cross-section
`C_l(s,t)`,

\[
 \Phi_l=\int_{C_l(s,t)}\mathbf B\cdot d\mathbf S,
 \qquad
 A_l(s,t)=\int_{C_l(s,t)}dA.
\]

Only in the infinitesimal, transversely uniform limit does this reduce to
`A_l|B_l|=Phi_l` using the centerline field. A scalar `Phi_l` and a centerline
therefore determine a thin-tube area scale but **not** a finite tube's shape,
observer overlap, or shock-patch overlap.

Every finite-measure first-release line uses a
`traced-reference-footprint`: a checksummed, oriented, non-self-intersecting
triangulated patch transverse to the field at `trace_time_s`. Its boundary and
interior quadrature vertices are traced with the same one-sided field
authority as the centerline, producing a nonfolded swept tube map
`X_l(s,xi_1,xi_2,t)` and a positive volume Jacobian `J_l`. The exporter stores
the cross-section boundary/connectivity or an equivalent lossless swept-cell
representation at every required node. It verifies the surface-flux integral,
positive Jacobian, geometric closure, and nonoverlap with every other tube in
the same measure group. The represented volume of any subset `D` is then

\[
 V_l(D,t)=\int\!\!\int\!\!\int
  \mathbf 1_D\!\left(\mathbf X_l(s,\boldsymbol\xi,t)\right)
  J_l(s,\boldsymbol\xi,t)\,ds\,d^2\boldsymbol\xi.
\]

This is evaluated by conservative intersection of swept tube cells with `D`,
not by testing only the centerline. The export must carry positive `Phi_l`
whenever absolute source normalization or count-based sampling is requested.
A characteristic line has no footprint or absolute flux and supports only
relative or per-unit-flux quantities. A one-line thin-tube approximation must
be identified as asymptotic and converge under footprint subdivision; it is
not presented as an exact reduction of a transversely structured finite tube.

For an infinitesimal flux tube crossing the shock transversely, a local shock
area associated with the tube is

\[
 dA_{\rm sh}=\frac{d\Phi}{|\mathbf B\cdot\hat{\mathbf n}_{\rm sh}|}.
\]

This relation is the infinitesimal transverse-crossing limit and is retained
as a verification identity. It becomes ill-conditioned at a tangency. An
absolute finite-tube export instead clips the physical shock patch against the
explicit swept tube and integrates the source density over that overlap;
surface area and tube-assigned rate close conservatively across neighboring
nonoverlapping footprints. It must not cap the infinitesimal denominator or
infer a finite overlap from the centerline. Characteristic-line calculations
may still report the local spectrum at a tangent root, but not an absolute
tube-integrated source.

For a zero-width characteristic line, the one-dimensional source location is
obtained from the exact centerline intersection

\[
 F(\mathbf x_l(s_{l,k}(t),t),t)=0,
\]

where `k` indexes all roots. For a finite tube, the authoritative geometry is
instead every connected component of

\[
 F(\mathbf X_l(s,\boldsymbol\xi,t),t)=0
\]

inside the swept footprint. The exporter conservatively clips that surface,
projects its distributed physical source onto the 1-D `s` measure, and retains
a component even when the centerline has no root. It must retain zero, one, or
multiple roots/components rather than select the nearest one. Each record
exports its `s` support/centroid, surface/component identity, area, and the
subelement quadrature records for normal, normal speed, `theta_Bn`, Mach
numbers, compression, spectrum, activity, and physical rate. Source-weighted
means may be emitted as diagnostics, but the injected line source remains the
mixture of subelement spectra rather than a spectrum reconstructed from mean
shock parameters. Time continuation preserves each identity over a topology-
constant interval and event-locates its creation, split, merge, annihilation,
or departure from the tube.

A geometric intersection record and an injecting source are independent
facts. Every record stores its source-envelope domain status, local shock/source eligibility, detailed reason,
all concurrent rejection flags plus an optional primary reason, independent
source-measure status, source-envelope factor, diagnostic-only flag, assigned physical rate, and
activity-interval identity. A line-level connection summary retains counts and
first time/radius for intersections inside the source domain, source-active
records, valid but inactive records, invalid upstream/jump states, tangency-measure
rejections, and post-termination diagnostic records. `No intersection within the
history horizon` and `history incomplete` remain distinct; neither is silently
reported as a permanent physical disconnection.

For a set of field lines, each line is an independent transport corridor
unless the input explicitly defines a quadrature over a finite bundle of
magnetic flux. Copying the complete 3-D shock rate onto every exported line
would multiply the physical source and is forbidden. A multi-line ensemble
must use one of the following declared measures:

- `characteristic`: each line reports intensity per unit area or per unit
  magnetic flux and is not summed as a global population;
- `flux-tube`: every line owns a positive `Phi_l`, an explicit nonoverlapping
  traced footprint, and its conservatively clipped fraction of the shock
  source; or
- `quadrature`: each line owns one nonoverlapping traced footprint cell and
  positive `Phi_l` within a named finite measure group. The group flux is
  `Phi_group=sum_l Phi_l`, and the derived, non-input quadrature weight is
  `Phi_l/Phi_group`; the footprints, fluxes, and source rates all close over
  the same group domain.

A `characteristic` line has no absolute `Phi_l`, `A_l`, or represented volume.
It may report the local distribution or intensity per unit magnetic flux, but
it cannot use a count-based residence/crossing observer estimator. Any
observer accumulated from weighted model particles requires a `flux-tube` or
`quadrature` measure with positive absolute flux and explicit finite geometry.
Dimensionless quadrature weights alone never define a volume.

Perpendicular diffusion is not representable by independent 1-D lines. It is
disabled for strict 1-D/3-D parity. A future coupled multi-line model may add
an explicit, conservative exchange operator between neighboring lines, but it
must not reinterpret independent line hopping as the 3-D perpendicular
diffusion operator.

### 11.5 Required contents of a field-line bundle

Bundle schema 3 is the first schema that carries lossless finite-footprint and
time-resolved overlap geometry; an older area-only bundle cannot be promoted
implicitly. The exchange bundle is immutable and self-describing. It contains one
manifest, one or more line records, optional time-dependent snapshots,
candidate-front intersection and shock/source-state histories, observer
mappings, and checksums. Every line node must
provide at least:

| Group | Required quantities |
|---|---|
| Geometry | stable line/node IDs, `s`, `x`, line-grid velocity, outward tangent, magnetic polarity, composite-footpoint sector, SCS-construction sector, PFSS and composite topology, magnetic region, interface side, transition-sheet state/distance/identity, source longitude `alpha`, `A_phi`, `J_phi`, fold margin, and mapping validity |
| Magnetic field | vector `B`, `B_abs`, `d ln(B_abs)/ds`, `div(b_hat)`, focusing length, signed HCS/interface distances, and discontinuity identities |
| Solar wind | vector `u`, `u_s`, mass density, pressure, temperature, `div(u)`, `b_hat b_hat:grad(u)` |
| Turbulence | `w_out`, `w_in`, direction basis, correlation length, spectral parameters, and validity flags; `w_+/-` may be derived |
| Tube measure | `Phi_l`, finite-footprint asset identity, cross-section boundary/swept-cell connectivity, volume Jacobian or conservative cell volumes, `A_l`, measure-group ID and flux, and derived quadrature weight when applicable |
| Observer mapping | stable observer ID; mapping times and topology-event times; stable overlap-component IDs; all connected support intervals; closest `s`, separation, tolerance, conservatively intersected tube/observer volume and volume exposure; detector position/velocity; estimator validity and coverage |
| Numerics | recommended coefficient-variation scale, segment length, and snapshot validity interval |
| Provenance | schema, units, frame, epoch, all active thermodynamic indices (`gamma_w`, `gamma_ad`, and optional `gamma_c`); raw rotation rate/convention, resolved sidereal rate/vector/axis, conversion-asset checksum; provider/current-sheet/topology/mapping and moving-front generations/hashes; source population/release/re-entry semantics; termination history, configuration fingerprint, and source checksums |

The manifest also carries every source application's compiled species entry:
stable semantic ID, source-local compiled slot, chemical symbol, mass, signed
charge, optional validated integer nucleon count, source normalization, and
transport role. `srcSEP` resolves the stable ID and physical identity against
its own immutable AMPS molecular table before allocation; it must not assume
that an integer slot has the same meaning in two separately compiled
executables. Duplicate stable IDs/slots, missing records, or physically
inconsistent identities are fatal, while duplicate chemical symbols are legal
for distinct roles and slots. For a matched parity run, every compiled
transport species must be present and mapped exactly once.

The source-history record serializes
`population_semantics=net-first-passage-upstream-released`,
`release_model=empirical-net-first-passage-release`, the release-boundary
discriminant, finite reference-surface distance/geometry/generation, focused
first-passage and no-escape policies, re-entry policy, all source-boundary
ledger terms, the signed-distance convention and event ordering from Section
10.3, and the candidate-front/shock generation identity used for each
interval. `srcSEP` consumes these values; its input has no local section that
can reinterpret the released population, move its reference surface, reverse
the magnetic-orientation test, or change the return policy.

The manifest contains a machine-readable data dictionary for every tabular
field: canonical name, SI unit, scalar/vector/category type, applicability,
companion validity field, and serialization sentinel. Nullable values have at
least the states `valid`, `inapplicable`, and `invalid`; more specific typed
causes may be added without collapsing those states. In memory, `NodeState`
and observer/source records carry the typed state, not a magic floating-point
value. Writers emit the finite sentinel only in the numeric file member, and
readers consult the companion state; equality comparison with the sentinel is
never an API.

Energy- or species-dependent diffusion coefficients should normally be
evaluated in both applications by the shared coefficient registry from the
exported primitive turbulence and plasma state. A tabulated coefficient may
be exported only with explicit species, energy, pitch-angle, units, and
interpolation axes. A single unlabeled `lambda_parallel` column is
insufficient.

Node coordinates are strictly increasing in `s`, begin on or outside the solar
absorbing boundary, and extend beyond every requested observer. Positive
scalars are interpolated with a positivity-preserving method; unit vectors are
renormalized after interpolation; and no interpolant may span an invalid
snapshot, `R_i`, `R_scs`, an HCS, separatrix, shock discontinuity, mapping fold,
or any derivative-side change. Every production-exported open line remains in
one unique composite-footpoint magnetic sector and maintains positive
clearance from `S_tr`. Touching the exclusion, a transition null, a typed kink,
or a composite/SCS-sector mismatch is a fatal production-export error; the
line is not truncated or moved to a nearby seed. A requested line rejected for
this reason is `topology-inapplicable`, not evidence of physical magnetic
disconnection. The manifest records the topology policy, minimum clearance,
rejection event, and D9 hash.
`srcSEP` must fail rather than extrapolate beyond either the spatial or
temporal validity range.

### 11.6 Multiple observers and energy-spectrum sampling

Observers are repeated, stable-ID records rather than a single global point.
Each record owns its frame/ephemeris, collection geometry, sampling cadence,
species selection, energy-bin edges, optional pitch-angle bins, and field-line
association. A 3-D observer samples a declared physical volume or surface; a
1-D observer uses an exported `LineObserverMapping`, never an inferred
`(line_id,s)` guess. A nearest-cell or nearest-line substitution outside the
configured tolerance is forbidden.

For the first-release 1-D residence estimator, the exporter conservatively
intersects the actual observer volume with the explicit swept finite tube. Let
`A_overlap,l,o(s,t)` denote the integral of the tube-coordinate volume
Jacobian over only that part of the transverse footprint inside observer `o`.
Then

\[
 V_{{\rm overlap},l,o}(t)=
 \int A_{{\rm overlap},l,o}(s,t)\,ds.
\]

At every mapping time the bundle stores stable IDs for each connected overlap
component, its nonzero-`s` support, closest centerline coordinate and
separation, conservative overlap volume, and geometric error bound. It never
substitutes `A_l` over the interval in which only the centerline is inside the
observer; that would count whole cross-sections at partial entry and exit. The
overlap, not a nominal center alone, defines the estimator. Schema-5 strict
1-D parity supports this volume-residence mapping. A 3-D disk or spherical
crossing observer may be retained for 3-D output, but a 1-D crossing estimator
is rejected unless a future bundle carries a separately qualified aperture--
tube overlap and geometric-factor operator.

Whenever observer and line have relative motion, the exporter evaluates a
cadence-resolved overlap table over the complete run horizon. Creation,
tangency, splitting, merging, and annihilation of overlap components are
event-located to the configured mapping tolerance. Each component has a stable
ID only on one topology-constant half-open time segment; interpolation is
bounded within that segment and forbidden across a topology event. Each valid
record is found from the physical tube/volume intersection, not nearest-node
selection. `srcSEP` samples only where the complete interpolated geometry
satisfies the connection and geometric-error tolerances. A requested run with
any uncovered sampling interval fails before allocation. If one line cannot
remain connected, the event needs multiple explicitly weighted line records or
a future topology-evolving bundle; the exporter never slides a static line to
follow the observer.

All reported energy, speed, arrival/look direction, and pitch-angle bins are
defined in the observer rest frame. An ephemeris supplies both position and
velocity; a fixed observer is stationary in its declared coordinate frame,
whose velocity relative to the transport frame is still transformed
explicitly. Before binning, the particle four-momentum is Lorentz-boosted to
the detector frame. Directional or pitch-angle sampling also transforms the
local electromagnetic field into that frame before constructing the field
direction; the analytic provider supplies `B` and the ideal-MHD electric field
`E=-u cross B` in the declared inertial frame. A plasma-frame particle speed or pitch angle is never substituted
for an instrument-frame observable.

Energy channels are defined by **edges**, not ambiguous nominal centers. The
supported constructions are linear, logarithmic, or an explicit checked file:

\[
 E_{k+1/2}=E_{\min}+k\frac{E_{\max}-E_{\min}}{N_E}
\]

or

\[
 E_{k+1/2}=E_{\min}
 \left(\frac{E_{\max}}{E_{\min}}\right)^{k/N_E},
 \qquad k=0,\ldots,N_E.
\]

The record states whether `E` is kinetic energy per particle or per nucleon.
Conversion uses the compiled mass and, for the latter, the validated integer
nucleon count from Section 2.4.

Observer cadence and accumulation intervals are measured in detector
acquisition-clock SI seconds. An ephemeris record supplies the monotone mapping
between that clock and run coordinate time; the identity map is recorded for a
fixed observer when applicable. Let `I_o` be one acquisition interval and
`V_repr,o(t)` the physical 3-D sampling volume or, for a 1-D line, the
conservative tube/observer overlap volume above. Define the volume exposure

\[
 \mathcal X_{V,o}=\int_{I_o}V_{{\rm repr},o}(t(\tau_o))\,d\tau_o.
\]

For a volume-residence estimator, the focused directional differential
intensity is

\[
 \widehat j_{kl}=
 \frac{\displaystyle\sum_a W_a
       \int_{I_o}\!\mathbf 1_{a,kl,o}(\tau_o)
       v_{{\rm det},a}(\tau_o)\,d\tau_o}
      {2\pi\mathcal X_{V,o}\,\Delta E_k\,\Delta\mu_l}
 \quad [\mathrm{m^{-2}\,s^{-1}\,sr^{-1}\,J^{-1}}].
\]

For an isotropic Parker sample, omit the pitch bin and replace the angular
measure `2*pi*Delta_mu` by `4*pi`. The familiar `V*Delta t` denominator is
used only when `V_repr,o` is proven constant over the interval. A 1-D observer
uses the identical exposure estimator with the time-resolved overlap above.

For focused bins whose union covers all of `mu in [-1,1]` without gaps and
whose gyrophase coverage is complete, the omnidirectional intensity is

\[
 \widehat J_k=2\pi\sum_l\widehat j_{kl}\Delta\mu_l
\]

for focused bins and `J_hat_k=4*pi*j_hat_k` for the explicitly isotropic Parker
distribution; its units omit steradians. A partial pitch-angle set or look cone
may instead report `accepted-solid-angle-integrated` intensity over its exact
acceptance, but it must not be labelled omnidirectional. A focused
omnidirectional value inferred from partial coverage requires a future named,
validated inversion operator. When output energy is per eV rather than per
joule, the conversion is applied to the bin width and declared in the column
metadata.

For a surface-crossing estimator, every accepted oriented crossing contributes
its physical weight once to `N_kl`. The declared aperture and look bin define

\[
 G_l=\int_A\int_{\Omega_l}
 |\hat{\mathbf v}\cdot\hat{\mathbf n}|\,d\Omega\,dA
 \quad [\mathrm{m^2\,sr}],
\]

For a time-dependent aperture, orientation, or angular acceptance, define

\[
 \mathcal X_{G,l}=\int_{I_o}G_l(\tau_o)\,d\tau_o,
 \qquad
 \widehat j_{kl}=\frac{N_{kl}}
 {\Delta E_k\,\mathcal X_{G,l}}.
\]

Only a proven constant `G_l` reduces this exposure to `G_l*Delta t`.

For a moving detector, velocity and crossing orientation are evaluated
relative to the detector surface. A disk uses its single configured normal. A
spherical surface uses the local outward normal at each crossing and a
declared inward, outward, or both-senses rule; it never reuses one global
normal. Geometry radius determines disk/sphere area, which is derived rather
than supplied independently. The geometry kernel computes and records `G_l`;
a point detector or zero-area surface is invalid. Volume residence and surface
crossing are different estimators and cannot be combined in one bin.
Instrument-response convolution is performed only after the raw physical
intensity has been stored; an online response is a separately configured,
checksummed linear operator.

Every bin carries represented particle number, effective sample count,
physical exposure/measure, value, uncertainty estimate, and validity. If no
particles contribute, count and intensity are finite zero and
`particles_present=0`; undefined ratios such as anisotropy are emitted as the
manifest sentinel `-1.7976931348623157e308` with `value_valid=0`. Empty sampling volumes and
invalid background cells can therefore be distinguished without `NaN`.

---

## 12. Boundary and initial conditions

### 12.1 Solar internal boundary

AMPS registers an internal sphere at

\[
  r=R_\odot.
\]

Particle trajectories intersecting this sphere are absorbed unless a future
explicit albedo model is selected. There is no specular reflection by default.
The magnetic and wind models are evaluated only on the physical side of the
boundary. Piston and candidate-front patches below the sphere are clipped.

### 12.2 Inner guard shell and source-support radius

There is exactly one inner particle-absorption surface: the registered sphere
`r=R_sun` in Section 12.1. `R_in` is the minimum radius at which shock-source
physics is qualified, not another boundary. Cells intersecting the guard shell
`R_sun<r<R_in` remain active and receive the complete finite PFSS, wind or
closed-plasma, turbulence, time-step, and categorical state. They permit an
inward particle to precipitate continuously to the solar sphere; they never
receive shock-source injection and are excluded from low-coronal validation
claims.

PFSS is mathematically defined from `R_sun`, and the global wind solution is
continued inward from its explicit base-reference sphere, so no constant
supersonic wind or uninitialized state is imposed in the guard shell. The
candidate-front dome may intersect the solar sphere geometrically, but its
apex at the initial epoch and every source-supporting patch must satisfy
`r>=R_in`. A proposed active source patch below `R_in` is retained as a
diagnostic candidate and assigned the typed reason `BelowQualifiedSourceRadius`.

### 12.3 Internal magnetic and topology interfaces

`R_i`, `R_scs`, the pure-SCS/Parker HCS, the transition sector surface
`S_tr`, and the open/closed separatrix are distinct internal model interfaces.
The ideal verification form at `R_i` has
continuous radial flux and the typed tangential surface-current jump from
Section 6.4; production replaces it with the resolved divergence-free overlap
shell. At `R_scs`, the SCS field and rotating-footpoint extension are fully
continuous, with one-sided derivatives where necessary. Pure-HCS and
separatrix crossings follow the mover-specific event rules in Section 6.8.
The ordinary ideal-HCS sign-flip rule is forbidden at `S_tr`; the selected
transition topology policy supplies either a qualified operator or the common
exclusion boundary. No quantity is assumed continuous unless its governing
interface condition requires continuity.

### 12.4 Outer boundary

The exact physical boundary is the Sun-centered sphere `r=R_dom`, and its
particle condition is escape. Cartesian AMR blocks wholly outside the sphere
are inactive; intersected blocks remain active and particle segments are
clipped/event-located at the exact sphere rather than a staircase of block
faces. An escaped particle is removed after its final boundary-crossing
contribution is recorded. The analytic background uses an outward
characteristic continuation up to the sphere. The Cartesian box wholly
encloses it, and every requested observer collection region lies strictly
inside it.

### 12.5 Lateral boundary of an active Parker corridor

AMPS may deactivate blocks outside a configured corridor. `full-domain` keeps
all physical blocks. `swept-field-line-tube` constructs the union of finite
tubes around the named line set over every exported geometry/source snapshot.
The corridor half-width is a constant or a tabulated function of line coordinate; an
additional geometric buffer covers field-line interpolation, shock-source
placement, and the maximum registered particle displacement used by the
qualification run.

A block is active when its closed bounding volume intersects the buffered
swept-tube volume, not when its center happens to lie inside. The mask then
adds every AMR ancestor, required child, face/edge/corner neighbor used by a
stencil, and ghost-halo owner. Finally, a graph check requires a face-connected
active path from the solar guard shell through every source-supporting
intersection to every observer collection region. This conservative geometric
closure is the normative hole-free criterion; visualization morphology is
only a diagnostic.

With `kappa_perp=0` **and** `drift_model=none`, the corridor is a computational
optimization provided no field-aligned step can leave it. A particle crossing
the lateral corridor surface under the only first-release boundary policy,
`escape`, is removed and
recorded as a lateral loss. Reflecting boundaries are forbidden because they
artificially confine particles. Any nonzero perpendicular diffusion or future
cross-field drift makes the corridor width a physical truncation. Schema 5
permits that case only in the qualified analytic-verification branch and
requires a width/buffer convergence sequence whose observer confidence
intervals and lateral-loss fraction meet preregistered bounds. A future
production capability must pass the same convergence gate. Until the
separatrix-crossing operator is qualified, any reachable open/closed interface
causes rejection.

### 12.6 Initial conditions

At `t=t0`:

1. PFSS and both topology maps are complete; when selected, finite-SCS
   coefficients, the constrained unsigned boundary, sector map, HCS, and
   resolved PFSS/SCS transition are also complete. Signed-potential
   reconstruction, composite-footpoint labeling, D9, and the selected
   transition-topology/exclusion mask finish before wind, source, observer,
   export, or particle setup, while the no-SCS
   verification branch records those products as typed inapplicable;
2. composite inner tubes, target-speed data, the selected resolution-invariant
   mass normalization, optional plasma-sheet modification, and an initial
   transonic-wind iterate are complete;
3. the wind/Parker-longitude-map iteration has converged to one fold-free,
   flux- and mass-conservative composite open state over the requested domain;
4. the closed hydrostatic provider and final open solution give every physical
   point exactly one thermodynamic authority;
5. open- and closed-field turbulence are complete under their explicitly
   selected policies;
6. cell-centered and vertex background values and categorical interface
   fields have been populated and halo
   exchanged;
7. the candidate-front surface and any explicitly enabled independent piston
   have been built;
8. shock admissibility and activation-history diagnostics are complete; zero
   initial active area is valid when the initial requirement is `none`;
9. every requested observer and, when export is enabled, field line has been
   traced/mapped, oriented, sampled, and checked against its source and
   energy-bin requirements; and
10. particle weights and time steps are initialized for every compiled AMPS
   species.

Initialization output and any field-line bundle are generated only after all
ten conditions pass. Bundle publication is transactional: an incomplete
line set, invalid shock history, or failed checksum leaves no apparently valid
final bundle.

### 12.7 Field-aligned `srcSEP` boundaries

Every imported line defines its own interval `[s_min,s_max]`. The inner end is
absorbing and must lie outside or on the application particle-absorption
surface. The outer end is absorbing/escape and must lie beyond the most distant
observer assigned to that line. An observer is represented by a stable
`LineObserverMapping` containing topology-segmented, time-resolved physical
overlap components, support intervals, closest coordinate, separation,
conservative represented volume/exposure, event times, and validity;
`srcSEP` must not infer it again from a nearest 3-D node.

The request registry, not the successful-export list, owns identity. Every
requested stable line/tube ID appears in the bundle manifest in original
order. A failed trace or unsupported interval remains as a typed rejection
record containing all reasons and the last valid support; later accepted lines
are not renumbered, and a rejected ID is never recycled for different
geometry.

If the observer does not lie on the selected line, the exporter reports the
minimum perpendicular separation. It may accept the association only when it
is below the explicitly configured connection tolerance. The actual mapped
interval(s), closest point, separation, tolerance, and coverage are retained
in the manifest. This tolerance
is a geometrical mapping uncertainty, not perpendicular particle diffusion.

---

## 13. Deriving model inputs from observations

Every observationally derived parameter must include source, instrument,
processing version, epoch, coordinate transformation, uncertainty, and file
checksum. A fitted value without those fields is incomplete.

### 13.1 Photospheric magnetic field

| Model input | Observation | Processing |
|---|---|---|
| `B_r(R_sun)` | GONG or SDO/HMI magnetogram | calibration, radial-field conversion, pole fill, flux balance, spherical-harmonic transform |
| solar orientation | ephemeris/SPICE or authoritative solar coordinates | transform to model frame at event epoch |
| `R_b`, `R_i`, `R_scs` | EUV coronal-hole boundaries, streamer/HCS morphology, open and unsigned flux | scan a joint ensemble; report interface kink and open-flux sensitivity; do not tune only to the SEP result |
| absolute open flux | rotation-matched in-situ radial magnetic field with uncertainty and folding correction | compare composite `Phi_open` on nested exterior spheres; optionally derive one positive magnetogram scale before PFSS/SCS construction |

The baseline analytic case instead declares harmonic coefficients directly and
does not claim event-specific topology.

For every event profile define

\[
 \Phi_{\rm open}^{\rm model}(r)=r^2\int_{4\pi}|B_r(r,\theta,\phi)|\,d\Omega
\]

on nested spheres outside `R_scs`. The corresponding observational product for
samples at positions `r(t)` is

\[
 \Phi_{\rm open}^{\rm obs}=4\pi
 \left\langle r(t)^2|B_r^{\rm in\ situ}(t)|\right\rangle,
\]

formed over the same Carrington rotation or registered time interval. Its
asset must declare spacecraft, cadence, gaps, averaging, uncertainty, sector
folding/disconnection correction, and coordinate conversion. The fixed-radius expression
`4*pi*r_obs^2*mean(|B_r|)` is legal only when the reference asset states that
every sample has already been normalized to the single declared comparison
radius. The comparison ratio and uncertainty form D6. Signed flux must close to zero, while unsigned
flux must remain invariant across the model's nested exterior spheres within
numerical tolerance. Coronal-hole agreement and open-flux magnitude are joint
constraints; matching either alone is insufficient [54]. Every D6 reference
also has exactly one data-use role: `construction`, `qualification`, or
`withheld-validation`. Byte-identical observations or overlapping samples may
not appear under two roles through differently named files.

`compare-only` leaves the lower boundary unchanged and requires D6 to pass.
`magnetogram-scale` uses an explicit two-pass background preparation. Pass A
builds the complete unscaled PFSS--transition--SCS--Parker background and raw
D6 product. The registered reference then derives one positive scalar. Pass B
discards the provisional state, applies that scalar **once** to the
flux-balanced photospheric coefficients before PFSS filtering, and rebuilds
and revalidates topology, SCS, composite field, wind, turbulence, and shock
inputs from the scaled authority. Only Pass B may be published. The unscaled
and scaled products are both retained, while no already-scaled intermediate is
ever fed back into calibration. The scale is

\[
 s_B=\frac{\Phi_{\rm open}^{\rm obs}}
          {\Phi_{\rm open}^{\rm model,unscaled}(r_{\rm cmp})},
\]

where `r_cmp` is deterministic: for `sample-position` it is the first element
of the canonically sorted `nested_sphere_radii_m`; for
`pre-normalized-to-comparison-radius` it is `comparison_radius_m`, which must
equal one member of that list. Nested-sphere closure must make the choice
immaterial within tolerance. The resulting `s_B` must lie inside its
preregistered admissible interval; otherwise the configuration fails rather
than clipping the factor. The scale
changes `B`, Alfvén speed, plasma beta, shock Mach numbers, and source
eligibility even when it leaves ideal field-line topology unchanged; all of
those consequences belong to the uncertainty campaign. Varying `R_b`, `R_i`,
or `R_scs` is not amplitude calibration because it can change topology. Such
radius variations remain a joint coronal-hole/streamer/open-flux ensemble and
are never selected by agreement with the SEP output.

For `magnetogram-scale`, the D6 reference used in the scale equation is a
**construction** asset. Agreement of Pass B with that same reference is a
closure/audit check and contributes no second likelihood factor or
qualification score. D6 may independently qualify or weight the tuple only
with a preregistered held-out product whose samples and processing lineage are
disjoint from the construction asset. In `compare-only`, a D6 asset may have a
qualification role because it did not set `s_B`; a withheld product remains
untouched until final validation.

The Alfvén-speed response must be recomputed from the coupled wind closure,
not described by an unconditional factor-of-two rule. From

\[
 v_A=\frac{|\mathbf B|}{\sqrt{\mu_0\rho}},
\]

a pure field rescaling gives `v_A -> s_B v_A` only when `rho` is held fixed.
If the structured wind instead holds mass loading per magnetic flux `eta_m`
and speed fixed, continuity gives `rho=eta_m |B|/u_s`, so
`rho -> s_B rho` and `v_A -> sqrt(s_B) v_A`. If speed, temperature, or mass
loading is re-solved, neither shortcut is authoritative; Pass B recomputes the
complete state and D1/D2 from first principles. Changing `R_b`, `R_i`, or
`R_scs` can alter both connectivity **and** low-coronal field magnitude through
the boundary-value solution. It therefore cannot be characterized as a
topology-only route to more open flux or substituted for amplitude scaling.

When shock-formation observations are used, open-flux calibration is therefore
not selected one axis at a time. The checksummed formation-height asset owns
the complete Cartesian product of magnetogram realization, scale mode/value,
radius tuple, density/wind member, and front member. D6 and coronal-hole/HCS
topology remain magnetic constraints, while D1/D2 and the independent event-
specific type-II/EUV likelihood constrain the resulting characteristic speeds
and activation height. Their joint preregistered rule may reject or weight a
tuple; it may not alter the scale factor inside a tuple, reuse a construction
D6 asset as qualification evidence, or use any SEP result.

### 13.2 Density and temperature

Candidate constraints include:

- EUV differential-emission-measure analysis in the low corona;
- polarized-brightness inversion from coronagraph data;
- spectroscopic density-sensitive line ratios;
- Doppler-dimming constraints on outflow speed;
- radio plasma frequency for a type-II source; and
- in-situ mass flux after propagation along the selected open tube.

The event-nominal kinematic product combines an inner electron-density or
mass-density reconstruction with an outer velocity reconstruction and species
temperatures, all registered against the topology and epoch of the magnetic
background. The inner/outer overlap `[r_a,r_b]` must be selected because both
products have support there; it is not an undocumented universal `3--5
R_sun` default. A candidate overlap is scanned with the products' native
resolution and uncertainty, then frozen before the SEP result is examined.
No accepted member may cross an unresolved data transition, topology event,
magnetic interface, field-line-ID change, or unsupported interval.

The normalization that fixes `eta_m` is observationally independent of a
lone speed. Acceptable authorities are a direct mass-per-flux estimate,
co-located density--velocity--field measurements/reconstructions, or radial
mass flux accompanied by a co-located mapped `B_r`. Each product states
whether it was used for construction, qualification, or withheld validation;
the same measurements cannot fill two roles. Electron-density products carry
the ion abundance and charge-state model used to obtain mass density, including
alpha/proton uncertainty and covariance. Velocity assets declare radial or
field-aligned component, inertial or corotating reference frame, a deterministic
consumer selector, and whether their
abscissa is heliocentric radius or oriented field-line arc length. A
field-aligned or arc-length asset additionally carries stable field-line IDs
and the exact background fingerprint; it is not a relabelled radial profile.

D7 reports fast/slow and tube-by-tube density, speed, temperature, and overlap
compatibility at 1.1, 2, and 5 solar radii when those radii are genuinely
covered, plus mass flux at
`solar_wind.mass_flux_reference_radius_m` for every structured branch. A
second target-radius mass flux is reported only for a finite-radius
target-speed polytrope; an asymptotic target and an empirical profile do not
invent such a radius. The comparisons use polarized-brightness/density
reconstructions, Doppler-dimming or tracked-flow speeds, and in-situ mass flux
with their full uncertainties [53, 55, 56]. The D7 coverage census is weighted
by represented open magnetic flux and open area, potential source
incident-number/energy flux, observer-footprint exposure, and requested export
support--never by the raw number of numerical tubes. A profile with no
observational coverage is a sensitivity input, not an event-grade closure.
Event qualification requires nonoverlapping data-use roles, complete required
consumer support, covariance or a named ensemble, and preregistered density,
speed, mass-flux, overlap-mismatch, momentum-residual, and coverage metrics.
An uncovered required source region, finite observer footprint, or exported
line interval fails the event profile rather than being filled by
extrapolation.

Closed-corona density and temperature constraints are kept separate from
open-wind mass flux. Streamer-belt density, temperature, speed, and angular
width constrain the optional empirical plasma-sheet closure. The Leblanc
corona-to-1-AU density profile [23] is an independent open-field diagnostic,
not a replacement for the composition-consistent tube solution.

The electron plasma frequency relation is

\[
 f_{pe}=\frac{1}{2\pi}\sqrt{\frac{n_e e^2}{m_e\epsilon_0}}
 \simeq 8.98\times10^3\sqrt{n_e[\mathrm{cm}^{-3}]}\ \mathrm{Hz}.
\]

Fundamental-versus-harmonic identification and the density model dominate the
uncertainty of a radio-derived height. Type-II band splitting must not be used
as a compression measurement without documenting its interpretation.

The preferred type-II constraint therefore preserves the calibrated dynamic-
spectrum lane in **frequency--time space**, including channel/time covariance,
lane-picking uncertainty, gaps, clock uncertainty, and the preregistered
fundamental/harmonic mixture. For every density/wind candidate, the provider
evaluates `n_e` on the candidate front, predicts `f_pe` and `2*f_pe`, and forms
that member's likelihood jointly with the EUV/front geometry and projection
covariance. It does not first turn the radio lane into one model-independent
height curve. If only a published height--time product is available, its
provenance must identify the density model used in the inversion: the product
may qualify only the matching density member unless the original radio data
and full cross-model covariance are supplied for a joint re-evaluation. The
same density-derived height product is never applied as independent evidence
to every alternative density member.

### 13.3 CME ejecta and shock geometry

Use EUV images for the earliest lateral expansion and multi-viewpoint
coronagraph images when available. The ejecta may be fitted with the Graduated
Cylindrical Shell model [20], while the faint outer front is fitted separately
as an ellipsoid [9]. The fit should provide:

- propagation longitude and latitude;
- tilt;
- center location;
- three semiaxes;
- front/piston separation; and
- covariance or an ensemble of acceptable fits.

Height-time and axis-time samples should be differentiated only after fitting
a smooth trajectory. Image cadence and feature-selection uncertainty must be
propagated into speed and acceleration. A two-point finite-difference speed is
not adequate for deciding whether a front at `1.01 R_sun` contains a locally
admissible fast-shock patch.

A PyThea/Kwon-compatible preprocessing adapter may accept center-plus-principal
axes or apex-plus-aspect-ratio observational coordinates, but it converts once
to the internal radial-principal-axis convention of Section 8.2. The original
parameters, frame, convention, covariance/ensemble, and conversion checksum are
retained. Center translation, each semiaxis, orientation, apex, fixed-angle
normal speeds, and the solar footprint are separate fitted/validated products;
an apex height-time curve alone does not determine flank speeds.
The PyThea release or commit, fit options, source-image identities, and exported
fit checksum are mandatory provenance [50].

### 13.4 Shock compression and magnetic field

White-light excess brightness can constrain the density enhancement when the
line-of-sight geometry and upstream density are modeled. Combined white-light,
EUV/UV spectroscopy, and Rankine--Hugoniot inversion can constrain compression,
heating, and coronal field strength [11, 12]. These products are validation data
for the stand-alone shock solver; they should not all be reused as independent
inputs in the same validation comparison.

### 13.5 Turbulence and scattering

Possible constraints are nonthermal line widths, density-fluctuation spectra,
radio scintillation, PSP/Solar Orbiter measurements farther out, and SEP
transport fits. A mean free path inferred by fitting the same SEP profile used
for validation is a calibrated event parameter, not independent validation.
Calibration and validation intervals must be separated.

### 13.6 SEP observations

Observer comparisons require instrument response, background subtraction,
energy-bin integration, species selection, dead-time/saturation handling, and
time standards. Model differential intensity is

\[
 j(E)=p^2f(p)
\]

under the usual isotropic phase-space convention, but the exact conversion
must match the distribution definition used by the mover and sampler.

---

## 18. Known limitations and required language for results

1. **The candidate front is prescribed, not generated.** The model determines
   patch by patch whether prescribed kinematics constitute an admissible fast
   shock. It does not predict eruption onset, front acceleration, or
   shock--piston standoff from the MHD force balance. A front initialized at
   `1.01--1.05 R_sun` may validly contain no shock patches.
2. **PFSS and SCS are magnetostatic constructions.** PFSS contains no
   active-region current or free energy, while the finite-shell SCS changes
   topology through a prescribed potential continuation and discrete polarity
   restoration. Neither represents CME-driven opening or reconnection.
3. **Solenoidality does not guarantee transition-sheet topology.** The
   vector-potential overlap is divergence free but does not, by that fact
   alone, preserve the inherited SCS HCS as a magnetic flux surface. Inside
   the overlap it is `S_tr`, a transition sector surface. The first-release
   production model excludes a convergence-qualified clearance around it from
   sources, particle transport, observer-feeding paths, and field-line export;
   excluded flux/area/rate and runtime losses are reported against
   preregistered physical-measure caps. Mask convergence does not make a large
   excluded source or observer fraction acceptable. A cap exceedance marks
   the event `not-event-grade`; repeated exceedance is evidence that the
   common-flux-surface or finite-thickness branch must be advanced. Results do not
   describe SEP transport through that neighborhood. A common-flux-surface or
   finite-thickness branch may replace the exclusion only after its separate
   topology and event-operator gates pass. Outside the overlap, the ideal HCS
   remains a zero-thickness sector discontinuity; HCS drift requires Stage 11A.
4. **The open wind is an effective closure.** A polytropic or WSA-targeted
   transonic tube does not calculate coronal heating, heat conduction,
   radiation, wave pressure, stream interaction, or CIR formation from first
   principles. The empirical-kinematic event branch combines an
   observation-constrained inner density with an outer velocity using one
   mass-per-flux authority; a speed alone cannot determine that authority.
   Its `C2` overlap and measured mismatch establish compatibility, not a
   solution of the heating problem. The branch enforces continuity but not the
   momentum equation; its measured residual bounds, rather than removes, this
   limitation. A radial velocity projection is valid only on outward segments
   that pass the projection guard, while a field-aligned/arc-length product is
   valid only on the exact fingerprinted lines for which it was produced.
   Missing required source, observer, or export coverage ends event-nominal
   validity. A longitude-map fold likewise marks the end of either branch's
   validity. The present pointwise and D7 gates do not constitute a continuous
   observational plausibility envelope between comparison radii; an
   event-specific, frame-explicit envelope is a future qualification, not an
   implied universal speed or monotonicity limit.
5. **The closed corona is hydrostatic and static.** It omits loop heating,
   evaporation, siphon flows, reconnection, and dynamic opening. The sharp
   separatrix or prescribed finite-width layer is a controlled baseline
   idealization. The default analytic composite is explicitly
   `diagnostic-kinematic`; it is not a stationary tangential discontinuity.
   Passing the hard normal-field/relative-flow gates and a bounded traction or
   volume-residual budget makes an event state a quantified approximation,
   not a global 3-D MHD equilibrium. Only a solver-produced or imported state
   that passes the full stationary-TD contract may use that equilibrium label.
6. **The heliospheric plasma sheet is empirical.** Its base density enhancement
   and width do not impose transverse MHD force balance or predict radial
   sheet broadening.
7. **The single-fluid shock EOS is limited.** The first release uses isotropic
   `gamma_ad=5/3`; it omits anisotropic pressure, heat flux, pickup ions,
   electron/ion nonequilibration, and kinetic shock structure. The critical-
   Mach classification inherits the domain of its versioned table. Budgeted
   exclusion of out-of-domain patches is a declared loss of source coverage,
   not an extrapolated criticality estimate.
8. **Moving-source mode has no global downstream CME.** Its preferred event
   branch describes an empirical net first-passage population at a finite
   physical upstream reference surface, not gross shock emission, a resolved
   sheath, magnetic cloud, downstream turbulence, or an accelerator-derived
   free-escape spectrum. The release fraction, conditional joint
   momentum--pitch distribution, and reference distance must be calibrated or
   propagated as distinct uncertainties. The region between the front and
   reference surface is unresolved source physics. A later post-reference
   return to the front is an absorbing transport loss until Stage 11B exists
   and must pass represented-number and birth-energy-normalized caps in
   immutable source cohorts. The conormal branch has no physical front loss
   and must not report a fabricated zero fraction.
   Schema 5 has no automatic physical inference for `L_ref`; it is a fixed
   calibration-owned length whose uncertainty must be propagated. Its
   diffusion-depth Péclet diagnostic is not a return probability. A
   momentum-dependent fixed-Péclet family and any return-to-release renewal
   kernel are future source models requiring new surface, conservation, and
   anti-double-counting contracts. Simple re-emission of returned particles is
   forbidden.
   The conormal no-through-flow boundary is a Parker sensitivity closure, not
   physical shock reflection. Shock-adjacent absorption is verification-only:
   its surviving Parker population tends to zero with placement distance and
   must never be used for event inference.
9. **Ellipsoidal kinematics remain an approximation.** Independent
   center/radial/lateral histories are more realistic than an apex-only law,
   but a fixed-orientation ellipsoid cannot represent arbitrary rotation,
   center-path deflection, deformation, or nonconvex fronts. Center deflection
   and body rotation are independent future histories and must not substitute
   for one another. Geometry inputs used in a fit are not validation outputs.
10. **Turbulence and mean free path are model dependent.** A prescribed wave
    energy does not uniquely determine pitch-angle scattering, and a fitted
    mean free path does not validate the wave model. Open and closed closures
    have different domains of applicability. Passing the wave-force gate only
    bounds the omitted field-aligned wave-pressure acceleration; it does not
    validate an empirical spectrum or every omitted turbulent heating term.
    Released particles do not amplify waves or modify their own scattering in
    schema 5. A prescribed post-reference scattering enhancement is at most a
    sensitivity proxy; a self-generated-wave claim requires a coupled
    spectral wave solver and particle--wave conservation.
11. **A radial source envelope is not shock weakening.** The hard/tapered
    termination acts only on injection. It cannot be interpreted as a
    prediction that the geometric front or MHD shock disappears at that
    radius.
12. **Disabled cells restrict transport.** Active-corridor convergence is
    mandatory, especially when perpendicular diffusion is enabled; no
    observer-feeding field line or source patch may intersect an inactive
    block.
13. **The three coupling radii are model parameters.** `R_b`, `R_i`, and
    `R_scs` have distinct mathematical roles. No `R_scs`, including
    `2.5 R_sun`, is production valid unless its realized harmonic attenuation
    and exterior latitude-flatness gates pass. Schema 5 starts radial Parker
    winding only at `R_scs`; a smaller independent winding radius requires a
    future generalized pushforward provider. Changing a radius may alter both
    topology and low-coronal field magnitude, so it is not a topology-only
    substitute for magnetogram scaling. Formation-height validation uses the
    complete joint candidate product rather than a favored one-dimensional
    scan.
14. **The low-coronal source is not a 1-AU shock model.** The baseline source
    ends near `20 R_sun`; extending shock acceleration through the heliosphere
    requires a separately validated propagation, weakening, and sheath model.
15. **Field-line connection is model-relative.** A root can be geometric yet
    sub-fast, subcritical, outside the source envelope, or terminated. The
    documented status and history horizon must accompany every statement that
    an observer is “connected” or “disconnected.” A line rejected for touching
    the transition clearance is `topology-inapplicable`, not observational
    evidence of disconnection, and cannot be replaced by a nearby line.
16. **Independent 1-D lines have no perpendicular coupling.** Multiple
    `srcSEP` lines sample different magnetic connections but cannot reproduce
    cross-field diffusion or drift between them. Sector/HCS drift is absent
    from the stand-alone baseline. Gradient/curvature drift is also absent in
    both production movers; adding it requires a frame- and energy-consistent
    operator with guiding-center and interface validity guards.
17. **The validated transport range ends at 1 AU.** A longer exported line is
    an API capability, not evidence that the wind, turbulence, and scattering
    closures are valid to 5 AU.
18. **The baseline has no flare-associated particle source.** It can test a
    shock-only hypothesis but cannot uniquely attribute early particles to or
    exclude a flare contribution. A future flare-associated source is an
    independently ledgered attribution sensitivity whose observational
    proxies do not directly determine an escaping-ion rate or spectrum.

A publication should call this a *stand-alone analytic composite
PFSS--SCS--Parker background with an observation-constrained kinematic
candidate front and locally admissible MHD-shock source*. It should not call
the result a self-consistent CME MHD simulation.

---

## 19. Supporting references

1. Parker, E. N. (1958), “Dynamics of the Interplanetary Gas and Magnetic
   Fields,” *Astrophysical Journal*, 128, 664--676,
   <https://doi.org/10.1086/146579>.
2. Altschuler, M. D., and G. Newkirk (1969), “Magnetic Fields and the Structure
   of the Solar Corona. I: Methods of Calculating Coronal Fields,” *Solar
   Physics*, 9, 131--149, <https://doi.org/10.1007/BF00145734>.
3. Schatten, K. H., J. M. Wilcox, and N. F. Ness (1969), “A Model of
   Interplanetary and Coronal Magnetic Fields,” *Solar Physics*, 6, 442--455,
   <https://doi.org/10.1007/BF00146478>.
4. National Solar Observatory, “Potential-Field Source-Surface Models,”
   <https://nso.edu/data/nisp-data/pfss/>.
5. Kopp, R. A., and T. E. Holzer (1976), “Dynamics of Coronal Hole Regions. I.
   Steady Polytropic Flows with Multiple Critical Points,” *Solar Physics*, 49,
   43--56, <https://doi.org/10.1007/BF00221484>.
6. Cranmer, S. R., A. A. van Ballegooijen, and R. J. Edgar (2007),
   “Self-consistent Coronal Heating and Solar Wind Acceleration from
   Anisotropic Magnetohydrodynamic Turbulence,” *Astrophysical Journal
   Supplement Series*, 171, 520--551, <https://doi.org/10.1086/518001>.
7. Zhou, Y., and W. H. Matthaeus (1990), “Transport and Turbulence Modeling of
   Solar Wind Fluctuations,” *Journal of Geophysical Research*, 95,
   10291--10311, <https://doi.org/10.1029/JA095iA07p10291>.
8. Jokipii, J. R. (1966), “Cosmic-Ray Propagation. I. Charged Particles in a
   Random Magnetic Field,” *Astrophysical Journal*, 146, 480,
   <https://doi.org/10.1086/148912>.
9. Kwon, R.-Y., J. Zhang, and O. Olmedo (2014), “New Insights into the Physical
   Nature of Coronal Mass Ejections and Associated Shock Waves within the
   Framework of the Three-dimensional Structure,” *Astrophysical Journal*,
   794, 148, <https://doi.org/10.1088/0004-637X/794/2/148>.
10. Kwon, R.-Y., and A. Vourlidas (2018), “The Density Compression Ratio of
    Shock Fronts Associated with Coronal Mass Ejections,” *Journal of Space
    Weather and Space Climate*, 8, A08,
    <https://doi.org/10.1051/swsc/2017045>.
11. Bemporad, A., and S. Mancuso (2010), “First Complete Determination of
    Plasma Physical Parameters across a Coronal Mass Ejection-driven Shock,”
    *Astrophysical Journal*, 720, 130--143,
    <https://doi.org/10.1088/0004-637X/720/1/130>.
12. Bemporad, A. (2016), “Measuring Coronal Magnetic Fields with Remote
    Sensing Observations of Shock Waves,” *Frontiers in Astronomy and Space
    Sciences*, 3, 17, <https://doi.org/10.3389/fspas.2016.00017>.
13. Edmiston, J. P., and C. F. Kennel (1984), “A Parametric Survey of the First
    Critical Mach Number for a Fast MHD Shock,” *Journal of Plasma Physics*,
    32, 429--441, <https://doi.org/10.1017/S002237780000218X>.
14. Gopalswamy, N., et al. (2013), “Height of Shock Formation in the Solar
    Corona Inferred from Observations of Type II Radio Bursts and Coronal Mass
    Ejections,” *Advances in Space Research*, 51, 1981--1989,
    <https://doi.org/10.1016/j.asr.2013.01.006>.
15. Blandford, R. D., and J. P. Ostriker (1978), “Particle Acceleration by
    Astrophysical Shocks,” *Astrophysical Journal Letters*, 221, L29--L32,
    <https://doi.org/10.1086/182658>.
16. Drury, L. O'C. (1983), “An Introduction to the Theory of Diffusive Shock
    Acceleration of Energetic Particles in Tenuous Plasmas,” *Reports on
    Progress in Physics*, 46, 973--1027,
    <https://doi.org/10.1088/0034-4885/46/8/002>.
17. Parker, E. N. (1965), “The Passage of Energetic Charged Particles through
    Interplanetary Space,” *Planetary and Space Science*, 13, 9--49,
    <https://doi.org/10.1016/0032-0633(65)90131-5>.
18. Roelof, E. C. (1969), “Propagation of Solar Cosmic Rays in the
    Interplanetary Magnetic Field,” in *Lectures in High-Energy Astrophysics*,
    NASA SP-199, NASA NTRS document 19690020281,
    <https://ntrs.nasa.gov/citations/19690020281>.
19. Skilling, J. (1971), “Cosmic Rays in the Galaxy: Convection or Diffusion?”
    *Astrophysical Journal*, 170, 265,
    <https://doi.org/10.1086/151210>.
20. Thernisien, A. (2011), “Implementation of the Graduated Cylindrical Shell
    Model for the Three-dimensional Reconstruction of Coronal Mass Ejections,”
    *Astrophysical Journal Supplement Series*, 194, 33,
    <https://doi.org/10.1088/0067-0049/194/2/33>.
21. Gopalswamy, N., et al. (2013), “The First Ground Level Enhancement Event of
    Solar Cycle 24: Direct Observation of Shock Formation and Particle Release
    Heights,” *Astrophysical Journal Letters*, 765, L30,
    <https://doi.org/10.1088/2041-8205/765/2/L30>.
22. Lario, D., et al. (2014), “The Solar Energetic Particle Event on 2013 April
    11: An Investigation of Its Solar Origin and Longitudinal Spread,”
    *Astrophysical Journal*, 797, 8,
    <https://doi.org/10.1088/0004-637X/797/1/8>.
23. Leblanc, Y., G. A. Dulk, and J.-L. Bougeret (1998), “Tracing the Electron
    Density from the Corona to 1 AU,” *Solar Physics*, 183, 165--180,
    <https://doi.org/10.1023/A:1005049730506>.
24. Lee, M. A. (1983), “Coupled Hydromagnetic Wave Excitation and Ion
    Acceleration at Interplanetary Traveling Shocks,” *Journal of Geophysical
    Research*, 88, 6109--6120,
    <https://doi.org/10.1029/JA088iA08p06109>.
25. Vainio, R., and T. Laitinen (2007), “Monte Carlo Simulations of Coronal
    Diffusive Shock Acceleration in Self-generated Turbulence,”
    *Astrophysical Journal*, 658, 622--630,
    <https://doi.org/10.1086/510284>.
26. Vršnak, B., et al. (2013), “Propagation of Interplanetary Coronal Mass
    Ejections: The Drag-Based Model,” *Solar Physics*, 285, 295--315,
    <https://doi.org/10.1007/s11207-012-0035-4>.
27. Žic, T., B. Vršnak, and M. Temmer (2015), “Heliospheric Propagation of
    Coronal Mass Ejections: Drag-Based Model Fitting,” *Astrophysical Journal
    Supplement Series*, 218, 32,
    <https://doi.org/10.1088/0067-0049/218/2/32>.
28. Rollett, T., et al. (2016), “ElEvoHI: A Novel CME Prediction Tool for
    Heliospheric Imaging Combining an Elliptical Front with Drag-Based Model
    Fitting,” *Astrophysical Journal*, 824, 131,
    <https://doi.org/10.3847/0004-637X/824/2/131>.
29. Shi, C., M. Velli, S. D. Bale, V. Réville, M. Maksimović, and J.-B. Dakeyo
    (2022), “Acceleration of Polytropic Solar Wind: Parker Solar Probe
    Observation and One-Dimensional Model,” *Physics of Plasmas*, 29(12),
    122901, <https://doi.org/10.1063/5.0124703>.
30. Mikić, Z., J. A. Linker, D. D. Schnack, R. Lionello, and A. Tarditi
    (1999), “Magnetohydrodynamic Modeling of the Global Solar Corona,”
    *Physics of Plasmas*, 6(5), 2217--2224,
    <https://doi.org/10.1063/1.873474>.
31. Rosner, R., W. H. Tucker, and G. S. Vaiana (1978), “Dynamics of the
    Quiescent Solar Corona,” *Astrophysical Journal*, 220, 643--665,
    <https://doi.org/10.1086/155949>.
32. Rouillard, A. P., et al. (2016), “Deriving the Properties of Coronal
    Pressure Fronts in 3D: Application to the 2012 May 17 Ground Level
    Enhancement,” *Astrophysical Journal*, 833(1), 45,
    <https://doi.org/10.3847/1538-4357/833/1/45>.
33. Kouloumvakos, A., A. P. Rouillard, Y. Wu, R. Vainio, A. Vourlidas,
    I. Plotnikov, A. Afanasiev, and H. Önel (2019), “Connecting the Properties
    of Coronal Shock Waves with Those of Solar Energetic Particles,”
    *Astrophysical Journal*, 876(1), 80,
    <https://doi.org/10.3847/1538-4357/ab15d7>.
34. Schatten, K. H. (1971), “Current Sheet Magnetic Model for the Solar
    Corona,” *Cosmic Electrodynamics*, 2, 232--245.
35. McGregor, S. L., W. J. Hughes, C. N. Arge, and M. J. Owens (2008),
    “Analysis of the Magnetic Field Discontinuity at the Potential Field Source
    Surface and Schatten Current Sheet Interface in the Wang--Sheeley--Arge
    Model,” *Journal of Geophysical Research: Space Physics*, 113(A8), A08112,
    <https://doi.org/10.1029/2007JA012330>.
36. Smith, E. J., and A. Balogh (1995), “Ulysses Observations of the Radial
    Magnetic Field,” *Geophysical Research Letters*, 22(23), 3317--3320,
    <https://doi.org/10.1029/95GL02826>.
37. Winterhalter, D., E. J. Smith, M. E. Burton, N. Murphy, and D. J. McComas
    (1994), “The Heliospheric Plasma Sheet,” *Journal of Geophysical Research:
    Space Physics*, 99(A4), 6667--6680,
    <https://doi.org/10.1029/93JA03481>.
38. Huang, J., et al. (2023), “Parker Solar Probe Observations of High Plasma
    Beta Solar Wind from the Streamer Belt,” *Astrophysical Journal Supplement
    Series*, 265(2), 47, <https://doi.org/10.3847/1538-4365/acbcd2>.
39. Arge, C. N., and V. J. Pizzo (2000), “Improvement in the Prediction of
    Solar Wind Conditions Using Near-Real-Time Solar Magnetic Field Updates,”
    *Journal of Geophysical Research: Space Physics*, 105(A5), 10465--10479,
    <https://doi.org/10.1029/1999JA000262>.
40. McGregor, S. L., W. J. Hughes, C. N. Arge, M. J. Owens, and D. Odstrcil
    (2011), “The Distribution of Solar Wind Speeds during Solar Minimum:
    Calibration for Numerical Solar Wind Modeling Constraints on the Source of
    the Slow Solar Wind,” *Journal of Geophysical Research: Space Physics*,
    116(A3), A03101, <https://doi.org/10.1029/2010JA015881>.
41. Gary, G. A. (2001), “Plasma Beta above a Solar Active Region: Rethinking
    the Paradigm,” *Solar Physics*, 203(1), 71--86,
    <https://doi.org/10.1023/A:1012722021820>.
42. Tóth, G., B. van der Holst, and Z. Huang (2011), “Obtaining Potential Field
    Solutions with Spherical Harmonics and Finite Differences,”
    *Astrophysical Journal*, 732(2), 102,
    <https://doi.org/10.1088/0004-637X/732/2/102>.
43. Kennel, C. F., and J. P. Edmiston (1988), “Switch-On Shocks,” *Journal of
    Geophysical Research: Space Physics*, 93(A10), 11363--11373,
    <https://doi.org/10.1029/JA093iA10p11363>.
44. Feng, H. Q., C. C. Lin, J. K. Chao, D. J. Wu, L. H. Lyu, and L. C. Lee
    (2009), “Observations of an Interplanetary Switch-On Shock Driven by a
    Magnetic Cloud,” *Geophysical Research Letters*, 36(7), L07106,
    <https://doi.org/10.1029/2009GL037354>.
45. Lang, J. T., R. D. Strauss, N. E. Engelbrecht, J. P. van den Berg,
    N. Dresing, D. Ruffolo, and R. Bandyopadhyay (2024), “A Detailed Survey of
    the Parallel Mean Free Path of Solar Energetic Particle Protons and
    Electrons,” *Astrophysical Journal*, 971(1), 105,
    <https://doi.org/10.3847/1538-4357/ad55c3>.
46. Subashchandar, N. S. M., et al. (2025), “Parallel and Perpendicular
    Diffusion of Energetic Particles in the Near-Sun Solar Wind Observed by
    Parker Solar Probe,” *Astrophysical Journal Letters*, 991(2), L30,
    <https://doi.org/10.3847/2041-8213/ae063f>.
47. Battarbee, M., S. Dalla, and M. S. Marsh (2017), “Solar Energetic Particle
    Transport Near a Heliospheric Current Sheet,” *Astrophysical Journal*,
    836(1), 138, <https://doi.org/10.3847/1538-4357/836/1/138>.
48. Battarbee, M., S. Dalla, and M. S. Marsh (2018), “Modeling Solar Energetic
    Particle Transport near a Wavy Heliospheric Current Sheet,”
    *Astrophysical Journal*, 854(1), 23,
    <https://doi.org/10.3847/1538-4357/aaa3fa>.
49. Mann, G., A. Klassen, H. Aurass, and H.-T. Classen (2003), “Formation and
    Development of Shock Waves in the Solar Corona and the Near-Sun
    Interplanetary Space,” *Astronomy & Astrophysics*, 400(1), 329--336,
    <https://doi.org/10.1051/0004-6361:20021593>.
50. Kouloumvakos, A., L. Rodríguez-García, J. Gieseler, D. Price,
    A. Vourlidas, and R. Vainio (2022), “PyThea: An Open-Source Software
    Package to Perform 3D Reconstruction of Coronal Mass Ejections and Shock
    Waves,” *Frontiers in Astronomy and Space Sciences*, 9, 974137,
    <https://doi.org/10.3389/fspas.2022.974137>.
51. Schatten, K. H. (1971), “Current Sheet Magnetic Model for the Solar
    Corona,” NASA-TM-X-65496 / X-692-71-132, NASA NTRS document 19710013957,
    <https://ntrs.nasa.gov/citations/19710013957>.
52. Caprioli, D., E. Amato, and P. Blasi (2010), “Non-linear Diffusive Shock
    Acceleration with Free Escape Boundary,” *Astroparticle Physics*, 33,
    307--311, <https://doi.org/10.1016/j.astropartphys.2010.03.001>.
53. Leer, E., and T. E. Holzer (1980), “Energy Addition in the Solar Wind,”
    *Journal of Geophysical Research*, 85(A9), 4681--4688,
    <https://doi.org/10.1029/JA085iA09p04681>.
54. Linker, J. A., et al. (2017), “The Open Flux Problem,” *Astrophysical
    Journal*, 848, 70, <https://doi.org/10.3847/1538-4357/aa8a70>.
55. Sheeley, N. R., Jr., et al. (1997), “Measurements of Flow Speeds in the
    Corona Between 2 and 30 Solar Radii,” *Astrophysical Journal*, 484,
    472--478, <https://doi.org/10.1086/304338>.
56. Guhathakurta, M., T. E. Holzer, and R. M. MacQueen (1996), “The
    Large-Scale Density Structure of the Solar Corona and the Heliospheric
    Current Sheet,” *Astrophysical Journal*, 458, 817--831,
    <https://doi.org/10.1086/176860>.
57. Dalla, S., M. S. Marsh, J. Kelly, and T. Laitinen (2013), “Solar Energetic
    Particle Drifts in the Parker Spiral,” *Journal of Geophysical Research:
    Space Physics*, 118(10), 5979--5985,
    <https://doi.org/10.1002/jgra.50589>.
58. Marsh, M. S., S. Dalla, J. Kelly, and T. Laitinen (2013), “Drift-Induced
    Perpendicular Transport of Solar Energetic Particles,” *Astrophysical
    Journal*, 774(1), 4, <https://doi.org/10.1088/0004-637X/774/1/4>.
59. Reames, D. V., and C. K. Ng (1998), “Streaming-Limited Intensities of Solar
    Energetic Particles,” *Astrophysical Journal*, 504, 1002--1005,
    <https://doi.org/10.1086/306124>.
60. Ng, C. K., D. V. Reames, and A. J. Tylka (2003), “Modeling
    Shock-Accelerated Solar Energetic Particles Coupled to Interplanetary
    Alfvén Waves,” *Astrophysical Journal*, 591(1), 461--485,
    <https://doi.org/10.1086/375293>.
61. Cheng, L., M. Zhang, D. Lario, L. A. Balmaceda, R.-Y. Kwon, and
    C. M. S. Cohen (2023), “Simulation of the Solar Energetic Particle Event
    on 2020 May 29 Observed by Parker Solar Probe,” *Astrophysical Journal*,
    943(2), 134, <https://doi.org/10.3847/1538-4357/acac21>.
62. Zhuang, B., N. Lugaz, and D. Lario (2022), “Widespread 1--2 MeV Energetic
    Particles Associated with Slow and Narrow Coronal Mass Ejections: Parker
    Solar Probe and STEREO Measurements,” *Astrophysical Journal*, 925(1), 96,
    <https://doi.org/10.3847/1538-4357/ac3af2>.
63. Palmerio, E., E. K. J. Kilpua, O. Witasse, D. Barnes, B. Sánchez-Cano,
    A. J. Weiss, et al. (2021), “CME Magnetic Structure and IMF Preconditioning
    Affecting SEP Transport,” *Space Weather*, 19(4), e2020SW002654,
    <https://doi.org/10.1029/2020SW002654>.
64. Kay, C., M. Opher, and R. M. Evans (2013), “Forecasting a Coronal Mass
    Ejection's Altered Trajectory: ForeCAT,” *Astrophysical Journal*, 775(1),
    5, <https://doi.org/10.1088/0004-637X/775/1/5>.
65. Kay, C., M. Opher, R. C. Colaninno, and A. Vourlidas (2016), “Using
    ForeCAT Deflections and Rotations to Constrain the Early Evolution of
    CMEs,” *Astrophysical Journal*, 827(1), 70,
    <https://doi.org/10.3847/0004-637X/827/1/70>.

---
