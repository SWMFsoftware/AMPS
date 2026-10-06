# Reduced corona–SWCME shock-front model

**Version:** 1.1 — 2026-10-04  
**Purpose:** a detailed physical and implementation specification for a prescribed, finite, three-dimensional shock front moving through a coronal and heliospheric ambient plasma and magnetic field.  
**Initial application:** `srcSEP3D`, with zero SEP particles. Shared algorithms belong under `src/models`; `srcSEP` integration and particle modeling follow separately.  
**Suggested repository location:** `src/models/sep_corona_swcme/docs/shock-front/model.md`.

This is a specification for the reduced mode. It is not a claim that the current source implements or passes every requirement below. Proposed configuration keys, API names, test identifiers, and example paths must be reconciled with the actual repository before implementation. Existing full-CME design documents and the unfinished sheath/ejecta implementation are retained separately.

**Revision:** incorporates the qualified corrections from the shock-front model review. In particular, geometric support is distinguished from accepted shock extent, body-drag proxy validity is tested, ambient sensitivity and flux-consistent calibration are required, and useful surface diagnostics are specified. The analytical fixture parameters remain unchanged. Driver/standoff dynamics, general front dynamics, and nonlinear formation remain separately qualified extensions.

## 1. Scientific objective and meaning of the model

The model answers the following questions throughout a prescribed evolution from the low corona to at least 1 AU:

1. Where is the finite disturbance front, and how fast does it move along its local normal?
2. What ambient magnetic field and plasma does each supported part of the front encounter?
3. Which parts are compatible with an outward fast shock under the selected jump model?
4. What are the local fast-mode Mach number, magnetic obliquity, density compression, and immediate downstream jump state?
5. Which ambient magnetic field lines intersect the supported shock, and when do those connections appear or disappear?
6. Does the front pass a specified observer, and what are the local front and jump parameters there?

The front trajectory is supplied by a smooth analytical history, a resolved observational history, or a documented outer propagation law. The ambient provider supplies plasma and magnetic data independently. Their combination defines local shock diagnostics. The trajectory is not obtained by solving the global momentum and energy equations of an erupting CME.

The word **launch** in this mode means initializing and accelerating a prescribed front near the Sun. It does not mean generating an instability, calculating an erupting flux rope, or proving that an initially smooth wave has steepened into a shock. A prescribed front can be sub-fast on some or all of its surface; those points remain geometric front points and have no accepted fast-shock jump.

This approach is useful for shock geometry, magnetic connectivity, and subsequent explicitly reduced SEP source/transport studies. It does not replace a complete spatial CME plasma solution for applications that require the sheath or ejecta interior. Observational front reconstruction combined with coronal field modeling provides a relevant precedent [R1].

### 1.1 The three represented objects

| Object | Represented quantities | Physical meaning |
| --- | --- | --- |
| Ambient reference | `rho`, `p`, `T`, `U`, `B`, EOS metadata, spatial/time validity | The prescribed corona/solar wind in the absence of the modeled disturbance |
| Supported front | Position, implicit geometry, normal, normal speed, finite extent, shape history | A prescribed propagating surface; not every point is necessarily a shock |
| Local shock record | Upstream state, admissibility, Mach numbers, obliquity, compression, immediate downstream state, jump residuals | One-sided jump limits at an accepted local fast shock |

The model has no sheath thickness, material contact, ejecta density, ejecta magnetic field, or downstream-volume velocity law. It therefore has no sheath inventory or contact-leakage budget to close. Such quantities must not be fabricated from a geometric contact fraction or local jump compression.

### 1.2 Explicit capability boundaries

An ambient reference can be evaluated at a mathematically supported point even if that point lies geometrically behind the prescribed front. In that case it remains an **undisturbed reference value**. It is not the actual post-CME plasma at that location.

The immediate downstream state is attached to the front record. It cannot be extended into an arbitrary downstream volume by copying it along rays, relaxing its velocity, or multiplying ambient fields by a compression ratio.

A consumer requesting a physical downstream volume receives a typed unsupported-capability result. A plotting consumer can request the ambient reference over the whole mesh, provided the output records this role and masks the interior as having no modeled CME plasma state.

### 1.3 Reading and qualification guide

| Purpose | Sections |
| --- | --- |
| Physical scope, frames, and shared authority | 1–4 |
| Ambient plasma, magnetic reconstruction, and calibration | 5–7 |
| Front geometry, trajectory validity, and handoff | 8–10 |
| Local shocks, diagnostic limitations, and connectivity | 11–13 |
| Proposed input/output and native execution contracts | 14–17 |
| Reproducible fixtures and required evidence | 18–20 |
| Implementation order, optional extensions, and references | 21–24 |

A conditional calculation answers “what shock properties would this prescribed front have in this ambient?” It becomes an event prediction only after the trajectory, plasma, field, and uncertainty have been independently assessed. Reaching 1 AU geometrically is a different gate from retaining an accepted shock there.

## 2. Assumptions and approximations

### 2.1 Core assumptions

| Assumption | Consequence and validation obligation |
| --- | --- |
| Nonrelativistic bulk plasma | All front and wind speeds must be well below the speed of light |
| One-fluid, scalar-pressure jump model | A consistent total plasma pressure and declared adiabatic index are required; pressure anisotropy and separate electron/ion shock heating are not solved |
| Locally thin shock | Jump conditions relate one-sided limits; the microscopic shock ramp is unresolved |
| Local surface planarity | Shock thickness and the upstream sampling stencil must be small relative to curvature and ambient variation lengths |
| Prescribed front kinematics | Front shape, acceleration, and deflection are inputs or reduced propagation outputs; global CME dynamics are not inferred |
| Prescribed ambient corona/wind | The eruption does not alter the ambient provider in this mode |
| Ideal magnetic jump conditions | Normal magnetic flux and tangential electric field are continuous across an accepted ideal-MHD shock |
| Test-background regime | No SEP pressure, particle current, wave growth, or particle feedback modifies the front or ambient |
| Finite supported surface | Unsupported flanks, artificial cap closures, and solar-interior points do not become shocks |
| Explicit input validity | Histories and magnetic maps are used only over their declared intervals and domains |

Using ideal-MHD **local jump relations** does not introduce a global MHD solver. It supplies the conservation constraints needed to calculate a physically admissible local compression and downstream limit.

### 2.2 Effects deliberately deferred

- Self-consistent CME eruption, reconnection, and flux-rope dynamics.
- Pressure-driven shock/contact standoff and downstream momentum evolution.
- Sheath mass accumulation, lateral flow, material maps, and magnetic draping.
- Ejecta expansion thermodynamics and internal magnetic topology.
- Eruption-induced opening or reconnection of ambient coronal field lines.
- CME–CME interaction and magnetic preconditioning, unless explicitly supplied as part of a separate ambient model.
- Collisionless ramp structure, cross-shock potential, electron/ion heating partition, anisotropic pressure, and nonthermal jump corrections.
- SEP injection, scattering, acceleration, transport, and wave feedback.

These omissions limit which physical conclusions can be drawn. In particular, an ambient connectivity calculation cannot establish the connectivity through the CME-disturbed downstream region.

### 2.3 Three levels of evidence

1. **Numerical verification:** geometry derivatives, field evaluation, local jump equations, event location, and host/MPI storage reproduce independent references.
2. **Reduced-model physical qualification:** the selected assumptions are appropriate for the stated use, inputs are compatible, and sensitivity and observational discrepancies are quantified.
3. **Event-specific validation:** the actual magnetic map, front history, plasma parameters, and observer geometry reproduce independent measurements within declared uncertainty.

A synthetic dipole fixture can pass level 1 without establishing levels 2 or 3. Passing this reduced model does not complete the deferred BG3D-4 sheath stage.

## 3. Units, coordinates, frames, and signs

### 3.1 Canonical units and constants

All resolved internal quantities use SI. Human-readable inputs can use explicit units if the parser supports them, but the event manifest stores the resolved SI values.

| Quantity | Internal unit |
| --- | --- |
| Position, radius, semiaxis | m |
| Time from reference epoch | s |
| Velocity and front normal speed | m s^-1 |
| Acceleration | m s^-2 |
| Mass density | kg m^-3 |
| Number density | m^-3 |
| Pressure | Pa |
| Temperature | K |
| Magnetic field | T |
| Field gradient | T m^-1 |
| Angles | rad |
| Drag coefficient `Gamma` | m^-1 |

Reference constants for examples are

\[
R_\odot=6.957\times10^8\ {\rm m},\qquad
1\ {\rm AU}=149\,597\,870\,700\ {\rm m}.
\]

Use the project's canonical constants and record their values in each manifest. Examples must not create a second inconsistent constants table. The permeability `mu0`, Boltzmann constant, species masses, and solar gravitational parameter are likewise recorded.

### 3.2 Spatial and time conventions

The canonical spatial frame is heliocentric inertial, **HCI**, with the Sun at the origin. The precise HCI realization, reference equinox/axis convention, and transformation version must be named. Photospheric maps supplied in Carrington coordinates require an explicit epoch-dependent rotation into HCI; Carrington longitude is not interchangeable with HCI longitude.

For spherical coordinates, `r` is heliocentric radius, `theta` is colatitude measured from the positive solar-rotation axis, and `phi` increases in the right-handed azimuthal direction. Cartesian vectors are preferred at poles and during shock solves.

Store the absolute reference epoch and its time scale. If the input timestamp is UTC, convert elapsed times consistently using the declared clock implementation. A sequence of nominal UTC timestamps must not silently ignore a leap-second crossing. The integration clock itself uses a monotone elapsed-SI-second coordinate.

All published plasma velocities and front speeds are in the HCI/Sun frame unless explicitly marked otherwise. A velocity measured in a rotating frame must be transformed, including the rotational velocity term, before use in a jump calculation.

### 3.3 Surface and shock signs

Let the level-set function `g(x,t)` be negative on the geometrical inner side and positive on the outer side of a supported front. Then

\[
\mathbf n=\frac{\nabla g}{|\nabla g|}
\]

points from the immediate downstream side toward the upstream side, generally away from the CME body. The local front normal speed is positive for motion along `n`.

For an outward forward shock,

\[
w_1=V_n-\mathbf U_1\cdot\mathbf n>0.
\]

The shock-frame upstream normal velocity is `u1n=-w1`. This distinction prevents sign mistakes between positive inflow magnitude and signed velocity.

The normal to a tilted or curved front is generally different from the heliocentric radial direction. Apex radial speed, local radial intersection speed, surface parameterization speed, and local normal speed are distinct outputs.

## 4. Model architecture and authority

### 4.1 Shared components

The reduced model comprises:

1. **Resolved event:** immutable parameters, assets, coordinate transformations, uncertainty description, and checksums.
2. **Ambient provider:** canonical coronal and heliospheric reference queries.
3. **Front evolution provider:** one continuous geometry/history authority from initial epoch to the endpoint.
4. **Local jump evaluator:** the maintained shared SWCME shock solver, after interface and regression verification.
5. **Connectivity/event evaluator:** field-line/front intersections and observer passage.
6. **Host adapter:** `srcSEP3D` storage, MPI epoch publication, masks, and diagnostics.

Shared launch/background functionality stays under `src/models/sep_coronal_cme` and/or the existing composite package `src/models/sep_corona_swcme`, according to the repository's established ownership. The SWCME jump solver remains under `src/models/swcme`. Avoid duplicating shared equations inside `srcSEP3D` or creating another independent corona implementation.

### 4.2 One continuous front authority

Use one composite identity and one immutable front epoch throughout the run. A change from coronal history evaluation to an outer SWCME propagation kernel is a phase change inside this authority. It is not a switch to another native background authority or a fresh initialization of radius, direction, width, speed, or time.

The ambient magnetic authority is also distinct from the front's propagation law. Entering the SWCME propagation phase must not replace a spatial coronal field with an unrelated radial/Parker field merely because a different front kernel is now active.

### 4.3 Separation from particles and the full CME mode

No source, injection spectrum, particle distribution, turbulence obligation, sheath parcel map, or contact asset is required to resolve or prepare a reduced front epoch. Existing richer epochs can remain for later work, but this mode must not depend on their particle-only or sheath-only prerequisites.

Do not alter generic PIC behavior for unrelated applications. Use the same maintained runtime-provider/native-storage mechanism already used for SEP3D–SWCME coupling. Native `DATAFILE` storage describes a buffer interface; it does not mean that the physical provider is a scheduled file reader.

## 5. Ambient plasma: what is needed and what is prescribed

### 5.1 Canonical upstream state

At a front point, the upstream provider returns

\[
\mathcal P_1=(\rho_1,p_1,\mathbf U_1,\mathbf B_1,
\gamma,\text{composition},\text{frame},\text{generation}).
\]

Temperature is derived from, or checked against, the declared EOS. Magnetic field alone is insufficient: Mach number and compression also depend on density, pressure, flow, and the shock-relative speed.

For a fully ionized hydrogen example,

\[
\rho\simeq m_p n_p,\qquad
p=n_p k_B T_p+n_e k_B T_e,\qquad n_e=n_p.
\]

If `Tp=Te=T`, then `p=2 np kB T`. Using `p=np kB T` while describing `T` as a common electron/ion temperature underestimates total pressure. Helium or other compositions require the appropriate mass density and total particle pressure. EOS composition is independent of the compiled energetic-particle species list.

The default local shock EOS can use a single ideal gas with `gamma=5/3`, with scalar total pressure and internal-energy density `p/(gamma-1)`. A different gamma is an explicit physical input. Do not tune it just to produce a desired compression.

### 5.2 Ambient-model choices

| Choice | Appropriate use | Required declaration |
| --- | --- | --- |
| Uniform analytical state | Independent geometry and jump verification | Synthetic; not a realistic low-coronal event |
| Canonical analytical corona/wind | Reduced production studies over its calibrated domain | Plasma profile equations, field/flow compatibility, EOS, and limitations |
| Empirical density/temperature with reconstructed field | Event diagnostics | Data provenance, uncertainty, interpolation, and lack of full force balance if applicable |
| Field-aligned flux-tube wind | Improved open-field ambient | Tube area, critical-point solution, mass flux, and open/closed routing |

Reuse the canonical provider when it meets these contracts. This specification does not authorize replacing it with a new wind solver unnecessarily.

### 5.3 Radial isothermal Parker reference

A useful analytical wind reference assumes steady radial flow, spherical area expansion, constant isothermal pressure coefficient `ciso^2=p/rho`, and solar gravity. Then

\[
\rho U r^2=\mathcal M,
\]

\[
\left(U-\frac{c_{\rm iso}^2}{U}\right)\frac{dU}{dr}
=\frac{2c_{\rm iso}^2}{r}-\frac{GM_\odot}{r^2}.
\]

The critical radius is `rc=GM_sun/(2 ciso^2)`. The transonic solution obeys

\[
y-\ln y=4\ln(r/r_c)+4r_c/r-3,\qquad
y=(U/c_{\rm iso})^2.
\]

Choose the subsonic branch below `rc`, the sonic limit at `rc`, and the supersonic branch above it. The isothermal wind sound parameter differs from the adiabatic wave speed used at the shock:

\[
c_{s,1}^2=\gamma p_1/\rho_1.
\]

Maintaining an isothermal wind implies an ambient heating/energy prescription. It does not make the shock isothermal. This radial reference must not be described as a complete wind solution through a nonradial PFSS closed-field region.

### 5.4 Optional field-aligned open-tube reference

For arclength `s` along an open tube, cross-sectional area `A(s)`, field magnitude `B(s)>0`, and isothermal coefficient `ciso`, an elementary prescribed-tube model has

\[
AB=\Phi_{\rm tube},\qquad \rho U_\parallel A=\dot M_{\rm tube},
\]

\[
\left(U_\parallel-\frac{c_{\rm iso}^2}{U_\parallel}\right)
\frac{dU_\parallel}{ds}
=c_{\rm iso}^2\frac{d\ln A}{ds}-\frac{d\Phi_g}{ds},
\qquad \Phi_g=-GM_\odot/r.
\]

Select the physical outward direction independently of magnetic polarity. A negative radial magnetic polarity does not imply a wind directed toward the Sun. Multiple critical points, nulls, and separatrices require an explicit supported solution strategy; a failed solve is not repaired by a speed or density floor.

Closed loops can use a declared hydrostatic reference with no through-flow. For constant `ciso`,

\[
p(s)=p(s_0)\exp[-(\Phi_g(s)-\Phi_g(s_0))/c_{\rm iso}^2].
\]

These constructions can satisfy along-field force balance. They do not by themselves establish transverse force balance for arbitrarily assigned tube pressures in a potential magnetic field.

### 5.5 Ambient consistency diagnostics

Record positivity and finite-value checks, EOS consistency, mass-flux behavior where applicable, plasma beta, field/flow alignment, and the ambient's stated equilibrium or kinematic role. If the model claims a stationary ideal magnetic field, examine the induction condition rather than combining arbitrary velocity and field profiles:

\[
\partial_t\mathbf B-\nabla\times(\mathbf U\times\mathbf B)=0.
\]

A prescribed ambient approximation need not be a fully force-balanced global solution, but its discrepancy and validity range must be explicit. Local RH closure cannot correct a poor ambient density or field normalization.

### 5.6 Temperature, density, and wind sensitivity

An isothermal Parker profile is a useful controlled reference, but a temperature calibrated in the low corona must not automatically be used unchanged at 1 AU. In pure fully ionized hydrogen with `Tp=Te=T` and shock `gamma=5/3`,

\[
c_s=\sqrt{2\gamma k_BT/m_p}.
\]

It is approximately 196.27 km/s at 1.4 MK and 52.45 km/s at 0.1 MK. This change does not imply a universal factor-of-four error in fast Mach number: the Alfvén speed, field orientation, upstream wind, and shock-relative speed also matter.

For event work, accept validated empirical density/temperature histories or a qualified nonisothermal wind model. Define their open/closed-field routing, normalization, interpolation, extrapolation limits, and uncertainties. Open-tube empirical wind prescriptions are admissible reduced inputs when their kinematic role is declared; they do not establish global MHD force balance.

A possible polytropic radial wind uses

\[
p=K\rho^{\gamma_{\rm eff}},\qquad
c_{\rm eff}^2=\gamma_{\rm eff}p/\rho,\qquad
\rho Ur^2=\mathrm{constant},
\]

\[
\left(U-\frac{c_{\rm eff}^2}{U}\right)\frac{dU}{dr}
=\frac{2c_{\rm eff}^2}{r}-\frac{GM_\odot}{r^2}.
\]

Its effective wind exponent is distinct from the local shock adiabatic index. It requires its own boundary conditions and supported transonic solution. Substituting a declining temperature into an unchanged isothermal wind is an empirical kinematic profile, not this polytropic solution; record the momentum/energy discrepancy accordingly.

At minimum, publish `rho`, total `p`, `T`, `U`, `cs`, `vA`, `cf`, beta, and fast-speed margin along the apex and selected observer connections. Identify the location of extrema instead of assuming a universal low-coronal Alfvén-speed peak. Vary density, temperature, wind, and magnetic amplitude jointly when their calibration is correlated. Section 18.5 supplies a reproducible orientation-aware sensitivity example.

## 6. Coronal magnetic field

### 6.1 Recommended baseline: magnetogram-driven PFSS

The potential-field source-surface model solves [R2]

\[
\nabla\cdot\mathbf B=0,\qquad
\nabla\times\mathbf B=0,\qquad
\mathbf B=-\nabla\Psi,\qquad \nabla^2\Psi=0
\]

between `R_sun` and a declared source-surface radius `Rss`. Its boundary conditions are

\[
B_r(R_\odot,\theta,\phi)=B_{r,\rm map}(\theta,\phi),
\]

\[
B_\theta(R_{ss})=B_\phi(R_{ss})=0.
\]

`Rss` is a model parameter, not a universal observed surface. The commonly used 2.5 solar radii is an example requiring sensitivity assessment. The magnetic source surface is distinct from the front's coronal/SWCME handoff radius, which may be much larger.

PFSS provides a useful ambient reconstruction with open and closed topology. It omits coronal currents and the eruptive distortion caused by the modeled CME. A reconstruction driven by an event's magnetic map is more informative than a synthetic dipole, but is still an approximation with input and topology uncertainty.

### 6.2 Harmonic representation and normalization

Use a real orthonormal harmonic basis `Y_lm`, with all phase and cosine/sine conventions declared. Let

\[
B_{r,\rm map}=\sum_{l=1}^{l_{\max}}\sum_m b_{lm}Y_{lm}.
\]

The coefficients `b_lm` are photospheric **radial-field amplitudes in tesla**, not potential amplitudes. Defining `z=r/R_sun` and `q=Rss/R_sun`, an explicit radial representation is

\[
\Psi_{lm}(r,\theta,\phi)
=-\frac{b_{lm}R_\odot}{D_l}
\left[z^l-q^{2l+1}z^{-(l+1)}\right]Y_{lm},
\]

\[
D_l=l+(l+1)q^{2l+1}.
\]

Therefore,

\[
B_{r,lm}=\frac{b_{lm}}{D_l}
\left[lz^{l-1}+(l+1)q^{2l+1}z^{-(l+2)}\right]Y_{lm},
\]

\[
B_{\theta,lm}=\frac{b_{lm}}{zD_l}
\left[z^l-q^{2l+1}z^{-(l+1)}\right]\partial_\theta Y_{lm},
\]

\[
B_{\phi,lm}=\frac{b_{lm}}{zD_l\sin\theta}
\left[z^l-q^{2l+1}z^{-(l+1)}\right]\partial_\phi Y_{lm}.
\]

These formulas make both boundary conditions explicit. At high degree, evaluate rescaled forms to avoid overflow in powers of `q`. At poles, use analytical Cartesian limits or a regular chart, rather than replacing `sin(theta)` by an arbitrary floor.

For the conventional orthonormal `Y_10=sqrt(3/(4 pi)) cos(theta)`, a pure dipole with polar photospheric field `Bpole` has

\[
b_{10}=B_{\rm pole}\sqrt{4\pi/3}.
\]

This is a verification field, not an event-specific reconstruction. Independent harmonic and flux tests are available in [R3].

### 6.3 Magnetic-map preprocessing

The resolved magnetic asset must retain:

- Observatory/instrument, product, observation time range, map version, and checksum.
- Whether the supplied map is radial field, line-of-sight field, or a derived product.
- Coordinate system, longitude convention, latitude sampling, and solid-angle weights.
- The transformation from the map frame to HCI at the reference epoch.
- Harmonic degree, normalization, filtering/apodization, and treatment of poorly observed regions.
- Signed flux, unsigned flux, reconstruction RMS, removed mean field if any, and the normalization/calibration operation.

Reject a significant monopole by default. If an explicit mean-removal policy is selected, record the removed component and its effect. A balanced map can have nonzero unsigned open flux; this does not authorize adding a global magnetic monopole.

A single explicitly documented field-amplitude calibration can be applied to the retained raw coefficients. Do not repeatedly scale already calibrated coefficients. Calibrating `|B|` at 1 AU and calibrating its radial component are different operations because spiral winding contributes to `|B|`.

### 6.4 Field derivatives, nulls, and topology

Return the Cartesian field, its valid derivative representation if needed, field-source identity, and reconstruction generation. Numerical divergence tests must differentiate the same field interpolation actually used by production queries.

At a genuine magnetic null, `B` can be a valid zero vector while `B/|B|`, field-line direction, and magnetic obliquity are undefined. Return explicit validity flags for these diagnostics. Do not insert a fictitious magnetic field floor to force a direction. The local fluid-wave calculation can approach its hydrodynamic limit, but magnetic connectivity through the null is unresolved.

Open/closed classification uses integration with declared stopping conditions. Distinguish a closed loop from a trace that hit a null, exited the available domain, exceeded a step limit, or failed numerically. Those outcomes have different scientific meanings.

### 6.5 Flux-consistent amplitude calibration

Treat low-coronal magnetic strength, open-field area, and heliospheric radial/unsigned open flux as joint constraints. The known open-flux discrepancy motivates an uncertainty study, rather than a universal factor-of-two rescaling [R10]. A global constant multiplier of a solenoidal field preserves its divergence and field-line topology, but also changes its low-coronal strength; report the effect at both ends of the model.

Independent scalar multipliers below and above a radius generally violate the magnetic contract. If `div(B)=0`, then

\[
\nabla\cdot[f(r)\mathbf B]=f'(r)B_r.
\]

A jump in the multiplier also creates a discontinuity of normal magnetic flux. A smooth radial blend is therefore not a repair when `Br` is nonzero. Extra calibration freedom requires a new flux/potential/boundary construction, checked for divergence and interface flux continuity, rather than two unrelated amplitude factors.

Changing `Rss` modifies the open-field distribution and unsigned open flux. It is a physical model variation with connectivity consequences, not just an amplitude knob. Preserve the raw map, fitted operations, calibration observables, uncertainties, and independent validation constraints in the manifest. Field strengths inferred from shock standoff or radio band splitting carry their own assumptions and are not independent validation if they supplied the fit.

## 7. Connecting the ambient field to the heliosphere

### 7.1 Required matching properties

Use the existing canonical corona/heliosphere matching construction when independently validated. The composite must preserve:

- Magnetic polarity and signed flux.
- Field-line routing and the source-footpoint coordinate convention.
- Continuous field values where continuity is part of the selected boundary contract.
- Solenoidality in each smooth region and normal-field continuity at interfaces.
- The claimed relation between field evolution and ambient velocity.

Componentwise interpolation between two unrelated vector fields generally does not preserve `div(B)=0`. A smooth transition must use a validated flux/vector-potential construction or a derived mapping, and must be checked independently. A current sheet or separatrix requires one-sided evaluation and explicit topology handling.

### 7.2 An explicit solenoidal exterior reference

The following is a useful derivation and verification option; it is not a claim about the current canonical provider. Assume:

1. A radial wind speed `U_r(r)>0` depending only on radius outside `Rss`.
2. Rigid rotation of the source map with angular speed `Omega` about the solar axis.
3. No latitudinal flow and `B_theta=0` outside `Rss`.
4. An HCI azimuthal wind `U_phi=Omega Rss^2 sin(theta)/r`, a prescribed angular-momentum continuation of corotation at `Rss`.

Let `G(theta,phi0)=Rss^2 Br(Rss,theta,phi0,tref)` and define

\[
f(r)=\int_{R_{ss}}^r
\frac{\Omega[1-R_{ss}^2/r'^2]}{U_r(r')}\,dr',
\qquad
\phi_0=\phi-\Omega(t-t_{\rm ref})+f(r).
\]

Then set

\[
B_r=\frac{G(\theta,\phi_0)}{r^2},\qquad
B_\theta=0,\qquad
B_\phi=-r\sin\theta f'(r)B_r.
\]

The radial and azimuthal divergence terms cancel:

\[
\frac{1}{r^2}\partial_r(r^2B_r)
=\frac{G_{,\phi_0}f'}{r^2},
\]

\[
\frac{1}{r\sin\theta}\partial_\phi B_\phi
=-\frac{G_{,\phi_0}f'}{r^2}.
\]

Since `f'(Rss)=0`, the exterior matches the source surface's radial field and zero tangential field. At large radius it approaches the familiar spiral relation

\[
\frac{B_\phi}{B_r}\simeq-\frac{\Omega r\sin\theta}{U_r}.
\]

For constant `U_r=U0`,

\[
f(r)=\frac{\Omega}{U_0}
\left[r-2R_{ss}+\frac{R_{ss}^2}{r}\right].
\]

This construction explicitly accounts for longitude structure; an arbitrary longitude-dependent radial field combined with an axisymmetric spiral formula can violate solenoidality. The assumed azimuthal wind is kinematic, not a solved magnetic-torque model. Verify induction with the stated rotating map and velocity if using this reference as a stationary-in-corotation ideal field.

Field values are continuous at `Rss`, but all spatial derivatives need not be. Do not claim a `C1` magnetic transition without checking it or implementing an additional validated transition layer. Connectivity studies using PFSS and heliospheric spirals have observational precedents [R1].

### 7.3 What qualifies as realistic coronal magnetism

For an event study, the model must use an appropriate observed magnetic map and report sensitivity to map age, harmonic resolution, field normalization, and `Rss`. Validate field-line connectivity and, where possible, magnetic strength or open-flux constraints independently.

The term realistic applies to an assessed ambient reconstruction. It does not mean that PFSS supplies coronal free energy, eruption-generated current sheets, or the post-shock/ejecta magnetic field. Synthetic fixtures and observed reconstructions must have different manifest labels.

## 8. Front geometry, normals, and speeds

### 8.1 Implicit surface formulation

A supported surface is defined by

\[
g(\mathbf x,t)=0,
\qquad \mathbf x\in\mathcal S_{\rm supported}(t).
\]

Differentiating along a point moving with the surface gives

\[
\partial_tg+\mathbf V_s\cdot\nabla g=0.
\]

Consequently,

\[
\boxed{V_n=-\frac{\partial_tg}{|\nabla g|}}.
\]

Only the normal component is geometrically determined. Tangential motion of surface labels can be a reparameterization and is not plasma motion. In the jump solve, choose surface-frame velocity `Vn n`; this is sufficient because tangential Galilean choices do not change the physical jump.

A small dimensionless level-set residual does not by itself bound position error. Also report a distance-like residual `|g|/|grad(g)|` near regular surface points.

### 8.2 Finite triaxial ellipsoid

Let `C(t)` be the center, `a_i(t)>0` the semiaxes, and `R(t)` a proper orthogonal rotation mapping body coordinates into HCI. Define

\[
Q=R\,\operatorname{diag}(a_1^{-2},a_2^{-2},a_3^{-2})R^T,
\qquad \mathbf y=\mathbf x-\mathbf C,
\]

\[
g=\mathbf y^TQ\mathbf y-1.
\]

At fixed Cartesian position,

\[
\nabla g=2Q\mathbf y,
\]

\[
\partial_tg=\mathbf y^T\dot Q\mathbf y
-2\dot{\mathbf C}^{\,T}Q\mathbf y,
\]

\[
\boxed{
V_n=\frac{2\dot{\mathbf C}^{\,T}Q\mathbf y
-\mathbf y^T\dot Q\mathbf y}{2|Q\mathbf y|}.
}
\]

`dot(Q)` includes both expansion and changing attitude. Omitting the rotation derivative yields incorrect normal speeds for a rotating ellipsoid.

The full mathematical ellipsoid is not automatically the physical supported front. Define an explicit leading patch or finite extent, and exclude artificial rear closures and points inside the Sun. Surface quadrature uses the actual Cartesian area element

\[
dA=|\partial_q\mathbf X\times\partial_p\mathbf X|\,dq\,dp.
\]

Use regular charts or a triangulated surface for pole-safe area integration.

Time-dependent semiaxes and attitude already permit changes within the ellipsoid family. They do not predict arbitrary ambient-driven deformation or corrugation. Neither an ellipsoid nor an SSE shape is a universally preferred event model; select the family that fits the observed front and disclose unresolved shape discrepancy.

### 8.3 Finite SSE-inspired front

A simple self-similar reference follows the circular-front geometry of SSE [R4], extended here by an explicit axisymmetric spherical leading surface. The three-dimensional extension and its finite support are model assumptions.

Let `R_a(t)` be the apex heliocentric distance, `d` a unit propagation direction, and `lambda` initially a constant angular half-width with `0<lambda<pi/2`. Define

\[
c(t)=\frac{R_a(t)}{1+\sin\lambda},\qquad
a(t)=\frac{R_a(t)\sin\lambda}{1+\sin\lambda},
\qquad \mathbf C=c\mathbf d.
\]

The enclosing sphere has radius `a` and center `C`. For a heliocentric ray `e` making angle `psi=acos(e dot d)`, the supported leading intersection is

\[
r_+(\psi,t)=c\cos\psi+
\sqrt{a^2-c^2\sin^2\psi},\qquad 0\le\psi\le\lambda.
\]

The smaller root is a rear intersection and is not included in the leading-front support. The radicand and angular condition are checked with controlled roundoff handling; a genuinely negative radicand means no intersection.

For fixed direction and width, `r_+=R_a K(psi,lambda)`. The radial intersection speed `dot(r_+)=dot(R_a) K` is finite at a fixed tangent ray; angular derivatives and tracking a moving tangency are ill-conditioned there. Evaluate normals and normal speeds from the Cartesian implicit surface, not an angular difference across the support boundary. With `s=sin(lambda)`,

\[
\boxed{V_n=\dot R_a\frac{\mathbf d\cdot\mathbf n+s}{1+s}}.
\]

At the tangent flank, `d dot n=-s` and `e dot n=0`, giving `Vn=0`. A purely radial wind also has `U dot n=0` there, so `w1=0`. Nonradial wind does not generally have that limit. **The geometric half-width is not the shock half-width.** Local fast-mode and RH tests determine the accepted shock patch; changing the width merely to remove sub-fast samples would change the physical input.

As an analytical reference only, set wind to zero and take uniform, orientation-independent `cf`. If `Ma=dot(R_a)/cf`, the super-fast angular boundary satisfies

\[
\frac{\sqrt{s^2-\sin^2\psi}
\left(\cos\psi+\sqrt{s^2-\sin^2\psi}\right)}{s(1+s)}
=\frac{1}{M_a}.
\]

The boundary half-angle in degrees is:

| Geometric half-width | `Ma=1.5` | `Ma=2` | `Ma=3` | `Ma=5` |
| --- | ---: | ---: | ---: | ---: |
| 30° | 19.1 | 23.3 | 26.6 | 28.6 |
| 45° | 26.1 | 32.4 | 37.8 | 41.6 |
| 60° | 31.1 | 39.2 | 46.8 | 52.7 |

These values are not a formula for a structured magnetized corona, where `cf` depends on position and field orientation. Report the actual accepted footprint, which can be asymmetric or disconnected.

The `lambda -> 0` limit degenerates to a point-like front and is not a regular three-dimensional shock surface. The `lambda -> pi/2` limit is a Sun-tangent sphere with `c=a=R_a/2`; it is not a heliocentric sphere of radius `R_a`. An independent Sun-centered spherical reference is a different geometry.

SSE geometry prescribes an expanding shape. It does not establish a plasma density distribution, material contact, or a dynamically solved shock. Smooth width/direction histories are a bounded extension of the same family. Their complete derivatives are

\[
\dot c=\frac{\dot R_a}{1+s}
-\frac{R_a\cos\lambda\,\dot\lambda}{(1+s)^2},\qquad
\dot a=\frac{\dot R_a s}{1+s}
+\frac{R_a\cos\lambda\,\dot\lambda}{(1+s)^2},
\]

\[
\boxed{
V_n=\dot R_a\frac{\mathbf d\cdot\mathbf n+s}{1+s}
+\frac{R_a\cos\lambda\,\dot\lambda(1-\mathbf d\cdot\mathbf n)}
{(1+s)^2}
+c\,\dot{\mathbf d}\cdot\mathbf n .
}
\]

Require a unit direction with \(\mathbf d\cdot\dot{\mathbf d}=0\), differentiable histories, and `0<lambda(t)<pi/2` throughout the validity interval. With fixed direction but changing width, the tangent value is `Vn=R_a cos(lambda) dot(lambda)/(1+sin(lambda))`; the zero-speed edge statement above applies to fixed width and direction. Match width, direction, and their derivatives through handoff, or reject an outer kernel that cannot represent them.

### 8.4 Support and near-Sun clipping

For each sample, distinguish:

- Mathematical surface point.
- Supported leading-front point.
- Point above the declared physical inner radius.
- Valid ambient query.
- Accepted fast-shock jump.

A finite cap's perimeter is a support boundary, not a new side shock or an impermeable wall. Edge-normal and crossing diagnostics must not be obtained by differentiating a discontinuous support mask. Keep an edge flag and use the parent surface's one-sided geometric derivatives.

At initial epochs a mathematical surface can intersect the photosphere. Only the supported part outside the physical inner boundary is evaluated. This is prescribed near-Sun geometry, rather than a proof of physical shock creation at the solar surface.

## 9. Prescribed coronal motion and outer propagation

### 9.1 History ownership and interpolation

A complete history supplies center, semiaxes, direction/attitude, and their derivatives, or supplies apex distance and parameters from which those quantities are derived. The history has a closed validity interval and a documented interpolation rule.

Tabulated values must be strictly time-ordered and finite. Duplicate or contradictory timestamps are rejected. Positive semiaxes, positive supported radius, and any selected monotonic-outward condition are checked between knots, not just at the knots. Interpolation overshoot can violate these conditions even when every tabulated value is valid.

Prefer a derivative-consistent Hermite or constrained smooth representation. Rotation interpolation must preserve proper orthogonality and resolve quaternion sign ambiguity; if normal speed requires angular velocity, provide a smooth orientation history and its derivative. Differentiating a discontinuous attitude sequence is unsupported.

No extrapolation is allowed by default. Any selected extrapolation law is a physical input with a finite authorized interval, not an implicit consequence of using the final tabulated speed.

### 9.2 An analytical acceleration fixture

The following history is a reproducible prescribed launch profile. It is not a calculated Lorentz-force eruption. For elapsed time `t` from initialization, acceleration duration `tau_a>0`, initial apex radius `R0`, and speeds `V0` and `Vf`, let

\[
s=t/\tau_a,\qquad \Delta V=V_f-V_0,
\]

\[
h(s)=10s^3-15s^4+6s^5,\qquad 0\le s\le1.
\]

During the pulse,

\[
V_a(t)=V_0+\Delta Vh(s),
\]

\[
R_a(t)=R_0+V_0\tau_as
+\Delta V\tau_a\left(\frac52s^4-3s^5+s^6\right),
\]

\[
\dot V_a(t)=\frac{30\Delta V}{\tau_a}s^2(1-s)^2.
\]

At `t=tau_a`,

\[
R_a=R_0+\frac{\tau_a}{2}(V_0+V_f),\qquad V_a=V_f.
\]

If the next prescribed coronal segment is constant-speed propagation, use this endpoint and `Vf`. The acceleration vanishes at both pulse endpoints. The profile is valid from its specified start; it does not prescribe a pre-eruption state.

If semiaxes scale with `R_a`, their derivatives must scale consistently. If width, direction, or aspect ratio changes independently, those histories are additional explicit inputs.

### 9.3 Outer SWCME propagation

Map the matched front state into the existing SWCME outer propagation kernel whenever it can represent the selected geometry. Reuse the kernel; do not create a second set of front kinematics merely to avoid an adapter.

If the kernel represents only a restricted shape family, accept only compatible histories or implement and qualify a documented generalization. Do not silently project a rotating triaxial front onto a fixed-width sphere while retaining the old geometry label.

One possible outer kinematic reference is quadratic drag,

\[
\frac{dR_a}{dt}=V_a,\qquad
\frac{dV_a}{dt}=-\Gamma(V_a-U_w)|V_a-U_w|,
\]

where `Gamma>=0` has units `m^-1`. For constant `Gamma>0` and wind `Uw`, with matched values at `th` and `Delta t=t-th`, define `z0=Vh-Uw`. Then

\[
V_a=U_w+\frac{z_0}{1+\Gamma|z_0|\Delta t},
\]

\[
R_a=R_h+U_w\Delta t
+\frac{\operatorname{sgn}(z_0)}{\Gamma}
\ln[1+\Gamma|z_0|\Delta t].
\]

Use the smooth limits for `z0=0` and `Gamma=0`, and `log1p` for small arguments. A radius-dependent or time-dependent wind/drag coefficient requires a numerical trajectory solve and independent timestep convergence.

Drag-based CME propagation has relevant precedents [R5, R6]. Applying a CME-body drag law to a **shock front** is an additional proxy approximation. It must be named, calibrated against front observations where available, and assessed separately from body-motion prediction. A local RH solution does not prove that this trajectory is dynamically generated by a CME driver.

The drag wind input should use the canonical upstream wind when the law is intended to represent interaction with that wind. If it is an independently fitted effective value, label it as such and record the discrepancy; it cannot silently replace the wind used for Mach-number calculation.

### 9.4 Shock-front proxy validity

For a freely propagating weak outward fast shock, the local normal propagation limit is

\[
V_n\longrightarrow\mathbf U_1\cdot\mathbf n+c_{f,1}
\quad\text{as}\quad M_{f,1}\longrightarrow1^+.
\]

The body-drag reference instead approaches its trajectory wind. It therefore does not guarantee a persistent shock, even if position and velocity remain numerically smooth. A front point with `Vn<U dot n` is overtaken by the local flow; a point with positive inflow can still be sub-fast. Export

\[
q_f=V_n-\mathbf U_1\cdot\mathbf n-c_{f,1}
\]

alongside the local classification. Positive `qf` is a super-fast candidate, not proof that a supported RH branch has been accepted.

Use the canonical wind for a physical wind-dependent trajectory by default. A fitted effective wind remains an explicitly labeled proxy input. Export its difference from the canonical apex wind every epoch. A difference is not automatically a failure: evaluate the actual fast-speed margin and declared model validity. Never alter the canonical wind, clip speed, or floor Mach number to rescue the drag trajectory.

The geometric fixture profile uses `report_only` validity: complete the specified trajectory and report where it ceases to support a shock. An event profile can require declared shock coverage at selected times/observers. Failure of that physical coverage requirement rejects the qualification; it must not be converted into a successful shock arrival merely because the radius endpoint was reached. A proxy's loss of admissibility is not evidence that the actual event's shock disappeared.

A constant-target variant approaching `S0=U_n+cf` can be written with the same analytical drag formulas by replacing `Uw` by `S0` and using positive excess `delta0=Vh-S0`. It is a new phenomenological law, not a correction already justified by CME-body drag. In a one-dimensional radial apex reference with variable `S(R,t)=Ur+cf`,

\[
\delta=V-S(R,t),\qquad
\dot\delta=\dot V-\partial_tS-V\partial_RS.
\]

Choosing `dot(delta)=-Gamma delta |delta|` would therefore require

\[
\dot V=\partial_tS+V\partial_RS-\Gamma\delta|\delta|.
\]

This optional research closure requires its own physical justification, numerical integration, convergence, and event calibration. A value of `S` sampled at the current apex cannot be substituted into a constant-target closed form. Neither variant supplies global driver/shock dynamics or guarantees flank admissibility.

## 10. Continuous corona–SWCME handoff

### 10.1 Distinct interfaces

There are two independent transition concepts:

1. **Ambient magnetic/plasma matching**, for example near a PFSS source surface.
2. **Front evolution handoff**, for example at an apex radius near 20 solar radii.

They do not have to occur at the same radius. Changing the front propagation kernel is not a reason to change the upstream plasma or field model.

### 10.2 Matched state

Find the handoff time as a bracketed event root,

\[
R_a(t_h)-R_h=0.
\]

For a nonmonotone history, declare whether the first outward crossing is intended. Store the resolved crossing time and the complete matched state:

- Apex distance and speed.
- Center and center velocity.
- Semiaxes and their derivatives.
- Attitude and angular velocity where supported.
- Direction, width, aspect ratios, and support-mask definition.
- Ambient model identity and input fingerprints.
- History phase, generation, and physical-quantity identity such as `shock_front` rather than `ejecta_leading_edge`.

The outer kernel begins from this state. It does not reread an unrelated outer reference radius or reset the elapsed time to a fresh event origin.

### 10.3 Required continuity

At minimum, establish continuity of supported front position and normal speed across handoff. Compatible parameter derivatives provide a `C1` front history. Check complete normal and shape continuity over the surface, rather than only apex radius/speed.

Acceleration continuity is required only if it is part of the declared contract or used by downstream diagnostics. A sharp change in the trajectory law can preserve position and speed while changing acceleration. Publish this limitation instead of claiming a smoother match.

For the constant-speed-to-drag examples, the one-sided apex accelerations at handoff are 0 and `-Gamma(Vh-Uw)|Vh-Uw|=-7.2 m/s^2`. Record both limits and the jump. This is a legitimate `C1` analytical fixture, not a `C2` claim. A smooth transition can be required for an event history or derivative-based diagnostic, but must be qualified as a changed history rather than advertised as universally more physical.

If a transition interval is used, blend a differentiable **history**, not separately calculated position and speed. For two compatible scalar trajectories and a smooth weight `w(t)`,

\[
R=(1-w)R_c+wR_o,
\]

\[
V=(1-w)V_c+wV_o+\dot w(R_o-R_c),
\]

\[
\dot V=(1-w)\dot V_c+w\dot V_o
+2\dot w(V_o-V_c)+\ddot w(R_o-R_c).
\]

The weight-derivative terms are essential. Use a validated matrix/rotation construction for attitude; entrywise averaging of rotation matrices is not a proper rotation. Reject nonphysical shape overshoot or inward motion introduced by blending.

### 10.4 What must remain unchanged at handoff

The upstream provider remains the same composite ambient authority. The same shock solver, EOS, field-line tracer, and support definition remain in use. Changing phase must not cause discontinuous Mach number or compression when the complete physical input state is continuous.

Any apparent discontinuity is traced to input kinematics, an ambient interface, an actual change in admissibility, or a numerical defect. It is not hidden by a compression floor or independent interpolation of jump outputs.

## 11. Local shock physics

### 11.1 Wave speeds and magnetic obliquity

For a positive upstream density and pressure,

\[
c_s^2=\frac{\gamma p_1}{\rho_1},\qquad
v_A^2=\frac{|\mathbf B_1|^2}{\mu_0\rho_1},\qquad
v_{A,n}^2=\frac{(\mathbf B_1\cdot\mathbf n)^2}{\mu_0\rho_1}.
\]

The fast magnetosonic speed along the front normal is

\[
c_f^2=\frac12\left[c_s^2+v_A^2+
\sqrt{(c_s^2+v_A^2)^2-4c_s^2v_{A,n}^2}\right].
\]

The slow-mode speed uses the minus sign before the square root. A slightly negative discriminant consistent with arithmetic roundoff can be handled with a documented bound; a physically invalid negative discriminant is an error.

The orientation limits must be tested explicitly:

\[
c_{f,\parallel}=\max(c_s,v_A),\qquad
c_{f,\perp}=\sqrt{c_s^2+v_A^2}.
\]

A polar dipole field parallel to a radial apex normal uses the first formula. Using the perpendicular expression there overestimates `cf` and biases `Mf` downward. At a magnetic null the fast speed has the hydrodynamic limit `cs`; this does not make magnetic direction or obliquity defined.

The positive inflow and fast Mach number are

\[
w_1=V_n-\mathbf U_1\cdot\mathbf n,\qquad M_f=w_1/c_f.
\]

When `B1` has a resolved nonzero direction, publish both a polarity-aware cosine and an acute shock obliquity:

\[
\mu_{Bn}=\frac{\mathbf B_1\cdot\mathbf n}{|\mathbf B_1|},
\qquad
\theta_{Bn}=\cos^{-1}(|\mu_{Bn}|)\in[0,\pi/2].
\]

Use bounded roundoff handling when evaluating `acos`. The acute angle deliberately removes magnetic polarity; preserve `mu_Bn` for consumers that need it. At a null the obliquity is undefined.

Also publish plasma beta, sonic Mach number, total-field Alfvén Mach number, and normal-field Alfvén Mach number when their denominators are meaningful:

\[
\beta_1=\frac{2\mu_0p_1}{B_1^2},\qquad
M_s=\frac{w_1}{c_s},\qquad
M_A=\frac{w_1}{v_A},\qquad
M_{A,n}=\frac{w_1}{|v_{A,n}|}.
\]

These Mach numbers are not interchangeable. Neither a universal `Mf=3` acceleration threshold nor an injection efficiency is part of this background-only model.

### 11.2 Admissibility and physical interpretation

A candidate outward fast shock requires `w1>0` and upstream super-fast inflow. A maintained jump solver must additionally find the appropriate admissible compressive branch.

This is an **admissibility test for a prescribed discontinuity**, not a simulation of shock formation. A local condition `Mf>1` does not reconstruct the wave steepening that produced a real shock. Formation can occur after a compression wave has propagated away from its driver [R7].

Sub-fast samples remain valid geometric front samples. They have no accepted fast-shock downstream state. Do not create a shock by assigning a minimum compression. Do not automatically claim that a complete smooth compression-wave volume has been solved there.

### 11.3 Shock-frame variables

Choose a local frame translating at `Vn n`. For either side,

\[
\mathbf u_i=\mathbf U_i-V_n\mathbf n,
\qquad u_{i,n}=\mathbf u_i\cdot\mathbf n,
\]

\[
\mathbf u_{i,t}=\mathbf u_i-u_{i,n}\mathbf n,
\qquad B_{i,n}=\mathbf B_i\cdot\mathbf n,
\qquad \mathbf B_{i,t}=\mathbf B_i-B_{i,n}\mathbf n.
\]

For the forward-shock convention, both normal fluid velocities are negative. Tangential vector equations can be evaluated in Cartesian form, avoiding arbitrary basis-angle singularities.

Let `[F]=F2-F1`. The one-fluid isotropic ideal-MHD Rankine–Hugoniot relations are [R8]:

**Mass**

\[
[\rho u_n]=0.
\]

**Normal magnetic field**

\[
[B_n]=0.
\]

**Tangential electric field**

\[
[u_n\mathbf B_t-B_n\mathbf u_t]=0.
\]

**Normal momentum**

\[
\left[\rho u_n^2+p+\frac{B_t^2}{2\mu_0}\right]=0.
\]

The full normal Maxwell stress also contains `-Bn^2/(2 mu0)`, which cancels because `Bn` is continuous.

**Tangential momentum**

\[
\left[\rho u_n\mathbf u_t-\frac{B_n\mathbf B_t}{\mu_0}\right]=0.
\]

**Total energy flux**

\[
\left[
u_n\left(\frac12\rho u^2+\frac{\gamma}{\gamma-1}p
+\frac{B^2}{\mu_0}\right)
-\frac{B_n}{\mu_0}(\mathbf u\cdot\mathbf B)
\right]=0.
\]

Gravity, curvature, and slow ambient gradients do not enter the ideal infinitesimal jump directly. They matter to the upstream state and large-scale trajectory. There is no extra heat-source term at a shock in this closure; the energy-conserving dissipative jump raises entropy.

### 11.4 Compression and downstream reconstruction

Define the density compression

\[
X=\rho_2/\rho_1>1.
\]

Mass conservation gives

\[
u_{2,n}=u_{1,n}/X,
\qquad U_{2,n}=V_n-\frac{w_1}{X}.
\]

For a trial `X`, solve tangential momentum and electric-field continuity for `B2t` and `u2t`, then use normal momentum to obtain

\[
p_2=p_1+\rho_1w_1^2(1-1/X)
+\frac{B_{1,t}^2-B_{2,t}^2}{2\mu_0}.
\]

The remaining energy-flux equation determines the physical compression root. For a regular coplanar oblique case with nonzero `Bn`, the tangential relation can be written

\[
\mathbf B_{2,t}
=X\frac{M_{A,n}^2-1}{M_{A,n}^2-X}\mathbf B_{1,t}.
\]

This formula is a diagnostic derivation, not a robust universal numerical implementation: perpendicular, parallel, switch-on/off, and near-singular cases require dedicated branch handling. In particular,

\[
B_{2,n}=B_{1,n}
\]

for every accepted ideal jump, but `B2t=X B1t` applies to a perpendicular shock (`Bn=0`) and certain limits, not every oblique shock. Multiplying the full magnetic vector by `X` is incorrect when `Bn` is nonzero.

For a common one-temperature fully ionized composition,

\[
T_2/T_1=(p_2/p_1)/X.
\]

If the provider retains separate electron and ion temperatures, the jump determines total pressure only. It must not invent a heating partition without an additional declared model.

### 11.5 Accepted branch checks

An accepted result must pass:

- Positive finite density and pressure, finite velocities and fields.
- A nontrivial compressive root in the maintained solver's physical interval.
- Entropy increase, measured for constant composition/gamma by `(p/rho^gamma)2/(p/rho^gamma)1`.
- All independent RH residuals within declared tolerances.
- Fast-family/evolutionary branch checks appropriate to the regular or degenerate geometry.
- Conditioning and reconstruction checks near tangential singularities.

For the ordinary ideal-gas fast-shock family used here, the strong-shock hydrodynamic compression limit is `(gamma+1)/(gamma-1)`, equal to 4 at `gamma=5/3`. This limit is not a command to set all strong or coronal shocks to compression 4. Physics outside this one-fluid closure can require a different model.

Strict branch inequalities become degenerate at exact parallel/perpendicular limits. Use the maintained solver's verified limiting treatment rather than dividing by a vanishing normal or tangential field component. A valid conservative state on the wrong branch is not an accepted forward fast shock.

### 11.6 Independent reference examples

**Hydrodynamic/weak-parallel-field example.** With `gamma=5/3`, sonic Mach `Ms=3`, and a field too weak to alter the gas-dynamic branch,

\[
X=\frac{(\gamma+1)M_s^2}{(\gamma-1)M_s^2+2}=3,
\]

\[
p_2/p_1=\frac{2\gamma M_s^2-(\gamma-1)}{\gamma+1}=11.
\]

If `Vn=1200 km/s` and `U1n=300 km/s`, then `w1=900 km/s`, `U2n=900 km/s`, and `T2/T1=11/3`. A zero-field reference has undefined magnetic obliquity; this must not be encoded as a measured zero angle.

**Perpendicular finite-beta example.** For `gamma=5/3` and `Bn=0`, eliminating the RH variables gives

\[
M_A^2=\frac{X(X+5+5\beta_1)}{2(4-X)}.
\]

With `beta1=0.2` and `MA=2`, the physical compression is `X=2`. The exact tangential field ratio is 2 and the pressure ratio is 6. This is a finite-positive-pressure reference, avoiding an artificial cold-plasma floor. Its fast Mach number is `2/sqrt(1+gamma beta1/2)`, about 1.852.

Use these independent relations as tests of the shared production solver. They do not represent a complete CME event.

### 11.7 Additional surface diagnostics

**Magnetic compression.** For an accepted jump with resolved nonzero upstream field, publish

\[
C_B=\frac{|\mathbf B_2|}{|\mathbf B_1|}
=\frac{\sqrt{B_{n,1}^2+|\mathbf B_{t,2}|^2}}
{\sqrt{B_{n,1}^2+|\mathbf B_{t,1}|^2}}.
\]

It is generally different from density compression. A perpendicular ideal-MHD jump has `CB=X`; a pure parallel field can have `CB=1`. At an upstream null or a rejected jump, this ratio is absent with a reason flag.

**de Hoffmann–Teller frame.** First use the shock frame with velocity `Vn n` and signed `u1n=-w1`. When `B1n` is resolved and nonzero, the tangential boost into the HT frame is

\[
\mathbf v_{{\rm HT},t}
=\mathbf u_{1,t}-\frac{u_{1,n}}{B_{1,n}}\mathbf B_{1,t},
\qquad
\mathbf u_{1,{\rm HT}}=\frac{u_{1,n}}{B_{1,n}}\mathbf B_1.
\]

Thus `u1_HT cross B1=0` and

\[
|\mathbf u_{1,{\rm HT}}|
=\frac{w_1}{|\cos\theta_{Bn}|}.
\]

This incident-flow speed is not the HT boost speed. Only in the normal-incidence frame, where `u1t=0`, is

\[
|\mathbf v_{{\rm HT},t}|=w_1|\tan\theta_{Bn}|.
\]

Name the exported frame, signed vector, and magnitude explicitly. If needed, the HCI velocity of this HT frame is `Vn n+vHT_t`. These Galilean relations support frame checks and later connection diagnostics [R14]; they are not an injection prescription. At `B1n=0` with nonzero inflow, no finite boost of this form exists. Near perpendicularity, report conditioning and an unsupported diagnostic when the declared division/velocity validity bounds are exceeded. Do not impose an obliquity floor, corrupt the otherwise valid jump, or label `w1/cos(theta)` as the tangential boost.

**First critical Mach number.** A collisionless first-critical-Mach diagnostic is optional. It requires an explicitly selected, verified procedure/table, its definition of Mach number, EOS, beta, obliquity, interpolation, uncertainty, and supported range. The critical value depends on these parameters [R13]; a universal constant such as 2.7 is not an acceptable general implementation. If a table uses a different Mach convention, convert consistently or reject the comparison. With `criticality_model=none`, publish absent criticality, rather than “subcritical.” When enabled, report `Mf/Mc1` only under a compatible fast-Mach convention. Crossing unity neither changes ideal-MHD jump admissibility nor enables particle reflection, injection, or acceleration.

## 12. Numerical treatment and failure semantics

### 12.1 Upstream sampling

Evaluate the canonical ambient at the front position when the provider directly represents an undisturbed reference. If using a gridded ambient product with a one-sided upstream stencil, specify the positive-normal offset/stencil and its extrapolation to the front.

The offset must be small relative to ambient-gradient and curvature scales, yet sufficient to avoid interpolation through a discontinuity already present in an imported product. Refine it independently. Never sample an interior reference and label it an actual downstream state.

### 12.2 Surface classification statuses

The conceptual statuses below are stable semantic requirements, not existing enum names:

| Status | Meaning | Downstream state |
| --- | --- | --- |
| `OUTSIDE_FRONT_SUPPORT` | Ray/point outside the finite leading patch | Absent |
| `BELOW_PHYSICAL_INNER_BOUNDARY` | Mathematical point lies inside the excluded solar region | Absent |
| `AMBIENT_UNAVAILABLE` | Requested upstream state is outside provider support | Absent; candidate rejected if required |
| `NON_FORWARD_INFLOW` | `w1<=0` under the selected sign convention | Absent |
| `SUBFAST_FRONT` | Valid geometric front, no accepted upstream super-fast condition | Absent |
| `SOLVED_FAST_SHOCK` | Accepted compressive fast branch and conservative jump | Valid local one-sided limit |
| `NUMERICALLY_UNRESOLVED_WEAK_SHOCK` | Physical/numerical classification cannot resolve a small super-fast jump reliably | Absent; typed numerical failure |
| `WRONG_BRANCH` | Conservative candidate belongs to an unsupported family | Absent; failure |
| `INVALID_JUMP` | Positive/finite/conservation/conditioning test failed | Absent; failure |

Absence of a shock is a physical classification; inability to solve an identified super-fast jump is a numerical/model failure. They must not share an undifferentiated `no shock` flag.

### 12.3 Weak shocks

The trivial no-jump root `X=1` exists algebraically. A solver must not return it as the compressive solution merely because a weak physical root lies nearby. Use the maintained deflated/bracketed formulation and verified special limits rather than an unconstrained root search from an arbitrary initial guess.

Publish upstream Mach excess, bracket, root iteration count, conditioning, and conservative residuals. Classification thresholds must reflect input uncertainty and arithmetic resolution, not a desired minimum source strength. Never manufacture `Mf>1` or `X>1` with a numerical floor.

The earlier `SHOCK_SOLVER_FAILURE` reports motivate coverage of near-perpendicular, parallel, oblique, and numerically unresolved weak shocks. This document does not claim that all such cases in the user's latest local source have been fixed.

### 12.4 Residual normalization

For each conservation flux `F`, use a declared dimensional-plus-relative tolerance,

\[
|F_2-F_1|\le\epsilon_{F,\rm abs}
+\epsilon_{F,\rm rel}\,F_{\rm scale}.
\]

The scale is constructed from finite physical flux magnitudes, with documented units and behavior when a component is zero. Use vector norms for vector jumps and retain individual components for diagnosis. A division guard such as `1e-300` is not an independently justified physical tolerance.

Evaluate residuals from the reconstructed primitive states independently of the nonlinear solver's internal scalar residual. Zero normal field and zero tangential flow must be legitimate test cases.

### 12.5 Front discretization

The analytical surface is the geometric authority. A visualization/diagnostic mesh approximates it. Its vertices, triangle orientation, area weights, and error bounds must be reproducible.

Do not interpolate compression across a support edge, across an accepted/sub-fast boundary, or across a failed jump. Interpolate primitive/history inputs consistently and reevaluate the local jump where needed. A thin high-Mach patch must not disappear solely because a coarse angular mesh missed it; quantify angular convergence and resolved coverage.

At a magnetic separatrix or current sheet, apply the provider's declared sidedness. Do not average opposite polarities into a fictitious null and then use that averaged field to infer an artificially high Alfvén Mach number.

### 12.6 Evolution is a sequence of diagnostic epochs

This mode does not admit or advance shocked material parcels. Each epoch reevaluates front geometry, ambient queries, and local jumps. A patch can transition from accepted shock to sub-fast and later back without resetting the event history.

Retain earlier exported records and generation lineage. A sub-fast transition does not erase previously written diagnostics, and it does not require the unimplemented sheath model to remain evaluable. No claim is made about the subsequent spatial location of that earlier shocked plasma.

Keep trajectory/model-validity flags separate from pointwise shock status and numerical failure. For example, `PROXY_OUTSIDE_DECLARED_VALIDITY` can coexist with a conservatively solved local jump, while `NON_FORWARD_INFLOW` describes that specific front/ambient combination. A required scientific-coverage failure has its own qualification result. A numerical candidate rejection leaves the previous committed epoch, clock, and generation unchanged; it does not silently advance with missing required patches.

## 13. Magnetic connectivity and observer passage

### 13.1 Ambient field-line tracing

At a specified epoch, trace an ambient field line using

\[
\frac{d\mathbf x}{ds}=\pm\frac{\mathbf B(\mathbf x,t)}{|\mathbf B(\mathbf x,t)|}.
\]

This is an instantaneous geometric line, not a plasma trajectory or a particle orbit. Trace both signs as needed to establish connectivity. HCI/Carrington transformations and any rotating/retarded field pattern must use the same epoch as the front.

Use adaptive spatial steps, local error control, and event detection for the photosphere, outer boundary, nulls, and domain exits. Field-line errors near separatrices can produce large connectivity changes even when local field errors are small. Report sensitivity rather than claiming an exact footpoint.

### 13.2 Intersections with a finite front

Locate all relevant roots of

\[
g(\mathbf x(s),t)=0
\]

that lie on the supported patch and above the physical inner boundary. A sign-change search alone misses tangent intersections; use a distance/minimum test or another validated tangent-detection method. Deduplicate coincident roots with a physical distance criterion.

A field line can intersect the surface more than once. Retain the full intersection list and classify each point by local shock status. If a consumer selects an intersection, its rule must be explicit—for example, the first accepted intersection reached from a specified observer along the relevant propagation direction. Do not assume the first root or the shock apex is always the accelerator.

Store observer identity, field-line identity, trace direction, intersection position, signed polarity information, arc length, local normal speed, `Mf`, `theta_Bn`, compression, and validity flags. Connectivity is based on the ambient approximation; it does not include eruption-generated reconnection or downstream draping.

An observer-connected accepted shock intersection can be called a **cobpoint**, following established shock/particle connection studies [R15]. Keep the full root list and the stated selection rule; this term does not imply that particles are being modeled.

Optionally publish `d(theta_Bn)/dt` and `d(Mf)/dt` along a continuously tracked selected intersection branch. For a scalar diagnostic `F` at its position `x_c(t)`,

\[
\frac{dF(\mathbf x_c(t),t)}{dt}
=\partial_tF+\dot{\mathbf x}_c\cdot\nabla F.
\]

Changing ambient time, front geometry, and moving field-line intersections all enter. A fixed-surface-label derivative is not the connected-point derivative. Preserve branch identifiers; at root mergers, grazing events, selection changes, null encounters, or unresolved discontinuities, publish an absent derivative and one-sided records where meaningful. Never difference across a branch change to manufacture a finite acceleration/obliquity rate.

### 13.3 Fixed or moving observers

For an observer trajectory `x_obs(t)`, a passage candidate satisfies

\[
g(\mathbf x_{\rm obs}(t),t)=0.
\]

Use bracketed temporal event location with support checks. At a tangent encounter, distinguish grazing contact from a transverse crossing. Observer passage speed involves relative motion:

\[
V_{\rm rel,n}=V_n-\dot{\mathbf x}_{\rm obs}\cdot\mathbf n.
\]

Report both local `Vn` and the observer-relative value. If the observer is outside finite support, report a miss. Arrival of the apex at 1 AU is not equivalent to shock arrival at Earth.

### 13.4 Meaning of future SEP applications

The reduced outputs can drive a separately declared phenomenological release source and upstream transport. They cannot determine an SEP spectrum, injection efficiency, or time-integrated acceleration from compression alone.

Explicit DSA requires a treatment of repeated crossings, upstream/downstream scattering and flow, and relevant finite-scale approximations. The local jump limits can support a local planar acceleration module, but that module and its downstream extent/boundary conditions must be qualified separately [R9]. Do not enable it implicitly while qualifying the background-only mode.

## 14. Input contract

### 14.1 Resolution rules

Resolve the event before mesh allocation or production stepping. Reject unknown keys, duplicate definitions, unsupported enum values, invalid units, inconsistent intervals, missing assets, and contradictory model authorities.

Relative asset paths resolve against the file containing the reference, not the shell's current directory. Include the transitive asset contents and resolved physical settings in the event fingerprint. A renamed file with unchanged bytes can retain a content identity, but its provenance path remains in the manifest.

Every input is classified as active, explicitly inactive for this mode, or unsupported. An old sheath/contact setting must not silently influence reduced-mode geometry. A parser can reject it or retain it as a clearly inactive legacy field according to a documented migration policy.

### 14.2 Run and ambient inputs

The names below describe the proposed schema; they are not an assertion of current parser support.

| Input | Meaning / unit | Constraints |
| --- | --- | --- |
| `mode` | `shock_front_ambient` | Explicit reduced-mode selection |
| `reference_epoch`, `time_scale` | Absolute epoch and clock convention | Required for observed maps or observers |
| `coordinate_frame` | HCI realization | One resolved canonical frame |
| `start_time_s`, `end_time_s` | Valid integration interval | End greater than start and within every required history |
| `background_dt_s` | Maximum epoch publication interval | Positive; convergence and event resolution assessed |
| `particle_mode` | `disabled` | No injection, seed population, or particle advance in qualification |
| `ambient_profile` / `plasma_asset` | Canonical analytical or empirical upstream model | One authority with full plasma/EOS specification |
| `gamma` | Shock adiabatic index | Greater than 1; matches the jump EOS |
| `composition` | Mass and pressure interpretation | Explicit; not inferred from energetic species |
| `density_normalization` | Density at a named point/surface | Positive, with units and normalization location |
| `temperature_model` | Isothermal reference, empirical profile, or qualified wind closure | Total electron/ion pressure and EOS checked; validity/uncertainty declared |
| `gamma_eff` | Effective exponent for a selected polytropic ambient wind | Distinct from shock `gamma`; inactive for other models |
| `wind_model` | Velocity profile and frame | Open/closed routing and physical role declared |
| `ambient_time_model` | Static, rigidly rotating, or another supported history | No undocumented map interpolation |

### 14.3 Magnetic inputs

| Input | Meaning | Constraints |
| --- | --- | --- |
| `magnetic_model` | PFSS plus validated exterior, or another maintained provider | Synthetic and observed modes distinguished |
| `magnetogram_asset` | Map plus metadata | Observation dates, frame, units, checksum |
| `harmonics_asset` | Precomputed coefficients as an alternative | Basis/normalization and map lineage required |
| `source_surface_radius_m` | `Rss` | Greater than solar radius; sensitivity recorded |
| `maximum_degree` | Harmonic resolution | Nonnegative supported degree; production has no monopole-only field |
| `monopole_policy` | Reject or explicitly remove | Removed flux recorded |
| `filter_model` / parameters | Map filtering | Applied once and checksummed |
| `amplitude_calibration` | Optional global constant field rescaling or qualified flux construction | No independent radial scalar factors; low-corona and IMF constraints recorded |
| `rotation_rate_rad_s` | Pattern angular speed | Signed convention and frame consistent |
| `exterior_field_model` | Flux/spiral matching construction | Includes its velocity assumptions and supported domain |
| `field_trace_tolerances` | Adaptive tracing controls | Physical units; separate null and numerical failure handling |

### 14.4 Front and history inputs

| Input | Meaning | Constraints |
| --- | --- | --- |
| `geometry` | Finite triaxial ellipsoid or finite SSE-inspired leading surface | Supported mathematical family |
| `direction_hci` | Propagation unit vector | Nonzero; normalized with a validation bound |
| `half_width_rad` | Finite SSE support | Regular supported interval; limits handled explicitly |
| `width_model` / `width_history_asset` | Fixed width or derivative-consistent time history | Fixed and history authorities mutually exclusive; full interval checked |
| `direction_history_asset` | Optional unit-direction history and derivative | Does not compete with a fixed direction; unit/tangent derivative constraints |
| `center_history` | `C(t)` and derivatives | Complete time coverage |
| `semiaxis_history` | Three `a_i(t)` and derivatives | Positive throughout interval |
| `attitude_history` | Proper rotation and derivatives | No instantaneous jumps |
| `support_definition` | Leading patch / finite footprint | Artificial closures excluded |
| `initial_apex_radius_m` | `R0` | Supports at least a declared coronal surface patch |
| `initial_apex_speed_m_s` | `V0` | Explicit reference frame |
| `final_pulse_speed_m_s` | `Vf` for analytical pulse | Consistent with selected history |
| `acceleration_duration_s` | `tau_a` | Positive |
| `history_asset` | Tabulated alternative | Does not compete with another active history |
| `history_quantity` | `shock_front` or explicitly named proxy | Avoid body/front confusion |
| `extrapolation_policy` | Default reject | Supported law and finite authorized interval |

### 14.5 Handoff and jump inputs

| Input | Meaning | Constraints |
| --- | --- | --- |
| `handoff_apex_radius_m` | `Rh` | Event crosses it inside valid history |
| `handoff_policy` | First outward crossing or another supported rule | Deterministic |
| `outer_propagation_model` | Existing SWCME kernel or declared reference | Represents complete matched shape |
| `drag_gamma_m_inv` | `Gamma` when quadratic drag selected | Nonnegative; not an EOS gamma |
| `drag_wind_policy` | Canonical wind or fitted effective speed | Identity and discrepancy explicit |
| `effective_wind_speed_m_s` | Constant fitted trajectory wind | Active only for the fitted-constant policy; never replaces Mach-number wind |
| `trajectory_validity_policy` | `report_only` or `require_declared_shock_coverage` | Separate proxy validity, numerical failure, and scientific coverage |
| `required_shock_coverage_asset` | Optional time/surface/observer requirements | Required for the coverage policy; no implicit requirement on every geometric edge |
| `transition_interval_s` | Optional smooth trajectory transition | No independent position/speed blending |
| `jump_model` | Maintained isotropic ideal-MHD fast solver | Physical family and EOS explicit |
| `jump_abs_tolerances` | Dimensional conservation bounds | Complete units and scales |
| `jump_rel_tolerances` | Relative conservation bounds | Declared before qualification |
| `weak_shock_policy` | Typed unresolved result | No compression or Mach floor |
| `upstream_sampling_model` | Direct ambient or one-sided grid stencil | Spatial refinement specified |

### 14.6 Numerical and output inputs

Declare angular/front resolution, spatial trace tolerance, event-time tolerance, maximum trace length/steps, surface output cadence, observer list and ephemerides, output formats, native mesh mode, and the required validation profile.

Scientific thresholds—such as a Mach number used for a diagnostic area statistic—must be labeled as analysis thresholds. They are not substituted for physical shock admissibility or used to hide weak-shock numerical failures.

The output directory is explicit and each run uses a unique directory. Restart is allowed only when event, ambient, geometry, EOS, and numerical-contract fingerprints match, except for explicitly allowed output-only changes.

Additional proposed controls include `require_fast_shock_at_observer`, `emit_trajectory_validity`, `magnetic_compression`, `ht_frame` and its conditioning/velocity bounds, `criticality_model` with its optional checksummed table, and connected-branch derivative/selection controls. An enabled diagnostic lacking a required validity bound or asset is a resolution error; a sample outside that diagnostic's valid domain has an explicit absent value.

If ensembles are enabled, declare parameter distributions, correlations, normalized member weights, reproducible sampling, common surface/observer labels, and unknown-outcome treatment. A parameter range without a probability measure is a sensitivity scan, not a probability model. All these settings and assets enter the event fingerprint.

## 15. Output contract

### 15.1 Resolved manifest

Write one machine-readable manifest before production stepping, containing:

- Model/document/schema versions and actual executable/source identity.
- Every resolved active parameter and explicitly inactive legacy parameter.
- Asset checksums, magnetic-map provenance, calibration, and time/coordinate transforms.
- Ambient and front provider identities and validity domains.
- EOS/composition/constants, selected jump family, and tolerance contract.
- History quantity (front or body proxy), trajectory-validity/required-coverage policies, diagnostic frame definitions, optional criticality model, and ensemble measure.
- Domain mode, geometry/trace resolutions, MPI rank count, and deterministic partition rules.
- Exact command, environment/build configuration identity, and selected test/run profile.
- Declared omissions: no downstream volume, contact, sheath, ejecta, particles, or feedback.

An absent observation is marked absent. A synthetic fixture is labeled synthetic. The manifest must not imply that an illustrative configuration has been calibrated against an event.

### 15.2 Per-epoch metadata

At each committed epoch, retain elapsed and absolute time, composite generation, front generation, ambient generation, handoff phase, complete front parameters and derivatives, support area, status counts, and checksums/identity references.

Also retain the trajectory/canonical wind comparison, apex fast-speed margin, validity flags, and any declared coverage result. A handoff receipt includes one-sided accelerations and the chosen `C1` or `C2` contract; derivative continuity cannot be inferred from two distant output samples.

An ambient field that is static can retain its ambient generation while the front generation advances. Do not claim CME-driven field updates merely because the front moved. An ambient time history can advance its own generation, provided the composite epoch commits both states atomically.

### 15.3 Surface records

Surface records have stable sample identifiers and explicit area weights. Required values include:

| Group | Required values |
| --- | --- |
| Geometry | Time, generation, position, outward unit normal, local `Vn`, apex speed, area weight, support/edge flags |
| Upstream | `rho1`, `p1`, `U1`, `B1`, EOS identity, upstream validity |
| Diagnostics | `w1`, `cs1`, `vA1`, `cf1`, fast-speed margin, `Mf1`, `Ms1`, `MA1`, `MAn1` where defined, beta, signed cosine, acute obliquity |
| Accepted jump | `X`, `rho2`, `p2`, `U2`, `B2`, derived temperature/entropy information, downstream validity |
| Additional diagnostics | `CB` and validity; HT boost/incident-flow vectors with frame, conditioning, and validity; optional compatible `Mc1` and criticality |
| Numerical evidence | Status, root/bracket information, conditioning, each RH residual and its scale |
| Provenance | Ambient identity, front identity, event fingerprint, local sample identifier |

Downstream values are absent for unsupported, sub-fast, or failed samples. In JSON, use null plus validity/status information rather than nonstandard NaN/Infinity. In CSV, use explicit validity columns and documented absent-value encoding. Never serialize absent states as zeros that look like physical vacuum.

### 15.4 Surface-integrated diagnostic quantities

Use actual Cartesian area weights and convergence-tested quadrature. Examples are:

\[
A_{\rm shock}=\int_{\mathcal S_{\rm accepted}}dA,
\qquad
\langle X\rangle_A=
\frac{\int_{\mathcal S_{\rm accepted}}X\,dA}{A_{\rm shock}},
\]

\[
\dot M_{\rm swept,diagnostic}=
\int_{\mathcal S_{\rm accepted}}\rho_1w_1\,dA.
\]

The latter is local incoming mass flux integrated over the modeled shock surface. Its time integral is **not a predicted retained sheath mass**, because downstream storage, exits, and material transport are not modeled. Integrals exclude unsupported/failed samples, and report the excluded area rather than renormalizing silently.

When `A_shock=0`, average compression is undefined. Report a valid zero shock area and an absent average; do not divide by zero or assign a fake mean compression of 1.

Distinguish the super-fast candidate area from the accepted RH shock area. Report failed/unknown area separately, including its fraction of total supported area. For an axisymmetric fixture the shock extent may be a half-angle; for a structured ambient, export azimuth-dependent intervals or the actual patch boundary, including disconnected regions. A single maximum half-angle does not establish filled shock coverage.

### 15.5 Connectivity, observer, and native receipts

Write separate files for field-line traces, intersection records, observer passage events, and native storage evidence. Preserve multiple epochs; overwriting the final epoch cannot demonstrate handoff continuity.

Connected records retain branch identity, selection rule, `CB`, HT diagnostic validity, optional criticality, and any branch-consistent time derivatives. Observer records distinguish geometric hit, accepted shock hit, sub-fast hit, unsupported coverage, and numerical unknown. Ensemble outputs retain member status/weights and the measure used for every probability.

Native receipts include owned-cell counts, ghost counts, ambient-byte/value errors, generation agreement, particle count, and collective commit/rejection results. If a front is not stored in every cell, record how its shared analytical epoch or partitioned diagnostic samples are validated.

### 15.6 Visualization

Preferred outputs are CSV/JSON for exact diagnostics and a maintained surface/field format such as VTK or the project's Tecplot convention for visualization. The surface geometry and field-line coordinates must be in the same frame and units.

Useful views include:

- Front colored by compression, Mach number, obliquity, and solver status.
- Ambient field lines colored by polarity or topology, with marked front intersections.
- Cuts of ambient field magnitude/direction and plasma beta.
- Time series along fixed surface labels and observer-connected intersections.
- Pre/inside/post-handoff geometry and diagnostic comparisons.

Plot labels must say **ambient reference** for volume fields. No plot may label the unmodeled downstream interior as a solved CME magnetic field or sheath.

## 16. Epoch preparation and native execution

### 16.1 Candidate algorithm

The following is semantic pseudocode, not an existing function-name contract:

```text
resolve_event_and_assets()
verify_mode_capabilities_and_validity_intervals()
prepare_ambient_reference()
prepare_continuous_front_history()
resolve_handoff_event_and_matched_state()

for each requested epoch t:
    candidate = private immutable epoch
    candidate.ambient = ambient_reference_at(t)
    candidate.front = complete_front_geometry_at(t)
    evaluate trajectory validity, canonical/effective wind, and apex margin

    for each deterministic supported surface sample:
        evaluate Cartesian position, normal, normal speed, and area weight
        query canonical upstream ambient at the same epoch
        classify inflow and fast-mode admissibility
        solve and independently check the local jump when admissible
        evaluate enabled magnetic compression, HT, and criticality diagnostics
        retain explicit status and validity for every sample

    calculate requested connectivity and observer diagnostics
    track selected intersection branches and valid requested time derivatives
    validate candidate geometry, field, jump records, and required coverage
    prepare native ambient-reference storage and generation metadata
    collectively accept or reject the complete candidate
    if accepted:
        publish storage, halo, and epoch atomically
        append receipts and scientific outputs
    else:
        retain the previous committed epoch unchanged
        report a structured error; do not advance the production clock
```

Valid sub-fast classifications do not reject an otherwise supported epoch. Invalid upstream states, unresolved required jumps, corrupted assets, inconsistent histories, and unacceptable numerical coverage do reject it under the declared contract. Do not convert these failures to sub-fast success to keep the run moving.

Apply `report_only` versus required scientific shock coverage explicitly. A numerical fixture can pass its declared geometry/diagnostic assertions while recording sub-fast or non-forward samples. A production shock-arrival gate cannot pass by applying that fixture's weaker success criterion.

### 16.2 Transaction and MPI rules

Preparation is private. It must not mutate the installed ambient/front authority or native buffers before collective acceptance. Rank-local rejection is communicated before any rank enters a different collective path. A failed preparation leaves the previous state, generations, and clock coherent.

Use stable global surface/sample identifiers and deterministic partitions. MPI rank identity is not a scientific sample label. Aggregate area integrals with a declared reproducible/tolerance-bounded reduction; distinguish scientific rank invariance from bitwise summation identity.

The event and matched crossing metadata are resolved consistently, for example on a designated authority and broadcast with checksums. Different ranks must not independently select different crossing times or branch roots because of a race or inconsistent configuration.

### 16.3 Native background storage

Populate native storage through the maintained SEP3D runtime-provider adapter. Disable scheduled-file updates/interpolation for that provider using the existing supported mechanism. Do not invent an application-specific core coupling mode or modify generic PIC merely to add a new provider name.

Owner and received ghost buffers must represent the same committed **ambient reference** and generation. Front masks or diagnostics, if stored on the mesh, have their own declared layout and epoch identity. A ghost check compares actual received data with independently evaluated expected values and cannot pass by rereading the sender's preparation buffer.

All allocated mesh regions required by the host must receive valid storage. Writing ambient-reference values into cells geometrically behind the front is allowed only with explicit reference-role metadata/masking. A physical downstream query remains unsupported there.

### 16.4 No particles

The qualification decks disable all particle sources, preloaded seeds, injection, and particle advance. Verify the actual global particle count remains zero. Merely setting an injection rate to zero is insufficient if a restart or another source creates particles.

Host stepping can still publish front/ambient epochs, halo data, diagnostics, and output. It must not require particle-only source spectra or wave state to instantiate this provider. Preserve baseline particle implementation for later modes rather than deleting it from the repository.

### 16.5 Background-only restart

A reduced checkpoint stores event/asset fingerprints, committed time and generation, trajectory phase/matched state, ambient identity, and any numerical outer-trajectory state needed for deterministic continuation. Analytical histories may be reevaluated, but their resolved parameters and crossing identity must match.

Restart does not require sheath inventories or particle checkpoints. Its equivalence test compares later front geometry, ambient values, jump classifications, observer/connectivity records, and receipts with an uninterrupted run. Changed physical assets or EOS are rejected rather than treated as output-only changes.

## 17. Domain, timestep, and sampling choices

### 17.1 Full domain versus corridor

Domain selection belongs to the application input and changes computational coverage, not the underlying analytical front.

| Domain mode | Appropriate use | Limitation |
| --- | --- | --- |
| Full coronal/heliospheric mesh | Broad front/field visualization and global native coverage | Greater memory and halo cost |
| Refined Parker/open-field corridor | Native checks or selected observer-connected regions | Omits other flanks, closed loops, and off-corridor structures |
| Provider-only analytical diagnostics | Complete front sampling and field-line tracing without a global host mesh | Does not qualify native AMPS owner/ghost integration |

A corridor deck must identify the observer, field-line family, angular/depth margins, and outer coverage. It cannot claim full-front native coverage. A field-line tracer can use the analytical provider outside the allocated mesh only if this is explicitly supported and recorded.

The model's geometry, normals, and jump parameters are not derived from AMR cell boundaries. Refinement improves storage/visualization/interpolation tests, but analytical surface accuracy has its own resolution controls.

### 17.2 Timestep meaning

There is no plasma-evolution CFL condition because this mode does not evolve a volumetric fluid. Nevertheless, epoch cadence must resolve front acceleration, passage events, rapid changes in ambient conditions sampled by the front, and connectivity changes.

A useful cadence criterion is

\[
\Delta t\,|V_n|\ll L_{\rm ambient},\qquad
\Delta t\ll\tau_{\rm acceleration},
\]

with physically defined local variation scales and a tested error target. Treat small/zero velocity and field gradients safely. Event roots can be located between output epochs if the histories/provider are continuously queryable; otherwise coarse output cadence limits timing precision.

Separate:

- Front/ambient publication timestep.
- Numerical outer-trajectory integration timestep.
- Geometry/field-line spatial resolution.
- Output cadence.

Converge each independently. A 10- or 20-step smoke run is not a propagation-to-1-AU qualification unless its actual committed epochs cover that evolution.

### 17.3 Endpoint and observer coverage

Find the apex endpoint root `R_a(t)=1 AU` or the configured greater radius. Export its actual resolved time and state. Also determine supported observer passages independently.

The endpoint requirement is satisfied by actual valid epoch evaluation and required native evidence at that radius, not by a history that merely permits extrapolation there. A short handoff deck and a long 1-AU deck must have different names and manifests.

## 18. Commented configuration examples

**Schema status:** the examples in this section are proposed logical configuration decks. They illustrate the required inputs and comments. They are not drop-in files for the current strict SEP3D parser. Implementation must translate these settings into the maintained schema, add parser tests, and document the accepted spelling before distributing runnable decks. Referenced ambient/harmonic assets must also be supplied and validated.

### 18.1 Short handoff smoke case

This synthetic case begins just below the front handoff radius. It checks phase matching and native generation/halo behavior without waiting for a full coronal transit. It does not claim to simulate CME formation.

```ini
# PROPOSED logical schema: translate to the maintained parser before use.
# This is a synthetic zero-particle handoff fixture, not an observed event.

[run]
mode = shock_front_ambient
reference_epoch = 2026-10-04T00:00:00Z
time_scale = UTC_resolved_to_elapsed_SI_seconds
coordinate_frame = HCI
start_time_s = 0
end_time_s = 600
background_dt_s = 60
particle_mode = disabled

[ambient]
# One canonical plasma authority; the asset includes profile equations,
# normalization, open/closed routing, EOS, validity range, and provenance.
provider = canonical_corona_wind
plasma_asset = ../assets/synthetic_ambient.profile.json

[magnetic]
# SYNTHETIC field for verification. Replace with observed-map harmonics
# plus provenance and uncertainty for event-specific qualification.
model = pfss_with_validated_exterior
harmonics_asset = ../assets/synthetic_dipole.harmonics.csv
source_surface_radius_m = 1739250000
monopole_policy = reject
rotation_rate_rad_s = 2.86533e-6
exterior_field_model = canonical_flux_matched_spiral

[shock]
model = maintained_isotropic_fast_RH
gamma = 1.6666666666666667
weak_shock_policy = typed_unresolved_failure
upstream_sampling = direct_canonical_ambient

[geometry]
model = finite_sse_leading_surface
# 45-degree half-width in radians; this is an axisymmetric shape fixture.
width_model = fixed
half_width_rad = 0.7853981633974483
direction_hci = 1, 0, 0
support = leading_ray_root_only
physical_inner_radius_m = 695700000

[coronal_history]
model = constant_apex_speed
# 19.9 R_sun: t_h = (20 - 19.9) R_sun / 1000 km/s = 69.57 s.
initial_apex_radius_m = 13844430000
initial_apex_speed_m_s = 1000000
valid_start_s = 0
valid_end_s = 600
quantity = shock_front

[handoff]
apex_radius_m = 13914000000
crossing = first_outward
outer_model = quadratic_drag_reference_or_equivalent_SWCME_kernel
drag_gamma_m_inv = 2e-11
# This is an explicitly fitted effective trajectory wind in this fixture.
# It does not replace the canonical wind used for local Mach numbers.
drag_wind_policy = fitted_effective_constant
effective_wind_speed_m_s = 400000
continuity = position_and_normal_speed
# A body-drag shock proxy can become sub-fast; export that outcome.
# For a validated event, select require_declared_shock_coverage and provide
# its checksummed time/surface/observer requirements instead.
trajectory_validity_policy = report_only

[diagnostics]
magnetic_compression = true
ht_frame = false
# Enabling HT output also requires declared conditioning/velocity bounds.
# No general first-critical Mach number is inferred from a constant.
criticality_model = none

[domain]
# Full domain: native coverage over the declared mesh.
# Corridor: change mode and supply observer/field-family footprint settings;
# then label native coverage as corridor-only. Geometry remains analytical.
mode = full_domain

[output]
directory = test_output/reduced_front/handoff_smoke
retain_all_epochs = true
surface_formats = csv, vtk
receipt_format = json
volume_field_role = ambient_reference
emit_trajectory_validity = true
```

The minimum useful receipt set includes initialization, at least one pre-handoff epoch, the resolved crossing, and a post-handoff epoch. With output at 0, 60, and 120 s, the 60-s record is before handoff and the 120-s record is after it; the crossing record should be evaluated explicitly at 69.57 s. Acceleration continuity is not promised for a sharp switch to the drag law.

### 18.2 Low-corona-to-1-AU case

Use a separately named deck, for example `shock_front_corona_to_1au.in`, with the same explicit ambient/magnetic authorities and a history such as:

```ini
# PROPOSED replacements for the corresponding smoke-deck sections.
# Do not append a second copy of an existing section/key.

[run]
mode = shock_front_ambient
reference_epoch = 2026-10-04T00:00:00Z
time_scale = UTC_resolved_to_elapsed_SI_seconds
coordinate_frame = HCI
start_time_s = 0
end_time_s = 345600
# Coarse illustrative maximum: use a smaller validated cadence during
# the acceleration/low-coronal phase and demonstrate convergence.
background_dt_s = 60
particle_mode = disabled

[coronal_history]
model = quintic_velocity_pulse_then_constant_speed
initial_apex_radius_m = 800055000
initial_apex_speed_m_s = 100000
final_pulse_speed_m_s = 1000000
acceleration_duration_s = 600
valid_start_s = 0
valid_end_s = 345600
quantity = shock_front

[endpoint]
apex_radius_m = 149597870700
evaluate_exact_crossing = true
# A synthetic +X target is not the actual Earth ephemeris.
observer_id = synthetic_target_x
observer_position_hci_m = 149597870700, 0, 0
# This fixture tests geometry/epochs even if the proxy loses shock validity.
# An event shock-arrival qualification must enable and pass this requirement.
require_fast_shock_at_observer = false

[output]
directory = test_output/reduced_front/corona_to_1au
retain_all_epochs = true
surface_formats = csv, vtk
receipt_format = json
volume_field_role = ambient_reference
emit_trajectory_validity = true
```

Here `R0=1.15 R_sun`. After the pulse, `R_a(600 s)=1.130055e9 m` and the speed is `1.0e6 m/s`. For a sharp handoff at 20 solar radii and the constant-drag reference above, the handoff occurs at about 13383.945 s. The reference apex reaches 1 AU at about 203884.494 s, or 56.635 hours after initialization. These are analytical fixture calculations, not an event forecast or a verified AMPS run. A four-day interval therefore covers the reference endpoint; an implementation must still locate and publish the actual endpoint and demonstrate its supported surface coverage.

If the maintained SWCME kernel uses another propagation law or the drag parameters vary, these reference times change. The full run's endpoint gate uses the actual resolved model, not these illustrative numbers.

This deck qualifies the specified geometry/epoch calculation, not shock survival to 1 AU. Retain the reported sub-fast/non-forward outcomes instead of tuning the ambient to make the test pass. A separately calibrated event deck must supply its own observed or assessed plasma/magnetic inputs, trajectory validity interval, and required shock coverage. Do not relabel the synthetic case as a prediction.

### 18.3 Synthetic harmonic asset example

For an orthonormal real harmonic basis, a pure dipole can be represented schematically as:

```csv
# PROPOSED asset format; the actual importer must define comments and columns.
# basis=real_orthonormal; phase=Condon_Shortley; coefficient_role=photospheric_Br
# units=T; reference_frame=declared_solar_axis_frame; monopole=absent
degree,order,cosine_T,sine_T
1,0,0.0002046653415892977,0
```

This coefficient corresponds to a polar photospheric amplitude of `1e-4 T` under the basis definition in section 6.2. It is an analytic verification input. A realistic event input must use reconstructed coefficients from a suitable magnetic map and retain that map's metadata; relabeling the dipole as an observed event is unacceptable.

### 18.4 Required asset structure

A runnable packet contains the translated input deck, magnetic map or harmonics, plasma profile asset, observer history where applicable, and a resolved manifest. The plasma asset must define actual profile equations/data and normalizations; a filename alone is not a physical specification.

Before distributing a packet, validate it from at least two working directories to prove that relative paths resolve correctly. Include one malformed/unknown-key deck and one missing-asset case in the parser tests. No current test runner should continue referring to renamed example files without a regression test for the canonical path.

### 18.5 Reproducible polar ambient-sensitivity reference

This additional probe uses the same pulse/drag apex history as section 18.2, but **overrides propagation and the target to +Z along the synthetic dipole axis**. It is not the +X deck's field orientation. Set `Rss=2.5 R_sun`, pure hydrogen `np=ne`, `Tp=Te=1.4e6 K`, shock `gamma=5/3`, and `np(1 AU)=7e6 m^-3`. Use the transonic isothermal Parker solution from section 5.3 with `ciso^2=2 kB T/mp`, and conserve its radial mass flux:

\[
n_p(r)=n_{p,\rm AU}\,
\frac{U_{\rm AU}r_{\rm AU}^2}{U(r)r^2}.
\]

For polar photospheric amplitude `Bp`, `z=r/R_sun`, and `q=2.5`, the PFSS/exterior polar field is

\[
B_{\rm pole}(r)=
\begin{cases}
B_p(1+2q^3/z^3)/(1+2q^3),&1\le z\le q,\\
B_p\,3(q/z)^2/(1+2q^3),&z\ge q.
\end{cases}
\]

At this apex, `B` and `n` are parallel, so use `cf=max(cs,vA)`. For reproducibility use `GM_sun=1.32712440018e20 m^3/s^2`, `kB=1.380649e-23 J/K`, `mp=1.67262192369e-27 kg`, `mu0=1.25663706212e-6 H/m`, and section 3's solar radius/AU. These give `ciso=152.027 km/s`, `cs=196.266 km/s`, `rc=4.12683 R_sun`, and `U_AU=601.223 km/s`. With `Bp=1e-4 T` (1 G) or `1e-3 T` (10 G):

| Radius / `R_sun` | Apex speed km/s | Wind km/s | `vA`, 1 G km/s | `Mf`, 1 G | `Mf`, 10 G |
| --- | ---: | ---: | ---: | ---: | ---: |
| 2 | 1000.000 | 49.334 | 334.155 | 2.845 | 0.284 |
| 4 | 1000.000 | 147.282 | 275.812 | 3.092 | 0.309 |
| 10 | 1000.000 | 281.187 | 152.439 | 3.662 | 0.472 |
| 20 | 1000.000 | 369.771 | 87.404 | 3.211 | 0.721 |
| 50 | 872.850 | 470.216 | 39.425 | 2.051 | 1.021 |
| 100 | 734.549 | 536.128 | 21.049 | 1.011 | 0.943 |
| 215.03216 (1 AU) | 582.592 | 601.223 | 10.366 | Absent: non-forward | Absent: non-forward |

The signed ratio `(Va-U)/cf` at 1 AU is about -0.095 for both amplitudes. It is not a positive inflow Mach number; the correct status is `NON_FORWARD_INFLOW`. The fitted-wind proxy is overtaken by this canonical polar wind at about `199.9234 R_sun`. No physical speed or Mach floor is permitted.

The 10 G example is not sub-fast everywhere: the 50-solar-radius sample is slightly super-fast. It remains a candidate until the parallel/switch-on branch requirements and conservative RH residuals are checked. These synthetic amplitudes do not certify a realistic coronal event.

This recomputation corrects the review's polar-table orientation: its quoted values follow `sqrt(cs^2+vA^2)`, appropriate for a perpendicular field. Treat using that expression for this parallel probe as a negative control. The probe also demonstrates why a constant 400 km/s trajectory wind cannot silently stand in for the actual 601 km/s upstream wind at 1 AU. A different calibrated ambient could give different outcomes; this is a sensitivity reference, not a universal shock lifetime.

### 18.6 Width-history reference

For a separate bounded geometry test, prescribe fixed `d` and

\[
\lambda(t)=\lambda_0+(\lambda_1-\lambda_0)h(t/\tau_\lambda),
\quad 0\le t\le\tau_\lambda,
\]

using section 9.2's quintic `h`, `lambda0=pi/6`, `lambda1=pi/4`, and `tau_lambda=300 s`. Keep the endpoint widths constant outside the pulse within the declared history interval. The derivative is `dot(lambda)=(lambda1-lambda0)30 s^2(1-s)^2/tau_lambda`, with `s=t/tau_lambda`. It vanishes at the endpoints.

The proposed resolution selects `width_model=history` and a checksummed derivative-consistent asset, instead of the fixed-width authority. Evaluate the complete section 8.3 normal speed against independent Cartesian time differentiation, including changing tangency/support. Add a separate direction-rotation reference to isolate the `c dot(d) dot n` term. This fixture verifies the input geometry; it does not claim to compute pressure-driven lateral expansion.

## 19. Validation specification

### 19.1 Registry and evidence policy

The `RSHxx` identifiers below are proposed labels. Reconcile them with the maintained registry before implementation; do not rename existing tests or reuse an identifier with another meaning. A single aggregate PASS must retain the result of each required assertion.

Each test records command, actual source/executable identity, input fingerprint, resolution, tolerances, expected values/reference provenance, measured errors, duration, and PASS/FAIL/SKIP/ERROR. A missing prerequisite is SKIP or ERROR according to the registry policy, never PASS. Any skipped required gate keeps the relevant qualification open.

Numerical tolerances are declared before assessing a fixture and include dimensional scales. Event uncertainty is separate from arithmetic or discretization error. A larger observational uncertainty does not excuse a violation of exact local conservation equations.

### 19.2 Required tests

| Proposed ID | Case | Required evidence |
| --- | --- | --- |
| RSH01 | Input resolution and units | Exact SI conversion, checksum coverage, duplicate/unknown-key rejection, path resolution independent of cwd |
| RSH02 | Coordinate/time transforms | Known HCI/map-frame rotations, velocity transformation, time ordering, reference-epoch consistency |
| RSH03 | Translating plane | Analytic normal and normal speed; radial/apex substitutions fail negative control |
| RSH04 | Expanding Sun-centered sphere | `n=x/r`, `Vn=dR/dt`; area `4 pi R^2` under quadrature refinement |
| RSH05 | Translating/expanding ellipsoid | Complete analytic level-set derivatives and distance residual match independent differentiation |
| RSH06 | Rotating triaxial ellipsoid | Angular-velocity contribution present; omitting `dot(R)` is detected |
| RSH07 | Finite SSE geometry | Apex, far-ray root, tangent flank, support miss, width limits, normal/radial speed distinction, fixed-width zero-speed edge, uniform-reference shock half-angles |
| RSH08 | Acceleration history | `dR/dt=V`, `dV/dt=a`, endpoint values, positivity and between-knot behavior |
| RSH09 | Outer trajectory | Closed-form drag cases above/below/equal wind, zero drag, numerical trajectory convergence |
| RSH10 | PFSS harmonic basis | Independent harmonic references, normalization, source-surface radial condition, pole limits [R3] |
| RSH11 | Magnetic-map processing | Signed flux, raw/calibrated identities, reconstruction error, monopole policy, filtering applied once |
| RSH12 | Magnetic field/derivatives | Divergence and valid curl conditions; actual interpolation checked; nonaxisymmetric exterior negative control |
| RSH13 | Magnetic matching | Normal flux, field continuity contract, polarity, line mapping; distinguish gradient jumps from failures |
| RSH14 | Ambient plasma/EOS | Positive physical state, total electron/ion pressure, mass normalization, sound/Alfvén/fast speed including parallel/perpendicular limits; hydrostatic or wind references where claimed |
| RSH15 | Frame/sign invariance | Galilean change preserves `w1` and jump invariants; normal/side conventions handled consistently |
| RSH16 | Gas-dynamic limit | `Ms=3`, `gamma=5/3`: `X=3`, `p2/p1=11`; temperature and downstream normal speed exact |
| RSH17 | Perpendicular finite-beta jump | `beta1=0.2`, `MA=2`: `X=2`, field ratio 2, pressure ratio 6; full energy and entropy checks |
| RSH18 | Oblique/parallel jump family | Independent states/residuals and verified special branches; full-vector multiplication by X fails negative control |
| RSH19 | Weak and sub-fast cases | No artificial floor, correct absent/downstream validity, unresolved positive Mach excess is typed failure |
| RSH20 | Handoff | Complete pre/crossing/post geometry, no authority reset; normal speed/jumps continuous as inputs allow; one-sided acceleration receipt and separate optional C2 transition |
| RSH21 | Field-line/front intersections | Independent simple field references, all roots, tangent root, duplicates, finite-support rejection, polarity |
| RSH22 | Observer crossing | Known apex/flank/miss/graze cases, moving observer relative speed, event-time convergence |
| RSH23 | Surface integrals | Actual area weights; geometric, super-fast, accepted and unknown area distinguished; asymmetric/disconnected footprints; zero-area handling and thin-patch convergence |
| RSH24 | Epoch transactions | Injected preparation/RH/asset/rank-local failure leaves committed state and clock unchanged |
| RSH25 | Native owner and halo | Actual native owner readback and received ghost values/metadata match the ambient-reference epoch |
| RSH26 | Native rank invariance | One- and four-rank epoch, jump, integral, and connectivity agreement; stable sample identity |
| RSH27 | Particle exclusion | Actual zero particles throughout initialization and stepping; legacy/source paths cannot create seeds |
| RSH28 | Actual 1-AU endpoint | Valid committed endpoint, accepted/absent status coverage, observer results, field values, and native receipts |
| RSH29 | Background-only restart | Uninterrupted versus resumed later epochs; changed physical input rejected |
| RSH30 | Maintained regressions | Shared SWCME/corona, SEP3D, architecture and sanitizer checks; generic PIC/particle baseline behavior preserved |
| RSH31 | Ambient sensitivity | Reproduce the polar probe with parallel fast speed; temperature/density/field/wind variations retain EOS and normalization; report rather than floor sub-fast/non-forward states |
| RSH32 | Width/direction histories | Independent derivatives, tangent limits, interval constraints, complete matched derivatives; omitted width and direction terms fail separately |
| RSH36 | Additional diagnostics | Magnetic compression, HT vectors/frame invariance and electric-field cancellation; parallel/perpendicular/null conditioning; criticality reference/table tests when that option is enabled |
| RSH37 | Flux-consistent calibration | Global scaling preserves divergence/topology; independent radial scaling fails; flux continuity, raw/fitted lineage, and joint coronal/IMF constraints |
| RSH38 | Trajectory validity and coverage | Constant-body-drag loss of fast/forward inflow is reported; canonical/effective wind discrepancy recorded; required event coverage can fail without altering speeds or ambient |
| RSH39 | Connected-branch derivatives | Moving intersection chain rule, stable IDs, root mergers/grazes/selection breaks give absent or one-sided derivatives; fixed-label derivative fails negative control |
| RSH40 | Ensemble accounting | Declared/correlated measure, common labels, known no-shock versus numerical unknown, unconditional supported-shock coverage and unknown-weight bounds; no discarded failing members |

RSH33–RSH35 are reserved below for optional dynamics/formation extensions. Their absence does not block reduced S0–S7. Conversely, an extension cannot be labeled qualified merely because the required reduced tests pass.

| Deferred ID | Extension-specific gate |
| --- | --- |
| RSH33 | Driver/standoff closure: independent planar/nose references, evolving curvature, implicit speed/Mach coupling, validity loss near weak shocks, event comparison |
| RSH34 | General front dynamics/GSD: amplitude–area/angle closure, verified dimensional extension, characteristic propagation, topology/self-intersection handling and mesh convergence |
| RSH35 | Nonlinear formation: declared piston/wave initialization, amplitude-dependent characteristics, steepening reference, attenuation/geometry limits and independent onset comparison |

### 19.3 Independent references and negative controls

Reference tests exercise actual production kernels where applicable. Separate code that computes only a reference identity does not establish that the production query follows it.

Each important negative control has a designated gate:

| Deliberate error | Gate detecting it |
| --- | --- |
| Use apex speed instead of flank normal speed | RSH03, RSH07 |
| Omit attitude, width, or direction derivatives | RSH06, RSH32 |
| Choose the rear SSE root or close the cap with an artificial side shock | RSH07 |
| Swap colatitude/latitude or harmonic normalization | RSH02, RSH10 |
| Use perpendicular fast speed for the parallel polar probe | RSH14, RSH31 |
| Omit electron pressure in the common electron/ion-temperature EOS | RSH14 |
| Multiply the entire upstream field by density compression | RSH18, RSH36 |
| Average opposite polarities through an unresolved current sheet | RSH12, RSH19 |
| Apply independent radial magnetic amplitude factors | RSH13, RSH37 |
| Treat a numerical RH failure as sub-fast or silently exclude its area | RSH19, RSH23, RSH24 |
| Clip trajectory speed/Mach or substitute fitted wind in upstream diagnostics | RSH31, RSH38 |
| Blend speed without the position-weight derivative | RSH20 |
| Reinitialize outer propagation at an unrelated radius/time | RSH20 |
| Call the incident HT flow speed the tangential boost; ignore singular conditioning | RSH36 |
| Return a universal critical Mach value outside a model/table's domain | RSH36 when enabled |
| Difference across a connected-branch change or omit intersection motion | RSH39 |
| Drop numerical failures from ensemble weights | RSH40 |
| Reuse the prepared owner buffer as a received-ghost check | RSH25 |
| Publish an interior ambient reference as actual downstream CME plasma | RSH01, RSH30 |

Each must fail its intended assertion. A negative control that also alters an unrelated part of the test is not a useful causal check.

### 19.4 Refinement campaign

Refine independently:

1. Surface angular/chart resolution and quadrature.
2. Ambient interpolation/harmonic resolution, when applicable.
3. Field-line integration accuracy and root-location tolerance.
4. Trajectory integration step and history interpolation resolution.
5. Publication cadence and observer/handoff event sampling.
6. Native AMR storage resolution and MPI decomposition.

Use at least three levels where an asymptotic convergence claim is made. Compare norms, extrema, surface-integrated quantities, and classification/coverage changes. A near-threshold patch can change its binary status under refinement; report the affected area and uncertainty rather than demanding bitwise status invariance at every marginal point.

Numerical convergence need not be strictly monotonic at the roundoff limit, but an unresolved growing error cannot be advertised as convergence. A physics approximation that persists under refinement remains part of the model discrepancy.

### 19.5 Handoff smoke and long native campaign

Maintain separate translated decks:

- `examples/shock-front/handoff_smoke.in`: short initialization and pre/crossing/post-handoff receipts.
- `examples/shock-front/corona_to_1au.in`: prescribed low-coronal acceleration, complete outer propagation, and exact endpoint/observer receipts.
- `examples/shock-front/event_case.in`: optional observed map/history campaign, independently labeled and validated.

The following commands are **templates for after the reduced mode and these decks are implemented**. They use the already familiar native AMPS test interface; they are not evidence that the current binary supports the proposed mode.

```bash
mkdir -p test_output/reduced-front/handoff-1 test_output/reduced-front/handoff-4
mpiexec -n 1 ./amps --test-suite sep-corona --test-input srcSEP3D/examples/shock-front/handoff_smoke.in --test-steps 10 --expect-mpi-ranks 1 --test-json test_output/reduced-front/handoff-1/native.json --artifact-directory test_output/reduced-front/handoff-1/artifacts
mpiexec -n 4 ./amps --test-suite sep-corona --test-input srcSEP3D/examples/shock-front/handoff_smoke.in --test-steps 10 --expect-mpi-ranks 4 --test-json test_output/reduced-front/handoff-4/native.json --artifact-directory test_output/reduced-front/handoff-4/artifacts
```

For the long case, replace the input/output paths and select enough steps, using the **actual application timestep and termination semantics**, to reach the endpoint. Do not infer the number of native steps from `background_dt_s` if the host uses another step interval. The executable or harness must explicitly assert endpoint coverage. Use a valid allocation and the site's maintained MPI invocation.

Before **each native rebuild**, confirm the AMPS root and ensure no active build/test uses `build`, then remove only the root build directory with `rm -rf -- build`. Regenerate the selected `srcSEP3D` configuration and production hooks using the established site recipe, and compile with `-j16`. Do not delete source, site configuration, libraries, or the read-only Mars1 reference. Do not add `apply_fix.py`; modify maintained sources directly. No commit or push is part of this implementation instruction.

Portable tests, shared regressions, sanitizers, and native MPI evidence remain distinct. If this environment lacks MPI or a complete AMPS build, record the missing native evidence and let the user's system run it. Do not substitute a portable mock for a native halo receipt.

## 20. Uncertainty, calibration, and event validation

### 20.1 Major uncertainty sources

- Photospheric map age, missing regions, radial-field conversion, polar treatment, and instrumental calibration.
- PFSS source-surface radius, harmonic truncation/filtering, field-strength normalization, and open/closed topology.
- Ambient density, pressure/temperature, composition, and wind speed.
- Front shape, direction, extent, attitude, time interpolation, and observed-front identification.
- Shock-front versus CME-body distinction in an outer propagation law.
- Earlier CMEs or coronal restructuring absent from the ambient model.
- Observer ephemeris and field-line routing near separatrices.

Near a magnetic null or a fast-mode threshold, nonlinear sensitivity can dominate. Do not propagate only a small linear error bar when input perturbations change topology or shock classification.

### 20.2 Ensemble propagation

When uncertainties are supplied, resolve a reproducible ensemble of compatible inputs. Reevaluate front geometry, ambient states, jumps, and connectivity for each member. Preserve correlated input variations instead of perturbing all parameters independently by default.

Report distributions of arrival time, supported hit/miss fraction, shock area, connected shock parameters, and numerical/model rejection counts. A member that has no admissible shock is a scientific outcome; a member that fails numerically is a separate outcome. Do not remove either silently before calculating confidence intervals.

Probabilities require a declared measure, not just a list of parameter choices. With normalized weights `wi` and a common observer, ray, or specified surface-label mapping, one possible quantity is

\[
P(\mathrm{supported\ and\ superfast})
=\sum_i w_i\,\mathbf 1[
\mathrm{supported}_i\ \land\ w_{1,i}>0\ \land\ M_{f,i}>1].
\]

Only evaluable members contribute a known yes/no value. A valid support miss is a known negative for this joint event; a failed ambient/jump evaluation is not. Keep super-fast probability distinct from probability of an accepted RH shock. If known positive weight is `Wplus` and unresolved weight is `Wunknown`, report the bounds `Wplus <= P <= Wplus+Wunknown` and the unresolved fraction, rather than renormalizing the survivors. Invalid priors/assets are an ensemble setup error; this uncertainty ledger is not a way to accept malformed inputs.

Label any conditional probability, its conditioning set, and denominator explicitly. Report member-level deterministic statuses as the primary evidence. An average compression across members cannot replace the no-shock/accepted/unknown distribution. Arrival-time summaries are conditional on the stated kind of hit and must retain the complementary miss/unknown weight.

### 20.3 Observational checks

Useful independent comparisons include:

| Observable | Constrained model quantity | Interpretation limit |
| --- | --- | --- |
| EUV/coronagraph front contours | Shape, direction, extent, early kinematics | A bright front alone does not uniquely identify a shock |
| Type-II radio signatures | Evidence of a shock and density-related constraints | Radio source location and harmonic interpretation are additional assumptions |
| In-situ shock arrival and upstream/downstream jumps | Passage time, local normal/speed, compression, field jump | Only supported local passages are predicted; not the following ejecta time series |
| Coronal topology / open-flux constraints | Ambient field connectivity and amplitude | Does not validate eruptive downstream topology |
| SEP connection/release timing, in a later phase | Candidate connected shock regions | Does not determine injection/acceleration without a particle/source model |

Calibration observables and validation observables must be separated. Matching the same arrival time used to fit a drag coefficient is a fit, not an independent forecast test. Report approximation limits instead of adjusting field amplitude solely to force a desired Mach number.

A shock standoff measurement or radio band-split-derived compression/field estimate is useful only with its geometric, density, emission, and shock-model assumptions declared. A comparison derived through the same ambient or RH assumptions is a conditional consistency check, not an independent measurement of every inferred field/plasma quantity.

## 21. Bounded implementation sequence

This sequence implements the reduced mode without claiming completion of the earlier full-CME stages.

| Step | Implementation | Exit evidence |
| --- | --- | --- |
| S0 | Freeze reduced capabilities, asset/schema semantics, validity/coverage policies, and baseline exclusions | Parser/architecture tests; explicit unsupported downstream volume and scientific/numerical outcome separation |
| S1 | Reuse and audit canonical ambient magnetic/plasma queries | PFSS, EOS, orientation limits, field matching, sensitivity and calibration checks (RSH31, RSH37) |
| S2 | Reuse front geometry; implement complete normal-speed/support and selected width/direction histories | Plane, sphere, rotating ellipsoid, SSE edge/patch references; RSH32 for histories |
| S3 | Connect coronal history to maintained SWCME outer evolution | Full matched-state, one-sided acceleration receipts, pre/crossing/post references and proxy-validity gate RSH38 |
| S4 | Expose shared RH diagnostics and typed weak/sub-fast behavior | Conservative jumps, magnetic compression and HT references; RSH36 criticality subgate only if selected |
| S5 | Add exports, connectivity, observers, selected branch derivatives/ensembles, and background-only restart | Root/output/restart references; RSH39–RSH40 for the enabled capabilities |
| S6 | Integrate native SEP3D ambient-reference storage/front epochs | Actual one/four-rank owner/halo, zero-particle, transaction, and regression receipts |
| S7 | Complete actual 1-AU campaign and applicable event validation | Geometric endpoint and native coverage; explicit shock coverage outcome, ambient sensitivity and declared qualification level |

Reuse qualified functionality after checking its actual source and interface. Do not rebuild shared PFSS, RH, or geometry kernels merely to place them behind a new mode. The user's current local tree is authoritative for implementation; historical uploaded trees do not establish its present status.

Once the reduced SEP3D background works, `srcSEP` can consume the same shared provider with a thin application adapter and its own native validation. Particle work remains a separate later task. If the full sheath/ejecta mode is resumed, its reopened contact, inventory, momentum, induction, and interface requirements remain in force.

### 21.1 Required documentation and comments

Every new or materially changed algorithm documents:

- Its governing equation and which approximations it makes.
- SI units, coordinate frame, normal/sign convention, and field/temperature interpretation.
- Mathematical surface labels versus material labels.
- Applicable geometry, regularity, and supported parameter limits.
- Numerical tolerances, conditioning, branches, and typed failure behavior.
- Provider ownership, generation, preparation/commit, and MPI responsibilities.
- Reference versus actual downstream physical state.

The README links this specification, runnable translated decks, asset provenance, exact build/test commands, and limitations. The plan records implemented versus validated capabilities and exact evidence. Keep the full-CME roadmap intact and mark the reduced milestone separately.

### 21.2 Optional physical extensions requiring separate qualification

The following are research branches, not corrections that can be enabled implicitly to make a failing fixture pass:

| Extension | Physical work required | Separate validation |
| --- | --- | --- |
| Shock-aware outer excess-speed law | Justify section 9.4's target/excess evolution and ambient derivatives; identify supported regimes and fit independent front data | Trajectory integration, weak limits, event calibration; no automatic global-dynamics claim |
| Prescribed driver plus empirical standoff | Supply a CME-body/contact geometry history and a justified nose relation; couple shock position/speed to changing standoff | RSH33; independent body/shock observations and relation-validity limits |
| General surface dynamics / GSD | A surface representation plus a physically justified strength/area/angle closure; arbitrary mesh motion alone is insufficient | RSH34; the published MHD GSD reference is two-dimensional [R12], so a 3-D coronal extension needs new evidence |
| Nonlinear formation diagnostic | Define piston versus wave quantities, initial amplitudes/state, amplitude-dependent characteristic speeds, geometric dilution and supported dissipation/attenuation | RSH35; static-ambient linear rays alone cannot compute nonlinear steepening [R7] |
| Imported MHD ambient | A separate data/support/EOS/interpolation and divergence-preservation contract | Imported-data gates; optional to this analytical/reconstructed-ambient objective |

For the standoff branch, a commonly used bow-shock nose reference is [R11]

\[
\frac{\Delta}{R_c}
=0.81\,\frac{(\gamma-1)M^2+2}
{(\gamma+1)(M^2-1)}.
\]

Here `Rc` is the driver's nose radius of curvature; the Mach definition and applicability must match the selected reference. It is not a general per-flank MHD law. If `Rs=Rdriver+Delta` then `Vs=Vdriver+dot(Delta)`. Since `Delta=Rc F(M)` and `M` depends on shock speed, this is a coupled history/derivative problem, not necessarily a short local fixed-point iteration. The divergence as `M->1` signals loss of this quasi-steady relation's validity, not a demonstrated detachment/fade rule. No downstream volume is supplied by this extension.

For the general-front branch, distinguish weak characteristic speed from shock strength: with otherwise comparable normal upstream flow, a larger local `cf` gives faster weak fast-mode propagation. “All parts move faster in low-Alfvén-speed regions” is not a general wave-speed law. Strong-shock/driver response requires the actual closure and can behave differently.

For the formation branch, amplitude dependence is essential. For example, a one-dimensional right-moving ideal-gas simple wave has `c=c0+(gamma-1)u/2` and characteristic speed `u+c=c0+(gamma+1)u/2`. This is a limited steepening reference, not an MHD coronal prescription. Caustics of linear rays in a static ambient do not establish shock formation from an accelerating CME.

Keep one detailed specification with a reading guide. Repository READMEs can provide short user instructions and links; splitting physics and verification into separate documents is optional and must preserve consistent cross-references, schema versions, and scope.

## 22. Known limitations and responsible interpretation

1. The front is prescribed. Local admissibility and RH closure do not establish a global eruption or driven shock solution.
2. A finite triaxial or SSE surface is a low-dimensional shape approximation. It may miss real corrugation, asymmetric deflection, and ambient-induced distortion.
3. PFSS models ambient potential fields. It does not reconstruct all active-region currents or evolving eruption topology.
4. An ambient plasma profile can be empirically useful without satisfying global force balance. Its discrepancy is not removed by the jump solver.
5. Immediate downstream jump fields are not the magnetic field throughout the sheath or ejecta.
6. No retained sheath mass, contact impermeability, ejecta energetics, or CME magnetic flux budget is predicted.
7. A sub-fast geometric front is not a solved compression wave. A numerically unresolved super-fast jump is not a physical no-shock result.
8. Instantaneous ambient field-line connection is not a particle orbit or a guarantee of SEP escape to an observer.
9. A CME-body drag law used for shock-front motion is a named proxy approximation, with distinct validation needs.
10. Exact 1-AU apex evaluation does not ensure Earth lies within finite support, or that its local crossing has the apex parameters.
11. Native mesh/MPI correctness is separate from analytical correctness and event realism.
12. Passing this reduced qualification does not close the deferred full-CME background stages or enable particle physics automatically.
13. Fixed-width/direction SSE has zero normal speed at its tangent support edge; geometric extent therefore need not be shock extent.
14. Body-drag motion can become sub-fast or overtaken by the canonical wind. Smooth handoff and geometric arrival do not establish physical shock persistence.
15. Ambient temperature, density, open/closed routing, and flux normalization can dominate inferred Mach number. Synthetic normalization is not event calibration.
16. HT/criticality and connected-point derivatives have explicit conditioning/model/branch limits. Their absence is not evidence of zero particle acceleration.
17. Ensemble probabilities depend on declared measures and preserve unknown weight; they are not prior-independent confidence.

## 23. References and their specific roles

The specification's ownership, input/output, numerical-status, and qualification contracts are design requirements developed for this application. The papers below support the physical ingredients and comparisons; they do not prove this implementation correct or validate an arbitrary combination of histories.

**[R1] Rouillard, A. P., et al. (2016).** *Deriving the Properties of Coronal Pressure Fronts in 3D: Application to the 17 May 2012 Ground Level Enhancement.* The Astrophysical Journal, 833, 45. [DOI](https://doi.org/10.3847/1538-4357/833/1/45); [author manuscript](https://arxiv.org/abs/1605.05208). Role: precedent for combining finite reconstructed fronts with ambient field/plasma information to examine local shock properties and connectivity. The paper compares several ambient approaches; it is not evidence that PFSS captures a complete eruptive plasma field.

**[R2] Stansby, D., Yeates, A., and Badman, S. T. (2020).** *pfsspy: A Python package for potential field source surface modelling.* Journal of Open Source Software, 5(54), 2732. [Paper and DOI](https://joss.theoj.org/papers/10.21105/joss.02732). Role: PFSS equations, radial source-surface boundary, numerical magnetic reconstruction, and independent software reference. This specification does not require a Python runtime in the native provider.

**[R3] Stansby, D., and Verscharen, D. (2022).** *Test Problems for Potential Field Source Surface Extrapolations of Solar and Stellar Magnetic Fields.* [Author manuscript](https://arxiv.org/abs/2201.07783). Role: independent harmonic/flux references for verifying a PFSS implementation and tracing. Passing these tests establishes numerical behavior rather than event-specific realism.

**[R4] Möstl, C., and Davies, J. A. (2012).** *Speeds and arrival times of solar transients approximated by self-similar expanding circular fronts.* Solar Physics. [DOI](https://doi.org/10.1007/s11207-012-9978-8); [author manuscript](https://arxiv.org/abs/1202.1299). Role: SSE geometry and apex/flank arrival distinctions. The finite axisymmetric three-dimensional surface and support selection here are explicitly specified model extensions.

**[R5] Žic, T., Vršnak, B., and Temmer, M. (2015).** *Heliospheric Propagation of Coronal Mass Ejections: Drag-Based Model Fitting.* The Astrophysical Journal Supplement Series, 218, 32. [DOI](https://doi.org/10.1088/0067-0049/218/2/32); [author manuscript](https://arxiv.org/abs/1506.08582). Role: reduced outer propagation and fitting of trajectory parameters; not a derivation of a complete shock-driven plasma state.

**[R6] Dumbović, M., et al. (2021).** *Drag-Based Model (DBM) Tools for Forecast of Coronal Mass Ejection Arrival Time and Speed.* Frontiers in Astronomy and Space Sciences. [DOI](https://doi.org/10.3389/fspas.2021.639986); [author manuscript](https://arxiv.org/abs/2103.14292). Role: variations of drag/cone models and the limitation of using magnetic-body propagation as a shock-arrival proxy.

**[R7] Lulić, S., et al. (2013).** *Formation of Coronal Shock Waves.* Solar Physics. [Author manuscript](https://arxiv.org/abs/1303.2786). Role: distinction between piston expansion, propagating compression, and nonlinear shock formation. It supports the limitation that local super-fast classification is not a calculated formation history.

**[R8] Koval, A., and Szabo, A. (2008).** *Modified “Rankine–Hugoniot” shock fitting technique: Simultaneous solution for shock normal and speed.* Journal of Geophysical Research: Space Physics. [DOI and full text](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2008JA013337). Role: conservative local shock relations, frame/normal/speed consistency, and independent synthetic/observational jump checks. Shock fitting from measured states and the forward jump solve proposed here are different uses of the same conservation laws.

**[R9] Battarbee, M., Vainio, R., Laitinen, T., and Hietala, H. (2013).** *Injection of thermal and suprathermal seed particles into coronal shocks of varying obliquity.* Astronomy & Astrophysics, 558, A110. [DOI](https://doi.org/10.1051/0004-6361/201321348); [author manuscript](https://arxiv.org/abs/1309.2062). Role: the importance of obliquity and downstream treatment for a later particle-injection/acceleration model. It does not authorize enabling particle physics during reduced-background qualification.

**[R10] Linker, J. A., et al. (2017).** *The Open Flux Problem.* The Astrophysical Journal, 848, 70. [DOI](https://doi.org/10.3847/1538-4357/aa8a70); [author manuscript](https://arxiv.org/abs/1708.02342). Role: tension between modeled coronal-hole area and heliospheric open flux across magnetic maps/models; supports joint calibration and uncertainty assessment, not a universal amplitude correction.

**[R11] Farris, M. H., and Russell, C. T. (1994).** *Determining the standoff distance of the bow shock: Mach number dependence and use of models.* Journal of Geophysical Research, 99(A9), 17681–17689. [DOI](https://doi.org/10.1029/94JA01020). Role: empirical bow-shock nose standoff/Mach/curvature comparison for an optional driver branch. Its quasi-steady applicability is not evidence for an arbitrary evolving CME flank law.

**[R12] Mostert, W., Pullin, D. I., Samtaney, R., and Wheatley, V. (2017).** *Geometrical shock dynamics for magnetohydrodynamic fast shocks.* Journal of Fluid Mechanics, 811, R2. [DOI](https://doi.org/10.1017/jfm.2016.767); [publisher abstract](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/geometrical-shock-dynamics-for-magnetohydrodynamic-fast-shocks/2DC5A580D7CBBC64D563EFEC8FD25200). Role: a verified two-dimensional MHD area–Mach–angle formulation, relevant to a future dynamics branch. It does not supply a qualified 3-D coronal implementation.

**[R13] Edmiston, J. P., and Kennel, C. F. (1984).** *A parametric survey of the first critical Mach number for a fast MHD shock.* Journal of Plasma Physics, 32(3), 429–441. [DOI](https://doi.org/10.1017/S002237780000218X). Role: dependence of collisionless first criticality on beta, obliquity and EOS. A verified optional procedure/table is needed; ideal-MHD compression alone is not an injection model.

**[R14] Wilkinson, W. P., Humphreys, M., and Kilgallon, S. (2023).** *The Link Between Individual Velocities in the Incident Plasma Flow and Ion Populations Upstream of Collisionless Shocks.* Journal of Geophysical Research: Space Physics. [DOI and full text](https://doi.org/10.1029/2023JA031615). Role: explicit frame distinctions for incident flow and HT transformations. The vectors in section 11.7 are written for this specification's signed normal convention; singular-limit controls are implementation requirements.

**[R15] Heras, A. M., Sanahuja, B., Lario, D., Smith, Z. K., Detman, T., and Dryer, M. (1995).** *Three Low-Energy Particle Events: Modeling the Influence of the Parent Interplanetary Shock.* The Astrophysical Journal, 445, 497–508. [DOI](https://doi.org/10.1086/175714); [archived paper](https://adsabs.harvard.edu/full/1995ApJ...445..497H). Role: historical shock/observer magnetic-connection and cobpoint terminology. Only connectivity diagnostics are adopted here; its particle modeling is outside the reduced mode.

## 24. Specification checks performed for this document

The preparation checks, run while writing this specification, verify selected equations and illustrative numbers independently of AMPS:

- PFSS harmonic radial normalization and zero tangential field at the source surface.
- Cartesian ellipsoid normal-speed formula against direct time differentiation, including rotation.
- Finite SSE apex/ray support and the distinction between radial and normal speeds.
- Fixed-width SSE edge and analytical super-fast half-angles; changing-width/direction normal speeds against Cartesian differentiation.
- Quintic acceleration integral and endpoint matching.
- Constant-drag trajectory derivative and the illustrative handoff/1-AU event times.
- One-sided handoff acceleration, polar Parker/PFSS sensitivity values, parallel/perpendicular fast speeds, and the non-forward crossing of the fitted proxy.
- Nonaxisymmetric exterior magnetic divergence and induction for the explicitly stated reference assumptions.
- Hydrodynamic and finite-beta perpendicular RH examples using independently reconstructed conservation fluxes.
- HT boost versus incident-flow speed, electric-field cancellation, magnetic compression limits, flux-calibration negative controls, and connected-point derivative/ensemble-accounting identities.
- Markdown code/math delimiters, proposed INI examples, reference labels, and required capability declarations.

These are **document/equation checks**, not production source tests. No claim of native compilation, MPI execution, current local-source correction, or qualification follows from them. Implementation must supply the separate evidence in section 19.

For version 1.1, the original document/equation suite passed 104 checks and the review-correction suite passed 48 checks. Production RSH implementation tests, native compilation, and MPI qualification were not run as part of this document revision.

### 24.1 Review disposition

| Review topic | Incorporated correction | Qualification retained |
| --- | --- | --- |
| SSE edge and finite shock extent | Exact normal-speed limit, conditional half-angle table, patch output and width/direction histories | Zero edge speed applies to fixed width/direction; real shock extent requires local ambient/RH evaluation |
| Shock versus body drag | Weak fast-speed limit, wind discrepancy and proxy/coverage diagnostics | A shifted or variable target is a new closure, not a proven dynamical repair |
| Ambient Mach sensitivity | EOS/profile sensitivity and a reproducible polar example | Parallel fast speed corrects the review's table; its 10 G case is not sub-fast everywhere |
| Magnetic calibration | Joint coronal/IMF constraints and divergence/normal-flux checks | Independent radial scalar scaling is rejected |
| General front and standoff dynamics | Optional research branches and dedicated gates | Mesh freedom and per-patch drag do not supply a global shock closure; bow-shock nose relations have limits |
| Formation timing | Separate nonlinear-formation contract and test gate | Linear rays alone cannot establish nonlinear steepening |
| Additional diagnostics | Magnetic compression, correctly named HT vectors, optional criticality, connected-branch derivatives | Frame, null/conditioning, table domain and branch changes remain explicit |
| Handoff acceleration | One-sided -7.2 m/s² fixture jump and optional C2 test | C1 fixtures remain valid within their stated contract |
| Ambient improvements and document layout | Open/closed-aware profiles, uncertainty and reading/test guides | Imported MHD is optional; the requested detailed model remains one document |
