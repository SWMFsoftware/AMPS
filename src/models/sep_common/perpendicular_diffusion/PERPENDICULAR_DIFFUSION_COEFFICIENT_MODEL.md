# Perpendicular diffusion coefficient models for SEP and GCR transport

**Specification version:** 2.1  
**Revision date:** 8 October 2026  
**Companion:** `PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md`, specification 1.3  
**Purpose:** a mathematical specification for Codex to implement a reusable coefficient library, and a source-traceable foundation for publication methods sections. The implementation roadmap is the last section.

Specification 2.0 replaced 1.0 after the first scientific review. Revision 2.1 incorporates the supported findings of the review of version 2.0 and the primary-source checks recorded in Section 18.4. It fixes definitions, dimensions, limiting assumptions, spectrum conventions and the connection between coefficients and stochastic transport. It does not introduce a new physical closure. Algebraic reductions of published equations are identified as such. Implementation choices are requirements for the proposed software, not claims made by the cited authors.

There is no universal perpendicular coefficient for SEPs or GCRs. A prescription, an asymptotic turbulence closure, a running displacement model and a source-width model represent different observables and assumptions. The library MUST preserve those distinctions. It MUST not silently switch between them, add their coefficients, infer a turbulence partition, choose a calibration constant or impose a Bohm bound.

The withdrawn version 1.0 file's unverified observational numbers, partially garbled simulation tables, unattached-script claims and assertions about the current AMPS implementation are excluded from the normative specification. A claim about AMPS MUST be checked against a concrete checkout before it is reported as a code defect. No AMPS compilation, particle simulation or performance measurement was performed for this revision.

## Contents

1. [Evidence and scope](#1-evidence-and-scope)
2. [Observables, particle variables and units](#2-observables-particle-variables-and-units)
3. [Transport equation and tensor](#3-transport-equation-and-tensor)
4. [Turbulence spectra and provider contract](#4-turbulence-spectra-and-provider-contract)
5. [Field-line transport and particle limits](#5-field-line-transport-and-particle-limits)
6. [Prescribed SEP coefficients](#6-prescribed-sep-coefficients)
7. [NLGC and ENLGC](#7-nlgc-and-enlgc)
8. [Coupled NLGCE-N and NLGCE-F](#8-coupled-nlgce-n-and-nlgce-f)
9. [UNLT and implicit slab contribution](#9-unlt-and-implicit-slab-contribution)
10. [Field-line–particle decorrelation theory](#10-field-lineparticle-decorrelation-theory)
11. [Random ballistic decorrelation and composite closed form](#11-random-ballistic-decorrelation-and-composite-closed-form)
12. [Compound transport, retracing and pre-diffusive diagnostics](#12-compound-transport-retracing-and-pre-diffusive-diagnostics)
13. [Isotropic-turbulence fits and MHD scalings](#13-isotropic-turbulence-fits-and-mhd-scalings)
14. [GCR prescriptions and drift](#14-gcr-prescriptions-and-drift)
15. [Stochastic transport, gradients and geometry](#15-stochastic-transport-gradients-and-geometry)
16. [Library interface, solvers and tables](#16-library-interface-solvers-and-tables)
17. [Numerical verification and acceptance](#17-numerical-verification-and-acceptance)
18. [Publication reporting and revision record](#18-publication-reporting-and-revision-record)
19. [References and equation provenance](#19-references-and-equation-provenance)
20. [Implementation roadmap for Codex in AMPS](#20-implementation-roadmap-for-codex-in-amps)

## 1. Evidence and scope

### 1.1 Evidence categories

| Category | Meaning |
|---|---|
| Source equation | The cited primary paper supplies the equation. The source version and equation locator are given in Section 19. A preprint is identified when used. |
| Algebraic consequence | Follows from explicitly stated definitions or source equations; it is not an additional empirical law. |
| Numerical fixture | Recomputed from the specified mathematical model by the accompanying reference scripts. It tests transcription and implementation, not physical accuracy. |
| Calibration input | A parameter that MUST be supplied with a fit or publication reference. No value is inferred from missing information. |
| Source gate | An inherited variant lacks a sufficiently resolved formula, normalization, domain or parameter definition. Its production registry entry remains unavailable. |
| Engineering requirement | Proposed software behavior or numerical acceptance criterion. Its number is a chosen control, not a measured physical constant. |

The accompanying `PERPENDICULAR_DIFFUSION_COEFFICIENT_MODEL_DATA.zip` contains equation fixtures, reference calculations, the two NLGCE-F coefficient tables, source-check records and checksums. It does not contain AMPS validation results. The references support particular equations and assumptions; they are not a claim that every cited model is accurate throughout all SEP and GCR environments.


Requirements use MUST/REQUIRED for binding behavior, SHOULD for a recommendation whose alternative needs documented justification, and MAY for an explicit option. Imperative contract statements such as “Require”, “Reject”, “Use” and “Do not” are binding within their stated scope. Mathematical definitions, equations, source conventions, required inputs and model domains remain normative independently of keyword capitalization. These software requirements do not turn an engineering choice into a source-derived physical law.

### 1.2 Registry scope and readiness

“Specified” below means the equations and required parameters are stated here. It does not mean a production implementation exists or has passed scientific validation. A parameterized model can be complete while requiring user-supplied calibration inputs.

| Proposed identifier | Return kind | Geometry / input | Specification state |
|---|---|---|---|
| `constant_kappa_perp`, `constant_lambda_perp` | Local scalar coefficient | Explicit nonnegative constant | Specified, Section 6 |
| `ratio_kappa`, `power_law_perp` | Local scalar coefficient | Supplied parallel coefficient or normalization | Specified, Section 6 |
| `pitch_angle_perp`, `droge_lambda_perp`, `paradise_alpha` | Pitch-angle coefficient or explicitly averaged coefficient | Prescribed SEP model | Specified, Section 6 |
| `fl_slab`, `fl_2d`, `fl_composite` | Field-line coefficient in metres | Section 4 spectra | Specified closure, Section 5 |
| `flrw_particle` | Streaming-limit particle coefficient | Field-line coefficient and pitch-angle convention | Specified limit, Section 5 |
| `nlgc`, `enlgc_2d` | Asymptotic coefficient | Supplied parallel coefficient and spectra | Specified, Section 7 |
| `nlgce_n`, `nlgce_f_2014` | Coupled parallel–perpendicular pair | The specific Qin–Zhang spectra and conventions | Specified, Section 8 |
| `unlt`, `implicit_slab_exact_2016`, `implicit_slab_rational_2016` | Asymptotic coefficient | Parallel coefficient and named turbulence model | Specified, Section 9 |
| `flpd_complete` | Asymptotic coefficient | 2D spectrum and stated field-line closure | Specified, Section 10 |
| `rbd_bc` | Asymptotic coefficient | Parallel coefficient and spectral tensor | Specified RBD/BC kernel, Section 11 |
| `composite_closed_2019` | Asymptotic coefficient | Parallel coefficient, field-line coefficient and required tagged perpendicular length | Parameterized relation; named 2019 bend-over / 2022 integral reproductions in Section 11.2 |
| `compound_diffusive_lines`, `gcd_compound` | MSD and running coefficients | Parallel displacement statistics and static field lines | Specified diagnostic, Section 12 |
| `frozen_fieldline` | Persistent geometric realization | Field-line coefficient / resolved field | Transport-layer mechanism, Section 12 |
| `prediffusive_fit` | Conditional median displacement | Laitinen–Dalla fit parameters | Specified diagnostic; not a diffusion coefficient |
| `iso_fit_candia_roulet_2004` | Calibrated asymptotic pair | Isotropic random field plus ordered field | Specified named fit and declared software domain, Section 13.2 |
| `iso_fit_snodin_2016` | Calibrated asymptotic pair | Sharp-cutoff Kolmogorov spectrum; total-rms-field radius | Specified basic pair only for this named shape, Section 13.3 |
| `iso_kappa0_snodin_2016` | Isotropic scalar coefficient | Zero mean field and one named source spectral shape | Specified finite-range $\kappa_0$ fits; ordered-field extensions for other shapes remain gated |
| `iso_fit_kuhlen_2025` | Running coefficients and asymptotic pair | Isotropic turbulence; mandatory calibration data | Parameterized algorithm specified, Section 13 |
| `classical_scattering` | Coefficient pair and Hall coefficient | One scattering time | Specified relation, Section 13 |
| `yan_lazarian_ma4` | Order-of-magnitude estimate | Sub-Alfvénic MHD regime | Explicit scaling only, Section 13 |
| `nwu_ratio_polar`, `corti_ams02`, `helmod_ratio` | Two perpendicular eigenvalues | Explicit rigidity law and angular convention | Parameterized prescriptions, Section 14; exact paper presets require their source data |
| `drift_weak_scattering`, `drift_rigidity_reduction`, `drift_classical` | Signed antisymmetric coefficient | Explicit drift convention | Specified, Section 14 |
| `tabulated_perp` | Declared coefficient / pair | Validated, versioned table | Specified data contract, Section 16 |
| `qlt_perp` | Diagnostic | Named QLT orbit / resonance model | Source gate for a general particle closure; only Section 5 limits are specified |
| `unlt_fgr_perturbative_diagnostic` | Perturbative perpendicular-path diagnostic | Pure 2D; weak scattering and suppressed velocity diffusion | Specified (U10), with its small-correction domain |
| `unlt_fgr` | Full finite-gyroradius closure | Harmonics, resonances and source-specific assumptions | Source gate for an executable (U11) backend |
| `enlgc_slab_secant_diagnostic` | Running slab MSD-ratio estimator | Supplied parallel coefficient and elapsed time | Specified (N6); not an asymptotic coefficient |
| `nlgc_slab_kernel_diagnostic` | Mathematical slab-kernel coefficient | Pure slab; closure assumptions identified as physically inadequate | Specified mathematical (N2)/(N4) diagnostic; prohibited in a production diffusion mover |
| `drift_candia_roulet_2004` | Signed Hall coefficient | Ordered-field radius and source fit parameters | Specified formula (G10); separately required declared domain, Section 14.5 |
| `iso_fit_casse_2002`, `restricted_scattering` | Calibration-dependent fit | Source-specific amplitude and domain | Source gate for an absolute coefficient |
| WNLT, SOQLT, dynamic / compressive extensions | Source-specific result | Additional physics and spectra | Outside this executable specification; source gate |

NLGC, ENLGC, UNLT, the exact/rational implicit-slab models, FLPD, RBD and the composite closed form take a parallel coefficient or its consistently converted mean free path as input. NLGCE determines a pair; it MUST not consume an independently selected parallel model. A direct coefficient, its polynomial surrogate and its interpolation table are separate evaluation backends with separate error evidence.


Unless a limit is stated explicitly, particle-closure formulas assume $v>0$ and a finite positive supplied $\kappa_\parallel$. For NLGC/ENLGC, UNLT, the exact/rational implicit-slab models, FLPD, RBD/BC and the composite closed form, the bounds in Section 16.2 establish the continuous fixed-input limit $\kappa_\perp\to0$ as $\kappa_\parallel\to0$. With otherwise valid required inputs, these backends MUST return a successful, explicitly tagged zero coefficient at exactly zero supplied $\kappa_\parallel$, without forming $Q=v^2/(3\kappa_\parallel)$ or a log of zero. Geometry and source gates still apply; an invalid pure-slab production NLGC request is not made valid by this branch.

At that endpoint $\kappa_\perp/\kappa_\parallel$ is an undefined $0/0$ quotient. A separately derived limiting ratio MAY be reported as a limit, not as an evaluated ratio. This mathematical extension does not establish every microscopic assumption of a late-time closure at zero parallel transport. NLGCE owns its parallel output; this is not a rule for assigning an input zero to that pair. Field-line coefficients themselves do not require a particle speed.

## 2. Observables, particle variables and units

### 2.1 Per-axis and vector conventions

In homogeneous, zero-mean displacement statistics define

$$
M_x(t)=\langle[\Delta x(t)]^2\rangle,\qquad
d_x(t)=\frac12\frac{dM_x}{dt},\qquad
\bar d_x(t)=\frac{M_x(t)}{2t}.
\tag{O1}
$$

With a nonzero systematic displacement, use the centered covariance for diffusion and report the mean displacement separately. Do not interpret deterministic drift as symmetric diffusion.

For axisymmetric perpendicular transport,

$$
M_\perp=M_x+M_y=2M_x,\qquad
d_\perp=\frac14\frac{dM_\perp}{dt}=d_x,\qquad
\kappa_\perp=\lim_{t\to\infty}d_\perp(t),\qquad
\lambda_\perp=\frac{3\kappa_\perp}{v}.
\tag{O2}
$$

The library coefficient is per transverse axis, not the sum of two axes. The conversion $\lambda=3\kappa/v$ is a convention for the isotropic spatial diffusion limit. A pitch-angle-dependent $D_\perp(\mu)$ MUST be identified explicitly rather than treated as an already averaged $\kappa_\perp$.

If $M_x=A t^\alpha$, then

$$
d_x=\frac{\alpha A}{2}t^{\alpha-1}
=\alpha\bar d_x.
\tag{O3}
$$

The two estimators agree for linear MSD growth. They differ during ballistic, sub-diffusive or super-diffusive motion. A nonnegative time-dependent Markov coefficient chosen to match an MSD is $d_x$, not $\bar d_x$; matching one MSD does not reproduce the memory or conditional statistics of compound transport.

For a static line parameterized by the signed mean-field coordinate $z$, use

$$
M_{\rm FL}(z)=\langle[X(z)-X(0)]^2\rangle,\qquad
d_{\rm FL}(z)=\frac12\frac{dM_{\rm FL}}{d|z|},\qquad
M_{\rm FL}=2\kappa_{\rm FL}|z|\ \text{in its diffusive regime}.
\tag{O4}
$$

$\kappa_{\rm FL}$ has units of length. An arc-length definition is not automatically a mean-field-$z$ definition; its conversion MUST be justified for the geometry, especially at large turbulence amplitude.

### 2.2 Relativistic variables

Use total momentum magnitude $p$, rest mass $m$ and signed charge $q_e$. In SI,

$$
E_{\rm tot}=\sqrt{m^2c^4+p^2c^2},\quad
v=\frac{pc^2}{E_{\rm tot}},\quad
\gamma_{\rm rel}=\frac{E_{\rm tot}}{mc^2},\quad
\mathcal R=\frac{pc}{|q_e|},\quad
r_L=\frac{p}{|q_e|B_{\rm ref}},\quad
\Omega=\frac{|q_e|B_{\rm ref}}{\gamma_{\rm rel}m}=\frac{v}{r_L}.
\tag{O5}
$$

$r_L$ is the nominal radius formed with total momentum. The instantaneous perpendicular gyroradius is $r_L\sqrt{1-\mu^2}$. $B_{\rm ref}$ is $B_0$ for the ordered-field closures, but some isotropic fits use the total rms field. This tag is part of the model identity. Kinetic energy per nucleon, total kinetic energy and rigidity are separate inputs; conversions require the species mass and charge.

Require $m>0$, $p\ge0$ and $|q_e|>0$ for the charged-particle conversions, and $-1\le\mu\le1$ when pitch angle is supplied. The equality $\Omega=v/r_L$ applies for $p>0$; evaluate the mass-and-charge expression for $\Omega$ at $p=0$ instead of forming $0/0$. Closure-specific positive-radius requirements still apply.

### 2.3 Units and required tags

| Quantity | Internal unit / convention |
|---|---|
| $\kappa_\parallel,\kappa_{\perp i},D_\perp(\mu),\kappa_A$ | $\mathrm{m^2\,s^{-1}}$; $\kappa_A$ may be signed |
| $\lambda_\parallel,\lambda_{\perp i},\ell,L_U,L_c,\kappa_{\rm FL}$ | m, each with its definition |
| $B_0,\delta B$ | T; variance in $\mathrm{T^2}$ |
| Reduced magnetic spectrum $g(k)$ | $\mathrm{T^2\,m}$ |
| Area spectrum $S_2(k)$ | $\mathrm{T^2\,m^2}$ |
| $P_{ij}(\mathbf k)$ in $\int d^3k\,P_{ij}$ | $\mathrm{T^2\,m^3}$, including distributional factors |
| $\mathcal R$ | V; $1\,\mathrm{GV}=10^9\,\mathrm V$ |
| Angles | radians in calculations; degrees require an explicit conversion |
| $d_\perp$, $\bar d_\perp$ | Distinct running-result types |
| Conditional median, ensemble MSD, injection width | Distinct observable types |

Exact unit conversions used in the fixtures are $1\,\mathrm{au}=149\,597\,870\,700\,\mathrm m$, $c=299\,792\,458\,\mathrm{m\,s^{-1}}$, $1\,\mathrm{cm^2\,s^{-1}}=10^{-4}\,\mathrm{m^2\,s^{-1}}$ and $1\,\mathrm{nT}=10^{-9}\,\mathrm T$. The first two are definitions [R20,R21]. No proton mass or other measured constant is required for the dimensionless fixture tables in Section 17.

Require a positive speed for coefficient-to-mean-free-path conversion. Accept zero diffusion eigenvalues as exact zero. Models referencing an ordered field require $B_0>0$; an explicitly isotropic, zero-mean-field model can use a nonzero total rms field and does not require an arbitrary ordered direction.

### 2.4 Symbol and implementation-name map

Local source notation can be retained for traceability, but implementation variables MUST have an unambiguous meaning and unit. The following map resolves the principal collisions; subordinate indices and fit parameters are defined beside their equations.

| Symbol | Meaning | SI unit | Defining equation / canonical name |
|---|---|---|---|
| $Q_{\rm src}$ | Transport source term | $f$ per s | (T1), `transport_source` |
| $Q$ | Parallel scattering rate $v^2/(3\kappa_\parallel)$ | $\mathrm{s^{-1}}$ | (N2), `parallel_scattering_rate` |
| $s,q$ | Smooth-spectrum inertial and energy-range indices | dimensionless | (S4)/(S5), `inertial_index` / `energy_index` |
| $\nu_{\rm QZ}=s/2$ | Qin–Zhang source index | dimensionless | Section 4.2 / (E1), `qz_half_index` |
| $s_{\rm A}$ | Area-spectrum high-wavenumber inertial index | dimensionless | (S12)/(S13), `area_inertial_index` |
| $C(s),D(s,q)$ | Distinct spectral normalization functions | dimensionless | (S3), `slab_norm` / `two_d_norm` |
| $\kappa_{\rm FL},\kappa_s,\kappa_2$ | Field-line diffusion coefficients | m | (F2)–(F4), `fieldline_coefficient` and component names |
| $D_{\rm FL,iso}$ | Snodin isotropic zero-mean-field line coefficient | m | Section 13.3, `snodin_isotropic_fieldline_coefficient` |
| $D_\perp(\mu)$ | Pitch-angle-dependent particle coefficient | $\mathrm{m^2\,s^{-1}}$ | (P3), `pitch_perpendicular_kappa` |
| $\kappa_\parallel,\kappa_{\perp i},\kappa_A$ | Parallel, perpendicular-axis and signed Hall coefficients | $\mathrm{m^2\,s^{-1}}$ | (O2)/(T2)/(G8), explicit component names |
| $d_x,\bar d_x$ | Derivative and MSD-ratio estimators | $\mathrm{m^2\,s^{-1}}$ | (O1), `running_derivative` / `running_secant` |
| $M_x,\mathbf C$ | Uncentered axis MSD and centered displacement covariance | $\mathrm{m^2}$ | Section 2.1, `axis_msd` / `centered_covariance` |
| $z,\sigma$ | Signed mean-field coordinate and field-line arc length | m | (O4) and Section 13.4, separate coordinate tags |
| $\xi_x,\xi_y$ | Independent standard normal line increments | dimensionless | (C7), `line_normal_x` / `line_normal_y` |
| $\alpha_{\rm MSD}$ | MSD exponent (the $\alpha$ in O3) | dimensionless | (O3), `msd_exponent` |
| $\alpha_{\rm FLPD}$ | FLPD branch-consistency parameter | dimensionless | (D3)/(D9), `flpd_branch_parameter` |
| $\eta_\kappa,\eta_B,\epsilon^2$ | Coefficient ratio, isotropic turbulent-energy fraction, ordered-field variance ratio | dimensionless | (D2)/(I1)/(S2), separate variables |
| $A_{\rm FLPD}$ | Product of paths divided by $3\ell_2^2$ (the $A$ in D2) | dimensionless | (D2), `flpd_path_product` |
| $A_s,A_2$ | NLGCE decorrelation rates | $\mathrm{s^{-1}}$ | (E4), `nlgce_slab_rate` / `nlgce_two_d_rate` |
| $A_{\rm NRMHD}(k)$ | NRMHD potential spectrum (the $A(k)$ in U3) | $\mathrm{T^2\,m^4}$ | (U3), `nrmhd_potential_spectrum` |
| $A_K,C_K$ | Kuhlen parallel-fit amplitude (the $A$ in I8) and line-fit normalization | dimensionless | (I8)/(I11), separate fit parameters |
| $L_s=L_{c,s}$ | Source-defined slab correlation length | m | (S9)/(E2), `slab_correlation_length` |
| $\ell_2,L_U,L_\perp$ | Smooth 2D bend-over, ultra and source-defined integral scales | m | (S5)/(S9), distinct scale tags |
| $\ell_{\rm tc}$ | Heuristic transverse-complexity scale | m | (F7), `transverse_complexity_length` |
| $\ell_{\perp,\rm int}=L_\perp/2$ | Snodin 2022 total perpendicular integral scale for S5 | m | (B8), `snodin_total_perpendicular_integral_length` |
| $B_n,B_{\rm norm}$ | Model-specific normalization fields, never assumed numerically identical | T | (P7)/(G1)/(G5), record the selected reference field |
| $h_\theta,h_C$ | Inverse-angle transition parameters | $\mathrm{rad^{-1}}$ | (G4)/(G6), `polar_inverse_angle` |
| $\varepsilon_{ijk}$ | Levi-Civita tensor | dimensionless | Section 14.5, `levi_civita`; never turbulence amplitude |

Other appearances of $A$ as an MSD amplitude, $\xi$ as a resonance/kernel variable, $\beta=v/c$, and $h$ as a finite-difference step keep their local defining indices and units. The symbol table does not identify different source length conventions merely because their units coincide.

## 3. Transport equation and tensor

### 3.1 Parker equation and its scope

For a nearly isotropic phase-space distribution $f(\mathbf x,p,t)$ in the plasma frame,

$$
\frac{\partial f}{\partial t}
=\nabla\cdot(\mathbf K_s\nabla f)
-(\mathbf U+\mathbf v_d)\cdot\nabla f
+\frac13(\nabla\cdot\mathbf U)p\frac{\partial f}{\partial p}+Q_{\rm src},
\qquad
\mathbf v_d=\nabla\times(\kappa_A\mathbf b).
\tag{T1}
$$

Here $\mathbf b=\mathbf B_0/B_0$. The coefficient library supplies the diffusion and optional antisymmetric coefficients. The solver owns sources, boundaries, solar-wind convection, focusing, momentum changes, shocks and time direction. Equation (T1) is not a number-density Fokker–Planck equation in Cartesian phase-space volume; forward and backward stochastic formulations MUST be derived for the actual dependent variable and measure [R01].

In an orthonormal field-aligned frame $(\mathbf b,\mathbf e_1,\mathbf e_2)$,

$$
\mathbf K_s=\kappa_\parallel\mathbf b\mathbf b
+\kappa_{\perp1}\mathbf e_1\mathbf e_1
+\kappa_{\perp2}\mathbf e_2\mathbf e_2.
\tag{T2}
$$

For equal perpendicular eigenvalues,

$$
\mathbf K_s=\kappa_\perp\mathbf I
+(\kappa_\parallel-\kappa_\perp)\mathbf b\mathbf b.
\tag{T3}
$$

Nonnegative eigenvalues give a positive semidefinite tensor. It is positive definite only if every eigenvalue is strictly positive. With distinct perpendicular eigenvalues the orientation of $\mathbf e_1$ is physical input; a computational choice of an arbitrary perpendicular basis is sufficient only when they are equal.

### 3.2 Parker-spiral projection

For a field with no polar component, write $\mathbf b=\cos\psi\,\mathbf e_r-\sin\psi\,\mathbf e_\phi$ and choose the in-plane perpendicular direction $\mathbf e_n=\sin\psi\,\mathbf e_r+\cos\psi\,\mathbf e_\phi$. In the physical orthonormal spherical basis,

$$
\begin{aligned}
K_{rr}&=\kappa_\parallel\cos^2\psi+\kappa_{\perp r}\sin^2\psi,\\
K_{\phi\phi}&=\kappa_\parallel\sin^2\psi+\kappa_{\perp r}\cos^2\psi,\\
K_{r\phi}&=(\kappa_{\perp r}-\kappa_\parallel)\sin\psi\cos\psi,\\
K_{\theta\theta}&=\kappa_{\perp\theta}.
\end{aligned}
\tag{T4}
$$

The labels $\perp r$ and $\perp\theta$ name the in-plane and polar eigenvalues, not the radial and polar projections of a general tensor. For fields with $B_\theta\ne0$, assemble (T2) in Cartesian coordinates from the declared frame before projecting.

The full radial mean free path is

$$
\lambda_{\rm rad}=\frac{3K_{rr}}v
=\lambda_\parallel\cos^2\psi+\lambda_{\perp r}\sin^2\psi.
\tag{T5}
$$

The convention $\lambda_\parallel^r=\lambda_\parallel\cos^2\psi$ is only the parallel contribution to (T5). PARADISE's prescription in Section 6 uses that convention. Do not solve $\lambda_\parallel=\lambda_{\rm rad}/\cos^2\psi$ when the supplied radial quantity includes perpendicular transport.

For a unit shock normal $\mathbf n$, the projection is

$$
\kappa_n=\mathbf n\cdot\mathbf K_s\mathbf n
=\kappa_\parallel(\mathbf n\cdot\mathbf b)^2
+\kappa_\perp[1-(\mathbf n\cdot\mathbf b)^2]
\tag{T6}
$$

in the axisymmetric case. The shock model owns this projection and the acceleration calculation.

## 4. Turbulence spectra and provider contract

### 4.1 Geometry and magnetic variances

The magnetostatic composite geometry is

$$
\mathbf B=B_0\mathbf e_z+\mathbf b_s(z)+\mathbf b_2(x,y),\qquad
\mathbf b_2=\nabla\times[a(x,y)\mathbf e_z].
\tag{S1}
$$

The components are statistically independent and transverse, with axisymmetric transverse power unless stated otherwise. Define

$$
\delta B^2=\delta B_s^2+\delta B_2^2,\quad
f_s=\frac{\delta B_s^2}{\delta B^2},\quad
\epsilon^2=\frac{\delta B^2}{B_0^2},\quad
\delta B_{x,a}^2=\frac12\delta B_a^2\ (a=s,2).
\tag{S2}
$$

These are total transverse variances of each component. Isotropic three-dimensional turbulence is a different geometry: equal Cartesian component variances are $\delta B^2/3$, and there is no automatically defined slab/2D partition.

### 4.2 Smooth bend-over spectra and tensor normalization

For $s>1$, $q>-1$ and positive bend-over scales $\ell_s,\ell_2$ [R02,R03],

$$
C(s)=\frac{\Gamma(s/2)}{2\sqrt\pi\,\Gamma((s-1)/2)},\qquad
D(s,q)=\frac{\Gamma((s+q)/2)}{2\Gamma((s-1)/2)\Gamma((q+1)/2)},
\tag{S3}
$$

$$
g_s(k_\parallel)=\frac{C(s)}{2\pi}\delta B_s^2\ell_s
[1+(k_\parallel\ell_s)^2]^{-s/2},\quad -\infty<k_\parallel<\infty,
\tag{S4}
$$

$$
g_2(k_\perp)=\frac{2D(s,q)}\pi\delta B_2^2\ell_2
\frac{(k_\perp\ell_2)^q}{[1+(k_\perp\ell_2)^2]^{(s+q)/2}},\quad k_\perp\ge0.
\tag{S5}
$$

Use the distributional spectral tensors

$$
P_{xx}^s=\frac{g_s(k_\parallel)}{k_\perp}\delta(k_\perp),\qquad
P_{xx}^2=\frac{g_2(k_\perp)}{k_\perp}\delta(k_\parallel)
\left(1-\frac{k_x^2}{k_\perp^2}\right).
\tag{S6}
$$

The radial delta convention is fixed by

$$
\int d^3k\,P_{xx}^s=2\pi\int_{-\infty}^{\infty}g_s\,dk=\frac{\delta B_s^2}{2},\quad
\int d^3k\,P_{xx}^2=\pi\int_0^\infty g_2\,dk=\frac{\delta B_2^2}{2}.
\tag{S7}
$$

Implement the reduced one-dimensional integrals instead of a numerical delta at a quadrature node. $\nu=s/2$ in the older papers; $C(\nu)$ there equals $C(s)$ here. In particular the Qin–Zhang 2014 model uses $s=5/3$, $\nu=5/6$, and **flat energy-range spectra for both components**. Its 2D spectrum is not the $q=3$ spectrum used by many FLPD examples.

### 4.3 Spectral moments and lengths

For the 2D spectrum, its moment follows by the beta integral:

$$
I_m\equiv\int_0^\infty k^m g_2(k)\,dk
=\frac{D(s,q)}\pi\delta B_2^2\ell_2^{-m}
\frac{\Gamma((q+m+1)/2)\Gamma((s-m-1)/2)}{\Gamma((s+q)/2)}.
\tag{S8}
$$

This expression requires $q+m>-1$ and $s>m+1$. An unavailable divergent moment is a status, not a number. Thus $I_{-2}$ requires $q>1$, $I_{-1}$ requires $q>0$, and $I_1$ diverges at $s=5/3$ without an explicitly specified high-wavenumber cutoff. A numerical integration endpoint is not a physical dissipation scale.

Define the lengths used in [R02]:

$$
L_s=2\pi C(s)\ell_s,\quad
L_U=\sqrt{\frac{s-1}{q-1}}\ell_2,\quad
L_\perp=
\frac{2\Gamma(q/2)\Gamma(s/2)}{\Gamma((q+1)/2)\Gamma((s-1)/2)}\ell_2.
\tag{S9}
$$

$L_U$ exists for $q>1$ and $L_\perp$ for $q>0$. Bend-over scale, integral scale, outer scale and ultra-scale are distinct tagged quantities. Convert (S9) only for this spectral shape; a different shape requires its own integral definitions. When nonzero variance is absent, report a shape-defined length separately from an undefined correlation estimator.

The Kubo diagnostic used in the ordered-field literature is

$$
K=\frac{\ell_\parallel}{\ell_\perp}\frac{\delta B_x}{B_0}.
\tag{S10}
$$

It requires a declared pair of reference lengths and an amplitude convention. It is not interchangeable with a convention using total transverse rms amplitude. A composite spectrum has more than one plausible diagnostic length; none is selected silently.

### 4.4 Area-spectrum adapters

For an axisymmetric 2D area density normalized by
$\int_0^\infty 2\pi k S_2(k)\,dk=\delta B_2^2$, the consistent reduced-spectrum relation is

$$
g_2(k)=kS_2(k).
\tag{S11}
$$

Then (S7) gives the component variance $\delta B_2^2/2$. This fixes the factors of $k$ and $2\pi$ at the adapter boundary.

The single-break spectrum used by Chhiber et al. [R04] is

$$
S_2(k)=C_2\delta B_2^2\lambda_2^2
\begin{cases}
(k\lambda_2)^p,&k\lambda_2\le1,\\
(k\lambda_2)^{-s_{\rm A}-1},&k\lambda_2>1,
\end{cases}
\quad
C_2=\frac{(s_{\rm A}-1)(p+2)}{2\pi(p+s_{\rm A}+1)}.
\tag{S12}
$$

Here $s_{\rm A}>1$, $p>-2$, and $\lambda_2$ is the break scale. The low-energy-range exponent of the reduced $g_2$ is $q=p+1$; the high-$k$ index is $s=s_{\rm A}$. These shapes share exponents with (S5), not the same full spectrum or length conversion.

For an explicitly supplied three-range area spectrum, one can normalize it algebraically without assigning a physical default:

$$
S_2(k)=\frac{C_0\lambda_2\delta B_2^2}{2\pi k}
\begin{cases}
(\lambda_2/\lambda_o)^{-1}(\lambda_o k)^q,& k<1/\lambda_o,\\
(\lambda_2 k)^{-1},&1/\lambda_o\le k<1/\lambda_2,\\
(\lambda_2 k)^{-s_{\rm A}},& k\ge1/\lambda_2,
\end{cases}
\quad
C_0=\left[\frac1{q+1}+\ln\frac{\lambda_o}{\lambda_2}+\frac1{s_{\rm A}-1}\right]^{-1}.
\tag{S13}
$$

Require $\lambda_o\ge\lambda_2>0$, $q>-1$, $s_{\rm A}>1$. Equation (S13) is a normalized **supplied spectrum definition**, not a certified reproduction of a particular turbulence-transport paper or a recommended choice of outer scale. Use (S11) and compute its moments directly. $\lambda_2$ is a break, not automatically $L_\perp$.

### 4.5 Turbulence energy and partition accounting

If a provider defines $Z^2=\langle|\delta\mathbf v|^2+|\delta\mathbf b_A|^2\rangle$ with $\delta\mathbf b_A=\delta\mathbf B/\sqrt{\mu_0\rho}$ and $r_A=\langle\delta v^2\rangle/\langle\delta b_A^2\rangle$, then

$$
\delta B^2=\frac{\mu_0\rho Z^2}{1+r_A}.
\tag{S14}
$$

This identity depends on the provider's definition of $Z^2$. For Alfvén waves in kinetic–magnetic equipartition, if $w_\pm$ is the **total** wave energy density of each propagation population,

$$
\delta B_\pm^2=\mu_0 w_\pm.
\tag{S15}
$$

If instead $w$ is magnetic energy alone, $\delta B^2=2\mu_0w$. For unhalved Elsasser variables $\mathbf z^\pm=\delta\mathbf v\pm\delta\mathbf b_A$, a pure population has $w_\pm=\rho\langle|\mathbf z^\pm|^2\rangle/4$; a convention using half-Elsasser variables has $w_\pm=\rho\langle|\mathbf z^\pm/2|^2\rangle$. Record which convention is supplied.

A propagation-direction split is not a slab/2D geometry split. The provider MUST state whether its preexisting variance means slab variance or total variance. Adding a separately supplied 2D population preserves an already supplied slab variance. Repartitioning a fixed total variance gives $\delta B_s^2=f_s\delta B^2$ and changes a slab-QLT coefficient proportionally to $1/f_s$ with other inputs fixed; these are different operations. The frequently used partition $f_s=0.2$ is an illustrative model assumption, never an automatic default.

Required provider fields depend on the selected model: mean-field definition and vector; geometry; component variances; spectra and spectral-index conventions; tagged scales; physical spectral cutoffs if present; position, time and revision; energy conversion convention; and provenance. Missing 2D information MUST produce `missing_input` in a model requiring it. Preserve the original parallel-library normalization in the shared adapter.

## 5. Field-line transport and particle limits

### 5.1 Slab, 2D and composite field-line coefficients

For pure slab turbulence with the mean field along $z$, the correlation representation is

$$
M_{\rm FL}(z)=\frac2{B_0^2}\int_0^{|z|}(|z|-u)R_{xx}(u)\,du.
\tag{F1}
$$

The short-distance form is $M_{\rm FL}=\delta B_x^2 z^2/B_0^2$. For a general non-slab field, $R_{xx}$ in (F1) is sampled along the stochastic field line; replacing it by the straight-line Eulerian correlation requires a closure.

For (S4), the slab asymptotic coefficient is [R02,R03]

$$
\kappa_s=\pi C(s)\ell_s\frac{\delta B_s^2}{B_0^2}
=L_s\frac{\delta B_{x,s}^2}{B_0^2}.
\tag{F2}
$$

Under the diffusive-decorrelation field-line closure with a finite 2D ultra-scale,

$$
\kappa_2^2=\frac\pi{B_0^2}\int_0^\infty\frac{g_2(k)}{k^2}\,dk
=\frac{s-1}{2(q-1)}\ell_2^2\frac{\delta B_2^2}{B_0^2},\qquad q>1,
\tag{F3}
$$

and

$$
\kappa_{\rm FL}=\kappa_s+\frac{\kappa_2^2}{\kappa_{\rm FL}},\qquad
\kappa_{\rm FL}=\frac{\kappa_s+\sqrt{\kappa_s^2+4\kappa_2^2}}2.
\tag{F4}
$$

Equation (F4) is a closure, not an exact general property of field lines. It is not additive in the two component coefficients. It reduces to $\kappa_s$ with no 2D power and to $\kappa_2$ with no slab power. Return zero for vanishing turbulence without evaluating undefined fractions. A flat 2D spectrum has divergent (F3); a numerical cutoff MUST not conceal this divergence.

### 5.2 Streaming FLRW conversion

For unbackscattered streaming along a static diffusive line, at fixed pitch angle

$$
D_\perp(\mu)=v|\mu|\kappa_{\rm FL},\qquad
\kappa_\perp=\frac12\int_{-1}^1D_\perp(\mu)\,d\mu
=\frac v2\kappa_{\rm FL},\qquad
\lambda_\perp=\frac32\kappa_{\rm FL}.
\tag{F5}
$$

The factor $1/2$ is the isotropic $\langle|\mu|\rangle$. A beam with $\mu=1$ has $D_\perp=v\kappa_{\rm FL}$. An explicit phenomenological reduction $a_{\rm FL}$ multiplies (F5) if selected; it is distinct from the $a^2$ of NLGC. Diffusive parallel motion on the same static line produces compound sub-diffusion instead of (F5).

For a local longitudinal angular diffusion coefficient satisfying $\langle\Delta\phi^2\rangle=2D_\phi|z|$,

$$
\kappa_{\rm FL,\phi}=(r\sin\theta)^2D_\phi.
\tag{F6}
$$

The commonly written $r^2D_\phi$ is its equatorial form. This converts longitudinal arc displacement, not automatically the direction normal to a Parker spiral; a field-aligned projection MUST also be stated. Require the radius, colatitude, distance parameter and angular units at the adapter boundary.

### 5.3 Regime limits as diagnostics

The heuristic fluid and CLRR estimates are [R02,R18]

$$
\frac{\kappa_\perp}{\kappa_\parallel}\simeq\frac{\delta B_x^2}{B_0^2}
\quad\text{(ballistic field lines, diffusive parallel motion)},\qquad
\frac{\kappa_\perp}{\kappa_\parallel}\simeq\frac{\kappa_{\rm FL}^2}{\ell_{\rm tc}^2}
\quad\text{(CLRR assumptions)}.
\tag{F7}
$$

Here $\delta B_x^2$ is the Cartesian variance of the field in the selected heuristic, and $\ell_{\rm tc}$ is its transverse-complexity scale. R18 initially identifies this scale with the perpendicular bend-over length and explicitly allows another stated scale, such as an integral scale. No scale is inferred automatically. Equation (F7) is informational regime scaling, not a selectable numerical coefficient backend without that definition and the applicable assumptions. For FLPD specifically the variance in (D2)–(D5) is $\delta B_{x,2}^2$, not the sum of slab and 2D Cartesian variances.

The CLRR estimate in (F7) is not the general composite-FLPD short-path limit. At fixed moment $I_{-1}$, (D4) varies inversely with $\kappa_{\rm FL}^2$; (F7) at fixed $\ell_{\rm tc}$ varies directly with it. Section 10 gives the restricted pure-2D matching identity. These are regime statements, not a universal switch controlled only by $\lambda_\parallel$. No coefficient is changed at an inferred transverse-complexity threshold unless the selected model specifies that transition.

For Gaussian parallel motion with variance $2\kappa_\parallel t$ on a diffusive static line, the exact moment consequence is

$$
M_x(t)=4\kappa_{\rm FL}\sqrt{\frac{\kappa_\parallel t}{\pi}},\qquad
d_x(t)=\kappa_{\rm FL}\sqrt{\frac{\kappa_\parallel}{\pi t}},\qquad
\bar d_x(t)=2\kappa_{\rm FL}\sqrt{\frac{\kappa_\parallel}{\pi t}}.
\tag{F8}
$$

This is an anomalous MSD law with $\kappa_\perp=0$ asymptotically under continued attachment, not a finite perpendicular diffusion coefficient. Its $t\to0$ singular running coefficient shows that a late-time diffusive-parallel approximation is not an early-time particle model.

## 6. Prescribed SEP coefficients

### 6.1 Constants, ratios and power laws

The basic definitions are

$$
\kappa_\perp=\kappa_0,\qquad
\kappa_\perp=\frac{v\lambda_0}{3},\qquad
\kappa_\perp=\eta_\kappa\kappa_\parallel,
\tag{P1}
$$

with explicitly supplied nonnegative $\kappa_0$, $\lambda_0$ or $\eta_\kappa$. The ratio of mean free paths equals the ratio of coefficients only for the same particle speed and the same averaging convention. Do not invent a default ratio or a universal upper bound $\eta_\kappa\le1$.

A parameterized prescription can be defined by

$$
\kappa_{\perp i}=\kappa_{0i}\,\beta^{a_{vi}}
\left(\frac{\mathcal R}{\mathcal R_0}\right)^{a_{Ri}}
\left(\frac{B_n}{B_0}\right)^{a_{Bi}}
\left(\frac r{r_0}\right)^{a_{ri}},\qquad \beta=v/c.
\tag{P2}
$$

Every exponent, normalization and reference value is configuration input. A factor omitted by the selected prescription has exponent zero. Differentiation MUST account for the declared independent particle variable and the provider's actual spatial fields. No power law is inferred from “SEP” or “GCR.”

### 6.2 Pitch-angle normalization

For the three forms compared by Strauss and Fichtner [R05],

$$
D_\perp(\mu)=D_0\Phi(\mu),\qquad
\Phi(\mu)=1,\quad 2|\mu|,\quad \frac4\pi\sqrt{1-\mu^2},\qquad
\frac12\int_{-1}^1\Phi\,d\mu=1.
\tag{P3}
$$

They all give $\kappa_\perp=D_0$ under an isotropic pitch-angle average. This normalization does not assert equivalence for an anisotropic particle population. A focused solver evaluates the selected $D_\perp(\mu)$ at the particle's current $\mu$.

If a prescription instead uses $D_\perp=v\lambda_\perp^0 g(\mu)$, obtaining $\lambda_\perp^0=3\kappa_\perp/v$ requires $\langle g\rangle=1/3$. Thus $g=\mu^2$ and $g=(2/3)|\mu|$ have the required normalization; bare $|\mu|$ does not.

The reduced streaming prescription is

$$
D_\perp(\mu)=a_{\rm FL}v|\mu|\kappa_{\rm FL},\qquad
\kappa_\perp=\frac{a_{\rm FL}v}{2}\kappa_{\rm FL},\qquad
\lambda_\perp=\frac32a_{\rm FL}\kappa_{\rm FL}.
\tag{P4}
$$

$a_{\rm FL}$ is supplied explicitly. It MUST not be squared accidentally or equated to the NLGC $a^2$ parameter.

### 6.3 Dröge–Dresing form

The source prescription [R06] is

$$
\lambda_\perp(\mu,r)=\alpha_D\lambda_\parallel(r)
\left(\frac r{r_0}\right)^2\cos\psi(r)\sqrt{1-\mu^2},\qquad
D_\perp(\mu,r)=\frac v3\lambda_\perp(\mu,r).
\tag{P5}
$$

The averaged coefficient ratio is the algebraic consequence

$$
\frac{\kappa_\perp}{\kappa_\parallel}
=\frac\pi4\alpha_D\left(\frac r{r_0}\right)^2\cos\psi.
\tag{P6}
$$

$\alpha_D$ is the unaveraged prefactor. A Parker-spiral assumption and its wind speed and rotation rate belong to the geometry configuration. With $B_r\propto r^{-2}$ and $B=B_r/\cos\psi$, the factor in (P6) is proportional to $1/B$, but it does not equal unity at $r=r_0$ unless the spiral angle vanishes. Keep the event-specific source value of $\alpha_D$ out of the model's defaults.

### 6.4 PARADISE form

Using the convention $\lambda_\parallel^r=\lambda_\parallel b_r^2$, the isotropic prescription in Wijsen et al. [R07] is

$$
\kappa_\perp=\frac\pi{12}\alpha_P v
\frac{\lambda_\parallel^r}{b_r^2}\frac{B_{\rm norm}}{B_0}
=\frac\pi4\alpha_P\kappa_\parallel\frac{B_{\rm norm}}{B_0}.
\tag{P7}
$$

$B_{\rm norm}$ is a reference field, not the local mean-field symbol in a second role. Evaluate the second form if $\lambda_\parallel$ is available: it avoids an artificial $0/0$ near $b_r=0$. A supplied radial path including the perpendicular contribution in (T5) is incompatible with the first expression without additional information.

Prescribed injection longitude or latitude width is an input to the source $Q_{\rm src}$ in (T1), not a perpendicular coefficient. The conditional diagnostic in Section 12.4 MUST not be used to determine such a width.

## 7. NLGC and ENLGC

### 7.1 NLGC kernel and reduced integral

For the magnetostatic NLGC closure [R08,R09],

$$
\kappa_\perp=\frac{a^2v^2}{3B_0^2}
\int d^3k\,\frac{P_{xx}(\mathbf k)}
{v/\lambda_\parallel+\kappa_\parallel k_\parallel^2+\kappa_\perp k_\perp^2},
\quad \lambda_\parallel=3\kappa_\parallel/v.
\tag{N1}
$$

The NLGC $a^2$ is an explicitly configured closure factor; the original commonly used value is $1/3$. Choosing another source-specific factor defines that variant. No general strong-turbulence cutoff follows merely from the denominator of (N1). The velocity-variance restriction of RBD in Section 11 does not restrict the NLGC integral.

For (S4)–(S6), set $Q=v/\lambda_\parallel=v^2/(3\kappa_\parallel)$. The exact spectral reduction of (N1) is

$$
\kappa_\perp=\frac{a^2v^2}{3B_0^2}
\left[
4\pi\int_0^\infty\frac{g_s(k)}{Q+\kappa_\parallel k^2}\,dk
+\pi\int_0^\infty\frac{g_2(k)}{Q+\kappa_\perp k^2}\,dk
\right].
\tag{N2}
$$


This reduction fixes the factors for a two-sided slab spectrum and a radial 2D spectrum. The composite production backend retains the published slab term in (N2). Its mathematical pure-slab limit is

$$
\eta_{\kappa,\rm slab}^{\rm NLGC}
=\frac{a^2}{2}\frac{\delta B_s^2}{B_0^2}
\mathcal H\!\left(s,0,\frac{\lambda_\parallel^2}{3\ell_s^2}\right)>0
\quad(\lambda_\parallel>0,\ \delta B_s^2>0,\ a^2>0).
\tag{N8}
$$

This nonzero asymptotic prediction is spurious for continued attachment in pure magnetostatic slab turbulence, where (F8) is sub-diffusive. A production asymptotic `nlgc` evaluation MUST reject pure slab as outside its intended physical domain. The explicitly selected `nlgc_slab_kernel_diagnostic` MAY evaluate the mathematical kernel, but its output is barred from a production diffusion mover.

As 2D variance tends to zero at fixed slab variance, the composite formula approaches (N8). The eligibility gate does not repair its near-slab physical reliability and does not establish a quantitative admissible threshold in $f_s$. None is supplied here. Do not delete, taper or clip the slab term under the NLGC identity. ENLGC and UNLT retain their separately defined equations and assumptions.

An additional $\gamma(\mathbf k)$ in the denominator is a **named dynamical decorrelation model**, with dimensions $\mathrm{s^{-1}}$. This specification implements $\gamma=0$. A dynamic backend requires a source-specific rate and consistency with the parallel model and transport frame.

### 7.2 Hypergeometric reduction as an independent check

For $h_{s,q}(x)=x^q/(1+x^2)^{(s+q)/2}$ define

$$
\mathcal H(s,q,A)=4D(s,q)\int_0^\infty\frac{h_{s,q}(x)}{1+Ax^2}\,dx
=\frac{s-1}{s+q}\,{}_2F_1\left(1,\frac{q+1}2;\frac{s+q}2+1;1-A\right).
\tag{N3}
$$

Equation (N3) follows by the substitution $u=x^2/(1+x^2)$; $\mathcal H(s,q,0)=1$. It requires $s>1$, $q>-1$, $A\ge0$. With $\eta_\kappa=\kappa_\perp/\kappa_\parallel$,

$$
\eta_\kappa=\frac{a^2}{2}
\left[
\frac{\delta B_s^2}{B_0^2}\mathcal H\left(s,0,\frac{\lambda_\parallel^2}{3\ell_s^2}\right)
+\frac{\delta B_2^2}{B_0^2}\mathcal H\left(s,q,\frac{\lambda_\parallel\lambda_\perp}{3\ell_2^2}\right)
\right].
\tag{N4}
$$

This is an algebraic check of (N2), not an additional approximation. For extreme arguments a special-function implementation MUST be checked against quadrature or appropriate asymptotics rather than assumed stable.

### 7.3 ENLGC and flat-energy-range asymptotics

The asymptotic 2D contribution in ENLGC [R10] uses $a^2=1$ and, for the original flat 2D spectrum $q=0$, reads

$$
\kappa_\perp=\frac{2v^2C(s)\ell_2\delta B_2^2}{3B_0^2}
\int_0^\infty
\frac{[1+(k\ell_2)^2]^{-s/2}}{Q+\kappa_\perp k^2}\,dk.
\tag{N5}
$$

The pure-slab contribution to the published running MSD-ratio estimator is

$$
\bar d_x^s(t)=2\sqrt\pi C(s)\ell_s\frac{\delta B_s^2}{B_0^2}
\sqrt{\frac{\kappa_\parallel}{t}}.
\tag{N6}
$$

The registry identity is `enlgc_slab_secant_diagnostic`. Equation (N6) is the secant estimator in (F8) with $\kappa_{\rm FL}=\kappa_s$. The derivative coefficient is half this value. ENLGC's asymptotic scalar result uses the 2D term; the decaying slab contribution is reported through a separate running-result API.

Let $A=\lambda_\parallel\lambda_\perp/(3\ell_2^2)$. For (N5),

$$
\lambda_\perp\simeq\frac12\frac{\delta B_2^2}{B_0^2}\lambda_\parallel
\quad(A\ll1),\qquad
\lambda_\perp\simeq
\left[\sqrt3\pi C(s)\ell_2\frac{\delta B_2^2}{B_0^2}\right]^{2/3}
\lambda_\parallel^{1/3}
\quad(A\gg1).
\tag{N7}
$$

For a flat-2D NLGC variant, replace the factor inside the second bracket by $a^2$ times itself. The $\kappa_\parallel^{1/3}$ dependence is this asymptotic result, not the sensitivity of every NLGC model or every energy-range spectrum. A $q=3$ 2D spectrum generally has a different large-$A$ asymptotic behavior. Compute the full implicit derivative or solve the full integral for a general configuration.

## 8. Coupled NLGCE-N and NLGCE-F

### 8.1 Fixed spectrum and corrected parameters

Use exactly the magnetostatic model in Qin and Zhang (2014) [R11], also specified in the parallel companion. For each component use a one-sided total-variance density

$$
\mathcal P_a(k)=4C(s)\ell_a\delta B_a^2[1+(k\ell_a)^2]^{-s/2},\qquad
\int_0^\infty\mathcal P_a\,dk=\delta B_a^2,
\quad a=s,2,\quad s=5/3.
\tag{E1}
$$

For the slab component $k=|k_\parallel|$. For 2D it is the reduced radial wavenumber after angular integration, not a Cartesian slice. Both component spectra in this particular backend have $q=0$.

The corrected NLPA parameter is

$$
L_{c,s}=2\pi C(s)\ell_s,\qquad
\xi=\frac{2\pi r_L}{L_{c,s}\epsilon}
=\frac{r_L/\ell_s}{C(s)\epsilon},\qquad
a_x=\frac12\left[
\frac{f_s}{\dfrac{\xi}{1+\xi}\dfrac1\epsilon+\dfrac{\epsilon}{2\xi}}
\right]^{1/2}.
\tag{E2}
$$

The first term is $(\xi/(1+\xi))/\epsilon$, **not** $(\xi/(1+\xi))^{1/\epsilon}$. $L_{c,s}$ is the slab integral length; do not substitute a variance-weighted composite correlation length. $r_L$ uses $B_0$. The modified perpendicular closure is

$$
a'^2=\left[
\sqrt{\frac{\ell_2}{\ell_s}}\frac1{f_s}
+\frac4{3(1-f_s)}
\right]^{-1}.
\tag{E3}
$$

Require $0<f_s<1$, $\epsilon>0$, $r_L>0$ and positive scales. Source-specific zero-component limits are not obtained by evaluating singular fractions.

### 8.2 Coupled integral equations

Define

$$
A_s(k)=\frac{v^2}{3\kappa_\parallel}+k^2\kappa_\parallel,\qquad
A_2(k)=\frac{v^2}{3\kappa_\parallel}+k^2\kappa_\perp.
\tag{E4}
$$

The simultaneous pair is

$$
\kappa_\parallel^{-1}=\frac{3a_x\Omega^2}{v^2B_0^2}
\left[
\int_0^\infty\mathcal P_s(k)\frac{A_s(k)}{\Omega^2+A_s(k)^2}\,dk
+\int_0^\infty\mathcal P_2(k)\frac{A_2(k)}{\Omega^2+A_2(k)^2}\,dk
\right],
\tag{E5}
$$

The coefficient in (E5) is $a_x$ to the first power, although its definition (E2) contains a square root. It is not $a_x^2$ and not the unrelated NLGC factor $a^2$.


$$
\kappa_\perp=\frac{a'^2v^2}{6B_0^2}
\int_0^\infty\frac{\mathcal P_2(k)}{A_2(k)}\,dk.
\tag{E6}
$$

Equations (E5)–(E6) use the total-variance normalization in (E1); each Cartesian transverse variance is half that total. The tensor form of NLPA has $A/[\Omega^2+(A+\Gamma_{\rm dec})^2]$; the numerator remains $A$. This backend fixes $\Gamma_{\rm dec}=0$ and does not define a dynamic extension.

Solve positive logarithmic unknowns for the pair. Require each integral residual, quadrature convergence and agreement from multiple initial guesses. The delivered representative states A–G and actual starting-guess/refinement procedure are enumerated in Section 17.5; additional domain coverage belongs to Section 17.6 and MUST be listed before it is claimed. An NLGCE-F pair may be an initial guess but does not override the direct integral result. The shared paired backend MUST have one owner in the joint coefficient library; neither library SHOULD solve the pair independently or overwrite its parallel component.

### 8.3 Polynomial approximation and coefficient assets

With $\alpha\in\{\parallel,\perp\}$ and natural-log inputs

$$
x_1=\ln\frac{r_L}{\ell_s},\quad
x_2=\ln f_s,\quad
x_3=\ln\epsilon^2,\quad
x_4=\ln\frac{\ell_s}{\ell_2},
\tag{E7}
$$

the published polynomial is

$$
F_\alpha=\sum_{i=0}^{5}\sum_{j=0}^{3}\sum_{k=0}^{3}\sum_{l=0}^{2}
d^\alpha_{ijkl}x_1^ix_2^jx_3^kx_4^l,\qquad
\lambda_\alpha=\ell_s e^{F_\alpha},\qquad
\kappa_\alpha=\frac v3\lambda_\alpha.
\tag{E8}
$$

The companion data files `NLGCE_F_2014_parallel.csv` and `NLGCE_F_2014_perpendicular.csv` store the 288 coefficients for each component from source Tables 3 and 4. Their header is `j,k,l,d_i0,d_i1,d_i2,d_i3,d_i4,d_i5`. There are 48 unique $(j,k,l)$ rows per component. These are the same source-checked tables delivered with parallel specification 1.3, reused unchanged. Never reconstruct missing entries, refit a table silently or transpose an inferred row ordering.

The published input box is

$$
10^{-5}\le r_L/\ell_s\le6.3,\quad
10^{-3}\le f_s\le0.85,\quad
10^{-4}\le\epsilon^2\le10^2,\quad
1\le\ell_s/\ell_2\le10^3.
\tag{E9}
$$

An input inside this box is evaluable under the published fit definition; it does not establish a uniform approximation-error bound. The corrected parallel-companion audit finds state-dependent differences between the published polynomial and the corrected integral closure. Preserve both backends and their provenance; do not “correct” a fit result with an undocumented multiplier. No universal turbulence-only accuracy threshold is specified.

For independent variables $u_a$, logarithmic sensitivities follow from

$$
\frac{\partial\ln\lambda_\alpha}{\partial\ln u_a}
=\frac{\partial\ln\ell_s}{\partial\ln u_a}
+\sum_{n=1}^4\frac{\partial F_\alpha}{\partial x_n}
\frac{\partial x_n}{\partial\ln u_a},\qquad
\frac{\partial\ln\kappa_\alpha}{\partial\ln u_a}
=\frac{\partial\ln v}{\partial\ln u_a}
+\frac{\partial\ln\lambda_\alpha}{\partial\ln u_a}.
\tag{E10}
$$

Do not independently differentiate $p$, $v$, $r_L$ and rigidity and then add them as if they were independent particle inputs. Compare indexed coefficient content for scientific transcription; a raw CSV checksum and a checksum of canonically serialized rows serve different purposes.

## 9. UNLT and implicit slab contribution

### 9.1 Diffusive UNLT

The source kernel [R02,R09] is

$$
\kappa_\perp=\frac{a^2v^2}{3B_0^2}
\int d^3k\,\frac{P_{xx}(\mathbf k)}
{Q+\tfrac43\kappa_\perp k_\perp^2+
\dfrac{v^2k_\parallel^2}{3\kappa_\perp k_\perp^2}}.
\tag{U1}
$$

The final denominator contains $\kappa_\perp$, not $\kappa_\parallel$. For the 2D part $k_\parallel=0$ its kernel is

$$
\kappa_\perp=\frac{a^2v^2\pi}{3B_0^2}
\int_0^\infty\frac{g_2(k)}{Q+(4/3)\kappa_\perp k^2}\,dk.
\tag{U2}
$$

It can be checked with (N3) by $A=(4/3)\lambda_\parallel\lambda_\perp/(3\ell_2^2)$. Do not evaluate a $0/0$ at $k_\parallel=k_\perp=0$ as a finite slab contribution: use the distributional slab limit, whose asymptotic coefficient vanishes under the magnetostatic attachment assumptions. For pure 2D and $\lambda_\parallel\to\infty$, (U2) gives $\kappa_\perp=(a v/2)\kappa_2$, so the unmodified $a^2=1$ reproduces (F5).

### 9.2 NRMHD kernel cross-check

For the specific spectral tensor in Shalchi and Hussein [R09],

$$
P_{xx}(\mathbf k)=\frac{k_y^2 A(k_\perp)}{4\pi k_c}
\mathbf1_{|k_\parallel|\le k_c},\quad
A(k)=\frac{A_0}{[1+(k\ell_\perp)^2]^{7/3}},\quad
A_0=\frac89\ell_\perp^4\delta B^2.
\tag{U3}
$$

$k_c$ is a parallel cutoff wavenumber, not the Kubo number. The tensor normalizes to $\int P_{xx}d^3k=\delta B^2/2$. The parallel integral is analytic:

$$
\kappa_\perp=\frac{a^2v^2}{3B_0^2}\frac1{2k_c}
\int_0^\infty k^3 A(k)
\frac{\arctan[k_c\sqrt{V/U}]}{\sqrt{UV}}\,dk,
\tag{U4}
$$

with $U=Q+\kappa_\perp k^2$, $V=\kappa_\parallel$ for NLGC and $U=Q+(4/3)\kappa_\perp k^2$, $V=v^2/(3\kappa_\perp k^2)$ for UNLT. This is a spectral reduction, useful for comparing an independent direct two-dimensional quadrature with the reduced solver. Use an analytic $\arctan x/x\to1$ limit where necessary. Do not insert the 2D-only FLPD kernel into NRMHD and retain the same physical label.

### 9.3 Exact implicit-slab kernel

For the model of Shalchi (2016) [R12], use the **slab** field-line coefficient $\kappa_s$ in the compound contribution. It is not the composite coefficient (F4). The exact kernel is

$$
\kappa_\perp=\frac{v^2\pi}{3B_0^2}
\int_0^\infty\frac{g_2(k)}{Q+\kappa_\perp k^2}\,\mathcal K(\xi_k)\,dk,
\tag{U5}
$$

$$
\mathcal K(\xi)=1-\sqrt\pi\,\xi\,\operatorname{erfcx}(\xi),\qquad
\xi_k=\frac{\kappa_s\lambda_\parallel k^2}{\sqrt{3\pi}\sqrt{1+\lambda_\parallel\lambda_\perp k^2/3}}.
\tag{U6}
$$

The source MSD ansatz is $M_x=4\kappa_s\sqrt{\kappa_\parallel t/\pi}+2\kappa_\perp t$. This is why its field-line contribution is slab-specific. The limits are

$$
\mathcal K(0)=1,\qquad
\mathcal K(\xi)=\frac1{2\xi^2}-\frac3{4\xi^4}+\frac{15}{8\xi^6}+\cdots
\quad(\xi\to\infty).
\tag{U7}
$$

Use a scaled complementary error function and, at large arguments, a stable asymptotic or an equivalent positive-integral evaluation. Direct subtraction in (U6) can lose precision. $\mathcal K(\xi)$ contains $\sqrt\pi\,\xi$, not $\sqrt{\pi\xi}$.

The exact positive-integral identities used by the reference check are

$$
\mathcal K(\xi)=2\int_0^\infty u\,e^{-u^2-2\xi u}\,du
=\frac1{2\xi^2}\int_0^\infty t\,e^{-t-t^2/(4\xi^2)}\,dt\quad(\xi>0),
\tag{U12}
$$

where the first identity also applies at $\xi=0$. They follow by integration by parts and the substitution $t=2\xi u$; they change the evaluation method, not the kernel.

The rational approximation is a distinct backend,

$$
\mathcal K_{\rm rat}(\xi)=\frac1{1+2\xi^2},
\tag{U8}
$$

which gives

$$
\lambda_\perp=2D(s,q)\lambda_\parallel\frac{\delta B_2^2}{B_0^2}
\int_0^\infty
\frac{h_{s,q}(x)}{1+A x^2+Gx^4}\,dx,
\quad
A=\frac{\lambda_\parallel\lambda_\perp}{3\ell_2^2},\quad
G=\frac{2\lambda_\parallel^2\kappa_s^2}{3\pi\ell_2^4}.
\tag{U9}
$$

It is not an exact reformulation of (U5). Its error is tested separately in Section 17. Missing slab power sets $\kappa_s=0$ and (U5) reduces to its 2D NLGC-like integral with $a^2=1$.

### 9.4 Finite gyroradius: specified perturbative diagnostic

The pure-2D weak-scattering, suppressed-velocity-diffusion approximation of Shalchi (2015) [R13] is

$$
\lambda_\perp\simeq\frac32\kappa_2
\left[1-\frac{q-1}{8(s-1)}\left(\frac{r_L}{\ell_2}\right)^2\right].
\tag{U10}
$$

Use the registry identity `unlt_fgr_perturbative_diagnostic`. Require $q>1$, $s>1$ and a small correction. Positivity of the bracket is necessary but is not the validity condition of a perturbation expansion. It MUST not be extended to arbitrary $r_L/\ell_2$ or clipped to zero to hide breakdown.

The fuller pitch-angle equation employs

$$
D_\perp(\mu)=\frac{v^2\mu^2}{B_0^2}
\sum_{n=-\infty}^{\infty}\int d^3k\,P_{xx}(\mathbf k)J_n^2(W)
\frac{D_\perp k_\perp^2}
{(D_\perp k_\perp^2)^2+(v\mu k_\parallel+n\Omega)^2},\quad
W=r_L k_\perp\sqrt{1-\mu^2}.
\tag{U11}
$$

Equation (U11) is retained as the source formula, not as a complete backend contract: harmonics, singular limits, permitted spectra and suppressed-velocity-diffusion assumptions require a dedicated source and numerical audit before `unlt_fgr` is made available. Do not combine (U10) with FLPD as an automatic correction; that would define an additional model.

## 10. Field-line–particle decorrelation theory

### 10.1 Full 2D FLPD kernel

For the 2D FLPD expression in Shalchi (2021) [R02],

$$
\kappa_\perp=\frac{\pi v^2}{3B_0^2}\int_0^\infty
\frac{g_2(k)\,dk}
{Q+\kappa_\perp k^2+
\dfrac{v\kappa_{\rm FL}}{\sqrt{3\kappa_\perp}}k\sqrt{Q+\kappa_\perp k^2}}.
\tag{D1}
$$

There is no NLGC $a^2$ multiplier. The perpendicular part is the 2D spectrum; the source's composite application uses (F4) for $\kappa_{\rm FL}$ so that slab turbulence enters through field-line wandering and the supplied parallel coefficient. The denominator's $\kappa_\perp k^2$ factor is one, not UNLT's $4/3$.

For $x=k\ell_2$, $\eta_\kappa=\kappa_\perp/\kappa_\parallel=\lambda_\perp/\lambda_\parallel$, set

$$
A=\frac{\lambda_\parallel\lambda_\perp}{3\ell_2^2},\qquad
c_{\rm FL}=\frac{\kappa_{\rm FL}}{\ell_2}.
$$

Then the exact dimensionless reduction is

$$
\eta_\kappa=4D(s,q)\frac{\delta B_{x,2}^2}{B_0^2}
\int_0^\infty\frac{h_{s,q}(x)\,dx}
{1+Ax^2+c_{\rm FL}\eta_\kappa^{-1/2}x\sqrt{1+Ax^2}}.
\tag{D2}
$$

Use the **2D component** $\delta B_{x,2}^2=\delta B_2^2/2$ in the prefactor. Evaluate (D1) and (D2) independently at reference states. If 2D power vanishes, this specified asymptotic expression returns zero; it is not a description of the full running slab displacement. For the uncut spectrum, (F4) requires $q>1$.

### 10.2 General short-parallel-path limit

For $\lambda_\parallel/\ell_2\to0$ at a finite ratio, $A\to0$ and

$$
\eta_\kappa=4D(s,q)\frac{\delta B_{x,2}^2}{B_0^2}
\int_0^\infty\frac{h_{s,q}(x)}{1+\alpha_{\rm FLPD}x}\,dx,\qquad
\alpha_{\rm FLPD}=c_{\rm FL}/\sqrt{\eta_\kappa}.
\tag{D3}
$$

This remains implicit. Short $\lambda_\parallel$ alone does not imply the CLRR approximation. The fluid branch further requires $\alpha_{\rm FLPD}\ll1$ and gives $\eta_\kappa\simeq V_2$, where $V_2=\delta B_{x,2}^2/B_0^2$. In the opposite limit $\alpha_{\rm FLPD}\gg1$, let $I_{-1}=\int g_2(k)/k\,dk$; (D3) gives

$$
\eta_\kappa\simeq\left[\frac{\pi I_{-1}}{B_0^2\kappa_{\rm FL}}\right]^2.
\tag{D4}
$$

Only when $\kappa_{\rm FL}=\kappa_2=L_U\delta B_{x,2}/B_0$ can (D4) be written

$$
\eta_\kappa\simeq\frac{L_\perp^2}{4L_U^2}\frac{\delta B_{x,2}^2}{B_0^2}.
\tag{D5}
$$

This uses $\pi I_{-1}=\delta B_{x,2}^2L_\perp/2$ for (S5). A composite coefficient MUST retain (D4). Evaluate the branch parameter at each candidate asymptotic solution:

$$
\alpha_{\rm fluid}=\frac{\kappa_{\rm FL}B_0}{\ell_2\delta B_{x,2}},\qquad
\alpha_{\rm D4}=\frac{\kappa_{\rm FL}^2B_0^2}{\pi\ell_2 I_{-1}}.
\qquad
\begin{aligned}
\alpha_{\rm fluid}^{(2)}&=\frac{L_U}{\ell_2}
=\sqrt{\frac{s-1}{q-1}},\\
\alpha_{\rm D4}^{(2)}&=\frac{2L_U^2}{\ell_2L_\perp}
=\frac{\Gamma((q-1)/2)\Gamma((s+1)/2)}
{\Gamma(q/2)\Gamma(s/2)}.
\end{aligned}
\tag{D9}
$$

The superscript $(2)$ denotes pure 2D turbulence with positive power, $q>1$, and (F3). These dimensionless numbers are independent of fluctuation amplitude for that closure. They are consistency diagnostics, not supplied accuracy thresholds. The finite values in Section 17 show why neither extreme branch can be assumed solely from short $\lambda_\parallel$. A requested asymptotic approximation MUST report its evaluated $\alpha_{\rm FLPD}$ and the selected branch; automatic switching requires an explicitly chosen numerical accuracy criterion.

For pure 2D turbulence, set $\eta_\kappa=\zeta V_2$. Then (D3) reduces to

$$
\zeta=4D(s,q)\int_0^\infty
\frac{h_{s,q}(x)}
{1+(L_U/\ell_2)\zeta^{-1/2}x}\,dx.
\tag{D10}
$$

Solve this scalar equation for $\zeta>0$; no physical calibration parameter has been added. Equating the coefficient in (D5) with the algebraic form of the heuristic (F7) would require

$$
\ell_{\rm tc,match}=\frac{2L_U^2}{L_\perp}.
\tag{D11}
$$

Equation (D11) is an algebraic matching length for these two expressions, not a universal definition of transverse complexity and not permission to assign that length in another closure.

#### 10.2.1 Nonanalytic approach to the short-path limit

For fixed $\alpha>0$ define

$$
J(A,\alpha)=\int_0^\infty
\frac{h_{s,q}(x)}
{1+Ax^2+\alpha x\sqrt{1+Ax^2}}\,dx,\qquad
\left.\partial_A\frac{h_{s,q}(x)}
{1+Ax^2+\alpha x\sqrt{1+Ax^2}}\right|_{A=0}
=-\frac{h_{s,q}(x)(x^2+\alpha x^3/2)}{(1+\alpha x)^2}.
\tag{D12}
$$

For the uncut (S5) spectrum, the derivative integrand behaves as $-x^{1-s}/(2\alpha)$ at large $x$. Its integral diverges for $1<s\le2$; therefore a regular first-order expansion in $A$ is unavailable. For $1<s<2$, a tail rescaling $x=y/\sqrt A$ gives

$$
J(0,\alpha)-J(A,\alpha)\sim
\frac{A^{s/2}}{\alpha}
\int_0^\infty y^{-s-1}\left[1-\frac1{\sqrt{1+y^2}}\right]dy.
\tag{D13}
$$

Both endpoints of this integral converge in the stated range. In particular the leading correction for the uncut Kolmogorov spectrum is proportional to $A^{5/6}$, not $A$. The self-consistent ratio inherits this exponent: the derivative of the right side of (D3) with respect to $\eta_\kappa$ is strictly between zero and $1/2$ at a positive nondegenerate root. This is an algebraic convergence property of the specified kernel. A supplied finite physical spectral cutoff changes the expansion; a numerical quadrature endpoint MUST NOT be treated as that physical cutoff.

### 10.3 General long-parallel-path limit

Taking $Q\to0$ in (D1), retaining the chosen $\kappa_{\rm FL}$, gives

$$
\kappa_\perp^2+\frac{v\kappa_{\rm FL}}{\sqrt3}\kappa_\perp
-\frac{v^2}{3}\kappa_2^2=0,
\tag{D6}
$$

$$
\kappa_\perp\longrightarrow
\frac{v}{2\sqrt3}\left[\sqrt{\kappa_{\rm FL}^2+4\kappa_2^2}-\kappa_{\rm FL}\right].
\tag{D7}
$$

For numerical stability use the equivalent expression
$2v\kappa_2^2/[\sqrt3(\sqrt{\kappa_{\rm FL}^2+4\kappa_2^2}+\kappa_{\rm FL})]$ when $\kappa_2$ is small and $\kappa_{\rm FL}>0$. Return zero at zero turbulence rather than evaluating $0/0$.

Only for $\kappa_{\rm FL}=\kappa_2$ does this become

$$
\kappa_\perp\longrightarrow
\frac{\sqrt5-1}{\sqrt3}\frac v2\kappa_{\rm FL}.
\tag{D8}
$$

The factor in (D8) is not universal for a composite field-line coefficient. These limits are consequences of the specified kernel, not newly proposed transport laws. The inspected 2021 arXiv source prints $\sqrt3$ in its quadratic; the earlier description's claim of a $\sqrt2$ typo is withdrawn.

## 11. Random ballistic decorrelation and composite closed form

### 11.1 RBD with backtracking correction

For the RBD/BC model restated and tested in Snodin et al. [R14],

$$
\kappa_\perp=\frac{a^2v^2}{3B_0^2}\int d^3k\,P_{xx}(\mathbf k)
\sqrt{\frac\pi2}\frac{\operatorname{erfc}(\beta_{\mathbf k})}{\sqrt{V_{\mathbf k}}},\qquad
\beta_{\mathbf k}=\frac Q{\sqrt{2V_{\mathbf k}}},\quad
V_{\mathbf k}=\sum_i k_i^2\langle\widetilde v_i^2\rangle.
\tag{B1}
$$

The ordinary Gaussian ballistic time integral and the backtracking-corrected one differ [R28, Equations (2.6)–(2.10); R14, Equations (22)–(23)]:

$$
T_{\rm plain}(Q,V)=\int_0^\infty e^{-Qt-Vt^2/2}dt
=\sqrt{\frac{\pi}{2V}}e^{\beta^2}\operatorname{erfc}(\beta)
=\sqrt{\frac{\pi}{2V}}\operatorname{erfcx}(\beta),\qquad
T_{\rm BC}(Q,V)=e^{-\beta^2}T_{\rm plain}
=\sqrt{\frac{\pi}{2V}}\operatorname{erfc}(\beta),
\quad\beta=\frac Q{\sqrt{2V}}.
\tag{B9}
$$

Thus (B1) deliberately uses the corrected kernel; replacing erfc by erfcx in the `rbd_bc` backend removes the source's backtracking correction. An uncorrected diagnostic MUST use a distinct identifier. The supplied magnetostatic expression has the source's temporal-decorrelation rate set to zero.

For its axisymmetric transverse fluctuation model,

$$
\langle\widetilde v_x^2\rangle=\langle\widetilde v_y^2\rangle
=\frac{a^2v^2}{6}\frac{\delta B^2}{B_0^2},\qquad
\langle\widetilde v_z^2\rangle=\frac{v^2}{3}
\left(1-a^2\frac{\delta B^2}{B_0^2}\right).
\tag{B2}
$$

The variance in (B2) is the full transverse variance belonging to that closure. Require nonnegative variances. In particular $a^2\delta B^2/B_0^2\le1$; with $a^2=1/3$ this is $\delta B/B_0\le\sqrt3$. Handle zero-variance endpoints by the actual kernel limit or reject an unsupported endpoint explicitly. This restriction arises from (B2), not the ordinary NLGC denominator.

RBD is an explicit quadrature once $\kappa_\parallel$ is supplied. Do not allocate a $\kappa_\perp$ root solve. If $Q>0$, the corrected kernel tends to zero as $V_{\mathbf k}\to0$; the plain kernel tends to $1/Q$. Evaluate the selected limit stably. With zero turbulence return zero directly. Infinite $\lambda_\parallel$ and degenerate velocity variances require separate convergence analysis.

For a 2D area spectrum and $k_\parallel=0$, angular integration gives

$$
\kappa_\perp=\frac{a^2v^2}{6B_0^2}\sqrt{\frac\pi2}
\int_0^\infty
\frac{2\pi kS_2(k)\operatorname{erfc}[Q/(k\sqrt{2V_x})]}
{k\sqrt{V_x}}\,dk,\qquad V_x=\langle\widetilde v_x^2\rangle.
\tag{B3}
$$

For **pure 2D** (S12), $a^2=1/3$, $s_{\rm A}=5/3$, $p=2$, its long-path limit is

$$
\frac{\kappa_\perp}{v\lambda_2}
\longrightarrow
\frac{\sqrt\pi}{6}(2\pi C_2)
\left[\frac1{p+1}+\frac{1}{s_{\rm A}}\right]\frac{\delta B_2}{B_0}.
\tag{B4}
$$

For precisely these indices $C_2=2/(7\pi)$ and the dimensionless prefactor in (B4) is $4\sqrt\pi/45$. These are exact reductions of (S12) and (B4), not fitted numbers.

The pure-2D condition matters: if the velocity variance also includes slab power, retain the total variance in $V_x$ when evaluating (B3). Comparison of a closure with a published particle simulation assesses its physical error; it cannot set the tolerance for numerical quadrature.

### 11.2 Composite closed form

The composite approximation MUST take a positive perpendicular length $\ell$ with a tagged definition [R18,R14]. Its source profiles are distinct:

| Length profile | Required meaning of $\ell$ | Reproduction scope |
|---|---|---|
| `shalchi_2019_bend_over` | Perpendicular bend-over scale of the stated turbulence model | Original heuristic source [R18] |
| `snodin_2022_integral` | Source's perpendicular integral length below | Application tested in [R14] |
| `parameterized_tagged_length` | Supplied positive length and its explicit definition | Parameterized use; no exact source reproduction claim |

Snodin et al. define, for their NRMHD area amplitude $A_{\rm NRMHD}(k)$,

$$
\ell_{\perp,\rm int}=
\frac{\int_0^\infty k^2A_{\rm NRMHD}(k)dk}
{\int_0^\infty k^3A_{\rm NRMHD}(k)dk}
=\frac{I_{-1}}{I_0}=\frac{L_\perp}{2},
\quad I_j=\int_0^\infty k^jg_2(k)dk.
\tag{B8}
$$

The second equality applies when the reduced tensor normalization maps $g_2\propto k^3A_{\rm NRMHD}$; the last follows from (S9). This factor of two distinguishes that source's integral length from this document's $L_\perp$. For (S5), compute both moments with the supplied spectrum and retain the length-profile tag. Do not relabel a bend-over scale as an integral length.

With $\kappa_{\rm FL}$ supplied from the declared field-line closure, the relation is

$$
\kappa_{\rm FL}=\frac23\lambda_\perp+
\ell\sqrt{\frac{\lambda_\perp}{\lambda_\parallel}}.
\tag{B5}
$$

For positive $\lambda_\parallel$ this gives

$$
\lambda_\perp=\frac{9\ell^2}{16\lambda_\parallel}
\left[\sqrt{1+\frac{8\kappa_{\rm FL}\lambda_\parallel}{3\ell^2}}-1\right]^2
=\frac{4\kappa_{\rm FL}^2\lambda_\parallel/\ell^2}
[\sqrt{1+8\kappa_{\rm FL}\lambda_\parallel/(3\ell^2)}+1]^2.
\tag{B6}
$$

Use the rationalized form, including its continuous value zero at $\lambda_\parallel=0$. Its limits are

$$
\lambda_\perp\simeq\frac{\kappa_{\rm FL}^2}{\ell^2}\lambda_\parallel
\quad\left(\frac{\kappa_{\rm FL}\lambda_\parallel}{\ell^2}\ll1\right),\qquad
\lambda_\perp\to\frac32\kappa_{\rm FL}\quad(\lambda_\parallel\to\infty).
\tag{B7}
$$

This is the closed-form approximation, not a fluid-limit implementation of (F7). Comparing (B7) with (D8) gives their factor difference only in the pure/quasi-2D limit with the same field-line coefficient; for composite FLPD use (D7). Missing length metadata returns `missing_input`.

## 12. Compound transport, retracing and pre-diffusive diagnostics

### 12.1 General compound moment and flat-2D GCD

If particles remain on static field lines and their parallel displacement probability density is $p_\parallel(z,t)$, the compound moment assumption is [R15]

$$
M_x(t)=\int_{-\infty}^{\infty}M_{\rm FL}(z)\,p_\parallel(z,t)\,dz.
\tag{C1}
$$

It assumes that the conditional field-line wandering entering the integral is the supplied $M_{\rm FL}$. Correlations not represented by this closure cannot be inferred from an MSD. For a zero-mean Gaussian parallel density of variance $M_z(t)>0$,

$$
p_\parallel(z,t)=\frac{e^{-z^2/(2M_z)}}{\sqrt{2\pi M_z}},\qquad
\langle|z|^a\rangle=\frac{2^{a/2}\Gamma((a+1)/2)}{\sqrt\pi}M_z^{a/2},
\quad a>-1.
\tag{C2}
$$

Thus $M_{\rm FL}=2\kappa_{\rm FL}|z|$ gives (F8) when $M_z=2\kappa_\parallel t$.

For the **flat** energy-range 2D spectrum ($q=0$), the large-$|z|$ field-line result used by Shalchi and Kourakis [R03,R15] is

$$
M_{\rm FL}(z)=
\left[9C(s)\sqrt{\frac\pi2}\right]^{2/3}
\left(\frac{\delta B_2}{B_0}\right)^{4/3}
\ell_2^{2/3}|z|^{4/3}.
\tag{C3}
$$

Combining (C2) and (C3) gives the complete, dimensionally consistent expression

$$
M_x(t)=\alpha_G(s)
\left(\frac{\delta B_2}{B_0}\right)^{4/3}
\ell_2^{2/3}M_z(t)^{2/3},\qquad
\alpha_G(s)=\frac{\Gamma(7/6)}{\sqrt\pi}
\left[18C(s)\sqrt{\frac\pi2}\right]^{2/3}.
\tag{C4}
$$

This coefficient is the authors' explicit conference Equation (6), checked visually against the primary PDF [R15]. Its independent evaluation at $s=5/3$ is **1.01024378252518**, not the nearby source sentence's approximate value 0.5. The latter is inconsistent with both that printed formula and the direct Gaussian integral. This specification uses (C4), records the discrepancy and does not claim a new physical prefactor. The journal full text could not be retrieved directly; exact reproduction of a journal-specific alternative remains a source gate. The older guide's numerical estimate based on the rounded 0.5 is not a normative fixture.

The length in (C4) is $\ell_2^{2/3}$, not an isolated $\ell_2$. For diffusive parallel motion $M_z=2\kappa_\parallel t$, let $A_G=\alpha_G(\delta B_2/B_0)^{4/3}\ell_2^{2/3}(2\kappa_\parallel)^{2/3}$. Then

$$
M_x=A_G t^{2/3},\qquad
d_x=\frac{A_G}{3}t^{-1/3},\qquad
\bar d_x=\frac{A_G}{2}t^{-1/3},\qquad \kappa_\perp=0.
\tag{C5}
$$

The scattering-time average discussed by the authors uses their secant estimator and $t_c=\lambda_\parallel/v$. With $\lambda_\parallel=3\kappa_\parallel/v$, this equals the nominal $\tau_s=3\kappa_\parallel/v^2$ in Section 13; it is distinct from that section's fitted field-line decorrelation time $\tau_c$:

$$
\overline\kappa_{\perp,c}=\frac1{t_c}\int_0^{t_c}\bar d_x(t)\,dt,\qquad
\overline\lambda_{\perp,c}=\frac3v\overline\kappa_{\perp,c}
=\left(\frac32\right)^{4/3}\alpha_G(s)
\left(\frac{\delta B_2}{B_0}\right)^{4/3}
\ell_2^{2/3}\lambda_\parallel^{1/3}.
\tag{C6}
$$

Equation (C6) is an average of an extrapolated diffusive-parallel approximation, not a nonzero asymptotic perpendicular path and not a claim of correct ballistic behavior at $t=0$. The averaged derivative estimator is instead

$$
\frac1{t_c}\int_0^{t_c}d_x(t)dt
=\frac{M_x(t_c)}{2t_c}
=\frac{A_G}{2}t_c^{-1/3}
=\frac23\,\overline\kappa_{\perp,c}.
\tag{C12}
$$

Both averages use the same extrapolated law. The factor $2/3$ follows from their estimator definitions and MUST NOT be absorbed into a turbulence normalization. A complete finite-time particle model MUST supply the appropriate $p_\parallel$ and short-distance field-line behavior in (C1).

The flat-2D GCD asymptotic is distinct from the finite-$q>1$ field-line coefficient in (F3). Do not use (F3) to regularize (C3), or apply (C4) to a $q=3$ spectrum without changing the model definition.

### 12.2 Persistent frozen-line algorithm

The field-line meandering method of Laitinen et al. [R16] generates a line before propagating the particle. For a coefficient defined per mean-field coordinate $z$, as in (F2)–(F4), a line diffusion increment over a newly explored monotonic segment is

$$
\Delta X=\sqrt{2\kappa_{\rm FL}\Delta z}\,\xi_x,\quad
\Delta Y=\sqrt{2\kappa_{\rm FL}\Delta z}\,\xi_y,\quad
\xi_x,\xi_y\sim N(0,1),\quad\Delta z>0.
\tag{C7}
$$

Store a line realization $\mathbf X_{\rm FL}(z)$ and its persistent identifier. A particle at signed position $z(t)$ samples that same realization; it MUST retrace the same geometry when it reverses. Incremental generation is permissible on first exploration of a segment. Previously visited segments are reused. A nested Brownian bridge or a sufficiently resolved stored curve can make line refinement reproducible. Random values keyed only by the current particle timestep do not define a frozen line.

For a constant-coefficient Brownian line, the bridge at $z_a<z<z_b$, conditioned on stored endpoints, is

$$
X(z)\mid X_a,X_b\sim N\left(
X_a+\frac{z-z_a}{z_b-z_a}(X_b-X_a),\
2\kappa_{\rm FL}\frac{(z-z_a)(z_b-z)}{z_b-z_a}
\right).
\tag{C8}
$$

This is a numerical construction of the chosen line process, not added particle physics. For variable $\kappa_{\rm FL}$ along a prescribed line coordinate, use the accumulated variance coordinate $u(z)=2\int\kappa_{\rm FL}(z')dz'$ on each monotonic branch; if the coefficient depends on the realized transverse path, solve the declared field-line stochastic model instead of assuming this scalar time change.

Applying fresh independent $\sqrt{2\kappa_{\rm FL}|\Delta z|}$ kicks at **every particle step** is wrong for compound retracing. For Gaussian parallel increments of duration $\Delta t$,

$$
\mathbb E\sum_{n=1}^{t/\Delta t}|\Delta z_n|
=2\sqrt{\frac{\kappa_\parallel}{\pi}}\frac t{\sqrt{\Delta t}},\qquad
M_x^{\rm fresh}(t)=4\kappa_{\rm FL}\sqrt{\frac{\kappa_\parallel}{\pi}}
\frac t{\sqrt{\Delta t}}.
\tag{C9}
$$

It diverges on timestep refinement and fails (F8). Required tests include exact reuse on an out-and-back path, compound $t^{1/2}$ MSD growth, timestep refinement, line-resolution refinement, checkpoint/restart and ordering-independent line identifiers. The transport configuration MUST explicitly select either an independent line realization per particle or a shared realization keyed by a supplied line identifier. Particles assigned the same line identifier see the same stored geometry. Both policies can have the same single-particle marginal (C1); they produce different interparticle correlations and statistical uncertainties. Physical sharing or connectivity MUST NOT be inferred from the MSD formula.

The coordinate tag is mandatory. A Brownian coefficient per mean-field coordinate $z$ cannot be applied per physical arc length $\sigma$ without a supplied model and conversion. That conversion is not determined by this specification. Persistent line identifiers, refinement, random-number ownership and restart serialization belong to the transport layer, not the deterministic coefficient library.

If independent particle diffusion relative to the line is added, its coefficient MUST come from a model explicitly describing that relative motion. Adding a closure that already includes the same field-line wandering can double count transport. Injection spread belongs to the source, not to either noise mechanism.

### 12.3 Observable separation

Return types MUST distinguish: an unconditional particle MSD; a centered particle covariance; a running derivative; an MSD-ratio estimator; a field-line MSD in space; a displacement relative to the original line; and a conditional return-plane statistic. None is substituted for another by changing a unit label.

### 12.4 Laitinen–Dalla pre-diffusive fit

The source [R17] studies protons in a **uniform mean field** plus slab/2D turbulence. It does not use a Parker spiral for this diagnostic. It fits the **median of squared gyrocentre displacement in the injection plane for returning particles**. The return-plane construction suppresses direct meandering displacement; the statistic is conditional and is not the ensemble transverse MSD.

The phenomenological fit is

$$
m_{\rm ret}(t)=A_1\left(\frac t{T_L}\right)^{\alpha_m}
\frac{1+(t/t_1)^{\beta_m-\alpha_m}}
{1+(t/t_2)^{\beta_m-1}},\qquad T_L=2\pi/\Omega.
\tag{C10}
$$

Supply $A_1$ in $\mathrm{m^2}$; $t_1,t_2,T_L$ in seconds; and dimensionless exponents. For the intended ordered transitions require $0<t_1<t_2$, $\beta_m>\alpha_m$ and $\beta_m>1$. $\beta_m$ is the intermediate transition exponent, not an undefined normalization. $A_1$ is the early-regime amplitude at the gyroperiod, not exactly $m_{\rm ret}(T_L)$ unless both transition factors are negligible.

The asymptotic amplitude follows directly from (C10):

$$
m_{\rm ret}(t)\sim A_2\frac t{T_L},\qquad
A_2=A_1\frac{t_2^{\beta_m-1}}
{T_L^{\alpha_m-1}t_1^{\beta_m-\alpha_m}}.
\tag{C11}
$$

This resolves the earlier undefined $t_0$ by algebra from source Equation (6); it is not an inferred physical timescale. The source describes the fit as a phenomenological descriptor. Do not convert its median slope to $\kappa_\perp$, use it as Markov noise, or identify it with an injection width. There are no default fitted parameters. Any fit dataset MUST identify the return rule, gyrocentre definition, species, field, energy, sample selection and fit window.

## 13. Isotropic-turbulence fits and MHD scalings

### 13.1 Distinct amplitude, radius and length conventions

For an isotropic random field plus a uniform ordered field, define

$$
\sigma^2=\frac{\delta B^2}{B_0^2},\qquad
\eta_B=\frac{\delta B^2}{B_0^2+\delta B^2},\qquad
\sigma^2=\frac{\eta_B}{1-\eta_B},\qquad
B_{\rm rms}=\sqrt{B_0^2+\delta B^2}.
\tag{I1}
$$

$\eta_B$ is a turbulence-energy fraction; $\eta_\kappa$ is a diffusion-coefficient ratio. They MUST not share an implementation variable. Candia–Roulet uses an ordered-field radius; Snodin uses a total-rms-field radius. Fit inputs SHOULD carry the exact radius definition rather than a generic `rigidity` number.

For a finite-range power spectrum with index $g>1$ between $L_{\min}$ and $L_{\max}$, the Candia–Roulet correlation convention is [R19]

$$
L_c=\frac{L_{\max}}2\frac{g-1}{g}
\frac{1-(L_{\min}/L_{\max})^g}
{1-(L_{\min}/L_{\max})^{g-1}}.
\tag{I2}
$$

It approaches $L_{\max}/5$ for $g=5/3$ and a wide inertial range. Snodin instead uses

$$
l_c=\frac\pi2\frac{\int k^{-1}M(k)dk}{\int M(k)dk},\qquad
l_c=\frac{\pi(s-1)}{2s k_0},\quad k_0=2\pi/L,
\tag{I3}
$$

for a sharp low-$k$ cutoff and $M(k)\propto k^{-s}$ above it. At $s=5/3$, $l_c=L/10$. These correlation conventions differ; a factor-of-two conversion is justified only if the two underlying outer scales and spectral shapes actually correspond.

### 13.2 Candia–Roulet pair

For ultrarelativistic particles in the source setup [R19], with $\rho_C=r_L(B_0)/L_{\max}$,

$$
\kappa_\parallel=cL_{\max}\rho_C\frac{N_\parallel}{\sigma^2}
\sqrt{\left(\frac{\rho_C}{\rho_\parallel}\right)^{2(1-g)}
+\left(\frac{\rho_C}{\rho_\parallel}\right)^2},
\tag{I4}
$$

$$
\frac{\kappa_\perp}{\kappa_\parallel}
=N_\perp(\sigma^2)^{a_\perp}
\begin{cases}1,&\rho_C\le0.2,\\
(\rho_C/0.2)^{-2},&\rho_C>0.2.
\end{cases}
\tag{I5}
$$

The source Table 1 parameters are:

| Spectrum $g$ | $N_\parallel$ | $\rho_\parallel$ | $N_\perp$ | $a_\perp$ | $N_A$ |
|---|---:|---:|---:|---:|---:|
| Kraichnan, $3/2$ | 2.0 | 0.22 | 0.019 | 1.37 | 17.6 |
| Kolmogorov, $5/3$ | 1.7 | 0.20 | 0.025 | 1.36 | 14.9 |
| Bykov–Toptygin, $2$ | 1.4 | 0.16 | 0.020 | 1.38 | 14.2 |

The source states the parallel fit's formal rigidity range $0.01\le\rho_C\le1$ and good behavior up to approximately $\sigma^2=10$. The perpendicular fit was made with $\rho_C\ge0.03$ and $\sigma^2\ge1$ to avoid sub-diffusion. Thus the conservative joint-fit domain used for an initial implementation is $0.03\le\rho_C\le1$, $1\le\sigma^2\le10$; this is a declared software domain based on the common checked ranges, not a universal physical boundary. No switch to the paper's high-turbulence parallel parameters is implied. Outside-domain formula evaluations are diagnostic extrapolations with explicit flags.

The $c$ in (I4) belongs to the relativistic fit. Replacing it by $v$ for low-energy SEPs is an additional prescription and MUST have a different model identifier and justification. The source returns a pair; replacing its parallel component changes the perpendicular result and the calibrated model.

### 13.3 Snodin basic fit and spectral-shape scope

The full basic pair below is specified for the sharp-cutoff Kolmogorov spectrum in [R22]. Use $x_S=R_L(B_{\rm rms})/L$, $L=2\pi/k_0$, and the source particle speed $v_0$:

$$
\kappa_0=v_0L(a_1+a_2x_S),\qquad
\kappa_\parallel=v_0L\left[a_1+a_2x_S+
\frac13x_S^{1/3}\frac{1-\eta_B}{\eta_B}\right].
\tag{I6}
$$

$$
\kappa_\perp=\frac{\kappa_0}{1+\chi(1-\eta_B)/\eta_B}.
\tag{I7}
$$

The source supplies these finite-range **isotropic $\kappa_0$ fits**:

| Spectral shape for $\kappa_0$ | $a_1$ | $a_2$ | Source definition |
|---|---:|---:|---|
| Sharp-cutoff Kolmogorov | 0.0031 | 0.74 | $M(k)\propto k^{-5/3}$ above $k_0$ |
| Sharp-cutoff Kraichnan | 0.0019 | 0.76 | $M(k)\propto k^{-3/2}$ above $k_0$ |
| Peaked Kolmogorov | 0.0017 | 0.75 | Source Equation (8), $s=5/3$, $k_b=5k_0$ |

For the peaked shape,
$M(k)=M_0(k/k_b)^4/[1+(k/k_b)^2]^{s/2+2}$.
The source's $k_b=5k_0$ does not redefine $L$: the fitted radius remains $R_L/L$ with $L=2\pi/k_0$. A reproduction that generates the turbulent field requires its actual wavenumber support and normalization at the specified variance. $M_0$ is not an additional independent calibration amplitude once that variance has been fixed.

For the Kolmogorov pair, the quoted best-fit $\chi\simeq2.35$ is approximate. The source relates it to $4D_{\rm FL,iso}/l_c$, where $D_{\rm FL,iso}$ is its field-line diffusion coefficient at $B_0=0$ (source Equation (22)), with length units. This is distinct from $D(s,q)$, and is not the same field-line closure as (F2)–(F4). The named pair uses 2.35 as its stated approximate fit parameter.

Registry `iso_fit_snodin_2016` MUST restrict (I6)'s ordered-field addition and the selected (I7) calibration to the sharp-cutoff Kolmogorov pair. The other rows supply $\kappa_0$ at $\eta_B=1$ through `iso_kappa0_snodin_2016`. They do not, by themselves, authorize changing only $a_1,a_2$ and retaining the Kolmogorov ordered-field exponent or $\chi$. Ordered-field pair backends for those shapes remain `source_gate` until their complete source definition, calibration and domain are provided.

The linear $\kappa_0$ fit was made over approximately $0.004\lesssim x_S\lesssim0.05$. For the initial Kolmogorov pair this interval and $0.1\le\eta_B\le1$ are a declared conservative software domain based on sampled comparisons, not a source guarantee of uniform accuracy throughout a filled parameter box. More precise spectral-resolution or amplitude domains require the actual fit dataset. At $\eta_B=1$, return $\mathbf K=\kappa_0\mathbf I$ with no preferred parallel axis, for each supported isotropic fit. The finite intercept belongs to these finite-range fits and is not a universal $r_L\to0$ diffusion floor.

Improved source fits with additional powers of $R_L$ are separate variants. Neither the basic pair nor $\kappa_0$ backend may infer a transition at $l_c/2$ or blend variants without a complete selected definition.

### 13.4 Kuhlen–Mertsch–Phan running model

The 2025 analytical model [R23] combines a running parallel coefficient with a running field-line fit. It is a heuristic finite-time closure, not identical to the Gaussian convolution (C1). The field-line source definition is per arc length $\sigma$ (source Equation (3)); replacing $\sigma$ by the mean-field coordinate $z$ is its stated small-turbulence approximation. These coordinates MUST retain separate tags.

Let $x_K=r_g/L_c$ with the source's correlation convention. Part I of the source series [R27] uses (I2), with $g=5/3$ for its finite-band Kolmogorov spectrum; $L_c\simeq L_{\max}/5$ is its wide-inertial-range limit. It defines $r_g$ using the total rms field $B_{\rm rms}$; in SI, $r_g=p/(|q_e|B_{\rm rms})$. Dimensional restoration of the normalized parallel fit is

$$
\kappa_\parallel=vL_c A x_K^{1/3}
\left[1+\left(\frac{x_K}{\rho_*}\right)^{5/(3s_\kappa)}\right]^{s_\kappa}.
\tag{I8}
$$

Use that total-rms-field radius for the source parameter sets below. Any externally supplied calibration dataset MUST retain its own radius and correlation-length metadata. A radius formed from $B_0$ is not interchangeable with one formed from $B_{\rm rms}$.

The running form is

$$
\tau_s=\frac{3\kappa_\parallel}{v^2},\quad
d_\parallel(t)=\kappa_\parallel
\left[1+(t/\tau_s)^{-1/s_d}\right]^{-s_d},\quad s_d=\frac12,
\quad M_z(t)=2\int_0^td_\parallel(u)du.
\tag{I9}
$$

For $s_d=1/2$, the integral has the stable exact reduction

$$
d_\parallel(t)=\frac{\kappa_\parallel t}{\sqrt{t^2+\tau_s^2}},\qquad
M_z(t)=2\kappa_\parallel[\sqrt{t^2+\tau_s^2}-\tau_s]
=\frac{2\kappa_\parallel t^2}{\sqrt{t^2+\tau_s^2}+\tau_s}.
\tag{I10}
$$

It gives $M_z\simeq v^2t^2/3$ at early times and $M_z\simeq2\kappa_\parallel t$ at late times. The field-line fit is

$$
d_{\rm FL}(z)=C_Kz\frac{\delta B_x^2}{B_0^2}
\left[1+(z/z_1)^{(1-\gamma_K)/s_1}\right]^{-s_1}
\left[1+(z/z_2)^{-\gamma_K/s_2}\right]^{s_2},\quad
s_1=1.5,\quad s_2=0.2.
\tag{I11}
$$

For isotropic fluctuations $\delta B_x^2=\delta B^2/3$. $C_K$ is a field-line fit normalization, not the spectral function $C(s)$. Retain it: the small-$z$ coefficient is $C_K\delta B_x^2/B_0^2$ for $\gamma_K<0$. With $\gamma_K=0$ the second bracket is the constant $2^{s_2}$; do not discard it silently. This is the explicitly selected no-intermediate-subdiffusion variant, not the same fitted normalization.

The source parameters checked in Table 1 and the text following Equation (24) are:

| $\eta_B$ | $A$ | $\rho_*$ | $s_\kappa$ | $C_K$ | $z_1/L_c$ | $z_2/L_c$ |
|---:|---:|---:|---:|---:|---:|---:|
| 0.2 | 3.81 | 0.755 | 0.796 | 0.8 | 1.5 | 8 |
| 0.5 | 1.06 | 0.627 | 0.906 | 0.55 | 1.5 | 5 |
| 0.8 | 0.331 | 0.428 | 1.22 | 0.2 | 2.5 | 5.5 |

The fitted $C_K$ values are retained exactly as quoted. They MUST NOT be renormalized by the slab/2D convention (F1), set to one, or replaced by $1-\eta_B$; in particular the middle value is 0.55, not 0.5. The source's isotropic random field includes fluctuating $B_z$, whereas the slab/2D geometry in (F1) has fixed longitudinal field. Its exact arc-length tangent is $\mathbf B/|\mathbf B|$; the small-turbulence approximation leading to the printed $B_0^2$ denominator in (I11) is not an identity for arbitrary isotropic amplitude. Retain the published approximate closure and state its scope.

The running perpendicular relation is

$$
z_*(t)=\sqrt{M_z(t)},\qquad
d_\perp^{\rm att}(t)=\frac{d_{\rm FL}(z_*)}{z_*}d_\parallel(t),\qquad
M_x^{\rm att}(t)=2\int_0^{z_*}d_{\rm FL}(z)dz.
\tag{I12}
$$

Use the small-$z$ ratio limit at $t=0$ and return $d_\perp(0)=M_x(0)=0$. The asymptotic result is frozen at a **decorrelation time**, not the scattering time. The source uses the approximate condition $d_\perp(\tau_c)\sim L_{c,\perp}^2/(2\tau_c)$ with a calibrated transverse scale. A specified numerical realization of that fit solves

$$
2\tau_c d_\perp^{\rm att}(\tau_c)=L_{c,\perp}^2,\qquad
\kappa_\perp=d_\perp^{\rm att}(\tau_c),\qquad
d_\perp(t)=\begin{cases}d_\perp^{\rm att}(t),&t\le\tau_c,\\
\kappa_\perp,&t>\tau_c.
\end{cases}
\tag{I13}
$$

$\gamma_K<0$ and $L_{c,\perp}>0$ are **required calibration inputs** for the intermediate-subdiffusion model. The source explains that they are fitted to perpendicular particle data; it does not supply their numerical values in the retrieved parameter table. A complete published perpendicular preset therefore remains unavailable until those data are supplied. No guess or interpolation in $\eta_B$ is provided. Replacing the approximate relation by equality in (I13) is the declared numerical fit convention; the uncertainty is absorbed into the calibrated scale, not justified by a new physical argument. Record a root-selection rule and report multiple candidate roots if they occur; use the first upward crossing only when that convention belongs to the supplied calibration.

For continuity of the MSD after freezing, use
$M_x(t)=M_x^{\rm att}(\tau_c)+2\kappa_\perp(t-\tau_c)$ for $t>\tau_c$. Choosing $d_\perp(\tau_s)$ would implement a different model. Without calibrated $\gamma_K$ and $L_{c,\perp}$, return `missing_calibration` rather than a perpendicular coefficient.

### 13.5 Classical scattering and limited scaling fits

For a single isotropic velocity-scattering time $\tau$, define $x=\Omega\tau=\lambda_\parallel/r_L$ and $\kappa_\parallel=v^2\tau/3$. The classical relation is

$$
\kappa_\perp=\frac{\kappa_\parallel}{1+x^2},\qquad
\kappa_A=\operatorname{sgn}(q_e)\frac{\kappa_\parallel x}{1+x^2},\qquad
\kappa_A^{\rm ws}=\operatorname{sgn}(q_e)\frac{vr_L}{3},\qquad
\frac{\kappa_A}{\kappa_A^{\rm ws}}=\frac{x^2}{1+x^2}.
\tag{I14}
$$

This is a specific scattering model, not a universal ceiling on $\kappa_\perp$. Turbulent field-line transport does not in general obey its ratio. At $x\gg1$, $\kappa_\perp\simeq vr_L^2/(3\lambda_\parallel)$.

Casse et al. [R24] report a scaling $\kappa_\perp/\kappa_\parallel\propto\eta_B^{2.3\pm0.2}$ for a particular regime and discuss an approximate prefactor rather than one exact, globally applicable fit. Restricted-scattering power laws likewise need amplitudes, spectrum support and fit domains. These are useful comparisons but do not specify an absolute coefficient. Keep their named production backends gated until the normalization and dataset-specific domains are supplied; a user-calibrated generic power law can use (P2) with its own provenance.

With zero mean field and statistical isotropy,

$$
K_{xx}=K_{yy}=K_{zz}=\kappa_{\rm iso},\qquad
\operatorname{tr}\mathbf K=3\kappa_{\rm iso}.
\tag{I15}
$$

There is no distinguished parallel coefficient whose one-third is a Cartesian eigenvalue. A transient or finite-realization anisotropy needs its own measured covariance tensor.

### 13.6 MHD $M_A^4$ estimate

Define $M_A=\delta v/v_A$ at the source model's injection scale. With the same-scale Alfvén ratio $r_A$, $M_A=\sqrt{r_A}\,\delta B/B_0$; equality to $\delta B/B_0$ requires equipartition. For the sub-Alfvénic regime $M_A<1$ and $\lambda_\parallel<L$, Yan and Lazarian's scaling [R25] is

$$
\kappa_\perp\sim\kappa_\parallel M_A^4.
\tag{I16}
$$

A numerical estimate requires an explicit coefficient $c_{M4}$, so a parameterized approximation uses $\kappa_\perp=c_{M4}\kappa_\parallel M_A^4$ and labels itself `order_of_magnitude`. Its coefficient is not a calibrated universal default. The source's other regimes use different scales, including $l_A=L M_A^{-3}$ for super-Alfvénic turbulence. They MUST be separate named branches with their source assumptions. Critically balanced MHD scalings are not slab/2D closures.

For field-line superdiffusive relations written with $\sim$, numerical factors obtained by replacing $\sim$ with $=$ are not determined by the physics of the scaling. The earlier claim of an exact factor discrepancy in such a relation is withdrawn. No uncalibrated $1/27$ or $1/81$ is promoted to an exact library coefficient.

## 14. GCR prescriptions and drift

### 14.1 General rigidity-law building block

Two dimensionally defined smooth broken forms suffice for the parameterized modulation families [R01,R26]:

$$
G_N(\mathcal R;\mathcal R_0,\mathcal R_k,a,b,c_s)
=\left(\frac{\mathcal R}{\mathcal R_0}\right)^a
\left[\frac{(\mathcal R/\mathcal R_0)^{c_s}+(\mathcal R_k/\mathcal R_0)^{c_s}}
{1+(\mathcal R_k/\mathcal R_0)^{c_s}}\right]^{(b-a)/c_s},
\tag{G1}
$$

$$
G_C(\mathcal R;\mathcal R_k,a,b,s_R)
=\left(\frac{\mathcal R}{\mathcal R_k}\right)^a
\left[1+(\mathcal R/\mathcal R_k)^{s_R}\right]^{(b-a)/s_R}.
\tag{G2}
$$

Positive reference rigidities and smoothing exponents are required. $G_N(\mathcal R_0)=1$; $G_C(\mathcal R_k)=2^{(b-a)/s_R}$, so their normalization constants are not interchangeable. A coefficient has the form $K_0\beta(B_n/B_0)G$, with $K_0$ in SI and its units declared before evaluation. Distinct perpendicular slopes give distinct rigidity laws; they do not remain a constant ratio to a parallel law with different slopes.

For either form, the effective logarithmic rigidity slope is $a+(b-a)y/(1+y)$, with $y=(\mathcal R/\mathcal R_k)^{c_s}$ or $(\mathcal R/\mathcal R_k)^{s_R}$ respectively, plus any $\beta$ derivative required by the chosen particle coordinate.

### 14.2 NWU ratio with polar enhancement

The common ratio structure is

$$
\kappa_{\perp r}=\eta_r\kappa_\parallel,\qquad
\kappa_{\perp\theta}=\eta_\theta\kappa_\parallel f(\theta),\qquad
\theta_A=\min(\theta,\pi-\theta),
\tag{G3}
$$

with a fully specified folded-colatitude function

$$
f(\theta)=\frac{d+1}{2}-\frac{d-1}{2}
\tanh[h_\theta(\theta_A-\pi/2+\theta_F)].
\tag{G4}
$$

The inputs $\eta_r,\eta_\theta,d,\theta_F,h_\theta$ are required; $h_\theta$ is in $\mathrm{rad^{-1}}$. Potgieter's review gives the commonly used ratios 0.02 and a proton enhancement parameter $d=3$ [R01]. Its printed sign and degree notation require an explicitly recorded angular interpretation. The software MUST not infer radians by visual plausibility or silently choose one ambiguous sign convention. A named `folded_radian` variant of (G4), with an explicitly supplied width, is a complete parameterized prescription, not certification of an undocumented historical implementation.

The function tends toward $d$ and 1 in the respective saturated branches; at finite angles it is not exactly equal to those values. At the folded equator it is continuous but generally not differentiable because $\theta_A$ has a cusp. Use one-sided derivatives there or the Cartesian derivative strategy of Section 15. Do not claim global $C^1$ regularity. A smoothing alteration is a separately identified numerical/model variant, not a correction that may be inserted silently.

If a distinct perpendicular rigidity law is selected, define $\kappa_{\perp\theta}=f(\theta)$ times the **declared perpendicular normalization/law**, rather than simultaneously asserting equality to two different parallel and perpendicular laws. A source preset MUST identify the species, epoch, equation, angular convention and exact fitted parameter set.

### 14.3 Corti AMS-02 parameterization

The Corti family [R26] uses (G2) with separate parallel and perpendicular low/high rigidity slopes. The radial and polar perpendicular components share $a_\perp,b_\perp$ in this study:

$$
\kappa_\parallel=K_0\beta\frac{B_n}{B_0}G_C(\mathcal R;a_\parallel,b_\parallel),\quad
\kappa_{\perp r}=0.02K_0\beta\frac{B_n}{B_0}G_C(\mathcal R;a_\perp,b_\perp),
\tag{G5}
$$

$$
\kappa_{\perp\theta}=0.01u(\theta)K_0\beta\frac{B_n}{B_0}
G_C(\mathcal R;a_\perp,b_\perp),\qquad
u(\theta)=\frac32+\frac12\tanh\left[h_C\left(|\theta-\pi/2|-\frac{35\pi}{180}\right)\right].
\tag{G6}
$$

In (G6), $\theta$ and the converted $35^\circ=35\pi/180$ are in radians; $h_C>0$ is a REQUIRED input in $\mathrm{rad^{-1}}$. Source Equation (11) prints the multiplier 8 alongside $\pi/2$ and $35^\circ$, without establishing the multiplier's angular-unit interpretation. Revision 2.0's unqualified use of 8 with converted radians is therefore withdrawn as a certified historical reproduction. This parameterized formula is usable only with an explicitly supplied width and angular convention; the exact source preset remains `source_gate` until that convention is established.

Algebraically (G6)'s angular function equals (G4) with $d=2$, $h_\theta=h_C$ and $\theta_F=35\pi/180$, using $|\theta-\pi/2|=\pi/2-\theta_A$. This identity does not resolve the source's unit ambiguity. The source fixes $B_n=1\,\mathrm{nT}$, $\mathcal R_k=4.3\,\mathrm{GV}$ and $s_R=2.2$ in this particular study, and fits the slopes and normalization. Those study-specific numbers were checked in the primary PDF. They are not universal defaults. The angular limits of 1 and 2 are approached, not exactly attained at finite equatorial or polar angles. The absolute value gives the same equatorial derivative cusp issue as (G4).

Equations (G5)/(G6) abbreviate the common $\mathcal R_k,s_R$ arguments of (G2). Separate radial and polar perpendicular slopes would be a different parameterized variant; they MUST not be attributed to this source preset.

Require the resolved source angular convention, fitted slopes and $K_0$ from the intended time interval before exposing an exact AMS-02 preset. A coefficient with substituted generic slopes is only a parameterized family member. A numeric $K_0$ quoted in units of $6\times10^{20}\,\mathrm{cm^2\,s^{-1}}$ is multiplied by $6\times10^{16}\,\mathrm{m^2\,s^{-1}}$ exactly once.

### 14.4 HelMod-style ratio

A ratio-family implementation can use

$$
\kappa_{\perp r}=\rho_H\kappa_\parallel,\qquad
\kappa_{\perp\theta}=\rho_H\kappa_\parallel F_H(\theta).
\tag{G7}
$$

It requires $\rho_H$ and a complete supplied polar function or table. A source-specific HelMod preset is gated until its latitude coordinate, enhancement rule, parameter set and temporal model are source checked. A step function is a discontinuous coefficient requiring appropriate interface treatment in the transport solver; central differencing across it is not an exact physical gradient. No inherited secondhand enhancement factor is made a default.

### 14.5 Drift convention and reduction

Use the signed weak-scattering coefficient

$$
\kappa_A^{\rm ws}=\operatorname{sgn}(q_e)\frac{vr_L}{3}
=\frac{pv}{3q_eB_0},\qquad
\mathbf v_d=\nabla\times(\kappa_A\mathbf b).
\tag{G8}
$$

With this convention the charge sign resides in $\kappa_A$ and field polarity in $\mathbf b$. Do not include either sign a second time. For $\kappa_{A,ij}=\kappa_A\epsilon_{ijk}b_k$, $\nabla\cdot(\mathbf K_A\nabla f)=-\mathbf v_d\cdot\nabla f$, consistent with (T1).

An explicitly prescribed rigidity reduction is

$$
\kappa_A=f_A\kappa_A^{\rm ws},\qquad
f_A=K_{A0}\frac{(\mathcal R/\mathcal R_A)^2}{1+(\mathcal R/\mathcal R_A)^2},
\quad 0\le K_{A0}\le1,\quad\mathcal R_A>0.
\tag{G9}
$$

Both parameters are configuration inputs. The classical factor is (I14), not automatically the same as a turbulence suppression model. For the Candia–Roulet empirical Hall fit [R19], a separate backend is

$$
\kappa_A=\operatorname{sgn}(q_e)\frac{cr_L}{3}
\left[1+(\sigma^2/\sigma_0^2)^2\right]^{-1/2},\qquad
\sigma_0^2=N_A\begin{cases}\rho_C^{0.3},&\rho_C\le0.2,\\
1.9\rho_C^{0.7},&\rho_C>0.2.
\end{cases}
\tag{G10}
$$

Registry `drift_candia_roulet_2004` uses (G10), the $N_A$ amplitudes in Section 13.2, and the same relativistic and radius definitions. It MUST carry a separately declared Hall-fit domain and out-of-domain policy; it does not inherit the perpendicular-pair domain. The source's Hall comparisons sample $\rho_C=0.01,0.1,1$ (Figure 7), which does not establish a uniformly calibrated filled parameter box. An exact calibrated preset remains gated if its domain metadata are missing. A supplied software or exploratory domain MUST be labeled as such. The two approximate fitted branches need not join exactly at 0.2; do not smooth them without a separate identified prescription. Do not apply this Galactic synthetic-field fit to a heliospheric coefficient merely because both have the same units.

For a spatially varying suppression,

$$
\nabla\times(f_A\kappa_A^{\rm ws}\mathbf b)
=f_A\nabla\times(\kappa_A^{\rm ws}\mathbf b)
+\nabla f_A\times(\kappa_A^{\rm ws}\mathbf b).
\tag{G11}
$$

Multiplying an existing weak-scattering velocity by $f_A$ omits the second term. It is equivalent to the reduced Parker coefficient only for a uniform factor or when the neglected term vanishes. A caller can select that approximation explicitly, with provenance; it MUST not be described as the exact curl of the reduced coefficient.

Parker's pitch-averaged antisymmetric transport and explicit pitch-angle-dependent guiding-centre gradient/curvature drifts have different contracts. Do not apply one isotropic reduction indiscriminately to a focused drift or double count an antisymmetric coefficient and explicit drifts. Current-sheet treatment belongs to the solver. The coefficient library supplies optional $\kappa_A$, $f_A$ and their available derivatives, not an inferred current-sheet crossing rule.

Unresolved Burger–Visser/Minnie formula variants and band-limited turbulence drift reductions remain source gates; their undefined logarithm bases, scales or amplitudes are not filled by analogy.

## 15. Stochastic transport, gradients and geometry

### 15.1 Forward density and backward generator

For the Cartesian Itô process

$$
dX_i=a_i\,dt+B_{i\alpha}\,dW_\alpha,\qquad
\mathbf B\mathbf B^T=2\mathbf K_s,
\tag{SDE1}
$$

the forward **number-density** equation and backward generator are

$$
\partial_t n=-\partial_i(a_i n)+\partial_i\partial_j(K_{ij}n),\qquad
L\varphi=a_i\partial_i\varphi+K_{ij}\partial_i\partial_j\varphi.
\tag{SDE2}
$$

For pure conservative spatial diffusion $\partial_t n=\nabla\cdot(\mathbf K_s\nabla n)$, choose

$$
a_i=\partial_jK_{ij}=(\nabla\cdot\mathbf K_s)_i.
\tag{SDE3}
$$

For the same diffusion operator acting on a scalar test function $\varphi$, its backward generator is $(\nabla\cdot\mathbf K_s)\cdot\nabla\varphi+\mathbf K_s:\nabla\nabla\varphi$. By contrast, $a=0$ has backward operator $\mathbf K_s:\nabla\nabla\varphi$ and forward operator $\partial_i\partial_j(K_{ij}n)$. A constant test function is stationary for the former; it need not give uniform forward particle density. This removes the forward/backward confusion in the earlier guide.

In one dimension with reflecting zero-flux boundaries, positive $\kappa(x)$ and zero additional flow:

$$
\begin{array}{ll}
a=\kappa': & J=-\kappa n',\quad n_{\rm stat}=\text{constant},\\
a=0: & J=-(\kappa n)',\quad n_{\rm stat}\propto1/\kappa.
\end{array}
\tag{SDE4}
$$

These analytic stationary densities are verification targets; an unattached Monte Carlo density ratio is not. The boundaries, sampling measure, initialization, relaxation and timestep MUST be specified when testing them.

For the spatial part of (T1), a backward Kolmogorov representation uses $a=\nabla\cdot\mathbf K_s-\mathbf U-\mathbf v_d$ along with the appropriate momentum term and time dependence. A forward physical-particle number-density mover requires a separate derivation for its phase-space measure. Do not copy this backward drift into a forward mover on the strength of a shared coefficient name.

### 15.2 Full tensor divergence

For the axisymmetric tensor (T3), the product rule gives

$$
\nabla\cdot\mathbf K_s
=\nabla\kappa_\perp
+\mathbf b[\mathbf b\cdot\nabla(\kappa_\parallel-\kappa_\perp)]
+(\kappa_\parallel-\kappa_\perp)
\left[(\mathbf b\cdot\nabla)\mathbf b+\mathbf b\nabla\cdot\mathbf b\right].
\tag{SDE5}
$$

The full gradient $\nabla\kappa_\perp$ is required. Knowing only its derivative along $\mathbf b$ is insufficient. If a code replaces the first term by $\mathbf b(\mathbf b\cdot\nabla\kappa_\perp)$, its omission is exactly $(\mathbf I-\mathbf b\mathbf b)\nabla\kappa_\perp$. This is a conditional algebraic diagnosis, not a verified claim about the current AMPS source.

For two perpendicular eigenvalues, write $\Delta_\parallel=\kappa_\parallel-\kappa_{\perp2}$ and $\Delta_1=\kappa_{\perp1}-\kappa_{\perp2}$. Then

$$
\begin{aligned}
\nabla\cdot\mathbf K_s={}&\nabla\kappa_{\perp2}
+\mathbf b(\mathbf b\cdot\nabla\Delta_\parallel)
+\Delta_\parallel[(\mathbf b\cdot\nabla)\mathbf b+\mathbf b\nabla\cdot\mathbf b]\\
&+\mathbf e_1(\mathbf e_1\cdot\nabla\Delta_1)
+\Delta_1[(\mathbf e_1\cdot\nabla)\mathbf e_1+\mathbf e_1\nabla\cdot\mathbf e_1].
\end{aligned}
\tag{SDE6}
$$

Frame derivatives are part of the operator. $\nabla\cdot\mathbf B_0=0$ does not imply $\nabla\cdot\mathbf b=0$. For constant scalar $\kappa_\perp$, the perpendicular-only tensor still has

$$
\nabla\cdot[\kappa_\perp(\mathbf I-\mathbf b\mathbf b)]
=-\kappa_\perp[(\mathbf b\cdot\nabla)\mathbf b+\mathbf b\nabla\cdot\mathbf b].
\tag{SDE7}
$$

A focused transport mover MUST declare whether its transverse operator is conservative divergence form. It cannot decide whether to include (SDE7) solely from the fact that the scalar coefficient is constant. A focused pitch-angle coefficient replaces $\kappa_\perp$ by $D_\perp(\mu)$ for that operator, with spatial derivatives at fixed $\mu,p,t$.

### 15.3 Noise, frame and coordinate conversion

A Cartesian covariance realization of (T2) is

$$
\mathbf B=
\left[\sqrt{2\kappa_\parallel}\mathbf b,
\sqrt{2\kappa_{\perp1}}\mathbf e_1,
\sqrt{2\kappa_{\perp2}}\mathbf e_2\right],\qquad
\Delta\mathbf X_{\rm noise}=\mathbf B\sqrt{\Delta t}\,\boldsymbol\xi,
\quad\boldsymbol\xi\sim N(\mathbf0,\mathbf I).
\tag{SDE8}
$$

The three scalar normal variates are independent. An exactly zero eigenvalue contributes exactly zero noise. A focused solver with deterministic parallel streaming uses only the appropriate transverse columns and its own pitch-angle scattering; it MUST not add parallel spatial noise on top of the same resolved parallel dynamics.

For unequal eigenvalues, a physical reference vector $\mathbf r_1$ defines

$$
\mathbf e_1=\frac{\mathbf r_1-(\mathbf r_1\cdot\mathbf b)\mathbf b}
{\sqrt{|\mathbf r_1|^2-(\mathbf r_1\cdot\mathbf b)^2}},\qquad
\mathbf e_2=\mathbf b\times\mathbf e_1.
\tag{SDE9}
$$

Require a nonsingular denominator. At degeneracy, an alternative physically defined frame or a failure is needed. Replacing unequal eigenvalues by their average is a model alteration and is not an automatic fallback. A reference-frame change that permutes eigenvalues MUST also permute their gradients and random-stream association.

Physical orthonormal spherical components in (T4) are **not** covariance entries for coordinates $(r,\theta,\phi)$. If a coordinate SDE is used,

$$
K^{rr}=K_{rr},\quad K^{\theta\theta}=K_{\theta\theta}/r^2,\quad
K^{\phi\phi}=K_{\phi\phi}/(r^2\sin^2\theta),\quad
K^{r\phi}=K_{r\phi}/(r\sin\theta),
\tag{SDE10}
$$

with analogous factors for other cross terms and the Itô coordinate drift. Derive that drift by Itô's formula, including the spherical Jacobian; a Cholesky factor of (T4) alone does not supply it. The reference implementation SHOULD advance Cartesian positions to avoid unnecessary coordinate singularities and then transform outputs for analysis.

### 15.4 Derivatives, interfaces and timesteps

Coefficient derivatives are spatial derivatives at fixed particle state and time. At each neighboring position, query the same provider and reevaluate the complete selected coefficient, including its parallel dependency or coupled pair. Differencing only its explicit radial factor misses the implicit turbulence and parallel dependencies.

An independent Cartesian reference for tensor divergence is

$$
(\nabla\cdot\mathbf K_s)_i\simeq
\sum_{j=1}^3\frac{K_{ij}(\mathbf x+h_j\mathbf e_j)-K_{ij}(\mathbf x-h_j\mathbf e_j)}{2h_j}.
\tag{SDE11}
$$

For smooth fields it has second-order truncation error. Specify $h_j$ and refine it relative to interpolation and solver noise. Six neighbor evaluations plus the center are seven samples; that count is not a measured seven-to-three AMPS runtime ratio. At boundaries use an explicitly derived one-sided stencil or report unavailable derivatives. At discontinuous coefficient interfaces use a solver scheme consistent with the weak transport operator; an unresolved central difference is not a physical interface law.

Example engineering timestep controls are

$$
\Delta t\le f_D\frac{h^2}{2\kappa_{\max}},\qquad
\Delta t\le f_a\frac h{|\mathbf a|},\qquad
\kappa_{\max}=\max(\kappa_\parallel,\kappa_{\perp1},\kappa_{\perp2}),
\tag{SDE12}
$$

with user/solver-specified fractions and zero-rate handling. They MUST be supplemented by snapshot boundaries, coefficient variation scales, field-frame curvature, focused scattering and the solver's existing controls. Eigenvalue gradients alone miss frame-curvature terms in (SDE5)–(SDE7). A zero eigenvalue MUST not be divided into a nonzero gradient; use a well-defined tensor-level control and timestep-convergence evidence. No universal numerical fraction is prescribed here.

With $\kappa_\perp=0$, rank-one parallel noise is tangent to a curved field locally. A finite Euler step need not stay exactly on the curved field line. Verify convergence to the intended geometric process; exact field-line following at finite steps requires a line-coordinate method.

## 16. Library interface, solvers and tables

### 16.1 Host-neutral inputs and typed outputs

The coefficient library has no mesh, PIC, MPI, trajectory, random-number or source-injection ownership. Input consists of an immutable configuration, a local particle state, a background/turbulence sample and any explicitly supplied parallel dependency. All formulas use SI after a single adapter conversion.

The following is a proposed host-neutral C++20 declaration contract, not a compiled AMPS implementation. The standard-library includes and provider details are omitted; named application records are defined below. Provider spectra are the tagged contract of Section 4, not an unspecified inferred spectrum.

```cpp
using Vec3 = std::array<double,3>;
using Mat3 = std::array<std::array<double,3>,3>;

enum class Status {
  success, invalid_input, missing_input, missing_calibration,
  outside_model_domain, incompatible_geometry, inconsistent_pair,
  divergent_moment, integration_failed, nonlinear_solver_failed,
  derivative_unavailable, source_gate
};
enum class Observable {
  symmetric_coefficient, paired_coefficients, pitch_angle_coefficient,
  signed_hall_coefficient, particle_msd, centered_particle_covariance,
  field_line_coefficient, field_line_msd, relative_to_line_msd,
  conditional_return_median, perturbative_diagnostic
};
enum class Estimator { asymptotic, instantaneous_derivative, secant,
                       finite_window_average, moment_value, not_applicable };
enum class FrameKind { axisymmetric_ordered, oriented_unequal, isotropic };
enum class DomainState { inside_declared_domain, extrapolated, not_applicable };
enum class Quality { source_fit, parameterized_closure, scaling_estimate,
                     diagnostic };
enum class DiffusionRegime { normal_diffusion, no_normal_diffusion,
                             not_established, not_applicable };
enum class DependencyOwner { supplied_parallel, paired_backend, none };
enum class LineCoordinate { mean_field_z, physical_arc_length };

struct ParticleState {
  double mass_kg, charge_C, momentum_kg_m_s;
  std::optional<double> mu;
};
struct ParallelInput {
  double kappa_m2_s;
  std::string model_id, equation_version, sample_fingerprint;
};
struct TurbulenceSample {
  std::string geometry_id, spectrum_id, energy_convention;
  std::string provider_revision, sample_fingerprint;
  std::map<std::string,double> tagged_si_parameters;
  // Required keys/units are validated against the selected Section 4 schema.
  // General spectral tensors are immutable provider handles in the actual API.
};
struct LocalState {
  Vec3 position_m;
  double time_s;
  std::optional<Vec3> mean_field_T;
  TurbulenceSample turbulence;
  std::optional<ParallelInput> parallel_dependency;
};
struct PhysicalFrame {
  FrameKind kind;
  std::optional<Vec3> b, e1, e2;
  std::string perpendicular_axis_definition;
};
struct SpatialCoefficients {
  PhysicalFrame frame;
  std::optional<double> parallel_m2_s;
  std::optional<double> perp1_m2_s, perp2_m2_s;
  std::optional<double> isotropic_m2_s;
  DependencyOwner parallel_owner;
};
struct ParticleMoments {
  double age_s;
  Vec3 raw_msd_m2;
  std::optional<Vec3> mean_displacement_m;
  std::optional<Mat3> centered_covariance_m2;
  std::optional<Vec3> derivative_m2_s, secant_m2_s;
};
struct FieldLineMoments {
  double coordinate_m;
  LineCoordinate coordinate_kind;
  std::optional<double> msd_per_transverse_axis_m2;
  std::optional<double> derivative_half_m, secant_half_m;
};
struct NamedSensitivity {
  std::string quantity, independent_input, unit;
  double value;
  std::optional<double> estimated_absolute_error;
};
struct NamedGradient {
  std::string quantity, unit;
  Vec3 cartesian_gradient;
  std::optional<double> estimated_absolute_error;
};
struct Provenance {
  std::string requested_model, used_model, equation_version, source_version;
  std::optional<std::string> fallback_reason;
  std::string configuration_fingerprint, sample_fingerprint;
  std::optional<std::string> calibration_id, table_checksum;
  std::map<std::string,std::string> conventions;
  DependencyOwner dependency_owner;
  std::optional<std::string> parallel_model_id;
};
struct NumericalDiagnostics {
  std::string backend, integration_method, root_method;
  std::optional<double> absolute_error_estimate, relative_error_estimate;
  std::optional<double> residual_norm, jacobian_condition_estimate;
  std::map<std::string,double> tolerances_and_tail_controls;
  std::vector<std::string> approximations, unavailable_derivatives;
  std::optional<std::size_t> iterations;
};
struct ModelResult {
  Status status;
  Observable observable;
  Estimator estimator;
  DomainState domain;
  Quality quality;
  DiffusionRegime diffusion_regime;
  std::optional<SpatialCoefficients> coefficients;
  std::optional<double> pitch_angle_m2_s;
  std::optional<double> field_line_m;
  std::optional<LineCoordinate> field_line_coordinate;
  std::optional<ParticleMoments> particle_moments;
  std::optional<FieldLineMoments> field_line_moments;
  std::optional<Vec3> relative_to_line_msd_m2;
  std::optional<double> conditional_median_m2, signed_hall_m2_s;
  std::vector<NamedSensitivity> sensitivities;
  std::vector<NamedGradient> gradients;
  Provenance provenance;
  NumericalDiagnostics numerical;
};
```

The implementation MUST enforce the following invariants. Optional quantities mean unavailable unless supplied by the selected model; absence is distinct from physical zero. This sketch uses one result record for exposition; a tagged union with the same invariants is permitted.

| Contract case | REQUIRED behavior |
|---|---|
| Axisymmetric ordered tensor | Finite unit $\mathbf b$ from a nonzero mean field; equal perpendicular eigenvalues. A perpendicular-only evaluation does not fabricate a parallel value. |
| Unequal perpendicular eigenvalues | Return physically defined $\mathbf b,\mathbf e_1,\mathbf e_2$ and axis definition, or an equivalent complete Cartesian tensor. Validate normalization, orthogonality and handedness. |
| Reference-vector construction | Require a supplied physical reference and a nonzero perpendicular projection; reject a parallel or unavailable reference. An arbitrary axis changes the model when eigenvalues differ. |
| Zero-mean-field isotropic tensor | Return `FrameKind::isotropic` and $\kappa_{\rm iso}$, with no axes. Never divide by a zero field norm. A supplied zero vector can represent zero field; a missing vector means unavailable. |
| Paired backend | Both parallel and perpendicular outputs have the same owner, configuration and sample. Never mix NLGCE-F perpendicular with an unrelated parallel result. |
| Diagnostic statistics | Preserve particle/field-line/relative/conditional observables and estimator tags. Field-line derivative is per line coordinate, not per time. A conditional median is not a covariance or noise amplitude. |
| Valid anomalous transport | Use `Status::success` with `DiffusionRegime::no_normal_diffusion` when that limit is established. Retain valid finite-time moments; supply a zero asymptotic coefficient only where its limit is proved, as in (C5)/(F8). |
| Failure | Keep physical output optionals unavailable. Preserve diagnostics and requested identity; do not represent a failed root, missing calibration or source gate as zero. |
| Explicit extrapolation/fallback | Set the domain and quality tags; preserve requested and used identities and the reason. A fallback success cannot hide the original failure. |
| Source preset | Require all source-specific definitions, fitted data and domain metadata. A parameterized family member is labeled separately. |

Domain records MUST identify their boundaries and whether they are source-stated, conservative software limits, or supplied exploratory limits. A rejected outside-domain request has no physical outputs; an explicitly enabled extrapolation is a successful tagged evaluation. Numerical error estimates remain separate from fit uncertainty and the quality label.

Define named independent sensitivity coordinates per model, for example $(p,B_0,\delta B_s^2,\delta B_2^2,\ell_s,\ell_2)$ at fixed spectral indices. Do not mix dependent particle variables or derived lengths in an undocumented array. Each sensitivity declares its output/input units. Particle-coefficient spatial gradients have units $\mathrm{m\,s^{-1}}$; field-line-coefficient gradients are dimensionless. Missing sensitivities or gradients are explicitly unavailable and do not invalidate an otherwise valid scalar evaluation unless the caller requested derivatives as a required output. Full tensor divergence also requires the frame derivatives of Section 15.

### 16.2 Numerical integrals and implicit roots

For a one-dimensional integral, use dimensionless wavenumber followed by $x=e^y$ if helpful:

$$
\int_0^\infty f(x)dx=\int_{-\infty}^{\infty}e^y f(e^y)dy.
\tag{NUM1}
$$

Control both tails from the actual power laws or convergence checks. In expressions such as $(1+e^{2y})^{-a}$ use stable logarithmic forms. Finite $y$ limits are numerical quadrature controls, not physical cutoffs. Separate a missing divergent spectral moment from a finite closure integral with additional suppressing factors.

For $\kappa_\perp=R(\kappa_\perp;\mathbf u)$, solve in $y=\ln[\kappa_\perp/(v\ell)]$ or in log coefficient ratio. A convenient dimensionless residual is

$$
r(y)=\ln R(v\ell e^y;\mathbf u)-\ln(v\ell e^y).
\tag{NUM2}
$$

Bracket a positive root using the model's known limits; do not assume every general kernel has a positive root. For finite $Q>0$ and convergent spectra, the specified kernels obey:

| Model / geometry | Upper bound |
|---|---|
| Composite NLGC (N1)/(N2) | $\kappa_\perp\le a^2\kappa_\parallel(\delta B_s^2+\delta B_2^2)/(2B_0^2)$ |
| 2D ENLGC (N5) | $\kappa_\perp\le\kappa_\parallel\delta B_2^2/(2B_0^2)$ |
| UNLT with transverse slab plus 2D (U2) | $\kappa_\perp\le a^2\kappa_\parallel\delta B_2^2/(2B_0^2)$ |
| General UNLT (U1) | $\kappa_\perp\le a^2\kappa_\parallel\int P_{xx}d^3k/B_0^2$ |
| Exact/rational implicit slab (U7)/(U9) | $\kappa_\perp\le\kappa_\parallel\delta B_2^2/(2B_0^2)$ |
| FLPD (D1)/(D2) | $\kappa_\perp\le\kappa_\parallel\delta B_2^2/(2B_0^2)$ |
| RBD/BC (B1) | $\kappa_\perp\le a^2\kappa_\parallel\int P_{xx}d^3k/B_0^2$ |
| Coupled NLGCE perpendicular equation (E4) | $\kappa_\perp\le a'^2\kappa_\parallel\delta B_2^2/(2B_0^2)$, with both coefficients outputs |
| Closed form (B6) | $\lambda_\perp\le\min(\kappa_{\rm FL}^2\lambda_\parallel/\ell^2,3\kappa_{\rm FL}/2)$ |

The 2D-only UNLT bound follows from that geometry; it is not a rule for arbitrary 3D spectra. For RBD, $T_{\rm BC}\le T_{\rm plain}\le1/Q$ follows from (B9). These bounds help bracketing and invariant checks; violating one is an evaluation error, not a reason to clip the result. Evaluate the continuous zero-parallel-input branches in Section 1.2 before constructing $Q$ or logarithms. Do not form $0/0$ ratios. The paired NLGCE bound constrains its output; it does not create an external parallel input.

For scalar implicit sensitivity, if $F(\kappa,\mathbf u)=\kappa-R(\kappa;\mathbf u)$,

$$
\frac{\partial\kappa}{\partial u_a}
=-\frac{\partial F/\partial u_a}{\partial F/\partial\kappa}.
\tag{NUM3}
$$

For a coupled residual $\mathbf F(\mathbf y,\mathbf u)=0$,

$$
\frac{\partial\mathbf y}{\partial u_a}
=-\left(\frac{\partial\mathbf F}{\partial\mathbf y}\right)^{-1}
\frac{\partial\mathbf F}{\partial u_a}.
\tag{NUM4}
$$

Report an ill-conditioned or unavailable derivative instead of fabricating a sensitivity from an asymptotic exponent. Differentiation under an integral requires its convergent model domain. Coefficient and gradient error budgets MUST both include integral, root and interpolation errors.

### 16.3 Tables and caching

A supplied table header MUST include: canonical model and equation version; quantity and per-axis convention; independent axes and units; spectrum type and normalization; radius and scale definitions; provider provenance; any parallel dependency or paired-model identity; calibration inputs; numerical construction method; domain; interpolation and extrapolation policies; value ordering; and content checksum.

Use logarithms only for strictly positive axes and values. Exact zero coefficients need a separate zero mask/branch or a specified linear scheme. Natural-log multilinear interpolation is a reproducible baseline; it is continuous but its derivative is generally discontinuous at cell boundaries. A named derivative-aware interpolant may be used only after positivity, overshoot, boundary and gradient-error tests. Do not assert that a generic monotone cubic or a folded polar formula is globally $C^1$.

For signed Hall tables, the header MUST additionally declare how signed values, zeros and sign changes are interpolated. Logging a negative or zero $\kappa_A$ is invalid. Interpolation of signed values on declared coordinates is a permissible engineering scheme; positive coordinates can use logarithmic spacing. Positive-magnitude interpolation is permitted only within a fixed-sign branch with an explicit zero/crossing policy. Keep charge-sign branches discrete. If no signed-value scheme is supplied, return `missing_input` rather than choosing one implicitly. Derivative and zero-crossing tests are required for the selected scheme.

Validate at held-out interior states, faces, corners, regime transitions and wherever the production solver differentiates the table. Compare values **and** derivatives with the direct reference. Record a maximum observed interpolation discrepancy, the validation sample definition and any untested domain; none is a proof of a global error bound.

A deterministic cache key MUST cover every value-relevant input: model/configuration version, particle state, spatial sample, time, provider revisions, local spectra/variances/scales and the parallel dependency or coupled-pair identity. Provider revision alone is insufficient when values vary spatially within that revision. Exact caching is distinct from momentum or space binning; binned reuse introduces an approximation whose errors MUST be measured and recorded. Never substitute a cached state across a snapshot boundary without a declared time rule.

### 16.4 Errors, provenance and configuration examples

Missing turbulence or calibration, an unsupported geometry, a failed root and out-of-domain input are distinguishable failures. The default behavior is to fail the requested model evaluation. An explicitly configured fallback returns both the requested and used model identities and its reason. It MUST not conceal a failed closure by silently returning a constant ratio.

Fingerprint canonical configuration fields individually, including units after conversion, source version, all active parameters, energy and length conventions, angular convention, table checksum, error controls and any approximation/fallback. Do not hash struct memory, unordered maps or only the coefficient's numeric value. Runtime sample fingerprints and configuration fingerprints serve different purposes. A checksum verifies unchanged data bytes; it does not prove a physical model.

Example **model-only** configurations below use dimensionless fixture states, not calibrated heliospheric profiles. The host MUST supply the background, species and the named parallel backend.

```toml
[perpendicular]
model = "flpd_complete"
backend = "direct"
parallel_input = "supplied"
out_of_domain = "fail"
field_line_closure = "dd_composite"

[turbulence]
geometry = "slab_plus_2d"
spectrum = "shalchi_composite"
inertial_index_s = 1.6666666666666667
energy_index_q = 3.0
variance_definition = "separate_component_variances"
length_definition = "bend_over"
```

```toml
[coefficients]
paired_model = "nlgce_f_2014"
out_of_domain = "fail"
parallel_table = "NLGCE_F_2014_parallel.csv"
perpendicular_table = "NLGCE_F_2014_perpendicular.csv"
```

These are proposed schema examples; they are not asserted to be accepted by today's AMPS parser. A real run MUST include every required dimensional sample and record the source-checked parameter set. The first example selects a perpendicular-only closure with an externally supplied parallel dependency. The second selects one backend that owns and supplies both coefficients. They are alternative configurations.

## 17. Numerical verification and acceptance

### 17.1 Reproducible scope

The reference script and machine-readable fixtures in the companion archive recompute the numerical tables below. The paired NLGCE script and its input tables are also included. These scripts are verification artifacts, not a production coefficient library. Source checking, mathematical consistency, numerical evaluation, transport-solver correctness and agreement with physical data are separate levels of evidence.

All numerical tolerances in the scripts are engineering choices. A comparison of a theory with a particle simulation, even agreement within a factor of two, MUST not be reused as a coefficient-solver tolerance. No particle simulation, AMPS test, fitted value for an unpublished parameter or measured runtime is represented by these fixtures.

Digits printed for equation fixtures describe numerical evaluation and rounding, not the physical precision of a closure or calibration. Run `python3 reference_verification.py` from the extracted data directory to recompute them and check the checksums. The revision 2.1 additions are independently recomputed by `python3 review_additions_verification.py`, using mpmath at the declared precision in addition to NumPy/SciPy. Run `python3 document_fixtures.py --document /path/to/PERPENDICULAR_DIFFUSION_COEFFICIENT_MODEL.md` to compare every generated table with its stored JSON, and `python3 document_audit.py` for structure and reference checks. The additional script also checks the printed $C(s)$, $\pi C(s)$, $L_s/\ell_s$ rounding and both quoted GCD coefficients. Default verification commands are read-only; fixture/report generation is an explicit maintenance operation.

### 17.2 Spectrum and field-line references

For $s=5/3$, the recomputed constants are $C(s)=0.118862354635443$, $\pi C(s)=0.373417100111093$ and $L_s/\ell_s=0.746834200222187$.

<!-- BEGIN GENERATED TABLE_LENGTHS -->
| $q$ | $L_U/\ell_2$ | $L_\perp/\ell_2$ | $L_\perp^2/(4L_U^2)$ |
| --- | --- | --- | --- |
| 1.5 | 1.15470054 | 1.13931016 | 0.243380181 |
| 2 | 0.816496581 | 0.950898837 | 0.339078224 |
| 3 | 0.577350269 | 0.7468342 | 0.418320992 |
<!-- END GENERATED TABLE_LENGTHS -->

For the next table, $B_0=1$, $\ell_s=\ell_2=1$ and $s=5/3$. All field-line coefficients are expressed in the common length unit. The state parameters are synthetic equation-verification inputs.

<!-- BEGIN GENERATED TABLE_FIELDLINES -->
| State | $q$ | $f_s$ | $\delta B^2/B_0^2$ | $\kappa_s$ | $\kappa_2$ | $\kappa_{\rm FL}$ |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 3 | 0 | 1 | 0 | 0.40824829 | 0.40824829 |
| 2 | 2 | 0.2 | 1 | 0.07468342 | 0.516397779 | 0.555087854 |
| 3 | 3 | 0.2 | 1 | 0.07468342 | 0.365148372 | 0.404394481 |
| 4 | 3 | 0.5 | 1 | 0.18670855 | 0.288675135 | 0.396748992 |
| 5 | 3 | 0.2 | 0.25 | 0.018670855 | 0.182574186 | 0.192148128 |
| 6 | 3 | 0.2 | 2 | 0.14936684 | 0.516397779 | 0.596453753 |
<!-- END GENERATED TABLE_FIELDLINES -->

### 17.3 FLPD references and limits

These values solve (D2), using the same states and (F4). Each cell is the per-axis ratio $\eta_\kappa=\lambda_\perp/\lambda_\parallel$. They are calculated values, not points digitized from a published simulation. The detailed JSON supplies more digits and the solver controls.

<!-- BEGIN GENERATED TABLE_FLPD -->
| State | $\lambda_\parallel/\ell_2=0.01$ | $\lambda_\parallel/\ell_2=0.1$ | $\lambda_\parallel/\ell_2=1$ | $\lambda_\parallel/\ell_2=10$ | $\lambda_\parallel/\ell_2=100$ | $\lambda_\parallel/\ell_2=1000$ |
| --- | --- | --- | --- | --- | --- | --- |
| 1 | 0.0905395152 | 0.0898707818 | 0.0756809986 | 0.0264989959 | 0.00391965196 | 0.000428639342 |
| 2 | 0.0572929984 | 0.0570843491 | 0.0517588756 | 0.0232140803 | 0.00411123254 | 0.000492334409 |
| 3 | 0.0647091203 | 0.0643472985 | 0.0559862612 | 0.0213919699 | 0.00329820227 | 0.00036461766 |
| 4 | 0.0311980426 | 0.0311025403 | 0.0285376566 | 0.0131438263 | 0.00224754814 | 0.000255546904 |
| 5 | 0.0171297787 | 0.0170951161 | 0.0160869009 | 0.00839921755 | 0.00157340928 | 0.000184024632 |
| 6 | 0.123304695 | 0.122191362 | 0.100082793 | 0.0327401962 | 0.00468474379 | 0.00050760348 |
<!-- END GENERATED TABLE_FLPD -->

The older table's $f_s=0.5$, $\lambda_\parallel/\ell_2=1000$ cell is corrected from 0.0002556 to **0.000255546904** at the precision shown here. Refining both log-wavenumber limits and quadrature tolerances reproduces the retained values. Independent dimensional (D1) and dimensionless (D2) evaluations are compared in the reference script.

The factor in (D8) is $0.713644179546180$. For the pure-2D, $q=3$, unit-variance state, (D7) gives $\lambda_\perp/\ell_2\to0.437016024448821$; the finite-path solution at $\lambda_\parallel/\ell_2=10^5$ is approximately $\eta_\kappa=4.36842679763362\times10^{-6}$. For the $q=3$, $f_s=0.2$, unit-variance composite state, the correct long-path limit is $\lambda_\perp/\ell_2\to0.372730281507756$.

The pure-2D short-path equation (D10) gives the following amplitude-independent results at $s=5/3$. Its ratios use positive $V_2$; neither a zero-power quotient nor a fitted turbulence amplitude is used.

<!-- BEGIN GENERATED TABLE_FLPD_SHORT -->
| $q$ | $\alpha_{\rm fluid}^{(2)}$ | $\alpha_{\rm D4}^{(2)}$ | $\zeta$ | Root $\alpha_{\rm FLPD}$ |
| --- | --- | --- | --- | --- |
| 1.5 | 1.15470053838 | 2.34059764404 | 0.12335251494 | 3.28772409187 |
| 2 | 0.816496580928 | 1.40218210533 | 0.155723479914 | 2.06908017903 |
| 3 | 0.57735026919 | 0.89265685271 | 0.181115665089 | 1.35662983904 |
<!-- END GENERATED TABLE_FLPD_SHORT -->

For $q=3$, $\eta_\kappa\to0.181115665089330\,V_2$ and $\alpha_{\rm FLPD}\to1.35662983904420$. These finite values support the implicit (D3), not either extreme asymptotic assumption. The additional fixtures cross-check (D10) with full (D2) at a chosen small path and verify the derivative bound used in Section 10.2.1.

The approach to (D3) for the uncut spectrum is nonanalytic. The following tests fix $\alpha$ at the $q=3$ short-path root and compute (D12)/(D13) by high-precision quadrature. The local exponent compares consecutive $A$ values; the first row has no preceding point. These are fixed-$\alpha$ integral checks, not full self-consistent finite-$A$ roots or claimed universal achieved error bounds.

<!-- BEGIN GENERATED TABLE_FLPD_TAIL -->
| $A$ | $J(0,\alpha)-J(A,\alpha)$ | Difference / $A^{5/6}$ | Local log exponent |
| --- | --- | --- | --- |
| 1e-08 | 2.52891516285e-07 | 1.17381843811 | unavailable |
| 1e-10 | 5.61371505253e-09 | 1.20943824491 | 0.826841947971 |
| 1e-12 | 1.22597150444e-10 | 1.22597150444 | 0.830384994275 |
<!-- END GENERATED TABLE_FLPD_TAIL -->

For pure 2D, $s=5/3$, $q=3$ and unit total variance, the CLRR simplification (D5) gives $\eta_\kappa\simeq0.209160495983068$. It does not describe the general short-path FLPD value $\eta_\kappa=0.090539515152693$ at $\lambda_\parallel/\ell_2=0.01$. This comparison demonstrates the additional $\alpha_{\rm FLPD}\gg1$ assumption. At these indices (D9) instead gives $\alpha_{\rm D4}^{(2)}=0.892656852710187$; the branch does not satisfy its own large-parameter assumption.

### 17.4 ENLGC, RBD, closed-form and kernel references

For (N5), $s=5/3$, $q=0$, $\delta B_2^2/B_0^2=0.8$ and $\ell_2=1$:

<!-- BEGIN GENERATED TABLE_ENLGC -->
| $\lambda_\parallel/\ell_2$ | $\lambda_\perp/\ell_2$ |
| --- | --- |
| 0.01 | 0.00391879938 |
| 0.1 | 0.0363550466 |
| 1 | 0.263129509 |
| 10 | 1.07757473 |
| 1000 | 6.36321559 |
<!-- END GENERATED TABLE_ENLGC -->

For (B3), use pure 2D, (S12), $a^2=1/3$, $s_{\rm A}=5/3$, $p=2$, $B_0=v=\lambda_2=1$. The numerical values and the exact long-path coefficient are regenerated, rather than mixing a rounded coefficient with a separately rounded limit:

<!-- BEGIN GENERATED TABLE_RBD -->
| $\delta B_2^2/B_0^2$ | $\lambda_\parallel/\lambda_2$ | $\kappa_\perp/(v\lambda_2)$ |
| --- | --- | --- |
| 0.1 | 100 | 0.0448755479 |
| 0.1 | 1000 | 0.049322237 |
| 0.4 | 100 | 0.0946606068 |
| 0.4 | 1000 | 0.0991443145 |
| 0.8 | 100 | 0.135927206 |
| 0.1 | $\infty$ | 0.0498221441 |
| 0.4 | $\infty$ | 0.0996442883 |
| 0.8 | $\infty$ | 0.140918304 |
<!-- END GENERATED TABLE_RBD -->

For (B6), the chosen input is $\kappa_{\rm FL}/\ell=0.1$:

<!-- BEGIN GENERATED TABLE_CLOSED -->
| $\lambda_\parallel/\ell$ | $\lambda_\perp/\ell$ |
| --- | --- |
| 0.1 | 0.000986884822 |
| 1 | 0.00885427379 |
| 10 | 0.0470789008 |
| 100 | 0.102075998 |
| 10000 | 0.144301936 |
<!-- END GENERATED TABLE_CLOSED -->

For the exact kernel (U6), the rational error means $\mathcal K_{\rm rat}/\mathcal K-1$:

<!-- BEGIN GENERATED TABLE_IMPLICIT_KERNEL -->
| $\xi$ | $\mathcal K(\xi)$ | Rational relative error |
| --- | --- | --- |
| 0.01 | 0.982473702 | 0.0176354214 |
| 0.1 | 0.841107137 | 0.165597239 |
| 0.5 | 0.454358639 | 0.467269705 |
| 1 | 0.242127844 | 0.376683194 |
| 2 | 0.0946459 | 0.173966448 |
| 5 | 0.0189056927 | 0.0371396311 |
| 10 | 0.00492681218 | 0.00980597613 |
<!-- END GENERATED TABLE_IMPLICIT_KERNEL -->

At large arguments the rational kernel is an asymptotic approximation; at $\xi=2$ its relative discrepancy is much larger than the quadrature tolerance. The exact and rational implicit-slab coefficients MUST therefore remain distinct backend identities.

The source integral-length profile in (B8), evaluated for the same smooth shape at $s=5/3$, is:

<!-- BEGIN GENERATED TABLE_CLOSED_LENGTHS -->
| $q$ | $\ell_{\perp,\rm int}/\ell_2$ |
| --- | --- |
| 1.5 | 0.569655077937 |
| 2 | 0.475449418542 |
| 3 | 0.373417100111 |
<!-- END GENERATED TABLE_CLOSED_LENGTHS -->

The direct moment integral and the gamma-function result agree under the numerical controls recorded in the added fixtures. This table converts definitions; it does not recalibrate a closed-form model.

For the same finite RBD variance/path family, independent evaluation of (B9) gives the following values (including the additional $0.8,1000$ state):

<!-- BEGIN GENERATED TABLE_RBD_COMPARE -->
| $\delta B_2^2/B_0^2$ | $\lambda_\parallel/\lambda_2$ | RBD/BC | Plain RBD |
| --- | --- | --- | --- |
| 0.1 | 100 | 0.0448755479307 | 0.0453317062419 |
| 0.1 | 1000 | 0.0493222370064 | 0.0493280794052 |
| 0.4 | 100 | 0.0946606068342 | 0.0949188156801 |
| 0.4 | 1000 | 0.0991443144773 | 0.0991472948213 |
| 0.8 | 100 | 0.135927205591 | 0.136117432316 |
| 0.8 | 1000 | 0.14041831775 | 0.140420438879 |
<!-- END GENERATED TABLE_RBD_COMPARE -->

The corrected and plain kernels share their long-path limit but differ at finite $Q$. The added fixtures independently integrate the Gaussian time transform before applying its backtracking factor. The plain numbers are an equation diagnostic, not an additional production model.

### 17.5 Additional equation checks

The companion calculations also check: spectrum normalization and permitted moments; (N3) against direct quadrature; (N2) against (N4); (U2) against its reduced integral; (U3) normalization and (U4) against nested quadrature; (U6) against (U12); GCD's coefficient against a direct Gaussian moment integral; (I10) against integration of (I9); the chain rule in (I12); (C10) and (C11) asymptotics with declared artificial fit inputs; classical Hall identities; tensor symmetry/eigenvalues/covariance; and analytic tensor divergence against Cartesian differences in a smooth curved field.

Parker geometry examples use **chosen inputs**, not fitted solar constants: $\Omega_{\rm ex}=2.87\times10^{-6}\,\mathrm{s^{-1}}$, $V_{\rm ex}=400\,000\,\mathrm{m\,s^{-1}}$, equatorial angle $\tan\psi=\Omega_{\rm ex}r/V_{\rm ex}$, and the exact au conversion. The reference values at 1 au are $\psi=47.0265300605^\circ$, $\cos^2\psi=0.464659869178404$ and $\cos\psi=0.681659643207960$. At $r=r_0=1\,\mathrm{au}$ with an explicitly chosen $\alpha_D=0.3$, (P6) gives 0.160612269551308, demonstrating the distinction between a prefactor and an averaged ratio. None of these example choices is a recommended default.

The archived paired NLGCE states A–G are now explicit here. They use $v=\ell_s=B_0=1$; $r=r_L/\ell_s$, $\epsilon^2=\delta B^2/B_0^2$, and the length ratio is $\ell_s/\ell_2$. These are synthetic equation fixtures, not observational samples.

<!-- BEGIN GENERATED TABLE_NLGCE_STATES -->
| State | $r$ | $f_s$ | $\epsilon^2$ | $\ell_s/\ell_2$ |
| --- | --- | --- | --- | --- |
| A | 0.035848041610664981 | 0.20000000000000001 | 1 | 10 |
| B | 0.01 | 0.20000000000000001 | 0.25 | 10 |
| C | 0.10000000000000001 | 0.5 | 1 | 1 |
| D | 0.01 | 0.20000000000000001 | 0.01 | 10 |
| E | 0.001 | 0.20000000000000001 | 0.001 | 10 |
| F | 0.01 | 0.20000000000000001 | 0.0001 | 10 |
| G | 0.01 | 0.20000000000000001 | 0.050000000000000003 | 10 |
<!-- END GENERATED TABLE_NLGCE_STATES -->

The actual procedure in `paired_reference_verification.py` initializes both nonlinear paths from the polynomial pair, repeats with both starting paths multiplied by 3, and compares quadrature controls $|y|\le40$, tolerance $10^{-10}$ against $|y|\le45$, tolerance $3\times10^{-11}$. Residuals are formed in log coefficients; the stored report records the results for these seven states. This tests two starts per state and does not establish root uniqueness or full fit-box coverage. Production tests MUST additionally enumerate faces, corners and transition states of the selected domain rather than treating A–G as exhaustive.

For $s=5/3$, the secant-average shape coefficient in (C6) is $(3/2)^{4/3}\alpha_G=1.73466066946129$. Equation (C12) gives the derivative-average coefficient by multiplying this value by $2/3$; it is a different estimator, not a different fitted normalization.

### 17.6 Production acceptance suite

The roadmap's implementation MUST add tests beyond the delivered algebraic reference checks:

1. **Contracts:** unit conversions; missing source angular, length and domain metadata; missing calibration; geometry/radius/scale mismatches; zero-mean isotropic output without axes; unequal-axis physical orientation; failed roots and unavailable outputs; success with no normal diffusion; and distinct observable/estimator routing. Other Snodin spectral shapes MUST fail an unsupported ordered-field pair request.
2. **Numerics:** independently reduced integrals; refined tails/tolerances; every Section 16.2 bound and zero-parallel branch without $Q$ or $0/0$; (D9)/(D10) consistency; the nonanalytic uncut FLPD limit; corrected/plain RBD separation; tagged closed-form lengths; positive-root residuals; multiple starts for NLGCE; enumerated boundary/corner states; and named sensitivities against independent finite differences.
3. **Assets:** all 576 NLGCE coefficients by explicit index; raw-file digests; polynomial results using direct summation and Horner; paired ownership; $\epsilon\ne1$ cases that detect the erroneous exponent in $a_x$.
4. **Geometry:** covariance $\mathbf B\mathbf B^T=2\mathbf K_s$, rotations, semidefinite zeros, unequal-axis frame orientation, tensor-divergence comparisons and spherical conversion factors.
5. **Transport:** homogeneous Gaussian covariance growth, zero/perpendicular limits, the analytic stationary densities (SDE4), boundary flux, forward/backward operator consistency and curved-field timestep convergence.
6. **Memory:** frozen-line out/back retracing, line-sharing, timestep and line-resolution convergence, compound $t^{1/2}$ scaling and restart/order reproducibility.
7. **Tables and cache:** held-out values and gradients, edges/corners, signed drift interpolation and zero crossings, discrete charge branches, explicit interpolation policy, zero handling, provider/time/momentum invalidation and quantified binned-cache error.

Use an analytic manufactured solution or deterministic operator check where available. Stochastic acceptance needs specified seeds, sample counts, uncertainty estimates, burn-in/fit windows and timestep refinement. Its statistical tolerance MUST be set from that design; no sample counts, achieved errors or density ratios are invented here.

## 18. Publication reporting and revision record

### 18.1 Required reporting

A methods section MUST state the transport equation and dependent variable; per-axis coefficient convention; mean field and turbulence geometry; spectral tensor and its normalization; component variances and energy conversion; length/radius definitions; particle species and particle variable; the coefficient closure and every free parameter; the parallel dependency or paired backend; interpolation and domain policy; tensor frame and full divergence; stochastic time direction and boundary treatment; numerical tolerances and demonstrated convergence; and the provenance of any fitted preset or data table.

For a running model additionally identify the observable, averaging ensemble, estimator, time window and transition/decorrelation assumption. For a conditional median state its selection explicitly. For frozen lines describe line generation, resolution, retracing and any independently modeled motion relative to a line. An observed longitude width does not uniquely determine $\kappa_\perp$; source extent, field-line geometry and finite-time motion MUST be handled by their selected models.

Do not describe numerical agreement between a transcription and its own integral solver as physical validation. Report source-version discrepancies, fit/calibration uncertainty and any approximations actually used. Only add an observational parameter range after checking its primary measurement, estimator, species and units.

### 18.2 Suggested methods wording

The following is a parameterized drafting aid, not a claim that a simulation was performed:

> We represent spatial diffusion by a symmetric tensor with one parallel and two perpendicular eigenvalues in a specified local magnetic frame. Perpendicular coefficients are defined per transverse axis, with mean free paths given by $\lambda_i=3\kappa_i/v$ in the isotropic diffusion limit. We evaluate the selected coefficient model using its stated magnetic spectral normalization, component variances and scale definitions. Nonlinear closures use the parallel coefficient from the declared parallel model; a coupled NLGCE backend supplies both components from the same state. The transport solver applies the divergence of the full spatially varying tensor and a stochastic covariance equal to twice that tensor. The model identifier, parameter values, source version, numerical backend, domain and convergence evidence are reported with the run.

Replace the generic choices by the actual closure and actual measured validation before publication. If using FLPD, cite (D1)/(D2), specify (F4) or the alternative source-defined field-line input, and use the composite limit (D7) when appropriate. If using NLGCE-F, identify both source tables and describe its state-dependent surrogate qualification. If using a conditional or running diagnostic, rewrite the diffusion-methods sentence accordingly.

### 18.3 Corrections carried from the review and additional source checks

| Earlier issue | Revision 2.0 resolution |
|---|---|
| Running derivative equated to MSD/$2t$ | Separate types and (O1)–(O3), (F8), (N6), (C5) |
| Dimensionally wrong GCD length and missing $\alpha$ | (C2)–(C4); recovered explicit primary prefactor and checked Gaussian integral |
| GCD average presented as an asymptotic path | Scattering-time secant average labeled in (C6) |
| Rounded GCD source approximation treated as accurate | Printed-formula discrepancy recorded; fixture uses (C4) |
| Reduced FLRW factor applied universally | Composite (D7); pure/quasi-2D specialization (D8) |
| Short-path FLPD treated as CLRR automatically | Implicit hybrid (D3), additional limit (D4)/(D5) |
| Claimed $\sqrt2$ typo in retrieved FLPD quadratic | Withdrawn; inspected source has $\sqrt3$ |
| Implicit-slab coefficient confused with composite FLRW | Slab-only $\kappa_s$ required in (U5)–(U9) |
| Area spectrum and reduced spectrum mixed | $g_2=kS_2$ in (S11), with both variance integrals checked |
| Corti preset given independent radial/polar rigidity slopes | Common perpendicular law in (G5)/(G6), as specified in the source |
| Kuhlen asymptotic coefficient evaluated at scattering time | Calibrated decorrelation condition (I13); missing parameters are required inputs |
| Pre-diffusive fit called Parker MSD/injection spread | Uniform-field conditional median; recovered $\beta_m$ and (C11) |
| Fresh field-line kicks on every particle step | Persistent realization, bridge/retracing and divergence counterexample (C9) |
| Forward density and backward generator confused | Explicit (SDE2)–(SDE4) and solver contract |
| RBD restrictions assigned to NLGC; RBD called implicit | Restriction tied to (B2); direct quadrature in (B1) |
| Radial mean free path omitted perpendicular term | Full (T5) distinguished from $\lambda_\parallel^r$ |
| Missing corrected NLGCE $a_x$ and independent parallel pairing | (E2) and one paired owner; all table coefficients reused by index |
| Slab variance repartitioned when adding 2D power | Provider semantics and energy accounting in Section 4.5 |
| Multiplying drift velocity called full reduced drift | Product-rule term (G11) and explicit approximation labeling |
| Divergent moments, arbitrary bounds and isotropic tensor factor | Moment conditions (S8), no universal ceiling, corrected (I15) |
| AMPS defects, measured errors and effort estimates asserted without assets | Removed from evidence; checkout audit is the first roadmap step |

### 18.4 Supported version 2.0 review findings incorporated in 2.1

The attached review was checked against the specified equations, the primary sources and independent calculations. The table records every numbered finding and the separate numerical-rounding correction. It states what was incorporated; unsupported extrapolations in the review are not adopted.

| Review item | Incorporated resolution and limits |
|---|---|
| Numerical-rounding correction | $\pi C(5/3)$ rounds to 0.373417100111093; high-precision computation and explicit in-text audit added. |
| A1 | (D9) gives general and pure-2D branch-consistency parameters; applicable only to the declared closure. |
| A2 | (D10), the short-path fixtures and (D12)/(D13) replace an unsupported $O(A)$ residual assumption with the uncut-tail analysis. |
| A3 | (D11) is an algebraic matching identity; it is not a universal complexity scale. |
| A4 | Section 1.2 states continuous zero-input branches; Section 16.2 gives all bounds, geometry scope and paired-output ownership. |
| A5 | Keep the correct erfc RBD/BC kernel; (B9) and original-author [R28] explain the backtracking step. |
| A6 | Require angular width $h_C$ and a declared unit convention; historical Corti multiplier interpretation remains gated. |
| A7 | Retain quoted $C_K$, distinguish arc length and mean-field coordinate, and prohibit inferred renormalization. |
| A8 | Explicit derivative average (C12), first power of $a_x$ in (E5), and exact RBD prefactor $4\sqrt\pi/45$. |
| B1 | Tag original bend-over and later integral reproductions separately; (B8) gives $\ell_{\perp,\rm int}=L_\perp/2$ for the stated mapping. |
| B2 | State Snodin $L=2\pi/k_0$, peaked $k_b=5k_0$, finite-range $\kappa_0$ variants, and the Kolmogorov-only pair scope. No independent fit amplitude is inferred from field-generation normalization. |
| B3 | Define source isotropic field-line $D_{\rm FL,iso}$, with length units, in $\chi=4D_{\rm FL,iso}/l_c$. |
| B4 | Give the heuristic's required transverse-complexity length and variance in (F7); retain its informational status without a selected definition. |
| B5 | Give (G10) its own registry identity and required declared domain; sampled rigidities do not define an author-calibrated box. |
| B6 | Selectable diagnostics have stable identities and routing; (U10) is separate from the gated full (U11). |
| B7 | Preserve the published composite NLGC slab term; show its mathematical pure-slab limit and prohibit production use at that known failure. No taper or near-boundary threshold is added. |
| B8 | Use the coefficient's coordinate $z$ and require a transport-layer independent/shared line policy; sharing changes correlations, not necessarily the single-particle compound moment. |
| B9 | Show existing paired states A–G and actual two-start/refinement procedure; broader production coverage remains required. |
| B10 | State (S13) as an explicitly supplied normalized form without reliance on withdrawn text. |
| B11 | State $r=r_0=1\,\mathrm{au}$ in the (P6) example. |
| B12 | Include implicit-slab and closed-form consumers in the parallel-dependency list. |
| C1 | Require a physical frame for unequal perpendicular eigenvalues and reject an invalid reference projection. |
| C2 | Represent zero-mean-field isotropy without a preferred axis or division by zero. |
| C3 | Separate observable, estimator and geometry; include covariance, field-line and relative-line statistics and signed Hall results. |
| C4 | Define provenance, numerical diagnostics, domain, quality, identities, ownership and output availability. |
| C5 | Choose success plus a diffusion-regime tag for valid anomalous moments; proved zero limits differ from unavailable coefficients. |
| C6 | Signed table interpolation requires a declared scheme and explicit zero/sign-crossing policy. |
| C7 | Correct configuration ownership: external parallel dependency in the first example; paired backend in the second. |
| D1 | Add the symbol/unit/implementation-name map and separate area index, inverse-rate/source notation, field-line lengths and FLPD parameters. |
| D2 | Define normative keywords, capitalize binding auxiliaries, and keep equations/domains and imperative requirements binding. |

Remaining missing information is recorded as a required input or source gate. In particular this revision assigns no Corti angular width, missing Kuhlen perpendicular calibration, unspecified Hall calibration domain, unsupported Snodin ordered-field variant, line-sharing policy, physical transverse-complexity definition, or AMPS checkout fact by assumption. Synthetic numerical inputs are labeled and are not offered as defaults.

## 19. References and equation provenance

The entries below are the sources used for the actual equations in this revision. References are deliberately scoped; the earlier broad bibliography is not a substitute for a formula audit. Use the specified source version for reproduction. Preprint equation numbering can differ from the journal version.

| Key | Source and inspected material | Equations / role in this specification |
|---|---|---|
| R01 | Potgieter (2013), *Solar Modulation of Cosmic Rays*, [arXiv:1306.4421v1](https://arxiv.org/abs/1306.4421v1), [DOI:10.12942/lrsp-2013-3](https://doi.org/10.12942/lrsp-2013-3). Diffusion/drift discussion and Equations (23)–(26). | Parker equation, normalized modulation family, perpendicular ratios and unresolved angular preset convention |
| R02 | Shalchi (2021), *Perpendicular Diffusion of Energetic Particles: A Complete Analytical Theory*, [arXiv:2109.07574v1](https://arxiv.org/abs/2109.07574v1), [DOI:10.3847/1538-4357/ac2363](https://doi.org/10.3847/1538-4357/ac2363). Definitions, spectra, field-line coefficients and Equations (85)–(98). | (S3)–(S9), (F2)–(F7), (U1), (D1)–(D8); (D9)–(D13) are explicitly derived consistency/tail consequences, not new source closures |
| R03 | Shalchi and Kourakis (2007), *Analytical description of stochastic field-line wandering in magnetic turbulence*, [arXiv:astro-ph/0703366v2](https://arxiv.org/abs/astro-ph/0703366v2), [DOI:10.1063/1.2776905](https://doi.org/10.1063/1.2776905). Field-line correlation closure and flat-spectrum asymptotic. | Field-line normalization, (F1)–(F4), (C3) |
| R04 | Chhiber et al. (2017), *Cosmic-Ray Diffusion Coefficients throughout the Inner Heliosphere from a Global Solar Wind Simulation*, [arXiv:1703.10322v2](https://arxiv.org/abs/1703.10322v2), [DOI:10.3847/1538-4365/aa74d2](https://doi.org/10.3847/1538-4365/aa74d2). Spectra and RBD evaluation. | (S12), area-spectrum normalization and (B3) |
| R05 | Strauss and Fichtner (2015), *On Aspects Pertaining to the Perpendicular Diffusion of Solar Energetic Particles*, [arXiv:1804.03689v1](https://arxiv.org/abs/1804.03689v1), [DOI:10.1088/0004-637X/801/1/29](https://doi.org/10.1088/0004-637X/801/1/29). Pitch-angle prescriptions. | (P3) and averaging convention |
| R06 | Dresing et al. (2012), *The Large Longitudinal Spread of Solar Energetic Particles During the 17 January 2010 Solar Event*, [arXiv:1206.1520v1](https://arxiv.org/abs/1206.1520v1), [DOI:10.1007/s11207-012-0049-y](https://doi.org/10.1007/s11207-012-0049-y). Monte Carlo prescription. | (P5); (P6) is its pitch-angle average |
| R07 | Wijsen et al. (2019), *Modelling three-dimensional transport of solar energetic protons in a corotating interaction region generated with EUHFORIA*, [arXiv:1901.09596v1](https://arxiv.org/abs/1901.09596v1), [DOI:10.1051/0004-6361/201833958](https://doi.org/10.1051/0004-6361/201833958). Perpendicular tensor prescription. | (P7), parallel-contribution radial convention |
| R08 | Matthaeus et al. (2003), *Nonlinear Collisionless Perpendicular Diffusion of Charged Particles*, [DOI:10.1086/376613](https://doi.org/10.1086/376613). Foundational NLGC attribution; source equation also inspected in R09 and R11. | (N1) attribution; source equation checked in R09 and R11 |
| R09 | Shalchi and Hussein (2014), *Perpendicular Diffusion of Energetic Particles in Noisy Reduced Magnetohydrodynamic Turbulence*, [arXiv:1409.2470v1](https://arxiv.org/abs/1409.2470v1), [DOI:10.1088/0004-637X/794/1/56](https://doi.org/10.1088/0004-637X/794/1/56). Equations (1), (3), (4) and NRMHD reduction. | (N1), (U1), (U3)/(U4); explicit normalization checks |
| R10 | Shalchi (2006), *Extended nonlinear guiding center theory of perpendicular diffusion*, [DOI:10.1051/0004-6361:20065465](https://doi.org/10.1051/0004-6361:20065465). Original ENLGC attribution; source equations restated in inspected R12. | (N5)/(N6); (N7) limiting assumptions follow from the integral |
| R11 | Qin and Zhang (2014), *The Modification of the Nonlinear Guiding Center Theory*, [arXiv:1401.1950v2](https://arxiv.org/abs/1401.1950v2), [DOI:10.1088/0004-637X/787/1/12](https://doi.org/10.1088/0004-637X/787/1/12). Equations (1)–(8), Table 2 and coefficient Tables 3/4. | Complete paired model (E1)–(E9), corrected $a_x$, table assets and fit box |
| R12 | Shalchi (2016), *The Implicit Contribution of Slab Modes to the Perpendicular Diffusion Coefficient of Particles Interacting with Two-component Turbulence*, [arXiv:1609.05227v1](https://arxiv.org/abs/1609.05227v1), [DOI:10.3847/0004-637X/830/2/130](https://doi.org/10.3847/0004-637X/830/2/130). Equations (30)–(33), (41)–(43), (55)–(58). | Slab-specific field-line coefficient, exact and rational kernels (U5)–(U9), ENLGC restatement |
| R13 | Shalchi (2015), *Finite gyroradius corrections in the theory of perpendicular diffusion 1. Suppressed velocity diffusion*, [arXiv:1506.07169v1](https://arxiv.org/abs/1506.07169v1), [DOI:10.1016/j.asr.2015.06.018](https://doi.org/10.1016/j.asr.2015.06.018). Suppressed-velocity-diffusion equations and expansion. | (U10)/(U11); source gate for a full harmonic backend |
| R14 | Snodin, Jitsuk, Ruffolo and Matthaeus (2022), *Energetic Particle Perpendicular Diffusion: Simulations and Theory in Noisy Reduced Magnetohydrodynamic Turbulence*, [arXiv:2205.09225v1](https://arxiv.org/abs/2205.09225v1), [DOI:10.3847/1538-4357/ac6e6d](https://doi.org/10.3847/1538-4357/ac6e6d). Equations (18)–(27). | (B1)/(B2), corrected/plain kernels, source integral length (Equation 8), composite closed-form application |
| R15 | Shalchi and Kourakis, *Nonlinear Field Line Random Walk and Generalized Compound Diffusion of Charged Particles*, ICRC 2007 proceedings published 2008, volume 1, pp. 405–408; [primary conference PDF](https://indico.nucleares.unam.mx/event/4/session/83/contribution/45/material/paper/0.pdf). Equations (2)–(6), (8)/(9), visually inspected. Related journal article: *A new theory for perpendicular transport of cosmic rays* (2007), [DOI:10.1051/0004-6361:20077260](https://doi.org/10.1051/0004-6361:20077260). | (C1)–(C6); recovered prefactor, explicit formula versus approximate-number discrepancy; conference version is the inspected primary equation source |
| R16 | Laitinen, Kopp, Effenberger, Dalla and Marsh (2016), *Solar energetic particle access to distant longitudes through turbulent field-line meandering*, [arXiv:1508.03164v2](https://arxiv.org/abs/1508.03164v2), [DOI:10.1051/0004-6361/201527801](https://doi.org/10.1051/0004-6361/201527801). Field-line method and Appendix. | Persistent-line mechanism; line is generated before particle propagation |
| R17 | Laitinen and Dalla (2017), *Energetic particle transport across the mean magnetic field: before diffusion*, [arXiv:1611.05347v1](https://arxiv.org/abs/1611.05347v1), [DOI:10.3847/1538-4357/834/2/127](https://doi.org/10.3847/1538-4357/834/2/127). Return-plane definition, Equation (6) and fit discussion. | Uniform-field conditional median (C10); (C11) derived directly from Equation (6) |
| R18 | Shalchi (2019), *Heuristic Description of Perpendicular Diffusion of Energetic Particles in Astrophysical Plasmas*, [arXiv:1908.00694v1](https://arxiv.org/abs/1908.00694v1), [DOI:10.3847/2041-8213/ab379d](https://doi.org/10.3847/2041-8213/ab379d). Heuristic regime estimates; composite expression also inspected in R14. | Regime interpretation, (F7) and attribution of (B5)/(B6) |
| R19 | Candia and Roulet (2004), *Diffusion and drift of cosmic rays in highly turbulent magnetic fields*, [arXiv:astro-ph/0408054v1](https://arxiv.org/abs/astro-ph/0408054v1), [DOI:10.1088/1475-7516/2004/10/007](https://doi.org/10.1088/1475-7516/2004/10/007). Equations (17)–(19), (22)/(23), Table 1 and domain discussion. | (I2), (I4)/(I5), (G10); exact parameter transcription and radius convention |
| R20 | IAU Resolution B2 (2012), *Re-definition of the astronomical unit of length*, [resolution text hosted by SYRTE, Observatoire de Paris](https://syrte.obspm.fr/IAU_resolutions/IAUResol_2012_0.html). | Exact au definition |
| R21 | BIPM, *The International System of Units (SI Brochure)*, [official publication](https://www.bipm.org/en/publications/si-brochure). | Exact speed of light and SI conversions |
| R22 | Snodin, Shukurov, Sarson, Bushby and Rodrigues (2016), *Global diffusion of cosmic rays in random magnetic fields*, [arXiv:1509.03766v2](https://arxiv.org/abs/1509.03766v2), [DOI:10.1093/mnras/stw217](https://doi.org/10.1093/mnras/stw217). Equations (7)/(8), (16), (20), (23), (26) and fit ranges. | (I3), Kolmogorov (I6)/(I7), shape-specific isotropic fits, peaked $k_b=5k_0$, isotropic field-line coefficient and total-field radius |
| R23 | Kuhlen, Mertsch and Phan (2025), *Diffusion of Relativistic Charged Particles and Field Lines in Isotropic Turbulence. II. Analytical Models*, [arXiv:2211.05882v3](https://arxiv.org/abs/2211.05882v3), [DOI:10.3847/1538-4357/adee94](https://doi.org/10.3847/1538-4357/adee94). Equations (4)/(5), (21)–(24), Table 1 and calibration discussion. | (I8)–(I13); exact running-parallel integral derived here; missing perpendicular calibration remains explicit |
| R24 | Casse, Lemoine and Pelletier (2002), *Transport of cosmic rays in chaotic magnetic fields*, [arXiv:astro-ph/0109223v1](https://arxiv.org/abs/astro-ph/0109223v1), [DOI:10.1103/PhysRevD.65.023002](https://doi.org/10.1103/PhysRevD.65.023002). Equation (19) and adjacent prefactor discussion. | Scaling comparison and source gate for an absolute named fit |
| R25 | Yan and Lazarian (2008), *Cosmic Ray Propagation: Nonlinear Diffusion Parallel and Perpendicular to Mean Magnetic Field*, [arXiv:0710.2617v2](https://arxiv.org/abs/0710.2617v2), [DOI:10.1086/524771](https://doi.org/10.1086/524771). MHD perpendicular transport regimes. | (I16), regime and order-of-magnitude qualification |
| R26 | Corti et al. (2019), *Numerical modeling of galactic cosmic ray proton and helium observed by AMS-02 during the solar maximum of Solar Cycle 24*, [arXiv:1810.09640v3](https://arxiv.org/abs/1810.09640v3), [DOI:10.3847/1538-4357/aafac4](https://doi.org/10.3847/1538-4357/aafac4). Equations (9)–(11) and fixed-parameter discussion. | (G2), (G5), parameterized (G6), unresolved multiplier units in source Equation (11), and checked study-specific constants |
| R27 | Kuhlen, Mertsch and Phan (2025), *Diffusion of Relativistic Charged Particles and Field Lines in Isotropic Turbulence. I. Numerical Simulations*, [arXiv:2211.05881v3](https://arxiv.org/abs/2211.05881v3), [DOI:10.3847/1538-4357/adee9a](https://doi.org/10.3847/1538-4357/adee9a). Reduced-radius definition and Equation (7). | Total-rms-field radius and finite-band correlation convention for (I8) |
| R28 | Ruffolo, Jitsuk, Pianpanit, Snodin, Matthaeus and Chuychai (2015), *Random ballistic interpretation of the nonlinear guiding center theory of perpendicular transport*, PoS(ICRC2015)197, [primary author conference PDF](https://pos.sissa.it/236/197/pdf). Equations (2.6)–(2.10). | Ordinary Gaussian time transform, erfcx, and explicit backtracking factor giving erfc in (B9); magnetostatic specialization |

The standard covariance, product-rule, Gaussian-moment, Brownian-bridge and Itô-generator equations are mathematical consequences of their stated definitions. Source records in the archive distinguish direct PDF inspection from an equation restated in another primary paper. No reference is used to justify an unsupported model parameter.

## 20. Implementation roadmap for Codex in AMPS

### 20.1 Placement and required checkout audit

Place reusable model physics in the existing shared coefficient layer **if confirmed by the current checkout**, preferably beside the parallel coefficient implementation under `src/models/sep_common/`. Keep the 3D mover adapter and tensor/line state in the corresponding SEP application, expected to be `srcSEP3D/`. These are proposed locations based on the supplied description; the exact tree, registry names and build manifests MUST be verified before edits. No directory creation or current-code defect is assumed from this document alone.

Start by reading the checkout's `AGENTS.md` instructions, locating the parallel registry and coefficient types, inspecting the turbulence-provider normalization and tracing the Parker/focused mover call paths. Record the commit, enabled applications, build system, source manifests, test harness and whether the mover represents a forward density or a backward scalar equation. Verify the present perpendicular modes and gradients directly. Repeat any inherited code claim as a new checkout-specific audit rather than citing the earlier file's alleged source inspection.

The physics layer MUST compile without AMPS particle, mesh or MPI dependencies. Reuse established unit, status, configuration and registry types where suitable. Extend the shared paired NLGCE backend instead of creating a second copy in a perpendicular library.

### 20.2 Proposed file ownership

| Proposed component | Responsibility | Host coupling |
|---|---|---|
| `sep_perpendicular_types.h` | Typed inputs/results, observable and domain descriptors | None |
| `sep_perpendicular_registry.{h,cpp}` | Identifiers, parameter schemas, validation and canonical configuration identity | None |
| `sep_perpendicular_spectra.{h,cpp}` | Normalized spectra, moments and tagged length conversions | None |
| `sep_perpendicular_fieldline.{h,cpp}` | Slab/2D/composite coefficients and explicit particle-limit conversions | None |
| `sep_perpendicular_prescribed.{h,cpp}` | Constants, ratios, pitch-angle and parameterized modulation families | None |
| `sep_perpendicular_closures.{h,cpp}` | NLGC, ENLGC, UNLT, implicit-slab, FLPD, RBD and composite relation | None |
| Existing or shared `sep_nlgce_*` backend | One paired NLGCE-N/F implementation used by both libraries | None |
| `sep_perpendicular_running.{h,cpp}` | Compound moments, conditional diagnostic and calibrated Kuhlen algorithm | None; returns statistics, not particle trajectories |
| `sep_perpendicular_fits.{h,cpp}` | Named isotropic fits, classical relation and qualified MHD estimate | None |
| `sep_perpendicular_table.{h,cpp}` | Table schema, loading, interpolation, checksums and validation metadata | None |
| SEP `turbulence/coefficient_bridge.*` | One-time conversion and local provider/dependency assembly | Provider and species coupling |
| SEP `transport/perpendicular_transport.*` | Tensor frame, full divergence, covariance and mover contract | Spatial derivatives and transport |
| SEP field-line state component | Persistent line realization, retracing, checkpoint state | Trajectory/random-stream ownership |
| Existing test locations | Meaningful coefficient, asset and mover tests | Follow checkout harness |

These names are proposals. Keep module boundaries even if the checkout uses different filenames. Reuse an existing checked digest implementation rather than adding a bespoke cryptographic primitive solely for this task. Update actual build manifests and layering tests with each added source.

### 20.3 Ordered implementation phases

**Phase 0 — audit and freeze conventions.** Complete the checkout audit in Section 20.1; reconcile parallel specification 1.3 and this specification; decide mover/observable compatibility; define source, radius, length, energy and angle tags. Resolve or preserve as explicit gates the closed-form length profile, Corti angular convention, Hall domain, Snodin spectral-shape scope and required calibration data; never infer their missing values. Output a concrete convention map and tests for exact unit conversions and existing accepted modes. Preserve or explicitly version any intentional compatibility behavior; do not promise bitwise equivalence without measuring it.

**Phase 1 — types, registry and provider bridge.** Implement Section 16's observable/estimator separation, physical frame kinds, success/regime/availability invariants, numerical diagnostics, quality/domain tags, dependency ownership, parameter validation and complete provenance. Provide stable diagnostic identifiers and bar them from coefficient movers. Add the explicit component-variance and spectral geometry block to the provider adapter. Resolve whether existing wave data are slab-only or total; verify the energy convention before conversion. No 2D fraction is inferred. Establish how paired NLGCE output supplies the parallel result. Tests exercise every missing-input and mismatch path before production formulas are enabled.

**Phase 2 — spectra and inexpensive coefficients.** Implement (S3)–(S12), supported normalized supplied spectra, moment convergence checks, tagged length conversions including (B8), separate bend-over/integral profiles, (F2)–(F5), constants, ratios, power laws and normalized pitch-angle forms. Add Dröge/PARADISE prescriptions with correct averaging and radial definitions. Import both NLGCE polynomial tables unchanged and evaluate the pair, with index and checksum validation. Exit only after these equations and the dimensional references pass; no mover wiring depends on an unimplemented registry entry.

**Phase 3 — direct closures and paired nonlinear solver.** Implement shared dimensionless quadrature and residual diagnostics, then NLGC/ENLGC, UNLT, implicit-slab exact/rational variants, FLPD, explicit RBD and the stable composite closed form. Implement or reuse the corrected paired NLGCE-N integral solver. Preserve distinct spectrum conventions and continuous zero-input branches, enforce every model-specific bound, and retain the published NLGC slab term with its production-geometry gate. Add (D9)/(D10) diagnostics and nonanalytic-tail convergence checks; do not hard-code a guessed asymptotic-switch threshold. Keep erfc for RBD/BC and test (B9) separately. Check dimensional versus reduced integrals, limiting assumptions, root residuals, quadrature refinement and NLGCE multiple starts using explicit A–G inputs plus enumerated application-domain boundaries. Analytic or independently verified implicit sensitivities are added where needed; a flat-spectrum asymptotic exponent is not a universal derivative.

**Phase 4 — tensor and mover integration.** Assemble Cartesian tensors and physical frames. Implement (SDE5)/(SDE6) and compare against (SDE11). Wire the already evaluated local coefficient/pair into the correct Parker or focused operator, with matching covariance and time direction. Audit the complete perpendicular gradient and frame terms. Keep source injection width and explicit drift choices in their correct layers. Check (SDE4), homogeneous covariance, zero eigenvalues, unequal axes, curved-field convergence and snapshot boundaries. Choose engineering timestep controls from the measured convergence evidence, not from an unverified earlier error percentage.

**Phase 5 — tables and efficient evaluation.** Generate tables from the direct backends only after Phase 3 is accepted. Choose axes from the actual application domain; record all conventions and dependencies. Require explicit signed-Hall interpolation, zero/crossing treatment and discrete charge branches. Validate held-out coefficients and gradients at interiors, boundaries, corners and transitions. Enable the selected interpolation only within its tested domain. Add complete cache invalidation; quantify approximations from binning separately. Profile the actual application before choosing table resolution, memory layout or parallel execution. Report measured cost; do not infer runtime from stencil sample count alone.

**Phase 6 — persistent lines and finite-time results.** Implement compound/GCD diagnostics in the statistics API. Add the frozen-line mechanism to the trajectory layer, including the coefficient's signed mean-field coordinate (or a separately supplied arc-length conversion), explicit independent/shared line assignment, reproducible refinement, revisiting and restart. Run the out/back and (F8)/(C9) timestep tests before enabling it for SEP events. Add the Laitinen–Dalla conditional-median descriptor without exposing it as a coefficient. Implement the Kuhlen algorithm only with explicitly supplied calibrated $\gamma_K,L_{c,\perp}$ and root convention; its partial source table alone cannot determine a perpendicular preset.

**Phase 7 — calibrated fits, modulation presets and scientific release.** Implement the Candia–Roulet pair, sharp-cutoff Kolmogorov Snodin pair, the three stated zero-mean-field Snodin $\kappa_0$ variants, and classical relations with their declared domains. Other Snodin ordered-field pairs remain gated. The Hall fit requires its own domain metadata. Add parameterized NWU/Corti families and optional drift coefficients; preserve (G11). A historical preset is admitted only when the source data, sign/angle interpretation, species and epoch are resolved. The $M_A^4$ estimate retains its scaling label and explicit normalization. Publication methods are generated from actual configuration/provenance plus measured validation; do not claim tests that were merely planned.

Phases are dependency ordered: registry/types and providers precede model wiring; direct solves precede generated tables; operator validation precedes production use; persistent-line state precedes field-line-following event models. Steps within a phase may be developed independently only when these interfaces and prerequisites are satisfied.

### 20.4 Source gates for later extensions

An exact Corti preset requires resolution of the angular multiplier convention; no numeric $h_C$ is assigned here. Hall source presets require independently documented domain metadata. Non-Kolmogorov Snodin ordered-field pairs require complete variant-specific definitions and calibration; their isotropic $\kappa_0$ fits do not close that gate.

For a general particle QLT backend, full finite-gyroradius harmonic model, WNLT/SOQLT, dynamic decorrelation, compressive turbulence, Casse absolute fit, restricted-scattering amplitudes, a specific HelMod polar preset or an unresolved drift fit, first produce a source-checked equation/convention/domain record. Required formulas include all resonance factors, harmonic sums and truncation controls where relevant. Reuse neither a dimensionally similar formula nor an inferred parameter as a substitute. Implement a separately versioned model only after the gate is closed.

The remaining journal-specific GCD check is an equation-version comparison; the fully specified Gaussian-convolution backend already follows the inspected conference equations and independent algebra. A user-requested reproduction of the source's inconsistent rounded prefactor MUST be separately identified as such and cannot overwrite the verified formula.

### 20.5 Definition of done

A model is ready for a scientific release when its exact source formula, units, independent inputs, observable and domain are documented; all source conventions and calibration inputs are supplied; zero, missing and diagnostic outputs have distinct tested semantics; normalization and limiting tests pass; an independent mathematical evaluation agrees within the stated numerical budget; derivatives required by the mover are verified; every failure path is tested; and configuration/provenance includes every relevant parameter and asset checksum. A model consumed by a mover additionally needs operator, timestep, boundary and statistical verification in that mover.

Deliver the shared library, adapters, explicit configuration examples, table assets, meaningful tests, source audit and publication-ready methods description tied to the implemented commit. Record which backends remain gated. Do not estimate scientific completion from lines of code, assert current-tree locations without inspection, or report unrun AMPS checks as passed.
