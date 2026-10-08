# Parallel diffusion coefficient models for SEP and GCR transport

**File:** PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md  
**Specification version:** 1.3  
**Revision date:** 8 October 2026  
**Literature cutoff:** 7 October 2026; principal review interval: 2006–2026.  
**Purpose:** a mathematical reference for implementing a reusable library of parallel spatial diffusion coefficients, and a source for the methods sections of publications.

This document describes coefficients for solar energetic particles (SEPs) and galactic cosmic rays (GCRs) transported with the Parker equation. It covers the seven families identified in the preceding review: constant mean free paths, SEP power laws, empirical broken rigidity laws, quasi-linear scattering, weakly nonlinear/resonance-broadened scattering, nonlinear parallel closures, and coefficients evaluated from transported turbulence. It also specifies a Bohm comparison model and adapters for supplied wave spectra and coefficient tables.

The library computes a **local transport coefficient**. It does not determine particle injection, shock geometry, magnetic connectivity, solar-wind evolution, or the turbulence evolution by itself. The same mathematical coefficient can be evaluated in a Parker solver or used to normalize a pitch-angle operator in a focused-transport solver, but these are different transport descriptions.

Published theory, algebraic consequences of the stated conventions, and proposed software requirements are distinguished throughout. Statements headed **Implementation requirement** are design decisions for the future library. Illustrative configurations and verification points are mathematical examples, not calibrated physical predictions.

**Numerical qualification:** NLGCE-F is the published polynomial approximation to the corrected NLGCE-N integral closure. Revision 1.2 corrects a transcription error in Equation (45): the source multiplies by 1/epsilon; earlier versions incorrectly exponentiated. All affected nonlinear benchmarks and discrepancy statistics have been regenerated. The 576 polynomial coefficients are unchanged and individually match the primary-source tables. Membership in the published input box remains distinct from a certified approximation-error bound.

The document remains a single implementation and publication reference. Sections 17–18 supply publication reporting and citations; Section 19 is a dated bibliography-snapshot appendix; Sections 20–22 supply coefficient data, BibTeX additions, and reproducibility assets. Section 23 gives the staged Codex implementation roadmap. The snapshot audit is not a physical assumption of any model.

## Contents

1. [Transport equation and scope](#1-transport-equation-and-scope)
2. [Variables, units, and validity](#2-variables-units-and-validity)
3. [Model inventory](#3-model-inventory)
4. [Constant coefficients](#4-constant-coefficients)
5. [SEP power-law models](#5-sep-power-law-models)
6. [Smooth broken rigidity models](#6-smooth-broken-rigidity-models)
7. [Pitch-angle conversion and normalization](#7-pitch-angle-conversion-and-normalization)
8. [Slab quasi-linear theory](#8-slab-quasi-linear-theory)
9. [Resonance broadening and weakly nonlinear theory](#9-resonance-broadening-and-weakly-nonlinear-theory)
10. [NLPA and coupled nonlinear models](#10-nlpa-and-coupled-nonlinear-models)
11. [NLGCE-F polynomial model](#11-nlgce-f-polynomial-model)
12. [Turbulence and wave-spectrum adapters](#12-turbulence-and-wave-spectrum-adapters)
13. [Bohm and supplied-table models](#13-bohm-and-supplied-table-models)
14. [Library interface and numerical rules](#14-library-interface-and-numerical-rules)
15. [Verification requirements and numerical references](#15-verification-requirements-and-numerical-references)
16. [Model selection and scientific interpretation](#16-model-selection-and-scientific-interpretation)
17. [Publication methods text](#17-publication-methods-text)
18. [References](#18-references)
19. [Bibliography comparison](#19-bibliography-comparison)
20. [NLGCE-F coefficient tables](#20-nlgce-f-coefficient-tables)
21. [BibTeX additions](#21-bibtex-additions)
22. [Reproducibility assets and revision record](#22-reproducibility-assets-and-revision-record)
23. [Codex implementation roadmap](#23-codex-implementation-roadmap)

## 1. Transport equation and scope

### 1.1 Parker equation

Let f(x,p,t) be the nearly isotropic phase-space distribution in the local plasma frame. The transport equation is

$$
\frac{\partial f}{\partial t}
=
\nabla\cdot\left(\mathbf K_s\cdot\nabla f\right)
-\left(\mathbf U+\mathbf v_d\right)\cdot\nabla f
+\frac{1}{3}\left(\nabla\cdot\mathbf U\right)
p\frac{\partial f}{\partial p}
+Q.
\tag{1}
$$

Here U is the convection velocity used by the selected transport formulation, v_d represents gradient, curvature, and current-sheet drifts, and Q is a source. Equation (1) omits momentum diffusion and other optional terms. Its classical foundation is Parker (1965) [R01]. The coefficient library supplies the parallel eigenvalue of the symmetric spatial diffusion tensor,

$$
\mathbf K_s
=
\kappa_\parallel\,\mathbf b\mathbf b
+\kappa_{\perp 1}\,\mathbf e_1\mathbf e_1
+\kappa_{\perp 2}\,\mathbf e_2\mathbf e_2,
\qquad
\mathbf b=\frac{\mathbf B_0}{B_0}.
\tag{2}
$$

The vectors (b,e1,e2) form a local orthonormal basis. The drift contribution is handled separately; it must not be added to the positive symmetric diffusion tensor as a third spatial diffusion process.

For an axisymmetric perpendicular coefficient, Equation (2) reduces to

$$
\mathbf K_s
=
\kappa_\perp\mathbf I+
(\kappa_\parallel-\kappa_\perp)\mathbf b\mathbf b.
\tag{3}
$$

Perpendicular diffusion enters this document where it is needed to define a nonlinear parallel closure. Otherwise its prescription belongs to a separate model.

### 1.2 Parallel, radial, and shock-normal coefficients

The transport mean free path is defined by

$$
\boxed{\kappa_\parallel=\frac{v\lambda_\parallel}{3}}.
\tag{4}
$$

It is a scattering measure, not the distance to the next large-angle collision unless a particular collision model establishes that interpretation.

In a Parker spiral with angle psi between the mean field and the radial direction,

$$
\kappa_{rr}
=
\kappa_\parallel\cos^2\psi+
\kappa_{\perp,r}\sin^2\psi,
\qquad
\lambda_r=\frac{3\kappa_{rr}}{v}.
\tag{5}
$$

For a unit shock normal n and axisymmetric perpendicular diffusion,

$$
\kappa_n
=
\mathbf n\cdot\mathbf K_s\cdot\mathbf n
=
\kappa_\parallel\cos^2\theta_{Bn}
+\kappa_\perp\sin^2\theta_{Bn}.
\tag{6}
$$

Consequently, a constant radial mean free path and a constant parallel mean free path are different prescriptions. A paper or input file must identify which is used. A parallel coefficient should be evaluated first and then projected through the tensor. Potgieter et al. (2014) explicitly distinguish these projections [R14].

**Implementation requirement:** no model accepts an unlabeled “mean free path.” A radial input requires a geometry adapter and a specified perpendicular coefficient; it must not be interpreted automatically as lambda_parallel.

### 1.3 Momentum changes and acceleration scope

The compression term in Equation (1), together with spatial transport across a shock, can represent first-order diffusive shock acceleration without momentum diffusion. A model that also includes stochastic acceleration or reacceleration requires additional scattering physics. A commonly used isotropic momentum-diffusion operator is

$$
\left.\frac{\partial f}{\partial t}\right|_{D_{pp}}
=\frac{1}{p^2}\frac{\partial}{\partial p}
\left(p^2D_{pp}\frac{\partial f}{\partial p}\right),
\qquad [D_{pp}]={\rm (kg\,m\,s^{-1})^2\,s^{-1}}.
\tag{1a}
$$

Here f is the same phase-space density as in Equation (1); an equation evolved for a momentum-weighted distribution requires the corresponding transformed operator. Marcowith and Kirk (1999) distinguish first-order shock acceleration from the optional second-order terms [R31].

No universal D_pp can be inferred from kappa_parallel alone. A compatible extension must specify scattering-center velocities, directional wave populations, the scattering frame, and any mixed pitch-angle/momentum coefficients. A relation derived for balanced Alfvén waves must not be reused for arbitrary imbalance or static turbulence. A future momentum-diffusion module is a separate extension of this spatial-coefficient library.

## 2. Variables, units, and validity

### 2.1 Particle variables

Use the total particle rest mass m, signed charge q=Ze, total kinetic energy T, and momentum magnitude p. In SI units,

$$
pc=\sqrt{T(T+2mc^2)},\qquad
\gamma=1+\frac{T}{mc^2},
\qquad
v=\frac{pc^2}{T+mc^2},\qquad
\beta=\frac{v}{c}.
\tag{7}
$$

For ions, energy per nucleon must be converted to total particle energy before Equation (7) is evaluated. The charge state Z is supplied independently of the mass number A.

The positive rigidity, relativistic gyrofrequency magnitude, and maximum gyroradius are

$$
\mathcal R=\frac{pc}{|q|},\qquad
\Omega=\frac{|q|B_0}{\gamma m},\qquad
r_L=\frac{p}{|q|B_0}
=\frac{\mathcal R}{cB_0}
=\frac{v}{\Omega}.
\tag{8}
$$

The maximum gyroradius means the radius at a pitch angle of 90 degrees; the actual perpendicular gyroradius is r_L sqrt(1-mu^2). Rigidity is in volts, not momentum units. Gyrofrequency is in radians per second, and wavenumber k below is in radians per metre.

Charge sign is needed for drift and polarized/directional wave interactions. The balanced, magnetostatic scalar models in this document depend on |q|.

### 2.2 Required unit conventions

| Quantity | Symbol | Internal SI unit |
|---|---|---|
| Parallel spatial diffusion coefficient | kappa_parallel | m² s⁻¹ |
| Parallel mean free path | lambda_parallel | m |
| Pitch-angle diffusion coefficient | D_mu_mu | s⁻¹ |
| Mean magnetic-field magnitude | B0 | T |
| Magnetic fluctuation variance | delta B² | T² |
| Spectral bend-over lengths | ell_s, ell_2 | m |
| Slab correlation length in NLPA/NLGCE | L_c,s; abbreviated L_c in Equations (28),(45) | m; fixed by Equation (28) |
| Other provider correlation length | Provider-specific named length | m; covariance, direction, and integration convention required |
| Particle rigidity | mathcal R | V |
| Momentum magnitude | p | kg m s⁻¹ |
| One-sided magnetic wavenumber spectrum | P(k) | T² m |
| Spectral angular frequency and width | omega, Gamma | s⁻¹ |

The astronomical-unit value follows [IAU Resolution B2 (2012)](https://iauarchive.eso.org/static/resolutions/IAU2012_English.pdf). SI prefixes and unit powers determine the remaining conversions; see the [BIPM SI Brochure](https://www.bipm.org/en/publications/si-brochure). The exact conversions are

$$
1\ {\rm AU}=149\,597\,870\,700\ {\rm m},\quad
1\ {\rm GV}=10^9\ {\rm V},\quad
1\ {\rm nT}=10^{-9}\ {\rm T},\quad
1\ {\rm cm^2\,s^{-1}}=10^{-4}\ {\rm m^2\,s^{-1}}.
\tag{9}
$$

A published normalization of 10^22 cm² s⁻¹ therefore equals 10^18 m² s⁻¹. A variance in nT² requires a factor of 10^-18 to obtain T².

**Implementation requirement:** input conversion occurs once at the boundary of the library. Stored configuration values and evaluation results carry explicit units; an unexplained bare number is never assumed to be an AU, GV, or nT value.

### 2.3 Local diffusion limit

The usual conversion between D_mu_mu and kappa_parallel assumes rapid pitch-angle relaxation compared with changes in the distribution and background. Useful diagnostics are

$$
\tau_{\rm sc}\sim\frac{\lambda_\parallel}{v},
\qquad
{\rm Kn}_f=\frac{\lambda_\parallel}{L_f},
\qquad
{\rm Kn}_B=\frac{\lambda_\parallel}{|L_B|},
\qquad
L_B=-\left(\mathbf b\cdot\nabla\ln B_0\right)^{-1}.
\tag{10}
$$

L_f is an estimated parallel gradient length of the distribution. The local, nearly isotropic approximation is best supported when these Knudsen numbers are small and the background changes on times much longer than tau_sc. They are diagnostics, not universal sharp thresholds.

Early SEP arrivals, weak scattering, strong focusing near the Sun, and some shock precursors can require focused transport. A positive finite coefficient does not by itself establish the validity of a Parker description. Shalchi and Klippenstein (2025) further show that parallel diffusion in a nonuniform field depends on the precise transport definition [R26]; a focusing correction cannot be inserted as a universal multiplier.

**Implementation requirement:** validity flags are separate from numerical success. In particular, a small integration residual must not clear a diffusion-limit warning.

For models using a mean field and a magnetic variance, define the total perturbation level

$$
\epsilon^2=\frac{\delta B_s^2+\delta B_2^2}{B_0^2}.
$$

Return a strong-turbulence diagnostic when epsilon>=1. This is an explicit operating diagnostic, not a universal rejection threshold: QLT requires a weak-perturbation justification that may fail below that threshold, while the published NLGCE-F fitting box includes stronger turbulence. The field-aligned tensor remains mathematically defined for B0>0; its physical interpretation depends on the chosen mean-field scale and transport closure. Record that mean-field definition. B0=0 is a distinct failure of the decomposition.

For NLGCE-F, distinguish an unbounded surrogate-error diagnostic from numerical evaluation success. The corrected audit in Section 11.4 does not support the former epsilon²<0.1 warning rule, which was based on a mistranscribed a_x. No turbulence-only accuracy threshold is imposed. Application-specific error assessment must use the actual four-dimensional state and the correct closure.

### 2.4 Symbol and implementation-name map

The following map disambiguates symbols that otherwise have similar typography. A code field must identify the physical quantity rather than reproduce an ambiguous single-letter name.

| Symbol | Meaning and scope | Suggested implementation name |
|---|---|---|
| mathcal R | Positive particle rigidity, Equation (8) | rigidity_V |
| r_* = r_L/ell_s | Dimensionless gyroradius used in the QLT expressions | gyroradius_to_slab_scale |
| r, x_vector | Heliocentric radius and spatial position | radius_m, position_m |
| x_1,...,x_4 | Natural-log inputs of NLGCE-F | fit_log_inputs |
| x=k ell_a | Dimensionless spectral integration variable | dimensionless_wavenumber |
| s, nu=s/2 | Magnetic inertial index and smooth-spectrum shape index | inertial_index, spectrum_nu |
| q_mu | Prescribed pitch-angle shape index, Equations (22)–(24) | pitch_angle_shape_index |
| s_arc | Field-line arc length, Equation (63) | arclength_m |
| A_tube | Flux-tube cross-sectional area | tube_area_m2 |
| A(k), A_s(k), A_2(k) | NLPA decorrelation denominators, with units s^-1 | nlpa_rate_s_inv |
| Gamma(nu) | Euler gamma function | gamma_function |
| Gamma_dec(k) | Temporal decorrelation rate in the resonance and NLPA kernels | decorrelation_rate_s_inv |
| mathscr R_k(omega) | Resonance function, with units s | resonance_kernel_s |
| Delta(k) | Gaussian frequency-width parameter | gaussian_width_s_inv |
| a, b, h | Low/high slopes and smoothness in Equation (17) only | low_slope, high_slope, smoothness |
| a in Equation (14) | Kinetic-energy exponent | kinetic_energy_exponent |
| a_x, a², a'^2 | Distinct NLPA/NLGC closure parameters | nlpa_ax, nlgc_a_squared, modified_a_squared |
| h_mu | Additional scattering near mu=0 | scattering_90deg_parameter |
| C(nu), C_k(t) | Spectrum normalization constant and temporal decorrelation function | spectrum_normalization, temporal_correlation |
| mu, mu_0 | Pitch-angle cosine and vacuum permeability | pitch_angle_cosine, vacuum_permeability |
| Z, Z² | Ion charge state and separately defined turbulence moment | charge_state, turbulence_moment_Z2 |
| sigma_D, f_s | Normalized residual energy and slab magnetic-variance fraction | residual_energy_fraction, slab_fraction |

The subscripted decorrelation rate Gamma_dec is distinct from the gamma function. The flux-tube area and arc length are distinct from the NLPA rate and spectral index. In prose, a plain rigidity abbreviation R always refers to mathcal R; the spectral ratio is r_*.

## 3. Model inventory

| Family | Proposed stable identifier | Direct output | Additional physical input |
|---|---|---|---|
| Constant parallel mean free path | constant_lambda | lambda, kappa | Specified positive lambda0 |
| Constant diffusion coefficient | constant_kappa | kappa, lambda | Specified positive kappa0 |
| SEP separable power law | power_law_lambda | lambda, kappa | Reference scales and exponents |
| Empirical smooth broken rigidity law | broken_rigidity_kappa | kappa, lambda | Amplitude, slopes, break, field scaling |
| Slab QLT, supplied or analytic spectrum | qlt_slab_spectrum | D_mu_mu, lambda, kappa | Slab magnetic spectrum and B0 |
| Slab QLT, inertial-range approximation | qlt_slab_inertial | lambda, kappa | Slab variance, bend-over length, slope |
| Prescribed mean free path with regularized pitch-angle shape | prescribed_lambda_mu_shape | D_mu_mu, lambda, kappa | Target lambda, slope, 90-degree parameter |
| Slab resonance broadening | broadened_slab | D_mu_mu, lambda, kappa | Spectrum and explicit broadening kernel |
| NLPA with externally supplied perpendicular diffusion | nlpa_given_perp | kappa_parallel, lambda_parallel | Slab/2D spectra, kappa_perp |
| Coupled original NLGC–NLPA | nlgc_e | Parallel and perpendicular pair | Slab/2D state |
| Coupled modified NLGC–NLPA | nlgce_n | Parallel and perpendicular pair | Slab/2D state |
| Polynomial fit of the preceding closure | nlgce_f_2014 | Parallel and perpendicular pair | Four dimensionless variables |
| Transported turbulence | turbulence_adapter | Output of selected closure | Supplied turbulence moments |
| Self-generated or prescribed wave spectrum | wave_spectrum_adapter | Output of selected scattering closure | Supplied wave spectrum and conventions |
| Bohm comparison | bohm | lambda, kappa | B0 and eta_B |
| Supplied numerical coefficient grid | tabulated_parallel | lambda or kappa | Validated grid, axes, units |

A turbulence adapter and a wave-spectrum adapter are **providers plus closures**, not new microscopic scattering laws. The library must record the underlying closure as part of their model identity.

“NLGC” and “UNLT” alone name perpendicular transport theories. They do not determine a parallel coefficient without additional input. “NLGCE-F” names a particular fit, not a generic polynomial regression. Matthaeus et al. (2003) [R05], Qin (2007, corrected 2013) [R08,R09], and Qin and Zhang (2014) [R10] establish these distinctions.

### 3.1 Implemented scope and paper-specific extensions

The displayed formulas define the proposed v1 library scope; they do not imply that every published nonlinear scattering theory is included.

| Capability | Scope in this specification | Requirements for a later extension |
|---|---|---|
| Regularized SEP pitch-angle prescription | Explicitly supplied in Equations (22),(23), with Dröge-style attribution [R30] | Report shape, amplitude/target-lambda mode, and spatial/energy dependence |
| Second-order QLT (SOQLT) | Literature context [R27]; no complete paper-specific evaluator specified | Select the exact theory and orbit-broadening approximation, normalize its spectral tensor, and validate its 90-degree and limiting behavior |
| Complete composite WNLT | Literature context [R06]; Section 9 supplies explicitly defined slab kernels | Specify the coupled orbit correlations, tensor geometry, and perpendicular dynamics of the chosen published theory |
| Directional/dynamical wave scattering and D_pp | Provider requirements and acceleration scope only | Supply compatible wave-frame and cross-diffusion physics; use a separate model identifier |
| Non-axisymmetric perpendicular tensor | Permitted in Equation (2), assembled by the transport solver | Supply both perpendicular coefficients and the spatial derivatives of all basis vectors |

Adding a Gaussian resonance kernel does not by itself implement SOQLT, and renaming a scalar coefficient does not implement a focused-transport solver.

## 4. Constant coefficients

### 4.1 Constant parallel mean free path

For a supplied lambda0>0,

$$
\lambda_\parallel=\lambda_0,
\qquad
\kappa_\parallel=\frac{v\lambda_0}{3}.
\tag{11}
$$

A constant mean free path produces a speed-dependent diffusion coefficient. Nonrelativistically, kappa is proportional to sqrt(T); relativistically it approaches c lambda0/3.

This model is useful for controlled comparisons and for separating changes in particle speed from assumed changes in scattering. It does not infer turbulence strength or a universal SEP scattering length.

### 4.2 Constant parallel diffusion coefficient

For a supplied kappa0>0,

$$
\kappa_\parallel=\kappa_0,
\qquad
\lambda_\parallel=\frac{3\kappa_0}{v}.
\tag{12}
$$

Here the inferred mean free path increases at low particle speed. Equations (11) and (12) therefore must have different identifiers and different input fields.

**Implementation requirement:** no position or time dependence is added to a constant model. A time-varying normalization is a separately configured provider. Spatial derivatives vanish when the particle momentum is held fixed.

## 5. SEP power-law models

### 5.1 Rigidity power law

A fully normalized separable prescription is

$$
\lambda_\parallel(\mathbf x,\mathcal R,t)
=
\lambda_0
\left(\frac{r}{r_0}\right)^\alpha
\left(\frac{\mathcal R}{\mathcal R_0}\right)^\delta
\left(\frac{B_0}{B_{\rm ref}}\right)^{-\eta}
g_t(t)\,g_{\rm region}(\mathbf x,t),
\qquad
\kappa_\parallel=\frac{v\lambda_\parallel}{3}.
\tag{13}
$$

The radial, field, time, and region factors are optional, and each equals one when disabled. lambda0 refers to the state where every ratio and factor is one.

| Parameter | Meaning | Requirement |
|---|---|---|
| lambda0 | Reference parallel mean free path | Positive |
| r0, R0, B_ref | Normalization scales | Positive and explicitly dimensioned |
| alpha | Explicit radial exponent | Real, justified for the chosen domain |
| delta | Mean-free-path rigidity exponent | Real |
| eta | Explicit mean-field exponent | Real; zero disables this factor |
| g_t, g_region | Supplied dimensionless factors | Positive and identified |

The frequently used low-rigidity value delta=1/3 corresponds to the inertial-range scaling of magnetostatic Kolmogorov slab QLT under fixed turbulence quantities. It is not a universal empirical result at every energy or heliocentric distance.

The factors in Equation (13) are an empirical prescription. Including both r^alpha and B0^-eta is allowed only when their combined dependence is intended; the Parker field already varies with radius. The same caution applies when a separately evolved turbulence level is used.

### 5.2 Energy and speed variants

If reproducing a model formulated in kinetic energy, implement a distinct option,

$$
\lambda_\parallel
=
\lambda_0
\left(\frac{T}{T_0}\right)^a
\left(\frac{r}{r_0}\right)^\alpha.
\tag{14}
$$

A speed power law is likewise distinct,

$$
\lambda_\parallel
=
\lambda_0
\left(\frac{v}{v_0}\right)^a_v
\left(\frac{r}{r_0}\right)^\alpha.
\tag{15}
$$

For a fixed species in the nonrelativistic limit, T is proportional to R², so an energy exponent a corresponds approximately to a rigidity exponent 2a. In the ultrarelativistic limit, T is approximately proportional to R. Equation (14) therefore cannot be converted to a single rigidity exponent across both regimes.

At fixed background, Equation (13) gives

$$
\frac{\partial\ln\lambda_\parallel}{\partial\ln\mathcal R}=\delta,
\qquad
\frac{\partial\ln\kappa_\parallel}{\partial\ln\mathcal R}
=
\delta+\frac{1}{\gamma^2}.
\tag{16}
$$

This speed contribution is important for low-energy ions.

**Implementation requirement:** require an explicit independent variable: rigidity, total kinetic energy, energy per nucleon, or speed. Convert species variables before evaluating the power law.

## 6. Smooth broken rigidity models

### 6.1 Dimensionally normalized formulation

A common GCR prescription uses a smooth transition between two rigidity slopes. An implementation form equivalent to the speed-factored, reference-normalized form of Potgieter et al. (2014) is

$$
\begin{aligned}
H(\mathcal R)
&=
\left(\frac{\mathcal R}{\mathcal R_0}\right)^a
\left[
\frac{(\mathcal R/\mathcal R_0)^h+
      (\mathcal R_b/\mathcal R_0)^h}
     {1+(\mathcal R_b/\mathcal R_0)^h}
\right]^{(b-a)/h},\\
\kappa_\parallel
&=
K_\star(t)\,\beta
\left(\frac{B_{\rm ref}}{B_0}\right)^\eta
H(\mathcal R)\,g_r(r)\,g_{\rm region}(\mathbf x,t).
\end{aligned}
\tag{17}
$$

K_star has units of diffusion coefficient, a and b are the low- and high-rigidity slopes after removing beta, R_b is the transition rigidity, and h>0 controls transition sharpness. The conventional field factor has eta=1; allowing other eta values is an explicitly configured generalization.

Because H(R0)=1, the actual coefficient at the reference field, position, and rigidity is K_star beta(R0), not K_star. With the displayed beta factor,

$$
\lambda_\parallel
=
\frac{3K_\star(t)}{c}
\left(\frac{B_{\rm ref}}{B_0}\right)^\eta
H(\mathcal R)\,g_r\,g_{\rm region}.
\tag{18}
$$

This convention allows the same mean free path for different species at equal rigidity while retaining their different speeds in kappa. If a user instead supplies the actual kappa at R0, the input converter must divide that value by beta(R0) to obtain K_star.

Potgieter et al. (2014), Vos and Potgieter (2015), and Corti et al. (2019) use this class of phenomenological coefficients in Parker-equation modulation studies [R14,R15,R18]. Their fitted parameters and normalization conventions are study-specific.

### 6.2 Limits and derivatives

For R much smaller than R_b, H is proportional to R^a; for R much larger than R_b, it is proportional to R^b. If a=b, the model is exactly a single power law.

Define w=[1+(R_b/R)^h]^-1. Then

$$
\frac{\partial\ln\lambda_\parallel}{\partial\ln\mathcal R}
=
a+(b-a)w,
\qquad
\frac{\partial\ln\kappa_\parallel}{\partial\ln\mathcal R}
=
a+(b-a)w+\frac{1}{\gamma^2}.
\tag{19}
$$

Thus the full low-rigidity kappa slope includes beta. At a nonrelativistic fixed species it approaches a+1, not a.

**Implementation requirement:** compute log(H) with a log-sum-exp evaluation of the terms R^h and R_b^h. This avoids overflow and cancellation at extreme rigidity ratios. Require K_star>0, R_b>0, R0>0, and h>0. No solar-cycle parameter table is assumed.

### 6.3 Time dependence and transient reductions

K_star, the slopes, and R_b can be supplied as functions of time. A local time sequence is not automatically a physically propagated heliospheric history. A transport/background provider must specify whether a value is simultaneous, convected, retarded, or part of a steady-state sequence.

A positive empirical reduction factor can represent stronger scattering in a sheath or other region, but it must be named and documented. Do not describe such a factor as self-generated turbulence unless a wave-growth calculation supplies it.

## 7. Pitch-angle conversion and normalization

### 7.1 Diffusion-limit relation

For a gyrotropic pitch-angle operator of the form d/dmu[D_mu_mu df/dmu], the homogeneous, weak-focusing diffusion limit is

$$
\boxed{
\lambda_\parallel
=
\frac{3v}{8}
\int_{-1}^{1}
\frac{(1-\mu^2)^2}{D_{\mu\mu}(\mu)}\,d\mu,
\qquad
\kappa_\parallel
=
\frac{v^2}{8}
\int_{-1}^{1}
\frac{(1-\mu^2)^2}{D_{\mu\mu}(\mu)}\,d\mu.
}
\tag{20}
$$

These relations connect scattering calculations with the Parker coefficient; see the QLT foundations [R02,R03,R04]. mu is the cosine of the particle pitch angle relative to the local mean field.

For even D_mu_mu, integrate over 0<mu<1 and multiply by two. Do not apply this symmetry to a polarized or directional wave model without checking it.

### 7.2 Exact isotropic-scattering check

For

$$
D_{\mu\mu}=D_0(1-\mu^2),
$$

Equation (20) gives

$$
\lambda_\parallel=\frac{v}{2D_0},
\qquad
\kappa_\parallel=\frac{v^2}{6D_0}.
\tag{21}
$$

This is an exact normalization check, independent of any turbulence spectrum.

### 7.3 Regularized pitch-angle shape

A common practical pitch-angle shape can be specified as

$$
D_{\mu\mu}
=
D_0(1-\mu^2)\left(|\mu|^{q_\mu-1}+h_\mu\right),
\qquad h_\mu\geq0.
\tag{22}
$$

Here q_mu is a scattering-shape spectral index and h_mu represents additional scattering near 90 degrees. The regularized shape is used in SEP focused-transport work, including Dröge et al. (2014) [R30]; q_mu and h_mu correspond to their shape index and 90-degree parameter, with amplitude differences absorbed in D0. This is a phenomenological prescription, not a complete nonlinear turbulence theory. The amplitude may depend on position and particle energy through a separately specified provider.

To enforce a supplied target mean free path, define

$$
I(q_\mu,h_\mu)
=
\int_{-1}^{1}
\frac{1-\mu^2}{|\mu|^{q_\mu-1}+h_\mu}\,d\mu,
\qquad
D_0=\frac{3v}{8\lambda_{\rm target}}I(q_\mu,h_\mu).
\tag{23}
$$

Then Equation (20) returns lambda_target exactly. Holding D0 fixed while changing h_mu changes the mean free path. The target-lambda and fixed-amplitude modes must therefore be distinct.

For h_mu=0 and 1<q_mu<2,

$$
I(q_\mu,0)=\frac{4}{(2-q_\mu)(4-q_\mu)}.
\tag{24}
$$

A divergent integral is a physical failure of the selected local diffusion closure, not a quadrature error to hide with an arbitrary mu cutoff.

**Implementation requirement:** at mu=±1, use the limiting integrand or open-interval quadrature rather than evaluating a literal 0/0. For the shape in Equation (22), integrate the reduced expression in Equation (23). A quadrature endpoint rule must not introduce a scattering floor.

**Implementation requirement:** a Parker solver uses kappa derived from the pitch-angle operator; a focused solver evolves D_mu_mu. Do not add a separate parallel spatial diffusion operator to a focused solver already representing the same scattering through its pitch-angle dynamics.

## 8. Slab quasi-linear theory

### 8.1 Assumptions and spectrum convention

The fully specified baseline here assumes:

- a locally uniform mean field B0;
- magnetostatic, transverse, axisymmetric slab fluctuations;
- no magnetic helicity or directional wave imbalance;
- weak perturbations for the interpretation as QLT;
- particle speed large enough that neglected wave frequencies are unimportant.

Let P_s(k), for k>=0, be the **one-sided spectrum of the total transverse slab magnetic variance**,

$$
\int_0^\infty P_s(k)\,dk=\delta B_s^2,
\qquad
\langle\delta B_x^2\rangle
=
\langle\delta B_y^2\rangle
=
\frac{\delta B_s^2}{2}.
\tag{25}
$$

Equivalently, the signed one-dimensional spectrum for one transverse component is S_xx(k)=P_s(|k|)/4. These definitions fix otherwise ambiguous factors of two and pi.

The gyroresonant wavenumber and pitch-angle coefficient are

$$
k_{\rm res}=\frac{\Omega}{v|\mu|}
=\frac{1}{r_L|\mu|},
\qquad
\boxed{
D_{\mu\mu}^{\rm QLT}
=
\frac{\pi\Omega^2}{4B_0^2v|\mu|}
(1-\mu^2)P_s(k_{\rm res}).
}
\tag{26}
$$

Equation (26) uses the exact spectrum convention in Equation (25). A formula copied from a paper with a per-component, two-sided, or differently normalized tensor spectrum requires conversion before use. The factor cannot be treated as an adjustable calibration constant.

The resonant QLT basis is Jokipii (1966), with heliospheric spectral and dynamical refinements discussed by Bieber et al. (1994), Teufel and Schlickeiser (2003), and Zank et al. (1998) [R02,R03,R04,R07].

### 8.2 Smooth bend-over spectrum

An analytic spectrum convenient for both QLT and the nonlinear closures is

$$
\begin{aligned}
P_s(k)
&=
4C(\nu)\,\ell_s\,\delta B_s^2
\left[1+(k\ell_s)^2\right]^{-\nu},\\
\nu&=\frac{s}{2},
\qquad
C(\nu)=
\frac{\Gamma(\nu)}
{2\sqrt{\pi}\,\Gamma(\nu-\tfrac12)},
\qquad s>1.
\end{aligned}
\tag{27}
$$

The integral in Equation (25) is exactly delta B_s². ell_s is the spectral bend-over length. It is not automatically an integral correlation length.

For the integral correlation-length definition based on the normalized even spatial covariance,

$$
L_c=
\frac{1}{\delta B_s^2}
\int_0^\infty
\langle\delta\mathbf B_s(0)\cdot\delta\mathbf B_s(z)\rangle\,dz
=
2\pi C(\nu)\ell_s.
\tag{28}
$$

Equation (28) applies to the spectrum in Equation (27), not to every turbulence model.

For r_*=r_L/ell_s and epsilon_s²=delta B_s²/B0², substitution into Equation (26) gives

$$
D_{\mu\mu}
=
\pi C(\nu)\epsilon_s^2
\frac{v}{\ell_s}r_\ast^{s-2}
(1-\mu^2)|\mu|^{s-1}
(1+r_\ast^2\mu^2)^{-s/2}.
\tag{29}
$$

The equivalent nonsingular expression away from mu=0 is obtained directly from Equation (26); Equation (29) makes the low-mu behavior explicit.

For 1<s<2, the exact mean free path for this unbroken high-wavenumber spectrum is

$$
\boxed{
\lambda_\parallel
=
\frac{3\ell_s}{4\pi C(s/2)\epsilon_s^2}
r_\ast^{2-s}
\int_0^1
(1-\mu^2)\mu^{1-s}
(1+r_\ast^2\mu^2)^{s/2}\,d\mu.
}
\tag{30}
$$

This equation is an algebraic consequence of Equations (20), (25), and (27). It defines the qlt_slab_spectrum analytic backend without relying on a schematic proportionality.

### 8.3 Inertial-range approximation

When r_*<<1 and the resonant interval remains within the inertial range,

$$
\boxed{
\lambda_\parallel^{\rm inertial}
=
\frac{3}
{2\pi C(s/2)(2-s)(4-s)}
\ell_s\frac{B_0^2}{\delta B_s^2}
\left(\frac{r_L}{\ell_s}\right)^{2-s},
\qquad 1<s<2.
}
\tag{31}
$$

For s=5/3,

$$
\lambda_\parallel^{\rm inertial}
\simeq
5.16465750497\,
\ell_s\,\frac{B_0^2}{\delta B_s^2}
\left(\frac{r_L}{\ell_s}\right)^{1/3}.
\tag{32}
$$

The displayed numerical prefactor is rounded; compute it from Equation (31) for numerical evaluation. It is specific to Equation (27) and the bend-over length ell_s. A formula written with another spectrum or another correlation-length convention can have a different prefactor without contradicting the r_\ast^(1/3) scaling.

At fixed rigidity and fixed slab variance, Equation (31) implies

$$
\lambda_\parallel
\propto
\ell_s^{s-1} B_0^s\,\mathcal R^{\,2-s}
(\delta B_s^2)^{-1}.
\tag{33}
$$

A statement that lambda always increases when B0 decreases is therefore incorrect without specifying how the fluctuation variance also changes. At fixed fractional variance epsilon_s², the B0 dependence becomes B0^(s-2).

For r_*>>1, Equation (30) has the leading behavior

$$
\lambda_\parallel
\sim
\frac{3\ell_s}{16\pi C(s/2)\epsilon_s^2}r_\ast^2.
\tag{34}
$$

The inertial approximation must not be used for this high-rigidity limit.

### 8.4 Multirange and measured spectra

A normalized three-range spectrum can be defined with x=k ell_s and x_d=k_d ell_s>1:

$$
\begin{aligned}
\Phi(x)&=
\begin{cases}
x^{q_E},&0<x\leq1,\\
x^{-s},&1<x\leq x_d,\\
x_d^{s_d-s}x^{-s_d},&x>x_d,
\end{cases}\\
P_s(k)&=\frac{\delta B_s^2\ell_s}{J}\Phi(k\ell_s),\\
J&=
\frac{1}{q_E+1}
+\frac{x_d^{1-s}-1}{1-s}
+\frac{x_d^{1-s}}{s_d-1}.
\end{aligned}
\tag{35}
$$

Here q_E>-1, s_d>1, and the displayed expression assumes s!=1. For s=1, the middle term is ln(x_d). This is a proposed explicit spectrum option; its parameters are not claimed to reproduce a particular publication unless set accordingly.

For supplied spectra, require positive k, nonnegative P, an identified normalization, and an explicit tail policy. A temporal spacecraft PSD is not automatically a slab wavenumber spectrum. Under a valid frozen-flow mapping with a known sampling velocity U_T along the relevant direction,

$$
k=\frac{2\pi f}{U_T},
\qquad
P_s(k)=\frac{U_T}{2\pi}
S_f\!\left(\frac{U_Tk}{2\pi}\right).
\tag{36}
$$

Here f is in cycles per second and S_f is a one-sided frequency PSD. An anisotropic geometry needs the appropriate sampling projection; a single measured PSD generally cannot identify slab and 2D power uniquely.

### 8.5 The 90-degree issue

For an unlimited inertial spectrum, D_mu_mu is proportional to |mu|^(s-1) near zero. Although D_mu_mu vanishes there, Equation (20) remains integrable when s<2. That mathematical convergence does not remove the physical limitations of resonant QLT near 90 degrees.

A dissipation spectrum with s_d>=2 makes the asymptotic inverse-scattering integral diverge. A hard upper-wavenumber cutoff leaves an interval with no resonant scattering and also gives no finite diffusive mean free path under strict magnetostatic QLT.

**Implementation requirement:** evaluate convergence from the actual high-k spectrum. Do not silently exclude |mu|<mu_min and report the resulting finite coefficient as the QLT value. Select a documented broadening/regularization model when a finite coefficient is required.

### 8.6 Numerical evaluation

For Equation (30), use

$$
\mu=z^{1/(2-s)},\qquad 0\leq z\leq1.
\tag{37}
$$

This cancels the integrable mu^(1-s) endpoint singularity. For s=5/3, mu=z³ and the transformed integrand is

$$
3(1-z^6)(1+r_\ast^2z^6)^{5/6}.
$$

For general spectra, partition the mu integral at the values corresponding to spectral breaks, evaluate the logarithm of positive spectra when appropriate, and distinguish missing spectrum coverage from a genuine zero in wave power.

## 9. Resonance broadening and weakly nonlinear theory

### 9.1 Physical idea and precise scope

Sharp QLT resonance assumes an unperturbed orbit and long-lived fluctuations. Finite wave lifetime, random sweeping, orbit perturbations, and displacement across a transversely structured field broaden that resonance. Weakly nonlinear theory (WNLT) treats orbit decorrelation and parallel/perpendicular transport together; the original composite-turbulence theory is Shalchi et al. (2004) [R06]. Shalchi (2009) treats transitions between quasilinear and isotropic pitch-angle scattering [R11].

There is no single universal “WNLT coefficient.” The following equations define explicit **slab resonance-broadened models** for implementation. They belong to that broader physical family, but they must not be labeled the complete paper-specific WNLT theory. A separate exact WNLT implementation requires its three-dimensional spectral tensor and orbit-correlation closure, including its perpendicular dynamics.

### 9.2 General resonance function

For a real decorrelation factor C_k(t) with C_k(0)=1, define

$$
\mathscr R_k(\omega)
=
{\rm Re}\int_0^\infty C_k(t)e^{i\omega t}\,dt.
\tag{38}
$$

Under the balanced, magnetostatic slab conventions of Equation (25),

$$
\boxed{
D_{\mu\mu}^{\rm broad}
=
\frac{\Omega^2(1-\mu^2)}{4B_0^2}
\int_0^\infty P_s(k)
\left[
\mathscr R_k(kv\mu-\Omega)
+\mathscr R_k(kv\mu+\Omega)
\right]\,dk.
}
\tag{39}
$$

Equation (39) and the kernels below are fully defined. Their amplitude normalization is fixed by their QLT limit. For C_k(t)=1, the resonance function is pi delta(omega), and Equation (39) reduces to Equation (26).

### 9.3 Exponential decorrelation: Lorentzian resonance

For C_k(t)=exp[-Gamma_dec(k)t],

$$
\mathscr R_k(\omega)
=
\frac{\Gamma_{\rm dec}(k)}{\omega^2+\Gamma_{\rm dec}(k)^2}.
\tag{40}
$$

Gamma_dec has units s^-1, and the integral of this kernel over all omega is pi. A completely specified minimal option is Gamma_dec(k)=Gamma_dec0>0. Another explicitly configured phenomenological option is

$$
\Gamma_{\rm dec}(k)=\Gamma_{{\rm dec},0}+u_{\rm dec}\,k,
\tag{41}
$$

where u_dec is a nonnegative decorrelation speed. These choices are implementation variants, not universal measurements of solar-wind decorrelation.

### 9.4 Gaussian decorrelation

For C_k(t)=exp[-Delta(k)^2 t²/4],

$$
\mathscr R_k(\omega)
=
\frac{\sqrt{\pi}}{\Delta(k)}
\exp\left[-\frac{\omega^2}{\Delta(k)^2}\right],
\qquad \Delta(k)>0.
\tag{42}
$$

Its integrated area is also pi. The minimal option is a supplied positive constant Delta0. The definition of Delta must be preserved; it is not interchangeable with the standard deviation of the Gaussian frequency profile.

Finite widths can provide nonzero D_mu_mu at mu=0. The resulting mean free path is still computed with Equation (20), using both kernels in Equation (39).

### 9.5 More complete nonlinear decorrelation

For a three-dimensional mode, an orbit-decorrelation approximation may introduce

$$
\Gamma_{\rm orb}(\mathbf k)
=
\kappa_\perp k_\perp^2+\kappa_\parallel k_\parallel^2.
\tag{43}
$$

This illustrates why parallel and perpendicular coefficients can become coupled. It does not turn the one-dimensional slab spectrum in Equation (39) into a complete 3D theory. Gyroharmonics, tensor polarization, and the assumed correlation functions must also be specified.

**Implementation requirement:** the v1 broadened_slab model accepts only its declared slab spectrum and kernel. It must reject a request to add an unprovided 2D/3D contribution. For small widths, integrate around resonance centers adaptively or switch to the analytically defined zero-width QLT branch; evaluating a very narrow peak on a fixed grid is inadequate.

### 9.6 Directional waves

For propagating Alfvén waves, the mismatch contains k_parallel v mu-omega_wave+n Omega, rather than k_parallel v mu+n Omega alone. For low-energy ions with v comparable to the Alfvén speed, that change and the scattering-frame physics matter.

Equations (26) and (39) must not be advertised as arbitrary imbalanced Alfvén-wave scattering formulas. A polarized, directional implementation needs separate spectral components and the appropriate electric/magnetic response factors.

## 10. NLPA and coupled nonlinear models

### 10.1 Published closure and corrected parameter

Qin (2007) derives an implicit nonlinear parallel approximation (NLPA), with perpendicular diffusion as input [R08]. Its corrected parameter is given by Qin (2013) and restated by Qin and Zhang (2014) [R09,R10]. The latter paper combines NLPA with a modified perpendicular closure.

Use kappa_z=kappa_parallel and kappa_x=kappa_perp. With a spectral tensor normalized so that the integral of S_xx over d³k is the x-component magnetic variance, define

$$
A(\mathbf k)
=
\frac{v}{\lambda_\parallel}
+k_\perp^2\kappa_x+k_\parallel^2\kappa_z.
$$

The published NLPA relation is

$$
\boxed{
\kappa_z^{-1}
=
6a_x\left(\frac{\Omega}{v}\right)^2
\int d^3k\,\frac{S_{xx}(\mathbf k)}{B_0^2}
\frac{A(\mathbf k)}
{\Omega^2+[A(\mathbf k)+\Gamma_{\rm dec}(\mathbf k)]^2}.
}
\tag{44}
$$

The numerator is A, not A+Gamma_dec. This matters when reproducing a dynamical-turbulence version. The magnetostatic specification below sets Gamma_dec=0.

Define

$$
\epsilon=\frac{\sqrt{\delta B_s^2+\delta B_2^2}}{B_0},
\quad
f_s=\frac{\delta B_s^2}{\delta B_s^2+\delta B_2^2},
\quad
\widetilde r=\frac{2\pi r_L}{L_c},
\quad
\xi=\frac{\widetilde r}{\epsilon}.
$$

The corrected coefficient is

$$
\boxed{
a_x=
\frac12
\left[
\frac{f_s}
{\frac{\xi}{1+\xi}\frac{1}{\epsilon}
+\frac{\epsilon}{2\xi}}
\right]^{1/2}.
}
\tag{45}
$$

**Implementation requirement:** use Equation (45), not the uncorrected parameter in the original 2007 printing. NLPA is not required to recover QLT in a supposed weak-turbulence limit; that is not a valid acceptance test for this particular closure [R10].

**Correlation-length requirement:** in Equation (45), L_c means the slab integral correlation length L_c,s=2 pi C(nu) ell_s of Equation (28). This is the definition stated by Qin and Zhang (2014), immediately after their Equation (6), even when the turbulence contains a 2D component [R10]. Total variance defines epsilon, whereas the slab spectrum defines L_c. Thus

$$
\xi=\frac{r_L/\ell_s}{C(\nu)\epsilon}.
\tag{45a}
$$

L_c is not an independently adjustable composite length for this backend. A variance-weighted combination of slab and 2D lengths, or substitution of ell_2, defines a different model and must not retain the published NLPA/NLGCE identifier.

Equation (45) uses the multiplicative factor (xi/(1+xi)) times (1/epsilon), followed by epsilon/(2xi), as printed in Qin and Zhang (2014), Equation (6), restating the 2013 correction [R09,R10]. The first factor is not raised to the power 1/epsilon. Both occurrences use epsilon, not epsilon². Implement the denominator as (xi/(1+xi))/epsilon + epsilon/(2xi). Revision 1.2 corrects the exponentiation error in earlier versions of this specification and regenerates all affected numerical results; it does not refit the published parameter.

### 10.2 Fully specified two-component reduction

For each component a=s,2, use

$$
P_a(k)
=
4C(\nu)\ell_a\delta B_a^2
[1+(k\ell_a)^2]^{-\nu},
\qquad
\int_0^\infty P_a(k)\,dk=\delta B_a^2.
\tag{46}
$$

For the slab component k means |k_parallel|. For the 2D component it means k_perp, after the angular integration of the axisymmetric transverse tensor. Thus P_2 is a reduced radial spectral density, not a Cartesian one-dimensional slice.

Set nu=5/6 for the published 2014 fit and magnetostatic coupled model. Define

$$
A_s(k)=\frac{v^2}{3\kappa_z}+k^2\kappa_z,
\qquad
A_2(k)=\frac{v^2}{3\kappa_z}+k^2\kappa_x.
\tag{47}
$$

Equation (44) then becomes

$$
\boxed{
\kappa_z^{-1}
=
\frac{3a_x\Omega^2}{v^2B_0^2}
\left[
\int_0^\infty
P_s(k)\frac{A_s(k)}{\Omega^2+A_s(k)^2}\,dk
+
\int_0^\infty
P_2(k)\frac{A_2(k)}{\Omega^2+A_2(k)^2}\,dk
\right].
}
\tag{48}
$$

The factors in Equation (48) follow from the one-sided, total-variance convention: each transverse component has half the total variance. This reduction provides an implementation check against the tensor equation.

For nlpa_given_perp, supply a positive kappa_x from an independent prescription and solve Equation (48) for kappa_z. Report which perpendicular model was supplied.

### 10.3 Original NLGC–NLPA pair

For the original coupled NLGC-E option,

$$
\kappa_x=
\frac{a^2v^2}{6B_0^2}
\left[
\int_0^\infty\frac{P_s(k)}{A_s(k)}\,dk+
\int_0^\infty\frac{P_2(k)}{A_2(k)}\,dk
\right],
\qquad a^2=\frac13,
\tag{49}
$$

together with Equation (48). This is the magnetostatic axisymmetric two-component form of the NLGC relation combined with NLPA [R05,R08].

### 10.4 Modified coupled pair: NLGCE-N

For the modified perpendicular closure, define

$$
\boxed{
a'^2=
\left[
\sqrt{\frac{\ell_2}{\ell_s}}\frac{1}{f_s}
+\frac{4}{3}\frac{1}{1-f_s}
\right]^{-1}.
}
\tag{50}
$$

The NLGCE-N pair is Equation (48) and

$$
\boxed{
\kappa_x=
\frac{a'^2v^2}{6B_0^2}
\int_0^\infty
\frac{P_2(k)}{A_2(k)}\,dk.
}
\tag{51}
$$

Only 2D power appears directly in Equation (51); slab power still affects kappa_x through the coupled parallel coefficient and a'^2. This is the modified model of Qin and Zhang (2014) [R10]. Its polynomial approximation is specified in Section 11.

Use lambda_parallel=3 kappa_z/v and lambda_perp=3 kappa_x/v. A textual occurrence of kappa_zz in the source's definition of the perpendicular mean free path is a typographical error; the perpendicular component is kappa_xx.

### 10.5 Numerical solution

**Implementation requirements:**

1. Evaluate all input ratios and the corrected a_x before solving.
2. Solve in logarithmic unknowns y_z=ln(kappa_z/(v ell_s)) and y_x=ln(kappa_x/(v ell_s)) so that both coefficients remain positive.
3. Use dimensionless integration variables x=k ell_a. Equations (48)–(51) then contain finite, normalized component spectra.
4. Solve the two residuals simultaneously. A pair of separate evaluations using stale values of the other coefficient is not a converged coupled solution.
5. A proposed numerical acceptance criterion is a maximum absolute logarithmic residual below 10^-8, with tighter quadrature errors than the requested residual.
6. Use continuation from a nearby validated state, or an in-domain NLGCE-F value, as an initial guess. Check more than one positive initial guess at representative points.
7. Record iteration count, residual, quadrature status, and any branch dependence.
8. Do not return the last iterate as a successful result when the solver fails.

The intended ordinary domain has B0>0, v>0, both bend-over lengths positive, nonzero total magnetic variance, and 0<f_s<1. Pure-component limits require an explicitly developed limit evaluator; the displayed mixed-component formulas must not be evaluated by dividing by zero or inserting an arbitrary small fraction.

For general nu, L_c must be recomputed from Equation (28). The empirical parameters and accuracy assessment of the 2014 fit apply to its stated spectrum; changing the energy-range shape or spectral index creates a different closure.

An NLGCE-F starting value is only an initial iterate. Its proximity to the nonlinear solution must be checked for the actual four-dimensional state. The acceptance decision comes from the integral residuals and quadrature convergence. A successful solve reproduces this published closure; it does not independently establish its physical validity or agreement with the published polynomial.

## 11. NLGCE-F polynomial model

### 11.1 Exact input variables and polynomial

The published NLGCE-F model fits the magnetostatic NLGCE-N solution in four logarithmic variables [R10]:

$$
x_1=\ln\frac{r_L}{\ell_s},
\quad
x_2=\ln f_s,
\quad
x_3=\ln\frac{\delta B_s^2+\delta B_2^2}{B_0^2},
\quad
x_4=\ln\frac{\ell_s}{\ell_2}.
\tag{52}
$$

For alpha=parallel or perpendicular,

$$
\boxed{
F_\alpha=
\ln\frac{\lambda_\alpha}{\ell_s}
=
\sum_{i=0}^{5}\sum_{j=0}^{3}\sum_{k=0}^{3}\sum_{l=0}^{2}
d^\alpha_{ijkl}\,x_1^i x_2^j x_3^k x_4^l,
\qquad
\lambda_\alpha=\ell_s e^{F_\alpha},
\qquad
\kappa_\alpha=\frac{v\lambda_\alpha}{3}.
}
\tag{53}
$$

All logarithms are natural. Each coefficient array contains 6×4×4×3=288 values. The complete parallel and perpendicular arrays are included in Section 20, preserving the published decimal precision.

The variables must not be replaced with log(r_L/L_c), log(epsilon), or log(ell_2/ell_s). Such substitutions are different polynomials. The turbulence input is particularly easy to confuse: x3 uses epsilon², whereas Equation (45) uses epsilon.

### 11.2 Published fit domain

| Input ratio | Lower bound | Upper bound |
|---|---:|---:|
| r_L / ell_s | 10⁻⁵ | 6.3 |
| f_s | 10⁻³ | 0.85 |
| Total delta B² / B0² | 10⁻⁴ | 10² |
| ell_s / ell_2 | 1 | 10³ |

These are the ranges reported in Qin and Zhang (2014), Table 2 [R10]. The fit is not valid for pure 2D turbulence. Its finite tabulated domain also excludes a pure slab endpoint.

This table specifies the published input box, not a uniformly validated approximation-error domain. Section 11.4 demonstrates large discrepancies inside the box. Its numerical bounds are retained as published; neither those bounds nor a turbulence-bin median supply an error bound for a particular state.

**Implementation requirement:** check all four original ratios before taking logarithms. The default out-of-domain behavior is an explicit status, not polynomial extrapolation. A user-selected fallback to a nonlinear solver is a different evaluated backend and must be reported as such.

### 11.3 Polynomial evaluation and derivatives

Use nested Horner evaluation in the order l, k, j, i; retain float64 precision or better. The table layout in Section 20 lists j,k,l as row indices and the six i values as columns. Parse that layout explicitly instead of assuming the storage order of a C or Fortran array.

The logarithmic derivatives are obtained by differentiating the polynomial. For example,

$$
\frac{\partial F_\alpha}{\partial x_1}
=
\sum_{i=1}^{5}\sum_{j,k,l}
i\,d^\alpha_{ijkl}x_1^{i-1}x_2^j x_3^k x_4^l.
\tag{54}
$$

A useful spatial derivative identity is

$$
\nabla\ln\lambda_\alpha
=
\nabla\ln\ell_s
+\sum_{a=1}^4
\frac{\partial F_\alpha}{\partial x_a}\nabla x_a.
\tag{55}
$$

The x_a gradients include the local B0, turbulence, and length-scale gradients. At fixed turbulence and field,

$$
\frac{\partial\ln\kappa_\parallel}{\partial\ln\mathcal R}
=
\frac{1}{\gamma^2}
+\frac{\partial F_\parallel}{\partial x_1}.
\tag{56}
$$

**Implementation requirement:** return both fitted eigenvalues together, even when the caller requests only the parallel coefficient. Discarding the fitted perpendicular output and substituting a different coefficient does not change the number returned by the polynomial, but it changes the physical closure used by the transport calculation.

### 11.4 Fit error and provenance

The polynomial is a published surrogate for a nonlinear integral model, not an exact identity. Distinguish numerical reproduction of each backend, discrepancy between the backends, closure agreement with particle simulations/observations, and validity of the Parker diffusion limit.

For alpha=parallel or perpendicular, define the measured discrepancy

$$
e_\alpha=
\left|\frac{\lambda_\alpha^{\rm F}}{\lambda_\alpha^{\rm N}}-1\right|,
\tag{53a}
$$

Here F is Equation (53) with the printed coefficients and N solves Equations (48),(50),(51), using the corrected multiplicative a_x in Equation (45). The quantity measures disagreement with the closure; it is not observational error.

**Correction to earlier versions:** the former calculation exponentiated xi/(1+xi) by 1/epsilon instead of multiplying by 1/epsilon. Its weak-turbulence discrepancy claims and numerical tables were therefore calculations of a different, mistranscribed closure. They are replaced by the regenerated results below. This transcription error does not demonstrate an inconsistency in the published Qin models.

The corrected numerical audit uses the same 300 states drawn uniformly in the four natural logarithms of Equation (52), within the published box. The generator is NumPy PCG64 through default_rng(20261008), with a 300 by 4 draw in the order r_L/ell_s, f_s, epsilon², ell_s/ell_2. The complete input/output CSV is authoritative for reproducing the sample.

| Total epsilon² bin | States | Median parallel discrepancy | Parallel share >25% | Median perpendicular discrepancy | Perpendicular share >25% |
|---|---:|---:|---:|---:|---:|
| 10^-4 to <10^-3 | 48 | 18.2% | 29.2% | 10.0% | 31.2% |
| 10^-3 to <10^-2 | 59 | 16.2% | 32.2% | 11.1% | 28.8% |
| 10^-2 to <0.1 | 49 | 10.2% | 12.2% | 8.6% | 10.2% |
| 0.1 to <1 | 42 | 8.6% | 16.7% | 7.7% | 11.9% |
| 1 to <10 | 48 | 13.2% | 14.6% | 6.5% | 4.2% |
| 10 to 100 | 54 | 10.7% | 22.2% | 5.0% | 3.7% |

All 300 integral solves met a maximum absolute logarithmic residual below 10^-8. Component log-wavenumber quadrature used numerical intervals [-40,40]. At the three largest parallel discrepancies, extending the intervals to [-55,55], tightening quadrature, and reducing initial mean free paths by a factor of 10^4 changed both coefficients by less than 10^-12 relative. The measured changes are recorded in the companion summary. These finite intervals are refinement controls, not physical spectrum cutoffs.

The table is a statistic of a defined log-uniform sample, not a probability distribution of heliospheric states or a local error bound. The largest discrepancies occur at particular four-dimensional inputs; medians conceal those outliers. The complete dataset permits analysis of their locations. A low median does not certify all states in a turbulence bin.

For point F, r_L/ell_s=0.01, f_s=0.2, epsilon²=10^-4, and ell_s/ell_2=10, the corrected nonlinear parallel mean free path is 87716.3439308 ell_s and the polynomial gives 95850.8828433 ell_s. Their signed relative difference is +9.274%, rather than the enormous discrepancy obtained from the former exponentiation error. The other explicit checks are in Section 15.4.

The attached review correctly motivated checks of source conventions, coefficients and more input states. Its large weak-turbulence numbers cannot be adopted as results for the correctly transcribed closure. The slab L_c definition and all 576 printed coefficients agree with the source [R09,R10]. Correcting a transcription is distinct from adjusting a_x to force agreement with the polynomial.

**Implementation requirements:** keep nlgce_n and nlgce_f_2014 as distinct backend identifiers; preserve their published expressions and coefficient tables; record the evaluated backend and coefficient-set digest; expose the general surrogate diagnostics in Section 2.3. Check the multiplicative a_x with a weak-turbulence regression point, since epsilon=1 makes multiplication and exponentiation indistinguishable. Any caller-selected fallback must record both requested and evaluated models.

Qin and Shen (2017) and Shen and Qin (2018) apply NLGCE-F to GCR modulation [R16,R24]. Those applications provide provenance; they do not establish a uniform error bound over the entire four-dimensional box.

## 12. Turbulence and wave-spectrum adapters

### 12.1 Separation of provider and scattering closure

A turbulence-transport calculation can supply local field magnitude, density, fluctuation variance, slab fraction, spectral lengths, and spectral slopes. A closure then maps these quantities to D_mu_mu or kappa_parallel.

The sequence is

**background/turbulence state → explicitly normalized spectrum or moments → selected scattering closure → mean free path and spatial coefficient.**

This separation allows the same library to operate with analytic Parker backgrounds, observational drivers, a reduced solar-wind model, or a global simulation. An MHD solver is not required by the coefficient evaluator.

Engelbrecht and Burger (2013a,b) use turbulence transport to obtain GCR diffusion coefficients [R12,R13]. Chhiber et al. (2017) evaluate spatially varying heliospheric coefficients using global background and turbulence data [R17]. Zhao et al. (2018) examine how solar-cycle changes in turbulence and the mean field affect the inferred coefficients [R23]. Wijsen et al. (2023) provide a recent SEP application of turbulence-dependent mean free paths [R19].

### 12.2 Converting turbulence moments

Suppose the supplied model defines an Elsässer fluctuation energy

$$
Z^2=
\langle\delta u^2\rangle+
\langle\delta b_A^2\rangle,
\qquad
\delta\mathbf b_A=
\frac{\delta\mathbf B}{\sqrt{\mu_0\rho}},
$$

and normalized residual energy

$$
\sigma_D=
\frac{\langle\delta u^2\rangle-\langle\delta b_A^2\rangle}
     {\langle\delta u^2\rangle+\langle\delta b_A^2\rangle}.
$$

Then

$$
\boxed{
\delta B^2=
\mu_0\rho\,\frac{1-\sigma_D}{2}Z^2,
\qquad
\delta B_s^2=f_s\delta B^2,
\qquad
\delta B_2^2=(1-f_s)\delta B^2.
}
\tag{57}
$$

Equation (57) is an algebraic conversion for the displayed definitions. Some providers define Z² differently; their conversion must be derived rather than guessed.

Let E_±=⟨delta z^±·delta z^±⟩ for delta z^±=delta u±delta b_A. Then

$$
E_++E_-=2\left(\langle\delta u^2\rangle+
\langle\delta b_A^2\rangle\right)=2Z^2.
\tag{57a}
$$

The specification therefore uses the half-sum Elsässer convention. An unhalved Elsässer sum differs by a factor of two; neither convention can be selected from the field name alone.

| Provider quantity M_prov, in m²/s² | Named convention | chi_Z in M_prov=chi_Z Z² | Conversion to canonical Z² |
|---|---|---:|---|
| ⟨delta u²⟩+⟨delta b_A²⟩ | kinetic_plus_magnetic_variance | 1 | M_prov |
| (E_++E_-)/2 | half_elsasser_sum | 1 | M_prov |
| E_++E_- | elsasser_sum | 2 | M_prov/2 |
| (⟨delta u²⟩+⟨delta b_A²⟩)/2 | specific_total_fluctuation_energy | 1/2 | 2 M_prov |

For a declared provider multiplier, the direct conversion is

$$
\delta B^2=\mu_0\rho\,\frac{1-\sigma_D}{2\chi_Z}M_{\rm prov}.
\tag{57b}
$$

The residual-energy sign also needs a convention: kinetic-minus-magnetic residual energy is sigma_D as defined here; magnetic-minus-kinetic residual energy must be negated before applying Equation (57). These conversions presume compatible averaging and the same density normalization for delta b_A. They do not reconstruct directional wave spectra or cross helicity from a total energy alone.

**Implementation requirement:** a moment adapter accepts a named energy convention and residual-energy convention, converts once to the canonical variables, and retains the original convention in provenance. Reject an unknown convention rather than assigning a factor of two by guesswork.

Require rho>0, Z²>=0, -1<=sigma_D<=1, and 0<=f_s<=1. The provider must identify whether rho is total mass density or a particular species density.

For a magnetic fluctuation energy density E_B and a total equipartition Alfvén-wave energy density E_wave,

$$
\delta B^2=2\mu_0 E_B,
\qquad
\delta B^2=\mu_0 E_{\rm wave}
\quad\text{under the stated magnetic/kinetic equipartition convention}.
\tag{58}
$$

The factor of two must not be inferred from a field simply named “wave energy.”

### 12.3 Radial scaling as a derived result

For illustration, let

$$
B_0\propto r^{-m_B},\qquad
\delta B_s^2\propto r^{-m_s},\qquad
\ell_s\propto r^{m_\ell}.
$$

Equation (31), at fixed rigidity, then gives

$$
\lambda_\parallel\propto
r^{\,m_\ell(s-1)-m_Bs+m_s}.
\tag{59}
$$

This explains why a single assumed radial power law need not remain valid from the low corona to the outer heliosphere. The mean field, spectral scale, and turbulent variance can change regime independently.

### 12.4 Supplied and self-generated wave spectra

If a wave solver supplies magnetic spectral power, the adapter evaluates the selected QLT or broadened scattering closure on the current spectrum. If it supplies wave energy rather than magnetic power, first apply its declared magnetic-energy convention.

A schematic wave-energy transport equation is

$$
\frac{\partial W_\pm}{\partial t}
+\nabla\cdot[(\mathbf U\pm v_A\mathbf b)W_\pm]
=
2\gamma_\pm W_\pm
-2\Gamma_{{\rm d},\pm}W_\pm
+\mathcal C_\pm[W]+Q_\pm.
\tag{60}
$$

Here the displayed growth and damping rates are amplitude rates, hence the factors of two in the energy equation. This equation illustrates the provider's responsibilities; it does not specify a unique cascade or particle-growth calculation.

Vainio and Laitinen (2007), Afanasiev and Vainio (2013), and Afanasiev et al. (2015) treat particle transport and wave growth self-consistently [R20,R21,R22]. Their wave dynamics must not be replaced by a static function of radius while retaining the label “self-consistent.”

**Implementation requirements:**

- The coefficient library reads wave state; it does not update particle streaming or grow waves.
- The state identifies propagation direction, polarization, units, normalization, frame, and update time.
- The balanced magnetostatic limit can use the combined slab spectrum of Section 8 when justified.
- An imbalanced/propagating extension requires a compatible scattering-center convection and energy-change model.
- Additional momentum diffusion or wave-frame energy changes are not determined by kappa_parallel alone.

## 13. Bohm and supplied-table models

### 13.1 Bohm comparison

Define the mean-field Bohm prescription by

$$
\lambda_\parallel=\eta_B r_L,
\qquad
\kappa_\parallel=\eta_B\frac{vr_L}{3},
\qquad \eta_B>0.
\tag{61}
$$

eta_B=1 is the conventional Bohm comparison value. Some applications use an effective total field instead of B0; that is a separate option whose field definition must be reported.

This is a useful strong-scattering comparison, especially near shocks, but it is not a universal solar-wind coefficient or a rigorous universal lower bound. Hussein and Shalchi (2014) examine the conditions associated with Bohm scaling [R29].

**Implementation requirement:** never clamp other models to kappa_B automatically. If a caller explicitly requests a bound, preserve the unclamped coefficient and report that a numerical/phenomenological bound was applied.

### 13.2 Supplied coefficient tables

A table model must identify whether the stored quantity is lambda_parallel or kappa_parallel and provide the independent variables and units. A rigidity table is not an energy-per-nucleon table.

Use linear interpolation of the logarithm of a positive coefficient versus the logarithms of positive axes. For time, use the explicitly supplied linear/step rule rather than taking its logarithm. Multilinear interpolation is the minimal reproducible multidimensional option; higher-order interpolation requires a named method and positivity checks.

Require an explicit policy at every boundary. The default is “out of domain.” Avoid automatic extrapolation. Values from a table generated with one turbulence spectrum must retain that spectrum's provenance even when the runtime uses a different background provider.

## 14. Library interface and numerical rules

### 14.1 Proposed data structures

The following is a C++20 interface sketch, not a completed library or a mandate to use a particular language. Optional quantities are represented explicitly, including outputs unavailable after failure. Spectrum handles stand for separately specified provider objects.

~~~cpp
#include <array>
#include <cstdint>
#include <optional>
#include <span>
#include <string>

using Vec3 = std::array<double, 3>;
using Mat3 = std::array<Vec3, 3>;

enum class Status {
  success,
  invalid_particle,
  invalid_background,
  missing_input,
  outside_model_domain,
  infinite_mean_free_path,
  integration_failed,
  nonlinear_solver_failed,
  inconsistent_spectrum
};

enum class Diagnostic : std::uint32_t {
  diffusion_limit_concern = 1u << 0,
  strong_turbulence = 1u << 1,
  qlt_weak_perturbation_concern = 1u << 2,
  surrogate_error_unbounded = 1u << 3,
  fallback_used = 1u << 5,
  user_bound_applied = 1u << 6,
  derivative_unavailable = 1u << 7
};

struct ParticleState {
  double mass_kg;
  double charge_C;
  double momentum_kg_m_s;
};

struct TurbulenceState {
  std::optional<double> density_kg_m3;
  std::optional<double> slab_variance_T2;
  std::optional<double> two_d_variance_T2;
  std::optional<double> slab_bendover_length_m;
  std::optional<double> two_d_bendover_length_m;
  std::optional<double> inertial_index;
  // Canonical moments after the declared provider conversion.
  std::optional<double> canonical_Z2_m2_s2;
  std::optional<double> canonical_sigma_D;
  std::optional<Vec3> grad_density_kg_m4;
  std::optional<Vec3> grad_slab_variance_T2_per_m;
  std::optional<Vec3> grad_two_d_variance_T2_per_m;
  std::optional<Vec3> grad_slab_length_m_per_m;
  std::optional<Vec3> grad_two_d_length_m_per_m;
  std::optional<Vec3> grad_inertial_index_per_m;
  std::optional<Vec3> grad_Z2_m_s2;
  std::optional<Vec3> grad_sigma_D_per_m;
  std::optional<std::uint64_t> slab_spectrum_handle;
  std::optional<std::uint64_t> directional_wave_spectrum_handle;
};

struct LocalState {
  double time_s;
  Vec3 position_m;
  std::optional<Vec3> mean_B_T;
  // Matrix convention: dB_dx[i][j] = partial B_i / partial x_j.
  std::optional<Mat3> dB_dx_T_per_m;
  std::optional<TurbulenceState> turbulence;
  std::uint64_t background_revision;
  std::uint64_t turbulence_revision;
};

struct PerpendicularPair {
  double kappa_perp_m2_s;
  double lambda_perp_m;
};

struct Provenance {
  std::string requested_model_id;
  std::string evaluated_model_id;
  std::string specification_version;
  std::string configuration_sha256;
  std::optional<std::string> coefficient_set_sha256;
  std::string input_moment_convention;
  std::uint64_t background_revision;
  std::uint64_t turbulence_revision;
};

struct ParallelResult {
  Status status;
  std::optional<double> kappa_parallel_m2_s;
  std::optional<double> lambda_parallel_m;
  std::optional<PerpendicularPair> perpendicular;
  std::optional<Vec3> grad_kappa_parallel_m_s;
  std::optional<Vec3> grad_kappa_perp_m_s;
  std::uint32_t diagnostic_mask;
  Provenance provenance;
  // An optional D_mu_mu evaluator is a separate typed interface.
};

struct ModelConfiguration;  // Defined by the selected model's schema.

ParallelResult evaluate_parallel(
    const ParticleState&,
    const LocalState&,
    const ModelConfiguration&);

void evaluate_parallel_batch(
    std::span<const ParticleState> particles,
    std::span<const LocalState> states,
    const ModelConfiguration& configuration,
    std::span<ParallelResult> results);
~~~

Particle-derived quantities are computed centrally from Equation (7) or the equivalent momentum formulation. A model must not recompute v or rigidity using a different relativistic approximation.

For NLGCE-F and nonlinear pairs, include the internally consistent perpendicular value. For empirical and QLT-only models, its absence is explicit.

**Batch requirements:** all three spans have the same length N; configuration is shared by the batch; result i corresponds to particle/state i. Empty batches are valid. A shape mismatch is a configuration error detected before evaluation. Each point retains its own status, diagnostics, provenance, and optional outputs; one invalid state must not invalidate unrelated points. There is no implicit broadcasting. A separate convenience adapter may explicitly broadcast a shared state.

The scalar and batch paths must evaluate the same model to the same declared tolerance. Batching permits vectorized polynomial evaluation and grouped quadrature, but does not relax accuracy checks or silently change the backend. Parallel evaluation must use immutable provider snapshots or independently synchronized caches, with the revision rules of Section 14.4.

Missing gradients remain absent. A genuinely uniform provider can supply explicit zero derivatives. A successful scalar coefficient evaluation is not a statement that a caller has enough derivatives to assemble a transport operator.

### 14.2 Evaluation sequence

1. Validate the particle state and compute p, gamma, v, and rigidity.
2. Validate only the background fields required by the selected model; compute Omega and r_L when the selected model needs the magnetic field. A constant mean-free-path or constant-coefficient model does not require a magnetic-field input.
3. Convert the provider's fluctuation conventions into the declared spectral/moment conventions.
4. Check the analytical and numerical model domain.
5. Evaluate the coefficient and any explicitly requested derivatives.
6. Verify lambda=3 kappa/v and positivity.
7. Attach numerical status, validity diagnostics, and provenance.
8. Let the transport solver assemble/project the tensor.

A requested neutral particle or p<=0 is outside the ordinary charged-particle coefficient API. A zero turbulence variance represents absence of scattering in the corresponding model; it must not be silently replaced by a small positive variance. At B0=0 the field-aligned decomposition is undefined.

### 14.3 Derivatives needed by transport solvers

For smooth diffusion-only transport in local Cartesian coordinates, a stochastic representation involves the diffusion drift

$$
dx_i=
\frac{\partial K_{ij}}{\partial x_j}\,dt
+L_{ia}\,dW_a,
\qquad
\mathbf L\mathbf L^T=2\mathbf K_s.
\tag{62}
$$

This statement is limited to the diffusion operator; the full forward/backward Parker SDE depends on the evolved distribution variable, convection, coordinates, and momentum terms.

For the axisymmetric tensor in Equation (3), define Delta kappa=kappa_parallel-kappa_perp. Its Cartesian divergence is

$$
\frac{\partial K_{ij}}{\partial x_j}
=\partial_i\kappa_\perp
+b_i b_j\partial_j(\Delta\kappa)
+\Delta\kappa\left[
b_j\partial_j b_i+b_i\partial_j b_j
\right].
\tag{62a}
$$

Thus the full drift requires both coefficient gradients and derivatives of the field direction. The latter can be derived from a supplied mean-field Jacobian:

$$
\partial_j b_i=
\frac{\delta_{ik}-b_i b_k}{B_0}\,\partial_j B_{0k}.
\tag{62b}
$$

Repeated Cartesian indices are summed, and dB_dx[i][j] means partial_j B_0i. Equations (62a),(62b) are algebraic consequences of Equation (3). The background provider supplies the field Jacobian; the coefficient evaluator supplies gradients of the coefficients it determines; the transport solver assembles the tensor and its divergence. If perpendicular diffusion comes from a separate model, that model must supply its value and gradient. For Equation (2) with unequal perpendicular coefficients, derivatives of both perpendicular basis vectors are also needed.

**Implementation requirement:** an operation requesting a complete tensor divergence must report missing required derivatives. It must not silently replace absent field-direction or perpendicular-coefficient derivatives with zero.

For a flux tube of cross-sectional area A_tube(s_arc), the parallel operator on a scalar f is

$$
\nabla\cdot(\kappa_\parallel\mathbf b\mathbf b\cdot\nabla f)
=
\frac{1}{A_{\rm tube}}\frac{\partial}{\partial s_{\rm arc}}
\left(A_{\rm tube}\kappa_\parallel
\frac{\partial f}{\partial s_{\rm arc}}\right).
\tag{63}
$$

Its first-derivative coefficient is d kappa_parallel/ds_arc+kappa_parallel d ln A_tube/ds_arc. Therefore a formula containing only d kappa_parallel/ds_arc is insufficient for a varying-area tube. This geometry belongs to the transport operator, not to an artificial “focusing correction” hidden in the coefficient.

**Implementation requirement:** analytic derivatives are preferred for explicit models. If derivatives are computed numerically, evaluate the whole declared state-dependent model and use provider-consistent gradients. Report derivative discontinuities at region boundaries, input-table knots, or bounds. Do not finite-difference across different background revisions.

### 14.4 Numerical tolerances and cache rules

Recommended initial numerical settings are float64 arithmetic, a relative integral tolerance of 10^-8 for general pitch-angle integrals, and tighter tolerances for reference calculations. A posteriori quadrature estimates and convergence under refinement determine acceptance.

Cache keys include species mass and charge, particle momentum, model/configuration version, and every relevant field/turbulence revision. A change in B0 alone changes r_L and the spectral resonance; reusing an old coefficient after a background update is incorrect.

A supplied gradient is optional. Its absence must not be interpreted as a zero gradient unless the provider explicitly declares uniform state.

### 14.5 Failure and fallback policies

| Situation | Required default behavior |
|---|---|
| Missing slab fraction or spectral length | missing_input; do not infer an undocumented value |
| B0<=0 for a field-aligned model | invalid_background |
| Zero resonant power causing divergent Equation (20) | infinite_mean_free_path |
| NLGCE-F ratio outside its fit domain | outside_model_domain |
| Failed coupled nonlinear solve | nonlinear_solver_failed |
| Incorrect spectral normalization | inconsistent_spectrum |
| Strong turbulence with otherwise admissible model inputs | Preserve the coefficient and return model-dependent validity diagnostics |
| In-box NLGCE-F with no certified local error bound | Return the exact polynomial and surrogate_error_unbounded diagnostic |
| Missing requested derivative | Preserve absence explicitly; complete derivative-dependent operation cannot report success |
| User-enabled alternate backend | Return actual backend identity and fallback flag |
| User-enabled lower/upper bound | Preserve raw result and bound diagnostics |

The solver, rather than this library, decides how to handle a ballistic/infinite-mean-free-path regime. A finite artificial replacement can affect SEP arrival times and shock acceleration; it is a modeling decision that must be recorded.

### 14.6 Configuration templates

The first three blocks are parameterized templates, not executable configurations. Values in angle brackets must be supplied from the chosen physical setup or a cited calibration; arbitrary numerical defaults are deliberately omitted. The NLGCE-F block selects the published backend and explicit evaluation policies.

The companion configuration_examples.json preserves these templates, with the same non-normative designation. Any later observational preset must identify its particle population, energy/rigidity interval, region, field geometry, and source. Palmer-consensus discussions [R03] provide context for selected heliospheric observations; their mean-free-path ranges are not universal bounds or automatic unit-validation cutoffs.

~~~json
{
  "model": "power_law_lambda",
  "lambda0": {
    "value": "<lambda0_AU>",
    "unit": "AU"
  },
  "rigidity0": {
    "value": "<rigidity0_GV>",
    "unit": "GV"
  },
  "radius0": {
    "value": "<radius0_AU>",
    "unit": "AU"
  },
  "rigidity_exponent": "<delta>",
  "radial_exponent": "<alpha>",
  "field_exponent": "<eta>",
  "normalization_factors": {
    "time": "unity",
    "region": "unity"
  }
}
~~~

~~~json
{
  "model": "broken_rigidity_kappa",
  "K_star": {
    "value": "<K_star_m2_per_s>",
    "unit": "m2/s"
  },
  "B_ref": {
    "value": "<B_ref_nT>",
    "unit": "nT"
  },
  "rigidity0": {
    "value": "<rigidity0_GV>",
    "unit": "GV"
  },
  "break_rigidity": {
    "value": "<break_rigidity_GV>",
    "unit": "GV"
  },
  "low_slope": "<a>",
  "high_slope": "<b>",
  "smoothness": "<h>",
  "field_exponent": "<eta>",
  "normalization_convention": "speed_factored_K_star"
}
~~~

~~~json
{
  "model": "qlt_slab_spectrum",
  "spectrum": {
    "form": "smooth_bendover",
    "slab_variance": {
      "value": "<slab_variance_nT2>",
      "unit": "nT2"
    },
    "bendover_length": {
      "value": "<slab_bendover_length_AU>",
      "unit": "AU"
    },
    "inertial_index": "<s>",
    "normalization": "one_sided_total_transverse_variance",
    "wavenumber_unit": "rad/m",
    "high_k_policy": "specified_analytic_tail"
  },
  "pitch_angle_regularization": "none"
}
~~~

~~~json
{
  "model": "nlgce_f_2014",
  "spectrum_provider": "local_slab_2d_state",
  "coefficient_set": "Qin_Zhang_2014_Tables_3_4",
  "logarithm": "natural",
  "out_of_domain": "error",
  "return_perpendicular_pair": true,
  "report_surrogate_diagnostics": true,
  "inside_domain_backend_switch": "disabled"
}
~~~

## 15. Verification requirements and numerical references

### 15.1 Minimal verification groups

These tests are intended for the future implementation; this document does not claim that the library has already been coded.

| Group | Required evidence |
|---|---|
| Units and species | Correct relativistic speed, rigidity, gyrofrequency, and SI conversions; ions and electrons |
| Constant/power laws | Equation (4); reference normalization; expected slopes including beta |
| Broken rigidity | H(R0)=1; a=b limit; low/high slopes; derivative Equation (19) |
| Pitch-angle integral | Exact isotropic result; target-lambda recovery for Equation (23) |
| Spectrum normalization | Integrated variance; per-component and one-/two-sided conversions |
| QLT | Equation (30), inertial limit, high-r_* limit, divergent-tail detection |
| Broadening | Kernel area pi; agreement with QLT as width tends to zero where both are regular |
| Nonlinear pair | Positive converged roots, dimensionless residuals, independent starting values |
| NLGCE-F | All 576 values and checksums, exact input ratios, natural logs, bounds, derivatives, source-parameter and small-epsilon regression fixtures |
| Geometry | Tensor projections; radial/parallel distinction; shock-normal coefficient; field-direction terms in Equation (62a) |
| Batch and optional inputs | Scalar/batch agreement; preserved per-point status; absent versus explicit-zero fields/derivatives |
| Moment conventions | Half-sum/unhalved Elsässer conversion; residual-energy sign; identical magnetic variance after conversion |
| Validity diagnostics | Strong-turbulence and surrogate-error flags remain separate from numerical success |
| State updates | Coefficients change when field or turbulence revisions change |
| Invalid states | Missing quantities, null field, absent scattering, fit-domain violations |

Match a coefficient table to its published fit at high numerical precision, but do not demand that its values equal the underlying nonlinear model at that same precision.

### 15.2 Relativistic/unit example

For a proton of total kinetic energy 10 MeV, using m_p=1.67262192369×10^-27 kg, c=299792458 m/s, and e=1.602176634×10^-19 C:

| Quantity | Reference value |
|---|---:|
| gamma | 1.01065788925 |
| beta | 0.144844004126 |
| v | 4.34231400235×10⁷ m/s |
| Rigidity | 1.37351526250×10⁸ V |
| lambda_parallel=0.1 AU | 1.495978707×10¹⁰ m (exact for the defined input) |
| kappa_parallel | 2.16533642887×10¹⁷ m²/s |

The proton mass is the [CODATA 2018 value archived by NIST](https://physics.nist.gov/cuu/Constants/ArchiveASCII/allascii_2018.txt), deliberately fixed for this regression example. The values of c and e are exact SI defining constants. These are calculated reference values for the defined 10 MeV and 0.1 AU inputs, not measurements or calibrated model parameters. The table rounds computed quantities to twelve significant digits; benchmark_points.json retains numerical reference values for comparisons at the declared tolerances. A later implementation using a different documented proton mass must regenerate the affected references.

### 15.3 Spectrum and pitch-angle constants

For the turbulence spectrum, s=5/3; for the regularized pitch-angle shape, q_mu=5/3:

| Quantity | Reference value |
|---|---:|
| C(5/6) | 0.118862354635 |
| L_c / ell_s | 0.746834200222 |
| Inertial prefactor in Equation (31) | 5.16465750497 |
| I(5/3,0.01) in Equation (23) | 4.2719856556 |
| D0 lambda_target / v for that shape | 1.60199462085 |

These computed values were checked by gamma-function evaluation and numerical quadrature; the regularized-shape integral was also checked with an independent closed-form expression. Displayed computed values are rounded to twelve significant digits. The explicitly defined indices and h_mu=0.01 are verification inputs, not inferred physical parameters.

For the exact smooth-spectrum QLT expression, Equation (30), with epsilon_s²=0.04:

| r_L / ell_s | lambda_parallel / ell_s |
|---:|---:|
| 0.001 | 12.9116445901 |
| 0.01 | 27.8174715425 |
| 1 | 137.17827805 |
| 10 | 1472.16829464 |

No 90-degree regularization or dissipation cutoff is applied to this table. It is the mathematical analytic-spectrum model with s=5/3 continued to arbitrarily high k.


### 15.4 NLGCE-F and nonlinear reference points

Both backends use nu=5/6, zero dynamical decorrelation and total epsilon². The nonlinear backend uses the multiplicative a_x of Equation (45) and fixes L_c=2 pi C(nu) ell_s. Points D–F test smaller turbulence levels; G is intermediate. These are defined numerical verification inputs, not measured event parameters.

| Point | r_L/ell_s | f_s | Total epsilon² | ell_s/ell_2 |
|---|---:|---:|---:|---:|
| A | 0.035848041610664981 | 0.2 | 1 | 10 |
| B | 0.01 | 0.2 | 0.25 | 10 |
| C | 0.1 | 0.5 | 1 | 1 |
| D | 0.01 | 0.2 | 0.01 | 10 |
| E | 0.001 | 0.2 | 0.001 | 10 |
| F | 0.01 | 0.2 | 0.0001 | 10 |
| G | 0.01 | 0.2 | 0.05 | 10 |

Point A corresponds to r_L/L_c=0.048 under Equation (28). E is an independently defined state, not a reconstruction of the review's incompletely specified M4 point.

| Point | Nonlinear lambda_parallel/ell_s | Fitted lambda_parallel/ell_s | Nonlinear lambda_perp/ell_s | Fitted lambda_perp/ell_s |
|---|---:|---:|---:|---:|
| A | 1.48799747985 | 1.57155490251 | 0.0543701824316 | 0.0523957778321 |
| B | 2.75584216326 | 3.20651578525 | 0.0261431510994 | 0.026524985621 |
| C | 2.15782150015 | 2.18084256045 | 0.0803778936919 | 0.0768736859119 |
| D | 107.219513777 | 97.8649318552 | 0.0128619686613 | 0.0116877324293 |
| E | 249.130915739 | 221.678269737 | 0.00360439951257 | 0.0030214818976 |
| F | 87716.3439308 | 95850.8828433 | 0.00603380578468 | 0.00682965096051 |
| G | 13.0325174575 | 15.5761048615 | 0.016749170539 | 0.016869816459 |

| Point | a_x | a'^2 | Signed parallel difference, 100(F/N-1) | Signed perpendicular difference, 100(F/N-1) |
|---|---:|---:|---:|---:|
| A | 0.162668319438 | 0.307900211697 | +5.615% | -3.631% |
| B | 0.167891402775 | 0.307900211697 | +16.353% | +1.461% |
| C | 0.344832514338 | 0.214285714286 | +1.067% | -4.360% |
| D | 0.103935583745 | 0.307900211697 | -8.725% | -9.130% |
| E | 0.0863571642202 | 0.307900211697 | -11.019% | -16.172% |
| F | 0.0236522188682 | 0.307900211697 | +9.274% | +13.190% |
| G | 0.181382697692 | 0.307900211697 | +19.517% | +0.720% |

These replace all affected nonlinear values in previous specification versions. A and C have epsilon=1 and are unchanged by the a_x correction; they cannot detect the multiplication-versus-exponentiation mistake. B and D–G must be included in source-parameter regression checks. Polynomial values are unchanged.

The polynomial uses the printed decimal coefficients and is independently checked with fifty-digit decimal arithmetic. The integral results solve Equations (48),(50),(51) in positive logarithmic unknowns, with maximum absolute logarithmic residual below 10^-8. A second bounded-wavenumber quadrature, k ell_a=z/(1-z), uses logarithmic kappa unknowns and an alternate initial pair; its measured differences are supplied in the companion summary.

For D–F, increasing the log-wavenumber range, tightening quadrature and using an initial guess 1000 times smaller reproduces both coefficients to better than 10^-12 relative. These are checks of the selected roots, not a proof of global uniqueness or physical accuracy.

**Regression rule:** compare each implementation against its own backend fixture and its declared numerical tolerance. Numerical agreement of the two backends is a separate approximation check. A nonlinear evaluator must solve the correct equations and cannot report the polynomial as a converged integral solution. In particular, verify the multiplicative factor in Equation (45) at epsilon different from one.

The printed derived values are rounded to twelve significant digits; benchmark_points.json retains full numerical fixtures and comparison tolerances. fit_error_samples.csv and fit_error_summary.json contain the regenerated 300-state audit. These are numerical reproducibility assets, not independent observational validation.

### 15.5 Meaningful transport checks

Once a solver uses the library, coefficient verification should be supplemented by transport checks that test its actual integration:

- homogeneous diffusion should give parallel variance 2 kappa_parallel t;
- a manufactured varying-coefficient solution should test the diffusion drift/divergence term;
- a Parker-spiral comparison should test Equation (5), including perpendicular diffusion;
- a slowly varying focused-transport solution should approach the corresponding Parker result as the diffusion limit improves;
- a shock test should use the normal coefficient on each side, rather than a scalar parallel value regardless of obliquity.

Observed SEP or GCR spectra constrain a combination of scattering, sources/boundaries, connectivity, convection, and other transport processes. Matching one flux profile cannot uniquely validate the coefficient in isolation.

## 16. Model selection and scientific interpretation

For controlled SEP comparisons, constant lambda and explicit power laws provide transparent assumptions. For long-term GCR modulation fits, a smooth broken rigidity model permits direct comparison with PAMELA/AMS-era work. For a background with credible turbulence moments, QLT and NLGCE provide a route from local fluctuation properties to coefficients.

A broader synthesis is provided by Engelbrecht et al. (2022) [R25]. The last two decades retain all of these approaches. Nonlinear closures did not eliminate empirical laws, and empirical fits did not establish a universal turbulence spectrum. The useful distinction is what data determine the coefficient and which assumptions remain free.

| Scientific setting | Model family to compare | Main limitation to report |
|---|---|---|
| SEP parameter study | Constant lambda and rigidity power laws | Assumed scattering and early anisotropic transport |
| GCR solar-cycle modulation | Smooth broken rigidity laws | Fitted rigidity/time dependence and boundary spectrum |
| Background with slab turbulence | QLT spectrum and broadening alternatives | Slab allocation, high-k power, 90-degree treatment |
| Composite slab/2D turbulence | NLPA/NLGCE-N and NLGCE-F | Spectral geometry, nonlinear closure, fit error/domain |
| Shock-wave feedback available | Wave provider plus compatible scattering closure | Wave growth, polarization, frame, and momentum coupling |
| Strong-scattering comparison | Bohm prescription | Assumed field choice and lack of universal-bound status |

Recent work reinforces limitations relevant to a reusable library. Shalchi and Klippenstein (2025) discuss nonuniform-field definitions [R26]; Shalchi (2026) develops pitch-angle scattering beyond standard QLT [R27]. These papers should guide future paper-specific extensions. Their full theories are not represented by renaming the generic kernels in Section 9.

Alonso Guzmán et al. (2025) provide a recent observational study of turbulence geometry in stream interaction regions, relevant to choosing slab/2D inputs [R28]. A fixed slab fraction may be a useful assumption, but it should remain an explicit model parameter rather than a universal property of every heliospheric region.

A staged coding order is: explicit constant/power/broken laws and unit handling; normalized pitch-angle integrals and QLT; separate NLGCE-F and nonlinear backends with the complete data and small-epsilon source-transcription fixtures; broadened spectra and provider adapters. The polynomial remains computationally convenient, but that convenience does not establish equivalence with the integral backend in the discrepant regimes. This ordering is a proposed development sequence, not a statement that the scientific theories must be used in that order.

## 17. Publication methods text

The paragraphs below are original methods prose. Use only paragraphs for models actually activated, retain the equations and references, and replace the descriptive parameter choices with the values used in the study.

The document identifiers [Rxx] are internal cross-references. Replace them with the canonical manuscript citation keys below, or a verified equivalent key in the study's bibliography. Do not cite a filename-derived year or reuse an entry with conflicting DOI/title/journal metadata.

### 17.0 Citation conversion and parameter reporting

| Reference | Canonical suggested manuscript key |
|---|---|
| R01 | Parker-1965-PSS |
| R02 | Jokipii-1966-AJ |
| R03 | Bieber-1994-AJ-Palmer |
| R04 | Teufel-2003-AA |
| R05 | Matthaeus-2003-AJ |
| R06 | Shalchi-2004-AJ-WNLT |
| R07 | Zank-1998-JGR |
| R08 | Qin-2007-AJ |
| R09 | Qin-2013-AJ-Erratum |
| R10 | Qin-2014-AJ-NLGCE |
| R11 | Shalchi-2009-AA |
| R12 | Engelbrecht-2013-AJ-772 |
| R13 | Engelbrecht-2013-AJ-779 |
| R14 | Potgieter-2014-SP-391 |
| R15 | Vos-2015-AJ |
| R16 | Qin-2017-AJ |
| R17 | Chhiber-2017-AJSS |
| R18 | Corti-2019-AJ |
| R19 | Wijsen-2023-JGR |
| R20 | Vainio-2007-AJ |
| R21 | Afanasiev-2013-ApJS |
| R22 | Afanasiev-2015-AA |
| R23 | Zhao-2018-AJ |
| R24 | Shen-2018-AJ |
| R25 | Engelbrecht-2022-SSRv |
| R26 | Shalchi-2025-AJ |
| R27 | Shalchi-2026-AJ |
| R28 | Alonso-Guzman-2025-JGR |
| R29 | Hussein-2014-AJ-Bohm |
| R30 | Droege-2014-JGR |
| R31 | Marcowith-1999-arXiv |

The R03 key denotes the verified 1994 paper, not the conflicting supplied Bieber entry. R22 uses the verified 2015 journal identity; its 2016 arXiv version is a legitimate alias.

Complete the applicable reporting rows before using the methods prose in a publication. Reference values below are numerical verification states, not publication parameter defaults.

| Reporting item | Values and conventions to supply |
|---|---|
| Transport formulation | Parker or focused transport; evolved distribution; numerical frame and independent variables |
| Particle/background state | Species mass and charge state; energy/rigidity interval; B0 and its averaging scale; provider and revisions |
| Constant/power model | Quantity prescribed; lambda0 or kappa0; reference scales; exponents; spatial/time factors |
| Broken rigidity model | K_star convention and units; low/high slopes; smoothness; break rigidity; field/time factors |
| Pitch-angle shape | q_mu, h_mu, and target-lambda versus fixed-amplitude mode; amplitude provider and dependencies |
| Spectrum/QLT/broadening | One-/two-sided and component normalization; variance; spectral geometry; lengths/slopes/tails; decorrelation kernel/width |
| Nonlinear pair | Exact backend, spectrum, epsilon², f_s, ell_s/ell_2, r_L/ell_s, slab L_c definition, solver residual and convergence checks |
| NLGCE-F | Tables 3/4 and their digests; actual input ranges; out-of-box policy; in-box surrogate diagnostics and any application-specific comparison |
| Turbulence moments | Named Elsässer/energy convention, residual-energy sign, density normalization, and conversion used |
| Tensor/SDE geometry | Perpendicular prescription; coefficient and field-direction derivatives; handling of unavailable derivatives |
| Limitations | Diffusion-limit and strong-turbulence diagnostics; numerical reproduction versus closure/observational validation |

### 17.1 General library description

We describe nearly isotropic energetic-particle transport using the Parker equation [R01]. The symmetric spatial diffusion tensor is expressed in a local field-aligned basis, with parallel coefficient kappa_parallel and separately specified perpendicular coefficients. The parallel mean free path is defined by lambda_parallel=3 kappa_parallel/v, where the particle speed is calculated relativistically from its total kinetic energy and rest mass. Coefficients are evaluated using the local background and turbulence state. The parallel, radial, and shock-normal coefficients are distinguished through tensor projection rather than treated as interchangeable scalar quantities.

### 17.2 Prescribed scattering model

The prescribed parallel mean free path is represented by a normalized power law in rigidity and heliocentric radius, with optional explicitly specified dependence on the mean field and time. The spatial diffusion coefficient follows from kappa_parallel=v lambda_parallel/3. The parameterization describes assumed scattering rather than a self-consistent turbulence calculation. Its rigidity exponent refers to the mean free path; the diffusion coefficient includes the additional rigidity dependence of particle speed.

### 17.3 Broken rigidity model

For the empirical modulation model, kappa_parallel is proportional to beta times a smooth broken rigidity function and an explicitly defined field-strength factor. The low- and high-rigidity indices describe the coefficient after removing beta. The normalization is expressed in SI units and its reference state is specified. This prescription follows the class of Parker-equation modulation coefficients employed by Potgieter et al. (2014), Vos and Potgieter (2015), and Corti et al. (2019) [R14,R15,R18]. Its parameter values are determined for the study rather than adopted as universal scattering constants.

### 17.4 QLT spectrum model

The slab scattering calculation uses a one-sided magnetic spectrum normalized to the total transverse slab variance. The pitch-angle diffusion coefficient is evaluated at the gyroresonant wavenumber using Equation (26), and the parallel mean free path is obtained from Equation (20). The spectral bend-over length, inertial index, and treatment of the high-wavenumber range are specified independently. The adopted treatment of pitch angles near 90 degrees is reported explicitly. The inertial-range approximation is used only when its resonant and gyroradius conditions are satisfied [R02,R03,R04,R07].

### 17.5 Nonlinear and fitted models

For the integral composite-turbulence calculation, the slab and two-dimensional magnetic variances and bend-over lengths define the local turbulence state. We solve the corrected NLPA relation coupled to the modified perpendicular relation [R08,R09,R10], using the slab correlation length L_c=2 pi C(nu) ell_s. Total magnetic variance defines epsilon². The pair is accepted according to its integral residual and quadrature convergence, with the numerical settings and input ranges reported for the study.

For calculations using NLGCE-F, we evaluate the printed Qin and Zhang (2014) polynomial in the four dimensionless variables of Equation (52), using both coefficient tables and their published input restrictions [R10]. This is the polynomial backend, whose approximation quality is distinct from numerical reproduction of its coefficients. The approximation discrepancy is assessed against the correctly transcribed integral closure at the actual input states; the published box does not itself guarantee uniform approximation accuracy. The study must report its actual parameter ranges, surrogate diagnostics, and any comparison performed for those states. Backend changes, bounds, or additional calibration are identified explicitly.

These two paragraphs describe different computational options and must be selected or qualified to match the calculation actually performed. No universal percentage error or equivalence of the backends is implied.

### 17.6 Transported turbulence and limitations

When turbulence moments or wave spectra evolve, the coefficient evaluator uses their declared magnetic-energy conventions and current local values. The turbulence provider and microscopic scattering closure are identified separately. This follows the distinction between turbulence evolution and energetic-particle scattering used in turbulence-based heliospheric studies [R12,R13,R17,R19,R23]. The Parker approximation is interpreted as a diffusion-limit description; periods with pronounced anisotropy or rapid focusing are evaluated with appropriate validity diagnostics rather than assumed valid solely because the computed coefficient is finite.

### 17.7 Prescribed pitch-angle operator

When focused transport is evolved, the pitch-angle operator uses the regularized shape D_mu_mu=D0(1-mu²)(|mu|^(q_mu-1)+h_mu), consistent with the class of SEP prescriptions used by Dröge et al. (2014) [R30]. In target-mean-free-path mode, D0 is set by Equation (23), so the integral definition of lambda_parallel recovers the prescribed target. In fixed-amplitude mode, the mean free path is derived from the chosen D0 and changes when the shape parameters change. The selected mode, shape parameters, and spatial/energy dependence must be reported; the focused solver does not add a second parallel spatial operator representing the same scattering.

## 18. References

References are listed below with persistent DOI links where available. They include pre-2006 foundations needed to define the models precisely and representative 2006–2026 developments/applications. This is a focused model bibliography, not an exhaustive census of every SEP/GCR modeling paper.

Each reference also gives a suggested BibTeX citation key for translating the document's [Rxx] identifiers into manuscript citations. For existing entries, retain the matching key in the supplied bibliography. The Bieber key is a proposed corrected identity and is not supplied as an automatic addition.

**[R01]** Parker, E. N. (1965). **The passage of energetic charged particles through interplanetary space.** *Planetary and Space Science*, **13**, 9–49. [DOI](https://doi.org/10.1016/0032-0633(65)90131-5). Role: Foundation of Parker transport. Suggested citation key: **Parker-1965-PSS**.

**[R02]** Jokipii, J. R. (1966). **Cosmic-Ray Propagation. I. Charged Particles in a Random Magnetic Field.** *The Astrophysical Journal*, **146**, 480. [DOI](https://doi.org/10.1086/148912). Role: Resonant quasi-linear scattering. Suggested citation key: **Jokipii-1966-AJ**.

**[R03]** Bieber, J. W.; Matthaeus, W. H.; Smith, C. W.; Wanner, W.; Kallenrode, M.-B.; Wibberenz, G. (1994). **Proton and Electron Mean Free Paths: The Palmer Consensus Revisited.** *The Astrophysical Journal*, **420**, 294. [DOI](https://doi.org/10.1086/173559). Role: Slab/composite spectra and heliospheric mean free paths. Suggested citation key: **Bieber-1994-AJ-Palmer**.

**[R04]** Teufel, A.; Schlickeiser, R. (2003). **Analytic calculation of the parallel mean free path of heliospheric cosmic rays. II. Dynamical magnetic slab turbulence and random sweeping slab turbulence with finite wave power at small wavenumbers.** *Astronomy & Astrophysics*, **397**, 15–25. [DOI](https://doi.org/10.1051/0004-6361:20021471). Role: Dynamical slab scattering and mean-free-path calculations. Suggested citation key: **Teufel-2003-AA**.

**[R05]** Matthaeus, W. H.; Qin, G.; Bieber, J. W.; Zank, G. P. (2003). **Nonlinear Collisionless Perpendicular Diffusion of Charged Particles.** *The Astrophysical Journal*, **590**, L53–L56. [DOI](https://doi.org/10.1086/376613). Role: Perpendicular NLGC closure needed by nonlinear parallel pairs. Suggested citation key: **Matthaeus-2003-AJ**.

**[R06]** Shalchi, A.; Bieber, J. W.; Matthaeus, W. H.; Qin, G. (2004). **Nonlinear Parallel and Perpendicular Diffusion of Charged Cosmic Rays in Weak Turbulence.** *The Astrophysical Journal*, **616**, 617–629. [DOI](https://doi.org/10.1086/424839). Role: Original weakly nonlinear composite-turbulence theory. Suggested citation key: **Shalchi-2004-AJ-WNLT**.

**[R07]** Zank, G. P.; Matthaeus, W. H.; Bieber, J. W.; Moraal, H. (1998). **The Radial and Latitudinal Dependence of the Cosmic Ray Diffusion Tensor in the Heliosphere.** *Journal of Geophysical Research*, **103**, 2085–2097. [DOI](https://doi.org/10.1029/97JA03013). Role: Turbulence-dependent spatial coefficients and spectral scaling. Suggested citation key: **Zank-1998-JGR**.

**[R08]** Qin, G. (2007). **Nonlinear Parallel Diffusion of Charged Particles: Extension to the Nonlinear Guiding Center Theory.** *The Astrophysical Journal*, **656**, 217–221. [DOI](https://doi.org/10.1086/510510). Role: NLPA and original coupled NLGC-E. Suggested citation key: **Qin-2007-AJ**.

**[R09]** Qin, G. (2013). **Erratum: Nonlinear Parallel Diffusion of Charged Particles: Extension to the Nonlinear Guiding Center Theory (2007, ApJ, 656, 217).** *The Astrophysical Journal*, **774**, 91. [DOI](https://doi.org/10.1088/0004-637X/774/1/91). Role: Required correction to the NLPA parameter. Suggested citation key: **Qin-2013-AJ-Erratum**.

**[R10]** Qin, G.; Zhang, L.-H. (2014). **The Modification of the Nonlinear Guiding Center Theory.** *The Astrophysical Journal*, **787**, 12. [DOI](https://doi.org/10.1088/0004-637X/787/1/12); [author preprint](https://arxiv.org/abs/1401.1950). Role: NLGCE-N, NLGCE-F, coefficients, and parameter domain. Suggested citation key: **Qin-2014-AJ-NLGCE**.

**[R11]** Shalchi, A. (2009). **Analytical description of nonlinear cosmic ray scattering: isotropic and quasilinear regimes of pitch-angle diffusion.** *Astronomy & Astrophysics*, **507**, 589–597. [DOI](https://doi.org/10.1051/0004-6361/200912755). Role: Nonlinear pitch-angle scattering and regime transitions. Suggested citation key: **Shalchi-2009-AA**.

**[R12]** Engelbrecht, N. E.; Burger, R. A. (2013). **An Ab Initio Model for Cosmic-Ray Modulation.** *The Astrophysical Journal*, **772**, 46. [DOI](https://doi.org/10.1088/0004-637X/772/1/46). Role: Turbulence-transport-derived GCR coefficients. Suggested citation key: **Engelbrecht-2013-AJ-772**.

**[R13]** Engelbrecht, N. E.; Burger, R. A. (2013). **An Ab Initio Model for the Modulation of Galactic Cosmic-Ray Electrons.** *The Astrophysical Journal*, **779**, 158. [DOI](https://doi.org/10.1088/0004-637X/779/2/158). Role: Electron scattering, dissipation-range and dynamical-turbulence sensitivity. Suggested citation key: **Engelbrecht-2013-AJ-779**.

**[R14]** Potgieter, M. S.; Vos, E. E.; Boezio, M.; De Simone, N.; Di Felice, V.; Formato, V. (2014). **Modulation of Galactic Protons in the Heliosphere During the Unusual Solar Minimum of 2006 to 2009.** *Solar Physics*, **289**, 391–406. [DOI](https://doi.org/10.1007/s11207-013-0324-6); [author preprint](https://arxiv.org/abs/1302.1284). Role: Smooth broken rigidity diffusion in a Parker model. Suggested citation key: **Potgieter-2014-SP-391**.

**[R15]** Vos, E. E.; Potgieter, M. S. (2015). **New Modeling of Galactic Proton Modulation During the Minimum of Solar Cycle 23/24.** *The Astrophysical Journal*, **815**, 119. [DOI](https://doi.org/10.1088/0004-637X/815/2/119). Role: PAMELA/Voyager-era empirical GCR modulation. Suggested citation key: **Vos-2015-AJ**.

**[R16]** Qin, G.; Shen, Z.-N. (2017). **Modulation of Galactic Cosmic Rays in the Inner Heliosphere, Comparing with PAMELA Measurements.** *The Astrophysical Journal*, **846**, 56. [DOI](https://doi.org/10.3847/1538-4357/aa83ad); [author preprint](https://arxiv.org/abs/1705.04847). Role: Application of NLGCE-F to GCR modulation. Suggested citation key: **Qin-2017-AJ**.

**[R17]** Chhiber, R.; Subedi, P.; Usmanov, A. V.; Matthaeus, W. H.; Ruffolo, D.; Goldstein, M. L.; Parashar, T. N. (2017). **Cosmic-Ray Diffusion Coefficients throughout the Inner Heliosphere from a Global Solar Wind Simulation.** *The Astrophysical Journal Supplement Series*, **230**, 21. [DOI](https://doi.org/10.3847/1538-4365/aa74d2); [author preprint](https://arxiv.org/abs/1703.10322). Role: Global evaluation of turbulence-dependent local mean free paths. Suggested citation key: **Chhiber-2017-AJSS**.

**[R18]** Corti, C.; Potgieter, M. S.; Bindi, V.; Consolandi, C.; Light, C.; Palermo, M.; Popkow, A. (2019). **Numerical Modeling of Galactic Cosmic-Ray Proton and Helium Observed by AMS-02 during the Solar Maximum of Solar Cycle 24.** *The Astrophysical Journal*, **871**, 253. [DOI](https://doi.org/10.3847/1538-4357/aafac4); [author preprint](https://arxiv.org/abs/1810.09640). Role: Broken rigidity coefficients and species-dependent speeds in Parker transport. Suggested citation key: **Corti-2019-AJ**.

**[R19]** Wijsen, N.; Li, G.; Ding, Z.; Lario, D.; Poedts, S.; Filwett, R. J.; Allen, R. C.; Dayeh, M. A. (2023). **On the Seed Population of Solar Energetic Particles in the Inner Heliosphere.** *Journal of Geophysical Research: Space Physics*, **128**, e2022JA031203. [DOI](https://doi.org/10.1029/2022JA031203). Role: SEP application with turbulence-dependent parallel mean free paths. Suggested citation key: **Wijsen-2023-JGR**.

**[R20]** Vainio, R.; Laitinen, T. (2007). **Monte Carlo Simulations of Coronal Diffusive Shock Acceleration in Self-generated Turbulence.** *The Astrophysical Journal*, **658**, 622–630. [DOI](https://doi.org/10.1086/510284). Role: SEP scattering with self-generated waves. Suggested citation key: **Vainio-2007-AJ**.

**[R21]** Afanasiev, A.; Vainio, R. (2013). **Monte Carlo Simulation Model of Energetic Proton Transport Through Self-Generated Alfvén Waves.** *The Astrophysical Journal Supplement Series*, **207**, 29. [DOI](https://doi.org/10.1088/0067-0049/207/2/29). Role: Particle–wave coupling and scattering-provider requirements. Suggested citation key: **Afanasiev-2013-ApJS**.

**[R22]** Afanasiev, A.; Battarbee, M.; Vainio, R. (2015). **Self-consistent Monte Carlo Simulations of Proton Acceleration in Coronal Shocks: Effect of Anisotropic Pitch-angle Scattering of Particles.** *Astronomy & Astrophysics*, **584**, A81. [DOI](https://doi.org/10.1051/0004-6361/201526750). Role: Anisotropic scattering in self-consistent shock acceleration. Suggested citation key: **Afanasiev-2015-AA**.

**[R23]** Zhao, L.-L.; Adhikari, L.; Zank, G. P.; Hu, Q.; Feng, X. S. (2018). **Influence of the Solar Cycle on Turbulence Properties and Cosmic-Ray Diffusion.** *The Astrophysical Journal*, **856**, 94. [DOI](https://doi.org/10.3847/1538-4357/aab362). Role: Observed turbulence dependence of GCR diffusion. Suggested citation key: **Zhao-2018-AJ**.

**[R24]** Shen, Z.-N.; Qin, G. (2018). **Modulation of Galactic Cosmic Rays in the Inner Heliosphere over Solar Cycles.** *The Astrophysical Journal*, **854**, 137. [DOI](https://doi.org/10.3847/1538-4357/aaab64); [author preprint](https://arxiv.org/abs/1709.08017). Role: Time-dependent turbulence-based GCR modulation. Suggested citation key: **Shen-2018-AJ**.

**[R25]** Engelbrecht, N. E.; Effenberger, F.; Florinski, V.; Potgieter, M. S.; Ruffolo, D.; Chhiber, R.; Usmanov, A. V.; Rankin, J. S.; Els, P. L. (2022). **Theory of Cosmic Ray Transport in the Heliosphere.** *Space Science Reviews*, **218**, 33. [DOI](https://doi.org/10.1007/s11214-022-00896-1). Role: Modern synthesis and context; not the source of the fitted coefficient arrays. Suggested citation key: **Engelbrecht-2022-SSRv**.

**[R26]** Shalchi, A.; Klippenstein, B. (2025). **The Parallel Diffusion Coefficient of Energetic Particles in a Spatially Varying Magnetic Field.** *The Astrophysical Journal*, **990**, 128. [DOI](https://doi.org/10.3847/1538-4357/adea53). Role: Limits of local nonuniform-field diffusion definitions. Suggested citation key: **Shalchi-2025-AJ**.

**[R27]** Shalchi, A. (2026). **Pitch-angle Scattering of Energetic Particles Interacting with MHD Turbulence. I. Analytical Theory.** *The Astrophysical Journal*, **1004**, 125. [DOI](https://doi.org/10.3847/1538-4357/ae70fc). Role: Recent development beyond standard QLT scattering. Suggested citation key: **Shalchi-2026-AJ**.

**[R28]** Alonso Guzmán, J. G.; Ghanbari, K.; Florinski, V. A.; Leske, R. A.; Zhao, L.-L.; Zhu, X.; et al. (2025). **Superposed Epoch Analysis of Stream Interaction Regions at 1 au During Solar Minimum With Turbulence Geometry Decomposition: Implications for Galactic Cosmic Ray Transport.** *Journal of Geophysical Research: Space Physics*, **130**, e2024JA033567. [DOI](https://doi.org/10.1029/2024JA033567). Role: Recent observational constraints on turbulence geometry. Suggested citation key: **Alonso-Guzman-2025-JGR**.

**[R29]** Hussein, M.; Shalchi, A. (2014). **Detailed Numerical Investigation of the Bohm Limit in Cosmic Ray Diffusion Theory.** *The Astrophysical Journal*, **785**, 31. [DOI](https://doi.org/10.1088/0004-637X/785/1/31). Role: Validity of Bohm scaling. Suggested citation key: **Hussein-2014-AJ-Bohm**.

**[R30]** Dröge, W.; Kartavykh, Y. Y.; Dresing, N.; Heber, B.; Klassen, A. (2014). **Wide longitudinal distribution of interplanetary electrons following the 7 February 2010 solar event: Observations and transport modeling.** *Journal of Geophysical Research: Space Physics*, **119**, 6074–6094. [DOI](https://doi.org/10.1002/2014JA019933). Role: Regularized pitch-angle prescriptions and focused SEP transport. Suggested citation key: **Droege-2014-JGR**.

**[R31]** Marcowith, A.; Kirk, J. G. (1999). **Computation of diffusive shock acceleration using stochastic differential equations.** [Author preprint](https://arxiv.org/abs/astro-ph/9905176). Role: First-order shock acceleration and the distinction from optional stochastic momentum-diffusion terms. Suggested citation key: **Marcowith-1999-arXiv**.

## 19. Bibliography comparison

**Snapshot appendix:** this audit describes the supplied file at the stated cutoff; it does not claim to track later edits. Publication citations must be verified against the primary identities in Section 18. This appendix is separable from the model specification and methods prose.

The comparison uses the supplied bibliography(20261008-011051).bib, containing 6,861 parsed bibliographic records. Its SHA-256 digest is:

~~~text
b9a1b85db18ee6a382e8bfb44b562d9069f62986c39665d7b50b295eafa829e5
~~~

“Not found” means no matching record was identified in this supplied snapshot after DOI, normalized title, and bibliographic-identity checks. It does not imply that the paper is absent from every version of the user's bibliography. Matching titles with a missing DOI are reported as present. Conflicting candidate metadata are identified separately and are not treated as safe new additions.

### 19.1 Relevant papers not found in the supplied bibliography

**12 references from 2006–2026 were not found** in the supplied snapshot:

| Reference | Paper | Why it is needed |
|---|---|---|
| R08 | Qin, G. (2007) — Nonlinear Parallel Diffusion of Charged Particles: Extension to the Nonlinear Guiding Center Theory | NLPA and original coupled NLGC-E |
| R09 | Qin, G. (2013) — Erratum: Nonlinear Parallel Diffusion of Charged Particles: Extension to the Nonlinear Guiding Center Theory (2007, ApJ, 656, 217) | Required correction to the NLPA parameter |
| R10 | Qin, G. et al. (2014) — The Modification of the Nonlinear Guiding Center Theory | NLGCE-N, NLGCE-F, coefficients, and parameter domain |
| R12 | Engelbrecht, N. E. et al. (2013) — An Ab Initio Model for Cosmic-Ray Modulation | Turbulence-transport-derived GCR coefficients |
| R13 | Engelbrecht, N. E. et al. (2013) — An Ab Initio Model for the Modulation of Galactic Cosmic-Ray Electrons | Electron scattering, dissipation-range and dynamical-turbulence sensitivity |
| R14 | Potgieter, M. S. et al. (2014) — Modulation of Galactic Protons in the Heliosphere During the Unusual Solar Minimum of 2006 to 2009 | Smooth broken rigidity diffusion in a Parker model |
| R15 | Vos, E. E. et al. (2015) — New Modeling of Galactic Proton Modulation During the Minimum of Solar Cycle 23/24 | PAMELA/Voyager-era empirical GCR modulation |
| R17 | Chhiber, R. et al. (2017) — Cosmic-Ray Diffusion Coefficients throughout the Inner Heliosphere from a Global Solar Wind Simulation | Global evaluation of turbulence-dependent local mean free paths |
| R18 | Corti, C. et al. (2019) — Numerical Modeling of Galactic Cosmic-Ray Proton and Helium Observed by AMS-02 during the Solar Maximum of Solar Cycle 24 | Broken rigidity coefficients and species-dependent speeds in Parker transport |
| R28 | Alonso Guzmán, J. G. et al. (2025) — Superposed Epoch Analysis of Stream Interaction Regions at 1 au During Solar Minimum With Turbulence Geometry Decomposition: Implications for Galactic Cosmic Ray Transport | Recent observational constraints on turbulence geometry |
| R29 | Hussein, M. et al. (2014) — Detailed Numerical Investigation of the Bohm Limit in Cosmic Ray Diffusion Theory | Validity of Bohm scaling |
| R30 | Dröge et al. (2014) — Wide longitudinal distribution of interplanetary electrons following the 7 February 2010 solar event | Explicit pitch-angle prescription attribution |

The following **pre-2006 foundations** are additional missing references; they are listed separately from the requested twenty-year interval:

| Reference | Paper | Why it is needed |
|---|---|---|
| R04 | Teufel, A. et al. (2003) — Analytic calculation of the parallel mean free path of heliospheric cosmic rays. II. Dynamical magnetic slab turbulence and random sweeping slab turbulence with finite wave power at small wavenumbers | Dynamical slab scattering and mean-free-path calculations |
| R06 | Shalchi, A. et al. (2004) — Nonlinear Parallel and Perpendicular Diffusion of Charged Cosmic Rays in Weak Turbulence | Original weakly nonlinear composite-turbulence theory |

### 19.2 Found references and metadata issues

| Reference | Status in the supplied file | Matching citation keys |
|---|---|---|
| R01 | present | parker65; Parker-1965-PSS |
| R02 | present | Jokipii-1966-AJ |
| R03 | conflicting candidate | Bieber-1994-AJ |
| R05 | present | Matthaeus-2003-AJ |
| R07 | present | Zank-1998-JGR |
| R11 | present | Shalchi-2009-AA |
| R16 | present | Qin-2017-AJ |
| R19 | present | Wijsen-2023-arXiv; Wijsen-2023-JGR |
| R20 | present | Vainio-2007-AJ |
| R21 | present | Afanasiev-2013-ApJS |
| R22 | verified journal entry and arXiv alias; conflicting/unresolved title candidates excluded | Afanasiev-2015-AA; Afanasiev-2016-arXiv |
| R23 | present | Zhao-2018-AJ |
| R24 | present | Shen-2018-AJ |
| R25 | present | Engelbrecht-2022-SSRv |
| R26 | present | Shalchi-2025-AJ |
| R27 | present; DOI and June 2026 publication verified | Shalchi-2026-AJ |
| R31 | present as author preprint | Marcowith-1999-arXiv |

The supplied bibliography has a conflicting Bieber-1994-AJ record: its title and three-author list correspond to the 1996 two-dimensional-turbulence paper, while its year, journal, volume, and pages resemble the 1994 mean-free-path paper. The correct identities are:

- Bieber et al. (1994), **Proton and Electron Mean Free Paths: The Palmer Consensus Revisited**, ApJ 420, 294–306, DOI 10.1086/173559 [R03].
- Bieber, Wanner, and Matthaeus (1996), **Dominant Two-dimensional Solar Wind Turbulence with Implications for Cosmic Ray Transport**, JGR 101, 2511–2522, DOI 10.1029/95JA02588.

The DOI 10.1086/173566 in that supplied record is not the DOI of the 1994 paper cited here. Review the conflicting entry before adding either intended paper.

The supplied Shalchi-2025-AJ record has the correct DOI and authors, but gives pages 1–14; the published citation uses article **128** in ApJ **990** [R26]. The Wijsen-2023-JGR record corresponds to DOI **10.1029/2022JA031203** [R19]; a separate arXiv-titled record in the supplied bibliography carries a different DOI and should be checked before reuse.

For R22, Afanasiev-2015-AA has the correct journal identity: A&A 584, A81 (2015), DOI 10.1051/0004-6361/201526750. The author preprint arXiv:1603.08857 was submitted in 2016 and explicitly identifies that same 2015 journal publication. It is a valid alias, not a different 2016 journal paper.

Afanasiev-2018-AA-618 is a conflicting candidate: its title matches the 2015 study, but A&A 618, A114 (2018) identifies Plaschke et al., First observations of magnetic holes deep within the coma of a comet, DOI 10.1051/0004-6361/201833300. Do not treat that entry as a verified match until corrected. Afanasiev-misc is an unresolved title-only candidate; prefer the verified canonical key for manuscript citations.

R27 is a real publication with DOI 10.3847/1538-4357/ae70fc, online publication on 9 June 2026 and print publication on 10 June 2026. Its date is within this document's literature cutoff.

These observations do not modify the supplied bibliography. Section 21 provides suggested entries only for unambiguous missing references.

## 20. NLGCE-F coefficient tables

The data below reproduce the numerical coefficients in Qin and Zhang (2014), Tables 3 and 4 [R10], as retrieved from arXiv:1401.1950v2. Each row contains j,k,l followed by d_0jkl through d_5jkl. Use these as **numeric data**, preserving the indices and scientific notation. Both tables have 48 rows and 288 coefficients.

No coefficient has been re-fit. All 576 values were compared individually against the printed primary-source tables, rather than inferred from three benchmark evaluations. The companion CSV files preserve the exact indexed decimal data; their SHA-256 checksums are recorded in Section 22 and SHA256SUMS. The mathematical checks in Section 15 use these values. Cite the original paper when using this fit.

### 20.1 Parallel coefficients

~~~csv
j,k,l,d_i0,d_i1,d_i2,d_i3,d_i4,d_i5
0,0,0,0.23553875E+01,0.13029025E+01,0.13427579E+00,-0.20846185E-01,-0.43999458E-02,-0.19732612E-03
0,0,1,-0.66243965E-02,0.77768482E-01,0.41205730E-01,0.90071283E-02,0.90300755E-03,0.33110445E-04
0,0,2,0.21922888E-02,-0.39684652E-02,-0.17124410E-02,-0.94109440E-04,-0.11972348E-06,-0.45759483E-07
0,1,0,-0.14247343E+01,-0.17294922E-01,0.50330910E-01,0.22555091E-02,-0.54744446E-03,-0.36329187E-04
0,1,1,0.29375712E-01,-0.10262387E-01,-0.11404086E-01,-0.16183987E-02,-0.55537309E-04,0.92677294E-06
0,1,2,-0.10823704E-01,-0.31550348E-02,0.10669635E-02,0.34794916E-03,0.30278040E-04,0.82627464E-06
0,2,0,0.96583711E-01,-0.12322827E-01,-0.14806351E-01,-0.19674998E-02,-0.47832394E-04,0.22621783E-05
0,2,1,-0.16641267E-01,-0.38108147E-02,0.15277614E-02,0.53316833E-03,0.49867082E-04,0.14719990E-05
0,2,2,0.24760055E-02,0.75282162E-03,-0.23891975E-03,-0.11721492E-03,-0.13940055E-04,-0.51898963E-06
0,3,0,0.11152796E-01,0.31889965E-03,-0.23482006E-02,-0.47383322E-03,-0.30607870E-04,-0.57224705E-06
0,3,1,-0.17297051E-02,-0.54967106E-03,0.19481087E-03,0.72844701E-04,0.72128259E-05,0.22722512E-06
0,3,2,0.28333003E-03,0.10946972E-03,-0.23059951E-04,-0.12751296E-04,-0.15280551E-05,-0.56788270E-07
1,0,0,-0.14786874E+00,0.35413822E+00,-0.23571521E-01,-0.29883137E-01,-0.41108264E-02,-0.16630590E-03
1,0,1,-0.33315220E+00,-0.19457101E+00,0.48492938E-01,0.19439119E-01,0.20217803E-02,0.68237390E-04
1,0,2,-0.70191117E-02,0.16863137E-01,-0.42473573E-02,-0.20407538E-02,-0.22445321E-03,-0.77342338E-05
1,1,0,-0.12645333E+00,-0.57381385E-01,0.19660793E-01,0.81465626E-02,0.93001667E-03,0.34019143E-04
1,1,1,-0.74499928E-01,0.21886029E-01,0.78026833E-02,-0.90778341E-03,-0.27852845E-03,-0.13910525E-04
1,1,2,0.88176328E-02,-0.24607326E-02,-0.16399476E-02,-0.65717348E-04,0.18575402E-04,0.12449253E-05
1,2,0,0.32714354E-01,-0.96871348E-02,-0.88372387E-02,-0.10864091E-02,0.19497311E-06,0.31551668E-05
1,2,1,-0.87371981E-02,0.23907528E-02,-0.84712003E-03,-0.58296444E-03,-0.80491321E-04,-0.32954014E-05
1,2,2,0.20518127E-02,0.25128221E-03,-0.21217125E-04,0.70656292E-07,0.91585353E-06,0.56694614E-07
1,3,0,0.38644947E-02,0.50051530E-03,-0.11014859E-02,-0.31141696E-03,-0.28640974E-04,-0.86437844E-06
1,3,1,0.99272926E-04,0.56214724E-04,-0.18481718E-03,-0.55454476E-04,-0.54975212E-05,-0.18100547E-06
1,3,2,0.87113797E-04,0.41654101E-04,0.12554768E-04,0.57474926E-06,-0.11142497E-06,-0.82310051E-08
2,0,0,0.12726928E+00,0.60497042E-01,-0.11952481E-01,-0.71946142E-02,-0.90296336E-03,-0.35119674E-04
2,0,1,-0.12984315E+00,-0.42695185E-01,0.14716482E-01,0.47345284E-02,0.44192205E-03,0.13655324E-04
2,0,2,0.76251665E-02,0.42709723E-02,-0.15275348E-02,-0.53455363E-03,-0.50940514E-04,-0.15623067E-05
2,1,0,-0.31863123E-01,-0.19815411E-01,0.52020996E-02,0.28355253E-02,0.35739999E-03,0.13896490E-04
2,1,1,-0.16244229E-01,0.93012640E-02,0.16637012E-02,-0.68382140E-03,-0.13628576E-03,-0.62450050E-05
2,1,2,0.27835739E-02,-0.97187943E-03,-0.47272017E-03,0.79149520E-05,0.95241782E-05,0.53642837E-06
2,2,0,0.46502480E-02,-0.32396413E-02,-0.15041494E-02,-0.75760531E-04,0.16825701E-04,0.12521577E-05
2,2,1,0.21825930E-03,0.12337052E-02,-0.50023564E-03,-0.24279710E-03,-0.29501141E-04,-0.11178305E-05
2,2,2,0.33644525E-03,0.26880995E-04,0.20519435E-04,0.76226242E-05,0.90872206E-06,0.33814742E-07
2,3,0,0.58192239E-03,0.61367478E-04,-0.17658762E-03,-0.52589587E-04,-0.51889990E-05,-0.16995503E-06
2,3,1,0.26078224E-03,0.54482841E-04,-0.84274690E-04,-0.24125560E-04,-0.22392916E-05,-0.68572929E-07
2,3,2,-0.31587447E-05,0.88101241E-05,0.77940884E-05,0.12632533E-05,0.62626899E-07,0.32419583E-09
3,0,0,0.94294495E-02,0.37121738E-02,-0.96294274E-03,-0.49755041E-03,-0.60494434E-04,-0.23180297E-05
3,0,1,-0.10376127E-01,-0.26461773E-02,0.11205302E-02,0.31769572E-03,0.27403681E-04,0.78509192E-06
3,0,2,0.86310176E-03,0.27169634E-03,-0.12254550E-03,-0.35939188E-04,-0.30618079E-05,-0.83621815E-07
3,1,0,-0.25905794E-02,-0.16267199E-02,0.44903253E-03,0.24949716E-03,0.31811253E-04,0.12472213E-05
3,1,1,-0.93864940E-03,0.78226701E-03,0.65698746E-04,-0.76664572E-04,-0.13261908E-04,-0.58433454E-06
3,1,2,0.20464969E-03,-0.79798598E-04,-0.31340467E-04,0.26263592E-05,0.96852467E-06,0.50155474E-07
3,2,0,0.17937967E-03,-0.26284951E-03,-0.77522141E-04,0.32019542E-05,0.19907745E-05,0.11330622E-06
3,2,1,0.14379785E-03,0.11350385E-03,-0.52325278E-04,-0.22296259E-04,-0.25628334E-05,-0.93609758E-07
3,2,2,0.14439318E-04,0.13040934E-05,0.28224994E-05,0.82772388E-06,0.83029680E-07,0.27175665E-08
3,3,0,0.26685171E-04,0.27702134E-05,-0.91601637E-05,-0.29156552E-05,-0.30693882E-06,-0.10715202E-07
3,3,1,0.30442128E-04,0.55124161E-05,-0.77587641E-05,-0.21593223E-05,-0.19452234E-06,-0.57567222E-08
3,3,2,-0.13835817E-05,0.65511747E-06,0.75405027E-06,0.12656149E-06,0.65449343E-08,0.47775400E-10
~~~

### 20.2 Perpendicular coefficients

The perpendicular table is included because NLGCE-F represents a consistent fitted pair, even when the main library request concerns parallel diffusion.

~~~csv
j,k,l,d_i0,d_i1,d_i2,d_i3,d_i4,d_i5
0,0,0,-0.20782397E+01,0.80542173E+00,0.43509413E-01,-0.21765985E-01,-0.34880893E-02,-0.14567092E-03
0,0,1,-0.46386191E+00,-0.98573919E-01,0.55008025E-02,0.60611018E-02,0.81158994E-03,0.32584425E-04
0,0,2,-0.20027545E-01,0.77767273E-02,-0.17837106E-03,-0.34531413E-03,-0.48546308E-04,-0.20284890E-05
0,1,0,0.32932663E-01,0.55078257E-01,0.57859672E-02,-0.98263229E-02,-0.15658018E-02,-0.64502352E-04
0,1,1,0.79801317E-01,-0.16056809E-01,-0.43059925E-02,0.15297052E-02,0.29837887E-03,0.13184386E-04
0,1,2,-0.96094000E-02,-0.70750559E-03,0.71313130E-04,-0.14638662E-03,-0.27273546E-04,-0.12478895E-05
0,2,0,0.45254433E-01,-0.15503612E-03,-0.58120004E-02,-0.95055489E-03,-0.41356187E-04,0.75481060E-07
0,2,1,-0.12156453E-01,-0.14601653E-02,0.16977702E-02,0.48775342E-03,0.45538505E-04,0.14225158E-05
0,2,2,0.14840409E-02,0.42407064E-03,-0.45749343E-04,-0.31128841E-04,-0.36845885E-05,-0.13458210E-06
0,3,0,0.32974684E-02,-0.60751030E-03,-0.81092420E-03,-0.65981454E-04,0.46104859E-05,0.42766278E-06
0,3,1,-0.10773557E-02,-0.29570722E-03,0.43219766E-04,0.13131432E-05,-0.13053938E-05,-0.84558084E-07
0,3,2,0.12302544E-03,0.56007217E-04,0.10139352E-04,0.24540202E-05,0.33361447E-06,0.14524518E-07
1,0,0,-0.11101067E+01,0.31228772E+00,-0.13177282E-01,-0.22255326E-01,-0.30794494E-02,-0.12445252E-03
1,0,1,-0.27352296E+00,-0.13632322E+00,0.95260785E-02,0.80307990E-02,0.96089118E-03,0.35198094E-04
1,0,2,0.12818726E-01,0.96450258E-02,-0.23937203E-02,-0.11687264E-02,-0.13273061E-03,-0.47474106E-05
1,1,0,-0.10207300E+00,0.11063794E-01,0.55260570E-02,-0.14631952E-02,-0.28667936E-03,-0.12195614E-04
1,1,1,-0.89940234E-02,0.15794357E-01,0.10450265E-01,0.23173513E-02,0.20829148E-03,0.66426911E-05
1,1,2,0.32828567E-02,-0.10346441E-02,-0.10776541E-02,-0.22747122E-03,-0.19402946E-04,-0.60184518E-06
1,2,0,0.16729966E-01,-0.86926828E-02,-0.47449408E-02,-0.25011753E-03,0.52349008E-04,0.39216750E-05
1,2,1,-0.68420124E-02,0.15559548E-02,-0.95408753E-03,-0.59427782E-03,-0.82159092E-04,-0.33850261E-05
1,2,2,0.80345403E-03,-0.54456722E-04,0.18658915E-03,0.84258499E-04,0.10778557E-04,0.42789250E-06
1,3,0,0.22208796E-02,-0.58431586E-03,-0.48924660E-03,-0.27928743E-04,0.51958515E-05,0.39174790E-06
1,3,1,-0.38128951E-04,0.13477913E-03,-0.25264265E-03,-0.10906523E-03,-0.13641539E-04,-0.53623767E-06
1,3,2,0.87033985E-05,-0.14763804E-04,0.25957223E-04,0.11060270E-04,0.13830492E-05,0.54549206E-07
2,0,0,-0.35795231E+00,0.65313508E-01,-0.50976726E-02,-0.52759476E-02,-0.70390662E-03,-0.28026723E-04
2,0,1,-0.80001226E-01,-0.34222505E-01,0.15934024E-02,0.15341728E-02,0.17475526E-03,0.60919963E-05
2,0,2,0.87639890E-02,0.34113860E-02,-0.18169113E-03,-0.16835979E-03,-0.18690148E-04,-0.62524682E-06
2,1,0,-0.26504203E-01,-0.19821713E-02,0.26199589E-02,0.52733957E-03,0.46965191E-04,0.16425008E-05
2,1,1,0.31124751E-03,0.71374069E-02,0.19284097E-02,0.75522804E-04,-0.14569881E-04,-0.96859723E-06
2,1,2,0.55387018E-03,-0.72916289E-03,-0.22397964E-03,-0.51781606E-06,0.30976122E-05,0.17280133E-06
2,2,0,0.25359898E-02,-0.29959443E-02,-0.89093469E-03,0.56975030E-04,0.25922842E-04,0.14212540E-05
2,2,1,-0.77368564E-03,0.92068559E-03,-0.29143432E-03,-0.18699531E-03,-0.25019943E-04,-0.10043456E-05
2,2,2,0.16753641E-03,-0.63376050E-04,0.32826799E-04,0.18838799E-04,0.24833739E-05,0.98856519E-07
2,3,0,0.37522544E-03,-0.17359581E-03,-0.84952992E-04,0.16008488E-05,0.19086792E-05,0.11081771E-06
2,3,1,0.77690104E-04,0.52930691E-04,-0.71571141E-04,-0.28796227E-04,-0.34496071E-05,-0.13183005E-06
2,3,2,-0.33219791E-05,-0.36269908E-05,0.64832962E-05,0.24768957E-05,0.29318259E-06,0.11173399E-07
3,0,0,-0.26370968E-01,0.45853683E-02,-0.42857969E-03,-0.38478482E-03,-0.50492454E-04,-0.19949740E-05
3,0,1,-0.59333101E-02,-0.24985017E-02,0.83613754E-04,0.93963534E-04,0.10360895E-04,0.34793936E-06
3,0,2,0.82864537E-03,0.29126396E-03,0.91218775E-05,-0.58268328E-05,-0.64542109E-06,-0.18840624E-07
3,1,0,-0.21245015E-02,-0.24321040E-03,0.27708187E-03,0.76882914E-04,0.81206123E-05,0.30009780E-06
3,1,1,0.21504225E-03,0.58749970E-03,0.68725126E-04,-0.23819124E-04,-0.45232048E-05,-0.20090275E-06
3,1,2,0.13663304E-04,-0.67148914E-04,-0.80040661E-05,0.38261029E-05,0.69827212E-06,0.30722829E-07
3,2,0,0.86363838E-04,-0.25412269E-03,-0.47306183E-04,0.10707095E-04,0.25867542E-05,0.12762594E-06
3,2,1,0.69686411E-05,0.95111050E-04,-0.24784177E-04,-0.15707547E-04,-0.20542771E-05,-0.81082411E-07
3,2,2,0.86898718E-05,-0.74726285E-05,0.15429839E-05,0.12006276E-05,0.16251525E-06,0.64782353E-08
3,3,0,0.18940167E-04,-0.14961780E-04,-0.47045372E-05,0.55656921E-06,0.17796779E-06,0.91749454E-08
3,3,1,0.10317043E-04,0.54645836E-05,-0.52425736E-05,-0.20926819E-05,-0.24595118E-06,-0.92493584E-08
3,3,2,-0.50524289E-06,-0.31899297E-06,0.39480448E-06,0.14383851E-06,0.16356627E-07,0.60444987E-09
~~~

## 21. BibTeX additions

The following entries are proposed additions for references not found unambiguously in the supplied snapshot. Existing records are not duplicated here, and the conflicting Bieber record is excluded pending correction. Pre-2006 foundational additions are identified in Section 19.

~~~bibtex
@article{Teufel-2003-AA,
  author = {Teufel, A. and Schlickeiser, R.},
  title = {{Analytic calculation of the parallel mean free path of heliospheric cosmic rays. II. Dynamical magnetic slab turbulence and random sweeping slab turbulence with finite wave power at small wavenumbers}},
  journal = {Astronomy \& Astrophysics},
  year = {2003},
  volume = {397},
  pages = {15--25},
  doi = {10.1051/0004-6361:20021471},
  note = {Suggested addition; used as R04 in PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md},
}

@article{Shalchi-2004-AJ-WNLT,
  author = {Shalchi, A. and Bieber, J. W. and Matthaeus, W. H. and Qin, G.},
  title = {{Nonlinear Parallel and Perpendicular Diffusion of Charged Cosmic Rays in Weak Turbulence}},
  journal = {The Astrophysical Journal},
  year = {2004},
  volume = {616},
  pages = {617--629},
  doi = {10.1086/424839},
  note = {Suggested addition; used as R06 in PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md},
}

@article{Qin-2007-AJ,
  author = {Qin, G.},
  title = {{Nonlinear Parallel Diffusion of Charged Particles: Extension to the Nonlinear Guiding Center Theory}},
  journal = {The Astrophysical Journal},
  year = {2007},
  volume = {656},
  pages = {217--221},
  doi = {10.1086/510510},
  note = {Suggested addition; used as R08 in PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md},
}

@article{Qin-2013-AJ-Erratum,
  author = {Qin, G.},
  title = {{Erratum: Nonlinear Parallel Diffusion of Charged Particles: Extension to the Nonlinear Guiding Center Theory (2007, ApJ, 656, 217)}},
  journal = {The Astrophysical Journal},
  year = {2013},
  volume = {774},
  pages = {91},
  doi = {10.1088/0004-637X/774/1/91},
  note = {Suggested addition; used as R09 in PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md},
}

@article{Qin-2014-AJ-NLGCE,
  author = {Qin, G. and Zhang, L.-H.},
  title = {{The Modification of the Nonlinear Guiding Center Theory}},
  journal = {The Astrophysical Journal},
  year = {2014},
  volume = {787},
  pages = {12},
  doi = {10.1088/0004-637X/787/1/12},
  url = {https://arxiv.org/abs/1401.1950},
  note = {Suggested addition; used as R10 in PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md},
}

@article{Engelbrecht-2013-AJ-772,
  author = {Engelbrecht, N. E. and Burger, R. A.},
  title = {{An Ab Initio Model for Cosmic-Ray Modulation}},
  journal = {The Astrophysical Journal},
  year = {2013},
  volume = {772},
  pages = {46},
  doi = {10.1088/0004-637X/772/1/46},
  note = {Suggested addition; used as R12 in PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md},
}

@article{Engelbrecht-2013-AJ-779,
  author = {Engelbrecht, N. E. and Burger, R. A.},
  title = {{An Ab Initio Model for the Modulation of Galactic Cosmic-Ray Electrons}},
  journal = {The Astrophysical Journal},
  year = {2013},
  volume = {779},
  pages = {158},
  doi = {10.1088/0004-637X/779/2/158},
  note = {Suggested addition; used as R13 in PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md},
}

@article{Potgieter-2014-SP-391,
  author = {Potgieter, M. S. and Vos, E. E. and Boezio, M. and De Simone, N. and Di Felice, V. and Formato, V.},
  title = {{Modulation of Galactic Protons in the Heliosphere During the Unusual Solar Minimum of 2006 to 2009}},
  journal = {Solar Physics},
  year = {2014},
  volume = {289},
  pages = {391--406},
  doi = {10.1007/s11207-013-0324-6},
  url = {https://arxiv.org/abs/1302.1284},
  note = {Suggested addition; used as R14 in PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md},
}

@article{Vos-2015-AJ,
  author = {Vos, E. E. and Potgieter, M. S.},
  title = {{New Modeling of Galactic Proton Modulation During the Minimum of Solar Cycle 23/24}},
  journal = {The Astrophysical Journal},
  year = {2015},
  volume = {815},
  pages = {119},
  doi = {10.1088/0004-637X/815/2/119},
  note = {Suggested addition; used as R15 in PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md},
}

@article{Chhiber-2017-AJSS,
  author = {Chhiber, R. and Subedi, P. and Usmanov, A. V. and Matthaeus, W. H. and Ruffolo, D. and Goldstein, M. L. and Parashar, T. N.},
  title = {{Cosmic-Ray Diffusion Coefficients throughout the Inner Heliosphere from a Global Solar Wind Simulation}},
  journal = {The Astrophysical Journal Supplement Series},
  year = {2017},
  volume = {230},
  pages = {21},
  doi = {10.3847/1538-4365/aa74d2},
  url = {https://arxiv.org/abs/1703.10322},
  note = {Suggested addition; used as R17 in PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md},
}

@article{Corti-2019-AJ,
  author = {Corti, C. and Potgieter, M. S. and Bindi, V. and Consolandi, C. and Light, C. and Palermo, M. and Popkow, A.},
  title = {{Numerical Modeling of Galactic Cosmic-Ray Proton and Helium Observed by AMS-02 during the Solar Maximum of Solar Cycle 24}},
  journal = {The Astrophysical Journal},
  year = {2019},
  volume = {871},
  pages = {253},
  doi = {10.3847/1538-4357/aafac4},
  url = {https://arxiv.org/abs/1810.09640},
  note = {Suggested addition; used as R18 in PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md},
}

@article{Alonso-Guzman-2025-JGR,
  author = {{Alonso Guzm{\'a}n}, J. G. and Ghanbari, K. and Florinski, V. A. and Leske, R. A. and Zhao, L.-L. and Zhu, X. and others},
  title = {{Superposed Epoch Analysis of Stream Interaction Regions at 1 au During Solar Minimum With Turbulence Geometry Decomposition: Implications for Galactic Cosmic Ray Transport}},
  journal = {Journal of Geophysical Research: Space Physics},
  year = {2025},
  volume = {130},
  pages = {e2024JA033567},
  doi = {10.1029/2024JA033567},
  note = {Suggested addition; used as R28 in PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md},
}

@article{Hussein-2014-AJ-Bohm,
  author = {Hussein, M. and Shalchi, A.},
  title = {{Detailed Numerical Investigation of the Bohm Limit in Cosmic Ray Diffusion Theory}},
  journal = {The Astrophysical Journal},
  year = {2014},
  volume = {785},
  pages = {31},
  doi = {10.1088/0004-637X/785/1/31},
  note = {Suggested addition; used as R29 in PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md},
}
~~~



### 21.1 Additional citation for the reviewed revision

The following paper was not found in the same supplied bibliography snapshot. R31 already has the verified author-preprint key Marcowith-1999-arXiv and is not duplicated here.

~~~bibtex
@article{Droege-2014-JGR,
  author = {Dr{\"o}ge, W. and Kartavykh, Y. Y. and Dresing, N. and Heber, B. and Klassen, A.},
  title = {Wide longitudinal distribution of interplanetary electrons following the 7 February 2010 solar event: Observations and transport modeling},
  journal = {Journal of Geophysical Research: Space Physics},
  year = {2014},
  volume = {119},
  pages = {6074--6094},
  doi = {10.1002/2014JA019933},
  url = {https://doi.org/10.1002/2014JA019933}
}
~~~

## 22. Reproducibility assets and revision record

### 22.1 Companion data contract

The downloadable PARALLEL_DIFFUSION_COEFFICIENT_MODEL_DATA.zip extracts to parallel_diffusion_model_data/. This is a numerical reference bundle for a future implementation; it is not the diffusion-coefficient library itself.

| File | Contents and required use |
|---|---|
| NLGCE_F_2014_parallel.csv | 48 rows, 288 coefficients; j,k,l followed by i=0,...,5 |
| NLGCE_F_2014_perpendicular.csv | Same indexing for the perpendicular polynomial |
| benchmark_points.json | Species/unit, pitch-angle, QLT, and A–G nonlinear/polynomial reference values |
| fit_error_samples.csv | All 300 input states, both backend pairs, errors, a_x, and nonlinear residuals |
| fit_error_summary.json | Generator, spectrum/length conventions, solver settings, environment, and binned statistics |
| configuration_examples.json | Parameterized, non-normative templates from Section 14.6; no invented physical defaults |
| reference_verification.py | Standalone numerical verification and optional audit reproduction |
| README.md | Data provenance, commands, dependencies, conventions, and interpretation |
| SHA256SUMS | Digests of the other companion files |

CSV coefficient rows have the exact header j,k,l,d_i0,d_i1,d_i2,d_i3,d_i4,d_i5. Construct d[i,j,k,l] explicitly; do not infer a Fortran memory layout from the row order. There are no duplicate or missing tuples in 0<=j<=3, 0<=k<=3, 0<=l<=2.

The coefficient checksums are:

~~~text
7cdc5ab9cddda0c7295a25bfcaff4ba3422002da12ac1b4660a98dcc64eaa762  NLGCE_F_2014_parallel.csv
7371b6557f40c5efa7a48460a5f1b2971374a2c43d0eb3e209f0a8e62670356a  NLGCE_F_2014_perpendicular.csv
~~~

Verify the digests before importing fixtures. A checksum protects the delivered transcription from later changes; source-level comparison, indexed parsing, and backend-specific numerical checks are separate verification steps.

From the extracted directory, run:

~~~bash
python3 reference_verification.py
python3 reference_verification.py --audit
~~~

The first command verifies checksums and deterministic reference values. The second also evaluates all supplied audit inputs, compares their recorded results, and prints the six-bin summary. Dependencies are NumPy and SciPy; the original evaluation environment and tolerances are recorded in fit_error_summary.json. Different compatible library versions may change the last floating-point digits.

### 22.2 Earlier review revisions and their correction

The review revisions clarified slab-only L_c; numerical and physical validity diagnostics; symbol and Elsässer conventions; field-direction derivatives; batching and optional inputs; publication reporting; extension scope; and bibliography matching. The 576 polynomial coefficients were correctly transcribed.

However, the earlier Equation (45) incorrectly used exponentiation where the source uses multiplication. Consequently the former nonlinear tables, weak-turbulence discrepancy claims and turbulence-only warning rule are superseded by revision 1.2. Equation (45), all affected tables and the companion calculation now use the published multiplicative form. Original equation numbers (1)–(63) are retained; supplementary identities use letter suffixes.

### 22.3 Numerical-value audit and revision 1.2

Every quantitative value must be identifiable as a published value, a mathematical constant/conversion, a calculated result at explicitly defined inputs, or a proposed numerical control. A numerical test input is not an observed turbulence state or a calibrated parameter. These categories must remain explicit when adapting this document for a publication or implementation.

| Numerical content | Provenance and check |
|---|---|
| NLGCE-F fitting bounds, Section 11.2 | Qin and Zhang (2014), Table 2 [R10]; all four original ratios and both bounds retained |
| Coefficient tables, Section 20 | Qin and Zhang (2014), Tables 3 and 4 [R10]; all 576 decimal values individually compared with a fresh retrieval of the source tables, with no discrepancy |
| Bibliographic numbers in Sections 18 and 21 | All 30 DOI-backed records checked against publisher-deposited Crossref metadata for year, volume and indexed page/article number; R31 is explicitly cited as its author preprint |
| AU and SI conversions, Section 2.2 | IAU Resolution B2 (2012) and SI definitions/prefixes; the AU value is exact |
| Proton example, Section 15.2 | Fixed CODATA 2018 proton mass and exact c, e, AU; relativistic arithmetic independently checked with sixty-digit decimal arithmetic |
| Spectrum, pitch-angle and QLT values, Section 15.3 and Equation (32) | Calculated from the displayed equations at defined inputs; quadrature, gamma-function normalization and the closed-form pitch-angle check agree |
| A–G values and signed discrepancies, Section 15.4 | Original verification calculations, not values claimed to have been published or measured; both nonlinear integration methods and independent high-precision polynomial evaluation reproduce them |
| Counts, medians and discrepancy shares, Section 11.4 | Recomputed from all 300 retained input/output rows; the draw, seed and weighting are specified, and the complete dataset is supplied |
| Refinement and alternate-seed statements | Measured comparisons recorded in the companion JSON files; these statements describe numerical convergence, not physical accuracy |
| Configuration templates, Section 14.6 | Physical parameter values are placeholders requiring caller input or cited calibration |
| Tolerances, numerical intervals and diagnostic thresholds | Proposed or recorded computational controls; none is presented as a measured universal physical boundary or certified surrogate-error bound |

Computed values in printed benchmark tables are rounded to twelve significant digits; percentages retain the stated rounding. Machine-readable values are numerical fixtures with explicit comparison tolerances, not assertions that every stored floating-point digit has physical significance. Numerical convergence and arithmetic checks verify the stated model evaluation; approximation quality and physical validity require their separate assessments.

Revision 1.2 corrects the a_x transcription, regenerates all affected nonlinear benchmarks and error statistics, removes the unsupported turbulence-only diagnostic, and adds explicit constant sources, independent numerical checks and parameterized configuration templates. It also corrects the companion README's erratum DOI to 10.1088/0004-637X/774/1/91; the R09 reference and BibTeX in the main document already used that correct DOI. Published coefficients are unchanged; the nonlinear parameter is restored to its published expression rather than calibrated.

## 23. Codex implementation roadmap

### 23.1 Implementation contract and scope

This roadmap turns the mathematical specification into independently reviewable implementation stages. It is a development plan for Codex, not a statement that the library already exists. Sections 1–15 determine the model equations, units, domains, output contracts, and numerical checks; the steps below determine how to build and verify them. Revision 1.3 adds this roadmap without changing the equations, numerical tables, or companion data audited in revision 1.2.

The first release implements the fully specified models in Section 3, together with the explicitly defined provider adapters. Its reusable core evaluates local coefficients and scattering functions. Background evolution, particle transport, mesh ownership, shock detection, wave growth, and MPI communication remain responsibilities of the application. The standalone core must be usable without an MHD simulation or a particular transport solver.

Before choosing source paths, dependencies, or a build system, Codex must inspect the actual target repository and its applicable AGENTS.md instructions. The C++20 sketch in Section 14.1 is a starting point for the data contract; adapt its syntax to the repository's supported language standard. Do not introduce a compiler requirement merely to reproduce the sketch. The file paths below are proposed organizational names, to be mapped to the actual repository during PD00.

Keep the published integral closures, the published polynomial, and phenomenological prescriptions identifiable throughout implementation. A numerical failure is not a reason to change a published parameter, insert a scattering floor, or return another backend under the requested identifier. Physical parameter values must be supplied by the caller or by a separately documented calibration. The configuration templates in Section 14.6 deliberately leave those values unspecified.

Complete SOQLT, complete composite WNLT, arbitrary directional/propagating-wave scattering, D_pp, and non-axisymmetric perpendicular dynamics are later extensions as described in Section 3.1. The first release may define extension interfaces and reject unsupported requests, but it must not advertise these theories as implemented. First-order shock acceleration through an existing transport solver does not by itself require a D_pp module; see Section 1.3.

### 23.2 Proposed module organization

Use the repository's existing conventions for public headers, implementation files, tests, and installation. The following layout expresses module responsibilities, not a claim that these paths currently exist.

| Proposed component | Responsibility and dependency boundary |
|---|---|
| include/parallel_diffusion/ | Public particle, local-state, configuration, result, provider, and evaluator interfaces; document SI units and optional quantities |
| src/parallel_diffusion/core/ | Central relativistic conversions, validation, status/diagnostic handling, configuration validation, provenance, and model registry |
| src/parallel_diffusion/models/ | Explicit laws, QLT, prescribed pitch-angle shape, slab broadening, NLPA/NLGC pairs, NLGCE-F, and coefficient tables |
| src/parallel_diffusion/numerics/ | Adaptive quadrature, positive/logarithmic nonlinear solves, convergence reporting, and derivative utilities shared by models |
| src/parallel_diffusion/providers/ | Spectrum and turbulence-moment conversion; immutable provider snapshots; no background or wave evolution |
| data/parallel_diffusion/ | Exact coefficient CSV files, backend fixtures, audit inputs, checksums, and provenance imported from Section 22 |
| tests/parallel_diffusion/ | Standalone model/reference tests and an aggregate runner with machine-readable results |
| adapters/parallel_diffusion/ | Application-specific background and transport bindings; dependency direction is from adapter to core |
| examples/parallel_diffusion/ | Small standalone evaluation programs and parameterized configuration templates |
| README.md and NUMERICS.md | Build/use instructions, implemented model inventory, normalization, numerical controls, and limitations |
| IMPLEMENTATION_STATUS.md | Stage state, acceptance evidence, unresolved work, and a continuation handoff for Codex |

If the host project uses a shared src/models hierarchy, place the reusable implementation there and map the public headers and tests accordingly. Existing SEP or GCR applications should call this shared implementation through adapters. Avoid independent copies of the same equations or coefficient tables in different applications.

The public configuration must be model-specific: required parameters, permitted provider types, independent-variable conventions, derivative requests, boundary policies, and fallback choices must be validated before evaluation. A model registry should expose which outputs and inputs each backend supports. An optional D_mu_mu evaluator must have a typed interface and lifetime rules; an unavailable scattering function must not be represented by a dummy successful evaluator.

### 23.3 Stage order and dependencies

Each stage ends with an implemented artifact, focused verification, and an update to IMPLEMENTATION_STATUS.md. Dependencies below are acceptance prerequisites. Where a module would benefit from an optional later facility, implement the simpler valid path first and record the later connection explicitly.

| Stage | Deliverable | Required earlier stages | Acceptance evidence |
|---|---|---|---|
| PD00 | Repository map, scope manifest, and verified reference-data import | None | Actual build conventions recorded; companion digests and indexed coefficient layout checked |
| PD01 | Shared API, SI conversion, validation, registry, and provenance | PD00 | Particle/unit fixtures and missing-input/status behavior pass |
| PD02 | Constant, power-law, broken-rigidity, and Bohm models | PD01 | Normalization, slopes, limits, and selective field requirements pass |
| PD03 | Canonical spectra, provider contracts, and shared quadrature | PD01 | Integrated variances, sidedness/component conversions, and refinement checks pass |
| PD04 | Pitch-angle conversion, regularized shape, and slab QLT | PD02, PD03 | Exact normalization, supplied QLT fixtures, endpoint handling, and divergence detection pass |
| PD05 | Exact published NLGCE-F polynomial and derivatives | PD01, PD03 | Coefficient digests, backend fixtures, input bounds, and polynomial derivatives pass |
| PD06 | NLPA, original NLGC-E, and modified NLGCE-N integral solves | PD03, PD04 | Correct a_x fixtures, positive roots, residuals, quadrature refinement, and backend-specific checks pass |
| PD07 | Lorentzian and Gaussian slab resonance broadening | PD03, PD04 | Kernel normalization, resolved resonance integration, and regular zero-width QLT limit pass |
| PD08 | Turbulence/wave adapters and supplied coefficient tables | PD02, PD03, PD04, PD05, PD06, PD07 | Convention conversions, closure dispatch, interpolation policies, and provenance pass |
| PD09 | Coefficient gradients and transport derivative contract | PD02, PD04, PD05, PD06, PD07, PD08 | Analytical/independent derivative checks and missing-derivative behavior pass |
| PD10 | Batch evaluation, immutable state, and correct caches | PD09 | Scalar/batch agreement, per-point failures, concurrency, and state-change invalidation pass |
| PD11 | Host application bindings and transport verification | PD10 | Real consumer uses shared core; applicable transport checks in Section 15.5 pass |
| PD12 | Release documentation, installation, and final evidence | PD11 | Supported inventory, reproducible tests, examples, and completion checklist agree with the delivered code |

PD05 does not depend on a nonlinear solver: the polynomial is fully specified by its own data. PD06 may subsequently use an accepted PD05 value as an initial guess, but its converged result must be accepted by integral residuals. If PD05 is unavailable or its input box is exceeded, PD06 must still have an independently usable positive initialization strategy. The stage numbers are task identifiers, not time estimates or performance targets.

### 23.4 Detailed implementation stages

#### PD00 — Inspect the repository and freeze the reference contract

- Read the applicable repository instructions and identify the supported compiler/language standard, build targets, dependency policy, test runner, installation rules, and application entry points. Record the chosen core and adapter paths before creating files.
- Map every stable identifier in Section 3 to its equations, required inputs, outputs, and proposed source/test module. Mark paper-specific extensions from Section 3.1 as deferred. There must be no unexplained identifier whose implementation silently delegates to a different physical model.
- Import the Section 22 data bundle without rewriting its decimal coefficient literals. Verify SHA256SUMS and the coefficient digests listed in Section 22.1. Parse j,k,l and the six i columns explicitly, checking unique indices and complete array coverage.
- Record this specification revision and the reference bundle's own metadata/digests separately. A data schema version is not the same quantity as a document revision or software release version.
- Establish a standalone verification target and the stage/test manifest. The supplied Python reference verifier is independent numerical evidence and may support cross-checks; it is not the C++ library or a production runtime dependency by default.

**Acceptance gate:** the repository map and model manifest are concrete, the delivered data are intact, and the reference checks can be run with their documented dependencies. No physical default, fit coefficient, reference value, or performance claim has been invented.

#### PD01 — Implement shared state, units, results, and model selection

- Implement the particle and local-state contracts of Section 14.1 using the host project's type conventions. Represent missing field, turbulence, perpendicular input, derivative, and failed output values explicitly.
- Centralize relativistic particle calculations and unit conversion from Sections 2.1–2.2. Retain signed charge where relevant while using its magnitude in the specified balanced scalar closures. An energy-per-nucleon input requires sufficient species information; do not infer a mass number from a generic charge or an unrelated species field.
- Validate only inputs required by the selected model. Missing B0 must not prevent evaluation of constant_lambda or constant_kappa. An invalid charged-particle state remains an error even if an explicit coefficient was supplied.
- Separate numerical status from physical validity diagnostics. Define the semantics of non-finite inputs, overflow, absent scattering, missing requested outputs, unknown model/configuration options, and explicitly selected fallback or bounds in the actual API.
- Implement reproducible configuration identity and provenance: requested/evaluated model, specification/software identity, coefficient-set digest where relevant, original provider convention, and state revisions. Preserve enough information to report raw and bounded outputs when a caller explicitly selects a bound.

**Acceptance gate:** the units/species example in Section 15.2 agrees with the full-precision fixture at its declared tolerance, and negative tests exercise the status distinctions in Section 14.5. Successful coefficient evaluation must not imply that missing gradients or unavailable outputs were computed.

#### PD02 — Implement explicit coefficient prescriptions

- Implement constant_lambda and constant_kappa as different models using Equations (11),(12). Their configuration schemas and speed dependence must remain distinct.
- Implement power_law_lambda with an explicit independent-variable option for rigidity, total kinetic energy, energy per nucleon, or speed. Apply the separately defined variants of Section 5, including enabled radial, field, time, and region factors. Disabled factors equal unity and introduce no additional input requirement.
- Implement broken_rigidity_kappa with the speed-factored K_star convention of Equations (17),(18), stable logarithmic evaluation, and an explicit converter if the caller instead supplies the actual reference coefficient. Do not silently treat K_star as kappa at the reference rigidity.
- Implement bohm from Equation (61) with a required field definition and explicitly configured eta_B. A total/effective-field variant must identify that convention; the Bohm value must not automatically bound other backends.
- Add model comments identifying the defining equations and units, and expose the analytical parameter/rigidity derivatives needed for later gradient support.

**Acceptance gate:** reference-state normalization, Equation (4), the a=b broken-law limit, the beta contribution to kappa slopes, and equal-rigidity species behavior all pass. Field-independent configurations operate without a field provider. Example physical parameters remain caller-supplied; unit tests may use explicitly labeled mathematical inputs.

#### PD03 — Implement spectrum conventions and shared numerical infrastructure

- Define a canonical one-sided total-transverse magnetic spectrum satisfying Equation (25). Spectrum objects must expose units, variance, wavenumber convention, valid coverage, tail policy, spectral breaks, and source identity.
- Implement the normalized smooth bend-over spectrum and the explicit multirange option in Equations (27),(35), including the s=1 middle-range limit. Component/two-sided spectra must be converted through named adapters before they reach a closure.
- Implement supplied-spectrum interpolation and coverage queries without confusing an unmeasured spectral interval with a physically zero spectrum. Frozen-flow frequency-to-wavenumber conversion is allowed only when its sampling velocity/projection and assumptions are supplied as in Equation (36).
- Implement reusable adaptive integration with error estimates, endpoint transformations, spectral-break partitioning, and refinement reporting. Quadrature ranges and transformations are numerical controls; a finite working integration interval must not be reinterpreted as a physical spectrum cutoff.
- Specify immutable provider snapshots and handle lifetimes so that one evaluation sees one coherent field/turbulence/spectrum state. These contracts will also be used by batch evaluation and cache invalidation.

**Acceptance gate:** analytic and numerical integrals reproduce the declared variance; one-/two-sided and component conversions give the same canonical spectrum; refinement resolves representative endpoint and tail integrals. Incorrect normalization, insufficient coverage, and invalid spectrum parameters produce explicit failures.

#### PD04 — Implement pitch-angle normalization and slab QLT

- Implement Equation (20) as a shared conversion routine with explicitly declared symmetry. For even magnetostatic models use the reduced interval; preserve the full interval for a future closure whose D_mu_mu is not even.
- Implement prescribed_lambda_mu_shape with separate target-lambda and fixed-amplitude modes. Use Equations (22),(23), reduced endpoint-safe integrands, and the Equation (24) check where applicable. Recompute D0 in target mode when a shape parameter changes; do not leave a stale amplitude.
- Implement qlt_slab_spectrum using Equation (26). Provide the smooth-spectrum Equation (30) route with the Equation (37) transformation and a general supplied-spectrum route with break-aware quadrature.
- Implement qlt_slab_inertial as the explicit Equation (31) approximation. Keep its identity and assumptions distinct from the full smooth spectrum, and compute its prefactor from the defining expression rather than the rounded display value.
- Determine convergence from the actual high-k behavior. Test hard cutoffs, dissipation tails, and genuine zero resonant power. An arbitrary mu exclusion interval or numerical scattering floor must not convert a divergent QLT model into successful finite diffusion.

**Acceptance gate:** the isotropic normalization, regularized-shape fixture, exact QLT fixtures, and appropriate inertial/high-rigidity limits pass. Missing spectral coverage is distinguishable from infinite_mean_free_path. The focused-transport consumer can obtain the same declared D_mu_mu used to compute the Parker coefficient.

#### PD05 — Implement the exact NLGCE-F polynomial

- Implement the nlgce_f_2014 stable identifier as the published 2014 backend, with its own validated configuration and result identity.
- Import both published coefficient arrays with explicit d[i,j,k,l] indexing and exact source decimal literals. Check their digests/coverage during data generation or build verification; do not hand-edit coefficients to improve agreement with another backend.
- Compute the original dimensionless ratios in Equation (52), validate all published bounds before logarithms, and evaluate Equation (53) with natural logs and the nested Horner convention in Section 11.3. Return both eigenvalues and Equation (4) conversions together.
- Differentiate the polynomial analytically with respect to all four log inputs. Provide the chain-rule quantities required by Equations (54)–(56); preserve the distinction between total epsilon squared in x3 and epsilon in the nonlinear a_x calculation.
- Implement explicit out-of-domain, overflow, missing-input, and caller-selected fallback behavior. In-domain evaluation returns the published polynomial and surrogate_error_unbounded; an input-box check is not a certified local error estimate.
- Keep coefficient-set identity in the result and report both requested and evaluated backends if a caller activates a different fallback. Do not introduce an automatic weak-turbulence threshold based on superseded review calculations.

**Acceptance gate:** both fitted outputs at A–G agree with their own full-precision fixtures. Derivatives agree with independent interior checks, each original ratio is checked at both boundaries and outside them, and a deliberate indexing/logarithm substitution is caught by the tests. Agreement with the nonlinear integral solution is reported separately from polynomial reproduction.

#### PD06 — Implement NLPA and the coupled integral closures

- Implement the explicit magnetostatic slab/2D reduction in Section 10. Use total magnetic variance for epsilon and the slab-only correlation length from Equation (28). Compute the Equation (45) denominator as (xi/(1+xi))/epsilon + epsilon/(2xi), with xi from Equation (45a).
- Implement nlpa_given_perp by taking an explicitly supplied positive perpendicular coefficient and solving Equation (48). Preserve the supplied perpendicular model/value/revision in provenance; if its derivatives are later requested, require that provider to supply them.
- Implement nlgc_e with Equations (48),(49), and nlgce_n with Equations (48),(50),(51), as distinct backend selections. The modified perpendicular equation contains only its specified direct 2D integral. Keep a_x, a squared, and a-prime squared separate in names and code.
- Use dimensionless spectral integrals and positive logarithmic unknowns as required by Section 10.5. Solve coupled residuals simultaneously, monitor quadrature error, and return iteration/residual information. A last iterate or a polynomial initial guess cannot count as a converged integral evaluation.
- Validate representative roots with alternate positive initial guesses and quadrature refinement. Use continuation or an in-domain accepted NLGCE-F value as an optional initial guess. Keep the independent initialization route available and report branch dependence if it is observed.
- Restrict published-fit reproduction to its specified spectrum and index. A generalized spectral closure requires a distinct documented configuration/model identity and its own validation. Pure-component endpoint evaluators remain unsupported until their limits have been explicitly developed.

**Acceptance gate:** nlgce_n reproduces both nonlinear outputs at A–G, satisfies the residual criterion in Section 10.5, and survives the stated refinement/initialization checks. B and D–G explicitly test the multiplicative a_x at epsilon different from unity. For nlpa_given_perp and nlgc_e, create independent residual/refinement checks and numerical references derived from their own equations; the NLGCE-N pair is not a reference for those different models. A weak-turbulence QLT limit is not an acceptance requirement for this NLPA closure.

#### PD07 — Implement the specified slab broadening kernels

- Implement broadened_slab from Equation (39), retaining both resonances and the canonical spectrum amplitude. Support the positive-width Lorentzian options of Equations (40),(41) and the explicitly defined Gaussian width of Equation (42).
- Resolve narrow resonance peaks adaptively in wavenumber and estimate the nested integral error before accepting the pitch-angle conversion. A fixed grid that misses a narrow peak cannot pass merely because its computed integral is small.
- Define an exact zero-width QLT branch. Use any small-width approximation only under a stated numerical acceptance check and where the QLT comparison is regular; finite broadening must not be discarded in a regime where it changes convergence.
- Treat decorrelation widths/speeds as supplied model inputs. Unsupported 2D/3D contributions, directional wave physics, or a request for complete WNLT/SOQLT must be rejected or identified as an unimplemented extension.

**Acceptance gate:** each frequency kernel has area pi under its stated width definition, the balanced regular zero-width comparison approaches QLT, and resonance refinement is demonstrated. These checks establish the declared slab model; the output/model documentation must retain that scope.

#### PD08 — Implement provider adapters and coefficient tables

- Implement turbulence_adapter as moment conversion plus an explicitly selected accepted closure. Convert the named energy and residual-energy conventions in Equations (57)–(58) once, retain the original convention, and reject unknown conventions or inadmissible density/variance inputs.
- Test equivalent half-sum, unhalved Elsasser, and specific-energy representations against the same canonical magnetic variance. A supplied total moment does not determine directional spectra, cross helicity, slab fraction, or lengths; those required quantities need separate inputs.
- Implement wave_spectrum_adapter with direction, polarization, frame, normalization, update time, and provider revision. The balanced magnetostatic reduction may dispatch to accepted QLT/broadened slab closures when justified. Other directional/propagating requests require a separately specified compatible scattering closure; retaining a provider handle does not implement that missing physics.
- Implement tabulated_parallel with an explicit stored quantity, axes, units, and boundary policy. Follow Section 13.2 for log interpolation of positive spatial/particle axes and coefficients, and for declared linear/step time handling. Convert a stored lambda to kappa with the current particle speed or vice versa.
- Preserve original table-generation/spectrum provenance, identify interpolation knots and step discontinuities, and avoid automatic extrapolation or undocumented cross-closure calibration.

**Acceptance gate:** equivalent moment conventions give identical coefficients after conversion; provider updates reach the selected backend; malformed/underspecified waves fail explicitly; table knots, interpolation, time rules, and boundary failures are reproducible. Every adapter result identifies its underlying evaluated closure.

#### PD09 — Implement requested gradients and transport contracts

- Supply analytical spatial/rigidity derivatives for explicit prescriptions and NLGCE-F, chaining through all enabled field, turbulence, length, region, and time providers as appropriate. Spatial gradients are evaluated at fixed particle momentum unless the API explicitly states another derivative.
- For integral closures, use a validated implicit derivative or provider-consistent differentiation of the full converged evaluation. Account for the externally supplied perpendicular coefficient in nlpa_given_perp. Report missing/failed derivatives separately from a successfully computed scalar value.
- Implement table derivatives only where the interpolation rule defines them. Return a documented one-sided convention or an unavailable/discontinuous indication at knots, step boundaries, region changes, and active bounds; do not silently report a smooth derivative there.
- Require the field-direction Jacobian and perpendicular-coefficient gradient needed by Equation (62a). The core supplies gradients of coefficients it determines; the transport adapter/consumer assembles and projects the tensor and its divergence using its own geometry.
- Keep radial and shock-normal projections distinct from a parallel eigenvalue, and retain the area derivative in the flux-tube operator of Equation (63). A non-axisymmetric perpendicular tensor additionally needs its basis derivatives and remains an explicit consumer capability.

**Acceptance gate:** smooth analytic cases agree with independent derivative checks under refinement, and the field-Jacobian convention reproduces Equations (62a),(62b). Missing gradients remain absent; a genuinely uniform provider can supply explicit zeros. A request for a complete derivative-dependent operator fails if a required derivative is missing.

#### PD10 — Implement batch evaluation and state-safe reuse

- Implement the equal-length scalar/batch contract of Section 14.1, including an empty batch, shape validation before execution, and per-point status/optional outputs. A convenience broadcast operation must be separately named and explicit.
- Use immutable provider snapshots or properly synchronized caches. Each evaluation must use a coherent state; numerical differentiation must not combine different background revisions.
- Key reused results by species, momentum, model/configuration, coefficient set, derivative request, provider identity/revision, and every locally relevant sampled input. Include position/time or an equivalent complete local-state identity whenever they affect the result. A single global background revision is insufficient to distinguish different positions in the same background.
- Include spectrum, kernel, turbulence, and supplied-perpendicular revisions where used. If a provider cannot supply a reliable revision/state identity, disable reuse for its changing inputs rather than assume they are constant.
- Optimize only after reference and scalar/batch equivalence pass. Vectorized Horner evaluation, grouped quadrature, and accepted continuation seeds may improve cost, but do not change the requested backend, convergence gate, or failure semantics. Keep MPI operations in the host adapter, outside the standalone core.

**Acceptance gate:** scalar and batch outputs agree at their declared tolerances; an invalid point leaves unrelated results intact; repeated concurrent calls are deterministic within the declared numerical tolerance. Changes in field, turbulence, position, time, species, momentum, spectrum, or supplied perpendicular state cannot return a stale cached result. Record measured optimization evidence without inventing a throughput target.

#### PD11 — Bind the shared core to actual transport consumers

- Implement the host background adapter using the repository paths selected in PD00. Convert host units and moment conventions at that boundary, sample a coherent local state, and pass application revisions and available gradients.
- Connect the existing Parker consumer to kappa and its required tensor derivatives. Connect any focused-transport consumer to its D_mu_mu interface when the selected model supplies one. An eigenvalue-only model needs a separately declared pitch-angle prescription before a focused solver can use it.
- Prevent duplicate parallel diffusion from representing the same scattering both through pitch-angle dynamics and an additional spatial diffusion operator. Keep injection, convection, drift, shock geometry, and wave growth in their existing application modules.
- For a host with multiple SEP/GCR applications, verify that they call the same shared coefficient implementation and data. Per-application configuration differences must remain explicit; do not resolve them by copying or silently altering a backend.
- Run the applicable transport checks in Section 15.5: homogeneous diffusion, a manufactured varying-coefficient solution, tensor projection in a Parker spiral, focused-to-Parker comparison where supported, and oblique-shock normal diffusion where supported. Define convergence/error criteria from the solver's own discretization and sampling; these tests are separate from coefficient arithmetic checks.

**Acceptance gate:** a real consumer builds and exercises the shared core, and the applicable transport tests demonstrate correct coefficient/operator coupling. If the implementation checkout contains no host solver, deliver the standalone core and adapter contract with PD11 explicitly blocked pending a named host checkout; do not claim transport integration from a coefficient-only example.

#### PD12 — Document, package, and close the implementation

- Install/export the public interface and required data using the host project's build conventions. Provide a standalone evaluation example, an aggregate verification command, and application examples only for bindings actually implemented.
- Document supported model IDs, units, normalization, required inputs, numerical controls, derivative availability, fallback/bound semantics, provider lifetime, and thread-safety. Comment each scientific implementation with its defining equation and reference; explain non-obvious endpoint transformations and convergence decisions.
- Preserve the Section 22 fixture provenance and add provenance for newly calculated references. Label test states as mathematical inputs. A calibrated observational preset requires its source, particle population, domain, and physical interpretation; do not turn placeholder templates into unexplained defaults.
- Verify that documentation, model registry, available outputs, installation, and tests agree. Deferred physics must remain listed as deferred. Remove temporary debug paths and ensure examples do not depend on a developer's private absolute paths.
- Complete the final checklist in Section 23.7 and update the continuation record with the actual delivered revision, commands, outcomes, and outstanding limitations.

**Acceptance gate:** a clean supported build can use the installed library and run the documented verification command; the supported test manifest passes; installation preserves required data identity. Release notes report implemented behavior and actual checks, with unavailable host/extensions explicitly identified.

### 23.5 Numerical acceptance and test reporting

Import the existing machine-readable references; do not transcribe rounded publication tables back into tests. The following controls already occur in the specification or companion fixtures and are restated here for implementation planning. They are numerical controls, not measured physical quantities or guarantees of model accuracy.

| Check | Required comparison/acceptance |
|---|---|
| Polynomial outputs against their fixtures | Relative tolerance 10^-11 from benchmark_points.json |
| Modified nonlinear outputs against their fixtures | Relative tolerance 5×10^-8 from benchmark_points.json |
| Supplied units/constants fixtures | Relative tolerance 10^-12 from benchmark_points.json, with the same explicitly fixed constants |
| Coupled nonlinear convergence | Maximum absolute logarithmic residual below 10^-8, plus quadrature convergence; Section 10.5 |
| General pitch-angle quadrature | Initial recommended relative tolerance 10^-8; actual acceptance also uses error estimates/refinement; Section 14.4 |
| Source coefficient import | Exact decimal/index coverage and digests for both arrays; Sections 20 and 22.1 |
| Polynomial versus nonlinear closure | Report the measured discrepancy; do not substitute a numerical-regression tolerance or require backend equality |
| New derivative, interpolation, or transport references | Declare the method, inputs, comparison rule, and convergence evidence before claiming a pass; no invented universal physical cutoff |

The nonlinear residual tolerance and the fixture comparison tolerance serve different purposes: meeting one does not establish the other. An absolute tolerance, if needed near a mathematically zero quantity, must be declared for that test and justified by its normalization. Do not loosen a tolerance merely to accept a failing implementation. If a different physical-constant set or a scientifically distinct model is intentionally selected, it needs separately identified references and provenance rather than an unexplained update of these fixtures.

The A–G fixtures are the normal deterministic nonlinear/polynomial regression set. Use the retained 300-state dataset for a reproducibility audit or numerical-release check when appropriate; do not make that full audit run inside every coefficient call. Newly implemented NLGC-E and externally supplied-perpendicular cases need independently generated references from their defining equations. Existing NLGCE-N values cannot be reused as if those different closures were identical.

Provide one aggregate runner with readable PASS, FAIL, and SKIP records and a machine-readable result file. Each record identifies the test, model, fixture/source, requested controls, measured discrepancy/residual, and failure reason when applicable. A required supported test failure gives a failing process exit status. SKIP is permitted only for an explicitly unavailable optional dependency, host consumer, or deferred extension and must state why; a skipped required test cannot satisfy a stage's acceptance gate. Tests for individual models should also be runnable separately so Codex can verify an edited component without repeatedly running unrelated transport simulations.

### 23.6 Codex execution and continuation record

At the start of implementation, create IMPLEMENTATION_STATUS.md in the agreed library/project location. It must persist enough evidence for a later Codex session to resume without inventing prior progress.

| Recorded field | Required content |
|---|---|
| Input contract | Specification revision, reference-data identity/digests, and any approved scientific departures |
| Repository decisions | Actual source/header/adapter paths, language standard, build targets, and numerical dependencies |
| Stage state | PD00–PD12 marked pending, in progress, blocked, or complete; never mark a stage complete without its gate evidence |
| Supported inventory | Each implemented stable ID, its capabilities, and explicitly deferred/unsupported options |
| Acceptance evidence | Exact commands, exit statuses, test report paths, backend-specific tolerances, and actual outcomes |
| Change identity | Actual commit/revision when available; otherwise an explicit uncommitted working-tree state |
| Known gaps | Missing derivatives, host bindings, unsupported limits, numerical failures, and incomplete documentation with concrete reasons |
| Next action | The first unfinished prerequisite and the files/commands needed to continue |

Codex should implement the first unfinished stage whose dependencies have passed, run its focused checks, inspect failures, and update this record. Broaden testing when a change affects shared behavior or leaves an unresolved numerical concern. Keep scientific changes reviewable: a new spectrum convention, altered closure, changed parameter, or new physical preset requires an explicit change to the model description and its references. Routine software refactoring must preserve the backend fixtures and data identity.

An unsupported or partially implemented model must fail selection explicitly. A stub returning an arbitrary positive coefficient, a silently substituted backend, or a test that only repeats the implementation formula is insufficient evidence. Meaningful verification uses independently calculated references, exact limiting identities, convergence/refinement, failure behavior, and actual consumer integration where available.

### 23.7 Completion checklist and handoff instruction

The fully specified standalone coefficient library is complete only when its declared supported inventory has passed the applicable gates; host integration is complete only for consumers that have been built and tested. The release checklist is:

- Every Section 3 identifier within the first-release scope has a validated schema, implementation, documented domain, and meaningful verification. A missing planned backend remains unfinished work; it must not disappear from the scope merely by being omitted from the registry. The provider adapters expose their evaluated closure.
- All supplied coefficient literals and reproducibility assets retain their identity. Polynomial and nonlinear backends reproduce their own fixtures; the corrected multiplicative a_x and slab-only L_c are covered by regression tests.
- Required SI/species conversions, optional inputs, per-point failures, numerical convergence, gradients, provenance, and state revisions behave as documented.
- Supplied spectra/tables have explicit normalization, coverage/tail/boundary rules, and cannot acquire undocumented floors, extrapolation, calibration, or universal bounds.
- The standalone core builds independently of background/transport ownership, and each delivered host binding uses that core. Derivative-dependent operations require all needed coefficient and geometric derivatives.
- The test runner distinguishes numerical failure, validity diagnostics, and unavailable extensions. The roadmap/status record describes actual results and any unresolved gate.
- Public comments, README, numerical documentation, examples, installation, and publication reporting describe the same implemented behavior.

**Codex handoff:** use this document and its verified companion bundle as the scientific contract; inspect the target repository; start at PD00; implement and verify the stages in dependency order; maintain IMPLEMENTATION_STATUS.md; and finish with the changed files, supported model inventory, actual build/test evidence, and any explicitly blocked or deferred work. Implement the library without filling missing physical inputs with invented numbers or changing the published closures to match one another.
