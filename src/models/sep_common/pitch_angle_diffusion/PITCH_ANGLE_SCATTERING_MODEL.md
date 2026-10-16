# Pitch-angle scattering models for SEP and GCR focused-transport particle pushers

**File:** PITCH_ANGLE_SCATTERING_MODEL_AMPS_v1_2.md  
**Specification version:** 1.2  
**Date:** 13 October 2026  
**Literature cutoff:** 12 October 2026; principal review interval 2006–2026, with pre-2006 foundations labelled.  
**Companion documents:** MEAN_FREE_PATH_MODEL.md (version 1.2) and PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md (revision 1.4).  
**Purpose:** a model-library specification for pitch-angle scattering of solar energetic particles (SEPs) and galactic cosmic rays (GCRs) in focused-transport particle pushers (stochastic-differential-equation, Monte Carlo and guiding-centre codes), and a source for the methods sections of publications.

This document collects the pitch-angle diffusion coefficients D_μμ that have been derived or used for heliospheric SEP and GCR transport, mainly in work published between 2006 and 2026, together with what a particle pusher needs in order to apply them. It covers the focused transport equations that pushers integrate (Section 2), the core definitions and conventions (Section 3), empirical and phenomenological D_μμ shapes (Section 5), quasilinear D_μμ for magnetostatic slab turbulence (Section 6), dynamical turbulence and the dissipation range (Section 7), wave-based D_μμ with propagation direction, helicity, wave frames and dispersive waves (Section 8), anisotropic Alfvénic and compressible turbulence (Section 9), nonlinear and second-order theories (Section 10), D_μμ measured in test-particle simulations and the estimators used (Section 11), adiabatic focusing and the D_μμ–λ∥ relation (Section 12), D_μμ computed from measured spectra (Section 13), self-generated waves and shock-associated scattering (Section 14), pitch-angle-resolved GCR models (Section 15), and the numerical formulation in a pusher: the Itô stochastic differential equation, importance sampling, integration schemes, the boundaries at μ = ±1, the point μ = 0, time-step rules, alternative kernels and frames (Section 16). Each formula is given as printed in its source, with the source location and the verification level of the reading (Section 1.3). Results computed in this review are labelled "derived" and were checked numerically; evaluations of published formulas at published inputs are labelled as such and are not presented as published results.

Where sources disagree, contain apparent typographical errors, or could not be read, the document records the discrepancy or the gap and does not fill it by assumption. Section 19 lists 128 printed discrepancies and reading conflicts, 67 decisions that the sources leave to the user or the host code, the information that was sought and not found, and 15 notes on the companion specifications. Several printed items in code descriptions must not be copied into an implementation, among them a forward pitch-angle SDE without the drift ∂D_μμ/∂μ, a Wiener increment written as dt·Λ, a Milstein update without its leading μ_t, and an importance-sampling drift with an undefined constant and no weight (Section 16.10); printed forms that are dimensionally inconsistent or not real on part of the pitch-angle interval are records, not runtime coefficients (Section 4.1). Section 17 specifies the library interface, Section 20 states the factual consequences of the sources for model selection, and Section 26, the final section, orders the implementation work. Version 1.1 incorporates the confirmed findings of an external review of version 1.0 (Section 25.6); among them, the diffusion-approximation coefficient of Effenberger & Litvinenko (2014) is distinguished from the long-time spatial variance that a focused-transport pusher produces (Section 12.3), and the stationary pitch-angle exponential is attributed to the correct equation.

The version 1.1 author describes the intended companion bundle PITCH_ANGLE_SCATTERING_MODEL_DATA.zip as containing the model registry (172 model IDs with source, equation, specification status and runtime state), every transcribed parameter set with provenance, 89 verification fixtures with 652 stored values and their generator, the derivation records with the scripts that check every derived statement, a self-checking script and an independent set of mathematical checks (Section 24). PITCH_ANGLE_SCATTERING_MODEL_additions.bib contains BibTeX entries only for the 44 cited references that are not in the supplied bibliography (Section 23).

**Implementation revision 1.2.** The prior physical corrections are retained. Sections 16.12 and 17 now define the scalar/tensor sampled-measure drift, state-dependent clock conversion, frame limits, trajectory streaming, adaptive noise, wave-channel bookkeeping and typed failure contract. Section 25.7 gives the current review findings and evidence; Section 26 provides 20 detailed work packages and 21 acceptance cases. This is an implementation specification, not certification that the AMPS library or all inventory models have already been implemented.

**Data availability.** The updated description refers to a v1.1 canonical archive, but the available uploaded archive is v1.0 (168 IDs / 83 fixtures / 567 values). The intended v1.1 data counts and author-reported checks are retained as historical claims and are not re-certified here. The accompanying v1.2 implementation supplement contains separate mathematical checks, a compilable declaration and the detailed roadmap; it does not relabel or reconstruct the canonical v1.1 archive.

## Contents

1. [Scope, companion specifications, and conventions](#1-scope-companion-specifications-and-conventions)
2. [The focused transport equation used by particle codes](#2-the-focused-transport-equation-used-by-particle-codes)
3. [Core definitions and conventions](#3-core-definitions-and-conventions)
4. [Model inventory and runtime states](#4-model-inventory-and-runtime-states)
5. [Empirical and phenomenological D_μμ shapes](#5-empirical-and-phenomenological-d_μμ-shapes)
6. [Quasilinear D_μμ for magnetostatic slab turbulence](#6-quasilinear-d_μμ-for-magnetostatic-slab-turbulence)
7. [Dynamical turbulence and the dissipation range](#7-dynamical-turbulence-and-the-dissipation-range)
8. [Wave-based D_μμ: propagation direction, helicity, wave frames and dispersive waves](#8-wave-based-d_μμ-propagation-direction-helicity-wave-frames-and-dispersive-waves)
9. [Anisotropic Alfvénic and compressible (fast-mode) turbulence](#9-anisotropic-alfvénic-and-compressible-fast-mode-turbulence)
10. [Nonlinear and second-order theories of D_μμ](#10-nonlinear-and-second-order-theories-of-d_μμ)
11. [D_μμ from test-particle simulations and how it is estimated](#11-d_μμ-from-test-particle-simulations-and-how-it-is-estimated)
12. [Adiabatic focusing and the D_μμ–λ∥ relation](#12-adiabatic-focusing-and-the-d_μμλ-relation)
13. [D_μμ computed from measured or prescribed spectra](#13-d_μμ-computed-from-measured-or-prescribed-spectra)
14. [Self-generated waves and shock-associated scattering](#14-self-generated-waves-and-shock-associated-scattering)
15. [Pitch-angle-resolved GCR transport models](#15-pitch-angle-resolved-gcr-transport-models)
16. [Pitch-angle scattering in a particle pusher: numerical formulation](#16-pitch-angle-scattering-in-a-particle-pusher-numerical-formulation)
17. [Library interface](#17-library-interface)
18. [Verification fixtures](#18-verification-fixtures)
19. [Discrepancy, decision and gap register](#19-discrepancy-decision-and-gap-register)
20. [Model selection](#20-model-selection)
21. [Publication methods text](#21-publication-methods-text)
22. [References](#22-references)
23. [Bibliography comparison](#23-bibliography-comparison)
24. [Companion data](#24-companion-data)
25. [Verification and revision record](#25-verification-and-revision-record)
26. [Codex implementation roadmap (final section)](#26-codex-implementation-roadmap-final-section)

## 1. Scope, companion specifications, and conventions

### 1.1 Scope

This specification describes the pitch-angle scattering models that SEP and GCR transport codes use, mainly in work published between 2006 and 2026, and the numerical rules for applying them in a particle pusher that integrates a focused transport equation. Pre-2006 papers appear only where a model in use is defined through them; they are labelled as foundations.

The document is a model library: each model has an identifier (Section 1.4), a source, the formula as printed, the definition of every symbol as stated in the source (frame of μ, units, spectrum normalization, operator convention), the published parameter values with their locations, and, where the source states them, the domain, singular points, boundary treatment and time-step rules. A model whose source leaves something open is recorded with the open item rather than completed (Section 19).

It does not cover:

- perpendicular diffusion, except where a pitch-angle-resolved model prints a μ-dependent perpendicular coefficient (Section 15.2);
- spatial diffusion closures, which are specified in the companion documents (Section 1.2);
- full-orbit particle tracing in synthetic turbulence as a transport method; such simulations appear only as sources of D_μμ and of estimators (Section 11);
- momentum diffusion, except where it accompanies pitch-angle scattering in a wave frame (Sections 8, 9 and 14).

### 1.2 Companion specifications

Two companion documents are referred to by section and equation number and are not repeated:

- **MFP v1.2** — *Mean-free-path models for SEP and GCR transport in the heliosphere*, MEAN_FREE_PATH_MODEL.md, version 1.2 (9 October 2026). It specifies the λ definitions (λ∥, λ_r, λ_rr; its Section 1.3), the SEP presets that fix λ (its Section 5), the pitch-angle shapes and their exact normalization to λ∥ (its Section 6, Eq. 11 there being Eq. (23) here), the closed-form QLT λ∥ (its Section 7), shock-region and GCR prescriptions (its Sections 8–9), turbulence inputs (its Section 11) and the Itô SDE for pure pitch-angle scattering (its Section 21.5, Eq. 46). References to "MFP v1.2 D-n" and "MFP v1.2 U-n" are items of its register (its Section 14).
- **PARALLEL rev. 1.4** — *Parallel diffusion coefficient models for SEP and GCR transport*, PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md, revision 1.4 (8 October 2026). It specifies kinematics and units (its Section 2), the diffusion-limit relation between D_μμ and κ∥ (its Section 7.1, Eq. 20), slab QLT and the spectrum conventions (its Section 8), nonlinear closures (its Sections 9–11) and turbulence and wave adapters (its Section 12).

Where this document and a companion disagree, the disagreement is recorded in Section 19.4 and is not resolved silently. An implementation pins the exact revisions of both companion documents and checks every cited section and equation number against them before implementing a delegated equation; a delegated path whose companion text is not available is reported as `BLOCKED_DEPENDENCY` (Section 4.1) rather than reconstructed.

Kinematics and unit conversions follow PARALLEL rev. 1.4, Section 2 (pc = [T(T + 2mc²)]^{1/2}, v = pc²/(T + mc²), rigidity 𝓡 = pc/|q|, r_L = 𝓡/(cB₀)), and the constants of MFP v1.2, Section 2 (exact c, e and AU; CODATA 2018 masses). The fixtures of Section 18 state the constants they use.

### 1.3 Status labels

| Label | Meaning |
|---|---|
| Full text ×n | Read in the full text; n readings (fetches, possibly of the same copy or the same text extraction) agreed (n = 1 means not double-checked). Readings of the same extraction confirm the transcription of that extraction, not the typeset original; only "two renderings" below confirms an item in two distinct versions |
| One rendering / two renderings | The number of distinct versions read (for example the arXiv PDF and the ar5iv HTML, or the arXiv and publisher versions); "two renderings agree" means the item was confirmed in both |
| Abstract only | Only the abstract or metadata was accessible |
| Secondary | Taken from a named citing paper |
| As printed | Transcribed exactly, including apparent typographical errors; not corrected |
| Not legible / structure only | The rendering was garbled; only the legible part or the structure of the equation is given, and nothing is reconstructed |
| Summarized reading; tool summary, tool assembly, tool-reconstructed | The item was read through a text-extraction tool that summarized it, assembled it from fragments or reconstructed it; it is given as structure only and is not verified as printed |
| Not found | Not located in any accessible source; nothing is substituted |
| Derived | An algebraic or numerical consequence computed in this review, not a published statement; every derived statement was checked numerically (Section 25) |
| Evaluation of printed formula | A published formula evaluated at published inputs; not a value the authors published |
| **Implementation requirement** | A design rule for the library that follows from the sources; not a physical statement |

Every number in the companion data carries one of these labels (Section 24). "The readings" are the records made when the sources were read for this review (Section 25.2); "not recorded in the readings" means that an item was not captured, not that the source lacks it.

### 1.4 Identifiers and cross-references

- **Equations** of this document are numbered and cited with parentheses, "Eq. (n)". Equations of the companion documents are written "MFP v1.2 Eq. n" and "PARALLEL rev. 1.4 Eq. n". A source's own equation numbers are written without parentheses, "Eq. n", normally after the source name (for example "Zhang et al. 2023, Eq. 6"); inside a paragraph about one source, "Eq. n" refers to that source.
- **References** are cited as [Pn] with the number of first citation; Section 22 lists them, and the companion file `source_key_map.json` maps each [Pn] to its source key (for example `Wijsen2019`), DOI and arXiv identifier.
- **Model IDs** carry a prefix by topic: FTE- (Section 2), PA- (Section 5; the MFP v1.2 shape IDs PA-QFORM, PA-EPS, PA-ISO, PA-EPREM, PA-KOLMO, PA-DROGE-VA and PA-AMPS-I…V are reused for the same shapes), QLT- (6), DYN- (7), WAVE- (8), CMP- (9), NL- (10), SIM- (11), FOC- (12), MEAS- (13), SGW- (14), GCR-PA- (15), and SDE-, MC- and VAL- for the numerical kernels and validation tests of Section 16. Section 4 lists them.
- **Register items** are D-n (printed discrepancies), U-n (decisions for the user or host code) and C-n (notes on the companion specifications), all in Section 19. Items of the companion registers always carry their prefix ("MFP v1.2 D-2").
- **Fixtures** are F-<prefix>-nn (Section 18). Check labels are written with their derivation record (for example "DR-06 C4"; also S8-n and S9-n of DR-08) and are distinct from the register items C-n. Monte Carlo runs are MCR-1 to MCR-7 and MCR-1b (labelled MC-1 to MC-7 inside DR-16) and are distinct from the model IDs MC-ROT, MC-ROT-HEMI, MC-REDRAW and MC-SPHERE-FLORINSKI13. The derivation records DR-02, DR-05, DR-06, DR-08, DR-10, DR-14 and DR-16 of the companion bundle (Section 24) contain the scripts and their output.

### 1.5 Notation used throughout

Each section defines its symbols as its sources do. The following are common to all sections:

| Symbol | Meaning |
|---|---|
| μ | cosine of the pitch angle, measured relative to the mean magnetic field **B** unless the source states otherwise (Zhang et al. 2023 measure it relative to the outward field direction; U-5) |
| v, p, γ | particle speed, momentum and Lorentz factor in the frame in which μ is defined (Sections 2.7 and 16.8 list the frames) |
| D_μμ | the coefficient of ∂/∂μ(·∂f/∂μ) in Eq. (21), unit s⁻¹; D′_μμ = ∂D_μμ/∂μ |
| ν | scattering frequency of isotropic scattering, Eq. (24) |
| λ∥, κ∥ | parallel mean free path and parallel diffusion coefficient, Eq. (23) |
| L | focusing length; most sources use 1/L = −b̂·∇ln B (Section 2.1; the second printed form of He & Wan 2015 Eq. 3 has the opposite sign, D-1) |
| ξ | dimensionless focusing strength λ₀/L, with λ₀ the mean free path without focusing (Sections 12 and 16.9) |
| q, H, h, ε | shape and gap parameters of the empirical forms (Section 5); the source's own symbol is kept |

### 1.6 Normative SI layer for executable formulas

The sources print their formulas in Gaussian units (Sections 6, 7, 9, 10 and most of 8 and 14) or in SI (Sections 8.3, 13 and 14 in part), and several print no unit system. The library evaluates every executable formula in SI; a published formula is a transcription, not the runtime expression, and enters the runtime only through the adapter named in its model record. The adapters are the following statements, each derived from the definitions and checked (F-PUSH-19 for the kinematics and the clock):

| Quantity | Normative SI representation in the library |
|---|---|
| Momentum, energy | p in kg m s⁻¹; kinetic energy T and total energy in J; eV, MeV and "energy per nucleon" are input and output conversions only. For an ion, T is the kinetic energy of the ion; an input per nucleon is converted with the stated mass number. p = [T(T + 2mc²)]^{1/2}/c, v = pc²/(T + mc²), γ = 1 + T/(mc²) (PARALLEL rev. 1.4, Section 2) |
| Gyrofrequency | signed Ω = qB₀/(γm) in rad s⁻¹ with the signed charge q and the relativistic mass γm; a resonance condition that uses |Ω| or omits γ says so (D-99, U-28, U-47); the Gaussian Ω = qB₀/(γmc) maps onto this value |
| Larmor scale | r_L = p/(|q|B₀) = v/|Ω| (the "maximal Larmor radius"); the gyroradius of a particle with pitch cosine μ is r_L(1 − μ²)^{1/2}; the Gaussian pc/(|q|B₀) is the same quantity |
| Rigidity | 𝓡 = pc/|q| in V; r_L = 𝓡/(cB₀); a rigidity scaling is written as (𝓡/𝓡_ref)^a with 𝓡_ref stated (D-21, U-9) |
| Spectral coordinates | k in rad m⁻¹; spacecraft-frame frequency f in Hz, ω = 2πf; every factor 2π in a k–f conversion is declared (Section 13.1) |
| Magnetic spectra | a one-sided transverse power spectral density P in T² m with ∫₀^∞P dk = the transverse variance in T²; the slab normalization δB² = 8π∫₀^∞g dk∥ of Section 6.1 defines g in the units of P/(8π), and the combinations g/B₀² and P/B₀² (length) are the same in Gaussian and SI units; sidedness and component content are declared per spectrum (Sections 6.7, 13.9) |
| Wave energy | per-k wave energy density in J m⁻², per-ln k in J m⁻³; the magnetic energy density of a fluctuation is δB²/(2μ₀), and the total energy of Alfvén waves under the equipartition the source assumes is δB²/μ₀ (the Gaussian B²/(8π) and B²/(4π)); the factor two between them is conditional on that assumption and is not a universal adapter (Section 14.2) |
| Transport coefficients | D_μμ in s⁻¹; λ∥ in m; κ∥ in m² s⁻¹; a positive physical elapsed time in s, with the clock of the coefficient declared (Eq. (173)) |
| Species | particle mass and signed charge per species; electrons and ions are tested with the signed Ω (U-28) |

The fixtures of Section 18 state the constants they use (exact c, e and au; CODATA 2018 masses and μ₀).

## 2. The focused transport equation used by particle codes

This section collects, as printed, the focused transport equations (FTE) that the reviewed particle codes integrate, with the frame of position (x, z or s), pitch-angle cosine μ and momentum p as each source states it, and the operator convention of the scattering term. Coefficients are in Sections 5–15 and schemes in Section 16; here only the direction of integration is recorded.

Most sources print Eq. (21) on the right, or on the left with a minus sign (the same operator); Ruffolo (1995) and EPREM use other factors (Sections 2.3, 2.4). The contraction is $\vec{ab}:\vec{cd}=a_ib_jc_jd_i$, so $\hat b\hat b:\nabla\vec V=b_ib_j\partial_jV_i$ [P1]. "Not stated" (n.s.) means the text read does not state it.

### 2.1 One-dimensional form along a field line

Effenberger & Litvinenko (2014), Eq. 1; Full text ×2 (ar5iv and arXiv PDF agree) [P2]:

$$
\frac{\partial f_0}{\partial t}+\mu v\frac{\partial f_0}{\partial z}+\frac{v}{2L}(1-\mu^2)\frac{\partial f_0}{\partial\mu}=\frac{\partial}{\partial\mu}\left(D_{\mu\mu}\frac{\partial f_0}{\partial\mu}\right)
\tag{1}
$$

Here f₀ is the gyrotropic phase-space distribution, z the distance along the mean field, v the constant speed (no solar wind, no energy change) and L = −B/(∂B/∂z) a constant focusing length; the operator is Eq. (21). The paper uses D_μμ = D₀(1−μ²) (Eq. 2; Eq. (24) with ν = 2D₀) and the linear density along the flux tube f = exp(z/L)f₀ (Eq. 3). Its Itô SDEs are "completely equivalent to the Fokker-Planck equation for the linear density" and are integrated forward in time with reflecting boundaries at μ = ±1 (Section 16). κ∥(ξ), ξ = λ₀/L (Eq. 8), is in Section 12.

The same equation is printed by Agüeda et al. (2014, Eq. 2; ×1; with a source; "convection and adiabatic deceleration neglected") [P3] and, for GCR protons, by Strauss et al. (2022, Eqs. 2–3; ×2) [P4]. Qin et al. (2005, Eq. 1) introduce it "in solar wind frame"; their extraction is garbled, and only the structure, with focusing written −[v(1−μ²)/(2B)](dB/dz)∂f/∂μ, is legible [P5].

**Focusing length as printed.** All rows that can be checked give a positive focusing rate where B decreases in the direction μ = +1, except the second printed form of He & Wan (2015, Eq. 3), whose first form agrees with the other rows; two rows cannot be checked.

| Definition as printed | Sources |
|---|---|
| L = −B/(∂B/∂z), or 1/L = −(1/B)dB/dz, or L = B/(−∂B/∂z); μ̇ contains +v(1−μ²)/(2L) | Effenberger & Litvinenko 2014; Agüeda et al. 2014; Ruffolo 1995 (def., Eq. 3); Strauss et al. 2022, Eq. 3; Qin et al. 2005, Eq. 8 (tool rendering, "= (1/A)(dA/dz)"); Kartavykh et al. 2015; Laitinen et al. 2016; Merten & Aerdker 2024; Tenishev et al. 2022, Eq. 4 |
| L⁻¹ = ∇·b̂, with −∂/∂μ[(1−μ²)vf/(2L)] on the right | Strauss & Fichtner 2015, Eq. 1; Strauss et al. 2017, Eq. 2; Strauss & le Roux 2019, Eq. 19; van den Berg et al. 2026, Sect. 2 |
| no L; −v b̂·∇ln B written in dμ/dt | Zhang et al. 2023, Eq. 3; Kubo et al. 2015, Eq. 3; EPREM (Young & Poduval 2025, Eq. A1) |
| L = (b̂·∇ln B₀)⁻¹, with +(1−μ²)/2·[−v/L + …]∂f/∂μ on the left | Qin & Wu 2020, Eq. 2 |
| text "L = (z·∇lnB)⁻¹" (∇ lost); the legible term in Eq. 3 is −(v/B)∂B/∂z | He et al. 2011 |
| dμ/dt = (1−μ²)/2·[−(v/B)∂B/∂z + …] = (1−μ²)/2·[v/L + …] with L = (z·∇ln B)⁻¹: **the two printed forms differ in sign** (D-1) | He & Wan 2015, Eq. 3 (two renderings); He & Wan 2019 and He et al. 2017 print the v/L form with the same L (whether they also print the first form was not recorded) |
| definition not captured; −(1−μ²)/2·[−v/L + …]∂f/∂μ on the right | Wang et al. 2014, Eq. 1 |
| "given by Eq. 22" (appendix not read) | van den Berg et al. 2020, Eq. 9 |

Sources: [P6, P7, P8, P9, P10, P11, P12, P13, P14, P15, P16, P17, P18, P19, P20, P21, P22, P23].

*Derived (check F-FTE-01).* If ∇·B = 0, ∇·b̂ = −b̂·∇ln B, so L⁻¹ = ∇·b̂ equals −(1/B)∂B/∂s. For a Parker spiral it is V(2V²+w²)/[r(V²+w²)^{3/2}], w = Ωr sin θ, the focusing length of MFP v1.2 Section 6.9. At r = 1 au, θ = 90°, V = 400 km s⁻¹, Ω = 2.86 × 10⁻⁶ rad s⁻¹ this gives L = 0.998552110968506 au; He & Wan's printed definition of L gives −0.998552110968506 au, so v/L with that L reverses the focusing term of the first form of their own Eq. 3.

**Flux-tube density form.** Qin et al. (2005, Eq. 7) rewrite the equation for g ≡ A(z)f, A ∝ 1/B the flux-tube area, because "we want to simulate the stochastic particle transport in a time-forward manner". Only its structure is legible (tool rendering): a forward Fokker–Planck equation with streaming in z, μ drift v(1−μ²)/(2L) + ∂D/∂μ and μ diffusion ∂²(Dg)/∂μ²; the μ of the streaming coefficient was lost ("vg"). van den Berg et al. (2020, Eq. 9; ×2: PDF and abs-page LaTeX) print the "Roelof 1969 form, without advection or energy losses", s the position along the field:

$$
\frac{\partial f}{\partial t}+\frac{\partial}{\partial s}\left[\mu vf\right]+\frac{\partial}{\partial\mu}\left[\frac{(1-\mu^2)v}{2L(s)}f\right]=\frac{\partial}{\partial\mu}\left[D_{\mu\mu}\frac{\partial f}{\partial\mu}\right]
\tag{2}
$$

Lang et al. (2024, Eq. 1; ∂/∂s_p) and Merten & Aerdker (2024, Eq. 1; with a source S) print the same form [P24, P9].

*Derived (check F-FTE-04).* If F ≡ Af₀ with d ln A/ds = 1/L, Eq. (2) for F is A times Eq. (1) for f₀, for any L(s). If f in Eq. (2) is the phase-space density, the two equations differ by μvf/L. The readings of the three 1-D printings do not state which density is meant. Strauss & Fichtner (2015) and Strauss et al. (2017) write the streaming term as the 3-D divergence −∇·(μv b̂ f) with L⁻¹ = ∇·b̂; that form is identical to Eq. (1) for the phase-space density (a special case of Eq. (13)).

**Implementation requirement.** A time-forward particle code samples a density obeying a Fokker–Planck equation in its own coordinates; along a field line this is the linear density Af, A ∝ 1/B [P5, P2]. The phase-space density is the sampled linear density divided by A (up to a constant). Time-backward codes integrate the advective form for f (Section 16).

### 2.2 Three-dimensional forms in the solar-wind frame

**Zhang et al. (2023)** [P15]. Verification: Eq. 1 Full text ×2 (ar5iv, PDF); Eq. 3 Full text ×1 (ar5iv); Eqs. 2 and 4 were relayed by the reading tool in summarized form.

$$
\frac{\partial f}{\partial t}-\nabla\cdot\boldsymbol{\kappa}_{\perp}\cdot\nabla f+\left(v\mu\hat{\mathbf{b}}_{o}+\mathbf{V}+\mathbf{V}_{d}\right)\cdot\nabla f-\frac{\partial}{\partial\mu}D_{\mu\mu}\frac{\partial f}{\partial\mu}+\frac{d\mu}{dt}\frac{\partial f}{\partial\mu}+\frac{dp}{dt}\frac{\partial f}{\partial p}=Q_{0}(t,\mathbf{x},\mu,p)
\tag{3}
$$

$$
\frac{d\mu}{dt}=-\frac{\left(1-\mu^{2}\right)v}{2}\hat{\mathbf{b}}_{o}\cdot\nabla\ln B+\frac{\mu\left(1-\mu^{2}\right)}{2}\left(\nabla\cdot\mathbf{V}-3\hat{\mathbf{b}}\hat{\mathbf{b}}:\nabla\mathbf{V}\right)
\tag{4}
$$

Eq. 4 (summarized rendering): dp/dt = −[(1−μ²)/2·(∇·V − b̂b̂:∇V) + μ² b̂b̂:∇V] p. Eq. 2 gives V_d; it is available only in summarized rendering and is not reproduced here.

Here p is "particle rigidity (proportional to momentum for a given charged particle species)" and μ the cosine of the pitch angle "to the outward magnetic field line". b̂_o is the outward unit vector, B = B b̂ = ±B b̂_o; this choice keeps f continuous across polarity reversals. V is the plasma velocity, κ⊥ the perpendicular tensor and Q₀ the seed source. The frame of p and μ is not stated; "We compute particle trajectories in the coordinates corotating with the sun" was read in the PDF only. Solution: time-backward SDEs with importance sampling.

**He et al. (2011); He & Wan (2015)** [P20]. He & Wan (2015), Eqs. 1–3: Full text ×1 (ar5iv, two fetches); Eq. 3 confirmed in a second rendering (arXiv PDF) by a later check.

$$
\frac{\partial f}{\partial t}+\mu v\frac{\partial f}{\partial z}+\mathbf V^{sw}\cdot\nabla f+\frac{dp}{dt}\frac{\partial f}{\partial p}+\frac{d\mu}{dt}\frac{\partial f}{\partial\mu}-\frac{\partial}{\partial\mu}\left(D_{\mu\mu}\frac{\partial f}{\partial\mu}\right)-\frac{\partial}{\partial x}\left(\kappa_{xx}\frac{\partial f}{\partial x}\right)-\frac{\partial}{\partial y}\left(\kappa_{yy}\frac{\partial f}{\partial y}\right)=Q(\mathbf x,p,t)
\tag{5}
$$

$$
\frac{dp}{dt}=-p\left[\frac{1-\mu^{2}}{2}\left(\frac{\partial V^{sw}_{x}}{\partial x}+\frac{\partial V^{sw}_{y}}{\partial y}\right)+\mu^{2}\frac{\partial V^{sw}_{z}}{\partial z}\right],\qquad
\frac{d\mu}{dt}=\frac{1-\mu^{2}}{2}\left[-\frac{v}{B}\frac{\partial B}{\partial z}+\mu\left(\cdots\right)\right]=\frac{1-\mu^{2}}{2}\left[\frac{v}{L}+\mu\left(\frac{\partial V^{sw}_{x}}{\partial x}+\frac{\partial V^{sw}_{y}}{\partial y}-2\frac{\partial V^{sw}_{z}}{\partial z}\right)\right]
\tag{6}
$$

L = (z·∇ln B)⁻¹ as printed; (⋯) in the first form stands for the same velocity-gradient bracket as in the second. With the printed L the two forms of dμ/dt differ in the sign of the focusing term (D-1). z is along the average Parker field and x, y are perpendicular to it. Drifts are neglected, and the frame is not stated. He et al. (2011, Eqs. 1–3) give the same equation with right-hand side 0 (the source is in their Eq. 7); their extraction is garbled [P19]. He & Wan (2019) and He et al. (2017) repeat the form [P21, P22]. All use time-backward SDEs.

*Derived (check F-FTE-02).* With z along b̂, ∂V_x/∂x + ∂V_y/∂y = ∇·V − b̂b̂:∇V and ∂V_z/∂z = b̂b̂:∇V. The dp/dt and the μ(…) term of Eq. (6) are then those of Zhang et al. (2023). The focusing term agrees with the first printed form, −(v/B)∂B/∂z = −v b̂·∇ln B, and with the second form only if v/L is read as −v b̂·∇ln B, the opposite of the printed L.

**Wijsen et al. (2019), PARADISE** [P25]. Eq. 1: Full text ×2 (ar5iv, PDF). Eqs. 2–5: Full text ×1 (ar5iv).

$$
\frac{\partial f}{\partial t}+\frac{d\mathbf{x}}{dt}\cdot\nabla f+\frac{d\mu}{dt}\frac{\partial f}{\partial\mu}+\frac{dp}{dt}\frac{\partial f}{\partial p}=\frac{\partial}{\partial\mu}\left(D_{\mu\mu}\frac{\partial f}{\partial\mu}\right)+\nabla\cdot\left(\boldsymbol{\kappa}_{\perp}\cdot\nabla f\right)
\tag{7}
$$

$$
\begin{aligned}
\frac{d\mathbf x}{dt}&=\mathbf{V}_{\rm sw}+\mathbf{V}_{d}+\mu v\mathbf{b},\qquad
\frac{d\mu}{dt}=\frac{1-\mu^{2}}{2}\left(v\nabla\cdot\mathbf{b}+\mu\nabla\cdot\mathbf{V}_{\rm sw}-3\mu\,\mathbf{b}\mathbf{b}:\nabla\mathbf{V}_{\rm sw}-\frac{2}{v}\mathbf{b}\cdot\frac{d\mathbf{V}_{\rm sw}}{dt}\right),\\
\frac{dp}{dt}&=\left(\frac{1-3\mu^{2}}{2}(\mathbf{b}\mathbf{b}:\nabla\mathbf{V}_{\rm sw})-\frac{1-\mu^{2}}{2}\nabla\cdot\mathbf{V}_{\rm sw}-\frac{\mu}{v}\mathbf{b}\cdot\frac{d\mathbf{V}_{\rm sw}}{dt}\right)p,\qquad
\mathbf{V}_{d}=\frac{vp}{QB}\left[\frac{1-\mu^{2}}{2}\left((\nabla\times\mathbf{b})_{\parallel}+\frac{\mathbf{b}\times\nabla B}{B}\right)+\mu^{2}(\nabla\times\mathbf{b})_{\perp}\right]
\end{aligned}
\tag{8}
$$

Frames: "the spatial coordinate is measured in the fixed inertial frame of the observer", "whereas the momentum and pitch-angle cosine are measured in a reference frame co-moving with the solar wind" (also Husidic et al. 2024 [P26]). The cross-field term "is not part of the standard FTE". The d/dt in b·dV_sw/dt is not defined in the portion read. Integration is time-forward: position by an Itô SDE, scattering by the small-angle method of Agüeda et al.; the standard μ SDE (Eq. 13) is printed but not used.

*Derived (check F-FTE-02).* With b·dV_sw/dt = 0 and ∇·b = −b·∇ln B, the coefficients in Eq. (8) equal those of Zhang et al. (2023).

**Further printings of these coefficients.**
- Laitinen et al. (2016, Eq. 1; Full text ×2) [P8] use the same coefficients with L = −B/(∂B/∂s) and ∇·κ̂∇f. They omit the dV_sw/dt terms and the drifts, and do not state the frame.
- Kubo et al. (2015, Eqs. 1–3; ×1) [P16] use the same coefficients in index notation; "the particle position and momentum are measured in the co-rotating frame of the Sun and the solar wind frame, respectively".
- Wang & Qin (2015): "x is the position in a non-rotating heliographic coordinate system", p, μ, v "in the solar wind frame" [P27]; Wang et al. (2014): μ "in the plasma reference frame".

### 2.3 Conservative and mixed-frame forms

**Ruffolo (1995), foundation** [P6]. Eqs. 3 and 7–9: Full text ×1 (ar5iv); scattering term: Full text ×2 (PDF). Mixed frames: z (arc length along the field) and t fixed; μ′, p′, v′ in the local solar-wind frame; ψ(z) the angle between r̂ and ẑ; 1/L(z) ≡ −(1/B)dB/dz.

$$
\dot z=\mu'v'+\left(1-\mu'^2\frac{v'^2}{c^2}\right)v_{\rm sw}\sec\psi,\quad
\dot\mu'=\frac{v'}{2L(z)}\left[1+\mu'\frac{v_{\rm sw}}{v'}\sec\psi-\mu'\frac{v_{\rm sw}v'}{c^2}\sec\psi\right](1-\mu'^2),\quad
\dot p'=-p'v_{\rm sw}\left[\frac{\sec\psi}{2L(z)}(1-\mu'^2)+\cos\psi\frac{d}{dr}\sec\psi\,\mu'^2\right]
\tag{9}
$$

These are Eqs. 8, 9 and 7; Eq. 10 adds differential convection (tool summary). Eq. 11, a conservative equation for F = d³N/(dz dμ dp) (N particles in a flux tube), was assembled by the tool; only its scattering term is verified, ∂/∂μ (φ(μ)/2) ∂/∂μ, "where φ(μ) is the pitch-angle scattering coefficient" (a first-order E′/E factor reported inside it is not verified). Hence φ = 2D_μμ. Ruffolo's λ ≡ 3D/v, D = (v²/4)∫(1−μ²)²/φ dμ, is Eq. (23) with D_μμ = φ/2; for φ = A|μ|^{q−1}(1−μ²) (Eq. 12) with q = 1.5, both give λ = 2.4 v/A (check F-FTE-06). Solution: finite differences, zero μ flux at μ = ±1.

**The primer: van den Berg, Strauss & Effenberger (2020)** [P1]. Eq. 17: Full text ×2. Eq. 18: Full text ×1 (PDF). The frames are stated as: "the GC position $\vec x$ and time are measured in the observer's frame" and "the momentum p and pitch-cosine μ are measured in the SW frame".

$$
\frac{\partial f}{\partial t}+\frac{\partial}{\partial x_i}\left[\frac{dx_i}{dt}f\right]+\frac{\partial}{\partial p}\left[\frac{dp}{dt}f\right]+\frac{\partial}{\partial\mu}\left[\frac{d\mu}{dt}f\right]=\frac{\partial}{\partial\mu}\left[D_{\mu\mu}\frac{\partial f}{\partial\mu}+D^{\perp}_{\mu j}\frac{\partial f}{\partial x_j}\right]+\frac{\partial}{\partial x_i}\left[D^{\perp}_{i\mu}\frac{\partial f}{\partial\mu}+D^{\perp}_{ij}\frac{\partial f}{\partial x_j}\right]
\tag{10}
$$

$$
\begin{aligned}
\frac{d\vec x}{dt}&=\mu v\,\hat{\vec b}+\vec v_{\rm sw}+\vec v_d,\qquad
\frac{dp}{dt}=p\left\{\frac{1-3\mu^2}{2}\hat{\vec b}\hat{\vec b}:\nabla\vec v_{\rm sw}-\frac{1-\mu^2}{2}\nabla\cdot\vec v_{\rm sw}-\frac{\mu}{v}\hat{\vec b}\cdot\left[\frac{\partial\vec v_{\rm sw}}{\partial t}+(\vec v_{\rm sw}\cdot\nabla)\vec v_{\rm sw}\right]\right\},\\
\frac{d\mu}{dt}&=\frac{1-\mu^2}{2}\left\{v\nabla\cdot\hat{\vec b}+\mu\nabla\cdot\vec v_{\rm sw}-3\mu\hat{\vec b}\hat{\vec b}:\nabla\vec v_{\rm sw}-\frac{2}{v}\hat{\vec b}\cdot\left[\frac{\partial\vec v_{\rm sw}}{\partial t}+(\vec v_{\rm sw}\cdot\nabla)\vec v_{\rm sw}\right]\right\}
\end{aligned}
\tag{11}
$$

Eq. 18 also prints a three-term drift v_d:
- (μp/qB) b̂×[∂b̂/∂t + (v_sw·∇)b̂];
- (pv/qB)[μ²(∇×b̂)_⊥ + (1−μ²)/2·b̂×∇B/B];
- (m/qB) v_sw×[∂(b̂/B)/∂t + (μv b̂ + v_sw)·∇(b̂/B)].

The mixed terms $D^\perp_{i\mu}=D^\perp_{\mu j}$ are "usually neglected". Strauss et al. (2020) print the same structure and frames; their extraction is garbled [P28].

**PARADISE/PARASOL intensity form.** Afanasiev et al. (2025), Eq. 37; Full text ×1 (arXiv v1, two fetches) [P29]:

$$
\frac{\partial j}{\partial t}+\nabla\cdot\left(\frac{d\mathbf{x}}{dt}j\right)+\frac{\partial}{\partial\mu}\left(\frac{d\mu}{dt}j\right)+\frac{\partial}{\partial p}\left(\frac{dp}{dt}j\right)=\frac{\partial}{\partial\mu}\left(D_{\mu\mu}\frac{\partial j}{\partial\mu}\right)+S
\tag{12}
$$

Here j = p²f is the "differential particle intensity" and S is non-zero only at the shock. Eq. 38 is dx/dt = V_sw + μv b̂ (no drifts, no κ⊥: "we only incorporate a pitch-angle diffusion process"). Eqs. 39–40 repeat Wijsen et al. (2019, Eqs. 3–4), except that Eq. 40 prints μ/ν (Greek nu) for μ/v. The frame of μ and p is not stated in the section read.

*Derived (checks F-FTE-03, F-FTE-08).* For the flow of Eq. (8) or Eq. (17) without drifts (v = pc²/E, E total energy; dV/dt read as the fluid derivative), the phase-space divergence is

$$
\nabla\cdot\dot{\mathbf x}+\frac{\partial\dot\mu}{\partial\mu}+\frac{1}{p^2}\frac{\partial(p^2\dot p)}{\partial p}=-\frac{\mu v}{c^2}\,\hat{\mathbf b}\cdot\left[\frac{\partial\mathbf V}{\partial t}+(\mathbf V\cdot\nabla)\mathbf V\right]
\tag{13}
$$

Hence (i) conservative and advective printings coincide for sets without an acceleration term (Zhang et al. 2023; He & Wan 2015; Laitinen et al. 2016), and otherwise differ by [Eq. (13)]·f, as Wijsen et al. (2019, advective) and Afanasiev et al. (2025, conservative) do; (ii) Eq. (12) for j = p²f is p² times Eq. (16); (iii) Eq. (10) has ∂/∂p[(dp/dt)f] without the p² Jacobian, so it agrees with the other forms if f is per unit p and differs by 2(ṗ/p)f if f is per unit momentum volume; the primer does not state the normalization in the portion read.

**Frame transformations in place of coefficients.** The Dröge-group code traces particles "in two reference frames" (solar wind, with static scattering centres; Sun-fixed, for advection and focusing) with Lorentz transformations and prints no FTE, attributing it to Isenberg (1997) [P7, P30]. Zuo et al. (2013) put p, μ and E "in the local solar wind plasma frame" and x "in a fixed frame at the shock" [P31].

### 2.4 Lagrangian co-moving forms

**EPREM.** Young & Poduval (2025), Eq. A1; Full text ×2 (arXiv HTML, PDF) [P17]:

$$
\frac{df_s}{dt}+v\mu\hat{b}\cdot\nabla f_s+\frac{(1-\mu^2)}{2}\left[-v\hat{b}\cdot\nabla\ln B-\frac{2}{v}\hat{b}\cdot\frac{d\vec{V}}{dt}+\mu\frac{d\ln(n^2/B^3)}{dt}\right]\frac{\partial f_s}{\partial\mu}+\left[-\frac{\mu}{v}\hat{b}\cdot\frac{d\vec{V}}{dt}+\mu^2\frac{d\ln(n/B)}{dt}+\frac{(1-\mu^2)}{2}\frac{d\ln B}{dt}\right]\frac{\partial f_s}{\partial\ln p}=\frac{\partial}{\partial\mu}\left(\frac{D_{\mu\mu}}{2}\frac{\partial f_s}{\partial\mu}\right)
\tag{14}
$$

The terms are labelled convection, streaming, adiabatic focusing, cooling and pitch-angle scattering; d/dt is the convective (Lagrangian) derivative in "a frame co-moving with the solar wind plasma"; n is the plasma density; the form is "based on the approach described in (Kóta et al. 2005)". Linker et al. (2019, ×1) and Kozarev et al. (2013, Eq. 1, ×2) multiply df/dt by (1 − V·e_b vμ/c²), which Eq. A1 lacks; Kozarev et al. (Eq. 1, two renderings in a later check) also print +v ê_b·∇ln B in the focusing bracket, where Eq. A1 has −v b̂·∇ln B (D-123); Linker et al. write −v e_b∇ln B without the dot, as printed [P32, P33]. Schwadron et al. (2010, ¶27; partly garbled) has the scattering term ∂/∂μ[(D_μμ/2)∂f/∂μ] [P34]. The foundation (Kóta et al. 2005) uses the frame "co-moving with the solar wind plasma"; its scattering factor is illegible [P35]. Thus D_μμ(EPREM) = 2 × the coefficient of Eq. (21); with Kozarev et al. (2013, Eq. 2) this returns 2λ∥ (MFP v1.2 Section 6.6, MFP v1.2 D-2 and U-8). EPREM is a grid code on co-moving nodes.

*Derived (check F-FTE-02).* With the ideal-MHD relations d ln n/dt = −∇·V and d ln B/dt = b̂b̂:∇V − ∇·V (verified against the Cauchy solution of a frozen-in field), μ d ln(n²/B³)/dt = μ(∇·V − 3b̂b̂:∇V) and μ² d ln(n/B)/dt + (1−μ²)/2·d ln B/dt = (1−3μ²)/2·b̂b̂:∇V − (1−μ²)/2·∇·V. The ∂f/∂μ and ∂f/∂ln p coefficients of Eq. (14) are then those of Eq. (8). The EPREM papers read do not state these relations.

**AMPS.** Tenishev et al. (2022), Eq. 1; Full text ×1 (arXiv PDF, two fetches) [P10]:

$$
\frac{\partial f}{\partial t}+\mu v\frac{\partial f}{\partial z}+\mathbf{V}^{sw}\cdot\nabla f+\frac{dp}{dt}\frac{\partial f}{\partial p}+\frac{d\mu}{dt}\frac{\partial f}{\partial\mu}-\frac{\partial}{\partial\mu}\left(D_{\mu\mu}\frac{\partial f}{\partial\mu}\right)=Q(\mathbf{x},\mathbf{p},t)
\tag{15}
$$

"The SEP transport can be conveniently described by the kinetic equation in Lagrangian coordinates", with field-line points moving "with the local solar wind velocity"; the frames of μ and p are not stated. The deterministic parts of the increments (Eqs. 2–4, "e.g., He et al., 2011") are dz = μv dt, dp = −p[(1−μ²)/2·(∂V^sw_x/∂x + ∂V^sw_y/∂y) + μ²∂V^sw_z/∂z]dt and the μ drift (1−μ²)/2·[v/L + μ(∂V^sw_x/∂x + ∂V^sw_y/∂y − 2∂V^sw_z/∂z)]dt with L(z) = B(z)/(−∂B/∂z): the coefficients of Eq. (6) with the focusing sign of Eq. (1). The stochastic part of Eq. 4 (no ∂D_μμ/∂μ drift; dW_μ = dt·Λ as extracted) is treated in Section 16. Not stated: whether the Lagrangian field-line motion replaces the V^sw·∇f term (Eq. 2 advances z by μv dt only), and along which axes ∂V_x/∂x, ∂V_y/∂y, ∂V_z/∂z are taken.

**M-FLAMPA kinetic description.** Borovikov et al. (2019, preprint) write f "in the co-moving frame of reference, which moves with the local plasma velocity". Scattering is written per Alfvén-wave branch σ in its own frame: (δf/δt)_scat = Σ_σ ∂/∂μ_σ(D^σ_μμ ∂f/∂μ_σ) (Eq. 4.29; Section 8) [P36].

### 2.5 The GCR equation of van den Berg, Engelbrecht & Strauss (2026)

Eqs. 1–4: Full text ×2 (A&A HTML alt-text LaTeX and A&A PDF); Eqs. 11, 36, 37: Full text ×1 (raw PDF text layer) [P14].

$$
\frac{\partial f}{\partial t}+\nabla\cdot(\vec v_{\rm gc}f)+\frac{1}{p^2}\frac{\partial}{\partial p}\left[p^2\frac{dp}{dt}f\right]+\frac{\partial}{\partial\mu}\left[\frac{d\mu}{dt}f\right]=\frac{\partial}{\partial\mu}\left[D_{\mu\mu}\frac{\partial f}{\partial\mu}\right]+\nabla\cdot(\mathsf D_\perp\cdot\nabla f)
\tag{16}
$$

$$
\vec v_{\rm gc}=\mu v\,\hat b+\vec u+\vec v_{\rm d}(\mu),\quad
\frac{dp}{dt}=p\left\{\frac{1-3\mu^2}{2}\hat b\hat b:\nabla\vec u-\frac{1-\mu^2}{2}\nabla\cdot\vec u-\frac{\mu}{v}\hat b\cdot\left[\frac{\partial\vec u}{\partial t}+(\vec u\cdot\nabla)\vec u\right]\right\},\quad
\frac{d\mu}{dt}=\frac{1-\mu^2}{2}\left\{v\nabla\cdot\hat b+\mu\nabla\cdot\vec u-3\mu\hat b\hat b:\nabla\vec u-\frac{2}{v}\hat b\cdot\left[\frac{\partial\vec u}{\partial t}+(\vec u\cdot\nabla)\vec u\right]\right\}
\tag{17}
$$

(Eqs. 2a–2c.) f(r, p, μ, t) is the gyrotropic number density per unit momentum volume; p and μ are "measured in a non-relativistic flow frame" (the frame of r is not stated); v is the speed in that frame, u the flow velocity, v_d(μ) the gyrophase-averaged drift; D⊥ = D⊥(μ)(I − b̂b̂), mixed coefficients neglected; L⁻¹ = ∇·b̂ ("so-called focusing length"), and near isotropy requires λ∥ ≪ |L|. Eqs. 3a–3e: F₀ = ½∫f dμ; V_d = ½∫v_d dμ; ½∫(dp/dt)dμ = −(p/3)∇·u (check F-FTE-05); κ∥ = (v²/8)∫(1−μ²)²/D_μμ dμ (Eq. (23)); κ⊥ = ½∫D⊥ dμ. Eq. 4 is the Parker equation with κ = κ∥b̂b̂ + κ⊥(I − b̂b̂). Eqs. 11 and 37:

$$
L'^{-1}=\nabla\cdot\hat b-\frac{2}{v^2}\hat b\cdot\left[\frac{\partial\vec u}{\partial t}+(\vec u\cdot\nabla)\vec u\right],\qquad
\mathcal{L}=-\frac{\mu v}{c^2}\hat b\cdot\left[\frac{\partial\vec u}{\partial t}+(\vec u\cdot\nabla)\vec u\right]
\tag{18}
$$

The PDF text layer renders the symbol of Eq. 37 as "Ł"; it is written $\mathcal L$ here. The extra term in L′ is second order in v_sw/v; L′ enters the analytic pitch-angle distributions (Eqs. 45–48; Section 12). Eq. 36 is the time-backward form in spherical coordinates; its legible μ part is (∂D_μμ/∂μ − dμ/dt)∂f/∂μ + ½·2D_μμ∂²f/∂μ², and its first term contains (μv/c²)b̂·[∂u/∂t + (u·∇)u]f.

*Derived (check F-FTE-03).* Eq. (13) equals $\mathcal L$, so the advective form of Eq. (16) reads ∂f/∂t = … − $\mathcal L$f. By the Feynman–Kac formula, a time-backward SDE for this equation carries the weight exp(−∫$\mathcal L$dτ) (derived). Their Eq. 15 as read (one rendering) has −$\mathcal L$P̃ with −∂P̃/∂τ on the left, and Eq. 18 as read has d ln w = +$\mathcal L$dτ; the sign relation between Eqs. 15, 18 and 37 is therefore not established (D-125). The relayed raw text of Eq. 36 shows this term with the opposite sign (single reading of a garbled layer; not resolved). Integration: time-backward Itô SDEs, "mathematically equivalent to a time-forward approach in the stationary limit" (Sections 15, 16).

Other pitch-angle-resolved GCR forms are Eq. (1) [P4] and Qin & Wu (2020, Eq. 2; ×2) [P18], whose streaming, convection, drift, κ⊥, cooling and focusing coefficients, with L = (b̂·∇ln B₀)⁻¹, are those of Zhang et al. (2023).

### 2.6 Forms written in a wave frame

Scattering coefficients in wave frames are covered in Sections 8 and 14; only the printed transport equations are given here.

**SOLPACS.** Afanasiev et al. (2023), Eq. 1; Full text ×2 (ar5iv, PDF); also Afanasiev et al. (2025), Eq. 12 [P37]:

$$
\frac{\partial f}{\partial t}+\left[v\mu+(1-M_A)V_A\right]\frac{\partial f}{\partial x}=\frac{\partial}{\partial\mu}\left(D_{\mu\mu}\frac{\partial f}{\partial\mu}\right)
\tag{19}
$$

x is the distance along the field from the shock (x > 0 upstream); v and μ are "as measured in the rest frame of Alfvén waves propagating at speed V_A wrt. the ambient plasma"; M_A = u₁/V_A in the de Hoffmann–Teller frame (Eq. 3). There is no focusing and no energy change.

**Ng & Reames (2008)** [P38]. NASA preprint; the raw strings were checked in two fetches, but the published PDF is garbled. Eq. 1 is a conservative equation for F = (B₀/B)(P/P₀)³f with a 2×2 diffusion matrix of D_μμ, P⁻¹D_μP and P⁻²D_PP. The layout of the matrix is not legible. The legible Eqs. 2–4 are:

$$
\dot x=\chi(\mu v+W_+)-V_{sh},\qquad
\dot\mu=(1-\mu^2)\left[r^{-1}\chi(v+\mu W_+)-v^{-1}\frac{dW_+}{dt}\right],\qquad
\dot\varepsilon=\frac{\dot P}{P}=-r^{-1}\chi W_+(1-\mu^2)-\mu v^{-1}\frac{dW_+}{dt}
\tag{20}
$$

Here χ = (1 + μvW₊/c²)⁻¹; W₊ = V_sw + V_A "relative to the Sun"; P, v, μ are in the local outward-wave frame, x is the distance from the shock and t is in the Sun's inertial frame; dW₊/dt ≡ ∂_tW₊ + ẋ∂_xW₊; B = B₀(r₀/r)², r = x + r_sh(t). Inward modes enter through V_σ = −2V_A (Sections 8, 14).

*Derived (check F-FTE-07).* With χ = 1 and ∇·b̂ = 2/r, Eq. (20) equals Eq. (17) with u → W₊r̂. The term −μ(1−μ²)∂_rW₊ of Eq. (17) is supplied in Eq. (20) by the path derivative dW₊/dt = ∂_tW₊ + W₊∂_rW₊ + μv∂_rW₊.

**Implementation requirement.** The time derivative of a frame velocity uses its source's definition: the fluid derivative ∂/∂t + V·∇ (van den Berg et al. 2020, 2026), the derivative along the particle path (Ng & Reames 2008), or not defined in the portion read (Wijsen et al. 2019; Afanasiev et al. 2025). Mixing them changes dμ/dt by μ(1−μ²)∂_∥W.

### 2.7 Summary

Abbreviations: F focusing; C energy change; V convection; D drifts; P perpendicular diffusion; op = Eq. (21); n.s. = not stated.

| Source, equation | x frame | μ, p frame | Scattering term | Terms | Verification |
|---|---|---|---|---|---|
| Effenberger & Litvinenko 2014, Eq. 1; Strauss et al. 2022, Eq. 2; Agüeda et al. 2014, Eq. 2 | field line, static | — | op | F | ×2; ×2; ×1 |
| Qin et al. 2005, Eqs. 1, 7 | field line | "solar wind frame" | op (legible) | F | ×1, garbled |
| van den Berg et al. 2020, Eq. 9; Lang et al. 2024; Merten & Aerdker 2024 | field line | n.s. | op | F (flux form) | ×2; ×1; ×1 |
| Strauss & Fichtner 2015; Strauss et al. 2017 | n.s. | n.s. | op | F, P | ×1–2 |
| Ruffolo 1995, Eqs. 7–11 | fixed | solar wind | ∂_μ[(φ/2)∂_μ] | F, C, V | ×1; Eq. 11 assembled |
| van den Berg et al. 2020, Eqs. 17–18; Strauss et al. 2020 | observer | solar wind | op + mixed D⊥ | F, C, V, D, P | ×2/×1; garbled |
| Zhang et al. 2023, Eqs. 1–4 | corotating (PDF only) | μ to outward field; p n.s. | op | F, C, V, D, P | ×2 (Eq. 1), ×1 |
| He & Wan 2015, Eqs. 1–3 (He et al. 2011, 2017; He & Wan 2019) | n.s. | n.s. | op | F (sign flagged), C, V, P | ×1 |
| Qin & Wu 2020, Eq. 2; Wang et al. 2014, Eq. 1 | n.s. | n.s.; μ "plasma frame" (Wang) | op | F, C, V, P; D (Qin & Wu) | ×2; ×1 |
| Kubo et al. 2015, Eqs. 1–3 | co-rotating | p solar wind; μ n.s. | op | F, C, V | ×1 |
| Laitinen et al. 2016, Eq. 1 | n.s. | n.s. | op | F, C, V, P | ×2 |
| Wijsen et al. 2019, Eqs. 1–5 | fixed inertial | solar wind | op | F, C, V, D, P | ×2 (Eq. 1), ×1 |
| Afanasiev et al. 2025, Eqs. 37–40 | n.s. | n.s. | op on j = p²f | F, C, V | ×1 |
| van den Berg et al. 2026, Eqs. 1–2 | n.s. | flow frame | op | F, C, V, D, P | ×2 |
| EPREM: Young & Poduval 2025; Kozarev et al. 2013; Linker et al. 2019 | co-moving nodes | co-moving | ∂_μ[(D_μμ/2)∂_μ] | F, C, V (in d/dt) | ×2; ×2; ×1 |
| AMPS: Tenishev et al. 2022, Eq. 1 | Lagrangian field line | n.s. | op | F, C, V | ×1 |
| SOLPACS: Afanasiev et al. 2023, Eq. 1 | shock | Alfvén-wave frame | op | V | ×2 |
| Ng & Reames 2008, Eqs. 1–4 | from shock; t Sun frame | outward-wave frame | matrix (not legible) | F, C, V | ×1 (preprint) |

**Direction of integration (as stated).**
- **Time-forward particle codes:** Effenberger & Litvinenko 2014 (Milstein); Qin et al. 2005; Kartavykh et al. 2015; PARADISE (Wijsen et al. 2019; Husidic et al. 2024, "integrated forward in time"); Zhang & Zhao 2017 [P39]; Zuo et al. 2013; Merten & Aerdker 2024; Tenishev et al. 2022 (increments per simulation time step).
- **Time-backward SDEs:** Zhang et al. 2009 (abstract: "in a fixed reference frame in which any spacecraft is roughly stationary") [P40]; Qin et al. 2006 (abstract) [P41]; He et al. 2011, 2017; He & Wan 2015, 2019; Wang et al. 2014; Qin & Wu 2020; Kubo et al. 2015; Zhang et al. 2023; Ding et al. 2022 [P42]; van den Berg et al. 2026.
- **SDEs, direction not stated:** Laitinen et al. 2016; Dröge et al. 2010 (abstract) [P43].
- **Grid solvers:** Ruffolo 1995; Strauss & Fichtner 2015; Strauss et al. 2017, 2022; the primer's main text; EPREM; Ng & Reames 2008.

Strauss & Effenberger (2017) note that "the time backward integration approach is mostly followed", while SEP modelling has used "the time forward SDE formalism" [P44].

**Implementation requirement.** Every FTE coefficient set records (1) the printed frame of x, μ and p, or "not stated"; (2) whether μ is relative to B or to the outward direction [P15]; (3) the operator factor (1, or ½ for Ruffolo's φ and for EPREM); (4) the density advanced (f, Af, p²f or F); (5) the definition of d/dt in frame-velocity terms. No set is converted silently; conversions use Eq. (13) or the relations stated in this section.

## 3. Core definitions and conventions

This section fixes the operator, the definition of D_μμ, the relation to λ∥ and κ∥, isotropic scattering and the Itô stochastic differential equation. Every other section refers to these five equations. The relations are standard results; their derivations are in the companion documents where cited, and the consequences stated here are labelled derived.

### 3.1 The pitch-angle operator

For a gyrotropic distribution f(**x**, p, μ, t), pitch-angle scattering enters the focused transport equation (Section 2) as

$$
\left.\frac{\partial f}{\partial t}\right|_{\rm scat}=\frac{\partial}{\partial\mu}\left(D_{\mu\mu}\,\frac{\partial f}{\partial\mu}\right),\qquad -1\le\mu\le1 .
\tag{21}
$$

Most sources print this term on the right-hand side, or on the left-hand side with a minus sign, which is the same operator (Section 2). D_μμ ≥ 0 has the unit s⁻¹. The flux in μ is −D_μμ∂f/∂μ; with zero flux at μ = ±1, ∫₋₁¹f dμ is conserved, an isotropic f is stationary, and ∫₋₁¹f²dμ does not increase, since its rate of change is −2∫₋₁¹D_μμ(∂f/∂μ)²dμ ≤ 0 (derived). Shapes that contain the factor (1 − μ²) give zero flux at μ = ±1 for any f with bounded ∂f/∂μ; printed variants whose coefficient does not vanish at μ = ±1 need the condition imposed explicitly (Section 5.5; MFP v1.2, Section 21.5). A coefficient printed with a sign that becomes negative for part of the μ range is not a diffusion coefficient there (for example Ding et al. [P42] with vμ in place of v|μ|; D-36, U-19).

### 3.2 Definition through the pitch-angle increments

The coefficient of Eq. (21) is

$$
D_{\mu\mu}=\lim_{\Delta t\to0}\frac{\langle(\Delta\mu)^2\rangle}{2\,\Delta t},
\tag{22}
$$

where Δμ is the change of μ in a time Δt that is short compared with the time scale on which D_μμ changes μ by order one, and the average is over particles with the same initial μ. For a Markov process with this second moment and the mean increment lim ⟨Δμ⟩/Δt = ∂D_μμ/∂μ, the forward Kolmogorov (Fokker–Planck) equation ∂P/∂t = −∂_μ(D′_μμP) + ∂²_μ(D_μμP) equals Eq. (21) term by term (derived; the product rule gives ∂²_μ(DP) − ∂_μ(D′P) = ∂_μ(D∂_μP)). The mean increment is therefore part of the definition: a process with the second moment of Eq. (22) and zero mean increment realizes a different operator.

Some papers define the coefficient through ⟨(Δμ)²⟩/Δt without the factor 2; a value from such a paper is twice Eq. (22) (Riordan & Pe'er [P45], Section 11.1; D-80, U-41). Estimators applied to simulations at finite Δt are biased (Section 11).

### 3.3 Operator conventions in the literature

| Convention | Printed scattering term | Relation to D_μμ of Eq. (21) | Sources |
|---|---|---|---|
| Standard | ∂_μ(D_μμ∂_μf) | identical | most sources (Sections 2, 5, 6) |
| Ruffolo | ∂_μ[(φ/2)∂_μF] | φ = 2D_μμ | [P6] (Sections 2.3, 16.1) |
| EPREM | ∂_μ[(D^E_μμ/2)∂_μf] | D^E_μμ = 2D_μμ | [P33, P17, P32] (Sections 2.4, 5.4; MFP v1.2, Section 6.6 and MFP v1.2 D-2) |
| Spherical form of isotropic scattering | D∂_μ[(1 − μ²)∂_μf] + D/(1 − μ²)∂²f/∂φ² | D_μμ = D(1 − μ²), ν = 2D | [P46] (Section 16.7) |
| Increment definition without the factor 2 | ⟨(Δμ)²⟩/Δt | twice D_μμ | [P45] (Section 11.1) |

The EPREM isotropic coefficient is printed as D^E_μμ = (1 − μ²)v/(2λ∥); the coefficient of Eq. (21) is then (1 − μ²)v/(4λ∥), for which Eq. (23) returns 2λ∥ (derived; Section 5.4; MFP v1.2 D-2). In a pusher, the drift and the noise both follow the coefficient of ∂_μ(·∂_μf) (Eq. (151)).

**Implementation requirement.** Every D_μμ model stores the operator convention of its source. The library converts to the coefficient of Eq. (21) once, at the model boundary, and never changes a printed coefficient by a factor of two silently.

### 3.4 The parallel mean free path

In the diffusion limit of the focused transport equation without focusing, the parallel diffusion coefficient and mean free path are (PARALLEL rev. 1.4, Section 7.1, Eq. 20; MFP v1.2, Section 6.1, Eq. 11)

$$
\lambda_\parallel=\frac{3v}{8}\int_{-1}^{1}\frac{(1-\mu^2)^2}{D_{\mu\mu}(\mu)}\,d\mu,\qquad
\kappa_\parallel=\frac{v\lambda_\parallel}{3}=\frac{v^2}{8}\int_{-1}^{1}\frac{(1-\mu^2)^2}{D_{\mu\mu}(\mu)}\,d\mu .
\tag{23}
$$

D_μμ here is the coefficient of Eq. (21). The relation does not assume that D_μμ is even in μ; for an even D_μμ the integral is twice the integral over 0 < μ < 1, and this symmetry is not applied to directional or polarized wave models without checking it (PARALLEL rev. 1.4, Section 7.1). For a shape D_μμ = D₀g(μ) the amplitude follows from λ∥ exactly: D₀ = (3v/(8λ∥))∫₋₁¹(1 − μ²)²/g dμ; the closed forms for the q-form and the ε-form are in Section 5.6.

Consequences (derived):

- If D_μμ ∝ |μ|^α near μ = 0, the integral converges for α < 1 and diverges for α ≥ 1. The q-form with H = 0 and 1 < q < 2 (α = q − 1) gives a finite λ∥; the ε-form with ε = 0 (α = 1) gives a logarithmically divergent λ∥; a dissipation-range spectrum k^{−p} gives D_μμ ∝ |μ|^{p−1} near μ = 0 in magnetostatic QLT and a divergent λ∥ for p ≥ 2 (Section 7.7).
- A cut-off at |μ| = μ_min (Chen et al. [P47]) defines a different quantity; it is a normalization device, not a scattering model (Sections 13.2 and 16.5).
- With focusing, the spatial spreading that a pusher produces differs from Eq. (23); the published relations are in Section 12, and none of them is applied as a default (U-46). The published relations are closure coefficients; the long-time centred variance of the constant-L isotropic process is Eq. (115), which differs from them (Section 12.3).

### 3.5 Isotropic scattering

Isotropic scattering with frequency ν is

$$
D_{\mu\mu}=\frac{\nu}{2}\left(1-\mu^2\right),\qquad \lambda_\parallel=\frac{v}{\nu},\qquad \kappa_\parallel=\frac{v^2}{3\nu}.
\tag{24}
$$

The sources write the same model as D₀(1 − μ²) with λ₀ = v/(2D₀) [P1, P2], as D_I(1 − μ²) with D_I = v²/(6κ∥) [P14], and as a constant scattering frequency ν₀ with λ∥ = v/ν₀ [P48] (Sections 5.1 and 15.7). Properties used for validation (derived; Section 16.9):

- The Legendre polynomials P_l(μ) are the eigenfunctions of Eq. (21) with eigenvalues −νl(l + 1)/2.
- Hence ⟨μ⟩(t) = ⟨μ⟩₀e^{−νt} and ⟨μ²⟩(t) = 1/3 + (⟨μ²⟩₀ − 1/3)e^{−3νt} (Eq. (174)).
- From an isotropic start, ⟨μ(0)μ(t)⟩ = e^{−νt}/3, so that the running spatial coefficient tends to κ∥ = v²/(3ν) (Eq. (175)).
- Poisson-timed isotropic redraws with rate r have the same ⟨μ⟩(t) and κ∥ for r = ν, but ⟨μ²⟩ relaxes at the rate r instead of 3ν; they do not realize Eq. (21) (Eq. (171)).

### 3.6 The Itô stochastic differential equation

A pusher realizes Eq. (21), together with the deterministic pitch-angle rate dμ/dt of the transport equation (focusing, adiabatic and frame terms; Section 2), by the Itô stochastic differential equation

$$
d\mu=\left(\frac{d\mu}{dt}+\frac{\partial D_{\mu\mu}}{\partial\mu}\right)dt+\sqrt{2D_{\mu\mu}}\;dW ,
\tag{25}
$$

where W is a standard Wiener process (dW has mean zero and variance dt) and all coefficients are evaluated at the start of the step (Itô calculus). For dμ/dt = 0 this is MFP v1.2 Eq. 46. Its properties, stated in Section 16:

- In forward time, the density of pseudo-particles obeys Eq. (150), in which the deterministic term is in conservative form ∂_μ[(dμ/dt)P]. A forward SDE therefore samples the density for which the transport equation is written in conservative form, for example the flux-tube density Af in one dimension, and the configuration states which density is sampled (Section 16.1; U-2).
- In backward time the deterministic rates change sign while ∂D_μμ/∂μ and the noise do not (Eq. (152)); the backward SDE samples f itself (Section 16.1).
- Codes that write the operator with φ/2 or D^E_μμ/2 convert both the drift and the noise (Eq. (151)).
- The exact process does not reach μ = ±1 when D_μμ vanishes linearly there, but discrete steps do; a boundary rule or a kernel that stays in [−1, 1] is needed (Section 16.4).
- For the q-form and the |μ|^{2/3} shape the drift ∂D_μμ/∂μ is unbounded at μ = 0 for every H ≥ 0, and for the ε-form it is discontinuous there (Eqs. (163) and (164); Section 16.5).

**Implementation requirement.** A pitch-angle step that claims to realize Eq. (21) reproduces the first two moments of Eq. (25) to first order in Δt, including the drift ∂D_μμ/∂μ, and its convergence with the step size is demonstrated against the isotropic results of Section 3.5 (Section 16.11).

## 4. Model inventory and runtime states

The inventory lists what the literature provides. Whether a configured run can execute is a separate question, answered by the runtime state and the admissibility fields of Section 4.1. The companion file `model_registry.json` holds, for each of the 172 model IDs, its section, equation labels, source keys, specification status, runtime state, what the model provides (`provides`), its physical admissibility (`admissibility`), its stated domain where the source restricts it (`domain`), and a note naming the open D-/U-items. Not every ID is a pitch-angle coefficient: the field `provides` separates D_μμ laws from transport-equation coefficient sets, spectra and spectrum adapters, resonance or velocity-space kernels, growth rates and wave providers, coupling schemes, focusing relations and spatial coefficients, estimators and benchmarks, pusher kernels and validation references, so that only a D_μμ law is invoked where a pitch-angle coefficient is required.

### 4.1 Runtime states

A configured run is in exactly one of these states. The registry records the state that applies when the inputs named in the source are supplied; the state of an actual run is evaluated per configuration.

| Runtime state | Definition | Examples |
|---|---|---|
| `READY_EXPLICIT_INPUTS` | The law is specified, and every input, convention and unit it needs is supplied by the configuration or the source | PA-ISO, PA-QFORM and PA-EPS with their parameters and an amplitude mode; QLT-SLAB with a spectrum whose normalization is declared; SDE-EM, SDE-MIL and MC-ROT; GCR-PA-VDB26-ANISO |
| `READY_PUBLISHED_VARIANT` | Reproducible only as an explicitly named printed variant, for the reproduction of the printed source or as a diagnostic; never a default, and not a physical model where its admissibility is "unestablished" | PA-QFORM-LANG-PRINTED (D-32), PA-EPS-PACHECO-PRINTED (D-33), FTE-SW-HEWAN15-PRINTED (D-1), WAVE-SL19-ND-PRINTED (D-48), GCR-PA-QINWU20-PRINTED (D-13), FOC-SH11-QUOTED |
| `REQUIRES_EXTERNAL_INPUT` | The formula is known, but a value it needs is not printed and must be supplied | WAVE-TWOFRAME-VV05 (weights ω_±, U-29); SGW-COUPLE-CSA (particle time step) |
| `REQUIRES_USER_DECISION` | The formula and inputs exist, but a decision listed in Section 19.2 is open | MEAS-CHEN24 (U-23, U-47), FTE-SW-HEWAN15 (U-1), PA-EPS-ETA (U-13), SDE-IS-ZHANG23 (U-58), MC-REDRAW (U-63), GCR-PA-VDB26-BSDE (U-55) |
| `REQUIRES_SOURCE_OR_CODE_AUDIT` | The equation, or what a host code implements, is not established; an open audit item blocks the model before any user decision on it is taken | PA-AMPS-I…V (garbled table: D-28 to D-31), SDE-AMPS-PRINTED (D-10, U-65), FTE-SW-ZHANG23 (equations read only in a summarized rendering; frames not stated), MEAS-ZHANG23 (D-41), DYN-NADT-SH06 (D-46), SDE-IS-HEWAN15-PRINTED (Section 16.2), WAVE-SPA22-BW (units of the prefactor, D-127), SGW-COUPLE-VK01 (backward-wave sign, D-126) |
| `REFERENCE_DATA_ONLY` | Published simulation values, a printed step or a printed form kept for comparison or demonstration, not a predictive closure | SIM-VDB24-C1…C3, SIM-RIORDAN19 (D-80), SIM-XUYAN13, SDE-MIL-EL14-PRINTED (Section 16.10), DYN-STRAUSS17-PRINTED (dimensionally inconsistent as printed, D-44) |

A model that needs an equation of a companion specification (Section 1.2) whose pinned text is not available to the implementer is not assigned one of these states: the implementation reports it as `BLOCKED_DEPENDENCY`. An unknown or incomplete configuration fails with a structured error that names the missing parameter, the source item or the D-/U-identifier.

**Admissibility (version 1.1).** The runtime state records whether the inputs of a model are specified. Whether a printed form is a physical coefficient at all is a separate question, recorded in the field `admissibility` of the registry with the four values "admissible", "admissible on a stated domain or regime only", "inadmissible as printed" and "unestablished". A printed form is inadmissible when it is dimensionally inconsistent (DYN-STRAUSS17-PRINTED, D-44; the prefactor of WAVE-SPA22-BW as printed, D-127), when it is negative or not real on part of the pitch-angle interval it is applied to (the printed μ^{s−1} of GCR-PA-QINWU20-PRINTED for μ < 0, D-13; the printed μ^{11/2} of Eq. (77) for μ < 0, Section 9.1), or when its sign contradicts the declared definition of the focusing length (the second printed form of He & Wan 2015 Eq. 3, D-1). It is unestablished when the printed form is consistent but the source gives no evidence that the printed behaviour is intended (the endpoint behaviour of the (1−μ)² and |μ|/(1−|μ|) variants, D-32, D-33). An inadmissible form is never evaluated as a coefficient of Eq. (21), whatever its runtime state and whatever the user configures; it stays in the registry so that the discrepancy can be cited and, where useful, demonstrated. An unestablished form is evaluated only under its own ID for the reproduction of the printed source or as a diagnostic. A model whose note names an open audit item (D-126, D-127) is in `REQUIRES_SOURCE_OR_CODE_AUDIT` even where a user decision on another aspect of it is also open: the audit blocks first. A configured runtime is executable only when the runtime state, the admissibility, the domain of the coefficient (field `domain`; Sections 9.1 and 16.5), the requested output (field `provides`) and the coefficient–kernel compatibility of Section 16.11 all pass.

Counts reported for the intended version 1.1 `model_registry.json` (matching archive not supplied to this review; Section 25.7):

| Runtime state | Model IDs |
|---|---|
| `READY_EXPLICIT_INPUTS` | 86 |
| `READY_PUBLISHED_VARIANT` | 6 |
| `REQUIRES_EXTERNAL_INPUT` | 2 |
| `REQUIRES_USER_DECISION` | 31 |
| `REQUIRES_SOURCE_OR_CODE_AUDIT` | 35 |
| `REFERENCE_DATA_ONLY` | 12 |
| total | 172 |

By admissibility (Section 4.1): admissible: 154; admissible on a stated domain or regime only: 7; inadmissible as printed: 6; unestablished: 5.

By what the ID provides (field `provides`): 73 IDs are D_μμ laws or pitch-angle tensors (including presets, printed variants and records); the other 99 are transport-equation coefficient sets, spectra and spectrum adapters, frame transformations, resonance or velocity-space kernels, spatial coefficients, focusing relations, growth rates and wave providers, coupling schemes, simulation benchmarks and estimators, pusher kernels, step and boundary rules, and validation references. Capability dispatch is typed: a scalar D_μμ law and a full (μ,p) tensor use separate interfaces (Section 17). Neither a string prefix nor a non-scattering inventory ID qualifies a provider for a kernel.

### 4.2 Inventory by section

Each model ID is listed with the subsection that defines it, the equations of this document that it implements (— where the subsection states the model in words or by reference to a source equation) and its runtime state.

**Section 2 — Focused transport equation coefficient sets**

| Model ID | Section | Equations | Runtime state |
|---|---|---|---|
| FTE-1D-EL14 | 2.1 | (1) | `READY_EXPLICIT_INPUTS` |
| FTE-1D-TUBE | 2.1 | (2) | `REQUIRES_USER_DECISION` |
| FTE-SW-ZHANG23 | 2.2 | (3), (4) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| FTE-SW-HEWAN15 | 2.2 | (5), (6) | `REQUIRES_USER_DECISION` |
| FTE-SW-HEWAN15-PRINTED | 2.2 | (5), (6) | `READY_PUBLISHED_VARIANT` |
| FTE-MIX-WIJSEN19 | 2.2 | (7), (8) | `REQUIRES_USER_DECISION` |
| FTE-MIX-RUFFOLO95 | 2.3 | (9) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| FTE-MIX-VDB20 | 2.3 | (10), (11) | `REQUIRES_USER_DECISION` |
| FTE-CONS-PARADISE25 | 2.3 | (12) | `REQUIRES_USER_DECISION` |
| FTE-DIV | 2.3 | (13) | `READY_EXPLICIT_INPUTS` |
| FTE-LAG-EPREM | 2.4 | (14) | `READY_EXPLICIT_INPUTS` |
| FTE-LAG-AMPS22 | 2.4 | (15) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| FTE-GCR-VDB26 | 2.5 | (16), (17), (18) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| FTE-WAVE-SOLPACS23 | 2.6 | (19) | `READY_EXPLICIT_INPUTS` |
| FTE-WAVE-NG08 | 2.6 | (20) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |

**Section 5 — Empirical and phenomenological shapes**

| Model ID | Section | Equations | Runtime state |
|---|---|---|---|
| PA-ISO | 5.1 | (24) | `READY_EXPLICIT_INPUTS` |
| PA-QFORM | 5.2 | (26), (27), (37), (38) | `READY_EXPLICIT_INPUTS` |
| PA-QFORM-HE11 | 5.2 | (26) | `REQUIRES_USER_DECISION` |
| PA-QFORM-WQ15 | 5.2 | (26) | `REQUIRES_USER_DECISION` |
| PA-QFORM-Z23 | 5.2 | (26) | `REQUIRES_USER_DECISION` |
| PA-QFORM-W19 | 5.2 | (26) | `REQUIRES_USER_DECISION` |
| PA-QFORM-QW15 | 5.2 | (28) | `REQUIRES_USER_DECISION` |
| PA-EPS | 5.3 | (29), (39) | `READY_EXPLICIT_INPUTS` |
| PA-EPS-W19 | 5.3 | (30) | `REQUIRES_USER_DECISION` |
| PA-EPS-ETA | 5.3 | (31) | `REQUIRES_USER_DECISION` |
| PA-DROGE-VA | 5.4 | (32), (41) | `REQUIRES_USER_DECISION` |
| PA-KOLMO | 5.4 | (33) | `READY_EXPLICIT_INPUTS` |
| PA-AMPS-I | 5.4 | (33) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| PA-AMPS-II | 5.4 | — | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| PA-AMPS-III | 5.4 | (43) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| PA-AMPS-IV | 5.4 | (26) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| PA-AMPS-V | 5.4, 6.3 | (45) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| PA-EPREM | 5.4 | (34) | `READY_EXPLICIT_INPUTS` |
| PA-QFORM-LANG-PRINTED | 5.5 | (35) | `READY_PUBLISHED_VARIANT` |
| PA-EPS-PACHECO-PRINTED | 5.5 | (36), (40) | `READY_PUBLISHED_VARIANT` |

**Section 6 — Quasilinear slab D_μμ**

| Model ID | Section | Equations | Runtime state |
|---|---|---|---|
| QLT-SLAB | 6.1 | (43) | `READY_EXPLICIT_INPUTS` |
| QLT-SLAB-HEL | 6.1 | (42) | `REQUIRES_USER_DECISION` |
| QLT-SLAB-PL | 6.3 | (44), (45) | `READY_EXPLICIT_INPUTS` |
| QLT-SLAB-TS03 | 6.4 | (46) | `READY_EXPLICIT_INPUTS` |
| QLT-SLAB-TS02 | 6.4 | (46) | `READY_EXPLICIT_INPUTS` |
| QLT-SLAB-BENDOVER | 6.5 | (47) | `READY_EXPLICIT_INPUTS` |

**Section 7 — Dynamical turbulence and dissipation range**

| Model ID | Section | Equations | Runtime state |
|---|---|---|---|
| DYN-TS-DT | 7.1 | (48), (50), (51) | `READY_EXPLICIT_INPUTS` |
| DYN-TS-RS | 7.1 | (49), (50), (51) | `READY_EXPLICIT_INPUTS` |
| DYN-STRAUSS17-PRINTED | 7.4 | (52) | `REFERENCE_DATA_ONLY` |
| DYN-SLS08 (alias NL-NLD-SH08) | 7.5 | (53) | `READY_EXPLICIT_INPUTS` |
| DYN-NADT-SH06 (alias NL-NADT-SH06) | 7.6, 10.5 | — | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| DYN-KD-SLR19 | 7.7 | — | `READY_EXPLICIT_INPUTS` |

**Section 8 — Wave-based D_μμ**

| Model ID | Section | Equations | Runtime state |
|---|---|---|---|
| WAVE-SL19-DT | 8.1 | (56) | `REQUIRES_USER_DECISION` |
| WAVE-SL19-ND-PRINTED | 8.1 | (57) | `READY_PUBLISHED_VARIANT` |
| WAVE-SPLIT-SPA22 | 8.1 | (58) | `READY_EXPLICIT_INPUTS` |
| WAVE-XFORM-LORENTZ | 8.2 | (60) | `READY_EXPLICIT_INPUTS` |
| WAVE-XFORM-GAL | 8.2 | (61) | `READY_EXPLICIT_INPUTS` |
| WAVE-DE | 8.2 | (59), (62) | `READY_EXPLICIT_INPUTS` |
| WAVE-TWOFRAME-VV05 | 8.2 | — | `REQUIRES_EXTERNAL_INPUT` |
| WAVE-A2-VV05 | 8.2 | (64) | `READY_EXPLICIT_INPUTS` |
| WAVE-BOR19-BRANCH | 8.3 | (66) | `READY_EXPLICIT_INPUTS` |
| WAVE-BOR19-PF | 8.3 | (67) | `READY_EXPLICIT_INPUTS` |
| WAVE-VLF03 | 8.3 | — | `READY_EXPLICIT_INPUTS` |
| WAVE-NR08 | 8.3 | — | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| WAVE-SPA22-GEN | 8.4 | (70) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| WAVE-SPA22-BW | 8.4 | (71) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| WAVE-DISP-SPA22 | 8.4 | (68), (72) | `READY_EXPLICIT_INPUTS` |

**Section 9 — Anisotropic Alfvénic and compressible turbulence**

| Model ID | Section | Equations | Runtime state |
|---|---|---|---|
| CMP-CH00-DXX | 9.1 | (74) | `READY_EXPLICIT_INPUTS` |
| CMP-CH00-KAPPA16 | 9.1 | (75) | `READY_EXPLICIT_INPUTS` |
| CMP-CH00-KAPPA17 | 9.1 | (75) | `READY_EXPLICIT_INPUTS` |
| CMP-YL02-ALF | 9.1 | (76) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| CMP-YL04-ALF-PRINTED | 9.1 | — | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| CMP-XL20-ALF | 9.1 | (77) | `READY_EXPLICIT_INPUTS` |
| CMP-XL20-SLOW | 9.1 | — | `READY_EXPLICIT_INPUTS` |
| CMP-XL20-FAST-G | 9.2 | (78) | `READY_EXPLICIT_INPUTS` |
| CMP-YL04-FAST-HALO | 9.2 | — | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| CMP-XL18-FAST-TTD-LIN | 9.2 | (79) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| CMP-YL08-RES (alias NL-YL08-RES) | 9.3, 10.5 | (81), (101) | `REQUIRES_USER_DECISION` |
| CMP-YL08-TTD | 9.3 | (82) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| CMP-YL08-TTD-SMALL | 9.3 | (83) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| CMP-XL18-SLOW-TTD | 9.4 | (85) | `READY_EXPLICIT_INPUTS` |
| CMP-XL18-FAST-TTD-B | 9.4 | (86) | `READY_EXPLICIT_INPUTS` |
| CMP-XL20-TRAP | 9.4 | (87) | `READY_EXPLICIT_INPUTS` |
| CMP-LYNN12-REXP | 9.5 | (88) | `READY_EXPLICIT_INPUTS` |
| CMP-LYNN12-RGAUSS | 9.5 | (88) | `READY_EXPLICIT_INPUTS` |
| SIM-XUYAN13 | 9.5, 11.5 | — | `REFERENCE_DATA_ONLY` |

**Section 10 — Nonlinear and second-order theories**

| Model ID | Section | Equations | Runtime state |
|---|---|---|---|
| NL-SOQLT-SH09-GAUSS | 10.1 | (89), (91), (90) | `READY_EXPLICIT_INPUTS` |
| NL-SOQLT-SH09-STEP | 10.1 | (89), (92) | `READY_EXPLICIT_INPUTS` |
| NL-SOQLT-SH09-CF | 10.3 | (94) | `READY_EXPLICIT_INPUTS` |
| NL-SOQLT-SH09-ISO | 10.3 | (93) | `READY_EXPLICIT_INPUTS` |
| NL-SOQLT-SH09-LAM | 10.3 | (95), (96) | `REQUIRES_USER_DECISION` |
| NL-SOQLT-SUB25 | 10.4 | (98), (99) | `READY_EXPLICIT_INPUTS` |
| NL-SOQLT-SUB25-MU0 | 10.4 | (100) | `READY_EXPLICIT_INPUTS` |
| NL-WNLT-SH04 | 10.5 | — | `REQUIRES_SOURCE_OR_CODE_AUDIT` |

**Section 11 — Simulation estimators and benchmarks**

| Model ID | Section | Equations | Runtime state |
|---|---|---|---|
| SIM-M5-ISO | 11.1 | (102) | `READY_EXPLICIT_INPUTS` |
| SIM-TAUTZ13II-ISO | 11.2 | (103), (104), (105) | `READY_EXPLICIT_INPUTS` |
| SIM-M2B-EXACT | 11.3 | (106) | `READY_EXPLICIT_INPUTS` |
| SIM-TAUTZ13-TDQLT | 11.4 | (107), (108) | `READY_EXPLICIT_INPUTS` |
| SIM-VDB24-C1 | 11.5, 11.6 | — | `REFERENCE_DATA_ONLY` |
| SIM-VDB24-C2 | 11.5, 11.6 | — | `REFERENCE_DATA_ONLY` |
| SIM-VDB24-C3 | 11.5, 11.6 | — | `REFERENCE_DATA_ONLY` |
| SIM-IVAS16-MC | 11.5 | — | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| SIM-RIORDAN19 | 11.5 | — | `REFERENCE_DATA_ONLY` |
| SIM-DUNDOVIC20 | 11.5 | — | `REFERENCE_DATA_ONLY` |
| SIM-LANGE13 | 11.5 | — | `REFERENCE_DATA_ONLY` |
| SIM-WEIDL15 | 11.5 | — | `REFERENCE_DATA_ONLY` |
| SIM-TAUTZ13-T1 | 11.6 | — | `REFERENCE_DATA_ONLY` |
| SIM-TAUTZ13II-F3 | 11.6 | — | `REFERENCE_DATA_ONLY` |

**Section 12 — Focusing relations**

| Model ID | Section | Equations | Runtime state |
|---|---|---|---|
| FOC-LITV12A-WEAK | 12.2 | (110) | `READY_EXPLICIT_INPUTS` |
| FOC-SH11-QUOTED | 12.2 | (111) | `READY_PUBLISHED_VARIANT` |
| FOC-EL14 | 12.3 | (112), (114) | `READY_EXPLICIT_INPUTS` |
| FOC-WQ18-A | 12.4 | (116) | `READY_EXPLICIT_INPUTS` |
| FOC-WQ18-SERIES | 12.4 | (117) | `READY_EXPLICIT_INPUTS` |
| FOC-T12-QL | 12.5 | (118) | `REQUIRES_USER_DECISION` |
| FOC-HW13 (alias MFP-FOCUS-HW13, FOC-T12-NL) | 12.5 | (118) | `READY_EXPLICIT_INPUTS` |
| FOC-TSD14-D0 | 12.6 | (120) | `READY_EXPLICIT_INPUTS` |

**Section 13 — Measured or prescribed spectra**

| Model ID | Section | Equations | Runtime state |
|---|---|---|---|
| MEAS-TAYLOR | 13.1 | (121) | `READY_EXPLICIT_INPUTS` |
| MEAS-CHEN24 | 13.2 | (122) | `REQUIRES_USER_DECISION` |
| MEAS-WIJSEN19 | 13.3 | (123) | `REQUIRES_USER_DECISION` |
| MEAS-LAITINEN16 | 13.4 | (124) | `REQUIRES_USER_DECISION` |
| MEAS-ZHANG23 | 13.5 | — | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| MEAS-DING22 | 13.6 | (125) | `REQUIRES_USER_DECISION` |
| MEAS-AFANASIEV15 | 13.7, 8.3 | (126), (65) | `REQUIRES_USER_DECISION` |
| MEAS-SOLPACS-GRID | 13.7 | — | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| MEAS-SUBASH25-TAYLOR | 13.8 | — | `REQUIRES_SOURCE_OR_CODE_AUDIT` |

**Section 14 — Self-generated waves and shocks**

| Model ID | Section | Equations | Runtime state |
|---|---|---|---|
| SGW-GROWTH-INT | 14.3 | (128) | `READY_EXPLICIT_INPUTS` |
| SGW-GROWTH-STREAM | 14.3 | (129) | `READY_EXPLICIT_INPUTS` |
| SGW-GROWTH-VK01 | 14.3 | (130) | `READY_EXPLICIT_INPUTS` |
| SGW-GROWTH-SLR19 | 14.3 | (131) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| SGW-WAVE-CSA | 14.4 | (132) | `READY_EXPLICIT_INPUTS` |
| SGW-WAVE-SHOCK | 14.4 | (133) | `READY_EXPLICIT_INPUTS` |
| SGW-BELL-A15 | 14.5 | (134) | `READY_EXPLICIT_INPUTS` |
| SGW-COUPLE-VK01 | 14.6 | (135), (136) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| SGW-COUPLE-CSA | 14.6 | — | `REQUIRES_EXTERNAL_INPUT` |
| SGW-COUPLE-SOLPACS | 14.6 | — | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| SGW-COUPLE-SANDROOS14 | 14.6 | (137) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| SGW-COUPLE-SLR19 | 14.6 | (138) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |

**Section 15 — Pitch-angle-resolved GCR models**

| Model ID | Section | Equations | Runtime state |
|---|---|---|---|
| GCR-PA-VDB26-ISO | 15.2 | (139), (140) | `READY_EXPLICIT_INPUTS` |
| GCR-PA-VDB26-ANISO | 15.2 | (139), (140) | `READY_EXPLICIT_INPUTS` |
| GCR-PA-VDB26-DPERP | 15.2 | (141) | `READY_EXPLICIT_INPUTS` |
| GCR-PA-VDB26-BSDE | 15.2 | (142), (143), (144) | `REQUIRES_USER_DECISION` |
| GCR-PA-VDB26-PAD | 15.2 | (145) | `READY_EXPLICIT_INPUTS` |
| GCR-PA-QINWU20 | 15.3 | (146) | `REQUIRES_USER_DECISION` |
| GCR-PA-QINWU20-PRINTED | 15.3 | (146) | `READY_PUBLISHED_VARIANT` |
| GCR-PA-STRAUSS22 | 15.4 | (26) | `READY_EXPLICIT_INPUTS` |
| GCR-PA-KUBO15 | 15.5 | (147), (26) | `REQUIRES_USER_DECISION` |
| GCR-PA-FLORINSKI13 | 15.6 | (148) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| GCR-PA-ZHANG14-LISM | 15.6 | — | `READY_EXPLICIT_INPUTS` |
| GCR-PA-EL14-BENCH | 15.7 | (149) | `READY_EXPLICIT_INPUTS` |

**Section 16 — Numerical kernels and validation targets**

| Model ID | Section | Equations | Runtime state |
|---|---|---|---|
| SDE-ITO-FWD | 16.1 | (25), (150) | `READY_EXPLICIT_INPUTS` |
| SDE-ITO-BWD | 16.1 | (152) | `READY_EXPLICIT_INPUTS` |
| SDE-CONV-RUFFOLO | 16.1 | (151) | `READY_EXPLICIT_INPUTS` |
| SDE-CONV-EPREM | 16.1 | (151) | `READY_EXPLICIT_INPUTS` |
| SDE-IS-ZHANG23 | 16.2 | (153), (154) | `REQUIRES_USER_DECISION` |
| SDE-IS-HEWAN15-PRINTED | 16.2 | (155) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| SDE-W-VDB26 | 16.2 | (156) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| SDE-EM | 16.3 | (157) | `READY_EXPLICIT_INPUTS` |
| SDE-MIL | 16.3 | (159) | `READY_EXPLICIT_INPUTS` |
| SDE-MIL-EL14-PRINTED (alias FOC-EL14-MILSTEIN-PRINTED) | 16.3, 12.3 | (158), (113) | `REFERENCE_DATA_ONLY` |
| SDE-PC-SANDROOS14 | 16.3 | (160) | `REQUIRES_USER_DECISION` |
| SDE-BC-MIRROR | 16.4 | (161), (162) | `READY_EXPLICIT_INPUTS` |
| SDE-DT-VDB26 | 16.6 | (165) | `READY_EXPLICIT_INPUTS` |
| SDE-DT-SANDROOS14 | 16.6 | (166) | `REQUIRES_USER_DECISION` |
| SDE-DT-SE17 | 16.6 | — | `READY_EXPLICIT_INPUTS` |
| MC-ROT | 16.7 | (167), (168), (169) | `READY_EXPLICIT_INPUTS` |
| MC-ROT-HEMI | 16.7 | (170) | `REQUIRES_USER_DECISION` |
| MC-REDRAW | 16.7 | (171) | `REQUIRES_USER_DECISION` |
| MC-SPHERE-FLORINSKI13 | 16.7 | (168) | `REQUIRES_SOURCE_OR_CODE_AUDIT` |
| VAL-PUSH-MOM | 16.9 | (174) | `READY_EXPLICIT_INPUTS` |
| VAL-PUSH-KAPPA | 16.9 | (175) | `READY_EXPLICIT_INPUTS` |
| VAL-PUSH-FOC | 16.9 | (176) | `READY_EXPLICIT_INPUTS` |
| VAL-PUSH-EST | 16.9 | — | `READY_EXPLICIT_INPUTS` |
| VAL-PUSH-KVAR | 16.9, 12.3 | (115) | `READY_EXPLICIT_INPUTS` |
| VAL-PUSH-FV | 16.9 | (178) | `READY_EXPLICIT_INPUTS` |
| VAL-PUSH-LEG | 16.9 | (177) | `READY_EXPLICIT_INPUTS` |
| VAL-PUSH-MMS | 16.9, 12.3 | (1) | `READY_EXPLICIT_INPUTS` |
| SDE-AMPS-PRINTED | 16.10 | — | `REQUIRES_SOURCE_OR_CODE_AUDIT` |

**Implementation requirement.** The registry is the single list of implementable items. An ID absent from it is not implemented under this specification; a merged ID (field `aliases`) is implemented once.

## 5. Empirical and phenomenological D_μμ shapes

Most focused-transport codes do not compute D_μμ from a turbulence spectrum. They prescribe a shape in μ and fix its amplitude from a prescribed λ∥ (or λ_r = λ∥cos²ψ), or they write the amplitude with an explicit rigidity factor. This section gives each published shape as printed, its operator convention, the parameter values in use with their locations, and the exact relation between amplitude and λ∥ through Eq. (23). MFP v1.2 Section 6 treats the same shapes; the values in Section 5.6 agree with its fixtures F-PA-01 to F-PA-05. Shapes computed from a turbulence spectrum are in Section 6.

**Conventions.** Unless a source states otherwise, D_μμ is the coefficient of the standard operator, Eq. (21). The frame of μ is stated per source where printed. Throughout, I(q,H) = ∫₀¹(1−μ²)/(μ^{q−1}+H) dμ. PARALLEL rev. 1.4 (Section 7.3, Eq. 23) defines the same integral over [−1, 1], which is twice this one.

### 5.1 Isotropic scattering (PA-ISO)

The isotropic coefficient, Eq. (24), is printed as D^iso_μμ = D₀(1−μ²) by [P1] (Eq. 11; Full text ×2, preprint), for "very strong turbulence", and as D^iso_μμ(μ) = D_I(1−μ²) by [P14] (Eq. 12a; PDF text layer, read twice) with D_I = v²/(6κ∥) (Eq. 13). The Barcelona papers write it as a constant scattering frequency ν = ν₀ with λ∥ = v/ν₀ [P48] (Full text ×2), and fit it as one option among ε-form runs [P49] (Full text ×1). [P5] state that the q-form "is reduced to the isotropic scattering model" for q = 1 (foundation; Full text ×1, equations garbled).

From Eq. (23), D₀(1−μ²) gives λ∥ = v/(2D₀) and κ∥ = v²/(6D₀), consistent with the printed D_I = v²/(6κ∥) (derived; Section 5.6). The process D^iso_μμ = ε(ν/2)(1−μ²) of [P25] (Eq. 14) is a component of the ε-form (Section 5.3), not a separate model. Full-orbit codes that redraw the velocity direction at Poisson-distributed times have no D_μμ (MFP v1.2, Section 6.6). EPREM writes an isotropic coefficient with a different operator factor; it is PA-EPREM (Section 5.4).

### 5.2 The regularized power law (PA-QFORM)

$$
D_{\mu\mu}=D_0\,(1-\mu^2)\left(|\mu|^{q-1}+H\right),\qquad 1<q<2,\quad H\ge 0 .
\tag{26}
$$

**Sources and symbols.** The form is attributed to [P50] (foundation; not read, known only through citing papers). It is printed in this form by [P1] (Eq. 10, "D^BW_μμ"; Full text ×2), [P4] (Eq. 19; ×2), [P11] (D_μμ,0(r,φ)(1−μ²){|μ|^{q−1}+H}; ×2, preprint; the arXiv PDF numbers it Eq. 2) and [P51] (Eq. 8, amplitude D(z); ×1 in the research pass, ×3 in MFP v1.2). [P52] print it with D₀ = ν₀/2 (Eq. 4; ×2). Here q is the spectral index of the slab (inertial-range) turbulence, and H is a constant that gives finite scattering through μ = 0, described as "an arbitrary (in terms of its value) correction to describe the inclusion of dynamical effects" [P1], "chosen in an ad hoc fashion" [P11], and as standing for non-resonant and nonlinear effects [P52, P7]. The Qin/Zhang group call the constant h (h₀ in [P15]); [P5] introduce it "to help particles transport through μ = 0 in simulation". With H = 0 the shape is the QLT result D^QLT_μμ = D₀(1−μ²)|μ|^{q−1} [P1] (Eq. 12), which "has the known problem of a resonance gap at μ = 0". The frame of μ is the solar-wind frame in [P5], [P27], [P53], [P1] ("p and μ … measured in the SW frame") and [P25]. [P7] scatter in the solar-wind frame with static fluctuations, and [P23] use the "plasma reference frame". In [P15] μ is the cosine "to the outward magnetic field line". No frame is stated in [P24], [P54], [P11], [P8] or [P19].

**Normalization.** Inserting Eq. (26) in Eq. (23) gives (derived; checked in Section 5.6)

$$
D_0=\frac{3v}{4\lambda_\parallel}\,I(q,H),\qquad I(q,0)=\frac{2}{(2-q)(4-q)} .
\tag{27}
$$

For H = 0 this is the printed ν₀ = 6v/[2λ∥(4−q)(2−q)] of [P52] (Eq. 2) and of the thesis [P55] (Eq. 3.8). For H > 0 no source prints a closed form; the sources that normalize to λ∥ obtain D₀ numerically ([P24]: "The level of scattering, D0, is obtained from λ∥"). Closed forms are in Section 5.6. [P8] print the equivalent λ∥ = (3v/4)∫₀¹(1−μ²)²/D_μμ dμ (Eq. 4) for this even shape. D_μμ has units of s⁻¹, so D₀, κ₀ and ν₀ are rates when no rigidity factor is written.

**Amplitude variants.** Several codes write the amplitude with explicit rigidity, momentum or speed factors. With D_μμ = A(1−μ²)(|μ|^{q−1}+H), Eq. (23) gives λ∥ = 3vI(q,H)/(4A) for any A that does not depend on μ (derived). The printed amplitudes and the dependences they imply:

| ID | Source (verification) | Printed amplitude A | λ implied by Eq. (23) (derived) |
|---|---|---|---|
| PA-QFORM-HE11 | [P19] Eq. 6 (×1, garbled; same form printed cleanly in [P20] Eq. 7, ×2); [P56] Eq. 4 (×1); [P21] Eq. 7 (R_d in place of R) | D^r_μμ ≡ D_μμ/cos²ψ = D₀vR^{−1/3}(\|μ\|^{q−1}+h)(1−μ²); D₀ "a constant"; R "the particle rigidity" (units not stated) | λ_r = λ∥cos²ψ = 3I(q,h)R^{1/3}/(4D₀); v cancels; λ_r is independent of position at constant D₀ |
| PA-QFORM-WQ15 | [P27] Eq. 3 (×1); [P16] Eq. 6 (×2); AMPS Type IV | D₀vp^{q−2}{\|μ\|^{q−1}+h}(1−μ²); p, μ in the solar-wind frame [Wang & Qin]; p in the solar-wind frame [Kubo] | λ∥ = 3I(q,h)p^{2−q}/(4D₀) |
| PA-QFORM-Z23 | [P15] Eq. 35 (×2, arXiv PDF) | D₀(x)p^{q−2}(1−μ²)(\|μ\|^{q−1}+h₀), "rigidity p is in the unit of GV"; Eq. 37 prints λ∥ = λ∥0(x)p^{2−q} | λ∥ = [3vI(q,h₀)/(4D₀(x))]p^{2−q}, so λ∥0 contains v (flagged below) |
| PA-QFORM-W19 | [P25] Eq. 7 (×3) | (Cπ/2)(\|Q\|/m)^{2−q}B^{−q}v^{q−1}; C from P = Ck^{−q} in the QLT Eq. 6 | λ∥ = [3I(q,H)/(2πC)](mv/\|Q\|)^{2−q}B^q; C is fixed by λ^r∥ = b_r²λ∥ (Eq. 11) |
| PA-QFORM | [P8] Eq. 6 (×2); [P7] Eq. 4.2 (×2); [P57] Eq. 2 (×1) | κ₀, which "contains the µ-independent terms in Eq. 2" (the QLT D_μμ) [Laitinen]; κ₀(r,E) [Kartavykh]; D^r_μμ ≡ D_μμ/cos²ψ = κ₀(r,R){\|μ\|^{q−1}+H}(1−μ²) [Dröge] | Eq. (27) with D₀ → κ₀ (and λ → λ_r for the D^r form) |
| PA-QFORM-QW15 | [P53] Eqs. 2–3 (×2); [P23] Eqs. 3–4 (×1); [P18] (×2), [P58]; AMPS Type V | D₀v(R_Lk_min)^{s−2}(μ^{s−1}+h)(1−μ²), D₀ = (δB_slab/B₀)²[π(s−1)/(4s)]k_min, R_L = pc/(\|q\|B₀) [QinWang]; [P18]: (δb_slab/B₀)²[π(s−1)/(4s)](v/l_slab)(R_L/l_slab)^{s−2}(μ^{s−1}+h)(1−μ²) | Eq. (28) |

With k_min = 1/l_slab the [P18] amplitude equals that of [P53]. For this slab form,

$$
\lambda_\parallel=\frac{3s\,I(s,h)}{\pi(s-1)}\left(\frac{B_0}{\delta B_{\rm slab}}\right)^2\frac{(R_Lk_{\min})^{2-s}}{k_{\min}}
\tag{28}
$$

(derived). At h = 0 this is MFP v1.2 Eq. 15, which MFP v1.2 (Section 6.4) identifies with the R_Lk_min ≪ 1 limit of the slab result of [P59] (Section 6). For [P25] Eq. 7 at H = 0, the amplitude follows from their Eq. 6 with P = Ck^{−q} and Ω = |Q|B/m (derived; check in derivation record DR-05). The notes do not record whether m is the rest mass or the relativistic mass.

**Published values (as printed).**

| Symbol | Value | Source and location |
|---|---|---|
| q | 5/3 | [P11] (text at the D_μμ equation); [P51] (text at Eq. 8); [P1] (Sect. 3.3.3, event application; Fig. 5); [P4] (text at Eq. 19); [P54] (text at Eq. 3); [P19]; [P20]; [P56]; [P27]; [P16]; [P15] ("Kolmogorov spectrum slope q = 5/3"); [P5] |
| s | 5/3 | [P24] (text at Eq. 3); [P53]; [P23]; [P18] |
| q | 1.67 | [P7] (text at Eq. 4.2); [P57] (figures); [P52] (Sect. 3.2) |
| q | 1.66 / 1.41 | [P52] (Sect. 3.1, cases A/B, H = 0); 1.66 also [P3] |
| q | 3/2 | [P1] (Fig. 5, QLT shape, H = 0, "Kraichnan") |
| q | 1.3–1.9, average 1.63 | measured slopes quoted by [P52] (Sect. 2) from Kunow et al. (1991) |
| H | 0.05 | [P11]; [P51] ("following Dröge et al. 2010"); [P1] (Sect. 3.3.3; Fig. 5); [P4]; [P24]; [P54]; [P57] (with q = 1.67); [P43] (secondary via [P25] and [P51]) |
| H | 0.1 | [P8] (text at Eq. 6) |
| H | 0.5 | [P7] (text at Eq. 4.2: "consistent with pitch angle distributions which are strongly forward peaked") |
| H | 0.01–0.50 (range studied) | [P52] (Sect. 3.2, Fig. 4) |
| h | 0.2 | [P5] (one of two main models); [P19]; [P20]; [P56]; [P16] ("[Qin et al., 2006]"; [P41] abstract only) |
| h₀ | 0.2 | [P15] (text at Eq. 35: "not very sensitive to h0 unless h0 ≪ 0.05") |
| h | 0.01 | [P5] (second model); [P27]; [P23]; [P53] ("we set h=0.01, and kmin=32 AU-1"); [P60] (k_min = 33 AU⁻¹, partly garbled); [P18] (Tables 1–2); AMPS Types IV, V |

**Printed issues.** (i) [P53], [P23], [P18] and [P58] print μ^{s−1}, not |μ|^{s−1}; for non-integer s this is undefined for μ < 0 unless |μ| is meant. (ii) [P23] and AMPS Type V print R^{s−2} with R a Larmor radius, which is dimensionally consistent only if R means R_Lk_min (MFP v1.2 D-6). (iii) [P53] use k_min = 32 AU⁻¹; [P60] and AMPS Type V print 33 AU⁻¹. (iv) [P15] Eq. 35 has no factor v, while Eq. 37 prints λ∥ = λ∥0(x)p^{2−q}; by Eq. (23) λ∥0 ∝ v/D₀(x), and the paper does not say where the speed dependence sits. (v) [P25] never give q a numerical value. [P8] give q no number; their slab spectrum has slope 5/3 (Eq. 10). (vi) The braces of [P5] Eq. 6 are not legible, and no rigidity factor is legible there. (vii) [P52] use q = 1.66 (Sect. 3.1) and q = 1.67 (Sect. 3.2) for nearly the same case.

**Statements for implementers in the sources.** [P11] (finite differences) impose zero flux through the μ = ±1 faces "because pitch-angle diffusion and focusing vanish" there, and average D_μμ over adjacent cells at interfaces. [P25] write of the μ = 0 singularity of Eq. 7 that "This pole needs to be treated carefully in order to obtain correct pitch-angle distributions", and that the Itô μ-SDE (their Eq. 13; cf. Eq. (25)) "requires the use of a very small time step to obtain correct PADs for values of μ near zero". They therefore adopt the ε-form with a different scheme (Section 5.3). The derivative ∂D_μμ/∂μ of Eq. (26) diverges as |μ|^{q−2} at μ = 0 for every H ≥ 0 (MFP v1.2, Section 21.5).

### 5.3 The ε-form (PA-EPS)

$$
D_{\mu\mu}=\frac{\nu_0}{2}\left(\frac{|\mu|}{1+|\mu|}+\varepsilon\right)(1-\mu^2),\qquad
\lambda_\parallel=\frac{v}{\nu_0}\,\varphi(\varepsilon).
\tag{29}
$$

**Sources.** The form originates in [P61] (journal paper not read). The first author's thesis [P55] (Ch. 3; Full text ×1, text extraction) gives ν(μ) = ν₀|μ|/(1+|μ|) (Eq. 3.19) and ν(μ) = ν₀(|μ|/(1+|μ|) + ε) (Eq. 3.21, ε symbol dropped by the extraction) with D_μμ = (ν/2)(1−μ²) (Eq. 3.7). Printed further by [P49] (×1), [P48] (Eq. 1; ×2), [P52] (Eq. 5; ×2), [P62] (Eq. 2, written D₀(…); ×1), [P29] (Eq. 41, D₀(…); ×2, preprint numbering) and [P14] (Eq. 12b, D_A(1−μ²)[|μ|/(1+|μ|)+ε]; PDF text layer, read twice). The λ relation and its small-ε asymptote φ(ε) = (3/2)ln(1/ε) are [P48] Eqs. 3–4. The thesis closed form of φ (Eq. 3.22) is not legible.

**Frame and meaning.** "In the solar wind frame, the pitch-angle diffusion coefficient can be expressed as D_μμ = ν(1−μ²)/2" [P49]; the thesis performs scattering "in the local solar wind frame". With ε = 0, ν vanishes linearly at μ = 0, which "guarantees that particles are not able to cross μ = 0 by scattering"; ε > 0 restores finite scattering through μ = 0, and a large ε approaches isotropic scattering ([P55]; [P49]: "from quasi-isotropic (ε ≥ 1) to fully anisotropic (ε = 0)"). For ε = 0.10 and 0.01, "Dμμ reaches its maximum value around μ ∼ ±0.5" [P55].

**PARADISE amplitude (PA-EPS-W19).** [P25] (Eq. 8; ×3) print

$$
D_{\mu\mu}=\frac{\nu_0}{2}\left(\frac{|Q|}{m}\right)^{2-q}B^{-q}v^{q-1}\left(\frac{|\mu|}{1+|\mu|}+\varepsilon\right)(1-\mu^2),
\tag{30}
$$

with ν₀ fixed through λ^r∥ = (3vb_r²/8)∫₋₁¹(1−μ²)²/D_μμ dμ (Eq. 11) and λ^r∥ = 0.3 AU "for 4 MeV protons", constant in r. Equation (29) then gives λ∥ = φ(ε)(mv/|Q|)^{2−q}B^q/ν₀ (derived). They split the coefficient into D^iso_μμ = ε(ν/2)(1−μ²) (Eq. 14; ε read in the ar5iv rendering, lost in the PDF extraction) and D^aniso_μμ = (ν/2)(1−μ²)|μ|/(1+|μ|) (Eq. 15), with ν = ν₀(|Q|/m)^{2−q}B^{−q}v^{q−1}, and treat scattering as "two isotropic scattering processes". [P29] write the ε-form without the rigidity factor and set the background λ⁰ = (0.1 au)(R/R_ref)^{2−q₀} with R_ref = 43 MV and q₀ = 5/3.

**Hemispheric transformation (PA-EPS-ETA).** For ε = 0 the thesis maps each hemisphere onto an isotropic problem ([P55], Eq. 3.20, as legible):

$$
\eta=\begin{cases}2\left(\mu-\tfrac12\right),&\mu\ge 0\\[2pt] 2\left(\mu+\tfrac12\right),&\mu<0\end{cases},\qquad
D_{\eta\eta}=\left(\frac{d\eta}{d\mu}\right)^2D_{\mu\mu}=\frac12\,\nu_0\,(1-\eta^2).
\tag{31}
$$

The back-transformation is not spelled out in the extracted text. With ε = 0, D_μμ = (ν₀/2)|μ|(1−|μ|), and Eq. (31) holds exactly. Because η is linear in μ within each hemisphere, ∂_μ(D_μμ∂_μf) = ∂_η(D_ηη∂_ηf) there. The point μ = 0 maps to η = ±1, where D_ηη = 0 (derived; check in derivation record DR-05). Section 16 gives the inverse map and the use of Eq. (31) as a scattering kernel, including the ε term.

**Published values.**

| ε | Source and location | Stated origin |
|---|---|---|
| 0.01 | [P49]; [P48] (inside r_bs); [P63] (Eq. 4, printed variant, Section 5.5); [P55] (Figs. 3.3–3.5) | equivalent to q = 1.66, H = 0 [P52] (Sect. 3.1, case A) |
| 0.10 | [P49]; [P55] | equivalent to q = 1.41, H = 0 (case B) |
| 0.05 | [P55] (Sect. 3.7.2) | to match Lario et al. (2007) |
| 0.045 | [P52] (Sect. 3.2) | least-squares equivalent of q = 1.67, H = 0.05 |
| 0.045 | [P14] (text at Eq. 13; GCR protons) | "parametrises dynamical effects" |
| 0.048 | [P25] (Sect. 4); [P62]; [P29] (after Eq. 41) | least-squares fit of Eq. 8 to Eq. 7 "for H = 0.05" (Wijsen; q not printed); "bridging the resonance gap" (Afanasiev, citing Klimas & Sandri 1971) |

The ε–(q, H) correspondences come from "minimizing the sum of the squared differences" between the two shapes [P52]. Neither the μ domain nor the weighting is stated, and ε(q, H) is given only as curves (their Figs. 2 and 4). An evaluation in this review under one assumed criterion (the same λ∥ for both shapes, unweighted squared difference of D_μμ over 0 ≤ μ ≤ 1) gives ε = 0.0100 (q = 1.66, H = 0), 0.106 (q = 1.41, H = 0), 0.0484 (q = 1.67, H = 0.05) and 0.0490 (q = 5/3, H = 0.05). These are not presets; the printed values are those in the table.

**Printed issues.** (i) [P48] state that for ε = 0.01, "ν(μ = 1) = 3.5υ/λ". Since ν(1)λ∥/v = φ(ε)(½ + ε), this is 3.523 with the asymptotic φ and 3.605 with the exact φ (derived; Section 5.6). (ii) [P29] use the symbol q₀ both for the Kolmogorov index in λ⁰ and for the exponent of their foreshock fit (their Eq. 21). (iii) The ε-form of [P29] and [P62] carries no rigidity factor, whereas [P25] Eq. 8 does. When the amplitude is fixed from λ∥ at each energy the two printings coincide; they differ when one ν₀ is used for all energies (derived).

**Statements for implementers in the sources.** [P25] describe Eq. 8 as bounded on μ ∈ [−1, 1], in contrast with the pole of Eq. 7. [P14] chose it "because its derivative is well behaved around µ = 0". D_μμ(0) = εν₀/2 > 0, and ∂D_μμ/∂μ jumps from −ν₀/2 to +ν₀/2 at μ = 0 (derived). [P25] use the fixed-δt small-angle rotation of [P55], stated for δt ≪ ν⁻¹ (their Eq. 16), instead of Eq. (25) (Section 16).

### 5.4 Other published shapes

**Dröge V_A form (PA-DROGE-VA).** [P64] (Eq. 8; Full text ×2, preprint) print, as "derived from a dynamical quasi-linear theory (Dröge 2003)",

$$
D_{\mu\mu}(\mu)=\frac{3v}{2(4-q)(2-q)\lambda_\parallel}\left(\sqrt{\mu^2+(V_A/v)^2}\right)^{q-1}(1-\mu^2),\qquad 1<q<2 .
\tag{32}
$$

The original [P65] was read as abstract only. Here q "corresponds to the spectral index of the power spectrum density of magnetic field fluctuations", V_A = 100 km s⁻¹ ("≪ v"), the prior on q is uniform on [1.1, 1.9], and λ∥ = ξv/Ω_n (Eq. 9; Ω_n undefined, MFP v1.2 D-29). The pitch-angle term of their Eq. 4 is ∂_μD_μμ∂_μ[(1 − μvV^c_sw/c²)f]; it acts on (1 − μvV^c_sw/c²)f, not on f, and v is "measured in the local solar wind frame". At V_A = 0 the prefactor is Eq. (27) with H = 0. For V_A > 0 the λ∥ in the prefactor is not the value returned by Eq. (23) (Section 5.6). The paper prints no integral definition of λ∥.

**M-FLAMPA Kolmogorov shape (PA-KOLMO).** [P36] (Full text ×2, preprint) derive, for I₋ + I₊ = I_Ck^{−5/3} (Eq. 6.1),

$$
D_{\mu\mu}=\frac{v}{\lambda_{\mu\mu}}(1-\mu^2)|\mu|^{2/3},\qquad \lambda_{\mu\mu}=\frac{4}{\pi}\frac{B^2/\mu_0}{I_C}\,r_L^{1/3}
\tag{33}
$$

(Eq. 6.2), in the plasma (co-moving) frame; their scattering operator is the wave-frame sum Σ_σ∂_{μσ}[D^σ_μμ∂_{μσ}f] (Eq. 4.29). This is Eq. (26) with q = 5/3, H = 0 and D₀ = v/λ_μμ, so Eq. (27) gives λ_xx/λ_μμ = 27/14, in agreement with the printed λ_xx = (54/7π)[(B²/μ₀)/I_C]r_L^{1/3} (Eq. 6.6). The approximate coefficients 0.9 (Eq. 6.10) and 0.5 (Eq. 6.11) are inconsistent with the stated convention k₀⁻¹ = L_max/2π (MFP v1.2 D-7; 81/(7π(2π)^{2/3}) = 1.0817, not "≈ 0.92"). M-FLAMPA itself solves the Parker equation [P66]; the shape is used to derive D∥.

**AMPS types (PA-AMPS-I…V).** [P10] (Table 1; arXiv only, partly garbled extraction) list the types below. The extraction strips absolute-value bars in places, and the readings differ: MFP v1.2 (Section 6.7) transcribes μ^{2/3}, μ^{q−1} and μ^{s−1}, the research notes |μ|. The |μ| below is therefore not verified.

| Type | Attributed to | Legible printed form and parameters | Printed issues |
|---|---|---|---|
| I | [P36] | (v/λ_μμ)(1−μ²)\|μ\|^{2/3}; λ_μμ = 0.5[(B²/μ₀)/(w₋+w₊)](L²_max r_L0)^{1/3}(pc/1 GeV)^{1/3}; r_L0 the Larmor radius at 1 GeV/c, L_max = 2πk_L⁻¹ | coefficient 0.5 (MFP v1.2 D-7) |
| II | [P67] | D₀[1/(1+ε)](1 − μV_A/v)²(…) + D₀[ε/(1+ε)](1 + μV_A/v)²(…), D₀ = (π/8)A²Ω²l_b(1−μ²); ε the backward/forward wave-energy ratio, ε = 1; V_A ≈ 38 km/s; l_b = 0.03 r⁻¹ AU upstream; δB⊥² = 0.1B₀²(1 AU/r) | the factors (…) are not legible; not fully specified |
| III | [P68]; Jokipii (1966); Zhao & Li (2014) | (π/4)(1−μ²)Ω₀kP^slab(k)/B², P(k) = A_βλ_c(δB)²/[1+(kλ_c)^β], ∫P dk = (δB)²; β = 5/3, λ_c = 10⁹ m, (δB/B₀)² = 0.05 at 1 AU | spectrum-based (Section 6); normalization symbol A_β vs A_s |
| IV | "Wang and Qin, 2004" | D₀vp^{q−2}(\|μ\|^{q−1}+h)(1−μ²); q = 5/3, h = 0.01 | year as printed; the matching paper found is [P27] |
| V | [P60] | (δB_slab/B₀)²[π(s−1)/(4s)]k_min vR^{s−2}(\|μ\|^{s−1}+h)(1−μ²); s = 5/3, h = 0.01, "k_min = 1/l_slab = 33 AU⁻¹" | R^{s−2} with R a Larmor radius (MFP v1.2 D-6) |

The paper does not state which type produced its figures. Its μ-update (Eq. 4) has no ∂D_μμ/∂μ term and is printed with dW_μ = dt·Λ (Section 16).

**EPREM (PA-EPREM).** [P33] (Full text ×2, preprint) print the scattering term and coefficient

$$
\frac{\partial}{\partial\mu}\left(\frac{D_{\mu\mu}}{2}\frac{\partial f}{\partial\mu}\right),\qquad
D_{\mu\mu}=\frac{(1-\mu^2)\,v}{2\lambda_\parallel},
\tag{34}
$$

(Eqs. 1–2), with λ∥ = λ₀(pc/1 GeV)^{1/3}(R/1 AU)^{2/3} and λ₀ = 0.05 AU (Eq. 3). The same operator is printed by [P69] (Eq. 1; ×1) and by [P17] (Eq. A1; ×2 renderings), who call D_μμ "the pitch-angle diffusion tensor" and work in "a frame co-moving with the solar wind plasma". [P34] Eq. 2 is not legible (a factor (R₁/r) with exponent ±3/2 of unreadable sign; MFP v1.2 Section 14.3). The coefficient of Eq. (21) is (1−μ²)v/(4λ∥), for which Eq. (23) returns 2λ∥ (derived; MFP v1.2 D-2, U-8).

### 5.5 Printed variants that are not the standard form

**Lang/Lavasa (PA-QFORM-LANG-PRINTED).** [P24] (Eq. 3; two arXiv renderings; publisher version not checked) and [P54] (Eq. 3; publisher and arXiv HTML agree) print

$$
D_{\mu\mu}=D_0\left(|\mu|^{q-1}+H\right)(1-\mu)^2
\tag{35}
$$

(Lang write s for q). Both use q = s = 5/3 and H = 0.05, with D₀ "obtained from λ∥" through the standard Eq. (23) (their Eq. 2). Lang et al. attribute the form to [P50, P70, P71]; those papers were not read here, but every other printing of this shape found in this review has (1−μ²). Taken literally (derived; Section 5.6): (i) D_μμ(μ)/D_μμ(−μ) = [(1−μ)/(1+μ)]², which is 1/9 at μ = 0.5; (ii) D_μμ(1) = 0 but D_μμ(−1) = 4D₀(1+H), so the no-flux condition at μ = −1 does not follow from the coefficient; (iii) the normalization stays finite, D₀ = (3v/4λ∥)I⁺(q,H) with I⁺(q,H) = ∫₀¹(1+μ²)/(μ^{q−1}+H)dμ. Then I⁺(q,0)/I(q,0) = 3 − q, and at q = 5/3, H = 0.05 the printed form needs a D₀ 1.473479 times that of Eq. (26) for the same λ∥ (MFP v1.2 D-3, U-5).

**Pacheco (PA-EPS-PACHECO-PRINTED).** [P63] (Eq. 4; publisher HTML and arXiv PDF agree) print

$$
D_{\mu\mu}=\frac{\nu_0}{2}\left(\frac{|\mu|}{1-|\mu|}+\varepsilon\right)(1-\mu^2)
=\frac{\nu_0}{2}\left[|\mu|(1+|\mu|)+\varepsilon(1-\mu^2)\right]
\tag{36}
$$

(the second form is derived), with ε = 0.01, λ_r logarithmically spaced between 0.01 and 0.5 AU and constant along the field line. The paper cites [P52], which print (1+|μ|). Taken literally (derived): (i) the bracket diverges at |μ| = 1, but the product is finite, with D_μμ(±1) = ν₀; (ii) for ε < 1, D_μμ increases monotonically with |μ| and is largest at |μ| = 1, whereas Eq. (29) peaks near |μ| ≈ 0.5; (iii) no-flux conditions at μ = ±1 must be imposed explicitly; (iv) λ∥ = (v/ν₀)φ_P(ε), with φ_P(0.01) = 5.287422 against φ(0.01) = 7.069533, and φ − φ_P → 2 as ε → 0 (MFP v1.2 D-4).

### 5.6 Exact normalization to λ∥

Every relation below is derived from Eq. (23) and was checked against direct quadrature at 40 digits (mpmath; scripts and outputs in derivation record DR-05). The values agree with MFP v1.2 fixtures F-PA-01 to F-PA-05.

**q-form.** By the substitution t = μ^{q−1} (Euler integral), and by a Pfaff transformation,

$$
I(q,H)=\frac1H\left[{}_2F_1\!\left(1,\tfrac{1}{q-1};\tfrac{q}{q-1};-\tfrac1H\right)-\tfrac13\,{}_2F_1\!\left(1,\tfrac{3}{q-1};\tfrac{q+2}{q-1};-\tfrac1H\right)\right]
=\frac{1}{1+H}\left[{}_2F_1\!\left(1,1;\tfrac{q}{q-1};\tfrac{1}{1+H}\right)-\tfrac13\,{}_2F_1\!\left(1,1;\tfrac{q+2}{q-1};\tfrac{1}{1+H}\right)\right].
\tag{37}
$$

The first form needs analytic continuation for H < 1. The second has its argument in (0, 1] and converges for all H ≥ 0, and at H = 0 it reduces to I(q,0) = 2/[(2−q)(4−q)] (derived). For q = 5/3 (substitution μ = t³):

$$
I\!\left(\tfrac53,H\right)=3\left[\frac67+\frac H5-\frac{H^2}{3}+H^3-\frac{(H+H^4)\arctan\!\left(1/\sqrt H\right)}{\sqrt H}\right],
\tag{38}
$$

with I(5/3, H) = 18/7 − (3π/2)√H + O(H) for small H. In double precision, Eq. (38) had relative errors ≤ 2.3 × 10⁻¹⁶ at the tested points H = 10⁻⁸, 10⁻⁴, 0.05 and 0.5. Cancellation raises the error to 1.5 × 10⁻¹² at H = 10 and 9 × 10⁻⁸ at H = 100 (derived). The printed variant of Section 5.5 has I⁺(q,H) = Eq. (37) with the minus signs before ⅓ replaced by plus signs, and I⁺(5/3,H) = 3[8/7 − H/5 + H²/3 − H³ − (H − H⁴)arctan(1/√H)/√H].

**ε-form.** With a = 1 + ε, polynomial division of the integrand gives

$$
\varphi(\varepsilon)=\frac32\int_0^1\frac{(1-\mu)(1+\mu)^2}{a\mu+\varepsilon}\,d\mu
=\frac{3}{2a^4}\left[\frac16+2\varepsilon+\frac52\varepsilon^2+\frac23\varepsilon^3+(1+2\varepsilon)\ln\frac{1+2\varepsilon}{\varepsilon}\right],
\tag{39}
$$

equivalently (3/2)[2/(3a) − 1/(2a²) + ε/a³ + (1+2ε)a⁻⁴ln((1+2ε)/ε)]. The bracket over a⁴ is the factor printed by [P14] (Eq. 13), D_A = v²/[4κ∥(1+ε)⁴][…], since ν₀ = 2D_A and κ∥ = vλ∥/3. That equation was read from a single rendering, and this check shows the reading is consistent with Eq. (23). For small ε, φ = (3/2)ln(1/ε) + 1/4 + O(ε ln ε), so the asymptote of [P48] (Eq. 4) omits the constant 1/4. In double precision, Eq. (39) agreed with 40-digit values to within a few units in the last place (relative errors below 7 × 10⁻¹⁶) at nine tested points from ε = 10⁻¹² to 10³; the exact figure depends on the points and the implementation. For the Pacheco printing, with b = 1 − ε and c = ε/b (0 < ε < 1),

$$
\varphi_P(\varepsilon)=\frac32\int_0^1\frac{(1-\mu)^2(1+\mu)}{b\mu+\varepsilon}\,d\mu
=\frac32\left[\frac{c^2+c/2-7/6}{b}+\frac{1-2\varepsilon}{b^4}\ln\frac1\varepsilon\right].
\tag{40}
$$

**Dröge V_A form.** With a = V_A/v, the ratio of the λ∥ returned by Eq. (23) to the λ∥ in the prefactor of Eq. (32) is

$$
\frac{\lambda_\parallel^{\rm int}}{\lambda_\parallel}=\frac{(4-q)(2-q)}{2}\,a^{1-q}\left[{}_2F_1\!\left(\tfrac{q-1}{2},\tfrac12;\tfrac32;-\tfrac1{a^2}\right)-\tfrac13\,{}_2F_1\!\left(\tfrac{q-1}{2},\tfrac32;\tfrac52;-\tfrac1{a^2}\right)\right],
\tag{41}
$$

For a → 0 the ratio is 1 + [(4−q)(2−q)/2][√π Γ((q−2)/2)/(2Γ((q−1)/2))]a^{2−q} + o(a^{2−q}) (derived). The coefficient is negative for 1 < q < 2, so the deficit decays only as (V_A/v)^{2−q}.

**Values** (to 7 digits; 15-digit values in derivation record DR-05):

| Shape | Parameters | Quantity | Value |
|---|---|---|---|
| PA-ISO | — | D₀λ∥/v | 1/2 |
| PA-QFORM | q = 5/3; H = 0, 0.01, 0.05, 0.1, 0.2, 0.5 | D₀λ∥/v = (3/4)I | 27/14 = 1.928571; 1.601995; 1.269780; 1.067701; 0.8398832; 0.5374346 |
| PA-QFORM | q = 1.67, H = 0.05 / 0.5 / 0 | D₀λ∥/v | 1.275869; 0.5383346; 1.950839 |
| PA-QFORM | H = 0; q = 1.66 / 1.63 / 1.5 / 1.41 | D₀λ∥/v | 1.885370; 1.710571; 1.2; 0.9816112 |
| PA-KOLMO | — | λ_xx/λ_μμ | 27/14 |
| PA-EPS | ε = 0.01, 0.045, 0.048, 0.05, 0.10 | φ = ν₀λ∥/v | 7.069533; 4.699191; 4.597262; 4.532860; 3.456958 |
| PA-EPS | ε = 0.01 | ν(1)λ∥/v exact / asymptotic | 3.605462 / 3.522955 |
| PA-EPS | ε = 0.045 ([P14]) | D_Aκ∥/v² | 0.7831984 |
| PA-QFORM-LANG-PRINTED | q = 5/3, H = 0.05 | D₀λ∥/v; ratio to PA-QFORM | 1.870995; 1.473479 |
| PA-EPS-PACHECO-PRINTED | ε = 0.01 / 0.048 | φ_P; ratio φ_P/φ | 5.287422, 0.7479167 / 3.218437, 0.7000770 |
| PA-DROGE-VA | q = 5/3; V_A/v = 0.001, 0.01, 0.05, 0.1, 0.3 | λ^int/λ | 0.9128698; 0.8123213; 0.6798308; 0.5985116; 0.4364452 |
| PA-DROGE-VA | V_A = 100 km/s, q = 5/3; protons 1, 10, 100 MeV (illustrative) | λ^int/λ | 0.8315350; 0.8849406; 0.9198294 |
| PA-EPREM | — | λ from Eq. (23) / printed λ∥ | 2 |

### 5.7 Summary of shapes

| ID | Printed form | Operator | Normalization | Parameter values in use | Sources |
|---|---|---|---|---|---|
| PA-ISO | D₀(1−μ²); (ν/2)(1−μ²); D_I(1−μ²) | standard | λ∥ = v/(2D₀); D_I = v²/(6κ∥) | — | [P1, P14, P48] |
| PA-QFORM | D₀(1−μ²)(\|μ\|^{q−1}+H), Eq. (26) | standard | Eq. (27), (37), (38) | q = 5/3, 1.67, 1.66, 1.41, 3/2; H = 0.05, 0.1, 0.5; h = 0.2, 0.01; h₀ = 0.2 | Section 5.2 |
| PA-QFORM-HE11, -WQ15, -Z23, -W19, -QW15 | D₀vR^{−1/3} (D^r form), D₀vp^{q−2}, D₀(x)p^{q−2}, (Cπ/2)(\|Q\|/m)^{2−q}B^{−q}v^{q−1}, slab D₀v(R_Lk_min)^{s−2} | standard | λ∥ = 3vI/(4A); Eq. (28) | as PA-QFORM | [P19, P27, P16, P15, P25, P53, P18] |
| PA-EPS | (ν₀/2)(\|μ\|/(1+\|μ\|)+ε)(1−μ²), Eq. (29) | standard; solar-wind frame | λ∥ = (v/ν₀)φ(ε), Eq. (39) | ε = 0.01, 0.10, 0.05, 0.045, 0.048 | Section 5.3 |
| PA-EPS-W19 | Eq. (30) | standard | λ∥ = φ(mv/\|Q\|)^{2−q}B^q/ν₀ | ε = 0.048; λ^r∥ = 0.3 AU at 4 MeV | [P25] |
| PA-EPS-ETA | Eq. (31), ε = 0 | standard, mapped to η | D_ηη = ½ν₀(1−η²) per hemisphere | — | [P55] |
| PA-DROGE-VA | Eq. (32) | standard, acting on (1 − μvV/c²)f | prefactor λ∥ ≠ integral λ∥ for V_A > 0, Eq. (41) | V_A = 100 km/s; q ∈ [1.1, 1.9] | [P64, P65] |
| PA-KOLMO | (v/λ_μμ)(1−μ²)\|μ\|^{2/3}, Eq. (33) | wave-frame sum (Eq. 4.29) | λ_xx/λ_μμ = 27/14 | q = 5/3, H = 0 | [P36] |
| PA-AMPS-I…V | Table in Section 5.4 | standard (Eq. 1) | not printed | Types IV, V: h = 0.01 | [P10] |
| PA-EPREM | (1−μ²)v/(2λ∥), Eq. (34) | ∂_μ[(D_μμ/2)∂_μf] | standard λ = 2λ∥ | λ₀ = 0.05 AU at pc = 1 GeV, 1 AU | [P33, P69, P17] |
| PA-QFORM-LANG-PRINTED | D₀(\|μ\|^{q−1}+H)(1−μ)², Eq. (35) | standard | D₀ = (3v/4λ∥)I⁺ | q = 5/3, H = 0.05 | [P24, P54] |
| PA-EPS-PACHECO-PRINTED | (ν₀/2)(\|μ\|/(1−\|μ\|)+ε)(1−μ²), Eq. (36) | standard | λ∥ = (v/ν₀)φ_P, Eq. (40) | ε = 0.01 | [P63] |

**Implementation requirement.** A shape object stores the printed form, the source's symbol for the gap parameter (H, h, h₀ or ε), the operator convention, and the amplitude mode: from λ∥, from λ_r with ψ, from a printed rigidity or momentum scaling (with the units of R or p, which [P19] and [P27] do not state), or from turbulence. The sources state that D₀, κ₀, C or ν₀ is fixed by λ∥ or λ_r through Eq. (23); the library evaluates that relation with the closed forms of Section 5.6 or by quadrature and checks the two against each other for every configured parameter set.

**Implementation requirement.** H, h, h₀ and ε are not converted into one another. The published correspondences (Section 5.3) depend on an unstated fitting criterion, so a conversion would be a new model, not a transcription.

**Implementation requirement.** The printed variants PA-QFORM-LANG-PRINTED and PA-EPS-PACHECO-PRINTED exist only under their own IDs. Because their coefficients do not vanish at μ = −1 (and at μ = +1 for Pacheco), a solver running them imposes the no-flux condition there explicitly. PA-EPS-PACHECO-PRINTED is evaluated in its product form, Eq. (36), never through the divergent bracket.

**Implementation requirement.** PA-DROGE-VA reports whether its λ∥ is the prefactor parameter or the value returned by Eq. (23). PA-EPREM uses D_μμ/2 as the coefficient of Eq. (21) in both drift and noise. Printings with μ^{s−1} are not evaluated for μ < 0 unless the configuration states that |μ|^{s−1} is meant.

## 6. Quasilinear D_μμ for magnetostatic slab turbulence

This section gives the resonant quasilinear (QLT) pitch-angle coefficient for magnetostatic slab turbulence and its explicit forms for the spectra that SEP codes use. D_μμ is the coefficient of the operator Eq. (21), defined as in Eq. (22). Dynamical turbulence and the dissipation range are in Section 7, propagating waves in Section 8, nonlinear broadening in Section 10, and measured spectra in Section 13. The companion PARALLEL rev. 1.4 (Section 8) treats the same physics with a one-sided total spectrum P_s; MFP v1.2 (Section 6.8) lists the printed forms and Section 7 the closed-form λ∥. All formulas below hold for one particle species; Ω carries the charge sign where stated.

### 6.1 Resonance condition and the general slab coefficient

**Source.** Teufel & Schlickeiser (2002) [P72], foundation (pre-2006). Full text ×2 (A&A HTML ALT text and PDF text layer) for Eqs. 2, 7, 24, 27–28, 30, 33, 34, 37; Eq. 25 only partly legible (ALT text truncated). The same coefficient is printed by Strauss et al. (2017) [P12] (Eqs. 17–19, Full text ×2) and Shalchi et al. (2009) [P73] (Eqs. 6–8, Full text ×1; the same form for |μ| ≫ δB/B₀ is their Eq. 29, fully legible). Engelbrecht et al. (2022) [P74] attribute it to Jokipii (1966) [P75] (secondary; the 1966 paper was not read).

**Set-up as printed.** B₀ = B₀e_z; magnetostatic fluctuations δE = 0, δB ≠ 0 (TS2002 Eqs. 4–5); all wave vectors parallel or antiparallel to B₀. μ = p∥/p along B₀ (Eq. 8); no plasma frame is named. Ω = eB₀/(mcγ) in Gaussian units (Eq. 7), with e the signed charge; R_L = v/Ω = pc/(|q|B₀) (Eq. 39, as printed; the two forms agree when Ω is taken as |Ω|, which is how R_L enters k_res below). D_μμ = Re∫₀^∞ dξ ⟨μ̇(t)μ̇*(t+ξ)⟩ (Eq. 2), which is the normalization of Eq. (22). With the slab tensor of Eq. 33, TS2002 print

$$
D_{\mu\mu}=\frac{\pi\Omega^2(1-\mu^2)}{B_0^2}\int_{-\infty}^{+\infty}dk_\parallel\,g(k_\parallel)\Big[(1-\sigma)\,\mathcal R(n=-1)+(1+\sigma)\,\mathcal R(n=+1)\Big],
\qquad
\mathcal R(n)={\rm Re}\int_0^\infty d\xi\,e^{-i(k_\parallel v_\parallel+n\Omega)\xi}F(\xi,k_\parallel),
\tag{42}
$$

(Eqs. 26, 34). F is the dynamical correlation function (Section 7); magnetostatic turbulence has F = 1. Then 𝓡(n) = πδ(k∥vμ + nΩ), as printed by Strauss et al. (2017, Eq. 17) and Shalchi et al. (2009, Eq. 7). The resonance is at |k∥| = |Ω|/(v|μ|) = 1/(R_L|μ|) (Strauss et al. 2017, Eq. 19; Engelbrecht et al. 2022, text of Eq. 11).

**The coefficient.** For σ independent of k∥ (Section 6.2) the delta functions give

$$
D_{\mu\mu}(\mu)=\frac{2\pi^2\,\Omega^2\,(1-\mu^2)}{B_0^2\,v\,|\mu|}\;g\!\left(k_{\rm res}\right),
\qquad
k_{\rm res}=\frac{|\Omega|}{v|\mu|}=\frac{1}{R_L|\mu|},
\qquad
(\delta B)^2=8\pi\int_0^\infty g(k_\parallel)\,dk_\parallel .
\tag{43}
$$

Status: as printed by Strauss et al. (2017, Eq. 18) in the form 2π²v(1−μ²)g(k_res)/(|μ|B₀²r_L²), here rewritten with v/r_L = Ω (derived rewriting); derived from Eq. (42) with F = 1 (DR-06 check C1).

**Spectrum convention.** g(k∥) = g(|k∥|) is even. The printed normalization is TS2002 Eq. 37 and TS2003 Eq. 5 [P59] (Full text ×2): (δB)² = Σ_{m=1}^{3}⟨δB_m²⟩ = 8π∫₀^∞ g dk∥, i.e. the total variance of the two transverse components (P₃₃ = 0 for slab), equal to 8π times the one-sided integral of g. Strauss et al. (2017, Eq. 11) call g the "omni-directional one-sided slab spectrum" with the same 8π. The formula involves g only through g/B₀², whose dimension is length; it is therefore the same in Gaussian and SI units as long as δB² and B₀² are expressed in the same units (derived).

Other printings of the spectrum map onto g as follows (derived; each relation reproduces the variance δB²):

| Spectrum as supplied | Its normalization | Relation to g |
|---|---|---|
| One-sided total P_s(k), k > 0 (PARALLEL rev. 1.4, Eq. 25) | ∫₀^∞ P_s dk = δB² | P_s = 8πg |
| Two-sided total I(k), k real | ∫₋∞^∞ I dk = δB² | I = 4πg |
| Two-sided, one transverse component S_xx(k) | ∫₋∞^∞ S_xx dk = ⟨δB_x²⟩ = δB²/2 | S_xx = 2πg |
| One-sided, one transverse component P_xx(k) | ∫₀^∞ P_xx dk = ⟨δB_x²⟩ | P_xx = 4πg |

With P_s, Eq. (43) is PARALLEL rev. 1.4 Eq. 26, D_μμ = πΩ²(1−μ²)P_s(k_res)/(4B₀²v|μ|). The spectra of TS2002, TS2003, Strauss et al. (2017, Eq. 12), Shalchi, Lazarian & Schlickeiser (2007) [P76] (Eqs. 4–5), Ding et al. (2022) [P42] (Eqs. 9–10), Shalchi et al. (2009, Eqs. 34–35) and Tautz & Shalchi (2013) [P77] (Eq. 4) all satisfy 8π∫₀^∞g dk = δB² (derived; checked numerically, DR-06 C2), except where cut-offs are added (Section 6.5).

**Implementation requirement.** A supplied spectrum is converted to g (or P_s) by a named rule before Eq. (43) is evaluated (PARALLEL rev. 1.4, Section 8.7). The prefactor is never rescaled to match a variance.

### 6.2 Polarization and helicity terms as printed

TS2002 define δB_L = (δB_x + iδB_y)/√2 and δB_R = (δB_x − iδB_y)/√2 (Eq. 11; visible parts legible) and print μ̇ = (iΩ/(√2B₀))√(1−μ²)[δB_R e^{iΦ} − δB_L e^{−iΦ}] (Eq. 12). The general coefficient, Eq. 25, contains Σ_n ∫d³k [J²_{n+1}P_RR + J²_{n−1}P_LL − J_{n+1}J_{n−1}(P_RL e^{2iΨ} + P_LR e^{−2iΨ})]𝓡 with argument k⊥v⊥/Ω (only the prefactor Ω²(1−μ²)/(2B₀²) and Σ_n∫d³k are legible in two renderings; the bracket is a single reading). The slab tensor is P_lm = g(k∥)[δ(k⊥)/k⊥][δ_lm + iσε_lm3] for l, m ∈ {1, 2} and 0 for l, m = 3 (Eq. 32), so that P_RR = (1−σ)gδ(k⊥)/k⊥, P_LL = (1+σ)gδ(k⊥)/k⊥, P_RL = P_LR = 0 (Eq. 33). For k⊥ → 0 only n = −1 couples to P_RR and n = +1 to P_LL, which gives Eq. (42).

**σ is not defined in TS2002 and the word "helicity" does not appear in the text read; σ is absent from the parameter list (Eq. 64).** TS2002 Eq. 35 is printed without σ, and Eq. 44 states D_μμ(−μ) = D_μμ(μ).

Derived (DR-06 C1, C14):

1. If σ does not depend on k∥ and g is even, the σ terms cancel exactly when Eq. (42) is folded onto k∥ > 0, for any F that depends on |k∥|. This reproduces Eq. 35 without σ and, for F = 1, Eq. (43).
2. For a real field, P_lm(−k) = P_lm(k)* (reality condition; the tensor is also Hermitian, P_lm = P_ml*). With g even this requires σ(−k∥) = −σ(k∥) for the form of Eq. 32.
3. For σ(k∥) = σ_m sgn(k∥) and F = 1, the coefficient is Eq. (43) multiplied by [1 − σ_m sgn(Ωμ)]. It is then not even in μ, unlike Eq. 44.

Which σ TS2002 intend is not stated. Spectra resolved by polarization and propagation direction are treated in Section 8.

**Implementation requirement.** The slab models of this section are specified for σ = 0. A helicity option requires user decision U-17 (Section 19.2).

### 6.3 Power-law spectrum: the q-form

For g = g₀|k∥|^{−s} at the resonant wavenumber, Eq. (43) gives (derived)

$$
D_{\mu\mu}=D_0\,(1-\mu^2)\,|\mu|^{s-1},
\qquad
D_0=\frac{2\pi^2 g_0\,|\Omega|^{2-s}\,v^{s-1}}{B_0^2}.
\tag{44}
$$

It holds for all μ whose k_res lies in the power-law range. It is the q-form Eq. (26) with q = s and H = 0. With the TS2003 normalization in its printed approximate form, (δB)² ≈ 8πs g₀k_min^{1−s}/(s−1) (Eq. 5, valid for k_min ≪ k_d), g₀ = (s−1)(δB)²k_min^{s−1}/(8πs) and (derived)

$$
D_0=\frac{\pi(s-1)}{4s}\left(\frac{\delta B}{B_0}\right)^2 v\,k_{\min}\,R^{\,s-2},
\qquad R=R_Lk_{\min}.
\tag{45}
$$

With the TS2002 spectrum (g = 0 below k_min, Eq. 37) the factor 1/s is absent. With the dissipation range kept exactly (Strauss et al. 2017, Eq. 12), D₀ is multiplied by N = s/[s + ((s−p)/(p−1))(k_min/k_d)^{s−1}] (derived). The identity Eq. (45) is fixture F-QLT-01.

**Printed forms.** Qin & Wang (2015) [P53] print D_μμ = D₀v(R_Lk_min)^{s−2}(μ^{s−1} + h)(1−μ²) with D₀ = (δB_slab/B₀)²[π(s−1)/(4s)]k_min, "adopted from Teufel & Schlickeiser 2003" (Eqs. 2–3, Full text ×2). For h = 0 and μ^{s−1} read as |μ|^{s−1}, this is Eq. (45) exactly (derived). R_L is "the maximum particle Larmor radius" and s = 5/3, k_min = 32 AU⁻¹, h = 0.01. Wang et al. (2014) [P23] (Eqs. 3–4) and AMPS Type V (Tenishev et al. 2022 [P10], Table 1, attributed to "Qin et al. 2013"; table read from a partly garbled extraction) print R^{s−2} with R a Larmor radius. That is dimensionally consistent only if R means R_Lk_min (MFP v1.2, D-6). AMPS Type V prints "k_min = 1/l_slab = 33 AU⁻¹" and h = 0.01. The term h is not part of QLT (Section 5). For h = 0, λ∥ is MFP v1.2 Eq. 15.

### 6.4 The TS2003 spectrum: piecewise coefficient

TS2003 (Sect. 2; PDF text layer, exponents resolved in the second reading) use g = g₀k_min^{−s} for |k∥| ≤ k_min, g₀|k∥|^{−s} for k_min ≤ |k∥| ≤ k_d, and g₀k_d^{p−s}|k∥|^{−p} for |k∥| ≥ k_d, with 1 < s < 2, 2 < p, k_min ≪ k_d. Inserting it in Eq. (43) gives, with R = R_Lk_min and Q = R_Lk_d (TS2003 Eq. 7) (derived),

$$
D_{\mu\mu}=D_{\rm TS}\,(1-\mu^2)\times
\begin{cases}
R^{-s}\,|\mu|^{-1}, & |\mu|\ge 1/R,\\
|\mu|^{s-1}, & 1/Q\le|\mu|<1/R,\\
Q^{\,p-s}\,|\mu|^{p-1}, & |\mu|<1/Q,
\end{cases}
\qquad D_{\rm TS}=\frac{2\pi^2 g_0\,|\Omega|^{2-s}v^{s-1}}{B_0^2}.
\tag{46}
$$

With the approximate normalization of TS2003 Eq. 5, D_TS equals D₀ of Eq. (45); with the exact one, it is multiplied by N of Section 6.3. The pieces join continuously at |μ| = 1/R and 1/Q, but ∂D_μμ/∂μ jumps there (derived). The first branch exists only for R > 1. Fixture F-QLT-02 gives values.

Checks (DR-06 C4, C6b): (i) the TS2003 dimensionless DT and RS coefficients (Section 7.2) tend to Eq. (46) as a → ∞ (α → 0); this agrees with rows (ii), (iii) and the μ-dependent part of row (iv) of TS2003 Table 1. (ii) Without a dissipation range (Q → ∞), Eq. (23) applied to Eq. (46) reproduces MFP v1.2 Eq. 45 to 40 digits (fixture F-QLT-03). (iii) With p ≥ 2, λ∥ diverges (Section 6.6).

**TS2002 spectrum.** With g = 0 below k_min and R > 1, D_μμ = 0 for |μ| > 1/R, so λ∥ is infinite in magnetostatic QLT (derived); the finite TS2002 high-rigidity λ∥ (K ∝ a², their Tables 3–4) comes only from the broadening of Section 7.

### 6.5 Bend-over spectrum

The spectrum g = [C(ν)/(2π)] l δB²_slab (1 + k∥²l²)^{−ν}, with C(ν) = Γ(ν)/[2√π Γ(ν − ½)] and inertial index s = 2ν, is printed by Shalchi, Lazarian & Schlickeiser (2007, Eqs. 4–5; Full text ×2) and Ding et al. (2022, Eqs. 9–10; Full text ×2). Tautz & Shalchi (2013, Eq. 4) write it with exponent −s/2 and C(s) = Γ(s/2)/[2√πΓ((s−1)/2)], the same function with ν = s/2. With r* = R_L/l, Eq. (43) gives (derived)

$$
D_{\mu\mu}=\pi C(\nu)\left(\frac{\delta B_{\rm slab}}{B_0}\right)^2\frac{v}{l}\;r_*^{\,s-2}\,(1-\mu^2)\,|\mu|^{s-1}\left(1+r_*^2\mu^2\right)^{-s/2}.
\tag{47}
$$

This is PARALLEL rev. 1.4 Eq. 29 (agreement checked to 40 digits, DR-06 C5). For r*|μ| ≪ 1 it is Eq. (44) with D₀ = πC(ν)(δB_slab/B₀)²(v/l)r*^{s−2}; for r*|μ| ≫ 1 it tends to πC(ν)(δB_slab/B₀)²(v/l)r*^{−2}(1−μ²)/|μ|. The exact λ∥ is MFP v1.2 Eq. 44 (fixture F-QLT-04) and PARALLEL rev. 1.4 Eq. 30.

Shalchi et al. (2009, Eqs. 34–35) add an energy-range index q: g = [D(s,q)/(2π)]δB²l_slab|k∥l_slab|^q/[1 + (k∥l_slab)²]^{(s+q)/2}, D(s,q) = Γ((s+q)/2)/[2Γ((s−1)/2)Γ((q+1)/2)] ("correctly normalized for q > −1 and s > 1"); 8π∫₀^∞g dk = δB² (derived check). Shalchi, Lazarian & Schlickeiser (2007, Eq. 4) and Tautz & Shalchi (2013, Eq. 4) also add g = 0 below k_min and a dissipation range (k_d/k∥)^p. With these cut-offs 8π∫g dk differs from δB²_slab (derived), which Shalchi et al. (2007) call "the turbulence strength", not the variance.

### 6.6 The 90° problem in QLT

Statements as printed:

| Source | Statement | Treatment in the source |
|---|---|---|
| Afanasiev et al. 2015 [P78] | QLT "does not work close to μ = 0" | linear gap filler below μ_min = Ω/(vk_max) (Section 13.7) |
| Chen et al. 2024 [P47] | κ∥ "blows up" as μ_min → 0; resonance near μ = 0 needs high-frequency waves that contribute little | κ∥ integral cut at μ_min = 0.05 (Section 13.2) |
| Shalchi et al. 2009 | QLT valid for \|μ\| > δB/B₀; non-resonant regime below | second-order QLT (Section 10) |
| Wijsen et al. 2019 [P25] | H "added to describe the scattering through μ = 0"; "dD_μμ/dμ has a pole at μ = 0" | H or ε-form (Section 5) |
| Strauss et al. 2017 | damping broadens the resonance and gives finite D_μμ at μ ≈ 0, which observations require | dynamical QLT (Section 7.4) |
| Zhang et al. 2023 [P15] | h₀ describes "the enhancement of scattering through μ = 0 … through either non-resonant scattering or non-linear effects" | h₀ = 0.2 (Section 5) |
| TS2002, TS2003 | not discussed | dynamical broadening gives finite D_μμ(0) for α > 0 (Section 7.3) |

Derived consequences of Eq. (43) (DR-06 C4):

1. For a power law with 1 < s < 2, D_μμ(0) = 0 but Eq. (23) converges, since (1−μ²)/D_μμ ∝ |μ|^{1−s}.
2. With a dissipation range of index p ≥ 2 (TS2003, Strauss et al. 2017: p = 3 and 2.6), D_μμ ∝ |μ|^{p−1} for |μ| < 1/Q and λ∥ diverges; the λ∥ integral cut at μ_c grows as μ_c^{2−p} for p > 2 and as ln(1/μ_c) for p = 2. PARALLEL rev. 1.4 (Section 8.5) states the same.
3. A hard upper cut-off k_max gives D_μμ = 0 for |μ| < 1/(R_Lk_max), so λ∥ is infinite.
4. With the TS2002 spectrum and R > 1, D_μμ = 0 for |μ| > 1/R (Section 6.4).

**Implementation requirement.** A magnetostatic slab model reports an infinite λ∥ where the integral diverges; it does not insert a μ cut-off or a floor silently (PARALLEL rev. 1.4, Section 8.5). For the Itô drift of Eq. (25), ∂D_μμ/∂μ of Eq. (44) diverges as |μ|^{s−2} at μ = 0, and Eqs. (46) and (47) have kinks or steep gradients where k_res crosses a spectral break. Time-step and grid convergence must be demonstrated for these shapes (MFP v1.2, Section 21.5).

### 6.7 Map of printed prefactors onto Eq. (43)

Each printed resonant coefficient below equals Eq. (43) when its spectrum is related to g as stated in the last column. This is algebra by this review (derived; each row checked numerically to 40 digits, DR-06 C10); it does not correct any source. k_res = |Ω|/(v|μ|) throughout.

| Source (eq.), status | Printed coefficient | Printed spectrum normalization | Agrees with Eq. (43) if |
|---|---|---|---|
| TS2002 Eq. 34 with F = 1; TS2003 | Eq. (42) | (δB)² = 8π∫₀^∞g dk | identical (σ independent of k) |
| Strauss et al. 2017 Eq. 18, Full text ×2 | 2π²v(1−μ²)g(k_res)/(\|μ\|B₀²r_L²) | Eq. 11: δB²_slab = 8π∫₀^∞g dk | identical (v/r_L = Ω) |
| Shalchi et al. 2009 Eq. 8, Full text ×1 | 2π²v(1−μ²)g/(\|μ\|B₀²R_L²) | Eqs. 34–35, 8π convention (derived) | identical |
| Ding et al. 2022 Eq. 8, Full text ×2 | 2π²Ω²(1−μ²)g^slab(k∥)/(B²vμ), "vμ" as printed | Eqs. 9–10, 8π convention (derived) | vμ is read as v\|μ\| |
| Wijsen et al. 2019 Eq. 6, Full text ×2 | (π/2)[(1−μ²)/(v\|μ\|)](Ω/B)²P(k_res) | not stated; P = Ck^{−q} | P = 4πg, i.e. ∫₋∞^∞P dk = δB² |
| Afanasiev et al. 2015 Eq. 2, Full text ×2 | (π/2)Ω\|k_res\|I_w,res(1−μ²)/B² (denominator: Section 13.7) | Eq. 5: δB² = ∫₋∞^∞I_w dk, sign of k = polarization | I_w = 4πg; consistent with the printed Eq. 5 for even I_w |
| Laitinen et al. 2016 Eq. 2, Full text ×2 | πΩ²(1−μ²)S∥(−(r_Lμ)⁻¹)/(v\|μ\|B²) | δB² = 2∫S dk, limits not printed | S∥ = 2πg if the integral runs over all k; if it runs over k > 0, S∥ = 4πg and Eq. 2 is twice Eq. (43) |
| Engelbrecht et al. 2022 Eq. 11, Full text ×2 | πv(1−μ²)S_xx(k_res)/(B₀²R_L²\|μ\|) | "power spectrum of b_x", not normalized in the text | S_xx = 2πg, the two-sided spectrum of one component |
| Chen et al. 2024 Eq. 3, Full text ×2 | (π/4)Ω₀(1−μ²)f_res P(f_res)/B₀² | not stated; P from B_N only | f_res P(f_res) = k_res P_s(k_res) with P_s = 8πg (Section 13.2) |
| Zhang et al. 2023 Eq. 36, not legible | reading A: π²Ω²(1−μ²)W⊥/(B²v\|μ\|); reading B: (π/2)Ω²(1−μ²)W⊥/(B²v\|μ\|) | not stated | A: W⊥ = 2g; B: W⊥ = 4πg (Section 13.5) |
| AMPS Type III (Tenishev et al. 2022, Table 1), partly garbled | (π/4)(1−μ²)Ω₀kP^slab(k)/B² | ∫P dk = (δB)² | P = 8πg if the integral runs over k > 0 |
| PARALLEL rev. 1.4 Eq. 26 | πΩ²(1−μ²)P_s/(4B₀²v\|μ\|) | ∫₀^∞P_s dk = δB_s² | P_s = 8πg |

A second, independent reading reaches the same relations for the rows it covers.

**Implementation requirement.** Every printed prefactor in this table is implemented together with its spectrum relation; a provider spectrum whose normalization is not stated (Wijsen, Engelbrecht, Chen, Zhang) needs user decisions U-10, U-22, U-23 and U-24 (Section 19.2). Model IDs: QLT-SLAB (Eq. (43), any converted spectrum), QLT-SLAB-PL (Eqs. (44), (45); with PA-AMPS-V and PA-QFORM reused for the printed q-forms), QLT-SLAB-TS03 and QLT-SLAB-TS02 (Section 6.4), QLT-SLAB-BENDOVER (Eq. (47)), QLT-SLAB-HEL (Section 6.2; requires U-17).

## 7. Dynamical turbulence and the dissipation range

In dynamical turbulence the correlation function F of Eq. (42) decays in time. The resonance then has a finite width, and D_μμ(0) becomes finite. This section covers the slab models of Teufel & Schlickeiser (2002, 2003) with the damping (DT) and random-sweeping (RS) correlation functions, the damped form of Strauss et al. (2017), the nonlinear-damping model of Shalchi, Lazarian & Schlickeiser (2007), the NADT model of Shalchi et al. (2006) as far as it is legible, and the dissipation-range coefficient used for electrons. Models whose resonance contains the wave frequency with separate forward and backward spectra are in Section 8. PARALLEL rev. 1.4 (Section 9) gives the generic Lorentzian and Gaussian kernels; MFP v1.2 (Section 7) gives the λ∥ closed forms built on these coefficients.

### 7.1 Damping and random-sweeping correlation functions

**Source.** Teufel & Schlickeiser (2002) [P72] and (2003) [P59], foundation (pre-2006). Full text ×2 (A&A ALT text and PDF text layer) for TS2002 Eqs. 24, 27, 28, 30, 35, 36; visible parts of Eqs. 26, 29, 31 agree; TS2003 Eqs. 3–4 from the PDF text layer, consistent with the visible ALT parts. The correlation functions are attributed to Bieber et al. (1994) [P70] (secondary; not read).

**As printed.** P_αβ(k, ξ) = P⁰_αβ(k)F(ξ, k) (Eq. 24). Damping model F₁ = e^{−ξ/q_D} (Eq. 27); random-sweeping model F₂ = e^{−(ξ/q_D)²} (Eq. 30), "with the same q_D"; q_D = 1/(αv_A|k∥|), 0 ≤ α ≤ 1 (Eq. 28). The resonance functions of Eq. (42) are 𝓡₁ = q_D/[1 + q_D²(k∥v∥ + nΩ)²] (Eq. 29) and 𝓡₂ = (√π/2)q_D e^{−(k∥v∥+nΩ)²q_D²/4} (Eq. 31). The coefficients are

$$
D_{\mu\mu}^{\rm (DT)}=\frac{2\pi\Omega^2(1-\mu^2)}{B_0^2}\int_0^\infty dk_\parallel\,g(|k_\parallel|)\,q_D\left[\frac{1}{1+q_D^2(k_\parallel v_\parallel-\Omega)^2}+\frac{1}{1+q_D^2(k_\parallel v_\parallel+\Omega)^2}\right],
\tag{48}
$$

$$
D_{\mu\mu}^{\rm (RS)}=\frac{\pi^{3/2}\Omega^2(1-\mu^2)}{B_0^2}\int_0^\infty dk_\parallel\,g(|k_\parallel|)\,q_D\left[e^{-(k_\parallel v_\parallel+\Omega)^2q_D^2/4}+e^{-(k_\parallel v_\parallel-\Omega)^2q_D^2/4}\right]
\tag{49}
$$

(TS2002 Eqs. 35–36; TS2003 Eqs. 3–4), with v∥ = vμ and g normalized as in Eq. (43).

**Conventions as stated.** The wave frequency does not enter the resonance; v_A enters only through q_D, and no dispersion relation is stated. μ is measured along B₀. D_μμ(−μ) = D_μμ(μ) (TS2002 Eq. 44; TS2003 Eq. 12). The final expressions "do not depend of the charge sign" (TS2002; also stated in TS2003). TS2002 also use the symbol α for κ⊥/κ∥ in their Sect. 1 (two meanings, as printed).

**Derived (DR-06 C1, C6b).** Eqs. (48) and (49) follow from Eq. (42) with F₁ and F₂ when σ does not depend on k∥ (Section 6.2). They are PARALLEL rev. 1.4 Eq. 39 with the Lorentzian kernel of width Γ_dec = αv_A|k∥| (its Eq. 40) and the Gaussian kernel with Δ = 2αv_A|k∥| (its Eq. 42). For α → 0 both tend to Eq. (43).

### 7.2 Spectrum insertion and dimensionless forms

TS2002 use the spectrum with g = 0 below k_min, normalized by Eq. 37; TS2003 use the spectrum flat below k_min (Section 6.4), normalized by Eq. 5. Both keep the approximate normalizations (δB)² ≈ 8πg₀k_min^{1−s}/(s−1) and ≈ 8πsg₀k_min^{1−s}/(s−1). The variables are (TS2002 Eq. 39; TS2003 Eq. 7; Full text ×2) ε = v_A/v, a = 1/(αε) = v/(αv_A), b = a/2, R_L = v/Ω = pc/(B₀|q|), R = R_Lk_min, Q = R_Lk_d. TS2003 print

$$
D_{\mu\mu}^{\rm (DT)}=\frac{(s-1)\,v\,(1-\mu^2)\,a}{4s\,k_{\min}R_L^2}\left(\frac{\delta B}{B_0}\right)^2 I(\mu),
\qquad
D_{\mu\mu}^{\rm (RS)}=\frac{\sqrt\pi\,(s-1)\,v\,(1-\mu^2)\,b}{4s\,k_{\min}R_L^2}\left(\frac{\delta B}{B_0}\right)^2 I(\mu),
\tag{50}
$$

$$
I=A+\left(\frac RQ\right)^{s}[C-B]+D,\quad
A=\int_0^1\!dx\,x^{s-1}K\!\left(\frac xR\right),\
B=\int_0^1\!dx\,x^{s-1}K\!\left(\frac xQ\right),\
C=\int_0^1\!dx\,x^{p-1}K\!\left(\frac xQ\right),\
D=\int_1^\infty\!\frac{dx}{x}\,K\!\left(\frac xR\right),
\tag{51}
$$

with K(y) = 1/[1 + a²(μ−y)²] + 1/[1 + a²(μ+y)²] for DT and K(y) = e^{−b²(μ−y)²} + e^{−b²(μ+y)²} for RS (TS2003 Eqs. 13–15, 17–19). Status: prefactors and I, Full text ×2; the brackets are complete in the ALT text of TS2002 Eq. 47 and partly truncated in TS2003. TS2002 print the same with I = A + (R/Q)^s[C − B] (no D term) and prefactors without the factor s (Eqs. 41, 43, 45–47, 50–52; Full text ×2). Intermediate dimensional forms are TS2002 Eqs. 38, 40, 42 and TS2003 Eqs. 6, 8, 10.

**Printed inconsistency.** TS2003 Eq. 11, the dimensionless RS form, prints the denominator 4k_minR_L² without s (ALT text and PDF agree), while Eqs. 9, 10, 13, 17 and 27 carry s. Derived: Eq. 4 with the normalization of Eq. 5 gives the s (checked; without s the coefficient is larger by the factor s).

**Derived checks (DR-06 C6).** Eq. (50) equals the k-integrals Eqs. (48), (49) with the approximate TS2003 normalization to 40 digits at three parameter points (values: fixture F-DYN-01). With the exact normalization of Strauss et al. (2017, Eq. 12) both are multiplied by N of Section 6.3. In the variable x the kernel peaks at x = μR (A, D) and x = μQ (B, C) with a relative width 1/(a|μ|) (DT) or 1/(b|μ|) (RS).

### 7.3 Pitch-angle asymptotics

TS2003 print the following approximations of I(μ) (Tables 1–2; Full text ×2; inequality directions from the HTML tables), with f₁ = 2/(p−2) + 2/(2−s), f₂ = π/sin(πs/2) and γ ≈ 0.577 (Euler's constant by reading):

| Row | Conditions | I(μ), DT | Conditions | I(μ), RS |
|---|---|---|---|---|
| (i) | μR ≫ 1, μQ ≫ 1, aμ ≪ 1, a/R ≪ 1 | 2/s − 2 ln(a/R) | μR ≫ 1, μQ ≫ 1, bμ ≪ 1, b/R ≪ 1 | 2/s − γ − 2 ln(b/R) |
| (ii) | μR ≫ 1, μQ ≫ 1, aμ ≫ 1 | π/(aμ) | μR ≫ 1, μQ ≫ 1, bμ ≫ 1 | √π/(bμ) |
| (iii) | μR ≪ 1, μQ ≫ 1, aμ ≫ 1 | π(R^s/a)μ^{s−1} | μR ≪ 1, b/R ≪ 1, μQ ≫ 1 | 2/s − γ − 2 ln(b/R) |
| (iv) | μR ≪ 1, μQ ≪ 1, aμ ≫ 1 | f₁R^sQ^{2−s}/a² + π(R^sQ^{p−s}/a)μ^{p−1} | μR ≪ 1, b/R ≪ 1, μQ ≪ 1, b/Q ≪ 1 | 2/s − γ − 2 ln(b/R) |
| (v) | μR ≪ 1, μQ ≫ 1, aμ ≪ 1, a/R ≫ 1 | f₂R^s/a^s | μR ≪ 1, b/R ≫ 1, bμ ≪ 1, μQ ≫ 1 | Γ(s/2)R^s/b^s |
| (vi) | μR ≪ 1, μQ ≪ 1, aμ ≪ 1, a/R ≫ 1, a/Q ≫ 1 | f₁R^sQ^{2−s}/a² | μR ≪ 1, b/R ≫ 1, bμ ≪ 1, μQ ≪ 1, b/Q ≪ 1 | Γ(s/2)R^s/b^s |
| (vii) | as (vi) with a/Q ≪ 1 | f₂R^s/a^s | μR ≪ 1, b/R ≫ 1, bμ ≪ 1, μQ ≪ 1, b/Q ≫ 1 | Γ(p/2)Q^{p−s}R^s/b^p |
| (viii) | μR ≪ 1, μQ ≫ 1, aμ ≪ 1, a/R ≪ 1 | 2/s − 2 ln(a/R) | μR ≪ 1, b/R ≫ 1, bμ ≫ 1, μQ ≫ 1 | √π(R^s/b)μ^{s−1} |
| (ix) | μR ≪ 1, μQ ≪ 1, aμ ≪ 1, a/R ≪ 1, a/Q ≪ 1 | 2/s − 2 ln(a/R) | μR ≪ 1, b/R ≫ 1, bμ ≫ 1, μQ ≪ 1, b/Q ≫ 1 | √π(R^sQ^{p−s}/b)μ^{p−1} |

The D integral alone is approximated by TS2003 Eqs. 16 and 20 (Full text ×2): DT, D ≈ −2 ln(a/R) (μR ≪ 1, a/R ≪ 1, or μR ≫ 1, aμ ≪ 1), R²/a² (μR ≪ 1, a/R ≫ 1), π/(aμ) (μR ≫ 1, aμ ≫ 1); RS, D ≈ −γ − 2 ln(b/R), (R²/b²)e^{−b²/R²}, √π/(bμ) in the same cases.

TS2002 print the asymptotics of A, B, C (Eqs. 49, 54) and Tables 1–2 (Full text ×2); without the D term they differ where the flat range matters, e.g. (2/s)/(1 + a²μ²) for μR ≫ 1, μQ ≫ 1 and 2/s where TS2003 have 2/s − 2 ln(a/R). TS2002 Table 2 writes √πΓ(s)/[2^{s−1}Γ((s+1)/2)] where TS2003 write Γ(s/2); the two are equal by the duplication formula (derived, checked).

**Reading of the tables (derived; DR-06 C6b, C7).**

1. Row (iii) inserted in Eq. (50) gives Eq. (45), independent of a: the magnetostatic q-form holds where the resonance is narrow (aμ ≫ 1) and inside the inertial range. Rows (ii) and the μ-dependent part of (iv) likewise give the flat and dissipation branches of Eq. (46).
2. Rows (v)–(ix) of DT and (iii)–(vii) of RS do not depend on μ (RS row (i) is also μ-independent but applies only for μR ≫ 1), so D_μμ ∝ (1−μ²) near μ = 0 and D_μμ(0) is finite for α > 0 (fixture F-DYN-01 gives I(0) at one point).
3. Spot checks of the exact I(μ) against nine rows give ratios that tend to 1 as the inequalities are sharpened. Row (vi) of DT also needs R ≪ Q; its leading correction is of relative order (R/Q)^{2−s}.

The mean-free-path forms built on these coefficients are, as printed (Full text ×2), λ/λ₀ = [3/(s−1)]R²K/(k_min a) (TS2002 Eq. 58, DT), [3/(√π(s−1))]R²K/(b k_min) (Eq. 62, RS), [3s/(s−1)]R²K/(k_min a) (TS2003 Eq. 24) and [3s/(√π(s−1))]R²K/(b k_min) (Eq. 28), with λ₀ = (B₀/δB)² and K = ∫₀¹(1−μ²)/I dμ. They follow from Eq. (23) (derived). The tabulated K and the resulting λ∥ are in MFP v1.2, Section 7.

### 7.4 Damped QLT of Strauss et al. (2017)

**Source.** Strauss, Dresing & Engelbrecht (2017) [P12]. Eqs. 9–18 Full text ×2 (ar5iv and arXiv PDF); Eqs. 19–22 arXiv PDF only (Full text ×1). The typeset ApJ version and the LaTeX source could not be read.

**As printed.** Γ(k∥, t) = exp(iωt − γt) (Eq. 9); ω = jV_Ak∥ with j = +1 forward and −1 backward (Eq. 10). Spectrum: TS2003 shape with dissipation range, δB²_slab = 8π∫₀^∞g^slab dk∥ (Eq. 11) and the exact g₀ of Eq. 12. Then

D_μμ = [2πv²(1−μ²)/(B₀²r_L²)][1 − μω/(k∥v)]²∫₀^∞g^slab[𝓡_{n=+1} + 𝓡_{n=−1}]dk∥ (Eq. 13), with 𝓡_{n=±1} = γ/[γ + (vμk∥ − ω ± Ω)²] (Eq. 14); [1 − μω/(k∥v)]² = [1 − μV_A/v]² ≈ 1 (Eq. 15); 𝓡_{n=±1} ≈ γ/[γ + (vμk∥ ± Ω)²] (Eq. 16); 𝓡 → πδ(vμk∥ ± Ω) for γ → 0 (Eq. 17); Eqs. 18–19 are Eq. (43). The form used in the paper is

$$
D_{\mu\mu}^{DT}=\frac{2\pi v^2(1-\mu^2)}{B_0^2r_L^2}\int_0^\infty g^{\rm slab}(k_\parallel)\left[\frac{\gamma}{\gamma+(v\mu k_\parallel+\Omega)^2}+\frac{\gamma}{\gamma+(v\mu k_\parallel-\Omega)^2}\right]dk_\parallel ,
\qquad \gamma=\alpha\omega=\alpha V_Ak_\parallel
\tag{52}
$$

(Eqs. 20–21), α ∈ [0, 1], α = 1 unless stated. μ is "the cosine of the pitch angle" (frame not stated); Ω is the "particle gyro-frequency" (no expression printed).

**Printed inconsistency (as printed in Eqs. 14, 16 and 20).** The denominator is γ + (…)², not γ² + (…)². γ is a rate and (…)² a squared rate, so the printed kernel is dimensionally inconsistent. Derived (DR-06 C8): (i) ∫γ/(γ + x²)dx = π√γ, which tends to zero as γ → 0, so the printed kernel does not have the limit stated in Eq. 17; γ/(γ² + x²) integrates to π. (ii) With γ² in the denominator, Eq. (52) is identical to Eq. (48) with q_D = 1/γ (checked to 40 digits). (iii) The printed form changes its value when the unit of time is changed. For this reason the printed form is not an executable coefficient (Section 4.1); the record DYN-STRAUSS17-PRINTED exists so that the discrepancy can be cited and demonstrated, not so that it can be run.

**Parameters (Table 1, Full text ×1).** s = 5/3; p = 2.6 (Smith et al. 2006); k_min = 35 AU⁻¹ (to match Weygand et al. 2011); k_d = 2π(a + bΩ_i)/V_sw with a = 0.2 Hz, b = 1.76 (Leamon et al. 2000); δB²(1 AU) = 13.2 nT² (Bieber et al. 1994) ∝ r^{−2.4}; δB²_slab = 0.2δB²; α = 1. Text: B = 5 nT, V_sw = 400 km s⁻¹ and n = 5 cm⁻³ (∝ r^{−2}) at Earth; λ∥ ∼ 0.7 AU for 100 keV electrons at Earth. The symbol a is used both in the k_d fit and for the perpendicular-diffusion parameter.

### 7.5 Nonlinear damping of slab modes (Shalchi, Lazarian & Schlickeiser 2007)

**Source.** [P76], arXiv PDF and ar5iv. Eqs. 4–9, 12, 13 Full text ×2; Eqs. 17–29 Full text ×1 (PDF; Eq. 20 confirmed by a second, targeted reading).

**As printed.** Γ^slab = e^{iωt−γ_ct} (Eq. 6); γ_c = v_A√(k∥/L) (Eq. 7, after Lazarian & Beresnyak 2006; L the outer scale); ω = jv_Ak∥, "we simply assume that all waves propagate forward and we therefore set j = +1" (Eq. 9). The coefficient is πΩ²(1−μ²)B₀⁻²∫₋∞^∞dk∥ g^slab Σ_{n=±1}R_n (Eq. 17) with R_n = v_A√(k∥/L)/[v_A²k∥/L + (k∥v∥ + nΩ − v_Ak∥)²] (Eq. 19). Eq. 20 folds it onto k∥ > 0 with the prefactor 2πΩ²(1−μ²)/B₀² and is **printed with the sum over n = +1 only**; Eqs. 17, 21 and 24 use n = ±1. The dimensionless form is

$$
\tilde D_{\mu\mu}=\frac{C(\nu)(1-\mu^2)}{R}\frac{\delta B_{\rm slab}^2}{B_0^2}\int_0^\infty dx\,h(x)\sum_{n=\pm1}\frac{\epsilon R\sqrt{x\eta}}{(\epsilon R\sqrt{x\eta})^2+[xR(\mu-\epsilon)+n]^2},
\qquad \tilde D_{\mu\mu}=\frac{D_{\mu\mu}l_{\rm slab}}{v}
\tag{53}
$$

(Eqs. 24–26), with x = l_slab k∥, R = R_L/l_slab = v/(Ωl_slab), ε = v_A/v, η = l_slab/L. The spectrum is Eq. 4 with Eqs. 22–23: h = 0 for x < x_min, (1+x²)^{−ν} for x_min ≤ x < ξ, (1+ξ²)^{−ν}(ξ/x)^p for x ≥ ξ, ξ = l_slab k_d; C(ν) as in Section 6.5. Eqs. 27–29 express ε through R and the rest rigidity R₀ = (0.511 MV or 938 MV)/(l_slab B₀) and state the validity range.

**Derived (DR-06 C9).** With ω odd in k∥ and γ_c a function of |k∥|, Eq. 17 equals Eq. 20 with n = ±1. Eqs. 21 and 24 are exactly consistent with that form (checked to 40 digits). One reading flagged the numerator of Eq. 21 (v_A√(k∥/L), not divided by Ω) against Eq. 24; the factor Ω is absorbed in the prefactor of Eq. 21, so the two agree. Eq. 20 with n = +1 alone gives a different value (1.2% of the n = ±1 value at the test point). The frame of μ is not stated; since the resonance contains the wave frequency, μ is presumably measured in the plasma frame (an inference of this review, not a source statement). D_μμ is not even in μ because of (μ − ε).

**Parameters (Table 2).** 2ν = 5/3, p = 3, v_A = 33.5 km s⁻¹, l_slab = 0.030 AU, k_d = 3 × 10⁶ AU⁻¹ ("k_slab" in the table), B₀ = 4.12 nT, δB_slab/B₀ = 1; Fig. 3 uses x_min = η = 0.1 and 0.01.

### 7.6 NADT model of Shalchi et al. (2006), as legible

**Source.** [P79]; one rendering, an OCR'd scan in which the symbol font is lost. Only structure is legible; no prefactor can be transcribed.

Legible items: slab + 2D turbulence (Eq. 2); Γ_slab = e^{−γ_slab t}e^{iωt} (Eq. 7; no minus glyph is visible in the second exponent, while Table 1 lists the plasma-wave model as e^{−iv_Ak∥t}); ω = jv_Ak∥, j = +1 (Eq. 9); γ_slab = εZ/l_2D, ε "of order 1" (Eq. 8), Z = √2 v_A δB_2D/B₀ (Eq. 12), as extracted. The slab resonance function (Eq. 22) has the structure rate/[rate² + (k∥v∥ + nΩ − v_Ak∥)²], with the rate squared, and the slab coefficient (Eq. 23) is ∝ (1−μ²)/B₀² times ∫₀^∞dk∥ g^slab Σ_{n=±1}R_n. The 2D part is non-resonant: R_n^{2D} = γ_2D/[γ_2D² + (nΩ)²], weighted by n²J_n²(W)/W² with W = k⊥R_L√(1−μ²) (Eqs. 31–34), and Σ_n[n²/(n² + y²)]J_n²(x) ≈ (x²/2)/(x² + 2y² + 2) (Eq. 35). The prefactors of Eqs. 20, 23, 27, 28 and 34 contain unidentifiable glyphs. Table 2 (OCR): "2ν = 5/3", p = 3, v_A = 33.5 km s⁻¹, l_slab = 0.030 AU, l_2D = 0.1 l_slab, B₀ = 4.12 nT, δB/B₀ = 1, δB²_slab = 0.2δB².

### 7.7 The dissipation range and electrons

**TS2003.** Beyond k_d the spectrum falls as k∥^{−p} (p = 3 in TS2003 Eq. 50). In magnetostatic QLT this gives D_μμ ∝ |μ|^{p−1} for |μ| < 1/Q (Eq. (46)), and λ∥ diverges for p ≥ 2 (Section 6.6). With dynamical broadening, rows (iv) and (vi) of Table 1 (DT) and rows (vii) and (ix) of Table 2 (RS) contain Q and p. For μQ ≪ 1, aμ ≪ 1 and a/Q ≫ 1 (row vi, DT), Eq. (50) gives D_μμ = [(s−1)v/(4sk_minR_L²)](δB/B₀)²f₁R^sQ^{2−s}(1−μ²)/a, which is proportional to 1/a = αv_A/v (derived). The scattering of particles that resonate only in the dissipation range is then set by the broadening parameter. The electron λ∥ built on these rows (TS2003 Tables 3–4; Engelbrecht et al. 2022 [P74], Eqs. 38–39, extraction ambiguous) is in MFP v1.2, Section 7.2 (Eqs. 22–23; the DT printing differences are MFP v1.2 D-1).

**Strauss et al. (2017).** p = 2.6 and the onset k_d = 2π(a + bΩ_i)/V_sw (Section 7.4), used in Eq. (52).

**Strauss & le Roux (2019)** [P13] (dissipation-range part only; wave terms are in Section 8). Full text ×1 for Eqs. 15–17. The onset is set by thermal damping: k_d^{LH,RH} ≈ |Ω|_{p,e}/(V_A + 3v_th^{p,e}) (Eq. 15), left-hand waves being damped by thermal protons and right-hand waves by thermal electrons; ω^{LH,RH} ∈ |Ω|_{p,e} ± 3k∥v_th (Eq. 16); T_{p,e} = T⁰(r/r₀)^{−χ} (Eq. 17). Parameters: T_p⁰ **printed as "30^5 K"**, T_e⁰ = 10⁵ K, χ_p = 4/3, χ_e = 1/3, α = 1; the slab spectrum shape is "exactly the same as used by Strauss et al. (2017)" except for k_d.

**Implementation requirements.**

1. DYN-TS-DT and DYN-TS-RS evaluate Eqs. (48), (49) or the equivalent Eq. (50). They declare the spectrum (TS2002 or TS2003 shape) and the normalization (the approximate one of the source, which reproduces its tables, or the exact Eq. 12 of Strauss et al. 2017) as named options. The prefactor of TS2003 Eq. 11 is not used.
2. α > 0 and v_A > 0 are required for broadening. α → 0 is the magnetostatic model of Section 6.4, with λ∥ possibly infinite.
3. The narrow peaks at x = μR and x = μQ (relative width 1/(a|μ|)) are integrated with breakpoints or adaptively (PARALLEL rev. 1.4, Section 9.5). For a < ∞ the coefficient is smooth in μ, so the drift ∂D_μμ/∂μ of Eq. (25) is finite (derived); it is computed from the same integral, not by differencing a coarse table.
4. Eq. (52) is not a runtime coefficient: a kernel that is dimensionally inconsistent as printed has no value independent of the unit of time, and a configuration choice cannot supply the missing consistency (admissibility rule, Section 4.1). DYN-STRAUSS17-PRINTED is kept as a record (`REFERENCE_DATA_ONLY`, admissibility "inadmissible as printed") and may be evaluated only in a diagnostic that demonstrates its unit dependence (U-25, resolved). The γ²-form is DYN-TS-DT.
5. DYN-SLS08 implements Eq. (53) with n = ±1 and records that μ is not even-symmetric. DYN-NADT-SH06 is not fully specified. DYN-KD-SLR19 supplies k_d from Eq. 15 of Strauss & le Roux (2019).

## 8. Wave-based D_μμ: propagation direction, helicity, wave frames and dispersive waves

Section 6 treats magnetostatic slab QLT, in which the plasma, wave and fixed frames coincide. Here the waves propagate. The resonance then contains the wave frequency, forward and backward waves and the two polarizations carry separate spectra, and μ and v belong to a definite frame. The balanced limit and the canonical spectrum are those of PARALLEL rev. 1.4, Section 8.1 (Eqs. 25–26); its Sections 9.6 and 12.4 state what a directional wave provider must supply. Wave growth and transport are in Section 14.

The *plasma frame* is the rest frame of the bulk flow. A *wave frame* moves along **b** with the phase speed of the waves of one direction: +V_A ("forward", "outward", j = +1, +) or −V_A ("backward", "inward", j = −1, −) for Alfvén waves. Primed or subscripted quantities (p′, μ′, v_σ, μ_σ, μ̃) are wave-frame quantities.

### 8.1 Forward and backward Alfvén waves, the resonance condition and helicity weights

**Source.** Strauss & le Roux (2019) [P13]. Eqs. 5–13 Full text ×2 (arXiv PDF and ar5iv). Eqs. 24–29 Full text ×2 (arXiv PDF only; two fetches agree). Among the sources read, it prints the four-component bookkeeping most explicitly.

The slab spectrum is split into four components g^slab_{n,j}(k∥): j = +1 forward and j = −1 backward propagating; n = −1 right-handed and n = +1 left-handed polarization (conventions quoted in ar5iv). The normalization is one-sided in k∥ with prefactor 2π:

$$
\delta B_{slab}^{2}=2\pi\int_{0}^{\infty}\sum_{j=\pm1}\sum_{n=\pm1}g^{slab}_{n,j}(k_{\parallel})\,dk_{\parallel}
\tag{54}
$$

(their Eq. 9). The cross helicity and the magnetic helicity are defined as ratios (Eqs. 10–11, Full text ×2):

$$
\sigma_{c}(k_{\parallel})=\frac{\sum_{n}g_{n,j=+1}-\sum_{n}g_{n,j=-1}}{\sum_{n}g_{n,j=+1}+\sum_{n}g_{n,j=-1}},
\qquad
\sigma_{m}(k_{\parallel})=\frac{\sum_{j}g_{n=+1,j}-\sum_{j}g_{n=-1,j}}{\sum_{j}g_{n=+1,j}+\sum_{j}g_{n=-1,j}} .
\tag{55}
$$

The dynamical QLT coefficient (Eqs. 24, 25 and 7) is

$$
D_{\mu\mu}=\frac{\pi v^{2}(1-\mu^{2})}{B_{0}^{2}r_{L}^{2}}\int_{0}^{\infty}dk_{\parallel}\sum_{n=\pm1}\sum_{j=\pm1}\mathcal{R}_{n,j}\,g^{slab}_{n,j},
\qquad
\mathcal{R}_{n,j}=\frac{\tau^{-1}}{(\tau^{-1})^{2}+(v\mu jk_{\parallel}-\omega+n\Omega)^{2}},
\qquad
\tau^{-1}=\alpha|\omega|=\alpha|V_Ak_{\parallel}| .
\tag{56}
$$

Here ω = V_A k∥ is the wave frequency in the bulk plasma frame (their App. B), α = 1 unless stated, V_A ≪ v is assumed, and r_L is called the "maximal Larmor radius". Ω is not defined explicitly in the main text; App. B defines Ω_p and |Ω_e| = ΛΩ_p. The frame of μ: one reading of the paper records that μ is taken in the flow (plasma) frame (their Sect. 5); another reading found no explicit statement. In the negligible-damping limit, R_{n,j} ≈ πδ(vμjk∥ − ω + nΩ) (Eq. 27) and

$$
D_{\mu\mu}^{ND}=\frac{\pi^{2}v(1-\mu^{2})}{B_{0}^{2}r_{L}^{2}|\mu|}\sum_{n}\sum_{j}g^{n,j}_{slab}(k_{\parallel}^{res}),
\qquad
k_{\parallel}^{res}=\frac{n\Omega}{V_A-v\mu j}
\tag{57}
$$

(Eqs. 28–29; the indices are printed as g^{n,j}_{slab}). Unlike Strauss et al. (2017, Eq. 14) [P12], which prints γ/(γ + x²) (Section 7.4), Eq. 25 has the squared rate in the denominator. σ_c(r) = exp(−r/r₀) with r₀ = 1 AU (Eq. 18). The spectrum shape and the thermal dissipation-range onsets are in Section 7.7.

**Derived results (checked numerically; derivation record DR-08, checks S8-1 to S8-4).**

1. *Resonant components.* For Ω > 0, k^res > 0 selects (n, j) = (−1, +1), (+1, −1) for μ > V_A/v; (+1, +1), (+1, −1) for |μ| < V_A/v; and (+1, +1), (−1, −1) for μ < −V_A/v. The selected polarizations thus depend on the sign convention for Ω, which the paper leaves open for electrons.
2. *Jacobian.* The δ-limit of Eq. 25 gives πg(k^res)/|vμj − V_A|; Eq. 28 has v|μ| instead. For V_A/v = 0.05 and g_{n,j} ∝ k^{−5/3} with equal component amplitudes, they differ by 1.04 % at μ = 0.4 and by a factor 2.78 at μ = 0.02 (derived). Eq. 24 with α = 10⁻⁶ converges to the exact-Jacobian form (fixture F-WAVE-04).
3. *Balanced limit.* For V_A = 0 and g_{n,j} = P_s/(8π) (P_s of PARALLEL rev. 1.4, Eq. 25), Eq. 28 equals PARALLEL rev. 1.4, Eq. 26, and Eq. 24 equals Strauss et al. (2017), Eq. 13 term by term for the same Lorentzian.
4. *Split.* Eq. (55) gives two ratios for four components. g_{n,j} = (G/4)(1 + jσ_c)(1 + nσ_m) satisfies both, but the split used by Strauss & le Roux is not recorded in the readings (Decision U-26).

**Printed split with helicity weights.** Spanier, Schreiner & Schlickeiser (2022) [P80] (Eqs. 25–28, Full text ×1, arXiv PDF) give an explicit factorized form. It applies to a spectrum that is two-sided in k∥, where positive k∥ "describes forward moving waves with positive phase speed" (ar5iv):

$$
g^{b}_{L,R}=\frac{1-H_{L,R}}{2}g_{L,R}(k)\Theta(-k),\quad
g^{f}_{L,R}=\frac{1+H_{L,R}}{2}g_{L,R}(k)\Theta(k),\quad
g_{L}=\frac{1+\sigma(k)}{2}g_{tot},\quad
g_{R}=\frac{1-\sigma(k)}{2}g_{tot},
\tag{58}
$$

with (δB)² = 2π∫_{−∞}^{∞}dk∥(g_L + g_R) (Eq. 23, Full text ×2). H is defined per polarization (H_L ≠ H_R allowed). By the printed definitions, positive cross helicity means forward dominance and positive magnetic helicity means left-hand dominance in both papers. Schlickeiser (1989) [P81], the usual origin of the (1 ∓ H)(1 ∓ σ) weights, was not read.

**Other printed bookkeepings.**
- Ng & Reames (2008) [P38] use four modes σ ∈ {R+, L+, R−, L−}, with ambient I_R− = I_L− = 0.05 I_R+. Ng, Reames & Tylka (2012) [P82] use 0.1 I_R+.
- Afanasiev et al. (2015) [P78] use one direction (outward), with the sign of k giving the circular polarization.
- Sandroos et al. (2014) [P83] use positive λ for L waves and negative λ for R waves, with (1/μ₀)∫I^±_{L,R}(λ)dλ = C^±_{L,R}B₀²/(2μ₀) (Eq. 11).
- Chandran (2000) [P84] and Yan & Lazarian (2002) [P85] use C_ij = σM_ij (Section 9).

The Dröge (2003) form [P65] (MFP v1.2, Section 6.5, `PA-DROGE-VA`; Secondary) is the only empirical closed shape containing V_A/v found here. Vainio (2000) [P86] (Abstract only) reports an approximately constant proton mean free path for |H_c| ∼ 1 of left-hand waves; its formulas were not obtained.

### 8.2 Wave-frame pitch angle and the transformation to the plasma frame

**Printed frame statements.**
- Afanasiev et al. (2015, Eq. 1; Full text ×2): "υ and μ are the particle speed and pitch-angle cosine as measured in the wave rest frame".
- Afanasiev et al. (2023, Eq. 1) [P37] and (2025, Eq. 12) [P29]: the rest frame of Alfvén waves moving at V_A relative to the plasma. The shock-frame advection speed is vμ + (1 − M_A)V_A.
- Vainio (2003) [P87]: p and μ are "measured in the wave frame".
- Ng & Reames (2008): P, v and μ are in the "local outward wave frame" (W₊ = V_sw + V_A relative to the Sun). Inward modes have V_σ = −2V_A relative to it.
- Borovikov et al. (2019) [P36] (Full text ×1): "v_σ = v − V_σ**b**", f(p∥) → f(p∥ + m_iV_σ), and the operator Σ_σ∂/∂μ_σ(D^σ_μμ∂f/∂μ_σ) (Eq. 4.29).

**Energy change per scattering.** Vainio & Kocharov (2001) [P88] (Full text ×2 for the unnumbered relations) scatter isotropically in each wave frame, with D^±_μμ = ½(1 − μ²)ν_±. They state the plasma-frame energy loss caused by a wave-frame scattering as

$$
-\Delta E=\mp V_A\,p_\pm\,\Delta\mu_\pm ,
\tag{59}
$$

where p_± and μ_± are measured in the forward (+) or backward (−) wave frame. Vainio (2003) and Battarbee et al. (2011, Eq. 7) [P89] print ΔE_w = −V_A pΔμ for the outward waves. Sandroos et al. (2014) print ΔU = mV_AṼΔμ̃. The wave-energy side is in Section 14.3.

**Transformations (Derived).** For a wave frame moving with velocity V_w**b** relative to the plasma (V_w = ±V_A for Alfvén waves), the Lorentz transformation of a particle with wave-frame momentum p′, energy E′ = (p′²c² + m²c⁴)^{1/2} and cosine μ′ is

$$
p_{\parallel}=\gamma_w\left(p'\mu'+\frac{V_wE'}{c^2}\right),\quad
p_{\perp}=p'\sqrt{1-\mu'^2},\quad
E=\gamma_w\left(E'+V_wp'\mu'\right),\quad
\mu=\frac{p_\parallel}{(p_\parallel^2+p_\perp^2)^{1/2}},
\tag{60}
$$

with γ_w = (1 − V_w²/c²)^{−1/2}. The inverse uses −V_w. The Galilean limit (γ_w → 1, velocities) is

$$
\mu'=\frac{\mu v-V_w}{\sqrt{v^2-2vV_w\mu+V_w^2}},\qquad v'=\sqrt{v^2-2vV_w\mu+V_w^2},
\tag{61}
$$

which is identical to Schreiner & Spanier (2015, Eq. 2.9) [P90] with V_w = v_ph (fixture F-WAVE-01). Elastic scattering in the wave frame keeps p′ and E′ fixed. Because E is linear in μ′, Eq. (60) gives, for any finite step Δμ′,

$$
\Delta E=\gamma_w V_w\,p'\,\Delta\mu' .
\tag{62}
$$

This is Eq. (59) with γ_w = 1. The Galilean energy ½m|**v**′ + V_w**b**|² gives the printed form exactly. The relativistic relative difference is γ_w − 1 (5.0004 × 10⁻⁵ for V_w = 0.01c; fixture F-WAVE-02).

**Jacobian (Derived).** At fixed p′ and E′,

$$
\frac{\partial\mu}{\partial\mu'}=\gamma_w\frac{p'}{p}\left(1-\frac{\mu V_w}{v}\right),
\qquad
\left(\frac{\partial\mu}{\partial\mu'}\right)^{2}(1-\mu'^2)=\gamma_w^{2}(1-\mu^{2})\left(1-\frac{\mu V_w}{v}\right)^{2},
\tag{63}
$$

where μ, p and v are plasma-frame quantities (fixture F-WAVE-03). A wave-frame coefficient ∝ (1 − μ′²) thus carries the factor (1 − μV_w/v)² in the plasma frame; with ω/k∥ = V_w this is the [1 − μω/(k∥v)]² of Strauss et al. (2017, Eq. 13) and Spanier et al. (2022, Eq. 22). The identity maps one step; it says nothing about rates or the cross terms D_μp.

**Two printed shortcuts.** Schreiner & Spanier (2015) print μ′ = μ − ω/(kv_α) as the form "used in prior work" and state that it is not applicable to dispersive waves (Eq. 2.4). The first-order expansion of Eq. (61) in ε = V_w/v is μ′ = μ − ε(1 − μ²) + O(ε²) (Derived). The two agree only at μ = 0; Eq. (61) keeps μ′ = ±1 at |μ| = 1. Their Eq. 2.10 evaluates γ′ from the Galilean speed.

**Monte Carlo procedures as printed.** The Vainio & Kocharov (2001) kernel and the Sandroos et al. (2014) step are in Section 14.6; for Sandroos et al. the frame-selection rule and transformation formulas were not found. Virtanen & Vainio (2005, App. A; Full text ×1) [P91]: "the particle is first scattered in the forward-wave frame, and then immediately again in the backward-wave frame", with weights ω_± ∈ [0, 1], ω₊ + ω₋ = 1, from the amplified-wave ratios (formula not given); "Scatterings are elastic in the wave frame"; time step printed as "Δt = a ν_e(γ)", a ≲ 0.01 (as extracted).

**Momentum diffusion (only as printed).** For waves streaming parallel and antiparallel at equal intensities, Virtanen & Vainio (2005, Eq. 15; Full text ×2), citing Schlickeiser (1989), give

$$
A_2(x,p)\simeq\frac{2\pi\Omega_e^{2-q}k_0^{q}I_0}{q(q+2)B_0^{2}}\,\frac{p^{2}v_A^{2}}{v^{3-q}},
\tag{64}
$$

with A₂ ∝ p^q for v ≈ c, p ≫ m_ec (Eq. 16; consistent with Eq. 15 via Ω_e = Ω_{e,0}/γ, check S8-9). Their Monte Carlo scheme obtains the energy change from the frame transformations alone. Ng & Reames (2008, Eq. 1) retain D_μP and D_PP (matrix not legible). Chandran (2000) and Yan & Lazarian (2002, 2004) print D_μp, D_pp for anisotropic turbulence (Section 9).

### 8.3 D_μμ in a wave frame

**Printed forms.** Afanasiev et al. (2015, Eqs. 1, 2, 5; Full text ×2, ar5iv and A&A HTML) give

$$
k_{res}=\frac{\Omega}{\upsilon\mu},\qquad
D_{\mu\mu}(\mu)=\frac{\pi}{2}\Omega\frac{|k_{res}|I_{w,res}}{B^{2}}(1-\mu^{2}),\qquad
\delta B^{2}=\int_{-\infty}^{+\infty}I_w(k)\,dk ,
\tag{65}
$$

with Ω = γ⁻¹Ω₀, B the mean field and μ, υ in the wave rest frame (SI). The reading of the denominator (B² or R²) is discussed in Section 13.7. Afanasiev et al. (2023, Eq. 4; 2025, Eq. 14) print (π/2)Ω₀|k_r|I(k_r)/(γB²)(1 − μ²) with k_r = Ω₀/(γvμ).

Borovikov et al. (2019, Eqs. 4.30, 4.31; Full text ×1, SI) give, per branch σ (σ = 1…4 enumerates the signs of V_σ and g_σ = ±1),

$$
D^{\sigma}_{\mu\mu}=\frac{\pi\omega_{ci}^{2}}{2(B^{2}/\mu_0)}(1-\mu_\sigma^{2})\int_0^\infty dk\,I_\sigma(k)\,\delta(|k|\mu_\sigma v_\sigma-g_\sigma\omega_{ci})
=\frac{\pi\omega_{ci}}{2(B^{2}/\mu_0)}(1-\mu_\sigma^{2})\frac{\omega_{ci}}{v_\sigma|\mu_\sigma|}I_\sigma\!\left(\frac{\omega_{ci}}{v_\sigma|\mu_\sigma|}\right),
\tag{66}
$$

which is nonzero for the two branches with g_σ = sign(μ_σ). The normalization is (δB)²/μ₀ = Σ_σ∫₀^∞d|k| I_σ(|k|). Neglecting the difference between μ and μ_σ and taking v_σ ≫ V_A, they obtain the plasma-frame form (Eq. 5.3)

$$
D_{\mu\mu}=\frac{\pi\omega_{ci}}{4(B^{2}/\mu_0)}(1-\mu^{2})\frac{\omega_{ci}}{v|\mu|}\left[I_-\!\left(\frac{\omega_{ci}}{v|\mu|}\right)+I_+\!\left(\frac{\omega_{ci}}{v|\mu|}\right)\right].
\tag{67}
$$

The wave-frame forms of Vainio (2003, Eq. 2), Battarbee et al. (2011, Eqs. 8–9), Vainio & Kocharov (2001, Eq. 1) and Afanasiev et al. (2015, Eq. 4) are tabulated in Section 14.2, with their mutual conversions (equal if W = I_w/(4π), and U = 2kW for Vainio & Kocharov; Decision U-30). Three further forms are:

| Source (verification) | Printed D_μμ | Frame of μ, v | Spectrum as printed | Resonance as printed |
|---|---|---|---|---|
| Vainio, Laitinen & Fichtner 2003 (×1) [P92] | (π/4)Ω(1−μ²)f_rP(f_r,r)/B² | not stated; f fixed-frame | P(f) per unit fixed-frame frequency | f_r = (Ω/2π)V/(v\|μ\|), V = u + v_A |
| Sandroos et al. 2014, Eq. 29 (×1) | (π/2)(Ω/B₀²)(1−μ̃²)\|k_res\|I(r,k_res) | parallel or antiparallel wave rest frame | I(λ)dλ = I(k)dk; δB²/μ₀ = (1/μ₀)∫I(λ)dλ | λ_r = (2π/Ω)Ṽμ̃ (Eq. 31) |
| Ng & Reames 2008, Eq. 6 (×1, preprint; Σ and ∫ lost in extraction) | Σ_σ(v²/(4P²))∫dk R^σ_μμ I_σ | outward-wave frame | ⟨δB·δB⟩ = Σ_σ∫₀^∞dk I_σ (Ng, Reames & Tylka 2003, Eq. 18 [P93]) | k = B/[P(μ − V_σ/v)] |

R^σ_μμ of Ng & Reames (2008) "includes minimal broadening of the quasilinear resonance condition" and is not printed in the Letter. Ng, Reames & Tylka (2003, Eq. 7) show the factors (1 − μ²), (1 − μV/v)² and a δ-function of (kμv − kV ∓ Ω). Their broadened form (App. A) was not readable. The saturation cap is discussed in Section 14.4.

**Conversions (Derived; checks S8-5 to S8-8, fixture F-WAVE-05).** The reference is PARALLEL rev. 1.4, Eq. 26, D = πΩ²(1 − μ²)P_s(k_res)/(4B₀²v|μ|), with P_s the one-sided total transverse variance spectrum.
- *Afanasiev et al. 2015, Eq. 2* equals PARALLEL rev. 1.4 Eq. 26 if I_w(k) = I_w(−k) = P_s(|k|)/2. This is the I_w = 4πg of Section 6.7. Only one sign of k resonates at a given μ. Their Eq. 4 equals Eq. 2 with I_w(k_res(μ))/|μ| replaced by its |μ| = 1 value. For I ∝ k^{−5/3} its ratio to Eq. 2 is 1.842 at μ = 0.4.
- *Borovikov et al. 2019.* Eq. (66), summed over the two resonant branches, gives Eq. (67) if each polarization carries half of I_± (this is the origin of π/4 versus π/2). With P_s = μ₀(I₊ + I₋) it equals PARALLEL rev. 1.4 Eq. 26.
- *Vainio, Laitinen & Fichtner 2003* equals PARALLEL rev. 1.4 Eq. 26 with P(f) = (2π/V)P_s(2πf/V). Battarbee et al. (2011, Eq. 2), ν = π²f_ci f_res P(f_res)/B² with f_ci = qB/(2πm_ic) and f_res = f_ci(u_sw + v_A)/(γv), equals it at |μ| = 1 only for γ = 1. As printed it is larger by the factor γ, because f_ci carries no γ while Ω is relativistic (S8-7).
- *Resonances.* Ng & Reames' k = B/[P(μ − V_σ/v)] is k(vμ − V_σ) = Ω with Ω = QeB/(γmc). This is the argument of the Borovikov δ-function written with v_σμ_σ = vμ − V_σ (S8-8). The sign of k_r differs between Vainio (2003) and Battarbee et al. (2011) (−ω_c/(vμ)) and Afanasiev et al. (+Ω/(γvμ)) (Section 14.2).

**Behaviour at μ = 0 as printed.** Afanasiev et al. (2015) state that QLT "does not work close to μ = 0" and fill the gap linearly (Section 14.6). Sandroos et al. (2014) use μ̃_min = 0.025. Vainio (2003) cites Vainio (2000) for "scattering near μ = 0, which can not be due to outward-propagating Alfvén waves in quasilinear theory". Ng, Reames & Tylka (2003) use a resonance-broadening parameter 0.10 (symbol not legible).

### 8.4 Dispersive waves and electrons

**Spanier, Schreiner & Schlickeiser (2022)** [P80] (arXiv numbering; Eqs. 18, 20–23 Full text ×2, others ×1). Cold electron–proton plasma, parallel waves:

$$
\frac{k_\parallel^2c^2}{\omega_{L,R}^2}=1-\frac{\omega_{pi}^2}{\omega(\omega\mp\Omega_i)}-\frac{\omega_{pe}^2}{\omega(\omega\pm\Omega_e)},\qquad
k^2=-\frac{(M+1)y_{L,R}^2}{(y_{L,R}\mp1)(y_{L,R}\pm M)},
\tag{68}
$$

(Eqs. 3 and 9; upper sign L), with y = ω/Ω_i, k = k∥/k_c, k_c = Ω_i/v_A (Eq. 8), M = m_p/m_e = 1836. The R branches are ω_R = V_A|k∥| (|k∥| < k_c), Ω_i + V_A²k∥²/Ω_i (k_c ≤ |k∥| ≤ M^{1/2}k_c) and Ω_e(1 − Mk_c²/k∥²) (beyond) (Eq. 16), interpolated as y_R ≈ |k|(1 + |k|) (Eq. 17). The resonance is

$$
v\mu k_\parallel-\omega_{R,L}(k_\parallel)\mp\frac{Q\Omega}{\gamma}=0,\qquad
\Omega\left(\frac{\mu k}{\epsilon}-y_{R,L}(k)\mp S_i\right)=0,\qquad
S_i=\frac{Q|Z_i|m_p}{m\gamma}=\frac{1}{\gamma}\begin{cases}1&\text{protons}\\-M&\text{electrons}\end{cases}
\tag{69}
$$

(Eqs. 18, 20, 21), with ε = v_A/v (Eq. 19), Ω = |q|B₀/(mc) and Ω′ = QΩ/γ. As printed, the prefactor of Eq. 22 below carries Ω without γ while the resonance uses Ω′; relative to a coefficient with the relativistic gyrofrequency in the prefactor (Eq. (43)), this is a factor γ² (derived; D-54). The coefficient (Eq. 22, Full text ×2) is

$$
D_{\mu\mu}(\mu)=\frac{\pi^2\Omega^2(1-\mu^2)}{B_0^2}\int_{-\infty}^{\infty}dk_\parallel\left(g_R\,\delta(vk_\parallel-\omega_R-\Omega')\left(1-\frac{\mu\omega_R}{k_\parallel v}\right)^2+g_L\,\delta(vk_\parallel-\omega_L+\Omega')\left(1-\frac{\mu\omega_L}{k_\parallel v}\right)^2\right).
\tag{70}
$$

It is printed with δ(vk∥ − …) without μ in both renderings, whereas the dimensionless Eq. 24 contains kμ/ε. This is a printed inconsistency (Section 19, D-54). The helicity-resolved four-term form is unnumbered in the extraction. In it, the L-terms' δ-functions contain y_R as extracted (flag). For backward whistlers (H_R = −1, σ = −1) the paper prints (Eq. 38)

$$
D_{\mu\mu}=\frac{\pi^2\Omega^2k_c(1-\mu^2)}{\Omega_i^2B_0^2}\,\Theta[k_r-k_0]\,\Theta[M^{1/2}-k_r]\,\frac{g_{tot}(k_r)}{|dy_R/dk+\mu/\epsilon|_{k_r}}\left(1+\frac{\epsilon\mu y_R(k_r)}{k_r}\right)^2 ,
\tag{71}
$$

with roots k_{1,2} = ½(±(μ²/ε² − 4S_i)^{1/2} − μ/ε), forward root k_r = ½((μ²/ε² − 4S_i)^{1/2} + μ/ε) and reality condition μ² ≥ 4S_iε² (Eqs. 33–34, 39). Derived from Eq. (70) with the Jacobian of the δ-function (v_A = Ω_i/k_c), the prefactor is π²Ω²k_c/(Ω_iB₀²): the printed Ω_i² makes D_μμ dimensionless unless it is expressed in units of Ω_i, which the text read does not state (D-127). The interaction limits are μ_R(x) = β_A(M − (1 + x²)^{1/2})/x, μ_L(x) = β_A(M + (1 + x²)^{1/2})/x and μ₀(x) = β_AM^{1/2}((1 + x²)^{1/2} − 1)/x, with x = p/(m_ec) and β_A = v_A/c. The whistler interaction range is −μ₀ ≤ μ ≤ μ_R (Eqs. 43–46). The text says the model fails at μ = 0 and for 10 MeV electrons; cyclotron-damped ranges are excluded and damping broadening is ignored. The PIC comparison parameters (runs S1/S2, fitted s = 3.1, σ = 0) are in the companion data (Section 24).

*Derived checks (S8-10, fixture F-WAVE-06).*
- The dimensionless Eq. 9 is Eq. 3 with the vacuum term 1 dropped, using ω_pe² = Mω_pi² and c²/v_A² = ω_pi²/Ω_i². The relative difference is 1.1 × 10⁻⁷ at y = 10 for c²/v_A² = 10⁸.
- The printed roots solve k² ± kμ/ε + S_i = 0, i.e. the resonance with y_R = k².
- The backward-wave resonance kμ/ε + y_R(k) + S_i = 0 (the form whose k-derivative is the denominator of Eq. 38), with S_i = −M/γ, gives μ = μ_R at k = 1 (y_R = 1) and μ = −μ₀ at k = M^{1/2} (y_R = k² = M). μ_L is the L-type resonance kμ/ε − y_L + S_i = 0 at k = y_L = 1.
- The interpolation Eq. 17 is approximate: at k = 3, Eq. 9 gives y_R = 9.855, Eq. 17 gives 12 and the Eq. 16 branch gives 10.

**Schreiner & Spanier (2015)** (ICRC; Full text ×1). The resonance of order n is ω − k∥μv_α = nΩ_α/γ_α (Eq. 2.1). The frame transformation is Eq. (61) (their Eq. 2.9). The dispersion relation is

$$
|k^\pm|=\frac{\omega}{c}\sqrt{1-\frac{\omega_p^2}{(\omega\pm\Omega_e)(\omega\pm\Omega_p)}}
\tag{72}
$$

(Eq. 2.11), with k⁺ the whistler/R branch, k⁻ the Alfvén/L branch and Ω_e "negative by definition". Derived (S8-11, fixture F-WAVE-07): with charge neutrality, Eq. 2.11 is the exact two-species cold-plasma relation (Eq. (68) with the vacuum term) if ω_p² = ω_pe² + ω_pp²; the notes do not record the source's definition of ω_p. Parameters: B₀ = 0.20 G; ω_p = 6.4 × 10⁷ (model) and 1.5 × 10⁸ rad s⁻¹ (simulation); ∼100 keV electrons. Since each resonant dispersive wave has its own phase speed, a dispersive spectrum has no single wave frame; the paper prescribes none for a whole spectrum.

**Strauss & le Roux (2019): four modes with dispersion.** The thermal dissipation-range onsets k_d^{LH,RH} (their Eq. 15) are in Section 7.7. The cold-plasma relations (App. B, PDF only) are k²_{∥,LH} ≈ (ω²/V_A²)Ω_p/(Ω_p − ω) (B15) and k²_{∥,RH} ≈ (ω²/c²){|Ω_e|/(|Ω_e| − ω) + ω_e²/[(|Ω_e| − ω)(Ω_p + ω)]} (B16). Derived (S8-12): (B15) is the y ≪ M limit of the L branch of Eq. (68). (B16) differs from its R branch by relative O(1/M) and O(V_A²/c²) (5.4 × 10⁻⁴ for V_A²/c² = 10⁻⁸).

### 8.5 What an implementer needs

| Model ID (equation) | Frame of μ, v | Conserved in one scattering | Plasma-frame change | Spectrum input | Fully specified? |
|---|---|---|---|---|---|
| MEAS-AFANASIEV15 (65; Section 13.7) | rest frame of outward waves | υ (wave frame) | Eq. (60); not printed | I_w two-sided, sign of k = polarization | yes, except gap filler b, b₁ |
| SGW-COUPLE-VK01 (Section 14.6) | each ± wave frame | p_± | Eq. (59) | U_± per ln k (U-30) | yes, isotropic ν_± |
| WAVE-TWOFRAME-VV05 | forward, then backward wave frame | speed in each frame ("elastic") | Lorentz transformation | weights ω_± (formula not given) | no |
| WAVE-BOR19-BRANCH (66), WAVE-BOR19-PF (67) | branch frame v_σ = v − V_σ**b**; Eq. 5.3 plasma frame | v_σ | not printed | I_σ per branch, SI | yes; Eq. 5.3 approximate |
| WAVE-SL19-DT (56), WAVE-SL19-ND-PRINTED (57) | plasma frame (one reading) | not stated | not stated | g_{n,j}, 2π one-sided | no (U-26; sign of Ω; r_L) |
| WAVE-VLF03 | not stated; f fixed-frame | not stated | not stated | P(f) per unit fixed-frame frequency | yes, except the frame of μ (not stated; declared by the configuration) |
| WAVE-NR08 | outward-wave frame | not stated | D_μP, D_PP retained (not legible) | I_σ one-sided per mode | no (R^σ_μμ not printed) |
| WAVE-SPA22-BW (71), WAVE-SPA22-GEN (70) | μ relative to B₀; frame not stated | not stated | factors (1 ∓ μω/(kv))² | g_tot, H_{L,R}, σ (WAVE-SPLIT-SPA22) | Eq. 38 yes; general form has printed issues |
| WAVE-XFORM-LORENTZ (60), WAVE-XFORM-GAL (61), WAVE-DE (62) | — | p′, E′ | exact | — | yes (derived kernels) |

**Implementation requirement.** Every wave-based D_μμ model declares the frame of μ and v, and the pusher applies the coefficient in that frame. Afanasiev et al. (2015, 2023, 2025), Vainio (2003), Ng & Reames (2008) and Borovikov et al. (2019) print wave-frame variables; applying their D_μμ to plasma-frame μ is a different model.

**Implementation requirement.** A wave-frame scattering keeps the particle speed in that frame (Vainio & Kocharov 2001; Virtanen & Vainio 2005). The plasma-frame energy change is then Eq. (62). A pusher that transforms explicitly must not also add Eq. (64) or another D_pp for the same waves; that would count the energy change twice (follows from Eq. (62)).

**Implementation requirement.** With forward and backward waves, each direction's D_μμ acts in its own frame. Vainio & Kocharov (2001) scatter in both frames every step; Virtanen & Vainio (2005) scatter forward then backward; the Sandroos et al. (2014) rule was not found. Order and weights are configuration, not inference.

**Implementation requirement.** A wave spectrum carries direction, polarization, sidedness and normalization. The conversions of Section 8.3 hold only under their stated equal-polarization or equipartition assumptions; without them the adapter refuses rather than rescales (PARALLEL rev. 1.4, Sections 8.7 and 12.4). For Strauss & le Roux (2019) the four-component split is an input (U-26), and the use of Eq. 28 or of the exact-Jacobian δ-limit is recorded. A model that applies one frame to a dispersive spectrum says so, since Schreiner & Spanier (2015) state that the single-frame shift Eq. 2.4 does not apply to dispersive waves.

## 9. Anisotropic Alfvénic and compressible (fast-mode) turbulence

The sources of this section replace the slab spectra of Sections 6 and 8 by three-dimensional spectra: Alfvénic and pseudo-Alfvén/slow turbulence with Goldreich–Sridhar (GS95) anisotropy, and isotropic fast modes. The coefficients contain gyroharmonic sums with Bessel functions of w = k⊥v⊥/Ω and an n = 0 transit-time-damping (TTD) term. μ (or ξ = cosθ) is measured relative to **B**₀ (relative to the local field in the Xu & Yan 2013 simulations, Section 9.5), the operator is the standard Fokker–Planck form (with D_μp, D_pp where printed), and the resonance is written in the plasma frame. The printed parameter sets are interstellar or in code units (Section 9.6).

### 9.1 Goldreich–Sridhar Alfvénic turbulence: suppressed gyroresonance

**Chandran (2000)** [P84], PRL 85, 4656 (Full text ×1, arXiv preprint). The paper gives QLT pitch-angle (ξ = cosθ) and momentum coefficients for shear-Alfvén (s) and pseudo-Alfvén (p) modes. The general integrals are its Eqs. 6–7. The spectrum (Eqs. 8–9) is

$$
M^s(k_\perp,k_\parallel,\tau)=\frac{L^3}{6\pi}k_\perp^{-10/3}l^{-1/3}\,g\!\left(\frac{k_\parallel}{k_\perp^{2/3}l^{-1/3}}\right)e^{-|\tau|/\tau_k},\qquad
\tau_k=\frac{l}{v_A}(k_\perp l)^{-2/3},
\tag{73}
$$

with g(x) = 1 for |x| < 1 and 0 otherwise. Also M^p = M^s = K^s = K^p and C^s = C^p = σM^s, with σ the cross helicity. The inline normalization is ∫_{l⁻¹}^∞k⊥dk⊥∫dk∥M^s(k⊥,k∥,0)B₀²/4 = L³B₀²/(8π). The lower limit l⁻¹ here differs from l′⁻¹ in the validity statement (as printed). For sinθ ≫ ε^{1/2} (Eqs. 11–12),

$$
D^s_{\xi\xi}=\frac{v}{l}\left[\frac{2\epsilon^{3/2}|\cos\theta|^{11/2}}{13\sin\theta}\sum_{n=1}^\infty n^{-9/2}-\frac{\delta\ln\epsilon}{3}\sin^2\theta\right],\quad
D^p_{\xi\xi}=\frac{v}{l}\left\{\frac{2\epsilon^{3/2}|\cos\theta|^{7/2}\sin\theta}{13}\sum n^{-9/2}-\frac{\delta\ln\epsilon}{6}\frac{\sin^4\theta}{\cos^2\theta}\left[1-\frac{v_A}{v_\parallel}\arctan\frac{v_\parallel}{v_A}\right]\right\},
\tag{74}
$$

with ε = v/(lΩ), δ = v_A/v and scattering frequency ν = 2D_ξξ/(1 − ξ²). D_ξp and D_pp are obtained by multiplying by −σpδ and by p²δ². For sinθ ≪ ε^{1/2}, D^s_ξξ = (v/l)(πθ²/8) and ν → (π/4)(v/l) (Eq. 13). The parallel diffusion coefficient is

$$
\kappa_\parallel=\frac{v^2}{8}\int_{-1}^{1}d\xi\,\frac{(1-\xi^2)^2}{D_{\xi\xi}},\qquad
\kappa_\parallel^{QLT}=vl(-\delta\ln\epsilon)^{-1}\left(\frac52-\frac{3\pi}{4}\right),\qquad
\kappa_\parallel^{QLT}=vl\,c_1(-\delta\ln\epsilon)^{-5/11}\epsilon^{-9/11},\ c_1\simeq0.88
\tag{75}
$$

(Eqs. 15–17). Eq. 16 holds for ε^{3/2} ≪ (−lnε)δ; the condition for Eq. 17 is not recorded in the readings. The text says that the n = 0 terms dominate at low energy, that QLT is unreliable in the trapping regime and that Eq. 16 is a lower bound on κ. Parameters (ε = 2.2 × 10⁻⁹E_GeV for protons; interstellar B₀, l, v_A) are in the companion data (Section 24).

*Derived checks (S9-1 to S9-4, fixture F-CMP-01).*
- The printed normalization integral equals L³B₀²/(8π), and Eq. 13 gives ν → (π/4)v/l.
- Eq. 16 is Eq. 15 with only the two n = 0 terms of Eq. (74) (arctan bracket set to 1): (1/8)∫₋₁¹6ξ²(1 − ξ²)/(1 + ξ²)dξ = 5/2 − 3π/4.
- Eq. 17 results from the pseudo-Alfvén gyroresonant term (∝ |ξ|^{7/2}) competing with the pseudo-Alfvén n = 0 term (∝ ξ⁻²) near ξ = 0: c₁ = ¼[(2π/11)/sin(6π/11)](2ζ(9/2)/13)^{−6/11}6^{5/11} = 0.87835, matching the printed 0.88. Quadrature of the reduced two-term integral (1/8)∫₋₁¹ξ²/(A|ξ|^{11/2} + B)dξ, which omits the factors (1 − ξ²)² and sin θ, approaches it (ratio 0.99737 at ε = 10⁻², δ = 10⁻⁹; 0.99992 at ε = 10⁻³, δ = 10⁻¹⁴); with those factors kept the ratios are lower (about 0.963 and 0.997 for the two pseudo-Alfvén terms; derived).

**Yan & Lazarian (2002, 2004) and Xu & Lazarian (2020): Alfvénic gyroresonance.** Yan & Lazarian (2002) [P85] (PRL; Full text ×1) use ∂_tf = ∂_μ(D_μμ∂_μf + D_μp∂_pf) + p⁻²∂_p[p²(D_μp∂_μf + D_pp∂_pf)], tensors normalized by B₀², V_AB₀ and V_A² (Eq. 1), and K_ij = C_aI_ijk⊥^{−10/3}exp(−L^{1/3}k∥/k⊥^{2/3}), C_a = L^{−1/3}/(6π), M_ij = K_ij, C_ij = σM_ij (Eq. 3). For a pitch angle α "not close to 0" (Eq. 8):

$$
[D_{\mu\mu},D_{\mu p},D_{pp}]=\frac{v^{2.5}\cos\alpha^{5.5}}{2\Omega^{1.5}L^{2.5}\sin\alpha}\,\Gamma[6.5,k_{max}^{-2/3}k_{res}L^{1/3}]\,[1,\sigma mV_A,m^2V_A^2],
\tag{76}
$$

with k_res = Ω/v∥, k_max = Ω₀/v_th, k_min = L⁻¹, m = γm_H; Γ[·,·] is not defined in the readings. Yan & Lazarian (2004, Eq. 12) [P94] print v^{2.5}μ^{5.5}/(Ω^{1.5}L^{2.5}(1 − μ²)^{0.5})Γ[6.5, k_max^{−2/3}k∥,res L^{1/3}][1, m²V_A²], without the factor ½ (as printed). Xu & Lazarian (2020) [P95] (Full text ×1) use I_A = C_Ak⊥^{−10/3}exp(−L^{1/3}k∥/k⊥^{2/3}), C_A = δB_A²L^{−1/3}/(6π), normalized to δB²/2 (Eqs. 5–6; L → LM_A⁻³ if super-Alfvénic), and give

$$
D_{\mu\mu,QLT,A}\approx\frac23\,8^{13/2}e^{-8}\,\frac{\delta B_A^2}{B_0^2}\left(\frac{v}{L\Omega}\right)^{3/2}\frac{v}{L}(1-\mu^2)^{-1/2}\mu^{11/2}
\tag{77}
$$

(Eq. 10; prefactor 165.820). For slow modes: (2/3)8^{13/2}e^{−8}(δB_s²/B₀²)(v/(LΩ))^{3/2}(v/L)(1 − μ²)^{1/2}μ^{7/2} (Eq. 17); the slow spectrum exponent is printed as k∥/k∥^{2/3} (XL2020 and XL2018, Eq. 13; as printed).

**Domain of the printed asymptotes (derived).** As printed, Eq. (77) and the slow-mode and fast-mode forms contain μ^{11/2}, μ^{7/2} and μ^{1/2}, which are not real for μ < 0, and Eq. (77) contains (1 − μ²)^{−1/2}, which diverges at μ = ±1; Yan & Lazarian (2002) state Eq. (76) for a pitch angle "not close to 0". They are asymptotic estimates for 0 < μ < 1 away from the endpoints, not full-domain coefficients: an endpoint divergence of D_μμ is not a scattering law, and extending an asymptote to μ = ±1 or joining asymptotes across regimes is not supported by the sources. A runtime evaluation over −1 ≤ μ ≤ 1 requires an explicit symmetry rule (|μ| in place of μ, which the sources do not print) and an endpoint treatment, both configured as U-66; outside the configured domain the coefficient returns a domain diagnostic (registry field `domain`, Section 4.1). The same applies to the regime limits of the Chandran (2000) forms and to the fast-mode halo and TTD forms of Sections 9.2–9.3. Derived (S9-6): with k⊥ ≥ L⁻¹, C = δB²L^{−1/3}/(6π) gives ∫d³k I = δB²/2 only over k∥ ≥ 0, and δB² over all k∥; the domain is not recorded (U-34).

### 9.2 Fast modes: gyroresonance and transit-time damping

*Spectra.* Yan & Lazarian (2002, Eqs. 4–5): low β, [M, C, K]_ij = (L^{−1/2}/(8π))J_ijk^{−7/2}[cos²θ, σcosθ, 1] with J_ij = k_ik_j/k⊥² and cosθ = k∥/k; high β, (L^{−1/2}/(8π))sin²θJ_ijk^{−7/2}[cos²θ/β, σcosθ/β^{1/2}, 1], whose prefactor Yan & Lazarian (2004, Eq. 10) print as 1/(2π). Xu & Lazarian (2018) [P96] and (2020) use I_f = C_fk^{−7/2}, C_f = δB_f²L^{−1/2}/(16π), with δB_f²/2 = ∫d³k I_f (XL2018 Eqs. 29–31; XL2020 Eqs. 19–20); this holds over the full sphere k ≥ L⁻¹, as does XL2018 App. B, I_f = C_fk^{−11/3}, C_f = δB_f²L^{−2/3}/(12π) (Derived, S9-6). ω_k = V_phk (XL2018 Eq. 34) and ω_tur,f = (V_Lf²/V_ph)L^{−1/2}k^{1/2} (Eq. 39). Yan & Lazarian (2004, Eq. 3) print τ_k = (k/L)^{−1/2}V_ph/V² with V_ph "equal to Alfven and sound velocity for high and low beta plasma, respectively" (apparently inverted; as printed).

*Gyroresonance.* The general form for slow and fast modes is D_μμ = C_μ∫d³k(k∥²/k²)[J₁′(x)]²I(**k**)R(**k**), with C_μ = (1 − μ²)Ω²/B₀² and x = k⊥v⊥/Ω (XL2018 Eqs. 2–3, 73; XL2020 Eq. 12). In the QLT limit Xu & Lazarian (2020) give

$$
D_{\mu\mu,QLT,f}\approx\frac{\pi}{56}\frac{\delta B_f^2}{B_0^2}\left(\frac{v}{L\Omega}\right)^{1/2}\Omega(1-\mu^2)\mu^{1/2}
\tag{78}
$$

(Eq. 22). Xu & Lazarian (2018, Eq. 75), (2π²/7)(C_f/B₀²)Ω^{1/2}v^{1/2}(1 − μ²)μ^{1/2}, is identical to it with their Eq. 30 (Derived, S9-7). For the Galactic halo, Yan & Lazarian (2004, Eq. 14) print D_μμ = (π(Ωvμ)^{0.5}(1 − μ²)/(2L^{0.5}))(1/7)(1 − [(k⊥,c/k∥,res)² + 1]^{−7/4}) (quasi-parallel; quasi-perpendicular addition Eq. 15), with k⊥,c set by damping (Yan & Lazarian 2002, Eq. 9; 2004, App. A).

*TTD.* In the linear (δ-function) limit, Xu & Lazarian (2018, Eqs. 36–37) give

$$
D_{\mu\mu,fL,T}=2H^*\pi^2\frac{C_f}{B_0^2}\left(l_{min}^{-1/2}-L^{-1/2}\right)v\,\frac{(1-\mu^2)^2}{\mu}\left(\frac{V_{ph}}{v_\parallel}\right)^2\left[1-\left(\frac{V_{ph}}{v_\parallel}\right)^2\right],
\tag{79}
$$

with H* step functions (Eq. 37; their form was not transcribed) and D ∝ (1 − μ²)²/μ³ for v∥ ≫ V_ph (Eq. 38). Yan & Lazarian (2004) evaluate TTD only numerically, with the Breit–Wigner factor τ_k/(τ_k² + (k∥v∥ − ω)²) (App. B, PDF only). The linear resonance is at v∥ = V_ph; XL2018 Eq. 5 states the condition v ≥ V_ph.

### 9.3 Resonance broadening by magnetic mirroring (Yan & Lazarian 2008)

**Source.** Yan & Lazarian (2008) [P97], ApJ 673, 942. Eqs. 1–9 Full text ×2 (ar5iv and arXiv PDF), except Eq. 9 (below). Eqs. 10–22 Full text ×1 (PDF OCR, partly garbled). The paper broadens the resonance with a Gaussian, from the spread of v∥ caused by δB∥. The spread is printed as

$$
\frac{\Delta v_\parallel}{v_\perp}=\frac{\langle(B-B_0)^2\rangle^{1/4}}{B_0^{1/2}}\simeq\left[\frac{\langle\delta B_\parallel^2\rangle}{B_0^2}+o\!\left(\frac{\langle(\delta B_\perp)^2\rangle^2}{B_0^4}\right)\right]^{1/4}
$$

(Eq. 1; the first form equals (⟨(B − B₀)²⟩/B₀²)^{1/4} and is dimensionless, although one reading flagged it as dimensionally odd), with σ_z² = v⊥²(⟨δB∥²⟩/B₀²)^{1/2}t² (Eq. 4). The coefficient (Eq. 5) is

$$
D_{\mu\mu}=\frac{\Omega^2(1-\mu^2)}{B_0^2}\int d^3k\sum_{n=0}^{\infty}R_n(k_\parallel v_\parallel-\omega\pm n\Omega)\left[I^A(\mathbf k)\frac{n^2J_n^2(w)}{w^2}+\frac{k_\parallel^2}{k^2}J_n'^2(w)I^M(\mathbf k)\right],
\tag{80}
$$

with w = k⊥v⊥/Ω = k⊥LR(1 − μ²)^{1/2}. R is printed as "R = v/(Ωl) is the dimensionless rigidity of the CRs, L is the injection scale" (l versus L as printed). I^M is the magnetosonic spectrum, "in our case at hand … the pseudo-Alfven modes". The resonance function is printed as

$$
R_n=\Re\int_0^\infty dt\,e^{i(k_\parallel v_\parallel+n\Omega-\omega)t-\frac12k_\parallel^2v_\perp^2t^2(\langle\delta B_\parallel^2\rangle/B_0^2)^{1/2}}
=\frac{\sqrt\pi}{|k_\parallel\Delta v_\parallel|}\exp\!\left[-\frac{(k_\parallel v\mu-\omega+n\Omega)^2}{k_\parallel^2\Delta v_\parallel^2}\right]
\simeq\frac{\sqrt\pi}{|k_\parallel|v_\perp M_A^{1/2}}\exp\!\left[-\frac{(k_\parallel v\mu-\omega+n\Omega)^2}{k_\parallel^2v^2(1-\mu^2)M_A}\right]
\tag{81}
$$

(Eq. 6), with M_A = δV/v_A = δB/B₀ and Δμ = Δv∥/v. The GS95 spectrum is I^A = I^S = (L^{−1/3}M_A^{4/3}/(6π))exp(−L^{1/3}|k∥|/(M_A^{4/3}k⊥^{2/3})) (Eq. 7); the factor k⊥^{−10/3} appears in the integrands. The TTD coefficient of pseudo-Alfvén modes in strong turbulence (Eq. 8) is

$$
D^T_{\mu\mu}=\frac{v\sqrt\pi(1-\mu^2)}{6\pi LR^2}\int_1^{k_{max}L}d^3x\,\frac{x_\parallel x_\perp^{-10/3}M_A^{4/3}}{(x_\perp^2+x_\parallel^2)\Delta\mu}J_1^2(w)\exp\!\left[-\frac{x_\parallel}{x_\perp^{2/3}M_A^{4/3}}-\frac{(\mu-v_A/v)^2}{\Delta\mu^2}\right],\qquad \mathbf x=\mathbf kL .
\tag{82}
$$

For small rigidity (Eq. 9) the two renderings disagree:

$$
D^T_{\mu\mu}\approx\frac{\sqrt\pi M_A^{7/2}v}{16L}(1-\mu^2)^{3/2}\left[-E_1(qx_\parallel)-e^{-qx_\parallel}\right]_1^{x_{\parallel,max}}\exp\!\left[-\frac{(\mu-v_A/v)^2}{\Delta\mu^2}\right]\quad\text{(ar5iv)},
\tag{83}
$$

with q = (x⊥,max M_A²)^{−2/3} and E₁(x) = ∫₁^∞dt e^{−xt}/t. The arXiv PDF text layer (OCR) gives the bracket as [−E₁(qx_k) − e^{−qx_k}/x_{k,max}]. The two copies differ, and the difference is not resolved from the source.

*Derived (S9-8 to S9-10, fixture F-CMP-03).*
- **Eq. 6 is internally inconsistent.** The integral, with exponent −½k∥²Δv∥²t², equals (π/2)^{1/2}/(|k∥|Δv∥)·exp[−ω̃²/(2k∥²Δv∥²)] (ω̃ the mismatch); the printed closed form is narrower by √2. Both have area π. Xu & Lazarian (2018, Eq. 7) print the self-consistent form. As printed, Eqs. 6, 8, 9 use Δμ² = (1 − μ²)M_A; the Xu & Lazarian convention corresponds to 2Δμ² (U-36).
- **Eq. 8** follows from Eqs. 5–7 for n = 0 with Ω²L/v = v/(R²L).
- **Small-R reduction of Eq. 8** (J₁² ≈ w²/4, x∥ ≪ x⊥, d³x = 2πx⊥dx⊥dx∥ on [1, max], Δμ = (1 − μ²)^{1/2}M_A^{1/2}) gives the ar5iv bracket with exactly the printed q, M_A^{7/2} and (1 − μ²)^{3/2}, plus a boundary term (3/(2c²))[E₁(cx∥) + e^{−cx∥}]₁^{x∥,max} (c = M_A^{−4/3}) from x⊥ = 1, and the prefactor √π/8 instead of the printed √π/16. It supports the ar5iv bracket but does not settle the factor 2 or the boundary term; the domain of Eq. 8 is not printed precisely.

The remaining forms are Full text ×1 and partly garbled. Eq. 10 (gyroresonance, strong turbulence) has an uncertain x∥ placement. Eq. 13 (weak-turbulence TTD, R ≪ 1) reads D^T_μμ = (v√π/(8L))(1 − μ²)^{3/2}M_A^{3/2}[tan⁻¹(M_A⁻²) − π/4]exp[−(μ − v_A/v)²/((1 − μ²)M_A)]. Eqs. 14–15 (fast modes) contain Δμ_k, undefined in the extracted text; Eq. 16 contains an undefined r. The fast-mode cutoff is k_cL = [4M_A⁴m_i cos²θ/(πm_eβsin⁴θ)]exp[2m_e/(βm_i cos²θ)] (Eqs. 17–18, prefactor unclear), and λ∥/L = (3/4)∫₀¹dμ v(1 − μ²)²/[(D^T_μμ + D^G_μμ)L] (Eq. 22). The text states that TTD covers 90°, with the resonance peak at α → 90° at k∥,res ∼ Ω/Δv. Table 1 (halo, WIM) is in the companion data (Section 24); its WIM β = 0.077 conflicts with β ≃ 0.11 in the text (as printed).

### 9.4 Resonance-broadened TTD and mirror trapping (Xu & Lazarian 2018, 2020)

**Xu & Lazarian (2018)**, ApJ 868, 36 (Full text ×1). The TTD coefficient is D_μμ,T = C_μ∫d³k(k∥²/k²)[J₀′(x)]²I(**k**)R(**k**) (Eqs. 1–3), with linear resonance R_L = πδ(ω_k − v∥k∥) (Eq. 4) and spread Δv∥ ≈ v⊥(⟨δB∥²⟩/B₀²)^{1/4} (Eq. 6). The broadened resonance function is

$$
R_B=\Re\int_0^\infty dt\,e^{i(\omega_k-v_\parallel k_\parallel)t-\frac12(\mathbf k\cdot\delta\mathbf r(t))^2}=\frac{\sqrt{2\pi}}{2(\Delta v_\parallel k_\parallel+\omega_{tur})}\exp\!\left[-\frac{(\omega_k-v_\parallel k_\parallel)^2}{2(\Delta v_\parallel k_\parallel+\omega_{tur})^2}\right]
\tag{84}
$$

(Eq. 7; limits Eqs. 9 and 11). For slow modes, with C_s = δB_s²L^{−1/3}/(6π) (Eq. 14) and ω_k = V_phk∥ (Eq. 18),

$$
D_{\mu\mu,sL,T}=\frac{\pi^2}{2}\frac{C_s}{B_0^2}L^{-2/3}\ln\frac{L}{l_{\perp,min}}\,v^2(1-\mu^2)^2\,\delta(v_\parallel-V_{ph}),
\quad
D_{\mu\mu,s1,T}\approx\frac{\sqrt2}{4}\pi^{3/2}\frac{C_s}{B_0^2}\left(\frac{\langle\delta B_\parallel^2\rangle}{B_0^2}\right)^{-1/4}L^{-2/3}\ln\frac{L}{l_{\perp,min}}\,v(1-\mu^2)^{3/2}\exp\!\left[-\frac{(V_{ph}-v_\parallel)^2}{2\Delta v_\parallel^2}\right]
\tag{85}
$$

(Eqs. 19 and 25). The text says D_μμ,s1,T is approximately constant at small μ, i.e. finite at 90°. For fast modes the broadened form is (Eq. 42)

$$
D_{\mu\mu,f1,T}\approx\frac{\sqrt2}{4}\pi^{3/2}\frac{C_f}{B_0^2}\left(\frac{\langle\delta B_\parallel^2\rangle}{B_0^2}\right)^{-1/4}\left(l_{min}^{-1/2}-L^{-1/2}\right)v(1-\mu^2)^{3/2}\exp\!\left[-\frac{v_\parallel^2}{2\Delta v_\parallel^2}\right],
\tag{86}
$$

and Eq. 45 is the decorrelation-dominated fast TTD. *Derived (S9-11, fixture F-CMP-04):* Eq. 25 is Eq. 19 with δ(v∥ − V_ph) replaced by a normalized Gaussian of standard deviation Δv∥ = v(1 − μ²)^{1/2}(⟨δB∥²⟩/B₀²)^{1/4}. The interstellar parameters (Fig. 1; Tables 2–3) are in the companion data (Section 24).

**Xu & Lazarian (2020)**, ApJ 894, 63 (Full text ×1). Gyroresonance is as in Sections 9.1–9.2, with R_L = πδ(ω_k − v∥k∥ + Ω) and the wave frequency "negligible for relativistic particles" (Eq. 4). The loss cone is μ_lc² = b_k/(B₀ + b_k) ≈ b_k/B₀ (Eqs. 23–24). The trapping rates are Γ_t,f = (v/(2L))(δB_f⁴/B₀⁴)(1 − μ²)/μ⁷ (Eq. 29) and Γ_t,s = (v/(2L))(δB_s²/B₀²)(1 − μ²)/μ³ (Eq. 36). The scattering rate is Γ_s = (1/μ²)⟨(Δμ)²⟩/Δt = 2D_μμ/μ² (Eq. 40). The cutoff below which trapping dominates, and the resulting parallel diffusion coefficient, are

$$
\mu_{cf,f}\approx\left[\frac{14}{\pi}\frac{\delta B_f^2}{B_0^2}\left(\frac{v}{L\Omega}\right)^{1/2}\right]^{2/11},\qquad
D_{\parallel,f,t}=\frac{v^2}{4}\int_{\mu_{cf,f}}^1d\mu\,\frac{(1-\mu^2)^2}{D_{\mu\mu,QLT,f}}=\frac{28}{5\pi}\frac{B_0^2}{\delta B_f^2}\left(\frac{v}{L\Omega}\right)^{-1/2}\frac{v^2}{\Omega}\left[4-\sqrt{\mu_{cf,f}}\,(5-\mu_{cf,f}^2)\right]
\tag{87}
$$

(Eqs. 42–44). *Derived (S9-12, fixture F-CMP-02):* μ_cf,f is the μ at which 2D_μμ,QLT,f/μ² (Eq. (78)) equals Γ_t,f. The closed form of Eq. 44 equals the integral of Eq. 43. The text says mirror trapping "removes the 90° singularity" in the parallel diffusion coefficient. No D_μμ value at μ = 0 is printed. Eq. 46 has B_s² in a denominator (flagged as printed). Appendix (A1)–(A4) repeats the broadened slow-mode TTD with the kernel of Eq. (84) (ω_tur = 0).

### 9.5 Simulation tests

**Xu & Yan (2013)** [P98], ApJ 779, 140 (Full text ×1; Eq. 5 Full text ×2). Test particles (1000 per run, Bulirsch–Stoer) move in 512³ isothermal compressible MHD turbulence with injection scale 0.4 of the box, M_A = 0.19–0.73 (ten values) and 1.5, and r_L = 0.01 and 0.03 of the box. No analytic D_μμ is printed. The estimator is D_μμ = ⟨(μ − μ₀)²⟩/(2t) (Eq. 5; the "short time interval" is not given numerically), and λ∥/L = (3/4)∫₀¹dμ₀u(1 − μ₀²)²/(D_μμL) (Eq. 6). The pitch angle is "the angle between the particle's velocity vector and the local magnetic field direction", i.e. a local-field μ. As printed: no resonance gap at 90° (the small gap "disappears because of the resonance broadening"); mirror interaction with compressible modes dominates at large pitch angles; D⊥ ∝ M_A^{4.21±0.75} (Eq. 12). Table values are in the companion data (`simulation_benchmarks.json`, SIM-XUYAN13).

**Lynn, Parrish, Quataert & Chandran (2012)** [P99], ApJ 758, 78 (Full text ×1). The paper calibrates TTD resonance broadening against test-particle simulations in slow/Alfvénic GS turbulence. It gives a diffusion coefficient in v∥, not D_μμ, and its μ is not the pitch-angle cosine. The resonance functions are R(**k**) = Re∫₀^∞dt e^{i(ω(**k**)−v∥k∥)t} (Eq. 4) and

$$
R_{exp}=\frac{\gamma v_A(k_\perp L)^{2/3}(2\pi)^{1/3}}{\gamma^2v_A^2(2\pi)^{2/3}(k_\perp L)^{4/3}/L+k_\parallel^2L(v_\parallel\pm v_p)^2},\qquad
R_{gauss}=\frac{\pi^{1/6}L}{\gamma v_A(4k_\perp L)^{2/3}}\exp\!\left[-\frac{k_\parallel^2L^2(v_\parallel\pm v_p)^2}{4\gamma^2v_A^2(2\pi)^{2/3}(k_\perp L)^{4/3}}\right]
\tag{88}
$$

(Eqs. 8 and 12), with v_p = 0.81c_s for β = 1 slow modes. The δ-limit (Eq. 16) is D∥ = (π²/(24L))(v_L/v_A)²v⊥⁴log(L/L_min)δ(v∥ ± v_p). The calibrated forms are Eqs. 10 and 13, with constants C₁, C₂ and γ = 1 (exponential) or 0.5 (Gaussian) in Table 2 (companion data); Eq. 13 uses 1/u₊³ in both terms (as printed).

### 9.6 Applicability statements and implementation rules

All parameter sets printed by the Section 9 sources, as recorded in the readings, are interstellar (ISM phases, Galactic halo, TeV–PeV cosmic rays) or in code units (Xu & Yan 2013). None of them, in the parts read, applies the coefficients to the inner heliosphere or to SEPs (Not found). The only heliospheric statement on fast modes found here is Vainio (2000) [P86], Abstract only: in a mixture of slab-mode and isotropic fast-mode waves, agreement with observations requires |H_c| ∼ 1.

| Model ID | Equation | Fully specified? |
|---|---|---|
| CMP-CH00-DXX, CMP-CH00-KAPPA16, CMP-CH00-KAPPA17 | (74) with Eq. 13; (75) | yes (asymptotic forms) |
| CMP-YL02-ALF, CMP-YL04-ALF-PRINTED | (76); Yan & Lazarian 2004 Eq. 12 | no (Γ[·,·] undefined in the readings; factor ½ differs) |
| CMP-XL20-ALF, CMP-XL20-SLOW, CMP-XL20-FAST-G | (77); XL2020 Eq. 17; (78) | yes |
| CMP-YL04-FAST-HALO | Yan & Lazarian 2004 Eq. 14 | no (k⊥,c from damping) |
| CMP-XL18-SLOW-TTD, CMP-XL18-FAST-TTD-B | (85), (86) | yes |
| CMP-XL18-FAST-TTD-LIN | (79) | no (H* not transcribed) |
| CMP-XL20-TRAP | (87) | yes |
| CMP-YL08-RES, CMP-YL08-TTD, CMP-YL08-TTD-SMALL | (81), (82), (83) | no (U-36; domain of Eq. 8; Eq. 9 renderings disagree) |
| CMP-LYNN12-REXP, CMP-LYNN12-RGAUSS | (88), for D in v∥ | yes, as kernels |

**Implementation requirement.** μ is relative to **B**₀ in the analytic forms and relative to the local field in Xu & Yan (2013); a benchmark comparison uses one definition or reports the mismatch.

**Implementation requirement.** A Gaussian resonance width is stored with its convention: exp(−x²/Δ²) (Yan & Lazarian 2008, Eq. 6) or exp(−x²/(2Δ²)) (Xu & Lazarian 2018, Eq. 7); they differ by √2 and each paper's closed forms use their own. PARALLEL rev. 1.4, Section 9.4 defines Δ through C_k(t) = exp(−Δ²t²/4).

**Implementation requirement.** Each spectrum declares its normalization (B₀²-normalized tensors in Yan & Lazarian 2002, 2004; δB²/2 = ∫d³k I in Xu & Lazarian 2018, 2020) and its k∥ domain (U-34). The integral forms (80), (82) and (84) are tabulated in (μ, p) with the harmonic truncation stated (Chandran 2000 checked with |n| ≤ 10).

## 10. Nonlinear and second-order theories of D_μμ

Quasilinear theory (QLT, Section 6) assumes unperturbed helical orbits. For magnetostatic slab turbulence with inertial-range index s > 1 its sharp resonance gives D_μμ ∝ |μ|^{s−1}, so D_μμ(0) = 0. The theories in this section broaden the resonance and give a finite D_μμ at μ = 0. Dynamical-turbulence resonance functions (damping, random sweeping, waves) are in Section 7. Fast-mode and transit-time-damping forms are in Section 9. PARALLEL rev. 1.4, Sections 9–11, defines generic broadened slab kernels (its Eqs. 38–42) and their λ∥ fixtures (its Section 15.8). Those kernels are explicitly not the paper-specific theories reviewed here.

Notation. Unless stated otherwise, the operator is Eq. (21), μ is the cosine of the pitch angle relative to **B**₀, and the turbulence is magnetostatic, so the plasma and fixed frames coincide. R_L = v/Ω = pc/(qB₀) (Gaussian units) and R = R_L/l_slab [P73, Eq. 4]. This review writes **b ≡ δB/B₀** as a shorthand. In Shalchi et al. (2009) δB is the turbulence strength of a pure slab model; whether it is the slab-only rms is not stated in the text read.

### 10.1 Slab spectrum and general form (Shalchi et al. 2009)

Source: Shalchi, Škoda, Tautz & Schlickeiser (2009) [P73].
- Eqs. 1, 2, 10–15, 29, 32 and 49 are legible in the aanda.org alt-text (Full text ×1; Eqs. 10–12 Full text ×2, see below).
- Eqs. 9, 20–24, 33, 35 and 45–47 are complete only in the raw text layer of the A&A PDF (Full text ×1 per equation). Eqs. 10–12 agree between the two renderings.

The general slab form is Eq. 6. With the QLT resonance function R_± = πδ(vμk ± Ω) (Eq. 7), Eq. 6 gives Eq. 8: D_μμ = [2π²v(1−μ²)/(|μ|B₀²R_L²)]g^slab[k∥ = (|μ|R_L)⁻¹].

$$
D_{\mu\mu}=\frac{2\pi v^2(1-\mu^2)}{B_0^2R_L^2}\int_0^\infty dk\,g^{\rm slab}(k)\left[R_+(k)+R_-(k)\right]
\tag{89}
$$

$$
g^{\rm slab}(k_\parallel)=\frac{D(s,q)}{2\pi}\,\delta B^2\,l_{\rm slab}\,\frac{|k_\parallel l_{\rm slab}|^q}{[1+(k_\parallel l_{\rm slab})^2]^{(s+q)/2}},\qquad
D(s,q)=\frac{\Gamma\!\left(\frac{s+q}{2}\right)}{2\Gamma\!\left(\frac{s-1}{2}\right)\Gamma\!\left(\frac{q+1}{2}\right)}
\tag{90}
$$

Eq. (90) is Eqs. 34–35. The paper says the spectrum is "correctly normalized for q > −1 and s > 1".

The following were derived in this review and checked with mpmath:
- **Normalization.** The normalization is 8π∫₀^∞ g^slab dk = δB².
- **Conversion to PARALLEL.** The one-sided total-variance spectrum of PARALLEL rev. 1.4 (its Eq. 25) is P_s = 8πg^slab. With this conversion, Eq. 8 is PARALLEL Eq. 26 and Eq. (43).
- **Equality of constants.** D(s, 0) equals C(s) = Γ(s/2)/[2√πΓ((s−1)/2)] [P77, Eq. 4]. It also equals C(ν) with ν = s/2 [P76, Eq. 5].

### 10.2 Second-order resonance function (Shalchi 2005, via Shalchi et al. 2009)

Shalchi et al. (2009, Sect. 3.1) "employ the second order theory of Shalchi (2005)" [P100]. Footnote 2 says the resonance function was obtained by combining a second-order approach with a late-time approximation and a 90° approximation. Shalchi (2005) itself is not accessible; it is known secondarily through Shalchi et al. (2009) and Subashchandar et al. (2025).

Eq. 9, Eq. 10 and Eqs. 11–12 as printed:

$$
R^{\rm SOQLT}_\pm(k)=\frac{\sqrt\pi}{vk}\frac{B_0}{\delta B}\exp\!\left[-\frac{(v\mu k\pm\Omega)^2B_0^2}{(vk\,\delta B)^2}\right]
=\frac{\sqrt\pi}{\sigma}\,e^{-\beta_\pm^2/\sigma^2},\quad \sigma=\frac{vk\,\delta B}{B_0},\quad \beta_\pm=v\mu k\pm\Omega
\tag{91}
$$

- **Area.** ∫dβ_± R_± = π (Eq. 13), the same area as the QLT δ-function.
- **Step-function approximation.** The closed forms below rest on the approximation of Eq. 15:

$$
R^{\rm SOQLT}_\pm(k)\approx\frac{\pi}{2\sigma}\,H(\sigma-|\beta_\pm|)
\tag{92}
$$

Eqs. 20–24 carry Heaviside functions, i.e. they use Eq. (92). With Γ = B₀/(R_LδB) (Eq. 23) and k_± = R_L⁻¹(|μ| ± b)⁻¹ (consistent with Eqs. 27 and 38), Eq. 24 gives D_μμ = [2πv²(1−μ²)/(B₀²R_L²)]·I. Here:
- non-resonant regime: I = (π/2Ω)Γ[∫_{−k₋}^∞ + ∫_{k₊}^∞](dk/k)g^slab for |μ| < b;
- quasi-resonant regime: I = (π/2Ω)Γ∫_{k₊}^{k₋}(dk/k)g^slab for |μ| > b.

Eq. 21 contains a k₀ that is not defined in the lines read.

### 10.3 Closed forms, μ = 0, and the strong-turbulence limit

Eq. 33 gives the value at μ = 0. Eq. 47 is the closed form for the spectrum Eq. (90) in the inertial-range approximation used in Sect. 5 of the paper.

$$
D=D_{\mu\mu}(\mu=0)=\frac{2\pi^2\Omega}{B_0^2}\,\Gamma\int_\Gamma^\infty\frac{dk}{k}\,g^{\rm slab}(k)
\tag{93}
$$

$$
D_{\mu\mu}=\frac{\pi D(s,q)}{2s}(1-\mu^2)\frac{v}{l_{\rm slab}}R^{s-2}\frac{\delta B}{B_0}\sum_{n=\pm1}{\rm sign}\!\left(\frac{\delta B}{B_0}+n|\mu|\right)\left||\mu|+n\frac{\delta B}{B_0}\right|^s
\tag{94}
$$

Eq. 47 was read once, in the PDF raw layer. That reading agrees with an earlier extraction that had been marked "reconstructed".

For δB ≫ B₀ the coefficient is isotropic, D_μμ = (1−μ²)D (Eq. 32). In that limit the paper states that "the previously used assumption of isotropic pitch-angle diffusion is indeed correct for sufficiently strong turbulence". The strong-turbulence mean free path (Eq. 49, alt-text) is

$$
\lambda_\parallel=\frac{s}{2\pi D(s,q)}\,l_{\rm slab}R^{2-s}\left(\frac{B_0}{\delta B}\right)^{s+1}
\tag{95}
$$

The weak-turbulence result (Eq. 51) was read once, in a single HTML fetch:

$$
\lambda_\parallel=\frac{3\,l_{\rm slab}R^{2-s}}{2\pi D(s,q)}\left[\frac{1}{(2-s)(4-s)}\left(\frac{B_0}{\delta B}\right)^2+\frac{s(2-s)-1}{2(2-s)}\left(\frac{B_0}{\delta B}\right)^s\right]
\tag{96}
$$

Validity as stated: |μ| > b is the quasi-resonant (QLT-like) regime and |μ| < b the non-resonant regime. Waves enter through |μ| → |μ − jε|, with ε = v_A/v and j = ±1 (Sect. 7). Published inputs: R = 0.0363 and b = 0.05 (Fig. 3); s = 5/3 (Fig. 5).

Each of the following was derived in this review and checked with mpmath (derivation record DR-10):
1. **Inertial-range evaluation.** Eq. (94) is the evaluation of Eq. 24 with g^slab ∝ k^{−s}, i.e. Eq. 24 in the inertial range. The paper's own Eq. 46 is the intermediate step. The approximation therefore requires k_± l_slab ≫ 1, i.e. R(|μ| + b) ≪ 1.
2. **QLT limit.** For b ≪ |μ|, Eq. (94) tends to the QLT inertial form πD(s,q)(1−μ²)(v/l_slab)R^{s−2}b²|μ|^{s−1}. The ratio is 0.9999999985 at μ = 0.5, b = 10⁻⁴.
3. **Value at μ = 0.** At μ = 0, Eq. (94) equals Eq. (93) with the inertial spectrum:
   - D = (πD(s,q)/s)(v/l_slab)R^{s−2}b^{s+1}.
   - Eq. (23) then gives λ∥ = v/(2D) for D_μμ = (1−μ²)D, which is Eq. (95).
   - This is Eq. (24) with ν = 2D.
4. **Bendover spectrum (q = 0).** For q = 0 and the full bendover spectrum, Eq. (93) has the closed form below. The Gaussian resonance (Eq. (91) in Eq. (89)) gives the second line instead (U is Tricomi's function).

$$
\frac{D^{\rm step}l_{\rm slab}}{v}=\frac{\pi D(s,0)}{s}b^{s+1}R^{s-2}\,{}_2F_1\!\left(\tfrac s2,\tfrac s2;1+\tfrac s2;-(Rb)^2\right),\qquad
\frac{D^{\rm Gauss}l_{\rm slab}}{v}=\sqrt\pi\,D(s,0)\,\Gamma\!\left(\tfrac s2\right)\frac{b}{R^2}\,U\!\left(\tfrac s2,1,(Rb)^{-2}\right)
\tag{97}
$$

   In the inertial range, Rb → 0, the ratio Gauss/step tends to sΓ(s/2)/√π. That is 1.06141647 for s = 5/3. At R = b = 1 the ratio is 0.96384. The two resonance forms are therefore not interchangeable at μ = 0.
5. **Mean free path from Eq. 47.** The small-b expansion of λ∥ from Eq. (94) inserted in Eq. (23) has the leading term of Eq. (96). Its (B₀/δB)^s coefficient is −0.7737038288 for s = 5/3, not −2/3. The printed coefficient (s(2−s)−1)/(2(2−s)) is reproduced exactly by a two-regime approximation: D(0)(1−μ²) for |μ| < b and QLT for |μ| > b. Whether the paper derived Eq. 51 this way was not read.

   At R = 0.0363, b = 0.05, s = 5/3 and q = 0, λ∥/l_slab is (F-NL-03):
   - 532.3186 from Eq. 47 through Eq. (23);
   - 553.3579 from Eq. (96);
   - 684.0225 from the QLT inertial form;
   - 2177.744 from Eq. (95).

### 10.4 SOQLT with time-dependent orbit broadening (Subashchandar et al. 2025)

Source: Subashchandar et al. (2025) [P101]. Full text ×2: HTML and PDF agree on main-text Eqs. 1–3 and appendix Eqs. 32–40; the appendix has its own numbering. The paper applies SOQLT to Parker Solar Probe spectra. Main text, Eq. 1:

$$
D^{(2)}_{\mu\mu}=\frac{4\pi\Omega^2(1-\mu^2)}{B_0^2}\int_0^\infty dk_\parallel\,g^{\rm slab}(k_\parallel)\int_0^\infty dt\,\cos(\Omega t)\,e^{ik_\parallel v\mu t}\,e^{-k_\parallel^2\sigma_z^2(t)/2}
\tag{98}
$$

App. Eqs. 32 and 34, with β_± = k∥vμ ± Ω:

$$
\sigma_z^2(t)=\frac{4\pi\Omega^2v^2(1-\mu^2)}{B_0^2}\int_0^\infty dk\,g(k)\left[M_+(t)+M_-(t)\right],\qquad
M_\pm(t)=\frac{1-\cos\beta_\pm t}{\beta_\pm^4}-\frac{\sin(\beta_\pm t)}{\beta_\pm^3}\,t+\frac{t^2}{2\beta_\pm^2}
\tag{99}
$$

Background:
- **Guiding-centre distribution.** It is assumed Gaussian: Γ = e^{ik⟨z⟩}e^{−k²σ_z²/2} with ⟨z⟩ = vμt (App. Eqs. 26 and 29).
- **Slab tensor.** P_ij = g(k∥)δ(k⊥)/k⊥ with ⟨δB_iδB_j⟩ = ½⟨δB²⟩δ_ij, so ⟨δB_xδB_x*⟩ = 4π∫g dk (App. Eq. 22).
- **Long-time limit at μ = 0.** App. Eq. 36 prints σ_z² → (v²/2)b²t².
- **Closed form near μ = 0.** It is App. Eq. 40, with K_± = ½Re∫₀^∞e^{i(±Ω+kvμ)t}e^{−k²σ_z²/2}dt (App. Eq. 38):

$$
K^{(2)}_\pm=\frac12\,\frac{\sqrt\pi}{vk\,(\delta B/B_0)}\exp\!\left[-\left(\frac{\mu\pm\Omega/(kv)}{\delta B/B_0}\right)^2\right]
\tag{100}
$$

Printed inconsistencies:
- App. Eq. 35 prints lim M₊ + M₋ = 1/(Ωt)² at μ = 0 and t ≫ 1/Ω. That is the same in both renderings.
- App. Eqs. 12–13 carry Ω signs that the text extraction shows as inconsistent with Eq. 14.
- The spectral indices are written "−q = 1.5" and "−5/3" while q^S = 3/2 and q^2D = 5/3, so the sign convention is inconsistent.
- κ⊥/κ∥ is given as "10⁻⁴–10⁻³", as "three orders" and as "three to four orders" in different places.

The following were derived and checked in this review:
- **Dominant term at μ = 0.** It is t²/Ω²: (M₊ + M₋)Ω²/t² = 0.99835 at Ωt = 10³. With 8π∫g dk = δB² this reproduces App. Eq. 36, whereas the printed Eq. 35 does not.
- **Gaussian integral.** App. Eq. 40 is the exact Gaussian integral of App. Eq. 39.
- **Relation to Shalchi et al. (2009).** Eq. (100) equals ½ × Eq. 11 of Shalchi et al. (2009), i.e. K = R/2. With cos(Ωt)e^{ik∥vμt} = ½Σ_±e^{i(kvμ±Ω)t}, the real part of Eq. (98) (printed without Re; App. Eq. 38 has Re) is then identical to Eq. (89). The two papers therefore use the same normalization; their QLT limits coincide with Eq. 8.

Inputs and results:
- **Spectrum input.** The measured-spectrum input (main Eq. 4) is given in Section 13. The conversion from P_S to g^slab is not stated in the text read.
- **Mean free path.** λ∥ uses main Eq. 2, κ∥ = (v²/8)∫(1−μ²)²/D^(2)_μμ dμ, which is Eq. (23).
- **Parameter ranges.** δB/B₀ ≈ 0.2–0.6; 10–60 R_⊙; protons 500 keV–1 GeV; a QLT comparison cut-off μ_min = 0.05 or 0.01.
- **Fitted results.** κ∥ ∝ r^{1.1} (500 keV) to r^{1.5} (1 GeV) and κ∥ ∝ (δB/B₀)^{−2.13}. These are also listed in MFP v1.2 Sections 4.5 and 7.6.

### 10.5 Other nonlinear broadenings (legible parts only)

| Theory | Resonance function as printed | Width | Status here |
|---|---|---|---|
| Weakly nonlinear theory (WNLT) [P102] | not read | perpendicular-motion broadening; per the book chapter abstract it "mainly increases D_μμ" and shortens λ∥ [P103] | abstract only; no formula. PARALLEL rev. 1.4, Section 9.1, states that its kernels are not WNLT |
| NADT [P79] | Eq. 22: R_n^slab = γ_slab/[γ²_slab + (k∥μv + nΩ − v_Ak∥)²]; 2D part Eq. 33: R_n^2D = γ_2D/[γ²_2D + (nΩ)²] | γ_slab = εZ/l_2D (Eq. 8, ε "of order 1"); Z = √2 v_AδB_2D/B₀ (Eq. 12) | OCR with Greek glyphs lost; structure only (the squared width in Eq. 22 is legible at structure level). Not implementable as printed. Parameters: MFP v1.2 Section 7.6 |
| Nonlinear damping [P76] | Eq. 19: v_A√(k∥/L)/[v_A²k∥/L + (k∥v∥ + nΩ − v_Ak∥)²] | γ_c = v_A√(k∥/L) (Eq. 7) | Full text ×2 for Eqs. 4–9, 12–13. Eq. 20 is printed with "Σ_{n=+1}", while Eqs. 17 and 21 use n = ±1. Dynamical forms: Section 7 |
| Yan & Lazarian [P97] | Eq. (101) below | Δv∥ from the magnetic-mirror (Völk) estimate, Eq. 1 | Eqs. 1–9 Full text ×2; Sect. 3 onward PDF only. Fast-mode and TTD forms: Section 9 |
| SOQLT, isotropic turbulence [P104, App. A.6b] | not legible | ξ "= v²(δB/B₀)^{2/3}" as extracted, dimensionally suspect | reconstructed by the text extraction; not used |

Yan & Lazarian (2008), Eq. 6, in the plasma frame (ω is the wave frequency, M_A = δV/v_A = δB/B₀):

$$
R_n=\frac{\sqrt\pi}{|k_\parallel\Delta v_\parallel|}\exp\!\left[-\frac{(k_\parallel v\mu-\omega+n\Omega)^2}{k_\parallel^2\Delta v_\parallel^2}\right]\simeq\frac{\sqrt\pi}{|k_\parallel|v_\perp M_A^{1/2}}\exp\!\left[-\frac{(k_\parallel v\mu-\omega+n\Omega)^2}{k_\parallel^2v^2(1-\mu^2)M_A}\right]
\tag{101}
$$

Details:
- **Gaussian form.** Eq. (101) has the same Gaussian form as Eq. (91). The width differs:
  - Yan & Lazarian: σ = |k∥|v⊥M_A^{1/2}, from compressive δB∥ (Eq. 4: σ_z² = v⊥²(⟨δB∥²⟩/B₀²)^{1/2}t²);
  - Shalchi et al. (2009): σ = vkδB/B₀, from transverse slab fluctuations.
- **Width in μ.** For ω = 0 the width in μ is Δμ² = (1−μ²)M_A. This is the reading of Eq. 6 recorded for this review; it is consistent with the text's sharp-resonance limit μ ≫ Δμ.
- **Printed anomalies.** The bracket of Eq. 9 differs between the ar5iv and PDF renderings. The plasma β of the warm ionized medium (WIM) is 0.077 in Table 1 and ≃ 0.11 in the text.

### 10.6 What the theories predict at μ = 0

| Source | Statement on D_μμ(0), as printed or as derived from the printed formula |
|---|---|
| QLT, magnetostatic slab (Section 6) | D_μμ ∝ \|μ\|^{s−1} → 0 for 1 < s < 2 in the inertial range |
| Shalchi et al. 2009, Eqs. 32–33 | finite, Eq. (93); the isotropic regime holds "for sufficiently strong turbulence" |
| Shalchi et al. 2009, Eq. 47 at μ = 0 | (πD(s,q)/s)(v/l_slab)R^{s−2}(δB/B₀)^{s+1} (derived, item 3 of Section 10.3) |
| Subashchandar et al. 2025, App. (40) | finite Gaussian resonance at μ → 0 |
| Yan & Lazarian 2008 | finite: the Gaussian width Δμ stays nonzero at μ = 0. The text states that TTD covers 90°, with the resonance peak at α → 90° at k∥,res ~ Ω/Δv |
| NADT; nonlinear damping | Lorentzian of nonzero width: nonzero at μ = 0 (derived from the printed structure) |
| Tautz et al. 2013a, simulations [P104] | isotropic turbulence: neither QLT nor SOQLT describes 90° scattering well. SOQLT overestimates D_μμ(0), which "compensates" values that are too low near μ ≈ 0, so λ∥ agrees |
| Qin & Shalchi 2009 (secondary via [P74]) | SOQLT reproduces simulated 90° scattering for weak slab turbulence [P105] |
| van den Berg et al. 2024, Case 3 | strong slab turbulence: isotropic for \|μ\| > 0.7, more scattering near μ = 0 [P106] |
| Dundovic et al. 2020 | QLT agreement at lag 10⁻²λ∥/c; at 10⁻¹λ∥/c a discrepancy near μ ~ 0, "attributed to resonance broadening" [P107] |
| Xu & Yan 2013, MHD simulations | no resonance gap at 90°, attributed to resonance broadening and mirror interaction [P98] |

**Implementation requirement.** A SOQLT D_μμ object records which resonance function it uses:
- Gaussian, Eq. (91);
- step, Eq. (92);
- the inertial closed form, Eq. (94).

At μ = 0 the Gaussian-to-step ratio is sΓ(s/2)/√π (1.0614 for s = 5/3) in the inertial limit Rb → 0, and it varies between about 0.956 (Rb ≈ 3) and 2/√π = 1.128 (Rb → ∞) outside it (derived from Eq. (97) and its step analogue; F-NL-05), so none may be substituted for another silently. Eq. (94) is evaluated only with its spectrum Eq. (90) in the inertial-range domain R(|μ| + δB/B₀) ≪ 1, and a request outside it returns a domain diagnostic. The spectrum is supplied in the 8π convention of Eq. (90), or converted by P_s = 8πg^slab.

**Implementation requirement.** NL-SOQLT-SUB25 evaluates σ_z²(t) from App. Eqs. 32 and 34. It does not use App. Eq. 35 as printed.

## 11. D_μμ from test-particle simulations and how it is estimated

This section collects the published estimators of D_μμ from particle trajectories, their numerical settings, the published simulated D_μμ usable as benchmarks, and the regimes where the estimators fail. It also states which estimators are consistent with Eq. (22). That is the property needed to validate a particle pusher against a prescribed D_μμ.

### 11.1 Estimators as printed

Labels M0–M5 are those of van den Berg, Els & Engelbrecht (2024) [P106]. That paper was read in Full text ×2 (arXiv HTML and PDF agree on all equations cited here) up to Sect. 4.4; its conclusions were not accessible.

| Label | Printed estimator | Sources |
|---|---|---|
| M0 | D_μμ(μ) = ∫₋₂² [(Δμ)²/(2Δt)] P(Δμ\|μ, Δt) d(Δμ) = ⟨(Δμ)²/(2Δt)⟩ | [P106, Eqs. 4a–4b] |
| M1a | D_μμ(μ(0), t) = ⟨[μ(t) − μ(0)]²/(2t)⟩ | [P106, Eq. 6]; [P104, Eq. 3]; [P108, Eq. 6, printed without ⟨ ⟩ and with lim t → ∞]; [P109, Eq. 24]; [P98, Eq. 5]; [P110] (with lim t → ∞) |
| M1b | D_μμ(μ(0), t) = d/dt⟨[μ(t) − μ(0)]²/2⟩ | [P106, Eq. 7]; [P111, Eq. 3] |
| M2a | D_μμ(μ(0)) = ∫₀^∞⟨μ̇(0)μ̇(τ)⟩dτ (Taylor–Green–Kubo) | [P106, Eq. 8]; [P104, Eq. 2a]; [P111, Eq. 2] |
| M2b | T_μ(μ(0)) = [1/C_μμ(μ(0), 0)]∫₀^∞C_μμ(μ(0), τ)dτ with C_μμ = ⟨μ(0)μ(τ)⟩, then D_μμ = (1−μ²)/(2T_μ) | [P106, Eqs. 11, 12, 14] (new method) |
| M3 | D_μμ(μ(0), t) = ⟨μ̇(t)[μ(t) − μ(0)]⟩ | [P106, Eq. 15]; ⟨μ̇(0)Δμ(t)⟩ in [P112, Eq. 34] |
| M4a / MIIb | D_μμ(μ, t) = [∫₋₁^μ ∂f(η, t)/∂t dη][∂f(μ, t)/∂μ]⁻¹ | [P106, Eq. 19]; [P108, Eq. 11] (= −j_μ); [P111, Eq. 24] |
| M4b / MIIa | fit of ∂_tf = D′∂_μf + D∂²_μf on a μ grid (central differences, tridiagonal system) | [P106, Eq. 20b]; [P108, Eqs. 9–10] (Savitzky–Golay derivatives) |
| M5 | stationary injection at μ_S with absorbing boundaries μ_L, μ_R, Eq. (102) | [P106, Eqs. 21–22, 26–27], after Kaiser et al. (1978) and Sun et al. (2016) [P113] |
| corr. | D_μμ = −½(1−μ²) d/dt ln⟨μ(t)μ(0)⟩ | [P111, Eq. 28]: arXiv PDF only; the A&A HTML numbering ends at Eq. 27 |
| R19 | D_μμ,MSD = ⟨(Δμ)²/Δt⟩, **no factor ½** | [P45, Eq. 11], one rendering |
| D20 | "histogram of Δμ with respect to μ(t)", no formula printed | [P107] |

The M5 relations as printed (Eqs. 21–22):

$$
0=\frac{d}{d\mu}\!\left[D_{\mu\mu}\frac{dF}{d\mu}\right]+\dot N_0\,\delta(\mu-\mu_S),\qquad
D_{\mu\mu}(\mu)=\begin{cases}\dot N_L\,(dF/d\mu)^{-1}, & \mu_L<\mu<\mu_S\\ -\dot N_R\,(dF/d\mu)^{-1}, & \mu_S<\mu<\mu_R\end{cases}
\tag{102}
$$

Symbols:
- F(μ_L) = F(μ_R) = 0.
- Ṅ_{L,R} are the escape rates at the two boundaries.
- F̂(μ_m) = Ṅ₀N_mδt/(N_pϖ_μ) (Eq. 26) and N̂_{L/R} = Ṅ₀N_{L/R}/N_p (Eq. 27).
- N_m counts particles in bin m after every step; N_p is the number injected; ϖ_μ is the bin width.
- Ṅ₀ cancels out of Eq. 22.

The paper says the Sun et al. (2016) form has a sign typo; that form is known only via van den Berg et al. (2024). The signs of Eq. (102) agree with the flux direction: the flux is −D_μμdF/dμ = −Ṅ_L left of μ_S and +Ṅ_R right of it (derived).

**Consistency with Eq. (22) and Eq. (21).** Each of the following was derived in this review and checked with mpmath (derivation record DR-10):
- **M1a.** At finite t, M1a is Eq. (22) evaluated at Δt = t, conditional on μ(0). For a differentiable (full-orbit) trajectory Δμ = μ̇Δt + O(Δt²), so ⟨(Δμ)²⟩/(2Δt) → 0 as Δt → 0 (the ballistic regime of Section 11.5); the limit of Eq. (22) is taken over a lag longer than the correlation and gyration times and shorter than the pitch-angle relaxation time (the admissibility window of Tautz et al. 2013a, Section 11.6), and a running estimator of an orbit code is not the finite-lag reference for a Markov pusher with a prescribed D_μμ, whose reference is the transition solution of its own operator (Section 16.9).
- **M1b and M3.** M1b has the same Δt → 0 limit. M3 equals M1b identically, because d(Δμ)²/dt = 2μ̇Δμ.
- **Inverting the operator.** M4a/MIIb, M4b/MIIa, Tautz (2013b) Eq. 24 and M5 invert the operator Eq. (21) directly. They are consistent with it by construction whenever the trajectories obey that operator.
- **M2b.** M2b reproduces D_μμ exactly only for D_μμ ∝ (1−μ²) (Section 11.3).
- **R19.** As printed, R19 is twice Eq. (22).
- **M0 (Eq. 4b).** M0 pools the pairs (μ, Δt). van den Berg et al. (2024) report that it gives a spurious "isotropic" shape and values that are too large in their composite case.

### 11.2 Isotropic-scattering reference and the bias of the running estimators (Tautz 2013b)

Source: Tautz (2013b) [P111]. Eqs. 1–27 were read; Eqs. 22–27 are confirmed in the A&A HTML and the arXiv PDF.
- **Claim.** The paper shows that the mean-square-displacement definition of D_μμ is flawed for bounded μ because time homogeneity fails.
- **Recommendation.** It recommends a definition based on the pitch-angle correlation function.

For D_μμ = (1−μ²)D (Eq. 8) and f(μ, 0) = δ(μ − μ₀), the solution is f = ½ + Σ_{k≥1}(k + ½)P_k(μ₀)P_k(μ)e^{−k(k+1)Dt} (Eq. 11). It gives ⟨μ⟩(t) = μ₀e^{−2Dt} (Eq. 14) and

$$
\langle(\Delta\mu)^2\rangle=\tfrac13\left(1-e^{-6Dt}\right)+\mu_0^2\left(1+e^{-6Dt}\right)-2\mu_0^2e^{-2Dt}
\tag{103}
$$

$$
D_{\mu\mu}=\tfrac12\frac{d}{dt}\langle(\Delta\mu)^2\rangle=(1-3\mu_0^2)D\,e^{-6Dt}+2\mu_0^2D\,e^{-2Dt}
\tag{104}
$$

Eq. (103) is Tautz Eq. 18, which prints the last exponent as e^{−2dt} (two renderings; read here as e^{−2Dt}, which follows from Eq. 11 and agrees with Eqs. 14 and 22; D-122), and Eq. (104) is Tautz Eq. 22. The limits:
- long times: ⅓ + μ₀² (Eq. 19);
- average over μ₀: ⅔(1 − e^{−2Dt}) (Eq. 20);
- short times: ⟨(Δμ)²⟩ = 2(1−μ₀²)Dt + 𝒪(t²) (Eq. 21).

Eq. 23, D_μμ ≈ (1−μ₀²)D, is printed with the condition "t ≪ D". The intended reading appears to be t ≪ 1/D; it is recorded as printed. The velocity correlation is ⟨v_z(t)v_z(0)⟩ = (v²/3)exp[−2∫₀^tD(t′)dt′] (Eq. 25), and ∫₀^tD dt′ = −½ln(3⟨v_z(t)v_z(0)⟩/v²) (Eq. 26).

Derived in this review and checked with mpmath:
- **Recovery of the printed equations.** Eqs. 14, 18, 20 and 22 follow from Eq. 11 by quadrature to 20 digits.
- **Short-time expansions.** For the running estimators they are

$$
\frac{\langle(\Delta\mu)^2\rangle}{2t}=(1-\mu_0^2)D+(7\mu_0^2-3)D^2t+\mathcal O(t^2),\qquad
\tfrac12\frac{d}{dt}\langle(\Delta\mu)^2\rangle=(1-\mu_0^2)D+(14\mu_0^2-6)D^2t+\mathcal O(t^2).
\tag{105}
$$

  The relative bias of M1a is (7μ₀² − 3)Dt/(1−μ₀²). It changes sign at μ₀² = 3/7, and at fixed t it grows without bound as |μ₀| → 1. Fixture values:
  - At μ₀ = 0.5, M1a/[(1−μ₀²)D] = 0.98357, 0.85489 and 0.34364 for Dt = 0.01, 0.1 and 1.
  - At μ₀ = 0.9 the same ratio is 1.13664, 2.06815 and 2.43493.
- **Correlation estimator.** For this D_μμ the averaged correlation is ⟨μ(t)μ(0)⟩ = ⅓e^{−2Dt}. The estimator "corr." therefore returns (1−μ²)D at every t, without bias.

### 11.3 Exact relations for M2b and M5 under Eq. (21) (derived)

**M2b.** Suppose μ obeys Eq. (25) with an even D_μμ and no focusing. Then u(μ₀) = ∫₀^∞⟨μ(τ)|μ₀⟩dτ solves the backward equation (D_μμu′)′ = −μ with zero flux at ±1 and odd u. Hence D_μμu′ = (1−μ²)/2, and with C_μμ(μ₀, 0) = μ₀² the exact relation and its inversion are

$$
T_\mu(\mu_0)=\frac{1}{\mu_0}\int_0^{\mu_0}\frac{1-\eta^2}{2D_{\mu\mu}(\eta)}\,d\eta,\qquad
D_{\mu\mu}(\mu)=\frac{1-\mu^2}{2\,d[\mu T_\mu(\mu)]/d\mu}.
\tag{106}
$$

Comparison with the printed Eq. 14:
- **Where they agree.** Eq. 14, T_μ = (1−μ²)/(2D_μμ), coincides with Eq. (106) only when (1−μ²)/D_μμ is constant.
- **QLT shape.** For D₀(1−μ²)|μ|^{q−1} the exact T_μ is larger than Eq. 14 by the factor 1/(2−q), which is 3 for q = 5/3.
- **Regularized shape.** For D₀(1−μ²)(|μ|^{2/3} + 0.05) the factor is 1.6986, 2.0556 and 2.1776 at μ₀ = 0.1, 0.5 and 0.9.
- **Independent check.** A finite-volume solution of the backward equation agrees to 0.1 %.

The paper reports that Eq. 14 "had to be multiplied by 2 to match QLT" in its weak-turbulence cases (Sects. 4.1–4.2) but not in the strong case (Sect. 4.3). In the strong case the simulated coefficient is close to isotropic for |μ| > 0.7. The derived relation does not resolve the printed factor 2.

**M5 with isotropic D_μμ.** For D_μμ = D₀(1−μ²) the stationary problem Eq. (102) has the solution F = Ṅ_L[artanh μ − artanh μ_L]/D₀ for μ < μ_S, and correspondingly on the right. The escape fraction is Ṅ_L/Ṅ₀ = [artanh μ_R − artanh μ_S]/[artanh μ_R − artanh μ_L]. Eq. 22 applied to this F returns D₀(1−μ²) exactly. With μ_{L,R} = ∓(1 − ϖ_μ/2), ϖ_μ = 0.05 and μ_S = 0.5, the fractions are Ṅ_L/Ṅ₀ = 0.374284769407670 and Ṅ_R/Ṅ₀ = 0.625715230592330. These are the values of fixture F-SIM-03 (model SIM-M5-ISO).

### 11.4 Time-dependent QLT for the running estimator (Tautz et al. 2013a)

Source: Tautz, Dosch, Effenberger, Fichtner & Kopp (2013a) [P104]. Eqs. 1–8 agree in ar5iv and the A&A HTML. Table 1 was read from the A&A PDF only. Appendix A is garbled; only (A.3) is legible.

With τ = Ωt, R = γv/(Ωℓ₀) and x = ℓ₀k, Eqs. 7–8 and A.3 are:

$$
\langle(\Delta\mu(\tau))^2\rangle=4\pi(1-\mu^2)\left(\frac{\delta B}{B_0}\right)^2\int_0^\infty dk\,G(k)\left[\mathcal K_++\mathcal K_-\right],\qquad
\mathcal K_\pm=\frac{1-\cos[(xR\mu\pm1)\tau]}{(xR\mu\pm1)^2}
\tag{107}
$$

$$
D^{\rm QLT}_{\mu\mu}(\mu)=\frac{2\pi^2}{|\mu|R}(1-\mu^2)\left(\frac{\delta B}{B_0}\right)^2G\!\left(k=\frac{1}{|\mu|R}\right)
\tag{108}
$$

Both are in units of Ω. The normalization of G in (A.2) was reconstructed by the text extraction and is not used here.

Derived and checked in this review:
- **Long-time limit.** With k in units of 1/ℓ₀, ⟨(Δμ)²⟩/(2τ) from Eq. (107) tends to Eq. (108) as τ → ∞, independently of the normalization of G. The ratio is 1.00089 at τ = 10⁴ for μ = 0.5, R = 0.1.
- **Fixture at τ = 10².** For G ∝ (1+x²)^{−5/6}, the ratio at τ = 10² is 1.08841010334818 (μ = 0.5, R = 0.1), 1.00754873626350 (μ = 0.5, R = 1) and 1.15872934147989 (μ = 0.2, R = 0.1). τ = 10² is the slab evaluation time of the paper. A pusher's M1a at the same τ should therefore be compared with Eq. (107), not with Eq. (108).

### 11.5 Published numerical settings

| Source | Turbulence and level | Particles, pusher, step | Estimator window, bins |
|---|---|---|---|
| van den Berg et al. 2024, Case 1 (re-run of Tautz et al. 2013a) | slab, δB²/B₀² = 10⁻³, l_b = 0.03 au, B₀ = 4 nT, R_L/l_b = 0.1 (≈ 144 MeV) | ≥ 8 × 10⁴ per run; Vay (2008) pusher, δt = P_gyro/64 | total 6 × 10⁴ λ_c/v, ϖ_t = 0.04 λ_c/v; M1a after ≈ 16 gyrations; 41 μ bins (ϖ_μ ≈ 0.05) |
| Case 2 (re-run of Qin & Shalchi 2009) | 20 % slab + 80 % 2D, δB²/B₀² = 9 × 10⁻⁴, l_b^2D/l_b^slab = 0.1, R_L/l_b = 0.05 (≈ 38 MeV) | as above | 10⁵ λ_c/v, ϖ_t = 64δt |
| Case 3 (re-run of Sun et al. 2016) | slab, δB²/B₀² = 1, l_b = 0.01 au, B₀ = 5 nT, R_L/l_b = 0.756 (1 GeV); gyrofrequency from the total field | 1.6 × 10⁵ | 10³ λ_c/v, ϖ_t = 6.25 × 10⁻² λ_c/v; M2b fitted only up to 10 gyrations; M5 \|μ_S\| = 0.866 with 2¹⁵ steps per gyration |
| Common settings of van den Berg et al. 2024 | spectra ∝ [1 + (l_bk)²]^{−5/6}, λ_c ≈ 0.747 l_b; slab N_z = 2²³, L_z = 10 au; 2D N_x = N_y = 2¹⁴ | 20 runs (10–25 realisations) | M5 boundaries at ±(1 − ϖ_μ/2), dF/dμ by linear regression, weak slab μ_S = ±0.5 (assignment possibly reversed in the extraction); M2b fit T_μ(1 − e^{−t/T_μ}); M1b 1.5× cumulative-standard-deviation criterion |
| Tautz et al. 2013a | slab and isotropic; δB/B₀ = 10⁻², 10⁻¹·⁵, 10⁻¹; R = 10⁻², 10⁻¹, 1, 10 | number not stated; N_m = 512 modes, k_minℓ₀ = 10⁻⁴, k_maxℓ₀ = 10⁴; vt_max/ℓ₀ = 1–10² (the bracketed values of the Section 11.6 table) | slab D_μμ at τ = 10²; diffusive after t ≃ 5Ω⁻¹, admissible for t ≥ 50Ω⁻¹; D_μμ decays ≈ τ^{−0.14} after τ = 10² |
| Ivascenko et al. 2016 [P108] | GISMO incompressible MHD, (δB/B₀)² ≈ 10⁻³, B₀ = 0.174 G, v_A = 1.2 × 10⁸ cm/s, 256³ | ≥ 10⁶ protons, v = 1.21 × 10¹⁰ cm/s; step not given | MI shown over 1–30 gyroperiods, figures after 10; initial f = a(μ ± 1)² + c; frame of μ not stated in the readings |
| Lange et al. 2013 [P109] | driven MHD with amplified modes, δB/B₀ ~ 10⁻²–10⁻³ | 10⁵ protons, Boris | "between 10 and 30 gyration periods" |
| Xu & Yan 2013 | compressible MHD 512³, M_A = 0.19–0.73 and 1.5 | 1000 per run, Bulirsch–Stoer | "short time interval" (not quantified); μ relative to the local field |
| Riordan & Pe'er 2019 [P45] | slab + 2D, s = 5/3, ℓ = 1, r_g = 2π (code units) | v = 0.1c | Δt = 20t_g (t_g ≪ Δt ≪ 1/⟨D_μμ⟩); t_max ≥ 2t_s |
| Dundovic et al. 2020 | slab, isotropic, 2D; slab δB/B₀ = 0.01 | Boris (no difference from Cash–Karp RK) | lags 10⁻² and 10⁻¹ λ∥/c, window T = t_max/5, t_max ≃ 10λ∥/c |
| Weidl et al. 2015 [P110] | imbalanced MHD 512³, σ_c = 0, 0.9 | — | not reached in the reading |

Printed inconsistencies in these settings:
- Dundovic et al. (2020): the slab grid is given as N_z = 2²⁵, L_z = 400 kpc, l_slab = 10⁻⁴L_z in one place and as 2³⁰, 384 kpc, 10⁻⁶L_z in another; the 2D box ratio is given as 100 and as 10³.
- Lange et al. (2013): L_scale ≈ 10⁸ cm in Eq. 29 and 3.4 × 10⁸ cm in the text.

Further inconsistencies are listed in Section 19.1.

### 11.6 Published simulated D_μμ usable as benchmarks

- **SIM-TAUTZ13-T1.** Tautz et al. 2013a, Table 1 (A&A PDF): χ² between simulated and theoretical slab D_μμ. Entries are QLT/SOQLT χ², with vt_max/ℓ₀ in brackets.

  | δB/B₀ | R = 10⁻² | R = 10⁻¹ | R = 1 | R = 10 |
  |---|---|---|---|---|
  | 10⁻² | 24.47 / 24.33 (1) | 6.55 / 6.43 (10) | 2.53 / 28.81 (100) | — |
  | 10⁻¹·⁵ | 16.74 / 16.60 (1) | 10.23 / 10.37 (10) | 13.55 / 27.71 (100) | 2.24 / no SOQLT entry (100) |

  Q = 0.9057–1. The A&A HTML text gives SOQLT χ² = 14.05 for (10⁻¹·⁵, R = 1), whereas the PDF table gives 27.71; the conflict is not resolved. The conclusion as printed is "excellent agreement with QLT" for slab turbulence.
- **SIM-TAUTZ13II-F3.** In Tautz 2013b, Fig. 3, isotropic scattering at R = 0.1 with a linear fit gives D = 1.2324 × 10⁻³ Ω. Fig. 4 uses the slab model D_μμ = (1−μ²)|μ|^{s−1}D (Eq. 27) with s = 5/3 and μ₀ = 0.5.
- **SIM-VDB24-C1…C3.** The three configurations of Section 11.5. The results are in figures only. In Case 2, M0 and M3 give κ∥^{Dμμ}/κ∥^{QLT} = 0.11. None of the nine methods reproduced the Qin & Shalchi (2009) coefficients [P105]; the paper states it is "unclear how Qin & Shalchi constructed their PADCs".
- **Qin & Shalchi 2009 and 2014** [P105, P114]. Full texts were not accessible. Secondary information: composite 20/80 slab/2D turbulence in 2009; 9 μ bins (ϖ_μ ≈ 0.2) in 2014.
- **Test-particle λ∥ tables** (Hussein et al. 2015 [P115]) are given in MFP v1.2, Section 7.6.
- **Analytic references for a pusher.** SIM-TAUTZ13II-ISO, SIM-M2B-EXACT, SIM-M5-ISO and SIM-TAUTZ13-TDQLT (Sections 11.2–11.4).

### 11.7 Where the estimators fail

- **Strong turbulence, δB/B₀ ≈ 1.**
  - Ivascenko et al. (2016) state that MI "fails at big wave amplitudes".
  - van den Berg et al. (2024) find that M1a and M1b break down because the mean-square displacement saturates too fast. M2a and M3 oscillate; M4a and M4b are very noisy (about two gyrations available).
  - M2b and M5 agree with each other in Case 3.
  - Engelbrecht et al. (2022) state that at strong turbulence pitch angles fluctuate too much to define D_μμ [P74].
- **Plateau of the MSD estimator.** Riordan & Pe'er (2019) find D_μμ ∝ δB² at weak turbulence, flattening near δB/B₀ ≈ 10⁻¹. They call the plateau an artefact of measuring a bounded quantity over too long a time. Their diagnostic is ⟨(Δμ)²⟩ ∝ Δt^b with b = d ln D_μμ/d ln Δt + 1. As read, D_μμ ∝ Δt², Δt⁰ and Δt⁻¹ in the ballistic, diffusive and tail regimes (Eqs. 13–14). With b = d ln D_μμ/d ln Δt + 1 and ballistic ⟨(Δμ)²⟩ ∝ Δt², the ballistic D_μμ ∝ Δt (derived), so the Δt² of the reading refers to ⟨(Δμ)²⟩ or is a misprint (D-124).
- **Late times.** Tautz et al. (2013a) state that all Fokker–Planck coefficients tend to zero at large t, dominated by a 1/t term. Eq. (103) gives the bounded form (⅓ + μ₀²)/(2t).
- **Validity of QLT and SOQLT.**
  - Dundovic et al. (2020) find QLT agreement for δB/B₀ ≲ 0.05 and SOQLT better for 0.05 ≲ δB/B₀ ≲ 0.5 (Sect. V). Slab δB/B₀ = 0.05 needs ≳ 10⁴Ω⁻¹ to reach the plateau.
  - Reichherzer et al. (2020) find QLT valid only for b/B ≲ 0.05 [P116]. That result concerns κ∥, not D_μμ.
- **Focusing fields.** In a focusing field the TGK (M3) estimator includes a deterministic focusing term and can become negative (Section 12.6) [P112].
- **Dynamical MHD fields.** The readings do not record the frame of μ for Ivascenko et al. (2016) and Lange et al. (2013); Xu & Yan (2013) measure μ relative to the local field.

**Implementation requirement.** A pusher validated against a prescribed D_μμ (Eq. (25)) is tested with an estimator consistent with Eq. (22). The estimators that meet this condition are:
- The conditional increment ⟨(Δμ)²|μ₀⟩/(2Δt), over a lag Δt that is small compared with 1/D_μμ. For isotropic D_μμ its expected value is Eq. (103)/(2Δt), with the first-order bias of Eq. (105). The bins next to μ = ±1 need correspondingly smaller lags.
- M5 or M4a/MIIb, which invert Eq. (21) exactly (SIM-M5-ISO).

M2b is used only with the inversion of Eq. (106), not with the printed Eq. 14, unless D_μμ ∝ (1−μ²). Every reported estimate states the estimator, the lag or time window, the μ-bin width and the frame of μ. An estimator printed without the factor ½ (Riordan & Pe'er 2019, Eq. 11) is converted explicitly and never mixed silently with Eq. (22).

## 12. Adiabatic focusing and the D_μμ–λ∥ relation

Eq. (23) is the diffusion limit of the focused transport equation without focusing. With a focusing length L the relation between D_μμ and the parallel diffusion coefficient changes. This section collects the published forms. They are different published results, obtained for different definitions and settings where the sources say so. The section reproduces and compares them; it does not decide between them. The companion documents treat this point:
- PARALLEL rev. 1.4, Section 2.3, states that a focusing correction "cannot be inserted as a universal multiplier".
- MFP v1.2, Section 6.9, supplies He & Wan (2013) only as the explicit option MFP-FOCUS-HW13.

### 12.1 Setting and notation

All sources in this section use the one-dimensional focused transport equation Eq. (1) with L = −B/(∂B/∂z). L > 0 where the field weakens along z, as stated by Wang & Qin (2018) and by Effenberger & Litvinenko (2014). Wang & Qin (2018, Eq. 1) write the focusing term inside the μ-derivative: ∂_μ[D_μμ∂_μf − (v/2L)(1−μ²)f].

Notation of this review:
- λ₀ is λ∥ from Eq. (23); κ₀ = vλ₀/3; ξ = λ₀/L.
- For isotropic scattering D_μμ = D₀(1−μ²), λ₀ = v/(2D₀) [P117, Eq. 4].
- For this D_μμ, Wang & Qin's ξ = v/(2DL) is the same quantity.

### 12.2 Itô equations and weak focusing (Litvinenko 2012a)

Source: Litvinenko (2012a), Full text ×1 [P117]. Eqs. 8–9:

$$
dz=\mu v\,dt,\qquad d\mu=\left[\frac{v}{2L}(1-\mu^2)+\frac{\partial D_{\mu\mu}}{\partial\mu}\right]dt+\sqrt{2D_{\mu\mu}}\,dW(t)
\tag{109}
$$

Here ⟨dW⟩ = 0 and ⟨dW²⟩ = dt. This is Eq. (25) with the focusing drift.

For isotropic scattering, Eqs. 11–13 give:
- dz = [v²/(4D₀L)](1−μ²)dt + v√((1−μ²)/(2D₀)) dW;
- κ∥ = [v²/(4D₀)]⟨1−μ²⟩.

The moment equations give ⟨μ²⟩ = ⅓[1 + (1/30)v²/(D₀²L²) + ⋯] (Eq. 17) and κ∥ = [v²/(6D₀)][1 − (1/60)v²/(D₀²L²) + ⋯] (Eq. 18). In terms of ξ:

$$
\frac{\lambda}{\lambda_0}=1-\frac{1}{15}\left(\frac{\lambda_0}{L}\right)^2+\cdots
\tag{110}
$$

Eq. (110) is Eq. 5. Litvinenko attributes it to Earl (1981), Beeck & Wibberenz (1986), Bieber & Burger (1990) and Kóta (2000). The stated validity is λ₀/L ≪ 1 and ⟨μ⟩ ≈ v/(6D₀L) ≪ 1.

Litvinenko (2012a, Eq. 6) quotes Shalchi (2011) [P118] as

$$
\frac{\lambda}{\lambda_0}=1-\frac{2}{5}\left(\frac{\lambda_0}{L}\right)^2+\cdots
\tag{111}
$$

He states that Shalchi (2011) is in error because Shalchi "erroneously postulated the complete loss of memory of the initial angular distribution". Shalchi (2011) itself was not accessible, so Eq. (111) is secondary.

Derived and checked in this review: with v²/(D₀²L²) = 4ξ², Eqs. 17 and 18 are the second-order terms of the stationary distribution f ∝ e^{ξμ} of Section 12.3, and Eq. 18 is Eq. (110).

### 12.3 Closed form for isotropic scattering (Litvinenko 2012b; Effenberger & Litvinenko 2014)

Sources:
- Litvinenko (2012b), Full text ×1 [P119]. Its Eq. 19 is garbled, and Eqs. 17–18 and 21 were "recovered from garbled text".
- Effenberger & Litvinenko (2014), Full text ×2 for Eqs. 8, 12–13 and 28–30 [P2].

The isotropic results are Eq. 20 of Litvinenko (2012b), κ/(vL) = coth ξ − 1/ξ, and Eq. 8 of Effenberger & Litvinenko (2014):

$$
\kappa_\parallel=\lambda_0v\left(\frac{\coth\xi}{\xi}-\frac{1}{\xi^2}\right),\qquad \xi=\frac{\lambda_0}{L},\quad\lambda_0=\frac{v}{2D_0}
\tag{112}
$$

Here κ∥ is the coefficient of the diffusion approximation for the isotropic phase-space density F₀ = ½∫f₀dμ: ∂F₀/∂t − u∂F₀/∂z = κ∥∂²F₀/∂z², with u = κ∥/L (source Eq. 5). The isotropic linear density is F = e^{z/L}F₀ (source Eq. 6); its equivalent diffusion-approximation equation is ∂F/∂t + u∂F/∂z = κ∥∂²F/∂z². The opposite displayed advection signs arise from the density transformation, not from opposite physical particle streams. The limits are κ∥ → κ∥,0 = v²/(6D₀) as ξ → 0 and κ = vL as ξ → ∞ [P119, Eq. 22].

The dimensionless telegraph coefficients are κ∥ = coth ξ/ξ − 1/ξ² and τ = tanh ξ/ξ (Eqs. 12–13); w = √(κ∥/τ) → 1/√3 as ξ → 0 (Effenberger & Litvinenko 2014, Eq. 21; one reading). The anisotropic generalization of Litvinenko (2012b) is:
- D_μμ = D₀(1−μ²)|μ|^{q−1}[1 + σ sgn(eB₀μ)] (Eq. 4);
- λ₀ = [3/(2(2−q)(4−q))]v/D₀ (Eq. 24);
- κ = [1/(1−σ²)]v²/(6D₀) (Litvinenko 2012b, Eq. 21, recovered from a garbled rendering).

Derived in this review: Eq. (23) with the σ-term gives the extra factor 1/(1−σ²). At q = 1 this matches the recovered Eq. 21.

The Milstein step of Effenberger & Litvinenko (2014, Eq. 30) is printed as follows in both renderings, with time in units of 1/(2D₀) and ε_t ~ N(0, 1):

$$
\mu_{t+\Delta t}=\left[\tfrac12\xi(1-\mu_t^2)-\mu_t\right]\Delta t+\sqrt{\Delta t\,(1-\mu_t^2)}\,\varepsilon_t-\tfrac12\mu_t\Delta t\left(\varepsilon_t^2-1\right)
\tag{113}
$$

As printed, the leading "μ_t +" term is missing. The paper uses reflecting boundaries at μ = ±1, 10⁷ particles and bin width 0.1, with ξ = 0 and 1.5 and λ₀ = 0.2 AU in the appendix example. Δt is not given.

Derived and checked in this review:
- **Nondimensional drift.** Eq. (109) with D₀(1−μ²), in units t′ = 2D₀t and z′ = z/λ₀, gives the drift ½ξ(1−μ²) − μ, the noise √(1−μ²) and the Milstein term −½μ(ε²−1)Δt. These are the increments of Eq. (113).
- **Monte Carlo check.** With μ_t restored (2 × 10⁵ particles, Δt = 10⁻³, reflection at ±1), the stationary ⟨μ⟩ agrees with coth ξ − 1/ξ to within 1 % at ξ = 0.5 and 1.5.
- **Stationary distribution.** The zero-flux stationary density of the Fokker–Planck equation for the linear (flux-tube) density f = e^{z/L}f₀ (Effenberger & Litvinenko Eq. 3), equivalently of the process Eq. (109) with reflection at μ = ±1, is uniform in z and proportional to e^{ξμ} in μ (zero angular flux: (v/2L)(1−μ²)f − D_μμ∂_μf = 0). The corresponding stationary solution of Eq. (1) for the phase-space density is f₀ ∝ e^{−z/L}e^{ξμ}; a z-independent e^{ξμ} is not a solution of Eq. (1), in which it leaves the residual (vμ/L)e^{ξμ} = 2D₀ξμ e^{ξμ} (derived, F-FOC-05; version 1.0 of this document called e^{ξμ} the stationary distribution of Eq. (1), Section 25.6). With the stationary angular density ρ₊(μ) = ξe^{ξμ}/(2 sinh ξ):
  - ⟨μ⟩ = coth ξ − 1/ξ, so Litvinenko's κ = vL⟨μ⟩ (2012b, Eq. 17) is Eq. (112).
  - (v²/4D₀)⟨1−μ²⟩ (2012a, Eq. 13) gives the same function, because the stationary mean of Litvinenko's Eq. 11 is ⟨μ⟩ = (ξ/2)⟨1−μ²⟩.
- **Small-ξ series.** Eq. (112) has the series

$$
\frac{\kappa_\parallel}{\kappa_0}=3\left(\frac{\coth\xi}{\xi}-\frac{1}{\xi^2}\right)=\sum_{n\ge1}\frac{3\cdot2^{2n}B_{2n}}{(2n)!}\xi^{2n-2}
=1-\frac{\xi^2}{15}+\frac{2\xi^4}{315}-\frac{\xi^6}{1575}+\frac{2\xi^8}{31185}-\frac{1382\,\xi^{10}}{212837625}+\cdots
\tag{114}
$$

  The series converges for |ξ| < π, and its second-order term is Eq. (110).

**What Eq. (112) is, and what the process Eq. (109) produces (derived).** Eq. (112) is the coefficient of the diffusion approximation: Effenberger & Litvinenko (2014) derive it, together with the telegraph coefficients, as an approximation to the full transport equation and compare both against its numerical solution. Its advection speed u = κ∥/L = v(coth ξ − 1/ξ) is the exact stationary mean velocity v⟨μ⟩ of the process (derived from ρ₊). Its κ∥ is not the long-time growth rate of the centred spatial variance of that process. For the constant-L, constant-v process Eq. (109) with D_μμ = D₀(1−μ²) (ν = 2D₀, ξ = v/(νL)), the long-time limit of Var[z(t)]/(2t) follows from the stationary generator (Poisson equation of the pitch-angle process; Kubo form):

$$
\kappa_{\rm var}\equiv\lim_{t\to\infty}\frac{{\rm Var}[z(t)]}{2t}=\int_{-1}^{1}\frac{H(\mu)^2}{D_{\mu\mu}(\mu)\,\rho_+(\mu)}\,d\mu,\qquad
H(\mu)=\int_{-1}^{\mu}v\,(y-\langle\mu\rangle)\,\rho_+(y)\,dy,\qquad
\rho_+(\mu)=\frac{\xi\,e^{\xi\mu}}{2\sinh\xi}
\tag{115}
$$

(H(±1) = 0; for ξ ≠ 0, H(μ) = v[e^{ξμ}(μ − coth ξ) + e^{−ξ}(1 + coth ξ)]/(2 sinh ξ) in closed form). Eq. (115) is the integral ∫D_μμ(φ′)²ρ₊dμ of the solution φ of ρ₊⁻¹∂_μ(D_μμρ₊∂_μφ) = −v(μ − ⟨μ⟩), i.e. the time integral of the stationary autocovariance of vμ. It contains no closure. With v = 1, ν = 1 (fixture F-FOC-04, quadrature class; checked by a conservative finite-volume solution of the Poisson equation with 1000 and 4000 cells, which gives 0.178275493 and 0.178275762 at ξ = 1.5, and by a Monte Carlo run of Eq. (109), MCR-8):

| ξ | κ_var, Eq. (115) | κ∥, Eq. (112) | κ_var/κ∥ | ⟨μ⟩ = coth ξ − 1/ξ |
|---|---|---|---|---|
| 0 | 0.333333333333333 | 0.333333333333333 | 1 | 0 |
| 0.1 | 0.332298461861175 | 0.333111322539896 | 0.997560 | 0.033311132254 |
| 1 | 0.248076740609470 | 0.313035285499331 | 0.792488 | 0.313035285499 |
| 1.5 | 0.178275780337345 | 0.292083150877230 | 0.610360 | 0.438124726316 |
| 5 | 0.010324829484838 | 0.160018160796404 | 0.064523 | 0.800090803982 |

Both tend to v²/(3ν) as ξ → 0. A pusher that integrates Eq. (109) correctly and measures the centred spatial variance reproduces κ_var, not Eq. (112); at ξ = 1.5 the two differ by 39 %. The raw second moment ⟨z²⟩/(2t) contains the term v²⟨μ⟩²t/2 and has no finite limit when ⟨μ⟩ ≠ 0. Wang & Qin (2019) [P120] state that the identification of the parallel diffusion coefficient with the displacement variance holds only for isotropic transport equations without cross terms or without a convection term, and derive a modified relation from the Fokker–Planck equation (Abstract only; their equations were not read). Consequences:

- Eq. (112) and Eq. (115) are different observables of the same process. A Monte Carlo value is compared only with the reference defined by the same observable and the same initial and boundary problem (Section 16.9; Section 20; PS08 and PS15 of Section 26).
- Backward trajectories. The backward SDE Eq. (152) with dμ/dt = v(1−μ²)/(2L) has the drift −v(1−μ²)/(2L) + D′_μμ; with reflection at μ = ±1 its stationary angular density is ρ₋(μ) = ξe^{−ξμ}/(2 sinh ξ) (zero flux, derived, F-FOC-05). ρ₋ is a property of the sampling process; the physical distribution is reconstructed from the boundary, initial and source values along the trajectories (Section 16.1), so a histogram of backward pitch cosines is not compared with ρ₊.
- Manufactured stationary solution. f₀(z, μ) = e^{−z/L}e^{ξμ} with its exact boundary data is a stationary solution of Eq. (1) (residual 0, derived, F-FOC-05) and serves a forward and a backward reconstruction test with the same reference (PS15).
- Sign of L. Under z → −z, μ → −μ, L → −L, Eq. (109) with an even D_μμ maps onto itself: ⟨μ⟩ changes sign, while κ_var and Eq. (112) are even in ξ (derived from Eq. (115); F-FOC-04 stores ξ = −1.5). The focusing strength is therefore carried as the signed inverse length 1/L (zero for a uniform field), and a negative ξ is admissible for even D_μμ (Section 12.7).

### 12.4 Wang & Qin (2018): κ∥^a and the additional term

Source: Wang & Qin (2018) [P121]. Sect. II was read in ar5iv and Sects. III onward in the arXiv PDF, one rendering each. Eqs. 10 and 39:

$$
M(\mu)=\frac{v}{2L}\int_{-1}^{\mu}\frac{1-\nu^2}{D_{\mu\mu}(\nu)}\,d\nu,\qquad
\kappa_\parallel^a=vL\,\frac{\int_{-1}^{1}\mu\,e^{M(\mu)}d\mu}{\int_{-1}^{1}e^{M(\mu)}d\mu}
\tag{116}
$$

They write κ∥^a = κ∥0 + T (Eq. 41) and state that this agrees with Beeck & Wibberenz (1986), Litvinenko (2012a) and He & Schlickeiser (2014) [P122]. T is attributed to Shalchi & Danos (2013) [P123]. Their full coefficient is κ∥ = κ∥^a + η_{0,2,0} = κ∥0 + T′ (Eqs. 54–56).

For isotropic scattering:
- M = ξ(μ + 1) and κ∥0 = v²/(6D);
- T = κ∥0S with S = −ξ²/15 + 2ξ⁴/315 + ⋯ (Eqs. B2–B3);
- η_{0,2,0} ≈ (v²/(30D))ξ² (Eqs. B16–B17).

Eqs. 62–63 then give

$$
\kappa_\parallel\approx\kappa_{\parallel0}\left(1+\frac{2}{15}\xi^2\right)
\tag{117}
$$

Its ξ² term has the opposite sign to that of Eq. (110). It refers to the full κ∥, which includes η_{0,2,0}, whereas κ∥^a alone has the −ξ²/15 term. The text extraction shows apparent printed typos: "ρ^{M(μ)}" for e^{M(μ)} in Eqs. 23/(25), "κ_κ0", and "νL" in Eq. A5.

Derived and checked in this review:
- **Isotropic case of Eq. 39.** Eq. (116) with M = ξ(μ + 1) equals Eq. (112), i.e. κ∥^a = vL⟨μ⟩ with f ∝ e^{M}. The quadrature agrees to 25 digits at ξ = 0.3 and 1.5.
- **Printed series.** S is Eq. (114) minus 1.
- **Sum of the corrections.** 1 + S + (6/30)ξ² gives the coefficient 2/15 of Eq. (117).

The same stationary-anisotropy form appears in two other sources:
- λ∥ = 3L∫μe^G dμ/∫e^G dμ in van den Berg, Strauss & Effenberger (2020, Eq. 14) [P1]. G is defined in an appendix that was not read.
- f̃_iso = (ξ csch ξ/2)e^{ξμ} with first-order anisotropy A₁ = 3(coth ξ − 1/ξ) in van den Berg, Engelbrecht & Strauss (2026, Eqs. 47–48) [P14], where ξ = λ∥/L′ and G = ξ(1 + μ); one rendering.

### 12.5 Tautz, Dosch & Lerche (2012) and He & Wan (2013)

Source: Tautz, Dosch & Lerche (2012), aanda.org Full text ×1 [P124]. The paper prints a "quasilinear" focusing result (Eq. 9) and a "nonlinear" result (Eq. 10):

$$
\lambda'_\parallel=\lambda_{\parallel,0}\left[1+\frac{4-s}{3}\,\frac{\lambda_{\parallel,0}}{L}\right]^{-1},\qquad
\lambda''_\parallel=3\lambda_{\parallel,0}\left(\frac{L}{\lambda_{\parallel,0}}\right)^3\left(\frac{\lambda_{\parallel,0}}{L}-\tanh\frac{\lambda_{\parallel,0}}{L}\right)
\tag{118}
$$

λ∥,0 is given by Eq. (23), printed as Eq. 6 of the source. The symbol s in Eq. 9 is not defined in the reading available here.

The simulations use the running κ∥(t) = ⟨(Δr∥)²⟩/(2t) with λ∥ = 3κ∥/v (Eq. 14). That equation was re-arranged by the text extraction, and whether the mean displacement is subtracted was not read. The paper finds that λ∥ increases in almost all cases, contradicting both Eqs. 9 and 10. It also finds that diffusion is asymmetric, depending on the initial direction, and attributes this to coherent motion toward the weaker field.

He & Wan (2013, Eqs. 2–3) [P125] use λ∥ = (3L³/λ∥,0²)[λ∥,0/L − tanh(λ∥,0/L)] with a Parker-field L. This is MFP v1.2, Section 6.9, Eq. 18, with the stable series in its Section 21.4, Eq. 43.

Derived and checked in this review:
- **Identity.** The He & Wan form is identical to λ″∥ of Eq. (118).
- **Series.** It is 3(ξ − tanh ξ)/ξ³ = 1 − 2ξ²/5 + 17ξ⁴/105 − ⋯, whose ξ² coefficient equals that of Eq. (111). The sources read here do not state the origin of Tautz et al. (2012) Eq. 10. Litvinenko (2012a) states that the −2/5 coefficient is in error.
- **Large ξ.** For ξ → ∞, λ″∥/λ∥,0 → 3/ξ², whereas Eq. (112) gives κ → vL.

### 12.6 Focusing contribution to D_μμ (Tautz, Shalchi & Dosch 2014)

Source: Tautz, Shalchi & Dosch (2014) [P112]. Eq. 13 appears in ar5iv and the PDF; Eqs. 14–34 were read from the PDF only.

The field is B_z = B₀e^{−z/L} with B_{x,y} = B₀(x, y)/(2L)e^{−z/L} (Eqs. 3a–c), so L⁻¹ = −(1/B)∂B/∂s (Eq. 2). The unperturbed focused orbit is (Eqs. 7a, 8a, 9a):

$$
v_\parallel=v\tanh\!\left(\frac{vt}{2L}+{\rm artanh}\,\mu_0\right),\quad
v_\perp=v\,{\rm sech}\!\left(\frac{vt}{2L}+{\rm artanh}\,\mu_0\right),\quad
z=2L\ln\cosh\!\left(\frac{vt}{2L}+{\rm artanh}\,\mu_0\right)+L\ln(1-\mu_0^2)+z_0
\tag{119}
$$

Using D_μμ = ∫₀^∞⟨μ̇(t)μ̇(0)⟩dt = ½d⟨(Δμ)²⟩/dt (Eqs. 11a–b) and ∫₀^∞(1 − μ(t)²)dt = (2L/v)(1 − μ₀) (Eq. 12), the paper obtains the large-scale contribution (Eq. 13) and its mean free path (Eq. 14):

$$
D^{(0)}_{\mu\mu}(\mu_0)=\frac{v}{2L}(1-\mu_0^2)(1-\mu_0),\qquad \lambda_\parallel^{(0)}=\frac{3v}{8}\int_{-1}^{1}\frac{(1-\mu_0^2)^2}{D^{(0)}_{\mu\mu}}d\mu_0=\frac{3L}{2}
\tag{120}
$$

Further printed results:
- The turbulent part, Eq. 26, is a double integral over the focused orbit, with the spectrum of Eqs. 23–24. Its L → ∞ limit, Eq. 29a, has units that look inconsistent as extracted; the variables may be normalized, so it is to be checked against the original.
- The estimator used in the simulations is ⟨μ̇(0)Δμ(t)⟩ (Eq. 34).
- Parameters: R = 1 with L/ℓ₀ = 0.02–1000 (Figs. 3–5); R = 1, L = ℓ₀ and 2 × 10⁷ particles (Fig. 6).
- D_μμ becomes negative for μ₀ ∈ (−1, 0) when L ≪ R_L, and only for Ωt ≳ 2.
- With "the modified relation", λ∥^(0) = 3L/2 is reduced by a factor of about 0.8 (as extracted; the relation is not reproduced here).

Derived and checked in this review:
- **Eqs. 12–14.** They follow from Eq. (119): μ̇(t) = (v/2L)sech²(·), so ∫₀^∞μ̇(t)μ̇(0)dt = μ̇(0)[1 − μ₀].
- **Meaning of D^(0).** D^(0) is the correlation integral of the deterministic focusing drift of Eq. (1). It is not a diffusion term of Eq. (21).

**Implementation requirement.** A D_μμ estimated from trajectories in a focusing field states whether the deterministic focusing contribution is included. Tautz, Shalchi & Dosch (2014) show that the correlation estimator includes it (Eq. (120) for the unperturbed orbit).

### 12.7 Numerical comparison (derived)

Values of κ∥/κ₀ or λ∥/λ₀ for isotropic scattering, evaluated with mpmath at 40 digits:

| ξ | Eq. (112) (E&L 2014; Litvinenko 2012b; Wang & Qin κ^a) | Eq. (118) λ″ (= He & Wan 2013) | Eq. (110) | Eq. (111) | Eq. (117) | Eq. (118) λ′, s = 5/3 (illustrative) |
|---|---|---|---|---|---|---|
| 0.1 | 0.999333967619688 | 0.996016125132549 | 0.999333 | 0.996 | 1.001333 | 0.927835051546392 |
| 0.3 | 0.994050969884083 | 0.965265283156566 | 0.994 | 0.964 | 1.012 | 0.810810810810811 |
| 0.5 | 0.983720482431917 | 0.909188225759766 | 0.983333 | 0.9 | 1.033333 | 0.72 |
| 1 | 0.939105856497994 | 0.715217532132705 | 0.933333 | 0.6 | 1.133333 | 0.5625 |
| 1.5 | 0.876249452631690 | 0.528757107871230 | 0.85 | 0.1 | 1.3 | 0.461538461538462 |
| 2 | 0.805972081091322 | 0.388489657471569 | 0.733333 | −0.6 | 1.533333 | 0.391304347826087 |
| 3 | 0.671636489980356 | 0.222771694034808 | 0.4 | −2.6 | 2.2 | 0.3 |

How to read the table:
- **Truncated forms.** Eq. (110), Eq. (111) and Eq. (117) are printed as truncated series for ξ ≪ 1. Their values at ξ ≳ 1 show only where the truncation departs from the closed forms.
- **Differences in definition.** The columns differ in what is defined:
  - Eq. (112) is the κ∥ of the diffusion–advection equation for F₀ with u = κ∥/L.
  - Eq. (117) includes Wang & Qin's additional term η_{0,2,0}.
  - Eq. (118) is stated by Tautz et al. (2012) as quasilinear (λ′) and nonlinear (λ″) results. Their simulations contradict both.
  - Shalchi (2011) is known only through Litvinenko (2012a), who states it is wrong.
- **Pure-focusing scale.** λ∥^(0)/λ₀ = 3/(2ξ) from Eq. (120) refers to focusing without scattering and is not a correction factor.

**Implementation requirement.** Each focusing relation is a separately named option with its source and definition: FOC-EL14, FOC-WQ18-A, FOC-WQ18-SERIES, FOC-LITV12A-WEAK, FOC-SH11-QUOTED, FOC-T12-QL, and FOC-HW13 (= MFP-FOCUS-HW13 = FOC-T12-NL). None is applied as a default multiplier when D_μμ is converted to λ∥ (MFP v1.2 Section 6.9; PARALLEL rev. 1.4 Section 2.3).

**Implementation requirement.** Eq. (112) is evaluated in double precision from Eq. (114) through ξ⁸ for |ξ| ≤ 0.1 and from the closed form above that. This choice was derived and checked: with the closed form written as 3(coth ξ/ξ − 1/ξ²) for κ∥/κ∥,0, the largest relative error on 10⁵ logarithmically spaced points of 10⁻⁸ ≤ ξ ≤ 10³ was below 10⁻¹³ against 40-digit values; with the switch at 0.05 or 0.2 it was below 5 × 10⁻¹³ and 7 × 10⁻¹³ (the figures depend on the sampling and the implementation). The function is even in ξ (derived). The sources read do not discuss L < 0. For an even D_μμ the reflection z → −z, μ → −μ, L → −L maps the process onto itself (Section 12.3), so a negative ξ is evaluated as |ξ| and the sign of ⟨μ⟩ and of u = κ∥/L follows the sign of L; for a D_μμ that is not even in μ (directional or helical wave models), a negative ξ returns a domain diagnostic unless the configuration admits it (U-44).

## 13. D_μμ computed from measured or prescribed spectra

Several SEP models compute D_μμ from a magnetic spectrum that is measured, fitted, or prescribed by a turbulence or wave model. Most of them use the resonant slab coefficient of Section 6, with printed forms that differ through the spectrum convention (mapped onto Eq. (43) in Section 6.7); the exceptions are the gyrofrequency without γ of Chen et al. (Section 13.2), the printed vμ of Ding et al. (13.6), the wave-frame coefficient of Afanasiev et al. (13.7), the second-order coefficient of Subashchandar et al. (13.8) and the illegible prefactor of Zhang et al. (13.5). This section records each printed form with its definitions and states the information an adapter needs. The companions cover related ground: MFP v1.2 Section 6.8 (printed forms), PARALLEL rev. 1.4 Sections 8.4 and 8.7 (supplied spectra and conversion rules) and Section 12 (adapters).

### 13.1 Frequency spectra and the Taylor hypothesis

A spacecraft measures a frequency power spectral density P_f(f). Under frozen-in flow with sampling velocity U_T along the relevant direction, k = 2πf/U_T and P_f df = P_k dk (PARALLEL rev. 1.4, Eq. 36). Then (derived; DR-06 C11)

$$
f\,P_f(f)=k\,P_k(k),\qquad k=\frac{2\pi f}{U_T},
\tag{121}
$$

so a printed f_res P(f_res) equals k_res P_k(k_res), and the sidedness and component content of P_f carry over to P_k unchanged. The sampling velocity differs between sources. Chen et al. (2024) use f_res = k_res V_SW/2π with V_SW the solar-wind speed (Section 13.2). Subashchandar et al. (2025) use k∥ = 2πf/(V_sw cos ψ) for M_A > 5 and 2πf/(V_sw cos ψ + V_A) for M_A ≤ 5 (Section 13.8). Chen et al. state that the hypothesis requires V_SW ≫ V_A and may fail near the Sun (paraphrase).

### 13.2 Chen et al. (2024): PSP spectra and a μ cut-off

**Source.** Chen, Giacalone, Guo & Klein (2024) [P47]; Full text ×2 (ar5iv and arXiv PDF) for Eqs. 2, 3 and the inline definitions.

**As printed.**

$$
D_{\mu\mu}=\frac{\pi}{4}\,\Omega_0\,(1-\mu^2)\,\frac{f_{\rm res}\,P(f_{\rm res})}{B_0^2},
\qquad
k_{\rm res}=\left|\frac{\Omega_0}{v\mu}\right|,\quad
f_{\rm res}=\frac{k_{\rm res}V_{SW}}{2\pi},\quad
\Omega_0=\frac{qB_0}{mc}
\tag{122}
$$

(Eq. 3 and the unnumbered definitions), and κ∥ = (v²/4)∫_{μ_min}^{1}(1−μ²)²dμ/D_μμ (Eq. 2). v is the particle velocity; μ "the cosine of the pitch angle" (frame not stated); B₀ the magnitude of the mean field; P(f_res) "the value of spectral density at the resonant frequency". P is computed from the RTN normal component B_N only, chosen because it is "always perpendicular to the mean magnetic field (Parker spiral)". **Not stated:** whether P is one- or two-sided, its units, the PSD estimator and the data cadence. Ω₀ is printed without a Lorentz factor. The cut-off is μ*_min = 0.05 ("dashed line in Figure 1(b)"), because κ∥ "blows up" as μ_min → 0 (90° problem). The frequency band runs from the resonant frequency of 1 GeV protons at μ = 1 to that of 100 keV protons at μ = μ_min, where the spectral indices are −3/2 to −5/3; there is no high-frequency extrapolation. The authors state: "Equation 3 is only applicable to magnetostatic fluctuations and does not include the dynamical effects". The κ∥(r, E) fit, Eq. 4, is the model SEP-CHEN24 of MFP v1.2 (Sections 4.5 and 5.2).

**Derived (DR-06 C10, C11; fixtures F-MEAS-01, F-MEAS-02).**

1. With Eq. (121), Eq. (122) becomes (π/4)Ω₀²(1−μ²)P_k(k_res)/(B₀²v|μ|). It equals Eq. (43) if P_k is the one-sided spectrum of the total transverse variance, P_k = P_s = 8πg. If P is the one-sided PSD of B_N alone and the turbulence is axisymmetric, P = P_s/2 and Eq. (122) is half of Eq. (43), so λ∥ is twice as large.
2. For a power-law PSD of index s, Ω₀ without γ gives a D_μμ larger by γ^{2−s} than Ω₀/γ would. For s = 5/3 this is 1.00004 at 100 keV, 1.0035 at 10 MeV and 1.274 at 1 GeV (protons).
3. For a pure k^{−5/3} spectrum the cut at μ_min = 0.05 retains 0.57035 of the full κ∥ integral, and the integral converges as μ_min → 0. It diverges only if the spectrum at the resonant frequencies of small μ is at least as steep as f^{−2} (Section 6.6, item 2); the text read does not show which case applies to the measured spectra.

**Implementation requirement.** MEAS-CHEN24 needs the user decisions on the sidedness and component content of P (U-23) and on Ω₀ (U-47). μ_min is a parameter of the κ∥ integral only. The source gives no prescription for D_μμ at |μ| < μ_min in a transport code, and none is supplied here.

### 13.3 Wijsen et al. (2019): PARADISE

**Source.** [P25]; Full text ×2 for Eqs. 6–8. The equations are

$$
D_{\mu\mu}=\frac{\pi}{2}\,\frac{1-\mu^2}{v|\mu|}\left(\frac{\Omega}{B}\right)^2P\!\left(k=\frac{\Omega}{|\mu|v}\right),
\tag{123}
$$

(Eq. 6, cited to "Jokipii 1966; Jaekel & Schlickeiser 1992"), with P "the power spectrum of the magnetic turbulence" and k "the parallel wave number". With P = Ck^{−q} the paper writes D_μμ = (Cπ/2)(|Q|/m)^{2−q}B^{−q}v^{q−1}(|μ|^{q−1} + H)(1−μ²) (Eq. 7), where H is "added to describe the scattering through μ = 0" and "caution is needed since dD_μμ/dμ has a pole at μ = 0". p and μ are measured in the frame co-moving with the solar wind (Eq. 1). C is fixed by a radial mean free path (Eqs. 11–12; Section 5). The normalization and sidedness of P are not stated, q is not given a number, and whether m is the rest mass is not stated.

**Derived.** Eq. 6 with P = Ck^{−q} and Ω = |Q|B/m gives Eq. 7 with H = 0 (checked, DR-06 C12). Eq. 6 equals Eq. (43) if P = 4πg, i.e. ∫₋∞^∞P dk = δB² (Section 6.7).

### 13.4 Laitinen et al. (2016)

**Source.** [P8]; Full text ×2 for Eqs. 2, 6, 11, 12. Eq. 2, cited to "Jokipii 1966; Matthaeus et al. 2003":

$$
D_{\mu\mu}=\frac{\pi\Omega^2}{v|\mu|B^2}(1-\mu^2)\,S_\parallel\!\left(-(r_L\mu)^{-1}\right),
\tag{124}
$$

with Ω the Larmor frequency, r_L the Larmor radius, and S∥ the slab spectrum. The spectrum is S(k) = S⊥(k⊥)δ(k∥) + S∥(k∥)δ(k⊥) (Eq. 10), with S∥ = C∥L∥δB∥² for L₀⁻¹ < k∥ < L∥⁻¹ and C∥L∥δB∥²(L∥k∥)^{−5/3} for k∥ ≥ L∥⁻¹ (Eq. 11), normalized by "δB²∥,⊥ = 2∫S∥,⊥(k)dk" with no integration limits printed. The paper uses the shape D_μμ = κ₀(|μ|^{q−1} + H)(1−μ²) (Eq. 6) with κ₀ containing "the µ-independent terms in Eq. 2" and H = 0.1. The spectrum is evolved radially with a WKB factor (Eqs. 13–14). Parameters: L∥ = L⊥ = 0.007 AU; L₀(r) = r; slab:2D = 20%:80%; δB²/B² = 0.04 at 1 AU; amplitudes normalized to λ∥ = 0.3 AU for 10 MeV protons at 1 AU; V_sw = 400 km s⁻¹; v_A0 = 30 km s⁻¹ and B₀ = 5 nT at 1 AU. The value of C∥ was not found in the text read.

Eq. 2 evaluates S∥ at a negative wavenumber for μ > 0, while Eq. 11 defines S∥ for positive k∥; an even extension is implied but not stated. Eq. (124) equals Eq. (43) if the printed normalization integrates over all k (S∥ = 2πg); if it integrates over k > 0, S∥ = 4πg and Eq. 2 is twice Eq. (43) (derived, Section 6.7).

### 13.5 Zhang et al. (2023)

**Source.** [P15], arXiv PDF only. The shape used in the model is the q-form D_μμ = D₀(x)p^{q−2}(1−μ²)(|μ|^{q−1} + h₀) with q = 5/3, h₀ = 0.2 and p the rigidity in GV (Eq. 35; Section 5). Eq. 36 motivates it from QLT. **Eq. 36 is not legible.** One extraction gave the raw string "D^QL_μμ = π / 2Ω2 (1 − μ 2) / B2v|μ| W⊥(k_res) with k_res = Ω / v|μ|", whose fraction grouping cannot be determined. Another reading of the same PDF gave π²Ω²(1−μ²)W⊥(k_res)/(B²v|μ|), and MFP v1.2 (Section 6.8) transcribes that reading. The two candidate readings, π²Ω² and (π/2)Ω², differ by 2π. W⊥ is "the spectral power density of the transverse magnetic field fluctuations as a function of wavenumber" (its argument is printed as α⊥, a notation clash); "forward and backward propagating fluctuations have equal power density"; no normalization is given. μ is the cosine of the pitch angle to the outward magnetic field line. Section 6.7 gives the spectrum relation for each reading. MEAS-ZHANG23 is not fully specified.

### 13.6 Ding et al. (2022): iPATH

**Source.** [P42]; Full text ×2 (arXiv PDF and ar5iv) for Eqs. 8–13.

$$
D_{\mu\mu}=\frac{2\pi^2\Omega^2(1-\mu^2)}{B^2\,v\mu}\,g^{\rm slab}(k_\parallel),
\qquad k_\parallel=\Omega(v|\mu|)^{-1},
\qquad \Omega=\frac{eB}{\gamma m},
\tag{125}
$$

(Eq. 8, **printed with vμ in the denominator, not v|μ|**, in both renderings), with g^slab = [C(ν)/(2π)]l_slab δB²_slab(1 + k²l²_slab)^{−ν} "given by Shalchi 2009" (Eq. 9), C(ν) = Γ(ν)/[2√πΓ(ν − ½)] (Eq. 10), and κ∥ = (v²/8)∫₋₁¹(1−μ²)²/D_μμ dμ (Eq. 11). B is not defined explicitly. Parameters: "the inertial range spectral index s = 2ν = 5/3"; l_slab = 10⁹ m; δB²/B² = 0.15 at 1 au (the slab share is not stated); δB² ∝ r^{γ} with γ = −3.5 (γ is also the Lorentz factor); l ∝ r^{1.0}. As printed, D_μμ < 0 for μ < 0 (derived). With v|μ| the coefficient is Eq. (47) with the 8π normalization (derived, Section 6.5).

### 13.7 Afanasiev et al. (2015): wave-frame coefficient and SOLPACS grid

**Source.** [P78]; Full text ×2 (arXiv PDF and ar5iv) for Eqs. 1–11, except the denominator of Eq. 2; Eq. 14 and the grid formulas PDF only.

**As printed.** k_res = Ω/(υμ), with "υ and μ the particle speed and pitch-angle cosine as measured in the wave rest frame", Ω = γ⁻¹Ω₀ (Eq. 1);

$$
D_{\mu\mu}(\mu)=\frac{\pi}{2}\,\Omega\,\frac{|k_{\rm res}|\,I_{w,{\rm res}}}{B^2}\,(1-\mu^2),
\qquad
\delta B^2=\int_{-\infty}^{+\infty}I_w(k)\,dk
\tag{126}
$$

(Eqs. 2 and 5), where "positive and negative values of k correspond to the waves with opposite circular polarizations", B is "the mean magnetic field" and I_w,res is the wave intensity at k = k_res. The denominator of Eq. 2 is B² in the ar5iv rendering of the TeX source and R² in the PDF text extraction; MFP v1.2 records R² as printed in the arXiv version (MFP v1.2 D-30; this document D-37). The B² reading agrees with the definition of B in the text. The earlier code used the simplified resonance k_res = Ω/υ (Eq. 3) and D_μμ = ½ν(1−μ²), ν = πΩ²I_w,res/(υB²) (Eq. 4). λ = (3υ/8)∫₋₁¹(1−μ²)²/D_μμ dμ (Eq. 11). The text states that QLT "does not work close to μ = 0".

**SOLPACS grid (Sect. 3 of the source; q_j, the interpolant, μ_j, ν(μ), ν₀,j and the gap filler: two renderings, PDF text and ar5iv, re-read for version 1.1; Eq. 14: one rendering).** The spectrum is stored on a grid k_j with local indices q_j = ln(I_w,j+1/I_w,j)/ln(k_j+1/k_j) and, "for k ∈ (k_j, k_{j+1})", I_w(k) = I_w,j(k/k_j)^{−q_j}, both as printed; μ_j = Ω/(υk_j); ν(μ) = ν₀,j|μ|^{q_j−1} with ν₀,j = (πΩ/B²)(υ/Ω)^{q_j−1}I_w,j|k_j|^{q_j} (the source parenthesizes only υ/Ω; the grouping of πΩ/B² is this review's reading); λ = (3υ/4)Σ_j∫_{μ_{j+1}}^{μ_j}(1−μ²)/(ν₀,j|μ|^{q_j−1})dμ (Eq. 14). For |μ| < μ_min = Ω/(υk_max) a linear gap filler ν̃(μ) = bμ + b₁ is used; how b and b₁ are fixed was not found in the text read.

**Printed inconsistency (D-128).** The two printed definitions contradict each other: with q_j = ln(I_w,j+1/I_w,j)/ln(k_j+1/k_j) the interpolant that passes through both nodes is I_w,j(k/k_j)^{+q_j}, whereas the printed I_w,j(k/k_j)^{−q_j} reproduces the second node only if q_j = −ln(I_w,j+1/I_w,j)/ln(k_j+1/k_j) (a positive decay index for a decreasing spectrum). Derived (F-MEAS-03): for k_j = 1, k_{j+1} = 2, I_w,j = 1, I_w,j+1 = 1/4 the printed index is −2 and the printed interpolant gives I_w(k_{j+1}) = 4, a factor 16 above the node value; with the decay-index convention q_j = 2 it gives 1/4. The subsequent ν(μ) ∝ |μ|^{q_j−1}, ν₀,j ∝ |k_j|^{q_j} and Eq. 14 are consistent only with the decay-index convention (D_μμ ∝ |k_res|I_w(k_res) ∝ |μ|^{q−1} for I_w ∝ k^{−q}, Section 6.3). Which of the two printed expressions carries the sign error is not established by the text; an implementation stores the signed logarithmic slope α_j = ln(I_w,j+1/I_w,j)/ln(k_j+1/k_j) and uses I_w(k) = I_w,j(k/k_j)^{α_j}, or q_j = −α_j with the printed exponent −q_j, and never the two printed expressions together. Interpolation is done on a positive |k| grid for each signed wave population separately (ln k is not defined for k < 0). Required checks of a spectrum provider: both nodes of every interval are reproduced, a pure power law is reproduced throughout each interval, a decreasing spectrum stays decreasing, and the integrated variance Eq. (126) is conserved by the interpolation to the stated accuracy (PS13).

**Derived (DR-06 C10, C13).** Eq. 2 with the power law I_w,j(k/k_j)^{−q_j} (decay-index convention, q_j > 0 for a decreasing spectrum) gives D_μμ = ½ν₀,j|μ|^{q_j−1}(1−μ²) with exactly the ν₀,j grouping above, and Eq. 14 follows from Eq. 11 for D_μμ = ν(1−μ²)/2 if the k_j grid covers both signs of μ (with positive k_j only, the sum gives half of (3v/4)∫₋₁¹; derived). This check assumed the power-law form and therefore did not expose D-128 (Section 25.6). With Eq. 5, I_w = 4πg for an even I_w, so Eq. (126) is Eq. (43) as printed (the B² reading). Parameters: n₀ = 3.6 × 10¹² m⁻³, B₀ = 3.4 × 10⁻⁵ T, k_max = 10⁻² m⁻¹ (∼ Ω₀/V_A); initial spectrum ∝ k^{−3/2}, normalized so that λ = 1 R☉ for 100 keV protons.

### 13.8 Subashchandar et al. (2025): spectrum adapter

**Source.** [P101]; Full text ×2 (arXiv HTML and PDF) for the main-text Eqs. 1–3 and appendix Eqs. 32–40; the fit model, Eq. 4, and App. Eq. 22 Full text ×1. The second-order QLT coefficient itself is in Section 10. The spectrum input is a fit of the measured PSD, fP_SUM(f) = 2C_S(2πf/(V_sw cos ψ))^{1−q^S} + 2C_2D(2πf/(V_sw sin ψ))^{1−q^{2D}} (Eq. 4), with P_S(k∥) = 2C_S k∥^{−q^S} and P_2D(k⊥) = 2C_2D k⊥^{−q^{2D}}. The slab part satisfies Eq. (121) with U_T = V_sw cos ψ (derived). The indices are fixed by assumption, q^S = 1.5 and q^{2D} = 5/3 (the text writes "−q = 1.5" and "−5/3"; sign convention inconsistent as printed). App. Eq. 22 prints ⟨δB_x(t)δB_x*(0)⟩ = 4π∫₀^∞dk∥ g(k∥)Γ(k∥, t), which at t = 0 is δB²/2 under the 8π convention (derived). The relation between P_S and g, and the sidedness of P_SUM, were not found in the text read. QLT comparisons use μ_min = 0.05 or 0.01.

### 13.9 Adapter information

| Source | Spectrum variable | Sampling velocity | Component | Normalization as printed | Frame of μ | Ω |
|---|---|---|---|---|---|---|
| Chen 2024 | f (Hz) | \|V_SW\| | B_N only | not stated | not stated | qB₀/mc, no γ printed |
| Subashchandar 2025 | f, fitted | V_sw cos ψ (+V_A if M_A ≤ 5) | P_SUM, slab + 2D | not stated (2C_S) | not stated | B₀q/(γm) (App.) |
| Wijsen 2019 | k∥ | — | not stated | not stated | solar-wind frame | via \|Q\|/m |
| Laitinen 2016 | k∥ | — | slab S∥ | 2∫S dk, limits not printed | not stated (μ relative to the local mean field) | Larmor frequency |
| Zhang 2023 | k | — | transverse W⊥ | not stated | not captured (μ relative to the outward field line) | not legible |
| Ding 2022 | k∥ | — | slab g | 8π (derived) | not stated | eB/(γm) |
| Afanasiev 2015 | k, sign = polarization | — | wave intensity | ∫₋∞^∞I_w dk = δB² | wave rest frame | Ω₀/γ |
| AMPS Type III | k | — | P^slab | ∫P dk = δB², limits not legible | not stated | Ω₀ |

**Implementation requirements.**

1. A spectrum adapter declares the spectral variable (f in Hz or k in rad m⁻¹), the sampling velocity and angle, the component content (one component or trace), the sidedness and the units. It converts to P_s or g by a named rule (Eq. (121); PARALLEL rev. 1.4, Section 8.7) before any printed prefactor is applied.
2. A single-component PSD is converted to the total transverse variance only under a declared axisymmetry assumption. A single measured PSD does not separate slab and 2D power (PARALLEL rev. 1.4, Section 8.4).
3. Where a source states no normalization or sidedness (Chen, Wijsen, Laitinen, Zhang, Engelbrecht, Subashchandar), the factor is a user decision (U-10, U-21 to U-24, U-40). The library does not choose one.
4. The frame of μ is part of the model. A wave-frame coefficient (Afanasiev) is used with the frame transformation of Section 8.
5. If k_res lies outside the supplied band, the adapter applies a declared tail policy. Missing coverage is not treated as zero power (PARALLEL rev. 1.4, Section 8.4).

Model IDs: MEAS-TAYLOR, MEAS-CHEN24, MEAS-WIJSEN19, MEAS-LAITINEN16, MEAS-ZHANG23, MEAS-DING22, MEAS-AFANASIEV15, MEAS-SOLPACS-GRID, MEAS-SUBASH25-TAYLOR (derivation record DR-06).

## 14. Self-generated waves and shock-associated scattering

Near shocks and in intense events, particles generate part of the Alfvén waves that scatter them. A model of this needs four parts:

- a D_μμ written in terms of the resonant wave intensity;
- a growth rate written in terms of the particle distribution;
- a wave transport equation;
- a rule that couples particle steps to wave updates.

This section gives each part as printed in the sources. Section 8 treats the wave-frame D_μμ forms; they appear here only to fix the spectrum convention of each growth rate. Mean-free-path prescriptions near shocks (Bohm, the steady-state λ of Afanasiev et al. 2015, PARASOL, M-FLAMPA) are in MFP v1.2, Section 8, and are not repeated. As PARALLEL rev. 1.4, Section 12.4, specifies, the coefficient library reads the wave state but does not grow the waves.

### 14.1 Sources, frames and verification

The frames are:

- the plasma frame;
- the wave frame, i.e. the rest frame of the waves of one propagation direction, moving at ±V_A along B relative to the plasma;
- the shock frame;
- the fixed frame.

All sources scatter particles elastically in a wave frame. In the plasma frame the scattering changes the particle energy, and the waves of that direction gain or lose that energy. Vainio & Kocharov [P88] print −ΔE = ∓V_A p_± Δμ_± for a scattering Δμ_± in the ± wave frame (unnumbered, Full text ×2).

| Source | Lineage | Used here | μ, p frame | Verification |
|---|---|---|---|---|
| [P88] (pre-2006) | Monte Carlo, ± waves | ν_±, Γ_±, kernel, wave update | ± wave frame | Full text ×2; Eq. A.2 complete in one rendering |
| [P87] (pre-2006) | analytic | energy exchange, Γ(k) | wave frame | Full text ×2 (Eqs. 2, 5; legible parts of 1, 3, 4, 6) |
| [P126, P89] | CSA | resonance, Γ_w, wave equation, coupling | wave frame; tracking in shock frame | Full text ×2 (Battarbee2011 Eqs. 3–5 ×1) |
| [P78, P37, P29] | SOLPACS; PARASOL | full-resonance D_μμ, steady state, shock-frame equations | Alfvén-wave rest frame | Full text ×2 (Afanasiev2015 Eq. 14 ×1) |
| [P127] | semi-analytic foreshock | λ(x, v), validity limits | — | Full text ×1 |
| [P93] (pre-2006); [P38]; [P82] | Ng–Reames finite differences | four modes, wave action | outward-wave frame | 2003 ×1, garbled; 2008 preprint, several fetches (Section 14.4); 2012 ×1 |
| [P83] | particles + wave packets | SDE step, energy-conserving growth | wave rest frames | Full text; Eqs. 27–43 PDF only (2 fetches) |
| [P91] (pre-2006) | Monte Carlo, relativistic shocks | two-frame scattering | scatterer frame | Full text ×2 (App. A ×1) |
| [P13] | finite differences, electrons | γ_{n,j}, spectrum update | not stated | Eqs. 4–22 ×2; Eqs. 24–39 PDF only (2 fetches) |
| [P36] | M-FLAMPA theory | branch D^σ_μμ, γ_σ | branch rest frame | PDF only; fraction bars inferred |

### 14.2 Spectrum conventions shared by D_μμ and the growth rate

A growth rate is consistent only with the D_μμ of the same paper, because the two share one spectrum. The printed forms below are wave-frame coefficients unless the table of Section 14.1 says the frame is not stated, and use Gaussian units unless SI is stated.

| Source | Printed D_μμ | Spectrum | Resonance |
|---|---|---|---|
| Vainio2003 Eq. 2 | (π/4) ω_c (1−μ²) \|k_r\| W(k_r)/U_B, U_B = B²/8π, ω_c = ω_cp/γ | W dk = wave energy density in [k, k+dk]; sidedness not stated | k_r = −ω_c/(vμ) |
| Battarbee2011 Eqs. 8, 9 | 2π² ω_ci (1−μ²)\|k_r\|W(k_r)/B² = (π/2)Qω_cp(1−μ²)\|k_r\|W/(Aγm_pn_pv_A²) | W = parallel wave energy per unit k | f_res = f_ci(u_sw + v_A)/(γv), μ dropped (Eq. 1) |
| Afanasiev2015 Eq. 2; Afanasiev2023 Eq. 4 (SI) | (π/2) Ω \|k_res\| I_w,res (1−μ²)/B², Ω = Ω_0/γ | δB² = ∫_{−∞}^{∞} I_w dk; sign of k = polarization (Eq. 5) | k_res = Ω/(vμ) |
| Afanasiev2015 Eqs. 3–4 (CSA) | (ν/2)(1−μ²), ν = πΩ²I_w,res/(vB²) | as Eq. 2 | k_res = Ω/v |
| VainioKocharov2001 Eq. 1 | ½(1−μ²)ν_±, ν_± = (π/4)ΩU_±/U_B | U_± = total wave energy per logarithmic bandwidth | fixed k = mΩ_p/p_0 |
| StraussLeRoux2019 Eqs. 24–25 | (πv²(1−μ²)/(B_0²r_L²))∫_0^∞dk Σ_{n,j}R_{n,j}g_{n,j} | δB² = 2π∫_0^∞Σ_{j,n}g_{n,j}dk (Eq. 9) | Lorentzian, τ⁻¹ = α\|V_Ak\| |
| Borovikov2019 Eq. 4.31 (SI) | (πω_ci/(2B²/μ_0))(1−μ_σ²)(ω_ci/(v_σ\|μ_σ\|))I_σ(ω_ci/(v_σ\|μ_σ\|)) | (δB)²/μ_0 = Σ_σ∫_0^∞d\|k\|I_σ | δ(\|k\|μ_σv_σ − g_σω_ci) |

Derived relations (checked numerically; see derivation record DR-14):

1. Vainio2003 Eq. 2 equals Battarbee2011 Eqs. 8 and 9, using U_B = B²/8π and v_A² = B²/(4πm_pn_p).
2. These equal Afanasiev2015 Eq. 2 if I = 4πW. That holds if W is the total wave energy at equipartition and W and I_w have the same k-sidedness.
3. At |μ| = 1, the D^±_μμ of Vainio & Kocharov (fixed resonant k) equals Vainio2003 Eq. 2 if U = 2kW, and half of it if U = kW (U-30). The authors' footnote says that the ν of Vainio et al. (2000) "is defined as having half the value of" their ν_±.

No source states the conversion between these conventions, so relations 2 and 3 are conditional. The resonance sign differs between papers: k_r = −ω_c/(vμ) with δ(k + ω_c/(vμ)) in Vainio2003 and Battarbee2011, but k_r = +Ω_0/(γvμ) with δ(k − k_r) in Afanasiev2015/2023/2025. The difference matters when I(k) is stored with signed k.

### 14.3 Energy exchange and growth rates

**Energy route.** Vainio2003 (Eq. 1) and Battarbee2011 (Eq. 7), Full text ×2:

$$
\frac{dU_w}{dt}=-\int d^3p\;V_A\,p\,\frac{\langle\Delta\mu\rangle}{\Delta t}\,f(\mathbf r,p,t),\qquad \frac{\langle\Delta\mu\rangle}{\Delta t}=\frac{\partial D_{\mu\mu}}{\partial\mu},
\tag{127}
$$

with ΔE_w = −V_A p Δμ per scattering. The Vainio2003 PDF rendering shows Δμ/Δt, and the HTML shows ⟨Δμ⟩/Δt.

**Integral growth rate** (Vainio2003 Eq. 3):

$$
\Gamma(k)=\frac{1}{W}\frac{dW}{dt}=\frac{\pi}{2}\,\frac{\omega_{cp}}{n V_A}\int d^3p\;v\,(1-\mu^2)\,|k|\,\delta\!\left(k+\frac{\omega_c}{v\mu}\right)\frac{\partial f}{\partial\mu},
\tag{128}
$$

where:

- ω_cp = eB/(m_pc), ω_c = ω_cp/γ;
- n is the background proton density;
- p, μ and v are measured in the wave frame.

Printed variants:

- Battarbee2011 Eq. 10 gives the per-species form, with Qω_cp, n_p(r), ω_ci = Qω_cp/(Aγ), d³v and F.
- Afanasiev2023 (Eq. 5) and Afanasiev2025 (Eq. 15) print δ(k − k_r) with k_r = Ω_0/(γvμ) and Ω_0 = eB/m_p (SI). They cite Vainio (2003).

Afanasiev2023 adds that "appreciable wave growth at a given point in space requires a finite number of particles to pass through that point". It prints no damping or saturation term. Model ID SGW-GROWTH-INT.

**Streaming form.** Vainio2003 integrates by parts in μ and replaces δ(k + ω_c/(vμ)) by ½δ(|k| − ω_c/v). The paper calls this an order-of-magnitude approximation (Eq. 4). The result is Eq. 5:

$$
\Gamma(\mathbf r,p,t)=\frac{\pi}{2}\,\frac{\omega_{cp}\,p\,S(\mathbf r,p,t)}{n V_A},\qquad S=2\pi\int_{-1}^{1}d\mu\;v\,\mu\,p^2 f,\qquad p=\frac{m_p\,\omega_{cp}}{|k|}.
\tag{129}
$$

S is the net streaming per unit area, momentum and time. The paper then treats S as a fixed-frame quantity.

Printed variants:

- Battarbee2011 Eq. 12, per species: Γ_{w,i} = (π/2)Q_iω_cp v_r S_v/(n_pv_A). Here S_v = 2πv²∫dμ vμF "in the Alfvén-wave frame", evaluated at v_r = ω_ci/|k|.
- Battarbee2010, frequency form: σ = π²f_cp p S_p/(n v_A), with p = m_pVf_cp/f.

These forms drop the μ-dependence of the resonance. Model ID SGW-GROWTH-STREAM.

**Two directions** (VainioKocharov2001 Eq. 2, Full text ×2; fixed resonant k):

$$
\Gamma_\pm=\pm\frac{1}{U_\pm}\int d^3p_\pm\,\nu_\pm V_A\,p_\pm\,\mu_\pm f_\pm=\pm\frac{\pi\,\Omega_p\,S_\pm}{2\,n_0V_A},\qquad S_\pm=\int d^3p_\pm\,v_\pm\mu_\pm f_\pm ,
\tag{130}
$$

where S_± is the particle flux in the ± wave frame. Eq. 3 integrates Γ_± along x ∓ V_A t′. Model ID SGW-GROWTH-VK01.

**Helicity-resolved rate** (StraussLeRoux2019 Eq. 37 = A12; PDF, 2 fetches). The rate applies to forward (j = +1) and backward (j = −1) waves of polarization n (−1 right-handed, +1 left-handed). It assumes the mono-energetic distribution f = p_0 f(r, μ, t) δ(p − p_0).

$$
\gamma_{n,j}=j\,\frac{\pi^2}{4}\,\frac{p_0^4}{k^{\rm res}_\parallel}\,\frac{\Omega^2}{\mathcal B}\,\frac{V_A}{v}\,(1-\mu_{\rm res}^2)\,\frac{\partial f}{\partial\mu}\Big|_{\mu_{\rm res}},\qquad \mathcal B=\frac{B_0^2}{8\pi},
\tag{131}
$$

with the definitions:

- k^res_∥ = nΩ/(V_A − vμj) (Eq. 29);
- μ_res = (ω − nΩ)/(jvk^res_∥) (A13);
- ω = V_Ak_∥ in the plasma frame.

γ is an amplitude rate, (1/g)∂g/∂t ≈ 2γ (Eq. 8). Eqs. 25 and 29 do not state which species Ω refers to. Model ID SGW-GROWTH-SLR19.

**Not implementable from the texts read.**

- Borovikov2019, Eqs. 4.15–4.17 and 5.6 (γ_σ for four branches). The legible Eq. 4.17 is γ_σ = 2π∫p²dp dμ K_σ v(1−μ²)[∂f/∂p + (1/(m_iV_σ) − μ/p)∂f/∂μ], with K_σ ∝ δ(|k|(μv − V_σ) − g_σω_ci). The reader inferred the constants of K_σ and the fraction bars of Eq. 5.6.
- Ng–Reames γ_σ (NgReamesTylka2003 Eq. 19; NgReames2008 Eq. 7). The factor G_f = P³[E(1 − μV_σ/v)]⁻²[∂_μf − (V_σ/v)(μ∂_μf − …)] is only partly legible. The 2003 Appendix B derivation was not reached.

### 14.4 Wave transport, damping and saturation

| Source | Printed wave evolution | Transport |
|---|---|---|
| Vainio2003 Eq. 6 | W = W_0 exp∫_0^tΓdt′ | none ("wave propagation explicitly neglected") |
| VainioKocharov2001 App. A | every Δx/V_A, U_±(X_i) ← U_±(X_{i∓1}) | one-cell shift; U_+ = U_0 at inflow; outflow lost |
| Battarbee2011 Eq. 6 | Eq. (132) | WKB at V = u_sw + v_A, diffusion in frequency |
| Afanasiev2023 Eq. 2; Afanasiev2025 Eq. 13 | Eq. (133) | advection in the shock frame |
| NgReames2008 Eq. 5 | ∂_tΨ_σ + ∂_x[(W_σ − V_sh)Ψ_σ] + η∂_η(V^σ_ηΨ_σ) = γ_σΨ_σ | WKB in x and η = k/B; Ψ_σ = 2I_σ η\|W_σ\|/(η_0V_A), V^σ_η = 2W_σ/r − dW_σ/dr |
| Sandroos2014 Eqs. 6–7 | packets at V_p ± V_A; d log λ^±/dt = d(V_p ± V_A)/dr | Lagrangian; shock transmission Eqs. 16–17 |
| StraussLeRoux2019 Eq. 38 | g(t) = g(t_0)exp(2∫γdt′) | none; cascade and dissipation of the new waves neglected |

CSA (Battarbee2011 Eq. 6, Full text ×2; model ID SGW-WAVE-CSA):

$$
\frac{\partial\tilde P}{\partial t}+V\frac{\partial\tilde P}{\partial r}=\Gamma_w\tilde P+\frac{\partial}{\partial f}\!\left(D_{ff}\frac{\partial\tilde P}{\partial f}\right),\qquad \tilde P=\frac{V^2}{B\,v_A}P,\qquad V=u_{sw}+v_A,\qquad D_{ff}=\frac{V}{r_\oplus}f^{8/3}f_b^{-2/3}.
\tag{132}
$$

Here P(f) is the power per unit fixed-frame frequency, f_b = 1 mHz and r_⊕ = 1 AU. D_ff is ad hoc: it makes unenhanced spectra approach a Kolmogorov form above f_b. Battarbee2010 prints the advection speed as a lower-case v, a symbol it does not define (as printed).

SOLPACS (Afanasiev2023 Eqs. 1–3; Afanasiev2025 Eqs. 12–13; Full text ×2; model ID SGW-WAVE-SHOCK):

$$
\frac{\partial f}{\partial t}+\left[v\mu+(1-M_A)V_A\right]\frac{\partial f}{\partial x}=\frac{\partial}{\partial\mu}\!\left(D_{\mu\mu}\frac{\partial f}{\partial\mu}\right),\qquad
\frac{\partial I}{\partial t}+(1-M_A)V_A\frac{\partial I}{\partial x}=\Gamma I ,
\tag{133}
$$

where:

- x is the distance along B from the shock, with x > 0 upstream;
- v and μ are measured "in the rest frame of Alfvén waves propagating at speed V_A wrt. the ambient plasma";
- M_A = u_1/V_A = M*_A/cos θ_Bn (Eq. 3), with u_1 the upstream bulk speed in the de Hoffmann–Teller frame;
- Γ is Eq. (128) with δ(k − k_r).

**Damping, floors and caps as printed.**

- Ng, Reames & Tylka (2003): a phenomenological rate of change of the ambient waves (0.01 hr⁻¹), and an upper cutoff κ_c with k⁻⁴ beyond it (Eq. 28, garbled). They set SEP-driven growth to zero downstream.
- Vainio & Kocharov: a minimum wave-energy floor.
- Virtanen & Vainio, Eq. 23: k′I(k′) = B_0² above γ′_0 (a Bohm-type cap as printed).
- Strauss & le Roux: dynamical damping τ⁻¹ = α|V_Ak| (α = 1), only inside the resonance function.
- Sandroos et al.: downstream waves are "very efficiently damped by particles".

**Saturation cap: the readings disagree.**

- The reading behind this section (three fetches of the NgReames2008 preprint) found no explicit cap I ≤ B²/(3πk), and none of the other texts in this section prints one.
- A one-fetch reading of the same preprint, and MFP v1.2, Section 8.2, record I_σ(k) ≤ I_sat(k) = B²/(3πk).

The disagreement is not resolved here, and the cap has no model ID (Section 19, C-7 and U-50).

### 14.5 Steady-state waves upstream of a parallel shock

Afanasiev2015 Eq. 7 (Full text ×2) gives a Bell-type steady state (after Bell 1978, as cited; not read):

$$
I_w(x,k)=\frac{1}{3\pi}\,\frac{\Omega B^2|k|^{-3}}{(u_1-V_A)(x+x_0)} .
\tag{134}
$$

The companion equations in Afanasiev2015 are:

- f = f_sx_0/(x + x_0) (Eq. 6);
- f_s = (σε_injn_0/(4πp_inj³))(p/p_inj)^{−σ} (Eq. 8), with σ = 3r_c/(r_c − 1) and r_c = V_1/V_2 the scattering-centre compression ratio;
- x_0 = (2/(πε_injσ))(V_A/Ω_0)(p/p_inj)^{σ−3} (Eq. 9);
- p_c^{σ−3} = p_inj^{σ−3}[1 + (π/2)(ε_injσ(M_A − 1)/r_c)Ω_0t_s] (Eq. 13).

The symbols u_1 and Ω denote the shock-frame wind speed and Ω_0/γ. The mean free path of Eq. 12, λ = 3(u_1 − V_A)(x + x_0)/v, is SHOCK-AFANASIEV15 of MFP v1.2, Section 8.2, Eq. 27. Model ID SGW-BELL-A15 denotes Eq. (134) used as a wave provider.

Vainio et al. [P127] (Eqs. 4–7) print the same structure. They add:

- the validity condition λ/L ≪ 3(u_1 − V_A)/v (Eq. 9), with L = −B/(∂B/∂R) = R/2;
- λ(x, E) = (3(u_1 − V_A)/v)(x + x_0)e^{x/x_feb} (Eq. 22), where x_feb is given by Eq. 24 (x_1 is undefined in the visible text);
- the remark that the CSA code neglects μ in the resonance k = mΩ/(pμ).

Derived (checked; fixture F-SGW-03):

1. Inserting Eq. (134) at k = Ω/v into the CSA rate (Afanasiev2015 Eq. 4) gives ν = v²/[3(u_1 − V_A)(x + x_0)] and λ = v/ν = Eq. 12 (this re-checks MFP F-SH-02).
2. Inserting it into the full-resonance Eq. 2 gives D_μμ = v²μ²(1−μ²)/[6(u_1 − V_A)(x + x_0)]. This D_μμ vanishes like μ², so the integral of Eq. (23) diverges. Excluding |μ| < δ gives λ_δ = (3/2)[(1 − δ)²/δ] × 3(u_1 − V_A)(x + x_0)/v. Eq. 12 therefore belongs to the pitch-angle-independent resonance. With the full resonance, λ depends on the gap treatment.
3. The printed SOLPACS values n_0 = 3.6 × 10¹² m⁻³ and B_0 = 3.4 × 10⁻⁵ T give Ω_0/V_A = 8.33 × 10⁻³ m⁻¹. This is consistent with the printed k_max = 10⁻² m⁻¹, described as ≈ Ω_0/V_A.

### 14.6 Coupling particle steps and wave updates

**(a) Vainio & Kocharov (2001), Appendix A** (model ID SGW-COUPLE-VK01). The step order below is this review's reading of the appendix.

1. Particles carry weights w_j normalized to the injection rate. The grid has N = 150 cells.
2. Move each particle: x_j ← x_j + v_jμ_jδt.
3. Scatter in each wave frame in turn: Lorentz-transform to the + frame, apply Eq. (135) below, and transform back; then repeat in the − frame. ϑ² is exponentially distributed with mean 2ν_±δt, and φ is uniform.
4. Update the wave energy from the mean plasma-frame energy loss by Eq. (136) below (Eq. A.2; complete only in the PDF rendering). The paper says this is more stable than summing individual particle changes.
5. Apply the wave-energy floor.
6. Shift the waves by one cell every Δx/V_A.

$$
\mu_{j\pm}\leftarrow\mu_{j\pm}\cos\vartheta_{j\pm}+\left(1-\mu_{j\pm}^2\right)^{1/2}\sin\vartheta_{j\pm}\cos\varphi_{j\pm},
\tag{135}
$$

$$
U_\pm(X_i)\leftarrow U_\pm(X_i)+\frac{1}{A\,\Delta x}\sum_{x_j\in I_i}w_j\,\nu_\pm(X_j)\,\delta t\,V_A\,p_{j\pm}\,\mu_{j\pm}.
\tag{136}
$$

As printed, the update has the same sign for both directions; with Eq. (130) the backward-wave term would carry the opposite sign unless μ_{j−} is measured along the backward-wave direction, which the text read does not state (D-126).

Time step: ⟨ϑ²_±⟩ = 2ν_±δt < 0.2 for all particles.

Derived (fixture F-SGW-01): to first order, the kernel gives ⟨Δμ⟩ = −μνδt and ⟨(Δμ)²⟩ = (1−μ²)νδt. These are D_μμ = (ν/2)(1−μ²) (Eq. (24)) with the Itô drift ∂D_μμ/∂μ, as the paper states. At ⟨ϑ²⟩ = 0.2, the exact first moment is 3.3 % below its first-order value for every μ; the second moment depends on μ: 2.1 % below at μ = 0.6 (F-SGW-01), 12 % below at μ = 0, and above the first-order value for |μ| > 0.636 (derived from Eq. (169)). The kernel rotates the velocity vector, so |μ| ≤ 1 holds without a boundary rule.

**(b) CSA** (model ID SGW-COUPLE-CSA).

- Particles are tracked in a shock-attached frame on a logarithmic grid that extends to 300 R_⊙ ahead of the shock.
- Streaming is "tracked via particle movement between simulation cells, weighted to take into account the movement of the grid" (Battarbee2010). Each crossing of a cell boundary is weighted by (v_wμ_w)/(v_sμ_s), where w is the wave frame and s the shock frame (Battarbee2011).
- Amplification is computed at cell boundaries over 1.6 ms intervals.
- Wave power is convected by "a Van Leer flux-limited Lax-Wendroff scheme" (Battarbee2010).
- Downstream, particles are followed for 2λ and then returned with probability P_ret = ((v_w − u_2)/(v_w + u_2))² (Eq. 3). A run stops when the cumulative return probability falls below 0.1 %.
- No particle time step is printed.

**(c) SOLPACS** (Afanasiev2015 Sect. 3; model ID SGW-COUPLE-SOLPACS, partly specified).

- The spectrum is stored at nodes k_j, with I_w(k) = I_w,j(k/k_j)^{−q_j} between nodes, and μ_j = Ω/(υk_j).
- ν(μ) = ν_{0,j}|μ|^{q_j−1}, with ν_{0,j} = (πΩ/B²)(υ/Ω)^{q_j−1}I_w,j|k_j|^{q_j} (×1). λ follows from Eq. 14 (×1).
- Near μ = 0 the resonant k exceeds k_max ≈ Ω_0/V_A. The gap is filled linearly: ν̃(μ) = bμ + b_1 on (−μ_min, μ_min), with μ_min = Ω/(υk_max). The fill is printed in μ, not |μ|, and b and b_1 are not given.
- ν(μ) is evaluated from the current spectrum "at each time step". The method "requires extremely small time steps" near μ = 0.
- The growth numerics are in Afanasiev & Vainio [P128], which was not accessible.

**(d) Ng–Reames.** Not implementable from the texts read.

- Finite differences with operator splitting, locally one-dimensional.
- Amplification at every coupling step; wave transport on a larger step; a dynamic coupling interval limited by growth and ∂F/∂μ (2003, garbled).
- μ_i = 0.05(i + ½), i = −20 … 19 (2008).
- Transport in the outward-wave frame, with inward waves at V_σ = −2V_A and resonance k = B/[P(μ − V_σ/v)]. D_μP and D_PP are retained, but the matrix of Eq. 1 is not legible.

**(e) Sandroos et al. (2014)** (model ID SGW-COUPLE-SANDROOS14). The order of operations per step, as reported:

1. Choose the global time step.
2. Map wave packets to the mesh and convert to intensity.
3. Scatter particles in the wave rest frames, accumulating growth factors.
4. Apply growth.
5. Move particles, then packets.
6. Apply boundaries, injection and particle splitting.

The pitch-angle step (Eqs. 32–36, PDF, 2 fetches) is

$$
\tilde\mu_{\rm trial}=\tilde\mu_t+a\,dt+b\sqrt{dt}\,N,\qquad
\tilde\mu_{t+1}=\tilde\mu_t+\left[\theta a_1+(1-\theta)a_0\right]dt+\left[\eta\, b(\tilde\mu_t)+(1-\eta)\,b(\tilde\mu_{\rm trial})\right]\sqrt{dt}\,N,
\tag{137}
$$

with:

- a = ∂D_μμ/∂μ and b² = 2D_μμ;
- a_0 = (1−η)a(μ̃_t) and a_1 = (1−η)a(μ̃_trial), as extracted (Eqs. 34–35);
- the first term on the right of Eq. 36 is printed μ̃₀ and read as μ̃_t (D-120); the same scheme is Eq. (160) of Section 16.3;
- θ = 0 and η = 0.5;
- N a unit normal deviate.

Further printed details:

- Time step (Eq. 39): dt_max = min(Δμ_max/a_max, Δμ_max²/b_max²), with N_sub = dt_sim/dt_max ≤ 10.
- Gap: μ̃_min = 0.025 in the resonant wavelength λ_r = (2π/Ω)Ṽμ̃ (Eq. 31).
- Reflection: printed as "μ̃ → 1 − μ̃" (as reported).
- Energy transfer:
  - each scattering transfers ΔU = mV_AṼΔμ̃;
  - a macroparticle of weight W (Eq. 40) changes each fully crossed wavelength cell by ΔU_w,l = −ΔuΔλ_l, with Δu = WmV_AΩ/(2π) (Eq. 41);
  - this change is smoothed twice per step with D = 0.25 (Eq. 42);
  - particles are split by a criterion R < 0.01 (Eq. 43, as reported).
- The text read does not give the rule for choosing the wave frame, nor the frame transformations.

Derived (fixture F-SGW-02): suppose the same N is used in both stages. The mean step is then (1−η)(a + b∂b/∂μ)dt = 2(1−η)a dt + O(dt²), since b∂b/∂μ = a. The variance is b²dt + O(dt²). With η = 0.5 the mean is a dt, as Eq. (25) requires. With η = 0.3 it would be 1.4a dt. The extracted form is therefore consistent only at the published η.

**(f) Strauss & le Roux (2019)** (model ID SGW-COUPLE-SLR19). At each step:

1. Compute ∂f/∂μ numerically.
2. Evaluate γ_{n,j} from Eq. (131).
3. Update the spectrum with the source's Eq. 39, given below as Eq. (138).
4. Recompute the transport coefficients for the next step.

$$
g^{\rm slab}_{n,j}(k^{\rm res},t_0+\Delta t)=g^{\rm slab}_{n,j}(k^{\rm res},t_0)\,\exp\!\left[2\gamma_{n,j}(k^{\rm res},t_0)\,\Delta t\right].
\tag{138}
$$

**(g) Virtanen & Vainio (2005), Appendix A** (PDF only).

- The particle "is first scattered in the forward-wave frame, and then immediately again in the backward-wave frame", after a Lorentz transformation to each frame.
- The two frames are weighted by ω_±, with ω_+ + ω_− = 1, computed from the amplified-wave ratios. The formula for ω_± is not given.
- The kernel is the same as Eq. (135) (their Eq. 24).
- The time step is printed as "Δt = a ν_e(γ)" with a ≲ 0.01 (as extracted).

### 14.7 Fitted shock-region scattering without wave evolution

PARASOL [P29] does not evolve waves in its pusher (PARADISE). Self-generated waves enter only through a time-independent upstream λ fitted to SOLPACS runs, and through the shock source. These are SHOCK-PARASOL in MFP v1.2, Section 8.3. The pusher uses the ε-form (Eq. (29)) with ε = 0.048 and D_0 fixed from λ∥ (Eqs. 41–42). It solves the conservative focused transport equation for j = p²f (Eq. 37). Vainio2014 Eq. 22 (Section 14.5) is another closed-form λ of this kind.

### 14.8 Implementation requirements

**Implementation requirement.** A self-generated-wave model stores its growth rate together with the D_μμ that shares its spectrum. It also stores the spectrum convention:

- energy or magnetic intensity;
- one- or two-sided;
- per k or per ln k;
- the unit system;
- the sign convention for k.

Conversions between the conventions of Section 14.2 need a declared assumption. Relations 2 and 3 there are not defaults.

**Implementation requirement.** Scattering by one wave population is done in that population's rest frame, and the plasma-frame energy change goes to the same population. When both propagation directions are present, the order and weighting of the two frames are configuration items. The published choices are: both frames every step (Vainio & Kocharov), sequential with unprinted weights (Virtanen & Vainio), or an unstated rule (Sandroos et al.).

**Implementation requirement.** A full-resonance D_μμ from an evolving spectrum needs an explicit μ ≈ 0 treatment. The published treatments are:

- the SOLPACS linear fill (b and b_1 not printed);
- μ̃_min = 0.025 (Sandroos et al.);
- ε (PARADISE);
- broadening (Ng–Reames; form not legible).

The treatment and its parameters are configuration inputs, and none is a default.

**Implementation requirement.** Floors, caps and phenomenological damping are separate, named contributions. A saturation cap I ≤ B²/(3πk) is applied only when configured.

**Implementation requirement.** Each coupling scheme is stored with its printed time-step rule:

- ⟨ϑ²⟩ < 0.2 (Vainio & Kocharov);
- 1.6 ms growth intervals (CSA);
- Eq. 39 with N_sub ≤ 10 (Sandroos et al.).

The coefficient layer reads the current wave snapshot (MFP v1.2, Section 21.6).

## 15. Pitch-angle-resolved GCR transport models

### 15.1 Models found

Most heliospheric GCR modulation codes solve the Parker equation (MFP v1.2, Section 9). The pitch-angle-resolved GCR models found are few:

- a 3-D focused-transport modulation model compared with the Parker equation [P14];
- a 3-D Forbush-decrease model [P18];
- a 1-D model of near-Sun GCR anisotropy [P4];
- local heliopause models [P46, P129];
- a TeV model in which D_μμ shapes only the interstellar boundary distribution [P130].

The GLE model of Kubo et al. [P16] and the benchmark of Effenberger & Litvinenko [P2] are included because GCR implementations use them.

| Model | Particles, method | D_μμ form | Amplitude fixed by | Near μ = 0 | At μ = ±1 |
|---|---|---|---|---|---|
| vandenBerg2026 | GCR p, 3-D, steady, backward SDE | D_I(1−μ²) or D_A(1−μ²)(\|μ\|/(1+\|μ\|)+ε) | κ∥(P), Eqs. 3d, 13 | ε = 0.045 | reflection (Eq. 23) |
| QinWu2020 | GCR p, 3-D, backward SDE | slab QLT × (μ^{s−1}+h)(1−μ²) | turbulence σ, l_slab | h = 0.01 | not stated |
| Strauss2022 | GCR p, 1-D, finite differences | D_0(1−μ²)(\|μ\|^{q−1}+H) | λ∥ | H = 0.05 | zero flux (code) |
| Kubo2015 | GLE p, 1-D, backward SDE | D_0vp^{q−2}(\|μ\|^{q−1}+h)(1−μ²) | λ_rr | h = 0.2 | not stated |
| Florinski2013 | 1 MeV p, gyrophase-resolved | D(1−μ²) (sphere) | rate D | — | random walk on sphere |
| StraussFichtner2014 | GCR, 2-D, backward SDE | not printed | λ∥ | — | not stated |
| Zhang2014 | TeV p, LISM boundary only | D_0(1−μ²) | λ∥ | — | — |
| EffenbergerLitvinenko2014 | 1-D benchmark | D_0(1−μ²) | λ_0 = v/(2D_0) | — | reflecting |

### 15.2 van den Berg, Engelbrecht & Strauss (2026)

**Source and verification.**

- A&A 712, A204. Sections 1–5 were read in the PDF text layer, Sections 1–2 also in the HTML.
- Eqs. 1–4: two renderings agree. Section 3 equations: one rendering, read twice independently ("×2"), unless marked ×1.
- Appendices B (SDE verification) and C (anisotropic distribution) were not reached.
- The arXiv identifier is not confirmed.

**Transport equation (Eq. 1).** The 3-D focused transport equation in conservative form is

∂f/∂t + ∇·(v_gc f) + p⁻²∂_p[p²(dp/dt)f] + ∂_μ[(dμ/dt)f] = ∂_μ(D_μμ∂_μf) + ∇·(D_⊥·∇f).

- v_gc = μvb̂ + u + v_d(μ); dp/dt and dμ/dt are given by Eqs. 2b–2c.
- p and μ are measured in the flow frame.
- D_⊥ = D_⊥(μ)(I − b̂b̂). The mixed coefficients D_μ⊥ are neglected.
- The runs are for solar minimum, steady state, without drifts.
- Normalization (Eqs. 3d–3e): κ∥ = (v²/8)∫(1−μ²)²/D_μμ dμ (Eq. (23)) and κ⊥ = ½∫D_⊥(μ)dμ.
- Focusing (Eq. 11, ×2): L⁻¹ = ∇·b̂ and L′⁻¹ = ∇·b̂ − (2/v²)b̂·[∂u/∂t + (u·∇)u].

**Pitch-angle diffusion** (Eqs. 12a–12b, ×2):

$$
D^{\rm iso}_{\mu\mu}(\mu)=D_I\,(1-\mu^2),\qquad
D^{\rm aniso}_{\mu\mu}(\mu)=D_A\,(1-\mu^2)\left(\frac{|\mu|}{1+|\mu|}+\epsilon\right),\qquad \epsilon=0.045 .
\tag{139}
$$

What the paper states about the two forms:

- The paper attributes Eq. 12b to Agüeda et al. (2008) [P61] and uses it "because its derivative is well behaved around µ = 0".
- "ϵ = 0.045 parametrises dynamical effects". Eq. 12b "is consistent with dynamical QLT".
- Eq. 12a "would be expected only in strong turbulence or in magnetostatic QLT".

Eq. 12b is the ε-form (Eq. (29)) with D_A = ν_0/2. Model IDs: GCR-PA-VDB26-ISO and GCR-PA-VDB26-ANISO.

**Amplitudes** (Eq. 13, ×2; both readings give the logarithm's argument as (1+2ε)/ε, but the fraction layout is not visible):

$$
D_I=\frac{v^2}{6\kappa_\parallel},\qquad
D_A=\frac{v^2}{4\kappa_\parallel(1+\epsilon)^4}\left[\frac16+2\epsilon+\frac52\epsilon^2+\frac23\epsilon^3+(1+2\epsilon)\ln\frac{1+2\epsilon}{\epsilon}\right],\qquad \epsilon\neq0 .
\tag{140}
$$

Derived (fixture F-GCR-PA-01):

- Eq. 3d with Eq. 12b gives D_A = (v²/4κ∥)∫_0^1(1−μ)(1+μ)²/[(1+ε)μ + ε]dμ.
- The residue of the integrand at μ = −ε/(1+ε) is (1+2ε)/(1+ε)⁴, which fixes the coefficient of the logarithm.
- Quadrature and Eq. (140) agree to better than 10⁻³² (relative) for ε = 0.01–0.5, which confirms the reading.
- At ε = 0.045: D_Aκ∥/v² = 0.783198422199096 and D_A/D_I = 4.69919053319457.
- D_A/D_I is the φ(ε) of MFP v1.2, Eq. 14. So Eq. (140) is a closed form of φ, and it reproduces MFP F-PA-02 (φ(0.048) = 4.597262, φ(0.01) = 7.069533).

**κ∥ and κ⊥** (×1, partly tool-reconstructed).

- Eq. 9: the legible fragments ("κ∥ = κ0 vB0/(cB)", "(P/P_r)", "γ + (P_k/P_r)γ", "δ−ρ") are consistent with the NWU law of MFP v1.2, Section 9.2, Eq. 33. The mapping is: κ_0 is a value (not a unit), with a → ρ, b → δ, c → γ and β = v/c.
- Parameters: κ_0 = 1.3 × 10²³ cm² s⁻¹ (printed "= 2.1 au² h⁻¹"; derived 2.0912), B_0 = 1 nT, P_r = 1 GV, P_k = 4 GV, ρ = 0.43, δ = 1.95, γ = 3.
- At P = P_r the bracket is 1 for any grouping. Fixture F-GCR-PA-02 uses this point: λ∥ = 0.195855 au at B = 4.44 nT.
- Eq. 10: κ⊥r = χκ∥ with χ = 0.02, and κ⊥θ = κ⊥r[(ς+1)/2 + ((ς−1)/2)tanh((|θ − 90°| − α_⊙ − θ_s)/θ_t)]. Values: ς = 3, α_⊙ = 13.5°, θ_s = 15°, θ_t = "45°/2π" (as reported).

**Perpendicular pitch-angle dependence** (Eq. 14, ×2; model ID GCR-PA-VDB26-DPERP):

$$
D_\perp(\mu)=h(\mu)\,\kappa_\perp,\qquad h=1\ \text{with}\ D^{\rm iso}_{\mu\mu},\qquad h=3\mu^2\ \text{with}\ D^{\rm aniso}_{\mu\mu}
\tag{141}
$$

The form h = 3μ² represents "non-linear theories". Derived: ½∫h dμ = 1 in both cases, so Eq. 3e returns κ⊥.

**Backward SDE** (model ID GCR-PA-VDB26-BSDE).

- Eq. 15 (×1): the backward Kolmogorov equation, −∂P̃/∂τ = a·∇_qP̃ + ½C:∇_q∇_qP̃ − ŁP̃.
- Eq. 16: dQ = a dτ + B·dW in backward time, with C = BB^T.
- Eq. 17: Euler–Maruyama, Q_{k+1} = Q_k + aΔτ_k + B·Λ√Δτ_k. The √ glyph was visible in the raw text layer in two readings; a third reading did not show it.

For μ (Eq. 38e and text, ×2; consistent with Eq. 36, ×1):

$$
a_\mu=\frac{\partial D_{\mu\mu}}{\partial\mu}-\frac{d\mu}{dt},\qquad B_{\mu\mu}=\sqrt{2D_{\mu\mu}} .
\tag{142}
$$

- Weights (×1): w ← w exp(ŁΔτ) (Eq. 18), with Ł = −(μv/c²)b̂·[∂u/∂t + (u·∇)u] (Eq. 37). The sign relation of this weight to Eq. 15 is not established (D-125).
- Start (×1): M(0) "uniformly distributed on [−1; 1]" and w(0) = 1.

Boundaries:

$$
M\leftarrow{\rm sign}(M)\,(2-|M|)\ \text{if}\ |M|>1\ \text{(Eq. 23)};\qquad R\leftarrow2r_\odot-R,\ M\leftarrow|M|\ \text{if}\ R<r_\odot\ \text{(Eq. 20)} .
\tag{143}
$$

- The inner boundary is r_⊙ = 0.005 au.
- At the outer boundary (90 au), GCRs "enter the modulation cavity isotropically … regardless of their pitch angles".
- The raw text of Eq. 23 reads "M ← sign(M) 2 − M"; no parentheses or bars are visible. Derived: sign(M)(2 − |M|) = 2sign(M) − M for M ≠ 0, so both readings that keep |M| ≤ 1 coincide. The rule is the mirror M → ±2 − M, valid for |M| < 3.

**Adaptive step** (Eqs. 19, 24, 25; the first entries of Eq. 24 ×2, the rest ×1):

$$
\Delta\tau=\min_i\left\{\frac{l_i}{|a_i|},\ \frac{l_i^2}{\sum_jB_{ij}^2}\right\},\qquad l_\mu=\ell\sqrt{1-\mu^2},\qquad \ell=2^{-6},
\tag{144}
$$

- i runs over r, θ, φ, p, μ, with no diffusive entry for p.
- l_r = ℓ min{L_r; λ_rr}, l_θ = ℓλ_θθ/r, l_φ = ℓ min{L_φ; λ_φφ}/(r sin θ), l_p = ℓp.
- Derived: for D^iso with dμ/dt = 0, the diffusive entry is ℓ²/(2D_I), independent of μ. But l_μ/|a_μ| = ℓ√(1−μ²)/(2D_I|μ|) → 0 as |μ| → 1.
- No lower bound on Δτ is stated.

**Setup** (×1).

- 90 r × 22 θ × 30 energy × 20 μ bins.
- 150 (TPE) and 300 (FTE) pseudo-particles per bin.
- Three-point smoothing.
- LIS of Vos (2011), not transcribed.
- B_⊕ = 4.44 nT; V_sw = 400 km s⁻¹ at the equator, doubled over the poles.

**Analytic distribution** (Eqs. 45–48, ×2; definition of ξ ×1; model ID GCR-PA-VDB26-PAD):

$$
\tilde f(\mu)=\frac{e^{G(\mu)}}{\int_{-1}^{1}e^{G(\eta)}\,d\eta},\qquad
G(\mu)=\frac{v}{2L'}\int_{-1}^{\mu}\frac{1-\eta^2}{D_{\mu\mu}(\eta)}\,d\eta,\qquad
\tilde f_{\rm iso}=\frac{\xi\,{\rm csch}\,\xi}{2}\,e^{\xi\mu},\qquad
A_1^{\rm iso}=3\left(\coth\xi-\frac1\xi\right),
\tag{145}
$$

with ξ = λ∥/L′ and "G(µ) = ξ(1 + µ)". Derived (fixture F-GCR-PA-03): Eq. 46 with D_I = v²/(6κ∥) gives G = ξ(1+μ). The resulting f̃ and A_1 = 3∫μf̃dμ match Eqs. 47–48, and A_1^iso → ξ for ξ ≪ 1. The anisotropic case (Appendix C) was not read.

**Result** (abstract). With equal κ∥ and κ⊥ in both equations, the Parker equation overestimates the proton intensity at Earth by about 30 % at low energies and by about 40 % over the poles. FTE intensities are nearly the same for Eqs. 12a and 12b at equal κ∥.

### 15.3 Qin & Wu (2020): Forbush decrease

**Source and verification.**

- arXiv:2010.07662 v2, "Submitted to ApJ"; publication not confirmed.
- D_μμ, λ∥ and the transport equation: Full text ×2. The equation numbers differ between versions (D_μμ is Eq. 15 in the PDF and Eq. 17 in ar5iv).
- Details are shared with paper I [P58].

**Transport equation.** Eq. 2 is the non-conservative 3-D focused transport equation with convection, cooling and perpendicular diffusion. It is solved by the "time-backward Markov stochastic process method" of Zhang [P131], with L = (b̂·∇ln B_0)⁻¹.

**Pitch-angle diffusion, as printed:**

$$
D_{\mu\mu}(\mu)=\left(\frac{\delta b_{\rm slab}}{B_0}\right)^2\frac{\pi(s-1)}{4s}\,\frac{v}{l_{\rm slab}}\left(\frac{R_L}{l_{\rm slab}}\right)^{s-2}\left(\mu^{s-1}+h\right)(1-\mu^2),
\tag{146}
$$

- R_L = pc/(|q|B_0), with B_0 the local field; s = 5/3.
- h = 0.01 was "introduced to model the non-linear effect of pitch-angle diffusion at μ = 0".
- The LaTeX source has \mu^{s-1} without absolute-value bars, here and in paper I (as printed). For s = 5/3 the real reading (μ²)^{1/3} equals |μ|^{2/3}. But a floating-point power of a negative base is undefined, so using |μ| is a user decision.
- λ∥ = (3v/8)∫(1−μ²)²/D_μμ dμ (Eq. 17).
- κ⊥ (Eq. 16, ×1) = (v/3)[(δb_2D/B_0)²√(3π)((s−1)/(2s))(Γ(s/2+1)/Γ(s/2+1/2))l_2D]^{2/3}λ∥^{1/3}(I − b̂b̂). It therefore depends on the λ∥ from Eq. (146).

**Parameters** (×1).

- Turbulence: σ = √(δb_slab² + δb_2D²)/B_0, with (δb_slab/B_0)_i = (√5/5)σ_i and (δb_2D/B_0)_i = (2√5/5)σ_i.
- Regions: σ = 0.3 (solar wind), 1.6 (sheath, i.e. δb > B_0; QLT is used as printed) and 0.1 (magnetic cloud). The sheath value lies outside the weak-turbulence assumption of quasilinear theory; the source applies the QLT form there, and this document records that application without taking it as evidence that the QLT form holds at δb > B_0 (an empirical extrapolation of the source; registry field `admissibility`).
- Lengths: l_slab = 0.025 au, l_2D = 0.0096 au.
- Wind and field: V_sw = 450 km s⁻¹; B = 5 nT at 1 au (radial component 3.62 nT).
- Domain: R_in = 0.05 au, R_out = 10 au.
- Source (Eq. 3): f_s = f_0p_0^{2.6}/[p(m_0²c² + p²)^{1.8}].
- Oulu effective energy: 6.54 GeV.
- Pseudo-particle number, time step and μ boundary are not stated.

**Derived** (fixture F-GCR-PA-04).

- With |μ|, Eq. (146) is Eq. (26) with q = s and H = h. Its amplitude equals the Qin & Wang (2015) form (MFP v1.2, Section 6.4) with k_min = 1/l_slab.
- Hence λ∥ = 3s I(s,h) l_slab (R_L/l_slab)^{2−s}/[π(s−1)(δb_slab/B_0)²], with I(5/3, 0.01) = 2.135993.
- At equal R_L/l_slab, the sheath λ∥ is (0.3/1.6)² = 0.0352 times the solar-wind value.

Model IDs: GCR-PA-QINWU20 (|μ|) and GCR-PA-QINWU20-PRINTED.

### 15.4 Strauss, van den Berg & Rankin (2022)

**Source and verification.** Full text ×2 for Eqs. 2, 3, 6, 7, 19 and for q and H.

**Model.** Steady 1-D focused transport of GCR protons along a field line, of the form of Eq. (1) with 1/L = −(1/B)dB/dz. It is used to estimate anisotropies near the Sun.

- D_μμ = D_0(1−μ²)(|μ|^{q−1} + H) (Eq. 19), with q = 5/3 and H = 0.05 (citing Beeck & Wibberenz 1986). This is Eq. (26).
- D_0 follows from the prescribed λ∥ by Eq. 6 (MFP v1.2, Eq. 13).
- λ∥ = 0.1 and 1 au, energy-independent; 100 MeV protons.
- Boundaries: absorbing at z = 0.05 au; isotropic at z = 1.2 au.

**Weak-focusing results** (λ∥ ≪ L).

- Eq. 10: ∂F/∂t + (κ∥/L)∂F/∂z = κ∥∂²F/∂z².
- Eq. 17: F(z)B(z) = const.
- Eq. 18: A = 3S/(vF) = −λ∥/L(z). This diverges as L → 0, and the authors note that the weak-focusing condition is not met for all parameters.

**Solver.** The finite-difference code of van den Berg et al. (2020) [P1]. In its current public version (code, ×1), μ = ±1 are zero-flux faces and q = 1.67, H = 0.05 are hard-coded. Its comments say H = 0.005 (D-112).

Model ID GCR-PA-STRAUSS22 (a preset of PA-QFORM).

### 15.5 Kubo, Kataoka & Sato (2015)

**Source and verification.** The paper models GLE protons above 100 MeV, not GCR. Full text ×2 for Eqs. 6, 7, 19 and for q and h.

**D_μμ.** D_μμ = D_0vp^{q−2}(|μ|^{q−1} + h)(1−μ²) (Eq. 6), with q = 5/3 and h = 0.2 (citing Qin et al. 2006). Then λ∥ ∝ p^{2−q}, independent of v (MFP v1.2, Section 6.4; preset SEP-KUBO15).

Eq. 7 is Eq. (23), and Eq. 8 sets λ_rr = λ∥cos²Φ, with λ_rr "assumed to be constant across the entire interplanetary space" and scanned over 0.2–2.0 AU (0.8 AU at 1 GV in Fig. 4). With D₀ constant, λ∥ does not depend on position, so a constant λ_rr requires D₀ ∝ cos²Φ (derived). How D₀ is obtained from λ_rr is not printed (Section 19, D-22).

**Backward SDEs** (Eqs. 17–21):

$$
ds=-(\mu v+V_s)\,d\tau,\qquad dp=-\dot p\,d\tau,\qquad d\mu=-\left(\dot\mu-\partial_\mu D_{\mu\mu}\right)d\tau+\sqrt{2D_{\mu\mu}}\,dw,\qquad d\tau=-dt,
\tag{147}
$$

- dw is Gaussian with variance dτ.
- ṗ and μ̇ are given by Eqs. 15–16 for a Parker spiral.
- The μ drift has the same sign structure as a_μ in Eq. (142).
- Particles are terminated at 0.05 and 80 AU. The solution (Eq. 22) averages the boundary values.
- The μ = ±1 treatment, particle number and time step are not stated.

Model ID GCR-PA-KUBO15.

### 15.6 Heliopause and interstellar medium

**Florinski (2013).** The equations are from the single-author arXiv preprint [P46]. A four-author Letter with a shorter title appeared later (ApJL 776, L37, 2013); its abstract reports a pancake rather than a loss-cone distribution, and its equations were not read. Eqs. 1–2 and parameters: Full text ×2. The model treats gyrophase-resolved 1 MeV protons crossing from the heliosheath into the "magnetic highway":

$$
\frac{\partial f}{\partial t}+v\sqrt{1-\mu^2}\cos\phi\,\frac{\partial f}{\partial x}+v\mu\frac{\partial f}{\partial z}-\Omega\frac{\partial f}{\partial\phi}=D\frac{\partial}{\partial\mu}\!\left[(1-\mu^2)\frac{\partial f}{\partial\mu}\right]+\frac{D}{1-\mu^2}\frac{\partial^2f}{\partial\phi^2}
\tag{148}
$$

- D is a rate; λ∥ = vτ with τ = (2D)⁻¹.
- Numerics: the "random walk on a sphere technique" (Ellison et al. 1990, not read), "that avoids the singularity at μ = ±1". For small angles it is "identical to a two-dimensional random walk on a plane". The step distribution is not given.
- Boundaries: isotropic at the heliosheath side; escape at z_max = 32 AU.
- Parameters: D_HS = 10⁻³ s⁻¹; D_MH = 2 × 10⁻⁶ s⁻¹; B = 3 and 4 μG; λ∥ = 23 AU; r_g = 0.0024 AU.

Derived (fixture F-GCR-PA-05): the right-hand side is D times the Laplacian on the sphere. Its μ part is Eq. (21) with D_μμ = D(1−μ²). Eq. (24) then gives λ∥ = v/(2D), which is 23.11 AU (D_MH) and 0.0462 AU (D_HS) for 1 MeV protons. The printed r_g = 0.0024 AU corresponds to 4 μG.

Model ID GCR-PA-FLORINSKI13 (not fully specified). Florinski et al. [P132] use isotropic scattering with D_0 ≈ 10⁻⁴ s⁻¹ in the inner heliosheath and a value at least 100 times smaller in the LISM. Their Eq. 3.2 is not legible.

**Strauss & Fichtner (2014)** (×1). A 2-D model: ∂f/∂t = −vμ∂_yf + ∂_μ(D_μμ∂_μf) + ∂_x(D_⊥∂_xf), solved with backward SDEs.

- D_⊥^FLRW ∼ |μ| and D_⊥^SCAT ∼ √(1−μ²) (Eqs. 2–3).
- D_μμ is not printed; its μ-dependence "follows Dröge et al. (2010)".
- λ∥ ≈ 1 AU inside the heliopause and ≈ 1000 AU outside; v = 1 AU h⁻¹.

The model is not fully specified.

**Zhang, Zuo & Pogorelov (2014)** (×1). TeV protons are traced backward on scatter-free full orbits. D_μμ enters only the LISM distribution:

- Eq. 6: δf = −(v/2)(∂F/∂z)∫_0^μ(1−μ²)/D_μμ dμ.
- For D_μμ = D_0(1−μ²) this becomes Eq. 7, δf = −λ∥μ∂F/∂z, with λ∥ = (3v/8)∫(1−μ²)²/D_μμ dμ.
- Derived check: both expressions give λ∥ = v/(2D_0).
- Parameters: B_ISM = 3 μG, κ∥ = 10²⁸ cm² s⁻¹, κ⊥/κ∥ = 0.001.

Model ID GCR-PA-ZHANG14-LISM.

### 15.7 Benchmark: Effenberger & Litvinenko (2014)

**Source and verification.** Eq. (1) with D_μμ = D_0(1−μ²). Eqs. 1, 8, 28–30 and the boundary rule: Full text ×2.

$$
\kappa_\parallel=\lambda_0v\left(\frac{\coth\xi}{\xi}-\frac{1}{\xi^2}\right),\qquad \xi=\frac{\lambda_0}{L},\qquad \lambda_0=\frac{v}{2D_0},\qquad \kappa_{\parallel,0}=\frac{v^2}{6D_0}
\tag{149}
$$

- SDEs: dz = μv dt and dμ = [(v/2L)(1−μ²) − 2D_0μ]dt + √(2D_0(1−μ²))dW (Eq. (25)).
- The Milstein step (Eq. 30) is printed without the leading μ_t (as printed).
- "We use reflecting boundaries at µ = ±1 to conserve the probability."
- Runs: 10⁷ particles; ξ = 0 and 1.5.

Derived (fixture F-GCR-PA-03): κ∥/κ∥,0 = A_1^iso(ξ)/ξ, with A_1^iso from Eq. (145). For example, κ∥/(λ_0v) = 0.292083 at ξ = 1.5. Model ID GCR-PA-EL14-BENCH. Eq. (149) is the coefficient of the diffusion approximation; the long-time centred spatial variance of the same SDE is Eq. (115), 0.178276 λ_0v at ξ = 1.5 (Section 12.3), and a pusher's variance is compared with that value, not with Eq. (149).

### 15.8 Models found not to be pitch-angle-resolved

- **Parker-equation models:** Luo et al. [P133]; Shen & Qin [P134]; Kóta & Jokipii [P135] (abstract); the NWU and HelMod families.
- **Full-orbit Forbush-decrease models without scattering:** [P136, P137, P138]. Lara et al. add white noise to B.
- **Sectored heliosheath field** with the pitch angle fixed at π/2 [P139].
- **SDE review:** [P44] contains no pitch-angle equations.
- **Read only as abstracts:** [P140, P141, P142] describe pitch-angle-dependent GCR transport, but no D_μμ is recorded.

### 15.9 Implementation requirements

**Implementation requirement.** A GCR pitch-angle preset stores:

- the D_μμ shape and its parameter (ε, h or H);
- the quantity that fixes its amplitude (a κ∥(P) model, turbulence inputs, or λ∥);
- κ⊥ and D_⊥(μ).

In vandenBerg2026, h(μ) is tied to the D_μμ shape.

**Implementation requirement.** An amplitude fixed by κ∥ or λ∥ is computed exactly from Eq. (23), using Eq. (140), MFP v1.2 Eq. 13, or quadrature. It is never computed from the asymptote of φ.

**Implementation requirement.** A backward-SDE kernel:

- uses its printed sign convention (Eqs. (142) and (147));
- records its μ boundary rule; where the source states none (QinWu2020, Kubo2015), the rule is a configuration item;
- if it uses Eq. (144), which vanishes at |μ| = 1 for D^iso, states a floor on Δτ or shows that |M| = 1 is not reached.

**Implementation requirement.** The printed μ^{s−1} of QinWu2020 runs only as GCR-PA-QINWU20-PRINTED.

**Implementation requirement.** Eq. (145) and Eq. (149) are reference results for validating a focused-transport SDE with constant L′: Eq. (145) for the forward stationary pitch-angle distribution and its first anisotropy, Eq. (149) for the mean streaming speed u = κ∥/L; the centred spatial variance is validated against Eq. (115) (VAL-PUSH-KVAR), and a backward-SDE histogram against ρ₋ of Section 12.3, not against Eq. (145).

## 16. Pitch-angle scattering in a particle pusher: numerical formulation

This section states how the operator Eq. (21) is carried into a pseudo-particle code: the stochastic differential equation (SDE) and its sign rules in forward and backward time, importance sampling and weights, integration schemes, the boundary at μ = ±1, the point μ = 0, time-step statements, kernels that are not Itô SDEs, frames, validation, and printed items that must not be copied. It extends MFP v1.2, Section 21.5 (Eq. 46 there is the pure-scattering case of Eq. (25)) and uses the shapes of Section 5 and MFP v1.2, Sections 6.2–6.6. Every item marked derived was checked numerically; the scripts, outputs and seeds are in derivation record DR-16 (fixtures F-PUSH-01 to F-PUSH-18, Monte Carlo runs MCR-1 to MCR-7). Monte Carlo (MC) runs use numpy 2.5.3 (PCG64 generator, fixed seeds); exact numbers use mpmath at 40 digits.

Notation. A(μ) is the total Itô drift, b(μ) the noise amplitude, ζ a standard normal variate (ξ is reserved for the focusing strength), Δt (Δτ, ds) the forward (backward) step, ν the frequency of Eq. (24); D′_μμ = ∂D_μμ/∂μ; "mirror rule" is Eq. (161).

### 16.1 The Itô SDE for the pitch-angle operator

**Forward form and its Kolmogorov equation.** For the Itô SDE dx = a dt + b dW, Strauss & Effenberger [P44] (Eqs. 2, 8, 9; Full text ×2 for Eqs. 2, 8, 9) print the time-forward Kolmogorov (Fokker–Planck) equation ∂ρ/∂t = −Σ∂_i(ã_iρ) + ½ΣΣ∂_i∂_j(C̃_ijρ) with C = b·b, and give b = √(2κ) for a diffusion coefficient κ (text after their Eq. 48). With A = dμ/dt + D′_μμ and b = √(2D_μμ), Eq. (25) therefore transports the density P(μ, t) of pseudo-particles by (derived)

$$
\frac{\partial P}{\partial t}=-\frac{\partial}{\partial\mu}\Big(\frac{d\mu}{dt}\,P\Big)+\frac{\partial}{\partial\mu}\Big(D_{\mu\mu}\frac{\partial P}{\partial\mu}\Big).
\tag{150}
$$

The pitch-angle term is Eq. (21) exactly; the deterministic term is in conservative (flux) form. A forward SDE therefore samples the density whose transport equation has ∂/∂μ[(dμ/dt)P], not (dμ/dt)∂f/∂μ. The sources handle this as follows. Effenberger & Litvinenko [P2] (Full text ×2) write the focused transport equation for f₀ with (v/2L)(1−μ²)∂f₀/∂μ (Eq. (1)) and integrate the SDE for the linear density f = e^{z/L}f₀ (their Eq. 3), "which is completely equivalent to the Fokker-Planck equation for the linear density". Qin et al. [P5] simulate g = A(z)f with A the flux-tube area (Eq. 7; equations partly not legible, structure only). Merten & Aerdker [P9] (Eq. 1) and van den Berg et al. [P1] (Eq. 9) print the conservative form directly. Derived check (F-PUSH-14): for constant L, F = e^{z/L}f₀ satisfies ∂_t F + ∂_z(μvF) + ∂_μ[(1−μ²)vF/(2L)] = ∂_μ(D_μμ∂_μF) exactly when f₀ satisfies Eq. (1) (residual 0 at 40 digits).

Forward Itô pitch-angle SDEs as printed:

| Source | Equation | Drift | Noise | Verification |
|---|---|---|---|---|
| [P2] | Eq. 28 | (v/2L)(1−μ²) − 2D₀μ | √(2D₀(1−μ²)) | Full text ×2 |
| [P7] | Eq. 4.1 | υ/(2L)(1−μ²) + ∂D_μμ/∂μ | √(2D_μμ) (radical lost) | Full text ×1 |
| [P143] | Eq. 3.2 | (v/2L)(1−μ²) + ∂D_μμ/∂μ | √(2D_μμ) (radical extent lost) | Full text ×1 |
| [P25] | Eq. 13, stated and not used | dμ/dt + ∂D_μμ/∂μ | √(2D_μμ) | Full text ×2 |
| [P9] | Eq. 19 | (v/2L)(1−μ²) + ∂D_μμ/∂μ | √(2D_μμ) | Full text ×1 |
| [P83] | Eq. 32, wave rest frame, scattering only | a = ∂_μD_μμ | b² = 2D_μμ | Full text ×1 |
| [P5] | Eqs. 9–10 | v(1−μ²)/(2L) + ∂D/∂μ | √(2D) | structure only (garbled) |
| [P8] | Eq. A.3 | ½(v/L + μR)(1−μ²) + ∂D_μμ/∂μ | √(2D_μμ dt) W | not legible as printed; structure only |
| [P10] | Eq. 4 | focusing and flow terms, no ∂D_μμ/∂μ | √(2D_μμ) dW_μ(dt), dW = dt·Λ | as printed; see 16.10 |

Effenberger & Litvinenko ("the Ito calculus") and Kartavykh et al. name the Itô calculus and Merten & Aerdker "Itô's Lemma"; the others do not name it; no source read uses a Stratonovich form.

**Operator conventions.** If a code writes the scattering term with a different factor, the Itô drift and noise follow the coefficient of ∂_μ(·∂_μf) (derived from Eq. (150); MFP v1.2, Section 21.5):

$$
\frac{\partial}{\partial\mu}\Big(\frac{\varphi}{2}\frac{\partial F}{\partial\mu}\Big):\; A_{\rm scat}=\tfrac12\varphi',\; b=\sqrt{\varphi};\qquad
\frac{\partial}{\partial\mu}\Big(\frac{D^{\rm E}_{\mu\mu}}{2}\frac{\partial f}{\partial\mu}\Big):\; A_{\rm scat}=\tfrac12 D^{\rm E\,\prime}_{\mu\mu},\; b=\sqrt{D^{\rm E}_{\mu\mu}} .
\tag{151}
$$

The first is Ruffolo's scattering coefficient φ(μ) [P6] (term ∂/∂μ (φ(μ)/2) ∂/∂μ read verbatim twice in the PDF; φ = 2D_μμ); the second is the EPREM form [P17, P32] (Young et al. Eq. A1, two renderings; Linker et al. one rendering; partly legible in Schwadron et al. [P34]). EPREM is a grid code; its factor-of-two question (MFP v1.2 Section 6.6 and MFP v1.2 D-2) is unchanged here.

**Time-backward form.** Zhang et al. [P15] (Eq. 6; ar5iv only, two fetches agree) print

$$
d\mu(s)=\Big[-\frac{d\mu}{dt}+\frac{\partial D_{\mu\mu}}{\partial\mu}\Big]ds+\sqrt{2D_{\mu\mu}}\,dw(s),
\tag{152}
$$

with s the backward running time and dw(s) "a Wiener process as a function of backward running time s"; the spatial and momentum equations (Eqs. 5, 7) likewise carry −vμb̂_o − V − V_d and −dp/dt, while the diffusive drift ∇·κ⊥ keeps its sign. The same sign rule is printed by Kubo et al. [P16] (Eqs. 17–21; Eq. 19 in two renderings: dμ = −(μ̇ − ∂_μD_μμ)dτ + √(2D_μμ)dw, dτ = −dt, p(dw; dτ) = (2πdτ)^{−1/2}exp(−dw²/2dτ)), by van den Berg et al. [P14] (Eq. 38e a_μ = ∂D_μμ/∂μ − dμ/dt and B_μμ = √(2D_μμ); PDF text layer, one rendering, two passes agree) and by He & Wan [P20] (unnumbered; see 16.2). The rule is: deterministic rates change sign, ∂D_μμ/∂μ and the noise do not. Strauss & Effenberger print the backward Kolmogorov equation (Eq. 7). Derived check: for the advective form ∂_t f + (dμ/dt)∂_μ f − ∂_μ(D_μμ∂_μ f) = Q, the μ-part −(dμ/dt)∂_μ f + D′∂_μ f + D∂²_μ f is the generator a∂_μ + ½b²∂²_μ of Eq. (152) (product rule, F-PUSH-18e), so the backward SDE samples f itself, not a flux-tube density. MCR-7 confirms this for dμ/dt = 2(1 − μ²), D_μμ = ½(1 − μ²), f(0, μ) = μ, t = 0.5: the backward SDE gives f(t, 0.3) = −0.4327 ± 0.0004 against −0.4324 from a finite-difference solution of the advective form, while the same SDE with the forward sign of dμ/dt gives 0.661. Zhang et al. (Eq. 8) evaluate f as the mean of the source integral along the trajectory plus the boundary or initial value at the first exit point; van den Berg et al. start pseudo-particles with "pitch cosine M(0) = µf (uniformly distributed on [−1; 1]), and weight w(0) = wf = 1" and state that "The time-backward approach described here is mathematically equivalent to a time-forward approach in the stationary limit".

μ reference direction: Zhang et al. measure μ "to the outward magnetic field line" b̂_o, so that f stays continuous across polarity reversals; an implementation coupling to such a code must convert μ when B = −B b̂_o.

### 16.2 Importance sampling and weights

**Zhang et al. (2023).** Zhang et al. [P15] (ar5iv, one rendering) substitute f = (1 + μ/a)u with a > 1 "constant tuning parameter" (text before their Eq. 9, which is the equation for u). The transformed equation has the same form for u with a modified rate and a loss term:

$$
\frac{d\mu'}{dt}=\frac{d\mu}{dt}-\frac{2D_{\mu\mu}}{a+\mu},\qquad
c=\frac{1}{a+\mu}\Big(\frac{d\mu}{dt}-\frac{\partial D_{\mu\mu}}{\partial\mu}\Big),
\tag{153}
$$

(Eqs. 10–11), and the solution (Eq. 12)

$$
f(t,\mathbf x,\mu,p)=\Big(1+\frac{\mu}{a}\Big)\Big\langle\int_0^t\frac{Q_0(t-s,\mathbf x(s),p(s),\mu'(s))}{1+\mu'(s)/a}\,
\exp\!\Big(-\!\int_0^s c\,(t-s_1,\mathbf x(s_1),p(s_1),\mu'(s_1))\,ds_1\Big)ds\Big\rangle .
\tag{154}
$$

Their Eq. 13 is the backward SDE for μ′, dμ′(s) = [−dμ′/dt + ∂D_μμ/∂μ]ds + √(2D_μμ)dw(s) (ar5iv, two of three readings; one reading did not find it). It is the form that Eq. (152) gives with dμ′/dt. Derived check (F-PUSH-13): inserting f = (1 + μ/a)u in the μ-part of the transport equation and dividing by 1 + μ/a gives (dμ′/dt)∂_μu − ∂_μ(D_μμ∂_μu) + cu with Eq. (153) exactly (residual ≤ 5 × 10⁻⁴¹ at three μ for a = 1.7 and arbitrary smooth D_μμ, dμ/dt, u). In backward time the extra drift is +2D_μμ/(a + μ) > 0 (derived from the sign rule): it pushes μ toward +1, and through the −vμb̂_o term of their Eq. 5 this moves backward trajectories toward the Sun. Strauss & Effenberger describe the method of Qin et al. (2005) and Zhang et al. [P40] as adding "an artificial convective velocity ... to advect the pseudo-particles towards the inner boundary" (Secondary, via [P44]; Zhang 2009: Abstract only). Zhang et al. do not list the value of a ("we do not explicitly list the value of a used in each of our calculations"); they also split trajectories in two "after a certain amount of time", halving the weight, recursively (split time not given).

**He & Wan (2015) and successors.** He & Wan [P20] (ar5iv, one rendering), He et al. [P22] and He & Wan [P21] print the backward pitch-angle SDE

$$
d\mu=\sqrt{2D_{\mu\mu}}\,dW_\mu(s)-\frac{1-\mu^2}{2}\Big[\frac{V}{L}+\mu\Big(\frac{\partial V^{sw}_x}{\partial x}+\frac{\partial V^{sw}_y}{\partial y}-2\frac{\partial V^{sw}_z}{\partial z}\Big)\Big]ds+\Big(\frac{\partial D_{\mu\mu}}{\partial\mu}+\frac{2D_{\mu\mu}}{M+\mu}\Big)ds
\tag{155}
$$

(as printed; V "pseudo-velocity", L = (ẑ·∇ln B)⁻¹). M is not defined on the page and no weight or loss factor is printed. The term has the form of the backward drift of Eq. (153) with M in place of a (derived comparison), but without c and without the factors of Eq. (154) the printed system does not sample f. He et al. (2017) state that the SDEs are "numerically solved with an Euler scheme as usual" and that trajectories "exit the system when the trajectories hit the physical boundaries for the first time".

**van den Berg et al. (2026).** The weight obeys (Eq. 18, PDF text layer, one rendering)

$$
w(\tau_k+\Delta\tau_k)\approx w(\tau_k)\,e^{\mathcal L(\mathbf Q(\tau_k),\tau_k)\,\Delta\tau_k},\qquad
\mathcal L=-\frac{\mu v}{c^2}\,\hat{\mathbf b}\cdot\Big[\frac{\partial\mathbf u}{\partial t}+(\mathbf u\cdot\nabla)\mathbf u\Big]
\tag{156}
$$

(𝓛, printed Ł, from Eq. 37), the loss term of their backward transport equation (Eq. 15); it is not an importance-sampling weight. Its sign relative to Eqs. 15 and 37 is not established: the Feynman–Kac weight for the advective form of Eq. (16) is exp(−∫𝓛dτ) (derived; D-125). The spectrum at the observer is "the weighted average over all pseudo-particles" (LIS weighting form not read).

### 16.3 Integration schemes

**Euler–Maruyama.** Strauss & Effenberger (Eqs. 4–5): x(t + Δt) = x(t) + aΔt + bΔW with ΔW = √Δt Λ, Λ standard normal (summarized reading). In μ:

$$
\mu_{n+1}=\mu_n+A(\mu_n)\,\Delta t+\sqrt{2D_{\mu\mu}(\mu_n)\,\Delta t}\;\zeta_n .
\tag{157}
$$

Merten & Aerdker print this as Eq. 24 (μ_{n+1} = μ_n + A_μh + B_μη_μ√h), with s_{n+1} = s_n + vμ_nh (Eq. 25) and B_μ² = 2D_μμ (Eq. B.34); van den Berg et al. use it in backward time (Eq. 17, "Λ(τk) are n independent normally distributed pseudo-random numbers"); Wijsen et al. use it for the spatial SDE (Eq. 12) and the deterministic Euler method for dμ/dt and dp/dt. One-step moments (derived, F-PUSH-04a/04b): E[Δμ] = AΔt and E[Δμ²] = A²Δt² + 2D_μμΔt.

**Milstein.** Effenberger & Litvinenko print (Eqs. 29–30, dimensionless: distance in λ₀, speed in v, time in λ₀/v)

$$
z_{t+\Delta t}=z_t+\mu_t\Delta t,\qquad
\mu_{t+\Delta t}=\Big[\tfrac12\xi(1-\mu_t^2)-\mu_t\Big]\Delta t+\sqrt{\Delta t(1-\mu_t^2)}\,\epsilon_t-\tfrac12\mu_t\Delta t(\epsilon_t^2-1)
\tag{158}
$$

as printed in both renderings: the right-hand side lacks the leading μ_t (Section 16.10). For b = √(2D_μμ), b b′ = D′_μμ wherever D_μμ > 0 (derived; F-PUSH-03 gives 0.607355341732574536 for both sides at μ = 0.3, q = 5/3, H = 0.05, D₀ = 1), so the Milstein step for Eq. (25) is (derived)

$$
\mu_{n+1}=\mu_n+\Big(\frac{d\mu}{dt}+D'_{\mu\mu}\Big)\Delta t+\sqrt{2D_{\mu\mu}\Delta t}\;\zeta_n+\tfrac12 D'_{\mu\mu}\,\Delta t\,(\zeta_n^2-1),
\tag{159}
$$

all coefficients at μ_n. With λ₀ = v/(2D₀) = 1 and v = 1 (so D₀ = ½ and v/L = ξ), Eq. (159) for D_μμ = D₀(1 − μ²) reproduces every term of Eq. (158) and adds μ_t (derived). One-step moments (derived, F-PUSH-04a/04b): E[Δμ] = AΔt, E[Δμ²] = A²Δt² + 2D_μμΔt + ½D′²Δt². The Milstein correction contains D′_μμ and is therefore singular at μ = 0 for the q-form (Section 16.5).

**Predictor–corrector of Sandroos et al. (2014).** As reported from the PDF (one rendering, two fetches; wave-frame tilde dropped): trial value μ_tr = μ_t + a(μ_t)dt + b(μ_t)√dt N (Eq. 33), a₀ = (1 − η)a(μ_t), a₁ = (1 − η)a(μ_tr) (Eqs. 34–35), and

$$
\mu_{t+1}=\mu_t+\big(\theta a_1+(1-\theta)a_0\big)dt+\big(\eta\,b(\mu_t)+(1-\eta)\,b(\mu_{\rm tr})\big)\sqrt{dt}\,N,\qquad \theta=0,\ \eta=0.5
\tag{160}
$$

(Eq. 36; the first term on the right is printed μ̃₀, read here as μ̃_t because the text calls μ̃_t and μ̃_{t+1} the initial and new pitch-angle cosines, D-120), with a = ∂_μD_μμ, b² = 2D_μμ (Eq. 32). Derived remark: the corrected drift (1 − η)a equals a − ηbb′ because bb′ = a for pure scattering; with η = ½ the noise is the average of the start and trial values. In MCR-1 (Section 16.6) this scheme converges to the Itô moments.

**Other statements.** Kartavykh et al. [P7] describe the Dröge-group method as "a time-forward stochastic process with" "an implicit scheme" with Lorentz transformations between frames (scheme not specified in the text read); Zuo et al. [P31] mention Runge–Kutta (sentence truncated). Van den Berg's thesis [P144] (Abstract only) reports that "different numerical integration schemes, of increasing accuracy, do not yield any improvement" for the SEP problem studied; its Appendix B (Euler–Maruyama, Milstein, order-1.5 scheme, variable step) was not read. Wijsen's thesis [P145] (Chapter 1 only) and van den Berg et al. (2026) call Euler–Maruyama schemes "unconditionally stable" for sufficiently regular (bounded) coefficients.

**Operator splitting and order of sub-steps (as stated).**

| Source | Order within a step |
|---|---|
| [P6] (finite differences) | half μ-step, p-step, z-shift, half μ-step; Crank–Nicolson in μ |
| [P55] | radial step δr with focusing-updated μ and elapsed time; transform (v, μ) to the local solar-wind frame; scatter; transform back; record at observer (tool summary of Sect. 3.6) |
| [P25] | Euler for dμ/dt and dp/dt; Euler–Maruyama for x; scattering by the Àgueda kernel (Section 16.7) |
| [P88] (pre-2006) | move x ← x + vμδt; scatter in both wave frames; update wave energy (reader's ordering) |
| [P83] | time step; waves to mesh; scatter in wave frames; apply wave growth; move particles, then waves |
| EPREM [P34, P32] | grid code with operator splitting and "a large number of sub-cycles (for time-step stability)" |
| [P15], [P20], [P16], [P9] | not stated (Merten: one Euler–Maruyama step for (s, μ)) |

No source read compares splitting orders or quantifies the splitting error.

### 16.4 The boundary at μ = ±1

**Does the exact process reach μ = ±1?** For Eq. (25) with pure scattering, the scale density of the diffusion is exp(−∫A/D dμ) = 1/D_μμ (derived). If D_μμ vanishes linearly at μ = ±1, as the shapes with the factor (1 − μ²) do, ∫dμ/D_μμ diverges logarithmically at the boundary (F-PUSH-18d: ∫_{0.5}^{1−δ}dμ/(1 − μ²) = 4.40, 9.01, 13.61 for δ = 10⁻⁴, 10⁻⁸, 10⁻¹²), so the exact process does not reach μ = ±1 (derived; Feller's test). This answers, for the exact process, Merten & Aerdker's footnote 9 on whether boundaries are needed as h → 0. The discrete step does leave [−1, 1]: for isotropic scattering, an Euler–Maruyama step from 1 − μ ≈ νΔt exits when ζ exceeds a threshold that tends to √2 (F-PUSH-10: 1.3784, 1.4107, 1.41418 for νΔt = 0.1, 0.01, 10⁻⁴), the largest exit probability per step is therefore 1 − Φ(√(2 − νΔt)) (0.0840392 at νΔt = 0.1, F-PUSH-10), which tends to 1 − Φ(√2) = 0.0786496 as Δt → 0, so a fixed fraction of the particles in that layer leaves per step at any step size (derived). The fraction of all particles reflected per step is proportional to Δt (MCR-1: 2.7 × 10⁻² at νΔt = 0.1, 2.6 × 10⁻³ at 0.01 for Euler–Maruyama).

**Reflection rules as printed.**

$$
x<x_{\min}:\ x\to-2x_{\min}-x=-2-x,\qquad x>x_{\max}:\ x\to 2x_{\max}-x=2-x
\tag{161}
$$

(Strauss & Effenberger Eqs. 25–26 as printed; with x_min = −1 the middle expression −2x_min − x of Eq. 25 equals 2 − x, whereas the reflection about x_min is 2x_min − x = −2 − x, the right-hand expression (derived; D-121); generic von Neumann boundary ∂f/∂x = 0 at x = ±1, "similarly to hard sphere reflection off a planar surface"; the review does not apply it to μ). Merten & Aerdker (Appendix B.3.1, read in two fetches) print, with m = μ_n + Δμ,

$$
\mu_{n+1}=m\ \ (|m|\le1),\qquad \mu_{n+1}=\operatorname{sign}(m)-\operatorname{sign}(m)\,\big|m-\operatorname{sign}(m)\big|\ \ (|m|>1),
\tag{162}
$$

with the example μ = 0.9, Δμ = 0.15 → 0.95, and reject periodic boundaries (−0.95) as unphysical; "reflective boundaries are applied during the propagation step". Eq. (162) equals Eq. (161) (derived: for m > 1 it gives 2 − m). Van den Berg et al. (Eq. 23, text layer "M ← sign(M) 2 − M if |M| > 1") is read as M ← sign(M)(2 − |M|), which is the same rule; brackets and inner bars are not visible, so the reading is unverified. Their inner spatial boundary is reflecting with M ← |M| (Eq. 20). Effenberger & Litvinenko: "We use reflecting boundaries at µ=±1 to conserve the probability." Sandroos et al. are reported to reflect by "μ̃ → 1 − μ̃" (Section 16.10). Zhang et al., the Qin/He papers, Kubo et al., Laitinen et al., Kartavykh et al., Tenishev et al. and the EPREM papers state no rule.

A single application of Eq. (161) returns a value in [−1, 1] only if |m| ≤ 3. With the q-form drift at small |μ| an Euler–Maruyama step can exceed this (MCR-2 recorded 15 such values at Δt = 0.1/D₀ and 2 at 0.03/D₀ for q = 5/3, H = 0.05); repeated application (folding) is then needed (derived).

**Grid codes.** Ruffolo (1995) requires the μ-flux to vanish as μ → ±1; Strauss & Fichtner [P11] (Eq. A11): "the flux through cell face i = N + 1/2 (located at μ = 1) is always zero"; the SEP_propagator code of the primer [P146] places μ = ±1 on cell faces with zero diffusive flux there.

**Kernels without reflection.** For the random rotation of Eq. (168), |μ′| ≤ 1 for every ϑ and φ (derived by the Cauchy–Schwarz inequality; μ′ is the field-aligned component of a rotated unit vector), so no boundary rule is needed. Florinski [P46] uses "the random walk on a sphere technique" (Ellison et al. 1990) "that avoids the singularity at μ = ±1". Van den Berg et al. shrink the μ step near the boundary through l_μ = ℓ√(1 − μ²) (Eq. (165)).

### 16.5 The point μ = 0

**Drift of the regularized shapes (derived).** For the q-form Eq. (26), D_μμ = D₀(1 − μ²)(|μ|^{q−1} + H),

$$
\frac{\partial D_{\mu\mu}}{\partial\mu}=D_0\Big[-2\mu\big(|\mu|^{q-1}+H\big)+(q-1)(1-\mu^2)|\mu|^{q-2}\operatorname{sgn}\mu\Big],
\tag{163}
$$

whose second term diverges at μ = 0 for 1 < q < 2 and every H ≥ 0 (H enters only the regular first term). The singular term is odd, so the drift points away from μ = 0 on both sides; it is integrable. For H = 0 also D_μμ(0) = 0, so the noise vanishes where the drift diverges. For the ε-form Eq. (29), D_μμ = D₀(1 − μ²)(|μ|/(1 + |μ|) + ε),

$$
\frac{\partial D_{\mu\mu}}{\partial\mu}=D_0\Big[-2\mu\Big(\frac{|\mu|}{1+|\mu|}+\varepsilon\Big)+\frac{(1-\mu^2)\operatorname{sgn}\mu}{(1+|\mu|)^2}\Big],
\tag{164}
$$

which is bounded and jumps from −D₀ to +D₀ at μ = 0 for every ε. Values (F-PUSH-01, -02; D₀ = 1; checked against numerical derivatives):

| μ | q = 5/3, H = 0.05: D′ | its singular term | q = 5/3, H = 0: D′ | ε = 0.048: D′ |
|---|---|---|---|---|
| 0.1 | 1.36883820162041 | 1.42192689542104 | — | — |
| 0.01 | 3.09215479871956 | 3.09408311648628 | 3.09315479871956 | — |
| 0.001 | 6.66654 | 6.66666 | 6.66664 | — |
| 10⁻¹² | — | — | — | 0.999999999997904 |
| 0.3 | — | — | — | 0.3712 |

The same holds for the M-FLAMPA shape |μ|^{2/3} (MFP v1.2, Sections 6.7, 21.5).

**Crossing of μ = 0 (derived).** With scale density 1/D_μμ, the two hemispheres communicate under the exact process if and only if ∫dμ/D_μμ is finite across μ = 0. If D_μμ ∝ |μ|^α near 0, this holds for α < 1 and fails for α ≥ 1 (F-PUSH-16: for α = 2/3 the integral from δ to 0.5 tends to 2.48001376922875; for α = 1 it grows as ln(1/δ)). Hence the strict q-form (H = 0, 1 < q < 2) still lets particles cross μ = 0, while Àgueda's ν(μ) = ν₀|μ|/(1 + |μ|) with ε = 0 does not; the thesis [P55] states that the linear dependence "guarantees that particles are not able to cross μ = 0 by scattering".

**Step-size consequence (derived).** Krülls & Achterberg's condition ⟨Δx⟩² ≪ ⟨Δx²⟩ (quoted by Strauss & Effenberger) applied to μ gives Δt ≪ 2D_μμ/D′²_μμ. For the q-form this bound tends to zero as μ → 0 for every H (≈ 2H|μ|^{2(2−q)}/((q − 1)²D₀) for H > 0); for the ε-form it tends to 2ε/D₀ (F-PUSH-11, D₀ = 1):

| μ | q = 5/3, H = 0.05 | q = 5/3, H = 0 | ε = 0.048 |
|---|---|---|---|
| 0.1 | 0.2805 | 0.2244 | 0.4403 |
| 0.01 | 0.02017 | 0.009702 | 0.1208 |
| 0.001 | 0.002700 | 0.000450 | 0.09841 |
| 0.0001 | 0.000506 | 0.0000209 | 0.09624 |

**What the sources do.**

- Regularization parameters H, h, ε: values are listed in Section 5 (for example h = 0.2 [P15, P16, P19], h = 0.01 [P10], H = 0.05 [P11], H = 0.5 [P7], ε = 0.048 [P25], ε = 0.045 [P14]). None of them removes the singularity of Eq. (163). Zhang et al.: "The result of our calculation is not very sensitive to h0 unless h0 ≪ 0.05".
- Wijsen et al. [P25] do not integrate Eq. 13 because it "requires the use of a very small time step to obtain correct PADs for values of μ near zero", and use the kernel of Section 16.7 with two isotropic processes.
- Van den Berg et al. chose the ε-form "because its derivative is well behaved around µ = 0" (bounded, though discontinuous, by Eq. (164)).
- Afanasiev et al. [P78] fill the resonance gap |μ| < μ_min = Ω/(vk_max) with ν̃(μ) = bμ + b₁, printed linear in μ (Sections 13.7, 14.6); separately they state that the method "requires extremely small time steps" near μ = 0.
- Sandroos et al. use μ̃_min = 0.025 for the resonant wavelength near μ̃ = 0.
- Ruffolo (1995) uses an effective coefficient φ_eff after Ng & Wong (1979) (form not read).
- Chen et al. [P47] cut the κ∥ integral at μ_min = 0.05; this is a normalization device, not a pusher rule.

**Numerical demonstration (derived, MCR-2).** Euler–Maruyama with the mirror rule from a uniform start, D₀ = 1, t = 1, N = 2 × 10⁵: the exact distribution stays uniform, so the ratios below should be 1 (standard errors 0.010 for the first column and 0.022 for the second).

| D₀Δt | q = 5/3, H = 0.05 | q = 5/3, H = 0 | ε = 0.048 |
|---|---|---|---|
| 0.1 | 0.830 / 0.852 / 0.730 | 0.762 / 0.737 / 0.727 | 0.898 / 0.873 / 0.687 |
| 0.03 | 0.816 / 0.823 / 0.723 | 0.727 / 0.674 / 0.723 | 0.943 / 0.919 / 0.821 |
| 0.01 | 0.855 / 0.794 / 0.877 | 0.789 / 0.689 / 0.883 | 0.970 / 0.920 / 0.932 |
| 0.003 | 0.946 / 0.903 / 0.960 | 0.882 / 0.689 / 0.956 | 1.002 / 0.989 / 0.989 |
| 0.001 | 0.984 / 0.922 / 0.987 | 0.951 / 0.791 / 0.994 | 0.995 / 0.970 / 0.999 |
| 0.0003 | 0.993 / 0.956 / 0.987 | 0.957 / 0.910 / 1.001 | 0.995 / 0.988 / 0.997 |
| 0.0001 | 0.996 / 0.961 / 1.008 | 1.013 / 0.943 / 1.006 | 0.997 / 1.034 / 1.002 |

Entries are P(|μ| < 0.05)/0.05, P(|μ| < 0.01)/0.01 and P(|μ| > 0.95)/0.05. Euler–Maruyama depletes the region around μ = 0, where the drift of Eq. (163) pushes particles away, and the layer next to μ = ±1. At D₀Δt = 0.01 the deficit in |μ| < 0.01 is 21% (H = 0.05), 31% (H = 0) and 8% (ε-form); the ε-form is within two standard errors in all three columns from D₀Δt = 0.003, while the q-form with H = 0 is still 6% (2.6 standard errors) low in |μ| < 0.01 at 10⁻⁴. This quantifies, for these shapes, the behaviour reported by Wijsen et al.; it is not a statement about any published code.

### 16.6 Time-step statements

| Source | Statement as printed | Applies to | Verification |
|---|---|---|---|
| [P44] Eq. 51 | Δt(x) = min{l/\|dκ/dx\|, l²/b²(x)}, l "shorter than the shortest spatial structure"; example l = 0.1L; Krülls & Achterberg ⟨Δx⟩² ≪ ⟨Δx²⟩ | spatial diffusion | Full text ×1 |
| [P14] Eqs. 19, 24, 25e | Δτ_k = min_i{l_i/\|a_i\|; l_i²/Σ_jB_ij²}; μ entries l_μ/\|a_μ\| and l_μ²/B²_μμ; l_μ = ℓ√(1 − µ²), ℓ = 2⁻⁶ | backward μ-SDE | one rendering (text layer); Eq. 19 and the μ entries of Eq. 24 two passes, the rest of Eq. 24 one |
| [P83] Eq. 39 | dt_max = min(Δμ_max/a_max, Δμ_max²/b_max²); N_sub = dt_sim/dt_max, target N_sub ≤ 10 | μ-SDE in wave frame | one rendering; Δμ_max value not found |
| [P146] (code) | Δt = min(0.9Δz/v, 0.9Δμ/B_max, 0.9Δμ²/(2D_μμ,max)), B_max = v/(2L_max), Δμ = 2/49 | explicit FD grid | code read once |
| [P55] Eqs. 3.27–3.31 | δμ ≤ u_c/v (3.27); vδt/(4L) ≤ u_c/(2v) (3.28); radial step δr = min[(u_c/v)max(r₀, 0.3), 0.1λ_r](\|μ₀\| + u_c/2v)cos ψ (3.31) | focusing/frame step | one rendering; the relation symbol of the time-step condition for Eq. 3.17 ("δt [?] ν⁻¹") is not legible |
| [P25] | δt ≪ ν⁻¹ (ar5iv alttext) | random rotation | Full text ×1 for this item |
| [P88] (pre-2006) | ⟨ϑ²_±⟩ = 2ν_±δt < 0.2 for all particles | random rotation | Full text ×2 |
| [P91] (pre-2006) | "Δt = a ν_e(γ)" as extracted, a ≲ 0.01 | random rotation | one rendering; see 16.10 |
| EPREM [P34, P17] | macro steps "tens of seconds", sub-steps "as small as 1 s"; tDel = 0.01 day | grid code | as printed |
| [P2], [P15], [P20], [P16], [P7], [P9] (pitch-angle tests) | Δt not stated | — | — |

In μ, Eqs. 19, 24 and 25e of van den Berg et al. read (derived substitution with a_μ and B_μμ of Section 16.1)

$$
\Delta\tau_\mu=\min\Big\{\frac{\ell\sqrt{1-\mu^2}}{\big|\partial D_{\mu\mu}/\partial\mu-d\mu/dt\big|}\,;\ \frac{\ell^2(1-\mu^2)}{2D_{\mu\mu}}\Big\},\qquad \ell=2^{-6},
\tag{165}
$$

and Sandroos et al. print

$$
dt_{\max}=\min\Big(\frac{\Delta\mu_{\max}}{a_{\max}},\frac{\Delta\mu_{\max}^2}{b_{\max}^2}\Big),\qquad a=\partial_\mu D_{\mu\mu},\ b^2=2D_{\mu\mu}.
\tag{166}
$$

These two are the only step rules for a pitch-angle Itô SDE found that are written through D_μμ and its derivative. Neither reading reports a test of the rule near μ = 0; for the q-form the first entry of Eq. (165) tends to zero as |μ|^{2−q} (derived from Eq. (163)), while van den Berg et al. use the ε-form, for which it stays finite. At μ = ±1 exactly the first entry of Eq. (165) is zero for every shape with D′_μμ(±1) ≠ 0 (the ε-form gives D′_μμ(1) = −D₀(1 + 2ε); the second entry stays finite, ℓ²/(2D₀(½ + ε)) for the ε-form), so a beam injected at μ = ±1 receives Δτ = 0 and does not advance (derived, F-PUSH-19). Neither source prints a lower bound on Δτ or an entrance treatment for such an injection (U-55); a derivative-based step rule is therefore not a complete solver for every coefficient and initial condition, and a step floor imposed to avoid the stall must not defeat the resolution required near μ = 0 (Section 16.5). When the step is adaptive, the trajectory is sampled at common physical times or with residence-time weights: counting the end points of adaptive steps with equal weight over-represents the regions where the steps are short and produces an artificial stationary pitch-angle distribution (derived: the occupation measure of a Markov process is the time spent, not the number of steps).

**Measured Euler–Maruyama error for isotropic scattering (derived, MCR-1).** D_μμ = (ν/2)(1 − μ²), μ₀ = 0.5, νt = 1, N = 10⁶, mirror rule for the SDE schemes; exact ⟨μ⟩ = 0.183940, ⟨μ²⟩ = 0.329184; standard errors 0.00054 and 0.00030. Errors (measured − exact):

| νΔt | Euler–Maruyama ⟨μ⟩ / ⟨μ²⟩ (reflected per step) | Milstein ⟨μ⟩ / ⟨μ²⟩ | Sandroos PC ⟨μ⟩ / ⟨μ²⟩ | rotation kernel ⟨μ⟩ / ⟨μ²⟩ |
|---|---|---|---|---|
| 0.1 | −0.0289 / −0.0189 (2.7 × 10⁻²) | −0.0106 / +0.0161 | −0.0055 / −0.0111 | −0.0040 / −0.0004 |
| 0.05 | −0.0155 / −0.0113 (1.4 × 10⁻²) | −0.0041 / +0.0087 | −0.0039 / −0.0073 | −0.0013 / +0.0002 |
| 0.02 | −0.0062 / −0.0053 (5.2 × 10⁻³) | −0.0015 / +0.0038 | −0.0024 / −0.0047 | −0.0015 / +0.0001 |
| 0.01 | −0.0027 / −0.0027 (2.6 × 10⁻³) | −0.0003 / +0.0016 | −0.0019 / −0.0030 | −0.0003 / +0.0001 |
| 0.005 | −0.0024 / −0.0010 (1.3 × 10⁻³) | −0.0007 / +0.0006 | −0.0009 / −0.0014 | −0.0008 / +0.0005 |
| 0.002 | +0.0000 / −0.0007 (5.0 × 10⁻⁴) | +0.0003 / +0.0007 | −0.0007 / −0.0008 | +0.0002 / +0.0005 |
| 0.001 | +0.0001 / +0.0002 (2.5 × 10⁻⁴) | +0.0004 / +0.0001 | −0.0009 / −0.0005 | +0.0003 / −0.0002 |

The Euler–Maruyama errors fall roughly in proportion to Δt and reach the statistical level (5 × 10⁻⁴) near νΔt = 0.002. Without reflection the ⟨μ⟩ error would be μ₀[(1 − νΔt)^n − e^{−νt}] = −0.0096 at νΔt = 0.1 (F-PUSH-18c), so about two thirds of the measured −0.0289 comes from the boundary. The exact rotation-kernel errors are −0.0032 / +0.0002 at νδt = 0.1 and below 3.1 × 10⁻⁴ for νδt ≤ 0.01 (F-PUSH-18a); the measured values agree within about two standard errors. The Milstein scheme reduces the ⟨μ⟩ error and the number of reflections but gives a ⟨μ²⟩ error of opposite sign and similar size. Running κ(t) at νΔt = 0.01 (MCR-3, N = 2 × 10⁵, v = ν = 1): Euler–Maruyama 0.3146 ± 0.0010 at νt = 50 against 0.32667 from Eq. (175) (3.7% low), the rotation kernel 0.3273 ± 0.0010.

Without a boundary, Euler–Maruyama would give ⟨μ⟩ = μ₀(1 − νΔt)^n and a stationary ⟨μ²⟩ = 1/(3 − νΔt); Milstein 1/(3 − 1.5νΔt) (derived recursions, F-PUSH-09). With the mirror rule from a uniform start (νt = 3; MCR-1b) the measured stationary bias of ⟨μ²⟩ is:

| νΔt | Euler–Maruyama (no-boundary value) | Milstein (no-boundary value) | rotation kernel |
|---|---|---|---|
| 0.1 | −0.0223 (+0.0115) | +0.0155 (+0.0175) | −0.0003 |
| 0.05 | −0.0134 (+0.0056) | +0.0089 (+0.0085) | +0.0006 |
| 0.02 | −0.0065 (+0.0022) | +0.0037 (+0.0034) | −0.0005 |
| 0.01 | −0.0034 (+0.0011) | +0.0014 (+0.0017) | +0.0002 |
| 0.005 | −0.0022 (+0.0006) | +0.0010 (+0.0008) | −0.0001 |

(N = 10⁶, standard error 0.0003. The Euler–Maruyama entry at νΔt = 0.01 is the mean of two independent runs time-averaged over 3 ≤ νt ≤ 13 or longer, −0.0033 ± 0.0001 and −0.0034 ± 0.0002; Section 25.3.) With reflection the Euler–Maruyama bias changes sign and is about three times the no-boundary value: 2.7, 3.0 and 3.4 times at νΔt = 0.02, 0.01 and 0.005 in time-averaged runs (2.9, 3.0 and 4.0 from the single-time table entries, whose standard error is 0.0003); the Milstein bias stays close to its no-boundary value because few particles are reflected; the rotation kernel has no stationary bias (Section 16.7).

### 16.7 Alternative kernels

**Small-angle random rotation (Àgueda 2008; Wijsen et al. 2019).** For isotropic scattering with frequency ν, the scattering-angle distribution after δt is (thesis Eq. 3.17; Wijsen et al. Eq. 16, which prints the same density in dϑ², "after a certain time δt≪ν⁻¹")

$$
F(\vartheta,\phi,\delta t)\,d\Omega=\frac{1}{2\pi\nu\delta t}\exp\!\Big(-\frac{\vartheta^2}{2\nu\delta t}\Big)\vartheta\,d\vartheta\,d\phi ,
\tag{167}
$$

so ϑ² is drawn from an exponential distribution with mean 2νδt and φ uniformly (thesis: "we pick the values of ϑ² and φ in each scattering process from exponential and uniform distributions, respectively"), and

$$
\mu_{\rm new}=\mu_{\rm old}\cos\vartheta+\sqrt{1-\mu_{\rm old}^2}\,\sin\vartheta\cos\varphi
\tag{168}
$$

(Wijsen et al. Eq. 17, two renderings agree; thesis Eq. 3.18 with prime and root garbled; the journal paper of Àgueda et al. [P61] was not read; Vainio & Kocharov Eq. A.1 and Virtanen & Vainio Eq. 24 print the same update). Equivalence with Eq. (21) (derived): with ⟨ϑ²⟩ = 2νδt, ⟨Δμ⟩ = −νμδt and ⟨Δμ²⟩ = ν(1 − μ²)δt to first order, i.e. D_μμ = (ν/2)(1 − μ²) (Eq. (24)) with the Itô drift ∂D_μμ/∂μ = −νμ, as Vainio & Kocharov state (⟨Δμ_±⟩/Δt = ∂D_μμ/∂μ = −μ_±ν_±). Exact one-step moments (derived; m = 2νδt, F the Dawson function):

$$
C_1\equiv\langle\cos\vartheta\rangle=1-\sqrt m\,F\!\big(\tfrac{\sqrt m}{2}\big),\quad
C_2\equiv\langle\cos2\vartheta\rangle=1-2\sqrt m\,F(\sqrt m),\quad
\langle\Delta\mu\rangle=\mu(C_1-1),\quad
\langle\Delta\mu^2\rangle=\mu^2\big(\tfrac{1+C_2}{2}-2C_1+1\big)+\tfrac{(1-\mu^2)(1-C_2)}{4}.
\tag{169}
$$

At ν = 1, δt = 0.01, μ = 0.3: ⟨Δμ⟩ = −0.0029900199714602886, ⟨Δμ²⟩ = 0.0089974528863896675 (F-PUSH-05; a 2-D quadrature agrees to 15 digits), against first-order −0.003 and 0.0091. By the addition theorem the kernel multiplies ⟨P_l(μ)⟩ by ⟨P_l(cos ϑ)⟩ per step, so the isotropic distribution is stationary exactly and ⟨μ⟩ decays with the effective rate −ln C₁/δt = ν(1 + νδt/6 + …) (1.0173641 for νδt = 0.1, 1.0016734 for 0.01; derived). One-step MC (MCR-6, N = 10⁷, νδt = 0.01) agrees with Eq. (169) within 1.2 standard errors at μ = 0.3 and 0.8.

**Hemispheric mapping for μ-dependent ν.** For ν(μ) = ν₀|μ|/(1 + |μ|) (thesis Eq. 3.19) the thesis maps each hemisphere onto [−1, 1] (Eq. 3.20): η = 2(μ − ½) for μ ≥ 0, η = 2(μ + ½) for μ < 0, with

$$
D_{\eta\eta}=\Big(\frac{d\eta}{d\mu}\Big)^2D_{\mu\mu}=\tfrac12\nu_0(1-\eta^2),
\tag{170}
$$

so scattering in η is isotropic with rate ν₀ (derived check F-PUSH-12 at three μ; because dη/dμ = 2 is constant, the operator ∂_μ(D_μμ∂_μ·) maps to ∂_η(D_ηη∂_η·) without extra drift). The back-transformation is not spelled out in the extracted thesis text; the inverse of Eq. (170) is μ = (η + 1)/2 (μ ≥ 0) and μ = (η − 1)/2 (μ < 0) (derived). The ε term (thesis Eq. 3.21; Wijsen et al. Eq. 8) is handled, according to Wijsen et al., "by means of two isotropic scattering processes": D^iso_μμ = ε(ν/2)(1 − μ²) (Eq. 14, ar5iv alttext read twice; ε lost in the PDF extraction) and D^aniso_μμ = (ν/2)(1 − μ²)|μ|/(1 + |μ|) (Eq. 15); their sum is their Eq. 8, Eq. (30) here (derived). The order of the two sub-steps is not stated. MCR-6 (ε = 0.048, ν₀ = 1, N = 10⁷; rotation with rate εν₀, then rotation in η with rate ν₀ in the current hemisphere, then the inverse map) reproduces ∂D_μμ/∂μ and D_μμ of the ε-form to first order (ν₀δt = 0.01: drift within 1.5 standard errors, ⟨Δμ²⟩/(2δt) 1.4% and 0.2% below D_μμ at μ = 0.3 and 0.8).

**Poisson-timed isotropic redraws.** Marsh et al. [P147] (one rendering) implement isotropic scattering in the solar-wind frame with events "Poisson distributed in time with an average scattering time" t_scat = λ/v₀ (v₀ the initial particle speed); at each event "the direction of the particle's velocity vector is randomly reassigned". Kelly et al. [P148] use λ/v "where v is the particle's velocity in a frame of reference which is stationary with respect to the solar wind" and transform to that frame, redraw, and transform back. With redraw rate r = 1/t_scat the generator is r(⟨·⟩_iso − ·), so (derived)

$$
\langle\mu\rangle(t)=\langle\mu\rangle_0e^{-rt},\qquad \langle\mu^2\rangle(t)=\tfrac13+\big(\langle\mu^2\rangle_0-\tfrac13\big)e^{-rt},\qquad \kappa_\parallel=\frac{v^2}{3r}=\frac{v\lambda}{3}.
\tag{171}
$$

The first moment and κ∥ coincide with isotropic diffusion with ν = r (MFP v1.2 Section 6.6), but ⟨μ²⟩ relaxes at rate r instead of 3ν; the redraw is a jump process and not a realization of Eq. (21). MCR-4 (r = 1, μ₀ = 0.5, t = 1, N = 10⁶): ⟨μ²⟩ = 0.3022 ± 0.0002 against 0.30268 from Eq. (171) and 0.32918 for diffusion; κ(10) = 0.2991 ± 0.0009 against 0.30000.

**Random walk on a sphere.** Florinski [P46] (arXiv version, two renderings agree) writes isotropic scattering as Sf = D∂_μ[(1 − μ²)∂_μf] + D/(1 − μ²)∂²f/∂φ² (Eq. 2) with λ∥ = vτ, τ = (2D)⁻¹, and integrates it by the random walk on a sphere; "For small scattering angles the process is identical to a two-dimensional random walk on a plane." In the notation of Eq. (167) this is ν = 2D (derived: λ∥ = v/ν).

### 16.8 Frames

**Where scattering is done (as stated).**

| Frame of μ (and v) in which D_μμ is applied | Sources |
|---|---|
| local solar-wind frame, propagation in a fixed or corotating frame | [P55] ("The scattering is performed in the local solar wind frame, where the magnetic fluctuations are assumed to be static"); [P7]; [P148]; [P147] |
| mixed coordinates: x, t in the observer's frame, p, μ in the solar-wind frame | [P6]; [P25]; [P1]; [P16] (position in the corotating frame, momentum in the solar-wind frame) |
| "non-relativistic flow frame" for p and μ | [P14] |
| Alfvén-wave rest frame(s) | [P78]; [P83]; [P88] (both wave frames each step); [P91] (forward-wave frame, then backward-wave frame) |
| Lagrangian frame co-moving with the plasma (grid) | EPREM [P34, P17, P32]; [P10] (field lines advected with the solar wind) |
| spatial coordinates corotating; frame of μ, p not extracted | [P15] (Zhang et al. 2009 abstract: "a fixed reference frame in which any spacecraft is roughly stationary") |

**Transformations printed.** Àgueda's thesis (Eq. 3.14, one rendering, partly garbled) transforms from the corotating to the solar-wind frame with γ′ legible as γ_u(γ − u_cμv/c²) and μ′ legible as [?]γ_uγ′(μv − u_c)/v′ (prefactor not legible), γ_u = (1 − u_c²/c²)^{−1/2}, u_c = u r̂ − Ωr φ̂ = u_c ẑ (Eq. 3.15) and u_c = u sec ψ(z) (Eq. 3.16). The standard Lorentz transformation of (v, μ) to a frame moving with speed u along the field is (derived; round trip exact at 40 digits, F-PUSH-15)

$$
\gamma'=\gamma_u\gamma\Big(1-\frac{u v\mu}{c^2}\Big),\qquad \mu'=\frac{\gamma_u\gamma\,(v\mu-u)}{\gamma' v'},\qquad \mu'\approx\mu-\frac{u}{v}(1-\mu^2)\ \ (u\ll v),
\tag{172}
$$

so the legible γ′ of Eq. 3.14 lacks a factor γ on its second term relative to this form; whether this is in the print or in the extraction is not known. Kartavykh et al. describe "proper Lorentz transformations between a stationary frame in which the Sun is at rest" and frames moving with the local solar wind, and energy changes that arise from them (no explicit dp/dt term). Ruffolo (1995, Eqs. 7–9) gives the mixed-frame rates ṗ′, ż and μ̇′ with first-order v_sw/v terms and multiplies the scattering term by (1 − μv_sw v sec ψ/c²) (E′/E; tool assembly of Eq. 11); Minoshima et al. [P64] print the operator acting on (1 − μvV^c_sw/c²)f. For wave-frame scattering, the plasma-frame energy change per scattering, −ΔE = ∓V_A p_±Δμ_± [P88], and the wave-frame D_μμ forms are in Sections 8 and 14.

**Clock of the coefficient (derived).** Transforming the particle alone does not transform a rate. A D_μμ defined per unit time of a frame moving with speed u along the field is applied with that frame's elapsed time along the particle path, which for a particle of speed v and pitch cosine μ in the first frame is

$$
dt'=\gamma_u\Big(1-\frac{u\,v\,\mu}{c^2}\Big)dt,\qquad \gamma_u=(1-u^2/c^2)^{-1/2}
\tag{173}
$$

(from t′ = γ_u(t − uz/c²) with z = vμt; F-PUSH-19). The sources that scatter in a wave or solar-wind frame (table above) do not state which clock their D_μμ refers to, except where the mixed-frame transport equations carry the factor explicitly (Ruffolo's (1 − μv_sw v sec ψ/c²), Minoshima et al.'s (1 − μvV^c_sw/c²), Sections 2.3 and 16.8). Each configured coefficient therefore declares its clock (host time, plasma-frame time, wave-frame time, or a mixed-frame time already contained in the source equation) and the rate conversion is applied once; a first-order mixed-frame operator is labelled by its retained order and not as an exact Lorentz-transformed scattering rate. A wave-frame operator that is elastic in its frame, combined with the transformations above, generates the plasma-frame energy change; adding a separate D_pp for the same wave population counts that change twice (WAVE-A2-VV05, Section 8.2). Where D_μp and D_pp are carried, the tensor is positive semidefinite: D_μμ ≥ 0, D_pp ≥ 0, D_μp² ≤ D_μμD_pp.

### 16.9 Validation of a pusher

**Analytic moment tests (derived).** For pure scattering with D_μμ = (ν/2)(1 − μ²), d⟨μ⟩/dt = ⟨D′⟩ = −ν⟨μ⟩ and d⟨μ²⟩/dt = ⟨2μD′ + 2D⟩ = ν − 3ν⟨μ²⟩, so

$$
\langle\mu\rangle(t)=\langle\mu\rangle_0e^{-\nu t},\qquad
\langle\mu^2\rangle(t)=\tfrac13+\big(\langle\mu^2\rangle_0-\tfrac13\big)e^{-3\nu t}.
\tag{174}
$$

F-PUSH-06 compares these with the Legendre-series solution of Tautz [P111] (Eq. 11, with D = ν/2) at μ₀ = 0.5, νt = 1: both give 0.18393972058572116 and 0.32918441096934467. Tautz's Eqs. 14 and 18 are the same results with ν = 2D. The running spatial coefficient from an isotropic stationary start is (derived from ⟨μ(0)μ(t)⟩ = e^{−νt}/3)

$$
\kappa(t)=\frac{\langle\Delta z^2\rangle}{2t}=\frac{v^2}{3\nu}\Big[1-\frac{1-e^{-\nu t}}{\nu t}\Big]\ \longrightarrow\ \frac{v^2}{3\nu}=\frac{v^2}{6D_0},
\tag{175}
$$

the limit printed by Effenberger & Litvinenko (Eq. 10) and Merten & Aerdker (Eq. 38) (PARALLEL rev. 1.4, Section 7.2) (values 0.1226265, 0.3000015, 0.3266667 at νt = 1, 10, 50 for v = ν = 1; F-PUSH-07). At νt = 50 κ(t) is still 2% below the limit.

**Focusing benchmarks.** Effenberger & Litvinenko print κ∥ = λ₀v(coth ξ/ξ − 1/ξ²), ξ = λ₀/L (Eq. 8), and compare SDE densities with diffusion (Eq. 7) and telegraph (Eqs. 15–19) solutions. For constant focusing with isotropic scattering, the stationary μ-distribution of Eq. (150) is ∝ e^{ξμ} with ξ = v/(2LD₀) (derived), the form of van den Berg et al. Eqs. 47–48 (f̃_iso = (ξ csch ξ/2)e^{ξμ}, A₁ = 3(coth ξ − 1/ξ)); hence

$$
\langle\mu\rangle_{\rm stat}=\coth\xi-\frac1\xi .
\tag{176}
$$

For Merten & Aerdker's test (v/(2L) = 5, D₀ = 1, so ξ = 5) this is 0.80009080398201938, whereas the zero of the Itô drift, 5(1 − μ²) − 2μ = 0, is μ = 0.81980390271855697, 2.5% higher (F-PUSH-08). Merten & Aerdker (Sect. 3.5.1, one rendering) describe the two fixed points of their SDE with focusing ("the previous stable one is pushed to positive values of μ", with "an unstable fixed point at μ < 0") and state that "the drift velocity is determined by the stable pitch-angle: v_drift = μ₁c₀"; μ₁ is not printed numerically and c₀ is not defined. This description fits the zero of the deterministic drift (0.8198 here), not the stationary mean of the process (0.8001); the identification is this review's reading. MCR-5 (mirror rule, uniform start, N = 2 × 10⁵, t = 3) gives ⟨μ⟩ = 0.7874, 0.7956, 0.7975, 0.7996 (Euler–Maruyama) and 0.7981, 0.8000, 0.7999, 0.7996 (Milstein) at Δt = 0.01, 0.003, 0.001, 0.0003 (standard error 0.0005): Eq. (176), not the drift zero.

**Which focused-transport observable is compared with which reference (derived, Section 12.3).** For constant focusing with isotropic scattering the forward process has three distinct references: the stationary angular density ρ₊ ∝ e^{ξμ} (Eq. (176), Eq. (145)), the mean streaming speed v⟨μ⟩ = v(coth ξ − 1/ξ) (equal to u = κ∥/L of Eq. (112)), and the long-time centred spatial variance Var[z]/(2t) → κ_var of Eq. (115), which is smaller than Eq. (112) for ξ ≠ 0 (0.248 against 0.313 for v = ν = 1 at ξ = 1). A backward run samples ρ₋ ∝ e^{−ξμ} in its pitch cosines and reconstructs f from the end-point values; its histogram is compared with ρ₋, its reconstructed solution with the manufactured f₀ = e^{−z/L}e^{ξμ} (F-FOC-05). MCR-8 (Euler–Maruyama, mirror rule, N = 10⁵, Δt = 0.004, uniform start, v = ν = 1, ξ = 1.5; group 13 of the physics-check script): the slope of Var[z(t)]/2 between t = 20 and t = 40 is 0.1776 ± 0.0012 (standard error from 20 particle batches), against κ_var = 0.17828 and Eq. (112) = 0.29208; the mean streaming speed in the same run is 0.4359 ± 0.0008 against coth ξ − 1/ξ = 0.4381. The same scheme in backward time (drift −v(1−μ²)/(2L) + D′_μμ, t = 8) gives ⟨μ⟩ = −0.4344 ± 0.0015 against −0.4381, the sign of ρ₋. Validation references VAL-PUSH-KVAR and VAL-PUSH-MMS.

**Time-dependent rate (derived).** For an imposed isotropic rate ν(t) and no other angular dynamics, the operator (ν(t)/2)∂_μ[(1−μ²)∂_μ·] has the Legendre polynomials as eigenfunctions with eigenvalues −ν(t)l(l+1)/2, so

$$
\langle P_l(\mu)\rangle(t)=\langle P_l(\mu)\rangle(0)\,\exp\Big[-\frac{l(l+1)}{2}\int_0^t\nu(t')\,dt'\Big]
\tag{177}
$$

(F-PUSH-20; Eq. (174) is the case l = 1, 2 with constant ν). It is an exact test for a host that updates the scattering rate between steps: the moments follow the integrated rate, with no contribution from a stale snapshot or from an extra kick at the update (VAL-PUSH-LEG).

**Derivative-free reference operator (derived).** A conservative finite-volume discretization of Eq. (21) on cells of width Δμ_i gives a continuous-time Markov chain with transition rates Q_{i→j} = T_{ij}/Δμ_i between neighbouring cells, where the transmissibility between the cell centres μ_i and μ_j is computed from the reciprocal coefficient,

$$
T_{ij}=\Big[\int_{\mu_i}^{\mu_j}\frac{d\mu}{D_{\mu\mu}(\mu)}\Big]^{-1},
\tag{178}
$$

not by sampling D_μμ at the interface. The chain conserves the number of particles, its detailed balance preserves the uniform density (zero-flux ends), and its continuum limit is Eq. (21). For D_μμ ∝ |μ|^α near μ = 0 the integral across μ = 0 is finite for α < 1 and infinite for α ≥ 1 (the same condition as for Eq. (23)): the q-form with H = 0 and 1 < q < 2 therefore communicates across μ = 0 with a finite transmissibility (F-PUSH-21: T = 0.452241 D₀ between μ = −0.05 and μ = 0.05 for q = 5/3, H = 0), whereas an interface value D_μμ(0) = 0 would disconnect the hemispheres; the ε-form with ε = 0 (α = 1) is disconnected in both treatments, consistent with Section 16.5. The chain needs no derivative of D_μμ and is the reference for coefficients whose drift is unbounded or discontinuous (Section 16.11); its angular-grid convergence is demonstrated, since it remains an approximation (VAL-PUSH-FV). The physics-check script of Section 24 uses this discretization for the Poisson equation behind Eq. (115).

**Estimator consistency with Eq. (22).** Ivascenko et al. [P108] (Eq. 6, two renderings) define D_μμ = lim (Δμ)²/(2Δt), the normalization of Eq. (22). One-step moments of a pusher are biased at O(Δt): Euler–Maruyama gives ⟨Δμ²⟩/(2Δt) = D_μμ + A²Δt/2 exactly, the rotation kernel 0.449873 instead of 0.455 at ν = 1, νδt = 0.01, μ = 0.3 (F-PUSH-04a/04b, -05). The running estimator is biased at late times: for D_μμ = D(1 − μ²), Tautz (Eq. 22) gives ½d⟨(Δμ)²⟩/dt = (1 − 3μ₀²)De^{−6Dt} + 2μ₀²De^{−2Dt}, which equals (1 − μ₀²)D only for small t (Eq. 23, printed condition "t ≪ D"). Van den Berg et al. [P106] compare nine estimators and find the mean-square-displacement estimators valid only in weak turbulence and within a time window; for their Eq. 14, T_μ = (1 − μ²)/(2D_μμ), the isotropic process gives T_μ = 1/ν = 1/(2D₀) exactly (derived from Eq. (174)); the factor 2 that they applied in weak turbulence is therefore not needed for isotropic scattering (Section 16.10). Ivascenko et al. tested their estimators with a Monte Carlo code whose input is the ε-form (ν₀, ε not found).

**Validation reported in the literature.**

| Source | Test |
|---|---|
| [P2] | SDE densities (10⁷ particles) vs analytic diffusion and telegraph solutions, ξ = 0 and 1.5 |
| [P11] | finite-difference code reproduces the SDE results of Effenberger & Litvinenko "very accurately" (Fig. 12) |
| [P55] | steady-state PADs vs analytic Fokker–Planck solution (Eq. 3.24) for isotropic and ε = 0.10 scattering; σ⁻¹ flux scaling; mean solar-wind-frame energy vs Eq. 3.41; radial indices vs the finite-difference code of Lario et al. (2007) |
| [P9] | isotropization of a μ₀ = 1 injection: "After a short time (∼2 s) the pitch angle diffusion leads to an isotropic distribution"; constant focusing (v/2L = 5, D₀ = 1); κ∥ = v²/(6D₀) asymptotics (no convergence number in the text read) |
| [P44] | SDE vs finite differences for spatial diffusion: agreement only for small enough Δt or with Eq. 51 |
| [P7] | stationary-shock spectrum approaches E^{−5/4} for R = 3; agreement up to 0.5 AU upstream |
| [P31] | simulated spectral index "always essentially the same as the prediction" γ = 3s/(s − 1) |
| [P6] | v_sw = 0 limit vs Ruffolo (1991); deceleration-only test vs τ_d |
| [P14], [P1] | SDE verification in Appendix B of each: Not found (not reachable); the primer code's "SDE_benchmark" folder compares with Dröge et al. (2010) (folder not readable) |

### 16.10 Printed issues in code descriptions that must not be copied

| Item | Location | As printed | Issue |
|---|---|---|---|
| No ∂D_μμ/∂μ drift | [P10] Eq. 4 | dμ = [(1 − μ²)/2 [v/L + μ(…)]]dt + √(2D_μμ)dW_μ(dt) | with the operator ∂_μ(D_μμ∂_μf) of their Eq. 1 the Itô drift needs ∂D_μμ/∂μ (Eq. (150)); whether the code includes it is not stated |
| Wiener increment | [P10] | "dW_μ(dt) = dt · Λ", Λ Gaussian (read twice, one rendering) | a Wiener increment scales as √dt·Λ |
| Missing μ_t | [P2] Eq. 30 | μ_{t+Δt} = [½ξ(1−μ_t²) − μ_t]Δt + … | the leading μ_t is absent in both renderings (Eq. (158) vs (159)) |
| Undefined M, no weight | [P20] (also [P22], [P21]) | +2D_μμ/(M + μ) in the backward μ drift | M undefined; no loss rate or weight printed (compare Eqs. (153), (154)) |
| Sign of the focusing term | [P20] Eq. 3 | −(v/B)∂B/∂z = v/L with L = (ẑ·∇ln B)⁻¹ | the two printed forms differ in sign; with the printed L, v/L is negative for a field decreasing along +z; the backward SDE, Eq. (155), uses V/L with the same L |
| Unit of D₀ | [P9] Sect. 3.5 | `D_0=1\,\mathrm{N}^2\,\mathrm{s}` | D_μμ has unit s⁻¹ |
| Reflection | [P83] | "μ̃ → 1 − μ̃" (as reported, one rendering) | not a mirror about μ = ±1 (Eq. (161) gives 2 − μ̃) |
| Time step | [P91] App. A | "Δt = a ν_e(γ)", a ≲ 0.01 (extraction) | dimensionally Δt = a/ν_e is expected; extraction-level |
| Eq. 14 factor | [P106] | Eq. 14 "had to be multiplied by 2 to match QLT" in weak but not strong turbulence | unresolved; for isotropic scattering Eq. 14 holds without a factor (Section 16.9) |
| Condition | [P111] Eq. 23 | "t ≪ D" | presumably t ≪ 1/D |
| Lorentz factor | [P55] Eq. 3.14 | γ′ = γ_u(γ − u_cμv/c²) (as legible) | standard form has γ on the second term (Eq. (172)); prefactor of μ′ not legible |

### 16.11 Implementation requirements

**Implementation requirement (operator and calculus).** A pitch-angle step that claims to realize Eq. (21) uses the Itô drift dμ/dt + ∂D_μμ/∂μ and the noise √(2D_μμ) in forward time, or the backward rule of Eq. (152), with D_μμ the coefficient of ∂_μ(·∂_μf); a code whose operator carries φ/2 or D_μμ/2 converts by Eq. (151) and states it. Tenishev et al. Eq. 4 and Effenberger & Litvinenko Eq. 30 are not implemented as printed (Section 16.10).

**Implementation requirement (sampled density).** The configuration states which density the forward pseudo-particles sample as an explicit coordinate measure composed of a spatial factor (d³x, or the flux-tube measure A(s)ds with A ∝ 1/B) and a momentum factor (dp dμ with or without the gyrotropic p², or the Ruffolo F convention), together with the weight or Jacobian that converts it to f (Eq. (150); Effenberger & Litvinenko Eq. 3; Qin et al. Eq. 7): in a static tube with compatible frames the number of particles is dN = 2πA(s)p²f ds dp dμ, and in three dimensions dN = 2πp²f d³x dp dμ; a mixed-frame model carries the measure and rate factors of its own transport equation (Section 2.3). The macroparticle weights, cell volumes or segment lengths and the sampled observables of the host represent that measure (Section 26, PS06b). A backward SDE samples f; its pitch-cosine histogram is a property of the sampling process (Section 12.3).

**Implementation requirement (coefficient–kernel compatibility).** Readiness is a property of a (coefficient, kernel) pair, not of the coefficient alone: a positive D_μμ(0) (H > 0) does not bound the drift of the q-form for 1 < q < 2, which Euler–Maruyama and Milstein (through the same derivative) both use (Section 16.5). A pair is released only after the tests of PS05 (Section 26) for that pair: convergence of the pitch-angle distribution near μ = 0 with the step, occupation of |μ| < 0.01, decay of several Legendre moments (Eq. (177)) and the finite-time spatial variance, each against the derivative-free reference of Eq. (178) or an exact solution, at a predeclared error budget. The coefficient interface therefore returns, besides D_μμ and the analytic one-sided derivatives, the regularity at μ = 0 and the domain of validity; a derivative-free kernel may be paired with an appropriate coefficient that has no finite derivative, after its pair tests pass; random rotation is restricted to its matching isotropic or source-specific hemispheric law, whereas the finite-volume reference realizes the declared scalar flux operator. An unbounded derivative is not replaced by zero, clamped or smoothed unless the modified coefficient is declared as a separate model and its convergence toward the intended operator is demonstrated.

**Implementation requirement (importance sampling).** Importance sampling of the Zhang et al. type carries the drift, the loss rate c and the factors of Eq. (154) together and records a. The He & Wan drift 2D_μμ/(M + μ) is not used without a defined M and a weight (Section 16.2).

**Implementation requirement (boundary).** Euler–Maruyama, Milstein and predictor–corrector steps apply a reflecting rule at μ = ±1 (Effenberger & Litvinenko; Merten & Aerdker; van den Berg et al. 2026), implemented as Eq. (161) applied until the value is in [−1, 1]; the rule and the reflected fraction per step are logged. The random-rotation and redraw kernels need no rule. A shape with D_μμ(±1) ≠ 0 (the printed (1 − μ)² variant, Section 5) needs the no-flux condition imposed explicitly (MFP v1.2, Section 21.5).

**Implementation requirement (μ = 0).** For the q-form and the |μ|^{2/3} shape the drift of Eq. (163) is unbounded at μ = 0 for every H; for the ε-form it is discontinuous there. The library evaluates Eq. (163) without division by zero, states its treatment of μ = 0 for every run (step-size rule, kernel or regularization, with its source), and demonstrates convergence of the μ-distribution near μ = 0 with the step size (Section 16.5, MCR-2). It does not insert an undocumented H, h or ε.

**Implementation requirement (time step).** The step rule is declared with its source: Eq. (165) (with ℓ), Eq. (166) (with Δμ_max), a fixed Δt with a convergence study, or δt ≪ ν⁻¹ for the rotation kernel. Convergence of ⟨μ⟩ and ⟨μ²⟩ to Eq. (174) and of κ(t) to Eq. (175) is shown at the chosen step (Section 16.6 gives the expected size of the errors).

**Implementation requirement (frames).** D_μμ is evaluated with μ and v in the frame in which its source defines it (Section 16.8, and Sections 8 and 14 for wave frames); the transformation to that frame and back is stated.

**Implementation requirement (kernels are distinct models).** The random rotation (Eq. (168)), the hemispheric two-process split (Eq. (170)) and the Poisson redraw (Eq. (171)) are separate kernels with their own IDs. The redraw is not a pitch-angle diffusion and is validated by Eq. (171), not by Eq. (174).

**Implementation requirement (estimators).** D_μμ measured from a pusher uses the normalization of Eq. (22) and accounts for the O(Δt) one-step bias and the late-time bias of the running estimator (Section 16.9).

### 16.12 Normative execution contract added in version 1.2

This subsection specifies the operator realized by an executable model. It supplements the source transcriptions without altering them. Equations N1–N7 are derived from the operators and transformations already given in this document. A source-specific approximation remains subject to its own domain. Physical admissibility, numerical convergence and observational validation are separate requirements.

**Scalar scattering and ownership of transport processes.** In a local frame in which scattering is elastic, the library realizes only

$$
\mathcal C f=\partial_\mu(D_{\mu\mu}\partial_\mu f),\qquad
d\mu=D'_{\mu\mu}\,dt_c+\sqrt{2D_{\mu\mu}}\,dW_c .
\tag{N1}
$$

The clock $t_c$ is the coefficient's declared clock. The host supplies streaming, deterministic focusing, flow gradients, adiabatic momentum change, spatial drift and shock crossings through an audited transport composer. For the static-tube reference, $ds=v\mu\,dt$ and $a_\mu=v(1-\mu^2)/(2L)$ are added once. Forward and backward deterministic terms follow the selected transport equation; backward integration does not change the sign of the diffusion term. A population already undergoing pitch-angle scattering and field-aligned streaming receives no additional parallel spatial Brownian displacement representing the same scattering. A full-orbit pusher that already resolves the magnetic mirror force receives no duplicate focusing drift.

**Full momentum-space diffusion.** A source providing $D_{\mu p}$ and $D_{pp}$ cannot be implemented through a scalar coefficient interface. For the symmetric tensor $D_{p\mu}=D_{\mu p}$ and the gyrotropic phase-space density, the operator of Section 9.1 is

$$
\mathcal C f=\partial_\mu(D_{\mu\mu}\partial_\mu f+D_{\mu p}\partial_p f)
+\frac{1}{p^2}\partial_p\left[p^2(D_{p\mu}\partial_\mu f+D_{pp}\partial_p f)\right].
\tag{N2}
$$

For forward particles sampling $P=p^2f$ with respect to $dp\,d\mu$, at fixed spatial coordinates and with no other transport terms, the Itô drift is

$$
A_\mu=\partial_\mu D_{\mu\mu}+\partial_pD_{\mu p}+\frac{2D_{\mu p}}p,
\qquad
A_p=\partial_\mu D_{p\mu}+\partial_pD_{pp}+\frac{2D_{pp}}p,
\qquad BB^{\mathsf T}=2D .
\tag{N3}
$$

The two stochastic increments are correlated. Omitting the $2D/p$ terms samples another measure; drawing independent noises generally discards the cross coefficient. More generally, a pure scattering operator $M^{-1}\partial_i(MD_{ij}\partial_j f)$ yields $A_i=\partial_jD_{ij}+D_{ij}\partial_j\ln M$ for the forward coordinate density $P=Mf$. A time-dependent or mixed-frame transport measure also has contributions from the complete transport equation and is not obtained by substituting $M$ into this static identity alone. A backward reconstruction uses the backward generator of that complete equation, its terminal/source data and any required weight.

All tensor entries and all derivatives needed by the chosen kernel are returned together. In SI, $D_{\mu\mu}$ has units s$^{-1}$, $D_{\mu p}$ momentum/time, and $D_{pp}$ momentum$^2$/time. Positive semidefiniteness requires $D_{\mu\mu}\ge0$, $D_{pp}\ge0$ and $D_{\mu p}^2\le D_{\mu\mu}D_{pp}$. A rank-one tensor is valid. If $D_{\mu\mu}=0$, positivity requires $D_{\mu p}=0$; the factorization must handle this case without division by zero. Check positivity after nondimensionalizing both coordinates, since directly comparing eigenvalues with different physical units is invalid. Only a documented roundoff correction at the factorization scale is allowed; it cannot repair a physically indefinite tensor.

Angular boundaries obey the no-normal-flux condition of the full tensor operator. If $D_{\mu\mu}$ vanishes at an endpoint, PSD forces $D_{\mu p}$ to vanish there as well. If a model instead has nonzero endpoint/cross diffusion, scalar mirror folding of μ alone is not automatically a reflecting realization of its tensor operator. Derive the reflected boundary dynamics from that full generator and verify the coupled flux condition; momentum cutoffs likewise have their own specified boundary. Scalar boundary rules cannot be copied to a tensor kernel without this check.

For $\eta=\ln(p/p_*)$, $p>0$, Itô's formula gives $A_\eta=A_p/p-D_{pp}/p^2$, $D_{\mu\eta}=D_{\mu p}/p$ and $D_{\eta\eta}=D_{pp}/p^2$. Logarithmic momentum can preserve positivity, but its extra drift is compulsory. The value $p=0$ has no pitch-angle direction and is outside this coordinate system. A momentum cutoff or thermalization boundary is an explicit application boundary, not a numerical floor.

**Domains and coefficient status.** The four admissibility values of Sections 4.1 and 17 are distinct. `AdmissibleOnDomain` is executable when its domain predicate passes and the requested transport problem is completely covered. It is not rejected merely because its enum differs from `Admissible`. A formula valid only at $\mu>0$ or away from $|\mu|=1$ is a valid local evaluator on that domain, but is not a complete scattering law on $[-1,1]$. Completing it requires a separately documented physical law and matching conditions; replacing fractional powers by absolute values, endpoint clipping or a fitted bridge is a model choice with its own provenance. A failed domain predicate returns no numerical coefficient.

A prescribed finite band and unavailable measured power are different cases. A physically declared truncated spectrum may produce zero scattering or an internal barrier; missing spectral coverage returns `MissingCoverage`. Mean-free-path integration returns `Finite`, `Divergent`, `Unavailable` or `Nonconverged`, with an error estimate when finite. It never substitutes a large finite number for infinity. The tag is `lambda_parallel_unfocused`: Eq. (23) does not include focusing or replace a coupled momentum-space closure. A radial mean free path, a focusing closure or an observational fit with another convention is not accepted as this tag without its named conversion.

**Singular points and the reference operator.** For a connected scalar divergence operator with an integrable internal resistance, Eq. (178) selects transmission through the singular point. The internal extension and continuity/flux conditions are part of the model record; the formula for $D$ alone need not specify an extension of a degenerate differential operator. For the q-form with $H=0$, $1<q<2$, retain the finite cross-zero conductance. For the pure-scattering ε-form with ε=0, preserve the disconnected hemispheres. A deterministic focusing term may move particles across zero even when the pure scattering operator cannot: the barrier test must state which operator it tests.

On a nonuniform grid, store cell widths, centers, transmissibilities, both outgoing rates and the diagonal. With rows of $Q$ denoting origin states, $Q_{ij}=T_{ij}/\Delta\mu_i$, $Q_{ii}=-\sum_{j\ne i}Q_{ij}$; the stationary probability masses are proportional to the widths, not equal per bin. The angular continuum limit is demonstrated with at least three refinements. A continuous-time chain is a discretized reference, not an exact continuum pusher. For constant focusing and isotropic scattering, a useful second reference has stationary mass $m_i=\int_{\text{cell }i}\rho_+(\mu)d\mu$ and conductance $T^\rho_{ij}=[\int_{\mu_i}^{\mu_j}(D\rho_+)^{-1}d\mu]^{-1}$, with $Q_{ij}=T^\rho_{ij}/m_i$. This follows from writing the angular flux as $-D\rho_+\partial_\mu(P/\rho_+)$ and independently checks Eq. (115).

**Streaming needs a trajectory integral.** Even an exact angular endpoint sampler does not determine $\Delta s=v\int\mu(t)dt$. For an isotropic stationary start in uniform physics, $v=\nu=1$ and $\Delta t=1$, the exact second moment of displacement is $2e^{-1}/3=0.245252961$, whereas streaming with a single stationary endpoint gives $1/3$. That 35.9% excess is a finite-step transport error despite perfect endpoint angular statistics. Stream along the subcycled trajectory; for the finite-state chain, integrate the dwell times between angular jumps. Resolve cell/segment exits, field rotation, shocks and background updates within the same interval. Compare finite-time displacement distributions and Eq. (175), including the ballistic limit, rather than only the asymptotic mean free path.

**Adaptive Brownian increments.** Select an accepted substep from known state and field information before drawing its noise whenever possible. If a step is shortened after its Brownian increment has been drawn, preserve that increment and its remaining path. For a parent interval $h$, split fraction $0<\theta<1$ and independent $Z\sim N(0,1)$,

$$
\Delta W_1=\theta\Delta W+\sqrt{\theta(1-\theta)h}\,Z,
\qquad \Delta W_2=\Delta W-\Delta W_1 .
\tag{N4}
$$

Discarding a rejected kick and resampling independent increments can bias the accepted noise. The bridge law is not a complete adaptive solver by itself: keep the unused increments, their physical intervals and their ordering, including across restart and migration. Multi-dimensional tensor noise needs one bridge per independent driving Wiener process. Rackauckas & Nie (2017), DOI 10.3934/dcdsb.2017133, provides a primary numerical-method reference for rejection sampling with memory. A geometry crossing located from a stochastic trial trajectory likewise requires a documented event algorithm; a blind reject-and-redraw is not authorized by the mirror boundary rule. Reflection is a discretization of the prescribed angular boundary, with measured bias and convergence, not an error-rejection mechanism.

**Frames, clock and state.** Store momentum and pitch cosine as authoritative in a named local frame, and derive $v$ and $\gamma$ from the same momentum and species. For kinetic energy use the stable identity $T=p^2c^2/(E+mc^2)$, $E=\sqrt{m^2c^4+p^2c^2}$. Use signed $\Omega=qB/(\gamma m)$ in charge-sensitive resonance terms and positive $r_L=p/(|q|B)$ as the maximal radius. A formula that omits $\gamma$ is executable only if its nonrelativistic approximation is declared and validated on the requested energies, or a separately identified corrected law is supplied. A source's omission is not permission to use $\Omega_0$ for relativistic GCRs.

An exact nonparallel Lorentz boost requires the full momentum vector, including gyrophase. The pair $(p,\mu)$ alone does not supply that vector, and a general boost does not preserve gyrotropy about the same axis. A 3-D gyrotropic AMPS adapter therefore keeps $(p,\mu)$ in the local plasma frame and uses a declared mixed-frame transport equation, unless the host actually supplies the extra vector information. Do not repeatedly sample an artificial gyrophase to claim an exact scalar observer-frame transformation. The scalar boost of Eq. (172) is limited to the along-field transformation for which it was derived. Apply the declared clock factor once; an already mixed-frame rate is not boosted a second time.

**A state-dependent clock multiplies the generator, not a newly differentiated coefficient.** For a coefficient-frame Itô process $dX=a_c(X)dt_c+B_c(X)dW_c$, with $dt_c=g(X)dt_h>0$, the path-parameter time change gives

$$
\mathcal L_h=g\mathcal L_c,\qquad a_h=g\,a_c,\qquad B_hB_h^{\mathsf T}=g\,B_cB_c^{\mathsf T}.
\tag{N7}
$$

For scalar pure scattering this means $A_h=gD'_c$ and $D_h=gD_c$. Using $(gD_c)'$ as the drift inserts the extra $g'D_c$ and changes the process. The Lorentz clock factor of Eq. (173) depends on the particle state, so this distinction matters. Transform the entire generator and state with Itô's formula, or scatter in the declared coefficient frame using its clock and a converged host-time subdivision; do not rescale D and then differentiate it as if the host collision measure were unchanged. A pure angular time change alone has stationary host-time density proportional to $\rho_c/g$, a useful mathematical test; the physical distribution of a full frame-transformed FTE additionally depends on its coordinate measure, spatial terms and boundary data. The static tensor drift Eq. (N3) applies in its stated frame/measure/clock; an arbitrary subsequent state-dependent time change is governed by Eq. (N7), not by reusing Eq. (N3) on rescaled entries without transforming the measure.

**Four wave populations.** Cross helicity and magnetic helicity are marginal ratios; they do not fix the correlation between direction and polarization. Let $j,n=\pm1$, $G=\sum g_{n,j}$, and introduce the explicitly supplied joint moment $\chi=\sum jn\,g_{n,j}/G$. The general four nonnegative weights satisfying Eq. (55) are

$$
g_{n,j}=\frac G4(1+j\sigma_c+n\sigma_m+jn\chi),\qquad
-1+|\sigma_c+\sigma_m|\le\chi\le1-|\sigma_c-\sigma_m| .
\tag{N5}
$$

The factorized split is the additional assumption $\chi=\sigma_c\sigma_m$. Supply the four spectra directly, or supply a named joint-split prescription and its domain; do not infer it from two helicities. At $G=0$, the ratios are undefined and no split is needed. Direction is defined relative to the signed magnetic unit vector, with an explicit mapping from outward-field conventions. Polarization conventions and the sign of the charge are source specific. The $V_A\to0$ balanced limit must reproduce the canonical slab coefficient with the same total transverse variance. Sharp resonances retain the actual Jacobian $|v\mu j-V_A|$; the replacement $v|\mu|$ is a separately domain-limited approximation near the Alfvénic resonances (Section 8.1).

**Particle–wave conservation and spectral channels.** An elastic angular kick at fixed $p'$ in a wave frame moving at signed along-field speed $u$ gives, exactly for that nondispersive frame transformation,

$$
\Delta E_{\rm particle}=\gamma_u u p'\Delta\mu',\qquad
\Delta p_{\parallel,\rm particle}=\gamma_u p'\Delta\mu',\qquad
\Delta E_{\rm particle}=u\Delta p_{\parallel,\rm particle} .
\tag{N6}
$$

Return the complete particle state and both exchanges in the plasma frame. In a closed nondispersive particle–wave test, deposit their negatives into the participating population, accounting for the host's statistical weight, volume or tube measure, and the chosen wave energy definition. Direction alone is an insufficient population key. Each exchange has an event ID, direction, polarization convention, population ID, spectral-bin ID and immutable spectrum epoch. For broad resonances, define the normalized spectral deposition weights from the same interaction kernel used in the scattering operator. For multiple populations, a resolved operator split or channel-resolved noise representation must define each contribution; a summed tensor alone does not identify a unique per-bin transfer. The bin contributions sum to the recorded particle change. Dispersive waves require their own consistent energy/momentum closure and are not assigned a single elastic Alfvén frame by convenience.

Event deposition and an independent growth update representing the same transfer are alternatives. Include damping, cascade fluxes, escape and any plasma/background work explicitly in the appropriate energy and momentum ledgers. A wave floor cannot silently erase an exchange while retaining the particle kick. Prevent inadmissible depletion through a conservative, converged coupling procedure, or represent an external reservoir explicitly. Positivity and roundoff energy cancellation are checked separately from convergence of the growth-rate limit. A backward pseudo-particle used to reconstruct a GCR distribution is not a physical energy-carrying macroparticle and cannot feed event-wise wave growth through its unweighted histogram.

**Numerical and background validity.** Return structured failures for $B=0$, nonunit/undefined magnetic direction, incompatible particle kinematics, unresolved source conventions, invalid approximation domains and nonconverged integrals. Do not floor $B$ to make a gyrotropic model executable near a magnetic null. Report the gyroradius relative to the physical background gradient and curvature scales; the gyrotropic/guiding-centre ordering depends on these scales, not simply on the AMR cell size. Fields, spectra and their gradients come from one immutable time snapshot, with background, turbulence and topology versions. The host subdivides an interval at a published update and does not silently evaluate half of a kick using another epoch.

## 17. Library interface

Version 1.2 replaces the illustrative interface of version 1.1 with a type-separated contract. The companion implementation supplement contains `pitch_angle_contract.hpp`, which compiles in C++17. It declares interfaces and result records; it does not implement physics or bind to an uninspected AMPS checkout. Map its names and error handling to the audited host compiler and coding style at PS00. Source transcription and fixture IDs remain unchanged.

The module has six layers: read-only registry and provenance; validated SI coefficients; canonical spectra and interaction channels; stochastic/reference kernels; a transport composer for the selected equation; and AMPS cell/field-line adapters. Pure coefficient evaluation has no mesh mutation, RNG consumption, file I/O or wave deposition. The coefficient core is shared between SEP and GCR applications; their boundary data, transport equations and temporal sampling differ.

| Object | Mandatory contract |
|---|---|
| `GyroState` / `Kinematics` | Species, authoritative SI momentum, μ and local frame ID; derive $v$, $\gamma$ and both energies consistently; $p>0$ |
| `FourMomentum` | Complete vector state for an exact nonparallel boost; no substitution from scalar μ without gyrophase |
| `ScalarCoefficient` | Standard $D_{\mu\mu}$ and both analytic one-sided derivatives where finite, frame/clock, domain and kernel-pair readiness; virtual destructor |
| `TensorCoefficient` | All symmetric diffusion entries, six derivatives and PSD diagnostics; separate type from scalar coefficients; Eq. (N3) drift and correlated noise |
| `Result<T>` / `Diagnostic` | A value or a named failure with D-/U-items; unavailable data, physical refusal and nonconverged numerics have different codes; no NaN-as-success |
| `LambdaValue` | Tagged unfocused parallel mean free path, `Finite/Divergent/Unavailable/Nonconverged`, numerical error and provenance |
| `FrameContract` | Unique local frame, reference axis, clock, and retained approximation order; named source conversion applied exactly once |
| `LocalContext` / `SnapshotKey` | SI field, unit direction, signed inverse focusing length, flow, Alfvén speed and immutable background/turbulence/topology epochs |
| `SpectrumView` | Canonical one-sided total transverse power on positive $|k|$ per identified population; source quantity, units, sidedness, component content and named conversion retained; band, tail policy, variance and wave-energy convention |
| `WavePopulation` / `Exchange` | Direction relative to $\hat b$, source polarization convention, frame, population and spectral-bin IDs, event ID, spectrum epoch and plasma-frame particle energy/momentum changes |
| `StepRequest` / `StepResult` | Requested and actually consumed physical host time, coefficient-clock interval, updated full gyrotropic state, pending random path, reflections, weight increment and exchanges; zero-time success forbidden |
| `FvGenerator` | Edges, centers, stationary masses, reciprocal-integral conductances, both outgoing rates and diagonal; orientation of Q, internal extension and quadrature error stated |
| `RunRecord` | Specification/data/code digests, supported model–kernel–geometry set, equation/measure/frame/clock, decisions, RNG contract, snapshot epochs, ownership of each transport process, observable definitions and error-budget ID |

Readiness is evaluated at three points. At construction, validate the registry, units, required inputs, equation provenance and domain parameters. Before a run, require a released model–kernel–geometry combination and compatible complete transport equation. At evaluation, recheck the local domain and snapshot coverage. A formerly valid configuration may leave its physical domain during propagation; the caller then follows a configured, physically justified alternative or stops that path with its diagnostic. It cannot silently turn the failed coefficient into zero.

The admissibility rule is explicit: `Admissible` and `AdmissibleOnDomain` may run when their predicates and all other gates pass; `InadmissibleAsPrinted` cannot run as scattering; `Unestablished` is diagnostic/reproduction only. Open source audits block before configuration choices. Capability checking uses typed `ScalarDmumu` or `MuMomentumTensor`, never a string-prefix test. Readiness in the source registry means that the equation and inputs are specified; it does not establish a production-certified kernel or host adapter.

Preserve the existing source IDs and presets. Add implementation metadata through a versioned manifest or registry extension, with a documented mapping to all four admissibility values, the capability enum and domain predicate IDs. This supplement does not synthesize or certify the missing v1.1 registry. Unknown enum values or missing required metadata are schema errors. Choose any default only in an explicit application preset; a parameter-loading utility must not choose $H$, ε, a wave split, a turbulence tail or a frame convention.

Scalar kernels include Euler–Maruyama, Milstein, isotropic random rotation, the source-specific hemispheric rotation and a finite-volume-chain reference. A predictor–corrector is added only when its full stochastic calculus and source choices are resolved. Full tensor and elastic-wave kernels are separate capabilities with their own readiness. Poisson redraw remains a jump-process model and cannot advertise realization of Eq. (21). The library returns a positive timestep limit or a structured `Nonprogress` result; the host never accepts a stalled zero-duration kick.

For production performance, validate configuration and select model dispatch outside the particle loop. Avoid per-particle string parsing, registry searches, virtual ownership allocation or construction of spectrum functions. Build immutable, bounded caches only where profiling justifies them. Cache keys include model/configuration, species/energy grid, field normalization and all snapshot epochs; tabulation has a documented approximation error and returns a derivative consistent with the same interpolant. A tensor table preserves PSD through an appropriate matrix representation or convex interpolation of PSD values, rather than interpolating its entries in unrelated log spaces. Check accuracy at table knots, between knots, near singularities and across band boundaries. Do not tabulate every AMR cell by default or cache a topology pointer across regridding.

At the host boundary, keep AMPS particle weights, IDs, physical time, momentum/frame state, remaining interval and RNG/bridge state together during migration and restart. Spectral exchange event IDs prevent duplicate deposition. A counter-based per-particle RNG is an option, not an assumption about the host. If exact reproducibility across ranks is promised, define its reduction and stream ordering; otherwise require statistical invariance and deterministic conservation under rank changes. No public repository listing substitutes for inspection of the actual target checkout and its transport/pusher build paths.

## 18. Verification fixtures

**Version/evidence note.** This section records the fixtures intended for the canonical v1.1 archive. That matching archive is not available in the present review; the separate implementation supplement has its own checks/results and does not supply or certify these canonical fixtures (Section 25.7).

The fixtures are reference values against which an implementation of this specification is tested. The generator `fixtures.py` recomputes each one from the formula stated in the section named in its reference, in mpmath with 40-digit arithmetic (`mp.dps = 40`). It writes `benchmark_points.json` (both at the root of the companion bundle, Section 24), which contains 89 fixtures with 652 stored values (version 1.1 added F-FOC-04, F-FOC-05, F-MEAS-03 and F-PUSH-19 to F-PUSH-21). The run takes under a minute; without an output path it writes `benchmark_points.regenerated.json` and compares it with the shipped file, which it never overwrites unless told to. Each value is stored to 15 significant digits. Each fixture records its inputs, its outputs, the section and equation labels it checks, and its kind: "derived" (70 fixtures), "evaluation of printed formula" (15) or "published-check" (4: a value computed from printed formulas and inputs, for comparison with a number that the source prints). No stored value is a published result unless its key says "printed". The generator also checks every fixture by an independent route, for example quadrature against a closed form, a numerical derivative against an analytic one, an equivalent printed form, or a second integration variable or tail split. It prints these checks without storing them. An implementation reproduces a value to a relative tolerance of 10⁻¹² if the value is algebraic or a closed form (including hypergeometric, Dawson, Tricomi and Bernoulli forms). The tolerance is 10⁻⁹ if the stored value was obtained by quadrature or numerical differentiation, which the table marks "quadrature class" (also stated in the fixture's note). For the few exact values that are zero, the comparison is absolute. Values obtained by Monte Carlo are not fixtures. Where a Monte Carlo run of Section 12 or 16 estimates an exact quantity, that quantity is the fixture: the stationary ⟨μ⟩ (F-FOC-01, F-PUSH-08), the moments of isotropic scattering (F-PUSH-06), the running κ(t) (F-PUSH-07), the rotation-kernel moments and errors (F-PUSH-05, F-PUSH-18a), the redraw moments (F-PUSH-17) and the centred-variance diffusivity with focusing (F-FOC-04, MCR-8). The fixture IDs carry the section prefixes and are listed in section order. SIM-M5-ISO (Section 11) is F-SIM-03. F-PUSH-04 of Section 16 consists of F-PUSH-04a and F-PUSH-04b, and F-PUSH-18 consists of F-PUSH-18a to F-PUSH-18e.

| ID | Checks | Key values |
|---|---|---|
| F-FTE-01 | Parker-spiral focusing length, 1/L = ∇·b̂ = −b̂·∇ln B; He & Wan (2015) printed sign. Eq. (1) | L, θ = 90° [au]: 0.998552110968506; 1/L, θ = 90° [m⁻¹]: 6.69427969641464e-12; L, θ = 45° [au]: 0.766336968848022; He & Wan definition, θ = 90° [au]: -0.998552110968506 |
| F-FTE-02 | Zhang, Wijsen, He & Wan and EPREM coefficient sets coincide (b̂·dV/dt = 0). Eqs. (4), (6), (8), (14) | μ̇ [s⁻¹]: 0.0004315712898; μ̇, He & Wan printed L [s⁻¹]: -0.0004315287102; ṗ/p [s⁻¹]: -1.44642888888889e-06; d ln B/dt [s⁻¹]: -2.91111111111111e-06 |
| F-FTE-03 | Phase-space divergence of the GCR flow equals Ł (1 GeV proton). Eqs. (13), (18) | v [m s⁻¹]: 262326089.435118; divergence = Ł [s⁻¹]: -8.75631414342931e-12; without b̂·Du/Dt: 0 |
| F-FTE-04 | Flux-tube form: Eq. (2) for Af₀ equals A × Eq. (1); difference −μvf₀/L for f₀ | A × residual of (1): -1.32574954852767; −μvf₀/L: -0.200303446219271 |
| F-FTE-05 | ½∫(dp/dt)dμ = −(p/3)∇·u (vandenBerg2026 Eq. 3c). Eq. (17) | p [kg m s⁻¹]: 9.06411099030581e-19; ½∫ṗ dμ [kg m s⁻²]: -6.04274066020387e-25 |
| F-FTE-06 | Ruffolo convention φ = 2D_μμ (q = 1.5). Eqs. (9), (23) | λ (three routes): 2.4 |
| F-FTE-07 | Ng & Reames (2008) Eqs. 2–4 (χ = 1) equal Eq. (17); path versus fluid derivative. Eq. (20) | μ̇ [s⁻¹]: 0.0118448817226667; ε̇ = ṗ/p [s⁻¹]: -0.00026963232; change of μ̇ [s⁻¹]: 3.4452e-05 |
| F-FTE-08 | Primer momentum term without the p² Jacobian = Ł − 2ṗ/p. Eqs. (10), (13) | 2ṗ/p [s⁻¹]: -1.45502287229613e-06; divergence [s⁻¹]: 1.45501411598199e-06 |
| F-EMP-01 | PA-QFORM normalization D₀λ∥/v = (3/4)I(q,H) for 13 (q, H). Eqs. (27), (37) | q = 5/3, H = 0: 1.92857142857143; q = 5/3, H = 0.05: 1.26978038591672; q = 5/3, H = 0.2: 0.839883195924871; q = 1.67, H = 0.05: 1.27586869391809; q = 1.41, H = 0: 0.981611151102677 |
| F-EMP-02 | Closed forms of I(q,H) agree; small-H behaviour. Eqs. (37), (38) | I(5/3, 0.05): 1.69304051455562; I(5/3, 0.2): 1.11984426123316; (I(5/3,10⁻¹²) − 18/7)/10⁻⁶: -4.71238538038469 |
| F-EMP-03 | PA-EPS normalization φ(ε) = ν₀λ∥/v. Eqs. (29), (39) | φ(0.01): 7.06953301508509; φ(0.045): 4.69919053319457; φ(0.048): 4.597261773672; φ(0.1): 3.45695783731876; D_Aκ∥/v² (ε = 0.045): 0.783198422199096 |
| F-EMP-04 | Agueda et al. (2010) ν(1)λ∥/v at ε = 0.01 (printed 3.5). Eq. (39) | exact φ: 3.6054618376934; asymptotic φ: 3.52295519228089 |
| F-EMP-05 | PA-QFORM-LANG-PRINTED normalization through I⁺(q,H). Eq. (35) | I⁺(5/3, 0.05): 2.49465990914552; D₀λ∥/v: 1.87099493185914; ratio to PA-QFORM: 1.47347915640418; ratio, H = 0: 1.33333333333333 |
| F-EMP-06 | PA-EPS-PACHECO-PRINTED normalization φ_P(ε). Eqs. (36), (40) | φ_P(0.01): 5.28742168458777; φ_P/φ (0.01): 0.747916683224392; φ_P(0.048): 3.2184370720209; φ − φ_P, ε = 10⁻¹²: 1.99999999984021 |
| F-EMP-07 | PA-DROGE-VA: integral λ∥ over prefactor λ∥, 15 (q, V_A/v) and three illustrative protons. Eqs. (32), (41) | q = 5/3, a = 0.001: 0.912869772632396; q = 5/3, a = 0.1: 0.5985115566934; q = 1.61, a = 0.3: 0.499412675363873; 1 MeV proton, V_A = 100 km/s: 0.831534993748934 |
| F-EMP-08 | PA-ISO, PA-KOLMO (27/14) and PA-EPREM (factor 2) normalizations. Eqs. (33), (34) | λ_xx/λ_μμ: 1.92857142857143; 81/(7π(2π)^{2/3}): 1.08172581454567; PA-EPREM λ/λ_printed: 2 |
| F-EMP-09 | Hemispheric map, ε = 0: D_ηη = 4D_μμ = ½ν₀(1 − η²) at five μ. Eq. (31) | η (μ = 0.2): -0.6; D_ηη (μ = 0.2): 0.32; η (μ = −0.75): -0.5; D_ηη (μ = −0.75): 0.375 |
| F-EMP-10 | Wijsen (2019) Eq. 6 with P = Ck^{−q} gives Eq. 7 at H = 0. Eq. (123) | Eq. 6 = Eq. 7 (q = 5/3): 2.10387238079303; ratio Eq. 6/Eq. 7: 1 |
| F-EMP-11 | Amplitude variants: He form (v cancels); slab form. Eq. (28) | He form λ_rR^{−1/3}: 0.839883195924871; slab form λ∥ [au]: 7.39652737834412 |
| F-EMP-12 | PA-EPS derivative at μ = 0: one-sided limits ±ν₀/2. Eqs. (29), (164) | ∂D/∂μ(0+): 0.5; ∂D/∂μ(0−): -0.5 |
| F-QLT-01 | Power-law q-form amplitude. Eqs. (44), (45) | D₀/(vk_min): 0.0583279255226493; D(μ = 0.3)/(vk_min): 0.023786584848757 |
| F-QLT-02 | Piecewise TS2003 coefficient. Eq. (46) | R = 0.5, μ = 0.05: 0.0535861775221503; R = 0.5, μ = 0.0005: 0.000989539431228299; R = 3, μ = 0.9: 0.00736916795286495; N, R = 0.5: 1.00252620774442 |
| F-QLT-03 | λ∥ from Eq. (46), Q → ∞, equals MFP v1.2 Eq. 45 | Λ, R = 0.1: 28.4939411143996; Λ, R = 1: 6.13883351925882; Λ, R = 10: 6.5647150946574 |
| F-QLT-04 | Bend-over spectrum: λ∥ from Eq. (47) equals MFP v1.2 Eq. 44 | r* = 0.1: 2.39875787016309; r* = 1: 5.4871311220144; r* = 10: 58.8867317856521; C(5/6): 0.118862354635443 |
| F-DYN-01 | TS2003 I(μ), damping (DT) and random sweeping (RS). Eqs. (50), (51) (quadrature class) | I_DT (0.01, 100, 30, 0.3): 2.40475763954238e-05; I_RS (0.01, 100, 30, 0): 5.62863643137267e-06; I_DT (0.5, 50, 5, 0.05): 0.0860888169104338; I_RS (3, 300, 10, 0.7): 0.528953096891988 |
| F-WAVE-01 | Galilean plasma-to-wave-frame transformation. Eq. (61) | v′: 0.986154146165801; μ′ (V_w = +0.05): 0.253510063281697; μ′ (V_w = −0.05): 0.344447481913586 |
| F-WAVE-02 | Lorentz transformation; ΔE = γ_wV_wp′Δμ′ exactly. Eqs. (60), (62) | μ before: 0.318048556868861; E before: 1.15649041515326; ΔE: -0.00115475827773666; γ_w − 1: 5.00037503125273e-05 |
| F-WAVE-03 | Jacobian ∂μ/∂μ′ at fixed p′. Eq. (63) | μ: 0.650256088762355; ∂μ/∂μ′: 0.547321787185016; γ_w²(1 − μ²)(1 − μV_w/v)²: 0.287578693178304 |
| F-WAVE-04 | Strauss & le Roux (2019) Eq. 28 versus the exact-Jacobian δ-limit. Eqs. (56), (57) | Eq. 28, μ = 0.4: 0.0361269791952787; exact, μ = 0.4: 0.0357534624836031; ratio, μ = 0.02: 2.77515502787715 |
| F-WAVE-05 | PARALLEL Eq. 26 = Afanasiev Eq. 2 = Borovikov Eq. 5.3. Eqs. (65), (67) | D_μμ, μ = 0.4: 0.035815896660891; Afanasiev Eq. 4/Eq. 2: 1.84201574932019 |
| F-WAVE-06 | Spanier et al. (2022) R-mode root and interaction limits. Eqs. (68), (71) | y_R(3): 9.85469363727099; μ_R: 3.6697639320225; μ_L: 3.6742360679775; μ₀: 0.0101151753935137 |
| F-WAVE-07 | Schreiner & Spanier (2015) Eq. 2.11 = two-species cold plasma. Eq. (72) | n²(+), ω = 0.3: 1540.55104110358; n²(−), ω = 50: -38.7558811435497 |
| F-CMP-01 | Chandran (2000) κ∥ constants (printed c₁ ≃ 0.88). Eq. (75) | 5/2 − 3π/4: 0.143805509807655; ζ(9/2): 1.05470751076145; c₁: 0.878345252059827 |
| F-CMP-02 | Xu & Lazarian (2020) fast-mode QLT, trapping cutoff and D∥. Eqs. (78), (87) | Eq. 22, μ = 0.45: 0.0519827704407154; μ_cf,f: 0.711978783273768; D∥,f,t: 0.859486946224829 |
| F-CMP-03 | Yan & Lazarian (2008) Eq. 6: printed closed form versus its integral. Eq. (81) | printed form: 1.48224756043472; integral: 1.25302095049645 |
| F-CMP-04 | Xu & Lazarian (2018) broadened slow-mode TTD, Eq. 25 = Eq. 19 with a Gaussian. Eq. (85) | D_μμ,s1,T: 0.0577017972165909 |
| F-NL-01 | Spectrum constant D(s,q); 8π∫g dk = δB². Eq. (90) | D(5/3, 0): 0.118862354635443; D(5/3, 3): 0.222222222222222 |
| F-NL-02 | Shalchi et al. (2009) Eq. 47, D_μμl_slab/v at R = 0.0363, b = 0.05. Eq. (94) | μ = 0: 0.000229595438570042; μ = 0.05: 0.000363548890611402; μ = 0.5: 0.00133161538587155 |
| F-NL-03 | λ∥/l_slab from Eqs. 47, 49, 51 and QLT. Eqs. (95), (96), (23) (Eq. 47 value: quadrature class) | Eq. 47: 532.318586751641; Eq. 49: 2177.74361334913; Eq. 51: 553.357891639988; QLT: 684.022508440936 |
| F-NL-04 | D(μ = 0), full bend-over spectrum: step and Gaussian resonance. Eq. (97) | step, (0.0363, 0.05): 0.000229595152078543; Gauss, (0.0363, 0.05): 0.000243695823331647; step, (1, 1): 0.169330084252815; Gauss, (1, 1): 0.163207121516267 |
| F-NL-05 | Inertial-range constants: sΓ(s/2)/√π; (B₀/δB)^s coefficient of λ∥ from Eq. 47 (printed −2/3) (coefficient: quadrature class) | sΓ(s/2)/√π: 1.06141647386329; coefficient: -0.773703828796142 |
| F-NL-06 | Subashchandar et al. (2025) App. Eq. 34 at μ = 0, (M₊ + M₋)Ω²/t². Eq. (99) | Ωt = 10: 1.1455856527594; Ωt = 100: 1.01015484904774; Ωt = 1000: 0.998347116160783 |
| F-SIM-01 | Isotropic running estimators M1a, M1b over (1 − μ₀²)D. Eqs. (103), (104), (105) | M1a, μ₀ = 0.5, Dt = 0.1: 0.854891025243379; M1b, μ₀ = 0.5, Dt = 0.1: 0.72875771408333; M1a, μ₀ = 0.9, Dt = 1: 2.43492556059886 |
| F-SIM-02 | M2b exact T_μ versus printed Eq. 14, D = (1 − μ²)(\|μ\|^{2/3} + 0.05). Eq. (106) | T_exact, μ₀ = 0.5: 1.51159291940223; ratio, μ₀ = 0.1: 1.69860295890624; ratio, μ₀ = 0.5: 2.05564702996714; ratio, μ₀ = 0.9: 2.17759682206488 |
| F-SIM-03 | M5 with isotropic D_μμ (SIM-M5-ISO). Eq. (102) | Ṅ_L/Ṅ₀: 0.37428476940767; Ṅ_R/Ṅ₀: 0.62571523059233; FD₀/Ṅ₀, μ = 0: 0.817708890949728 |
| F-SIM-04 | Tautz et al. (2013a) running ⟨(Δμ)²⟩/(2τ) over QLT at τ = 100. Eqs. (107), (108) (quadrature class) | μ = 0.5, R = 0.1: 1.08841010334818; μ = 0.5, R = 1: 1.0075487362635; μ = 0.2, R = 0.1: 1.15872934147989 |
| F-FOC-01 | Focusing functions for isotropic scattering at nine ξ. Eqs. (112), (118), (110), (111), (117) | EL, ξ = 1: 0.939105856497994; HW, ξ = 1: 0.715217532132705; T12 Eq. 9, ξ = 1: 0.5625; ⟨μ⟩_stat, ξ = 1.5: 0.438124726315845; EL, ξ = 10: 0.270000001236692 |
| F-FOC-02 | Series coefficients of EL (3·2^{2n}B_{2n}/(2n)!) and HW. Eq. (114) | EL ξ¹⁰ coefficient: -6.49321284241919e-06; HW ξ⁴ coefficient: 0.161904761904762; first omitted EL term, ξ = 0.1: -6.49321284241919e-16 |
| F-FOC-03 | Tautz, Shalchi & Dosch (2014) Eqs. 12–14 from the focused orbit. Eqs. (119), (120) | ∫(1 − μ²)dt: 0.98; ∫μ̇(t)μ̇(0)dt: 0.455; λ∥^(0): 1.05 |
| F-FOC-04 | Exact long-time centred spatial variance of the constant-focusing isotropic process against the diffusion-approximation closure, six ξ including ξ = −1.5. Eqs. (115), (112) (quadrature class) | κ_var, ξ = 1: 0.248076740609470; κ∥ closure, ξ = 1: 0.313035285499331; ratio, ξ = 1.5: 0.610359686280837; κ_var, ξ = 5: 0.0103248294848381; κ_var, ξ = −1.5: 0.178275780337345 |
| F-FOC-05 | Stationary solutions of Eq. (1) with constant focusing: residual of the z-independent e^{ξμ}, the manufactured solution e^{−z/L}e^{ξμ}, forward and backward angular densities. Eqs. (1), (109), (152) (quadrature class) | residual of e^{ξμ} at μ = 0.3: 0.174275136409242; residual of the manufactured solution: 0; ρ₊(0.3): 0.557400071771623; ρ₋(0.3): 0.412932129377716; ⟨μ⟩ backward: −0.163953413738653 |
| F-MEAS-01 | Chen et al. (2024) Eq. 3 for a power-law PSD, 10 MeV proton. Eq. (122) | Ω₀ [s⁻¹]: 4.78941657797182; D, μ = 0.5, printed Ω₀ [s⁻¹]: 0.0193870952818243; D, μ = 0.5, Ω₀/γ [s⁻¹]: 0.0193187054653726; γ^{1/3}, 1 GeV: 1.27358709044433 |
| F-MEAS-02 | κ∥ integral retained by the cut μ_min = 0.05 for k^{−5/3} | retained fraction: 0.570349826471065 |
| F-MEAS-03 | SOLPACS grid: the printed index and the printed interpolant do not reproduce the second node (D-128); decreasing, increasing and flat intervals | printed q_j (decreasing): −2; printed interpolant at k_{j+1}: 4; ratio to the node: 16; decay-index convention: 0.25 |
| F-SGW-01 | Vainio & Kocharov (2001) kernel at μ = 0.6: exact over first-order moments. Eqs. (135), (169) | mean, ⟨ϑ²⟩ = 0.2: 0.967323914389265; second moment, ⟨ϑ²⟩ = 0.2: 0.978645660700087 |
| F-SGW-02 | Sandroos et al. (2014) step as extracted: mean and variance at dt = 10⁻⁶. Eq. (137) (quadrature class) | a(0.4): 0.142865056995694; mean/(a dt), η = 0.5: 0.999999761499488; mean/(a dt), η = 0.3: 1.39999966609928 |
| F-SGW-03 | Afanasiev (2015) steady state with the full resonance; SOLPACS Ω₀/V_A. Eq. (134) | λ_δ/(3C/v), C = (u₁ − V_A)(x + x₀) (δ = 0.01): 147.015; V_A [m s⁻¹]: 390862.352905799; Ω₀/V_A [m⁻¹]: 0.008332353445679 |
| F-SGW-04 | Normalization identities of the wave-frame D_μμ forms (Section 14.2) | common value: 50.3745016466493; VK01/V03, U = kW: 0.5 |
| F-GCR-PA-01 | vandenBerg2026 Eq. 13 closed form of D_A (reading check). Eq. (140) | D_Aκ∥/v², ε = 0.045: 0.783198422199096; D_A/D_I, ε = 0.045: 4.69919053319457; D_Iκ∥/v²: 0.166666666666667 |
| F-GCR-PA-02 | vandenBerg2026 at P = P_r = 1 GV, B = 4.44 nT (printed κ₀ = 2.1 au² h⁻¹). Eq. (140) | λ∥ [au]: 0.195855387282196; D_A [s⁻¹]: 0.0175320604468709; κ₀ [au² h⁻¹]: 2.09119739377523 |
| F-GCR-PA-03 | Analytic PAD (vandenBerg2026 Eqs. 47–48) and κ∥(ξ) of Effenberger & Litvinenko. Eqs. (145), (149) | A₁, ξ = 1.5: 1.31437417894754; f̃(0.3), ξ = 1.5: 0.552409471806625; κ∥/(λ₀v), ξ = 1.5: 0.29208315087723 |
| F-GCR-PA-04 | Qin & Wu (2020) D_μμ with \|μ\|. Eq. (146) | I(5/3, 0.01): 2.13599282779915; λ∥/l_slab, σ = 0.3: 283.294847460892; λ∥/l_slab, σ = 1.6: 9.95958448104698 |
| F-GCR-PA-05 | Florinski (2013): λ∥ = v/(2D), r_g, 1 MeV protons (printed 23 AU, 0.0024 AU). Eq. (148) | λ∥, D = 2 × 10⁻⁶ s⁻¹ [au]: 23.1120764197665; r_g, 4 μG [au]: 0.00241539952081483 |
| F-PUSH-01 | q-form drift and its singular term near μ = 0. Eq. (163) | D′, μ = 0.1, H = 0.05: 1.36883820162041; singular term, μ = 0.1: 1.42192689542104; D′, μ = 0.01, H = 0: 3.09315479871956 |
| F-PUSH-02 | ε-form drift; one-sided limits ±D₀. Eq. (164) | D′, μ = 10⁻¹²: 0.999999999997904; D′, μ = 0.3: 0.3712 |
| F-PUSH-03 | Milstein identity bb′ = ∂D_μμ/∂μ, b = √(2D_μμ). Eq. (159) | bb′, μ = 0.3: 0.607355341732575 |
| F-PUSH-04a | Euler–Maruyama and Milstein one-step moments, isotropic. Eqs. (157), (159) | E[Δμ]: -0.003; EM E[Δμ²]: 0.009109; Milstein E[Δμ²]: 0.0091135 |
| F-PUSH-04b | As F-PUSH-04a, q-form with focusing. Eqs. (157), (159), (163) | E[Δμ]: 0.00141758653013895; EM E[Δμ²]: 0.000754640315094283; Milstein E[Δμ²]: 0.000755079849345032 |
| F-PUSH-05 | Random-rotation kernel: exact one-step moments, effective rate. Eq. (169) | C₁(0.02): 0.990033266761799; C₂(0.02): 0.960529090939652; E[Δμ], μ = 0.3: -0.00299001997146029; E[Δμ²], μ = 0.3: 0.00899745288638967; −ln C₁/δt, νδt = 0.1: 1.01736414473968 |
| F-PUSH-06 | Isotropic moment decay, μ₀ = 0.5, νt = 1. Eq. (174) | ⟨μ⟩: 0.183939720585721; ⟨μ²⟩: 0.329184410969345 |
| F-PUSH-07 | Running κ(t), v = ν = 1. Eq. (175) | νt = 1: 0.122626480390481; νt = 10: 0.300001513330992; νt = 50: 0.326666666666667 |
| F-PUSH-08 | Constant focusing, ξ = 5: stationary ⟨μ⟩ versus the Itô drift zero. Eq. (176) | ⟨μ⟩_stat: 0.800090803982019; drift zero: 0.819803902718557 |
| F-PUSH-09 | No-boundary stationary ⟨μ²⟩ of Euler–Maruyama and Milstein (Section 16.6) | EM, νΔt = 0.1: 0.344827586206897; Milstein, νΔt = 0.1: 0.350877192982456 |
| F-PUSH-10 | Maximal per-step Euler–Maruyama exit probability; threshold √(2 − νΔt) (Section 16.4) | threshold, νΔt = 0.1: 1.37840487520902; P_max, νΔt = 0.1: 0.0840391595174851; 1 − Φ(√2): 0.0786496035251426 |
| F-PUSH-11 | Krülls–Achterberg scale 2D_μμ/(∂D_μμ/∂μ)² near μ = 0. Eqs. (163), (164) | q-form H = 0.05, μ = 0.001: 0.00270009990282157; q-form H = 0, μ = 0.001: 0.000450003150018; ε-form, μ = 0.001: 0.0984100064104579 |
| F-PUSH-12 | Hemispheric η map, ε = 0. Eq. (170) | η (μ = 0.7): 0.4; 4D_μμ (μ = 0.7): 0.42 |
| F-PUSH-13 | Importance-sampling identity L[f]/(1 + μ/a) = L′[u]. Eqs. (153), (154) (quadrature class) | μ = 0.1: 1.90031420310832; μ = 0.85: 3.82402885342927 |
| F-PUSH-14 | Conservative focusing for F = e^{z/L}f₀ equals e^{z/L} × advective form. Eq. (150) (quadrature class) | both sides: -0.185128772105901 |
| F-PUSH-15 | Lorentz transformation of (v, μ). Eq. (172) | γ′: 1.04803534443389; v′: 0.29927693604669; μ′: 0.394378291750548 |
| F-PUSH-16 | Crossing of μ = 0: ∫dμ/D_μμ for \|μ\|^α(1 − μ²) (Section 16.5) | α = 2/3, δ = 0: 2.48001376922875; α = 1, δ = 10⁻¹²: 27.0817149715945 |
| F-PUSH-17 | Poisson-timed isotropic redraws: moments and κ. Eq. (171) | ⟨μ²⟩(1): 0.302676713235713; κ(T = 10): 0.300001513330992 |
| F-PUSH-18a | Rotation kernel: exact errors at νt = 1 from μ₀ = 0.5 (Sections 16.6–16.7) | err ⟨μ⟩, νδt = 0.1: -0.00316638558422697; err ⟨μ²⟩, νδt = 0.1: 0.000210281242324494; err ⟨μ⟩, νδt = 0.01: -0.000307540507421492 |
| F-PUSH-18b | Addition theorem E[μ′²] = 1/3 + (μ² − 1/3)⟨P₂(cos ϑ)⟩. Eq. (169) | E[μ′²], μ = 0.3: 0.0972034409035135; E[Δμ²]/(2δt), μ = 0.3: 0.449872644319483 |
| F-PUSH-18c | Euler–Maruyama without boundary: error of ⟨μ⟩ (Section 16.6) | νΔt = 0.1: -0.00960050053572116; νΔt = 0.01: -0.000923549949106408 |
| F-PUSH-18d | μ = ±1 not reached: ∫_{0.5}^{1−δ}dμ/(1 − μ²) diverges (Section 16.4) | δ = 10⁻⁴: 4.40241263130899; δ = 10⁻¹²: 13.6127780039099 |
| F-PUSH-18e | Product rule ∂_μ(D∂_μf) = D′f′ + Df″ used by the backward SDE. Eq. (152) (quadrature class) | value at μ = 0.4: -2.78272515195437 |
| F-PUSH-19 | Step rule Eq. (165) at μ = ±1 for the ε-form; clock factor of Eq. (173); SI kinematics of a 100 MeV proton (Section 1.6) | Δτ at μ = ±1: 0; second entry at μ = 1: 0.000223982224770642; dt′/dt (u = 400 km s⁻¹, v = c/2, μ = 0.7): 0.999533899972278; rigidity: 444583420.329749 V; Ω at 5 nT: 0.432812934579193 s⁻¹ |
| F-PUSH-20 | Legendre-moment decay under ν(t) = 1 + 0.5 sin 3t to t = 2. Eq. (177) | ∫ν dt: 2.00663828555827; l = 1: 0.134439864291891; l = 2: 0.00242987647398007; l = 3: 5.90429967880181e-06 |
| F-PUSH-21 | Transmissibility across μ = 0 of the finite-volume reference operator for the q-form, H = 0, |μ| ≤ 0.05. Eq. (178) (quadrature class) | ∫dμ/D: 2.211209399028; T: 0.452241203587312 |

## 19. Discrepancy, decision and gap register

This section collects four lists: printed discrepancies and reading conflicts found in the sources (19.1, IDs D-n); decisions that the sources leave to the user or the host code (19.2, IDs U-n); information that was sought and not found (19.3); and findings that concern the companion specifications MFP v1.2 and PARALLEL rev. 1.4 (19.4, IDs C-n). Each item is recorded as printed or as read and is not resolved by assumption. Identifiers of the companion registers always carry their prefix ("MFP v1.2 D-30", "MFP v1.2 U-8"); an entry here that repeats an MFP v1.2 item says so. "Where" lists the sections of this document that state or use the item. D-1 to D-119, U-1 to U-65 and C-1 to C-15 are ordered by the first of them; D-120 to D-127 were added by the final source and review checks (Section 25.3); D-128, U-66 and U-67 were added in version 1.1 after the external review (Section 25.6), which also settled U-25 and U-44 and extended U-2 and U-55. "Handling" restates what those sections do. Fixture labels F-… refer to Section 18; check labels (C1–C14, S8-…, S9-…) refer to the derivation records DR-02 to DR-16 of the companion bundle (Section 24).

### 19.1 Printed discrepancies and inconsistencies

| ID | Item | Where | Handling |
|---|---|---|---|
| D-1 | He & Wan 2015 Eq. 3 (two renderings) prints dμ/dt = (1−μ²)/2·[−(v/B)∂B/∂z + …] = (1−μ²)/2·[v/L + …] with L = (ẑ·∇ln B)⁻¹. With this L, v/L < 0 where B decreases along +z, so the second form has the opposite sign to the first, which agrees with every other source (including He et al. 2011 and Tenishev et al. 2022, L = B/(−∂B/∂z)). He & Wan 2019 and He et al. 2017 print the v/L form with the same L, and the backward SDE of Section 16.2 uses V/L. At 1 au in a Parker field the printed definition gives L = −0.998552 au (F-FTE-01). | 1.5, 2.1, 2.2, 2.7, 16.2, 16.10 | Coefficient set implemented with the first printed form, −v b̂·∇ln B; the v/L form with the printed L only as variant FTE-SW-HEWAN15-PRINTED (U-1). |
| D-2 | van den Berg et al. 2020 Eq. 9, Lang et al. 2024 Eq. 1 and Merten & Aerdker 2024 Eq. 1 print the conservative 1-D form ∂_s[μvf] + ∂_μ[(1−μ²)vf/(2L)]. It equals Eq. (1) only if f is the flux-tube (linear) density Af; for the phase-space density the two differ by μvf/L (derived, F-FTE-04). None of the three states which density is meant. | 2.1 | Density is a configuration item (U-2). |
| D-3 | Afanasiev et al. 2025 Eq. 40 (arXiv v1): the acceleration term of dp/dt is read as μ/ν (Greek nu) in one reading and μ/v in another. Wijsen et al. 2019 Eq. 4 and Afanasiev et al. 2025 Eq. 39 use the particle speed v. | 2.3 | Read as μ/v; recorded as printed. |
| D-4 | van den Berg et al. 2020 Eq. 17 prints ∂/∂p[(dp/dt)f]; van den Berg et al. 2026 Eq. 1 prints (1/p²)∂/∂p[p²(dp/dt)f] with f "per unit momentum volume". For f per unit momentum volume the two differ by 2(ṗ/p)f (derived, F-FTE-08). The 2020 normalization of f is not stated in the portion read. | 2.3 | Eq. 17 is not used with f per d³p without the Jacobian (U-2). |
| D-5 | Wijsen et al. 2019 Eq. 1 (advective, for f) and Afanasiev et al. 2025 Eq. 37 (conservative, for j = p²f) use the same coefficients. They are not algebraically identical when b̂·dV_sw/dt ≠ 0: they differ by −(μv/c²)(b̂·dV_sw/dt)j (derived, F-FTE-03). | 2.3 | Recorded; each preset implements one stated form. |
| D-6 | Linker et al. 2019 print the focusing term as −v e_b∇ln B (no dot); Young & Poduval 2025 Eq. A1 print −v b̂·∇ln B. | 2.4 | Read as e_b·∇ln B. |
| D-7 | EPREM: Kozarev et al. 2013 Eq. 1 and Linker et al. 2019 multiply df/dt by (1 − V·e_b vμ/c²); Young & Poduval 2025 Eq. A1 do not. Same code, different printings. | 2.4 | Variant flag on the EPREM coefficient set (U-3). |
| D-8 | EPREM operator ∂_μ[(D_μμ/2)∂_μf] with D_μμ = (1−μ²)v/(2λ∥) (Kozarev et al. 2013 Eqs. 1–2; Young & Poduval 2025 Eq. A1; Schwadron et al. 2010): Eq. (23) returns 2λ∥ (derived). Same item as MFP v1.2 D-2. | 2.4, 5.4, 16.1 | Operator factor recorded per coefficient set; the EPREM-compatible mode is MFP v1.2 U-8. |
| D-9 | Tenishev et al. 2022 Eqs. 1–2: Eq. 1 contains V^sw·∇f, transport is "in Lagrangian coordinates", and Eq. 2 is dz = μv dt. Whether the field-line motion replaces the V^sw·∇f term, and along which axes ∂V_x/∂x, ∂V_y/∂y, ∂V_z/∂z are taken, is not stated. | 2.4 | U-6. |
| D-10 | Tenishev et al. 2022 (arXiv v1) Eq. 4: dμ = [(1−μ²)/2·(v/L + μ(…))]dt + √(2D_μμ)dW_μ(dt). The drift lacks the ∂D_μμ/∂μ that the Itô SDE of their Eq. 1 operator needs (Eq. (25)). The text after Eq. 4 gives "dW_μ(dt) = dt · Λ", Λ Gaussian; a Wiener increment scales as √dt. | 2.4, 5.4, 16.10, 16.11 | Not implemented as printed; recorded as a printed variant only. Whether the code contains the drift is U-65. |
| D-11 | van den Berg et al. 2026 Eq. 36 (time-backward form; raw PDF text layer, one reading): the term (μv/c²)b̂·[∂u/∂t + (u·∇)u]f appears with the sign opposite to that implied by their Eqs. 15, 18 and 37 and by the derived phase-space divergence (F-FTE-03). | 2.5 | Not resolved; Eqs. 15, 18 and 37 are followed (Section 15.2); Eq. 36 is to be re-read on the typeset page. |
| D-12 | Ng & Reames 2008: the 2×2 diffusion matrix of Eq. 1 and the placement of c, e³, g_σ and V_A in the growth rate, Eq. 7, are not legible (NASA preprint, two fetches; published PDF garbled). | 2.6, 8.2, 8.5, 14.3, 14.6 | The Ng–Reames model is not implementable from the texts read; WAVE-NR08 not fully specified. |
| D-13 | μ^{s−1} printed without absolute value: Qin & Wang 2015 Eq. 2; Wang et al. 2014 Eq. 3; Qin & Wu 2020 Eq. 15 (PDF v2) / Eq. 17 (ar5iv; LaTeX \mu^{s-1}); Wu & Qin 2020 Eq. 10 / 12. For non-integer s a floating-point power of μ < 0 is undefined; for s = 5/3 the real reading (μ²)^{1/3} equals \|μ\|^{2/3}. | 5.2, 6.3, 15.3, 15.9 | Printed form only as variant (GCR-PA-QINWU20-PRINTED); \|μ\|^{s−1} only when the configuration says so (U-7). |
| D-14 | Wang et al. 2014 Eq. 3 and AMPS Type V (Tenishev et al. 2022, Table 1) print R^{s−2} with R a Larmor radius, which is dimensionally consistent only if R means R_Lk_min. Same item as MFP v1.2 D-6. | 5.2, 5.4, 6.3 | Qin & Wang 2015 form used; it equals the TS2003 power-law form of Section 6.3 (derived). |
| D-15 | k_min = 33 AU⁻¹ (Qin et al. 2013, partly garbled; AMPS Type V, exponent corrupted in the extraction) versus 32 AU⁻¹ (Qin & Wang 2015). Different papers of the same group; not an error of one paper. | 5.2, 5.4, 6.3 | Both stored as separate presets. |
| D-16 | Zhang et al. 2023 Eq. 35 prints D₀(x)p^{q−2}(1−μ²)(\|μ\|^{q−1}+h₀) with no factor v, and Eq. 37 prints λ∥ = λ∥0(x)p^{2−q} (arXiv PDF, two fetches). Eq. (23) gives λ∥ = [3vI(q,h₀)/(4D₀(x))]p^{2−q}, so λ∥0 must contain v (or D₀ must be ∝ v); the paper does not say which. | 5.2 | Recorded; amplitude mode "D₀(x) fixed from λ_r at each (x, p)" (U-8); companion note C-2. |
| D-17 | Wijsen et al. 2019 (arXiv and publisher): ε = 0.048 is "obtained from minimizing the sum of the squared differences between (8) and (7) for H = 0.05", but q is never given a numerical value, so the ε–(q, H) correspondence cannot be reproduced. | 5.2, 5.3, 13.3 | q stored as "not printed" (U-12). |
| D-18 | Laitinen et al. 2016 Eq. 6: q is "the spectral index of the slab turbulence", with no number; the slope 5/3 is only implied by Eq. 10. | 5.2 | q = 5/3 stored with the flag "implied by Eq. 10". |
| D-19 | Qin et al. 2005 Eq. 6 (garbled repository PDF): the grouping of h with \|μ\|^{q−1} is not legible, and no rigidity factor is legible. | 5.2 | Foundation only; later printings used. |
| D-20 | Àgueda & Vainio 2013: q = 1.66 (Sect. 3.1, H = 0, ε = 0.01, case A) and q = 1.67 (Sect. 3.2, H = 0.05, ε = 0.045) for nearly the same case. | 5.2 | Both values stored with their ε. |
| D-21 | Units not stated: R in He et al. 2011 Eq. 6, He & Wan 2015 Eq. 7 and He 2015 Eq. 4 (D₀vR^{−1/3}); p in Wang & Qin 2015 Eq. 3 and Kubo et al. 2015 Eq. 6 (D₀vp^{q−2}). D₀ values are not portable without them. | 5.2 | Units stated by the configuration (U-9). |
| D-22 | Kubo et al. 2015 Eqs. 6–8: D₀vp^{q−2}(…) with λ_rr = λ∥cos²Φ "constant across the entire interplanetary space". With a constant D₀, λ∥ is independent of position, so λ_rr cannot be constant; D₀ must vary as cos²Φ, and the procedure is not printed. | 5.2 (amplitude table), 15.5 | Recorded in 15.5; amplitude mode "from λ_r with ψ". |
| D-23 | Àgueda et al. 2010, text after Eq. 5: "ν(μ = 1) = 3.5υ/λ" for ε = 0.01. Since ν(1)λ∥/v = φ(ε)(½ + ε), this is 3.523 with the asymptotic φ = (3/2)ln(1/ε) of their Eq. 4 and 3.605 with the exact φ (derived, F-EMP-04). | 5.3 | Kept as printed; the library uses the exact φ. |
| D-24 | Afanasiev et al. 2025 (arXiv numbering): q₀ denotes both the Kolmogorov index in λ⁰ and the exponent of the foreshock fit (Eq. 21). The ε-form of their Eq. 41 (also Wijsen et al. 2022) has no rigidity factor, unlike Wijsen et al. 2019 Eq. 8; the two coincide only when the amplitude is fixed from λ∥ at each energy (derived). | 5.3 | Two amplitude modes (PA-EPS, PA-EPS-W19). |
| D-25 | Minoshima et al. 2026 Eq. 4: the pitch-angle term is ∂_μD_μμ∂_μ[(1 − μvV^c_sw/c²)f]; the operator acts on (1 − μvV/c²)f, not on f. D_μμ of their Eq. 8 enters with this operator. | 5.4 | PA-DROGE-VA carries the operator flag; standard for V = 0. |
| D-26 | Minoshima et al. 2026 Eq. 8 (Dröge V_A form): λ∥ appears in the prefactor, and no integral definition is printed. For V_A > 0, Eq. (23) returns a smaller value (ratio 0.913 at V_A/v = 10⁻³, q = 5/3; derived). | 5.4, 5.6 | Meaning of the parameter is U-14; MFP v1.2 Section 6.5 already requires it to be reported. |
| D-27 | Minoshima et al. 2026 Eq. 9 prints λ∥ = ξv/Ω_n; Ω_n is not defined (the text gives ν_p = Ω_p/2π). Same item as MFP v1.2 D-29. | 5.4 | The gyrofrequency is stated by the configuration. |
| D-28 | Borovikov et al. 2019 state k₀⁻¹ = L_max/2π and evaluate 81/(7π(2π)^{2/3}) "≈ 0.92"; the expression equals 1.0817, so the printed coefficients 0.9 (Eq. 6.10) and 0.5 (Eq. 6.11; AMPS Type I) are inconsistent with the stated convention. Same item as MFP v1.2 D-7. | 5.4 | Exact expressions used. |
| D-29 | AMPS Table 1 (Tenishev et al. 2022, arXiv only, partly garbled extraction): MFP v1.2 Section 6.7 transcribes μ^{2/3}, μ^{q−1}, μ^{s−1}; the readings read \|μ\| with bars "stripped in places". The two readings of the same table disagree. | 5.4, 6.3 | Not resolved; \|μ\| in Section 5.4 is marked not verified; the typeset table is to be read before PA-AMPS-I, IV and V are implemented (companion note C-3). |
| D-30 | AMPS Type II (Tenishev et al. 2022, Table 1): D₀[1/(1+ε)](1 − μV_A/v)²(…) + D₀[ε/(1+ε)](1 + μV_A/v)²(…). One reading gives the factors (…) as \|vμ ∓ V_A\|^{2/3}/(\|vμ ∓ V_A\|^{5/3} + (Ωl_b)^{5/3}); the other found them not legible. | 5.4 | PA-AMPS-II not implemented until read on the typeset page (U-15). |
| D-31 | AMPS Type IV (Tenishev et al. 2022, Table 1) is attributed to "Wang and Qin, 2004"; the matching paper found is Wang & Qin 2015 (arXiv 2013). | 5.4 | Year as printed. |
| D-32 | Lang et al. 2024 Eq. 3 and Lavasa et al. 2026 Eq. 3 print D₀(\|μ\|^{q−1}+H)(1−μ)²; every other printing of this shape has (1−μ²). Taken literally (derived): D_μμ(μ)/D_μμ(−μ) = [(1−μ)/(1+μ)]², D_μμ(−1) = 4D₀(1+H) ≠ 0, and at q = 5/3, H = 0.05 the same λ∥ needs a D₀ 1.473479 times larger. Same item as MFP v1.2 D-3. | 5.5, 5.6, 16.11 | Standard (1−μ²) form implemented; printed form only as variant PA-QFORM-LANG-PRINTED with an explicit no-flux condition at μ = −1 (MFP v1.2 U-5). |
| D-33 | Pacheco et al. 2019 Eq. 4 print \|μ\|/(1−\|μ\|) in the ε-form; Àgueda & Vainio 2013, which they cite, print \|μ\|/(1+\|μ\|). Taken literally (derived): D_μμ(±1) = ν₀, D_μμ increases monotonically with \|μ\| for ε < 1, and φ_P(0.01) = 5.287422 against φ(0.01) = 7.069533. Same item as MFP v1.2 D-4. | 5.5, 5.6 | Standard \|μ\|/(1+\|μ\|) implemented; printed form only as variant PA-EPS-PACHECO-PRINTED (U-16). |
| D-34 | TS2002 Eqs. 32–35, 44: σ in P_RR = (1−σ)g, P_LL = (1+σ)g is never defined; Eq. 35 is printed without σ, and Eq. 44 states D_μμ(−μ) = D_μμ(μ). This is consistent only if σ does not depend on k∥. Hermiticity of a real field requires σ(−k∥) = −σ(k∥) for the printed tensor form, and with σ = σ_m sgn(k∥) the magnetostatic D_μμ is multiplied by 1 − σ_m sgn(Ωμ) and is not even (derived, DR-06 C1, C14). | 6.2, 7.1 | Slab models specified for σ = 0; helicity only after U-17. |
| D-35 | Reading disagreements in the readings: TS2003 Eq. 5 integrand (one reading: possibly k g(k)); a possible leading minus in Shalchi et al. 2006 Eq. 20; charge-sign independence stated in TS2003 only. The second reading has g(k∥) (two renderings), finds no minus sign, and finds the statement also in TS2002. | Sections 6–7 (readings used) | Second reading used. |
| D-36 | Ding et al. 2022 Eq. 8 prints B²vμ in the denominator (two renderings; a further single-fetch reading did not establish whether \|μ\| is printed), while the same sentence gives k∥ = Ω(v\|μ\|)⁻¹. As printed, D_μμ < 0 for μ < 0 (derived). | 3.1, 6.7, 13.6 | Read as v\|μ\| (U-19). |
| D-37 | Afanasiev et al. 2015 Eq. 2: the denominator is B² in the ar5iv rendering of the TeX source (one reading also reads B² in the A&A journal HTML) and R² in the PDF text extraction; the text defines B as the mean field. MFP v1.2 D-30 records R² "as printed in the arXiv version". | 6.7, 8.3, 13.7 | The B² reading is used in the relations of Sections 6.7 and 13.7; the choice is U-20; companion note C-11. |
| D-38 | Laitinen et al. 2016: δB² = 2∫S dk is printed without integration limits; Eq. 2 evaluates S∥ at −(r_Lμ)⁻¹, whereas Eq. 11 defines S∥ for k∥ > 0. The coefficient differs by a factor 2 between the two readings of the limits (derived, DR-06 C10). | 6.7, 13.4, 13.9 | U-21. |
| D-39 | Engelbrecht et al. 2022: S_xx in Eq. 11 is "the power spectrum of b_x" and is not normalized in the text. The DT prefactor of their Eq. 39 is 3s/(√π(s−1)) as extracted, whereas TS2003 Eq. 24 (DT) has no √π; the extraction of Eq. 39 is flagged ambiguous. | 6.7, 13.9 | Recorded; the λ∥ forms are in MFP v1.2 Section 7.2; the S_xx convention is U-22. |
| D-40 | Chen et al. 2024 Eq. 3: P(f_res) is computed from B_N; its sidedness, units and estimator are not stated. "One-sided PSD of B_N" and "one-sided total transverse PSD" differ by a factor 2 under axisymmetry (derived, DR-06 C10, C11). | 6.7, 13.2, 13.9 | U-23. |
| D-41 | Zhang et al. 2023 Eq. 36 is not legible. One extraction gives the raw string "π / 2Ω2 (1 − μ 2) / B2v\|μ\| W⊥(k_res)", whose fraction grouping cannot be determined; a second fetch gives π²Ω²(1−μ²)W⊥(k_res)/(B²v\|μ\|). The candidate readings π²Ω² and (π/2)Ω² differ by 2π. The argument of W⊥ is printed as α⊥. | 6.7, 13.5, 13.9 | MEAS-ZHANG23 not fully specified (U-24); companion note C-4. |
| D-42 | TS2003 Eq. 11 (dimensionless RS form) prints √π(s−1)v(1−μ²)b/(4k_minR_L²), without s (ALT text and PDF agree), while Eqs. 9, 10, 13, 17 and 27 carry s. Eq. 4 with the normalization of Eq. 5 gives s; without it the coefficient is larger by the factor s (derived, DR-06 C6). | 7.2 | TS2003 Eq. 17 used; the prefactor of Eq. 11 is not used. |
| D-43 | TS2002 Eqs. 58, 62 and TS2003 Eq. 24 (λ/λ₀): one reading of the notes has a (or b) in the numerator; the second reading (complete ALT text) has it in the denominator together with k_min, which follows from Eqs. 57/23 and 55 (derived). | 7.3 | Second reading used. |
| D-44 | Strauss et al. 2017 Eqs. 14, 16, 20 print the Lorentzian as γ/[γ + (…)²] (ar5iv and PDF for Eqs. 14, 16; PDF only for Eq. 20). γ is a rate and (…)² a squared rate, so the kernel is dimensionally inconsistent; ∫γ/(γ + x²)dx = π√γ → 0, which contradicts the printed limit Eq. 17. With γ² the form equals TS2002 Eq. 35 exactly (derived, DR-06 C8). | 7.4, 8.1 | DYN-TS-DT is the default; DYN-STRAUSS17-PRINTED only as an explicit variant (U-25). |
| D-45 | Shalchi, Lazarian & Schlickeiser 2007 Eq. 20 prints Σ_{n=+1}, whereas Eqs. 17, 21 and 24 use n = ±1. n = +1 alone gives 1.24 % of the n = ±1 value at the test point (DR-06 C9). | 7.5, 10.5 | n = ±1 used. |
| D-46 | Shalchi et al. 2006 Eq. 7 (OCR): e^{+iωt} (no minus glyph visible) versus e^{−iv_Ak∥t} for the same model in Table 1. | 7.6 | NADT not implemented (DYN-NADT-SH06 not fully specified). |
| D-47 | Strauss & le Roux 2019, parameter list: T_p⁰ printed as "30^5 K", which is not a valid magnitude. | 7.7 | Data row flagged; no default. |
| D-48 | Strauss & le Roux 2019 Eqs. 27–29: D^ND ∝ 1/(v\|μ\|) with k^res = nΩ/(V_A − vμj). The δ-limit of their Eq. 25 gives 1/\|vμj − V_A\|; the printed form keeps V_A in k^res but drops it from the Jacobian. For V_A/v = 0.05 the two differ by 1.04 % at μ = 0.4 and by a factor 2.78 at μ = 0.02 (S8-2). | 8.1, 8.5 | Implemented as printed (WAVE-SL19-ND-PRINTED); the exact-Jacobian δ-limit is a separate derived variant (U-27). |
| D-49 | Strauss & le Roux 2019 Eq. 28 prints the indices as g^{n,j}_{slab}; elsewhere g^{slab}_{n,j}. | 8.1 | Notation only. |
| D-50 | Strauss & le Roux 2019, frame of μ: one reading records "μ is taken in the flow frame (Sec. 5)"; another found no explicit statement. | 8.1, 8.5 | Both readings stated; Sect. 5 of the paper is to be re-read. |
| D-51 | Virtanen & Vainio 2005, App. A (pre-2006; extraction): "Δt = a ν_e(γ)", a ≲ 0.01. Dimensionally Δt = a/ν_e is expected. | 8.2, 14.6, 16.6, 16.10 | As printed; extraction-level; not used. |
| D-52 | Battarbee et al. 2011 Eqs. 1–2: ν = π²f_ci f_res P(f_res)/B² with f_ci = qB/(2πm_ic), no γ. With Ω = 2πf_ci/γ the Vainio, Laitinen & Fichtner 2003 form is smaller by γ, and Battarbee et al. Eq. 8 uses ω_ci = Qω_cp/(Aγ). As printed, Eq. 2 exceeds the Vainio, Laitinen & Fichtner value by γ at \|μ\| = 1 (S8-7). | 8.3 | U-31. |
| D-53 | Resonance sign: k_r = −ω_c/(vμ) with δ(k + ω_c/(vμ)) in Vainio 2003 Eqs. 2–3 and Battarbee et al. 2011 Eq. 10; k_r = +Ω₀/(γvμ) with δ(k − k_r) in Afanasiev et al. 2015 Eq. 1, 2023 Eq. 5 and 2025 Eqs. 15–16. Opposite conventions for k, not an error within either paper. | 8.3, 14.2, 14.8 | The sign convention is stored with any signed-k spectrum (U-32). |
| D-54 | Spanier et al. 2022 Eq. 22 (Full text ×2) prints δ(vk∥ − ω_R − Ω′) and δ(vk∥ − ω_L + Ω′) without μ, whereas Eq. 18 has vμk∥ and the dimensionless Eq. 24 has kμ/ε. | 8.4 | Eqs. 24 and 38 implemented; Eq. 22 flagged as printed. |
| D-55 | Spanier et al. 2022, unnumbered helicity-resolved form (PDF only): the L-terms contain δ(±kμ/ε + y_R(k) − S_i), where y_L is expected by analogy with Eq. 24. This may be an extraction error. | 8.4 | Not implemented; Eqs. 24–28 and 38 used. |
| D-56 | Spanier et al. 2022, author order: Spanier, Schreiner, Schlickeiser in the library index; Schlickeiser, Spanier, Schreiner in the MDPI citation block. | references (source of 8.1, 8.4) | To be checked with Crossref; reference record PARTIAL. |
| D-57 | Chandran 2000, after Eq. 9: the normalization integral runs from l⁻¹, whereas the validity statement has (l′)⁻¹ < k⊥. | 9.1 | As printed; the normalization evaluates exactly (S9-1). |
| D-58 | Yan & Lazarian 2002 Eq. 8 has the prefactor 1/(2Ω^{1.5}L^{2.5}); Yan & Lazarian 2004 Eq. 12 has 1/(Ω^{1.5}L^{2.5}) for the same quantity. | 9.1, 9.6 | Two IDs (CMP-YL02-ALF, CMP-YL04-ALF-PRINTED). |
| D-59 | Xu & Lazarian 2018 Eq. 13 and 2020 Eq. 13: the slow-mode spectrum is printed with exp(−L^{1/3}k∥/k∥^{2/3}); k⊥^{2/3} is expected by analogy with Xu & Lazarian 2020 Eq. 5. | 9.1 | As printed. |
| D-60 | Xu & Lazarian 2018 Eq. 14 and 2020 Eq. 6 versus 2018 Eq. 30: C = δB²L^{−1/3}/(6π) is "normalized to δB²/2", and C_f = δB_f²L^{−1/2}/(16π). With k⊥ ≥ L⁻¹, the Alfvén/slow C gives δB²/2 only over k∥ ≥ 0 (δB² over all k∥); the fast-mode C_f gives δB_f²/2 over the full sphere (derived, S9-6). | 9.1, 9.2, 9.6 | The k∥ domain is U-34. |
| D-61 | Yan & Lazarian 2002 Eq. 1: ⟨v_iv_i*⟩/V_A² = δK_ij, with the index i repeated. | 9.1 (source) | Notation only. |
| D-62 | High-β fast-mode spectrum prefactor: L^{−1/2}/(8π) in Yan & Lazarian 2002 Eq. 5, 1/(2π) in Yan & Lazarian 2004 Eq. 10; a factor 4. | 9.2 | Both recorded. |
| D-63 | Yan & Lazarian 2004 Eq. 3: V_ph is "equal to Alfven and sound velocity for high and low beta plasma, respectively", inverted relative to the rest of the paper. | 9.2 | As printed. |
| D-64 | Yan & Lazarian 2004 Eq. 11: the D_pp row is duplicated (in the print or in the extraction). | 9.2 (source) | Not implemented. |
| D-65 | Yan & Lazarian 2008 Eq. 6 (Full text ×2): Re∫₀^∞dt e^{iω̃t − ½k∥²Δv∥²t²} is printed as (√π/\|k∥Δv∥\|)exp[−ω̃²/(k∥²Δv∥²)]. The integral equals (π/2)^{1/2}/(\|k∥\|Δv∥)·exp[−ω̃²/(2k∥²Δv∥²)]; the printed closed form is narrower by √2 (both have area π). At k∥²Δv∥² = 0.7, ω̃ = 0.5: printed 1.48225, integral 1.25302 (S9-8). | 9.3 | Closed form implemented as printed under CMP-YL08-RES with the convention recorded (U-36); companion note C-15. |
| D-66 | Yan & Lazarian 2008 Eq. 9 (small-rigidity TTD): ar5iv bracket [−E₁(qx∥) − e^{−qx∥}]₁^{x∥,max}; PDF OCR bracket [−E₁(qx_k) − e^{−qx_k}/x_{k,max}]. A small-R reduction of their Eq. 8 reproduces the ar5iv bracket, q, M_A^{7/2} and (1 − μ²)^{3/2}, with prefactor √π/8 (printed √π/16) and an extra x⊥ = 1 boundary term (derived, S9-10). | 9.3, 9.6, 10.5 | CMP-YL08-TTD-SMALL blocked until the typeset page is read (U-37). |
| D-67 | Yan & Lazarian 2008, after Eq. 5: "R = v/(Ωl) is the dimensionless rigidity of the CRs, L is the injection scale" (l versus L). | 9.3 | R = v/(ΩL) used (S9-9 reproduces Eq. 8 with L). |
| D-68 | Yan & Lazarian 2008: WIM β = 0.077 in Table 1 and β ≃ 0.11 in Sect. 4.2. | 9.3, 10.5 | Both recorded; the companion-data row is flagged. |
| D-69 | Yan & Lazarian 2008 Eqs. 14–16: Δμ_k (Eqs. 14–15) and r (Eq. 16) are not defined in the extracted text. | 9.3 | Not implemented. |
| D-70 | Xu & Lazarian 2020 Eq. 46: B_s² in a denominator is not a defined symbol. | 9.4 | As printed. |
| D-71 | Lynn et al. 2012 Eq. 13 (calibrated D∥): 1/u₊³ in both terms; u₋ is expected in the second. | 9.5 | As printed. |
| D-72 | Shalchi et al. 2009 Eq. 21 (PDF raw layer): k₀ appears without a definition in the lines read. | 10.2 | Eq. 24 used instead. |
| D-73 | Subashchandar et al. 2025 App. Eq. 35 (both renderings): lim_{t≫1/Ω, μ=0} M₊ + M₋ = 1/(Ωt)². The dominant term of App. Eq. 34 is t²/Ω², and with 8π∫g dk = δB² only t²/Ω² gives App. Eq. 36, σ_z² = (v²/2)(δB/B₀)²t² (derived, F-NL-06). | 10.4 | App. Eqs. 32–34 are evaluated directly; Eq. 35 is not used. |
| D-74 | Subashchandar et al. 2025 App. Eqs. 12–13: the signs of the gyration terms are flagged by the extractor as inconsistent with Eq. 14. | 10.4 | Not used (Eq. 15 for μ̇ is consistent). |
| D-75 | Subashchandar et al. 2025, text near Eq. 4: the indices are written "−q = 1.5" and "−5/3", with P_S = 2C_Sk∥^{−q^S}, while q^S = 1.5 and q^{2D} = 5/3 are used as positive indices. | 10.4, 13.8 | q^S = 3/2, q^{2D} = 5/3 as positive indices of k^{−q}, as in Eq. 4; to be confirmed by the user (U-40). |
| D-76 | Subashchandar et al. 2025: κ⊥/κ∥ is given as "10⁻⁴–10⁻³", as "three orders" and as "three to four orders" in different places. | 10.4 | Not used numerically. |
| D-77 | Tautz et al. 2013a, App. A (A.2, A.4, A.6b): garbled; ξ = v²(δB/B₀)^{2/3} "as extracted" is dimensionally suspect and was reconstructed by the extractor. | 10.5 | Not used; only A.3. |
| D-78 | Weidl et al. 2015 (extraction): "μ = cos(v_z/v)"; presumably μ = v_z/v. | 11.1, 11.5 (source) | Extraction error assumed; flagged. |
| D-79 | Tautz 2013b Eq. 28 is present in the arXiv PDF only; the A&A HTML numbering ends at Eq. 27. | 11.1 | Cited as arXiv-only. |
| D-80 | Riordan & Pe'er 2019 Eq. 11 (one rendering): D_μμ,MSD = ⟨(Δμ)²/Δt⟩, without the ½ of Eq. (22); it is twice that normalization. | 3.2, 11.1, 11.7 | Converted explicitly; never mixed silently with Eq. (22) (U-41). |
| D-81 | Sun et al. 2016: van den Berg et al. 2024 state that the M5 formula of Sun et al. contains a sign typo; the formula was not read (secondary). | 11.1 | The M5 form of Section 11.1, whose signs agree with the flux direction, is used. |
| D-82 | Ivascenko et al. 2016, journal: one reading gives "ApJ Supplement" / "ApJS (volume not printed)", another ApJ 833, 223; the DOI prefix 10.3847/1538-4357 is that of ApJ. | references (source of 11.1, 11.5, 11.7, 16.9) | ApJ (from the DOI); reference record PARTIAL. |
| D-83 | van den Berg et al. 2024, DOAJ record: issue date "Jan 2024" for vol. 977, which is December 2024. | references (source of Section 11) | Metadata flag. |
| D-84 | Tautz 2013b Eq. 23: D_μμ ≈ (1−μ₀²)D is printed with the condition "t ≪ D"; dimensionally t ≪ 1/D. | 11.2, 16.9, 16.10 | As printed. |
| D-85 | van den Berg et al. 2024, Sects. 4.1–4.3: Eq. 14 "had to be multiplied by 2 to match QLT" in the weak-turbulence cases but not in the strong case. Derived: Eq. 14 is exact only for D_μμ ∝ (1−μ²), where it needs no factor; the exact relation of Section 11.3 gives 1/(2−q) = 3 for the pure QLT shape and 1.70–2.18 for H = 0.05 (F-SIM-02). | 11.3, 11.7, 16.9, 16.10 | The printed factor 2 is not adopted; the exact relation of Section 11.3 is used (U-42). |
| D-86 | van den Berg et al. 2024, M5 for weak slab turbulence: the μ_S = 0.5 run is used for μ < 0 and the −0.5 run for μ > 0; the extractor flags the assignment as possibly reversed. | 11.5 | Recorded as extracted. |
| D-87 | Dundovic et al. 2020: the slab grid is given as N_z = 2²⁵, L_z = 400 kpc, l_slab = 10⁻⁴L_z in one place and as 2³⁰, 384 kpc, 10⁻⁶L_z in another; the 2D box ratio is given as 100 and as 10³. | 11.5 | Both recorded. |
| D-88 | Lange et al. 2013: L_scale ≈ 10⁸ cm in Eq. 29 and 3.4 × 10⁸ cm in the text; the sign of ω differs between Eq. 22 and Eq. B.2. | 11.5 | Both recorded. |
| D-89 | Riordan & Pe'er 2019: the α̃ peak is "about 3" in the text and about 2.5 by the Fig. 6 caption; α̃ ≈ 0.7 (caption) versus 0.5 < α̃ < 1 (conclusions). | 11.5 (source) | Both recorded. |
| D-90 | Tautz et al. 2013a: SOQLT χ² for (δB/B₀ = 10⁻¹·⁵, R = 1) is 14.05 in the A&A HTML text and 27.71 in the PDF Table 1. | 11.6 | Both in the companion data, flagged; not resolved. |
| D-91 | Tautz et al. 2013a, Table 1: no SOQLT entry for (10⁻¹·⁵, R = 10). | 11.6 | Recorded as missing. |
| D-92 | Focusing corrections to λ∥: ξ² coefficients −1/15 (Litvinenko 2012a Eq. 5), −2/5 (Tautz et al. 2012 Eq. 10 = He & Wan 2013 = MFP v1.2 Eq. 18; also Shalchi 2011 as quoted by Litvinenko) and +2/15 (Wang & Qin 2018 Eqs. 62–63). They rest on different definitions (Section 12.7); Litvinenko 2012a states that the −2/5 coefficient is in error. | 12.2, 12.4, 12.5, 12.7 | Separate named options; no default (U-46); companion note C-5. |
| D-93 | Litvinenko 2012b Eq. 24: λ₀ = [3/(2(2−q)(4−q))]v/D₀, although D_μμ of Eq. 4 contains σ. Eq. (23) with Eq. 4 gives the extra factor 1/(1−σ²), consistent with the recovered Eq. 21 at q = 1 (derived). | 12.3 | Eq. 24 read as σ = 0; flagged. |
| D-94 | Effenberger & Litvinenko 2014 Eq. 30 (both renderings): μ_{t+Δt} = [½ξ(1−μ_t²) − μ_t]Δt + √(Δt(1−μ_t²))ε_t − ½μ_tΔt(ε_t² − 1). The leading μ_t is absent; the printed right side is the increment, and with μ_t the equation is the Milstein step (derived). | 12.3, 15.7, 16.3, 16.10, 16.11 | The Milstein step of Section 16.3 is implemented; the printed form is kept as a record only (SDE-MIL-EL14-PRINTED, FOC-EL14-MILSTEIN-PRINTED). |
| D-95 | Effenberger & Litvinenko 2014, equation numbering: ar5iv labels dz = μv dt as Eq. 27; the PDF extraction reports Eq. 27 missing. | 12.3 | Eqs. 28–30 cited. |
| D-96 | Wang & Qin 2018 Eqs. 23/25 and A5: typos flagged by the extractor ("ρ^{M(μ)}", "κ_κ0", "νL"). | 12.4 | e^{M(μ)} used. |
| D-97 | Tautz et al. 2012 Eq. 9: s is not defined in the reading available. | 12.5, 12.7 | U-45; Section 12.7 uses s = 5/3 only as an illustration. |
| D-98 | Tautz, Shalchi & Dosch 2014 Eq. 29a: D(s,q)[π/(v\|μ\|)](δB/\|B\|)²(1−μ²)g(1/(R\|μ\|)); the units look inconsistent as extracted (the variables may be normalized). | 12.6 | Not used. |
| D-99 | Chen et al. 2024, inline definitions: Ω₀ = qB₀/mc, without a Lorentz factor. For a power-law PSD of index s, D_μμ is then larger by γ^{2−s} than with Ω₀/γ (1.274 for 1 GeV protons, s = 5/3; F-MEAS-01). | 13.2, 13.9 | U-47. |
| D-100 | Afanasiev et al. 2015, Sect. 3 (inline): the gap below μ_min is filled with ν̃(μ) = bμ + b₁ on (−μ_min, μ_min), printed in μ, not \|μ\|; b and b₁ are not given. The odd term makes ν̃ asymmetric, and possibly negative for μ < 0, unless b and b₁ are constrained. | 13.7, 14.6 | Recorded as printed; the parameters are U-49. |
| D-101 | Borovikov et al. 2019 Eqs. 4.16, 4.27, 5.6 (PDF only): raw strings with lost fraction bars; the placement was inferred by the reader. | 14.1, 14.3 | Not implemented. |
| D-102 | Vainio 2003 Eq. 1: ⟨Δμ⟩/Δt in the HTML, Δμ/Δt in the PDF rendering. | 14.3 | ⟨Δμ⟩/Δt used (consistent with Battarbee et al. 2011 Eq. 7). |
| D-103 | Battarbee et al. 2010, Sect. II.1: ∂P̃/∂t + v ∂P̃/∂r = … with a lower-case v that is not defined there; Battarbee et al. 2011 Eq. 6 has V = u_sw + v_A. | 14.4 | Battarbee et al. 2011 Eq. 6 used; recorded as printed. |
| D-104 | Ng, Reames & Tylka 2003, Table 1: ambient intensity "2.5 × 10[?]3" MeV cm⁻²; the sign of the exponent is lost. | 14.4 (value not used) | Not used. |
| D-105 | Sandroos et al. 2014 Eqs. 34–35 (as extracted): a₀ = (1 − η)a(μ̃_t), a₁ = (1 − η)a(μ̃_trial). With the same N in both stages the mean drift is 2(1 − η)a, correct only for η = 0.5 (derived, F-SGW-02). | 14.6, 16.3 | Accepted only at the published θ = 0, η = 0.5 (U-51). |
| D-106 | Sandroos et al. 2014, after Eq. 36 (PDF, as reported, one rendering): "μ̃ → 1 − μ̃" at \|μ̃\| = 1. The mirror reflection at μ̃ = 1 is 2 − μ̃; 1 − μ̃ maps 1.01 to −0.01. | 14.6, 16.4, 16.10 | Flagged; not implemented without the original; the Monte Carlo checks of Section 16 use the mirror rule (U-51). |
| D-107 | Afanasiev et al. 2025 Eqs. 29–30 (arXiv v1 only): Δx(E) = C₁ exp(E/ΔE₁)^{α₁} and q(E) = C₂ exp(E/ΔE₂)^{α₂}; the exponent may apply to E/ΔE or to the exponential. | 14.7 (fit not transcribed) | Not resolved; partial input to MFP v1.2 U-6 (companion note C-9). |
| D-108 | van den Berg et al. 2026 Eq. 9 (κ∥): tool-reconstructed, single reading; only fragments are legible, and the bracket grouping is not established. | 15.2 | At P = P_r the value is independent of the grouping (F-GCR-PA-02); for P ≠ P_r the grouping is U-53. |
| D-109 | van den Berg et al. 2026, Sect. 3.1 (Eq. 10): θ_t = "45°/2π" (as reported, one reading); about 7.16° if read literally. | 15.2 | U-54. |
| D-110 | van den Berg et al. 2026 Eq. 23 (text layer): "M ← sign(M) 2 − M if \|M\| > 1"; parentheses and bars are not visible. Every reading that keeps \|M\| ≤ 1 equals 2 sign(M) − M (derived). | 15.2, 16.4 | Treated as the mirror rule; flagged unverified. |
| D-111 | Qin & Wu 2020 and Wu & Qin 2020: equation numbers differ between arXiv versions (D_μμ is Eq. 15 in PDF v2 and Eq. 17 in ar5iv; Eq. 10 versus 12 in paper I). | 15.3 | Both numbers cited. |
| D-112 | SEP_propagator code (van den Berg et al. 2020; used by Strauss et al. 2022; master branch, one reading): the code has `(ABS(MU(j)))**(1.67 - 1.) + 0.05`, so q = 1.67, whereas the papers print q = 5/3 and H = 0.05; the comments mention H = 0.005; the header says "M = 99 and N = 200", the code has M = 49, N = 100. | 15.4 | Recorded only; the paper values define presets (U-57). |
| D-113 | Florinski et al. 2015 (ICRC) Eq. 3.2: the text layer is garbled; a tool reading gives D = D₀(B/B₀)Q/(γA). | 15.6 | Not used. |
| D-114 | Strauss, Fichtner & le Roux 2015 (ICRC): D⊥(μ) = (v r_L/2)(1−μ²)δB_z²/B₀², while δB_x, δB_y ≫ δB_z is assumed. | not cited in a section | Listed for completeness. |
| D-115 | He & Wan 2015 (also He et al. 2017, He & Wan 2019), backward μ SDE: +2D_μμ/(M + μ) in the drift, with M not defined and no loss rate or weight printed. The form is that of the Zhang et al. importance-sampling drift with a → M, but without the loss rate c and the (1 + μ/a) factors the system does not sample f. | 16.2, 16.10, 16.11 | Not used without a defined M and a weight (SDE-IS-HEWAN15-PRINTED blocked). |
| D-116 | Zhang et al. 2023 Eq. 13 (ar5iv): one reading records dμ′(s) = [−dμ′/dt + ∂D_μμ/∂μ]ds + √(2D_μμ)dw(s); another reading of the same rendering lists Eqs. 1–12 only and states that the μ′ SDE was not printed in the text read; a third reading confirms Eq. 13 as above. | 16.2 | Resolved: Eq. 13 is printed (two of three readings); it is the form that the backward rule of Section 16.1 gives. |
| D-117 | Àgueda 2008 thesis Eq. 3.14 (as legible): γ′ = γ_u(γ − u_cμv/c²). The standard form γ′ = γ_uγ(1 − u_cvμ/c²) has γ on the second term; the prefactor of μ′ is not legible. | 16.8, 16.10 | The derived transformation of Section 16.8 is used and stated. |
| D-118 | Merten & Aerdker 2024, constant-focusing test (Sect. 3.5.1): "the drift velocity is determined by the stable pitch-angle: v_drift = μ₁c₀", where the stable pitch angle is described as a fixed point of the SDE drift; μ₁ is not printed numerically. If μ₁ is the zero of the Itô drift (0.819804 for v/(2L) = 5, D₀ = 1), it differs by 2.5 % from the stationary mean coth ξ − 1/ξ = 0.800091 (derived, F-PUSH-08); the Monte Carlo run MCR-5 converges to the latter. | 16.9 | Re-read: the description fits the drift zero; Section 16.9 uses the stationary mean as the test value (U-64). |
| D-119 | Merten & Aerdker 2024, Sect. 3.5: `D_0=1\,\mathrm{N}^2\,\mathrm{s}`; D_μμ has the unit s⁻¹. | 16.9, 16.10 | Treated as D₀ = 1 s⁻¹ in the test; flagged. |
| D-120 | Sandroos et al. 2014 Eq. 36 (PDF text layer, two passes) prints μ̃₀ as the first term on the right-hand side; the next sentence calls μ̃_t and μ̃_{t+1} the initial and new pitch-angle cosines. | 14.6, 16.3 | Read as μ̃_t in Eq. (160); recorded as printed. |
| D-121 | Strauss & Effenberger 2017 Eq. 25: "x < x_min : x → −2x_min − x = −2 − x". With x_min = −1 the middle expression equals 2 − x; the reflection about x_min is 2x_min − x = −2 − x, the right-hand expression (derived). | 16.4 | The right-hand expressions of Eqs. 25–26 are used (Eq. (161)). |
| D-122 | Tautz 2013b Eq. 18 prints the last exponent as e^{−2dt} (lower-case d; two renderings). Eq. 11 gives e^{−2Dt}, consistent with Eqs. 14 and 22. | 11.2 | Read as e^{−2Dt} in Eq. (103); recorded as printed. |
| D-123 | Kozarev et al. 2013 Eq. 1 (arXiv PDF and ar5iv, one checking pass) prints +v ê_b·∇ln B in the focusing bracket of the ∂f/∂μ coefficient on the left-hand side; Young & Poduval 2025 Eq. A1 prints −v b̂·∇ln B, the sign of every other source (Section 2.1). Same code, different printings. | 2.4 | Eq. (14) follows Young & Poduval; the Kozarev sign is recorded as printed and is not offered as a variant. |
| D-124 | van den Berg et al. 2024 Eqs. 13–14, as read (summarized reading): D_μμ ∝ Δt², Δt⁰ and Δt⁻¹ in the ballistic, diffusive and tail regimes, with ⟨(Δμ)²⟩ ∝ Δt^b and b = d ln D_μμ/d ln Δt + 1. Ballistic ⟨(Δμ)²⟩ ∝ Δt² gives D_μμ ∝ Δt (derived), so "Δt²" refers to ⟨(Δμ)²⟩ or is a misprint. | 11.7 | Recorded; not used for a model. |
| D-125 | van den Berg et al. 2026: Eq. 15 (as read: −∂P̃/∂τ = a·∇P̃ + ½C:∇∇P̃ − ŁP̃) and Eq. 18 (w ← w e^{ŁΔτ}) with Ł of Eq. 37. The phase-space divergence of Eq. (16) equals Ł (F-FTE-03), so the advective form has −Łf, for which the Feynman–Kac weight of a backward SDE is exp(−∫Ł dτ) (derived). The sign relation between Eqs. 15, 18 and 37 is not established; Eq. 36 (D-11) shows the term with the opposite sign in one garbled reading. | 2.5, 15.2, 16.2 | SDE-W-VDB26 not executable until the typeset Eqs. 15, 18 and 37 are re-read; the term is of order (v/c)(u/v)²(v/L) and is omitted unless U-3 includes it. |
| D-126 | Vainio & Kocharov 2001 Eq. A.2 (complete in the PDF rendering only) updates U_+ and U_− with the same sign of ν V_A p μ δt, whereas Γ_± = ±(1/U_±)∫ν V_A p μ f (Eq. (130)) and −ΔE = ∓V_A pΔμ imply opposite signs, unless μ_{j−} is measured along the backward-wave direction, which the text read does not state. | 14.3, 14.6 | Recorded as printed; SGW-COUPLE-VK01 needs U-52 and this sign before use. |
| D-127 | Spanier et al. 2022 Eq. 38 (backward whistlers): prefactor π²Ω²k_c/(Ω_i²B₀²). Derived from Eq. 22 with the δ-function Jacobian (v_A = Ω_i/k_c), the prefactor is π²Ω²k_c/(Ω_iB₀²); with Ω_i² the coefficient is dimensionless unless it is expressed in units of Ω_i, which the text read does not state. Eq. 22 also carries Ω without γ in the prefactor (Section 8.4). | 8.4 | WAVE-SPA22-BW is in REQUIRES_SOURCE_OR_CODE_AUDIT until this unit is established; the audit blocks before U-33 (Section 4.1). |
| D-128 | Afanasiev, Battarbee & Vainio 2015, Sect. 3 (two renderings, PDF text and ar5iv): the local spectral index is printed as q_j = ln(I_w,j+1/I_w,j)/ln(k_j+1/k_j) and the interpolant as I_w(k) = I_w,j(k/k_j)^{−q_j}. The two are inconsistent: the printed interpolant reproduces the second node only if q_j carries the opposite sign (a positive decay index for a decreasing spectrum); with the printed q_j the node value is missed by (k_j+1/k_j)^{2q_j} (a factor 16 for k_j+1/k_j = 2, I_w,j+1/I_w,j = 1/4; F-MEAS-03). The following ν(μ) ∝ \|μ\|^{q_j−1}, ν₀,j ∝ \|k_j\|^{q_j} and Eq. 14 are consistent with the decay-index convention only. Which expression carries the sign error is not established. Found by the external review of version 1.0 (Section 25.6). | 13.7 | An implementation stores the signed slope α_j and interpolates with (k/k_j)^{α_j}, or q_j = −α_j with the printed exponent, never both printed expressions together; MEAS-SOLPACS-GRID stays in REQUIRES_SOURCE_OR_CODE_AUDIT (also U-49). |

### 19.2 Decisions required from the user or host code

| ID | Decision |
|---|---|
| U-1 | He & Wan 2015, 2019 and He et al. 2017 coefficient sets: the first printed form of He & Wan 2015 Eq. 3, −v b̂·∇ln B, which all other sources share (the default of Section 2.2), or in addition the v/L form with the printed L as a named variant (D-1). |
| U-2 | Density advanced or reported by each preset, declared as a composed measure (spatial factor d³x or A ds; momentum factor with or without p²; or Ruffolo's F = d³N/(dz dμ dp)) with its Jacobian to f (Section 16.11): phase-space f (per d³x d³p), linear density Af along a flux tube, j = p²f (per unit p), or their product. This concerns van den Berg et al. 2020 Eqs. 9 and 17, Lang et al. 2024 Eq. 1 and Merten & Aerdker 2024 Eq. 1 (D-2, D-4), and the host's macroparticle weights and cell or segment measures (PS06b). |
| U-3 | Whether relativistic or second-order frame terms are carried: Ł of van den Berg et al. 2026 Eq. 37, the EPREM prefactor (1 − V·e_b vμ/c²) (D-7), Ruffolo's v_sw v/c² terms, and L′ of van den Berg et al. 2026 Eq. 11. |
| U-4 | Definition of d/dt in the frame-velocity terms of the Wijsen et al. 2019 and Afanasiev et al. 2025 coefficient sets: the fluid derivative (as in van den Berg et al. 2020, 2026) or the derivative along the particle path (as in Ng & Reames 2008). The choice changes dμ/dt by μ(1−μ²)∂∥W (derived). |
| U-5 | Whether μ is measured relative to **B** or to the outward field direction (Zhang et al. 2023); it matters for D_μμ that is not even in μ (wave-frame forms). |
| U-6 | AMPS (Tenishev et al. 2022): frames of μ and p; whether the V^sw·∇f term of Eq. 1 is represented by the Lagrangian field-line motion; the axes of ∂V_x/∂x, ∂V_y/∂y, ∂V_z/∂z in Eqs. 3–4 (D-9). |
| U-7 | Whether the printings with μ^{s−1} (Qin & Wang 2015, Wang et al. 2014, Qin & Wu 2020, Wu & Qin 2020) are evaluated as \|μ\|^{s−1}. The configuration says so explicitly; the printed form runs only as GCR-PA-QINWU20-PRINTED (D-13). |
| U-8 | Zhang et al. 2023: whether the speed dependence sits in D₀(x) or in λ∥0(x) when one D₀(x) is used for all energies (D-16). |
| U-9 | Units of R (He-type amplitudes) and of p (Wang–Qin-type amplitudes), needed to make D₀ values portable (D-21). |
| U-10 | Wijsen et al. 2019 Eqs. 6–8: whether m is the rest mass or γm (Ω = \|Q\|B/m reproduces Eq. 7 from Eq. 6; the choice affects the energy scaling of relativistic particles), and the sidedness and normalization of P in Eq. 6 (a two-sided total spectrum gives P = 4πg). |
| U-11 | Whether the library offers any ε ↔ (q, H) conversion; the sources do not state the domain or the weighting of their least-squares fits (Section 5.3). |
| U-12 | Which q to pair with ε = 0.048 if a q-form equivalent of PARADISE is requested; q is not printed (D-17). |
| U-13 | Whether the Àgueda hemispheric transformation (Section 5.3; kernel in Section 16.7) is offered as a solver option for ε = 0. |
| U-14 | Meaning of the λ∥ parameter of PA-DROGE-VA: the prefactor value (as printed) or the value of Eq. (23), which requires inverting the ratio of Section 5.6 (D-26). |
| U-15 | Whether PA-AMPS-II is implemented after its factors are read on the typeset page (D-30). |
| U-16 | Whether PA-EPS-PACHECO-PRINTED is offered at all (the analogue of MFP v1.2 U-5 for the Lang/Lavasa variant; D-33). |
| U-17 | Helicity in slab QLT: σ = 0 only (the default of Sections 6–7); σ(k∥) = σ_m sgn(k∥), the Hermitian-consistent form, which makes D_μμ ∝ 1 − σ_m sgn(Ωμ) and breaks TS2002 Eq. 44; or the (n, j)-resolved spectra of Section 8 (D-34). |
| U-18 | Normalization of the TS spectra: approximate (TS2002 Eq. 37, TS2003 Eq. 5, "≈"; reproduces the TS tables) or exact (Strauss et al. 2017 Eq. 12; factor N = s/[s + ((s−p)/(p−1))(k_min/k_d)^{s−1}], e.g. N = 1.00253 at R = 0.5, Q = 1000). Both are named options; which one is the default. |
| U-19 | Ding et al. 2022 Eq. 8: v\|μ\| or the printed vμ (negative D_μμ for μ < 0) (D-36). |
| U-20 | Afanasiev et al. 2015 Eq. 2: B² (TeX rendering; consistent with the text) or R² (PDF extraction; MFP v1.2 D-30) (D-37). |
| U-21 | Laitinen et al. 2016: normalization integral over all k (S∥ = 2πg; D_μμ equals Eq. (43)) or over k > 0 (S∥ = 4πg; D_μμ is twice Eq. (43)) (D-38). |
| U-22 | Engelbrecht et al. 2022 Eq. 11, if implemented: whether S_xx is the two-sided spectrum of one component (S_xx = 2πg) (D-39). |
| U-23 | Chen et al. 2024 spectrum: P is the one-sided PSD of B_N (then D_μμ = ½ × Eq. (43) with P_s = 2P under axisymmetry, and λ∥ is twice as large) or the one-sided total transverse PSD (D_μμ = Eq. (43)) (D-40). |
| U-24 | Zhang et al. 2023 Eq. 36: reading A (π²Ω², W⊥ = 2g), reading B ((π/2)Ω², W⊥ = 4πg), or not implemented (Eq. 35 only) (D-41). |
| U-25 | Strauss et al. 2017 printed kernel. Resolved in version 1.1: the printed kernel is dimensionally inconsistent and is not a runtime coefficient under the admissibility rule of Section 4.1; DYN-STRAUSS17-PRINTED is a record (REFERENCE_DATA_ONLY) that may be evaluated only in a diagnostic demonstrating its unit dependence. DYN-TS-DT equals the γ²-form (D-44). Listed for traceability. |
| U-26 | Four-component split for Strauss & le Roux 2019, whose Eqs. 10–11 define σ_c and σ_m only as ratios: the factorized g_{n,j} = (G/4)(1 + jσ_c)(1 + nσ_m) (derived; satisfies both), the Spanier et al. 2022 form (Eqs. 25–28, which allows H_L ≠ H_R), or a user-supplied table. |
| U-27 | Strauss & le Roux 2019 in the negligible-damping limit: the printed Eq. 28 or the exact-Jacobian δ-limit (D-48). |
| U-28 | Sign convention for Ω of electrons in Strauss & le Roux 2019; it decides which polarization resonates (check S8-1). |
| U-29 | Coexisting forward and backward waves: scatter in both wave frames every step (Vainio & Kocharov 2001) or sequentially with weights ω_± (Virtanen & Vainio 2005; the ω_± formula is not printed). |
| U-30 | Conversions between printed spectrum conventions (Section 14.2, relations 2–3): I = 4πW (equipartition, same sidedness) and U = 2kW (two-sided W). For Vainio & Kocharov 2001, U = k[W(k) + W(−k)] reproduces Vainio 2003 Eq. 2 at \|μ\| = 1, and U = kW(k) gives half of it. The default is no conversion. |
| U-31 | Battarbee et al. 2011 Eq. 2: f_ci as printed (no γ) or f_ci/γ, which matches Vainio, Laitinen & Fichtner 2003 (D-52). |
| U-32 | Sign convention of signed k (polarization) in stored spectra (D-53). |
| U-33 | Spanier et al. 2022: the exact root y_R(k) of Eq. 9 or the printed approximations, Eq. 16 or 17 (at k = 3: 9.8547, 10 and 12). |
| U-34 | k∥ integration domain (k∥ ≥ 0 or all k∥) for the Alfvénic and slow-mode normalizations of Xu & Lazarian 2018, 2020 (D-60). |
| U-35 | Γ[6.5, x] in Yan & Lazarian 2002 Eq. 8 and 2004 Eq. 12 is not defined in the notes; the upper incomplete gamma function is the natural reading but is not confirmed. |
| U-36 | Gaussian-width convention for Yan & Lazarian 2008: the printed closed form (Δμ² = (1 − μ²)M_A) or the value implied by their integral (2Δμ²) (D-65). |
| U-37 | Yan & Lazarian 2008 Eq. 9: the ar5iv or the PDF bracket, and the prefactor (√π/16 printed; √π/8 from the reduction) (D-66). |
| U-38 | Default SOQLT resonance form for NL-SOQLT-SH09: Gaussian (Shalchi et al. 2009 Eq. 11), step (Eq. 15) or the closed form (Eq. 47). At μ = 0 their ratio is sΓ(s/2)/√π (1.0614 for s = 5/3) in the inertial limit and between about 0.956 and 1.128 outside it (Section 10.6). This document sets no default. |
| U-39 | λ∥ for NL-SOQLT-SH09-CF: Shalchi et al. 2009 Eq. 51 as printed, or Eq. (23) applied to Eq. 47; they give 553.358 l and 532.319 l at R = 0.0363, b = 0.05. |
| U-40 | The mapping from the Subashchandar et al. 2025 P_S(k∥) to g^slab (not stated) and the sign convention of q^S (D-75). |
| U-41 | Whether numbers from Riordan & Pe'er 2019 are divided by 2 before comparison with Eq. (22); the single rendering shows no ½ (D-80). |
| U-42 | Whether the printed factor 2 of van den Berg et al. 2024 Eq. 14 is reproduced as a variant or replaced by the exact relation of Section 11.3 (D-85). |
| U-43 | Acceptance estimator and tolerance for pusher validation. Candidates: one-step or few-step conditional increments with the bias correction of Section 11.2; M5 (SIM-M5-ISO); M2b with the exact relation of Section 11.3. |
| U-44 | Domain policy for L < 0 (ξ < 0). The sources do not discuss L < 0. Settled in version 1.1 for even D_μμ: the reflection z → −z, μ → −μ, L → −L maps the process onto itself (Section 12.3), so FOC-EL14 and Eq. (115) are evaluated at \|ξ\| and the sign of ⟨μ⟩ follows the sign of L, with 1/L carried as a signed inverse length. Still open for D_μμ that are not even in μ (directional or helical wave models): a domain diagnostic unless the configuration admits it. |
| U-45 | Value and definition of s in Tautz et al. 2012 Eq. 9 (D-97). |
| U-46 | Whether, and which, focusing relation is applied when converting λ∥ ↔ D_μμ in a focusing field: FOC-EL14, FOC-WQ18-A, FOC-WQ18-SERIES, FOC-HW13 or none (D-92). |
| U-47 | Chen et al. 2024 gyrofrequency: as printed (Ω₀ = qB₀/mc) or relativistic (Ω₀/γ) (D-99). |
| U-48 | Scattering below μ_min with measured spectra: Chen et al. 2024 print none. Options are the gap fillers of Section 5, the dynamical models of Section 7 or the nonlinear models of Section 10. |
| U-49 | SOLPACS gap-filling parameters b and b₁, and whether the linear fill is in μ or in \|μ\| (D-100). |
| U-50 | Whether to offer the Ng & Reames 2008 cap I ≤ B²/(3πk), and in which spectrum normalization; the readings disagree on whether the cap is printed (C-7). |
| U-51 | Sandroos et al. 2014: whether η ≠ 0.5 is allowed (the extracted a₀, a₁ are valid only at η = 0.5), and the reflection rule at μ̃ = ±1 (printed "1 − μ̃") (D-105, D-106). |
| U-52 | Whether to reproduce the order of steps in Vainio & Kocharov 2001, App. A, as read here (the paper gives no explicit order), and the value of the wave-energy floor (not printed). |
| U-53 | van den Berg et al. 2026 Eq. 9: whether the MFP v1.2 Eq. 33 grouping is adopted for P ≠ P_r (tool-reconstructed, single reading) (D-108). |
| U-54 | van den Berg et al. 2026, θ_t = "45°/2π": the literal value of about 7.16° or another reading (D-109). |
| U-55 | van den Berg et al. 2026: a lower bound on Δτ near \|μ\| = 1 and an entrance treatment for an injection at μ = ±1; the printed step rule gives Δτ → 0 there for D^iso, and exactly Δτ = 0 at μ = ±1 for every shape with D′_μμ(±1) ≠ 0 (F-PUSH-19), so an injected beam at μ = ±1 does not advance. Neither source prints a floor or an entrance rule (Section 16.6). |
| U-56 | The μ = ±1 boundary rule for Qin & Wu 2020 and Kubo et al. 2015; neither paper states one. |
| U-57 | Whether the SEP_propagator code choices (zero-flux μ faces; q = 1.67 instead of 5/3) define the Strauss et al. 2022 preset, or only the paper values do (D-112). |
| U-58 | Importance-sampling parameter a, which Zhang et al. 2023 do not list: set by the user (a > 1) for each run. |
| U-59 | Sub-step order for focusing, cooling and scattering (no source quantifies the splitting error), and for the two isotropic processes of the PARADISE split (order not stated). |
| U-60 | Reflection rule at μ = ±1. The sources use the mirror rule; a single application fails for \|m\| > 3 (q-form drift, large Δt). Settled by the implementation requirement of Section 16.11: the mirror rule applied until the value is inside (folding), with a logged count. The alternatives (step rejection, re-sampling) have no source; listed for traceability. |
| U-61 | Treatment of μ = 0 for the q-form, chosen per preset and recorded: a fixed small Δt with a convergence demonstration (Wijsen et al. 2019 found it impractical); the adaptive step of van den Berg et al. 2026 (it goes to zero at μ = 0 for the q-form); the ε-form instead (van den Berg et al. 2026; a different model); the rotation kernel with the hemispheric split (only for the Àgueda ν(μ); U-13); linear gap filling (Afanasiev et al. 2015; U-49). |
| U-62 | Step-length parameters: Δμ_max of Sandroos et al. 2014 Eq. 39 (value not found), and ℓ of the van den Berg et al. 2026 step rule when it is used with D_μμ shapes other than the ε-form. |
| U-63 | Speed in the redraw time λ/v: Marsh et al. 2013 use the initial speed v₀, Kelly et al. the solar-wind-frame speed. |
| U-64 | Constant-focusing validation value: settled in Section 16.9 and Section 26.3 (configuration C) as the stationary mean coth ξ − 1/ξ; Merten & Aerdker compare with the drift zero (D-118). Listed for traceability. |
| U-65 | AMPS pusher: model the code as printed (Tenishev et al. 2022 Eq. 4) or as the standard Itô step; this needs the code (D-10). |
| U-66 | Symmetry extension and endpoint treatment of the printed Alfvénic, slow-mode and fast-mode asymptotes of Section 9.1 (Eqs. (76), (77), (78) and the slow-mode form), which are not real for μ < 0 as printed and, for Eq. (77), diverge at μ = ±1: the interval on which they are evaluated, whether \|μ\| replaces μ, and what is returned outside (a domain diagnostic by default). The sources print neither. |
| U-67 | Observable of each focused-transport validation run (Section 16.9; PS08, PS15, configuration E): forward stationary angular density ρ₊, mean streaming speed, centred spatial variance (Eq. (115)), backward trajectory histogram ρ₋, or the reconstructed solution of a stated boundary-value problem (manufactured f₀ = e^{−z/L}e^{ξμ}); each run names the one it reports and the reference it is compared with. |

### 19.3 Information not found

The following were sought and not found in the accessible sources; nothing in this document substitutes for them.

**FTE and SDE methods.**
- van den Berg et al. 2020 (primer): Appendices A–F, including B (SDE solver), D (G(μ)) and F (the focusing length, "Eq. 22"); every rendering reached ends before the appendices. The SDE_benchmark folder of the primer code was not readable.
- van den Berg et al. 2026: Appendix B (SDE verification); the expression of the drift v_d(μ); the frame of r.
- Wijsen 2020 thesis, Chapter 3 (SDE integration, boundary conditions, time step); van den Berg 2019 thesis, Appendix B (schemes, variable step).
- Full texts not read (abstract only or not accessible): Zhang, Qin & Rassoul 2009; Qin, Zhang & Dwyer 2006; Dröge et al. 2010 (the FTE as printed); Roelof 1969; Skilling 1971; Isenberg 1997; le Roux & Webb 2009, 2012 (the moving-plasma FTE); Kopp et al. 2012; Ellison et al. 1990 (random walk on a sphere; step distribution); Ng & Wong 1979 (φ_eff); Àgueda et al. 2008 (ApJ); Kloeden & Platen.
- Kóta et al. 2005 (AIP Conf. Proc. 781), the version cited by EPREM: not read. The ICRC 2005 paper that was read has a garbled equation and an illegible scattering factor.
- Equations not legible or not captured: Qin et al. 2005 Eq. 7 (the μ of the streaming term is lost); He et al. 2011 Eqs. 2–3; Ruffolo 1995 Eqs. 10–11 beyond the scattering term, Eq. 13, and the DOI; the matrix layout of Ng & Reames 2008 Eq. 1 (published version garbled); the definition of L in Wang et al. 2014.
- Frames not stated or not captured: He & Wan 2015; He et al. 2011; Laitinen et al. 2016; μ and p in Zhang et al. 2023 (the "corotating" statement was read only in the PDF).
- The definition of d/dt in Wijsen et al. 2019 Eqs. 3–4 and Afanasiev et al. 2025 Eqs. 39–40.
- Values not printed: the importance-sampling parameter a of Zhang et al. 2023; Δμ_max of Sandroos et al. 2014; Δt of Effenberger & Litvinenko 2014; the step of Merten & Aerdker 2024 for the pitch-angle tests; the order of the two scattering sub-steps in PARADISE; the back-transformation of the Àgueda η map (derived in this review).
- Dedicated studies of the step-size convergence of pitch-angle SDEs near μ = 0 or μ = ±1 (Itô versus Stratonovich, Euler–Maruyama versus Milstein): none located.

**D_μμ theory.**
- TS2002 Eq. 25 in full (ALT text truncated; bracket from a single PDF reading) and whether Eq. 22 has a complex conjugate; the definition of σ in TS2002 and TS2003.
- Strauss et al. 2017: the typeset ApJ version and the LaTeX source (robots-disallowed); Eq. 20 was read once.
- Shalchi et al. 2006 (NADT): the typeset version; the symbol glyphs of Eqs. 20, 23, 27, 28 and 34 (OCR with Greek glyphs lost).
- Dröge 2003 (dynamical QLT): abstract only; the V_A form is known only through Minoshima et al. 2026. Engelbrecht 2018 and 2019 dissipation-range papers: not read. Jokipii 1966: not read (metadata unconfirmed). Schlickeiser 1989: not read, so the origin of the (1 ∓ H)(1 ∓ σ) weights and of A₂ in Virtanen & Vainio 2005 Eq. 15 is not checked.
- Nonlinear theories: Shalchi 2005 (Phys. Plasmas, the SOQLT foundation; not open access; article number 052905 unconfirmed); Shalchi et al. 2004 (WNLT, ApJ 616, 617; the WNLT D_μμ was not found); Shalchi 2009 book, Chapters 5 (WNLT) and 6 (SOQLT), abstracts only; Shalchi 2011, Shalchi & Danos 2013 and He & Schlickeiser 2014 (no copies; titles and DOIs not found); Tautz, Shalchi & Schlickeiser 2007 (isotropic-turbulence SOQLT resonance; not searched).
- Focusing: Litvinenko 2012b Eq. 19 is garbled, and Eqs. 17–18 and 21 are only "recovered"; Tautz et al. 2012, the definition of s in Eq. 9 and the exact printed form of Eq. 14, including whether the mean displacement is subtracted; Tautz, Shalchi & Dosch 2014, the "modified relation" that reduces λ∥^(0) by about 0.8.
- Measured and prescribed spectra: Chen et al. 2024, units and sidedness of P(f), PSD estimator and cadence (the OSTI copy was not read); Wijsen et al. 2019, normalization of P, value of q, and whether m is the rest mass; Laitinen et al. 2016, the value of C∥ and the integration limits of the normalization; Zhang et al. 2023, the normalization of W⊥ and a legible Eq. 36; Subashchandar et al. 2025, the relation between P_S and g^slab and the sidedness of P_SUM.
- Compressible turbulence: the validity condition of Chandran 2000 Eq. 17; the step functions H* of Xu & Lazarian 2018 Eq. 37. No application of the Section 9 sources to the inner heliosphere or to SEPs was found in the parts read; only Vainio 2000 (abstract) makes a heliospheric statement, about fast-mode mixtures. "Lynn et al. 2014" was not identified; Lynn et al. 2012 is used.

**Simulations.**
- Qin & Shalchi 2009, 2014: full texts not accessible; their simulated D_μμ values, fits and set-ups are not recorded.
- Sun et al. 2016: abstract only; the M5 formula with its reported sign typo was not read.
- Tautz et al. 2013a: Appendix A is garbled (only A.3 legible); the number of particles is not stated in any rendering read.
- van den Berg et al. 2024: the Summary and Conclusions were not reached (extraction truncated in Sect. 4.4); their QLT Eq. C8 was not captured.
- Weidl et al. 2015: the estimator time window and the μ binning were not reached.
- Ivascenko et al. 2016: ν₀ and ε of the Monte Carlo test input; the time step; the number of μ bins.
- Riordan & Pe'er 2019: volume and article number.

**SEP shapes.**
- Beeck & Wibberenz 1986, the origin of the q–H form: not read (title and initials unconfirmed).
- Àgueda et al. 2008 (ApJ 675, 1601): not readable; its Appendix A φ(ε) is cited by Àgueda et al. 2010. The thesis version (Eq. 3.22) is garbled; φ(ε) here is derived and cross-checked against van den Berg et al. 2026 Eq. 13.
- Dröge et al. 2010: abstract only; the q used with H = 0.05 is not known. Dröge & Kartavykh 2009 and Dröge et al. 2014, 2016, 2018: not read.
- The origin of h = 0.2 beyond Qin et al. 2005 (Qin, Zhang & Dwyer 2006, abstract only).
- The least-squares criterion (domain, weighting) of Àgueda & Vainio 2013 and Wijsen et al. 2019; the ε(q) and ε(q, H) curves of Àgueda & Vainio 2013 (Figs. 2 and 4) were not digitized.
- AMPS: the published version of Tenishev et al. 2022 and its typeset Table 1; le Roux & Webb 2009 and Hu et al. 2017, the sources of Types II and III (abstract only or secondary); Qin et al. 2013, the source of Type V (not read).
- Lang et al. 2024, publisher version (to check the (1−μ)² printing); the Lavasa et al. 2026 publisher HTML agrees with arXiv.
- Schwadron et al. 2010 Eq. 2: the sign of the 3/2 exponent.
- Kunow et al. 1991 (q range 1.3–1.9, average 1.63): known only as quoted.
- Ruffolo-group SEP fits with φ(μ) shapes; no observationally tabulated D_μμ(μ) was found.
- Names cited only through other papers, without reference records: Kunow et al. 1991, Lario et al. 2007, Klimas & Sandri 1971, Jokipii 1966, Zhao & Li 2014.

**Waves and shocks.**
- Vainio 2000 (ApJS 131, 519): abstract only; the wave-frame D_μμ, D_μp and D_pp for forward and backward slab waves with cross helicity were not obtained.
- Ng & Reames 2008: the explicit R^σ_μμ is not printed in the Letter. Ng, Reames & Tylka 2003, Appendices A (broadened resonance function) and B (growth rate): not readable. Ng et al. 1999 (GRL): scanned.
- Virtanen & Vainio 2005: the explicit formula for the weights ω_±.
- Sandroos et al. 2014: the rule for choosing the parallel or antiparallel wave frame for each scattering; the explicit frame transformations; the journal version.
- Strauss & le Roux 2019: the four-component split actually used; the species of Ω in Eqs. 25 and 29; the meaning of "maximal Larmor radius".
- Schreiner & Spanier 2015: the definition of ω_p in Eq. 2.11.
- Afanasiev et al. 2015: how the gap-filler coefficients b and b₁ are fixed. Afanasiev & Vainio 2013 (ApJS 207, 29), the SOLPACS growth numerics: access blocked.
- Vainio & Laitinen 2007, 2008: abstract only. Afanasiev, Vainio & Kocharov 2014; Battarbee 2013 PhD thesis: not read.
- Lee 1983, Lee 2005, Bell 1978 (steady-state theory): cited only.
- PATH self-excited-wave formulas: Zank et al. 2000; Rice et al. 2003; Li et al. 2003, 2005; Gordon et al. 1999.
- Berezhko-school self-generated-wave models: only Petukhova et al. 2017 was found, with prescribed waves; it is not cited.

**GCR.**
- van den Berg et al. 2026: Appendix C (anisotropic distribution and A₁); the LIS formula (Eq. 5); the energy range of the bins; the arXiv identifier (a listing gives 2606.09298, which is absent from the OpenAlex record).
- Zhang 1999: full text.
- Kóta & Jokipii 2017; Florinski et al. 2015 (ApJ) and 2024: abstracts only.
- Engelbrecht et al. 2022, Sections 7–8.
- Qin & Wu 2020: the published version; the pseudo-particle number, time step and μ boundary.
- Strauss & Fichtner 2014: the D_μμ amplitude (figure only) and the Dröge et al. 2010 D_μμ it uses.
- No pitch-angle-resolved model of the GCR anisotropy at 1 au other than van den Berg et al. 2026, and no Forbush-decrease model with pitch-angle scattering other than Qin & Wu 2020, was found.

### 19.4 Notes on the companion specifications

These items concern MFP v1.2 or PARALLEL rev. 1.4. They are recorded, not resolved; a change to a companion is for its editor.

| ID | Companion location | Finding | Evidence |
|---|---|---|---|
| C-1 | MFP v1.2, Section 6.2 (table of H values) and Section 6.3 | Additions, not disagreements: H = 0.5 (Kartavykh et al. 2015, ICRC), H = 0.05 with q = 1.67 (Dröge & Kartavykh 2008, ICRC), h = 0.01 (Qin & Wu 2020) and ε = 0.045 (van den Berg et al. 2026) are not in the MFP v1.2 tables. The ε–(q, H) values of MFP v1.2 Section 6.3 agree with the sources as printed. An exploratory evaluation in Section 5.3, under one assumed criterion, suggests that ε = 0.048 (Wijsen et al. 2019) and ε = 0.045 (Àgueda & Vainio 2013) need not come from the same criterion. | Values as printed (Sections 5.2–5.3); exploratory evaluation not a preset. |
| C-2 | MFP v1.2, Section 6.4 (also the SEP-ZHANG23 row of Section 5.2) | Zhang et al. 2023 is listed with "λ∥ = λ∥0(x)p^{2−q} (printed)" without the note that Eq. 35 has no speed factor, so that λ∥0 must contain v (D-16, U-8). | Eqs. 35 and 37: arXiv PDF, two fetches; the missing factor is derived from Eq. (23). |
| C-3 | MFP v1.2, Section 6.7 (AMPS table) | MFP v1.2 transcribes μ^{2/3}, μ^{q−1} and μ^{s−1}; the readings read \|μ\| with bars stripped in places (D-29). | Both readings come from the same partly garbled arXiv extraction; the typeset table was not read. |
| C-4 | MFP v1.2, Section 6.8 (Zhang et al. 2023 Eq. 36 row) | MFP v1.2 transcribes π²Ω²(1−μ²)W⊥(k_res)/(B²v\|μ\|) as printed. This review finds Eq. 36 not legible: the readings π²Ω² and (π/2)Ω² differ by 2π (D-41, U-24). This review suggests marking the MFP v1.2 row "not legible". | Two fetches of the arXiv PDF that disagree; one raw string with undeterminable fraction grouping. |
| C-5 | MFP v1.2, Section 6.9 (MFP-FOCUS-HW13, Eq. 18) | MFP v1.2 calls He & Wan 2013 "a published isotropic-scattering result". This review finds the He & Wan function identical to Tautz et al. 2012 Eq. 10; its ξ² coefficient (−2/5) is the one that Litvinenko 2012a states is in error for Shalchi 2011 (D-92). Suggested: MFP v1.2 Section 6.9 cross-references Section 12.5; no formula change. | Identity derived and checked numerically; Litvinenko 2012a Full text ×1; Shalchi 2011 not accessible. |
| C-6 | MFP v1.2, Section 7.6 | MFP v1.2 states that "No new closed-form λ∥ was found in the 2006–2026 nonlinear literature that was accessible". Shalchi et al. 2009 print closed-form λ∥ for the strong-turbulence limit (Eq. 49) and the weak-turbulence limit (Eq. 51); Section 10.3 records both. Suggested: MFP v1.2 Section 7.6 cross-references Section 10.3. | Eq. 49: aanda.org alt-text, Full text ×1; Eq. 51: one HTML fetch. |
| C-7 | MFP v1.2, Section 8.2, and its NgReames2008 reference record | MFP v1.2 states that Ng & Reames 2008 "cap the wave intensity at I_σ(k) ≤ I_sat(k) = B²/(3πk)"; the record's role is "Wave saturation cap and ambient λ". One reading (one fetch of the NASA preprint) lists the same cap. Another reading (three fetches of the same preprint; the published PDF is garbled) found no explicit cap of this or a similar form. The readings disagree. Section 14.4 reports both and gives the cap no model ID (U-50). A resolution needs one targeted re-read of the preprint (NTRS 20080038680; search "I_sat", "saturat", "3π"). No normalization comparison was attempted: Ng & Reames normalize ⟨δB·δB⟩ = Σ_σ∫₀^∞I_σ dk (one-sided, per mode), Afanasiev et al. 2015 use a two-sided I_w. | One fetch versus three fetches of the same preprint. |
| C-8 | MFP v1.2, Section 8.2 (agreements) | Agree with MFP v1.2: λ = 0.23 AU at 1 MeV and 0.68 AU at 100 MeV for an ambient I_σ ∝ k^{−1.5} (Ng & Reames 2008); I_R− = I_L− = 0.1 I_R+ (Ng, Reames & Tylka 2012); the Afanasiev et al. 2015 initial spectrum ∝ k^{−3/2} normalized to λ = 1 R☉ at 100 keV; n₀, B₀, V_s, E_inj and r_c = 4; Afanasiev et al. 2015 Eq. 7 inserted in Eq. 4 gives Eq. 12 (re-checked, F-SGW-03). Ng & Reames 2008 themselves use I_R− = I_L− = 0.05 I_R+; this is a different paper, not a conflict. | Two sets of notes agree on 0.05; Eq. 12 re-derived. |
| C-9 | MFP v1.2, U-6 and Section 8.3 (PARASOL Λ(E), Δx(E)) | The notes contain: Eq. 18, Λ(E) = K(E/E_inj)^{β_λ}exp[(E/E_b)^{δ_λ}]; Eq. 19, X_λ = [1 + (x/Δx)^q]^{1/q}; from arXiv v1 only, K = 0.22V_A/(ε_injΩ₀) (Eq. 28), Δx(E) = C₁exp(E/ΔE₁)^{α₁} with C₁ = 0.35((M_A/(ε_inj x_ref))(V_A/Ω₀))^{0.71} R☉, α₁ = 0.83 − 0.02E_b/E_ref, ΔE₁ = 1.64(E_b/E_ref)^{1.07} MeV, and q(E) = C₂exp(E/ΔE₂)^{α₂} with C₂ = 0.66, α₂ = 0.47, ΔE₂ = 9.45(E_b/E_ref)^{0.82} MeV; β_λ = 0.2, δ_λ = 1.4, x_ref = 1 R☉, E_ref = 1 MeV. The exponent placement in Δx and q is ambiguous (D-107), and E_b is not defined in the notes. This is partial input for MFP v1.2 U-6; the companion specification decides whether it suffices. | Eqs. 18–19: two renderings; Eqs. 28–30 and the parameter values: arXiv v1 only. |
| C-10 | MFP v1.2, Section 9.2 (NWU form, Eq. 33) | van den Berg et al. 2026 Eq. 9 is an instance of the NWU form with a numerical κ₀ (1.3 × 10²³ cm² s⁻¹, without 1/3, β included, B₀ = 1 nT, P_r = 1 GV). It is not in MFP v1.2 because it was published later; the companion data of this document hold it (`selfgenerated_wave_and_gcr_sets.json`). | Eq. 9 tool-reconstructed, single reading (D-108); at P = P_r the value does not depend on the grouping. |
| C-11 | MFP v1.2, register D-30 (and Section 6.8) | D-30 records R² in the denominator of Afanasiev et al. 2015 Eq. 2 "as printed in the arXiv version" (journal version not checked). This review reads B² in the ar5iv rendering of the TeX source; one reading also records Eqs. 1–5 with B² in the A&A journal HTML; R² appears only in the PDF text extraction, which another reading calls an extraction artifact. One reading suggests amending D-30 to "PDF text extraction shows R²; TeX and journal HTML show B²"; another suggests closing D-30 as an extraction artifact; this document keeps the choice open (U-20). Not resolved here (D-37). | ar5iv (TeX): two readings; journal HTML: one reading; R²: PDF text extraction only. |
| C-12 | MFP v1.2, register items re-confirmed | MFP v1.2 D-2, D-3, D-4, D-6, D-7 and D-29 were re-confirmed by the readings used here; they appear as D-8, D-32, D-33, D-14, D-28 and D-27. MFP v1.2 U-5 and U-8 are referred to and not repeated. | As stated in those entries. |
| C-13 | MFP v1.2, Section 21.5 | MFP v1.2 states that under pure scattering "no boundary condition is needed because D_μμ(±1) = 0". This holds for the exact process (Feller's test, Section 16.4) but not for a discrete step: an Euler–Maruyama step from the boundary layer leaves [−1, 1] with probability up to 1 − Φ(√(2 − νΔt)) per step (0.084 at νΔt = 0.1, tending to 0.0786 as Δt → 0; F-PUSH-10). Suggested addition to MFP v1.2 Section 21.5: "for the exact process; discrete SDE steps need the reflection of Section 16.4". | Derived and checked numerically (F-PUSH-10). |
| C-14 | PARALLEL rev. 1.4, Section 7.3, Eqs. 23–24 | PARALLEL rev. 1.4 defines I(q_μ, h_μ) over [−1, 1], with D₀ = (3v/8λ)I and I(q, 0) = 4/[(2−q)(4−q)]. MFP v1.2 Eq. 13 and Section 5 define I over [0, 1], with D₀ = (3v/4λ)I and I(q, 0) = 2/[(2−q)(4−q)]. The D₀ is the same; the notation differs by a factor 2. The conventions paragraph of Section 5 states this once. | Text of both documents; PARALLEL Eqs. 23–24 re-read for this register. |
| C-15 | PARALLEL rev. 1.4, Section 9.4, Eq. 42 | For C_k(t) = exp(−Δ²t²/4), Eq. 42 gives (√π/Δ)exp(−ω²/Δ²). With Δ = \|k∥\|Δv∥ it reproduces the printed closed form of Yan & Lazarian 2008 Eq. 6, not their integral, whose exponent −½k∥²Δv∥²t² corresponds to Δ = √2\|k∥\|Δv∥ (D-65). No change to PARALLEL rev. 1.4 is needed; a mapping statement is suggested. | Derived (check S9-8). |

## 20. Model selection

These are factual consequences of the sources and of the derived checks of this review; they are not recommendations of one model over another.

- **Reproducing a published SEP or GCR simulation.** The paper's shape, its gap parameter under the paper's own symbol (H, h, h₀ or ε), its amplitude mode (from λ∥, from λ_r with ψ, from a rigidity or momentum scaling, or from turbulence) and its operator convention are all needed (Section 5; `dmumu_shape_presets.json`). The gap parameters are not interchangeable: the published correspondences between ε and (q, H) rest on unstated fitting criteria (Section 5.3, U-11). The same λ∥ with different shapes gives different D_μμ(μ) and different pitch-angle distributions.
- **Turbulence-based D_μμ.** Magnetostatic slab QLT (Section 6) gives D_μμ ∝ |μ|^{s−1} → 0 at μ = 0 for an inertial range with 1 < s < 2, with a finite λ∥; a dissipation range with p ≥ 2 or a hard upper cut-off k_max makes λ∥ infinite (Section 6.6). The models that give a nonzero D_μμ(0) are separate models with their own assumptions and inputs: dynamical turbulence (Section 7.3, finite for α > 0), second-order and nonlinear theories (Sections 10.3 and 10.6), resonance broadening and transit-time damping in compressible turbulence (Section 9), and the empirical gap parameters (Section 5). Test-particle simulations agree with SOQLT near μ = 0 in some regimes and with neither QLT nor SOQLT in others (Section 10.6); no source read establishes one of these models as generally valid.
- **Wave-frame scattering.** In wave-frame models μ is measured in the frame of the forward or backward waves, and each scattering changes the plasma-frame energy (Sections 8.2, 14.1 and 16.8); models that scatter in the solar-wind frame have no such term. The relative weight of forward and backward waves, their helicity and the sidedness of the stored spectrum are inputs, several of which the sources do not print (U-26 to U-30, U-32).
- **Measured spectra.** Converting a measured power spectral density into D_μμ requires the sidedness and component content of the spectrum, the frame and the Taylor-hypothesis speed, and, for Chen et al. (2024), whether the gyrofrequency includes γ (Section 13; U-10, U-21, U-23, U-24, U-40, U-47). The μ_min cut-off of Chen et al. applies to the κ∥ integral only; below μ_min the source prints no D_μμ (U-48).
- **GCR.** Few pitch-angle-resolved GCR models were found (Section 15): van den Berg et al. (2026) use the ε-form with ε = 0.045 and an amplitude from κ∥, Strauss et al. (2022) the q-form with H = 0.05, Qin & Wu (2020) a slab-QLT q-form printed with μ^{s−1} (D-13, U-7). Most GCR modulation codes solve the Parker equation and need only κ∥ (MFP v1.2, Section 9).
- **Pusher numerics.** For isotropic scattering the Euler–Maruyama errors of ⟨μ⟩ and ⟨μ²⟩ fall roughly in proportion to Δt and reach 5 × 10⁻⁴ near νΔt = 0.002; with reflection at μ = ±1, the stationary ⟨μ²⟩ bias has the opposite sign to the no-boundary value and about three times its size (2.7 to 3.4 in time-averaged runs) for 0.005 ≤ νΔt ≤ 0.02, and the running κ(t) at νt = 50 is 3.7 % below the exact κ(50) at νΔt = 0.01. The random-rotation kernel has no stationary bias but realizes only isotropic scattering, or the ε-form through the hemispheric mapping and an isotropic second process (Sections 16.6 and 16.7). Near μ = 0, Euler–Maruyama with the q-form depletes the distribution: at D₀Δt = 0.01 by 21% (H = 0.05) and 31% (H = 0) in |μ| < 0.01, against 8% for the ε-form; the ε-form is within two standard errors from D₀Δt = 0.003, while the q-form with H = 0 is still 6% low at 10⁻⁴ (Section 16.5). These numbers are derived for the stated shapes and are not statements about any published code.
- **Focusing.** With focusing, the spatial spreading produced by a pusher differs from Eq. (23), and the published relations of Section 12 are closure coefficients, not the quantity a pusher measures: for isotropic scattering the diffusion-approximation closed form κ∥ = λ₀v(coth ξ/ξ − 1/ξ²) gives 0.939 of the unfocused value at ξ = 1, whereas the long-time centred spatial variance of the same process gives 0.744 (Eq. (115), Section 12.3; 0.535 against 0.876 at ξ = 1.5). The mean streaming speed v(coth ξ − 1/ξ) is the same in both. A pusher result is compared with the reference of the same observable; the other published relations differ from both, and none is applied as a default (U-46).
- **D_μμ from simulations.** The estimators of Section 11 differ in normalization and in their bias at finite lag; a simulated D_μμ is compared with Eq. (22) only after the estimator's normalization and lag are stated (Sections 11.2, 11.7; U-41, U-42, U-43).
- **Code descriptions.** Some printed pusher equations cannot be implemented as printed (Section 16.10). Reproducing such a code requires its source code, not the printed equation (U-65).

## 21. Publication methods text

The following paragraphs are templates. Bracketed items are filled from the run record (Section 17), and every equation cited is given with its number in this document and its original source. The templates cite this document's reference identifiers (Section 22).

**Prescribed q-form.** "Pitch-angle scattering was described by the term ∂/∂μ(D_μμ ∂f/∂μ) of the focused transport equation, with D_μμ = D₀(1 − μ²)(|μ|^{q−1} + H), q = [ ] and H = [ ], following [source of the shape and parameters]. The amplitude D₀ was fixed by the prescribed parallel mean free path through λ∥ = (3v/8)∫₋₁¹(1 − μ²)²/D_μμ dμ, evaluated [with the closed form of Eq. (38) for q = 5/3 / with Eq. (37) / by quadrature], with λ∥ = [ ] at [energy and position]." If the source quantity is λ_r, add "λ_r = λ∥cos²ψ with ψ the local field angle (MFP v1.2, Section 1.3)".

**Prescribed ε-form.** "Pitch-angle scattering in the solar-wind frame used D_μμ = (ν₀/2)(|μ|/(1 + |μ|) + ε)(1 − μ²) with ε = [ ] [P55, P25]; ν₀ was fixed by λ∥ = (v/ν₀)φ(ε) with φ given by Eq. (39)."

**Slab quasilinear coefficient.** "D_μμ was computed from magnetostatic slab quasilinear theory, Eq. (43), for a slab spectrum g(k∥) normalized as δB²_slab = 8π∫₀^∞g dk∥ [P72, P59], with B₀ = [ ], δB²_slab = [ ], s = [ ], k_min = [ ] [and dissipation range k_d = [ ], p = [ ]]. The spectrum [name] was converted to g by [rule of Section 6.7 or 13.9]." Name the treatment of μ ≈ 0 (dynamical broadening, second-order theory or gap parameter, with its source).

**Itô pusher.** "The pitch-angle evolution was integrated as the Itô stochastic differential equation dμ = (dμ/dt + ∂D_μμ/∂μ)dt + (2D_μμ)^{1/2}dW [P44, P2], with the [Euler–Maruyama / Milstein] scheme and a time step [Δt = [ ] / chosen by the rule of [P14] with ℓ = [ ]]. Values of μ outside [−1, 1] were reflected, μ → 2 − μ for μ > 1 and μ → −2 − μ for μ < −1 [P44, P9]; the reflection was applied repeatedly until the value was inside the interval. The forward SDE samples the [flux-tube density / phase-space density], and [weights]. Convergence was verified against ⟨μ⟩(t) = μ₀e^{−νt}, ⟨μ²⟩(t) = 1/3 + (μ₀² − 1/3)e^{−3νt} and κ∥ = v²/(3ν) for isotropic scattering [and, for constant focusing, against the stationary ⟨μ⟩ = coth ξ − 1/ξ, the mean streaming speed v(coth ξ − 1/ξ) and the centred spatial variance Var[z]/(2t) → κ_var of Eq. (115); the diffusion-approximation coefficient of Effenberger & Litvinenko (2014) was not used as the reference for the variance], and the pitch-angle distribution near μ = 0 was checked for convergence with Δt. [The scattering rate was applied per [host / plasma-frame / wave-frame] time, and the forward pseudo-particles sampled [the flux-tube density A(s)p²f / f] (Section 16.11).]"

**Time-backward integration.** "The transport equation was solved with time-backward stochastic differential equations [P15, P14], in which the deterministic rates change sign while ∂D_μμ/∂μ and the noise term keep their sign; [importance sampling with f = (1 + μ/a)u, a = [ ], and the loss rate of [P15] / no importance sampling]."

**Random-rotation kernel.** "Isotropic scattering with frequency ν was applied as small-angle random rotations [P55, P25]: in each step δt the squared deflection ϑ² was drawn from an exponential distribution with mean 2νδt and the azimuth φ uniformly, and μ′ = μcosϑ + (1 − μ²)^{1/2}sinϑ cosφ, with νδt = [ ]." For the ε-form add "using the hemispheric mapping of [P55] for the |μ|/(1 + |μ|) part and an isotropic process of rate εν₀, in the order [ ]".

**Wave-frame scattering.** "Particles were scattered in the rest frames of the forward and backward Alfvén waves [P88] with D_μμ^± = [ ]; the particle was transformed to the wave frame, scattered, and transformed back [Lorentz / Galilean transformation], so that each scattering changes the plasma-frame energy. The wave spectra [I_±(k), sidedness, helicity] were [prescribed / evolved by [model of Section 14]]."

**Reporting checklist.** A methods section that uses this library states: the model IDs; the operator convention and the conversion applied (Section 3.3); the frame and reference direction of μ; the gap parameter under its source symbol; the amplitude mode and the λ quantity (λ∥, λ_r or λ_rr); for spectra, the printed quantity, sidedness, component content and conversion rule; the kernel, time direction, boundary rule, step rule and step parameter; the treatment of μ = 0; the sampled measure and weights (Section 16.11); the clock the scattering rate refers to (Eq. (173)); the observable reported and the reference it was compared with (Section 16.9); the number of pseudo-particles and the random-number generator and seed; and the open D-/U-items that apply to the run.

## 22. References

References are numbered [Pn] in order of first citation. Each entry gives the role of the paper here, the reading level of the content used (Section 1.3) and the citation key: the existing key when the paper is in the supplied bibliography, otherwise the proposed key of the additions file (Section 23). Bibliographic fields were checked against publisher, Crossref, OpenAlex or arXiv records (Section 25.3); an entry marked "Metadata not fully confirmed" has fields that could not be confirmed and are checked before a publication cites it.

**[P1]** van den Berg, J. P.; Strauss, R. D.; Effenberger, F. (2020). **A Primer on Focused Solar Energetic Particle Transport: Basic Physics and Recent Modelling Results.** *Space Science Reviews*, **216**(8), 146. [DOI](https://doi.org/10.1007/s11214-020-00771-x); [arXiv:2012.07570](https://arxiv.org/abs/2012.07570). Role: 1-D conservative FTE (Eq. 9), mixed-frame FTE with mixed perpendicular terms (Eq. 17) and coefficients (Eq. 18); contraction convention. Reading: full text. Key in the supplied bibliography: **vanDenBerg-2020-arXiv**. Note: authors: the published byline reads 'Jabus van den Berg, Du Toit Strauss, Frederic Effenberger' (hence 'Strauss, D. T.' in some indexes); the record's initials J. P. and R. D. are the authors' standard initials (OpenAlex display name J.P. van den Berg).

**[P2]** Effenberger, F.; Litvinenko, Y. E. (2014). **The Diffusion Approximation versus the Telegraph Equation for Modeling Solar Energetic Particle Transport with Adiabatic Focusing. I. Isotropic Pitch-angle Scattering.** *The Astrophysical Journal*, **783**(1), 15. [DOI](https://doi.org/10.1088/0004-637X/783/1/15); [arXiv:1410.1225](https://arxiv.org/abs/1410.1225). Role: 1-D FTE with constant focusing length (its Eq. 1 is Eq. (1) here); linear density f = exp(z/L) f0; time-forward Ito SDE. Reading: full text. Key in the supplied bibliography: **Effenberger-2014-AJ**.

**[P3]** Agueda, N.; Klein, K.-L.; Vilmer, N.; Rodríguez-Gasén, R.; Malandraki, O. E.; Papaioannou, A.; Subirà, M.; Sanahuja, B.; et al. (2014). **Release timescales of solar energetic particles in the low corona.** *Astronomy & Astrophysics*, **570**, A5. [DOI](https://doi.org/10.1051/0004-6361/201423549). Role: 1-D FTE with source (Eq. 2), L = -B/(dB/dz). Reading: full text. Proposed key (additions file): **Agueda-2014-AA**.

**[P4]** Strauss, R. D.; van den Berg, J. P.; Rankin, J. S. (2022). **Cosmic-Ray Transport near the Sun.** *The Astrophysical Journal*, **928**(1), 22. [DOI](https://doi.org/10.3847/1538-4357/ac582a); [arXiv:2202.11320](https://arxiv.org/abs/2202.11320). Role: 1-D FTE for GCR protons (Eqs. 2-3), 1/L = -(1/B) dB/dz. Reading: full text. Key in the supplied bibliography: **Strauss-2022-AJ**.

**[P5]** Qin, G.; Zhang, M.; Dwyer, J. R.; Rassoul, H. K.; Mason, G. M. (2005). **The Model Dependence of Solar Energetic Particle Mean Free Paths under Weak Scattering.** *The Astrophysical Journal*, **627**(1), 562–566. [DOI](https://doi.org/10.1086/430136). Role: 1-D FTE "in solar wind frame" and flux-tube density g = A f for time-forward SDEs (Eqs. 1, 7, 8; equations garbled). Reading: full text. Pre-2006 foundation. Key in the supplied bibliography: **Qin-2005-AJ**.

**[P6]** Ruffolo, D. (1995). **Effect of Adiabatic Deceleration on the Focused Transport of Solar Cosmic Rays.** *The Astrophysical Journal*, **442**, 861. [DOI](https://doi.org/10.1086/175489); [arXiv:astro-ph/9408056](https://arxiv.org/abs/astro-ph/9408056). Role: mixed-frame FTE coefficients (Eqs. 7-9); phi = 2 D_mumu operator convention. Reading: full text. Pre-2006 foundation. Key in the supplied bibliography: **Ruffolo-1995-AJ**. Note: pages: first page 861 confirmed; last page not given by OpenAlex (record keeps first page only).

**[P7]** Kartavykh, Y.; Dröge, W.; Gedalin, M. (2015). **Simulation of energetic particle interactions with a shock wave in a focused transport model.** *Proceedings of the 34th International Cosmic Ray Conference (ICRC2015), PoS(ICRC2015)205*, **236**, 205. [DOI](https://doi.org/10.22323/1.236.0205). Role: κ₀(r,E)(|μ|^{q-1}+H)(1−μ²) with q = 1.67, H = 0.5 (Eq. 4.2). Reading: full text. Proposed key (additions file): **Kartavykh-2015-ICRC**. Note: year given as the conference year (ICRC2015); the PoS volume 236 was published 2016-08-18.

**[P8]** Laitinen, T.; Kopp, A.; Effenberger, F.; Dalla, S.; Marsh, M. S. (2016). **Solar energetic particle access to distant longitudes through turbulent field-line meandering.** *Astronomy & Astrophysics*, **591**, A18. [DOI](https://doi.org/10.1051/0004-6361/201527801); [arXiv:1508.03164](https://arxiv.org/abs/1508.03164). Role: 3-D FTE with focusing, cooling, convection, perpendicular diffusion (Eq. 1). Reading: full text. Key in the supplied bibliography: **Laitinen-2016-AA**.

**[P9]** Merten, L.; Aerdker, S. (2024). **Modeling Cosmic-Ray Transport: A CRPropa based stochastic differential equation solver.** *arXiv preprint*. [arXiv:2410.01472](https://arxiv.org/abs/2410.01472). Role: 1-D conservative FTE (Eq. 1), L(s) = -B/(dB/ds); time-forward Ito SDE. Reading: full text. Key in the supplied bibliography: **Merten-2025-CPC**.

**[P10]** Tenishev, V.; Zhao, L.; Sokolov, I. (2022). **Application of the Monte Carlo Method in Modeling Transport and Acceleration of Solar Energetic Particles.** *arXiv preprint*. [arXiv:2209.09346](https://arxiv.org/abs/2209.09346). Role: AMPS Lagrangian FTE (Eq. 1) and printed increments (Eqs. 2-4). Reading: full text. Key in the supplied bibliography: **Tenishev-2022-arXiv**.

**[P11]** Strauss, R. D.; Fichtner, H. (2015). **On Aspects Pertaining to the Perpendicular Diffusion of Solar Energetic Particles.** *The Astrophysical Journal*, **801**(1), 29. [DOI](https://doi.org/10.1088/0004-637X/801/1/29); [arXiv:1804.03689](https://arxiv.org/abs/1804.03689). Role: conservative FTE with L^-1 = div b (Eq. 1); grid solver. Reading: full text. Key in the supplied bibliography: **Strauss-2015-AJ**.

**[P12]** Strauss, R. D.; Dresing, N.; Engelbrecht, N. E. (2017). **Perpendicular Diffusion of Solar Energetic Particles: Model Results and Implications for Electrons.** *The Astrophysical Journal*, **837**(1), 43. [DOI](https://doi.org/10.3847/1538-4357/aa5df5); [arXiv:1804.03693](https://arxiv.org/abs/1804.03693). Role: conservative FTE with L^-1 = div b (Eqs. 1-2); grid solver. Reading: full text. Key in the supplied bibliography: **Strauss-2017-AJ-837**.

**[P13]** Strauss, R. D.; le Roux, J. A. (2019). **Solar Energetic Particle Propagation in Wave Turbulence and the Possibility of Wave Generation.** *The Astrophysical Journal*, **872**(2), 125. [DOI](https://doi.org/10.3847/1538-4357/aafe02); [arXiv:1902.02038](https://arxiv.org/abs/1902.02038). Role: FTE with L^-1 = div b (Eq. 19). Reading: full text. Key in the supplied bibliography: **Strauss-2019-AJ**.

**[P14]** van den Berg, J. P.; Engelbrecht, N. E.; Strauss, R. D. (2026). **Are the Parker and focused transport equations equivalent for galactic cosmic ray modulation?.** *Astronomy & Astrophysics*, **712**, A204. [DOI](https://doi.org/10.1051/0004-6361/202660355); [arXiv:2606.09298](https://arxiv.org/abs/2606.09298). Role: GCR FTE (Eqs. 1-4), L prime (Eq. 11), backward form and weight term (Eqs. 36-37). Reading: full text. Key in the supplied bibliography: **vandenBerg-2026-AA**. Note: arXiv id confirmed through OpenAlex's record of the arXiv DOI 10.48550/arXiv.2606.09298 (title and authors match); the arXiv abstract page itself was refused by the fetch proxy (HTTP 429) and export.arxiv.org was not reachable (robots).

**[P15]** Zhang, M.; Cheng, L.; Zhang, J.; Riley, P.; Kwon, R.-Y.; Lario, D.; Balmaceda, L.; Pogorelov, N. V. (2023). **A Data-driven, Physics-based Transport Model of Solar Energetic Particles Accelerated by Coronal Mass Ejection Shocks Propagating through the Solar Coronal and Heliospheric Magnetic Fields.** *The Astrophysical Journal Supplement Series*, **266**(2), 35. [DOI](https://doi.org/10.3847/1538-4365/accb8e); [arXiv:2212.07259](https://arxiv.org/abs/2212.07259). Role: 3-D FTE (Eqs. 1-4), mu relative to outward field; time-backward SDEs. Reading: full text. Key in the supplied bibliography: **Zhang-2023-AJSS**.

**[P16]** Kubo, Y.; Kataoka, R.; Sato, T. (2015). **Interplanetary particle transport simulation for warning system for aviation exposure to solar energetic particles.** *Earth, Planets and Space*, **67**, 117. [DOI](https://doi.org/10.1186/s40623-015-0260-9); [arXiv:1506.00825](https://arxiv.org/abs/1506.00825). Role: FTE in index notation with frame statement; time-backward SDEs. Reading: full text. Key in the supplied bibliography: **Kubo-2015-EPS**.

**[P17]** Young, M. A.; Poduval, B. (2025). **Quantifying the Effects of Parameters in Widespread SEP Events with EPREM.** *arXiv preprint*. [arXiv:2509.16401](https://arxiv.org/abs/2509.16401). Role: EPREM Lagrangian FTE (Eq. A1). Reading: full text. Proposed key (additions file): **Young-2025-arXiv**. Note: no journal publication found (no journal-ref or DOI on arXiv up to v4 of 2026-03-17; OpenAlex lists only arXiv and Zenodo locations).

**[P18]** Qin, G.; Wu, S.-S. (2020). **Magnetic Cloud and Sheath in the Ground-Level Enhancement Event of 2000 July 14. II. Effects on the Forbush Decrease.** *arXiv preprint*. [arXiv:2010.07662](https://arxiv.org/abs/2010.07662). Role: 3-D GCR FTE of the Qin/Zhang family (Eq. 2); time-backward method. Reading: full text. Key in the supplied bibliography: **Qin-2021-AJ**. Note: published as The Astrophysical Journal 908(2), 236 (2021), DOI 10.3847/1538-4357/abd77c (confirmed by OpenAlex and the arXiv DOI field); record kept as the arXiv preprint.

**[P19]** He, H.-Q.; Qin, G.; Zhang, M. (2011). **Propagation of Solar Energetic Particles in Three-Dimensional Interplanetary Magnetic Fields: In View of Characteristics of Sources.** *The Astrophysical Journal*, **734**(2), 74. [DOI](https://doi.org/10.1088/0004-637X/734/2/74). Role: 3-D FTE of the Qin/Zhang family (garbled extraction); focusing term -(v/B)dB/dz. Reading: full text. Key in the supplied bibliography: **He-2011-AJ**.

**[P20]** He, H.-Q.; Wan, W. (2015). **Numerical Study of the Longitudinally Asymmetric Distribution of Solar Energetic Particles in the Heliosphere.** *The Astrophysical Journal Supplement Series*, **218**(2), 17. [DOI](https://doi.org/10.1088/0067-0049/218/2/17); [arXiv:1502.02683](https://arxiv.org/abs/1502.02683). Role: 3-D FTE Eqs. 1-3 with focusing-length sign as printed; time-backward SDEs. Reading: full text. Key in the supplied bibliography: **He-2015-AJ-218**.

**[P21]** He, H.-Q.; Wan, W. (2019). **Propagation of Solar Energetic Particles in the Outer Heliosphere: Interplay between Scattering and Adiabatic Focusing.** *The Astrophysical Journal Letters*, **885**(2), L28. [DOI](https://doi.org/10.3847/2041-8213/ab50bd); [arXiv:1911.02787](https://arxiv.org/abs/1911.02787). Role: repeats the He & Wan 2015 FTE form; time-backward SDEs. Reading: full text. Key in the supplied bibliography: **He-2019-AJL, He-2019-arXiv**.

**[P22]** He, H.-Q.; Zhou, G.; Wan, W. (2017). **Propagation of Solar Energetic Particles in Three-dimensional Interplanetary Magnetic Fields: Radial Dependence of Peak Intensities.** *arXiv preprint*. [arXiv:1706.05443](https://arxiv.org/abs/1706.05443). Role: repeats the He & Wan 2015 FTE form; time-backward SDEs (Euler). Reading: full text. Key in the supplied bibliography: **He-2017-AJ**. Note: published as He, H.-Q.; Zhou, G.; Wan, W. 2017, The Astrophysical Journal 842(2), 71, DOI 10.3847/1538-4357/aa7574.

**[P23]** Wang, Y.; Qin, G.; Zhang, M.; Dalla, S. (2014). **A Numerical Simulation of Solar Energetic Particle Dropouts during Impulsive Events.** *The Astrophysical Journal*, **789**(2), 157. [DOI](https://doi.org/10.1088/0004-637X/789/2/157); [arXiv:1403.6544](https://arxiv.org/abs/1403.6544). Role: 3-D FTE (Eq. 1); mu frame; time-backward method. Reading: full text. Key in the supplied bibliography: **Wang-2014-AJ**.

**[P24]** Lang, J. T.; Strauss, R. D.; Engelbrecht, N. E.; van den Berg, J. P.; Dresing, N.; Ruffolo, D.; Bandyopadhyay, R. (2024). **A Detailed Survey of the Parallel Mean Free Path of Solar Energetic Particle Protons and Electrons.** *The Astrophysical Journal*, **971**(1), 105. [DOI](https://doi.org/10.3847/1538-4357/ad55c3); [arXiv:2406.05765](https://arxiv.org/abs/2406.05765). Role: 1-D conservative FTE printing (Eq. 1). Reading: full text. Key in the supplied bibliography: **Lang-2024-arXiv**.

**[P25]** Wijsen, N.; Aran, A.; Pomoell, J.; Poedts, S. (2019). **Modelling three-dimensional transport of solar energetic protons in a corotating interaction region generated with EUHFORIA.** *Astronomy & Astrophysics*, **622**, A28. [DOI](https://doi.org/10.1051/0004-6361/201833958); [arXiv:1901.09596](https://arxiv.org/abs/1901.09596). Role: PARADISE mixed-frame FTE and coefficients (Eqs. 1-5). Reading: full text. Key in the supplied bibliography: **Wijsen-2019-AA-A28**.

**[P26]** Husidic, E.; Wijsen, N.; Baratashvili, T.; Poedts, S.; Vainio, R. (2024). **Energetic particle acceleration and transport with the novel Icarus + PARADISE model.** *Journal of Space Weather and Space Climate*, **14**, 11. [DOI](https://doi.org/10.1051/swsc/2024009); [arXiv:2404.12164](https://arxiv.org/abs/2404.12164). Role: PARADISE frame statement and forward-in-time Ito integration. Reading: full text. Key in the supplied bibliography: **Husidic-2024-JSWSC**.

**[P27]** Wang, Y.; Qin, G. (2015). **Estimation of the Release Time of Solar Energetic Particles Near the Sun.** *The Astrophysical Journal*, **799**(1), 111. [DOI](https://doi.org/10.1088/0004-637X/799/1/111); [arXiv:1311.7469](https://arxiv.org/abs/1311.7469). Role: frame statement for x, p, mu. Reading: full text. Key in the supplied bibliography: **Wang-2015-AJ**.

**[P28]** Strauss, R. D.; van den Berg, J. P.; Steyn, P. J.; Effenberger, F. J.; Wijsen, N.; Laitinen, T.; le Roux, J. A. (2020). **Perpendicular diffusion of solar energetic particles: When is the diffusion approximation valid?.** *Journal of Physics: Conference Series*, **1620**(1), 012021. [DOI](https://doi.org/10.1088/1742-6596/1620/1/012021). Role: conservative mixed-frame FTE with mixed perpendicular terms (Eq. 1, garbled); frame statement. Reading: full text. Key in the supplied bibliography: **Strauss-2020-JP**.

**[P29]** Afanasiev, A.; Wijsen, N.; Vainio, R. (2025). **Towards advanced forecasting of solar energetic particle events with the PARASOL model.** *Journal of Space Weather and Space Climate*, **15**, 3. [DOI](https://doi.org/10.1051/swsc/2024039); [arXiv:2412.11852](https://arxiv.org/abs/2412.11852). Role: PARADISE conservative intensity-form FTE (Eqs. 37-40); SOLPACS Eq. 12. Reading: full text. Key in the supplied bibliography: **Afanasiev-2025-JSWSC**.

**[P30]** Isenberg, P. A. (1997). **A hemispherical model of anisotropic interstellar pickup ions.** *Journal of Geophysical Research*, **102**(A3), 4719–4724. [DOI](https://doi.org/10.1029/96JA03671). Role: FTE source cited by Kartavykh et al. 2015; scattering in the solar-wind frame (abstract). Reading: abstract only. Pre-2006 foundation. Key in the supplied bibliography: **Isenberg-1997-JGR, isenberg97**.

**[P31]** Zuo, P.; Zhang, M.; Rassoul, H. K. (2013). **Acceleration of Low-Energy Ions at Parallel Shocks with a Focused Transport Model.** *The Astrophysical Journal*, **767**(1), 6. [DOI](https://doi.org/10.1088/0004-637X/767/1/6). Role: frame statement (p, mu, E in plasma frame; x fixed at shock); time-forward SDEs. Reading: full text. Key in the supplied bibliography: **Zuo-2013-AJ**.

**[P32]** Linker, J. A.; Caplan, R. M.; Schwadron, N.; Gorby, M.; Downs, C.; Török, T.; Lionello, R.; Wijaya, J. (2019). **Coupled MHD-Focused Transport Simulations for Modeling Solar Particle Events.** *arXiv preprint*. [arXiv:1905.05299](https://arxiv.org/abs/1905.05299). Role: EPREM FTE printing with relativistic prefactor. Reading: full text. Key in the supplied bibliography: **Linker-2019-JPCS**. Note: published as Linker, J. A.; Caplan, R. M.; Schwadron, N.; Gorby, M.; Downs, C.; Török, T.; Lionello, R.; Wijaya, J. 2019, Journal of Physics: Conference Series 1225(1), 012007, DOI 10.1088/1742-6596/1225/1/012007.

**[P33]** Kozarev, K. A.; Evans, R. M.; Schwadron, N. A.; Dayeh, M. A.; Opher, M.; Korreck, K. E.; van der Holst, B. (2013). **Global Numerical Modeling of Energetic Proton Acceleration in a Coronal Mass Ejection Traveling through the Solar Corona.** *The Astrophysical Journal*, **778**(1), 43. [DOI](https://doi.org/10.1088/0004-637X/778/1/43); [arXiv:1406.2377](https://arxiv.org/abs/1406.2377). Role: EPREM FTE with (1 - V.e_b v mu/c^2) prefactor and D_mumu/2 operator. Reading: full text. Key in the supplied bibliography: **Kozarev2013**.

**[P34]** Schwadron, N. A.; Townsend, L.; Kozarev, K.; Dayeh, M. A.; Cucinotta, F.; Desai, M.; Golightly, M.; Hassler, D.; et al. (2010). **Earth-Moon-Mars Radiation Environment Module framework.** *Space Weather*, **8**(1), S00E02. [DOI](https://doi.org/10.1029/2009SW000523). Role: EPREM scattering operator (D_mumu/2) (para. 27, partly garbled). Reading: full text. Key in the supplied bibliography: **Schwadron-2010-SW**.

**[P35]** Kóta, J.; Manchester, W. B.; Gombosi, T. I. (2005). **SEP Acceleration at Realistic CMEs: Two Sites of Acceleration?.** *Proceedings of the 29th International Cosmic Ray Conference (Pune)*, **1**, 125–128. Role: foundation of the EPREM co-moving (Lagrangian) FTE; equation garbled. Reading: full text. Pre-2006 foundation. Proposed key (additions file): **Kota-2005-ICRC**.

**[P36]** Borovikov, D.; Sokolov, I. V.; Huang, Z.; Roussev, I. I.; Gombosi, T. I. (2019). **Toward quantitative model for simulation and forecast of solar energetic particle production during gradual events - II: kinetic description of SEP.** *arXiv preprint*. [arXiv:1911.10165](https://arxiv.org/abs/1911.10165). Role: M-FLAMPA co-moving frame and per-branch wave-frame scattering operator (Eq. 4.29). Reading: full text. Key in the supplied bibliography: **Borovikov-2019-arXiv**.

**[P37]** Afanasiev, A.; Vainio, R.; Trotta, D.; Nyberg, S.; Talebpour Sheshvan, N.; Hietala, H.; Dresing, N. (2023). **Self-consistent modeling of the energetic storm particle event of November 10, 2012.** *Astronomy & Astrophysics*, **679**, A111. [DOI](https://doi.org/10.1051/0004-6361/202346220); [arXiv:2310.00993](https://arxiv.org/abs/2310.00993). Role: SOLPACS wave-frame transport equation (Eq. 1). Reading: full text. Key in the supplied bibliography: **Afanasiev-2023-AA**.

**[P38]** Ng, C. K.; Reames, D. V. (2008). **Shock Acceleration of Solar Energetic Protons: The First 10 Minutes.** *The Astrophysical Journal Letters*, **686**(2), L123–L126. [DOI](https://doi.org/10.1086/592996). Role: outward-wave-frame transport coefficients (Eqs. 1-4). Reading: full text. Key in the supplied bibliography: **Ng-2008-AJL**.

**[P39]** Zhang, M.; Zhao, L. (2017). **Precipitation and Release of Solar Energetic Particles from the Solar Coronal Magnetic Field.** *The Astrophysical Journal*, **846**(2), 107. [DOI](https://doi.org/10.3847/1538-4357/aa86a8). Role: corotating spatial frame; time-forward SDEs. Reading: full text. Key in the supplied bibliography: **Zhang-2017-AJ**.

**[P40]** Zhang, M.; Qin, G.; Rassoul, H. K. (2009). **Propagation of Solar Energetic Particles in Three-Dimensional Interplanetary Magnetic Fields.** *The Astrophysical Journal*, **692**(1), 109–132. [DOI](https://doi.org/10.1088/0004-637X/692/1/109). Role: time-backward stochastic solution in a fixed frame (abstract). Reading: abstract only. Key in the supplied bibliography: **Zhang-2009-AJ**.

**[P41]** Qin, G.; Zhang, M.; Dwyer, J. R. (2006). **Effect of adiabatic cooling on the fitted parallel mean free path of solar energetic particles.** *Journal of Geophysical Research: Space Physics*, **111**(A8), A08101. [DOI](https://doi.org/10.1029/2005JA011512). Role: adiabatic cooling added to the FTE; time-backward stochastic integration (abstract). Reading: abstract only. Key in the supplied bibliography: **Qin-2006-JGR**. Metadata not fully confirmed: article number A08101 not confirmed by a primary source (OpenAlex gives none; it is taken from the Kubo et al. 2015 reference list); all other fields confirmed.

**[P42]** Ding, Z.; Wijsen, N.; Li, G.; Poedts, S. (2022). **Modeling the 2020 November 29 solar energetic particle event using EUHFORIA and iPATH models.** *Astronomy & Astrophysics*, **668**, A71. [DOI](https://doi.org/10.1051/0004-6361/202244732); [arXiv:2210.16967](https://arxiv.org/abs/2210.16967). Role: iPATH time-backward SDE solver (direction of integration). Reading: full text. Key in the supplied bibliography: **Ding-2022-AA, Ding-2022-arXiv**.

**[P43]** Dröge, W.; Kartavykh, Y. Y.; Klecker, B.; Kovaltsov, G. A. (2010). **Anisotropic Three-Dimensional Focused Transport of Solar Energetic Particles in the Inner Heliosphere.** *The Astrophysical Journal*, **709**(2), 912–919. [DOI](https://doi.org/10.1088/0004-637X/709/2/912). Role: 3-D focused transport with SDEs (abstract; direction not stated). Reading: abstract only. Key in the supplied bibliography: **Droge-2010-AJ**.

**[P44]** Strauss, R. D.; Effenberger, F. (2017). **A Hitch-hiker's Guide to Stochastic Differential Equations: Solution Methods for Energetic Particle Transport in Space Physics and Astrophysics.** *Space Science Reviews*, **212**(1), 151–192. [DOI](https://doi.org/10.1007/s11214-017-0351-y); [arXiv:1703.06192](https://arxiv.org/abs/1703.06192). Role: statement on forward versus backward SDE integration. Reading: full text. Key in the supplied bibliography: **Strauss-2017-SSR**.

**[P45]** Riordan, J. D.; Pe'er, A. (2019). **Pitch-angle Diffusion and Bohm-type Approximations in Diffusive Shock Acceleration.** *The Astrophysical Journal*, **873**(1), 13. [DOI](https://doi.org/10.3847/1538-4357/aaffd2); [arXiv:1810.11817](https://arxiv.org/abs/1810.11817). Role: MSD estimator without ½ (Eq. 11), Δt = 20 t_g window, plateau artefact. Reading: full text. Proposed key (additions file): **Riordan-2019-AJ**.

**[P46]** Florinski, V. (2013). **Entering the "magnetic highway": energetic particle anisotropies at the heliospheric boundary.** *arXiv preprint*. [arXiv:1302.2179](https://arxiv.org/abs/1302.2179). Role: gyrophase-resolved heliopause model; random walk on a sphere; λ∥ = v/(2D). Reading: full text. Proposed key (additions file): **Florinski-2013-arXiv**. Note: single-author preprint (posted 2013-02-08) from which the equations were read; a four-author Letter with the shorter title "Energetic particle anisotropies at the heliospheric boundary" (Florinski, Jokipii, Alouani-Bibi & le Roux, ApJL 776, L37, 2013, DOI 10.1088/2041-8205/776/2/L37) reports a pancake rather than a loss-cone distribution in its abstract; its equations were not read.

**[P47]** Chen, X.; Giacalone, J.; Guo, F.; Klein, K. G. (2024). **Parallel Diffusion Coefficient of Energetic Charged Particles in the Inner Heliosphere from the Turbulent Magnetic Fields Measured by Parker Solar Probe.** *The Astrophysical Journal*, **965**(1), 61. [DOI](https://doi.org/10.3847/1538-4357/ad33c3); [arXiv:2403.08141](https://arxiv.org/abs/2403.08141). Role: D_μμ from PSP frequency spectra (Eq. 3), μ_min = 0.05. Reading: full text. Key in the supplied bibliography: **Chen-2024-AJ-965**.

**[P48]** Agueda, N.; Vainio, R.; Lario, D.; Sanahuja, B. (2010). **Solar near-relativistic electron observations as a proof of a back-scatter region beyond 1 AU during the 2000 February 18 event.** *Astronomy & Astrophysics*, **519**, A36. [DOI](https://doi.org/10.1051/0004-6361/200913963). Role: ε-form Eq. 1, λ∥ = (v/ν₀)φ(ε) and asymptote Eqs. 3–4; ε = 0.01; ν(1) = 3.5v/λ. Reading: full text. Proposed key (additions file): **Agueda-2010-AA**.

**[P49]** Agueda, N.; Lario, D.; Vainio, R.; Sanahuja, B.; Kilpua, E.; Pohjolainen, S. (2009). **Modeling solar near-relativistic electron events: Insights into solar injection and interplanetary transport conditions.** *Astronomy & Astrophysics*, **507**(2), 981–993. [DOI](https://doi.org/10.1051/0004-6361/200912224). Role: ε-form in the solar-wind frame; ε = 0.10, 0.01 and isotropic fits. Reading: full text. Proposed key (additions file): **Agueda-2009-AA**.

**[P50]** Beeck, J.; Wibberenz, G. (1986). **Pitch angle distributions of solar energetic particles and the local scattering properties of the interplanetary medium.** *The Astrophysical Journal*, **311**, 437. [DOI](https://doi.org/10.1086/164784). Role: attributed origin of the q-form with H (not read). Reading: metadata only. Pre-2006 foundation. Key in the supplied bibliography: **Beeck-1986-AJ**. Note: pages field holds the first page only (437), as in both sources; the last page was not available (IOP and ADS refuse automated fetches).

**[P51]** Strauss, R. D.; Ogunjobi, O.; Moraal, H.; McCracken, K. G.; Caballero-Lopez, R. A. (2017). **On the Pulse Shape of Ground-Level Enhancements.** *Solar Physics*, **292**(4), 51. [DOI](https://doi.org/10.1007/s11207-017-1086-3); [arXiv:1703.05906](https://arxiv.org/abs/1703.05906). Role: q-form D(z)(1−μ²){|μ|^{q-1}+H} (Eq. 8), q = 5/3, H = 0.05. Reading: full text. Key in the supplied bibliography: **Strauss-2017-SP**.

**[P52]** Agueda, N.; Vainio, R. (2013). **On the parametrization of the energetic-particle pitch-angle diffusion coefficient.** *Journal of Space Weather and Space Climate*, **3**, A10. [DOI](https://doi.org/10.1051/swsc/2013034). Role: q-form with ν₀ normalization (Eqs. 2–4), ε-form (Eq. 5), ε–(q,H) least-squares equivalences. Reading: full text. Proposed key (additions file): **Agueda-2013-JSWSC**.

**[P53]** Qin, G.; Wang, Y. (2015). **Simulations of a Gradual Solar Energetic Particle Event Observed by Helios 1, Helios 2, and IMP 8.** *The Astrophysical Journal*, **809**(2), 177. [DOI](https://doi.org/10.1088/0004-637X/809/2/177); [arXiv:1505.02974](https://arxiv.org/abs/1505.02974). Role: slab-QLT amplitude D₀v(R_L k_min)^{s-2}, h = 0.01, k_min = 32 AU⁻¹ (Eqs. 2–3). Reading: full text. Key in the supplied bibliography: **Qin-2015-AJ**.

**[P54]** Lavasa, E.; Lang, J. T.; Papaioannou, A.; Strauss, R. D.; Mallios, S. A.; Hillaris, A.; Kouloumvakos, A.; Anastasiadis, A.; et al. (2026). **Multi-spacecraft constraints on relativistic solar energetic particle transport in the widespread 28 October 2021 event.** *Astronomy & Astrophysics*, **707**, A12. [DOI](https://doi.org/10.1051/0004-6361/202558094); [arXiv:2603.09839](https://arxiv.org/abs/2603.09839). Role: printed (1−μ)² variant of the q-form (Eq. 3), q = 5/3, H = 0.05. Reading: full text. Key in the supplied bibliography: **Lavasa-2026-AA**.

**[P55]** Àgueda Costafreda, N. (2008). **Near-relativistic electron events. Monte Carlo simulations of solar injection and interplanetary transport.** *PhD thesis, Universitat de Barcelona (hdl:2445/35265)*. Role: ε-form origin (Eqs. 3.19, 3.21), H = 0 normalization (Eq. 3.8), hemispheric η transformation (Eq. 3.20). Reading: full text. Proposed key (additions file): **AguedaCostafreda-2008-PHD**.

**[P56]** He, H.-Q. (2015). **Perpendicular Diffusion in the Transport of Solar Energetic Particles from Unconnected Sources: The Counter-streaming Particle Beams Revisited.** *The Astrophysical Journal*, **814**(2), 157. [DOI](https://doi.org/10.1088/0004-637X/814/2/157); [arXiv:1512.00027](https://arxiv.org/abs/1512.00027). Role: D^r_μμ rigidity form with h = 0.2 for electrons. Reading: full text. Key in the supplied bibliography: **He-2015-AJ-814**.

**[P57]** Dröge, W.; Kartavykh, Y. Y. (2008). **Testing transport theories with solar energetic particles.** *Proceedings of the 30th International Cosmic Ray Conference (Mérida, 2007)*, **1 (SH)**, 159–162. Role: D^r_μμ = κ₀(r,R){|μ|^{q-1}+H}(1−μ²) with q = 1.67, H = 0.05 (Eq. 2). Reading: full text. Proposed key (additions file): **Droge-2008-ICRC**.

**[P58]** Wu, S.-S.; Qin, G. (2020). **Magnetic Cloud and Sheath in the Ground-Level Enhancement Event of 2000 July 14. I. Effects on the Solar Energetic Particles.** *arXiv preprint*. [arXiv:2008.06427](https://arxiv.org/abs/2008.06427). Role: same D_μμ definition as QinWu2020 (μ^{s-1}). Reading: full text. Key in the supplied bibliography: **Wu-2020-AJ**. Note: published as Wu, S.-S.; Qin, G. 2020, The Astrophysical Journal 904(2), 151, DOI 10.3847/1538-4357/abc0f2 (confirmed by OpenAlex and the arXiv DOI field); record kept as the arXiv preprint. The former container wording '(accepted draft, The Astrophysical Journal)' is not supported by the arXiv page (no comments line), hence container set to 'arXiv preprint'.

**[P59]** Teufel, A.; Schlickeiser, R. (2003). **Analytic calculation of the parallel mean free path of heliospheric cosmic rays. II. Dynamical magnetic slab turbulence and random sweeping slab turbulence with finite wave power at small wavenumbers.** *Astronomy & Astrophysics*, **397**(1), 15–25. [DOI](https://doi.org/10.1051/0004-6361:20021471). Role: slab λ∥ form whose R_L k_min ≪ 1 limit equals Eq. (28) at h = 0 (via MFP v1.2 Sect. 6.4). Reading: full text. Pre-2006 foundation. Proposed key (additions file): **Teufel-2003-AA**.

**[P60]** Qin, G.; Wang, Y.; Zhang, M.; Dalla, S. (2013). **Transport of Solar Energetic Particles Accelerated by ICME Shocks: Reproducing the Reservoir Phenomenon.** *The Astrophysical Journal*, **766**(2), 74. [DOI](https://doi.org/10.1088/0004-637X/766/2/74). Role: h = 0.01, k_min = 33 AU⁻¹ (AMPS Type V attribution). Reading: full text. Key in the supplied bibliography: **Qin-2013-AJ**.

**[P61]** Agueda, N.; Vainio, R.; Lario, D.; Sanahuja, B. (2008). **Injection and Interplanetary Transport of Near-Relativistic Electrons: Modeling the Impulsive Event on 2000 May 1.** *The Astrophysical Journal*, **675**(2), 1601–1613. [DOI](https://doi.org/10.1086/527527). Role: origin of the ε-form (journal paper not read). Reading: secondary. Proposed key (additions file): **Agueda-2008-AJ**.

**[P62]** Wijsen, N.; Aran, A.; Scolini, C.; Lario, D.; Afanasiev, A.; Vainio, R.; Sanahuja, B.; Pomoell, J.; et al. (2022). **Observation-based modelling of the energetic storm particle event of 14 July 2012.** *Astronomy & Astrophysics*, **659**, A187. [DOI](https://doi.org/10.1051/0004-6361/202142698); [arXiv:2201.06454](https://arxiv.org/abs/2201.06454). Role: ε-form written with D₀ (Eq. 2), ε = 0.048. Reading: full text. Key in the supplied bibliography: **Wijsen-2022-AA**.

**[P63]** Pacheco, D.; Agueda, N.; Aran, A.; Heber, B.; Lario, D. (2019). **Full inversion of solar relativistic electron events measured by the Helios spacecraft.** *Astronomy & Astrophysics*, **624**, A3. [DOI](https://doi.org/10.1051/0004-6361/201834520); [arXiv:1902.06602](https://arxiv.org/abs/1902.06602). Role: printed |μ|/(1−|μ|) variant of the ε-form (Eq. 4), ε = 0.01. Reading: full text. Key in the supplied bibliography: **Pacheco-2019-AA**.

**[P64]** Minoshima, T.; Miyoshi, Y.; Murakami, G.; Pinto, M.; Schmid, D.; Matsuoka, A.; Baumjohann, W.; Fischer, D.; et al. (2026). **A study of solar energetic particle transport on 30 March 2022 using multi-spacecraft data assimilation.** *Earth, Planets and Space*, **78**(1), 64. [DOI](https://doi.org/10.1186/s40623-026-02389-9); [arXiv:2602.00765](https://arxiv.org/abs/2602.00765). Role: Dröge (2003) V_A form as printed (Eq. 8), V_A = 100 km/s. Reading: full text. Key in the supplied bibliography: **Minoshima-2026-EPS**.

**[P65]** Dröge, W. (2003). **Solar Particle Transport in a Dynamical Quasi-linear Theory.** *The Astrophysical Journal*, **589**(2), 1027–1039. [DOI](https://doi.org/10.1086/374812). Role: attributed origin of the V_A shape (secondary via Minoshima 2026). Reading: abstract only. Pre-2006 foundation. Key in the supplied bibliography: **Droge-2003-AJ**.

**[P66]** Liu, W.; Sokolov, I. V.; Zhao, L.; Gombosi, T. I.; Sachdeva, N.; Chen, X.; Tóth, G.; Lario, D.; et al. (2025). **Physics-based Simulation of the 2013 April 11 Solar Energetic Particle Event.** *The Astrophysical Journal*, **985**(1), 82. [DOI](https://doi.org/10.3847/1538-4357/adc4e3); [arXiv:2412.07581](https://arxiv.org/abs/2412.07581). Role: M-FLAMPA solves the Parker equation; shape used for D∥. Reading: full text. Key in the supplied bibliography: **Liu-2025-AJ**.

**[P67]** le Roux, J. A.; Webb, G. M. (2009). **Time-Dependent Acceleration of Interstellar Pickup Ions at the Heliospheric Termination Shock Using a Focused Transport Approach.** *The Astrophysical Journal*, **693**(1), 534–551. [DOI](https://doi.org/10.1088/0004-637X/693/1/534). Role: attributed source of AMPS Type II (secondary via Tenishev 2022). Reading: abstract only. Key in the supplied bibliography: **Roux-2009-AJ**.

**[P68]** Hu, J.; Li, G.; Ao, X.; Zank, G. P.; Verkhoglyadova, O. (2017). **Modeling Particle Acceleration and Transport at a 2-D CME-Driven Shock.** *Journal of Geophysical Research: Space Physics*, **122**(11). [DOI](https://doi.org/10.1002/2017JA024077). Role: attributed source of AMPS Type III (secondary via Tenishev 2022). Reading: secondary. Key in the supplied bibliography: **Hu-2017-JGR**. Metadata not fully confirmed: page range not confirmed: neither OpenAlex nor NTRS gives it and the Wiley pages could not be read.

**[P69]** Young, M. A.; Schwadron, N. A.; Gorby, M.; Linker, J.; Caplan, R. M.; Downs, C.; Török, T.; Riley, P.; et al. (2021). **Energetic Proton Propagation and Acceleration Simulated for the Bastille Day Event of 2000 July 14.** *The Astrophysical Journal*, **909**(2), 160. [DOI](https://doi.org/10.3847/1538-4357/abdf5f); [arXiv:2012.09078](https://arxiv.org/abs/2012.09078). Role: EPREM operator with D_μμ/2 (Eq. 1). Reading: full text. Key in the supplied bibliography: **Young-2021-AJ**.

**[P70]** Bieber, J. W.; Matthaeus, W. H.; Smith, C. W.; Wanner, W.; Kallenrode, M.-B.; Wibberenz, G. (1994). **Proton and Electron Mean Free Paths: The Palmer Consensus Revisited.** *The Astrophysical Journal*, **420**(1), 294–306. [DOI](https://doi.org/10.1086/173559). Role: cited by Lang et al. (2024) for the q-form (not read). Reading: metadata only. Pre-2006 foundation. Proposed key (additions file): **Bieber-1994-AJ-420**. Metadata not fully confirmed: last page (306) not confirmed: OpenAlex and NTRS give only the first page 294 or none; issue 1 from the NTRS record only.

**[P71]** Dröge, W. (2000). **The Rigidity Dependence of Solar Particle Scattering Mean Free Paths.** *The Astrophysical Journal*, **537**(2), 1073–1079. [DOI](https://doi.org/10.1086/309080). Role: cited by Lang et al. (2024) and Agueda et al. (2010) for the q-form / finite scattering through μ = 0 (not read). Reading: metadata only. Pre-2006 foundation. Proposed key (additions file): **Droge-2000-AJ**.

**[P72]** Teufel, A.; Schlickeiser, R. (2002). **Analytic calculation of the parallel mean free path of heliospheric cosmic rays. I. Dynamical magnetic slab turbulence and random sweeping slab turbulence.** *Astronomy & Astrophysics*, **393**(2), 703–715. [DOI](https://doi.org/10.1051/0004-6361:20021046). Role: DT/RS slab D_μμ (Eqs. 34-36, 41-54), 8π normalization, σ tensor. Reading: full text. Pre-2006 foundation. Proposed key (additions file): **Teufel-2002-AA**.

**[P73]** Shalchi, A.; Škoda, T.; Tautz, R. C.; Schlickeiser, R. (2009). **Analytical description of nonlinear cosmic ray scattering: isotropic and quasilinear regimes of pitch-angle diffusion.** *Astronomy & Astrophysics*, **507**(2), 589–597. [DOI](https://doi.org/10.1051/0004-6361/200912755). Role: Slab QLT D_μμ (Eqs. 6-8, 29), normalized spectrum D(s,q) (Eqs. 34-35), QLT validity |μ| > δB/B₀. Reading: full text. Key in the supplied bibliography: **Shalchi-2009-AA**.

**[P74]** Engelbrecht, N. E.; Effenberger, F.; Florinski, V.; Potgieter, M. S.; Ruffolo, D.; Chhiber, R.; Usmanov, A. V.; Rankin, J. S.; et al. (2022). **Theory of Cosmic Ray Transport in the Heliosphere.** *Space Science Reviews*, **218**(4), 33. [DOI](https://doi.org/10.1007/s11214-022-00896-1). Role: Slab QLT D_μμ with S_xx (Eq. 11); electron λ∥ forms. Reading: full text. Key in the supplied bibliography: **Engelbrecht-2022-SSRv**.

**[P75]** Jokipii, J. R. (1966). **Cosmic-ray propagation. I. Charged particles in a random magnetic field.** *The Astrophysical Journal*, **146**, 480. [DOI](https://doi.org/10.1086/148912). Role: Original magnetostatic slab QLT (secondary via Engelbrecht et al. 2022). Reading: secondary. Pre-2006 foundation. Key in the supplied bibliography: **Jokipii-1966-AJ**.

**[P76]** Shalchi, A.; Lazarian, A.; Schlickeiser, R. (2007). **Non-linear damping of slab modes and cosmic ray transport.** *Monthly Notices of the Royal Astronomical Society*, **383**(2), 803–808. [DOI](https://doi.org/10.1111/j.1365-2966.2007.12590.x); [arXiv:0710.5418](https://arxiv.org/abs/0710.5418). Role: Nonlinear-damping slab D_μμ (Eqs. 17-24), bend-over spectrum. Reading: full text. Key in the supplied bibliography: **Shalchi-2008-MNRAS**. Note: year: the publisher gives December 2007 (OpenAlex 2007); the tag Shalchi2008 reflects a 2008 citation that was not confirmed in any source read.

**[P77]** Tautz, R. C.; Shalchi, A. (2013). **Simulated energetic particle transport in the interplanetary space: The Palmer consensus revisited.** *Journal of Geophysical Research: Space Physics*, **118**(2), 642–647. [DOI](https://doi.org/10.1002/jgra.50155); [arXiv:1301.7162](https://arxiv.org/abs/1301.7162). Role: Bend-over slab spectrum with C(s) (Eq. 4). Reading: full text. Key in the supplied bibliography: **Tautz-2013-JGR, Tautz-2013-JGRSP**.

**[P78]** Afanasiev, A.; Battarbee, M.; Vainio, R. (2015). **Self-consistent Monte Carlo simulations of proton acceleration in coronal shocks: Effect of anisotropic pitch-angle scattering of particles.** *Astronomy & Astrophysics*, **584**, A81. [DOI](https://doi.org/10.1051/0004-6361/201526750); [arXiv:1603.08857](https://arxiv.org/abs/1603.08857). Role: Wave-frame resonant D_μμ (Eq. 2), I_w normalization (Eq. 5), SOLPACS grid. Reading: full text. Key in the supplied bibliography: **Afanasiev-2015-AA**.

**[P79]** Shalchi, A.; Bieber, J. W.; Matthaeus, W. H.; Schlickeiser, R. (2006). **Parallel and Perpendicular Transport of Heliospheric Cosmic Rays in an Improved Dynamical Turbulence Model.** *The Astrophysical Journal*, **642**(1), 230–243. [DOI](https://doi.org/10.1086/500728). Role: NADT slab and 2D D_μμ (structure only, OCR). Reading: full text. Proposed key (additions file): **Shalchi-2006-AJ**.

**[P80]** Spanier, F.; Schreiner, C.; Schlickeiser, R. (2022). **Determining Pitch-Angle Diffusion Coefficients for Electrons in Whistler Turbulence.** *Physics*, **4**(1), 80–103. [DOI](https://doi.org/10.3390/physics4010008); [arXiv:2112.03999](https://arxiv.org/abs/2112.03999). Role: whistler/Alfvén D_μμ with cross helicity H_{L,R} and magnetic helicity σ; cold e–p dispersion. Reading: full text. Key in the supplied bibliography: **Spanier-2022-Physics**.

**[P81]** Schlickeiser, R. (1989). **Cosmic-Ray Transport and Acceleration. I. Derivation of the Kinetic Equation and Application to Cosmic Rays in Static Cold Media.** *The Astrophysical Journal*, **336**, 243–263. [DOI](https://doi.org/10.1086/167009). Role: cited origin of helicity weights and of A_2 (not read). Reading: metadata only. Pre-2006 foundation. Key in the supplied bibliography: **Schlickeiser-1989-AJ**. Metadata not fully confirmed: pages: last page 263 not confirmed (OpenAlex gives first page 243 only; ADS not reachable).

**[P82]** Ng, C. K.; Reames, D. V.; Tylka, A. J. (2012). **Solar Energetic Particles: Shock Acceleration and Transport through Self-Amplified Waves.** *AIP Conference Proceedings*, **1436**, 212–218. [DOI](https://doi.org/10.1063/1.4723610). Role: inward/outward wave initialization I_R− = I_L− = 0.1 I_R+. Reading: full text. Key in the supplied bibliography: **Ng-2012-AIP**. Metadata not fully confirmed: volume (1436) not confirmed: OpenAlex gives no volume and the AIP pages could not be read; all other fields confirmed.

**[P83]** Sandroos, A.; Li, G.; Zhao, L. (2014). **First- and Second-order Fermi Acceleration at Parallel Shocks.** *arXiv preprint*. [arXiv:1408.3022](https://arxiv.org/abs/1408.3022). Role: D_μμ in parallel/antiparallel wave rest frames (Eqs. 29–31); L/R helicity labels. Reading: full text. Proposed key (additions file): **Sandroos-2014-arXiv**. Note: no journal publication found (no journal-ref on arXiv; OpenAlex lists only the arXiv preprint).

**[P84]** Chandran, B. D. G. (2000). **Scattering of energetic particles by anisotropic magnetohydrodynamic turbulence with a Goldreich-Sridhar power spectrum.** *Physical Review Letters*, **85**(22), 4656–4659. [DOI](https://doi.org/10.1103/PhysRevLett.85.4656); [arXiv:astro-ph/0008498](https://arxiv.org/abs/astro-ph/0008498). Role: GS95 shear/pseudo-Alfvén D_ξξ and κ∥ (Eqs. 6–17). Reading: full text. Pre-2006 foundation. Key in the supplied bibliography: **Chandran-2000-PRL**.

**[P85]** Yan, H.; Lazarian, A. (2002). **Scattering of Cosmic Rays by Magnetohydrodynamic Interstellar Turbulence.** *Physical Review Letters*, **89**(28), 281102. [DOI](https://doi.org/10.1103/PhysRevLett.89.281102); [arXiv:astro-ph/0205285](https://arxiv.org/abs/astro-ph/0205285). Role: Alfvénic gyroresonance approximation (Eq. 8); fast-mode tensors (Eqs. 4–5). Reading: full text. Pre-2006 foundation. Key in the supplied bibliography: **Yan-2002-PRL**.

**[P86]** Vainio, R. (2000). **Charged-Particle Resonance Conditions and Transport Coefficients in Slab-Mode Waves.** *The Astrophysical Journal Supplement Series*, **131**(2), 519–529. [DOI](https://doi.org/10.1086/317372). Role: abstract statements on cross helicity and fast-mode mixtures. Reading: abstract only. Pre-2006 foundation. Key in the supplied bibliography: **Vainio-2000-AJSS**.

**[P87]** Vainio, R. (2003). **On the generation of Alfvén waves by solar energetic particles.** *Astronomy & Astrophysics*, **406**(2), 735–740. [DOI](https://doi.org/10.1051/0004-6361:20030822). Role: wave-frame D_μμ (Eq. 2); ΔE_w = −V_A pΔμ. Reading: full text. Pre-2006 foundation. Key in the supplied bibliography: **Vainio-2003-AA-406**.

**[P88]** Vainio, R.; Kocharov, L. (2001). **Proton transport through self-generated waves in impulsive flares.** *Astronomy & Astrophysics*, **375**(1), 251–259. [DOI](https://doi.org/10.1051/0004-6361:20010846). Role: isotropic scattering in each ± wave frame; −ΔE = ∓V_A p_± Δμ_±. Reading: full text. Pre-2006 foundation. Proposed key (additions file): **Vainio-2001-AA**.

**[P89]** Battarbee, M.; Laitinen, T.; Vainio, R. (2011). **Heavy-ion acceleration and self-generated waves in coronal shocks.** *Astronomy & Astrophysics*, **535**, A34. [DOI](https://doi.org/10.1051/0004-6361/201117507); [arXiv:1303.4340](https://arxiv.org/abs/1303.4340). Role: CSA resonance and ν (Eqs. 1–2); ΔE_w (Eq. 7). Reading: full text. Key in the supplied bibliography: **Battarbee-2013-AA-A34**.

**[P90]** Schreiner, C.; Spanier, F. (2015). **Effects of dispersive wave modes on charged particles transport.** *Proceedings of the 34th International Cosmic Ray Conference (ICRC2015), The Hague, PoS(ICRC2015)*, **236**, 177. [DOI](https://doi.org/10.22323/1.236.0177); [arXiv:1512.04692](https://arxiv.org/abs/1512.04692). Role: wave-frame μ for dispersive waves (Eqs. 2.4, 2.9–2.11). Reading: full text. Proposed key (additions file): **Schreiner-2015-ICRC**. Note: year given as the conference year (ICRC2015); the PoS volume 236 was published 2016-08-18.

**[P91]** Virtanen, J. J. P.; Vainio, R. (2005). **Stochastic acceleration in relativistic parallel shocks.** *The Astrophysical Journal*, **621**(1), 313–323. [DOI](https://doi.org/10.1086/427324); [arXiv:astro-ph/0411184](https://arxiv.org/abs/astro-ph/0411184). Role: sequential forward/backward wave-frame scattering; momentum diffusion A_2 (Eq. 15). Reading: full text. Pre-2006 foundation. Proposed key (additions file): **Virtanen-2005-AJ**.

**[P92]** Vainio, R.; Laitinen, T.; Fichtner, H. (2003). **A simple analytical expression for the power spectrum of cascading Alfvén waves in the solar wind.** *Astronomy & Astrophysics*, **407**(2), 713–723. [DOI](https://doi.org/10.1051/0004-6361:20030914). Role: D_μμ from a fixed-frame frequency spectrum, f_r = (Ω/2π)V/(v|μ|). Reading: full text. Pre-2006 foundation. Proposed key (additions file): **Vainio-2003-AA**.

**[P93]** Ng, C. K.; Reames, D. V.; Tylka, A. J. (2003). **Modeling shock-accelerated solar energetic particles coupled to interplanetary Alfvén waves.** *The Astrophysical Journal*, **591**(1), 461–485. [DOI](https://doi.org/10.1086/375293). Role: legible factors of the resonance function (Eq. 7); intensity normalization (Eq. 18). Reading: full text. Pre-2006 foundation. Key in the supplied bibliography: **Ng-2003-AJ**.

**[P94]** Yan, H.; Lazarian, A. (2004). **Cosmic-Ray Scattering and Streaming in Compressible Magnetohydrodynamic Turbulence.** *The Astrophysical Journal*, **614**, 757–769. [DOI](https://doi.org/10.1086/423733); [arXiv:astro-ph/0408172](https://arxiv.org/abs/astro-ph/0408172). Role: fast-mode spectra and gyroresonance (Eqs. 9–15). Reading: full text. Pre-2006 foundation. Key in the supplied bibliography: **Yan-2004-AJ**. Note: title hyphenation 'Cosmic-Ray' is supported only by a truncated INSPIRE title entry ('COSMIC-RAY'); the arXiv title has 'Cosmic Ray'. Issue number not available from the sources read (record leaves it blank).

**[P95]** Xu, S.; Lazarian, A. (2020). **Trapping of Cosmic Rays in MHD Turbulence.** *The Astrophysical Journal*, **894**(1), 63. [DOI](https://doi.org/10.3847/1538-4357/ab8465); [arXiv:2003.12593](https://arxiv.org/abs/2003.12593). Role: QLT gyroresonance approximations; mirror-trapping cutoff μ_cf (Eqs. 10–44). Reading: full text. Key in the supplied bibliography: **Xu-2020-AJ**.

**[P96]** Xu, S.; Lazarian, A. (2018). **Resonance-broadened Transit Time Damping of Particles in MHD Turbulence.** *The Astrophysical Journal*, **868**(1), 36. [DOI](https://doi.org/10.3847/1538-4357/aae840); [arXiv:1810.07726](https://arxiv.org/abs/1810.07726). Role: broadened TTD for slow and fast modes; fast gyroresonance (Eqs. 1–75). Reading: full text. Proposed key (additions file): **Xu-2018-AJ**.

**[P97]** Yan, H.; Lazarian, A. (2008). **Cosmic-Ray Propagation: Nonlinear Diffusion Parallel and Perpendicular to Mean Magnetic Field.** *The Astrophysical Journal*, **673**(2), 942–953. [DOI](https://doi.org/10.1086/524771); [arXiv:0710.2617](https://arxiv.org/abs/0710.2617). Role: Gaussian resonance broadening; broadened TTD and gyroresonance (Eqs. 1–22). Reading: full text. Key in the supplied bibliography: **Yan-2008-AJ**.

**[P98]** Xu, S.; Yan, H. (2013). **Cosmic-Ray Parallel and Perpendicular Transport in Turbulent Magnetic Fields.** *The Astrophysical Journal*, **779**(2), 140. [DOI](https://doi.org/10.1088/0004-637X/779/2/140); [arXiv:1307.1346](https://arxiv.org/abs/1307.1346). Role: test-particle D_μμ estimator in compressible MHD turbulence. Reading: full text. Key in the supplied bibliography: **Xu-2013-AJ**.

**[P99]** Lynn, J. W.; Parrish, I. J.; Quataert, E.; Chandran, B. D. G. (2012). **Resonance Broadening and Heating of Charged Particles in Magnetohydrodynamic Turbulence.** *The Astrophysical Journal*, **758**(2), 78. [DOI](https://doi.org/10.1088/0004-637X/758/2/78); [arXiv:1204.0155](https://arxiv.org/abs/1204.0155). Role: simulation-calibrated TTD resonance functions (D in v∥). Reading: full text. Key in the supplied bibliography: **Lynn-2012-AJ**.

**[P100]** Shalchi, A. (2005). **Second-order quasilinear theory of cosmic ray transport.** *Physics of Plasmas*, **12**(5). [DOI](https://doi.org/10.1063/1.1895805). Role: origin of the SOQLT resonance function (secondary via Shalchi et al. 2009). Reading: secondary. Pre-2006 foundation. Proposed key (additions file): **Shalchi-2005-PoP**. Metadata not fully confirmed: pages: article number (052905 per the task brief) not confirmed; OpenAlex has no page/article number and the AIP and ADS pages could not be fetched.

**[P101]** Subashchandar, N. S. M.; Zhao, L.; Shalchi, A.; Zank, G. P.; le Roux, J. A.; Li, H.; Zhu, X.; Silwal, A.; et al. (2025). **Parallel and Perpendicular Diffusion of Energetic Particles in the Near-Sun Solar Wind Observed by Parker Solar Probe.** *The Astrophysical Journal Letters*, **991**(2), L30. [DOI](https://doi.org/10.3847/2041-8213/ae063f); [arXiv:2509.10648](https://arxiv.org/abs/2509.10648). Role: Taylor mapping and fitted PSD model for SOQLT input (Eq. 4). Reading: full text. Key in the supplied bibliography: **Subashchandar-2025-AJL**.

**[P102]** Shalchi, A.; Bieber, J. W.; Matthaeus, W. H.; Qin, G. (2004). **Nonlinear Parallel and Perpendicular Diffusion of Charged Cosmic Rays in Weak Turbulence.** *The Astrophysical Journal*, **616**(1), 617–629. [DOI](https://doi.org/10.1086/424839). Role: weakly nonlinear theory (secondary via the Shalchi 2009 book chapter abstract). Reading: secondary. Pre-2006 foundation. Proposed key (additions file): **Shalchi-2004-AJ-616**.

**[P103]** Shalchi, A. (2009). **Nonlinear Cosmic Ray Diffusion Theories.** *Astrophysics and Space Science Library*, **362**, XIII+199. [DOI](https://doi.org/10.1007/978-3-642-00309-7). Role: WNLT chapter abstract (broadening increases D_μμ). Reading: abstract only. Key in the supplied bibliography: **Shalchi-2009**. Metadata not fully confirmed: volume (series no. 362) and pages (XIII+199) not confirmed: OpenAlex gives no series volume or page count and the Springer page was rate limited.

**[P104]** Tautz, R. C.; Dosch, A.; Effenberger, F.; Fichtner, H.; Kopp, A. (2013). **Pitch-angle scattering in magnetostatic turbulence. I. Test-particle simulations and the validity of analytical results.** *Astronomy & Astrophysics*, **558**, A147. [DOI](https://doi.org/10.1051/0004-6361/201322142); [arXiv:1309.7836](https://arxiv.org/abs/1309.7836). Role: running estimator (Eq. 3), time-dependent QLT (Eqs. 7-8), Table 1 χ². Reading: full text. Key in the supplied bibliography: **Tautz-2013-AA-A147**.

**[P105]** Qin, G.; Shalchi, A. (2009). **Pitch-Angle Diffusion Coefficients of Charged Particles from Computer Simulations.** *The Astrophysical Journal*, **707**(1), 61–66. [DOI](https://doi.org/10.1088/0004-637X/707/1/61). Role: simulated D_μμ vs QLT/SOQLT/WNLT (abstract; set-up secondary via van den Berg et al. 2024). Reading: abstract only. Key in the supplied bibliography: **Qin-2009-AJ**.

**[P106]** van den Berg, J. P.; Els, P. L.; Engelbrecht, N. E. (2024). **An Evaluation of Different Numerical Methods to Calculate the Pitch-angle Diffusion Coefficient from Full-orbit Simulations: Disentangling a Rope of Sand.** *The Astrophysical Journal*, **977**(2), 174. [DOI](https://doi.org/10.3847/1538-4357/ad8b4a); [arXiv:2410.21001](https://arxiv.org/abs/2410.21001). Role: nine D_μμ estimators M0-M5, benchmark cases 1-3. Reading: full text. Key in the supplied bibliography: **vandenBerg-2024-AJ**.

**[P107]** Dundovic, A.; Pezzi, O.; Blasi, P.; Evoli, C.; Matthaeus, W. H. (2020). **Novel aspects of cosmic ray diffusion in synthetic magnetic turbulence.** *Physical Review D*, **102**(10), 103016. [DOI](https://doi.org/10.1103/PhysRevD.102.103016); [arXiv:2007.09142](https://arxiv.org/abs/2007.09142). Role: histogram D_μμ at lags 10⁻² and 10⁻¹ λ∥/c; QLT/SOQLT validity ranges. Reading: full text. Key in the supplied bibliography: **Dundovic-2020-PRD**.

**[P108]** Ivascenko, A.; Lange, S.; Spanier, F.; Vainio, R. (2016). **Determining pitch-angle diffusion coefficients from test particle simulations.** *The Astrophysical Journal*, **833**(2), 223. [DOI](https://doi.org/10.3847/1538-4357/833/2/223); [arXiv:1610.03240](https://arxiv.org/abs/1610.03240). Role: estimators MI, MIIa, MIIb (Eqs. 6, 9-11); failure at δB/B₀ ≈ 1. Reading: full text. Key in the supplied bibliography: **Ivascenko-2016-AJ**.

**[P109]** Lange, S.; Spanier, F.; Battarbee, M.; Vainio, R.; Laitinen, T. (2013). **Particle scattering in turbulent plasmas with amplified wave modes.** *Astronomy & Astrophysics*, **553**, A129. [DOI](https://doi.org/10.1051/0004-6361/201220804); [arXiv:1303.7463](https://arxiv.org/abs/1303.7463). Role: running estimator (Eq. 24), 10-30 gyroperiod window. Reading: full text. Key in the supplied bibliography: **Lange-2013-AA**.

**[P110]** Weidl, M. S.; Jenko, F.; Teaca, B.; Schlickeiser, R. (2015). **Cosmic-ray pitch-angle scattering in imbalanced MHD turbulence simulations.** *The Astrophysical Journal*, **811**(1), 8. [DOI](https://doi.org/10.1088/0004-637X/811/1/8); [arXiv:1508.03156](https://arxiv.org/abs/1508.03156). Role: running-estimator definition in MHD simulations (reading truncated before the estimator settings). Reading: full text. Key in the supplied bibliography: **Weidl-2015-AJ**.

**[P111]** Tautz, R. C. (2013). **Pitch-angle scattering in magnetostatic turbulence. II. Analytical considerations and pitch-angle isotropization.** *Astronomy & Astrophysics*, **558**, A148. [DOI](https://doi.org/10.1051/0004-6361/201322143); [arXiv:1309.7838](https://arxiv.org/abs/1309.7838). Role: isotropic reference solution (Eqs. 11-22), running-estimator bias, correlation estimator (Eqs. 25-28). Reading: full text. Key in the supplied bibliography: **Tautz-2013-AA-A148**.

**[P112]** Tautz, R. C.; Shalchi, A.; Dosch, A. (2014). **Pitch-angle scattering of energetic particles with adiabatic focusing.** *The Astrophysical Journal*, **794**(2), 138. [DOI](https://doi.org/10.1088/0004-637X/794/2/138); [arXiv:1408.6947](https://arxiv.org/abs/1408.6947). Role: focused orbit, focusing contribution D^(0) (Eq. 13), estimator ⟨μ̇(0)Δμ(t)⟩ (Eq. 34). Reading: full text. Key in the supplied bibliography: **Tautz-2014-AJ**.

**[P113]** Sun, P.; Jokipii, J. R.; Giacalone, J. (2016). **Pitch-angle scattering of energetic charged particles in nearly constant magnitude magnetic turbulence.** *The Astrophysical Journal*, **827**(1), 16. [DOI](https://doi.org/10.3847/0004-637X/827/1/16). Role: origin of the stationary absorbing-boundary method M5 (secondary via van den Berg et al. 2024). Reading: abstract only. Proposed key (additions file): **Sun-2016-AJ**.

**[P114]** Qin, G.; Shalchi, A. (2014). **Detailed numerical investigation of 90° scattering of energetic particles interacting with magnetic turbulence.** *Physics of Plasmas*, **21**(4). [DOI](https://doi.org/10.1063/1.4873895). Role: 9 μ bins (secondary via van den Berg et al. 2024). Reading: secondary. Proposed key (additions file): **Qin-2014-PoP**. Metadata not fully confirmed: pages: article number (042906 per the task brief) not confirmed; OpenAlex has no page/article number and the AIP and ADS pages could not be fetched.

**[P115]** Hussein, M.; Tautz, R. C.; Shalchi, A. (2015). **The influence of different turbulence models on the diffusion coefficients of energetic particles.** *Journal of Geophysical Research: Space Physics*, **120**(6), 4095–4111. [DOI](https://doi.org/10.1002/2015JA021060); [arXiv:1505.05099](https://arxiv.org/abs/1505.05099). Role: test-particle λ∥ tables (via MFP v1.2 Section 7.6). Reading: full text. Key in the supplied bibliography: **Hussein-2015-JGR, Hussein-2015-JGRSP**.

**[P116]** Reichherzer, P.; Becker Tjus, J.; Zweibel, E. G.; Merten, L.; Pueschel, M. J. (2020). **Turbulence-level dependence of cosmic ray parallel diffusion.** *Monthly Notices of the Royal Astronomical Society*, **498**(4), 5051–5064. [DOI](https://doi.org/10.1093/mnras/staa2533); [arXiv:1910.07528](https://arxiv.org/abs/1910.07528). Role: QLT validity b/B ≲ 0.05 (κ∥ simulations). Reading: full text. Key in the supplied bibliography: **Reichherzer-2020-MNRAS**.

**[P117]** Litvinenko, Y. E. (2012). **The parallel diffusion coefficient and adiabatic focusing of cosmic-ray particles.** *The Astrophysical Journal*, **745**(1), 62. [DOI](https://doi.org/10.1088/0004-637X/745/1/62). Role: Itô SDE with focusing (Eqs. 8-13), weak-focusing λ/λ₀ (Eq. 5), statement on Shalchi (2011). Reading: full text. Proposed key (additions file): **Litvinenko-2012-AJ**.

**[P118]** Shalchi, A. (2011). **Improved Analytical Description of Parallel Diffusion with Adiabatic Focusing.** *The Astrophysical Journal*, **728**(2), 113. [DOI](https://doi.org/10.1088/0004-637X/728/2/113). Role: λ/λ₀ = 1 − (2/5)(λ₀/L)² (secondary via Litvinenko 2012a Eq. 6). Reading: secondary. Proposed key (additions file): **Shalchi-2011-AJ**. Note: title capitalization converted from the all-caps IOP metadata to title case.

**[P119]** Litvinenko, Y. E. (2012). **Effects of non-isotropic scattering, magnetic helicity, and adiabatic focusing on diffusive transport of solar energetic particles.** *The Astrophysical Journal*, **752**(1), 16. [DOI](https://doi.org/10.1088/0004-637X/752/1/16). Role: κ/(vL) = coth ξ − 1/ξ (Eq. 20); helicity form (Eqs. 4, 21, 24). Reading: full text. Proposed key (additions file): **Litvinenko-2012-AJ-752**.

**[P120]** Wang, J. F.; Qin, G. (2019). **The Diffusion Coefficient with Displacement Variance of Energetic Particles Caused by Adiabatic Focusing.** *The Astrophysical Journal*, **886**, 89. [DOI](https://doi.org/10.3847/1538-4357/ab505e); [arXiv:1909.11547](https://arxiv.org/abs/1909.11547). Role: Distinction between the transport-equation (closure) diffusion coefficient and the displacement-variance definition under focusing (Section 12.3; Abstract only). Reading: abstract only. Key in the supplied bibliography: **Wang-2019-AJ**. Note: journal, volume and article number from the supplied bibliography record (Wang-2019-AJ) and the arXiv abstract page (DOI); abstract read on arXiv.

**[P121]** Wang, J. F.; Qin, G. (2018). **Parallel and Perpendicular Diffusion Coefficients of Energetic Charged Particles with Adiabatic Focusing.** *The Astrophysical Journal*, **868**(2), 139. [DOI](https://doi.org/10.3847/1538-4357/aae927); [arXiv:1707.08058](https://arxiv.org/abs/1707.08058). Role: κ∥^a (Eq. 39) and full κ∥ ≈ κ∥0(1 + 2ξ²/15) (Eqs. 62-63). Reading: full text. Key in the supplied bibliography: **Wang-2018-AJ**.

**[P122]** He, H.-Q.; Schlickeiser, R. (2014). **Modification of the Parallel Scattering Mean Free Path of Cosmic Rays in the Presence of Adiabatic Focusing.** *The Astrophysical Journal*, **792**(2), 85. [DOI](https://doi.org/10.1088/0004-637X/792/2/85). Role: agreement with κ∥^a (secondary via Wang & Qin 2018). Reading: secondary. Proposed key (additions file): **He-2014-AJ**.

**[P123]** Shalchi, A.; Danos, R. J. (2013). **On the Different Analytical Results Obtained for the Parallel Diffusion Coefficient of Cosmic Particles with Adiabatic Focusing.** *The Astrophysical Journal*, **765**(2), 153. [DOI](https://doi.org/10.1088/0004-637X/765/2/153). Role: correction term T (secondary via Wang & Qin 2018). Reading: secondary. Proposed key (additions file): **Shalchi-2013-AJ**. Note: title capitalization converted from the all-caps IOP metadata to title case.

**[P124]** Tautz, R. C.; Dosch, A.; Lerche, I. (2012). **Simulating cosmic-ray transport with adiabatic focusing.** *Astronomy & Astrophysics*, **545**, A149. [DOI](https://doi.org/10.1051/0004-6361/201219636). Role: quasilinear and nonlinear focusing λ∥ (Eqs. 9-10); simulation finding. Reading: full text. Key in the supplied bibliography: **Tautz-2012-AA**.

**[P125]** He, H.-Q.; Wan, W. (2013). **The dependence of the parallel and perpendicular mean free paths on the rigidity of the solar energetic particles: theoretical model versus observations.** *Astronomy & Astrophysics*, **557**, A57. [DOI](https://doi.org/10.1051/0004-6361/201321025). Role: tanh focusing correction of λ∥ (via MFP v1.2 Section 6.9). Reading: full text. Proposed key (additions file): **He-2013-AA**.

**[P126]** Battarbee, M.; Laitinen, T.; Vainio, R.; Agueda, N. (2010). **Acceleration of Energetic Particles Through Self-Generated Waves in a Decelerating Coronal Shock.** *AIP Conference Proceedings (Solar Wind 12)*, **1216**, 84–87. [DOI](https://doi.org/10.1063/1.3395969); [arXiv:1303.4334](https://arxiv.org/abs/1303.4334). Role: CSA wave equation, frequency growth rate, streaming bookkeeping, Van Leer advection. Reading: full text. Key in the supplied bibliography: **Battarbee-2010-AIP**.

**[P127]** Vainio, R.; Pönni, A.; Battarbee, M.; Koskinen, H. E. J.; Afanasiev, A.; Laitinen, T. (2014). **A semi-analytical foreshock model for energetic storm particle events inside 1 AU.** *Journal of Space Weather and Space Climate*, **4**, A08. [DOI](https://doi.org/10.1051/swsc/2014005). Role: Bell-type foreshock λ(x, E), validity limits. Reading: full text. Key in the supplied bibliography: **Vainio-2014-JSWSC**.

**[P128]** Afanasiev, A.; Vainio, R. (2013). **Monte Carlo Simulation Model of Energetic Proton Transport Through Self-Generated Alfvén Waves.** *The Astrophysical Journal Supplement Series*, **207**(2), 29. [DOI](https://doi.org/10.1088/0067-0049/207/2/29). Role: SOLPACS wave-growth numerics (not accessible). Reading: metadata only. Key in the supplied bibliography: **Afanasiev-2013-ApJS**.

**[P129]** Strauss, R. D.; Fichtner, H. (2014). **Cosmic ray anisotropies near the heliopause.** *Astronomy & Astrophysics*, **572**, L3. [DOI](https://doi.org/10.1051/0004-6361/201424842). Role: 2-D pitch-angle-dependent heliopause model; D_⊥(μ) forms. Reading: full text. Proposed key (additions file): **Strauss-2014-AA**.

**[P130]** Zhang, M.; Zuo, P.; Pogorelov, N. (2014). **Heliospheric Influence on the Anisotropy of TeV Cosmic Rays.** *The Astrophysical Journal*, **790**(1), 5. [DOI](https://doi.org/10.1088/0004-637X/790/1/5). Role: D_μμ-based LISM anisotropy for TeV protons. Reading: full text. Proposed key (additions file): **Zhang-2014-AJ**.

**[P131]** Zhang, M. (1999). **A Markov Stochastic Process Theory of Cosmic-Ray Modulation.** *The Astrophysical Journal*, **513**(1), 409–420. [DOI](https://doi.org/10.1086/306857). Role: time-backward Markov method used by Qin & Wu (2020). Reading: abstract only. Pre-2006 foundation. Key in the supplied bibliography: **Zhang-1999-AJ**.

**[P132]** Florinski, V.; Stone, E. C.; Cummings, A. C.; le Roux, J. A. (2015). **Cosmic rays beyond the boundary of the heliosphere.** *Proceedings of Science (ICRC2015)*, **236**, 214. [DOI](https://doi.org/10.22323/1.236.0214). Role: isotropic scattering rate in the heliosheath and LISM. Reading: full text. Proposed key (additions file): **Florinski-2015-ICRC**. Note: year given as the conference year (ICRC2015); the PoS volume 236 was published 2016-08-18.

**[P133]** Luo, X.; Potgieter, M. S.; Zhang, M.; Pogorelov, N. V.; Feng, X.; Strauss, R. D. (2016). **A Numerical Simulation of Cosmic Ray Modulation near the Heliopause. II. Some Physical Insights.** *The Astrophysical Journal*, **826**(2), 182. [DOI](https://doi.org/10.3847/0004-637X/826/2/182). Role: Parker-equation SDE modulation model (not pitch-angle-resolved). Reading: full text. Proposed key (additions file): **Luo-2016-AJ**.

**[P134]** Shen, Z.-N.; Qin, G. (2018). **Modulation of Galactic Cosmic Rays in the Inner Heliosphere over Solar Cycles.** *The Astrophysical Journal*, **854**(2), 137. [DOI](https://doi.org/10.3847/1538-4357/aaab64); [arXiv:1709.08017](https://arxiv.org/abs/1709.08017). Role: Parker-equation GCR modulation (listed as not pitch-angle-resolved); read as arXiv:1709.08017. Reading: full text. Key in the supplied bibliography: **Shen-2018-AJ**.

**[P135]** Kóta, J.; Jokipii, J. R. (2014). **Are cosmic rays modulated beyond the heliopause?.** *The Astrophysical Journal*, **782**(1), 24. [DOI](https://doi.org/10.1088/0004-637X/782/1/24). Role: Parker-equation model beyond the heliopause (not pitch-angle-resolved). Reading: abstract only. Proposed key (additions file): **Kota-2014-AJ**.

**[P136]** Kubo, Y.; Shimazu, H. (2010). **Effect of Finite Larmor Radius on the Cosmic Ray Penetration into an Interplanetary Magnetic Flux Rope.** *arXiv preprint*. [arXiv:1007.2714](https://arxiv.org/abs/1007.2714). Role: full-orbit flux-rope model without scattering. Reading: full text. Proposed key (additions file): **Kubo-2010-arXiv**. Note: published as Kubo, Y.; Shimazu, H. 2010, 'Effect of Finite Larmor Radius on Cosmic-Ray Penetration into an Interplanetary Magnetic Flux Rope', The Astrophysical Journal 720(1), 853-861, DOI 10.1088/0004-637X/720/1/853.

**[P137]** Benella, S.; Laurenza, M.; Vainio, R.; Grimani, C.; Consolini, G.; Hu, Q.; Afanasiev, A. (2020). **A New Method to Model Magnetic Cloud-driven Forbush Decreases: The 2016 August 2 Event.** *The Astrophysical Journal*, **901**(1), 21. [DOI](https://doi.org/10.3847/1538-4357/abac59). Role: full-orbit Forbush-decrease model without scattering. Reading: full text. Proposed key (additions file): **Benella-2020-AJ**.

**[P138]** Lara, A.; Borgazzi, A.; Guennam, E.; Niembro, T.; Arunbabu, K. P. (2024). **Interaction of Cosmic Rays With Magnetic Flux Ropes.** *Journal of Geophysical Research: Space Physics*, **129**(8), e2024JA032478. [DOI](https://doi.org/10.1029/2024JA032478). Role: full-orbit model with white noise added to B (no scattering operator). Reading: full text. Proposed key (additions file): **Lara-2024-JGR**. Metadata not fully confirmed: article number e2024JA032478 not read in any source (it follows the AGU eLocator pattern of the DOI, but neither OpenAlex nor a publisher page shows it); second author's middle initial not confirmed (OpenAlex display name 'Andrea Inés Borgazzi', raw byline 'A. Borgazzi').

**[P139]** Florinski, V.; Alouani-Bibi, F.; Kóta, J.; Guo, X. (2012). **Cosmic-ray diffusion in a sectored magnetic field in the distant heliosheath.** *The Astrophysical Journal*, **754**(1), 31. [DOI](https://doi.org/10.1088/0004-637X/754/1/31); [arXiv:1302.2927](https://arxiv.org/abs/1302.2927). Role: gyrophase scattering at fixed pitch angle π/2 (not pitch-angle-resolved). Reading: full text. Key in the supplied bibliography: **Florinski-2012-AJ**.

**[P140]** Kóta, J.; Jokipii, J. R. (2017). **Transient Cosmic-ray Events beyond the Heliopause: Interpreting Voyager-1 Observations.** *The Astrophysical Journal*, **839**(2), 126. [DOI](https://doi.org/10.3847/1538-4357/aa6c64). Role: pitch-angle-dependent energisation beyond the heliopause (abstract only). Reading: abstract only. Proposed key (additions file): **Kota-2017-AJ**.

**[P141]** Florinski, V.; Stone, E. C.; Cummings, A. C.; le Roux, J. A. (2015). **Energetic Particle Anisotropies at the Heliospheric Boundary. II. Transient Features and Rigidity Dependence.** *The Astrophysical Journal*, **803**(1), 47. [DOI](https://doi.org/10.1088/0004-637X/803/1/47). Role: pitch-angle-dependent heliopause model (abstract only). Reading: abstract only. Key in the supplied bibliography: **Florinski-2015-AJ**.

**[P142]** Florinski, V.; Guzman, J. A.; Kleimann, J.; Baliukin, I.; Ghanbari, K.; Turner, D.; Zieger, B.; Kóta, J.; et al. (2024). **Magnetic Trapping of Galactic Cosmic Rays in the Outer Heliosheath and Their Preferential Entry into the Heliosphere.** *The Astrophysical Journal*, **961**(2), 244. [DOI](https://doi.org/10.3847/1538-4357/ad0b15). Role: GCR transport with pitch-angle scattering in the outer heliosheath (abstract only). Reading: abstract only. Key in the supplied bibliography: **Florinski-2024-AJ**. Metadata not fully confirmed: second author's surname not confirmed: the record has 'Guzman, J. A.', but his name variants indicate the compound surname 'Alonso Guzman' (i.e. 'Alonso Guzman, J.'); how the journal splits the name was not read (Crossref not consulted). All other fields confirmed.

**[P143]** Dröge, W.; Kartavykh, Y.; Klassen, A.; Dresing, N.; Lario, D. (2015). **Multi-spacecraft observations and transport modeling of energetic electrons for a series of solar particle events in August 2010.** *Proceedings of the 34th International Cosmic Ray Conference (ICRC2015), PoS(ICRC2015)*, **236**, 206. [DOI](https://doi.org/10.22323/1.236.0206). Role: Forward Itô μ-SDE (Eq. 3.2). Reading: full text. Key in the supplied bibliography: **Droge-2016-ICRC**. Note: year given as the conference year (ICRC2015); the PoS volume 236 was published 2016-08-18.

**[P144]** van den Berg, J. P. (2019). **Modelling of solar energetic particles by stochastic processes.** *MSc dissertation, North-West University (Potchefstroom), hdl:10394/33908*. Role: Statement that higher-order schemes give no improvement. Reading: abstract only. Proposed key (additions file): **vandenBerg-2019-misc**. Metadata not fully confirmed: year: 2019 not confirmed (the repository's Date field is blank and the full metadata view could not be fetched).

**[P145]** Wijsen, N. (2020). **PARADISE: a model for energetic particle transport in the solar wind.** *PhD thesis, Universitat de Barcelona and KU Leuven, hdl:2445/178025*. Role: Euler–Maruyama stability statement (Chapter 1 only read). Reading: full text. Key in the supplied bibliography: **Wijsen-2020-PHD**.

**[P146]** Strauss, R. D. (n.d.). **SEP_propagator (GitHub repository github.com/RDStrauss/SEP_propagator, master branch).** *Computer code (GitHub)*. Role: FD reference code of the primer: CFL step and zero-flux μ faces. Reading: full text. Proposed key (additions file): **Strauss-nd-Code**. Metadata not fully confirmed: year: not determinable (no release, tag or dated citation; commit history not reachable); version used not identified.

**[P147]** Marsh, M. S.; Dalla, S.; Kelly, J.; Laitinen, T. (2013). **Drift-induced Perpendicular Transport of Solar Energetic Particles.** *The Astrophysical Journal*, **774**(1), 4. [DOI](https://doi.org/10.1088/0004-637X/774/1/4); [arXiv:1307.1585](https://arxiv.org/abs/1307.1585). Role: Poisson-timed isotropic redraw, mean time λ/v0. Reading: full text. Key in the supplied bibliography: **Marsh-2013-AJ**.

**[P148]** Kelly, J.; Dalla, S.; Laitinen, T. (2012). **Cross-Field Transport of Solar Energetic Particles in a Large-Scale Fluctuating Magnetic Field.** *The Astrophysical Journal*, **750**(1), 47. [DOI](https://doi.org/10.1088/0004-637X/750/1/47); [arXiv:1202.6010](https://arxiv.org/abs/1202.6010). Role: Poisson-timed isotropic redraw in solar-wind frame. Reading: full text. Key in the supplied bibliography: **Kelly-2012-arXiv**.

## 23. Bibliography comparison

**Snapshot appendix.** This audit describes the supplied file at the date of this document. It is separable from the model specification.

The comparison uses the supplied bibliography (bibliography.bib), containing 7186 parsed records. Its SHA-256 digest is:

~~~text
1d26e45fbdb1eb4c4efd7ef4fbd7e34b80a79f35c754fa48bb68f1c5827f2fc9
~~~

A reference counts as present when its DOI, its arXiv identifier or its normalized title matches a record; title matches with a conflicting year or DOI were checked individually (Section 23.3). "Not found" means no matching record in this snapshot. The additions file contains BibTeX entries only for references not found.

### 23.1 Cited references not found in the supplied bibliography

**34 references from 2006–2026 were not found.** Their entries are in PITCH_ANGLE_SCATTERING_MODEL_additions.bib:

| Reference | Paper | Proposed key | Role |
|---|---|---|---|
| P3 | Agueda et al. (2014) — Release timescales of solar energetic particles in the low corona | Agueda-2014-AA | 1-D FTE with source (Eq. 2), L = -B/(dB/dz) |
| P7 | Kartavykh et al. (2015) — Simulation of energetic particle interactions with a shock wave in a focused transport model | Kartavykh-2015-ICRC | κ₀(r,E)(\|μ\|^{q-1}+H)(1−μ²) with q = 1.67, H = 0.5 (Eq. 4.2) |
| P17 | Young & Poduval (2025) — Quantifying the Effects of Parameters in Widespread SEP Events with EPREM | Young-2025-arXiv | EPREM Lagrangian FTE (Eq. A1) |
| P45 | Riordan & Pe'er (2019) — Pitch-angle Diffusion and Bohm-type Approximations in Diffusive Shock Acceleration | Riordan-2019-AJ | MSD estimator without ½ (Eq. 11), Δt = 20 t_g window, plateau artefact |
| P46 | Florinski (2013) — Entering the "magnetic highway": energetic particle anisotropies at the heliospheric boundary | Florinski-2013-arXiv | gyrophase-resolved heliopause model; random walk on a sphere; λ∥ = v/(2D) |
| P48 | Agueda et al. (2010) — Solar near-relativistic electron observations as a proof of a back-scatter region beyond 1 AU during the 2000 February 18 event | Agueda-2010-AA | ε-form Eq. 1, λ∥ = (v/ν₀)φ(ε) and asymptote Eqs. 3–4; ε = 0.01; ν(1) = 3.5v/λ |
| P49 | Agueda et al. (2009) — Modeling solar near-relativistic electron events: Insights into solar injection and interplanetary transport conditions | Agueda-2009-AA | ε-form in the solar-wind frame; ε = 0.10, 0.01 and isotropic fits |
| P52 | Agueda & Vainio (2013) — On the parametrization of the energetic-particle pitch-angle diffusion coefficient | Agueda-2013-JSWSC | q-form with ν₀ normalization (Eqs. 2–4), ε-form (Eq. 5), ε–(q,H) least-squares equivalences |
| P55 | Àgueda Costafreda (2008) — Near-relativistic electron events. Monte Carlo simulations of solar injection and interplanetary transport | AguedaCostafreda-2008-PHD | ε-form origin (Eqs. 3.19, 3.21), H = 0 normalization (Eq. 3.8), hemispheric η transformation (Eq. 3.20) |
| P57 | Dröge & Kartavykh (2008) — Testing transport theories with solar energetic particles | Droge-2008-ICRC | D^r_μμ = κ₀(r,R){\|μ\|^{q-1}+H}(1−μ²) with q = 1.67, H = 0.05 (Eq. 2) |
| P61 | Agueda et al. (2008) — Injection and Interplanetary Transport of Near-Relativistic Electrons: Modeling the Impulsive Event on 2000 May 1 | Agueda-2008-AJ | origin of the ε-form (journal paper not read) |
| P79 | Shalchi et al. (2006) — Parallel and Perpendicular Transport of Heliospheric Cosmic Rays in an Improved Dynamical Turbulence Model | Shalchi-2006-AJ | NADT slab and 2D D_μμ (structure only, OCR) |
| P83 | Sandroos et al. (2014) — First- and Second-order Fermi Acceleration at Parallel Shocks | Sandroos-2014-arXiv | D_μμ in parallel/antiparallel wave rest frames (Eqs. 29–31); L/R helicity labels |
| P90 | Schreiner & Spanier (2015) — Effects of dispersive wave modes on charged particles transport | Schreiner-2015-ICRC | wave-frame μ for dispersive waves (Eqs. 2.4, 2.9–2.11) |
| P96 | Xu & Lazarian (2018) — Resonance-broadened Transit Time Damping of Particles in MHD Turbulence | Xu-2018-AJ | broadened TTD for slow and fast modes; fast gyroresonance (Eqs. 1–75) |
| P113 | Sun et al. (2016) — Pitch-angle scattering of energetic charged particles in nearly constant magnitude magnetic turbulence | Sun-2016-AJ | origin of the stationary absorbing-boundary method M5 (secondary via van den Berg et al. 2024) |
| P114 | Qin & Shalchi (2014) — Detailed numerical investigation of 90° scattering of energetic particles interacting with magnetic turbulence | Qin-2014-PoP | 9 μ bins (secondary via van den Berg et al. 2024) |
| P117 | Litvinenko (2012) — The parallel diffusion coefficient and adiabatic focusing of cosmic-ray particles | Litvinenko-2012-AJ | Itô SDE with focusing (Eqs. 8-13), weak-focusing λ/λ₀ (Eq. 5), statement on Shalchi (2011) |
| P118 | Shalchi (2011) — Improved Analytical Description of Parallel Diffusion with Adiabatic Focusing | Shalchi-2011-AJ | λ/λ₀ = 1 − (2/5)(λ₀/L)² (secondary via Litvinenko 2012a Eq. 6) |
| P119 | Litvinenko (2012) — Effects of non-isotropic scattering, magnetic helicity, and adiabatic focusing on diffusive transport of solar energetic particles | Litvinenko-2012-AJ-752 | κ/(vL) = coth ξ − 1/ξ (Eq. 20); helicity form (Eqs. 4, 21, 24) |
| P122 | He & Schlickeiser (2014) — Modification of the Parallel Scattering Mean Free Path of Cosmic Rays in the Presence of Adiabatic Focusing | He-2014-AJ | agreement with κ∥^a (secondary via Wang & Qin 2018) |
| P123 | Shalchi & Danos (2013) — On the Different Analytical Results Obtained for the Parallel Diffusion Coefficient of Cosmic Particles with Adiabatic Focusing | Shalchi-2013-AJ | correction term T (secondary via Wang & Qin 2018) |
| P125 | He & Wan (2013) — The dependence of the parallel and perpendicular mean free paths on the rigidity of the solar energetic particles: theoretical model versus observations | He-2013-AA | tanh focusing correction of λ∥ (via MFP v1.2 Section 6.9) |
| P129 | Strauss & Fichtner (2014) — Cosmic ray anisotropies near the heliopause | Strauss-2014-AA | 2-D pitch-angle-dependent heliopause model; D_⊥(μ) forms |
| P130 | Zhang et al. (2014) — Heliospheric Influence on the Anisotropy of TeV Cosmic Rays | Zhang-2014-AJ | D_μμ-based LISM anisotropy for TeV protons |
| P132 | Florinski et al. (2015) — Cosmic rays beyond the boundary of the heliosphere | Florinski-2015-ICRC | isotropic scattering rate in the heliosheath and LISM |
| P133 | Luo et al. (2016) — A Numerical Simulation of Cosmic Ray Modulation near the Heliopause. II. Some Physical Insights | Luo-2016-AJ | Parker-equation SDE modulation model (not pitch-angle-resolved) |
| P135 | Kóta & Jokipii (2014) — Are cosmic rays modulated beyond the heliopause? | Kota-2014-AJ | Parker-equation model beyond the heliopause (not pitch-angle-resolved) |
| P136 | Kubo & Shimazu (2010) — Effect of Finite Larmor Radius on the Cosmic Ray Penetration into an Interplanetary Magnetic Flux Rope | Kubo-2010-arXiv | full-orbit flux-rope model without scattering |
| P137 | Benella et al. (2020) — A New Method to Model Magnetic Cloud-driven Forbush Decreases: The 2016 August 2 Event | Benella-2020-AJ | full-orbit Forbush-decrease model without scattering |
| P138 | Lara et al. (2024) — Interaction of Cosmic Rays With Magnetic Flux Ropes | Lara-2024-JGR | full-orbit model with white noise added to B (no scattering operator) |
| P140 | Kóta & Jokipii (2017) — Transient Cosmic-ray Events beyond the Heliopause: Interpreting Voyager-1 Observations | Kota-2017-AJ | pitch-angle-dependent energisation beyond the heliopause (abstract only) |
| P144 | van den Berg (2019) — Modelling of solar energetic particles by stochastic processes | vandenBerg-2019-misc | Statement that higher-order schemes give no improvement |
| P146 | Strauss (n.d.) — SEP_propagator (GitHub repository github.com/RDStrauss/SEP_propagator, master branch) | Strauss-nd-Code | FD reference code of the primer: CFL step and zero-flux μ faces |

The following **10 pre-2006 foundations** are also not in the supplied bibliography and are included in the additions file:

| Reference | Paper | Proposed key | Role |
|---|---|---|---|
| P35 | Kóta et al. (2005) — SEP Acceleration at Realistic CMEs: Two Sites of Acceleration? | Kota-2005-ICRC | foundation of the EPREM co-moving (Lagrangian) FTE; equation garbled |
| P59 | Teufel & Schlickeiser (2003) — Analytic calculation of the parallel mean free path of heliospheric cosmic rays. II. Dynamical magnetic slab turbulence and random sweeping slab turbulence with finite wave power at small wavenumbers | Teufel-2003-AA | slab λ∥ form whose R_L k_min ≪ 1 limit equals Eq. (28) at h = 0 (via MFP v1.2 Sect. 6.4) |
| P70 | Bieber et al. (1994) — Proton and Electron Mean Free Paths: The Palmer Consensus Revisited | Bieber-1994-AJ-420 | cited by Lang et al. (2024) for the q-form (not read) |
| P71 | Dröge (2000) — The Rigidity Dependence of Solar Particle Scattering Mean Free Paths | Droge-2000-AJ | cited by Lang et al. (2024) and Agueda et al. (2010) for the q-form / finite scattering through μ = 0 (not read) |
| P72 | Teufel & Schlickeiser (2002) — Analytic calculation of the parallel mean free path of heliospheric cosmic rays. I. Dynamical magnetic slab turbulence and random sweeping slab turbulence | Teufel-2002-AA | DT/RS slab D_μμ (Eqs. 34-36, 41-54), 8π normalization, σ tensor |
| P88 | Vainio & Kocharov (2001) — Proton transport through self-generated waves in impulsive flares | Vainio-2001-AA | isotropic scattering in each ± wave frame; −ΔE = ∓V_A p_± Δμ_± |
| P91 | Virtanen & Vainio (2005) — Stochastic acceleration in relativistic parallel shocks | Virtanen-2005-AJ | sequential forward/backward wave-frame scattering; momentum diffusion A_2 (Eq. 15) |
| P92 | Vainio et al. (2003) — A simple analytical expression for the power spectrum of cascading Alfvén waves in the solar wind | Vainio-2003-AA | D_μμ from a fixed-frame frequency spectrum, f_r = (Ω/2π)V/(v\|μ\|) |
| P100 | Shalchi (2005) — Second-order quasilinear theory of cosmic ray transport | Shalchi-2005-PoP | origin of the SOQLT resonance function (secondary via Shalchi et al. 2009) |
| P102 | Shalchi et al. (2004) — Nonlinear Parallel and Perpendicular Diffusion of Charged Cosmic Rays in Weak Turbulence | Shalchi-2004-AJ-616 | weakly nonlinear theory (secondary via the Shalchi 2009 book chapter abstract) |

### 23.2 Cited references found in the supplied bibliography

| Reference | Matching key(s) | Notes |
|---|---|---|
| P1 | vanDenBerg-2020-arXiv | matched by title and authors; the library record has journal 'arXiv' and DOI 10.1007/s11214-020-00755-2, whereas Crossref gives 10.1007/s11214-020-00771-x for this paper (Section 23.3) |
| P2 | Effenberger-2014-AJ | — |
| P4 | Strauss-2022-AJ | — |
| P5 | Qin-2005-AJ | matched by title, volume and pages; the library DOI 10.1086/430201 belongs to the next article in the same volume (OpenAlex), the correct DOI is 10.1086/430136 (Section 23.3) |
| P6 | Ruffolo-1995-AJ | Ruffolo-1995-AJ: library entry has no DOI |
| P8 | Laitinen-2016-AA | Laitinen-2016-AA: library entry has no DOI |
| P9 | Merten-2025-CPC | matched by arXiv identifier; the library record is the published version (Computer Physics Communications 311, 109542, 2025); this review read the arXiv version 2410.01472 |
| P10 | Tenishev-2022-arXiv | — |
| P11 | Strauss-2015-AJ | Strauss-2015-AJ: library entry has no DOI |
| P12 | Strauss-2017-AJ-837 | — |
| P13 | Strauss-2019-AJ | — |
| P14 | vandenBerg-2026-AA | — |
| P15 | Zhang-2023-AJSS | — |
| P16 | Kubo-2015-EPS | matched by title; a second record, Kubo-2015-PS, has the journal name in its year field; neither has a DOI (Section 23.3) |
| P18 | Qin-2021-AJ | the library record is the published version (ApJ 908, 236, 2021, DOI 10.3847/1538-4357/abd77c, confirmed in OpenAlex); this review read the arXiv preprint 2010.07662 |
| P19 | He-2011-AJ | He-2011-AJ: library entry has no DOI |
| P20 | He-2015-AJ-218 | He-2015-AJ-218: library entry has no DOI |
| P21 | He-2019-AJL, He-2019-arXiv | two records of the same paper; neither has a DOI, and He-2019-AJL has no article number |
| P22 | He-2017-AJ | — |
| P23 | Wang-2014-AJ | Wang-2014-AJ: library entry has no DOI |
| P24 | Lang-2024-arXiv | Lang-2024-arXiv: library entry has no DOI |
| P25 | Wijsen-2019-AA-A28 | Wijsen-2019-AA-A28: library entry has no DOI |
| P26 | Husidic-2024-JSWSC | — |
| P27 | Wang-2015-AJ | — |
| P28 | Strauss-2020-JP | Strauss-2020-JP: library entry has no DOI |
| P29 | Afanasiev-2025-JSWSC | — |
| P30 | Isenberg-1997-JGR, isenberg97 | two records of the same paper, both without DOI |
| P31 | Zuo-2013-AJ | — |
| P32 | Linker-2019-JPCS | — |
| P33 | Kozarev2013 | — |
| P34 | Schwadron-2010-SW | Schwadron-2010-SW: library entry has no DOI |
| P36 | Borovikov-2019-arXiv | — |
| P37 | Afanasiev-2023-AA | Afanasiev-2023-AA: library entry has no DOI |
| P38 | Ng-2008-AJL | matched by title, authors, volume and page; the library DOI 10.1086/592998 belongs to another article (Crossref), the correct DOI is 10.1086/592996; a second record, Ng-2008-NASA, has the same title (Section 23.3) |
| P39 | Zhang-2017-AJ | Zhang-2017-AJ: library entry has no DOI |
| P40 | Zhang-2009-AJ | Zhang-2009-AJ: library entry has no DOI |
| P41 | Qin-2006-JGR | Qin-2006-JGR: library entry has no DOI |
| P42 | Ding-2022-AA, Ding-2022-arXiv | Ding-2022-arXiv: library entry has no DOI |
| P43 | Droge-2010-AJ | Droge-2010-AJ: library entry has no DOI |
| P44 | Strauss-2017-SSR | Strauss-2017-SSR: library entry has no DOI |
| P47 | Chen-2024-AJ-965 | Chen-2024-AJ-965: library entry has no DOI |
| P50 | Beeck-1986-AJ | Beeck-1986-AJ: library entry has no DOI |
| P51 | Strauss-2017-SP | — |
| P53 | Qin-2015-AJ | Qin-2015-AJ: library entry has no DOI |
| P54 | Lavasa-2026-AA | — |
| P56 | He-2015-AJ-814 | — |
| P58 | Wu-2020-AJ | — |
| P60 | Qin-2013-AJ | — |
| P62 | Wijsen-2022-AA | — |
| P63 | Pacheco-2019-AA | — |
| P64 | Minoshima-2026-EPS | — |
| P65 | Droge-2003-AJ | — |
| P66 | Liu-2025-AJ | — |
| P67 | Roux-2009-AJ | Roux-2009-AJ: library entry has no DOI |
| P68 | Hu-2017-JGR | Hu-2017-JGR: library entry has no DOI |
| P69 | Young-2021-AJ | Young-2021-AJ: library entry has no DOI |
| P73 | Shalchi-2009-AA | — |
| P74 | Engelbrecht-2022-SSRv | — |
| P75 | Jokipii-1966-AJ | Jokipii-1966-AJ: library entry has no DOI |
| P76 | Shalchi-2008-MNRAS | matched by DOI; the library gives year 2008, the publisher December 2007 (MNRAS 383, 803) |
| P77 | Tautz-2013-JGR, Tautz-2013-JGRSP | — |
| P78 | Afanasiev-2015-AA | matched by DOI; three further records carry the same title: Afanasiev-2016-arXiv (the arXiv version, year 2016), Afanasiev-2018-AA-618 (volume 618, A114 with this title) and Afanasiev-misc (no year) (Section 23.3) |
| P80 | Spanier-2022-Physics | — |
| P81 | Schlickeiser-1989-AJ | — |
| P82 | Ng-2012-AIP | — |
| P84 | Chandran-2000-PRL | — |
| P85 | Yan-2002-PRL | — |
| P86 | Vainio-2000-AJSS | — |
| P87 | Vainio-2003-AA-406 | — |
| P89 | Battarbee-2013-AA-A34 | same article (A&A 535, A34); the library record has year 2013, no volume and the arXiv DOI 10.48550/arxiv.1303.4340 (Section 23.3) |
| P93 | Ng-2003-AJ | Ng-2003-AJ: library entry has no DOI |
| P94 | Yan-2004-AJ | — |
| P95 | Xu-2020-AJ | — |
| P97 | Yan-2008-AJ | — |
| P98 | Xu-2013-AJ | — |
| P99 | Lynn-2012-AJ | — |
| P101 | Subashchandar-2025-AJL | — |
| P103 | Shalchi-2009 | Shalchi-2009: library entry has no DOI |
| P104 | Tautz-2013-AA-A147 | — |
| P105 | Qin-2009-AJ | Qin-2009-AJ: library entry has no DOI |
| P106 | vandenBerg-2024-AJ | — |
| P107 | Dundovic-2020-PRD | — |
| P108 | Ivascenko-2016-AJ | — |
| P109 | Lange-2013-AA | — |
| P110 | Weidl-2015-AJ | — |
| P111 | Tautz-2013-AA-A148 | — |
| P112 | Tautz-2014-AJ | — |
| P115 | Hussein-2015-JGR, Hussein-2015-JGRSP | — |
| P116 | Reichherzer-2020-MNRAS | — |
| P120 | Wang-2019-AJ | — |
| P121 | Wang-2018-AJ | — |
| P124 | Tautz-2012-AA | — |
| P126 | Battarbee-2010-AIP | — |
| P127 | Vainio-2014-JSWSC | — |
| P128 | Afanasiev-2013-ApJS | — |
| P131 | Zhang-1999-AJ | Zhang-1999-AJ: library entry has no DOI |
| P134 | Shen-2018-AJ | Shen-2018-AJ: library entry has no DOI |
| P139 | Florinski-2012-AJ | same paper: the arXiv record 1302.2927 (posted 12 February 2013) gives DOI 10.1088/0004-637X/754/1/31, the published version in ApJ 754, 31 (2012), which is the library record |
| P141 | Florinski-2015-AJ | Florinski-2015-AJ: library entry has no DOI |
| P142 | Florinski-2024-AJ | — |
| P143 | Droge-2016-ICRC | matched by DOI; the library gives the PoS publication year 2016, this list the conference year 2015 |
| P145 | Wijsen-2020-PHD | — |
| P147 | Marsh-2013-AJ | — |
| P148 | Kelly-2012-arXiv | Kelly-2012-arXiv: library entry has no DOI |

### 23.3 Records in the supplied bibliography that need attention

These observations concern records of the supplied bibliography that this review matched or nearly matched; they were checked against Crossref, OpenAlex, arXiv or publisher pages during the bibliographic check (Section 25.3). The supplied file is not modified.

**Wrong or conflicting DOIs.**

| Record | Problem | Value found |
|---|---|---|
| `Ng-2008-AJL` | DOI 10.1086/592998 belongs to another article (an American Naturalist paper, Crossref) | Ng & Reames 2008, ApJL 686, L123: DOI 10.1086/592996 (Crossref) |
| `Ng-2008-NASA` | same paper as `Ng-2008-AJL`; its DOI field 10.48550/arxiv.20080038680 is not an arXiv DOI (the number is a NASA report identifier) | — |
| `Qin-2005-AJ` | DOI 10.1086/430201 belongs to the next article of ApJ 627 (OpenAlex) | Qin et al. 2005, ApJ 627, 562: DOI 10.1086/430136 (OpenAlex) |
| `vanDenBerg-2020-arXiv` | journal "arXiv" and DOI 10.1007/s11214-020-00755-2 | the Primer is Space Sci. Rev. 216, 146 (2020), DOI 10.1007/s11214-020-00771-x (Crossref); what 10.1007/s11214-020-00755-2 refers to was not checked (lookup refused) |
| `Battarbee-2013-AA-A34` | year 2013, no volume, DOI 10.48550/arxiv.1303.4340 | Battarbee, Laitinen & Vainio 2011, A&A 535, A34, DOI 10.1051/0004-6361/201117507 |
| `Battarbee-2013-AJ` | carries the title of Battarbee et al. 2010 (AIP Conf. Proc., DOI 10.1063/1.3395969) with year 2013 and DOI 10.48550/arxiv.1303.4334 | identity not checked |
| `Afanasiev-2018-AA-618` | title of Afanasiev, Battarbee & Vainio 2015 (A&A 584, A81) with volume 618, A114 | title and volume belong to different papers; not resolved here |
| `Afanasiev-2024-JSWSC` | same title as Afanasiev, Wijsen & Vainio 2025 (JSWSC 15, 3), year 2024, DOI "n/a" | JSWSC 15, 3 (2025), DOI 10.1051/swsc/2024039 |
| `Shalchi-2008-MNRAS` | year 2008 | MNRAS 383, 803 is dated December 2007 by the publisher and OpenAlex |

**Malformed fields and duplicates.**

- `Kubo-2015-PS` has the journal name ("Planets and Space") in its year field and duplicates `Kubo-2015-EPS`; neither has the DOI 10.1186/s40623-015-0260-9.
- Duplicate records of one paper: `Hussein-2015-JGR` and `Hussein-2015-JGRSP`; `Tautz-2013-JGR` and `Tautz-2013-JGRSP`; `Ding-2022-AA` and `Ding-2022-arXiv`; `He-2019-AJL` and `He-2019-arXiv`; `Strauss-2019-AJ` and `Strauss-2019-arXiv`; `Isenberg-1997-JGR` and `isenberg97`; `Afanasiev-2015-AA`, `Afanasiev-2016-arXiv` and `Afanasiev-misc`.
- `Merten-2025-CPC` and `Qin-2021-AJ` are the published versions of preprints read here (Section 23.2); `Florinski-2013-AJL` is the later four-author Letter, not the single-author preprint whose equations this review used.

**Records without a DOI for which the DOI was confirmed here.**

| Record | DOI |
|---|---|
| `Afanasiev-2023-AA` | 10.1051/0004-6361/202346220 |
| `Beeck-1986-AJ` | 10.1086/164784 |
| `Chen-2024-AJ-965` | 10.3847/1538-4357/ad33c3 |
| `Droge-2010-AJ` | 10.1088/0004-637X/709/2/912 |
| `Florinski-2015-AJ` | 10.1088/0004-637X/803/1/47 |
| `He-2011-AJ` | 10.1088/0004-637X/734/2/74 |
| `He-2015-AJ-218` | 10.1088/0067-0049/218/2/17 |
| `He-2019-AJL` | 10.3847/2041-8213/ab50bd |
| `Isenberg-1997-JGR` | 10.1029/96JA03671 |
| `Jokipii-1966-AJ` | 10.1086/148912 |
| `Kelly-2012-arXiv` | 10.1088/0004-637X/750/1/47 |
| `Kubo-2015-EPS` | 10.1186/s40623-015-0260-9 |
| `Laitinen-2016-AA` | 10.1051/0004-6361/201527801 |
| `Lang-2024-arXiv` | 10.3847/1538-4357/ad55c3 |
| `Ng-2003-AJ` | 10.1086/375293 |
| `Qin-2009-AJ` | 10.1088/0004-637X/707/1/61 |
| `Qin-2015-AJ` | 10.1088/0004-637X/809/2/177 |
| `Roux-2009-AJ` | 10.1088/0004-637X/693/1/534 |
| `Ruffolo-1995-AJ` | 10.1086/175489 |
| `Schwadron-2010-SW` | 10.1029/2009SW000523 |
| `Shen-2018-AJ` | 10.3847/1538-4357/aaab64 |
| `Strauss-2015-AJ` | 10.1088/0004-637X/801/1/29 |
| `Strauss-2017-SSR` | 10.1007/s11214-017-0351-y |
| `Strauss-2020-JP` | 10.1088/1742-6596/1620/1/012021 |
| `Wang-2014-AJ` | 10.1088/0004-637X/789/2/157 |
| `Wijsen-2019-AA-A28` | 10.1051/0004-6361/201833958 |
| `Young-2021-AJ` | 10.3847/1538-4357/abdf5f |
| `Zhang-1999-AJ` | 10.1086/306857 |
| `Zhang-2009-AJ` | 10.1088/0004-637X/692/1/109 |
| `Zhang-2017-AJ` | 10.3847/1538-4357/aa86a8 |
| `isenberg97` | 10.1029/96JA03671 |

## 24. Companion data

**Version note.** The table below describes the author-reported canonical v1.1 bundle. The currently available bundle is v1.0; the v1.2 implementation supplement is a separate evidence/package layer. Pin and verify each archive independently before using its registry or equation labels (Section 25.7).

The bundle PITCH_ANGLE_SCATTERING_MODEL_DATA.zip has this layout:

| Path | Content |
|---|---|
| `README.md` | layout, field definitions, provenance and verification labels, parsing rules, regeneration instructions |
| `parameters/model_registry.json` | 172 model IDs with section, equation labels, source keys, specification status, runtime state (Section 4.1), what the ID provides (`provides`), physical admissibility (`admissibility`), the stated domain where the source restricts it (`domain`), register items and the derivation records that list them |
| `parameters/dmumu_shape_presets.json` | 55 objects: the empirical shapes of Section 5 and the GCR presets of Section 15 with their printed form, operator convention, normalization, parameter values as printed, source location, verification and flags |
| `parameters/qlt_wave_parameter_sets.json` | 42 parameter sets printed with the quasilinear, dynamical, wave, compressible, nonlinear and measured-spectrum models of Sections 6–10 and 13 |
| `parameters/simulation_benchmarks.json` | 17 simulation set-ups and published simulated values of Sections 9.5, 11 and 12 |
| `parameters/selfgenerated_wave_and_gcr_sets.json` | 21 parameter sets of Sections 14 and 15 |
| `parameters/numerical_methods.json` | 34 printed numerical settings of Sections 2, 14 and 16: time-step rules and constants, particle counts, grids, reflection rules and schemes |
| `benchmark_points.json` | the 89 fixtures of Section 18 (652 stored values) with inputs, outputs, kind, reference and tolerance class |
| `fixtures.py` | the generator of `benchmark_points.json` (mpmath, 40 digits); it also runs an independent cross-check for every fixture and prints the results. Run without arguments it writes `benchmark_points.regenerated.json` next to itself and compares it with the shipped `benchmark_points.json`, which stays unchanged; an explicit output path is required to write anywhere else (version 1.1; the version 1.0 generator overwrote the shipped file by default) |
| `source_key_map.json` | source key → [Pn], DOI and arXiv identifier, and whether the reference is in the supplied bibliography; equation labels → numbers |
| `register_map.json` | the local item IDs used inside the derivation records → the IDs of Section 19, and the items that were dropped, with the reason |
| `derivation_records/DR-02.md` … `DR-16.md` | the working records of the review, one per group of sections: printed discrepancies, decisions, not-found lists, fixture computations, derived results with the scripts that check them and their output, and the data rows from which the parameter files were built (DR-02: Section 2; DR-05: Section 5; DR-06: Sections 6, 7, 13; DR-08: Sections 8, 9; DR-10: Sections 10–12; DR-14: Sections 14, 15; DR-16: Section 16) |
| `reference_verification.py` | checks every fixture value against an expectation computed without the generator, the provenance fields and runtime states of the data files, the equation and register cross-references, and SHA256SUMS |
| `independent_physics_checks.py` | mathematical checks of the identities used by this specification by numerical quadrature, numerical differentiation and seeded stochastic simulation; it does not check transcriptions |
| `SHA256SUMS` | digests of all other files |

All transcribed numbers are stored as strings exactly as printed, including "≤", "≥", "~", ranges and uncertainties, so that no rounding or reinterpretation enters the data. Every row or object names its source key (mapped to [Pn] by `source_key_map.json`), its location in the source and its verification label, and lists the register items (D-n, U-n of Section 19; MFP v1.2 items with their prefix) that apply to it. The derivation records are reproduced as written during the review, under a header that explains their local identifiers; where a record and this document differ, this document governs.

**Implementation requirement.** A loader converts a printed string to a number only when no interpretation is needed: a plain decimal or scientific-notation number, or an exact rational literal of the form integer/integer (for example "5/3") in a field declared dimensionless. Any field containing a qualifier, range, uncertainty, unit-only normalization or symbol is surfaced to the caller rather than converted (MFP v1.2, Section 21.3).

## 25. Verification and revision record

### 25.1 Version

| Version | Date | Content |
|---|---|---|
| 1.0 | 13 October 2026 | First issue. Literature searched to 12 October 2026. |
| 1.1 | 13 October 2026 | Revision after an external physical review of version 1.0 (Section 25.6): correction of the stationary-distribution statement of Section 12.3 and of the focusing statement of Section 20; the exact long-time variance of the constant-focusing process (Eq. (115)) with fixtures and checks; the printed sign inconsistency of the SOLPACS spectrum interpolation (D-128); admissibility, `provides` and `domain` fields in the registry with the rule of Section 4.1 (DYN-STRAUSS17-PRINTED to `REFERENCE_DATA_ONLY`; WAVE-SPA22-BW and SGW-COUPLE-VK01 to `REQUIRES_SOURCE_OR_CODE_AUDIT`); the SI layer of Section 1.6; the clock of a coefficient (Eq. (173)); the step-rule stall at μ = ±1; the finite-volume reference operator (Eq. (178)) and the time-dependent-rate test (Eq. (177)); the composed sampled measure and coefficient–kernel readiness in Sections 16.11 and 17; the host-integration stage PS06b and revised gates in Section 26; six new fixtures; `fixtures.py` no longer overwrites the shipped fixtures by default. No source transcription of version 1.0 was changed except the re-reading of the SOLPACS grid formulas (Section 13.7). |
| 1.2 | 13 October 2026 | Implementation review of v1.1: full scalar/tensor sampled-measure and boundary contracts; state-dependent clock/generator correction; trajectory streaming and retained adaptive noise; frame/gyrophase limits; four-population joint split; channel-resolved energy/momentum bookkeeping; typed interfaces and domain-qualified readiness; clarification of phase-space versus linear-density advection; 20 detailed AMPS work packages and 21 acceptance cases. Canonical source equation numbers/model IDs and data are not regenerated. Separate supplementary mathematical checks and compile-only contract evidence; matching v1.1 data and actual AMPS integration remain unverified. |

### 25.2 How the content was obtained

1. **Readings.** The sources were located by literature search and read in the accessible renderings (arXiv PDF and HTML, ar5iv, publisher pages, NASA ADS, repositories). Each reading records the equations with their source numbers, the definition of every symbol as stated, the parameter values with their locations, apparent printed inconsistencies, and the verification level of Section 1.3. There are six topical sets: transport equations and SDE methods; quasilinear theory and waves; nonlinear theories and simulations; empirical SEP shapes; GCR models; waves and shocks. Two follow-up sets re-read the transport equations and the D_μμ forms.
2. **Sections.** Sections 2 and 5–16 were drafted from the readings, one group of sections at a time. Each group has a working record (DR-02 to DR-16, Section 24) that contains the scripts checking every derived statement: mpmath at 40 digits, and numpy Monte Carlo with fixed seeds. Sections 1, 3, 4, 17 and 20–26 were written last, from the drafted sections.
3. **Register.** Section 19 was compiled from the candidate items of the working records. Duplicates were merged, and items withdrawn by their own checks were dropped with the reason (`register_map.json`).
4. **Data and fixtures.** The parameter files were built from the data rows of the working records and compared with the section texts. `fixtures.py` generates the 89 fixtures at 40 digits and runs an independent cross-check for each: 380 cross-checks, none failed. The largest relative difference is 6.7 × 10⁻⁵, in the finite-volume cross-check of F-FOC-04 (400 angular cells, tolerance 10⁻⁴); the next largest are 4.6 × 10⁻⁷ (a limit evaluated at μ = 1 − 10⁻⁶, F-PUSH-19) and 1.8 × 10⁻⁸ (F-NL-05); all others are below 2 × 10⁻¹².

### 25.3 Independent checks of version 1.0

Each check was made by a separate checker that had not produced the content it checked.

1. **Second reading of the sources.** Two checkers re-read the pusher-critical and convention-critical equations of 20 sources against the transcriptions. The sources were: Effenberger & Litvinenko 2014; Merten & Aerdker 2024; Zhang et al. 2023; Strauss & Effenberger 2017; Wijsen et al. 2019; Sandroos et al. 2014; van den Berg et al. 2026; Tenishev et al. 2022; Kozarev et al. 2013 and Young & Poduval 2025; Lang et al. 2024 and Lavasa et al. 2026; Pacheco et al. 2019; He & Wan 2015; Ivascenko et al. 2016; Tautz 2013b; Chen et al. 2024; Kubo et al. 2015; Riordan & Pe'er 2019; Florinski 2013; Vainio & Kocharov 2001. Most transcriptions were confirmed. The following were changed as a result:
   - He & Wan 2015 Eq. 3 prints two forms of the focusing term with opposite signs. D-1, U-1 and Sections 2.1, 2.2 and 16.10 were revised; the earlier text attributed only the second form to the paper.
   - New register items: D-120 (Sandroos et al. Eq. 36 prints μ̃₀), D-121 (the middle expression of Strauss & Effenberger Eq. 25), D-122 (the exponent of Tautz 2013b Eq. 18) and D-123 (the focusing sign of Kozarev et al. Eq. 1).
   - Zhang et al. 2023 Eq. 13 was confirmed (D-116 resolved).
   - The Merten & Aerdker fixed-point statement was re-read (D-118, Section 16.9).
   - A quotation was re-attributed to the Àgueda thesis (Section 16.7).
   - Florinski 2013 was recorded as the single-author preprint actually read (Section 15.6; reference list).
   - The frame of μ for Ivascenko et al. and Lange et al. is now "not stated" (Section 11).
   - Several verification levels were aligned.
   - A third reading of van den Berg et al. 2026 Eq. 17 did not show the √ glyph that two earlier readings recorded (Section 15.2).
2. **Fixtures recomputed.** Two checkers recomputed all 567 stored fixture values from the formulas of this document. They did not use the generator. All values agree within the tolerances of Section 18, after two kinds of fix:
   - One generator error was corrected. F-SGW-02 stored the raw second moment of the predictor–corrector step as its "variance"; the corrected values are 0.999997549880043 (η = 0.5) and 0.99999658130875 (η = 0.3).
   - Several faults in the checkers themselves were fixed. The recomputation is shipped as `reference_verification.py` (Section 24).
3. **Independent mathematical checks.** A third checker wrote `independent_physics_checks.py`. It checks the identities and the Monte Carlo statements of Sections 3, 5, 6, 12 and 16 by quadrature, numerical differentiation, finite-volume solutions of the pitch-angle equation and seeded Monte Carlo. Two statements were found false and corrected:
   - **Section 16.4.** The largest exit probability per step is 1 − Φ(√(2 − νΔt)), not 1 − Φ(√2).
   - **Section 16.6.** The Euler–Maruyama reflection bias is about three times the no-boundary value, not four. The table entry at νΔt = 0.01 is now −0.0034, the mean of two independent time-averaged runs (−0.0033 ± 0.0001 and −0.0034 ± 0.0002). The first run had given −0.0043, 3.3 standard errors away.

   Three error statements for double-precision evaluation were also restated (Sections 5.6 and 12.7).
4. **Review of the text.** Four review passes covered Sections 1–21, 24 and 26. They led to these corrections:
   - TS2003 asymptotic rows (Section 7.3).
   - The range of the SOQLT Gaussian-to-step ratio (Section 10.6, U-38).
   - The ballistic exponent of van den Berg et al. 2024 (D-124).
   - The residue in Section 15.2.
   - The μ dependence of the rotation-kernel second moment (Section 14.6).
   - The bound in C-13.
   - A rounding in the Section 16.6 table.
   - Fixture and equation cross-references.
   - The notation for source equation numbers (Section 1.4).
   - Three open questions were recorded: the sign of the van den Berg et al. 2026 weight (D-125), the sign of the Vainio & Kocharov wave update (D-126) and the units of Spanier et al. Eq. 38 (D-127).
   - The runtime states of SGW-COUPLE-SLR19 and SDE-W-VDB26 were changed to `REQUIRES_SOURCE_OR_CODE_AUDIT`.
5. **Bibliographic metadata.** All 147 references were checked against OpenAlex, arXiv abstract pages and publisher pages. Crossref was used where it answered; some lookups were refused for rate limits. The result was 91 confirmed, 43 corrected and 13 partly confirmed. One of the partial records, Florinski 2013, was then replaced by the preprint actually read. Twelve records keep a metadata note in Section 22. Corrections include missing co-authors (Qin et al. 2005), years (Shalchi, Lazarian & Schlickeiser 2007; Florinski et al. 2012), titles, issue numbers and DOIs.
6. **Bibliography comparison.** Section 23 compares the cited references with the supplied bibliography. Problems found in the supplied records are listed in Section 23.3.

### 25.4 What rests on a single reading or was not verified

- Items labelled Full text ×1, "one rendering" or "summarized reading" in Sections 2–16, in particular:
  - the Section 3 equations of van den Berg et al. 2026 (PDF text layer);
  - Sandroos et al. 2014 Eqs. 27–43 (PDF only);
  - Tenishev et al. 2022 (arXiv v1 only);
  - the Zhang et al. 2023 equations (ar5iv);
  - Ruffolo 1995 Eq. 11 (assembled by the extraction tool);
  - the partly garbled tables and equations listed in Section 19.
- Papers read as abstract only or not reached (Section 19.3).
- The open printed discrepancies and decisions of Sections 19.1 and 19.2. These include D-125 (weight sign in van den Berg et al. 2026), D-126 (wave-update sign in Vainio & Kocharov 2001) and D-127 (units of Spanier et al. Eq. 38), which need the typeset pages.
- The twelve references with partly confirmed metadata (notes in Section 22).
- References into the companion specifications (MFP v1.2, PARALLEL rev. 1.4). They were checked against the companion files at writing time but are not rechecked by the scripts of the bundle.

### 25.5 Companion bundle check

**Historical v1.1 author report, retained for provenance.** The following results were stated in the uploaded description. They were not repeated against a matching v1.1 archive in this review, because only the v1.0 archive was supplied. See Section 25.7 for the current evidence boundary.

The bundle PITCH_ANGLE_SCATTERING_MODEL_DATA.zip (23 files) of version 1.1 was extracted to an empty folder and checked there:

- `reference_verification.py`: PASS. All 652 stored fixture values of the 89 fixtures were recomputed independently and agree within the tolerances of Section 18. The six data files pass the provenance, source-key, model-ID, runtime-state, section, equation-label and register-identifier checks (D-1 to D-128, U-1 to U-67, C-1 to C-15), and all 22 digests of `SHA256SUMS` match. Run time is about 22 s.
- `independent_physics_checks.py`: 274 of 274 checks PASS (253 of version 1.0 and the 21 checks of group 13 added in version 1.1). Run time is 291 s. The Monte Carlo run MCR-8 quoted in Section 16.9 is this run.
- The SHA-256 digest of the file `SHA256SUMS` is f75ca823a10bb129d415c0816612e070ece34bbe5ff630c9cabbaed3c741b9a0.
- Version 1.0 (SHA-256 of `SHA256SUMS` aba757e4c214bca1987c5f1a4014359315b8bc8e576a2fd6d39fee9f72dcfd13) had 83 fixtures with 567 values and 253 checks; the version 1.0 fixtures are unchanged in version 1.1 (the six new fixtures were appended).

### 25.6 External review of version 1.0 and the changes of version 1.1

Version 1.0 (SHA-256 of the specification a95bfee7aa60621d7305aa467e021e586b41fe344a169b5e5e762b984b56ac70, of the data bundle 56991409de20142cca6abf29f6651258ed7eb20ca34f7085c1178550c5a64b6e) was reviewed externally against the supplied bundle; the reviewer ran `reference_verification.py` and `independent_physics_checks.py` (all 567 fixture values and 253 checks passed) and reported findings with their own calculations. Each finding was re-derived here before any change was made; the list below records which findings were confirmed and what was changed, and which were not adopted and why.

**Confirmed and adopted.**

1. *Spectrum interpolation sign (Section 13.7).* The two printed definitions of the SOLPACS local index and interpolant are inconsistent. The source was re-read in two renderings (PDF text and ar5iv), which both print q_j = ln(I_w,j+1/I_w,j)/ln(k_j+1/k_j) with the exponent −q_j; the inconsistency is therefore in the source as printed, not in the version 1.0 transcription, and version 1.0 had reproduced it without comment. Recorded as D-128 with the node test F-MEAS-03; the derived check DR-06 C13 had assumed the power-law form and could not expose it. The spectrum provider stores one convention and reproduces both nodes of every interval (PS13).
2. *Stationary distribution (Section 12.3).* Version 1.0 called e^{ξμ} "the zero-flux stationary distribution of Eq. (1)". A z-independent e^{ξμ} is not a solution of that equation (residual (vμ/L)e^{ξμ}, F-FOC-05); it is the stationary density of the Fokker–Planck equation for the linear density, i.e. of the forward process with reflection, and the corresponding solution of Eq. (1) is f₀ ∝ e^{−z/L}e^{ξμ}. The backward process samples e^{−ξμ}. Section 12.3, 15.7, 16.9 and the roadmap (PS15, configuration E, U-67) now state which observable is compared with which reference, and the manufactured solution is a fixture (VAL-PUSH-MMS).
3. *Focusing closure versus measured variance (Sections 12.3, 20, 15.7).* Version 1.0 presented Eq. (112) as "the κ∥ produced by a pusher". Eq. (112) is the diffusion-approximation coefficient of Effenberger & Litvinenko; the long-time centred spatial variance of the constant-L isotropic process is Eq. (115), re-derived here from the Poisson equation of the stationary generator and checked by 40-digit quadrature, a finite-volume Poisson solve (0.178275493 and 0.178275762 with 1000 and 4000 cells against 0.178275780 at ξ = 1.5) and a Monte Carlo run of the SDE (MCR-8). The two coincide only as ξ → 0; at ξ = 1.5 the variance is 39 % below the closure. The mean streaming speed is the same in both. Wang & Qin (2019) [P120] is cited for the distinction (Abstract only). Fixture F-FOC-04; validation reference VAL-PUSH-KVAR; the sign symmetry of ξ (U-44) was derived at the same time.
4. *Inadmissible printed forms (Sections 4.1, 7.4, 7.7, 17).* A runtime state recorded whether inputs are specified, not whether the printed form is a physical coefficient; a dimensionally inconsistent kernel (D-44) was executable as a named variant. The registry now carries an `admissibility` field and Section 4.1 the rule that an inadmissible form is never evaluated as a coefficient of Eq. (21); DYN-STRAUSS17-PRINTED is a record (`REFERENCE_DATA_ONLY`, U-25 resolved), and the models with open audit items D-126 and D-127 (SGW-COUPLE-VK01, WAVE-SPA22-BW) are in `REQUIRES_SOURCE_OR_CODE_AUDIT`, the audit blocking before any user decision. The printed endpoint variants (D-32, D-33) are "unestablished" and run only under their own IDs for reproduction or diagnosis, as before.
5. *Registry capabilities (Section 4).* IDs that are not pitch-angle coefficients (resonance kernels for diffusion in v∥, spatial coefficients, growth rates, focusing closures, kernels, references) were counted with the D_μμ laws. The field `provides` now separates them; 73 of the 172 IDs are D_μμ laws or pitch-angle tensors.
6. *Domains of asymptotic coefficients (Section 9.1).* The printed Alfvénic, slow-mode and fast-mode asymptotes are not real for μ < 0 and Eq. (77) diverges at μ = ±1; version 1.0 listed them as `READY_EXPLICIT_INPUTS` without a domain. The registry field `domain` and U-66 (symmetry extension and endpoint treatment, not printed by the sources) were added; the QLT application in the sheath at δb > B₀ (Section 15.4) is labelled an application of the source outside the weak-turbulence assumption.
7. *Coefficient–kernel compatibility (Sections 16.5, 16.11, 17, PS03–PS05).* Readiness is now a property of the (coefficient, kernel) pair; `dD_dmu` may be absent, and a derivative-free finite-volume reference with reciprocal-coefficient transmissibility (Eq. (178), F-PUSH-21) is the reference for coefficients with unbounded or discontinuous drift. The statement that a finite reciprocal-coefficient integral across μ = 0 (α < 1) connects the hemispheres, while an interface value D_μμ(0) = 0 would disconnect them, was derived and is consistent with Section 16.5.
8. *Step rule at μ = ±1 (Section 16.6).* Eq. (165) returns exactly zero at μ = ±1 for every shape with D′_μμ(±1) ≠ 0 (F-PUSH-19); U-55 was extended, and residence-time weighting of adaptive trajectories was added as a requirement.
9. *Clocks, SI layer and interface (Sections 1.6, 16.8, 16.11, 17).* The rate conversion dt′ = γ_u(1 − uvμ/c²)dt (Eq. (173)), the declared clock of each coefficient, the composed sampled measure in place of the single enumeration `SampledDensity` (which could not express A(s)p²f), a `PitchStep` that returns the full particle state and the plasma-frame energy exchange per wave population, the positive-semidefiniteness of the (μ, p) diffusion tensor and a normative SI table were added. The statement that total Alfvén-wave energy is δB²/μ₀ only under the equipartition assumption, against δB²/(2μ₀) for the magnetic part, is in Section 1.6.
10. *Host integration (Section 26).* Version 1.0 asked the implementer to inspect the host but stated no acceptance tests for a three-dimensional adaptive mesh or for segmented field lines. Stage PS06b (geometry, time and sampling adapters), the time-dependent-rate test Eq. (177) (F-PUSH-20, VAL-PUSH-LEG), the rule against duplicating a transport process (parallel spatial diffusion from the companion κ∥ library for a population already scattered in μ; focusing added to a full-orbit pusher) and the conservation tests for particle–wave coupling were added, with D-126 blocking PS14 until audited.
11. *Estimators (Section 11.1).* For a differentiable orbit the literal Δt → 0 limit of ⟨(Δμ)²⟩/(2Δt) is zero; the estimator window and the distinction between an orbit-code running estimator and the transition solution of a Markov pusher were added.
12. *Test policy and data package (Sections 18, 24, 26).* `fixtures.py` no longer overwrites the shipped `benchmark_points.json` by default (PS01); error budgets and statistical bounds are fixed before a test (PS16); Section 1.3 states that repeated readings of one text extraction do not confirm the typeset original.

**Not adopted, with the reason.**

- The reviewer's proposed replacement stage sequence (R00–R13) was not substituted for PS00–PS16, because the PS identifiers are referenced throughout the document and the data; its content was merged into the existing stages and PS06b, and a correspondence is given in Section 26.4.
- The suggested first production subset (isotropic and ε-form scattering first, the q-form later) is an ordering of the implementation work and is left to the user; this document does not rank the models (Section 20). The compatibility tests that the suggestion rests on are required for every pair regardless of order.
- The suggestion to add a document hash and a typeset equation capture for every published formula is not met here: the readings used text renderings (Section 25.2), and this is stated rather than claimed.

**Unchanged.** The reviewer confirmed the central operator and SDE, the normalization Eq. (23) for even and non-even D_μμ, the isotropic relations, the exact normalizations of the q- and ε-forms, the factor-of-two adapters, the distinction between angular diffusion and Poisson redraws, declared spectral sidedness and the rule against counting wave-frame energy changes twice. No source transcription of version 1.0 was changed by this revision, except that the SOLPACS grid formulas were re-read in a second rendering (Section 13.7).

### 25.7 Implementation review of version 1.1 and evidence of version 1.2

Reviewed input: `PITCH_ANGLE_SCATTERING_MODEL (1).md`, specification v1.1, dated 13 October 2026. Revised deliverable: v1.2. This is a scientific specification/implementation-contract revision, not a completed AMPS library or a certification of every literature model.

The update adopts the earlier corrections to spectrum interpolation, focused stationary densities, the distinction between a diffusion closure and spatial variance, dimensionally invalid printed forms, singular-drift kernel compatibility, SI clocks and host geometry. Those corrections should be retained. The main remaining improvements are the equations and interfaces that determine what a pusher actually samples, and explicit gates for implementing them in AMPS.

| Priority | Remaining issue in v1.1 | v1.2 improvement |
|---|---|---|
| Required before affected release | Tensor capability advertised through D and D′ scalar accessors | Separate scalar/tensor interfaces; full operator, p²-measure drift, cross noise and log-p correction; Sections 16.12/17 |
| Required before frame-clock release | Clock conversion stated, but state-dependent generator/drift conversion not explicit | Eq. (N7): multiply the whole generator; do not differentiate gD as an unchanged-measure coefficient |
| Required before 3-D exact-boost claims | μ,p-only state cannot supply a nonparallel vector boost | Keep gyrotropic state in the local plasma frame, or require full vector/gyrophase and the correct transformed equation |
| Required before finite-step transport release | Angular endpoint tests do not determine the displacement integral | Trajectory streaming/dwell integration and finite-time spatial tests; PS04a |
| Required before adaptive rejection | A rejected stochastic kick can be discarded and redrawn | Brownian bridge with retained unused intervals; checkpoint and migration ownership |
| Required before tensor boundary release | Scalar mirror rules do not automatically implement mixed tensor no-flux | Full angular/momentum boundary contract; coupled flux and conservation tests |
| Required before domain-limited release | Section 4 lists only three admissibility labels, while the enum/counts use four; interface rule rejects every non-Admissible value | Explicit fourth value; domain-qualified execution permitted only with complete problem coverage |
| Required before registry-backed certification | Available archive remains v1.0, while the description claims v1.1 data/fixture checks | Version pinning; historical claims distinguished from checks actually run; no regenerated canonical archive claimed |
| Required before wave feedback | Direction-only exchange key cannot identify polarization, resonance channel and bin | Event/population/bin/epoch keys, energy and momentum transfers, normalized deposition from the same kernel |
| Required before four-population wave model | Two helicity ratios leave the joint split undetermined | Explicit joint moment χ and nonnegative-weight bounds, or four supplied spectra |
| Required before normalization/table release | Missing band, physical zero, divergent λ and failed quadrature can be conflated | Typed status results, canonical spectra, tagged unfocused λ∥ and consistent value/derivative interpolation |
| Required before AMPS handoff | Broad stage table lacks concrete work packages and measurable budgets | 20 work packages, 21 acceptance cases, three numerical release scopes, commands/ledger requirements and lifecycle tests |

The tensor drift follows directly from the physical operator for gyrotropic f. For particles sampling P=p²f, it includes ∂pDμp+2Dμp/p in the angular drift and ∂μDpμ+∂pDpp+2Dpp/p in momentum drift. Its noise covariance is 2D. This is the measure-specific connection missing from the scalar sketch; merely returning a positive tensor is insufficient.

A further mathematical issue is the state-dependent clock factor g=dt_c/dt_h. For the original coefficient-frame generator L_c, the host path-clock generator is gL_c. A scalar source drift D′ therefore becomes gD′; replacing it with (gD)′ introduces g′D. The revised contract distinguishes a time change from a new divergence operator and requires transformation of the full frame/measure equation.

The prescribed first release uses individually certified scalar models and AMPS adapters. Isotropic scattering provides exact moment and finite-time displacement references; q-form and ε-form kernels must pass their own singularity/occupation/transport tests. Spectral models and wave feedback follow separate gates. A mathematically admissible prescribed diffusion law is not evidence that it uniquely describes measured SEP/GCR turbulence: physical approximation domains and observational comparisons remain model-specific.

The supplementary checks independently verify tensor sampled-density identities and covariance, Lorentz mass shell/elastic energy exchange, helicity joint-weight bounds, conservative nonuniform finite volumes, focused Poisson references, manufactured solutions, trajectory-integral distinctions, Brownian bridges, spectrum interpolation, optical-time moments and clock-induced drift. `verification_results.json` records 262/262 passing checks. `pitch_angle_contract.hpp` compiles with C++17 warnings treated as errors. These checks do not run an AMPS pusher, implement the declared functions, or certify the missing v1.1 archive.

The available archive was inspected: its README states v1.0, its registry contains 168 IDs, and its fixture file contains 83 cases. Its `SHA256SUMS` digest matches the v1.0 digest quoted in v1.1. The v1.1 claims of 172 IDs, 89 fixtures and 652 values are therefore historical author-reported evidence awaiting the matching archive. The canonical v1.0 files were not overwritten or relabeled.

The revised model is ready to guide implementation through explicit contracts and tests. Production readiness still requires the actual AMPS checkout audit, matching source data, implemented kernels/adapters and the advertised release's convergence campaign. Source audits D-125/D-126/D-127 and delegated companion equations remain capability-specific blockers.

Primary references checked during this review:

1. [AMPS official repository](https://github.com/SWMFsoftware/AMPS): public project layout only; the target checkout/pusher was not audited.
2. [Rackauckas & Nie (2017)](https://www.aimsciences.org/article/doi/10.3934/dcdsb.2017133), DOI 10.3934/dcdsb.2017133: rejection sampling with memory for adaptive SDE integration.
3. [Effenberger & Litvinenko (2014), author-hosted publisher PDF](https://researchcommons.waikato.ac.nz/bitstreams/7c28952e-22d1-4f50-b560-8e04162a8d0f/download), DOI 10.1088/0004-637X/783/1/15: distinction between the full focused-transport equation and diffusion/telegraph approximations.
4. [Impact of evolving turbulence on diffusion coefficients, Appendix A1](https://academic.oup.com/mnras/article/547/4/stag393/8499610), DOI 10.1093/mnras/stag393: an accessible primary example of the full gyrotropic momentum/pitch-angle scattering operator. The SDE/time-change identities here are derived and independently checked, not quoted as that paper's numerical results.

All runtime equations added in Sections 16.12/17 are labeled as derived contracts. The review does not claim a fresh typeset-source audit of all 147 historical references.

## 26. Codex implementation roadmap (final section)

This roadmap preserves PS00–PS16 and PS06b so existing cross-references remain useful. It adds PS02a for the tensor contract and PS04a for trajectory/time integration. Sections 16.12 and 17 specify executable behavior; the source transcriptions, their discrepancy register and approximation domains remain the scientific record. Each stage below is a work package with inputs, deliverables and a falsifiable exit gate. An unavailable source or reference is a blocker for the affected capability, never a passing test.

### 26.1 Operating instructions and release scope

**Implementation objective.** Build one mesh-independent pitch-angle model library and thin adapters for AMPS field-line and 3-D transport applications. First release a validated prescribed-scattering subset; add turbulence-derived scattering and particle–wave coupling through separate gates. This is an engineering order intended to expose convention and pusher errors early. It does not select the physical model for a scientific application or rank published theories by validity.

| Release | Supported scope | Required evidence |
|---|---|---|
| A: prescribed elastic scattering | PA-ISO and individually certified PA-EPS/PA-QFORM pairs; constant or supplied λ∥; selected forward transport and, where implemented, backward reconstruction; no wave evolution | Core, singularity, time/streaming, geometry, frame and manufactured-solution gates; explicit supported-pair list |
| B: prescribed turbulence | Canonical spectra; selected magnetostatic/dynamical slab and nonlinear laws on certified domains | Release A plus spectrum normalization, resonance/quadrature, PSD/domain and table-convergence evidence |
| C: particle–wave interaction | Identified wave populations, elastic-frame or full-tensor kernel, particle energy changes, selected self-generated-wave/shock coupling | Release B plus channel-resolved transfers, closed-box energy/momentum, positivity and growth-limit tests |
| Application validation | A named SEP or GCR setup with its source/initial/boundary data | Numerical release for its exact feature set, independent transport comparison and, if claimed, comparison with observations |

PA-ISO with a bounded small-angle rotation is a useful first forward implementation because it has exact angular and spatial references and no internal cusp. PA-QFORM is not paired with a derivative-based kernel by default when $1<q<2$. Use the derivative-free reference and release only a specifically certified implementation. PA-EPS with ε>0 is a useful subsequent anisotropic case, with its discontinuous drift treated explicitly. A coefficient may be evaluable before any production kernel is released for it.

**Host audit before edits.** The public AMPS repository lists `srcSEP`, `srcSEP3D`, `src` and `test`, but the actual target checkout and revision are not supplied with this description. PS00 must inspect that checkout, its instruction files, build rules and particle layout. No API, model directory, input keyword or compiler upgrade below is asserted to exist in that checkout. Reuse shared model infrastructure when present. Keep the physics in one location and avoid separate copies in the SEP and GCR pushers.

**Proposed module boundaries, to be mapped at PS00.** Relative paths are design labels until the checkout is inspected.

| Logical module | Proposed files under the approved shared model directory | Responsibility |
|---|---|---|
| Contract and registry | `pitch_angle/contract.*`, `registry.*`, `configuration.*` | SI types, diagnostics, domain/readiness, provenance and version checks |
| Prescribed shapes | `pitch_angle/shapes.*`, `normalization.*` | D, analytic derivatives, amplitude and tagged λ∥ |
| Spectral physics | `pitch_angle/spectrum.*`, `resonance.*`, `slab_qlt.*`, `dynamic_qlt.*` | Canonical variance, wave populations and converged integrals |
| Numerical kernels | `pitch_angle/kernels.*`, `fv_reference.*`, `tensor_kernel.*`, `frames.*` | Pure scattering updates and independent reference solver |
| Transport composition | `pitch_angle/transport.*`, `step_scheduler.*` | Selected complete FTE, streaming/focusing/momentum terms and physical interval ownership |
| Host adapters | Existing SEP/SEP3D adapter locations selected at PS00 | AMPS state, cell/segment sampling, snapshot, weights, migration/restart and wave deposition |
| Verification | Existing test system plus `pitch_angle` cases and frozen data | Independent references, stage ledger, machine-readable results and release matrix |

Prefer immutable configuration and snapshot views, explicit errors and model dispatch selected outside the particle loop. The C++17 contract in the supplement is a compile-checked declaration, not an AMPS implementation. If the host uses another language standard, adapt optional/result types locally rather than changing the entire project solely to adopt the sketch.

**Stage records.** Maintain a ledger with `NOT_STARTED`, `IN_PROGRESS`, `PASS`, `FAIL`, `BLOCKED` and `SKIP`. Each entry includes code revision, specification/data digests, configured model/pair/domain, changed files, exact command, reference identity, measured errors and uncertainties, output locations and outstanding D-/U-items. Record skipped optional features in the released capability matrix. Baselines and canonical fixtures remain read-only; regeneration goes to a separate path.

**Test ownership.** Separate the coefficient implementation from the reference computation. Compare analytic moments or an independently assembled flux operator, rather than running the same coefficient helper on both sides of a supposed validation. Frozen fixtures check conventions and algebra; they do not establish physical validity of a turbulence approximation. The supplementary Python checks are mathematical evidence for this specification and are not evidence of successful AMPS integration.

### 26.2 Stage dependencies and deliverables

| Stage | Depends on | Main deliverable | Release gate |
|---|---|---|---|
| PS00 | Target checkout | Host process/state/time inventory and frozen baseline | Every μ/p/time/weight update mapped; build and baseline reproducible |
| PS01 | PS00 | Versioned data loaders and provenance/schema validator | Matching archive passes; unsafe or incomplete configurations refused |
| PS02 | PS01 | Scalar standard operator, kinematics and λ∥ API | Convention and normalization references pass |
| PS02a | PS02; tensor capability requested | Tensor operator/measure contract | Eq. (N2) equals the sampled-density operator using Eq. (N3); PSD and covariance pass |
| PS03 | PS02 | Prescribed coefficient evaluators and presets | Round-trip normalization, zero/endpoints and domain checks pass |
| PS04 | PS03 | Scalar kernels and conservative angular reference | Number, stationary density and angular convergence pass |
| PS04a | PS04 | Physical-time scheduler and trajectory streaming | Finite-time spatial references, Brownian bridges and nonprogress handling pass |
| PS05 | PS04, PS04a | Certified coefficient–kernel pairs | Cusp/barrier/occupation/transport error budgets pass |
| PS06 | PS04a | Alternative kernels as separate capabilities | Rotation, hemispheric split and redraw each match their own operator |
| PS06b | PS00, PS05; PS06 when used | Cell and field-line adapters | Geometry/time/measure/migration/restart cases pass |
| PS07 | PS06b; PS02a when used | Selected FTE, frames and reconstruction weights | Forward/backward manufactured problems and frame/clock references pass |
| PS08 | PS04a, PS07 | Focusing diagnostics and closure separation | Stationary density, drift and Eq. (115) pass independently |
| PS09 | PS02; pinned spectral dependencies | Canonical spectra and slab QLT | Variance, prefactor, resonance and finite-band tests pass |
| PS10 | PS09 | Selected dynamic/dissipation models | Width, narrow-resonance and limiting cases converge |
| PS11 | PS07, PS09; PS02a where needed | Wave populations and energy-changing scattering | Charge/direction/polarization and energy/frame tests pass |
| PS12 | PS09; PS02a where needed | Selected compressible/nonlinear models | Full integral/asymptote distinctions and certified domains pass |
| PS13 | PS09, PS04a | Measured-spectrum adapters and estimators | Node/variance/coverage and finite-lag tests pass |
| PS14 | PS11, PS06b; host wave/shock capability | Conservative wave coupling | Closed-box, positivity and growth-limit tests pass |
| PS15 | PS07, PS08; features needed by chosen GCR law | GCR application adapter and boundary problems | Reconstruction and independent application reference pass |
| PS16 | All stages needed by advertised release | Regression, performance and methods handoff | Supported-pair matrix and reproducible evidence complete |

Dependency rules apply per feature. PS02a, PS09–PS14 or a backward pusher are not prerequisites for a purely scalar forward Release A. PS16 is run for each release's selected subset, not deferred until all 172 inventory records are executable. A missing companion equation blocks only the path that delegates to it, with that dependency identified.

### 26.3 Detailed stage work packages

#### PS00 — Host baseline and process audit

Inputs: actual AMPS checkout/revision, local instructions, build configuration and the intended SEP/GCR transport setup. First build and run the existing baseline. Use repository search to locate all μ, momentum, velocity, focusing, parallel displacement and statistical-weight updates, local timestep selectors, sampling callbacks and RNG calls. A starting search is `rg -n 'Dmumu|D_mu|pitch|Scattering|focusing|FieldLine|TimeStep|Random'` in the audited source paths; inspect callers, not only names.

Create a table for each pusher: equation solved, forward/backward direction, sampled density, coordinates/frame, clock, weight meaning, deterministic terms, stochastic terms, angular boundary rule, geometry exits and restart state. Map every process to one owner. Identify whether focusing is resolved by an orbit/ guiding-centre mover or inserted as a drift, and whether a κ∥ module already supplies spatial noise. Document when particles see MHD/turbulence updates, whether field-line segments move, and how local clocks are synchronized. Inspect existing test/build targets and compiler capabilities before choosing files.

Deliverables: baseline command/configuration/output; process-ownership matrix; host-state conversion contract; approved module locations; stage ledger and test inventory. Gate: the baseline repeats, all relevant update sites are accounted for, and no physical term can be added twice. Unmapped behavior is documented as a blocker for that adapter. Do not diagnose a host defect from the printed Tenishev equation alone.

#### PS01 — Data, versions, schemas and provenance

Inputs: this specification, the exact companion archive and any delegated MFP/PARALLEL texts. Confirm archive version and digest before loading presets. The currently available archive is v1.0: 168 IDs, 83 fixtures, 567 values, and SHA-256 of `SHA256SUMS` `aba757e4c214bca1987c5f1a4014359315b8bc8e576a2fd6d39fee9f72dcfd13`. The v1.1 description claims 172/89/652 and digest `f75ca823a10bb129d415c0816612e070ece34bbe5ff630c9cabbaed3c741b9a0`; that matching archive was not supplied for this review. Do not pair its claims with the v1.0 files.

Implement read-only loading and strict schemas for all four admissibility states, capabilities, domain/extension, clocks, units, parameter qualifiers and register items. Plain literals and exact dimensionless rationals may be parsed; ranges, uncertainties and qualitative statements remain structured source data. Qualify registry readiness as source completeness, and keep kernel/host certification separately. Supply a versioned migration only after inspecting the actual v1.1 registry; do not infer shifted equation labels from v1.0.

Gate: matching archive checksums and all supplied verification scripts pass; explicit-output fixture regeneration agrees within its declared comparison policy and leaves original fixtures unchanged; schema rejects inadmissible formulas, wrong capability, unknown conventions, open source audits and unresolved units. `AdmissibleOnDomain` succeeds inside its domain and fails outside. Deliver machine-readable validation and a data manifest. Archive certification remains BLOCKED until the matching data are available; independent normalized core checks can still proceed with their own declared provenance.

#### PS02 — Scalar operator, SI and normalization

Implement Eq. (21) as the sole scalar convention. Convert Ruffolo/EPREM coefficients once at model construction. Create authoritative particle kinematics from species and SI momentum; validate μ and field direction without silently repairing grossly invalid input. Derive signed Ω, positive maximal r_L, speed, γ and stable kinetic energy from that state. Implement PA-ISO with ν=2D₀, λ∥=v/ν and κ₀=vλ∥/3 only in the unfocused elastic diffusion limit.

Implement Eq. (23) with endpoint/internal singularity analysis and typed finite/divergent/unavailable/nonconverged results. A zero scattering amplitude has divergent λ∥; it is a ballistic case, not an exception to hide with a floor. Use scaling to avoid subtractive cancellation in closed forms and quadrature split points at nonsmooth locations. Keep λ tags distinct from radial prescriptions and focusing closures.

Gate: canonical factors of two, isotropic identities, both μ signs, endpoint limits, low/high-momentum kinematics and unit conversions pass; quadrature reproduces closed-form λ∥ with the budget in Section 26.5. Reject an incompatible λ tag. Deliver evaluator tests and a coefficient grid with units and convention metadata. Do not require absent companion equations when the selected formula is complete in this document.

#### PS02a — Tensor diffusion and sampled measure

This work package is required before any full $(\mu,p)$ tensor capability is advertised. Implement the operator of Eq. (N2), six derivatives and Eq. (N3); factor the PSD tensor to produce correlated noise, including zero/rank-one cases. If using log momentum, implement its Itô correction and explicit momentum-domain boundaries. Define the tensor's coordinates, frame and clock; apply normalization conversions to every entry consistently.

Compare the analytical physical operator with the forward coordinate-density operator acting on several nonconstant manufactured functions. A constant f under pure closed scattering is stationary in the physical measure; its sampled momentum density carries p². Verify the covariance $BB^{\mathsf T}=2D$ after coordinate scaling, reject indefinite inputs and prove that the implementation does not drop mixed derivatives. Independent scalar and momentum noise is permitted only when $D_{\mu p}=0$.

Gate: pointwise manufactured-operator residuals, first/second/cross increment moments and PSD checks pass; momentum boundaries conserve or remove number exactly as specified. For a released tensor pusher, compare a two-coordinate solution against an independent conservative solver and sweep Δt and both grids. The supplement verifies the mathematical contract, but production tensor integration remains unreleased until these implementation tests run.

#### PS03 — Prescribed shapes and presets

Implement PA-ISO, PA-QFORM, PA-EPS and selected aliases/variants through source-tagged shape records. Return D, its regularity and actual analytic one-sided derivatives where finite. Separate coefficient evaluation from kernel readiness. Check ranges of q, H, ε and amplitude; source conventions such as H versus h are not interchangeable. A preset records whether normalization is from λ∥, ν, D₀ or another source-specific amplitude, including conversion provenance.

Use closed-form normalization where available and an independently evaluated integral as the check. Include H→0, ε→0, q→2 and endpoint limiting cases with the correct finite/divergent classification. Anisotropic/directional coefficients need both hemispheres integrated explicitly. Do not repair an unbounded drift by setting its value at zero to zero. A deliberately smoothed shape has its own model identity, smoothing scale and convergence study.

Gate: amplitude→λ∥→amplitude round trips pass; positivity and real-valuedness hold on the advertised domain; derivatives agree with smooth-point numerical checks and known one-sided limits; cusp and barrier classifications are correct. Printed (1−μ)² or endpoint-singular variants remain diagnostic until their physical behavior is established. Deliver source-specific evaluator and normalization tables; no SDE release is implied by this stage.

#### PS04 — Scalar kernels and independent angular reference

Implement Euler–Maruyama and Milstein for compatible regular coefficients, including $D'$ and the standard $\sqrt{2D\Delta t}$ noise. For a selected transport composer, add the correct deterministic drift once; keep pure scattering tests separate. Implement repeated mirror folding for the configured no-flux angular boundary and log reflections. Resolve endpoint entrance rather than allowing the van-den-Berg timestep rule to stall at $|\mu|=1$.

Build the independent finite-volume generator using reciprocal-integral transmissibilities, nonuniform widths and the declared internal extension. Store rate orientation and diagonal explicitly. Verify the stochastic semigroup, detailed balance, uniform stationary density and count conservation. Use tridiagonal/matrix-action methods appropriate to the reference problem; do not rely on a dense production matrix. If sampling a chain, use valid continuous-time jumps or converged uniformization and integrate dwell times when streaming.

Gate: Eq. (174) and several Legendre modes agree under timestep refinement; the flux reference converges under angular refinement, including irregular grids, cusp transmission and pure-scattering ε=0 disconnection. Start beams at μ=±1 and at/around zero. Invalid derivative pairs are refused. Deliver kernel limits, boundary-bias measurements and reference-solver outputs. A reflection fraction alone is not a convergence certificate.

#### PS04a — Physical time, streaming and adaptive randomness

Define one owner of each particle's elapsed physical interval. At the beginning of a substep compute the next limit from scattering accuracy, streaming/cell or segment exit, deterministic focusing, momentum change, background update and shock crossing. Consume a strictly positive interval no larger than the earliest applicable limit and advance to common physical output times. If the coefficient-clock rate differs from the host clock, apply the declared factor once and record both durations.

Stream using the evolving pitch-angle trajectory. Compare endpoint, subcycled/trapezoidal and independently integrated references to quantify the spatial integration error; a bounded angular sampler does not authorize an arbitrarily long streaming step. For chain references, integrate each residence interval. If reducing a drawn stochastic trial, implement Eq. (N4), retain the unused path and preserve it across any retry or checkpoint. Draw-before-reject logic must be audited explicitly.

Gate: finite-time isotropic displacement, ballistic support $|\Delta s|\le v\Delta t$ in the static no-flow reference, Eq. (175), Brownian split covariance and exact time-dependent-rate moments pass. A forced update at an interior time produces the correct optical time. Endpoint beams progress, and failed steps return a structured nonprogress/error result without a silent kick or consumed time. Common-time and residence-time diagnostics agree where they represent the same observable.

Test state-dependent clock factors independently: multiply the whole generator according to Eq. (N7). A rescaled coefficient $gD$ followed by differentiation produces an erroneous extra drift $g'D$. Use a simple positive $g(\mu)=1+0.3\mu$ mathematical case and its stationary density proportional to $1/g$, then the actual Lorentz/mixed-frame transform required by the application.

#### PS05 — Zero pitch cosine and pair certification

For every advertised shape–kernel pair, specify the mathematical extension at μ=0 and any numerical approximation. Cover PA-QFORM for representative $1<q<2$ with H=0 and H>0, PA-EPS for ε=0 and ε>0, and any requested source-specific cusp. Run from uniform, zero-centered, near-zero and beamed starts. Compare with the independent angular reference at common times rather than only plotting the final stationary density.

Measure probability in $|\mu|<0.01$, exchange between hemispheres, several Legendre moments, finite-time angular distributions and spatial variance. Sweep timestep and angular reference grids independently. For a regularized candidate, also sweep its smoothing parameter independently; fixed smoothing can match a numerical implementation while changing the physical operator. Do not assume the isotropic Legendre decay law for a nonisotropic coefficient: its reference is that coefficient's operator.

Gate: confidence-bounded discrepancies fit the predeclared budgets, and reference error is smaller than the implementation budget. Pair certification records exact q/H/ε/energy/domain ranges and step-rule limits. A failing or unsupported pair stays unreleased while its coefficient evaluator remains usable for diagnostics. A derivative-free grid kernel may be released only with its quantization and transport errors documented; it is not called an exact continuum kernel.

#### PS06 — Rotation, hemispheric split and redraw

Implement MC-ROT as the bounded small-angle isotropic approximation of Eq. (168), with its optical-time condition, source angle distribution and evaluation time documented. Validate at μ=0, near ±1 and exactly at ±1. Build MC-ROT-HEMI only for the matching ε-form with the correct inverse map and explicit process order. Variable rates and operator splitting require their own error sweep.

Implement MC-REDRAW separately as Eq. (171), using its declared speed/rate and residence-time dynamics. Its Legendre spectrum differs from angular diffusion, so matching λ∥ or the first moment does not establish equivalence. Do not use redraws as a replacement for D_μμ near zero.

If the Sandroos predictor–corrector is requested, add it as a separate source-specific kernel after U-51 and its evaluation times/calculus are resolved. Verify conditional first/second moments and the correct finite-time operator independently; a predictor–corrector label alone does not establish an Itô update or a higher order. Extend the kernel enum/certification matrix only after that gate.

Gate: rotation one-step and finite-time moments converge to their scalar operator, the hemispheric kernel matches the independently assembled ε-form operator, and redraws match the jump-process reference, including higher modes. Release each kernel under its own ID/capability. Record step bias, reflected/entrance behavior, variable-rate handling and cost. An application can use a kernel only when its specific coefficient pairing has passed PS05.

#### PS06b — AMPS geometry, weights, snapshots and lifecycle

Implement cell and field-line adapters to the same immutable context and coefficient interface. Map AMPS macroparticle state and weight to the stated number measure, and record cell volume or the tube cross-section times segment length. For a static tube, verify $A(s)B(s)$ is constant. Segment splitting changes geometry representation, not represented particle number or the physical scattering clock. A plasma-comoving/moving segment needs the full selected transport equation, not the static-tube measure substituted without a Jacobian.

Split subcycles at cell/segment transitions, shocks and background publication times. Obtain B, its orientation, focusing gradient and turbulence from consistent snapshots. Cache by all relevant physics and topology versions; invalidate on AMR refinement/coarsening, segment splitting/merging, repartitioning and wave updates. Preserve particle ID, species, momentum frame, weights, physical time, unconsumed interval, pending bridge increments and undeposited exchange records during migration or checkpointing. Define a deposition transaction so an event cannot be applied twice after restart.

Gate: AP-GEOM cases of Section 26.5 pass for uniform and radial fields, irregular cell/segment sizes, background-rate changes, forced migration and mid-subcycle restarts. Equivalent 1-D/3-D physics produces matching normalized ensemble observables, not necessarily identical individual trajectories. Coupled regridding conserves represented integrated wave energy with explicit external work. Report guiding-centre ordering in physical scales. Deliver conversion code, lifecycle serialization tests and process/snapshot ownership documentation.

#### PS07 — Complete FTE, frames and forward/backward reconstruction

Select exactly one supported transport coefficient set for each application. Derive the generator, weight and observable mapping from that equation and the sampled measure; do not assemble drift terms from papers with incompatible conventions. Implement only delegated equations whose pinned companion texts are available. Handle signed inverse focusing length and the declared μ reference consistently with magnetic polarity and flow gradients.

Implement along-field Galilean/Lorentz transformations on full authoritative state and the selected clock. An exact nonparallel transform uses a vector state; a 3-D gyrotropic application otherwise remains in the local plasma frame with the justified mixed-frame equation. Test charge, boost sign, $u\to0$, $p\ll mc$ and relativistic momenta. Never transform only μ while keeping incompatible p, v or γ.

For backward reconstruction, specify terminal/outer-boundary/source data, stopping criteria and any importance-sampling drift and loss/weight terms together. Test importance sampling against an unweighted reference for the same solution, including effective sample size. Source-printed weights with unresolved sign or constants remain blocked.

Gate: boost round trips and mass-shell/clock identities pass; manufactured stationary $f_0=e^{-s/L}e^{\xi\mu}$ with exact boundary data is reconstructed by the forward and backward implementations within their separate error budgets. A backward trajectory histogram is compared with ρ−; the reconstructed physical f is compared with the manufactured f. Deliver the complete generator and observable derivation alongside code, not just a list of coefficients.

#### PS08 — Focusing, stable references and physical observables

Keep FOC-EL14 and other published relations as explicitly named closure outputs. Do not multiply a scattering amplitude or measured variance by a focusing factor by default. Implement diagnostics for stationary angular density, mean streaming velocity, finite-time centered variance and its long-time slope as different observables.

For the constant-L isotropic process, compare the mean with $v(\coth\xi-1/\xi)$ and variance with Eq. (115), not Eq. (112). Use the independent reversible angular generator and Poisson solve as the main variance reference; Monte Carlo confirms the transport pusher. Check ξ=0, ±0.1, ±1.5 and a stronger-focus case such as ±5. For small ξ use series, and for large |ξ| use scaled/logarithmic stationary weights; direct `exp`/`sinh` ratios can overflow. The stationary density may concentrate in an endpoint layer, which the reference grid must resolve.

Gate: normalization, stationary mean, variance and sign symmetry pass independently; finite-time drift is removed before taking variance; the reference grid and quadrature are converged before measuring SDE error. Eq. (115) is used only for its constant-L/constant-v/constant-ν reference; varying fields require the matching transport/reference problem. Deliver labeled diagnostic plots and a test-harness check that rejects comparison of a variance observable with a closure coefficient. Complete Release A's focusing evidence here before advertising a focused pusher.

#### PS09 — Canonical spectra and magnetostatic slab QLT

Implement a canonical one-sided total transverse PSD $P_s(k)$ in T² m on positive $k$, with $\int P_s\,dk=\delta B_\perp^2$ and $g=P_s/(8\pi)$. Source adapters retain their sidedness/component conventions and explicit conversion, not an empirical prefactor chosen to fit λ∥. Parse the complete band and tail policy. A one-component observation is not automatically total transverse or slab power.

Implement selected QLT-SLAB/PL/TS02/TS03/BENDOVER forms using their declared spectra, relativistic Ω and charge/reference-axis conventions. Validate the canonical coefficient directly and through every supported source adapter. Cover resonance outside the band, exactly on knots and near zero; never claim a physically meaningful zero where the spectrum is unavailable. Carry turbulence amplitude and the weak-turbulence/order assumptions with the model.

Gate: normalization integrals, pure power-law slopes/prefactors, kinematic scaling, factor-of-two/8π conversions and balanced limits pass against independent quadrature. λ∥ divergence from a magnetostatic gap is reported explicitly. Demonstrate grid/quadrature and any table interpolation error. A physically prescribed truncation and a missing measurement trigger different results. Deliver the first certified spectral model/adapter combinations with their documented domains.

#### PS10 — Dynamic turbulence and dissipation range

Implement selected admissible Lorentzian/Gaussian resonance kernels with declared decorrelation rates and units. DYN-STRAUSS17-PRINTED remains a nonexecutable record. Narrow resonances require a quadrature transformation or adaptive refinement around every relevant resonant root; simply increasing a uniform k grid is not a universal solution. Compute D and any required derivative from the same mathematical integral with compatible error control.

Sweep broadening width toward the sharp limit, and compare with the exact resonance Jacobian before using any small-$V_A/v$ approximation. Treat α=0 as a distinct distributional limit rather than substituting it into a numerically singular Lorentzian. Include electron/ion charge signs, relativistic gyrofrequency, dissipation break and spectral band coverage. Finite broadening can yield finite D(0)/λ∥ when the endpoint and spectrum conditions allow; magnetostatic λ∥ may diverge for the same high-k index.

Gate: resolved-width quadrature converges; Gaussian/Lorentzian normalization and sharp-limit comparisons pass; finite/different divergent limits are returned correctly. Derivative error is controlled independently of value error if it feeds a kernel. Declare a dimensionless admissibility predicate for each approximation; vague “much smaller” statements require a certified parameter range or remain a diagnostic evaluator. Deliver model-domain tables and independently reproduced fixtures.

#### PS11 — Wave populations, tensor/event kernels and clocks

Supply four nonnegative direction/polarization spectra or a named joint prescription satisfying Eq. (N5). Declare j relative to signed B, the source's polarization basis, harmonic/charge conventions, frame velocity, clock and wave-energy definition. Test outward versus AlongB mappings through a polarity change without guessing a polarization sign convention. Balanced $V_A\to0$ returns the same canonical variance and coefficient as PS09.

Choose one implementation route per population: elastic scattering in the appropriate nondispersive wave frame with a state/clock transform, or a consistent plasma-frame tensor. Do not add the second route's momentum diffusion to the first route's energy-changing kicks. Resolve each interacting channel before returning spectral energy/momentum transfers. Multi-population operator splitting is a finite-step approximation; sweep step and order. A full summed tensor needs an explicit channel decomposition before it can drive per-bin wave evolution.

Gate: particle invariant energy in the chosen wave frame, Eq. (N6), host/plasma mass shell, clock positivity and both signs of u pass; positive/negative charge resonances and exact-Jacobian limits reproduce their reference. Verify transfers sum to the particle change, include empty populations and rank-one tensors, and test finite and relativistic energies. Deliver channel identifiers, retained-order documentation and a certified wave-kernel matrix. Dispersive models remain separate until their energy/momentum closure is verified.

#### PS12 — Compressible, anisotropic and nonlinear theories

Implement only the selected complete formulas or full integral laws for which source, normalization and inputs are established. Keep Gaussian, step and inertial closed-form SOQLT resonance prescriptions separate. A local asymptote at positive μ or away from endpoints does not become a complete law by an absolute-value substitution. Tensor sources pass PS02a and cannot silently discard D_μp/D_pp.

Tabulate full integrals with their harmonic truncation and geometry assumptions, controlling value, derivative and PSD errors. Where a closed form is an asymptote of an available full integral, compare them over a parameter sweep and document the range achieving the application error budget. Where there is no verified full law, retain the source's stated asymptotic domain and keep whole-interval transport unsupported until a physical completion is specified. Strong turbulence does not make weak-QLT use admissible merely because the integral remains positive.

Gate: model-specific fixtures, resonance choice, truncation convergence and approximation-domain checks pass. Compare angular shape and transport observables, not λ∥ alone: a matching integral can hide a wrong D(0). Report source comparisons as within-source evidence, and orbit/MHD benchmarks as separate physical validation with their estimator window. Deliver evaluators and the finite domain of each supported pusher pairing; block unread or unresolved expressions.

#### PS13 — Measured spectra and simulation estimators

Implement frequency-to-wavenumber adapters with the sampling velocity/angle, 2π conversion, Taylor/order assumptions, sidedness and component/slab decomposition declared. Treat multi-component or low-Alfvén-Mach data only through an applicable model. Preserve measurement uncertainties separately from numerical quadrature/interpolation errors. Support missing bands explicitly.

For positive spectral nodes, implement a single signed log slope α and $P=P_j(k/k_j)^\alpha$ or a positive decay index −α consistently. Test both nodes, monotonicity and an exact interval integral, including the α=−1 logarithmic case. Literal zero nodes need a separately documented nonnegative interpolant; log zero is invalid, while an absent datum must not become zero. If the supplied spectrum stores integrated bin powers, reconstruct an interpolant whose bin integral matches them rather than renormalizing point samples silently.

Validate estimators using exact/flux transition solutions for the prescribed Markov process. Orbit estimates require a coarse-graining window longer than correlation/gyration scales and shorter than relaxation; the infinitesimal limit of a smooth orbit is zero. Match initial conditioning, lag, boundaries and normalization. Gate: node/integral/coverage errors and finite-lag biases are controlled; extrapolated tails and inferred slab fractions are named assumptions. Deliver observable/lag metadata and an inference plan using PAD shape/anisotropy and timing in addition to λ∥ when parameters are calibrated.

#### PS14 — Conservative particle–wave and shock coupling

Implement feedback only for an interaction closure with a consistent coefficient, wave convention, phase-space measure, clock and channel deposition. A prescribed phenomenological D does not determine a unique wave-growth spectrum. Choose event deposition or a compatible ensemble growth update for the same exchange, and document the choice. Backward sampling trajectories cannot be treated as physical energy-bearing particles.

Accumulate per-event weighted particle energy and parallel momentum changes and the negative wave changes in the same plasma frame; map to spectral-bin energy using the exact host measure. Define transactions across restart, migration and wave publication. Account for damping, spectral cascade boundary fluxes, escape, floors/caps and external/plasma work. Prevent an unphysical energy floor from deleting a transfer. Coupling can need its own substeps and positive conservative scheme; angular timestep accuracy alone does not control wave depletion.

Gate: a homogeneous closed-box nondispersive test conserves total energy and parallel interaction momentum, has nonnegative spectra, and matches the small-step growth limit. Sweep particle angular step, coupling step, spectral resolution and particle count independently. Inspect conservation per population/bin before combining totals. A shock test separates frame work and escaped energy from scattering exchange. D-126 blocks SGW-COUPLE-VK01; unresolved D-125/D-127 or incompatible source conventions block their affected capabilities. Deliver complete ledgers and a feedback-valid model matrix.

#### PS15 — GCR configurations and boundary reconstruction

Implement GCR presets through the same local coefficient library and the audited application equation. State plasma/observer frame, energy/rigidity range, charge, local turbulence inputs, polarity and perpendicular-diffusion model. GCR propagation near weak fields or large gyroradii still requires the physical ordering test; a mesh refinement cannot supply missing gyrotropic validity. A source using nonrelativistic Ω₀ at relativistic energies requires a domain restriction or a separately sourced corrected model.

For a backward application, freeze the outer-boundary distribution, sources, stopping surfaces, energy mapping and weights. Reconstruct physical intensities from f with its actual momentum convention and detector/solid-angle normalization. Test a uniform-field/no-flow analytic problem, the manufactured focusing problem and then an independently solved application problem. A forward GCR population uses forward number weights rather than the backward reconstruction convention.

Gate: known boundary values, positivity, normalization and forward/backward manufactured reconstruction pass at multiple resolutions. Compare the observable identified by each source; a ρ− sampling histogram is not tested against a forward PAD. D-125 prevents use of the unresolved van-den-Berg weight. An observational comparison is claimed only with actual data, uncertainty and calibration/validation separation. Deliver runnable GCR examples only for certified feature combinations.

#### PS16 — Regression campaign, performance and handoff

Run the baseline and full matrix for the advertised release, including refusal paths. Freeze tolerances/seeds before evaluating results. Sweep Δt, angular grid, spatial/segment resolution, coupling step and particle count only where they affect the selected model. Capture errors at common physical output times. Use independent references and separate discretization, sampling, finite-time/transient and approximation errors.

Profile representative homogeneous and varying-background workloads before optimizing: coefficient integral/table cost, subcycles per particle, memory per snapshot/cache and deposition/RNG overhead. Validate the physics again after any approximation or interpolation change. Check thread/MPI safety, deterministic conserved totals, configured reproducibility and bounded cache growth. Do not impose an invented speedup target before measuring the host baseline.

Gate: no failed required test, no changed threshold or seed used to obtain a pass, and no blocked optional path counted as validated. Deliver the release matrix, code/data/specification revisions, reproducible commands, test artifacts, open limitations and methods text. Include one independent SEP transport comparison and one GCR comparison only when those applications are advertised. If observational data were not used, describe the result as numerical/physical-model verification rather than observational validation.

### 26.4 Correspondence with the previous roadmap

PS00–PS16 and PS06b retain their earlier meanings. PS02a makes the previously scalar-only tensor contract explicit; PS04a separates angular integration from physical-time, trajectory and stochastic event handling. The earlier external-review sequence R00–R13 still maps as: R00→PS00; R01→PS01/PS02/PS07; R02→PS03; R03→PS04; R04→PS05; R05→PS06b/PS04a; R06→PS07/PS08/Release A/PS16; R07→PS09; R08→PS10/PS12; R09→PS11/PS02a; R10→PS14; R11→PS13; R12→PS15; R13→PS16. The added work packages are implementation gates, not new published model IDs or additions to the unavailable canonical fixture archive.

### 26.5 Acceptance tests, references and initial error budgets

The numbers below are proposed **engineering acceptance budgets**, not published physical constants or observational tolerances. They provide an executable starting contract for Release A. A demanding application may choose stricter budgets before testing. Do not loosen a budget after observing a failure. Use dimensionless reference scales near zeros rather than relative error alone.

| Quantity | Initial acceptance budget | Reference/uncertainty requirement |
|---|---|---|
| Algebraic conversions, boosts and closed forms | Scaled absolute 10⁻¹⁴ plus relative 10⁻¹², with stable low-energy forms | Independently evaluated identity; tolerance scaled to kinematic dynamic range |
| Smooth quadrature and normalization round trips | Relative 10⁻⁸ plus scaled absolute 10⁻¹² | Repeat with tighter quadrature; numerical reference error ≤ one quarter of budget |
| Generator number/detailed-balance residual | 10⁻¹² relative to generator-rate and mass scales | Deterministic residual; zero-flux or declared boundary fluxes accounted |
| Mean Legendre moments, modes 1–4 | Absolute 5×10⁻³ | Confidence-bounded error and independent reference-grid error included |
| Probability in $|\mu|<0.01$ | Absolute 5×10⁻⁴ | Binomial or batch-based uncertainty; correct weighting and bin integration |
| Finite-time displacement variance | Relative 1%, with a declared ballistic/near-zero absolute scale | Eq. (175) or independent joint/reference process; time, grid and sampling errors separated |
| Focusing mean and variance | Mean: absolute 5×10⁻³; variance: relative 1% | Converged burn-in/time window and Eq. (115), not the closure |
| Frame-elastic particle energy | Scaled roundoff-level residual | Stable kinetic energy; scalar/vector frame domain respected |
| Event energy/momentum deposition | Scaled roundoff residual before global reduction | Bin sum equals particle exchange; global ledger includes summation and explicit sinks |
| Coupled growth/application observables | Budget declared for the selected problem before running | Independent convergence/reference; no universal observational error assigned here |

For an observable Y, choose its budget $\epsilon_Y$ and family-wise confidence level before the run. Require $|\widehat Y-Y_{ref}|+z\,SE+\epsilon_{ref}\le\epsilon_Y$, where $\epsilon_{ref}$ bounds reference error and z accounts for the planned number of comparisons. A four-standard-error rule with an arbitrarily noisy sample is insufficient: design the particle/batch count so $zSE$ uses at most one quarter of the budget. Estimate uncertainty with independent particle batches or valid correlation-aware blocks; do not treat successive samples of the same trajectory as independent. A fixed seed does not replace uncertainty estimation. Store pilot design separately and freeze production seeds/counts before the acceptance campaign.

Use at least three timestep and, where applicable, angular/spatial grid levels. Compare the finest two plus the independent continuum reference; formal order is asserted only where the scheme's regularity assumptions hold. In singular cases measure convergence directly. For shared Brownian refinements, use a consistent bridge path; do not reuse unrelated increments and call it a strong convergence test. A stationary histogram alone can miss incorrect relaxation and streaming.

| Case ID | Setup and measurements | Pass condition and stages |
|---|---|---|
| AP-ISO-01 | ν=v=1; starts μ₀=0, 0.5, ±1 and isotropic; common times νt=0.1,1,5; modes 1–4 | Eq. (174)/(177) and exact initial data; fixed/adaptive timestep refinement; PS04/PS06 |
| AP-ISO-02 | Isotropic stationary start; νt=0.1,1,10,50; record mean, centered displacement variance and support | Eq. (175) including its finite-time correction; no endpoint-only streaming; PS04a |
| AP-ZERO-01 | q=5/3, H=0 and a supplied H>0; zero/near-zero and uniform starts | Connected flux extension, occupation and hemisphere exchange versus converged reference; PS05 |
| AP-ZERO-02 | ε=0 and a specified positive ε; pure scattering, then separately add focusing | Correct pure-scattering disconnection and focused transport behavior; drift jump handled; PS05 |
| AP-FV-01 | Nonuniform angular grids, then at least three refinements | Q row sums, stationary mass proportional to widths, positive semigroup and converged moments; PS04 |
| AP-TIME-01 | Piecewise rate with a known update inside a host interval; e.g. ν=0.4 for 0.3 s then 2.1 for 0.7 s | Exact optical time 1.59 and Eq. (177); no stale or duplicate kick; PS04a/PS06b |
| AP-RNG-01 | Split drawn Wiener increments, force rejected/refined trials and restart with a pending remainder | Eq. (N4), correct disjoint variances/covariance, retained parent sum and repeatable interval consumption; PS04a/PS06b |
| AP-FOC-01 | Isotropic scattering, ξ=0,±0.1,±1.5,±5; stationary angular and displacement diagnostics | Forward ρ+, mean, Eq. (115), sign symmetry and converged transient removal; PS08 |
| AP-FOC-02 | Finite s-domain, exact boundary/initial data for $f_0=e^{-s/L}e^{\xi\mu}$ | Physical solution reconstructed forward/backward; backward histogram follows ρ−; PS07/PS15 |
| AP-TENSOR-01 | Nonconstant PSD manufactured tensors; regular, rank-one, zero-cross and indefinite cases | Eq. (N2)/(N3), correlated covariance, log-p correction, invalid-input refusal; PS02a |
| AP-FRAME-01 | ±u, low/relativistic p, endpoint/interior μ; exact parallel and permitted vector boosts | Mass shell, round trip, positive declared clock, wave-frame elastic invariant and Eq. (N6); PS07/PS11 |
| AP-CLOCK-01 | Pure scalar angular time change with $g=1+0.3\mu>0$, then selected physical frame transform | Eq. (N7), no spurious $g'D$ drift; mathematical host-clock density proportional to $1/g$; PS04a/PS07 |
| AP-SPEC-01 | Positive pure-power-law nodes including slope −1; declared zero nodes and missing band | Both nodes, monotonicity, analytic interval variance, explicit zero/missing distinction; PS09/PS13 |
| AP-WAVE-01 | Four spectra with equal helicity marginals but different χ, plus empty/imbalanced populations | Weights satisfy Eq. (N5), resonance selection differs correctly, balanced limit and channel sums pass; PS11 |
| AP-GEOM-01 | Uniform B/ν on equivalent field-line and 3-D meshes with nonuniform resolution | Same coefficients and common-time ensemble observables; physical optical time unchanged; PS06b |
| AP-GEOM-02 | Radial $B=B_0(r_0/r)^2\hat r$ on r>0, with $1/L=2/r$ and $A\propto r^2$ | Divergence-free background away from origin; correct tube/volume measure and converged matching distributions; PS06b |
| AP-GEOM-03 | Split/merge segments, refine/coarsen cells and force migration/repartitioning | Number, normalized sampling, time and statistics invariant within declared discretization budgets; PS06b |
| AP-GEOM-04 | Checkpoint mid-subcycle with pending noise/exchanges; restart and vary rank count | No lost/duplicated time, kick or deposition; preserved frame/weight/epoch/RNG; promised reproducibility or statistical equivalence; PS06b |
| AP-ENERGY-01 | Frozen homogeneous plasma-frame scattering, then separate single-wave elastic test | Energy invariant in the configured elastic frame; plasma exchange/bin ledger closes; PS11 |
| AP-ENERGY-02 | Closed particle–wave box and topology change; separate explicit damping/escape run | Energy/momentum, positive spectrum and growth limit; external terms in ledger, no double counting; PS14 |
| AP-REFUSE-01 | Wrong capability/λ tag, unknown frame, B=0, inadmissible source, domain exit, missing power, zero step | Typed failure names parameter/model/source/register; no successful value or silent fallback; PS01/PS16 |

References for nonisotropic moments and spatial spread are assembled for the same coefficient, extension, initial/boundary problem and sample measure. The isotropic analytical exponent is not assigned to a q-form or ε-form. A 1-D/3-D comparison uses the same approximation order, frame and physical background, rather than two differently defined implementations.

### 26.6 Minimum reproducible configurations

**A — initial prescribed-scattering release.** PA-ISO, ν=v=1 in normalized verification units, no focusing/flow/momentum change, source-defined number measure. Compare compatible EM, Milstein and bounded rotation with fixed νΔt=0.01,0.003,0.001,0.0003 as initial study values. Output Legendre modes, endpoint entrance, finite-time displacement and bias/confidence bounds. The accepted production step is chosen from convergence, not assumed from these initial values. SI application examples use declared physical species/energies and rates.

**B — nonisotropic pair release.** PA-QFORM at q=5/3 with H=0 and a separately supplied positive H, or PA-EPS with ε=0 and a supplied ε>0. Declare the internal extension, source normalization and kernel explicitly. Compare against angular grids 200,400,800 or tighter until reference error is below budget. Release no derivative-based q-form path merely because these coefficient grids can be evaluated.

**C — focused transport.** v=ν=1 and ξ=±1.5: stationary mean ±0.438124726316, centered variance slope 0.178275780337345; Eq. (112)'s closure value is 0.292083150877230 and is a different output. Add ξ=0 and ±5. On a finite physical reconstruction domain use exact manufactured boundaries and compare f, not an unweighted histogram.

**D — spectral coefficient.** Canonical slab PSD with declared total variance, break/cutoff indices, band/tails and amplitude. Output integrated variance, resonant population/band coverage, D(μ), regularity, λ∥ status and certified weak-turbulence/domain predicates. A gap-filling prescription is a named physical/model choice rather than an implementation detail.

**E — wave-frame energy change.** One nondispersive population with signed u and known p′, elastic angular kick, full back-transform and identified spectral deposition. Output particle/plasma/wave energy and parallel momentum, both elapsed clocks and event/bin IDs. Extend to four populations only after the channel and split contract is specified.

**F — host equivalence and lifecycle.** Use identical immutable uniform physics in cell and field-line adapters, then the radial background of AP-GEOM-02. Force rate updates, AMR/segment changes, migrations and a mid-subcycle restart. Compare common-time observables and integrated ledgers. This is the first AMPS integration acceptance case, not a re-run of local coefficient fixtures.

### 26.7 Error handling and physical alternatives

| Failure | Required behavior | Permitted continuation |
|---|---|---|
| Unread/ambiguous source, incompatible dimensions or signs | Block affected model; name D-/U-item | Use another explicitly selected, independently established model |
| Local approximation/domain exit | No coefficient value from that evaluator | Configured physical composite/fallback with provenance, or stop the affected path |
| Missing measured spectral coverage | `MissingCoverage`, preserve missing status | Supplied tail/other data with its named assumptions |
| Physically truncated spectrum | Evaluate actual declared power and report zero/barrier/divergent λ if applicable | Continue with the specified transport operator, not an artificial scattering floor |
| Unbounded/discontinuous drift incompatible with kernel | `IncompatibleKernel` | A certified derivative-free or correctly transformed/regularized model pairing |
| Nonconverged resonance/quadrature/table | `Nonconverged` with error/work limit | Refine within the declared domain or use a verified evaluator |
| B=0 or violated gyrotropic ordering | Physical diagnostic, no B floor | A justified host transport/orbit model that covers that regime |
| Zero elapsed step / invalid state | `Nonprogress`/`InvalidInput`; no successful kick | Correct the step/state or stop with the recorded cause |
| Wave depletion/floor violation | Preserve transfer ledger, do not erase energy | Conservative coupling refinement or explicit physical reservoir |
| Data/specification mismatch | `VersionMismatch` for the affected registry path | Supply matching archive or an explicitly scoped independent-data implementation |

### 26.8 Completion and handoff

Deliver a final table with each stage's status, code/data/specification revision, command, measured error and uncertainty, open item and next action. Deliver the supported model–kernel–geometry–frame–domain matrix and refuse all combinations outside it. Keep source completeness, physical admissibility, kernel verification, AMPS integration and observational evidence as distinct columns.

Release A is complete when its selected prescribed models, time/geometry/measure adapters and FTE reconstruction gates pass; it need not advertise wave evolution or every literature record. Release B/C and application validation each require their additional evidence. No “all models ready” statement follows from the 172 inventory IDs or the 652 intended fixtures. A fully described model still requires an implemented and tested pusher on the actual AMPS checkout.

The handoff package contains the API and configuration reference; source/derived equation mapping; immutable data manifest; reproduction commands; tests, results and plots; performance measurements; refusal examples; open scientific assumptions; and methods text populated from the run record. The implementing agent should work in small stage-complete changes, retain the baseline and avoid tuning thresholds or source formulas to obtain a pass.

**Present review boundary.** This version makes the scientific execution contract and roadmap concrete and supplies independent mathematical checks and a compilable declaration. The actual AMPS implementation, matching v1.1 companion archive, delegated companion texts and end-to-end numerical/observational validation remain evidence that must be supplied or produced during these stages. They are not represented as completed by this specification revision.
