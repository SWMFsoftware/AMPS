# Mean-free-path models for SEP and GCR transport in the heliosphere

**File:** MEAN_FREE_PATH_MODEL.md  
**Specification version:** 1.2  
**Date:** 9 October 2026 (version 1.0: 8 October 2026; changes in Section 20.1)  
**Literature cutoff:** 7 October 2026; principal review interval 2006–2026, with pre-2006 foundations labelled.  
**Companion documents:** PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md (revision 1.4) and the perpendicular diffusion specification.  
**Purpose:** a reference for implementing mean-free-path prescriptions in SEP and GCR transport codes, and a source for the methods sections of publications.

This document collects the mean-free-path models that heliospheric SEP and GCR transport codes actually use: observational constraints and event-fitted values, empirical presets of SEP codes, pitch-angle diffusion shapes and their exact normalization, closed-form turbulence-based λ∥, shock-region prescriptions, GCR modulation forms with their published parameter sets, and the turbulence inputs these forms need. Each formula is given as printed in its source, with the source location and the verification level of the reading (Section 1.2). Algebraic consequences computed in this review are labelled "derived"; evaluations of published formulas at published inputs are labelled as such and are not presented as published results.

Where sources disagree, contain apparent typographical errors, or could not be read, this document records the discrepancy or the gap (Section 14) and does not fill it by assumption. Several items require a decision by the user or by the host code (Section 14.2); in particular, no λ formula was found for the target-code mode "Tenishev2005AIAA", and the mode "Chen2024AA" most plausibly refers to Chen et al. (2024, ApJ 965, 61), which is an inference. Section 21 states the implementation contracts that follow from the equations without adding physics (runtime states, units, parsing, domains, stable numerical evaluation and validation layers), and Section 22 orders the implementation work.

The companion bundle MEAN_FREE_PATH_MODEL_DATA.zip contains every transcribed table and parameter set with provenance, 34 groups of verification fixtures, a self-checking script and a separate set of independent mathematical checks (Section 19). MEAN_FREE_PATH_MODEL_additions.bib contains BibTeX entries only for references that are not in the supplied bibliography (Section 18).

## Contents

1. [Scope, companion specifications, and conventions](#1-scope-companion-specifications-and-conventions)
2. [Variables, units, and kinematics](#2-variables-units-and-kinematics)
3. [Model inventory](#3-model-inventory)
4. [Observational constraints at and inside 1 au](#4-observational-constraints-at-and-inside-1-au)
5. [Empirical SEP prescriptions used by transport codes](#5-empirical-sep-prescriptions-used-by-transport-codes)
6. [Pitch-angle diffusion shapes and the λ∥ normalization](#6-pitch-angle-diffusion-shapes-and-the-λ-normalization)
7. [Closed-form turbulence-based parallel mean free paths](#7-closed-form-turbulence-based-parallel-mean-free-paths)
8. [Shock-region and self-generated-wave prescriptions](#8-shock-region-and-self-generated-wave-prescriptions)
9. [GCR modulation prescriptions and published parameter sets](#9-gcr-modulation-prescriptions-and-published-parameter-sets)
10. [Coupling to perpendicular and radial transport](#10-coupling-to-perpendicular-and-radial-transport)
11. [Turbulence inputs](#11-turbulence-inputs)
12. [Library interface additions](#12-library-interface-additions)
13. [Verification fixtures](#13-verification-fixtures)
14. [Discrepancy, decision and gap register](#14-discrepancy-decision-and-gap-register)
15. [Model selection](#15-model-selection)
16. [Publication methods text](#16-publication-methods-text)
17. [References](#17-references)
18. [Bibliography comparison](#18-bibliography-comparison)
19. [Companion data](#19-companion-data)
20. [Verification and revision record](#20-verification-and-revision-record)
21. [Implementation contracts, numerical safeguards and validation design](#21-implementation-contracts-numerical-safeguards-and-validation-design)
22. [Codex implementation roadmap (final section)](#22-codex-implementation-roadmap-final-section)

## 1. Scope, companion specifications, and conventions

### 1.1 Scope

This specification describes the mean-free-path prescriptions that SEP and GCR transport codes have used in the heliosphere, mainly in work published between 2006 and 2026. It covers:

- the observational constraints that event analyses report (Section 4);
- empirical λ presets of SEP transport and acceleration codes (Section 5);
- pitch-angle diffusion shapes and their exact normalization to λ∥ (Section 6);
- closed-form turbulence-based λ∥ (Section 7);
- shock-region prescriptions (Section 8);
- GCR modulation prescriptions with their published parameter sets (Section 9);
- the turbulence inputs these forms require (Section 11).

Pre-2006 papers appear only where a recent model is defined through them, and are labelled as foundations.

The document is self-contained as a catalogue and review of mean-free-path models, but not as an implementation of every closure: it deliberately delegates several equations to two companion specifications and does not repeat them:

- **PARALLEL** — *Parallel diffusion coefficient models for SEP and GCR transport*, PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md, revision 1.4 (8 October 2026). It specifies the Parker transport context, kinematics and units (its Sections 1–2), the generic rigidity power law (Section 5), the smooth broken rigidity law (Section 6), pitch-angle normalization (Section 7), slab QLT and spectrum conventions (Section 8), nonlinear closures including NLGCE-F (Sections 9–11), turbulence adapters (Section 12) and the Bohm model (Section 13). References here to "PARALLEL rev. 1.4, Section n" point to that document.
- **PERPENDICULAR** — the companion perpendicular diffusion specification. It is cited here only by topic (FLRW, NLGC/UNLT, constant ratios, Dröge-type gyroradius scaling), not by section number.

An implementation pins the exact revisions of both companion documents it uses and checks every section and equation number cited here against those revisions before it implements a delegated equation; a delegated path whose companion text is not available to the implementer is blocked rather than reconstructed (Section 21.10).

### 1.2 Status labels

| Label | Meaning |
|---|---|
| Full text ×n | Read in the full text; n independent readings (fetches, possibly of the same copy) agreed, counted over versions 1.0 and 1.2 (n = 1 means not double-checked) |
| Abstract only | Only the abstract or metadata was accessible |
| Secondary | Taken from a named citing paper |
| As printed | Transcribed exactly, including apparent typographical errors; not corrected |
| Not found | Not located in any accessible source; nothing is substituted |
| Derived | An algebraic or numerical consequence computed in this review, not a published statement |
| Evaluation of printed formula | A published formula evaluated at published inputs; not a value the authors published |
| **Implementation requirement** | A design decision for the library, not a physical statement |

Every number in the companion data carries one of these labels (Section 19).

### 1.3 Mean-free-path definitions

Four different quantities are called "mean free path" in this literature.

$$
\lambda_\parallel=\frac{3\kappa_\parallel}{v},\qquad
\lambda_\perp=\frac{3\kappa_\perp}{v},
\tag{1}
$$

$$
\lambda_r\equiv\lambda_\parallel\cos^2\psi\quad\text{(SEP focused-transport convention)},
\tag{2}
$$

$$
\lambda_{rr}\equiv\frac{3\kappa_{rr}}{v}=\lambda_\parallel\cos^2\psi+\lambda_\perp\sin^2\psi\quad\text{(tensor projection)},
\tag{3}
$$

where ψ is the angle between the mean field and the radial direction. Equation (2) is used by Agueda et al. (2010, Eq. 2), He et al. (2011), Wijsen et al. (2019, Eq. 11), Strauss et al. (2017, Eq. 10), Kubo et al. (2015, Eq. 8) and Zhang et al. (2023); Eq. (3) is used by Chhiber et al. (2017, Eq. 2), Zhao et al. (2018, Eq. 15) and the GCR tensor (Eq. 32). Strauss et al. (2017, GLE) and Kubo et al. (2015) write the SEP-convention quantity of Eq. (2) as λ_rr; the symbol therefore does not identify the definition. The two definitions coincide only when λ⊥ sin²ψ is negligible. The fourth quantity is the isotropic-scattering length λ of full-orbit codes (Section 6.6), for which κ = vλ/3.

A spatially constant λ_r implies a λ∥ that grows with r through 1/cos²ψ; a constant λ∥ implies a λ_r that decreases. These are different physical prescriptions (PARALLEL rev. 1.4, Section 1.2 makes the same point for κ).

**Implementation requirement.** Every λ value entering or leaving the library carries a tag from {λ∥, λ_r, λ_rr, λ_iso, λ_unspecified}. "λ_unspecified" is used where a source does not say (SEP-PATH09 and SEP-SOLPENCO05); converting it requires a user decision (Section 14, U-11).

## 2. Variables, units, and kinematics

Kinematics, units and the exact unit conversions follow PARALLEL rev. 1.4, Section 2 (its Eqs. 7–9): pc = [T(T + 2mc²)]^{1/2}, v = pc²/(T + mc²), rigidity 𝓡 = pc/|q| in volts, r_L = 𝓡/(cB₀), 1 AU = 149 597 870 700 m. Two additional points recur in this review.

**Momentum versus rigidity normalization.** Several presets normalize by pc = 1 GeV while their text says "1 GV" (Kozarev et al. 2013 state both). Since

$$
\frac{pc}{1\ \mathrm{GeV}}=|Z|\,\frac{\mathcal R}{1\ \mathrm{GV}},
\tag{4}
$$

the two coincide numerically for every singly charged species (|Z| = 1: protons, electrons, positrons, singly charged ions), whatever its mass; for an alpha particle pc = 1 GeV is 𝓡 = 0.5 GV (fixture F-KIN-01). At fixed kinetic energy, by contrast, the rigidity depends on the mass. Zhang et al. (2023) give momentum "in GV", i.e. momentum per unit charge. **Implementation requirement:** each preset stores whether its variable is pc, 𝓡, or kinetic energy (total or per nucleon), and the library never substitutes one for another.

**Constants used in the fixtures.** Exact c, e and AU; CODATA 2018 masses m_p = 1.67262192369 × 10⁻²⁷ kg, m_e = 9.1093837015 × 10⁻³¹ kg and m_α = 6.6446573357 × 10⁻²⁷ kg (as in PARALLEL rev. 1.4); and, where a paper gives a length in solar radii, the IAU 2015 nominal R☉ = 6.957 × 10⁸ m. The papers do not state which solar radius they used, so conversions from R☉ are labelled derived. A 1 MeV proton has 𝓡 = 43.3306 MV, and a proton with pc = 1 GeV has T = 432.988 MeV.

## 3. Model inventory

The inventory lists what the literature provides. Whether a given configured run can execute is a separate question, answered by the runtime state of Section 21.1, which is evaluated per configured run; the last column gives the state when the inputs named in the source are supplied.

| ID | Section | Quantity returned | Implementable from this document? | Runtime state (Section 21.1) |
|---|---|---|---|---|
| OBS-* (event data) | 4 | λ∥ or λ_r benchmarks | Data only | REFERENCE_DATA_ONLY |
| SEP-PATH09, SEP-EPREM10, SEP-EPREM13, SEP-MFLAMPA19, SEP-MFLAMPA25, SEP-SOFIE24, SEP-ZHANG23, SEP-PARASOL25, SEP-SOLPENCO05, SEP-SPARX15, SEP-MARSH13, SEP-HE11, SEP-WANGQIN15, SEP-KUBO15, SEP-PARADISE19, SEP-STRAUSS17G, SEP-STRAUSS15, SEP-DROGE16P | 5.2 | λ (tag per preset) | Yes, with the tag and normalization recorded; λ_unspecified presets need U-11 | READY_EXPLICIT_INPUTS, except SEP-PATH09 and SEP-SOLPENCO05 (REQUIRES_USER_DECISION, U-11) and SEP-EPREM10 away from 1 AU (REQUIRES_SOURCE_OR_CODE_AUDIT: radial factor not legible, Section 5.2) |
| SEP-LAITINEN16, SEP-LAITINEN18 | 5.2, 11.4 | λ∥ via QLT over a turbulence model | Yes, through the QLT adapter of PARALLEL rev. 1.4 with Eq. (41) | READY_EXPLICIT_INPUTS when the pinned PARALLEL QLT adapter is available; otherwise the implementation reports BLOCKED_DEPENDENCY (Section 21.10) |
| SEP-MINOSHIMA26 | 5.2 | λ∥ = ξv/Ω_n (Ω_n undefined in the source) | Yes, with ξ user-supplied and the gyrofrequency stated (D-29) | REQUIRES_SOURCE_OR_CODE_AUDIT (D-29) |
| SEP-CHEN24 | 4.5, 5.2 | κ∥ (published); λ∥ (derived) | κ∥ yes; λ∥ = 3κ∥/v once the species is chosen (U-2) | READY_PUBLISHED_KAPPA_ONLY |
| Tenishev2005AIAA | 5.4 | — | No: no formula found (U-4) | REQUIRES_SOURCE_OR_CODE_AUDIT (U-4) |
| PA-QFORM, PA-EPS, PA-DROGE-VA, PA-ISO, PA-KOLMO (M-FLAMPA), PA-AMPS-I…V | 6 | D_μμ(μ) and its amplitude from λ∥ | Yes; PA-QFORM-LANG-PRINTED and PA-EPS-PACHECO-PRINTED only as explicit variants | READY_EXPLICIT_INPUTS; the printed variants and PA-EPREM are READY_PUBLISHED_VARIANT (D-2, D-3, D-4) |
| MFP-FOCUS-HW13 | 6.9 | focusing-corrected λ∥ | Yes, explicit option | READY_EXPLICIT_INPUTS |
| QLT-TS03-P | 7.2 | λ∥ (ions) | Yes | READY_EXPLICIT_INPUTS with named turbulence inputs (U-9); otherwise REQUIRES_EXTERNAL_INPUT |
| QLT-TS03-E-RS | 7.2 | λ∥ (electrons) | Yes | as QLT-TS03-P |
| QLT-TS03-E-DT | 7.2 | λ∥ (electrons) | Only with a variant flag (U-7) | READY_PUBLISHED_VARIANT (D-1) |
| QLT-ZANK98 | 7.5 | λ∥ | Yes, with the variance convention stated (U-12) | READY_EXPLICIT_INPUTS with the variance convention named; otherwise REQUIRES_USER_DECISION (U-12) |
| SHOCK-BOHM, SHOCK-AFANASIEV15, SHOCK-PARASOL, SHOCK-MFLAMPA | 8 | λ or κ near shocks | Bohm, Afanasiev and M-FLAMPA yes; PARASOL needs U-6 | Bohm, Afanasiev, M-FLAMPA: READY_EXPLICIT_INPUTS with host shock coordinates (Section 21.7); PARASOL: REQUIRES_SOURCE_OR_CODE_AUDIT (U-6) |
| GCR-NWU14 and variants, GCR-CORTI19, GCR-LUO19, GCR-HELMOD17, GCR-HELMOD19, GCR-BOBIK12, GCR-STRAUSS11, GCR-EFFENBERGER12, GCR-WANG19, GCR-TOMASSETTI17, GCR-PERUGIA21, GCR-PERUGIA25, GCR-JIANG23, GCR-DUAN25, GCR-QINSHEN17, GCR-EB13 | 9 | K∥ (or λ∥) | Yes, with normalization recorded; GCR-LUO19 form unknown (data only) | READY_EXPLICIT_INPUTS with a named published parameter set and, where the field B(r) is built from Eq. (42) with a published 1 AU value, the U-13 normalization stated (HelMod also with its activity index and the D-21 numeric convention; GCR-QINSHEN17 and GCR-EB13 with their turbulence inputs, GCR-QINSHEN17 through the pinned PARALLEL NLGCE-F, otherwise BLOCKED_DEPENDENCY); NWU variants with the D-9 printed exponents REQUIRES_USER_DECISION (U-10); GCR-LUO19 REFERENCE_DATA_ONLY; GCR-BOBIK12 REQUIRES_SOURCE_OR_CODE_AUDIT (units not stated); GCR-JIANG23 REQUIRES_USER_DECISION (K₀ unit printed as cm⁻² s⁻¹ sr⁻¹ GeV⁻¹); GCR-WANG19, GCR-PERUGIA21, GCR-PERUGIA25 and GCR-DUAN25 REQUIRES_EXTERNAL_INPUT (time-dependent values printed only in figures), GCR-DUAN25 then READY_PUBLISHED_VARIANT (D-33) |
| TURB-* | 11 | turbulence inputs | Adapters and published values (data only) | REFERENCE_DATA_ONLY |

Every model ID maps to a source location, and to its runtime state, in `parameters/model_registry.json` (Section 19); the data groups OBS-* and TURB-* are the files `observations/event_fitted_mfp.csv` and `turbulence/published_turbulence_values.csv`, whose rows carry their own source, location and verification. The host-code label "Chen2024AA" is not a model ID; it is mapped by U-3 (Section 5.4).

## 4. Observational constraints at and inside 1 au

This section records what event analyses actually report. It is a constraint and benchmark source, not a model. Every value is tagged with the quantity it refers to (λ∥ or λ_r), because the two are mixed in the literature (Section 1.3). Complete transcriptions are in the companion file `observations/event_fitted_mfp.csv` (Section 19).

### 4.1 The Palmer consensus

The primary statement, from the abstract of Palmer (1982) [M01], is:

> "A consensus is found: at 1 AU, λ∥ = 0.08–0.3 AU over a wide range of rigidity, R = 5 × 10⁻⁴ to 5 GV."

The printed range 5 × 10⁻⁴–5 GV is identical to the 0.5–5000 MV quoted by later authors. The same abstract states that scatter-free events are those "where λ∥ ≳ 1 AU", that "K⊥r / K∥ < 0.1 at 1 AU", and that "a reasonable mean is K⊥r/β = 10^21 cm² s⁻¹". The quantity is the **parallel** mean free path at **1 AU**; the word "radial" does not appear in the abstract. Status: abstract only (the full text was not accessible).

Later restatements agree on the parallel quantity and the 0.5–5000 MV range: Tautz and Shalchi (2013) [M02], Chhiber et al. (2017) [M03], Reames (2013) [M04], Lavasa et al. (2026) [M05], Minoshima et al. (2026) [M06], and Subashchandar et al. (2025) [M07]. Subashchandar et al. argue that the band should not be extrapolated to the near-Sun region. Shalchi et al. (2006) [M08] convert Palmer's mean perpendicular coefficient to λ⊥ ≈ 0.0067 AU and give 0.02 ≤ λ⊥/λ∥ ≤ 0.083.

Bieber et al. (1994) [M09] state in their abstract that "the mean free path of cosmic-ray electrons and protons may be fundamentally different at low to intermediate (less than 50 MV/c) rigidities", that "the mean free path of 1.4 MV/c electrons is often similar to that of 187 MV/c protons", and that "'consensus' ideas about cosmic-ray mean free paths may require drastic revision". Whether their tabulated values are radial or parallel mean free paths was not determined (full text not accessible). Engelbrecht et al. (2022) [M10] summarize Bieber et al. (1994) as finding the consensus "applicable to electron parallel MFPs at rigidities below 25 MV", and report that parallel MFPs "span a range of about two orders of magnitude … with values often larger than the Palmer (1982) consensus range". Chen et al. (2024) [M11] use the phrase "Palmer consensus" for a different statement: that QLT predicts κ∥ about an order of magnitude below observations.

**Specification statement.** A consensus check is λ∥(1 au) ∈ [0.08, 0.3] au for 0.5 MV ≤ P ≤ 5 GV. It is a historical band, not a bound. Lang et al. (2024) report that "Proton pMFPs reported on here remain mostly well above the Palmer consensus range" [M12].

### 4.2 Definitions used by the event-fitting papers

| Group | Fitted quantity | Relation to λ∥ | Pitch-angle model |
|---|---|---|---|
| Agueda et al. 2008–2014; Agueda & Lario 2016; Pacheco et al. 2019 [M13, M14, M15, M16] | λ_r, constant along the field line and in energy within a fit | λ_r = λ∥ cos²ψ (Agueda et al. 2010, Eq. 2) | ε-form or q-form (Section 6.3) |
| Lang et al. 2024; Lavasa et al. 2026 [M12, M05] | λ∥ at the observer | λ∥ = 3κ∥/v; λ_r = λ∥ cos²ψ used only to convert literature λ_r, with ψ = 45° | q-form with H = 0.05, printed factor (1−μ)² (Section 6.2) |
| Dröge et al. 2016 [M17] | λ∥ "normalized to a distance of 1 au" | — | not read (abstract only) |
| Chen et al. 2024 [M11] | κ∥ from measured PSD (QLT), not an event fit | none given | QLT, Eq. (8) |

Converting a published λ_r into λ∥ requires the field angle ψ that the authors used. Lang et al. (2024, footnote 8) use ψ = 45°, so cos²ψ = 0.5 and λ∥ = 2λ_r. Agueda et al. and Pacheco et al. fit λ_r directly, and the per-event ψ used in their simulations is not part of the tabulated results. A converted value is therefore labelled as derived.

### 4.3 Event-fitted values

All values below are transcriptions. Formal per-fit uncertainties are not given in any of these tables; the papers report goodness of fit instead, and the grids on which λ was searched set the effective resolution.

**Helios electrons, 0.31–0.94 au (Pacheco et al. 2019, Tables 2–3) [M16].** Fifteen fits of 0.3–0.8 MeV nominal electrons give λ_r = 0.020–0.270 AU. λ_r was searched on a logarithmic grid from 0.01 to 0.5 AU. The paper states: "We find no dependence of the radial mean free path on the radial distance." Its D_μμ uses ε = 0.01 (Section 6.3). Table 4 of the same paper compares with Kallenrode et al. (1992a), Kallenrode (1993), and Agueda & Lario (2016); those values are secondary here.

**Near-relativistic electrons at 1 au (Agueda et al. 2014, Table 3) [M14].** Seven events (1999–2004), λ_r per energy channel from ACE/EPAM (62–312 keV) and Wind/3DP (50–230 keV). The authors summarize the range as 0.12–0.44 AU (Sect. 5.1) and state that "the electron radial mean free path is rigidity dependent in the range from 0.3 MV to 0.5 MV", with "increasing values of λr toward smaller rigidities". That range matches the ACE "Total" column (0.12–0.44 AU); the individual channel values span 0.10–0.52 AU. The D_μμ has q = 1.66, and perpendicular diffusion is neglected.

**Single-event electron spectrum (Lang et al. 2024, Table 3) [M12].** Wind/3DP electrons, 20 January 2022: λ∥ = 0.28, 0.08, 0.37, 0.23, 0.18, 0.36 au at 0.04, 0.07, 0.11, 0.18, 0.31, 0.52 MeV, with R² = 0.97, 0.79, 0.86, 0.96, 0.96, 0.65. The paper's stated acceptance criterion is R² ≥ 0.80; two of the six printed rows are below it. Lang et al. fit 15 events between 1998 and 2022 (electrons 0.02–4.0 MeV, 0.14–4.48 MV; protons 1.30–130 MeV, 49.4–510.7 MV); all other fitted values are only plotted (their Fig. 8).

**GLE 73, 28 October 2021 (Lavasa et al. 2026, Tables 2–3) [M05].** For STEREO-A, electrons from 0.05 to 3.3 MeV (0.232–3.78 MV) give λ∥ = 0.07–0.17 AU, and protons from 14.3 to 77.5 MeV (164–389 MV) give λ∥ = 0.26–0.38 AU. R² is 96–99% for electrons and 77–95% for protons. Table 3 adds λ⊥ from multi-observer fits, with λ⊥/λ∥ = 0.9–3% for electrons and 4.5–10% for protons. Two printed energy–rigidity pairs are internally inconsistent and are transcribed as printed and flagged (Section 14):

- Table 2 prints "0.94" MeV between the 0.08 and 0.115 MeV rows, with rigidity 0.324 MV. An electron of 0.94 MeV has P = 1.358 MV; the printed 0.324 MV corresponds to 0.094 MeV.
- Table 3 prints 2.5 MeV with 2.46 MV. Table 2 pairs 2.46 MV with 2.0 MeV, and 2.5 MeV corresponds to 2.967 MV.

These are kinematic checks with P = [E(E + 2m_ec²)]^{1/2} (Section 2); they do not establish which printed number is wrong. All other rows agree with the printed rigidity to within rounding.

**Ranges quoted only in abstracts.** Dröge et al. (2016) [M17]: λ∥ = 0.15–0.6 au and λ⊥ = 0.005–0.01 au, normalized to 1 au, for August 2010 electron events. Agueda et al. (2008) [M18]: λ_r = 0.9 AU for the 2000 May 1 event. Dröge (2000) [M19]: earlier events "exhibit mean free paths in the range of 0.02 to 0.5 AU" and "a uniform shape of the functional form of the rigidity dependence, which varies only in absolute height for different events, can explain all observations".

**Spatially structured fits.** Agueda et al. (2010) [M13] fit the 2000 February 18 event with λ_r = 3.2 AU inside r = 1.2 AU and λ_r = 0.2 AU beyond it. Minoshima et al. (2026) [M06], using data assimilation for 30 March 2022 (BepiColombo at 0.6 AU, STEREO-A at 1.0 AU), report λ∥ that "decreases over time and reaches roughly 0.5-1.0 AU at STEREO-A during the decay phase", with medians of about 0.5 AU at 1.5 MeV and 1.0 AU at 5.9 MeV for t = 12–23 h. They give uncertainties only as posterior distributions.

**Model inputs, not fits.** Battarbee et al. (2018) [M20] use λ = 0.3 au for protons. Dalla et al. (2020) [M21] use constant λ = 0.1 and 0.5 AU, and 1.0 and 0.3 AU for GLE 71, stating "There is no consensus within the literature about the degree of scattering experienced by GLE energy protons." Houeibib et al. (2025) [M22] hold λ∥ constant at 0.1–1 AU.

### 4.4 Rigidity dependence

Engelbrecht et al. (2022) [M10] summarize the observations as electron parallel MFPs that are rigidity-independent at low rigidity, and proton parallel MFPs that at higher rigidities "appear to display" a "~P^{1/3} dependence, expected from magnetostatic quasilinear theory (QLT)". Lang et al. (2024) find that "For most of the proton pMFPs no clear rigidity dependence can be discerned" and that "For electrons at the lowest rigidities, pMFPs display an increase with decreasing rigidity". Agueda et al. (2014) find the same low-rigidity electron trend in λ_r. Lavasa et al. (2026) find the opposite trend for 0.20–0.45 MV electrons in one event (λ∥ rising with rigidity) and λ∥ rising with rigidity for 180–400 MV protons. Chhiber et al. (2017) quote Bieber et al. (1994) power indices of 0.2–0.56 for 10–10³ MV; that is secondary. Engelbrecht et al. (2022) report that Dröge & Kartavykh (2009) [M23] found observed electron pitch-angle distributions inconsistent with dynamical-QLT predictions (secondary). Tan et al. (2011) [M24] relate the change from scatter-free transport of lower-energy electrons to diffusive transport of higher-energy electrons (~25–500 keV) to the level and breaks of the magnetic power spectrum (abstract only; no λ values in the abstract).

The theoretical regimes listed by Lang et al. (2024, Sect. 4) are ~P² at high rigidity, ~P^{1/3} at intermediate rigidity (magnetostatic QLT), ~P^{2−s} in the inertial range, and ~P^{2−p} for random-sweeping low-rigidity electrons, where s and p are the inertial and dissipation-range indices (Section 7.2). The low-rigidity electron behaviour is attributed to resonance with the slab dissipation range. Strauss et al. (2020) [M25] use (their equation is numbered (9) in the arXiv PDF and (13) in the ar5iv rendering) the dissipation-range onset

$$
k_d=\frac{2\pi}{V_{sw}}\left(a+b\,\Omega_i\right),\qquad a=0.2\ \mathrm{Hz},\quad b=1.76,
\tag{5}
$$

after Leamon et al. (2000), and note it may vary by a factor of about 5. Ω_i is the ion cyclotron frequency in the source's notation. Strauss et al. define Ω_i as "the proton cyclotron frequency" without stating whether it is angular or cyclic; the Leamon et al. (2000) fit with the same a and b is written with Ω_ci/2π (Section 11.6).

### 4.5 Radial dependence

**Chen et al. (2024) PSP fit [M11].** From PSP orbits 5–13 (r = 0.062–0.8 AU), Chen, Giacalone, Guo and Klein compute κ∥ from QLT with measured spectra and fit

$$
\kappa_\parallel=(5.16\pm1.22)\times10^{18}\;
r^{\,1.17\pm0.08}\,E^{\,0.71\pm0.02}\ \ \mathrm{cm^2\,s^{-1}},
\tag{6}
$$

"where r and E are in the units of AU and keV respectively". The stated domain is "within about 0.1-0.8AU for the energetic particles of the energy between 100keV-1GeV". The resonance calculation is for protons; the species is not stated for Eq. (6) itself. The paper gives no λ∥ and no κ = vλ/3 relation. The underlying definitions are

$$
\kappa_\parallel=\frac{v^2}{4}\int_{\mu_{\min}}^{1}\frac{(1-\mu^2)^2}{D_{\mu\mu}}\,d\mu,
\qquad \mu_{\min}=0.05,
\tag{7}
$$

$$
D_{\mu\mu}=\frac{\pi}{4}\,\Omega_0\,(1-\mu^2)\,\frac{f_{\rm res}P(f_{\rm res})}{B_0^2},
\qquad f_{\rm res}=\frac{k_{\rm res}V_{SW}}{2\pi},\quad k_{\rm res}=\left|\frac{\Omega_0}{v\mu}\right|,
\tag{8}
$$

with Ω₀ = qB₀/mc and P(f) the measured spectrum of B_N. For a D_μμ even in μ, Eq. (7) equals the symmetric form (v²/8)∫₋₁¹ of Eq. (11), restricted to |μ| ≥ 0.05. The abstract's phrase "increases exponentially" describes a power-law fit.

**Other radial evidence.**

- Subashchandar et al. (2025) [M07] find, for PSP at about 0.06–0.3 AU with SOQLT, that κ∥ for 500 keV protons "follows approximately r^{1.1}" and for 1 GeV protons "increases more steeply as r^{1.5}".
- Pacheco et al. (2019) find no dependence of λ_r on r for Helios events (Section 4.3).
- Dröge et al. (2016) state that MFPs "can vary significantly not only as a function of radial distance, but also of heliospheric longitude".
- Minoshima et al. (2026) take λ∥ ∝ |B|⁻¹ ∝ r^b with b = 2 for r ≪ A and b = 1 for r ≫ A (A = 1 AU), quote Palmer (1982) as giving b = 0.9–1.8 over 1.0–5.0 AU, and quote b = 1.17 over 0.1–0.8 AU from Chen et al. (2024). The Palmer value is secondary.
- Zhong, Wang and Qin (2024) [M26] model λ_r as a power function of r for 13–64 MeV protons; the exponent was not available (abstract only).
- Cao, Wang and Guo (2025) [M27] use κ = κ₀(E/E₀)^β with E₀ = 1 MeV inside a peak-flux model; this is not a λ fit.

### 4.6 Use in a library

**Implementation requirement.** Observational values are distributed as benchmark data with their definition (λ∥ or λ_r), observer distance, species, energy or rigidity, and source table. A test that compares a model with them must state which conversion (if any) was applied, including ψ. The library must not convert λ_r to λ∥ silently.

The Chen et al. (2024) fit is supplied as the separate model `SEP-CHEN24` (Sections 5.2–5.3), with its published units and domain. The corresponding λ∥ follows from the definition λ∥ = 3κ∥/v (Eq. 1) once the particle speed is known; the species is not stated by the source (protons are implied), so it must be chosen explicitly (Section 14, U-2), and whether κ∥ or the derived λ∥ is reported is the caller's output choice (U-1).

## 5. Empirical SEP prescriptions used by transport codes

SEP transport and acceleration codes almost always prescribe the ambient scattering length by a separable power law, normalized at one rigidity (or momentum) and one heliocentric distance. This section records each code's published choice exactly as printed, then states how the choices map onto one generic form. The generic form is the rigidity power law of the parallel specification (PARALLEL rev. 1.4, Section 5.1, Eq. 13); it is not repeated here.

### 5.1 Generic separable form and the three normalization choices

The codes below fit

$$
\Lambda(r,X)=\Lambda_0\left(\frac{X}{X_0}\right)^{a}\left(\frac{r}{r_0}\right)^{b},
\tag{9}
$$

where Λ is one of three different quantities and X is one of three different momentum variables:

| Choice | Options | Consequence |
|---|---|---|
| Quantity Λ | λ∥; λ_r = λ∥ cos²ψ; "λ" not stated as either | λ_r constant in r implies λ∥ = λ_r/cos²ψ, which grows with r through ψ |
| Variable X | rigidity P (or R) = pc/(Ze); momentum pc; kinetic energy E | P/1 GV and pc/1 GeV coincide only for \|Z\| = 1 (Section 2, Eq. 4) |
| Reference | (X₀, r₀), usually 1 GV or 1 GeV and 1 au | λ at another energy follows only from the stated exponent |

**Implementation requirement.** A preset records all three choices. A value such as "0.3 au" is not portable between codes unless Λ, X, X₀ and r₀ are also carried.

### 5.2 Published code presets

Values are transcriptions; "as printed" marks wording that is ambiguous in the source. Machine-readable copies with locations are in `parameters/sep_code_presets.json` (Section 19).

| ID | Code / source | Printed form | Published parameters | Quantity | Status |
|---|---|---|---|---|---|
| SEP-PATH09 | PATH; Verkhoglyadova et al. 2009, Eq. (2) [M28] | λ = λ₀(pc/1 GeV)^α (r/1 AU)^β | λ₀ = 0.8 AU, α = 1/3, β = 2/3 ("chosen … in this study") | "λ"; parallel or radial not stated | Full text, 2 fetches |
| SEP-EPREM10 | EPREM; Schwadron et al. 2010 [M29] | λ∥ ∝ P^{1/3}, scaled by "a nominal value at 1 GV rigidity"; Eq. (2) (Sect. 2.2): D_μμ = (R₁/r)^{±3/2}(1−μ²)v/(2λ₀) with λ₀ "the parallel mean free path at R₁ = 1 AU" (sign of the exponent not legible) | 0.05 AU at 1 GV ("does a reasonable job with the event onsets") | λ∥ | Full text, 3 fetches (two of them copies of the journal PDF); the radial factor of Eq. (2) is not restated for the Sect. 3 SEP runs, and the sign of its 3/2 exponent could not be read (Section 14.3) |
| SEP-EPREM13 | EPREM; Kozarev et al. 2013, Eq. (3) [M30] | λ∥ = λ₀(pc/1 GeV)^{1/3}(R/1 AU)^{2/3} | λ₀ = 0.05 AU "at 1 AU and 1 GV proton rigidity" | λ∥ | Full text, 2 fetches; text says 1 GV, equation uses pc = 1 GeV |
| SEP-MFLAMPA19 | M-FLAMPA; Borovikov et al. 2019, Eq. (6.7) [M31] | λ_xx = λ₀(R/1 AU)(pc/1 GeV)^{1/3} | "λ₀ ∼ 0.1 ÷ 0.4 AU is a free parameter" | λ∥ (λ_xx along B) | Full text, 2 fetches (preprint) |
| SEP-MFLAMPA25 | M-FLAMPA; Liu et al. 2025, Eqs. (14)–(16) [M32] | λ∥ = λ₀(r/1 au)(pc/1 GeV)^{1/3}; D∥ = λ∥v/3 | λ₀ = 0.3 au (Table 1), "the mean free path for 1 GeV particles at 1 au" | λ∥ | Full text, 2 fetches |
| SEP-SOFIE24 | SOFIE; Zhao et al. 2024 [M33] | constant mean free path upstream of the shock | 0.3 AU in all nine events; sensitivity runs 0.05, 0.3, 1 AU (2013 Apr 11 event, Fig. 7) | λ∥: "the same constant parallel mean free path in all SEP events" (Sect. 4.5); elsewhere only "mean free path" | Full text (arXiv v1), 3 fetches |
| SEP-ZHANG23 | Zhang et al. 2023, Eqs. (35), (37) [M34] | λ_r constant in space; λ∥ = λ∥0(x) p^{2−q} | λ_r = 200 R☉ at 1 GV (also 20 R☉ in Fig. 7); q = 5/3; h₀ = 0.2 | λ_r | Full text, 2 fetches |
| SEP-PARASOL25 | PARASOL; Afanasiev et al. 2025, their Eq. 43 [M35] | λ⁰ = (0.1 au)(R/R_ref)^{2−q₀} | R_ref = 43 MV ("rigidity of a 1 MeV proton"), q₀ = 5/3 | λ∥ (background) | Full text, 2 fetches (preprint) |
| SEP-SOLPENCO05 | SOLPENCO; Aran et al. 2005 [M36] | λ ∝ P^{1/2} | 0.2 AU or 0.8 AU for 0.5 MeV protons | "proton mean free path" | Full text, 4 fetches (pre-2006 foundation) |
| SEP-SPARX15 | SPARX; Marsh et al. 2015 [M37] | isotropic scattering of pitch angle and gyrophase in the solar-wind frame; Poisson-distributed scattering times "determined by a prescribed value of mean free path" | λ = 0.3 AU (assumed for the database runs, Sect. 3.2); no energy or radial dependence stated | isotropic scattering length | Full text, 4 fetches (including the published version) |
| SEP-MARSH13 | Marsh et al. 2013 [M38] | isotropic scattering in the solar-wind frame at Poisson-distributed times with mean λ/v₀ (Sect. 2.2) | λ = 0.3, 1, 10 AU | isotropic scattering length | Full text, 3 fetches |
| SEP-HE11 | He, Qin & Zhang 2011, Eqs. (4)–(6) [M39] | λ_r constant; D^r_μμ = D₀vR^{−1/3}(\|μ\|^{q−1}+h)(1−μ²) | λ_r = 0.3 AU for 50 MeV protons ("corresponding to λ = 0.67 AU at 1 AU"); q = 5/3; h = 0.2 | λ_r | Full text, 1 fetch (equation reconstructed from a garbled PDF) |
| SEP-WANGQIN15 | Wang & Qin 2015, Eq. (3) [M40] | D_μμ = D₀vp^{q−2}{\|μ\|^{q−1}+h}(1−μ²) | λ∥ = 0.126 AU or 0.3 AU for 10 MeV protons at 1 AU; q = 5/3; h = 0.01 | λ∥ | Full text, 1 fetch |
| SEP-KUBO15 | Kubo et al. 2015, Eqs. (6), (8) [M41] | λ_rr = λ∥cos²Φ constant | test runs λ_rr = 0.2, 0.8, 1.4, 2.0 AU at 142 MeV (Fig. 2); 0.8 AU at 1 GV (Fig. 4); 2012 Jan 27 event: 0.2–2.0 AU in 0.3 AU increments compared with GOES data (no reference energy stated in that sentence; no single best-fit value); q = 5/3; h = 0.2 | λ_r | Full text, 3 fetches |
| SEP-PARADISE19 | PARADISE; Wijsen et al. 2019, Eqs. (8), (11) [M42] | constant λ^r∥ | 0.3 AU for 4 MeV protons; ε = 0.048 | λ_r | Full text, 2 fetches |
| SEP-LAITINEN16 | Laitinen et al. 2016 [M43] | λ∥ from QLT over a radially evolving spectrum (Section 11.4) | normalized to λ∥ = 0.3 AU for a 10 MeV proton at 1 AU; H = 0.1 | λ∥ | Full text, 2 fetches |
| SEP-LAITINEN18 | Laitinen et al. 2018 [M44] | as above | λ∥ = 0.1, 0.3, 1 AU for 10 MeV protons at 1 AU; only 10 MeV protons are simulated and no rigidity scaling is printed | λ∥ | Full text, 4 fetches (including the published version) |
| SEP-STRAUSS17G | Strauss et al. 2017 (GLE), Eqs. (10)–(11) [M45] | λ_rr = λ₀[r/r₀]^α, λ_rr = λ∥cos²Ψ | r₀ = 1 AU; α ∈ {−2,−1,0,1,2} explored (Appendix B); λ₀ free; 2 GV protons | λ_r | Full text, 3 fetches; equation numbers follow the arXiv PDF (the ar5iv rendering numbers them one lower) |
| SEP-STRAUSS15 | Strauss & Fichtner 2015 [M46] | constant λ∥ | λ∥ = 1 AU (Sect. III), 0.5 AU (Sect. IV); ~85 keV electrons; H = 0.05 | λ∥ | Full text (preprint), 3 fetches |
| SEP-DROGE16P | Dröge et al. 2016 (ICRC 2015 proceedings) [M47] | λ∥ = λ_r/cos²ψ(r) | 2010 Aug 7: λ_r = 0.12 AU, α = 0.08 (STEREO-B sector: 0.5 AU, α = 0.02); 2010 Aug 18: λ_r = 0.22/0.06/0.12 AU (STEREO-A/ACE/STEREO-B sectors); 65–105 keV electrons | λ_r | Full text, 4 fetches; the paper says the Aug 7 values were "assumed" and the STEREO-B-sector values "determined" |
| SEP-MINOSHIMA26 | Minoshima et al. 2026, Eq. (9) [M06] | λ∥ = ξ v/Ω_n as printed (Ω_n undefined; D-29), ξ spatially uniform | ξ assimilated; V_A = 100 km/s; ν_p = Ω_p/2π = 0.063 Hz | λ∥ | Full text, 2 fetches (preprint) |
| SEP-CHEN24 | Chen et al. 2024, Eq. (4) [M11] | κ∥ = (5.16±1.22)×10¹⁸ r^{1.17±0.08} E^{0.71±0.02} cm² s⁻¹ | r in AU, E in keV; 0.1–0.8 AU; 100 keV–1 GeV | κ∥ | Full text, 3 fetches |

Two further variants are recorded without a preset because the source does not fix a value: Young et al. (2021) write κ∥ = κ∥0(R_g/R_g0)^χ for EPREM [M48], with no λ₀ in the accessible text; and the CCMC M-FLAMPA page lets users set the far-upstream coefficient with no default.

### 5.3 Normalization details by preset

**SEP-PATH09.** The transport module uses "the Monte Carlo approach of Li et al. (2003, 2005)" with pitch-angle scattering, focusing and cooling; no D_μμ is printed and perpendicular diffusion is not included. Borovikov et al. (2019) attribute the (R/1 AU)(pc/1 GeV)^{1/3} form to Li et al. (2003) and the (R/1 AU)^{2/3} dependence to Zank et al. (2007); both attributions are secondary.

**SEP-EPREM13.** Kozarev et al. (2013) print Eq. (2), D_μμ = (1−μ²)v/(2λ∥), entering the transport equation as ∂_μ[(D_μμ/2)∂_μ f]. Section 6.6 shows that, if this is the operator, the standard relation returns 2λ∥. The version used "does not include perpendicular (cross-field) diffusion".

**SEP-MFLAMPA25.** The equivalent energy form printed by Liu et al. (2025, Eq. 16) is

$$
D_\parallel=\frac13\,c\,\lambda_0\left(\frac{r}{1\,\mathrm{au}}\right)
\left[\frac{E_k(E_k+2E_{p0})}{(1\,\mathrm{GeV})^2}\right]^{1/6}
\left[\frac{E_k(E_k+2E_{p0})}{(E_k+E_{p0})^2}\right]^{1/2}.
\tag{10}
$$

The first bracket is (pc/1 GeV)^{1/3} and the second is v/c, so Eq. (10) is identical to D∥ = vλ∥/3 with the preset form (algebraic check, fixture F-SEP-03). Liu et al. state that λ₀ = 0.3 au "is consistent with the results of Chen et al. (2024)", and they compare the dependence adopted in M-FLAMPA with the Chen et al. formula (their Eq. 17). Self-generated waves are not included in this paper.

**SEP-ZHANG23.** D₀(x) is set so that λ_r = λ∥cos²ψ is constant ("We follow Bieber (1994) to set the radial mean free path λ_r to be constant"). Momentum p is in GV, i.e. the code variable is momentum per charge; for protons it is rigidity. With 1 R☉ = 6.957 × 10⁸ m (IAU nominal) and 1 au = 1.495978707 × 10¹¹ m, 200 R☉ = 0.930093 au; this conversion is ours, not the paper's.

**SEP-PARASOL25.** The rigidity of a 1 MeV proton is 43.33 MV (Section 2), so the printed R_ref = 43 MV is rounded. Using 43 rather than 43.33 MV changes λ⁰ by the factor (43.33/43)^{1/3} = 1.0026 at fixed R (fixture F-SEP-05). The upstream and downstream modifications are in Section 8.3.

**SEP-HE11 family.** The printed pairs give the cos²ψ the authors used at 1 AU, as ratios of their own numbers: He et al. 2011 and He & Wan 2015 (λ_r = 0.3 AU, λ∥ = 0.67 AU) imply 0.448; He 2015 (0.15 AU, 0.34 AU) implies 0.441; He & Wan 2019 (0.25 AU, 0.5 AU and 0.28 AU, 0.56 AU) imply 0.500. The papers state the solar-wind speed (400 km/s in He & Wan 2015) but the rotation rate and resulting ψ were not recorded here; the ratios are therefore descriptive only [M39, M49, M50, M51].

**SEP-LAITINEN16.** The normalization fixes λ∥ only at (1 AU, 10 MeV); other radii and energies follow from the QLT integral over the prescribed turbulence (Section 11.4). There is no closed-form λ∥(R, r).

**SEP-SPARX15 and SEP-MARSH13.** Marsh et al. (2013, Sect. 2.2) print the mean time between scatterings as λ/v₀; Marsh et al. (2015) state only that the Poisson-distributed scattering times are determined by the prescribed mean free path. Isotropic redistribution of the velocity direction at the mean rate v/λ gives the hard-sphere result κ = vλ/3 (derived), so λ is the transport mean free path. The papers do not label λ as parallel or radial. Cross-field transport in these full-orbit codes comes from drifts.

**SEP-MINOSHIMA26.** Equation (9) is printed with Ω_n, which the paper does not define; the text describes λ∥ as "proportional to the local cyclotron radius", which implies λ∥ ∝ v/|B| if Ω_n is the local proton gyrofrequency (not stated; D-29). The authors take λ∥ ∝ |B|⁻¹ ∝ r^b with b = 2 for r ≪ A and b = 1 for r ≫ A, A = 1 AU.

**SEP-CHEN24.** See Section 4.5 for the defining equations. The published quantity is κ∥. The corresponding λ∥ follows from the definition λ∥ = 3κ∥/v (Eq. 1), which is not an additional physical assumption; it does require the particle speed, and therefore the species, which the paper does not state (Section 14, U-2). Whether a run reports κ∥ or the derived λ∥ is the caller's output choice (U-1). Fixture F-SEP-06 tabulates κ∥ and the derived proton λ∥.

### 5.4 Target-code modes "Chen2024AA" and "Tenishev2005AIAA"

**Chen2024AA.** No Chen et al. (2024) paper in A&A presenting an SEP mean-free-path or diffusion model was found. The searches of A&A were rate-limited, so this absence is not exhaustive. The Chen et al. (2024) model cited by SEP-modelling groups for this purpose is SEP-CHEN24 (ApJ 965, 61), cited as such by Liu et al. (2025) and Minoshima et al. (2026). Mapping the code label to SEP-CHEN24 is an inference that the code owner must confirm (Section 14, U-3).

**Tenishev2005AIAA.** The only matching publication found is Tenishev & Combi (2005), "Monte-Carlo Model for Dust/Gas Interaction in Rarefied Flows", AIAA-2005-4832 [M52]. Tenishev, Zhao & Sokolov (2022) cite it in Sect. 3.3.1 for the SEP adaptation of AMPS [M53]. The full text of the AIAA paper was not readable, and no mean-free-path formula attributed to it was found in any accessible source. This specification therefore defines no `Tenishev2005AIAA` equation. The code's own implementation, or the author, is the only source (Section 14, U-4). The D_μμ types that Tenishev et al. (2022) list for AMPS are in Section 6.7.

### 5.5 Reference values side by side

The following illustrates how different the presets are at a common point. It is evaluated by fixture F-SEP-01 with the stated conversions, for a proton at 1 au and pc = 1 GeV (P = 1 GV):

| Preset | Value at 1 au, pc = 1 GeV | Conversion applied |
|---|---|---|
| SEP-PATH09 | λ = 0.8 au | none |
| SEP-EPREM13 | λ∥ = 0.05 au | none |
| SEP-MFLAMPA25 | λ∥ = 0.3 au | none |
| SEP-ZHANG23 | λ_r = 0.930093 au | 200 R☉ in au; λ∥ requires ψ |

The values are not comparable without ψ (for SEP-ZHANG23) and without knowing whether SEP-PATH09's λ is parallel or radial. A comparison of SEP-CHEN24 with SEP-MFLAMPA25 at several energies is in fixture F-SEP-06.

## 6. Pitch-angle diffusion shapes and the λ∥ normalization

Focused-transport codes need D_μμ(μ), not only λ∥. Most of them prescribe a shape and fix its amplitude from a prescribed λ∥ (or λ_r). This section collects the shapes in use and gives the exact amplitude for each, so a library can convert between a λ preset (Section 5) and a D_μμ.

### 6.1 Normalization integral

The relation used by every code in this review is

$$
\lambda_\parallel=\frac{3v}{8}\int_{-1}^{1}\frac{(1-\mu^2)^2}{D_{\mu\mu}(\mu)}\,d\mu,
\qquad
\kappa_\parallel=\frac{v\lambda_\parallel}{3}=\frac{v^2}{8}\int_{-1}^{1}\frac{(1-\mu^2)^2}{D_{\mu\mu}(\mu)}\,d\mu .
\tag{11}
$$

It is the diffusion limit of the focused-transport equation with D_μμ defined as the coefficient in ∂_μ(D_μμ ∂_μ f) (PARALLEL rev. 1.4, Section 7.1). Equivalent printings occur: (3v/4)∫₀¹ for even D_μμ (Laitinen et al. 2016, Eq. 4; Shalchi, Yan & Lazarian 2005, Eqs. 35–39 [M54]), (v²/4)∫_{μ_min}^1 with a cut-off (Chen et al. 2024, Eq. 2), and D_zz = v²⟨(1−μ²)²/(4D_μμ)⟩_μ with ⟨·⟩_μ = ½∫₋₁¹ (Borovikov et al. 2019, Eq. 3.7). These agree with Eq. (11) except for the deliberate cut-off in Chen et al.

If a code writes its operator with a different factor (for example ∂_μ[(D_μμ/2)∂_μ f]), the D_μμ inserted in Eq. (11) must be the coefficient of ∂_μ(·∂_μ f), i.e. D_μμ/2 in that example (Section 6.6). Section 21.5 gives the Itô stochastic differential equation that reproduces the operator ∂_μ(D_μμ∂_μ f) in a particle solver (Eq. 46).

**Implementation requirement.** A shape object stores (i) the printed functional form, (ii) the operator convention, and (iii) whether its amplitude is set from λ∥, from λ_r (with ψ), or from turbulence. Conversion between amplitude and λ uses Eq. (11) evaluated exactly or by the closed forms below.

### 6.2 Regularized power-law shape (q-form)

The most common shape, after Beeck & Wibberenz (1986), is

$$
D_{\mu\mu}=D_0\,(1-\mu^2)\left(|\mu|^{q-1}+H\right),\qquad 1<q<2,\ H\ge 0 .
\tag{12}
$$

Inserting Eq. (12) in Eq. (11) gives exactly

$$
D_0=\frac{3v}{4\lambda_\parallel}\,I(q,H),
\qquad
I(q,H)=\int_0^1\frac{1-\mu^2}{\mu^{q-1}+H}\,d\mu,
\qquad
I(q,0)=\frac{2}{(2-q)(4-q)} .
\tag{13}
$$

For H = 0 this is D₀ = 3v/[2λ∥(2−q)(4−q)], which is the printed normalization of Agueda & Vainio (2013, Eq. 2: ν₀ = 6v/[2λ∥(4−q)(2−q)] with D_μμ = (ν₀/2)…) and of Minoshima et al. (2026, Eq. 8 at V_A = 0). For q = 5/3, I(5/3, 0) = 18/7. Values of I(5/3, H) for the published H are fixture F-PA-01:

| H | Used by | Status |
|---|---|---|
| 0.01 | Wang & Qin 2015 ("h = 0.01 is chosen for non-linear effect of pitch angle diffusion at μ = 0 in the solar wind"); Wang et al. 2014; Qin & Wang 2015; AMPS Type V [M40, M55, M56, M53] | Full text |
| 0.05 | Strauss & Fichtner 2015 ("chosen in an ad hoc fashion"); Strauss et al. 2017 (GLE) [M45]; Lang et al. 2024; Lavasa et al. 2026; Dröge et al. 2010 [M57] (secondary, via Wijsen et al. 2019 and Strauss et al. 2017) | Full text except Dröge |
| 0.1 | Laitinen et al. 2016 ("enhances scattering between pitch angle hemispheres") | Full text |
| 0.2 | He et al. 2011, He & Wan 2015, He 2015; Kubo et al. 2015 ("[Qin et al., 2006]"); Zhang et al. 2023 (h₀ = 0.2; "not very sensitive to h₀ unless h₀ ≪ 0.05") | Full text |

At a fixed λ∥, a larger H lowers I and therefore lowers D₀; the H values above are not interchangeable.

**Printed factor (1−μ)².** Lang et al. (2024, Eq. 3) print D_μμ = D₀(|μ|^{s−1}+H)(1−μ)² with s = 5/3 "the inertial range (Kolmogorov) spectral index of the slab turbulence component", and Lavasa et al. (2026, Eq. 3) print D₀(|μ|^{q−1}+H)(1−μ)² with q = 5/3. The conventional factor, used by the sources they cite, is (1−μ²). With (1−μ)² the coefficient is not even in μ and vanishes only at μ = 1. This specification transcribes the printing and implements only (1−μ²) under the ID `PA-QFORM`; using the printed (1−μ)² requires the explicit variant `PA-QFORM-LANG-PRINTED`, which is a user decision (Section 14, U-5). With the printed factor, D_μμ(−1) = 4D₀(1 + H) ≠ 0, so the pitch-angle flux D_μμ∂_μf does not vanish at μ = −1 by itself; a solver that runs the variant must impose the no-flux condition at μ = −1 explicitly (derived). The normalization integral stays finite, since (1−μ²)²/(1−μ)² = (1+μ)².

### 6.3 ε-form (Agueda et al.; PARADISE)

$$
D_{\mu\mu}=\frac{\nu_0}{2}\left(\frac{|\mu|}{1+|\mu|}+\varepsilon\right)(1-\mu^2),
\qquad
\lambda_\parallel=\frac{v}{\nu_0}\,\varphi(\varepsilon),
\qquad
\varphi(\varepsilon)=\frac32\int_0^1\frac{(1-\mu^2)(1+\mu)}{\mu+\varepsilon(1+\mu)}\,d\mu .
\tag{14}
$$

The definition of φ follows from Eq. (11); Agueda et al. (2010, Eqs. 1–4) print λ∥ = (v/ν₀)φ(ε) and the asymptote φ(ε) = (3/2)ln(1/ε) for ε ≪ 1 [M13, M58]. The asymptote is only the leading term; at ε = 0.048 and 0.01 it differs from the exact φ (fixture F-PA-02). **Implementation requirement:** evaluate φ exactly.

Values in use:

- ε = 0.048 in PARADISE (Wijsen et al. 2019, Sect. 4; Afanasiev et al. 2025, Eq. 41). Wijsen et al. obtained it by least-squares matching of the ε-form to the q-form "for H = 0.05". PARADISE writes the amplitude as ν₀(|Q|/m)^{2−q}B^{−q}v^{q−1}/2 (Wijsen et al. 2019, Eq. 8), so λ∥ ∝ (mv/|Q|)^{2−q}B^q at fixed ν₀.
- ε = 0.01 in Pacheco et al. (2019, Eq. 4), which prints the first term as |μ|/(1−|μ|), not |μ|/(1+|μ|). With (1−|μ|), D_μμ = (ν₀/2)[|μ|(1+|μ|) + ε(1−μ²)], which does not vanish at |μ| = 1, unlike every other shape in this section. Transcribed as printed and flagged (Section 14, D-4); the ID `PA-EPS` implements |μ|/(1+|μ|) only.

### 6.4 Rigidity-scaled amplitudes

Several codes write the amplitude with an explicit momentum or rigidity factor so that one D₀ covers all energies:

| Source | Printed form | Implied λ∥ scaling |
|---|---|---|
| He et al. 2011 Eq. (6); He & Wan 2015 Eq. (7) | D^r_μμ = D_μμ/cos²ψ = D₀vR^{−1/3}(\|μ\|^{q−1}+h)(1−μ²) | λ∥ ∝ R^{1/3} |
| Wang & Qin 2015 Eq. (3); Kubo et al. 2015 Eq. (6) | D_μμ = D₀vp^{q−2}(\|μ\|^{q−1}+h)(1−μ²) | λ∥ ∝ p^{2−q} |
| Zhang et al. 2023 Eqs. (35), (37) | D_μμ = D₀(x)p^{q−2}(1−μ²)(\|μ\|^{q−1}+h₀), p in GV | λ∥ = λ∥0(x)p^{2−q} (printed) |
| Qin & Wang 2015 Eqs. (2)–(3) | D_μμ = D₀v(R_Lk_min)^{s−2}(μ^{s−1}+h)(1−μ²), D₀ = (δB_slab/B₀)²π(s−1)k_min/(4s) | slab-QLT inertial range |
| Wang et al. 2014 Eqs. (3)–(4); AMPS Type V | same, printed with R^{s−2} and "R = pc/(\|q\|B₀) is the maximum particle Larmor radius" | see below |
| Wijsen et al. 2019 Eq. (7) | D_μμ = (Cπ/2)(\|Q\|/m)^{2−q}B^{−q}v^{q−1}(\|μ\|^{q−1}+H)(1−μ²) | QLT inertial range at fixed C |

In the He et al. form, v·R^{−1/3} gives λ∥ ∝ R^{1/3} from Eq. (11) (v cancels); in the Wang & Qin form, λ∥ ∝ p^{2−q}, which equals R^{1/3} for protons at q = 5/3. The hidden scaling matters when a single-energy λ value (Section 5.2) is used at other energies.

**Dimensional check.** In Qin & Wang (2015) the factor (R_Lk_min)^{s−2} is dimensionless and D₀ has units of k_min, so D_μμ has units of s⁻¹. Wang et al. (2014) and AMPS Type V print R^{s−2} with R defined as a Larmor radius; that product has units m^{−1/3} s⁻¹ for s = 5/3 and is dimensionally consistent only if R means R_Lk_min. This specification implements the Qin & Wang (2015) form (Section 14, D-6).

**Consistency with closed-form QLT (algebraic check).** With h = 0, the Qin & Wang (2015) form inserted in Eq. (11) gives

$$
\lambda_\parallel=\frac{6s}{\pi(s-1)(2-s)(4-s)}\left(\frac{B_0}{\delta B_{\rm slab}}\right)^2\frac{(R_Lk_{\min})^{2-s}}{k_{\min}},
\tag{15}
$$

which is exactly the R = R_Lk_min ≪ 1 limit of the Teufel & Schlickeiser (2003) proton form, Eq. (21) in Section 7.2. The factor s in the numerator comes from the TS2003 spectrum that is flat below k_min (Section 7.3).

### 6.5 Dynamical-turbulence shape (Dröge 2003 form)

Minoshima et al. (2026, Eq. 8) print, attributing it to Dröge (2003),

$$
D_{\mu\mu}(\mu)=\frac{3v}{2(4-q)(2-q)\lambda_\parallel}\left(\sqrt{\mu^2+(V_A/v)^2}\right)^{q-1}(1-\mu^2),\qquad 1<q<2 .
\tag{16}
$$

For V_A = 0 this is the q-form with H = 0 and exact normalization. For V_A > 0 the term (V_A/v) fills the resonance gap, and Eq. (11) applied to Eq. (16) returns a value smaller than the nominal λ∥ in the prefactor; the ratio depends on V_A/v and q (fixture F-PA-03). **Implementation requirement:** `PA-DROGE-VA` must report whether λ∥ is the prefactor parameter or the integral value. The Dröge (2003) original was not accessible.

### 6.6 Isotropic shapes

**Isotropic pitch-angle diffusion.** D_μμ = (ν/2)(1−μ²) gives λ∥ = v/ν from Eq. (11).

**EPREM printing.** Kozarev et al. (2013, Eq. 2) print D_μμ = (1−μ²)v/(2λ∥) and the operator ∂_μ[(D_μμ/2)∂_μ f]; Young et al. (2021) print the same operator form. The coefficient of ∂_μ(·∂_μ f) is then (1−μ²)v/(4λ∥) = (ν/2)(1−μ²) with ν = v/(2λ∥), and Eq. (11) returns 2λ∥ (derivation in this review; fixture F-PA-04). Either the code's λ∥ is half the standard transport mean free path, or the extraction lost a factor. Not resolved (Section 14, D-2). **Implementation requirement:** an EPREM-compatible mode must state which of the two it reproduces.

**Hard-sphere (isotropic redistribution).** Full-orbit codes (Marsh et al. 2013, 2015; Kelly et al. 2012 [M59]; Battarbee et al. 2018) redistribute the velocity direction isotropically at Poisson-distributed times; Marsh et al. (2013) print the mean time λ/v₀, and Marsh et al. (2015) state that the times are determined by the prescribed mean free path. The spatial diffusion coefficient is κ = vλ/3 (derived), so λ is the transport mean free path. This is not a D_μμ model and needs no shape.

### 6.7 Kolmogorov shape without gap filling (M-FLAMPA) and AMPS types

Borovikov et al. (2019, Eqs. 5.3, 6.1–6.6) derive, for an Alfvén-wave spectrum I₋ + I₊ = I_C k^{−5/3},

$$
D_{\mu\mu}=\frac{v}{\lambda_{\mu\mu}}(1-\mu^2)|\mu|^{2/3},
\qquad
\lambda_{\mu\mu}=\frac{4}{\pi}\frac{B^2/\mu_0}{I_C}r_L^{1/3},
\qquad
\lambda_{xx}=\frac{54}{7\pi}\frac{B^2/\mu_0}{I_C}r_L^{1/3}.
\tag{17}
$$

Inserting the shape in Eq. (11) gives λ_xx/λ_μμ = 27/14 exactly, in agreement with the printed (54/7π)/(4/π) (algebraic check, fixture F-PA-05). This is the q-form with q = 5/3, H = 0. If the spectrum is normalized so that ∫_{k₀}^∞(I₋+I₊)dk = (δB)²/μ₀, then I_C = (2/3)(δB)²k₀^{2/3}/μ₀ and λ_xx becomes the downstream form printed by Liu et al. (2025, Eq. 18), with coefficient 81/(7π) (Section 8.4). Borovikov et al. also print approximate forms with coefficients 0.9 (λ_xx) and 0.5 (λ_μμ) in terms of the wave energy density w₋+w₊ and L_max (Eqs. 6.10–6.11). They state the convention, k₀⁻¹ = L_max(R)/2π with L_max(R) ∼ 0.03R, and evaluate "the numerical factor 81/(7π(2π)^{2/3}) ≈ 0.92". That factor equals 1.0817 (0.92 is close to its reciprocal, 0.9244), and the corresponding λ_μμ factor is 0.5609 (fixture F-PA-05); the printed 0.9 and 0.5 are therefore not consistent with the stated convention (Section 14, D-7). This specification uses the exact expressions.

**AMPS D_μμ types.** Tenishev, Zhao & Sokolov (2022, Table 1) list the shapes implemented in AMPS. The table was read from OCR-quality text; |μ| versus μ and the exponent placement in Type I were not double-checked.

| Type | Attributed to | Printed form and parameters |
|---|---|---|
| I | Borovikov et al. 2019 | D_μμ = (v/λ_μμ)(1−μ²)μ^{2/3}; λ_μμ = 0.5[(B²/μ₀)/(w₋+w₊)](L²_max r_L0)^{1/3}(pc/1 GeV)^{1/3} |
| II | le Roux & Webb 2009 | D₀ = (π/8)A²Ω²l_b(1−μ²); l_b = 0.03 r⁻¹ AU upstream; ε = 1; ⟨δB⊥²⟩ = 0.1B₀²(1 AU/r) |
| III | Hu et al. 2017; Jokipii 1966; Zhao & Li 2014 | D_μμ = (π/4)(1−μ²)Ω₀kP^slab(k)/B²; P(k) = A_β λ_c(δB)²/[1+(kλ_c)^β] with ∫₀^∞P dk = (δB)²; λ_c = 10⁹ m; "typical values at 1 au: k_L = 2.0×10⁻⁷ m⁻¹, k_R = 1.0×10⁻¹⁰ m⁻¹"; (δB/B₀)² = 0.05 |
| IV | Wang & Qin 2004 | D_μμ = D₀vp^{q−2}(μ^{q−1}+h)(1−μ²); q = 5/3 |
| V | Qin et al. 2013 | D_μμ = (δB_slab/B₀)²[π(s−1)/(4s)]k_min vR^{s−2}(μ^{s−1}+h)(1−μ²); s = 5/3; "k_min = 1/l_slab = 33AU^1" (as printed); h = 0.01 |

The extraction did not show which type produced the paper's figures.

### 6.8 QLT D_μμ from a supplied spectrum

Codes that compute D_μμ from a spectrum use magnetostatic slab QLT; the parallel specification gives its full treatment (PARALLEL rev. 1.4, Section 8). Printed forms differ only through the spectrum convention:

| Source | Printed D_μμ | Spectrum convention |
|---|---|---|
| Chen et al. 2024 Eq. (3) | (π/4)Ω₀(1−μ²)f_res P(f_res)/B₀² | measured frequency PSD of B_N, Taylor hypothesis |
| Ding et al. 2022 Eq. (8) [M60] | 2π²Ω²(1−μ²)g^slab(k∥)/(B²vμ), k∥ = Ω/(v\|μ\|) | g^slab = [C(ν)/2π]l_slab δB²_slab(1+k²l²_slab)^{−ν}, C(ν) = Γ(ν)/[2√πΓ(ν−½)] |
| Laitinen et al. 2016 Eq. (2) | πΩ²(1−μ²)S∥(−(r_Lμ)⁻¹)/(v\|μ\|B²) | S∥ as defined in the paper |
| Zhang et al. 2023 Eq. (36) | π²Ω²(1−μ²)W⊥(k_res)/(B²v\|μ\|) | W⊥ as defined in the paper |
| Wijsen et al. 2019 Eq. (6) | (π/2)(1−μ²)(Ω/B)²P(k=Ω/(\|μ\|v))/(v\|μ\|), P = Ck^{−q} | one-sided power law |
| Afanasiev et al. 2015 Eq. (2) [M61] | (π/2)Ω\|k_res\|I_w,res(1−μ²)/R² as printed in the arXiv version, followed by "where B is the mean magnetic field" (B² intended; D-30), k_res = Ω/(vμ) | δB² = ∫₋∞^∞ I_w dk |
| Strauss et al. 2017 Eq. (20) [M62] | damped QLT with Lorentzian resonances, γ = αV_Ak∥ (Eq. 21) | g^slab of the paper |

**Implementation requirement.** A spectrum adapter must carry its normalization (one- or two-sided, per component or total, wavenumber or frequency) and convert it before using any of these prefactors. PARALLEL rev. 1.4, Section 8.7 gives the conversion rules.

### 6.9 Focusing-corrected parallel mean free path

He & Wan (2013, Eqs. 2–3) use [M63]

$$
\lambda_\parallel=\frac{3L^3}{\lambda_{\parallel,0}^2}\left[\frac{\lambda_{\parallel,0}}{L}-\tanh\!\left(\frac{\lambda_{\parallel,0}}{L}\right)\right],
\qquad
L(r,\theta,V)=\frac{r\,(V^2+\Omega^2r^2\sin^2\theta)^{3/2}}{V\,(2V^2+\Omega^2r^2\sin^2\theta)},
\tag{18}
$$

where L is the Parker-field focusing length (Ω here is the solar rotation rate). For λ∥,0 ≪ L, λ∥ → λ∥,0 (x − tanh x ≈ x³/3). In floating-point arithmetic the difference x − tanh x cancels for small x = λ∥,0/L; Section 21.4 gives the equivalent series (Eq. 43), which changes only the evaluation. This is a published isotropic-scattering result; Shalchi and co-workers show that focusing corrections depend on the transport definition (PARALLEL rev. 1.4, Section 2.3), so it is supplied only as the explicit option `MFP-FOCUS-HW13`.

## 7. Closed-form turbulence-based parallel mean free paths

Two closed forms account for nearly all turbulence-driven λ∥ in recent SEP and GCR modelling: the Teufel & Schlickeiser (2003) asymptotics, assembled into one continuous expression by the North-West University group, and the Zank et al. (1998) fit to slab QLT. Both are magnetostatic or weakly dynamical slab results. They require turbulence inputs (Section 11) and differ in spectrum conventions, which this section states explicitly.

### 7.1 Spectrum conventions

| Convention | Slab spectrum | Variance normalization | Used by |
|---|---|---|---|
| TS2003 / EB2013 / Lang 2024 | g = g₀k_min^{−s} (\|k∥\| ≤ k_min); g₀\|k∥\|^{−s} (k_min ≤ \|k∥\| ≤ k_d); g₀k_d^{p−s}\|k∥\|^{−p} (\|k∥\| ≥ k_d) | (δB)² = 8π∫₀^∞g dk∥, total over components | Teufel & Schlickeiser 2003 [M64]; Engelbrecht & Burger 2013b [M65]; Lang et al. 2024 |
| TS2002 | as above, but g = 0 for \|k∥\| ≤ k_min | same | Teufel & Schlickeiser 2002 [M66] |
| Bendover (Shalchi school; iPATH) | g = [C(ν)/2π] l δB²_slab (1+k²l²)^{−ν}, s = 2ν, C(ν) = Γ(ν)/[2√π Γ(ν−½)] | 8π∫₀^∞g dk = δB²_slab | Shalchi et al. 2006, 2008 [M08, M67]; Ding et al. 2022 |
| Zank 1998 fit | bendover spectrum, written with the correlation length λ_s = 2πC(ν)l | per-component variance B²_x,slab (Zank 1998), or total ⟨b_s²⟩ (Chhiber 2017) | Section 7.5 |

Lang et al. (2024, Eq. 9) give the normalization of the TS2003 spectrum including the dissipation range,

$$
g_0=\frac{\delta B_{\rm slab}^2\,k_{\min}^{s-1}(s-1)}{8\pi}\left[s+\frac{s-p}{p-1}\left(\frac{k_{\min}}{k_d}\right)^{s-1}\right]^{-1},
\tag{19}
$$

(the same expression is Engelbrecht & Burger 2013b, Eq. 2), and relate k_min to the slab correlation length by (Lang et al. 2024, Eq. 13, attributed to Engelbrecht's 2012 thesis)

$$
\lambda_{\rm slab}=\frac{\pi(s-1)}{2k_{\min}}\left[s+\frac{s-p}{p-1}\left(\frac{k_{\min}}{k_d}\right)^{s-1}\right]^{-1}
\;\xrightarrow{\;k_{\min}/k_d\to0\;}\;
\frac{\pi(s-1)}{2s\,k_{\min}} .
\tag{20}
$$

### 7.2 Teufel & Schlickeiser (2003) based continuous form

Definitions (TS2003, Eq. 7; Lang et al. 2024): R_L = P/(B₀c) with P = pc/|q| in volts (SI); R = R_Lk_min; Q = R_Lk_d; a = v/(α_D V_A), b = a/2; 0 ≤ α_D ≤ 1 is the dynamical parameter (DT: exp[−α_D|k∥|V_A|t|]; RS: exp[−α_D²k∥²V_A²t²]).

**Protons (and other ions).** Engelbrecht & Burger (2013a, 2014) attribute the continuous construction to Burger et al. (2008) [M68, M69], whose original equation was not read. Lang et al. (2024, Eq. 10 with 𝒦 = 0), Engelbrecht & Burger (2014, Eq. 9) [M70] and TS2003 (Eq. 24 with the R ≪ 1 and high-rigidity table entries) give the same expression:

$$
\lambda_\parallel=\frac{3s\,R_L^2k_{\min}}{4\pi(s-1)}\left(\frac{B_0}{\delta B_{\rm slab}}\right)^2
\left[1+\frac{8}{(2-s)(4-s)}\,\frac{1}{R^{s}}\right].
\tag{21}
$$

The three printings are algebraically identical (EB2014 writes the bracket as [1/(4√π) + 2R^{−s}/(√π(2−s)(4−s))] with prefactor 3s/(√π(s−1))·R²/k_m). Equation (21) is a sum of two asymptotes, not a uniformly valid QLT result:

- R ≪ 1: λ∥ → [6s/(π(s−1)(2−s)(4−s))](B₀/δB_slab)²R^{2−s}/k_min, i.e. ∝ P^{1/3} for s = 5/3. This is Eq. (15).
- R ≫ 1: λ∥ → [3s k_min R_L²/(4π(s−1))](B₀/δB_slab)², i.e. ∝ P² (TS2003, Eqs. 31–34; for RS the paper states it is "exactly the same result as for the DT-model").

Because Eq. (21) adds asymptotes, it is not exact at R ∼ 1: compared with the exact magnetostatic slab QLT for the same spectrum (flat below k_min, k^{−s} above, no dissipation range), it is too large by up to 26% (exact/form = 0.7938 at R = 3.03 for s = 5/3; derived, fixture F-QLT-07). That exact integral has a closed form, Eq. (45) in Section 21.4, which is a validation reference and not a replacement for the published Eq. (21). The dynamical parameter does not appear in Eq. (21); TS2003 Table 3 shows this regime holds for R ≪ 1 ≪ a ≪ Q. Engelbrecht & Burger (2014) state that dissipation-range effects are neglected in this form.

**Electrons.** Lang et al. (2024, Eqs. 10–12) add a dissipation-range term:

$$
\lambda_\parallel=\frac{3s\,R_L^2k_{\min}}{4\pi(s-1)}\left(\frac{B_0}{\delta B_{\rm slab}}\right)^2
\left[1+\frac{8}{(2-s)(4-s)}\frac{1}{R^{s}}+\frac{4\,\mathcal K}{Q^{\,p-s}R^{s}}\right],
\tag{22}
$$

$$
\mathcal K_{\rm RS}=\left(\frac{\sqrt\pi}{\Gamma(p/2)}+\frac{1}{p-2}\right)\left(\frac a2\right)^{p-2},
\qquad
\mathcal K_{\rm DT}={}_2F_1\!\left(1,\frac{1}{p-1};\frac{p}{p-1};-\frac{a}{f_1Q}\right)\frac{a}{f_1},
\qquad
f_1=\frac{2(p-s)}{\pi(p-2)(2-s)} .
\tag{23}
$$

Lang et al. note that p > 2 is required for convergence. The f₁ of Lang et al. equals the TS f₁ = 2/(p−2) + 2/(2−s) divided by π, so −a/(f₁Q) is the TS argument −πa/(f₁^TS Q).

- **RS term.** Equation (23) is algebraically identical to TS2003 Table 4 (R ≪ 1 ≪ Q ≪ b) and to Engelbrecht & Burger (2013b, Eq. 3), using b = a/2; both print 1/Γ(p/2). TS2002 Table 4 prints the equivalent 2^{p−1}Γ((p+1)/2)/[√πΓ(p)], which equals 1/Γ(p/2) by the duplication formula.
- **DT term (unresolved printing difference).** The three sources differ:

| Source | Q-power multiplying 𝒦 | ₂F₁ argument |
|---|---|---|
| TS2002, TS2003 Table 3 (R ≪ 1 ≪ Q ≪ a) | Q^{3−s} | −πa/(f₁^TS Q) |
| Engelbrecht & Burger 2013b, Eq. (4) | Q^{p−s} | −πa/(f₁^TS Q^{p−2}) |
| Lang et al. 2024, Eqs. (10), (12) | Q^{p−s} | −a/(f₁Q) (= −πa/(f₁^TS Q)) |

All three coincide for p = 3, the TS reference value. They differ for the fitted p = 3.57 of Lang et al. (Table 6). No published erratum was found, and the typeset ApJ version could not be checked (Section 14, D-1). **Implementation requirement:** the DT electron term is implemented only with an explicit variant flag (`TS-DT-Q3MS`, `TS-DT-EB13`, `TS-DT-LANG24`); there is no default. Fixture F-QLT-04 evaluates all three at the Lang et al. Table 6 parameters, labelled as evaluations of printed formulas, not as published values.

### 7.3 TS2002 versus TS2003

TS2002 sets g = 0 below k_min; TS2003 makes the spectrum flat there. For the same δB²_slab this raises the inertial-range λ∥ by the factor s. With the Bieber et al. (1994) parameter set used in both papers — B₀ = 4.12 nT, k_min = 10⁻¹⁰ m⁻¹, k_d = 2 × 10⁻⁵ m⁻¹, s = 5/3, p = 3, v_A = 33.5 km/s, α = 1 (TS2003, Eq. 50) — the printed DT proton approximations are (every entry is λ/λ₀):

| Paper | Mid-rigidity | High rigidity | Low rigidity |
|---|---|---|---|
| TS2002, Eq. (65); Appendix C, Eq. (C.2) | ≈ 0.0106 AU (r/MV)^{1/3}, 10⁻¹ MV ≪ r ≪ 10⁴ MV | Eq. (C.2): ≈ 1.96 × 10⁻⁶ AU (r/MV)², r ≫ 10⁴ MV | ≈ 0.0062 AU (r/MV), printed without an exponent (first power), r ≪ 10⁻¹ MV |
| TS2003, Eq. (51) | ≈ 0.018 AU (r/MV)^{1/3} | ≈ 2.62 × 10⁻¹⁰ AU (r/MV)², r ≫ 10⁴ MV | ≈ 0.010 AU (r/MV), r ≪ 10⁻¹ MV |

Here λ₀ = (B₀/δB)² is dimensionless, so the coefficients are λ(δB/B₀)². The R ≪ 1 term of Eq. (21) reproduces 0.01775 AU (TS2003) and, without the factor s, 0.01065 AU (TS2002); the R ≫ 1 term reproduces 2.615 × 10⁻¹⁰ AU; R = 1 occurs at r = 12 351 MV, printed as "r < 1.23 × 10⁴ MV" (fixture F-QLT-01). The TS2002 high-rigidity coefficient 1.96 × 10⁻⁶ AU (r/MV)² belongs to the spectrum with no power below k_min and was not checked against a derivation here; TS2002 evaluate it with δB ≈ 0.33B₀ (λ₀ = 10) as λ ≃ 2 × 10⁻⁵ AU (r/MV)² (Eq. C.3), state that the resulting λ > 2000 AU is "unphysical", and attribute it to the sharp cut-off at k_min. The electron approximations of TS2002 Eq. (67) and TS2003 Eqs. (52), (54) were not legible in the accessible extraction and are not transcribed.

### 7.4 Published parameter sets for Eqs. (21)–(23)

| Set | B₀ (nT) | δB²_slab (nT²) | s | p | k_d (km⁻¹) | k_min (km⁻¹) | V_A (km/s) | α_D | Source status |
|---|---|---|---|---|---|---|---|---|---|
| TS-BIEBER94 | 4.12 | not fixed (results printed as λ/λ₀) | 5/3 | 3 | 2 × 10⁻² | 10⁻⁷ | 33.5 | 1 | TS2003 Eq. 50, printed as k_d = 2 × 10⁻⁵ m⁻¹ and k_min = 10⁻¹⁰ m⁻¹ |
| LANG24-T5 nominal 1 au | 6.27 ± 2.30 | 3.18 ± 2.15 | 1.69 ± 0.04 | 2.61 (+0.96/−0.60) | (4.28 ± 4.16) × 10⁻³ | (2.28 ± 1.22) × 10⁻⁷ | 56.77 ± 22.59 | — | Lang Table 5 |
| LANG24-T6-DT (electrons) | 3.97 | 5.33 | 1.65 | 3.57 | 8.44 × 10⁻³ | 2.28 × 10⁻⁷ | 56.77 | 0.5 | Lang Table 6 |
| LANG24-T6-RS (electrons) | 3.97 | 5.33 | 1.65 | 2.61 | 8.44 × 10⁻³ | 3.50 × 10⁻⁷ | 79.36 | 1.0 | Lang Table 6 |
| LANG24-T6-P (protons) | 3.97 | 1.03 | 1.65 | 3.57 | 0.12 × 10⁻³ | 2.28 × 10⁻⁷ | 79.36 | 0.5 | Lang Table 6 |
| EB13b | 5 (at Earth) | model | 5/3 | 2.6 (tested 2.3, 2.6, 2.8) | from Leamon et al. 2000 | model | model | 0–1 | EB2013b text |

Lang et al. Table 5 also gives δB² = 15.91 ± 10.74 nT² with a 20/80 slab/2D split, λ_slab = (2.81 ± 1.51) × 10⁶ km and λ_2D = (1.10 ± 0.49) × 10⁶ km. Table 5 prints the unit of B₀ as "m nT" (apparent typo) and states that the s and p ranges are 1σ while the other ranges are 2σ. Tables 5 and 6 were confirmed by an independent re-reading (Section 20).

### 7.5 Zank et al. (1998) fit

As printed in the ICRC 1999 paper of the same group, which reproduces Zank et al. (1998) [M71, M72]:

$$
\lambda_\parallel=3.1371\,\frac{B^{5/3}}{B^2_{x,\rm slab}}\left(\frac{P}{c}\right)^{1/3}\ell_{\rm slab}^{2/3}
\left[1+\frac{7A/9}{(q+1/3)(q+7/3)}\right],
\tag{24}
$$

$$
A=(1+s^2)^{5/6}-1,\qquad
q=\frac{5s^2/3}{1+s^2-(1+s^2)^{1/6}},\qquad
s\equiv0.746834\,\frac{R_L}{\ell_{\rm slab}} .
\tag{25}
$$

B²_x,slab is the variance of one slab component. Chhiber et al. (2017, Eqs. 4–7) and Perri et al. [M73] print the coefficient 6.2742 with the total slab variance ⟨b_s²⟩ and the slab correlation length λ_s in place of ℓ_slab; Perri et al. state the units: "B in nT, ⟨b_s²⟩ in nT², P in V, c in m/s, λ_s in m", giving λ∥ in km. The validity "at rigidities ranging from from 10 MV to 10 GV" is stated by Chhiber et al. and Perri et al.; the JGR original was not read. Here s and q are auxiliary symbols of the fit, not the spectral index or the TS quantity Q.

**Origin of the constants (derivation in this review; fixture F-QLT-02).** For the bendover spectrum with ν = 5/6, C(5/6) = 0.1188622 and the correlation length is λ_s = 2πC(5/6)l = 0.7468343 l, which is the constant in s. The exact inertial-range (R_L ≪ l) slab-QLT result for this spectrum is λ∥ = [3/(2πC(2−s)(4−s))](B/δB_slab)²R_L^{1/3}l^{2/3} with s = 5/3 here the spectral index; rewritten with λ_s, the coefficient is 6.27421 for the total variance δB²_slab and 3.13710 for the per-component variance δB²_slab/2. Thus 3.1371 and 6.2742 are the same fit in two variance conventions, provided λ_s (Chhiber) and ℓ_slab (Zank) both denote the correlation length 2πCl. The bracket in Eq. (24) reproduces the exact magnetostatic QLT integral for this spectrum to within 2.2%: the largest deviation is +2.16% at R_L/l = 5.787, and the ratio returns to 1.0014 at R_L/l = 100 (fixture F-QLT-02 tabulates it). The exact integral itself equals the inertial-range result times a Gauss hypergeometric function, Eq. (44) in Section 21.4 (derived), which gives an independent check of these numbers.

**Le Roux et al. (1999) coefficient.** Quenby & Webber (2015) [M74] quote λ∥ = 2.433(B²/b²_rsl)(P/(cB))^{1/3}λ_sl^{2/3} × F, "derived by Le Roux et al. (1999) for a composite field model", with F "very close to unity". Since (B²)(P/(cB))^{1/3} = B^{5/3}(P/c)^{1/3}, 2.433 can be compared directly with 3.1371; it is 0.7756 times smaller. It could not be reconciled with the conventions above and is not implemented (Section 14, D-8).

**Published evaluations.** Chhiber et al. (2017, Table 1): λ∥ = 0.29 AU (standard case), 0.21 AU (doubled turbulence) and 0.40 AU (halved turbulence) for a 100 MeV (445 MV) proton at 1 AU in the ecliptic. Perri et al.: λ∥ from 2 × 10⁻² AU at 1 MV to 200 AU at 10⁵ MV at 1 AU (case D1), scaling as P^{0.33} below ≈ 2 × 10³ MV and P^{1.31} above. These quoted end values and slopes are not mutually consistent (0.02 AU × 2000^{0.33} × 50^{1.31} ≈ 41 AU, not 200 AU; derived check), and the paper gives 0.05–0.2 AU (Sec. 3.1) and 0.02–0.2 AU (Sec. 3.2) for 100 MeV (Section 14, D-31).

### 7.6 Nonlinear and simulation results usable as benchmarks

No new closed-form λ∥ was found in the 2006–2026 nonlinear literature that was accessible; nonlinear D_μμ (NADT, nonlinear damping, SOQLT) is integrated numerically through Eq. (11), as in the parallel specification (PARALLEL rev. 1.4, Sections 9–11).

- **Shalchi et al. (2006)**, NADT model, Table 2 at 1 AU: s = 5/3, p = 3, v_A = 33.5 km/s, l_slab = 0.030 AU, l_2D = 0.1 l_slab, k_slab = k_2D = 3 × 10⁶ AU⁻¹, B₀ = 4.12 nT, δB/B₀ = 1, slab fraction 0.2, 2D fraction 0.8. The paper reports λ∥ within 0.08–0.3 AU for 0.5–5000 MV but tabulates no λ values. Shalchi, Lazarian & Schlickeiser (2008) find the slab-only nonlinear-damping λ∥ at medium rigidity about 5 times smaller than the composite result of Shalchi et al. (2006).
- **Hussein, Tautz & Shalchi (2015)** [M75], test-particle λ∥ for magnetostatic turbulence with δB = B₀, s = 5/3:
  - Table 3 (isotropic turbulence), λ∥/l₀ at R_L/l₀ = 0.01, 0.05, 0.1, 0.5, 1.0, 5.0, 10.0: 0.51, 1.05, 1.35, 3.55, 7.1, 82.0, 330.
  - Table 5 (NRMHD), λ∥/l⊥ at R_L/l⊥ = 0.01, 0.05, 0.1, 1.0, 5.0, 10: 1.05, 2.1, 3.2, 31, 700, 1875.
  - Table 2 (slab/2D, 20/80): λ∥/l_slab = 0.90, 1.88, 2.7, 8.5, 15.0, 87.0, 255. The R_L/l_slab column headers were lost in extraction; they are not assumed here, so Table 2 is not distributed as a benchmark.
- **Subashchandar et al. (2025)**: SOQLT with PSP spectra gives κ∥ ∝ (δB/B₀)^{−2.13} for δB/B₀ ≈ 0.2–0.6, plus the radial scalings of Section 4.5. No tabulated λ.

### 7.7 Dissipation-range onset

The electron λ∥ is controlled by k_d (Lang et al. 2024: "the electron pMFP is extremely sensitive to the onset and shape of the dissipation range"). Published prescriptions:

- Strauss et al. (2020; numbered (9) in the arXiv PDF and (13) in the ar5iv rendering), Eq. (5) here: k_d = (2π/V_sw)(a + bΩ_i) with a = 0.2 Hz and b = 1.76, after Leamon et al. (2000).
- Engelbrecht & Strauss (2018, Eq. 17) [M76], as extracted: k_d,LH ≈ Ω_p/(V_A + 3v_th^p), with 1 AU values B = 5 nT, n = 7 cm⁻³, T_p = 0.95 × 10⁵ K, T_e = 1.46 × 10⁵ K, V_A = 41 km/s, v_th^p = 40 km/s (Table 1). Their other onset estimates (Eqs. 18–21) were not legible in the extraction.
- Engelbrecht & Burger (2013b, Table 2): 1 AU breakpoint frequencies of 0.456, 0.452, 0.298 and 0.178 Hz for the Leamon et al. (2000) models.

## 8. Shock-region and self-generated-wave prescriptions

Near a CME-driven shock the ambient prescriptions of Section 5 are replaced or capped. The codes reviewed use four approaches: a Bohm limit, a steady-state self-generated-wave solution, fits to self-consistent wave simulations, and a turbulence-energy-based λ downstream. None of these is a single universal formula.

### 8.1 Bohm limit

Zhang et al. (2023, Sect. 2.4) state: "We choose the Bohm limit for it or κ₁ = vp/(3qB₁)" for the upstream coefficient, with "typically κ₂ ≪ κ₁" downstream, and "We have only partially implemented this feature by applying the Bohm diffusion limit for calculating shock acceleration of particle sources." The Bohm form is

$$
\kappa_{\rm Bohm}=\frac{v\,r_g}{3},\qquad r_g=\frac{p}{|q|B},\qquad \lambda_{\rm Bohm}=r_g .
\tag{26}
$$

It is the same comparison model as PARALLEL rev. 1.4, Section 13.1. Kozarev & Schwadron (2016) [M77] use the hard-sphere tensor κ = κ∥cos²θ_BN + κ⊥sin²θ_BN with κ⊥ = κ∥/[1 + (λ∥/r_g)²] in an analytic shock model; that model is not EPREM.

### 8.2 Steady-state self-generated waves upstream of a shock

Afanasiev, Battarbee & Vainio (2015, Eqs. 7 and 12) give the steady-state upstream wave intensity and mean free path for isotropic scattering in the shock frame:

$$
I_w(x,k)=\frac{1}{3\pi}\,\frac{\Omega B^2|k|^{-3}}{(u_1-V_A)(x+x_0)},
\qquad
\lambda(x,v)=\frac{3(u_1-V_A)}{v}\,(x+x_0),
\tag{27}
$$

where x is the distance upstream, u₁ the upstream flow speed in the shock frame and x₀ a constant. Inserting I_w at k_res = Ω/v into their isotropic rate ν = πΩ²I_w,res/(vB²) (Eqs. 3–4) gives ν = v²/[3(u₁−V_A)(x+x₀)] and hence λ = v/ν, which is the second expression (consistency check in this review). Equivalently κ = vλ/3 = (u₁−V_A)(x+x₀). Their simulations start from outward Alfvén waves ∝ k^{−3/2} "calculated assuming the initial mean free path λ=1R_⊙ for 100 keV protons", with n₀ = 3.6 × 10¹² m⁻³, B₀ = 3.4 × 10⁻⁵ T, V_s = 1500 km/s, E_inj = 18.1 keV and a scattering-centre compression ratio of 4. Battarbee et al. (2010) [M78] set "The ambient 100 keV proton mean free path … to λ₀ = 1R_⊙ at r₀ = 1.5R_⊙".

Ng & Reames ("Shock Acceleration of Solar Energetic Protons: The First 10 Minutes") [M79] cap the wave intensity at I_σ(k) ≤ I_sat(k) = B²/(3πk) and state that a "Prescribed ambient Iσ ∝ k^{−1.5} gives λ = 0.23 AU at 1 MeV (0.68 AU at 100 MeV)". A pure R^{2−q} scaling with q = 1.5 from the 1 MeV value would give 0.74 AU at 100 MeV; the reason for the difference is not stated in the accessible text. Ng, Reames & Tylka (2012) [M80] initialize inward waves as I_R− = I_L− = 0.1 I_R+. The wave-intensity normalization of Ng & Reames was not compared with that of Afanasiev et al., so the two caps are not converted into each other here.

### 8.3 PARASOL shock-region fits

Afanasiev, Wijsen & Vainio (2025) fit SOLPACS self-consistent simulations for a shock at 1 au and feed the result to PARADISE:

$$
\text{upstream:}\quad \lambda_\parallel=\min\!\left[\lambda^0,\ \lambda^{(P)}\right],\qquad
\lambda^{(P)}=\Lambda_0\left[1+\left(\frac{x}{\Delta x_0}\right)^{q_0}\right]^{1/q_0},
\tag{28}
$$

$$
\text{downstream:}\quad \lambda_\parallel=\min\!\left[\lambda^0,\ \Lambda_0\,e^{-d/d_0}\right],\qquad d_0=1\,R_\odot,
\tag{29}
$$

with the background λ⁰ = (0.1 au)(R/R_ref)^{2−q₀}, R_ref = 43 MV, q₀ = 5/3 (SEP-PARASOL25; their Eqs. 21, 43, 44). Related fitted relations printed in the paper are K = 0.22 V_A/(ε_inj Ω₀) (their Eq. 28), C₂ = 0.66 and α₂ = 0.47 (their Eqs. 34–35). Λ₀ and Δx₀ are the paper's Λ(E) (their Eq. 18) and Δx(E) (their Eq. 29) evaluated for a shock at 1 au. The scattering-centre compression ratio is r_sc = r(1 − M_A⁻¹). Free parameters are the injection efficiency ε_inj (example 5 × 10⁻⁴) and Δt (example 10 h). The expressions for Λ(E) and Δx(E) were not transcribed in this review, so `SHOCK-PARASOL` requires them to be supplied from the source (Section 14, U-6).

### 8.4 M-FLAMPA downstream and numerical floor

Liu et al. (2025, Eqs. 18–20) use, downstream of the shock,

$$
\lambda_\parallel=\frac{81}{7\pi}\left(\frac{B}{\delta B}\right)^2\frac{r_{L0}^{1/3}}{k_0^{2/3}}\left(\frac{pc}{1\,\mathrm{GeV}}\right)^{1/3},
\qquad r_{L0}=\frac{1\,\mathrm{GeV}}{c\,e\,B},\qquad k_0=\frac{2\pi}{L_{\max}(r)},\quad L_{\max}(r)=0.4\,r,
\tag{30}
$$

$$
D_\parallel=\max\{D_\parallel,\ D_{\min}\},\qquad D_{\min}=0.1\,R_s\times10^5\ \mathrm{m\,s^{-1}} .
\tag{31}
$$

δB is taken from the AWSoM-R Alfvén-wave intensity ("a Kolmogorov spectrum with an index of -5/3", Zhao et al. 2024). D_min is a numerical floor "as used in Sokolov et al. (2004) and Borovikov et al. (2018)", compensating for the shock width eroded by the mesh. Borovikov et al. (2019) give instead L_max(R) ∼ 0.03R, and a trial k₀ = const ∼ 0.1/R_S (Borovikov et al. 2018); these are different published choices, not a conflict to be resolved. Equation (30) follows from Eq. (17) when ∫_{k₀}^∞ I_C k^{−5/3}dk = (δB)²/μ₀ (Section 6.7). Borovikov et al. (2019, Eq. 6.9) print the same expression; Liu et al. (2025) is used as the reference printing.

### 8.5 Other codes

- **PATH** (Verkhoglyadova et al. 2009): near the shock "A self-consistent treatment of this interaction yields a value of the parallel diffusion coefficient" from "the quasi-linear formulation by Gordon et al. (1999)"; no explicit κ or I(k) and no Bohm limit are printed. The foundations (Zank et al. 2000; Rice et al. 2003; Li et al. 2003, 2005) were not accessible.
- **iPATH** (Hu et al. 2017; Ding et al. 2020, 2022) [M81, M82, M60]: parallel and perpendicular coefficients from QLT and NLGC. Ding et al. (2020, Eq. 8) use κ∥ ∼ v^{4/3}B^{5/3}r^{(2/3)α+3} and κ⊥ ∼ v^{10/9}B^{−7/9}r^{(8/9)α−1} with l_slab ∼ l_2D ∼ r^α, α = 0.8, δB² ∼ r⁻³, κ⊥/κ∥ = 0.0099 at 1 AU for 1 MeV protons and δB²/B² = 0.5 at 1 AU. Ding et al. (2022, Eq. 13) use κ∥ ∼ v^{4/3}B^{5/3}r^{(2/3)a_slab−γ}, κ⊥ ∼ v^{10/9}B^{−7/9}r^{γ/3+2a_slab/9+2a_2D/3}, with γ = −3.5, a_slab = a_2D = 1.0, κ⊥/κ∥ = 0.03 for 10 MeV protons at 1 au, δB²/B² = 0.15 at 1 au, l_slab = 10⁹ m. These exponents follow from slab QLT (λ∥ ∝ (B/δB)²l^{2/3}r_L^{1/3}) and NLGC (κ⊥ ∝ [v(δB_2D/B)²l_2D]^{2/3}κ∥^{1/3}) (algebraic check in this review). The shock-region coefficient of Hu et al. (2017) was not accessible.
- **EPREM** (Kozarev et al. 2013): the paper describes no modification of λ∥ near the shock.
- **SEPMOD** (Luhmann et al. 2010) [M83]: "The subsequent transport in the CORHEL interplanetary fields is scatter-free", with an option for "Monte-Carlo-like pitch angle scattering"; no λ value.
- **SOLPENCO** (Aran et al. 2005): an optional foreshock region confines particles ahead of the shock; no λ formula for it in the accessible text.
- **Subashchandar et al. (2025)** measure upstream rise times at a PSP shock and write κ_rr = W₁V_shΔt and κ∥ = κ_rr/cos²θ_Bn (inline). With their numbers (V_sh = 556 km/s, W₁ = 286 km/s, θ_Bn ≈ 66°, Δt = 5.2 min at 75 keV and 12 min at 750 keV) these give κ∥ ≈ 3.0 × 10¹⁸ and ≈ 6.9 × 10¹⁸ cm² s⁻¹ (fixture F-SH-03; the paper only plots them). These are foreshock values and must not be used as ambient constraints.

## 9. GCR modulation prescriptions and published parameter sets

GCR modulation codes specify K∥ (usually written K or κ, in cm² s⁻¹ or AU² s⁻¹) rather than λ∥, with λ∥ = 3K∥/v. The functional forms are few; the normalization conventions are many. This section gives each form exactly as printed, the published parameter sets, and the normalization facts needed to use them. Machine-readable tables are in `parameters/` (Section 19).

### 9.1 Tensor projection

Potgieter (2013, Eq. 22) [M84] gives the heliocentric tensor elements used by all of these codes:

$$
K_{rr}=K_\parallel\cos^2\psi+K_{\perp r}\sin^2\psi,\quad
K_{\theta\theta}=K_{\perp\theta},\quad
K_{\phi\phi}=K_{\perp r}\cos^2\psi+K_\parallel\sin^2\psi,\quad
K_{\phi r}=K_{r\phi}=(K_{\perp r}-K_\parallel)\cos\psi\sin\psi .
\tag{32}
$$

The SEP quantity λ_r = λ∥cos²ψ (Section 4.2) is the first term of K_rr only.

### 9.2 NWU smooth broken power law

Potgieter et al. (2014, Eq. 5) [M85]:

$$
K_\parallel=(K_\parallel)_0\,\beta\left(\frac{B_0}{B}\right)\left(\frac{P}{P_0}\right)^{a}
\left[\frac{(P/P_0)^{c}+(P_k/P_0)^{c}}{1+(P_k/P_0)^{c}}\right]^{\frac{b-a}{c}},
\qquad P_0=1\ \mathrm{GV},\ B_0=1\ \mathrm{nT}.
\tag{33}
$$

"(K∥)₀ is a constant in units of 10²² cm².s⁻¹"; no numerical value is printed. b = 1.95 and c = 3.0 ("determines the smoothness of the transition"); a and P_k vary with time (Table 1). Exact properties (algebra in this review):

- At P = P₀ the bracket equals 1, so K∥(1 GV) = (K∥)₀β(B₀/B). For P ≪ P_k, K∥ ∝ βP^a; for P ≫ P_k, K∥ ∝ βP^b.
- Since v = βc, λ∥ = 3K∥/v = [3(K∥)₀/c](B₀/B)G(P), where G is the rigidity factor of Eq. (33). λ∥ does not depend on β, and λ∥(P₁)/λ∥(P₂) = G(P₁)/G(P₂) at fixed B.

This is the smooth broken rigidity model of the parallel specification (PARALLEL rev. 1.4, Section 6.1), which states the dimensionally normalized form and its derivatives. Variants in the same group:

| Paper | Printed exponent symbols | Printed exponent of the bracket | Notes |
|---|---|---|---|
| Potgieter 2013, Eq. 23 | a, b, c | (b−a)/c | B_n/B_m with B_n = 1 nT; b = 1.95, c = 3.0 |
| Potgieter et al. 2015 (electrons), Eq. 9 [M86] | a, b (b₁), c | (b₁−a)/c, reconstructed from a garbled PDF | b = 1.55, c = 3.5; K_⊥r has its own high-rigidity slope b₂ |
| Vos & Potgieter 2015, Eqs. 3 and 5 [M87] (κ∥ = κ₀βFG in Eq. 3; G(P) in Eq. 5) | a, b, c, P_k | not legible (minus sign lost in extraction) | separate b for λ∥ and λ⊥ (Table 2) |
| Vos & Potgieter 2016, Eq. 4 [M88] | a, c | printed (c−a)/c, confirmed by a re-fetch | b is mentioned in the text but absent from Eq. (4) |
| Aslam et al. 2019, Eq. 14 [M89]; Aslam et al. 2023, Eq. 5 [M90] | c₁, c₂∥, c₃ | (c₂∥−c₁)/c₃ | (K∥)₀ unit 10²² cm² s⁻¹ |
| Aslam et al. 2021, Eq. 5 [M91] | c₁, c₂∥, c₃ | printed "2∥ − 1" over c₃ (Eq. 5) and (c₂⊥ − c₃)/c₃ (Eq. 6) | (K∥)₀ unit 6 × 10²⁰ cm² s⁻¹ |
| Raath 2015 (MSc thesis), Eq. 3.31 [M92] | a, b, c | (b−a)/c | unit 6 × 10²⁰ cm² s⁻¹; P′₀ = 1 GV, B₀ = 1 nT; the prefactor symbol βc is Raath's notation for v/c |

The printed exponents of Vos & Potgieter (2016) and Aslam et al. (2021) are transcribed as printed. Using (b−a)/c for them is an interpretation that a configuration must state explicitly (Section 14, D-9).

**Corti et al. (2019, Eq. 9)** [M93] normalize rigidity to the break:

$$
k_\parallel=k_\parallel^0\,\beta\,\frac{1\,\mathrm{nT}}{B}\left(\frac{R}{R_k}\right)^{a}\left[1+\left(\frac{R}{R_k}\right)^{s}\right]^{\frac{b-a}{s}},
\qquad R_k=4.3\ \mathrm{GV},\ s=2.2,
\tag{34}
$$

with k∥⁰ "in units of 6×10²⁰ cm²/s", k = vλ/3, k⁰_⊥,r = 0.02, k⁰_⊥,θ = 0.01 multiplied by the polar function u(θ) = 3/2 + ½tanh[8(|θ − π/2| − 35°)] (their Eqs. 10–11; perpendicular slopes a⊥, b⊥ fitted separately), R_A = 0.55 GV and k_A⁰ = 1. With c = s and the same a, b, Eqs. (33) and (34) have the same rigidity shape; the normalizations are related exactly by

$$
k_\parallel^0=(K_\parallel)_0\,\frac{(R_k/1\,\mathrm{GV})^{b}}{\left[1+(R_k/1\,\mathrm{GV})^{s}\right]^{(b-a)/s}}
\tag{35}
$$

when both are expressed in the same unit (identity derived in this review; fixture F-GCR-03). The two published unit choices (10²² versus 6 × 10²⁰ cm² s⁻¹) must be converted separately.

**Luo et al. (2019, Eq. 14)** [M94] use a double power law with coefficients "in units of 10²⁰ cm²/s"; the exponent placement, including the power of B, was not legible, so their Table 1 is distributed with the form flagged as unknown.

### 9.3 Published NWU-family parameter sets

**Protons, PAMELA 2006–2009 (Potgieter et al. 2014, Table 1).** Only values that changed are tabulated; b = 1.95 and c = 3.0 throughout.

| Parameter | 2006 | 2007 | 2008 | 2009 |
|---|---|---|---|---|
| α (tilt, °) | 15.7 | 14.0 | 14.3 | 10.0 |
| B at Earth (nT) | 5.05 | 4.50 | 4.25 | 3.94 |
| λ∥ at Earth, 100 MV (AU) | 0.04 | 0.06 | 0.09 | 0.12 |
| a | 0.56 | 0.48 | 0.39 | 0.28 |
| P_k (GV) | 4.0 | 4.0 | 4.0 | 4.2 |
| P_A0 (GV) | 1/√10 | 1/√10 | 1/√10 | 1/√40 |

Section 5 of the final arXiv version (v3) states that at 1 GV λ∥ "increased by a factor of ~2.3, from ~0.13 AU in 2006, to ~0.3 AU in 2009"; versions v1 and v2 printed "~23" and "~30 AU". Evaluating Eq. (33) with the Table 1 parameters gives λ∥(1 GV)/λ∥(100 MV) = 1/G(0.1 GV), so the tabulated 100 MV values imply λ∥(1 GV) = 0.146 AU (2006) and 0.230 AU (2009), a factor 1.57 (derived, fixture F-GCR-02). The approximate values in the text differ from these (Section 14, D-10).

**Protons, PAMELA 2006–2009 (Vos & Potgieter 2015, Table 2; "MFP values at 1.0 GV at the Earth").**

| Parameter | 2006e | 2007e | 2008e | 2009e |
|---|---|---|---|---|
| λ∥ (AU) | 0.742 | 0.888 | 0.970 | 1.204 |
| a | 0.91 | 0.88 | 0.86 | 0.80 |
| b for λ∥ | 2.10 | 2.10 | 2.10 | 2.10 |
| b for λ⊥ | 1.58 | 1.58 | 1.58 | 1.58 |
| c | 2.60 | 2.40 | 3.0 | 2.2 |
| P_k (GV) | 4.00 | 4.05 | 4.08 | 4.30 |

Their Table 1 gives α, B_e and r_TS for seven half-year periods (CSV). The constants B₀ and P₀ are illegible in the 2015 extraction ("=B 10 nT and =P 10 GV"); the 2016 paper by the same authors prints B₀ = 1 nT and P₀ = 1 GV. A normalizing parameter with the same symbol in two publications is not automatically the same coefficient (D-19).

**Electrons, PAMELA 2006–2009 (Potgieter et al. 2015, Table 1; "At Earth at 100 MeV").** The λ rows are printed with headers "λ∥ (AU) × 10⁻¹" and "λ⊥r (AU) × 10⁻³" and values 1.01–1.35 and 2.85–3.78; whether the numbers are to be multiplied by 10⁻¹ and 10⁻³ is not stated, and they are distributed as printed. Other values: a = 0.0, b₁ = 1.55, b₂ = 1.24, c = 3.5, P_k = 0.33–0.34 GV, d = 1.2.

**Positrons, PAMELA 2006–2009 (Aslam et al. 2019, Table 1; λ∥ at Earth at 1 GV).** λ∥ = 0.438, 0.465, 0.506, 0.539, 0.557, 0.574, 0.593 AU (2006b–2009b); c₁ = 0, c₂∥ = 2.25, c₂⊥ = 1.688, c₃ = 2.50 (2.70 for 2009a–b), P_k = 0.565–0.585 GV, K_A0 = P_A0 = 0.90, d_⊥θ = 6.

**Values in running text (Aslam et al. 2021, 2023).**

| Paper | Species | (K∥)₀ printed | Unit as printed |
|---|---|---|---|
| Aslam et al. 2021 | e⁻ and e⁺ | 32.13 (2006), 34.29 (2009), 31.38 (2015) | 6 × 10²⁰ cm² s⁻¹ |
| Aslam et al. 2023 | protons | 49.0 (May 2011), 26.49 (end of reversal), 30.52 (May 2015) | 10²² cm² s⁻¹ |
| Aslam et al. 2023 | electrons | 29.84 (May 2011), 21.09 (Mar–Apr 2014), 25.21 (May 2015) | 10²² cm² s⁻¹ |

For electrons, Aslam et al. (2023) give P_k = 0.50–0.60 GV in Sec. 6.4 and "P_k = 0.4 GV" in Sec. 7 (Section 14, D-20); for protons, P_k varies "between 2.50 to 4.0 GV". The 2021 and 2023 electron values are numerically similar under units that differ by a factor 16.7; the texts do not resolve this (Section 14, D-11). The values must be stored with their unit and never compared as bare numbers. Slopes and drift parameters printed in the same texts are in the CSV.

**AMS-02 protons per Bartels rotation (Corti et al. 2019, Table 2, A < 0, partial).** Eight rows (BR 2426–2445) and two rows of Table 3 (A > 0) were transcribed; the full Table 3 was not obtained. k∥⁰ = 110–130 (unit 6 × 10²⁰ cm² s⁻¹); λ⊥ at 1 GV and 5 GV at Earth is in AU (Table 2, footnote e). The table marks fits at the lower or upper edge of the parameter grid with a symbol; these are transcribed as "≤" and "≥" with a flag.

**AMS-02 transient events (Luo et al. 2019, Table 1).** Seven Forbush-decrease and GMIR rows with a, b, c and κ₀∥ = 1105–1615, κ₀⊥θ, κ₀⊥r (unit 10²⁰ cm² s⁻¹); P_k = 4.3 GV, a⊥ = 0.8, b⊥ = 1.0, c⊥ = 2.2.

### 9.4 HelMod

Boschini et al. (2018, ASR, Eq. 5) [M95] and Della Torre et al. (2016, Eq. 2) [M96]:

$$
K_\parallel=\frac{\beta}{3}\,K_0\left[\frac{P}{1\,\mathrm{GV}}+g_{\rm low}\right]\left(1+\frac{r}{1\,\mathrm{AU}}\right),\qquad P>1\ \mathrm{GV},
\tag{36}
$$

with K₀ in AU² GV⁻¹ s⁻¹ (Boschini 2018 ASR). As printed, the bracket is dimensionless while K₀ carries GV⁻¹. The numerical-use convention is K∥ [AU² s⁻¹] = (β/3)·K₀·(P/GV + g_low)·(1 + r/AU) with K₀ inserted as the printed number; a dimensionally typed implementation names this convention and keeps the printed unit as metadata, and does not convert K₀ as though the bracket carried an extra factor of GV (Section 14, D-21). HelMod v4 (Boschini et al. 2019, Eq. 2) [M97] replaces (1 + r/1 AU) by (R_c + R/1 AU), R_c "tuned to describe radial GCR intensity gradients in the inner heliosphere".

K₀ at low activity (tilt below about 50°) is a cubic in the smoothed sunspot number, and at high activity an exponential in the McMurdo neutron-monitor count rate (Boschini 2018 ASR, Eqs. 6–7, Tables 1–2):

$$
K_0^{\rm SSN}=c_0+c_1\,{\rm SSN}+c_2\,{\rm SSN}^2+c_3\,{\rm SSN}^3,
\qquad
K_0^{\rm NMCR}=p_0\exp\!\left(p_1\,{\rm NMCR}+p_2\,{\rm NMCR}^2\right).
\tag{37}
$$

| Polarity / phase | c₀ | c₁ | c₂ | c₃ | rms |
|---|---|---|---|---|---|
| A<0 ascending | 0.0003059 | −2.51 × 10⁻⁶ | 1.284 × 10⁻⁸ | −2.838 × 10⁻¹¹ | 0.1097 |
| A<0 descending | 0.0002876 | −3.715 × 10⁻⁶ | 2.534 × 10⁻⁸ | −5.689 × 10⁻¹¹ | 0.1400 |
| A>0 ascending | 0.0002262 | −5.058 × 10⁻⁷ | – | – | 0.1153 |
| A>0 descending | 0.0002267 | −7.118 × 10⁻⁷ | – | – | 0.1607 |

| Station | p₀ | p₁ | p₂ | rms |
|---|---|---|---|---|
| MCMU | 0.003753 | −0.04791 | 0.0001365 | 0.100 |
| OULU | 0.001354 | −0.10070 | 0.0007697 | 0.094 |

The SSN validity is stated as 2.2 ≤ SSN ≤ 266.9 in the text and about 10 ≲ SSN ≲ 165 in the Fig. 1 caption; the NMCR fit range is about 210 ≲ NMC ≲ 300; only McMurdo was used. g_low has maximum 0.3 at low activity, decreasing to 0 at high activity, by "an empirical function" not shown. ρ = K⊥/K∥ = 0.06. Della Torre et al. (2016) give g_low = 0.2 (low activity) and ρ ≈ 0.058; Boschini et al. (2017, ApJ) [M98] give g_low = 0 at maximum and 0.3 at minimum. The arithmetic conversion 1 AU² s⁻¹ = 2.2380 × 10²⁶ cm² s⁻¹ gives c₀(A<0 asc.) ≈ 6.846 × 10²² cm² s⁻¹ GV⁻¹ (fixture F-GCR-05).

**HelMod v4 transition function** (Boschini et al. 2019, Eq. 21, Table 2): F(α_t) = (F_min + F_max)/2 − [(F_min − F_max)/2]tanh[(α_t − α₀)/s], α_t the tilt angle of the heliospheric neutral sheet.

| Parameter | F_min | F_max | α₀ (asc) | s (asc) | α₀ (desc) | s (desc) |
|---|---|---|---|---|---|---|
| P₀,d (GV) | 0.5 | 4 | 73 | 1 | 63 | 10 |
| P₀,NS (GV) | 0.5 | SSN*/50 | 73 | 1 | 63 | 10 |
| K_c, q>0, A>0 | 3 | 1 | 40 | 18 | 53 | 5 |
| K_c, q>0, A<0 | 1 | 1 | – | – | – | – |
| K_c, q<0, A<0 | 3 | 1 | 40 | 18 | 53 | 5 |
| K_c, q<0, A>0 | 0.7 | 1 | 47 | 5.8 | 58.4 | 5.8 |
| g_low (e⁺/e⁻) | 0.4 | 0 | 67 | 20 | 45 | 10 |
| g_low (p, ions) | 0.5 | 0 | 60 | 9 | 45 | 10 |
| R_c | 4 | 1 | 60 | 9 | 45 | 10 |

The e⁺/e⁻ row label was read as "g_low, e⁺/e⁻" in one fetch and "K_c e⁺; e⁻" in another (identical values); it is flagged. ρ ≈ 0.065 (protons, ions), 0.05 (e±). The "Boschini et al. 2018a" that v4 cites for K₀(t) (its Eqs. 6–7) is, by the v4 reference list, Boschini et al. (2018, ASR 62, 2859) [M95]; its Eqs. 6–7 are Eq. (37) (equation numbers read in arXiv:1704.03733v1).

**Bobik et al. (2012)** [M99], the earlier form: K∥ ≈ βk₁(r,t)K_P(P,t)[B⊕/(3B)] (Eq. 6), K_P ≈ P above a threshold of 0.4–1.015 GV and about constant below, K_F = c₁ + c₂SSN⁻¹ + c₃SSN + c₄SSN² (Eq. 13) for 10 ≲ SSN ≲ 165, with Table 1:

| | A<0 asc | A<0 desc | A>0 asc | A>0 desc |
|---|---|---|---|---|
| c₁ | +0.0001686 | +8.872 × 10⁻⁵ | +2.39708 × 10⁻⁴ | +2.28037 × 10⁻⁴ |
| c₂ | +0.001488 | +0.001874 | — | — |
| c₃ | — | — | −8.28987 × 10⁻⁷ | −1.00984 × 10⁻⁶ |
| c₄ | −3.164 × 10⁻⁹ | — | — | — |

The units of K₀ and K_F are not stated. ρ_k = 0.05; K⊥θ = 10K⊥r for θ ≲ 30° or ≳ 150°. These coefficients belong to a different parameterization from Eq. (37) and must not be mixed with it.

### 9.5 Other GCR prescriptions

| ID | Source | Printed form | Published values | Status |
|---|---|---|---|---|
| GCR-STRAUSS11 | Strauss et al. 2011, Eq. 4 [M100] | λ∥ = λ₀(P/P₀)(1 + r/r₀) for P ≥ P₀; λ₀(1 + r/r₀) for P < P₀ | λ₀ = 0.15 AU (VALUE), P₀ = 1 GV, r₀ = 1 AU; κ⊥θ = κ⊥r = 0.01κ∥; electrons | Full text ×2; "≥" missing in extraction; κ relation printed "κ∥ = v/3λ∥" |
| GCR-EFFENBERGER12 | Effenberger et al. 2012, Eqs. 22–24 [M101] | κ∥ = κ∥0β(p/p₀)(B_e/B)^{a∥}; κ⊥1 = κ⊥0β(p/p₀)(B_e/B)^{a⊥}; κ⊥2 = ξκ⊥1 | κ∥0 = 0.9 × 10²² cm²/s (VALUE), p₀ = 1 GeV/c, a∥ = 0.75, a⊥ = 0.97, κ⊥0 = 0.1κ∥0, ξ = 2, κ_A = 0 | Full text ×2 |
| GCR-WANG19 | Wang et al. 2019, Eq. 4 [M102] | K∥ = (1/30)kβB_E/B for R < 0.1 GV; (1/3)kβR B_E/B for R ≥ 0.1 GV | k = k₀ · 3.6 × 10²² cm²/s; k₀ = (B_c/⟨B_E⟩)^n, n ≈ 2; K⊥ = 0.02K∥ | Full text ×2; B_c only in figures |
| GCR-TOMASSETTI17 | Tomassetti 2017, Sec. IV [M103] | K∥ = κ⁰ · 10²² βℛ/GV/(3B/B₀) cm² s⁻¹; κ⁰(t) = aφ_d⁻¹(t) + b | B₀ ≅ 3.4 nT "sets the HMF at r₀ = 1 AU" (final version; v1 printed "3.4 nT AU²"); K⊥ ≅ 0.02K∥; Table II (below) | Full text ×2 |
| GCR-PERUGIA21 | Fiandrini et al. 2021, Eq. 8 [M104] | K∥ = K₀(β/3)(R/R₀)^a/(B/B₀)·[((R/R₀)^h + (R_k/R₀)^h)/(1 + (R_k/R₀)^h)]^{(b−a)/h} | K₀ "in units of 10²³ cm²s⁻¹" (v1; "of the order of" in v2); R₀ = 1 GV; R_k = 3 GV; h = 3; B₀ at Earth; ξ⊥r = ξ⊥θ = 0.02; R_A = 0.5 GV | Full text ×3; results differ between versions |
| GCR-PERUGIA25 | Tomassetti et al. 2025, Eq. 5 [M105] | as GCR-PERUGIA21 | "h≡0.01"; R_k = 3 GV; K₀ "of the order of 10²³"; R_A ≅ 0.5 GV | Full text ×2 |
| GCR-JIANG23 | Jiang et al. 2023, Eq. 8 [M106] | K∥ = K₀β(B₀/\|B\|)(R/R₀)^a[((R/R₀)^m + (R_k/R₀)^m)/(1 + (R_k/R₀)^m)]^{(b−a)/m} | R_k = 3 GV, m = 3.0, R₀ = 1 GV, R_A = 0.5 GV; Sec. 4.3 averages K₀ = 0.3 × 10²³ (unit printed as cm⁻² s⁻¹ sr⁻¹ GeV⁻¹), a = 1.80, b = 0.989, B₀ = 5 nT, α = 70° | Full text ×2 |
| GCR-DUAN25 | Duan et al. 2025, Eqs. 3.5–3.6 [M107] | K∥ = K₀βk₁k₂, k₁ = B_eq/B, k₂ = (R/R_k)^a[1 + (R/R_k)^{(b−a)/c}]^c; prose: "a for rigidity below the turnover point R_k and b for rigidity above it" | K₀ "a constant with units 10²⁰ cm² s⁻¹"; c = 2.2 (fixed); K⊥,r = 0.02K∥; K⊥,θ = (2 + tanh[8\|θ − 90°\| − 280°]) × K⊥,r (Eq. 3.5, as printed) | Full text ×5; a, b, R_k and K₀ only as plotted points (Fig. 2); slope order depends on the sign of b − a (D-33) |
| GCR-QINSHEN17 | Qin & Shen 2017 [M108]; Shen & Qin 2018 [M109] | NLGCE-F λ∥ (PARALLEL rev. 1.4, Section 11) with δB = δB_1AU R^S(1 + sin²θ)/2, S = −1.56 + 0.09 ln(α/α_c), α_c = 1° | λ_slab = 0.02r; E_slab/E_total = 0.2; λ_slab/λ_2D = 10.0 (2017) or 2.6 (2018) | Full text ×1 (2017), ×2 (2018) |
| GCR-EB13 | Engelbrecht & Burger 2013a, Eq. 16; 2014, Eq. 9 | Eq. (21) | k_m = 1/λ_sl (as extracted); NLGC a² = 1/3; λ_out = 12.5λ_c,2D | Full text ×2 |

Notes on the table:

- **Wang et al. 2019.** The two branches meet at R = 0.1 GV when R is in GV, since (1/3)(0.1) = 1/30; this is arithmetic, not a statement of the paper.
- **Tomassetti 2017.** The final arXiv version (v5) prints "B₀ ≅ 3.4 nT sets the HMF at r₀ = 1 AU"; v1 printed "B₀ ≅ 3.4 nT AU²" (Section 14, D-12). Table II (12 years, January 2005–January 2017):

  | NM station | a (MV) | b | χ² |
  |---|---|---|---|
  | NEWK | 697.5 ± 4.5 | 0.39 ± 0.01 | 2434.12 |
  | OULU | 680.3 ± 4.8 | 0.40 ± 0.01 | 2393.75 |
  | APTY | 709.8 ± 4.6 | 0.42 ± 0.01 | 2335.46 |
  | JUNG | 689.9 ± 4.5 | 0.38 ± 0.01 | 2426.4 |

- **Perugia models.** Fiandrini et al. (2021) report, in v2, "average value of a = 1.21±0.06", b = 0.74 ± 0.03 during 2006–2009, a maximum average b of 1.3 ± 0.07, and "λk is found to range between 0.05 AU and 0.3 AU, depending on solar activity"; v1 reports different numbers. K₀(t), a(t) and b(t) are only in figures. h = 3 (2021) and h ≡ 0.01 (2025) give very different transitions at R_k, so parameters from one paper cannot be used with the other's h. The polar function g(θ) of Fiandrini et al. (Eq. 9) is printed with "θ_A + π/2"; the NWU original has θ_A − 90° (Section 14, D-13). Tomassetti et al. (2023, Eq. 3) [M110] print λ∥ = K₀(B₀/B)(R₀/R)^a[((R/R₀)^h + (R_k/R₀)^h)/(1 + (R_k/R₀)^h)]^{(b−a)/h}, with K∥ = βcλ∥/3 given separately, R₀ ≡ 1 GV and K₀ "in units of 10²³ cm² s⁻¹", which is dimensionally inconsistent for a length. The factor (R₀/R)^a is inverted relative to Fiandrini et al. (2021, Eq. 8); as printed it gives λ∥ ∝ R^{−a} below R_k, whereas the same paper says a and b "set the slopes of the rigidity dependence of λ∥ below and above Rk, respectively" (Section 14, D-23).
- **Duan et al. 2025.** With x = R/R_k and c > 0, the logarithmic slope of the printed k₂ is a + (b − a)y/(1 + y) with y = x^{(b−a)/c} (Eq. 47, derived). For b > a, y → 0 at x ≪ 1 and y → ∞ at x ≫ 1, so the low- and high-rigidity slopes are a and b, as the paper's prose says and as in Eq. (34). For b < a the limits are reversed: the printed k₂ has slope b below R_k and a above it, contrary to the prose, whereas Eq. (34) gives a below and b above for either sign of b − a. Whether the fitted sets of Duan et al. have b > a cannot be read from the text (values only in Fig. 2). The printed form is therefore reproduced as printed, under a variant name, and both orderings are tested (Section 14, D-33; fixtures F-GCR-08 and F-GCR-09). Even for b > a the transition shape differs from Eq. (34); c in Duan et al. is not the c of Eq. (33), and the two forms are not interchangeable.
- **Engelbrecht & Burger 2013a.** The extraction defines k_m = 1/λ_sl ("slab bendover"), whereas Lang et al. (2024, Eq. 13) relate k_min to λ_slab through Eq. (20). Both are published choices; a configuration must name which it uses.

### 9.6 Normalization conventions

| Model | 1/3 inside the normalization? | Normalization symbol: UNIT or VALUE | Rigidity reference | Field reference |
|---|---|---|---|---|
| NWU 2013–2023 (Eq. 33) | no (β, no 1/3) | (K∥)₀: UNIT 10²² cm² s⁻¹ (Potgieter 2013, 2014, 2015; Aslam 2019, 2023); UNIT 6 × 10²⁰ (Raath 2015; Aslam 2021) | P₀ = 1 GV | B₀ = 1 nT |
| Vos & Potgieter 2015–2016 | no | κ∥0 "in units of cm² s⁻¹", no value | P₀ = 1 GV (2016) | B₀ = 1 nT (2016) |
| Corti 2019 (Eq. 34) | no | k∥⁰: UNIT 6 × 10²⁰ cm² s⁻¹ | R_k = 4.3 GV | 1 nT |
| Luo 2019 | no | κ₀: UNIT 10²⁰ cm² s⁻¹ | P_k = 4.3 GV | 1 nT (B power illegible) |
| HelMod (Eq. 36) | yes (β/3) | K₀: VALUE in AU² GV⁻¹ s⁻¹ | 1 GV | none |
| Bobik 2012 | inside B⊕/(3B) | k₁, K₀: unit not stated | K_P ≈ P | B⊕ ≈ 5 nT |
| Strauss 2011 | λ given | λ₀: VALUE 0.15 AU | P₀ = 1 GV | none |
| Effenberger 2012 | no | κ∥0: VALUE 0.9 × 10²² cm² s⁻¹ | p₀ = 1 GeV/c | B_e at 1 AU |
| Wang 2019 | yes | k = k₀ · 3.6 × 10²² cm² s⁻¹, k₀ dimensionless | R in GV | B_E at Earth |
| Tomassetti 2017 | yes | κ⁰ dimensionless × 10²² cm² s⁻¹ | GV | B₀ ≅ 3.4 nT at r₀ = 1 AU (final version; v1 printed "nT AU²", D-12) |
| Perugia 2021–2025 | yes (β/3) | K₀: UNIT 10²³ cm² s⁻¹ (v1) / "order of" | R₀ = 1 GV, R_k = 3 GV | B₀ at Earth |
| Jiang 2023 | no | K₀: "order of" 10²³ cm² s⁻¹ | R₀ = 1 GV, R_k = 3 GV | B₀ at Earth (5 nT) |
| Duan 2025 | no | K₀: UNIT 10²⁰ cm² s⁻¹ | R_k (free) | B_eq near Earth |

**Implementation requirement.** Every GCR preset stores (i) whether 1/3 is included, (ii) the unit of the normalization constant as printed (UNIT) or its value (VALUE), (iii) the rigidity reference, and (iv) the field reference. A tabulated λ∥(P_ref) at Earth is stored as the primary published quantity when the paper gives no normalization value; converting it to (K∥)₀ uses λ∥ = 3K∥/v and the paper's own B at Earth, and is labelled derived.

Published perpendicular-to-parallel ratios used with these forms: 0.02 radial and 0.02 polar (×f(θ)) in Potgieter 2013–2014, Aslam 2019–2023, Fiandrini 2021 and Jiang 2023; 0.02 radial and 0.01 polar times a polar-enhancement function in Vos & Potgieter 2015 (h(θ), Eqs. 9–10), Vos & Potgieter 2016 (H⊥θ, Eqs. 8–9) and Corti 2019 (u(θ), Eqs. 10–11); 0.02 radial and 0.01 polar ("scaled to 2% and 1% of κ∥") in Di Felice et al. 2017 [M111]; 0.01 in Strauss 2011; 0.1 with ξ = 2 in Effenberger 2012; 0.05 with a polar factor 10 in Bobik 2012; 0.058, 0.06 and 0.065 (p) or 0.05 (e±) in HelMod 2016–2019; 0.02 in Wang 2019 and Tomassetti 2017; 0.02 radial and 0.02 × (2 + tanh[8|θ − 90°| − 280°]) polar in Duan 2025 (Eq. 3.5). From Vos & Potgieter (2015) onward the NWU perpendicular coefficient has its own high-rigidity slope (b⊥ = 1.58 versus b∥ = 2.10; b₂ = 1.24 versus b₁ = 1.55; c₂⊥ = 0.75c₂∥), so K⊥ is not always a constant multiple of K∥. The perpendicular specification covers these closures.

## 10. Coupling to perpendicular and radial transport

Perpendicular closures are specified in the PERPENDICULAR companion document. This section records only the forms that SEP codes pair with the λ presets of Section 5, because they change how a fitted λ∥ or λ_r must be interpreted.

**Dröge-type gyroradius scaling.** Dröge et al. (2016, ICRC, Eq. 3.4) use

$$
\Lambda_\perp(\mu)=\frac{3}{v}D_\perp=\alpha\,\lambda_\parallel(r)\left(\frac{r}{1\,\mathrm{AU}}\right)^2\cos\psi(r)\sqrt{1-\mu^2},
\tag{38}
$$

where "α is a parameter to model the relative contributions of parallel and perpendicular diffusion" (2010 Aug 7: α = 0.08 "assumed", and α = 0.02 "determined" for the sector surrounding STEREO-B; 2010 Aug 18: Λ⊥ = 0.17 · 0.06 · Σ(r,ψ,μ) as printed, with no separate α stated). Wijsen et al. (2019, Eq. 10) use, in PARADISE,

$$
\kappa_\perp=\frac{\pi}{12}\,\frac{\alpha\,v\,\lambda^r_\parallel}{b_r^2}\,\frac{B_0}{B},\qquad \alpha=10^{-4},
\tag{39}
$$

with λ^r∥/b_r² = λ∥. Averaging Eq. (38) over μ with κ⊥ = ½∫₋₁¹(v/3)Λ⊥ dμ and ⟨(1−μ²)^{1/2}⟩ = π/4 gives κ⊥ = (π/12)αvλ∥(r/1 AU)²cos ψ, which equals Eq. (39) when B₀/B = (r/1 AU)²cos ψ, i.e. for a Parker field with B₀ the radial field at 1 AU (derivation in this review; fixture F-PERP-01). Wijsen et al. write "Similarly to Dröge et al. (2010), we assume a perpendicular mean free path which scales with the gyro-radius of the particle", but choose B₀ = max_{r=1 AU}B = 9.7 nT, the maximum field strength at 1 AU (their Sect. 5.2), rather than the radial field. The two prescriptions therefore have the same form but different normalizations, and their α values are not interchangeable.

**Other pairings recorded in Section 5.**

| Code / paper | Perpendicular prescription | Values |
|---|---|---|
| He et al. 2011; He & Wan 2015 | constant λ_x = λ_y | 0.006 AU (50 MeV protons); 0.03 AU (40 keV electrons, He 2015) |
| Wang & Qin 2015 | κ⊥ ∝ (p/1 GeV c⁻¹)^{1/3}, set by ratio | κ⊥/κ∥ = 0.01 in the ecliptic at 1 AU |
| Wang et al. 2014 | NLGC ratio scan | κ⊥/κ∥ = 0, 10⁻⁵, 5 × 10⁻⁵, 10⁻⁴ |
| Qin & Wang 2015 | NLGC | (δB_2D)²/(δB_slab)² = 4, l_2D = 3.1 × 10⁻³ AU; fitted κ⊥ ∼ 1–3% of κ∥ at 1 AU |
| Strauss & Fichtner 2015 | λ⊥ = ηλ∥; D⊥ forms constant, FLRW 2D⊥,0\|μ\|, scattering (4/π)D⊥,0(1−μ²)^{1/2}; κ⊥ = ½∫D⊥dμ | η = 0.02, 0.1 |
| Strauss et al. 2017 [M62] | FLRW D⊥ = av\|μ\|κ_FL, λ⊥ = (3/2v)∫D⊥dμ | a = 0.2 (electrons); a = 1/√3 ≈ 0.6, "the generally used value … used mainly for protons" |
| Laitinen et al. 2016, 2018 | NLGC, a_NLGC = 1/√3, plus field-line random walk D_FL | — |
| Zhang et al. 2023 | κ⊥ = (v/2V)α⊥κ_gd0B₀/B (Eq. 38); "v/V is the ratio of particle to solar wind plasma speed", B₀ the field on the solar surface on the same field line | κ_gd0 = 3.4 × 10¹³ cm² s⁻¹; α⊥ = 0.37 (λ_r = 200 R☉) and 0.074 (λ_r = 20 R☉) |
| Lavasa et al. 2026 | D⊥(μ) = 2D⊥,0\|μ\| (FLRW), κ⊥ = ½∫D⊥dμ, λ⊥ = 3κ⊥/v | fitted (Table 3) |
| iPATH (Ding et al. 2020, 2022) | NLGC | κ⊥/κ∥ = 0.0099 (1 MeV) and 0.03 (10 MeV) at 1 AU |
| Kozarev & Schwadron 2016 | κ⊥ = κ∥/[1 + (λ∥/r_g)²] | analytic shock model |
| EPREM (Kozarev 2013), M-FLAMPA, SOFIE, PATH 2009, PARASOL, SEPMOD, Kubo 2015 | none | — |

**Published ratios.** Palmer (1982): "K⊥r / K∥ < 0.1 at 1 AU". Shalchi et al. (2006): 0.02 ≤ λ⊥/λ∥ ≤ 0.083. Tautz & Shalchi (2013), test particles with δB/B₀ = 0.5: average λ⊥/λ∥ = 0.0379 ± 0.0093 (electrons) and 0.0597 ± 0.0426 (protons). He & Wan (2013): 0.001–0.2 across events. Dröge et al. (2016, ApJ 826, 134, abstract) [M17]: λ⊥ = 0.005–0.01 au with λ∥ = 0.15–0.6 au. Lavasa et al. (2026): 0.9–3% (electrons), 4.5–10% (protons). Strauss et al. (2017): 0.001–0.03 for a = 1/10.

**Implementation requirement.** A λ_r preset combined with a nonzero perpendicular coefficient does not return λ_rr of Eq. (3) unless λ⊥ sin²ψ is added explicitly. The library reports λ∥, λ⊥, λ_r and λ_rr separately and never relabels one as another.

## 11. Turbulence inputs

The closed forms of Section 7 and the QLT shapes of Section 6.8 need turbulence quantities. This section records how each quantity is defined in the literature and which values have been published, so that an input can be traced. Published values are not defaults: a configuration supplies its own inputs and states their source. Transcriptions are in `turbulence/published_turbulence_values.csv` (Section 19).

### 11.1 Inputs required by each closure

| Closure | Required inputs |
|---|---|
| TS2003 / NWU (Eqs. 21–23) | B₀; δB²_slab; k_min (or λ_slab with Eq. 20, or k_m = 1/λ_sl); s; for electrons p, k_d, V_A, α_D and the DT/RS choice |
| Zank fit (Eq. 24) | B; per-component or total slab variance (state which); slab correlation length ℓ_slab or λ_s |
| QLT with supplied spectrum (Section 6.8) | the spectrum and its normalization convention |
| Qin & Wang amplitude (Eq. 15) | (δB_slab/B₀)²; k_min; s; h |
| All | particle rigidity P, speed v, and for λ_r the field angle ψ |

### 11.2 Slab fraction

"Slab fraction" has three distinct definitions in this literature, and they are not the same quantity:

1. **Band-limited magnetic-energy fraction** from ratio or anisotropy tests on spacecraft data, written with spectral amplitudes (C_s, C_2 or A_S, A_2D): MacBride et al. 2010 [M112], r = C_s/(C_s + C_2); Fa & He 2026 [M113], r = A_S/(A_S + A_2D); Zhao et al. 2022 [M114]; Cheng et al. 2025 [M115].
2. **Magnetic-variance ratio** ⟨b_s²⟩/⟨b_2²⟩: Chhiber et al. 2017, Eq. (20), = 20/80 = 0.25 (the left side is typeset ⟨b_s²⟩/⟨b_s²⟩); Zhao et al. 2018 [M116], 80:20 and 60:40 applied to ⟨b²⟩.
3. **Elsässer-energy partition** Z²:W² or E_2D:E_slab, which includes kinetic energy: Oughton et al. 2011 [M117]; Engelbrecht & Burger 2013a; Adhikari et al. 2017, 2020 [M118, M119].

Under a constant Alfvén ratio applied identically to both components (Engelbrecht & Burger 2013a; Wiengarten et al. 2016, Eq. 36 [M120]) definitions 2 and 3 give the same ratio; otherwise they do not (algebraic consequence of the cited formulas).

**Published values.**

| Source | Location / data | Definition | Value |
|---|---|---|---|
| Bieber, Wanner & Matthaeus 1996 [M121] | Helios; abstract only | 1 | "a dominant 2D (∼85% by energy) component" |
| MacBride et al. 2010 | Helios 1, 0.3–1 au | 1 | 27 ± 12% slab (all data); Table 3 fits 0.27 (all), 0.33 (fast), 0.16 (slow); averages 0.25 ± 0.02, 0.30 ± 0.05, 0.25 ± 0.04; "no statistically significant spatial evolution" |
| Leamon et al. 2000 [M122] | 1 au | 1 | inertial range ≈ 80% 2D / 20% slab (citing Bieber 1996); dissipation range "approximately 50% two-dimensional and 50% slab energy" |
| Fa & He 2026 | Wind 1995–2023, 0.0026–0.0547 Hz | 1 | r_ratio = 0.27, r_aniso = 0.13, mean 0.20; yearly r_ratio 0.17 (2018)–0.39 (1999), r_aniso 0.00 (1996)–0.24 (2002); r_ratio = 0.19 S_N^{0.1}, r_aniso = 0.07 S_N^{0.2} |
| Zhao et al. 2022 | PSP orbits 1–7, 0.01–0.1 Hz | 1 | C₂/C_s within 0.3 au from 0.11 to 0.75 (slab 57–90%); 0.3–0.6 au mostly 1.50–1.70 (slab 37–40%), 0.67 in orbit 7; "about 60%–80%" slab within 0.3 au |
| Cheng et al. 2025 | PSP encounters 1–19, < 0.3 au | 1 | C₂/C_s = 0.35 (26%:74%) coronal-hole wind; 0.83 (45%:55%) streamer wind |
| Subashchandar et al. 2025 | PSP | spectral fit | slab inertial-range spectra "approximately one to two orders of magnitude (i.e., 10–100 times) higher in power than the corresponding 2D turbulence spectrum" |
| Chhiber et al. 2017 | model | 2 | slab/2D = 0.25 |
| Zhao et al. 2018 | model on OMNI | 2 | 80:20 and 60:40 |
| Oughton et al. 2011; Wiengarten et al. 2016 | model boundary | 3 | Z² = 1500, W² = 150 km² s⁻² (≈ 91:9); Wiengarten: "90%-10% partitioning for Z²-W²" |
| Zank et al. 2017, 2020 [M123, M124]; Adhikari et al. 2017, 2020 | model | 3 | 80:20 |

Near the Sun, Zhao et al. (2022), Cheng et al. (2025) and Subashchandar et al. (2025) all find slab-dominated fluctuations, with different methods (ratio tests versus a fitted composite spectrum). At 1 au, the long-term means are 0.13–0.27 (Fa & He 2026) and the yearly values 0.00–0.39; MacBride et al. (2010) give 0.16–0.33 by wind type. A slab fraction of 0.2 is therefore a convention, not a measurement.

### 11.3 Fluctuation variance and its radial dependence

No directly fitted observational law δB² ∝ r^{−3.5} was found in a primary source. Chhiber (2022) [M125] uses the asymptotic exponents of Zank, Matthaeus & Smith (1996) [M126]: with δb² ∝ r^{−α}, α = 3 for WKB without dissipation or mixing, α = 3.5 for "dissipation … without mixing", and α = 4 with mixing and dissipation; for a radial mean field ∝ r⁻², δb/B₀ ∝ r^{2−α/2}. Chhiber (2022) finds δb/B₀ "well approximated by a r^{1/4} power law", "with a value ∼0.8 near 1 AU and ∼0.5 at 0.14 AU", which α = 3.5 matches. The same paper quotes turbulence energy "from ∼3 nT² at 200 R☉ to ∼50 nT² at 30 R☉" and B₀ "∼60 nT at 30 R☉"; with δb/B₀ ≈ 0.5 these magnitudes are not mutually consistent under a direct reading (0.5² × 60² = 900 nT²), so neither is used as a normalization here (Section 14, D-14).

Model outputs: Adhikari et al. (2020), 0.17–1 au, quasi-2D ∝ r^{−3.12}, slab ∝ r^{−2.97}, total ∝ r^{−3.1} (Conclusions: r^{−3.14}, r^{−2.99}, r^{−3.1}); Adhikari et al. (2017), 2D ⟨b²⟩ ∝ r^{−2.91} (1–10 au) and r^{−2.30} (10–75 au). Observed 1 au levels vary with the solar cycle: Zhao et al. (2018) find total ⟨b²⟩ in 2003 "almost four times larger than in 2009" and "1.5 times larger than in 2015", with Table 1 (inward IMF) ⟨z+²⟩ = 1500.31, 363.99, 635.72 km² s⁻², ⟨z−²⟩ = 629.98, 198.01, 326.57 km² s⁻², E_D = −387.01, −110.83, −208.52 km² s⁻² and |B| = 8.61, 4.71, 7.65 nT for 2003, 2009 and 2015. Smith et al. (2001, Sec. 3.1) [M127] write of the 1 au Omnitape data that Z² values of "200 to 400 km² s⁻² are most typical". Subashchandar et al. (2025) report δB/B₀ ≈ 0.2–0.6 for PSP.

**Conversion to δB².** Breech et al. 2008 [M128], Oughton et al. 2011, Wiengarten et al. 2016, Usmanov et al. 2018 [M129] and Chhiber et al. 2017 use Z² = ⟨v² + b²⟩ = ⟨z+² + z−²⟩/2 ("twice the energy per unit mass"), with b in Alfvén units. With the Alfvén ratio r_A = ⟨v²⟩/⟨b²⟩,

$$
\delta B^2=\frac{\mu_0\rho Z^2}{r_A+1}\ \ (\mathrm{SI}),\qquad \langle B'^2\rangle=\frac{4\pi\rho Z^2}{r_A+1}\ \ (\mathrm{Gaussian}),
\tag{40}
$$

(Wiengarten et al. 2016, Eq. 36, applied to Z² for 2D and W² for slab; Engelbrecht & Burger 2013a, Eq. 15; Chhiber et al. 2017, Eq. 19). r_A = 1/2 corresponds to σ_D = (r_A − 1)/(r_A + 1) = −1/3. The only readable rendering of Breech et al. (2008) shows "σ_D = 1/3" next to r_A = 1/2, but that rendering drops minus signs throughout, so the printed sign is not established; their own relation r_A = (1 + σ_D)/(1 − σ_D) and Usmanov et al. (2018, Table 1) give −1/3 (Section 14, D-15). Wiengarten et al. (2016) print the SI factor μ₀ρ in Eq. 36 while normalizing their Elsässer variables with √(4πρ); a configuration uses one unit system throughout. Chhiber et al. (2017) use r_A = 1 for 1–45 R☉ and 1/2 beyond. NI-MHD outputs (Zank et al. 2017; Adhikari et al. 2017, 2020) are per-component ⟨z±²⟩ and E_D; their conversion to ⟨b²⟩ was not transcribable here, so a data set using them must state its conversion.

### 11.4 Radially evolving prescriptions used by SEP codes

- **Laitinen et al. (2016)** (SEP-LAITINEN16): δB²/B² = 0.04 at 1 AU; slab:2D = 20%:80%; breakpoints L∥ = L⊥ = 0.007 AU; outer scale L₀(r) = r; radial evolution (Eqs. 13–14)

  $$
  W(r)=\left(\frac{r_0}{r}\right)^3\left(\frac{V_{sw,r_0}+v_{A0}}{V_{sw,r_0}+\frac{r_0}{r}v_{A0}}\right)^2,
  \qquad V_{sw,r_0}=400\ \mathrm{km/s},\ v_{A0}=30\ \mathrm{km/s},\ r_0=1\ \mathrm{AU},
  \tag{41}
  $$

  with the amplitude normalized so that λ∥ = 0.3 AU for a 10 MeV proton at 1 AU. The resulting particles propagate "essentially without scattering from the Sun to 0.4 AU".
- **Strauss et al. (2017)** [M62] (Table 1): δB²_slab = 0.2δB²; s = 5/3; p = 2.6; k_min = 35 AU⁻¹; δB²(1 AU) = 13.2 nT² ∝ r^{−2.4}; 2D spectral index q = 7. This is the only published source in this review for a 1 AU value of 13.2 nT².
- **Ding et al. (2020, 2022)** (iPATH; Section 8.5): δB² ∼ r⁻³ or r^{−3.5}; l_slab ∼ l_2D ∼ r^{0.8} or r^{1.0}; δB²/B² = 0.5 or 0.15 at 1 au; l_slab = 10⁹ m.
- **Qin & Shen (2017), Shen & Qin (2018)** (GCR-QINSHEN17): δB = δB_1AU R^S(1 + sin²θ)/2 with S = −1.56 + 0.09 ln(α/α_c); λ_slab = 0.02r.
- **Effenberger et al. (2012)** print an alternative QLT λ∥ with l_slab = 0.03ρ^{0.5} and δB²_slab = B_e²ρ^{−2.15}, from a garbled extraction; it is not implemented.

### 11.5 Correlation lengths

Published "correlation lengths" use incompatible definitions:

| Definition | Source |
|---|---|
| e-folding: R_C(τ_e) = 1/e, λ_C = V_SWτ_e | Cuesta et al. 2022 [M130]; Adhikari et al. 2017 |
| integral scale: ∫₀^∞⟨b(0)·b(r)⟩/⟨b²⟩dr | Ruiz et al. 2014 [M131] |
| larger decay constant of C = β e^{−r/λ_CS} + (1 − β)e^{−r/λ₂} | Weygand et al. 2011 [M132] |
| spectral bendover: k_m = 1/λ_sl | Engelbrecht & Burger 2013a; Wiengarten et al. 2016 |
| bendover-spectrum correlation length λ_s = 2πC(ν)l | Zank fit (Section 7.5) |
| f_mid, 50% of fluctuation power | Subashchandar et al. 2025 |

The integral and e-folding scales coincide only for an exponential correlation function. Published values:

- **Cuesta et al. (2022), Table 1 (ACE/WIND, 10⁶ km):** all wind λ_C^∥ = 1.98 ± 0.09, λ_C^⊥ = 1.55 ± 0.02 (ratio 1.28 ± 0.06); V_SW < 450 km/s 1.98 ± 0.11 and 1.67 ± 0.02 (1.19 ± 0.07); V_SW > 600 km/s 1.78 ± 0.31 and 1.05 ± 0.05 (1.70 ± 0.31). Table 2, λ_C ∝ R^α: α = 0.97 ± 0.04 (R < 0.30 au), 0.29 ± 0.01 (0.30–1.0 au), 0.27 ± 0.01 (R > 1 au); for λ_C^∥ 1.03 ± 0.06, 0.28 ± 0.03, 0.64 ± 0.12; for λ_C^⊥ 0.61 ± 0.06, 0.27 ± 0.02, 0.23 ± 0.01.
- **Ruiz et al. (2014):** λ(D) = 0.89(D/1 AU)^{0.43} × 10⁶ km over 0.3–5 au (no exponent error reported); Table 1 mean/median 0.85/0.67 (Helios, 0.3–0.7 au), 1.12/0.97 (ACE), 2.03/1.85 (Ulysses, 3–5.3 au), in 10⁶ km.
- **Weygand et al. (2011), Table 1 (10⁶ km):** slow wind parallel 2.8 ± 0.8, perpendicular 1.1 ± 0.1 (ratio 2.55 ± 0.76); fast wind parallel 1.0 ± 0.2, perpendicular 1.4 ± 0.5 (ratio 0.71 ± 0.29).
- **Zhao et al. (2018):** average λ_s ≈ 0.88 × 10⁶ km (OMNI); λ_s in 2003 twice that in 2009.
- **Lang et al. (2024), Table 5:** λ_slab = (2.81 ± 1.51) × 10⁶ km, λ_2D = (1.10 ± 0.49) × 10⁶ km at 1 au.
- **Model conventions:** λ_s = 2λ_2D (Chhiber et al. 2017; Zhao et al. 2018); in Chhiber et al. λ ∝ B^{−1/2} in the inner region; Adhikari et al. (2020) model correlation lengths ∝ r^{1.12} (2D) and r^{1.11} (slab).

### 11.6 Spectral indices and dissipation range

- Inertial range: α_B ≈ −3/2 at 0.17 au to ≈ −5/3 at 0.6 au for 10⁻²–10⁻¹ Hz (Chen et al. 2020 [M133]); q = 1.55 ± 0.07 (coronal-hole wind) and 1.59 ± 0.10 (streamer wind) within 0.3 au (Cheng et al. 2025); Subashchandar et al. (2025) assume a slab index 1.5 near the Sun. Theory (Zank et al. 2020): 2D Kolmogorov in k⊥; slab −5/3, −3/2 or −2 depending on the time-scale regime.
- Dissipation-range onset at 1 au: "a few tenths of a hertz" (Leamon et al. 1998, abstract) and "≥0.3 Hz" (Smith et al. 2006, abstract). Leamon et al. (2000, Table 1) fit the break as a + b(X/2π): for X = Ω_ci, a = 0.200, b = 1.760, χ² = 2.93; for X = k_res·V_sw, a = 0.274, b = 0.360; for X = k_ii·V_sw, a = 0.152, b = 0.451, the last "only marginally better than the gyrofrequency fit" (OCR-quality extraction). The Strauss et al. (2020) relation, Eq. (5) in this document, uses the same a and b; it agrees with this fit if its Ω_i denotes Ω_ci/2π.
- Ion-scale index ≈ −2.8 (Alexandrova et al. 2009, abstract) [M134]. Smith et al. (2006) report "a broad range of power-law indexes" in the dissipation range, steeper for greater cascade rates; their numbers were not accessible.
- The dissipation-range index p and onset k_d are therefore explicit inputs with no universal value; the published fits used with Eq. (22) are in Section 7.4.

### 11.7 Turbulence-transport model boundary values

| Model | Location | Values as printed |
|---|---|---|
| Breech et al. 2008 | 0.3 AU, ecliptic (Fig. 4, Case 1) | Z² = 3000 (km/s)², λ = 0.008 AU, σ_c = 0.60, T = 1.8 × 10⁵ K |
| Breech et al. 2008 | 0.3 AU, θ = 75° (Fig. 5, Case 1) | Z² = 6000 (km/s)², λ = 0.03 AU, σ_c = 0.80, T = 8.5 × 10⁵ K |
| Oughton et al. 2011 | 0.3 AU (Fig. 1) | Z² = 1500, W² = 150 km² s⁻², ℓ = l = 0.008 AU, l∥ = 0.036 AU, σ_c = σ̃_c = 0.6, T printed "1.6 × 10⁶ K" |
| Oughton et al. 2011 | 1 AU (§3.2) | ℓ = l = 0.014 AU, l∥ = 0.037 AU, σ_c = σ̃_c = 0.5 |
| Wiengarten et al. 2016 | 0.3 au, validation | Z² = 1500, W² = 150 km² s⁻², σ_c = 0.6, ℓ = λ = 0.008 AU, λ∥ = 0.036 AU, T = 1.6 × 10⁵ K |
| Wiengarten et al. 2016 | 0.3 au, low latitude | Z² = 900, W² = 90, σ_c = 0.4, ℓ = λ = 0.012 AU, λ∥ = 0.03 AU, T = 3.0 × 10⁵ K |
| Wiengarten et al. 2016 | 0.3 au, high latitude | Z² = 5000, W² = 500, σ_c = 0.6, ℓ = λ = 0.018 AU, λ∥ = 0.03 AU, T = 1.5 × 10⁶ K |
| Engelbrecht & Burger 2013a | 0.3 AU, ecliptic / polar | Z² = 1250 / 1600 km² s⁻²; W² = 350 / 3000; λ = l = 0.004 / 0.015 AU; λ_c,s = 0.011 / 0.011 AU; σ_c = σ̃_c = 0.6 / 0.8; T = 2 × 10⁵ / 1.6 × 10⁶ K (Table 1); B_E = 5 nT, n = 7 cm⁻³ (text) |
| Smith et al. 2001 | 1 AU (Sec. 3.2; alternates Sec. 3.2 and App. B) | Z² = 350 km² s⁻² (alternates 250, 400, 650), Λ = 0.03 AU, T = 60 000 K |
| Adhikari et al. 2017 | 1 au | 2D ⟨z∞+²⟩ = 4000, ⟨z∞−²⟩ = 800, E_D∞ = −100 km² s⁻²; slab ⟨z*+²⟩ = 1000, ⟨z*−²⟩ = 200, E_D* = −25 km² s⁻²; T = 8 × 10⁴ K |
| Adhikari et al. 2020 | 0.165 au | quasi-2D ⟨z⁺²⟩ = 9338.4, ⟨z⁻²⟩ = 952.4, E_D = −112.48 km² s⁻²; slab 2334.6, 238.1, −28.12 km² s⁻²; T = 1.75 × 10⁵ K |

The Oughton et al. (2011) temperature "1.6 × 10⁶ K" (Fig. 1 caption, re-read in a copy of the published paper) conflicts with the 1.6 × 10⁵ K listed by Wiengarten et al. (2016) for the same validation case, whose other boundary values are identical (Section 14, D-16). Correlation-length quantities L (km³ s⁻²) of the NI-MHD models are in the CSV.

### 11.8 Mean field

Raath et al. (2016, Eqs. 2–4) [M135] give the Parker field used by the NWU codes:

$$
\mathbf B=B_0\left(\frac{r_0}{r}\right)^2(\mathbf e_r-\tan\psi\,\mathbf e_\phi),\qquad
\tan\psi=\frac{\Omega(r-r_\odot)\sin\theta}{V_{sw}},\qquad
B=B_0\left(\frac{r_0}{r}\right)^2\sqrt{1+\tan^2\psi},
\tag{42}
$$

with r₀ = 1 AU, r☉ = 0.005 AU, Ω = 2.66 × 10⁻⁶ rad s⁻¹, V_sw ≈ 400 km/s (equator) and 800 km/s (poles), and field values 5.05 nT (2006), 4.50, 4.25, 3.94 nT (2007–2009), labelled B_e in their Table 1. Raath et al. describe B₀ as "the magnitude of the HMF at r₀ = 1 AU (i.e. at Earth)", but in Eq. (42) B₀ is the amplitude of the radial component at r₀: the magnitude there is B₀√(1 + tan²ψ(r₀)) = |B_r(r₀)|/|cos ψ(r₀)|, which is 1.4071B₀ in the ecliptic for V_sw = 400 km/s (derived). Inserting a published field magnitude at Earth as B₀, or a radial amplitude as the magnitude, changes the field by that factor; the choice must be configured explicitly (Section 14, D-34 and U-13). They also give the Smith–Bieber modification (Eq. 31, b = 20r☉, B_T(b)/B_R(b) ≈ −0.02) and the Jokipii–Kóta modification (Eqs. 29–30, δ_m = 8.7 × 10⁻⁵, δ = 0.002 near the poles); Parker (1958) and Smith & Bieber (1991) are secondary here. Other normalizations in use: B₀ = 5 nT (Engelbrecht & Burger 2013a), 5.5 nT (Zhao et al. 2018 example), and fields taken directly from MHD (Chhiber et al. 2017) or OMNI (Zhao et al. 2018, tan ψ = B_Y/B_X).

## 12. Library interface additions

The PARALLEL rev. 1.4 interface (its Section 14) evaluates κ∥ from a model and a local state. The mean-free-path layer specified here sits on top of it and adds four things: tagged λ quantities, presets with explicit normalization, pitch-angle shapes with exact amplitude conversion, and published parameter sets loaded from the companion data.

```cpp
enum class LambdaKind { Parallel, RadialSEP, RadialTensor, IsotropicScattering, Unspecified };
enum class MomentumVariable { Rigidity_V, MomentumPc_eV, KineticTotal_eV, KineticPerNucleon_eV };
enum class NormalizationKind { Value, UnitOnly };        // Section 9.6
enum class RuntimeState { ReadyExplicitInputs, ReadyPublishedVariant, ReadyPublishedKappaOnly,
                          RequiresExternalInput, RequiresUserDecision, RequiresSourceOrCodeAudit,
                          ReferenceDataOnly };            // Section 21.1, evaluated per configured run

struct LambdaValue {             // every returned or supplied length is tagged
  double metres;
  LambdaKind kind;
  std::optional<double> psi_rad; // required to convert RadialSEP <-> Parallel
};

struct PowerLawPreset {          // Eq. (9)
  std::string id;                // e.g. "SEP-MFLAMPA25"
  LambdaKind kind;
  double lambda0_m;
  MomentumVariable x_var; double x0;      // e.g. pc, 1 GeV
  double a;                      // momentum exponent
  double r0_m; double b;         // radial exponent
  std::string source_location;   // e.g. "Liu et al. 2025, Eq. (15)"
};

enum class PitchShape { QForm, QFormLangPrinted, EpsForm, EpsFormPachecoPrinted,
                        DrogeVA, Isotropic, Kolmogorov23, AMPS_I, AMPS_II, AMPS_III, AMPS_IV, AMPS_V };
enum class OperatorConvention { Standard /* d_mu(D d_mu f) */, HalfD /* d_mu((D/2) d_mu f) */ };

struct PitchAngleModel {
  PitchShape shape;
  double q = 5.0/3.0, H = 0.0, eps = 0.0, VA_over_v = 0.0;
  OperatorConvention op = OperatorConvention::Standard;
};

// D_mumu amplitude from a tagged lambda: exact Eq. (11) or closed forms of Section 6
double amplitude_from_lambda(const PitchAngleModel&, const LambdaValue&, double v);
LambdaValue lambda_from_amplitude(const PitchAngleModel&, double amplitude, double v);

enum class TSDissipationVariant { None, RS, DT_TS_Q3MS, DT_EB13, DT_LANG24 };   // Section 7.2, U-7
struct TSInputs { double B0_T, dB2_slab_T2, s, p, kmin_per_m, kd_per_m, VA_m_s, alphaD; TSDissipationVariant variant; };
LambdaValue lambda_ts2003(const TSInputs&, double rigidity_V, double v);         // Eqs. (21)-(23)

enum class ZankVariance { PerComponent, Total };                                 // U-12
LambdaValue lambda_zank98(double B_T, double dB2_slab_T2, ZankVariance, double lslab_corr_m, double rigidity_V);
```

**Implementation requirements.**

1. `LambdaKind::Unspecified` values (SEP-PATH09, SEP-SOLPENCO05) cannot be converted to κ∥ without a configured interpretation (U-11). The library refuses rather than guessing.
2. A GCR preset carries `NormalizationKind`, the printed unit, whether 1/3 is included, and the rigidity and field references (Section 9.6). Loading a published (K∥)₀ with a different unit from the preset is an error.
3. `QFormLangPrinted`, `EpsFormPachecoPrinted` and every DT variant exist only to reproduce papers as printed; they are never defaults.
4. Every published parameter set is loaded from the companion data by source key and location; the library does not hard-code transcribed numbers.
5. Derived quantities (for example an implied (K∥)₀ from a tabulated λ∥) are returned with a `derived` flag.

## 13. Verification fixtures

All fixtures are in `benchmark_points.json` (34 groups: the 32 groups of version 1.0, whose values are unchanged, plus F-GCR-09 and F-NUM-01) and are recomputed by `reference_verification.py` (Section 19). They are computed with 40-digit arithmetic and reported to 15 significant digits. A float64 implementation should reproduce algebraic fixtures to a relative 10⁻¹² and quadrature fixtures to a relative 10⁻⁹; tolerances are relative to the reference value, except that reference values that are zero or are themselves differences at the 40-digit level (keys `relative_difference`, `rel_diff`) are checked against absolute bounds. Values given with fewer digits are checked with the relative tolerances 10⁻⁴ (slopes from numerical differentiation, truncation-term estimates) and 10⁻⁷ (extremum locations). Values obtained by numerical quadrature (F-PA-05 numeric, F-QLT-07, the F-NUM-01 quadrature entries) are accurate to about 5 × 10⁻¹⁵, i.e. in the last printed digit.

| ID | Checks | Key values |
|---|---|---|
| F-KIN-01 | kinematics | 1 MeV proton 𝓡 = 43.3306378480744 MV; pc = 1 GeV proton T = 432.988102839397 MeV; alpha with pc = 1 GeV has 𝓡 = 0.5 GV; electron 0.094 MeV → 0.323889 MV, 0.94 MeV → 1.358042 MV, 2.0 MeV → 2.458454 MV, 2.5 MeV → 2.967321 MV |
| F-SEP-01 | presets at 1 au, pc = 1 GeV | 200 R☉ = 0.930093452192432 au |
| F-SEP-02 | λ_r ↔ λ∥ | factor 2 at ψ = 45°; implied cos²ψ 0.447761 (0.3/0.67), 0.441176 (0.15/0.34), 0.5 (He & Wan 2019); Parker example (tan ψ = Ω(r − r☉)/V with Ω = 2.66 × 10⁻⁶ rad/s, V = 400 km/s, r = 1 AU, r☉ = 0.005 AU, θ = 90°) ψ = 44.7078°, cos²ψ = 0.505100 |
| F-SEP-03 | Liu et al. Eq. (16) = vλ∥/3 | relative difference ≤ 10⁻⁴⁰ at 0.1 MeV–1 GeV; λ∥(1 au) = 0.0717824, 0.105371, 0.154786, 0.228967, 0.357767 au at 0.1, 1, 10, 100, 1000 MeV |
| F-SEP-05 | PARASOL | (43.3306/43)^{1/3} = 1.00255654 |
| F-SEP-06 | Chen et al. 2024 | κ∥ grid r = 0.1–0.8 au, E = 0.1–1000 MeV; derived proton λ∥ = 3κ∥/v (species choice U-2), e.g. 0.00420445 au (0.1 au, 100 keV), 0.0448555 au (0.5 au, 1 MeV), 0.126981 au (0.8 au, 10 MeV), 0.220292 au (0.8 au, 100 MeV) |
| F-PA-01 | I(5/3, H) | H = 0: 18/7; 0.01: 2.135993; 0.05: 1.693041; 0.1: 1.423601; 0.2: 1.119844 |
| F-PA-02 | φ(ε) | ε = 0.048: 4.597262 (asymptote 4.554831); ε = 0.01: 7.069533 (asymptote 6.907755) |
| F-PA-03 | Dröge V_A form | λ_integral/λ_nominal = 0.912870, 0.812321, 0.598512, 0.436445 at V_A/v = 0.001, 0.01, 0.1, 0.3 |
| F-PA-04 | EPREM operator | effective λ/λ∥ = 2 |
| F-PA-05 | M-FLAMPA shape | λ_xx/λ_μμ = 27/14; 81/(7π) = 3.683300; derived approximate coefficients 1.081726 and 0.560895 versus printed 0.9 and 0.5 |
| F-PA-06 | He & Wan focusing | λ/λ₀ = 0.999960, 0.996016, 0.715218, 0.222772, 0.0270000 at λ₀/L = 0.01, 0.1, 1, 3, 10 |
| F-QLT-01 | TS coefficients | 0.0106514 AU (TS2002 mid), 0.0177523 AU (TS2003 mid), 2.61511 × 10⁻¹⁰ AU (high), R = 1 at 12 351.45 MV |
| F-QLT-02 | Zank constants and accuracy | C(5/6) = 0.118862354635443; 2πC = 0.746834200222187; 5.16465750496608; 6.27420534071733 (total); 3.13710267035867 (per component); bracket/exact = 1.0000000, 1.0000004, 1.0000274, 1.0017973, 1.0156918, 1.0181874, 1.0014405 at R_L/l = 0.01, 0.1, 0.3, 1, 3, 10, 100; maximum 1.0216122 at R_L/l = 5.7868 |
| F-QLT-03 | Eq. (21), Lang Table 6 protons | 0.494342, 0.740432, 1.118621 au at 1, 10, 100 MeV |
| F-QLT-04 | Eq. (22), Lang Table 6 | e.g. 0.04 MeV: DT-LANG24 0.436753, DT-TS 0.538683, DT-EB13 0.473401, RS 0.639548 au; 0.52 MeV: 0.148798, 0.380851, 0.210156, 0.447411 au |
| F-QLT-05 | Lang Eq. 13 with Table 5 | λ_slab = 2.813926 × 10⁶ km (printed 2.81 × 10⁶ km) |
| F-QLT-07 | Eq. (21) versus exact slab QLT (TS2003 spectrum, s = 5/3) | exact/form = 0.997910, 0.911392, 0.793832, 0.900401 at R = 0.1, 1, 3, 10; minimum 0.793821 at R = 3.0284 |
| F-QLT-06 | identities | duplication identity and f₁(Lang) = f₁(TS)/π at p = 2.61, 3, 3.57 |
| F-SH-01 | Bohm | 1, 10, 100 MeV protons in 5 nT: r_g = 2.89071 × 10⁷, 9.16311 × 10⁷, 2.96594 × 10⁸ m |
| F-SH-02 | Afanasiev steady state | ratio exactly 1 |
| F-SH-03 | Subashchandar inline relations | κ∥ = 2.99895 × 10¹⁸ and 6.92065 × 10¹⁸ cm² s⁻¹ |
| F-SH-04 | M-FLAMPA floor | D_min = 6.957 × 10¹² m² s⁻¹ with R☉ = 6.957 × 10⁸ m |
| F-GCR-01 | NWU G(P) | G(1 GV) = 1; limiting slopes a and b |
| F-GCR-02 | Potgieter 2014 | implied λ∥(1 GV) = 0.146277, 0.182578, 0.222710, 0.230366 au for 2006–2009 |
| F-GCR-03 | Corti ↔ NWU | k∥⁰/(K∥)₀ = 3.16037378 for a = 0.8, b = 1.7, s = 2.2, R_k = 4.3 GV |
| F-GCR-04 | Vos & Potgieter 2015 | implied (K∥)₀ = 54.9078, 57.8794, 59.5989, 70.3766 × 10²² cm² s⁻¹ for 2006e–2009e with B_e = 4.95, 4.36, 4.11, 3.91 nT (Table 1); derived and conditional on B₀ = 1 nT, P₀ = 1 GV taken from the 2016 paper (D-19) |
| F-GCR-05 | HelMod | 1 AU² s⁻¹ = 2.23795229 × 10²⁶ cm² s⁻¹; c₀(A<0 asc) = 6.84590 × 10²² cm² s⁻¹ GV⁻¹ |
| F-GCR-06 | Tomassetti κ⁰ | NEWK: 2.13375, 1.55250, 1.08750 at φ = 400, 600, 1000 MV |
| F-GCR-07 | Strauss 2011; Wang 2019 | Strauss λ∥ = 0.3 au at (≤ 1 GV, 1 AU); Wang branches equal 1/30 at 0.1 GV |
| F-GCR-08 | Duan versus Corti, b > a | same limiting slopes; values at R_k differ (4.594793 versus 1.327849 for a = 0.8, b = 1.7, c = s = 2.2, R_k = 4.3 GV) |
| F-GCR-09 | Duan versus Corti, b < a (D-33) | a = 1.7, b = 0.8, c = s = 2.2, R_k = 4.3 GV (illustrative): Duan slopes 0.80001 (R = 10⁻¹² GV) and 1.7 (10¹⁴ GV), i.e. b below R_k and a above; Corti 1.7 and 0.8; Eq. (47) agrees with numerical differentiation |
| F-NUM-01 | stable numerics and closed forms (Section 21.4) | He–Wan λ/λ₀ = 0.9999999999996 (x = 10⁻⁶), 0.999960001618982 (x = 0.01); Zank q = 2.0 (s = 10⁻⁸), 1.99998333449065 (s = 0.01); ₂F₁ form of the Zank exact integral equals the quadrature (1.06243852893212 at x = 1); TS2003 closed form exact/form = 72/79 = 0.911392405063291 at R = 1, minimum 0.793820697646137 at R = 3.0284426 |
| F-PERP-01 | Dröge/PARADISE identity | ⟨(1−μ²)^{1/2}⟩ = π/4 |

The observational CSV also carries the kinematic consistency of the Lavasa et al. (2026) energy–rigidity pairs (F-KIN-01).

## 14. Discrepancy, decision and gap register

### 14.1 Printed discrepancies and inconsistencies

Each item is transcribed as printed and is not resolved by assumption.

| ID | Item | Where | Handling |
|---|---|---|---|
| D-1 | DT electron dissipation term: Q^{3−s} (TS2002/2003) versus Q^{p−s} (EB2013b, Lang 2024); ₂F₁ argument with Q (TS, Lang) versus Q^{p−2} (EB2013b) | Section 7.2 | variant flag, no default (U-7); identical at p = 3 |
| D-2 | EPREM D_μμ with operator ∂_μ[(D_μμ/2)∂_μf] implies effective λ = 2λ∥ | Section 6.6 | `OperatorConvention`; U-8 |
| D-3 | Lang et al. 2024 and Lavasa et al. 2026 print (1−μ)² in D_μμ | Section 6.2 | standard (1−μ²) implemented; printed form only as variant (U-5) |
| D-4 | Pacheco et al. 2019 print \|μ\|/(1−\|μ\|) in the ε-form | Section 6.3 | standard \|μ\|/(1+\|μ\|) implemented |
| D-5 | Lavasa et al. 2026: "0.94" MeV with 0.324 MV (Table 2); 2.5 MeV with 2.46 MV (Table 3) | Section 4.3 | flagged in CSV |
| D-6 | Wang et al. 2014 and AMPS Type V print R^{s−2} with R a Larmor radius (dimensionally inconsistent unless R = R_Lk_min) | Section 6.4 | Qin & Wang 2015 form implemented |
| D-7 | Borovikov et al. 2019 state k₀⁻¹ = L_max/2π and evaluate 81/(7π(2π)^{2/3}) "≈ 0.92"; the expression equals 1.0817 (λ_μμ analogue 0.5609), so the printed 0.9 and 0.5 are inconsistent with the stated convention | Sections 6.7, 8.4 | exact expressions used |
| D-8 | Le Roux et al. 1999 coefficient 2.433 (via Quenby & Webber) unreconciled with 3.1371 | Section 7.5 | not implemented |
| D-9 | Printed bracket exponents: Vos & Potgieter 2016 "(c−a)/c"; Aslam et al. 2021 "2∥ − 1" and "(c₂⊥ − c₃)/c₃" | Section 9.2 | interpretation must be configured explicitly |
| D-10 | Potgieter et al. 2014 Sec. 5 (arXiv v3): λ∥(1 GV) "increased by a factor of ~2.3, from ~0.13 AU in 2006, to ~0.3 AU in 2009" (v1–v2: "~23", "~30 AU") | Section 9.3 | Table 1 with the paper's Eq. 5 (Eq. 33 here) implies 0.146 → 0.230 AU (F-GCR-02) |
| D-11 | Aslam et al. 2021 (unit 6 × 10²⁰) and 2023 (unit 10²²) electron (K∥)₀ numerically similar | Section 9.3 | stored with units; not compared |
| D-12 | Tomassetti 2017: v1 prints "B₀ ≅ 3.4 nT AU²"; the final version (v5) prints "B₀ ≅ 3.4 nT sets the HMF at r₀ = 1 AU" | Section 9.5 | final version used |
| D-13 | Fiandrini et al. 2021 g(θ) printed with "θ_A + π/2" (NWU: θ_A − 90°) | Section 9.5 | as printed |
| D-14 | Chhiber 2022 turbulence-energy and δb/B₀ magnitudes mutually inconsistent under a direct reading | Section 11.3 | not used as normalization |
| D-15 | Breech et al. 2008: the only readable rendering shows σ_D = 1/3 with r_A = 1/2 (which requires −1/3), but it drops minus signs throughout; the printed sign is not established | Section 11.3 | r_A = 1/2 is the quantity used; the sign is to be read on the typeset page |
| D-16 | Oughton et al. 2011 boundary T "1.6 × 10⁶ K" (Fig. 1 caption, re-read in a copy of the published paper) versus Wiengarten et al. 2016 "1.6 × 10⁵ K" for the same case | Section 11.7 | both as printed |
| D-17 | Kozarev et al. 2013: λ₀ "at 1 AU and 1 GV" with pc/1 GeV in the equation | Section 5.2 | numerically identical for every singly charged species (\|Z\| = 1); differs for \|Z\| > 1 |
| D-18 | Potgieter et al. 2015 electron table header "λ∥ (AU) × 10⁻¹" ambiguous | Section 9.3 | values as printed |
| D-19 | Vos & Potgieter 2015 B₀, P₀ illegible ("=B 10 nT and =P 10 GV") | Section 9.3 | 2016 values (1 nT, 1 GV) noted; used only in the conditional fixture F-GCR-04 |
| D-20 | Aslam et al. 2023 electron P_k: 0.50–0.60 GV (Sec. 6.4) versus 0.4 GV (Sec. 7) | Section 9.3 | both recorded |
| D-21 | HelMod SSN validity: 2.2–266.9 (text) versus about 10–165 (Fig. 1 caption); K₀ unit carries GV⁻¹ while the bracket [P/(1 GV) + g_low] is dimensionless | Section 9.4 | both ranges recorded; the numerical-use convention (K₀ inserted as the printed number, result in AU² s⁻¹) is named in the configuration, the printed unit is kept as metadata, and K₀ is not rescaled as though the bracket carried a factor of GV |
| D-22 | Fiandrini et al. 2021 results differ between arXiv v1 and v2 | Section 9.5 | v2 quoted, v1 in data |
| D-23 | Tomassetti et al. 2023 (Eq. 3) print λ∥ with K₀ "in units of 10²³ cm² s⁻¹" and with the rigidity factor (R₀/R)^a, inverted relative to Fiandrini et al. 2021 and to their own statement that a is the slope below R_k | Section 9.5 | as printed; not implemented |
| D-24 | Strauss et al. 2011 print "κ∥ = v/3λ∥" | Section 9.5 | as printed; κ = vλ/3 is the conventional reading |
| D-25 | Chen et al. 2024 abstract "increases exponentially" for a power-law fit | Section 4.5 | power law implemented |
| D-26 | Strauss et al. 2020: the k_d equation is numbered (9) in the arXiv PDF and (13) in the ar5iv rendering | Sections 4.4, 7.7 | both numbers given |
| D-27 | Engelbrecht & Burger 2013a k_m = 1/λ_sl versus Lang et al. 2024 Eq. 13 | Sections 7.1, 9.5 | both recorded as choices |
| D-28 | Repository records of Potgieter et al. 2014 (5 authors) and Aslam et al. 2021 (5 authors) are incomplete; the journal versions list 6 and 7 authors | references | journal author lists used |
| D-29 | Minoshima et al. 2026 Eq. (9) prints λ∥ = ξv/Ω_n; Ω_n is not defined (the text gives ν_p = Ω_p/2π) | Section 5.2 | gyrofrequency must be stated by the configuration |
| D-30 | Afanasiev et al. 2015 Eq. (2) (arXiv version) prints R² in the denominator; the text defines B as the mean field | Section 6.8 | B² used; journal version not checked |
| D-31 | Perri et al. 2020: quoted λ∥ end values (0.02 AU at 1 MV, 200 AU at 10⁵ MV) and slopes (P^{0.33}, P^{1.31}) are mutually inconsistent; 100 MeV values 0.05–0.2 AU (Sec. 3.1) versus 0.02–0.2 AU (Sec. 3.2) | Section 7.5 | not used as benchmarks |
| D-32 | Lang et al. 2024 Table 5 prints the B₀ unit as "m nT" | Section 7.4 | nT assumed by the unit context; flagged |
| D-33 | Duan et al. 2025 print k₂ = (R/R_k)^a[1 + (R/R_k)^{(b−a)/c}]^c and state that a applies below R_k and b above; the printed form has these limiting slopes only for b > a, and the reverse for b < a (Eq. 47); the fitted a, b are given only in a figure | Sections 9.5, 21.8 | printed form behind a variant name; both orderings tested (F-GCR-08, F-GCR-09); the equation is not rewritten to match the prose |
| D-34 | Raath et al. 2016 call B₀ "the magnitude of the HMF at r₀ = 1 AU", but in their field (Eq. 42) B₀ is the radial-component amplitude at r₀; the magnitude there is B₀/\|cos ψ(r₀)\| (1.4071B₀ in the ecliptic for 400 km/s); the tabulated values are labelled B_e | Section 11.8 | normalization configured explicitly (U-13); never inferred |

### 14.2 Decisions required from the user or host code

| ID | Decision |
|---|---|
| U-1 | Whether a run reports SEP-CHEN24 as the published κ∥ or as the derived λ∥ = 3κ∥/v (Eq. 1; the conversion is a definition, but it needs the species, U-2) |
| U-2 | Which species SEP-CHEN24 applies to (protons implied, not stated) |
| U-3 | Whether the code label "Chen2024AA" refers to Chen et al. 2024, ApJ 965, 61 (inference) |
| U-4 | What the code mode "Tenishev2005AIAA" implements; no formula found in any accessible source |
| U-5 | Whether to offer the printed (1−μ)² variant of Lang/Lavasa |
| U-6 | Supplying the PARASOL Λ(E) (their Eq. 18) and Δx(E) (their Eq. 29) expressions from the source |
| U-7 | Which DT electron printing (if any) to support |
| U-8 | Whether an EPREM-compatible mode reproduces λ∥ or 2λ∥ |
| U-9 | Default turbulence inputs: none are provided; a configuration must name its source |
| U-10 | Interpretation of the D-9 exponents |
| U-11 | Interpretation of λ_unspecified presets (SEP-PATH09, SEP-SOLPENCO05) |
| U-12 | Variance convention of inputs to QLT-ZANK98 |
| U-13 | Whether a published 1 AU field value normalizes the field magnitude or the radial amplitude B₀ of Eq. (42) (D-34) |

### 14.3 Information not found

The following were sought and not found in accessible sources; nothing in this document substitutes for them.

- **SEP codes:** the sign of the 3/2 exponent of the radial factor in Schwadron et al. (2010, Eq. 2; minus signs lost in both readable copies); Zhang, Qin & Rassoul 2009 parameter values; full texts of Dröge et al. 2010, 2014, 2016 (ApJ/JGR); Hu et al. 2017 shock-region κ; PATH foundations (Zank et al. 2000; Rice et al. 2003; Li et al. 2003, 2005); Ng, Reames & Tylka 2003; Sokolov et al. 2004; Luhmann et al. 2007; Young et al. 2021 λ₀; SaRoN model description; Qin et al. 2006 (source of h = 0.2); Ruffolo-group φ(μ).
- **Observations:** full texts of Palmer 1982 and Bieber et al. 1994 (including whether Bieber's values are radial or parallel); Dresing et al. 2023 transport section; event fits by Tan et al. 2013, Kartavykh, Laitinen, Wei, Strauss 2017, Rodríguez-García; Zhong et al. 2024 radial exponent; a machine-readable archive of event-fitted λ.
- **Theory:** Shalchi 2026 (ApJ 1004, 125); Shalchi & Klippenstein focusing λ∥; Engelbrecht 2015, 2019; Burger et al. 2008 original equation; Zank et al. 1998 JGR full text; Le Roux et al. 1999; Pei et al. 2010 λ∥ formula; Dröge 2003; electron approximations TS2002 Eq. 67, TS2003 Eqs. 52, 54 (illegible).
- **GCR:** Song et al. 2021; Wang et al. 2022 (PRD 106, 063006) content; Moloto & Engelbrecht papers; Kopp et al. 2012; Zhao/Zhu/Luo/Shen 2017–2021 papers; Jokipii & Kopriva 1979, Jokipii & Thomas 1981, Langner et al. 2003, Heber et al. 2006 forms; Corti et al. 2019 complete Table 3; Bobik et al. 2012 units; Aslam et al. 2021, 2023 c₃ and P_k values; Perugia K₀(t), a(t), b(t), Wang 2019 B_c(t) and Duan 2025 a, b, R_k, K₀ (figures only); HelMod-4 2023 K₀–SSN coefficients.
- **Turbulence:** the typeset sign of σ_D in Breech et al. 2008 (D-15); Bieber et al. 1996 full text; Matthaeus et al. 1990 numbers; Alonso Guzmán et al. 2025 values; Hamilton et al. 2008 and Dasso et al. 2005 fractions; dissipation-range index distributions of Leamon et al. 1998 and Smith et al. 2006; Zank et al. 1996 and 2017 equations; Shiota et al. 2017; primary Parker 1958 and Smith & Bieber 1991 forms.

The web-search budget of the review was exhausted before these could be pursued further.

## 15. Model selection

These are factual consequences of the sources, not recommendations of one model over another.

- **Reproducing a published SEP simulation:** use the paper's preset (Section 5.2) with its quantity tag, its D_μμ shape and H/ε (Section 6), and its perpendicular pairing (Section 10). The same headline value (for example 0.3 AU) means different things in SEP-MFLAMPA25 (λ∥ at 1 GeV, ∝ r), SEP-SOFIE24 (constant λ∥, called parallel once in the source), SEP-PARADISE19 (λ_r at 4 MeV) and SEP-LAITINEN16 (λ∥ at 10 MeV from turbulence).
- **Comparing with event fits:** compare λ_r with λ_r (Agueda, Pacheco) and λ∥ with λ∥ (Lang, Lavasa), converting only with a stated ψ.
- **Turbulence-driven λ∥:** Eq. (21) for ions and Eq. (24) are the two closed forms in use; both require the spectrum convention of Section 7.1. Electron λ∥ is sensitive to k_d and p, and the DT form has the unresolved printing D-1.
- **GCR modulation:** each published parameter set is valid only with its own formula variant, unit, B at Earth, and companion parameters (tilt, drift reduction, heliopause and termination-shock positions). Parameters from different groups, or from Perugia papers with different h, are not interchangeable.
- **Near-Sun use:** Chen et al. (2024) is the only fitted κ∥(r, E) for 0.1–0.8 au; the Palmer band is a 1 au statement and Subashchandar et al. (2025) argue against extrapolating it inward.

## 16. Publication methods text

The following paragraphs are templates; bracketed items are filled from the configuration. They cite this document's reference identifiers (Section 17).

**Prescribed SEP mean free path.** "The parallel mean free path was prescribed as λ∥ = λ₀(r/1 au)(pc/1 GeV)^{1/3} with λ₀ = [value] au, following Liu et al. (2025) [M32]; the corresponding parallel diffusion coefficient is κ∥ = vλ∥/3. Pitch-angle scattering used D_μμ = D₀(1−μ²)(|μ|^{q−1}+H) with q = 5/3 and H = [value], with D₀ fixed exactly by the λ∥ integral (Agueda & Vainio 2013 [M58])."

**Radial mean free path.** "We prescribe a spatially constant radial mean free path λ_r = λ∥cos²ψ = [value] au at [energy] (convention of Agueda et al. 2010 [M13]); λ∥ therefore varies along the field line as 1/cos²ψ."

**Turbulence-based λ∥.** "The proton parallel mean free path was computed from the slab-turbulence expression of Teufel & Schlickeiser (2003) [M64] in the continuous form used by Engelbrecht & Burger (2014) [M70] and Lang et al. (2024) [M12], with B₀ = [ ] nT, δB²_slab = [ ] nT², s = [ ], k_min = [ ] km⁻¹." For electrons add the dissipation-range term and name the printed variant used.

**Zank fit.** "λ∥ was evaluated with the slab-QLT fit of Zank et al. (1998) [M71] in the form of Chhiber et al. (2017) [M03], using the total slab variance and the slab correlation length λ_s = [ ] (coefficient 6.2742)."

**GCR modulation.** "The parallel diffusion coefficient followed Potgieter et al. (2014) [M85], K∥ = (K∥)₀β(B₀/B)(P/P₀)^a[((P/P₀)^c + (P_k/P₀)^c)/(1 + (P_k/P₀)^c)]^{(b−a)/c}, with (K∥)₀ = [ ] × 10²² cm² s⁻¹, P₀ = 1 GV, B₀ = 1 nT, a = [ ], b = [ ], c = [ ], P_k = [ ] GV."

**Observational comparison.** "Model λ∥ at 1 au is compared with the Palmer (1982) [M01] consensus range 0.08–0.3 AU for 0.5 MV–5 GV and with event-fitted values of [source] ([λ∥ or λ_r])."

## 17. References

References are numbered [Mxx] in order of first citation. Each entry gives the reading level for the content used here (Section 1.2) and the citation key: the existing key when the paper is in the supplied bibliography, otherwise the proposed key of the additions file (Section 18). Entries marked with a metadata note have bibliographic fields that were not all confirmed against the publisher record; they must be checked before publication citation (Section 20 records what was verified).

**[M01]** Palmer, I. D. (1982). **Transport coefficients of low-energy cosmic rays in interplanetary space.** *Reviews of Geophysics*, **20**(2), 335–351. [DOI](https://doi.org/10.1029/RG020i002p00335). Role: Palmer consensus for λ∥ at 1 AU. Reading: abstract only. Pre-2006 foundation. Proposed citation key (additions file): **Palmer-1982-RG**.

**[M02]** Tautz, R. C.; Shalchi, A. (2013). **Simulated energetic particle transport in the interplanetary space: The Palmer consensus revisited.** *Journal of Geophysical Research: Space Physics*, **118**(2), 642–647. [DOI](https://doi.org/10.1002/jgra.50155); [arXiv:1301.7162](https://arxiv.org/abs/1301.7162). Role: Test-particle λ versus Palmer band. Reading: full text. Citation key in the supplied bibliography: **Tautz-2013-JGRSP**.

**[M03]** Chhiber, R.; Subedi, P.; Usmanov, A. V.; Matthaeus, W. H.; Ruffolo, D.; Goldstein, M. L.; Parashar, T. N. (2017). **Cosmic-Ray Diffusion Coefficients throughout the Inner Heliosphere from a Global Solar Wind Simulation.** *The Astrophysical Journal Supplement Series*, **230**(2), 21. [DOI](https://doi.org/10.3847/1538-4365/aa74d2); [arXiv:1703.10322](https://arxiv.org/abs/1703.10322). Role: Turbulence-transport-derived λ∥, λ⊥ in the inner heliosphere. Reading: full text. Citation key in the supplied bibliography: **Chhiber-2017-AJSS**.

**[M04]** Reames, D. V. (2013). **The Two Sources of Solar Energetic Particles.** *Space Science Reviews*, **175**(1-4), 53–92. [DOI](https://doi.org/10.1007/s11214-013-9958-9); [arXiv:1306.3608](https://arxiv.org/abs/1306.3608). Role: Restatement of the Palmer consensus. Reading: full text. Citation key in the supplied bibliography: **Reames-2013-SSR**.

**[M05]** Lavasa, E.; Lang, J. T.; Papaioannou, A.; Strauss, R. D.; Mallios, S. A.; Hillaris, A.; Kouloumvakos, A.; Anastasiadis, A.; Daglis, I. A. (2026). **Multi-spacecraft constraints on relativistic solar energetic particle transport in the widespread 28 October 2021 event.** *Astronomy & Astrophysics*, **707**, A12. [DOI](https://doi.org/10.1051/0004-6361/202558094); [arXiv:2603.09839](https://arxiv.org/abs/2603.09839). Role: GLE 73 λ∥ and λ⊥ versus rigidity. Reading: full text. Citation key in the supplied bibliography: **Lavasa-2026-AA**.

**[M06]** Minoshima, T.; Miyoshi, Y.; Murakami, G.; Pinto, M.; Schmid, D.; Matsuoka, A.; Baumjohann, W.; Fischer, D.; Iwai, K.; Imada, S. (2026). **A study of solar energetic particle transport on 30 March 2022 using multi-spacecraft data assimilation.** *Earth, Planets and Space*, **78**(1), 64. [DOI](https://doi.org/10.1186/s40623-026-02389-9); [arXiv:2602.00765](https://arxiv.org/abs/2602.00765). Role: Data-assimilated λ∥ = ξ v/Ω_p. Reading: full text. Citation key in the supplied bibliography: **Minoshima-2026-EPS**.

**[M07]** Subashchandar, N. S. M.; Zhao, L.; Shalchi, A.; Zank, G. P.; le Roux, J. A.; Li, H.; Zhu, X.; Silwal, A.; Alonso Guzman, J. (2025). **Parallel and Perpendicular Diffusion of Energetic Particles in the Near-Sun Solar Wind Observed by Parker Solar Probe.** *The Astrophysical Journal Letters*, **991**(2), L30. [DOI](https://doi.org/10.3847/2041-8213/ae063f); [arXiv:2509.10648](https://arxiv.org/abs/2509.10648). Role: PSP SOQLT/UNLT radial scalings. Reading: full text. Citation key in the supplied bibliography: **Subashchandar-2025-AJL**.

**[M08]** Shalchi, A.; Bieber, J. W.; Matthaeus, W. H.; Schlickeiser, R. (2006). **Parallel and Perpendicular Transport of Heliospheric Cosmic Rays in an Improved Dynamical Turbulence Model.** *The Astrophysical Journal*, **642**(1), 230–243. [DOI](https://doi.org/10.1086/500728). Role: Dynamical turbulence λ∥, λ⊥ with heliospheric parameters. Reading: full text. Proposed citation key (additions file): **Shalchi-2006-AJ**.

**[M09]** Bieber, J. W.; Matthaeus, W. H.; Smith, C. W.; Wanner, W.; Kallenrode, M.-B.; Wibberenz, G. (1994). **Proton and Electron Mean Free Paths: The Palmer Consensus Revisited.** *The Astrophysical Journal*, **420**, 294–306. [DOI](https://doi.org/10.1086/173559). Role: Electron/proton mean-free-path discrepancy; slab/2D composite. Reading: abstract only. Pre-2006 foundation. Proposed citation key (additions file): **Bieber-1994-AJ-420**.

**[M10]** Engelbrecht, N. E.; Effenberger, F.; Florinski, V.; Potgieter, M. S.; Ruffolo, D.; Chhiber, R.; Usmanov, A. V.; Rankin, J. S.; Els, P. L. (2022). **Theory of Cosmic Ray Transport in the Heliosphere.** *Space Science Reviews*, **218**(4), 33. [DOI](https://doi.org/10.1007/s11214-022-00896-1). Role: Review of observed and theoretical mean free paths. Reading: full text. Citation key in the supplied bibliography: **Engelbrecht-2022-SSRv**.

**[M11]** Chen, X.; Giacalone, J.; Guo, F.; Klein, K. G. (2024). **Parallel Diffusion Coefficient of Energetic Charged Particles in the Inner Heliosphere from the Turbulent Magnetic Fields Measured by Parker Solar Probe.** *The Astrophysical Journal*, **965**(1), 61. [DOI](https://doi.org/10.3847/1538-4357/ad33c3); [arXiv:2403.08141](https://arxiv.org/abs/2403.08141). Role: PSP κ∥(r,E) power-law fit. Reading: full text. Citation key in the supplied bibliography: **Chen-2024-AJ-965**.

**[M12]** Lang, J. T.; Strauss, R. D.; Engelbrecht, N. E.; van den Berg, J. P.; Dresing, N.; Ruffolo, D.; Bandyopadhyay, R. (2024). **A Detailed Survey of the Parallel Mean Free Path of Solar Energetic Particle Protons and Electrons.** *The Astrophysical Journal*, **971**(1), 105. [DOI](https://doi.org/10.3847/1538-4357/ad55c3); [arXiv:2406.05765](https://arxiv.org/abs/2406.05765). Role: λ∥ survey of 15 events; TS2003 closed form with fitted turbulence parameters. Reading: full text. Citation key in the supplied bibliography: **Lang-2024-arXiv**.

**[M13]** Agueda, N.; Vainio, R.; Lario, D.; Sanahuja, B. (2010). **Solar near-relativistic electron observations as a proof of a back-scatter region beyond 1 AU during the 2000 February 18 event.** *Astronomy & Astrophysics*, **519**, A36. [DOI](https://doi.org/10.1051/0004-6361/200913963). Role: λ_r definition, ε-form D_μμ, back-scatter region. Reading: full text. Proposed citation key (additions file): **Agueda-2010-AA**.

**[M14]** Agueda, N.; Klein, K.-L.; Vilmer, N.; Rodríguez-Gasén, R.; Malandraki, O. E.; Papaioannou, A.; Subirà, M.; Sanahuja, B.; Valtonen, E.; Dröge, W.; Nindos, A.; Heber, B.; Braune, S.; Usoskin, I. G.; Heynderickx, D.; Talew, E.; Vainio, R. (2014). **Release timescales of solar energetic particles in the low corona.** *Astronomy & Astrophysics*, **570**, A5. [DOI](https://doi.org/10.1051/0004-6361/201423549). Role: Event-fitted electron λ_r (Table 3). Reading: full text. Proposed citation key (additions file): **Agueda-2014-AA**.

**[M15]** Agueda, N.; Lario, D. (2016). **Release History and Transport Parameters of Relativistic Solar Electrons Inferred from Near-the-Sun In Situ Observations.** *The Astrophysical Journal*, **829**(2), 131. [DOI](https://doi.org/10.3847/0004-637X/829/2/131). Role: Helios event-fitted λ_r. Reading: abstract only. Proposed citation key (additions file): **Agueda-2016-AJ**.

**[M16]** Pacheco, D.; Agueda, N.; Aran, A.; Heber, B.; Lario, D. (2019). **Full inversion of solar relativistic electron events measured by the Helios spacecraft.** *Astronomy & Astrophysics*, **624**, A3. [DOI](https://doi.org/10.1051/0004-6361/201834520); [arXiv:1902.06602](https://arxiv.org/abs/1902.06602). Role: Helios λ_r versus r (Tables 2-4). Reading: full text. Proposed citation key (additions file): **Pacheco-2019-AA**.

**[M17]** Dröge, W.; Kartavykh, Y. Y.; Dresing, N.; Klassen, A. (2016). **Multi-Spacecraft Observations and Transport Modeling of Energetic Electrons for a Series of Solar Particle Events in August 2010.** *The Astrophysical Journal*, **826**(2), 134. [DOI](https://doi.org/10.3847/0004-637X/826/2/134). Role: λ∥ and λ⊥ normalized to 1 au. Reading: abstract only. Proposed citation key (additions file): **Droge-2016-AJ**.

**[M18]** Agueda, N.; Vainio, R.; Lario, D.; Sanahuja, B. (2008). **Injection and Interplanetary Transport of Near-Relativistic Electrons: Modeling the Impulsive Event on 2000 May 1.** *The Astrophysical Journal*, **675**(2), 1601–1613. [DOI](https://doi.org/10.1086/527527). Role: Event-fitted λ_r. Reading: abstract only. Proposed citation key (additions file): **Agueda-2008-AJ**.

**[M19]** Dröge, W. (2000). **The Rigidity Dependence of Solar Particle Scattering Mean Free Paths.** *The Astrophysical Journal*, **537**(2), 1073–1079. [DOI](https://doi.org/10.1086/309080). Role: Rigidity dependence of event-fitted λ. Reading: abstract only. Pre-2006 foundation. Proposed citation key (additions file): **Droge-2000-AJ**.

**[M20]** Battarbee, M.; Guo, J.; Dalla, S.; Wimmer-Schweingruber, R.; Swalwell, B.; Lawrence, D. J. (2018). **Multi-spacecraft observations and transport simulations of solar energetic particles for the May 17th 2012 event.** *Astronomy & Astrophysics*, **612**, A116. [DOI](https://doi.org/10.1051/0004-6361/201731451); [arXiv:1706.08458](https://arxiv.org/abs/1706.08458). Role: Constant λ = 0.3 au full-orbit model input. Reading: full text. Proposed citation key (additions file): **Battarbee-2018-AA**.

**[M21]** Dalla, S.; de Nolfo, G. A.; Bruno, A.; Giacalone, J.; Laitinen, T.; Thomas, S.; Battarbee, M.; Marsh, M. S. (2020). **3D propagation of relativistic solar protons through interplanetary space.** *Astronomy & Astrophysics*, **639**, A105. [DOI](https://doi.org/10.1051/0004-6361/201937338); [arXiv:2002.00929](https://arxiv.org/abs/2002.00929). Role: Constant λ for GLE protons. Reading: full text. Citation key in the supplied bibliography: **Dalla-2020-AA**.

**[M22]** Houeibib, A.; Pantellini, F.; Griton, L. (2025). **Dynamics of energetic electrons scattered in the solar wind: Magnetohydrodynamics and test-particle simulations.** *Astronomy & Astrophysics*, **694**, A211. [DOI](https://doi.org/10.1051/0004-6361/202451436); [arXiv:2403.06706](https://arxiv.org/abs/2403.06706). Role: Constant λ∥ 0.1-1 AU in electron PAD modelling. Reading: full text. Citation key in the supplied bibliography: **Houeibib-2025-AA**.

**[M23]** Dröge, W.; Kartavykh, Y. Y. (2009). **Testing Transport Theories with Solar Energetic Particles.** *The Astrophysical Journal*, **693**(1), 69–74. [DOI](https://doi.org/10.1088/0004-637X/693/1/69). Role: Dynamical-QLT pitch-angle distributions vs observations. Reading: secondary. Citation key in the supplied bibliography: **Droge-2009-AJ**.

**[M24]** Tan, L. C.; Reames, D. V.; Ng, C. K.; Shao, X.; Wang, L. (2011). **What Causes Scatter-Free Transport of Non-Relativistic Solar Electrons?.** *The Astrophysical Journal*, **728**(2), 133. [DOI](https://doi.org/10.1088/0004-637X/728/2/133). Role: Scatter-free to diffusive electron transition. Reading: abstract only. Proposed citation key (additions file): **Tan-2011-AJ**.

**[M25]** Strauss, R. D.; Dresing, N.; Kollhoff, A.; Brüdern, M. (2020). **On the shape of SEP electron spectra: The role of interplanetary transport.** *The Astrophysical Journal*, **897**(1), 24. [DOI](https://doi.org/10.3847/1538-4357/ab91b0); [arXiv:2005.03486](https://arxiv.org/abs/2005.03486). Role: Dissipation-range onset for electron λ∥. Reading: full text. Proposed citation key (additions file): **Strauss-2020-AJ**.

**[M26]** Zhong, Y.; Wang, Y.; Qin, G. (2024). **The Mean Free Path of 13-64 MeV Protons Derived from Statistical Results of Solar Energetic Particle Events.** *The Astrophysical Journal*, **974**(2), 228. [DOI](https://doi.org/10.3847/1538-4357/ad70aa). Role: Radial power-law λ_r (exponent not obtained). Reading: abstract only. Citation key in the supplied bibliography: **Zhong-2024-AJ-974**.

**[M27]** Cao, Y.; Wang, Y.; Guo, J. (2025). **Radial dependence of solar energetic particle peak fluxes and fluences.** *Astronomy & Astrophysics*, **695**, A25. [DOI](https://doi.org/10.1051/0004-6361/202452591). Role: κ = κ0(E/E0)^β in a peak-flux model. Reading: full text. Citation key in the supplied bibliography: **Cao-2025-AA**.

**[M28]** Verkhoglyadova, O. P.; Li, G.; Zank, G. P.; Hu, Q.; Mewaldt, R. A. (2009). **Using the Path Code for Modeling Gradual SEP Events in the Inner Heliosphere.** *The Astrophysical Journal*, **693**(1), 894–900. [DOI](https://doi.org/10.1088/0004-637X/693/1/894). Role: PATH ambient λ power law. Reading: full text. Citation key in the supplied bibliography: **Verkhoglyadova-2009-AJ**.

**[M29]** Schwadron, N. A.; Townsend, L.; Kozarev, K.; Dayeh, M. A.; Cucinotta, F.; Desai, M.; Golightly, M.; Hassler, D.; Hatcher, R.; Kim, M.-Y.; Posner, A.; PourArsalan, M.; Spence, H. E.; Squier, R. K. (2010). **Earth-Moon-Mars Radiation Environment Module framework.** *Space Weather*, **8**(1), S00E02. [DOI](https://doi.org/10.1029/2009SW000523). Role: EPREM λ∥ normalization. Reading: full text. Citation key in the supplied bibliography: **Schwadron-2010-SW**. Metadata note: unconfirmed field(s): pages.

**[M30]** Kozarev, K. A.; Evans, R. M.; Schwadron, N. A.; Dayeh, M. A.; Opher, M.; Korreck, K. E.; van der Holst, B. (2013). **Global Numerical Modeling of Energetic Proton Acceleration in a Coronal Mass Ejection Traveling through the Solar Corona.** *The Astrophysical Journal*, **778**(1), 43. [DOI](https://doi.org/10.1088/0004-637X/778/1/43); [arXiv:1406.2377](https://arxiv.org/abs/1406.2377). Role: EPREM D_μμ and λ∥ power law. Reading: full text. Citation key in the supplied bibliography: **Kozarev2013**.

**[M31]** Borovikov, D.; Sokolov, I. V.; Huang, Z.; Roussev, I. I.; Gombosi, T. I. (2019). **Toward quantitative model for simulation and forecast of solar energetic particle production during gradual events - II: kinetic description of SEP.** *arXiv preprint*. [arXiv:1911.10165](https://arxiv.org/abs/1911.10165). Role: M-FLAMPA D_μμ shape and λ_xx forms. Reading: full text. Citation key in the supplied bibliography: **Borovikov-2019-arXiv**.

**[M32]** Liu, W.; Sokolov, I. V.; Zhao, L.; Gombosi, T. I.; Sachdeva, N.; Chen, X.; Tóth, G.; Lario, D.; Manchester IV, W. B.; Whitman, K.; Cohen, C. M. S.; Bruno, A.; Mays, M. L.; Bain, H. M. (2025). **Physics-based Simulation of the 2013 April 11 Solar Energetic Particle Event.** *The Astrophysical Journal*, **985**(1), 82. [DOI](https://doi.org/10.3847/1538-4357/adc4e3); [arXiv:2412.07581](https://arxiv.org/abs/2412.07581). Role: M-FLAMPA upstream and downstream λ∥. Reading: full text. Citation key in the supplied bibliography: **Liu-2025-AJ**.

**[M33]** Zhao, L.; Sokolov, I.; Gombosi, T.; Lario, D.; Whitman, K.; Huang, Z.; Toth, G.; Manchester, W.; van der Holst, B.; Sachdeva, N.; Liu, W. (2024). **Solar Wind with Field Lines and Energetic Particles (SOFIE) Model: Application to Historical Solar Energetic Particle Events.** *Space Weather*, **22**(9), e2023SW003729. [DOI](https://doi.org/10.1029/2023SW003729); [arXiv:2309.16903](https://arxiv.org/abs/2309.16903). Role: SOFIE constant upstream λ. Reading: full text. Citation key in the supplied bibliography: **Zhao-2024-SW**.

**[M34]** Zhang, M.; Cheng, L.; Zhang, J.; Riley, P.; Kwon, R.-Y.; Lario, D.; Balmaceda, L.; Pogorelov, N. V. (2023). **A Data-driven, Physics-based Transport Model of Solar Energetic Particles Accelerated by Coronal Mass Ejection Shocks Propagating through the Solar Coronal and Heliospheric Magnetic Fields.** *The Astrophysical Journal Supplement Series*, **266**(2), 35. [DOI](https://doi.org/10.3847/1538-4365/accb8e); [arXiv:2212.07259](https://arxiv.org/abs/2212.07259). Role: Constant λ_r, h0 = 0.2, Bohm shock limit. Reading: full text. Citation key in the supplied bibliography: **Zhang-2023-AJSS**.

**[M35]** Afanasiev, A.; Wijsen, N.; Vainio, R. (2025). **Towards advanced forecasting of solar energetic particle events with the PARASOL model.** *Journal of Space Weather and Space Climate*, **15**, 3. [DOI](https://doi.org/10.1051/swsc/2024039); [arXiv:2412.11852](https://arxiv.org/abs/2412.11852). Role: PARASOL background and shock-region λ∥. Reading: full text. Citation key in the supplied bibliography: **Afanasiev-2025-JSWSC, Afanasiev-2024-JSWSC**.

**[M36]** Aran, A.; Sanahuja, B.; Lario, D. (2005). **Fluxes and fluences of SEP events derived from SOLPENCO.** *Annales Geophysicae*, **23**(9), 3047–3053. [DOI](https://doi.org/10.5194/angeo-23-3047-2005). Role: SOLPENCO λ ∝ P^1/2. Reading: full text. Pre-2006 foundation. Citation key in the supplied bibliography: **Aran-2005-AG**.

**[M37]** Marsh, M. S.; Dalla, S.; Dierckxsens, M.; Laitinen, T.; Crosby, N. B. (2015). **SPARX: A modeling system for Solar Energetic Particle Radiation Space Weather forecasting.** *Space Weather*, **13**(6), 386–394. [DOI](https://doi.org/10.1002/2014SW001120). Role: Isotropic constant λ default. Reading: full text. Citation key in the supplied bibliography: **Marsh-2015-SW**.

**[M38]** Marsh, M. S.; Dalla, S.; Kelly, J.; Laitinen, T. (2013). **Drift-induced Perpendicular Transport of Solar Energetic Particles.** *The Astrophysical Journal*, **774**(1), 4. [DOI](https://doi.org/10.1088/0004-637X/774/1/4); [arXiv:1307.1585](https://arxiv.org/abs/1307.1585). Role: Isotropic scattering full-orbit model. Reading: full text. Proposed citation key (additions file): **Marsh-2013-AJ**.

**[M39]** He, H.-Q.; Qin, G.; Zhang, M. (2011). **Propagation of Solar Energetic Particles in Three-Dimensional Interplanetary Magnetic Fields: In View of Characteristics of Sources.** *The Astrophysical Journal*, **734**(2), 74. [DOI](https://doi.org/10.1088/0004-637X/734/2/74). Role: Constant λ_r, h = 0.2. Reading: full text. Citation key in the supplied bibliography: **He-2011-AJ**.

**[M40]** Wang, Y.; Qin, G. (2015). **Estimation of the Release Time of Solar Energetic Particles Near the Sun.** *The Astrophysical Journal*, **799**(1), 111. [DOI](https://doi.org/10.1088/0004-637X/799/1/111); [arXiv:1311.7469](https://arxiv.org/abs/1311.7469). Role: h = 0.01 D_μμ, λ∥ at 10 MeV. Reading: full text. Citation key in the supplied bibliography: **Wang-2004-arXiv**.

**[M41]** Kubo, Y.; Kataoka, R.; Sato, T. (2015). **Interplanetary particle transport simulation for warning system for aviation exposure to solar energetic particles.** *Earth, Planets and Space*, **67**, 117. [DOI](https://doi.org/10.1186/s40623-015-0260-9); [arXiv:1506.00825](https://arxiv.org/abs/1506.00825). Role: Fitted constant λ_rr at 142 MeV. Reading: full text. Citation key in the supplied bibliography: **Kubo-2015-EPS, Kubo-2015-PS**.

**[M42]** Wijsen, N.; Aran, A.; Pomoell, J.; Poedts, S. (2019). **Modelling three-dimensional transport of solar energetic protons in a corotating interaction region generated with EUHFORIA.** *Astronomy & Astrophysics*, **622**, A28. [DOI](https://doi.org/10.1051/0004-6361/201833958); [arXiv:1901.09596](https://arxiv.org/abs/1901.09596). Role: PARADISE ε-form D_μμ and κ⊥. Reading: full text. Citation key in the supplied bibliography: **Wijsen-2019-AA-A28**.

**[M43]** Laitinen, T.; Kopp, A.; Effenberger, F.; Dalla, S.; Marsh, M. S. (2016). **Solar energetic particle access to distant longitudes through turbulent field-line meandering.** *Astronomy & Astrophysics*, **591**, A18. [DOI](https://doi.org/10.1051/0004-6361/201527801). Role: Turbulence-normalized λ∥, H = 0.1. Reading: full text. Citation key in the supplied bibliography: **Laitinen-2016-AA**.

**[M44]** Laitinen, T.; Effenberger, F.; Kopp, A.; Dalla, S. (2018). **The effect of turbulence strength on meandering field lines and Solar Energetic Particle event extents.** *Journal of Space Weather and Space Climate*, **8**, A13. [DOI](https://doi.org/10.1051/swsc/2018001). Role: λ∥ normalization variants. Reading: full text. Citation key in the supplied bibliography: **Laitinen-2018-JSWSC, Laitinen-2018-arXiv**.

**[M45]** Strauss, R. D.; Ogunjobi, O.; Moraal, H.; McCracken, K. G.; Caballero-Lopez, R. A. (2017). **On the Pulse Shape of Ground-Level Enhancements.** *Solar Physics*, **292**(4), 51. [DOI](https://doi.org/10.1007/s11207-017-1086-3); [arXiv:1703.05906](https://arxiv.org/abs/1703.05906). Role: Radial power-law λ_rr. Reading: full text. Proposed citation key (additions file): **Strauss-2017-SP**.

**[M46]** Strauss, R. D.; Fichtner, H. (2015). **On Aspects Pertaining to the Perpendicular Diffusion of Solar Energetic Particles.** *The Astrophysical Journal*, **801**(1), 29. [DOI](https://doi.org/10.1088/0004-637X/801/1/29); [arXiv:1804.03689](https://arxiv.org/abs/1804.03689). Role: H = 0.05 q-form; constant λ∥. Reading: full text. Citation key in the supplied bibliography: **Strauss-2015-AJ**.

**[M47]** Dröge, W.; Kartavykh, Y.; Klassen, A.; Dresing, N.; Lario, D. (2016). **Multi-spacecraft observations and transport modeling of energetic electrons for a series of solar particle events in August 2010.** *Proceedings of the 34th International Cosmic Ray Conference (ICRC2015), PoS(ICRC2015)*, **236**, 206. [DOI](https://doi.org/10.22323/1.236.0206). Role: Dröge SDE model, Λ⊥ scaling, fitted λ_r. Reading: full text. Proposed citation key (additions file): **Droge-2016-ICRC**.

**[M48]** Young, M. A.; Schwadron, N. A.; Gorby, M.; Linker, J.; Caplan, R. M.; Downs, C.; Török, T.; Riley, P.; Lionello, R.; Titov, V.; Mewaldt, R. A.; Cohen, C. M. S. (2021). **Energetic Proton Propagation and Acceleration Simulated for the Bastille Day Event of 2000 July 14.** *The Astrophysical Journal*, **909**(2), 160. [DOI](https://doi.org/10.3847/1538-4357/abdf5f); [arXiv:2012.09078](https://arxiv.org/abs/2012.09078). Role: EPREM κ∥ rigidity form. Reading: full text (partial). Citation key in the supplied bibliography: **Young-2021-AJ**.

**[M49]** He, H.-Q.; Wan, W. (2015). **Numerical Study of the Longitudinally Asymmetric Distribution of Solar Energetic Particles in the Heliosphere.** *The Astrophysical Journal Supplement Series*, **218**(2), 17. [DOI](https://doi.org/10.1088/0067-0049/218/2/17); [arXiv:1502.02683](https://arxiv.org/abs/1502.02683). Role: Constant λ_r, h = 0.2. Reading: full text. Citation key in the supplied bibliography: **He-2015-AJ-218**.

**[M50]** He, H.-Q. (2015). **Perpendicular Diffusion in the Transport of Solar Energetic Particles from Unconnected Sources: The Counter-streaming Particle Beams Revisited.** *The Astrophysical Journal*, **814**(2), 157. [DOI](https://doi.org/10.1088/0004-637X/814/2/157); [arXiv:1512.00027](https://arxiv.org/abs/1512.00027). Role: Constant λ_r for electrons. Reading: full text. Citation key in the supplied bibliography: **He-2015-AJ-814**.

**[M51]** He, H.-Q.; Wan, W. (2019). **Propagation of Solar Energetic Particles in the Outer Heliosphere: Interplay between Scattering and Adiabatic Focusing.** *The Astrophysical Journal Letters*, **885**(2), L28. [DOI](https://doi.org/10.3847/2041-8213/ab50bd); [arXiv:1911.02787](https://arxiv.org/abs/1911.02787). Role: λ_r/λ∥ pairs. Reading: full text. Citation key in the supplied bibliography: **He-2019-AJL, He-2019-arXiv**.

**[M52]** Tenishev, V. M.; Combi, M. R. (2005). **Monte-Carlo Model for Dust/Gas Interaction in Rarefied Flows.** *38th AIAA Thermophysics Conference, AIAA-2005-4832*. [DOI](https://doi.org/10.2514/6.2005-4832); [link](https://hdl.handle.net/2027.42/76515). Role: Origin of the code label Tenishev2005AIAA (no λ formula found). Reading: metadata only. Pre-2006 foundation. Citation key in the supplied bibliography: **Tenishev-2005-AIAA-4832**.

**[M53]** Tenishev, V.; Zhao, L.; Sokolov, I. (2022). **Application of the Monte Carlo Method in Modeling Transport and Acceleration of Solar Energetic Particles.** *arXiv preprint*. [arXiv:2209.09346](https://arxiv.org/abs/2209.09346). Role: AMPS D_μμ types I-V. Reading: full text. Citation key in the supplied bibliography: **Tenishev-2022-arXiv**.

**[M54]** Shalchi, A.; Yan, H.; Lazarian, A. (2005). **Spurious contribution to cosmic ray scattering calculations.** *Monthly Notices of the Royal Astronomical Society*, **356**(3), 1064–1070. [DOI](https://doi.org/10.1111/j.1365-2966.2004.08531.x); [arXiv:astro-ph/0411074](https://arxiv.org/abs/astro-ph/0411074). Role: Magnetostatic slab D_μμ and λ∥ asymptote. Reading: full text. Pre-2006 foundation. Proposed citation key (additions file): **Shalchi-2005-MNRAS**.

**[M55]** Wang, Y.; Qin, G.; Zhang, M.; Dalla, S. (2014). **A Numerical Simulation of Solar Energetic Particle Dropouts during Impulsive Events.** *The Astrophysical Journal*, **789**(2), 157. [DOI](https://doi.org/10.1088/0004-637X/789/2/157); [arXiv:1403.6544](https://arxiv.org/abs/1403.6544). Role: Slab-QLT-amplitude D_μμ with h = 0.01. Reading: full text. Citation key in the supplied bibliography: **Wang-2014-AJ**.

**[M56]** Qin, G.; Wang, Y. (2015). **Simulations of a Gradual Solar Energetic Particle Event Observed by Helios 1, Helios 2, and IMP 8.** *The Astrophysical Journal*, **809**(2), 177. [DOI](https://doi.org/10.1088/0004-637X/809/2/177); [arXiv:1505.02974](https://arxiv.org/abs/1505.02974). Role: (R_L k_min)^{s-2} D_μμ and NLGC κ⊥. Reading: full text. Citation key in the supplied bibliography: **Qin-2015-AJ**.

**[M57]** Dröge, W.; Kartavykh, Y. Y.; Klecker, B.; Kovaltsov, G. A. (2010). **Anisotropic Three-Dimensional Focused Transport of Solar Energetic Particles in the Inner Heliosphere.** *The Astrophysical Journal*, **709**(2), 912–919. [DOI](https://doi.org/10.1088/0004-637X/709/2/912). Role: 3D focused transport with λ⊥ scaling. Reading: abstract only. Citation key in the supplied bibliography: **Droge-2010-AJ**.

**[M58]** Agueda, N.; Vainio, R. (2013). **On the parametrization of the energetic-particle pitch-angle diffusion coefficient.** *Journal of Space Weather and Space Climate*, **3**, A10. [DOI](https://doi.org/10.1051/swsc/2013034). Role: q-form and ε-form D_μμ normalizations. Reading: full text. Proposed citation key (additions file): **Agueda-2013-JSWSC**.

**[M59]** Kelly, J.; Dalla, S.; Laitinen, T. (2012). **Cross-Field Transport of Solar Energetic Particles in a Large-Scale Fluctuating Magnetic Field.** *The Astrophysical Journal*, **750**(1), 47. [DOI](https://doi.org/10.1088/0004-637X/750/1/47); [arXiv:1202.6010](https://arxiv.org/abs/1202.6010). Role: Constant-λ full-orbit scattering. Reading: full text. Citation key in the supplied bibliography: **Kelly-2012-arXiv**.

**[M60]** Ding, Z.; Wijsen, N.; Li, G.; Poedts, S. (2022). **Modeling the 2020 November 29 solar energetic particle event using EUHFORIA and iPATH models.** *Astronomy & Astrophysics*, **668**, A71. [DOI](https://doi.org/10.1051/0004-6361/202244732); [arXiv:2210.16967](https://arxiv.org/abs/2210.16967). Role: iPATH QLT/NLGC coefficients and scalings. Reading: full text. Citation key in the supplied bibliography: **Ding-2022-AA, Ding-2022-arXiv**.

**[M61]** Afanasiev, A.; Battarbee, M.; Vainio, R. (2015). **Self-consistent Monte Carlo simulations of proton acceleration in coronal shocks: Effect of anisotropic pitch-angle scattering of particles.** *Astronomy & Astrophysics*, **584**, A81. [DOI](https://doi.org/10.1051/0004-6361/201526750); [arXiv:1603.08857](https://arxiv.org/abs/1603.08857). Role: Self-generated-wave D_μμ and steady-state λ(x,v). Reading: full text. Citation key in the supplied bibliography: **Afanasiev-2015-AA, Afanasiev-2016-arXiv, Afanasiev-misc**.

**[M62]** Strauss, R. D.; Dresing, N.; Engelbrecht, N. E. (2017). **Perpendicular Diffusion of Solar Energetic Particles: Model Results and Implications for Electrons.** *The Astrophysical Journal*, **837**(1), 43. [DOI](https://doi.org/10.3847/1538-4357/aa5df5); [arXiv:1804.03693](https://arxiv.org/abs/1804.03693). Role: Damped QLT D_μμ and FLRW for electrons. Reading: full text. Citation key in the supplied bibliography: **Strauss-2017-AJ-837**.

**[M63]** He, H.-Q.; Wan, W. (2013). **The dependence of the parallel and perpendicular mean free paths on the rigidity of the solar energetic particles: theoretical model versus observations.** *Astronomy & Astrophysics*, **557**, A57. [DOI](https://doi.org/10.1051/0004-6361/201321025). Role: Focusing-corrected λ∥. Reading: full text. Proposed citation key (additions file): **He-2013-AA**.

**[M64]** Teufel, A.; Schlickeiser, R. (2003). **Analytic calculation of the parallel mean free path of heliospheric cosmic rays. II. Dynamical magnetic slab turbulence and random sweeping slab turbulence with finite wave power at small wavenumbers.** *Astronomy & Astrophysics*, **397**(1), 15–25. [DOI](https://doi.org/10.1051/0004-6361:20021471). Role: Slab QLT λ∥ asymptotics used by the NWU constructions. Reading: full text. Pre-2006 foundation. Proposed citation key (additions file): **Teufel-2003-AA**.

**[M65]** Engelbrecht, N. E.; Burger, R. A. (2013). **An Ab Initio Model for the Modulation of Galactic Cosmic-Ray Electrons.** *The Astrophysical Journal*, **779**(2), 158. [DOI](https://doi.org/10.1088/0004-637X/779/2/158). Role: Electron RS/DT λ∥ constructions. Reading: full text. Citation key in the supplied bibliography: **Engelbrecht-2013-AJ-779**.

**[M66]** Teufel, A.; Schlickeiser, R. (2002). **Analytic calculation of the parallel mean free path of heliospheric cosmic rays. I. Dynamical magnetic slab turbulence and random sweeping slab turbulence.** *Astronomy & Astrophysics*, **393**(2), 703–715. [DOI](https://doi.org/10.1051/0004-6361:20021046). Role: Slab QLT λ∥ asymptotics (g = 0 below k_min). Reading: full text. Pre-2006 foundation. Proposed citation key (additions file): **Teufel-2002-AA**.

**[M67]** Shalchi, A.; Lazarian, A.; Schlickeiser, R. (2008). **Non-linear damping of slab modes and cosmic ray transport.** *Monthly Notices of the Royal Astronomical Society*, **383**(2), 803–808. [DOI](https://doi.org/10.1111/j.1365-2966.2007.12590.x); [arXiv:0710.5418](https://arxiv.org/abs/0710.5418). Role: Nonlinear damping slab λ∥. Reading: full text. Proposed citation key (additions file): **Shalchi-2008-MNRAS**.

**[M68]** Engelbrecht, N. E.; Burger, R. A. (2013). **An Ab Initio Model for Cosmic-Ray Modulation.** *The Astrophysical Journal*, **772**(1), 46. [DOI](https://doi.org/10.1088/0004-637X/772/1/46). Role: Turbulence-driven GCR λ∥, λ⊥. Reading: full text. Citation key in the supplied bibliography: **Engelbrecht-2013-AJ-772**.

**[M69]** Burger, R. A.; Krüger, T. P. J.; Hitge, M.; Engelbrecht, N. E. (2008). **A Fisk-Parker Hybrid Heliospheric Magnetic Field with a Solar-Cycle Dependence.** *The Astrophysical Journal*, **674**(1), 511–519. [DOI](https://doi.org/10.1086/525039). Role: Origin of the continuous proton λ∥ construction. Reading: metadata only. Proposed citation key (additions file): **Burger-2008-AJ**.

**[M70]** Engelbrecht, N. E.; Burger, R. A. (2014). **Cosmic-Ray Modulation: an Ab Initio Approach.** *Brazilian Journal of Physics*, **44**(5), 512–519. [DOI](https://doi.org/10.1007/s13538-014-0241-7). Role: Proton λ∥ construction. Reading: full text. Proposed citation key (additions file): **Engelbrecht-2014-BJP**.

**[M71]** Zank, G. P.; Matthaeus, W. H.; Bieber, J. W.; Moraal, H. (1998). **The radial and latitudinal dependence of the cosmic ray diffusion tensor in the heliosphere.** *Journal of Geophysical Research*, **103**(A2), 2085–2097. [DOI](https://doi.org/10.1029/97JA03013). Role: Slab QLT fit for λ∥. Reading: secondary. Pre-2006 foundation. Citation key in the supplied bibliography: **Zank-1998-JGR**.

**[M72]** Zank, G. P.; Le Roux, J. A.; Matthaeus, W. H.; Moraal, H. (1999). **Solar Wind Turbulence, Diffusion Coefficients, and Cosmic Ray Modulation.** *Proceedings of the 26th International Cosmic Ray Conference, paper SH.3.1.11*, **7**, 41–44. [link](https://dpnc.unige.ch/ams/OLD_WebSite/ICRC-99/root/vol7/s3_1_11.pdf). Role: Printed form of the Zank 1998 fit. Reading: full text. Pre-2006 foundation. Proposed citation key (additions file): **Zank-1999-ICRC**.

**[M73]** Perri, B.; Brun, A. S.; Strugarek, A.; Réville, V. (2020). **Impact of solar magnetic field amplitude and geometry on cosmic rays diffusion coefficients in the inner heliosphere.** *Journal of Space Weather and Space Climate*, **10**, 55. [DOI](https://doi.org/10.1051/swsc/2020057); [arXiv:2010.01880](https://arxiv.org/abs/2010.01880). Role: Zank fit with units. Reading: full text. Citation key in the supplied bibliography: **Perri-2020-JSWSC**.

**[M74]** Quenby, J. J.; Webber, W. R. (2015). **Transient Heliosheath Modulation.** *Monthly Notices of the Royal Astronomical Society*, **453**(2), 1297–1304. [DOI](https://doi.org/10.1093/mnras/stv1482); [arXiv:1409.1105](https://arxiv.org/abs/1409.1105). Role: Le Roux et al. 1999 coefficient. Reading: full text. Proposed citation key (additions file): **Quenby-2015-MNRAS**.

**[M75]** Hussein, M.; Tautz, R. C.; Shalchi, A. (2015). **The influence of different turbulence models on the diffusion coefficients of energetic particles.** *Journal of Geophysical Research: Space Physics*, **120**(6), 4095–4111. [DOI](https://doi.org/10.1002/2015JA021060); [arXiv:1505.05099](https://arxiv.org/abs/1505.05099). Role: Test-particle λ∥ benchmarks. Reading: full text. Citation key in the supplied bibliography: **Hussein-2015-JGRSP**.

**[M76]** Engelbrecht, N. E.; Strauss, R. D. (2018). **A Tractable Estimate for the Dissipation Range Onset Wavenumber Throughout the Heliosphere.** *The Astrophysical Journal*, **856**(2), 159. [DOI](https://doi.org/10.3847/1538-4357/aab495). Role: Dissipation-range onset k_d. Reading: full text. Citation key in the supplied bibliography: **Engelbrecht-2018-AJ**.

**[M77]** Kozarev, K. A.; Schwadron, N. A. (2016). **A Data-Driven Analytic Model for Proton Acceleration by Large-Scale Solar Coronal Shocks.** *The Astrophysical Journal*, **831**(2), 120. [DOI](https://doi.org/10.3847/0004-637X/831/2/120); [arXiv:1608.00240](https://arxiv.org/abs/1608.00240). Role: Hard-sphere κ tensor at shocks. Reading: full text. Citation key in the supplied bibliography: **Kozarev-2016-AJ**.

**[M78]** Battarbee, M.; Laitinen, T.; Vainio, R.; Agueda, N. (2010). **Acceleration of Energetic Particles Through Self-Generated Waves in a Decelerating Coronal Shock.** *AIP Conference Proceedings (Solar Wind 12)*, **1216**, 84–87. [DOI](https://doi.org/10.1063/1.3395969); [arXiv:1303.4334](https://arxiv.org/abs/1303.4334). Role: Ambient coronal λ0 for self-generated-wave DSA. Reading: full text. Citation key in the supplied bibliography: **Battarbee-2010-AIP**.

**[M79]** Ng, C. K.; Reames, D. V. (2008). **Shock Acceleration of Solar Energetic Protons: The First 10 Minutes.** *The Astrophysical Journal Letters*, **686**, L123–L126. [DOI](https://doi.org/10.1086/592996); [link](https://ntrs.nasa.gov/citations/20080038680). Role: Wave saturation cap and ambient λ. Reading: full text. Citation key in the supplied bibliography: **Ng-2008-AJL, Ng-2008-NASA**.

**[M80]** Ng, C. K.; Reames, D. V.; Tylka, A. J. (2012). **Solar Energetic Particles: Shock Acceleration and Transport through Self-Amplified Waves.** *AIP Conference Proceedings*, **1436**, 212–218. [DOI](https://doi.org/10.1063/1.4723610). Role: Self-amplified wave transport. Reading: full text. Citation key in the supplied bibliography: **Ng-2012-AIP**. Metadata note: unconfirmed field(s): volume.

**[M81]** Hu, J.; Li, G.; Ao, X.; Zank, G. P.; Verkhoglyadova, O. (2017). **Modeling Particle Acceleration and Transport at a 2-D CME-Driven Shock.** *Journal of Geophysical Research: Space Physics*, **122**(11). [DOI](https://doi.org/10.1002/2017JA024077). Role: iPATH shock module. Reading: abstract only. Citation key in the supplied bibliography: **Hu-2017-JGR**. Metadata note: unconfirmed field(s): pages.

**[M82]** Ding, Z.-Y.; Li, G.; Hu, J.-X.; Fu, S. (2020). **Modeling the 2017 September 10 solar energetic particle event using the iPATH model.** *Research in Astronomy and Astrophysics*, **20**(9), 145. [DOI](https://doi.org/10.1088/1674-4527/20/9/145); [arXiv:2005.02326](https://arxiv.org/abs/2005.02326). Role: iPATH κ scalings. Reading: full text. Citation key in the supplied bibliography: **Ding-2020-RAA, Ding-2020-arXiv**.

**[M83]** Luhmann, J. G.; Ledvina, S. A.; Odstrcil, D.; Owens, M. J.; Zhao, X. P.; Liu, Y.; Riley, P. (2010). **Cone model-based SEP event calculations for applications to multipoint observations.** *Advances in Space Research*, **46**(1), 1–21. [DOI](https://doi.org/10.1016/j.asr.2010.03.011). Role: Scatter-free SEPMOD transport. Reading: full text. Citation key in the supplied bibliography: **Luhmann-2010-ASR**.

**[M84]** Potgieter, M. S. (2013). **Solar Modulation of Cosmic Rays.** *Living Reviews in Solar Physics*, **10**, 3. [DOI](https://doi.org/10.12942/lrsp-2013-3); [arXiv:1306.4421](https://arxiv.org/abs/1306.4421). Role: Diffusion tensor and NWU broken power law. Reading: full text. Citation key in the supplied bibliography: **Potgieter-2013-LRSP**.

**[M85]** Potgieter, M. S.; Vos, E. E.; Boezio, M.; De Simone, N.; Di Felice, V.; Formato, V. (2014). **Modulation of Galactic Protons in the Heliosphere During the Unusual Solar Minimum of 2006 to 2009.** *Solar Physics*, **289**(1), 391–406. [DOI](https://doi.org/10.1007/s11207-013-0324-6); [arXiv:1302.1284](https://arxiv.org/abs/1302.1284). Role: NWU K∥ form and 2006-2009 proton parameters. Reading: full text. Proposed citation key (additions file): **Potgieter-2014-SP**.

**[M86]** Potgieter, M. S.; Vos, E. E.; Munini, R.; Boezio, M.; Di Felice, V. (2015). **Modulation of Galactic Electrons in the Heliosphere during the Unusual Solar Minimum of 2006-2009: A Modeling Approach.** *The Astrophysical Journal*, **810**(2), 141. [DOI](https://doi.org/10.1088/0004-637X/810/2/141). Role: Electron parameters 2006-2009. Reading: full text. Proposed citation key (additions file): **Potgieter-2015-AJ**.

**[M87]** Vos, E. E.; Potgieter, M. S. (2015). **New Modeling of Galactic Proton Modulation during the Minimum of Solar Cycle 23/24.** *The Astrophysical Journal*, **815**(2), 119. [DOI](https://doi.org/10.1088/0004-637X/815/2/119). Role: Proton λ∥(1 GV) and shape parameters 2006-2009. Reading: full text. Citation key in the supplied bibliography: **Vos-2015-AJ**.

**[M88]** Vos, E. E.; Potgieter, M. S. (2016). **Global gradients for cosmic-ray protons in the heliosphere during the solar minimum of cycle 23/24.** *Solar Physics*, **291**(7), 2181–2195. [DOI](https://doi.org/10.1007/s11207-016-0945-7); [arXiv:1608.01688](https://arxiv.org/abs/1608.01688). Role: Printed exponent variant. Reading: full text. Proposed citation key (additions file): **Vos-2016-SP**.

**[M89]** Aslam, O. P. M.; Bisschoff, D.; Potgieter, M. S.; Boezio, M.; Munini, R. (2019). **Modeling of Heliospheric Modulation of Cosmic-Ray Positrons in a Very Quiet Heliosphere.** *The Astrophysical Journal*, **873**(1), 70. [DOI](https://doi.org/10.3847/1538-4357/ab05e6); [arXiv:1811.10710](https://arxiv.org/abs/1811.10710). Role: Positron parameters 2006-2009. Reading: full text. Proposed citation key (additions file): **Aslam-2019-AJ**.

**[M90]** Aslam, O. P. M.; Luo, X.; Potgieter, M. S.; Ngobeni, M. D.; Song, X. (2023). **Unfolding Drift Effects for Cosmic Rays over the Period of the Sun's Magnetic Field Reversal.** *The Astrophysical Journal*, **947**(2), 72. [DOI](https://doi.org/10.3847/1538-4357/acc24a); [arXiv:2212.13397](https://arxiv.org/abs/2212.13397). Role: (K∥)0 values in 10^22 units across the reversal. Reading: full text. Proposed citation key (additions file): **Aslam-2023-AJ**.

**[M91]** Aslam, O. P. M.; Bisschoff, D.; Ngobeni, M. D.; Potgieter, M. S.; Munini, R.; Boezio, M.; Mikhailov, V. V. (2021). **Time and Charge-sign Dependence of the Heliospheric Modulation of Cosmic Rays.** *The Astrophysical Journal*, **909**(2), 215. [DOI](https://doi.org/10.3847/1538-4357/abdd35); [arXiv:2011.02052](https://arxiv.org/abs/2011.02052). Role: (K∥)0 values in 6×10^20 units. Reading: full text. Proposed citation key (additions file): **Aslam-2021-AJ**.

**[M92]** Raath, J. L. (2015). **A comparative study of cosmic ray modulation models.** *MSc (Magister Scientiae) dissertation, North-West University*. [link](https://hdl.handle.net/10394/15516). Role: 6×10^20 unit convention. Reading: full text. Proposed citation key (additions file): **Raath-2015-MSc**.

**[M93]** Corti, C.; Potgieter, M. S.; Bindi, V.; Consolandi, C.; Light, C.; Palermo, M.; Popkow, A. (2019). **Numerical Modeling of Galactic Cosmic-Ray Proton and Helium Observed by AMS-02 during the Solar Maximum of Solar Cycle 24.** *The Astrophysical Journal*, **871**(2), 253. [DOI](https://doi.org/10.3847/1538-4357/aafac4); [arXiv:1810.09640](https://arxiv.org/abs/1810.09640). Role: R/R_k broken power law and per-BR parameters. Reading: full text. Citation key in the supplied bibliography: **Corti-2019-AJ**.

**[M94]** Luo, X.; Potgieter, M. S.; Bindi, V.; Zhang, M.; Feng, X. (2019). **A Numerical Study of Cosmic Proton Modulation Using AMS-02 Observations.** *The Astrophysical Journal*, **878**(1), 6. [DOI](https://doi.org/10.3847/1538-4357/ab1b2a). Role: Transient-event parameters (10^20 units). Reading: full text. Proposed citation key (additions file): **Luo-2019-AJ**.

**[M95]** Boschini, M. J.; Della Torre, S.; Gervasi, M.; La Vacca, G.; Rancoita, P. G. (2018). **Propagation of Cosmic Rays in Heliosphere: The HELMOD Model.** *Advances in Space Research*, **62**, 2859–2879. [DOI](https://doi.org/10.1016/j.asr.2017.04.017); [arXiv:1704.03733](https://arxiv.org/abs/1704.03733). Role: HelMod K0(SSN), K0(NMCR) coefficients. Reading: full text. Citation key in the supplied bibliography: **Boschini-2018-ASR**.

**[M96]** Della Torre, S.; Gervasi, M.; Grandi, D.; Jóhannesson, G.; La Vacca, G.; Masi, N.; Moskalenko, I. V.; Orlando, E.; Porter, T. A.; Quadrani, L.; Rancoita, P. G.; Rozza, D. (2016). **HelMod: A Comprehensive Treatment of the Cosmic Ray Transport Through the Heliosphere.** *XXV European Cosmic Ray Symposium (ECRS 2016) Proceedings, eConf C16-09-04.3*. [arXiv:1612.08445](https://arxiv.org/abs/1612.08445). Role: HelMod K∥ form. Reading: full text. Proposed citation key (additions file): **DellaTorre-2016-ECRS**.

**[M97]** Boschini, M. J.; Della Torre, S.; Gervasi, M.; La Vacca, G.; Rancoita, P. G. (2019). **The HelMod Model in the Works for Inner and Outer Heliosphere: from AMS to Voyager Probes Observations.** *Advances in Space Research*, **64**(12), 2459–2476. [DOI](https://doi.org/10.1016/j.asr.2019.04.007); [arXiv:1903.07501](https://arxiv.org/abs/1903.07501). Role: HelMod v4 K∥ and transition function. Reading: full text. Proposed citation key (additions file): **Boschini-2019-ASR**.

**[M98]** Boschini, M. J.; Della Torre, S.; Gervasi, M.; Grandi, D.; Jóhannesson, G.; Kachelriess, M.; La Vacca, G.; Masi, N.; Moskalenko, I. V.; Orlando, E.; Ostapchenko, S. S.; Pensotti, S.; Porter, T. A.; Quadrani, L.; Rancoita, P. G.; Rozza, D.; Tacconi, M. (2017). **Solution of Heliospheric Propagation: Unveiling the Local Interstellar Spectra of Cosmic-ray Species.** *The Astrophysical Journal*, **840**(2), 115. [DOI](https://doi.org/10.3847/1538-4357/aa6e4f); [arXiv:1704.06337](https://arxiv.org/abs/1704.06337). Role: HelMod g_low and ρ. Reading: full text. Citation key in the supplied bibliography: **Boschini-2017-AJ**.

**[M99]** Bobik, P.; Boella, G.; Boschini, M. J.; Consolandi, C.; Della Torre, S.; Gervasi, M.; Grandi, D.; Kudela, K.; Pensotti, S.; Rancoita, P. G.; Tacconi, M. (2012). **Systematic Investigation of Solar Modulation of Galactic Protons for Solar Cycle 23 Using a Monte Carlo Approach with Particle Drift Effects and Latitudinal Dependence.** *The Astrophysical Journal*, **745**(2), 132. [DOI](https://doi.org/10.1088/0004-637X/745/2/132); [arXiv:1110.4315](https://arxiv.org/abs/1110.4315). Role: Early HelMod K∥ and K_F(SSN). Reading: full text. Citation key in the supplied bibliography: **Bobik-2012-AJ**.

**[M100]** Strauss, R. D.; Potgieter, M. S.; Büsching, I.; Kopp, A. (2011). **Modeling the Modulation of Galactic and Jovian Electrons by Stochastic Processes.** *The Astrophysical Journal*, **735**(2), 83. [DOI](https://doi.org/10.1088/0004-637X/735/2/83). Role: λ0 = 0.15 AU electron prescription. Reading: full text. Citation key in the supplied bibliography: **Strauss-2011-AJ**.

**[M101]** Effenberger, F.; Fichtner, H.; Scherer, K.; Barra, S.; Kleimann, J.; Strauss, R. D. (2012). **A Generalized Diffusion Tensor for Fully Anisotropic Diffusion of Energetic Particles in the Heliospheric Magnetic Field.** *The Astrophysical Journal*, **750**(2), 108. [DOI](https://doi.org/10.1088/0004-637X/750/2/108); [arXiv:1202.6319](https://arxiv.org/abs/1202.6319). Role: κ∥0 = 0.9×10^22 cm²/s prescription. Reading: full text. Proposed citation key (additions file): **Effenberger-2012-AJ**.

**[M102]** Wang, B.-B.; Bi, X.-J.; Fang, K.; Lin, S.-J.; Yin, P.-F. (2019). **Time-dependent solar modulation of cosmic rays from solar minimum to solar maximum.** *Physical Review D*, **100**(6), 063006. [DOI](https://doi.org/10.1103/PhysRevD.100.063006); [arXiv:1904.03747](https://arxiv.org/abs/1904.03747). Role: Two-branch K∥. Reading: full text. Proposed citation key (additions file): **Wang-2019-PRD**.

**[M103]** Tomassetti, N. (2017). **Solar and nuclear physics uncertainties in cosmic-ray propagation.** *Physical Review D*, **96**(10), 103005. [DOI](https://doi.org/10.1103/PhysRevD.96.103005); [arXiv:1707.06917](https://arxiv.org/abs/1707.06917). Role: κ0 = a/φ + b normalization. Reading: full text. Proposed citation key (additions file): **Tomassetti-2017-PRD**.

**[M104]** Fiandrini, E.; Tomassetti, N.; Bertucci, B.; Donnini, F.; Graziani, M.; Khiali, B.; Reina Conde, A. (2021). **Numerical modeling of cosmic rays in the heliosphere: Analysis of proton data from AMS-02 and PAMELA.** *Physical Review D*, **104**(2), 023012. [DOI](https://doi.org/10.1103/PhysRevD.104.023012); [arXiv:2010.08649](https://arxiv.org/abs/2010.08649). Role: Perugia broken power law. Reading: full text. Citation key in the supplied bibliography: **Fiandrini-2021-PRD**.

**[M105]** Tomassetti, N.; Bertucci, B.; Fiandrini, E.; Khiali, B. (2025). **Propagation Times and Energy Losses of Cosmic Protons and Antiprotons in Interplanetary Space.** *Galaxies*, **13**(2), 23. [DOI](https://doi.org/10.3390/galaxies13020023); [arXiv:2503.14025](https://arxiv.org/abs/2503.14025). Role: Perugia h ≡ 0.01 variant. Reading: full text. Proposed citation key (additions file): **Tomassetti-2025-Galaxies**.

**[M106]** Jiang, J.; Lin, S.; Yang, L. (2023). **A New Scenario of Solar Modulation Model during the Polarity Reversing.** *The Astrophysical Journal*, **957**(2), 72. [DOI](https://doi.org/10.3847/1538-4357/acf719); [arXiv:2303.04460](https://arxiv.org/abs/2303.04460). Role: Broken power law with averages. Reading: full text. Proposed citation key (additions file): **Jiang-2023-AJ**.

**[M107]** Duan, K.-K.; Wang, X.; Li, W.-H.; Xu, Z.-H.; Tsai, Y.-L. S.; Fan, Y.-Z. (2025). **Scrutinizing the impact of the solar modulation on AMS-02 antiproton excess.** *Journal of Cosmology and Astroparticle Physics*, **2025**(10), 049. [DOI](https://doi.org/10.1088/1475-7516/2025/10/049); [arXiv:2506.13352](https://arxiv.org/abs/2506.13352). Role: Alternative broken power law. Reading: full text. Proposed citation key (additions file): **Duan-2025-JCAP**.

**[M108]** Qin, G.; Shen, Z.-N. (2017). **Modulation of Galactic Cosmic Rays in the Inner Heliosphere, Comparing with PAMELA Measurements.** *The Astrophysical Journal*, **846**(1), 56. [DOI](https://doi.org/10.3847/1538-4357/aa83ad); [arXiv:1705.04847](https://arxiv.org/abs/1705.04847). Role: NLGCE-F in a modulation model. Reading: full text. Citation key in the supplied bibliography: **Qin-2017-AJ**.

**[M109]** Shen, Z.-N.; Qin, G. (2018). **Modulation of Galactic Cosmic Rays in the Inner Heliosphere over Solar Cycles.** *The Astrophysical Journal*, **854**(2), 137. [DOI](https://doi.org/10.3847/1538-4357/aaab64); [arXiv:1709.08017](https://arxiv.org/abs/1709.08017). Role: NLGCE-F modulation over cycles. Reading: full text. Citation key in the supplied bibliography: **Shen-2018-AJ**.

**[M110]** Tomassetti, N.; Bertucci, B.; Donnini, F.; Graziani, M.; Fiandrini, E.; Khiali, B.; Reina Conde, A. (2023). **Data driven analysis of cosmic rays in the heliosphere: diffusion of cosmic protons.** *Rendiconti Lincei. Scienze Fisiche e Naturali*, **34**(2), 333–338. [DOI](https://doi.org/10.1007/s12210-023-01149-1); [arXiv:2303.12239](https://arxiv.org/abs/2303.12239). Role: Perugia λ∥ printing. Reading: full text. Proposed citation key (additions file): **Tomassetti-2023-RL**.

**[M111]** Di Felice, V.; Munini, R.; Vos, E. E.; Potgieter, M. S. (2017). **New Evidence for Charge-sign-dependent Modulation during the Solar Minimum of 2006 to 2009.** *The Astrophysical Journal*, **834**(1), 89. [DOI](https://doi.org/10.3847/1538-4357/834/1/89); [arXiv:1608.01301](https://arxiv.org/abs/1608.01301). Role: Perpendicular ratios and drift reduction. Reading: full text. Proposed citation key (additions file): **DiFelice-2017-AJ**.

**[M112]** MacBride, B. T.; Smith, C. W.; Vasquez, B. J. (2010). **Inertial-range anisotropies in the solar wind from 0.3 to 1 AU: Helios 1 observations.** *Journal of Geophysical Research*, **115**(A7), A07105. [DOI](https://doi.org/10.1029/2009JA014939). Role: Slab energy fraction 0.3-1 au. Reading: full text. Proposed citation key (additions file): **MacBride-2010-JGR**. Metadata note: unconfirmed field(s): pages.

**[M113]** Fa, Z.; He, H.-Q. (2026). **Solar-cycle variability of composite geometry in solar wind turbulence.** *Astronomy & Astrophysics*, **712**, A86. [DOI](https://doi.org/10.1051/0004-6361/202556095); [arXiv:2505.12870](https://arxiv.org/abs/2505.12870). Role: Slab fraction 1995-2023. Reading: full text. Proposed citation key (additions file): **Fa-2026-AA**.

**[M114]** Zhao, L.-L.; Zank, G. P.; Adhikari, L.; Nakanotani, M. (2022). **Inertial-range Magnetic-fluctuation Anisotropy Observed from Parker Solar Probe's First Seven Orbits.** *The Astrophysical Journal Letters*, **924**(1), L5. [DOI](https://doi.org/10.3847/2041-8213/ac4415); [arXiv:2112.01711](https://arxiv.org/abs/2112.01711). Role: Near-Sun slab fraction. Reading: full text. Proposed citation key (additions file): **Zhao-2022-AJL**.

**[M115]** Cheng, W.; Xiong, M.; Jiao, Y.; Ran, H.; Yang, L.; Hu, H.; Wang, R. (2025). **Inertial-range Turbulence Anisotropy of the Young Solar Wind from Different Source Regions.** *The Astrophysical Journal Letters*, **988**(1), L15. [DOI](https://doi.org/10.3847/2041-8213/adeb8a); [arXiv:2507.04288](https://arxiv.org/abs/2507.04288). Role: Near-Sun slab fraction by source region. Reading: full text. Proposed citation key (additions file): **Cheng-2025-AJL**.

**[M116]** Zhao, L.-L.; Adhikari, L.; Zank, G. P.; Hu, Q.; Feng, X. S. (2018). **Influence of the Solar Cycle on Turbulence Properties and Cosmic-Ray Diffusion.** *The Astrophysical Journal*, **856**(2), 94. [DOI](https://doi.org/10.3847/1538-4357/aab362). Role: Solar-cycle turbulence and λ∥ at 1 au. Reading: full text. Citation key in the supplied bibliography: **Zhao-2018-AJ**.

**[M117]** Oughton, S.; Matthaeus, W. H.; Smith, C. W.; Breech, B.; Isenberg, P. A. (2011). **Transport of solar wind fluctuations: A two-component model.** *Journal of Geophysical Research*, **116**(A8), A08105. [DOI](https://doi.org/10.1029/2010JA016365). Role: Two-component turbulence transport. Reading: full text. Citation key in the supplied bibliography: **Oughton-2011-JGR**. Metadata note: unconfirmed field(s): pages.

**[M118]** Adhikari, L.; Zank, G. P.; Hunana, P.; Shiota, D.; Bruno, R.; Hu, Q.; Telloni, D. (2017). **II. Transport of Nearly Incompressible Magnetohydrodynamic Turbulence from 1 to 75 au.** *The Astrophysical Journal*, **841**(2), 85. [DOI](https://doi.org/10.3847/1538-4357/aa6f5d). Role: NI-MHD boundary values and radial exponents. Reading: full text. Proposed citation key (additions file): **Adhikari-2017-AJ**.

**[M119]** Adhikari, L.; Zank, G. P.; Zhao, L.-L.; Kasper, J. C.; Korreck, K. E.; Stevens, M.; Case, A. W.; Whittlesey, P.; Larson, D.; Livi, R.; Klein, K. G. (2020). **Turbulence Transport Modeling and First Orbit Parker Solar Probe (PSP) Observations.** *The Astrophysical Journal Supplement Series*, **246**(2), 38. [DOI](https://doi.org/10.3847/1538-4365/ab5852). Role: Inner-heliosphere turbulence model. Reading: full text. Proposed citation key (additions file): **Adhikari-2020-AJSS**.

**[M120]** Wiengarten, T.; Oughton, S.; Engelbrecht, N. E.; Fichtner, H.; Kleimann, J.; Scherer, K. (2016). **A Generalized Two-component Model of Solar Wind Turbulence and ab initio Diffusion Mean-free Paths and Drift Lengthscales of Cosmic Rays.** *The Astrophysical Journal*, **833**(1), 17. [DOI](https://doi.org/10.3847/0004-637X/833/1/17); [arXiv:1609.08271](https://arxiv.org/abs/1609.08271). Role: Two-component model and δB² conversion. Reading: full text. Proposed citation key (additions file): **Wiengarten-2016-AJ**.

**[M121]** Bieber, J. W.; Wanner, W.; Matthaeus, W. H. (1996). **Dominant two-dimensional solar wind turbulence with implications for cosmic ray transport.** *Journal of Geophysical Research*, **101**(A2), 2511–2522. [DOI](https://doi.org/10.1029/95JA02588). Role: Slab/2D composite geometry. Reading: abstract only. Pre-2006 foundation. Citation key in the supplied bibliography: **Bieber-1996-JGR**.

**[M122]** Leamon, R. J.; Matthaeus, W. H.; Smith, C. W.; Zank, G. P.; Mullan, D. J.; Oughton, S. (2000). **MHD-driven Kinetic Dissipation in the Solar Wind and Corona.** *The Astrophysical Journal*, **537**(2), 1054–1062. [DOI](https://doi.org/10.1086/309059). Role: Dissipation-range onset fits. Reading: full text. Pre-2006 foundation. Proposed citation key (additions file): **Leamon-2000-AJ**.

**[M123]** Zank, G. P.; Adhikari, L.; Hunana, P.; Shiota, D.; Bruno, R.; Telloni, D. (2017). **Theory and Transport of Nearly Incompressible Magnetohydrodynamic Turbulence.** *The Astrophysical Journal*, **835**(2), 147. [DOI](https://doi.org/10.3847/1538-4357/835/2/147). Role: NI-MHD turbulence transport. Reading: full text. Proposed citation key (additions file): **Zank-2017-AJ**.

**[M124]** Zank, G. P.; Nakanotani, M.; Zhao, L.-L.; Adhikari, L.; Telloni, D. (2020). **Spectral Anisotropy in 2D plus Slab Magnetohydrodynamic Turbulence in the Solar Wind and Upper Corona.** *The Astrophysical Journal*, **900**(2), 115. [DOI](https://doi.org/10.3847/1538-4357/abad30). Role: Slab and 2D spectral indices. Reading: full text. Proposed citation key (additions file): **Zank-2020-AJ**.

**[M125]** Chhiber, R. (2022). **Anisotropic Magnetic Turbulence in the Inner Heliosphere—Radial Evolution of Distributions Observed by Parker Solar Probe.** *The Astrophysical Journal*, **939**(1), 33. [DOI](https://doi.org/10.3847/1538-4357/ac9386); [arXiv:2205.14096](https://arxiv.org/abs/2205.14096). Role: δb/B0 radial scaling. Reading: full text. Proposed citation key (additions file): **Chhiber-2022-AJ**.

**[M126]** Zank, G. P.; Matthaeus, W. H.; Smith, C. W. (1996). **Evolution of turbulent magnetic fluctuation power with heliospheric distance.** *Journal of Geophysical Research*, **101**(A8), 17093–17107. [DOI](https://doi.org/10.1029/96JA01275). Role: Radial exponents of δB². Reading: abstract only. Pre-2006 foundation. Proposed citation key (additions file): **Zank-1996-JGR**.

**[M127]** Smith, C. W.; Matthaeus, W. H.; Zank, G. P.; Ness, N. F.; Oughton, S.; Richardson, J. D. (2001). **Heating of the low-latitude solar wind by dissipation of turbulent magnetic fluctuations.** *Journal of Geophysical Research*, **106**(A5), 8253–8272. [DOI](https://doi.org/10.1029/2000JA000366). Role: 1 AU turbulence energy and boundary values. Reading: full text. Pre-2006 foundation. Citation key in the supplied bibliography: **Smith-2001-JGR**.

**[M128]** Breech, B.; Matthaeus, W. H.; Minnie, J.; Bieber, J. W.; Oughton, S.; Smith, C. W.; Isenberg, P. A. (2008). **Turbulence transport throughout the heliosphere.** *Journal of Geophysical Research*, **113**(A8), A08105. [DOI](https://doi.org/10.1029/2007JA012711). Role: Single-component turbulence transport. Reading: full text. Citation key in the supplied bibliography: **Breech-2008-JGR**. Metadata note: unconfirmed field(s): pages.

**[M129]** Usmanov, A. V.; Matthaeus, W. H.; Goldstein, M. L.; Chhiber, R. (2018). **The Steady Global Corona and Solar Wind: A Three-dimensional MHD Simulation with Turbulence Transport and Heating.** *The Astrophysical Journal*, **865**(1), 25. [DOI](https://doi.org/10.3847/1538-4357/aad687). Role: Global turbulence model. Reading: full text. Proposed citation key (additions file): **Usmanov-2018-AJ**.

**[M130]** Cuesta, M. E.; Chhiber, R.; Roy, S.; Goodwill, J.; Pecora, F.; Jarosik, J.; Matthaeus, W. H.; Parashar, T. N.; Bandyopadhyay, R. (2022). **Isotropization and Evolution of Energy-containing Eddies in Solar Wind Turbulence: Parker Solar Probe, Helios 1, ACE, WIND, and Voyager 1.** *The Astrophysical Journal Letters*, **932**(1), L11. [DOI](https://doi.org/10.3847/2041-8213/ac73fd); [arXiv:2205.00526](https://arxiv.org/abs/2205.00526). Role: Correlation lengths and radial exponents. Reading: full text. Proposed citation key (additions file): **Cuesta-2022-AJL**.

**[M131]** Ruiz, M. E.; Dasso, S.; Matthaeus, W. H.; Weygand, J. M. (2014). **Characterization of the Turbulent Magnetic Integral Length in the Solar Wind: From 0.3 to 5 Astronomical Units.** *Solar Physics*, **289**(10), 3917–3933. [DOI](https://doi.org/10.1007/s11207-014-0531-9); [arXiv:1404.2826](https://arxiv.org/abs/1404.2826). Role: Integral scale versus distance. Reading: full text. Proposed citation key (additions file): **Ruiz-2014-SP**.

**[M132]** Weygand, J. M.; Matthaeus, W. H.; Dasso, S.; Kivelson, M. G. (2011). **Correlation and Taylor scale variability in the interplanetary magnetic field fluctuations as a function of solar wind speed.** *Journal of Geophysical Research*, **116**(A8), A08102. [DOI](https://doi.org/10.1029/2011JA016621). Role: Parallel/perpendicular correlation scales. Reading: full text. Proposed citation key (additions file): **Weygand-2011-JGR**. Metadata note: unconfirmed field(s): pages.

**[M133]** Chen, C. H. K.; Bale, S. D.; Bonnell, J. W.; Borovikov, D.; Bowen, T. A.; Burgess, D.; Case, A. W.; Chandran, B. D. G.; Dudok de Wit, T.; Goetz, K.; Harvey, P. R.; Kasper, J. C.; Klein, K. G.; Korreck, K. E.; Larson, D.; Livi, R.; MacDowall, R. J.; Malaspina, D. M.; Mallet, A.; McManus, M. D.; Moncuquet, M.; Pulupa, M.; Stevens, M. L.; Whittlesey, P. (2020). **The Evolution and Role of Solar Wind Turbulence in the Inner Heliosphere.** *The Astrophysical Journal Supplement Series*, **246**(2), 53. [DOI](https://doi.org/10.3847/1538-4365/ab60a3); [arXiv:1912.02348](https://arxiv.org/abs/1912.02348). Role: Inner-heliosphere spectral index. Reading: full text. Citation key in the supplied bibliography: **Chen-2019-AJSS**.

**[M134]** Alexandrova, O.; Saur, J.; Lacombe, C.; Mangeney, A.; Mitchell, J.; Schwartz, S. J.; Robert, P. (2009). **Universality of Solar-Wind Turbulent Spectrum from MHD to Electron Scales.** *Physical Review Letters*, **103**(16), 165003. [DOI](https://doi.org/10.1103/PhysRevLett.103.165003). Role: Ion-scale spectral index. Reading: abstract only. Proposed citation key (additions file): **Alexandrova-2009-PRL**.

**[M135]** Raath, J. L.; Potgieter, M. S.; Strauss, R. D.; Kopp, A. (2016). **The effects of magnetic field modifications on the solar modulation of cosmic rays with a SDE-based model.** *Advances in Space Research*, **57**, 1965–1977. [DOI](https://doi.org/10.1016/j.asr.2016.01.017); [arXiv:1506.07305](https://arxiv.org/abs/1506.07305). Role: Parker, Smith-Bieber and Jokipii-Kóta field models. Reading: full text. Citation key in the supplied bibliography: **Raath-2016-ASR**.

## 18. Bibliography comparison

**Snapshot appendix.** This audit describes the supplied file at the date of this document. It is separable from the model specification.

The comparison uses the supplied bibliography (bibliography.bib), containing 6893 parsed records. Its SHA-256 digest is:

~~~text
592b51b08ebe0e1635f26dc6545bec98c5805d7c1be2d93173408521d9b7ffd2
~~~

A reference counts as present when its DOI, its arXiv identifier, or its normalized title matches a record; candidate matches by first author and year were then checked by title. "Not found" means no matching record in this snapshot; it does not imply absence from every version of the bibliography. Conflicting records are listed in Section 18.3 and are not treated as safe additions.

### 18.1 Cited papers not found in the supplied bibliography

**52 references from 2006–2026 were not found.** Their BibTeX entries are in MEAN_FREE_PATH_MODEL_additions.bib:

| Reference | Paper | Proposed key | Role |
|---|---|---|---|
| M08 | Shalchi et al. (2006) — Parallel and Perpendicular Transport of Heliospheric Cosmic Rays in an Improved Dynamical Turbulence Model | Shalchi-2006-AJ | Dynamical turbulence λ∥, λ⊥ with heliospheric parameters |
| M13 | Agueda et al. (2010) — Solar near-relativistic electron observations as a proof of a back-scatter region beyond 1 AU during the 2000 February 18 event | Agueda-2010-AA | λ_r definition, ε-form D_μμ, back-scatter region |
| M14 | Agueda et al. (2014) — Release timescales of solar energetic particles in the low corona | Agueda-2014-AA | Event-fitted electron λ_r (Table 3) |
| M15 | Agueda et al. (2016) — Release History and Transport Parameters of Relativistic Solar Electrons Inferred from Near-the-Sun In Situ Observations | Agueda-2016-AJ | Helios event-fitted λ_r |
| M16 | Pacheco et al. (2019) — Full inversion of solar relativistic electron events measured by the Helios spacecraft | Pacheco-2019-AA | Helios λ_r versus r (Tables 2-4) |
| M17 | Dröge et al. (2016) — Multi-Spacecraft Observations and Transport Modeling of Energetic Electrons for a Series of Solar Particle Events in August 2010 | Droge-2016-AJ | λ∥ and λ⊥ normalized to 1 au |
| M18 | Agueda et al. (2008) — Injection and Interplanetary Transport of Near-Relativistic Electrons: Modeling the Impulsive Event on 2000 May 1 | Agueda-2008-AJ | Event-fitted λ_r |
| M20 | Battarbee et al. (2018) — Multi-spacecraft observations and transport simulations of solar energetic particles for the May 17th 2012 event | Battarbee-2018-AA | Constant λ = 0.3 au full-orbit model input |
| M24 | Tan et al. (2011) — What Causes Scatter-Free Transport of Non-Relativistic Solar Electrons? | Tan-2011-AJ | Scatter-free to diffusive electron transition |
| M25 | Strauss et al. (2020) — On the shape of SEP electron spectra: The role of interplanetary transport | Strauss-2020-AJ | Dissipation-range onset for electron λ∥ |
| M38 | Marsh et al. (2013) — Drift-induced Perpendicular Transport of Solar Energetic Particles | Marsh-2013-AJ | Isotropic scattering full-orbit model |
| M45 | Strauss et al. (2017) — On the Pulse Shape of Ground-Level Enhancements | Strauss-2017-SP | Radial power-law λ_rr |
| M47 | Dröge et al. (2016) — Multi-spacecraft observations and transport modeling of energetic electrons for a series of solar particle events in August 2010 | Droge-2016-ICRC | Dröge SDE model, Λ⊥ scaling, fitted λ_r |
| M58 | Agueda et al. (2013) — On the parametrization of the energetic-particle pitch-angle diffusion coefficient | Agueda-2013-JSWSC | q-form and ε-form D_μμ normalizations |
| M63 | He et al. (2013) — The dependence of the parallel and perpendicular mean free paths on the rigidity of the solar energetic particles: theoretical model versus observations | He-2013-AA | Focusing-corrected λ∥ |
| M67 | Shalchi et al. (2008) — Non-linear damping of slab modes and cosmic ray transport | Shalchi-2008-MNRAS | Nonlinear damping slab λ∥ |
| M69 | Burger et al. (2008) — A Fisk-Parker Hybrid Heliospheric Magnetic Field with a Solar-Cycle Dependence | Burger-2008-AJ | Origin of the continuous proton λ∥ construction |
| M70 | Engelbrecht et al. (2014) — Cosmic-Ray Modulation: an Ab Initio Approach | Engelbrecht-2014-BJP | Proton λ∥ construction |
| M74 | Quenby et al. (2015) — Transient Heliosheath Modulation | Quenby-2015-MNRAS | Le Roux et al. 1999 coefficient |
| M85 | Potgieter et al. (2014) — Modulation of Galactic Protons in the Heliosphere During the Unusual Solar Minimum of 2006 to 2009 | Potgieter-2014-SP | NWU K∥ form and 2006-2009 proton parameters |
| M86 | Potgieter et al. (2015) — Modulation of Galactic Electrons in the Heliosphere during the Unusual Solar Minimum of 2006-2009: A Modeling Approach | Potgieter-2015-AJ | Electron parameters 2006-2009 |
| M88 | Vos et al. (2016) — Global gradients for cosmic-ray protons in the heliosphere during the solar minimum of cycle 23/24 | Vos-2016-SP | Printed exponent variant |
| M89 | Aslam et al. (2019) — Modeling of Heliospheric Modulation of Cosmic-Ray Positrons in a Very Quiet Heliosphere | Aslam-2019-AJ | Positron parameters 2006-2009 |
| M90 | Aslam et al. (2023) — Unfolding Drift Effects for Cosmic Rays over the Period of the Sun's Magnetic Field Reversal | Aslam-2023-AJ | (K∥)0 values in 10^22 units across the reversal |
| M91 | Aslam et al. (2021) — Time and Charge-sign Dependence of the Heliospheric Modulation of Cosmic Rays | Aslam-2021-AJ | (K∥)0 values in 6×10^20 units |
| M92 | Raath (2015) — A comparative study of cosmic ray modulation models | Raath-2015-MSc | 6×10^20 unit convention |
| M94 | Luo et al. (2019) — A Numerical Study of Cosmic Proton Modulation Using AMS-02 Observations | Luo-2019-AJ | Transient-event parameters (10^20 units) |
| M96 | Della Torre et al. (2016) — HelMod: A Comprehensive Treatment of the Cosmic Ray Transport Through the Heliosphere | DellaTorre-2016-ECRS | HelMod K∥ form |
| M97 | Boschini et al. (2019) — The HelMod Model in the Works for Inner and Outer Heliosphere: from AMS to Voyager Probes Observations | Boschini-2019-ASR | HelMod v4 K∥ and transition function |
| M101 | Effenberger et al. (2012) — A Generalized Diffusion Tensor for Fully Anisotropic Diffusion of Energetic Particles in the Heliospheric Magnetic Field | Effenberger-2012-AJ | κ∥0 = 0.9×10^22 cm²/s prescription |
| M102 | Wang et al. (2019) — Time-dependent solar modulation of cosmic rays from solar minimum to solar maximum | Wang-2019-PRD | Two-branch K∥ |
| M103 | Tomassetti (2017) — Solar and nuclear physics uncertainties in cosmic-ray propagation | Tomassetti-2017-PRD | κ0 = a/φ + b normalization |
| M105 | Tomassetti et al. (2025) — Propagation Times and Energy Losses of Cosmic Protons and Antiprotons in Interplanetary Space | Tomassetti-2025-Galaxies | Perugia h ≡ 0.01 variant |
| M106 | Jiang et al. (2023) — A New Scenario of Solar Modulation Model during the Polarity Reversing | Jiang-2023-AJ | Broken power law with averages |
| M107 | Duan et al. (2025) — Scrutinizing the impact of the solar modulation on AMS-02 antiproton excess | Duan-2025-JCAP | Alternative broken power law |
| M110 | Tomassetti et al. (2023) — Data driven analysis of cosmic rays in the heliosphere: diffusion of cosmic protons | Tomassetti-2023-RL | Perugia λ∥ printing |
| M111 | Di Felice et al. (2017) — New Evidence for Charge-sign-dependent Modulation during the Solar Minimum of 2006 to 2009 | DiFelice-2017-AJ | Perpendicular ratios and drift reduction |
| M112 | MacBride et al. (2010) — Inertial-range anisotropies in the solar wind from 0.3 to 1 AU: Helios 1 observations | MacBride-2010-JGR | Slab energy fraction 0.3-1 au |
| M113 | Fa et al. (2026) — Solar-cycle variability of composite geometry in solar wind turbulence | Fa-2026-AA | Slab fraction 1995-2023 |
| M114 | Zhao et al. (2022) — Inertial-range Magnetic-fluctuation Anisotropy Observed from Parker Solar Probe's First Seven Orbits | Zhao-2022-AJL | Near-Sun slab fraction |
| M115 | Cheng et al. (2025) — Inertial-range Turbulence Anisotropy of the Young Solar Wind from Different Source Regions | Cheng-2025-AJL | Near-Sun slab fraction by source region |
| M118 | Adhikari et al. (2017) — II. Transport of Nearly Incompressible Magnetohydrodynamic Turbulence from 1 to 75 au | Adhikari-2017-AJ | NI-MHD boundary values and radial exponents |
| M119 | Adhikari et al. (2020) — Turbulence Transport Modeling and First Orbit Parker Solar Probe (PSP) Observations | Adhikari-2020-AJSS | Inner-heliosphere turbulence model |
| M120 | Wiengarten et al. (2016) — A Generalized Two-component Model of Solar Wind Turbulence and ab initio Diffusion Mean-free Paths and Drift Lengthscales of Cosmic Rays | Wiengarten-2016-AJ | Two-component model and δB² conversion |
| M123 | Zank et al. (2017) — Theory and Transport of Nearly Incompressible Magnetohydrodynamic Turbulence | Zank-2017-AJ | NI-MHD turbulence transport |
| M124 | Zank et al. (2020) — Spectral Anisotropy in 2D plus Slab Magnetohydrodynamic Turbulence in the Solar Wind and Upper Corona | Zank-2020-AJ | Slab and 2D spectral indices |
| M125 | Chhiber (2022) — Anisotropic Magnetic Turbulence in the Inner Heliosphere—Radial Evolution of Distributions Observed by Parker Solar Probe | Chhiber-2022-AJ | δb/B0 radial scaling |
| M129 | Usmanov et al. (2018) — The Steady Global Corona and Solar Wind: A Three-dimensional MHD Simulation with Turbulence Transport and Heating | Usmanov-2018-AJ | Global turbulence model |
| M130 | Cuesta et al. (2022) — Isotropization and Evolution of Energy-containing Eddies in Solar Wind Turbulence: Parker Solar Probe, Helios 1, ACE, WIND, and Voyager 1 | Cuesta-2022-AJL | Correlation lengths and radial exponents |
| M131 | Ruiz et al. (2014) — Characterization of the Turbulent Magnetic Integral Length in the Solar Wind: From 0.3 to 5 Astronomical Units | Ruiz-2014-SP | Integral scale versus distance |
| M132 | Weygand et al. (2011) — Correlation and Taylor scale variability in the interplanetary magnetic field fluctuations as a function of solar wind speed | Weygand-2011-JGR | Parallel/perpendicular correlation scales |
| M134 | Alexandrova et al. (2009) — Universality of Solar-Wind Turbulent Spectrum from MHD to Electron Scales | Alexandrova-2009-PRL | Ion-scale spectral index |

The following **9 pre-2006 foundations** are also not in the supplied bibliography and are included in the additions file:

| Reference | Paper | Proposed key | Role |
|---|---|---|---|
| M01 | Palmer (1982) — Transport coefficients of low-energy cosmic rays in interplanetary space | Palmer-1982-RG | Palmer consensus for λ∥ at 1 AU |
| M09 | Bieber et al. (1994) — Proton and Electron Mean Free Paths: The Palmer Consensus Revisited | Bieber-1994-AJ-420 | Electron/proton mean-free-path discrepancy; slab/2D composite |
| M19 | Dröge (2000) — The Rigidity Dependence of Solar Particle Scattering Mean Free Paths | Droge-2000-AJ | Rigidity dependence of event-fitted λ |
| M54 | Shalchi et al. (2005) — Spurious contribution to cosmic ray scattering calculations | Shalchi-2005-MNRAS | Magnetostatic slab D_μμ and λ∥ asymptote |
| M64 | Teufel et al. (2003) — Analytic calculation of the parallel mean free path of heliospheric cosmic rays. II. Dynamical magnetic slab turbulence and random sweeping slab turbulence with finite wave power at small wavenumbers | Teufel-2003-AA | Slab QLT λ∥ asymptotics used by the NWU constructions |
| M66 | Teufel et al. (2002) — Analytic calculation of the parallel mean free path of heliospheric cosmic rays. I. Dynamical magnetic slab turbulence and random sweeping slab turbulence | Teufel-2002-AA | Slab QLT λ∥ asymptotics (g = 0 below k_min) |
| M72 | Zank et al. (1999) — Solar Wind Turbulence, Diffusion Coefficients, and Cosmic Ray Modulation | Zank-1999-ICRC | Printed form of the Zank 1998 fit |
| M122 | Leamon et al. (2000) — MHD-driven Kinetic Dissipation in the Solar Wind and Corona | Leamon-2000-AJ | Dissipation-range onset fits |
| M126 | Zank et al. (1996) — Evolution of turbulent magnetic fluctuation power with heliospheric distance | Zank-1996-JGR | Radial exponents of δB² |

### 18.2 Cited papers found in the supplied bibliography

| Reference | Matching key(s) | Notes |
|---|---|---|
| M02 | Tautz-2013-JGRSP | — |
| M03 | Chhiber-2017-AJSS | — |
| M04 | Reames-2013-SSR | Reames-2013-SSR: library entry has no DOI |
| M05 | Lavasa-2026-AA | — |
| M06 | Minoshima-2026-EPS | matched by title; library record supplies EPS 78 and DOI |
| M07 | Subashchandar-2025-AJL | — |
| M10 | Engelbrecht-2022-SSRv | — |
| M11 | Chen-2024-AJ-965 | Chen-2024-AJ-965: library entry has no DOI |
| M12 | Lang-2024-arXiv | only the arXiv version is in the library; the published ApJ 971, 105 record is missing |
| M21 | Dalla-2020-AA | library record lacks volume and DOI |
| M22 | Houeibib-2025-AA | matched by first author, year, journal and article (A211); research notes did not record the title |
| M23 | Droge-2009-AJ | Droge-2009-AJ: library entry has no DOI |
| M26 | Zhong-2024-AJ-974 | — |
| M27 | Cao-2025-AA | — |
| M28 | Verkhoglyadova-2009-AJ | Verkhoglyadova-2009-AJ: library entry has no DOI |
| M29 | Schwadron-2010-SW | Schwadron-2010-SW: library entry has no DOI |
| M30 | Kozarev2013 | — |
| M31 | Borovikov-2019-arXiv | — |
| M32 | Liu-2025-AJ | — |
| M33 | Zhao-2024-SW | — |
| M34 | Zhang-2023-AJSS | — |
| M35 | Afanasiev-2025-JSWSC, Afanasiev-2024-JSWSC | Afanasiev-2024-JSWSC duplicates the record with doi "n/a" and no volume |
| M36 | Aran-2005-AG | — |
| M37 | Marsh-2015-SW | Marsh-2015-SW: library entry has no DOI |
| M39 | He-2011-AJ | He-2011-AJ: library entry has no DOI |
| M40 | Wang-2004-arXiv | library record has year 2004 and journal "arXiv:1311.7469v4" (published ApJ 799, 111, 2015; no DOI in the record) |
| M41 | Kubo-2015-EPS, Kubo-2015-PS | Kubo-2015-PS has the journal split into the year field |
| M42 | Wijsen-2019-AA-A28 | Wijsen-2019-AA-A28: library entry has no DOI |
| M43 | Laitinen-2016-AA | Laitinen-2016-AA: library entry has no DOI |
| M44 | Laitinen-2018-JSWSC, Laitinen-2018-arXiv | duplicate arXiv record |
| M46 | Strauss-2015-AJ | Strauss-2015-AJ: library entry has no DOI |
| M48 | Young-2021-AJ | Young-2021-AJ: library entry has no DOI |
| M49 | He-2015-AJ-218 | He-2015-AJ-218: library entry has no DOI |
| M50 | He-2015-AJ-814 | — |
| M51 | He-2019-AJL, He-2019-arXiv | duplicate arXiv record; neither has a DOI |
| M52 | Tenishev-2005-AIAA-4832 | Tenishev-2005-AIAA-4832: library entry has no DOI |
| M53 | Tenishev-2022-arXiv | — |
| M55 | Wang-2014-AJ | Wang-2014-AJ: library entry has no DOI |
| M56 | Qin-2015-AJ | Qin-2015-AJ: library entry has no DOI |
| M57 | Droge-2010-AJ | Droge-2010-AJ: library entry has no DOI |
| M59 | Kelly-2012-arXiv | only the arXiv version is in the library; published ApJ 750, 47 |
| M60 | Ding-2022-AA, Ding-2022-arXiv | duplicate arXiv record |
| M61 | Afanasiev-2015-AA, Afanasiev-2016-arXiv, Afanasiev-misc | duplicates: arXiv copy (2016) and an undated misc record; Afanasiev-2018-AA-618 has the same title with A&A 618, A114 (Section 18.3) |
| M62 | Strauss-2017-AJ-837 | — |
| M65 | Engelbrecht-2013-AJ-779 | — |
| M68 | Engelbrecht-2013-AJ-772 | — |
| M71 | Zank-1998-JGR | — |
| M73 | Perri-2020-JSWSC | matched by title; library record supplies JSWSC 10, 55 and DOI |
| M75 | Hussein-2015-JGRSP | — |
| M76 | Engelbrecht-2018-AJ | — |
| M77 | Kozarev-2016-AJ | — |
| M78 | Battarbee-2010-AIP | Battarbee-2013-AJ duplicates this title with ApJ 658, 622-639 and an arXiv DOI (Section 18.3) |
| M79 | Ng-2008-AJL, Ng-2008-NASA | Ng-2008-NASA duplicates the paper with an invalid DOI 10.48550/arxiv.20080038680 |
| M80 | Ng-2012-AIP | — |
| M81 | Hu-2017-JGR | Hu-2017-JGR: library entry has no DOI |
| M82 | Ding-2020-RAA, Ding-2020-arXiv | duplicate arXiv record |
| M83 | Luhmann-2010-ASR | — |
| M84 | Potgieter-2013-LRSP | — |
| M87 | Vos-2015-AJ | — |
| M93 | Corti-2019-AJ | — |
| M95 | Boschini-2018-ASR | matched by DOI |
| M98 | Boschini-2017-AJ | — |
| M99 | Bobik-2012-AJ | Bobik-2012-AJ: library entry has no DOI |
| M100 | Strauss-2011-AJ | — |
| M104 | Fiandrini-2021-PRD | — |
| M108 | Qin-2017-AJ | Qin-2017-AJ: library entry has no DOI |
| M109 | Shen-2018-AJ | Shen-2018-AJ: library entry has no DOI |
| M116 | Zhao-2018-AJ | Zhao-2018-AJ: library entry has no DOI |
| M117 | Oughton-2011-JGR | — |
| M121 | Bieber-1996-JGR | Bieber-1994-AJ carries this title with ApJ 420 metadata (Section 18.3) |
| M127 | Smith-2001-JGR | — |
| M128 | Breech-2008-JGR | Breech-2008-JGR: library entry has no DOI |
| M133 | Chen-2019-AJSS | library record has year 2019 and no volume, pages or DOI (published ApJS 246, 53, 2020) |
| M135 | Raath-2016-ASR | matched by authors, year and arXiv record; title in the library differs from the paraphrase in the research notes; library entry has no DOI |

### 18.3 Records in the supplied bibliography that need attention

These records were encountered while matching. They are reported for correction in the bibliography; this document does not edit the supplied file.

| Record | Issue |
|---|---|
| Bieber-1994-AJ | Title, authors and JGR identity of Bieber, Wanner & Matthaeus (1996) combined with ApJ 420, 294–306 and DOI 10.1086/173566. The paper Bieber et al. (1994), "Proton and Electron Mean Free Paths: The Palmer Consensus Revisited", ApJ 420, 294, DOI 10.1086/173559 [M09], is therefore not in the bibliography under its own title; it is supplied in the additions file under the key Bieber-1994-AJ-420. |
| Wang-2004-arXiv | Year 2004 for arXiv:1311.7469 (Wang & Qin, published ApJ 799, 111, 2015). |
| Battarbee-2013-AJ | Title of Battarbee et al. (2010, AIP Conf. Proc.) with ApJ 658, 622–639 (the volume and pages of Vainio & Laitinen 2007) and DOI 10.48550/arxiv.1303.4334. |
| Afanasiev-2018-AA-618 | Same title as Afanasiev et al. (2015, A&A 584, A81) but A&A 618, A114 (2018); one of the two identities is wrong. |
| Afanasiev-misc | Undated duplicate of Afanasiev et al. (2015). |
| Afanasiev-2024-JSWSC | Duplicate of Afanasiev-2025-JSWSC with doi "n/a". |
| Kubo-2015-PS | Journal "Earth, Planets and Space" split into the year field ("Planets and Space"). |
| Ng-2008-AJL | DOI 10.1086/592998 belongs to an unrelated paper; the DOI of Ng & Reames (2008, ApJL 686, L123–L126) is 10.1086/592996. |
| Ng-2008-NASA | Duplicate of Ng-2008-AJL with the invalid DOI 10.48550/arxiv.20080038680. |
| Chen-2019-AJSS | Year 2019 and no volume, pages or DOI for Chen et al. (2020, ApJS 246, 53). |
| Lang-2024-arXiv | arXiv record only; the published ApJ 971, 105 is missing. |
| Kelly-2012-arXiv | arXiv record only; the published ApJ 750, 47 (DOI 10.1088/0004-637X/750/1/47) is missing. |
| Dalla-2020-AA | No volume or DOI (A&A 639, A105; DOI 10.1051/0004-6361/201937338). |
| Raath-2016-ASR | No DOI (10.1016/j.asr.2016.01.017). |
| Kubo-2015-EPS | No DOI (10.1186/s40623-015-0260-9). |
| Minoshima-2026-EPS | No article number (64). |
| He-2019-arXiv / He-2019-AJL; Laitinen-2018-arXiv / Laitinen-2018-JSWSC; Ding-2022-arXiv / Ding-2022-AA; Ding-2020-arXiv / Ding-2020-RAA; Afanasiev-2016-arXiv / Afanasiev-2015-AA | arXiv duplicates of published papers. |
| Several matched records (for example Chen-2024-AJ-965, Droge-2010-AJ, He-2011-AJ, Wijsen-2019-AA-A28, Laitinen-2016-AA, Strauss-2015-AJ, Bobik-2012-AJ, Breech-2008-JGR, Zhao-2018-AJ) | No DOI field; matched by title. |

## 19. Companion data

The bundle MEAN_FREE_PATH_MODEL_DATA.zip has this layout:

| Path | Content |
|---|---|
| `README.md` | layout, column definitions, provenance and verification labels, parsing rules, regeneration instructions |
| `observations/event_fitted_mfp.csv` | 114 rows: Pacheco et al. 2019 Tables 2–4, Agueda et al. 2014 Table 3, Lang et al. 2024 Table 3, Lavasa et al. 2026 Tables 2–3; quantity tag (lambda_r or lambda_par), energy, rigidity, distance, verification, flags |
| `parameters/sep_code_presets.json` | 22 SEP presets of Section 5.2 with form, quantity tag, variable, parameters as printed, source location, verification and runtime state (Section 21.1) |
| `parameters/gcr_nwu_family_parameter_sets.csv` | 461 rows: Potgieter et al. 2014 Table 1; Vos & Potgieter 2015 Tables 1–2; Potgieter et al. 2015 Table 1; Aslam et al. 2019 Table 1; Aslam et al. 2021, 2023 text values; Corti et al. 2019 Tables 2–3 (partial); Luo et al. 2019 Table 1 |
| `parameters/gcr_helmod.json` | Boschini et al. 2018 (ASR) Tables 1–2 with the D-21 numerical-use convention; Boschini et al. 2019 Table 2; Bobik et al. 2012 Table 1; scalar values of Della Torre et al. 2016 and Boschini et al. 2017 (ApJ) |
| `parameters/gcr_other_models.json` | Strauss 2011, Effenberger 2012, Wang 2019, Tomassetti 2017 (Table II), Perugia 2021/2025 (with the Tomassetti et al. 2023 printing, D-23), Jiang 2023, Duan 2025 (D-33), Qin & Shen 2017/Shen & Qin 2018, Engelbrecht & Burger 2013a |
| `parameters/qlt_parameter_sets.json` | Section 7.4 sets (TS-BIEBER94, LANG24-T5, LANG24-T6-DT/RS/P, EB13b), the TS2002 and TS2003 printed approximations (Section 7.3), Shalchi et al. 2005 Eq. 27 and 2006 Table 2, Hussein et al. 2015 Tables 3 and 5 (Table 2 stored but not usable as a benchmark), Strauss et al. 2017 Table 1, and the Chhiber et al. 2017 and Zhao et al. 2018 λ∥ evaluations |
| `parameters/model_registry.json` | every model ID with section, equation, source, implementation status and runtime state |
| `turbulence/published_turbulence_values.csv` | 139 rows of Section 11 values with definition, condition, verification and flags |
| `benchmark_points.json` | 34 fixture groups of Section 13 |
| `source_key_map.json` | source key → [Mxx], citation key, DOI, arXiv identifier; equation labels → numbers |
| `reference_verification.py` | compares every numerical fixture value with an expectation computed without the generator (derived values recomputed, with closed forms where they exist; row coordinates and printed constants compared with the values stated here), with the tolerances of Section 13; checks the provenance fields of every data file, the runtime states and SHA256SUMS |
| `independent_physics_checks.py` | mathematical checks of the identities used by the specification against numerical quadrature, numerical differentiation and a seeded stochastic simulation (Section 21.9); it does not check transcriptions |
| `SHA256SUMS` | digests of all other files |

All transcribed numbers are stored as strings exactly as printed, including "≤", "≥", "~", ranges and uncertainties, so that no rounding or reinterpretation enters the data. Every row or object names its source key (mapped to the [Mxx] identifiers of Section 17 by `source_key_map.json`), its location in the source, and its verification label. The label "SIGN NOT LEGIBLE" marks a value whose sign could not be read in any accessible rendering.

**Implementation requirement.** A loader converts a printed string to a number only when no interpretation is needed: a plain decimal or scientific-notation number, or an exact rational literal of the form integer/integer (for example "1/3", "2/3", "5/3") in a field declared dimensionless (an exponent or index). Any field containing a qualifier, range, uncertainty or unit-only normalization is surfaced to the caller rather than converted (Section 21.3).

## 20. Verification and revision record

### 20.1 Versions

- **Version 1.0, 8 October 2026:** first issue.
- **Version 1.1:** an external review of version 1.0, consisting of a findings document and a revised text. Its findings were checked against the sources and the mathematics before any was adopted.
- **Version 1.2, 9 October 2026:** this version. It adopts the findings of version 1.1 that were confirmed, corrects those that were not, and adds the results of a second reading of the sources (Section 20.4).

**Adopted from version 1.1 (each confirmed independently):**

| Item | Where |
|---|---|
| The Duan et al. (2025) printed k₂ has slopes a below and b above R_k only for b > a (D-33); slope identity Eq. (47); fixture F-GCR-09 for b < a | Sections 9.5, 14.1, 21.8; F-GCR-08 relabelled |
| λ∥ = 3κ∥/v for SEP-CHEN24 is the definition Eq. (1), not an assumption; only the species (U-2) is open, and U-1 becomes the caller's output choice | Sections 3, 4.5, 5.3, 14.2 |
| pc/1 GeV and 𝓡/1 GV coincide for every singly charged species, not only protons and electrons | Section 2; D-17 |
| HelMod K₀: the numerical-use convention is named explicitly (merged into D-21 rather than a second register entry) | Sections 9.4, 14.1, 21.8 |
| The B₀ of the Parker field is the radial-component amplitude; a field magnitude at Earth is not interchangeable with it (D-34, U-13), confirmed against Raath et al. (2016), whose wording calls B₀ the magnitude | Sections 11.8, 14, 21.2 |
| Runtime states per configured run, with the added state REQUIRES_USER_DECISION | Sections 3, 12, 21.1; `model_registry.json`, `sep_code_presets.json` |
| Exact rational literals ("1/3") are parsed only in fields declared dimensionless | Sections 19, 21.3 |
| The version 1.0 verifier compared values with the tolerance 10⁻¹² × max(1, \|value\|), which is absolute, not relative, for values below 1; version 1.2 uses relative tolerances for all nonzero reference values, and absolute bounds only for zero references and for values that are themselves 40-digit differences | Sections 13, 21.9; `reference_verification.py` |
| Implementation contracts: typed units, parsing, domains, stable evaluation (He & Wan series, expm1/log1p in the Zank fit), the Itô form of the pitch-angle operator and its first-moment test, conservation and boundary tests for the printed (1−μ)² variant, turbulence and shock interfaces, validation layers L0–L3, blocker policy | Section 21 |
| The document is a self-contained catalogue but delegates equations to the companion specifications; a parameter with the same symbol in two publications is not automatically the same coefficient | Sections 1.1, 9.3 |
| Exact slab-QLT integrals as validation references: the Gauss hypergeometric form for the bendover spectrum and the piecewise closed form for the TS2003 spectrum | Section 21.4, Eqs. (44), (45); F-NUM-01 |
| A second script of independent mathematical checks | `independent_physics_checks.py` |
| Roadmap MF00–MF16 as the final section, replacing MF00–MF10 | Section 22 |

**Not adopted, or corrected before adoption:**

- Version 1.1 lists MEAN_FREE_PATH_MODEL_additions.bib as a blocker because the file was not among the attachments it reviewed; the file was issued with version 1.0 as a separate attachment and is issued again with version 1.2, so the blocker row is omitted.
- Version 1.1 states that the companion specifications were not attached. That describes the material available to the reviewer, not this specification; it is replaced by the requirement to pin the companion revisions and check every cited section and equation number against them (Sections 1.1, 21.10).
- A separate register entry for HelMod would duplicate D-21. The register identifiers of version 1.1 map as follows: its D-33 (Duan) is kept; its D-34 (HelMod unit) is merged into D-21; its D-35 (verifier tolerance) and D-36 (rational literals) are resolved in the verifier and in Sections 19 and 21.3 and have no register entry; its D-37 (Parker B₀) is D-34 here.
- Version 1.1 retained the version 1.0 statement that fixture regeneration agreed "to a relative 10⁻¹²" while the verifier's tolerance was absolute for small values; the statement now matches the verifier.
- Repository paths named in the version 1.1 roadmap were not verified in any repository and are not named; MF00 locates the host layout by inspection.
- The He & Wan series of version 1.1 (through x⁶) is extended by one term, given with its truncation bound and a switching rule that is tested in double precision.

**Other corrections in version 1.2:**

- From the second reading of the sources (Section 20.4): the SEP-EPREM10 radial factor (present in Schwadron et al. 2010, Eq. 2, with an illegible exponent sign); SEP-SOFIE24 is a parallel mean free path (stated once, Sect. 4.5), so it leaves U-11; the SEP-SPARX15 and SEP-MARSH13 forms and the SEP-SPARX15 wording; the SEP-KUBO15 values (scanned, not fitted, with three distinct reference conditions); the TS2002 printed approximations, which are λ/λ₀ with ≪ ranges and a first-power low-rigidity branch; the Dröge et al. (2016) α values and the Strauss et al. (2017) quotation in Section 10; the polar functions accompanying the 0.01 and 0.02 perpendicular ratios of Vos & Potgieter (2015, 2016), Corti et al. (2019) and Duan et al. (2025); the Tomassetti et al. (2023) printed form (D-23); D-15 (the Breech et al. sign is not legible in any accessible rendering); the identification of "Boschini et al. 2018a" with Boschini et al. (2018, ASR); four turbulence-table locations and wordings (Smith et al. 2001 twice, Breech et al. 2008 Fig. 5, Usmanov et al. 2018).
- SEP-SOLPENCO05 was "implementable" in the version 1.0 registry while Sections 12 and 14 listed it under U-11; the registry now agrees with the text.
- The F-GCR-02 fixture string quoting Potgieter et al. (2014) used the v1–v2 wording; it now quotes the final version, as Section 9.3 does. No numerical fixture value of version 1.0 changed.
- Raath (2015) is an MSc dissertation (Magister Scientiae); the metadata note is removed.
- Added in version 1.2 (not in version 1.1): the explicit no-flux condition of the printed (1−μ)² variant, D_μμ(−1) = 4D₀(1 + H) ≠ 0; the second-moment decay test; the statement that the Itô drift is singular or discontinuous at μ = 0 for the q-, M-FLAMPA and ε-shapes for every H; the s-dependence of the Eq. (21) error; checks of the full prefactors of Eqs. (44) and (45) against Eq. (11).
- Wording and labels: the Raath (2015) row (no longer "reconstructed"; P′₀, B₀ and the βc notation added); the SEP-LAITINEN18 variable (no rigidity scaling printed); the Marsh et al. wording also in Section 6.6; the SEP-DROGE16P values described as "assumed" and "determined"; the Alexandrova et al. (2009) value written as printed (~k^−2.8); the Oughton et al. (2011) location (Fig. 1 caption); the Tomassetti (2017) field reference in Section 9.6 (final version); the inputs of the F-SEP-02 Parker example (r☉ = 0.005 AU, θ = 90°) and R_k = 4.3 GV for F-GCR-08/09 in Section 13; the PARASOL source equation written "their Eq. 43" to avoid confusion with Eq. (43); the scope of `model_registry.json` (model IDs; the OBS-* and TURB-* groups are the data files); the verification counts, which now give the total number of agreeing readings over versions 1.0 and 1.2 (Section 1.2).
- Table cells containing |…| are escaped so that tables render correctly.

### 20.2 How the content was obtained

The literature was reviewed in five topic areas (observational constraints, SEP code prescriptions, closed-form theory, GCR modulation prescriptions, turbulence inputs). All sources were read through a web-fetching tool that returns model-generated renderings of PDF or HTML pages; this is why every item carries a fetch count, and why critical numbers were read from two renderings where possible. Publisher sites that refused automated access (IOPscience, Wiley/AGU, ADS) were replaced by arXiv, ar5iv, institutional repositories and publisher table pages. The web-search budget of the version 1.0 review was exhausted; the items in Section 14.3 could not be pursued further.

### 20.3 Independent checks of version 1.0

Five independent checks were run on the drafted version 1.0. None of the checkers had seen the drafting notes or the fixture generator.

| Check | Scope | Result | Action |
|---|---|---|---|
| Observations and turbulence tables | 277 items: every row of `event_fitted_mfp.csv`, Lang et al. 2024 Tables 5–6 and Eqs. 2, 3, 9–13, and `published_turbulence_values.csv` | 249 match, 2 mismatch, 26 not checked | Corrected the inverted Subashchandar et al. (2025) slab/2D ratio (and the conflict it had created in Section 11.2) and the Lang Eq. (3) exponent symbol; added three Kallenrode (1993) values and four γ values that had been omitted; the Lavasa Table 3 rows are now confirmed in two renderings |
| GCR tables and equations | 711 items: every cell of the GCR CSV and JSON files and the Section 9 equations | 678 match, 26 mismatch (5 distinct issues), 7 not checked; every numeric table cell matched | Corrected the Corti λ⊥ unit (AU), the Potgieter et al. 2014 sentence (final arXiv version), the Tomassetti B₀ wording (final version), the Vos & Potgieter equation numbers, and the Potgieter et al. 2015 row-label wording |
| SEP and theory equations | 164 items: forms, coefficients, parameter values and quotations of Sections 4–8 and 11 | 145 match, 14 mismatch, 5 not checked; no coefficient or parameter value was wrong | Corrected non-verbatim quotations (Kozarev, Liu, TS2003, Laitinen, Wijsen), notation (Lang Eq. 3, Minoshima Ω_n, Afanasiev R²), the attribution of the h = 0.01 quotation, the Strauss et al. 2020 equation number, and three factual statements (Borovikov convention, Zhang V and B₀, Wijsen B₀) |
| Bibliographic metadata | all 135 references against Crossref, OpenAlex, arXiv, INSPIRE and publisher pages | 82 confirmed, 45 corrected (mostly missing volume, issue, pages or DOI; ten preprints found published), 8 partial | Applied all corrections, including one wrong DOI (Ng & Reames 2008) and two titles; the remaining unconfirmed fields are listed as metadata notes in Section 17 |
| Mathematics and consistency | independent recomputation of 52 derived claims and of every number in Section 13; equation numbering, cross-references, citations, units | all derived identities and printed fixture values reproduced; 5 errors (one numerical, four reference errors) and 11 warnings | Corrected the Zank-bracket accuracy statement (maximum deviation 2.16% at R_L/l = 5.787, previously stated as 1.8%), all broken references, and the ambiguous citations; added fixture F-QLT-07 quantifying the accuracy of Eq. (21) |

### 20.4 Second reading of the sources for version 1.2

The items that version 1.0 listed as resting on a single reading, and the remaining parameter sets and equations, were read again from the sources by four independent checkers that had not seen the drafting notes.

| Check | Scope | Result | Action |
|---|---|---|---|
| SEP presets and observational rows | 109 items: presets EPREM10, SOFIE24, SOLPENCO05, SPARX15, MARSH13, KUBO15, STRAUSS17G, STRAUSS15, DROGE16P, LAITINEN18 (JSON and Section 5.2 rows); the He et al. λ_r/λ∥ pairs; the Kallenrode values of Pacheco et al. (2019, Table 4) and the γ values of Agueda et al. (2014, Table 3) | 98 match, 10 mismatch (4 presets, each in the JSON and the table), 1 not checked (an interpretive label, reworded) | Corrections listed in Section 20.1; SOLPENCO05, MARSH13, STRAUSS17G, STRAUSS15 and DROGE16P confirmed; all observational rows confirmed |
| Turbulence values | 26 rows of 11 sources not re-read for version 1.0 | 18 match, 4 mismatch (locations and wording; no value wrong), 4 not checked; two of these (Engelbrecht & Burger 2013a Table 1) were then read in the copy of the published paper found by the GCR check and match | Locations and wording corrected; the Oughton et al. (2011) temperature 1.6 × 10⁶ K confirmed in a copy of the published paper (D-16); Leamon et al. (2000) Table 1 read a second time; the Breech et al. (2008) sign is not legible (D-15) |
| GCR items | 48 items: Raath (2015) Eq. 3.31 and units, Engelbrecht & Burger (2013a) Eqs. 15–16, Tomassetti et al. (2023), perpendicular ratios, HelMod tilt definitions, Effenberger et al. (2012), Duan et al. (2025) | 40 match, 7 mismatch, 1 not checked (whether Duan's fitted b exceeds a; values only in a figure) | Corrections listed in Section 20.1; Effenberger et al. (2012) fully confirmed, B_e being "the magnetic field strength at r_e = 1 AU" |
| QLT parameter sets and equations | 88 items: TS2002 approximations, Shalchi et al. (2005, 2006), Hussein et al. (2015), Strauss et al. (2017), Chhiber et al. (2017), Zhao et al. (2018), Eq. (40) and its sources, the Raath et al. (2016) field, Eq. (38) and the Section 10 ratios | 81 match, 5 mismatch, 2 not checked (Engelbrecht & Burger 2013a, then read by the GCR check and matching) | TS2002 entries, Dröge α values and the Strauss quotation corrected; the B₀ wording of Raath et al. recorded as D-34; the derivation of Section 10 recomputed and confirmed |

### 20.5 What rests on a single reading or was not verified

- Zank et al. (2017) 80:20 partition: rests on the version 1.0 reading; no copy was accessible for the second reading.
- The sign of σ_D in Breech et al. (2008) and of the 3/2 exponent in Schwadron et al. (2010, Eq. 2): minus signs are lost in every accessible rendering.
- Only one accessible copy (in some cases read in both versions): Leamon et al. (2000), Oughton et al. (2011), Smith et al. (2001), Breech et al. (2008), Usmanov et al. (2018), Engelbrecht & Strauss (2018; Eq. 17 layout not read directly), Engelbrecht & Burger (2013a, b), He et al. (2011), Raath (2015) and Vos & Potgieter (2015); Zhao et al. (2024) and Strauss & Fichtner (2015) were read as preprints only; Bieber et al. (1996) is abstract only.
- AGU citation numbers of Schwadron 2010, Breech 2008, Oughton 2011, MacBride 2010 and Weygand 2011, the page range of Hu et al. 2017, and the AIP volume of Ng, Reames & Tylka 2012 are unconfirmed.

### 20.6 Companion bundle check

The issued bundle was extracted into an empty directory and both scripts were run in isolated mode. `reference_verification.py`: 5 542 checks passed and none failed. They compare all 548 numerical values of the 34 fixture groups with expectations computed without the generator, at the tolerances of Section 13 (445 derived values recomputed, 84 row coordinates and 19 printed or illustrative constants), and check the digests and the structure, runtime states and provenance fields of every data file. `independent_physics_checks.py`: 38 checks passed and none failed. These results test mathematics, structure and integrity; they are not a verification of the transcriptions (Sections 20.3–20.4) and not an end-to-end transport run. The SHA-256 digests of the issued bundle are:

~~~text
d37a3763400f64d216e336c049cf83732ad641360ec02d4e1acf8b49694c1706  README.md
2a461634a2c6f8a55885efe4875e89b97885e8c95fd18c56d6ce37ffa329f528  benchmark_points.json
582e70d4d718b7bb9f53340da41fef3623826a7b5cacdaf619cf9f7a2cb64ecc  independent_physics_checks.py
63e94ed726cb2adfd438a73aac7904928e22a4e4ddf2c0847ee4d20dd46228c9  observations/event_fitted_mfp.csv
77bc3c86eba90d14c10ca33928f4e17ffdd486c909bc76cf47a1a6a381039c91  parameters/gcr_helmod.json
6d7c66442bf96a2a9d4cbec49484e5970b550025f327db2c1113d903289af90f  parameters/gcr_nwu_family_parameter_sets.csv
f399ec3f815b1cd939fa98a9c7bce08f10f7be5c4d6ebc2af65abecda346d0b3  parameters/gcr_other_models.json
e308a3cd8908c5aacb7bd2ed0b7c33a8d36d0e71a1cef4a93bcf68da7684c4aa  parameters/model_registry.json
7b6a9424303de51861df029855346bfc96957ccbd9d6e8174480e0e3afd6c360  parameters/qlt_parameter_sets.json
20e186840d2d729de0a800b85a50ba5f833a23edaf1e656151fb73b5b007f818  parameters/sep_code_presets.json
1836610ec61080f066fdb5d35de7db4af00e321f4553b726c138acd7d2d3eaaf  reference_verification.py
8415bb8831b3e041158dfb6a4dededa40db7f3d09bd4f2d00c2bd40996dd6f42  source_key_map.json
5e3b7daa3db9f38cd274a3c5543d6829f604a758766d47d304ca93784f410130  turbulence/published_turbulence_values.csv
~~~

## 21. Implementation contracts, numerical safeguards and validation design

This section states engineering and mathematical requirements that follow from Sections 1–19; it adds no scattering physics. Equations (1)–(42), the transcriptions and the published parameter sets keep their documented provenance. The equations introduced here are exact rewritings or consequences of those equations and are labelled derived. A choice that would change the physical closure is a configuration decision, never an implementation convenience. Delegated equations (PARALLEL rev. 1.4, PERPENDICULAR) are implemented only against pinned revisions of those documents (Section 21.10).

### 21.1 Runtime states

The inventory of Section 3 lists what the literature provides; it is not a promise that every entry can run without further information. A configured run is in exactly one of these states:

| Runtime state | Definition | Examples |
|---|---|---|
| `READY_EXPLICIT_INPUTS` | The law is specified, and every input, convention and unit it needs is supplied by the configuration or the source | tagged SEP power-law presets, q-form and ε-form shapes, isotropic scattering, Bohm, Eq. (21) with named turbulence inputs, the Zank fit with a named variance convention, a GCR form with a named published parameter set |
| `READY_PUBLISHED_VARIANT` | A formula is reproducible only under an explicitly named printed variant or operator convention | DT electron printings (D-1), the EPREM half-D operator (D-2), the printed (1−μ)² factor (D-3), the printed ε-form (D-4), the printed Duan law (D-33) |
| `READY_PUBLISHED_KAPPA_ONLY` | The published quantity is a diffusion coefficient; a mean free path needs the particle speed and hence a species choice | SEP-CHEN24: κ∥ is the published output, λ∥ = 3κ∥/v is derived after U-2 |
| `REQUIRES_EXTERNAL_INPUT` | The formula is known, but required time-dependent, turbulence or shock inputs are absent | Eq. (21) without slab parameters (U-9); GCR models whose time-dependent parameters are printed only in figures (Wang 2019, Perugia 2021/2025, Duan 2025) |
| `REQUIRES_USER_DECISION` | The formula and inputs exist, but a decision listed in Section 14.2 is open | λ_unspecified presets (U-11), Zank without a variance convention (U-12), NWU variants with the D-9 printed exponents (U-10), a field normalization of Eq. (42) (U-13), the K₀ unit of Jiang 2023 |
| `REQUIRES_SOURCE_OR_CODE_AUDIT` | The equation, or what a host-code mode implements, is not established | `Tenishev2005AIAA` (U-4), the `Chen2024AA` mapping (U-3), SHOCK-PARASOL (U-6), SEP-MINOSHIMA26 (D-29), GCR-BOBIK12 units, SEP-EPREM10 away from 1 AU |
| `REFERENCE_DATA_ONLY` | An observational fit or a published comparison value, not a predictive closure | event-fitted λ (Section 4), the Palmer band, GCR-LUO19 (form not legible), published turbulence values |

The state is evaluated per configured run, not once per model ID: a publication's tabulated value, a model family and an executable configuration are different things. `model_registry.json` and `sep_code_presets.json` record the state that applies when the inputs named in the source are supplied, with a note where it depends on the configuration. An unknown or incomplete configuration fails with a structured error that names the missing parameter, source item, or D-/U-identifier. A normalization stored as a UNIT (Section 9.6) is never used as a VALUE. A model that needs a delegated companion equation (for example the PARALLEL QLT adapter or NLGCE-F) whose pinned text is not available is not assigned one of these states: the implementation reports it as `BLOCKED_DEPENDENCY` (Section 21.10).

### 21.2 Typed quantities and units

The evaluation path keeps these quantities distinct: species (rest mass, charge number Z with its sign), total kinetic energy, kinetic energy per nucleon, momentum pc, rigidity, position and time, `LambdaKind`, κ∥, κ⊥ and D_μμ. Internally lengths are in metres, speeds in m s⁻¹, fields in tesla, wavenumbers in m⁻¹, diffusion coefficients in m² s⁻¹, D_μμ in s⁻¹ and rigidities in volts. Input and output may use AU, R☉, GV, MeV, nT or cm² s⁻¹, through explicit and tested conversions.

- For |Z| > 0, pc and rigidity obey Eq. (4). At fixed pc every singly charged species has the same rigidity; at fixed kinetic energy the rigidity depends on the mass.
- A speed of zero is an invalid argument of λ = 3κ/v and of the amplitude ↔ λ conversions. Speed, energy and rigidity must be positive and finite; no small positive value is substituted unless a numerical-policy setting says so.
- The isotropic-scattering length λ_iso of a full-orbit code (Section 6.6) is not converted into λ∥ by changing its tag.
- Eq. (2) is applied only on request and only with the local field angle ψ. With a nonzero perpendicular coefficient the radial tensor element is Eq. (3), not Eq. (2). A radial-to-parallel conversion is rejected when cos²ψ = 0 or when ψ is missing or not finite.
- Field magnitude |B|, the radial component |B_r| and the component normal to a surface are different quantities. For the unmodified Parker field of Eq. (42), |B(r₀)| = |B_r(r₀)|/|cos ψ(r₀)| (derived). A published field value at 1 AU is used as the B₀ of Eq. (42) or as the magnitude only as the configuration states (D-34, U-13).
- The relativistic gyrofrequency qB/(γm) and the nonrelativistic cyclotron frequency qB/m differ for γ ≠ 1, and angular (rad s⁻¹) and cyclic (Hz) frequencies differ by 2π. The undefined Ω_n of D-29 is not mapped to either automatically.
- QLT expressions printed in Gaussian units (for example Ω = qB/(mc), 4πρ) are converted by unit system, not copied into an SI implementation.

A returned value carries at least: model ID, formula variant, quantity kind, value in SI, species, local time or snapshot ID, source key and location, original unit, whether it is derived, the active D-/U-decisions, and domain warnings. A host may hold this in its existing types; published plots must be reproducible from it.

### 21.3 Parsing the companion data

The bundle is a verified transcription archive, not a runtime configuration. A loader keeps the raw strings (`value_as_printed` and the JSON strings) for audit, and may build a separately versioned numeric view only from strings that parse without interpretation.

| Printed form | Permitted handling | Not permitted |
|---|---|---|
| Plain decimal or scientific notation (`0.3`, `2.5e-3`) | parse, with the dimension declared by the field | inferring a unit from nearby text |
| Exact rational literal integer/integer (`1/3`, `2/3`, `5/3`) | parse as an exact rational, only in a field declared dimensionless (exponents, spectral indices); keep the source string | reading it as a date, or rejecting a published exponent while reporting the preset as runnable |
| Number with an explicit unit (`0.3 AU`, `5 nT`) | parse through a field-specific unit schema | dropping the unit, or confusing cm with m |
| Range, ± uncertainty, approximate value, inequality, UNIT normalization | keep as a tagged structure until an explicit selection is made | replacing a range by its midpoint, dropping an uncertainty, treating a UNIT as a VALUE |
| Undefined symbol, incomplete equation, conflicting printed variants | reject the configuration with the D-/U-identifier | substituting a nearby published formula |

The loader records the bundle version and its SHA-256 digest and the files a run used. A change of reference data never changes default physics silently. Schema validation requires, for every transcribed row or object, a nonempty source key, location and verification label, using each file's own field names: in the CSV files `source`, `location` and `verification`; in `sep_code_presets.json` `source`, `location` and `verification`; in `gcr_other_models.json` and `qlt_parameter_sets.json` `source` (which may carry the location, for example "TS2003 Eq. (51)") with `location` where present and `status`; in `gcr_helmod.json` the object key as source with the `equation` or `table` text as location and `status`; and in `model_registry.json` `source`, `equation` and `status`. Model and fixture IDs are unique. The observational CSV is benchmark data; it is not a source of code defaults or tuning targets without an explicit calibration workflow.

### 21.4 Domains and stable evaluation

Every model has a declared domain and an out-of-domain policy (`error`, or `warn-and-evaluate` for sensitivity studies). Nothing in the literature reviewed here authorizes silent clamping of energy, rigidity, distance or mean free path.

| Model or operation | Minimum checks |
|---|---|
| q-form, Eq. (12) | 1 < q < 2 and H ≥ 0 as printed; D₀ > 0; −1 ≤ μ ≤ 1. I(q, 0) is finite for q < 2 although D_μμ(0) = 0 when H = 0 (derived) |
| ε-form, Eq. (14) | ε > 0 (φ diverges logarithmically as ε → 0); the printed Pacheco variant is separate (D-4); normalization checked by re-integrating Eq. (11) |
| Dröge V_A form, Eq. (16) | the λ in its prefactor and the λ given by Eq. (11) differ when V_A/v > 0 (F-PA-03) |
| Eqs. (21)–(23) | B₀, δB²_slab, k_min > 0; k_d > k_min where used; 1 < s < 2; p > 2 for the dissipation term (Lang et al.); Eq. (21) is a sum of asymptotes and exceeds the exact integral, Eq. (45), by an amount that depends on s: up to 26% for s = 5/3 (exact/form = 0.7938 at R = 3.03; F-QLT-07), 27% for s = 1.65 and 35% for s = 1.5 (derived); the comparison is made at the configured s |
| Dynamical electron terms | V_A > 0 and α_D > 0 where a = v/(α_DV_A); the magnetostatic limit α_D → 0 is a separate limiting case that this document does not specify |
| Zank fit, Eqs. (24)–(25) | B, slab variance, correlation length and rigidity positive; variance convention named (U-12); stated validity 10 MV–10 GV (Chhiber et al., Perri et al.) |
| SEP-CHEN24, Eq. (6) | fit domain about 0.1–0.8 AU and 100 keV–1 GeV as stated; no automatic extrapolation; the printed coefficient uncertainties are not treated as independent without a published covariance |
| Shock prescriptions, Eqs. (27)–(31) | side of the shock and sign of the distance checked; λ > 0 and I_w > 0 in Eq. (27) require (u₁ − V_A)(x + x₀) > 0; the physical upstream case u₁ > V_A, x + x₀ > 0 is required and checked |
| GCR broken power laws | positive normalization, rigidity reference and field; evaluation in logarithms for extreme rigidity ratios; UNIT and VALUE normalizations kept apart |

**He & Wan focusing correction.** With x = λ∥,0/L, Eq. (18) is

$$
\frac{\lambda_\parallel}{\lambda_{\parallel,0}}=\frac{3(x-\tanh x)}{x^3}
=1-\frac{2x^2}{5}+\frac{17x^4}{105}-\frac{62x^6}{945}+\frac{1382x^8}{51975}-\frac{21844x^{10}}{2027025}+\dots
\tag{43}
$$

(derived from the Taylor series of tanh; convergent for |x| < π/2). The difference x − tanh x cancels in floating point: its relative rounding error is about 3ε/x² for machine precision ε, i.e. complete loss of accuracy in double precision near x = 10⁻⁸. The series through x⁸ has a truncation error below 21844x¹⁰/2027025, which is 1.1 × 10⁻¹⁵ at |x| = 0.05 (alternating series with decreasing terms), while the direct formula is accurate to about 3 × 10⁻¹³ there; evaluating the series for |x| ≤ 0.05 and the direct formula above is therefore accurate to better than 10⁻¹² for all x up to about 10¹⁰⁰, beyond which x³ overflows in double precision and 3(1 − tanh x/x)/x² is used instead (derived; fixture F-NUM-01). This changes only the evaluation, not the published closure.

**Zank auxiliary functions.** In Eq. (25), compute A = expm1((5/6) log1p(s²)) and the denominator of q as s² − expm1((1/6) log1p(s²)), and use the limit q → 2 as s → 0 (derived) instead of evaluating 0/0. The bracket of Eq. (24) tends to 1 as s → 0 (F-NUM-01).

**Exact slab-QLT integrals (validation references).** For the bendover spectrum of Section 7.1 with index ν (s = 2ν, 1/2 < ν < 1) and magnetostatic slab QLT, the substitution t = μ² turns Eq. (11) into Euler's integral for the Gauss hypergeometric function (derived):

$$
\lambda_\parallel^{\rm exact}=\frac{3}{2\pi C(\nu)(2-s)(4-s)}\left(\frac{B}{\delta B_{\rm slab}}\right)^2R_L^{\,2-s}\,l^{\,s-1}\;
{}_2F_1\!\left(-\nu,\,1-\nu;\,3-\nu;\,-\frac{R_L^2}{l^2}\right),
\tag{44}
$$

with δB²_slab the total slab variance and l the bendover length (λ_s = 2πC(ν)l). The prefactor is the inertial-range limit, so the hypergeometric factor is the ratio of the exact result to that limit; for ν = 5/6 it is the integral that the bracket of Eq. (24) approximates (F-QLT-02). For large R_L/l the Pfaff transformation ₂F₁(a, b; c; z) = (1 − z)^{−a}₂F₁(a, c − b; c; z/(z − 1)) gives the convergent form (1 + x²)^ν ₂F₁(−ν, 2; 3 − ν; x²/(1 + x²)) with x = R_L/l.

For the TS2003 spectrum without dissipation range (Section 7.1: flat below k_min, k^{−s} above) the same integral splits at μ = 1/R and gives (derived)

$$
\lambda_\parallel^{\rm exact}=\frac{3s\,R_L^2k_{\min}}{\pi(s-1)}\left(\frac{B_0}{\delta B_{\rm slab}}\right)^2\tilde J(R),\qquad
\tilde J(R)=
\begin{cases}
\dfrac{2R^{-s}}{(2-s)(4-s)}, & R\le 1,\\[1.5ex]
\dfrac14+\left(\dfrac{1}{2-s}-\dfrac12\right)R^{-2}+\left(\dfrac14-\dfrac{1}{4-s}\right)R^{-4}, & R>1,
\end{cases}
\tag{45}
$$

with R = R_Lk_min. Eq. (21) is the same expression with J̃ replaced by 1/4 + 2R^{−s}/[(2−s)(4−s)], the sum of the two limits. For s = 5/3 the ratio exact/Eq. (21) is 72/79 = 0.911392 at R = 1 and has its minimum 0.793821 at R = 3.0284 (F-QLT-07, F-NUM-01). Equations (44) and (45) are validation references; they do not replace the published Eqs. (24) and (21).

### 21.5 Pitch-angle operator in particle solvers

For the operator ∂_μ(D_μμ∂_μf), the Itô stochastic differential equation whose forward (Fokker–Planck) equation is that operator is (derived)

$$
d\mu=\frac{\partial D_{\mu\mu}}{\partial\mu}\,dt+\sqrt{2D_{\mu\mu}}\;dW_t ,
\tag{46}
$$

with W_t a Wiener process. If the host writes its operator as ∂_μ[(D_μμ/2)∂_μf] (the EPREM form, Section 6.6), the effective coefficient is D_μμ/2, and both the drift and the noise of Eq. (46) use it. Other discretizations and boundary algorithms are admissible only if they reproduce the same operator, keep −1 ≤ μ ≤ 1, conserve the particle number under no-flux boundaries, and converge to the analytic limits below. An implementation never changes the printed coefficient by a factor of two silently.

**Analytic checks (derived).** For D_μμ = (ν₀/2)(1 − μ²), Eq. (11) gives λ∥ = v/ν₀ and κ∥ = v²/(3ν₀). Under pure scattering, ⟨μ(t)⟩ = ⟨μ(0)⟩e^{−ν₀t} and ⟨μ²(t)⟩ = 1/3 + (⟨μ²(0)⟩ − 1/3)e^{−3ν₀t}; no boundary condition is needed because D_μμ(±1) = 0. For the printed (1−μ)² variant (Section 6.2), D_μμ(−1) = 4D₀(1 + H) ≠ 0 and the no-flux condition at μ = −1 must be imposed explicitly. For the q-form (1 < q < 2) and the M-FLAMPA |μ|^{2/3} shape, ∂_μD_μμ contains (q − 1)D₀(1 − μ²)|μ|^{q−2} sgn μ, which diverges (integrably) at μ = 0 for every H ≥ 0; for the ε-form it is discontinuous at μ = 0; and at H = 0 also D_μμ(0) = 0 (derived). The drift of Eq. (46) is therefore singular or discontinuous at μ = 0 for these shapes, so time-step and grid convergence must be demonstrated and the numerical treatment stated for every H. An undocumented positive H is not inserted, and it would not remove the drift singularity. The isotropic direction-reset model of full-orbit codes (Section 6.6) is not a pitch-angle SDE and has its own statistical tests.

### 21.6 Turbulence and wave-spectrum provider

A provider gives the closure a timestamped, immutable local snapshot with what the selected closure needs: at least B and, where required, density, Alfvén speed, slab variance, slab/2D split with its definition (Section 11.2), correlation or bendover length with its convention, spectral indices, dissipation onset, and a resolved spectrum if the closure uses one. Missing items stay missing; the library does not manufacture them from the mean field.

A wavenumber-resolved provider also states: the wavenumber grid and its unit (rad m⁻¹ or cycles m⁻¹), the representation of positive and negative k, one- or two-sided spectral density, magnetic variance or energy-density normalization, propagation direction and polarization if provided, the frame of any frequency, its time-interpolation policy, and the resonance convention of PARALLEL rev. 1.4, Section 8. Converting Elsässer energy to magnetic variance needs the density and the Alfvén ratio, or another documented closure (Eq. 40). An evolving wave spectrum, a prescribed static spectrum and an empirical D_μμ shape are alternative providers; they are not combined multiplicatively by default.

Update order: record the snapshot ID, evaluate the transport coefficient for that snapshot, and publish the coefficient snapshot before the particle step uses it. Interpolation between snapshots follows the host's background scheduler if it has a validated one, and otherwise an interpolation policy declared in the configuration; it is never inferred from a time column in the data. Where a source does not define wave growth, saturation or coupling to particle streaming, the feature is "not specified" and is not added.

### 21.7 Shock-region interface

A shock-region closure needs a shock locator with its epoch, normal direction, signed normal distance, front speed and upstream/downstream classification, plus any shock or wave quantity its published formula uses. In Eqs. (27)–(29), x ≥ 0 is the upstream and d ≥ 0 the downstream distance along the chosen shock-normal convention; a heliocentric distance, an unsigned distance to a surface, or a distance to a CME centre is not a substitute.

The runtime keeps four contributions separate: ambient scattering, shock-region scattering, self-generated-wave scattering, and numerical floors. Their priority and blending are not determined by the sources reviewed here and are set explicitly by the host configuration. A triangulated shock surface supplies geometry only, not a wave spectrum or an injection efficiency. SHOCK-PARASOL stays blocked until U-6 is resolved. The M-FLAMPA floor, Eq. (31), is a numerical floor on κ∥, not a physical Bohm bound; the corresponding λ floor is 3D_min/v (derived), which depends on energy and species.

### 21.8 GCR-specific requirements

A published GCR parameter set may require solar-activity indices, magnetic polarity, tilt angle and an epoch to be selected; a time-labelled parameter is never applied to another epoch silently, and a missing time-dependent value is a configuration error, not a reason to take a nearby year. Perpendicular and drift parameters in Section 9 are companion inputs, not outputs of a parallel-MFP form. A solver that needs the full tensor obtains every component and applies Eq. (32), including the off-diagonal terms, rather than building a scalar radial coefficient from λ_r.

**Duan et al. (2025) printed law.** With x = R/R_k and c > 0, the logarithmic slope of k₂ (Section 9.5) is (derived)

$$
\frac{d\ln k_2}{d\ln x}=a+(b-a)\,\frac{y}{1+y},\qquad y=x^{(b-a)/c}.
\tag{47}
$$

For b > a the low- and high-rigidity slopes are a and b; for b < a they are b and a. The paper's prose states a below R_k and b above without this condition. The printed function is reproduced only under its variant name, both orderings are tested (F-GCR-08, F-GCR-09), and the equation is not rewritten to match the prose (D-33).

**Other printed normalizations.** HelMod's K₀ is used with the numerical-use convention of D-21: the printed number is inserted and the result is in AU² s⁻¹, the printed unit is kept as metadata, and K₀ is not rescaled as though the dimensionless bracket carried a factor of GV. The Tomassetti et al. (2023) printing of λ∥ (D-23) and the Jiang et al. (2023) K₀ unit are not used without an explicit configuration decision. The Parker field normalization follows D-34 and U-13.

### 21.9 Validation layers and the shipped checks

The 34 fixture groups test equations and numbers; they do not show that a transport code reproduces SEP time profiles, GCR modulation, wave–particle feedback or shock acceleration. Four layers are kept apart:

1. **L0 — provenance, schema and algebra:** the fixtures of Section 13, unit checks, source keys, parsing rules and error paths, and the independent mathematical checks shipped with the bundle; both orderings for Duan (F-GCR-08, F-GCR-09) and both variance conventions for Zank.
2. **L1 — numerical operators:** scattering statistics of Eq. (46), amplitude ↔ λ round trips for the q-, ε- and Dröge forms, conservation at the μ boundaries, small-argument stability (F-NUM-01), comparison of Eqs. (21) and (24) with Eqs. (45) and (44), spectral normalization, and convergence in time step, quadrature and grid.
3. **L2 — runtime coupling:** changes of B, turbulence, species, energy and time; background interpolation; deterministic provider selection; signed shock distances; provenance; coefficients independent of the process decomposition; controlled refusal for missing input. At least one fixed-background case precedes time-dependent tests.
4. **L3 — science comparisons:** published parameter sets over their stated domains; model–observation comparisons of like quantities (λ∥ with λ∥, λ_r with λ_r, the same energy definition); transport experiments where reference simulations or observations exist. These are validation evidence, not a licence to tune parameters until they agree.

Algebraic L0 checks use a relative tolerance of 10⁻¹² where the conditioning allows, quadrature checks 10⁻⁹ (Section 13). Near-singular integrals, very small values, Monte Carlo statistics and transport observables need tolerances that are justified and recorded separately (convergence tests, sampling uncertainty, independent solutions). A failed criterion is not loosened to report a pass. Each test reports PASS, FAIL or SKIP with a reason, and SKIP is never counted as PASS.

The bundle contains two scripts. `reference_verification.py` compares every numerical value in `benchmark_points.json` with an expected value that it computes without the generator: derived values are recomputed with Eqs. (44) and (45), the closed form of I(5/3, H) and hypergeometric forms of the Dröge integral in place of numerical quadrature, and Eq. (47) and its NWU analogue in place of numerical differentiation, with the GCR inputs read from the data files; the row coordinates (energies, radii, rigidities) and the printed or exact constants are compared with the values stated in this document. It uses the tolerances of Section 13 and checks the provenance fields of every data file, the runtime states and the digests. `independent_physics_checks.py` compares the identities of this document with numerical quadrature, numerical differentiation and a seeded simulation of Eq. (46): the q-form normalization, the EPREM factor, Eqs. (44) and (45) including their prefactors (by integrating Eq. (11) with the QLT D_μμ of Section 6.8 and the normalized spectra), Eqs. (43) and (47), the moment decay of Section 21.5, the tensor projection of Eq. (32), the Parker-field magnitude, and the other derived identities of Sections 6–10. Neither script verifies that a transcription matches its source; that was done by reading the sources (Section 20).

### 21.10 Blockers

| Blocker | Required before the feature can be reported as implemented |
|---|---|
| U-3, U-4 | the host-code behaviour of `Chen2024AA` and `Tenishev2005AIAA` audited and documented; existing behaviour preserved until then |
| D-1/U-7, D-2/U-8, D-3/U-5, D-4 | each printed variant kept separate; no implicit correction, operator-factor change or default |
| U-6 | PARASOL's Λ(E) and Δx(E) obtained from the source and verified |
| D-29 | an unambiguous definition of Ω_n (angular or cyclic, relativistic or not) |
| U-9, U-12 | named turbulence inputs and the slab-variance convention |
| U-2 | the species of SEP-CHEN24 before any λ∥ is derived |
| U-10, U-11, U-13 | the interpretation of the D-9 printed exponents; the interpretation of λ_unspecified presets; the field normalization of Eq. (42) |
| D-33 | GCR-DUAN25 carries its printed formula as a named variant; its slopes are never described as a below and b above without the condition b > a |
| Companion specifications | the exact PARALLEL rev. 1.4 and PERPENDICULAR revisions used are pinned, and every section and equation number cited here is checked against them; a delegated path whose companion text is not available to the implementer is reported as `BLOCKED_DEPENDENCY` and left unimplemented |
| Items not found or not legible (Section 14.3) | not implemented and listed as gaps, including the sign of the Schwadron et al. (2010) radial exponent, the typeset σ_D of Breech et al. (2008) and the Duan et al. (2025) fitted values |

A stage may contain explicitly blocked optional features, but the handoff distinguishes them from features that are implemented and validated, for both the implementation record and later methods sections.

## 22. Codex implementation roadmap (final section)

This is the implementation plan for Codex. It replaces the MF00–MF10 table of version 1.0 and keeps its objectives and fixture identifiers; the stages extend the PARALLEL rev. 1.4 roadmap (its Section 23, stages PD00–PD12, of which PD01 is the shared interface and SI conversion, PD02 the power-law, broken-rigidity and Bohm models, PD03 the spectra and provider contracts, PD04 pitch-angle conversion and slab QLT, PD05 NLGCE-F, and PD08 the turbulence and wave adapters). Sections 1–21 are the authoritative reference. Physical equations are not simplified or modified to make a test pass, and missing closure parameters, wave physics, turbulence spectra or literature numbers are not invented: a blocked source-specific feature is preferable to a working but undocumented substitute.

### 22.1 Operating instructions

**Host repository.** Determine the repository layout first. Locate, by inspection, the shared particle-transport physics and every consumer of λ, κ or D_μμ that the repository contains (for example field-line, three-dimensional or GCR solvers); do not create or assume a directory path, or a solver type, before locating it. A shared location for physical formulas is preferable to duplicated code. Reuse the PARALLEL and PERPENDICULAR interfaces only after checking their revisions and symbols (Section 21.10). Do not change tested physics before its baseline behaviour is recorded.

**Project records.** At MF00, create or update a repository-local instruction file for this task, a stage-and-status ledger, and a test inventory for the mean-free-path work, following the repository's existing naming pattern if there is one. Record the specification version, the companion-data digest, the dependency revisions, the source-code revision, the active D-/U-decisions, and for each stage its status (`NOT_STARTED`, `IN_PROGRESS`, `PASS`, `FAIL`, `BLOCKED` or `SKIP`, with the reason). Update the ledger when a stage starts and ends. A completed stage records the command that reproduces its tests, the files changed, the tests run, the results and the known limitations.

**Tests.** Keep the host's existing tests and their thresholds. If the repository has a global test runner, add the new tests to it without relaxing existing gates. Each new test has its own runner, reference data, machine-readable result, log and, where useful, plots; the global summary shows progress and the final PASS/FAIL/SKIP/ERROR counts. An unavailable comparison is not reported as passing. A unit test of a local closure is distinguished from validation of a transport simulation.

**Documentation.** Code comments name the source equation or the derived identity implemented, the source key, units, domain and the relevant D-/U-items. Configuration examples name the source and the formula variant, not only a code name. Every parameter with a choice of unit or normalization has a visible input field or an unambiguous source-specific preset. Regression output records the active parameter values so that a methods section can be written from it.

### 22.2 Stages

| Stage | Work | Depends on | Exit gate |
|---|---|---|---|
| MF00 — baseline and source audit | Inspect the host repository: existing λ, κ and D_μμ implementations, the input parser, scattering types, and, if present, GCR code, turbulence provider and shock interface; record what the modes "Chen2024AA" and "Tenishev2005AIAA" compute; freeze the baseline configuration and test log | PD00 | written inventory; existing tests captured; U-3 and U-4 status recorded; original behaviour reproducible |
| MF01 — data and provenance | Extract and checksum the bundle; run `reference_verification.py` and `independent_physics_checks.py`; read-only source-key loader and the parser of Section 21.3; record the bundle digest | MF00 | both scripts pass (34 fixture groups); qualifiers preserved; rational literals, units and rejected forms tested; no physics configuration changed |
| MF02 — species and SI conversions | Typed species, mass, charge, kinetic energy (total, per nucleon), pc, rigidity, speed, γ, Larmor radius, SI boundaries | MF01, PD01 | F-KIN-01; protons, electrons, alpha particles and other charge states; round trips; invalid inputs rejected |
| MF03 — λ kinds and geometry | Tags λ∥, λ_r, λ_rr, λ_iso, unspecified; field angle from the host field; Eqs. (1)–(3) without implicit reinterpretation | MF02 | F-SEP-02, F-PERP-01; zero and nonzero perpendicular coefficient; ψ → π/2 rejected; field magnitude versus radial component (D-34) |
| MF04 — prescribed SEP closures | Eq. (9) and the Section 5.2 presets whose runtime state is READY; momentum variable and λ kind required; rational-exponent parsing; source-normalized records | MF02, MF03, PD02 | F-SEP-01, F-SEP-03, F-SEP-05; grids in energy, distance, species and charge; λ_unspecified refused (U-11) |
| MF05 — pitch-angle diffusion | q-form, ε-form, isotropic and the supported printed variants; amplitude ↔ λ by Eq. (11); `OperatorConvention`; Dröge nominal versus effective λ; He & Wan option with Eq. (43); stochastic and PDE interfaces | MF03, MF04, PD04 | F-PA-01 to F-PA-06, F-NUM-01; moment decay and λ = v/ν₀ (Section 21.5); μ-domain, no-flux and conservation tests; EPREM factor; convergence plots |
| MF06 — turbulence and spectrum providers | Static and timestamped providers; slab variance, energy units, spectral normalization and correlation-length conventions; wavenumber-resolved spectra only where the pinned PARALLEL physics exists | MF02, PD03; PD08 only for wavenumber-resolved providers | normalization integrals; Eq. (40) with a declared Alfvén ratio; invalid sources and units rejected; unchanged inputs give unchanged coefficients |
| MF07 — ion turbulence-based λ∥ | Eq. (21) with Eq. (45) as a validation reference (not a replacement); Eq. (24) with the variance convention and Eq. (44) as reference | MF06 | F-QLT-01, -02, -03, -05, -06, -07, F-NUM-01; limits and transition; small-argument stability; domain warnings |
| MF08 — electron turbulence-based λ∥ | Eq. (22) RS branch; each DT printing only under its own ID; p, k_d, V_A, α_D and domain required; no DT default | MF06, MF07 | F-QLT-04, -06; variation of p and k_d; the DT differences reproduced and labelled; α_D = 0 rejected |
| MF09 — shock-region λ and κ | Bohm, Eq. (27), Eqs. (30)–(31); tagged shock coordinates and side; ambient, shock and floor kept separate; SHOCK-PARASOL blocked until U-6 | MF03–MF08, PD02 (Bohm); host shock interface if present | F-SH-01 to F-SH-04; side and signed-distance tests; floor reported separately; documented refusal when shock information is missing |
| MF10 — GCR forms | Eqs. (33), (34) and the other fully specified Section 9 forms with unit-specific presets; K∥ and λ∥ kept distinct; tensor components by Eq. (32) | MF02, MF03, PD02; PD05 for GCR-QINSHEN17 | F-GCR-01 to F-GCR-09; field and rigidity-reference conversions, with the U-13 field normalization stated; λ ↔ K round trips; UNIT versus VALUE; Duan with b > a and b < a; HelMod D-21 convention; D-9 exponents only after U-10; no inferred units |
| MF11 — source-specific conditional models | SEP-CHEN24 κ∥ as published, with derived λ∥ only for a stated species; Minoshima Ω_n audit; gating of incomplete SEP and GCR prescriptions | MF00, MF02, MF04, MF10 | F-SEP-06 with domain gates; no `Tenishev2005AIAA` code path without documented code evidence; incomplete models listed BLOCKED, not PASS |
| MF12 — host transport integration | One coefficient interface for the consumers found at MF00; D_μμ integrated with each solver's operator convention; dimensional basis and tensor behaviour documented | MF03–MF11 as applicable | where both a field-aligned and a three-dimensional consumer exist, fixed-background runs of both give the same local λ and κ at the same state (otherwise recorded as not applicable, with the reason); baseline tests pass; coefficients independent of the process decomposition |
| MF13 — static scientific comparisons | Observational and parameter data loaded as reference only; like-for-like comparisons (λ kind, species, energy definition); curves and tables of each supported model over its stated domain | MF04–MF12 | model-versus-reference plots for SEP electrons, protons and GCR with source keys in the metadata; the Palmer band is not an absolute pass/fail criterion |
| MF14 — time-dependent coupling | Coefficients updated from frozen snapshots at declared epochs when the host provides a time-dependent background or wave spectrum; only already implemented and validated wave providers | MF06, MF09, MF12 | fixed-snapshot invariance; interpolation and stale-snapshot tests; reproducible runs; provenance of spectrum, κ and λ retained |
| MF15 — regression campaign | Per-test runners and global summary; precision and convergence sweeps; error-path tests; at least one end-to-end SEP case and one GCR case where a reference exists | MF01–MF14 | L0–L3 results reported separately (Section 21.9); logs and plots; counts by status; no thresholds weakened; out-of-domain behaviour documented |
| MF16 — handoff | Examples, comments, interface documentation and methods templates (Section 16) updated with the models, versions and units actually chosen; final validation matrix and gap register | all attempted stages | reproducibility package with commands, source revisions, digests, measured results, source mappings and the open D-/U-items |

Stage dependencies are feature-aware: a blocked optional feature (for example a DT electron printing) does not block the release of an independently validated simple SEP power law, but the inventory is then not described as fully implemented. If a required PARALLEL or PERPENDICULAR equation is not available, the module is recorded as `BLOCKED_DEPENDENCY` and left unimplemented.

### 22.3 Minimum runnable configurations

**A — prescribed SEP preset.** A proton energy grid and a position grid; a named λ∥ preset whose runtime state is READY_EXPLICIT_INPUTS; a selected q-form or ε-form with its parameter and operator convention. Output λ∥, κ∥, D_μμ(μ), the local field angle, λ_r, λ_rr only if κ⊥ is available, source labels and validity flags; plot λ∥(r, E) and D_μμ(μ) and check the normalization integral numerically.

**B — static turbulence-based closure.** Supplied B, slab variance, spectral or bendover scale and all model indices; compare Eq. (21) and Eq. (24) where their assumptions and domains overlap, without claiming equality, and each with its exact reference, Eqs. (45) and (44). For electrons, compare RS and the selected DT printings at their published parameters, labelling the spread as a source discrepancy, not an uncertainty band.

**C — shock and evolving background.** A documented shock provider with front and signed normal coordinates and any closure-specific turbulence quantities; ambient-only, shock-region and post-shock tests. Output the side classification, the active prescription, the numerical floor and the coefficient snapshots. No wave spectrum, injection efficiency or PARASOL function is inferred from shock speed or compression.

**D — GCR parameter set.** A named source, year or Bartels rotation, species and epoch; a rigidity grid; the background field and every solar-activity input the source requires. Output K∥(P), λ∥(P), and the tensor only if the perpendicular and drift inputs are available; show the rigidity break and the normalization unit. A missing time-dependent coefficient is a configuration error.

For each configuration, save (1) machine-readable input with every parameter and quantity tag, (2) machine-readable output with units and model IDs, (3) diagnostic plots, (4) test results and tolerances, and (5) a provenance and dependency manifest. Command lines are generated from the host executable's actual input grammar after MF00; no unverified option appears in an example.

### 22.4 Completion and release rules

1. Record the baseline before editing, and rerun the relevant baseline and new tests after every stage. A stage is PASS only when its gates pass without modified thresholds.
2. Every runtime formula maps to one equation of this document and its source, to an identity labelled derived here, or to an explicitly named non-default printed variant. Parameter changes are visible in the configuration and in model reports.
3. Missing information stays missing: for U-4, U-6, D-29, unresolved printed discrepancies and unread source items, the result is BLOCKED, not a guessed implementation.
4. Equation fixtures, operator convergence tests, coupled-run tests and observational comparisons are separate evidence with separate acceptance criteria.
5. Keep a final table: stage, status, source revision, files modified, command, PASS/FAIL/SKIP/ERROR, measured error, open item, next action.
6. Close with an accurate statement of which SEP, electron, ion, shock and GCR closures are ready for production, which are source-specific variants, and which are not implemented. A passing test is not an observational validation unless that comparison was actually made.

The library reproduces the selected and documented scattering physics; it does not choose the physics on the user's behalf. A defensible partial implementation with explicit blockers is preferable to an undocumented complete one.
