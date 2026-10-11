# Polynomial-Spectral Model of SEP-Excited Alfvén Waves on a Prescribed Turbulent Background

Revision 3, 10 October 2026, native SI units. This revision incorporates findings R01–R18 of the second physical review (SI edition): the finite-frequency electrodynamic factor r_σ² in the scattering rate and in every matched growth and action source (R01–R03), corrected maintenance algebra and physical excess budgets (R04–R05), finite-frequency validity thresholds and an enforcement policy (R06), a reproducible background profile with field-aligned derivatives (R07), restriction of nonlinear transfer and fast-shock classification (R08–R09), total extensive wave action as the primary transported state with corrected drift and source terms (R10–R11), conservative remapping (R12), space-time flux registers on 3D AMR and the native AMPS interfaces (R13–R16), finite-step particle–wave coupling (R17) and native SI quantities with scalable caches (R18). Section 16 is the SI/AMPS implementation roadmap. Every formula marked (checked) was verified by code as recorded in Section 15; checks of the previous revision that the corrections invalidate are withdrawn there.

## Summary

The model represents the ratio of the evolved parallel-wave spectrum to a reference spectrum, g = ln(I/I₀), as a Legendre polynomial in a magnetic-field-scaled logarithmic wavenumber coordinate, for two propagation families σ = ± and two helicity channels h = ±. The transported numerical state is the total extensive wave action of each channel; g is recovered from it to evaluate spectra. The representation guarantees a positive spectrum; accuracy is established only by convergence of scattering rates, κ∥, energy exchange and SEP observables.

- **One scattering tensor for particles and waves.** Scattering is exactly elastic in the frame of each wave family: D_μμ carries the electrodynamic factor r_σ² of a propagating Alfvén wave, D_μp = a_σD_μμ and D_pp = a_σ²D_μμ with a_σ = V_σp/(v r_σ). The growth rate and the action source are derived from this tensor by the energy identity with the same finite-frequency resonance, interaction window and species list; the exchanged energy is booked once per accepted interaction.
- **The split is exact when the reference solves the SEP-free equation.** A dynamically evolved SEP-free reference (policy A) is the first implementation. A prescribed reference (policy B) carries a residual E₀ that acts as a maintenance forcing with rate τ = (S + E₀)/I₀; it is admissible only when τ ≥ 0 or the step and domain are restricted.
- **Growth needs no wavenumber grid.** The projection of the windowed growth rate onto P_n equals a particle moment; for Monte Carlo particles and Parker solvers the moment is a contraction of per-cell kernels with pitch-angle Legendre moments or with the derived diffusive anisotropy.
- **WKB transport leaves g invariant and conserves action.** Coefficient vectors are advected at u + σV_A; in Lagrangian variables the spectral drift is −D ln(Δs B)/Dt − σ∂_s(V_AB)/B.
- **Boundaries are classified by speed relative to the actual surface**; a material outer end of a Lagrangian line is inflow for the inward family.
- **Conservative numerics.** Total action moments are transported with one oriented flux per shared face, remapped by summing extensive moments, and recovered by moment matching; wave energy is a separately controlled observable with an explicit WKB work term.
- **Nonlinear transfer is restricted.** C1 (relaxation toward a maintained background by shearing in perpendicular turbulence) is retained with a declared regime and calibration; C2 is disabled in the strict parallel slab baseline, because purely parallel transverse Alfvén waves do not interact in incompressible MHD.

Assumptions that must accompany any result: positive ions; slab reduction with the wavevector kept along B; nondispersive finite-frequency Alfvén resonance over a declared (p, μ, position, family) region with a specified treatment outside it; quasilinear scattering with δB²/B² ≲ 0.1 and |Γ| ≪ κV_A; a prescribed open-system MHD background with its residual reported; reflection omitted only where a non-WKB benchmark shows an acceptable error.

## 0. Open decisions register and stop-and-ask protocol

The physics in Sections 1–15 is fixed. The items below are choices the model leaves to the owner of the intended application; none has a silent default. An implementing agent (Codex or Claude Code) must not infer, approximate or pick any of them. It stops and asks, as specified in the protocol after the table, and keeps every stage that depends on an unresolved item BLOCKED.

| ID | Decision | Options the model admits | Affects | Status |
| --- | --- | --- | --- | --- |
| D01 | Reference policy | A: evolved SEP-free reference with the same discrete operator (recommended first); B: prescribed I₀ with E₀ computed and τ = (S+E₀)/I₀ admissible | Sections 5, 8.4, 9; S02 | OPEN |
| D02 | Energy state for the extensive moments | (a) an additional intrinsic-energy moment with its own closure parameter; (b) energy-constrained source projection preserving realizability; (c) action-only state with the energy defect reported and bounded | Sections 7, 9, 10; S03, S04 | OPEN, no default |
| D03 | Closure outside the accepted interaction region | (i) a fully specified fallback tensor for the background with its frame, resonance or broadening kernel and energy rule; (ii) block predictions where the region is excluded | Section 6; S02, event runs | OPEN |
| D04 | Accepted region for the intended event | Minimum energy per family, pitch-angle extent, species, position range, as a (p, μ, z, σ) map | Section 6; S00, S12 | OPEN, event-specific |
| D05 | Nonlinear transfer C1 | Enabled or disabled; turbulence regime; ε source; Γ_T form; C_T calibration reference | Section 8.4; S02, S12 | OPEN; C_T has no derived value |
| D06 | Reflection of the excess | Omitted after V09 passes for the declared profile and event, or extension E02 | Section 4; S12 | OPEN until V09 is run |
| D07 | Shock closure | Injection law, acceleration treatment, wave jump or one-way transmission (Vainio & Schlickeiser 1998) or none | Section 11; S12 | OPEN |
| D08 | Particle solver per host | Full focused transport; or the diffusive closure of Section 3.3 where its ordering is tested | Section 3; S05, S09 | OPEN per host |
| D09 | Frame, polarity and μ axis | μ against B or against t_out; outward tangent; host w± to (σ, h) mapping; both polarities | Sections 1, 2; S00 | OPEN |
| D10 | Reference spectrum reconstruction | Spectral shape φ(k), support and break, slab fraction f_slab, helicity split, family imbalance, from host w± or an empirical model | Sections 5, 14; S02 | OPEN |
| D11 | Tube normalization | Φ_line or seed area and its provenance for each AMPS line; characteristic-only lines | Section 10.2; S06 | OPEN |
| D12 | Plasma composition | Ion species, densities and charges entering d_i, V_A and the validity thresholds | Sections 1, 6; S02 | OPEN |
| D13 | Species feeding back on the waves | Which species enter the growth sum; others booked as external terms | Section 8.2; S01 | OPEN |
| D14 | Stochastic operator in AMPS | Itô SDE with the drifts of Section 12, or an event-based operator with demonstrated generator equivalence | Section 12; S09 | OPEN |
| D15 | Time synchronization | Wave interval; common step or AMR subcycling; particle local-step policy | Section 12; S10 | OPEN |
| D16 | Time-order target | Accept first order (backward Euler transport) or upgrade subsolvers for a second-order claim | Section 12; S04, S10 | OPEN |
| D17 | Release observables and accuracy limits | Which SEP spectra, arrival times, anisotropies, wave power and exchanged energy, with numerical error limits | Section 15; S00, S13 | OPEN; H4 is BLOCKED until set |
| D18 | Initial numerical settings | x interval, window ramp, N, Q, L, shells, μ cells, used only as starting points for convergence studies | Section 7; S03, S13 | OPEN, convergence decides |

**Stop-and-ask protocol for the implementing agent**

1. Before any implementation (stage S00), copy this register into `docs/sep_wave/decisions.md` with one line per item: the chosen option, who approved it and when, or the word OPEN.
2. An OPEN item blocks every stage, test and physical run that its "Affects" column names. Stages that do not depend on it may proceed.
3. When work reaches an OPEN item, or when an input, closure, convention or tolerance is needed that neither this document nor the register supplies, the agent stops that line of work and writes a question to `docs/sep_wave/questions.md` stating the item ID, the options with their consequences, and what it would do under each; it does not proceed on the affected stages until a person answers. It never substitutes a plausible value, a taper, a cap or a disabled term in place of the answer.
4. A recommendation in this document (for example policy A in D01) is still a decision to confirm, not a default to apply.
5. Any answer that changes physics is recorded as a dated entry in the register and in `docs/CHANGELOG.md`, and the affected gates are re-run.

## 1. Scope, conventions and units

The model covers SEP transport along field lines upstream of a fast-mode shock, or from a fixed inner boundary before a shock exists. All physical quantities are in SI; dimensionless variables (x, ξ, μ, g, p̂) are declared explicitly.

| Quantity | SI unit or normalization |
| --- | --- |
| Magnetic field B; charge q; rest mass m | T; C; kg |
| Position, arc length z, time, speeds u, v, V_A | m; s; m/s |
| Momentum p; total energy E; kinetic energy T | kg m/s; J; J |
| Number density n; mass density ρ | m⁻³; kg/m³ |
| Gyrofrequency; Alfvén speed; resonance scale | Ω = qB/(γm); V_A = B/√(μ₀ρ); k* = q_pB/p* |
| Spectral coordinate κ = \|k\|, k* | m⁻¹; x, ξ, g dimensionless |
| Magnetic power spectrum I_σh(k) | T² m; ∫I_σh dk over the channel = δB²_σh in T² |
| Wave energy per unit x; forcing S, E₀ | κI/μ₀ in J/m³; T² m/s |
| Wave action per unit x; extensive moments M_jn | 𝒜 = I/(μ₀V_A) in J s/m³; J s |
| Segment geometry | area A in m²; volume in m³; flux Φ_line in Wb |
| Distribution f; rates Γ, D_μμ | m⁻³ (kg m/s)⁻³; s⁻¹ |
| D_μp; D_pp; κ∥ | kg m/s²; kg² m²/s³; m²/s |

Constants: μ₀ = 1.25663706127×10⁻⁶ H/m, q_p = 1.602176634×10⁻¹⁹ C, c = 299792458 m/s (CODATA 2022). Particle number n = 2π∫p²f dp dμ. If a code stores p̂ = p/p_ref, then f̂ = p_ref³f, dp = p_ref dp̂ and ∂f/∂p = p_ref⁻⁴∂f̂/∂p̂; every kernel and measure is transformed once at a named unit boundary. Host Alfvén wave energy densities satisfy δB_σ² = μ₀w_σ when w_σ is the total (magnetic plus kinetic) energy density; slab and helicity fractions are applied once.

| Assumption | Condition it needs | Where it fails |
| --- | --- | --- |
| Positive ions | q > 0 so that Ω > 0 and ln(q/q_p) is defined | Electrons need a separate helicity and polarity contract |
| Slab (parallel) Alfvén waves carry the resonance; the quasi-2D background is a spectator | Slab fraction known (about 20% at 1 AU) | Oblique or compressive modes near quasi-parallel shocks |
| The wavevector is kept along B | Model restriction; the discarded perpendicular Hamiltonian component is diagnosed | Strongly curved or sheared fields |
| WKB propagation | κL_A ≫ 1 and background change slow against κV_A; non-WKB benchmark V09 | Localized gradients; the Alfvén-point neighbourhood for inward waves |
| Non-dispersive finite-frequency Alfvén branch | κ_r d_i = V_A/(γ\|vμ − σV_A\|) ≤ 0.1 at every accepted interaction, with d_i for the declared composition | Low energies and small \|μ\| near the Sun (Section 6) |
| Quasilinear theory, random phases | Resonant δB²/B² ≲ 0.1 over the occupied support; \|Γ\|/(κV_A) ≪ 1 | Strongly amplified foreshocks |
| Reference spectrum | Policy A evolved SEP-free; or policy B with E₀ reported and τ admissible | Preceding events that amplified the resonant range |
| Particle transport order | First order in u/v and u/c; Du/Dt neglected; stated along the profile and near shocks | Strongly accelerating or sheared flow |
| Fast-mode shock on the connected field line | u_n,upstream > c_f(θ_Bn) with c_f² = ½[V_A² + c_s² + √((V_A² + c_s²)² − 4V_A²c_s²cos²θ)] | Low corona before shock formation; weak flanks |

Sign conventions: σ = +1 propagates away from the Sun in the plasma frame, σ = −1 toward it; the outward tangent t_out of the line is defined independently of polarity, and whether μ is measured against B or against t_out is declared together with the family labels, resonant helicity and host w± mapping. The helicity channel h = ±1 is the one resonant with particles of sign(vμ − σV_A) = h, which for positive ions is magnetic helicity −h. Coriolis terms are neglected for the transverse perturbations.

## 2. Geometry and notation

All quantities live on one field line with arc length z measured along t_out. The 1D formulas are written in the frame where the background flow is field-aligned; Section 10 gives the Lagrangian and 3D forms.

- **Background.** Field-aligned flow u, field B, tube area A = Φ_line/|B|, density ρ with ρuA = const, V_A = B/√(μ₀ρ), focusing length L_B = −B/(∂B/∂z).
- **Spectrum.** I_σh(z, k, t) on signed k, with |k| = κ; the channel integral ∫I_σh dk = δB²_σh. The plasma-frame wave energy density of a channel is δB²_σh/μ₀ (magnetic plus kinetic); per unit x = ln κ the energy is κI/μ₀ and the action is 𝒜 = I/(μ₀V_A).
- **Particles.** f(z, p, μ, t) with momentum in the local plasma frame; Ω = qB/(γm).

Resonance and the log coordinate, with V_σ = σV_A and r_σ = 1 − μV_σ/v (checked):

$$
k_r(\mu,p)=\frac{\Omega}{v\mu-V_\sigma},\qquad \kappa_r=\frac{\Omega}{|v\mu-V_\sigma|},\qquad x=\ln\frac{\kappa}{k_*},\qquad k_*=\frac{q_pB}{p_*}
$$

$$
x_r(\mu,p)=\ln\!\left[\frac{q}{q_p}\,\frac{p_*}{p\,|\mu-V_\sigma/v|}\right]
$$

B cancels from x_r; the map from (μ, p) to the polynomial coordinate depends on the cell through V_A/v and is evaluated per cell. Equal-rigidity species share resonances only in the zero-frequency limit; with finite frequency each species keeps its own V_A/v in x_r.

## 3. Particle side

### 3.1 Focused transport equation

To first order in u/v and u/c, with plasma-frame momentum and the flow acceleration Du/Dt neglected (Roelof 1969; Ruffolo 1995; Isenberg 1997):

$$
\frac{\partial f}{\partial t}+(u+\mu v)\frac{\partial f}{\partial z}+\dot\mu\,\frac{\partial f}{\partial\mu}+\dot p\,\frac{\partial f}{\partial p}=\sum_{\sigma=\pm}\mathcal{L}_\sigma f+Q
$$

$$
\dot\mu=\frac{1-\mu^2}{2}\left[\frac{v}{L_B}+\mu\left(\nabla\!\cdot\!\mathbf u-3\,\hat b\hat b\!:\!\nabla\mathbf u\right)\right],\qquad \dot p=-p\left[\frac{1-\mu^2}{2}\left(\nabla\!\cdot\!\mathbf u-\hat b\hat b\!:\!\nabla\mathbf u\right)+\mu^2\,\hat b\hat b\!:\!\nabla\mathbf u\right]
$$

In the flux tube, b̂b̂:∇u = ∂u/∂z and ∇·u = ∂u/∂z + u/L_B. Limits the implementation must reproduce: static focusing changes μ and not p; homogeneous isotropic compression gives ṗ = −(p/3)∇·u, i.e. p ∝ ρ^(1/3). The neglected ordering is quantified along the supplied profile and near shocks (test P5); where it is not small the host code's complete focused-transport coefficients are used. Q is the injection source with its normalization, angular shape, species, history and frame stated as inputs.

### 3.2 Scattering tensor, exactly elastic in the wave frame

Scattering by family σ is elastic in the frame moving at V_σ along B. For V_A ≪ c the scattering surfaces in plasma-frame variables are E − V_σpμ = const; their tangent gives the diffusion direction and the rank-one tensor (checked):

$$
a_\sigma=\frac{V_\sigma\,p}{v\,r_\sigma},\qquad \mathcal{X}_\sigma f=\frac{\partial f}{\partial\mu}+a_\sigma\frac{\partial f}{\partial p},\qquad D_{\mu p,\sigma}=a_\sigma D_\sigma,\quad D_{pp,\sigma}=a_\sigma^2D_\sigma
$$

$$
\mathcal{L}_\sigma f=\frac{\partial}{\partial\mu}\Big(D_\sigma\,\mathcal{X}_\sigma f\Big)+\frac{1}{p^2}\frac{\partial}{\partial p}\Big(p^2\,a_\sigma\,D_\sigma\,\mathcal{X}_\sigma f\Big)
$$

The quasilinear pitch-angle coefficient for a propagating parallel Alfvén wave contains the electrodynamic factor r_σ² in addition to the finite-frequency resonance (Schlickeiser 1989; Ng et al. 2003, Appendix B; Thomas & Pfrommer 2019). It follows from the magnetostatic rate in the wave frame transformed to the plasma frame: (∂μ/∂μ′)² = (v′/v)²r_σ² and 1 − μ′² = (v/v′)²(1 − μ²), so

$$
D_\sigma(\mu,p)=\frac{\pi}{2}\,\Omega\,(1-\mu^2)\,\frac{\kappa_r\,I_{\sigma h}(k_r)}{B^2}\,r_\sigma^2,\qquad h=\mathrm{sign}(v\mu-V_\sigma)
$$

For an unpolarized one-sided magnetic spectrum I⁽¹⁾ the prefactor becomes π/4 and r_σ² remains; for V_A/v → 0 the magnetostatic rate is recovered. At V_A/v = 0.1 and μ = ±0.9, r_σ² = 0.83 and 1.19 for outward waves, so the factor acts at the same order as the other finite-frequency terms retained here. The previous revision's rate without r_σ² is withdrawn. A distribution that depends only on E − V_σpμ satisfies 𝒳_σf = 0 and is not scattered; the operator conserves particle number; and because v a_σ = V_σ(p + μa_σ) = V_σp/r_σ, the energy change equals V_σ times the parallel-momentum change exactly when momentum-boundary fluxes vanish or are booked (checked):

$$
\frac{d}{dt}\int E f\,d^3p\,\Big|_\sigma=V_\sigma\frac{d}{dt}\int p\mu f\,d^3p\,\Big|_\sigma=-2\pi V_\sigma\int p^2\,\frac{p}{r_\sigma}\,D_\sigma\,\mathcal{X}_\sigma f\,dp\,d\mu=-\frac{1}{\mu_0}\int\Gamma_\sigma I_\sigma\,dk
$$

The total pitch-angle coefficient is D = D₊ + D₋. The accepted-interaction policy of Section 6 is applied before any resonance division; interactions at vμ = V_σ are rejected, not evaluated and multiplied by zero.

### 3.3 Diffusive (Parker) closure

A host code carrying only f₀(p) is coupled through the leading diffusive balance. With f = f₀ + h, zero pitch-angle average of h, D = D₊ + D₋ and Ā = (D₊a₊ + D₋a₋)/D (checked; Skilling 1975; Borovikov et al. 2019, Section 3):

$$
\frac{\partial h}{\partial\mu}=-\frac{v(1-\mu^2)}{2D}\frac{\partial f_0}{\partial z}-\bar A\,\frac{\partial f_0}{\partial p},\qquad
J=-\kappa_\parallel\frac{\partial f_0}{\partial z}-\frac{pU_{sc}}{3}\frac{\partial f_0}{\partial p},\qquad \kappa_\parallel=\frac{v^2}{8}\int_{-1}^{1}\frac{(1-\mu^2)^2}{D}d\mu,\qquad U_{sc}=\frac{3v}{4p}\int_{-1}^{1}(1-\mu^2)\bar A\,d\mu
$$

$$
\frac{\partial f_0}{\partial t}+u\frac{\partial f_0}{\partial z}+\frac{1}{A}\frac{\partial(AJ)}{\partial z}-\frac{p}{3}\frac{1}{A}\frac{\partial(Au)}{\partial z}\frac{\partial f_0}{\partial p}=\frac{1}{p^2}\frac{\partial}{\partial p}\Big\{p^2\Big[-\frac{pU_{sc}}{3}\frac{\partial f_0}{\partial z}+D_{pp,\mathrm{eff}}\frac{\partial f_0}{\partial p}\Big]\Big\}+Q_0,\qquad D_{pp,\mathrm{eff}}=\frac{1}{2}\int_{-1}^{1}\frac{D_+D_-(a_+-a_-)^2}{D}d\mu
$$

For momentum-independent U_sc this is Parker's equation convected at u + U_sc; D_pp,eff ≥ 0 vanishes for a single family. These are leading diffusive results valid under a tested ordering of the anisotropy (test P5), not a replacement for focused transport at arbitrary anisotropy. The growth kernels take the derived anisotropy 𝒳_σf ≃ −v(1−μ²)/(2D)∂_zf₀ + (a_σ − Ā)∂_pf₀ separately per helicity channel.

## 4. Wave kinetic equation for the total spectrum

Wave action per unit volume and unit k, N_σh = I_σh/(μ₀κV_A), is conserved along rays absent sources (Dewar 1970). In conservative form for a time-dependent background and a deforming tube:

$$
\frac{\partial(AN)}{\partial t}+\frac{\partial(Ac_\sigma N)}{\partial z}+\frac{\partial(\dot k\,AN)}{\partial k}=\frac{A}{\mu_0\kappa V_A}\Big(\Gamma I+\mathcal{N}[I]+\mathcal{R}[I]+S\Big),\qquad c_\sigma=u+\sigma V_A,\qquad \dot k=-k\frac{\partial c_\sigma}{\partial z}
$$

The ray equations are Hamilton's equations and hold for time-dependent backgrounds. For a steady background ω = κc_σ is conserved along a ray and A c_σ²E_ω/V_A = const (Jacques 1977).

**Growth and damping by particles.** The energy identity of Section 3.2 evaluated per mode gives the intensity growth rate of channel (σ, h), with the SI prefactor μ₀π²V_A/B² (checked: closes the energy budget to 2×10⁻⁶ for a relativistic anisotropic f and helicity-asymmetric spectrum at V_A/v up to 0.1):

$$
\Gamma_{\sigma h}(k_h)=\sigma\,\frac{\mu_0\pi^2V_A}{B^2}\sum_{\mathrm{species}}\int dp\;p^3\,(1-\mu_k^2)\,\frac{\Omega^2\,r_\sigma(\mu_k)}{\kappa\,v}\;\mathcal{X}_\sigma f(p,\mu_k),\qquad k_h=h\kappa,\quad \mu_k=\frac{\Omega/k_h+V_\sigma}{v},\quad|\mu_k|\le1
$$

For a Compton–Getting power law f ∝ p^(−α)(1 + αμu_str/v) it reduces, to first order in u_str/v and V_A/v, to Γ₊ = (π/2)Ω_i(n_>(p_k)/n_i)(u_str/V_A − 1)(α−3)/(α−2), an intensity rate (Kulsrud & Pearce 1969; Kulsrud & Cesarsky 1971; Lee 1983; Ng et al. 2003). Tests of this limit alone cannot detect the r_σ² factor.

**Windowed rates.** With a species-, momentum- and pitch-angle-dependent validity window W (Section 6), the factorization W(k)Γ(k) is not valid. Define Γ_W by inserting W into each species contribution of the integrand above, and Γ_{1−W} with 1 − W. Then the particles scatter on I_eff = WI + (1 − W)I₀, the evolved spectrum receives the energy source density IΓ_W/μ₀ per unit k (J m⁻² s⁻¹), and the reference reservoir receives I₀Γ_{1−W}/μ₀; their sum balances the particle exchange with I_eff, including the ramp.

**Other damping.** Thermal-ion cyclotron absorption is not included; accepted interactions are restricted to κ_r d_i ≤ 0.1 for the declared composition.

**Reflection.** The linear non-WKB flux-tube equations with z_± = δv ∓ δB/√(μ₀ρ) (checked; Heinemann & Olbert 1980; Velli 1993):

$$
\frac{\partial z_\sigma}{\partial t}+c_\sigma\frac{\partial z_\sigma}{\partial z}-\frac{c_{-\sigma}}{4}\frac{\partial\ln\rho}{\partial z}z_\sigma-\frac{c_{-\sigma}}{2}\frac{\partial\ln V_A}{\partial z}z_{-\sigma}=0,\qquad
\varepsilon_{R,\sigma}=\frac{|c_\sigma|}{4\kappa V_A}\left|\frac{\partial\ln V_A}{\partial z}\right|
$$

ε_R is the amplitude ratio of the local asymptotic particular solution for κL_A ≫ 1; ε_R² is a local energy estimate and not a bound on the reflected flux over a finite domain. The excess is evolved without reflection only where the non-WKB benchmark V09 (declared profile, incoming data, a localized gradient, the Alfvén-point neighbourhood, opposite-family effects on drift and transfer) shows an acceptable error. The conditional profile below (Section 14) uses the field-aligned derivative ∂/∂z = cos ψ ∂/∂r and u∥ = u_r sec ψ; the previous table used the radial derivative and is withdrawn.

| r | u_r (km/s) | u∥ (km/s) | V_A (km/s) | L_A,∥ (R☉) | ε_R, 1 MeV | ε_R, 100 MeV | ε_R, 1 GeV | κ_r d_i, 0.1 MeV | κ_r d_i, 1 MeV |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 3 R☉ | 74 | 75 | 886 | 8.33 | 3.5×10⁻⁷ | 3.8×10⁻⁶ | 1.5×10⁻⁵ | 0.25 | 0.068 |
| 10 R☉ | 348 | 349 | 574 | 11.82 | 4.2×10⁻⁶ | 4.5×10⁻⁵ | 1.7×10⁻⁴ | 0.15 | 0.043 |
| 64.5 R☉ | 399 | 417 | 99 | 73.83 | 8.9×10⁻⁵ | 9.2×10⁻⁴ | 3.5×10⁻³ | 0.023 | 0.007 |
| 215 R☉ | 400 | 568 | 41 | 616.96 | 2.5×10⁻⁴ | 2.6×10⁻³ | 9.9×10⁻³ | 0.009 | 0.003 |

ε_R and κ_r d_i are for σ = +, μ = 1, protons, with the finite-frequency resonance. 215 R☉ = 1.4958×10¹¹ m is the draft's outer radius, not exactly 1 AU.

## 5. The split: reference spectrum, residual and the ratio equation

Write the wave operator as 𝓛_wI = ΓI + 𝒩[I] + ℛ[I] + S. For a reference I₀ define the residual E₀ = 𝓛_wI₀ − 𝒩[I₀] − ℛ[I₀] − S (T² m/s). Policy A: I₀ is evolved SEP-free with the same discrete operator, geometry and quadrature as I, so E₀ = 0 within numerical error and g = 0 is a discrete equilibrium. Policy B: I₀ is prescribed and E₀ is the computed forcing that maintains it. With the common forcing S + E₀ in both equations, I = I₀e^g obeys exactly

$$
\frac{\partial g}{\partial t}+c_\sigma\frac{\partial g}{\partial z}+\dot x_\sigma\frac{\partial g}{\partial x}=\Gamma_W+\frac{\mathcal{N}[I]-e^{g}\mathcal{N}[I_0]}{I}+\frac{\mathcal{R}[I]-e^{g}\mathcal{R}[I_0]}{I}+\frac{S+E_0}{I_0}\left(e^{-g}-1\right)
$$

$$
\dot x_\sigma=-\frac{\partial c_\sigma}{\partial z}-c_\sigma\frac{\partial\ln B}{\partial z}-\frac{\partial\ln B}{\partial t}\ \text{(Eulerian)},\qquad
\dot x_\sigma=-\frac{D\ln(\Delta s\,B)}{Dt}-\frac{\sigma}{B}\frac{\partial(V_AB)}{\partial s}\ \text{(Lagrangian)}
$$

- The amplitude terms cancel because I and I₀ obey the same transport operator; the drift forms were checked symbolically against each other, with D ln Δs/Dt = b̂·(∇u)·b̂ and Δs B ∝ B²/ρ in ideal MHD.
- The maintenance term relaxes g toward 0 at the rate τ = (S + E₀)/I₀ only when τ ≥ 0. For an arbitrary reconstructed background τ can be negative, in which case the affine update can drive I to zero in finite time and the logarithm cannot preserve positivity; such a case requires an admissible maintenance prescription, a restricted step and domain with a stated meaning, or policy A.
- The ratio-equation terms are not the physical excess budget. Subtracting the total and reference intensity equations, which share S + E₀, gives for C1 the physical excess equation 𝓛_wδI = Γ_WI − Γ_TδI with δI = I − I₀; the C1 contribution to the physical excess energy is −Γ_TδI, and E₀ is a reservoir power on the total and reference accounts, not an additional excess source.
- Reflection of the excess is dropped consistently from the g equation and from E₀; a reference model that retains non-WKB reflection contributes it to E₀ under policy B.

## 6. Validity domain and its enforcement

Every accepted interaction (species, p, μ, cell, family) satisfies the inequalities below; interactions outside them are not evaluated by the resonant tensor.

| Condition | Inequality | Protons at V_A = 8.86×10⁵ m/s, κ_r d_i ≤ 0.1 |
| --- | --- | --- |
| Non-dispersive branch | γ\|vμ − σV_A\| ≥ 10 V_A, i.e. κ_r d_i ≤ 0.1 with d_i = V_A/Ω_p for the declared composition | μ = 1: T ≥ 7.942×10⁻¹⁴ J (0.4957 MeV) for σ = +, 5.316×10⁻¹⁴ J (0.3318 MeV) for σ = −; μ = 0.05: 2.995×10⁻¹¹ J (186.9 MeV) and 1.941×10⁻¹¹ J (121.2 MeV) |
| Quasilinear amplitude | ∫κI dx/B² ≤ 0.1 summed over the four channels on the occupied support; \|Γ\| ≪ κV_A | Run-time diagnostic |
| Transport ordering | u/v, u/c and the neglected Du/Dt terms small along the profile | Run-time diagnostic |

The thresholds are family-dependent and exchange under μ → −μ; when both families are accepted the stricter value applies. They establish only the dispersion condition: at T = 0.5 MeV and μ = 0.05, κ_r d_i = 2.23 (σ = +) and 0.64 (σ = −), so the narrow margin at μ = 1 does not admit all pitch angles. The zero-frequency ordering V_A/(v|μ|) ≪ 1 is no longer required by the tensor; an additional pitch-angle cutoff is a numerical operating choice whose neglected error is stated, not a proof.

**Outside the accepted region** the model supplies a declared closure, never an implicit extrapolation: (i) a fully specified fallback tensor for the background, with its wave frame or scattering mechanism, resonance or broadening kernel and energy-transfer rule, booked against the reference reservoir — a purely plasma-frame pitch-angle closure exchanges no energy locally; or (ii) the run leaves its declared domain and predictions there are blocked. Transported excess at excluded κ is booked as "excess outside the accepted region"; whether it remains a valid mode with another dispersion law, enters an unresolved reservoir or invalidates the domain is stated per case. The 43% small-|μ| share of κ∥ for an unbroken magnetostatic Kolmogorov slab is a sensitivity illustration only; κ∥ sensitivity is recomputed with the full tensor and the actual background band and closure.

**Window.** The accepted region is encoded as a smooth C¹ window W(p, μ, species, cell) whose ramp (width about 0.5 in x_r) is inside the species integrals of Γ_W and of the particle exchange (Section 4). Angular quadrature splits at μ = V_σ/v and at window edges.

## 7. Polynomial representation and recovery

Each channel's g is a Legendre series on one global interval, g = Σ_{n≤N} a_nP_n(ξ), ξ = (2x − x_min − x_max)/Δx. Uniform weight in log κ, diagonal mass matrix Δx δ_nm/(2n+1), closed-form derivative and endpoint values P_n(±1) = (±1)ⁿ make it the natural basis; monomials are ill-conditioned beyond N ≈ 6.

- **Interval.** Global [x_min, x_max] because k* ∝ B; outside it g = 0. The basis interval and the background generation are part of the state metadata. If regridding changes the reconstructed B while physical κ is retained, x_new = x_old − ln(B_new/B_old) with a spectral remap and endpoint accounting; this coordinate change is distinct from the ray-equation drift and is not applied twice.
- **State.** The transported state is the total extensive action M_jn = ∫_cell dV ∫𝒜P_n dx of each channel (Section 10); a_n and kernel tables are rebuildable caches. g is recovered by moment matching on the exponential closure (Levermore 1996; Junk 1998), with residual, condition number, quadrature refinement and realizability of the implied positive total spectrum recorded; an infeasible target triggers the conservative fallback and the ledger, never a clipped negative spectrum.
- **Energy.** Wave energy ∫κV_A𝒜 dx is not a finite combination of Legendre action moments; it is carried as an independently controlled energy state, or enforced by a constrained source projection, with the projection defect reported (Section 9).
- **Degree and quadrature.** N ≈ 8–12 and Q are starting choices; the accepted values follow from convergence of D, κ∥, the exchange rate and SEP observables (test V12). The same quadrature nodes are used in moments, budgets and recovery.
- **Evaluating D_σ.** x_r depends on the cell through V_A/v; the table P_n(ξ_r(μ_i, p_j)) and the window are evaluated per cell or interpolated on a documented V_A grid keyed by every dependence (species, V_A/v, support, window, background generation). D_σ(μ_i, p_j) = (π/2)Ω_j(1−μ_i²)κ_rI₀[1 + W(e^g − 1)]r_σ²/B², with g from one matrix–vector product; Monte Carlo particles evaluate the recurrence at their own ξ_r.

## 8. Coefficient and moment dynamics

### 8.1 Growth moments

The projection of the windowed growth rate onto P_n equals a particle moment with the corrected weight Ω(1−μ²)r_σ (checked to 10⁻⁶ against direct quadrature):

$$
\mathcal{G}_{\sigma hn}=\int P_n\,\Gamma_W\,dx=\sigma\,\frac{\mu_0\pi^2V_A}{B^2}\sum_{\mathrm{species}}\int dp\,p^3\!\!\int_{\mathrm{channel}\ h}\!\!d\mu\;\Omega\,(1-\mu^2)\,r_\sigma\;\mathcal{X}_\sigma f\;P_n\big(\xi_r\big)\,W(p,\mu,x_r)
$$

For a representation by coefficients, Δa_n = Δt(2n+1)𝒢_n/Δx is exact for frozen f; for the extensive state the matching action source is given in Section 10.

### 8.2 Growth from pitch-angle moments: Monte Carlo and Parker solvers

With f = Σ_ℓ f_ℓ(p)P_ℓ(μ) on (−1, 1) the moments are a contraction with per-cell kernels (checked):

$$
\mathcal{G}_{\sigma hn}=\sum_{\ell=0}^{L}\int dp\Big[K^{\mu}_{n\ell}f_\ell+K^{p}_{n\ell}\frac{\partial f_\ell}{\partial p}\Big],\quad
K^{\mu}_{n\ell}=\sigma\frac{\mu_0\pi^2V_A}{B^2}p^3\!\!\int\!\Omega(1-\mu^2)r_\sigma P_nW\,P_\ell'\,d\mu,\quad
K^{p}_{n\ell}=\sigma\frac{\mu_0\pi^2V_A}{B^2}p^3\!\!\int\!\Omega(1-\mu^2)r_\sigma a_\sigma P_nW\,P_\ell\,d\mu
$$

- Kernels are computed per cell by Gauss quadrature split at μ = V_σ/v and window edges, or cached on a documented grid; scalar B, charge and rate prefactors are factored out of reusable angular maps. Full per-cell kernel tables are not stored (Section 14).
- Monte Carlo accumulation: per cell and momentum shell, f_ℓ(p_j) = (2ℓ+1)/2 × Σ_iw_iP_ℓ(μ_i) / (2πV_cell∫_shell p²dp) with w_i the physical particle number; ∂_pf_ℓ by differences between shells; N_eff = (Σw)²/Σw² and the weighted variance per cell, shell and channel are reported; near Γ = 0 an absolute uncertainty and a finite correlation-time bound are used. In a test (40 shells, L = 3, 2×10⁶ samples) 𝒢₀ of a streaming power law was recovered to 4%, 𝒢₁ and 𝒢₃ to 1.5%, with error ∝ N⁻¹ᐟ².
- Parker solvers insert the derived anisotropy of Section 3.3; the two helicity channels differ through the odd-in-μ parts of a_σ − Ā and of the resonance map and are kept separate.
- The per-particle estimator obtained by integrating by parts onto P_n(ξ_r) has a factor 1/(μ − V_σ/v) and heavy-tailed variance; it is not used.
- A species in the growth sum feeds back on the waves; a species omitted from it is a declared external budget term.

### 8.3 Spectral drift

In coefficient form the drift is −ẋ_σΣ_mK_nma_m − |ẋ_σ|P_n(ξ_in)g(ξ_in) with K_nm = 2 for m > n, m + n odd (checked), and the upwind penalty at the inflow edge. The extensive-moment form is given in Section 10.

### 8.4 Nonlinear transfer and background relaxation

**C1 — relaxation toward the maintained background.** Parallel waves are sheared by a declared background of perpendicular eddies and leave slab resonance (Shebalin et al. 1983; Goldreich & Sridhar 1995; Yan & Lazarian 2002; Farmer & Goldreich 2004). With 𝒩 = −Γ_TI on the total spectrum and policy B, E₀ = 𝓛_wI₀ + Γ_TI₀ − S and the ratio equation contains (checked: the example 𝓛_wI₀/I₀ = 0.2, Γ_T = 0.7, S/I₀ = 0.1 s⁻¹ gives τ = 0.9 s⁻¹)

$$
\mathcal{K}^{(\mathrm{C1})}=-\tau_{C1}\left(1-e^{-g}\right),\qquad \tau_{C1}=\frac{S+E_0}{I_0}=\Gamma_T+\frac{\mathcal{L}_wI_0}{I_0},\qquad \Gamma_T\simeq C_T\left(\frac{\varepsilon\,\kappa}{V_A}\right)^{1/2}
$$

The physical excess damping is −Γ_TδI (Section 5). Γ_T needs the outer scale, cascade rate per unit mass ε, amplitude and counter-propagating population of the declared regime (Lazarian 2016); C_T absorbs the amplitude-versus-intensity convention and is a calibration constant to be fixed against that regime, and the applicability of C1 is tested before event use (gate P5). Pointwise, y = e^g obeys dy/dt = Γy − τ(y − 1); the exact step for frozen rates is y_new = e^{(Γ−τ)Δt}y + τ(e^{(Γ−τ)Δt} − 1)/(Γ − τ), evaluated with expm1 and as y + τΔt at Γ = τ.

**C2 — disabled in the baseline.** For transverse incompressible fluctuations depending only on the parallel coordinate, the Elsässer nonlinear term (z_{−σ}·∇)z_σ vanishes; counter-propagating purely parallel waves do not cascade along k. The transformed local-transfer law of the previous revision is algebraically correct but is phenomenology that implies perpendicular structure (Chandran & Hollweg 2009). It may be enabled only as extension E03 with a stated mechanism, scales, energy recipient, calibration and validation, and never together with C1 as two descriptions of the same transfer.

## 9. Energy and action budgets of the open system

Let e = κV_A be the intrinsic frequency. Multiplying the action equation by e gives the energy equation with the WKB work term:

$$
\frac{\partial(e\mathcal{A})}{\partial t}+\nabla\!\cdot\!(\mathbf v_\sigma e\mathcal{A})+\frac{\partial}{\partial x}(\dot x_\sigma e\mathcal{A})=e\,S_{\mathcal{A}}+\mathcal{A}\,D_{\mathrm{ray}}e,\qquad D_{\mathrm{ray}}=\partial_t+\mathbf v_\sigma\!\cdot\!\nabla+\dot x_\sigma\partial_x
$$

Budgets are computed for the total and for the reference with the same conservative operator and then subtracted; action fluxes are used for action and e-weighted fluxes for energy. Per channel the accounts are: particle exchange with the evolved spectrum (IΓ_W/μ₀, equal to the weighted particle energy increments of the accepted interactions); particle exchange with the reference reservoir (I₀Γ_{1−W}/μ₀) and with any fallback tensor; C1 transfer (−Γ_TδI, signed); imposed forcing S + E₀ on the total and reference accounts; WKB work 𝒜D_ray e; physical boundary fluxes q_f e𝒜 and spectral endpoint fluxes ẋ_σe𝒜; and numerical defects (recovery projection, remap, fallback, floors) reported separately. The residual is normalized by the sum of absolute terms plus a declared floor. Expansion and spectral drift change wave energy while conserving action, so boundary escape and transfer are not the only sinks.

## 10. Conservative discretization on field-line segments and AMR meshes

The primary numerical state is the total extensive wave action of each channel on the physical control volume; coefficient advection is retained only as a reference and fallback. A first-order upwind update of g at Courant 0.1 on two periodic cells with I = 1 and 100 loses 36% of the action in one step (checked); monitoring that defect does not make the update conservative.

Within the parallel-ray reduction the conservative equation is

$$
\frac{\partial\mathcal{A}}{\partial t}+\nabla\!\cdot\!(\mathbf v_\sigma\mathcal{A})+\frac{\partial}{\partial x}(\dot x_\sigma\mathcal{A})=S_{\mathcal{A}},\qquad \mathbf v_\sigma=\mathbf u+\sigma V_A\,\mathbf t_{\mathrm{out}}
$$

integrated over moving control volumes with the relative face speed q_f = (v_σ − w_face)·n̂_f. For moment n the interior drift coupling, with F_x,n the outward spectral-endpoint flux, is (checked: a packet moving toward increasing x increases M₁)

$$
\frac{dM_n}{dt}=-F_{x,n}+\frac{2\dot x_\sigma}{\Delta x}\sum_{m<n,\ m+n\ {\mathrm{odd}}}(2m+1)M_m+S_n-\sum_f\Phi_{fn}
$$

The particle action source per unit arc length of a segment, from the energy exchange 2πp²v a_σD_σ𝒳_σf divided by the intrinsic frequency κ_rV_A (checked):

$$
S_{n,\mathrm{line}}=\sigma A\sum_{\mathrm{species}}\int 2\pi p^3\,dp\int_{\mathrm{channel}\ h}\frac{D_\sigma\,\mathcal{X}_\sigma f}{\kappa_r\,r_\sigma}\,P_n(\xi_r)\,d\mu
$$

where D_σ is the evolved-spectrum part including W; the reference-reservoir and fallback parts are booked separately; for 3D cells A dz is replaced by the physical volume measure.

### 10.1 Face fluxes and the shared reconstruction

One oriented, time-integrated flux per physical face or fine subface, Φ_fn = ∫dt∫_f dS q_f∫𝒜_upP_n dx, is applied with opposite signs to the two neighbours. A physical action spectrum is reconstructed and upwinded at the shared face from the recovered g (positivity-preserving, with a conservative first-order fallback that both neighbours use if either falls back), then its moments are integrated. The discrete geometry update and the mesh-velocity flux satisfy the geometric conservation law, which is tested, not assumed. Explicit transport obeys Δt ≤ V_j/Σ_f|q_f|A_f with the velocities used by the flux; the stability of the multidimensional update and the spectral drift is demonstrated separately.

### 10.2 AMPS field-line segments

A segment is the interval between two vertices with stored data, particle association, ownership and MPI exchange. Its control volume is V_j = ∫A dz with A(z) = Φ_line/|B(z)| from the owned magnetic flux and the host's segment-volume hook (Simpson quadrature in the inspected implementation); particle weights and wave moments use the same volume and the same Φ_line normalization, declared with its provenance. A characteristic-only line without a finite tube normalization needs a separate contract. Endpoint fluxes are F_{j+1/2,n} = A(z_{j+1/2}) q_{j+1/2}∫𝒜_upP_n dx with q = (v_σ − w_endpoint)·t_line, which reduces to σV_A only at a truly material endpoint; a line re-traced from successive MHD snapshots is not automatically material and is handled by conservative overlap at the retracing event. Geometry, area, B and background epochs are kept coherent; the frozen-in identity D ln(ΔsB)/Dt = 2D ln B/Dt − D ln ρ/Dt is a consistency check, never a repair.

- **Insertion, removal, split, merge, re-tracing.** Each owned segment carries one wave state keyed by a stable identity. Splitting an interval along one tube keeps Φ_line (serial subdivision); splitting a bundle into separate tubes partitions Φ_line with positive fractions summing to the parent. For a length split or merge the total action over the physical overlap is conserved and each new spectrum is recovered with its own I₀ and V_A.
- **Distributed implicit sweep.** The one-way bidiagonal sweep (a_j^{n+1} = (a_j^n + λa_up^{n+1})/(1 + λ), λ = V_AΔt/Δs_j) is valid for a single ordered chain in the intensive form; across MPI ranks the upstream updated state is communicated in sweep order or the assembled system is solved; cyclic lines and sign-changing relative speeds need their own solver. The conservative moment solver has its own matrix from volumes, endpoint areas and speeds; the convex-combination property of the scalar sweep is not transferred to it without derivation.
- **Moving shock.** Retiring a whole segment is first order in cell size; a conservative partial-segment or cut-volume operation books action at the actual crossing time with the same swept volume for particles. Boundary classification uses normal relative velocities (Section 11).

### 10.3 3D AMR meshes

Owned cells carry all channel moments, the energy state and the reference; ghost copies are support data and never enter global budgets, which are computed on active leaves only. Face fluxes use the physical volume, face area and face quadrature of v_σ, B direction and action; coarse–fine interfaces accumulate one oriented, time-integrated flux over all fine subfaces and all fine substeps, Φ_coarse,n = Σ_f Σ_r Φ_fn,r, in a flux register per family, helicity and moment, and the coarse conserved state is corrected at synchronization before recovery. The native implementation starts with a common wave time step; fine-level subcycling is enabled only after the common-step conservation tests pass. Restriction transfers eight children to a common physical κ basis before summing extensive moments and energy; prolongation distributes the extensive state with positive weights summing to the parent and reconstructs from the children's actual background measures. The reduced drift ẋ_σ = −[b̂·(∇u)·b̂ + σ b̂·∇V_A] − (∂_t + v_σ·∇)ln B is conditional on the parallel-ray restriction; the alignment angle from short vector-ray integrations is diagnosed and the accepted error stated. Transverse numerical spreading on a Cartesian mesh is measured under rotation and refinement and reported separately from physical cross-field transport; where it exceeds the observable's tolerance the field-line representation is used.

### 10.4 Conservative remapping

A length- or volume-weighted average of coefficients averages g and loses action (e^{⟨g⟩} ≤ ⟨e^g⟩); a nodal log-mean conserves only if the new cell's background action equals the sum of the old weights, which an independently reconstructed I₀ and V_A do not guarantee (two cells with I₀ = 1, 3 and V_A = 1, 2 merged into a parent with I₀ = 2, V_A = 1.5 gain 6.67% with the draft denominator; checked). The remap used here sums the extensive moments of the old cells (or distributes them on a split), then recovers g against the new cell's own C₀(x) = ∫_cell I₀/(μ₀V_A) dV; the recovery projection error is reported. Zero padding is exact only for raising the degree on the same interval; lowering the degree, changing elements or transferring between different spectral partitions uses conservative restriction, common interface quadrature and a realizability check.

## 11. Boundary and initial conditions

Data are required on a boundary for each characteristic entering the domain, judged by q_σ,b = (v_σ − w_b)·n̂ < 0 with n̂ the outward normal and w_b the boundary velocity. Fluxes use (u + σV_At_out − w)·n̂ directly and never divide by a vanishing b̂·n̂ at a tangential intersection.

| Boundary | Relative direction | Treatment |
| --- | --- | --- |
| Fast shock moving outward: u_n,upstream > c_f(θ_Bn) in the shock frame | Both families leave into the shock | No upstream data from downstream; the shock's own injection, acceleration and wave jump or transmission law is a separately specified closure |
| Fixed outer edge beyond the Alfvén point | u ± V_A > 0: both leave | None |
| Material outer edge moving with u | +V_A leaves; −V_A enters | Incoming g₋ or action from the exterior (0 if no SEPs beyond) |
| Material inner edge moving with u | +V_A enters; −V_A leaves | Incoming g₊ (0 unless sources below) |
| Fixed inner edge | By sign of u ± V_A | Data for each incoming characteristic |
| x edges, for drift | By sign of ẋ | g = 0 (zero action) imposed at inflow; outflow flux booked |

Initial condition g = 0. For a shock-fitted coordinate with fixed z_out, L = z_out − z_sh and ζ = (z − z_sh)/L, dζ/dt = [c_σ − (1 − ζ)V_sh,z]/L. Cut cells, partial segments and retirement converge with cell size (test V17). Both polarity choices are tested.

## 12. Finite-step coupling and time integration

The finite-step wave gain must equal the actual weighted particle energy loss of the accepted interactions; a frozen Γ from a mid-step f advanced analytically does not guarantee this. Each accepted interaction produces one record (species, weight, physical time, owner, family, helicity, resonant mode or broadening distribution, window or fallback weight, ΔE and Δp∥ in the plasma frame), and the matched wave energy and action are deposited once: ΔE_wave = −ΔE_particle,scatt and ΔM_n = ΔE_wave P_n(ξ_r)/(κ_rV_A) for an interaction allocated to the evolved spectrum, with the same physically derived weights for reservoir or broadened allocations. If a step would remove more wave energy than the positive spectrum can supply, the coupled step is reduced or solved implicitly; clipping one side breaks conservation. A moment-based mean-field growth estimator may be used with a time-centred coupled solve and an energy-constrained source recovery, and its finite-step transfer is compared with the particle increments.

**Stochastic operator.** The Itô drifts for one family, from the conservative operator in the measure p²dp dμ (checked against the operator), are

$$
b_\mu=\frac{\partial D}{\partial\mu}+\frac{1}{p^2}\frac{\partial}{\partial p}\big(p^2aD\big),\qquad b_p=\frac{\partial(aD)}{\partial\mu}+\frac{1}{p^2}\frac{\partial}{\partial p}\big(p^2a^2D\big),\qquad \text{noise}=\sqrt{2D}\,(1,\,a)
$$

with an independent Wiener process per family; an isotropization event at rate v/λ does not reproduce the μ-dependent tensor unless generator equivalence is demonstrated.

**Composition and order.** Growth and maintenance are one exact nodal step with τ_C1 (Section 8.4). A symmetric composition (particles ½ · transport ½ · drift ½ · [growth + maintenance] · drift ½ · transport ½ · particles ½) is second order only when every subsolver is; backward Euler transport inside it leaves the scheme first order. The achieved order of the complete coupled method is measured from three refinements (gate C2, test V13). Step control: max|Δg| ≤ 0.1; the explicit CFL where used; |ẋ|Δt ≤ Δx/(2N²) for explicit drift; sources accumulated over the common wave synchronization interval with local particle times and physical weights, including particles that migrate between owners.

## 13. Runtime diagnostics

| Diagnostic | What is computed | Action |
| --- | --- | --- |
| Quasilinear amplitude | Resonant δB²/B² summed over four channels on the occupied support, per log band; \|Γ\|/(κV_A) | Flag; predictions in the region are extrapolations |
| Accepted-interaction map | κ_r d_i at every contributing (species, p, μ, cell, family); fractions of exchange and of κ∥ outside tolerance | Use the declared fallback or block predictions there |
| WKB and reflection | ε_R at the longest occupied wavelength; background time scale vs κV_A | Compare with V09 |
| Budgets | Total, reference and excess action and energy accounts of Section 9; residual normalized by the sum of absolute terms plus floor | Reduce Δt, raise N or Q, or revise closures |
| Geometry | Discrete GCL; D ln(ΔsB)/Dt vs 2D ln B/Dt − D ln ρ/Dt | Host advection or retracing suspect |
| Monte Carlo | N_eff per cell, shell, channel; weighted variance; absolute error near zero growth | More particles or longer window |
| Recovery | Residual, conditioning, realizability of the total spectrum, fallback count | Fallback with shared flux; book floors |
| Reservoirs | Power of E₀, C1 transfer, fallback tensor, omitted species, floors | Report |
| Ray alignment | Perpendicular Hamiltonian k̇ component vs κ D_v b̂ + κ̇ b̂ | State accepted error; leave the slab model when exceeded |
| Nonresonant modes | (U_SEP/U_B)(V_s/v) from the injected current | Outside scope if not small |

## 14. Inputs, storage and references used in the test problem

**Inputs stated for every run:** B normalization, polarity, spiral geometry, A(z) and ∂/∂z along the field; the vector flow, its frame, density model, composition and ρuA consistency; background spectrum normalization, support, break, helicity fraction, imbalance, scaling and initial state; turbulence outer scale, ε, imbalance and regime; shock normal, trajectory, speed history and connectivity; injection amplitude, units, species, spectral and angular shape, start and stop; momentum endpoints and units, angular representation, x interval, window definition, boundary motion and exterior waves.

**Conditional profile of Section 4.** Leblanc et al. (1998) density n = 3.3×10¹¹/r̂² + 4.1×10¹²/r̂⁴ + 8.0×10¹³/r̂⁶ m⁻³ with r̂ = r/R☉, R☉ = 6.957×10⁸ m, pure protons; radial mass flux nu_rr² = const with u_r = 4.0×10⁵ m/s at r̂ = 215; B_r = 3.5×10⁻⁹(215/r̂)² T with B_φ/B_r = −Ω_sun r/u_r, Ω_sun = 2.7×10⁻⁶ s⁻¹, equatorial, corotating frame with u∥ = u_r sec ψ so that ρu∥/|B| is constant; field-aligned derivatives ∂/∂z = cos ψ ∂/∂r. A script reproduces the table. The same interpolation and derivative operators are used in focusing, drift, E₀ and the benchmark.

The energy-containing scale need not lie below all SEP resonances: at 5 nT, κ_r = 4.9×10⁻⁵ km⁻¹ (4.9×10⁻⁸ m⁻¹) for 0.5 MeV and 1.4×10⁻⁶ km⁻¹ for 500 MeV protons at μ = 1, while the correlation-scale break at 1 AU is near 10⁻⁶ km⁻¹ = 10⁻⁹ m⁻¹; the overlap is checked for each background.

**Storage.** At 5×10⁵ cells, N = 12, L = 3, 40 shells and 64 angular nodes, four coefficient vectors need about 0.21 GB; full per-cell two-component growth kernels would need 67 GB and resonance tables for two families 266 GB. Scalar prefactors are factored out of reusable dimensionless angular maps; caches are keyed by species, V_A/v, support, window and background generation and invalidated after refinement or background updates; interpolation error in D, κ∥ and the signed exchange rate is validated.

**Extensions requiring separate gates:** E01 dispersive and small-pitch-angle physics (full tensor, matched source, re-run P1 P2 P4 C1 C4); E02 reflection or empirical turbulent maintenance; E03 nonlinear cascade (why the strict slab is extended, modes, invariant, calibration); E04 higher order and extra species or charge signs.

## 15. Verification record and references

| Relation | Check | Result |
| --- | --- | --- |
| Scattering rate with r_σ² (Section 3.2) | Galilean transform of the wave-frame magnetostatic rate | (∂μ/∂μ′)²(1−μ′²)/(1−μ²) = r_σ² |
| Wave-frame equilibrium | f(E − V_σpμ), relativistic, V_A = 0.05c | \|𝒳f\| < 2×10⁻⁸ of \|∂_μf\| |
| Energy–momentum relation | Anisotropic relativistic f, helicity-asymmetric spectrum | ΔE = V_σΔp∥ to round-off |
| Matched growth rate (Section 4), SI prefactor | Same f, both signed branches, V_A/v up to 0.1 | Wave gain = particle loss to 2×10⁻⁶ |
| Growth moment with Ω(1−μ²)r_σ (Section 8.1) | n = 0, 1, 3 vs direct quadrature of ∫P_nΓ dx | 10⁻⁶ |
| Action source with 1/r_σ (Section 10) | Particle-side sum vs ∫ΓI/(μ₀V_A) P_n dx | 2×10⁻⁶ |
| Diffusive closure (Section 3.3) | Kolmogorov slab D, two families | J to round-off; D_pp,eff ≥ 0, = 0 for one family |
| Kulsrud–Pearce limit | CG power law, α = 4.5 | Leading order to 7 digits; exact ∂_pf adds 5% at V_A/v ≈ 0.1 |
| Pitch-angle-moment kernels | Exact f_ℓ | 10⁻¹⁰ |
| Monte Carlo moment estimator | 40 shells, L = 3, 2×10⁴–2×10⁶ particles | 𝒢₀ within 4%, 𝒢₁, 𝒢₃ within 1.5%; ∝ N⁻¹ᐟ² |
| Elsässer, WKB, g equation, drift forms, C1, C2 algebra | Symbolic | Zero residual |
| Legendre K_nm | Quadrature | 2×10⁻¹³ |
| Extensive drift sign (Section 10) | Gaussian packet, ẋ > 0 | dM₁/dt = +0.396 numeric vs +0.396 formula |
| τ_C1 (Section 8.4) | 𝓛_wI₀/I₀ = 0.2, Γ_T = 0.7, S/I₀ = 0.1 s⁻¹ | τ = 0.9 s⁻¹ |
| Finite-frequency thresholds (Section 6) | V_A = 8.86×10⁵ m/s, m_pc² = 1.503277616×10⁻¹⁰ J | 0.4957/0.3318 MeV at μ = 1; 186.9/121.2 MeV at μ = 0.05 |
| Profile reconstruction (Section 4) | Field-aligned derivative, u∥ | L_A,∥ = 8.33, 11.82, 73.83, 616.96 R☉ |
| g-upwind non-conservation | Two cells, I = 1 and 100, Courant 0.1 | 35.96% action loss |
| Log-mean remap with mismatched background | I₀ = 1, 3; V_A = 1, 2; parent 2, 1.5 | 6.67% gain with the old denominator; conservative g_new = −0.0645 |
| Sequential growth/relaxation substeps | Δt = 0.1 … 0.0125 | Error ratios → 2: first order |

Withdrawn: the previous revision's energy-identity, growth-rate and threshold checks that used the rate without r_σ² (they verified algebraic matching of an incorrect rate), its reflection table with radial derivatives, and its C2 energy-conservation claim as physics.

Gates adopted from the review (Sections 19–20 there): P1 tensor; P2 window; P3 background; P4 validity; P5 physical omissions; N1 action transport; N2 drift; N3 recovery; N4 remap; N5 3D AMR; N6 segments; C1 coupled box; C2 time clocks; C3 units; C4 stochastic operator; H1 restart; H2 native 3D; H3 native line; H4 physical event. The older V01–V20 labels of Revision 2 map onto these gates and are superseded.

### References

Entries marked † were checked against the publisher or arXiv record; the others are standard citations to be confirmed before publication.

**Quasilinear theory, scattering tensor and focused transport**

- Beeck, J., & Wibberenz, G. 1986, ApJ, 311, 437
- Hasselmann, K., & Wibberenz, G. 1968, Z. Geophys., 34, 353
- Isenberg, P. A. 1997, JGR, 102, 4719
- Jokipii, J. R. 1966, ApJ, 146, 480
- Roelof, E. C. 1969, in Lectures in High-Energy Astrophysics, NASA SP-199, 111
- Ruffolo, D. 1995, ApJ, 442, 861 †
- Schlickeiser, R. 1989, ApJ, 336, 243
- Schlickeiser, R. 2002, Cosmic Ray Astrophysics (Springer)
- Shalchi, A. 2009, Nonlinear Cosmic Ray Diffusion Theories (Springer)
- Skilling, J. 1975, MNRAS, 172, 557
- Thomas, T., & Pfrommer, C. 2019, MNRAS, 485, 2977 † — Eqs. 53–56: matched tensor with finite-frequency factors

**Wave growth by streaming particles and coupled SEP–wave models**

- Afanasiev, A., & Vainio, R. 2013, ApJS, 207, 29
- Afanasiev, A., Battarbee, M., & Vainio, R. 2015, A&A, 584, A81
- Battarbee, M., Laitinen, T., & Vainio, R. 2011, A&A, 535, A34 †
- Bell, A. R. 1978, MNRAS, 182, 147
- Bell, A. R. 2004, MNRAS, 353, 550
- Gordon, B. E., Lee, M. A., Möbius, E., & Trattner, K. J. 1999, JGR, 104, 28263
- Kulsrud, R., & Pearce, W. P. 1969, ApJ, 156, 445
- Kulsrud, R. M., & Cesarsky, C. J. 1971, ApL, 8, 189
- Lee, M. A. 1983, JGR, 88, 6109
- Li, G., Zank, G. P., & Rice, W. K. M. 2003, JGR, 108, 1082
- Ng, C. K., & Reames, D. V. 1994, ApJ, 424, 1032
- Ng, C. K., Reames, D. V., & Tylka, A. J. 2003, ApJ, 591, 461 † — Appendix A resonance kernel; Appendix B tensor and energy exchange
- Reames, D. V., & Ng, C. K. 1998, ApJ, 504, 1002
- Reames, D. V. 2013, Space Sci. Rev., 175, 53
- Rice, W. K. M., Zank, G. P., & Li, G. 2003, JGR, 108, 1369 †
- Vainio, R. 2003, A&A, 406, 735 †
- Vainio, R., & Schlickeiser, R. 1998, A&A, 331, 793; 1999, A&A, 343, 303
- Wentzel, D. G. 1974, ARA&A, 12, 71
- Zank, G. P., Rice, W. K. M., & Wu, C. C. 2000, JGR, 105, 25079

**Wave transport, WKB, reflection and turbulence**

- Bieber, J. W., Wanner, W., & Matthaeus, W. H. 1996, JGR, 101, 2511
- Bruno, R., & Carbone, V. 2013, Living Rev. Solar Phys., 10, 2
- Chandran, B. D. G., & Hollweg, J. V. 2009, ApJ, 707, 1659 † — Eq. 8: Elsässer nonlinear term; perpendicular-scale dissipation
- Dewar, R. L. 1970, Phys. Fluids, 13, 2710
- Dobrowolny, M., Mangeney, A., & Veltri, P. 1980, PRL, 45, 144
- Farmer, A. J., & Goldreich, P. 2004, ApJ, 604, 671 †
- Goldreich, P., & Sridhar, S. 1995, ApJ, 438, 763
- Heinemann, M., & Olbert, S. 1980, JGR, 85, 1311 †
- Jacques, S. A. 1977, ApJ, 215, 942 †
- Lazarian, A. 2016, ApJ, 833, 131 †
- Shebalin, J. V., Matthaeus, W. H., & Montgomery, D. 1983, J. Plasma Phys., 29, 525
- Tu, C.-Y., & Marsch, E. 1995, Space Sci. Rev., 73, 1
- Velli, M. 1993, A&A, 270, 304
- Yan, H., & Lazarian, A. 2002, PRL, 89, 281102
- Zank, G. P., Matthaeus, W. H., & Smith, C. W. 1996, JGR, 101, 17093
- Zhou, Y., & Matthaeus, W. H. 1990, JGR, 95, 14881

**Background solar wind and host codes**

- Borovikov, D., Sokolov, I. V., Roussev, I. I., Taktakishvili, A., & Gombosi, T. I. 2018, ApJ, 864, 88 †
- Borovikov, D., Sokolov, I. V., Huang, Z., Roussev, I. I., & Gombosi, T. I. 2019, arXiv:1911.10165 † — Section 3 and Appendix A: diffusive balance and flow-acceleration terms
- Leblanc, Y., Dulk, G. A., & Bougeret, J.-L. 1998, Sol. Phys., 183, 165 †
- Sokolov, I. V., Roussev, I. I., Gombosi, T. I., et al. 2004, ApJ, 616, L171
- SWMFsoftware AMPS public source, commit 50c8789eaaa75a23bc0760ab313d49f8ef4534e6 † — segment storage and exchange, flux-tube geometry, SI field-line data contracts, srcSEP3D documentation
- Tenishev, V., Shou, Y., Borovikov, D., et al. 2021, JGR Space Physics, 126, e2020JA028242
- Tenishev, V., Zhao, L., & Sokolov, I. 2022, arXiv:2209.09346 †
- van der Holst, B., Sokolov, I. V., Meng, X., et al. 2014, ApJ, 782, 81
- Zhao, L., Sokolov, I., Gombosi, T., et al. 2024, Space Weather, 22, e2023SW003729 †

**Numerical methods and constants**

- AMReX documentation, AmrCore and flux registers † — general space-time refluxing principle, not a dependency
- Hesthaven, J. S., & Warburton, T. 2008, Nodal Discontinuous Galerkin Methods (Springer)
- Junk, M. 1998, J. Stat. Phys., 93, 1143
- Levermore, C. D. 1996, J. Stat. Phys., 83, 1021
- NIST CODATA 2022 fundamental constants †

## 16. Codex implementation roadmap for AMPS in SI units

This roadmap replaces the Revision 2 roadmap; the physical corrections R01–R18 are incorporated in Sections 1–15. The accepted target is a shared SI particle–wave model with two native AMPS geometry adapters: individual physical field-line segments in 1D and active leaf control volumes on a 3D AMR mesh. Conservative total wave action is the primary numerical state. The Python verification model and the AMPS C++ implementation must implement the same approved equations. A Fortran or M-FLAMPA port is a separate task and does not satisfy an AMPS acceptance gate.

### Implementation contract

| Item | Required contract |
| --- | --- |
| Physical scope | Begin with positive ions, four family–helicity channels, the corrected finite-frequency nondispersive parallel Alfvén tensor, and the full focused particle equation. State the ordering of every omitted term. Use the diffusive closure only where its ordering is tested. Resolve field polarity, pitch-angle axis and host w± mapping together. |
| Canonical units | Use B in T, q in C, mass in kg, momentum in kg m/s, particle energy in J, time in s, length in m, physical κ in m⁻¹, mass density in kg/m³ and I per physical κ in T² m. Gyrofrequency is qB/(γm); Alfvén speed is B/√(μ₀ρ). Input conversions occur once at a named boundary. |
| Spectrum and state | Energy density per logarithmic κ is κI/μ₀ in J/m³. Action density per logarithmic κ is I/(μ₀V_A) in J s/m³. Extensive action moments are in J s. An independently controlled extensive energy state is in J. Recover g only for evaluating a spectrum; never transport or average g as the conserved state. |
| Particle measure | For gyrotropic f, n = 2π∫p² f dp dμ. If p_hat = p/p_ref and f_hat = p_ref³ f_SI, carry dp = p_ref dp_hat and ∂f_SI/∂p = p_ref⁻⁴ ∂f_hat/∂p_hat. Particle weights and segment or cell volumes must reproduce the same SI number density. |
| Reference policy | Implement a dynamically evolved SEP-free reference first, with the same conservative geometry and numerical operators as the total spectrum. An empirical reference requires a separately gated maintenance contract, including nonzero imposed forcing and admissibility of signed rates. |
| Excluded interactions | A validity window belongs inside the species interaction and growth integrals. Define evolved-wave and reference-reservoir transfers separately. A window does not supply a physical dispersion relation or a missing small-\|μ\| scattering model. Until those are specified, affected physical runs remain blocked; bounded verification cases may proceed within their declared domain. |
| Nonlinear and omitted physics | Keep C2 disabled in the strict parallel incompressible slab baseline. Enabling a cascade requires a stated physical mechanism, additional scales or modes and independent validation. C1 and omitted reflection require regime tests before event use; an arbitrary coefficient or a small reflected energy fraction is not sufficient evidence. |
| Geometry and ownership | 3D cells use physical active volume, face area and owned-leaf state. A finite line has declared magnetic flux in Wb, area Φ_line/\|B\| in m² and segment volume ∫A dz in m³. A characteristic-only line cannot silently acquire a finite tube normalization. Geometry, backgrounds, particle sources and wave state share an epoch contract. |

### Codex execution rules

0. Follow the stop-and-ask protocol of Section 0. Every item of the open decisions register D01–D18 is OPEN until a person records a choice; an OPEN item blocks the stages it names, and a needed input that neither the description nor the register supplies is a question, never an inference.

1. In S00, read repository instructions, pin the branch revision and record verified native build and test commands. The proposed paths and test IDs below require final repository locations before implementation.

2. Correct the description before making it an oracle. Map R01–R18 to equations, domains, units, primary sources and tests. Retain independent integral checks; Python-to-C++ parity alone cannot establish correctness.

3. Make one reviewable implementation change per stage or clearly identified substage. Run the affected regressions and prerequisite gates. A stage is complete only after its evidence report and machine-readable result manifest are committed with the implementation revision.

4. Use PASS, FAIL or BLOCKED for every required test. A skipped test is BLOCKED. On failure, preserve the failing case and diagnose it; do not loosen a tolerance, change the physics, disable a required channel or replace an independent oracle merely to make a gate pass.

5. Preserve one oriented space-time flux per shared interface and one accepted particle-interaction record per transfer. A limiter, failed recovery, regrid or rejected coupled step must retain the conservation ledger and give an explicit status.

6. Keep unresolved physical closures as explicit blockers or separately disabled extensions. Independent stages may continue when their prerequisites pass. No downstream physical acceptance can bypass a failed prerequisite.

7. Record every field in the evidence table, including achieved errors, convergence and confidence intervals. Report numerical error separately from approximation and calibration error.

8. Use small deterministic fixtures for algebra and conservation, controlled ensemble tests for Monte Carlo statistics, and native host tests for ownership, geometry and MPI behavior. All active acceptance gates must run on the actual AMPS branch before release.

### Module responsibilities

Proposed placement is a shared module under src/models/sep_common/sep_wave, thin adapters under srcSEP and srcSEP3D, and independent verification and regression directories following the target branch conventions. Keep geometry types out of the shared physics kernels. Do not duplicate resonance, energy conversion or window logic between the two adapters.

| Module | Responsibility |
| --- | --- |
| Shared SI contract | units, species, signs, frames, spectral measures, background snapshots and epoch validation |
| Shared particle physics | finite-frequency resonance, admissible interaction map, full diffusion tensor, matched growth and reservoir transfer |
| Shared wave numerics | extensive moments, independent energy control, realizable recovery, spectral drift, reference evolution and paired source ledger |
| Native 1D adapter | physical tube geometry, endpoint fluxes, segment identities, conservative overlap, line MPI exchange and distributed solve |
| Native 3D adapter | owned and ghost states, common face fluxes, flux registers, active-leaf restriction, AMR transfer, rebalance and restart |
| Verification tools | independent p–μ and p–κ integrals, direct κ-grid reference, deterministic focused transport reference and regression manifests |

### Acceptance tolerances

These criteria apply to the controlled verification fixtures. Production tolerances are fixed against the selected physical observables in S00. All dimensioned comparisons include an explicit SI scale and an absolute floor near zero.

| Test class | Criterion |
| --- | --- |
| Analytic algebra and SI scaling | For well-conditioned binary64 fixtures, require scaled error ≤ 10⁻¹². Use an explicitly dimensioned absolute floor near a zero result. Compare the correct tensor, invariants, source normalization and unit transformations, not only stored constants. |
| Integral exchange fixture | For the archived relativistic unequal-helicity fixture of review Section 21, require evolved energy, reservoir energy and action-source discrepancies ≤ 10⁻⁶ at 192 Gauss points. Repeat 96, 192 and at least one finer quadrature and show convergence. New fixtures receive tolerances from their independently converged reference. |
| Discrete conservation | For a closed or fully booked control-volume test, define residual R_Q after subtracting boundary flux, imposed sources and, for energy, physical work. Require \|R_Q\| ≤ 64 γ_N Q_scale, where γ_N = N ε/(1−N ε), ε is machine precision, N is a documented conservative count of accumulated arithmetic operations and N ε < 0.01. Q_scale includes initial absolute content and absolute exchanged contributions; use a declared SI floor for an empty state. Measure recovery or truncation error separately. |
| Recovery and kernel parity | Use a 10⁻¹¹ scaled residual target for realizable, well-conditioned moment-and-energy recovery fixtures and 10⁻¹² for regular Python-to-C++ kernel fixtures. Report conditioning, quadrature refinement and absolute floors. An infeasible spectrum must trigger the documented conservative fallback or rollback, rather than a looser acceptance threshold. |
| Space and time accuracy | Use at least three refinements with ratio 2 in a smooth asymptotic case and p_obs = log₂(error_h/error_h/2). A first-order design must give 0.85 ≤ p_obs ≤ 1.15 in that case. A claimed second-order design must give 1.85 ≤ p_obs ≤ 2.15. Discontinuities, limiter activation and stochastic noise require a separate stated error measure; they are not used to certify a smooth formal order. |
| Stochastic operator | Predeclare ensemble size, seeds, observables and multiple-comparison policy. Use at least 32 independent batches for the controlled acceptance fixture, 99% confidence intervals with family-wise error controlled at 1%, and at least three particle-count levels to separate bias from sampling error. Booked finite particle–wave energy exchange still satisfies the deterministic conservation bound per realization. |
| Physical event accuracy | Before the release comparison, set numerical error limits for the selected particle spectra, arrival times, anisotropy, wave power and energy exchange in the case manifest. Converge κ, p, μ, space and time separately. This gate is BLOCKED until the required observables and accuracy limits are fixed; an arbitrary N = 12 or a single grid match is not acceptance. |

### Stage dependencies and gate map

Execute a stage only after its listed prerequisites pass. Dependency order controls acceptance; shared components may be developed independently only within already approved contracts. Gate IDs P1–P5, N1–N6, C1–C4 and H1–H4 refer to review Sections 19 and 20. Gate C1 is the coupled-box test and gate C2 is the time-clock test; the turbulent maintenance and cascade closures named C1 and C2 are separate model terms. The revised stage test IDs are distinct from the older V01–V20 verification claims.

| Stage | Prerequisites | Acceptance gates |
| --- | --- | --- |
| S00 Approve the physical specification and host contract | None | R01–R18 disposition and the prerequisites for P1–P5 C3 H2 H3 |
| S01 Implement the corrected SI particle and wave kernels | S00 | P1 P2 C3 and R01–R03 R18 |
| S02 Implement background evolution and physical validity policies | S00 S01 | P3 P4 and the C1 regime part of P5 R04–R07 |
| S03 Build the conservative spectral state and recovery | S01 S02 | N2 N3 and R10–R12 |
| S04 Verify conservative fixed geometry transport and source steps | S03 | N1 N2 N3 and R04 R10 R11 R17 |
| S05 Verify the focused particle equation and deterministic coupling | S01 S02 S04 | P1 P2 C1 and the diffusive comparison part of P5 |
| S06 Implement native AMPS individual field line segments | S00 S03 S04 | N4 N6 H3 and R15 R16 |
| S07 Implement native 3D transport before refinement | S00 S03 S04 | N1 H2 and R13 R14 |
| S08 Add 3D AMR restriction remapping and reflux | S07 | N4 N5 H2 and R12–R14 |
| S09 Implement native Monte Carlo scattering and paired transfer | S05 S06 S08 | C1 C3 C4 and R17 |
| S10 Synchronize local particle clocks and optional wave subcycling | S08 S09 | C2 N5 and R13 R17 |
| S11 Complete restart rebalance and diagnostic persistence | S06 S08 S10 | H1 H2 H3 and R13–R18 |
| S12 Verify moving geometry shocks and approximation domains | S05 S06 S08 S10 S11 | P5 N6 H3 and R07–R09 R13 R16 |
| S13 Accept a converged native physical case and production configuration | S01–S12 for all features used in the release | H4 and all applicable P N C H gates |

### S00 Approve the physical specification and host contract

**Prerequisites:** None. **Acceptance mapping:** R01–R18 disposition and the prerequisites for P1–P5 C3 H2 H3.

**Deliverable:** An approved SI equation set, scope and closure decision register, and a pinned native AMPS integration plan.

#### Implementation work

- Resolve the open decisions register D01–D18 of Section 0 with the owner before any kernel is written; record choices in `docs/sep_wave/decisions.md`. Unresolved items stay OPEN and BLOCK their dependent stages.
- Audit the actual branch, repository instructions, native particle clocks, field-line volume hooks, mesh ownership, MPI exchange and restart paths. Record verified build and run commands and retain a reproducible baseline before adding waves.
- Apply R01–R18 to the model description. Declare positive-ion support, four channels, frame and polarity mapping, reference policy, wave-energy normalization and total action state. Select how independently controlled wave energy is represented.
- Specify the accepted dispersion range, singular-channel behavior and excluded-interaction reservoir or fallback tensor. For the first bounded tests, reject out-of-domain inputs explicitly. Decide which omissions can be justified for the intended event and which require an extension.
- Define shared kernel and geometry adapter interfaces, stable identities, owned state, source clocks, background epochs and output schema. Set case-specific error targets before results are available.

#### Required tests

| Test | Required result |
| --- | --- |
| S00-T00 | The decisions register exists, every D01–D18 entry has an approved choice or OPEN, and every OPEN entry is mapped to the stages it blocks; no stage dependent on an OPEN entry is marked in progress. |
| S00-T01 | Build and run the unmodified native 1D and 3D baseline using the commands recorded from the branch; preserve logs, revisions and configuration. |
| S00-T02 | Dimensional and sign audit of every stored quantity, source and flux; reject an intentionally mismatched unit or frame declaration. |
| S00-T03 | Every required review finding has a correction and a test owner, or an explicit physical blocker. No undocumented inferred closure is accepted. |

**Exit gate:** The corrected specification and integration contracts are complete for the controlled baseline. Unresolved excluded-range, cascade, reflection or shock closures block affected production cases, even when bounded numerical stages can proceed.

### S01 Implement the corrected SI particle and wave kernels

**Prerequisites:** S00. **Acceptance mapping:** P1 P2 C3 and R01–R03 R18.

**Deliverable:** Independent Python reference kernels and stored SI fixtures for resonance, scattering, growth, action source and reservoir transfer.

#### Implementation work

- Implement relativistic velocity, SI gyrofrequency, finite-frequency resonance, channel selection and r_sigma. Include r_sigma² in D_μμ, the matched rank-one D_μp and D_pp, and the corrected SI μ₀π² growth normalization.
- Place the species-, momentum- and pitch-angle-dependent window inside the particle and growth integrals. Split evolved-wave and reference-reservoir contributions; retain the required 1/r_sigma factor in the action source.
- Build independent p–μ particle-energy integrals and p–κ wave-energy integrals with separate mappings and quadratures. Book physical momentum boundaries. Factor cacheable tables only after direct kernels agree.
- Implement the leading diffusive closure as a separately labeled approximation; retain its single-family zero effective momentum diffusion limit and species-dependent finite-frequency resonance.

#### Required tests

| Test | Required result |
| --- | --- |
| S01-T01 | Both families and helicities: tensor symmetry, nonnegative eigenvalues, wave-frame equilibrium and the energy–momentum identity, including finite V_A/v and boundary terms. |
| S01-T02 | Unequal helicities, relativistic particles and a p–μ window ramp: independent evolved plus reservoir energy exchange and action source meet the converged integral criteria. |
| S01-T03 | Physical SI and independently dimensionless fixtures agree to the algebra target for Ω, κ_r, all tensor components, growth, energy and action; test normalized momentum Jacobians. |
| S01-T04 | Recover the unpolarized magnetostatic π/4 limit and the one-family D_pp,eff = 0 limit; verify finite-frequency species dependence and polarity reversal. |

**Exit gate:** P1 P2 and the kernel part of C3 pass with archived inputs and achieved errors. Direct integral disagreement blocks all particle-coupled stages; copying one formula into both integrals is not independent verification.

### S02 Implement background evolution and physical validity policies

**Prerequisites:** S00 S01. **Acceptance mapping:** P3 P4 and the C1 regime part of P5 R04–R07.

**Deliverable:** Epoch-consistent SI background snapshots, a dynamic SEP-free reference and explicit validity and reservoir policies.

#### Implementation work

- Use a single immutable background generation for B, density, flow, gradients, I₀ and windows over each declared bracket. Validate units, frame, polarity and availability before evaluation.
- Implement the dynamic reference with the same transport and quadrature as the total spectrum. Then, if requested, implement empirical maintenance with E₀ = L_w I₀ + Γ_T I₀ − S and τ = (S+E₀)/I₀; do not subtract S twice.
- Test positivity and admissibility for signed maintenance rates. Use relaxation language only for an admissible nonnegative rate. Reconstruct host Alfvén energy with δB² = μ₀ w when w is total Alfvén energy.
- Evaluate validity before a singular resonance division. Record accepted and excluded particle ranges by species and family; do not interpret tapering as proof of nondispersive validity. Validate the declared multi-ion prescription separately.

#### Required tests

| Test | Required result |
| --- | --- |
| S02-T01 | With SEP exchange disabled, a total spectrum identical to its dynamic reference preserves g = 0 under their common discrete update, including nonuniform background and nonzero boundary transport. |
| S02-T02 | Empirical policy: nonzero S, nonzero L_w I₀, positive and negative admissible rates; verify total, reference and excess budgets without double-counting maintenance. |
| S02-T03 | Reproduce the four SI threshold roots in review Section 6 and their family swap under μ reversal; singular-channel and invalid-composition cases produce explicit status. |
| S02-T04 | Reproduce the conditional SI profile with the declared velocity frame and field-aligned derivative; compare C1 applicability with its assumed perpendicular turbulence regime. |

**Exit gate:** P3 P4 pass for the chosen policy. Cases needing an unspecified fallback tensor, dispersive branch or inadmissible forcing remain BLOCKED; no automatic scattering or energy-loss patch is permitted.

### S03 Build the conservative spectral state and recovery

**Prerequisites:** S01 S02. **Acceptance mapping:** N2 N3 and R10–R12.

**Deliverable:** Mandatory extensive total-action moments, independent energy control, spectral quadrature and realizable spectrum recovery.

#### Implementation work

- Implement Legendre polynomials, derivatives and quadrature with explicit physical κ to local x and ξ mappings. Keep the basis interval and background generation in the state metadata.
- Store extensive action moments for all four channels and the selected independent energy state. Treat recovered g and kernel tables as rebuildable caches. Include a spectrum compatible with the chosen moment-and-energy constraints.
- Implement conditioned recovery with positivity or realizability checks, quadrature refinement and a conservative fallback representation. If constraints are infeasible, preserve the extensive state and record the failure; never clip a negative recovered spectrum without booking the effect.
- Derive moment sources from the corrected interaction map and spectral drift from the conservative action equation. For constant positive x drift, the lower-moment coupling has the positive sign given in review Section 10.

#### Required tests

| Test | Required result |
| --- | --- |
| S03-T01 | Legendre orthogonality, recurrence and derivative matrix meet the algebra target on regular fixtures; polynomial projection and quadrature convergence are recorded. |
| S03-T02 | Recover realizable spectra spanning the declared dynamic range; match action moments and energy to the recovery target and report condition number and quadrature error. |
| S03-T03 | Positive x drift increases M₁ for the archived Gaussian packet; reverse drift, nonconstant coefficients and endpoint fluxes are also tested. |
| S03-T04 | Force an ill-conditioned and an infeasible recovery. The fallback or rollback preserves the conservation ledger and exposes status; finite action moments alone are not credited with exact energy conservation. |

**Exit gate:** N2 N3 pass. The transported state and its energy contract are selected before either geometry adapter is written. An optional log-coefficient transport or unbooked recovery projection fails this gate.

### S04 Verify conservative fixed geometry transport and source steps

**Prerequisites:** S03. **Acceptance mapping:** N1 N2 N3 and R04 R10 R11 R17.

**Deliverable:** A direct κ-grid reference and moment-based finite-volume wave solver on fixed physical control volumes.

#### Implementation work

- Use one oriented interface flux and apply it with opposite signs to neighboring extensive states. Assemble spatial transport, conservative spectral drift, physical boundary fluxes and the independent energy-work budget.
- Begin with a validated first-order monotone spatial and time scheme. Determine its actual spatial and spectral CFL limits for the full update; add higher order only as a separately tested extension.
- Implement a stable frozen-coefficient source update using expm1 and the continuous Γ = τ limit. Evolve total and reference states consistently; projection into the conservative representation must satisfy the declared energy constraint.
- Book spectral endpoint escape and WKB energy work separately from particle transfer and maintenance. Keep C2 off. Build a κ-grid reference whose refinement is independent of the Legendre representation.

#### Required tests

| Test | Required result |
| --- | --- |
| S04-T01 | Periodic smooth and high-contrast packets conserve action within the discrete bound. The dimensionless I = 1 and 100 two-cell fixture must avoid the draft's 35.96% loss at Courant 0.1. |
| S04-T02 | Gaussian spectral translation with both signs, physical κ mapping and open endpoints agrees with the independently converged κ-grid solver; all escaped action and energy are booked. |
| S04-T03 | Frozen growth and maintenance reproduce the analytic pointwise solution, including Γ = τ, τ = 0 and a small exponent. Quantify any moment projection error separately. |
| S04-T04 | Three space and time refinements certify the actual first-order design. An implicit update is tested for its particular matrix, boundaries and coefficient signs; positivity is not inferred from an unrelated scalar backward-Euler case. |

**Exit gate:** N1–N3 pass in fixed geometry with separate action and energy budgets. A claim of second-order splitting remains blocked until the complete composed algorithm, including all subsolvers, demonstrates that order.

### S05 Verify the focused particle equation and deterministic coupling

**Prerequisites:** S01 S02 S04. **Acceptance mapping:** P1 P2 C1 and the diffusive comparison part of P5.

**Deliverable:** A deterministic focused-transport reference and a finite-step particle–wave coupling solve with explicit reservoirs.

#### Implementation work

- Implement the approved focused particle equation, including focusing, flow terms, compression, momentum boundaries and the corrected full tensor. Keep the leading diffusive solver separate and document its applicability.
- For each accepted time interval, compute weighted particle energy change and the corresponding evolved-wave or reference-reservoir transfer using the same channel and interaction map.
- Use a coupled or energy-constrained source step that reconciles actual finite particle increments with the conservative wave state. Reject and reduce the interval if the requested transfer cannot be supplied by an admissible spectrum.
- Construct growth, damping and nearly canceling exchange fixtures. Separate interaction energy from external particle injection, escape, compression work and any prescribed background forcing.

#### Required tests

| Test | Required result |
| --- | --- |
| S05-T01 | Pure static focusing, with scattering and flow work disabled, preserves particle number and kinetic energy when boundary fluxes vanish. Homogeneous isotropic compression gives p proportional to ρ^(1/3) in its controlled limit. |
| S05-T02 | Increasing scattering under a controlled ordering recovers the diffusive flux and effective momentum diffusion; outside that ordering the focused solver remains the reference. |
| S05-T03 | Closed box: actual particle energy plus evolved-wave and reservoir energy satisfies the discrete budget per family and helicity for growth, damping, window ramps and near-zero net exchange. |
| S05-T04 | A deliberately excessive damping or energy-transfer request triggers whole-step rollback or a conservative admissible solve; a one-sided particle or wave clip fails. |

**Exit gate:** The deterministic form of C1 passes and its source ledger is ready for native adapters. Comparing Γ computed from a distribution with a nominal wave update is insufficient unless the accepted finite-step energy change also closes.

### S06 Implement native AMPS individual field line segments

**Prerequisites:** S00 S03 S04. **Acceptance mapping:** N4 N6 H3 and R15 R16.

**Deliverable:** A native C++ 1D adapter using actual AMPS segment storage, physical volumes, stable identities and conservative endpoint fluxes.

#### Implementation work

- Bind the shared state to individual segments using the actual owned storage and MPI exchange contracts. Use stable line and segment identities rather than transient pointers or local indices.
- Read finite Φ_line, B and geometry from the SI exchange contract. Use A = Φ_line/|B| and the host's physical segment-volume hook consistently with particle weights. Declare seed-area or flux normalization and its provenance; the legacy π m² registry default is not a measured tube area.
- Compute each endpoint's relative velocity from its actual motion. Use A q times the action integral; q reduces to σV_A only at a material endpoint. Re-traced geometry is not automatically material.
- Implement conservative overlap for segment insertion, deletion, split, merge and re-tracing, including changes of x coordinates and reference weights. Serial subdivision keeps the same flux; dividing a bundle into separate tubes partitions the parent flux.
- Implement distributed ordering for an implicit chain solve. Cyclic lines and sign-changing speeds require their appropriate solve; independent local MPI sweeps are not a global upstream solution.

#### Required tests

| Test | Required result |
| --- | --- |
| S06-T01 | Native nonuniform B, area and segment length: physical volumes agree with the host geometry contract and particle density normalization; uniform physical states and variable-area packets satisfy their correct analytic budgets. |
| S06-T02 | Fixed, material and moving endpoints select incoming data using actual q. Test inward and outward families and polarity reversal. |
| S06-T03 | Pure geometric splitting, merging, insertion and re-tracing conserve action and the declared energy state; evolving backgrounds book physical work separately. Unequal backgrounds use the new physical reference denominator, not copied g. |
| S06-T04 | One distributed line on at least two MPI ranks matches the serial conservative solution; periodic and sign-changing-speed cases exercise the appropriate distributed solver. |

**Exit gate:** N4 N6 and the native segment part of H3 pass in the actual AMPS executable. A continuous line reference or a Python segment harness does not replace tests on individual host segments.

### S07 Implement native 3D transport before refinement

**Prerequisites:** S00 S03 S04. **Acceptance mapping:** N1 H2 and R13 R14.

**Deliverable:** A native AMPS 3D adapter on a single level with owned/ghost state, background epochs and shared physical face fluxes.

#### Implementation work

- Allocate all channel moments, independent energy state and reference state on owned cells. Specify ghost packing, deterministic halo exchange and cache invalidation; never include ghost copies in global physical budgets.
- Use physical cell volumes and face quadrature for transport velocity, B direction and action flux. Share the oriented flux between neighbors. Confirm the adopted parallel-ray reduction remains valid in curved geometry against the vector ray equation.
- Use background fields from one coherent epoch or bracket. Deposit native particle sources only on the owning active cell and bind source records to the wave interval.
- Compare the C++ shared kernels with stored Python fixtures. Run transport on a 3D grid with variable field orientation before adding refinement.

#### Required tests

| Test | Required result |
| --- | --- |
| S07-T01 | Native 3D packets aligned with a grid axis, oblique to it and in a curved field conserve action with booked escape and energy work; compare physical observables under three spatial refinements. |
| S07-T02 | Changing MPI decomposition preserves ownership, shared fluxes and global active-cell budgets within the deterministic bound. |
| S07-T03 | Halo and background epoch mismatch fixtures reject stale states or rebuild caches explicitly; restart-independent kernel fixtures meet C++ parity targets. |
| S07-T04 | Transverse numerical packet spreading decreases under refinement and is reported separately from physical cross-field diffusion; verify the declared ray-alignment approximation. |

**Exit gate:** The single-level native 3D portion of H2 passes. Prescribed scattering-coefficient bridges or a 2D reference harness are not evidence of a self-consistent 3D wave solver.

### S08 Add 3D AMR restriction remapping and reflux

**Prerequisites:** S07. **Acceptance mapping:** N4 N5 H2 and R12–R14.

**Deliverable:** Conservative 3D AMR evolution with common time steps, eight-child transfer, active-leaf budgets and face-flux registers.

#### Implementation work

- Implement coarse-to-fine conservative distribution and eight-child restriction. Transfer child moments to a common physical κ basis before summing extensive moments and energy; recover only afterward. Remap κ support and reference weights when B or basis intervals change.
- At every coarse–fine interface, accumulate one oriented time-integrated action and energy flux from all fine subfaces. Reflux every channel and moment plus the independent energy budget. Start with one common wave time step.
- Use active-leaf masks in sources and diagnostics, excluding covered coarse cells. Preserve data under refinement, derefinement, MPI rebalance and ownership transfer; invalidate geometry-dependent caches.
- Use a common physical κ quadrature or conservative mortar for differing local spectral representations. Degree reduction and fallback retain the flux ledger and quantify the stated projection error.

#### Required tests

| Test | Required result |
| --- | --- |
| S08-T01 | Oblique 3D packet crosses fine–coarse boundaries in both directions and each relevant face orientation. Verify fine-subface sums, coarse corrections and active-leaf budgets. |
| S08-T02 | Eight children with unequal background and V_A restrict to the exact summed extensive state after conservative basis transfer. The archived two-cell 6.67% artificial action gain must be absent. |
| S08-T03 | Repeated refine/derefine and MPI rebalance with high spectral contrast preserve action and controlled energy; test physical κ shifts, changed degree and conservative fallback. |
| S08-T04 | Compare common-step AMR with an independently refined uniform 3D mesh and multiple MPI decompositions. Record convergence of particle-relevant wave power, not only the global total. |

**Exit gate:** N4 N5 pass for common-step native AMR. Subcycling is deliberately deferred to S10; it cannot be enabled before the common-step flux register and conservative remap tests pass.

### S09 Implement native Monte Carlo scattering and paired transfer

**Prerequisites:** S05 S06 S08. **Acceptance mapping:** C1 C3 C4 and R17.

**Deliverable:** A native AMPS stochastic operator reproducing the approved tensor and a shared accepted-interaction energy ledger in both geometries.

#### Implementation work

- For the selected Itô formulation, include drift terms required by the p² dp dμ phase-space measure and all tensor derivatives. Use independent noise for each wave family and preserve the approved boundary behavior.
- Do not substitute an event rate based only on a mean free path for the full μ-dependent operator without demonstrating generator equivalence. Include finite V_A/v, anisotropic and unequal-helicity spectra.
- Write one accepted record containing species, physical weight, accepted time, owner, family, helicity, physical mode, window/fallback classification and particle energy and momentum increments. Book weighted wave or reservoir energy exactly once.
- Deposit the matching extensive action increment using intrinsic frequency κV_A and the selected energy-constrained source representation. Roll back the particle and wave update together if an admissible spectrum cannot supply the exchange.

#### Required tests

| Test | Required result |
| --- | --- |
| S09-T01 | Weak drift and diffusion moments of the native algorithm agree with the tensor generator over predeclared ensembles. Test wave-frame equilibrium, single-family invariant error under time-step refinement, and both families separately and together. |
| S09-T02 | Every realized closed-box interaction satisfies the paired energy budget within the discrete bound; include window ramps, reservoir transfer, growth, damping and particle momentum boundaries. |
| S09-T03 | Ensemble mean growth and pitch-angle statistics agree with the deterministic reference within the predeclared confidence policy. Three particle-count levels and time-step refinement separate sampling error from bias. |
| S09-T04 | Run the same controlled physical scattering case through native segments and native 3D cells with matched physical volumes and weights; compare SI density, energy transfer and tensor statistics. |

**Exit gate:** C1 C3 C4 pass in both native geometries. A confidence interval does not excuse a deterministic energy-ledger defect, and a scalar λ comparison does not satisfy the stochastic tensor gate.

### S10 Synchronize local particle clocks and optional wave subcycling

**Prerequisites:** S08 S09. **Acceptance mapping:** C2 N5 and R13 R17.

**Deliverable:** A documented wave-interval schedule, source synchronization and, if enabled, conservative AMR wave subcycling.

#### Implementation work

- Define the accepted wave interval and accumulate local particle interactions by their accepted times. Carry ownership and partial interval records across cell, segment and MPI migration without duplication.
- Use time brackets for backgrounds and fine ghost states at every substep. Freeze or predict coefficients only at the declared scheme order; reject inconsistent epochs.
- If wave subcycling is required, accumulate coarse flux as the sum of all fine subfaces and all fine substeps over the coarse interval, then reflux and synchronize conservative states.
- Measure the time order of the complete coupled algorithm. First-order backward Euler inside a symmetric composition does not establish second-order accuracy; upgrade each necessary subsolver before making that claim.

#### Required tests

| Test | Required result |
| --- | --- |
| S10-T01 | Particles with unequal local steps cross several owners within one wave interval. Audit interaction IDs and energy ledgers for omitted or repeated transfers, including rejected steps. |
| S10-T02 | Fine/coarse time-step ratio 2, repeated substeps and time-dependent speeds: verify the complete space-time fine flux sum and post-reflux budgets against common-step AMR. |
| S10-T03 | Three synchronization-step refinements for an analytic source problem, a packet and the coupled box report time order for particle moments, wave energy and budget residual. |
| S10-T04 | MPI migration with partial source accumulators and a background bracket change reproduces the serial ledger; test both segment and 3D ownership transitions. |

**Exit gate:** C2 passes. If subcycling is enabled, its N5 tests are mandatory. If it is disabled, mark the optional feature as excluded from the release configuration and retain common-step AMR; do not label its unrun tests PASS.

### S11 Complete restart rebalance and diagnostic persistence

**Prerequisites:** S06 S08 S10. **Acceptance mapping:** H1 H2 H3 and R13–R18.

**Deliverable:** Versioned restart and MPI migration records covering the full conservative state, geometry, reference, clocks and stochastic state.

#### Implementation work

- Serialize owned moments, independent energy and reference state, physical spectral mappings, tube flux or active mesh geometry, epochs, source accumulators, accepted interaction identities and RNG state or deterministic counter keys.
- Rebuild caches after restore using complete dependency keys: species, B, V_A/v, support, window, reference and geometry generations. V_A alone is not a sufficient kernel key unless every other dependence has been factored and fixed.
- Support repartition or rebalance without treating ghost or duplicate wave copies as additional physical content. Reduce Monte Carlo source records once; do not sum replicated physical wave states as sources.
- Persist action and energy ledgers, boundary escape, WKB work, external/reference forcing, projection defects, fallback counts and domain-exclusion diagnostics. Fail clearly on unsupported checkpoint schema or unit metadata.

#### Required tests

| Test | Required result |
| --- | --- |
| S11-T01 | Uninterrupted versus restarted deterministic runs match the conservative state and budgets within the declared bound, including a restart with partial source accumulation. |
| S11-T02 | Restart and rebalance to a different MPI decomposition preserve stable segment or cell ownership, physical tube flux, active masks and wave source identities. |
| S11-T03 | With the same deterministic RNG contract, matched-layout stochastic restart reproduces the realized trajectory. For a changed layout, require exact ledger conservation and the predeclared statistical equivalence test unless bitwise trajectory preservation is part of the contract. |
| S11-T04 | A changed background, basis, species or window invalidates the correct caches; corrupted or incompatible checkpoint metadata produces explicit failure. |

**Exit gate:** H1 passes and both adapters preserve their physical state under ownership changes. A restart containing g alone or omitting partial particle-source clocks is rejected.

### S12 Verify moving geometry shocks and approximation domains

**Prerequisites:** S05 S06 S08 S10 S11. **Acceptance mapping:** P5 N6 H3 and R07–R09 R13 R16.

**Deliverable:** Conservative moving segment and boundary handling with documented shock and ray assumptions and event-specific approximation tests.

#### Implementation work

- Use actual interface or endpoint velocities and physical swept volumes. Apply a discrete geometric conservation law with particle and wave volumes synchronized; use overlap for re-tracing rather than assuming material motion.
- Represent shock-cut partial cells or partial segments consistently for particles, waves and source deposition. Classify a fast shock using normal relative speed and the oblique fast magnetosonic speed; do not use field-line intersection speed as the MHD classification.
- Specify the adopted particle injection, acceleration and wave boundary or jump treatment, including energy sources and frame. An unspecified shock transmission or compression law blocks the physical shock case.
- Test omitted reflection against a non-WKB reference in the event's declared regime, measuring opposite-family scattering and drift effects. Test ray alignment and focused-flow ordering where the background curves or accelerates.

#### Required tests

| Test | Required result |
| --- | --- |
| S12-T01 | A uniform physical state on a deforming control volume satisfies the discrete GCL; moving endpoints and re-traced segments preserve booked action and energy without duplicated spectral drift. |
| S12-T02 | A shock traverses segment interiors and 3D cells with shared swept volume and source timing. Grazing field orientation avoids division by a vanishing field-normal component. |
| S12-T03 | Fixed, material and shock boundaries classify incoming families using relative normal velocity. Test both polarity choices and verify particle and wave escape and injection budgets. |
| S12-T04 | Reflection, C1 applicability, ray alignment and focused/diffusive ordering errors are quantified against the release observables. A failed approximation test requires the corresponding physical extension or a narrower declared domain. |

**Exit gate:** Moving-geometry gates pass, and P5 is satisfied for the intended event. Algorithmic shock transport alone does not validate an unspecified shock-wave physical closure or justify omitting a consequential reflected family.

### S13 Accept a converged native physical case and production configuration

**Prerequisites:** S01–S12 for all features used in the release. **Acceptance mapping:** H4 and all applicable P N C H gates.

**Deliverable:** A reproducible native AMPS acceptance dossier, converged physical comparisons and a measured production configuration.

#### Implementation work

- Converge the independent reference in κ, p, μ, space and time before comparing the reduced spectral representation. Vary polynomial degree, spectral intervals, recovery quadrature and interaction-table resolution separately.
- Compare selected SEP spectra, arrival times, anisotropy, wave power, effective scattering and exchanged energy in native 1D segments and native 3D AMR. Keep injection amplitude, frame, background policy and observer definition fixed during numerical comparisons.
- Measure memory and runtime on the actual branch and hardware. Start from direct or factorized kernels; introduce caches only with validated accuracy and invalidation. Avoid allocating full cell-by-channel-by-species tables without a measured need and memory budget.
- Publish per-gate status, achieved errors and confidence intervals, enabled physics and declared approximation domain. Complete the release only when every required gate passes for that configuration.

#### Required tests

| Test | Required result |
| --- | --- |
| S13-T01 | Independent reference convergence and native 1D/3D comparisons meet the observer-specific limits fixed in S00; all open-system budget terms and numerical defects are reported. |
| S13-T02 | Adaptation, MPI decomposition, restart and the permitted time-step hierarchy do not alter accepted observables beyond their declared numerical error limits. |
| S13-T03 | Cached and direct kernels agree for tensor components, λ and signed growth over the full supported parameter set, including cache invalidation events and near-zero growth. |
| S13-T04 | Production memory and runtime remain within the predeclared budget; all active gate manifests identify the exact code, inputs, configuration and hardware. |

**Exit gate:** H4 and every applicable prerequisite pass. Physical validity is stated for the tested domain and enabled model only. Disabled extensions and blocked physical regimes are excluded explicitly from the release claims.

### Physical extensions requiring separate gates

The baseline does not authorize an unspecified physical closure. Implement an extension when the intended case requires it, then re-run all affected gates. Disabled optional extensions are named in the release configuration; a required extension cannot be counted as completed by disabling it.

| Extension | Physical requirement and revalidation |
| --- | --- |
| E01 Dispersive and small pitch-angle physics | Specify dispersion or broadening, polarization, the full particle tensor and matched wave or reservoir energy source. Re-run P1 P2 P4 C1 C4 and the event comparison. A ν = v/λ substitution or a taper alone cannot pass. |
| E02 Reflection or empirical turbulent maintenance | Add the omitted mode coupling or declared perpendicular turbulence model when the S12 approximation test fails. Preserve channel energy and action bookkeeping, reference forcing and physical work; re-run P3 P5 and native coupling gates. |
| E03 Nonlinear cascade | State why the strict parallel incompressible slab assumption is being extended, the new scales or modes and the transfer operator's invariant. Calibrate against a declared independent physical reference; verify the appropriate energy flux and re-run P5 N2 N3 H4. Keep C2 disabled until this gate passes. |
| E04 Higher order and extra species | Higher order requires a complete scheme convergence gate; extra charge signs require a new helicity and polarity contract. Additional ions require species-dependent resonance and validity checks. Re-run all affected kernels, sources, caches and native acceptance cases. |

### SI test case and production prerequisites

The draft test case may be retained only as an initial computational fixture. Convert its inputs once to SI: the spectral break 10⁻⁶ km⁻¹ is 10⁻⁹ m⁻¹, shock speed 1500 km/s is 1.5×10⁶ m/s, and proton energies 0.5–500 MeV are 8.01088317×10⁻¹⁴–8.01088317×10⁻¹¹ J. The fractional magnetic power 0.02 remains dimensionless. Forty momentum shells, 64 pitch-angle cells and degrees 8 and 12 are starting resolutions, not convergence evidence.

The finite-frequency thresholds show that the lower-energy part of this fixture is not valid across all pitch angles in a nondispersive near-Sun model. State the accepted physical region and the matched excluded-range treatment before using the fixture for physical acceptance. A spectral window cannot make the entire 0.5–500 MeV population physically admissible. The chosen density and Parker geometry also require their declared velocity frame and composition.

Set the injection amplitude, species weights, magnetic flux or 3D cell normalization, shock energy source, observer positions, background policy and permitted physical closures in the parameter file. C_T = 1 is an initial calibration choice, not a derived constant. A streaming plateau test is meaningful only in a declared regime where self-generated scattering controls that observable.

### Evidence required for every gate

| Record | Required contents |
| --- | --- |
| Identification | stage and test ID, PASS/FAIL/BLOCKED, code revision, input and reference hashes, compiler and flags, verified native command and exit code |
| Physical configuration | SI units, species, momentum convention, frames and polarity, active channels, dispersion/window policy, reservoirs, reference policy and enabled omissions or extensions |
| Numerical configuration | mesh levels and active leaves or stable segment geometry, volumes and tube flux, time-step hierarchy, polynomial basis and intervals, quadrature, limiter and recovery policy, MPI layout |
| Results | achieved absolute and scaled errors, residual scale and operation count, observed order and refinement data, stochastic ensemble design and confidence bounds, failed or excluded inputs |
| State continuity | background and geometry epochs, restart schema, source interval and accepted record counts, RNG contract, adaptation or ownership changes |
| Gate decision | required criterion, achieved result, prerequisite status, reviewer disposition and a reproducible failure artifact when a gate does not pass |

Use docs/sep_wave/stage_Sxx.md and a machine-readable stage_Sxx.json as proposed evidence locations, adjusted once to the native repository conventions in S00. An example test label is S08-T01, distinct from the older V01–V20 claims. The stage overview maps the revised work to the review's P, N, C and H acceptance gates. The release definition of done is the corrected SI description, every required gate passed on the actual AMPS branch, reproducible evidence, conservative restart and geometry handling, and numerical errors within the declared physical-observable limits.
