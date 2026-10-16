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
