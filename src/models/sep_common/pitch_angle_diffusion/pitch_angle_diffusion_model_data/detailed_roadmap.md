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
