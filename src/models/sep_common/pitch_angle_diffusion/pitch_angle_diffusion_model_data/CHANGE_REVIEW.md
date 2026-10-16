# Physical and AMPS implementation review of the updated model

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
