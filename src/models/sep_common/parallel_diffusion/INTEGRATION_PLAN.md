# Parker mover integration plan

This plan records how the shared parallel-diffusion selector will enter the
existing Parker paths. The standalone revision-1.4 backends and batch API are
implemented and qualified; D14 and PD11 are explicitly deferred by the user's
2026-10-08 direction. It is a plan, not a claim that either application
currently consumes the new library.

## Shared build and configuration boundary

1. Promote `parallel_diffusion.o` into the canonical `sep_common.a` membership
   only after updating both applications' exact shared-member audits and the
   source manifest in one change. Until then the subdirectory builds a
   standalone archive, preventing an unqualified backend from silently
   replacing the seven established kernels.
2. Add one application adapter per host. Each adapter converts existing host
   species/background records to `ParticleState` and one coherent
   `LocalState`. Unit conversion happens there exactly once. PIC, MPI, mesh,
   and application headers remain outside this library.
3. After D14 is resumed, add a `[parallel_diffusion]` section to the `srcSEP3D` INI parser and a
   `ParallelDiffusion on` block to the legacy `srcSEP` parser. Both collect
   source-located string assignments, convert supported external units if the
   application syntax permits them, then call `BuildConfiguration` and
   `SetActiveConfiguration`. Unknown, duplicate, missing, and inactive-model
   keys fail before particle initialization. Preserve all existing input
   choices as explicit compatibility mappings; do not silently reinterpret an
   old coefficient name as a scientifically different new backend.
4. Freeze the active configuration before mover threads start and publish the
   stable model ID plus configuration fingerprint in startup/output metadata.

## `srcSEP` consumer

The existing Parker core already depends on
`SEP::Transport::SpatialDiffusionProvider`; `AdvanceParker` in
`srcSEP/util/sep_parker_core.cpp` consumes its finite `kappaParallelM2PerS`,
`dKappaParallelDsMPerS`, and provenance. The binding therefore belongs in the
PIC adapter, not in the stochastic step:

1. Add a parallel-diffusion-backed branch to
   `PICSpatialDiffusionProvider::Evaluate` in
   `srcSEP/coefficient_providers.cpp`. Construct total momentum and signed
   charge from the sampled particle/species, and sample position, mean field,
   factors, and background revision from the same field-line state.
2. Evaluate the shared function at the centre and through the existing
   refinement stencil used by `PICSpatialDiffusionProvider`. Until PD09
   provides a complete analytic spatial gradient, retain that coherent
   provider-level stencil; never replace an absent gradient with zero.
3. Translate shared statuses to the existing `StatusCode` and `ValueState`
   without converting unsupported, missing-input, or infinite-MFP states into
   finite coefficients. Preserve the shared provenance string/fingerprint in
   `SpatialDiffusionSample::provenance`.
4. Leave focused `D_mumu` and event-MFP movers on their existing providers
   unless a selected future backend supplies the matching typed output. Do not
   apply spatial diffusion in addition to pitch-angle dynamics for the same
   scattering representation.
5. Extend the parser/registry tests and run the existing PARK01--PARK07,
   coefficient, source-ownership, and controlled scientific suites, followed
   by the roadmap's homogeneous and manufactured varying-coefficient tests.

## `srcSEP3D` consumer

`ResolveLocalTransportImpl` in `srcSEP3D/main_lib.cpp` currently builds a
`CoefficientSelection`, calls `Turbulence::EvaluateLocalScattering`, assigns
`LocalTransportRecord::kappaParallelM2PerS`, and obtains the parallel gradient
from two neighbor evaluations. The shared binding should preserve that
separation:

1. Add a new explicit `SpatialDiffusionModel` selection in
   `srcSEP3D/runtime/run_configuration.h`; do not overload
   `CorrelationMeanFreePath` or `PitchAngleIntegral` with new meanings.
2. In `Turbulence::EvaluateLocalScattering`, dispatch that selection through a
   small adapter that builds `ParticleState`/`LocalState` from the already
   coherent background, turbulence, position, species, charge, momentum, and
   generation arguments. Return the shared finite kappa and provenance.
3. Keep `ResolveLocalTransportImpl`'s bounded neighbor stencil for
   `dKappaParallelDsMPerS` until PD09 supplies and validates every required
   local/provider derivative. The Parker mover in
   `srcSEP3D/transport/parker_transport.cpp` remains coefficient-agnostic.
4. For `constant_ratio` perpendicular diffusion, continue deriving
   `kappa_perp` from the selected shared parallel value exactly once. Other
   perpendicular closures and field-direction derivatives retain their own
   existing ownership.
5. Extend configuration fingerprints, restart compatibility, publication,
   parser tests, COEF3D tests, transport tests, and one-/four-rank native
   evidence. Verify the same shared configuration produces the same scalar
   coefficient in both applications at identical SI states.

## Required transport acceptance before activation

- homogeneous diffusion gives parallel variance `2*kappa_parallel*t`;
- the existing manufactured varying-coefficient cases recover the Itô drift
  using the shared coefficient at every stencil point;
- Parker-spiral tensor projection retains the existing perpendicular choice;
- an oblique-shock check uses the normal projection, not raw parallel kappa;
- scalar results agree across `srcSEP` and `srcSEP3D` for identical particle
  and local-state fixtures; and
- configuration/restart provenance identifies requested/evaluated backend,
  parameters, background revision, and any future explicit fallback.

Application binding is PD11 and remains pending until the standalone roadmap
and D14 host contract are authorized. No current source deck selects this
library. The latest approved direction deliberately leaves this work for a
later task; no srcSEP3D source was changed by the standalone completion.
