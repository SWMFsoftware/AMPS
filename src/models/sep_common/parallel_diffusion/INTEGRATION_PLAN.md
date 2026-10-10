# Parker mover integration record

This record distinguishes the implemented srcSEP3D binding from the still-open
srcSEP and native-qualification work. The standalone revision-1.4 backends and
batch API remain shared and host-neutral. D13--D16 were approved for srcSEP3D
on 2026-10-08; those decisions do not silently apply to srcSEP.

## Shared build and configuration boundary

1. srcSEP3D builds the canonical standalone `libparallel_diffusion.a` and
   flattens its two objects into `mainlib.a`; the seven-member `sep_common.a`
   remains unchanged until both application ownership audits are updated.
2. Add one application adapter per host. Each adapter converts existing host
   species/background records to `ParticleState` and one coherent
   `LocalState`. Unit conversion happens there exactly once. PIC, MPI, mesh,
   and application headers remain outside this library.
3. srcSEP3D schema 5 now implements `[parallel_diffusion]`. The future
   `ParallelDiffusion on` block in the legacy `srcSEP` parser will collect
   source-located string assignments, convert supported external units if the
   application syntax permits them, then call `BuildConfiguration` and
   `SetActiveConfiguration`. Unknown, duplicate, missing, and inactive-model
   keys fail before particle initialization. Preserve all existing input
   choices as explicit compatibility mappings; do not silently reinterpret an
   old coefficient name as a scientifically different new backend.
4. srcSEP3D freezes the active configuration during serial Runtime setup and
   includes the stable model ID, complete validated parameters, and library
   fingerprint in its physics/restart manifest.

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

`ResolveLocalTransportImpl` in `srcSEP3D/main_lib.cpp` builds a
`CoefficientSelection`, calls `Turbulence::EvaluateLocalScattering`, assigns
`LocalTransportRecord::kappaParallelM2PerS`, and obtains the parallel gradient
from two neighbor evaluations. The shared binding should preserve that
separation:

1. Implemented: a new explicit `SpatialDiffusionModel` selection in
   `srcSEP3D/runtime/run_configuration.h`; do not overload
   `CorrelationMeanFreePath` or `PitchAngleIntegral` with new meanings.
2. Implemented: `Turbulence::EvaluateLocalScattering` dispatches through a
   small adapter that builds `ParticleState`/`LocalState` from the already
   coherent background, turbulence, position, species, charge, momentum, and
   generation arguments. Return the shared finite kappa and provenance.
3. Implemented: `ResolveLocalTransportImpl` keeps its bounded neighbor stencil for
   `dKappaParallelDsMPerS` until PD09 supplies and validates every required
   local/provider derivative. The Parker mover in
   `srcSEP3D/transport/parker_transport.cpp` remains coefficient-agnostic.
4. Implemented: `constant_ratio` perpendicular diffusion continues deriving
   `kappa_perp` from the selected shared parallel value exactly once. Other
   perpendicular closures and field-direction derivatives retain their own
   existing ownership.
5. Implemented component scope: configuration/restart fingerprints plus
   `CFG3D17` and `COEF3D08`. Still required: native one-/four-rank evidence and
   cross-application scalar parity after srcSEP obtains its own approved binding.

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

## Implemented srcSEP seam

1. Input: schema 4 of the srcSEP `--input` INI (`util/sep_initialization.cpp`)
   adds an optional `[parallel_diffusion]` section whose raw body lines and line
   numbers are stored uninterpreted. Schemas 1-3 reject the section.
2. Selection: `--spatial-diffusion-provider parallel-diffusion-library`
   (`SpatialKind::ParallelDiffusionLibrary` in the canonical coefficient
   registry). The registry admits it only for the `parker` mover; fte-dmumu and
   fte-mfp (including `mean-free-path=from-spatial`) are rejected.
3. Installation: `adapters/parallel_diffusion_adapter.cpp` (`Configure`, called
   from `main.cpp` before `Initialization::Install` and before AMPS/MPI
   initialization) requires the section exactly when the provider is selected,
   calls `ParseSection`, applies `CheckHostInputAvailability` with all inputs
   unavailable, installs the model, and copies the model ID and SHA-256
   fingerprint into the startup fingerprint.
4. Evaluation: `PICSpatialDiffusionProvider` (`coefficient_providers.cpp`)
   samples species mass/charge, `|p|` from the Parker speed, the segment's
   Cartesian position (srcSEP's heliocentric frame, Sun at (0,0,0); a field
   line need not start at the origin), the snapshot B vector and generation,
   and the step
   epoch, and calls `ActiveParallelDiffusion` through the adapter. The existing
   refined arc-length stencil supplies `d(kappa)/ds`.
5. Evidence: PIC-free component tests `make -C srcSEP
   test-parallel-diffusion-binding-unit` (PDB01-PDB05) and COEF01. Natively,
   the native srcSEP build (`./Config.pl -application=test/sep_parker_spiral__field_line`,
   `make -j`) compiles and links the PIC sampling wrapper and the `main.cpp`
   wiring; one-rank `--initialization-only` runs install the model and reject a
   missing section or a missing provider as specified. A transport step that
   evaluates the coefficient natively has not run: the shipped example aborts in
   the Parker mover's plasma-density check at step 1 with or without this
   provider, so no native transport or MPI evidence exists for it yet.

PD11 is partial: both application bindings are implemented and
component-tested, while the native acceptance bullets above (including
cross-application scalar parity) remain pending. No maintained production deck
selects srcSEP3D schema 5 or the srcSEP library provider yet; the srcSEP
example `srcSEP/examples/sep_parker_mesh_parallel_diffusion.in` demonstrates
syntax only.
