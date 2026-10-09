# Perpendicular diffusion implementation status

The standalone revision-2.1 deterministic library is implemented and tested
through the coefficient/statistics portion of roadmap Phases 0-3, the tensor
assembly portion of Phase 4, the explicitly bounded one-dimensional scalar
table portion of Phase 5, the deterministic statistics portion of Phase 6,
and the fully specified parameterized models in Phase 7.

Implemented:

- typed observables, estimators, frames, domain/quality/regime tags, dependency
  ownership, provenance, status, registry, strict model schemas, transactional
  active configuration, dispatch pointer, and batch evaluation;
- smooth spectra, moments and length conversions; slab/2D/composite field-line
  closures; prescribed SEP coefficients;
- NLGC, ENLGC, UNLT, exact/rational implicit-slab, FLPD, corrected RBD/BC,
  stable composite closed form, and shared-owner NLGCE-N/F pairs;
- compound/GCD/slab-running statistics and the conditional pre-diffusive fit;
- calibrated-input Kuhlen algorithm, specified isotropic fits, classical and
  M_A^4 relations, parameterized GCR families, signed drift coefficients, and
  a checked scalar table backend;
- independent compiled fixtures and the supplied source/equation/data audits.

Intentionally not implemented in this library:

- source-gated physics listed in README.md;
- persistent frozen-line trajectory state;
- multidimensional or signed-Hall tables without an application domain,
  interpolation contract, and validation dataset;
- application provider/input adapters, coefficient gradients, full tensor
  divergence, forward/backward mover selection, random increments, timestep
  policy, caching, MPI synchronization, or native application qualification.

No source-gated model falls back to another closure. No unspecified physical
constant, turbulence partition, angular width, fit amplitude, spectral cutoff,
or Kuhlen calibration is installed.

