# Step 10 change manifest: unified coefficient providers

## Implemented

- Added `util/sep_coefficient_registry.*` with canonical spatial, pitch-angle,
  and mean-free-path names, SI units, parameter schemas, parsers, and one active
  validated configuration.
- Added `coefficient_providers.*` as the only production PIC bridge. Parker,
  `fte-dmumu`, and `fte-mfp` now share source/background validation, provenance,
  snapshot identity, and invalid-value policy.
- Centralized the supported isotropic closures:
  `kappa_parallel=v*lambda_parallel/3` and
  `Dmumu=(v/(2*lambda_parallel))(1-mu^2)`, including inverse conversions.
- Registered analytical QLT, QLT1, Tenishev-2005, Chen-2024, configured Dmumu,
  self-consistent, and SWMF-backed paths through explicit interfaces.
- Rejected the spatial-from-MFP/MFP-from-spatial recursion, direct analytical
  MFP models labelled as self-consistent/SWMF, mutable/mislabeled SWMF state,
  and self-consistent relabelling of imported read-only SWMF state.
- Added CLI selectors for coefficient source, all three provider classes, and
  invalid policy. Startup output records each canonical choice.
- Added `COEF01`–`COEF05`, `COEF-SOURCE`, and
  `make test-coefficients-unit`.
- Step 11 consumes these coefficient contracts only through derived views of
  the single authoritative turbulence state; it adds neither a mover nor a
  competing coefficient owner.

## Compatibility

No-argument defaults preserve the prior srcSEP selection: prescribed source,
spatial diffusion integrated from configured Dmumu, configured pitch-angle
diffusion, Tenishev-2005 mean free path, and fail-on-invalid behavior. Provider
names may be changed with either `--option value` or `--option=value` syntax.
