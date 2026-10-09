# Perpendicular diffusion coefficient library

This directory contains a host-neutral C++17 implementation of
`PERPENDICULAR_DIFFUSION_COEFFICIENT_MODEL.md`, revision 2.1. The library
evaluates perpendicular spatial coefficients, paired parallel/perpendicular
closures, signed Hall coefficients, field-line coefficients, and explicitly
tagged finite-time diagnostics. It does not own a mesh, MPI state, particle
storage, random streams, a turbulence evolution model, or a Parker mover.

All dimensional API values are SI. Momentum is total particle momentum in
kg m s^-1, charge is signed coulombs, magnetic field is tesla, length is
metres, time is seconds, field-line diffusion is metres, and particle spatial
diffusion is m^2 s^-1. Variances in `TurbulenceSample` are total magnetic
variances, not per-Cartesian-component values. The implementation applies the
factors of one half required by the specification; callers must not pre-divide
them.

## Delivered boundary

The deterministic coefficient/statistics library implements every revision
2.1 registry entry whose equations and required inputs are complete. Forty-one
identifiers are executable, including diagnostics that are deliberately not
mover eligible. Five stable identifiers fail with `StatusCode::SourceGate`:

- `frozen_fieldline` is persistent trajectory state, not a deterministic
  coefficient call;
- `qlt_perp` lacks a complete general orbit/resonance closure;
- `unlt_fgr` lacks the harmonics, singular-limit, spectrum, and velocity-
  diffusion contract needed to execute Equation (U11);
- `iso_fit_casse_2002` and `restricted_scattering` lack an absolute
  normalization and complete calibrated domain.

An exact historical Corti preset is also not claimed. `corti_ams02` is the
specified parameterized family and requires the caller's explicit
`parameterized_folded_radian` convention, inverse-radian width, slopes,
normalization, and calibration identity. The Kuhlen backend is executable only
with caller-supplied `gamma_K`, transverse correlation length, calibration
identity, and `first_upward_crossing` root convention; none is inferred from
the partial source table.

Application integration is not part of this standalone change. Neither
`srcSEP` nor `srcSEP3D` calls the library yet. `INTEGRATION_PLAN.md` records the
remaining provider, gradient, tensor-divergence, mover, and MPI qualification
work. A standalone PASS is not Parker-equation or AMPS qualification.

## Public API and manager

Include `perpendicular_diffusion.h` and use
`SEP::PerpendicularDiffusion`.

`BuildConfiguration(model, assignments, &configuration)` is the model parser
manager. The application retains ownership of file syntax and source-location
diagnostics. Each selected model has its own strict schema: duplicate,
unknown, missing, malformed, nonfinite, or inactive-model keys fail closed.
The output is written only after complete validation, so a failed parse cannot
partially replace a working configuration.

`SetActiveConfiguration` and `ConfigureActiveModel` transactionally install a
validated configuration and the public dispatch pointer:

```cpp
extern ModelFunction ActiveModelFunction;
```

Configure it during serial startup before mover threads exist. Direct
`Evaluate` and equal-shape `EvaluateBatch` calls are reentrant; concurrent
reconfiguration is outside the contract. The library has no result cache.

Example parser-neutral configuration:

```cpp
using namespace SEP::PerpendicularDiffusion;

Status status = ConfigureActiveModel(
    "flpd_complete",
    {{"relative_tolerance", "1e-10"},
     {"maximum_refinements", "20"},
     {"maximum_iterations", "120"}});
```

Those numbers are numerical controls, not physical parameters or spectral
cutoffs. The local call must separately supply the particle, mean field,
component turbulence variances/scales/indices, and a coherent parallel
dependency. Model identifiers are case-insensitive; parameter names and
enumerated values are exact strings.

`ModelResult` separates status from optional outputs. A failed model has no
coefficient. Successful anomalous compound transport uses
`DiffusionRegime::NoNormalDiffusion` and retains its finite-time moment; it is
not relabeled as a nonzero asymptotic coefficient. Field-line coefficients,
pitch-angle coefficients, conditional medians, particle moments, symmetric
spatial eigenvalues, and signed Hall coefficients occupy different members.

Every successful result records the requested/evaluated identity, equation
version, configuration fingerprint, sample identity, quality, domain state,
dependency owner, and numerical method where applicable. The current
configuration fingerprint is a canonical FNV-1a-64 identity for restart and
comparison, not a cryptographic asset digest. Supplied tables and data retain
their separately declared checksums.

## Runtime state and invariants

`ParticleState` contains positive rest mass, nonzero signed charge, positive
total momentum, and optional pitch-angle cosine. Relativistic speed and
rigidity are computed centrally without a species default.

`LocalState::meanFieldT` defines the ordered-field direction and magnitude.
Ordered results require a nonzero field. Statistical zero-mean isotropic
results use `FrameKind::Isotropic` and fabricate no axis. Equal transverse
eigenvalues may use any deterministic perpendicular factorization because the
tensor is rotation invariant. Unequal eigenvalues require caller-supplied,
right-handed physical `perpendicularAxis1/2`; the library rejects missing or
nonorthogonal axes.

`ParallelInput` is a value plus model, equation, and sample identity. A
perpendicular-only closure will not accept an unlabeled scalar. If the
parallel and turbulence sample fingerprints differ, evaluation fails with
`InconsistentPair`. NLGCE-N and NLGCE-F do not consume `ParallelInput`; the
existing parallel library is the sole owner of the coupled solve and returns
both coefficients.

The smooth slab/2D provider contract uses total transverse component
variances `slabVarianceT2` and `twoDVarianceT2`, bend-over scales, `s>1`, and
`q>-1`. No slab fraction, 2D fraction, cutoff, integral scale, or outer-scale
conversion is inferred. Isotropic fits consume `totalVarianceT2` and their
own explicitly named scale. A missing optional is unavailable, never physical
zero.

## Implemented models and configuration keys

No physical parameter has a production default. Numerical implicit models
may additionally accept `relative_tolerance`, `maximum_refinements`, and
`maximum_iterations`.

| Model(s) | Equations and required configuration |
| --- | --- |
| `constant_kappa_perp` | P1; `kappa_perp_m2_per_s` |
| `constant_lambda_perp` | P1; `lambda_perp_m` |
| `ratio_kappa` | P1; `eta_kappa`; supplied parallel dependency |
| `power_law_perp` | P2; `kappa0_m2_per_s`, all four reference values and exponents |
| `pitch_angle_perp` | P3; `D0_m2_per_s`; shape `isotropic`, `abs_mu`, or `sqrt_one_minus_mu2` |
| `droge_lambda_perp` | P5-P6; `alpha_D`, `radius0_m`; runtime mu, radius, spiral cosine, parallel dependency |
| `paradise_alpha` | P7; `alpha_P`, `field_reference_T`; supplied parallel dependency |
| `fl_slab`, `fl_2d`, `fl_composite` | F2-F4; smooth spectrum runtime state; F3/F4 require q>1 |
| `flrw_particle` | F5/P4; `a_FL`; mode `local_mu` or `isotropic_average`; identified field-line dependency |
| `nlgc`, `nlgc_slab_kernel_diagnostic` | N1-N4/N8; `a_squared`; production model rejects pure slab, diagnostic requires it |
| `enlgc_2d` | N5; exactly flat q=0 2D spectrum |
| `nlgce_n`, `nlgce_f_2014` | E1-E9; no physical keys; optional numerical controls; pair delegated to the shared parallel backend |
| `unlt` | U1-U2; `a_squared`; transverse slab contributes zero in the stated distributional limit |
| `implicit_slab_exact_2016`, `implicit_slab_rational_2016` | U5-U9; exact positive/asymptotic kernel or separately identified rational kernel |
| `flpd_complete` | D1-D2; smooth 2D spectrum with q>1 and F4 composite field-line closure |
| `rbd_bc` | B1-B4/B9; `a_squared`, `area_energy_index`, `area_inertial_index`; pure-2D S12 implementation |
| `composite_closed_2019` | B5-B8; `perpendicular_length_m`; profile tag `shalchi_2019_bend_over`, `snodin_2022_integral`, or `parameterized_tagged_length` |
| `compound_diffusive_lines` | F8; `age_s`; identified field-line and parallel dependencies |
| `gcd_compound` | C3-C6; `age_s`; q=0 2D state and parallel dependency |
| `enlgc_slab_secant_diagnostic` | N6; `age_s`; slab spectrum and parallel dependency |
| `prediffusive_fit` | C10-C11; age, amplitude, gyroperiod, transition times/exponents, and `calibration_id` |
| `iso_fit_candia_roulet_2004` | I4-I5; `outer_scale_m`, spectrum row, and domain policy |
| `iso_fit_snodin_2016` | I6-I7; sharp-cutoff Kolmogorov pair, `outer_scale_m`, domain policy |
| `iso_kappa0_snodin_2016` | I6 zero-mean variants; outer scale, named spectrum, domain policy |
| `iso_fit_kuhlen_2025` | I8-I13; all fit parameters, negative `gamma_K`, calibrated transverse length, `calibration_id`, and root rule |
| `classical_scattering` | I14; supplied parallel coefficient; returns perpendicular pair and signed Hall value |
| `yan_lazarian_ma4` | I16; `c_M4`, sub-Alfvénic Mach number, injection scale; scaling-quality tag |
| `nwu_ratio_polar` | G3-G4; both ratios, polar enhancement, fold angle/width, `folded_radian` convention |
| `corti_ams02` | G2/G5-G6 parameterized family; normalization, field/break references, slopes, smoothing, angular width/convention, calibration identity |
| `helmod_ratio` | G7; `rho_H`, comma-separated colatitude/factor arrays, polar-function identity, and explicit `linear` interpolation |
| `drift_weak_scattering` | G8; no physical configuration; runtime particle and ordered field |
| `drift_rigidity_reduction` | G9; `K_A0`, `rigidity_A_V` |
| `drift_classical` | I14; supplied parallel coefficient |
| `drift_candia_roulet_2004` | G10; outer scale, named spectrum, explicit Hall-domain bounds/policy/identity |
| `tabulated_perp` | one-dimensional Section 16.3 scalar table; axis, SI axis/value arrays, generation identity, checksum, `reject` boundary policy, and explicit `linear`/`linear_explicit_zero` or `log_log`/`strictly_positive` interpolation policy |
| `unlt_fgr_perturbative_diagnostic` | U10; `maximum_relative_correction`; q>1 pure-2D perturbative diagnostic |

Domain policies for source fits are `reject` or `tag_extrapolated`. An
extrapolated result is explicitly tagged and is not represented as in-domain.
The current table backend is a scalar one-dimensional table. Its interpolation
and zero behavior are mandatory configuration, not hidden defaults. Multidimensional
tables, signed Hall interpolation, derivative-aware interpolation, and cache
ownership require the application domain and validation dataset and are not
invented here.

## Numerical methods

Closure integrals are reduced analytically to dimensionless one-dimensional
forms. Adaptive Simpson quadrature operates on `x=exp(y)` in sixteen pieces;
log-tail bounds expand through 24, 32, 40, and 48. These are convergence
controls, not physical wavenumber cutoffs. Positive implicit ratios use a
model-specific analytic upper bound and bracketed bisection in log space.
Failure to converge returns `IntegrationFailed` or `NonlinearSolverFailed`;
no clipping or fallback coefficient is installed.

The exact implicit-slab kernel uses direct `erfcx` algebra only where stable
and the U7 positive asymptotic series at large argument. The RBD/BC backend
uses `erfc`, not `erfcx`, preserving the stated backtracking correction.
FLPD solves D2 directly and evaluates F4 without an asymptotic branch switch.
All applicable externally supplied zero-parallel branches return an exact
successful zero before constructing Q or a logarithm, after geometry and
required-input validation.

## Example input mapping

The shared manager is intentionally parser-neutral. A future application deck
might map this section:

```text
[perpendicular_diffusion]
model = flpd_complete
relative_tolerance = 1e-10
maximum_refinements = 20
maximum_iterations = 120
```

to:

```cpp
ConfigureActiveModel("flpd_complete",
    {{"relative_tolerance", "1e-10"},
     {"maximum_refinements", "20"},
     {"maximum_iterations", "120"}});
```

This illustrates the library keys only. It is not accepted AMPS syntax until
an application parser explicitly implements and tests that mapping. Runtime
turbulence and parallel dependencies come from coherent local providers, not
from global model parameters.

## Build and verification

From this directory:

```sh
make -f makefile
make -f makefile verify
make -f makefile clean
```

`verify` builds both shared coefficient libraries, runs the compiled C++
contract/equation tests, verifies all 576 NLGCE coefficient literals and
paired fixtures through the supplied scripts, recomputes the independent
spectral/closure fixtures, checks revision-2.1 review additions, compares all
generated specification tables, and audits document structure. The review-
addition script additionally requires Python `mpmath`; `verify` prints an
explicit SKIP when it is unavailable. In a release environment where that
audit is mandatory, install NumPy, SciPy, and mpmath and run:

```sh
make -f makefile verify-high-precision
```

The compiled tests cover strict parsing/source gates, spectral constants and
moment divergence, F4, prescribed coefficients, tensor assembly, U6/U8,
FLPD, ENLGC, RBD/BC, every tested zero branch, B6, anomalous F8 routing,
NLGCE pair ownership, and table boundaries. The Python artifacts explicitly
state that they are mathematical/source-transcription checks, not physical or
AMPS validation.

## Known qualification limits

- There is no srcSEP/srcSEP3D provider or input-parser binding, mover call,
  tensor-divergence implementation, stochastic covariance test, native build,
  MPI run, or event validation in this library change.
- Spatial gradients and implicit sensitivities remain unavailable and are
  reported as such in numerical diagnostics; a host must not treat absence as
  zero. Full divergence additionally needs frame derivatives owned by the
  background/transport layer.
- Persistent frozen-line realization, retracing, random-stream ownership, and
  restart serialization belong to the trajectory layer and remain gated.
- Parameterized fits are equation implementations, not endorsements or new
  calibrations. Physical validation against measurements or particle
  simulations has not been performed here.
