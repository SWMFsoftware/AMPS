# Parallel diffusion coefficient library

This directory contains the dependency-free C++17 implementation of
`PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md`, revision 1.4. It evaluates the
parallel eigenvalue of the symmetric spatial diffusion tensor used by a Parker
transport equation. It does not project that eigenvalue into radial or shock-
normal directions, assemble the full tensor, add drifts, update turbulence, or
infer momentum diffusion.

All public values are SI. Momentum is total particle momentum in kg m s^-1,
charge is signed coulombs, rigidity is positive volts, magnetic field is
tesla, lengths are metres, and diffusion coefficients are m^2 s^-1. Host
adapters must convert external units exactly once before entering this API.

## Delivered boundary

The standalone implementation now spans PD00–PD10. The parser-neutral manager,
model-specific schemas, dispatch pointer, scalar and batch evaluation,
pitch-angle interface, spectral closures, nonlinear closures, adapters, and
table evaluator are present. The supplied revision-1.4 data bundle is retained
under `parallel_diffusion_model_data/`, and the compiled NLGCE-F arrays are
generated without changing their decimal literals or index order.

The approved srcSEP3D portion of PD11 is implemented. Schema 5 calls the same
model-specific parser through `[parallel_diffusion]`, freezes its SHA-256
configuration identity, installs the active function during serial Runtime
configuration, and evaluates it through the existing Parker coefficient path.
The Parker core remains coefficient-agnostic and retains its coherent
field-aligned gradient stencil. The separate srcSEP binding and native one-/
four-rank transport qualification remain pending; standalone or component-test
success must not be reported as MPI qualification.

The srcSEP3D syntax is:

```ini
[run]
schema_version = 5
transport = parker

[transport]
spatial_diffusion_model = parallel-diffusion-library

[parallel_diffusion]
model = constant_kappa
kappa_parallel_m2_per_s = 1.0e18
```

The example value demonstrates syntax, not a calibrated production choice.
All parameter names are case-sensitive and all numeric text is suffix-free SI.
The section is required exactly when the library selector is active. Legacy
schema-4 selectors retain their original meaning.

srcSEP3D currently admits the library models whose runtime inputs it can supply
without reinterpretation: `constant_lambda`, `constant_kappa`, supported
variants of `power_law_lambda` and `broken_rigidity_kappa`, mean-field `bohm`,
`prescribed_lambda_mu_shape`, and tables without an energy-per-nucleon axis.
Selections requiring nucleon count, slab/2D variance, spectral bend-over
lengths, an effective Bohm field, external region/time factors, or an external
perpendicular closure fail explicitly. Total variance is not relabeled slab
variance, and the existing correlation length is not relabeled a bend-over
length.

The paper-specific extensions explicitly excluded by Section 3.1 remain out
of scope: complete SOQLT, complete composite WNLT, arbitrary directional wave
scattering, momentum diffusion, and a non-axisymmetric perpendicular tensor.
The implemented `broadened_slab` model is exactly the stated slab model, not a
claim to implement those broader theories.

## Public API and manager

Include `parallel_diffusion.h` and use `SEP::ParallelDiffusion`.

`BuildConfiguration(model, assignments, &configuration)` is the semantic
parser manager. The application retains ownership of file syntax and source
locations. The manager selects the model-specific reader, rejects duplicate,
unknown, missing, malformed, and inactive-model keys, validates a temporary
typed configuration, and publishes it only on success. Parameter keys are
case-sensitive; stable model IDs and enumerated values are case-insensitive.
Numeric text must be finite and contain no unit suffix.

Numerically integrated model readers accept `relative_tolerance` and
`maximum_refinements`; nonlinear readers also accept `maximum_iterations`.
These controls govern convergence only and never alter a spectrum, insert a
physical cutoff, or replace a failed closure. Keys not consumed by the
selected model are rejected, so parameters cannot silently survive a model
change with a different meaning.

`SetActiveConfiguration` and `ConfigureActiveModel` install a validated
configuration transactionally. They update this required public dispatch
pointer only after validation succeeds:

```cpp
extern ModelFunction ActiveModelFunction;
```

Configure it once during serial startup. Reconfiguration concurrent with
evaluation is outside the contract. Direct `Evaluate` calls are reentrant;
`EvaluateBatch` evaluates equal-length particle/state arrays with per-point
statuses. The library deliberately uses no result cache, so a revision or
local-state change cannot return a stale coefficient. A future cache must obey
the complete key rules in specification Section 14.4.

Every result carries a SHA-256 configuration fingerprint computed from a
canonical serialization of only the selected model's active schema. Dataset
digests, requested/evaluated model IDs, and provider revisions remain separate
provenance fields.

Typical parser-neutral setup:

```cpp
using namespace SEP::ParallelDiffusion;

Status configured = ConfigureActiveModel(
    "power_law_lambda",
    {{"lambda0_m", "1.495978707e10"},
     {"independent_variable", "rigidity"},
     {"rigidity0_V", "1.0e9"},
     {"independent_exponent", "0.3333333333333333"},
     {"radius0_m", "1.495978707e11"},
     {"radial_exponent", "0.5"}});
```

These numbers illustrate syntax only; they are not installed defaults or a
recommended calibration.

`ParallelResult` returns status separately from optional outputs. On scalar
success, `lambdaParallelM` and `kappaParallelM2PerS` are both present and obey
`kappa=v*lambda/3`. Coupled and fitted nonlinear models also return their
internally consistent perpendicular pair. Missing derivatives remain absent
and set `DerivativeUnavailable`; absence is never interpreted as zero.

`EvaluatePitchAngleDiffusion(mu, ...)` returns D_mu_mu in s^-1 only for
`qlt_slab_spectrum`, `prescribed_lambda_mu_shape`, `broadened_slab`, and an
adapter whose selected underlying closure supplies D_mu_mu. Eigenvalue-only
models return `UnsupportedModel` through this separate interface.

## Runtime state

`ParticleState` requires positive total rest mass and momentum and a finite,
nonzero signed charge. `nucleonCount` is required only when an energy-per-
nucleon variable is selected; it is never inferred from mass or charge.

`LocalState::meanFieldT` is the resolved mean field that defines the local
field-aligned basis. Canonical turbulence variances are total two-component
magnetic variances in T^2. Bend-over lengths are spectral bend-over lengths,
not silently substituted integral correlation lengths. Nonlinear models need
both positive slab and 2D variances and lengths. Pure-component limit formulas
are not invented.

Provider generations are copied into provenance. Optional Jacobians and
gradients must describe the same immutable snapshot. Explicit power laws,
broken laws, Bohm, and NLGCE-F return spatial gradients when every gradient
needed by their selected factors is available. Otherwise the scalar result is
preserved and the gradient is absent.

## Implemented models and parameters

No physical parameter has a production default. Numerical tolerances default
to the initial controls in Sections 10.5 and 14.4 and may be tightened.

### `constant_lambda`

Equation (11). Required key: `lambda_parallel_m` (positive). No background
input is required. The spatial gradient is exactly zero.

### `constant_kappa`

Equation (12). Required key: `kappa_parallel_m2_per_s` (positive). No
background input is required. This is distinct from constant mean free path.

### `power_law_lambda`

Equations (13)–(16). Required keys are `lambda0_m`,
`independent_variable`, its matching reference, and `independent_exponent`.
The variable and reference pairs are:

| Variable | Reference key |
| --- | --- |
| `rigidity` | `rigidity0_V` |
| `total_kinetic_energy` | `kinetic_energy0_J` |
| `energy_per_nucleon` | `energy_per_nucleon0_J` |
| `speed` | `speed0_m_per_s` |

Supplying both `radius0_m` and `radial_exponent` enables the heliocentric
radial factor. Supplying both `field_reference_T` and `field_exponent` enables
the field factor, except that an exact zero exponent is canonical unity.
`use_time_factor` and `use_region_factor` explicitly require their positive
runtime factors when true.

### `broken_rigidity_kappa`

Equations (17)–(19). Required keys are `K_star_m2_per_s`, `rigidity0_V`,
`break_rigidity_V`, `low_slope`, `high_slope`, and `smoothness`. Optional
paired field keys are `field_reference_T` and `field_exponent`.
`use_radial_factor` and `use_region_factor` require already evaluated positive
runtime factors. `K_star` is the speed-factored normalization, not the actual
coefficient at the reference rigidity. Use `KStarFromReferenceKappa` for the
explicit species-dependent conversion.

### `bohm`

Equation (61). Required keys are `eta_B` and `field_definition`, which is
`mean_field` or `effective_field`. This comparison model never clamps another
backend.

### `prescribed_lambda_mu_shape`

Equations (20)–(24). Required keys are `amplitude_mode`, `q_mu`, and `h_mu`.
`target_lambda` mode additionally requires `target_lambda_m`; `fixed_amplitude`
mode requires `D0_per_s`. The target mode recomputes D0 from the shape
integral. A divergent unregularized integral reports `InfiniteMeanFreePath`;
no pitch-angle cutoff or scattering floor is inserted.

### `qlt_slab_spectrum`

Equations (25)–(30) and (35), with the endpoint transformation in Equation
(37). Required `spectrum_form` is `smooth_bendover`, `multirange`, or
`supplied_log_log`. The smooth form reads B, slab variance, bend-over length,
and inertial index from the runtime snapshot.

The normalized `multirange` form additionally requires
`energy_range_index` (`q_E>-1`), `dissipation_index` (`s_d>1`), and
`dissipation_wavenumber_rad_per_m` (positive `k_d`). The runtime snapshot
still supplies the slab variance, bend-over length `ell_s`, and inertial index
`s`; evaluation requires `k_d*ell_s>1`. The normalization uses the exact
three-range integral in Equation (35), including `ln(k_d*ell_s)` when `s=1`.
For strict magnetostatic QLT, `s_d>=2` returns `InfiniteMeanFreePath` because
the actual high-wavenumber inverse-scattering integral diverges. The library
does not hide that divergence behind a finite integration cutoff.

A supplied spectrum additionally requires comma-separated
`spectrum_k_rad_per_m`, `spectrum_power_T2_m`, `spectrum_source_identity`, and
`spectrum_declared_variance_T2`, plus explicit `low_k_policy` and
`high_k_policy`. Each policy is `out_of_domain`,
`zero`, or `power_law`; power-law policies require `<policy>_index`.

The supplied spectrum is already canonical one-sided total-transverse power.
Use these named Section 8.7 API conversions before building it:

- `ConvertOneSidedComponentsToCanonical` sums independently supplied
  transverse components;
- `ConvertTwoSidedComponentsToCanonical` explicitly sums both signs of both
  transverse components, without assuming evenness or axisymmetry;
- `ConvertEvenTwoSidedTotalToCanonical` folds an explicitly even signed total
  spectrum;
- `ConvertOneSidedCyclesPerMToCanonical` applies the cycles-to-radians
  coordinate Jacobian;
- `ConvertFrozenFlowFrequencyToCanonical` applies Equation (36) only when a
  positive sampling-velocity projection and nonempty assumption identity are
  supplied; and
- the two `ConvertQinZhang*ComponentToCanonical` functions apply the stated
  factor four to the source's slab and reduced-radial 2D conventions.

These functions preserve real zeros and reject negative or nonfinite power.
They do not rescale an unidentified convention to force agreement with a
variance. Zero resonant power is not confused with missing coverage. Before
evaluation, the log-log segments and explicit tails are integrated and
checked against the declared variance; divergent or inconsistent tails return
`InconsistentSpectrum`.

The following is a parser-neutral assignment example, not an approved
srcSEP3D or srcSEP input-file block. It shows every model parameter for the
multirange option; the four runtime quantities named in the comments must
come from one coherent local provider snapshot.

```text
model = qlt_slab_spectrum
spectrum_form = multirange
energy_range_index = 0.0
dissipation_index = 1.5
dissipation_wavenumber_rad_per_m = 1.0e-8
relative_tolerance = 1.0e-8
maximum_refinements = 18

# Runtime provider, not model-block parameters:
# mean_B_T, slab_variance_T2, slab_bendover_length_m, inertial_index
```

The numbers are labeled mathematical syntax inputs only. They are not a
physical calibration, and no application installs them as defaults. The
approved srcSEP3D schema-5 binding uses the surrounding section shown above
and accepts suffix-free SI text only. The separate srcSEP parser syntax and
any conversion at that future host boundary remain undecided.

### `qlt_slab_inertial`

Equation (31). It has no physical configuration keys; B, slab variance,
bend-over length, and index are runtime state. Its identity remains separate
from exact Equation (30), and the library does not claim that a caller's r-star
is in the inertial approximation's validity range.

### `broadened_slab`

Equations (39)–(42). It accepts the same spectrum schema as full QLT plus
`kernel` and `width0_per_s`. Kernel values are `lorentzian_constant`,
`lorentzian_linear`, and `gaussian`; the linear Lorentzian also requires
`decorrelation_speed_m_per_s`. An exact zero width selects the exact QLT
branch while retaining the requested broadened-slab model identity. Positive widths use nested,
refinement-checked wavenumber and pitch-angle quadrature with both resonances.

### `nlpa_given_perp`, `nlgc_e`, and `nlgce_n`

Equations (45)–(51). Their physical turbulence inputs are runtime state;
optional parser keys only tighten numerical controls. `nlpa_given_perp`
additionally requires a positive runtime perpendicular coefficient plus its
model identity and revision. All solves use positive logarithmic unknowns,
dimensionless spectral integrals, simultaneous coupled residuals, and the
maximum absolute logarithmic residual gate of 1e-8. Failed solves never return
the last iterate as success.

### `nlgce_f_2014`

Equations (52)–(56). It uses the fixed
`Qin_Zhang_2014_Tables_3_4` coefficient set, natural logarithms, published
input box, and both 288-value arrays. Optional `coefficient_set` may only name
that exact set. Both eigenvalues and the analytic rigidity derivative are
returned. In-box success carries `SurrogateErrorUnbounded` because membership
in the published box is not a certified local fit-error bound. Provenance
contains both audited CSV SHA-256 values.

### `turbulence_adapter`

Equations (57) and (57b). Required keys are `energy_convention`,
`residual_energy_convention`, `moment_source`, `underlying_closure`, and
`vacuum_permeability_H_per_m`. The permeability is explicit because the
post-2019 SI value is measured and the specification does not select a CODATA
release; the library does not silently install the former exact SI value.
`moment_source=input_parameters` additionally requires
`provider_moment_m2_per_s2`, `residual_energy`, and `slab_fraction` in the
model block; `moment_source=local_state` requires those three values from each
immutable runtime snapshot, so transported provider updates reach the closure.
Named energy conventions are
`kinetic_plus_magnetic_variance`, `half_elsasser_sum`, `elsasser_sum`, and
`specific_total_fluctuation_energy`; residual sign is
`kinetic_minus_magnetic` or `magnetic_minus_kinetic`. Density remains a
runtime input. The selected QLT or broadened closure's keys are also required.
The adapter never infers slab fraction, spectral lengths, or directionality.

### `wave_spectrum_adapter`

Revision 1.4 implements only the explicitly admissible balanced,
transverse-axisymmetric, plasma-frame reduction. Required declarations are
`propagation=balanced`, `polarization=transverse_axisymmetric`, `frame=plasma`,
and `underlying_closure`, followed by that QLT/broadened closure's keys.
Directional, propagating, or imbalanced requests fail rather than being
relabeled as magnetostatic scattering.

### `tabulated_parallel`

Required keys are `stored_quantity` (`lambda_parallel` or `kappa_parallel`),
comma-separated `axes`, one `axis_N_values_SI` list per axis,
`coefficient_values_SI`, and `generation_identity`. Supported axes are
`rigidity`, `total_kinetic_energy`, `energy_per_nucleon`, `speed`,
`heliocentric_radius`, `time`, and `mean_field_magnitude`. Values are flattened
row-major with the final axis fastest. Positive axes and coefficients use
multilinear log interpolation. A time axis additionally requires `time_rule`
equal to `linear` (linear interpolation of the positive coefficient after the
other axes) or `step_previous`; time is never logarithmically transformed.
Boundaries are closed, but no extrapolation is performed outside them. Knot
derivatives remain unavailable.

Example parser-neutral table:

```text
model = tabulated_parallel
stored_quantity = lambda_parallel
axes = rigidity,time
axis_0_values_SI = 1e8,1e9,1e10
axis_1_values_SI = 0,86400
time_rule = linear
coefficient_values_SI = <six positive row-major values>
generation_identity = <dataset name and checksum>
```

## Status and diagnostics

Important failures include `MissingInput`, `InvalidBackground`,
`OutsideModelDomain`, `InfiniteMeanFreePath`, `IntegrationFailed`,
`NonlinearSolverFailed`, `InconsistentSpectrum`, and `InvalidConfiguration`.
No failure is converted to zero or a finite fallback. The default policies are
no fallback, no bound, and no extrapolation.

Diagnostics do not turn a successful scalar into a failure.
`QltWeakPerturbationConcern`, `SurrogateErrorUnbounded`, and
`DerivativeUnavailable` identify limitations requiring consumer judgment.
The host decides how to handle ballistic regimes and diagnostics.

## Build and verification

From this directory:

```sh
make -f makefile clean
make -f makefile verify
make -f makefile CXXFLAGS='-O2 -std=c++17 -Wall -Wextra -Wpedantic -Werror' verify
(cd parallel_diffusion_model_data && sha256sum -c SHA256SUMS)
python3 parallel_diffusion_model_data/reference_verification.py
python3 parallel_diffusion_model_data/reference_verification.py --audit
python3 parallel_diffusion_model_data/reference_verification.py --broadened
```

`verify` builds `libparallel_diffusion.a`, runs the original analytical and
parser tests, and runs fixture-driven advanced tests. Expected advanced values
are read from `benchmark_points.json`; they are not copied into the test.
Reports are written to `build/test-report.json` and
`build/advanced-test-report.json`.

The instrumented API selfcheck used for the final standalone audit can be
reproduced without changing the normal build:

```sh
g++ -std=c++17 -O1 -g -Wall -Wextra -Wpedantic -Werror \
  -fsanitize=address,undefined -fno-omit-frame-pointer -pthread \
  parallel_diffusion.cpp parallel_diffusion_advanced.cpp \
  test_parallel_diffusion_advanced_driver.cpp \
  -o /tmp/parallel_diffusion_advanced_sanitized
ASAN_OPTIONS=detect_leaks=0 \
  /tmp/parallel_diffusion_advanced_sanitized selfcheck
```

`reference_verification.py` requires NumPy and SciPy. The default library and
tests require only the C++17 and Python standard libraries. The generated
archive and reports are not source files.

## Thread safety and numerical limits

Direct evaluation uses only caller-owned immutable inputs and is reentrant.
The active global configuration is startup-only. No MPI or PIC types enter the
library. Adaptive quadrature uses float64 and reports failure when its error
gate is not met. The nonlinear solver reports its iteration count and final
maximum logarithmic residual.

Spatial derivatives are returned only where the declared inputs determine
them. Integral-closure spatial derivatives, complete tensor divergence,
field-direction derivatives, table-knot policies, and consumer geometry are
not fabricated. The srcSEP/srcSEP3D adapters must retain their coherent
neighbour stencil wherever a required analytic derivative is absent.

See `IMPLEMENTATION_STATUS.md` for exact evidence and open gates, and
`INTEGRATION_PLAN.md` for the implemented srcSEP3D seam and remaining srcSEP/
native qualification work.
