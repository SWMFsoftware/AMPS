# Parallel diffusion coefficient library

This directory is the shared, application-independent home of parallel
diffusion coefficients intended for the Parker equation in `srcSEP` and
`srcSEP3D`. The scientific contract is
[`PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md`](PARALLEL_DIFFUSION_COEFFICIENT_MODEL.md),
revision 1.3. The implementation uses C++17, SI units, and no PIC, MPI, mesh,
or application headers.

The current code is the data-independent PD01/PD02 release slice. It provides
the common particle/local-state API, validated model registry, parser-facing
configuration bridge, active function pointer, relativistic conversions, and
the five explicit analytical prescriptions from roadmap stage PD02. It is not
yet wired into either application mover. Later spectrum, QLT, nonlinear,
polynomial, table, gradient, batch, and provider stages remain unfinished.

## Physical scope and conventions

The returned quantity is the parallel eigenvalue of the symmetric spatial
diffusion tensor, not a diffusion coefficient in an arbitrary coordinate
direction. Every successful result satisfies

```text
kappa_parallel = v * lambda_parallel / 3.
```

A host that needs radial diffusion, shock-normal diffusion, or a Cartesian
tensor must combine this value with its magnetic-field direction and its
separately selected perpendicular coefficient. Particle drift is an
antisymmetric/advective contribution and remains outside this library. The
library also does not infer momentum diffusion `D_pp`; doing so would require
additional wave-frame, directional-wave, and scattering-center physics.

`ParticleState::momentumKgMPerS` is the positive magnitude of total particle
momentum, not kinetic energy, momentum per nucleon, or one Cartesian
component. `massKg` is total rest mass and `chargeC` is signed charge. The
implemented scalar laws use the magnitude of charge in
`rigidity = p*c/|q|`, but the sign is retained for future closures where wave
direction or polarization can matter. Neutral particles, non-positive mass,
and non-positive momentum are outside this charged-particle API.

Kinematics are exact relativistic quantities computed centrally:

```text
gamma = hypot(m*c, p)/(m*c)
beta  = p/hypot(m*c, p)
v     = c*beta
R     = p*c/|q|
T     = (p*c)^2 / (sqrt((m*c^2)^2 + (p*c)^2) + m*c^2).
```

The last form is algebraically equal to total energy minus rest energy but is
numerically safer for nonrelativistic particles. Energy per nucleon is formed
only when the caller supplies `ParticleState::nucleonCount`; the code never
infers mass number from rest mass or charge state.

`LocalState::positionM` is heliocentric Cartesian position. An adapter for a
translated computational mesh must subtract the configured solar origin
before constructing it. `meanFieldT` is the resolved mean field defining
`B0`; it is not automatically an RMS or total field. The separately named
`effectiveFieldMagnitudeT` is used only when the Bohm configuration explicitly
selects that convention. All scalar inputs and outputs are SI. Host parsers
must convert AU, GV, nT, MeV, or cgs coefficient values exactly once at the
application boundary.

The Parker diffusion approximation presumes a nearly isotropic distribution
with scattering rapid relative to distribution/background evolution. A
positive finite coefficient establishes numerical evaluation only; it does
not establish small Knudsen number, weak focusing, or validity during early
anisotropic SEP arrival and shock-precursor regimes. The present PD02 inputs
do not contain the distribution/background length scales needed to compute
those diagnostics, so the library does not manufacture them.

## Implemented models

No physical parameter below has a library default. Every scale or exponent
used by an implemented model is supplied explicitly by a typed caller or the
parser bridge.

| Stable model ID | Implemented equation | Required configuration |
| --- | --- | --- |
| `constant_lambda` | Equation (11), `kappa=v*lambda/3` | `lambda_parallel_m > 0` |
| `constant_kappa` | Equation (12), `lambda=3*kappa/v` | `kappa_parallel_m2_per_s > 0` |
| `power_law_lambda` | Equations (13)--(15) | `lambda0_m`, `independent_variable`, its SI reference key, and `independent_exponent`; optional paired radial/field parameters and optional supplied time/region factors |
| `broken_rigidity_kappa` | Equations (17)--(19) | `K_star_m2_per_s`, `rigidity0_V`, `break_rigidity_V`, `low_slope`, `high_slope`, `smoothness`; optional paired field parameters and optional supplied radial/region factors |
| `bohm` | Equation (61) | `eta_B > 0` and explicit `field_definition = mean_field` or `effective_field` |

The equations evaluated by these entries are:

```text
constant_lambda:
  lambda = lambda0
  kappa  = v*lambda0/3

constant_kappa:
  kappa  = kappa0
  lambda = 3*kappa0/v

power_law_lambda:
  lambda = lambda0 * (X/X0)^a * (r/r0)^alpha
                   * (B0/Bref)^(-eta) * g_time * g_region
  kappa  = v*lambda/3

broken_rigidity_kappa:
  H(R)   = (R/R0)^a * [((R/R0)^h + (Rb/R0)^h)
                       /(1 + (Rb/R0)^h)]^((b-a)/h)
  kappa  = K_star*beta*(Bref/B0)^eta*H(R)*g_radial*g_region
  lambda = 3*K_star/c*(Bref/B0)^eta*H(R)*g_radial*g_region

bohm:
  r_L    = p/(|q|*B)
  lambda = eta_B*r_L
  kappa  = v*lambda/3.
```

For `power_law_lambda`, `X` is exactly the configured independent variable;
the library does not convert a published energy exponent into a rigidity
exponent. At fixed species, such a conversion changes between the
nonrelativistic and ultrarelativistic limits.

`power_law_lambda.independent_variable` is one of `rigidity`,
`total_kinetic_energy`, `energy_per_nucleon`, or `speed`. Its corresponding
reference key is, respectively, `rigidity0_V`, `kinetic_energy0_J`,
`energy_per_nucleon0_J`, or `speed0_m_per_s`. Energy per nucleon additionally
requires an explicit positive particle nucleon count; charge and mass are
never used to guess it.

The optional power-law radial factor is enabled only by supplying both
`radius0_m` and `radial_exponent`. The optional magnetic factor is enabled only
by supplying both `field_reference_T` and `field_exponent`. `use_time_factor`
and `use_region_factor` are `true`/`false`; when enabled, the local-state
provider must supply a positive dimensionless value. Disabled factors are
exactly unity and create no background requirement. Consistent with Equation
(13), `field_exponent = 0` also makes the field factor exactly unity and does
not require a local magnetic field.

The broken law uses the speed-factored convention in Equation (17):
`kappa(R0)=K_star*beta(R0)` at the reference field and factors, not `K_star`.
`KStarFromReferenceKappa` is the explicit conversion for a caller whose input
is the actual reference coefficient. A radial law is not specified by
Equation (17), so enabling `use_radial_factor` requires the application to
supply the already defined positive `radialFactor`; the library does not
invent one.

For Bohm scaling, `mean_field` uses the magnitude of `LocalState::meanFieldT`.
`effective_field` uses `effectiveFieldMagnitudeT`. The field choice is retained
in the configuration identity. Bohm is a requested comparison model only; it
does not clamp or bound other models.

### Parameter ownership and runtime requirements

The parser bridge uses exact, case-sensitive parameter keys. Model IDs and
enumerated values are case-insensitive. Numeric text must be finite and fully
consumed; unit suffixes are not accepted by this SI-only bridge.

| Configuration key | Model | Units/domain | Runtime `LocalState` requirement |
| --- | --- | --- | --- |
| `lambda_parallel_m` | constant lambda | m, positive | none |
| `kappa_parallel_m2_per_s` | constant kappa | m² s⁻¹, positive | none |
| `lambda0_m` | power law | m, positive | depends on enabled factors below |
| `independent_variable` | power law | `rigidity`, `total_kinetic_energy`, `energy_per_nucleon`, or `speed` | selects exactly one reference key below |
| `rigidity0_V` | power law or broken law | V, positive | none |
| `kinetic_energy0_J` | power law | J, positive | none |
| `energy_per_nucleon0_J` | power law | J per nucleon, positive | positive particle `nucleonCount` |
| `speed0_m_per_s` | power law | m s⁻¹, positive | none |
| `independent_exponent` | power law | dimensionless, finite | none |
| `radius0_m`, `radial_exponent` | power law | m positive; exponent finite | finite positive norm of heliocentric `positionM` |
| `field_reference_T`, `field_exponent` | power or broken law | T positive when active; exponent finite | finite nonzero `meanFieldT` when exponent is nonzero |
| `use_time_factor` | power law | explicit Boolean | positive `timeFactor` when true |
| `use_region_factor` | power or broken law | explicit Boolean | positive `regionFactor` when true |
| `K_star_m2_per_s` | broken law | m² s⁻¹, positive | none beyond selected factors |
| `break_rigidity_V` | broken law | V, positive | none |
| `low_slope`, `high_slope` | broken law | dimensionless, finite | none |
| `smoothness` | broken law | dimensionless, positive | none |
| `use_radial_factor` | broken law | explicit Boolean | positive pre-evaluated `radialFactor` when true |
| `eta_B` | Bohm | dimensionless, positive | selected field below |
| `field_definition` | Bohm | `mean_field` or `effective_field` | matching positive field input |

For the broken-rigidity law, the application owns the physical definition of
`radialFactor` and `regionFactor`; this library only validates and multiplies
them. For the power law, the radial factor is explicitly computed from
heliocentric radius. `timeFactor`, `regionFactor`, and broken-law
`radialFactor` must be positive because the multiplicative implementation is
evaluated in logarithmic form. The provider must document discontinuities and
whether a time variation is simultaneous, convected, or retarded.
Omitted optional Boolean selectors are disabled. Required scale, exponent, and
field-definition keys are never filled from physical defaults.

## Reserved but unavailable models

The registry contains every stable identifier in Section 3 so a misspelled,
deferred, and implemented model remain distinguishable. Selection of the
following identifiers currently returns `UnsupportedModel` rather than a
substitute coefficient:

- `qlt_slab_spectrum`, `qlt_slab_inertial`, and
  `prescribed_lambda_mu_shape` (PD03/PD04);
- `nlgce_f_2014` (PD05);
- `nlpa_given_perp`, `nlgc_e`, and `nlgce_n` (PD06);
- `broadened_slab` (PD07); and
- `turbulence_adapter`, `wave_spectrum_adapter`, and `tabulated_parallel`
  (PD08).

Complete SOQLT, complete composite WNLT, arbitrary directional wave
scattering, momentum diffusion, and non-axisymmetric perpendicular dynamics
are later extensions in the specification and are not registered as
implemented models.

## API and active dispatch

Include `parallel_diffusion.h` and use namespace
`SEP::ParallelDiffusion`. A direct, side-effect-free evaluation is:

```cpp
ModelConfiguration configuration;
configuration.model = ModelId::ConstantLambda;
configuration.constantLambda.lambdaParallelM = 1.495978707e10;

ParticleState particle;
particle.massKg = 1.67262192369e-27;
particle.chargeC = 1.602176634e-19;
particle.momentumKgMPerS = /* total SI momentum supplied by the caller */;

LocalState local;
ParallelResult result = Evaluate(particle, local, configuration);
```

`ParallelResult` distinguishes status from optional values. A successful
explicit model returns both `lambdaParallelM` and `kappaParallelM2PerS`, their
analytic logarithmic rigidity slopes, provenance, and a configuration
fingerprint. Constant models also return an explicit zero spatial gradient.
Spatial gradients of non-constant laws remain absent and carry the
`DerivativeUnavailable` diagnostic until PD09 connects coherent background
gradients. A successful scalar value must not be read as a successful complete
Parker drift.

The public function pointer is:

```cpp
extern ModelFunction ActiveModelFunction;
```

Application parsers must not write it directly. They call
`ConfigureActiveModel(model_id, parameters)` or construct a typed
`ModelConfiguration` and call `SetActiveConfiguration`. Validation is
transactional: configuration and pointer change only after the complete model
schema passes, so rejected input preserves the previous selection. After
initialization, movers call `EvaluateActive(particle, local)`.

The active configuration is a process-global startup choice. Configure it
once during serial initialization, before MPI worker or OpenMP mover activity;
do not reconfigure it concurrently with evaluation. Direct `Evaluate` calls
are reentrant because their configuration is an explicit value.

`BuildConfiguration` is the preferred parser boundary when an application
needs to retain the typed value rather than install global dispatch. It parses
into a temporary candidate and writes the caller's output only after the full
schema validates. `ConfigureActiveModel` adds the installation step. A failed
call to either path does not leave a partially updated typed configuration;
a failed active selection also preserves the prior function pointer.

The public pointer is deliberately visible to satisfy the application
dispatch requirement, but its lifetime is static and ownership remains with
the library. Do not delete it, point it at a function with incompatible
parameter semantics, or mutate it independently of the active configuration.
There is no synchronization around active reconfiguration. Configure exactly
once during serial startup, then treat the active configuration and pointer as
immutable throughout MPI/OpenMP particle advancement.

### Result, derivative, and failure contract

On `Success`, both `lambdaParallelM` and `kappaParallelM2PerS` are present,
finite, positive, and related by the exact computed particle speed. On any
failure, those optionals remain absent. Failure is never encoded as a zero,
NaN, infinity, clamp, or fallback coefficient.

`dLnLambdaDLnRigidity` and `dLnKappaDLnRigidity` mean logarithmic derivatives
at fixed species and fixed `LocalState`. They include the exact relativistic
identity

```text
d ln(v) / d ln(R) = 1/gamma^2.
```

Consequently, constant lambda has a nonzero kappa rigidity slope, constant
kappa has a nonzero lambda slope of the opposite sign, and the empirical laws
include their configured independent-variable slope plus the speed term where
appropriate. `gradKappaParallelMPerS` is a Cartesian derivative at fixed
particle momentum and has units m s⁻¹. Constant models return an explicit
zero vector. Non-constant PD02 models leave it absent and set
`DerivativeUnavailable`, because their background factors do not yet carry
coherent gradients. A mover must retain its provider-consistent stencil or
reject a complete derivative-dependent operation; it must not replace absence
with zero.

| Status | Meaning in the current release |
| --- | --- |
| `InvalidParticle` | Mass/momentum is non-positive, charge is zero/non-finite, or relativistic conversion is unrepresentable. |
| `InvalidBackground` | A selected local value exists but is non-finite, zero, or negative where the physical decomposition requires positivity. |
| `MissingInput` | A selected runtime factor, field, nucleon count, or required parser key is absent. |
| `OutsideModelDomain` | Positive output overflows/underflows the representable finite domain; later models also use this for stated fit domains. |
| `InvalidConfiguration` | A value, pair of keys, Boolean, duplicate key, or unknown key violates the selected schema. |
| `UnsupportedModel` | The stable identifier is unknown or is registered but its roadmap backend is not implemented. |

`InfiniteMeanFreePath`, `IntegrationFailed`, `NonlinearSolverFailed`, and
`InconsistentSpectrum` are reserved for later scattering, quadrature,
nonlinear, and spectrum stages. Numerical status and diagnostic bits are
separate: future backends may return a successful finite coefficient together
with a physical-validity warning. The current explicit models set only
`DerivativeUnavailable`; they cannot diagnose the diffusion limit without
additional host-provided scales.

## Parser bridge and input examples

`BuildConfiguration` and `ConfigureActiveModel` accept an exact model ID and a
vector of `InputParameter{name,value}`. They reject duplicate keys, unknown
keys, missing paired parameters, malformed numbers, non-finite values, and
unsupported models. Unit conversions belong at the application input boundary;
the bridge accepts only the SI-labelled keys documented above.

The following is the planned `srcSEP3D` INI spelling. It documents the exact
library keys, but it is not accepted by the application until the PD11 parser
binding in `INTEGRATION_PLAN.md` is implemented:

```ini
[parallel_diffusion]
model = power_law_lambda
lambda0_m = 1.495978707e10
independent_variable = rigidity
rigidity0_V = 1.0e9
independent_exponent = 0.3333333333333333
radius0_m = 1.495978707e11
radial_exponent = 0.5
field_reference_T = 5.0e-9
field_exponent = 1.0
use_time_factor = false
use_region_factor = false
```

The equivalent planned legacy `srcSEP` block is:

```text
ParallelDiffusion on
model = power_law_lambda
lambda0_m = 1.495978707e10
independent_variable = rigidity
rigidity0_V = 1.0e9
independent_exponent = 0.3333333333333333
radius0_m = 1.495978707e11
radial_exponent = 0.5
field_reference_T = 5.0e-9
field_exponent = 1.0
use_time_factor = false
use_region_factor = false
```

These numbers demonstrate syntax and SI units only. They are not calibrated
defaults and are not asserted to describe any SEP population or event.

Other complete parameter examples are:

```ini
[parallel_diffusion]
model = constant_lambda
lambda_parallel_m = 1.495978707e10
```

```ini
[parallel_diffusion]
model = constant_kappa
kappa_parallel_m2_per_s = 1.0e18
```

```ini
[parallel_diffusion]
model = broken_rigidity_kappa
K_star_m2_per_s = 1.0e18
rigidity0_V = 1.0e9
break_rigidity_V = 3.0e9
low_slope = 0.3
high_slope = 1.8
smoothness = 2.0
field_reference_T = 5.0e-9
field_exponent = 1.0
use_radial_factor = false
use_region_factor = false
```

```ini
[parallel_diffusion]
model = bohm
eta_B = 1.0
field_definition = mean_field
```

## Numerical implementation

Particle speed, gamma, kinetic energy, and rigidity are calculated once from
total momentum with the exact relativistic relations in Equations (7),(8).
The kinetic-energy expression is algebraically rationalized to avoid
nonrelativistic cancellation. Multiplicative laws are evaluated in logarithmic
form. The broken-rigidity transition uses log-sum-exp rather than direct powers
of rigidity. No scattering floor, coefficient clamp, extrapolation, fallback,
or inferred turbulence quantity is present.

The final exponential is still checked for a finite positive result, so the
logarithmic formulation improves intermediate stability without pretending
that an arbitrarily large physical coefficient fits in binary64. Magnetic
vector magnitudes use nested `hypot`, and the broken-law derivative evaluates
its logistic transition with a bounded exponent after reaching the
double-precision asymptote. These transformations are algebraic; they do not
change the configured model or add a numerical regularization parameter.

The configuration fingerprint is a deterministic 64-bit FNV-1a identity of
the active, canonical parameter values. It is deliberately named a
fingerprint, not SHA-256. Later coefficient data sets retain their required
SHA-256 identities separately.

## Build and tests

From this directory:

```sh
make -f makefile verify
```

This builds `libparallel_diffusion.a`, runs the aggregate C++ test program, and
writes `build/test-report.json`. Individual test groups are currently in one
small executable; every record prints `PASS` or `FAIL` and records its model,
fixture/identity, and reason in the JSON report. A required failure gives a
nonzero process exit status.

The PD01/PD02 checks cover relativistic units/species, both constant laws,
reference normalization, separable power scaling, missing nucleon data,
`H(R0)=1`, the exact `a=b` limit, the beta contribution to rigidity slope,
equal-rigidity species behavior, Bohm scaling, selective field requirements,
`K_star` conversion, duplicate/invalid input, transactional pointer selection,
and explicit rejection of unavailable models. Mathematical test inputs are
not observational calibrations.

Coefficient expectations are constructed outside the selected production
model evaluator, while the shared kinematic conversion they sometimes consume
has its own independent fixture. The 10 MeV fixture constructs momentum from
kinetic energy in long-double test arithmetic; constant and Bohm checks apply
their defining identities directly; power-law checks use analytically chosen
ratios; and the broken-law derivative is checked with refined centered
differences plus Richardson cancellation. The tests therefore detect mistakes
in normalization, beta ownership, and derivatives instead of merely
reproducing implementation arithmetic.

## Qualification status and limitations

The companion `parallel_diffusion_model_data/` bundle named in Section 22 is
absent from this checkout. Consequently PD00 cannot verify its SHA-256 sums,
the full-precision `benchmark_points.json` fixture is unavailable, and the
coefficient tables/audit inputs required by PD05--PD08 cannot be imported.
Tests use independent algebraic identities and the printed Section 15.2
values at a tolerance compatible with their twelve significant digits; they
do not claim the unavailable full-precision-fixture gate.

See `IMPLEMENTATION_STATUS.md` for exact stage state and evidence,
`INTEGRATION_PLAN.md` for the planned `srcSEP`/`srcSEP3D` Parker bindings, and
`info_request.md` for the missing assets, scientific choices, provider
conventions, and calibration guidance still required for later stages.
