# Pitch-angle diffusion coefficient library (Release A)

Host-neutral C++17 implementation of the prescribed-scattering subset
("Release A", roadmap Section 26.1) of
[`PITCH_ANGLE_SCATTERING_MODEL.md`](PITCH_ANGLE_SCATTERING_MODEL.md),
specification v1.2 (13 October 2026; SHA-256
`1a7ecf1a599ea54a24e5e0da061a0d71c00ed75f4b74e3b561262d724ea96b92`, the
digest the supplement's `review_manifest.json` gives for the reviewed v1.2).
It supplies the pitch-angle diffusion coefficient D_μμ, its derivative and a
scattering step limit to a focused-transport particle pusher, and the pure
scattering kernels and reference operator needed to verify such a pusher.
srcSEP3D's focused-diffusion mover uses it through
`pitch_angle_diffusion_model = pitch-angle-library` (see
`srcSEP3D/README.md`, "Shared pitch-angle diffusion"); the srcSEP binding is
deferred.

The stage ledger with exact commands, measured errors and open items is
[`IMPLEMENTATION_STATUS.md`](IMPLEMENTATION_STATUS.md).  It distinguishes what
is implemented from what is independently validated.

Like the parallel-diffusion library next to it, this library has no PIC, MPI
or application dependency. Every value crossing its boundary is SI.

## Contents

| File | Purpose |
|---|---|
| `pitch_angle_diffusion.h/.cpp` | registry, kinematics, shapes, λ∥ normalization, step rules, parser, active model and dispatch pointers |
| `pitch_angle_certification.cpp`, `pitch_angle_certified_pairs.inc` | the coefficient-kernel certification matrix (PS05 evidence) and `CheckKernelPair` |
| `pitch_angle_kernels.h/.cpp` | Euler-Maruyama, Milstein and random-rotation kicks; mirror folding; finite-volume reference operator, its exact propagator and Green-Kubo variance |
| `pitch_angle_quadrature.h` | internal adaptive Gauss-Kronrod quadrature |
| `pitch_angle_registry.inc`, `generate_registry.py` | the 172-ID inventory generated from specification Section 4.2 |
| `test_pitch_angle_diffusion.cpp` | deterministic unit tests PAD01-PAD12 |
| `pitch_angle_driver.cpp`, `reference_checks.py` | independent reference checks PAR01-PAR11 |
| `test_pitch_angle_montecarlo.cpp` | statistical acceptance campaign (AP-* cases of Section 26.5) |
| `run_all_tests.py` | global runner of all of the above |
| `pitch_angle_diffusion_model_data/` | the v1.2 implementation supplement (roadmap, contract header, `check_physics.py`); not a data archive |

## What is implemented

### Operator and convention

Every coefficient is the D_μμ of the standard operator (Eq. 21)

```
(df/dt)_scat = d/dmu ( D_mumu(mu) df/dmu ),    -1 <= mu <= 1,   D_mumu >= 0 [s^-1]
```

and a forward Itô pusher realizes it with the drift dD_μμ/dμ and the noise
√(2 D_μμ) (Eq. 25), so both are returned together.  A source printed with
another operator factor is converted once, inside the model (Section 3.3).
The library never adds streaming, focusing, flow or adiabatic terms; those
belong to the host's transport composer (Section 16.12).

### Models (shapes D = D0 g(μ), g = (1 − μ²) h(μ))

| `model` | h(μ) | Parameters | D0 from λ_law (Eq. 23, Section 5.6) | Behaviour at μ = 0 |
|---|---|---|---|---|
| `PA-ISO` | 1 | — | D0 = v/(2λ) | regular |
| `PA-QFORM` | \|μ\|^(q−1) + H | `q` (1 < q < 2), `gap_h` (H ≥ 0) | D0 = (3v/4λ) I(q,H) | unbounded drift (Eq. 163) |
| `PA-EPS` | \|μ\|/(1+\|μ\|) + ε | `epsilon` (ε > 0) | ν0 = vφ(ε)/λ, D0 = ν0/2 | drift jump −D0 → +D0 (Eq. 164) |
| `PA-EPREM` | 1 | — | printed D^E = (1−μ²)v/(2λ) with operator D^E/2: standard D0 = v/(4λ); Eq. (23) returns 2λ | regular |

* I(q,H) = ∫₀¹(1−μ²)/(μ^(q−1)+H) dμ: exact 2/((2−q)(4−q)) at H = 0, Eq. (38)
  at q = 5/3, otherwise the second form of Eq. (37) (2F1(1,1;c;z),
  z = 1/(1+H), summed with a rigorous geometric tail bound).
* φ(ε): the closed form Eq. (39).
* Each closed form is checked against an independent adaptive quadrature at
  configuration time.  If they differ by more than 1e-8 relative, the
  configuration fails with `Nonconverged` (Section 5.7).
* ε = 0 is refused for a λ-normalized PA-EPS, because λ∥ diverges as
  (3/2)ln(1/ε). It is accepted only by the amplitude-level `EvaluateShape` and
  the finite-volume reference, where it disconnects the hemispheres.
* H, h, h₀ and ε are never converted into one another (Section 5.7).

### The λ∥ amplitude law (user decision 2026-10-13)

```
lambda_law = lambda_ref (R/R_ref)^a (r/r_ref)^b,   R = p c/|q| [V],  r = |x - x_sun| [m]
```

The law is evaluated per particle from the pusher's state: mass, signed
charge, momentum and heliocentric position. The returned
`lambdaParallelUnfocusedM` is the tagged Eq. (23) value of the actual
coefficient: equal to `lambdaLawM`, except twice it for PA-EPREM. No focusing
relation is applied (U-46).

### Scattering step rules (Section 16.6; user decisions 2026-10-13)

The rule is selected in the input; none has a default.

| `step_rule` | Δt | Parameters |
|---|---|---|
| `rate-fraction` | f/ν_max | `rate_fraction` f |
| `van-den-berg-2026` | max( min{ ℓ√(1−μ²)/\|D′+dμ/dt\|, ℓ²(1−μ²)/(2D) }, f/ν_max ) | `ell` ℓ (source: 2⁻⁶), `rate_fraction` f (floor) |
| `sandroos-2014` | min( Δμ_max/a_max, Δμ_max²/b_max² ) | `delta_mu_max` (0 < Δμ_max ≤ 2) |

* ν_max = sup_μ 2D/(1−μ²) = 2 D0 sup h is the largest local scattering
  frequency.
* a_max = sup_μ\|D′\| and b_max² = sup_μ 2D, both at the particle's p and r
  (the printed "max").
* The van den Berg rule is written for forward time: its drift entry uses the
  forward Itô drift D′ + dμ/dt, while the source prints the backward one.
* Without the floor, Eq. (165) is exactly 0 at \|μ\| = 1 whenever D′(±1) ≠ 0,
  so an endpoint beam would never advance (U-55). The rate-fraction floor
  prevents that.
* A rule never returns a zero-length success.

### Coefficient-kernel readiness (Section 16.11)

Release is a property of a (model, kernel) pair and of the step it runs with.
`CheckKernelPair(configuration, kernel)` accepts only rows of the
certification matrix (`pitch_angle_certified_pairs.inc`). The step is
measured rule-independently as the largest dimensionless step the configured
rule can produce, ν_max Δt (`WorstCaseScaledStep`):

* `rate-fraction`: f;
* `van-den-berg-2026`: max(f, ℓ² sup h / inf h);
* `sandroos-2014`: min(Δμ_max · 2 sup h / sup\|g′\|, Δμ_max² sup h / sup g).

The matrix and its evidence are in `IMPLEMENTATION_STATUS.md` (stage PS05).
In summary:

* **PA-QFORM** with Euler-Maruyama or Milstein is always refused, because its
  drift is unbounded at μ = 0 for every H (user decision).
* **Euler-Maruyama with mirror folding** is not certified. It fails the
  endpoint-beam case AP-ISO-01 at every tested step: the reflection bias of
  Section 16.4 decreases too slowly with Δt.
* **Milstein** is certified for PA-ISO and PA-EPREM (the same angular shape,
  h = 1) at ν_max Δt ≤ 0.001: every Milstein case passes there (worst budget
  use 50.9%, AP-FOC-01 ξ = −5). At 0.003 it fails one focusing case (ξ = +5,
  101% of the 1% variance budget) and at 0.01 the endpoint-beam case.
* **PA-EPS with Milstein** is not certified. In AP-ZERO-02 the drift jump at
  μ = 0 (D′ = ±D0) makes the scheme over-populate |μ| < 0.01: from a uniform
  start the exact band probability 0.01 is measured as 0.0140, 0.0119 and
  0.0107 at ν_max Δt = 0.01, 0.003 and 0.001 for ε = 0.01 (budget 5e-4). An
  independent numpy Milstein reproduces the value (0.013991 ± 1.6e-5), so this
  is the scheme's bias, not a coding error. ε = 1 is inconclusive: the
  particle counts of the concentrated starts were too small for the design
  check. PA-EPS remains evaluable (coefficient, λ∥, finite-volume reference).
* The random rotation passes AP-ISO-01/02 as a library kernel for the
  isotropic law, but no focusing case was run with it and no host uses it, so
  it has no matrix row.
* The finite-volume chain is a reference operator, not a pusher kernel.

### Frame and clock (user decision 2026-10-13)

* The coefficient is evaluated with p and μ in the local plasma
  (solar-wind) frame, the frame in which the host stores them.
* It is applied per unit host time, retaining order 0 in U/c: the exact clock
  factor of Eq. (173), γ_u(1 − u v μ/c²), differs from one by about 10⁻³ at
  U = 400 km/s and v ~ c, and is not applied (`FrameContract`).
* A host that needs it must multiply the whole generator (Eq. N7). It must not
  rescale D and then differentiate.

### Registry and refusals

* `Registry()` holds all 172 inventory IDs of Section 4.2: ID, aliases,
  section, equations and runtime state, generated by `generate_registry.py`
  and checked against the specification digest.
* A registered ID outside Release A is refused with `UnsupportedModel`, naming
  its section and runtime state (for example
  `QLT-SLAB (Section 6.1, READY_EXPLICIT_INPUTS) is registered but outside
  Release A`).
* An ID absent from the inventory is `InvalidConfiguration`.
* No other model is ever substituted.

## Input block

A host passes the body of its block (srcSEP3D:
`#begin physics: pitch_angle_diffusion`) to `ParseSection` with physical line
numbers:

```text
#begin physics: pitch_angle_diffusion
  model = PA-ISO                        ! inventory ID, case-insensitive
  lambda_reference_m = 4.48793612e10    ! lambda_ref [m]
  lambda_reference_rigidity_v = 1e9     ! R_ref [V]
  lambda_rigidity_exponent = 1/3        ! a (exact integer ratio allowed)
  lambda_reference_radius_m = 1.495978707e11   ! r_ref [m]
  lambda_radial_exponent = 2/3          ! b
  step_rule = rate-fraction
  rate_fraction = 0.001                 ! certified limit for PA-ISO/Milstein
#end physics
```

PA-EPS takes `epsilon` (ε > 0), PA-QFORM takes `q` and `gap_h`; PA-ISO and
PA-EPREM take no shape key. Parsing a model does not release it with a
kernel: `CheckKernelPair` does that.
`van-den-berg-2026` takes `rate_fraction` and `ell`; `sandroos-2014` takes
`delta_mu_max` only.

* **Comments:** `#` and `!` start a comment.
* **Keys:** they are case-sensitive and name SI units.
* **Values:** decimal SI numbers. Dimensionless keys also accept an exact
  integer ratio `a/b` (Section 26.3 PS01).
* **Failures:** unknown, duplicate, missing and inactive keys fail, and the
  diagnostic begins with `line N:`.

## API sketch

```cpp
#include "pitch_angle_diffusion/pitch_angle_diffusion.h"
namespace PA = SEP::PitchAngleDiffusion;

PA::ParsedSection parsed;                               // serial initialization
PA::Status s = PA::ConfigureActiveModelFromSection(lines, &parsed);
s = PA::CheckKernelPair(parsed.configuration, PA::Kernel::Milstein);

// Per particle substep (pure, reentrant):
PA::ParticleState particle{mass_kg, charge_C, momentum_kg_m_s};
PA::LocalState local;  local.heliocentricPositionM = {x, y, z};
PA::CoefficientResult c = PA::ActivePitchAngleDiffusion(mu, particle, local);
if (!c.status.ok()) { /* typed failure; never substitute a value */ }
double D = c.dMuMuPerS, Dprime = *c.dDmuMuDmuPerS;      // D' absent => no finite drift
PA::StepLimitResult lim = PA::ActiveScatteringStepLimit(mu, dmu_dt, particle, local);
```

The kernels in `pitch_angle_kernels.h` take the random variates as arguments,
so the host keeps RNG ownership and restart state:

```cpp
SEP::PitchAngleDiffusion::Kernels::MilsteinKick(mu, D, Dprime, dmu_dt, dt, zeta)
```

## Build and test

```bash
cd src/models/sep_common/pitch_angle_diffusion
make                       # libpitch_angle_diffusion.a
make tests                 # unit, reference-driver and acceptance binaries
python3 run_all_tests.py   # everything (reuses an unchanged acceptance run)
python3 run_all_tests.py --rerun-acceptance   # force the campaign (52 cases, ~2.3 h
                                              # on 16 idle cores; 3.9 h measured on a
                                              # shared node)
python3 run_all_tests.py --skip-acceptance    # fast stages only
```

`reference_checks.py` and the supplement's `check_physics.py` need numpy,
scipy and mpmath (`python3 -m pip install --user mpmath`). A missing
dependency is reported as SKIP, never as PASS. srcSEP3D builds the library
through its makefile (`pitch-angle-diffusion` target) and flattens the three
objects into `mainlib.a`. `srcSEP3D/test/run_all_test.sh` runs the global
runner as step `pitch-angle-library`.

### Validation criteria

| Stage | Tests | Reference | Tolerance |
|---|---|---|---|
| Unit | PAD01-PAD12 | Section 18 fixtures (F-EMP, F-PUSH, F-FOC, F-GCR-PA), Section 5.6 table, exact identities, tanh-sinh quadrature, finite differences | closed forms 1e-12; quadrature class 1e-9; 7-digit table values 5e-7 |
| Reference | PAR01-PAR11 | mpmath (40 digits; 2F1 and substituted quadrature cross-checked), mean-free-path and parallel-diffusion companion bundles, supplement `check_physics.py` (261 checks), registry digest, compile check of the C++17 contract | 1e-12 (1e-9 transmissibility, 1e-10 Eq. 23) |
| Acceptance | AP-ISO-01/02, AP-TIME-01, AP-ZERO-02, AP-FOC-01 | exact solutions (Eqs. 174, 175, 177, 176, F-FOC-04) and the finite-volume operator (fine grid, error from coarse/fine difference) | Section 26.5 budgets: Legendre and hemisphere 5e-3, occupation of \|μ\| < 0.01 5e-4, displacement and focusing variance 1%, focusing mean 5e-3; criterion \|err\| + z·SE + ε_ref ≤ ε and z·SE ≤ ε/4 (family-wise 99%, Bonferroni) |

## Known limitations and open items

* **PS01: canonical data archive.** It is BLOCKED.
  `PITCH_ANGLE_SCATTERING_MODEL_DATA` (registry with admissibility/capability
  per ID, 89 fixtures) is not available. Only the Section 4.2 fields are
  compiled. The tests use the Section 18 values, the companion bundles and
  independent references instead.
* **Release B and C** (spectral QLT, dynamical, nonlinear, wave, compressible,
  measured-spectrum and self-generated-wave models) are registered and refused.
* **Not implemented:** tensor (μ,p) coefficients, frame transformations and
  clock factors, backward pushers and importance sampling, the hemispheric and
  redraw kernels (U-13, U-63), and a Brownian bridge. No adaptive rejection is
  used, so no bridge is needed.
* **PA-QFORM** is evaluable, with D, D′ away from 0, λ∥ and finite-volume
  transmissibility, but no pusher pair is released.
* **Euler-Maruyama** is not certified, see above.
* **PA-EPS** has no released pusher pair (AP-ZERO-02, see above). A
  derivative-free or jump-aware kernel, or Milstein at a smaller step with a
  correctly sized campaign, would need new evidence.
* **AP-ZERO-02 design defect.** The particle counts of its concentrated
  ("zero", "near-zero") starts were sized for a band probability near 0.01,
  but it is much larger there, so those comparisons fail the z·SE ≤ ε/4
  design check. The counts were frozen and not changed after the run; the
  PA-EPS verdict rests on the uniform start, whose sizing is correct.
* **Host lifecycle evidence** (AP-GEOM-01..04: migration, AMR,
  restart-mid-subcycle with the library active) is not produced by this
  library. srcSEP3D's component tests cover the bridge, the kick equality and
  the composer. Native evidence is listed in `IMPLEMENTATION_STATUS.md`.
* **Fixed step for observables.** With the van den Berg rule the step depends
  on μ. srcSEP3D samples at common AMPS step times, not at substep endpoints,
  as Section 16.6 requires.
* **No observational validation.** The prescribed shapes are mathematically
  admissible laws. Their agreement with measured SEP or GCR scattering is a
  property of each published parameterization, not established here.
