# Controlled analytical component validation

This document describes the common registry cases used to test individual
transport movers and turbulence operators in deliberately simplified physical
environments. The linked `amps` CLI and the source-only sanitizer runners call
the same descriptor callbacks in `util/sep_mover_validation.cpp` and
`util/sep_turbulence_validation.cpp`; the expected solution is therefore not a
separate Make-only implementation.

## Selecting tests

In a configured AMPS checkout, inspect and run the catalog with:

```sh
../amps --list-tests
../amps --test PARK04
../amps --test FTED08 --test-json pitch-mode.json
../amps --test TURB21 --test TURB22 --test TURB23 --test-junit turbulence.xml
../amps --test-group parker
../amps --test-group fte-dmumu
../amps --test-group fte-mfp
../amps --test-group turbulence
```

Explicit group selection includes extended statistical cases. `--all-tests`
continues to select only descriptors marked `routine`; `PARK02`, `PARK07`,
`FTED01`, `FTED08`, `FTEM01`, `FTEM02`, and `FTEM08` are `extended`.
Source-only equivalents are:

```sh
make test-parker-unit
make test-fte-dmumu-unit
make test-fte-mfp-unit
make test-turbulence-core-unit
```

All focused runners compile with C++11, strict warnings, AddressSanitizer, and
UndefinedBehaviorSanitizer. Their JSON and JUnit reports are created only in a
disposable build directory.

## Parker cases

| ID | Controlled equation or boundary | Acceptance quantity |
|---|---|---|
| `PARK01` | (ds/dt=U), `kappa=0` | (s-s_0=U\Delta t) |
| `PARK02` | (ds=\sqrt{2\kappa},dW) | mean zero and variance (2\kappa\Delta t) |
| `PARK03` | Itô gradient drift | (\Delta s=(d\kappa/ds)\Delta t), including sign |
| `PARK04` | (dp/dt=-p\,\nabla\cdot U/3) | (p=p_0e^{-(\nabla\cdot U)t/3}) |
| `PARK05` | outward step through a finite interval | explicit absorbing-boundary status |
| `PARK06` | (ds/dt=0.4s) | error reduction toward (s=e^{0.4t}) |
| `PARK07` | drift-free Brownian motion on `[-1,1]` | exact upper-exit probability and mean first-passage time |

The analytical references use scalar arithmetic outside `AdvanceParker`. The
statistical case records its fixed keyed seed, ensemble size, mean error, and
variance error.

For `PARK07`, the particle starts at `x=0.2 m` with
`kappa=0.25 m^2 s^-1`. The independent absorbing-boundary references are

```text
P(X_T=1) = (x-a)/(b-a),
E[T] = (x-a)(b-x)/(2 kappa).
```

The test follows 20,000 fixed-seed paths with `dt=0.002 s`, applies the shared
absorbing-coordinate kernel after each production Parker increment, and records
the completion fraction. Its statistical limits use five standard errors plus
an explicit allowance for discrete end-of-step boundary detection.

## Coefficient-driven focused-transport cases

`FTED01`–`FTED07` cover constant-diffusion Itô moments, derivative drift,
magnetic focusing, combined streaming/cooling, reflective pitch boundaries,
QLT normalization and turbulence identity, and deterministic split
refinement. `FTED08` adds a continuum eigenmode benchmark.

For

```text
D_mumu = D0 (1-mu^2),
f(mu,0) = 0.5 [1 + a P2(mu)],
```

the pitch-angle Fokker–Planck operator has Legendre eigenfunctions and

```text
<P2>(t) = (a/5) exp(-6 D0 t).
```

`FTED08` independently rejection-samples the initial distribution, advances
40,000 production particles with seed 308, and evaluates the final mode with a
five-standard-error acceptance limit. The reference uses only Legendre
orthogonality and the analytical eigenvalue; it does not call a production
flux or reflection routine.

## Event-driven mean-free-path cases

| IDs | Controlled reference |
|---|---|
| `FTEM01` | Exponential waiting-time mean (1/\nu) |
| `FTEM02` | Poisson event-count mean (\nu t) |
| `FTEM03` | Infinite-mean-free-path characteristic (\Delta s=(U+v\mu)t) |
| `FTEM04` | Exact conservation of speed in the selected Alfvén-wave frame |
| `FTEM05` | Focusing convergence under interval refinement |
| `FTEM06` | Exact adiabatic momentum characteristic |
| `FTEM07` | Combined focusing/streaming/cooling refinement |
| `FTEM08` | Persistent-flight MSD and long-time spatial-diffusion limit |

For isotropic pitch resets occurring as a Poisson process of rate
`nu=v/lambda`, the stationary parallel-velocity correlation and exact MSD are

```text
<v_parallel(0)v_parallel(t)> = (v^2/3) exp(-nu t),
<Delta s^2> = [2 v^2/(3 nu^2)] [nu t - 1 + exp(-nu t)].
```

`FTEM08` evaluates this finite-time telegraph-like persistent-flight solution
at `nu*t=0.2`, then evaluates the same formula and the asymptotic
`<Delta s^2>=2 kappa_parallel t`, `kappa_parallel=v*lambda/3`, at `nu*t=50`.
The Alfvén speed is zero so wave-frame scattering is an exact isotropic reset
in the plasma frame; 40,000 fixed-seed histories are used.

## Automatic refinement order

All refinement cases use the common registry utility
`EstimateRefinementOrder`:

```text
p = log(e_coarse/e_fine) / log(h_coarse/h_fine).
```

`PARK06`, `FTED07`, `FTEM05`, `FTEM07`, and `TURB21` report `p` and the
coarse/fine error ratio in terminal, JSON, and JUnit evidence. The helper
rejects nonfinite or nonpositive errors and invalid resolution ordering. The
minimum expected orders are `0.85`, `1.80`, `1.80`, `0.85`, and `0.85`,
respectively, reflecting each scheme's formal first- or second-order behavior
with a small pre-asymptotic allowance.

## Turbulence cases

`TURB02`–`TURB20` are now registered individually. They cover particle-wave
ledger closure, units, initialization, spectral projection, uniform advection,
CFL subcycling, boundary validation, reflection, cascade/dissipation, growth,
coefficient closure, shock injection, operator order, positivity limiting,
conservative remap, restart, deterministic histories, mover/source separation,
and common standalone/coupled evolution. `TURBOWN01` retains the former focused
source-ownership check under a nonconflicting ID because the production catalog
already used `TURB01` for the 1-AU magnetic-pressure closure.

### TURB21: translated nonuniform profile

The periodic analytical solution is

```text
E(s,0) = 1 + 0.25 sin(2 pi s),
E(s,t) = 1 + 0.25 sin[2 pi (s-Ut)].
```

The test initializes exact finite-volume cell averages, not center samples,
and evolves the production first-order upwind operator at CFL 0.4 to `t=0.25`
with 32 and 64 cells. It requires decreasing L1 error, observed refinement
order `p>=0.85`, fine-grid L1 error no larger than
`0.025`, and conservation to `2e-13 J`. This detects incorrect propagation
speed or phase and numerical diffusion that the uniform `TURB06` invariant
cannot detect.

### TURB22: time-dependent growth history

The prescribed rate is

```text
G(t) = 0.30 + 0.20 t W,
E(t) = E0 + 0.30 t + 0.10 t^2.
```

At every 0.1 s interval the test supplies the exact integral of `G(t)` through
the production pending particle-wave source and compares every resulting state,
not only the final point. The maximum time-history error must be at most
`2e-13 J`.

This separation is intentional. `TURB22` proves that the turbulence update
applies a known time-dependent energy source once, accumulates it at the
correct time, and reports it consistently. It
does not test how a particle distribution is converted into `G(t)`. The native
`growth_rate_validation_test.cpp` comparison exercises that QLT/coupling path
and therefore remains an AMPS/PIC/MPI integration test rather than a
dependency-light analytical unit test.

### TURB23: controlled particle-wave total energy

Two deterministic wave-frame scattering cases exercise the outward-particle
minus branch and inward-particle plus branch. For each event the test computes
the relativistic particle kinetic-energy change independently, selects a
macro-particle weight giving `|Delta E_particle|=0.25 J`, and applies

```text
Delta E_wave = -Delta E_particle.
```

The production turbulence update must preserve weighted particle-plus-wave
energy and report the same signed increment in `particleExchangeJ`, both within
`2e-12 J`. The controlled exchange is smaller than the initial energy of each
branch, so a positivity-limiter activation fails the case instead of hiding an
incorrect sign or branch assignment.

## Validation gate boundary

The source-only analytical command is `make test-controlled-analytical`.
Native AMPS, real SWMF, and observational validation remain separate commands
and cannot borrow status from these controlled cases. See
`validation/README.md` for manifest schemas, checksum requirements, held-out
event policy, and the exact `test-native-amps-validation`,
`test-swmf-validation`, and `test-observational-validation` invocations.

## Interpretation and limitations

These tests establish equation-level behavior in controlled environments. They
do not establish accuracy for arbitrary spatially varying coefficients, a full
AMPS particle adapter, native MPI decomposition, a coupled SWMF run, or an
observed SEP event. Those remain separate integration and scientific-validation
gates. Statistical descriptors use fixed reproducible seeds; an eventual
campaign-level seed panel would provide a stronger assessment of rare
statistical failures.
