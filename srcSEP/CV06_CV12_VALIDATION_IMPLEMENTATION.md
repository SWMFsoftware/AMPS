# CV06–CV12 linked-validation implementation

## Execution and evidence contract

Every case is a native descriptor in `component_tests.cpp` and is linked into
the srcSEP/AMPS application through `advanced_validation_models.o`. The public
runner first requires the selected executable to advertise the requested ID,
then invokes `amps --test CVxx --test-input ... --test-output-dir ...`. A native
PASS and nonempty raw model CSV are required before an independent Python
reference is run. No standalone harness result is accepted as scientific
evidence.

Each case preserves the resolved unit-bearing JSON, strict native argument
manifest, raw model/reference/solution CSVs, native and aggregate JSON/JUnit,
command log, executable/input/output hashes, and PNG/EPS figures. Case-local
READMEs state the equations, sweeps, gates, and limitations.

## CV06: Legendre pitch-angle diffusion

The production focused `Dmumu` mover evolves separately initialized modes
`l=1..6` for `D_mumu=D0(1-mu^2)`. All measured modes 0–6 are saved. The scorer
fits decay rates against `l(l+1)D0`, reports leakage and normalization in
standard-error units, calculates a timestep refinement diagnostic, and retains
the mover's boundary-reflection count. The finite ensemble makes this a weak
stochastic convergence gate; larger publication campaigns should increase the
registered particle and seed counts.

## CV07: telegraph transport

Production keyed exponential waiting times drive a symmetric two-velocity
persistent flight. Three rates and early/transition/late normalized times save
front mass, binned profiles, event count, MSD, kurtosis, and causal-support
violations. Independent closed forms check the unscattered fronts, MSD, Poisson
events, and the late `v^2/(2nu)` diffusion limit.

## CV08: drift-diffusion first passage

The Parker core advances constant drift and diffusion to an absorbing boundary.
Every particle retains arrival/censor status, interpolated arrival time, and
overshoot. Ten seeds in two Peclet regimes support per-seed KS tests against the
inverse-Gaussian CDF, moment uncertainty, overshoot refinement, and a wrong-
diffusion negative control.

## CV09: planar DSA

A controlled linked shock-cycle model samples constant return/escape and gain
for compression ratios 2, 3, and 4. An O(N) ranked survival fit estimates the
phase-space index and compares it with `3r/(r-1)`; acceleration time, particle
accounting, and diffusion-length resolution are checked independently. This is
not presented as validation of a fully resolved PIC shock geometry.

## CV10: nonuniform turbulence advection

The production spectral turbulence state transports both `U+V_A` and `U-V_A`
branches for fixed-area, expanding-area, and coarse-to-fine remapped grids.
Cell-integrated sinusoidal characteristics provide the reference. Three spatial
resolutions measure L2 order, while all bins and scenarios contribute to wave-
action conservation and positivity gates.

## CV11: resonant growth and damping

A one-hot spectral state receives constant growth, constant damping, exact
cancellation, and sign-changing source histories through the production energy
ledger. The sinusoidal rate uses midpoint quadrature and is compared with an
exact time integral at three steps, supplying a real second-order test. Inactive
bins, positivity, and zero-rate invariance are retained as explicit metrics.

## CV12: closed particle-wave exchange

Production wave-frame scattering determines the particle energy change; the
turbulence ledger receives the opposite transfer in the selected resonant
branch. Particle-only, wave-only, coupled, and suppressed-deposition runs span
three timesteps and three population sizes. Coupled total energy and per-step
ledger residual must close, uncoupled states must remain fixed, and the negative
control must be detected.

## Commands

```sh
# Strict source, Python, registry, and documentation-independent smoke gate.
make test-cv06-cv12-unit

# Run the actual linked application and generate the complete evidence set.
make test-cv06-cv12-unit SEP_EXECUTABLE=/absolute/path/to/amps

# Equivalent explicit user-facing campaign selection.
python3 test/run_tests.py --amps /absolute/path/to/amps \
  --validation-case CV06 --validation-case CV07 \
  --validation-case CV08 --validation-case CV09 \
  --validation-case CV10 --validation-case CV11 \
  --validation-case CV12 --output-dir /absolute/path/to/evidence
```
