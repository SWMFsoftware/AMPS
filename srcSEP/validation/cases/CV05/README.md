# CV05 — magnetic focusing in a prescribed field gradient

## Purpose and isolated physics

CV05 verifies the deterministic focusing terms of the full gyrotropic
focused-transport mover. It prescribes `B(s)=B0 exp(g s)`, hence constant
`g=d ln|B|/ds`, while setting `D_mumu`, plasma flow, velocity gradients,
momentum change, sources, losses, and particle-wave coupling to zero.

Both signs of `g` are run with pitch-angle cosines near zero, near both bounds,
and exactly at `mu=+-1`. This exposes field-line orientation mistakes,
singular endpoint handling, sign errors, and leakage outside the physical
pitch-angle interval. Positive quadrature weights represent a smooth angular
distribution `f(mu)=1+0.4 mu` for moment checks.

## Exact characteristic and gate

For constant particle speed `v`,
`dmu/dt=-(v g/2)(1-mu^2)` and `ds/dt=v mu`. For interior angles,
`mu=tanh[atanh(mu0)-v g t/2]`; integrating that expression gives the reference
position. The `mu=+-1` solutions are evaluated by their finite limiting
characteristics. `reference_solution.py` implements these formulas in a
separate process.

The linked mover runs 20/10/5 s timesteps to 1000 s. Acceptance covers position
relative to the path length, pitch angle, momentum conservation, the adiabatic
magnetic-moment proxy `(1-mu^2)/B`, the first angular moment, second-order
the first two angular moments, total weight, second-order trajectory refinement,
and pitch-angle bounds. Every saved time participates in the characteristic
comparison. Values and tolerances are in
`input.json`.

## Run and evidence

```sh
python3 test/run_tests.py --amps /absolute/path/to/amps \
  --validation-case CV05 --output-dir test_output/CV05
```

The supplied binary must list and execute CV05. Artifacts include its native
JSON/JUnit and `CV05_model.csv`, the independent `CV05_reference.csv`, standard
solution CSV, executable/input/output hashes, exact command log, and PNG/EPS
pitch-angle comparison with residuals. Registry selection without the case
input/output paths reports prerequisite SKIP rather than inventing defaults.

## Interpretation and limitations

Errors that reverse with gradient sign indicate a focusing-sign convention;
large errors only near `mu=+-1` indicate endpoint treatment; invariant drift
with good `mu` but bad position indicates inconsistent streaming/focusing
coupling. CV05 is a prescribed smooth-gradient verification and does not cover
pitch-angle scattering, magnetic mirrors with turning points, discontinuous
fields, SWMF interpolation, or observational agreement.
