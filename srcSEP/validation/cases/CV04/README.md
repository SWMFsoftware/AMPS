# CV04 — adiabatic momentum change in expanding solar wind

## Purpose and isolated physics

CV04 isolates the Parker-mover momentum law `dp/dt=-(p/3) div(U)`. Spatial
diffusion, pitch-angle physics, sources, losses, and wave coupling are disabled.
The sweep covers protons and alpha particles at 0.1, 10, and 1000 MeV per
nucleon, so unit, mass-number, and relativistic energy-conversion errors cannot
hide in one nonrelativistic example.

Two backgrounds are exercised. Constant positive divergence has the exact
solution `p=p0 exp[-div(U)t/3]`. Constant-speed spherical wind has
`div(U)=2U/r`, `r=r0+Ut`, and `p=p0(r/r0)^(-2/3)`. The native spherical case
evaluates divergence at the radial midpoint and is expected to converge at
second order; the constant frozen coefficient is integrated exponentially.

## Reference, numerics, and acceptance

`reference_solution.py` evaluates both characteristics and its own
relativistic momentum/kinetic-energy transformations. Both implementations use
the algebraically stable rationalized expression for low-energy conversion,
avoiding cancellation in `sqrt(1+x^2)-1`.

The linked model uses three timesteps per scenario (20/10/5 s for constant
divergence and 200/100/50 s for spherical flow) with ten saved intervals. Gates
cover constant and spherical momentum error, spherical refinement order,
relativistic energy error on the finest realizations, and statistical-weight
conservation. A weighted momentum power law is transported with each energy;
the fitted final momentum slope must retain its configured index. Thresholds
and SI constants are explicit in `input.json`.

## Run and evidence

```sh
python3 test/run_tests.py --amps /absolute/path/to/amps \
  --validation-case CV04 --output-dir test_output/CV04
```

The linked registry writes `CV04_model.csv` plus native JSON/JUnit. The case
adds `CV04_reference.csv`, `CV04_solution.csv`, hashes/logs, and the proton
spherical energy comparison in PNG and EPS. A missing/stale executable fails
before reference generation; direct native invocation without reviewed case
paths reports SKIP.

## Interpretation and limitations

Species-dependent momentum error suggests an incorrect total mass or
per-nucleon conversion. Correct momentum but wrong energy implicates the
relativistic output transform. Loss of second-order spherical convergence
indicates coefficient evaluation at the wrong radius. CV04 validates adiabatic
characteristics only; it does not validate a self-consistent solar-wind model,
shock acceleration, spectra after transport, or spacecraft observations.
