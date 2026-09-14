# CV02 — constant-coefficient spatial diffusion Green function

## Purpose and isolated physics

CV02 verifies the stochastic spatial-diffusion term in the production Parker
mover. Solar-wind advection, focusing, momentum change, sources, losses, and
particle-wave coupling are zero. A failure therefore points to the
`sqrt(2 kappa dt)` increment, random-stream use, normalization, or statistics
rather than an interacting heliospheric process.

The linked application advances a point packet with `kappa_parallel=10^12
m^2/s` on a long uniform line. At 200 s its standard deviation is only 20 Mm;
the packet center is 25 standard deviations from either boundary, satisfying
the plan's minimum five-sigma isolation criterion.

## Reference and acceptance

The independent solution is

`f(s,t) = [4 pi kappa t]^(-1/2) exp[-(s-s0)^2/(4 kappa t)]`,

with mean `s0`, variance `2 kappa t`, skewness zero, and kurtosis three.
`reference_solution.py` evaluates exact error-function integrals over the same
80 finite bins as the model; it does not sample the density at bin centers.

Ten independent seeds, particle counts 1,000/4,000/16,000, and timesteps
20/10/5 s expose Monte Carlo scaling and timestep dependence. The gate fits
kappa from the variance-time slope, measures the mean in standard errors,
compares the pooled final profile in L1/L2/L-infinity, checks skewness,
kurtosis and escaped weight, and requires a wrong `1.5*kappa` negative control
to be clearly worse. Thresholds are in
`input.json`, not hidden in code.

## Run and evidence

```sh
python3 test/run_tests.py --amps /absolute/path/to/amps \
  --validation-case CV02 --output-dir test_output/CV02
```

The binary must advertise `CV02`; the numerical stage is
`amps --test CV02 --test-input ... --test-output-dir ...`. Outputs include the
native JSON/JUnit and `CV02_model.csv`, independent `CV02_reference.csv`, the
standard `CV02_solution.csv`, provenance hashes, logs, and
`CV02_comparison.png/.eps`. JSON metrics are authoritative; figures are review
aids. A direct registry run without test paths reports SKIP because it lacks a
reviewed configuration.

## Interpretation and limitations

A biased mean suggests unintended drift; incorrect variance slope suggests a
factor-of-two or unit error; a profile-only error suggests bin normalization or
non-Gaussian increments. CV02 is controlled numerical verification, not an
observational or SWMF-coupling claim. Statistical thresholds apply only to the
documented ensemble and must not be loosened after observing a failed run.
