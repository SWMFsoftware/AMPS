# CV07 — telegraph transport and diffusion limit

CV07 validates finite-speed persistent transport before its diffusion limit.
Particles move at `+v` or `-v`; production exponential waiting-time sampling
switches direction with Poisson rate `nu`. The exact causal support is
`|s|<=vt`, the unscattered front mass is `exp(-nu t)`, and
`<s^2>=v^2[t/nu-(1-exp(-2nu t))/(2nu^2)]`. At late time the effective spatial
diffusion coefficient approaches `v^2/(2nu)`.

Three switching rates and four values of `nu*t` span the ballistic-front,
transition, and diffusion regimes. Ten seeds retain uncertainty estimates;
profiles, front counts, events, MSD, kurtosis, and support violations are saved
before scoring. The independent reference calculates the closed-form moments
and front weights. Gates use binomial/seed standard errors, require zero mass
outside the causal cone, and verify both event statistics and the late-time
diffusion coefficient. The plot shows MSD across the regime sweep.

Run `python3 test/run_tests.py --amps /path/to/amps --validation-case CV07
--output-dir test_output/CV07`.
