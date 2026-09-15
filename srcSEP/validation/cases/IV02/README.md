# IV02 — scattering–focusing equilibrium

IV02 couples focusing and pitch diffusion in one production mover. With
`D_mumu=D0(1-mu^2)` and a constant focusing strength, the zero-flux stationary
distribution is proportional to `exp(xi*mu)`. Ratios `xi=0.2,1,2` cover a
decade, while isotropic and beam-like initial populations test ergodicity.

The model stores normalized 40-bin PDFs and first moments for three timesteps.
The case-local reference analytically integrates the exponential over each
production bin. Gates cover L1 error, first moment, normalization, and agreement
between initial conditions. This validates the combined operator convention;
CV06 remains the sharper isolated `D_mumu` eigenmode test.

Run: `python3 test/run_tests.py --amps /path/to/amps --validation-case IV02
--output-dir test_output/IV02`.
