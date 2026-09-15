# IV03 — manufactured spatial–pitch–momentum transport

IV03 defines a positive separable trigonometric/exponential distribution in
`s`, `mu`, `ln p`, and time with smooth advection, divergence, pitch diffusion,
and spatial diffusion coefficients. The exact manufactured source is the
analytic transport residual. The linked C++ evaluator reconstructs the same
residual from centered derivatives at three coordinate increments.

The saved table exposes the exact/numerical residual at every point rather than
only a final norm. Acceptance requires a second-order L2 trend, a bounded
finest-level error, and positivity. This implementation validates derivative
and source assembly for the available field-line transport terms; it is not
presented as a deterministic replacement for production Monte Carlo evolution.

Run: `python3 test/run_tests.py --amps /path/to/amps --validation-case IV03
--output-dir test_output/IV03`.
