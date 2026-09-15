# CV11 — time-dependent resonant wave growth and damping

CV11 isolates the production turbulence source ledger on a one-cell, one-hot
spectral state. For prescribed net rate `Gamma-gamma`, the independent answer
is `W(t)=W0 exp(2 integral[Gamma-gamma]dt)`. Constant growth, constant damping,
exact cancellation, and a sign-changing sinusoidal history cover positive,
negative, zero, and time-dependent rates.

Only the registered resonant bin is energized. Constant-rate increments are
constructed exactly; the sinusoidal rate uses midpoint quadrature so three
time steps exhibit a genuine second-order convergence trend against the exact
integral. Acceptance checks active-bin histories, temporal order, inactive-bin
contamination, cancellation, ledger closure, and positivity. This controlled
rate test does not validate a particular physical growth-rate closure; it
validates how a supplied rate changes the production wave state.

Run `python3 test/run_tests.py --amps /path/to/amps --validation-case CV11
--output-dir test_output/CV11`.
