# IV04 — moving field-line grid conservation and GCL

IV04 initializes a uniform cell-integrated state and calls the production
conservative overlap remap for rigid translation, stretched, compressed, and
sinusoidal segment distributions at 16, 32, and 64 cells. In arclength
coordinates a rigid translation leaves every cell measure unchanged, which is
the discrete geometric-conservation-law free-stream control. Segment lengths are normalized
to the same physical domain and validated as positive before scoring.

For a uniform physical density, the exact target cell integral equals its new
segment length. Therefore any deviation is a direct free-stream/geometric-
conservation failure. The case also checks total integral, production remap
ledger residual, and minimum segment length. Raw node positions and mapped
profiles support diagnosis of phase, coverage, or tangled-cell errors.

Run: `python3 test/run_tests.py --amps /path/to/amps --validation-case IV04
--output-dir test_output/IV04`.
