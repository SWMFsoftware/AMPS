# U11 — Thermal inertia and eclipse response

The complete executable contract—including purpose, production entry points,
configuration, units, frames/time, oracle, controls, procedure, metrics,
acceptance criteria, failure modes, data/hashes, predecessor gates, artifacts,
and status semantics—is in [reference/acceptance.json](reference/acceptance.json).

Run from the AMPS repository root:

```bash
python3 srcMoon/test/U11_thermal_inertia/test.py
```

Exit codes are `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. A skipped
capability is intentionally not converted into a passing source-preflight result.

The November 12 radiative-equilibrium case initializes every depth at 250 K,
applies exactly balancing absorbed and emitted surface flux, and advances the
same production column solver at 1, 60 and 3600 s. Uniform temperature and
zero conductive gradient are an analytical steady solution; each result must
remain within 1e-8 K. This tests a thermal invariant, without observational
inputs, spin-up assumptions or a claim of transient time-step convergence.
