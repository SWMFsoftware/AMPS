# U07 — Surface-refinement selector unit verification

The complete executable contract—including purpose, production entry points,
configuration, units, frames/time, oracle, controls, procedure, metrics,
acceptance criteria, failure modes, data/hashes, predecessor gates, artifacts,
and status semantics—is in [reference/acceptance.json](reference/acceptance.json).

The probe preserves the historical no-input callback check and also invokes
the production M1 configuration kernels directly. It covers the six approved
SZA knots, an interpolation point, the terminator and body-fixed polar minimum
selectors, and the radial near/far endpoints. Linked AMR realization is I12;
observable convergence is I11. The configured radial checks cover the 10 km
near-surface target through 200 km, the 40 km science-region endpoint at
10,000 km altitude, and the separate 17,371 km target at the legacy
400-lunar-radius boundary. This prevents the near-Moon qualification scale
from being silently imposed throughout the tail domain.

The same probe evaluates the committed QUICK targets—20 km near-surface,
80 km at the science-region endpoint, 34,742 km at the outer boundary, and
10 km terminator/polar refinement—while retaining the FULL physical extents.
This verifies selector arithmetic only; the approximately 20-minute runtime
estimate and realized linked QUICK mesh still require an I12 wiring run.

Run from the AMPS repository root:

```bash
python3 srcMoon/test/U07_surface_refinement/test.py
```

Exit codes are `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. A skipped
capability is intentionally not converted into a passing source-preflight result.
