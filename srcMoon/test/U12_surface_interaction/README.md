# U12 — Sticking, temperature, residence, and re-emission

The complete executable contract—including purpose, production entry points,
configuration, units, frames/time, oracle, controls, procedure, metrics,
acceptance criteria, failure modes, data/hashes, predecessor gates, artifacts,
and status semantics—is in [reference/acceptance.json](reference/acceptance.json).

Run from the AMPS repository root:

```bash
python3 srcMoon/test/U12_surface_interaction/test.py
```

Exit codes are `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. A skipped
capability is intentionally not converted into a passing source-preflight result.

The accommodation probe evaluates both supported input interpretations at
different incident and thermal speeds. In particular, it independently
computes the kinetic-energy identity in Tenishev et al. (2013), Eq. 14, and
checks its zero- and unity-accommodation limits. This deterministic U-test does
not replace I18's linked stochastic angular/speed-distribution and inventory
budget gates.
