# U03 — Lunar gravity kernel and reference trajectories

The complete executable contract—including purpose, production entry points,
configuration, units, frames/time, oracle, controls, procedure, metrics,
acceptance criteria, failure modes, data/hashes, predecessor gates, artifacts,
and status semantics—is in [reference/acceptance.json](reference/acceptance.json).

Run from the AMPS repository root:

```bash
python3 srcMoon/test/U03_lunar_gravity/test.py
```

The probe calls `Moon::OrbitalDynamics::AddLunarPointMassAcceleration`, the
same production helper dispatched by `Moon::TotalParticleAcceleration`. It is
isolated deliberately so the closed-form lunar point-mass check remains valid
in both SPICE/orbit-on and no-SPICE regression builds; Sun/Earth and rotating
terms belong to U04 and are not silently included in the U03 oracle.

Exit codes are `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. A skipped
capability is intentionally not converted into a passing source-preflight result.
