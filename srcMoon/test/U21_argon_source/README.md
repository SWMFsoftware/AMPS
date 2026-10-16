# U21 — Radiogenic argon geography and transients

The complete executable contract—including purpose, production entry points,
configuration, units, frames/time, oracle, controls, procedure, metrics,
acceptance criteria, failure modes, data/hashes, predecessor gates, artifacts,
and status semantics—is in [reference/acceptance.json](reference/acceptance.json).

Run from the AMPS repository root:

```bash
python3 srcMoon/test/U21_argon_source/test.py
```

Exit codes are `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. A skipped
capability is intentionally not converted into a passing source-preflight result.

The runner now executes the compiled species-specific uniform-source kernel,
checks area-weighted rate conservation at the literature-selected
`3.8e20 atoms s^-1` lower bound, verifies the maintained input receipt, and
hashes the local D10 metadata/proxy. The final U21 result remains `SKIPPED`:
D10 covers only the South Pole-Aitken region and declares
`D10 DATA INCOMPLETE`, while U21 also requires potassium-weighted global and
transient linked injection tests. The passing nested `local_verification` in
`result.json` is kernel evidence only.
