# U15 — Multi-species photoionization tables

The complete executable contract—including purpose, production entry points,
configuration, units, frames/time, oracle, controls, procedure, metrics,
acceptance criteria, failure modes, data/hashes, predecessor gates, artifacts,
and status semantics—is in [reference/acceptance.json](reference/acceptance.json).

Run from the AMPS repository root:

```bash
python3 srcMoon/test/U15_multispecies_photoionization/test.py
```

Exit codes are `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. A skipped
capability is intentionally not converted into a passing source-preflight result.

The global runner now invokes supporting photon branch statistics, conservative
reaction primitives and runtime-input controls through this existing owner.
`result.json` preserves each child status, runner, result path and log under
`supporting_checks`; their counts are separate from the U15 gate. The actual
production registry is sampled over 12 parent/solar families, with 100000 draws
per family, frozen seed/stream construction and 99% familywise binomial bounds.
CSV/PNG/EPS comparisons and data/source hashes remain in the photon child.
See [the photon selection contract](../common/m3_rates/PHOTON_SELECTION_CONTRACT.md).

Use a new artifact directory without reconfiguring the application:

```sh
/usr/bin/python3 -B srcMoon/test/run_all_test.py --skip-build --only U15 \
  --output-dir test_output/srcMoon/U15-support-NEW
```

Passing these components does not resolve the original U25 predecessor or
qualify species/state-to-particle mapping, molecular recoil or linked I22.
The complete U15 remains SKIPPED until its unchanged full contract is met;
a child FAIL/ERROR is retained as that outcome, never hidden by a SKIP.
