# U17 — Charge-exchange kernels

The complete executable contract—including purpose, production entry points,
configuration, units, frames/time, oracle, controls, procedure, metrics,
acceptance criteria, failure modes, data/hashes, predecessor gates, artifacts,
and status semantics—is in [reference/acceptance.json](reference/acceptance.json).

Run from the AMPS repository root:

```bash
python3 srcMoon/test/U17_charge_exchange/test.py
```

Exit codes are `0=PASS`, `1=FAIL`, `2=ERROR`, and `77=SKIPPED`. A skipped
capability is intentionally not converted into a passing source-preflight result.

The registered owner now runs the maintained D04/CEX rate-registry controls as
supporting checks. They exercise actual production code against independently
evaluated published He coefficients and Na source knots, verify units and
strict energy domains, and retain the exact source PDF/data hashes. The
explicit CEX source directory comes from `examples/m3_reaction_rates.in`;
missing authoritative sources remain SKIPPED, not regenerated approximations.

```sh
/usr/bin/python3 -B srcMoon/test/run_all_test.py --skip-build --only U17 \
  --output-dir test_output/srcMoon/U17-support-NEW
```

Child runners, logs, results and counts are preserved in `supporting_checks`.
A passing total-cross-section/monoenergetic-rate service does not satisfy the
unchanged full U17 contract for differential/state-resolved products and
supported bath/projectile sampling. The full gate remains SKIPPED while those
capabilities are missing; evaluated child FAIL/ERROR remains visible. See
[the rate contract](../common/m3_rates/CONTRACT.md) and
[the data proposal](../M3_ION_SURFACE_AND_CEX_DATA_PROPOSAL.md).
