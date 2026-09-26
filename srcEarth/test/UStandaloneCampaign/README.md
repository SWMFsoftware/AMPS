# UStandaloneCampaign — Roadmap Step 8 validation

This dependency-free suite validates the standalone event/observation campaign
layer built on the unchanged Step-7 particle solver. It exercises orchestration,
data provenance, exposure averaging, observation comparison, restart, and
release evidence. It does **not** replace the existing trajectory, cutoff,
flux, spectrum, field, or observational physics tests.

Run from the repository root (or any directory):

```bash
./srcEarth/test/UStandaloneCampaign/run_test.sh
```

The script runs `test_campaign.py` and the independent
`test_step8_source_contract.py`. It requires Python 3 only; no AMPS executable,
MPI installation, network access, or mutable external data is used.

## Reference solutions and gates

The tests use four kinds of independent references rather than checking merely
that a file exists:

| Area | Reference or invariant | Failure condition |
|---|---|---|
| Cadence folding | The analytic mean of the linear series 1, 2, 3 over a ten-minute window is exactly 2 | Wrong integral, endpoint handling, duplicate cadence, nonfinite input, or any extrapolation |
| Tecplot extraction | A hand-written table has ratio 6/3 = 2, exact variables/order, and one row | Schema/order drift, row-count drift, missing column, zero denominator, nonfinite or inconsistent bounds |
| Observation metrics | Identical modeled/observed directional ratio gives log-RMSE 0 and correct sign fraction 1 | Missing comparison, unit mismatch, incomplete window coverage, or incorrect metric |
| Event execution | A deterministic fake Step-7 executable emits three epochs, numerical evidence, and canonical predictions | Missing artifacts, unresolved support above 0.01, energy/angular change above 0.02, or absent observation comparisons |

The adapter test also parses the committed observation references used by C9
(PAMELA), C10 (POES/MetOp), and C19 (GOES EPEAD). Those are real frozen tables,
not generated fixtures. A compact REPT row checks its released canonical schema.
The existing C9/C10/C19 scripts remain the authoritative full observation-facing
physics tests and their thresholds are unchanged.

## Negative tests

The suite deliberately proves that the following cannot pass:

- a resource whose SHA-256 changes, lacks provenance, or uses a URL;
- platform-specific normalization or an unregistered validation exclusion;
- an unresolved-support threshold above 0.01, or missing two-percent energy or
  angular convergence gates in a validation/holdout campaign;
- incompatible units, uncovered observation windows, invalid detector-response
  bounds, or a missing prediction/comparison;
- a previous PASS whose output artifact has been modified (that run is executed
  again while intact runs are restart-skipped);
- a dry run presented as scientific PASS; and
- incomplete O1/O2 evidence or an unfrozen/unverified O4 holdout manifest.

The O1/O2/O4 release-gate fixture has a complete PASS reference and then removes
an O2 comparison to confirm fail-closed behavior. No test edits existing C/F or
U-F thresholds, expected states, references, mover settings, or `last pass:`
records. The new test-list entry starts with an empty `last pass:` value until it
has actually passed in the target repository.
