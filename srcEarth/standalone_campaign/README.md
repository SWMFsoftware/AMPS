# Standalone event campaigns — Roadmap Step 8

This directory turns the Step-7 standalone cutoff/flux/spectrum executable into a
reproducible event workflow. It does not contain a second particle model. Every
trajectory, access value, spectrum, density, flux, and detector rate still comes from
the production Step 2–7 C++ kernels.

Existing gates remain unchanged. Step 8 adds an independent fail-closed campaign
gate; it does not modify an existing P/F command, tolerance, reference, mover,
trace limit, expected result, or `last pass:` value.

## Physics scope

The released Phase-1 workflow evaluates instantaneous or quasi-static magnetic-field
snapshots. The field is frozen during each backward-trajectory batch and the boundary
distribution, ephemeris, attitude, response, and output are associated with the same
authoritative epoch. The campaign layer performs exposure-window averaging after the
snapshot calculations. It does not imply time-dependent particle acceleration,
long-duration trapping/loss physics, or a time-dependent electric characteristic.

The comparison policy is `SHARED_EVENT_BOUNDARY_NO_PLATFORM_SCALE`: every platform in
an event uses the same independently established boundary population. The runner
rejects a manifest that requests platform-specific multiplicative normalization.

## Files

- `campaign.py` validates the manifest, verifies every SHA-256 digest, creates a
  content-addressed local cache, renders strict input templates, parses termination
  evidence, evaluates declared numerical gates, and verifies restart artifacts.
- `adapters.py` normalizes PAMELA cutoff intervals, POES/MetOp MEPED cutoff
  boundaries, GOES EPEAD physical east/west ratios, and Van Allen Probes REPT proton
  spectra. It also performs exact trapezoidal averaging of the declared piecewise
  linear model series over each observation accumulation window.
- `run_campaign.py` creates one directory per field model and epoch, streams AMPS
  output to the console and `amps.log`, records progress, resumes only hash-verified
  PASS runs, extracts declared prediction columns/ratios from AMPS products, and
  writes machine-readable numerical and observation summaries.
- `release_gate.py` consumes the completed O1 and O2 summaries and preflights the
  untouched O4 holdout manifest. It reports PASS only when all required runs,
  numerical gates, and comparisons are present and successful.

All Python code uses the standard library and is compatible with direct execution
from an AMPS checkout.

## Manifest contract

The schema identifier is `earth-standalone-campaign/v1`. A valid manifest must state:

| Section | Required content |
|---|---|
| Identity/event | Portable campaign ID, explicit event ID, UTC start/end, and field cadence; the interval must contain an integer number of cadences. |
| Physics | Field models and requested cutoff/access/spectrum/density/flux/detector products. |
| References | Input template, driver, boundary spectrum, ephemeris, attitude, instrument response, and observations. Every local file has a provenance statement and fixed SHA-256 digest. URLs are forbidden during scoring. |
| Instruments | Adapter, observation and response resource IDs, physical quantity, exact units, cadence, and whether comparison is required. |
| Execution | Tokenized command, MPI/thread counts, expected artifacts, canonical prediction file, and Tecplot extractor with exact variable order, row count, operation, identity, and units. |
| Gates | Every campaign retains unresolved support `<= 0.01`; validation/holdout campaigns also retain energy and angular relative change `<= 0.02`. Missing or weaker evidence is rejected before execution. |
| Exclusions | An explicit list, including an empty list when nothing is excluded. Validation/holdout exclusions must have a reason and be preregistered. |
| Blinding | Campaign role and a literal `frozen` flag. `VALIDATION` and `HOLDOUT` manifests cannot be unfrozen. |

Resource acquisition is deliberately outside the runner. Download and document an
official product once, record its checksum and calibration/quality policy, then run
the campaign only against that immutable local file. A changed byte blocks preflight.

## Prediction extraction

AMPS Step 7 writes Tecplot products. `prediction_extractors` converts reviewed output
columns into the canonical comparison table. `COLUMN` extracts one flux, spectrum, or
rate column. `RATIO` divides two declared columns, for example physical east and west
detector rates. A glob must match exactly one artifact; a missing/ambiguous file,
unknown/reordered variable, changed row count, invalid row, zero denominator,
nonfinite value, or inconsistent uncertainty bounds fails the run.

The canonical prediction fields are:

```text
utc,instrument_id,platform,channel,direction,quantity,units,value,lower,upper
```

An externally produced canonical file is also accepted, but a required comparison
cannot pass when the file or matching rows are absent. This is intentional: a run that
does not produce a comparison is a failed validation run, not an unscored pass.
Units are part of the comparison identity, so incompatible unit strings never match
implicitly.

## Cadence and observation gates

The model series is linearly interpolated between frozen snapshots and integrated over
the exact observation start/end interval with the trapezoidal rule. The routine never
extrapolates. A window without complete model coverage is counted as missing.

Applicable roadmap gates are evaluated without post-hoc adjustment:

- spectrum/count amplitude: median absolute log10 ratio no larger than 0.30;
- per-instrument factor-of-two success fraction at least 60% and pooled event
  factor-of-two fraction at least 70%;
- LEO cutoff latitude: MAE no larger than 2 degrees AACGM and absolute bias no larger
  than 1 degree;
- directional ratio: log10 RMSE no larger than 0.30 and correct asymmetry sign in at
  least 80% of quality-controlled samples; and
- when model bounds are provided, 90% interval coverage between 80% and 98%.

Every quality-controlled observation remaining after the frozen adapter policy must
have model coverage (`comparison_fraction = 1`). A missing window is not removed from
the denominator and cannot be hidden as a successful partial comparison.

Energy/angular, mover/time-step, mesh, backend, and parallel-convergence evidence is
not fabricated by this runner. Production manifests declare the corresponding frozen
JSON evidence gates. Missing evidence fails.

## Running

Preflight every file, adapter, cadence, and manifest rule without AMPS:

```bash
python3 srcEarth/standalone_campaign/run_campaign.py \
  --manifest srcEarth/examples/standalone_step8_campaign/campaign.json \
  --output-dir test_output/standalone_step8 \
  --validate-only
```

Render every input and command without launching AMPS:

```bash
python3 srcEarth/standalone_campaign/run_campaign.py \
  --manifest srcEarth/examples/standalone_step8_campaign/campaign.json \
  --amps ./amps \
  --output-dir test_output/standalone_step8 \
  --dry-run
```

Execute or safely resume:

```bash
python3 srcEarth/standalone_campaign/run_campaign.py \
  --manifest srcEarth/examples/standalone_step8_campaign/campaign.json \
  --amps ./amps \
  --output-dir test_output/standalone_step8 \
  --restart
```

`--restart` skips a run only when its prior status is PASS, its manifest/input/
executable fingerprint is unchanged, and every recorded artifact still matches its
digest. Failed, incomplete, modified, or stale runs execute again.

After the full O1/O2 campaigns and the frozen O4 holdout package exist, aggregate the
Step-8 release decision without rerunning or reinterpreting a campaign:

```bash
python3 srcEarth/standalone_campaign/release_gate.py \
  --o1-summary results/O1/campaign_summary.json \
  --o2-summary results/O2/campaign_summary.json \
  --holdout-manifest campaigns/O4/campaign.json \
  --output results/step8_release_gate.json
```

The release command exits 2 on incomplete/failed evidence. O1 and O2 must be frozen
`VALIDATION` campaigns with real executed PASS runs, successful required numerical
and observation gates, and the shared no-platform-scale policy. O4 must be an
untouched frozen `HOLDOUT` manifest whose local resources all match their hashes.

## Result tree

```text
campaign_preflight.json
campaign_summary.json
progress.jsonl
termination_summary.json
convergence_summary.json
observation_summary.json
campaign_predictions.csv
cache/<sha256>/<reference-file>
observations/<instrument-id>.csv
runs/<model>/<epoch>/amps.in
runs/<model>/<epoch>/amps.log
runs/<model>/<epoch>/run_status.json
comparisons/<model>/<instrument-id>/comparison.csv
comparisons/<model>/<instrument-id>/summary.json
```

`campaign_summary.json` is the top-level status. A required AMPS run, artifact,
numerical gate, or observation comparison failure makes it `FAIL` and the runner exits
with code 2. `--dry-run` reports `NOT_RUN`, never PASS.

## O1/O2 release use

- O1 (13–15 December 2006) uses the existing C9 PAMELA and C10 archive-derived
  POES/MetOp assets, plus independently documented boundary, response, ephemeris, and
  attitude files in its frozen manifest.
- O2 (17 May 2012 GLE71) uses the existing C19 GOES-13/15 physical east/west mapping,
  response provenance, event geometry, and unresolved-support diagnostics.

The code is ready to orchestrate these events, but the small example is a workflow
smoke case, not evidence that the full O1/O2 release ensembles passed. Those claims
require the complete linked AMPS executable, production resources, numerical
convergence ensembles, and archived comparison summaries.

## Tests

From the AMPS repository root:

```bash
./srcEarth/test/UStandaloneCampaign/run_test.sh
```

The suite uses analytic exposure averages, hand-derived Tecplot ratios, the committed
C9/C10/C19 observation assets, and a deterministic fake Step-7 executable. It checks
correct PASS behavior, immutable hash/provenance rejection, no-network scoring,
no-platform-normalization policy, all four adapters, exact units and schemas,
hash-verified restart/re-execution, no extrapolation, required comparison coverage,
dry-run status, O1/O2/O4 aggregation, and negative numerical cases that must fail at
the unchanged 0.01 unresolved and 0.02 energy/angular gates. See
`test/UStandaloneCampaign/README.md` for the test-by-test reference specification.
