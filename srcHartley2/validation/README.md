# Stage 3 scientific decision gate

This directory contains the permanent decision gate for determining whether
the Stage 3 systematic-residual scheme gives adequate minor-daughter
statistics without implementing individual particle corrections in Stage 4.
It separates three operations that must not be conflated:

- `stage3_science_policy.json` freezes the approved scientific thresholds;
- `run_stage3_science_gate.py` creates and resumes immutable AMPS evidence;
- `evaluate_stage3_science_gate.py` reads a complete archive and assigns the
  scientific outcome.

The gate uses one frozen executable for both residual modes. The executable's
`--chemistry-residual-sampling bernoulli|systematic` option is applied before
application block storage is requested, so Bernoulli runs reserve no
systematic residual state while both arms retain identical compiled physics,
chemistry, mesh, and population policies.

## Approved experiment

The policy is intentionally validated against constants in
`stage3_gate_common.py`; replacing the JSON with a newly hashed, relaxed
policy cannot silently weaken the gate. The approved matrix is:

- four fixed independent seeds in both Bernoulli and Systematic modes;
- 12 MPI processes and one OpenMP thread per run;
- 10,000 iterations of burn-in followed by three complete 10,000-iteration
  accepted windows, with a hard limit of 40,000 iterations;
- every accepted window must have a maximum species-profile L1 change of at
  most 5%, valid reaction/species production closure, direct-source agreement
  within 1%, and no emergency request or executed emergency adaptation;
- exact fully contained cells from 10 km through 100,000 km altitude above
  the 6 km nucleus;
- common adaptive-bin support chosen by the minimum exposure over all eight
  runs, with at least 100,000 particle-iterations per bin;
- C2, C2H, CN, C2H5, and HCO evaluated with the unlocked `CHEMICAL_ONLY`
  50--200 particles-per-cell block policy;
- Systematic CV at most 25% in at least 90% of qualified bins, no qualified
  bin above 50% CV, raw zero-cell fraction at most 5%, no empty adaptive bin,
  and at least 90% of populated complete blocks inside the population target;
- Bernoulli/Systematic means agreeing within the larger of 5% of their scale
  and two combined standard errors; and
- no more than 10% regression in the Systematic ensemble's log-density RMSE
  relative to an explicitly accepted Rubin baseline.

The first statistically valid evaluation writes
`evaluation/rubin_candidate_baseline.json` with `accepted=false`. Review that
file and the comparison products before copying it to a durable source file
and changing only `accepted` to `true`. Baseline acceptance is a scientific
review decision; neither runner nor evaluator promotes it automatically.

## Build and preflight

The canonical case is `input/halley1d.input`. Remove exactly the repository
`build/` directory before compiling:

```sh
cd /data/vtenishe/COMET1D/AMPS
rm -rf /data/vtenishe/COMET1D/AMPS/build
make -j8 EXTRACOMPILEROPTIONS=-O3
python3 srcHartley2/validation/run_stage3_science_gate.py plan \
  validation_runs/stage3_science_gate_YYYYMMDD \
  --executable amps \
  --input input/halley1d.input \
  --pic-global build/pic/picGlobal.dfn
```

`plan` is read-only and prints the exact hashes, environment, ranks, and
alternating order. Inspect that output before launching.

## Run, recover, and evaluate

Start a new empty archive only once:

```sh
python3 srcHartley2/validation/run_stage3_science_gate.py run \
  validation_runs/stage3_science_gate_YYYYMMDD \
  --executable amps \
  --input input/halley1d.input \
  --pic-global build/pic/picGlobal.dfn
```

The runner executes cases sequentially so they do not contend for the same
host memory. It atomically records an attempt before launch, stores stdout and
stderr in that attempt, validates all four complete diagnostic windows, and
hashes every required artifact. It never uses AMPS particle restart data.
After a reboot or signal, inspect and resume with:

```sh
python3 srcHartley2/validation/run_stage3_science_gate.py status \
  validation_runs/stage3_science_gate_YYYYMMDD
python3 srcHartley2/validation/run_stage3_science_gate.py resume \
  validation_runs/stage3_science_gate_YYYYMMDD
```

Completed cases are re-hashed and skipped. A failed or interrupted case gets
a new numbered attempt starting at iteration zero; partial `PT` output is
preserved for diagnosis but never reused as science evidence.

After all eight cases are operationally complete:

```sh
python3 srcHartley2/validation/evaluate_stage3_science_gate.py \
  validation_runs/stage3_science_gate_YYYYMMDD \
  --rubin-baseline path/to/accepted_stage3_rubin_baseline.json
```

Before reading any science value, the evaluator independently re-hashes the
frozen executable/input/`picGlobal.dfn` and every file listed by each
`artifacts.sha256`; post-run edits or damage are `INVALID_EVIDENCE`. It then
writes JSON, CSV, and Markdown under `evaluation/`. Its outcome and exit
status are:

| Outcome | Exit | Meaning |
| --- | ---: | --- |
| `PASS_SKIP_STAGE4` | 0 | statistical gates pass with no detected bias; retain Stage 3 |
| `FAIL_CONSIDER_STAGE4` | 10 | unbiased evidence is statistically insufficient; evaluate Stage 4 |
| `INCONCLUSIVE_RERUN` | 11 | convergence, emergency, closure, mode-agreement, or baseline gate is unresolved |
| `INVALID_EVIDENCE` | 12 | archive, provenance, geometry, or schema is invalid |

The operational runner separately uses exit 0 for success, 2 for invalid
policy/input, 3 for provenance/archive preflight failure, 4 for AMPS or
artifact failure, and 5 for interruption. Operational completion is not a
scientific pass.

## Focused tests

These tests do not launch the long eight-run ensemble:

```sh
srcHartley2/tests/run_stage3_accumulator_tests.sh
srcHartley2/tests/run_stage3_parser_tests.sh
srcHartley2/tests/run_stage3_science_gate_tests.sh
PYTHONDONTWRITEBYTECODE=1 python3 srcHartley2/tests/test_stage3_adaptive_bins.py
PYTHONDONTWRITEBYTECODE=1 python3 srcHartley2/tests/test_stage3_diagnostic_semantics.py
```

Do not merge to `master` from this workflow. Retain the complete archive,
accepted baseline, exact hashes, outcome summary, performance/memory record,
and comparison plots for the later user-approved integration review.
