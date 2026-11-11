# Comet validation and scientific decision gates

## Stage 10 scratch directories in the integrated AMPS tree

Run the focused receipt and merge-readiness tests from the AMPS root:

```sh
bash srcHartley2/tests/run_stage10_tests.sh
```

Stage 10 creates unique `stage10-readiness-selftest-*` and
`stage10-all-stage-result-*` temporary subdirectories directly beneath the
AMPS root, located relative to the test source rather than the current working
directory. Neither test uses `/tmp`, `TMPDIR`, or a platform fallback for these
fixtures. The AMPS root must be writable; an unavailable root produces an
error rather than redirecting the fixture elsewhere. Concurrent invocations
receive independent directories.

The readiness test initializes a separate Git repository at
`AMPS/stage10-readiness-selftest-*/repository/.git` and makes synthetic commits
there to exercise ancestry, changed-file classification, and rejected evidence.
These commits contain test inputs, not test results. The test validates every
fixture Git target, clears inherited Git path/configuration overrides, and
disables external hooks and commit signing. Git state in the enclosing AMPS
repository and in the source/reference checkout remains read-only.

Both scratch directories are automatically removed on completion or ordinary
Python exceptions. An uncatchable termination such as `SIGKILL` or a machine
crash can leave a directory behind; inspect the exact printed scratch path
before removing that individual directory. No test removes an existing
directory or changes Git ignore rules. Durable suite logs/JSON remain the
responsibility of the outer runner and are separate from this scratch data.

The application, inputs and shared comet libraries have now been imported.
The portable all-stage runner passes 45 jobs across 14 groups; its archived
runtime arms recheck frozen evidence rather than executing a new native
simulation. The subsequent all-stage runner passes 46 jobs across 14 groups,
with fresh Stage-8/9 runs at 1, 4, and 12 MPI ranks. Earlier fresh Stage-4–7
matrices also pass and remain in their own evidence directories. Fresh I3,
canonical Halley, and radial-mesh checks pass. The earlier complete-core campaign
preserved all 785 baseline passing cases and added six passes, with its original
two failures, one timeout error and 240 skips unchanged. A later user run
(`run-20261110T092615Z-924820`) exposed two additional configuration-dependent
failures in native exosphere/dust fixtures inheriting global SPICE selection.
Those runners are now repaired and checked with independent starting states.
The subsequent full 1,035-case core campaign preserves all 785 original passes,
restores both failing cases, and adds one harness pass (792 PASS, original
2 FAIL / 240 SKIP / 1 ERROR unchanged). The earlier final-tree legacy
Moon/fast-wave comparisons have zero-byte diffs. The two additional 3-D build
exceptions were approved and applied on 2026-11-10; the clean native build
passes with both site-default and optimized compiler flags. The optimized
four-rank nightly run reached the local 1,800-second smoke limit during
initial cut-cell volume construction, before any particle iteration. That
runtime remains unvalidated. The subsequent read-only performance investigation
identifies repeated full-surface ray searches; see
[the diagnosis and proposed validation plan](../../src/test/GEOMETRY_PERFORMANCE_INVESTIGATION.md).
No production geometry optimization is included in these build fixes. See
[the integration ledger](../../COMET_CORE_INTEGRATION_STATUS.md) for exact
commands, frozen binaries, approved compatibility exceptions and evidence.
These implementation checks do not establish long-run science acceptance.

Focused verification on 2026-11-09 passed the eight receipt-writer cases and
the readiness mutation suite with `TMPDIR`, `TEMP`, and `TMP` pointing to
nonexistent directories. Two concurrent readiness runs, launched from a
different working directory with invalid inherited `GIT_DIR`, `GIT_WORK_TREE`,
and `GIT_INDEX_FILE`, also passed and removed their distinct AMPS-local scratch
directories. SHA-256 inventories of every enclosing AMPS `.git` file were
identical before and after execution. These checks exercise the new location,
cleanup, concurrency, and Git-environment isolation rather than native AMPS
physics or completion of the full integration.

## Stage 3 scientific decision gate

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

The canonical case is `input/halley1d.input`. Run these commands from the AMPS
repository root; all paths are relative to it. Remove exactly the root
`build/` directory before compiling:

```sh
rm -rf -- build
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

## Stage-9 balanced analog/enhanced campaign

`run_stage9_science_campaign.py` manages immutable balanced runs from one
frozen executable and physical input. `plan` authenticates and copies the
executable, input, `picGlobal.dfn`, generated chemistry table, and generated
Stage-6/7/8 policies before any launch. `run` accepts only a pristine plan;
`resume` retains already accepted cases but creates a fresh attempt for an
interrupted or failed case instead of reusing partial `PT` state. `status` is
read-only.

Each completed case is reduced by `stage9_run_summary.py`. The reducer reads
the AMPS convergence and Stage-0 chemistry reports, Stage-9 checkpoint JSON/
CSV and per-cell carrier tables, profiles, and the transported-collision log.
It refuses missing, ambiguous, altered, nonfinite, or internally inconsistent
evidence; it never fills a missing science value by interpolation or a default
assumption.

Planning requires a resolved copy of
`stage9_science_policy.template.json`. The distributed template deliberately
uses `DECISION_REQUIRED` for the minimum daughter carrier-ESS fraction,
variance-per-walltime improvement, required resolved-bin fraction, and flux-
closure choice. Flux closure must be explicitly either
`independent_seed_ci_contains_zero` with a documented null tolerance or
`relative_tolerance` with a positive numeric tolerance. Synthetic fixture
values are not production choices.

Run the dependency-free mutation suites with:

```bash
srcHartley2/tests/run_stage9_science_gate_tests.sh
srcHartley2/validation/run_stage9_tooling_tests.sh
```

The first currently covers 15 gate mutations and the second covers 9 archive,
summary, recovery, and publication cases. Their PASS validates tooling only;
the Halley claim requires a complete checksum-protected manifest with at least
four independent seeds per arm.

## Stage-9 ensemble Rubin comparison products

After the complete balanced Stage-9 campaign passes the strict science gate,
publish the independent-seed ensemble comparison with:

```bash
python3 srcHartley2/validation/generate_stage9_rubin_comparison.py \
  validation_runs/STAGE9_CAMPAIGN/RUN_MANIFEST.json \
  --output-directory validation_runs/STAGE9_CAMPAIGN/rubin2011_comparison
```

The exporter calls the checksum-protected science evaluator itself. It refuses
an incomplete or altered summary, differing physical provenance, a failed
convergence/closure/ESS/variance gate, and any attempt to overwrite an
existing product. The five products share the basename
`stage9_rubin2011_ensemble_comparison` and have `png`, `eps`, `pdf`, `csv`, and
`json` extensions. Model density is converted from m^-3 to cm^-3. Every model
curve is the mean over the documented independent seeds; dashed curves are
the confidence bounds calculated with the policy-supplied multiplier. The
committed Rubin vector digitization is overlaid only for species present in
the Stage-9 gate policy, and every available published variant is identified
explicitly.

The JSON records seed identities, summary and manifest checksums, physical
input identity, convergence state, complete post-burn-in window exposure,
wall/MPI cost, paper-data and source-PDF provenance, and the exact plotted
values. The long-form CSV distinguishes `ANALOG`, `FORCED`, and `RUBIN2011`
rows, so the comparison can be reproduced without reading a raster plot.
Synthetic tooling tests are run with:

```bash
srcHartley2/validation/run_stage9_comparison_tests.sh
```

Those tests validate authentication, fail-closed behavior, unit conversion,
and output formats. Their synthetic PASS is not comet science evidence.

## Stage-10 focused-test receipt and merge-readiness bundle

Stage 10 separates four facts that must not be collapsed into one PASS:

1. focused Stage 0--10 implementation tests completed;
2. Stage-9 long-run science evidence was supplied and accepted;
3. the resulting exact diff is complete enough for user review; and
4. the user explicitly authorized an operation affecting `master`.

Only the first three can be evaluated here. The fourth is always serialized
as `NOT_GRANTED`; none of these programs runs a mutating Git command.

Run the focused suite with an explicit new log root. In addition to
`summary.txt`, it writes `all_stage_result.json` with schema
`amps.comet.all-stage-test-result.v1`:

```bash
srcHartley2/tests/run_all_stage_tests.sh \
  --log-root /new/absolute/path/to/stage10-all-stage-tests
```

The receipt contains the exact HEAD commit/tree, branch state, every test and
group result, log checksum, aggregate counts, archived/live mode, and
interpretation. A receipt created without a real accepted Stage-9 manifest
must explicitly say focused implementation passed but science did not.

The policy and tool are:

```text
stage10_release_policy.json
stage10_merge_readiness.py
```

Use `stage10_merge_readiness.py --help` and the subcommand help as the
authoritative CLI. `collect` requires an explicit policy-listed target ref,
full target/base/development commits, the versioned policy, the focused-suite
receipt and its SHA-256, and a brand-new output directory. It optionally
authenticates real Stage-9 campaign/comparison evidence. It does not infer a
target from the spelling `master`. `verify --bundle DIRECTORY` recomputes
every bundle checksum and validates the recorded status and authorization.

The report schema is `amps.comet.stage10-merge-readiness.v1`. Without accepted
real science it reports `FOCUSED_ONLY_NOT_SCIENCE_READY`; with both gates it
may report `READY_FOR_USER_REVIEW`. Both results retain
`merge_authorization=NOT_GRANTED`. Dirty or untracked reviewable source,
unclassified diff paths, altered receipts/logs, missing evidence, commit/tree
mismatches, reused output, and checksum changes fail closed. Explicit
generated/runtime untracked exclusions are summarized by policy rule; a
custom policy cannot weaken the canonical classifications.

Run the Stage-10 mutation tests with:

```bash
srcHartley2/tests/run_stage10_tests.sh
```

The test fixtures use temporary repositories/evidence only. They prove the
tooling contracts and cannot constitute comet science or permission to merge.
The manual semantic integration order, core/application ownership boundaries,
exclusions, and post-integration matrix are documented in
`../MASTER_MERGE.md`.
