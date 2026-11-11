# Full-volume Corona–SWCME global test owner

Updated 2026-11-10. `run_test.py` owns
`CORE-MODEL-SEP-CORONA-SWCME-FULL-001`. It imports the exact 46-row maintained
software registry; it does not duplicate the model's physical oracles or turn
a software PASS into BG3D-4/full-CME qualification.

[Global procedure](../../../../README.md) ·
[All CME cases](../README.md) ·
[Model APIs](../../../../../models/sep_corona_swcme/README_CROSS_RAY.md) ·
[Local binary commands](../../../../../models/sep_corona_swcme/test/README.md).

## Commands and selection

From the AMPS root, select the complete owner:

```sh
python3 -B src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --test CORE-MODEL-SEP-CORONA-SWCME-FULL-001 \
  --output-dir test_output/cswc-full-check --build-jobs 16
```

Use a fresh output directory. To validate all registered CORE suites:

```sh
python3 -B src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --match 'CORE-*' --jobs 8 --build-jobs 16 \
  --output-dir test_output/core-validation
```

Use `--match 'CORE-*'` without `--all` for core-only selection; `--all` selects
the complete registry, including applications. `--config` is mandatory. On
another machine, start from `test/qualification-machine/template.json` and
set real tools, paths and capabilities. No compiler/kernel/data path is baked
into this adapter.

A report subset can be requested explicitly:

```sh
python3 -B src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --case CSWC-FULL-CSWC0628-XRAY-MATERIAL-NONLINEAR \
  --output-dir test_output/cswc-material-case
```

This owner's case filter applies AFTER complete local-registry execution and
validation. It does not make the selected binary a quick one-case run. For
quick development, build/run the exact standalone target in the local test
README; the global owner remains the full import/registration check.

## Actual execution and concurrency

1. Resolve the supplied machine context and exact requested case registry.
2. Copy the real shared-model source/dependency trees and required thin
   adapters into `work_dir/isolated-repo`, excluding generated build objects,
   caches and test output. Fixtures and parser headers are copied from the
   actual source tree; no mock physics replaces an include dependency.
3. Build the 23 registered C++ executables through the local Makefile with
   configured parent make parallelism. A shared coronal archive is built by
   its maintained dependency recipe.
4. Invoke each binary sequentially in Makefile order. Print `[CASE] RUN` with
   each owned GLOBAL case ID before starting that binary. A binary with
   several cases still runs as one process; its RUN lines are not separate
   concurrent jobs or individual timing receipts.
5. Execute the maintained public-header/archive architecture audit.
6. Reconcile the complete ordered globalized registry, then apply an optional
   selected-case subset and write the typed result.

The current manifest marks this owner `exclusive=true`: it runs without other
core suites in the same stream. `--jobs 8` controls eligible nonexclusive
suites, not the binaries within this owner. `--streams N` independently bounds
registered core/application streams. Native checkout-mutating builds across
streams remain serialized by the global procedure. These model-library builds
use private trees and are not native AMPS builds or G1 native-run evidence.

## Portable make/compiler settings

The adapter consumes `tools.make` and `tools.python` as argv vectors. The
model Makefiles default to g++; a model compiler override can be carried by a
site make wrapper or make assignment in a layered overlay, for example:

```json
{
  "schema_version": 1,
  "tools": {
    "make": ["/usr/bin/make", "CXX=/opt/toolchain/bin/g++"]
  },
  "paths": {
    "output_root": "test_output/core-results",
    "temporary_root": "test_output/core-work"
  },
  "execution": {"build_jobs": 16}
}
```

Pass the overlay with a second `--config` after the base machine profile.
Relative paths resolve through the global configuration mechanism. Setting
`tools.cxx` alone does not override the CXX variable inside these maintained
Makefiles. The architecture script currently invokes `g++` and `nm` from PATH,
and the coronal dependency recipe explicitly requests `-j16`; those are
current portability boundaries, not a claim of arbitrary-toolchain
qualification. A site must provide them or retain the resulting setup ERROR.

## Budgets and truthful failure import

| Operation | Current ceiling / result |
| --- | --- |
| Exact-registry make invocation | 3600 s; build failure is owner-level ERROR |
| Ordinary local binary | 600 s |
| `test_bg3d4_piston_refinement` | 180 s; timeout is ERROR for its owned ID, with no inferred PASS from partial levels |
| Architecture audit | 300 s |
| Whole owner | Manifest timeout 3600 s |

Parent machine timeout settings do not silently change the hard-coded local
budgets. The refinement's initial cells, all supported rays, 3600-s epoch,
CFL/refinement levels, single-branch condition and physical/convergence limits
are unchanged. Increasing a diagnostic execution budget is not qualification.

Import rules deliberately fail closed:

- Exit 0 must publish every owned PASS marker in exact order, with no omitted
  or foreign marker.
- A typed nonzero result must retain the complete owned registry and agree
  with status/exit mapping: PASS/SKIP=0, FAIL=1, ERROR=2, using the maximum
  status exit for mixed rows.
- A narrow maintained legacy form, exit 1 plus explicit recognized FAIL
  prose, becomes FAIL for that binary's owned identities. A signal, unknown
  nonzero result or unresolved setup problem remains ERROR.
- A handled per-binary timeout preserves that binary's ERROR rows and allows
  later binaries to execute.
- Missing/duplicate/reordered/partial markers are owner-level protocol ERROR.
  Requested IDs remain visible as ERROR rather than disappearing. In
  particular, a multi-case binary aborting after partial output has not
  established the remaining assertions; those are not fabricated PASS rows.

The combined transcript is written after complete owner import. On an early
build/protocol abort it may be absent/incomplete; use the global owner error
message/log and retained work directory. This is a current diagnostic-retention
limitation, not evidence that all physical cases executed.

## Registry and contributor contract

`LOCAL_IDS`, `BINARY_CASES`, the model Makefile and central
`case_registries.json` entry `T11-CSWC-FULL` must agree. The original owner has
15 rows; new work adds 31 model rows in ten binaries (29 transport/material/
contact cases plus contact-derivative and shell controls), giving 46 including
`CSWC-FULL-ARCHCSWC01`. The separate coverage-governance case adds the 32nd new
global case in the harness.

A contributor/agent adding a case must update its independent oracle, local
binary/markers, exact adapter/central registry, ownership, local/global README
and generated catalog in the same change. Run:

```sh
python3 -B src/test/generate_test_catalog.py
python3 -B src/test/generate_test_catalog.py --check
python3 -B src/test/validate_traceability.py
```

Catalog and ownership consistency do not supply missing physical evidence.
Never delete an expensive/failing identity, rewrite a fixture to the calculated
answer, downgrade a required gate or upgrade partial coverage to make green.

## Reports and current evidence

The global output directory contains:

```text
amps-core-tests.json                    complete combined machine-readable result
amps-core-tests.xml                     JUnit result
streams/core.log                        suite progress and stream summary
streams/core.json                       core stream report
logs/CORE-MODEL-SEP-CORONA-SWCME-FULL-001.log
artifacts/CORE-MODEL-SEP-CORONA-SWCME-FULL-001/full-volume-local.txt
```

The JSON suite row records `log`, `artifact_dir`, `work_dir`, evidence paths
and elapsed time. `individual_results` retains exact case IDs/status/messages;
join its `suite` to the suite row to locate that owner's log. The effective
configuration/context records actual tools, paths, build jobs and selection.
A final terminal summary lists every case and all nonpasses; the transcript
supplies local numeric diagnostics and source-specific failure details.

The complete current receipt is repository-local:
`test_output/remaining-failures-repair-20261110/core-current/amps-core-tests.json`.
It reports **1067 = 825 PASS / 1 FAIL / 240 SKIP / 1 ERROR**, 50 suites,
1958.293 s. This owner reports 45 PASS/1 ERROR; its sole ERROR is the original
CSWC0628 refinement ceiling. All 792 old PASS identities and all 32 new global
cases remain PASS; CMBGU05 changes FAIL to PASS. Strict complete coverage is
the separate T12 FAIL with 26 incomplete components. Full-CME qualification,
actual lateral-map/thermal/work closure, common moving geometry and later
native one-/four-rank propagation remain open.
