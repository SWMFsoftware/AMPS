# Native initialization readiness fix — 2026-09-30

This package fixes `SCCM3D01` rejecting a prepared, physically inactive shock
before the configured `event.valid_from` time. The supplied example starts the
source at 3600 s while initialization runs at 0 s. The canonical SWCME provider
deliberately returns `active=false`, `generation=0` before that boundary.

## Changes

- `srcSEP3D/main_lib.cpp`: accept a prepared inactive shock during native state
  capture. Require valid status, time coverage, provider/configuration identity
  and finite state fields. A physically active shock still needs a positive
  generation and radius. The activation time and transport/source code are
  unchanged.
- `src/models/sep_coronal_cme/src/runtime_integration.cpp`: retain the strict
  ten-stage ledger check while reporting each missing stage, the completed and
  expected masks, and any unexpected bits.
- `src/models/sep_coronal_cme/test/tests_stage8.cpp`: extend `INIT3D01` with
  independent missing-stage, multiple-stage and unknown-bit negative controls.

The existing shock/background provenance test is outside this narrow fix. In
particular, this change does not install the new coronal model as the production
background provider or fix the separate first-step SIGSEGV.

## Install

Extract this archive from the AMPS root, where `srcSEP3D/` and `src/` reside:

```bash
tar -xzf '/path/to/srcSEP3D_coronal_cme_native_tests_20260930.tar(1).gz'
```

This is the complete original uploaded package with the three source/test
files above updated. It includes the original partial `sep_common` additions;
it is intended to overlay an existing AMPS checkout, not replace the complete
AMPS source tree or its shared models. Existing files not in the archive are
retained. It contains no generated build directories or executables.

Alternatively, apply only `srcSEP3D/initialization_readiness_fix.patch` to the
matching original source tree from its AMPS root. Do not apply the patch after
installing the already-updated sources.

Rebuild the shared archive and application so neither uses the old diagnostic
or readiness predicate:

```bash
make -C src/models/sep_coronal_cme -j8 lib
make clean
make -j8
```

Use the same configured compiler/MPI environment as your existing AMPS build.

## Verify initialization without stepping

Run from the AMPS root:

```bash
mpiexec -n 4 ./amps --test SCCM3D01 --test-input srcSEP3D/examples/sep3d_analytic_parker.in --test-steps 0 --expect-mpi-ranks 4 --test-json test_output/sep3d-init/SCCM3D01.json --artifact-directory test_output/sep3d-init/artifacts
```

If the inactive-shock readiness condition was the only missing stage, the
captured mask becomes 1023 (`0x3ff`) and SCCM3D01 passes. Any genuinely missing
stage continues to fail with a message such as:

```text
initialization output requested before all ten stages: missing=Shock; completed_mask=0x3ef; expected_mask=0x3ff
```

Inspect the captured state with:

```bash
cat test_output/sep3d-init/artifacts/coronal-cme-application-state.txt
```

## Verification performed when preparing this package

- `make -j8 test` in `src/models/sep_coronal_cme`: 196/196 tests passed, including
  the extended INIT3D01 negative controls; adapter compilation also completed.
- An isolated C++ regression harness compiled the actual updated readiness
  expression, the actual ShockState record and its Covers implementation:
  16/16 positive/negative checks passed. The harness supplied only an empty
  patch-record type to avoid unrelated, missing SWCME build dependencies; it
  did not simulate a native AMPS execution.
- The complete linked AMPS executable could not be compiled or run in this
  environment because the upload omits the AMPS core, complete SWCME tree and
  existing shared SEP kernel build files. The zero-step native command above
  must therefore be verified on the user's configured AMPS installation.

This note describes the first initialization-readiness fix. The cumulative
package also includes `ACTIVE_REGION_READINESS_FIX.md`, which documents the
subsequent full-domain/corridor bookkeeping and allocation-evidence correction.
Other findings in the earlier application review are not implemented here.
