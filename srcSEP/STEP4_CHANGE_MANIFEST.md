# Step 4 change manifest

> Step 14 supersession: the documented transition-alias interval has ended.
> Only `parker`, `fte-dmumu`, and `fte-mfp` now parse; see
> `MIGRATION_MANIFEST.md` for replacements.

Step 4 introduces the three-mover public production API on top of the completed
Step 3 SI geometry/source normalization.

## New files

- `PRODUCTION_MOVER_API.md`
- `production_mover_runtime.cpp`
- `util/sep_production_mover.h`
- `util/sep_production_mover.cpp`
- `test/run_step4_tests.sh`
- `test/step4/test_production_mover_cli.cpp`
- `STEP4_CHANGE_MANIFEST.md`

## Principal modified files

- `sep.h` routes PIC callbacks through the validating adapter;
- `util/sep_cli.h` and `util/sep_cli.cpp` expose only canonical production
  choices; Step 14 removed the temporary transition aliases;
- `main.cpp` implements pre-initialization `--list-movers`, prints complete
  mover/coefficient metadata, and uses capabilities for turbulence dispatch;
- `makefile`, `test/run_step1_tests.sh`, `README.md`, `test/README.md`, and
  `CURRENT_CLI_CONTRACT.md` include the new source and tests.

## Quick verification

```sh
make test-mover-api-unit
make test-cli-unit
```

Generated objects, archives, sanitizer files, coverage files, and test output
must not be committed.
