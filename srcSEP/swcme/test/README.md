# SWCME validation

This directory contains the standalone SWCME validation executable. It uses
the repository's existing Make-based build approach and calls production
SWCME interfaces directly.

Run the validation suite from this directory with:

```sh
make test
```

Run one registered validation by ID with:

```sh
./output/test_swcme --test CFG01
```

Directory roles:

- `core/`: the test runner and shared validation helpers.
- `1d/`: validations of the production 1-D SWCME model.
- `3d/`: validations of the production 3-D SWCME model.
- `reference/`: reviewed reference data used by comparison tests.
- `profiles/`: inputs describing validation sampling profiles.
- `output/`: generated executables and test results; ignored by Git.

All tests are linked into the single `output/test_swcme` executable.
