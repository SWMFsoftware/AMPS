# SRC3D06

This is a canonical Stage 9 release test defined in Section 17.2 of `model.md`. Its physics, numerical contract, rejection cases, and acceptance criteria are documented in `docs/STAGE9_MOVING_SHOCK_SOURCE.md`. The executable assertion is registered in `test/tests_stage9.cpp`.

Run from this directory with `python3 test.py`. The launcher delegates to the global registry, so the individual test and cumulative stage gate execute exactly the same implementation.
