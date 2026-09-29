# SHK3D07

This is a canonical Stage 7 release test defined in Section 17.2 of `model.md`. Its physics, numerical contract, rejection cases, and acceptance criteria are documented in `docs/STAGE7_TRANSACTIONAL_SHOCK_PROVIDER.md`. The executable assertion is registered in `test/tests_stage7.cpp`.

Run from this directory with `python3 test.py`. The launcher delegates to the global registry, so the individual test and cumulative stage gate execute exactly the same implementation.
