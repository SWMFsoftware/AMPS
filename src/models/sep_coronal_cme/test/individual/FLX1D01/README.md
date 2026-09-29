# FLX1D01

This is a canonical Stage 10 release test defined in Section 17.2 of `model.md`. Its physics, numerical contract, rejection cases, and acceptance criteria are documented in `docs/STAGE10_FIELD_LINE_EXCHANGE.md`. The executable assertion is registered in `test/tests_stage10.cpp`.

Run from this directory with `python3 test.py`. The launcher delegates to the global registry, so the individual test and cumulative stage gate execute exactly the same implementation.
