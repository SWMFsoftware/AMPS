# POP3D01

This is a canonical Stage 8 release test defined in Section 17.2 of `model.md`. Its physics, numerical contract, rejection cases, and acceptance criteria are documented in `docs/STAGE8_RUNTIME_MESH_AND_OBSERVERS.md`. The executable assertion is registered in `test/tests_stage8.cpp`.

Run from this directory with `python3 test.py`. The launcher delegates to the global registry, so the individual test and cumulative stage gate execute exactly the same implementation.
