# ARCHSCCM01

compiles the shared library against only C++17 and approved `sep_common` headers, scans its dependency graph for forbidden AMPS/PIC, MPI, Tecplot, application, or `swcme` includes, and verifies that both applications link through adapters instead of compiling private copies of the shared physics sources.

Run from any directory with `python3 test.py`. The launcher delegates to the global registry, so individual and cumulative gates execute the identical implementation.
