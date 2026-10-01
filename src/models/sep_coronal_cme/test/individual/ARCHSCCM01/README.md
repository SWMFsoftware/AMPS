# ARCHSCCM01

Checks the shared C++ source/header boundary for forbidden AMPS/PIC, MPI,
Tecplot, application and `swcme` includes, rejects upward dependencies from
neutral `sep_common` sources, scans the built archive's undefined symbols, and
compiles an external C++17 consumer of the public headers. Adapter compilation
is a separate `make check-adapters` check included by the aggregate driver.

The text scan selects regular C/C++ source/header files by suffix. Compiled
`sep_*.o`, archives, shared libraries, dependency files and source-named
directories are excluded from text decoding. This fixes the UTF-8 traceback
seen when AMPS leaves compiled objects beside `sep_common` sources. The archive
symbol audit remains active; excluding object bytes from text decoding does
not relax the model's dependency boundary. An invalid UTF-8 source file still
fails explicitly with its path.

Five isolated regressions run inside this same canonical gate: adjacent binary
artifacts/directories, nested/alternate-suffix forbidden includes, upward
neutral dependencies, invalid source encoding, and valid Unicode comments.
They are followed by the real archive/public-header checks. They do not add
canonical IDs or change the release suite's test count.

Run from any directory with `python3 test.py`. The launcher delegates to the global registry, so individual and cumulative gates execute the identical implementation.

From the AMPS root, rerun only this gate with:

```sh
python3 src/models/sep_coronal_cme/test/run_tests.py --test ARCHSCCM01
```

The Python audit update needs no AMPS executable rebuild. Its existing shared
library archive must be present; `make -C src/models/sep_coronal_cme lib` builds
it if necessary.
