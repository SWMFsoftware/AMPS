# Shared chemical database software checks

This is the single global runner for the shared C++/Python data-access service.
It performs no AMPS configuration, particle simulation or observational
validation. Read [CONTRACT.md](CONTRACT.md) for the preimplementation inventory,
scientific boundaries and exact software acceptance rules.

From the repository root, choose a new output directory:

```bash
/usr/bin/python3 -B test/chemical_database/run_all_tests.py \
  --output-dir test_output/chemical-database/my-new-run \
  --database-lock /data/vtenishe/CHEMICAL_DATABASE/database-lock.json
```

The optional `--database-lock` selects a real prepared release for additional
Python/native exact-source agreement and actual installed CLI/C++ SDK checks.
The suite copies that release (never hardlinks), relocates the copy, runs its
own `bin/chemdb` and compiles all six shared CPPs from its own `sdk/` tree
without AMPS checkout includes/libraries. `installed-deployment/result.json`
records those deployed source/executable/lock identities. The original release
and canonical sources remain read-only. Omit the lock to run the native and
Python administrative metadata/integrity controls; dependent real-release
checks report SKIPPED. Use the actual
installation lock filename written by the builder. No latest-release search,
implicit download or alternative source is used. Python3 and a C++17 compiler
(`g++`, or the compiler argv named by `CXX`) are required. No MPI/PIC/SPICE is used.
The CORE adapter explicitly passes `--protocol-only` for the site-independent
software selection. It preserves its configured compiler argv through
`--cxx-json`; there is no site-data prerequisite hidden in that CORE case.

The native controls compile the actual maintained `Runtime.cpp` and `DataIO.cpp`
and call them through a thin read-only adapter. Administrative fixtures contain
only labelled parser strings, never observations or numerical chemistry data.
They verify exact SHA-256 binding, published empty/`abc` SHA vectors,
schema/codec ownership, ordered identifiers, quoted/newline CSV strings and
exact padded column labels, capability and permission refusal, unavailable
metadata, resource bounds, path and symlink confinement, immutable shared
lifetime and concurrent readers, null/missing-output guards, failed-reload
rollback, row counts and relocation. Compile-time type identities prove that
the historical Moon namespaces alias the shared services. Real installed and
canonical source bytes are read-only. Altered fixtures remain inside this run.

`result.json` contains every evaluated test or fixture outcome and the executed,
PASS, FAIL, SKIPPED and ERROR totals. Every failed/error entry names the test,
runner and exact log path. `unittest.log`, `build.json`, compiler stdout/stderr
and `queries/query-*/{stdout.log,stderr.log,result.json}` preserve the actual
commands and responses. Prior outputs are never overwritten. Exit codes are
0=PASS,1=FAIL,2=ERROR,77=SKIPPED. A prerequisite fixture outcome is counted as
such; it does not pretend its dependent methods were executed.

When adding a new software test, add its method to `test_protocol.py` or a
maintained `test_*.py` module in this directory; the runner discovers its
`unittest.TestCase` methods automatically and binds their source hashes. Keep native
adapters thin and update the contract before implementation. Existing Moon
U15/U16/U17 numerical/statistical runners retain their original independent
oracles and acceptance criteria; the shared-service move changes their source
link lists and hashes, not their physical scope. No permanent U/I ID is added.

The November 19 source-acquisition checks cover strict packet admission,
NASA CEA fixed-width records and domains, IAEA native printed fit coefficients,
ExoMol source projections and TU Wien native CSV members. Administrative parser
fixtures are labelled as such. Real-source tests read the explicit selectors
below; absent inputs are SKIPPED and never substituted with generated science:

```bash
AMPS_CEA_SOURCE=/path/to/acquisition/raw/nasa-cea-data-thermo.inp \
AMPS_EXOMOL_PACKET=/path/to/exomol-native \
AMPS_SURFACE_SOURCE_ZIP=/path/to/As_published_in_Comm_Earth_Env.zip \
/usr/bin/python3 -B test/chemical_database/run_all_tests.py \
  --database-lock /data/vtenishe/CHEMICAL_DATABASE/snapshots/amps-chemistry-20261119-r2/database-lock.json \
  --output-dir test_output/chemical_database/NEW_REAL_RELEASE_RUN
```

[`verify_acquired_release.py`](ACQUIRED_RELEASE_QA.md) separately checks the actually published new tables
against their acquired originals, preserving a result/manifest and per-check
logs. Invoke its `--help` for the required packet and lock arguments. This is
source reproduction and deployed API verification, not linked I21–I24 or
observational validation. Both old and new locks can be tested independently;
no test changes a selected source, threshold or default installation.

The r3 follow-up adds `test_native_chianti.py`, `test_native_apollo_tpd.py`,
`test_native_workbooks.py`, `test_model_definitions.py` and
`test_reduced_models.py` to the same automatic discovery. The first three keep
their actual-source classes marked `requires_authoritative_source=True`:
missing originals are SKIPPED in the full suite and explicitly outside the
CORE protocol selection. Generic parser controls are labelled administrative,
and analytic/conservation tests exercise the actual shared C++ helpers. Model
metadata controls prove old-lock compatibility and that definitions grant no
physical capability. No lunar U/I ID or acceptance threshold is changed.

Select CHIANTI originals with `AMPS_CHIANTI_LEVELS=/path/to/o_1.elvlc` and
`AMPS_CHIANTI_RADIATIVE=/path/to/o_1.wgfa`, and Apollo originals with
`AMPS_APOLLO_TPD_SOURCE_DIR=/path/to/jones2020-native`; the existing CEA/ExoMol
selectors above remain applicable. Workbook real-source controls require the
exact archived ZIPs documented in `NATIVE_WORKBOOK_SOURCES.md`. The global
runner records all actual selections and distinguishes source access from
physical agreement. Choose r3's own lock to exercise its contents.

[`verify_physical_model_release.py`](PHYSICAL_MODEL_RELEASE_QA.md) checks an
explicit r3-style release and selected packet family through Python and C++,
including exact rows, recipe members, held-source refusal and model-registry
agreement. Use `--help`; supply the compiled shared `catalogue-probe` with
`--native-probe`, the actual lock and repeated `--source-packet` arguments.
Its output directory must be new. This is shared source/API verification, not
formal linked I21–I24 qualification or observational validation.

[REDUCED_MODEL_NUMERICS.md](REDUCED_MODEL_NUMERICS.md) records the independent
Decimal/Fraction oracles, retained development failures and the proof behind
the binary floating-point zero-COM certificate. Mathematical roundoff criteria
are distinct from physical approximation or observational agreement limits.
