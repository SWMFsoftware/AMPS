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

When adding a new software test, add its method to `test_protocol.py` or
explicitly register its TestCase in this runner's loaded module. Keep native
adapters thin and update the contract before implementation. Existing Moon
U15/U16/U17 numerical/statistical runners retain their original independent
oracles and acceptance criteria; the shared-service move changes their source
link lists and hashes, not their physical scope. No permanent U/I ID is added.
