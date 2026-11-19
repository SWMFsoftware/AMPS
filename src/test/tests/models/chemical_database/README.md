# Shared chemical database CORE adapter

`CORE-CHEMICAL-DATABASE-001` delegates its
`CHEMICAL-DATABASE-PROTOCOL-001` case to the maintained
[`test/chemical_database` suite](../../../../../test/chemical_database/README.md)
with explicit `--protocol-only`. The contract was recorded before source/test
changes in [CONTRACT.md](../../../../../test/chemical_database/CONTRACT.md).

The case compiles actual shared Runtime/DataIO and checks namespace aliases,
known SHA vectors, independent hashlib, exact quoted/padded CSV labels/rows,
hash/schema/path/capability/rights refusals, missingness, resource limits,
immutable lifetime/concurrent readers, transactional refusal, relocation,
Python API/CLI and runner-result classification. Administrative fixtures are
labelled parser/invariant strings, never observations. No site database,
network, PIC, MPI, SPICE or application configuration is required.

Run from the repository root with the machine profile used by the full CORE
campaign and a new evidence directory:

```sh
python3 -B src/test/run_tests.py --config test/qualification-machine/passion.json \
  --match CORE-CHEMICAL-DATABASE-001 --output-dir test_output/core-chemical-db
```

The adapter preserves the configured compiler argv (including wrappers/flags)
and the standard global timeout/process-group ownership. Child `result.json`,
`unittest.log`, build/executable/source hashes and per-query logs are retained
under the case artifact directory's `shared-protocol`. Missing dependencies
remain SKIP, evaluated disagreement FAIL, corrupt/missing/contradictory results
or build/runtime failures ERROR. The full exact child registry/counts must
match its declared status/exit; an available PASS cannot hide a required SKIP.

Real-release row/source agreement and installed SDK relocation remain separate
explicit `--database-lock` work; physical model validity, CUDA, source-rights
holds, kinetic completeness and formal application gates are not established
by this case. Component coverage therefore remains `partial`; the existing
strict incomplete-coverage failure is preserved. No Moon U/I ID is added.
