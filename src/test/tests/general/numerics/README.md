# General numerical and CRC tests

`CORE-GENERAL-NUMERICS-001` retains the ten existing array/ownership/resize/
managed-stack/specfunc/quadrature identities and adds four independent CRC
contracts. The runner mirrors production general headers and supplies a private
empty generated configuration; it never rewrites the application's config.
Compiler, MPI launcher, rank flag and output/work paths come from `--config`.

| New case | What / how / why |
|---|---|
| `PERF-CRC-STREAM-001` | Compare all-byte, unaligned and long pseudorandom streams with an independent bitwise polynomial recurrence, retaining full unsigned-long upper bits |
| `PERF-CRC-SPLIT-001` | Compare many split boundaries, scalar templates, clear and empty appends with the same byte stream |
| `PERF-CRC-COPY-001` | Copy an initialized CRC object, append independently and verify that neither accumulator aliases the other |
| `PERF-CRC-THREAD-001` | Concurrent host construction/update checks shared immutable-table initialization and per-object accumulator independence |

Implementation/reference: `test_general_numerics.cpp`; production class:
`src/general/specfunc.h`. The reference constructs each byte's polynomial
contribution bit-by-bit, without using production tables or slicing code.
On LP64 the legacy return value is not equivalent to a masked standard CRC-32;
the test explicitly seeds and checks the full accumulator. Device execution
and unavailable alternative architectures remain separate qualification work.

```sh
python3 -B src/test/run_tests.py \
  --config test/qualification-machine/passion.json \
  --match CORE-GENERAL-NUMERICS-001 \
  --output-dir test_output/general-numerics
```

Use a fresh directory and a site overlay for local paths. Each case executes
in an isolated process and needs one exact result plus a successful exit;
the array-bounds identity additionally retains its two-rank MPI gather checks.
Logs are preserved as `general-numerics-output.txt`. All 14 cases passed the
targeted performance-preservation run. The new native
`PIC-PBUFFER-SIGNATURE-001` separately verifies canonical compiled particle
bytes, list inclusion/exclusion and source immutability through the actual PIC
API. New permanent identities also require manifest/catalog/README updates.
Keep their suite/case ownership current in `component_inventory.json` and
`source_traceability.json`; compiled case sources must be listed under the
latter's `case_registry_sources`. Agent contributors must run traceability
and generated-catalog checks after registration changes, preserving these
independent polynomial and numerical oracles.
