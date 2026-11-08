# T12 portability, sanitizer, backend, and guardrail matrix

This directory implements a bounded part of Stage T12.  It is deliberately
not a claim that T12 or complete AMPS coverage is closed.  The global adapter
compiles a dependency-light executable against the production
`src/general/array_1d.h` ownership/indexing implementation and runs it through
compiler, sanitizer, OpenMP, leak, and performance profiles.  It also reuses
the maintained C++/Fortran ABI fixtures and the two-rank `array_3d` gather
fixture.  Expected values are independent closed-form checksums or explicit
ABI constants; no production result is copied into an oracle.

## Exact cases

| Case | Contract and result policy |
|---|---|
| `T12-CXX-GNU-DEBUG-001` | Configured MPI C++ wrapper, `-O0 -g`, production-array ownership and checksum. A compile-time guard requires `__GNUC__` and rejects Clang's compatibility definition. A present compiler failure is `FAIL`. |
| `T12-CXX-GNU-OPTIMIZED-001` | Same guarded GNU contract under `-O3 -DNDEBUG`; its executable is reused by bounded OpenMP/performance/leak cases. |
| `T12-CXX-CLANG-001` | Clang through the OpenMPI wrapper. Missing Clang/wrapper is a typed `SKIP`; a present compiler failure is `FAIL`. |
| `T12-FORTRAN-ABI-001` | Both maintained legacy-underscore and ISO-C-binding C++/Fortran ABI controls. C++ translation units use the configured C++ driver; the selected Fortran driver owns the final link and therefore supplies its own GNU, Flang, or Intel runtime rather than a hard-coded `libgfortran`. Missing Fortran is `SKIP`. |
| `T12-SANITIZER-ASAN-001` | Production-array ownership under AddressSanitizer. LeakSanitizer is disabled because ptrace-supervised runners cannot execute it reliably; the independent Valgrind case owns leak detection. |
| `T12-SANITIZER-UBSAN-001` | Same fixture under UndefinedBehaviorSanitizer with recovery disabled. |
| `T12-SANITIZER-TSAN-001` | A real four-thread OpenMP fill/reduction of the production array under ThreadSanitizer, requiring its exact four-thread PASS marker. Missing TSan support files or a known host address-map incompatibility is `SKIP`; a data-race report is `FAIL`. This is not an MPI+TSan qualification. |
| `T12-VALGRIND-LEAK-001` | Full leak inventory of the optimized ownership fixture; definite, indirect, or possible loss fails, while process-lifetime runtime-library reachable bytes remain recorded but are not treated as a production-array leak. Missing Valgrind is `SKIP`. |
| `T12-OPENMP-INVARIANCE-001` | Static 1/2/4-thread fills and reductions have the same independently known `N(N-1)/2` checksum. |
| `T12-MPI-OPENMPI-001` | The maintained two-rank production `array_3d` gather identity/rejection contract with the explicit OpenMPI wrapper/launcher pair. |
| `T12-MPI-MPICH-001` | The same source and oracle with the explicit MPICH pair. Missing implementations are separate typed `SKIP`s. |
| `T12-CUDA-COMPILE-001` | NVCC compiles a kernel including production `constants.h`; this is header/device-compiler evidence only. |
| `T12-CUDA-RUNTIME-001` | The constants kernel runs on a device. Compiler present but no usable driver/device is `SKIP`, never a CPU PASS. |
| `T12-PERFORMANCE-ARRAY-001` | One-million-element production-array fill/checksum, eight repetitions, with a broad 30 s fixture guardrail. |
| `T12-PERFORMANCE-IO-001` | Two-megabyte production-array binary round trip with checksum and a broad 15 s fixture guardrail. It is not the AMPS restart codec. |
| `T12-COVERAGE-STRICT-001` | Runs `--coverage --require-complete-coverage` and reports its real exit status. The current incomplete inventory is therefore `FAIL`, keeping T12 open. |

The performance cases separate the correctness invariant from elapsed time.
Their generous ceilings detect catastrophic local regressions or unusable
filesystems, not small benchmark variation.  They do **not** satisfy the T12
requirements for configured mesh construction, particle push, collision,
halo exchange, and restart-I/O performance; those remain open until each has a
production-linked, machine-calibrated budget.

## Portability and machine configuration

Required baseline commands (`python`, `cxx`, `mpi_cxx`, `make`, `mpiexec`) are
resolved by the global test configuration. Optional tools first honor complete
argv vectors in the machine profile, then conventional executable fallbacks:
`clang_cxx` (`clang++`), `valgrind` (`valgrind`), `nvcc` (`nvcc`),
`openmpi_cxx`/`openmpi_exec` (`mpicxx.openmpi`/`mpiexec.openmpi`), and
`mpich_cxx`/`mpich_exec` (`mpicxx.mpich`/`mpiexec.mpich`). The existing
`fortran` vector similarly precedes common Fortran driver names. No vector is
passed through a shell, so wrappers and their arguments retain exact argv
boundaries. Optional absence is a stable per-case `SKIP`; hard-coded absolute
paths are forbidden.

The runner rejects a context whose owning `expected_cases` differs in content
or order from its complete 16-case registry before applying an optional
`selected_cases` subset. This makes an omitted manifest identity a protocol
ERROR instead of a smaller green run. A Valgrind timeout is localized to an
`ERROR` for its own case, allowing the remaining identities to execute and
appear in the receipt.

All outputs and compilation products live below the global suite work and
artifact directories.  A private mirror supplies the otherwise generated
`.general.conf`, so the suite neither reads a stale application selection nor
changes repository configuration.  It never invokes `Config.pl`, removes
root `build`, or performs a native AMPS rebuild.

After registration in the global manifest, run:

```sh
python3 src/test/run_tests.py --test CORE-T12-PORTABILITY-001 \
  --output-dir /tmp/amps-t12-portability --jobs 1
```

A nonzero final exit is expected while `T12-COVERAGE-STRICT-001` remains red.
Inspect the final individual table and
`artifacts/CORE-T12-PORTABILITY-001/t12-portability-matrix.txt`; do not relabel
that expected open gate as a suite PASS.

## Known qualification limits

- No native AMPS application is rebuilt by this suite.
- Clang, TSan, CUDA runtime, Valgrind, Fortran, and the second MPI
  implementation may be unavailable and remain visible SKIPs.
- Sanitizer coverage is the production general-array boundary, not all PIC,
  mesh, model, MPI, or application code.
- The OpenMP case checks deterministic disjoint writes/reductions; it does not
  qualify concurrent particle lists or hybrid MPI+OpenMP transport.
- CUDA compilation/runtime does not qualify a CUDA particle mover or backend.
- The strict coverage gate remains authoritative; the bounded green cases
  cannot upgrade partial/planned inventory entries.
