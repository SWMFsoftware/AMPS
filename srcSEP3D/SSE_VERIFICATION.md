# Finite SSE and weak-flank fix verification — 2026-10-02

The cumulative source overlay includes finite SSE support and the exact
180-s weak-shock runtime correction. See [SSE_WEAK_SHOCK_FIX.md](SSE_WEAK_SHOCK_FIX.md)
for the reproduced failure, algorithm changes and target-system rebuild/run.

| Check | Result | Scope |
| --- | --- | --- |
| Portable standalone application suite | 148 passed; 0 failed/errors/skipped | Existing layers plus nine finite-SSE tests |
| Focused final SSE suite | 9 passed | Geometry, fields, particle/source handoff, restart, exact failed coordinates through ten steps, actual weak layer, and genuine unresolved-layer rejection |
| AddressSanitizer + UndefinedBehaviorSanitizer | 24 passed | SSE3D, SHK3D, RST3D, SWBG3D groups; includes both new runtime regressions |
| Canonical shock regressions | 14 test functions; 2,553 assertions passed | SHK01–12, SHK15, SHK16; includes the unchanged independent oblique/weak references and 100,000 deterministic physical stress inputs |
| Weak random population | All 25,000 solved | Excesses 1e-10..1e-6 now resolve conserved evolutionary near-identity jumps |
| Shared canonical 1-D/3-D integration suite | 4 passed | Configuration, one-dimensional, three-dimensional and cross-dimensional integration |
| Independent failed-ray reference | Passed at 80 decimal digits | Direct conserved-flux bisection, independent of the production cubic |
| Production integration audit | 11 passed; 1 skipped | Static integration and compiled header/build probes; configured native build unavailable |
| Receive-mask fixture | All 13 checks passed | Six face masks select current data, six stale receives reject, one unmasked receive passes |
| Runtime DATAFILE ownership fixture | 84 assertions passed | 42 checks each with temporal interpolation off/on |

The optimized portable build uses C++17 and
`-Wall -Wextra -Wpedantic -Werror -O2`. Sanitizer builds add
`-fsanitize=address,undefined -fno-omit-frame-pointer -g`; runs use
`ASAN_OPTIONS=detect_leaks=0:halt_on_error=1` and
`UBSAN_OPTIONS=halt_on_error=1`. LeakSanitizer is not included in the claim.
After the full suite and sanitizer run, only comments, documentation and the
SSE3D08 reference literal changed (one rounding unit, within the unchanged
5e-13 comparison tolerance); the optimized focused suite passed again.

The canonical shock checks invoke the actual production test functions and
shared model archive. The 100,000 stress inputs are fixed-seed physical states,
not native MPI mesh cells. Evidence includes the complete application log,
focused/sanitizer JSON, canonical shock log, exact before/after reproduction,
and independently verified reference output.

## Native limitation

There is no configured NASA AMPS build or HPE MPT environment here. The rebuilt
four-rank ten-step sep-corona test has **not** been executed. The receive-mask
fixture compiles production packing/capture functions but does not simulate
MPI transport. All native MPI transport, decomposition and publication checks
must still be confirmed with the target command in SSE_WEAK_SHOCK_FIX.md.

## Source package

All archive paths are relative to the AMPS root. It contains directly modified
sources, comments, READMEs, the complete input, tests and evidence; no installer,
object file or binary is included. `FINITE_SSE_SOURCE_MANIFEST.sha256` records
all packaged file contents. The cumulative overlay preserves the earlier
polar Rankine-Hugoniot, runtime DATAFILE ownership and receive-mask corrections.
Passing these checks verifies numerical/software contracts; the physical
model remains the prescribed finite front described in SSE.md.
