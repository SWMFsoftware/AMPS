# Native build support

The native exosphere, dust charging/injection and ECSIM stencil fixtures have
an explicit global no-SPICE contract. Their adapters call Config.pl with
`-spice-path=nospice` in addition to application selection. An exosphere input
block's `SPICE=off` alone does not clear a saved global library installation.
The shared `src/test/lib/native_configuration.py` helper checks the selected
application/input, saved mode, effective make `SPICE`, `SEARCH_C` and `FLAGCC`,
then corroborates the actual clean-build compiler/archive commands. It uses a
private dependency-free inspection target, not `make -n amps`; recursive make
commands can execute even during a make dry run.

Each owner retains a `*-configuration.json` receipt and configure/probe log,
including initial selection, exact command, observed options, check outcomes
and configuration hashes. Build identity points to that receipt. A mismatch
is ERROR before native execution and preserves the requested case IDs. The
helper is called inside each owner's existing exclusive lock/preflight; it
does not change global runner defaults or SPICE-dependent application owners.
`HARNESS-NATIVE-CONFIG-ISOLATION-001` checks valid and injected-invalid control
flows using isolated synthetic tools; it does not qualify the physics cases.
Those cases require real clean native builds and their independent oracles.

The global runner's native owners use their application adapters (for example
`../tests/pic/native_core/run_test.py`) and the exclusive checkout gate.
Legacy `make test_Moon -j` and `make test_fast-wave -j` instead
call `prepare_legacy_native.py` immediately before `Config.pl`. Run these make
targets separately from a global campaign and separately from each other:

```sh
make -j test_Moon
make test_fast-wave -j
```

The legacy preflight requires the repository root as its working directory,
inspects Linux `/proc` for competing builds/runs in that checkout, refuses a
symlink at `build`, and removes exactly the generated `build/` directory.
Ancestor make processes belong to the current invocation and are excluded
from the conflict scan. Configuration then selects the application and
regenerates its production hooks; the AMPS compilation uses `-j16`, even when
the outer make has unlimited parallelism. Configuration, compilation and test
execution must remain serialized against other native owners. This helper
is an idle-state check, not an interprocess locking service for arbitrary
external commands; an operator must not launch a competing build after the
check. It preserves site configuration and libraries and performs no Git
operation.

The Linux process guard is deliberately conservative: it treats a competing
`make` process below the checkout as a possible native build even when the
target is documentation-only. Schedule `make PDF` and other make targets
before or after a native campaign. A refused preflight is an execution ERROR,
not physics evidence; retain its first log and rerun the blocked owner from
an idle root. Do not disable this guard to make a report green.

The recipes are maintained in `MakefileTest/Table`; `Makefile.test` is the
site-generated counterpart. The fast-wave AMPS compilation propagates build
failure instead of continuing with an older executable. Historical reference
comparison recipes can still mask a numerical failure, so integration
validation must inspect every generated `.diff` file and child log in
addition to the make exit status. Application SPICE/input/reference paths
remain site settings; this helper does not replace them.

For fast-wave only, the native compile recipe forwards `FAST_WAVE_CXXFLAGS`
through `FLAGC_EXTRA`. The root Makefile supplies `-O1` for GNU host C/C++
when no explicit optimization level is already configured; existing site
levels, NVCC and other compiler families retain their policy. Disable this
fallback with `make test_fast-wave -j FAST_WAVE_CXXFLAGS=`. Site files and
Fortran flags are preserved. The measured complete command improves from
119.04 s to 87.08 s on passion, with all six sampled files byte-identical and
the original comparison bounds unchanged. See
[`test/srcFastWave/README.md`](../../../test/srcFastWave/README.md) for exact
commands, timings, validation scope and evidence locations.

Optional individual-function CPU sampling and MPI call-site measurement for
this fixed fixture are documented in
[`test/srcFastWave/profiling/README.md`](../../../test/srcFastWave/profiling/README.md).
They reuse a frozen validated binary and repeat its unchanged checks in
private directories; compiling the diagnostic interposer is not an AMPS
rebuild. The resulting optimization proposal retains all native/core gates.

The subsequent algorithm implementation is documented in
[`test/srcFastWave/STATUS.md`](../../../test/srcFastWave/STATUS.md).
Its final native regression takes 79.387 s total / 26.326 s simulation on
passion and preserves six exact full-file hashes. The native particle owner
adds `PIC-PBUFFER-SIGNATURE-001`: its compiled particle-record fixture checks
canonical CloneParticle-equivalent bytes, an independent bitwise full-width
CRC, list-field inclusion/exclusion, incremental state and source immutability.
It runs as part of `CORE-PIC-NATIVE-001`; no separate configuration or
production signature mode is introduced. MPI transport, stencil and sparse
unit fixtures are separate nonexclusive owners under
[`tests/pic/performance`](../tests/pic/performance/README.md); the core registry,
component inventory, source traceability and generated catalog all register
the additions. Native builds still serialize globally against other owners.
