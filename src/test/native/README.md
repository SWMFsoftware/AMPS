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

The recipes are maintained in `MakefileTest/Table`; `Makefile.test` is the
site-generated counterpart. The fast-wave AMPS compilation propagates build
failure instead of continuing with an older executable. Historical reference
comparison recipes can still mask a numerical failure, so integration
validation must inspect every generated `.diff` file and child log in
addition to the make exit status. Application SPICE/input/reference paths
remain site settings; this helper does not replace them.
