# Shared global test-runner support

Updated 2026-11-10. These modules provide configuration, exact registration,
execution/result protocol and bounded native-configuration verification. They
are infrastructure helpers, not physics oracles.

| Module | Responsibility |
| --- | --- |
| `configuration.py` | Layer default/site/environment/CLI settings, validate tools/capabilities/paths and resolve portable argv. |
| `manifest.py` | Expand exact case registries, validate ownership/inventory and retain incomplete components. |
| `runner_support.py` | Load suite contexts, resolve commands, run/terminate process groups, validate typed case/result records and preserve selected IDs on infrastructure ERROR. |
| `native_configuration.py` | Explicitly select and verify no-SPICE mode for fixtures with no target/frame/kernel contract, retaining configuration and real-build receipts before MPI. |

[Global procedure](../README.md) · [Harness controls](../tests/harness/self/README.md) ·
[Current evidence](../STATUS.md).

## Native no-SPICE configuration repair

A component input such as `SPICE=off` does not remove a globally saved SPICE
installation. `.amps.conf`, `Makefile.local`, expanded make variables and the
actual recursive compiler commands can otherwise disagree. The exosphere,
dust-charging, dust-injection and ECSIM-stencil fixtures do not declare a
SPICE target/frame/kernel contract, so their adapters explicitly select a
no-SPICE build and verify it before launching the physical test.

`configure_no_spice(context, application, input_deck, label)` must be called
inside the adapter's existing exclusive checkout lock, after the required
root/idle/generated-build-removal preflight. It does not acquire that lock or
remove `build` itself. Its sequence is:

1. Record the preselection and persist a STARTED JSON receipt.
2. Invoke the established application/configuration workflow with explicit
   `-spice-path=nospice`; never edit `.amps.conf` directly.
3. Read saved application/input/SPICE choices with LAST-assignment semantics.
   Verify the staged input bytes equal the requested real input deck.
4. Record SHA-256 identities of available configuration/make files.
5. Add a private dependency-free make target in the suite work directory and
   inspect actual expanded SPICE, SEARCH_C and FLAGCC from the real Makefile.
   `make -n amps` is not used: recursive make recipes can still execute in a
   nominal dry run.
6. Require saved/effective no-SPICE mode and the correct compiler/search
   definitions; persist PASS or durable ERROR evidence.

The compiler-option parser respects ORDER and both split/combined -D/-U
spellings. These headers use `#ifdef`: `-D_NO_SPICE_CALLS_=0` still defines the
macro, and a later `-U_NO_SPICE_CALLS_` removes it. A substring search would
incorrectly accept some contradictory flag sequences.

After the actual clean top-level build, call
`verify_no_spice_build(configuration, build_stdout)` before any MPI execution.
It requires nonempty C++ compile-command evidence, an active no-SPICE define
in EVERY such command, and absence of `-lcspice` or `cspice.a` commands. An
empty compiler inventory or a recursive override is missing/contradictory
evidence and rejects; it is not accepted because the intended probe passed.

`NativeConfiguration.save()` atomically replaces the JSON receipt through a
temporary file. `reject()` saves status ERROR and raises
`NativeConfigurationError` carrying evidence paths. Callers retain these
paths in their global result. Configuration failure is distinguished from a
failed physical source/charging/sampling invariant.

## Required native build contract

The helper does not replace the repository's native preflight. Before every
actual native rebuild, the owner confirms the AMPS root and idle state,
removes only the generated root `build`, regenerates the selected application
and production hooks through the established workflow, then uses top-level
`MAKEFLAGS=-j16 make amps`. Preserve site libraries, source/input data and
`Makefile.conf`. Selection changes occur through the configuration program.
No direct `.amps.conf` edit, real-repository Git mutation or source-tree reset
is part of this mechanism.

Machine settings come from the supplied profile and suite context; tools are
argv vectors rather than shell strings. The no-SPICE policy is specific to
fixtures without a declared SPICE contract. It does not disable SPICE for
applications that require real targets, frames and kernels.

## Evidence and independent controls

Each caller's artifact directory contains `LABEL-configuration.json` and
`.txt`: expected/saved choices, actual make expansion, exact command/return,
configuration hashes, build-command count and rejection checks. The suite's
native build/run logs and binary/deck identities remain separate evidence.

`HARNESS-NATIVE-CONFIG-ISOLATION-001` executes the maintained adapters with
private synthetic tools: valid enabled/disabled starts must reach their
intended no-SPICE native path, while ignored selection, stale flags,
contradictory make/compiler values or unexpected links must reject before MPI.
The retained repair campaign also runs actual exosphere/dust fixtures from
both starting modes, reverse order and explicit SPICE reselection. These
controls validate this contract, not every configuration option or application.

The completed core campaign preserves every previously passing case, including
all repaired native owners. The 1035-case repair baseline had 792 PASS;
the enlarged 1067-case result has 825 PASS with only the original refinement
ERROR and strict incomplete-coverage FAIL remaining. See STATUS.md for exact
commands, reports and scope. The earlier seven additions relative to the
1028-case preintegration baseline are the five comet-library cases, generic
INP-05 and this native-configuration case; the subsequent repair adds 32 more.

## Result/coverage rules for contributors and agents

Use the maintained `case`, `write_result` and `error_result` interfaces.
`error_result` preserves every selected identity as ERROR when execution/setup
cannot establish its assertion; it does not assert that those cases ran.
Never reduce expected membership to suppress a nonpass. Missing or contradictory
context, stale IDs and status/exit mismatches remain visible.

Complete coverage additionally requires source/configuration/catalog ownership
through `validate_traceability.audit_repository()` and actual independent
scientific evidence. The strict CLI distinguishes incomplete FAIL=1 from
contradictory/setup ERROR=2. A helper or harness PASS cannot promote partial
production components. Every new case must be registered, documented and
added to the generated TEST_CATALOG.md in the same change.
