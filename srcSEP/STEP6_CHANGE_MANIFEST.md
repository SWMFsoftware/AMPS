# Step 6 change manifest: common field-line transport core

## Implemented

- Added `util/sep_transport_common.*` with explicit status values, shared
  particle-state validation, strict relativistic SI conversion, the common
  focusing convention, exact plasma-frame cooling, boundary policies, named
  composable timestep limits and diagnostics, and keyed random streams.
- Added `transport_common.*` as the sole PIC-facing load, physical field-line
  advance, immutable local-background read, commit, and segment-attachment
  boundary.
- Refactored `parker_mover.cpp`, `focused_transport_dmumu.cpp`, and the interim
  mean-free-path shell to enter through the common boundary. Step 9 subsequently
  replaced that shell with canonical `focused_transport_mfp.cpp`.
- Replaced the legacy forward-Euler cooling/floor behavior in `cooling.*` with a
  status-returning adapter over the exact common kernel. The old scalar function
  remains only as a compatibility wrapper and reports invalid input as NaN.
- Added `CORE01`–`CORE07` plus the `CORE-SOURCE` production-shell check and the
  `make test-transport-common-unit` target.

## Removed duplication and hidden policy

The canonical shells no longer own separate particle validation, raw coordinate
arithmetic, segment insertion, relativistic conversion, cooling formulas,
pitch-angle reflection, or global random draws. A substep below the declared
minimum is an error rather than an implicit clamp.

## Remaining native gate

This archive does not contain the enclosing `Makefile.conf`, PIC headers, MPI
configuration, or linked executable. All three PIC-facing mover shells must
therefore still be compiled together and exercised in a complete field-line
AMPS checkout before Step 6 native integration is signed off.
