# Steps 10–12 validation report

Validation date: 2026-09-13 (UTC)

## Completed source-only gates

- Step 10: `COEF01`–`COEF05` and `COEF-SOURCE` passed.
- Step 11: `TURB01`–`TURB20`, `TURB-CLI`, and `TURB-SOURCE` passed.
- Step 12: `PAR01`–`PAR05` and `PAR-SOURCE` passed.
- Regression: every focused script from `run_step1_tests.sh` through
  `run_step12_tests.sh` passed sequentially.
- Compiler policy: the Step 10–12 dependency-light binaries compile as C++11
  with `-Wall -Wextra -Werror -pedantic`, AddressSanitizer, and
  UndefinedBehaviorSanitizer.
- Source audit: no object, archive, shared-library, Python-cache, or test-binary
  build products are included in the prepared source tree.
- Patch audit: `git diff --no-index --check` against the supplied baseline
  reports no whitespace errors after the final correction.

## Reproducibility evidence policy

PAR02 and PAR03 feed the same 200 physical contributions through different
worker and synthetic MPI-rank layouts. The canonical reduced arrays produce an
identical 64-bit evidence hash. PAR04 proves that drawing from a diagnostic
purpose stream does not change the mover-purpose stream. This is bitwise
evidence for the dependency-light reduction core; it is not a claim about an
uncontrolled third-party MPI collective tree.

## Native-build boundary

`make -j2 lib` cannot run from this attachment because the enclosing AMPS file
`../../Makefile.conf` is absent. Consequently, this report does not claim a
linked PIC/AMPS executable, real OpenMP/MPI process-count comparison, live SWMF
handoff/restart, or long-duration scientific/observational validation. Those
remain explicit native integration gates on a complete AMPS checkout.
