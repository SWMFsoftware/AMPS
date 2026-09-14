# WP11--WP20 validation report

Date: 2026-09-14 (UTC)

## Scope

This report records verification performed on the source-only delivery after
implementing WP11 through WP20. It distinguishes dependency-light numerical
evidence from integration and external-data evidence that cannot be produced by
an unattached srcSEP source tree.

## Focused work-package gate

Command:

```sh
make test-wp11-wp20-unit
```

Result: **PASS**.

| Case | Controlled assertion | Result |
|---|---|---|
| WP11 | zero-flux bounded Milstein diffusion preserves domain, symmetry, isotropy, and endpoint behavior | PASS |
| WP12 | named tolerances, step-doubling estimator, limiter histogram, and accepted/rejected accounting | PASS |
| WP13 | coefficient responds only to the selected immutable source view | PASS |
| WP14 | configured constant Dmumu and zero derivative, including invalid inputs | PASS |
| WP15 | Jokipii analytic derivative matches finite-difference oracle with finite endpoint limits | PASS |
| WP16 | Florinskiy branch mirror symmetry, output assignment, and bounded derivative | PASS |
| WP17 | adaptive analytical kappa integral and typed resonance-gap behavior | PASS |
| WP18 | proton/alpha/electron scaling, source normalization, and energy convention | PASS |
| WP19 | named-scale validation, reachability, fingerprinting, and no hard-coded sampling field | PASS |
| WP20 | exhaustive provider/mover policy matrix and finite-to-ballistic type boundary | PASS |

The runner compiled with C++11, `-Wall -Wextra -Wpedantic -Werror`,
AddressSanitizer, and UndefinedBehaviorSanitizer. Leak detection was disabled in
the same manner as the existing source-only runners because the managed test
environment does not expose the process interfaces LeakSanitizer requires.

## Regression gates

The following existing suites were rerun after the implementation:

- Step 1 CLI/registry/reporting: PASS;
- Step 2 background snapshot and single-clock contract: PASS;
- Step 4 exact three-mover API: PASS;
- Step 5 field-line-only scope: PASS;
- WP01--WP10 production contracts: PASS;
- `make test-sanitizer` (Steps 3 and 6--13 plus WP11--WP20): PASS;
- Step 14 documentation, source hygiene, strict warnings, and compiler static
  analysis: PASS;
- Step 15 numerical, cross-mover, cross-model, and SWCME source-only campaign:
  PASS for included cases.

The Step 15 release campaign correctly reported **INCOMPLETE** because this
delivery did not include a checksum-verified real SWMF replay or held-out
spacecraft observations. This is an expected fail-closed evidence status, not a
software-test failure.

## Native validation boundary

The source archive does not contain the enclosing AMPS-generated headers,
configured `Makefile.conf`, MPI runtime, or linked application. Consequently no
native PIC-adapter build or native mover smoke result is claimed here. Complete
that gate in the configured host application with:

```sh
make WARNINGS='-Wall -Wextra -Wpedantic -Werror' lib amps
make test-native-amps-validation SEP_EXECUTABLE=/path/to/amps
```

Real SWMF and observational evidence remain separate commands documented in
`validation/README.md`; neither inherits PASS from the source-only suite.
