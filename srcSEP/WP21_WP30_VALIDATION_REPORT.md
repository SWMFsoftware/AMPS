# WP21–WP30 validation report

## Executed source-only gate

Command: `./test/run_wp21_wp30_tests.sh`

Result: **PASS**. The strict C++11 `-Wall -Wextra -Wpedantic -Werror` build with
AddressSanitizer and UndefinedBehaviorSanitizer completed, and all ten work
package assertions plus the production-wiring/static-safety gate passed.

| Work package | Controlled evidence | Status |
|---|---|---|
| WP21 | knot values, bitwise timestep partition, serialized restart | PASS |
| WP22 | outside-start double crossing, orientation, outward policy, degenerate rejection | PASS |
| WP23 | normalized power-law quadrature and index-one limit | PASS |
| WP24 | identical semantic key reproducibility and swapped-field separation | PASS |
| WP25 | relative-normal source value, exact branch closure, zero-relative-speed limit | PASS |
| WP26 | `A|B|=Phi`, conservative refinement, versioned roundtrip | PASS |
| WP27 | relativistic local-field Larmor reference and superluminal rejection | PASS |
| WP28 | exact edge/underflow/overflow semantics, effective sample size, density/flux normalization | PASS |
| WP29 | transactional completion, duplicate/path rejection, corruption detection | PASS |
| WP30 | layer precedence, immutable serialization fingerprint, restart mismatch rejection | PASS |

The runner also confirmed production source references to the sphere kernel,
queued turbulence ledger, normalized inverse CDF, field-line flux registry,
strict invalid-particle path, transactional writer, and startup fingerprint.

## Separate gates not executed in this source-only workspace

| Evidence class | Status | Requirement |
|---|---|---|
| Native AMPS/PIC strict-warning build and registered runtime tests | BLOCKED | enclosing AMPS `Makefile.conf`, generated PIC headers, MPI, linked executable |
| MPI/OpenMP decomposition and long-run turbulence/source ledger | BLOCKED | native executable and multi-rank runtime |
| Real SWMF replay | BLOCKED | checksum-authenticated SWMF manifest/data |
| Held-out observational validation | BLOCKED | authenticated event products and instrument forward operators |

No blocked gate is reported as passed by the dependency-free test.

## Regression gates

The following pre-existing source-only gates were re-run after WP21–WP30:

- `make test-sanitizer`: PASS, including Steps 3 and 6–13, WP11–WP20, and the
  new WP21–WP30 suite;
- `./test/run_wp01_wp10_tests.sh`: PASS;
- `./test/run_wp11_wp20_tests.sh`: PASS;
- `./test/run_step1_tests.sh`: PASS (its printed `HIDDEN01 FAIL` is the
  intentional fixture proving hidden callback failures are detected; the
  encompassing `HIDDEN PASS` and process exit status are authoritative); and
- `make test-documentation-unit`: PASS, including docs, public-surface,
  archive-hygiene, strict-warning, static-analyzer, and sanitizer assertions.
