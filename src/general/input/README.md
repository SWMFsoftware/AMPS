# Generic AMPS runtime-input reader

`amps_input.h` is the header-only reader of the generic AMPS runtime input
format (`amps.in`, `-input FILE`).  The normative specification — grammar,
block kinds and owners, the srcSEP/srcSEP3D layouts, the model files and the
migration from the retired formats — is [`../../INPUT_FILE.md`](../../INPUT_FILE.md).

## What the reader does and does not do

* It parses the grammar only: `#begin KIND[: NAME]` / `#end KIND` blocks,
  `KEY = VALUE` assignments, `!` comments, `\` continuation and `#include`
  (relative to the including file, cycle-checked, at most 64 deep).
* It enforces the structural rules (assignments inside blocks, unique keys
  per block, unique sibling blocks, block-balanced files, no unknown
  directives) and reports every violation as `file:line: reason` plus the
  offending line.
* It never interprets a key.  Each block is validated by its owner (the AMPS
  core, an application, or a model library under `src/models`).  Owner
  helpers: `RequireOnlyBlocks` (unknown blocks are fatal),
  `RequireNoAssignments`, `Flatten` (dotted model vocabularies) and
  `ToCanonicalText` (deterministic expanded text for broadcast/fingerprints).

It is C++17 and uses only the standard library, so model libraries that must
stay free of PIC/MPI/application headers can use it by adding
`-I<AMPS>/src/general/input`; their architecture checks allow exactly
`#include "amps_input.h"`.

## Users

| Code | Use |
|---|---|
| `srcMoon/MoonInput.cpp` | flat `application: moon` runtime block and Moon-owned value validation |
| `srcSEP3D/runtime/configuration_io.cpp` | complete srcSEP3D input (`kSrcSep3dLayout`) |
| `srcSEP/util/sep_initialization.cpp` | complete srcSEP input (`kSrcSepLayout`) |
| `src/models/sep_corona_swcme/shock_front/configuration.cpp` | reduced shock-front event (`physics: shock_front`) |
| `src/models/sep_corona_swcme/src/cme_event.cpp` | CME background event (`physics: cme_background`) |
| `src/models/sep_coronal_cme/src/configuration_parser.cpp` | coronal CME schema 5 (`physics: sep_coronal_cme`) |

## Test

```bash
make -C src/general/input test
```

`test_amps_input.cpp` checks, with hand-written fixtures: nested and named
blocks, comments, continuation and locations (INP-01); fourteen structural
errors with their `file:line` (INP-02); include resolution, include inside an
open block, the balanced-file rule, cycles and missing files (INP-03); and the
canonical-text round trip, flattening and the unknown-block policy (INP-04).
It prints `[INP-xx] PASS|FAIL message` lines and exits nonzero on failure; it
is the `amps-input` step of `srcSEP/test/run_all_test.sh` and
`srcSEP3D/test/run_all_test.sh`.
