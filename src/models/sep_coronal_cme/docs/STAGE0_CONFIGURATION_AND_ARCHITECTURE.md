# Stage 0: configuration and architecture contract

Stage 0 implements no magnetic or plasma provider. Its purpose is to make an
invalid physical choice impossible to reinterpret silently before later stages
allocate a field, a mesh, or particles.

The deck is written in the generic AMPS runtime input format
(`src/INPUT_FILE.md`) and read by the header-only, standard-library-only
`src/general/input/amps_input.h` (the model makefile adds
`-I../../general/input`; `tools/check_architecture.py` allows exactly
`#include "amps_input.h"` as the one AMPS include). The complete configuration
is one `#begin physics: sep_coronal_cme` ... `#end physics` block whose child
blocks are the schema sections named as in the registry (`#begin pfss`; nested
blocks for dotted section names such as `shock_geometry.center_distance`;
`#begin observer: earth` for an identified `observer.ID` section).
Assignments directly in the model block and dotted keys are rejected, and
comments use `!`. Only `run.schema_version = 5` is accepted: the former
pass-through of schema-1--4 bytes to the frozen application parser was
removed, because those schemas are `srcSEP3D` application decks that
`srcSEP3D` now reads itself. Schema 5 is parsed using the generated normative
registry. Unknown, duplicate, incomplete, unit-bearing, non-finite, and
inconsistent branch values are errors. The Stage-0 fixture
`test/data/schema5_stage0.in` is regenerated in this format by
`tools/generate_stage0_fixture.py`.

Numerical keys encode their SI dimension in the name and accept bare SI
values. The typed record materializes Stage 0--2 authorities, while the full
registry freezes later selectors. An unavailable selection returns typed
`NotImplemented`; it never substitutes a different physical model.

Restart identity is SHA-256 over sorted, length-prefixed normalized records.
Path spelling and comments are excluded, but content checksums, selectors,
radii, frames, EOS, interface policy, campaign seed, and numerical authorities
remain included. A mismatch fails before any provider is published.

The library depends only on C++17, neutral `sep_common` status records, and
the header-only generic input reader `amps_input.h` (itself standard-library
only, with no PIC or MPI). AMPS/PIC, MPI, Tecplot, applications, and `swcme`
are forbidden dependencies.
`ARCHSCCM01` checks source includes, archive symbols, and an external public-
header consumer compilation.
