# Stage 0: configuration and architecture contract

Stage 0 implements no magnetic or plasma provider. Its purpose is to make an
invalid physical choice impossible to reinterpret silently before later stages
allocate a field, a mesh, or particles.

The parser first lexes assignments and locates the unique
`run.schema_version`. Schemas 1--4 are returned byte-for-byte to their frozen
application parser. Schema 5 is then parsed using the generated normative
registry. Unknown, duplicate, incomplete, unit-bearing, non-finite, and
inconsistent branch values are errors.

Numerical keys encode their SI dimension in the name and accept bare SI
values. The typed record materializes Stage 0--2 authorities, while the full
registry freezes later selectors. An unavailable selection returns typed
`NotImplemented`; it never substitutes a different physical model.

Restart identity is SHA-256 over sorted, length-prefixed normalized records.
Path spelling and comments are excluded, but content checksums, selectors,
radii, frames, EOS, interface policy, campaign seed, and numerical authorities
remain included. A mismatch fails before any provider is published.

The library depends only on C++17 and neutral `sep_common` status records.
AMPS/PIC, MPI, Tecplot, applications, and `swcme` are forbidden dependencies.
`ARCHSCCM01` checks source includes, archive symbols, and an external public-
header consumer compilation.
