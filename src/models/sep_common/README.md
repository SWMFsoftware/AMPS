# sep_common

`sep_common` owns the dependency-free SEP kernels shared by the one-dimensional
`srcSEP` model and the three-dimensional `srcSEP3D` model. It is an AMPS model
sibling, not an SWCME submodule: most of its APIs are transport, coefficient,
injection, species, background-snapshot, or test-infrastructure APIs that are
also valid for Parker-spiral and other non-SWCME backgrounds.

## Source-distribution ownership

`SOURCE_MANIFEST.json` declares this directory as the sole source owner for
the seven shared SEP kernels. `build/` and `sep_common.a` are generated and
must not appear in a source release. The AMPS-level B01 hygiene gate checks
that rule from strict allowlists shared with both applications and SWCME.
Generated and retired paths take precedence, and an unclassified path fails
the release. Its isolated self-test inserts stale objects, caches, output,
retired, and unclassified controls; its deterministic archive mode then
reopens and verifies the exact member set. The archive is rebuilt from these
sources in every clean extraction; a preexisting archive is never accepted as
evidence that the current sources compile.

## Ownership boundary

The directory contains the canonical source and header for each shared kernel:

| Kernel | Responsibility |
| --- | --- |
| `sep_transport_common` | Relativistic conversions, timestep selection, keyed random streams, and transport status types |
| `sep_coefficient_physics` | Gyro/rigidity, pitch-angle diffusion, mean-free-path, and spatial-diffusion calculations |
| `sep_coefficient_registry` | Named coefficient-provider selection and configuration fingerprints |
| `sep_background_snapshot` | Immutable prepared background data exchanged at model boundaries |
| `sep_injection_spectrum` | Dependency-free injection-spectrum functions |
| `sep_species_source` | Species/abundance/source properties |
| `sep_test_registry` | Shared test descriptors, result records, and JSON/JUnit writers |

SWCME remains a separate sibling because it supplies one possible solar-wind,
CME, and shock-background provider. `sep_common` must never include `pic.h`,
`mpi.h`, or an SWCME header. Coupling adapters belong in `srcSEP`, `srcSEP3D`,
or `swcme`, depending on which side owns the API.

The former copies under `srcSEP/util` were removed. Do not add forwarding
`.cpp` files there: compiling the same definitions under application-specific
flags defeats the shared-library boundary and can create duplicate symbols.

## Build and verification

From `AMPS/src/models/sep_common`:

```sh
make
make verify
make clean
```

The build produces `sep_common.a` and keeps intermediate objects under
`build/`. Both applications include this directory for headers. Production
archives consume the objects built here; dependency-light sanitizer tests may
compile the canonical sources directly with instrumentation, but must never
maintain private source copies.

The numerical build deliberately excludes `-ffast-math`. `srcSEP3D` test
`UTIL02` freezes representative kernel results byte-for-byte and will detect a
floating-point flag or implementation change.

Both applications audit this archive without inspecting each other's trees.
For `srcSEP`, `make test-sep-common-ownership-unit` verifies the seven-member
archive, unique strong definitions, canonical-header consumer link, and
source-versus-`build/main` path equivalence.  `srcSEP3D` performs its own
archive and byte-exact frozen-kernel checks.  The AMPS-level source-package
gate is the only place that may inspect both independent applications at once.
