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

The [`parallel_diffusion`](parallel_diffusion/README.md) library is a
dependency-free sibling within this ownership boundary. Its standalone
revision-1.4 PD00–PD10 archive implements all 16 first-release identifiers,
including spectral, pitch-angle, nonlinear, adapter, table and batch paths.
It remains built in its own subdirectory. srcSEP3D now consumes that
canonical archive through its approved Parker-only D13--D16 binding and
flattens the two objects into the AMPS application archive; srcSEP remains
unbound. This keeps the established seven-member `sep_common.a` and both
applications' existing exact-member audits unchanged until promotion is
deliberately coordinated. The application incorporation record, component-test
evidence, derivative limitations and native-qualification boundary are
recorded beside the library.

The [`perpendicular_diffusion`](perpendicular_diffusion/README.md) library is
the corresponding revision-2.1 host-neutral coefficient/statistics sibling.
It implements every fully specified registry identity, delegates coupled
NLGCE pairs to the parallel library's single numerical owner, and preserves
incomplete literature closures as typed source gates. It is also built and
verified in its own directory: neither application consumes it yet, so its
standalone PASS does not alter the established seven-member `sep_common.a` or
claim Parker-mover/native/MPI qualification. Its integration boundary and
remaining transport work are recorded beside the library.

The [`mean_free_path`](mean_free_path/README.md) library is the host-neutral
focused-transport coefficient layer for `MEAN_FREE_PATH_MODEL.md` version 1.2.
It supplies typed mean-free-path kinds, pitch-angle shapes, SEP and GCR
prescriptions, shock closures, strict parser-neutral selection, provenance and
source-incomplete gates. Its current delivery is standalone: `srcSEP3D` does
not call it yet, and its
[`INTEGRATION_PLAN.md`](mean_free_path/INTEGRATION_PLAN.md) records the future
focused-mover binding. Standalone coefficient tests are not particle-mover,
native, MPI or event validation.

The [`pitch_angle_diffusion`](pitch_angle_diffusion/README.md) library is the
host-neutral Release A implementation of `PITCH_ANGLE_SCATTERING_MODEL.md`
v1.2: prescribed D_μμ shapes (PA-ISO, PA-QFORM, PA-EPS, PA-EPREM) normalized
per particle to a λ∥ law, three scattering step rules, the pure scattering
kernels with an independent finite-volume reference operator, and a
coefficient-kernel certification matrix derived from its statistical
acceptance campaign.  srcSEP3D's focused-diffusion mover consumes it through
`pitch_angle_diffusion_model = pitch-angle-library` and flattens its three
objects into the application archive; srcSEP remains unbound.  Its stage
ledger (`IMPLEMENTATION_STATUS.md`) separates library, component and native
evidence; the canonical data archive (roadmap PS01) is not available.

## Shared runtime-input vocabulary (`sep_input_layout.h`)

`sep_input_layout.h` is the single table of the shared SEP runtime-input
vocabulary, content schema 6 (normative specification:
[`../../INPUT_FILE.md`](../../INPUT_FILE.md) section 5).  It lists every block
(`kBlocks`) and every application-owned key (`kKeys`) that srcSEP or srcSEP3D
understands, each marked with the application(s) that implement it
(`kSrcSep`, `kSrcSep3d`, `kBoth`), and supplies the helpers that encode the
decisions D1-D5: `kSchemaVersion = 6` (D1), the typed
`NotSupportedMessage()` "... is a srcSEP3D setting that is not supported by
srcSEP (shared SEP input schema 6, src/INPUT_FILE.md section 5)" diagnostic
(D3), `SwcmeOwnedParkerKey()` for the seven Parker ambient keys owned by
`physics: swcme`, `SwcmeSourceKeys()` and `SynthesizeSwcmeSource()` (exact
round-trip decimal `source.*` assignments in J/kg/si from `source: shock`,
D4), and `SwcmeThreeDimensionalKeys()` (SWCME3D-only keys that srcSEP's
SWCME1D rejects).

The header is header-only and depends only on the C++ standard library, so a
parser can validate a complete input before AMPS, MPI or a mesh exists; it is
not part of `sep_common.a` and does not change the seven-member archive.
Consumers: `srcSEP3D/runtime/configuration_io.cpp` and
`srcSEP/util/sep_initialization.cpp` walk their documents against the table,
and the srcSEP3D stage-1 test `CFG3D18` cross-checks
`SwcmeThreeDimensionalKeys()` against
`swcme::input3d::detail::geometry_keys()`.  Validation: srcSEP3D stage-1
`CFG3D` group (including `CFG3D01` and `CFG3D18`) and srcSEP `INIT01`-`INIT03`
pass; the full native regressions with the schema-6 decks are pending.

## Model siblings and boundaries

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
