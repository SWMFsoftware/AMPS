# sep_common migration manifest

This change replaces the former top-level `sep_util` build wrapper with a
source-owning AMPS model sibling:

```text
AMPS/src/models/sep_common
```

## Canonical moves

The following `.h`/`.cpp` pairs moved from `srcSEP/util` to `sep_common`:

- `sep_transport_common`
- `sep_coefficient_physics`
- `sep_coefficient_registry`
- `sep_background_snapshot`
- `sep_test_registry`
- `sep_injection_spectrum`
- `sep_species_source`

`sep_common/makefile` now compiles these sources locally and produces
`sep_common.a`. The srcSEP and srcSEP3D makefiles insert those already-built
objects into their application `mainlib.a`, preserving AMPS's existing final
link contract without requiring another top-level linker argument.

There are deliberately no forwarding `.cpp` files. Dependency-light srcSEP
tests may compile the canonical sources directly with sanitizer flags, but no
application owns a private source copy. Includes use the explicit sibling
include path supplied by each makefile/test runner.

## SWCME relocation

The private `srcSEP/swcme` implementation and its temporary forwarding headers
were removed. srcSEP now includes public SWCME headers through `SWCME_DIR`
(normally `AMPS/src/models/swcme`), and its SWCME integration runner uses that
canonical tree. There is no application-local fallback or compatibility
directory. SWCME remains separate from `sep_common` because it is one
background/shock provider, not a prerequisite of general SEP transport.

## Additional stale-source cleanup

The uploaded srcSEP tree still contained `mover.cpp`, `fte_mover.cpp`, and
`fte_mover_dmumu.cpp` even though its Step 14 manifest, README, makefile, and
tests identify them as deleted legacy implementations. They were unlinked and
are excluded from this clean tree.

## Build compatibility correction

`srcSEP3D/core/sep3d_types.h` used a namespaced constant named `Pi`. AMPS
`general/constants.h` defines `Pi` as a global macro, so namespaces could not
prevent preprocessor substitution. The constant is now `kPi`, and test
`BLDL3D04` compiles the header after defining the AMPS macro to protect the
`cell_centered_linear_interpolation_cpp.cpp` include path.

The srcSEP3D makefile is copied by AMPS to `AMPS/build/main`. Its AMPS and
`sep_common` paths are now resolved as absolute paths from the active makefile
location instead of from the build process working directory. `BLDL3D05`
protects both the source and copied layouts.

## Compatibility policy

There is no `sep_util` or SWCME compatibility tree. Detached users can set
`SEP_COMMON_DIR` and `SWCME_DIR` while migrating; installed paths should use
the documented AMPS-root application layout and `src/models` shared-library
layout. Any application-local shared-model source or forwarding header is an
ownership error.
