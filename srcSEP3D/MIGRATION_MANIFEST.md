# srcSEP3D Migration Manifest

This manifest records what Phase R0 removed, retained, or replaced. It is an
evidence document: completion claims must correspond to an executable test or
an explicitly identified external build gate.

## Phase R0 — Rebaseline the production tree

### Deleted production files

| File | Disposition | Reason |
|---|---|---|
| `SEP3D.cpp` | deleted | contained the axisymmetric y=0 mover and legacy global energy sampler; neither belongs in a three-dimensional transport application |
| `amps/amps_sampling.h` | deleted | declaration-only placeholder with no production implementation; retaining it implied sampling existed when it did not |
| `core/sep3d_energy_distribution.{h,cpp}` | deleted | unintegrated replacement kernel had no registered acceptance test or production consumer; Phase O will add sampling as one complete implementation |
| `test/stage1` | removed from deliverable | generated executable from the uploaded archive, not source; rebuilt locally by the runner |

### Retained and rewritten files

| File | R0 state |
|---|---|
| `main_lib.cpp` | retains only AMPS/Exosphere link hooks and explicit execution guards; wedge mesh, fixed resolution, Maxwellian prepopulation, placeholder sampling, and placeholder outputs are removed |
| `main.cpp` | clean AMPS driver source retained for the `main.a` archive; strict-warning clean by construction |
| `SEP3D.h` | declares only current application hooks and includes the AMPS mover-status adapter; no retired mover or sampler declarations |
| `core/sep3d_types.h` | replaces integer-valued mover return codes with semantic `ParticleMotionOutcome` |
| `makefile` | truthful production manifest, enclosing-AMPS production build, archive symbol audit, AMPS-independent test targets, and absolute path discovery valid from both `srcSEP3D` and copied `build/main` locations |
| `test/run_tests.py` | rebuilt around the selector vocabulary used by `srcSEP` and structured JSON/JUnit evidence |

### Mover return-code correction

The prototype core assigned local integers to mover results. They did not match
the AMPS definitions in `pic.h`: AMPS uses deleted-on-face `0`, left-domain
`2`, and motion-finished `3`. A model extension also occupied a value used by
AMPS. This is corrected as follows:

1. L0 returns `Core::ParticleMotionOutcome`, whose enumerators are semantic and
   are never used as AMPS integers.
2. `amps/amps_mover_status.h` is the only conversion boundary.
3. The adapter maps `Advanced` to `_PARTICLE_MOTION_FINISHED_` and
   `LeftDomain` to `_PARTICLE_LEFT_THE_DOMAIN_`.
4. `static_assert` binds the adapter to the actual AMPS ABI at production
   compile time.
5. `BLDL3D03` also reads the supplied AMPS `pic.h` and records the observed
   values in the test evidence.

Invalid background, underflow, and internal error statuses are not mapped to a
valid AMPS particle result. They must be handled before the adapter is called.

### Production execution policy

R0 does not provide a simulation. Calling `amps_init_mesh`, `amps_init`,
`amps_time_step`, or `localResolution` terminates with an explicit message.
This is safer and more truthful than running the previous wedge or silently
substituting placeholder values. The guard is removed only when R2 and the
mesh phase supply a validated Runtime and production domain.

### R0 evidence

| Gate | Implementation | Source-only behavior |
|---|---|---|
| `BLDL3D01` | delegates through `make strict-production` to the enclosing top-level `make amps`, then audits `build/main` archives | SKIP when a real `Makefile.conf` is unavailable |
| `BLDL3D02` | production-tree and makefile-manifest inspection | executable without AMPS |
| `BLDL3D03` | adapter contract plus real `pic.h` macro inspection | executable with an AMPS source tree; SKIP if `pic.h` is unavailable |
| `BLDL3D04` | compile-only AMPS `Pi`-macro collision regression | executable without a configured AMPS build |
| `BLDL3D05` | source/copy makefile path-resolution fixture | executable without AMPS libraries; invokes both makefiles from an unrelated working directory |
| `LAY01`, `LAY02`, `BLD01` | three-layer AMPS-dependency guard | executable without AMPS/MPI |

## Earlier Step 1–4 material: corrected status

The uploaded tree described Steps 1–4 as complete. R0 corrects that claim:

- the directory layout and standalone registry are useful and retained;
- the original production source had not actually retired the axisymmetric
  mover until this R0 change;
- `sep_common` now owns the canonical shared sources and archive beside both
  applications; the former `srcSEP/util` copies are absent;
- SWCME now has one canonical `src/models/swcme` tree and srcSEP consumes its public
  headers from that location; the former `srcSEP/swcme` duplicate is absent;
- a configured enclosing AMPS executable link is still required before the
  production-build portion of R1 can be declared complete;
- no complete three-dimensional coupled background snapshot owner exists yet.

No later physics phase should use the old Step 1–4 completion labels as release
evidence. The current runner reports the gates that are actually executable.

## AMPS macro compatibility correction

`SEP3D::Core::Const::Pi` was renamed to `kPi`. AMPS exposes `Pi` as a global
macro in `general/constants.h`; after `pic.h` includes the application header,
the preprocessor rewrote the namespaced declaration and caused
`cell_centered_linear_interpolation_cpp.cpp` to fail with “expected
unqualified-id before numeric constant.” `BLDL3D04` now reproduces this include
order in a compile-only regression.

## Copied-build makefile path correction

AMPS copies `srcSEP3D` into `AMPS/build` and renames that application directory
to `AMPS/build/main`. A literal `../Makefile.conf` is therefore wrong in the
copied tree: it resolves to `AMPS/build/Makefile.conf` instead of the canonical
`AMPS/Makefile.conf`.

The makefile now obtains the absolute directory containing the active makefile
from `MAKEFILE_LIST`, identifies whether it is the source or copied location,
and resolves `AMPS_ROOT` once. `AMPS_CONFIG`, `SEP_COMMON_DIR`, the common
archive, and all common-object paths are absolute descendants of that root.
This logic does not depend on the current working directory. `BLDL3D05`
constructs both layouts, invokes each makefile from a third directory, and
requires identical canonical paths.

## Enclosing production-build correction

The application makefile is designed to compile production objects after AMPS
has copied the application to `AMPS/build/main`. A direct production submake in
`AMPS/srcSEP3D` lacks the include flags for generated headers such as
`AMPS/build/pic/pic.h`; such a compile is not equivalent to the production
build. `strict-production` and `BLDL3D01` therefore delegate to the top-level
`make amps` target and audit the resulting absolute `build/main` archive paths.

## Forbidden production identifiers and operations

`BLDL3D02` rejects the following in production sources:

- `Mover_Axisymmetric_SecondOrder`
- `TotalParticleAcceleration`
- `GlobalEnergyDistribution`
- `inject_particle_onto_field_line`
- `CMPI_channel`
- the former wedge bounds `8.760e+08` and `9.445e+08`
- `PrepopulateDomain`, placeholder mesh output, or mesh-file creation in the
  R0 application boundary

Documentation and tests may name retired identifiers when explaining or
detecting them; production headers and translation units may not.
