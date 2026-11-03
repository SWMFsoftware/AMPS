# Mars Hot-Oxygen Corona Development Status

Last updated: 2026-11-03 (America/Chicago)

## Objective

Develop and validate the AMPS Mars hot-oxygen corona application in `srcMars/`.
The specific physics or implementation change for the current effort has not yet
been defined; this file is the persistent handoff and development log.

## Working Constraints

- All reads, writes, builds, tests, and generated artifacts must remain inside
  the current AMPS repository (`/data/vtenishe/Mars-test/AMPS`) or its
  subdirectories.
- Preserve unrelated and pre-existing user changes.
- Update this file after each material code, configuration, build, test, or
  design change.
- Record exact verification commands and outcomes, including failures and known
  limitations.

## Repository Baseline

- Branch: `master`
- Starting commit observed: `2369abe`
- Baseline inspection date: 2026-11-03
- Build state at baseline: not tested before repair work.

The following changes already existed when tracking began and are treated as
user-owned/inherited work.

### Modified tracked files

- `srcMars/BackgroundAtmosphere.cpp`
- `srcMars/DifferentialCrossSection.h`
- `srcMars/MTGCM.h`
- `srcMars/Mars.cpp`
- `srcMars/Mars.dfn`
- `srcMars/Mars.h`
- `srcMars/MarsBackgroundAtmosphereFox.h`
- `srcMars/m-gitm.h`
- `srcMars/main.cpp`
- `srcMars/main_lib.cpp`

At baseline these tracked changes totaled approximately 540 insertions and 504
deletions. They have not been comprehensively reviewed; this session addressed
only the blockers encountered while compiling the current oxygen-only case.

### Untracked files and directories

- `srcMars/DifferentialCrossSection_O-CO2-Edependent.h`
- `srcMars/DifferentialCrossSection_O-CO2.h`
- `srcMars/DifferentialCrossSection_O-CO2_singleCollE_2_5eV.h`
- `srcMars/NewOCO2-2.5eV-v03.h`
- `srcMars/NewOCO2-allE-total-smoothed.h`
- `srcMars/NewOCO2-table_E.h`
- `srcMars/newestsources/` (alternate source tree containing Mars application
  sources and a makefile)

## Current State

- Workspace scope and status-tracking requirements are established.
- Existing `srcMars` modifications were inventoried before build repairs began.
- The configured application is oxygen-only (`SpeciesList=O`; generated
  `_O_SPEC_=0`).
- A clean `make -j` completes and produces the root-level `amps` executable.
- The generated executable is a dynamically linked, non-stripped, 64-bit x86-64
  ELF file. At the final check it was 2,258,968 bytes.
- The 26 M-GITM path-construction `-Wformat-overflow` warnings and both obsolete
  configuration warnings have been removed.
- No simulation/runtime validation has been performed yet.

## Build Procedure

From the repository root:

```bash
rm -rf /data/vtenishe/Mars-test/AMPS/build
make -j
```

The clean removal is required before compilation so that AMPS regenerates the
working source tree from the current configuration and `srcMars` sources.

## Changes Made During Build Repair

### `mars.pl`

- Added the current repository directory to Perl's module search path using
  `Cwd` and `use lib cwd`.
- This allows plain `make -j` to locate the local `ampsConfigLib.pm`; an external
  `PERL5LIB=.` environment override is no longer required.
- The root-level `mars.pl` was already present but is currently untracked by Git.

### `srcMars/Mars.cpp`

- Updated volume-injection associated-data access from the obsolete
  `PIC::Mesh::cDataCenterNode::LocalParticleVolumeInjectionRateOffset` owner to
  `PIC::Mesh::cDataCenterNode_static_data::LocalParticleVolumeInjectionRateOffset`.
- Compiled hot-carbon integration cases only when `_C_SPEC_ >= 0`.
- Compiled the legacy hot-carbon particle-production implementation only when
  carbon is configured; the oxygen-only build receives a zero-return stub.
- Restored pointer member access (`->`) for the current AMPS mesh API.

### `srcMars/main.cpp`

- Restored pointer member access (`PIC::Mesh::mesh->...`) throughout mesh setup,
  allocation, distribution, diagnostics, and sampling code. The current AMPS
  API declares `PIC::Mesh::mesh` as a pointer.

## Known Warnings and Limitations

- OpenMPI repeatedly prints `opal_ifinit: socket() failed with errno=1` in this
  restricted environment for compiler-wrapper invocations. Compilation and
  linking still complete with exit 0; this is environment noise rather than a
  Mars source warning.
- The build verifies compilation and linkage only; it does not establish that
  the model initializes successfully with its input datasets or produces
  physically valid hot-oxygen corona results.

## Warning Remediation

The warning fix was implemented and verified on 2026-11-03:

1. Added `newMars::GetMGITMInputFilePath()` in `Mars.h`. It constructs the full
   filename as a `std::string` once per loading routine.
2. Replaced all 26 active fixed-buffer `sprintf` calls with reuse of that path:
   5 sites in `BackgroundAtmosphere.cpp`, 12 in `Mars.h`, and 9 in `main.cpp`.
3. Hardened both dataset readers (`MTGCM.h` and `m-gitm.h`) by rejecting null or
   overlong filenames before copying into their 2000-byte `DataFileName` arrays.
   Checked `strlen` plus `memcpy` avoids both overflow and silent truncation.
4. Removed deprecated `WorkingSourceDirectory=srcTemp` from `mars.input`.
5. Replaced the obsolete `#UserDefinitions` block with supported explicit
   definitions in `#General` and `#Mesh`. A generated-file comparison caught
   that simply removing the obsolete block changed PIC and mesh modes to off,
   despite the configurator's claim that the block had no effect. The final
   generated configuration again has PIC user definitions enabled and mesh user
   definitions on, matching the pre-fix behavior.

The required clean `make -j` completes without the former compiler or
configuration warnings. The OpenMPI socket diagnostic remains environment-only.

## Verification Log

| Date | Command/check | Result |
| --- | --- | --- |
| 2026-11-03 | `git status --short -- srcMars` | Confirmed inherited modified and untracked files listed above. |
| 2026-11-03 | `git diff --stat -- srcMars` | Confirmed 10 modified tracked files; no interpretation or validation performed. |
| 2026-11-03 | Remove `build/`; `make -j` | Failed during configuration because `mars.pl` could not locate local `ampsConfigLib.pm`. |
| 2026-11-03 | Remove `build/`; `PERL5LIB=. make -j` | Reached Mars compilation; exposed obsolete center-node offset ownership and unguarded hot-carbon references. |
| 2026-11-03 | Remove `build/`; `PERL5LIB=. make -j` after first repair | Reached `main.cpp`; exposed object-style access to pointer-valued `PIC::Mesh::mesh`. |
| 2026-11-03 | Remove `build/`; `PERL5LIB=. make -j` after mesh repair | Passed compilation and final link; exit 0. |
| 2026-11-03 | Remove `build/`; plain `make -j` after `mars.pl` repair | Passed clean configuration, compilation, and final link; exit 0. |
| 2026-11-03 | `file amps` and `stat` | Verified executable 64-bit x86-64 ELF, 2,258,744 bytes. |
| 2026-11-03 | Remove `build/`; `make -s -j` warning investigation | Passed configuration, compilation, and final link; exit 0. Found one repeated C++ warning class (26 unsafe path-formatting sites), two obsolete configuration constructs, and one restricted-environment OpenMPI diagnostic. |
| 2026-11-03 | Remove `build/`; `make -j` after initial warning fix | Passed with exit 0 and removed the compiler/configuration warnings. Generated-mode comparison found that direct removal of the obsolete block had disabled PIC and mesh user-definition modes. |
| 2026-11-03 | Remove `build/`; `make -j` after supported mode definitions | Passed configuration, compilation, and final link with exit 0. Former `-Wformat-overflow` and obsolete-configuration warnings were absent; PIC user definitions remained enabled and mesh user definitions remained on. The OpenMPI environment diagnostic remained. |
| 2026-11-03 | `file amps` and `stat` after warning fix | Verified executable 64-bit x86-64 ELF, 2,258,968 bytes. |

## Open Questions

- What concrete physical model, input data, defect, or simulation case should be
  addressed first?
- Is `srcMars/newestsources/` a reference snapshot, a migration target, or the
  intended active implementation?
- What runtime case and physical diagnostics define acceptance for the current
  hot-oxygen corona work?

## Development Log

### 2026-11-03

- Created this status file at the user's request.
- Established the repository boundary and captured the pre-existing `srcMars`
  working-tree state.
- Removed and regenerated `build/` for each full compilation attempt.
- Repaired Perl module discovery, two AMPS API compatibility issues, and
  oxygen-only conditional compilation.
- Completed a clean parallel build with the exact requested `make -j` command.
- Produced the root-level `amps` executable; did not launch a simulation.
- Repeated a clean quiet parallel build to classify warnings. Documented the
  exact warning sources and a prioritized repair proposal; no warning fixes were
  applied during this diagnostic pass.
- Implemented centralized safe M-GITM path construction and checked dataset
  filename storage, removing all 26 format-overflow warnings.
- Removed obsolete configuration constructs while preserving the previously
  generated PIC and mesh user-definition modes through supported definitions.
- Completed the required clean parallel build after the final fix; only the
  restricted-environment OpenMPI socket diagnostic remains.

