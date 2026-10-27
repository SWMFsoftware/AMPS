# Hartley2 current-core compatibility status

Date: 2026-10-26

## Status

| Check | Result |
| --- | --- |
| Application preprocessing | Pass |
| Clean C++17/MPI compile and link | Pass, with legacy warnings |
| Runtime initialization | Not verified; required Hartley2 SPICE kernels are unavailable locally |
| Required AMPS core changes | None |

The `build` directory was removed before each compilation. The final clean build
produced the `amps` executable successfully.

## Application changes

### `input/hartley2.pl`

- Added the current working directory to Perl's module search path. Current Perl
  versions no longer include `.` in `@INC`, which prevented the AMPS
  preprocessing modules from loading.
- Replaced the hard-coded `srcTemp` working-source directory with `$ARGV[1]`,
  which is the directory supplied by the current AMPS preprocessor.

### `srcHartley2/CometData.cpp`

- Replaced the removed
  `PIC::Mesh::cDataCenterNode::totalAssociatedDataLength` member with
  `PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength`.
- This follows the current declaration in `src/pic/pic.h`; no compatibility
  member was added to the core.

### `srcHartley2/main.cpp`

- Updated the surface-intersection callback to accept the owning AMR tree node.
  The current `PIC::Mover::fProcessTriangleCutFaceIntersection` callback has five
  arguments and the mover supplies that node when it invokes the callback.
- Replaced the removed single-fluid Tecplot setters with a
  `PIC::CPLR::DATAFILE::TECPLOT::cIonFluidDescriptor`, then registered that
  descriptor in `IonFluidDescriptorTable`. This is the current multi-fluid
  interface and matches its use in other current AMPS applications.

### `srcHartley2/Comet.cpp`

- Initialized the dust-injection result to `false` and returned it outside the
  dust compile-time block. The function now returns a defined value when dust is
  disabled.

## Core assessment

No AMPS core source modification is required for Hartley2. The current core
already exposes every interface needed by the application:

- `PIC::Mesh::cDataCenterNode_static_data::totalAssociatedDataLength` owns the
  center-node static-data length.
- `PIC::Mover::fProcessTriangleCutFaceIntersection` includes the AMR node in its
  callback contract.
- `PIC::CPLR::DATAFILE::TECPLOT::cIonFluidDescriptor` and
  `IonFluidDescriptorTable` describe the loaded ion-fluid variables.

Adding aliases for the removed interfaces would preserve obsolete semantics.
In particular, a four-argument surface callback would discard the node passed
by the mover, and the removed Tecplot setters cannot describe the current
multi-fluid layout. Keeping the migration in `srcHartley2` makes the application
conform to the core contract.

Possible core improvements are optional and do not block Hartley2:

- Report an invalid SPICE toolkit or kernel path earlier during configuration.
- Provide a checked path-formatting helper for legacy application code that uses
  fixed-size buffers and `sprintf`.
- If many more legacy applications need migration, provide documented adapter
  examples for the current mover callback and Tecplot fluid descriptor APIs.

## Verification performed

The local SPICE toolkit was configured from
`/home/vtenishe/SPICE/cspice`. The compatibility work was then checked with:

```sh
./Config.pl -spice-path=/home/vtenishe/SPICE/cspice
rm -rf build
make -j8
perl -c input/hartley2.pl
git diff --check
```

The clean compile and link completed successfully. The resulting `amps` file is
a 64-bit ELF executable, and its dynamic library dependencies resolve locally.

The build still reports legacy warnings, including deprecated C++ `register`
qualifiers, one `%ld`/`int` format mismatch, integer narrowing, and possible
fixed-buffer truncation in `sprintf` calls. These warnings did not prevent the
build and are separate from the current-core API compatibility errors fixed
above.

## Runtime blockers and remaining configuration

An end-to-end Hartley2 run was not attempted because the mission kernel set
listed in `input/hartley2.input` is incomplete on this machine. Missing files
include the Hartley 2 ephemeris, frame, spacecraft-clock, and attitude kernels,
such as:

- `dif_epoch_nav_v1.bsp`
- `dif_dixi_nav_v1.bsp`
- `di_v18.tf`
- `di_hartley_2_v01.tpc`
- `dif_sclkscet_00084_science.tsc`
- `dif_sc_101104_highrate_v2.bc`

The actual Hartley2 kernel root must be supplied through `Config.pl` before a
runtime test. No unrelated Deep Impact kernels were substituted for the missing
files.

`srcHartley2/main.cpp` also contains a user-specific absolute path to a Tecplot
input file in the optional Lorentz-force/data-import branch. That path must be
made configurable before enabling that branch; the current Hartley2 input keeps
dust disabled.
