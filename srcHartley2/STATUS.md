# Hartley2 1D/3D compatibility status

Date: 2026-10-27

## Current status

| Check | Result |
| --- | --- |
| Existing 3D application preprocessing | Pass |
| Existing 3D clean C++17/MPI compile and link | Pass, with legacy warnings |
| New spherical DIM1 clean C++17/MPI compile and link | Pass with the normal AMPS model libraries enabled |
| DIM1 radial mesh runtime | Pass on one MPI process |
| DIM1 geometric-profile unit test | Pass |
| Full 3D Hartley2 runtime | Not verified; the required mission SPICE kernels are unavailable locally |
| Full DIM1 coma simulation | Not implemented; the current DIM1 target generates and saves the radial mesh |

The `build` directory was removed before every compilation and recompilation.
The final repository configuration was restored to `hartley2.input`, and the
final clean build is the existing 3D Hartley2 executable.

## Spherical radial mesh

`input/hartley2.radial.input` is the DIM1 mesh configuration. The Hartley2
parser now accepts these values in its block:

- `Dimension=1` or `Dimension=3`;
- `NucleusRadius`, the inner radial coordinate in metres;
- exactly one of `DomainRadialExtent` (distance above the surface) and
  `DomainOuterRadius` (radius from the comet centre);
- `RadialCellNumber`, the requested number of cells in the ideal geometric
  profile; and
- `NucleusSurfaceResolution`, the width of the ideal first cell at the
  nucleus surface.

The two domain-size keywords are separate because “domain size” can mean
either radial extent or outer radius. The parser rejects a DIM1 input that does
not specify exactly one convention.

For surface width `dx0`, requested count `N`, radial extent `L`, and
expansion ratio `q`, `RadialMeshProfile.cpp` solves

```text
dx_i = dx0 q^i
L = dx0 (q^N - 1) / (q - 1).
```

`q=1` is used for the uniform limiting case. The input is rejected when
`L < N*dx0`, because satisfying that case would require cell widths to decrease
away from the nucleus instead of expand. The resolution callback returns the
width of the ideal geometric cell containing a requested radius.

The AMPS mesh is a dyadic block AMR tree: a refined block is split in half, and
all cells within a block have one width. It therefore approximates a general
geometric sequence conservatively; it cannot represent arbitrary `q` exactly.
The DIM1 driver reports both the requested ideal count and the number of cells
actually allocated by AMR. For the validation values in
`input/hartley2.radial.input`, the results are:

```text
nucleus radius                 500 m
outer radius                 10000 m
requested surface width         50 m
requested ideal cells            50
geometric width ratio       1.046969034872
ideal outer-cell width        473.9451 m
allocated AMR blocks                12
allocated AMR cells                 96
smallest allocated cell width  37.109375 m
```

An exact geometric grid with exactly `N` physical cells would require a
nonuniform radial coordinate mapping or a new non-AMR mesh representation.
That is a separate core design from the dyadic AMR refinement implemented
here.

## Hartley2 application changes

### `input/hartley2.pl`

- Added `Dimension=1|3` selection while preserving DIM3 as the default.
- Added parsing and validation for all radial mesh parameters.
- Configured DIM1 for spherical symmetry, one radial mesh dimension, no cut
  cells, and no NASTRAN boundary user data.
- Included the common exosphere surface-data definition required by current
  core boundary class declarations.
- Selected the implemented first-order spherical particle mover for DIM1. The
  second-order spherical mover path still contains unfinished runtime logic.
- Replaced the 3D Hartley2 object lists with `RadialMesh`,
  `RadialMeshProfile`, and `RadialMain` for a DIM1 target.
- Removed the 3D `Comet.h` include from the generated DIM1 `pic.h`.

### `srcHartley2/RadialMeshProfile.cpp` and `RadialMesh.h`

- Added validation and a bounded numerical solve for the geometric expansion
  ratio.
- Added ideal cell width, boundary, outer width, and local-resolution queries.
- Kept the geometric calculation independent of AMPS so it can be unit tested
  directly.

### `srcHartley2/RadialMesh.cpp`

- Added the AMPS mesh-resolution callback and input-backed radial parameters.
- Added diagnostics for the resolved domain, geometric profile, and allocated
  AMR cell count.

### `srcHartley2/RadialMain.cpp`

- Added a minimal DIM1 executable that initializes AMPS, builds the radial
  domain from the supplied nucleus radius to the resolved outer radius,
  distributes the mesh, and saves `hartley2.radial.mesh.msh`.

### `srcHartley2/makefile`

- Added the radial mesh/profile objects to the source object list.
- Permits the input parser to select the radial-only library and main program
  for `Dimension=1`.

### `srcHartley2/test_radial_mesh.cpp`

- Checks the solved ratio, final boundary, uniform limit, consecutive-width
  ratios, and invalid inputs.

The previously completed current-core migration remains in place: the parser
uses the current preprocessing directory, Hartley2 uses the current centre-node
static data length, its triangle callback accepts the owning AMR node, Tecplot
ion fluids use `cIonFluidDescriptor`, and dust injection returns a defined
value when dust is disabled.

## AMPS core changes required by DIM1

### `src/meshAMR/meshAMRgeneric.h`

- Renamed the computed face-neighbour lookup to `findNeibNodeFace` so it does
  not collide with stored face-neighbour arrays in DIM1 and DIM2.
- Implemented computed DIM1 face-neighbour lookup and DIM1 neighbour refinement
  limits.
- Retained the DIM1/DIM2 stored face-neighbour arrays needed by existing core
  paths when global 3D connection storage is disabled.
- Added DIM1 and DIM2 index bounds to centre-node signatures and associated-data
  save/load operations.
- Accepted planar or spherical symmetry for DIM1 and planar or axial symmetry
  for DIM2 during mesh initialization.
- Implemented spherical DIM1 face areas as `4*pi*r^2` and corrected the DIM1
  face-area multiplier array length.
- Made mesh-refinement diagnostics iterate over the compiled dimension, which
  removes out-of-bounds Y/Z reads in DIM1.

The existing cell-measure path already computes a spherical shell volume as
`4*pi*(r_max^3-r_min^3)/3`; no change was needed there.

### `src/pic/pic.h`, `pic_parallel.cpp`, and `parallel/init_send_halo.cpp`

- Replaced a hard-coded `cMeshAMR3d` particle-table type with `cAmpsMesh`.
- Replaced obsolete dimension-specific nested AMR node-ID types with the
  current global `cAMRnodeID`.
- Updated halo construction to use the renamed computed face-neighbour lookup.

### `src/pic/pic_mover.cpp`

- Made body-of-rotation boundary declarations available to dimension-neutral
  switch code.
- Initialized the DIM1 boundary-face index in the second-order mover and
  replaced stale `x`/`v` identifiers with that mover's `xInit`/`vMiddle`
  variables.
- Corrected three spherical rotation expressions that copied the Y coordinate
  into the Z component.

The radial configuration deliberately selects the first-order mover because
these compile fixes do not complete the second-order spherical algorithm.

### `src/pic/pic_field_line.cpp`

- Kept physical Parker-spiral position, field, and direction vectors at three
  components even when the computational mesh dimension is one.

## Verification performed

The DIM1 executable was built with all standard AMPS model libraries enabled,
so the radial input does not leave persistent model switches that can break a
later 3D build. Verification used:

```sh
rm -rf build
make -j8                         # Dimension=1 configuration
/data/vtenishe/COMET1D/AMPS/amps # run from /tmp, one MPI process

rm -rf build
make -j8                         # restored Dimension=3 configuration

g++ -std=c++17 -Wall -Wextra -pedantic \
  srcHartley2/test_radial_mesh.cpp srcHartley2/RadialMeshProfile.cpp \
  -o /tmp/hartley2_radial_mesh_test
/tmp/hartley2_radial_mesh_test

perl -c input/hartley2.pl
git diff --check
```

The DIM1 runtime generated a 4,049-byte mesh file and completed mesh
distribution. Direct execution inside the filesystem sandbox could not open
Open MPI's local communication socket; the same executable completed when run
with normal local socket access.

Both builds retain legacy warnings in old Hartley2 code, including deprecated
`register` qualifiers, `%ld`/`int` format mismatches, narrowing conversions,
and possible fixed-buffer `sprintf` overflows. These warnings did not prevent
either link.

## Remaining work for a DIM1 coma simulation

The delivered DIM1 path is a verified radial mesh generator. Turning it into a
coma simulation still requires physical inputs that were not supplied:

- nucleus production rate or surface flux by species;
- source temperature and bulk radial velocity distribution;
- outer-boundary particle behavior and loss processes;
- timestep and particle-weight policy for the requested physical regime; and
- requested sampling/output variables.

Those values should not be inferred from the existing irregular 3D Hartley2
surface model because spherical symmetry changes the source representation.

The full 3D runtime also remains blocked by missing Hartley2 ephemeris, frame,
spacecraft-clock, and attitude SPICE kernels. No unrelated mission kernels were
substituted.
