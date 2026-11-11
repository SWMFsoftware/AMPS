# AMPS core mesh tests

This configured AMPS application validates the native AMR tree, its adaptation
metadata, cell/corner interpolation, and deterministic cut-cell geometry.  It
is intentionally separate from `srcCoreUnitTests`: mesh dimensionality is a
compile-time property, so the same source is configured and built in three
profiles without weakening the global uniqueness of test IDs.

## Profiles and configuration

| Application | Compile-time macro | Geometry | Registered scope |
|---|---:|---|---|
| `test/core-mesh-1d` | `DIM=1` | planar line, unit length | binary topology and adaptation |
| `test/core-mesh-2d` | `DIM=2` | planar square, unit area | quadtree topology and adaptation |
| `test/core-mesh-3d` | `DIM=3` | Cartesian cube, unit volume | octree topology/adaptation, interpolation, cut cells |

The decks are `input/test/core-mesh-{1d,2d,3d}.input`.  They set `DIM` with
`defineLocation(pic/picGlobal.dfn)`; no manual header edit is part of this test.
All profiles use two cells and one ghost cell in each represented block
direction.  The 3-D deck enables internal-surface data types for the synthetic
cut-cell fixture, but installs no exosphere source, target, SPICE geometry, or
external NASTRAN file.

Each executable uses a deterministic unit domain.  The root is split once and
its all-lower child is split once more.  In `D` dimensions the independent
oracles are therefore

```
branching factor = 2^D
all nodes        = 1 + 2*2^D
leaves           = 2*2^D - 1
leaf volume      = 2^(-D*level)
sum leaf volume  = 1
```

Bounds are dyadic (`0`, `1/4`, `1/2`, `1`) and exactly representable in IEEE
binary arithmetic.  Topology checks therefore use exact comparisons; only
accumulated volume and interpolation equations use roundoff-scaled tolerances.

## Stable case IDs

Every dimension registers the following profile-qualified IDs (replace `D`
with `1D`, `2D`, or `3D`):

- `CORE-MESH-D-TOPOLOGY-UNIFORM-001`
- `CORE-MESH-D-TOPOLOGY-NONUNIFORM-001`
- `CORE-MESH-D-SEARCH-001`
- `CORE-MESH-D-PARENT-CHILD-001`
- `CORE-MESH-D-NEIGHBORS-001`
- `CORE-MESH-D-CHECKSUM-001`
- `CORE-MESH-D-ADAPT-REFINE-001`
- `CORE-MESH-D-ADAPT-COARSEN-001`
- `CORE-MESH-D-ADAPT-DATA-001`
- `CORE-MESH-D-ADAPT-VOLUME-001`
- `CORE-MESH-D-ADAPT-COUNTERS-001`
- `CORE-MESH-D-ADAPT-LOAD-001`
- `CORE-MESH-D-ADAPT-MAX-LEVEL-001`
- `CORE-MESH-D-ADAPT-CYCLE-001`

The 3-D profile additionally registers:

- `CORE-MESH-3D-INTERPOLATION-UNIFORM-001`
- `CORE-MESH-3D-INTERPOLATION-BOUNDARY-001`
- `CORE-MESH-3D-INTERPOLATION-AMR-001`
- `CORE-MESH-3D-INTERPOLATION-LEGACY-MASK-001`
- `CORE-MESH-3D-INTERPOLATION-DUPLICATE-001`
- `CORE-MESH-3D-INTERPOLATION-PERMUTATION-001`
- `CORE-MESH-3D-INTERPOLATION-CACHE-001`
- `CORE-MESH-3D-INTERPOLATION-INPUT-IMMUTABLE-001`
- `MESH-CUTCELL-TRIANGLE-METRICS-001`
- `MESH-CUTCELL-RAY-001`
- `MESH-CUTCELL-TANGENT-DEGENERATE-001`
- `MESH-CUTCELL-CLOSEST-001`
- `MESH-CUTCELL-PLANE-VOLUME-001`
- `MESH-CUTCELL-SPHERE-VOLUME-001`
- `MESH-CUTCELL-SPHERE-AREA-CONVERGENCE-001`
- `MESH-CUTCELL-NASTRAN-ROUNDTRIP-001`

Coarsen, data-transfer, and split/coarsen-cycle registrations are maintained in
their dedicated test object because those cases exercise an explicit production
restriction/prolongation policy; they must never be replaced by no-op PASS
receipts.

## Numerical contracts

The topology cases validate physical search at every leaf midpoint, reciprocal
same-level face connectivity, parent/child refinement and integer geometry,
unique production AMR IDs, and an allocation-order-independent canonical
checksum. The CHECKSUM case also writes and reads an actual native mesh file
in a private mesh, then checks the exact unit-cube split, stable node IDs,
parent links, midpoint ownership and reciprocal face geometry. This covers
lower-dimensional face-table initialization and serialization, which an
in-memory hash cannot exercise. The temporary file belongs to the case's
isolated run directory and is removed after a successful check. In every
dimension the uniform case also creates a private,
topology-only mesh and exercises production `buildMesh` with bounded
resolution, requiring the root split, the expected generation/maximum level,
and a unique `2^D`-node bottom list. Adaptation cases validate accepted split
counts, generation counters, leaf-only load callback application, and typed maximum-level
rejection with no generation change.

Corner interpolation is checked against the affine field

```
f(x,y,z) = 1.25 + 2*x - 3*y + 0.5*z,
```

which trilinear interpolation must reproduce exactly up to floating-point
roundoff.  Boundary tests exercise the half-open final-cell convention.  The
AMR row test requires nonnegative, unit-sum support from both sides of a 2:1
interface.  The legacy-mask test independently enumerates all 256 choices of
refined root octants. For each topology it probes all eight source octants near
the seven central face/edge/corner combinations: 14,336 production rows in
total. Every row must have current leaf support, nonnegative weights, unit sum,
and affine coordinate reproduction. Duplicate and permutation cases validate
logical-row algebra, while cache replay must invalidate after adjacent
refinement. The input-immutability case guards the production contract that an
interpolation coordinate is read-only.

The exhaustive masks exposed three production selection defects. A coarse cell
surrounded by fine support must use the `<=` level selector; every applicable
face/edge/corner neighbor slot must participate before choosing the coarsest
candidate; and the edge-1/edge-3 geometric predicates must reference their
actual transverse coordinates. The tests retain their affine oracle and tight
roundoff tolerance; these defects were fixed in production rather than hidden
by weaker acceptance.

The cut-cell cases use independent Euclidean and integral oracles.  A 3-4-5
right triangle fixes area, winding, and unit normal; direct ray/plane equations
fix hit parameters and spatial `EPS` behavior.  Parallel, near-tangent,
zero-length, and zero-area inputs are required to return a typed miss without a
NaN. A closed outward-wound box first places its near face at the dyadic
`x=0.25` m plane, giving exactly 5.0 m^3 in the `[-1,1]^3` m cell, and then at
the non-dyadic `x=0.30` m plane, giving 5.2 m^3 within the requested relative
error. Production intersection tests use the open block interior so mere
contact with a block face is not counted as a cut. Recursive depth follows the
interface-dominated leading error `O(A*h)`. Translated and
scaled spherical exclusions are checked against `4*pi*R^3/3`, including a
complete enclosing block and a single octant.  Planar sphere facets at 8x16
and 16x32 angular resolution must show the expected second-order area-error
reduction toward `4*pi*R^2`.

The NASTRAN fixture and its format contract are documented in
`fixtures/README.md`.  The reader treats sparse, out-of-order `GRID*` IDs as
keys.  The writer emits deterministic consecutive IDs and promises a
read/write/read round trip of coordinates, connectivity, winding, and face
attributes—not comments, sparse input IDs, record order, element IDs, or
per-file provenance. This scope is 3-D and one-rank; it does not qualify open,
nonmanifold, self-intersecting, large production, or distributed surfaces.
Detailed fixture syntax is in `fixtures/README.md`.

## Build and run

Native builds must follow the repository-wide preflight in `AGENTS.md`.  From
the confirmed repository root, after confirming no build/test process is
active, remove only `build`, configure exactly one profile, and build through
the top-level workflow with 16 jobs.  For example:

```
rm -rf -- build
./Config.pl -application=test/core-mesh-3d
./Config.pl -spice-path=nospice
env MAKEFLAGS="-j16" make amps
```

The unified runner is the normal entry point. It performs the same clean
preflight separately for all three dimensional applications, verifies the
effective generated `DIM`, and launches each case in a fresh one-rank MPI
process:

```
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-MESH-NATIVE-001 --order normal \
  --output-dir /tmp/amps-t3-qualified-normal
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-MESH-NATIVE-001 --order reverse \
  --output-dir /tmp/amps-t3-qualified-reverse
python3 src/test/run_tests.py --config test/qualification-machine/passion.json \
  --test CORE-MESH-NATIVE-001 --order shuffle \
  --shuffle-seed 1729 --output-dir /tmp/amps-t3-qualified-shuffle-1729
```

For direct executable diagnosis, use the configured AMPS
binary's standard location and the structured CLI:

```
mpiexec -n 1 ./amps --list
mpiexec -n 1 ./amps --case CORE-MESH-3D-INTERPOLATION-AMR-001
```

`--filter`, repeated `--case`, `--order normal|reverse|shuffle`, and `--seed`
have the same semantics as the core native harness.  Mesh-mutating acceptance
cases should be launched one case per process by the outer runner; `--list`
does not construct a mesh.

## Implemented versus validated

The source fixtures and independent analytic oracles are implemented and
qualified for the bounded T3 scope. The first normal/reverse/seed-1729
campaigns each passed 58/58, but an independent audit subsequently corrected
one additional lower-dimensional array-bound defect, so those receipts remain
bring-up evidence. The three corrected-tree commands above each pass 58/58.
Their per-profile executable hashes agree across all orders: 1-D
`54e79ddf1c34e2910a059a23e38232bea6322ab7a1ee6d31e35204acd8dbf4e0`, 2-D
`8d5311b99bfba3812599336b2cdbb4fbcd53ab16a3e625122b6db44a30b4aad9`, and
3-D `96f13a4be39ed221b7a252bb8313db3371ed49b2fd464ea1b5d09a545bf42e34`.
Source inspection or a syntax-only compile is not a native PASS. The tests do not qualify MPI
decomposition,
axisymmetric/spherical cell measures, arbitrary application-associated data,
or an unconfigured interpolation backend. The legacy `OutputLayerNodeMap`
path remains hard-coded for 3-D and is not executed or qualified by the 1-D/2-D
profiles. The cut-cell volume/parity oracle
requires a closed, consistently wound, manifold triangulation; this stage does
not qualify self-intersections, open/nonmanifold surfaces, or multi-rank
surface ownership.

The three profile decks match the repository's existing `*.input` ignore rule.
Codex did not modify `.gitignore`; a human preparing a commit must explicitly
force-add `input/test/core-mesh-{1d,2d,3d}.input`.

The 1-D/2-D decks enable `InternalBoundary` only to generate the shared type
graph required by AMPS headers; the test applications register no physical
internal boundary or surface source.
