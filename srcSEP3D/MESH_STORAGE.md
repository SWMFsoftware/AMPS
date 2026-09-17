# Phase M Mesh and Storage

Phase M defines one mesh-resolution and storage contract for both the
AMPS-independent verifier and the production AMPS boundary. The implementation
is in `mesh/mesh_model.{h,cpp}`; the production adapter is in `main_lib.cpp`.

## Domain contract

The numerical mesh is a Cartesian cube. `DomainBounds` retains both the cube
and the physical heliocentric shell:

- Earth preset outer radius: 1 AU;
- Mars preset outer radius: 1.666 AU;
- inner radius: supplied by immutable run configuration;
- Cartesian bounds: `[-r_outer,+r_outer]` on every axis.

The inner sphere is a physical boundary, not a reason to change the Cartesian
allocation. Production background filling excludes cell centers inside the
inner radius; their allocated storage remains zero-initialized until a later
boundary phase defines its behavior.

## Resolution law

The radial request is

\[
h_r(r)=\operatorname{clip}\left(h_{\min}\frac{\max(r,r_{\rm in})}
{r_{\rm in}},h_{\min},h_{\rm bg}\right).
\]

It equals the minimum size on the inner sphere, doubles when radius doubles,
and joins the background size continuously.

The optional Parker-tube centerline has a configured longitude and colatitude
at the inner radius. Its longitude winds by

\[
\Delta\phi=-p\,\Omega_\odot(r-r_{\rm in})/V_{\rm sw},
\]

where `p` is ±1. `TubeDistanceM` uses
`r*atan2(|r_hat × t_hat|, r_hat·t_hat)`, which remains first-order accurate for
nearly coincident directions. Inside the core the tube requests `h_tube`;
between core and shoulder radii it linearly joins `h_bg`. The final request is
the minimum of radial, tube, and background requests, clipped to declared
bounds.

## Standalone octree and ownership

`StandaloneOctree::Build` recursively refines blocks whose cell width exceeds
the smallest request sampled at their center/corners. It then enforces 2:1
face-neighbor balance. Leaves are sorted deterministically by spatial minimum
and level before receiving:

- `globalLeaf`;
- `firstCell`;
- round-robin `ownerRank`.

`CellStorage` proves that only the owner may write a global cell. This model is
not a second production mesh; it is an AMPS-free oracle for the resolution,
balance, identity, storage, and memory contracts.

The memory estimate is

\[
N_{cell}(B_{static}+B_{sample})+N_{leaf}B_{block-overhead}.
\]

A build that exceeds the configured budget fails before allocating the
standalone storage image.

## Frozen cell layout

Offsets are bytes relative to the srcSEP3D static-data base assigned by AMPS.
The default layout contains 17 doubles (136 bytes):

| Order | Field | Count | Units |
|---:|---|---:|---|
| 1 | magnetic field `B` | 3 | T |
| 2 | bulk velocity `U` | 3 | m/s |
| 3 | number density | 1 | m⁻³ |
| 4 | velocity divergence | 1 | s⁻¹ |
| 5 | temperature | 1 | K |
| 6 | pressure | 1 | Pa |
| 7 | Alfvén speed | 1 | m/s |
| 8 | `div(b_hat)` | 1 | m⁻¹ |
| 9 | focusing length | 1 | m |
| 10 | field-line curvature | 3 | m⁻¹ |
| 11 | field-aligned strain | 1 | s⁻¹ |

Optional fields follow in this exact order: magnetic gradient (9 doubles),
velocity gradient (9 doubles), and—for SWMF turbulence authority—directional
magnetic wave variances along/against `+B` (2 doubles, T²). Sampling bytes are
tracked separately because AMPS duplicates/switches sampling buffers according
to its sampling configuration.

## Production initialization order

`amps_init_mesh()` performs this order:

1. require an immutable configured `Runtime`;
2. validate the Phase-M resolution configuration;
3. initialize AMPS allocation registries;
4. register srcSEP3D static/sampling byte callbacks;
5. let AMPS freeze its complete center-node layout;
6. build the tree with the tested `localResolution()` function;
7. partition the tree and create owner lists;
8. allocate blocks and initialize cell measures;
9. bind the exact frozen `StorageLayout` to `Runtime`.

Changing or appending fields after step 5 is a layout error. Background filling
iterates `DomainBlockDecomposition::BlockTable`, so each rank writes only cells
in blocks assigned to that rank.

## Gradient reconstruction

`ReconstructScalarGradient` and `ReconstructVectorGradient` solve weighted
least-squares normal equations using actual neighbor displacement vectors.
Mixed neighbor distances are therefore valid at coarse/fine AMR interfaces.
A stencil without three-dimensional rank returns an error; it never invents a
zero gradient.

## Evidence

| IDs | Contract |
|---|---|
| `MSH3D01` | one million points remain between resolution floor/background bounds |
| `MSH3D02` | radial surface, transition, and octave closed forms |
| `MSH3D03–04` | Parker-tube centerline and convergence |
| `MSH3D05–06` | shoulder balance, live negative control, rotation invariance |
| `MSH3D07` | five octrees, exact histograms/memory, owner-only storage |
| `MSH3D08` | Earth/Mars preset extents |
| `MSH3D09` | coarse/fine linear exactness and rank-deficient rejection |

Run `test/run_tests.py --suite phase-m --rebuild`.
