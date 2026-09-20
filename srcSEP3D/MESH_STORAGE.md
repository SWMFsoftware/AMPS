# Phase M Mesh and Storage

Phase M defines one mesh-resolution and storage contract for both the
AMPS-independent verifier and the production AMPS boundary. The implementation
is in `mesh/mesh_model.{h,cpp}`; the production adapter is in `main_lib.cpp`.

## Domain contract

The numerical mesh is a Cartesian cube. `DomainBounds` retains both the cube
and the physical heliocentric shell:

- solar preset outer radius: 0.30 AU;
- one-AU (`earth` input alias) preset outer radius: 1 AU;
- Mars preset outer radius: 1.666 AU;
- inner radius: supplied by immutable run configuration;
- Cartesian bounds: `origin+[-r_outer,+r_outer]` on every axis.

There is no numeric automatic-radius sentinel. `OuterRadiusMode::Preset`
resolves one of the named values before fingerprinting, while
`OuterRadiusMode::Explicit` uses the supplied positive `outerRadiusM`.
Observer, shock, and reference locations are checked against those resolved
bounds before AMPS initialization.

The inner sphere is a physical boundary, not a reason to change the Cartesian
allocation. Production background filling excludes cell centers inside the
inner radius; their allocated storage remains zero-initialized until a later
boundary phase defines its behavior.

## Resolution law

The near-Sun request is a named interpolation. With
\(s=\operatorname{clip}[(r-r_{in})/(r_{transition}-r_{in}),0,1]\),

\[
h_r(r)=h_{surface}+P(s)(h_{global}-h_{surface}),
\]

where `P` is `linear`, `power-law`, or `smoothstep`, optionally raised to the
declared positive exponent. It equals the surface target at the inner sphere
and joins the global target continuously at the declared transition radius.

The optional Parker-tube centerline has a configured longitude and colatitude
at the inner radius. Its longitude winds by

\[
\Delta\phi=-\Omega_\odot(r-r_{\rm in})/V_{\rm sw}.
\]

Magnetic polarity is absent from geometry: reversing polarity changes the
analytic field and pitch orientation, not the refined tube. `TubeDistanceM` uses
`r*atan2(|r_hat × t_hat|, r_hat·t_hat)`, which remains first-order accurate for
nearly coincident directions. The tube radius is declared at a reference
heliocentric distance and is either physically constant or proportional to
radius (constant angular width). A named transverse profile joins `h_tube` on
the centerline to `h_global` exactly at the tube boundary. The final request is
the finer of radial and tube requests, clipped to declared bounds.

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

The former application-storage-only estimate has been replaced by a whole-run
planning model:

\[
M_{resident}=N_{cell}(B_{base-cell}+B_{static}+B_{sample}+B_{node}
+N_{particle/cell}B_{particle})
+N_{leaf}(B_{block}+B_{app-overhead}),
\]

plus per-block communication buffers, a configurable halo fraction, and a
configurable safety margin. Every ABI-dependent coefficient is explicit in
the immutable configuration. `BuildRefinementPreflight` samples limiting
surfaces without building the full tree, reports requested resolution extrema
and estimated blocks by level, and rejects an estimate above the configured
budget. Native integration evidence remains responsible for comparing this
conservative plan to actual AMPS peak memory.

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
8. for schema 3, write the final distributed tree with AMPS'
   `outputMeshTECPLOT` and write the finite Parker centreline once on rank zero;
9. allocate blocks and initialize cell measures;
10. bind the exact frozen `StorageLayout` to `Runtime`.

Changing or appending fields after step 5 is a layout error. Background filling
iterates `DomainBlockDecomposition::BlockTable`, so each rank writes only cells
in blocks assigned to that rank.

The two schema-3 paths are explicit input fields. The centreline writer builds
and validates the entire ordered curve before opening its output and records
arc length, Cartesian position, heliocentric radius, and requested cell size in
metres. A write failure is fatal; these are initialization products, not
best-effort diagnostics.

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
| `MSH3D02` | linear radial surface, midpoint, and transition closed forms |
| `MSH3D03–04` | polarity-independent Parker-tube centerline and convergence |
| `MSH3D05–06` | composite tube balance, live negative control, rotation invariance |
| `MSH3D07` | five octrees, exact histograms/memory, owner-only storage |
| `MSH3D08` | Earth/Mars preset extents |
| `MSH3D09` | coarse/fine linear exactness and rank-deficient rejection |
| `MSH3D10–11` | finite-line arc length/origin invariance and Tecplot initialization output |
| `CFG3D03–05` | normalized domains, shared Parker geometry, composite preflight and whole-run memory |

Run `test/run_tests.py --suite phase-m --rebuild`.
