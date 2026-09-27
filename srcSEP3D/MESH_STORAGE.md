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
at the inner radius. It is the integral curve of the same SWCME Parker field
that initializes the background. For source radius \(r_0\), its longitude is

\[
\Delta\phi=-\frac{\Omega_\odot}{V_{\rm sw}}
\left[(r-r_0)-r_0\ln\left(\frac{r}{r_0}\right)\right].
\]

This follows from
\(B_\phi/B_r=-\Omega_\odot(r-r_0)\sin\theta/V_{\rm sw}\) and
\(r\sin\theta\,d\phi/dr=B_\phi/B_r\).  The older
`-Omega*(r-r0)/V` angle was not an integral curve of the initialized field and
could displace the 1-AU mask by more than its requested width.

Magnetic polarity is absent from geometry: reversing polarity changes the
analytic field and pitch orientation, not the refined tube. `TubeDistanceM` uses
`r*atan2(|r_hat × t_hat|, r_hat·t_hat)`, which remains first-order accurate for
nearly coincident directions. The tube radius is declared at a reference
heliocentric distance and is either physically constant or proportional to
radius (constant angular width). A named transverse profile joins `h_tube` on
the centerline to `h_global` exactly at the tube boundary. The final request is
the finer of radial and tube requests, clipped to declared bounds.

## Active Parker transport corridor

AMPS represents activity with
`cTreeNodeAMR::IsUsedInCalculationFlag` on a complete leaf block. Its public
`SetTreeNodeActiveUseFlag` gathers selected AMR node IDs from every MPI rank,
broadcasts the common flag change, and suppresses block allocation and load for
inactive leaves. Movers that enter an inactive leaf follow AMPS' ordinary
left-domain path. There is no public per-finite-volume-cell allocation switch;
therefore `[mesh.active_region]` intentionally describes a block mask even
though every cell in the block is disabled together.

`mode = full-domain` preserves historical behavior. For
`mode = parker-tube`, the physical active radius is

\[
R_a(r)=R_{a,ref}
\quad\hbox{or}\quad
R_a(r)=R_{a,ref}\,r/r_{a,ref-radius},
\]

for `physical-constant` and `constant-angular-width`, respectively. The active
radius and reference are independent of the refinement radius, but validation
requires the active tube to contain the refined tube at both endpoints and
therefore everywhere along the finite active line for the supported constant
and linearly radius-scaled width laws.

The active-mask planner uses the *finite* `[parker_spiral]` line, clipped at the
physical outer sphere.  It subdivides the exact analytic curve independently
of diagnostic `point_count`, expands every chord by the physical tube radius
plus a rigorous curve-to-chord bound, and evaluates the exact Euclidean
segment-to-axis-aligned-box distance.  Consequently a curve that crosses a
leaf is retained even when it misses the leaf centre and all eight corners.
Malformed boxes remain active so AMPS' structural validator, rather than the
pruning pass, supplies the authoritative diagnostic.

`buffer_blocks=N` now means exactly `N` complete AMR touching-neighbour layers
(faces, edges, and corners), using AMPS' native coarse/fine neighbour links.
It is never converted to a candidate leaf's physical diagonal.  This is
important at resolution transitions: diagonal scaling gave coarse leaves an
oversized halo while rejecting intervening fine leaves, producing the holes
and detached rectangular islands visible in active-mesh plots.

After halo dilation, the planner flood-fills inactive leaves from the
Cartesian boundary and promotes any bounded inactive cavity to a safety halo.
It then requires one face-connected active component from the source leaf to
the finite-line endpoint.  A disconnected core, detached halo, or missing
endpoint stops initialization before any AMPS block storage is allocated.

The replicated leaf list is partitioned deterministically by ordinal modulo
MPI rank before calling `SetTreeNodeActiveUseFlag`; every node ID is submitted
exactly once. The call occurs after `buildMesh()` and before parallel load
measurement, distribution-list creation, or block allocation. Configuration
also rejects fixed observers whose collection spheres do not intersect the
physical tube and rejects moving/field-connected observers with a static mask.
Shock source patches outside the active leaves receive explicit disconnected
ledger rows; the exact per-species macro count is apportioned only over
connected physical patches.

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
The default layout contains 19 doubles (152 bytes): 17 background primitives
plus two directional turbulence variances.

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
| 12 | `deltaB_plus_squared`, `deltaB_minus_squared` | 2 | T² |

The actual canonical allocation order is the 17 background primitives,
optional magnetic gradient (9 doubles), optional velocity gradient (9
doubles), and then the two wave variances. The variances are allocated for both
prescribed and SWMF turbulence, so the initialization file always exposes the
state used by scattering. SI energy densities are derived for output as
`w_plus/minus=deltaB_plus/minus_squared/mu0`; their sum is emitted as
`turbulence_wave_energy_density_J_per_m3`. The public six-column group (total,
plus, and minus variance followed by total, plus, and minus energy density) is
mandatory and is tested independently of AMPS. These derived quantities are
not duplicate mutable cell fields. Sampling bytes are tracked separately because AMPS
duplicates/switches sampling buffers according to its sampling configuration.

## Production initialization order

`amps_init_mesh()` performs this order:

1. require an immutable configured `Runtime`;
2. validate the Phase-M resolution configuration;
3. initialize AMPS allocation registries;
4. register srcSEP3D static/sampling byte callbacks;
5. let AMPS freeze its complete center-node layout;
6. build the tree with the tested `localResolution()` function;
7. apply the active-use mask through AMPS' public node-ID synchronization API;
8. install an active-only load measure, partition the retained tree, and
   create owner lists;
9. for schema 3 or 4, write the final distributed tree with AMPS'
   `outputMeshTECPLOT` and write the finite Parker centreline once on rank zero;
10. allocate active blocks and initialize cell measures;
11. bind the exact frozen `StorageLayout` to `Runtime`.

Changing or appending fields after step 5 is a layout error. Background filling
iterates `DomainBlockDecomposition::BlockTable`, so each rank writes only cells
in blocks assigned to that rank.

The two schema-3/schema-4 paths are explicit input fields. The centreline writer builds
and validates the entire ordered curve before opening its output and records
arc length, Cartesian position, heliocentric radius, and requested cell size in
metres. A write failure is fatal; these are initialization products, not
best-effort diagnostics.

`amps_init()` then performs the data-bearing sequence:

1. install positive finite global and block-local time step/weight for every
   species compiled from `SpeciesList`;
2. build and validate the complete owner-local background snapshot;
3. zero the application and native DATAFILE records so Cartesian padding is a
   finite, deterministic interpolation source;
4. after preparing the selected turbulence provider, prescribe background and
   directional turbulence together to every owner-local physical center node,
   read back the two variance values, and collectively verify prescribed
   turbulence is positive in every physical cell;
5. exchange associated-data halos (and, for relativistic GCA, generate and
   exchange its neighbor-dependent derived fields);
6. publish/install runtime and mover state, including an optional restart; and
7. verify the active snapshot and every species numerical value, then call
   `outputMeshDataTECPLOT` for `sep3d-initialization-data.dat` as the final
   operation of `amps_init()`.

The separation between steps 9–10 of `amps_init_mesh()` and this sequence is
intentional: `outputMeshTECPLOT` needs only the finalized octree, whereas
`outputMeshDataTECPLOT` must not observe the native buffer before Parker/SWCME
and turbulence installation have completed on every rank.

`outputMeshDataTECPLOT` represents a FEBRICK at its vertices. AMPS creates a
temporary center node for each vertex and calls the registered center-node
interpolators before printing it. Static bytes requested by an application are
not included by AMPS' built-in interpolation. srcSEP3D consequently registers
`InterpolateInitializationCellData` alongside its print callbacks before layout
freeze. It interpolates all `cellAssociatedBytes / sizeof(double)` values as a
single state vector, so the two turbulence variances cannot remain zero while
the native DATAFILE plasma/IMF fields are populated. Cartesian padding is
explicitly zeroed before the physical-cell fill and is identified in output by
`background_valid=0`.

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
| `MSH3D10–11` | exact equal-arc finite line/origin invariance and Tecplot initialization output |
| `MSH3D12` | conservative finite capsule/box intersection and full-domain identity |
| `MSH3D13` | analytic curve derivative agrees with the Parker tangent and arc-length inversion round-trips |
| `MSH3D14` | whole-octree core coverage, exact topological halo depth, pruning, cavity elimination, and connectivity |
| `CFG3D03–05` | normalized domains, shared Parker geometry, composite preflight and whole-run memory |

Run `test/run_tests.py --suite phase-m --rebuild`.
