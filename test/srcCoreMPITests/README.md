# AMPS T10 collective MPI fixture

This configured application validates AMPS' native distributed AMR/PIC paths.
It is not a physical model and registers no source, plasma, surface, or
external coupling. Every case is collective: all MPI ranks construct the same
replicated AMR tree, enter the same communication operations, publish a local
assertion record, and participate in the final failure reduction.

## Compile profiles

- `input/test/core-mpi.input` creates a 3-D nonperiodic tree with one refined
  root child. The 15 used leaves provide both same-level and 2:1 coarse/fine
  contacts while remaining small enough for repeated 1/2/4-rank runs.
- `input/test/core-mpi-periodic.input` enables AMPS' compile-time periodic
  boundary mode and uses the production periodic mesh/pair manager.

Both profiles use two-cell blocks, one ghost cell, H2O as a neutral record,
species-global time step/weight storage, no SPICE, and no model modules. The
fixture reserves two doubles on every center and corner node before associated
data storage is frozen. Those fields make owner-to-ghost and block-migration
checks independent of the built-in sampling payload.

## Independent oracles and ownership conventions

The AMR tree is replicated; each leaf's `Thread` is the global ownership
authority. After domain-boundary ghosts are installed, only
`ParallelNodesDistributionList[ThisThread]` is a traversable owner load list.
Remote leaves reuse the intrusive next/previous links in this rank's received-
ghost lists. Tests therefore derive the global owner map from the replicated
tree and separately require the local load list to contain every and only
local leaves.

Geometry identities and expected values use dyadic bounds/refinement indices,
not production communication output. Center coordinates use
`x_min + (i+1/2) Delta x`; triangle ray hits use barycentric inequalities;
particle totals/momenta/energy use direct SI definitions; and rank-invariant
checksums are commutative XOR/FNV receipts independent of list order. Floating
reductions use a relative `2e-14` gate, while integer counts, ownership,
particle identity, associated-data values, and serialized checksums are exact.

## Cases

| Family | Cases | Production path |
|---|---:|---|
| Ownership/geometry | 4 | replicated owner map, local load list, received ghosts, weighted repartition, distributed cut-cell ray geometry |
| Halo/serialization | 4 | center-node Replace/Add, send halo and same/coarse-fine masks, block associated data plus particle pack/unpack |
| Migration/boundaries | 4 | ghost-to-owner exchange, block/particle repartition, outer deletion, periodic ghost-to-real transfer |
| Portability/restart/invariants | 4 | CPU-list/online discovery, same and changed decomposition restart, global number/weight/mass/momentum/energy/rate checksum |

The production block send-mask callback is invoked only for an actual
`(owned block,destination)` relation, matching its AMR caller contract. Its
historical no-neighbor fallback creates a basic mask and is not used as a
neighbor-discovery predicate.

## Run

From `/home/vtenishe/Mars2/AMPS`:

```sh
python3 src/test/run_tests.py --test CORE-MESH-MPI-OWNERSHIP-001 \
  --output-dir /tmp/amps-t10
```

Aliases select bounded families without adding another build owner:

```sh
python3 src/test/run_tests.py --test CORE-PIC-MPI-HALO-001
python3 src/test/run_tests.py --test CORE-PIC-MPI-MIGRATION-001
python3 src/test/run_tests.py --test CORE-PIC-MPI-RANK-INVARIANCE-001
```

The adapter enforces the required clean native preflight before each selected
compile profile: confirm repository/process state, remove only `build`, select
the application through `Config.pl`, disable SPICE through `Config.pl`, and
compile through the top-level workflow with `MAKEFLAGS=-j16`. It retains each
executable and generated configuration identity before running cases.

`native-t10-mpi-output.txt` contains every build/run transcript;
`t10-matrix-receipts.json` records the rank-count receipts; the global JSON and
JUnit reports contain one aggregated result per stable case. A case passes
only if exactly ranks `0..N-1` publish local PASS at each of N=1,2,4 and rank
zero publishes one matching collective PASS summary.

## Qualified boundary and remaining work

T10 qualifies CPU domain-boundary-layer communication for this bounded 3-D
Cartesian fixture and fixed-weight neutral record. It does not claim CUDA,
hybrid/OpenMP list concurrency, scheduler-authorized affinity mutation,
distributed source injection, ECSIM's consumption of remote halo values in a
multi-block solve, mesh-topology file restart, non-Cartesian meshes, or a
production application's external callback state.
