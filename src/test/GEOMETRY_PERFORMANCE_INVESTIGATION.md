# Hartley2 initial cut-cell volume cost

Date: 2026-11-10. Scope: investigate the slow native 3-D initialization exposed
after the approved Hartley2 build fixes. This report proposes an optimization;
the production geometry algorithm has not been changed by this investigation.

## Observed results

The two application fixes compile and link successfully through the clean
native workflow. Both builds use the unchanged canonical `input/hartley2.input`,
DIM3, species-global weights, individual corrections ON, radial backend OFF,
and brightness OFF.

| Attempt | Result |
| --- | --- |
| Site-default native build, `OPT3=-g` | PASS, 63.932 s |
| Site-default four-rank runtime | Intentionally interrupted after 587.993 s during volume construction; zero particle steps |
| Native build with make-only `OPT3=-O2 -g` | PASS, 97.007 s |
| Optimized four-rank runtime | Local 1,800 s timeout; measured 1,800.231 s; zero particle steps; no completed statistics file |

The local runtime receipt labels the timed-out attempt `FAIL` and records
`failure="native runtime exceeded 1800 seconds"`, with no normal exit code.
This is an incomplete runtime validation, not a measured physical rejection.
The 1,800 s ceiling is an execution safeguard introduced by the local smoke
wrapper; it is not a physics, accuracy, or convergence threshold. All owned MPI
processes had exited when checked after the timeout.

The observed mesh has 3,858 leaf blocks and a 32,040-triangle surface. All four
ranks reported changing processed-block counts during initialization. A brief
GDB sample found this active path:

```text
InitCellMeasureBlock
  -> CutCell::GetRemainedBlockVolume (recursive subdivision, level <= 5)
     -> PointInTriangulatedComputationalDomain
        -> cTriangleFace::RayIntersection
```

There is no completed 100-iteration 3-D run, sampled flyby result, or scientific
acceptance associated with either attempt.

## Provenance and reproducibility

All evidence below is local generated output under
`test_output/comet-integration-20261109/`, excluded from the source import:

- `native/hartley3d-final-20261110/build.json` and `build.log`;
- `hartley3d-live/result.json`, `interruption.json`, and
  `initial-volume-stack.log`;
- `native/hartley3d-optimized-20261110/build.json`, `build.log`, and
  `selection-audit.json`;
- `hartley3d-live-optimized-20261110/result.json`, `command.json`, and `run.log`;
- `provenance/hartley3d-volume-path-audit.json`;
- `provenance/final-validation-commands.json`.

The build drivers record the exact root/idle check, removal of only `build/`,
Config.pl regeneration, and top-level `MAKEFLAGS=-j16` compilation. The optimized
build passes `OPT3=-O2 -g` on the make command line. The site `Makefile.conf`,
canonical deck, generated physics flags, and generated application sources
have identical checksums across the two builds.

The optimized executable SHA-256 is
`a8762dbfbf724417da30b496534126fc4bbe80e333fe202f7ffbdf11d9201b59`.

## Identified bottleneck

[PointInTriangulatedComputationalDomain](../meshAMR/meshAMRcutcell.cpp)
casts one fixed-direction ray and invokes the existing triangle intersection
predicate on every face. It then sorts intersection parameters, merges shared
crossings using the existing EPS/roundoff tolerance, and evaluates parity.
One point query therefore performs 32,040 triangle tests in this case.

The volume integrator recursively subdivides cells intersecting the surface.
It classifies the center of every nonintersecting child and every terminal cut
child against the complete surface. The native call requests relative error
0.01; its calculated depth is capped by the compiled maximum of 5. A fully
refined octree at this depth has

```text
8^5 = 32,768 terminal cells
32,768 * 32,040 = 1,049,886,720 triangle tests
```

The latter is a worst-case bound for one original cell, not an observed query
count. Early exits reduce it, but the linear full-surface search still repeats
at many distinct adaptive child centers. An allocation of capacity `nSurface`
for the intersection-parameter vector also occurs per point; its contribution
has not been isolated from the ray-test cost.

The per-block volume assembly routine is byte-identical to pre-integration
master `f6c6340b59d71173cfa825b2b242ee6440d341e9`. Its SHA-256 is
`1be9744322220c549e5488cf07347efdddc3e86fc4ad69407340b067d731c6f5`.
The only integrated cut-cell translation-unit code change is the escape-length
loop bound from `3` to `_MESH_DIMENSION_`, which remains `3` in this build.
Thus the sampled slow kernel is retained master geometry code. No claim is
made about the runtime of an unexecuted pre-integration Hartley2 application.

## Recommended implementation

Use a conservative spatial hierarchy, such as a triangle bounding-volume tree,
to enumerate candidates for the existing fixed ray. Keep the current
`RayIntersection` predicate, hit parameters, sorting, crossing de-duplication,
parity convention, subdivision, terminal reconstruction, and accuracy limits.
The intended change is candidate search work, with unchanged geometric
decisions.

Important implementation boundaries:

1. Build and reuse the hierarchy for one explicit immutable surface epoch.
   Pass a query context through recursive integration. Pointer/count equality
   alone cannot establish freshness when face coordinates are edited in place.
2. Bounds must conservatively include the predicate's spatial EPS,
   barycentric tolerance, origin-hit clamping, and roundoff behavior. Ambiguous,
   nonfinite, or unsuitable bounds fall back to the existing full scan.
3. Preserve compatibility for explicit caller-supplied surfaces as well as the
   process-global surface; retain existing public overloads.
4. Define context lifetime for sequential calls, mesh adaptation, geometry
   changes, MPI processes, and OpenMP tasks. Query scratch belongs to a call or
   worker; read-only acceleration data can be shared only while geometry is
   immutable. Account for the current projected-coordinate helper's write to
   the stored Gram determinant when reasoning about thread safety.
5. Retain the physical mesh, surface, EPS, refinement ceiling and requested
   error. A coarser problem would not demonstrate an optimization of this run.

Caching point answers alone is not the leading recommendation: the current
global classifications are predominantly at distinct recursive child centers.
Reusing the surface search structure addresses the full-surface scan instead.

## Required validation before accepting an optimization

- Retain the linear implementation as a differential oracle during development.
  Compare parity, hit parameters and volumes for deterministic probes on both
  synthetic surfaces and the actual 32,040-face surface. Use the same compiler
  flags when checking exact reproduction.
- Preserve independent analytic tests for closed planar cuts, translated/scaled
  spheres and volume complements. Differential agreement alone cannot prove
  physical correctness.
- Exercise origins on faces/edges/vertices, near-parallel rays, duplicate/shared
  crossings, degenerate/sliver triangles, translations/scales, zero/positive
  EPS, and explicit surfaces unrelated to the global mesh.
- Test surface replacement and in-place geometry changes between query epochs,
  plus supported threaded use and MPI owner/ghost volume agreement.
- Measure triangle-test counts, hierarchy build/rebuild work, memory, and wall
  time. Freeze a deterministic operation-count regression fixture so machine
  speed is not the only performance criterion. No speedup is claimed yet.
- Register new tests in the mesh native owner, manifest, traceability records,
  README and `TEST_CATALOG.md`, following the existing agent/catalog rule.
- Run DIM1/DIM2/DIM3 mesh tests, distributed mesh tests, the full core registry
  with exact baseline identity comparison, and the required Moon/fast-wave and
  Hartley2 application checks after any production geometry change.
- Rebuild canonical Hartley2 cleanly and complete all 100 particle iterations
  on four ranks. Retain the exact input, executable, cache provenance, logs and
  microscopic statistics. Keep this compatibility result separate from
  Hartley Stage-9 long-run science and CME qualification gates.

Investigation outcome: the build fixes are validated; the native 3-D runtime
remains open. Spatial-search acceleration is the proposed next implementation,
and has not been added to production in this investigation.
