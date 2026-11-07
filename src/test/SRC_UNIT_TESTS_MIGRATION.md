# `test/srcUnitTests` analysis and migration ledger

This ledger maps the legacy GoogleTest application to the portable core-test
registry. A row marked **ported** has an independent oracle in
`test/srcCoreUnitTests`; it does not mean that the entire owning production
component is qualified. Rows marked **open** remain visible work and must not
be inferred from a related passing case.

## Source-by-source mapping

| Legacy source and behavior | Registry replacement | State and analysis |
|---|---|---|
| `pbuffer.cpp`: initialization, allocate/delete counters | `PIC-PBUFFER-INIT-001` | **Ported.** Uses an isolated 64-record production buffer and exact allocation accounting. |
| `pbuffer.cpp`: position, velocity, species, clone | `PIC-PBUFFER-STATE-001` | **Ported.** Asymmetric values expose swapped dimensions/offsets. |
| `pbuffer.cpp`: intrusive-list insertion/removal | `PIC-PBUFFER-LIST-001` | **Ported.** Exact head/next invariants. |
| `pbuffer.cpp`: optional data requested before initialization | `PHYS-IDF-STORAGE-001` plus configured particle layout | **Ported for real LB optional fields.** A future generic layout-only profile should also request anonymous fields of multiple alignments. |
| `pbuffer.cpp`: individual statistical weight | none | **Open.** The active fixture deliberately compiles `ParticleWeightCorrectionMode=off`; this needs a separate mutually exclusive generated configuration. |
| `pbuffer.cpp`: full-buffer failure | none | **Open.** The legacy test itself is commented out. Add a subprocess test because the production fatal path must not terminate the suite runner. |
| `distribution.cpp`: Maxwellian means, speed, temperature | `PHYS-DISTRIBUTION-MAXWELLIAN-001` | **Ported.** Tests nonzero drift, all component means, analytic mean thermal speed, and recovered temperature with a fixed seed. |
| `idf.cpp`: storage access | `PHYS-IDF-STORAGE-001` | **Ported.** Checks configured offsets and every H2O vibrational mode. |
| `idf.cpp`: rotational equilibrium | `PHYS-IDF-ROTATION-001` | **Ported.** Uses the gamma/equipartition mean and a monatomic zero-mode boundary. |
| `idf.cpp`: vibrational equilibrium | `PHYS-IDF-VIBRATION-001` | **Ported.** Uses the independent harmonic-oscillator law per mode. |
| `idf.cpp`: microcanonical partition | `PHYS-IDF-PARTITION-001` | **Ported and corrected.** Resets accumulators for every degree pair and tests the beta mean plus every-draw conservation. |
| `idf.cpp`: RT transition | `PHYS-IDF-RT-001`, `PHYS-IDF-RT-RATE-001` | **Ported and corrected.** Tests the two-stage beta law, exact energy conservation, molecular/monatomic branch, and analytic `1/Z_rot=1/5` occurrence rate. The legacy molecular/monatomic assertion reads element 2 of a two-element array and is undefined behavior. |
| `collisions.cpp`: hard-sphere pair kernel | `PHYS-COLLISION-CONSERVATION-001`, `PHYS-COLLISION-ISOTROPY-001` | **Ported.** Unequal masses, both rotation branches, exact invariants, and solid-angle moments. |
| `collisions.cpp`: populated-cell heat bath | `PHYS-3D-COLLISION-RELAXATION-001`, `PHYS-3D-COLLISION-RATE-001` | **Ported for production NTC.** A real allocated AMR cell holds 160 H2O particles for RT equilibration and 240+240 H2O/O particles for the independent `n0*n1*sigma*<g>` rate. Energy/momentum and sampled rate are asserted. The other legacy collision selectors remain open. |
| `split_merge.cpp`: fixed/variable-weight split and merge | none | **Open.** Requires the individual-weight generated profile, then must test weight, species, momentum, energy, cell-list integrity, and deterministic limits. |
| `main_lib.cpp`: application mesh, particle movement, sphere interaction | `PHYS-MOVER-*`, `PHYS-3D-*` | **Ported as isolated assertions.** The fixture builds an eight-owner-block 3-D mesh, registers a sphere, allocates real particles, checks Simple/Boris free flight, analytic specular event time, NTC/LB relaxation, collision rate, and photolytic parent survival. |
| `sampling.h` helpers | no production case | **Not migrated.** The helpers only accumulated legacy test statistics; the new cases use explicit analytic estimators and confidence bounds. |

## Defects intentionally not copied

- `ParticleTestBase.h` can deallocate the Linux `aligned_alloc` buffer through
  a C++ array-delete path. The new scopes use `free` for that exact allocation
  mode and restore every process-global pointer/counter.
- `idf.cpp::RTTransitionTest` accesses `rotationalEnergyFraction[2]` although
  the array length is two. The replacement directly asserts the monatomic
  element at index zero and total conservation.
- Several legacy statistical accumulators live across nested parameter loops,
  so later expectations include earlier samples. Replacement cases reset all
  estimators locally and set their random seed inside the case.
- The legacy executable used to discard `RUN_ALL_TESTS()` and return success.
  It now propagates GoogleTest's status, while the new registry independently
  enforces case status, suite status, and process-exit agreement.

## Next individual-component suites

Priority additions that were absent or incomplete in `srcUnitTests` are:

1. Charged mover convergence in uniform `E`, uniform `B`, and crossed `E x B`,
   including gyro phase/radius, energy, observed order, block/refinement
   crossings, outer-boundary loss, absorbing surfaces, and moving surfaces.
2. The remaining collision-cell algorithms (`mf_simplified`,
   `ntc_simplified`, `mf_Yinsi`, and `mf`) against the same rate/conservation
   fixture; unequal time steps/weights; zero-relative-speed and majorant/
   collision-limiting branches.
3. A separate individual-weight build for split/merge and reaction-daughter
   weight bookkeeping, followed by injection normalization and sampling
   normalization tests.
4. Mesh topology, refine/coarsen conservation, cut-cell volume/face geometry,
   interpolation, restart round trip, and deterministic owner/ghost tests.
   The first variants can run at one rank; MPI ownership requires later
   one/two/four-rank evidence.
5. Species/mass/charge tables, units and round trips, RNG stream
   reproducibility/independence, special functions, quadrature, invalid input,
   fatal diagnostics in subprocesses, and serialization corruption tests.
6. Photochemical three-body angular laws and daughter insertion/weight
   conservation; charge exchange, electron impact, chemical solver stiffness/
   positivity/conservation, surface adsorption/desorption, sputtering, dust,
   ECSIM discrete operators, and gyro guiding-center validity boundaries.

These items are also tracked by `component_inventory.json`; this document is
the detailed legacy-to-registry provenance rather than a claim of completion.

One production discrepancy is intentionally explicit: the current LB
`RedistributeEnergy` path uses a hard-coded `Zrot=5`, although the input
generator emits `RotationZnumber` and LB initialization validates it. The
current occurrence test qualifies that implemented closure. It does not claim
that changing `RotZnum` changes the rate; defining the physically correct
cross-species combination requires a separate reviewed production change.
