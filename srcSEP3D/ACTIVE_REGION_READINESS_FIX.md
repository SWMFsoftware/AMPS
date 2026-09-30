# Active-region plan readiness and corridor compatibility — 2026-09-30

This cumulative source package fixes `SCCM3D04` falsely rejecting a verified
full-domain plan when no complete leaf needs deactivation. It retains the
earlier inactive-shock readiness fix described in
`INITIALIZATION_READINESS_FIX.md`.

## Production behavior

The active-region installed flag is set only after every leaf's AMPS active-use
flag matches the computed plan. Whether the plan actually removes blocks is
recorded separately. Allocation verification now runs for both full-domain and
corridor modes, and stops initialization for missing owner blocks, resident
inactive blocks or counts inconsistent with the plan.

The existing corridor geometry is unchanged:

- `mode=parker-tube` retains a conservative finite corridor around the selected
  analytic Parker field line, plus the configured topological halo and any
  bounded cavities filled by the safety policy.
- Whole leaves entirely inside the physical Sun are excluded even if selected
  by the corridor or halo. Surface-intersecting leaves retain AMPS cut-cell
  treatment; fully interior cells retain zero physical volume.
- Source-to-endpoint connectivity, observer coverage, refinement/active-width
  consistency and source-patch connectivity checks remain in force.
- A corridor that happens to retain all available leaves is still a valid plan
  and emits the existing no-pruning warning. No outside leaf is activated just
  to make a test pass, and no arbitrary leaf is removed to force pruning.

The provided corridor remains static and tied to the analytic Parker curve.
Arbitrary PFSS/imported or evolving field-line corridors require their own
provider/geometry and reactivation contracts; this patch does not implement
those future capabilities or fix the separate first-step crash.

## Native evidence

SCCM3D04 now reports the particular failure and requires a registered solar
boundary, a verified installed plan, verified allocation, a supported mode,
nonzero active-leaf count and allocated-block count equal to the active plan.
Pruning evidence must match the inactive count. Solar-interior leaves must be
included in the excluded count; in full-domain mode they are the only excluded
leaves. Wide corridors do not fail merely because pruning is zero.

The existing native JSON schema is retained with an additive `active_region`
object. The text state artifact also includes mode, installation, pruning,
allocation verification and the active/inactive/solar-interior counts. These
fields enter the native MPI identity fingerprint.

## Install and rebuild

Extract the complete archive into the AMPS root, preserving its directory
structure. It overlays the original submitted package and includes the prior
fix. The archive is not a complete AMPS checkout: keep your existing AMPS,
SWCME and sep_common build files.

From the AMPS root, rebuild using the same compiler/MPI configuration:

```bash
make -C src/models/sep_coronal_cme -j8 lib
make clean
make -j8
```

The changed `NativeApplicationState` header requires rebuilding all its callers.
The included `srcSEP3D/initialization_readiness_fix.patch` is a cumulative patch
against the original upload, including both fixes. Use it as an alternative to
installing the updated sources, not after extracting those sources.

## Verify the future corridor configuration

The existing complete `sep3d_analytic_parker_active_tube.in` example selects
`parker-tube`, an active width of 0.05 AU at 1 AU, and one block halo. Its
observers are aligned with the configured field-line corridor. Run all registered coupled
native initialization checks from the AMPS root, without particle stepping:

```bash
mpiexec -n 4 ./amps --test-suite sep-corona --test-input srcSEP3D/examples/sep3d_analytic_parker_active_tube.in --test-steps 0 --expect-mpi-ranks 4 --test-json test_output/corridor-init/native.json --artifact-directory test_output/corridor-init/artifacts
```

For the unpruned reference, use the same command with
`sep3d_analytic_parker.in` and a separate output directory.

The static mask operates at complete-leaf resolution: the computational
corridor includes whole retained blocks and its halo rather than a sharp
per-cell radial cutoff at the tube width. The physical solar surface remains
the internal absorption boundary; the configurable transport/source shell is
separate and unchanged.

## Portable regressions and validation limits

```bash
python3 srcSEP3D/test/run_native_boundary_regression.py
```

This compiles the actual native evaluator and shared geometry. Its 45 checks
cover full-domain zero pruning, solar-only pruning, normal and wide corridors,
missing installation/allocation, incorrect counts, inconsistent pruning,
unsupported modes, suite selection and CLI conflicts, solar-exclusion evidence,
real finite-line core/halo/pruning
geometry, and cut-cell/exterior preservation. The runner also verifies the new
JSON evidence round-trip. Those checks passed when preparing this package.

The existing 19 SEP3D Python tests also passed. The earlier shared-model suite
passed 196/196 after the retained inactive-shock fix; its sources did not change
in this update. Full AMPS compilation, native MPI initialization and particle
stepping remain to be verified on the user's configured AMPS installation.
