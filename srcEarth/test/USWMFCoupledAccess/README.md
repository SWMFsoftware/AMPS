# Step 10: SWMF-coupled cutoff and directional-access validation

This suite validates the roadmap Step 10 production contract: a complete SWMF receive
is evaluated at the configured physical cadence, one immutable Step-9 field snapshot is
used for the entire cutoff/access calculation, and the resulting product is directly
comparable with standalone `FIELD_MODEL SWMF_SNAPSHOT` replay of that exact state.

The forward Monte Carlo density and sphere-flux samplers are **not** a Step-10
reference. They solve a different problem. Step 10 compares the live coupled Mode3D
backtracer with the same Mode3D backtracer reading the exported frozen field.

## Fast numerical/source suite

From the AMPS repository root:

```bash
./srcEarth/test/USWMFCoupledAccess/run_test.sh
```

The runner is independent of MPI/SWMF and returns nonzero on any failure. Every
component and the runner itself write an explicit `RESULT: PASS` or `RESULT: FAIL`.

| ID | Test | Independent reference or failure condition |
|---|---|---|
| S10-U01 | Cadence state machine | Fixed 0/30/60-s schedule; duplicate skip, stale-time rejection, and fresh-restart evaluation |
| S10-U02 | Product identity | Fixed nine-decimal time/snapshot suffix; restart identity and content/time sensitivity |
| S10-U03 | Shue parameters | Published analytic Shue coefficients at Pdyn=2 nPa and IMF Bz=-5 nT, with fixed subsolar/terminator radii |
| S10-U04 | Escape classification | Manufactured BOX, Shue nose, Y-face data-limit, and X-tail-cap crossings with exact fractions |
| S10-U05 | Artifact manifest | Required snapshot/content/mesh/time fields, quasi-static label, nonempty unique artifact list, and result/suffix-identity rejection |
| S10-U06 | Replay comparator | Exact manufactured live/replay Tecplot product |
| S10-U07 | Comparator gate | A changed cutoff fails the default exact gate; tolerance works only when explicitly supplied |
| S10-U08 | Provenance rejection | Different/missing field identity, changed zone topology, or changed resolved Shue surface fails before numeric comparison |
| S10-U09 | Legacy-consumer compatibility | The real C9 and C10 access/penumbra readers consume production-ordered, validated Step-10 `AUXDATA` and still recover the manufactured reference rows |
| S10-U10 | Fail-closed consumer gate | Malformed `AUXDATA` and inconsistent `access_state`/`allowed`/`unresolved` flags remain hard failures in both C9 and C10 |
| S10-SOURCE | Production wiring | Collective cadence, post-success commit, final content-derived suffix, output inventory, active boundary policy, and unchanged F4/C9/C19 gates |

S10-U03 and S10-U04 are numerical reference-solution tests, not string/source checks.
They exercise the same dependency-free functions called by the production tracer.

S10-U09 closes a producer/consumer regression exposed by the C9/C10 campaign: Step 10
correctly added provenance records to the C++ artifacts, but the older observational
readers initially attempted to parse those records as floating-point rows.  The repair
keeps the provenance and recognizes only syntactically valid `AUXDATA KEY="value"`
records. Unknown text, malformed metadata, missing variables, wrong row widths, and
all existing observational thresholds continue to fail.

## Outer-boundary meaning

`BOUNDARY_TYPE BOX` preserves the historical behavior: every computational face is a
physical escape surface. `BOUNDARY_TYPE SHUE` uses

```text
r(theta) = r0 [2/(1+cos(theta))]^alpha
```

plus `DOMAIN_X_MIN` as the nightside tail cap. `SHUE_R0 AUTO` and
`SHUE_ALPHA AUTO` are resolved from configured `PDYN` [nPa] and `IMF_BZ` [nT]; numeric
`SHUE_R0` is in Earth radii and numeric `SHUE_ALPHA` is dimensionless.

For `SHUE`, reaching YMAX/YMIN/ZMAX/ZMIN/XMAX while still inside the physical
magnetopause is `INVALID_FIELD`, not `OUTER_BOUNDARY_ALLOWED`. Likewise, a missing AMR
leaf or interpolation row inside the computational domain is unavailable mesh data.
This distinction is essential: unavailable MHD coverage must not increase calculated
particle access.

## Linked live-versus-replay acceptance

The complete Step-10 acceptance test requires a configured AMPS/SWMF tree. For every
profile below, keep all existing C/F thresholds unchanged:

1. Run a live coupled cutoff/access callback with `SWMF_SNAPSHOT_EXPORT T`.
2. Replay the exported CSV in standalone Mode3D using
   `examples/standalone_step9_swmf_replay.in.template` and the identical product input.
3. Compare the corresponding artifacts:

   ```bash
   python3 srcEarth/test/USWMFCoupledAccess/compare_cutoff_access.py \
     --live  coupled/cutoff_3d_points.swmf_t..._sid....dat \
     --replay replay/cutoff_3d_points....dat \
     --json  comparison.json
   ```

   The default gate is exact. If a campaign has a preregistered nonzero parity
   tolerance, pass it explicitly with `--rtol` and/or `--atol`; this script never
   chooses or expands a tolerance on its own. Zone declarations must also match
   exactly, and SHUE artifacts must carry identical resolved r0, alpha, and tail cap.

4. Repeat for `POINTS`, `TRAJECTORY`, and `SHELLS`, including `DIRECT_ACCESS` for
   `A(E,Omega)` and one scalar cutoff search.
5. Repeat the same snapshot with 1 MPI rank × 1 thread, 2 ranks × 8 threads, and the
   release-default 8 ranks × 16 threads. Include STATIC, DYNAMIC, and BLOCK_CYCLIC
   scheduler profiles where the input supports them.
6. Restart at a saved epoch. The exported snapshot/content IDs, product suffix, numeric
   artifact, and replay comparison must reproduce; a duplicate callback in the same
   process must be skipped collectively.
7. Deliberately test an incomplete receive, stale clock, missing AMR row, invalid Shue
   configuration, and root output failure. Each requested epoch must leave an explicit
   `swmf_product_status*.json` with `FAILED`, never a passing or partial manifest.

The live calculation writes `swmf_cutoff_access_manifest*.json` only after every
recorded cutoff/access file was successfully closed and verified nonempty. Each product
also contains `AUXDATA` for snapshot ID, absolute epoch, mesh revision, content
fingerprint, and outer-boundary policy.

Step 10 remains an instantaneous/quasi-static access calculation: B is frozen over the
trajectory batch. Passing this suite does not validate time-dependent characteristics,
electric acceleration, long-duration trapping across snapshots, or Step-11 coupled
flux/spectrum acceptance.
