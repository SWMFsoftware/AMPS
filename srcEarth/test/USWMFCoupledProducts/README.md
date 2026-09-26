# Step 11: SWMF-coupled flux and spectrum validation

This suite validates the Roadmap Step-11 implementation: the live SWMF callback folds
directional access from the frozen Step-9 field with the same Step-6 boundary spectrum
and product integrator used by standalone Mode3D. It produces local differential
spectra, density, omnidirectional and one-way planar integral flux, configured energy
channels, optional detector-response rates, conservative lower/upper bounds, and a
trajectory-termination summary.

The reference path is the backward-characteristic Mode3D solver replaying the exported
SWMF snapshot. `3d_forward/Density3D.cpp` and `SphereFlux3D.cpp` are forward
Monte-Carlo population samplers and are **not** reference solutions for this step.

## Fast tests

From the AMPS repository root:

```bash
./srcEarth/test/USWMFCoupledProducts/run_test.sh
```

The fast suite requires only a C++17 compiler and Python 3. It returns nonzero on the
first failure and prints `RESULT: PASS` only after all components pass.

| ID | Test | Independent reference or failure condition |
|---|---|---|
| S11-U01 | Coupled-control identity | Deterministic spectrum-map order; fingerprint changes for a detector factor, channel bound, or spacecraft position |
| S11-U02 | Product manifest | Exact field/spectrum epoch, suffix identity, complete artifact roles, termination closure, and unchanged 1% unresolved gate |
| S11-U03 | Open/blocked products | Closed-form constant-spectrum local spectrum, omnidirectional/planar flux, clipped channel, top-hat rate, and zero-access solution |
| S11-U04 | Unresolved bounds | Closed-form 0.49/0.50/0.51 access bounds; unresolved access remains an interval instead of becoming zero |
| S11-U05 | Complete-set replay | Four manufactured density/spectrum/flux/termination artifacts agree exactly |
| S11-U06 | Comparator threshold | A changed numeric value fails the default exact comparison; a tolerance is used only when explicitly supplied |
| S11-U07 | Provenance/schema rejection | Missing spectrum time or configured detector column fails even when remaining numbers match |
| S11-U08 | Release gate rejection | A manifest above its retained unresolved tolerance cannot pass |
| S11-SOURCE | Production wiring | Shared Step-6 kernel, ephemeris selection, collective controls, close-before-manifest ordering, no forward-MC reference, and unchanged F4/C9/C19 entries |

S11-U03 and S11-U04 are numerical reference-solution tests. They do not merely compare
two copies of implementation output. The expected integrals are closed-form constants
stored in the test, while the calculation uses the production `BoundaryProducts`
kernel.

## Product transaction

At each accepted SWMF/PT epoch, the bridge performs this transaction:

1. Assemble and freeze one content-identified B/u snapshot.
2. Set the product epoch to the snapshot's absolute UTC. For `TRAJECTORY`, retain only
   ephemeris samples matching that UTC within 1 ms and remap location-qualified
   apertures.
3. Collectively compare the complete spectrum/channel/response/observation control
   fingerprint across MPI ranks.
4. Call `Earth::Mode3D::RunDensityAndFlux()` directly. This is the same Step-6
   backward-characteristic integrator used by standalone Mode3D.
5. Flush and close every numeric artifact before adding it to the per-call inventory.
6. Require spectrum, density, flux, and termination roles; require termination counts
   to sum to sampled trajectories; apply the configured `DS_UNRESOLVED_TOL` unchanged.
7. Write `swmf_flux_spectrum_manifest*.json`, release the snapshot, and only then write
   `swmf_product_status*.json` with `PASS`.

A failed write, missing role, stale/mismatched epoch, inconsistent rank control,
termination-count defect, or unresolved-gate failure leaves the prewritten status as
`FAILED`. The bridge forces the termination file on for the coupled transaction, but
does not alter trajectory classification, numerical resolution, or any tolerance.

Every numeric artifact includes snapshot ID, absolute UTC, mesh revision, field-content
fingerprint, physical outer-boundary policy, boundary-spectrum evaluation/active-table
epochs, energy basis and units, temporal interpolation status, spectrum uncertainty,
product/channel/response/observation fingerprints, and unresolved accounting.

## Live SWMF versus offline replay (I-F03, I-F07)

The linked acceptance test needs an AMPS/SWMF build and cannot be simulated by the
portable runner:

1. Run the live coupled template with `SWMF_SNAPSHOT_EXPORT T` and
   `CALC_TARGET CUTOFF_RIGIDITY+DENSITY_SPECTRUM`.
2. Replay the exported CSV with standalone `FIELD_MODEL SWMF_SNAPSHOT`, using the same
   species, locations, energy grid, spectrum table, channels, responses, boundary
   policy, and unresolved settings.
3. Compare the complete set:

   ```bash
   python3 srcEarth/test/USWMFCoupledProducts/compare_flux_products.py \
     --live-manifest coupled/swmf_flux_spectrum_manifest.swmf_t..._sid....json \
     --replay-dir replay \
     --json coupled_vs_replay.json
   ```

The default is exact (`--rtol 0 --atol 0`). A preregistered campaign may pass a
nonzero tolerance explicitly, but the comparator never chooses or expands one. It also
requires all product schemas, zone topology, provenance, channel columns, and detector
columns to match.

For I-F03, repeat using an analytic dipole loaded through both gridless and Mode3D with
identical states/directions. Kernel quantities must agree to `1e-10`; energy-integrated
products must agree within 2%, and mesh-backed products within 5% at the declared mesh
resolution. No result may be omitted from the comparison report.

For I-F07, include a time-dependent boundary table and a timestamped trajectory. The
field epoch, boundary-spectrum evaluation epoch, selected/interpolated table rows,
ephemeris sample, product timestamp, and response schema must be identical by
construction. A missing matching trajectory sample is a hard failure.

## Convergence and parallel matrix (F6, F7, F13, I-F06)

Run the same frozen snapshot and product input across the following matrix, retaining
the plan's existing gates:

- Energy-grid and angular-grid refinement (F6): successive integrated density/flux
  changes must be at most 2%; response-supported bins and directional detector products
  must be within 5%.
- Gridless versus Mode3D product parity (F7): kernel values within `1e-10`, integrated
  values within 2%, and declared-resolution mesh products within 5%.
- Parallel reproducibility (F13/I-F06): 1 MPI rank × 1 thread, 2 × 8, and 8 × 16;
  STATIC, DYNAMIC, and BLOCK_CYCLIC schedulers. Deterministic arrays must agree within
  `1e-10`, reductions within `1e-7`, with identical termination-category counts and
  complete artifact comparisons.

These are pass/fail requirements, not suggested tolerances. This implementation does
not modify the corresponding existing test gates.

## Response closure and observations (F17, I-F10, O3)

For I-F10, run synthetic delta-like and top-hat spectra/responses with frozen instrument
schemas. Verify exact fold values, all configured output columns, response fingerprints,
and nonempty comparison records. Missing comparisons fail the test.

For F17, run at least two configured instrument/channel schemas simultaneously. Require
the direct product integral, energy-channel integral, and response-folded rate to close
against independent quadrature. The response-weighted unresolved contribution must be
at most 1%, and the lower/upper unresolved width must be at most 5%. Do not replace an
unresolved interval with zero or drop unsupported energy bins.

O3 remains an observational validation campaign, not a result claimed by this portable
suite. Use the September 2017 interval with GOES-16 SGPS and Van Allen Probes REPT after
recording data versions, energy/response definitions, time averaging, ephemeris/frame
conversion, background/dead-time processing, and uncertainty handling. Compare modeled
response-folded rates and local spectra at matching epochs, report coverage and missing
data explicitly, and retain all preregistered observational gates.

## Scope limitation

Step 11 is instantaneous/quasi-static and magnetic-only: B is frozen for a trajectory
batch and `STATIC_MAGNETIC` mapping is used. Passing these tests does not validate
electric acceleration, time-dependent characteristics, forward population evolution,
or long-duration trapping across field snapshots.
