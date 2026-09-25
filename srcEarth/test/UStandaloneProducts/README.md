# UStandaloneProducts — strict Roadmap Step 7 contract suite

Run from the AMPS repository root:

```bash
./srcEarth/test/UStandaloneProducts/run_test.sh
```

The suite compiles the production `util/StandaloneProductContract.h` with C++11
warnings promoted to errors, then audits its production wiring. It has no AMPS, MPI,
SPICE, Geopack, empirical-field, or SWMF link dependency.

## What the suite validates

| Gate | Check | Reference or exact expectation |
|---|---|---|
| S7-C01 | Model aliases and roles | Exact canonical table for DIPOLE, IGRF, T96, T01, T05/TS05, TA15N/B, and TA16; `NONE` is explicitly validation-only; T89 is rejected. |
| S7-C02 | Driver mapping | Exact wrapper inputs: T01 adds G1–G3, T05 adds W1–W6, and TA15 consumes only PDYN/BY/BZ/XIND. |
| S7-C03 | Native units and values | nT, nPa, km/s, cm^-3, and dimensionless aliases; Pa where nPa is required, partial numeric tokens, NaN, and infinity must fail. |
| S7-C04 | Product/domain selection | Cutoff-only, flux/spectrum-only, and combined products across POINTS, TRAJECTORY, and SHELLS; unknown partial tokens fail. |
| S7-C05 / I-F07 contract | Snapshot synchronization | Field, driver record, boundary spectrum, ephemeris, and output epoch must match exactly; missing snapshot identity, driver coverage, Geopack initialization, or field validity fails. |
| S7-C06 | Manifest | Exact schema fields, JSON escaping, snapshot ID, selected products, validation status, and `NONE` labeling. |
| S7-W01 | Production wiring | Both gridless and Mode3D consume the shared contract, emit manifests, and enforce immutable snapshot identity. |
| S7-W02 | SWMF isolation | T01/TA15 headers remain inside the non-SWMF compile guard; Step 7 adds no empirical-field call to coupled builds. |

These are startup/contract tests, not substitutes for physical reference solutions.
The numerical references remain deliberately independent:

- UFluxNumerics, UFieldProvider, UTrajectoryCore, UDirectionalAccess, and
  UBoundaryProducts test the shared numerical kernels against analytic solutions.
- F1/F2/F15/F16 test exact zero-field spectra, flux, and density.
- C1–C7/C11/C14 and F3–F5/F11/F12 test dipole/IGRF/phenomenological fields,
  transmission, cross-solver behavior, and reconstruction.
- C8/C9/C10/C19 retain their directional and observation-facing gates.

No existing expected result, tolerance, mover, trace limit, input, reference file, or
`last pass:` record is modified by this suite.
